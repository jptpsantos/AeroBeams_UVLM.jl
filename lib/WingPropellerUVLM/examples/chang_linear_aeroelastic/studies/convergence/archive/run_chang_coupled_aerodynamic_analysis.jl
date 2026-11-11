# Rigid aerodynamic analysis of the combined Chang wing--propeller geometry.
#
# This is the inexpensive first stage of the UVLM convergence workflow.  The
# wing is held rigid, the propeller rotates at constant advance ratio, and the
# free wakes are marched until periodic loads are obtained.  Wing lift and
# propeller thrust are reduced from the same dimensional Imperial near-field
# nodal loads used by the production aeroelastic driver.
#
# Run from the repository root with
#
#   julia --project=lib/WingPropellerUVLM lib/WingPropellerUVLM/examples/chang_linear_aeroelastic/studies/convergence/run_chang_coupled_aerodynamic_analysis.jl
#
# Every numerical control can be overridden from PowerShell, for example:
#
#   $env:CHANG_AERO_INTERACTION = "true"
#   $env:CHANG_AERO_FCORE_SEGMENT_FACTOR = "0.25"
#   $env:CHANG_AERO_SIMULATED_REVOLUTIONS = "4"

using Dates
using DelimitedFiles
using LinearAlgebra
using Printf
using Statistics
using StaticArrays

const CONVERGENCE_DIR = @__DIR__
const EXAMPLE_DIR = normpath(joinpath(CONVERGENCE_DIR,"..","..",".."))

using WingPropellerUVLM:
    Freestream,
    Reference,
    Uniform,
    copy_surfaces_to_previous!,
    get_nodal_properties_chang,
    grid_to_surface_panels,
    imperial_nodal_forces,
    imperial_nodal_positions,
    initialize_bohnisch_uvlm_system,
    near_field_forces!,
    propagate_system!,
    update_propeller_grids!,
    update_system_surfaces!

include(joinpath(@__DIR__,"chang_aerodynamic_study.jl"))
using .ChangAerodynamicStudy

"""Return lift in the stability-axis vertical direction."""
function stability_lift(bodyForce,angleOfAttackRad)
    sine,cosine = sincos(angleOfAttackRad)
    return -sine * bodyForce[1] + cosine * bodyForce[3]
end

function sum_surface_force(nodalForce)
    total = zero(first(nodalForce))
    for force in nodalForce
        total += force
    end
    return total
end

function propeller_loads(nodal_forces,nodalPositions,propellerSurfaceIndices,hub)
    totalForce = SVector(0.0,0.0,0.0)
    totalMoment = SVector(0.0,0.0,0.0)
    for surface_index in propellerSurfaceIndices
        for vertexIndex in eachindex(nodal_forces[surface_index])
            force = nodal_forces[surface_index][vertexIndex]
            position = nodalPositions[surface_index][vertexIndex]
            totalForce += force
            totalMoment += cross(position - hub,force)
        end
    end
    return totalForce,totalMoment
end

"""March the rigid interacting wing--propeller system and return periodic loads."""
function simulate_coupled_aerodynamics(options; verbose=true)
    source_sha256 = ChangAerodynamicStudy.source_fingerprint()
    steps_per_revolution = validate_options(options)
    omega = options.reference_rpm * 2pi / 60 *
            options.flow_speed_mps / options.reference_speed_mps
    rpm = omega * 60 / (2pi)
    dt = deg2rad(options.azimuth_step_deg) / omega
    total_steps = options.simulated_revolutions * steps_per_revolution
    maximum_wake_rows = max(
        1,
        ceil(Int,options.retained_wake_revolutions * steps_per_revolution)
    )

    root_chord = options.wing_root_chord_m
    tip_chord = options.wing_tip_chord_m
    span = options.wing_span_m
    span_nodes = collect(range(0.0,span,length=options.wing_span_panels + 1))
    chordDistribution = collect(range(
        root_chord,
        tip_chord,
        length=options.wing_span_panels + 1
    ))
    xle_tip = 0.5 * (root_chord - tip_chord)
    xle_distribution = collect(range(
        0.0,
        xle_tip,
        length=options.wing_span_panels + 1
    ))
    propeller_span_position = options.propeller_attachment_eta * span
    propellerAttachNode = argmin(abs.(span_nodes .- propeller_span_position))

    _,_,twists_at_nodes = get_nodal_properties_chang(options.propeller_radial_panels)
    blade_twists = deg2rad.(
        twists_at_nodes .- 90.0 .+ options.collective_pitch_offset_deg,
    )
    fcore = if isnothing(options.core_radius_m)
        (chordValue,segmentLength) -> max(
            options.finite_core_segment_factor * segmentLength,
            options.finite_core_chord_factor * chordValue
        )
    else
        radius = Float64(options.core_radius_m)
        (chordValue,segmentLength) -> radius
    end
    reference = Reference(
        span * (root_chord + tip_chord) / 2,
        0.5 * (root_chord + tip_chord),
        span,
        SVector(0.3 * root_chord,0.0,0.0),
        options.flow_speed_mps,
        options.air_density_kgpm3
    )
    freestream = Freestream(
        options.flow_speed_mps,
        deg2rad(options.angle_of_attack_deg),
        deg2rad(options.sideslip_deg),
        SVector(0.0,0.0,0.0)
    )
    hubCenter = SVector(-options.pylon_length_m,0.0,0.0)

    uvlm = initialize_bohnisch_uvlm_system(
        xle=[0.0,xle_tip],
        yle=[0.0,span],
        zle=[0.0,0.0],
        chord_geo=[root_chord,tip_chord],
        theta_geo=[0.0,0.0],
        phi_geo=[0.0,0.0],
        ns_wing=options.wing_span_panels,
        nc_wing=options.wing_chord_panels,
        mirror_wing=false,
        spacing_s_wing=Uniform(),
        spacing_c_wing=Uniform(),
        R_prop=options.propeller_radius_m,
        c_prop=options.propeller_chord_m,
        ns_prop=options.propeller_radial_panels,
        nc_prop=options.propeller_chord_panels,
        blade_twists_prop=blade_twists,
        Nb_prop=options.propeller_blades,
        Npropellers=1,
        span_nodes=span_nodes,
        prop_attach_nodes=[propellerAttachNode],
        propeller_span_positions=[propeller_span_position],
        chord=chordDistribution,
        xle_distribution=xle_distribution,
        ref=reference,
        symmetric_wing=options.wing_symmetric,
        fs=freestream,
        dt=fill(dt,total_steps),
        nnodes=length(span_nodes),
        prop_pivot_offset_from_ea_A=SVector(0.0,0.0,0.0),
        hub_center_prop_A=hubCenter,
        fcore=fcore,
        elastic_axis_fraction=options.elastic_axis_fraction,
        maximum_wake_rows_wing=maximum_wake_rows,
        maximum_wake_rows_propeller=maximum_wake_rows,
        verbose=verbose
    )

    system = uvlm.system
    active_wake_rows = zeros(Int,length(system.surfaces))
    propellerSurfaceIndices = uvlm.prop_surface_indices[1]
    zeroNodes = zeros(length(span_nodes))
    zeroPropellerCoordinates = zeros(2)

    wing_cl = Vector{Float64}(undef,total_steps)
    propeller_ct = Vector{Float64}(undef,total_steps)
    wing_force_history = Matrix{Float64}(undef,total_steps,3)
    propeller_force_history = Matrix{Float64}(undef,total_steps,3)
    propeller_torque = Vector{Float64}(undef,total_steps)

    dynamic_pressure = 0.5 * options.air_density_kgpm3 * options.flow_speed_mps^2
    wingArea = span * (root_chord + tip_chord) / 2
    rotationRateHz = omega / (2pi)
    diameter = 2 * options.propeller_radius_m
    ctDenominator = options.air_density_kgpm3 * rotationRateHz^2 * diameter^4
    cqDenominator = options.air_density_kgpm3 * rotationRateHz^2 * diameter^5
    hub = SVector(
        uvlm.ea_x_aero[1],
        uvlm.attach_node_y[1],
        0.0
    ) + hubCenter

    verbose && println("Rigid coupled aerodynamic analysis")
    coreDescription =
        isnothing(options.core_radius_m) ?
        "max($(options.finite_core_segment_factor) ds, $(options.finite_core_chord_factor) c)" :
        "fixed $(options.core_radius_m) m"
    verbose && println("  Core radius: $coreDescription")
    verbose && @printf(
        "  V=%.3f m/s, RPM=%.3f, J=%.6f, grids wing=%dx%d prop=%dx%d\n",
        options.flow_speed_mps,
        rpm,
        options.flow_speed_mps / (rotationRateHz * diameter),
        options.wing_span_panels,
        options.wing_chord_panels,
        options.propeller_radial_panels,
        options.propeller_chord_panels,
    )
    verbose && @printf(
        "  interaction=%s, wake=%.3f rev\n",
        options.interaction_on,
        options.retained_wake_revolutions,
    )
    verbose && flush(stdout)

    core_diagnostics = [
        begin
            radii = [panel.core_size for panel in surface]
            widths = [norm(panel.rtr - panel.rtl) for panel in surface]
            chords = [panel.chord for panel in surface]
            (; minimum_m=minimum(radii),maximum_m=maximum(radii),
                max_radius_over_span_edge=maximum(radii ./ widths),
                max_radius_over_panel_chord=maximum(radii ./ chords))
        end for surface in system.surfaces
    ]

    for step = 1:total_steps
        copy_surfaces_to_previous!(system,length(system.surfaces))
        time_s = step * dt
        update_propeller_grids!(
            uvlm.grids_prop_current,
            uvlm.T_pivot_A_current,
            uvlm.grids_prop_ref,
            zeroPropellerCoordinates,
            [propellerAttachNode],
            uvlm.ea_x_aero,
            uvlm.attach_node_y,
            zeroNodes,
            zeroNodes,
            zeroNodes,
            zeroNodes,
            SVector(0.0,0.0,0.0),
            hubCenter,
            omega,
            time_s,
            1,
            options.propeller_blades
        )
        update_system_surfaces!(
            system,
            uvlm.grid_wing_init,
            uvlm.grids_prop_current,
            uvlm.ratio_wing,
            1,
            options.propeller_blades;
            fcore
        )

        propagate_system!(
            system,
            freestream,
            dt;
            additional_velocity=nothing,
            repeated_points=uvlm.repeated_points,
            nwake=active_wake_rows,
            eta=options.wake_relaxation,
            calculate_influence_matrix=true,
            near_field_analysis=true,
            near_field_force_function=near_field_forces!,
            derivatives=false,
            interaction_id=uvlm.surface_interaction_id,
            interaction=options.interaction_on,
            advance_wake=true
        )
        active_wake_rows .= min.(active_wake_rows .+ 1,maximum_wake_rows)

        nodal_forces = imperial_nodal_forces(system)
        nodalPositions = imperial_nodal_positions(system)
        wingForce = sum_surface_force(nodal_forces[1])
        propForce,propMoment = propeller_loads(
            nodal_forces,
            nodalPositions,
            propellerSurfaceIndices,
            hub
        )
        wing_force_history[step,:] .= wingForce
        propeller_force_history[step,:] .= propForce
        propeller_torque[step] = propMoment[1]
        wing_cl[step] = stability_lift(wingForce,deg2rad(options.angle_of_attack_deg)) /
                        (dynamic_pressure * wingArea)
        propeller_ct[step] = -propForce[1] / ctDenominator

        all(isfinite,(wing_cl[step],propeller_ct[step],propeller_torque[step])) ||
            error("Nonfinite aerodynamic coefficients at step $step")

        if verbose && step % steps_per_revolution == 0
            revolution = step ÷ steps_per_revolution
            indices = (step-steps_per_revolution+1):step
            @printf(
                "  revolution %d/%d: mean wing CL=%+.7f, propeller CT=%+.7f\n",
                revolution,
                options.simulated_revolutions,
                mean(wing_cl[indices]),
                mean(propeller_ct[indices]),
            )
            flush(stdout)
        end
    end

    revolution_mean_cl = [
        mean(wing_cl[((rev-1)*steps_per_revolution+1):(rev*steps_per_revolution)])
        for rev = 1:options.simulated_revolutions
    ]
    revolution_mean_ct = [
        mean(propeller_ct[((rev-1)*steps_per_revolution+1):(rev*steps_per_revolution)])
        for rev = 1:options.simulated_revolutions
    ]
    averageSteps = options.averaged_revolutions * steps_per_revolution
    averageIndices = (total_steps-averageSteps+1):total_steps

    return (;
        options,
        core_diagnostics,
        source_sha256,
        omega_radps=omega,
        rpm,
        advance_ratio=options.flow_speed_mps / (rotationRateHz * diameter),
        time_step_s=dt,
        steps_per_revolution,
        maximum_wake_rows,
        time_s=collect(1:total_steps) .* dt,
        azimuth_deg=(collect(1:total_steps) .* options.azimuth_step_deg) .% 360,
        wing_cl,
        propeller_ct,
        wing_force_history,
        propeller_force_history,
        propeller_torque,
        propeller_cq=propeller_torque ./ cqDenominator,
        mean_wing_cl=mean(wing_cl[averageIndices]),
        mean_propeller_ct=mean(propeller_ct[averageIndices]),
        std_wing_cl=std(wing_cl[averageIndices]),
        std_propeller_ct=std(propeller_ct[averageIndices]),
        revolution_mean_cl,
        revolution_mean_ct,
        last_revolution_cl_change=revolution_mean_cl[end] - revolution_mean_cl[end-1],
        last_revolution_ct_change=revolution_mean_ct[end] - revolution_mean_ct[end-1]
    )
end

function write_results(result,output_directory; plot_results=true)
    mkpath(output_directory)
    history_path = joinpath(output_directory,"coupled_aerodynamic_history.csv")
    revolution_path = joinpath(output_directory,"coupled_aerodynamic_revolutions.csv")
    summary_path = joinpath(output_directory,"coupled_aerodynamic_summary.txt")
    plot_path = joinpath(output_directory,"coupled_aerodynamic_coefficients.png")

    open(history_path,"w") do stream
        println(
            stream,
            "time_s,azimuth_deg,wing_CL,propeller_CT,propeller_CQ," *
            "wing_Fx_N,wing_Fy_N,wing_Fz_N,propeller_Fx_N,propeller_Fy_N,propeller_Fz_N"
        )
        writedlm(
            stream,
            hcat(
                result.time_s,
                result.azimuth_deg,
                result.wing_cl,
                result.propeller_ct,
                result.propeller_cq,
                result.wing_force_history,
                result.propeller_force_history
            ),
            ','
        )
    end
    open(revolution_path,"w") do stream
        println(stream,"revolution,mean_wing_CL,mean_propeller_CT")
        writedlm(
            stream,
            hcat(
                collect(1:result.options.simulated_revolutions),
                result.revolution_mean_cl,
                result.revolution_mean_ct
            ),
            ','
        )
    end
    open(summary_path,"w") do stream
        options = result.options
        println(stream,"Rigid coupled Chang wing--propeller aerodynamic analysis")
        println(stream,"Generated: $(Dates.format(now(),"yyyy-mm-dd HH:MM:SS"))")
        println(stream,"Core radius by surface (m and mesh ratios): $(result.core_diagnostics)")
        println(stream,"Wake shedding fraction eta: $(options.wake_relaxation)")
        println(stream,"Interaction enabled: $(options.interaction_on)")
        println(stream,"Wing image symmetry across y=0: $(options.wing_symmetric); blade symmetry: false")
        println(stream,"Fixed core radius (m; nothing uses factors): $(options.core_radius_m)")
        @printf(stream,"Flow speed: %.8f m/s\n",options.flow_speed_mps)
        @printf(stream,"Angle of attack: %.8f deg\n",options.angle_of_attack_deg)
        @printf(stream,"RPM: %.8f\n",result.rpm)
        @printf(stream,"Advance ratio J: %.10f\n",result.advance_ratio)
        @printf(stream,"Wing mesh: %d span x %d chord panels\n",
            options.wing_span_panels,options.wing_chord_panels)
        @printf(stream,"Propeller mesh: %d radial x %d chord panels per blade\n",
            options.propeller_radial_panels,options.propeller_chord_panels)
        @printf(stream,"Azimuth step: %.8f deg; dt: %.12g s\n",
            options.azimuth_step_deg,result.time_step_s)
        @printf(stream,"Simulated/averaged revolutions: %d/%d\n",
            options.simulated_revolutions,options.averaged_revolutions)
        @printf(stream,"Retained wake: %.8f revolutions (%d rows)\n",
            options.retained_wake_revolutions,result.maximum_wake_rows)
        if isnothing(options.core_radius_m)
            @printf(stream,"Corrected finite core: max(%.8g ds, %.8g c)\n",
                options.finite_core_segment_factor,options.finite_core_chord_factor)
        else
            @printf(stream,"Corrected finite core: fixed %.8g m\n",options.core_radius_m)
        end
        @printf(stream,"Mean wing CL: %+.12e\n",result.mean_wing_cl)
        @printf(stream,"Wing CL standard deviation: %.12e\n",result.std_wing_cl)
        @printf(stream,"Mean propeller CT: %+.12e\n",result.mean_propeller_ct)
        @printf(stream,"Propeller CT standard deviation: %.12e\n",result.std_propeller_ct)
        @printf(stream,"Last two revolution mean-CL change: %+.12e\n",
            result.last_revolution_cl_change)
        @printf(stream,"Last two revolution mean-CT change: %+.12e\n",
            result.last_revolution_ct_change)
        println(stream,"Mean wing CL by revolution: $(result.revolution_mean_cl)")
        println(stream,"Mean propeller CT by revolution: $(result.revolution_mean_ct)")
    end

    write_metadata(output_directory,result.options; expected_fingerprint=result.source_sha256)
    plot_results || return (; history_path,revolution_path,summary_path,plot_path="")
    @eval using Plots
    return Base.invokelatest(write_coefficient_plot,result,output_directory,
        history_path,revolution_path,summary_path,plot_path)
end

function write_coefficient_plot(result,output_directory,history_path,revolution_path,summary_path,plot_path)
    firstAverageStep = length(result.time_s) -
                       result.options.averaged_revolutions * result.steps_per_revolution + 1
    plotIndices = firstAverageStep:length(result.time_s)
    coefficientPlot = plot(
        result.time_s[plotIndices] .* result.rpm ./ 60,
        result.wing_cl[plotIndices];
        xlabel="Rotor revolutions",
        ylabel="Coefficient",
        label="wing CL",
        linewidth=2,
        framestyle=:box,
        gridalpha=0.25
    )
    plot!(
        coefficientPlot,
        result.time_s[plotIndices] .* result.rpm ./ 60,
        result.propeller_ct[plotIndices];
        label="propeller CT",
        linewidth=2
    )
    savefig(coefficientPlot,plot_path)
    return (; history_path,revolution_path,summary_path,plot_path)
end

function main()
    options = options_from_environment()
    result = simulate_coupled_aerodynamics(options)
    runStamp = Dates.format(now(),"yyyymmdd_HHMMSS")
    output_directory = abspath(get(
        ENV,
        "CHANG_AERO_OUTPUT_DIR",
        joinpath(EXAMPLE_DIR,"output","coupled_aerodynamic_$runStamp")
    ))
    paths = write_results(result,output_directory;
        plot_results=aero_bool("CHANG_AERO_PLOT_RESULTS",true))

    println("\nFinal-window aerodynamic coefficients (check periodicity in the sweep)")
    @printf("  wing CL       = %+.10f\n",result.mean_wing_cl)
    @printf("  propeller CT  = %+.10f\n",result.mean_propeller_ct)
    @printf("  delta CL/rev  = %+.3e\n",result.last_revolution_cl_change)
    @printf("  delta CT/rev  = %+.3e\n",result.last_revolution_ct_change)
    println("  summary: $(paths.summary_path)")
    println("  history: $(paths.history_path)")
    println("  plot:    $(paths.plot_path)")
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
