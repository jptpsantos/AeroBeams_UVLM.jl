# One integration through the reversal. Reuse the independent driver's helpers.
if !isdefined(@__MODULE__, :pazy_lco_aligned_timestep)
    include(joinpath(@__DIR__, "PazyLinearLCOSweep.jl"))
end

function pazy_continuous_airspeed(t; initial_speed, peak_speed, upsweep_duration, downsweep_duration)
    if t <= upsweep_duration
        return pazy_linear_airspeed(t; start_speed=initial_speed,
            end_speed=peak_speed,duration=upsweep_duration)
    end
    return pazy_linear_airspeed(t-upsweep_duration; start_speed=peak_speed,
        end_speed=initial_speed,duration=downsweep_duration)
end

function pazy_continuous_response_row(t,speed,structure,angle_of_attack,gauge_y,gauge_z)
    root = structure.model.elements[1].compStates
    twist = wingtip_twist_degrees(structure.model)
    strain = 1e6*(root.γ[1]+gauge_z*root.κ[2]-gauge_y*root.κ[3])
    return [t,speed,strain,root.γ[1],root.κ[2],root.κ[3],
        -structure.model.elements[end].nodalStates.u_n2[1],twist,rad2deg(angle_of_attack)+twist]
end

function pazy_continuous_branch(history,indices,name)
    data = history[indices,:]
    return (; name,indices,time=data[:,1],branch_time=data[:,1].-data[1,1],
        airspeed_history=data[:,2],root_strain_ue=data[:,3],
        tip_out_of_plane=data[:,7],tip_twist_degrees=data[:,8],wingtip_aoa_degrees=data[:,9])
end

"""Integrate initial -> peak -> initial with one static initialization and clock.

The wing, velocities, circulation and wake evolve through the reversal without
reinitialization, a hold, or an additional pulse. Both CSV branches share the
same turning-point sample. `time_s` is global; `branch_time_s` is local.
Set `fixed_time_step` to match an independent run exactly. Completed upsweep
data and a turnaround checkpoint are saved before starting the downsweep.
"""
function run_pazy_continuous_LCO(; initial_speed=35.0,peak_speed=52.0,
    upsweep_duration=40.0,downsweep_duration=20.0,
    density=1.225,angle_of_attack=deg2rad(5.0),
    chordwise_panels=10,spanwise_panels=20,time_step_chords=0.25,
    fixed_time_step=nothing,wake_length_chords=5.0,gauge_y=0.0,gauge_z=2.25e-3/2,
    coupling_scheme=:strong,coupling_maximum_iterations=20,coupling_relaxation=0.5,
    coupling_geometry_tolerance=1e-4,coupling_load_tolerance=1e-4,
    progress_frequency=100,equilibrium_options=(;),
    output_directory=joinpath(@__DIR__,"output","LCO_continuous_sweeps"))

    @assert all(isfinite,(initial_speed,peak_speed,upsweep_duration,downsweep_duration,
        time_step_chords,wake_length_chords,gauge_y,gauge_z))
    @assert 0 < initial_speed < peak_speed
    @assert upsweep_duration > 0 && downsweep_duration > 0
    @assert time_step_chords > 0 && wake_length_chords > 0
    _,_,chord,_ = AeroBeams.geometrical_properties_Pazy()
    dt = pazy_lco_aligned_timestep(upsweep_duration,downsweep_duration,
        time_step_chords*chord/peak_speed;fixed_time_step)
    turnaround_index = round(Int,upsweep_duration/dt)+1
    n_steps = round(Int,(upsweep_duration+downsweep_duration)/dt)
    @assert turnaround_index >= 2 && n_steps+1 > turnaround_index
    maximum_wake_rows = ceil(Int,wake_length_chords*chord/(initial_speed*dt))
    static_settings = merge((; maximum_aerodynamic_time=0.20,minimum_aerodynamic_time=0.05,
        aerodynamic_window_time=0.01,aerodynamic_force_tolerance=1e-4,
        aerodynamic_circulation_tolerance=1e-4,maximum_static_iterations=100,
        static_relaxation=0.3,geometry_tolerance=1e-5,load_tolerance=1e-4,
        consecutive_equilibrium_iterations=2,newton_maximum_iterations=50,
        newton_absolute_tolerance=1e-8,newton_relative_tolerance=1e-8,
        newton_display_iterations=false,progress_frequency),equilibrium_options)
    mkpath(output_directory)
    settings_path = joinpath(output_directory,"pazy_linear_LCO_settings.txt")
    open(settings_path,"w") do io
        println(io,"Continuous simulation: one coupled static initialization, one dynamic integration, no state reset at reversal.")
        println(io,"Upsweep: $initial_speed -> $peak_speed m/s in $upsweep_duration s")
        println(io,"Downsweep: $peak_speed -> $initial_speed m/s in $downsweep_duration s")
        println(io,"Root incidence [deg]: ",rad2deg(angle_of_attack))
        println(io,"dt [s]: $dt; wake rows: $maximum_wake_rows; mesh: $chordwise_panels x $spanwise_panels")
        println(io,"Gauge y/z [m]: $gauge_y / $gauge_z; no pulse; gravity zero; first-element midpoint strain.")
        println(io,"time_s is global in all CSVs; branch_time_s is local. The turning sample belongs to both exported branches.")
        println(io,"Static controls: ",static_settings)
    end
    equilibrium = run_pazy_wing_uvlm_static(; static_settings...,airspeed=initial_speed,
        density,angle_of_attack,chordwise_panels,spanwise_panels,symmetric_wing=true,
        maximum_wake_rows,time_step_chords=dt*initial_speed/chord,
        core_radius=1e-3,wake_shedding_fraction=0.1,aerodynamic_interaction=true)
    writedlm(joinpath(output_directory,"pazy_linear_LCO_initialization.csv"),
        vcat(permutedims(equilibrium.history_columns),equilibrium.history),',')
    open(settings_path,"a") do io
        println(io,"Initialization: $(equilibrium.termination_reason); converged=$(equilibrium.converged); r_geometry=$(equilibrium.geometry_residual); r_load=$(equilibrium.load_residual)")
    end
    equilibrium.converged || error("Cannot start continuous sweep: $(equilibrium.termination_reason)")
    profile(t) = pazy_continuous_airspeed(t;initial_speed,peak_speed,upsweep_duration,downsweep_duration)
    header = ["time_s","airspeed_mps","root_longitudinal_strain_ue","root_axial_strain",
        "root_curvature_2_per_m","root_curvature_3_per_m","tip_out_of_plane_m","tip_twist_deg","wingtip_aoa_deg"]
    history = zeros(n_steps+1,length(header))
    history[1,:] .= pazy_continuous_response_row(0.0,initial_speed,
        equilibrium.structural_problem,angle_of_attack,gauge_y,gauge_z)
    turnaround = Ref{Any}(nothing)
    first_downsweep_step = Ref{Any}(nothing)
    function write_branch(indices,name)
        data = history[indices,:]
        writedlm(joinpath(output_directory,"pazy_linear_LCO_$name.csv"),
            vcat(permutedims(vcat(header,["branch_time_s"])),
                hcat(data,data[:,1].-data[1,1])),',')
    end
    function observe(i,t,structure,aerodynamic)
        history[i,:] .= pazy_continuous_response_row(t,profile(t),structure,
            angle_of_attack,gauge_y,gauge_z)
        if i == turnaround_index
            turnaround[] = (; index=i,time=t,structural_state=copy(structure.x),
                maximum_velocity=maximum(norm(vcat(e.states.V,e.states.Ω)) for e in structure.model.elements),
                aerodynamic_snapshot=UVLM.snapshot_uvlm(aerodynamic.state.system),
                aerodynamic_time=aerodynamic.state.timeNow,step_index=aerodynamic.state.stepIndex)
            write_branch(1:i,"upsweep")
            writedlm(joinpath(output_directory,"pazy_linear_LCO_continuous_history.csv"),
                vcat(permutedims(header),history[1:i,:]),',')
            marker_header = vcat(header,["maximum_section_velocity_norm","aerodynamic_time_s",
                "aerodynamic_step_index","retained_wake_rows","circulation_norm"])
            marker = vcat(vec(history[i,:]),[turnaround[].maximum_velocity,
                aerodynamic.state.timeNow,aerodynamic.state.stepIndex,
                aerodynamic.state.system.nwake[1],norm(aerodynamic.state.system.Γ)])
            writedlm(joinpath(output_directory,"pazy_linear_LCO_turnaround.csv"),
                vcat(permutedims(marker_header),permutedims(marker)),',')
            figure = pazy_LCO_strain_plot(((;airspeed=history[1:i,2],
                strain=history[1:i,3],label="Upsweep",color=:blue),))
            Plots.savefig(figure,joinpath(output_directory,"pazy_linear_LCO_root_strain_vs_airspeed.png"))
            open(settings_path,"a") do io
                println(io,"Upsweep completed at global t=$t s. Saved live turnaround state: max velocity=$(turnaround[].maximum_velocity); wake rows=$(aerodynamic.state.system.nwake[1]); aerodynamic step=$(aerodynamic.state.stepIndex). Continuing the same integration.")
            end
            println("Continuous upsweep saved at t=$t s; proceeding through reversal with accepted beam/wake state.")
            flush(stdout)
        elseif i == turnaround_index+1
            first_downsweep_step[] = (;time=t,aerodynamic_time=aerodynamic.state.timeNow,
                step_index=aerodynamic.state.stepIndex,nwake=copy(aerodynamic.state.system.nwake))
        end
    end
    result = run_pazy_wing_uvlm(; airspeed=peak_speed,airspeed_profile=profile,
        initial_equilibrium=equilibrium,accepted_step_callback=observe,
        density,angle_of_attack,sideslip=0.0,initial_airspeed_fraction=1.0,
        airspeed_ramp_duration=0.0,duration=upsweep_duration+downsweep_duration,
        settling_time=0.0,perturbation_amplitude=0.0,perturbation_duration=0.1,
        chordwise_panels,spanwise_panels,symmetric_wing=true,
        time_step_chords=dt*peak_speed/chord,maximum_wake_rows,
        core_radius=1e-3,wake_shedding_fraction=0.1,aerodynamic_interaction=true,
        save_uvlm_history=false,uvlm_save_frequency=1,newton_maximum_iterations=20,
        newton_absolute_tolerance=1e-5,newton_relative_tolerance=1e-5,
        newton_display_iterations=false,newton_always_update_jacobian=true,
        coupling_scheme,coupling_maximum_iterations,coupling_relaxation,
        coupling_geometry_tolerance,coupling_load_tolerance,coupling_display_iterations=false,
        animation_frames=2,save_animation_history=false,save_structural_history=false,
        save_aerodynamic_load_history=false,progress_frequency)
    write_branch(turnaround_index:size(history,1),"downsweep")
    writedlm(joinpath(output_directory,"pazy_linear_LCO_continuous_history.csv"),
        vcat(permutedims(header),history),',')
    upsweep = pazy_continuous_branch(history,1:turnaround_index,"upsweep")
    downsweep = pazy_continuous_branch(history,turnaround_index:size(history,1),"downsweep")
    figure = pazy_LCO_strain_plot((
        (;airspeed=upsweep.airspeed_history,strain=upsweep.root_strain_ue,label="Upsweep",color=:blue),
        (;airspeed=downsweep.airspeed_history,strain=downsweep.root_strain_ue,label="Downsweep",color=:red)))
    Plots.savefig(figure,joinpath(output_directory,"pazy_linear_LCO_root_strain_vs_airspeed.png"))
    trace = Plots.plot(result.time,history[:,3];xlabel="Global time [s]",
        ylabel="Root strain [microstrain]",label="Strain",color=:blue,linewidth=0.6,
        right_margin=8Plots.mm)
    Plots.vline!(trace,[upsweep_duration];label="Speed reversal",color=:gray,linestyle=:dash)
    Plots.plot!(Plots.twinx(trace),result.time,result.airspeed_history;
        ylabel="Airspeed [m/s]",label=false,color=:orangered,grid=false,linewidth=1.2)
    Plots.savefig(trace,joinpath(output_directory,"pazy_linear_LCO_root_strain_history.png"))
    open(settings_path,"a") do io
        println(io,"Continuous sweep completed at t=$(last(result.time)) s; accepted steps=$(length(result.time)-1); final wake rows=$(result.aerodynamic_problem.state.system.nwake[1]).")
        println(io,"Maximum accepted coupling residuals: geometry=$(maximum(result.coupling_geometry_residual)), load=$(maximum(result.coupling_load_residual)).")
    end
    println("Continuous LCO outputs saved in: $output_directory")
    return (;result,equilibrium,upsweep,downsweep,turnaround=turnaround[],
        first_downsweep_step=first_downsweep_step[],dt,maximum_wake_rows,output_directory)
end
