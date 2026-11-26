# Steady Pazy response for swept-wing native AeroBeams and UVLM comparisons.
ENV["GKSwstype"] = get(ENV,"GKSwstype","100")
include(joinpath(@__DIR__,"PazyWingUVLMStatic.jl"))
include(joinpath(@__DIR__,"PazySteadyNativeReference.jl"))
using DelimitedFiles
import Plots

const SWEPT_PITCH_RANGE = (0,1,3,5,7,10)
const SWEPT_ANGLES_DEG = (10,20)
const SWEPT_AIRSPEEDS = [1e-3; collect(10.0:10.0:70.0)]
const SWEPT_STUDY_DIRECTORY = joinpath(@__DIR__,"output","swept_steady_benchmark_vlm_def")
const SWEPT_UVLM_OPTIONS = (;
    chordwise_panels=4,spanwise_panels=15,symmetric_wing=true,maximum_wake_rows=40,
    core_radius=1e-6,wake_shedding_fraction=0.1,aerodynamic_interaction=true,
    aerodynamic_force_tolerance=1e-4,aerodynamic_circulation_tolerance=1e-4,
    maximum_static_iterations=100,static_relaxation=1.0,geometry_tolerance=1e-5,
    load_tolerance=1e-4,load_residual_force_floor=1e-6,
    consecutive_equilibrium_iterations=2,newton_maximum_iterations=50,
    newton_absolute_tolerance=1e-8,newton_relative_tolerance=1e-8,
    newton_display_iterations=false,progress_frequency=1000,
    relaxation_mode=:anderson,anderson_depth=4,geometry_cap_chords=2.0,wake_reuse=true)

function create_swept_pazy_model(airspeed,pitch_angle,sweep_angle; hasTipCorrection=true)
    return AeroBeams.create_Pazy(; mesh=:nominal,sweepAngle=sweep_angle,
        upright=true,pitchAngle=pitch_angle,airspeed,withWingtipRod=true,
        beamWingtipRod=false,sweptTipMassConfig=:TE,
        hasTipCorrection,tipLossType="VLM-def",tipLossFunctionIsAirspeedDependent=hasTipCorrection,
        aeroSolver=AeroBeams.Indicial(),
        airfoil=deepcopy(AeroBeams.FlatPlate))
end

function create_swept_pazy_uvlm_structure(sweep_angle)
    return create_pazy_uvlm_structure(
        create_swept_pazy_model(0.0,0.0,sweep_angle; hasTipCorrection=false)...;
        parameters=(; mesh=:nominal,sweepAngle=sweep_angle,sweptTipMassConfig=:TE))
end

function swept_te_strain(model)
    wing = model.beams[1]
    nodes = wing.discretization.normalizedNodalPositions
    midpoints = (nodes[1:end-1] .+ nodes[2:end]) ./ 2
    gauge = argmin(abs.(midpoints .- 0.40))
    strain = wing.elements[gauge].state.midpoint
    # Gauge coordinates on the upper surface at the aft edge of the spar.
    return 1e6 * (strain.γ[1] - strain.κ[2] * (2.23e-3 / 2) + strain.κ[3] * (50e-3 / 2))
end

function swept_uvlm_tip_aoa(problem)
    system = problem.state.system
    grid = system.grids[1]
    chord = normalize(grid[:,end,end] - grid[:,1,end])
    tipSpan = grid[:,1,end] + grid[:,end,end] - grid[:,1,end-1] - grid[:,end,end-1]
    normal = normalize(cross(chord,tipSpan))
    velocity = UVLM.freestream_velocity(system.freestream[])
    return atand(dot(velocity,normal),dot(velocity,chord))
end

function solve_swept_native_point(airspeed,pitchAngle,sweepAngle,x0)
    model,nElem,span = create_swept_pazy_model(airspeed,pitchAngle,sweepAngle)
    problem = AeroBeams.create_SteadyProblem(; model,x0)
    AeroBeams.solve!(problem)
    problem.systemSolver.loadStepState.convergedFinalSolution ||
        error("Native Pazy solve failed at sweep $(rad2deg(sweepAngle)) deg, pitch $(rad2deg(pitchAngle)) deg, $airspeed m/s")
    verify_pazy_vlm_def(model,airspeed,pitchAngle,sweepAngle)
    tip = model.beams[1].elements[end].state.nodal
    aoa = problem.loadstepHistory.aeroVariables[end][nElem].unsteadyKin.αₑ
    return (; oop=-100 * tip.u2[1] / span,twist=AeroBeams.twist_angle(tip.p2_b),
        strain=swept_te_strain(model),aoa=rad2deg(aoa),x=copy(problem.systemArrays.x))
end

function run_pazy_swept_steady_pitch_range(; outputDirectory=SWEPT_STUDY_DIRECTORY,
    sweepAngles=deg2rad.(collect(SWEPT_ANGLES_DEG)),pitches=SWEPT_PITCH_RANGE,
    airspeeds=SWEPT_AIRSPEEDS,uvlmOptions=SWEPT_UVLM_OPTIONS)

    @assert !isempty(sweepAngles) && all(isfinite,sweepAngles) "Use finite sweep angles"
    @assert all(angle -> any(supported -> isapprox(rad2deg(angle),supported; atol=1e-8),
        SWEPT_ANGLES_DEG),sweepAngles) "Pazy swept wings support 10° and 20°"
    @assert !isempty(airspeeds) && all(isfinite,airspeeds) && first(airspeeds) > 0 &&
        all(>(0),diff(airspeeds)) "Use finite, positive, increasing airspeeds"

    atmosphere = AeroBeams.standard_atmosphere(0)
    density = atmosphere.ρ
    mkpath(outputDirectory)
    columns = ["sweep_deg","pitch_deg","airspeed_m_s","dynamic_pressure_pa",
        "density_kg_m3","native_oop_pct","native_twist_deg","uvlm_oop_pct",
        "uvlm_twist_deg","uvlm_converged","uvlm_static_iterations","uvlm_wake_steps",
        "uvlm_geometry_residual","uvlm_load_residual","uvlm_wake_force_residual",
        "uvlm_wake_circulation_residual","uvlm_seconds","modeled_pitch_deg",
        "native_te_strain_ue","uvlm_te_strain_ue","native_tip_aoa_deg","uvlm_tip_aoa_deg"]
    failures = 0
    open(joinpath(outputDirectory,"steady_pitch_range.csv"),"w") do io
        writedlm(io,permutedims(columns),',')
        flush(io)
        for sweep_angle in sweepAngles
            sweep_deg = rad2deg(sweep_angle)
            _,span,chord,_,_ = AeroBeams.geometrical_properties_Pazy(; mesh=:nominal,
                Λ=sweep_angle)
            for pitch in pitches
                angle = deg2rad(pitch + (sweep_deg ≈ 20 ? 0.4 : 0.0))
                native_guess = Float64[]
                previous = nothing
                previous_speed = 0.0
                pazy = create_swept_pazy_uvlm_structure(sweep_angle)
                for airspeed in airspeeds
                    println("Pazy sweep $sweep_deg deg, pitch $pitch deg, airspeed $airspeed m/s")
                    flush(stdout)
                    native = solve_swept_native_point(airspeed,angle,sweep_angle,native_guess)
                    native_guess = native.x

                    target = isnothing(previous) ? nothing :
                        previous.loads .* (airspeed / previous_speed)^2
                    time_step_chords = get(uvlmOptions,:time_step_chords,
                        1 / uvlmOptions.chordwise_panels)
                    dt = time_step_chords * chord / airspeed
                    elapsed = @elapsed uvlm = run_pazy_wing_uvlm_static(;
                        uvlmOptions...,density,airspeed,angle_of_attack=angle,
                        sweepAngle=sweep_angle,pazy_structure=pazy,warm_start=previous,
                        warm_start_target_loads=target,maximum_aerodynamic_time=1500dt,
                        minimum_aerodynamic_time=0.0,aerodynamic_window_time=8dt)
                    # Map from the unpitched UVLM frame into the wind tunnel frame.
                    tipDisplacement = uvlm.tip_displacement_A
                    uvlm_oop = uvlm.converged ?
                        100 * (-cos(angle) * tipDisplacement[1] + sin(angle) * tipDisplacement[2]) / span : NaN
                    uvlm_twist = uvlm.converged ? uvlm.tip_twist_degrees : NaN
                    uvlm_strain = uvlm.converged ? swept_te_strain(uvlm.structural_problem.model) : NaN
                    uvlm_aoa = uvlm.converged ? swept_uvlm_tip_aoa(uvlm.aerodynamic_problem) : NaN
                    writedlm(io,permutedims([sweep_deg,pitch,airspeed,
                        0.5 * density * airspeed^2,density,native.oop,native.twist,
                        uvlm_oop,uvlm_twist,Int(uvlm.converged),size(uvlm.history,1),
                        uvlm.total_aerodynamic_steps,uvlm.geometry_residual,uvlm.load_residual,
                        uvlm.history[end,12],uvlm.history[end,13],elapsed,rad2deg(angle),
                        native.strain,uvlm_strain,native.aoa,uvlm_aoa]),',')
                    flush(io)
                    if uvlm.converged
                        previous = merge(uvlm.warm_start,(aero=nothing,))
                        previous_speed = airspeed
                    else
                        failures += 1
                        @warn uvlm.termination_reason sweep_deg pitch airspeed
                        pazy = create_swept_pazy_uvlm_structure(sweep_angle)
                        previous = nothing
                    end
                end
            end
        end
    end
    return failures
end

function plot_pazy_swept_steady_pitch_range(outputDirectory=SWEPT_STUDY_DIRECTORY)
    values,_ = readdlm(joinpath(outputDirectory,"steady_pitch_range.csv"),',',Float64;
        header=true)
    Plots.gr()
    common = (; framestyle=:axes,grid=true,gridalpha=0.15,guidefontsize=11,
        tickfontsize=9,left_margin=5Plots.mm,right_margin=3Plots.mm,
        bottom_margin=2Plots.mm,xlabel="Airspeed [m/s]",xlims=(0,70),xticks=0:10:70)
    colors = ("#800080","#0000ff","#00566d","#6c8c00","#ff8500","#ff0000")
    labels = ("Tip OOP disp. [% semispan]","Strains @ 40% span (με)","Tip angle of attack [deg]")
    captionLabels = ("tip OOP displacement.","strains at spar TE (approximately 40% span).",
        "tip angle of attack.")
    nativeColumns = (6,19,21)
    uvlmColumns = (8,20,22)
    panels = Plots.Plot[]
    for response in eachindex(labels)
        captions = Plots.Plot[]
        for (column,sweep) in enumerate(SWEPT_ANGLES_DEG)
            rows = values[isapprox.(values[:,1],sweep; atol=1e-8),:]
            panel = Plots.plot(; common...,ylabel=labels[response],
                legend=response == 1 ? :topleft : false,legendfontsize=8)
            Plots.plot!(panel,[NaN],[NaN]; color=:black,lw=1.3,label="AeroBeams (VLM-def)")
            Plots.plot!(panel,[NaN],[NaN]; color=:black,lw=1.7,ls=:dash,label="AeroBeams + UVLM")
            Plots.scatter!(panel,[NaN],[NaN]; color=:black,ms=3,msw=0,label="Exp.")
            for (pitch,color) in zip(SWEPT_PITCH_RANGE,colors)
                branch = rows[rows[:,2] .== pitch,:]
                Plots.plot!(panel,branch[:,3],branch[:,nativeColumns[response]];
                    color,lw=1.3,label="θ = $(pitch)°")
                Plots.plot!(panel,branch[:,3],branch[:,uvlmColumns[response]];
                    color,lw=1.7,ls=:dash,label=false)
                # Repository strain series refer to the root gauge. Use the
                # corresponding measurements only for displacement and tip AoA.
                if response != 2
                    quantity = response == 1 ? "tipdisp" : "tipaoa"
                    file = joinpath(pkgdir(AeroBeams),"test","referenceData","sweptPazyS$sweep",
                        "steady_$(quantity)_aoa$pitch.txt")
                    experimental = readdlm(file)
                    scale = response == 1 ? 100.0 : 1.0
                    Plots.scatter!(panel,experimental[1,:],scale .* experimental[2,:];
                        color,ms=2.6,msw=0,label=false)
                end
            end
            push!(panels,panel)
            letter = 'a' + (2 * (response - 1) + column - 1)
            caption = Plots.plot(xlims=(0,1),ylims=(0,1),axis=false,ticks=false,grid=false,
                framestyle=:none,margin=0Plots.mm,legend=false)
            Plots.annotate!(caption,0.5,0.5,
                Plots.text("($letter) S$sweep, $(captionLabels[response])",10,"Times Bold"))
            push!(captions,caption)
        end
        append!(panels,captions)
    end
    layout = Plots.grid(6,2; heights=[0.29,0.043,0.29,0.043,0.29,0.044])
    figure = Plots.plot(panels...; layout,size=(1080,1030),dpi=200)
    stem = joinpath(outputDirectory,"pazy_s10_s20_steady_pitch_range")
    Plots.savefig(figure,stem * ".png")
    Plots.savefig(figure,stem * ".pdf")
    return figure
end
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    if ARGS == ["--plot-only"]
        plot_pazy_swept_steady_pitch_range()
    elseif isempty(ARGS)
        LinearAlgebra.BLAS.set_num_threads(1)
        failures = run_pazy_swept_steady_pitch_range()
        plot_pazy_swept_steady_pitch_range()
        iszero(failures) || error("$failures UVLM conditions did not converge")
    else
        error("Usage: run_pazy_swept_steady_pitch_range.jl [--plot-only]")
    end
end
