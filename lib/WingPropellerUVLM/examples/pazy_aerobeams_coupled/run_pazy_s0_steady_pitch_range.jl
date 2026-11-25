# Steady S0 beam-rod benchmark, with native strip aerodynamics and UVLM.
ENV["GKSwstype"] = get(ENV,"GKSwstype","100")
include(joinpath(@__DIR__,"PazyWingUVLMStatic.jl"))
using DelimitedFiles
import Plots

const S0_PITCH_RANGE = (3,5,7)
const S0_STUDY_DIRECTORY = joinpath(@__DIR__,"output","s0_steady_pitch_range")
const S0_MODEL_OPTIONS = (;
    mesh=:nominal,beamWingtipRod=true,nElementsWingtipRod=12,upright=true,
    sweepAngle=0.0,withWingtipRod=true,rootTorsionCorr=true,sweptTipMassConfig=:baseline)
const S0_UVLM_OPTIONS = (;
    chordwise_panels=4,spanwise_panels=15,symmetric_wing=true,maximum_wake_rows=40,
    core_radius=1e-6,wake_shedding_fraction=0.1,aerodynamic_interaction=true,
    aerodynamic_force_tolerance=1e-4,aerodynamic_circulation_tolerance=1e-4,
    maximum_static_iterations=100,static_relaxation=1.0,geometry_tolerance=1e-5,
    load_tolerance=1e-4,load_residual_force_floor=1e-6,
    consecutive_equilibrium_iterations=2,newton_maximum_iterations=50,
    newton_absolute_tolerance=1e-8,newton_relative_tolerance=1e-8,
    newton_display_iterations=false,progress_frequency=1000,
    relaxation_mode=:anderson,anderson_depth=4,geometry_cap_chords=2.0,wake_reuse=true)

# Avin et al. test 1 used a 10 g balance mass 40 mm behind the trailing edge.
function create_s0_benchmark_model(airspeed,pitchAngle; hasTipCorrection=true)
    _,_,chord,sparFraction = AeroBeams.geometrical_properties_Pazy(mesh=S0_MODEL_OPTIONS.mesh)
    tipOffset = [0.0,-chord * (1 - sparFraction) - 0.04,0.0]
    return AeroBeams.create_Pazy(; S0_MODEL_OPTIONS...,pitchAngle,airspeed,
        hasTipCorrection,tipLossType="Exponential",aeroSolver=AeroBeams.Indicial(),
        airfoil=deepcopy(AeroBeams.FlatPlate),
        tipMasses=[0.01],ηtipMasses=[tipOffset],tipMassesInertias=[zeros(3,3)])
end

function run_s0_steady_pitch_range(; outputDirectory=S0_STUDY_DIRECTORY,
    pitches=S0_PITCH_RANGE,airspeeds=[1e-3; collect(1.0:60.0)],
    uvlmOptions=S0_UVLM_OPTIONS)

    @assert !isempty(airspeeds) && all(isfinite,airspeeds) && first(airspeeds) > 0 &&
        all(>(0),diff(airspeeds)) "Use finite, positive, increasing airspeeds"
    atmosphere = AeroBeams.standard_atmosphere(0)
    density = atmosphere.ρ
    _,semispan,chord = AeroBeams.geometrical_properties_Pazy(mesh=S0_MODEL_OPTIONS.mesh)
    mkpath(outputDirectory)
    columns = ["pitch_deg","airspeed_m_s","dynamic_pressure_pa","density_kg_m3",
        "native_oop_pct","native_twist_deg","uvlm_oop_pct","uvlm_twist_deg",
        "uvlm_converged","uvlm_static_iterations","uvlm_wake_steps",
        "uvlm_geometry_residual","uvlm_load_residual","uvlm_wake_force_residual",
        "uvlm_wake_circulation_residual","uvlm_seconds"]
    failures = 0
    open(joinpath(outputDirectory,"steady_pitch_range.csv"),"w") do io
        writedlm(io,permutedims(columns),',')
        flush(io)
        for pitch in pitches
            angle = deg2rad(pitch)
            nativeGuess = Float64[]
            previous = nothing
            previousSpeed = 0.0
            # In the upright S0 model gravity is parallel to the pitch axis.
            # Keep the root unpitched and rotate the UVLM freestream by angle.
            pazy = create_pazy_uvlm_structure(
                create_s0_benchmark_model(0.0,0.0; hasTipCorrection=false)...;
                parameters=(; S0_MODEL_OPTIONS...,test1_balance_mass_kg=0.01))
            for airspeed in airspeeds
                println("S0 pitch $pitch deg, airspeed $airspeed m/s")
                flush(stdout)
                model = first(create_s0_benchmark_model(airspeed,angle))
                native = AeroBeams.create_SteadyProblem(model=model,x0=nativeGuess)
                AeroBeams.solve!(native)
                native.systemSolver.loadStepState.convergedFinalSolution ||
                    error("Native S0 solve failed at $pitch deg, $airspeed m/s")
                nativeGuess = copy(native.systemArrays.x)
                tip = model.beams[1].elements[end].state.nodal
                nativeOOP = 100 * tip.u2_b[3] / semispan
                nativeTwist = AeroBeams.twist_angle(tip.p2_b)

                # Only the structural state is continued across airspeeds.
                # Each condition develops a new wake at its own speed and dt.
                target = isnothing(previous) ? nothing :
                    previous.loads .* (airspeed / previousSpeed)^2
                timeStepChords = get(uvlmOptions,:time_step_chords,1 / uvlmOptions.chordwise_panels)
                dt = timeStepChords * chord / airspeed
                elapsed = @elapsed uvlm = run_pazy_wing_uvlm_static(;
                    uvlmOptions...,density,airspeed,angle_of_attack=angle,pazy_structure=pazy,
                    warm_start=previous,warm_start_target_loads=target,
                    maximum_aerodynamic_time=1500dt,minimum_aerodynamic_time=0.0,
                    aerodynamic_window_time=8dt)
                # NaN breaks the curve at any failed equilibrium.
                uvlmOOP = uvlm.converged ? 100 * uvlm.tip_out_of_plane / semispan : NaN
                uvlmTwist = uvlm.converged ? uvlm.tip_twist_degrees : NaN
                writedlm(io,permutedims([pitch,airspeed,0.5 * density * airspeed^2,
                    density,nativeOOP,nativeTwist,uvlmOOP,uvlmTwist,Int(uvlm.converged),
                    size(uvlm.history,1),uvlm.total_aerodynamic_steps,uvlm.geometry_residual,
                    uvlm.load_residual,uvlm.history[end,12],uvlm.history[end,13],elapsed]),',')
                flush(io)
                if uvlm.converged
                    previous = merge(uvlm.warm_start,(aero=nothing,))
                    previousSpeed = airspeed
                else
                    failures += 1
                    @warn uvlm.termination_reason pitch airspeed
                    pazy = create_pazy_uvlm_structure(
                        create_s0_benchmark_model(0.0,0.0; hasTipCorrection=false)...;
                        parameters=pazy.parameters)
                    previous = nothing
                end
            end
        end
    end
    return failures
end

function plot_s0_steady_pitch_range(outputDirectory=S0_STUDY_DIRECTORY)
    values,_ = readdlm(joinpath(outputDirectory,"steady_pitch_range.csv"),',',Float64;
        header=true)
    _,semispan = AeroBeams.geometrical_properties_Pazy(mesh=S0_MODEL_OPTIONS.mesh)
    dataDirectory = joinpath(pkgdir(AeroBeams),"test","referenceData","Pazy")
    Plots.gr()
    common = (;
        framestyle=:axes,grid=true,gridalpha=0.15,guidefontsize=12,
        tickfontsize=9,left_margin=4Plots.mm,bottom_margin=2Plots.mm)
    maxOOP = maximum(filter(isfinite,vec(values[:,[5,7]])))
    tipOOP = Plots.plot(; common...,xlabel="Dynamic pressure [Pa]",
        ylabel="Tip OOP disp. [% semispan]",xlims=(0,2210),xticks=0:500:2000,
        ylims=(0,max(55.0,maxOOP + 2)),yticks=0:10:max(50,ceil(Int,maxOOP)),
        legend=:topleft,legendfontsize=8)
    finiteTwist = filter(isfinite,vec(values[:,[6,8]]))
    minTwist = min(-0.2,minimum(finiteTwist) - 0.1)
    maxTwist = max(2.9,maximum(finiteTwist) + 0.15)
    tipTwist = Plots.plot(; common...,xlabel="Airspeed [m/s]",ylabel="Tip twist [deg]",
        xlims=(0,60),xticks=0:10:60,ylims=(minTwist,maxTwist),
        yticks=ceil(2minTwist)/2:0.5:maxTwist,legend=false)

    Plots.plot!(tipOOP,[NaN],[NaN]; color=:black,lw=1.35,label="AeroBeams")
    Plots.plot!(tipOOP,[NaN],[NaN]; color=:black,lw=1.65,ls=:dash,label="AeroBeams + UVLM")
    Plots.scatter!(tipOOP,[NaN],[NaN]; color=:black,markershape=:diamond,
        ms=2.5,msw=0,label="Exp. - Avin et al. (2022)")
    for (pitch,color) in zip(S0_PITCH_RANGE,("#800080","#008000","#ff0000"))
        rows = values[values[:,1] .== pitch,:]
        Plots.plot!(tipOOP,rows[:,3],rows[:,5]; color,lw=1.35,label="αᵣ = $pitch°")
        Plots.plot!(tipOOP,rows[:,3],rows[:,7]; color,lw=1.65,ls=:dash,label=false)
        Plots.plot!(tipTwist,rows[:,2],rows[:,6]; color,lw=1.35,label=false)
        Plots.plot!(tipTwist,rows[:,2],rows[:,8]; color,lw=1.65,ls=:dash,label=false)

        oopData = readdlm(joinpath(dataDirectory,"tip_u3Vsq_rootPitch$(pitch)_ExpTest1.txt"))
        Plots.scatter!(tipOOP,oopData[1,:],100 .* oopData[2,:] ./ semispan;
            color,markershape=:diamond,ms=1.9,msw=0,label=false)
        # AeroBeams distributes the wind-tunnel-derived twist series at 5° and 7°.
        if pitch in (5,7)
            twistData = readdlm(joinpath(dataDirectory,"tip_thetaVsU_rootPitch$(pitch)_Exp.txt"))
            Plots.scatter!(tipTwist,twistData[1,:],twistData[2,:];
                color,markershape=:diamond,ms=1.9,msw=0,label=false)
        end
    end

    layout = Plots.@layout [grid(1,2){0.83h}; grid(1,2){0.08h}; caption]
    figure = Plots.plot(tipOOP,tipTwist,
        s0_caption("(a) Tip OOP displacement.",11),s0_caption("(b) Tip twist.",11),
        s0_caption("S0 beam-rod model: steady response over the pitch range.",12);
        layout,size=(1080,465),dpi=200)
    stem = joinpath(outputDirectory,"pazy_s0_steady_pitch_range")
    Plots.savefig(figure,stem * ".png")
    Plots.savefig(figure,stem * ".pdf")
    return figure
end

function s0_caption(label,size)
    panel = Plots.plot(xlims=(0,1),ylims=(0,1),axis=false,ticks=false,grid=false,
        framestyle=:none,margin=0Plots.mm,legend=false)
    Plots.annotate!(panel,0.5,0.5,Plots.text(label,size,"Times Bold",:center))
    return panel
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    if ARGS == ["--plot-only"]
        plot_s0_steady_pitch_range()
    elseif isempty(ARGS)
        LinearAlgebra.BLAS.set_num_threads(1)
        failures = run_s0_steady_pitch_range()
        plot_s0_steady_pitch_range()
        iszero(failures) || error("$failures UVLM conditions did not converge")
    else
        error("Usage: run_pazy_s0_steady_pitch_range.jl [--plot-only]")
    end
end
