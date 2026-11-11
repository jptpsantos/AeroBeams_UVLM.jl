# Compare two nonlinear, time-domain simulations of the same Pazy wing:
#
#   1. AeroBeams structure + WingPropellerUVLM aerodynamics
#   2. AeroBeams structure + AeroBeams Indicial strip aerodynamics
#
# Both simulations use the same structural model, time vector, speed ramp and
# tip-force pulse.  The root incidence is applied once: through the UVLM
# operating point in case 1 and through the Pazy root pitch in case 2.

include(joinpath(@__DIR__,"PazyWingUVLMCoupling.jl"))
using DelimitedFiles

# ---------------------------------------------------------------------------
# Controls: edit this section, then run this file.
# ---------------------------------------------------------------------------

# Flight condition
airspeed = 55.0                       # target airspeed [m/s]
density = 1.225                       # air density [kg/m^3]
angle_of_attack = deg2rad(2.0)        # root incidence [rad]
sideslip = 0.0                        # strip comparison requires zero sideslip

# Identical smooth startup for both aerodynamic models
initial_airspeed_fraction = 0.05      # positive because UVLM cannot start at U=0
airspeed_ramp_duration = 0.5          # [s]
duration = 3.0                        # [s]
settling_time = 1.0                   # tip pulse begins here [s]
perturbation_amplitude = 0.01         # tip F1A pulse amplitude [N]
perturbation_duration = 0.04          # [s]
equilibrium_average_duration = 0.20   # pre-pulse interval used for mean removal [s]

# UVLM discretization and wake
chordwise_panels = 4
spanwise_panels = 15
symmetric_wing = true
wake_length_chords = 10.0
maximum_wake_rows = ceil(Int,wake_length_chords * chordwise_panels)
core_radius = 1e-3                    # [m]
wake_shedding_fraction = 0.1
aerodynamic_interaction = true

# The comparison only needs the wingtip histories, so UVLM field histories are
# disabled.  Set this to true if circulation histories are also required.
save_uvlm_history = false
uvlm_save_frequency = 10

# AeroBeams Indicial strip-theory controls
strip_has_tip_correction = true
strip_tip_loss_type = "Exponential"  # alternatives: "VLM-undef" or "VLM-def"

# Newton-Raphson controls used by both structural simulations
newton_maximum_iterations = 50
newton_absolute_tolerance = 1e-8
newton_relative_tolerance = 1e-8
newton_display_iterations = false
newton_always_update_jacobian = false
progress_frequency = 100

# Strong coupling controls for the UVLM branch of the comparison
coupling_maximum_iterations = 20
coupling_relaxation = 0.3
coupling_geometry_tolerance = 1e-5
coupling_load_tolerance = 1e-3
coupling_display_iterations = false

output_directory = joinpath(@__DIR__,"output","aero_comparison")

# Numerical agreement between the two curves is a model comparison, not a
# UVLM mesh-convergence proof.  For convergence, repeat this file with finer
# (chordwise_panels, spanwise_panels), for example (4,15), (6,20), (8,30).

# ---------------------------------------------------------------------------
# Analysis
# ---------------------------------------------------------------------------

"""Smooth speed law used in both time-domain simulations."""
function comparison_airspeed(t,target,initialFraction,ramp_duration)
    ramp_duration == 0 && return target
    fraction = clamp(t / ramp_duration,0.0,1.0)
    smoothstep = fraction^2 * (3 - 2 * fraction)
    return target * (initialFraction + (1 - initialFraction) * smoothstep)
end

"""The same compact, smooth tip pulse used by the UVLM coupling."""
function comparison_tip_force(t,start_time,pulse_duration,amplitude)
    phase = (t - start_time) / pulse_duration
    return 0 < phase < 1 ? amplitude * sin(2pi * phase) * sinpi(phase)^2 : 0.0
end

function comparison_newton_solver()
    return AeroBeams.create_NewtonRaphson(
        maximumIterations=newton_maximum_iterations,
        absoluteTolerance=newton_absolute_tolerance,
        relativeTolerance=newton_relative_tolerance,
        displayStatus=newton_display_iterations,
        alwaysUpdateJacobian=newton_always_update_jacobian,
        minConvRateAeroJacUpdate=1.2,
        minConvRateJacUpdate=1.2
    )
end

"""Run the nonlinear AeroBeams model with its Indicial strip aerodynamics."""
function run_strip_comparison(time_vector)
    n_elements,span,_,_ = AeroBeams.geometrical_properties_Pazy()

    # create_Pazy replaces this dummy beam with the actual Pazy beam.
    dummyBeam = AeroBeams.create_Beam(
        length=span,
        nElements=n_elements,
        S=[AeroBeams.isotropic_stiffness_matrix()]
    )
    tipPulse = AeroBeams.create_BC(
        name="comparison tip pulse",
        beam=dummyBeam,
        node=n_elements + 1,
        types=["F1A"],
        values=[t -> comparison_tip_force(
            t,settling_time,perturbation_duration,perturbation_amplitude)]
    )
    speed = t -> comparison_airspeed(
        t,airspeed,initial_airspeed_fraction,airspeed_ramp_duration)

    model,_,_,_,_ = AeroBeams.create_Pazy(
        aeroSolver=AeroBeams.Indicial(),
        upright=true,
        θ=angle_of_attack,
        airspeed=speed,
        g=0.0,
        withSkin=true,
        sweepStructuralCorrections=false,
        GAy=1e16,
        GAz=1e16,
        hasTipCorrection=strip_has_tip_correction,
        tipLossType=strip_tip_loss_type,
        tipLossFunctionIsAirspeedDependent=true,
        additionalBCs=[tipPulse]
    )
    model.atmosphere.ρ = density

    problem = AeroBeams.create_DynamicProblem(
        model=model,
        timeVector=time_vector,
        systemSolver=comparison_newton_solver(),
        trackingTimeSteps=true,
        trackingFrequency=1,
        displayProgress=true,
        displayFrequency=progress_frequency,
        saveInitialSolution=true
    )
    AeroBeams.solve!(problem)
    length(problem.savedTimeVector) == length(time_vector) ||
        error("AeroBeams strip solution stopped before the requested final time")

    tipBending = [-states[n_elements].u_n2[1]
                  for states in problem.nodalStatesOverTime]
    tip_twist = [
        begin
            rotation,_ = AeroBeams.rotation_tensor_WM(states[n_elements].p_n2_b)
            chordDirection = rotation * [0.0; 1.0; 0.0]
            asind(clamp(chordDirection[3],-1.0,1.0))
        end for states in problem.nodalStatesOverTime
    ]

    return (
        problem=problem,
        time=copy(problem.savedTimeVector),
        tip_bending_displacement=tipBending,
        tip_twist_degrees=tip_twist
    )
end

mean_of(values) = sum(values) / length(values)
rms_of(values) = sqrt(sum(abs2,values) / length(values))

function run_aero_comparison()
    @assert airspeed > 0 && density > 0 && duration > 0
    @assert iszero(sideslip) "AeroBeams' Pazy strip model in this comparison has no sideslip input"
    @assert 0 < initial_airspeed_fraction <= 1
    @assert 0 <= airspeed_ramp_duration <= settling_time < duration
    @assert 0 < perturbation_duration && settling_time + perturbation_duration < duration
    @assert 0 < equilibrium_average_duration <= settling_time
    @assert chordwise_panels > 0 && spanwise_panels > 0
    @assert maximum_wake_rows > 0 && progress_frequency > 0

    println("\n=== 1/2: UVLM + nonlinear AeroBeams structure ===")
    uvlm = run_pazy_wing_uvlm(;
        airspeed,density,angle_of_attack,sideslip,
        initial_airspeed_fraction,airspeed_ramp_duration,
        duration,settling_time,
        chordwise_panels,spanwise_panels,symmetric_wing,
        maximum_wake_rows,core_radius,wake_shedding_fraction,
        aerodynamic_interaction,save_uvlm_history,uvlm_save_frequency,
        newton_maximum_iterations,newton_absolute_tolerance,
        newton_relative_tolerance,newton_display_iterations,
        newton_always_update_jacobian,
        perturbation_amplitude,perturbation_duration,
        animation_frames=2,
        progress_frequency,
        coupling_maximum_iterations,
        coupling_relaxation,
        coupling_geometry_tolerance,
        coupling_load_tolerance,
        coupling_display_iterations
    )

    println("\n=== 2/2: AeroBeams Indicial strip aerodynamics ===")
    strip = run_strip_comparison(uvlm.time)
    @assert strip.time ≈ uvlm.time

    # Remove the mean immediately before the pulse to compare perturbations
    # around each aerodynamic model's own deformed operating condition.
    equilibriumIndices = findall(t ->
            settling_time - equilibrium_average_duration <= t <= settling_time,
        uvlm.time)
    isempty(equilibriumIndices) && error("No samples in equilibrium averaging interval")

    uvlmBendingMean = mean_of(uvlm.tip_bending_displacement[equilibriumIndices])
    stripBendingMean = mean_of(strip.tip_bending_displacement[equilibriumIndices])
    uvlmTwistMean = mean_of(uvlm.tip_twist_degrees[equilibriumIndices])
    stripTwistMean = mean_of(strip.tip_twist_degrees[equilibriumIndices])

    uvlmBendingDelta = uvlm.tip_bending_displacement .- uvlmBendingMean
    stripBendingDelta = strip.tip_bending_displacement .- stripBendingMean
    uvlmTwistDelta = uvlm.tip_twist_degrees .- uvlmTwistMean
    stripTwistDelta = strip.tip_twist_degrees .- stripTwistMean

    responseIndices = findall(t ->
            t >= settling_time + perturbation_duration,uvlm.time)
    bendingRMSError = rms_of(
        uvlmBendingDelta[responseIndices] .- stripBendingDelta[responseIndices])
    twistRMSError = rms_of(
        uvlmTwistDelta[responseIndices] .- stripTwistDelta[responseIndices])

    mkpath(output_directory)
    data_path = joinpath(output_directory,"pazy_aero_comparison.csv")
    data = hcat(
        uvlm.time,
        uvlm.airspeed_history,
        uvlm.tip_bending_displacement,
        strip.tip_bending_displacement,
        uvlm.tip_twist_degrees,
        strip.tip_twist_degrees,
        uvlmBendingDelta,
        stripBendingDelta,
        uvlmTwistDelta,
        stripTwistDelta
    )
    header = permutedims([
        "time_s","airspeed_m_per_s","uvlm_bending_m","strip_bending_m",
        "uvlm_twist_deg","strip_twist_deg","uvlm_delta_bending_m",
        "strip_delta_bending_m","uvlm_delta_twist_deg","strip_delta_twist_deg"
    ])
    writedlm(data_path,vcat(header,data),',')

    plots = AeroBeams.Plots
    p1 = plots.plot(
        uvlm.time,uvlm.tip_bending_displacement,
        label="UVLM",color=:blue,linewidth=2,
        ylabel="Tip bending [m]",title="Absolute histories"
    )
    plots.plot!(p1,strip.time,strip.tip_bending_displacement,
        label="Indicial strip",color=:orange,linewidth=2,linestyle=:dash)
    p2 = plots.plot(
        uvlm.time,uvlm.tip_twist_degrees,
        label="UVLM",color=:blue,linewidth=2,
        xlabel="Time [s]",ylabel="Tip twist [deg]"
    )
    plots.plot!(p2,strip.time,strip.tip_twist_degrees,
        label="Indicial strip",color=:orange,linewidth=2,linestyle=:dash)

    responseTime = uvlm.time[responseIndices] .-
                   (settling_time + perturbation_duration)
    p3 = plots.plot(
        responseTime,uvlmBendingDelta[responseIndices],
        label="UVLM",color=:blue,linewidth=2,
        ylabel="Δ bending [m]",title="Post-pulse, equilibrium removed"
    )
    plots.plot!(p3,responseTime,stripBendingDelta[responseIndices],
        label="Indicial strip",color=:orange,linewidth=2,linestyle=:dash)
    p4 = plots.plot(
        responseTime,uvlmTwistDelta[responseIndices],
        label="UVLM",color=:blue,linewidth=2,
        xlabel="Time after pulse [s]",ylabel="Δ twist [deg]"
    )
    plots.plot!(p4,responseTime,stripTwistDelta[responseIndices],
        label="Indicial strip",color=:orange,linewidth=2,linestyle=:dash)

    figure = plots.plot(p1,p2,p3,p4;
        layout=(2,2),size=(1200,750),margin=5 * AeroBeams.Measures.mm)
    figurePath = joinpath(output_directory,"pazy_aero_comparison.png")
    plots.savefig(figure,figurePath)

    println("\n=== Comparison complete ===")
    println("UVLM equilibrium bending [m]: ",uvlmBendingMean)
    println("Strip equilibrium bending [m]: ",stripBendingMean)
    println("UVLM equilibrium twist [deg]: ",uvlmTwistMean)
    println("Strip equilibrium twist [deg]: ",stripTwistMean)
    println("Post-pulse bending RMS difference [m]: ",bendingRMSError)
    println("Post-pulse twist RMS difference [deg]: ",twistRMSError)
    println("Plot: ",figurePath)
    println("Data: ",data_path)

    return (
        uvlm=uvlm,
        strip=strip,
        equilibrium=(
            uvlm_bending=uvlmBendingMean,
            strip_bending=stripBendingMean,
            uvlm_twist=uvlmTwistMean,
            strip_twist=stripTwistMean
        ),
        rms_difference=(bending=bendingRMSError,twist=twistRMSError),
        plot_path=figurePath,
        data_path=data_path
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    comparison = run_aero_comparison()
end
