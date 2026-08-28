# Imperial-load windmilling trim for the isolated rigid Chang propeller.

import Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..", "..")))

using DelimitedFiles
using Printf
using WingPropellerUVLM

include(joinpath(@__DIR__, "chang_case.jl"))
include(joinpath(@__DIR__, "chang_windmilling_trim.jl"))

env_float(name, default) = parse(Float64, get(ENV, name, string(default)))
env_int(name, default) = parse(Int, get(ENV, name, string(default)))

options = ChangWindmillingTrimOptions(
    flow_speed_mps = env_float("CHANG_TRIM_SPEED_MPS", PROPELLER_CONFIG.trim_speed_mps),
    air_density_kgpm3 = env_float("CHANG_TRIM_AIR_DENSITY", 1.225),
    angle_of_attack_deg = env_float("CHANG_TRIM_ALPHA_DEG", SIMULATION_CONFIG.angle_of_attack_deg),
    sideslip_deg = env_float("CHANG_TRIM_BETA_DEG", SIMULATION_CONFIG.sideslip_deg),
    propeller_radius_m = PROPELLER_CONFIG.radius_m,
    propeller_chord_m = PROPELLER_CONFIG.chord_m,
    blade_count = PROPELLER_CONFIG.blades,
    radial_panels = env_int("CHANG_TRIM_RADIAL_PANELS", PROPELLER_CONFIG.radial_panels),
    chordwise_panels = env_int("CHANG_TRIM_CHORDWISE_PANELS", PROPELLER_CONFIG.chordwise_panels),
    blade_pitch_75_deg = env_float("CHANG_TRIM_BETA75_DEG", 0.0),
    azimuth_step_deg = env_float("CHANG_TRIM_AZIMUTH_STEP_DEG", SIMULATION_CONFIG.azimuth_step_deg),
    simulated_revolutions = env_int("CHANG_TRIM_SIMULATED_REVOLUTIONS", 4),
    averaged_revolutions = env_int("CHANG_TRIM_AVERAGED_REVOLUTIONS", 1),
    retained_wake_revolutions = env_int("CHANG_TRIM_RETAINED_WAKE_REVOLUTIONS", 3),
    wake_relaxation = env_float("CHANG_TRIM_WAKE_RELAXATION", 0.1),
    vortex_core_span_fraction = env_float("CHANG_TRIM_FCORE_SPAN_FACTOR", 0.5),
    vortex_core_chord_fraction = env_float("CHANG_TRIM_FCORE_CHORD_FACTOR", 0.3),
    initial_rpm = env_float("CHANG_TRIM_INITIAL_RPM", PROPELLER_CONFIG.rotation_rpm),
    rpm_bracket = (
        env_float("CHANG_TRIM_MIN_RPM", 0.85 * PROPELLER_CONFIG.rotation_rpm),
        env_float("CHANG_TRIM_MAX_RPM", 1.15 * PROPELLER_CONFIG.rotation_rpm),
    ),
    torque_tolerance_nm = env_float("CHANG_TRIM_TORQUE_TOLERANCE_NM", 0.25),
    rpm_tolerance = env_float("CHANG_TRIM_RPM_TOLERANCE", 0.25),
    maximum_root_iterations = env_int("CHANG_TRIM_MAX_ITERATIONS", 8),
)

println("Imperial Chang windmilling trim")
@printf(
    "  V=%.3f m/s, alpha=%.3f deg, beta75=%.3f deg, grid=%dx%d, azimuth step=%.3f deg\n",
    options.flow_speed_mps,
    options.angle_of_attack_deg,
    options.blade_pitch_75_deg,
    options.radial_panels,
    options.chordwise_panels,
    options.azimuth_step_deg,
)
trim_solution = trim_chang_windmilling_rpm(options)
trim = trim_solution.trim

output_directory = normpath(get(
    ENV,
    "CHANG_TRIM_OUTPUT_DIR",
    joinpath(@__DIR__, "output", "chang_windmilling_trim"),
))
mkpath(output_directory)

evaluation_path = joinpath(output_directory, "imperial_windmilling_rpm_evaluations.csv")
evaluation_rows = [
    permutedims([
        result.rpm,
        result.mean_torque_nm,
        result.mean_thrust_n,
        result.torque_coefficient,
        result.thrust_coefficient,
        result.advance_ratio,
        result.inflow_ratio,
        result.periodic_torque_change_nm,
        result.torque_standard_deviation_nm,
    ])
    for result in trim_solution.evaluations
]
evaluation_matrix = reduce(vcat, evaluation_rows)
open(evaluation_path, "w") do io
    println(io, "rpm,mean_torque_nm,mean_thrust_n,CQ,CT,J,mu,last_revolution_torque_change_nm,torque_std_nm")
    writedlm(io, evaluation_matrix, ',')
end

history_path = joinpath(output_directory, "imperial_windmilling_trim_history.csv")
azimuth_degrees = (
    collect(1:length(trim.torque_history_nm)) .* options.azimuth_step_deg
) .% 360
open(history_path, "w") do io
    println(io, "step,azimuth_deg,torque_nm,thrust_n")
    writedlm(
        io,
        hcat(
            collect(1:length(trim.torque_history_nm)),
            azimuth_degrees,
            trim.torque_history_nm,
            trim.thrust_history_n,
        ),
        ',',
    )
end

summary_path = joinpath(output_directory, "imperial_windmilling_trim_summary.txt")
open(summary_path, "w") do io
    println(io, "Chang isolated-propeller windmilling trim (Imperial near-field loads)")
    @printf(io, "Flow speed: %.8f m/s\n", options.flow_speed_mps)
    @printf(io, "Air density: %.8f kg/m^3\n", options.air_density_kgpm3)
    @printf(io, "Angle of attack: %.8f deg\n", options.angle_of_attack_deg)
    @printf(io, "Blade pitch at 75%% radius: %.8f deg\n", options.blade_pitch_75_deg)
    @printf(io, "Grid: %d radial x %d chordwise panels per blade\n", options.radial_panels, options.chordwise_panels)
    @printf(io, "Azimuth step: %.8f deg\n", options.azimuth_step_deg)
    @printf(io, "Simulated/averaged/retained-wake revolutions: %d/%d/%d\n",
        options.simulated_revolutions, options.averaged_revolutions, options.retained_wake_revolutions)
    @printf(io, "Initial comparison RPM: %.8f\n", trim_solution.initial.rpm)
    @printf(io, "Imperial torque at initial RPM: %+.10f N m\n", trim_solution.initial.mean_torque_nm)
    @printf(io, "Imperial windmilling RPM: %.8f\n", trim.rpm)
    @printf(io, "Residual mean torque: %+.10f N m\n", trim.mean_torque_nm)
    @printf(io, "Mean thrust: %+.10f N\n", trim.mean_thrust_n)
    @printf(io, "Advance ratio J: %.10f\n", trim.advance_ratio)
    @printf(io, "Inflow ratio mu=J/pi: %.10f\n", trim.inflow_ratio)
    @printf(io, "CQ: %+.12e\n", trim.torque_coefficient)
    @printf(io, "CT: %+.12e\n", trim.thrust_coefficient)
    @printf(io, "Last-revolution torque standard deviation: %.10f N m\n", trim.torque_standard_deviation_nm)
    @printf(io, "Last two revolution mean-torque change: %+.10f N m\n", trim.periodic_torque_change_nm)
    println(io, "Mean torque by revolution (N m): $(trim.revolution_mean_torque_nm)")
end

println("\nTrim result")
@printf("  initial RPM:   %.4f (Imperial mean Q=%+.4f N m)\n",
    trim_solution.initial.rpm, trim_solution.initial.mean_torque_nm)
@printf("  Imperial RPM:  %.4f\n", trim.rpm)
@printf("  residual Q:    %+.6f N m\n", trim.mean_torque_nm)
@printf("  J:             %.8f\n", trim.advance_ratio)
@printf("  mu=J/pi:       %.8f\n", trim.inflow_ratio)
println("  summary:       $summary_path")
println("  evaluations:   $evaluation_path")
println("  trim history:  $history_path")
