# Run zero-torque windmilling and 650 N thrusting trims for the Chang propeller.

import Pkg
const EXAMPLE_DIR = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(normpath(joinpath(EXAMPLE_DIR, "..", "..")))

using DelimitedFiles
using Printf
using WingPropellerUVLM

include(joinpath(EXAMPLE_DIR, "src", "ChangAeroelastic.jl"))
using .ChangAeroelastic: load_chang_configuration
config = load_chang_configuration()
include(joinpath(EXAMPLE_DIR, "src", "chang_propeller_trim.jl"))

env_float(name, default) = parse(Float64, get(ENV, name, string(default)))
env_int(name, default) = parse(Int, get(ENV, name, string(default)))
env_symbol(name, default) = Symbol(lowercase(strip(get(ENV, name, string(default)))))

options = ChangWindmillingTrimOptions(
    flow_speed_mps = env_float("CHANG_TRIM_SPEED_MPS", 65.0),
    air_density_kgpm3 = env_float("CHANG_TRIM_AIR_DENSITY", 1.225),
    angle_of_attack_deg = env_float("CHANG_TRIM_ALPHA_DEG", 3.0),
    sideslip_deg = env_float("CHANG_TRIM_BETA_DEG", config.simulation.sideslip_deg),
    propeller_radius_m = config.propeller.radius_m,
    propeller_chord_m = config.propeller.chord_m,
    blade_count = config.propeller.blades,
    radial_panels = env_int("CHANG_TRIM_RADIAL_PANELS", 10),
    chordwise_panels = env_int("CHANG_TRIM_CHORDWISE_PANELS", 10),
    collective_pitch_offset_deg = env_float("CHANG_TRIM_COLLECTIVE_OFFSET_DEG", 0.0),
    azimuth_step_deg = env_float("CHANG_TRIM_AZIMUTH_STEP_DEG", 2.5),
    simulated_revolutions = env_int("CHANG_TRIM_SIMULATED_REVOLUTIONS", 4),
    averaged_revolutions = env_int("CHANG_TRIM_AVERAGED_REVOLUTIONS", 1),
    retained_wake_revolutions = env_int("CHANG_TRIM_RETAINED_WAKE_REVOLUTIONS", 1),
    wake_relaxation = env_float("CHANG_TRIM_WAKE_RELAXATION", 0.1),
    core_radius_m = ChangAeroelastic.environment_optional_number(
        Float64, "CHANG_TRIM_CORE_RADIUS_M", 0.001,
    ),
    vortex_core_span_fraction = env_float("CHANG_TRIM_FCORE_SPAN_FACTOR", 0.0),
    vortex_core_chord_fraction = env_float("CHANG_TRIM_FCORE_CHORD_FACTOR", 0.001),
    near_field_force_model = env_symbol(
        "CHANG_TRIM_NEAR_FIELD_FORCE_MODEL",
        config.simulation.near_field_force_model,
    ),
    initial_rpm = env_float("CHANG_TRIM_INITIAL_RPM", config.propeller.rotation_rpm),
    rpm_bracket = (
        env_float("CHANG_TRIM_MIN_RPM", 0.85 * config.propeller.rotation_rpm),
        env_float("CHANG_TRIM_MAX_RPM", 1.15 * config.propeller.rotation_rpm),
    ),
    torque_tolerance_nm = env_float("CHANG_TRIM_TORQUE_TOLERANCE_NM", 0.25),
    thrust_tolerance_n = env_float("CHANG_TRIM_THRUST_TOLERANCE_N", 1.0),
    rpm_tolerance = env_float("CHANG_TRIM_RPM_TOLERANCE", 0.25),
    maximum_root_iterations = env_int("CHANG_TRIM_MAX_ITERATIONS", 8),
)

target_thrust_n = env_float("CHANG_THRUST_TARGET_N", 650.0)
output_directory = normpath(get(
    ENV,
    "CHANG_TRIM_OUTPUT_DIR",
    joinpath(EXAMPLE_DIR, "output", "chang_propeller_trims"),
))
mkpath(output_directory)

function write_trim_outputs(output_directory, label, solution; target_thrust_n = nothing)
    trim = solution.trim
    evaluation_path = joinpath(output_directory, "$(label)_evaluations.csv")
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
            result.periodic_thrust_change_n,
            result.torque_standard_deviation_nm,
            result.thrust_standard_deviation_n,
        ]) for result in solution.evaluations
    ]
    evaluation_matrix = reduce(vcat, evaluation_rows)
    open(evaluation_path, "w") do io
        println(io, "rpm,mean_torque_nm,mean_thrust_n,CQ,CT,J,mu,"
            * "last_revolution_torque_change_nm,last_revolution_thrust_change_n,"
            * "torque_std_nm,thrust_std_n")
        writedlm(io, evaluation_matrix, ',')
    end

    history_path = joinpath(output_directory, "$(label)_history.csv")
    azimuth_degrees = (collect(1:length(trim.torque_history_nm)) .*
        options.azimuth_step_deg) .% 360
    open(history_path, "w") do io
        println(io, "step,azimuth_deg,torque_nm,thrust_n")
        writedlm(io, hcat(
            collect(1:length(trim.torque_history_nm)),
            azimuth_degrees,
            trim.torque_history_nm,
            trim.thrust_history_n,
        ), ',')
    end

    summary_path = joinpath(output_directory, "$(label)_summary.txt")
    open(summary_path, "w") do io
        println(io, "Chang isolated-propeller trim: $label")
        println(io, "Near-field force model: $(options.near_field_force_model)")
        println(io, "Finite-core radius: $(chang_trim_core_description(options))")
        @printf(io, "Flow speed: %.8f m/s\n", options.flow_speed_mps)
        @printf(io, "Angle of attack: %.8f deg\n", options.angle_of_attack_deg)
        @printf(io, "Blade angle at 75%% radius: %.8f deg\n", chang_blade_angle_75_deg(options))
        @printf(io, "Grid: %d radial x %d chordwise panels per blade\n",
            options.radial_panels, options.chordwise_panels)
        @printf(io, "RPM: %.8f\n", trim.rpm)
        @printf(io, "Mean torque: %+.10f N m\n", trim.mean_torque_nm)
        @printf(io, "Mean thrust: %+.10f N\n", trim.mean_thrust_n)
        if isnothing(target_thrust_n)
            @printf(io, "Torque residual: %+.10f N m\n", trim.mean_torque_nm)
        else
            @printf(io, "Target thrust: %.10f N\n", target_thrust_n)
            @printf(io, "Thrust residual: %+.10f N\n", trim.mean_thrust_n - target_thrust_n)
        end
        @printf(io, "J: %.10f\n", trim.advance_ratio)
        @printf(io, "mu=J/pi: %.10f\n", trim.inflow_ratio)
        @printf(io, "CQ: %+.12e\n", trim.torque_coefficient)
        @printf(io, "CT: %+.12e\n", trim.thrust_coefficient)
    end
    return (; evaluation_path, history_path, summary_path)
end

println("Chang propeller trims")
@printf("  V=%.3f m/s, alpha=%.3f deg, target thrust=%.3f N\n",
    options.flow_speed_mps, options.angle_of_attack_deg, target_thrust_n)
println("  Running windmilling trim (Q=0)...")
windmilling_solution = trim_chang_windmilling_rpm(options)
println("  Running thrusting trim (T=$(target_thrust_n) N)...")
thrusting_solution = trim_chang_thrusting_rpm(
    options;
    target_thrust_n = target_thrust_n,
)

windmilling_paths = write_trim_outputs(output_directory, "windmilling", windmilling_solution)
thrusting_paths = write_trim_outputs(
    output_directory,
    "thrusting_$(Int(round(target_thrust_n)))N",
    thrusting_solution;
    target_thrust_n = target_thrust_n,
)

summary_path = joinpath(output_directory, "trim_comparison_summary.txt")
open(summary_path, "w") do io
    println(io, "Chang propeller trim comparison")
    @printf(io, "Windmilling: RPM=%.8f, Q=%+.10f N m, T=%+.10f N\n",
        windmilling_solution.trim.rpm,
        windmilling_solution.trim.mean_torque_nm,
        windmilling_solution.trim.mean_thrust_n)
    @printf(io, "Thrusting: target T=%.8f N, RPM=%.8f, Q=%+.10f N m, T=%+.10f N\n",
        target_thrust_n,
        thrusting_solution.trim.rpm,
        thrusting_solution.trim.mean_torque_nm,
        thrusting_solution.trim.mean_thrust_n)
end

println("\nTrim results")
@printf("  windmilling: RPM=%.4f, Q=%+.6f N m, T=%+.6f N\n",
    windmilling_solution.trim.rpm,
    windmilling_solution.trim.mean_torque_nm,
    windmilling_solution.trim.mean_thrust_n)
@printf("  thrusting:   RPM=%.4f, Q=%+.6f N m, T=%+.6f N\n",
    thrusting_solution.trim.rpm,
    thrusting_solution.trim.mean_torque_nm,
    thrusting_solution.trim.mean_thrust_n)
println("  comparison:   $summary_path")
println("  windmilling:  $(windmilling_paths.summary_path)")
println("  thrusting:    $(thrusting_paths.summary_path)")