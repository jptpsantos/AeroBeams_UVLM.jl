# Isolated-propeller trim solver shared by the windmilling and powered studies.

using LinearAlgebra
using Printf
using StaticArrays
using Statistics

using WingPropellerUVLM:
    Freestream,
    Reference,
    RotationMatrix,
    System,
    copy_surfaces_to_previous!,
    generate_propeller_blades_grid,
    get_nodal_properties_chang,
    get_twist_deg_interp,
    grid_to_surface_panels,
    imperial_nodal_forces,
    imperial_nodal_positions,
    legacy_imperial_segment_forces!,
    near_field_forces!,
    propagate_system!,
    repeated_trailing_edge_points

"""Numerical and physical controls for isolated Chang propeller trims."""
Base.@kwdef struct ChangWindmillingTrimOptions
    flow_speed_mps::Float64 = 65.0
    air_density_kgpm3::Float64 = 1.225
    angle_of_attack_deg::Float64 = 3.0
    sideslip_deg::Float64 = 0.0
    propeller_radius_m::Float64 = 1.15
    propeller_chord_m::Float64 = 0.197
    blade_count::Int = 4
    radial_panels::Int = 20
    chordwise_panels::Int = 20
    # Added uniformly to the complete Chang twist distribution. Leave at zero
    # to use the blade geometry defined in BladeGeometry.jl unchanged.
    collective_pitch_offset_deg::Float64 = 0.0
    azimuth_step_deg::Float64 = 5.0
    simulated_revolutions::Int = 4
    averaged_revolutions::Int = 1
    retained_wake_revolutions::Int = 3
    wake_relaxation::Float64 = 0.1
    vortex_core_span_fraction::Float64 = 0.5
    vortex_core_chord_fraction::Float64 = 0.0
    near_field_force_model::Symbol = :imperial
    initial_rpm::Float64 = 1207.96
    rpm_bracket::NTuple{2,Float64} = (1200, 1400)
    torque_tolerance_nm::Float64 = 0.25
    thrust_tolerance_n::Float64 = 1.0
    rpm_tolerance::Float64 = 0.25
    maximum_root_iterations::Int = 8
end

function validate_chang_windmilling_options(options::ChangWindmillingTrimOptions)
    options.flow_speed_mps > 0 || error("flow_speed_mps must be positive")
    options.air_density_kgpm3 > 0 || error("air_density_kgpm3 must be positive")
    options.propeller_radius_m > 0 || error("propeller_radius_m must be positive")
    options.propeller_chord_m > 0 || error("propeller_chord_m must be positive")
    options.blade_count > 0 || error("blade_count must be positive")
    options.radial_panels > 0 || error("radial_panels must be positive")
    options.chordwise_panels > 0 || error("chordwise_panels must be positive")
    isfinite(options.collective_pitch_offset_deg) ||
        error("collective_pitch_offset_deg must be finite")
    options.azimuth_step_deg > 0 || error("azimuth_step_deg must be positive")
    options.simulated_revolutions >= 2 || error("simulate at least two revolutions")
    1 <= options.averaged_revolutions < options.simulated_revolutions ||
        error("averaged_revolutions must be positive and smaller than simulated_revolutions")
    options.retained_wake_revolutions > 0 || error("retained_wake_revolutions must be positive")
    options.wake_relaxation >= 0 || error("wake_relaxation must be nonnegative")
    options.vortex_core_span_fraction >= 0 ||
        error("vortex_core_span_fraction must be nonnegative")
    options.vortex_core_chord_fraction >= 0 ||
        error("vortex_core_chord_fraction must be nonnegative")
    options.near_field_force_model in (:imperial, :legacy_imperial_segments) ||
        error(
            "near_field_force_model must be :imperial or " *
            ":legacy_imperial_segments",
        )
    options.initial_rpm > 0 || error("initial_rpm must be positive")
    0 < options.rpm_bracket[1] < options.rpm_bracket[2] ||
        error("rpm_bracket must contain two increasing positive values")
    options.torque_tolerance_nm > 0 || error("torque_tolerance_nm must be positive")
    options.thrust_tolerance_n > 0 || error("thrust_tolerance_n must be positive")
    options.rpm_tolerance > 0 || error("rpm_tolerance must be positive")
    options.maximum_root_iterations > 0 || error("maximum_root_iterations must be positive")

    steps_per_revolution = round(Int, 360 / options.azimuth_step_deg)
    isapprox(
        steps_per_revolution * options.azimuth_step_deg,
        360.0;
        atol = 100 * eps(Float64),
        rtol = 0.0,
    ) || error("azimuth_step_deg must divide 360 degrees exactly")
    return steps_per_revolution
end

"""Return the Chang blade angle at 75% radius, including any collective offset."""
chang_blade_angle_75_deg(options::ChangWindmillingTrimOptions) =
    get_twist_deg_interp(0.75) + options.collective_pitch_offset_deg

"""Return the near-field implementation selected by the trim options."""
function chang_trim_near_field_function(model::Symbol)
    model === :imperial && return near_field_forces!
    model === :legacy_imperial_segments && return legacy_imperial_segment_forces!
    error("Unsupported near-field force model: $model")
end

function rotate_chang_propeller_grids!(current_grids, reference_grids, omega, time)
    spin = RotationMatrix(-omega * time, 1)
    for blade_index in eachindex(reference_grids)
        reference_grid = reference_grids[blade_index]
        current_grid = current_grids[blade_index]
        for vertex_index in CartesianIndices((size(reference_grid, 2), size(reference_grid, 3)))
            i, j = Tuple(vertex_index)
            current_grid[:, i, j] = spin * SVector{3}(reference_grid[:, i, j])
        end
    end
    return current_grids
end

"""Return total rotor thrust and aerodynamic shaft torque from nodal loads."""
function chang_propeller_shaft_loads(system; hub_position = SVector(0.0, 0.0, 0.0))
    nodal_forces = imperial_nodal_forces(system)
    nodal_positions = imperial_nodal_positions(system)
    total_force = zero(hub_position)
    total_moment = zero(hub_position)

    for surface_index in eachindex(nodal_forces)
        for vertex_index in eachindex(nodal_forces[surface_index])
            force = nodal_forces[surface_index][vertex_index]
            position = nodal_positions[surface_index][vertex_index]
            total_force += force
            total_moment += cross(position - hub_position, force)
        end
    end

    # The UVLM body frame points +x downstream. Positive thrust is defined upstream.
    return (thrust_n = -total_force[1], torque_nm = total_moment[1])
end

"""
    simulate_chang_propeller_rpm(rpm, options; verbose=true)

March the isolated rigid Chang propeller with a free wake and average the
selected near-field shaft loads over complete final revolutions.
"""
function simulate_chang_propeller_rpm(
    rpm::Real,
    options::ChangWindmillingTrimOptions;
    verbose::Bool = true,
)
    steps_per_revolution = validate_chang_windmilling_options(options)
    rpm > 0 || error("rpm must be positive")

    omega = Float64(rpm) * 2pi / 60
    dt = deg2rad(options.azimuth_step_deg) / omega
    total_steps = options.simulated_revolutions * steps_per_revolution
    maximum_wake_rows = options.retained_wake_revolutions * steps_per_revolution
    near_field_force_function = chang_trim_near_field_function(
        options.near_field_force_model,
    )

    _, _, twists_at_nodes = get_nodal_properties_chang(options.radial_panels)
    blade_twists = deg2rad.(
        twists_at_nodes .- 90.0 .+ options.collective_pitch_offset_deg,
    )
    reference_grids = generate_propeller_blades_grid(
        options.propeller_radius_m,
        options.propeller_chord_m,
        options.radial_panels,
        options.chordwise_panels,
        blade_twists,
        options.blade_count,
    )
    current_grids = deepcopy(reference_grids)
    fcore = (chord, segment_width) -> max(
        options.vortex_core_span_fraction * segment_width,
        options.vortex_core_chord_fraction * chord,
    )
    initial_surfaces = [
        grid_to_surface_panels(grid; fcore = fcore)[3]
        for grid in current_grids
    ]
    system = System(initial_surfaces; nw = fill(maximum_wake_rows, options.blade_count))
    system.reference[] = Reference(
        pi * options.propeller_radius_m^2,
        options.propeller_chord_m,
        2 * options.propeller_radius_m,
        SVector(0.0, 0.0, 0.0),
        options.flow_speed_mps,
        options.air_density_kgpm3,
    )
    system.symmetric .= false
    system.surface_id .= 1:options.blade_count
    system.wake_finite_core .= true
    system.trailing_vortices .= false
    system.surfaces[:] = initial_surfaces
    system.previous_surfaces[:] = deepcopy(initial_surfaces)
    system.freestream[] = Freestream(
        options.flow_speed_mps,
        deg2rad(options.angle_of_attack_deg),
        deg2rad(options.sideslip_deg),
        SVector(0.0, 0.0, 0.0),
    )

    repeated_points = repeated_trailing_edge_points(system.surfaces)
    active_wake_rows = zeros(Int, options.blade_count)
    interaction_id = collect(1:options.blade_count)
    thrust_history = Vector{Float64}(undef, total_steps)
    torque_history = Vector{Float64}(undef, total_steps)

    verbose && @printf(
        "  RPM=%8.3f: %d steps, %d retained wake rows, model=%s\n",
        rpm,
        total_steps,
        maximum_wake_rows,
        string(options.near_field_force_model),
    )

    for step in 1:total_steps
        copy_surfaces_to_previous!(system, options.blade_count)
        time = step * dt
        rotate_chang_propeller_grids!(current_grids, reference_grids, omega, time)
        for blade_index in 1:options.blade_count
            system.surfaces[blade_index] = grid_to_surface_panels(
                current_grids[blade_index];
                fcore = fcore,
            )[3]
        end

        propagate_system!(
            system,
            system.freestream[],
            dt;
            additional_velocity = nothing,
            repeated_points = repeated_points,
            nwake = active_wake_rows,
            eta = options.wake_relaxation,
            calculate_influence_matrix = true,
            near_field_analysis = true,
            near_field_force_function = near_field_force_function,
            derivatives = false,
            interaction_id = interaction_id,
            interaction = true,
            advance_wake = true,
        )
        active_wake_rows .= min.(active_wake_rows .+ 1, maximum_wake_rows)

        loads = chang_propeller_shaft_loads(system)
        thrust_history[step] = loads.thrust_n
        torque_history[step] = loads.torque_nm
        if verbose && step % steps_per_revolution == 0
            completed_revolution = step ÷ steps_per_revolution
            revolution_range = (step - steps_per_revolution + 1):step
            @printf(
                "      revolution %d/%d: mean Q=%+.4f N m\n",
                completed_revolution,
                options.simulated_revolutions,
                mean(torque_history[revolution_range]),
            )
        end
    end

    revolution_mean_torque = [
        mean(torque_history[((rev - 1) * steps_per_revolution + 1):(rev * steps_per_revolution)])
        for rev in 1:options.simulated_revolutions
    ]
    revolution_mean_thrust = [
        mean(thrust_history[((rev - 1) * steps_per_revolution + 1):(rev * steps_per_revolution)])
        for rev in 1:options.simulated_revolutions
    ]
    average_steps = options.averaged_revolutions * steps_per_revolution
    average_range = (total_steps - average_steps + 1):total_steps
    mean_torque = mean(torque_history[average_range])
    mean_thrust = mean(thrust_history[average_range])
    torque_standard_deviation = std(torque_history[average_range])
    thrust_standard_deviation = std(thrust_history[average_range])
    n_revolutions_per_second = omega / (2pi)
    diameter = 2 * options.propeller_radius_m
    advance_ratio = options.flow_speed_mps / (n_revolutions_per_second * diameter)
    torque_coefficient = mean_torque / (
        options.air_density_kgpm3 * n_revolutions_per_second^2 * diameter^5
    )
    thrust_coefficient = mean_thrust / (
        options.air_density_kgpm3 * n_revolutions_per_second^2 * diameter^4
    )
    periodic_torque_change = revolution_mean_torque[end] - revolution_mean_torque[end - 1]
    periodic_thrust_change = revolution_mean_thrust[end] - revolution_mean_thrust[end - 1]

    return (
        near_field_force_model = options.near_field_force_model,
        rpm = Float64(rpm),
        omega_radps = omega,
        advance_ratio = advance_ratio,
        inflow_ratio = advance_ratio / pi,
        mean_thrust_n = mean_thrust,
        mean_torque_nm = mean_torque,
        thrust_coefficient = thrust_coefficient,
        torque_coefficient = torque_coefficient,
        torque_standard_deviation_nm = torque_standard_deviation,
        thrust_standard_deviation_n = thrust_standard_deviation,
        periodic_torque_change_nm = periodic_torque_change,
        periodic_thrust_change_n = periodic_thrust_change,
        revolution_mean_torque_nm = revolution_mean_torque,
        revolution_mean_thrust_n = revolution_mean_thrust,
        thrust_history_n = thrust_history,
        torque_history_nm = torque_history,
        steps_per_revolution = steps_per_revolution,
        time_step_s = dt,
    )
end

# Backward-compatible name used by earlier windmilling scripts.
simulate_chang_windmilling_rpm(args...; kwargs...) =
    simulate_chang_propeller_rpm(args...; kwargs...)

function _opposite_sign_or_zero(a, b)
    return iszero(a) || iszero(b) || signbit(a) != signbit(b)
end

"""Solve mean aerodynamic shaft torque = 0 with a safeguarded secant iteration."""
function trim_chang_windmilling_rpm(
    options::ChangWindmillingTrimOptions;
    verbose::Bool = true,
)
    validate_chang_windmilling_options(options)
    evaluations = Dict{Float64,NamedTuple}()

    function evaluate(rpm)
        key = round(Float64(rpm); digits = 8)
        result = get!(evaluations, key) do
            simulate_chang_propeller_rpm(key, options; verbose = verbose)
        end
        verbose && @printf(
            "    mean Q=%+10.4f N m, CQ=%+.6e, J=%.6f, rev drift=%+.4f N m\n",
            result.mean_torque_nm,
            result.torque_coefficient,
            result.advance_ratio,
            result.periodic_torque_change_nm,
        )
        return result
    end

    initial_result = evaluate(options.initial_rpm)
    if abs(initial_result.mean_torque_nm) <= options.torque_tolerance_nm
        return (trim = initial_result, initial = initial_result, evaluations = collect(values(evaluations)))
    end

    lower_rpm, upper_rpm = options.rpm_bracket
    lower_result = evaluate(lower_rpm)
    upper_result = evaluate(upper_rpm)
    expansion_count = 0
    while !_opposite_sign_or_zero(lower_result.mean_torque_nm, upper_result.mean_torque_nm)
        expansion_count += 1
        expansion_count <= 5 || error(
            "Unable to bracket zero torque after five RPM-bracket expansions",
        )
        lower_rpm = max(50.0, 0.75 * lower_rpm)
        upper_rpm = 1.25 * upper_rpm
        lower_result = evaluate(lower_rpm)
        upper_result = evaluate(upper_rpm)
    end

    if abs(lower_result.mean_torque_nm) <= options.torque_tolerance_nm
        ordered_evaluations = sort!(collect(values(evaluations)); by = result -> result.rpm)
        return (trim = lower_result, initial = initial_result, evaluations = ordered_evaluations)
    elseif abs(upper_result.mean_torque_nm) <= options.torque_tolerance_nm
        ordered_evaluations = sort!(collect(values(evaluations)); by = result -> result.rpm)
        return (trim = upper_result, initial = initial_result, evaluations = ordered_evaluations)
    end

    best_result = abs(lower_result.mean_torque_nm) <= abs(upper_result.mean_torque_nm) ?
        lower_result : upper_result
    for _ in 1:options.maximum_root_iterations
        torque_span = upper_result.mean_torque_nm - lower_result.mean_torque_nm
        trial_rpm = upper_rpm - upper_result.mean_torque_nm *
            (upper_rpm - lower_rpm) / torque_span
        bracket_width = upper_rpm - lower_rpm
        edge_margin = 0.1 * bracket_width
        if !isfinite(trial_rpm) ||
                trial_rpm <= lower_rpm + edge_margin || trial_rpm >= upper_rpm - edge_margin
            trial_rpm = (lower_rpm + upper_rpm) / 2
        end

        trial_result = evaluate(trial_rpm)
        if abs(trial_result.mean_torque_nm) < abs(best_result.mean_torque_nm)
            best_result = trial_result
        end
        if abs(trial_result.mean_torque_nm) <= options.torque_tolerance_nm ||
                bracket_width <= options.rpm_tolerance
            best_result = trial_result
            break
        end

        if _opposite_sign_or_zero(lower_result.mean_torque_nm, trial_result.mean_torque_nm)
            upper_rpm, upper_result = trial_rpm, trial_result
        else
            lower_rpm, lower_result = trial_rpm, trial_result
        end
    end

    abs(best_result.mean_torque_nm) <= options.torque_tolerance_nm || @warn(
        "Windmilling trim stopped above the requested torque tolerance",
        residual_torque_nm = best_result.mean_torque_nm,
        torque_tolerance_nm = options.torque_tolerance_nm,
    )
    ordered_evaluations = sort!(collect(values(evaluations)); by = result -> result.rpm)
    return (trim = best_result, initial = initial_result, evaluations = ordered_evaluations)
end

"""
    trim_chang_thrusting_rpm(options; target_thrust_n, verbose=true)

Solve `mean_thrust_n = target_thrust_n` by varying RPM. Unlike windmilling
trim, the resulting aerodynamic shaft torque is generally nonzero and must be
balanced by the motor.
"""
function trim_chang_thrusting_rpm(
    options::ChangWindmillingTrimOptions;
    target_thrust_n::Real,
    verbose::Bool = true,
)
    validate_chang_windmilling_options(options)
    isfinite(target_thrust_n) && target_thrust_n > 0 ||
        error("target_thrust_n must be finite and positive")
    target_thrust = Float64(target_thrust_n)
    evaluations = Dict{Float64,NamedTuple}()

    function evaluate(rpm)
        key = round(Float64(rpm); digits = 8)
        result = get!(evaluations, key) do
            simulate_chang_propeller_rpm(key, options; verbose = verbose)
        end
        residual = result.mean_thrust_n - target_thrust
        verbose && @printf(
            "    mean T=%+10.4f N, residual=%+10.4f N, CT=%+.6e, J=%.6f, rev drift=%+.4f N\n",
            result.mean_thrust_n,
            residual,
            result.thrust_coefficient,
            result.advance_ratio,
            result.periodic_thrust_change_n,
        )
        return result
    end

    thrust_residual(result) = result.mean_thrust_n - target_thrust
    initial_result = evaluate(options.initial_rpm)
    if abs(thrust_residual(initial_result)) <= options.thrust_tolerance_n
        return (
            trim = initial_result,
            initial = initial_result,
            target_thrust_n = target_thrust,
            residual_thrust_n = thrust_residual(initial_result),
            evaluations = collect(values(evaluations)),
        )
    end

    lower_rpm, upper_rpm = options.rpm_bracket
    lower_result = evaluate(lower_rpm)
    upper_result = evaluate(upper_rpm)
    lower_residual = thrust_residual(lower_result)
    upper_residual = thrust_residual(upper_result)
    expansion_count = 0
    while !_opposite_sign_or_zero(lower_residual, upper_residual)
        expansion_count += 1
        expansion_count <= 5 || error(
            "Unable to bracket the target thrust after five RPM-bracket expansions",
        )
        lower_rpm = max(50.0, 0.75 * lower_rpm)
        upper_rpm = 1.25 * upper_rpm
        lower_result = evaluate(lower_rpm)
        upper_result = evaluate(upper_rpm)
        lower_residual = thrust_residual(lower_result)
        upper_residual = thrust_residual(upper_result)
    end

    if abs(lower_residual) <= options.thrust_tolerance_n
        ordered_evaluations = sort!(collect(values(evaluations)); by = result -> result.rpm)
        return (
            trim = lower_result,
            initial = initial_result,
            target_thrust_n = target_thrust,
            residual_thrust_n = lower_residual,
            evaluations = ordered_evaluations,
        )
    elseif abs(upper_residual) <= options.thrust_tolerance_n
        ordered_evaluations = sort!(collect(values(evaluations)); by = result -> result.rpm)
        return (
            trim = upper_result,
            initial = initial_result,
            target_thrust_n = target_thrust,
            residual_thrust_n = upper_residual,
            evaluations = ordered_evaluations,
        )
    end

    best_result = abs(lower_residual) <= abs(upper_residual) ?
        lower_result : upper_result
    for _ in 1:options.maximum_root_iterations
        residual_span = upper_residual - lower_residual
        trial_rpm = upper_rpm - upper_residual *
            (upper_rpm - lower_rpm) / residual_span
        bracket_width = upper_rpm - lower_rpm
        edge_margin = 0.1 * bracket_width
        if !isfinite(trial_rpm) ||
                trial_rpm <= lower_rpm + edge_margin ||
                trial_rpm >= upper_rpm - edge_margin
            trial_rpm = (lower_rpm + upper_rpm) / 2
        end

        trial_result = evaluate(trial_rpm)
        trial_residual = thrust_residual(trial_result)
        if abs(trial_residual) < abs(thrust_residual(best_result))
            best_result = trial_result
        end
        if abs(trial_residual) <= options.thrust_tolerance_n ||
                bracket_width <= options.rpm_tolerance
            best_result = trial_result
            break
        end

        if _opposite_sign_or_zero(lower_residual, trial_residual)
            upper_rpm, upper_result, upper_residual =
                trial_rpm, trial_result, trial_residual
        else
            lower_rpm, lower_result, lower_residual =
                trial_rpm, trial_result, trial_residual
        end
    end

    final_residual = thrust_residual(best_result)
    abs(final_residual) <= options.thrust_tolerance_n || @warn(
        "Thrusting trim stopped above the requested thrust tolerance",
        residual_thrust_n = final_residual,
        thrust_tolerance_n = options.thrust_tolerance_n,
    )
    ordered_evaluations = sort!(collect(values(evaluations)); by = result -> result.rpm)
    return (
        trim = best_result,
        initial = initial_result,
        target_thrust_n = target_thrust,
        residual_thrust_n = final_residual,
        evaluations = ordered_evaluations,
    )
end
