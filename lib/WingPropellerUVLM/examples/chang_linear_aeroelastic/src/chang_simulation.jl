# Time integration and runtime checks for the Chang aeroelastic example.
#
# The entry-point file keeps the physical setup visible. This file contains the
# step-level bookkeeping: structural history, trim-load subtraction, coupling
# diagnostics, accepted-wake updates, optional animation snapshots, and aborts.

# ### Runtime option groups

"""Read plotting and animation controls shared by interactive and batch runs."""
function chang_visualization_options()
    options = (;
        plot_results = environment_flag("CHANG_PLOT_RESULTS", true),
        animate_wake = environment_flag("CHANG_ANIMATE_WAKE", false),
        animation_stride = environment_number(Int, "CHANG_ANIMATION_STRIDE", 5),
        animation_fps = environment_number(Int, "CHANG_ANIMATION_FPS", 15),
        plot_time_limit_s = environment_number(
            Float64,
            "CHANG_PLOT_END_TIME_S",
            5.0,
        ),
    )
    options.animation_stride > 0 || error("CHANG_ANIMATION_STRIDE must be positive")
    options.animation_fps > 0 || error("CHANG_ANIMATION_FPS must be positive")
    return options
end

"""Read the UVLM regularization, load-arm, and retained-wake controls."""
function chang_aerodynamic_options(wing_chord_panels::Int, pylon_length::Real)
    segment_core_factor = environment_number(
        Float64,
        "CHANG_FCORE_SEGMENT_FACTOR",
        0.01,
    )
    chord_core_factor = environment_number(
        Float64,
        "CHANG_FCORE_CHORD_FACTOR",
        0.01,
    )
    hub_load_arm_factor = environment_number(
        Float64,
        "CHANG_HUB_LOAD_ARM_FACTOR",
        0.5,
    )
    maximum_wake_rows_wing = environment_number(
        Int,
        "CHANG_WAKE_ROWS_WING",
        10 * wing_chord_panels,
    )
    maximum_wake_rows_propeller = environment_number(
        Int,
        "CHANG_WAKE_ROWS_PROPELLER",
        72,
    )

    isfinite(segment_core_factor) && segment_core_factor >= 0.0 || error(
        "CHANG_FCORE_SEGMENT_FACTOR must be finite and nonnegative",
    )
    isfinite(chord_core_factor) && chord_core_factor >= 0.0 || error(
        "CHANG_FCORE_CHORD_FACTOR must be finite and nonnegative",
    )
    hub_load_arm_factor >= 0.0 || error(
        "CHANG_HUB_LOAD_ARM_FACTOR must be nonnegative",
    )
    maximum_wake_rows_wing >= 0 || error(
        "CHANG_WAKE_ROWS_WING must be nonnegative",
    )
    maximum_wake_rows_propeller >= 0 || error(
        "CHANG_WAKE_ROWS_PROPELLER must be nonnegative",
    )

    finite_core = (chord_length, segment_length) -> max(
        segment_core_factor * segment_length,
        chord_core_factor * chord_length,
    )
    return (;
        segment_core_factor,
        chord_core_factor,
        finite_core,
        elastic_axis_fraction = 0.30,
        propeller_pivot_offset_A = SVector(0.0, 0.0, 0.0),
        physical_hub_center_A = SVector(-pylon_length, 0.0, 0.0),
        load_center_A = SVector(-hub_load_arm_factor * pylon_length, 0.0, 0.0),
        maximum_wake_rows_wing,
        maximum_wake_rows_propeller,
    )
end

"""Read the pre-impulse trim window and smooth pitch-pulse controls."""
function chang_excitation_options(
    angular_speed::Real,
    number_of_propellers::Int,
    propeller_indices,
)
    revolution_period = 2pi / abs(angular_speed)
    trim_revolutions = environment_number(
        Float64,
        "CHANG_TRIM_REVOLUTIONS",
        10.0,
    )
    trim_average_revolutions = environment_number(
        Float64,
        "CHANG_TRIM_AVERAGE_REVOLUTIONS",
        1.0,
    )
    start_time = environment_number(
        Float64,
        "CHANG_IMPULSE_START_S",
        trim_revolutions * revolution_period,
    )
    duration = environment_number(Float64, "CHANG_IMPULSE_DURATION_S", 0.15)

    trim_revolutions > 0.0 || error("CHANG_TRIM_REVOLUTIONS must be positive")
    trim_average_revolutions > 0.0 || error(
        "CHANG_TRIM_AVERAGE_REVOLUTIONS must be positive",
    )
    start_time >= 0.0 || error("CHANG_IMPULSE_START_S must be nonnegative")
    duration > 0.0 || error("CHANG_IMPULSE_DURATION_S must be positive")

    return (;
        propeller_indices,
        number_of_propellers,
        magnitude = environment_number(
            Float64,
            "CHANG_IMPULSE_MAGNITUDE",
            1000.0,
        ),
        start_time,
        duration,
        trim_average_revolutions,
        trim_average_start_time = max(
            0.0,
            start_time - trim_average_revolutions * revolution_period,
        ),
    )
end

"""Read generalized-alpha, partitioned-coupling, and safety controls."""
function chang_integration_options(
    structural;
    reference,
    freestream_speed::Real,
    reference_area::Real,
    reference_chord::Real,
    propeller_radius::Real,
    wing_node_count::Int,
    dofs_per_node::Int,
)
    generalized_alpha = generalized_alpha_parameters(
        environment_number(Float64, "CHANG_GA_RHO_INF", 1.0),
    )
    coupling = PartitionedCouplingOptions(
        maximum_iterations = environment_number(
            Int,
            "CHANG_COUPLING_MAX_ITER",
            10,
        ),
        state_tolerance = environment_number(
            Float64,
            "CHANG_COUPLING_TOL_U",
            1.0e-5,
        ),
        load_tolerance = environment_number(
            Float64,
            "CHANG_COUPLING_TOL_F",
            1.0e-2,
        ),
        equilibrium_tolerance = environment_number(
            Float64,
            "CHANG_COUPLING_TOL_EQ",
            1.0e-10,
        ),
        coupled_equilibrium_tolerance = environment_number(
            Float64,
            "CHANG_COUPLING_TOL_COUPLED_EQ",
            1.0e-4,
        ),
        relaxation = environment_number(
            Float64,
            "CHANG_COUPLING_RELAXATION",
            1.0,
        ),
    )
    scales = build_chang_coupling_scales(
        structural;
        reference,
        freestream_speed,
        reference_area,
        reference_chord,
        propeller_radius,
        wing_node_count,
        dofs_per_node,
    )
    state_norm_limit = environment_number(
        Float64,
        "CHANG_STATE_ABORT_NORM",
        1.0e3,
    )
    propeller_angle_limit_deg = environment_number(
        Float64,
        "CHANG_PROP_ANGLE_ABORT_DEG",
        Inf,
    )
    state_norm_limit > 0.0 || error("CHANG_STATE_ABORT_NORM must be positive")
    propeller_angle_limit_deg > 0.0 || error(
        "CHANG_PROP_ANGLE_ABORT_DEG must be positive",
    )
    return (;
        generalized_alpha,
        coupling,
        state_scale = scales.state,
        load_scale = scales.load,
        state_norm_limit,
        propeller_angle_limit_deg,
    )
end

"""Print the numerical choices that govern the coupled response."""
function report_chang_solver_options(excitation, integration)
    parameters = integration.generalized_alpha
    coupling = integration.coupling
    println(
        "Partitioned generalized-alpha: rho_inf=$(parameters.rho_inf), " *
        "alpha_m=$(parameters.alpha_m), alpha_f=$(parameters.alpha_f), " *
        "gamma=$(parameters.gamma), beta=$(parameters.beta)",
    )
    println(
        "Coupling correction: max_iter=$(coupling.maximum_iterations), " *
        "tol_u=$(coupling.state_tolerance), tol_f=$(coupling.load_tolerance), " *
        "tol_linear=$(coupling.equilibrium_tolerance), " *
        "tol_coupled=$(coupling.coupled_equilibrium_tolerance), " *
        "relaxation=$(coupling.relaxation)",
    )
    println(
        "Trim baseline: average $(excitation.trim_average_revolutions) revolution(s), " *
        "from t=$(round(excitation.trim_average_start_time, digits = 4)) s " *
        "to t=$(round(excitation.start_time, digits = 4)) s",
    )
    return nothing
end

# ### Structural preparation

"""Print the structural diagnostics and reject an invalid reduced model."""
function report_and_validate_structural_model(
    structural;
    symmetry_tolerance::Real = 1.0e-12,
)
    diagnostics = chang_structural_diagnostics(structural)
    println(
        "Wing-only modal frequencies (Hz): " *
        "$(round.(diagnostics.wing_modal_frequencies_hz, digits = 4))",
    )
    println(
        "Structural discretization: mass=$(structural.mass_model), " *
        "stiffness=$(structural.stiffness_model)",
    )
    println(
        "Structural checks: min eig(M)=$(diagnostics.minimum_mass_eigenvalue), " *
        "min eig(K)=$(diagnostics.minimum_stiffness_eigenvalue), " *
        "symmetry(M/K)=($(diagnostics.mass_symmetry_error), " *
        "$(diagnostics.stiffness_symmetry_error))",
    )

    diagnostics.minimum_mass_eigenvalue > 0.0 || error(
        "The free structural mass matrix is not positive definite",
    )
    diagnostics.minimum_stiffness_eigenvalue > 0.0 || error(
        "The free structural stiffness matrix is not positive definite",
    )
    diagnostics.mass_symmetry_error <= symmetry_tolerance || error(
        "The free structural mass matrix is not symmetric",
    )
    diagnostics.stiffness_symmetry_error <= symmetry_tolerance || error(
        "The free structural stiffness matrix is not symmetric",
    )
    return diagnostics
end

"""Allocate displacement, velocity, and acceleration at every output time."""
function initialize_structural_history(structural, number_of_states::Int)
    displacement = Vector{Vector{Float64}}(undef, number_of_states)
    velocity = similar(displacement)
    acceleration = similar(displacement)

    displacement[1] = zeros(structural.ndof_free)
    velocity[1] = zeros(structural.ndof_free)
    acceleration[1] = isempty(structural.M) ? zeros(structural.ndof_free) :
        structural.M \ (
            -structural.C * velocity[1] - structural.K * displacement[1]
        )
    return (; displacement, velocity, acceleration)
end

"""Build dimensional reference scales for the partitioned residuals."""
function build_chang_coupling_scales(
    structural;
    reference,
    freestream_speed::Real,
    reference_area::Real,
    reference_chord::Real,
    propeller_radius::Real,
    wing_node_count::Int,
    dofs_per_node::Int,
)
    state = ones(structural.ndof_free)
    load = ones(structural.ndof_free)
    force = max(
        0.5 * reference.rho * freestream_speed^2 * reference_area,
        1.0,
    )
    moment = max(force * reference_chord, 1.0)

    # The root is clamped, so the reduced wing state contains nnodes - 1 nodes.
    for node_index in 1:(wing_node_count - 1)
        node_start = dofs_per_node * (node_index - 1)
        state[node_start .+ (1:3)] .= max(reference_chord, eps(Float64))
        state[node_start .+ (4:6)] .= 1.0
        load[node_start .+ (1:3)] .= force
        load[node_start .+ (4:6)] .= moment
    end

    propeller_moment = max(
        0.5 * reference.rho * freestream_speed^2 * pi * propeller_radius^3,
        1.0,
    )
    load[(structural.ndof_wing_free + 1):end] .= propeller_moment
    return (; state, load)
end

# ### Accepted-step helpers

"""Print the generalized trim load captured before the pitch impulse."""
function report_trim_load(
    trim_load,
    sample_count::Int,
    capture_time::Real,
    wing_dof_count::Int,
    number_of_propellers::Int,
)
    println(
        "\n=== Mean trim generalized load captured at " *
        "t=$(round(capture_time, digits = 4)) s from $sample_count samples ===",
    )
    for propeller_index in 1:number_of_propellers
        pitch_index = wing_dof_count + 2 * (propeller_index - 1) + 1
        yaw_index = pitch_index + 1
        println(
            "  P$propeller_index: pitch F0 = " *
            "$(round(trim_load[pitch_index], digits = 3)) N*m, " *
            "yaw F0 = $(round(trim_load[yaw_index], digits = 3)) N*m",
        )
    end
    println(
        "  |F0_wing| = $(round(norm(trim_load[1:wing_dof_count]), digits = 2)), " *
        "|F0_prop| = $(round(norm(trim_load[(wing_dof_count + 1):end]), digits = 2))",
    )
    return nothing
end

"""Advance and activate the wake exactly once after accepting a coupled step."""
function commit_chang_wake!(system, freestream, time_step, wake)
    advance_wake!(
        system,
        freestream,
        time_step;
        additional_velocity = nothing,
        repeated_points = wake.repeated_points,
        nwake = wake.active_rows,
        interaction_id = wake.interaction_ids,
        interaction = wake.interaction_on,
    )
    for surface_index in eachindex(wake.active_rows)
        wake.active_rows[surface_index] < wake.maximum_rows[surface_index] || continue
        wake.active_rows[surface_index] += 1
    end
    return nothing
end

"""Describe whether an accepted structural state exceeds a safety limit."""
function chang_state_limit_status(
    state,
    wing_dof_count::Int;
    state_norm_limit::Real,
    propeller_angle_limit_deg::Real,
)
    maximum_state = maximum(abs, state)
    propeller_angles_deg = rad2deg.(state[(wing_dof_count + 1):end])
    maximum_propeller_angle_deg = maximum(abs, propeller_angles_deg)
    nonfinite = any(value -> !isfinite(value), state)
    state_limit_exceeded = maximum_state > state_norm_limit
    angle_limit_exceeded = maximum_propeller_angle_deg > propeller_angle_limit_deg
    reason = nonfinite ? "non-finite structural state" :
        state_limit_exceeded ? "state limit $state_norm_limit" :
        angle_limit_exceeded ? "propeller-angle limit $propeller_angle_limit_deg deg" : ""
    return (;
        exceeded = nonfinite || state_limit_exceeded || angle_limit_exceeded,
        reason,
        maximum_state,
        maximum_propeller_angle_deg,
    )
end

# ### Coupled solution

"""
    solve_chang_aeroelastic!(system, structural; kwargs...)

March the coupled generalized-alpha/UVLM problem. Every aerodynamic trial is
evaluated from the same beginning-of-step snapshot. A failed coupled step is
restored and rejected; a converged step advances the free wake exactly once.
"""
function solve_chang_aeroelastic!(
    system,
    structural;
    time,
    time_steps,
    freestream_history,
    aerodynamic_load,
    wake,
    excitation,
    integration,
    animation,
)
    history = initialize_structural_history(structural, length(time))
    displacement = history.displacement
    velocity = history.velocity
    acceleration = history.acceleration

    pitch_dof_indices = [
        structural.ndof_wing_free + 2 * (index - 1) + 1
        for index in excitation.propeller_indices
    ]
    trim_load = zeros(structural.ndof_free)
    trim_load_accumulator = zeros(structural.ndof_free)
    trim_sample_count = 0
    have_trim_load = false
    impulse_was_reported = false

    coupling_iterations = fill(0, length(time_steps))
    coupling_state_residual = fill(NaN, length(time_steps))
    coupling_load_residual = fill(NaN, length(time_steps))
    coupling_equilibrium_residual = fill(NaN, length(time_steps))
    coupling_converged = fill(false, length(time_steps))
    accepted_perturbation_load = zeros(structural.ndof_free)
    last_accepted_step = length(time_steps)

    animation_surface_history = Vector{typeof(system.surfaces)}()
    animation_wake_history = Vector{typeof(system.wakes)}()
    animation_active_wake_rows_history = Vector{Vector{Int}}()
    animation_time_history = Float64[]
    if animation.enabled
        animation.record_frame(
            animation_surface_history,
            animation_wake_history,
            animation_active_wake_rows_history,
            animation_time_history,
            time[1],
        )
    end

    println("Starting coupled aeroelastic simulation...")
    for step in eachindex(time_steps)
        # All fixed-point trials start from this exact accepted aerodynamic state.
        copy_surfaces_to_previous!(system, wake.surface_count)
        aerodynamic_snapshot = snapshot_uvlm(system)
        time_step = time_steps[step]

        external_load_n = smooth_hann_pulse_load(
            time[step],
            structural.ndof_free,
            pitch_dof_indices;
            magnitude = excitation.magnitude,
            start_time = excitation.start_time,
            duration = excitation.duration,
        )
        external_load_np1 = smooth_hann_pulse_load(
            time[step + 1],
            structural.ndof_free,
            pitch_dof_indices;
            magnitude = excitation.magnitude,
            start_time = excitation.start_time,
            duration = excitation.duration,
        )
        if !impulse_was_reported &&
            (maximum(abs, external_load_n) > 0.0 || maximum(abs, external_load_np1) > 0.0)
            println("<<<<<<<<<<< Applying Smooth Pitch Impulse >>>>>>>>>>>")
            impulse_was_reported = true
        end

        # Converge the structural state and UVLM load while holding the wake fixed.
        full_aerodynamic_load = zeros(structural.ndof_free)
        step_solution = partitioned_generalized_alpha_step(
            structural.M,
            structural.C,
            structural.K,
            displacement[step],
            velocity[step],
            acceleration[step],
            accepted_perturbation_load,
            external_load_np1,
            external_load_n,
            time_step,
            integration.generalized_alpha,
            state_guess -> begin
                full_aerodynamic_load .= aerodynamic_load(
                    aerodynamic_snapshot,
                    state_guess,
                    step,
                )
                return have_trim_load ? full_aerodynamic_load .- trim_load :
                    zeros(structural.ndof_free)
            end;
            options = integration.coupling,
            require_load_convergence = have_trim_load,
            state_scale = integration.state_scale,
            load_scale = integration.load_scale,
        )

        if !step_solution.converged
            restore_uvlm!(system, aerodynamic_snapshot)
            error(
                "Partitioned coupling failed at step $step " *
                "(t=$(time[step + 1]) s) after $(step_solution.iterations) iterations: " *
                "state_res=$(step_solution.state_residual), " *
                "load_res=$(step_solution.load_residual), " *
                "coupled_equilibrium_res=$(step_solution.coupled_equilibrium_residual), " *
                "linear_equilibrium_res=$(step_solution.equilibrium_residual). " *
                "The time step was not committed.",
            )
        end

        displacement[step + 1] = copy(step_solution.displacement)
        velocity[step + 1] = copy(step_solution.velocity)
        acceleration[step + 1] = copy(step_solution.acceleration)
        accepted_full_load = copy(full_aerodynamic_load)
        accepted_perturbation_load .= step_solution.trial_load

        # Average the periodic full load, then subtract it from later trials so
        # the structural solution contains aerodynamic perturbations only.
        accepted_time = time[step + 1]
        trim_sample = !have_trim_load &&
            accepted_time >= excitation.trim_average_start_time - eps(accepted_time) &&
            accepted_time <= excitation.start_time + eps(accepted_time)
        if trim_sample
            trim_load_accumulator .+= accepted_full_load
            trim_sample_count += 1
        end
        if !have_trim_load && accepted_time >= excitation.start_time - eps(accepted_time)
            trim_sample_count > 0 || error(
                "No aerodynamic samples were available for trim averaging",
            )
            trim_load .= trim_load_accumulator ./ trim_sample_count
            have_trim_load = true
            accepted_perturbation_load .= 0.0
            report_trim_load(
                trim_load,
                trim_sample_count,
                accepted_time,
                structural.ndof_wing_free,
                excitation.number_of_propellers,
            )
        end

        coupling_iterations[step] = step_solution.iterations
        coupling_state_residual[step] = step_solution.state_residual
        coupling_load_residual[step] = step_solution.load_residual
        coupling_equilibrium_residual[step] = step_solution.coupled_equilibrium_residual
        coupling_converged[step] = step_solution.converged

        # The converged circulation/load trial is already in `system`; only now
        # is its wake convected and the newly shed row activated.
        commit_chang_wake!(system, freestream_history[step], time_step, wake)

        if animation.enabled &&
            (step == 1 || step % animation.stride == 0 || step == length(time_steps))
            animation.record_frame(
                animation_surface_history,
                animation_wake_history,
                animation_active_wake_rows_history,
                animation_time_history,
                accepted_time,
            )
        end

        println(
            "Step $step/$(length(time_steps)) (t=$(round(accepted_time, digits = 4)) s): " *
            "$(step_solution.iterations) iterations, " *
            "state_res=$(round(step_solution.state_residual, sigdigits = 4)), " *
            "load_res=$(round(step_solution.load_residual, sigdigits = 4)), " *
            "coupled_eq_res=$(round(step_solution.coupled_equilibrium_residual, sigdigits = 4)), " *
            "linear_eq_res=$(round(step_solution.equilibrium_residual, sigdigits = 4))",
        )
        flush(stdout)

        limit = chang_state_limit_status(
            displacement[step + 1],
            structural.ndof_wing_free;
            state_norm_limit = integration.state_norm_limit,
            propeller_angle_limit_deg = integration.propeller_angle_limit_deg,
        )
        if limit.exceeded
            sample_range = max(1, step - 400):step
            pitch_envelope = [
                abs(displacement[index][pitch_dof_indices[1]])
                for index in sample_range
            ]
            sample_time = [time[index] for index in sample_range]
            positive_indices = findall(>(1.0e-12), pitch_envelope)
            message =
                "\n>>> ABORT at t=$(round(accepted_time, digits = 4)) s, " *
                "reason=$(limit.reason), " *
                "|U|=$(round(limit.maximum_state, sigdigits = 4)), " *
                "max propeller angle=" *
                "$(round(limit.maximum_propeller_angle_deg, sigdigits = 4)) deg"
            if length(positive_indices) > 5
                first_index, last_index = first(positive_indices), last(positive_indices)
                growth_rate = (
                    log(pitch_envelope[last_index]) - log(pitch_envelope[first_index])
                ) / (sample_time[last_index] - sample_time[first_index])
                message *= "; growth sigma(P1 pitch) approx " *
                    "$(round(growth_rate, digits = 4)) /s"
            end
            println(message)
            last_accepted_step = step
            break
        end

        if step in wake.saved_steps
            wake.surface_history[step] = [copy(surface) for surface in system.surfaces]
        end
    end
    println("Simulation finished.")

    if animation.enabled && last(animation_time_history) != time[last_accepted_step + 1]
        animation.record_frame(
            animation_surface_history,
            animation_wake_history,
            animation_active_wake_rows_history,
            animation_time_history,
            time[last_accepted_step + 1],
        )
    end

    return (;
        displacement_history = displacement,
        velocity_history = velocity,
        acceleration_history = acceleration,
        last_step = last_accepted_step,
        trim_load,
        coupling_iterations,
        coupling_state_residual,
        coupling_load_residual,
        coupling_equilibrium_residual,
        coupling_converged,
        animation_surface_history,
        animation_wake_history,
        animation_active_wake_rows_history,
        animation_time_history,
    )
end
