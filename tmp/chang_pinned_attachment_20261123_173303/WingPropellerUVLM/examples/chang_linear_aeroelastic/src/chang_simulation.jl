# Time integration and checks for the Chang aeroelastic example.
# Options below derive solver data from the explicit, resolved configuration.
# The solver owns structural history, trim subtraction, and accepted wake updates.

"""Create the solver's frame callback; disabled animation needs no plotting package."""
function chang_animation_options(system,active_rows,visualization)
    record_frame = if visualization.animate_wake
        (surfaces,wakes,rows,times,time) -> record_chang_animation_frame!(
            surfaces,wakes,rows,times,system,active_rows,time
        )
    else
        (arguments...) -> nothing
    end
    return (;
        enabled=visualization.animate_wake,
        stride=visualization.animation_stride,
        record_frame
    )
end


"""Build the finite-core law and load reference points for one configured case."""
function chang_aerodynamic_options(config)
    (; core_radius_m,segment_core_factor,chord_core_factor,hub_load_arm_factor,
        elastic_axis_fraction) = config.aerodynamic
    pylonLength = config.structural.pylon_length_m
    finite_core = if isnothing(core_radius_m)
        (chordLength,segmentLength) -> max(
            segment_core_factor * segmentLength,
            chord_core_factor * chordLength
        )
    else
        # Choose the law once; panel construction reuses this captured radius.
        radius = Float64(core_radius_m)
        (chordLength,segmentLength) -> radius
    end
    return (;
        core_radius_m,segment_core_factor,chord_core_factor,finite_core,elastic_axis_fraction,
        # The pylon pivot is attached directly to the wing elastic-axis node.
        # The physical rotor hub is one full pylon length ahead of that pivot;
        # the modal load point follows the Chang assumed-mode arm.
        propeller_pivot_offset_A=SVector(0.0,0.0,0.0),
        physical_hub_center_A=SVector(-pylonLength,0.0,0.0),
        load_center_A=SVector(-hub_load_arm_factor * pylonLength,0.0,0.0),
        maximum_wake_rows_wing=config.wake.maximum_rows_wing,
        maximum_wake_rows_propeller=config.wake.maximum_rows_propeller
    )
end

"""Derive the trim averaging window and the smooth pitch pulse from the case."""
function chang_excitation_options(config,parameters)
    excitation = config.excitation
    if parameters.Npropellers == 0
        # Apply a nodal moment to the wing itself, using the structural DOF basis.
        node = clamp(round(Int,excitation.wing_impulse_eta * parameters.Ne) + 1,
            2,parameters.nnodes)
        return (;
            dof_indices=[(node - 2) * parameters.ndof + excitation.wing_impulse_dof],
            number_of_propellers=0,
            magnitude=excitation.impulse_magnitude_nm,
            start_time=excitation.impulse_start_s,
            duration=excitation.impulse_duration_s,
            trim_average_duration_s=excitation.trim_average_duration_s,
            trim_average_start_time=max(0.0,
                excitation.impulse_start_s - excitation.trim_average_duration_s)
        )
    end
    revolutionPeriod = 2pi / abs(parameters.Ω)
    return (;
        propeller_indices=config.simulation.impulse_propeller_indices,
        number_of_propellers=parameters.Npropellers,
        magnitude=excitation.impulse_magnitude_nm,
        start_time=excitation.impulse_start_s,
        duration=excitation.impulse_duration_s,
        trim_average_revolutions=excitation.trim_average_revolutions,
        trim_average_start_time=max(
            0.0,
            excitation.impulse_start_s - excitation.trim_average_revolutions * revolutionPeriod
        )
    )
end

"""Build independent structural-integration/coupling options and residual scales."""
function chang_integration_options(config,structural,parameters)
    scales = build_chang_coupling_scales(
        structural;
        reference=parameters.ref,
        freestream_speed=parameters.Vinf,
        reference_area=parameters.Sref,
        reference_chord=parameters.cref,
        propeller_radius=parameters.R_prop,
        wing_node_count=parameters.nnodes,
        dofs_per_node=parameters.ndof
    )
    time_integrator = validate_structural_time_integrator(
        config.integration.time_integrator
    )
    coupling_scheme = validate_aeroelastic_coupling_scheme(config.coupling.scheme)
    newmark_beta = newmark_beta_parameters(config.integration.newmark_alpha)
    generalized_alpha = generalized_alpha_parameters(config.integration.rho_inf)
    structural_parameters = time_integrator == :newmark_beta ?
                            newmark_beta : generalized_alpha
    coupling = PartitionedCouplingOptions(
        maximum_iterations=config.coupling.maximum_iterations,
        state_tolerance=config.coupling.state_tolerance,
        load_tolerance=config.coupling.load_tolerance,
        equilibrium_tolerance=config.coupling.equilibrium_tolerance,
        coupled_equilibrium_tolerance=config.coupling.coupled_equilibrium_tolerance,
        relaxation=config.coupling.relaxation
    )
    return (;
        time_integrator,
        coupling_scheme,
        structural_parameters,
        newmark_beta,
        generalized_alpha,
        coupling,
        coupling_verbose=config.coupling.verbose,
        state_scale=scales.state,
        load_scale=scales.load,
        state_norm_limit=config.integration.state_norm_limit,
        propeller_angle_limit_deg=config.integration.propeller_angle_limit_deg
    )
end

"""Print the numerical choices that govern the coupled response."""
function report_chang_solver_options(excitation,integration)
    coupling = integration.coupling
    if integration.time_integrator == :newmark_beta
        parameters = integration.newmark_beta
        println(
            "Structural integrator: Newmark-beta, alpha=$(parameters.alpha), " *
            "gamma=$(parameters.gamma), beta=$(parameters.beta)"
        )
    else
        parameters = integration.generalized_alpha
        println(
            "Structural integrator: generalized-alpha, rho_inf=$(parameters.rho_inf), " *
            "alpha_m=$(parameters.alpha_m), alpha_f=$(parameters.alpha_f), " *
            "gamma=$(parameters.gamma), beta=$(parameters.beta)"
        )
    end
    println("Aeroelastic coupling: $(integration.coupling_scheme)")
    println(
        "Coupling correction: max_iter=$(coupling.maximum_iterations), " *
        "tol_u=$(coupling.state_tolerance), tol_f=$(coupling.load_tolerance), " *
        "tol_linear=$(coupling.equilibrium_tolerance), " *
        "tol_coupled=$(coupling.coupled_equilibrium_tolerance), " *
        "relaxation=$(coupling.relaxation), verbose=$(integration.coupling_verbose)"
    )
    averageWindow =
        hasproperty(excitation,:trim_average_duration_s) ?
        "$(excitation.trim_average_duration_s) s" :
        "$(excitation.trim_average_revolutions) revolution(s)"
    println(
        "Trim baseline: average $averageWindow, " *
        "from t=$(round(excitation.trim_average_start_time,digits = 4)) s " *
        "to t=$(round(excitation.start_time,digits = 4)) s"
    )
    return nothing
end

# ### Structural preparation

"""Print the structural diagnostics and reject an invalid reduced model."""
function report_and_validate_structural_model(
    structural;
    symmetry_tolerance::Real=1.0e-12
)
    diagnostics = chang_structural_diagnostics(structural)
    println(
        "Wing-only modal frequencies (Hz): " *
        "$(round.(diagnostics.wing_modal_frequencies_hz,digits = 4))"
    )
    println(
        "Structural discretization: mass=$(structural.mass_model), " *
        "stiffness=$(structural.stiffness_model)"
    )
    println(
        "Structural checks: min eig(M)=$(diagnostics.minimum_mass_eigenvalue), " *
        "min eig(K)=$(diagnostics.minimum_stiffness_eigenvalue), " *
        "symmetry(M/K)=($(diagnostics.mass_symmetry_error), " *
        "$(diagnostics.stiffness_symmetry_error))"
    )

    diagnostics.minimum_mass_eigenvalue > 0.0 || error(
        "The free structural mass matrix is not positive definite"
    )
    diagnostics.minimum_stiffness_eigenvalue > 0.0 || error(
        "The free structural stiffness matrix is not positive definite"
    )
    diagnostics.mass_symmetry_error <= symmetry_tolerance || error(
        "The free structural mass matrix is not symmetric"
    )
    diagnostics.stiffness_symmetry_error <= symmetry_tolerance || error(
        "The free structural stiffness matrix is not symmetric"
    )
    return diagnostics
end

"""Allocate the structural history and enforce equilibrium at the initial state."""
function initialize_structural_history(
    structural,
    nStates::Int;
    initial_aerodynamic_load=nothing,
    initial_external_load=nothing
)
    displacement = Vector{Vector{Float64}}(undef,nStates)
    velocity = similar(displacement)
    acceleration = similar(displacement)

    displacement[1] = zeros(structural.ndof_free)
    velocity[1] = zeros(structural.ndof_free)
    initialForce = zeros(structural.ndof_free)
    !isnothing(initial_aerodynamic_load) && (initialForce .+= initial_aerodynamic_load)
    !isnothing(initial_external_load) && (initialForce .+= initial_external_load)
    acceleration[1] =
        isempty(structural.M) ? zeros(structural.ndof_free) :
        structural.M \ (
            initialForce .- structural.C * velocity[1] .-
            structural.K * displacement[1]
        )
    return (; displacement,velocity,acceleration)
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
    dofs_per_node::Int
)
    state = ones(structural.ndof_free)
    load = ones(structural.ndof_free)
    force = max(
        0.5 * reference.rho * freestream_speed^2 * reference_area,
        1.0
    )
    moment = max(force * reference_chord,1.0)

    # The root is clamped, so the reduced wing state contains nnodes - 1 nodes.
    for nodeIndex = 1:(wing_node_count-1)
        nodeStart = dofs_per_node * (nodeIndex - 1)
        state[nodeStart.+(1:3)] .= max(reference_chord,eps(Float64))
        state[nodeStart.+(4:6)] .= 1.0
        load[nodeStart.+(1:3)] .= force
        load[nodeStart.+(4:6)] .= moment
    end

    propellerMoment = max(
        0.5 * reference.rho * freestream_speed^2 * pi * propeller_radius^3,
        1.0
    )
    load[(structural.ndof_wing_free+1):end] .= propellerMoment
    return (; state,load)
end

# ### Accepted-step helpers

"""Print the generalized trim load captured before the pitch impulse."""
function report_trim_load(
    trim_load,
    sample_count::Int,
    captureTime::Real,
    wing_dof_count::Int,
    number_of_propellers::Int
)
    println(
        "\n=== Mean trim generalized load captured at " *
        "t=$(round(captureTime,digits = 4)) s from $sample_count samples ==="
    )
    for propeller_index = 1:number_of_propellers
        pitch_index = wing_dof_count + 2 * (propeller_index - 1) + 1
        yawIndex = pitch_index + 1
        println(
            "  P$propeller_index: pitch F0 = " *
            "$(round(trim_load[pitch_index],digits = 3)) N*m, " *
            "yaw F0 = $(round(trim_load[yawIndex],digits = 3)) N*m"
        )
    end
    println(
        "  |F0_wing| = $(round(norm(trim_load[1:wing_dof_count]),digits = 2)), " *
        "|F0_prop| = $(round(norm(trim_load[(wing_dof_count + 1):end]),digits = 2))"
    )
    return nothing
end

"""Advance and activate the wake exactly once after accepting a coupled step."""
function commit_chang_wake!(system,freestream,time_step,wake)
    advance_wake!(
        system,
        freestream,
        time_step;
        additional_velocity=nothing,
        repeated_points=wake.repeated_points,
        nwake=wake.active_rows,
        interaction_id=wake.interaction_ids,
        interaction=wake.interaction_on
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
    propeller_angle_limit_deg::Real
)
    maximum_state = maximum(abs,state)
    propellerAnglesDeg = rad2deg.(state[(wing_dof_count+1):end])
    maximum_propeller_angle_deg = maximum(abs,propellerAnglesDeg; init=0.0)
    nonfinite = any(value -> !isfinite(value),state)
    stateLimitExceeded = maximum_state > state_norm_limit
    angleLimitExceeded = maximum_propeller_angle_deg > propeller_angle_limit_deg
    reason =
        nonfinite ? "non-finite structural state" :
        stateLimitExceeded ? "state limit $state_norm_limit" :
        angleLimitExceeded ? "propeller-angle limit $propeller_angle_limit_deg deg" : ""
    return (;
        exceeded=nonfinite || stateLimitExceeded || angleLimitExceeded,
        reason,
        maximum_state,
        maximum_propeller_angle_deg
    )
end

# ### Coupled solution

"""
    solve_chang_aeroelastic!(system, structural; kwargs...)

March the selected structural integrator and UVLM coupling scheme. Implicit
trials are evaluated from the same beginning-of-step snapshot. Loose coupling
evaluates one target-time aerodynamic load with the structural state lagged at
`t_n`. A failed step is restored and rejected; every accepted step advances the
free wake exactly once.
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
    animation
)
    # This is a perturbation formulation: before the trim baseline and pulse,
    # both structural right-hand-side components are zero. The initializer
    # solves M*a0 = F0 - C*v0 - K*u0 rather than assuming acceleration is zero.
    history = initialize_structural_history(structural,length(time))
    displacement = history.displacement
    velocity = history.velocity
    acceleration = history.acceleration

    impulseDOFIndices =
        hasproperty(excitation,:dof_indices) ? excitation.dof_indices :
        [
            structural.ndof_wing_free + 2 * (index - 1) + 1
            for index in excitation.propeller_indices
        ]
    trim_load = zeros(structural.ndof_free)
    trimLoadAccumulator = zeros(structural.ndof_free)
    trimSampleCount = 0
    haveTrimLoad = false
    impulseWasReported = false

    coupling_iterations = fill(0,length(time_steps))
    coupling_state_residual = fill(NaN,length(time_steps))
    coupling_load_residual = fill(NaN,length(time_steps))
    coupling_equilibrium_residual = fill(NaN,length(time_steps))
    coupling_converged = fill(false,length(time_steps))
    acceptedPerturbationLoad = zeros(structural.ndof_free)
    lastAcceptedStep = length(time_steps)

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
            time[1]
        )
    end

    println("Starting coupled aeroelastic simulation...")
    for step in eachindex(time_steps)
        # All fixed-point trials start from this exact accepted aerodynamic state.
        copy_surfaces_to_previous!(system,wake.surface_count)
        aerodynamic_snapshot = snapshot_uvlm(system)
        time_step = time_steps[step]

        external_load_n = smooth_hann_pulse_load(
            time[step],
            structural.ndof_free,
            impulseDOFIndices;
            magnitude=excitation.magnitude,
            start_time=excitation.start_time,
            duration=excitation.duration
        )
        external_load_np1 = smooth_hann_pulse_load(
            time[step+1],
            structural.ndof_free,
            impulseDOFIndices;
            magnitude=excitation.magnitude,
            start_time=excitation.start_time,
            duration=excitation.duration
        )
        if !impulseWasReported &&
           (maximum(abs,external_load_n) > 0.0 || maximum(abs,external_load_np1) > 0.0)
            println("<<<<<<<<<<< Applying Smooth Moment Impulse >>>>>>>>>>>")
            impulseWasReported = true
        end

        # The coupling scheme changes only the exchange sequence. Both paths
        # use the same target-time geometry/load operator, structural equation,
        # aerodynamic snapshot, and single accepted wake commit.
        fullAerodynamicLoad = zeros(structural.ndof_free)
        load_for_state =
            stateGuess -> begin
                fullAerodynamicLoad .= aerodynamic_load(
                    aerodynamic_snapshot,
                    stateGuess,
                    step
                )
                return haveTrimLoad ? fullAerodynamicLoad .- trim_load :
                       zeros(structural.ndof_free)
            end
        stepSolution = if integration.coupling_scheme == :loose_explicit
            # One UVLM evaluation at target time/rotor azimuth t[n+1], using
            # the known structural state U[n]. Thus the structural forcing is
            # F_aero,n+1^explicit(U[n]), with no additional one-step delay.
            explicit_load = load_for_state(displacement[step])
            loose_explicit_aeroelastic_step(
                integration.time_integrator,
                structural.M,
                structural.C,
                structural.K,
                displacement[step],
                velocity[step],
                acceleration[step],
                explicit_load,
                acceptedPerturbationLoad,
                external_load_np1,
                external_load_n,
                time_step,
                integration.structural_parameters;
                options=integration.coupling,
                load_scale=integration.load_scale
            )
        elseif integration.coupling_scheme == :implicit_predictor_corrector
            diagnostic =
                integration.coupling_verbose ?
                values -> println(
                    "  step=$step, t=$(time[step + 1]) s, coupling_iter=$(values.iteration), " *
                    "state_res=$(values.state_residual), load_res=$(values.load_residual), " *
                    "coupled_eq_res=$(values.coupled_equilibrium_residual)"
                ) : nothing
            partitioned_aeroelastic_step(
                integration.time_integrator,
                structural.M,
                structural.C,
                structural.K,
                displacement[step],
                velocity[step],
                acceleration[step],
                acceptedPerturbationLoad,
                external_load_np1,
                external_load_n,
                time_step,
                integration.structural_parameters,
                load_for_state;
                options=integration.coupling,
                require_load_convergence=haveTrimLoad,
                state_scale=integration.state_scale,
                load_scale=integration.load_scale,
                diagnostic
            )
        else
            # Configuration validation should make this unreachable, but keep
            # the time marcher safe for programmatically constructed options.
            error("Unknown aeroelastic coupling scheme: $(integration.coupling_scheme)")
        end

        if !stepSolution.converged
            restore_uvlm!(system,aerodynamic_snapshot)
            error(
                "$(integration.coupling_scheme) coupling failed at step $step " *
                "(t=$(time[step + 1]) s) after $(stepSolution.iterations) iterations: " *
                "state_res=$(stepSolution.state_residual), " *
                "load_res=$(stepSolution.load_residual), " *
                "coupled_equilibrium_res=$(stepSolution.coupled_equilibrium_residual), " *
                "linear_equilibrium_res=$(stepSolution.equilibrium_residual). " *
                "The time step was not committed."
            )
        end

        displacement[step+1] = copy(stepSolution.displacement)
        velocity[step+1] = copy(stepSolution.velocity)
        acceleration[step+1] = copy(stepSolution.acceleration)
        acceptedPerturbationLoad .= stepSolution.trial_load

        # Average the periodic full load, then subtract it from later trials so
        # the structural solution contains aerodynamic perturbations only.
        acceptedTime = time[step+1]
        trimSample =
            !haveTrimLoad &&
            acceptedTime >= excitation.trim_average_start_time - eps(acceptedTime) &&
            acceptedTime <= excitation.start_time + eps(acceptedTime)
        if trimSample
            trimLoadAccumulator .+= fullAerodynamicLoad
            trimSampleCount += 1
        end
        if !haveTrimLoad && acceptedTime >= excitation.start_time - eps(acceptedTime)
            trimSampleCount > 0 || error(
                "No aerodynamic samples were available for trim averaging"
            )
            trim_load .= trimLoadAccumulator ./ trimSampleCount
            haveTrimLoad = true
            acceptedPerturbationLoad .= 0.0
            report_trim_load(
                trim_load,
                trimSampleCount,
                acceptedTime,
                structural.ndof_wing_free,
                excitation.number_of_propellers
            )
        end

        coupling_iterations[step] = stepSolution.iterations
        coupling_state_residual[step] = stepSolution.state_residual
        coupling_load_residual[step] = stepSolution.load_residual
        coupling_equilibrium_residual[step] = stepSolution.coupled_equilibrium_residual
        coupling_converged[step] = stepSolution.converged

        # The converged circulation/load trial is already in `system`; only now
        # is its wake convected and the newly shed row activated.
        commit_chang_wake!(system,freestream_history[step],time_step,wake)

        if animation.enabled &&
           (step == 1 || step % animation.stride == 0 || step == length(time_steps))
            animation.record_frame(
                animation_surface_history,
                animation_wake_history,
                animation_active_wake_rows_history,
                animation_time_history,
                acceptedTime
            )
        end

        couplingDescription =
            integration.coupling_scheme == :loose_explicit ?
            "1 explicit exchange" : "$(stepSolution.iterations) implicit iterations"
        println(
            "Step $step/$(length(time_steps)) (t=$(round(acceptedTime,digits = 4)) s): " *
            "$couplingDescription, " *
            "state_res=$(round(stepSolution.state_residual,sigdigits = 4)), " *
            "load_res=$(round(stepSolution.load_residual,sigdigits = 4)), " *
            "coupled_eq_res=$(round(stepSolution.coupled_equilibrium_residual,sigdigits = 4)), " *
            "linear_eq_res=$(round(stepSolution.equilibrium_residual,sigdigits = 4))"
        )
        flush(stdout)

        if step in wake.saved_steps
            wake.surface_history[step] = [copy(surface) for surface in system.surfaces]
        end

        limit = chang_state_limit_status(
            displacement[step+1],
            structural.ndof_wing_free;
            state_norm_limit=integration.state_norm_limit,
            propeller_angle_limit_deg=integration.propeller_angle_limit_deg
        )
        if limit.exceeded
            sampleRange = max(1,step - 400):step
            pitchEnvelope = [
                abs(displacement[index][impulseDOFIndices[1]])
                for index in sampleRange
            ]
            sampleTime = [time[index] for index in sampleRange]
            positive_indices = findall(>(1.0e-12),pitchEnvelope)
            message =
                "\n>>> ABORT at t=$(round(acceptedTime,digits = 4)) s, " *
                "reason=$(limit.reason), " *
                "|U|=$(round(limit.maximum_state,sigdigits = 4)), " *
                "max propeller angle=" *
                "$(round(limit.maximum_propeller_angle_deg,sigdigits = 4)) deg"
            if length(positive_indices) > 5
                firstIndex,lastIndex = first(positive_indices),last(positive_indices)
                growth_rate =
                    (
                        log(pitchEnvelope[lastIndex]) - log(pitchEnvelope[firstIndex])
                    ) / (sampleTime[lastIndex] - sampleTime[firstIndex])
                message *= "; growth sigma(excited coordinate) approx " *
                           "$(round(growth_rate,digits = 4)) /s"
            end
            println(message)
            lastAcceptedStep = step
            break
        end
    end
    println("Simulation finished.")

    if animation.enabled && last(animation_time_history) != time[lastAcceptedStep+1]
        animation.record_frame(
            animation_surface_history,
            animation_wake_history,
            animation_active_wake_rows_history,
            animation_time_history,
            time[lastAcceptedStep+1]
        )
    end

    return (;
        time_integrator=integration.time_integrator,
        coupling_scheme=integration.coupling_scheme,
        displacement_history=displacement,
        velocity_history=velocity,
        acceleration_history=acceleration,
        last_step=lastAcceptedStep,
        trim_load,
        coupling_iterations,
        coupling_state_residual,
        coupling_load_residual,
        coupling_equilibrium_residual,
        coupling_converged,
        animation_surface_history,
        animation_wake_history,
        animation_active_wake_rows_history,
        animation_time_history
    )
end
