const STRUCTURAL_TIME_INTEGRATORS = (:newmark_beta, :generalized_alpha)
const AEROELASTIC_COUPLING_SCHEMES = (:loose_explicit, :implicit_predictor_corrector)

"""Reject an unsupported structural integration method with a user-facing error."""
function validate_structural_time_integrator(method::Symbol)
    method in STRUCTURAL_TIME_INTEGRATORS ||
        error("Unknown structural time integration method: $method")
    return method
end

"""Reject an unsupported aeroelastic coupling strategy with a user-facing error."""
function validate_aeroelastic_coupling_scheme(scheme::Symbol)
    scheme in AEROELASTIC_COUPLING_SCHEMES ||
        error("Unknown aeroelastic coupling scheme: $scheme")
    return scheme
end

"""
    newmark_beta_parameters(alpha=0.05)

Return the reference solver's dissipative Newmark-beta parameters. This is
intentionally not the classical average-acceleration pair: with `alpha=0.05`,
`gamma=0.55` and `beta=0.275625`.
"""
function newmark_beta_parameters(alpha::Real = 0.05)
    alpha_nb = Float64(alpha)
    isfinite(alpha_nb) || throw(ArgumentError("newmark alpha must be finite"))
    gamma = 0.5 + alpha_nb
    beta = 0.25 * (gamma + 0.5)^2
    gamma > 0.0 || throw(ArgumentError("newmark gamma must be positive"))
    beta > 0.0 || throw(ArgumentError("newmark beta must be positive"))
    return (; alpha = alpha_nb, alpha_nb, gamma, beta)
end

"""Calculate all Newmark coefficients from the current physical time step."""
function newmark_beta_coefficients(dt::Real, parameters)
    dt > 0 || throw(ArgumentError("The Newmark-beta time step must be positive"))
    beta = parameters.beta
    gamma = parameters.gamma
    a0 = 1.0 / (beta * dt^2)
    a1 = gamma / (beta * dt)
    a2 = 1.0 / (beta * dt)
    a3 = 1.0 / (2.0 * beta) - 1.0
    a4 = gamma / beta - 1.0
    a5 = (dt / 2.0) * (a4 - 1.0)
    a6 = dt * (1.0 - gamma)
    a7 = gamma * dt
    return (; a0, a1, a2, a3, a4, a5, a6, a7)
end

"""Recover the Newmark end-of-step velocity and acceleration from displacement."""
function newmark_beta_kinematics(
    displacement_np1,
    displacement_n,
    velocity_n,
    acceleration_n,
    dt::Real,
    parameters,
)
    coefficients = newmark_beta_coefficients(dt, parameters)
    (; a0, a2, a3, a6, a7) = coefficients
    acceleration = a0 .* (displacement_np1 .- displacement_n) .-
        a2 .* velocity_n .- a3 .* acceleration_n
    velocity = velocity_n .+ a6 .* acceleration_n .+ a7 .* acceleration
    return (; velocity, acceleration)
end

"""
    newmark_beta_corrector(M, C, K, u_n, v_n, a_n,
        force_np1, external_np1, dt, parameters)

Solve one Newmark-beta structural correction using the force at `n+1`.
`M`, `C`, and `K` are held fixed by the caller throughout any coupling
subiterations.
"""
function newmark_beta_corrector(
    M,
    C,
    K,
    displacement_n,
    velocity_n,
    acceleration_n,
    force_np1,
    external_np1,
    dt::Real,
    parameters;
    effective_stiffness = nothing,
    effective_factorization = nothing,
)
    coefficients = newmark_beta_coefficients(dt, parameters)
    (; a0, a1, a2, a3, a4, a5, a6, a7) = coefficients

    # Newmark rewrites M*a[n+1] + C*v[n+1] + K*u[n+1] = F[n+1]
    # as one displacement solve with the constant effective matrix below.
    if isnothing(effective_stiffness)
        effective_stiffness = K .+ a0 .* M .+ a1 .* C
    end
    total_force = force_np1 .+ external_np1
    effective_force = total_force .+
        M * (a0 .* displacement_n .+ a2 .* velocity_n .+ a3 .* acceleration_n) .+
        C * (a1 .* displacement_n .+ a4 .* velocity_n .+ a5 .* acceleration_n)

    effective_solver = isnothing(effective_factorization) ?
        effective_stiffness : effective_factorization
    displacement = effective_solver \ effective_force
    equilibrium_residual = norm(effective_stiffness * displacement .- effective_force) /
        max(norm(effective_force), 1.0)
    # Recover acceleration and velocity from the corrected displacement so the
    # three structural fields always satisfy the same Newmark kinematics.
    acceleration = a0 .* (displacement .- displacement_n) .-
        a2 .* velocity_n .- a3 .* acceleration_n
    velocity = velocity_n .+ a6 .* acceleration_n .+ a7 .* acceleration

    return (; displacement, velocity, acceleration, equilibrium_residual)
end

"""Normalized physical equilibrium residual for a Newmark end state."""
function newmark_beta_equilibrium_residual(
    M,
    C,
    K,
    displacement_np1,
    velocity_np1,
    acceleration_np1,
    force_np1,
    external_np1;
    load_scale = nothing,
)
    total_force = force_np1 .+ external_np1
    residual = M * acceleration_np1 .+ C * velocity_np1 .+
        K * displacement_np1 .- total_force
    return _scaled_relative_residual(residual, total_force, load_scale)
end

"""Return the selected integrator's parameters without coupling it to a scheme."""
function structural_integration_parameters(
    method::Symbol;
    newmark_alpha::Real = 0.05,
    generalized_alpha_rho_infinity::Real = 1.0,
)
    validate_structural_time_integrator(method)
    return method == :newmark_beta ? newmark_beta_parameters(newmark_alpha) :
        generalized_alpha_parameters(generalized_alpha_rho_infinity)
end

"""Predict the `n+1` state using the selected method's own beta and gamma."""
function structural_predictor(
    method::Symbol,
    displacement_n,
    velocity_n,
    acceleration_n,
    dt::Real,
    parameters,
)
    validate_structural_time_integrator(method)
    dt > 0 || throw(ArgumentError("The structural time step must be positive"))
    displacement = displacement_n .+ dt .* velocity_n .+
        dt^2 .* (0.5 - parameters.beta) .* acceleration_n
    velocity = velocity_n .+ dt .* (1.0 - parameters.gamma) .* acceleration_n
    acceleration = copy(acceleration_n)
    return (; displacement, velocity, acceleration)
end

function structural_kinematics(
    method::Symbol,
    displacement_np1,
    displacement_n,
    velocity_n,
    acceleration_n,
    dt::Real,
    parameters,
)
    validate_structural_time_integrator(method)
    if method == :newmark_beta
        return newmark_beta_kinematics(
            displacement_np1, displacement_n, velocity_n, acceleration_n, dt, parameters,
        )
    end
    return generalized_alpha_kinematics(
        displacement_np1, displacement_n, velocity_n, acceleration_n, dt, parameters,
    )
end

"""Build and factor the constant effective matrix for one physical step."""
function structural_effective_system(method::Symbol, M, C, K, dt::Real, parameters)
    validate_structural_time_integrator(method)
    dt > 0 || throw(ArgumentError("The structural time step must be positive"))
    if method == :newmark_beta
        coefficients = newmark_beta_coefficients(dt, parameters)
        matrix = K .+ coefficients.a0 .* M .+ coefficients.a1 .* C
    else
        beta = parameters.beta
        gamma = parameters.gamma
        alpha_m = parameters.alpha_m
        alpha_f = parameters.alpha_f
        a0 = 1.0 / (beta * dt^2)
        velocity_coefficient = gamma / (beta * dt)
        matrix = (1.0 - alpha_m) .* a0 .* M .+
            (1.0 - alpha_f) .* velocity_coefficient .* C .+
            (1.0 - alpha_f) .* K
    end
    return (; matrix, factorization = factorize(matrix))
end

function structural_corrector(
    method::Symbol,
    M,
    C,
    K,
    displacement_n,
    velocity_n,
    acceleration_n,
    force_np1,
    force_n,
    external_np1,
    external_n,
    dt::Real,
    parameters;
    effective_system = nothing,
)
    validate_structural_time_integrator(method)
    effective_keywords = isnothing(effective_system) ? NamedTuple() : (;
        effective_stiffness = effective_system.matrix,
        effective_factorization = effective_system.factorization,
    )
    if method == :newmark_beta
        return newmark_beta_corrector(
            M, C, K, displacement_n, velocity_n, acceleration_n,
            force_np1, external_np1, dt, parameters; effective_keywords...,
        )
    end
    return generalized_alpha_corrector(
        M, C, K, displacement_n, velocity_n, acceleration_n,
        force_np1, force_n, external_np1, external_n, dt, parameters;
        effective_keywords...,
    )
end

function structural_equilibrium_residual(
    method::Symbol,
    M,
    C,
    K,
    displacement_np1,
    velocity_np1,
    acceleration_np1,
    displacement_n,
    velocity_n,
    acceleration_n,
    force_np1,
    force_n,
    external_np1,
    external_n,
    parameters;
    load_scale = nothing,
)
    validate_structural_time_integrator(method)
    if method == :newmark_beta
        return newmark_beta_equilibrium_residual(
            M, C, K, displacement_np1, velocity_np1, acceleration_np1,
            force_np1, external_np1; load_scale,
        )
    end
    return generalized_alpha_equilibrium_residual(
        M, C, K, displacement_np1, velocity_np1, acceleration_np1,
        displacement_n, velocity_n, acceleration_n,
        force_np1, force_n, external_np1, external_n, parameters; load_scale,
    )
end

"""
    loose_explicit_aeroelastic_step(...)

Perform exactly one structural solve with an aerodynamic load that the caller
has evaluated exactly once from the lagged structural state. The UVLM state is
owned and committed by the caller, separately from this structural operation.
"""
function loose_explicit_aeroelastic_step(
    method::Symbol,
    M,
    C,
    K,
    displacement_n,
    velocity_n,
    acceleration_n,
    force_np1,
    force_n,
    external_np1,
    external_n,
    dt::Real,
    parameters;
    options::PartitionedCouplingOptions = PartitionedCouplingOptions(),
    load_scale = nothing,
)
    validate_partitioned_coupling_options(options)
    correction = structural_corrector(
        method, M, C, K, displacement_n, velocity_n, acceleration_n,
        force_np1, force_n, external_np1, external_n, dt, parameters,
    )
    coupled_equilibrium_residual = structural_equilibrium_residual(
        method, M, C, K,
        correction.displacement, correction.velocity, correction.acceleration,
        displacement_n, velocity_n, acceleration_n,
        force_np1, force_n, external_np1, external_n, parameters; load_scale,
    )
    converged = all(isfinite, correction.displacement) &&
        all(isfinite, correction.velocity) && all(isfinite, correction.acceleration) &&
        correction.equilibrium_residual <= options.equilibrium_tolerance &&
        coupled_equilibrium_residual <= options.coupled_equilibrium_tolerance
    return (;
        displacement = correction.displacement,
        velocity = correction.velocity,
        acceleration = correction.acceleration,
        trial_load = copy(force_np1),
        iterations = 1,
        state_residual = 0.0,
        load_residual = 0.0,
        equilibrium_residual = correction.equilibrium_residual,
        coupled_equilibrium_residual,
        converged,
    )
end

"""
    partitioned_aeroelastic_step(method, ..., load_at_state; ...)

Iterate aerodynamic loads and structural corrections at one physical `n+1`
time. `load_at_state` must restore the same beginning-of-step aerodynamic
snapshot before every call and must not advance the wake or physical time.
Previous physical structural states remain read-only until the returned state
is accepted by the caller.
"""
function partitioned_aeroelastic_step(
    method::Symbol,
    M,
    C,
    K,
    displacement_n,
    velocity_n,
    acceleration_n,
    load_n,
    external_load_np1,
    external_load_n,
    dt::Real,
    parameters,
    load_at_state;
    options::PartitionedCouplingOptions = PartitionedCouplingOptions(),
    require_load_convergence::Bool = true,
    state_scale = nothing,
    load_scale = nothing,
    diagnostic = nothing,
)
    validate_structural_time_integrator(method)
    validate_partitioned_coupling_options(options)
    _validate_residual_scale(state_scale, length(displacement_n), "state_scale")
    _validate_residual_scale(load_scale, length(load_n), "load_scale")
    effective_system = structural_effective_system(method, M, C, K, dt, parameters)

    # Start the coupling loop with a load-free structural prediction at n+1.
    predictor = structural_predictor(
        method, displacement_n, velocity_n, acceleration_n, dt, parameters,
    )
    state_guess = copy(predictor.displacement)
    accepted_displacement = copy(state_guess)
    accepted_velocity = copy(predictor.velocity)
    accepted_acceleration = copy(predictor.acceleration)
    accepted_load = copy(load_n)
    previous_load = nothing
    state_residual = Inf
    load_residual = Inf
    equilibrium_residual = Inf
    coupled_equilibrium_residual = Inf
    converged = false
    iterations = 0

    for iteration in 1:options.maximum_iterations
        iterations = iteration
        accepted_displacement .= state_guess
        # The callback maps this structural guess to UVLM geometry, solves one
        # aerodynamic trial from the n snapshot, and returns generalized loads.
        candidate_load = load_at_state(state_guess)
        length(candidate_load) == length(load_n) || throw(DimensionMismatch(
            "load_at_state returned $(length(candidate_load)) entries; expected $(length(load_n))",
        ))
        accepted_load .= candidate_load
        # Solve the structure with the current aerodynamic load at n+1.
        correction = structural_corrector(
            method, M, C, K, displacement_n, velocity_n, acceleration_n,
            candidate_load, load_n, external_load_np1, external_load_n, dt, parameters;
            effective_system,
        )
        equilibrium_residual = correction.equilibrium_residual

        accepted_kinematics = structural_kinematics(
            method, accepted_displacement, displacement_n, velocity_n,
            acceleration_n, dt, parameters,
        )
        accepted_velocity .= accepted_kinematics.velocity
        accepted_acceleration .= accepted_kinematics.acceleration
        state_residual = _scaled_relative_residual(
            correction.displacement .- accepted_displacement,
            accepted_displacement,
            state_scale,
        )
        load_residual = !require_load_convergence ? 0.0 :
            (isnothing(previous_load) ? Inf : _scaled_relative_residual(
                candidate_load .- previous_load,
                candidate_load,
                load_scale,
            ))
        coupled_equilibrium_residual = structural_equilibrium_residual(
            method, M, C, K,
            accepted_displacement, accepted_velocity, accepted_acceleration,
            displacement_n, velocity_n, acceleration_n,
            candidate_load, load_n, external_load_np1, external_load_n,
            parameters; load_scale,
        )
        load_is_converged = !require_load_convergence ||
            (iteration > 1 && load_residual <= options.load_tolerance)

        if !isnothing(diagnostic)
            diagnostic((;
                iteration,
                state_residual,
                load_residual,
                equilibrium_residual,
                coupled_equilibrium_residual,
            ))
        end

        # Accept only when the exchanged state/load and the structural equation
        # all meet their tolerances for this physical time step.
        if state_residual <= options.state_tolerance &&
            equilibrium_residual <= options.equilibrium_tolerance &&
            coupled_equilibrium_residual <= options.coupled_equilibrium_tolerance &&
            load_is_converged
            converged = true
            break
        end

        # Velocity and acceleration are not relaxed independently: for both
        # methods they are recovered from this relaxed displacement using the
        # method's exact kinematic relation at the next trial.
        state_guess .= options.relaxation .* correction.displacement .+
            (1.0 - options.relaxation) .* state_guess
        previous_load = copy(candidate_load)
    end

    return (;
        displacement = accepted_displacement,
        velocity = accepted_velocity,
        acceleration = accepted_acceleration,
        trial_load = accepted_load,
        iterations,
        state_residual,
        load_residual,
        equilibrium_residual,
        coupled_equilibrium_residual,
        converged,
    )
end
