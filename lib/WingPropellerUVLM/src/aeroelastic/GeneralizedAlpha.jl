"""
    generalized_alpha_parameters(rho_infinity)

Return the Chung-Hulbert generalized-alpha parameters for a second-order
structural system. `rho_infinity` is the desired high-frequency spectral
radius and must lie in `[0, 1]`.
"""
function generalized_alpha_parameters(rho_infinity::Real)
    rho = Float64(rho_infinity)
    0.0 <= rho <= 1.0 || throw(ArgumentError("rho_infinity must be in [0, 1]"))

    alpha_m = (2.0 * rho - 1.0) / (rho + 1.0)
    alpha_f = rho / (rho + 1.0)
    gamma = 0.5 - alpha_m + alpha_f
    beta = 0.25 * (1.0 - alpha_m + alpha_f)^2
    return (; rho_infinity = rho, rho_inf = rho, alpha_m, alpha_f, gamma, beta)
end

"""
    generalized_alpha_corrector(M, C, K, u_n, v_n, a_n,
        force_np1, force_n, external_np1, external_n, dt, parameters)

Perform one linear generalized-alpha structural correction. The returned
named tuple contains the displacement, velocity, acceleration, and normalized
effective-equilibrium residual at `n + 1`.
"""
function generalized_alpha_corrector(M, C, K, u_n, v_n, a_n,
    force_np1, force_n, external_np1, external_n, dt::Real, parameters)

    dt > 0 || throw(ArgumentError("The generalized-alpha time step must be positive"))
    beta = parameters.beta
    gamma = parameters.gamma
    alpha_m = parameters.alpha_m
    alpha_f = parameters.alpha_f

    a0 = 1.0 / (beta * dt^2)
    a2 = 1.0 / (beta * dt)
    a3 = 1.0 / (2.0 * beta) - 1.0
    a6 = dt * (1.0 - gamma)
    a7 = gamma * dt
    velocity_coefficient = a7 * a0

    constant_acceleration = -a0 .* u_n .- a2 .* v_n .- a3 .* a_n
    constant_velocity = v_n .+ a6 .* a_n .+ a7 .* constant_acceleration
    force_alpha = (1.0 - alpha_f) .* (force_np1 .+ external_np1) .+
        alpha_f .* (force_n .+ external_n)

    effective_stiffness = (1.0 - alpha_m) .* a0 .* M .+
        (1.0 - alpha_f) .* velocity_coefficient .* C .+
        (1.0 - alpha_f) .* K
    right_hand_side = force_alpha .-
        M * ((1.0 - alpha_m) .* constant_acceleration .+ alpha_m .* a_n) .-
        C * ((1.0 - alpha_f) .* constant_velocity .+ alpha_f .* v_n) .-
        K * (alpha_f .* u_n)

    displacement = effective_stiffness \ right_hand_side
    equilibrium_residual = norm(effective_stiffness * displacement .- right_hand_side) /
        max(norm(right_hand_side), 1.0)
    acceleration = a0 .* (displacement .- u_n) .- a2 .* v_n .- a3 .* a_n
    velocity = v_n .+ a6 .* a_n .+ a7 .* acceleration

    return (; displacement, velocity, acceleration, equilibrium_residual)
end

Base.@kwdef struct PartitionedCouplingOptions{T<:Real}
    maximum_iterations::Int = 10
    state_tolerance::T = 1.0e-5
    load_tolerance::T = 1.0e-2
    equilibrium_tolerance::T = 1.0e-10
    relaxation::T = 0.5
end

function validate_partitioned_coupling_options(options::PartitionedCouplingOptions)
    options.maximum_iterations > 0 ||
        throw(ArgumentError("maximum_iterations must be positive"))
    options.state_tolerance > 0 || throw(ArgumentError("state_tolerance must be positive"))
    options.load_tolerance > 0 || throw(ArgumentError("load_tolerance must be positive"))
    options.equilibrium_tolerance > 0 ||
        throw(ArgumentError("equilibrium_tolerance must be positive"))
    0.0 < options.relaxation <= 1.0 ||
        throw(ArgumentError("relaxation must lie in (0, 1]"))
    return options
end

"""
    partitioned_generalized_alpha_step(..., load_at_state;
        options, require_load_convergence)

Perform the fixed-point aeroelastic correction for one physical time step.
`load_at_state(state)` must evaluate and return the aerodynamic generalized
load at the guessed end-of-step structural state. Trial-state restoration and
commit semantics remain the responsibility of that callback and its caller.
"""
function partitioned_generalized_alpha_step(
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
)
    validate_partitioned_coupling_options(options)

    predictor = displacement_n .+ dt .* velocity_n .+
        dt^2 .* (0.5 - parameters.beta) .* acceleration_n
    state_guess = copy(predictor)
    corrected_displacement = copy(state_guess)
    corrected_velocity = copy(velocity_n)
    corrected_acceleration = copy(acceleration_n)
    load_guess = copy(load_n)
    state_residual = Inf
    load_residual = Inf
    equilibrium_residual = Inf
    converged = false
    iterations = 0

    for iteration in 1:options.maximum_iterations
        iterations = iteration
        candidate_load = load_at_state(state_guess)
        correction = generalized_alpha_corrector(
            M,
            C,
            K,
            displacement_n,
            velocity_n,
            acceleration_n,
            candidate_load,
            load_n,
            external_load_np1,
            external_load_n,
            dt,
            parameters,
        )
        corrected_displacement = correction.displacement
        corrected_velocity = correction.velocity
        corrected_acceleration = correction.acceleration
        equilibrium_residual = correction.equilibrium_residual

        state_residual = norm(corrected_displacement .- state_guess) /
            max(norm(corrected_displacement), 1.0)
        load_residual = norm(candidate_load .- load_guess) /
            max(norm(candidate_load), 1.0)
        load_is_converged = !require_load_convergence ||
            (iteration > 1 && load_residual <= options.load_tolerance)

        if state_residual <= options.state_tolerance &&
            equilibrium_residual <= options.equilibrium_tolerance &&
            load_is_converged
            converged = true
            state_guess .= corrected_displacement
            load_guess .= candidate_load
            break
        end

        state_guess .= options.relaxation .* corrected_displacement .+
            (1.0 - options.relaxation) .* state_guess
        load_guess .= candidate_load
    end

    return (;
        displacement = corrected_displacement,
        velocity = corrected_velocity,
        acceleration = corrected_acceleration,
        trial_load = load_guess,
        iterations,
        state_residual,
        load_residual,
        equilibrium_residual,
        converged,
    )
end
