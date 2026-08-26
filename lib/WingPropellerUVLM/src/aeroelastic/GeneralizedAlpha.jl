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
    generalized_alpha_kinematics(displacement_np1, displacement_n,
        velocity_n, acceleration_n, dt, parameters)

Recover the end-of-step velocity and acceleration associated with a prescribed
end-of-step displacement using the Newmark relations of generalized-alpha.
This is used when a partitioned iteration accepts the aerodynamic evaluation
state rather than the result of one additional structural solve.
"""
function generalized_alpha_kinematics(
    displacement_np1,
    displacement_n,
    velocity_n,
    acceleration_n,
    dt::Real,
    parameters,
)
    dt > 0 || throw(ArgumentError("The generalized-alpha time step must be positive"))
    beta = parameters.beta
    gamma = parameters.gamma
    a0 = 1.0 / (beta * dt^2)
    a2 = 1.0 / (beta * dt)
    a3 = 1.0 / (2.0 * beta) - 1.0
    acceleration = a0 .* (displacement_np1 .- displacement_n) .-
        a2 .* velocity_n .- a3 .* acceleration_n
    velocity = velocity_n .+ dt .* (1.0 - gamma) .* acceleration_n .+
        gamma .* dt .* acceleration
    return (; velocity, acceleration)
end

function _validate_residual_scale(scale, expected_length::Int, name::AbstractString)
    isnothing(scale) && return nothing
    length(scale) == expected_length || throw(DimensionMismatch(
        "$name must contain $expected_length entries",
    ))
    all(isfinite, scale) || throw(ArgumentError("$name must contain only finite values"))
    all(>(0), scale) || throw(ArgumentError("$name entries must be positive"))
    return scale
end

function _scaled_relative_residual(delta, reference, scale)
    validated_scale = _validate_residual_scale(scale, length(delta), "residual scale")
    if isnothing(validated_scale)
        return norm(delta, Inf) / max(norm(reference, Inf), 1.0)
    end
    scaled_delta = delta ./ validated_scale
    scaled_reference = reference ./ validated_scale
    return norm(scaled_delta, Inf) / max(norm(scaled_reference, Inf), 1.0)
end

"""
    generalized_alpha_equilibrium_residual(...; load_scale=nothing)

Return the normalized residual of the complete generalized-alpha equilibrium
equation for a proposed end-of-step state and load. Unlike the linear-solver
residual returned by `generalized_alpha_corrector`, this quantity detects a
mismatch between the structural state and the aerodynamic load used by the
partitioned iteration.
"""
function generalized_alpha_equilibrium_residual(
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
    alpha_m = parameters.alpha_m
    alpha_f = parameters.alpha_f
    acceleration_alpha = (1.0 - alpha_m) .* acceleration_np1 .+
        alpha_m .* acceleration_n
    velocity_alpha = (1.0 - alpha_f) .* velocity_np1 .+ alpha_f .* velocity_n
    displacement_alpha = (1.0 - alpha_f) .* displacement_np1 .+
        alpha_f .* displacement_n
    force_alpha = (1.0 - alpha_f) .* (force_np1 .+ external_np1) .+
        alpha_f .* (force_n .+ external_n)
    residual = M * acceleration_alpha .+ C * velocity_alpha .+
        K * displacement_alpha .- force_alpha
    return _scaled_relative_residual(residual, force_alpha, load_scale)
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
    kinematics = generalized_alpha_kinematics(displacement, u_n, v_n, a_n, dt, parameters)
    acceleration = kinematics.acceleration
    velocity = kinematics.velocity

    return (; displacement, velocity, acceleration, equilibrium_residual)
end

Base.@kwdef struct PartitionedCouplingOptions{T<:Real}
    maximum_iterations::Int = 10
    state_tolerance::T = 1.0e-5
    load_tolerance::T = 1.0e-2
    equilibrium_tolerance::T = 1.0e-10
    coupled_equilibrium_tolerance::T = 1.0e-5
    relaxation::T = 0.5
end

function validate_partitioned_coupling_options(options::PartitionedCouplingOptions)
    options.maximum_iterations > 0 ||
        throw(ArgumentError("maximum_iterations must be positive"))
    options.state_tolerance > 0 || throw(ArgumentError("state_tolerance must be positive"))
    options.load_tolerance > 0 || throw(ArgumentError("load_tolerance must be positive"))
    options.equilibrium_tolerance > 0 ||
        throw(ArgumentError("equilibrium_tolerance must be positive"))
    options.coupled_equilibrium_tolerance > 0 ||
        throw(ArgumentError("coupled_equilibrium_tolerance must be positive"))
    0.0 < options.relaxation <= 1.0 ||
        throw(ArgumentError("relaxation must lie in (0, 1]"))
    return options
end

"""
    partitioned_generalized_alpha_step(..., load_at_state;
        options, require_load_convergence)

Perform the fixed-point aeroelastic correction for one physical time step.
`load_at_state(state)` must restore the beginning-of-step aerodynamic snapshot,
evaluate the generalized load at `state`, and leave that trial as the current
aerodynamic state. On convergence, the returned displacement is exactly the
last state passed to the callback, so the structural and aerodynamic states
form one accepted fixed-point pair. If convergence fails, the caller must
restore its snapshot and must not advance the physical time step.
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
    state_scale = nothing,
    load_scale = nothing,
)
    validate_partitioned_coupling_options(options)
    _validate_residual_scale(state_scale, length(displacement_n), "state_scale")
    _validate_residual_scale(load_scale, length(load_n), "load_scale")

    predictor = displacement_n .+ dt .* velocity_n .+
        dt^2 .* (0.5 - parameters.beta) .* acceleration_n
    state_guess = copy(predictor)
    accepted_displacement = copy(state_guess)
    accepted_velocity = copy(velocity_n)
    accepted_acceleration = copy(acceleration_n)
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
        candidate_load = load_at_state(state_guess)
        length(candidate_load) == length(load_n) || throw(DimensionMismatch(
            "load_at_state returned $(length(candidate_load)) entries; expected $(length(load_n))",
        ))
        accepted_load .= candidate_load
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
        equilibrium_residual = correction.equilibrium_residual

        accepted_kinematics = generalized_alpha_kinematics(
            accepted_displacement,
            displacement_n,
            velocity_n,
            acceleration_n,
            dt,
            parameters,
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
        coupled_equilibrium_residual = generalized_alpha_equilibrium_residual(
            M,
            C,
            K,
            accepted_displacement,
            accepted_velocity,
            accepted_acceleration,
            displacement_n,
            velocity_n,
            acceleration_n,
            candidate_load,
            load_n,
            external_load_np1,
            external_load_n,
            parameters;
            load_scale,
        )
        load_is_converged = !require_load_convergence ||
            (iteration > 1 && load_residual <= options.load_tolerance)

        if state_residual <= options.state_tolerance &&
            equilibrium_residual <= options.equilibrium_tolerance &&
            coupled_equilibrium_residual <= options.coupled_equilibrium_tolerance &&
            load_is_converged
            converged = true
            break
        end

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
