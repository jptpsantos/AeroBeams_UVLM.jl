module PazyLinearUVLM

let examples = normpath(joinpath(@__DIR__,".."))
    examples in LOAD_PATH || pushfirst!(LOAD_PATH,examples)
end

using LinearAlgebra
import WingPropellerUVLM as UVLM

export PazyUVLMStateSpace,linearize_uvlm_step

struct PazyUVLMStateSpace
    A::Matrix{Float64}
    B::Matrix{Float64}
    C::Matrix{Float64}
    D::Matrix{Float64}
    state_reference::Vector{Float64}
    force_reference::Vector{Float64}
    time_step::Float64
    relative_state_drift::Float64
    relative_force_drift::Float64
end

function pack_state(system)
    state = copy(system.Γ)
    for surface in eachindex(system.wakes),row = 1:system.nwake[surface],
        span_panel in axes(system.wakes[surface],2)
        push!(state,UVLM.circulation_strength(system.wakes[surface][row,span_panel]))
    end
    return state
end

function unpack_state!(system,state)
    bound_count = length(system.Γ)
    wake_count = sum(system.nwake[surface] * size(system.wakes[surface],2)
        for surface in eachindex(system.wakes))
    length(state) == bound_count + wake_count || throw(DimensionMismatch(
        "Linear UVLM state length does not match bound and active wake panels"))
    copyto!(system.Γ,view(state,1:bound_count))
    index = bound_count
    for surface in eachindex(system.wakes),row = 1:system.nwake[surface],
        span_panel in axes(system.wakes[surface],2)
        index += 1
        panel = system.wakes[surface][row,span_panel]
        system.wakes[surface][row,span_panel] = UVLM.set_circulation_strength(
            panel,state[index])
    end
    return system
end

function pack_forces(loads)
    force_count = sum(3 * length(surface) for surface in loads.forces)
    forces = Vector{Float64}(undef,force_count)
    index = 0
    for surface in loads.forces,force in surface,component in force
        index += 1
        forces[index] = component
    end
    return forces
end

function advance_wake_circulation!(system)
    bound_offset = 0
    for surface in eachindex(system.wakes)
        rows = system.nwake[surface]
        span_panels = size(system.wakes[surface],2)
        panel_indices = LinearIndices(size(system.surfaces[surface]))
        surface_panels = size(system.surfaces[surface],1)

        for span_panel in 1:span_panels
            previous = [UVLM.circulation_strength(system.wakes[surface][row,span_panel])
                for row in 1:rows]
            for row in rows:-1:2
                panel = system.wakes[surface][row,span_panel]
                system.wakes[surface][row,span_panel] = UVLM.set_circulation_strength(
                    panel,previous[row-1])
            end
            panel = system.wakes[surface][1,span_panel]
            trailing_edge_index = bound_offset + panel_indices[surface_panels,span_panel]
            system.wakes[surface][1,span_panel] = UVLM.set_circulation_strength(
                panel,system.Γ[trailing_edge_index])
        end
        bound_offset += length(system.surfaces[surface])
    end
    return system
end

function validate_directions(directions,grids,label)
    directions isa AbstractVector || throw(ArgumentError("$label must be a vector of grid directions"))
    for direction in directions
        length(direction) == length(grids) || throw(DimensionMismatch(
            "$label surface count differs from the UVLM reference"))
        for surface in eachindex(grids)
            size(direction[surface]) == size(grids[surface]) || throw(DimensionMismatch(
                "$label grid direction does not match the UVLM reference"))
            all(isfinite,direction[surface]) || throw(ArgumentError(
                "$label contains nonfinite values"))
        end
    end
    return length(directions)
end

function scaled_directions(directions,scale)
    return [[scale .* grid for grid in direction] for direction in directions]
end

function evaluate_step!(workspace,reference_system,problem,time_step,state,
    displacement,velocity)

    UVLM.copy_system!(workspace,reference_system)
    unpack_state!(workspace,state)
    grids = [reference_system.grids[surface] .+ displacement[surface]
        for surface in eachindex(reference_system.grids)]
    UVLM.apply_trial_geometry!(workspace,grids,problem.aeroSolver)
    for surface in eachindex(grids)
        previous_grid = grids[surface] .- time_step .* velocity[surface]
        UVLM.update_surface_panels!(workspace.previous_surfaces[surface],previous_grid;
            ratios=workspace.ratios[surface],
            fcore=(chord,span) -> problem.aeroSolver.coreRadius)
    end

    nwake = copy(workspace.nwake)
    freestream = UVLM.freestream(problem.operatingPoint)
    UVLM.propagate_system!(workspace,freestream,time_step;
        additional_velocity=nothing,repeated_points=nothing,nwake,
        eta=problem.aeroSolver.wakeSheddingFraction,
        calculate_influence_matrix=true,near_field_analysis=true,derivatives=false,
        interaction_id=problem.model.interactionGroups,
        interaction=problem.aeroSolver.interaction,advance_wake=false)
    forces = pack_forces(UVLM.dimensional_loads(workspace))

    advance_wake_circulation!(workspace)
    return pack_state(workspace),forces
end

"""Finite-difference the discrete fixed-wake UVLM map around a settled state.

The aerodynamic state is bound circulation followed by active wake circulation.
Each displacement/velocity direction defines one reduced interface input. Wake
panel geometry is frozen at the reference shape while its circulation convects.
This returns an aerodynamic state-space block only; it is not yet a coupled
beam/UVLM flutter model.
"""
function linearize_uvlm_step(problem;
    displacement_directions,velocity_directions,
    finite_difference_step=1e-6,steady_tolerance=1e-3)

    isfinite(finite_difference_step) && finite_difference_step > 0 ||
        throw(ArgumentError("finite_difference_step must be positive and finite"))
    isfinite(steady_tolerance) && steady_tolerance > 0 ||
        throw(ArgumentError("steady_tolerance must be positive and finite"))

    reference_system = problem.state.system
    all(reference_system.nwake .== size.(reference_system.wakes,1)) ||
        throw(ArgumentError("Linearization requires a fully developed wake"))
    reference_system.near_field_analysis[] ||
        throw(ArgumentError("Linearization reference must contain near-field forces"))
    displacement_count = validate_directions(
        displacement_directions,reference_system.grids,"displacement_directions")
    velocity_count = validate_directions(
        velocity_directions,reference_system.grids,"velocity_directions")
    displacement_count == velocity_count || throw(DimensionMismatch(
        "Displacement and velocity direction counts must match"))
    time_step = problem.timeVector[2] - problem.timeVector[1]
    time_step > 0 || throw(ArgumentError("UVLM time step must be positive"))

    state_reference = pack_state(reference_system)
    force_reference = pack_forces(UVLM.dimensional_loads(reference_system))
    state_scale = max(problem.operatingPoint.airspeed * problem.model.reference.c,1e-12)
    force_scale = max(norm(force_reference),
        0.5 * problem.operatingPoint.density * problem.operatingPoint.airspeed^2 *
        problem.model.reference.S * 1e-12)
    zero_grids = [zeros(size(grid)) for grid in reference_system.grids]
    workspace = deepcopy(reference_system)
    next_reference,force_next = evaluate_step!(workspace,reference_system,problem,
        time_step,state_reference,zero_grids,zero_grids)
    state_drift = norm(next_reference - state_reference) /
        max(norm(state_reference),sqrt(length(state_reference)) * state_scale * 1e-12)
    force_drift = norm(force_next - force_reference) / force_scale
    state_drift <= steady_tolerance && force_drift <= steady_tolerance ||
        throw(ArgumentError("UVLM reference is not a one-step fixed point: " *
            "state drift=$state_drift, force drift=$force_drift"))

    state_count = length(state_reference)
    output_count = length(force_reference)
    input_count = 2 * displacement_count
    A = Matrix{Float64}(undef,state_count,state_count)
    B = Matrix{Float64}(undef,state_count,input_count)
    C = Matrix{Float64}(undef,output_count,state_count)
    D = Matrix{Float64}(undef,output_count,input_count)

    for column in 1:state_count
        step = finite_difference_step * max(abs(state_reference[column]),state_scale)
        state_plus,state_minus = copy(state_reference),copy(state_reference)
        state_plus[column] += step
        state_minus[column] -= step
        next_plus,force_plus = evaluate_step!(workspace,reference_system,problem,
            time_step,state_plus,zero_grids,zero_grids)
        next_minus,force_minus = evaluate_step!(workspace,reference_system,problem,
            time_step,state_minus,zero_grids,zero_grids)
        A[:,column] .= (next_plus .- next_minus) ./ (2 * step)
        C[:,column] .= (force_plus .- force_minus) ./ (2 * step)
    end

    for direction in 1:displacement_count
        displacement = scaled_directions(displacement_directions[direction],finite_difference_step)
        next_plus,force_plus = evaluate_step!(workspace,reference_system,problem,
            time_step,state_reference,displacement,zero_grids)
        next_minus,force_minus = evaluate_step!(workspace,reference_system,problem,
            time_step,state_reference,scaled_directions(displacement_directions[direction],-finite_difference_step),zero_grids)
        B[:,direction] .= (next_plus .- next_minus) ./ (2 * finite_difference_step)
        D[:,direction] .= (force_plus .- force_minus) ./ (2 * finite_difference_step)

        velocity = scaled_directions(velocity_directions[direction],finite_difference_step)
        next_plus,force_plus = evaluate_step!(workspace,reference_system,problem,
            time_step,state_reference,zero_grids,velocity)
        next_minus,force_minus = evaluate_step!(workspace,reference_system,problem,
            time_step,state_reference,zero_grids,
            scaled_directions(velocity_directions[direction],-finite_difference_step))
        column = displacement_count + direction
        B[:,column] .= (next_plus .- next_minus) ./ (2 * finite_difference_step)
        D[:,column] .= (force_plus .- force_minus) ./ (2 * finite_difference_step)
    end

    return PazyUVLMStateSpace(A,B,C,D,state_reference,force_reference,time_step,
        state_drift,force_drift)
end

end