module AeroBeamsCoupling

import AeroBeams as AB

export prepare_structural_problem!,set_nodal_wrench!,structural_checkpoint,restore_structural_checkpoint!

"""Prepare the v1 prescribed-data caches before a custom partitioned time loop."""
function prepare_structural_problem!(problem)
    AB.precompute_motion_basis_A!(problem)
    AB.precompute_BCs_history!(problem)
    AB.precompute_distributed_loads!(problem)
    AB.precompute_flap_deflection!(problem)
    AB.update_motion_basis_A!(problem,1)
    AB.update_BCs_state!(problem,1)
    return problem
end

"""Set a six-component dead wrench in basis A on the active nodal BC.

Update both prescribed functions and current state. A steady solve may rebuild
the prescribed histories; a partitioned trial instead consumes the current BC
state without reading a stale precomputed time sample.
"""
function set_nodal_wrench!(bc,wrench,time)
    length(wrench) == 6 || throw(DimensionMismatch("A nodal wrench has six components"))
    all(isfinite,wrench) || throw(ArgumentError("Nonfinite nodal wrench"))
    loadTypes = (:F1A,:F2A,:F3A,:M1A,:M2A,:M3A)
    values = Function[bc.definition.values...]
    for (i,type) in enumerate(loadTypes)
        index = findfirst(==(type),bc.definition.types)
        isnothing(index) && throw(ArgumentError("Nodal load BC is missing $type"))
        value = Float64(wrench[i])
        loadFunction = let value=value
            t -> value
        end
        values[index] = loadFunction
        bc.loadFuns.deadLoadsOnAofTime[i] = loadFunction
    end
    definition = bc.definition
    bc.definition = AB.BCDefinition(definition.name,definition.beam,definition.node,
        definition.types,values,definition.toBeTrimmed)
    AB.update_BC_state!(bc,Float64(time))
    return bc
end

"""Checkpoint mutable structural trial data without copying prescribed histories.

The current Pazy/Chang adapters do not support hinge-axis constraints. Keeping
the original model graph preserves the beam/element/BC aliases and avoids
copying all accepted solution history on every physical step.
"""
function structural_checkpoint(problem)
    isempty(problem.model.hingeAxisConstraints) || throw(ArgumentError(
        "Structural trial checkpoints require a model without hinge-axis constraints"))
    elements = problem.model.assemblyData.elements
    nodes = problem.model.assemblyData.specialNodes
    return (; systemArrays=deepcopy(problem.systemArrays),
        systemSolver=deepcopy(problem.systemSolver),timeState=deepcopy(problem.time.state),
        elementStates=[copy(e.state) for e in elements],
        elementCaches=[deepcopy(e.stateCache) for e in elements],
        nodeStates=[copy(node.state) for node in nodes],
        motionState=deepcopy(problem.model.kinematicsBasisA.state),
        BCs=[(; definition=bc.definition,loadFuns=deepcopy(bc.loadFuns),
            state=deepcopy(bc.state),stateCache=deepcopy(bc.stateCache)) for bc in problem.model.BCs])
end

function restore_structural_checkpoint!(problem,checkpoint)
    problem.systemArrays = checkpoint.systemArrays
    problem.systemSolver = checkpoint.systemSolver
    for field in fieldnames(typeof(problem.time.state))
        setfield!(problem.time.state,field,getfield(checkpoint.timeState,field))
    end
    for (element,state,cache) in zip(problem.model.assemblyData.elements,
        checkpoint.elementStates,checkpoint.elementCaches)
        element.state = state
        element.stateCache = cache
    end
    for (node,state) in zip(problem.model.assemblyData.specialNodes,checkpoint.nodeStates)
        node.state = state
    end
    for field in fieldnames(typeof(checkpoint.motionState))
        getfield(problem.model.kinematicsBasisA.state,field) .= getfield(checkpoint.motionState,field)
    end
    for (bc,saved) in zip(problem.model.BCs,checkpoint.BCs)
        # Keep the BC identities referenced by special nodes, and its active beam.
        bc.definition = saved.definition
        bc.loadFuns = saved.loadFuns
        bc.state = saved.state
        bc.stateCache = saved.stateCache
    end
    return problem
end

end
