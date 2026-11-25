"""Accepted aerodynamic state. The backend System remains a transitional storage facade."""
mutable struct UVLMState{S}
    system::S
    timeNow::Float64
    stepIndex::Int
end
"""Private trial System and scratch storage, isolated from the accepted state."""
mutable struct UVLMWorkspace{S}
    system::S
    commitSystem::S
    trialReady::Bool
    stepOpen::Bool
    timeEndTimeStep::Float64
end
UVLMWorkspace(system,trialReady,stepOpen,timeEndTimeStep) =
    UVLMWorkspace(system,deepcopy(system),trialReady,stepOpen,timeEndTimeStep)
mutable struct UVLMResults
    completed::Bool
    terminationReason::String
    savedTimeVector::Vector{Float64}
    circulationOverTime::Vector{Vector{Float64}}
    forceOverTime::Vector{Vector{Float64}}
    momentOverTime::Vector{Vector{Float64}}
    surfaceGeometryOverTime::Vector{Any}
end
UVLMResults() = UVLMResults(false,"initialized",Float64[],Vector{Float64}[],
    Vector{Float64}[],Vector{Float64}[],Any[])

mutable struct UVLMDynamicProblem{M,F,S}
    model::M
    aeroSolver::UVLMSolver
    operatingPoint::OperatingPoint
    timeVector::Vector{Float64}
    surfaceMotion::F
    state::UVLMState{S}
    workspace::UVLMWorkspace{S}
    trackingTimeSteps::Bool
    trackingFrequency::Int
    saveGeometry::Bool
    results::UVLMResults
end

function initialize_backend(model,aeroSolver,operatingPoint,maximumWakeRows,grids)
    # Allocate one aerodynamic System and initialize it at the accepted t_0 geometry.
    length(grids) == length(model.grids) || throw(DimensionMismatch("Surface count changed"))
    system = System(deepcopy(model.grids); nw=maximumWakeRows)
    system.reference[] = Reference(model.reference.S,model.reference.c,model.reference.b,
        model.reference.r,operatingPoint.airspeed,operatingPoint.density)
    system.freestream[] = freestream(operatingPoint)
    system.symmetric .= model.symmetric
    system.surface_id .= model.interactionGroups
    apply_trial_geometry!(system,grids,aeroSolver)
    copy_surfaces_to_previous!(system,length(system.surfaces))
    return system
end
function apply_trial_geometry!(system,grids,aeroSolver)
    length(grids) == length(system.grids) || throw(DimensionMismatch("Surface count changed"))
    for i in eachindex(grids)
        size(grids[i]) == size(system.grids[i]) || throw(DimensionMismatch("Grid topology changed"))
        all(isfinite,grids[i]) || throw(ArgumentError("Nonfinite trial geometry"))
        system.grids[i] .= grids[i]
        update_surface_panels!(system.surfaces[i],system.grids[i];
            ratios=system.ratios[i],fcore=(c,ds) -> aeroSolver.coreRadius)
    end
    return system
end

"""Create a fixed-grid or prescribed-motion aerodynamic analysis.
`surfaceMotion(model,time)` returns all vertex grids at that time, without
mutating model. No structural equations or automatic rotor motion are inferred.
Specify either timeVector or initialTime/Δt/finalTime. Initial wake is empty.
"""
function create_UVLMDynamicProblem(; model,aeroSolver=create_UVLMSolver(),operatingPoint,
    initialTime=0.0,Δt=nothing,finalTime=nothing,timeVector=nothing,
    maximumWakeRows=100,surfaceMotion=(model,t) -> model.grids,
    trackingTimeSteps=true,trackingFrequency=1,saveGeometry=false)
    trackingFrequency isa Integer && trackingFrequency > 0 || throw(ArgumentError("trackingFrequency must be positive"))
    if isnothing(timeVector)
        !isnothing(Δt) && !isnothing(finalTime) || throw(ArgumentError("Provide Δt and finalTime, or timeVector"))
        all(isfinite,(initialTime,Δt,finalTime)) && Δt > 0 && finalTime > initialTime ||
            throw(ArgumentError("Invalid time interval"))
        resolvedTimeVector = collect(Float64(initialTime):Float64(Δt):Float64(finalTime))
    else
        isnothing(Δt) && isnothing(finalTime) || throw(ArgumentError("timeVector conflicts with Δt/finalTime"))
        resolvedTimeVector = Float64.(timeVector)
    end
    length(resolvedTimeVector) >= 2 && all(isfinite,resolvedTimeVector) && all(>(0),diff(resolvedTimeVector)) ||
        throw(ArgumentError("At least two finite, increasing time samples are required"))
    nSurfaces = length(model.grids)
    wakeRows = maximumWakeRows isa Integer ? fill(Int(maximumWakeRows),nSurfaces) : Int.(maximumWakeRows)
    length(wakeRows) == nSurfaces && all(>(0),wakeRows) || throw(ArgumentError("Positive wake capacity required for every surface"))
    system = initialize_backend(model,aeroSolver,operatingPoint,wakeRows,surfaceMotion(model,first(resolvedTimeVector)))
    state = UVLMState(system,first(resolvedTimeVector),1)
    workspace = UVLMWorkspace(deepcopy(system),false,false,first(resolvedTimeVector))
    return UVLMDynamicProblem(model,aeroSolver,operatingPoint,resolvedTimeVector,surfaceMotion,
        state,workspace,trackingTimeSteps,trackingFrequency,saveGeometry,UVLMResults())
end

function begin_time_step!(problem::UVLMDynamicProblem)
    problem.workspace.stepOpen && throw(ArgumentError("A time step is already open"))
    problem.state.stepIndex < length(problem.timeVector) || throw(ArgumentError("No remaining time steps"))
    # Open one n -> n+1 transaction without changing the accepted state at n.
    problem.workspace.timeEndTimeStep = problem.timeVector[problem.state.stepIndex+1]
    problem.workspace.stepOpen = true
    problem.workspace.trialReady = false
    return problem
end

"""Evaluate a trial from the accepted state. This never commits wake history.
Returns dimensional vertex forces and their matching application points.
"""
function evaluate_trial!(problem::UVLMDynamicProblem; surfaces=nothing)
    workspace = problem.workspace
    workspace.stepOpen || throw(ArgumentError("Call begin_time_step! first"))
    workspace.trialReady = false
    # Every trial starts from the same accepted circulation, geometry, and wake.
    copy_system!(workspace.system,problem.state.system)
    system = workspace.system
    # UVLM uses the accepted and trial surfaces to obtain the surface velocity.
    copy_surfaces_to_previous!(system,length(system.surfaces))
    grids = isnothing(surfaces) ? problem.surfaceMotion(problem.model,workspace.timeEndTimeStep) : surfaces
    apply_trial_geometry!(system,grids,problem.aeroSolver)
    # Solve circulation and loads, but leave wake convection for commit_time_step!.
    propagate_system!(system,freestream(problem.operatingPoint),workspace.timeEndTimeStep - problem.state.timeNow;
        additional_velocity=nothing,repeated_points=nothing,
        nwake=system.nwake,eta=problem.aeroSolver.wakeSheddingFraction,
        calculate_influence_matrix=true,near_field_analysis=true,derivatives=false,
        interaction_id=problem.model.interactionGroups,interaction=problem.aeroSolver.interaction,advance_wake=false)
    all(isfinite,system.Γ) && all(isfinite,system.dΓdt) || error("Nonfinite aerodynamic trial")
    loads = dimensional_loads(system)
    all(f -> all(v -> all(isfinite,v),f),loads.forces) || error("Nonfinite aerodynamic forces")
    workspace.trialReady = true
    return loads
end

"""Dimensional vertex loads in backend coordinates; positions use the same force geometry."""
function dimensional_loads(system::System)
    return (; forces=imperial_nodal_forces(system),positions=imperial_nodal_positions(system))
end
dimensional_loads(problem::UVLMDynamicProblem) = dimensional_loads(problem.state.system)

function save_time_step!(problem::UVLMDynamicProblem)
    loads = dimensional_loads(problem)
    force = zeros(3)
    moment = zeros(3)
    for i in eachindex(loads.forces),j in eachindex(loads.forces[i])
        f = loads.forces[i][j]
        r = loads.positions[i][j] - problem.model.reference.r
        force .+= f
        moment .+= cross(r,f)
    end
    results = problem.results
    push!(results.savedTimeVector,problem.state.timeNow)
    push!(results.circulationOverTime,copy(problem.state.system.Γ))
    push!(results.forceOverTime,force)
    push!(results.momentOverTime,moment)
    problem.saveGeometry && push!(results.surfaceGeometryOverTime,deepcopy(problem.state.system.surfaces))
    return problem
end

function commit_time_step!(problem::UVLMDynamicProblem)
    workspace = problem.workspace
    workspace.stepOpen && workspace.trialReady || throw(ArgumentError("A successful uncommitted trial is required"))
    # The preallocated commit buffer keeps failures isolated from accepted/trial state.
    system = copy_system!(workspace.commitSystem,workspace.system)
    # A successful trial becomes physical only after its wake is advanced once.
    wakeRows = system.nwake
    advance_wake!(system,freestream(problem.operatingPoint),workspace.timeEndTimeStep - problem.state.timeNow;
        nwake=wakeRows,interaction_id=problem.model.interactionGroups,interaction=problem.aeroSolver.interaction)
    for i in eachindex(wakeRows)
        wakeRows[i] = min(wakeRows[i] + 1,size(system.wakes[i],1))
    end
    workspace.commitSystem = problem.state.system
    problem.state.system = system
    problem.state.timeNow = workspace.timeEndTimeStep
    problem.state.stepIndex += 1
    workspace.stepOpen = false
    workspace.trialReady = false
    if problem.trackingTimeSteps &&
       (mod(problem.state.stepIndex - 1,problem.trackingFrequency) == 0 || problem.state.stepIndex == length(problem.timeVector))
        save_time_step!(problem)
    end
    problem.results.completed = problem.state.stepIndex == length(problem.timeVector)
    problem.results.terminationReason = problem.results.completed ? "completed" : "running"
    return problem
end
function rollback_time_step!(problem::UVLMDynamicProblem)
    copy_system!(problem.workspace.system,problem.state.system)
    problem.workspace.stepOpen = false
    problem.workspace.trialReady = false
    problem.workspace.timeEndTimeStep = problem.state.timeNow
    return problem
end
function solve!(problem::UVLMDynamicProblem)
    problem.workspace.stepOpen && throw(ArgumentError("Finish or roll back the open trial before solve!"))
    try
        while problem.state.stepIndex < length(problem.timeVector)
            begin_time_step!(problem)
            evaluate_trial!(problem)
            commit_time_step!(problem)
        end
    catch e
        rollback_time_step!(problem)
        problem.results.terminationReason = sprint(showerror,e)
        rethrow()
    end
    return problem
end

mutable struct UVLMSteadyProblem{M}
    model::M
    aeroSolver::UVLMSolver
    operatingPoint::OperatingPoint
    system::Union{Nothing,System}
end
function create_UVLMSteadyProblem(; model,aeroSolver=create_UVLMSolver(),operatingPoint)
    aeroSolver.interaction || throw(ArgumentError("Steady facade currently requires full interaction"))
    return UVLMSteadyProblem(model,aeroSolver,operatingPoint,nothing)
end
function solve!(problem::UVLMSteadyProblem)
    model = problem.model
    operatingPoint = problem.operatingPoint
    reference = Reference(model.reference.S,model.reference.c,model.reference.b,model.reference.r,operatingPoint.airspeed,operatingPoint.density)
    problem.system = steady_analysis(deepcopy(model.grids),reference,freestream(operatingPoint); symmetric=model.symmetric,
        surface_id=model.interactionGroups,fcore=(c,ds) -> problem.aeroSolver.coreRadius,derivatives=false)
    return problem
end
