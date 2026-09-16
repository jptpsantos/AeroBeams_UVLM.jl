"""Accepted aerodynamic state. The backend System remains a transitional storage facade."""
mutable struct UVLMState{S}
    system::S
    timeNow::Float64
    stepIndex::Int
end
"""Private trial System and scratch storage, isolated from the accepted state."""
mutable struct UVLMWorkspace{S}
    system::S
    trialReady::Bool
    stepOpen::Bool
    timeEndTimeStep::Float64
end
mutable struct UVLMResults
    completed::Bool
    terminationReason::String
    savedTimeVector::Vector{Float64}
    circulationOverTime::Vector{Vector{Float64}}
    forceOverTime::Vector{Vector{Float64}}
    momentOverTime::Vector{Vector{Float64}}
    surfaceGeometryOverTime::Vector{Any}
end
UVLMResults()=UVLMResults(false,"initialized",Float64[],Vector{Float64}[],
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

function initialize_backend(model,solver,point,maximumWakeRows,grids)
    length(grids)==length(model.grids) || throw(DimensionMismatch("Surface count changed"))
    system=System(deepcopy(model.grids);nw=maximumWakeRows)
    system.reference[]=Reference(model.reference.S,model.reference.c,model.reference.b,
        model.reference.r,point.airspeed,point.density)
    system.freestream[]=freestream(point)
    system.symmetric .= model.symmetric
    system.surface_id .= model.interactionGroups
    apply_trial_geometry!(system,grids,solver)
    copy_surfaces_to_previous!(system,length(system.surfaces))
    return system
end
function apply_trial_geometry!(system,grids,solver)
    length(grids)==length(system.grids) || throw(DimensionMismatch("Surface count changed"))
    for i in eachindex(grids)
        size(grids[i])==size(system.grids[i]) || throw(DimensionMismatch("Grid topology changed"))
        all(isfinite,grids[i]) || throw(ArgumentError("Nonfinite trial geometry"))
        system.grids[i] .= grids[i]
        _,ratio,panels=grid_to_surface_panels(system.grids[i];fcore=(c,ds)->solver.coreRadius)
        system.ratios[i] .= ratio
        system.surfaces[i] .= panels
    end
    return system
end

"""Create a fixed-grid or prescribed-motion aerodynamic analysis.
`surfaceMotion(model,time)` returns all vertex grids at that time, without
mutating model. No structural equations or automatic rotor motion are inferred.
Specify either timeVector or initialTime/Δt/finalTime. Initial wake is empty.
"""
function create_UVLMDynamicProblem(;model,aeroSolver=create_UVLMSolver(),operatingPoint,
    initialTime=0.0,Δt=nothing,finalTime=nothing,timeVector=nothing,
    maximumWakeRows=100,surfaceMotion=(m,t)->m.grids,
    trackingTimeSteps=true,trackingFrequency=1,saveGeometry=false)
    trackingFrequency isa Integer && trackingFrequency>0 || throw(ArgumentError("trackingFrequency must be positive"))
    if isnothing(timeVector)
        !isnothing(Δt) && !isnothing(finalTime) || throw(ArgumentError("Provide Δt and finalTime, or timeVector"))
        all(isfinite,(initialTime,Δt,finalTime)) && Δt>0 && finalTime>initialTime ||
            throw(ArgumentError("Invalid time interval"))
        times=collect(Float64(initialTime):Float64(Δt):Float64(finalTime))
    else
        isnothing(Δt) && isnothing(finalTime) || throw(ArgumentError("timeVector conflicts with Δt/finalTime"))
        times=Float64.(timeVector)
    end
    length(times)>=2 && all(isfinite,times) && all(>(0),diff(times)) ||
        throw(ArgumentError("At least two finite, increasing time samples are required"))
    n=length(model.grids)
    rows=maximumWakeRows isa Integer ? fill(Int(maximumWakeRows),n) : Int.(maximumWakeRows)
    length(rows)==n && all(>(0),rows) || throw(ArgumentError("Positive wake capacity required for every surface"))
    sys=initialize_backend(model,aeroSolver,operatingPoint,rows,surfaceMotion(model,first(times)))
    state=UVLMState(sys,first(times),1)
    workspace=UVLMWorkspace(deepcopy(sys),false,false,first(times))
    return UVLMDynamicProblem(model,aeroSolver,operatingPoint,times,surfaceMotion,
        state,workspace,trackingTimeSteps,trackingFrequency,saveGeometry,UVLMResults())
end

function begin_time_step!(p::UVLMDynamicProblem)
    p.workspace.stepOpen && throw(ArgumentError("A time step is already open"))
    p.state.stepIndex<length(p.timeVector) || throw(ArgumentError("No remaining time steps"))
    p.workspace.timeEndTimeStep=p.timeVector[p.state.stepIndex+1]
    p.workspace.stepOpen=true
    p.workspace.trialReady=false
    return p
end

"""Evaluate a trial from the accepted state. This never commits wake history.
Returns dimensional vertex forces and their matching application points.
"""
function evaluate_trial!(p::UVLMDynamicProblem;surfaces=nothing)
    w=p.workspace
    w.stepOpen || throw(ArgumentError("Call begin_time_step! first"))
    w.trialReady=false
    w.system=deepcopy(p.state.system)
    sys=w.system
    copy_surfaces_to_previous!(sys,length(sys.surfaces))
    grids=isnothing(surfaces) ? p.surfaceMotion(p.model,w.timeEndTimeStep) : surfaces
    apply_trial_geometry!(sys,grids,p.aeroSolver)
    propagate_system!(sys,freestream(p.operatingPoint),w.timeEndTimeStep-p.state.timeNow;
        additional_velocity=nothing,repeated_points=repeated_trailing_edge_points(sys.surfaces),
        nwake=copy(p.state.system.nwake),eta=p.aeroSolver.wakeSheddingFraction,
        calculate_influence_matrix=true,near_field_analysis=true,derivatives=false,
        interaction_id=p.model.interactionGroups,interaction=p.aeroSolver.interaction,advance_wake=false)
    all(isfinite,sys.Γ) && all(isfinite,sys.dΓdt) || error("Nonfinite aerodynamic trial")
    loads=dimensional_loads(sys)
    all(f -> all(v->all(isfinite,v),f),loads.forces) || error("Nonfinite aerodynamic forces")
    w.trialReady=true
    return loads
end

"""Dimensional vertex loads in backend coordinates; positions use the same force geometry."""
function dimensional_loads(system::System)
    return (;forces=imperial_nodal_forces(system),positions=imperial_nodal_positions(system))
end
dimensional_loads(p::UVLMDynamicProblem)=dimensional_loads(p.state.system)

function save_time_step!(p::UVLMDynamicProblem)
    loads=dimensional_loads(p)
    force=zeros(3); moment=zeros(3)
    for i in eachindex(loads.forces), j in eachindex(loads.forces[i])
        f=loads.forces[i][j]; r=loads.positions[i][j]-p.model.reference.r
        force .+= f; moment .+= cross(r,f)
    end
    result=p.results
    push!(result.savedTimeVector,p.state.timeNow)
    push!(result.circulationOverTime,copy(p.state.system.Γ))
    push!(result.forceOverTime,force); push!(result.momentOverTime,moment)
    p.saveGeometry && push!(result.surfaceGeometryOverTime,deepcopy(p.state.system.surfaces))
    return p
end

function commit_time_step!(p::UVLMDynamicProblem)
    w=p.workspace
    w.stepOpen && w.trialReady || throw(ArgumentError("A successful uncommitted trial is required"))
    # Commit on a private copy: a wake-convection failure leaves the trial retryable.
    sys=deepcopy(w.system)
    rows=copy(p.state.system.nwake)
    advance_wake!(sys,freestream(p.operatingPoint),w.timeEndTimeStep-p.state.timeNow;
        nwake=rows,interaction_id=p.model.interactionGroups,interaction=p.aeroSolver.interaction)
    commit_wake_rows!(rows,size.(sys.wakes,1)); sys.nwake .= rows
    p.state.system=sys
    p.state.timeNow=w.timeEndTimeStep; p.state.stepIndex+=1
    w.stepOpen=false; w.trialReady=false
    if p.trackingTimeSteps && (mod(p.state.stepIndex-1,p.trackingFrequency)==0 || p.state.stepIndex==length(p.timeVector))
        save_time_step!(p)
    end
    p.results.completed=p.state.stepIndex==length(p.timeVector)
    p.results.terminationReason=p.results.completed ? "completed" : "running"
    return p
end
function rollback_time_step!(p::UVLMDynamicProblem)
    p.workspace.system=deepcopy(p.state.system)
    p.workspace.stepOpen=false; p.workspace.trialReady=false
    p.workspace.timeEndTimeStep=p.state.timeNow
    return p
end
function solve!(p::UVLMDynamicProblem)
    p.workspace.stepOpen && throw(ArgumentError("Finish or roll back the open trial before solve!"))
    try
        while p.state.stepIndex<length(p.timeVector)
            begin_time_step!(p); evaluate_trial!(p); commit_time_step!(p)
        end
    catch e
        rollback_time_step!(p)
        p.results.terminationReason=sprint(showerror,e)
        rethrow()
    end
    return p
end

mutable struct UVLMSteadyProblem{M}
    model::M
    aeroSolver::UVLMSolver
    operatingPoint::OperatingPoint
    system::Union{Nothing,System}
end
function create_UVLMSteadyProblem(;model,aeroSolver=create_UVLMSolver(),operatingPoint)
    aeroSolver.interaction || throw(ArgumentError("Steady facade currently requires full interaction"))
    return UVLMSteadyProblem(model,aeroSolver,operatingPoint,nothing)
end
function solve!(p::UVLMSteadyProblem)
    m=p.model; op=p.operatingPoint
    ref=Reference(m.reference.S,m.reference.c,m.reference.b,m.reference.r,op.airspeed,op.density)
    p.system=steady_analysis(deepcopy(m.grids),ref,freestream(op);symmetric=m.symmetric,
        surface_id=m.interactionGroups,fcore=(c,ds)->p.aeroSolver.coreRadius,derivatives=false)
    return p
end
