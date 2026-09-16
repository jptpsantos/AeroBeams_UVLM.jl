"""Reference aerodynamic vertex grids and coefficient references; treat as read-only.
Grids have size (3, nChordwisePanels+1, nSpanwisePanels+1), in metres.
"""
struct UVLMModel{G,R}
    name::String
    grids::G
    reference::R
    symmetric::Vector{Bool}
    interactionGroups::Vector{Int}
end
function create_UVLMModel(; surfaces,reference,name="",symmetric=false,interactionGroups=nothing)
    isempty(surfaces) && throw(ArgumentError("At least one vertex grid is required"))
    for g in surfaces
        ndims(g)==3 && size(g,1)==3 && size(g,2)>=2 && size(g,3)>=2 ||
            throw(ArgumentError("surfaces must contain (3,nChord+1,nSpan+1) vertex grids"))
        all(isfinite,g) || throw(ArgumentError("Geometry must be finite"))
    end
    reference isa Reference || throw(ArgumentError("reference must be a Reference"))
    all(x -> isfinite(x) && x>0,(reference.S,reference.c,reference.b,reference.V,reference.rho)) &&
        all(isfinite,reference.r) || throw(ArgumentError("Invalid reference quantities"))
    n=length(surfaces)
    sym=symmetric isa Bool ? fill(symmetric,n) : Bool.(symmetric)
    ids=isnothing(interactionGroups) ? collect(1:n) : Int.(interactionGroups)
    length(sym)==length(ids)==n || throw(DimensionMismatch("One symmetry/group entry per surface is required"))
    return UVLMModel(String(name),[Float64.(g) for g in surfaces],deepcopy(reference),sym,ids)
end

"""Aerodynamic algorithm controls. AIC assembly remains enabled for moving geometry."""
struct UVLMSolver
    coreRadius::Float64
    interaction::Bool
    wakeSheddingFraction::Float64
end
function create_UVLMSolver(;coreRadius=1e-3,interaction=true,wakeSheddingFraction=.1)
    isfinite(coreRadius) && coreRadius>0 || throw(ArgumentError("coreRadius must be positive"))
    isfinite(wakeSheddingFraction) && 0<=wakeSheddingFraction<=1 ||
        throw(ArgumentError("wakeSheddingFraction must lie in [0,1]"))
    return UVLMSolver(coreRadius,interaction,wakeSheddingFraction)
end

"""AeroBeams-style constructor for the existing linear partitioned solver options."""
function create_PartitionedCoupling(;maximumIterations=10,stateTolerance=1e-5,
    loadTolerance=1e-2,equilibriumTolerance=1e-10,coupledEquilibriumTolerance=1e-4,
    relaxation=1.0)
    options=PartitionedCouplingOptions(;maximum_iterations=maximumIterations,
        state_tolerance=stateTolerance,load_tolerance=loadTolerance,
        equilibrium_tolerance=equilibriumTolerance,
        coupled_equilibrium_tolerance=coupledEquilibriumTolerance,relaxation)
    validate_partitioned_coupling_options(options)
    return options
end

