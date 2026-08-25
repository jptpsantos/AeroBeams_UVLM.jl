# Coupling AeroBeams with a wing--propeller UVLM

## Purpose and status

This document is the implementation blueprint for coupling the geometrically
exact nonlinear structural formulation in AeroBeams with the wing--propeller
unsteady vortex-lattice method (UVLM) developed in
[`Wing_Propeller_UVLM`](https://github.com/jptpsantos/Wing_Propeller_UVLM).
It covers the software architecture, coordinate transformations, structural and
aerodynamic state ownership, conservative transfer of motion and loads,
propeller equations, nonlinear time marching, verification, and representative
Julia implementations.

The examples below follow AeroBeams conventions:

- Julia source files and multiple dispatch;
- `CamelCase` names for composite types;
- `snake_case` names for functions, with `!` on mutating functions;
- keyword constructors and explicit input validation;
- `@unpack` and `@pack!` where they improve readability;
- no solver state stored in module-level globals;
- comments that describe the physical operation, not only the syntax.

Code marked **existing call** uses an API already present in one of the two
repositories. Code marked **proposed interface** is the concrete interface to
add during the coupling work. The UVLM repository must first be converted into
a Julia package before `import WingPropellerUVLM` will work as shown.

The reference implementation inspected for this design is commit
[`12f58e73`](https://github.com/jptpsantos/Wing_Propeller_UVLM/commit/12f58e73fedc09c32c6fc4a2305fdc62f5a48872)
of `Wing_Propeller_UVLM` and AeroBeams version 0.8.1.

The UVLM source has now been imported into this repository as the internal
package `lib/WingPropellerUVLM`. Its `src/backend` directory contains the
modified VortexLattice implementation, `src/wing_propeller` contains the
reusable research-model helpers, and `src/UVLMState.jl` provides the
transactional snapshot/restore interface described below. The remaining
AeroBeams coupling hooks in this document are still proposed work.

## 1. Design decision

The UVLM must remain a global aerodynamic subsystem. It should not be inserted
into AeroBeams' existing element-local `AeroSurface` implementation.

An AeroBeams aerodynamic element owns a small local vector of aerodynamic
states, `χ`, and evaluates strip-theory loads independently for each beam
element. The wing--propeller UVLM instead owns a dense, global state shared by
the wing, every propeller blade, and every wake:

- current and previous surface panels;
- bound circulation `Γ` and its time derivatives;
- wing and propeller free wakes;
- active wake lengths;
- aerodynamic interaction identifiers;
- panel and segment forces.

Those states cannot be divided into independent beam-element states without
destroying the mutual aerodynamic interactions. The target implementation is
therefore a strongly partitioned Dirichlet--Neumann coupling:

```text
accepted AeroBeams state at t_n
             |
             v
  structural interface motion guess
             |
             v
  motion transfer and UVLM grid update
             |
             v
 trial UVLM solve restored from state at t_n
             |
             v
 conservative wing loads + propeller hub wrenches
             |
             v
 nonlinear AeroBeams solve from t_n to t_(n+1)
             |
             +------ relaxed outer iteration ------+

after convergence: commit structural state, Γ, and wakes exactly once
```

The UVLM load is fixed during each AeroBeams Newton solve. The UVLM is
reevaluated only in the outer coupling loop. Advancing the free wake from inside
`element_arrays!` or a Newton residual evaluation would advance physical time
once per Newton iteration and produce a path-dependent, incorrect solution.

## 2. Separation into three software layers

The preferred final organization has three packages or modules.

```text
AeroBeams.jl
|-- geometrically exact structure
|-- generic externally supplied nodal resultants
`-- public one-step dynamic API

WingPropellerUVLM.jl
|-- UVLM system creation
|-- wing and blade grids
|-- circulation and wake state
|-- trial propagation and rollback
`-- panel/segment force access

AeroBeamsUVLMCoupling.jl
|-- frame transformations
|-- interface mesh and interpolation stencils
|-- motion transfer
|-- load transfer
|-- spinning-rotor reactions
`-- coupled time-step controller
```

During initial research development, the coupling layer can live in this fork:

```text
src/
|-- AeroBeams.jl
|-- ExternalResultants.jl
|-- Problem.jl
|-- SystemSolver.jl
`-- Coupling/
    |-- Coupling.jl
    |-- CouplingFrames.jl
    |-- InterfaceMesh.jl
    |-- MotionTransfer.jl
    |-- LoadTransfer.jl
    |-- UVLMAdapter.jl
    |-- SpinningRotor.jl
    `-- CoupledDynamicProblem.jl

test/
`-- coupling/
    |-- frameTests.jl
    |-- motionTransferTests.jl
    |-- loadTransferTests.jl
    |-- uvlmRollbackTests.jl
    |-- structuralStepTests.jl
    `-- coupledRegressionTests.jl

examples/
`-- changWingPropellerUVLM.jl
```

Once the interfaces stabilize, `Coupling/` should become a small standalone
package. Raw UVLM source files should not be copied into `src/AeroBeams.jl` with
`include`. A package boundary provides dependency locking, namespace isolation,
precompilation, reproducible tests, and a stable API.

The current driver functions migrate as follows:

| Current implementation | Target owner and API |
|:--|:--|
| `initialize_bohnisch_uvlm_system` | `WingPropellerUVLM.create_UVLMAssembly` |
| `snapshot_aero_step` | `WingPropellerUVLM.snapshot` |
| `restore_aero_step!` | `WingPropellerUVLM.restore!` |
| `update_aero_geometry_for_state!` | `MotionTransfer.update_geometry!` |
| `assemble_structural_aero_load!` | `LoadTransfer.transfer_loads!` |
| `aero_load_for_state!` | `restore!` + `update_geometry!` + `advance_trial!` + `transfer_loads!` |
| `generalized_alpha_corrector` | `AeroBeams.solve_dynamic_step!` for the first nonlinear implementation |
| top-level time loop | `solve_coupled_step!` and a thin simulation driver |

## 3. State ownership and transaction rules

Each subsystem must have one owner for every mutable quantity.

| Quantity | Owner | Trial behavior | Commit behavior |
|:--|:--|:--|:--|
| AeroBeams mixed state `x` | AeroBeams | restored to the beginning of the step before each structural solve | retained after convergence |
| Element states and rates | AeroBeams | recomputed by the nonlinear solve | saved to history |
| Interface position/orientation guess | Coupler | relaxed by the outer iteration | discarded after extracting the accepted state |
| UVLM surfaces | UVLM | replaced by the trial deformed grids | retained for the accepted geometry |
| `Γ`, `dΓ`, `dΓdt` | UVLM | restored before every aerodynamic trial | advanced once |
| Wake panels | UVLM | restored before every aerodynamic trial | shed and convected once |
| External nodal resultants | Coupler/AeroBeams bridge | fixed during an inner structural Newton solve | overwritten at the next outer iteration |
| Rotor phase `φ` | Rotor controller | evaluated at `t_(n+1)` | advanced once |

The required invariant is:

> Every evaluation in an outer coupling iteration starts from the same accepted
> structural and aerodynamic states at `t_n`.

This is the same sound transaction principle already used by
`snapshot_aero_step`, `restore_aero_step!`, and `aero_load_for_state!` in the
current partitioned UVLM driver. The production implementation must extend the
snapshot to all mutable aerodynamic history, including surface history and
active wake counters.

## 4. Coordinate systems

### 4.1 Select one right-handed mapping

The existing UVLM uses a right-handed frame with approximately

- `U1 = X`: chord/freestream direction;
- `U2 = Y`: wing span direction;
- `U3 = Z`: upward.

The current Chang structural driver uses span, chord, and downward displacement
as its six-degree-of-freedom nodal convention. Preserve that convention in the
first AeroBeams implementation:

- `A1 = +U2`: span;
- `A2 = +U1`: chord/freestream;
- `A3 = -U3`: downward.

This is important. Merely swapping chord and span while retaining `A3 = +U3`
would produce a matrix with determinant `-1`, which is a reflection rather than
a valid rotation. The proper rotation mapping a vector resolved in basis `A`
to basis `U` is

```math
R_{UA} =
\begin{bmatrix}
0 & 1 & 0 \\
1 & 0 & 0 \\
0 & 0 & -1
\end{bmatrix},
\qquad
\boldsymbol{v}_U=R_{UA}\boldsymbol{v}_A.
```

The same proper rotation transforms polar vectors (positions and forces) and
axial vectors (rotations and moments).

### 4.2 Frame type

Create one immutable transformation object; do not scatter component swaps and
minus signs throughout geometry and load routines.

```julia
# src/Coupling/CouplingFrames.jl

struct CouplingFrames
    R_UA::Matrix{Float64}
    R_AU::Matrix{Float64}
end

function create_CouplingFrames(; R_UA::Matrix{Float64} =
    [0.0 1.0  0.0;
     1.0 0.0  0.0;
     0.0 0.0 -1.0])

    @assert size(R_UA) == (3,3)
    @assert isapprox(R_UA' * R_UA, Matrix(I,3,3); atol=1e-12)
    @assert isapprox(det(R_UA), 1.0; atol=1e-12)

    return CouplingFrames(R_UA, Matrix(R_UA'))
end

to_uvlm(frames::CouplingFrames, v_A::AbstractVector) = frames.R_UA * v_A
to_aerobeams(frames::CouplingFrames, v_U::AbstractVector) = frames.R_AU * v_U
```

The first tests should be explicit:

```julia
@test to_uvlm(frames,[1.0;0.0;0.0]) ≈ [0.0;1.0;0.0] # span
@test to_uvlm(frames,[0.0;1.0;0.0]) ≈ [1.0;0.0;0.0] # chord
@test to_uvlm(frames,[0.0;0.0;1.0]) ≈ [0.0;0.0;-1.0] # down
```

### 4.3 Rotation notation

Use the notation `R_XY` for a matrix that maps components from basis `Y` into
basis `X`:

```math
\boldsymbol{v}_X=R_{XY}\boldsymbol{v}_Y.
```

For a deformed AeroBeams section, `R_AB` maps section-basis components to model
basis `A`. At an element midpoint, the existing `element.RR0` is used in that
role by the residual implementation. At a node, construct it from the nodal
Wiener--Milenkovic parameters and the undeformed nodal frame:

```julia
function nodal_rotation_A(element::Element, side::Int)
    p_A = side == 1 ? element.nodalStates.p_n1 : element.nodalStates.p_n2
    R0_AB = side == 1 ? element.R0_n1 : element.R0_n2
    R_deformation = rotation_tensor_WM(p_A)[1]
    return R_deformation * R0_AB
end
```

This convention must be verified with rigid-rotation tests before aerodynamic
loads are enabled.

## 5. AeroBeams structural model

### 5.1 Wing sectional properties

AeroBeams expects one sectional stiffness matrix and one sectional inertia
matrix per element, or one constant matrix for the whole beam. Preserve the full
`6 x 6` matrices wherever possible. Do not reduce anisotropic coupling terms to
only `EA`, `GJ`, `EIy`, and `EIz` unless the reference model is truly uncoupled.

For an uncoupled preliminary model:

```julia
using AeroBeams

S_wing = [isotropic_stiffness_matrix(
    EA  = EA[e],
    GAy = GAy[e],
    GAz = GAz[e],
    GJ  = GJ[e],
    EIy = EIy[e],
    EIz = EIz[e],
) for e in eachindex(EA)]

I_wing = [inertia_matrix(
    ρA   = ρA[e],
    ρIy  = ρIy[e],
    ρIz  = ρIz[e],
    ρIyz = ρIyz[e],
    e2   = e2[e],
    e3   = e3[e],
) for e in eachindex(ρA)]

wing = create_Beam(
    name = "wing",
    length = spanLength,
    nElements = length(S_wing),
    normalizedNodalPositions = spanNodes ./ spanLength,
    S = S_wing,
    I = I_wing,
)
```

For the Chang data, confirm the placement and sign of `EIzy` through a unit
curvature test before populating off-diagonal entries. The structural matrix
ordering used by the old linear beam is not automatically identical to the
sectional strain/curvature ordering in AeroBeams.

Nodal lumped masses from the linear model should not be silently divided by an
element length. Either reconstruct the physical sectional mass distribution or
attach explicitly validated point inertias. The acceptance checks are total
mass, center of mass, inertia tensor, static flexibility, and modes.

### 5.2 Pylons as geometrically exact beams

The preferred nonlinear model represents each pylon as an AeroBeams beam
connected to its wing node. With the frame selected above, a pylon extending in
the negative UVLM `X` direction extends in negative AeroBeams `A2`. For the
current `E321` implementation, `p0=[-π/2;0;0]` rotates the pylon local `x1` axis
from `+A1` to `-A2`.

```julia
function create_pylon(; name, wing, wingNode, length, nElements,
    S, I, hubMass, hubInertia)

    pylon = create_Beam(
        name = name,
        length = length,
        nElements = nElements,
        S = length(S) == 1 ? S : copy(S),
        I = length(I) == 1 ? I : copy(I),
        rotationParametrization = "E321",
        p0 = [-π/2; 0.0; 0.0],
        connectedBeams = [wing],
        connectedNodesThis = [1],
        connectedNodesOther = [wingNode],
    )

    # PointInertia is attached to an element midpoint. Offset it from the last
    # element midpoint to place it at the physical pylon tip.
    lastElementLength = pylon.elements[end].Δℓ
    hub = PointInertia(
        elementID = nElements,
        η = [lastElementLength/2; 0.0; 0.0],
        mass = hubMass,
        inertiaMatrix = hubInertia,
    )
    add_point_inertias_to_beam!(pylon,inertias=[hub])

    return pylon
end
```

Create the model only after all beams and point inertias are finalized:

```julia
pylons = [create_pylon(
    name = "pylon$ip",
    wing = wing,
    wingNode = propellerWingNodes[ip],
    length = pylonLength[ip],
    nElements = nElemPylon[ip],
    S = S_pylon[ip],
    I = I_pylon[ip],
    hubMass = hubMass[ip],
    hubInertia = hubInertia[ip],
) for ip in eachindex(propellerWingNodes)]

clamp = create_BC(
    name = "wingRootClamp",
    beam = wing,
    node = 1,
    types = ["u1A","u2A","u3A","p1A","p2A","p3A"],
    values = zeros(6),
)

model = create_Model(
    name = "ChangWingPropeller",
    beams = [wing; pylons],
    BCs = [clamp],
    v_A = [0.0; airspeed; 0.0],
)
```

Do not retain the old two-degree-of-freedom pylon stiffness and inertia while
also using physical pylon beams. That would double-count pylon flexibility,
mass, and gyroscopic coupling.

## 6. Generic external-resultant interface in AeroBeams

The current distributed-load callbacks are evaluated as prescribed functions of
space and time. UVLM loads instead change after every outer coupling iteration.
AeroBeams needs a generic container for externally supplied nodal forces and
moments.

### 6.1 Data type

```julia
# src/ExternalResultants.jl -- proposed interface

@with_kw mutable struct ExternalNodalResultants
    forces_A::Matrix{Float64}
    moments_A::Matrix{Float64}
    forceEquations::Vector{Vector{Int}}
    momentEquations::Vector{Vector{Int}}
    active::Bool = true
end

# Place this model-dependent constructor in Model.jl, after Model is defined.
function create_ExternalNodalResultants(model::Model)
    nNodes = model.nNodesTotal
    forceEquations = [Int[] for _ in 1:nNodes]
    momentEquations = [Int[] for _ in 1:nNodes]

    for element in model.elements
        sides = (
            (element.nodesGlobalID[1],element.eqs_Fu1,element.eqs_Fp1),
            (element.nodesGlobalID[2],element.eqs_Fu2,element.eqs_Fp2),
        )
        for (node,eqs_Fu,eqs_Fp) in sides
            if isempty(forceEquations[node])
                forceEquations[node] = copy(eqs_Fu)
                momentEquations[node] = copy(eqs_Fp)
            end
        end
    end

    @assert all(!isempty,forceEquations)
    @assert all(!isempty,momentEquations)

    return ExternalNodalResultants(
        forces_A = zeros(3,nNodes),
        moments_A = zeros(3,nNodes),
        forceEquations = forceEquations,
        momentEquations = momentEquations,
    )
end

function clear!(loads::ExternalNodalResultants)
    fill!(loads.forces_A,0.0)
    fill!(loads.moments_A,0.0)
    return loads
end

function add_wrench!(loads::ExternalNodalResultants,node::Int,
    force_A::AbstractVector,moment_A::AbstractVector)

    @assert 1 <= node <= size(loads.forces_A,2)
    @assert length(force_A) == length(moment_A) == 3
    loads.forces_A[:,node] .+= force_A
    loads.moments_A[:,node] .+= moment_A
    return loads
end

function similar_external_loads(loads::ExternalNodalResultants)
    return ExternalNodalResultants(
        forces_A = zeros(size(loads.forces_A)),
        moments_A = zeros(size(loads.moments_A)),
        forceEquations = loads.forceEquations,
        momentEquations = loads.momentEquations,
        active = loads.active,
    )
end
```

The simple equation map above is appropriate for rigidly connected wing and
pylon nodes. Before supporting hinge moments, extend the map to distinguish the
separate rotational equilibrium equations used by AeroBeams at a hinge.

Add an optional field to `Model`:

```julia
externalNodalResultants::Union{Nothing,ExternalNodalResultants} = nothing
```

and a setter:

```julia
function set_external_nodal_resultants!(model::Model,
    loads::ExternalNodalResultants)

    @assert size(loads.forces_A) == (3,model.nNodesTotal)
    @assert size(loads.moments_A) == (3,model.nNodesTotal)
    model.externalNodalResultants = loads
    return model
end
```

### 6.2 Residual assembly hook

In `assemble_system_arrays!`, apply the supplied resultants after the regular
element and special-node residuals have been assembled:

```julia
function apply_external_nodal_resultants!(problem::Problem)
    loads = problem.model.externalNodalResultants
    if isnothing(loads) || !loads.active
        return nothing
    end

    forceScaling = problem.model.forceScaling
    for node in axes(loads.forces_A,2)
        problem.residual[loads.forceEquations[node]] .-=
            loads.forces_A[:,node] ./ forceScaling
        problem.residual[loads.momentEquations[node]] .-=
            loads.moments_A[:,node] ./ forceScaling
    end

    return nothing
end
```

Then add one call at the end of the existing assembly function:

```julia
for specialNode in specialNodes
    special_node_arrays!(problem,model,specialNode)
end

apply_external_nodal_resultants!(problem)
return problem.residual
```

The negative sign is consistent with the existing element residual, where an
applied positive nodal resultant is subtracted from the equilibrium equation.
The UVLM load is frozen during the inner Newton solve, so this first version has
no aerodynamic Jacobian contribution. Add a derivative only if a future
monolithic or quasi-monolithic method requires it.

Place the type and its model-independent `clear!`, `add_wrench!`, and
`similar_external_loads` methods in `ExternalResultants.jl`. Include that file
before `Model.jl` in `src/AeroBeams.jl`, because the `Model` field refers to the
type. Place `create_ExternalNodalResultants(model::Model)` and the setter in
`Model.jl` after `Model` has been defined. Export only the public constructor and
setter; residual-assembly helpers can remain internal. Recreate the equation map
after any call to `update_model!`, because that call can change global node and
equation indices.

## 7. Public single-step structural API

`solve_dynamic!` currently initializes a problem and owns the complete time
loop. The coupler needs control of one physical time step. Expose small wrappers
around the existing internal operations rather than duplicating the structural
residual or Newton solver.

```julia
# src/Problem.jl -- proposed interface

function begin_dynamic_step!(problem::DynamicProblem,timeNow::Real,Δt::Real)
    @assert Δt > 0
    @assert timeNow > problem.timeNow

    problem.timeBeginTimeStep = problem.timeNow
    problem.timeEndTimeStep = timeNow
    problem.timeNow = timeNow
    problem.Δt = Δt

    update_basis_A_orientation!(problem)
    for BC in problem.model.BCs
        update_BC_data!(BC,timeNow)
    end

    # These equivalent rates define AeroBeams' implicit trapezoidal update.
    get_equivalent_states_rates!(problem)
    update_BL_complementary_variables!(problem)

    return problem
end

function solve_dynamic_step!(problem::DynamicProblem;
    initialGuess::Union{Nothing,Vector{Float64}}=nothing)

    if !isnothing(initialGuess)
        @assert length(initialGuess) == length(problem.x)
        problem.x .= initialGuess
    end

    solve_time_step!(problem)
    return problem.systemSolver.convergedFinalSolution
end

function commit_dynamic_step!(problem::DynamicProblem; saveSolution::Bool=true)
    if saveSolution
        save_time_step_data!(problem,problem.timeNow)
    end
    return problem
end
```

Extract initialization from `solve_dynamic!` into a reusable public helper:

```julia
function initialize_dynamic_problem!(problem::DynamicProblem)
    if !problem.skipInitialStatesUpdate
        solve_initial_dynamic!(problem)
    end
    if problem.saveInitialSolution
        save_time_step_data!(problem,problem.timeNow)
    end
    return problem
end
```

The existing high-level solver can then call `initialize_dynamic_problem!`
followed by its normal fixed or adaptive time loop, so this refactor does not
change behavior for current AeroBeams users.

Expose a public structural snapshot wrapper around the existing `copy_state` and
`restore_state!` functions:

```julia
snapshot_dynamic_state(problem::DynamicProblem) = copy_state(problem)

function restore_dynamic_state!(problem::DynamicProblem,snapshot)
    restore_state!(problem,snapshot)
    return problem
end
```

For fixed-step coupling, call `begin_dynamic_step!` once, then restore the
beginning-of-step structural state before every outer iteration. The equivalent
rates calculated by `begin_dynamic_step!` remain those of `t_n`.

Adaptive time stepping requires a richer snapshot containing time variables,
the model basis orientation, boundary-condition state, equivalent rate arrays,
and saved-history lengths. Implement and test fixed steps first; then make a
rejected coupled step restore both solvers completely.

## 8. Packaging the UVLM

### 8.1 Target module

Replace top-level scripts and global variables with a module and configuration
objects:

```julia
module WingPropellerUVLM

using LinearAlgebra
using Parameters
using StaticArrays

include("VortexLatticeBackend.jl")
include("UVLMConfiguration.jl")
include("UVLMAssembly.jl")
include("UVLMSnapshot.jl")
include("Geometry.jl")
include("Loads.jl")

export UVLMConfiguration,
       UVLMAssembly,
       create_UVLMAssembly,
       snapshot,
       restore!,
       set_geometry!,
       advance_trial!,
       commit!,
       surface_loads

end
```

All values currently read from `Main` or captured as global arrays become fields
of `UVLMConfiguration` or `UVLMAssembly`.

```julia
@with_kw struct UVLMConfiguration
    nSpanWing::Int
    nChordWing::Int
    nSpanPropeller::Int
    nChordPropeller::Int
    nBlades::Int
    nPropellers::Int
    wakeRowsWing::Int
    wakeRowsPropeller::Int
    interaction::Bool = true
    interactionID::Vector{Int}
    coreRadius::Function
end

@with_kw mutable struct UVLMAssembly{TS,TF,TR,TG}
    system::TS
    configuration::UVLMConfiguration
    freestream::TF
    repeatedPoints::TR
    referenceGeometry::TG
    propellerSurfaceIndices::Vector{Vector{Int}}
    activeWakeRows::Vector{Int}
    maximumWakeRows::Vector{Int}
    time::Float64 = 0.0
    step::Int = 0
end

function get_freestream(aero::UVLMAssembly,time::Real,index::Int)
    if aero.freestream isa Function
        return aero.freestream(time)
    end
    @assert 1 <= index <= length(aero.freestream)
    return aero.freestream[index]
end
```

Parameterizing the backend-dependent fields preserves type stability without
making the coupling layer depend on the concrete internal `System` type.

### 8.2 Complete aerodynamic snapshot

```julia
struct UVLMSnapshot{TW,TΓ,TDΓ,TΓT,TS}
    wakes::TW
    Γ::TΓ
    dΓ::TDΓ
    dΓdt::TΓT
    surfaces::TS
    previousSurfaces::TS
    activeWakeRows::Vector{Int}
    time::Float64
    step::Int
end

function snapshot(aero::UVLMAssembly)
    system = aero.system
    dΓ = getfield(system,Symbol("dΓ"))
    return UVLMSnapshot(
        deepcopy(system.wakes),
        copy(getfield(system,Symbol("Γ"))),
        tuple((copy(value) for value in dΓ)...),
        copy(getfield(system,Symbol("dΓdt"))),
        deepcopy(system.surfaces),
        deepcopy(system.previous_surfaces),
        copy(aero.activeWakeRows),
        aero.time,
        aero.step,
    )
end

function restore!(aero::UVLMAssembly,snapshot::UVLMSnapshot)
    system = aero.system
    system.wakes .= deepcopy(snapshot.wakes)
    system.surfaces .= deepcopy(snapshot.surfaces)
    system.previous_surfaces .= deepcopy(snapshot.previousSurfaces)
    getfield(system,Symbol("Γ")) .= snapshot.Γ
    getfield(system,Symbol("dΓdt")) .= snapshot.dΓdt

    dΓ = getfield(system,Symbol("dΓ"))
    for i in eachindex(dΓ)
        dΓ[i] .= snapshot.dΓ[i]
    end

    aero.activeWakeRows .= snapshot.activeWakeRows
    aero.time = snapshot.time
    aero.step = snapshot.step
    return aero
end
```

Depending on the final backend, `system.wakes .= ...` may need to be replaced by
an element-wise assignment because wake matrices can change size. The rollback
test, not superficial field equality, is decisive: restore a snapshot, execute
the same trial twice, and require identical circulation, wake geometry, and
loads.

### 8.3 Trial propagation

The existing UVLM propagation call should be wrapped directly:

```julia
function advance_trial!(aero::UVLMAssembly,freestream,Δt::Real)
    @assert Δt > 0
    config = aero.configuration

    # Existing Wing_Propeller_UVLM call.
    propagate_system!(
        aero.system,
        freestream,
        Δt;
        additional_velocity = nothing,
        repeated_points = aero.repeatedPoints,
        nwake = aero.activeWakeRows,
        eta = 0.1,
        calculate_influence_matrix = true,
        near_field_analysis = true,
        derivatives = false,
        interaction_id = config.interactionID,
        interaction = config.interaction,
    )

    aero.time += Δt
    aero.step += 1
    return aero
end
```

`calculate_influence_matrix` must be `true` whenever structural deformation or
propeller azimuth changes the panel geometry. Reusing an influence matrix is
valid only after proving that the relevant geometry is unchanged.

`commit!` advances active wake counters once and establishes the accepted
surface as the previous surface for the following physical step:

```julia
function commit!(aero::UVLMAssembly)
    for i in eachindex(aero.activeWakeRows)
        maximumRows = aero.maximumWakeRows[i]
        aero.activeWakeRows[i] = min(aero.activeWakeRows[i]+1,maximumRows)
    end
    return aero
end
```

Adapt the maximum-wake lookup to the actual `System` fields. Keep this operation
out of `advance_trial!` so rejected trials never increase the active wake.

## 9. Structural interface mesh

The coupler should not assume that UVLM span stations always coincide with beam
nodes. Build the interpolation topology once.

```julia
struct SpanStencil
    node1::Int
    node2::Int
    N1::Float64
    N2::Float64
end

@with_kw struct WingInterfaceMesh
    spanCoordinates::Vector{Float64}
    stencils::Vector{SpanStencil}
    referenceOffsets_B::Vector{Matrix{Float64}}
end

@with_kw struct PropellerInterface
    pylonTipNode::Int
    rotorAxis_H::Vector{Float64} = [1.0;0.0;0.0]
    referenceBladeGrids_H::Vector{Array{Float64,3}}
    spinDirection::Float64 = -1.0
end
```

For an aerodynamic station at beam coordinate `s` in element
`[s1,s2]`, use

```julia
ξ = (s-s1)/(s2-s1)
stencil = SpanStencil(node1,node2,1-ξ,ξ)
```

The same stencil must be used by motion transfer and its transpose by load
transfer. This is the basis of discrete work conservation.

## 10. Extracting AeroBeams nodal kinematics

Only one copy of a shared global node should enter the interface arrays.

```julia
@with_kw mutable struct StructuralKinematics
    positions_A::Matrix{Float64}
    rotations_AB::Vector{Matrix{Float64}}
    assigned::BitVector
end

function extract_nodal_kinematics(model::Model)
    kinematics = StructuralKinematics(
        positions_A = zeros(3,model.nNodesTotal),
        rotations_AB = [Matrix(I,3,3) for _ in 1:model.nNodesTotal],
        assigned = falses(model.nNodesTotal),
    )

    for element in model.elements
        data = (
            (element.nodesGlobalID[1],element.r_n1,
             element.nodalStates.u_n1,element.nodalStates.p_n1,
             element.R0_n1),
            (element.nodesGlobalID[2],element.r_n2,
             element.nodalStates.u_n2,element.nodalStates.p_n2,
             element.R0_n2),
        )

        for (node,r0_A,u_A,p_A,R0_AB) in data
            if !kinematics.assigned[node]
                R_deformation = rotation_tensor_WM(p_A)[1]
                kinematics.positions_A[:,node] .= r0_A .+ u_A
                kinematics.rotations_AB[node] .= R_deformation * R0_AB
                kinematics.assigned[node] = true
            end
        end
    end

    @assert all(kinematics.assigned)
    return kinematics
end
```

Call `update_states!(problem)` or complete a residual assembly before extracting
kinematics so that `element.nodalStates` corresponds to `problem.x`.

## 11. Motion transfer

### 11.1 Wing grid

At a wing interface station, interpolate the reference-line translation with
the span stencil. Do not linearly interpolate the nine entries of two rotation
matrices. Use a quaternion/rotation-vector interpolation, or make UVLM stations
coincide with structural nodes in the first implementation.

For the aligned-grid first version:

```julia
function update_wing_grid!(grid_U,interface::WingInterfaceMesh,
    kinematics::StructuralKinematics,frames::CouplingFrames)

    for (j,stencil) in enumerate(interface.stencils)
        @assert stencil.node1 == stencil.node2 ||
            (stencil.N1 == 1.0 || stencil.N2 == 1.0)

        node = stencil.N1 == 1.0 ? stencil.node1 : stencil.node2
        r_A = kinematics.positions_A[:,node]
        R_AB = kinematics.rotations_AB[node]

        offsets_B = interface.referenceOffsets_B[j]
        for i in axes(offsets_B,2)
            grid_U[:,i,j] .= frames.R_UA * (r_A + R_AB*offsets_B[:,i])
        end
    end

    return grid_U
end
```

The offsets are measured from the structural reference line, normally the
elastic axis, to each chordwise UVLM vertex in the undeformed section basis.
This automatically includes bending, twist, sweep, dihedral, and large rigid
rotations without reconstructing Euler angles.

For nonmatching grids, use

```math
\boldsymbol r_U(s,c)=R_{UA}
\left[\boldsymbol r_A(s)+R_{AB}(s)\boldsymbol d_B(s,c)\right]
```

with `r_A(s)` interpolated by shape functions and `R_AB(s)` interpolated on
`SO(3)`.

### 11.2 Propeller grids

With a physical pylon, its tip supplies hub translation and whirl orientation.
Do not add separate pitch/yaw deformation angles from the old two-DOF model.
Only prescribed shaft rotation is superposed:

```math
R_{UH}(t)=R_{UA}R_{AH},
\qquad
R_{U\mathcal{B}_k}(t)=R_{UH}(t)R_{\mathrm{spin}}(\phi)R_{k0}.
```

```julia
function update_propeller_grids!(grids_U,propeller::PropellerInterface,
    kinematics::StructuralKinematics,frames::CouplingFrames,φ::Real)

    node = propeller.pylonTipNode
    rHub_U = frames.R_UA * kinematics.positions_A[:,node]
    R_UH = frames.R_UA * kinematics.rotations_AB[node]
    spinAngle = propeller.spinDirection*φ
    sine,cosine = sincos(spinAngle)
    R_spin = [1.0 0.0     0.0;
              0.0 cosine -sine;
              0.0 sine    cosine]

    for blade in eachindex(grids_U)
        referenceGrid_H = propeller.referenceBladeGrids_H[blade]
        for j in axes(referenceGrid_H,2), k in axes(referenceGrid_H,3)
            point_H = referenceGrid_H[:,j,k]
            grids_U[blade][:,j,k] .= rHub_U + R_UH*R_spin*point_H
        end
    end

    return grids_U
end
```

The sample `R_spin` call illustrates composition but must use the packaged UVLM
rotation helper or a dedicated axis-angle function after the hub-axis convention
is finalized. Test one positive azimuth increment against the current driver to
confirm the spin sign.

### 11.3 Surface history

Surface velocity in the UVLM is derived from current and previous panels. At the
beginning of a physical step:

1. retain the accepted `t_n` geometry as `previous_surfaces`;
2. replace only `surfaces` with the trial `t_(n+1)` geometry;
3. restore both arrays before another trial;
4. after convergence, leave the accepted `surfaces` in place.

Never update `previous_surfaces` for every outer coupling iteration.

## 12. Conservative aerodynamic load transfer

### 12.1 Wing resultants

For each aerodynamic span station, first sum all chordwise forces and take their
moment about the deformed beam reference line:

```math
\boldsymbol F_{U,j}=\sum_i \boldsymbol f_{U,ij},
```

```math
\boldsymbol M_{U,j}=\sum_i
\left[\boldsymbol m_{U,ij}+
(\boldsymbol r_{U,ij}-\boldsymbol r_{U,EA,j})
\times\boldsymbol f_{U,ij}\right].
```

Transform the complete wrench to `A` and distribute it with the transpose of
the motion interpolation:

```julia
function transfer_station_wrench!(loads::ExternalNodalResultants,
    stencil::SpanStencil,frames::CouplingFrames,
    force_U::AbstractVector,moment_U::AbstractVector)

    force_A = frames.R_AU * force_U
    moment_A = frames.R_AU * moment_U

    add_wrench!(loads,stencil.node1,
        stencil.N1*force_A,stencil.N1*moment_A)
    add_wrench!(loads,stencil.node2,
        stencil.N2*force_A,stencil.N2*moment_A)

    return loads
end
```

If a stencil collapses to one aligned node, set the other weight to zero or
store a single-node stencil without adding the same node twice.

The current UVLM implementation exposes span-segment, chord-segment, and
unsteady forces. Audit their definitions once during packaging and establish
one function that returns the total physical panel or vertex load. Do not sum
two arrays unless the UVLM force derivation confirms that they are additive and
non-overlapping.

### 12.2 Propeller hub wrench

Reduce all blade loads at the physical hub center:

```julia
function propeller_hub_wrench(panelPositions_U,panelForces_U,
    panelMoments_U,rHub_U)

    force_U = zeros(3)
    moment_U = zeros(3)
    for i in eachindex(panelForces_U)
        force_U .+= panelForces_U[i]
        moment_U .+= panelMoments_U[i] +
            cross(panelPositions_U[i]-rHub_U,panelForces_U[i])
    end
    return force_U,moment_U
end
```

Apply the full wrench to the pylon-tip global node:

```julia
force_A = frames.R_AU * forceHub_U
moment_A = frames.R_AU * momentHub_U
add_wrench!(loads,pylonTipNode,force_A,moment_A)
```

Do not also apply the same propeller wrench directly to the wing attachment.
The pylon internal resultants transmit it to the wing.

### 12.3 Conservation tests

For every transfer call, test

```math
\sum_n\boldsymbol F_n=\sum_a\boldsymbol f_a,
```

```math
\sum_n\left(\boldsymbol M_n+oldsymbol r_n\times\boldsymbol F_n\right)
=
\sum_a\left(\boldsymbol m_a+oldsymbol r_a\times\boldsymbol f_a\right),
```

and, for a random admissible virtual displacement,

```math
\delta\boldsymbol q_s^T\boldsymbol Q_s
=
\delta\boldsymbol x_a^T\boldsymbol f_a.
```

Force/moment conservation can pass while virtual work fails, so all three are
needed.

## 13. Propeller structural equations

### 13.1 Prescribed-speed rotor

Let `e_s` be the shaft direction, `Ω` its prescribed speed, and `φ` its phase:

```math
\dot\phi=\Omega,
\qquad
\phi_{n+1}=\phi_n+\Omega\Delta t.
```

If the rotor mass and complete non-spinning inertia tensor are already included
as an AeroBeams `PointInertia`, add only the angular momentum associated with
relative shaft spin:

```math
\boldsymbol h_{spin}=J_s\Omega\boldsymbol e_s.
```

The additional reaction moment on the support is based on

```math
\boldsymbol m_{spin}
=
\left.\frac{d\boldsymbol h_{spin}}{dt}\right|_A
=J_s\dot\Omega\boldsymbol e_s
+\boldsymbol\omega_h\times\boldsymbol h_{spin}.
```

For constant `Ω`, only the gyroscopic cross product remains. Do not add the
complete rigid-body inertia again in a `SpinningRotor` component.

An initial partitioned implementation may evaluate this moment from the latest
hub angular velocity and freeze it with the other external resultants during the
inner structural Newton solve:

```julia
@with_kw mutable struct SpinningRotor
    node::Int
    axialInertia::Float64
    speed::Function
    acceleration::Function = t -> 0.0
    axis_H::Vector{Float64} = [1.0;0.0;0.0]
    phase::Float64 = 0.0
end

function spinning_rotor_moment_A(rotor::SpinningRotor,R_AH,
    ωHub_A::AbstractVector,time::Real)

    Ω = rotor.speed(time)
    Ωdot = rotor.acceleration(time)
    axis_A = R_AH * rotor.axis_H
    hSpin_A = rotor.axialInertia * Ω * axis_A
    return rotor.axialInertia*Ωdot*axis_A + cross(ωHub_A,hSpin_A)
end

function commit_rotor_phase!(rotor::SpinningRotor,time_n::Real,time_np1::Real)
    Δt = time_np1-time_n
    @assert Δt > 0
    rotor.phase += Δt/2 * (rotor.speed(time_n)+rotor.speed(time_np1))
    return rotor
end
```

The nonlinear production version should assemble this support moment and its
consistent tangent as a special-node contribution. That lets AeroBeams Newton
iterations see the dependence on nodal rotation and angular velocity. The
partitioned version remains useful as a regression implementation.

### 13.2 Variable-speed rotor

If RPM becomes a dynamic state, add one scalar shaft equation per propeller:

```math
J_s\dot\Omega=
\tau_{motor}(t,\Omega)-Q_{aero}-Q_{loss}(\Omega),
\qquad \dot\phi=\Omega.
```

Here `Q_aero` is the component of the UVLM hub moment along the instantaneous
shaft axis. Define its sign from power:

```math
P_{aero}=Q_{aero}\Omega.
```

The rotor phase must be included in both structural and aerodynamic snapshots
so a rejected step cannot advance azimuth.

### 13.3 Baseline two-degree-of-freedom option

To reproduce the old linear solution before using a physical pylon, a custom
`WhirlPropeller` can retain pitch and yaw coordinates with the current
`M`, `C`, `K`, unbalance, and gyroscopic blocks. Treat this as a temporary
verification model. It requires extending AeroBeams global degrees of freedom
and residual assembly, whereas the physical pylon beam uses AeroBeams native
states. Do not combine the two representations.

### 13.4 Damping

AeroBeams 0.8.1 does not provide general structural damping. If the pylon mount
damping ratio in the old model is physically required, implement it explicitly
as a nodal/pylon constitutive component or a carefully documented Rayleigh
damping extension. Numerical dissipation in a time integrator is not a physical
replacement for pylon damping in whirl-flutter calculations.

## 14. Coupling data types

```julia
# src/Coupling/CoupledDynamicProblem.jl

@with_kw struct CouplingSolver
    maxIterations::Int = 25
    positionTolerance::Float64 = 1e-6
    rotationTolerance::Float64 = 1e-7
    loadTolerance::Float64 = 1e-5
    initialRelaxation::Float64 = 0.5
    minimumRelaxation::Float64 = 0.05
    maximumRelaxation::Float64 = 1.0
end

@with_kw mutable struct CouplingIteration
    iteration::Int = 0
    relaxation::Float64
    positionResidual::Float64 = Inf
    rotationResidual::Float64 = Inf
    loadResidual::Float64 = Inf
    converged::Bool = false
end

@with_kw mutable struct CoupledDynamicProblem{TP,TA,TM,TL,TS}
    structuralProblem::TP
    aerodynamicSystem::TA
    motionTransfer::TM
    loadTransfer::TL
    couplingSolver::TS
    externalLoads::ExternalNodalResultants
    rotors::Vector{SpinningRotor}
    timeVector::Vector{Float64}
    interfaceHistory::Vector{Any} = Any[]
    loadHistory::Vector{ExternalNodalResultants} = ExternalNodalResultants[]
    iterationHistory::Vector{CouplingIteration} = CouplingIteration[]
end
```

Replace the `Any` interface history with the final concrete interface-state type
once the motion mapper is implemented. It appears here only because that type
depends on the selected rotation representation.

Use a dedicated interface state instead of relaxing AeroBeams' complete mixed
vector `x`. That vector contains displacements, rotations, internal forces,
moments, velocities, and possibly aerodynamic states with different units.
Relaxing all of it componentwise can violate the structural compatibility
equations. The relaxed interface state is used only to generate the next UVLM
geometry; every accepted structural state remains a fully converged AeroBeams
solution.

## 15. Aitken relaxation

Form an interface residual from translation and a proper relative-rotation
measure. For two rotation matrices, use the rotation vector of
`R_candidate*R_guess'`; do not subtract Euler angles.

For a flattened, nondimensional interface residual `r_k`, dynamic Aitken
relaxation is

```julia
function aitken_relaxation(ωPrevious,rPrevious,rCurrent,solver::CouplingSolver)
    Δr = rCurrent-rPrevious
    denominator = dot(Δr,Δr)
    if denominator <= eps(Float64)
        return ωPrevious
    end

    ω = -ωPrevious * dot(rPrevious,Δr) / denominator
    return clamp(ω,solver.minimumRelaxation,solver.maximumRelaxation)
end
```

Nondimensionalize translations by a reference length and loads by reference
force and moment. Report position, rotation, and load residuals separately even
if a combined vector is used for Aitken relaxation.

## 16. Strongly coupled time-step algorithm

The following is the intended controller. Helper functions such as
`extract_interface_state`, `update_geometry!`, and `transfer_loads!` belong to
the coupling layer.

```julia
function solve_coupled_step!(coupled::CoupledDynamicProblem,timeIndex::Int)
    @unpack structuralProblem,aerodynamicSystem,couplingSolver = coupled

    time_n = coupled.timeVector[timeIndex-1]
    time_np1 = coupled.timeVector[timeIndex]
    Δt = time_np1-time_n

    structuralSnapshot = snapshot_dynamic_state(structuralProblem)
    aerodynamicSnapshot = snapshot(aerodynamicSystem)

    begin_dynamic_step!(structuralProblem,time_np1,Δt)

    interface_n = extract_interface_state(structuralProblem.model)
    interfaceGuess = predict_interface_state(coupled,interface_n,Δt)
    previousResidual = nothing
    previousLoads = nothing
    relaxation = couplingSolver.initialRelaxation

    iterationData = CouplingIteration(relaxation=relaxation)

    for iteration in 1:couplingSolver.maxIterations
        iterationData.iteration = iteration

        # All aerodynamic trials begin with the accepted history at t_n.
        restore!(aerodynamicSystem,aerodynamicSnapshot)
            update_geometry!(coupled.motionTransfer,aerodynamicSystem,
                interfaceGuess,time_np1)
            advance_trial!(aerodynamicSystem,
            get_freestream(aerodynamicSystem,time_np1,timeIndex),Δt)

        clear!(coupled.externalLoads)
        transfer_loads!(coupled.loadTransfer,coupled.externalLoads,
            aerodynamicSystem,interfaceGuess)
        add_rotor_reactions!(coupled.externalLoads,coupled.rotors,
            interfaceGuess,time_np1)
        set_external_nodal_resultants!(structuralProblem.model,
            coupled.externalLoads)

        # Every structural trial also begins at t_n. begin_dynamic_step! was
        # called once, so its trapezoidal equivalent rates remain prepared.
        restore_dynamic_state!(structuralProblem,structuralSnapshot)
        convergedStructure = solve_dynamic_step!(structuralProblem)
        convergedStructure || error(
            "AeroBeams Newton solve failed at t=$time_np1, iteration=$iteration")

        interfaceCandidate = extract_interface_state(structuralProblem.model)
        residual = interface_residual(interfaceGuess,interfaceCandidate)

        iterationData.positionResidual = position_residual(residual)
        iterationData.rotationResidual = rotation_residual(residual)
        iterationData.loadResidual = isnothing(previousLoads) ? Inf :
            relative_load_change(coupled.externalLoads,previousLoads)

        motionConverged =
            iterationData.positionResidual <= couplingSolver.positionTolerance &&
            iterationData.rotationResidual <= couplingSolver.rotationTolerance
        loadConverged = iteration > 1 &&
            iterationData.loadResidual <= couplingSolver.loadTolerance

        if motionConverged && loadConverged
            iterationData.converged = true

            # Recompute at the accepted structural geometry from the same
            # snapshot. This is the only aerodynamic state retained.
            restore!(aerodynamicSystem,aerodynamicSnapshot)
            update_geometry!(coupled.motionTransfer,aerodynamicSystem,
                interfaceCandidate,time_np1)
            advance_trial!(aerodynamicSystem,
                get_freestream(aerodynamicSystem,time_np1,timeIndex),Δt)

            committedLoads = similar_external_loads(coupled.externalLoads)
            transfer_loads!(coupled.loadTransfer,committedLoads,
                aerodynamicSystem,interfaceCandidate)
            add_rotor_reactions!(committedLoads,coupled.rotors,
                interfaceCandidate,time_np1)

            commitMismatch = relative_load_change(
                committedLoads,coupled.externalLoads)
            commitMismatch <= couplingSolver.loadTolerance || error(
                "Committed UVLM loads are inconsistent at t=$time_np1")

            coupled.externalLoads = committedLoads
            set_external_nodal_resultants!(structuralProblem.model,
                coupled.externalLoads)
            commit!(aerodynamicSystem)
            for rotor in coupled.rotors
                commit_rotor_phase!(rotor,time_n,time_np1)
            end
            commit_dynamic_step!(structuralProblem)
            push!(coupled.iterationHistory,deepcopy(iterationData))
            return true
        end

        if !isnothing(previousResidual)
            relaxation = aitken_relaxation(
                relaxation,previousResidual,residual,couplingSolver)
        end

        interfaceGuess = relax_interface(
            interfaceGuess,interfaceCandidate,relaxation)
        previousResidual = residual
        previousLoads = deepcopy(coupled.externalLoads)
        iterationData.relaxation = relaxation
    end

    # Reject the physical step without changing either subsystem.
    restore_dynamic_state!(structuralProblem,structuralSnapshot)
    restore!(aerodynamicSystem,aerodynamicSnapshot)
    push!(coupled.iterationHistory,deepcopy(iterationData))
    return false
end
```

The production version should return a typed status object rather than throw for
normal step rejection. Errors are shown above to make failure points explicit.
If adaptive stepping is enabled, a rejected step halves `Δt`, restores the full
coupled snapshot, and retries without advancing rotor phase or wake length.

## 17. Top-level initialization and call sequence

After the UVLM is packaged and the generic AeroBeams hooks are implemented, an
example should look like this:

```julia
using AeroBeams
using WingPropellerUVLM
using AeroBeamsUVLMCoupling

# Structural model
model,wing,pylons = create_chang_aerobeams_model(structuralData)

Δψ = deg2rad(2.5)
Δt = Δψ/abs(rotorSpeed)
timeVector = collect(0.0:Δt:finalTime)

newton = create_NewtonRaphson(
    maximumIterations = 50,
    relativeTolerance = 1e-8,
)

structuralProblem = create_DynamicProblem(
    model = model,
    systemSolver = newton,
    timeVector = timeVector,
    trackingTimeSteps = false,
    displayProgress = false,
)

# Solve consistent structural initial states using an AeroBeams public helper
# before handing control of subsequent steps to the coupled driver.
initialize_dynamic_problem!(structuralProblem)

# Global wing + all propeller blades
uvlm = create_UVLMAssembly(
    configuration = uvlmConfiguration,
    wingGeometry = wingGeometry,
    propellerGeometries = propellerGeometries,
    freestream = freestreamHistory,
)

frames = create_CouplingFrames()
motionTransfer = create_MotionTransfer(
    model = model,
    frames = frames,
    wingStations = wingStations,
    propellerTipNodes = [p.nodeRange[end] for p in pylons],
    referenceWingGrid = wingGeometry.grid,
    referenceBladeGrids = propellerGeometries,
)
loadTransfer = create_LoadTransfer(motionTransfer)
externalLoads = create_ExternalNodalResultants(model)
set_external_nodal_resultants!(model,externalLoads)

rotors = [SpinningRotor(
    node = pylon.nodeRange[end],
    axialInertia = rotorData[ip].axialInertia,
    speed = t -> rotorSpeed,
) for (ip,pylon) in enumerate(pylons)]

coupled = CoupledDynamicProblem(
    structuralProblem = structuralProblem,
    aerodynamicSystem = uvlm,
    motionTransfer = motionTransfer,
    loadTransfer = loadTransfer,
    couplingSolver = CouplingSolver(),
    externalLoads = externalLoads,
    rotors = rotors,
    timeVector = timeVector,
)

for timeIndex in 2:length(timeVector)
    accepted = solve_coupled_step!(coupled,timeIndex)
    accepted || error("Coupling failed at time index $timeIndex")
end
```

`initialize_dynamic_problem!` is another small public wrapper to extract from
the initialization portion of the current `solve_dynamic!`. It should establish
consistent initial structural states without entering AeroBeams' own time loop.

## 18. Steady equilibrium and initialization

A rotating propeller normally gives a periodic load rather than an instantaneous
steady load. A nonlinear static equilibrium should use a periodic aerodynamic
state and cycle-averaged resultants:

1. Hold the current structural geometry fixed.
2. Run the UVLM for enough rotor revolutions to remove wake startup transients.
3. Average loads over an integer number of converged revolutions.
4. Solve the AeroBeams steady nonlinear problem with those loads.
5. Update the UVLM geometry and repeat until deformation and averaged loads
   converge.

Use the resulting deformed structure, circulation, and periodic wake as the
initial state for a time-domain perturbation. Capturing a single instantaneous
load as `F0` does not guarantee a periodic trim and can contaminate the measured
modal damping with startup transients.

## 19. Stability analysis

AeroBeams' standard eigenproblem linearizes its structural equations and local
aerodynamic states. It does not automatically include a global UVLM free wake.
Use this development sequence:

1. Apply a small smooth structural or propeller-whirl perturbation.
2. Run the fully coupled time-domain system.
3. Identify modal frequency and damping from several response channels.
4. Repeat with smaller perturbation to verify linear behavior.
5. Later, finite-difference the complete discrete coupled map to obtain a
   state-transition Jacobian.
6. When the base solution is rotor-periodic, use a one-revolution monodromy
   matrix and Floquet multipliers.

Do not pass UVLM-coupled cases directly to `create_EigenProblem` and interpret
the output as coupled free-wake stability without this additional linearization.

## 20. Verification plan and acceptance criteria

### 20.1 Unit tests

| Test | Required result |
|:--|:--|
| Frame orthogonality | `R'R = I`, `det(R)=+1` |
| Basis directions | chord, span, and vertical mappings match the declared convention |
| Rigid translation | every aerodynamic point receives exactly the same translation |
| Rigid rotation | grid matches direct rotation of the undeformed geometry |
| Rotor azimuth | one positive increment matches the current UVLM driver |
| Force conservation | structural and aerodynamic total forces agree |
| Moment conservation | totals about an arbitrary common origin agree |
| Virtual work | structural and aerodynamic work agree for random admissible perturbations |
| UVLM rollback | repeated trials from one snapshot are bitwise equal where practical |
| Wake transaction | rejected trials do not change active wake length |

### 20.2 Structural regression

- Match total mass, center of mass, and inertia tensor.
- Match static tip displacement and twist in the small-load limit.
- Match the first wing and pylon modes of the current linear model.
- Reduce load amplitude and confirm convergence to the linear response.
- Increase load and demonstrate a genuine geometric-nonlinearity effect.

Suggested initial modal acceptance is within one or two percent after both
models use equivalent property distributions and boundary conditions. Tighten
this based on discretization studies rather than treating it as a universal
tolerance.

### 20.3 Aerodynamic regression

- Isolated rigid wing lift and moment.
- Isolated propeller thrust and torque coefficients.
- Wing--propeller interaction on and off.
- Circulation and wake geometry for one time step from a saved initial state.
- Surface-force decomposition without missing or double-counted terms.

### 20.4 Coupled limiting cases

- Zero freestream and zero rotor speed reproduce structural-only AeroBeams.
- Infinite structural stiffness reproduces rigid UVLM.
- Interaction disabled reproduces the sum of isolated components.
- Zero propeller load reproduces wing-only aeroelasticity.
- Zero wing load reproduces pylon/rotor structural dynamics.
- Zero rotor spin removes gyroscopic terms.
- Reversing spin reverses gyroscopic cross-coupling but not non-spinning mass.

### 20.5 Discretization and coupling convergence

Perform independent studies of:

- AeroBeams elements;
- wing chordwise and spanwise UVLM panels;
- propeller radial and chordwise panels;
- wake length;
- vortex-core model;
- rotor azimuth step, for example `10`, `5`, and `2.5` degrees;
- outer coupling tolerances;
- maximum coupling iterations.

Do not compensate for a nonconverged outer iteration by reducing only the
structural Newton tolerance. They solve different residuals.

## 21. Performance rules

Correctness takes precedence during the first coupled implementation. After the
verification ladder passes:

- preallocate grid, wrench, and interface-residual arrays;
- reuse span stencils and reference offsets;
- avoid `deepcopy` in inner loops by implementing explicit buffer copies;
- keep concrete parametric field types;
- profile before changing the UVLM algorithms;
- rebuild the influence matrix whenever geometry changes;
- consider aerodynamic subcycling only after a single-rate solution converges;
- do not differentiate through the mutable free-wake solver initially.

The propeller azimuth requirement will often impose a much smaller time step
than the structural dynamics. A later multirate algorithm may execute several
UVLM substeps per structural step, but it requires a new definition of load
averaging and a rollback transaction for all substeps. It should not be part of
the first coupled milestone.

## 22. Coding and review checklist

Before merging each stage:

- no executable work occurs at module load time;
- no physical parameter is read from `Main`;
- every public constructor validates dimensions, units, and coordinate basis;
- each mutating function returns its mutated object or `nothing` consistently;
- units are SI at the coupling boundary;
- snapshots contain all history needed for deterministic rollback;
- trial and commit operations have separate names;
- force and moment reference points are stated in docstrings;
- the structural and aerodynamic solvers can still run independently;
- tests do not depend on plotting or interactive state;
- regression data records repository commit, mesh, time step, and tolerances;
- the current generalized-`α` linear result remains archived as a reference,
  even though the first AeroBeams coupling uses its existing trapezoidal rule.

## 23. Implementation milestones

1. **Package the UVLM.** Move configuration and mutable state out of the current
   driver; add snapshot/restore determinism tests.
2. **Build the AeroBeams structure.** Match mass, static response, and modes with
   no aerodynamics.
3. **Add generic external resultants and one-step APIs.** Prove a prescribed
   nodal load gives the same result through the existing and new interfaces.
4. **Implement coordinate and motion transfer.** Pass rigid-body and large-
   rotation grid tests.
5. **Implement conservative load transfer.** Pass force, moment, and virtual-
   work tests.
6. **Run one-way cases.** AeroBeams motion to UVLM, then prescribed UVLM loads to
   AeroBeams, separately.
7. **Run strong two-way coupling.** Use fixed time steps, prescribed RPM,
   rollback, and Aitken relaxation.
8. **Add physical rotor effects.** Add relative-spin gyroscopic reaction without
   duplicating `PointInertia` contributions.
9. **Establish periodic trim.** Initialize the dynamic calculation from a
   deformed, periodic aerodynamic state.
10. **Perform stability and convergence studies.** Only after all limiting-case
    tests pass.

This ordering preserves a runnable, validated configuration at every milestone
and isolates structural, aerodynamic, transfer, and time-integration errors.
