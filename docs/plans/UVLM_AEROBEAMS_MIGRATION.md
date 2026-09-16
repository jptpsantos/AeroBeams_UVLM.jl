> Implementation update (15 September 2026): the user authorized the terminology migration. The working facade APIs and shared convergence package are described in [the self-contained analysis guide](../migration/UVLM_ANALYSIS_GUIDE.tex). That guide identifies implemented behavior and remaining items; the P0--P11 list below remains the broader roadmap, including the future native AeroBeams structural bridge.
**UVLM migration to AeroBeams terminology and structure — implementation plan**

Prepared 14 September 2026 against working tree based on `d0fe40bf7f58f174c49cac62c62e424e32b41586`, AeroBeams 0.8.1 and WingPropellerUVLM 0.1.0. This is a proposed implementation sequence. Names marked **proposed** below are design targets, not callable APIs today. No solver refactoring is performed by this document.

Adopt AeroBeams' model–problem–solver organization and constructor conventions, retain a standalone UVLM package, and put the actual structural coupling in a small bridge package. Complete the organization and compatibility work before replacing the Chang structural solver. Changes to physics, time integration, and numerical approximations receive their own validation gates.

The intended outcome is that a user builds aerodynamic and structural models, creates a problem with explicit solver and time settings, calls `solve!`, and inspects named results. The global circulation and wake remain aerodynamic state; AeroBeams owns the structural state and structural time integration in coupled analyses.

**1. Scope, existing evidence, and decisions**

The scope covers `lib/WingPropellerUVLM/src`, the reusable machinery currently under `examples/chang_linear_aeroelastic`, convergence studies, package environments, and the minimal AeroBeams interfaces needed for coupling. Preserve the backend's numerical formulas during the initial migration. Keep Chang geometry tables, fitted stiffnesses, modal assumptions, and experiment-specific defaults in benchmark builders.

Use these decisions as the implementation baseline:

| Decision | Selected direction | Reason |
|---|---|---|
| Package boundary | Keep `WingPropellerUVLM`; add proposed `lib/AeroBeamsUVLM` bridge | Each discipline remains independently usable; one explicit place owns the mapping |
| Dependency direction | Bridge depends on AeroBeams and WingPropellerUVLM; neither discipline depends on the bridge | Avoid a circular dependency and implicit loading through examples |
| Public naming | `UVLMModel`, `UVLMDynamicProblem`, `create_Type`, snake_case functions, camelCase keywords | Follow local AeroBeams conventions while avoiding ambiguous exported `Model`/`DynamicProblem` names |
| First structural integration | Partitioned coupling with frozen dimensional external loads during each inner structural Newton solve | Reuse the current structural residual and solver without claiming an unavailable global aerodynamic Jacobian |
| Time ownership | Standalone UVLM problem owns its clock; the contained AeroBeams `DynamicProblem` owns the coupled clock | Eliminate independent physical-step advancement |
| Aerodynamic state | One accepted global state, one trial state, explicit commit/rollback | Preserve nonlocal interactions and exactly one wake advancement per accepted step |
| Damping | Retain the user-selected moving-block reference method | Refactoring must preserve the chosen estimator and its diagnostics |
| Compatibility | Translate old inputs at public boundaries and retain old output readers | Existing studies and archived results remain interpretable |
| First coupled capability | Fixed-step dynamics, prescribed rotor speed, restricted supported attachments | Establish a small validated interface before adding adaptive stepping, hinges, or linearization |

The preceding [implementation review](../reviews/uvlm_aerobeams_2026-09-14/REVIEW.md) is historical evidence. Two findings have since changed: the aeroelastic entry now uses `MovingBlockDamping.jl`, and its nested `Convergence.run_aeroelastic(settings)` delegates to the entry's speed-override path. The recorded 124 targeted checks passed after those changes; they do not establish that all earlier findings are resolved. FFT-bin resolution, raw-positive-peak offset sensitivity, dependency isolation, termination evidence, and broad fingerprints still need explicit treatment.

The older [coupling blueprint](../src/wing-propeller-uvlm-coupling.md) contains useful derivations and prototype interfaces. Use this plan for sequence, package ownership, and current terminology. Review prototype code before promoting it: its `UVLMAssembly` target becomes the model/state/workspace separation below, and its external-resultant equation map needs explicit restrictions for hinges and constrained nodes. The existing [library guide](../src/wing-propeller-uvlm-library-guide.md) must be updated with each completed public API phase.

**2. Conventions to copy, with their actual meanings**

The local source is the authority for this checkout. AeroBeams distinguishes steady, trim, eigenvalue, and dynamic problems; its official project describes the same analysis categories. These names must preserve those distinctions. [Official AeroBeams overview](https://github.com/luizpancini/AeroBeams.jl#features).

| Local AeroBeams evidence | Rule for UVLM |
|---|---|
| [Model.jl](../../src/Model.jl), `create_Model`; [Beam.jl](../../src/Beam.jl), `create_Beam` | Model constructors collect physical inputs, validate them, and build derived topology |
| [Problem.jl](../../src/Problem.jl), `create_DynamicProblem` | Problem holds model, time grid, solver choice, current solve state, and tracked output |
| [SystemSolver.jl](../../src/SystemSolver.jl), `create_NewtonRaphson` | Algorithm and tolerance configuration belong in a solver object |
| [Element.jl](../../src/Element.jl), `ElementalStates` | Preserve the meaning and basis of `u,p,F,M,V,Ω,χ`; do not rename an unrelated vector to these symbols |
| [Utilities.jl](../../src/Utilities.jl), `rotation_tensor_WM` | Use AeroBeams' deformation-rotation reconstruction and associated tangent operators |
| [Problem.jl](../../src/Problem.jl), `trackingTimeSteps`, `trackingFrequency`, `savedTimeVector` | Expose predictable sampling and storage controls with corresponding time coordinates |
| [UnitsSystem.jl](../../src/UnitsSystem.jl), constructor docstring | Unit labels only affect plotting; implement actual input conversions at the bridge boundary |
| [AeroSolver.jl](../../src/AeroSolver.jl) and [AeroSurface.jl](../../src/AeroSurface.jl) | Existing solvers describe element-local aerodynamic models; a UVLM marker alone cannot provide the required global solve |

Use PascalCase for types/files, `create_Type(; kwargs...)` for public keyword constructors, snake_case for operations, and `!` when an operation changes arguments. Use camelCase for public descriptive fields and keywords, with established mathematical symbols where their physical meaning is clear. Private numerical kernels may retain short mathematical names and existing backend names during migration.

Each public type needs documentation of inputs, derived fields, mutable fields, units, coordinate basis, ownership, and validation. Use concrete or parameterized hot-path storage; copying broad `Real`/`Function` fields or adding `@pack!` everywhere is not a requirement for compatibility.

Do not immediately export a second unqualified `solve!` beside AeroBeams' function in examples. Use module-qualified calls. The bridge may explicitly extend `AeroBeams.solve!` for its own `CoupledDynamicProblem` type; do not replace existing AeroBeams methods or add methods on unrelated foreign types. Julia requires an explicit import or a qualified function name when extending another module's function. [Julia namespace rules](https://docs.julialang.org/en/v1/manual/modules/#Namespace-management).

**3. Proposed terminology map**

This table is the public migration contract. TOML/CSV spellings remain versioned independently from Julia field names.

| Current name or concept | Proposed public name / location | Semantic requirement |
|---|---|---|
| `System` | Transitional backend `System`; public `UVLMModel`, `UVLMState`, `UVLMWorkspace` | A model is not a mutable simulation workspace; do not alias all three to the old object |
| `initialize_bohnisch_uvlm_system` | `create_UVLMModel` + problem initialization | Bohnisch/Chang defaults move into named example builders |
| `ChangModel`, `build_chang_model` | Benchmark `create_ChangModel`; composed generic UVLM model | Preserve the original linear model as a regression backend |
| `run_chang` | Compatibility wrapper over construction, `solve!`, and explicit output | Preserve returned data through an adapter during transition |
| `steady_analysis!`, `unsteady_analysis!` | `solve!(::UVLMSteadyProblem)`, `solve!(::UVLMDynamicProblem)` | Separate analysis orchestration from one-step numerical kernels |
| `PartitionedCouplingOptions` | `PartitionedCoupling`, `create_PartitionedCoupling` | Keep structural Newton settings separate from outer coupling settings |
| `maximum_iterations`, `state_tolerance`, `load_tolerance` | `maximumIterations`, `stateTolerance`, `loadTolerance` | Preserve residual definitions, normalization, and units |
| `equilibrium_tolerance`, `coupled_equilibrium_tolerance` | `equilibriumTolerance`, `coupledEquilibriumTolerance` where applicable | A tolerance on the Chang matrix residual is not automatically an AeroBeams residual tolerance |
| `speed_mps`, `freestream_speed_mps` | `airspeed` in `OperatingPoint` | Scalar speed in m/s; keep the velocity vector/basis explicit separately |
| `rpm`, `rotation_rpm` | Input `rotationRPM`; resolved `angularSpeed` | Convert once using `angularSpeed = 2π*rotationRPM/60` in rad/s |
| Aeroelastic proportional RPM override | `rotorSpeedPolicy = :constantAdvanceRatio` | Preserve the original reference point; reject contradictory prescribed-speed inputs |
| `angle_of_attack_deg`, `sideslip_deg` | Resolved `angleOfAttack`, `sideslip` in radians | Legacy degree-valued input is converted, not relabeled |
| `azimuth_deg`, `azimuth_step_deg` | Resolved `azimuthStep` in radians | Derive Δt once from prescribed angular speed; later adaptive dt uses actual rotor phase increments |
| `spanwise_panels`, `ns` | `nSpanwisePanels` | Independent of AeroBeams `nElements` |
| `radial_panels` | `nRadialPanels` | Blade discretization, not wing span discretization |
| `chordwise_panels`, `nc` | `nChordwisePanels` | Preserve array order `(coordinate, chord, span/radius)` |
| `nw`, maximum wake rows | `maximumWakeRows` | Allocated capacity per surface |
| `nwake`, `iwake`, active wake rows | `activeWakeRows` in accepted/trial state | One counter owner; valid range `0:maximumWakeRows` |
| `wake_revolutions` | `wakeRevolutions` in a wake-retention policy | Conversion to rows depends on time/azimuth policy |
| `core_radius_m`, core factors | `coreRadius`, `segmentCoreFactor`, `chordCoreFactor` | Preserve distinction between fixed physical radius and a discretization-dependent rule |
| `interaction_on`, `surface_id` | `interactionGroups` / resolved interaction policy | `false` still retains within-propeller blade interaction |
| `attachment_eta`, nearest-node utilities | `attachments` with beam ID and normalized span position | Use interpolation for off-node attachments; preserve physical position when remeshing |
| `elastic_axis_fraction` | `normSparPos` only after confirming matching geometry convention | Both must refer to the same leading-edge chord fraction; retain an explicit offset otherwise |
| `surface_history` | Optional `surfaceGeometryOverTime` | Heavy geometry history is separate from response samples |
| Response `time`, output arrays | `savedTimeVector`, `responseTimeVector`, `responses` | Each history has an explicit sampling policy; no assumed identical lengths |
| `moving_block_lambda_per_s` | `growthRate` / serialized `growth_rate_per_s` | λ in s⁻¹; negative means decay |
| `frequency_hz`, `omega_radps` | `frequencyHz`, `angularFrequency` | AeroBeams eigen `frequencies` are angular frequencies in the local implementation |
| `damping_ratio_percent` | `dampingRatio` internally, explicit percent output | ζ = −λ/√(λ²+ω²); fraction and percent are different quantities |
| Current startup `trim_*` | `baselineRevolutions`, `baselineAverageRevolutions` | Current mean-load subtraction is baseline preparation, not a solved deformed trim |
| `validation_passed` | Separate `completed`, `converged`, `acceptedForStudy` | Numerical completion, solver convergence, and scientific acceptance are different outcomes |

Do not reinterpret backend `PanelProperties.gamma`, `Γ`, `dΓ`, and `dΓdt` solely from their names. Audit storage and normalization along circulation and force paths, document each unit, and keep conversions in one adapter. In particular, `dΓ` contains freestream derivative arrays and must not be silently renamed as a temporal circulation increment.

**4. Target packages, files, and responsibility**

Proposed layout; move files incrementally rather than creating empty abstractions for every future feature:

```text
lib/WingPropellerUVLM/
  src/WingPropellerUVLM.jl
  src/OperatingPoint.jl
  src/UVLMSurface.jl
  src/UVLMModel.jl
  src/UVLMSolver.jl
  src/UVLMState.jl
  src/UVLMWorkspace.jl
  src/UVLMProblem.jl
  src/UVLMResults.jl
  src/Core.jl
  src/Kinematics.jl
  src/Loads.jl
  src/Compatibility.jl
  src/backend/                       existing numerical kernels initially retained
  src/aeroelastic/                   existing linear benchmark solver retained
  examples/chang_linear_aeroelastic/ case data, builders, thin executable entries
  studies/Project.toml               reproducible analysis environment
  studies/src/Convergence.jl
  studies/src/MovingBlockDamping.jl
  studies/src/Provenance.jl
  test/                              package numerical and compatibility tests

lib/AeroBeamsUVLM/                    proposed independent bridge package
  Project.toml
  src/AeroBeamsUVLM.jl
  src/CouplingFrames.jl
  src/BeamUVLMCoupling.jl
  src/StructuralKinematics.jl
  src/LoadTransfer.jl
  src/PartitionedCoupling.jl
  src/CoupledDynamicProblem.jl
  src/SpinningRotor.jl               when prescribed-spin structural effects validated
  test/
  examples/

src/ExternalResultants.jl            proposed generic AeroBeams facility
src/Problem.jl                      generic structural step transaction hooks
src/SystemSolver.jl                 generic external-load assembly hook
docs/src/                           public API, migration guide, coupled examples
```

`UVLMModel` owns read-only topology, reference geometry, surface/rotor definitions, attachment-independent aerodynamic meshes, and coefficient reference definitions. Resolved operating-point data belongs to the problem when it changes between analyses; the model may contain a default but must not be silently mutated by a speed override.

`UVLMSolver` owns force-model selection, core policy, interaction policy, wake method, and applicable linear-solve/cache policies. Keep interpolation/grid discretization inputs with geometry, rather than scattering the same control across model and solver. Every setting must have one authoritative owner after resolution.

`UVLMState` owns bound/wake circulation, wake geometry, accepted surface history needed for motion/rate estimates, active wake counts, time-dependent flow state, and rotor phase where the aerodynamic subsystem owns prescribed motion. `UVLMWorkspace` owns AIC/RHS storage, factorization caches, temporary grids, velocities, force arrays, and transfer buffers. Results hold saved copies or documented immutable records, never accidental live aliases to workspaces.

The existing `System` initially remains a private backend facade referencing the new storage. Establish alias and rollback tests before physically splitting its fields; do not maintain two competing circulation arrays. Hide the facade only after all active callers use the new lifecycle. Preserve legacy advanced backend exports through compatibility wrappers until a documented breaking release.

The bridge owns structural/aerodynamic IDs, interpolation stencils, frame transforms, dimensional load buffers, coupling iteration records, and a reference to the real structural problem. AeroBeams owns structural equations, constraints, initial-state initialization, Newton iterations, and structural integration. A model cannot be shared concurrently by two mutable AeroBeams problems without an explicit isolation strategy.

**5. Public workflow to implement**

The following is illustrative **proposed API**, not code to execute now. Constructor details are finalized in phase P2 and documented together.

```julia
import WingPropellerUVLM as UVLM
import AeroBeams
import AeroBeamsUVLM as Bridge

aeroModel = UVLM.create_UVLMModel(; surfaces, rotors, reference)
aeroSolver = UVLM.create_UVLMSolver(; coreRadius, interactionGroups)
operatingPoint = UVLM.create_OperatingPoint(;
    airspeed = 85.0,
    referenceAirspeed = 65.0,
    referenceRotationRPM = referenceRPM,
    rotorSpeedPolicy = :constantAdvanceRatio,
)

structuralModel = AeroBeams.create_Model(; beams, BCs)
structuralProblem = AeroBeams.create_DynamicProblem(;
    model = structuralModel,
    systemSolver = AeroBeams.create_NewtonRaphson(),
    Δt, finalTime = 6.0,
    trackingTimeSteps = true, trackingFrequency = 1,
)
mapping = Bridge.create_BeamUVLMCoupling(;
    structuralModel, aerodynamicModel = aeroModel, attachments, frames,
)
problem = Bridge.create_CoupledDynamicProblem(;
    structuralProblem, aerodynamicModel = aeroModel, aerodynamicSolver = aeroSolver,
    operatingPoint, coupling = mapping,
    couplingSolver = Bridge.create_PartitionedCoupling(;
        maximumIterations = 10, stateTolerance = 1e-5, loadTolerance = 1e-2,
    ),
)
AeroBeams.solve!(problem)  # bridge-defined dispatch on its own wrapper type
results = problem.results
```

For standalone aerodynamics, use `UVLM.create_UVLMDynamicProblem(; model, aeroSolver, operatingPoint, Δt, finalTime, ...)` followed by `UVLM.solve!`. Coupled construction allocates UVLM state/workspace directly rather than attaching a second independently advancing aerodynamic problem.

`solve!` mutates the problem and returns a documented value consistently across the new UVLM/bridge APIs; prefer returning the problem for convenience. Do not change AeroBeams' existing return contract globally. Constructors must not run simulations, activate environments, write output, or load plotting backends. A separate writer persists results; legacy `run_chang` may retain its orchestration behavior through a wrapper.

**6. Kinematic, unit, and load contract**

Before coupling, produce a small convention document and executable tests for these invariants:

| Quantity | Required contract |
|---|---|
| Frames | Define `R_UA` by `v_U = R_UA*v_A`; test orthogonality and determinant +1. Define reference-frame translation and time dependence separately |
| Position | Include initial beam geometry/orientation, displacement, section rotation, chord offset, hub/pylon offset, and rotor phase in a single map |
| Rotation | Reconstruct AeroBeams Wiener–Milenkovic deformation rotations using existing utilities; initial `rotationParametrization` and deformation `p` are distinct |
| Velocity | Use physical translational/angular velocities and frame transport; include a moving basis' angular velocity and origin velocity once |
| Wake coordinates | Choose a documented physical frame. If using a moving frame, transform old wake positions and convection consistently; a time-varying matrix alone is insufficient |
| Loads | Structural interface carries dimensional N and N·m, paired with the positions at which the UVLM force formulation applies them |
| Moments | Preserve the moment reference point: translated moment is `m_new = m_old + (r_old-r_new) × f` |
| Structural scaling | Apply AeroBeams `forceScaling` exactly once at its residual boundary; never pre-scale loads and scale them again |
| Spin | Define positive shaft axis, phase, aerodynamic torque, and support reaction by one angular-momentum/power convention |
| Units | SI internally; convert degree/RPM/imperial benchmark inputs once. Plotting labels do not convert solver arrays |

For a surface motion map `x_a = g(q_s,t)`, verify `δx_a = H δq_s` and `Q_s = Hᵀ f_a`, hence `f_aᵀδx_a = Q_sᵀδq_s`. For rotational coordinates, a generalized force conjugate to `p` is not automatically a physical moment. Derive the appropriate tangent map before inserting a wrench into AeroBeams equilibrium equations. Use actual node/element reconstruction and IDs; never slice `DynamicProblem.x` in Chang's DOF order.

Do not equate structural and aerodynamic meshes. Precompute span/element attachment stencils, then update only deformation-dependent geometry/tangents each trial. Test nonmatching meshes, nonzero rotations, chord offsets, arbitrary rigid frame rotations, and off-node propeller attachments. Rebuild maps after topology/equation numbering changes such as `update_model!`.

**7. Solver integration and transaction contract**

Start from the existing `snapshot_uvlm`, `restore_uvlm!`, `advance_uvlm_trial!`, `advance_wake!`, and `commit_wake_rows!` behavior. The new API must make a noncommitting trial explicit; the existing trial helper defaults `advanceWake=true`, so callers cannot rely on its default for partitioned coupling.

Proposed operations and their mutation boundaries:

| Operation | Allowed mutation | Required postcondition |
|---|---|---|
| `initialize!` | Fresh problem state/workspace | Initial state, wake, phase, and time agree; no completed physical step |
| `begin_time_step!` | Step snapshot and step context | Accepted state captured once; current Δt and trial time fixed |
| `evaluate_trial!` | Trial state and scratch only | Repeat evaluation from identical inputs produces identical loads; accepted wake unchanged |
| `commit_time_step!` | Accepted state, wake, phase, completion status | One commit token per step; duplicate commit rejected |
| `rollback_time_step!` | Restore both disciplines and step bookkeeping | No rejected history entry, wake row, rotor phase increment, or advanced reference frame survives |
| `save_time_step!` | Selected result records only | Store accepted state with its actual time and sampling policy |

One coupled physical step proceeds as follows:

1. Save accepted structural/aerodynamic state and current external loads. AeroBeams defines the trial end time and Δt.
2. Prepare basis motion, boundary conditions, and equivalent state rates once for this attempted step. Preserve the beginning-of-step integration data across outer coupling iterations.
3. Reconstruct all coupled aerodynamic surfaces from the latest structural iterate. Restore accepted aerodynamic history before the aerodynamic trial, then apply that trial geometry using the established backend order.
4. Solve circulation and dimensional loads globally with wake advancement disabled. Transfer loads to a fresh external-resultant buffer.
5. Hold that buffer fixed in basis A during the inner AeroBeams structural solve. Structural Newton convergence is necessary but does not establish outer coupling convergence.
6. Recompute aerodynamic loads at the resulting structural iterate. Check state, load, and structural equilibrium with these updated loads. An accepted step must represent a consistent structural-state/load pair; avoid accepting solely because the lagged-load Newton solve converged.
7. Repeat with configured relaxation until all required criteria pass. On success commit aerodynamic history, wake and rotor phase exactly once, then save results. On failure restore both disciplines and report a structured failure; adaptive retry is enabled only after its separate gate.

For AeroBeams, add a generic `ExternalNodalResultants` facility and append its residual contribution after element/special-node assembly in `assemble_system_arrays!`. During the first partitioned implementation the loads are frozen through the inner solve; therefore the external UVLM contribution has no inner aerodynamic tangent. Keep existing gravity, structural loads, and any explicitly separate sectional aerodynamics. Reject overlapping sectional and UVLM coverage by default to prevent duplicate aerodynamic loads.

The earlier blueprint locates the external-load buffer on `Model`; retain that minimal integration point initially, document its single-run mutable ownership, and include it in the bridge snapshot. Do not introduce a module-global registry keyed by model or problem. The bridge rebuilds the equation map when topology changes. First support unconstrained translational/rotational equilibrium locations and tested rigid connections; reject unsupported hinge-side and prescribed-DOF cases explicitly until their resultant/reaction assembly is verified.

Expose generic structural step preparation/solve/accept/restore operations in `Problem.jl`, factoring existing code rather than duplicating `time_march!`. Inspect both `time_march!` and `adaptable_time_march!`. Current `copy_state`/`restore_state!` covers only a subset of state: inventory time/index values, basis orientation, rates and equivalent rates, constraints, BC caches, solver status, complementary aerodynamic variables, external loads, and output history. Save only after acceptance. Initialization and initial-rate correction must also receive the correct external loads and must not advance a physical step.

The current `solve!(::Problem)` branches on exact types. A new subtype alone does not obtain a solver. The proposed bridge wrapper contains a real `DynamicProblem`, and its explicit `AeroBeams.solve!(::CoupledDynamicProblem)` method orchestrates the generic step API.

For later monolithic coupling, `assemble_system_arrays!` already gathers element states before evaluating element contributions. A global aerodynamic evaluation would additionally require up-to-date rotations and rates for every coupled element, then nonlocal residual/Jacobian blocks across surfaces. Do not solve the global lattice inside every element's `aerodynamic_loads!`. Start derivative verification with finite differences at fixed accepted history; the current force derivative limitations and Float64-specific paths prevent an assumption of automatic AD support.

**8. Rotor structure and analysis semantics**

Keep prescribed speed with the requested proportional scaling as the first supported rotor case. Preserve the verified aerodynamic reference operating point in metadata when the aeroelastic speed changes. For multiple rotors, resolve speed/phase/rotation direction per rotor. Do not infer every rotor's phase from one global RPM.

Moving from the Chang linear model to AeroBeams changes the structural formulation. Map sectional stiffness/inertia distributions into beams, and compare structural-only mass, static compliance, modes, and attachment motion before comparing coupled responses. Case-fitted frequencies are calibration targets, not direct entries in an AeroBeams global stiffness matrix.

Rotor mass and nonspinning inertia may be represented using native inertial components. Relative shaft-spin angular momentum and the resulting support reaction still require explicit treatment; do not assume `PointInertia` automatically reproduces Chang's gyroscopic matrix. Check torque sign against angular-momentum balance and power, and avoid adding inertial contributions twice. Validate the zero-spin limit and the small-angle gyroscopic matrix against the retained linear benchmark before using the full coupled rotor case.

The local README states structural damping is not modeled by AeroBeams. A nonzero Chang pylon or stiffness damping setting must therefore be rejected with an actionable message until a documented physical damping component is implemented, or supported through an explicitly labeled legacy model. Numerical integrator damping cannot stand in for physical mount damping. Zero physical damping is the first parity case.

Do not label startup mean-load subtraction as `TrimProblem`. Introduce baseline preparation as its own operation. Later steady/trim support must state whether the rotor operating state is steady, azimuth-frozen, phase-averaged, or periodic. A rotor-periodic base state may need periodic/Floquet treatment rather than an ordinary autonomous eigenproblem; that is a separate numerical design task.

**9. Damping, convergence, and output migration**

Move the current moving-block module into the studies layer once provenance migration is ready. Remove `MovingBlockStudy` and the entry-local method extension only after the central study runner accepts an explicit estimator object, proposed `MovingBlockDampingMethod`, and every public entry reaches that same runner.

Preserve the reference settings and calculations: peak-based cut, 25–50% record sizing, one-sample shift, per-block mean removal, rectangular FFT, global dominant component, and linear log-amplitude slope. Preserve the optional fixed-duration mode, frequency-band acceptance check, minimum peaks/R², and CSV/PNG diagnostics. Rename controls through the input adapter; keep old method fingerprints and output fields readable.

Record `growthRate`, `frequencyHz`, `angularFrequency`, `frequencyResolutionHz`, `dampingRatio`, block duration/count/size, fit interval, R², and explicit validity reasons. Do not call FFT-bin spacing an uncertainty bound. Include off-bin signals, out-of-band signals, mode switching, offsets, noise, and fitting-window sensitivity in acceptance studies. Offset sensitivity retained to match the user's reference must remain documented; a detrended variant requires a separate method/version and comparison.

For fixed-step damping studies, keep dense response samples even when heavy geometry tracking is disabled or decimated. An adaptive accepted time grid is generally not uniform: either restrict moving-block analysis to fixed-step data initially, or implement an explicit validated resampling stage that records grid, interpolation, and effects on damping. Do not pass irregular samples directly to the FFT or silently fit only saved animation frames.

Keep aerodynamic and aeroelastic convergence independent. Add structural mesh convergence after decoupling meshes, wing force/moment metrics for wing aerodynamic convergence, and a final combined aeroelastic refinement check. CT/CQ alone cannot certify all wing loads, particularly when interaction is disabled. Separate output fields for `completed`, `solverConverged`, `metricValid`, `familyAccepted`, and `combinedAccepted`.

Define a proposed version-2 result schema with explicit units/bases, requested and actual operating point, requested/accepted time grid, termination reason, residual history, discretization, solver/method IDs, input hashes, dependency versions, and result hashes. Keep legacy snake_case serialized keys unless a versioned conversion is necessary; adopting camelCase Julia fields does not require breaking CSV consumers.

Existing source fingerprints include convergence, Chang, and package source plus the manifest. Moving numerical files will invalidate those fingerprints. Preserve old selections and their seals unchanged. Provide two explicit paths: read them for archival/postprocessing using the old schema, or rerun/verify under the migrated numerical source to issue a new selection. Matching rename-regression tests alone do not authorize rewriting an old verification seal.

New provenance should separate aerodynamic solver, structural model/solver, coupling, estimator, and output schema identities using conservative dependency manifests, not manually selected hashes that omit indirect dependencies. Record the loaded source identity at process start; use a fresh Julia process for provenance-bearing runs after edits. Key any case reuse by all inputs affecting that stage, including speed/RPM, initial state, excitation, time grid, and solver tolerances.

**10. File migration worklist**

Paths are relative to `lib/WingPropellerUVLM` unless otherwise specified.

| Existing files | Work to perform | Boundary to preserve |
|---|---|---|
| `src/WingPropellerUVLM.jl` | Introduce public includes/exports and compatibility layer incrementally | Package UUID/name and existing imports remain stable |
| `src/backend/system.jl`, `src/UVLMState.jl` | Inventory fields; build facade; separate state/workspace; formalize transaction | Array ownership, accepted history, force/circulation conventions |
| `src/backend/{panel,geometry,circulation,induced,analyses,wake}.jl` | Add adapters first; rename/move internal files only after numerical regression | Core kernels and wake chronology |
| `src/backend/{nearfield,nearfield_postprocessing,reference,freestream}.jl` | Document units; expose in-place dimensional load path | Dimensional loads versus coefficients; vortex locations and moments |
| `src/backend/{legacy_nearfield,nonlinear,rotors,vspgeom,farfield,stability,visualization}.jl` | Inventory legacy exports and mark support status; retain behind adapters | Do not silently remove less-used workflows or claim new validation |
| `src/wing_propeller/{Initialization,GridUtilities,Kinematics,BladeGeometry}.jl` | Separate generic geometry/motion from named tabulated cases | Attachment interpolation, grid ordering, blade twist and spin direction |
| `examples/.../src/{ChangAeroelastic,chang_configuration,chang_model_parameters}.jl` | Thin case builder + input adapter; generic settings become typed constructors | Explicit config isolation and proportional RPM behavior |
| `examples/.../src/{chang_structural_matrices,chang_structural_model}.jl` | Retain benchmark implementation; add separate AeroBeams builder later | No wholesale replacement of M/C/K by superficially similar beam inputs |
| `examples/.../src/{chang_workspaces,chang_uvlm_coupling,chang_simulation}.jl` | Extract workspaces and generic loop; retain Chang kinematic adapter | Virtual work, trial repeatability, accepted state/load consistency |
| `src/aeroelastic/{GeneralizedAlpha,Excitations}.jl` | Preserve linear benchmark API; separate excitation from solver | Generalized-alpha is not used to advance the AeroBeams structural state |
| `examples/.../src/chang_postprocessing.jl` and plotting/animation files | Structured status/results and explicit writers | Preserve response columns and optional visualization |
| Convergence entries, `src/Convergence.jl`, stage adapters, selection/options files | Consolidate runner; move reusable study orchestration | Selection verification, resolved input keys, case deduplication |
| Convergence `MovingBlockDamping.jl`, `moving_block_adapter.jl`, legacy metrics | Promote reference method; remove temporary wrapper after API replacement | User-selected algorithm and old diagnostic readers |
| Root `src/{AeroBeams,Model,Problem,SystemSolver}.jl` | Add generic external resultants and structural step API | Existing uncoupled analyses unchanged when no external loads are attached |
| Root `src/Core.jl` | Audit residual/load hooks and rotation reconstruction; add only necessary generic support | Constrained equations, units and force scaling |
| Both package projects, study/bridge projects, tests and docs | Explicit dependency graph; isolated tests; update executable examples | No fallback to undeclared packages in the default Julia environment |

**11. Implementation phases and acceptance gates**

Implement each phase in reviewable commits. P0–P5 deliver the requested public terminology/structure for standalone UVLM and its current Chang workflow. P6–P9 establish actual AeroBeams coupling. P10–P11 cover performance and broader analysis capabilities. The latter do not block accepting the naming/organization migration.

| Phase | Concrete deliverables | Prerequisites | Gate before proceeding |
|---|---|---|---|
| P0 — Freeze baseline | API/export inventory; resolved preset snapshots; archived small numerical fixtures; dependency audit; supported-feature matrix | None | Fresh-process tests and current smoke produce recorded results; exact current limitations documented |
| P1 — Environments and entry points | Declared study/test dependencies; headless tests; guarded script entry points; one public stage runner | P0 | `LOAD_PATH=["@","@stdlib"]` works in designated environments; include performs no solve/output; both speed entry paths agree |
| P2 — Naming and constructor facade | Types and terminology map; constructor validation; old-key conversion/deprecation; thin builders | P1 | Old/new inputs resolve identically; contradictory aliases error; old numerical fixtures match |
| P3 — State/problem/workspace separation | `UVLMModel`, problem types, state/workspace, facade ownership, explicit step operations | P2 | Repeat-trial, rollback, single-commit, array-alias, and independent-run isolation tests pass |
| P4 — Results and study layer | Named status/results; output/history policy; explicit estimator; central convergence runner | P3 | Current moving-block numerical/diagnostic fixtures preserved; incomplete/unconverged runs never accepted |
| P5 — Provenance and documentation | Versioned schemas; conservative stage identities; legacy readers; guide and example migration | P4 | Old selections remain unchanged/readable; stale caches rejected; new selections verify against new source/environment |
| P6 — Generic AeroBeams interfaces | External resultants; structural step prepare/solve/commit/restore; default no-load behavior | P0; P3 contract | Structural test suite passes; external load signs/scaling/reactions verified; forced retry restores all mutated fields |
| P7 — Geometry/load bridge | Separate bridge package; explicit frames; node/element maps; nonmatching meshes | P5, P6 | Rigid-motion, force/moment, virtual-work and moving-frame tests pass on declared supported connections |
| P8 — Coupled fixed-step dynamics | Wrapper problem; frozen-load inner solve; outer residual checks; initialization and failure records | P7 | Prescribed rigid motion, structural-only limit, repeat trials and one wake commit pass; a small flexible-wing benchmark converges |
| P9 — Chang rotor parity | Separate AeroBeams Chang builder; prescribed-spin inertia/reaction support; mode/response comparisons | P8 | Structural-only match established first; zero-spin and small-angle gyroscopic checks; coupled time/mesh/coupling convergence |
| P10 — Measured performance | History reductions, buffer reuse, cache rules, group solves, benchmark report | P3–P9 as relevant | Equivalent loads/states within fixture tolerances; repeatable warm-run allocation/time evidence |
| P11 — Advanced analyses | Adaptive coupling and resampling; steady/trim; nonlocal linearization; broader hinges/rotor models | P8/P9 and individual design review | Dedicated numerical validation for each feature; unsupported modes rejected until then |

The critical sequence is P0 → P1 → P2 → P3 → P4 → P5 → P7 → P8 → P9; P6 supplies the structural interface required by P7. Performance work may begin after its ownership dependencies are stable. Avoid a calendar estimate until P0 measures current environments, fixture runtime, and required structural interfaces.

**12. Validation matrix and concrete completion criteria**

| Validation group | Required cases / assertions |
|---|---|
| API equivalence | Current Chang preset and both convergence entries; old/new keyword forms; unknown/duplicate keys; defaults; speed override 65→85 m/s with RPM ×85/65 |
| Environments | UVLM alone; studies with FFTW; plotting enabled/disabled; bridge with local packages; no global environment fallback; Julia supported-version CI |
| Aerodynamic preservation | Existing panel/core/force tests; rigid wing, isolated propeller, interacting wing/propeller; Γ, dΓdt, dimensional loads, wake geometry, CT/CQ histories |
| Transactions | Two identical trials; A→B→A trial sequence; rejected coupling step; attempted double commit; active/capacity wake limits; phase rollback; snapshots unaffected by later mutation |
| Structural hook | Existing AeroBeams suite with no provider; known static force/moment; dynamic prescribed loads; reactions; special nodes; explicit rejection of unsupported connections |
| Transfers | Independent meshes, offsets, rigid rotations, nonzero Wiener–Milenkovic rotations, moving frame; force resultant, moment resultant, virtual work |
| Coupled limits | Zero aerodynamic loads; prescribed rigid structure; wing without propeller; no cross-group aerodynamic interaction; zero-spin rotor; small-angle coupled comparison |
| Completion | Requested discrete grid/sample count versus accepted grid; integration failure; outer-coupling failure; finite-state/load limits; no invalid damping acceptance |
| Damping | Current reference regression; decay/growth; off-bin/out-of-band signals; fixed-duration mode; offsets/noise/multiple modes; irregular-time rejection; CSV/TOML/plot consistency |
| Studies/provenance | Original selection unmodified; overridden point correctly nested; combined verification; stage-specific invalidation; changed excitation/solver/estimator prevents improper reuse |
| Performance | Warm numerical sections separately from compilation and I/O; allocated bytes; history memory growth; per-case isolation |

Use the current [smoke](../reviews/uvlm_aerobeams_2026-09-14/coupled_smoke.jl) and [summary](../reviews/uvlm_aerobeams_2026-09-14/smoke_summary.toml) as the starting Chang regression evidence. Its 38-step runs do not validate a six-second damping result. Preserve its virtual-work threshold of `1e-6` for that fixture; extend test coverage before applying the same threshold to other conditioning/scales.

For rename/facade phases, aim for identical resolved inputs and deterministic outputs. If floating-point evaluation order changes, define scaled absolute/relative tolerances from fixture magnitudes and demonstrate that deviations are far below study acceptance thresholds. Do not loosen tests simply to pass a new result. AeroBeams replacement is a formulation change: compare converged physical observables and documented modeling assumptions rather than requiring bitwise agreement with the Chang integrator.

Each completed phase must include migrated call sites, appropriate tests, updated API documentation, and a short record of remaining unsupported features. No phase is complete merely because a constructor exists or a single example starts running.

**13. Performance work tied to the new structure**

Start with optional geometry tracking and workspace ownership. Keep dense lightweight response histories for damping, while making full panels/wakes opt-in with a stride and memory estimate. The prior review estimated roughly 1.18 GiB of panel payload for one representative six-second case; confirm actual configurations rather than treating this as a universal measurement.

Then reuse force, geometry, snapshot, and interpolation buffers; cache factorizations only under explicit invalidation rules. The retained Chang generalized-alpha effective matrix may be cached for fixed M/C/K/Δt and parameters; its gyroscopic terms mean an SPD factorization cannot be assumed. AeroBeams' structural Jacobian follows its own update rules.

Cache aerodynamic work by geometry, interaction policy, core policy, and relevant flow/time dependencies, not panel count alone. Exploit disconnected AIC interaction groups only after verifying identical within-group coupling. Test motion flags so `calculateInfluenceMatrix=false` cannot falsely promise reuse while geometry is changing.

Measure geometry update, AIC assembly/factorization, force evaluation, wake convection, transfer, structural solve, snapshot/restore, estimator, and I/O separately. Reuse the moving-block FFT plan and buffers. Consider process-based case concurrency only with measured memory limits and isolated workspaces/output directories. Wake coarsening, FMM, subcycling, and new integration methods are later numerical changes, not part of terminology refactoring.

**14. Compatibility and retirement policy**

Accept old snake_case Julia keywords for a documented transition release. Translate once at construction, emit one actionable deprecation warning at the user boundary, and retain old names in archive readers. Reject supplying both old and new aliases unless the input adapter deliberately verifies equivalent values. Never emit warnings from panel loops or every time step.

Keep `WingPropellerUVLM` package identity and current example paths. Thin scripts call the new API and preserve familiar editable settings during transition. Remove the entry-local moving-block adapter only when the central runner implements equivalent behavior; do not leave two estimators selected by namespace accidentally.

Keep compatibility wrappers until the migrated examples, tests, docs, and legacy-reader fixtures all pass and a breaking-version policy has been stated. Removal is a final dedicated change, not part of moving files. Preserve archived source/result hashes and avoid reusing output directories with incompatible schemas.

The first implementation increment should deliver P0–P1: establish clean environments, capture current outputs, guard entry points, and consolidate runner semantics. Then implement the P2 constructor facade using the existing backend. This gives the user an AeroBeams-style workflow early while preserving a verified route back to the current UVLM calculations.
