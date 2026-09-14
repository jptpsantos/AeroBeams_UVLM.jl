**Implementation review: Chang convergence, WingPropellerUVLM, and AeroBeams integration — 14 September 2026**

**Follow-up:** The [full terminology and structure migration plan](../../plans/UVLM_AEROBEAMS_MIGRATION.md) defines the implementation sequence and acceptance gates. This review records the earlier snapshot: moving-block postprocessing and the aeroelastic entry's nested speed-override delegation were subsequently implemented, as noted in the plan.

The reviewed time-marching and load-transfer path has a sound foundation, but the current automated damping acceptance can report a misleading result. Correct the damping estimator and project isolation before using `accepted=true` as the basis for final mesh or flutter conclusions. The speed override works through the new entry-level API; the older internal API still silently uses the original speed.

This review concerns the current working tree, including uncommitted changes, at base commit `68b685be56e0f8800bcb6c632aca6ef7e6a8b5db`. The numerical source fingerprint exercised by the smoke tests is recorded in [smoke_summary.toml](smoke_summary.toml). Production code and existing aerodynamic result files were not modified during this review. New files in this directory contain diagnostics, test outputs, and this proposal.

The review follows the active convergence-to-solver call chain in detail. It also examines the local AeroBeams extension points and inventories adjacent code. It is not a proof of all exported UVLM methods, all parameter combinations, or agreement with experiments. In particular, the legacy nonlinear/rotor-import/stability utilities and plotting backends did not receive the same dynamic coverage as the active Chang solver.

**Evidence obtained**

All existing package assertions passed under Julia 1.11.5:

| Existing test group | Passed |
|---|---:|
| Aerodynamic convergence integrity | 41 |
| Chang fixed-core configuration and propagation | 73 |
| Fixed-core convergence study | 66 |
| Independent two-stage convergence | 35 |
| Independent aeroelastic configuration | 24 |
| WingPropellerUVLM core suite | 130 |
| Total | 369 |

The added [coupled_smoke.jl](coupled_smoke.jl) exercises 85 m/s with RPM scaled from 1212 RPM at 65 m/s. It runs 38 physical steps through baseline collection, excitation, and genuine structural/aerodynamic iterations. Both `interaction=false` and `interaction=true` completed their requested discrete time grids; the former needed at most three iterations and the latter four. These deliberately short, coarse cases demonstrate execution and state consistency, not converged aerodynamics or a meaningful damping estimate. The additional 13 assertions passed.

A frozen-force virtual-work test checked all 26 free DOFs with nonuniform wing rotations and nonzero propeller angles. Its maximum scaled discrepancy was `1.5162671857330032e-7`, below the `1e-6` threshold. This supports the current work-conjugate transfer, including the interpolated attachment, for the exercised geometry.

[reproduce_findings.jl](reproduce_findings.jl) and [diagnostics.toml](diagnostics.toml) contain independent synthetic counterexamples for the damping estimator. They intentionally record the current problematic behavior rather than treating it as a passing scientific requirement.

The sandbox required a temporary first Julia depot for Pkg's usage log. Julia emitted Windows temporary-pidfile cleanup warnings. PowerShell also wrapped stderr activation messages as `NativeCommandError` in [coupled_smoke.log](coupled_smoke.log); its test summaries and final TOML were nevertheless produced with all 13 assertions passing. The existing package-suite invocation completed with exit code zero.

**Correctness findings and fixes to prioritize**

1. **High priority: a mode outside the requested band can pass as an in-band mode.** In `studies/convergence/src/damping_metrics.jl`, lines 187–216 select the largest bin *inside* the band, without establishing that an actual mode exists there. Line 256 accepts based on the amplitude fit's R². A synthetic `exp(-0.15t) sin(2π·8t)` signal, analyzed over 3–6 Hz with the default one-second Hann window, is reported as approximately **5 Hz**, with lambda `-0.15008175 s⁻¹`, R² `0.9977275`, and `valid_for_convergence=true`. Leakage has a clean exponential envelope too. A good envelope fit alone does not establish modal identity. Add peak prominence/energy checks against the full usable spectrum, detect band-edge/leakage cases, track a consistent modal component, and test known out-of-band signals. Where modes overlap, a validated multimode fit is preferable to taking each block's largest bin.

2. **High priority: frequency quantization can conceal differences larger than the acceptance tolerance.** The same code stores the bin center as the measured frequency, then takes its median. For a one-second block, the bin spacing is approximately 1 Hz, yet stage 2 requests 0.1 Hz agreement. The synthetic 4.1 Hz and 4.4 Hz signals both report approximately **4.0 Hz** and valid lambda values within `0.00002 s⁻¹` of one another. Both channels constructed this way would pass the numerical comparison despite a true 0.3 Hz frequency difference. Use a validated sub-bin estimator or time-domain modal fit, record its uncertainty, and require that uncertainty to be compatible with the convergence tolerance. Increasing duration helps resolution, but simply zero-padding an FFT does not add information or establish accuracy. Test off-bin frequencies and mode switching, not only sinusoids aligned with bins.

3. **High priority for reproducibility: the examples/tests depend on undeclared FFTW and Plots.** `Convergence.jl:7` eagerly includes `damping_metrics.jl`, whose line 2 imports FFTW—even for the aerodynamic stage. `lib/WingPropellerUVLM/Project.toml` declares neither FFTW nor Plots; neither is present in that project's manifest. With `LOAD_PATH=["@","@stdlib"]`, both resolve to `nothing`. The successful suite currently benefits from the user's shared Julia environment. A clean installation can therefore repeat a package-not-found failure; plotting can also fail only after expensive simulations finish. Give the studies a declared, reproducible environment, or declare the required dependencies at the appropriate package boundary. Keep plotting optional and validate it before executing a study that requests plots. Add a clean-project CI run without the default global environment.

4. **Medium priority: two similarly named run functions have different speed semantics.** `chang_convergence_aeroelastic.jl:89` verifies the original selection and applies `speed_mps`; `src/aeroelastic_adapter.jl:108` loads the original selection and ignores that field. Consequently, the previously documented call `ChangConvergenceAeroelastic.Convergence.run_aeroelastic(settings)` can run the original speed while `settings.speed_mps` requests another one. The new `main()` path is correct. Consolidate these entry points into one public API, and have any retained compatibility wrapper delegate to it or reject an unsupported override explicitly. The entry-only implementation preserved existing numerical fingerprints, but that compatibility choice should not become the permanent API design.

5. **Medium priority: offset-dependent peak screening rejects otherwise identical dynamics.** `damping_metrics.jl:126` forms a positive threshold from the raw signal, before the per-block mean subtraction. The valid damped 4.1 Hz sinusoid becomes unavailable after subtracting a constant 2 from the response. Physical damping is unchanged. Detrend/center appropriately before peak qualification, with a documented treatment of equilibrium drift; do not let the sign of a response offset control validity. Include positive and negative offsets, trends, and decays near the noise floor in regression tests.

6. **Medium priority: completion checks are weaker than the solver's own completion data.** `completed_requested_window` accepts a history ending up to `1.01*final_step` before the requested end. For example, `0:0.01:0.99` satisfies a requested end time of 1.0 although the regular requested grid contains 1.0. The normal solver intentionally uses a range that may stop just below the requested time, so exact equality to that requested time is not the fix. Compare against the expected generated grid and its sample count. Also consume structured termination/coupling information instead of inferring success from time alone. `aeroelastic_case` currently discards `run_chang`'s return value and `damping_result` reads only time/pitch/yaw. Validate the actual completion status, all accepted coupling steps, and configured state limits before fitting.

7. **Medium priority: including the aerodynamic entry currently launches the full study.** Its final `PROGRAM_FILE` guard is commented out, so inclusion in a REPL or another script starts `main()` immediately. Tests now explicitly filter out that call. This is a substantial surprise for reusable library code and can write into the default output folder. Restore a guarded command-line entry and expose a separate explicit REPL `main()`/`solve!` call. Do not rely on AST filtering as the long-term import mechanism.

8. **Medium priority for large cases: surface history is retained even when animation is disabled.** `wing_propeller/Initialization.jl:124` requests every step; `chang_simulation.jl:499` copies every surface into that history. This is independent of the separate animation flag and its own history. The smoke test measured 107,352 bytes for only 12 panels and 38 steps. Each `SurfacePanel{Float64}` occupies 208 bytes. As an illustrative calculation, 700 panels, 1212 RPM, 5-degree increments, and six seconds retain approximately **1.18 GiB of panel payload alone** in this history. At 85 m/s with proportional RPM that becomes approximately 1.55 GiB. Add explicit history policy and stride controls; convergence runs usually need response histories and selected snapshots, not every panel at every step.

Additional interface edge cases found by inspection: configuration allows zero allocated wake rows, but `shed_wake!` indexes `wake[end,j]`; either disallow zero for shedding runs or implement its semantics. A response duration smaller than one generated time step leaves `dt` empty and initialization indexes `fs_vec[1]`; report an actionable validation error before allocation. `surface_motion=true` is hard-coded in `propagate_system!`, making `calculate_influence_matrix=false` ineffective for skipping assembly. These are not triggered by the supplied positive-wake, multi-step presets.

**Scientific acceptance limits, separate from implementation defects**

The aerodynamic selector correctly requires complete families, a stable final pair, periodic CT/CQ, agreement of all finer retained levels, and a combined candidate/finer check. It does not simply choose the last or finest setting. The combined check keeps the selected core fixed. Source and result hashes are verified, including the supporting histories.

However, CT/CQ measure propeller loads. With `interaction=false`, changing the wing grid does not test convergence of wing aerodynamics through those metrics. Stage 2 then freezes the wing span count, which also sets the structural mesh. Establish wing-load and structural-mesh convergence separately; treat the saved wing count as a baseline until that evidence exists. Even with interaction enabled, insensitive CT/CQ do not guarantee converged wing moments.

Stage 2 presently reports individual parameter families. It neither publishes a final verified aeroelastic selection nor automatically verifies all chosen refinements together. Add an explicit combined verification and a machine-readable overall outcome after the estimator is corrected. The current generic `report.md` does not itself certify that all families passed.

The startup operation is **mean aerodynamic load subtraction about the undeformed configuration**: structural response is suppressed before the pulse, then the mean is subtracted. It is not AeroBeams' steady deformed trim solution. Periodic rotor loads remaining after mean subtraction may still force the response. Check settling, pulse-amplitude independence, fitting-window sensitivity, and separation from rotor harmonics before interpreting lambda as a free modal decay/growth rate. Consider a phase-resolved baseline or a properly linearized perturbation formulation when justified by the intended analysis.

The generalized-alpha structural method alone does not establish second-order accuracy of the coupled algorithm. Wake positions use an explicit velocity-times-dt update, circulation rates use backward differences, and geometry velocities come from consecutive configurations. Measure the order and damping sensitivity of the complete coupled scheme. Near a stability crossing, a fixed `0.01 s⁻¹` acceptance tolerance can allow a change of sign; include a stability-sign/uncertainty criterion appropriate to that investigation.

**File-by-file assessment of the active path**

Paths below are relative to `lib/WingPropellerUVLM`, except the AeroBeams rows.

| File/module | What was checked | Assessment / next action |
|---|---|---|
| `Project.toml`, `src/WingPropellerUVLM.jl` | Package boundary, imports, exports | Core loads successfully; examples/tests lack declared FFTW/Plots. A number of case-specific and legacy helpers are exported as if generic. |
| `examples/.../chang_convergence_aerodynamic.jl` | Editable inputs, speed/RPM, activation, entry behavior | Speed variable is now correctly scoped. Explicit simulation call on include needs guarding. |
| `examples/.../chang_convergence_aeroelastic.jl` | Original verification, override, RPM propagation, recorded reference | Correct through entry-level `main`; unify duplicate runner semantics. Original seal is correctly nested rather than attached to the modified operating point. |
| `studies/convergence/src/Convergence.jl` | Control validation, case keys, reuse, error comparisons, reporting | Individual acceptance logic is coherent. Damping false positives propagate through it. Add overall outcome, completion evidence, and history policy. |
| `src/aerodynamic_selection.jl` under convergence | Selection, combined check, load-time integrity | Substantial verification is present. Keep checks when reorganizing APIs. This does not verify wing or damping convergence by itself. |
| `src/aerodynamic_options.jl` under convergence | Options, source/result hashes, periodicity, phase comparison | Phase-grid union handles differing sample grids. Fingerprint scope is broader than aerodynamic dependencies and does not prove that an already-loaded Julia module matches edited disk files. |
| `src/aerodynamic_backend.jl` under convergence | Geometry, spin, wake stepping, load normalization | Uses dimensional vertex loads and explicit configuration. AIC/force work and frequent allocation are optimization targets. Revolutions log CL/CT; adding CQ would improve monitoring. |
| `src/aeroelastic_adapter.jl` under convergence | Speed inputs to solver, wake rows, fit timing, case output | Effective physical inputs are consistently consumed by the lower-level path. Original-speed runner remains a trap; return/validate solver termination information. |
| `src/damping_metrics.jl` under convergence | Peaks, FFT, regression, frequency, completion | Highest-priority correction area; three synthetic counterexamples reproduced. |
| `examples/.../src/ChangAeroelastic.jl` | Model/workspace ownership, explicit config, solve/output orchestration | Useful starting separation. Promote reusable types and constructors into the package; keep case data in examples. |
| `chang_configuration.jl` | Parsing, explicit environment isolation, derived defaults, validation | `env=Dict()` correctly isolates convergence cases. Add time-grid/wake edge checks and validate coupling options fully before expensive model work. |
| `chang_model_parameters.jl` | Units, mass distribution, speed/RPM, rotor inertia, attachment, time grid | RPM scales once as intended; model/config values agree in the dynamic smoke. Split tabulated case data from generic model construction. |
| `chang_structural_matrices.jl`, `chang_structural_model.jl` | Mass symmetry, positive definiteness, gyroscopic signs, assembly | Existing tests plus smoke support matrix consistency; zero-damping C is skew-symmetric in the smoke. This linear six-DOF beam benchmark is not a substitute for AeroBeams' formulation. |
| `chang_workspaces.jl` | Per-run storage | Fresh workspaces avoid shared global run state. Several nominal buffers are still followed by fresh allocations inside hot functions. |
| `chang_uvlm_coupling.jl` | Basis mapping, Euler composition, attachment interpolation, virtual work | All-free-DOF frozen-load test passes for exercised case. Preserve this as a regression reference; replace kinematic mapping for AeroBeams rotation parameters. |
| `chang_simulation.jl` | Trial restore, accepted state/load pair, trim, commit, failures | Failed fixed-point steps restore state and throw; accepted trials commit wake once. Add structured termination/checkpointing and remove unconditional surface history. |
| `chang_postprocessing.jl` | Accepted range, residual columns, validation summary | Useful status already exists; convergence should consume it. `validation_passed` denotes computational completion, not physical validation or damping convergence. |
| `src/aeroelastic/GeneralizedAlpha.jl`, `Excitations.jl` | Parameter formulas, Newmark recovery, fixed-point criteria, pulse | Consistent with tested linear model; both load and coupled-equilibrium residuals are checked. Cache the fixed effective operator. |
| `src/UVLMState.jl` | Snapshot/restore and split wake commit | Existing repeat-trial/commit tests pass. Distinguish accepted state, trial state, scratch loads, and active wake count in the future API. |
| `src/wing_propeller/Initialization.jl` | Geometry, interaction IDs, symmetry, active rows, history | Physical attachment positions are interpolated, not snapped. History and empty-time/wake edges need work. Avoid initial geometry aliases in a general reusable workspace. |
| `src/wing_propeller/GridUtilities.jl`, `Kinematics.jl`, `BladeGeometry.jl` | Geometry transforms, interpolation, existing twist tests | Active Chang adapter uses its own six-DOF mapping. Exported `wing_kinematics_from_free_state` assumes older plunge/slope/torsion ordering and must not be used blindly for AeroBeams. |
| `src/backend/panel.jl`, `geometry.jl` | Panel representation, grid mapping, update API | Fixed Float64 entry points and fresh grid/panel construction limit direct AD and allocation efficiency. Geometry changes must invalidate appropriate AIC data. |
| `src/backend/induced.jl`, `circulation.jl` | Core kernel, interaction mask, circulation solve, geometry velocity | Kernel invariance/limits have tests. Disconnected interaction groups form block-diagonal AIC but are solved as one dense system. |
| `src/backend/analyses.jl`, `wake.jl`, `system.jl` | Circulation history, force evaluation, wake commit and storage | Active split-trial architecture is suitable for partitioned coupling. Make motion/cache policy effective and formalize wake-count ownership. |
| `src/backend/nearfield.jl`, `nearfield_postprocessing.jl` | Circulation convention, dimensional loads, node/moment transfer | Existing force/resultant tests and new virtual-work audit support active transfer. Reuse force/node buffers; keep coefficient conventions away from structural load interfaces. |
| `chang_propeller_trim.jl` | Adjacent trim entry and root-search structure | Separate from stage-2 mean-load subtraction; no new windmilling/thrust root-search validation performed here. |
| Plot/animation, legacy nearfield, nonlinear/rotor import, far-field/stability helpers | Dependency/export and call-path inspection | Legacy nearfield has existing tests; other adjacent functionality is not certified by the new smoke. Imperial near-field derivatives are explicitly unsupported. |
| AeroBeams `Model.jl`, `Beam.jl`, `Element.jl`, `Problem.jl` | Ownership, state meanings, constructors, solver flow | Use actual beam states and problem lifecycle, not Chang vector slices. See integration proposal below. |
| AeroBeams `Core.jl`, `Aerodynamics.jl`, `AeroSolver.jl`, `AeroSurface.jl`, `AeroProperties.jl`, `SystemSolver.jl`, `Utilities.jl` | Loads, derivatives, rotations, Newton and time-step boundaries | Requires a global UVLM coupling hook and consistent tangents; an `AeroSolver` subtype alone is insufficient. |

**Efficiency improvements, in implementation order**

| Priority | Proposed change | Why / verification required |
|---|---|---|
| 1 | Disable full surface history for convergence; add `trackingTimeSteps`/`trackingFrequency`-style controls | Measured allocation and analytic memory estimate above. Response CSVs and selected snapshots suffice for the current estimator. |
| 2 | Precompute and factor the generalized-alpha effective matrix for fixed M, C, K, dt and parameters | `GeneralizedAlpha.jl` rebuilds and solves it in every fixed-point iteration. Gyroscopic C makes the effective matrix generally nonsymmetric: use an appropriate factorization such as LU, not an assumed SPD Cholesky. Rebuild on speed, dt or model changes. |
| 3 | Reuse geometry, circulation-sign, segment-force, vertex-position and generalized-load buffers | `near_field_forces!`, `imperial_nodal_forces`, `imperial_nodal_positions`, and the adapter allocate within every trial. Use in-place APIs; verify resultant force, moment and virtual work after changing ownership. |
| 4 | Solve disconnected interaction groups independently | With `interaction=false`, wing and each propeller group are separate AIC blocks. Keep the blades of one propeller together. Compare against the dense reference; do not drop within-propeller interactions. |
| 5 | Introduce explicit AIC invalidation/caching rules | Rigid unchanged geometry can reuse applicable work. Rotor motion, deformation, core changes and shedding-location changes invalidate relevant terms. Never reuse a factorization merely because panel counts match. |
| 6 | Preallocate aerodynamic snapshots and separate immutable geometry from history | Current deepcopy snapshots copy allocated wakes each time; repeated trials restore full storage. Optimize only after checking which scratch arrays must be recomputed and whether inactive rows need copying. Preserve rollback tests. |
| 7 | Profile direct induced-velocity and wake convection loops | Their all-to-all work grows quickly with wake length. Start with allocations, cache-friendly access and disjoint receiver parallelism. Tree/FMM or wake coarsening is a later numerical change requiring accuracy studies. No speedup factor is established here. |
| 8 | Cache identical structural models and run unique cases in separate worker processes when memory permits | Chordwise/radial aerodynamic sweeps often share the same structural operators. Isolate mutable workspaces and output directories. Current `redirect_stdout`/`redirect_stderr` usage makes naive threaded case execution unsuitable. Control BLAS oversubscription. |
| 9 | Decouple logging/plotting from the solve and add completed-case metadata/checkpoints | Buffered logging and preflight plotting avoid late failures. Aeroelastic runs currently always recompute across sessions; any reuse must verify operating speed, structural inputs, pulse, tolerances, source and histories. |

Before optimizing, measure warm-run wall time and allocated bytes separately for geometry update, AIC assembly/factorization, force calculation, wake convection, structural correction, snapshot/restore and output. Compilation time must be reported separately. The table identifies code-backed opportunities; only the surface-storage cost was quantified here.

**AeroBeams-compatible terminology and ownership**

The local root package is AeroBeams 0.8.1. Its local source is the authority for this proposal; public documentation may describe another revision. The [official project](https://github.com/luizpancini/AeroBeams.jl) describes geometrically exact beam structures with steady, trim, eigen and dynamic analyses. Those distinctions should remain visible in the bridge.

Use AeroBeams' existing `create_Model`, `create_Beam`, `create_DynamicProblem`, `create_NewtonRaphson` and `solve!` conventions for structural setup. Introduce a small bridge that owns the UVLM state and the mapping; retain the low-level backend behind it. Proposed names below are design suggestions, not APIs already implemented:

| Current concept | Proposed public concept | Meaning to preserve |
|---|---|---|
| Case NamedTuples mixed with solver orchestration | `create_UVLMModel(...)` and explicit study configuration | Geometry/physical data separate from analysis settings |
| Large `System` containing mixed data | `UVLMSystem` initially; then `UVLMState` and `UVLMWorkspace` | Circulation/wake history versus scratch matrices and force buffers |
| Chang-only free-vector mapping | `BeamUVLMCoupling` | Beam/node/element IDs, physical attachments, frame transformations and interpolation |
| `run_chang` and second-order matrices | A bridge-owned `CoupledDynamicProblem` containing an AeroBeams `DynamicProblem` | AeroBeams owns structural equations and structural time integration |
| `PartitionedCouplingOptions` snake-case fields | `create_PartitionedCoupling` with `maximumIterations`, `stateTolerance`, `loadTolerance`, `relaxation` | Match public naming; retain compatibility wrappers during migration |
| `time`, `dt`, step indices | `timeVector`, `Δt`, `timeNow`, `timeBeginTimeStep`, `timeEndTimeStep` | One owner of physical time and accepted-step state |
| `surface_history`, response arrays | `savedTimeVector`, `aeroStatesOverTime`, response result types | Consistent tracking policy and predictable memory use |
| `ns`, `nc` and ambiguous wake lengths | `nSpanwisePanels`, `nChordwisePanels`, `nWakeRows`, `wakeRevolutions` | Keep allocated capacity distinct from active wake length |
| `speed_mps`, RPM and flow angles | Explicit operating-point object; SI internals, conversions at boundary | Rotor angular speed in rad/s; RPM only as a user input/output convention |

Do not imitate broad `Real` or `Function` fields merely for naming consistency. Use concrete or parameterized storage where practical, while retaining the recognizable constructor and lifecycle conventions.

There are five essential semantic differences:

1. **State vector:** AeroBeams elemental states include `u,p,F,M,V,Ω,χ`; `Problem.x` also contains formulation-specific and constrained DOFs. It is not the Chang sequence of six nodal displacement/rotation values followed by two propeller modal coordinates. Use element/node IDs and reconstruction APIs; never slice `x` using the Chang order.
2. **Rotations:** Chang composes Euler rotations. AeroBeams' deformation rotations use Wiener–Milenkovic parameters and rotation/tangent utilities (`Utilities.jl:439`, `Core.jl` rotation update). Initial beam orientation has its own parametrization. Reconstruct rotation matrices with the existing utilities and map velocities using their tangent operators. Renaming Euler angles to `p` is incorrect.
3. **Frames and loads:** AeroBeams displacements are resolved in A; sectional F/M/V/Ω have documented B-frame meanings; aerodynamic nodal resultants are assembled in A. Define an explicit proper rotation between UVLM coordinates and AeroBeams A. Transfer dimensional vertex forces with their matching vortex positions. Do not feed normalized `body_forces` coefficients or convention-adjusted roll/yaw coefficient signs into structural residuals. Apply AeroBeams `forceScaling` once at residual assembly.
4. **Time integrator:** AeroBeams currently computes state rates through `2/Δt*state - equivalent_rate` in `Core.jl:155`, within its own `DynamicProblem` solve. Replace the Chang structural correction with that structural solve when coupling. Do not advance both generalized-alpha and AeroBeams for one physical step.
5. **Modal output units:** `Problem.jl:823` stores AeroBeams eigen frequencies as angular frequencies; its `dampings` are real eigenvalue parts. Chang reports Hz and lambda in s⁻¹. Convert angular frequency by `2π`; lambda has the same growth/decay sign meaning and is not a damping ratio. If desired, compute `ζ=-λ/sqrt(λ²+ω²)` explicitly.

For kinematic mapping, let aerodynamic vertex positions be `x_a=g(q_s,t)`. The transfer must satisfy `δx_a=H δq_s` and `Q_s=Hᵀ f_a`, with the appropriate rotational-coordinate tangent or physical nodal-wrench interpretation at the AeroBeams residual boundary. Preserve `f_aᵀδx_a=Q_sᵀδq_s`. This is the invariant already supported by the Chang virtual-work test. Decouple structural and aerodynamic meshes using this mapping rather than equating panel rows to beam nodes.

**Where the bridge must enter the solver**

`Core.jl:366` currently computes aerodynamics per element. UVLM circulation and wake influence are nonlocal across surfaces, so solving the whole lattice separately inside every `aerodynamic_loads!` call would both repeat work and risk advancing the wake multiple times per Newton iteration. A new subtype of `AeroSolver` may identify a surface model, but the global state cannot be owned independently by each `AeroProperties` object.

The proposed lifecycle is:

1. At the physical-step boundary, save accepted aerodynamic and structural states.
2. During a structural/coupling trial, reconstruct every coupled surface from the current AeroBeams configuration and velocities.
3. Restore the accepted aerodynamic snapshot and solve one circulation/load trial without committing wake convection.
4. Scatter work-conjugate loads to the appropriate AeroBeams nodal/element residuals. Explicitly disable duplicate sectional aerodynamics on those same surfaces.
5. Solve/update structural equations and test state, load and coupled residual convergence.
6. On success, commit circulation/history and advance the wake once. On failure or adaptive-dt retry, restore both disciplines and recompute at the new dt.
7. Save configured histories using the common time vector and outcome record.

For an initial partitioned implementation, provide trial/commit/rollback hooks around `solve_time_step!` and `time_march!` (`Problem.jl:1318`). Repeated structural trials must not mutate accepted AeroBeams history or advance the global frame repeatedly. State-dependent UVLM loads cannot be supplied as only a precomputed time-only load function: `precompute_distributed_loads!` runs before the solve and cannot represent an unknown current deformation.

For a later monolithic Newton method, assemble nonlocal aerodynamic derivative blocks across connected elements. Existing per-element derivative machinery in `Core.jl:844` does not automatically contain those cross-element couplings. The UVLM code explicitly rejects Imperial near-field derivatives, and several geometry APIs require Float64, so directly turning on AeroBeams AD is not sufficient. Start with a verified finite-difference Jacobian at fixed accepted wake history, then consider analytic derivatives or AD-compatible storage. Derivative evaluations must have identical rollback semantics.

Also, `AeroBeams.solve!(problem::Problem)` uses exact type branches (`Problem.jl:634`). A new `Problem` subtype will not acquire a solver merely by subtyping it. Either provide an explicit dispatch method in the bridge or extend the existing lifecycle deliberately. A wrapper that contains a real `DynamicProblem` is a safer first integration boundary than pretending the Chang model is an AeroBeams `Model`.

**Proposed migration sequence and acceptance gates**

| Phase | Concrete deliverable | Acceptance gate |
|---|---|---|
| A — correctness | Robust damping estimator, one operating-point API, isolated example dependencies, explicit termination | Synthetic off-bin/out-of-band/offset/multimode tests; missing-dependency CI; no silent speed mismatch |
| B — package organization | Move reusable convergence/model/solver code from examples into package modules; retain thin entry files | Existing 369 tests, new smoke and virtual-work tests remain valid; document output schema migration |
| C — performance | History controls, reusable workspaces, cached structural factorization | Numerical equivalence, bounded allocations, warm-run profile; no changed circulation/force/wake results without an explicit numerical-method change |
| D — geometry/load bridge | `BeamUVLMCoupling` with independently chosen structural/aero meshes and explicit frames | Identity/rigid-motion tests; resultant force/moment; virtual work at nonzero rotations, nonmatching meshes and off-node attachments |
| E — dynamic bridge | AeroBeams-owned structure/time, global UVLM trial/commit protocol | Rigid prescribed motion; structural-only limit; small coupled benchmark; repeat-trial determinism; forced rejection and smaller-dt retry; exactly one wake commit |
| F — analysis parity | Distinct steady trim, dynamic perturbation, eigen/linearized analyses and overall convergence selection | Mode identity, frequency-unit checks, time/mesh/iteration convergence; combined refinement; suitable external benchmark |

For phase B, source fingerprints need an explicit migration policy. They presently hash all convergence and Chang source files, so even a damping-only or structural reporting change invalidates rigid aerodynamic provenance. Split fingerprints by the numerical dependencies of each stage and version the schemas; never silently rewrite existing verified seals. Disk hashes also need to be associated with the code actually loaded in a long-lived Julia session, or runs must require a fresh process after source edits.

A Julia extension/weak dependency or separate bridge package can keep AeroBeams optional for standalone UVLM users. Since this checkout currently nests UVLM under AeroBeams without a root dependency, choose one dependency direction deliberately; avoid having each package require the other. Keep Chang tabulated properties and modal pylon assumptions as a benchmark builder, not requirements imposed on all beam models.

**How to reproduce this review**

From the repository root:

```powershell
julia --startup-file=no --compiled-modules=existing --project=lib/WingPropellerUVLM lib/WingPropellerUVLM/test/runtests.jl
julia --startup-file=no --compiled-modules=existing --project=lib/WingPropellerUVLM docs/reviews/uvlm_aerobeams_2026-09-14/reproduce_findings.jl
julia --startup-file=no --compiled-modules=existing --project=lib/WingPropellerUVLM docs/reviews/uvlm_aerobeams_2026-09-14/coupled_smoke.jl
```

The first command uses the user's normal depot; in the restricted review environment it was invoked through a wrapper that prepended a temporary depot. The smoke script already does that. The diagnostics require the currently available shared FFTW package to reproduce the estimator, then explicitly inspect isolated-project dependency resolution. No full aerodynamic convergence sweep, six-second damping study, flutter boundary, or actual AeroBeams–UVLM coupled solve was run as part of this review.
