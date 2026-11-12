# AeroBeams upstream update assessment

Assessed on 11 November 2026. This is an assessment and migration procedure; the update has not been merged and no solver or case file has been edited.

**Recommendation.** Integrate AeroBeams 1.0 on a separate migration branch, preserve a runnable 0.8.1 baseline, and port the nonlinear structural adapters before using the updated branch for research results. The standalone UVLM and linear Chang implementations can remain structurally unchanged. A dependency update alone cannot perform this migration because this repository is itself an AeroBeams fork.

**Exact revisions inspected**

| Item | Revision / observation |
| --- | --- |
| Current branch | `update-aerobeams` |
| Current committed code | `55f06211` (`Phase analysis`) |
| Common ancestor with upstream | `69b2feeb9405f719b3abf349dc07894093c55700`, 2 July 2026 |
| Upstream fetched during assessment | `285efa818b70a610ab5d2fcfbece14340ec4bb72`, 10 November 2026 |
| Divergence from upstream | 38 commits only on this branch; 11 only on upstream |
| Local AeroBeams package | 0.8.1 |
| Upstream package | 1.0.0 |
| Julia invoked by the current shell | 1.11.5 |
| Upstream Julia compatibility declaration | `julia = "1.13"` |
| Existing working changes | `analyze_pitch_phase.jl` and `chang_case.jl` in `examples/chang_linear_aeroelastic/` |

The Git remotes are already correct: `origin` points to `jptpsantos/AeroBeams_UVLM.jl`, and `upstream` points to `luizpancini/AeroBeams.jl`. Fetching refreshed the cached upstream reference from `7751b374` to `285efa81`. The latest commit only disables scheduled CompatHelper; the substantial migration is the previously unmerged October `refactor for v1.0` work.

Authoritative comparison: [upstream commit](https://github.com/luizpancini/AeroBeams.jl/commit/285efa818b70a610ab5d2fcfbece14340ec4bb72), [upstream package requirements](https://github.com/luizpancini/AeroBeams.jl/blob/285efa818b70a610ab5d2fcfbece14340ec4bb72/Project.toml).

**What the update affects**

| Local area | Assessment | Required work |
| --- | --- | --- |
| `lib/WingPropellerUVLM/src/` | Low direct API exposure | Keep the UVLM physics and public API; run existing regressions under the new Julia/environment |
| `examples/chang_linear_aeroelastic/src/` | Low direct exposure | This uses the local linear structural matrices and UVLM, without importing AeroBeams; do not replace it with the nonlinear solver as part of an upstream update |
| `examples/chang_linear_aeroelastic/analyze_pitch_phase.jl` | No direct AeroBeams API dependency | The saved-history fit uses DelimitedFiles, LinearAlgebra, Statistics, Printf, Plots and ColorSchemes; preserve its input columns and phase definitions |
| `examples/pazy_aerobeams_coupled/` | High | Port construction, load injection, structural stepping, state access, static coupling, histories, animation, and diagnostics |
| `examples/chang_aerobeams_coupled/` | High | Port wing/pylon construction and the custom coupled residual, Jacobian, DOF maps, time state and rollback |
| `lib/WingPropellerUVLM/studies/` | Environment/provenance impact | Revalidate dependencies, regenerate its manifest only as needed, and begin new studies with the new Julia/source identity |
| Root `src/` and `test/` | Major upstream refactor | Integrate the v1 implementations and review the local optimization patches individually |
| Documentation/workflows/development scripts | Repository maintenance impact | Reconcile examples and paths, preserve useful local research material, and review fork CI/release behavior |

A textual inventory found qualified AeroBeams references in 32 custom Julia files, excluding the `backup`, `output`, and `archive` directories. This includes comments and diagnostic references; it is an exposure inventory, not a claim that every file needs edits. Indirect field access requires inspection as well.

The current saved linear Chang histories can still be analyzed with the current code and environment. Updating AeroBeams does not retroactively change those files or improve their model assumptions. Rerunning the nonlinear models after migration can change the response, and requires fresh validation.

**Merge conflicts and local AeroBeams patches**

Upstream changes all 24 files currently under `src/`. File-level three-way merge checks used the actual common ancestor, committed local files and fetched upstream files, without modifying the working tree or index. They predict these conflicts:

| File | Conflict blocks |
| --- | ---: |
| `src/AeroBeams.jl` | 1 |
| `src/Core.jl` | 8 |
| `src/Element.jl` | 1 |
| `src/Problem.jl` | 3 |
| `src/SystemSolver.jl` | 6 |
| `src/Utilities.jl` | 3 |
| `test/runtests.jl` | 1 |
| Total | 23 |

The overlapping changes in `README.md`, `docs/make.jl`, and `docs/src/index.md` merge cleanly in this check. This predicts textual conflicts; it does not establish a compiling or numerically compatible merge. In particular, custom `lib/` files can merge cleanly because upstream does not have them, while retaining invalid calls to the old AeroBeams API.

The six locally modified AeroBeams source files mostly contain efficiency and history-isolation work: removal of Revise from module loading; views and fused arithmetic; concretely dimensional distributed-load arrays; avoiding unused load histories; cheaper Newton checkpoints; reusable LU factorizations; guarding spring-Jacobian resets; copying saved state vectors; cached Gauss rules; and reduced matrix temporaries.

Several objectives now have upstream equivalents: Revise is removed, state/cache containers are reorganized and specialized, distributed-load allocation is conditional, and `SystemArrays.luCache` reuses sparse factorizations. Do not transplant the old cache/checkpoint machinery wholesale into v1. Compare behavior, measure the new implementation, and port only a demonstrated remaining need. Rewrite `test/efficiencyTests.jl` around valid v1 invariants before adding it back to the new upstream runner; its old fields and helper functions no longer exist.

The very large overall Git diff is partly repository cleanup: upstream deletes 6,035 tracked files under `dev/`, including generated material. Review those deletions if any local workflow uses them. That count does not mean thousands of custom UVLM source files need migration.

**Concrete API changes to port**

Source references below are to the pinned [v1 source directory](https://github.com/luizpancini/AeroBeams.jl/tree/285efa818b70a610ab5d2fcfbece14340ec4bb72/src).

| Existing usage | v1 representation / action |
| --- | --- |
| `problem.x`, `.residual`, `.jacobian`, `.inertia` | `problem.systemArrays.x`, `.residual`, `.jacobian`, `.inertia` |
| `problem.timeNow`, `.Δt` | `problem.time.state.timeNow`, `.Δt`; use `timeIndex` in the new time state |
| `problem.timeVector` | `problem.time.timeVector` |
| `problem.savedTimeVector` | `problem.solutionHistory.t` |
| `problem.nodalStatesOverTime` | `problem.solutionHistory.states[k][e].nodal` for a saved element state |
| `model.elements`, `.specialNodes`, `.nNodesTotal`, `.systemOrder`, `.forceScaling` | Corresponding fields in `model.assemblyData` |
| `element.states` and `.compStates` | `element.state.midpoint`; the new midpoint includes strains and momenta |
| `element.nodalStates` | `element.state.nodal` |
| `element.statesRates` | `element.state.rates` |
| `element.r_n1`, `.r_n2`, `.R0_n1`, `.R0_n2` | Geometry lives in `element.geometry` |
| Element global IDs / node IDs / DOFs | Inspect `element.assemblyData` and `element.assemblyIndices`; verify index shapes and ordering |
| Special-node equation and DOF maps | `specialNode.assemblyData`; prescribed flags live in `specialNode.flags` |
| `solver.convergedFinalSolution` | `solver.loadStepState.convergedFinalSolution` |
| `assemble_system_arrays!(problem)` | `fj!(problem)` is the new combined residual/Jacobian entry; also inspect `f!`, `j!` and the solver skip flags |
| `solve_NewtonRaphson!(problem)` | Nonlinear dispatch now goes through `solve_nonlinear_system!`; use the public solve workflow when the custom coupling permits it |
| `update_time_variables!(problem,i)` | `update_time_state!(problem,i)` |
| `update_basis_A_orientation!` | Reconcile with `update_motion_basis_A!`; the old `update_R_A_array` keyword is gone |
| `get_equivalent_states_rates!` | `update_equivalent_states_rates!` |
| `save_time_step_data!(problem,t)` | `save_time_step_data!(problem)`; history is preallocated and has its own index |
| `save_load_factor_data!` | Reconcile with `save_loadstep_data!` and `loadstepHistory` |
| `copy_state` / `restore_state!` | No same-named helpers in v1; design and test the coupled driver's accepted/trial checkpoint explicitly |

These are implementation directions, not a compatibility shim already present in the repository. Some transformations change lifecycle and storage semantics and cannot be implemented safely by global search/replace.

Four changes deserve particular attention:

1. **Rotation return value and initial orientation.** `rotation_tensor_WM(p)` now returns the rotation matrix directly, instead of `(R,...)`. Replace old `first(...)`, `[1]`, and tuple unpacking at AeroBeams interfaces. Selecting `[1]` now returns a scalar. `create_Beam` also removes the `rotationParametrization` keyword. Convert an existing Euler input explicitly, for example `p0=AB.transform_to_WM([pi/2,0.0,pi]; from=:E321)`, and verify that the resulting matrix equals the current reference frame. `initialPosition` becomes a new placement convention involving `firstNodeOffset`; confirm actual assembled nodal coordinates rather than assuming the meanings are identical. `aeroSurface` becomes `surface`.

2. **Boundary conditions and external loads.** `create_BC` accepts `conditions=Dict(:u1A=>0.0,...)`, replacing parallel string `types` and numeric `values` inputs. Definitions, functions, state and histories now live in separate containers; stored definition values are functions. `bc.values .= loads` and `update_BC_data!` cannot be retained. Introduce a narrow load adapter using the v1 BC state/functions, and explicitly refresh the active BC state for every trial load. The normal v1 time loop reads precomputed BC histories, which must not overwrite the newly evaluated UVLM trial load. `update_BC_state!(bc,timeNow::Float64)` evaluates functions, while the integer overload reads a history index; distinguish them deliberately. Test loads changing between successive Newton/FSI solves at the same physical time.

3. **Pazy sample data.** `nodal_positions_Pazy()` is removed. `geometrical_properties_Pazy()` now returns five values, adding the normalized node locations as the fifth. Update consumers to obtain `ne,span,chord,spar,nodes` from this helper. The old four-variable destructuring may ignore the extra value, but the separate old node-position call still fails. Inspect all Pazy sweep and plotting wrappers as well as the main model builder.

4. **Nonlinear Chang assembly and connections.** `BeamPylon.jl` passes `connectedBeams`, `connectedNodesThis`, and `connectedNodesOther` to `create_Beam`; these are removed. Model construction now accepts `connections`, expressed as beam-pair and node-pair tuples. Move the connection definition there and verify the resulting global IDs. Equation maps also change shape: v1 special-node force/moment equilibrium rows are fixed three-component index vectors. Do not retain the old flattening logic without checking it. The custom augmented residual/Jacobian and physical force scaling must still produce the same hub reactions and virtual work.

**Initialization, stepping and solver settings**

The public v1 `solve!` precomputes motion of basis A, BC histories, distributed loads and flap data, then initializes the current motion and BC state before dispatching the analysis. The existing custom drivers bypass much of the public solve function. Their current sequence of only `precompute_distributed_loads!` plus `solve_initial_dynamic!` is insufficient as a general v1 initialization procedure. Port the preparation sequence explicitly, while ensuring UVLM loads remain trial-dependent.

`create_DynamicProblem` removes `saveInitialSolution`. The new problem stores time and solution history separately. Saving the initial state and subsequent accepted steps must follow the new tracking/preallocation rules; enabling tracking only after construction is not an equivalent operation. Preserve the existing convention of forming equivalent structural rates once per physical step and advancing the aerodynamic wake once per accepted step.

| Existing Newton keyword | v1 keyword / treatment |
| --- | --- |
| `maximumIterations` | `maxIter` |
| `desiredIterations` | `desIter` |
| `absoluteTolerance` | `absTol` |
| `relativeTolerance` | `relTol` |
| `initialLoadFactor` | `initLoadFac` |
| `minimumLoadFactor` | `minLoadFac` |
| `maximumLoadFactorStep` | `maxLoadFacStep` |
| `minimumLoadFactorStep` | `minLoadFacStep` |
| `trackingLoadSteps` | `trackLoadSteps` |
| `alwaysUpdateJacobian=true` | Explicitly select/review `method=:Newton`; reuse behavior is controlled differently in `:NewtonChord` |

The default problem solver is now `create_NLsolver()`, rather than relying on the old Newton defaults. Preserve the intended solver explicitly. Defaults for iteration counts, load increments, error limits, aerodynamic convergence allowances and trim solution have changed. In v1's Newton loop the absolute norm is an infinity norm and acceptance checks both absolute and relative criteria; the old loop accepts either criterion. Matching tolerance numbers alone does not prove matching convergence behavior.

The nonlinear Chang driver also restores old clock fields and calls the removed `copy_state`/`restore_state!` helpers. Its new checkpoint must cover all mutable structural trial state, equivalent rates/caches, solver state, time, and external load state that affect a retry. A rollback test should compare a retried step with a clean run from the accepted state. Resetting the initial model state is not a substitute for restoring the accepted state.

**Environments, outputs and workflows**

Use Julia 1.13 or another version satisfying the upstream compatibility declaration before resolving the v1 environment. Lowering `compat.julia` locally does not establish compatibility. No migration runtime result was obtained with Julia 1.13 in this assessment.

Upstream adds ADTypes, ArnoldiMethod, Arpack, FunctionWrappers, KrylovKit, LineSearches, LinearMaps, NLsolve, Roots and StaticArrays as direct dependencies. It removes several old direct dependencies, including BenchmarkTools, Dierckx, GR, Interpolations, MAT, PyPlot, Revise and Statistics. It raises JLD2 from the 0.5 line to 0.6 and PrettyTables from 2 to 3, and tightens DSP. Keep packages explicitly required by local analysis code in the appropriate project, rather than relying on old transitive availability. WingPropellerUVLM already declares its own Interpolations dependency, so removing it from the root project is not itself a reason to remove it from UVLM.

Regenerate environments with Julia's package manager after updating their projects. Do not hand-edit manifests or run an indiscriminate `Pkg.update()` first. The root `Manifest.toml` exists locally but is untracked; the UVLM and study manifests are tracked. Preserve baseline copies and resolve the root, UVLM and study environments independently. Their manifests are distinct environments, not three copies of the AeroBeams version. The coupled examples currently stack checked-out packages through `LOAD_PATH`; verify `pathof(AeroBeams)` resolves to this checkout, not a registry-installed copy with the same UUID.

The Pazy visualization helper edits old beam aerodynamic fields, constructs old `AeroProperties`, and consumes old structural history arrays. It requires a substantive port; restoring simulation execution alone does not restore animation. Prefer v1-supported surface construction/attachment paths and preserve zero native aerodynamic forcing when UVLM supplies the loads.

Archive existing CSV/TXT histories and their case/source information. They remain readable by matching postprocessors. Treat `.jls` snapshots that contain old AeroBeams objects as version-bound: do not assume they deserialize into v1's different type layouts. Keep the old branch/environment for replay. New runs need new output directories and an AeroBeams revision/Julia/environment record. The existing linear-study source fingerprint includes Julia's version and UVLM/study sources, but does not fingerprint the root AeroBeams source; it is not sufficient evidence for nonlinear AeroBeams provenance.

Review upstream workflow changes as part of integration. CI now targets pushes to `main` and tags rather than all feature branches, so ensure migration PRs/custom checks cover the intended branch. The new registration workflow reacts to version changes on `main`; configure its behavior intentionally for this fork. Add UVLM/study/coupling checks separately: upstream's test runner does not know about them.

**Recommended integration procedure**

The following commands are instructions for the migration, not commands executed during this assessment. Start in the repository root. Commit the two existing working changes once they represent the state you want to preserve, or save them deliberately before proceeding. Do not use `git add .` to capture generated results.

```powershell
git status --short
git add -- lib/WingPropellerUVLM/examples/chang_linear_aeroelastic/analyze_pitch_phase.jl lib/WingPropellerUVLM/examples/chang_linear_aeroelastic/chang_case.jl
git commit -m "Preserve current Chang analysis before AeroBeams v1 migration"
git branch baseline-aerobeams-0.8.1
git switch -c aerobeams-v1-migration
git fetch upstream main
git merge --no-commit --no-ff 285efa818b70a610ab5d2fcfbece14340ec4bb72
```

The merge deliberately targets the assessed revision. Review newer upstream changes separately if the fetched branch has advanced. Resolve the seven expected conflicting files by understanding the v1 structures and the local patches. Preserve custom `lib/` code and scientifically useful local files. `git merge --abort` returns to the pre-merge migration branch while the merge is still in progress. Finish and commit a coherent integration after conflict resolution; the baseline branch remains available for comparisons.

Then work in reviewable steps: establish v1 loading and environments; port structural constructors and coordinate/state access; port external-load injection and stepping; port nonlinear Chang augmentation and rollback; restore histories/visualization/diagnostics; run regressions and scientific comparisons. Avoid combining a solver update with unrelated mesh, force-model, units, or case-default changes.

For the checked-out root package, use a fresh Julia process and verify the loaded source:

```powershell
julia +1.13 --project=. -e 'using Pkg; Pkg.resolve(); Pkg.instantiate()'
julia +1.13 --project=. -e 'using AeroBeams; println(pathof(AeroBeams)); println(VERSION)'
```

These commands assume a Juliaup `1.13` channel is installed. Use the corresponding Julia executable if managing versions another way. Resolve/instantiate UVLM and studies in their own projects as needed, preserving and reviewing their manifest diffs.

**Validation required before accepting new research results**

1. Run the new upstream structural, aerodynamic, eigenvalue, trim and dynamic tests. Adapt the local efficiency tests without weakening their scientific checks.
2. Run the standalone UVLM and study suites. Compare the linear Chang structural matrices, natural frequencies, short responses and saved-history analysis against the preserved baseline. Keep timestep, mesh, core model, RPM, force model, excitation and estimator settings fixed.
3. Run Pazy coordinate/load-transfer and virtual-work tests; verify that updating loads at the same time changes the structural residual. Compare dry structural modes, static deformation and short forced responses before long aeroelastic runs.
4. Validate both nonlinear Chang support configurations, including beam-pylon connections, DOF maps, hub reactions and gyroscopic coupling. Check residual/Jacobian consistency independently and exercise rejected-step rollback.
5. Check output sampling, saved-state isolation, final history frame, and animation/load vectors. Keep the output field names expected by existing postprocessors where possible.
6. Only then repeat representative Pazy stability/LCO points and nonlinear Chang growth/phase cases. Compare frequencies, growth rates, amplitudes and convergence behavior, and refine timestep/mesh when discrepancies persist. A passing implementation test does not establish a converged flutter boundary.

Existing entry points, to run after their code has been ported and environments are ready:

```powershell
julia +1.13 --project=. test/runtests.jl
julia +1.13 --project=lib/WingPropellerUVLM lib/WingPropellerUVLM/test/runtests.jl
julia +1.13 --project=lib/WingPropellerUVLM/studies lib/WingPropellerUVLM/studies/test/runtests.jl
julia +1.13 --project=lib/WingPropellerUVLM lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/test/runtests.jl --full
```

Also use the existing nonlinear Chang `test_reduced_rotor.jl`, `test_tip_rotor.jl`, `validation/test_coupled_solver.jl`, `validation/test_aerodynamic_interface.jl`, and `validation/compare_wing_modes.jl` as applicable. Some of those tests require the same v1 API migration as the implementation.

**Assessment limits and completed checks**

Completed: live upstream SHA verification, fetch, ancestry/divergence inspection, local/upstream source comparison, file-level three-way conflict prediction, qualified custom-call inventory, and inspection of constructors, BC caches, time integration lifecycle, output histories, workflows and package requirements. No merge, package resolution, Julia upgrade, simulation, or numerical regression was performed. The findings establish that a direct update breaks current nonlinear adapters; they do not establish the magnitude or direction of any change in validated flutter/LCO predictions.
