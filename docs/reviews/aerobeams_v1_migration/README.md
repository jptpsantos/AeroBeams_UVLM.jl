# AeroBeams 1.0 migration record

Migration date: 12 November 2026. The pre-update assessment is
[preserved separately](../aerobeams_update_2026-11-11.md).

**Branches and preserved baseline**

The baseline branch is `baseline-aerobeams-0.8.1`, pointing to the user's saved
commit `730e4af6`. It includes the saved Chang case and pitch-analysis changes.
The migration branch is `aerobeams-v1-migration`, integrating upstream commit
`285efa818b70a610ab5d2fcfbece14340ec4bb72` (AeroBeams 1.0.0).

The previously untracked root environment is preserved as
[`baseline-Manifest.toml`](baseline-Manifest.toml). Its original Julia version
is 1.11.5. The migration installs Julia 1.13.3 alongside it and sets a Juliaup
override for this repository directory. The global Julia default is unchanged.
An already-running Julia process retains its loaded code and runtime; restart
the IDE Julia session before loading the migrated code.
The local, ignored `.vscode/settings.json` now selects the shared examples
environment using `${workspaceFolder}`; its previous value is preserved in
[`baseline-vscode-settings.json`](baseline-vscode-settings.json).
[`migrated-root-Manifest.toml`](migrated-root-Manifest.toml) records the resolved
root validation environment, since the upstream package ignores its root manifest.

The upstream `src/` implementations are retained. The local performance tests
are ported to v1's native storage, sparse-factorization cache and copied saved
states, rather than carrying the old 0.8.1 solver/cache implementations into
the refactored core.

**Implementation changes**

- Pazy and nonlinear Chang construction now uses v1 beam geometry, boundary
  conditions, solver keywords, model connections and sample-data helpers.
- Structural arrays, element states, geometry/assembly data, solver convergence
  flags, time state and histories use the new v1 containers.
- `AeroBeamsCoupling.jl` prepares the prescribed-data caches needed by the
  custom time loops, updates the active nodal wrench for every aerodynamic
  trial, and checkpoints/restores the mutable state used by a rejected Chang
  step. Its checkpoint preserves the active model graph and avoids copying
  the entire accepted solution history.
- The augmented Chang rotor residual uses v1's single force/moment balance
  per shared node. Pylons attach to the already assembled wing through model
  connection pairs. The original force scaling and physical wrench conventions
  are retained.
- Chang's two exact half-turn reference orientations use explicit WM
  parameters. This preserves the original director matrices and avoids an
  upstream Euler/quaternion conversion that can evaluate `asin` slightly
  outside its domain through roundoff at a half-turn.
- Structural animation builds a display-only model with the surface attached
  at construction. It reads selected accepted states, includes the final
  frame, and leaves the simulation model and histories untouched. The v1
  supported backends are `:gr` and `:plotlyjs`; the structural animation
  default is `:gr`.
- The Pazy comparison plot reads v1's consolidated `damps.txt`, preserving its
  mode-by-speed, column-major layout. Windows directory paths are normalized
  without a trailing separator before containment/relative-path checks.
- The linear Chang numerical solver and saved-history phase-fit algorithms
  are retained. The dynamic case loader obtains the freshly included case
  function in the latest world before calling it, as required by newer Julia
  world-age semantics. Case values, public configuration keys and CSV columns
  are preserved.
- The root keeps an explicit Statistics dependency for existing analysis
  scripts. The study environment uses compatible FFTW/plotting dependencies.
- A custom CI workflow exercises UVLM, studies and both structural couplings.
  The upstream registration workflow is guarded to run only on the upstream
  repository.

The main interface mappings for future custom scripts are:

| AeroBeams 0.8.1 usage | AeroBeams 1.0 usage |
| --- | --- |
| `problem.x`, `residual`, `jacobian`, `inertia` | `problem.systemArrays.<field>` |
| `model.elements`, `specialNodes`, `forceScaling` | `model.assemblyData.<field>` |
| `element.states`, `nodalStates`, `statesRates` | `element.state.midpoint`, `.nodal`, `.rates` |
| Nodal fields such as `u_n1`, `p_n2_b` | `u1`, `p2_b` in the nodal state |
| Element DOF/equilibrium indices | `element.assemblyIndices.<field>` |
| `problem.timeNow`, `timeIndex`, `Δt` | `problem.time.state.<field>` |
| `problem.timeVector` | `problem.time.timeVector` |
| Saved nodal state arrays | `problem.solutionHistory.states[k][e].nodal` |
| Static saved load-factor states | `problem.loadstepHistory.states` |
| `create_BC(types=...,values=...)` | `create_BC(conditions=Dict(:u1A=>0.0,...))` |
| `maximumIterations`, `absoluteTolerance`, `relativeTolerance` | `maxIter`, `absTol`, `relTol` |
| `alwaysUpdateJacobian` | `method=:Newton` or `:NewtonChord` |
| `assemble_system_arrays!(problem)` | `fj!(problem)` for residual and Jacobian |
| `R,_ = rotation_tensor_WM(p)` | `R = rotation_tensor_WM(p)` |

Geometry, reference positions and assembly mappings are distinct v1
containers; a mechanical replacement of every old field with one common
prefix would be incorrect. Likewise, updating a load definition alone does
not update the active load consumed by a custom trial. These changes are
handled explicitly at the coupling interfaces.

**Environments**

The coupled examples now have one reproducible environment at
`lib/WingPropellerUVLM/examples/Project.toml`. Its sources refer to the two
packages in this checkout. Resolving them together exposes the CCBlade 0.2 /
ImplicitAD 0.3 dependency on ForwardDiff 0.10, which is incompatible with
AeroBeams 1.0's ForwardDiff 1 requirement. WingPropellerUVLM therefore allows
both CCBlade 0.2 and 1; the shared environment selects CCBlade 1 and ImplicitAD
1. Standalone environments retain their own compatible locked versions.
The shared environment constrains ReverseDiff to 1.17.3: the 1.18.9
`@grad_from_chainrules` macro generates duplicate `func` arguments when loading
ImplicitAD 1.1.1. The compatible release supports ForwardDiff 1 and loads both
packages without modifying third-party package sources.
The [CCBlade 1 package declaration](https://github.com/byuflowlab/CCBlade.jl/blob/master/Project.toml)
documents the modern differentiation support.

Use the shared environment for the nonlinear Pazy/Chang examples. Mixing two
independently resolved manifests through `LOAD_PATH` can otherwise select
incompatible versions of shared packages or their extensions.

From the repository root:

```powershell
julia +1.13 --project=lib/WingPropellerUVLM/examples -e 'using Pkg; Pkg.instantiate()'
julia +1.13 --project=lib/WingPropellerUVLM/examples lib/WingPropellerUVLM/examples/chang_aerobeams_coupled/run_chang_aerobeams.jl
julia +1.13 --project=lib/WingPropellerUVLM/examples lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/run_pazy_wing_uvlm.jl
```

First-time package precompilation can be slow on Windows. The validation can
also use `--compiled-modules=existing` to use available caches and load missing
ones from source, avoiding stalled third-party precompile background tasks.
This does not change solver settings or numerical acceptance criteria.

For the active saved-history analysis, use the shared examples environment
(or the root environment) and a fresh Julia 1.13 session. The source history
files remain unchanged. `verify_pitch_analysis.jl` replays the reviewer
analysis into a separate QA directory and compares its results with the saved CSV.

**Validation evidence**

Validation is being completed in the final resolved environments. Final
results and any upstream-only failures are recorded here before delivery.

Completed independent checks:

- Standalone UVLM regressions passed under Julia 1.13.3.
- The study/convergence/damping suite passed, including 13 executable-example
  and Chang-facade parity checks.
- All 40 migrated native efficiency/state-ownership checks passed.
- The first six dry-wing frequencies at 8, 16 and 32 elements matched the
  committed AeroBeams 0.8.1 log within approximately `1e-14` relative difference.
  The verification requires all three meshes and uses the archived values,
  rather than recomputing the expected values with v1.
- The physical nonlinear beam-pylon/rotor check passed all 11 assertions.
- Animation-path and stability-reference reader checks passed all 24 assertions.
- A synthetic structural animation and its final PNG were rendered and
  inspected. Source-model/state preservation checks passed. These are rendering
  QA artifacts, not new physical simulation results.

Relevant scripts and artifacts:

- [Baseline modal comparison](verify_baseline_modes.jl), [results](baseline-modes.log).
- [Animation QA](verify_animation.jl), [final frame](render_qa/pazy_synthetic_final.png).
- [Saved pitch-analysis replay](verify_pitch_analysis.jl), [results](pitch-analysis.log).
- [Isolated upstream Helios comparison](verify_upstream_helios.jl).
- [Remaining upstream suites](verify_upstream_remaining.jl).

The original upstream test runner reports two failures in the Helios payload /
stiffness reference comparisons (`test/eigenAeroelasticTests.jl:83–84`). The
isolated check reproduces both with unchanged upstream source, test helper,
reference data and tolerances (`atol=1e-6`, `rtol=1e-4`):

| Output | Actual/reference shape | Different NaN locations | Maximum finite absolute difference |
| --- | --- | --- | --- |
| Frequencies | `(3,612)` | 0 | `0.03905772252835504` |
| Dampings | `(3,612)` | 0 | `0.002071262830894094` |

The cause of these Windows/environment reference differences remains
undetermined. This is not evidence that the entire upstream suite passes.
No golden data or tolerance is changed to obtain a passing result. The native
solver implementation remains identical to upstream; these comparisons are
separate from the Pazy and Chang adapter regressions.

**Returning to the baseline**

Preserve the old manifest outside the checkout before switching branches,
because the migration report does not exist on the baseline branch:

```powershell
Copy-Item docs/reviews/aerobeams_v1_migration/baseline-Manifest.toml "$env:TEMP/AeroBeams-baseline-Manifest.toml"
Copy-Item docs/reviews/aerobeams_v1_migration/baseline-vscode-settings.json "$env:TEMP/AeroBeams-baseline-settings.json"
git switch baseline-aerobeams-0.8.1
Copy-Item "$env:TEMP/AeroBeams-baseline-Manifest.toml" Manifest.toml
Copy-Item "$env:TEMP/AeroBeams-baseline-settings.json" .vscode/settings.json
juliaup add 1.11.5
julia +1.11.5 --project=.
```

The explicit `+1.11.5` overrides this directory's Julia 1.13 selection.
Use the preserved Julia 1.11.5 runtime for an exact baseline replay. Switching
branches does not switch an already-running Julia session. CSV/TXT results
remain suitable for their matching postprocessors; old serialized `.jls`
objects retain the old AeroBeams type layouts and should be replayed on the
baseline environment.
