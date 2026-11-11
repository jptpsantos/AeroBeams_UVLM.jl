# Chang implementation review — 2026-11-10

The follow-up [notation pass](../../CODE_STYLE.md) covers the whole custom
package, including the backend, studies, Pazy and nonlinear Chang coupling.
See [its validation report](../../NOTATION_VALIDATION.md). The line references
below describe the earlier Chang correctness review, before the notation pass.

The review covers the primary and fast Chang runners, configuration/model
construction, structural matrices and integration, UVLM load transfer, wake
convection, history retention, and result output. AeroBeams `src/AeroBeams.jl`,
`src/SystemSolver.jl`, `src/Problem.jl`, and `src/Utilities.jl` supplied the style
reference: explicit numerical operations, small local corrections, existing
data structures, and straightforward loops. This is a focused review of the
Chang execution path, not an audit of every AeroBeams example or study.

| Classification | Finding and implemented correction | Source location |
|:--|:--|:--|
| A — Confirmed bug | Shared trailing-edge wakes with unequal active lengths could pass an invalid vertex index into the shorter sending wake. Reproduced a `BoundsError` in the backend and the same error wrapped by threaded tasks in the fast implementation. Use direct induction at the receiving coordinate past the sending wake's last vertex. | `../../src/backend/wake.jl:258`; `src/ChangRuntimeOptimizations.jl:80` |
| A — Confirmed bug | A scalar comparison rejected arrays of per-propeller wake capacities before the existing normalization and elementwise validation. Removed the redundant comparison. Valid arrays now work; negative entries and incorrect lengths still fail. | `../../src/wing_propeller/Initialization.jl:69` |
| A — Confirmed bug | A state/angle limit stop committed the final aerodynamic state but skipped its requested surface-history entry. Save the accepted geometry before checking the stop limit. | `src/chang_simulation.jl:569` |
| A — Confirmed workflow bug | The primary runner required undeclared CSV/DataFrames packages even when the optional organized exporter was unused. Use NamedTuple columns and the same plain CSV-writing approach as the existing result writer. | `run_chang_linear_aeroelastic.jl:23` |
| A — Confirmed workflow bug | Organized histories ignored the configured output directory. Write them beneath that directory's `multi_prop/` child. | `run_chang_linear_aeroelastic.jl:63` |
| A — Confirmed workflow bug | The benchmark could fail while writing its report on a clean checkout because `output/validation/` had not been created. Create the parent before allocating a fresh report directory. | `validation/benchmark_chang_runtime_optimizations.jl:75` |
| A — Confirmed test defect | The baseline package suite demanded bitwise equality of damping matrices obtained through equivalent rotor-speed expressions. The observed entries were `433.73425849987865` and `433.7342584998787`. Use `rtol=10eps(Float64)` for that assertion; leave the production formulas intact. | `../../test/general_structural_assembly.jl:166` |
| C/D — Simplification and allocation reduction | Remove the accepted full-load copy; the trim accumulator can consume the original buffer before it is modified again. This removes one vector copy per accepted step without changing arithmetic. Reuse the already computed propeller ordering in export metadata. | `src/chang_simulation.jl:517`; `run_chang_linear_aeroelastic.jl:77` |
| C/E — Simplicity and AeroBeams conventions | Keep explicit file-output loops, named columns, existing array updates, and local branches. Remove an obsolete commented finite-core assignment. Avoid activating an already active package project, which otherwise unnecessarily writes Julia's global usage log. | `run_chang_linear_aeroelastic.jl:6`; `../../src/wing_propeller/Initialization.jl:50` |

The wake correction changes the result only for the confirmed unequal-length
shared-edge case, which previously failed or could read an inactive row. It
evaluates the existing Biot–Savart operator at the actual receiving vertex.
The current Chang configuration does not exercise that exceptional topology.
Mass, stiffness, damping, force/moment transfer, integrators, coupling formulas,
mesh settings, time steps, core laws, and tolerances in the solvers are unchanged.

The organized export retains its filenames, column names/order, angular units,
and the outer return fields `csv_path`, `txt_path`, and `df`. Its `df` field is
now a NamedTuple of column vectors. Column access such as `df.time_s` works;
callers using DataFrame-specific indexing must explicitly convert it.

Validation uses Julia 1.11.5 and four Julia threads. Commands were run from the
repository root with `--startup-file=no --compiled-modules=existing
--threads=4 --project=lib/WingPropellerUVLM`. The compilation flag permits
unchanged cached dependencies while compiling modified sources in memory,
because this workspace cannot write Julia's global compilation cache.

| Executed check | Result |
|:--|:--|
| `lib/WingPropellerUVLM/test/runtests.jl` | 775 assertions passed; process exit code 0 |
| `lib/WingPropellerUVLM/test/chang_runtime_optimizations.jl` | 1,245 assertions passed; process exit code 0 |
| `validation/verify_chang_runner.jl` | 12 assertions passed; reduced primary simulation, output-directory handling, and independent one-/two-propeller export ordering checks |
| `validation/audit_chang_virtual_work.jl` | Passed at reference, uniform, and nonuniform deformations; maximum scaled work error `1.1484226014601973e-7`, below the existing `1e-6` threshold |
| `validation/benchmark_chang_runtime_optimizations.jl` | All three repetitions had exactly equal velocities and shed/convected panels; median serial/threaded times 6.810/2.413 s, or 2.82× for wake commit only |
| `git diff --check` | Passed for the reviewed paths |

The final three test scripts passed 2,032 assertions in total. Logs are saved
under [the review output directory](output/validation/code_review_20261110_2656b261/).
The [bounded benchmark report](output/validation/fast_runtime_benchmark_nsRaMH/report.md)
records the mesh, wake capacities, individual times, and measurement limits.

The runtime suite checks both structural integrators and coupling schemes,
exact serial/threaded circulation, wake and structural histories, symmetry,
interaction groups, shared edges, full-history retention, animation storage,
output preservation, and limit stops. The new independent shared-edge oracle
isolates the shorter sending wake and compares its velocity with direct
induction at the physical receiving point, so backend/fast agreement alone
cannot hide the original bug. The original runtime baseline passed 875
assertions before these regressions were added.

Two numerical questions remain for separate investigation (classification B):

- The shared-edge optimization assumes matching row indices are coincident
  vertices for rows present on both wakes. Initially coincident surfaces that
  subsequently move independently need a dedicated topology/kinematics check.
  The present fix addresses vertices downstream of a shorter wake, without
  redefining shared-edge exclusion. See `../../src/backend/wake.jl:219` and
  `src/ChangRuntimeOptimizations.jl:65`.
- A propeller pulse at time zero, or a trim window that contains no accepted
  time point, can pass configuration loading and then encounter the solver's
  explicit no-trim-samples error. Supporting those cases needs a defined
  initial-baseline/averaging policy; inventing a fallback would change the
  perturbation model. See `src/chang_configuration.jl:278` and
  `src/chang_simulation.jl:513`.

The reduced checks establish execution and regression consistency. They do
not establish mesh/time convergence or a flutter boundary, and the full
five-second production transient and PNG/GIF rendering were not run. The
AeroBeams root test suite was not run; no AeroBeams source files changed.
No new speedup is attributed to
the removed vector copy; the benchmark measures the existing threaded wake
implementation. Unrelated Pazy output moves/deletions and the user's case
parameters were preserved.
