# Repository-wide custom-code notation pass

Scope: every existing custom Julia file under `lib/WingPropellerUVLM`,
including both packages, the backend, structural routines, all examples,
studies, tests, validation tools and archived custom scripts. The user's
choice was to preserve existing public names and saved-file schemas.

The primary references were the local AeroBeams `Beam.jl`, `Model.jl`,
`Element.jl`, `Core.jl` and `SystemSolver.jl`. The resulting conventions and
physical notation are documented in [CODE_STYLE.md](CODE_STYLE.md), with
instructions for subsequent work in [AGENTS.md](AGENTS.md).

| Source audit | Result |
|:--|:--|
| Custom Julia files parsed | 211 |
| Files changed from the start of this pass | 200 |
| Already consistent files retained | 11 |
| Recorded descriptive internal-name mappings | 956 |
| Pure name/layout changes verified by Julia syntax-tree comparison | 204 files |
| Study-entry guards verified separately | 2 files |
| Diagnostic/test-only adjustments reviewed separately | 5 files |
| AeroBeams source/test/development files with unchanged SHA-256 hashes | 502 |
| Renaming collisions within a file | 0 |

The descriptive mappings use `camelCase` for internal names and keep
`snake_case` functions. Additional reviewed renames make the core objects
explicit (`problem`, `workspace`, `model`, `aeroSolver`, `operatingPoint`,
`nSurfaces`) and use `r_*` for propeller position vectors. Existing mathematical
subscripts, rotation definitions and public fields/keywords are retained.
The two old step-counter spellings that become `nSteps` occur in separate
files and represent the same counting concept.

The syntax comparisons ignore source-line metadata and apply the recorded
identifier mappings. They preserve numerical expression trees, literal
values, accumulation order and public record keys. No equations, units,
integration/coupling algorithms, solver tolerances, discretization choices,
case inputs or runtime dependencies were changed by the notation pass.

The separately reviewed changes are:

- `chang_convergence_aerodynamic.jl` and
  `chang_convergence_aeroelastic.jl` now use a script-entry guard. Direct CLI
  execution still runs the selected study; `include` loads definitions without
  starting a production sweep. This makes the existing import/dry-run tests
  executable. The first test attempt exposed and was interrupted at the
  previous unconditional sweep call.
- `audit_pazy_failure.jl` recognizes the formatted save/frame calls when
  instrumenting a diagnostic copy; its hooks do not change production solvers.
- `test_pazy_LCO_output.jl` matches driver source independently of comma
  spacing and sets explicit fixture angles, speeds and branches. The user's
  editable driver values are preserved.
- `test/moving_block_damping.jl` expects the existing legacy API's
  `ErrorException` on invalid records, and separately checks the validated
  metrics API's `DimensionMismatch`. Numerical FFT/regression expectations
  are unchanged. The wider tests exposed these pre-existing expectation
  mismatches; the implementation and reference routine both use `error()`.
- `studies/test/migration_smoke.jl` now sets its fixed 1,212 RPM / 65 m/s
  reference operating point explicitly. Its original angular-speed expectation
  is retained, independently of the user's editable study RPM.
- `test/chang_aeroelastic_preview.jl` explicitly supplies the five-family,
  four-second full-study planning fixture used by its original assertions.
  The user's current full-study duration and enabled families are preserved.

Tests run with Julia 1.11.5, four Julia threads, headless plotting and
`--startup-file=no --compiled-modules=existing`. Modified package sources
are compiled in memory because Julia's global compilation cache is read-only
in this workspace. No dependency installation was required.

| Regression checks | Result |
|:--|:--|
| Full studies test entry, including package tests, convergence, damping and migration | 1,192 assertions passed; exit code 0 |
| Chang runtime/fast-runner checks | 1,245 assertions passed, including exact serial/threaded histories and unequal shared-edge wake lengths |
| Nonlinear Chang/AeroBeams coupling validation | 119 assertions passed |
| Pazy 11-script integration batch: geometry, virtual work, Newton recovery, output, transients and bifurcation plots | 854 passed; 2 legacy exception-type expectations failed before their correction; no numerical assertions failed |
| Corrected Pazy moving-block/stability subset rerun | 76 assertions passed, including the corrected expectations and the new typed-metrics check |
| Parse/name/layout and protected-source audit | Passed |
| `git diff --check` | Passed |

The [raw evidence directory](examples/chang_linear_aeroelastic/output/validation/notation_20261110_76208cfe/)
contains the pre-pass custom-source snapshot, the complete
[211-file hash/rename inventory](examples/chang_linear_aeroelastic/output/validation/notation_20261110_76208cfe/notation_source_validation.toml),
source-audit logs and regression logs. The initial Pazy batch log deliberately
retains its two failed legacy exception expectations; the separate corrected
subset log records the successful rerun. Those failures were test assumptions,
and no physical/numerical assertion failed in that batch.

The normal AeroBeams root test suite and full production parameter sweeps are
outside this notation pass. The original AeroBeams code, public API spellings,
case values and unrelated simulation-output moves/deletions are preserved.
Old study selections remain subject to their existing source-provenance
checks; the checks were not weakened to accept a different source snapshot.
