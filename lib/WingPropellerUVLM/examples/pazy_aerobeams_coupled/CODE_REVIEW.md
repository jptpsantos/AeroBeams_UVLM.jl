# Pazy coupling review and test cleanup

This review covers the custom Pazy dynamic and static coupling, visualization,
study runners, and their tests. It follows the package's
[AeroBeams conventions](../../CODE_STYLE.md), preserving public keyword and
result names. Original AeroBeams source files and editable case settings are
outside the changes made in this review.

## Corrections

1. **Final UVLM animation state.** Frame selection previously used only the
   stride multiples. A run whose step count was not divisible by that stride
   omitted its final accepted state. The dynamic loop now records stride
   multiples and the final state. The startup/pulse regression uses a frame
   stride that leaves a remainder and checks the final timestamp.
2. **Animation output path.** A string-prefix check accepted a sibling directory
   whose name started with the repository name. The helper now checks the
   relative path components. Tests cover an internal path, a parent path,
   a similarly named sibling, and a path containing `..` components.
3. **Static load-transfer allocation.** The wake-settling loop now constructs
   `PazySpanMap` once and reuses its nodal-load matrix through `beam_loads!`.
   The same transfer operations are used as in the allocating helper.
   Convergence-window entries still receive independent copies, so buffer
   reuse cannot overwrite the previous wake-step loads. No speedup is claimed
   without a benchmark.

## Definitions and documentation

- Positions use `r`; AeroBeams rotation parameters use `p` and are interpreted
  with `rotation_tensor_WM`, separately from Euler-angle inputs.
- `A_TO_UVLM` converts AeroBeams vectors to UVLM axes. Loads use its transpose
  on return, and transferred moments include the moment arm crossed with force.
- Geometry and load transfer share the interpolation weights and vortex offsets.
  Existing finite-rotation virtual-work and force/moment balance tests remain.
- The README now describes the default strong coupling, the optional loose
  predictor, numeric boundary-condition updates, and the actual relation between
  timestep and nominal wake-row length. These documentation corrections do not
  change the coupling equations, tolerances, or user-selected conditions.

## Test organization

The twelve existing automated scripts moved into `test/`, with their include
paths adjusted. `test/runtests.jl` runs them in isolated modules, with a default
suite and `--full` integration coverage. The shared moving-block damping suite
remains in the package tests instead of being included again by the Pazy
stability test. The saved-result animation preview moved to `validation/`.
No simulation result folders or backups were removed.

See [test/README.md](test/README.md) for commands and opt-in physical, distributed,
and rendering checks. Test results for this review are recorded below after
verification. Short regression cases do not establish mesh, timestep, wake,
or flutter-boundary convergence for a production study.
