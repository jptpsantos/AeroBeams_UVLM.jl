# Pazy verification

From the repository root, with the AeroBeams and WingPropellerUVLM environments
and the shared `lib/WingPropellerUVLM/examples` environment instantiated under
Julia 1.13 or newer:

```powershell
julia +1.13 --project=lib/WingPropellerUVLM/examples lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/test/runtests.jl
julia +1.13 --project=lib/WingPropellerUVLM/examples lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/test/runtests.jl --full
```

The default suite checks coordinate and load transfer, virtual work, Newton
load updates, animation output paths, stability statistics, spectrograms,
output routing, bifurcation analysis and plots. It includes a small coupled
bifurcation smoke test. `--full` also runs the startup/pulse coupling check,
static equilibrium, independent and continuous LCO sweeps, and diagnostics.
These are implementation regressions, not flutter convergence studies.

Every test file can also run directly. Existing optional switches remain:

- `test_coupling.jl --dynamic`: coupled startup/pulse and final animation frame.
- `test_pazy_model_options.jl --dynamic`: native matrices, rod configurations,
  swept geometry and load conservation, plus short coupled runs with AePW3,
  a refined no-rod mesh, 20-degree sweep and a 10-degree explicit LE rod.
  Also checks static-to-dynamic initialization and rejects mismatched model
  settings. Add `--animations` to save swept-rod preview GIFs under `tmp/`.
- `test_pazy_stepped_bifurcation.jl --physical`: stepped physical continuation.
- `test_pazy_bifurcation.jl --parallel`: distributed sweep consistency.
- `test_pazy_stability_sweep.jl --plots`: stability plot smoke check.
- `test_static_acceleration.jl`: Anderson mixing, wake reuse, geometry cap and warm starts of the
  static solver and of the static sweep reach the same equilibria as the original fixed-relaxation
  algorithm with fewer UVLM steps. Included in `--full`.
- `test_static_spanwise_cn_alpha.jl --plots`: constant-incidence equilibrium,
  local-angle cn/alpha ratios, load conservation, optional derivative modes and convergence gates,
  with an optional plot smoke check. Included in `--full`.

The physical, parallel and plotting switches are opt-in even with `--full`; they can also be passed
to `runtests.jl`. Each script runs in a separate module within that process.
Shared moving-block damping tests live in the package's `test/` directory and
are covered by its study suite, rather than being duplicated here.

For a manual structural-animation preview from an existing serialized result:

```powershell
julia --project=lib/WingPropellerUVLM lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/validation/verify_pazy_structural_animation.jl path/to/result.jls
```

This preview uses the original `output/structural_style_preview/` location.
It is excluded from automated tests because it requires a saved result.
