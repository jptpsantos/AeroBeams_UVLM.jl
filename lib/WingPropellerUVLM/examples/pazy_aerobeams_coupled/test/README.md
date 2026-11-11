# Pazy verification

From the repository root, with the AeroBeams and WingPropellerUVLM environments
already instantiated:

```powershell
julia --project=lib/WingPropellerUVLM lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/test/runtests.jl
julia --project=lib/WingPropellerUVLM lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/test/runtests.jl --full
```

The default suite checks coordinate and load transfer, virtual work, Newton
load updates, animation output paths, stability statistics, spectrograms,
output routing, bifurcation analysis and plots. It includes a small coupled
bifurcation smoke test. `--full` also runs the startup/pulse coupling check,
static equilibrium, independent and continuous LCO sweeps, and diagnostics.
These are implementation regressions, not flutter convergence studies.

Every test file can also run directly. Existing optional switches remain:

- `test_coupling.jl --dynamic`: coupled startup/pulse and final animation frame.
- `test_pazy_stepped_bifurcation.jl --physical`: stepped physical continuation.
- `test_pazy_bifurcation.jl --parallel`: distributed sweep consistency.
- `test_pazy_stability_sweep.jl --plots`: stability plot smoke check.

The last three options are opt-in even with `--full`; they can also be passed
to `runtests.jl`. Each script runs in a separate module within that process.
Shared moving-block damping tests live in the package's `test/` directory and
are covered by its study suite, rather than being duplicated here.

For a manual structural-animation preview from an existing serialized result:

```powershell
julia --project=lib/WingPropellerUVLM lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/validation/verify_pazy_structural_animation.jl path/to/result.jls
```

This preview uses the original `output/structural_style_preview/` location.
It is excluded from automated tests because it requires a saved result.
