# Chang validation tools

These files are focused checks, not alternative production drivers.

| File | Check |
|:--|:--|
| `verify_chang_structural_model.jl` | Mass/stiffness integrity and published modal frequencies |
| `compare_chang_pylon_frequencies.jl` | Chang RCAS Tables 3/4/7 versus the retuned rigid and assumed-mode pylons; current defaults and a 30-element case at the published 6.25 m attachment |
| `audit_chang_virtual_work.jl` | All-coordinate virtual work at reference, uniform, and nonuniform wing rotations |
| `investigate_chang_wing_moment_transfer.jl` | Production load-map regression against geometry derivatives at three step sizes |
| `compare_imperial_legacy_loads.jl` | Corrected and retained legacy load reconstruction |
| `analyze_chang_damping.jl` | Simple exponential-envelope fit for a saved history |
| `conservative_lumped_modal/` | Evidence and regression tests for the adopted inertia remap |

The structural, virtual-work, investigation, and load-comparison scripts run
or construct a case themselves. The damping tool takes a history CSV path as
its first argument. See [the moment-transfer review](WING_MOMENT_TRANSFER_REVIEW.md)
for the Euler-axis projection and its validation.

Run the pylon frequency comparison from the repository root with the studies
environment, which provides plotting dependencies:

```powershell
julia --project=lib/WingPropellerUVLM/studies lib/WingPropellerUVLM/examples/chang_linear_aeroelastic/validation/compare_chang_pylon_frequencies.jl
```

Tables, full nonspinning spectra, model settings, and PNG/SVG comparisons are
written to `output/validation/chang_pylon_comparison/`. Missing pylon torsion and
higher pylon bending modes are marked as not retained; isolated-wing reference
frequencies are compared only to the bare-wing model.
