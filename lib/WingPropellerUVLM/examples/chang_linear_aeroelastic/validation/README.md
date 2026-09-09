# Chang validation tools

These files are focused checks, not alternative production drivers.

| File | Check |
|:--|:--|
| `verify_chang_structural_model.jl` | Mass/stiffness integrity and published modal frequencies |
| `audit_chang_virtual_work.jl` | All-coordinate virtual work at reference, uniform, and nonuniform wing rotations |
| `investigate_chang_wing_moment_transfer.jl` | Production load-map regression against geometry derivatives at three step sizes |
| `compare_imperial_legacy_loads.jl` | Corrected and retained legacy load reconstruction |
| `analyze_chang_damping.jl` | Simple exponential-envelope fit for a saved history |
| `conservative_lumped_modal/` | Evidence and regression tests for the adopted inertia remap |

The structural, virtual-work, investigation, and load-comparison scripts run
or construct a case themselves. The damping tool takes a history CSV path as
its first argument. See [the moment-transfer review](WING_MOMENT_TRANSFER_REVIEW.md)
for the Euler-axis projection and its validation.
