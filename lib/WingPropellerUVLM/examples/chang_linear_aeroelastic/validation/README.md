# Chang validation tools

These files are focused checks, not alternative production drivers.

| File | Check |
|:--|:--|
| `verify_chang_structural_model.jl` | Mass/stiffness integrity and published modal frequencies |
| `audit_chang_virtual_work.jl` | Work-conjugacy of aerodynamic load transfer |
| `compare_imperial_legacy_loads.jl` | Corrected and retained legacy load reconstruction |
| `analyze_chang_damping.jl` | Simple exponential-envelope fit for a saved history |
| `conservative_lumped_modal/` | Evidence and regression tests for the adopted inertia remap |

The first three scripts run or construct a case themselves. The damping tool
takes a history CSV path as its first argument.
