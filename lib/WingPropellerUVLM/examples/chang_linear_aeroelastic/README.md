# Chang linear aeroelastic example

This example couples the Chang linear wing--propeller structural model to the
free-wake UVLM solver. The source is arranged so the normal workflow starts
from two files: edit `chang_case.jl`, then run `run_chang_linear_aeroelastic.jl`.

## Directory map

| Path | Purpose |
|:--|:--|
| `chang_case.jl` | Editable geometry, mesh, operating point, and model choices |
| `run_chang_linear_aeroelastic.jl` | Primary partitioned aeroelastic simulation |
| `src/` | Structural assembly, UVLM coupling/integration, trim, plotting, and output helpers |
| `studies/` | Trim and airspeed study entry points |
| `studies/convergence/` | Rigid-aerodynamic and gated aeroelastic convergence workflow |
| `validation/` | Focused structural, load-transfer, and damping checks |
| `output/` | Generated results; everything except `.gitignore` is ignored |

The [library and coupling guide](../../../../docs/src/wing-propeller-uvlm-library-guide.md)
documents the coordinate systems, accepted-state UVLM transaction, structural
load transfer, and generalized-alpha coupling in detail.

The primary runner follows the AeroBeams example sequence: problem setup,
problem solution, and post-processing. Step-level iteration and wake
bookkeeping are isolated in `src/chang_simulation.jl` so the entry point stays
focused on the physical workflow.

## Run the primary case

From the repository root:

```powershell
julia --threads=auto lib/WingPropellerUVLM/examples/chang_linear_aeroelastic/run_chang_linear_aeroelastic.jl
```

The default case uses a 20-by-5 wing grid, a 5-by-5 grid per propeller blade,
85 m/s, and the corrected Imperial near-field loads. It writes a CSV history,
a text summary, and (by default) a PNG under `output/`. Wake animation is off
by default because retained geometry and GIF files can be large.

For a quick path and coupling check:

```powershell
$env:CHANG_END_TIME_S = "0.001"
$env:CHANG_WING_SPAN_PANELS = "4"
$env:CHANG_WING_CHORD_PANELS = "2"
$env:CHANG_PROP_RADIAL_PANELS = "2"
$env:CHANG_PROP_CHORD_PANELS = "2"
$env:CHANG_PLOT_RESULTS = "false"
$env:CHANG_ANIMATE_WAKE = "false"
julia --threads=1 lib/WingPropellerUVLM/examples/chang_linear_aeroelastic/run_chang_linear_aeroelastic.jl
```

## Configuration

Routine inputs live together in `chang_case.jl`. Automated studies can
override them without editing source:

| Variable | Case default | Meaning |
|:--|--:|:--|
| `CHANG_WING_SPAN_PANELS` | `20` | Wing spanwise panels and beam elements |
| `CHANG_WING_CHORD_PANELS` | `5` | Wing chordwise panels |
| `CHANG_PROP_RADIAL_PANELS` | `5` | Radial panels per blade |
| `CHANG_PROP_CHORD_PANELS` | `5` | Chordwise panels per blade |
| `CHANG_ROTATION_RPM` | `1217.6962` | RPM at the trim reference speed |
| `CHANG_TRIM_SPEED_MPS` | `65` | RPM reference speed |
| `CHANG_SPEED_MPS` | `85` | Aeroelastic flow speed |
| `CHANG_AOA_DEG` | `0` | Angle of attack |
| `CHANG_AZIMUTH_STEP_DEG` | `5` | Rotor increment and aerodynamic time step |
| `CHANG_END_TIME_S` | `5` | Requested duration |
| `CHANG_INTERACTION` | `false` | Wing--propeller aerodynamic interaction |
| `CHANG_NEAR_FIELD_FORCE_MODEL` | `imperial` | `imperial` or `legacy_imperial_segments` |
| `CHANG_PROP_MOMENT_PROJECTION` | `exact_virtual_work` | Exact or fixed-axis modal projection |
| `CHANG_PLOT_RESULTS` | `true` | Write/display the response PNG |
| `CHANG_ANIMATE_WAKE` | `false` | Record accepted states and write a GIF |
| `CHANG_OUTPUT_DIR` | `output/` | Result directory |

Finite-core, wake, impulse, generalized-alpha, coupling, abort, and animation
controls are collected by the named option helpers in `src/chang_simulation.jl`;
they all begin with `CHANG_`.

## Studies

| Task | Entry point |
|:--|:--|
| Zero-torque windmilling trim | `studies/run_chang_windmilling_trim.jl` |
| Powered thrust trim | `studies/run_chang_thrusting_trim.jl` |
| Airspeed stability sweep | `studies/run_chang_speed_sweep.jl` |
| UVLM and damping convergence | [`studies/convergence/README.md`](studies/convergence/README.md) |

For example:

```powershell
julia --project=lib/WingPropellerUVLM lib/WingPropellerUVLM/examples/chang_linear_aeroelastic/studies/run_chang_windmilling_trim.jl
```

The recommended convergence process is aerodynamic-first: establish rigid
wing/propeller `CL` and `CT` convergence, then run only the gated aeroelastic
damping cases selected from that result.

## Validation

Run the structural modal and matrix checks with:

```powershell
julia --project=lib/WingPropellerUVLM lib/WingPropellerUVLM/examples/chang_linear_aeroelastic/validation/verify_chang_structural_model.jl
```

`validation/audit_chang_virtual_work.jl` checks work-conjugate load transfer,
while `validation/conservative_lumped_modal/` retains the compact evidence for
the spatial-inertia remap used by the production model.

`VALIDATION_PASSED=true` means a requested time march completed with finite
states and converged partitioned steps. It is a runtime integrity check, not by
itself proof of grid convergence or agreement with a flutter boundary.
