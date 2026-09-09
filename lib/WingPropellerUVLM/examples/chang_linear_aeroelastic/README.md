# Chang linear aeroelastic example

This example couples the Chang linear wing--propeller structural model to the
free-wake UVLM solver. The source is arranged so the normal workflow starts
from two files: edit `chang_case.jl`, then run `run_chang_linear_aeroelastic.jl`.

## Directory map

| Path | Purpose |
|:--|:--|
| `chang_case.jl` | Editable physical inputs, numerical controls, and output settings |
| `run_chang_linear_aeroelastic.jl` | Primary partitioned aeroelastic simulation |
| `src/` | Structural assembly, UVLM coupling/integration, trim, plotting, and output helpers |
| `studies/` | Trim and airspeed study entry points |
| `studies/convergence/` | Rigid-aerodynamic and gated aeroelastic convergence workflow |
| `validation/` | Focused structural, load-transfer, and damping checks |
| `output/` | Generated results; everything except `.gitignore` is ignored |

The [library and coupling guide](../../../../docs/src/wing-propeller-uvlm-library-guide.md)
documents the coordinate systems, accepted-state UVLM transaction, structural
load transfer, and generalized-alpha coupling in detail.

## How to read the main script

`run_chang_linear_aeroelastic.jl` follows five numbered steps:

| Step | What it does | Where the details live |
|:--|:--|:--|
| 1. Case and output | Read physical inputs, select load models, enable requested plots | `chang_case.jl`, `src/chang_postprocessing.jl` |
| 2. Structure | Build and check the mass, damping/gyroscopic, and stiffness matrices | `src/chang_model_parameters.jl`, `src/chang_structural_model.jl` |
| 3. Aerodynamics | Initialize surfaces, finite cores, and wake storage | `src/chang_uvlm_coupling.jl` |
| 4. Simulation | Configure trim, impulse, and coupling; advance the response | `src/chang_simulation.jl` |
| 5. Results | Write CSV/summary files and the requested PNG/GIF | `src/chang_postprocessing.jl` |

The main objects follow the same flow: `structural` holds the reduced matrices,
`uvlm` holds aerodynamic storage, `solution` holds the state histories and
coupling diagnostics, and `results` holds extracted responses and file paths.
Within each time step, the solver converges motion and loads with the wake
fixed, then advances the accepted wake once.

`src/chang_workspaces.jl` binds the initialized arrays to the names used by the
geometry/load adapter and validation scripts. Include it after creating
`structural`, `aerodynamic_options`, and `uvlm`. These bindings share the same
arrays; the file does not construct another aerodynamic model.

For routine work, edit `chang_case.jl`. Its groups include core factors, wake
retention, trim/impulse timing, integration, coupling tolerances, output, and
structural damping. The option functions in `src/chang_simulation.jl` read those
defaults and apply environment overrides. Tabulated Chang wing mass/stiffness
distributions remain reference data in `src/chang_model_parameters.jl`.

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

Edit the plain values in the numbered groups in `chang_case.jl`:

| Group | What to adjust |
|:--|:--|
| `WING_DEFAULTS` | Wing dimensions and mesh |
| `PROPELLER_DEFAULTS` | Blade geometry, mesh, installation, RPM, and pitch offset |
| `SIMULATION_DEFAULTS` | Air density, speed, flow angles, duration, azimuth step, interaction, and load models |
| `AERODYNAMIC_DEFAULTS` | Finite-core factors, elastic axis, and modal load arm |
| `WAKE_DEFAULTS` | Retained wake row counts |
| `EXCITATION_DEFAULTS` | Trim revolutions, averaging window, and pitch impulse |
| `INTEGRATION_DEFAULTS` | Generalized-alpha damping and stop limits |
| `COUPLING_DEFAULTS` | Iteration limit, convergence tolerances, and relaxation |
| `OUTPUT_DEFAULTS` | Directory, filename label, plots, and animation |
| `STRUCTURAL_DEFAULTS` | Damping, pylon properties, and rotor mass/inertia inputs |

The existing `WING_CONFIG`, `PROPELLER_CONFIG`, and `SIMULATION_CONFIG` objects
are constructed at the bottom of the file. The runner and its helpers resolve
the other groups when preparing the model and runtime options.

Environment parsing, configuration types, and input checks live in
`src/chang_configuration.jl`. Automated studies can override the visible
defaults through the existing environment variables; a set environment value
takes precedence over the corresponding setting in the file:

| Variable | Case default | Meaning |
|:--|--:|:--|
| `CHANG_WING_SPAN_PANELS` | `20` | Wing spanwise panels and beam elements |
| `CHANG_WING_CHORD_PANELS` | `5` | Wing chordwise panels |
| `CHANG_PROP_RADIAL_PANELS` | `5` | Radial panels per blade |
| `CHANG_PROP_CHORD_PANELS` | `5` | Chordwise panels per blade |
| `CHANG_ROTATION_RPM` | `1217.6962` | RPM at the trim reference speed |
| `CHANG_TRIM_SPEED_MPS` | `65` | RPM reference speed |
| `CHANG_SPEED_MPS` | `85` | Aeroelastic flow speed |
| `CHANG_AOA_DEG` | `3` | Angle of attack |
| `CHANG_AZIMUTH_STEP_DEG` | `5` | Rotor increment and aerodynamic time step |
| `CHANG_END_TIME_S` | `5` | Requested duration |
| `CHANG_INTERACTION` | `false` | Wing--propeller aerodynamic interaction |
| `CHANG_NEAR_FIELD_FORCE_MODEL` | `imperial` | `imperial` or `legacy_imperial_segments` |
| `CHANG_PROP_MOMENT_PROJECTION` | `exact_virtual_work` | Exact or fixed-axis modal projection |
| `CHANG_PLOT_RESULTS` | `true` | Write/display the response PNG |
| `CHANG_ANIMATE_WAKE` | `false` | Record accepted states and write a GIF |
| `CHANG_OUTPUT_DIR` | `output/` | Result directory |

For example, set `interaction_on = true` in `SIMULATION_DEFAULTS` to include
aerodynamic influence between the wing and propellers (and between separate
propellers). With `false`, those groups are aerodynamically isolated, while
blades within each propeller still interact. Structural wing–propeller
coupling remains active in both cases. If `CHANG_INTERACTION` is set in the
environment, it takes precedence over this file setting.

For example, the following edits can all be made in `chang_case.jl`:

```julia
# Inside the corresponding *_DEFAULTS group:
segment_core_factor = 0.01,             # AERODYNAMIC_DEFAULTS
maximum_rows_wing = 120,                # WAKE_DEFAULTS
maximum_rows_propeller = 144,           # WAKE_DEFAULTS
rho_inf = 0.8,                          # INTEGRATION_DEFAULTS
state_tolerance = 1.0e-6,               # COUPLING_DEFAULTS
maximum_iterations = 20,                # COUPLING_DEFAULTS
```

`maximum_rows_wing = nothing` retains automatic sizing: the configured
`wing_rows_per_chord_panel` times the **active** wing chordwise panel count,
including mesh overrides. Set an integer to choose the row count explicitly.
Wake length is controlled by retained age (approximately row count times time
step); its spatial extent evolves with convection. Changing the azimuth step
also changes retained wake age at a fixed row count.

`impulse_start_s = nothing` starts the pulse after `trim_revolutions` at the
active rotor speed. Set a time in seconds for an explicit start. Output
`label = nothing` derives the filename label from the selected force model.

Existing numerical overrides such as `CHANG_FCORE_SEGMENT_FACTOR`,
`CHANG_WAKE_ROWS_WING`, `CHANG_GA_RHO_INF`, and `CHANG_COUPLING_TOL_U` still take
precedence over file settings. Independent trim/convergence study grids and
validation-specific inputs remain in their study scripts; the groups above
configure the primary aeroelastic runner.

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

Wing moments use the instantaneous axes of the implemented Euler rotation
order. Propeller wrenches use the interpolated attachment angles before their
generalized moments are distributed to wing nodes. `CHANG_PROP_MOMENT_PROJECTION`
selects only the propeller pitch/yaw modal projection.

`validation/audit_chang_virtual_work.jl` asserts work-conjugate load transfer for
all wing coordinates at the reference, uniformly rotated, and nonuniformly
deformed states. It also checks all propeller modal coordinates when their
projection is `exact_virtual_work`. The extended
`validation/investigate_chang_wing_moment_transfer.jl` checks eleven deformation
scenarios at three finite-difference step sizes; the
[moment-transfer review](validation/WING_MOMENT_TRANSFER_REVIEW.md) gives the
derivation. `validation/conservative_lumped_modal/` retains the evidence for
the spatial-inertia remap used by the production model.

`VALIDATION_PASSED=true` means a requested time march completed with finite
states and converged partitioned steps. It is a runtime integrity check, not by
itself proof of grid convergence or agreement with a flutter boundary.
