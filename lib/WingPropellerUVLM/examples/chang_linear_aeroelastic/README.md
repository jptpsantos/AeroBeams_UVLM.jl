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

The runner loads one complete `config`, calls `run_chang(config)`, and exposes
`chang_run` and `results` for inspection. All case settings, environment
overrides, and automatic values are resolved before model construction.

| Object | Contents |
|:--|:--|
| `chang_run.config` | Active geometry, flow, finite-core, wake, excitation, integration, coupling, output, and structural settings |
| `chang_run.model` | Configuration, derived parameters, structural matrices, and aerodynamic load model |
| `chang_run.workspace` | This run's mutable UVLM system, wake, geometry, and load buffers |
| `chang_run.solution` | Accepted structural histories and coupling diagnostics |
| `chang_run.results` | Extracted responses and output paths |

`src/ChangAeroelastic.jl` defines the example module and the build/solve/save
workflow. `src/chang_workspaces.jl` provides `build_chang_workspace(model)`,
which allocates fresh aerodynamic storage. The geometry and load-transfer
functions receive `model` and `workspace` explicitly. Loading the module
defines functions; it does not construct or run a case.

Within each time step, the solver converges motion and loads with the wake
fixed, then advances the accepted wake once. Numerical integration remains in
`src/chang_simulation.jl`. Tabulated Chang mass and stiffness distributions
remain reference data in `src/chang_model_parameters.jl`.

### Use the API or run multiple cases

With the WingPropellerUVLM project active, load the module once and build each
case from fresh defaults:

```julia
include("lib/WingPropellerUVLM/examples/chang_linear_aeroelastic/src/ChangAeroelastic.jl")
using .ChangAeroelastic

defaults = chang_case_defaults()
defaults = merge(defaults, (
    simulation = merge(defaults.simulation, (; freestream_speed_mps = 80.0)),
    coupling = merge(defaults.coupling, (; state_tolerance = 1e-6)),
))
config = load_chang_configuration(defaults; env = Dict{String,String}())
response = run_chang(config)
```

Omit `env` to apply the existing `CHANG_*` environment overrides. Resolve
settings again after changing inputs so automatic wake counts and impulse
timing use the new case. Environment changes after loading do not affect the
resolved configuration. Each run owns a copy of its settings and fresh
aerodynamic storage, so successive cases do not overwrite one another.

For an audit that only needs the model and geometry, use
`model = build_chang_model(config)` and
`workspace = build_chang_workspace(model)`. Then call, for example,
`update_aero_geometry_for_state!(model, workspace, state, time)`.
Treat model data as read-only after construction; build a new model when
changing the case. The primary runner reloads `chang_case.jl` on each
invocation, including when rerun in the same IDE session.

## Run the primary case

From the repository root:

```powershell
julia --threads=auto lib/WingPropellerUVLM/examples/chang_linear_aeroelastic/run_chang_linear_aeroelastic.jl
```

The default case uses a 20-by-10 wing grid, a 7-by-7 grid per propeller blade,
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
| `wing` | Wing dimensions and mesh |
| `propeller` | Blade geometry, mesh, installation, RPM, and pitch offset |
| `simulation` | Air density, speed, flow angles, duration, azimuth step, interaction, and load models |
| `aerodynamic` | Finite-core factors, elastic axis, and modal load arm |
| `wake` | Retained wake row counts |
| `excitation` | Trim revolutions, averaging window, and pitch impulse |
| `integration` | Generalized-alpha damping and stop limits |
| `coupling` | Iteration limit, convergence tolerances, and relaxation |
| `output` | Directory, filename label, plots, and animation |
| `structural` | Damping, pylon properties, and rotor mass/inertia inputs |

The file returns these groups through `chang_case_defaults()`. The runner uses
`load_chang_configuration()` to turn them into the complete active `config`.

Environment parsing, automatic values, and input checks live in
`src/chang_configuration.jl`. Automated studies can override the visible
defaults through the existing environment variables; a set environment value
takes precedence over the corresponding setting in the file:

| Variable | Case default | Meaning |
|:--|--:|:--|
| `CHANG_WING_SPAN_PANELS` | `20` | Wing spanwise panels and beam elements |
| `CHANG_WING_CHORD_PANELS` | `10` | Wing chordwise panels |
| `CHANG_PROP_RADIAL_PANELS` | `7` | Radial panels per blade |
| `CHANG_PROP_CHORD_PANELS` | `7` | Chordwise panels per blade |
| `CHANG_ROTATION_RPM` | `1212.0` | RPM at the trim reference speed |
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

For example, set `interaction_on = true` in `simulation` to include
aerodynamic influence between the wing and propellers (and between separate
propellers). With `false`, those groups are aerodynamically isolated, while
blades within each propeller still interact. Structural wing–propeller
coupling remains active in both cases. If `CHANG_INTERACTION` is set in the
environment, it takes precedence over this file setting.

For example, the following edits can all be made in `chang_case.jl`:

```julia
# Inside the corresponding group:
segment_core_factor = 0.01,             # aerodynamic
maximum_rows_wing = 120,                # wake
maximum_rows_propeller = 144,           # wake
rho_inf = 0.8,                          # integration
state_tolerance = 1.0e-6,               # coupling
maximum_iterations = 20,                # coupling
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

Run `validation/verify_chang_context.jl` to check configuration isolation,
automatic settings, and interleaved aerodynamic trials with different models
in one Julia session.

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
