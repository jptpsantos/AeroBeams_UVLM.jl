# Chang linear aeroelastic validation

This example ports the latest partitioned generalized-alpha Chang driver from
the local `UVLM_Chang` research folder into `WingPropellerUVLM`. The source
folder is treated as read-only. The copied case retains its wing, nodal mass,
sectional stiffness, flexible-pylon modal, rotor, and default time-integration
parameters.

The aerodynamic side uses the current package APIs:

- `initialize_bohnisch_uvlm_system` for the wing, blades, wakes, and interaction groups;
- `snapshot_uvlm` / `restore_uvlm!` for repeatable partitioned trials;
- `propagate_system!(...; advance_wake=false)` for repeated circulation/load trials;
- `advance_wake!` for the single accepted free-wake update;
- `near_field_forces!` for the direct Imperial chord/span/unsteady loads;
- `imperial_nodal_forces` and `imperial_nodal_positions` for structural load transfer.

Select the corrected direct Imperial port with
`near_field_force_model = :imperial` in `SimulationConfig`.
`legacy_nearfield.jl` now contains only the segment-load formulation identified
in the original source as following the Imperial College C++ UVLM
implementation; the older `PanelProperties` reconstruction and inactive trial
implementations were removed. Circulation, interaction masking, free-wake
propagation, and the partitioned structural time integrator are unchanged. The
force loop reuses the system force arrays and can distribute independent panels
over Julia threads without changing the segment equations. Its chordwise
segment velocities include the receiving bound surface and exclude only the
coincident filament, correcting the broader self-surface omission in the
original legacy block. To run that alternative implementation, set
`near_field_force_model = :legacy_imperial_segments` in `SimulationConfig`.

Propeller pitch/yaw moments can be projected in two ways with
`SimulationConfig.propeller_moment_projection`. The default
`:exact_virtual_work` option projects the aerodynamic modal wrench onto the
instantaneous axes implied by `R_y(pitch) * R_z(-yaw)` and the local wing
rotation. Use `:fixed_aero_axes` to reproduce the original small-angle
`(M_y, -M_z)` projection.

The [complete library and time-marching guide](../../../../docs/src/wing-propeller-uvlm-library-guide.md)
walks through every phase of the driver, the callback, source files, and public
functions. The run file and `chang_uvlm_coupling.jl` are also commented at the
transaction and coordinate/load-transfer boundaries.

Run the full Chang input from the repository root with:

```powershell
julia --threads=auto lib/WingPropellerUVLM/examples/chang_linear_aeroelastic/run_chang_linear_aeroelastic.jl
```

## Isolated-propeller trims

Both trim drivers use `SimulationConfig.near_field_force_model` by default.
Override it with `CHANG_TRIM_NEAR_FIELD_FORCE_MODEL=imperial` for the corrected
direct formulation or `legacy_imperial_segments` for the original-compatible
segment formulation.

The windmilling trim varies RPM until the mean aerodynamic shaft torque is
zero:

```powershell
$env:CHANG_TRIM_NEAR_FIELD_FORCE_MODEL = "imperial"
julia --project=lib/WingPropellerUVLM lib/WingPropellerUVLM/examples/chang_linear_aeroelastic/run_chang_windmilling_trim.jl
```

The powered trim varies RPM until the mean propeller thrust reaches the
positive value supplied with `CHANG_TRIM_TARGET_THRUST_N`. Its default target
is 650 N. The blade angle at 75% radius is taken directly from the Chang twist
distribution; no separate beta75 input is required.

```powershell
$env:CHANG_TRIM_NEAR_FIELD_FORCE_MODEL = "imperial"
$env:CHANG_TRIM_TARGET_THRUST_N = "650.0"
julia --project=lib/WingPropellerUVLM lib/WingPropellerUVLM/examples/chang_linear_aeroelastic/run_chang_thrusting_trim.jl
```

For a deliberate collective-pitch variation, use
`CHANG_TRIM_COLLECTIVE_OFFSET_DEG`. Leave it unset (the default is zero) for
the original Chang geometry.

Each driver marches one rigid four-bladed propeller with its free wake and
averages the selected near-field loads over complete final revolutions. It
writes all RPM evaluations, the final periodic history, and a summary under
`output/chang_windmilling_trim/` or `output/chang_thrusting_trim/`. The
windmilling RPM can replace
`PropellerConfig.rotation_rpm`; `chang_model_parameters.jl` then preserves its
advance ratio when the aeroelastic flow speed changes.

The principal convergence controls are `CHANG_TRIM_RADIAL_PANELS`,
`CHANG_TRIM_CHORDWISE_PANELS`, `CHANG_TRIM_AZIMUTH_STEP_DEG`,
`CHANG_TRIM_SIMULATED_REVOLUTIONS`, and
`CHANG_TRIM_RETAINED_WAKE_REVOLUTIONS`. Trim and aeroelastic calculations must
use matching blade grid, azimuth-step, vortex-core, and wake-length settings.
The root-solver controls are `CHANG_TRIM_INITIAL_RPM`, `CHANG_TRIM_MIN_RPM`,
`CHANG_TRIM_MAX_RPM`, `CHANG_TRIM_RPM_TOLERANCE`, and
`CHANG_TRIM_MAX_ITERATIONS`. Windmilling uses
`CHANG_TRIM_TORQUE_TOLERANCE_NM`; powered trim uses
`CHANG_TRIM_THRUST_TOLERANCE_N`.

For reference, a preliminary 5×5-mesh check at 65 m/s, 1.225 kg/m³, zero
incidence, a 5° azimuth step, and three retained wake revolutions placed the
Imperial zero-torque point near 1217.69 rpm. This value is not the final trim
for the configured 20-radial × 10-chordwise mesh; rerun the trim after any
mesh, azimuth-step, vortex-core, or wake-length change.

The default input retains the original three-second simulation. For a short
coupling smoke test, override only the run duration and impulse timing:

```powershell
$env:CHANG_END_TIME_S = "0.012"
$env:CHANG_IMPULSE_START_S = "0.002"
$env:CHANG_IMPULSE_DURATION_S = "0.003"
julia lib/WingPropellerUVLM/examples/chang_linear_aeroelastic/run_chang_linear_aeroelastic.jl
```

Available runtime overrides are `CHANG_END_TIME_S`, `CHANG_OUTPUT_DIR`,
`CHANG_OUTPUT_LABEL`, `CHANG_TRIM_REVOLUTIONS`,
`CHANG_TRIM_AVERAGE_REVOLUTIONS`, `CHANG_IMPULSE_MAGNITUDE`,
`CHANG_IMPULSE_START_S`, `CHANG_IMPULSE_DURATION_S`, `CHANG_GA_RHO_INF`,
`CHANG_COUPLING_MAX_ITER`, `CHANG_COUPLING_TOL_U`,
`CHANG_COUPLING_TOL_F`, `CHANG_COUPLING_TOL_EQ`,
`CHANG_COUPLING_TOL_COUPLED_EQ`, `CHANG_COUPLING_RELAXATION`,
`CHANG_STATE_ABORT_NORM`, and `CHANG_PROP_ANGLE_ABORT_DEG`. Plotting is enabled
by default; use `CHANG_PLOT_RESULTS=false` to disable it or
`CHANG_PLOT_END_TIME_S` to change the default five-second horizontal axis.

## 80--86 m/s divergence sweep

Run the standalone sweep driver from the repository root:

```powershell
julia --project=lib/WingPropellerUVLM lib/WingPropellerUVLM/examples/chang_linear_aeroelastic/run_chang_speed_sweep.jl
```

It runs 80, 81, ..., 86 m/s in separate Julia processes and stops after the
first divergent or failed case. The complete output of each speed is saved in
its own log, while one compact result is printed to the terminal. The partial
or complete sweep table is continuously saved as
`output/speed_sweep_80_86/speed_sweep_summary.csv`.

Divergence is declared when the run terminates before its requested end time,
the absolute pitch or yaw reaches 15 degrees, or a post-impulse envelope fit
has growth rate at least 0.02/s, fitted growth ratio at least 1.05, and
R-squared at least 0.5. The principal controls are
`CHANG_SWEEP_END_TIME_S`, `CHANG_SWEEP_FIT_START_S`,
`CHANG_SWEEP_DIVERGENCE_RATE_PER_S`, `CHANG_SWEEP_ENVELOPE_RATIO`,
`CHANG_SWEEP_MIN_FIT_R2`, and `CHANG_SWEEP_ABORT_ANGLE_DEG`. The start, stop,
and increment can be changed with `CHANG_SWEEP_SPEED_START_MPS`,
`CHANG_SWEEP_SPEED_STOP_MPS`, and `CHANG_SWEEP_SPEED_STEP_MPS`.

The sweep defaults to the current 30x5 wing grid, 10x5 propeller grid, legacy
Imperial segment forces, and exact virtual-work moment projection. For a fast
smoke test, reduce the mesh and duration, for example:

```powershell
$env:CHANG_SWEEP_SPEED_STOP_MPS = "80"
$env:CHANG_SWEEP_WING_SPAN_PANELS = "4"
$env:CHANG_SWEEP_PROP_RADIAL_PANELS = "2"
$env:CHANG_SWEEP_END_TIME_S = "0.08"
$env:CHANG_SWEEP_FIT_START_S = "0.04"
julia --project=lib/WingPropellerUVLM lib/WingPropellerUVLM/examples/chang_linear_aeroelastic/run_chang_speed_sweep.jl
```

Wake animation is controlled directly near the top of
`run_chang_linear_aeroelastic.jl`:

```julia
const ANIMATE_WAKE = true
```

Set it to `false` when the additional memory and rendering time are not
desired. Frame stride and playback rate can still be overridden at runtime:

```powershell
$env:CHANG_ANIMATION_STRIDE = "5"
$env:CHANG_ANIMATION_FPS = "15"
julia lib/WingPropellerUVLM/examples/chang_linear_aeroelastic/run_chang_linear_aeroelastic.jl
```

The driver then calls `animate_chang_wing_wake` from `chang_animation.jl` and
writes `<output-label>_wing_wake.gif`. Frames contain accepted structural/UVLM
states only; rejected partitioned trials are never shown. A smaller stride
gives smoother motion but requires more storage and rendering time. The three
plot axes use one common physical scale, so one metre has the same displayed
length in `x_A`, `y_A`, and `z_A`. For the current input, the modeled wing span
is 7.50 m while the propeller diameter is 2.30 m.

Choose which propellers receive the pitch impulse in `chang_case.jl`:

```julia
const SIMULATION_CONFIG = SimulationConfig(
    impulse_propeller_indices = [1, 2],
)
```

The listed indices must exist in `PropellerConfig.attachment_eta`.

Results are written to `output/` as a CSV history and a text validation
summary. With the default force model, the main input also saves
`chang_linear_imperial_uvlm_history.png`.
When wake animation is enabled it additionally saves
`chang_linear_imperial_uvlm_wing_wake.gif`.
For one propeller, it contains the wing-tip displacement and twist followed by
the propeller pitch and yaw. For two propellers, it contains the inboard pitch
and yaw followed by the outboard pitch and yaw. A successful run prints
`VALIDATION_PASSED=true`.

Generate a six-panel PNG of the structural response and coupling convergence
from the saved CSV with:

```powershell
julia lib/WingPropellerUVLM/examples/chang_linear_aeroelastic/plot_chang_time_histories.jl
```
