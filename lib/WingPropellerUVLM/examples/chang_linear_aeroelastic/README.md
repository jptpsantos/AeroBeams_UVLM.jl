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
- `imperial_nodal_forces` and `imperial_nodal_positions` for structural load transfer.

The [complete library and time-marching guide](../../../../docs/src/wing-propeller-uvlm-library-guide.md)
walks through every phase of the driver, the callback, source files, and public
functions. The run file and `chang_uvlm_coupling.jl` are also commented at the
transaction and coordinate/load-transfer boundaries.

Run the full Chang input from the repository root with:

```powershell
julia lib/WingPropellerUVLM/examples/chang_linear_aeroelastic/run_chang_linear_aeroelastic.jl
```

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
`CHANG_COUPLING_TOL_COUPLED_EQ`, and `CHANG_COUPLING_RELAXATION`. Plotting is
enabled by default; use `CHANG_PLOT_RESULTS=false` to disable it or
`CHANG_PLOT_END_TIME_S` to change the default five-second horizontal axis.

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
summary. The main input also saves `chang_linear_imperial_uvlm_history.png`.
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
