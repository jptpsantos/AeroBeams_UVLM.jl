# Chang linear aeroelastic validation

This example ports the latest partitioned generalized-alpha Chang driver from
the local `UVLM_Chang` research folder into `WingPropellerUVLM`. The source
folder is treated as read-only. The copied case retains its wing, nodal mass,
sectional stiffness, flexible-pylon modal, rotor, and default time-integration
parameters.

The aerodynamic side uses the current package APIs:

- `initialize_bohnisch_uvlm_system` for the wing, blades, wakes, and interaction groups;
- `snapshot_uvlm` / `restore_uvlm!` for repeatable partitioned trials;
- `propagate_system!` for the unsteady free-wake step;
- `imperial_nodal_forces` and `imperial_nodal_positions` for structural load transfer.

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

Available runtime overrides are `CHANG_END_TIME_S`,
`CHANG_IMPULSE_MAGNITUDE`, `CHANG_IMPULSE_START_S`,
`CHANG_IMPULSE_DURATION_S`, `CHANG_GA_RHO_INF`,
`CHANG_COUPLING_MAX_ITER`, `CHANG_COUPLING_TOL_U`,
`CHANG_COUPLING_TOL_F`, `CHANG_COUPLING_TOL_EQ`, and
`CHANG_COUPLING_RELAXATION`. Plotting is enabled by default; use
`CHANG_PLOT_RESULTS=false` to disable it or `CHANG_PLOT_END_TIME_S` to change
the default five-second horizontal axis.

Results are written to `output/` as a CSV history and a text validation
summary. The main input also saves `chang_linear_imperial_uvlm_history.png`.
For one propeller, it contains the wing-tip displacement and twist followed by
the propeller pitch and yaw. For two propellers, it contains the inboard pitch
and yaw followed by the outboard pitch and yaw. A successful run prints
`VALIDATION_PASSED=true`.

Generate a six-panel PNG of the structural response and coupling convergence
from the saved CSV with:

```powershell
julia lib/WingPropellerUVLM/examples/chang_linear_aeroelastic/plot_chang_time_histories.jl
```
