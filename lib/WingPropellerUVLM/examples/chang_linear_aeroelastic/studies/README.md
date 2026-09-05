# Chang parameter studies

These entry points reuse the production model without changing its source.
Each writes generated data below the example-level `output/` directory.

| Entry point | Purpose |
|:--|:--|
| `run_chang_windmilling_trim.jl` | Find the isolated-propeller zero-torque RPM |
| `run_chang_thrusting_trim.jl` | Find the RPM for a requested positive thrust |
| `run_chang_speed_sweep.jl` | March consecutive airspeeds and stop at divergence/failure |
| [`convergence/`](convergence/README.md) | Aerodynamic-first UVLM and gated damping convergence |

Run all files from the repository root with the local package project. For
example:

```powershell
julia --project=lib/WingPropellerUVLM lib/WingPropellerUVLM/examples/chang_linear_aeroelastic/studies/run_chang_speed_sweep.jl
```

Environment variables used only by a study have a task prefix such as
`CHANG_TRIM_*`, `CHANG_SWEEP_*`, `CHANG_AERO_SWEEP_*`, or `CHANG_AE_*`.
The child aeroelastic cases receive the ordinary production `CHANG_*`
configuration variables.
