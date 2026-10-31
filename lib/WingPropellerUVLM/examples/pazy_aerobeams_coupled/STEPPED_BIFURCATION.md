# Pazy bifurcation with state-preserving speed steps

This separate analysis follows the supplied propeller procedure: hold speed,
discard the transient, record late extrema, move to the next speed with the
accepted physical state, and reverse the speed sequence. It uses the existing
nonlinear AeroBeams/UVLM coupling without changing existing files.

From the repository root:

```powershell
# Validate the configuration; no simulations or output files:
julia --project=. lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/run_pazy_stepped_bifurcation.jl --check-config

# Run the configured study:
julia --project=. lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/run_pazy_stepped_bifurcation.jl
```

Configure the NEW `run_pazy_stepped_bifurcation.jl` file. Its defaults are 5 degrees,
40:1:50 m/s, 10 seconds per plateau, and extrema from the final two seconds.
There are 22 plateaus and 220 seconds of simulated time. The reversal speed is
held twice, once per direction, just as the propeller script solves its highest
speed at the end of the upsweep and again at the start of the downsweep.
`:down_up`, `:upsweep` and `:downsweep` are also available. The first speed of the
chosen sequence determines the single static initialization.

To match the propeller script's observation durations, set `plateau_time=120.0`
and `transient_time=115.0`. This can be expensive for UVLM; select durations using
settling and timestep sensitivity, rather than assuming the same seconds suit
both models. The existing solver retains scalar arrays for the entire single
integration, so memory use increases with total duration and number of speeds.

Programmatic use, without executing an entry point:

```julia
include("lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/PazySteppedBifurcation.jl")
using .PazySteppedBifurcation
result = run_pazy_stepped_bifurcation(
    U_min=40., U_max=50., dU=0.5,
    plateau_time=20., transient_time=15.,
    angle_of_attack=deg2rad(5.),
    sweep_order=:up_down,
    on_unconverged=:warn,
    save_time_histories=true)
println(result.output_directory)
```

Only a single initial tip pulse is applied (default coefficient 0.1 N, starting
at 0.1 s for 0.04 s). No pulse, structural reset, equilibrium solve, or wake reset
is introduced at later speeds or at reversal. The pulse must finish before the
first observation window and span at least eight timesteps. Set its coefficient
to zero for an unforced study. Root incidence is in radians.

The integration uses a common physical timestep sized at the maximum speed and
shortened to divide each plateau exactly. An explicit `fixed_time_step` must
divide `plateau_time` and respect the requested maximum `U*dt/chord`. Wake capacity
is sized at the minimum speed. A boundary endpoint belongs to the plateau just
completed; the next accepted step uses the new speed. The shared boundary sample
is excluded from the next plateau's observation window and optional history.

`extrema_method=:quadratic` fits three neighboring samples of each response and
solves the local fitted derivative equal to zero. This refines the peak time and
value without treating, for example, nodal velocity as strain rate. It is an
approximation from scalar samples, not the propeller ODE's dense solution.
`:sampled` preserves the original discrete turning-point convention. Neither
method fits across a speed change or observation-window boundary. Flat extrema
retain the last sample of the plateau. Summary envelopes include refined peaks
and raw window extrema; `sample_minimum` and `sample_maximum` preserve the latter.

Convergence uses the existing final-window diagnostic on raw samples. Every
status is retained. `on_unconverged=:warn` reports unverified responses and moves
on; `:continue` records them silently; `:error` saves the completed plateau and
stops. Holds are fixed in advance: there is no automatic extension. Increase
the hold/window if a periodic response has not settled. Modulated or nonperiodic
responses can also fail the cycle-stationarity diagnostic.

Outputs go into a fresh directory under `output/bifurcation_stepped`:

- `plateau_schedule.csv` and `settings.txt` describe the complete speed order.
- `bifurcation_summary.csv` and `bifurcation_local_extrema.csv` contain all
  completed plateaus, including direction and global/local time information.
- `upsweep/` and `downsweep/` contain separate compatible summary/peak CSVs.
- `plateau_NNNN_extrema.csv`, `_local_extrema.csv`, and optional `_history.csv`
  preserve each completed plateau. Aggregate CSVs update after every plateau.
- `plateau_boundaries.csv` records uninterrupted aerodynamic time/step counters,
  state norms and problem identities. Optional `save_boundary_states=true`
  returns physical snapshots in memory for diagnostics; these are not a disk
  restart implementation. `keep_solver_result=true` returns the full solver result.
- `failure.txt`, if present, explains a failure and lists the completed plateaus.
- `figures/plots_*/` contains PNG/PDF bifurcation figures and interpolated curves.
  AeroBeams incidence colors are retained; solid/open marks identify the upsweep,
  dashed/filled marks identify the downsweep, and crosses mark unverified samples.

CSV-only replotting, always into a new figure directory:

```powershell
julia --project=. lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/PazySteppedBifurcationPlots.jl path/to/stepped_run
```

Each direction is interpolated separately. Missing static equilibria remain NaN;
the oscillatory mean is never used as a static reference. The fourth panel is
geometric incidence (root incidence plus twist), not effective aerodynamic AoA.

A descending sweep should start on an established LCO if the goal is to follow
that attractor. A peak speed beyond an oscillatory region may let the motion
decay before reversal. Compare disturbance amplitudes, hold durations and speed
increments before interpreting a directional difference as settled hysteresis.
These time integrations sample attracting responses; they do not compute
unstable periodic branches or certify Hopf/fold locations.

Verification:

```powershell
julia --startup-file=no --compiled-modules=existing --project=. lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/test_pazy_stepped_bifurcation.jl --physical
```

Without `--physical`, only schedule, derivative-zero extraction and synthetic
statistics checks run. The optional millisecond physical tests check coupling,
state continuity, output preservation and plotting; they are not settled LCO studies.
