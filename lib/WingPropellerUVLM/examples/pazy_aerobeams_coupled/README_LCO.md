# Pazy LCO sweeps with UVLM

For independent **constant-speed** simulations and final-window minimum/maximum
bifurcation diagrams, use `run_pazy_bifurcation.jl`; see [BIFURCATION.md](BIFURCATION.md)
for configuration, CSV replotting, convergence checks, and optional process parallelism.

Run `run_pazy_LCO.jl` for the independent upsweep and downsweep analyses,
directly from the IDE or the repository root:

```powershell
julia --project=. lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/run_pazy_LCO.jl
```

## Continuous test (separate runner)

Use `run_pazy_continuous_LCO.jl` for the single continuous simulation:

```powershell
julia --project=. lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/run_pazy_continuous_LCO.jl
```

This separate runner reverses at 52 m/s. One dynamic integration runs 35 -> 52 m/s in 40 s and
52 -> 35 m/s in 20 s, at 5 degrees root incidence. It starts from one coupled
static equilibrium at 35 m/s. There is no second initialization, hold, state
reset, or tip pulse at reversal. Structural displacement and velocity,
circulation, aerodynamic clocks, and the evolving wake remain in the same run.

The continuous runner explicitly uses `comparison_time_step=0.00041207376120325537` s and
35 wake rows to match the latest completed independent experiment. Set this
value to `nothing` for automatic dt selection. Matching dt is important because
the earlier diagnostics showed strong timestep sensitivity.

Outputs are saved separately to `output/LCO_continuous_sweeps`. The completed
upsweep CSV, partial global history, strain-versus-speed figure, and
`pazy_linear_LCO_turnaround.csv` are saved at reversal while integration
continues. The full history, downsweep CSV, and combined figures are saved when
the simulation finishes. Branch CSVs share exactly one peak-speed sample;
`time_s` is global (0-40 s up, 40-60 s down), and `branch_time_s` starts at zero
for each branch. A single global time/airspeed/strain plot shows the reversal.

Set `PAZY_CONTINUOUS_LCO_RUN=false` before including the continuous runner to
load its functions without executing a simulation. The independent runner uses
its own `PAZY_LCO_RUN` switch.

To regenerate the continuous strain-versus-speed figure from its saved CSVs,
run `plot_pazy_continuous_LCO.jl` directly, or use:

```julia
include("lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/plot_pazy_continuous_LCO.jl")
plot_pazy_continuous_LCO()
```

For a matched independent reset at 52 m/s, load the independent runner with
execution disabled (as shown below), then call:

```julia
run_pazy_linear_LCO(; initial_speed=35.0, peak_speed=52.0,
    downsweep_initial_speed=52.0, upsweep_duration=40.0, downsweep_duration=20.0,
    fixed_time_step=0.00041207376120325537, angle_of_attack=deg2rad(5.0),
    output_directory=joinpath("lib", "WingPropellerUVLM", "examples",
        "pazy_aerobeams_coupled", "output", "LCO_independent_sweeps_matched"))
```

This comparison uses the continuous test's speed ranges, ramp rates, dt, wake
capacity, mesh, and convergence controls; its downsweep starts from its own
static equilibrium. The restored independent runner uses the same 35-52 m/s
range, with automatic timestep selection, and writes to
`output/LCO_independent_sweeps`.

A short physical regression compares the continuous upsweep against the
independent upsweep and checks the shared boundary sample, live velocity,
continuing aerodynamic clock/wake counters, fixed dt, and CSV outputs:

```powershell
julia --startup-file=no --compiled-modules=existing --project=. lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/test_pazy_continuous_LCO.jl
```

## Independent runs (original runner)

`run_pazy_LCO.jl` performs two independent simulations at 5 degrees root incidence:

| Branch | Linear speed ramp | Dynamic duration | Initial condition |
| --- | --- | --- | --- |
| Upsweep | 35 to 52 m/s (+0.425 m/s^2) | 40 s | Its own coupled static equilibrium at 35 m/s |
| Downsweep | 52 to 35 m/s (-0.85 m/s^2) | 20 s | Its own coupled static equilibrium at 52 m/s |

Both branches cover the same 35-52 m/s range. The downsweep starts at
`peak_speed` by default. Edit the speed and duration arguments at the bottom of
`run_pazy_LCO.jl` to change the ramps. The optional `downsweep_initial_speed`
argument allows a different downsweep starting speed when explicitly supplied.

Each simulation has a local time axis starting at zero. No structural state,
circulation, or wake passes from the upsweep to the downsweep. There is no
preparation hold. Independent equilibrium initialization is the chosen numerical
initial condition; the quoted paper passage specifies the ramp rates and durations
but does not fully specify its initial states.

Before each ramp, the existing static UVLM/beam solver iterates nonlinear beam
balance and develops a free wake on each fixed trial geometry. Both wake settling
and coupled force/geometry convergence must pass. An unconverged equilibrium
stops the branch before time marching and leaves its convergence history/status
on disk. A fixed waiting time is not substituted for equilibrium.

The dynamic problem copies the converged beam model and full state vector, with
zero structural velocity, and seeds UVLM with that branch's developed wake,
bound circulation, previous geometry and aerodynamic history. Both clocks are
reset to zero. Equilibrium and dynamic timesteps, wake capacities, operating
points, meshes and UVLM settings are checked for compatibility. The static
initialization is preserved separately for inspection.

Edit the parameter block at the bottom of `run_pazy_LCO.jl`, or load without
execution and call the driver:

```julia
ENV["PAZY_LCO_RUN"] = "false"
include("lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/run_pazy_LCO.jl")
delete!(ENV, "PAZY_LCO_RUN")

lco_results = run_pazy_linear_LCO(;
    initial_speed=35.0, peak_speed=52.0,
    upsweep_duration=40.0, downsweep_duration=20.0,
    angle_of_attack=deg2rad(5.0))

up = lco_results.upsweep
down = lco_results.downsweep
# up.result.time and down.result.time both start at zero.
# Each branch has .equilibrium, .result, .root_strain_ue and .wingtip_aoa_degrees.
```

Run the separate tip-force case directly from the IDE or repository root:

```powershell
julia --project=. lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/run_pazy_LCO_downsweep_pulse.jl
```

With `pulse_branches=:downsweep`, this runner performs the 35 -> 52 m/s upsweep in 40 s without a pulse,
then an independently initialized 52 -> 35 m/s downsweep in 20 s with a tip
pulse triggered at 49 m/s. Edit `pulse_speed`, `pulse_coefficient`, and
`pulse_duration` in that file. It saves to `output/LCO_downsweep_pulse` and
does not include or execute `run_pazy_LCO.jl`. Set
`PAZY_LCO_DOWNSWEEP_PULSE_RUN=false` before including it to load without running.

Set `pulse_branches=:both` to excite both independently initialized branches,
`:upsweep` to excite only the upsweep, or `:none` for an unforced control.
`run_branches` still selects the simulations to execute; selecting a pulse
does not start an otherwise skipped branch. For the ascending branch, edit
`upsweep_pulse_speed`, `upsweep_pulse_coefficient`, and `upsweep_pulse_duration`.
The initial diagnostic values are 45 m/s, 0.1 N and 0.04 s. The descending
branch retains `pulse_speed`, `pulse_coefficient`, and `pulse_duration`.
At 40 s ascending duration, the 45 m/s trigger occurs at 23.5294 s.

The output folder follows the pulse selection: `LCO_downsweep_pulse`,
`LCO_both_pulses`, `LCO_upsweep_pulse`, or `LCO_no_pulse`. Use the same selection
when replotting saved results, or explicitly set `output_directory`.
Repeating a case in the same directory replaces that case's outputs; use a
different directory when comparing mesh, duration, or force parameters.

Pulses provide a controlled initial disturbance, not a correction guaranteeing
an LCO. Compare identical speed histories, dt, mesh and wake settings first.
Inspect post-pulse growth/decay and test smaller coefficients. The bipolar
0.1 s waveform has a Fourier zero at 30 Hz; 0.04 s has stronger content near
the 28-31 Hz frequencies in the earlier diagnostics. Actual frequencies may
change with airspeed and discretization.

The same file has a `uvlm_options` block for chordwise/spanwise panel counts,
the maximum nondimensional timestep `time_step_chords`, optional
`fixed_time_step` [s], wake length [chords], optional `maximum_wake_rows`,
`core_radius` [m], `wake_shedding_fraction`, `symmetric_wing`, and
`aerodynamic_interaction`. These settings are applied to both the static
initialization and the dynamic calculation. An explicit wake-row count
overrides automatic sizing from wake length. A fixed timestep must obey the
nondimensional maximum and align both ramp endpoints. The resolved timestep,
wake capacity/age, and aerodynamic controls are recorded in the settings file.

At the top of the runner, select simulation and plotting independently:

```julia
run_branches = :both        # :both, :upsweep, or :downsweep
pulse_branches = :both      # :both, :upsweep, :downsweep, or :none
plot_only = false       # true: read saved CSVs without running either analysis
plot_branches = run_branches # or explicitly :both, :upsweep, :downsweep
```

When `plot_only=false`, `run_branches` selects which independent simulations
run. For example, `run_branches=:downsweep` skips the entire upsweep, including
its static initialization, and runs only the pulsed downsweep. The default
`:both` retains both analyses. `plot_branches` selects the newly computed
branches shown in the strain-versus-speed and time-history figures. Its default
follows `run_branches`; a plotting selection must include at least one computed
branch. Both configured ramps still determine timestep and wake sizing so that
single-branch runs can be compared consistently with the paired case. Skipped
branches return `nothing` in the result and produce no new initialization or
history CSV; existing files are preserved and never read for automatic plots.

When `plot_only=true`, only the strain-versus-speed
figure is regenerated from the selected saved CSVs; changing `uvlm_options`
does not alter existing histories. Single-branch figures have `_upsweep` or
`_downsweep` filename suffixes, preserving the combined figure.

For only a downsweep simulation and plot, use:

```julia
run_branches = :downsweep
plot_only = false
plot_branches = run_branches
```

To replot the saved results without rerunning this case:

```julia
include("lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/plot_pazy_LCO.jl")
plot_pazy_LCO(; plot_branches=:downsweep,
    output_directory=joinpath("lib", "WingPropellerUVLM", "examples",
    "pazy_aerobeams_coupled", "output", "LCO_downsweep_pulse"))
```

For custom calls, the driver accepts `upsweep_pulse_amplitude` [N],
`upsweep_pulse_start` [local upsweep s], `upsweep_pulse_duration` [s], and
`downsweep_pulse_amplitude` [N], `downsweep_pulse_start` [local downsweep s],
and `downsweep_pulse_duration` [s]. After loading the driver without execution
as above, this diagnostic example applies the pulse when the downsweep reaches
49 m/s:

```julia
downsweep_duration = 20.0
pulse_speed = 49.0
pulse_start = (52.0-pulse_speed)/(52.0-35.0)*downsweep_duration
pulsed_results = run_pazy_linear_LCO(;
    initial_speed=35.0, peak_speed=52.0,
    upsweep_duration=40.0, downsweep_duration,
    downsweep_pulse_amplitude=0.1,
    downsweep_pulse_start=pulse_start, # approximately 3.5294 s
    downsweep_pulse_duration=0.04,
    output_directory=joinpath("lib", "WingPropellerUVLM", "examples",
        "pazy_aerobeams_coupled", "output", "LCO_downsweep_pulse"))
```

The upsweep in the example above remains unforced. Both driver coefficients default to zero, preserving the
original independent analysis. For phase `(t-start)/duration` between zero and
one, the existing tip force along structural A1 is
`coefficient*sin(2pi*phase)*sinpi(phase)^2`; it is zero outside that interval.
This smooth bipolar pulse has zero net force impulse and a peak magnitude
about 0.6495 times the coefficient (about 0.065 N for a 0.1 N coefficient).
The speed ramp continues during the pulse. A nonzero pulse must fit entirely
inside its branch and span at least eight timesteps. Its applied history is
available as `pulsed_results.upsweep.tip_force_A1` or
`pulsed_results.downsweep.tip_force_A1` and the appended CSV column
`tip_force_A1_N`; settings record its parameters.

The 49 m/s trigger is a diagnostic choice, not an established instability
boundary for this model. To apply the pulse near the beginning instead, set
`downsweep_pulse_start=0.1`. A disturbance applied at 52 m/s may decay before
the wing reaches an unstable range. Compare with zero and smaller pulse
coefficients at identical numerical settings; inspect the response after the
pulse has ended before interpreting growth as self-excited oscillation.

`equilibrium_options` can override static convergence controls, for example
`(; maximum_static_iterations=150, maximum_aerodynamic_time=0.3)`. Defaults match
the standalone static analysis (100 outer iterations, 0.2 s maximum wake time
per fixed-geometry trial, force/circulation tolerances 1e-4, geometry tolerance
1e-5, and two consecutive converged outer iterations).

Both independent ramps and their static initializations use the same physical
timestep. Without `fixed_time_step`, it is sized at the highest speed of either
branch and slightly reduced to align each branch endpoint. An explicit
`fixed_time_step` must align both endpoints and obey the nominal maximum dt.
Simple fractional durations, such as `(60.0-35.0)/0.85`, are aligned using their
rational approximation. Duration ratios that would reduce the timestep below
half its nominal value are rejected before any simulation starts.
The default 10 x 20 UVLM mesh and five nominal chord lengths of wake are retained.
Static initialization adds computation, but its wake-development time is not
physical wing-motion time and is not added to either ramp duration.

Outputs in `output/LCO_independent_sweeps`:

- `pazy_linear_LCO_upsweep.csv` and `pazy_linear_LCO_downsweep.csv`: independent
  histories, each starting at time zero, with speed, root strain, raw root axial
  strain/curvatures, tip displacement, twist and incidence-plus-tip-twist AoA.
- `pazy_linear_LCO_upsweep_initialization.csv` and
  `pazy_linear_LCO_downsweep_initialization.csv`: static convergence histories.
- `pazy_linear_LCO_root_strain_vs_airspeed.png`: instantaneous root longitudinal
  strain [με] versus airspeed [m/s], blue upsweep / red downsweep, with a lower-left
  legend and no steady curve. All samples are connected in their original time
  order, retaining the oscillations rather than averaging them into an envelope.
- `pazy_linear_LCO_root_strain_history.png`: separate time plots for each branch.
- `pazy_linear_LCO_settings.txt`: controls and initialization convergence statuses.

`plot_pazy_LCO.jl` defaults to the independent results. Run it directly or use
the following to plot saved CSVs without rerunning the simulations:

```julia
include("lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/plot_pazy_LCO.jl")
plot_pazy_LCO()
```

If only one branch CSV exists, it is plotted and the missing branch is reported.
Run again after both branches finish to obtain the combined plot. To set the
reference figure's axis limits explicitly from the Julia REPL:

```julia
include("lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/plot_pazy_LCO.jl")
plot_pazy_LCO(; sweep_mode=:independent, xlims=(27.5, 52.5), ylims=(-2500, -500))
```

The default automatic axis limits show the full simulated strain range.

Independent runs have no concatenated physical time history. Continuous mode
writes the actual single-run history in its separate output directory.

The coupling retains near-root section values even with detailed histories
disabled. Surface strain uses local beam coordinates (+y forward, +z upper):

    strain_ue = 1e6 * (gamma_1 + gauge_z*kappa_2 - gauge_y*kappa_3)

Defaults `gauge_y=0`, `gauge_z=2.25e-3/2` place the gauge on the upper spar surface
at its chordwise center. The axial component of `gamma + kappa x offset` fixes
the sign: upper-surface compression under upward bending is negative. The first
element midpoint approximates the root station. Confirm experimental gauge
location and section conventions before quantitative comparison.

This changes the sweep and initialization procedure. The current beam coupling
still uses zero gravity and UVLM aerodynamics, while the paper includes gravity
and compares linear and dynamic-stall aerodynamic models. The dashed steady
reference curve requires separate converged static solutions over the speed
range, using the same strain reconstruction. A dynamic mean is not that curve.

The earlier staircase analysis remains available as `run_pazy_LCO()` after
loading with execution disabled, and writes to `output/LCO_sweep`.

## Constant-speed diagnostic

Run the separate pulse/sensitivity study without changing either ramp:

```powershell
julia --startup-file=no --compiled-modules=existing --project=. lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/run_pazy_LCO_diagnostic.jl
```

The default is 45 m/s at 5 degrees for 5 s. Four cases compare a small pulse,
an unforced control, half the timestep, and twice the retained wake length.
Each distinct speed/timestep/wake combination starts from its own checked static
equilibrium. The unforced control copies the same equilibrium as the baseline.
The pulse is the existing smooth, zero-integral tip-force waveform in the beam
coupling: coefficient 0.1 N, starting at 0.1 s, lasting 0.04 s. The coefficient
is a multiplier; the actual force extrema are smaller than 0.1 N.

Baseline dt and wake capacity match the default sweep (dt about 0.000475477 s,
30 wake rows). At 45 m/s this retains about 6.48 chords. Half dt uses 60 rows
to preserve the same wake age and physical length; double wake uses the baseline
dt with 60 rows. The 10 x 20 mesh and coupling tolerances remain the same.
The direct diagnostic process uses one BLAS thread to avoid overhead on small
dense systems. Existing Julia sessions are unaffected.

Results are written after each completed case to
`output/LCO_constant_speed_diagnostic`: per-case `history.csv`,
`initialization.csv`, and `windows.csv`, plus `diagnostic_summary.csv`,
`diagnostic_settings.txt`, `strain_histories.png`, and `window_amplitudes.png`.
The strain history plot subtracts each case's initial equilibrium value for
comparison; the CSV retains absolute surface strain and the applied pulse.
Window statistics exclude the pulse and remove a linear drift before computing
half range and RMS. The final four windows give a log-RMS slope in 1/s:
positive indicates growth over those windows, negative indicates decay.
These values alone do not establish a settled LCO or a converged numerical result.
Extend the run if the window amplitudes are still changing; compare refinements
over the same physical time and with the same pulse before drawing conclusions.

To select additional speeds or extend the duration from the Julia REPL:

```julia
include("lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/run_pazy_LCO_diagnostic.jl")
BLAS.set_num_threads(1)
run_pazy_LCO_diagnostic(; airspeeds=(44.0, 45.0, 47.0), duration=10.0,
    output_directory="output/pazy_LCO_diagnostic_extended")
```

The short diagnostic regression uses a coarse mesh and short durations only to
check force application, constant speed, equilibrium reuse, preserved wake age
under timestep refinement, coupling convergence, and exports:

```powershell
julia --startup-file=no --compiled-modules=existing --project=. lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/test_pazy_LCO_diagnostic.jl
```

To analyze the completed 45 m/s CSVs without rerunning the solver:

```powershell
julia --startup-file=no --compiled-modules=existing --project=. lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/analyze_pazy_LCO_diagnostic.jl
```

The analysis checks the stored strain/force histories and writes
`analysis_metrics.csv`, `analysis_comparison.png`, `analysis_late_responses.png`,
and `analysis_report.md` in the diagnostic output directory. The report discusses
the completed baseline, unforced, half-timestep, and double-wake study; its growth
and frequency estimates describe finite time windows, not converged eigenvalues
or guaranteed settled LCO amplitudes.

Regression checks exercise two independent static initializations and ramps,
verify copied beam/wake states and zero starting velocity, require independent
clocks, reject unconverged initialization, and check exports and plots:

```powershell
julia --startup-file=no --compiled-modules=existing --project=. lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/test_pazy_LCO.jl
```
