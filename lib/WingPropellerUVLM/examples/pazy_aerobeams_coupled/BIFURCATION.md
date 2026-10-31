# Independent time-domain bifurcation analysis

Run from the repository root:

```powershell
julia --project=. lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/run_pazy_bifurcation.jl
```

Edit the named tuple at the bottom of that file before a full study. The current
entry point uses 40:0.25:50 m/s, root incidence 5 degrees, 10 s at each fixed speed,
statistics over `t >= 8 s`, zero pulse, serial execution, and no saved time histories.
The function defaults use 35:1:52 m/s, 10 s with `t >= 8 s`, and a 0.1 N pulse.
The exact upper
endpoint is included even if the speed increment does not divide the interval.
This is a separate runner; existing ramp and continuous LCO scripts retain
their behavior. Including the new runner loads functions without executing.
`enabled=false` returns without constructing models or creating outputs.

```julia
include("lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/run_pazy_bifurcation.jl")
bifurcation = (
    enabled=true, U_min=40.0, U_max=48.0, dU=0.5,
    simulation_time=10.0, transient_time=8.0,
    angle_of_attack=deg2rad(5.0),
    responses=(:tip_out_of_plane_percent, :tip_twist_deg, :root_strain_ue),
    save_time_histories=false, plot_results=true,
    zoom_limits=(42.0, 47.0), plot_amplitude=true,
)
result = run_pazy_bifurcation(; bifurcation...)
println(result.summary_path)
```

## Initialization and existing solver

`PazyBifurcation.jl` calls `run_pazy_wing_uvlm_static` and
`run_pazy_wing_uvlm` at each prescribed speed. It reuses the geometrically
nonlinear beam, time integration, strong/loose coupling, UVLM, and free wake.
No equations or algorithms in those solvers were changed for this analysis.

The default `initialization=:equilibrium` independently solves the coupled
static state and develops its wake at each speed. Its convergence must pass
before dynamics starts. The dynamic solver copies that speed's equilibrium;
structural velocities start at zero. Every dynamic clock starts at zero and
freestream velocity remains constant, including the initial sample. No final
beam, circulation, or wake state passes between airspeeds. With
`initialization=:undeformed`, each call instead uses the original unloaded
startup and a fresh empty wake; a static reference is then unavailable.

The default disturbance is the **existing** bipolar tip-force pulse, coefficient
0.1 N, starting at 0.1 s and lasting 0.04 s. Its peak magnitude is approximately
0.6495 times the coefficient. The same settings apply independently at each
speed. Set `pulse_amplitude=0.0` for unforced runs. Nonzero pulses must finish
before `transient_time` and span at least eight integration steps.

At each speed the physical step follows `time_step_chords*chord/U`, shortened
slightly to land exactly at `simulation_time`. Wake capacity is recomputed
from `wake_length_chords`, or set explicitly with `maximum_wake_rows`.
Mesh, core radius, shedding fraction, symmetry and aerodynamic interaction
are shared by that speed's equilibrium and dynamics. Override numerical
controls using `dynamic_options=(coupling_scheme=:strong, ...)` and
`equilibrium_options=(maximum_static_iterations=150, ...)`. Only solver
control keys are accepted there, so state/speed overrides cannot accidentally
turn this into a continuation or ramp study.

## Scalar outputs and coordinate conventions

| Name | Unit | Definition |
| --- | --- | --- |
| `tip_out_of_plane_m` | m | Existing tip reference-node output, `-u_n2[1]` in AeroBeams A axes |
| `tip_out_of_plane_percent` | % semispan | The same displacement divided by semispan, times 100 |
| `tip_in_plane_m` | m | `-u_n2[2]`, positive downstream in the existing UVLM coordinate map |
| `tip_twist_deg` | deg | Existing `wingtip_twist_degrees` rotation reconstruction |
| `root_strain_ue` | microstrain | `1e6*(gamma_1 + gauge_z*kappa_2 - gauge_y*kappa_3)` at the first element midpoint |
| `tip_geometric_incidence_deg` | deg | Root incidence plus elastic tip twist |

Gauge offsets use local section axes (+y forward, +z upper surface); defaults
are y=0 and z=0.001125 m. The strain is a near-root proxy at the first element
midpoint. Tip displacement refers to the existing beam reference node, whereas
the papers report midchord displacement. Geometric incidence omits plunge and
wake effects and must not be called effective angle of attack.

## Extrema, convergence and interpretation

The driver retains raw samples satisfying `t >= transient_time`, with no
filtering, detrending, absolute value, or static subtraction. It saves their
minimum, maximum, arithmetic mean and half peak-to-peak amplitude. The actual
window bounds and sample count are included in the CSV. At least two samples
must remain. Optional convergence checks use these statuses:

| Status | Meaning |
| --- | --- |
| `stationary` | Entire final-window range is at most twice the absolute tolerance |
| `converged` | Minima and maxima across the last three complete peak-to-peak cycles agree within tolerance |
| `unconverged` | Enough cycles are present but their extrema still vary beyond tolerance |
| `insufficient_cycles` | Too few complete cycles for the check |
| `not_checked` | `check_convergence=false` |

For the cycle check, tolerance is `atol + rtol*mean_cycle_amplitude`, with
`rtol=0.02` and `convergence_cycles=3` by default. A constant static offset
does not inflate the tolerance. Absolute tolerances are in each response's
units (see `RESPONSE_DEFINITIONS`); override selected definitions using
`response_definitions=merge(PazyBifurcation.RESPONSE_DEFINITIONS,
(tip_twist_deg=(unit="deg", label="Tip twist", atol=0.005),))`.

All statuses are saved. Figures mark every result whose convergence was not
verified with red crosses. This simple check tests late-cycle stationarity;
it is not a proof of asymptotic stability and may flag quasiperiodic or chaotic
responses. Extend the run/window and check timestep, mesh and wake sensitivity
before interpreting a branch as settled. The mean oscillatory response is
saved separately from the independently calculated static equilibrium.

The summary retains one overall minimum and maximum per speed for the response
envelope. With `save_local_extrema=true` (default), a separate compact CSV saves
**every local minimum and maximum** within the same final time window, including
sample time and signed response value. Flat peaks count once, at their last
sample. Window endpoints do not count as local peaks; strictly stationary or
monotone windows yield no artificial oscillatory extrema. No smoothing, peak
threshold, or convergence-based removal is applied. The paper-style plots use
these local peaks, so modulated responses can produce multiple points per speed.
Neither extrema extraction nor the convergence diagnostic identifies Hopf points.
Independent runs sample the attraction basin selected by their initialization
and pulse. They do not trace both coexisting branches of a hysteresis loop.

## Saved data and replotting

Default output directories are unique timestamped folders under
`output/bifurcation`. An explicitly supplied directory must be empty, preventing
accidental mixing with an earlier run.

- `settings.txt`: speed grid, observation/convergence settings and model controls.
- `bifurcation_summary.csv`: one row per speed and selected response, with
  airspeed, variable, units, min, max, mean, amplitude, status, convergence flag,
  actual window bounds, sample count, static equilibrium and display label.
- `case_NNNN_extrema.csv`: compact per-speed results, also preserving completed
  cases if a later solve fails. Serial execution updates the combined summary
  after each speed; process execution writes it when all tasks complete.
- `bifurcation_local_extrema.csv` and `case_NNNN_local_extrema.csv`: every retained
  local extremum, with columns `airspeed_mps,response,unit,kind,time_s,value`.
  These do not require full saved time histories. Set `save_local_extrema=false`
  to disable them.
- `case_NNNN_history.csv`: selected scalar histories only, when
  `save_time_histories=true`. Full beam/load/lattice/animation histories are
  disabled in the Pazy adapter. One run's scalar arrays are released after
  extraction; complete solver objects are not accumulated across speeds.
- `bifurcation_RESPONSE.png`, optional `_zoom.png` and `_amplitude.png`.
- `pazy_bifurcation.png` / `.pdf`: the four-panel paper-style figure, plus
  `pazy_bifurcation_RESPONSE.png` / `.pdf` for each individual panel.

## Four-panel figure matching the supplied style

The 5-degree runner enables `plot_paper=true`. Panels show tip out-of-plane
displacement (% semispan), near-root longitudinal strain (microstrain), tip twist
(degrees), and **geometric tip incidence** (root incidence + elastic twist).
The fourth panel is explicitly labelled geometric incidence: the existing
UVLM output does not provide the reference figure's effective tip angle of
attack. It does not include plunge velocity or induced-flow effects.

The figure uses orange for 5 degrees, downward/upward triangles for all local
minima/maxima, a dashed curve for independently solved static equilibrium,
serif text, light grids, captions below each panel, and a two-column legend.
The combined figure is 1100 by 860 pixels at the base plotting size, exported
at 300 dpi. Plot axes follow the supplied reference exactly:

| Panel | Y range | Y tick step |
| --- | --- | --- |
| Tip OOP displacement | 15 to 35 % semispan | 5 |
| Root longitudinal strain | -5000 to 0 microstrain | 1000 |
| Tip twist | -15 to 15 degrees | 5 |
| Geometric tip incidence | -10 to 20 degrees | 5 |

The reference x axes cover 30 to 60 m/s with ticks every 5 m/s. The current
runner selects a 40 to 50 m/s view and 0.25 m/s simulation increments
(41 independent speeds). The image does not establish a
simulation increment, timestep, pulse, or transient duration. Those numerical
settings remain explicit in the runner. Configure the axes with:

```julia
paper_plot_options=(
    xlims=(30.0,60.0), xticks=30:5:60,
    ylims=((15.0,35.0), (-5000.0,0.0), (-15.0,15.0), (-10.0,20.0)),
    yticks=(15:5:35, -5000:1000:0, -15:5:15, -10:5:20),
    formats=("png","pdf"), save_panels=true,
    fontfamily="Times", captionfontfamily="Times Bold",
    guidefontsize=12, tickfontsize=11, legendfontsize=10, captionfontsize=11,
    linewidth=0.55, markersize=1.35, markerstrokewidth=0.4,
    figure_size=(1100,860), panel_size=(550,420), dpi=300, gridalpha=0.18,
)
```

These appearance controls are passed directly to `plot_pazy_bifurcation`;
no plotting-module edits are needed to tune them. `linewidth` controls the
dashed steady curves; extrema remain unconnected triangle markers. Marker
sizes and line widths are in plotting points. Font sizes are in points,
`figure_size`/`panel_size` set base (width,height), and `dpi` controls PNG
resolution. GR's built-in TeX renders the subscripted scientific labels
without an external LaTeX installation. The reference-style outputs begin
with `pazy_bifurcation`; the older `bifurcation_` outputs are diagnostics.

The same reference axes apply to CSV-only replotting by default. Fixed axes may
clip results outside the reference ranges; all raw extrema remain saved in CSV.
For automatic axes use `xlims=nothing, xticks=:auto,
ylims=(nothing,nothing,nothing,nothing), yticks=(:auto,:auto,:auto,:auto)`.

All convergence statuses remain in the summary and the existing diagnostic
plots. The paper-style figure plots raw local extrema without convergence-cross
overlays; inspect the diagnostic plots/CSV before interpreting a branch.
Hopf diamonds are drawn only for independently identified, explicitly supplied
points. The legend includes the reference's Hopf symbol key by default even
when no points have been supplied; `show_hopf_key=false` hides that unused key.
Supply actual airspeed and response values rather than treating an
amplitude threshold as a Hopf detection:

```julia
plot_pazy_bifurcation(result.summary_path;
    hopf_points=((angle_of_attack_deg=5.0, airspeed_mps=U_hopf,
        responses=(tip_out_of_plane_percent=u3_hopf, root_strain_ue=strain_hopf,
            tip_twist_deg=twist_hopf, tip_geometric_incidence_deg=incidence_hopf)),))
```

Regenerate the paper-style figures from CSV without loading an aeroelastic solver:

```powershell
julia --project=. lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/plot_pazy_bifurcation.jl path/to/bifurcation_summary.csv
```

The plotter reads root incidence from the run's `settings.txt`. Without that file,
supply `angle_of_attack_deg=5.0` programmatically. Old results without a local-peak
CSV still plot, with a warning and legend labels identifying their final-window
envelopes. Missing geometric incidence can be reconstructed exactly from the
saved tip twist and root incidence. Missing local peaks cannot be recovered from
summary statistics alone; rerun the analysis to save them.

Later, overlay other incidences by passing multiple summary paths:

```julia
include("lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/plot_pazy_bifurcation.jl")
plot_pazy_bifurcation([summary_3deg, summary_5deg, summary_7deg];
    output_directory="tmp/pazy_bifurcation_comparison")
```

Colors are blue/orange/yellow for 3/5/7 degrees. Individual files and the combined
figure are exported as PNG and PDF by default; SVG is also supported. Set
`save_panels=false` to export only the combined figure.

The existing envelope/amplitude plot API remains `plot_bifurcation(result.summary_path;
zoom_limits=(42.0,47.0), plot_amplitude=true)`. Static curves use only actual
equilibrium values; unavailable values are NaN and omitted.

## Optional process parallelism

```julia
using Distributed
ids = addprocs(2; exeflags=`--project=$(dirname(Base.active_project()))`)
try
    parallel_result = run_pazy_bifurcation(; bifurcation..., worker_ids=ids)
finally
    rmprocs(ids)
end
```

Workers load the adapter before tasks are serialized. Each process constructs
independent model and wake data and writes uniquely numbered case files. Workers
must share access to the source tree and output directory (local workers do).
Plotting takes place only on the main process. Empty `worker_ids` means serial.

## Wing–propeller cases and reusable analysis

`../bifurcation/TimeDomainBifurcation.jl` is a solver-independent module. Its
`run_bifurcation(simulate; ...)` calls `simulate(U)` once per speed, expecting
`(time, responses, static_equilibrium)` with a NamedTuple of scalar response
vectors and an optional NamedTuple of independently computed static values.
The callback must initialize every state afresh. This interface also accepts
additional structural generalized coordinates without modifying a solver.

For the existing general wing–propeller cases, use the supplied adapter:

```julia
include("lib/WingPropellerUVLM/examples/wing_propeller_aeroelastic/WingPropellerBifurcation.jl")
using .WingPropellerBifurcation: run_wing_propeller_bifurcation
case = WingPropellerBifurcation.WingPropellerAeroelastic.xu_case()
propeller_result = run_wing_propeller_bifurcation(case;
    U_min=60.0, U_max=110.0, dU=1.0,
    simulation_time=10.0, transient_time=8.0,
    responses=(:tip_out_of_plane_m, :tip_twist_rad,
        :propeller_1_pitch_rad, :propeller_1_yaw_rad))
```

This adapter copies the case and updates its physical operating speed before
calling the existing structural assembly and UVLM solver. Their existing speed
law recomputes angular speeds and gyroscopic terms for constant advance ratio;
fixed-RPM cases retain fixed RPM. Coordinates come from the existing
`full_beam_aerodynamic_kinematics` routine. Case integration/coupling settings
are retained. Use `solver_options` for the existing UVLM runner keywords,
including initial displacement/velocity, external load, finite core, interaction
and wake capacities. These runs use fresh empty wakes and the solver's default
zero structural initial conditions unless explicitly supplied. No static
equilibrium is inferred from the response mean. The existing general solver
retains full state histories during one run; this adapter releases them after
scalar extraction. Pazy itself has no propeller DOFs.

## Verification and reference scope

```powershell
julia --startup-file=no --compiled-modules=existing --project=. lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/test_pazy_bifurcation_plots.jl --keep-output
julia --startup-file=no --compiled-modules=existing --project=. lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/test_pazy_bifurcation.jl --parallel
julia --startup-file=no --compiled-modules=existing --project=. lib/WingPropellerUVLM/examples/wing_propeller_aeroelastic/validation/test_bifurcation.jl
```

Tests cover known periodic/stationary/growing signals, exclusion of startup
transients, CSV and plot round trips, short two-speed physical runs, a repeated
first-speed run after the second speed, and agreement between serial and
separate-process results. The millisecond physical runs test plumbing and
initialization, not settled LCO behavior. The full 10-second grid is deliberately
separate from these checks.
The CSV-only plotting test also verifies raw peak values/times, plateau handling,
multiple extrema per speed, stationary header-only peak files, legacy envelopes,
and combined/individual PNG/PDF exports. Its retained preview is synthetic and
is saved under `tmp/pazy_bifurcation_style_validation_*`, separate from studies.

The supplied Riso/Cesnik IFASD 2022 preprint, p. 13 and Fig. 11 (p. 15), states
the 10-second late-extrema procedure and describes static initialization followed
by a one-cosine vertical gust. The supplied SciTech 2024 paper, pp. 9–10/Fig. 6,
describes stepped speed sweeps and late-cycle extrema. The implementation follows
the user's explicitly requested independent-run protocol, using the existing
tip-force pulse rather than implementing the reference gust. These are reference
methodologies, not claims of identical models, forcing, or numerical results.

Validated in this workspace on 23 October 2026 with Julia 1.11.5: 75 new serial
analysis checks, 13 wing-propeller adapter checks, four separate-process/serial
comparison checks, and all 179 existing `test_pazy_LCO.jl` regression checks
passed. The generated figures were visually inspected. The parallel comparison
ran outside the execution sandbox after its worker-pipe startup restriction was
confirmed with a minimal probe. Full 10-second airspeed grids were not run.

## Revise errors while loading dependencies

### Unsupported bifurcation keywords

`got unsupported keyword arguments "save_local_extrema", "plot_paper", "paper_plot_options"`
means the adapter being called has the older function signature. The Pazy
adapter must accept those three keywords, forward `save_local_extrema` to the
shared driver, and invoke the paper-style plotter with `paper_plot_options`.
The entry point now checks its module API version and gives a restart message
when an older module is already retained in the IDE REPL.

After updating the source, start a fresh Julia process to replace any cached
older definitions. Check the runner's keyword compatibility without starting
simulations or creating result folders:

```powershell
julia --startup-file=no --project=. lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/run_pazy_bifurcation.jl --check-config
```

This checks the entry-point keyword dispatch with `enabled=false`; it does not
solve a structural model or verify numerical convergence. Omit `--check-config`
to run the configured study. Do not remove the plotting/extrema keywords to
work around a stale session, since they control the requested outputs.

### Dependency reload errors

Messages such as `Failed to revise PrettyTables`, `invalid redefinition of
constant PrettyTables.ColumnTable`, or `invalid redefinition of constant
Adapt.WrappedArray` identify a failure in the interactive session's live reload.
Later missing-type errors can follow the unsuccessful redefinitions. They do
not report an aerodynamic or structural convergence failure.

Restart the affected Julia REPL to rebuild its package state. On Julia 1.11,
Revise cannot replace changed struct definitions in place; see the
[Revise limitations](https://timholy.github.io/Revise.jl/stable/limitations/).
The exact trigger may be a package/source/environment change during an existing
session; the stack trace alone does not establish which change occurred.

For a long study, a fresh terminal process avoids the existing REPL state:

```powershell
julia --startup-file=no --project=. lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/run_pazy_bifurcation.jl
```

VS Code's Julia extension loads Revise automatically when `julia.useRevise`
is enabled (its default). To disable that automatic loading for new REPLs,
set `"julia.useRevise": false` in workspace settings, then restart the REPL.
Changing this setting does not repair a session that already has failed type
redefinitions. No solver tolerance or aerodynamic parameter adjustment is
indicated by these messages.

The coupling adapter prepends the local UVLM environment to `LOAD_PATH`. At
the 24 October check, the root manifest selected PrettyTables 2.4.0 / Adapt
3.7.2, while the UVLM manifest selected PrettyTables 3.4.8 / Adapt 4.7.0.
A fresh adapter load used the latter versions and succeeded. This is a reason
to keep a study in one fresh process rather than reuse a session across
environment changes; the version difference alone did not reproduce the
reported Revise failure in an isolated probe.

Read-only verification logs: `tmp/pazy_bifurcation_clean_load.log` records a
successful guarded entry-point load without Revise; `tmp/pazy_prettytables_clean_load.log`
records a successful table rendering with the root project's PrettyTables;
`tmp/pazy_revise_load_path.log` records the successful isolated Revise probe.
These checks load dependencies without launching the full bifurcation study.
