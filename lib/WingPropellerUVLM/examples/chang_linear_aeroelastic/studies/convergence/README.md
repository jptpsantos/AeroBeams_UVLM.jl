# Independent Chang convergence studies

The entries now expose camelCase settings through `create_study()`. The shared
implementation is in `lib/WingPropellerUVLM/studies/src`; local `src/` files
forward to it for compatibility and provide the explicit aeroelastic entry
wrapper. Run with the studies project. See the
[complete LaTeX guide](../../../../../../docs/migration/UVLM_ANALYSIS_GUIDE.tex)
and [studies setup](../../../../studies/README.md). The snake_case examples below
remain supported through the legacy `aerodynamic_study()` / `aeroelastic_study()`
settings accessors. Including the aeroelastic entry defines its functions
without launching a study. The aerodynamic preset still calls `main()` at the
bottom; comment that call when loading only its definitions in the REPL.

For the staged finite-core, aerodynamic and damping investigation, see
[the analysis plan](FINITE_CORE_ANALYSIS_PLAN.md). It distinguishes existing
automation from additional checks and proposed extensions.

Only two entry files need editing:

| File | Purpose |
|---|---|
| [chang_convergence_aerodynamic.jl](chang_convergence_aerodynamic.jl) | Define independent physical/numerical inputs; converge propeller CT and CQ |
| [chang_convergence_aeroelastic.jl](chang_convergence_aeroelastic.jl) | Set explicit inputs or read a verified aerodynamic selection; converge pitch/yaw damping |
| `src/` | Shared-package compatibility includes and the explicit aeroelastic driver |
| `archive/` | Previous drivers, their editable presets, and reference notes |

Neither new study reads `chang_case.jl` or `CHANG_*` environment settings.
They reuse the existing UVLM and Chang structural solvers and tabulated
wing/blade properties. Numerical source changes invalidate old verification;
editing the interactive case does not.

## 1. Aerodynamic study

Edit `chang_convergence_aerodynamic.jl`:

- `physical`: wing/blade geometry, **actual RPM at the specified speed**,
  density, angles, collective, symmetry, and aerodynamic interaction.
- `core_mode`: `:fixed` for metres, `:chord` for a fraction of full local chord,
  or `:segment` for a fraction of each spanwise/radial vortex edge.
- `nominal`: absolute wing/blade panel counts, core value, wake age, and azimuth increment.
- `sweeps`: explicit counts/values in refinement order. Mesh and wake levels
  increase; core and azimuth levels decrease. Each enabled family needs >=3 levels.
- `families`: which sweeps to execute. A partial study produces inspection
  results; publishing a selection requires all seven families.
- `simulated_revolutions`, `averaged_revolutions`, tolerances, and output settings.

Propeller panel counts are per blade. Symmetry mirrors only the wing and its
wake across y=0; all physical blades remain explicit. Both wakes retain the
specified age in revolutions, so row counts increase as the time step decreases.
Each family changes one parameter while the others keep their nominal values.
Fixed-core mode separates dimensional radius from mesh refinement; segment
mode changes the radius when its defining edge length changes.

The supplied `interaction=false` disables wing/propeller aerodynamic interaction.
In that configuration CT/CQ cannot establish convergence of the wing loads:
changing the wing mesh need not change either propeller coefficient. Enable
interaction when studying the coupled aerodynamics, and inspect wing loads
separately if wing aerodynamic convergence is required. Stage 2 checks the
resulting damping sensitivity.

Run from the repository root:

Direct execution activates `lib/WingPropellerUVLM/studies`. In the IDE/REPL,
activate that project before including an entry. The aeroelastic entry can also
activate it automatically when the studies package is unavailable (see below).

```powershell
julia --startup-file=no --project=lib/WingPropellerUVLM/studies lib/WingPropellerUVLM/examples/chang_linear_aeroelastic/studies/convergence/chang_convergence_aerodynamic.jl
```

`dry_run=true` writes only the proposed case matrix. The supplied preset uses
`false`, so the command runs simulations. Cases run sequentially with fresh
solver states, sharing compilation and reusing duplicate cases. Existing rigid
histories are reused only when their options, source and file hashes match.
The example values are starting points, not established converged settings.

Results go to `output/convergence_aerodynamic/` relative to the Chang example:

- `case_matrix.csv`: every numerical setting, physical time step, and wake rows.
- `convergence.csv`: CT/CQ means, waveform errors, previous-level changes,
  periodicity, validity and acceptance.
- `convergence.png`: CT/CQ against refinement level for each family.
- Per-case raw histories, logs, summaries and metadata. CL remains available
  in raw histories but does **not** control this CT/CQ study's acceptance.
- `study.toml` and `report.md`: frozen study inputs and selection status.

To plot existing results without running simulations:

```powershell
julia --startup-file=no --project=lib/WingPropellerUVLM/studies lib/WingPropellerUVLM/examples/chang_linear_aeroelastic/studies/convergence/plot_chang_convergence_aerodynamic.jl
```

An optional first argument selects another saved results directory. PNG and
SVG figures for mean CT, mean CQ and successive-level errors are written to
`plots/` inside that directory. Axes use the saved case matrix; `plot_notes.md`
records the saved operating point and acceptance status. Including the plotting
file defines `ChangConvergenceAerodynamicPlots.plot_saved_convergence(directory)`
without running it.

For each coefficient, error is the larger of the mean difference and the
phase-waveform RMS difference. Tolerance is the larger of the absolute
tolerance and the relative tolerance times the absolute reference mean.
Cycle-to-cycle periodicity has its own absolute RMS limit. A candidate needs
a complete family, valid signals, agreement of the last two levels, and
agreement of all subsequent levels with the final reference. A finest case
cannot pass merely by comparing against itself.

When all families pass, the driver selects the first accepted value from
each family. It then runs that **combined candidate** and a finer combined
mesh/wake/time-step case, keeping the selected core law/value. Only if both
are periodic and agree does it write a verified `aerodynamic_selection.toml`.
Otherwise the selection status remains `not_converged`, with a reason in the
report. Increase the settling duration or extend/refine the failing ranges
before starting stage 2. Selecting a numerical core plateau does not calibrate
a physical vortex radius. RPM and collective remain fixed; no retrim is performed.

## 2. Aeroelastic study

The aeroelastic entry also activates the local studies project when loaded with
`include` or from the IDE if `WingPropellerUVLMStudies` is unavailable in the
current environment. This is a repository-local package; `Pkg.add` is unnecessary.
Its execution guard runs a study only when the file is executed directly, so
`include` is safe for loading definitions for manual REPL use.

### Full explicit analysis and runtime planning

The default `runMode=:full` runs the editable nominal settings and sweeps with
`coupling.scheme=:loose_explicit`. Each physical step evaluates UVLM once at the
target time and rotor azimuth, using the known structural displacement `U[n]`,
then advances the structure and commits the wake once. The structural integrator
remains Newmark-beta: this is an explicit partitioned exchange between the
aerodynamic and structural solvers. `maximumIterations` does not affect this
scheme because there are no coupling subiterations. The equilibrium tolerances
still check the structural solve against its prescribed aerodynamic load.

The full profile uses 4 s instead of the previous 2 s: a saved 2 s case took about
10 hours and retained only six positive peaks, below `minimumPeaks=8`. A longer
record is still needed for damping estimation; 4 s supplies more cycles but
does not guarantee an acceptable fit. The four mesh families are retained, and
`azimuthStepDeg=[5.0, 2.5, 1.25]` adds a three-level time-step sweep to check
sensitivity to the lagged explicit exchange. The five families produce 15 case
entries and 11 unique simulations because their nominal cases coincide. Full
results go to `output/convergence_aeroelastic_explicit/`.

Select `--preview` for a very short coupled response. It uses
0.15 s, a 4-by-2 wing mesh, four physical blades with one chordwise panel and
radial counts `[2, 3, 4]`, a 10-degree azimuth increment, and 0.5 retained wake
revolutions (18 rows). At the supplied RPM this gives 110 time steps per case,
330 across the three cases, and 16/20/24 bound panels respectively. Plots are
disabled and results go to `output/convergence_aeroelastic_explicit_preview/`,
keeping the full-study directory separate. Preview damping validity is forced to
`false`: this short interval checks the solver and its response output, and
cannot establish convergence of damping at 3--6 Hz.

From the repository root, run the default full analysis, select the preview,
or inspect either plan without solving:

```powershell
julia --startup-file=no --project=lib/WingPropellerUVLM/studies lib/WingPropellerUVLM/examples/chang_linear_aeroelastic/studies/convergence/chang_convergence_aeroelastic.jl
julia --startup-file=no --project=lib/WingPropellerUVLM/studies lib/WingPropellerUVLM/examples/chang_linear_aeroelastic/studies/convergence/chang_convergence_aeroelastic.jl --preview
julia --startup-file=no --project=lib/WingPropellerUVLM/studies lib/WingPropellerUVLM/examples/chang_linear_aeroelastic/studies/convergence/chang_convergence_aeroelastic.jl --dry-run
julia --startup-file=no --project=lib/WingPropellerUVLM/studies lib/WingPropellerUVLM/examples/chang_linear_aeroelastic/studies/convergence/chang_convergence_aeroelastic.jl --preview --dry-run
```

`--full` explicitly selects the default full profile. `--dry-run` writes
`runtime_plan.csv` with time steps, step counts, panel counts and wake capacities.
These are workload indicators, not predictions of elapsed wall time.

For REPL planning, including the entry does not run it:

```julia
include("lib/WingPropellerUVLM/examples/chang_linear_aeroelastic/studies/convergence/chang_convergence_aeroelastic.jl")
study = ChangConvergenceAeroelastic.create_study(; runMode=:full, dryRun=true)
ChangConvergenceAeroelastic.WingPropellerUVLMStudies.solve!(study)
```

Duration and the azimuth increment determine the number of physical steps:
`dt = azimuthStepDeg / (6 * rpm)`. A smaller increment adds steps and also
increases retained wake rows at fixed wake age. More bound panels increase
influence-matrix assembly and circulation-solve work; more wake panels increase
induced-velocity and wake-convection work. Use the preview to check setup before
paying for long records and refine each numerical control in the full study.
Duplicate nominal cases are already reused in memory within one study.
Check the azimuth family before accepting damping convergence with explicit
coupling: the lagged exchange can change growth rates and frequencies, and
stability of the structural integrator alone does not establish stability of
the coupled response. Keep the equilibrium tolerances while refining the time
step. Later kernel profiling can assess
reusable snapshot buffers, structural factorization reuse, unused surface-history
copies and per-step logging. No runtime speedup is assumed from the plan alone.

Edit `chang_convergence_aeroelastic.jl`. The supplied `selectionFile = nothing`
uses the entry's own `physical`, `coreMode`, `nominal`, and `sweeps`, just as in
the aerodynamic study. Every nominal control needs an absolute value. An
incomplete aerodynamic study does not block this independent damping study.
Console output, `study.toml`, `aerodynamic_input.toml`, and `report.md` explicitly
record that aerodynamic convergence is **unverified**. No saved aerodynamic
selection is read or certified. A span-panel sweep rebuilds both aerodynamic
and structural meshes and is identified as such in the report.

To retain the strict two-stage workflow, use `runMode=:full` and set `selectionFile` to the verified
`aerodynamic_selection.toml` path. Physical inputs and the core law are then
inherited from that selection; set individual nominal values to `nothing` to
inherit its numerical controls. This path rechecks source hashes, every raw
aerodynamic history, all seven family acceptance decisions, and the combined
refinement evidence. An `incomplete` or `not_converged` selection is rejected.
Wing span remains fixed in this verified-selection workflow.

To investigate a different aeroelastic speed using the same aerodynamic
selection, set `airspeed = 85.0` (for example) in the aeroelastic entry file.
`nothing` retains the aerodynamic speed. RPM scales in proportion to speed,
preserving the advance ratio; selected meshes, wake revolutions and core settings
remain the starting values for the sweeps. Choose a separate `output_directory`
for each speed. The time step and automatic excitation/fitting times use the
scaled RPM.

Run the entry script normally, or call
`ChangConvergenceAeroelastic.run_aeroelastic(settings)` in the REPL after including
it. `ChangConvergenceAeroelastic.create_study(settings; overrides...)` and
`WingPropellerUVLMStudies.solve!(study)` also support both entry workflows.
The nested `Convergence.run_aeroelastic(settings)` call retains the shared
verified-selection workflow.
In verified-selection mode the source TOML is verified without modification.
For an override,
`aerodynamic_input.toml` records the effective physical inputs and nests the
verified original selection under `aerodynamic_reference`; `report.md` records
both speeds and RPMs. Aerodynamic verification remains specific to the original
speed. The explicit driver lives outside the shared numerical sources, preserving
the source fingerprint of an aerodynamic study already in progress. Its own
driver fingerprint is recorded separately in `study.toml`.

Set the independent structural/pylon properties, response duration and impulse,
integration/coupling settings, and damping fit controls in the second file.
The table-based Chang wing model remains implemented in the shared model source.
`frequency_band_hz` must correspond to the modes being studied; the 3–6 Hz
example is an initial choice, not automatic mode identification.

The aeroelastic entry uses the moving-block formulation supplied in
`moving_block_2.jl`. Its implementation is in `MovingBlockDamping.jl` alongside
the entry file. It cuts the response between selected positive peaks, chooses
an FFT block occupying 25–50% of that segment, shifts the block **one sample**,
removes each block's mean, and fits a straight line to log dominant spectral
amplitude versus block start time. The slope is lambda in s⁻¹; negative means
decay and positive means growth. The equivalent damping ratio is also saved.
It uses rectangular blocks, without the previous Hann window. Polynomial and
peak-finding dependencies are not needed for the linear regression and local
maximum operation. The reference file's demonstration and flutter-speed
interpolation are not executed.

Edit the `damping` settings in the aeroelastic entry:

```julia
block_size_mode = :record_fraction,
block_size = 512,
size_ratio_lb = 0.25,
size_ratio_ub = 0.50,
peak_from_start = 1,
peak_from_end = 0,
```

`:record_fraction` follows the reference's block sizing and one-sample shift.
For a fixed physical window across time-step refinements, use `:duration`;
only that mode uses `block_duration_s` and `block_overlap`. In either mode,
the dominant non-DC frequency is sought over the full spectrum. The configured
`frequency_band_hz` checks whether it remains in the expected band; the code
does not choose an in-band leakage component in place of the dominant mode.
Minimum peak count and R² remain acceptance checks for the study.

Each case saves `pitch_moving_block.csv` and `yaw_moving_block.csv`, containing
block start times, log amplitudes, fitted values and dominant frequencies.
With `make_plots=true`, corresponding PNGs show the analyzed segment and fit.
`damping_fit.toml` includes the FFT frequency spacing (`frequency_resolution_hz`),
the growth rate, frequency, angular frequency, damping ratio and fit quality.
Frequencies remain FFT-bin estimates as in the reference: a fine convergence
tolerance alone does not confer sub-bin accuracy. Positive-peak trimming also
retains the reference's offset sensitivity; inspect the recorded fit curves.

For standalone postprocessing, include `MovingBlockDamping.jl` and call
`MovingBlockDamping.moving_block_damping(t, x)` to obtain the reference-compatible
tuple `(lambda, frequency_hz, omega_radps)`. Use `moving_block_metrics(t, x)`
for the summary and diagnostic arrays. Inputs must be finite and uniformly sampled.

The entry-level postprocessing extension uses a dedicated study wrapper and
separate estimator fingerprints in `study.toml` and `damping_fit.toml`. The
shared aerodynamic solver sources are unchanged, preserving existing verified
aerodynamic selections. The legacy metrics in `src/` remain available to
archived/internal workflows; the two aeroelastic entry calls above select the
new method.

For each damping sweep, `nothing` selects the aerodynamic value and two finer
levels unless overridden by `nominal`. The aeroelastic entry uses the same
`nominal`, `sweeps`, and `families = collect(keys(sweeps))` structure as the
aerodynamic entry:

- `nominal`: fixed controls for the other parameters during each sweep.
  Explicit mode requires an absolute value for every control. In verified mode,
  each `nothing` inherits the aerodynamic selection; set an absolute value to
  override it. The core law is also inherited in that mode.
- `sweeps`: add or comment entries to choose which parameters to vary. Explicit
  lists need at least three strictly ordered levels: increasing for panels/wake,
  decreasing for core/azimuth. Their first value may differ from the nominal
  or aerodynamic selection. A `nothing` entry enables automatic refinement
  from the resolved nominal; comment the entry to disable that family.

For example, set `nominal.propellerRadialPanels` to 10 and use
`propellerRadialPanels = [5, 10, 15]` in `sweeps`. The radial family runs those
three counts; all other families keep 10 radial panels. The report records
the resolved nominal controls, and the original aerodynamic selection is
preserved when using the verified-selection workflow. In that workflow wing
span stays fixed at the selected count because it also defines the structural
elements. Explicit mode permits a span sweep, which changes both meshes.

```powershell
julia --startup-file=no --project=lib/WingPropellerUVLM/studies lib/WingPropellerUVLM/examples/chang_linear_aeroelastic/studies/convergence/chang_convergence_aeroelastic.jl --full
```

Full-study results go to `output/convergence_aeroelastic_explicit/`. Each case
includes the full response history, `resolved_configuration.toml`, and
`damping_fit.toml` with
fit quality, frequency, lambda, and damping ratio. The summary compares both
pitch/yaw lambda and frequency. Positive lambda means growth; negative lambda
means decay. FFT windows follow the configured record-fraction or
physical-duration mode.
Fits with insufficient peaks or poor quality remain indeterminate. Inspect
raw signals for beating, mode changes, and nonlinear response amplitudes.

Aeroelastic cases always run fresh; duplicate cases within one study are
reused in memory. The solver receives a complete explicit configuration,
without launching the interactive entry point. The current production
shedding fraction is 0.1; the handoff rejects other aerodynamic values.
This study checks damping sensitivity at one operating point. To establish
flutter-speed convergence, repeat near the stability crossing.

## Previous scripts and verification

Old filenames now live under `archive/`; their relative solver paths and
test imports have been updated. Existing simulation output directories are
preserved. The archive retains the older CT/CQ/CL and case-dependent workflows;
their instructions describe those versions, not the new two-file workflow.

Run package checks with `julia --project=lib/WingPropellerUVLM/studies lib/WingPropellerUVLM/test/runtests.jl`.
The independent workflow tests cover planning, CT/CQ-only acceptance,
selection provenance, rejected stale/tampered results, and explicit model configuration.

`validation/verify_independent_convergence.jl` (relative to the Chang example)
runs three small rigid cases and a short production aeroelastic case. This
checks solver integration and rejection of insufficient damping data; it does
not establish aerodynamic or damping convergence for the supplied study preset.
