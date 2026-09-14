# Independent Chang convergence studies

For the staged finite-core, aerodynamic and damping investigation, see
[the analysis plan](FINITE_CORE_ANALYSIS_PLAN.md). It distinguishes existing
automation from additional checks and proposed extensions.

Only two entry files need editing:

| File | Purpose |
|---|---|
| [chang_convergence_aerodynamic.jl](chang_convergence_aerodynamic.jl) | Define independent physical/numerical inputs; converge propeller CT and CQ |
| [chang_convergence_aeroelastic.jl](chang_convergence_aeroelastic.jl) | Read the verified aerodynamic selection; converge pitch/yaw damping |
| `src/` | Shared solver adapters, metrics, plots, and result verification |
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

Both entry scripts activate the local `lib/WingPropellerUVLM` project before
loading the solver, including when run or included from the IDE.

```powershell
julia --startup-file=no --project=lib/WingPropellerUVLM lib/WingPropellerUVLM/examples/chang_linear_aeroelastic/studies/convergence/chang_convergence_aerodynamic.jl
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

Edit `chang_convergence_aeroelastic.jl`. Its `selection_file` points to stage 1.
The physical wing/propeller/flow configuration is inherited from that saved
selection, so changing stage 1's script later does not silently change stage 2.
The handoff rechecks source hashes, every raw aerodynamic history, all family
acceptance decisions, and the combined refinement evidence.

To investigate a different aeroelastic speed using the same aerodynamic
selection, set `speed_mps = 85.0` (for example) in the aeroelastic entry file.
`nothing` retains the aerodynamic speed. RPM scales in proportion to speed,
preserving the advance ratio; selected meshes, wake revolutions and core settings
remain the starting values for the sweeps. Choose a separate `output_directory`
for each speed. The time step and automatic excitation/fitting times use the
scaled RPM.

Run the entry script normally, or call
`ChangConvergenceAeroelastic.run_aeroelastic(settings)` in the REPL after including
it. This entry-level function applies the override; the internal
`Convergence.run_aeroelastic` function uses the original operating point.
The source TOML is verified without modification. For an override,
`aerodynamic_input.toml` records the effective physical inputs and nests the
verified original selection under `aerodynamic_reference`; `report.md` records
both speeds and RPMs. Aerodynamic verification remains specific to the original
speed. This entry-file change preserves existing numerical source fingerprints.

Set the independent structural/pylon properties, response duration and impulse,
integration/coupling settings, and damping fit controls in the second file.
The table-based Chang wing model remains implemented in the shared model source.
`frequency_band_hz` must correspond to the modes being studied; the 3–6 Hz
example is an initial choice, not automatic mode identification.

For each damping sweep, `nothing` selects the aerodynamic value and two finer
levels. Alternatively, supply an absolute list, e.g. `prop_radial=[10,12,15]`
when 10 was selected. The list must begin at the aerodynamic selection.
Wing span stays fixed at the selected count because it also defines the
structural elements. Structural mesh convergence is a separate task.

```powershell
julia --startup-file=no --project=lib/WingPropellerUVLM lib/WingPropellerUVLM/examples/chang_linear_aeroelastic/studies/convergence/chang_convergence_aeroelastic.jl
```

Results go to `output/convergence_aeroelastic/`. Each case includes the full
response history, `resolved_configuration.toml`, and `damping_fit.toml` with
fit quality, frequency, lambda, and damping ratio. The summary compares both
pitch/yaw lambda and frequency. Positive lambda means growth; negative lambda
means decay. FFT windows use a fixed physical duration, rounded to a sample.
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

Run package checks with `julia --project=lib/WingPropellerUVLM lib/WingPropellerUVLM/test/runtests.jl`.
The independent workflow tests cover planning, CT/CQ-only acceptance,
selection provenance, rejected stale/tampered results, and explicit model configuration.

`validation/verify_independent_convergence.jl` (relative to the Chang example)
runs three small rigid cases and a short production aeroelastic case. This
checks solver integration and rejection of insufficient damping data; it does
not establish aerodynamic or damping convergence for the supplied study preset.
