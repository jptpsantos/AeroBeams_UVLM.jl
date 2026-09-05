# Chang two-stage UVLM convergence workflow

The convergence study is split into two fronts so that expensive damping
simulations are not used to diagnose basic aerodynamic discretization errors.

| File | Role |
|:--|:--|
| `run_chang_coupled_aerodynamic_sweep.jl` | Recommended rigid aerodynamic sweep |
| `run_chang_coupled_aerodynamic_analysis.jl` | Isolated child case used by that sweep |
| `run_chang_aeroelastic_sweep_from_aerodynamic.jl` | Gated damping study |
| `run_chang_uvlm_convergence_sweep.jl` | Legacy one-stage entry point |
| `chang_aeroelastic_convergence.jl` | Shared damping metrics and process/report helpers |
| `plot_chang_uvlm_convergence_report.jl` | Rebuild plots for a legacy one-stage result |

## Stage 1: rigid coupled aerodynamics

Run `run_chang_coupled_aerodynamic_sweep.jl`.  The wing is rigid, the
propeller rotates, and wing--propeller aerodynamic interaction is enabled.
Each family changes one quantity and compares periodic wing `CL` and
propeller `CT` and `CQ` with the reference. Both the child and reuse checker
read shared option definitions. The default operating point is 65 m/s,
3 degrees angle of attack, 1217.6962 RPM, with interaction enabled:

- wing span panels: 20, 30, 40;
- wing chord panels: 10, 20, 30;
- propeller radial panels: 10, 15, 20;
- propeller chord panels: 10, 15, 20;
- retained wake: 1, 2, 3 revolutions;
- core/full-chord factor: 0.04, 0.02, 0.01, 0.005; and
- propeller azimuth step: 5, 2.5, 1 degrees.

Nominal grids are wing 30x10 and blades 10x10, with two retained wake
revolutions, a 5-degree step, and `rc = 0.01 c`. Eight revolutions are
simulated and the last two averaged. These are starting settings for a study,
not validated converged settings.

The default `CHANG_AERO_SWEEP_CORE_MODE=chord` keeps radius independent of
span/radial mesh. `c` is the full local chord, not a panel's chord. The primary
aeroelastic runner retains its editable mixed law; the gated driver transfers
both selected factors. See [the finite-core review](FINITE_CORE_REVIEW.md)
for dimensional radii, evidence, and remaining modeling limitations.

Core acceptance means insensitivity over the tested interval, not physical
validation of the smallest core. For a separate mesh-linked experiment, use
`CHANG_AERO_SWEEP_CORE_MODE=segment`: default factors are 0.5, 0.25, 0.125,
0.0625 with nominal 0.25. That mode couples span/radial mesh and core changes.
Only the factor selected by the mode is active; the other is forced to zero.

Mean and waveform comparisons use `max(absolute_tolerance, relative_tolerance
* abs(reference_mean))`: defaults are 1% relative, `1e-4` absolute for `CL`,
and `1e-6` for `CT` and `CQ`. Near zero thrust/torque, the absolute floor matters.

Acceptance also requires:

1. A complete history with matching source, physical/numerical options,
   sample count, and time/azimuth grids.
2. Wake filling before the averaging and periodicity windows: simulate at least
   `ceil(retained_revolutions) + max(averaged_revolutions, 2) + 1` revolutions.
3. Phase-aligned loads that repeat. The worst RMS change over the final
   `max(averaged_revolutions, 2)` cycle pairs must satisfy absolute limits of
   `1e-4`, `1e-6`, and `1e-6` for `CL`, `CT`, and `CQ` respectively.
4. All family levels present, agreement of the final two levels, and agreement
   of the candidate and every finer level with the reference in mean and waveform.

The finest case agreeing with itself cannot establish convergence. Waveform
standard deviations measure deterministic fluctuations, not confidence intervals.
Increase duration when periodicity fails; extend refinement when the final pair
disagrees. All comparison and periodicity limits are recorded in the CSV.

The driver now streams each flushed child log to the terminal.  It rewrites
the CSV and Markdown report and refreshes the convergence PNG after every
completed case. TOML metadata records all controls, including duration,
averaging, density, sideslip, collective pitch, and shedding fraction. SHA-256
fingerprints check numerical source, Julia version, manifest, and result files.
Old outputs without this metadata are recomputed. Changing tolerances only
recomputes comparisons; changing simulated settings or source invalidates reuse.

The PNG shows the larger mean/waveform error divided by tolerance for each
coefficient; the dashed line marks one and green stars mark accepted cases.
Children run without plotting dependencies. Set
`CHANG_AERO_SWEEP_PLOT_RESULTS=false` to disable sweep plots, or
`CHANG_AERO_PLOT_RESULTS=false` for a standalone rigid case.

From the repository root in PowerShell:

```powershell
$env:CHANG_AERO_SWEEP_LIVE_LOG = "true"
$env:CHANG_AERO_SWEEP_REUSE_EXISTING = "true"
Remove-Item Env:CHANG_AERO_SWEEP_DRY_RUN -ErrorAction SilentlyContinue

julia --project=lib/WingPropellerUVLM `
  lib/WingPropellerUVLM/examples/chang_linear_aeroelastic/studies/convergence/run_chang_coupled_aerodynamic_sweep.jl
```

To resume a stopped sweep, set `CHANG_AERO_SWEEP_OUTPUT_DIR` to its existing
directory before invoking the same command.  The metadata audit prevents a
result computed with one grid from being silently reused under another grid.

Refinement tables can be extended without editing Julia:

```powershell
$env:CHANG_AERO_SWEEP_FAMILIES = "time_step,wake_length"
$env:CHANG_AERO_SWEEP_AZIMUTH_LEVELS = "5,2.5,1,0.5"
$env:CHANG_AERO_SWEEP_WAKE_LEVELS = "2,3,4"
$env:CHANG_AERO_SWEEP_SIMULATED_REVOLUTIONS = "12"
$env:CHANG_AERO_SWEEP_AVERAGED_REVOLUTIONS = "3"
$env:CHANG_AERO_SWEEP_DRY_RUN = "true"
```

All level lists have prefix `CHANG_AERO_SWEEP_`: `WING_SPAN_LEVELS`,
`WING_CHORD_LEVELS`, `PROP_RADIAL_LEVELS`, `PROP_CHORD_LEVELS`, `WAKE_LEVELS`,
`CORE_LEVELS`, and `AZIMUTH_LEVELS`. Supply at least three strictly ordered
values. Core and azimuth decrease; other coordinates increase. Azimuth must
divide 360 degrees. Nominal meshes use the `CHANG_AERO_*_PANELS` controls.

Tolerance suffixes under the same prefix are `CL_ABS_TOL`, `CT_ABS_TOL`,
`CQ_ABS_TOL`, `REL_TOL`, `PERIODIC_CL_TOL`, `PERIODIC_CT_TOL`, and
`PERIODIC_CQ_TOL`. The old name `CHANG_AERO_WAKE_RELAXATION` actually controls
the trailing-edge shedding fraction eta, which the production adapter fixes
to 0.1.

## Stage 2: gated aeroelastic damping

Run `run_chang_aeroelastic_sweep_from_aerodynamic.jl` only after every
aerodynamic family is complete.  This driver reads
`coupled_aerodynamic_convergence.csv` and refuses to proceed if:

- a family is incomplete (including all four default core levels), duplicated,
  or contains a failed case;
- a result file does not match the numerical settings recorded in the table;
- no level and all levels finer than it satisfy the `CL`/`CT` tolerances and
  revolution-to-revolution waveform and mean-drift limits; or
- the 30-panel wing fails the aerodynamic spanwise tolerance.

The last condition is required because the current production coupling uses
the same spanwise grid for the beam and aerodynamic surface.  Structural
modal convergence established 30 beam elements.  Using 40 aerodynamic span
panels with 30 beam elements would first require a conservative,
work-conjugate noncollocated mapping.

Start with a dry run, which only writes and prints the selected case matrix:

```powershell
$env:CHANG_AEROELASTIC_AERO_SUMMARY = "FULL_PATH_TO\coupled_aerodynamic_convergence.csv"
$env:CHANG_AE_DRY_RUN = "true"
$env:CHANG_AE_VALIDATION_STAGE = "combined"

julia --project=lib/WingPropellerUVLM `
  lib/WingPropellerUVLM/examples/chang_linear_aeroelastic/studies/convergence/run_chang_aeroelastic_sweep_from_aerodynamic.jl
```

After checking `aeroelastic_case_matrix.md`, run the two combined cases:

```powershell
$env:CHANG_AE_DRY_RUN = "false"
$env:CHANG_AE_OUTPUT_DIR = "FULL_PATH_TO\output\gated_aeroelastic_combined"

julia --project=lib/WingPropellerUVLM `
  lib/WingPropellerUVLM/examples/chang_linear_aeroelastic/studies/convergence/run_chang_aeroelastic_sweep_from_aerodynamic.jl
```

If the selected and combined-refined damping slopes differ materially, run
the attribution stage in a new output directory:

```powershell
$env:CHANG_AE_VALIDATION_STAGE = "attribution"
$env:CHANG_AE_OUTPUT_DIR = "FULL_PATH_TO\output\gated_aeroelastic_attribution"

julia --project=lib/WingPropellerUVLM `
  lib/WingPropellerUVLM/examples/chang_linear_aeroelastic/studies/convergence/run_chang_aeroelastic_sweep_from_aerodynamic.jl
```

The default aeroelastic controls are deliberately fixed across cases:

- full wing--propeller aerodynamic interaction;
- generalized-alpha `rho_inf = 1.0` (no high-frequency algorithmic damping);
- partitioned-coupling relaxation `1.0`;
- 6 s response with a 0.8--6 s damping interval;
- 1.5 s moving FFT blocks with 90% overlap and a Hann window; and
- tracking restricted to 3--6 Hz for the pitch/yaw mode.

The direct comparison metric is the moving-block exponential slope `lambda`
in 1/s.  Both pitch and yaw fits must have `R2 >= 0.8`; the default numerical
convergence tolerance is `0.01 1/s`.  The CSV and PNG are refreshed after
each completed aeroelastic case.

The gate transfers both selected core factors into the aeroelastic child. It
currently requires the rigid study to use the production density (1.225 kg/m^3),
zero sideslip/collective offset, and eta=0.1. Rigid options explicitly define
the Chang geometry; reflect primary case geometry edits in those options too.

## Verification

```powershell
julia --project=lib/WingPropellerUVLM lib/WingPropellerUVLM/test/runtests.jl
julia --project=lib/WingPropellerUVLM `
  lib/WingPropellerUVLM/examples/chang_linear_aeroelastic/validation/verify_chang_aerodynamic_gate.jl
```
