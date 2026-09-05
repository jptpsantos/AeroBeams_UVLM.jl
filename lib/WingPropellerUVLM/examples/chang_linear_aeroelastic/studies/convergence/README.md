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
propeller `CT` with the finest level:

- wing span panels: 20, 30, 40;
- wing chord panels: 5, 10, 20;
- propeller radial panels: 5, 10, 15;
- propeller chord panels: 5, 10, 15;
- retained wake: 0.5, 1, 2 revolutions;
- segment-core factor: 0.5, 0.25, 0.125; and
- propeller azimuth step: 10, 5, 2.5 degrees.

The finite-core family is a model-sensitivity check, not a claim that the
smallest core is automatically the most physical value.  The chord-based
core is held at zero because the corrected core law uses the local three-
dimensional vortex-segment length.

The default comparison uses a 1% relative tolerance with absolute floors of
`1e-4` for wing `CL` and `1e-6` for propeller `CT`.  Separate floors are used
because `CT` is only of order `1e-3` at this operating point; a `1e-4` thrust
floor would be too permissive.

The driver now streams each flushed child log to the terminal.  It rewrites
the CSV and Markdown report and refreshes the convergence PNG after every
completed case.  Existing results are reused only when their saved metadata
exactly match the requested mesh, wake, core, and time step.

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

## Stage 2: gated aeroelastic damping

Run `run_chang_aeroelastic_sweep_from_aerodynamic.jl` only after every
aerodynamic family is complete.  This driver reads
`coupled_aerodynamic_convergence.csv` and refuses to proceed if:

- a family is incomplete or contains a failed case;
- a result file does not match the numerical settings recorded in the table;
- no level and all levels finer than it satisfy the `CL`/`CT` tolerances and
  revolution-to-revolution drift limits; or
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
