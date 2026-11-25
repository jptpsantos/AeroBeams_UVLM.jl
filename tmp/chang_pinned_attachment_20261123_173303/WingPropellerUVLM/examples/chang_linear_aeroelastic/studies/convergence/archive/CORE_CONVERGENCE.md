# Fixed-radius convergence study

Use `run_chang_core_convergence.jl` to inspect the sensitivity of the Chang
wing–propeller model to a **core radius specified in metres**. The study reads
the geometry, operating point, collective pitch, and interaction switch from
`../../chang_case.jl`. It ignores unrelated `CHANG_*` environment overrides.
The existing factor-based sweep remains a separate workflow.

Wing image symmetry is taken from `wing.symmetric` in `chang_case.jl` in
both modes. Propeller blade symmetry stays disabled. Changing the switch
invalidates cached aerodynamic results. The separate factor-based scripts
use `CHANG_AERO_WING_SYMMETRIC`; their gated aeroelastic stage inherits that
value from the aerodynamic metadata.

## Run

1. Edit `core_convergence_case.jl`. Start with `mode = :aerodynamic`.
2. Run from the repository root:

   ```powershell
   julia --startup-file=no --project=lib/WingPropellerUVLM lib/WingPropellerUVLM/examples/chang_linear_aeroelastic/studies/convergence/run_chang_core_convergence.jl
   ```

3. The current preset enables all aerodynamic families with `dry_run = false`,
   so the command executes the full study.
4. To preview first, set `dry_run = true`. This writes `case_matrix.csv` without
   simulating. Check panel counts, radii, time steps, and wake rows there,
   then restore `dry_run = false` to execute.

These are editable starting values, not validated resolutions or core radii.
The current case radius is automatically included in the radius list and is
the nominal radius for the separate discretization families. For example,
to explore millimetre-scale radii independently of the current case value:

```julia
radii_m = [0.01, 0.003, 0.001, 0.0003], # largest to smallest
nominal_radius_m = 0.001,
families = [:core_mesh, :time_step, :wake_length],
```

Define absolute panel counts in the same settings file:

```julia
# Baseline for wake/time-step runs and unchanged mesh directions.
nominal_mesh = (wing_span=20, wing_chord=5, prop_radial=5, prop_chord=5),

# Independent one-direction refinement lists; lengths can differ.
panel_levels = (
    wing_span = [20, 30, 40],
    wing_chord = [5, 8, 12],
    prop_radial = [5, 10, 15],
    prop_chord = [3, 6, 9],
),

# Separate combined meshes used by the core-radius sweep.
core_meshes = [
    (wing_span=20, wing_chord=5,  prop_radial=5,  prop_chord=3),
    (wing_span=30, wing_chord=8,  prop_radial=10, prop_chord=6),
    (wing_span=40, wing_chord=12, prop_radial=15, prop_chord=9),
],
```

These three controls are independent. Editing `panel_levels` changes the
separate mesh families; edit `core_meshes` to change the meshes paired with
each radius. Counts are absolute and override the panel counts in
`chang_case.jl` for this study. Propeller counts are **per blade**.
The example above illustrates nonproportional counts, not a calibrated mesh.

Select families before launching if you want a smaller study. Three mesh
levels and at least three radii/time-step/wake levels are required in the
settings. Counts must be positive integers, with strictly increasing lists.
Combined meshes must refine at least one direction per row without coarsening
the others. Doubling both grid directions substantially increases solver cost. Cases
execute sequentially; identical cases shared by several families run once.
Rigid cases share one Julia session to avoid repeated compilation, with a
fresh solver state for each simulation.

## What is varied

| Family | Changed quantity | Held fixed |
|---|---|---|
| `core_mesh_1`, `_2`, `_3` | Complete radius sweep on each mesh | Physical geometry, time step, wake age |
| `mesh_at_core_1`, `_2`, ... | Mesh refinement at each fixed radius, reusing the same runs | Radius, physical geometry, time step, wake age |
| `wing_span` | Wing span panels | Radius and other numerical settings |
| `wing_chord` | Wing chord panels | Radius and other numerical settings |
| `prop_radial`, `prop_chord` | One blade-grid direction | Radius and other numerical settings |
| `time_step` | Azimuth increment | Radius, mesh, physical duration and retained wake age |
| `wake_length` | Retained age in rotor revolutions | Radius, mesh and time step |

Both wakes retain the requested age. Their row counts increase when the time
step decreases; a nonintegral row count is rounded up. This deliberately uses
a common wake-age convention rather than the separate row limits in the main
case. The actual values are recorded in the matrix.

The rigid stage supports the current single-propeller geometry and Imperial
near-field loads. Its wing is undeformed. It runs a fixed number of complete
revolutions and compares the final phase-averaged wing `CL`, propeller `CT`,
and propeller `CQ`. Startup duration must exceed wake filling plus the
comparison cycles; increase `simulated_revolutions` if periodicity fails.
RPM and collective are held fixed throughout: no thrust or windmilling
retrim is performed.

## Read the results

Results are written under `../../output/core_convergence_<mode>_<timestamp>/`
unless `output_directory` is set explicitly.

- `case_matrix.csv`: complete numerical settings, including radius in metres,
  panel counts, physical time step, and wake rows.
- `convergence.csv`: measured quantities, error relative to the last level,
  change from the preceding level, tolerances, periodicity, and acceptance.
- `core_sensitivity.png`: coefficients or fitted growth rates/frequencies
  versus radius, with a separate curve for each mesh.
- `convergence_errors.png`: each comparison error divided by its tolerance.
- `report.md` and `study_snapshot.txt`: interpretation, operating conditions,
  resolved settings, and source fingerprint.
- Each case directory: raw history, log, and aerodynamic summary/metadata.
  The rigid summary also records actual panel core radii and mesh-size ratios.
  Aeroelastic cases save `damping_fit.toml` with each channel's fit quality,
  growth rate, damping ratio, frequency, window sizes, and validity.

Tables are updated after each case, so partial results survive an interrupted
study. Plots are generated after the run. A row can finish successfully but
fail the periodicity or damping-fit check; examine `valid` and `accepted`.
The core plot shows completed simulations even when their periodicity fails.

For each aerodynamic quantity, the comparison error is the larger of the
mean difference and the phase-waveform RMS difference. Its tolerance is
`max(absolute_tolerance, relative_tolerance * abs(reference_mean))`.
Periodicity uses a separate absolute cycle-to-cycle RMS tolerance.

Acceptance requires all planned levels in that family to finish, valid
signals, agreement of the last two levels, and agreement of the candidate
and every subsequent level with the final reference. The finest result
cannot pass solely because its error against itself is zero. A failed family
does not prevent other families from being reported.

Acceptance is **within a family**, not a declaration that the whole model is
converged. In particular, a core-family plateau only measures sensitivity to
regularization. Check the `mesh_at_core_*` families at the same radius as
well: a flat core curve whose value moves with mesh refinement is insufficient. Use the
separate mesh/time/wake families at your candidate nominal radius, then extend
their levels until the final comparison settles. A smaller radius is not
automatically a better physical model.

For reuse, set `output_directory` to a consistent directory before running.
Rigid results are reused only when their options, numerical-source hashes,
and output hashes match; changed or incomplete files are recomputed. Changing
the case or numerical source invalidates reuse. Aeroelastic cases run fresh
on every launch. Keep studies with different operating points in separate
directories, and do not edit the case or source while a study is running.

## Inspect aeroelastic stability

After resolving basic aerodynamic sensitivity, set:

```julia
mode = :aeroelastic,
dry_run = false,
families = [:core_mesh, :time_step, :wake_length],
fit_start_s = nothing,       # after the pulse; defaults to at least 0.75 s
block_duration_s = 1.0,      # fixed physical window across time-step refinements
block_overlap = 0.9,
frequency_band_hz = (3.0, 6.0),
minimum_fit_r_squared = 0.8,
lambda_tolerance_per_s = 0.01,
frequency_tolerance_hz = 0.1,
```

This mode runs the production aeroelastic solver in a fresh Julia child for
each unique case. It preserves the main case's duration, excitation,
structural properties, force/moment choices, and integration/coupling settings.
The wing span mesh remains fixed at `nominal_mesh.wing_span` because it also
changes the structural discretization. This is also the structural element
count used for every aeroelastic case in this study. `wing_span` is omitted,
and the wing-span entries in `core_meshes` are overridden; combined mesh refinement changes wing
chord panels and both blade-grid directions. Structural mesh convergence is
a separate study.

The moving-block FFT method fits the logarithm of modal amplitude against
time. `lambda > 0` indicates growth and `lambda < 0` decay. Both pitch and yaw
growth rates and their tracked frequencies must agree for acceptance.
Choose the fit interval after the forcing ends, use a sufficiently long
response, and choose a frequency band appropriate to the mode of interest.
The example 3–6 Hz band is a starting assumption. Inspect raw pitch/yaw
histories for beating, mode changes, and departure from a small-amplitude
response. Agreement at one speed is not a flutter-speed convergence result;
repeat at speeds around the stability crossing.

## Verification

`test/core_convergence.jl` exercises planning, metadata, and convergence
decisions as part of the package tests. A bounded solver/report check is:

```powershell
julia --startup-file=no --project=lib/WingPropellerUVLM lib/WingPropellerUVLM/examples/chang_linear_aeroelastic/validation/verify_chang_core_convergence.jl
```

It runs small rigid cases, verifies histories, plots, and reuse, and checks
the aeroelastic child configuration and preview. Its deliberately coarse
results in `output/core_convergence_smoke` are software checks, not a
production convergence study.
