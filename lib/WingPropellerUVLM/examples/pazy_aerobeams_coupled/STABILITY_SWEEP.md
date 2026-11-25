# Pazy flutter and LCO boundaries

`run_pazy_flutter_sweep.jl` runs the fixed-speed equilibrium flutter sweep
from 0° through 10° in 1° steps by default. Run from the repository root.

```powershell
$sweep = 'lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/run_pazy_flutter_sweep.jl'
$results = 'output/Pazy Flutter/aoa5'
julia --startup-file=no --compiled-modules=existing --threads=10 --project=. $sweep --check-config --angle 5 --output-directory $results
# Fixed-speed flutter sweep at 5 deg.
julia --startup-file=no --compiled-modules=existing --threads=10 --project=. $sweep --angle 5 --output-directory $results
# Fixed-speed flutter sweep at all configured angles.
julia --startup-file=no --compiled-modules=existing --threads=10 --project=. $sweep --output-directory 'output/Pazy Flutter/full-grid'
```

Set `--output-directory` to the new directory for that run. Absolute paths
are supported, and relative paths are resolved from the current directory.
The runner creates any missing parent directories. Omit the option for a
fresh `output/flutter_sweep/run_*/` directory beside the runner. Julia's thread
count is set at startup. The runner uses one BLAS thread during the sweep.

Edit `pazy_flutter_configuration()` in the runner. Each entry in `slices` specifies
the incidence, numerical pilot speed, independent pulse test speeds, LCO seed
and LCO speed limits. The starting configuration explicitly uses the refined
structural mesh with `rootTorsionCorr=false` for both trim and dynamics.
`--equilibrium` explicitly selects the default sweep. `--lco` runs the stepped
LCO paths, and `--all` runs both sweeps. `--pilot` runs the standalone timestep
and wake check at 5°. Add `--angle DEG` to select one configured incidence for
any of these modes or to inspect its schedule with `--check-config`. Each
simulation invocation creates a new output directory. In an IDE Julia REPL,
include the runner and call `main_pazy_flutter_sweep(["--angle","5"])`.

The initial sweep uses a 0.2 ms timestep. At 5° and 45 m/s, this setting
reached trim in the saved diagnostic, while both 0.1 ms wakes exceeded the
settling limit. The equilibrium bracket checks still test 0.1 ms at their
endpoints. The full grid uses three time integrations per speed before
refinement.

The solver assembles independent aerodynamic matrix rows across Julia threads,
checks wake core sizes once per calculation, and skips induction from wakes
with exactly zero circulation. Strong coupling reuses an aerodynamic trial
when its grid matches exactly. Each successful case still uses both pulses
and the unforced trace, with the pulsed runs attempted first. These changes
preserve the configured timestep, wake length and coupling tolerances.

The 0°–7° airspeed grids use the historical sweeps as a starting point.
The 8°–10° grids and LCO seeds are exploratory. The 0° LCO paths initially
search up to 90 m/s. Extend the speed limit if oscillations persist there
and a higher-speed extinction bracket is needed.

- **Numerical pilot.** `--pilot` compares the two finer timesteps and a longer
  wake at the slice's `pilot_speed`. `pilot.csv` records agreement in signed
  modal growth, frequency, bending/twist signature, and trim.
- **Equilibrium flutter.** Each fixed-speed case starts from a converged trim
  and developed wake. The fitter subtracts the matching unforced response,
  compares two pulse amplitudes, and checks two post-pulse windows in tip
  bending and twist. It uses Hann-weighted fits at a common frequency to
  reduce leakage from other modes and compares the bending/twist amplitude
  and phase across windows.
  Positive growth denotes instability. The driver refines opposite-sign
  brackets and checks their endpoints with a smaller timestep and a longer
  wake. A speed estimate requires adjacent, resolved fits of the same mode.
- **LCO boundaries.** Two paths start at the slice's `lco_seed` after verifying
  a periodic response and step toward either speed limit. Two further paths start from a
  verified quiet response at the outer speeds and traverse the full range.
  Each path carries beam and wake state and receives one initial pulse.
  Classification compares two separate late windows for amplitude, frequency,
  and waveform repetition. In `--all`, a fixed-speed control at the LCO seed
  speed can raise the configured amplitude floors. A finer speed path reruns
  from the original preparation and adds points inside coarse brackets. All
  LCO transitions are reported as speed intervals.

Outputs in the run directory include `pilot.csv` in the pilot,
`equilibrium_cases.csv`, per-case `fits.csv` and pulse histories,
`lco_classification.csv` for each completed portion of a path, and
`boundaries.csv`. The boundary figure draws resolved intervals as solid lines,
provisional intervals as dotted lines, and accepted equilibrium speed estimates
as markers. Intervals with unverified mode identity or discretization remain in
the CSV. A failed case breaks a bracket. An unresolved point between valid
observations keeps the full interval and records NaN for the estimate.

The default fit windows cover 0.50–1.20 s and 1.20–1.90 s to give early bending
transients time to decay. The pulse coefficients are 0.05 and 0.10 N to make
small, decaying responses easier to measure. Fits still require proportional
amplitudes across the two pulses. Increase `duration` and move the two
`growth_options.windows` farther apart for weak growth, or move them earlier
if saturation starts in the fitting record. Check the individual fits and their
frequencies before adjusting `frequency_band`. Increase
`plateau_time` and `lco_options.tail_time` for slowly settling, low-frequency,
or modulated responses. Check the sensitivity of the resulting LCO brackets
to plateau duration, step size, timestep and mesh before using them as physical
boundaries.

## Legacy stability sweep

Open `run_pazy_stability_sweep.jl` and edit the settings at the top:

- `pitch_angles_deg` and `airspeeds_mps`: every listed airspeed is run at every
  listed pitch angle. The reference curves do not select simulation cases.
- The flight, UVLM and structural settings use the same names as the single-case run.
- `match_wake_panel_chord=true`: use a nominal wake-row length equal to one
  chordwise bound panel; `false` keeps `time_step_chords=0.25`.
- `response_to_fit`: `:twist` or `:bending`.
- With `response_to_fit=:twist`, the time history is
  `result.tip_twist_degrees` sampled at `result.time`. The original
  `moving_block_2.jl` is preserved in `studies/reference`, and its callable
  method is `studies/src/MovingBlockDamping.jl`.
- The moving-block calculation uses the full response and the original
  automatic block-size settings. `boundary_minimum_fit_r_squared` filters
  unreliable fits when constructing the stability boundary; it does not
  change the moving-block calculation.
- `midpoint_bracket_angles_deg`: angles with saturating responses whose
  unresolved boundaries are shown at the midpoint of their speed brackets.
  The CSV keeps both bracket endpoints; the midpoint is not a damping root.
- `run_simulations=false`: read saved histories instead of solving again.

The sweep starts at full airspeed (`initial_airspeed_fraction=1.0`,
`airspeed_ramp_duration=0.0`), like the single-case example. With zero tip-pulse
amplitude, the startup response supplies the excitation for the fit.

Run the file from the IDE. It performs four steps:

1. Solve every listed angle/speed combination independently.
2. Fit the peak-trimmed full response with `moving_block_metrics`.
3. Find sign changes in growth rate between adjacent valid speeds.
4. Shade the UVLM instability region and plot the earlier AeroBeams boundaries
   as lines alongside the experimental and literature data.

Positive growth means a growing response; positive damping means decay. The first
decay-to-growth crossing is onset; growth-to-decay is restabilization. Failed or
poorly fitted cases and nearly zero growth are excluded from interpolation.
If an unresolved case separates opposite valid signs, its neighboring speeds
are recorded as an uncertainty bracket in the CSV. The shaded region uses the
last confirmed growing speed as a conservative edge. Solver failures are never
bridged.

Set the speed range wide enough to include decaying cases on both sides of the
flutter region. The AeroBeams, UM/NAST, and SHARPy curves are used only in the
comparison plot.

## Output

In `output/stability_sweep`:

- `damping_summary.csv`: growth, damping, frequency, fit quality and status.
- `uvlm_boundary.csv`: onset/offset speeds from accepted damping fits only.
- `uvlm_boundary_updated.csv`: estimated onset/offset speeds, including
  provisional crossings anchored by an unresolved positive-growth fit. Its
  `quality` and endpoint columns show which estimates need confirmation. In
  particular, `bracket_midpoint_unresolved` is only a display estimate. The
  comparison plot uses these estimates when the sweep runs.
- `uvlm_boundary_brackets.csv`: speed intervals containing an unresolved crossing.
- `pazy_stability_boundary.png` and `.svg`: the shaded UVLM region, earlier
  AeroBeams third-mode onset/offset lines, UM/NAST and SHARPy curves, and
  experimental points. The legend labels the UVLM region as estimated when it
  includes unresolved fits.
- Each `aoa3_V50`-style folder: `response.csv` with both tip histories and airspeed,
  `moving_block.png` with the selected response and fit, or `failure.txt`.

## Edit the comparison plot

The last block of `run_pazy_stability_sweep.jl` creates `comparison_plot` and
sets its colors, line styles, legend labels, axis labels, font sizes, line and
marker sizes, axis limits and legend position. Experimental points
and literature curves follow the AeroBeams plot colors and markers; the earlier
AeroBeams eigenvalue boundaries share a dark blue from `ColorSchemes.darkrainbow`
and use distinct line styles. Edit and rerun only that block in
the Julia REPL to change the figure without repeating the cases or fits.
You can also edit the existing plot directly:

```julia
xlims!(comparison_plot, (35, 75))
ylims!(comparison_plot, (2, 7))
savefig(comparison_plot, joinpath(output_directory, "pazy_stability_boundary_zoom.png"))
```

In a new Julia session, set `run_simulations=false` and run the sweep script to
rebuild the plot from saved responses without solving the cases again.

There is no automatic restart, signature checking or case-management system.
Rerunning overwrites outputs. Use a different directory to preserve a study.
When using saved histories, you are responsible for selecting the correct cases.
Old files from the previous driver are not deleted; only the outputs listed above
are maintained by this simplified version.

## Interpretation

The reference plot reads AeroBeams' **cached, tracked third-mode eigenanalysis**
and the UM/NAST, SHARPy and experimental data used by
`test/examples/PazyWingFlutterPitchRange.jl`. It does not run a new strip-theory
simulation or linearize the UVLM.

The fit interval may contain both initial growth and a later saturated response.
Inspect the response and fit yourself: use
small-amplitude oscillations about a settled deformed configuration and verify
that the same mode is being fitted at adjacent speeds. The dominant-FFT estimator
is not a mode tracker. A finite-amplitude LCO with near-zero growth does not locate
flutter of the equilibrium.

These are provisional boundaries. Refine speed spacing, mesh, timestep, wake
length and fitting window before drawing conclusions. The existing solver and
the previously reported dynamic-Jacobian issue are unchanged.

