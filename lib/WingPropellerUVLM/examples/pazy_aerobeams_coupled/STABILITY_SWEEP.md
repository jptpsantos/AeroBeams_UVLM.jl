# Simple Pazy stability sweep

Open `run_pazy_stability_sweep.jl` and edit the settings at the top:

- `pitch_angles_deg` and `airspeeds_mps`: angles and increasing speeds to simulate.
- The flight, UVLM and structural settings use the same names as the single-case run.
- `response_to_fit`: `:twist` or `:bending`.
- `fit_window_s`: `(start_s, end_s)` for the single damping fit.
- `block_duration_s` and `minimum_fit_r_squared`: moving-block fit settings.
- `run_simulations=false`: read saved histories instead of solving again.

The sweep starts at full airspeed (`initial_airspeed_fraction=1.0`,
`airspeed_ramp_duration=0.0`), like the single-case example. With zero tip-pulse
amplitude, the startup response supplies the excitation for the fit.

Run the file from the IDE. It performs four steps:

1. Solve each angle/speed case independently.
2. Fit the selected window with the existing `moving_block_metrics`.
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

## Output

In `output/stability_sweep`:

- `damping_summary.csv`: growth, damping, frequency, fit quality and status.
- `uvlm_boundary.csv`: interpolated onset/offset speeds.
- `uvlm_boundary_brackets.csv`: speed intervals containing an unresolved crossing.
- `pazy_stability_boundary.png` and `.svg`: the shaded UVLM region, earlier
  AeroBeams third-mode onset/offset lines, UM/NAST and SHARPy curves, and
  experimental points. Unresolved crossings are not drawn as exact boundaries.
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

