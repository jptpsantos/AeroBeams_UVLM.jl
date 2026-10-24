# Pazy constant-speed diagnostic: completed results

All four cases completed at 45 m/s, 5 degrees root incidence, about 5 s physical time, with the same 10 x 20 surface mesh. No new simulations were run for this analysis.

**Main result: the present response is strongly timestep-sensitive. The small downsweep cannot be attributed solely to missing excitation, and the current results do not establish numerical convergence.**

| Case | dt [ms] | Initial strain [microstrain] | Final-window half range [microstrain] | Early RMS growth [1/s] | Tail RMS growth [1/s] | Late dominant frequency [Hz] |
|---|---:|---:|---:|---:|---:|---:|
| Pulse baseline | 0.475477 | -1540.77 | 1356.506 | 2.049 | 0.231 | 27.99 |
| Unforced | 0.475477 | -1540.77 | 109.266 | 2.262 | 1.964 | 29.99 |
| Half timestep | 0.237739 | -1383.48 | 1.064 | 0.297 | 0.294 | 30.99 |
| Double wake | 0.475477 | -1559.99 | 1308.948 | 2.228 | 0.000 | 27.99 |

Half range and RMS come from a linear-detrended 0.25 s window; the final complete window is 4.64-4.89 s. The early growth fit uses window centers between 0.75 and 2.5 s; it is a finite-response estimate, not an eigenvalue. Tail growth uses the last four windows and may include the approach to saturation. Frequencies are Hann FFT peaks over the last 1 s, with about 1 Hz frequency-bin spacing.

## Interpretation

1. **Excitation and finite time matter.** The unforced response grows from very small disturbances to about 109 microstrain by the final window. The pulse case reaches about 1357 microstrain over the same elapsed time. Both grow with the baseline timestep; the unforced case has not reached a plateau. Residual equilibrium mismatch/roundoff supplies disturbances in the nominally unforced run. This supports delayed growth as a contributor to the ramp result, but does not quantify its contribution to the downsweep.

2. **Timestep sensitivity dominates this comparison.** Halving dt reduces the final-window half range by about 99.92%, while reducing the magnitude of the initial strain by about 10.21%. Initial tip displacement changes from about 0.11990 m to 0.10818 m. The refined response is still growing slowly (about 0.29 1/s), so 1.06 microstrain is a transient value, not its settled LCO amplitude. Because the equilibrium changes too, the dynamic runs are linearizing/evolving about different configurations.

3. **A longer wake has a smaller effect here, but wake convergence is unproven.** Doubling retained wake length from 6.49 to 12.98 chords changes the final-window half range by about -3.51% and initial strain magnitude by about +1.25%. It also changes the growth trajectory. The double-wake case is nearly stationary over its final four RMS windows; this short interval alone does not establish a unique, fully settled LCO.

4. **Iteration convergence passed; discretization convergence did not.** Recorded geometry and load residuals remain within 1e-4, and the static initializations report convergence. Those tests measure iteration consistency at each chosen discretization. They do not remove the differences between timestep/wake discretizations.

## What this means for the sweeps

The old downsweep result remains a computed transient for its selected settings. These tests do not establish a plotting bug, a unique solver bug, or agreement with the reference's LCO branch. Raising the downsweep starting speed also reduces the shared timestep in the current driver, so a comparison against the old 52 m/s case changes both starting speed and numerical resolution. To isolate starting-speed effects, rerun both configurations at the same dt and physical wake length after identifying suitable discretization settings.

## Recommended order of further work

1. First investigate fixed-geometry aerodynamic/steady-load sensitivity to dt, without structural motion or a pulse. Check circulation and total loads together with wake geometry. The roughly 10% equilibrium-strain shift already exists before the pulse.

2. Inspect near-trailing-edge wake spacing and finite-core treatment. In the current code the shedding offset is eta*V*dt. At eta=0.1 and 45 m/s this is about 2.14 mm for baseline and 1.07 mm at half dt, compared with the fixed 1 mm core radius. This is a concrete numerical sensitivity to investigate, not a confirmed cause.

3. Then use an additional timestep refinement and controlled wake/core/mesh checks to establish trends in steady loads and small-disturbance growth. Keep physical wake length and pulse fixed; tighter coupling/Newton tolerances may also need checking for small-amplitude responses.

4. Once those quantities are sufficiently insensitive, extend the constant-speed runs as needed for amplitude stationarity and repeat speed ramps with controlled excitation/initial conditions. A 45 m/s check alone cannot identify the full flutter boundary or reproduce the reference diagram.

## Files and checks

[Comparison figure](analysis_comparison.png), [late response traces](analysis_late_responses.png), [computed metrics](analysis_metrics.csv), [original summary](diagnostic_summary.csv).

The analysis checked finite CSV values, uniform time grids, constant 45 m/s speed, reconstructed surface strain, and recorded coupling residuals. Baseline/unforced initial strains match exactly. Forced cases have matching pulse values at common sample times. Baseline and half dt retain the same nominal physical wake age/length (30 and 60 rows respectively); double wake uses 60 rows at baseline dt.

