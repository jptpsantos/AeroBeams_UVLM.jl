# Downsweep review, 23 October 2026

The inspected branch driver has no identified reversal, force-application, or
strain-export error. This is a scoped code/data audit, not proof that the
coupled aerodynamic discretization or physical model is converged. The saved
independent downsweep exhibits delayed, small oscillations; the current pulse
case changes several parameters and has not yet exported a response CSV.
No simulation settings or solver source were changed for this review.

## Which results are available

- `output/LCO_independent_sweeps/pazy_linear_LCO_downsweep.csv` now contains
  **52 -> 35 m/s in 20 s**, 42,064 samples. Its upsweep contains 35 -> 52 m/s
  in 40 s, 84,127 samples. Both have dt = 0.00047547726030002614 s.
  This is the completed unforced case, not the new pulse case.
- The older `latest_sweep_analysis.md` and `latest_sweep_window_analysis.csv`
  in that directory describe a previous 60 m/s run. Their conclusions and
  amplitudes must not be assigned to the present CSVs.
- `output/LCO_downsweep_pulse/` contains settings and a converged static
  initialization, but no completed dynamic response CSV or PNG at inspection.
  Missing output alone cannot distinguish a running, interrupted, or failed job.
- The continuous log has passed the reversal and reached approximately
  t = 48.58 s, U = 44.70 m/s; its output directory has no completed downsweep
  CSV yet. The available full-history file was written at the reversal.

## Fresh measurements of the completed independent case

The audit checks finite data, increasing time, strictly decreasing downsweep
speed, and strain reconstruction from the stored section strain and curvature.
The reconstruction error is zero for both branches. In 0.5 m/s bins, a quadratic
time trend is removed and the residual half peak-to-peak range is measured.
These are transient amplitude estimates, not settled LCO amplitudes or exact
stability boundaries. The plotted curves retain the original undetrended strain.

| Downsweep speed (m/s) | Branch time (s) | Detrended half range (microstrain) |
| --- | --- | --- |
| 49 | 3.5295 | 0.000109 |
| 48 | 4.7058 | 0.000615 |
| 47 | 5.8824 | 0.0356 |
| 46 | 7.0589 | 1.774 |
| 45 | 8.2353 | 26.06 |
| 44 | 9.4116 | 66.01 |
| 43 | 10.5884 | 51.03 |
| 42 | 11.7647 | 11.91 |
| 41 | 12.9411 | 1.972 |
| 40 | 14.1176 | 0.2004 |

The largest sampled downsweep bin is **66.07 microstrain at 43.5 m/s**;
the upsweep reaches **741.60 microstrain at 49 m/s**. The downsweep disturbance
grows by many orders of magnitude before becoming visible, then decays as the
speed falls further. This supports an excitation/history and finite-growth-time
explanation. It does not establish that a larger seed would reproduce the
experimental amplitude or lower termination speed.

## Code checks

- `PazyLinearLCOSweep.jl` sets the down branch to
  `(downsweep_initial_speed, initial_speed, downsweep_duration)`. Its ramp is
  `start + (end-start)*clamp(t/duration,0,1)`, with no reversed time axis.
- Each selected branch gets a separate converged coupled static initialization
  and developed wake. The dynamic problem copies that equilibrium. No upsweep
  oscillation or structural velocity is supplied to the independent downsweep.
- The custom speed profile overrides the startup ramp. `settling_time` sets
  the pulse onset here; it does not hold the speed constant until the pulse.
- The dynamic loop updates both the UVLM operating point and dimensional-load
  reference speed every step. Passing the peak speed as `airspeed` fixes dt;
  it does not keep the downsweep freestream at the peak speed.
- The pulse acts on `nodal_loads[1,end]`, the tip force along structural A1,
  which is the out-of-plane direction in this model. Each coupling iteration
  first replaces the load array with aerodynamic loads, then adds the pulse.
  The force is not accumulated over iterations. Boundary-condition caches are
  refreshed before the structural solve.
- All aerodynamic trials begin from the same accepted state. Wake convection
  is committed once per accepted physical step, not once per coupling iteration.
- The output uses `1e6*(gamma1 + z*kappa2 - y*kappa3)` for both branches.
  The plot uses instantaneous speed/strain pairs in their original time order.
  No absolute value, sorting, envelope substitution, or down-branch sign flip
  explains the discrepancy.

The existing short physical regression log records 151 passing independent
branch checks and 28 passing legacy sweep checks, including forced/unforced
comparisons, isolated branch execution, static-state copying, exported force,
and plotted series. That suite was reviewed, not rerun here. Fresh CSV checks
and pulse arithmetic/spectral checks were run for this audit; all passed.
Neither set of checks validates experimental flutter boundaries.

## Current pulse runner differs from the completed baseline

The current `run_pazy_LCO_downsweep_pulse.jl` selects only the downsweep:

| Parameter | Current value |
| --- | --- |
| Speed history | 52 -> 35 m/s in **10 s** |
| Deceleration | -1.7 m/s per second |
| Mesh | 5 chordwise x 16 spanwise panels |
| Requested minimum-speed wake length | 3 chords |
| Resolved wake capacity / age | 18 rows / 0.0085584 s |
| Pulse trigger | 49 m/s, t = 1.764705882 s |
| Pulse coefficient / duration | 1 N / 0.1 s |
| Actual peak force magnitude | about 0.6495 N |
| Speed interval during pulse | 49 -> 48.83 m/s |

The 10 s ramp provides half the time between any two speeds compared with the
20 s baseline. The coarser mesh and shorter wake also prevent treating this as
a pulse-only comparison with the original 10 x 20, 5-chord configuration.

### Pulse duration can suppress excitation near the observed frequency

For phase s = (t-t_start)/T in (0,1), the prescribed waveform is

`F(s) = A*sin(2*pi*s)*sin(pi*s)^2`

`     = A*(0.5*sin(2*pi*s) - 0.25*sin(4*pi*s))`.

Its net impulse is zero. A zero net impulse does not imply no vibration, but
this waveform has exact Fourier zeros at f = n/T for integers n >= 3.
For T = 0.1 s, **30 Hz is the first such zero**. The earlier 45 m/s diagnostics
reported dominant response frequencies around 28-31 Hz.

Numerical quadrature of `abs(integral F(t)*exp(-i*2*pi*f*t) dt)` gives:

| Duration | 28 Hz | 30 Hz | 31 Hz |
| --- | --- | --- | --- |
| 0.10 s | 0.0010685 | approximately zero | 0.0003055 |
| 0.04 s | 0.0100656 | 0.0099662 | 0.0098761 |

Units are N s for a 1 N coefficient. The 0.04 s diagnostic used a 0.1 N
coefficient, so its actual values are one tenth of the second row. These are
forcing spectral magnitudes, not calculated modal amplitudes. The frequency
at 49 m/s with the current mesh is not established; damping, coupling, and
time-varying speed can alter the response. Nevertheless, increasing the
coefficient of the 0.1 s pulse alone does not remove its 30 Hz spectral zero.
The current pulse spans about 210 timesteps, so under-sampling of the pulse
is not the evident problem.

## Causes to separate before interpreting the reference mismatch

1. **Initialization and growth time.** The independent down branch begins
   essentially quiet. Its small residual grows late, rather than starting
   from a developed oscillatory state. The faster current ramp reduces growth
   time further. Large-amplitude branch selection is a possible additional
   effect, but coexistence of solutions has not been demonstrated.
2. **Unresolved numerical sensitivity.** In the existing constant-45 m/s
   diagnostic, halving dt while retaining nominal wake age changed initial
   strain from -1540.77 to -1383.48 microstrain and early fitted growth from
   about 2.05 to 0.297 per second. Last-window amplitudes were about 1356.5
   versus 1.06 microstrain, with the refined case still growing. A passing
   coupling residual is insufficient to establish mesh/timestep convergence.
3. **Forcing spectrum.** The new 0.1 s force is not a broadband impulse and
   can weakly excite the relevant mode. This applies to the pending pulse
   case; it cannot explain the already completed unforced result.
4. **Reference protocol and model.** The earlier paper review records that
   the experimental 5-degree downsweep follows a 7-degree upsweep and a root
   pitch change to 5 degrees. A fresh 5-degree static start, or a fixed-5-degree
   round trip, does not reproduce that history. Gravity, gauge location,
   and aerodynamic formulation also differ. See `PAZY_LCO_PAPER_REVIEW.md`.
   A similar-looking upsweep alone does not resolve these differences.

## Controlled next comparisons

First use the original 20 s downsweep, 10 x 20 mesh and 5-chord wake, with
`fixed_time_step=0.00047547726030002614`, to match the completed independent
baseline. At 49 m/s compare no pulse with coefficients 0.01 and 0.1 N and
duration 0.04 s, in separate output directories. Retain the current 1 N,
0.1 s case as a different experiment. The pulse onset for a 20 s ramp is
3.529411765 s, calculated automatically by the runner.

Next compare a 40 s downsweep with the same dt, mesh, wake rule and pulse
trigger speed, to isolate rate effects. Use held-speed tests near 45 m/s to
measure growth/decay after forcing and determine whether the amplitude settles.
Do not label a pulse transient as an LCO just because it forms a band in the
strain-speed plot.

Resolve timestep and mesh/wake sensitivity before claiming quantitative
agreement. When reducing dt, increase wake rows proportionally to retain wake
age; test wake length and core/shedding controls separately. A controlled
continuous-versus-independent comparison must also match dt and wake capacity:
the currently running continuous case uses a different dt and wake capacity
from the completed independent case.

Reproduce the fresh read-only audit from the repository root with:

`julia --startup-file=no --project=. tmp/audit_pazy_downsweep.jl`

Captured output: `tmp/pazy_downsweep_audit_2026-10-23.log`.

## Later follow-up: pulse run completed, upsweep appears flat

The pulse-folder availability statement above describes the earlier inspection.
At the later inspection, the folder contains completed up/down CSVs and a
combined figure. The saved run differs from the currently edited runner:

| Parameter | Completed pulse-folder results | Runner at follow-up inspection |
| --- | --- | --- |
| Upsweep duration | 20 s | 40 s |
| Downsweep duration | 10 s | 20 s |
| Down pulse coefficient | 1 N | 0.1 N |
| Down pulse duration | 0.1 s | 0.1 s |
| Mesh / wake setting | 5 x 16 / 3 chords | 5 x 16 / 3 chords |

Both branches were run and plotted. The upsweep force column is identically
zero; detrended half ranges near 49, 50 and 51 m/s are approximately 0.000104,
0.000109 and 0.000121 microstrain. The largest upsweep residual is an initial
transient of about 0.267 microstrain near 35.5 m/s. The forced downsweep peaks
near 42.5 m/s at about 806.43 microstrain. Thus the missing upsweep band is a
computed response, not a missing plot series or force accidentally applied to
the wrong branch. Reproduction: `julia --startup-file=no --project=.
tmp/analyze_pazy_pulse_branches.jl` (one command).

Compared with the original oscillating upsweep, the completed new upsweep has
half the duration, a coarser mesh and a shorter wake. Its starting strain is
-897.74 microstrain versus -834.69 in the original run at the same 35 m/s.
This approximately 7.6% difference exists before either branch moves and shows
that the two aerodynamic discretizations already give different equilibria.
The relative contributions of changed growth rate, initial seed and traversal
time cannot be determined by comparing these two runs alone.

The driver now accepts optional `upsweep_pulse_amplitude`,
`upsweep_pulse_start`, and `upsweep_pulse_duration` in addition to the existing
downsweep controls. Both default to zero amplitude. The separate pulse runner
adds `pulse_branches=:both/:upsweep/:downsweep/:none`; it retains the existing
downsweep-only selection unless edited. Ascending diagnostic controls are
45 m/s, 0.1 N coefficient and 0.04 s duration. `:both` uses the separate
`output/LCO_both_pulses` folder. Existing user-selected numerical values and
descending pulse parameters were retained.

Applying a pulse to each branch tests controlled excitation, and can make
small-seed growth visible. It does not guarantee a sustained LCO, numerical
convergence, or reference agreement. Match discretization and ramp duration
before comparing seeded/unseeded runs. After the force ends, distinguish
growing, sustained and decaying motion, and check dependence on pulse amplitude.

Validation after adding the ascending controls: 164 independent-sweep checks
and 28 legacy-sweep checks passed (192 total). This includes a forced upsweep
that differs from its unforced control, an unchanged independently forced
downsweep when the upsweep pulse is enabled, pulse timing and exported forces,
and rejection of invalid pulse durations. Log: `tmp/pazy_LCO_both_pulses_test.log`.
No full production simulation was launched for this change.
