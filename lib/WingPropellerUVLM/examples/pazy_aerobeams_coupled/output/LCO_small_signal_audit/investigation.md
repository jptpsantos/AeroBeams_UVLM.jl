# Detailed implementation and small-disturbance investigation

This investigation tests mechanisms beyond differences between the two input
files. It does not assume that the changed settings are the sole explanation.
The production runners, structural solver and UVLM solver were not modified.

**Final finding, 24 October:** the latest completed unforced upsweeps in the
original and pulse-run directories are identical in every sample and recover
the original oscillation band. The displayed combined pulse-run PNG is older
than the new CSV and still represents the earlier calculation. The nine
constant-speed diagnostic runs find no small-disturbance cutoff, but do find
about 14% coupling-tolerance error in the tested short-time normalized RMS.
Details and limitations follow. A fresh reference-format view of the verified
up response is `latest_unforced_upsweep.png`; the controlled tolerance figure
is `normalized_pulse_response.png`.

## Original driver versus current driver

`tmp/pazy_driver_call_audit.jl` loads the original driver from Git HEAD
(`e687ef8d`) and the current driver into separate modules, replacing only the
static/dynamic solvers with recording stubs. With the original physical and
numerical inputs and a downsweep-only pulse in the current driver:

- All static initialization keyword values for the upsweep match.
- All 38 dynamic solver keywords match, excluding direct function-object
  identity for the speed closure; the closure values match at 10,001 times.
- Both upsweeps receive a zero force coefficient; only the current downsweep
  receives the specified nonzero coefficient.

This is a direct orchestration comparison, not a new full-length simulation.
The existing physical branch-isolation regressions separately verify matching
upsweep displacement/strain when only the downsweep pulse is changed.
Log: `tmp/pazy_driver_call_audit.log`.

The production coupling differs from that Git revision only by the optional
accepted-step callback. Neither production linear sweep activates it. The
inspected structural solver and UVLM sources have no working-tree diff from
that revision. This comparison does not identify what methods an existing
interactive Julia session loaded previously.

## State and load paths reviewed

- Each independent branch constructs its own static equilibrium and copies it
  into a new dynamic problem. There is no branch-to-branch velocity/wake reuse.
- `create_DynamicProblem` receives the equilibrium state vector; it does not
  call the unloaded dynamic initialization in this path.
- The speed profile is sampled on the fixed time vector. Both UVLM freestream
  and dimensional load reference speed are updated at each accepted step.
- The external tip force is added to freshly assigned aerodynamic loads on
  each coupling iteration, rather than accumulating over iterations. Active
  model boundary-condition caches are updated before structural solves.
- Equivalent structural rates are formed once per physical time step. Inner
  structural solves retain those rates while updating the trial state.
- UVLM trials copy the same accepted aerodynamic history. Wake convection and
  the physical clock advance once, at commit.
- The small-signal runs below assert that the source equilibrium structural
  vector and circulation remain unchanged after repeated dynamic runs.
- Newton applies a computed state increment before testing convergence; the
  inspected convergence condition does not explicitly zero a small increment.
- Both branches use the same strain reconstruction and the plot uses saved
  samples directly. No threshold in the plot hides the response.

No pulse-branch logic error was identified in these paths. Coupling and temporal
discretization errors can still affect the physical growth rates.

## Controlled constant-speed experiment

All cases use the same serialized coupled equilibrium at 45 m/s, 5 degrees,
5 x 16 UVLM panels, 18 wake rows, dt = 0.00047547726030002614 s, core radius
0.001 m, and shedding fraction 0.1. Duration is approximately 0.5002 s.
The bipolar force acts from 0.01 to 0.05 s. Coefficients are 0, 0.00001 and
0.01 N, a factor of 1,000 between the two forced cases. Each tolerance setting
has its own matching unforced trajectory from the identical equilibrium.

Measure `delta_strain = forced_strain - matching_unforced_strain`. This removes
the shared equilibrium/startup transient. Dividing by force coefficient checks
whether a tiny disturbance is treated differently from a larger one. Raw
histories, coupling iterations, residuals and maximum section speed are saved.

| Coupling tolerances | Newton tolerances | Pulse coefficient (N) | RMS delta strain / coefficient, 0.39-0.49 s (microstrain/N) | Mean coupling iterations in that window |
| --- | --- | --- | --- | --- |
| 1e-4 | 1e-5 | 1e-5 | 10.40689 | 2.00 |
| 1e-4 | 1e-5 | 0.01 | 10.40733 | 2.00 |
| 1e-7 | 1e-9 | 1e-5 | 10.14664 | 3.00 |
| 1e-7 | 1e-9 | 0.01 | 9.16297 | 7.39 |
| 1e-9 | 1e-11 | 1e-5 | 9.13428 | 9.31 |
| 1e-9 | 1e-11 | 0.01 | 9.12899 | 14.07 |

At the production tolerances, the normalized last-window RMS agrees to about
0.0042% between the two pulse sizes. Thus, in this test, a tiny disturbance is
not discarded or heavily suppressed relative to a larger one. Both use two
coupling iterations. This contradicts a simple small-amplitude cutoff as the
explanation for the absent band.

The intermediate 1e-7 tolerance changes the last-window response by about 2.5%
for the small pulse and 12.0% for the larger pulse. At 1e-9 coupling / 1e-11
Newton tolerances, normalized responses agree to about 0.058% in the last
window, or 0.0635% in a relative L2 comparison over 0.06-0.49 s. The small and
larger normalized waveforms differ by 0.00577% at the default tolerances and
8.45% at the intermediate tolerance. The intermediate adaptive stopping rule
was not sufficiently tight for the smaller perturbation. The mesh, dt,
equilibrium and forcing remain identical across all nine runs.

The default last-window RMS is approximately 14.0% higher than the most refined
result. The larger-pulse RMS changes by only about 0.37% from the intermediate
to the final refinement. This identifies a real numerical coupling error;
it does not establish convergence with respect to mesh, wake, or time step.
The tests do not support an implementation that rejects a tiny disturbance
while admitting the larger disturbance at the default tolerances.

A modest short-time difference can accumulate over a long speed sweep. These
0.5 s tests neither establish a converged flutter boundary nor prove that the
coupling tolerances are immaterial. Conversely, the force-proportional response
does not establish that the resulting discrete growth rate is physically exact.

Scripts/logs: `tmp/pazy_small_signal_audit.jl`,
`tmp/pazy_small_signal_audit.log`, `tmp/pazy_small_signal_refinement.log`.
The `--refine` option reuses `equilibrium.jls` and writes separate files.

## New 40-second unforced upsweep result

While this investigation continued, a new completed upsweep CSV appeared in
the pulse output directory. It spans 35 -> 52 m/s in 40 s and contains zero
tip force. Its initial strain is -897.746 microstrain, matching the earlier
5 x 16 / 18-row initialization, not the original -834.690 microstrain state.

The detrended half range is now 0.01664 microstrain near 47 m/s, 4.637 near
48 m/s, 41.997 near 49 m/s, and 11.463 near 50 m/s. The unforced upsweep is
therefore no longer flat. This is new evidence that traversal time matters,
although the amplitude remains much lower than the original 741.60 microstrain.

The response has been preserved as `upsweep_40s_snapshot.csv`, SHA256
`853368fd8a0ecf37ca942c4e1324b59a8c4c62739a8af0e3496004f6f407480e`.
The pulse folder's settings/initialization were subsequently overwritten by
another invocation indicating 10 x 20 / 30 rows, while this response CSV still
has the earlier initial state. Current editor values and that mutable settings
file must not be used as proof of this completed CSV's configuration.

## Limits of attribution

The controlled driver comparison rules out a changed upsweep call caused by
adding a downsweep pulse at matched inputs. The small-signal runs reject an
obvious finite-amplitude cutoff at the tested speed and pulse sizes. They do
identify coupling-tolerance dependence that warrants convergence assessment.
The newest 40 s unforced result also shows that the previously absent band was
not an unconditional consequence of using the pulse entry point.

The exact contributions of mesh, wake, timestep, coupling error, initialization
residual and available growth time over the full speed range are not yet
separated. Adding a pulse to both branches would not perform that separation.

## Completed matched production histories, inspected 24 October

Both saved upsweep files were updated by later invocations:

- `LCO_downsweep_pulse/pazy_linear_LCO_upsweep.csv`, saved at approximately
  06:49 local time on 24 October; its settings specify a 40 s, 10 x 20,
  30-row case with zero up pulse and nonzero down pulse.
- `LCO_independent_sweeps/pazy_linear_LCO_upsweep.csv`, saved at approximately
  08:21; its settings specify the matching unforced 40 s case, including zero
  down pulse. Their separate initialization timestamps are approximately
  17:39 and 19:10 on 23 October.

The two current files have **identical SHA256 hashes**,
`e1d6264a7859bed953d4e0bc9126ae0c5d006f3a1d7fbc5bda8edc38a9e4a2d2`.
All 84,127 samples and all ten columns compare exactly, including time, speed,
axial strain, both curvatures, reconstructed strain, displacement, twist,
wingtip incidence proxy, and identically zero up force. Maximum and RMS
differences are zero in every column.

The current files therefore supersede the earlier file comparison; they are
not the older nine-column original file with hash beginning b493ebf1. The
matched strain history recovers a half range of 741.596 microstrain near
49 m/s. The 40 s coarse-mesh history with about 42 microstrain remains preserved
in `upsweep_40s_snapshot.csv`, so the changed file does not erase that evidence.

Reproduction: `julia --startup-file=no --project=.
tmp/compare_matched_pazy_upsweeps.jl` (one command). Results and a frozen current
history are `matched_upsweep_comparison.txt` and `matched_upsweep_snapshot.csv`.

This matched result is stronger evidence than comparing input settings alone:
the original entry point and the downsweep-pulse entry point now give the same
saved up response at matching configuration. The presence of the down pulse
does not explain the earlier loss of the up band. The earlier different runs
still do not isolate mesh, wake and rate individually. The coupling refinement
shows an additional numerical error that merits attention before quantitative
flutter/LCO validation, without demonstrating a pulse-specific software bug.

The existing combined PNG in the pulse directory was last saved on 23 October,
before these new upsweep CSVs. It still shows the earlier nearly flat up branch.
The new down CSV has not replaced the old 10 s history at this inspection, so
that directory is not yet a completed matching up/down pair for the latest run.
