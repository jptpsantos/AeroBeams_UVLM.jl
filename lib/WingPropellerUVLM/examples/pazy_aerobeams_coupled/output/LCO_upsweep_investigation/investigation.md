# Why the unforced upsweep differs in the pulse experiment

This review compares actual saved histories, not current editor values. No
production simulation was run or changed. Measurements and source SHA256 hashes
are in `measurements.txt`; the reproducible analysis is
`tmp/investigate_pazy_upsweep.jl`, run from the repository root. The amplitude
figure uses quadratic-detrended half ranges in 0.5 m/s bins. Very small residuals
include numerical noise and trend error; they are not measured modal amplitudes.

## Findings

The downsweep pulse does not disable the upsweep. The shared driver completes
the upsweep, writes its history, and only then initializes the independent
downsweep. Each branch receives its own pulse parameters. Previously completed
regressions compare the unforced upsweep with and without a downsweep pulse at
identical settings and verify matching displacement and strain histories
(`test_pazy_LCO.jl`, forced/quiet comparison). The current source has the same
separation. These short tests verify branch isolation, not full LCO accuracy.

| Quantity | Original unforced upsweep | Unforced upsweep in pulse folder |
| --- | --- | --- |
| Saved duration | 40 s | 20 s |
| Ramp rate | 0.425 m/s per second | 0.85 m/s per second |
| Mesh used | 10 x 20 | 5 x 16 |
| Wake capacity | 30 rows | 18 rows |
| Wake length setting at minimum speed | 5 chords | 3 chords |
| dt | 0.000475477260300 s | 0.000475465956638 s |
| Initial strain at 35 m/s | -834.690 microstrain | -897.744 microstrain |
| Initial tip displacement | 0.0653325 m | 0.0702973 m |
| Initial tip twist | 0.901115 deg | 0.840604 deg |
| Half range near 47 m/s | 0.005959 microstrain | 0.00008293 microstrain |
| Half range near 48 m/s | 17.85 microstrain | 0.00009225 microstrain |
| Half range near 49 m/s | 741.60 microstrain | 0.00010399 microstrain |

Durations, dt, section states and amplitudes are read from the CSVs. Original
mesh/wake settings are from its runner/driver defaults and the earlier audit;
the pulse-run mesh/wake were captured in the earlier saved-settings inspection.
The old CSV format does not embed mesh metadata. The source hashes match those
used in the previous audit/selected comparison where available.

The starting strain differs by about 7.6% in magnitude. This difference exists
before the ramp or pulse, so growth time alone does not explain every change:
the aerodynamic discretizations already produce different static equilibria.
The timestep difference here is only about 0.0024%; this is not the earlier
half-timestep experiment. Mesh and wake were both changed, so neither can be
identified individually as the cause from this comparison.

## Why an unforced run can oscillate

Zero imposed tip pulse does not mean an exactly undisturbed trajectory. The
speed ramp moves the equilibrium, static convergence is finite, and numerical
disturbances exist. If a relevant perturbation is unstable, it can amplify.
The original upsweep's rapid growth near 47-49 m/s is consistent with that
mechanism; the source of the initial modal seed was not separately measured.

For a small perturbation with local growth rate sigma(U), a useful diagnostic
approximation is A(t) = A0 * exp(integral sigma(U(t)) dt). A twice-as-fast ramp
halves the accumulated exponent over the same speed interval if sigma(U) is
unchanged. It does not simply halve the final amplitude. The new mesh/wake can
also change sigma and the initial seed, so this relation is an explanation of
sensitivity, not a quantitative prediction for these different discretizations.

## Why the pulsed downsweep grows

The saved down branch lasts 10 s and includes a bipolar force with coefficient
approximately 1 N, peak magnitude 0.6495 N, applied over about t = 1.765-1.865 s
at 49-48.83 m/s. These are the values recovered from its force column, not the
0.1 N coefficient currently shown in the editor.

After the pulse has ended, the detrended half range grows from 4.89 microstrain
at 48 m/s to 18.82 at 47, 327.45 at 45, and 806.43 near 42.5 m/s. Thus the red
band is not confined to the time the external force acts. It is consistent
with amplification of the seeded motion during descent. The varying-speed
history alone does not establish a settled LCO or rule out finite-amplitude
branch selection. A fixed-speed study and convergence checks are needed for
those conclusions.

The down branch receives a controlled disturbance near the speed range of
interest. The up branch starts at 35 m/s without that disturbance and traverses
a different history. These conditions need not produce equal visible growth.

## The current folder mixes new metadata with old completed histories

At this inspection, `LCO_downsweep_pulse/pazy_linear_LCO_settings.txt` and its
upsweep-initialization CSV were updated at approximately 14:32 on 23 October.
The settings request 40/20 s and a 0.1 N down pulse. However, the response CSVs
are still the 20/10 s files saved at approximately 12:34 and 13:47; the combined
PNG is also from 13:47. Therefore that PNG is not evidence of the new 40 s
upsweep's result. The metadata was written by a newer invocation before new
completed response files appeared. This does not by itself establish whether
that invocation is still running, stopped, or failed.

## Comparisons that isolate the causes

Keep the upsweep unforced and retain the same explicit dt, solver controls and
gauge. Compare this matrix in separate output directories:

| Case | Mesh / wake setting | Up duration | Status |
| --- | --- | --- | --- |
| A | 10 x 20 / 5 chords | 40 s | Original completed oscillating case |
| B | 10 x 20 / 5 chords | 20 s | Not performed in this review; isolates rate against A |
| C | 5 x 16 / 3 chords | 40 s | Current requested setup, but no matching completed CSV at inspection |
| D | 5 x 16 / 3 chords | 20 s | Completed nearly flat case |

Use dt = 0.00047546595663750473 s for a strictly matched new matrix; it aligns
both 20 and 40 s durations. The existing A differs slightly in dt, as reported
above. Match wake capacity within comparisons intended to change only duration.
Changing mesh and wake separately would further distinguish their effects.

For a direct software-control comparison, use the original 10 x 20, 5-chord,
40/20 s setup and original dt = 0.00047547726030002614 s, leaving the up pulse
disabled. Adding only a down pulse should retain the original up calculation.
If a full matched repeat differs materially, investigate code version,
initialization, solver state and numerical reproducibility before attributing
the difference to physical hysteresis. Merely adding an up pulse would not
isolate the reason the previous unforced up calculation changed.
