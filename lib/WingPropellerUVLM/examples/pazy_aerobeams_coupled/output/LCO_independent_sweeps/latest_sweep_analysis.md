# Latest independent sweeps: 60 m/s downsweep

The completed files specify 35 -> 52 m/s in 40 s and an independently initialized
60 -> 35 m/s downsweep in 29.4117647 s. Both use dt = 0.00041207376120325537 s,
5 degrees root incidence, and the same strain reconstruction. Both static
initializations report convergence.

## Observations from the raw histories

- Upsweep: 97,071 samples; peak detrended half range about 1,015.39 microstrain
  in the 48.5 m/s bin.
- Downsweep: 71,376 samples; strictly decreasing speed, ending at 35 m/s. Strain
  matches the exported axial strain/curvature reconstruction. Around 43, 45,
  49, and 51 m/s, detrended half ranges are approximately 0.000346, 0.000184,
  0.000103, and 0.000129 microstrain, respectively. These are minute residual
  variations, not a developed LCO band.
- The largest downsweep residual bin is about 0.134 microstrain near 59.5 m/s,
  during the initial transient. These bin statistics remove a quadratic trend
  over 0.5 m/s; they are not frequency-filtered LCO amplitudes and can include
  remaining slow-trend error and numerical noise.
- The 60 m/s branch takes (60-52)/0.85 = 9.4118 s just to reach 52 m/s.

## Interpretation and limits

The plotted straight red curve represents the computed history accurately.
The driver resets this branch to a coupled static equilibrium with zero starting
structural velocity, without inheriting the previous oscillation, circulation,
or wake. Beginning higher adds a substantial interval before the lower-speed
oscillatory region seen in the upsweep. Decay of small disturbances during that
interval, followed by delayed growth from a very small seed, is consistent with
this result. The files alone do not prove the damping or its exact boundary.

The previous 45 m/s diagnostic supports excitation/growth-time sensitivity at
its baseline discretization, but also demonstrated strong timestep sensitivity:
halving dt changed both equilibrium strain and the measured growth rate. Thus
initial-condition history is a testable explanation, not a uniquely established
cause or a guarantee that a continuous sweep will reproduce the reference.

## Continuous test selected by the user

The new single integration reverses immediately at 52 m/s: 35 -> 52 over 40 s,
then 52 -> 35 over 20 s. It retains beam displacement, velocity, circulation,
and wake, and has no repeated static solve, pulse, or preparation hold.
It uses the exact dt and derived 35-row wake capacity of these independent
files, retaining the same mesh and solver tolerances. Results go to
`../LCO_continuous_sweeps`, preserving this experiment.

For a controlled reset/no-reset comparison, also use independent mode with
the downsweep starting at 52 m/s and the same fixed dt. Comparing the continuous
52 m/s reversal directly against this 60 m/s run changes both history and the
upper endpoint; that difference must be kept explicit. Neither comparison
removes the unresolved discretization sensitivity.

The numeric bin results are in `latest_sweep_window_analysis.csv`.
