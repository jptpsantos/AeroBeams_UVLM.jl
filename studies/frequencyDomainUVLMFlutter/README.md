# Flutter of the Pazy wing with frequency-domain UVLM aerodynamics

Flutter of the undeformed Pazy wing (zero root pitch angle, no gravity) with the aerodynamic forces of a frequency-domain unsteady vortex-lattice method (UVLM), in the basis of the vibration modes of AeroBeams (the same structural model of `create_Pazy`). It is a reference for the strip theory of AeroBeams, which has the loss of lift of the finite wing (tip loss function) as a correction, as the three-dimensional unsteady potential flow includes the effects of the finite wing on the unsteady loads.

- `modalData.jl`: modes of the model without airflow (frequencies, generalized masses from the strain energy of the mode, bending displacement and twist of the reference axis along the wing).
- `unsteadyVLM.jl`: vortex lattice (straight or swept planar wing, wall as a symmetry plane, planar wake of finite length), generalized aerodynamic forces for harmonic motion in the modes, strip theory (Theodorsen/Jones) for comparison, and the p-k solution.
- `flutterUVLM.jl`: script with the options at its top (sweep angle, tip mass configuration, panels, wake length), which prints the flutter onset and offset speeds of the UVLM and, optionally, of the strip theory of AeroBeams with and without the finite-wing lag correction (`Indicial(finiteWingLag=true)`).

## Validation

- Steady lift slope of the planform, `CL/(2πα)` with a wake of 10 chords and 5 × 15 panels: 0.7946 (Λ=0°), 0.7906 (10°), 0.7560 (20°), equal to the 4 digits to the UVLM data of the sweeps in `test/referenceData/PazyVLM` at the lowest airspeed and pitch angle.
- Unsteady lift and moment of the sections of a wing of very high aspect ratio in plunge and pitch: within 1% (magnitude) and 1° (phase) of Theodorsen's solution (with Jones' approximation of the Wagner function) for `k` from 0.05 to 0.3.
- Flutter speed of the straight wing at zero pitch angle, plain structure: 69.9 m/s (34.2 Hz, torsion branch) and 95 m/s (17 Hz, second bending branch) for 8 × 30 panels and a wake of 20 chords (5 × 15 panels: 70.6 m/s; 16 × 30: 69.2 m/s; wake of 2 chords: 67.7 m/s). The time-domain coupling of the UVLM of the fork `AeroBeams_UVLM.jl` with the same structure (10 × 20 panels, time step of 1/4 chord) has a neutral stability at 65.5 m/s and 34.6 Hz, Sharpy 64.5 m/s and ZAERO about 67 m/s.

## Zero pitch angle flutter, UVLM and AeroBeams (flat plate, `"VLM_tan_fixed"` tip loss; onset in m/s, frequency in Hz)

| Configuration | UVLM (modes of AeroBeams) | Experiment | Strip theory | Strip theory with `finiteWingLag` |
|---|---|---|---|---|
| S0 | 69.9 (34.2) | - | 99.5 (14.9, bending branch; no torsion flutter) | 71.4 (34.0) |
| S10 LE | 57.4 (33.7) | 56.0 (32.8) | 63.6 (33.0) | 59.8 (33.5) |
| S10 TE | 76.9 (17.4) | 78.1 (18.8) | 80.3 (15.2) | 76.2 (18.0) |
| S20 LE | 50.5 (32.1) | 52.3 (31.8) | 56.8 (31.6) | 54.6 (31.9) |
| S20 TE | 75.1 (16.6) | 76.7 (16.3), at 3° | 75.8 (14.2) | 72.7 (17.2) |

(Strip theory: flat plate, `"VLM"` tip loss of the undeformed wing, plain structure of the straight wing, rod modeled as in the former `dev/Pazy/S10/InertiaRod` and `dev/Pazy/S20/InertiaRod` scripts (now those of `dev/Pazy/S10` and `dev/Pazy/S20` with `beamWingtipRod = false`), i.e., lumped with the ad hoc inertia correction for the S10 and as a beam for the S20.)

The strip theory scales the whole circulatory response with the tip loss factor ϖ, which models the steady loss of lift, so that the initial (impulsive) circulatory load of the Wagner function, which is that of the two-dimensional wing, becomes ϖ times its value. `Indicial(finiteWingLag=true)` scales the amplitudes of the exponentials of the lag by `κ = (1 - (1-ΣA)/ϖ)/ΣA`, which keeps the initial load, with no parameters to calibrate.
