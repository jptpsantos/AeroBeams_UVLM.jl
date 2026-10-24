# Pazy LCO paper comparison and proposed investigation

Review date: 22 October 2026. This review combines the supplied papers, the current Julia implementation, and the completed diagnostic results. Recommendations and proposed acceptance thresholds below are engineering proposals, not additional instructions or specifications found in the papers. No production solver settings were changed and no additional aeroelastic simulations were launched for this review.

## Sources and comparison targets

1. L. G. P. dos Santos, C. E. S. Cesnik, and F. D. Marques, *Effects of Separated-Flow Nonlinearities on the Flutter Behavior of the Pazy Wing*, SciTech 2024, DOI: https://doi.org/10.2514/6.2024-0829. Supplied file: `SciTech_2024___Luiz__paper_.pdf`, 14 PDF pages. Relevant locations: Sections III-IV (pp. 2-5), Section V and Fig. 3/Table 2 (pp. 5-7), stepped bifurcation procedure (pp. 9-10), interpretation of stall and sweep rates (pp. 12-13).
2. A. Drachinsky, O. Avin, D. E. Raveh, Y. Ben-Shmuel, and M. Tur, *Flutter Tests of the Pazy Wing*, AIAA Journal 60(9), 5414-5421, DOI: https://doi.org/10.2514/1.J061717. Supplied file: `pazy.pdf`, eight scanned PDF pages, all visually inspected. Relevant locations: test sequence and gauge discussion, p. 5416; Figs. 4-6, p. 5417; Tables 1-2, p. 5418; harmonic content and conclusions, pp. 5419-5420.

The immediate target is the instantaneous compressive strain-versus-speed history in the experimental paper's Fig. 6, reproduced as Fig. 3(a) of the first paper. Fig. 3(b,c) of the first paper are numerical predictions using different aerodynamic formulations. They are useful model comparisons, not the experimental curves themselves.

## What the first procedure actually prescribes

At fixed root incidence 5 degrees, Section V prescribes these numerical ramps:

| Branch | Airspeed law, with local branch time in seconds | Duration | Rate |
| --- | --- | --- | --- |
| Upsweep | `U(t) = 35 + 0.425*t` | 40 s | +0.425 m/s^2 |
| Downsweep | `U(t) = 52 - 0.85*t` | 20 s | -0.85 m/s^2 |

The plotted ordinate retains the static bending contribution and the oscillation. Connect instantaneous `(U(t), strain(t))` samples in time order, with blue for upsweep and red for downsweep. Do not detrend, take absolute values, or replace the history by its RMS for this comparison. A separate detrended signal is useful for identifying growth, frequencies, and amplitude. The dashed steady reference may be omitted as requested; its omission does not change the required dynamic calculation.

The first paper also runs small airspeed steps, holds for several seconds, and retains extrema from late cycles after stabilization (pp. 9-10). That is the procedure for Figs. 6-7, which investigate bifurcations with reduced rate effects. It is a different numerical experiment from the continuously changing speed curves of Fig. 3. The paper does not give a single universal dwell duration or an automatic amplitude-stationarity criterion.

The first model is a geometrically exact strain-based beam in UM/NAST coupled to Peters finite-state inflow/strip aerodynamics, with a calibrated spanwise correction; the alternative adds a dynamic stall model. Section IV explicitly includes gravity for the upright wing. Therefore, matching the ramps in an AeroBeams/UVLM model applies the sweep protocol but does not reproduce the full UM/NAST formulation. Adding the paper's strip-theory tip correction to UVLM without a separate derivation would risk duplicating three-dimensional effects.

Most importantly, the supplied first paper does not completely specify the numerical ramp initial states, perturbation amplitudes, preparation holds, or whether every dynamic state was carried between its plotted branches. Those details cannot be inferred uniquely from the appearance of Fig. 3. A continuous 5-degree reversal is a useful controlled hypothesis, not a documented exact reconstruction of the authors' initialization.

## What the experimental curve represents

The experimental paper's Section IV and Figs. 3-4 establish the loading history:

- Test 1: increase speed at approximately 7 degrees, change root incidence from 7 to 5 degrees at the high-speed endpoint, then reduce speed at approximately 5 degrees.
- Test 2: increase speed at approximately 5 degrees, change incidence to approximately 3 degrees at high speed, then reduce speed.

Consequently, the blue 5-degree branch is from test 2 and the red 5-degree branch is from test 1. The red branch has experienced a preceding 7-degree history and a pitch change. Neither a fresh static 5-degree state at 52 m/s nor a single fixed-5-degree round trip reproduces that complete experimental history. The experimental speed rates vary in time, and actual incidence differs slightly from the nominal settings; the first paper's linear ramps approximate their average rates.

The experimental Table 1 provides clear comparison values:

| 5-degree branch | Entry into LCO | Exit from LCO |
| --- | --- | --- |
| Increasing speed | 43 m/s | 51 m/s |
| Decreasing speed | 48 m/s | 36 m/s |

For comparison, the first paper's Table 2 gives `(onset up, offset up, onset down, offset down)` of `(44, 52, 49, 39)` m/s with nonlinear aerodynamics and `(46, 51, 46, 39)` m/s with linear aerodynamics. These are model-dependent outcomes; do not impose them as UVLM stability boundaries.

The second paper suggests aerodynamic stall is a major nonlinearity, based partly on the twist excursions. The first paper finds hysteresis with both aerodynamic formulations, only modest improvements from its stall model, and smaller hysteresis regions when rate effects are reduced. It also notes that a previous VLM model predicted weaker instability and different bifurcation behavior. Together, these findings make missing stall a possible model limitation, not a demonstrated explanation for the nearly flat current downsweep. Likewise, a satisfactory-looking upsweep alone does not validate the other branch or the numerical discretization.

## Strain reconstruction and gauge matching

The implementation currently uses the local beam-section relation

```text
strain_microstrain(s,t) = 1e6 * [gamma_1(s,t) + z_g*kappa_2(s,t) - y_g*kappa_3(s,t)]
```

This is the axial component of the centerline strain plus `curvature cross gauge_offset`, with the code's local coordinates: spanwise axis 1, +y forward, +z upper surface. Use the actual section states in those axes. Large global wing displacement is already handled by the nonlinear beam model; adding a separate displacement-derived bending correction would double count it. For the present bending convention, the upper surface is compressive and the result is negative. Microstrain is dimensionless strain multiplied by one million; curvature is in 1/m and offsets must be in metres.

`PazyLinearLCOSweep.jl:33` applies this reconstruction. `PazyWingUVLMCoupling.jl` saves the first element's midpoint axial strain and curvatures, and the continuous helper uses the same first-element midpoint. Defaults are `y_g=0`, `z_g=0.001125 m`. With the present nonuniform beam mesh, that midpoint is approximately 19.1 mm from the root; it is an approximation to root strain, not a measurement exactly at the clamp.

Gauge correspondence is unresolved. The experimental paper explicitly places the Fig. 5 sample sensor 127 mm from the root; Fig. 6 does not separately establish its sensor's complete spanwise and section coordinates in the supplied text. The first paper calls its numerical output upper-surface root strain. Do not silently equate all three signals, or automatically assign 127 mm to Fig. 6. Confirm the Fig. 6 FBG channel, spanwise location, leading/trailing-side placement, spar thickness, and zeroing convention from the experimental data/metadata. The official workshop page is https://nescacademy.nasa.gov/workshops/AePW3/public/wg/largedeflection; it currently directs users to the working-group lead for repository access.

Proposed change: retain the existing root proxy, and add independently labeled outputs for specified gauge stations `(s_m,y_m,z_m)`, using a documented section reconstruction/interpolation consistent with the beam discretization. Export axial and both bending contributions separately. At minimum compare upper/lower spar faces and the two experimental sensor sides once their locations are known. Neither rescaling an offset nor changing the sign to match the plot is a valid substitute for establishing the gauge geometry.

The main strain-versus-speed plot must retain the mean compression. The experimental paper multiplies the Fig. 5 time trace by -1 for presentation, whereas its Fig. 6 shows negative compression. Its 15 Hz high-pass operation applies to the oscillatory displacement/twist analysis associated with Fig. 9; it is not a prescription to high-pass the Fig. 6 strain curves.

## Current implementation and evidence

`run_pazy_LCO.jl` has been restored to two independent ramps, 35 -> 52 m/s in 40 s and 52 -> 35 m/s in 20 s. `PazyLinearLCOSweep.jl` computes a coupled static equilibrium separately for each branch, starts with zero structural velocity, and applies no controlled pulse. The UVLM wake and circulation are seeded from that branch's static initialization. No evolving upsweep state is transferred to the independent downsweep.

`run_pazy_continuous_LCO.jl` instead starts once at 35 m/s and carries motion and aerodynamic memory through a fixed-5-degree reversal at 52 m/s. Its explicit dt is 0.00041207376120325537 s, with 35 retained wake rows. The restored independent default uses about 0.00047547726030002614 s, with 30 rows. Comparing those defaults directly changes timestep and wake discretization as well as initialization. A controlled history comparison must specify the same dt and wake retention rule in both.

There is also an output-provenance issue at the time of inspection. `output/LCO_independent_sweeps/pazy_linear_LCO_settings.txt` now says 52 -> 35 m/s in 20 s and has a later modification time than the saved branch CSVs. The saved downsweep CSV itself still starts at 60 m/s and ends at 35 m/s at t=29.41176470588235 s. Those files must not be interpreted as one completed restored 52 m/s run. The previously generated `latest_sweep_analysis.md` describes the 60 m/s histories, which is consistent with their contents. The continuous log had not reached its reversal when inspected, so it supplies no completed continuous downsweep result yet.

The completed 60 m/s downsweep history shows almost no oscillatory band. Its detrended half ranges near 43-51 m/s are below 0.001 microstrain, while the accompanying upsweep reaches approximately 1015 microstrain in a 48.5 m/s bin. These statistics come from the earlier stored analysis, which removes a quadratic trend over 0.5 m/s bins; they are not directly the raw plotted range or a frequency-filtered modal amplitude. Starting at 60 m/s adds 9.4118 s before reaching 52 m/s. Without a controlled disturbance, that additional history can strongly reduce the seed for subsequent growth. It does not prove the true stability boundary.

The completed 45 m/s, 5 s diagnostics are stronger evidence that several effects must be separated:

| Case | dt (ms) | Initial strain (microstrain) | Last complete 0.25 s window, detrended half range (microstrain) | Early log-RMS growth (1/s) |
| --- | --- | --- | --- | --- |
| Pulse baseline | 0.475477 | -1540.77 | 1356.51 | 2.049 |
| Unforced | 0.475477 | -1540.77 | 109.27 | 2.262 |
| Half dt, same nominal wake age | 0.237739 | -1383.48 | 1.06 | 0.297 |
| Double retained wake, baseline dt | 0.475477 | -1559.99 | 1308.95 | 2.228 |

Source: `output/LCO_constant_speed_diagnostic/analysis_metrics.csv` and its report. The final complete window is approximately 4.64-4.89 s. The pulse is the existing smooth tip-force waveform, with coefficient 0.1 N, duration 0.04 s, starting at 0.1 s; its actual extrema are smaller than its coefficient. Growth values are finite-window fits, not computed eigenvalues. The unforced and refined responses are still growing. In particular, 1.06 microstrain is not an established asymptotic amplitude or evidence of stability.

Halving dt changes the initial equilibrium strain by about 10.2% in magnitude before the dynamic pulse. This prevents attributing all the dynamic differences to a single change in instability growth. Iteration residuals passing tolerance establishes convergence of the chosen discrete calculation, not convergence with respect to dt, mesh, or wake treatment.

One concrete numerical mechanism to investigate is near-wake regularization. The nominal shedding distance is approximately `eta*U*dt`. At 45 m/s, eta=0.1 gives about 2.14 mm for the baseline and 1.07 mm for half dt, compared with a fixed 1 mm vortex core. This is a plausible sensitivity, not an identified bug. The backend also contains an unsteady circulation-derivative force term; the small downsweep by itself does not establish that added-mass/unsteady loading is absent.

## Proposed code modifications

Keep `run_pazy_LCO.jl` as the original independent baseline and the existing continuous runner separate. Add new options or experimental runners so the comparison remains reproducible.

| Priority | Files / capability | Proposed modification and purpose |
| --- | --- | --- |
| 1 | Both sweep drivers and plot loader | Use a unique case directory/run ID, immutable per-run manifest, completion status per branch, and CSV metadata checks for endpoints, duration and dt. Stop the plot/analysis from mixing new settings with old histories. Preserve existing results. |
| 1 | `PazyLinearLCOSweep.jl`, continuous runner | Explicitly record dt, retained wake age and nominal length versus speed, core radius, shedding fraction, mesh, initialization provenance, and pulse settings. Compare reset/no-reset cases with identical numerical settings. Fixed row count represents changing physical wake length as speed changes; match the entire retention rule. |
| 1 | Sweep initialization / forcing | Add explicit controlled pulse coefficient, start time (or target speed), and duration. Keep the default unforced. Add consistent dynamic restart support separately from `initial_equilibrium`; a live state is not a static equilibrium. |
| 1 | Restart support | Preserve complete structural integration state and history, accepted beam velocities, aerodynamic geometry history, circulation and its history, wake, operating point, and clocks or consistent rebasing. Store a serialized, versioned restart if reruns from reversal are required. The existing turnaround CSV is diagnostic metadata, not a complete restart file. |
| 2 | Coupling and both response exporters | Add physical gauge-station outputs and component decomposition as described above. Preserve the original first-element proxy with an explicit station label. |
| 2 | Accepted-step diagnostics / analysis | Export local strain/velocity amplitudes, dominant frequency, harmonic content, and growth/decay windows alongside raw strain. Define a common noise-aware onset/offset criterion for all runs; retain raw data for comparison. |
| 2 | Constant-speed / continuation driver | Hold at selected speeds and carry the full dynamic state to adjacent speed steps. Compute late-cycle minima/maxima and mean only after stationarity checks. This implements the first paper's separate bifurcation procedure. |
| 3 | Static and dynamic model construction | Expose gravity consistently in both solvers and validate its orientation in the upright-wing axes. The current coupling explicitly sets gravity to zero. Check modal properties, skin/inertia configuration and strain reference together before comparing with the paper. Do not tune damping to obtain a desired red band. |
| 3 | Experimental-history runner | Support measured `U(t)` and root-pitch kinematics for the 7 -> 5 degree preparation. A scalar freestream-angle jump alone does not reproduce a physically pitching root and its surface velocities/inertial effects. Build and validate those kinematics consistently with the UVLM mesh and beam boundary conditions. |

Implementation follow-up: the independent driver now accepts optional
`downsweep_pulse_amplitude`, `downsweep_pulse_start`, and
`downsweep_pulse_duration` controls, with the default coefficient zero. It
exports the force as `tip_force_A1_N`. See `README_LCO.md` for the B2 call.
Stationarity and full dynamic restart controls above remain proposed additions.

## Proposed test cases and how to interpret them

First resolve the steady-load timestep sensitivity at representative speeds, then use the resulting settings for the branch investigation. The existing continuous run still provides useful exploratory information at its present discretization.

| Case | Procedure | Question answered |
| --- | --- | --- |
| N0: steady-load refinement | Fixed geometry, then coupled static solutions at 35, 45 and 52 m/s. Compare dt, dt/2, dt/4 with proportional wake-row counts to retain the same nominal wake age. Separately inspect longer wake, core radius and shedding fraction, then surface mesh and solver tolerances. | Does the reference equilibrium and the aerodynamic load become insensitive to discretization? |
| B0: matched independent baseline | 35 -> 52 in 40 s; fresh static 52 -> 35 in 20 s; no pulse. Use exactly the continuous run's dt/wake/core/mesh controls for this exploratory comparison. Save in a new directory. | What does state reset alone change relative to B1? |
| B1: continuous history | Existing single 35 -> 52 -> 35 simulation, with no reset and no hold. | Does the actual motion surviving the upsweep and reversal materially alter the downsweep? Carrying state does not guarantee that a large LCO still exists at 52 m/s. |
| B2: controlled downsweep seeds | Independent 52 -> 35 in 20 s, with otherwise identical settings. Apply the same short pulse when reaching 49 m/s (local t=3.5294 s), optionally compare a pulse at 45 m/s (t=8.2353 s). Bracket small pulse coefficients, e.g. 0.01 and 0.1 N, and include the unforced control. | Is the quiet red branch caused by a tiny seed and finite available growth, or does the response remain decaying? Pulse timing must be matched by speed when changing sweep rates. |
| B3: constant-speed basin/growth tests | At 45 m/s initially, then around the observed upper and lower transitions, compare a quiet equilibrium plus a small pulse with a physically developed oscillatory state from a compatible saved trajectory. Extend the duration until amplitudes/growth are resolved. | Do tiny disturbances grow, and do finite-amplitude oscillations persist at the same speed? This separates instability from a possible coexistence of attractors. |
| B4: descending settled continuation | Start within a confirmed oscillatory range, decrease speed in proposed 0.5 m/s steps, refining to 0.25 m/s near loss of oscillation. Carry the full state; also continue upward from the same branch when useful. Retain late-cycle extrema. | Is there a finite-amplitude branch below the small-disturbance boundary, as in the first paper's bifurcation analysis? This is not the 20 s ramp comparison. |
| B5: rate sensitivity | Repeat full 35-52 m/s ramps with durations up/down 40/20, 80/40, and, if necessary, 160/80 s. Keep dt, wake retention, and seed protocol fixed. | How much apparent hysteresis comes from delayed growth/decay rather than settled coexistence? |
| B6: experimental-history reconstruction | Use measured test 1 speed and pitch history to prepare the 5-degree downsweep; use test 2 for the 5-degree upsweep. If data are unavailable, label a 7-degree preparation/pitch transition as an approximation, with explicit assumed timing. | Does the preceding pitch/speed history account for differences that a fixed-5-degree reversal cannot isolate? |

For B3-B4, a proposed starting stationarity check is agreement of detrended RMS and mean over at least three consecutive windows, each spanning at least 25 dominant cycles, with an RMS change below about 2%. Check extrema and harmonic content as well. These are proposed screening thresholds; near a transition or in a modulated/chaotic response, longer records and different diagnostics may be needed. Do not label a run stable solely because its amplitude is small at an arbitrary stopping time. Record cases that have not settled.

In the small-disturbance regime, a useful interpretation is

```text
A(t) approximately A0 * exp(integral sigma(U(t)) dt)
log[A(U)/A(U0)] approximately integral from U0 to U of sigma(v)/Udot dv
```

Here sigma is a local growth/decay rate, not yet known over the current speed range. A tiny seed can remain visually undetectable through part of an unstable region, and an existing large oscillation can persist transiently into a stable region. This explains why a continuous sweep, a seeded sweep and a settled continuation answer different questions.

Possible outcomes and their implications:

- B1 or B2 restores an oscillatory red band: initialization/excitation and finite growth time contribute. Confirm that the band survives numerical refinement before claiming physical agreement.
- A pulse grows at a held speed but not visibly during B2: available time/rate is central; use B5 and compare the accumulated growth.
- Small disturbances decay but a developed finite-amplitude state persists at the same fixed speed: coexistence is consistent with subcritical behavior. Repeat from multiple amplitudes and after refinement to distinguish it from a long transient.
- Quiet and developed states both decay at converged numerical settings: the model may have a different unstable region or bifurcation structure. The first paper itself documents differences between VLM and strip-model predictions.
- The branch changes substantially with dt, core or near-wake treatment: resolve the numerical sensitivity before diagnosing physical hysteresis.
- Mean strain differs while fluctuation frequency and normalized motion are plausible: investigate gauge station, section offset, gravity/reference convention and steady-load bias. Gauge changes alone do not create the missing structural oscillation.
- Persistent differences remain after those checks: reassess aerodynamic model fidelity and experimental history; evaluate nonlinear aerodynamics with appropriate Pazy-flow calibration rather than assuming it will fix the curve.

The recommended first dynamic comparison is B0 versus the already-running B1, using matched discretization, followed by B2/B3 if the downsweep remains quiet. N0 is required before treating any branch amplitudes or flutter boundaries as quantitatively reliable. The original independent runner and its 20 s downsweep should remain available throughout.

## Concrete configuration available with the existing API

The following prepares B0 from the repository root, without changing the original runner's defaults. This is an executable configuration example, not a simulation performed during this review. It starts a full independent analysis when its final call is executed, so it should be scheduled deliberately alongside the current long-running job.

```julia
ENV["PAZY_LCO_RUN"] = "false"
include("lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/run_pazy_LCO.jl")
delete!(ENV, "PAZY_LCO_RUN")

matched_independent = run_pazy_linear_LCO(;
    initial_speed=35.0, peak_speed=52.0, downsweep_initial_speed=52.0,
    upsweep_duration=40.0, downsweep_duration=20.0,
    angle_of_attack=deg2rad(5.0),
    fixed_time_step=0.00041207376120325537,
    wake_length_chords=5.0, chordwise_panels=10, spanwise_panels=20,
    output_directory=joinpath("lib", "WingPropellerUVLM", "examples",
        "pazy_aerobeams_coupled", "output", "LCO_paper_comparison",
        "independent_52_matched"))
```

The default gauge, core, shedding fraction, and coupling settings in this call match the existing continuous test. Matching them controls the comparison; it does not establish that they are converged. A repeat should use a new case directory to retain earlier outputs. B2 can use the new optional pulse controls. B3-B6 require the extensions or preparation procedures listed above and should not be represented as already implemented sweep modes.
