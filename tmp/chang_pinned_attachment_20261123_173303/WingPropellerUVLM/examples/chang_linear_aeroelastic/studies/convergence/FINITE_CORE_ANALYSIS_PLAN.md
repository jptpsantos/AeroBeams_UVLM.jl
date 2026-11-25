# Plan for finite-core, aerodynamic, and damping convergence

Prepared 2026-09-12. This document defines the proposed analysis; it does not
report completed production convergence runs. Numerical thresholds below are
initial engineering targets, not universal literature requirements.

The objective is to identify a documented range of core radii for which CT/CQ,
relevant wing/hub loads, and modal growth rates are insensitive to further
numerical refinement. Physical vortex-core calibration is a separate objective.
Retain the present smooth kernel during the initial study.

## 1. Freeze the physical problem and its scope

Record geometry, blade twist, density, speed, actual RPM, collective, symmetry,
interaction settings, load formulation, structural properties and excitation.
Do not retrim or change the physical problem between numerical cases.

The aerodynamic entry currently specifies:

| Setting | Current value |
|---|---|
| Speed | 65 m/s |
| Actual RPM | `1212 * 85 / 65`, approximately 1584.923 |
| Wing chord / blade chord | 1.8 m / 0.197 m |
| Propeller radius / blades | 1.15 m / 4 |
| Wing symmetry / aerodynamic interaction | true / false |
| Core law / nominal radius | fixed / 0.001 m |

Resolve whether the intended RPM is 1212 or 1584.923 before collecting final
results. Neither choice is inferred from speed automatically by this entry.
Wing image symmetry must correspond to the intended physical configuration;
all propeller blades remain explicit. Freeze the symmetry choice throughout.

Use the interaction-disabled problem as a baseline, then repeat the decisive
study with interaction enabled if wing-propeller aerodynamic interference is
part of the intended model. Keep their results in separate directories. CT/CQ
with interaction disabled cannot establish wing-load convergence.

Deliverable: a frozen operating-point record and a statement of whether the
claim concerns propeller coefficients, coupled wing loads, damping, or all three.

## 2. Complete the instrumentation before production runs

Keep the two current entry files. Add helpers inside `src/` rather than more
competing entry scripts. The current source fingerprint invalidates aerodynamic
evidence after numerical source changes, so finish required extensions first.

| Capability | Current status | Work needed for this plan |
|---|---|---|
| Absolute mesh counts; fixed/chord/segment core laws | Available | Use existing settings |
| Independent CT/CQ sweeps, phase errors and periodicity | Available | Retain current checks |
| Combined candidate and finer aerodynamic verification | Available | Retain current checks |
| Verified handoff to damping study | Available | Preserve verification |
| Crossed core x mesh x time studies | Not automated by current families | Add case generation and tables |
| Separate wing and propeller core controls | Not exposed by current shared callback | Extend geometry, configuration, metadata and handoff |
| Local wake/induction diagnostics | Not all summarized | Add inexpensive reductions and sparse snapshots |
| Wing-load and hub-moment convergence gates | Not part of CT/CQ acceptance | Add a separate assessment for coupled claims |
| Fit-window and excitation-amplitude studies | Not current refinement families | Add controlled repeated runs/postprocessing |
| Structural spanwise mesh convergence | Explicitly excluded from current damping sweeps | Separate structural study |

For independent source cores, retain a consistent core on a physical shared
edge and on its reflected copy. Verify cancellation of equal-circulation shared
edges. Propagate source-core identity into wakes, updates, trim, load evaluation,
snapshots and result metadata. Distinguish initial wake and bound cores only if
their roles are explicitly modelled. Do not silently modify the physical model.

## 3. Verify the kernel and solver interfaces

Run the package tests and the bounded independent-workflow validation. Review
the existing archived kernel audit before adding missing checks. Verify:

- Agreement with independent quadrature of the softened straight-line integral.
- Midpoint, endpoint and far-field behaviour at positive radius.
- Rotation/translation invariance and correct dimensional scaling.
- Reversal of circulation/orientation and segment-subdivision invariance at a
  common dimensional core.
- Shared-edge cancellation and matching finite/semi-infinite limits.
- Ring-centre induction against the matching regularized circular-line reference;
  off-axis/filament self-induction needs an appropriate additional reference.
- Core inheritance and preservation under wake motion and solver snapshots.

These checks establish implementation consistency. They do not calibrate a real
tip-vortex radius. Stop and resolve a failed kernel check before expensive sweeps.

## 4. Establish a settled, affordable baseline

Start with the present nominal mesh (wing 20 x 5; each blade 5 x 5), 1 mm fixed
core, 5-degree step, and two retained revolutions. Use `dry_run=true` to inspect
the matrix before execution. Measure warmed-up runtime per step and peak memory;
record compilation time separately. Estimate the expensive combined cases before
launching them. Do not extrapolate wall time solely from the number of cases.

Check duration separately: initially compare 8, 12 and 16 simulated revolutions
at the same retained wake age, averaging the last two revolutions of each run.
Extend if needed. At the selected duration, compare averaging the last two and
four revolutions provided both windows exclude transients. Duration changes are
separate configurations, not an existing sweep family.

Acceptance: CT/CQ means and phase waveforms cease drifting and the current
cycle-to-cycle RMS limits are met. Present defaults use absolute periodicity
limits of 1e-6 for each coefficient. Inspect their suitability and numerical
noise before production; freeze the criterion before accepting cases. Do not
increase the core or relax a criterion merely to make an unsteady wake pass.
Persistent aperiodic behaviour requires distinguishing numerical failure from
physical/statistical unsteadiness and revising the analysis method if necessary.

Use fresh rest states for independent cases. Reuse matching histories only with
the existing source/options/history validation. Test a refined/smaller-core case
for settling too; coarse-grid settling does not guarantee fine-grid settling.

## 5. Screen individual aerodynamic sensitivities

The current ranges are starting points, not known converged settings:

| Family | Levels |
|---|---|
| Wing span panels | 20, 40, 60 |
| Wing chord panels | 5, 10, 20 |
| Blade radial panels, per blade | 5, 10, 20 |
| Blade chord panels, per blade | 5, 10, 20 |
| Retained wake revolutions | 1, 2, 3 |
| Fixed core, metres | 0.01, 0.003, 0.001, 0.0003 |
| Azimuth increment, degrees | 5, 2.5, 1.25 |

This is 22 entries and 16 distinct parameter combinations. If all families pass,
the driver attempts two additional combined verification runs. These counts
exclude duration checks, repeated operating configurations and crossed studies.

Initially inspect mesh, then time step, then wake extent and core. Hold all
non-varied controls fixed within each family. Preserve physical retained wake
age as the step changes: `Nwake = ceil(360 * revolutions / azimuth_deg)`.
The total simulation duration must remain sufficient for every retained age.
For a propeller, `dt = azimuth_deg / (6 * RPM)`.

For coefficient q, use the existing assessment:

`E_q = max(abs(mean_q - mean_q_ref), RMS(phase_q - phase_q_ref))`

`T_q = max(absolute_tolerance_q, relative_tolerance * abs(mean_q_ref))`

Current starting tolerances are 1e-6 absolute and 1% relative. Near-zero torque
requires a justified absolute scale. Each family must be complete, its last two
levels must agree, and a selected candidate must agree with all subsequent
levels. Agreement against the finest case alone is insufficient.

If a family fails, extend its range or improve settling; do not accept the final
level by default. Re-centre nominal controls on promising settings and repeat
the necessary families, because coarse settings held fixed can mask sensitivity.
Preserve earlier campaigns rather than overwriting their interpretation.

## 6. Resolve core-mesh-time interactions

Use three core values bracketing the candidate region, two mesh configurations
(candidate and finer), and two time increments (candidate and half-step). This
produces 12 combinations at a fixed, adequately long wake extent. For example,
use 0.5, 1 and 2 mm only if screening identifies a plateau around 1 mm.
The decreasing core-family ordering in the current script would be 2, 1, 0.5 mm;
crossed cases require the planned extension.

Compare meshes at the same core/time step, time steps at the same mesh/core,
and cores at the same mesh/time step. Plot CT/CQ against radius with separate
curves for each discretization. Use the same wake age and settling criteria.

Record min/max and representative percentiles of core/chord and core/edge length,
non-self interaction distance/core, maximum induced velocity, wake segment
lengths, and wake geometry near the propeller, wing and truncation boundary.
Separate connected/self interactions from near encounters. Save sparse phase-
aligned wake snapshots; avoid dumping all pairwise interaction distances.

Pass only if the acceptable radius region persists under refinement. Otherwise
refine the affected mesh/time direction or investigate close encounters, load
evaluation and wake truncation. A core-sensitive result is an outcome to report,
not a reason to choose the radius that produces the preferred damping sign.

## 7. Separate component sensitivity and compare core laws

After separate source controls exist, keep the wing core fixed while sweeping
the blade core, then reverse the roles. Near the candidates, use three wing
radii x three propeller radii (nine combinations) if interaction effects require
it. Hold mesh/time/wake settings fixed and inspect wing/hub loads as well as CT/CQ.

Only after establishing this fixed-core baseline, compare chord and segment
laws. Match the initial dimensional core at an explicitly chosen reference
station, then report the actual radius distribution everywhere. For example,
1 mm corresponds to 0.005076 times the 0.197 m blade chord and 0.0005556 times
the 1.8 m wing chord. A common chord factor cannot reproduce both values.

For segment scaling, refine at fixed factor to examine that numerical method,
and separately adjust the factor to preserve a reference dimensional radius.
Distinguish those experiments. In chord mode the callback uses full local chord;
in segment mode it uses a spanwise/radial vortex edge. Neither is automatically
the wake's streamwise spacing.

Do not average answers from different laws or assume that equal numeric
parameters produce equal physical cores. Core-law comparison is model
sensitivity; within-law refinement is numerical convergence.

## 8. Verify the combined aerodynamic choice and save the handoff

Run the selected combined controls and a finer combined mesh/time/wake case at
the same fixed core. Use the current verified selection mechanism. Retain the
crossed-study evidence beside it; until the handoff is extended, the current
`verified` label certifies only the checks currently implemented, not every
additional assessment in this document.

For coupled wing-load claims, require a separate wing/hub-load assessment.
Finish any numerical-source changes before producing the final selection.
Save the operating condition, source/version hashes, dimensional core settings,
case matrix, raw histories, acceptance limits and combined verification evidence.

## 9. Establish that the damping estimate is identifiable

Load the verified selection in `chang_convergence_aeroelastic.jl`. First run a
single candidate response. Check trim completion before excitation and start
the fit after the impulse and its immediate transient. Match the trim/mean-load
state to the rigid aerodynamic calculation where applicable.

Current defaults are a six-second response, ten trim revolutions, a 0.15-second
1200 N m impulse, a one-second fit block and a 3-6 Hz frequency band. Inspect
the actual spectrum and structural modes before choosing the band: the configured
uncoupled pylon pitch/yaw frequencies are 7.97 Hz, which is not itself a prediction
of the coupled modal frequencies. A channel peak is not automatically one mode.
Track the same modes across cases and resolve beating/mode mixtures if present.

Compare the nominal impulse to half amplitude, normalized by impulse amplitude.
Use a smaller excitation if fitted rates depend on amplitude. Compare fit starts
and block durations that remain inside the free response; freeze an appropriate
physical window for step refinement. Extend the response (initially 6, 9, 12 s)
when weak damping or insufficient cycles prevents identification.

Current automatic checks require at least eight peaks and R-squared >=0.8.
Treat these as minimum screening criteria. Examine residuals and sensitivity to
fit windows; a high R-squared does not prove a single-mode response. Where useful,
compare the moving-block estimate with an independent modal/envelope estimate.

## 10. Conduct damping convergence

The current six families produce 18 entries (13 distinct controls): wing chord,
blade radial/chord counts, retained wake age, core and time step. They start at
the aerodynamic selection. `nothing` requests two additional refinement levels.
Choose explicit lists when the automatic levels are too costly or unsuitable.

Keep structural discretization fixed for these families. Independently check
coupling tolerances/iterations and structural integration so their errors do
not dominate the aerodynamic study. An adaptive stopping limit is not a damping
result; failed or truncated histories remain invalid.

Starting acceptance targets: changes in both tracked modal growth rates below
0.01 1/s and frequencies below 0.1 Hz, with valid fits and agreement between the
last two levels and the candidate. Tighten these targets if they cannot resolve
the stability margin. Reuse postprocessing for fit-window changes; rerun dynamics
only when a simulation parameter changes.

Repeat nearby core values with a refined mesh/time combination, including a
larger as well as a smaller core around the candidate. The present decreasing
core family alone does not perform this bracketing. After individual families
pass, run a combined finer aeroelastic verification; this is not currently an
automatic final check in the damping driver.

Report observed numerical/fit variation with the growth rate. If it is comparable
to the growth-rate magnitude, classify the sign as unresolved rather than
assigning stable/unstable. The code's lambda is a growth rate: positive grows,
negative decays. Do not confuse its sign with conventional positive damping ratio.

## 11. Optional extensions according to the intended claim

For complete aeroelastic mesh convergence, conduct a structural spanwise mesh
study: dry modal frequencies/shapes and mass/stiffness first, then coupled rates.
The current wing span panel count also defines structural elements and is held
fixed by the damping driver. A general full-model mesh claim needs this extra work.

For flutter-speed convergence, bracket a growth-rate zero crossing and refine
speed near it with candidate and finer discretizations and plausible core values.
Declare whether actual RPM, advance ratio or a trim schedule is held fixed across
speed. One operating point cannot establish flutter-boundary convergence.

For physical wake-core calibration, add an explicit initial core and age law
only after the numerical baseline. Validate against wake velocity/core measurements
or an appropriate higher-fidelity reference. Test time-step consistency of growth
and any stretching treatment. Keep this campaign separate from the constant-core
study; it changes the physical model.

## 12. Deliverables and run management

Use one output directory per frozen campaign. Retain raw CSV histories, study
TOML snapshots, logs, per-case runtime/memory, fit diagnostics and source hashes.
Keep failed cases visible. Stop escalation when new failures or unresolved
settling appear rather than spending the remaining budget on finer cases.

Final figures: CT/CQ versus actual control values; phase-waveform comparisons;
core-mesh-time plots; wing/hub-load comparisons where relevant; matched wake
snapshots; pitch/yaw response and fit residuals; growth rates and frequencies
versus numerical controls; and a table of accepted values and tested ranges.

The current plots use refinement level on the x-axis; the case matrix supplies
actual values. Add radius/log-axis and crossed-study plots for the final report.

Execution order: kernel checks -> settled baseline -> individual aerodynamic
sweeps -> crossed/core-component studies -> combined aerodynamic verification ->
identified damping response -> damping sweeps and combined verification -> optional
structural/flight-boundary/physical-core studies. Repeat earlier stages only when
changed settings, failed criteria or new evidence justify it.

References supporting the methodology, not prescribing this plan's tolerances:

- [VortexLattice.jl geometry/core callback](https://raw.githubusercontent.com/byuflowlab/VortexLattice.jl/master/src/geometry.jl).
- [OpenFAST/OLAF regularization and discretization guidelines](https://github.com/OpenFAST/openfast/blob/main/docs/source/user/aerodyn-olaf/RunningOLAF.rst).
- [QBlade vortex-core and wake formulation](https://docs.qblade.org/src/theory/aerodynamics/lifting_line/lifting_line.html).
- [Van Hoydonck and van Tooren: core-correction consistency](https://arxiv.org/abs/1204.2378).
- [Bhagwat and Leishman: generalized viscous vortex models](https://www.researchgate.net/publication/255470975_Generalized_viscous_vortex_model_for_application_to_free-vortex_wake_and_aeroacoustic_calculations).
