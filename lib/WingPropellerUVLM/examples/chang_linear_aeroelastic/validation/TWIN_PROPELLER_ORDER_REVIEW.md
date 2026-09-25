# Investigation of the apparent outboard-first flutter result

Investigated on 2026-09-24, using the current `chang_case.jl` and the saved
`output/chang_linear_imperial_uvlm_history.csv` (75 m/s, two propellers, 5 s).
Production code and the user's case settings were not changed. The saved summary
does not record every resolved setting, so current defaults cannot establish the
complete provenance of that historical run.

**Conclusion:** no swapped propeller indices or load-transfer sign error was found
in the checks below. The saved response does not establish reversed flutter onset:
both propellers decay. There are concrete analysis limitations and numerical
settings that should be resolved before interpreting an outboard-first result as
a physical instability. No converged flutter sweep was performed in this review.

## 1. The apparent damping ordering depends on the fitting window

Applying the speed-sweep helper's mean-subtracted absolute-peak/log-amplitude
regression separately to all four saved channels gives:

| Channel | Growth rate, 1–5 s (1/s) | Growth rate, 3–5 s (1/s) |
| --- | ---: | ---: |
| Inboard pitch | -0.5842 | -0.4375 |
| Inboard yaw | -0.5954 | -0.4623 |
| Outboard pitch | -0.4216 | -0.4915 |
| Outboard yaw | -0.4029 | -0.4736 |

All rates are negative. Over 1–5 s the outboard response decays more slowly;
over 3–5 s the inboard response decays more slowly. These are channel-envelope
fits, **not separately identified aeroelastic eigenvalues**. Pitch RMS decreases
from 2.156 to 0.359 degrees inboard and from 1.955 to 0.574 degrees outboard
between the 1–2 s and 4–5 s windows.

Kher's twin-propeller study explicitly discusses modal beating and the failure
of simple logarithmic decrement for non-monotonic peaks. Its configurations are
windmilling-trimmed. Its A–C positions are 83% and 42% semispan.
[Kher et al., author-uploaded paper](https://www.researchgate.net/publication/404059967_Impact_of_Spanwise_Propeller_Placement_and_Aerody-namic_Interaction_Effects_on_Whirl_Flutter_of_Wing-Twin-Propeller_Configurations).

The current equal, simultaneous pitch pulses excite both propellers. A five-second
channel envelope can contain several modes with different participation and phase.
Track both whirl modes across speed and check identification stability against
window length, model order, and excitation choice before assigning an onset order.

There is also a confirmed **single-propeller assumption in the analysis tools**:
`studies/run_chang_speed_sweep.jl:115–116` reads only propeller 1 pitch/yaw;
`validation/analyze_chang_damping.jl:17` reads only propeller 1 pitch. These helpers
cannot determine which of two propellers governs flutter. This does not alter
the simulated histories, but it limits any flutter conclusion obtained from them.

## 2. Wake truncation is the leading aerodynamic numerical suspect

Current resolved settings give:

| Quantity | Value |
| --- | ---: |
| Time step | 0.00029672236 s |
| Wing wake rows | 25 |
| Propeller wake rows | 72 (half a revolution) |
| Wing nominal wake age / convection length | 0.007418 s / 0.55635 m |
| Propeller nominal wake age / convection length | 0.021364 s / 1.60230 m |
| Physical hub-to-elastic-axis distance | 1.70688 m |
| Physical hub-to-wing-trailing-edge distance | 2.96688 m |

Lengths above are `U * rows * dt`, **not measured free-wake extents**. Induction,
blade chord, and wake deformation modify the actual positions. Nevertheless, the
propeller's nominal retained convection length is shorter than the hub-to-elastic-
axis distance and much shorter than the distance to the trailing edge. The wing
wake is only about 0.31 chord long.

`src/wing_propeller/Initialization.jl:154` sets `system.trailing_vortices .= false`.
Thus there is no semi-infinite continuation compensating for the short finite
wake. Truncation can alter the wing's unsteady loads and the aerodynamic feedback
on each rotor differently. Its effect on the flutter ordering has **not** been
measured here.

The automatic wing rule in `src/chang_configuration.jl:111` is also unsuitable
for an isolated time-step convergence study: `rows = 5 * chordwise_panels`,
independent of `dt`. Halving the azimuth step halves wing wake age and length
unless the row count is doubled. Refining chordwise panels changes wake length
too. Such studies otherwise mix time/spatial discretization with wake truncation.

First compare progressively longer physical wake ages (including 1, 2, and 3
propeller revolutions), extending the wing wake independently. Hold retained
physical wake age fixed when changing the time step. Convergence, rather than one
selected longer wake, must establish adequacy.

## 3. The current coupling checks do not measure aerodynamic lag error

`src/chang_simulation.jl:424` evaluates the next-step aerodynamic load using the
previous structural displacement, then advances the structure once. The wake
is committed from that lagged geometry. This is the configured `:loose_explicit`
scheme, not a converged implicit aeroelastic step.

`src/aeroelastic/StructuralTimeIntegration.jl:278–321` checks structural equilibrium
against the supplied explicit load and returns zero state/load residuals. Those
zeros and `all_coupling_steps_converged = true` in the saved summary therefore
do not demonstrate time-step independence or agreement between the load and the
updated geometry. Lag can bias damping near a small separation between critical
speeds. Its sign and size for the two modes remain unquantified.

Compare `:implicit_predictor_corrector` with the explicit result at the same
physical setup. Also halve `dt` while maintaining both wake ages. Keep the
structural integrator fixed during the coupling comparison; subsequently check
sensitivity to Newmark's configured numerical damping (`newmark_alpha = 0.05`).

## 4. The operating point and reference configuration need alignment

`chang_case.jl:31` prescribes 1217 rpm at 65 m/s, then
`src/chang_model_parameters.jl:331` scales RPM with speed at fixed advance ratio.
The main run's “trim” averages and subtracts loads at zero structural displacement
(`src/chang_simulation.jl:507`). It does not adjust RPM to zero mean shaft torque.
There is a separate windmilling trim script, but the main driver does not call it.

Santos's multi-propeller paper reports a windmilling trim of 1275 rpm at 65 m/s
and uses moving-block identification. Its illustrated UVLM B–F case places the
propellers at 83% and 18%, with excitation applied to the outboard propeller.
Its position labels differ from Kher's: the current 42%/83% pair is Kher A–C but
Santos B–D, not Santos B–F. These differences matter when selecting comparison
figures. [Santos et al., author-uploaded paper](https://www.researchgate.net/publication/405402273_Aeroelastic_Stability_Analysis_of_Wing-Multi-Propeller_Systems_Under_Aerodynamic_Interactions_Using_the_Unsteady_Vortex-Lattice_Method).

The prescribed 1217 rpm is about 4.55% below 1275 rpm. Copying a reference RPM
alone would not establish windmilling for the present aerodynamic discretization;
verify the mean torque for the actual model and operating condition. Changing
RPM alters both the aerodynamic and gyroscopic terms.

## 5. The actual installation differs from the displayed positions

`src/chang_model_parameters.jl:403` rounds attachments to beam nodes. On the
current 20-element mesh:

| Propeller | Requested eta / y | Actual eta / y |
| --- | --- | --- |
| Inboard P1 | 0.42 / 3.150 m | 0.40 / 3.000 m |
| Outboard P2 | 0.83 / 6.225 m | 0.85 / 6.375 m |

`src/wing_propeller/Initialization.jl:100` uses these same rounded nodes for the
aerodynamic geometry. This is **not an aero/structural attachment mismatch**;
the complete installation moves, while the exported labels retain the requested
values. The spacing increases from 3.075 to 3.375 m. Quantitative reference and
mesh comparisons must use the actual positions or a mesh containing the desired
attachment nodes. This discrepancy alone does not prove the reversed trend.

## Checks that passed

The accompanying `investigate_twin_propeller_order.jl` audit found:

- Reversing the two attachment entries and permuting the associated coordinates
  gives identical M, C, and K matrices (zero norm differences).
- M and K are symmetric; the undamped gyroscopic C is skew-symmetric. The
  minimum mass eigenvalue is positive (0.0784085).
- All 124 free-coordinate aerodynamic load transfers agree with independent
  central differences of the actual vortex-vertex geometry. Maximum scaled
  error is 1.83e-9 at zero state and 1.73e-9 at a nonuniform perturbed state.
  Frozen synthetic forces are used; this validates kinematics and load mapping,
  not the aerodynamic force law or wake convergence.
- The first structural/gyroscopic frequencies are 4.7551 and 4.8795 Hz, with
  outboard/inboard propeller-angle norm ratios 2.882 and 0.342 respectively.
  Their growth rates are numerically zero. These are **not aeroelastic poles**.
- Plotting and organized-history export sort the propellers by attachment eta;
  the inboard/outboard labels follow the correct coordinate ordering.

The current finite-core factors are 0.05, and the working-tree change raised
them from 0.01. Core sensitivity should be checked with mesh/wake/trim fixed.
The output summary does not record the factors, so this review does not attribute
the saved history to that particular edit.

## Reproduction and remaining evidence needed

```powershell
julia --startup-file=no --compiled-modules=existing --project=lib/WingPropellerUVLM lib/WingPropellerUVLM/examples/chang_linear_aeroelastic/validation/investigate_twin_propeller_order.jl
```

The audit builds current defaults with environment overrides disabled and reads
the saved history; it does not rerun the historical simulation. An optional first
argument selects another standard Chang history CSV. The existing history needs
to cover the 1–5 s windows used by this diagnostic.

To identify the cause rather than a plausible contributor, obtain both modal
growth-rate curves across speed after matching physical positions/trim and
converging wake extent. Repeat with implicit coupling and a smaller time step.
Compare inboard-only and outboard-only excitation to separate mode participation
from stability. None of the successful indexing/virtual-work checks establishes
which aeroelastic mode crosses zero first.
