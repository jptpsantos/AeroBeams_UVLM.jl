# Pazy coupled-run failure investigation

## Scope and status

Read-only audit of the production implementation, with separate diagnostic scripts and output directories. No production solver, runner, visualization, or existing result was changed by this investigation. This is not a certification of the aeroelastic model or its predicted LCOs.

The attached failure was reproduced twice at **t = 2.8893635 s**, including the initial Newton residual **67490.9572950132**. One repeat used single-threaded BLAS. The accepted histories end at t = 2.888869 s. Changing BLAS thread count did not remove the failure.

Conditions: U = 50 m/s, incidence = 3 degrees, chord = 0.0989 m, dt = 0.0004945 s = 0.25 c/U, 10 x 20 aerodynamic panels, 15 beam elements, 40 retained wake rows (10 nominal chords), core radius 0.001 m, loose coupling, original Newton controls. The reproduction runs target 3 seconds instead of the runner's current 6 seconds. Animation/history storage was reduced in the diagnostic only; the numerical coupling loop was retained.

## 1. Confirmed error: incomplete dynamic rotational Jacobian

`src/Core.jl`, `element_states_rates!`, defines

    pdot = (2/dt)*p - pdotEquiv

The dynamic force and moment residuals include

    length/2 * (Rdot*R0*P + R*R0*Pdot)
    length/2 * (Rdot*R0*H + R*R0*Hdot)

The rotational derivatives in `element_jacobian!` (lines 1114 and 1128 at audit time) include the derivative of Rdot with pdot fixed, but omit the additional chain-rule term arising from pdot's dependence on p. For each rotation component j, the missing force/moment terms are

    length/dt * R_pj * R0 * P
    length/dt * R_pj * R0 * H

These must receive the same assembly, hinge transformations, and force scaling as the existing terms. The discrete-time contribution belongs in the dynamic tangent, not in the generic partial-derivative utility or steady-state tangent.

Independent central differences of the full assembled residual confirm this at the saved failing state:

| Rotational direction | Original relative derivative error | With diagnostic chain-rule correction |
| --- | ---: | ---: |
| 1 | 0.057526 | 6.54e-11 |
| 2 | 0.113960 | 1.03e-10 |
| 3 | 0.022538 | 5.09e-11 |

Table uses perturbation h = 1e-5. Checks at h = 1e-6 and 1e-7 give the same original discrepancy, with corrected errors below 1.2e-8. Maximum element rotation-parameter norm is 0.302, so this result is not caused by crossing the rotation-parameter scaling threshold. Other state-block directional derivatives were also checked in `step_replay.txt`.

**Important limit:** the correction makes the tested tangent consistent; it is not proof that it eliminates the long-time growth. Undamped Newton using the corrected tangent still diverges from the already highly disturbed saved state. Newton robustness and the preceding coupling history must also be addressed.

## 2. Loose coupling accepts an increasingly inconsistent interface

The current runner selects `coupling_scheme = :loose`. In `PazyWingUVLMCoupling.jl`, the aerodynamic evaluation uses extrapolated geometry, the structure advances once, and the mismatch is accepted without iteration. The aerodynamic trial committed at the end of the step retains the predicted geometry; the structural solution generally has a different geometry.

This is the approximation made by this explicit partitioned scheme, not by the spatial force-transfer map. `coupling_relaxation` and the FSI convergence tolerances are not applied in loose mode. Consequently FSI=1 and r_F=n/a are expected, but do not establish interface convergence.

| Time interval (s) | Maximum interface mismatch / chord | Maximum beam angular speed (rad/s) |
| --- | ---: | ---: |
| 0–0.5 | 0.000427 | 8.75 |
| 1–1.5 | 0.000842 | 11.19 |
| 2–2.5 | 0.013128 | 54.25 |
| 2.5–2.8 | 0.045952 | 148.70 |
| 2.8–failure | 0.260458 | 1472.84 |

The last mismatch is approximately 25.8 mm for a 98.9 mm chord. Maximum aerodynamic surface-motion speed reaches 70.6 m/s. Thus the Newton failure is preceded by severe interface and high-rate growth; it is not an isolated harmless first-iteration residual spike.

An increment-based interface-work diagnostic also detects a growing difference between work evaluated on the predicted and solved surfaces. This is **not** a complete energy balance and cannot by itself quantify artificial energy injection or prove its sign.

A full strong-coupling comparison is still running. It has passed 1.4 s with tightly matched interface geometry, but has **not yet reached 2.889 s**. It is therefore not evidence yet that strong coupling cures this late failure or produces a stable LCO.

## 3. Newton recovery and convergence checks need attention

At the exact saved failing step, with the same previous state and aerodynamic loads:

| Replay | Outcome | Reassembled final residual norm |
| --- | --- | ---: |
| Original Newton settings | Fails | 3.21e7 |
| Refresh Jacobian every iteration | Reports convergence | 9.88e-4 |
| Zero applied loads, refreshed Jacobian | Reports convergence | 3.12e-3 |

Removing the applied loads changes the initial residual only from 67490.96 to 67491.84. This explains why external-load-factor reduction is ineffective here: the large imbalance is dominated by the evolved dynamic state, not the magnitude of the current external load alone.

`src/SystemSolver.jl` accepts an iteration when the absolute residual **OR** the relative state increment is small. The printed absolute residual precedes the most recent update. A small relative increment does not necessarily mean the updated residual meets the requested absolute tolerance. At t = 2.8676055 s, an accepted step has an independently reassembled residual of 1.79e-4 despite an absolute tolerance of 1e-6.

Refreshing the existing incomplete Jacobian is helpful in this replay, but is not a full correction. A fresh-Jacobian direct solve still has residual 7.56e-6 after 20 updates when actual residual convergence is required. A diagnostic backtracking experiment with the corrected tangent reduces the residual from 67490.96 to 3938.28, then stalls. Thus neither corrected undamped Newton nor this simple backtracking implementation recovers the already disturbed state in this test. The experiment is recorded in `dynamic_tangent_audit.txt`; it changes only the nonlinear solution procedure, not the physical time step or loads. Corrections must be tested from before the interface-growth onset, not judged solely by trying to salvage the last step.

## 4. Checks that passed, and what they establish

- Existing Pazy coupling regression: **35/35** checks passed. Includes coordinate orientation, elastic-axis offsets, resultant force/moment transfer, finite-rotation virtual work on noncoincident meshes, predictor/relaxed offset consistency, and updating active load BCs after a Newton retry.
- UVLM migration/core checks: **49/49** passed.
- Main UVLM checks: **167/167** passed when run separately from an unrelated preceding Chang assertion failure.
- The unmodified full package test entry point is **not all green**: `test/general_structural_assembly.jl:165` compares damping matrices with exact equality and fails on approximately 309.5662930277377 versus 309.56629302773774. This was not modified or hidden.
- A separate dry-structure run uses the same Pazy beam, dt and original Newton controls. A short smooth force/moment pulse is followed by free vibration, without UVLM. It completes **6067 steps to 3.0001315 s**. Free-vibration energy ranges from 0.000630132483650051 to 0.0006301344229224831 J; final/initial energy = **1.0000002685970426**. This small-amplitude check does not cover large-amplitude nonlinear dynamics, but does not support a generic dry-integrator runaway explanation.
- Logged bound-panel and wake-panel areas remain finite and nonzero through the last accepted coupled step. This excludes panel-area collapse in the logged data, not every possible wake proximity, self-intersection, or force-model issue.

Spatial virtual-work conservation at matching configurations must not be confused with time-discrete conservation across a lagged fluid–structure interface.

The inspected structural update uses the trapezoidal relation `qdot_new = 2*(q_new-q_old)/dt - qdot_old` on the mixed first-order states. This wrapper does not expose a Newmark-beta/alpha parameter. Introducing a dissipative integration scheme would require consistent changes to the state-rate updates and their derivatives, not merely changing a beta value in the run file.

## Recommended next implementation and verification

1. Correct the structural dynamic rotational Jacobian and add finite-difference regression tests at several deformed/moving states. Test steady/eigen tangents separately so the transient correction is not applied there.
2. Use strong partitioned coupling for the reference LCO investigation; verify both geometry and loads, and compare against loose coupling only as a controlled numerical experiment.
3. Refresh the tangent during verification. Independently check the residual after convergence using suitable block scaling. Add a safeguarded Newton strategy if needed; do not advance through failed solves or clip the physical response.
4. If a step is rejected, restore **both** fluid and structural accepted states before retrying. Existing catch/rollback handles the UVLM trial and stops; it is not a complete coupled adaptive-step implementation.
5. Repeat beyond the observed failure time, then perform separate time-step, aerodynamic mesh, retained-wake-length, core-size and FSI-tolerance studies. Keep physical wake length controlled when changing dt. Compare cycle amplitude, mean deformation, frequency and energy/work trends, not merely whether the run finishes.
6. Only then assess an alternative dissipative integrator. An alpha parameter cannot correct the missing Jacobian term or guarantee a physically correct LCO. Adding numerical damping before isolating these errors risks changing the inferred flutter/LCO response.

The confirmed defects and numerical risks must be addressed before attributing the late response to physical LCO behavior. Conversely, a nonlinear structure plus a free wake does not by itself establish that this exact model and operating point must saturate into the experimental LCO.

## Reproducible artifacts

All paths below are relative to this example directory:

- `audit_pazy_failure.jl`: isolated full-run reproduction and per-step diagnostics; accepts `loose`, `strong`, or `loose_freshjac`, followed by duration in seconds.
- `replay_pazy_failure.jl`: repeat the exact failed structural step and check state-block directional derivatives.
- `audit_pazy_dynamic_tangent.jl`: compare original and discrete-chain-rule-corrected tangents, plus local Newton experiments. No production functions are replaced.
- `audit_pazy_dry_structure.jl`: uncoupled beam control experiment.
- `summarize_pazy_failure.jl`: compare the saved histories and generate a diagnostic figure. The strong history/plot is partial until that run finishes.
- `output/failure_audit_20260927_143526_loose/`: first exact reproduction, failure checkpoint, `step_replay.txt`, `dynamic_tangent_audit.txt`, `comparison_summary.txt`, `coupling_comparison.png`.
- `output/failure_audit_20260927_144139_loose/`: second exact reproduction with single-threaded BLAS.
- `output/failure_audit_20260927_143903_strong/`: ongoing comparison.
- `output/dry_audit_20260927_162719/`: dry-beam history and summary.

Example (run from repository root):

```powershell
julia --project=. lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/audit_pazy_dynamic_tangent.jl lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/output/failure_audit_20260927_143526_loose/failure.jls
```
