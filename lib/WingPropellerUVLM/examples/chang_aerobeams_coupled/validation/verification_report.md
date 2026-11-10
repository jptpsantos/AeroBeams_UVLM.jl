# Chang nonlinear AeroBeams coupling: verification record

Date: 6 November 2026. Julia 1.11.5 on Windows.

All implementation changes are confined to the new
`examples/chang_aerobeams_coupled` directory. Existing analysis sources, package
projects and manifests were left unchanged. Generated results are stored only
in this example's `output/`.

## Executed checks

The assembled focused suite completed with **119 passed / 119 total**:

```powershell
julia --compiled-modules=existing --project=. lib/WingPropellerUVLM/examples/chang_aerobeams_coupled/validation/runtests.jl
```

See [verification.log](verification.log). It includes third-party ImplicitAD
method-declaration warnings; the Julia Test summary has no failed or errored
tests. PowerShell's stderr redirection formats these warnings as native error
records, so use the Test summary when reading this saved combined stream.

Coverage includes:

- Chang's complete installed mass and four gyroscopic blocks, their residual
  derivatives, spin reversal and analytical fixed-base whirl frequencies.
- Finite carrier rotation, moving-center-of-mass transport, Galilean invariance
  and mechanical power.
- Physical pylon-tip rotor inertia, angular-momentum derivatives and objectivity.
- Chang rigid-grid parity, finite-rotation virtual work, force/moment balance,
  shifted vortex force locations, different span meshes and aerodynamic rollback.
- Independent finite differences of the complete augmented beam/rotor residual.
- Converged nonlinear interacting steps with both reduced supports and beam
  pylons, plus rejection/restore of a deliberately failed step.
- Mean-loaded nonlinear initialization, developed-wake response and output.
- Growth-rate fitting against a signal with known exponential decay.

## Dry-wing comparison

`compare_wing_modes.jl` completed **9 passed / 9 total**, comparing 8, 16 and
32 elements. The model conserves the sum of the local spatial inertia blocks
and total wing mass, 177.3285166847 kg. The distributed AeroBeams reconstruction
and the original nodal-lumped Chang model are not identical discretizations.

At 32 elements:

| Mode, ordered by frequency | Existing Chang [Hz] | AeroBeams [Hz] | Difference |
|---|---:|---:|---:|
| 1 | 6.63652 | 6.59371 | −0.6450% |
| 2 | 20.81284 | 20.80096 | −0.0571% |
| 3 | 36.51431 | 36.42610 | −0.2416% |
| 4 | 42.90502 | 42.80883 | −0.2242% |
| 5 | 93.61634 | 93.96137 | +0.3686% |
| 6 | 112.71452 | 113.08562 | +0.3292% |

The maximum difference among these six frequencies is 0.645%. Higher-mode
differences are considerably larger on the coarse meshes. This is a dry-wing
structural comparison, not a coupled free-wake flutter calculation.

Full data: [wing_modal_comparison.csv](../output/wing_modal_comparison.csv).
Execution details: [modal_comparison.log](modal_comparison.log).

## Loaded interacting case

The deliberately coarse integration check uses the actual Chang distributions
with four wing beam elements, a 4×2 wing aerodynamic grid, four rigid blades
with 2×2 panels each, 65 m/s, 3° incidence, attachment at 50% span and fully
enabled wing/propeller aerodynamic interaction. Wake capacity is 12 rows per
surface. These settings exercise the algorithms; they are not a converged
production aerodynamic discretization.

Mean equilibrium converged in seven structural/aerodynamic outer iterations:

| Quantity | Result |
|---|---:|
| Initial mean wingtip z | 0.0301395864 m |
| Initial wingtip twist | −0.0434118° |
| Relative support pitch | +0.639145° |
| Relative support yaw | −0.358952° |
| Final mean iteration state change | 3.2153×10⁻⁷ |
| Final mean iteration load change | 3.9609×10⁻⁶ |
| Periodic aerodynamic initialization residual | 2.5126×10⁻⁵ |

The subsequent 14-step, 0.00946483 s developed-wake response retained all mean
loads. Its maximum residuals were:

| Residual | Observed | Configured limit |
|---|---:|---:|
| Geometry | 4.8726×10⁻⁹ | 1×10⁻⁷ |
| Scaled attachment velocity | 1.1631×10⁻⁷ | 1×10⁻⁶ |
| Aerodynamic loads | 9.1544×10⁻⁶ | 1×10⁻⁵ |
| Coupled mechanical equilibrium | 3.0339×10⁻⁸ | 1×10⁻⁶ |

Artifacts:

- [Time histories](../output/verification_mean_equilibrium_history.csv)
- [Response plot](../output/verification_mean_equilibrium_history.png)
- [Initialization iterations](../output/verification_mean_equilibrium_initialization.csv)
- [Resolved settings and residuals](../output/verification_mean_equilibrium_settings.txt)

## Interpretation limits

Both requested structural variations are implemented and exercise converged
coupled steps. The beam-pylon defaults match a static moment/rotation
compliance; distributed modal equivalence and measured rotor/nacelle properties
still require calibration. The two-DOF variation remains a reduced small-
relative-deformation support on a finite-rotation wing.

These results verify implementation and initialization. No converged flutter
speed, nonlinear branch continuation or Floquet stability boundary is claimed.
The longer paired-trajectory speed-sweep entry point is supplied, but its full
campaign was not run as part of these focused checks.

A stricter follow-up used eight wing elements / an 8×3 wing aerodynamic grid,
3×2 panels per blade, 36 wing wake rows and 72 propeller wake rows (one full
revolution at a 5° azimuth step). Its fixed-shape aerodynamic state did not
converge to the configured 0.1% periodic tolerance within twelve revolutions;
the reported residual was 0.0088319. The trim was correctly rejected. See
[longer_wake_trim_rejection.log](longer_wake_trim_rejection.log).

This rejection must not be interpreted as a validated trim or flutter result.
The longer response-check script uses an explicitly unloaded startup; its
results are distinct from the converged coarse mean-equilibrium case above.

## Longer perturbed startup responses

Both support formulations then completed **369 accepted timesteps** to
0.24946598 s with the longer wake, starting explicitly unloaded. A 2 N m smooth
pitch pulse acts from 0.05 to 0.07 s. All **14 additional checks passed**, and
the Julia process returned exit code 0. No periodic trim was assumed for these
runs. See [perturbed_responses.log](perturbed_responses.log).

| Quantity | Reduced support | Nonlinear beam pylon |
|---|---:|---:|
| Maximum wingtip z | 0.06421274 m | 0.06284696 m |
| Maximum absolute pitch | 1.34034° | 1.42152° |
| Maximum absolute yaw | 1.03021° | 1.19045° |
| Maximum FSI iterations | 13 | 14 |
| Maximum load residual | 9.9964×10⁻⁶ | 8.2081×10⁻⁶ |
| Maximum coupled mechanical residual | 1.1793×10⁻⁷ | 9.9876×10⁻⁷ |

Pitch/yaw are relative support coordinates for the reduced model and absolute
shaft-orientation angles in U for the beam-pylon model. Their numerical values
are not an equivalence test between the two support formulations. Pylon modal
calibration remains separate from these integration checks.

Artifacts:

- [Reduced-support response plot](../output/verification_reduced_history.png)
- [Reduced-support history](../output/verification_reduced_history.csv)
- [Beam-pylon response plot](../output/verification_beam_pylon_history.png)
- [Beam-pylon history](../output/verification_beam_pylon_history.csv)
- [Numerical response summary](../output/perturbed_response_summary.csv)
