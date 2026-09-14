# Finite-core selection for Chang wing-propeller stability

The [fixed-radius convergence driver](CORE_CONVERGENCE.md) implements a
radius-by-mesh study with separate time-step and wake-age refinements,
CSV histories, plots, and an optional aeroelastic growth-rate comparison.

The core must be selected together with spatial resolution, time step, and
wake extent. A core that produces stable time histories is not thereby
validated. First separate discretization error from core sensitivity, then
check the aerodynamic damping and flutter boundary. The numerical trial
values below are a proposed study, not calibrated physical vortex radii.

## Current implementation and measured settings

The main case now exposes `aerodynamic.core_radius_m` (see `chang_case.jl`
for the active value). For example, `1e-3` selects a fixed 1 mm trial radius
for all surfaces and their shed wakes.
This overrides the core factors and stays fixed during deformation and wake
convection. Set it to `nothing` to select the original factor-based rule:

```julia
rc = max(segment_core_factor * ds, chord_core_factor * c)
```

Here `ds` is the three-dimensional spanwise/radial edge length and `c` is
the full local chord. A panel stores one radius for its ring edges. Wake
panels inherit the trailing-edge panel radius and retain it during convection;
there is no wake-age diffusion model. Same-surface bound induction is
unregularized, while cross-surface bound induction and wake induction use the
configured core in the Chang setup.

For the factor-based settings audited on 2026-09-10 before adding the fixed
option (wing 20x5, each blade 5x5, segment factor 0.1, chord factor zero),
geometry construction after the twist correction gives:

| Surface | Initial radius | Radius / full chord | Largest radius / panel chord |
|---|---:|---:|---:|
| Wing | 37.50 mm | 2.083% | 0.1042 |
| Blade | 23.0005-23.1429 mm | 11.675-11.748% | 0.5891 |

The blade range comes from the twisted three-dimensional panel edges.
Radial refinement changes these radii when the segment factor is fixed.
These settings therefore mix changes in discretization and regularization.

The implemented segment kernel integrates the softened Biot-Savart integrand
`dl x r / (|r|^2 + rc^2)^(3/2)`. Its infinite straight-filament limit is

```text
v_theta(r) = Gamma/(2*pi) * r/(r^2 + rc^2).
```

Thus changing the radius changes induced velocities, and potentially load
phase and aerodynamic damping. It is not merely a numerical damping knob.
For illustration, at a distance of 10 mm from an infinite filament, radii of
23 mm and 1.97 mm retain approximately 16% and 96% of the singular velocity.
This is an illustration of kernel sensitivity, not a measured encounter in
the present wake. Finite-segment and curved-wake behavior also depends on
geometry.

At 85 m/s, the configured advance-ratio scaling gives 1584.923 RPM and a
5-degree step of approximately 0.000525788 s. The retained propeller wake is
one revolution (72 rows); the automatically sized wing wake is 50 rows,
or approximately 0.694 rotor revolutions. Halving the azimuth step while
keeping these row counts halves the retained wake age.

## What the references establish

- [OpenFAST OLAF theory, regularization and diffusion](https://openfast.readthedocs.io/en/main/source/user/aerodyn-olaf/OLAFTheory.html#regularization):
  smoothing a discrete vortex field and modeling physical viscous diffusion
  are different operations. A constant smoothing radius does not establish
  a physical diffusion model.
- [OpenFAST OLAF input definitions](https://openfast.readthedocs.io/en/main/source/user/aerodyn-olaf/InputFiles.html):
  bound and wake regularization have separate parameters, with several
  kernels and choices of dimensional, chord-based, or spacing-based radii.
  This supports making those choices explicit in a solver interface.
- [OpenFAST OLAF operating guidance](https://openfast.readthedocs.io/en/main/source/user/aerodyn-olaf/RunningOLAF.html):
  its current lifting-line guidance suggests spacing-based factors near 0.6
  and addresses time step, wake extent, and downstream regularization.
  The wind-turbine formulation differs from this panel UVLM; its factors
  and wake-length recommendations are not calibrated Chang settings.
- [Van Hoydonck, Gerritsma, and van Tooren, 2012](https://arxiv.org/abs/1204.2699),
  *On Core and Curvature Corrections used in Straight-Line Vortex Filament
  Methods*: vortex-ring benchmarks expose errors from local curvature and
  inappropriate core corrections. Apparent convergence against the finest
  computed case can be misleading. The relevant check is against an
  independent reference for the same regularized model.
- [Bhagwat and Leishman, 2002](https://www.researchgate.net/publication/255470975_Generalized_viscous_vortex_model_for_application_to_free-vortex_wake_and_aeroacoustic_calculations),
  *Generalized Viscous Vortex Model for Application to Free-Vortex Wake and
  Aeroacoustic Calculations*: wake core growth affects geometry, induced
  velocities, and loads; eddy-viscosity modeling is informed by vortex
  measurements. Rotor-scale empirical constants require justification.

## Recommended numerical study

1. **Freeze the physical case and correct geometry first.** The Chang twist
   samples are now located at segment midpoints. Record airspeed, RPM,
   collective, density, interaction setting, and structural model. Existing
   trim values and convergence results must be checked for this geometry.

2. **Run fixed-operating-point comparisons first.** Hold RPM and collective
   fixed across numerical variants so that numerical effects are identifiable.
   Afterwards, if the physical comparison is windmilling or a specified thrust,
   re-trim each retained variant to that same condition. Record the new RPM
   or collective; do not mix these two comparisons.

3. **Separate mesh and core refinement.** Use `core_radius_m` to hold the
   dimensional radius fixed. Alternatively, select `core_radius_m = nothing`,
   set the segment factor to zero, and use a fixed full-chord factor on the
   constant-chord geometry. Refine chordwise and spanwise/radial meshes at fixed core, then repeat
   for smaller cores. This first establishes the discretization limit of a
   regularized model; it does not prove an inviscid limit or physical validity.
   A spacing-based strategy remains possible, but must be studied as a joint
   mesh/core refinement rather than as a pure mesh study.

4. **Use more than one core on more than one mesh.** The existing chord-factor
   family `0.04, 0.02, 0.01, 0.005` is a useful diagnostic range. Include a larger
   core if needed to connect it to the current approximately 0.117-chord blade
   core. The nominal 0.01 is only a trial value; a smaller core may demand a
   finer mesh and time step. If the trend persists under refinement, report
   the remaining core sensitivity instead of declaring the smallest value
   converged.

5. **Refine time at fixed wake age/extent.** For example, compare 5, 2.5, and
   1.25 degrees, scaling wake row counts inversely with the angle increment.
   For one retained revolution these require 72, 144, and 288 rows. Test wake
   extent separately, for example 1, 2, and 3 revolutions as an initial study.
   Extend the range if the final pair disagrees. Neither this sequence nor
   OLAF's wind-turbine defaults establish the necessary length for this propeller.

6. **Screen aerodynamics before long aeroelastic runs.** Check periodic
   circulation, CL, CT, CQ, hub pitch/yaw moments, and phase-resolved waveforms.
   Monitor extreme induced velocities, close wake/surface encounters, and
   downstream wake distortion. Useful diagnostic ratios include core/chord,
   core/panel width, core/convected wake spacing, and encounter distance/core.
   None of these ratios has a universal acceptance threshold.

7. **Converge the stability quantities themselves.** Constant mean thrust or
   lift does not guarantee converged damping. A useful additional diagnostic
   is prescribed small pitch/yaw motion at the candidate mode frequency and
   evaluation of the incremental aerodynamic work, with baseline periodic
   forcing removed consistently:

   ```text
   W = integral(delta_M_pitch * pitch_rate + delta_M_yaw * yaw_rate, dt).
   ```

   Use compatible force/velocity coordinates and account for rotor-phase
   periodicity when averaging. With moments defined on the structure, positive
   net work adds structural energy. Then compare matched aeroelastic modes,
   decay/growth rates, frequencies, and bracketed flutter speeds. Absolute
   tolerances are needed near zero damping. Set a flutter-speed uncertainty
   target before the study; 1% is an example engineering target, not a standard.

Spatial refinement of the production wing also changes its beam mesh. Use
the rigid aerodynamic study first and keep the structural model fixed for
the final damping comparisons, or establish structural convergence separately.

## Trial radii through the existing interface

```julia
# Values inside the aerodynamic configuration for a diagnostic run.
core_radius_m = nothing
segment_core_factor = 0.0
chord_core_factor = 0.01
```

| Chord factor | Wing radius, c = 1.8 m | Blade radius, c = 0.197 m |
|---|---:|---:|
| 0.040 | 72 mm | 7.880 mm |
| 0.020 | 36 mm | 3.940 mm |
| 0.010 | 18 mm | 1.970 mm |
| 0.005 | 9 mm | 0.985 mm |

Both surfaces change in this table. Attribution to wing versus propeller
requires separate component settings; the current callback is shared by both.
The user-editable production factors were not changed by the twist
interpolation fix or by the addition of the fixed-radius option.

The existing rigid sweep already implements separate core, mesh, time, and
wake families. See [the workflow README](README.md) for environment overrides
and dry runs. It has its own physical defaults; align them with the case under
study. A rigid sweep passing its acceptance checks is only the first stage,
not validation of the flutter boundary.

## Physical wake-core modeling as a subsequent extension

For a physically calibrated wake, introduce distinct bound smoothing and
initial wake core, with separate wing/propeller settings. Use a consistent
radius on a shared vortex edge. Panel-specific radii can differ on twisted,
nonuniform, or deformed grids; equal opposing circulations on that edge then
need not cancel under different kernels.

One possible age-dependent model is

```text
rc(age)^2 = rc0^2 + 4 * alpha * nu_eff * age
```

Here age is in seconds, not rotor degrees. For the Lamb-Oseen convention
where core radius locates maximum swirl, `alpha = 1.25643`; `nu_eff` represents
the selected effective viscosity. This convention and coefficient must be
reconciled with the chosen kernel, rather than appended unchanged to a
different smoothing law. See the OLAF theory and Bhagwat-Leishman references.

Initial core and diffusion parameters should be checked against resolved CFD
or measured vortex profiles at appropriate wake ages and operating conditions.
If these data are unavailable, identify the core as a model assumption and
report its uncertainty in the flutter boundary. Increasing the core until
oscillations disappear does not establish physical stability.

If implemented in the partitioned solver, advance age/core from the accepted
wake state once per physical step; reset trial state during coupling iterations.
Any far-wake continuation or freezing should receive its own convergence check.

## Scope of the current verification

The twist correction passed 130 library checks (including 23 new twist checks)
and 41 aerodynamic convergence-integrity checks. Geometry-only diagnostics
produced the dimensions above. No production core sweep, re-trim, or flutter
boundary convergence is claimed by this note.

## Source comparison and additional kernel audit (2026-09-10)

The external implementations inspected were:

| Code | Revision / primary source | Core interpretation |
|---|---|---|
| Imperial UVLM | [`d8af34a`](https://github.com/ImperialCollegeLondon/UVLM/blob/d8af34a22baf1cddd38f1e362274c407637aab1c/include/biotsavart.h) | Hard exclusions; nonlinear `segment` and `UVLMlin` use different tests. Inspect the called routine before interpreting the input as a radius. |
| VortexLattice.jl | [`b89ecd0`](https://github.com/byuflowlab/VortexLattice.jl/blob/b89ecd089b3fac5737378afeecc73a47ed1006d7/src/induced.jl) | Smooth finite-core branch; the retrieved segment denominator loses its core term on the perpendicular-bisector plane. The local fork already corrects this. |
| PteraSoftware | [`0eed4d2`](https://github.com/camUrban/PteraSoftware/blob/0eed4d28117953a553ff5fa8a9f5d1c267fa649b/pterasoftware/_aerodynamics_functions.py) | Smooth denominator regularization with an initial radius and circulation-dependent wake-age growth. |

SHARPy exposes separate `vortex_radius` and `vortex_radius_wake_ind` settings;
its [default constant](https://github.com/ImperialCollegeLondon/sharpy/blob/main/sharpy/utils/constants.py)
is `1e-6`. These are exclusion tolerances, not interchangeable with the local
smooth core. Likewise, `near_field_force_model = :imperial` selects the local
load formulation and does not switch the Biot-Savart kernel to Imperial's.

Run the additional independent checks from the repository root:

```powershell
julia --startup-file=no --compiled-modules=existing --project=lib/WingPropellerUVLM lib/WingPropellerUVLM/examples/chang_linear_aeroelastic/studies/convergence/audit_core_kernels.jl
```

All 114 checks passed. The existing suite was also rerun: 41 convergence
integrity checks and 130 library checks passed. The new audit compares finite
segments with independent 100,000-point integration of the softened law at
interior, endpoint, exterior, and oblique locations. The maximum relative
error in these examples was `1.391e-9`. It also checks reversal, rotation,
length scaling, and the semi-infinite limit in two directions.

An independent curved-filament check evaluates a polygonal ring at its center
against the exact circular integral of the same softened kernel:

```text
v_z = Gamma * R^2 / (2 * (R^2 + rc^2)^(3/2)).
```

Error decreases over 16, 64, and 256 segments; at 256 it is below `5.02e-5`
relative for the three tested radii. This checks the center velocity, not
the difficult self-induced velocity on the curved filament itself.

For a unit straight segment, unit circulation, and radius 0.1, the audit
also reproduces the upstream bisector failure:

| Perpendicular distance | Local speed | Retrieved upstream speed |
|---|---:|---:|
| 0.01 | 0.154489 | 15.6034 |
| 0.001 | 0.0156048 | 156.064 |
| 0.0001 | 0.00156064 | 1560.64 |

Two limitations remain separate from the verified positive-core formula:

- `finite_core=true` with zero radius is singular at the filament. The
  unregularized branch also returns NaNs at segment interiors/endpoints.
  The Chang configuration rejects an all-zero core, and the callers have
  topological self-edge exclusions. A collapsed or intersecting geometry
  can still require additional detection; the kernel alone is not a guard.
- Shared edges can have different panel radii. For the corrected 5x5 Chang
  blade, the maximum chord-adjacent radius difference is approximately
  0.25662%. In a diagnostic with equal opposing unit circulations and an
  observation point 1 mm from each tested shared edge, the maximum residual
  is 0.52009% of the larger individual contribution. This is a constructed
  cancellation check, not a measured error in the production loads or
  flutter speed. A common edge radius and net edge circulation are the
  consistent treatment for a discrete vortex sheet.

The current mixed bound/wake regularization also needs a trailing-edge
cancellation check if the kernel policy is changed. Keep the singularity
guard, bound smoothing, and wake-core evolution explicit. No production
kernel or case setting was changed by this audit.
