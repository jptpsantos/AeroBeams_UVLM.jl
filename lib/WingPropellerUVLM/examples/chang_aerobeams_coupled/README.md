# Chang wing–propeller coupling with nonlinear AeroBeams

This example reads the existing Chang reference data and UVLM implementation,
and constructs a separate nonlinear AeroBeams 1.0 structural model. Use Julia
1.13 or newer and instantiate the shared `examples/Project.toml` environment
once; it resolves the checked-out AeroBeams and UVLM packages together.

## Run a case

From the repository root:

```powershell
julia +1.13 --project=lib/WingPropellerUVLM/examples lib/WingPropellerUVLM/examples/chang_aerobeams_coupled/run_chang_aerobeams.jl
```

Edit `chang_aerobeams_case.jl` for the physical configuration. Shared new solver,
initialization, pylon and output defaults are in
`src/CoupledSolver.jl::chang_aerobeams_defaults()`.

`reference` supplies Chang physical, mesh, flow and wake data. The legacy
`reference.integration`, `coupling`, `excitation` and `output` groups do not
control the new driver: use the top-level new `solver`, `initialization` and
`output` groups. In particular, `solver.impulse_*` defines the applied pulse.

For the nonlinear beam-pylon variation:

```powershell
julia +1.13 --project=lib/WingPropellerUVLM/examples lib/WingPropellerUVLM/examples/chang_aerobeams_coupled/run_chang_beam_pylon.jl
```

The example adds the shared environment to `LOAD_PATH` when launched from
another project. Running a case does not change any environment. In a restricted environment
where Julia cannot write its normal compiled cache, add
`--compiled-modules=existing` to the command. Initial compilation can take time.

The default 0.25 s case is a response demonstration, not a sufficient duration
for establishing a flutter boundary. Increase the duration and unperturbed
settling interval for damping identification.

## Two structural formulations

`solver.support=:reduced` retains two relative support coordinates per rotor,
pitch and yaw, on a fully nonlinear AeroBeams wing. Its objective reduced
kinetic energy preserves the current Chang reference mass blocks at zero
relative deformation, including both acceleration cross-blocks. A moving-CM
extension and SE(3) momentum transport retain finite wing-carrier motion.
Relative support deformation and intrinsic modal inertia remain a reduced,
small-relative-deformation approximation. The configured angle limit stops an
analysis outside its intended range.

`solver.support=:beam_pylon` replaces those two coordinates and their support
mass/stiffness with nonlinear AeroBeams pylon elements. A separate physical
rotor mass/inertia and spin momentum are applied at each physical tip. There is
no duplicate reduced support or `PointInertia` at the tip. The rotor follows
the pylon-tip orientation without independent bearing pitch/yaw coordinates.

The default beam-pylon properties need physical calibration:

- `EI_pitch=K_pitch*L` and `EI_yaw=K_yaw*L` match tip rotational compliance
  under a pure tip moment. They do not imply equal modal frequencies or
  effective inertias to Chang's calibrated reduced model.
- Pylon distributed mass uses the Chang mass-per-length datum. Its default
  sectional rotary inertia is an explicit configurable assumption.
- Rotor mass is blade count times blade mass. Its polar spin inertia uses
  Chang's calibrated `Ix_prop`; absent measured data, transverse inertia uses
  the thin axisymmetric disk assumption `J_transverse=Ix_prop/2`.
- Supply `pylon.rotor_inertia_tensor` or transverse inertia when measured
  rotor/nacelle data are available. `pylon.rotor_mass` and `spin_inertia`
  independently specify carried assembly mass and independently spinning polar
  inertia. Stationary nacelle mass is not invented.

The initial beam-pylon variation has no material damping. Nonzero reduced
support damping or global reference Rayleigh damping is rejected rather than
silently interpreted as a nonlinear beam constitutive law. The reduced
variation supports the reference relative viscous damping.

## Gyroscopic terms and coordinate conventions

All interface poses and physical wrenches use UVLM axes U=(chord aft, span, up).
Wing section directors b=(span, chord aft, down) map to U with

```text
R_U_b = [0 1 0; 1 0 0; 0 0 -1].
```

AeroBeams WM parameters are used through their actual rotation tensor. Parameter
rates are differentiated through that tensor to obtain physical angular
velocity and acceleration; they are never treated as Euler angles.

Chang support coordinates use `Rrel=Ry(theta)*Rz(-psi)` and prescribed blade
spin `Rx(-Omega*t)`. The added spin momentum is initially `[-Ix*Omega,0,0]`
in U. The reduced component retains wing–wing, wing–relative,
relative–wing and relative–relative gyroscopic blocks. The complete gyro matrix
is skew and does no mechanical work. In the beam-pylon variation the same
physics enters through the finite-rotation nodal moment

```text
M_inertia = J_world*alpha + omega × (J_world*omega + H_spin).
```

Ordinary AeroBeams inertia alone does not contain this prescribed independent
shaft spin. Setting aerodynamic RPM or adding a point inertia is insufficient.

## Nonlinear wing and the Chang reference

The wing is an actual AeroBeams mixed beam model with finite displacement and
rotation. UVLM supplies its aerodynamics; native `AeroSurface` loads are off.
Full 6×6 Chang spatial inertia blocks, including eccentric mass and products
of inertia, are conservatively reconstructed as element densities. Section
compliance is integrated across every tabulated stiffness break point.

This preserves the local mass data but differs from Chang's nodal lumping and
cubic Euler–Bernoulli element. The AeroBeams element has constant compliance
and uses a large finite shear stiffness. Its physical curvature convention
also requires transforming the sign of bending cross-stiffness. The legacy
Chang bending-slope convention should not be assumed to give identical matrix
entries under a single unexamined rotation sign map.

Run `validation/compare_wing_modes.jl` to export 8/16/32-element structural modal
comparisons. Differences are reported, not hidden by retuning the wing. These
are dry wing modes, not free-wake aeroelastic eigenvalues.

The current adapter supports a straight, constant-x wing elastic axis. A
different reference line must be implemented explicitly rather than silently
represented by the same straight beam. Wing nodes and pylon nodes have explicit
global maps; adding a pylon cannot make a wingtip output refer to its tip.

## Aerodynamic and structural load transfer

One global UVLM system contains the wing and every physical blade. Wing and
blade wakes are owned by that system. All interference is enabled by default;
wing image symmetry is disabled for the initially asymmetric installation.

The new interface preserves unique physical surface IDs separately from
interference groups. It also handles the facade's vector-symmetry constructor
locally, without changing the backend.

The initial adapter requires a fixed positive vortex core and the corrected
`:imperial` force model. Factor-only cores and legacy force formulas are
rejected explicitly rather than silently replaced.

Paired dimensional vortex-vertex forces and force-application positions are
used. The wing load map is the transpose of its spanwise motion interpolation,
including the quarter-panel vortex offsets. Each rotor's forces reduce to a
force and moment about its actual hub.

The reduced variation preserves Chang's distinct full physical hub offset and
half-length modal arm. Both its grid motion and generalized load projection
use the same effective modal center. The beam-pylon variation has no modal arm:
its complete wrench acts once at the actual pylon tip and is transmitted to
the wing by the beam assembly.

## Strong time coupling and initialization

The native nonlinear AeroBeams residual and analytic Jacobian are augmented
with the rotor component's attachment residual. Only the small attachment and
relative-coordinate block is centrally differentiated. Structural Newton
updates spin terms at each iterate. The aerodynamic load remains fixed during
that mechanical solve.

AeroBeams retains its first-order trapezoidal integration. Relative support
and nodal attachment parameter histories use compatible average-acceleration
Newmark updates. Every FSI trial uses the same accepted history. Configuration,
attachment velocity, aerodynamic load and final coupled mechanical residual
are checked before acceptance. Rotor azimuth comes from physical target time.
Circulation and wakes commit once; failed trials restore the accepted clocks
and structural/rotor state.

`initialization.mode=:mean_equilibrium` performs the following:

1. Hold structural geometry fixed and rotate/develop all wakes.
2. Check revolution-to-revolution mean loads, endpoint loads and circulation.
3. Solve the mean-loaded nonlinear wing/support equilibrium.
4. Repeat until loads and structural state converge.
5. Seed the dynamic system with that shape, developed wake and consistent phase.

Mean loads are retained throughout dynamics. No trim baseline is subtracted.
This is a mean equilibrium, not a fully coupled periodic structural orbit.
Allow subsequent dynamic settling before applying a stability perturbation.
Mean initialization requires an integer number of timesteps per revolution;
the default azimuth-based timestep satisfies that requirement.

A periodic aerodynamic state is a checked condition, not a guarantee. The
8×3 wing / 3×2 blade check with 72 propeller wake rows and a 5° azimuth increment
did not meet a 0.1% periodic tolerance in twelve revolutions. Its trim was
rejected; that tolerance was not relaxed. Use `mode=:unloaded` to study a
startup transient, or refine/settle the aerodynamic model before accepting
loaded periodic initialization. The verification record distinguishes those
two kinds of response.

`solver.speed_policy=:fixed_rpm` keeps configured RPM fixed when changing
airspeed. `:fixed_advance_ratio` retains the existing Chang speed/RPM scaling.
Prescribed-speed torque is external to the reduced shaft description; this
version does not solve motor, governor or windmilling shaft dynamics.

## Verification and stability study

Run all focused tests:

```powershell
julia --compiled-modules=existing --project=. lib/WingPropellerUVLM/examples/chang_aerobeams_coupled/validation/runtests.jl
```

The tests cover reference mass/gyro blocks, analytical dry whirl frequencies,
spin reversal, finite-rotation covariance, Galilean invariance, kinetic power,
physical tip inertia, frozen-force virtual work, force/moment balance, separate
span meshes, UVLM rollback, augmented Jacobian derivatives, converged dynamic
steps in both support formulations, and mean-loaded initialization.

`run_speed_sweep.jl` provides a separate longer study. It runs independently
trimmed baseline/perturbed pairs at each airspeed and fits RMS perturbation
growth. Initial phases must match. Fit quality is saved alongside the exponent;
a rotor's forced blade-passage response is not itself evidence of flutter.
Repeat with mesh, timestep, wake length, core and perturbation refinement before
interpreting a sign change as a flutter boundary. Floquet linearization and
nonlinear branch continuation are not implemented by this study helper.

`validation/run_verification_responses.jl` runs longer unloaded-start transients
in both structural formulations, retaining a complete revolution of propeller
wake and applying a smooth support-moment pulse. Those are numerical response
checks, not trim, damping or flutter-boundary measurements.

Generated CSV histories, settings and PNG responses go only to this example's
`output/`. The settings record the physical and numerical assumptions. CSV
pitch/yaw denote relative support coordinates for the reduced case and shaft
orientation in the fixed U frame for the beam-pylon case.

## Source map

| File | Responsibility |
|---|---|
| `src/AeroBeamsWing.jl` | Chang beam properties, explicit nodal maps, finite rotations and attachment kinematics |
| `src/ReducedRotor.jl` | Two-DOF objective moving-base support and all gyroscopic blocks |
| `src/BeamPylon.jl` | Connected nonlinear pylon assembly |
| `src/TipRotor.jl` | Physical tip inertia and independent spin momentum |
| `src/AerodynamicInterface.jl` | Multi-surface grids, conservative load transfer and UVLM construction |
| `src/CoupledSolver.jl` | Augmented Newton, strong FSI, initialization, histories and output |
| `src/Studies.jl` | Phase-matched perturbation growth study |
