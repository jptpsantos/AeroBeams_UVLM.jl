# Pazy UVLM / AeroBeams verification

Audit and correction dates: 2026-09-24 and 2026-09-25. The corrections change
only the Pazy example coupling and visualization; both solver backends remain
unchanged.

## Coordinate-system correction (2026-09-25)

The structural and wake animations exposed a physical convention error: the
UVLM chord was mirrored about the spar. AeroBeams' Pazy airfoil points forward
along positive structural y, whereas a UVLM chord fraction increases from the
leading edge to the trailing edge. The coupling now uses the right-handed map
`(x_U, y_U, z_U) = (-y_A, z_A, -x_A)` and constructs chord offsets as
`(spar_fraction - chord_fraction) * chord`.

This correction puts the physical leading edge at `+spar_fraction*chord` in
AeroBeams y, the trailing edge at `-(1-spar_fraction)*chord`, and maps positive
UVLM lift to the Pazy model's upward direction. The same corrected transform is
used for forces and moments, so kinematics and load transfer remain work
conjugate. New regression checks explicitly test the leading edge, trailing
edge, transform handedness, and freestream direction.

Results generated before this correction used the reversed chord and must be
rerun; changing only the visualization cannot repair those saved results.

## Conclusion

Both reproduced defects are corrected and the expanded coupling regression
tests pass. This verifies the specific transfer and recovery properties, not
the complete post-flutter response. Neither test establishes which defect
triggered the reported long-time crash or guarantees an LCO at 55 m/s. That
requires a converged long-duration calculation and matching the intended
benchmark configuration.

## Applied code changes

In `PazyWingUVLMCoupling.jl`:

- `wing_geometry` now returns `(grid, positions, vortex_offsets)`. The offsets
  belong to each individual structural node, are resolved in AeroBeams axes,
  and include the UVLM quarter-panel shift and physical trailing-edge row.
- `beam_loads(loads, vortex_offsets, weights)` now uses those offsets in the
  nodal moments. Translational force weights are unchanged. This preserves
  force, moment, and virtual work for the implemented kinematic mapping without
  requiring coincident meshes.
- The FSI loop predicts, relaxes, and accepts `vortex_offsets` with exactly the
  same coefficients as the aerodynamic grid. The final load evaluation uses
  offsets from the final structural geometry.
- Load BCs now contain six numeric values, not functions capturing an external
  array. The new `apply_pazy_nodal_loads!` writes them into `problem.model` and
  refreshes the BC data before each structural solve and after FSI acceptance.
  Its simple BC ordering assumption is explicit: clamp first, then nodal loads.
- Geometry and tip histories always read `structure.model`, including after
  initialization and Newton recovery. No changes to AeroBeams' Newton algorithm
  or its model-copy/retry behavior were needed.

In `verify_pazy_coupling.jl`, the original failing checks now use the corrected
transfer API and the same numeric-BC update helper as the real coupling. Added
checks cover finite rotations, force/moment balance, quarter-panel and trailing
edge offsets, FSI prediction/relaxation, and a second load solve after an actual
Newton recovery. Chordwise meshes of 1, 4, and 7 panels are tested.

This document records the original evidence below and the verification results.

## 1. Original defect: load transfer did not preserve virtual work

Locations: `wing_geometry` and `beam_loads` in `PazyWingUVLMCoupling.jl`.

The geometry interpolates each neighboring beam section's chord points. The old
moment transfer instead used the arm from each structural node to the final
interpolated aerodynamic point. These are not work-conjugate mappings.

Writing the geometry of an aerodynamic vortex vertex as

```text
x_a(j,k) = sum_n W(j,n) * (r_n + d_n(k)),
```

where `d_n(k)` is the rotated chord offset at structural node `n`, the compatible
transfer for spatial nodal forces and moments is

```text
F_n = sum_(j,k) W(j,n) * f_a(j,k)
M_n = sum_(j,k) W(j,n) * cross(d_n(k), f_a(j,k)).
```

The old code substituted `x_a(j,k) - r_n` for `d_n(k)`. It conserved total
force and moment, but that alone did not conserve virtual work. The compatible
offsets must correspond to the vortex vertices (including their quarter-panel
shift and trailing edge), not blindly to the physical wing-grid vertices.

Reproduction using the actual Pazy nodal locations and the current uniform
15-panel aerodynamic span mesh: apply a 1 N vertex force and virtually rotate
one neighboring section about its chord, with its beam position held fixed.
The aerodynamic points do not move, but the transferred beam moment does work:

```text
Aerodynamic work per unit virtual rotation:  0.0 N*m
Structural work per unit virtual rotation: -0.007902345442707792 N*m
Total force and moment balance:             PASS
Coincident-span-mesh control:               PASS
```

This permitted artificial interface energy exchange; its net sign over an actual
oscillation still needs measurement. It is not evidence that this numerical
value is the energy error in the user's simulation.

Applied remedy: make force transfer the transpose of the implemented motion
mapping. The corrected value in the test above is 0.0 N*m, matching aerodynamic
work. Coincident span stations are not required. Merely setting
`spanwise_panels=15` does not achieve coincidence: the structural Pazy stations
are nonuniform, while the current aerodynamic stations are uniform.

## 2. Original defect: Newton recovery disconnected the external load updates

Locations: `src/SystemSolver.jl:112,184,188` and
the load BC construction and FSI loop in `PazyWingUVLMCoupling.jl`.

Newton saves a deep copy of the model and can restore that copy after an
unsuccessful load step. The old coupling continued updating BCs and reading
geometry through its original local `model` reference. In addition, deep-copying
the old BC functions copied their captured `nodal_loads` array. Updating the
original array then no longer updated the restored structural model's loads.

The audit exercised a real retry in the same Newton routine, using a small
nonlinear cantilever with the coupling's captured-array BC pattern:

```text
Model replaced during recovery: true
Newton finally converged:        true
Updated external force:          31.0 N
Force read by restored model:    30.0 N
```

This is conditional on recovery; it does not explain growth that occurs before
any recovery. The reproduction uses a structural steady problem to isolate the
shared Newton routine, not the user's entire 55 m/s trajectory.

Applied remedy: always operate on `structure.model` and explicitly update numeric
BC values and cached data on that model each coupling iteration. The corrected
test confirms 31.0 N both in the active model BC and its special-node BC after
recovery, followed by a converged second solve with a changed structural state.

## Other checks and remaining limits

- The aerodynamic trial starts from a copy of accepted circulation, geometry,
  and wake. Wake advancement occurs at commit, not repeatedly inside FSI.
- Structural equivalent rates are formed once per physical step, not once per
  coupling iteration. The underlying AeroBeams update is trapezoidal on its
  first-order states (`Problem.jl:1510`, `Core.jl:155`); this Pazy wrapper does
  not implement a separate adjustable Newmark-beta algorithm.
- UVLM uses backward differences for surface motion and circulation rate
  (`backend/circulation.jl:8`, `backend/analyses.jl:808`). Converging FSI does
  not eliminate their time-discretization error. A time-step study is necessary
  for damping and LCO amplitude/frequency.
- The 2026-09-24 audit run used 4 by 15 aerodynamic panels, 40 retained wake rows, and
  `dt = 0.00044954545 s` at 55 m/s. The retained wake age is approximately
  0.01798 s. This audit does not establish wake-length or mesh convergence.
  The rule `U*dt=c/Nc` is a mesh-spacing choice, not a stability guarantee.
- The example has zero gravity and mirror symmetry. Match these assumptions,
  structural properties, angle of attack, and aerodynamic model before expecting
  quantitative agreement with a specific published LCO.

## Test commands

```powershell
julia --project=lib/WingPropellerUVLM lib/WingPropellerUVLM/test/runtests.jl
julia --project=lib/WingPropellerUVLM lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/test_coupling.jl --dynamic
julia --project=lib/WingPropellerUVLM lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/verify_pazy_coupling.jl
julia --project=lib/WingPropellerUVLM lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/test_static_uvlm.jl
```

Initial audit, before corrections:

- UVLM suite: 390 assertions passed.
- Existing Pazy checks: 18 geometry/motion assertions and 13 short dynamic
  assertions passed (97 steps, approximately 0.06 s at a final speed of 40 m/s).
- New audit: 7 assertions passed and 2 failed, reproducing the defects above.
  These failures were not suppressed or marked as expected passes.

After the 2026-09-25 coordinate correction:

- Expanded coupling verification: 35 assertions passed, including the explicit
  leading-edge, trailing-edge, handedness, and freestream-direction checks.
- Existing Pazy checks rerun: all 20 geometry/motion and 13 short dynamic
  assertions passed (97 accepted steps, approximately 0.06 s).
- Static UVLM/AeroBeams checks: all 22 assertions passed. At 15 m/s and 3
  degrees the corrected coarse-mesh smoke test converged to positive upward tip
  displacement, consistent with the Pazy and UVLM axis definitions.
- The earlier 55 m/s smoke test was run before the coordinate correction. It is
  retained only as historical information and is not validation of the corrected
  model; post-flutter cases must be rerun.

Next: compare long-duration responses with independently refined time step, wake,
mesh, and FSI tolerances. Diagnose cycle work and amplitude envelopes before
adding numerical damping to force a bounded response. No long-duration
55 m/s LCO validation or reproduction of the user's exact crash was completed
in this audit.
