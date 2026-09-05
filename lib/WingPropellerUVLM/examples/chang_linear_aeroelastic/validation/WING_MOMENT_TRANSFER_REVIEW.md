# Investigation: work-conjugate wing moments

The discrepancy is confirmed for the geometry currently implemented. The
structural-to-aerodynamic sign conversion is correct; the missing operation
is projection of spatial moments onto the instantaneous axes of the wing's
Euler coordinates. It occurs in both direct wing loads and propeller wrenches
distributed to wing attachment nodes.

This investigation adds a standalone diagnostic and validates an analytic
candidate. It does not change the production load mapping or structural model.

## Derivation for the actual rotation order

Let the structural wing rotations be:

- alpha: span rotation, structural DOF 4;
- beta: chord rotation, structural DOF 5;
- gamma: down-axis rotation, structural DOF 6.

Their aerodynamic components are theta_y = alpha, theta_x = beta,
theta_z = -gamma. Both the wing grid and propeller attachment use

```text
R = Rz(-gamma) Rx(beta) Ry(alpha).
```

For a point x = p + R r, a coordinate variation produces
delta x = a_i cross (x-p) delta q_i. Therefore

```text
Q_i = sum(F dot dx/dq_i) = a_i dot M,
M   = sum((x-p) cross F).
```

The axes expressed in the aerodynamic basis are

```text
a_span  = ( sin(gamma) cos(beta), cos(gamma) cos(beta), sin(beta) )
a_chord = ( cos(gamma),          -sin(gamma),           0         )
a_down  = ( 0,                   0,                  -1         ).
```

Consequently the work-conjugate generalized moments are

```text
Q_span  = sin(gamma) cos(beta) Mx + cos(gamma) cos(beta) My + sin(beta) Mz
Q_chord = cos(gamma) Mx - sin(gamma) My
Q_down  = -Mz.
```

The current code uses (My, Mx, -Mz). These expressions agree when beta and
gamma are zero, including pure spanwise torsion at nonzero alpha. A nonzero
rotation by itself does not imply an error. The outer and middle rotations
change the axes associated with the inner coordinates.

The third coordinate remains conjugate to -Mz exactly for this Euler order.
Do not rotate all moments into the full body frame: body-axis components
are also not generally the generalized moments conjugate to Euler angles.

Relevant sources:

- [Wing geometry](../../../src/wing_propeller/GridUtilities.jl), line 86:
  Rz(theta_z) Rx(theta_x) Ry(theta_y).
- [Chang geometry/load adapter](../src/chang_uvlm_coupling.jl), lines 91–99,
  145–147, 273–275, and 319–336.
- [Paired force application points](../../../src/backend/nearfield.jl),
  imperial_nodal_positions. The audit differentiates these vortex vertices,
  so it does not mix physical-grid points with vortex-grid force locations.

## Propeller attachment chain rule

For an attachment between two structural nodes,

```text
q_attachment = w_left q_left + w_right q_right
Q_node,rotation = w_node A(q_attachment)' M_about_pivot.
```

Evaluate the axes at the interpolated attachment angles, then distribute
the projected moments with the same interpolation weights. Using each
node's own axes on a weighted attachment moment would fail for nonuniform
rotations. Direct wing loads instead use each wing node's angles and moment
about its own deformed elastic axis.

The propeller's pitch/yaw modal projections already use their instantaneous
axes. Those components pass the current audit and need no analogous change
for this particular issue.

## Numerical evidence

The diagnostic explicitly reproduces the earlier review setup: 84 m/s,
interaction enabled, 4x2 wing panels, 2x2 panels per blade, and both core
factors 0.025. These are isolated process settings; the current editable
0.01 defaults were not changed. The short run supplies aerodynamic forces
at approximately 0.00953 s, before trim or structural excitation.

Nodal aerodynamic forces are frozen while differentiating the geometry.
The constructed reference propeller pitch/yaw angles are 4/3 degrees.
Central differences use steps of 1e-5, 1e-6, and 1e-7 (radians for rotation,
metres for translation).

For uniform 1-degree wing rotations, the original attachment result is:

| Right attachment, chord rotation | Generalized moment (N m) |
| --- | ---: |
| Current mapping | 20.7960458652 |
| Analytic candidate | 17.9386876274 |
| Geometry finite difference, h=1e-6 | 17.9386876227 |

The required correction is -2.8573582378 N m:

| Source | Correction (N m) |
| --- | ---: |
| Direct wing loads | -0.4440731223 |
| Propeller attachment wrench | -2.4132851155 |

About 84.5% of this component correction comes from the propeller attachment.
The reported 13.74% is the component's scaled discrepancy, with denominator
max(abs(current), abs(FD), 1 N m). It is not an overall force or response error.

The extended check includes all 26 free coordinates, not only the attachment
and propeller coordinates checked previously. At node 2, which has no
propeller attachment, a pure 1-degree down-axis rotation gives a spurious
chord-rotation generalized moment of 0.5069174933 N m in the current mapping.
The candidate and independent derivative both give zero to numerical
precision. At zero span rotation, a variation of the middle chord rotation
does not move a point lying on the chord axis; the nonzero Cartesian moment
therefore is not work-conjugate to this coordinate.

Eleven scenarios cover zero rotations, uniform amplitudes of 0.001 through
5 degrees, each pure rotational component, nonuniform wing rotations and
translations, and zero propeller modal angles. All 26 coordinates are checked
at all three FD steps: 858 component comparisons. The maximum candidate
scaled discrepancy is 5.28666e-8. Pure span torsion passes the existing mapping.
The diagnostic completed 72 projection assertions and two near-zero
tangent assertions.

Larger scaled errors also occur at other wing coordinates with very small
reference moments. Those values should be read with their dimensional
moments, not presented as percentages of the total aerodynamic load.

## Why this can matter near zero displacement

For small beta and gamma,

```text
Q_span,exact - Q_span,current   = gamma Mx + beta Mz + O(angle^2)
Q_chord,exact - Q_chord,current = -gamma My + O(angle^2).
```

If a reference aerodynamic moment M0 is nonzero, these are first-order
load terms. Subtracting a constant mean trim load does not remove a
state-dependent axis correction.

With forces frozen at the short-run loaded reference and all propeller
modal angles reset to zero, the maximum missing directional tangent is
325.641 N m/rad for a direction increasing all free wing rotations equally.
Halving the small perturbation reproduces the tangent. This is a load-map
directional derivative, not a modal stiffness, flutter-speed shift, or
damping prediction. The sampled forces are startup loads, not a converged
trim solution.

The implemented structure remains linear. Making the aerodynamic load map
work-conjugate to its existing finite geometry does not introduce a complete
geometrically nonlinear beam or prove that every prestress term required
for a particular linearization is present.

## Recommended correction and validation

For the existing finite geometry, apply the analytic moment projection at
every wing node and at each interpolated propeller attachment. Keep the
current force components, moment reference points, and propeller modal
axes. The candidate uses three scalar projections; numerical differentiation
is only a validation tool and need not enter the time-marching solver.

If the intended model instead uses strictly linearized geometry, derive
that displacement map explicitly and transfer loads through its transpose.
Keeping exact geometry with an unexplained fixed-axis moment map leaves
the virtual-work mismatch demonstrated here.

After applying either modeling choice, the production audit should assert
all wing and propeller coordinates at nonzero, nonuniform states. Then compare
identical coupled-response cases, recording full wing rotations and
generalized loads. The ordinary response CSV contains only tip torsion and
propeller angles; it cannot establish the size of this error throughout
the user's current full run. No damping or flutter change is quantified
by the frozen-load investigation.

## Reproduction

[Investigation script](investigate_chang_wing_moment_transfer.jl)

```powershell
julia --startup-file=no --compiled-modules=existing --project=lib/WingPropellerUVLM lib/WingPropellerUVLM/examples/chang_linear_aeroelastic/validation/investigate_chang_wing_moment_transfer.jl
```

The script calls the existing audit, compares the production mapping to the
candidate and independent geometry derivatives, and writes:

- [All component comparisons](../output/wing_moment_investigation/components.csv)
- [Scenario summaries](../output/wing_moment_investigation/summary.csv)
- [Omitted directional tangent](../output/wing_moment_investigation/omitted_preload_tangent.csv)

These checks validate the proposed projection. They must not be interpreted
as assertions that the unchanged production mapping is correct.
