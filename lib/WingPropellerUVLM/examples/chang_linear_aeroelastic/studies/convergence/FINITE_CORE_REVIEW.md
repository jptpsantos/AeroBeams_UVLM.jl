# Chang aerodynamic code and finite-core review

Scope: the Chang rigid and partitioned aeroelastic drivers, convergence
handoff, geometry/core assignment, induced-velocity kernels, and existing
structural/load-transfer validation paths. This is not a formal verification
of every routine in the repository.

## Errors corrected

| Finding | Consequence | Correction |
|:--|:--|:--|
| Rigid defaults were 65 m/s and 3 degrees; sweep expected 80 m/s and 0 degrees | Valid outputs could fail validation and stop the batch | Shared options |
| Reuse omitted duration, averaging, collective pitch, density, sideslip, shedding fraction, and source revision | Different problems could be compared or reused | Complete TOML metadata, source/result fingerprints, sample-count and time/phase checks |
| Means determined acceptance; finest case always agreed with itself | Changing oscillations or unresolved refinement could pass | Cycle waveform RMS, final-pair agreement, and all-finer-level checks |
| Fixed segment factor changed radius under span/radial refinement | Mesh and regularization errors were mixed | Fixed chord-based radius during ordinary mesh sweeps; separate core sensitivity |
| Downstream gate assumed three core levels and zero chord factor | Incomplete four-level study could pass, or chord core be discarded | Expected level counts and both core factors carried to children |
| Importing sweep executed `main()` unconditionally | Helper imports could launch a batch | Restore program-file guard |
| Standalone plots wrapped azimuth through zero | Curves connected unrelated cycle endpoints | Continuous rotor revolutions |
| Trailing kernel softened perpendicular distance but not endpoint distance | Inconsistency with long finite-segment limit near the core | Same softened endpoint distance in both kernels |

The last kernel issue affects callers enabling semi-infinite trailing
vortices. Chang sets `system.trailing_vortices .= false`, so that correction
does **not** explain Chang's time-step sensitivity. Earlier workspace
corrections to the finite-segment denominator and rotation-invariant edge
length were retained and tested.

## Meaning of the core

The radius assigned by `grid_to_surface_panels` is:

```text
rc = max(segment_factor * norm(rtr - rtl), chord_factor * c)
```

The edge length is the panel's three-dimensional spanwise/radial bound edge.
`c` is full local surface chord. One radius is stored per panel and used on
its ring edges. A shed wake panel inherits its trailing-edge surface panel's
radius; convection preserves it. There is no viscous core spreading or
stretching model here.

This radius is a numerical regularization parameter. Specifying it does not
by itself model viscous diffusion. OpenFAST's theory distinguishes those
effects and offers separate evolution models.
[OpenFAST OLAF theory](https://openfast.readthedocs.io/en/main/source/user/aerodyn-olaf/OLAFTheory.html)

The finite-segment kernel integrates the softened Biot-Savart integrand
`dl × r / (|r|^2 + rc^2)^(3/2)`. The trailing kernel now also uses
`sqrt(|r|^2 + rc^2)` in its endpoint term. Tests compare against independent
quadrature, the long-segment limit, circulation reversal, geometric scaling,
and the vanishing-core limit away from a filament.

## Choosing a radius for this study

Hold radius fixed while refining space and time first. The new sweep's
nominal `rc = 0.01 c` gives **0.018 m on the wing** and **0.00197 m on the
blades** for chords of 1.8 m and 0.197 m. These are explicit trial values,
not experimentally established viscous cores. The interactive aeroelastic
runner's editable mixed law remains available.

The chord-factor family 0.04, 0.02, 0.01, 0.005 spans wing radii of
72, 36, 18, 9 mm and blade radii of 7.88, 3.94, 1.97, 0.985 mm. Seek a range
where periodic CL/CT/CQ and their waveforms change less than the selected
tolerances. Repeat that core test on a finer mesh before using a value for
damping calculations. A trend toward smaller radius is insufficient evidence.

OpenFAST's lifting-line/free-wake guidance currently suggests a
span-spacing-based radius and factor near 0.6, with a broad anticipated range.
This supports testing mesh-based regularization, but its factor should not
be transferred directly to this panel UVLM and kernel.
[OpenFAST OLAF guidance](https://openfast.readthedocs.io/en/main/source/user/aerodyn-olaf/RunningOLAF.html)

Straight-segment representations of curved vortices can also retain core and
curvature errors under some refinement procedures. Reducing the core is not
a general cure for distorted wakes or poor convergence.
[Van Hoydonck, Gerritsma, and van Tooren, 2012](https://arxiv.org/abs/1204.2699)

Each rigid summary reports actual minimum/maximum radii by surface and maximum
ratios to span/radial edge and panel chord. Inspect these ratios for anisotropic
meshes. The separate `segment` mode explicitly couples mesh and core changes.

## Remaining limits

- **Wake truncation:** Chang discards old free-wake rows without semi-infinite
  continuation. Retained wake length remains an independent error source.
  Shedding fraction eta controls placement, not iterative wake relaxation.
- **Panel-level cores:** adjacent panels on nonuniform/tapered grids may
  assign different radii to a shared edge. A future edge-based model should
  make shared-edge induction consistent. This review does not redesign panel
  storage or add viscous diffusion.
- **Structural span coupling:** production span refinement also changes the
  beam. The rigid study isolates aerodynamics; the damping gate keeps the
  validated 30-element structural mesh.
- **Geometry:** rigid options explicitly define Chang geometry and do not
  import arbitrary edits to `chang_case.jl`. Align both before studying a
  modified physical model.
- **Reference accuracy:** agreement with a tested reference is not a formal
  uncertainty estimate. No Richardson/GCI extrapolation is claimed: several
  refinement ratios are unequal and an asymptotic regime is unproven.
- **Mean trim:** subtracting the mean leaves periodic aerodynamic forcing.
  Damping still needs an appropriate free-response interval and modal tracking;
  passing rigid coefficient checks does not validate flutter damping.

## Validation evidence

Automated checks pass: 107 solver tests, 41 sweep-integrity tests, and 11
aerodynamic-to-aeroelastic gate tests. The reduced sweep also reuses all three
unchanged cases successfully and regenerates the comparisons and plots.

Structural verification passes its published-frequency tolerances on 30
elements; the largest listed isolated-wing frequency error is about 1.50%.
Mass/stiffness positivity and symmetry checks pass. The virtual-work audit
at 4-degree propeller pitch and 3-degree yaw gives maximum scaled error
`5.76e-9`.

A deliberately coarse end-to-end sweep (wing 4x2, blades 2x2, 30/20/15-degree
steps, one retained wake revolution, five simulated revolutions) completes
all children. Mean CT values are approximately `-0.03021`, `-0.01379`, and
`-0.007080`: clearly unresolved. The revised checks reject these cases.
These are workflow checks, not production settings.

No full production-resolution aerodynamic sweep or long damping study is
claimed converged. The accompanying README gives commands for establishing
that evidence at the operating point of interest.
