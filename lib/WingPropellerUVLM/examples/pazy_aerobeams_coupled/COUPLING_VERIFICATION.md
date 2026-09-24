# Pazy UVLM / AeroBeams verification

Audit date: 2026-09-24. Existing simulation files and settings were not modified.

## Conclusion

The implementation cannot currently be considered verified for post-flutter
response. Two reproducible coupling defects were found. Neither test establishes
which defect triggers the reported long-time crash, or guarantees that fixing
them will produce an LCO at 55 m/s. That requires a corrected, converged
long-duration calculation and matching the intended benchmark configuration.

## 1. Load transfer does not preserve virtual work

Locations: `PazyWingUVLMCoupling.jl:342` (`wing_geometry`) and
`PazyWingUVLMCoupling.jl:376` (`beam_loads`), especially line 385.

The geometry interpolates each neighboring beam section's chord points. The
moment transfer instead uses the arm from each structural node to the final
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

The current code substitutes `x_a(j,k) - r_n` for `d_n(k)`. It conserves total
force and moment, but that alone does not conserve virtual work. The compatible
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

This permits artificial interface energy exchange; its net sign over an actual
oscillation still needs measurement. It is not evidence that this numerical
value is the energy error in the user's simulation.

Remedy: make force transfer the transpose of the implemented motion mapping.
Alternatively, start with genuinely coincident span stations. Merely setting
`spanwise_panels=15` does not achieve coincidence: the structural Pazy stations
are nonuniform, while the current aerodynamic stations are uniform.

## 2. Newton recovery disconnects the external load updates

Locations: `src/SystemSolver.jl:112,184,188` and
`PazyWingUVLMCoupling.jl:64,67,213,220`.

Newton saves a deep copy of the model and can restore that copy after an
unsuccessful load step. The coupling continues updating BCs and reading geometry
through its original local `model` reference. In addition, deep-copying the BC
functions copies their captured `nodal_loads` array. Updating the original array
then no longer updates the loads seen by the restored structural model.

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

Remedy: always operate on `structure.model`, and explicitly maintain/rebind the
live load data after recovery. Changing only the local model reference is
insufficient because the load-function captures have also been copied. A simple
alternative is explicitly updating numeric BC values on the active model each
coupling iteration, provided all BC bookkeeping is refreshed consistently.

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
- Current settings use 4 by 15 aerodynamic panels, 40 retained wake rows, and
  `dt = 0.00044954545 s` at 55 m/s. The retained wake age is approximately
  0.01798 s. This audit does not establish wake-length or mesh convergence.
  The rule `U*dt=c/Nc` is a mesh-spacing choice, not a stability guarantee.
- The example has zero gravity and mirror symmetry. Match these assumptions,
  structural properties, angle of attack, and aerodynamic model before expecting
  quantitative agreement with a specific published LCO.

## Tests executed

```powershell
julia --project=lib/WingPropellerUVLM lib/WingPropellerUVLM/test/runtests.jl
julia --project=lib/WingPropellerUVLM lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/test_coupling.jl --dynamic
julia --project=lib/WingPropellerUVLM lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/verify_pazy_coupling.jl
```

- UVLM suite: 390 assertions passed.
- Existing Pazy checks: 18 geometry/motion assertions and 13 short dynamic
  assertions passed (97 steps, approximately 0.06 s at a final speed of 40 m/s).
- New audit: 7 assertions passed and 2 failed, reproducing the defects above.
  Its nonzero exit code is intentional while those defects remain; failures
  have not been suppressed or marked as expected passes.

Recommended order: correct the two coupling defects; rerun these tests; then
compare long-duration responses with independently refined time step, wake,
mesh, and FSI tolerances. Diagnose cycle work and amplitude envelopes before
adding numerical damping to force a bounded response. No long-duration
55 m/s LCO validation or reproduction of the user's exact crash was completed
in this audit.
