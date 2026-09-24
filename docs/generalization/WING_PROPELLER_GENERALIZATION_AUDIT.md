# Wing-Propeller Generalization: Pre-Implementation Audit

This audit records the implementation state inspected before extracting the
configuration-independent structural and aeroelastic path. It covers the staged
working tree as found on 2026-09-21; those staged edits are treated as existing
user work.

## 1. Current execution path

The active high-level path is:

```text
chang_case_defaults
  -> load_chang_configuration
  -> build_chang_model_parameters
  -> assemble_structural_model
  -> build_chang_workspace / initialize_chang_uvlm
  -> solve_chang_aeroelastic!
  -> update_aero_geometry_for_state!
  -> assemble_structural_aero_load!
```

The Newmark-beta and generalized-alpha kernels are already configuration
independent in `src/aeroelastic`. The loose-explicit and implicit
predictor-corrector choices are also selected without a structural-case branch.
The configuration coupling is in the Chang example's model construction,
geometry adapter, and load adapter, not in those integration kernels.

## 2. Full-beam DOF convention

The implementation uses six structural coordinates at every wing node:

| Local index | Code name | Physical coordinate | Structural axis |
|---:|---|---|---|
| 1 | `u_span` | axial/span translation | +span |
| 2 | `v_chord` | in-plane/chord translation | +chord/aft |
| 3 | `w_vertical_down` | out-of-plane translation | +down |
| 4 | `theta_span` | torsion | rotation about +span |
| 5 | `theta_chord` | out-of-plane bending slope | rotation about +chord |
| 6 | `theta_vertical` | in-plane bending slope | rotation about +down |

The structural-to-aerodynamic basis conversion in
`chang_uvlm_coupling.jl` is

```text
(span, chord, down) -> (aero y, aero x, -aero z)
```

with the same axial mapping for rotations. The global coordinate order is all
wing-node coordinates first, followed by two coordinates per propeller:

```text
wing node 1 [1:6]
wing node 2 [1:6]
...
wing node N [1:6]
propeller 1 pitch
propeller 1 yaw
propeller 2 pitch
propeller 2 yaw
...
```

The root is currently clamped by removing the first six global coordinates.

### Bohnisch active-subspace mapping

The supplied Bohnisch element order is `[h, slope, alpha]` at each node. The
exact active embedding in the full beam is:

| General beam DOF | Physical quantity | Bohnisch treatment |
|---|---|---|
| 1, `u_span` | axial displacement | artificially stiff |
| 2, `v_chord` | in-plane displacement | artificially stiff |
| 3, `w_vertical_down` | out-of-plane bending displacement `h` | physical |
| 4, `theta_span` | torsion `alpha` | physical |
| 5, `theta_chord` | out-of-plane bending slope | physical |
| 6, `theta_vertical` | in-plane bending slope | artificially stiff through in-plane `EI` |

Thus the Bohnisch local element order `[h1,slope1,alpha1,h2,slope2,alpha2]`
maps to full-element indices `[3,5,4,9,11,10]`. This mapping follows the
implemented strain-displacement matrix and was not inferred from variable
names alone.

## 3. Chang-specific assumptions found

| File/function | Existing assumption | Required generalization |
|---|---|---|
| `chang_case.jl`, `chang_case_defaults` | Chang planform, mesh, propeller distribution, pylon modal data, advance-ratio law, and solver controls share one named case | move physical data to case files consumed by one interface |
| `chang_configuration.jl` | `CHANG_*` names and one Chang-only speed calculation | retain compatibility wrapper; resolve a physical speed-model option |
| `chang_model_parameters.jl`, `build_chang_model_parameters` | hard-coded Chang mass, inertia, CG, stiffness, and blade-angle tables | map the tables into the general stiffness and `:lumped_nodal` mass inputs |
| same | uniform beam mesh and one radius/chord/blade count/mesh shared by every propeller | accept general span nodes and per-propeller physical data |
| same | attachment specified by continuous `attachment_eta`, rounded afterward to a node | make `attachment_node` authoritative; derive every location from it |
| same | angular speed always preserves Chang's trim advance ratio | select `:fixed_omega` or `:constant_advance_ratio` per propeller |
| `chang_structural_matrices.jl`, `assemble_structural_matrices` | function is generically named in staged work but accepts Chang-shaped scalar propeller arguments and only a nodal mass path | extract one physical assembly with selectable wing mass model and per-propeller data |
| same | each propeller reuses identical `M_P`, `C_P`, `K_P`, inertia coupling, and `Omega` | assemble each propeller from its own fields |
| same | root BC is `7:end` | construct constrained/free coordinates from the six-DOF layout and configured root constraints |
| `chang_structural_model.jl` | diagnostics and modal classifier are Chang-named, though their mathematics is generic | provide generic diagnostics/modes and compatibility aliases |
| `chang_simulation.jl` | physical hub center and modal load center are derived from Chang pylon length | make aerodynamic hub the attachment node; keep any physical modal-coordinate coupling in the structural input, not as a second hub |
| `chang_uvlm_coupling.jl`, `initialize_chang_uvlm` | passes `propeller_span_positions` to the UVLM initializer | omit the independent position so node geometry is authoritative |
| same, `update_aero_geometry_for_state!` | staged work reads translation/rotation from one node, but the stored `attach_node_y` can still be the old continuous location | construct the pivot/hub from the selected node's reference and deformed position |
| same, `assemble_structural_aero_load!` | staged work transfers a wrench to one node, but still forms `cross(hub-pivot,F)` because hub and pivot are distinct | set hub equal to pivot/attachment and transfer the UVLM hub moment directly; retain only physical blade-to-hub reduction |
| `Initialization.jl`, `initialize_bohnisch_uvlm_system` | Bohnisch name, optional independent propeller span positions, and identical propeller geometry | add/use a configuration-neutral entry point and direct-node location; do not change the UVLM backend |
| `ChangAeroelastic.jl`, `ChangModel`, `run_chang` | model and analysis entry points are case-named | add a common model/analysis entry point; keep aliases for existing studies |
| validation scripts | several scripts still destructure the removed interpolation pair/weight fields | update tests to assert single-node support, hub colocation, and direct load support |

## 4. Actual Chang inertia formulation

The active Chang wing model is not a consistent distributed beam mass matrix.
The workbook supplies concentrated nodal masses, CG offsets, and inertia
tensors. `build_chang_model_parameters` converts those data to SI units, forms
six-by-six rigid spatial-inertia blocks, and conservatively remaps the blocks to
the active structural nodes by control volume. `assemble_structural_matrices`
then adds each complete block only to the diagonal block of its node.

Therefore the objective description "Chang is lumped" is correct for the
active wing discretization, with this qualification: each lump is a full rigid
spatial-inertia matrix, not merely three diagonal translational masses. It
contains translational mass, rotational inertia, first-moment coupling from the
CG offset, and inertia products. No wing element consistent mass is currently
assembled. Additional propeller/nacelle concentrated inertia enters separately
through `Gs_W`; the pitch/yaw inertia enters `Ms_P`; and their cross inertia
enters `Bs_P`/`Fs_W`.

This model must remain `:lumped_nodal`; replacing it by a distributed mass
model would not be a Chang regression.

## 5. Current stiffness formulation

The wing stiffness is already a full three-dimensional Euler-Bernoulli beam
representation. The four strain-resultant coordinates are axial strain,
in-plane curvature, out-of-plane curvature, and twist. The constitutive matrix
uses `EA`, `EIz`, `EIy`, optional `EIzy`, and `GJ`. Chang's tabulated properties
are integrated through every element and split at tabulated discontinuities.

The formulation has no transverse-shear coordinate and therefore no `GA`
input should be added. Artificial stiffness for reduced-reference cases belongs
only in `EA` and the excluded bending rigidity. Its scale must be configurable
and selected by a modal-convergence/conditioning study.

## 6. Current attachment and aerodynamic-location state

`HEAD` used two-node work-conjugate interpolation for concentrated structural
inertia, motion, and returned propeller wrench. The staged user changes already
replace those operations with a single attachment operator supported on one
node. However, `initialize_chang_uvlm` still passes the pre-rounding continuous
span positions to `initialize_bohnisch_uvlm_system`, which deliberately uses
those independent positions when supplied.

For the default 20-element Chang case, the inspected state is:

| Propeller | requested y [m] | selected node | node y [m] |
|---:|---:|---:|---:|
| 1 | 3.150 | 9 | 3.000 |
| 2 | 6.225 | 18 | 6.375 |

Consequently, direct structural attachment is only partially implemented: the
UVLM pivot can remain off-node. The original conclusion below that direct-node
attachment also required the aerodynamic hub to be colocated was incorrect.
For the Chang case, the selected node is the pylon pivot on the elastic axis;
the rotor hub remains one pylon length away. Its physical hub wrench must
therefore include `cross(hub-pivot,F)` when transferred to the wing node. This
distinction was restored after the generalization regression was identified.

Wing surface loads are a separate path. Chordwise panel forces are still
reduced at each aerodynamic span station and mapped to beam-node loads. That
wing aerodynamic-to-beam mapping is not part of the direct-propeller change.

## 7. Rotational speed, boundary conditions, and downstream solver

Chang currently computes

```text
mu = V_trim / (Omega_trim R)
Omega = V_inf / (mu R)
```

so its active law is constant advance ratio. There is no fixed-omega path in
the case model yet. The time step is then obtained from the configured azimuth
increment divided by `Omega`.

The first wing node is clamped and all appended propeller pitch/yaw coordinates
are free. This is currently implemented by a contiguous range rather than a
constraint list. The common Newmark-beta, generalized-alpha, loose-explicit,
and implicit predictor-corrector algorithms already operate only on reduced
`M`, `C`, and `K`; no case branch is needed downstream.

## 8. Units verified from assembly operations

| Quantity | Unit used by code |
|---|---|
| `EA` | N |
| `EIy`, `EIz`, `EIzy`, `GJ` | N m^2 |
| distributed mass | kg/m |
| distributed sectional torsional inertia `rhoJ` | kg m |
| nodal mass | kg |
| nodal/propeller inertia | kg m^2 |
| CG/reference offset | m |
| mass first moment | kg m |
| pitch/yaw stiffness | N m/rad (radian dimensionless numerically) |
| angular speed `Omega` | rad/s |
| damping ratio | dimensionless |
| stiffness-proportional damping coefficient | s |

## 9. Bohnisch distributed-inertia mapping

For the supplied reference matrix:

| Bohnisch quantity | Physical meaning | General representation |
|---|---|---|
| `rho_struct` | mass per unit span | `mass_per_length` [kg/m] in the consistent element mass |
| `rhoJ` | torsional mass moment of inertia per unit span | `torsional_inertia_per_length` [kg m] |
| `xalpha*sc` | signed chordwise CG offset from the beam reference/elastic axis | `cg_offset_chord` [m] |
| `rho_struct*xalpha*sc` | distributed first moment | bending-torsion cross block between `w_vertical_down`/`theta_chord` and `theta_span` |

With the supplied coordinates, `a=-0.30`, `e=-0.08`, `xalpha=0.22`, and
`xalpha*sc=+0.1375 m`: the CG is aft of the elastic axis when positive chord is
aft. Embedding the supplied six-by-six matrix into full element indices
`[3,5,4,9,11,10]` preserves its sign and all cross terms exactly. This is not
mathematically equivalent to Chang's nodal spatial-block mass and will be a
separate physical mass option, not a case-name branch.

## 10. Pre-generalization numerical baseline

For the staged default Chang case with 20 elements:

```text
free matrix size: 124 x 124
norm(M): 112.57948461239296
norm(C): 1536.3848889790759
norm(K): 1.9837309186277426e11
trace(M): 810.3556035919382
trace(C): 0.0
trace(K): 8.40046808500192e11
```

The first ten coupled undamped frequencies are:

```text
5.8286053255, 7.9285876657, 7.9543038276, 7.9677424283,
8.5817084495, 19.2160180311, 33.5507788576, 38.6435355067,
85.4548801612, 103.7172169672 Hz
```

These values establish the structural regression target for mapping Chang into
the new physical interface. Aerodynamic histories involving the now-prohibited
off-node hub are expected to change when exact colocation is enforced and must
be reported separately from the matrix regression.
