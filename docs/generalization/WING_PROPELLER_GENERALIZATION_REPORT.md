# Wing-Propeller Generalization Report

## Result

The structural and high-level aeroelastic path is now driven by physical case
data rather than a Chang/Xu/Bohnisch branch:

```text
case data
  -> full 3D beam stiffness
  -> selected physical wing mass model
  -> direct-node propeller/nacelle assembly
  -> common M, C, K
  -> common Newmark-beta/generalized-alpha integration
  -> common loose/implicit coupling
  -> configuration-neutral high-level UVLM adapter
  -> unchanged UVLM backend
```

The implementation is intentionally small. There is one structural assembly,
two explicit wing mass choices, one direct-node propeller assembly, and one
common analysis module. There is no configuration registry, type hierarchy,
or case-name conditional in the solver.

The pre-implementation findings are preserved separately in
`WING_PROPELLER_GENERALIZATION_AUDIT.md`.

## A. Chang-specific assumptions found and removed

| File/function | Assumption found | Generalization |
|---|---|---|
| `chang_case.jl`, `chang_case_defaults` | geometry, two propeller stations, operating law, and solver controls were one Chang-shaped input | Chang is mapped by `cases/chang_case.jl` into the same physical fields used by Xu and Bohnisch |
| `chang_model_parameters.jl`, `build_chang_model_parameters` | Chang tables directly fed the only assembly path | retained as a compatibility/data reader; its stiffness distributions and nodal spatial-inertia blocks are inputs to the general assembly |
| `chang_structural_matrices.jl` | scalar properties shared by all propellers and a nodal-only wing mass path | new `assemble_wing_propeller_structure` accepts per-propeller data and selects `:lumped_nodal` or `:consistent_distributed` by physics |
| same | root removal assumed the first six coordinates and every other coordinate was free | constraints are generated from configured wing node/local-DOF lists; appended propeller DOFs remain independent of mesh size |
| `chang_uvlm_coupling.jl`, initialization | requested continuous span station could override the selected structural node | the selected `attachment_node` is authoritative; the optional legacy station may only equal the node station |
| same, deformation | propeller attachment was historically interpolated between adjacent nodes | all translation and rotation are now read from one selected node |
| same, load transfer | a hub-to-pivot offset added an extra `r x F` | hub, pivot, and attachment are colocated; the physical UVLM moment about the hub is retained directly |
| `chang_simulation.jl` and configuration | `hub_load_arm_factor`, pylon-derived hub center, and modal load center represented separate points | the unused load-arm option was removed; the aerodynamic origin is the structural node; physical pylon effects remain in structural first moments/inertias |
| `Initialization.jl` | function name was Bohnisch-specific, propellers were geometrically identical, and an independent aerodynamic station was allowed | `initialize_wing_propeller_uvlm_system` is neutral and accepts independent per-propeller radius, chord, blade count, twist, panel counts, and wake limits |
| integration/coupling driver | the Chang example owned the only complete time-marching path | `run_time_domain_analysis` and `run_uvlm_time_domain_analysis` operate only on a physical case and common `M,C,K` |
| validation scripts | several scripts still expected node-pair/weight fields | direct-node assertions and virtual-work checks now use only `prop_attach_nodes` |

The validated Chang reader and named example remain for backward compatibility;
they are no longer the definition of the general solver.

## B. Final structural DOF convention

The actual six-DOF nodal order is:

| Local index | Coordinate | Physical meaning | Positive structural axis |
|---:|---|---|---|
| 1 | `u_span` | axial translation | span |
| 2 | `v_chord` | in-plane translation | chord/aft |
| 3 | `w_vertical_down` | out-of-plane translation | down |
| 4 | `theta_span` | torsion | about span |
| 5 | `theta_chord` | out-of-plane bending slope | about chord |
| 6 | `theta_vertical` | in-plane bending slope | about down |

The structural-to-aerodynamic transformation is

```text
(span, chord, down) -> (aero y, aero x, -aero z)
```

and is applied to translations and rotations. Global ordering is:

```text
all six-DOF wing nodes,
propeller 1 pitch, propeller 1 yaw,
propeller 2 pitch, propeller 2 yaw,
...
```

For the supplied Bohnisch order `[h,slope,alpha]`, the active full-beam element
indices are `[3,5,4,9,11,10]`:

| General DOF | Bohnisch treatment |
|---|---|
| axial translation, index 1 | finite high `EA` |
| in-plane translation/rotation, indices 2 and 6 | finite high in-plane `EI` |
| out-of-plane translation, index 3 | physical `EI` |
| torsion, index 4 | physical `GJ` |
| out-of-plane slope, index 5 | physical `EI` |

## C. General 3D beam stiffness

`general_beam_element_stiffness` integrates

```math
K_e = \int_{x_e}^{x_{e+1}} B(x)^T C(x) B(x)\,dx,
```

where the strain vector contains axial strain, in-plane curvature,
out-of-plane curvature, and twist. The constitutive data are `EA`,
`EI_in_plane`, `EI_out_of_plane`, optional `EI_coupling`, and `GJ`.
The existing Euler-Bernoulli formulation does not contain transverse shear,
so no unused `GA` parameter was introduced.

Every stiffness may be:

- one scalar;
- one value per element;
- one value per node;
- `(positions, values)` in metres; or
- `(eta, values)` in normalized span.

Properties are evaluated locally and stiffness integration is split at any
tabulated breakpoints.

For Bohnisch, excluded axial and in-plane motion uses configurable finite
values:

```text
EA = multiplier * EI_out_of_plane / span^2
EI_in_plane = multiplier * EI_out_of_plane
```

The selected multiplier is `1e3`, not an extreme sentinel. Section N records
the frequency/conditioning study that selected it.

## D. Wing mass formulations

### `:lumped_nodal`

The input is either a six-by-six spatial-inertia block at each node or the
physical nodal mass, CG offset, and inertia tensor. A block has the form

```math
M_n = \begin{bmatrix}
mI & -mS(r) \\
mS(r) & I_{CG}-mS(r)S(r)
\end{bmatrix}.
```

This retains translational mass, rotational inertia, mass first moments, and
products of inertia. It is the Chang model.

### `:consistent_distributed`

Each element assembles

```math
M_e=M_{translation}+M_{torsion}+M_{bending-torsion}.
```

It contains axial consistent mass, both Euler-Bernoulli transverse consistent
mass blocks, sectional torsional inertia, and the signed bending-torsion first
moment. It is the Bohnisch and Xu model. Scalar, element, nodal, or simple
tabulated spanwise inertia inputs use the same property interface.

Both branches return `M_wing`; downstream assembly and integration do not know
which branch produced it.

## E. Actual Chang inertia model

Chang is genuinely lumped at the active wing nodes, but each lump is a full
rigid spatial-inertia block. The workbook masses, CGs, inertia tensors, and
products of inertia are converted to SI units and conservatively remapped by
control volume. No consistent element wing mass is added.

The following remain separate:

- wing nodal spatial inertia: `M_wing`;
- concentrated nacelle/rotor inertia at the selected node: `M_attached`;
- pitch/yaw inertia: `M_propeller`;
- wing/propeller cross inertia: off-diagonal global mass blocks.

Thus Chang was not converted to the Bohnisch distributed formulation.

## F. Bohnisch inertia mapping

| Supplied quantity | Units | Physical meaning | General representation |
|---|---:|---|---|
| `rho_struct` | kg/m | mass per span | `mass_per_length` |
| `rhoJ` | kg m | torsional inertia per span | `torsional_inertia_per_length` |
| `xalpha*sc` | m | signed CG-to-elastic-axis chord offset | `cg_offset_chord` |
| `rho_struct*xalpha*sc` | kg | distributed first moment coefficient in the element integral | active bending-torsion cross block |

For the supplied coordinates, `a=-0.30`, `e=-0.08`, `xalpha=0.22`, and
`xalpha*sc=+0.1375 m`. Positive chord is aft, so the sign means that the CG is
aft of the elastic axis. This is embedded exactly in indices
`[3,5,4,9,11,10]`; it is not inferred from similarly named Chang fields.

## G. Bending-torsion coupling

`consistent_distributed_element_mass` inserts the full supplied coupling
matrix multiplied by

```math
\rho_{struct}l_e(x_\alpha s_c).
```

Numerical verification gives an active global mass relative error of
`1.727e-16` versus an independently assembled supplied Bohnisch matrix. In
particular, the `w_vertical_down`/`theta_span` terms are nonzero and are covered
by a unit test.

## H. Propeller attachment

The former representation had two neighboring structural nodes and two
weights. The final representation has one integer:

```julia
attachment_node = propeller.attachment_node
node_dofs = wing_node_dofs(attachment_node)
```

Every concentrated attached inertia, wing/propeller mass coupling, gyroscopic
coupling, motion extraction, and returned propeller wrench has support only on
that node. `attachment_operators` are explicit six-row Boolean-like operators,
which makes this support testable.

The general wing aerodynamic surface mapping remains separate. Surface
stations may interpolate to their enclosing beam element; propeller loads may
not.

## I. Aerodynamic propeller position

The neutral initializer derives span station, chord, and leading edge only
from `span_nodes[attachment_node]`. A legacy `propeller_span_positions` input
can only confirm this value; it cannot override it.

At every trial state, `direct_propeller_hub_positions` reads the translation of
the selected node. `general_uvlm_coupling.jl` then applies that node's rotation,
the propeller pitch/yaw coordinates, and prescribed shaft azimuth. There is no
adjacent-node interpolation and no independent aerodynamic origin.

## J. Force and moment transfer

Blade vertex forces are first reduced to a physical force and moment about the
hub. `colocated_hub_wrench` asserts that hub and attachment coordinates are
equal. The force and the already physical UVLM hub moment are then applied to
the selected structural node. No extra `cross(hub-attachment, force)` exists.

The Chang finite-difference virtual-work audit checked every free coordinate
for reference, uniform-rotation, and nonuniform-rotation states. Its maximum
scaled discrepancy was `9.25e-8` with a `1e-6` acceptance threshold.

## K. Rotational speed

`propeller_angular_speed` implements per propeller:

```text
:fixed_omega             -> omega_rad_s
:fixed_rpm               -> rpm * 2*pi/60
:constant_advance_ratio  -> pi*Vinf/(J*R)
```

Each propeller may select a different law and physical values. The resulting
speed is used by the gyroscopic `Ix*Omega` terms and aerodynamic shaft azimuth.

The reference cases now declare their intended physical laws directly:

- Chang and Xu use `:constant_advance_ratio`, so `Omega` changes with
  freestream speed while `J` remains fixed.
- Bohnisch uses `:fixed_rpm` with a default of 2500 RPM. Its configurable
  `propeller_beta75_deg` is added to a twist distribution whose value at
  `0.75R` is zero, making the resulting blade angle at `0.75R` equal to the
  requested `beta75_deg`.

## L. Chang regression

The general physical interface reproduces the former structural result exactly:

| Quantity | Relative error |
|---|---:|
| `M` | 0 |
| `C` | 0 |
| `K` | 0 |

The first ten coupled frequencies remain:

```text
5.828605, 7.928588, 7.954304, 7.967742, 8.581708,
19.216018, 33.550779, 38.643536, 85.454880, 103.717217 Hz
```

The 30-element isolated-wing comparison remains within 1.50% of all five
published modal targets, including `6.6305 Hz` for OOP1 and `42.8477 Hz` for T1.

A two-propeller, 33-step, direct-node UVLM response with a short pitch pulse
completed with finite history and all steps converged. Representative maxima
were:

| Response | Maximum absolute value |
|---|---:|
| tip displacement | 0.024293 m |
| tip twist | 1.91831 deg |
| propeller 1 pitch/yaw | 0.54191 / 0.09095 deg |
| propeller 2 pitch/yaw | 0.91758 / 0.09883 deg |

An entry-by-entry aerodynamic history match to the former implementation is
not an appropriate regression because the requested hub relocation deliberately
changes geometry and removes the former artificial moment arm. The structural
regression is exact; the new aerodynamic history is stored as a direct-node
validation artifact. No full flutter sweep was rerun because it was not a
readily bounded validation during this change.

## M. Xu verification

`cases/xu_case.jl` uses the same structural assembly, consistent distributed
mass, and direct-node propeller attachment. The published `1.6 m` propeller
station is made an exact node while retaining 20 elements.

| Mode set | General result [Hz] | Xu reference [Hz] |
|---|---|---|
| isolated wing | 2.8820, 16.6796, 18.0615, 46.5824, 50.5734 | 2.88, 16.68, 18.06, 46.59, 50.57 |
| coupled, zero RPM | 2.8480, 5.4784, 6.9978, 17.7453, 19.4879, 49.2815, 51.3924 | 2.85, 5.48, 7.00, 17.75, 19.48, 49.25, 51.37 |

The maximum discrepancies are below `0.015 Hz` and `0.04 Hz`, respectively.
The physical inputs were transcribed from Liu Xu, *Propeller-Wing Whirl
Flutter*, TU Delft (2020), Tables 9.1, 9.3, and 9.4.

## N. Bohnisch verification

The prompt did not provide a Bohnisch attachment station; the validation case
explicitly chooses the wing-tip node. It can be changed by editing only the
case's `attachment_node`.

The independently assembled reduced Bohnisch matrices and mapped active
full-beam matrices agree as follows:

| Matrix | Maximum absolute difference | Relative difference |
|---|---:|---:|
| mass | `3.997e-15` | `1.727e-16` |
| damping | `1.444e-8` | `3.121e-15` |
| stiffness | `1.431e-5` | `3.125e-15` |

The inactive-stiffness convergence study was:

| Multiplier | `cond(K)` | First five active state-space frequencies [Hz] | Max relative change from `1e5` |
|---:|---:|---|---:|
| `1e2` | `7.750e7` | 3.685735, 6.773117, 15.588069, 28.842747, 31.343079 | `3.217e-4` |
| `1e3` | `7.750e8` | 3.685771, 6.774760, 15.588565, 28.851532, 31.345186 | `2.609e-5` |
| `1e4` | `7.750e9` | 3.685774, 6.774921, 15.588609, 28.851982, 31.345273 | `2.362e-6` |
| `1e5` | `7.750e10` | 3.685775, 6.774937, 15.588613, 28.852028, 31.345278 | 0 |

`1e3` is the smallest tested value below the `1e-4` active-frequency criterion
and has 100 times better stiffness conditioning than `1e5`; it is the default.
Artificial high-frequency modes are deliberately excluded from this comparison.

## O. Final configuration interface

The three reference files are:

```text
examples/wing_propeller_aeroelastic/cases/chang_case.jl
examples/wing_propeller_aeroelastic/cases/xu_case.jl
examples/wing_propeller_aeroelastic/cases/bohnisch_case.jl
```

Routine run controls are centralized at the top of
`examples/wing_propeller_aeroelastic/run_case.jl`. The editable `settings`
tuple selects the reference case and overrides structural elements, independent
wing aerodynamic panels, propeller panels, attachment span fraction(s),
azimuth step or explicit time step, total duration, operating condition,
propeller speed law, time integrator, coupling scheme, interaction, finite-core
law, wake capacity, and inactive stiffness multiplier. The same block also
exposes the Newmark damping parameter, generalized-alpha spectral radius, and
the implicit-coupling iteration limit, tolerances, and relaxation. No reference
case file needs to be opened for these changes.

For a uniform `Ne`-element station definition, direct attachment is calculated
as `round(Int, eta*Ne)+1`. Xu additionally anchors that calculated node at the
published physical position, so its direct hub remains exactly at 1.6 m and
its reference modal validation is retained when using the default mesh.

A new configuration needs only the following shape (fields omitted here are
ordinary physical propeller coupling fields illustrated by the cases):

```julia
wing = (;
    span_nodes,
    stiffness = (;
        EA,
        EI_in_plane,
        EI_out_of_plane,
        EI_coupling = 0.0,
        GJ,
    ),
    mass_model = :consistent_distributed, # or :lumped_nodal
    inertia = (;
        mass_per_length,
        torsional_inertia_per_length,
        cg_offset_chord,
    ),
    damping = (; stiffness_coefficient = 0.0),
    boundary_conditions = (;
        constrained_nodes = [1],
        constrained_local_dofs = collect(1:6),
    ),
    geometry = (;
        root_chord,
        tip_chord,
        xle_root,
        xle_tip,
        elastic_axis_fraction,
    ),
)

propellers = [(;
    attachment_node,
    radius,
    blades,
    blade_chord,
    radial_panels,
    chordwise_panels,
    blade_twists,
    mass,
    spin_inertia,
    pitch_inertia,
    yaw_inertia,
    pitch_stiffness,
    yaw_stiffness,
    damping_ratio,
    speed_model = :fixed_rpm,
    rpm,
)]

case = (;
    wing,
    propellers,
    operating_condition = (; freestream_speed, air_density,
        angle_of_attack, sideslip),
    aerodynamic = (; spanwise_panels, chordwise_panels,
        mirror = false, symmetric = false),
    solver_options = (; time_integrator = :newmark_beta,
        newmark_alpha = 0.05,
        generalized_alpha_rho_infinity = 1.0,
        coupling_scheme = :implicit_predictor_corrector,
        coupling_maximum_iterations = 10,
        coupling_state_tolerance = 1.0e-5,
        coupling_load_tolerance = 1.0e-2,
        coupling_equilibrium_tolerance = 1.0e-10,
        coupling_coupled_equilibrium_tolerance = 1.0e-5,
        coupling_relaxation = 0.5),
)

structural = build_structural_model(case)
results = run_uvlm_time_domain_analysis(
    case;
    time,
    time_steps,
)
```

For a lumped wing, replace the distributed inertia tuple by
`spatial_inertia_blocks` or by nodal masses, CG offsets, and CG inertia tensors.
No solver function needs modification.

## Units

| Input | Unit |
|---|---|
| `EA` | N |
| both `EI`, `EI_coupling`, `GJ` | N m^2 |
| mass per length | kg/m |
| sectional torsional inertia per length | kg m |
| nodal/attached mass | kg |
| CG offset | m |
| mass first moment | kg m |
| `Ix`, pitch/yaw/nodal rotational inertia | kg m^2 |
| pitch/yaw stiffness | N m/rad |
| `Omega` | rad/s |
| RPM | rev/min |
| damping ratio | dimensionless |

## Verification executed

- Complete package suite: 388 assertions passed.
- Exact Chang general-interface `M,C,K` regression.
- 30-element Chang isolated-wing modal check.
- Direct attachment support, node motion, hub colocation, and no-offset wrench checks.
- Chang aerodynamic virtual-work audit, maximum scaled error `9.25e-8`.
- Xu isolated and coupled modal tables.
- Exact Bohnisch active-subspace matrix mapping.
- Bohnisch inactive-stiffness convergence and conditioning.
- Both time integrators crossed with both coupling algorithms on common `M,C,K`.
- A complete configuration-neutral UVLM time step, including transactional wake commit.

No file under `src/backend` was modified.
