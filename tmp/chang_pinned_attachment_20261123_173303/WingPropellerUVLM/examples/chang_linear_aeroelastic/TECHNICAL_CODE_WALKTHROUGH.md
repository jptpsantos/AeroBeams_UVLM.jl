# Chang wing--propeller aeroelastic solver: complete code walkthrough

This document describes the code that is actually executed by the Chang time-domain example. It is intended to support an independent physical, mathematical, and numerical audit, not just explain how to run the example.

The following labels are used throughout:

- **Verified from code**: follows directly from the active implementation, and where stated from an executed test.
- **Inferred model interpretation**: the mathematics makes this interpretation likely, but the source does not name it explicitly.
- **Concern**: behavior is confirmed in the code, but its physical intent or robustness needs review.

Paths in this document are relative to `examples/chang_linear_aeroelastic` unless they begin with `../../src`, in which case they refer to the `WingPropellerUVLM` package source.

## 1. Scope and the normal entry point

The normal executable is [`run_chang_linear_aeroelastic.jl`](run_chang_linear_aeroelastic.jl). It activates `lib/WingPropellerUVLM`, loads [`src/ChangAeroelastic.jl`](src/ChangAeroelastic.jl), obtains one resolved configuration, and calls `run_chang`.

The main path is:

```text
chang_case.jl: editable defaults
        |
        v
load_chang_configuration(): defaults + CHANG_* overrides + derived settings
        |
        v
build_chang_model()
        |-- build_chang_model_parameters(): geometry, reference data, Omega, time grid
        |-- assemble_chang_structural_model(): M, C, K and attachment operators
        `-- chang_aerodynamic_options(): finite core and reference/load points
        |
        v
build_chang_workspace()
        `-- initialize_chang_uvlm(): wing/blade panels, System, empty wake, buffers
        |
        v
solve_chang_aeroelastic!()
        |-- initialize U0, Udot0, Uddot0
        |-- for every physical step n -> n+1
        |     |-- save accepted aerodynamic state at n
        |     |-- evaluate full UVLM trial(s) at target time t[n+1]
        |     |-- map UVLM vertex loads to structural generalized loads
        |     |-- Newmark-beta or generalized-alpha structural correction
        |     |-- loose exchange once, or implicit fixed-point iterations
        |     `-- after acceptance, convect/shed the wake exactly once
        `-- retain structural histories and coupling diagnostics
        |
        v
write_chang_results(): CSV, summary, optional plot/animation
```

The two numerical choices are independent:

- `config.integration.time_integrator` selects `:newmark_beta` or `:generalized_alpha`.
- `config.coupling.scheme` selects `:loose_explicit` or `:implicit_predictor_corrector`.

Thus all four combinations use the same structural matrices, geometry mapping, UVLM operator, load mapping, trim logic, and wake commit. Only the structural quadrature and the number/order of aerodynamic--structural exchanges change.

## 2. Actual call hierarchy

### 2.1 Model construction and run control

```text
run_chang_example()
`-- load_chang_configuration()
    |-- chang_case_defaults()
    |-- environment_overrides()
    |-- chang_rotor_angular_speed()
    `-- validate_configuration()
`-- run_chang(config)
    |-- load_chang_visualization()
    `-- run_chang_case(config)
        |-- build_chang_model(config)
        |   |-- validate_configuration()
        |   |-- build_chang_model_parameters()
        |   |   |-- chang_spatial_inertia_block()
        |   |   |-- chang_control_volume_spatial_remap()
        |   |   |-- chang_properties_from_spatial_inertia()
        |   |   `-- get_nodal_properties_chang()
        |   |-- assemble_chang_structural_model()
        |   |   `-- assemble_chang_structural_matrices()
        |   `-- chang_aerodynamic_options()
        |-- report_and_validate_structural_model()
        |-- build_chang_workspace()
        |   `-- initialize_chang_uvlm()
        |       `-- initialize_bohnisch_uvlm_system()
        |           |-- wing_to_grid()
        |           |-- generate_propeller_blades_grid()
        |           |-- initialize_propeller_grids!()
        |           |-- grid_to_surface_panels()
        |           `-- System(...)
        |-- chang_wake_context()
        |-- chang_excitation_options()
        |-- chang_integration_options()
        |-- solve_chang_aeroelastic!()
        `-- write_chang_results()
```

### 2.2 One aerodynamic trial

```text
aero_load_for_state!(snapshot, Utrial, step)
|-- restore_uvlm!(system, snapshot)
|-- update_aero_geometry_for_state!(Utrial, t[n+1])
|   |-- generate_panel_grid_and_interpolate()       # wing
|   |-- RotationMatrix()                            # wing/whirl/spin
|   `-- grid_to_surface_panels()                    # replace current panels
|-- propagate_system!(..., advance_wake=false)
|   |-- get_surface_velocities!()
|   |-- update_wake_shedding_locations!()
|   |-- influence_coefficients!()
|   |-- update_trailing_edge_coefficients!()
|   |-- normal_velocity!()
|   |-- circulation!()                              # AIC * Gamma = w
|   |-- dGamma/dt = (Gamma[n+1] - Gamma[n]) / dt
|   `-- near_field_forces!()
|       |-- induced_velocity()
|       |-- imperial_segment_force()
|       `-- imperial_unsteady_panel_force()
`-- assemble_structural_aero_load!()
    |-- imperial_nodal_forces()
    |-- imperial_nodal_positions()
    `-- chang_wing_generalized_moment()
```

### 2.3 Structural and coupling call tree

```text
solve_chang_aeroelastic!()
|-- loose_explicit_aeroelastic_step()
|   `-- structural_corrector(method)
|       |-- newmark_beta_corrector(), or
|       `-- generalized_alpha_corrector()
|
`-- partitioned_aeroelastic_step()
    |-- structural_effective_system()               # factor once per step
    |-- structural_predictor()
    |-- repeatedly call aero_load_for_state!()
    |-- structural_corrector(method)
    |-- structural_kinematics(method)
    `-- structural_equilibrium_residual(method)

after successful loose/implicit solution:
`-- commit_chang_wake!()
    `-- advance_wake!()
        |-- get_wake_velocities!()
        `-- shed_wake!()
            |-- translate_wake!()
            `-- rowshift!()
```

## 3. File-by-file map

### 3.1 Example-owned files

| File | Main responsibility and important definitions | Called by / dependencies | Category |
|---|---|---|---|
| [`run_chang_linear_aeroelastic.jl`](run_chang_linear_aeroelastic.jl) | User entry point; `run_chang_example()` activates the package and starts one case. | Loads `ChangAeroelastic`; calls configuration and run functions. | Analysis control |
| [`chang_case.jl`](chang_case.jl) | `chang_case_defaults()` returns the editable nested `NamedTuple`: wing, propeller, flow, wake, excitation, integration, coupling, output, and structural data. | Read only by `chang_case_defaults()` in the module. | Input |
| [`src/ChangAeroelastic.jl`](src/ChangAeroelastic.jl) | Defines module `ChangAeroelastic`, immutable `ChangModel`, `build_chang_model`, `run_chang`, and `run_chang_case`; includes the active example helpers. | Top-level script and studies. Depends on exported package routines and all files below. | Architecture/control |
| [`src/chang_configuration.jl`](src/chang_configuration.jl) | Parses typed `CHANG_*` overrides, resolves automatic wake rows/output paths/impulse start, computes the fixed-advance-ratio rotor speed, validates settings. | `run_chang_example`, studies, `build_chang_model`. | Input/configuration |
| [`src/chang_model_parameters.jl`](src/chang_model_parameters.jl) | Converts Chang tables to SI; builds structural/aero geometry, conservative spatial-inertia remap, distributed stiffness data, pylon/rotor properties, attachment weights, reference quantities, `Omega`, `dt`, and `t`. | `build_chang_model`. Uses package interpolation, blade data, `Freestream`, `Reference`. | Physical model/data |
| [`src/chang_structural_matrices.jl`](src/chang_structural_matrices.jl) | Beam strain matrix and stiffness integration; nodal spatial masses; propeller/pylon modal, inertial and gyroscopic blocks; global `M`, `C`, `K`; clamped-root reduction. | Included by `chang_structural_model`. | Structural model |
| [`src/chang_structural_model.jl`](src/chang_structural_model.jl) | Passes parameter data to the matrix assembler; provides isolated-wing modal analysis and matrix diagnostics. | `build_chang_model`, run diagnostics, validation scripts. | Structural model/verification |
| [`src/chang_uvlm_coupling.jl`](src/chang_uvlm_coupling.jl) | Initializes the case-specific UVLM system; maps structural state to wing/blade panels; maps UVLM vertex forces back to wing and propeller generalized loads; defines the deterministic aerodynamic trial. | `build_chang_workspace`, time marcher. | Coordinates, mesh/load mapping |
| [`src/chang_workspaces.jl`](src/chang_workspaces.jl) | `build_chang_workspace()` allocates fresh mutable aerodynamic and load-mapping storage for a model. | `run_chang_case`. | State allocation |
| [`src/chang_simulation.jl`](src/chang_simulation.jl) | Converts settings into aerodynamic, excitation, integration and animation options; initializes structural history; performs the complete accepted-step loop; trim subtraction; wake commit; safety checks. | `run_chang_case`. Depends on coupling and package integrators. | Time marching/coupling |
| [`src/chang_postprocessing.jl`](src/chang_postprocessing.jl) | Extracts tip and propeller channels, writes history CSV and summary, invokes optional plotting and animation. | `run_chang_case`. | Output |
| [`src/chang_plotting.jl`](src/chang_plotting.jl) | Plots accepted displacement/angle histories. Loaded only if requested. | `load_chang_visualization`, postprocessing. | Visualization |
| [`src/chang_animation.jl`](src/chang_animation.jl) | Copies accepted surface/wake frames, computes display limits, and writes a wake GIF. Loaded only if requested. | Simulation callback and postprocessing. | Visualization |
| [`src/chang_propeller_trim.jl`](src/chang_propeller_trim.jl) | Separate rigid-propeller RPM trim utility. It simulates an isolated propeller and root-finds windmilling or thrusting RPM. | Trim study scripts only; **not called by the main coupled run**. | Separate propeller study |

The `studies/` and `validation/` trees construct configurations and call the same model/run functions. They do not provide a second production time marcher. `studies/convergence/archive/` contains historical drivers and is not on the normal call path.

### 3.2 Package files on the active path

| File | Active role | Important definitions |
|---|---|---|
| [`../../src/WingPropellerUVLM.jl`](../../src/WingPropellerUVLM.jl) | Package include/export hub. | Loads backend, wing/propeller helpers, transaction support, time integrators, excitation, and generic problem API. |
| [`../../src/wing_propeller/BladeGeometry.jl`](../../src/wing_propeller/BladeGeometry.jl) | Chang blade reference distribution. | `get_twist_deg_chang`, `get_nodal_properties_chang`. |
| [`../../src/wing_propeller/GridUtilities.jl`](../../src/wing_propeller/GridUtilities.jl) | Interpolation and physical-grid construction. | `linear_interpolate_1d`, `generate_panel_grid_and_interpolate`, `generate_propeller_blades_grid`. |
| [`../../src/wing_propeller/Kinematics.jl`](../../src/wing_propeller/Kinematics.jl) | Accepted/current geometry history helper. | `copy_surfaces_to_previous!` is active. The older generic wing/propeller update helpers are not used by the Chang adapter. |
| [`../../src/wing_propeller/Initialization.jl`](../../src/wing_propeller/Initialization.jl) | Allocates and initializes the combined wing/blade `System`, interaction groups, wake capacity, and load buffers. | `initialize_bohnisch_uvlm_system`. |
| [`../../src/backend/rotors.jl`](../../src/backend/rotors.jl) | Supplies `RotationMatrix(angle, axis)`. The package's generic rotor generator is not used here. | `RotationMatrix`. |
| [`../../src/backend/panel.jl`](../../src/backend/panel.jl) | Immutable aerodynamic primitive types and their geometric accessors. | `SurfacePanel`, `WakePanel`. |
| [`../../src/backend/geometry.jl`](../../src/backend/geometry.jl) | Wing physical grid and vortex-lattice panels. | `wing_to_grid`, `grid_to_surface_panels`, `repeated_trailing_edge_points`. |
| [`../../src/backend/reference.jl`](../../src/backend/reference.jl) | Dimensional reference area/chord/span/point/speed/density. | `Reference`. |
| [`../../src/backend/freestream.jl`](../../src/backend/freestream.jl) | Freestream vector and optional body-frame rotation velocity. | `Freestream`, `freestream_velocity`, `rotational_velocity`. |
| [`../../src/backend/system.jl`](../../src/backend/system.jl) | The mutable-by-contained-arrays UVLM state/storage object. | `System`, `PanelProperties`. |
| [`../../src/backend/induced.jl`](../../src/backend/induced.jl) | Regularized Biot--Savart kernels, bound/wake induced velocities, and AIC blocks. | `bound_induced_velocity`, `ring_induced_velocity`, `influence_coefficients!`, `induced_velocity`. |
| [`../../src/backend/circulation.jl`](../../src/backend/circulation.jl) | Mesh-relative velocity and the no-penetration right-hand side. | `get_surface_velocities!`, `normal_velocity!`, `circulation!`. |
| [`../../src/backend/analyses.jl`](../../src/backend/analyses.jl) | Orders one unsteady UVLM solve and optionally its wake update. | `propagate_system!`, `advance_wake!`. |
| [`../../src/backend/nearfield.jl`](../../src/backend/nearfield.jl) | Production Imperial-style dimensional bound-segment and unsteady panel loads; transfer to vortex vertices. | `near_field_forces!`, `imperial_nodal_forces`, `imperial_nodal_positions`. |
| [`../../src/backend/legacy_nearfield.jl`](../../src/backend/legacy_nearfield.jl) | Selectable compatibility force model. | `legacy_imperial_segment_forces!`; active only when configured. |
| [`../../src/backend/wake.jl`](../../src/backend/wake.jl) | Shedding location, free-wake velocity, convection, circulation assignment, and fixed-capacity row rotation. | `update_wake_shedding_locations!`, `get_wake_velocities!`, `translate_wake!`, `shed_wake!`. |
| [`../../src/UVLMState.jl`](../../src/UVLMState.jl) | Snapshot/restore transaction for repeated aerodynamic trials. | `UVLMSnapshot`, `snapshot_uvlm`, `restore_uvlm!`. |
| [`../../src/aeroelastic/StructuralTimeIntegration.jl`](../../src/aeroelastic/StructuralTimeIntegration.jl) | Integrator dispatch and independent loose/implicit coupling. | Newmark routines, `structural_*`, `loose_explicit_aeroelastic_step`, `partitioned_aeroelastic_step`. |
| [`../../src/aeroelastic/GeneralizedAlpha.jl`](../../src/aeroelastic/GeneralizedAlpha.jl) | Generalized-alpha equations, residual scaling, and coupling options. | `generalized_alpha_*`, `PartitionedCouplingOptions`. The older `partitioned_generalized_alpha_step` is retained but the Chang path uses the generic method-dispatched step. |
| [`../../src/aeroelastic/Excitations.jl`](../../src/aeroelastic/Excitations.jl) | Smooth generalized-force pulse. | `smooth_hann_pulse`, `smooth_hann_pulse_load`. |
| [`../../src/UVLMProblem.jl`](../../src/UVLMProblem.jl) | Generic aerodynamic-only transactional problem facade. | `UVLMDynamicProblem`, `begin_time_step!`, `evaluate_trial!`, `commit_time_step!`. **The Chang coupled solver does not call this facade**; it implements the analogous transaction explicitly so it can insert the structural fixed-point loop. |

Backend files for VSP import, nonlinear lifting-line analysis, far-field drag, stability derivatives, and VTK output are loaded by the package but are not called by this analysis.

## 4. Function behavior and data flow

### 4.1 Configuration and model creation

| Function | Inputs and return | Mutation/state effect | Role |
|---|---|---|---|
| `chang_case_defaults` | No inputs; returns a fresh nested settings tuple. | The vector fields are newly created. | Single editable source of normal case settings. |
| `environment_value`, `environment_flag`, `environment_number`, `environment_optional_number`, `environment_symbol` | Setting name/default and environment dictionary; return typed value. | None. | Strict parsing of overrides. |
| `environment_overrides` | One settings group plus binding table. | None. | Applies only declared overrides. |
| `chang_rotor_angular_speed` | Resolved propeller and simulation settings. | None. | Preserves the trim inflow ratio when airspeed changes. |
| `load_chang_configuration` | Defaults and optionally an environment dictionary. | Deep-copies defaults; does not retain `ENV`. | Resolves wake rows, impulse time, output path, damping reference, then validates. |
| `validate_configuration` | Complete configuration. | None. Throws on invalid input. | Guards geometry, meshes, force-model names, pulse, integrator and coupling options. |
| `build_chang_model_parameters` | Resolved configuration. | Prints diagnostics only. | Produces a large read-only named tuple used by structural and aerodynamic builders. |
| `build_chang_model` | Configuration. | Deep-copies it. | Creates `ChangModel(config, parameters, structural, aerodynamic_options, force_function)`. |
| `build_chang_workspace` | `ChangModel`. | Allocates all mutable arrays. | Makes separate runs from the same model independent. |

The rotor relationship is implemented as

\[
\mu=\frac{V_{\rm trim}}{\Omega_{\rm trim}R},\qquad
J=\pi\mu,\qquad
\Omega=\frac{V_\infty}{R\mu}.
\]

`rotation_rpm` is therefore a trim-speed datum, not a constant RPM for every case. `build_chang_model_parameters` recomputes `Omega` at the simulation airspeed, holding `mu` and `J` constant. The time step is

\[
\Delta t=\frac{\Delta\psi_{\rm azimuth}}{\Omega}.
\]

The time range is `range(0, end_time, step=dt)`, so the integrated final time is the last grid point not exceeding the requested end time.

### 4.2 Structural helpers

| Function | Operation | Mutation |
|---|---|---|
| `chang_spatial_inertia_block` | Forms the complete 6-by-6 spatial inertia about the beam reference axis from mass, CG inertia and CG offset. | None. |
| `chang_control_volume_spatial_remap` | Integrates a piecewise-linear spatial-inertia density over target nodal control volumes. Moves the constrained-root block to the first free node. | Creates new blocks; source untouched. |
| `chang_properties_from_spatial_inertia` | Recovers mass, offset and CG inertia for diagnostics. | None. |
| `chang_point_mass_matrix` | Forms the nodal translational/rotational mass matrix including first moments and parallel-axis shift. | None. |
| `chang_beam_strain_displacement_matrix` | Returns the 4-by-12 axial, two-bending-curvature and torsional-strain matrix `B`. | None. |
| `chang_beam_element_stiffness_matrix` | Integrates `B' C B` for constant `C` by two-point Gauss quadrature. | None. |
| `chang_constitutive_matrix_at_eta` | Interpolates `EA`, `EIz`, `EIy`, `EIzy`, `GJ`. | None. |
| `chang_distributed_beam_element_stiffness_matrix` | Splits an element at every tabulated property break and integrates `B' C(x) B` by three-point Gauss quadrature. | None. |
| `assemble_chang_structural_matrices` | Assembles wing and pylon/propeller block matrices, attachment operators, gyroscopic terms, Rayleigh damping, and removes the root DOFs. | Allocates and returns matrices. |
| `assemble_chang_structural_model` | Supplies the active parameter set to the matrix routine. | None. |
| `chang_wing_modal_analysis`, `chang_wing_modal_frequencies` | Solve the isolated fixed-root wing generalized eigenproblem and classify modal strain-energy families. | None. |
| `chang_structural_diagnostics` | Checks symmetry, eigenvalues and frequencies. | None. |

### 4.3 Geometry and mapping helpers

| Function | Exact action | Persistent effect? |
|---|---|---|
| `linear_interpolate_1d` | Piecewise-linear scalar/vector interpolation with endpoint clamping. | No. |
| `generate_panel_grid_and_interpolate` | At each span station, rotates a chord line about its elastic-axis point using `Rz*Rx*Ry`; returns a newly generated physical grid. | No. Despite its name, the current Chang mesh has one structural node per aerodynamic span station, so no separate span interpolation is performed. |
| `generate_propeller_blades_grid` | Creates one twisted constant-chord blade and rotates copies evenly about rotor `x`. | No. |
| `grid_to_surface_panels` | Converts physical panel-edge vertices to quarter-chord vortex rings, three-quarter-chord collocation points, normals, chord and finite-core radius. | Returns new immutable panels. |
| `update_aero_geometry_for_state!` | Reconstructs root DOFs, maps bases, generates the wing, composes attachment/whirl/spin rotations, generates every blade, and replaces `system.surfaces`. | Yes, trial geometry in `system` and workspace hub/axis arrays. Safe in an implicit iteration only because the snapshot is restored first. |
| `chang_wing_generalized_moment` | Projects an aerodynamic moment onto the instantaneous work-conjugate axes of the wing Euler sequence. | No. |
| `assemble_structural_aero_load!` | Converts dimensional vortex-vertex loads into the reduced structural generalized-force ordering. | Reuses/overwrites workspace force and moment buffers; does not advance physical state. |
| `aero_load_for_state!` | Restores the accepted aerodynamic snapshot, installs target geometry, solves circulation/loads with no wake advance, and returns generalized loads. | Leaves `system` holding the latest trial, deliberately. |

### 4.4 Time integration and run helpers

| Function | Role | Mutation |
|---|---|---|
| `chang_aerodynamic_options` | Chooses fixed or geometry-scaled core radius and defines pivot, physical hub, and modal load points. | None. |
| `chang_excitation_options` | Converts rotor speed to revolution period and defines the trim averaging interval. | None. |
| `chang_integration_options` | Independently validates integrator and coupling scheme; builds both parameter sets, coupling tolerances and dimensional residual scales. | None. |
| `initialize_structural_history` | Allocates histories; sets `U0=0`, `Udot0=0`, and computes `Uddot0=M\(F0-C Udot0-K U0)`. | New arrays only. |
| `build_chang_coupling_scales` | Creates translation/angle and force/moment reference scales for fixed-point residuals. | None. |
| `smooth_hann_pulse_load` | Applies `0.5(1-cos(2*pi*phase))` to selected propeller pitch coordinates. | Returns a new load vector. |
| `solve_chang_aeroelastic!` | Owns structural histories, trim averaging, coupling selection, accepted/rejected step semantics, wake commit and diagnostics. | This is the main persistent state-changing routine. |
| `commit_chang_wake!` | Convects/sheds once, then increments each active row count up to capacity. | Permanently changes accepted wake geometry/circulation and active counts. |
| `chang_state_limit_status` | Checks finite values and configured mixed-component bounds. | None. |
| `write_chang_results` | Extracts and writes accepted histories. | Creates output directory/files; solver state unchanged. |

## 5. Coordinate systems and signs

### 5.1 Structural wing frame `S`

**Verified from code:** each node is ordered

```text
[u_span, v_chord, w_vertical_down,
 theta_span, theta_chord, theta_vertical_down_axis]
```

- `X_S`: spanwise, root to tip.
- `Y_S`: chordwise, leading edge to trailing edge/downstream.
- `Z_S`: vertically downward.
- This is right-handed for the modeled side: span cross chord gives down.
- The beam element coordinate is the same span axis; element DOFs are the six components at the left node followed by the six at the right node.

### 5.2 Aerodynamic body frame `A`

- `x_A`: chordwise/downstream.
- `y_A`: spanwise, root to modeled tip.
- `z_A`: upward.
- The origin for the undeformed wing is its root leading-edge plane; the elastic axis is at `xle + elastic_axis_fraction*chord`.
- `Freestream` produces
  `Vinf*[cos(alpha)cos(beta), -sin(beta), sin(alpha)cos(beta)]`.

The basis transformation in `update_aero_geometry_for_state!` is

\[
[x_A,y_A,z_A]^T=[v_S,u_S,-w_S]^T,
\]

and the same permutation/sign is applied to the rotation components. The conjugate force mapping in `assemble_structural_aero_load!` is therefore

\[
[F_{u_S},F_{v_S},F_{w_S}]^T=[F_{y_A},F_{x_A},-F_{z_A}]^T.
\]

These are the principal `X` permutation and `Z` sign changes in the active path.

### 5.3 Wing rotation convention

At a span station the aerodynamic chord is transformed by

\[
R_W=R_z(\theta_{z_A})R_x(\theta_{x_A})R_y(\theta_{y_A}).
\]

The structural rotational coordinates are not simply projected onto fixed aerodynamic axes at finite angle. `chang_wing_generalized_moment` uses the virtual-work axes of this Euler chain:

- structural span rotation: `Rz*Rx*e_y`;
- structural chord rotation: `Rz*e_x`;
- structural vertical/down rotation: `-e_z`.

This is why the mapping needs the current `theta_x_A` and `theta_z_A`.

### 5.4 Propeller and blade frames

- The reference rotor axis is aerodynamic `x_A`.
- Blade 1 is constructed with radial direction initially along `+y_A`; `z_A` completes the disk-plane basis. Other blades are rotations about `x_A`.
- The blade chord at each radius is placed by the Chang twist table; the stored angle is `twist_deg - 90 deg + collective_offset`.
- The pivot is on the deformed wing elastic axis at the selected structural node.
- The physical rotor hub is one pylon length from that pivot along `-x_A`.
- The modal load point uses `hub_load_arm_factor * pylon_length`; the default Chang assumed-mode factor is `0.5`.
- Pylon/nacelle offsets are represented consistently in both the aerodynamic geometry/load map and the structural first moments and cross inertias.
- Propeller pitch is positive about the wing-rotated `+y_A` axis.
- Structural propeller yaw is positive about structural down, hence the aerodynamic yaw angle is its negative and its work-conjugate axis is `-R_W R_y(pitch)e_z`.
- Prescribed spin is `R_x(-Omega*t)`. Azimuth is computed from absolute target time; no mutable azimuth counter exists.

## 6. Structural model

### 6.1 State and beam theory

The reduced structural equation is

\[
M\ddot U+C\dot U+KU=F_{\rm aero}+F_{\rm external}.
\]

The code mapping is:

| Mathematical quantity | Code |
|---|---|
| `M`, `C`, `K` | `model.structural.M`, `.C`, `.K` |
| `U_n` | `displacement[step]` |
| `Udot_n` | `velocity[step]` |
| `Uddot_n` | `acceleration[step]` |
| aerodynamic perturbation load | `accepted_perturbation_load`, `candidate_load`, or `explicit_load` |
| external pulse | `external_load_n`, `external_load_np1` |

`chang_beam_strain_displacement_matrix` contains linear axial strain, linear torsion, and cubic-Hermite bending curvatures in two planes. There is no shear-strain row or independent shear DOF. **Inferred model interpretation:** this is a linear Euler--Bernoulli-type beam with axial/torsional deformation and coupled two-plane bending, rather than a Timoshenko beam.

The structure itself is linear: matrices do not depend on `U`. The aerodynamic geometry applies finite Euler rotations, so the overall formulation mixes a linear structural model with nonlinear kinematic placement. This is normally defensible only while the structural response remains in the intended small-deformation range.

### 6.2 Mass

The wing does not use a conventional distributed element mass matrix. Chang's nodal mass, CG, and full inertia tables are converted to 6-by-6 spatial-inertia blocks and conservatively remapped to the active mesh. These blocks are lumped at structural nodes.

The remap transfers mass, first moments, products of inertia, and rotational inertia together. Because the root node is removed by the clamp, its remapped block is deliberately added to the first free node. This preserves total inertia but is a modeling choice that locally changes its spanwise distribution near the root.

### 6.3 Stiffness

For every element,

\[
K_e=\int_{x_e}B(x)^T C(x)B(x)\,dx.
\]

`C` contains `EA`, `EIz`, `EIy`, the bending cross-coupling `EIzy`, and `GJ`. An element crossing a tabulated discontinuity is split at that station before Gauss integration. Element matrices are scattered into `Ks_W` using the standard two-node 12-DOF index window.

### 6.4 Propeller/pylon coordinates and coupling blocks

Each propeller adds two structural coordinates:

```text
[pitch_about_span, yaw_about_vertical]
```

The modal inertias are obtained from the configured stiffness/frequency pairs, for example `In_theta=K_theta/omega_theta^2`. The rotor axial inertia similarly comes from the twist pair. Pylon/rotor mass, first moments, and cross inertias form the coupling blocks.

The full matrices are assembled as

\[
M_g=\begin{bmatrix}M_W+G_W&F_W\\B_P&M_P\end{bmatrix},\quad
C_g=\begin{bmatrix}C_W&H_W\\D_P&C_P\end{bmatrix},\quad
K_g=\begin{bmatrix}K_W&0\\0&K_P\end{bmatrix}.
\]

`B_P` and `F_W` carry wing--pylon inertial coupling. `H_W`, `D_P`, the local propeller skew block, and the local wing block carry terms proportional to `Ix_prop*Omega`; these are gyroscopic and intentionally skew-symmetric where appropriate.

The attachment operator interpolates all six wing quantities between the two nodes bracketing the exact propeller station. Its transpose is used in the force direction, preserving discrete virtual work. The first six global DOFs are removed to clamp the root.

### 6.5 Damping

The propeller modal block uses viscous diagonal damping `2*zeta*I*omega` plus its gyroscopic skew block. The complete reduced model then adds stiffness-proportional Rayleigh damping

\[
C_R=\frac{2\xi_{\rm ref}}{\omega_{\rm ref}}K.
\]

No nonlinear or aerodynamic damping matrix is assembled. Aerodynamic damping emerges from the unsteady coupled forces.

## 7. Aerodynamic and propeller model

### 7.1 Initialization

`initialize_bohnisch_uvlm_system` creates one wing surface and one surface per blade. Surface order is

```text
1: wing
2 ... 1+Nb: propeller 1 blades
then propeller 2 blades, etc.
```

The wing receives optional image symmetry across `y=0`; blades are never mirrored. The modeled structural wing remains one root-to-tip cantilever. When symmetry is true, the mirror influences the aerodynamics but does not create a second structural beam.

Interaction IDs differ from surface IDs:

- every surface has a unique `surface_id`, used in core/self-segment logic;
- the wing has interaction group 1;
- all blades of propeller `ip` share group `ip+1`.

With `interaction_on=false`, wing--propeller and different-propeller AIC, wake induction, near-field induction and wake-convection contributions are omitted. Blades belonging to the same propeller still interact. Structural attachment coupling remains active.

Bound circulation starts at zero and every active wake-row count starts at zero. Wake arrays are allocated to maximum capacity but inactive entries are not read.

### 7.2 Panel and circulation equations

Each physical quadrilateral is represented by a vortex ring whose chordwise vertices lie at quarter-chord positions. Its no-penetration control point is at three-quarter chord. `grid_to_surface_panels` computes the normal and finite-core radius.

At every trial, `get_surface_velocities!` evaluates

\[
V_{\rm surface,relative}=\frac{r_n-r_{n+1}}{\Delta t}=-V_{\rm mesh}.
\]

This sign is correct because the boundary condition is written using fluid velocity relative to the moving surface. `normal_velocity!` forms

\[
w_i=-n_i\cdot\left(V_\infty+V_{\rm body\ rotation}
 +V_{\rm additional}-V_{\rm mesh}+V_{\rm wake}\right).
\]

`influence_coefficients!` fills the bound-vortex matrix from regularized Biot--Savart ring velocities, including enabled symmetry images and interaction masks. The solve is

\[
AIC\,\Gamma=w,
\]

implemented by `circulation!(Gamma,AIC,w)`. Immediately before the solve, `dGamma/dt` is seeded with `-Gamma_n`; after it,

\[
\dot\Gamma_{n+1}\approx\frac{\Gamma_{n+1}-\Gamma_n}{\Delta t}.
\]

Wing and blades use the same panel, influence and circulation implementation. Their differences are geometry, symmetry, interaction group, motion, and discretization.

### 7.3 Aerodynamic force equations

The active `:imperial` model converts the package ring-circulation sign to the Imperial convention and constructs spanwise/chordwise circulation jumps. For a vortex segment from `r1` to `r2`,

\[
f_{\rm segment}=\rho\,\Delta\Gamma\,V\times(r_2-r_1).
\]

The velocity includes freestream, optional body rotation/additional flow, segment motion, and induced velocity from allowed bound surfaces and active wakes. Adjacent self-segments are excluded from the bound induction at their own force evaluation points.

The unsteady panel term is

\[
f_{\rm unsteady}=-\rho A n\dot\Gamma.
\]

`imperial_nodal_forces` distributes half of every segment force to each endpoint and one quarter of every unsteady panel force to each corner. `imperial_nodal_positions` returns the matching vortex-lattice vertices. The resultant is conserved, and load-transfer moments use the same points at which these forces are defined.

## 8. Structural deformation to aerodynamic geometry

For every trial state:

1. `update_aero_geometry_for_state!` splits the reduced vector into free wing DOFs and propeller pitch/yaw DOFs.
2. It reinserts six zeros for the clamped root.
3. It transforms structural components to the aerodynamic basis.
4. For each span station it places the elastic-axis point at its undeformed position plus translation.
5. It rotates the complete local chord line by `Rz*Rx*Ry` about that point.
6. It interpolates attachment translations and rotations with the same weights used in structural assembly.
7. It evaluates pitch, yaw, and prescribed spin using absolute `t[n+1]`.
8. It transforms every reference blade-grid point to its current global position.
9. It regenerates new `SurfacePanel` matrices from the physical grids and replaces entries of `system.surfaces`.

The mesh is regenerated from the undeformed reference grids and the current trial state. It is not incrementally deformed from the previous trial, preventing accumulation of geometric drift during fixed-point iterations.

All six wing DOFs are used. Because `ns_wing == Ne`, the structural nodes and aerodynamic span stations coincide in the active model. If those discretizations were separated in the future, `generate_panel_grid_and_interpolate` would need real spanwise interpolation; its present direct indexing assumes equal counts.

## 9. Aerodynamic force and moment transfer

### 9.1 Wing

1. `imperial_nodal_forces(system)[1]` gives force at every wing vortex vertex.
2. The deformed elastic-axis position is constructed at every structural/aerodynamic span node.
3. Every vertex moment is computed explicitly as

   \[
   M_A=(r_{\rm vortex}-r_{EA})\times F_A.
   \]

4. Forces and moments are summed across chord at each span node.
5. Force components are permuted/signed into `[span,chord,down]`.
6. Moments are projected by `chang_wing_generalized_moment` onto the current Euler-coordinate virtual-work axes.
7. Root entries are removed with the same `free_dofs` convention as `M`, `C`, and `K`.

### 9.2 Propeller

For every propeller, forces from all of its blade vertices are summed. The code first forms

\[
F_P=\sum_a F_a,\qquad
M_{hub}=\sum_a(r_a-r_{hub})\times F_a.
\]

The physical hub and the elastic-axis pivot are distinct. The modal coordinates
receive

\[
M_{modal}=M_{hub}+(r_{load}-r_{pivot})\times F_P,
\]

while the wing attachment receives the hub wrench translated to the pivot,

\[
M_{pivot}=M_{hub}+(r_{hub}-r_{pivot})\times F_P.
\]

The modal moment is projected onto the instantaneous pitch/yaw axes (or fixed
aerodynamic axes if requested), while the complete pivot wrench is projected
into structural coordinates and applied only to the selected attachment node.

## 10. Initial conditions, aerodynamic startup, trim, and perturbation

The structural histories start with

```text
U[1]   = 0
Ud[1]  = 0
Udd[1] = M \ (F_initial - C*Ud[1] - K*U[1])
```

The current call supplies no initial aerodynamic or external load, so `Udd[1]=0`. There is no prior steady aerodynamic equilibrium solve.

The UVLM starts from undeformed geometry, zero bound circulation and an empty wake. During the configured trim period, the complete unsteady aerodynamic calculation still runs and its wake develops, but `load_for_state` returns zero to the structure. The structure is consequently held at its undeformed perturbation reference while full aerodynamic generalized loads are accumulated.

Over the final configured number of revolutions before the pulse, accepted full-load samples are averaged. At the impulse start, that mean becomes `trim_load`. Later structural forcing is

\[
F_{aero,perturbation}=F_{aero,full}-F_{trim}.
\]

This is a rigid-reference mean-load subtraction, not a static aeroelastic trim. The structure is not first deflected to an equilibrium under the mean aerodynamic load.

The selected propeller pitch DOFs receive a finite-duration Hann moment. Both `F_external,n` and `F_external,n+1` are calculated because generalized-alpha weights the two; Newmark uses the `n+1` value.

## 11. Newmark-beta as implemented

For the configured `alpha`, `newmark_beta_parameters` sets

\[
\gamma=0.5+\alpha,\qquad
\beta=\frac14(\gamma+0.5)^2.
\]

For the default `alpha=0.05`, `gamma=0.55` and `beta=0.275625`.

`newmark_beta_coefficients` forms the standard `a0...a7`. In `newmark_beta_corrector`:

\[
K_{eff}=K+a_0M+a_1C
\]

maps to `effective_stiffness = K .+ a0.*M .+ a1.*C`, and

\[
F_{eff}=F_{n+1}+M(a_0U_n+a_2\dot U_n+a_3\ddot U_n)
+C(a_1U_n+a_4\dot U_n+a_5\ddot U_n)
\]

maps to `total_force`, followed by `effective_force`. The solve and recovery are exactly

\[
U_{n+1}=K_{eff}^{-1}F_{eff},
\]

\[
\ddot U_{n+1}=a_0(U_{n+1}-U_n)-a_2\dot U_n-a_3\ddot U_n,
\]

\[
\dot U_{n+1}=\dot U_n+a_6\ddot U_n+a_7\ddot U_{n+1}.
\]

In an implicit aeroelastic step, `structural_effective_system` factorizes `Keff` once and all coupling corrections reuse it. `M`, `C`, `K`, `U_n`, `Udot_n`, and `Uddot_n` remain unchanged during those iterations.

`newmark_beta_equilibrium_residual` separately checks the physical end-state equation with the candidate aerodynamic load; this is distinct from the linear-solver residual of `Keff*U-F_eff`.

## 12. Generalized-alpha as implemented

`generalized_alpha_parameters(rho_inf)` uses the Chung--Hulbert second-order parameters

\[
\alpha_m=\frac{2\rho_\infty-1}{\rho_\infty+1},\quad
\alpha_f=\frac{\rho_\infty}{\rho_\infty+1},
\]

\[
\gamma=\frac12-\alpha_m+\alpha_f,\quad
\beta=\frac14(1-\alpha_m+\alpha_f)^2.
\]

The equilibrium enforced by `generalized_alpha_corrector` is

\[
M a_{n+1-\alpha_m}+C v_{n+1-\alpha_f}+K u_{n+1-\alpha_f}
=F_{n+1-\alpha_f},
\]

where the code defines, for example,

\[
a_{n+1-\alpha_m}=(1-\alpha_m)a_{n+1}+\alpha_m a_n
\]

and analogously for state and force with `alpha_f`. `force_alpha` includes both aerodynamic and external loads at `n` and `n+1`.

The effective matrix is

\[
K_{eff}=(1-\alpha_m)a_0M
+(1-\alpha_f)\frac{\gamma}{\beta\Delta t}C
+(1-\alpha_f)K.
\]

The corrector constructs the remaining known kinematic terms on the right-hand side, solves for `U[n+1]`, then uses `generalized_alpha_kinematics` to recover velocity and acceleration. As for Newmark, the implicit loop reuses one factorization and checks both the linear solve and full coupled equilibrium.

At `rho_inf=1`, the parameters reduce to `alpha_m=alpha_f=0.5`, `gamma=0.5`, `beta=0.25`; the weighted equation has no intended high-frequency algorithmic damping. Smaller `rho_inf` introduces controlled high-frequency dissipation.

## 13. Exact physical step: loose explicit coupling

Let the accepted state at the start be

```text
structural: U_n, Ud_n, Udd_n
aerodynamic: surfaces_n, Gamma_n, wake_n, active_rows_n
previous accepted aero perturbation load: F_aero,n
```

The actual sequence is:

1. `copy_surfaces_to_previous!` copies current accepted panels into `previous_surfaces`. These represent geometry at `t_n`.
2. `snapshot_uvlm` deep-copies wake, `Gamma_n`, circulation derivatives, current/previous surfaces, shedding locations, active rows and freestream.
3. External loads are evaluated at `t_n` and `t[n+1]`.
4. `load_for_state(U_n)` is called exactly once.
5. The snapshot is restored, so the trial begins from the accepted aerodynamic state at `n`.
6. Geometry is regenerated at absolute time `t[n+1]`, but using lagged structural displacement `U_n`. Rotor azimuth is therefore the target azimuth `-Omega*t[n+1]`; the elastic deformation is one exchange behind.
7. `propagate_system!(advance_wake=false)` computes relative mesh velocity from accepted geometry at `n` to this target-time/lagged-structure geometry, updates the shedding location, solves `Gamma_{n+1}^{explicit}`, computes `dGamma/dt`, and computes full loads. Wake panels are not convected or shed.
8. Before trim capture, the structural aerodynamic load is zero. After capture it is `F_full-F_trim`.
9. `loose_explicit_aeroelastic_step` performs one Newmark or generalized-alpha solve, producing `U[n+1]`, `Ud[n+1]`, `Udd[n+1]`.
10. The structural state is accepted. The aerodynamic system still corresponds to the last trial, whose deformation was `U_n`, not the just-computed `U[n+1]`.
11. `commit_chang_wake!` convects the active wake and sheds one row using the accepted trial circulation and geometry, then increments active rows.
12. Histories/diagnostics/optional animation are stored.

Therefore the loose scheme evaluates an aerodynamic operator labeled at target physical time/azimuth `n+1`, but with lagged structural displacement `U_n`. This is the scheme's explicit one-step aeroelastic lag. There is not an additional delayed reuse of an older aerodynamic load: a fresh UVLM solve is performed in each physical step.

## 14. Exact physical step: implicit predictor--corrector

Steps 1--3 above are identical. Then:

1. `structural_effective_system` builds/factorizes the selected method's constant matrix.
2. `structural_predictor` creates

   \[
   U_{n+1}^{(0)}=U_n+\Delta t\dot U_n+Delta t^2(1/2-\beta)\ddot U_n.
   \]

3. At iteration `k`, the current `state_guess=U[n+1]^(k)` is saved as the aerodynamic evaluation state.
4. `aero_load_for_state!` restores the **same** beginning-of-step snapshot.
5. It regenerates geometry at `t[n+1]` from `U[n+1]^(k)` and absolute rotor azimuth `-Omega*t[n+1]`.
6. It solves one circulation/load trial from `Gamma_n` and `wake_n`, with wake advancement disabled, returning `F_aero,n+1^(k)`.
7. The structural corrector solves for `U_corrected^(k+1)` using this load and the fixed previous physical state.
8. It evaluates:
   - scaled state residual between corrected displacement and the state passed to the aero callback;
   - scaled load residual between successive candidate loads (after trim exists);
   - effective linear-system residual;
   - complete Newmark or generalized-alpha coupled-equilibrium residual at the aerodynamic evaluation state.
9. If all required tolerances pass, it accepts the exact state that produced the current aerodynamic trial. Velocity and acceleration are recovered from that accepted displacement by the integrator's kinematic equations.
10. Otherwise only displacement is relaxed:

    \[
    U^{(k+1)}=\omega U_{corrected}^{(k+1)}+(1-\omega)U^{(k)}.
    \]

    Velocity and acceleration are not independently relaxed.
11. A failure after the maximum iterations restores the snapshot and rejects the physical step.
12. A successful step commits the wake once, outside the iteration.

Verified invariants during subiterations:

- `displacement[step]`, `velocity[step]`, and `acceleration[step]` are read-only.
- physical time is an argument; no time counter is incremented;
- rotor azimuth is recomputed from the same `t[n+1]` on every trial;
- every trial begins from identical `Gamma_n`, wake geometry, active rows and previous/current surfaces;
- `advance_wake=false` prevents convection/shedding during trials;
- the accepted wake is advanced exactly once after convergence.

## 15. Free-wake algorithm

### 15.1 Shedding location and circulation

During the circulation trial, `update_wake_shedding_locations!` puts the downstream edge of the new near-wake panel at

\[
r_{shed}=r_{TE}+\eta V_{TE,relative}\Delta t,
\]

with `eta=0.1` in the Chang call. This velocity includes freestream, configured body rotation, optional additional flow, and trailing-edge mesh-relative velocity. It does not include induced velocity for this initial separation.

On commit, `shed_wake!` creates a panel from these stored locations and assigns it the circulation of the corresponding trailing-edge bound panel. This is the code's Kelvin-like shedding rule; there is no separate global Kelvin-constraint equation.

### 15.2 Convection

`get_wake_velocities!` treats the newest shedding-location vertices with prescribed freestream/body/additional/trailing-edge-motion velocity. Older active wake vertices additionally receive induced velocity from allowed bound surfaces and active wakes, including self-induced wake velocity with finite-core regularization and enabled symmetry images.

`translate_wake!` moves active panels by explicit first-order position update `r_new=r_old+V*dt`. `rowshift!` rotates the fixed storage so the newly shed row becomes row 1.

### 15.3 Truncation and interaction

One row is activated per accepted physical step until the configured capacity is reached. Once full, the last/oldest row is overwritten before the array is shifted. Thus truncation is by fixed wake age/row count; there is no wake relaxation or far-wake roll-up model beyond explicit convection and finite-core regularization.

The wing and propeller wake capacities can differ. Interaction masks are applied consistently to AIC assembly, wake contribution to no-penetration, near-field induced velocity, and wake convection.

## 16. Persistent state and mutation

| State | Meaning | Updated by | Time level/persistence |
|---|---|---|---|
| `displacement`, `velocity`, `acceleration` | Accepted structural histories | `solve_chang_aeroelastic!` after convergence | Persistent `n`, `n+1` vectors |
| `accepted_perturbation_load` | Previous accepted aerodynamic perturbation load | Main time loop | Persistent; needed at `n` by generalized-alpha |
| `trim_load` | Mean full aerodynamic generalized load | Main time loop at trim capture | Persistent baseline |
| `system.surfaces` | Current bound-panel geometry | `update_aero_geometry_for_state!` | Trial during coupling; final trial persists after acceptance |
| `system.previous_surfaces` | Accepted geometry at start of step | `copy_surfaces_to_previous!`, restored by snapshot | Represents `n` during trials |
| `system.Γ` | Current bound circulation | `circulation!` inside propagation | Trial; final trial accepted |
| `system.dΓdt` | First-order bound-circulation rate | `propagate_system!` | Trial; final trial accepted |
| `system.wakes` | Wake panel vertices, core and circulation | `advance_wake!`/`shed_wake!` | Unchanged in trials, persistent after commit |
| `workspace.iwake` | Authoritative active rows used by the Chang solver | `commit_chang_wake!` | Persistent, capped |
| `system.nwake` | Backend copy of active rows | `propagate_system!`/`advance_wake!` | Synchronized before a solve/commit, but currently lags `iwake` immediately after a Chang commit; see the issue list |
| `system.wake_shedding_locations` | Target near-wake edge | aerodynamic trial | Trial; final value used by commit |
| `system.Vcp`, `Vh`, `Vv`, `Vte` | Relative surface-motion velocities | `get_surface_velocities!` | Scratch/trial; `Vte` reused at commit |
| `system.AIC`, `w`, `properties`, force buffers | UVLM scratch and current loads | propagation/near-field solve | Overwritten each trial |
| propeller azimuth | `-Ω*t_target` | computed in geometry mapping | Derived, not stored/incremented |
| hub/load points and pitch/yaw axes | Current trial propeller kinematics | geometry mapping | Workspace scratch used immediately by load mapping |

Important mutating functions:

- `copy_surfaces_to_previous!`: copies panel values into preallocated previous-surface matrices.
- `restore_uvlm!`: overwrites wake, circulation, surfaces, shedding locations, rows and freestream in place while retaining array identities.
- `update_aero_geometry_for_state!`: replaces current panel matrices and overwrites propeller kinematic buffers.
- `propagate_system!`: overwrites mesh velocities, AIC, RHS, circulation, circulation rate, panel properties and force buffers; optionally the wake.
- `near_field_forces!`: returns new force matrices and overwrites panel properties.
- `assemble_structural_aero_load!`: overwrites reusable nodal force/moment arrays.
- `advance_wake!`, `translate_wake!`, `shed_wake!`, `rowshift!`: permanently mutate wake storage when called at commit.

## 17. Output and postprocessing

The returned `solution` retains the complete accepted `U`, `Udot`, and `Uddot` histories plus trim load and coupling diagnostics. The workspace retains only the final aerodynamic system unless optional frames were recorded.

The history CSV writes:

- time;
- wing-tip vertical-down displacement;
- wing-tip span-axis twist;
- each propeller pitch/yaw angle;
- iteration count, state/load/coupled-equilibrium residual, and convergence flag.

The summary records model choices, completion/finite/convergence status and response maxima. The default CSV does **not** store velocity, acceleration, full aerodynamic generalized-load history, circulation history, or wake history. Those would require extending result storage, not merely postprocessing the existing CSV.

## 18. Exact complete algorithm

### 18.1 Initialization common to both schemes

1. Load editable defaults from `chang_case.jl`.
2. Apply typed environment overrides and resolve automatic settings.
3. Validate the complete case.
4. Convert Chang reference tables to SI and create the active meshes/distributions.
5. Compute constant-advance-ratio `Omega`, `dt`, time grid, freestream and reference data.
6. Remap spatial inertias and integrate distributed beam stiffness.
7. Assemble full wing/propeller `M_global`, `C_global`, `K_global`.
8. Remove six root DOFs, add Rayleigh damping, and validate the reduced matrices.
9. Create undeformed wing and blade grids.
10. Convert them to vortex-ring panels; allocate UVLM arrays and fixed wake capacity.
11. Set symmetry, interaction groups, zero circulation and zero active wake rows.
12. Allocate load-transfer and optional animation storage.
13. Select integrator parameters and coupling scheme independently.
14. Set `U0=Udot0=0`; calculate `Uddot0` from structural equilibrium.
15. Enter the physical time loop.

### 18.2 Loose explicit step

1. Copy accepted geometry to `previous_surfaces` and snapshot the accepted UVLM state.
2. Compute external pulse at `n` and `n+1`.
3. Restore the snapshot and build target-time geometry from `U_n`.
4. Solve target-time bound circulation and loads from the accepted wake/`Gamma_n`, without wake convection.
5. Map all wing/blade loads into one generalized vector; subtract trim if active.
6. Perform one selected-integrator structural solve.
7. Check linear and physical equilibrium residuals.
8. Accept structural state and current aerodynamic trial.
9. Update trim accumulation/capture if applicable.
10. Convect active wake panels and shed/activate one row exactly once.
11. Store diagnostics, optional accepted frame, and safety-check the state.

### 18.3 Implicit predictor--corrector step

1. Copy accepted geometry and snapshot as above.
2. Compute external loads and factor the effective structural matrix once.
3. Predict `U[n+1]^(0)`.
4. For each coupling iteration:
   1. restore the same snapshot;
   2. build geometry from `U[n+1]^(k)` at fixed `t[n+1]`/azimuth;
   3. solve circulation and loads with no wake advance;
   4. map/subtract trim;
   5. solve the structural corrector;
   6. compute state, load, linear and coupled-equilibrium residuals;
   7. accept if all required checks pass, otherwise relax displacement.
5. If no convergence, restore the snapshot and abort without committing the step.
6. If converged, store the state that generated the accepted aerodynamic trial.
7. Update trim logic.
8. Convect/shed/activate the wake exactly once.
9. Store diagnostics/frame and safety-check.

## 19. Implementation verification summary

### Confirmed correct

- The Newmark effective matrix, effective force, velocity and acceleration recovery match the stated Newmark equations term by term.
- Generalized-alpha parameters, weighted state/load levels, effective matrix and recovery equations are internally consistent.
- Integrator selection and coupling selection are independent.
- The implicit loop holds previous physical structural states fixed, restores the same aerodynamic snapshot for every trial, uses absolute target rotor time, and commits the wake once.
- Loose coupling performs one fresh target-time UVLM solve per physical step; its lag is specifically structural `U_n` in target-time geometry.
- Initial acceleration is computed from equilibrium rather than assumed independently.
- The structural and aerodynamic basis conversion is applied consistently to translations and forces.
- Wing and propeller moment arms use positions paired with the dimensional vortex-vertex forces.
- Propeller attachment kinematics and load transfer use the same single node; the dedicated structural verification reports zero attachment virtual-work error.
- The complete package suite passes 388 assertions across its test summaries.
- The four integrator/coupling smoke combinations and their refined Newmark cases pass 31/31 checks.
- The structural verification passes with positive reduced `M`/`K`, symmetry errors near machine precision, conservative total spatial inertia, and first reference modes within the reported comparison errors.
- The finite-difference load-transfer audit covers reference, uniform-rotation, and nonuniform-rotation states; its maximum scaled virtual-work discrepancy was approximately `9.25e-8`.

### Potential concerns or modeling assumptions

1. **Linear structure with finite-rotation aerodynamic geometry.** This is a mixed model, not a geometrically nonlinear beam. Establish a response-amplitude range where it remains valid.
2. **Rigid-reference aerodynamic trim.** The mean load is subtracted while the structure stays at zero; no static aeroelastic equilibrium is found. This should match the intended perturbation experiment.
3. **Root inertia relocation.** Moving the clamped-root control-volume block to the first free node preserves totals but changes local inertia distribution. Retain the current modal comparison as a regression check.
4. **Direct-node station resolution.** A requested continuous span fraction is rounded to a structural node in the Chang compatibility input. Refine or place the structural nodes explicitly when the exact installation station matters.
5. **Trim mean endpoint.** Sampling both ends of a nominal periodic window can duplicate a phase when the interval is exactly an integer revolution. The effect is small for many samples but can be removed with a half-open window or phase-aware average.
6. **Newmark input range.** The code checks only finite `alpha` and positive `gamma,beta`. The stated unconditional-stability family also assumes `alpha>=0`; negative values currently pass if the derived coefficients remain positive.
7. **`system.grids` is not kept current by the Chang adapter.** The active solve uses `system.surfaces`, so this does not affect propagation, but generic inspection of `system.grids` can see stale/uninitialized data.
8. **Wake circulation stretching rule.** `translate_wake` replaces panel `gamma` by `gamma*l_old/l_new`. If `gamma` is circulation, Kelvin conservation would normally retain it rather than scale by filament perimeter. The code behavior is verified, but the intended vorticity representation must be confirmed before calling it physically correct.
9. **First-order wake convection and finite truncation.** These are valid numerical choices but require time-step/core/wake-length convergence, especially for flutter prediction.
10. **Current history output is structural.** It is insufficient by itself to audit circulation, aerodynamic energy/work, or wake convergence over time.

### Definite implementation/robustness issues found

1. **Zero wake capacity is accepted but cannot be committed.** `src/chang_configuration.jl:238-239` accepts `maximum_rows_* == 0`, while `../../src/backend/wake.jl:798` writes `wake[end,j]`. A zero-row wake therefore fails on the first accepted step. Require strictly positive capacity, or explicitly bypass shedding for a zero-capacity surface.
2. **Repeated trailing-edge velocity reuse does not skip the outer receiver loop.** In `get_wake_velocities!` at `../../src/backend/wake.jl:321-333` and `:361-373`, `continue` is inside the loop over duplicate matches, so the copied velocity is immediately recomputed afterward. This is a control-flow bug for configurations with coincident trailing-edge vertices. Use a flag/helper and continue the enclosing receiver-vertex loop after a match.
3. **Runs shorter than one computed time step are not rejected.** `src/chang_model_parameters.jl:436` can create a one-state, zero-step range; `src/chang_postprocessing.jl:101-103` and `:201` subsequently assume at least one coupling diagnostic. Configuration validation should require `end_time >= dt` or postprocessing should explicitly support a zero-step run.
4. **The returned backend active-row field can lag the Chang counter.** `commit_chang_wake!` passes `wake.active_rows` into `advance_wake!`, then increments only `wake.active_rows` at `src/chang_simulation.jl:280-284`. Consequently `system.nwake` is one activation behind immediately after a commit and in a returned short-run workspace. The next propagation resynchronizes it, so the active time marcher uses the correct external count, but generic inspection/default backend calls can see stale state. Assign `system.nwake .= wake.active_rows` after incrementing.

The context-validation script had two stale hard-coded expectations after the case defaults changed (wake rows and default coupling scheme). They were changed to derive from the active defaults; this was a test-maintenance issue, not a solver defect.

### Verification executed for this review

| Check | Result | What it establishes |
|---|---:|---|
| `Pkg.test()` | 167/167 pass | Package integrators, UVLM transactions, near-field forces, wake history, grids and structural helpers remain regression-clean. |
| `verify_chang_context.jl` | 38/38 pass | Per-case configuration ownership, environment parsing, independent workspaces and deterministic interleaved trials. |
| `verify_chang_solver_combinations.jl` | 31/31 pass | Both integrators execute with both coupling schemes; implicit cases converge in the smoke problem. This is execution/regression evidence, not a flutter-convergence proof. |
| `verify_chang_structural_model.jl` | pass | Positive/symmetric reduced matrices, conservative mass remap, modal comparison and attachment virtual work. |
| `audit_chang_virtual_work.jl` | pass; max scaled error about `4.99e-8` | Finite-difference agreement of mapped generalized loads with aerodynamic virtual work for multiple nonzero states. |

### Recommended numerical verification tests

1. Single-DOF undamped/damped oscillator against the analytical solution for both integrators.
2. Newmark energy decay test for several `alpha` values and explicit rejection of unsupported negative `alpha`.
3. Generalized-alpha spectral-radius and high-frequency dissipation regression.
4. Structural free-vibration modal time histories against the generalized eigenproblem.
5. Static unit-load/virtual-work tests for every wing DOF and propeller pitch/yaw coordinate.
6. Finite-difference derivative of aerodynamic virtual work versus every geometry coordinate, including nonzero Euler angles.
7. Fixed-wing steady/unsteady UVLM benchmarks with wake disabled/enabled.
8. Isolated propeller thrust, torque, rotation-sign and azimuth-periodicity tests.
9. Wake convection test with a uniform prescribed velocity: exact expected vertex displacement and unchanged circulation.
10. Coincident trailing-edge test that exercises `repeated_points`.
11. Loose versus implicit comparison under decreasing `dt`; both should approach a common response.
12. Implicit tolerance/relaxation refinement and failed-step rollback test.
13. Independent refinement of wing panels, blade panels, azimuth step, core radius and retained wake age.
14. Aerodynamic-work/structural-energy balance around the force/moment mapping.
15. Flutter-speed convergence with time step, coupling tolerance, wake length and finite core varied separately.

## 20. Bottom line

The active architecture is a partitioned linear-structural/free-wake-UVLM solver with transactional aerodynamic trials. The structural solver sends a displacement guess to a case-specific geometry adapter; the shared UVLM backend returns dimensional vortex-vertex forces; the adapter maps those forces and moments into exactly the reduced structural coordinate ordering. Newmark-beta and generalized-alpha act only on the structural equation. Loose and implicit coupling determine how often that same aerodynamic load operator is evaluated before an accepted physical state. The wake is deliberately excluded from coupling subiterations and advanced once after acceptance.

The central time-level and mutation logic is coherent and is covered by passing package tests. The three definite robustness issues above and the wake-circulation stretching interpretation should be addressed or explicitly accepted before treating a flutter result as fully verified.
