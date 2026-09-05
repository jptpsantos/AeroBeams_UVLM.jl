# WingPropellerUVLM library guide

## 1. What the library is

`WingPropellerUVLM` is the in-repository Julia package that owns the global
wing--propeller unsteady vortex-lattice model. It provides:

- wing and propeller lattice generation;
- bound-circulation and free-wake evolution;
- optional aerodynamic interaction between lifting-surface groups;
- Imperial College C++ UVLM-compatible near-field forces;
- nodal aerodynamic forces and positions for structural coupling;
- rollback-safe aerodynamic trial steps; and
- a reusable partitioned generalized-alpha integrator for the current linear
  Chang validation model.

The package is located at `lib/WingPropellerUVLM`. It is deliberately separate
from the `AeroBeams` module because a UVLM is one global aerodynamic system:
the wing, all blades, their bound circulation, and all wakes influence the same
dense solve. It should not be split into AeroBeams' element-local aerodynamic
states.

The present status is important:

- the rigid free-wake UVLM example is runnable;
- the linear Chang structural model is coupled and runnable in the time domain;
- the force calculation uses the Imperial near-field formulation;
- the transactional APIs needed for nonlinear coupling exist; but
- the final adapter to AeroBeams' geometrically exact nonlinear dynamic solver
  is still an implementation step. The architecture for that step is described
  in [Coupling AeroBeams with a wing--propeller UVLM](wing-propeller-uvlm-coupling.md).

## 2. Run the package locally

From the repository root, instantiate and test the local package once:

```powershell
julia --project=lib/WingPropellerUVLM -e "using Pkg; Pkg.instantiate(); Pkg.test()"
```

Run the rigid rectangular-wing case:

```powershell
julia lib/WingPropellerUVLM/examples/rectangular_wing_free_wake.jl
```

Run the coupled Chang case:

```powershell
julia lib/WingPropellerUVLM/examples/chang_linear_aeroelastic/run_chang_linear_aeroelastic.jl
```

Both input files call `Pkg.activate` with a path based on `@__DIR__`, so VS
Code's **Run Julia File** command also uses the correct local project. Do not
use `Pkg.add("WingPropellerUVLM")`: this package is not obtained from the Julia
registry.

To use it from a separate Julia project, develop the local path once:

```julia-repl
pkg> develop path="path/to/AeroBeams_UVLM.jl/lib/WingPropellerUVLM"
```

Then use `using WingPropellerUVLM` normally in that project.

## 3. The mental model

There are four distinct data layers:

| Layer | Main object | Meaning |
|:--|:--|:--|
| Geometry | grids and `SurfacePanel` matrices | Current wing and blade vortex lattices |
| Aerodynamic state | `System` | Circulation, panel motion, loads, wakes, and work arrays |
| Transaction state | `UVLMSnapshot` | Beginning-of-step copy used to repeat or reject a trial |
| Structural state | `solution.*_history` in the Chang case | Free generalized displacement, velocity, and acceleration |

The `System` is mutable even though it is declared as an immutable Julia
`struct`: its fields are arrays, and the contents of those arrays change in
place. One `System` contains every aerodynamic surface. In the Chang case the
surface order is:

```text
surface 1       wing
next Nb         propeller 1 blades
next Nb         propeller 2 blades, if present
...
```

`system.nwake[i]` is the active wake-row count for surface `i`. The wake matrix
is allocated to its maximum size at initialization; only the first active rows
participate in a step.

## 4. Coordinate and DOF conventions in the Chang example

The UVLM aerodynamic frame is

```text
x_A = chord/freestream direction
y_A = wing span direction
z_A = upward
```

The Chang structural nodal order is

```text
[u_span, v_chord, w_down, theta_span, theta_chord, theta_vertical]
```

Therefore `chang_uvlm_coupling.jl` applies

```text
u_x_A     = v_chord
u_y_A     = u_span
u_z_A     = -w_down
theta_x_A = theta_chord
theta_y_A = theta_span
theta_z_A = -theta_vertical
```

The reverse mapping is applied to forces and moments:

```text
F_struct = [F_y_A, F_x_A, -F_z_A]
M_struct = [M_y_A, M_x_A, -M_z_A]
```

These swaps and signs form a proper right-handed rotation. They are not merely
plotting conventions; an incorrect sign changes aerodynamic work and can
create artificial damping or excitation.

Wing free DOFs exclude the clamped root. The global reduced state is

```text
[all free wing DOFs,
 propeller 1 pitch, propeller 1 yaw,
 propeller 2 pitch, propeller 2 yaw, ...]
```

## 5. One accepted UVLM time step

For a partitioned moving-geometry problem, separate repeatable aerodynamic
trials from the one accepted wake update:

```julia
copy_surfaces_to_previous!(system, length(system.surfaces))
snapshot = snapshot_uvlm(system)

# Repeat this block for each structural coupling guess.
restore_uvlm!(system, snapshot)
update_geometry_for_trial!(system, structural_guess)
propagate_system!(system, freestream, dt;
    additional_velocity = nothing,
    repeated_points = repeated_trailing_edge_points(system.surfaces),
    nwake = active_wake_rows,
    eta = 0.1,
    calculate_influence_matrix = true,
    near_field_analysis = true,
    derivatives = false,
    interaction_id = system.surface_id,
    interaction = true,
    advance_wake = false,
)

# Call only after the final trial passes every coupled convergence test.
advance_wake!(system, freestream, dt;
    additional_velocity = nothing,
    repeated_points = repeated_trailing_edge_points(system.surfaces),
    nwake = active_wake_rows,
    interaction_id = system.surface_id,
    interaction = true,
)

commit_wake_rows!(active_wake_rows, maximum_wake_rows)
system.nwake .= active_wake_rows
```

`propagate_system!` still advances the wake by default, so rigid and prescribed
one-pass analyses retain their original behavior. The `advance_wake=false`
option is specifically for repeated trials representing the same physical time
interval. `advance_wake!` uses the accepted circulation, trailing-edge motion,
and shedding locations already stored in `system`; it does not re-solve the
AIC system or recalculate loads.

Do not call `propagate_system!` again after accepting the final coupling trial.
That would solve circulation and loads twice. Call only `advance_wake!`, then
increment the active-row counter.

## 6. What `propagate_system!` updates

`propagate_system!` performs the following update in order:

1. Compute panel control-point, bound-segment, and trailing-edge velocities
   from `current_surfaces - previous_surfaces` divided by `dt`.
2. Store the current freestream and active wake-row counts.
3. Update the locations from which the new wake row will be shed.
4. Rebuild the aerodynamic influence coefficient matrix when the geometry is
   moving or `calculate_influence_matrix=true`.
5. Update trailing-edge influence coefficients.
6. Form the no-penetration right-hand side from freestream, surface motion,
   bound surfaces, and active wakes.
7. Solve the dense linear system for the new bound circulation `Gamma`.
8. Form `dGamma/dt = (Gamma[n+1] - Gamma[n])/dt`.
9. Compute Imperial near-field segment and unsteady forces when
   `near_field_analysis=true`.
10. If `advance_wake=true` (the default), call `advance_wake!` to compute wake
    convection velocities, convect existing panels, and shed the new row.

With `advance_wake=false`, steps 1--9 are still fully evaluated and provide the
same circulation and Imperial loads; only step 10 is deferred.

The keyword `interaction` controls cross-group influence. With
`interaction=true`, all surfaces interact. With `interaction=false`, a
receiving surface sees only sending surfaces with the same entry in
`interaction_id`. The Bohnisch initializer assigns group 1 to the wing, group
2 to all blades of propeller 1, group 3 to all blades of propeller 2, and so
on. Thus an isolated propeller retains interaction among its own blades while
wing--propeller and propeller--propeller interactions are disabled.

For the Imperial force path, use `derivatives=false`; freestream derivatives of
that formulation have not been implemented.

## 7. How the Chang coupled time step works

The principal driver is
`lib/WingPropellerUVLM/examples/chang_linear_aeroelastic/run_chang_linear_aeroelastic.jl`.
It presents the setup/solve/output sequence described here and delegates the
step-level algorithm to `src/chang_simulation.jl`.

### 7.1 Before the loop

The driver:

1. activates the local package environment;
2. reads the case and model parameters;
3. assembles root-constrained structural `M`, `C`, and `K` matrices;
4. asks `solve_chang_aeroelastic!` to initialize the structural histories;
5. creates the global UVLM `System`, reference grids, wakes, interaction groups,
   and load-transfer buffers;
6. defines a smooth pitch perturbation and a periodic trim-load averaging
   window; and
7. creates generalized-alpha and fixed-point coupling options.

`solution.displacement_history[it]`, `solution.velocity_history[it]`, and
`solution.acceleration_history[it]` are the accepted structural state at
`t[it]`. The solver also retains the accepted aerodynamic perturbation load at
that same time.

### 7.2 Beginning a physical step

At the top of step `it`, the driver executes:

```julia
copy_surfaces_to_previous!(system, wake.surface_count)
aerodynamic_snapshot = snapshot_uvlm(system)
```

The first line freezes the accepted geometry at `t[n]`, which is required to
calculate surface velocity. The second line records the complete aerodynamic
history from which every `t[n+1]` coupling trial must start.

### 7.3 The structural--aerodynamic callback

`partitioned_generalized_alpha_step` receives an anonymous callback. Its core
operation is:

```julia
full_aerodynamic_load .= aero_load_for_state!(
    system, aerodynamic_snapshot, state_guess, step
)
```

`aero_load_for_state!` is defined in `chang_uvlm_coupling.jl` and performs:

```text
restore_uvlm!(system, aerodynamic_snapshot)
        |
        v
update_aero_geometry_for_state!(system, state_guess, t[it+1])
        |
        v
propagate_system!(system, freestream[it], dt[it], ...; advance_wake=false)
        |
        v
assemble_structural_aero_load!(system, kinematics)
```

Restoring first makes circulation history and geometry identical for every
trial. Wake convection is disabled in the callback because all iterations
represent the same physical time interval.

### 7.4 Geometry update

`update_aero_geometry_for_state!`:

- reconstructs the clamped root in each wing DOF array;
- converts structural translations and rotations to the aerodynamic frame;
- interpolates the nodal motion to the wing lattice about the 30%-chord elastic
  axis;
- evaluates propeller pitch and yaw from the final structural DOFs;
- applies full wing-node rotation as `Rz * Rx * Ry`;
- applies propeller whirl as `Ry(pitch) * Rz(yaw)`;
- applies prescribed spin as `Rx(-Omega*t[n+1])`;
- tracks the pivot, physical hub, and modal load point separately; and
- replaces `system.surfaces` with panels made from the trial grids.

It changes only geometry. The subsequent `propagate_system!` call solves the
aerodynamics.

### 7.5 Imperial load transfer

After propagation, `imperial_nodal_forces(system)` returns a matrix of
dimensional body-frame forces for every surface. The matching
`imperial_nodal_positions(system)` returns the actual vortex vertices at which
those forces act.

For the wing, the adapter:

1. forms each vertex moment about the deformed elastic-axis node;
2. sums all chordwise forces and moments at a span station;
3. rotates them to the structural component order; and
4. removes the clamped-root entries.

For each propeller, it:

1. sums all blade vertex forces;
2. forms the aerodynamic moment about the physical hub;
3. shifts that wrench to the pylon modal load point for pitch/yaw generalized
   forces;
4. shifts it to the wing attachment pivot for the wing nodal load; and
5. appends `[pitch moment, yaw moment]` to the reduced load vector.

The returned vector has exactly the same ordering as the structural state.

### 7.6 Generalized-alpha fixed-point iteration

For a second-order system

```math
M\ddot{q} + C\dot{q} + Kq = F_\mathrm{aero} + F_\mathrm{external},
```

the parameters are computed from the requested high-frequency spectral radius
`rho_inf`:

```math
\alpha_m=\frac{2\rho_\infty-1}{\rho_\infty+1},\quad
\alpha_f=\frac{\rho_\infty}{\rho_\infty+1},\quad
\gamma=\frac12-\alpha_m+\alpha_f,\quad
\beta=\frac14(1-\alpha_m+\alpha_f)^2.
```

For each outer iteration, `partitioned_generalized_alpha_step`:

1. evaluates aerodynamic load at the current displacement guess;
2. solves one linear generalized-alpha structural corrector with that load;
3. compares the corrected displacement with the displacement at which the
   aerodynamic load was evaluated;
4. compares successive aerodynamic loads;
5. evaluates the complete coupled equilibrium residual;
6. accepts only if all enabled tolerances pass; otherwise
7. forms the next guess with fixed relaxation:

```math
q^{k+1}_\mathrm{guess}=omega q^{k}_\mathrm{corrected}
 +(1-\omega)q^{k}_\mathrm{guess}.
```

Translations, rotations, forces, and moments use separate physical scales in
the residuals. This avoids comparing unlike units in one raw vector norm.

The returned displacement is the last state evaluated by the aerodynamic
callback. Its velocity and acceleration are recomputed from that displacement
using the generalized-alpha/Newmark kinematic relations. Consequently the
returned structural state and the UVLM state left in `system` are a consistent
fixed-point pair.

If convergence fails, the solver restores `aerodynamic_snapshot` and terminates without
committing the step.

### 7.7 Trim baseline and perturbation load

The rotating wake produces periodic generalized loads at zero structural
perturbation. During startup the callback still runs the UVLM, but returns zero
load to the structure. The driver averages full aerodynamic loads over the
configured final trim revolutions:

```math
F_0 = \mathrm{mean}\left(F_\mathrm{aero}(q=0,t)\right).
```

After the baseline is captured, structural dynamics uses

```math
F_\mathrm{pert}(q,t)=F_\mathrm{aero}(q,t)-F_0.
```

This is a perturbation analysis about the mean periodic aerodynamic load; it is
not a static structural trim to a deformed equilibrium.

### 7.8 Accepting the step

On convergence, the driver stores the structural state and load, calls
`advance_wake!` exactly once with the accepted aerodynamic state, and then
increases `iwake` by one row per surface up to the allocated maximum. It does
not call `propagate_system!` again. At the beginning of the next loop, the
accepted `system.surfaces` become `previous_surfaces`, so panel-motion velocity
is based on accepted geometries only.

## 8. Imperial near-field force model

`src/backend/nearfield.jl` follows the Imperial College London C++ UVLM
operations named in `PROVENANCE.md`:

- dimensional Joukovski segment force
  `rho * delta_gamma * cross(velocity, r2-r1)`;
- unsteady panel force `-rho * area * normal * gamma_dot`;
- spanwise and chordwise circulation jumps;
- half-to-each-endpoint segment-force transfer; and
- quarter-to-each-corner unsteady transfer, with the Imperial trailing-edge
  rule.

Because the imported VortexLattice circulation direction is opposite to
Imperial's vertex traversal, the conversion
`gamma_imperial = -gamma_julia` is applied at the force boundary.

The following `System` fields hold the dimensional contributions:

| Field | Shape per surface | Contents |
|:--|:--|:--|
| `span_seg_forces` | `(nc+1, ns)` | Forces on spanwise vortex segments |
| `chord_seg_forces` | `(nc, ns+1)` | Forces on chordwise vortex segments |
| `unsteady_forces` | `(nc, ns)` | `gamma_dot` panel force |

Use `imperial_nodal_forces(system)` as the structural-coupling API rather than
reassembling these arrays in an input file.

## 9. Source-file map

### 9.1 Package entry and coupling infrastructure

| File | Responsibility and main functions |
|:--|:--|
| `src/WingPropellerUVLM.jl` | Declares the module, loads dependencies and source files in dependency order, and defines the exported public API. |
| `src/UVLMState.jl` | Transactional `UVLMSnapshot`, `snapshot_uvlm`, `restore_uvlm!`, `advance_uvlm_trial!`, and `commit_wake_rows!`. |
| `src/aeroelastic/GeneralizedAlpha.jl` | Generalized-alpha parameters, kinematics, structural corrector, coupled equilibrium residual, options, validation, and fixed-point step. |
| `src/aeroelastic/Excitations.jl` | Smooth Hann pulse scalar and generalized-load-vector helpers. |

### 9.2 Modified VortexLattice backend

| File | Responsibility and main functions |
|:--|:--|
| `src/backend/panel.jl` | `SurfacePanel`, `WakePanel`, and `TrefftzPanel` geometry; translation, rotation, reflection, normals, wake circulation, and Trefftz drag primitives. |
| `src/backend/geometry.jl` | Spacing types; wing/grid interpolation; grid-to-panel conversion; panel updates; lifting-line geometry; repeated trailing-edge detection. |
| `src/backend/system.jl` | `PanelProperties` and the global preallocated `System`, including circulation, surfaces, motion velocities, wakes, and force arrays. |
| `src/backend/reference.jl` | `Reference` quantities and the `Body`, `Stability`, and `Wind` output-frame marker types. |
| `src/backend/freestream.jl` | `Freestream`, frame transforms, translational/rotational velocity, derivatives, and trajectory conversion. |
| `src/backend/induced.jl` | Biot--Savart segment/ring velocities, surface/wake induced velocity, AIC assembly, and trailing-edge coefficient updates. |
| `src/backend/circulation.jl` | Surface-motion velocities, no-penetration right-hand side, circulation solve, and optional freestream derivatives. |
| `src/backend/wake.jl` | Shedding locations, wake-induced velocities, wake convection, row shifting, and shedding. |
| `src/backend/analyses.jl` | High-level steady/unsteady interfaces and the one-step `propagate_system!` algorithm. |
| `src/backend/nearfield.jl` | Production Imperial segmented, unsteady, and vertex-force calculation. |
| `src/backend/legacy_nearfield.jl` | Retained Imperial College-inspired legacy segment loads used by the Chang aeroelastic compatibility model. |
| `src/backend/nearfield_postprocessing.jl` | Generic body-force, force-history, lifting-line, derivative-reduction, and reference-frame postprocessing. |
| `src/backend/farfield.jl` | Trefftz-plane far-field induced drag. |
| `src/backend/stability.jl` | Body- and stability-axis derivatives from stored UVLM derivatives. |
| `src/backend/visualization.jl` | VTK/VTM/PVD output for surfaces, properties, wakes, and time histories. |
| `src/backend/rotors.jl` | Rotor/blade file readers, blade-lattice generation, `RotationMatrix`, mirroring, and airfoil selection. |
| `src/backend/vspgeom.jl` | OpenVSP degen-geometry reading and import. |
| `src/backend/nonlinear.jl` | Section data and optional nonlinear sectional VLM correction. This is aerodynamic nonlinearity, not AeroBeams structural nonlinearity. |

### 9.3 Wing--propeller helpers

| File | Responsibility and main functions |
|:--|:--|
| `src/wing_propeller/BladeGeometry.jl` | Tabulated chord/twist laws and nodal blade properties for the generic and Chang propellers. |
| `src/wing_propeller/GridUtilities.jl` | Linear interpolation, span attachment lookup, structurally deformed wing grids, and blade-grid generation. |
| `src/wing_propeller/Kinematics.jl` | Accepted-surface copying, reduced wing kinematics, initial/current propeller transformations, and in-place surface updates. |
| `src/wing_propeller/Initialization.jl` | `initialize_bohnisch_uvlm_system`, which constructs the complete wing/blade/wake system and reusable work buffers. |

### 9.4 Examples and validation files

| File | Responsibility |
|:--|:--|
| `examples/rectangular_wing_free_wake.jl` | Rigid unsteady wing, accepted-state stepping, VTK series, lift history, and mean spanwise lift. |
| `examples/README.md` | Rigid-case inputs, outputs, regression values, and environment overrides. |
| `examples/chang_linear_aeroelastic/chang_case.jl` | User-facing wing, propeller, and simulation configurations plus validation. |
| `examples/chang_linear_aeroelastic/run_chang_linear_aeroelastic.jl` | Readable primary workflow: case, structure, aerodynamics, excitation/coupling, solve, and output. |
| `examples/chang_linear_aeroelastic/src/chang_model_parameters.jl` | SI conversion and interpolation of Chang mass, inertia, stiffness, pylon, rotor, mesh, freestream, and time data. |
| `examples/chang_linear_aeroelastic/src/chang_structural_matrices.jl` | Point-mass and 3-D beam-element matrices and global assembly; `inertia_reference=:beam_axis` retains the historical comparison. |
| `examples/chang_linear_aeroelastic/src/chang_structural_model.jl` | High-level structural assembly, modal analysis/classification, and matrix diagnostics. |
| `examples/chang_linear_aeroelastic/src/chang_uvlm_coupling.jl` | Chang-specific motion transfer, load transfer, and rollback-safe aerodynamic callback. |
| `examples/chang_linear_aeroelastic/src/chang_simulation.jl` | Runtime options, coupled time marching, trim subtraction, accepted-wake commit, animation capture, and safety checks. |
| `examples/chang_linear_aeroelastic/src/chang_postprocessing.jl` | CSV/summary output, validation flags, history extraction, and optional plotting. |
| `examples/chang_linear_aeroelastic/src/chang_animation.jl` | Opt-in accepted-state recorder and 3-D GIF renderer. |
| `examples/chang_linear_aeroelastic/studies/` | Trim, airspeed, and two-stage convergence entry points. |
| `examples/chang_linear_aeroelastic/validation/` | Structural, virtual-work, force-model, damping, and inertia-remap checks. |

`Project.toml` defines package dependencies and compatibility. `test/runtests.jl`
contains package, force, snapshot/rollback, generalized-alpha, structural, and
example regression tests. `PROVENANCE.md`, `LICENSE`, and
`THIRD_PARTY_NOTICES.md` record the origin and licensing of imported code.

## 10. Public API reference

The functions below are exported by `WingPropellerUVLM.jl`. Functions not in
this list should be treated as backend implementation details unless a new
public interface and tests are deliberately added.

### 10.1 Geometry, panels, and rotors

| API | Use |
|:--|:--|
| `SurfacePanel`, `WakePanel`, `TrefftzPanel`, `Wake` | Core bound-surface, wake, and far-field representations. |
| `AbstractSpacing`, `Uniform`, `Sine`, `Cosine` | Select lattice point spacing. |
| `wing_to_grid` | Generate a lifting-surface vertex grid from leading-edge stations, chord, twist, and dihedral. |
| `grid_to_surface_panels` | Convert a `(3,nc+1,ns+1)` vertex grid to UVLM panels and control-point ratios. |
| `lifting_line_geometry`, `lifting_line_geometry!` | Obtain lifting-line coordinates and chords for sectional output. |
| `translate`, `translate!`, `rotate`, `rotate!`, `reflect`, `set_normal` | Transform panels/grids or control panel orientation. |
| `repeated_trailing_edge_points` | Detect shared trailing-edge vertices for consistent wake-velocity evaluation. |
| `generate_rotor` | Read rotor data and create blade grids/section data. |
| `RotationMatrix` | Create a 3-D axis rotation used by propeller kinematics. |
| `read_degengeom`, `import_vsp` | Import OpenVSP degen geometry. |

### 10.2 Reference, freestream, and system

| API | Use |
|:--|:--|
| `Reference(S,c,b,r,V,rho)` | Dimensional reference area, chord, span, point, velocity, and density. The five-argument form retains `rho=1`. |
| `Freestream(Vinf,alpha,beta,Omega)` | Freestream speed, angle of attack, sideslip, and body angular rate. Angles are radians. |
| `trajectory_to_freestream` | Convert a prescribed trajectory history to freestream states. |
| `AbstractFrame`, `Body`, `Stability`, `Wind` | Abstract output-frame interface and its body-, stability-, and wind-axis selectors. |
| `System(grids_or_surfaces; nw=...)` | Allocate the global UVLM state for all surfaces and maximum wake sizes. |
| `PanelProperties` | Per-panel circulation, velocity, and nondimensional force data. |
| `get_surface_properties` | Return `system.properties`. |

`Reference.rho` is the package's only aerodynamic density source. The
production Imperial near-field loads, legacy load routines, dynamic-pressure
normalization, and Trefftz-plane calculation all read it from the same
analysis reference. For the Chang case, `AIR_DENSITY = 1.225` is used only to
construct `ref`; subsequent solver, coupling-scale, and output calculations
read `ref.rho`.

### 10.3 Aerodynamic analyses and results

| API | Use |
|:--|:--|
| `steady_analysis`, `steady_analysis!` | Create/run a steady VLM solution. The bang form mutates a supplied `System`. |
| `unsteady_analysis`, `unsteady_analysis!` | High-level prescribed surface-history analysis. |
| `propagate_system!` | Solve one explicitly controlled UVLM step; use `advance_wake=false` for repeated coupling trials. |
| `advance_wake!` | Convect active wakes and shed one row from the accepted circulation/geometry without re-solving loads. |
| `imperial_nodal_forces` | Return dimensional vertex forces for structural transfer. |
| `imperial_nodal_positions` | Return the matching vortex-vertex positions. |
| `body_forces`, `body_forces_history` | Integrated coefficients/loads for one state or a stored history. |
| `spanwise_force_coefficients` | Spanwise coefficient distribution from panel properties. |
| `lifting_line_coefficients`, `lifting_line_coefficients!` | Sectional force and moment coefficients along lifting lines. |
| `far_field_drag` | Trefftz-plane induced drag. |
| `body_derivatives`, `stability_derivatives` | Aerodynamic derivatives when a derivative-capable analysis was run. |
| `write_vtk` | Write surface/wake state or histories as ParaView VTK/VTM/PVD output. |

### 10.4 Wing--propeller construction and motion

| API | Use |
|:--|:--|
| `linear_interpolate_1d` | Checked one-dimensional linear interpolation. |
| `span_position_to_node_index`, `span_positions_to_node_indices` | Map physical span stations to structural nodes. |
| `propeller_attachment_nodes_from_eta` | Map nondimensional attachment positions to structural nodes. |
| `generate_panel_grid_and_interpolate` | Deform a wing lattice by interpolated six-DOF nodal motion about an elastic axis. |
| `generate_aero_panel_grid_and_interpolate` | Related aerodynamic-grid interpolation helper retained for research inputs. |
| `generate_propeller_blades_grid` | Build reference blade grids from radius, chord, panel counts, twists, and blade count. |
| `get_chord_over_R`, `get_twist_deg`, `get_twist_deg_interp`, `get_nodal_properties` | Generic tabulated blade geometry. |
| `get_twist_deg_chang`, `get_nodal_properties_chang` | Chang-specific blade geometry. |
| `copy_surfaces_to_previous!` | Freeze accepted surface geometry before a time step. |
| `wing_kinematics_from_free_state` | Convert the earlier reduced three-DOF wing state to aerodynamic kinematics. |
| `initialize_propeller_grids!`, `update_propeller_grids!` | Build or update blade grids for attachment, wing motion, whirl, and spin. |
| `update_system_surfaces!` | Replace wing/blade panels in an existing `System`. |
| `initialize_bohnisch_uvlm_system` | Construct the complete Bohnisch/Chang wing--propeller system and work buffers. |

`initialize_bohnisch_uvlm_system` currently has a detailed keyword interface
because it preserves research-case inputs. Its named-tuple return exposes the
`system`, maximum/active wake vectors, surface indices, interaction groups,
reference/current blade grids, attachment geometry, and preallocated load
buffers. Use the call in `run_chang_linear_aeroelastic.jl` as the canonical
example and keep case-specific aliases outside the package.

### 10.5 Transactional time stepping

| API | Use |
|:--|:--|
| `UVLMSnapshot` | Deep copy of mutable history needed to repeat a time-step trial. |
| `snapshot_uvlm` | Capture accepted wakes, circulation, rates, surfaces, shedding points, wake counts, and freestream. |
| `restore_uvlm!` | Restore a snapshot in place while preserving array identities. |
| `advance_uvlm_trial!` | Restore a snapshot and run one UVLM trial; its `advanceWake` keyword defaults to `true` for compatibility. |
| `commit_wake_rows!` | Increment active wake lengths once after acceptance, capped at allocated maxima. |

### 10.6 Structural coupling and excitation

| API | Use |
|:--|:--|
| `generalized_alpha_parameters` | Compute second-order method parameters from `rho_infinity`. |
| `generalized_alpha_kinematics` | Recover end-step velocity and acceleration from an accepted displacement. |
| `generalized_alpha_corrector` | Solve one linear structural step for prescribed old/new aerodynamic and external loads. |
| `generalized_alpha_equilibrium_residual` | Check the complete interpolated structural equilibrium with the candidate aerodynamic load. |
| `PartitionedCouplingOptions` | Store iteration limit, state/load/equilibrium tolerances, and relaxation. |
| `partitioned_generalized_alpha_step` | Iterate structural correction and a user-supplied `load_at_state` callback to a fixed point. |
| `smooth_hann_pulse` | Unit-amplitude smooth finite-duration pulse. |
| `smooth_hann_pulse_load` | Put that pulse on selected generalized DOFs. |

### 10.7 Optional nonlinear aerodynamic correction

| API | Use |
|:--|:--|
| `SectionProperties` | Section grouping, area, airfoil, and circulation metadata. |
| `grid_to_sections` | Build section data from a grid and airfoil definitions. |
| `nonlinear_analysis!` | Apply the backend's nonlinear sectional aerodynamic correction. |

This optional aerodynamic correction is separate from the geometrically exact
nonlinear structural model in AeroBeams.

## 11. Chang runtime controls and outputs

The coupled input recognizes these environment variables:

| Variable | Default | Meaning |
|:--|:--|:--|
| `CHANG_END_TIME_S` | case `end_time_s` | Requested simulation duration |
| `CHANG_PLOT_RESULTS` | `true` | Load Plots.jl and create/display the history PNG |
| `CHANG_ANIMATE_WAKE` | `false` | Retain accepted wake states and create a GIF |
| `CHANG_OUTPUT_DIR` | example `output` directory | Output directory |
| `CHANG_OUTPUT_LABEL` | `chang_linear_imperial_uvlm` | Filename prefix |
| `CHANG_TRIM_REVOLUTIONS` | `10` | Wake-development revolutions before the perturbation |
| `CHANG_TRIM_AVERAGE_REVOLUTIONS` | `1` | Final revolutions used for mean baseline load |
| `CHANG_IMPULSE_MAGNITUDE` | `1000` | Pitch generalized-force pulse amplitude |
| `CHANG_IMPULSE_START_S` | end of trim | Explicit pulse start time |
| `CHANG_IMPULSE_DURATION_S` | `0.15` | Pulse duration |
| `CHANG_GA_RHO_INF` | `1.0` | High-frequency spectral radius |
| `CHANG_COUPLING_MAX_ITER` | `10` | Maximum fixed-point iterations per physical step |
| `CHANG_COUPLING_TOL_U` | `1e-5` | Scaled state residual tolerance |
| `CHANG_COUPLING_TOL_F` | `1e-2` | Scaled load-change tolerance |
| `CHANG_COUPLING_TOL_EQ` | `1e-10` | Linear corrector residual tolerance |
| `CHANG_COUPLING_TOL_COUPLED_EQ` | `1e-4` | Complete coupled equilibrium tolerance |
| `CHANG_COUPLING_RELAXATION` | `1.0` | Fixed displacement relaxation factor |
| `CHANG_PLOT_END_TIME_S` | `5` | Plot x-axis endpoint |
| `CHANG_ANIMATION_STRIDE` | `5` | Accepted physical steps between animation frames |
| `CHANG_ANIMATION_FPS` | `15` | GIF playback frame rate |

GIF generation is opt-in with `CHANG_ANIMATE_WAKE=true` because retaining every
selected surface/wake state can consume substantial memory and disk space.

A short, plot-free smoke run is:

```powershell
$env:CHANG_END_TIME_S = "0.012"
$env:CHANG_IMPULSE_START_S = "0.002"
$env:CHANG_IMPULSE_DURATION_S = "0.003"
$env:CHANG_PLOT_RESULTS = "false"
julia lib/WingPropellerUVLM/examples/chang_linear_aeroelastic/run_chang_linear_aeroelastic.jl
```

The driver writes `<label>_history.csv` and `<label>_summary.txt`; plotting adds
`<label>_history.png`. The CSV includes wing-tip response, each propeller's
pitch/yaw response, coupling iteration count, all fixed-point residuals, and a
convergence flag. `VALIDATION_PASSED=true` means the requested run completed,
all stored values were finite, and every accepted coupling step converged. It
does not by itself prove mesh or time-step convergence or agreement with an
experimental flutter boundary.

## 12. Extending the callback to AeroBeams

The UVLM side of the nonlinear coupling should preserve the same callback
contract:

```julia
function uvlm_load_at_aerobeams_state!(coupler, aero_snapshot, beam_state, time)
    restore_uvlm!(coupler.uvlm, aero_snapshot)
    update_geometry_from_aerobeams!(coupler, beam_state, time)
    propagate_system!(coupler.uvlm, coupler.freestream, coupler.dt;
        advance_wake=false, ...)
    return transfer_uvlm_loads_to_aerobeams(coupler, beam_state)
end

# After the outer aeroelastic iteration converges:
advance_wake!(coupler.uvlm, coupler.freestream, coupler.dt; ...)
commit_wake_rows!(coupler.active_wake_rows, coupler.maximum_wake_rows)
```

The changes are structural, not aerodynamic:

- obtain nodal positions and section rotations from the converged AeroBeams
  trial state rather than the linear Chang vector;
- use interpolation/rotation transfer that remains valid for large rotations;
- apply external nodal forces and moments through an AeroBeams-supported load
  interface;
- restore the AeroBeams beginning-of-step state before every outer trial;
- hold transferred UVLM loads fixed during an inner AeroBeams Newton solve; and
- commit both solvers only when motion, load, and nonlinear equilibrium
  residuals converge.

Never call the UVLM from each structural element residual or Newton evaluation.
Evaluate it at the outer coupling level, and never call `advance_wake!` inside
that iteration; the wake must advance according to physical time rather than
solver iteration count.

## 13. Common mistakes and diagnostics

| Symptom | Likely cause and check |
|:--|:--|
| `Package WingPropellerUVLM not found` | Run the example directly, activate `lib/WingPropellerUVLM`, or `develop` its local path. |
| Wake grows too quickly | Wake-row increment or propagation is inside the coupling iteration. Commit once after convergence. |
| Result changes with maximum coupling iterations even at tight tolerance | A trial is not restoring the same snapshot, or a mutable field is missing from the snapshot. |
| Large artificial damping/excitation | Check frame signs, rotation order, force/moment mapping, and moment reference points. |
| Loads are zero | Ensure `near_field_analysis=true` before calling `imperial_nodal_forces`. |
| Derivative error | Imperial near-field derivatives are unavailable; use `derivatives=false`. |
| Coupling stops as nonconverged | Reduce `dt`, lower relaxation, increase the iteration cap, then inspect state/load/coupled-equilibrium residuals separately. Do not simply disable a convergence test. |
| Correct total force but wrong torsion/whirl | Use `imperial_nodal_positions` and form moments about the deformed elastic axis, hub, pivot, and modal point consistently. |
| Interaction-off case still has blade--blade induction | This is expected: blades of one propeller share an interaction group. Cross-group influence is what is disabled. |
| VTK contains rejected geometries | Save histories only after the coupled step is accepted. |

## 14. Recommended verification sequence

Before using a new wing--propeller configuration for stability conclusions:

1. Run package tests.
2. Verify the structural matrices and natural modes without aerodynamics.
3. Run the rigid wing and each rigid propeller separately.
4. Compare interaction off with the sum of isolated components.
5. Enable interactions and check force/moment conservation.
6. Repeat one trial twice from one snapshot and confirm identical results.
7. Check that rejected iterations do not change active wake length.
8. Perform independent structural-mesh, UVLM-mesh, wake-length, vortex-core,
   azimuth-step, physical-time-step, and coupling-tolerance studies.
9. Reduce perturbation amplitude to verify a linear response before extracting
   damping or flutter onset.
10. For rotating systems, distinguish startup, rotor-periodic trim, and the
    later perturbation window in all response identification.

The linear Chang example is a useful regression and integration test. It is not
a substitute for the nonlinear AeroBeams limiting-case and convergence studies
described in the coupling blueprint.
