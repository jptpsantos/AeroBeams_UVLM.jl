# Nonlinear Chang wing–propeller analysis with UVLM and AeroBeams

**Technical implementation report**  
**Date:** 6 November 2026  
**Implementation:** `lib/WingPropellerUVLM/examples/chang_aerobeams_coupled`  
**Verification environment:** Julia 1.11.5 on Windows

## 1. Purpose and present status

This report explains the new code for a Chang wing–propeller system whose wing is represented by geometrically nonlinear AeroBeams elements and whose aerodynamics are evaluated by the existing unsteady vortex lattice method, WingPropellerUVLM. It describes the equations, their implementation, the load transfer, the solution procedure, the additions needed to combine the packages, and the verification completed so far.

Two structural models are available:

1. **Nonlinear wing with a two-coordinate propeller support**, selected with `solver.support=:reduced`. Each propeller has relative pitch and yaw coordinates. Its installed mass, moving-wing acceleration coupling, support stiffness and damping, and all four wing/support gyroscopic blocks are retained.
2. **Nonlinear wing with nonlinear beam pylons**, selected with `solver.support=:beam_pylon`. Pylon beam elements replace the reduced support. A physical rotor inertia and prescribed spin momentum act at each pylon tip. Rotor pitch and yaw then follow the beam-tip orientation.

Both models use explicit rotating blade surfaces in a common UVLM system. Wing–propeller aerodynamic interaction is enabled by default. The structural response retains the mean aerodynamic loading; it is not computed solely from deviations from a rigid aerodynamic baseline.

The implementation has passed component and coupled-step verification, a dry-wing modal comparison, and two short interacting transient calculations. A coarse mean-loaded equilibrium also converged. A more demanding fixed-shape wake-periodicity calculation did **not** meet its tolerance and was rejected. Consequently, the existing results establish implementation behavior and short-time integration, but do not establish a flutter speed, a nonlinear bifurcation diagram, or agreement with a Chang experiment.

### Contents

- [2. Scope of the additions](#2-scope-of-the-additions)
- [3. Code organization and execution](#3-code-organization-and-execution)
- [4. Frames, coordinates and signs](#4-frames-coordinates-and-signs)
- [5. Reuse of the Chang configuration](#5-reuse-of-the-chang-configuration)
- [6. Nonlinear wing construction](#6-nonlinear-wing-construction)
- [7. Two-coordinate support formulation](#7-two-coordinate-support-formulation)
- [8. Beam-pylon formulation](#8-beam-pylon-formulation)
- [9. Aerodynamic geometry](#9-aerodynamic-geometry)
- [10. Transfer of aerodynamic loads](#10-transfer-of-aerodynamic-loads)
- [11. Augmented mechanical equations](#11-augmented-mechanical-equations)
- [12. Time integration and strong coupling](#12-time-integration-and-strong-coupling)
- [13. Initialization and developed wakes](#13-initialization-and-developed-wakes)
- [14. Configuration and running cases](#14-configuration-and-running-cases)
- [15. Results and file formats](#15-results-and-file-formats)
- [16. Verification evidence](#16-verification-evidence)
- [17. Stability-study support](#17-stability-study-support)
- [18. Assumptions and current limits](#18-assumptions-and-current-limits)
- [19. Remaining work for a predictive study](#19-remaining-work-for-a-predictive-study)
- [20. Function reference and traceability](#20-function-reference-and-traceability)

## 2. Scope of the additions

### 2.1 What changed in the analysis framework

The implementation was added in a separate example directory. It reads the existing Chang configuration and physical distributions, imports AeroBeams and WingPropellerUVLM, and supplies a new structural adapter, propeller components, aerodynamic interface and coupled driver.

No changes were required in the existing AeroBeams beam formulation, WingPropellerUVLM source, Pazy analysis, original Chang analysis, convergence studies, package projects, or manifests. Their existing execution paths remain separate. This report also adds documentation only; it does not change those analyses.

“Modifications” in this report therefore means additions and local adaptations in the new example. It does not mean that all working-tree changes visible in the shared repository were produced by this implementation.

### 2.2 Why a new driver was needed

The existing UVLM–AeroBeams wing coupling supplies a basis for passing deformed wing geometry and aerodynamic loads. A propeller adds several requirements beyond a wing-only coupling:

| Requirement | Addition in this implementation |
|---|---|
| A rotor mounted on a moving nonlinear wing | Objective attachment kinematics and a separate mechanical component |
| Two relative support coordinates | Additional pitch/yaw unknowns and their time integration |
| Acceleration coupling | Full installed reference mass blocks and moving-carrier momentum transport |
| Prescribed shaft spin | Explicit spin angular momentum and gyroscopic terms |
| Wing plus rotating blades | One UVLM model containing the wing and every blade |
| Aerodynamic feedback in both directions | Interaction groups, updated surface grids and common wakes |
| Consistent generalized forces | Transpose wing interpolation and exact relative-coordinate virtual work |
| A flexible physical pylon | Connected AeroBeams beams and a physical tip-rotor component |
| A loaded initial condition | Fixed-shape periodic aerodynamic settling and nonlinear mean equilibrium |
| Reliable nonlinear time steps | Mechanical Newton iterations, strong FSI iterations and wake rollback |

Adding only a rotor point mass or setting an aerodynamic RPM would omit the independently spinning rotor's angular momentum. Adding only the two-by-two pitch/yaw gyro matrix would omit gyroscopic feedback through the moving wing attachment.

### 2.3 The main differences from the original Chang calculation

The original Chang module remains the source of reference physical data and reduced support matrices. The new model changes how the wing and the coupled time response are computed:

| Aspect | New example |
|---|---|
| Wing structure | Actual geometrically nonlinear AeroBeams mixed beam model |
| Wing mass representation | Distributed reconstruction of the full local Chang spatial inertia blocks |
| Wing section properties | Constant element compliance obtained by integrating the tabulated law |
| Support, first model | Chang reduced support extended to a finite-motion carrier |
| Support, second model | Connected nonlinear beams replacing the reduced support |
| Aerodynamic loads | Dimensional paired vortex forces and application points |
| Structural loading | Full mean and unsteady aerodynamic loads |
| Integration | Native AeroBeams first-order trapezoidal relations plus attachment/support Newmark reconstruction |
| Coupling | Strong partitioned geometry/load iterations with a final equilibrium check |
| Stability study | Phase-matched baseline/perturbed time responses and an RMS growth estimate |

The legacy `reference.integration`, `reference.coupling`, `reference.excitation` and `reference.output` settings are retained in the reused reference configuration object but do **not** select the algorithms in the new driver. The operative controls are the new top-level `solver`, `initialization`, `pylon` and `output` groups.

## 3. Code organization and execution

### 3.1 Files and responsibilities

| File | Responsibility |
|---|---|
| [chang_aerobeams_case.jl](chang_aerobeams_case.jl) | Constructs the editable case settings; does not execute a response |
| [run_chang_aerobeams.jl](run_chang_aerobeams.jl) | Runs the nonlinear-wing/reduced-support case |
| [run_chang_beam_pylon.jl](run_chang_beam_pylon.jl) | Runs the same case with beam pylons and a distinct output label |
| [run_speed_sweep.jl](run_speed_sweep.jl) | Sets a longer duration and runs paired trajectories at several speeds |
| [src/ChangAeroBeams.jl](src/ChangAeroBeams.jl) | Module entry point, imports, includes and public API |
| [src/AeroBeamsWing.jl](src/AeroBeamsWing.jl) | Wing properties, nonlinear model, node access, physical motion and nodal loads |
| [src/ReducedRotor.jl](src/ReducedRotor.jl) | Reduced installed mass, support elasticity/damping and gyro mechanics |
| [src/BeamPylon.jl](src/BeamPylon.jl) | Connected pylon beams, global node maps and physical hub poses |
| [src/TipRotor.jl](src/TipRotor.jl) | Rotor-only ordinary inertia and prescribed spin at a pylon tip |
| [src/AerodynamicInterface.jl](src/AerodynamicInterface.jl) | Surface geometry, UVLM setup and aerodynamic load maps |
| [src/CoupledSolver.jl](src/CoupledSolver.jl) | Defaults, problem storage, Newton/FSI solves, initialization and saving |
| [src/Studies.jl](src/Studies.jl) | Paired speed studies and perturbation-growth fitting |
| [validation/runtests.jl](validation/runtests.jl) | Focused verification-suite entry point |
| [validation/verification_report.md](validation/verification_report.md) | Compact record of completed verification |

Generated response data and plots are written under this example's `output/`. Its `.gitignore` excludes that generated directory from ordinary tracking. Links to those artifacts in this report require the local generated files to be present.

### 3.2 How the packages are loaded

`ChangAeroBeams.jl` places the repository root and local WingPropellerUVLM root on the process `LOAD_PATH`, imports the two packages, and includes the existing `ChangAeroelastic.jl` module under the alias `Chang`. The new source files are then included into `ChangAeroBeams`.

This is a local example-level extension. It does not install a package, edit an environment, or add a new global rotor type to the AeroBeams package.

### 3.3 Main call sequence

```text
run_chang_aerobeams.jl / run_chang_beam_pylon.jl
  include src/ChangAeroBeams.jl
  include chang_aerobeams_case.jl
  run_chang_aerobeams(case)
    build_chang_aerobeams(case)
      copy and resolve Chang reference data
      build nonlinear wing [and connected pylon beams]
      construct ReducedRotor [or TipRotor] components
      create wing + all-blade UVLM transaction problem
    solve_chang_aerobeams!(problem)
      initialize_mean_equilibrium! [if selected]
      for each physical time step
        snapshot accepted structural and support state
        begin UVLM transaction
        iterate aerodynamic trial / mechanical Newton solve
        check final geometry, velocity, loads and equilibrium
        commit one accepted wake step, or restore on failure
        record accepted response
    save_chang_aerobeams(problem) [if enabled]
```

### 3.4 What is stored in the problem

`ChangAeroBeamsProblem` owns the reference model, nonlinear wing/pylon model, vector of rotor components, aerodynamic problem, shared time vector, relative coordinates and rates, accepted nodal history, copied configuration, rotor phase offset, lifecycle flags and results dictionary.

The principal dimensions are:

| Quantity | Size or meaning |
|---|---|
| `wing.problem.x` | Native AeroBeams mixed unknown vector; not simply six displacements per node |
| `eta`, `eta_velocity`, `eta_acceleration` | `2*Npropellers` for reduced supports; empty for beam pylons |
| Augmented Newton vector `z` | `[x_AeroBeams; eta]` |
| Nodal load matrix | `6 × Nglobalnodes`, ordered force then physical moment |
| Interface node positions | `3 × Nnodes` |
| Node rotation data | One proper `3 × 3` section rotation per node |
| Aerodynamic surfaces | One wing plus `Npropellers*Nb_prop` explicit blades |

The builder deep-copies the supplied configuration before resolving it. Multiple runs can therefore use a common input configuration without changing its nested attachment arrays or other mutable data through the builder.

## 4. Frames, coordinates and signs

### 4.1 Global and section axes

All interface positions, physical velocities, forces and moments use the right-handed aerodynamic frame

$$
U=(x\text{ chord aft/downstream},\ y\text{ spanwise},\ z\text{ upward}).
$$

The wing section frame is

$$
b=(1\text{ spanwise},\ 2\text{ chord aft},\ 3\text{ downward}),
\qquad
R_0=R_{U\leftarrow b}=
\begin{bmatrix}0&1&0\\1&0&0\\0&0&-1\end{bmatrix}.
$$

The wing is created with the initial E321 orientation `[pi/2,0,pi]`; the code checks that its initial matrix equals this `R0`. The E321 values define the reference geometry. The evolving AeroBeams rotation unknowns are Wiener–Milenković (WM) parameters.

For a current WM vector `p`, the full wing section orientation is

$$
R_{\mathrm{section}}=R_{\mathrm{WM}}(p)R_0.
$$

The carrier deformation rotation needed by the rotor is

$$
D=R_{\mathrm{section}}R_0^T=R_{\mathrm{WM}}(p),
$$

which is initially the identity. Using the full section matrix as the initial rotor orientation would rotate the shaft into the wrong reference direction.

### 4.2 Relative support and spin conventions

For each reduced support,

$$
\eta=\begin{bmatrix}\theta\\\psi\end{bmatrix},\qquad
R_{\mathrm{rel}}=R_y(\theta)R_z(-\psi),\qquad
R_{\mathrm{spin}}=R_x(-\Omega t).
$$

Thus positive structural yaw is opposite a positive rotation about global `U-z` in the undeformed configuration. The signed shaft spin momentum is initially

$$
H_0=-h e_x,\qquad h=J_s\Omega.
$$

These conventions control the signs in the gyro matrices, relative aerodynamic moments and blade motion. Reversing only the structural `spin_scale` reverses the structural gyro terms. A physical reversal of blade rotation must also change the aerodynamic motion consistently; the case interface does not currently implement that complete operating variation automatically.

### 4.3 Physical angular velocity is not a WM parameter rate

`wm_spatial_angular_kinematics` differentiates the actual rotation along

$$
p(t)=p+t\dot p+\tfrac12t^2\ddot p.
$$

Using `ForwardDiff`, it obtains

$$
\omega=\operatorname{vee}(\dot R R^T),\qquad
\alpha=\operatorname{vee}(\ddot R R^T+\dot R\dot R^T).
$$

The second expression is the time derivative of the spatial angular velocity matrix. This conversion retains the finite-rotation terms required by the rotor mechanics. Treating `pdot` as `omega` would change the gyro and inertia loads at finite rotation.

### 4.4 Equilibrium signs and units

Interface wrenches are `[Fx,Fy,Fz,Mx,My,Mz]` in `U`, with force in N and moment in N m. Relative angles are in radians internally. CSV and plot angle outputs use degrees.

`rotor_mechanics` returns the **required inertial wrench on the left side of equilibrium**. Its negative is the reaction exerted on the supporting structure. The driver adds the returned wrench to the native mechanical residual and applies aerodynamic loads on the external-load side. This convention prevents adding the same rotor reaction twice.

## 5. Reuse of the Chang configuration

### 5.1 Data read from the original implementation

The new builder calls the existing `Chang.load_chang_configuration` and `Chang.build_chang_model`. These supply the geometry, tabulated structural distributions, blade geometry/twist, support coefficients, operating point and aerodynamic settings. The original Chang time-response solver is not called.

The principal reused definitions are in:

- [chang_case.jl](../chang_linear_aeroelastic/chang_case.jl): reference values and editable groups.
- [chang_configuration.jl](../chang_linear_aeroelastic/src/chang_configuration.jl): resolution and validation of those groups.
- [chang_model_parameters.jl](../chang_linear_aeroelastic/src/chang_model_parameters.jl): physical and derived parameters.
- [chang_structural_matrices.jl](../chang_linear_aeroelastic/src/chang_structural_matrices.jl): original reduced installed mass and gyro blocks.
- [ChangAeroelastic.jl](../chang_linear_aeroelastic/src/ChangAeroelastic.jl): original module entry point.

The resolver is called with an empty environment dictionary. The new example therefore uses the explicit supplied settings rather than silently taking legacy `CHANG_*` environment overrides.

### 5.2 Important derived support quantities

The inherited model uses pylon length `L`, pylon mass `m_p`, and rotor blade mass `m_r` to form

$$
m_P=m_p+m_r,\qquad m_p=\mu_p L,\qquad m_r=N_bm_b,
$$

$$
S_{\mathrm{root}}=-\left(m_rL+\frac{m_pL}{2}\right),\qquad
S_{\mathrm{modal}}=-\left(\frac{m_rL}{2}+\frac{m_pL}{6}\right),
$$

$$
I_{\mathrm{root}}=m_rL^2+\frac{m_pL^2}{3},\qquad
I_{\mathrm{cross}}=\frac{m_rL^2}{2}+\frac{m_pL^2}{8}.
$$

The calibrated relative inertias and polar spin inertia are

$$
I_{\theta}=\frac{K_{\theta}}{(2\pi f_{\theta})^2},\quad
I_{\psi}=\frac{K_{\psi}}{(2\pi f_{\psi})^2},\quad
J_s=I_{x,\mathrm{prop}}=\frac{K_{\mathrm{twist}}}{(2\pi f_{\mathrm{twist}})^2}.
$$

The last relation is a reference calibration. It is not a new geometric integration of the explicit UVLM blades. The new code preserves this distinction when using `Js` for the beam-tip rotor.

### 5.3 Fixed RPM and fixed advance ratio

The inherited parameter calculation scales `Omega` with flight speed when trim speed is held fixed. The new `solver.speed_policy` makes the operating choice explicit:

- `:fixed_rpm`: before resolving the copied reference, set its trim speed equal to the current freestream speed. The resulting `Omega` equals the requested RPM at each speed.
- `:fixed_advance_ratio`: preserve the inherited scaling, so `Omega` changes with freestream speed.

The choice affects aerodynamic loads and `h=Js*Omega`, and must be held consistent when interpreting a speed sweep.

### 5.4 Attachment stations

The inherited Chang parameters select the nearest structural node to each requested `attachment_eta`. The new interface uses that selected node's actual span station. It rejects an inconsistent attachment map instead of interpolating a mechanical rotor onto an unsupported location.

For the default 20-element uniform wing, the requested `eta=0.65` is a node. With another element count, the actual location can differ from the requested fraction. Check `reference.parameters.prop_attach_nodes` and `reference.parameters.span_nodes` when refining the mesh.

## 6. Nonlinear wing construction

### 6.1 Actual AeroBeams elements

`build_wing` constructs an AeroBeams `Beam`, fixes all six root coordinates, creates a numeric six-component load boundary condition at every node, and creates an AeroBeams `DynamicProblem`. The beam has `aeroSurface=nothing`: UVLM supplies the aerodynamic loading, so native sectional aerodynamics are not added a second time.

The model uses zero gravity and a stationary reference frame. Freestream motion belongs to the aerodynamic operating point. The wing is geometrically nonlinear with an elastic sectional law; this example does not add nonlinear material constitutive behavior.

The current geometric adapter requires a straight elastic axis with constant global `x`. It explicitly checks the leading-edge and chord distributions against that requirement. It does not silently represent a swept or curved elastic axis with the same straight beam.

### 6.2 Reconstructing distributed inertia

The legacy wing supplies a conservative local `6×6` spatial inertia block `B_i` at each source node. `_wing_control_lengths` assigns half-element widths at the ends and half the distance between adjacent neighbor stations internally. `_wing_element_inertias` forms nodal densities

$$
\mathcal I_i=B_i/\ell_i,
$$

interpolates those densities linearly between source stations, integrates them exactly over every target element, and divides by the target element length:

$$
\bar{\mathcal I}_e=\frac{1}{L_e}\int_e\mathcal I(y)\,dy.
$$

It symmetrizes each result, checks positive definiteness and verifies

$$
\sum_e L_e\bar{\mathcal I}_e=\sum_iB_i
$$

within a relative tolerance of `1e-12`. This retains the complete sum of the **local** blocks, including local eccentricity and products of inertia. It avoids reducing the wing to a scalar mass per unit length.

This conservation statement does not imply identical nodal lumping, identical global spatial inertia about every origin, identical spanwise mass moments, or identical modal matrices. Moving a distributed density between stations changes its spatial weighting. The dry-wing modal comparison is therefore a required check rather than an assumed equality.

### 6.3 Section stiffness and compliance

`_wing_element_stiffness` evaluates the inherited constitutive distribution and transforms it into the AeroBeams section convention. Its nonzero structure is

$$
S_{11}=EA,\quad S_{22}=S_{33}=10^{12},\quad S_{44}=GJ,
$$

$$
S_{55}=EI_y,\quad S_{66}=EI_z,\quad S_{56}=S_{65}=-EI_{zy}.
$$

The assignment of the two bending entries follows the source distribution and section axes in the code. The negative cross-stiffness follows the physical curvature mapping `kappa_2=-w_down''`, `kappa_3=v_chord''`. It should not be replaced by an unexamined equality between the original Hermite slope variables and physical rotation components.

AeroBeams uses constant compliance within each element. The adapter splits integration intervals at every tabulated material station and uses three-point Gauss quadrature to compute

$$
\bar C_e=\frac{1}{L_e}\int_eS(y)^{-1}\,dy,\qquad S_e=\bar C_e^{-1}.
$$

This preserves integrated compliance for constant section resultants. It is different from simply averaging stiffness entries. The finite shear stiffness approximates the Euler–Bernoulli limit and remains an explicit numerical approximation to investigate during refinement.

### 6.4 Nodes, rotations and loads

Every node has a named load BC, `Chang wing load node N`. This makes it a special node with an explicit location in the mixed unknown and equation maps.

`wing_node_parameters` reads the current node displacement and WM vector directly from `problem.x`. The fixed root is handled separately because its corresponding unknowns are reactions, not free displacements. `wing_node_poses` combines reference positions/directors with the current displacement/rotation. Thus aerodynamic trials use current Newton unknowns, rather than stale element output from the previous accepted state.

`set_wing_loads!` updates numeric physical force and moment BC values and invokes `AB.update_BC_data!`. Native AeroBeams assembly handles its own parameterization-specific load transformations. The coupling interface supplies physical spatial moments, not a guessed generalized moment conjugate to `p`.

When pylons are present, the code retains both `wing_node_count` and total `node_count`. Wing aerodynamics and wingtip output use the former; structural assembly and loading use the latter.

## 7. Two-coordinate support formulation

### 7.1 Why the moving attachment needs more than a two-by-two model

A fixed-base pitch/yaw support can be written as

$$
M_\eta\ddot\eta+(D_\eta+G_{\eta\eta})\dot\eta+K_\eta\eta=Q_\eta.
$$

On a flexible wing, the base translates and rotates. Rotor/support inertia produces wing loads, wing accelerations drive the relative coordinates, and the spin momentum couples both motions. `ReducedRotor` therefore operates on eight local quasi-velocities: three carrier translations, three carrier rotations and two relative-coordinate rates.

### 7.2 Reference installed mass

The component reconstructs each propeller's installed mass from the physical reference parameters, rather than extracting an attachment block that could already contain another coincident propeller. In the original structural axes `S=(span,downstream,down)`, it forms

$$
M_S=\begin{bmatrix}A&B^T\\B&M_\eta\end{bmatrix},
$$

where `A` contains the three translational masses, installed first moments and transverse carrier inertias. The source coupling block is

$$
B=\begin{bmatrix}
0&0&S_\theta&I_{\theta\alpha}&0&0\\
-S_\psi&0&0&0&0&I_{\psi\gamma}
\end{bmatrix}.
$$

The axes are transformed with

$$
Q=\operatorname{diag}(R_{U\leftarrow S},R_{U\leftarrow S},I_2),\qquad
M_0=QM_SQ^T,
$$

where `R_(U←S)=R0`. The two relative coordinates retain their inherited signs.

`M0` preserves the source Chang component, including its null carrier shaft-roll direction. For the reference parameters verified here, it is positive semidefinite as a standalone installed component, rather than a general positive-definite rigid rotor tensor. The constructor does not enforce positive semidefiniteness for arbitrary changes to the calibrated support frequencies or inertias; recheck the complete installed matrix after changing those data. The full wing still supplies its own structural inertia. The separate `TipRotor` model validates a physical positive-definite ordinary inertia tensor.

### 7.3 Objective extension to finite carrier motion

The reference mass implies a center of mass `c0` and relative CM mode matrix `C_eta`, obtained from its translation–rotation and translation–relative blocks. The implementation uses

$$
c(\eta)=c_0+C_\eta\eta,
$$

$$
L(\eta)=\begin{bmatrix}I&-[c(\eta)]_\times&C_\eta\end{bmatrix},\qquad
M(\eta)=M_0+m_P\left(L^TL-L_0^TL_0\right).
$$

Here `[a]_× b=a×b`. This changes the translational kinetic energy consistently with the moving CM while freezing the intrinsic reduced/modal inertia at the reference configuration. The code evaluates both shape derivatives `M_,theta` and `M_,psi` analytically.

For carrier deformation rotation `D`, spatial base velocity `v` and spatial angular velocity `omega`, define

$$
u=\begin{bmatrix}D^Tv\\D^T\omega\\\dot\eta\end{bmatrix},\qquad
T=\tfrac12u^TM(\eta)u,
$$

$$
\dot u=\begin{bmatrix}D^Ta-\omega_b\times v_b\\D^T\alpha\\\ddot\eta\end{bmatrix},\qquad
\dot M=\sum_iM_{,\eta_i}\dot\eta_i.
$$

The transport term in the first row is necessary because the components are expressed in a moving frame.

Let `mu=M*u` and `mu_dot=M*udot+Mdot*u`. With translational momentum `P=mu[1:3]` and angular momentum `Pi=mu[4:6]`, the ordinary inertial loads in carrier coordinates are

$$
f_b=\dot P+\omega_b\times P,
$$

$$
m_b=\dot\Pi+\omega_b\times\Pi+v_b\times P,
$$

$$
r_{\eta,i}=\dot\mu_{6+i}-\tfrac12u^TM_{,\eta_i}u
 +(K_\eta\eta+D_\eta\dot\eta)_i.
$$

The physical base wrench returned to the driver is `[D*f_b; D*m_b]`, after the gyro and optional gravity terms are included. These equations retain moving-base acceleration coupling and momentum transport while allowing finite wing motion.

The support remains a reduced small-relative-deformation model. Exact rotation matrices for blade motion and spin do not make its linear CM reconstruction and frozen intrinsic modal inertia a geometrically exact pylon. The default guard is `max(abs(eta)) <= 0.35 rad`.

### 7.4 Support stiffness and damping

The constructor obtains `K_eta` from the original `Ks_P` block. It takes the symmetric part of `Cs_P` as physical viscous damping, separating the antisymmetric spin terms that it reconstructs explicitly. The constructor can include the inherited stiffness-proportional contribution, but the coupled builder rejects nonzero global Chang Rayleigh damping because it is not implemented as a nonlinear beam material law.

In the unified coordinates `[wing base; relative support]`, spring potential depends on `eta` and viscous dissipation depends on `eta_dot`. Their generalized contributions therefore enter the relative equations. Adding the same mount spring moment as a second independent wing load would double count an internal interaction. Wing feedback already follows from the coupled inertial and external-load equations.

### 7.5 Full gyroscopic coupling

`rotor_spin_data` computes the exact carrier-frame spin momentum

$$
H(\eta)=-hR_y(\theta)R_z(-\psi)e_x
=\begin{bmatrix}
-h\cos\theta\cos\psi\\h\sin\psi\\h\sin\theta\cos\psi
\end{bmatrix}.
$$

Writing `H_eta=∂H/∂eta`, the eight-by-eight gyro operator has blocks

$$
G_{\omega\omega}=-[H]_\times,\quad
G_{\omega\eta}=H_\eta,\quad
G_{\eta\omega}=-H_\eta^T,
$$

$$
G_{\eta\eta}=h\cos\psi
\begin{bmatrix}0&1\\-1&0\end{bmatrix},
$$

with zero translation rows/columns. `G*u` is added to the ordinary inertia residual. These are the carrier/carrier, carrier/relative, relative/carrier and relative/relative blocks. At the undeformed reference, they reproduce the original Chang blocks after the stated axis transformation.

In particular,

$$
G_{\eta\eta}(0)=\begin{bmatrix}0&h\\-h&0\end{bmatrix}.
$$

Although a linear second-order equation groups this with a damping-like velocity matrix, it is gyroscopic coupling: `G^T=-G`, so `u^T G u=0`. It does no dissipative work. Physical support damping is stored separately in `D_eta`.

The verification checks all four blocks, their residual derivatives, sign reversal, finite-carrier objectivity, energy/power and fixed-base whirl frequencies. It does not infer correctness from the two-by-two relative block alone.

## 8. Beam-pylon formulation

### 8.1 Topology and replacement of the reduced support

`build_wing_with_beam_pylons` first builds the wing and then creates one AeroBeams pylon per propeller, connected at the selected physical wing node. The pylon reference axis points upstream along `-U-x`:

$$
(b_{p1},b_{p2},b_{p3})=(-e_x,e_y,-e_z).
$$

The initial E321 parameters are `[0,pi,0]`. A preliminary model resolves global node numbering before the final named load BCs and dynamic problem are created.

In this variation:

- There are no additional `eta` unknowns.
- The reduced support's installed mass, `K_eta` and `D_eta` are not assembled.
- The pylon beams supply their own distributed inertia and elasticity.
- A `TipRotor` at each physical tip supplies rotor-only ordinary inertia and spin.
- No duplicate AeroBeams `PointInertia` is added.

The connection is an assembled beam joint. The rotor carrier follows the tip; an additional flexible pitch/yaw bearing is not included in series with that beam.

### 8.2 Pylon section properties and their meaning

The pylon section matrix is diagonal:

$$
S_p=\operatorname{diag}(EA,GA,GA,GJ,EI_{\mathrm{pitch}},EI_{\mathrm{yaw}}).
$$

The default bending properties are

$$
EI_{\mathrm{pitch}}=K_\theta L,\qquad EI_{\mathrm{yaw}}=K_\psi L.
$$

For a cantilever under a pure tip moment, the tip rotation is `M*L/EI`; this choice matches the reduced static rotation under that loading. It does not match displacement under a tip force, distributed mode shape, modal inertia, dynamic frequency or whirl boundary automatically.

Default numerical assumptions are `EA=GA=1e10 N`, `GJ=1e8 N m²`, four elements per pylon and sectional rotary inertia per length `1e-3 kg m`. Pylon translational density uses the Chang mass-per-length datum. These properties are configurable and require measured or independently justified values for a physical pylon study.

The implemented beam-pylon variation has no structural material damping. A nonzero inherited reduced support damping ratio is rejected, as is global inherited stiffness-proportional damping.

### 8.3 Rotor ordinary inertia and independent spin

`TipRotor` places its CM at the hub. Its default rotor mass is `Nb*blade_mass`, excluding all pylon mass. For the reference four blades at 1.44 kg each, this is 5.76 kg. Stationary nacelle mass is not introduced without user data.

The default spin inertia is the calibrated `Ix_prop`, approximately `2.63146 kg m²` for the reference data. The ordinary carried tensor defaults to

$$
J=\operatorname{diag}(J_s,J_s/2,J_s/2),
$$

which is an axisymmetric thin-disk assumption. It is not a measured Chang transverse inertia. `rotor_inertia_tensor`, `rotor_transverse_inertia`, `rotor_mass` and `spin_inertia` allow separate carried assembly and independently spinning rotor properties. The constructor checks mass/inertia finiteness, symmetry, positive definiteness and principal-inertia triangle inequalities.

The ordinary tensor is carried by the tip frame. The default represents an axisymmetric mean rotor; the implementation does not rotate a nonaxisymmetric blade-resolved inertia tensor through azimuth. A supplied full tensor should be interpreted consistently as the effective carried tensor used by this model.

### 8.4 Tip equations

For incremental tip rotation `D`, define

$$
J_U=DJD^T,\qquad H_U=D[-h,0,0]^T.
$$

The component residual wrench is

$$
F=m(a-g),
$$

$$
M=J_U\alpha+\omega\times(J_U\omega+H_U).
$$

The first two rotational terms are ordinary rigid-body inertia; `omega×H_U` is the prescribed-spin gyro term. The latter is necessary even though the pylon already has nonlinear beam inertia. The component API supports gravity, but the current example driver passes zero gravity.

Aerodynamic rotor loads act once at the physical tip. The nonlinear beam carries their resulting force and moment to the wing. Applying an additional equivalent root wrench on top of the tip wrench would double count that aerodynamic load.

## 9. Aerodynamic geometry

### 9.1 One interacting multi-surface problem

`_ca_build_aerodynamic` creates a `UVLMModel`, operating point, solver and `UVLMDynamicProblem`. One surface represents the wing and each blade is a separate physical surface with its own wake. The blade mesh remains rigid within each blade; its pose changes with shaft spin and its supporting structure.

Interaction group 1 is the wing, group `ip+1` contains the blades of propeller `ip`. Physical surface IDs are unique even when several blades share one interaction group. The interaction switch comes from `reference.simulation.interaction_on`. Image symmetry is disabled by default for the asymmetric installation; explicit multiple propeller stations can be supplied.

The initial adapter uses the corrected `:imperial` near-field force model and, through the normal builder, a fixed positive core radius. Legacy force selection and factor-only core selection are not silently substituted. The aerodynamic solver retains the inherited UVLM assumptions; this extension does not add viscous separation or a stall model.

### 9.2 Local compatibility adaptations

Two facade behaviors required local handling:

1. The existing convenience constructor's vector symmetry broadcast can produce a `BitVector`, whereas its model field requires `Vector{Bool}`. The new adapter allocates through the scalar `symmetric=false` path, then sets the allocated wing/blade flags.
2. The facade initially uses interaction-group IDs as physical surface IDs. The adapter sets unique physical IDs in the accepted, trial and commit system buffers. The interaction groups remain available separately for interference control.

These adaptations operate on the new problem objects. No UVLM core file was edited. Preserving the two kinds of IDs matters because finite-core/source identity decisions and interaction-group masking serve different purposes.

### 9.3 Wing surface motion

The aerodynamic and structural span stations may differ. `_ca_wing_geometry` finds linear weights `w_ij` from structural nodes to each aerodynamic station. At chord fraction `f`, node `i` contributes the actual chord point

$$
x_i(f)=r_i+(f-f_{EA})c_i d_{2i},
\qquad
x_a(f)=\sum_iw_{ia}x_i(f),
$$

where `r_i` is the current elastic-axis position and `d_2i` is the current physical chord director. Finite rotations enter before interpolation; rotation matrices themselves are not linearly interpolated and then treated as exact orientations.

UVLM moves bound vortex rows to quarter-panel locations. The load map therefore stores the corresponding fractions

$$
f^v_k=0.75f_k+0.25f_{k+1}
$$

for all but the last row, which remains at the trailing edge. This keeps the motion/load pair consistent with the force locations produced by the backend.

### 9.4 Reduced-support blade motion: physical hub versus modal center

Let `r_A` be the wing attachment point, `d_h` the full physical hub offset, and `d_m` the inherited modal load arm. The default offsets correspond to `(-L,0,0)` and `(-L/2,0,0)`. Define

$$
r_C=r_A+D(d_h-d_m),\qquad
r_h=r_C+DR_{rel}d_m.
$$

For a blade-reference point `rho`, the actual point is

$$
x_b=r_C+DR_{rel}(d_m+R_{spin}\rho).
$$

At zero relative angle, the hub is at `r_A+D*d_h`, preserving the full physical pylon length. Relative motion rotates about the effective modal center `r_C`, preserving the reduced formulation's arm. Replacing both offsets with the same arm would alter either the physical geometry or the reduced-coordinate kinematics.

The half-length modal arm is part of the inherited Chang assumed-mode definition. Its structural first moments and cross-inertias retain those assumed-mode coefficients. Changing `hub_load_arm_factor` alone changes the aerodynamic motion/load map without recalibrating the structural reduced shape. Keep its default 0.5 for the reference regression, or derive and calibrate consistent structural support data when changing that assumption.

The instantaneous relative axes are

$$
a_\theta=De_y,\qquad a_\psi=-DR_y(\theta)e_z.
$$

### 9.5 Beam-pylon blade motion

For a pylon-tip hub pose `(r_h,D_h)`,

$$
x_b=r_h+D_hR_{spin}\rho.
$$

`pylon_hub_poses` obtains the physical tip position and the incremental WM tip rotation, initially the identity in `U`. There is no separate modal center, half-length arm or relative support rotation in this variant.

## 10. Transfer of aerodynamic loads

### 10.1 Dimensional paired forces and positions

The code uses `loads.forces` together with `loads.positions`, rather than assuming a force is applied at a mesh point chosen independently. Each rotor's force and moment about its physical hub are

$$
F_h=\sum_a f_a,\qquad M_h=\sum_a(x_a-r_h)\times f_a.
$$

Moments are computed from the actual force application locations. They include the effect of force offsets; they are not reconstructed from a prescribed sectional pitching coefficient alone.

### 10.2 Wing loads as the transpose of the motion map

For each aerodynamic wing force, the map applies

$$
F_i\mathrel{+}=w_{ia}f_a,\qquad
M_i\mathrel{+}=w_{ia}\big[(f^v_a-f_{EA})c_i d_{2i}\big]\times f_a.
$$

This is the transpose of the chord-point motion interpolation with respect to physical node translation and rotation. It satisfies virtual work

$$
\sum_af_a\cdot\delta x_a
=\sum_i(F_i\cdot\delta r_i+M_i\cdot\delta\phi_i).
$$

The nodal moment arm is the contributing node's own chord offset. Adding an extra aerodynamic-station-to-node spanwise arm would not be the transpose of this interpolation and would introduce spurious work. Global wrench balance follows from the weighted chord-point construction.

### 10.3 Reduced-support generalized loads

The carrier gets the equivalent physical wrench at `r_A`:

$$
F_A=F_h,\qquad M_A=M_h+(r_h-r_A)\times F_h.
$$

The relative-coordinate moment about the modal center is

$$
M_C=M_h+(r_h-r_C)\times F_h,
$$

$$
Q_\theta=a_\theta\cdot M_C,\qquad Q_\psi=a_\psi\cdot M_C.
$$

Both carrier and relative generalized loads are needed because the same blade point moves under both independent coordinate variations. They are not two applications of the same physical load. Double counting would instead occur if an additional full support reaction were applied on top of this unified coordinate residual.

### 10.4 Beam-pylon loads

`pylon_aerodynamic_loads` returns wing-only nodal aerodynamic loads, a `6×Npropellers` matrix of hub wrenches, and an empty relative-load vector. `_aero_load_map` inserts each hub wrench into the global pylon-tip node column. Internal beam equilibrium supplies the wing reaction.

The verification checks virtual work at finite rotations for both formulations, including the reduced modal center, the pylon-tip map, independent span meshes and the shifted vortex force locations.

## 11. Augmented mechanical equations

### 11.1 Residual assembly

At frozen aerodynamic loads, the unknown is

$$
z=\begin{bmatrix}x_{AB}\\\eta\end{bmatrix}.
$$

`coupled_mechanical_arrays!` first updates the AeroBeams nodal load BCs and calls `AB.assemble_system_arrays!`. It then adds the rotor correction:

$$
R(z)=\begin{bmatrix}R_{AB}(x_{AB};F_{aero})\\0\end{bmatrix}
+\begin{bmatrix}\sum_j\mathcal A_jW_{inertia,j}\\r_\eta-Q_\eta\end{bmatrix}.
$$

The native residual already contains the wing and, when selected, pylon beam mechanics. The correction contains only the additional installed rotor/support mechanics.

`_insert_wrench_residual!` uses explicit special-node equation maps and native `forceScaling`. Shared joint equilibrium rows are updated once; separate connected-element compatibility rows are not mistaken for another place to add the same wrench.

### 11.2 Jacobian and Newton solve

The Jacobian combines the native AeroBeams analytic matrix with a central finite difference of the small component correction. Perturbed columns are limited to the free attachment-node translation/rotation coordinates and the added relative coordinates. The finite-difference step is

$$
\delta z_j=\sqrt[3]{\epsilon_{mach}}\max(|z_j|,10^{-3}).
$$

This includes rotor acceleration, mass-shape, finite-rotation and gyro derivatives without differentiating the full UVLM solve inside mechanical Newton. The aerodynamic loading remains frozen during that inner solve.

`_solve_mechanics!` uses `delta_z=-(J\R)` and a halving line search with up to 12 trial lengths. It accepts a Newton solution only when the absolute infinity norm of the assembled scaled residual is at most `1e-8` by default. Failure to reduce the residual, nonfinite values, or exhaustion of 25 iterations raises an error.

The assembly uses a dense augmented matrix and linear solve. This is adequate for verification-size models but is a performance limitation for substantially larger meshes or many rotors. Sparse block assembly or component automatic differentiation is a future optimization; it is not present in the current example.

### 11.3 What is added to AeroBeams, and where

The rotor contribution is added by this external augmented driver. It is not registered as a new core AeroBeams constitutive object. Native beam assembly and state updates are reused, but `AB.solve!(wing.problem)` by itself is not the complete new wing–propeller simulation.

Similarly, an ordinary AeroBeams eigenproblem built from only `wing.problem.model` would not automatically include this external reduced support, its full gyro residual, or UVLM wake states. The dry-wing eigenproblem used for verification is intentionally only a dry-wing structural comparison.

## 12. Time integration and strong coupling

### 12.1 Time grid and kinematic reconstruction

The default time step comes from the inherited azimuth step. An explicit `solver.time_step_s` can override it. The driver uses a uniform grid ending at the largest whole step not exceeding the requested duration. Therefore the final recorded time can be slightly shorter than `end_time_s`.

Native AeroBeams element equations use their first-order trapezoidal equivalent-rate relations. For the additional support coordinates and attachment-node parameter history, the driver reconstructs average-acceleration Newmark rates, `gamma=0.5`, `beta=0.25`:

$$
\dot q_{n+1}=\frac{2}{\Delta t}(q_{n+1}-q_n)-\dot q_n,
$$

$$
\ddot q_{n+1}=\frac{4}{\Delta t^2}(q_{n+1}-q_n)
-\frac{4}{\Delta t}\dot q_n-\ddot q_n.
$$

Here `q` is a translation, WM parameter vector or relative angle as appropriate. WM rates are converted to physical angular quantities as described in Section 4. The added nodal history is held fixed throughout every trial at the same physical time step.

This mixed discretization is the method actually implemented. Selecting a legacy generalized-alpha option does not change it, and an exact fully discrete energy-conservation result is not asserted for the combined FSI algorithm.

### 12.2 A physical time step

For every target time, the driver:

1. Copies the accepted native unknowns, native element states/rates, structural clocks, support state and accepted nodal history.
2. Updates the native structural time variables and computes equivalent rates once from accepted history.
3. Predicts node translations, WM parameters and relative angles.
4. Calls `UVLM.begin_time_step!` once.
5. Evaluates aerodynamic trial loads at the current predicted/relaxed geometry.
6. Applies the prescribed support pitch pulse, if active.
7. Solves the frozen-load augmented mechanical equations with Newton iterations.
8. Compares aerodynamic geometry, attachment velocity and load changes.
9. Re-evaluates aerodynamics at the solved geometry and checks the full mechanical residual with these final loads.
10. If all checks pass, updates accepted structural/support history, commits UVLM once, and records the step. Otherwise, it relaxes the mixed state and repeats, up to the configured FSI limit.

This is a strong partitioned fixed-point solve with an inner nonlinear mechanical solve. It is not a single monolithic UVLM/beam Newton system.

### 12.3 Perturbation pulse

The applied moment is

$$
M_p(t)=M_0\sin^2\left(\pi\frac{t-t_0}{T_p}\right)
\quad\text{for }t_0<t<t_0+T_p,
$$

and zero otherwise. The default is 2 N m for 0.02 s, beginning at 0.15 s.

For reduced supports, it enters the selected relative pitch equation. For beam pylons, it is a physical moment at the selected tip along the current hub `y` direction. These are different support formulations and the pulse acts through their respective coordinates.

`hub_forces` and `hub_moments` in results are the aerodynamic hub loads. They do not include this imposed pulse or the inertial rotor reaction.

### 12.4 Acceptance criteria

The four default FSI checks are:

| Check | Implemented measure | Limit |
|---|---|---:|
| Geometry | Maximum componentwise aerodynamic-grid change divided by reference chord | `1e-7` |
| Attachment velocity | Change of base velocities scaled by `max(V,1)`, angular velocities and relative rates by `max(abs(Omega),1)` | `1e-6` |
| Loads | Infinity-norm change in concatenated forces and moments/chord, normalized by the larger current/previous norm or 1 | `1e-5` |
| Final equilibrium | Infinity norm of the full native-plus-component residual at the final aerodynamic trial | `1e-6` |

The geometry check covers all surfaces; the additional velocity check samples rotor attachment motion and relative rates. The load check includes all nodal loads and reduced generalized moments. The first load iteration has infinite change, so acceptance requires a load comparison with an earlier trial.

After provisional convergence, the final aerodynamic recomputation must independently satisfy the load tolerance and final equilibrium tolerance. An inner frozen-load Newton success alone is insufficient. The default relaxation of the entire mixed state is 0.6, with at most 30 coupling iterations.

The residual numbers are the implementation's scaled numerical measures. They should not be read directly as a force error in N or an angle error in radians.

### 12.5 Wake transactions and failure handling

The existing UVLM facade has accepted, trial and commit buffers. Repeated `evaluate_trial!` calls at one target time start from accepted aerodynamic history; FSI iterations do not each advance a new physical wake row. The final accepted trial is committed exactly once.

On failure, the catch block calls `rollback_time_step!`, restores native beam states/rates, unknowns and clocks, restores support states and nodal history, and raises the error. The driver does not silently accept a failed step or automatically loosen the tolerances.

A completed problem cannot be run again in place. A partially advanced response cannot be restarted at zero through this API. Build a fresh problem for another trajectory; persistent checkpoint/restart is not implemented.

## 13. Initialization and developed wakes

### 13.1 Two available modes

`initialization.mode=:unloaded` starts from the unloaded structural state and initial aerodynamic storage. It is useful for startup integration checks, but its response contains mean-load application and wake-development transients.

`initialization.mode=:mean_equilibrium` computes a mean-loaded nonlinear configuration and develops a periodic fixed-shape aerodynamic wake. It retains all aerodynamic loads during the subsequent response.

The second mode is a mean equilibrium, not a fully coupled periodic structural orbit. Periodic blade loading remains after initialization and can drive structural settling.

### 13.2 Fixed-shape aerodynamic settling

For a fixed structural shape and relative support state, `_settle_fixed_shape` constructs a fresh aerodynamic problem and rotates the blades over complete revolutions. The time step must divide a full revolution to within the code tolerance. This is checked explicitly.

It averages the mapped loads over each revolution and compares successive cycles using:

- Revolution-mean mapped load change.
- Mapped load change at the same endpoint phase.
- Bound-circulation change at the same endpoint phase.

The maximum of the three is the periodic residual. The default requires at least four revolutions, allows at most twenty, and uses `1e-3` tolerance. It does not directly compare every retained wake coordinate or every phase sample of the entire cycle; it is the implemented periodicity guard, not a proof of a periodic wake orbit in every state variable.

If the guard fails, initialization raises an error with the separate residual components. It does not accept a nonperiodic trim or automatically replace it with an unloaded state.

### 13.3 Nonlinear mean equilibrium

`initialize_mean_equilibrium!` temporarily sets native structural `Delta_t=Inf` and zeroes rate terms. It alternates:

1. Periodic fixed-shape aerodynamic settling at the current shape.
2. A static augmented mechanical solve driven by the revolution-mean loads.
3. A relaxed update of the mixed structural/support state.

The outer convergence checks are normalized mixed-state change and mean-load change, both defaulting to `1e-5`, with 15 maximum iterations and relaxation 0.6. The state norm includes the mixed unknown vector, not only physical wingtip displacement.

Once converged, the driver sets relative rates to zero and settles a final wake at the accepted shape. The final accepted aerodynamic system is copied into a new dynamic problem. `phase_time_offset` preserves the rotor azimuth corresponding to that developed wake, while the reported response time starts again at zero.

At a static carrier shape, the support/wing gyro velocity terms vanish. Spin is still present in blade motion and the developed wake. Gyroscopic mechanical coupling becomes active when the carrier or relative support moves during the response.

### 13.4 What the initialized response means

A mean-loaded starting configuration avoids applying the entire mean aerodynamic load suddenly. It does not remove periodic aerodynamic forcing. For damping identification, an unperturbed dynamically settled baseline and a phase-matched perturbed trajectory are still needed. A short no-pulse response can contain periodic forcing or residual settling even if the static mean iterations converged.

## 14. Configuration and running cases

### 14.1 Where to edit settings

The intended user-facing configuration is [chang_aerobeams_case.jl](chang_aerobeams_case.jl). It merges selected overrides into `chang_aerobeams_defaults()`. Julia named tuples are merged at each nested level; replacing an entire group with only one field would drop its other required fields.

| Group | Controls |
|---|---|
| `reference.wing` | Wing planform, aerodynamic panel counts, symmetry and inherited structural discretization |
| `reference.propeller` | Radius, chord, blades, blade mesh, RPM, requested attachment fractions and collective offset |
| `reference.structural` | Calibrated reduced-support inputs, pylon/blade mass data and damping ratios |
| `reference.simulation` | Air density, speed, incidence, sideslip, azimuth step, duration and aerodynamic interaction |
| `reference.aerodynamic` | Fixed core radius, elastic-axis fraction and reduced modal arm |
| `reference.wake` | Wing/blade wake capacities |
| `solver` | Support choice, time step, speed policy, Newton/FSI controls, initial relative state and pulse |
| `initialization` | Unloaded or mean-loaded start and settling/convergence controls |
| `pylon` | Beam-pylon sections and rotor-only tip properties |
| `output` | New example output path, label, plotting and native history storage |

The main builder currently uses the inherited wing element count, normally set by `reference.wing.spanwise_panels`. Lower-level wing/interface helpers support independent span meshes, and verification exercises them, but the public case does not yet expose a separate structural wing-count field.

The tabulated Chang wing mass and stiffness distributions come from the inherited model helpers linked in Section 5. They are not editable arrays exposed through `reference.structural` in this case interface.

### 14.2 Principal new defaults

| Setting | Default |
|---|---|
| Freestream speed / duration | 65 m/s / 0.25 s |
| Propeller station | Requested 65% span |
| Support | `:reduced` |
| Aerodynamic interaction / wing symmetry | Enabled / disabled |
| Near-field force / relative moment map | `:imperial` / exact virtual work in the new adapter |
| Time step | Inherited azimuth-based value |
| Speed policy | `:fixed_rpm` |
| Newton iterations / tolerance | 25 / `1e-8` |
| FSI iterations / relaxation | 30 / 0.6 |
| Geometry / velocity / load / final equilibrium tolerance | `1e-7` / `1e-6` / `1e-5` / `1e-6` |
| Relative angle guard | 0.35 rad |
| Initial relative coordinates/rates | Zero unless supplied |
| Pulse | Rotor 1 pitch, 2 N m, start 0.15 s, duration 0.02 s |
| Initialization | `:mean_equilibrium` |
| Mean iterations / tolerance / relaxation | 15 / `1e-5` / 0.6 |
| Aerodynamic settling revolutions / tolerance | Minimum 4, maximum 20 / `1e-3` |
| Pylon elements | 4 per pylon |
| Write results / plot / native history | Enabled / enabled / disabled |

These defaults specify a case to investigate; they are not a certification of mesh, wake or trim convergence. In particular, the stricter wake-periodicity check has rejected one refined setup, as documented in Section 16.

### 14.3 Commands

From the repository root:

```powershell
# Nonlinear wing with two-coordinate supports
julia --compiled-modules=existing --project=. lib/WingPropellerUVLM/examples/chang_aerobeams_coupled/run_chang_aerobeams.jl

# Nonlinear wing with beam pylons and physical tip rotors
julia --compiled-modules=existing --project=. lib/WingPropellerUVLM/examples/chang_aerobeams_coupled/run_chang_beam_pylon.jl
```

`--compiled-modules=existing` was used in this restricted verification environment to avoid writing new compiled caches outside the workspace. A normal writable Julia setup can omit it. No new dependency installation is performed by these runners.

The reduced runner uses the case as configured. The beam-pylon runner overrides the support choice and output label. Both call `run_chang_aerobeams`.

### 14.4 Programmatic use

```julia
include("lib/WingPropellerUVLM/examples/chang_aerobeams_coupled/src/ChangAeroBeams.jl")
using .ChangAeroBeams

cfg = chang_aerobeams_defaults()
cfg = merge(cfg, (
    solver=merge(cfg.solver, (support=:beam_pylon,)),
    output=merge(cfg.output, (label="my_beam_pylon_case",)),
))

problem = build_chang_aerobeams(cfg)
solve_chang_aerobeams!(problem)
paths = save_chang_aerobeams(problem)
```

The builder initializes native structural storage but does not run the coupled response. This separation permits inspection before solving. `run_chang_aerobeams(cfg)` combines build, solve and optional save.

To request an unloaded startup explicitly, merge `initialization=(... mode=:unloaded ...)` into the existing group. Do not interpret this mode as an accepted replacement for a failed periodic trim when reporting a loaded stability result.

## 15. Results and file formats

### 15.1 In-memory histories

`coupled.results` stores:

| Key | Content |
|---|---|
| `:time` | Response time starting at zero |
| `:node_positions` | All structural node positions, `3×Nnodes×Ntimes` |
| `:eta`, `:eta_velocity`, `:eta_acceleration` | Relative support history, empty first dimension for beam pylons |
| `:rotor_angles` | Two displayed shaft/support angles per rotor |
| `:tip_twist_rad` | Wingtip chord inclination from its finite deformation matrix |
| `:hub_forces`, `:hub_moments` | Aerodynamic resultants at physical rotor hubs |
| `:coupling_iterations` | FSI iterations at accepted steps |
| `:geometry_residual`, `:velocity_residual`, `:load_residual`, `:mechanical_residual` | Acceptance diagnostics |
| `:initialization`, `:mean_equilibrium` | Mean-init records/state/loads when that mode is used |
| `:paths` | Saved response paths when the combined run API saves outputs |

### 15.2 Angle definitions require care

For reduced supports, displayed rotor pitch/yaw are the relative coordinates `theta,psi`, measured with respect to the moving wing carrier. For beam pylons, the code takes `s=D_h*e_x` and reports

$$
\theta_{out}=\operatorname{atan2}(-s_z,s_x),\qquad
\psi_{out}=\operatorname{atan2}(-s_y,s_x).
$$

These are absolute shaft-direction angles in `U`. They are not identical coordinate definitions to the reduced output, and at finite combined rotation they should not be assumed to be the original sequential rotation parameters. A physical cross-model comparison should calculate the same absolute shaft vector in both models.

The output called “wingtip twist” is computed from the deformed tip chord `d=D_tip*e_x` as `atan2(d_z,d_x)`. It is a geometric chord inclination, not a raw WM component or a separately extracted local torsional strain. The current plots retain their short existing output label.

### 15.3 Saved artifacts

`save_chang_aerobeams` requires a completed run and writes:

- `<label>_history.csv`: time, physical wingtip position, chord inclination in degrees, FSI count, four residuals, and per-rotor displayed angles plus three aerodynamic hub force components.
- `<label>_settings.txt`: resolved reference, actual solver/initialization settings, wing/pylon metadata, phase offset, time step and maximum residuals.
- `<label>_initialization.csv`: mean-iteration history, if applicable.
- `<label>_history.png`: a four-panel plot of wingtip height, chord inclination, pitch and yaw, if plotting is enabled.

Hub moments and relative rates remain in memory and are not included in the current history CSV. `save_structural_history=true` additionally invokes native AeroBeams time-step storage; it does not create a persistent checkpoint/restart facility.

The current settings file does not serialize the entire `config.pylon` group or each `TipRotor`'s resolved mass, ordinary tensor and independent spin inertia. In particular, nondefault rotor/nacelle overrides cannot all be reconstructed from that file alone. Retain the input case/configuration alongside the outputs when using those overrides.

Some verification settings files were generated before the explanatory “actual integration” header was expanded. Their embedded legacy groups still do not control these simulations. The operative solver settings and this report identify the algorithms that were run.

## 16. Verification evidence

### 16.1 What was executed

The existing [verification record](validation/verification_report.md) and logs contain the actual executed results. This report summarizes those runs; simulations were not rerun merely to produce the documentation.

| Verification group | Result | Evidence |
|---|---:|---|
| Assembled component/interface/coupled suite | 119/119 passed | [verification.log](validation/verification.log) |
| Dry-wing data/modal comparison | 9/9 passed | [modal_comparison.log](validation/modal_comparison.log) |
| Two longer perturbed startup responses | 14/14 passed | [perturbed_responses.log](validation/perturbed_responses.log) |
| Stricter longer-wake periodic trim | Rejected for nonconvergence | [longer_wake_trim_rejection.log](validation/longer_wake_trim_rejection.log) |

The focused log includes third-party ImplicitAD method warnings. PowerShell's combined stderr formatting labels those warning records as native errors; the Julia Test summary reports no failed or errored tests. The later two-response process explicitly returned Julia exit code zero.

The corresponding reproduction commands, from the repository root, are:

```powershell
julia --compiled-modules=existing --project=. lib/WingPropellerUVLM/examples/chang_aerobeams_coupled/validation/runtests.jl
julia --compiled-modules=existing --project=. lib/WingPropellerUVLM/examples/chang_aerobeams_coupled/validation/compare_wing_modes.jl
julia --compiled-modules=existing --project=. lib/WingPropellerUVLM/examples/chang_aerobeams_coupled/validation/run_verification_responses.jl
```

The suite includes the coarse mean-loaded check. The last command reproduces the explicitly unloaded longer response checks; it does not convert the rejected longer-wake trim into a passing trim.

### 16.2 Scope of the focused checks

| Verification source | Property checked |
|---|---|
| [test_reduced_rotor.jl](test_reduced_rotor.jl) | Complete Chang mass/gyro blocks, finite-difference residual tangents, skewness, spin reversal, finite-motion invariance and power, dry whirl frequencies |
| [test_tip_rotor.jl](test_tip_rotor.jl) | Rotor/pylon mass ownership, valid inertia, ordinary/gyro tangents, finite angular-momentum derivative, power and objectivity |
| [test_aerodynamic_interface.jl](validation/test_aerodynamic_interface.jl) | Rigid Chang grid parity, finite-rotation virtual work, vortex force offsets, wrench balance, independent meshes, IDs and aerodynamic rollback |
| [test_coupled_solver.jl](validation/test_coupled_solver.jl) | Independent directional derivatives of the augmented residual, interacting response, known growth fit and failed-step restore |
| [check_pylon.jl](validation/check_pylon.jl) | Connected nonlinear pylon and converged physical tip-rotor response |
| [check_mean_equilibrium.jl](validation/check_mean_equilibrium.jl) | Loaded nonlinear mean equilibrium, developed-wake response and output |

These checks use physics identities, reference blocks and independent differences. They establish more than code execution alone, but they do not replace aerodynamic refinement or experimental validation.

### 16.3 Dry-wing mass and modal comparison

The reconstructed wing mass is approximately **177.3285166847 kg**. The integrated local spatial-block conservation error is at roundoff level in the executed cases.

The first six ordered dry-wing frequencies at 32 elements were:

| Ordered mode | Existing Chang [Hz] | AeroBeams [Hz] | Relative difference |
|---|---:|---:|---:|
| 1 | 6.63652 | 6.59371 | −0.6450% |
| 2 | 20.81284 | 20.80096 | −0.0571% |
| 3 | 36.51431 | 36.42610 | −0.2416% |
| 4 | 42.90502 | 42.80883 | −0.2242% |
| 5 | 93.61634 | 93.96137 | +0.3686% |
| 6 | 112.71452 | 113.08562 | +0.3292% |

The maximum absolute difference among those six ordered frequencies was 16.16% at 8 elements, 3.31% at 16 elements and 0.645% at 32 elements. Coarse high-mode discrepancies must therefore not be hidden by the close first frequency.

This comparison supports the section/property mapping and demonstrates mesh sensitivity. It compares ordered frequencies, not a modal assurance criterion or a proof of identical eigenvectors. It excludes installed rotors and UVLM wake dynamics. Full data are in [wing_modal_comparison.csv](output/wing_modal_comparison.csv).

### 16.4 Coarse mean-loaded interacting verification

The mean-equilibrium integration check used the actual Chang distributions with:

- Four wing beam elements and a `4×2` wing aerodynamic panel grid.
- Four blades, `2×2` panels per blade.
- 65 m/s, 3° incidence and an attachment at 50% span.
- Wing/propeller interaction enabled.
- Twelve wake rows per surface.
- 5° azimuth step, `Delta_t=0.00067605956 s`.
- At least three and at most ten settling revolutions, periodic tolerance `0.002`, and mean relaxation 1.

This deliberately coarse setup exercises the coupling algorithm; it is not an aerodynamically converged production case.

Mean equilibrium converged in seven outer iterations:

| Quantity | Value |
|---|---:|
| Initial loaded wingtip `z` | 0.0301395864 m |
| Initial chord inclination | −0.0434118° |
| Relative support pitch / yaw | +0.639145° / −0.358952° |
| Final normalized mixed-state change | `3.2153e-7` |
| Final mean-load change | `3.9609e-6` |
| Final fixed-shape periodic residual | `2.5126e-5` |
| Preserved aerodynamic phase-time offset | 0.146028865 s |

The subsequent no-pulse response completed 14 steps to 0.00946483 s, with maximum diagnostics:

| Measure | Observed maximum | Limit |
|---|---:|---:|
| Geometry | `4.8726e-9` | `1e-7` |
| Attachment velocity | `1.1631e-7` | `1e-6` |
| Loads | `9.1544e-6` | `1e-5` |
| Coupled equilibrium | `3.0339e-8` | `1e-6` |

![Coarse mean-equilibrium response: wingtip displacement, chord inclination and relative support angles.](output/verification_mean_equilibrium_history.png)

Data: [history](output/verification_mean_equilibrium_history.csv), [initialization](output/verification_mean_equilibrium_initialization.csv), [settings](output/verification_mean_equilibrium_settings.txt).

### 16.5 Longer-wake trim rejection

A stricter follow-up used eight wing elements, an `8×3` wing aerodynamic grid, `3×2` panels per blade, 36 wing wake rows and 72 blade wake rows. At 5° azimuth increments, 72 rows correspond to one revolution of blade wake.

The fixed-shape aerodynamic periodic residual was **0.0088319** after the allowed twelve revolutions, exceeding the configured `0.001` limit. Initialization raised an error; that case did not produce an accepted periodic trim.

This is an unresolved convergence result. It is not evidence of a validated trim, and it does not by itself prove a physical instability. The saved older error reports the maximum periodic measure rather than its individual components; the current error message reports mean-load, endpoint-load and circulation errors separately to assist further diagnosis.

### 16.6 Longer perturbed startup responses in both models

The two formulations subsequently completed unloaded-start transients with the longer-wake discretization. Each accepted **369 steps**, ending at **0.24946598 s**. A 2 N m pulse acted from 0.05 to 0.07 s. These runs used `mode=:unloaded`; they are distinct from the accepted coarse mean-loaded case and the rejected longer-wake trim.

| Quantity | Reduced support | Beam pylon |
|---|---:|---:|
| Initial wingtip `z` | 0 m | 0 m |
| Maximum wingtip `z` | 0.06421274 m | 0.06284696 m |
| Maximum absolute displayed pitch | 1.34034° | 1.42152° |
| Maximum absolute displayed yaw | 1.03021° | 1.19045° |
| Maximum FSI iterations | 13 | 14 |
| Maximum load residual | `9.9964e-6` | `8.2081e-6` |
| Maximum coupled equilibrium residual | `1.1793e-7` | `9.9876e-7` |

All 14 additional checks passed. The similar displacement scale demonstrates that both routes execute an interacting response; it is not a proof of equivalent pylon dynamics. In particular, displayed angle coordinates differ between the models as explained in Section 15.

![Reduced-support unloaded-start interacting response.](output/verification_reduced_history.png)

![Beam-pylon unloaded-start interacting response.](output/verification_beam_pylon_history.png)

Data: [reduced history](output/verification_reduced_history.csv), [beam-pylon history](output/verification_beam_pylon_history.csv), [numerical summary](output/perturbed_response_summary.csv).

The startup loading, developing wake and imposed pulse all contribute to these short histories. Their peaks must not be interpreted as a flutter boundary or a growth-rate measurement about a converged periodic state.

## 17. Stability-study support

### 17.1 Why paired trajectories are supplied

A rotating propeller creates periodic forcing even without an instability. A response envelope from a single loaded trajectory can therefore mix forcing, wake startup, mean-equilibrium error and perturbation growth.

`run_chang_aerobeams_speed_sweep` requires mean-equilibrium initialization. At every speed, it independently runs a zero-pulse baseline and a trajectory with the configured pulse. It checks that their initialization phase offsets agree to `1e-12`, then forms

$$
\delta y(t)=y_{perturbed}(t)-y_{baseline}(t).
$$

Both runs share the operating point and numerical settings. The code preserves the baseline/perturbed results in memory for independent checks instead of returning only a fitted number.

### 17.2 Growth estimator

`perturbation_growth_rate` computes nonoverlapping RMS windows after a selected start time,

$$
A_k=\sqrt{\frac1{N_k}\sum_{j\in k}\delta y_j^2},
$$

and fits

$$
\log A_k=a+\sigma t_k.
$$

It requires at least three samples per retained window and at least four resolved, non-negligible windows. The result contains `growth_rate_per_s`, `r_squared`, window centers and RMS amplitudes. A known decaying signal is used in verification.

The estimate is an observed signal growth exponent. It is not automatically a modal damping ratio because no modal frequency is extracted. Multiple participating modes, periodic coefficients, beating, fitting-window choice and nonlinear saturation can change the fitted value.

### 17.3 Supplied speed runner and current status

`run_speed_sweep.jl` requests speeds 50:10:90 m/s, five seconds per response, a pulse beginning at one second, and a fit starting at two seconds with 0.2 s RMS windows. That full campaign has **not** been executed as part of the recorded verification.

The implementation is not a Floquet eigenanalysis, periodic-orbit shooting method, automatic flutter root search or nonlinear branch continuation. A speed where the fitted exponent changes sign is a candidate requiring duration, amplitude, phase and discretization checks.

## 18. Assumptions and current limits

| Area | Current assumption or limitation | Consequence |
|---|---|---|
| Wing | Straight constant-`x` elastic axis; geometric nonlinearity with elastic sectional properties | Swept/curved reference axes need an explicit geometry extension |
| Wing discretization | Distributed reconstruction and finite high shear stiffness | Refine and compare dry modes; do not assume identical legacy matrices |
| Reduced support | Linear CM modes and frozen intrinsic modal inertia | Suitable initially for small relative deformation, despite finite wing rotation |
| Beam pylon | Default bending calibration matches pure tip-moment compliance only | Measured sections and independent modal calibration are needed |
| Rotor inertia | Default axisymmetric effective tensor and calibrated polar inertia | Measured rotor/nacelle data are needed for a predictive installation |
| Hub/CM | Coincident for `TipRotor` | A measured eccentric CM requires an additional spatial-inertia extension |
| Rotor speed | Prescribed constant spin under the selected speed policy | No motor, governor, spin-up acceleration, windmilling or shaft-speed DOF |
| Blades | Rigid aerodynamic surfaces; effective rotor structural inertia | No elastic blade modes or centrifugal blade-stiffness model |
| Pylon aerodynamics | `aeroSurface=nothing` and no explicit UVLM pylon surface | Pylon/nacelle aerodynamic forces are not included |
| Damping | Relative viscous damping in reduced model; no beam material damping | Beam-pylon damping must be added and calibrated if required |
| Gravity | Zero in the current coupled examples | Loaded sag under weight is not included, despite component gravity APIs |
| Aerodynamics | Inherited UVLM and corrected inviscid force formulation | No new viscous, separation or stall model |
| Wake | Finite retained wakes and code-defined periodic guard | Wake length and phase convergence remain necessary |
| Initialization | Mean static structure plus periodic fixed-shape aerodynamic wake | A fully coupled periodic baseline still needs dynamic settling |
| Time discretization | Mixed native trapezoidal and added Newmark reconstruction | Time-step refinement is required; legacy integrator settings do not apply |
| Mechanical solver | Dense augmented Jacobian, finite differences of component corrections | Larger systems need performance work |
| Stability | Paired time-response growth fit | No full coupled eigenvalue/Floquet solver or nonlinear continuation |
| Lifecycle | Fresh problem per run | No persistent checkpoint/restart |
| Verification | Component identities and short responses | No experimental validation or demonstrated flutter boundary yet |

## 19. Remaining work for a predictive study

The code now supplies both requested structural routes. A credible Chang wing–propeller stability comparison should proceed through the following stages.

1. **Fix the physical comparison definition.** Record whether speed sweeps hold RPM or advance ratio fixed; confirm actual attachment nodes, pitch setting, rotor mass, pylon geometry, and shaft direction. Decide which absolute shaft angles will be compared between formulations.
2. **Complete structural calibration.** Refine dry-wing modes and compare mode shapes. Compare fixed-base reduced whirl frequencies with an independent solution. For beam pylons, supply physical section properties and rotor/nacelle inertia, then compare installed dry modes and tip compliance. Do not tune only a single static rotation and call the models dynamically equivalent.
3. **Resolve aerodynamic periodicity.** At fixed structure, study azimuth step, wing/blade panels, core radius, wake length and revolution count. Use the separately reported mean, endpoint and circulation errors to diagnose the rejected trim. Retain strict rejection when the chosen periodicity criteria fail.
4. **Converge the loaded state.** Repeat nonlinear mean equilibrium with refined settings, then run a zero-pulse dynamic response long enough to characterize residual structural settling and periodic forcing. Establish repeatable shaft phase and loading.
5. **Check perturbation linearity.** Run phase-matched baseline/perturbed pairs at several sufficiently small pulse amplitudes. Check whether differences scale with amplitude before interpreting growth rates as local stability information.
6. **Perform the speed study.** Use consistent physical assumptions and identical comparison signals. Inspect several fit windows, RMS-window lengths, frequencies and response durations. Refine around candidate sign changes in growth.
7. **Establish numerical convergence.** Repeat candidate boundaries with finer time steps, panels, beam meshes and wakes, and tighter FSI/Newton tolerances where needed. Record shifts in the inferred boundary separately from fit uncertainty.
8. **Investigate nonlinear response.** Beyond a verified local boundary, increase run duration and examine whether oscillations decay, grow without bound in the valid model range, or saturate. Use a periodic-orbit/continuation method if a bifurcation branch is required. Respect the reduced support angle guard or use the calibrated beam-pylon model.
9. **Validate against an external reference.** Compare the same physical installation and operating convention with published/experimental Chang data or an independent model. The present report has not performed that external validation.

These are remaining analysis tasks, not completed results. The current short response demonstrations should be retained as integration checks alongside the more demanding campaign.

## 20. Function reference and traceability

The following locations identify the main implementation blocks in the source reviewed for this report. Line numbers are navigation aids and may change with later development.

| Source/function | Location | Role |
|---|---|---|
| Module entry point | [ChangAeroBeams.jl](src/ChangAeroBeams.jl#L1) | Package imports, legacy reference include and new module API |
| `_wing_element_inertias` | [AeroBeamsWing.jl](src/AeroBeamsWing.jl#L25) | Integrates reconstructed full local inertia densities |
| `_wing_element_stiffness` | [AeroBeamsWing.jl](src/AeroBeamsWing.jl#L63) | Integrates compliance with the physical curvature convention |
| `build_wing` | [AeroBeamsWing.jl](src/AeroBeamsWing.jl#L105) | Creates the nonlinear wing, clamp, load BCs and dynamic problem |
| `wing_node_parameters` / `wing_node_poses` | [AeroBeamsWing.jl](src/AeroBeamsWing.jl#L181) | Reads current mixed unknowns and computes physical poses |
| `set_wing_loads!` | [AeroBeamsWing.jl](src/AeroBeamsWing.jl#L241) | Updates physical nodal force/moment boundary conditions |
| `wm_spatial_angular_kinematics` | [AeroBeamsWing.jl](src/AeroBeamsWing.jl#L264) | Converts WM rates to physical angular velocity/acceleration |
| `wing_node_motion` | [AeroBeamsWing.jl](src/AeroBeamsWing.jl#L282) | Reconstructs trial node motion from accepted history |
| `ReducedRotor` constructor | [ReducedRotor.jl](src/ReducedRotor.jl#L42) | Builds the individual reference installed component |
| `rotor_mass_data(::ReducedRotor)` | [ReducedRotor.jl](src/ReducedRotor.jl#L81) | Objective moving-CM mass and shape derivatives |
| `rotor_spin_data(::ReducedRotor)` | [ReducedRotor.jl](src/ReducedRotor.jl#L97) | Exact relative spin momentum and all gyro blocks |
| `rotor_mechanics(::ReducedRotor)` | [ReducedRotor.jl](src/ReducedRotor.jl#L146) | Moving-carrier wrench and relative residual |
| `build_wing_with_beam_pylons` | [BeamPylon.jl](src/BeamPylon.jl#L8) | Pylon topology, section assumptions and global load maps |
| `pylon_hub_poses` | [BeamPylon.jl](src/BeamPylon.jl#L68) | Physical tips and incremental rotor orientations |
| `TipRotor` constructor | [TipRotor.jl](src/TipRotor.jl#L36) | Validates rotor-only mass, tensor and spin inertia |
| `rotor_mechanics(::TipRotor)` | [TipRotor.jl](src/TipRotor.jl#L113) | Ordinary and gyro tip wrench |
| `_ca_wing_geometry` | [AerodynamicInterface.jl](src/AerodynamicInterface.jl#L60) | Structural-node chord interpolation and vortex fractions |
| `coupled_grids` | [AerodynamicInterface.jl](src/AerodynamicInterface.jl#L112) | Reduced-support blade/wing geometry |
| `aerodynamic_loads` | [AerodynamicInterface.jl](src/AerodynamicInterface.jl#L184) | Reduced generalized loads and aerodynamic hub resultants |
| `_ca_build_aerodynamic` | [AerodynamicInterface.jl](src/AerodynamicInterface.jl#L240) | UVLM objects, local symmetry/ID compatibility handling |
| `pylon_coupled_grids` | [AerodynamicInterface.jl](src/AerodynamicInterface.jl#L293) | Actual beam-tip blade geometry |
| `pylon_aerodynamic_loads` | [AerodynamicInterface.jl](src/AerodynamicInterface.jl#L335) | Physical hub loads once at each pylon tip |
| `chang_aerobeams_defaults` | [CoupledSolver.jl](src/CoupledSolver.jl#L7) | New configuration groups and shared defaults |
| `build_chang_aerobeams` | [CoupledSolver.jl](src/CoupledSolver.jl#L61) | Copies/resolves inputs and constructs all subsystems |
| `_component_residual` | [CoupledSolver.jl](src/CoupledSolver.jl#L204) | Adds installed support/rotor mechanics at mapped nodes |
| `coupled_mechanical_arrays!` | [CoupledSolver.jl](src/CoupledSolver.jl#L242) | Native residual/Jacobian plus small component correction |
| `_solve_mechanics!` | [CoupledSolver.jl](src/CoupledSolver.jl#L272) | Newton solve and line search |
| `_settle_fixed_shape` | [CoupledSolver.jl](src/CoupledSolver.jl#L356) | Revolution averaging and aerodynamic periodicity guard |
| `initialize_mean_equilibrium!` | [CoupledSolver.jl](src/CoupledSolver.jl#L414) | Nonlinear mean loading and developed wake transfer |
| `_record!` | [CoupledSolver.jl](src/CoupledSolver.jl#L472) | Physical histories, angle definitions and diagnostics |
| `solve_chang_aerobeams!` | [CoupledSolver.jl](src/CoupledSolver.jl#L510) | Accepted-step loop, FSI iteration, commit and rollback |
| `save_chang_aerobeams` | [CoupledSolver.jl](src/CoupledSolver.jl#L650) | CSV, settings, initialization records and plots |
| `run_chang_aerobeams` | [CoupledSolver.jl](src/CoupledSolver.jl#L716) | Build/solve/save convenience API |
| `perturbation_growth_rate` | [Studies.jl](src/Studies.jl#L7) | RMS windows and logarithmic growth fit |
| `run_chang_aerobeams_speed_sweep` | [Studies.jl](src/Studies.jl#L38) | Phase-matched baseline/perturbed speed campaign |

The report is based on the local implementation and recorded verification artifacts linked above. It introduces no external experimental claims and no changes to existing analysis code.
