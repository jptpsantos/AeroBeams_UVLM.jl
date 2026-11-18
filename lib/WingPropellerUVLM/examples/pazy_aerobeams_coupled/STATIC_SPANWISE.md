# Static Pazy spanwise normal force and pitching moment

Run `run_pazy_wing_uvlm_static_sweep.jl` from the IDE, or from the repository root:

```powershell
julia --project=lib/WingPropellerUVLM/examples lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/run_pazy_wing_uvlm_static_sweep.jl
```

Edit these inputs at the top of the runner:

```julia
velocity_range = collect(10.0:10.0:30.0) # m/s, strictly positive
angles_of_attack_deg = [1e-6,3.0]        # first case approximates zero incidence
sweep_angles_deg = [0.0]                # available native data: 0, 10, 20 deg
slope_modes = (:local_ratio,)
```

Every velocity, incidence and sweep combination is solved independently.
The prescribed freestream incidence stays constant within that solve; no
perturbed angles or additional equilibrium solves are used by the default mode.
The script uses the existing `run_pazy_wing_uvlm_static`: nonlinear AeroBeams
`SteadyProblem` balance with externally applied UVLM loads. There is no
structural inertia or dynamic beam settling. At each static trial geometry,
a fresh UVLM free wake develops while the wing is held fixed. The method
requires both wake convergence and coupled load/geometry convergence.

The native geometry, structural matrices, rod options and coordinate mapping
match `run_pazy_wing_uvlm.jl`. UVLM x points downstream, y spans the wing and z
points up. The leading edge is upstream of the trailing edge, and the physical
chord is streamwise even with sweep. Positive native sweep is aft sweep; the
available native structural data do not support arbitrary or negative sweep.
Root pitch and gravity are zero. Incidence is applied once through the
freestream. `symmetric_wing=true` includes an image wing/wake at y=0;
`false` represents an isolated cantilever. Reported forces concern only the
explicitly modelled semispan.

## Normal-force coefficient

Each UVLM spanwise node produces

```text
q_inf = rho * U^2 / 2
cn_j = dot(F_j, n_j) / (q_inf * chord * dy0_j)
alpha_local_j = atan(dot(V_inf, n_j), dot(V_inf, t_chord_j))
cn_alpha_local_j = cn_j / alpha_local_j
```

`F_j` sums the dimensional vortex-vertex forces along the chord at station j.
The existing `imperial_nodal_forces` redistribution assigns half of each
segment force to each endpoint and one quarter of each residual unsteady
panel force to each corner, exactly as for structural load transfer. The
section force vectors conserve the full UVLM resultant. `n_j` is the unit
normal formed from the physical chord at that node and the centered span
tangent of the physical chord midpoints, with one-sided tangents at the root
and tip. It points up on the undeformed wing and follows the deformed shape.

`dy0_j` is the undeformed projected tributary width: a full panel width
`L*cos(sweep)/spanwise_panels` at interior nodes, and half that width at the
root and tip. These tributary areas sum to the reference semispan area; the
existing result key `strip_area_m2` holds these nodal tributary areas. The
fixed normalization supports comparison at the same material stations across
deformations. `eta = s/L` labels the spanwise nodes along the undeformed spar.
Twenty spanwise panels therefore produce 21 samples at `eta = 0:0.05:1`.
`s_reference_m` and `y_reference_m` give the undeformed spar distance and
projected span coordinate; separate columns give the physical chord midpoint
at each deformed node. Both plots use normalized span from 0 (root) to 1
(wingtip), with ticks every 0.1. Endpoint coefficients are averages over their
half-width tributary strips, not pointwise limits; the tip value need not be
zero. The vertical axes show sectional cn and cn_alpha, respectively.

The default `:local_ratio` mode uses the **local geometric incidence of the
deformed wing**. `t_chord_j` is the downstream unit vector along the physical
chord at node j; `n_j` is its deformed surface normal. The two-argument
`atan` gives signed incidence in that chord/normal plane, in radians. This
includes elastic twist and bending, so local incidence generally varies along
the wing even though the prescribed freestream angle stays fixed. For a flat,
undeformed wing with its physical chord streamwise, it equals the prescribed
angle, including for the native swept configurations.

The angle definition uses the freestream, excluding induced downwash. UVLM's
induced-flow effects remain in cn. Dividing by an incidence that includes
downwash would produce a different coefficient and is not this convention.
The coefficient continues to use freestream dynamic pressure and undeformed
projected tributary area; it is an apparent finite-wing coefficient, not a direct
replacement for a 2D airfoil input in a strip model.

`cn/alpha_local` is a **secant approximation**, not the tangent derivative
`dcn/dalpha`. Interpreting it as a slope assumes an approximately linear,
zero-intercept relation near the operating point. Division uses the local
angle in radians, so ratios have units **1/rad**. Multiply by `pi/180` to
convert to per degree. The cutoff for `abs(alpha_local)` is configurable
with `local_alpha_tolerance`. The runner uses `1e-10` rad for the
small-incidence case. Below or at the cutoff, or if the sectional flow
projection vanishes, the script writes NaN and sets `cn_alpha_local_valid=false`.

### Proxy for zero incidence

At exactly 0° on the unloaded symmetric wing, the local incidence and
normal force vanish, so the local ratio is 0/0. The runner uses a single
equilibrium at `1e-6` degrees (about `1.745e-8` rad) as the proxy for
zero. The tiny nonzero load is divided by the resulting local angle at
each station. The actual imposed incidence is recorded in the CSV,
plot and `aoa1p0em6deg` folder. The pitching moment uses the loads from that same
equilibrium and therefore remains small near zero incidence.

The runner uses a `load_residual_force_floor` of `1e-9` N to scale the
static load residual against the tiny nodal loads. With the solver API's
1 N default floor, a small absolute mismatch can pass the coupling
tolerance while the elastic twist is still changing by an appreciable
fraction of the prescribed incidence. The cutoff and residual settings
preserve the use of one equilibrium per flight condition.

Optional derivative modes remain available through the solver API for earlier
studies: `:equilibrium` re-solves the complete flexible-wing equilibrium at
`alpha +/- alpha_step`, and `:fixed_shape` independently settles wakes at those
incidences while holding the central equilibrium shape fixed. Both settle
separate wakes at constant airspeed and vary the freestream direction.
Use `alpha_step_deg` in the runner or `alpha_step` in radians in the solver
API. Both methods also work at exactly zero incidence.

### Perturbation slope with respect to local incidence

The existing derivative columns divide the change in cn by the imposed
incidence change `2 * alpha_step`. The perturbations now also report the
change in geometric local incidence at each material station. The local
perturbation slope in `spanwise_cn_alpha_perturbations.csv` is

```text
delta_cn[j] = cn_plus[j] - cn_minus[j]
delta_alpha_local_rad[j] = alpha_local_plus[j] - alpha_local_minus[j]
cn_alpha_local_perturbation_per_rad[j] = delta_cn[j] / delta_alpha_local_rad[j]
```

The input perturbation changes incidence across the entire wing. The
result describes the sectional response along that imposed change.
A derivative for a single independently perturbed strip would require
a local velocity perturbation or a local change in geometry. A local
normal velocity can enter the boundary-condition right-hand side with
the geometry held fixed, and the solved circulation can change across
the span. Perturbing each strip would yield a matrix of strip-to-strip
sensitivities.

Set `slope_modes = (:local_ratio,:fixed_shape,:equilibrium)` to compare
both perturbation modes with the secant. Use `:fixed_shape` for the
steady aerodynamic sensitivity on the central equilibrium geometry,
and `:equilibrium` to include the change in elastic deformation.

## Pitching-moment coefficient

`spanwise_cm.csv` contains nose-up-positive moment coefficients about each
section's deformed physical quarter chord and about the native spar. The
section pitch axis is the cross product of the local normal and downstream
chord unit vector. The extraction sums `(r_vortex - r_reference) × F_vertex`
over the chord at each station and divides its pitch component by
`q_inf * strip_area_m2 * chord`. The vortex vertices are the positions used
by the UVLM load transfer. The same reference area and chord apply to both
moment coefficients.

Use `case.section.cm_quarter_chord` for the moment about the physical
aerodynamic reference and `case.section.cm_spar` for the pitch component
about the beam attachment point. AeroBeams projects its section pitching
moment onto the beam axis for swept wings. Its attached-flow model already
adds the lever arm from `x_ac` to the spar, so a total spar moment supplied as
`cm₀` would count that arm twice. The two derivative modes also export both
moment derivatives per radian of imposed incidence. On a fixed shape with
section chord equal to the reference chord, `x_ac_over_c = 0.25 -
(dcm_quarter_chord/dalpha) / (dcn/dalpha)` where the normal-force slope is
sufficiently far from zero.
This is an apparent aerodynamic center for the imposed perturbation.

## Using the curves in a flutter model

- Use `cn_alpha_local_per_rad` where `cn_alpha_local_valid=true`. Match
  `sweep_deg`, `alpha_deg`, `U_m_s` and `eta` when selecting data, and
  select the `1e-6` degree profile as the proxy for zero incidence.
- The nominal 5 × 15 S0 run at 1° and 120 m/s reaches 12.03 rad⁻¹ at
  `eta ≈ 0.067` with local incidence 0.85°. The corresponding maxima over
  the S10 and S20 1° sweeps are 5.79 and 7.12 rad⁻¹. UVLM loads at each
  station respond to the circulation and geometry of the whole wing.
  The secant can include strong load redistribution as the equilibrium
  bends and twists. A fixed-shape central difference probes the steady
  response with the structure held fixed. A full-equilibrium derivative
  also includes compliance, which the eigenproblem already represents.
- If calibrating the AeroBeams tip-loss function, fit a **dimensionless
  replacement multiplier**. At low incidence on the straight wing the
  first approximation is the fixed-shape UVLM slope divided by the
  airfoil's actual AeroBeams slope at the same Mach number. Applying the
  raw per-radian slope directly as the multiplier over-scales the circulatory
  normal force and pitching moment. Multiplying a UVLM-calibrated slope by
  the stock tip-loss curve counts finite-span effects twice. For deformed
  and swept cases, compare dimensional loads on a common geometry. The
  models use different velocity projections, sweep factors and load axes.
  Fit or average over AeroBeams' span quadrature using normalized material
  arclength. The UVLM endpoint values represent half strips.
  At 3° and 50 m/s the saved S0 profile has an area-weighted mean cn of
  0.301. Evaluating the stock exponential factor with the same local angles
  and a Prandtl–Glauert corrected `2π` slope gives an estimated 0.337
  before re-trimming. If the reduced loading decreases the equilibrium
  deflection, the rectangular-wing hump mode can move to higher speeds.
  [Pazy workshop comparison of flutter branches](https://elib.dlr.de/202893/1/AIAA_SciTech_2024_LDWG_PazyWing_Paper.pdf)
- Compare static deflection and twist at each speed before comparing the
  eigenvalues. Match the structural mesh, matrix set, rod and root-stiffness
  options, tip-mass configuration, gravity, incidence convention and density.
  The stock S0 eigen example uses an exponential correction and a refined
  structural mesh. The archived UVLM curves used the nominal mesh, and the
  coupled static solver uses zero gravity. For the same upright incidence
  convention in `create_Pazy`, use `pitchAngle=0` and
  `flightDirection=[sin(alpha),cos(alpha),0]` with `alpha` in radians.
  Set `tipLossPitchAngle=alpha` when evaluating the stock correction by
  root incidence. Use a fine eigenvalue speed grid
  and track modes through both ends of each unstable branch. The S0 example
  lists a 3° experimental onset and offset at 49 and 58 m/s, respectively.
  A coarse grid can miss a short unstable interval. The `find_flutter`
  helper also requires adjacent nonzero frequencies and a strict damping
  sign change.
- Static derivatives measure the zero-frequency response of the settled
  wake to uniform incidence. Spanwise cross-coupling, unsteady phase and
  damping for a flexible eigenmode require a response to distributed motion
  at the modal frequency. The Pazy
  workshop study describes stability calculations about the deformed
  equilibrium with frequency-dependent aerodynamics.
  [Pazy Aeroelastic Prediction Workshop study](https://ntrs.nasa.gov/citations/20210024322)

## Outputs and convergence

Default base folder: `output/static_spanwise_cn_cm/`. Each prescribed angle
has a separate subfolder containing only that angle's cases. For example:

```text
output/static_spanwise_cn_cm/
    aoa3deg/
        spanwise_cn_alpha.csv
        spanwise_cn_alpha_perturbations.csv
        spanwise_cm.csv
        status.csv
        settings.txt
        case_1_central_convergence.csv
        spanwise_sweep..._alpha....png
        spanwise_cm_sweep..._alpha....png
    aoa5deg/
        ...
```

Decimal and negative angles use `p` and `m`: 2.5 degrees gives `aoa2p5deg`,
and -2.5 degrees gives `aoam2p5deg`. Velocity and sweep cases at the same
incidence share that incidence's folder; case numbering starts at 1 within
each folder and matches its CSV rows and convergence filenames.

- `spanwise_cn_alpha.csv`: one row per case and spanwise node, including cn,
  `alpha_local_deg`, `alpha_local_rad`, `cn_alpha_local_per_rad` and
  `cn_alpha_local_valid`, velocity, prescribed incidence, sweep, reference
  stations, deformed coordinates and wingtip displacement/twist. Earlier
  derivative columns are retained; they are NaN unless explicitly requested.
  `alpha_step_deg` is NaN when no derivative mode is requested.
- `spanwise_cm.csv` has the same case and station keys, both moment
  coefficients and their optional derivatives. The moment plot shows the
  two references separately.
- `spanwise_cn_alpha_perturbations.csv` has one row for each requested
  perturbation mode and spanwise station. It records `delta_cn`,
  `delta_alpha_local_rad`, their ratio and a validity flag. The ratio is
  NaN if either solve fails or the local incidence change has magnitude
  at or below the cutoff. With `:local_ratio` alone, the file has only a
  header.
- `status.csv`: static convergence flags and reasons. The default mode has one
  central solve per condition. Optional derivative diagnostics are also included
  when requested.
- `case_N_SOLVE_convergence.csv`: static iteration histories for central and
  perturbed flexible-wing solves. Time is cumulative aerodynamic settling time.
- `spanwise_sweep..._alpha....png`: cn and cn/alpha_local curves, grouped by
  incidence and sweep with separate curves for velocity. Only converged data
  are plotted.
- `settings.txt`: numerical controls, structural settings and conventions.

The CSVs are checkpointed after each completed case. Running again replaces
files within the matching angle folders; other angle folders remain available.
Choose a new base directory to preserve another study at the same incidence.
Unconverged cases export NaN cn, cm and ratios. Local ratios are also NaN at
near-zero local incidence or when that mode was not requested. Optional tangent
derivatives are NaN if either perturbed solve fails to converge or their mode
was not requested. Inspect validity and convergence flags to distinguish cases.

The result remains available as `static_sweep`. For example:

```julia
case = static_sweep.cases[1]
case.section.eta
case.section.cn
case.section.cm_quarter_chord
case.section.cm_spar
case.section.alpha_local  # radians
case.cn_alpha_local      # cn / alpha_local, 1/rad
case.cn_alpha_local_valid
case.cn_alpha_local_perturbations.fixed_shape  # when both fixed-shape wakes converge
case.diagnostics
```

Wake durations are specified as convective distances `U*t/chord` and converted
to seconds separately for each speed. `maximum_aerodynamic_chords` limits each
wake settling trial; it is not the physical wing simulation time or total
wall-clock time. The existing static solver also requires at least two retained
wake ages. Increase the maximum if `status.csv` reports unsettled wakes.

The runner's panel mesh and wake length are starting settings. Check
panel, timestep and wake-length refinement and tighter static/wake tolerances
for quantitative results. Ratios near zero incidence are sensitive to small
angle and load errors, even outside the numerical cutoff. UVLM assumes inviscid
attached flow; this script provides no stall model. Static equilibria do not
establish dynamic stability at the selected speed.

## Implementation checks

The rectangular diagnostic under
`output/static_spanwise_cn_cm/rectangular_diagnostic/` uses the nominal
structure. At 0° and 60 m/s on the archived 5 × 15 lattice,
`eta = 1/15` has a fixed-shape cn slope of 5.23 rad⁻¹ and spar moment
slope of 1.21 rad⁻¹. Both central coefficients are zero. The deflected
checks under `reduced_uvlm_mesh/` use a 3 × 8 lattice. At 3° and 50 m/s
the station `eta = 0.125` gives a secant of 5.77 rad⁻¹ and a fixed-shape
cn slope of 5.43 rad⁻¹. Its cm is 0.000962 about the quarter chord and
0.0694 about the spar. At 1° and 120 m/s the same station gives a secant
of 9.64 rad⁻¹ and a fixed-shape slope of 4.22 rad⁻¹. These runs check
extraction on the stated meshes.

```powershell
julia --project=lib/WingPropellerUVLM/examples lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/test/test_static_spanwise_cn_alpha.jl --plots
```

The small-mesh regressions check independently prescribed local twist angles,
known cn/alpha_local ratios, signed and near-zero incidence, and a single
constant-incidence equilibrium without perturbed solves. They also check
nodal force redistribution, pitching-moment sign and reference, half-width endpoint areas, 21-station output for
20 panels, force conservation on a deformed swept wing, the earlier optional
derivative modes, CSV/plot generation and convergence gates. Upstream-facing
streamwise chords are covered by the existing native geometry tests. The test is included in the Pazy
`test/runtests.jl --full` suite. These checks verify the procedure and do not
establish mesh convergence of the default study.
