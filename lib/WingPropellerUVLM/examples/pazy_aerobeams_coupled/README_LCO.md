# Independent Pazy LCO sweeps with UVLM

Run `run_pazy_LCO.jl` directly from the IDE or the repository root:

```powershell
julia --project=. lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/run_pazy_LCO.jl
```

The default case consists of two independent simulations at 5 degrees root incidence:

| Branch | Linear speed ramp | Dynamic duration | Initial condition |
| --- | --- | --- | --- |
| Upsweep | 35 to 52 m/s (+0.425 m/s^2) | 40 s | Its own coupled static equilibrium at 35 m/s |
| Downsweep | 52 to 35 m/s (-0.85 m/s^2) | 20 s | Its own coupled static equilibrium at 52 m/s |

Each simulation has a local time axis starting at zero. No structural state,
circulation, or wake passes from the upsweep to the downsweep. There is no
preparation hold. Independent equilibrium initialization is the chosen numerical
initial condition; the quoted paper passage specifies the ramp rates and durations
but does not fully specify its initial states.

Before each ramp, the existing static UVLM/beam solver iterates nonlinear beam
balance and develops a free wake on each fixed trial geometry. Both wake settling
and coupled force/geometry convergence must pass. An unconverged equilibrium
stops the branch before time marching and leaves its convergence history/status
on disk. A fixed waiting time is not substituted for equilibrium.

The dynamic problem copies the converged beam model and full state vector, with
zero structural velocity, and seeds UVLM with that branch's developed wake,
bound circulation, previous geometry and aerodynamic history. Both clocks are
reset to zero. Equilibrium and dynamic timesteps, wake capacities, operating
points, meshes and UVLM settings are checked for compatibility. The static
initialization is preserved separately for inspection.

Edit the parameter block at the bottom of `run_pazy_LCO.jl`, or load without
execution and call the driver:

```julia
ENV["PAZY_LCO_RUN"] = "false"
include("lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/run_pazy_LCO.jl")
delete!(ENV, "PAZY_LCO_RUN")

lco_results = run_pazy_linear_LCO(;
    initial_speed=35.0, peak_speed=52.0,
    upsweep_duration=40.0, downsweep_duration=20.0,
    angle_of_attack=deg2rad(5.0))

up = lco_results.upsweep
down = lco_results.downsweep
# up.result.time and down.result.time both start at zero.
# Each branch has .equilibrium, .result, .root_strain_ue and .wingtip_aoa_degrees.
```

`equilibrium_options` can override static convergence controls, for example
`(; maximum_static_iterations=150, maximum_aerodynamic_time=0.3)`. Defaults match
the standalone static analysis (100 outer iterations, 0.2 s maximum wake time
per fixed-geometry trial, force/circulation tolerances 1e-4, geometry tolerance
1e-5, and two consecutive converged outer iterations).

Both ramps and their static initializations use the same physical timestep,
sized at the maximum speed and slightly reduced to align each branch endpoint.
The default 10 x 20 UVLM mesh and five nominal chord lengths of wake are retained.
Static initialization adds computation, but its wake-development time is not
physical wing-motion time and is not added to either ramp duration.

Outputs in `output/LCO_independent_sweeps`:

- `pazy_linear_LCO_upsweep.csv` and `pazy_linear_LCO_downsweep.csv`: independent
  histories, each starting at time zero, with speed, root strain, raw root axial
  strain/curvatures, tip displacement, twist and incidence-plus-tip-twist AoA.
- `pazy_linear_LCO_upsweep_initialization.csv` and
  `pazy_linear_LCO_downsweep_initialization.csv`: static convergence histories.
- `pazy_linear_LCO_root_strain_vs_airspeed.png`: blue upsweep / red downsweep.
- `pazy_linear_LCO_root_strain_history.png`: separate time plots for each branch.
- `pazy_linear_LCO_settings.txt`: controls and initialization convergence statuses.

There is no concatenated physical time history. The new default output directory
keeps these independent runs separate from the previous connected-ramp outputs.

The coupling retains near-root section values even with detailed histories
disabled. Surface strain uses local beam coordinates (+y forward, +z upper):

    strain_ue = 1e6 * (gamma_1 + gauge_z*kappa_2 - gauge_y*kappa_3)

Defaults `gauge_y=0`, `gauge_z=2.25e-3/2` place the gauge on the upper spar surface
at its chordwise center. The axial component of `gamma + kappa x offset` fixes
the sign: upper-surface compression under upward bending is negative. The first
element midpoint approximates the root station. Confirm experimental gauge
location and section conventions before quantitative comparison.

This changes the sweep and initialization procedure. The current beam coupling
still uses zero gravity and UVLM aerodynamics, while the paper includes gravity
and compares linear and dynamic-stall aerodynamic models. The dashed steady
reference curve requires separate converged static solutions over the speed
range, using the same strain reconstruction. A dynamic mean is not that curve.

The earlier staircase analysis remains available as `run_pazy_LCO()` after
loading with execution disabled, and writes to `output/LCO_sweep`.

Regression checks exercise two independent static initializations and ramps,
verify copied beam/wake states and zero starting velocity, require independent
clocks, reject unconverged initialization, and check exports and plots:

```powershell
julia --startup-file=no --compiled-modules=existing --project=. lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/test_pazy_LCO.jl
```