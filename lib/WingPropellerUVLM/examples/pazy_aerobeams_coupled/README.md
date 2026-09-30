# Pazy wing: AeroBeams + UVLM

Start with `run_pazy_wing_uvlm.jl`. All inputs are ordinary variables in that
file: airspeed, angle of attack, chordwise/spanwise panel counts, wake length,
duration, tip pulse and animation frame count. Run the whole file to simulate,
generate both GIFs, and save the wingtip time-history plot.

The runner also exposes the UVLM core radius, wake-shedding fraction, wing
symmetry, sideslip, interaction flag and aerodynamic-history saving frequency.
`wake_length_chords = 20.0` retains approximately 20 wing chords of wake:
the runner calculates `maximum_wake_rows = ceil(Int, wake_length_chords *
chordwise_panels)`. With four chordwise panels this is 80 rows, or 1.978 m
nominal length. To use a fixed row count, assign that integer directly to
`maximum_wake_rows` in the runner. The wake starts empty, grows until this
capacity is reached, and then discards its oldest rows.

`core_radius` is in metres. `wake_shedding_fraction` controls the separation
between the trailing edge and the wake-shedding location, not the retained wake
length. Wake convection includes induced velocity. `aerodynamic_interaction`
controls influence between surface groups and has no effect for this single
wing; it does not disable wake roll-up. `symmetric_wing=true` includes an image
wing and wake across UVLM y=0. Set it to `false` for an isolated cantilever or
nonzero sideslip.

With `save_uvlm_history=true`, `result.aerodynamic_problem.results` contains
`savedTimeVector`, `circulationOverTime`, `forceOverTime`, and `momentOverTime`.
These are sampled at `uvlm_save_frequency` accepted steps, plus the final step;
the initial unloaded state is not included. Forces and moments use UVLM axes
and refer to the explicitly modelled wing (the symmetry image is not added to
the totals). Set `generate_animations=false` to skip GIF generation, or change
`animation_fps` and `output_directory` in the runner. The accepted UVLM loads
transferred to the beam are available in `result.aerodynamic_nodal_load_history`.

`progress_frequency = 1` prints every accepted coupled timestep. Increase it
to, for example, `100` to print only every 100 steps; the final step is always
printed.

The startup uses a cubic smooth airspeed ramp when
`initial_airspeed_fraction < 1` and `airspeed_ramp_duration > 0`, and then
holds full speed until the tip pulse starts at `settling_time`.
`result.airspeed_history` stores the applied sweep. UVLM requires a positive
speed, so `initial_airspeed_fraction` must be greater than zero.

The ramp and hold reduce the startup transient, but they do not guarantee that
coupled equilibrium has been reached. Verify that the tip motion and
aerodynamic loads are nearly constant before the pulse, and increase
`settling_time` when necessary.

With `newton_display_iterations = true`, AeroBeams prints each Newton iteration:
`i` is the iteration, `σ` is the load factor, `ϵ_abs` is the residual norm, and
`ϵ_rel` is the relative solution correction. The timestep is accepted only
when `convergedFinalSolution` is true; otherwise the example stops. The Newton
iteration limit and tolerances are also defined in the runner.

`PazyWingUVLMCoupling.jl` contains one main function, read from top to bottom:

1. Create the original 15-element nonlinear Pazy beam and clamp its root.
2. Set `dt = time_step_chords * chord / airspeed`, independently of the UVLM
   mesh. The main runner uses Ribeiro's baseline `time_step_chords = 0.25`.
3. Create the independent UVLM mesh and its interpolation weights.
4. March both solvers forward in the same time loop.

Inside that loop, UVLM predicts forces, AeroBeams advances the beam, and UVLM
recomputes the circulation on the deformed wing and advances the wake once.
The force predictor extrapolates the last two wing geometries to the new time,
so UVLM includes wing motion instead of incorrectly treating the wing as stationary.
The array `nodal_loads` connects them: each column contains the three forces
and three moments applied to one beam node. AeroBeams reads that array through
its boundary-condition functions.

Only three small helpers follow the main function:

- `wing_geometry`: beam positions and rotations -> UVLM grid.
- `beam_loads`: UVLM forces -> beam forces and moments, using `arm × force`.
- `save_wing_frame!`: copy the wing and wake for animation.

The aerodynamic root uses the exact prescribed clamp position and rotation.
It does not use the small residual displacements in the beam's recovered root
outputs: even a tiny offset from the symmetry plane can make the root vortex
nearly coincide with its image and produce nonfinite aerodynamic forces.

The simulation lasts approximately `duration`. The tip pulse starts at
`settling_time`;
the wing and wake continue from their existing state. This is a prescribed
settling interval, not an automatic equilibrium check. Check the tip history
to see whether the startup transient has decayed before interpreting the pulse
response. Neither the structure nor UVLM is linearized.

The coupling is explicit: the structural step uses predicted aerodynamic loads.
Mesh, timestep and wake-length convergence still need checking for quantitative
flutter results. The timestep sets the *nominal* wake-row length to one bound
panel chord; wing motion and induced velocity deform the wake. The final time
is rounded to the nearest whole step to retain this timestep.
During the speed ramp, `U(t)*dt` is smaller than a bound panel chord; the nominal
length matches only after the target speed is reached.

## Run

Run the focused geometry, finite-force and wake-commit checks with:

```powershell
julia --project=lib/WingPropellerUVLM lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/test_coupling.jl
```

Append `--dynamic` to also run a short coupled startup-and-pulse test. This
checks the implementation; it does not replace the 3-second run or a flutter
convergence study.

From the repository root, prepare the environments once:

```powershell
julia --project=. -e "using Pkg; Pkg.instantiate()"
julia --project=lib/WingPropellerUVLM -e "using Pkg; Pkg.instantiate()"
```

Then run:

```powershell
julia --project=lib/WingPropellerUVLM lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/run_pazy_wing_uvlm.jl
```

You can also run `run_pazy_wing_uvlm.jl` directly from the IDE while the
repository-root project is active. The example adds the local AeroBeams and
WingPropellerUVLM package directories to Julia's load path itself.

The runner adds the flight-condition suffix to all outputs, for example
`pazy_response_structure_aoa3_V50.gif`, `pazy_response_uvlm_wake_aoa3_V50.gif`,
and `pazy_tip_time_histories_aoa3_V50.png`, inside `output_directory`.

The structural GIF uses AeroBeams' own deformation renderer, physical scale
1, transparent blue NACA 0018 surfaces, an emphasized beam elastic axis, and
the native clamped-root symbols. The run file exposes:

```julia
structure_backend = :pyplot       # supports translucent 3D surfaces
structure_camera = (45, 45)       # azimuth, elevation [degrees]
structure_surface_alpha = 0.5     # opacity: smaller means more transparent
structure_axis_linewidth = 2.5
show_clamped_root = true
```

PyPlot is already an AeroBeams dependency. The wake GIF retains its own GR
renderer and `wake_camera`. The structural camera and semispan-based axis
limits are fixed across frames. Airfoil surfaces are attached only while
plotting; this does not enable AeroBeams aerodynamic loads or alter the
computed response.

On Windows, the structural renderer temporarily selects sequential NumPy MKL
to avoid the threaded-library crash observed in the 3D preview test. The
environment is restored after rendering; Julia's solver BLAS settings are
not changed. Use a fresh Julia session if PyPlot/NumPy was already initialized
with a different MKL threading layer.

The clamp display is independent of `show_force_vectors` and
`show_moment_vectors`, which remain off by default. If enabled, the structural
GIF shows lift and torsional-moment arrows; the wake GIF shows the full force
vectors. Arrow normalization is fixed over the simulation.

To regenerate just the structural GIF with an existing in-memory `result`,
reload `PazyWingUVLMVisualization.jl` and call:

```julia
save_structural_animation(result;
    output_path=joinpath(output_directory, "pazy_response_structure_$(case_suffix).gif"),
    camera=(45,45), surface_alpha=0.5, show_clamped_root=true)
```

This replaces only that GIF; it does not rerun the simulation. The tip plot
shows out-of-plane bending at the structural reference line and twist using
the rotated-local-chord definition of the AeroBeams Pazy examples.
