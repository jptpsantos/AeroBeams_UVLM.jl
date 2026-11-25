# Pazy wing: AeroBeams + UVLM

Use AeroBeams 1.0 with Julia 1.13 or newer. From the repository root,
instantiate `lib/WingPropellerUVLM/examples/Project.toml` once, then run:

```powershell
julia +1.13 --project=lib/WingPropellerUVLM/examples lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/run_pazy_wing_uvlm.jl
```

The shared environment keeps AeroBeams, UVLM and differentiation dependencies
compatible. Restart existing Julia sessions when switching from the 0.8.1 baseline.

Start with `run_pazy_wing_uvlm.jl`. All inputs are ordinary variables in that
file: airspeed, angle of attack, chordwise/spanwise panel counts, wake length,
duration, tip pulse, animation frame count and the native Pazy structural options.
Run the whole file to simulate,
generate both GIFs, and save the wingtip time-history plot.

For a **static velocity/incidence/sweep study with spanwise cn_alpha**, run
`run_pazy_wing_uvlm_static_sweep.jl`. Edit `velocity_range`,
`angles_of_attack_deg` and `sweep_angles_deg` at its top. By default it holds
the prescribed incidence constant during each static solve and computes
`cn/alpha_local` on the deformed equilibrium wing, without perturbations. It saves
spanwise CSV data, plots and convergence diagnostics in
`output/static_spanwise_cn_alpha/`, with a separate angle folder such as
`aoa3deg/` for each prescribed incidence. See [STATIC_SPANWISE.md](STATIC_SPANWISE.md)
for coefficient definitions, units, controls and validation.

## S0 steady pitch-range comparison

```powershell
julia +1.13 -t 8 --project=lib/WingPropellerUVLM/examples lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/run_pazy_s0_steady_pitch_range.jl
```

The runner solves the upright S0 beam-rod model at 3°, 5°, and 7°, at
0.001 m/s and at integer speeds from 1 to 60 m/s. Both methods use
the nominal structure, a 12-element tip rod, gravity, and the Test 1 balance
mass of 10 g at 40 mm aft of the trailing edge. AeroBeams uses its indicial
strip model with exponential tip loss. The UVLM solution uses a symmetric
4 × 15 panel grid and a 40-row free wake, with incidence imposed through
the freestream. Each airspeed starts with a fresh wake and can continue the
structural equilibrium from the preceding airspeed.

`output/s0_steady_pitch_range/steady_pitch_range.csv` records both solvers'
wingtip response, UVLM convergence residuals, and elapsed time for each case.
The two-panel PNG and PDF use solid lines for native AeroBeams, dashed lines
for AeroBeams + UVLM, and diamonds for the [wind-tunnel reference data][s0-data].
[Avin et al. (2022)][s0-avin] provide displacement data at all three pitches.
[Riso and Cesnik (2023)][s0-riso] estimated twist from the measured leading-
and trailing-edge displacement. The library includes those estimates at 5°
and 7°. Add `--plot-only` to regenerate the figure from the CSV.

[s0-data]: ../../../../test/referenceData/Pazy
[s0-avin]: https://doi.org/10.2514/1.J060621
[s0-riso]: https://doi.org/10.2514/6.2023-0759

The static iteration (`run_pazy_wing_uvlm_static`) accepts options that change only the path to the
equilibrium, not the acceptance tests (wake window, unrelaxed load residual, geometry change):

| Option | Default | Meaning |
| --- | --- | --- |
| `relaxation_mode` | `:fixed` | `:fixed`: loads ← (1-ω) loads + ω G(loads), ω = `static_relaxation`. `:anderson`: Anderson mixing of depth `anderson_depth` (default 4) with mixing parameter `static_relaxation` (use 1.0); needs no tuning and converges in far fewer iterations than any fixed ω. |
| `geometry_cap_chords` | `Inf` | Maximum grid change per load update in chords; the load step is halved (and also when the structural Newton solve fails) until it complies. Recommended 2 with `:anderson`, otherwise the first full step at high load can move the wing so far that its free wake never settles. |
| `wake_reuse` | `false` | Keep the wake and circulation between static iterations (the accepted UVLM state is moved onto the new shape without surface velocity) instead of starting from an empty wake. Each iteration then needs about a window of steps instead of at least two wake ages. |
| `warm_start`, `warm_start_target_loads` | `nothing` | Start from the `warm_start` field of a previous result (same structure and discretization), moving the loads to the target (e.g. previous loads scaled by the dynamic-pressure ratio). The wake is used only with `wake_reuse` and if the entry `aero` is not `nothing` (same time step, e.g. a perturbed incidence). |

`run_pazy_wing_uvlm_static_sweep(; warm_start=true)` starts each equilibrium of a sweep angle from the previous
one and the `:equilibrium`/`:fixed_shape` perturbations from the central one (with `wake_reuse=true` in
`static_options` also from its wake). On the Pazy S0 at U = 10:10:120 m/s and 7° this needed 11 times fewer UVLM steps than
fixed relaxation 0.3, and high-load cases where relaxation 1 fails converge in 10-20 iterations.

The dynamic runner and `run_pazy_wing_uvlm_static.jl` expose the same structural
keywords as the updated `AeroBeams.create_Pazy`:

| Setting | Default | Meaning |
| --- | --- | --- |
| `mesh` | `:nominal` | Native structural mesh; also accepts `:refined`. Independent of the UVLM mesh. |
| `sweepAngle` | `deg2rad(0.0)` | Available sweep angles: 0, 10 and 20 degrees, supplied in radians. |
| `useAePW3Matrices` | `false` | `true` selects the original AePW3 matrices and reference axis; straight wing and nominal mesh only. |
| `withWingtipRod` | `true` | `false` selects the straight wing without its rod, using the no-rod stiffness and inertia data. |
| `rootTorsionCorr` | `true` | Correct root torsional stiffness for the straight wing with rod; ignored for other configurations and AePW3 data. |
| `beamWingtipRod` | `false` | `true` models the flexible rod and rigid link explicitly; available for 0 and 10 degrees with the new matrices. Requires `withWingtipRod=true`. |
| `nElementsWingtipRod` | `12` | Explicit rod mesh; at least 6 elements for the straight wing and 4 for the 10-degree wing. |
| `sweptTipMassConfig` | `:baseline` | Swept-wing mass configuration: `:baseline`, `:LE` or `:TE`. |
| `sweptTipMassInertiaCorrection` | `false` | Native lumped-inertia correction for swept LE/TE cases. Must be false with an explicit rod. |

For example, set `sweepAngle=deg2rad(10.0)` and `sweptTipMassConfig=:TE`
in the runner for the swept TE case. Add `beamWingtipRod=true` to resolve
the rod flexibility. For a straight wing without the rod, keep
`sweepAngle=0.0`, `useAePW3Matrices=false`, `beamWingtipRod=false`, and set
`withWingtipRod=false`.

The shared constructor obtains the matrices, mesh, chord, reference axis and
rod connections from native AeroBeams. Its combination checks and missing-data
errors also apply here. The default now follows the new condensed matrices,
including the native root correction; select `useAePW3Matrices=true` explicitly
to use the original AePW3 data. These choices change the physical model and can
change flutter results.

All aerodynamic loads come from UVLM. Gravity and root pitch remain zero;
`angle_of_attack` sets incidence through the freestream. Swept lattices retain
the physical streamwise chord and use the projected span for reference area
and span. Explicit rod/link beams receive no aerodynamic panels. Tip motion
and the tip pulse refer to the main spar tip. Both solver results record the
structural settings in `pazy_parameters`; a supplied `initial_equilibrium`
must use the same settings as the dynamic call.

The wingtip twist time history is `result.tip_twist_degrees`, sampled at
`result.time`. With `analyze_tip_twist=true`, the runner applies the library's
moving-block damping function to the full signal with the original block defaults,
prints the growth rate, frequency and damping ratio, and displays and
saves `pazy_tip_twist_moving_block_aoa..._V....png`. The supplied
`moving_block_2.jl` is preserved unchanged in `studies/reference`; its
callable implementation is `studies/src/MovingBlockDamping.jl`.

After running the single-case script, the direct call in the Julia REPL is:

```julia
lambda, frequency, omega = moving_block_damping(result.time, result.tip_twist_degrees)
```

The function returns the same three values and displays the same two-panel
diagnostic plot as `moving_block_2.jl`.

The runner also exposes the UVLM core radius, wake-shedding fraction, wing
symmetry, sideslip, interaction flag and aerodynamic-history saving frequency.
`match_wake_panel_chord=true` sets `time_step_chords=1/chordwise_panels`, so the
nominal streamwise wake-row length `U*dt` equals one bound-panel chord. With
`false` (the default), the runner keeps `time_step_chords=0.25`. Wake motion
and induced velocity can change actual panel dimensions. The runner computes
`maximum_wake_rows = ceil(Int, wake_length_chords/time_step_chords)` to keep the
requested nominal wake length in either mode. To use a fixed row count, assign
that integer directly to
`maximum_wake_rows` in the runner. The wake starts empty, grows until this
capacity is reached, and then discards its oldest rows.

The static sweep script exposes the same `match_wake_panel_chord` option and
defaults it to `true`, preserving its existing `time_step_chords=1/chordwise_panels`
setting. Set it to `false` to use `time_step_chords=0.25`; its maximum wake-row
count is then recalculated from the selected step.

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

1. Create the selected native Pazy structure and clamp its root.
2. Set `dt = time_step_chords * chord / airspeed`. The runner uses 0.25 by
   default or `1/chordwise_panels` when wake and bound chord lengths are matched.
3. Create the independent UVLM mesh and its interpolation weights.
4. March both solvers forward in the same time loop.

The default `coupling_scheme=:strong` iterates the aerodynamic loads and beam
deformation within each physical time step until the interface geometry and
load tolerances are satisfied. The `:loose` option uses one aerodynamic trial
on the predicted wing geometry before the structural solve. Both schemes
include wing motion and advance the accepted wake once per physical step.

Each column of `nodal_loads` contains the three forces and three moments at
one beam node, expressed in AeroBeams axes. `apply_pazy_nodal_loads!` updates
the numeric boundary conditions on the active structural model, including
when AeroBeams retries a Newton solve.

The transfer helpers use reusable geometry and load buffers:

- `wing_geometry!`: beam positions and rotations -> UVLM grid.
- `beam_loads!`: UVLM forces -> beam forces and moments, using `arm × force`.
- `PazySpanMap`: cached nonzero interpolation weights for both transfers.
- `save_wing_frame!`: copy the accepted wing and wake for animation, including
  the final state even when the frame stride does not divide the step count.

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

Mesh, timestep and wake-length convergence still need checking for quantitative
flutter results. The nominal wake-row length at the target speed is
`time_step_chords * chord`; it equals a bound panel chord when
`time_step_chords = 1/chordwise_panels`. Wing motion and induced velocity deform
the wake. During the speed ramp, the convected length `U(t)*dt` is smaller.
The final time is rounded to the nearest whole step to retain the timestep.

## Run

Automated checks are grouped in [`test/`](test/README.md). Run the default suite
with `test/runtests.jl`, or append `--full` for static and dynamic integration
tests. Saved-result animation verification is in `validation/`.

Run the focused geometry, finite-force and wake-commit checks with:

```powershell
julia --project=lib/WingPropellerUVLM lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/test/test_coupling.jl
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
