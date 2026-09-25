# This file can be run directly from the IDE or from either local Julia project.
include(joinpath(@__DIR__, "PazyWingUVLMCoupling.jl"))
include(joinpath(@__DIR__, "PazyWingUVLMVisualization.jl"))

# Flight conditions
airspeed = 50.0                       # m/s
density = 1.225                       # kg/m^3
angle_of_attack = deg2rad(3.0)
sideslip = deg2rad(0.0)               # use symmetric_wing=false for nonzero sideslip

# Smooth startup: ramp from 5% to 100% airspeed, then hold before the tip pulse.
# A positive initial fraction is required because UVLM does not accept U=0.
initial_airspeed_fraction = 1.0
airspeed_ramp_duration = 0.0           # seconds; must not exceed settling_time

# UVLM mesh (independent of the 15-element AeroBeams mesh)
chordwise_panels = 5
spanwise_panels = 15
symmetric_wing = true                 # mirror wing and wake across UVLM y=0
                                     # false: isolated cantilever, without its image

# Retained wake length. One row is shed per timestep; old rows are discarded.
# Since U*dt = chord/Nc, retaining N rows gives nominal length N*chord/Nc.
wake_length_chords = 10.0             # nominal retained length / wing chord
@assert isfinite(wake_length_chords) && wake_length_chords > 0
maximum_wake_rows = ceil(Int, wake_length_chords * chordwise_panels)
# To specify the row count directly, replace the expression above with e.g. 80.
# The wake starts empty and grows to this limit. Induced velocities deform it.

# UVLM solver controls
core_radius = 1e-3                    # vortex core radius [m], strictly positive
wake_shedding_fraction = 0.1          # trailing-edge offset = fraction * V_relative * dt
                                     # allowed range: 0 to 1; this does not set wake length
aerodynamic_interaction = true       # influence between different surface groups
                                     # no effect with this single wing; NOT a roll-up switch

# Store UVLM circulation and total force/moment histories in result.
save_uvlm_history = true
uvlm_save_frequency = 1               # save every N accepted steps (final step always saved)

# AeroBeams Newton-Raphson controls
newton_maximum_iterations = 20
newton_absolute_tolerance = 1e-6
newton_relative_tolerance = 1e-6
newton_display_iterations = true       # print i, load factor and convergence errors
newton_always_update_jacobian = true

# Strong UVLM-AeroBeams coupling at every physical timestep. The geometry and
# aerodynamic loads must both converge before the wake is advanced.
coupling_maximum_iterations = 20
coupling_relaxation = 0.5              # smaller is more robust but slower
coupling_geometry_tolerance = 1e-5     # max interface-coordinate change / chord
coupling_load_tolerance = 1e-5         # relative nodal force/moment change
coupling_display_iterations = false    # true: print every inner FSI iteration

# Time: dt is calculated as chord / (chordwise_panels * airspeed).
duration = 5.0                        # total simulated time, approximately [s]
settling_time = 1.0                   # time at which the tip pulse starts [s]
perturbation_amplitude = 0.0         # tip-force multiplier [N]
perturbation_duration = 0.1          # pulse duration [s]
progress_frequency = 1                # 1: print every step; 100: every 100 steps

# Animation output (independent of UVLM force-history storage)
generate_animations = false
generate_time_history_plot = true
animation_time_step = 0.01              # simulated seconds between frames, as in Pazy gust
animation_frames = 150                  # fallback maximum if animation_time_step=nothing
animation_fps = 30                      # same requested FPS as the AeroBeams example
wake_camera = (45, 30)                  # high oblique view: (azimuth, elevation) [deg]
show_force_vectors = true               # green UVLM force arrows in both GIFs
force_vector_scale = 1.0                # 1.0: largest arrow is about span/10
output_directory = joinpath(@__DIR__, "output")

# Solve, then generate the AeroBeams deformation and UVLM wake GIFs.
result = run_pazy_wing_uvlm(;
    airspeed, density, angle_of_attack, sideslip,
    initial_airspeed_fraction, airspeed_ramp_duration,
    duration, settling_time,
    chordwise_panels, spanwise_panels, symmetric_wing, maximum_wake_rows,
    core_radius, wake_shedding_fraction, aerodynamic_interaction,
    save_uvlm_history, uvlm_save_frequency,
    newton_maximum_iterations, newton_absolute_tolerance,
    newton_relative_tolerance, newton_display_iterations,
    newton_always_update_jacobian,
    perturbation_amplitude, perturbation_duration, animation_frames,
    progress_frequency, animation_time_step,
    coupling_maximum_iterations, coupling_relaxation,
    coupling_geometry_tolerance, coupling_load_tolerance,
    coupling_display_iterations)
if generate_animations
    animations = save_animations(
        result;
        output_directory,
        fps=animation_fps,
        wake_camera,
        show_force_vectors,
        force_vector_scale,
    )
end
if generate_time_history_plot
    tip_time_history_plot = save_tip_time_histories(
        result;
        output_path=joinpath(output_directory, "pazy_tip_time_histories.png"),
    )
end
println("Final tip displacement [m]: ", last(result.tip_out_of_plane))
println("Final tip twist [deg]: ", last(result.tip_twist_degrees))
println("Maximum FSI iterations in one step: ", maximum(result.coupling_iterations))
# Applied startup sweep: result.airspeed_history
# Aerodynamic histories: result.aerodynamic_problem.results
# Fields: savedTimeVector, circulationOverTime, forceOverTime, momentOverTime.
# Forces/moments use UVLM axes and include only the explicitly modelled wing.
