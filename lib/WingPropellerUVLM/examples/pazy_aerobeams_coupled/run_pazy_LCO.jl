# Speed sweep for the Pazy wing LCO response. Run directly from the IDE or
# either local Julia project, like run_pazy_wing_uvlm.jl.
include(joinpath(@__DIR__, "PazyWingUVLMCoupling.jl"))
include(joinpath(@__DIR__, "PazyWingUVLMVisualization.jl"))
using DelimitedFiles

# Flight condition and speed sweep [m/s].
initial_speed = 45.0
final_speed = 60.0
speed_step = 5.0
density = 1.225
angle_of_attack = deg2rad(3.0)
sideslip = deg2rad(0.0)

@assert initial_speed > 0 && final_speed >= initial_speed
@assert speed_step > 0
speeds = collect(initial_speed:speed_step:final_speed)
if !isapprox(last(speeds), final_speed; atol=1e-10, rtol=0)
    push!(speeds, final_speed)
end

# Simulate 10 seconds at each speed. The short tip-force pulse excites the
# wing after the settling interval; it is then free to reach its LCO response.
duration_per_speed = 10.0
settling_time = 1.0
perturbation_amplitude = 0.0
perturbation_duration = 0.1
lco_window_duration = 2.0      # final interval used for the amplitude estimate

# UVLM mesh and wake settings from run_pazy_wing_uvlm.jl.
chordwise_panels = 10
spanwise_panels = 20
symmetric_wing = true
initial_airspeed_fraction = 1.0
airspeed_ramp_duration = 0.0
match_wake_panel_chord = false
time_step_chords = match_wake_panel_chord ? 1 / chordwise_panels : 0.25
wake_length_chords = 5.0
maximum_wake_rows = ceil(Int, wake_length_chords / time_step_chords)
core_radius = 1e-3
wake_shedding_fraction = 0.1
aerodynamic_interaction = true

# Solver controls from the base Pazy run.
save_uvlm_history = false
uvlm_save_frequency = 1
newton_maximum_iterations = 20
newton_absolute_tolerance = 1e-5
newton_relative_tolerance = 1e-5
newton_display_iterations = false
newton_always_update_jacobian = true
coupling_scheme = :strong
coupling_maximum_iterations = 20
coupling_relaxation = 0.5
coupling_geometry_tolerance = 1e-4
coupling_load_tolerance = 1e-4
coupling_display_iterations = false
animation_frames = 2
animation_time_step = 10.0
save_animation_history = false
progress_frequency = 100

@assert 0 < lco_window_duration <= duration_per_speed
output_directory = joinpath(@__DIR__, "output", "LCO_sweep")
mkpath(output_directory)

speed_histories = Vector{Vector{Float64}}()
aoa_histories = Vector{Vector{Float64}}()
response_amplitudes = Float64[]

for speed in speeds
    println("\nRunning Pazy LCO case at U = $speed m/s for $duration_per_speed s")
    result = run_pazy_wing_uvlm(;
        airspeed=speed, density, angle_of_attack, sideslip,
        initial_airspeed_fraction, airspeed_ramp_duration,
        duration=duration_per_speed, settling_time,
        perturbation_amplitude, perturbation_duration,
        chordwise_panels, spanwise_panels, symmetric_wing,
        time_step_chords, maximum_wake_rows,
        core_radius, wake_shedding_fraction, aerodynamic_interaction,
        save_uvlm_history, uvlm_save_frequency,
        newton_maximum_iterations, newton_absolute_tolerance,
        newton_relative_tolerance, newton_display_iterations,
        newton_always_update_jacobian, animation_frames, progress_frequency,
        animation_time_step, save_animation_history, coupling_scheme,
        coupling_maximum_iterations, coupling_relaxation,
        coupling_geometry_tolerance, coupling_load_tolerance,
        coupling_display_iterations)

    # The coupling helper reports elastic twist relative to the undeformed
    # wing. Add it to the fixed freestream AoA for the wingtip effective AoA.
    wingtip_aoa_degrees = rad2deg(angle_of_attack) .+ result.tip_twist_degrees
    push!(speed_histories, result.time)
    push!(aoa_histories, wingtip_aoa_degrees)

    window_start = duration_per_speed - lco_window_duration
    window_indices = findall(t -> t >= window_start, result.time)
    amplitude = (maximum(wingtip_aoa_degrees[window_indices]) -
        minimum(wingtip_aoa_degrees[window_indices])) / 2
    push!(response_amplitudes, amplitude)

    response = hcat(result.time, fill(speed, length(result.time)),
        result.tip_twist_degrees, wingtip_aoa_degrees)
    writedlm(joinpath(output_directory, "pazy_LCO_V$(speed).csv"),
        vcat(permutedims(["time_s", "airspeed_mps", "tip_twist_deg", "wingtip_aoa_deg"]), response), ',')
    println("Wingtip AoA final-window amplitude: $amplitude deg")
end

# Plot wingtip effective AoA histories for all analyzed speeds.
aoa_plot = Plots.plot(; xlabel="Time [s]", ylabel="Wingtip AoA [deg]",
    title="Pazy wingtip AoA response by airspeed", legend=:outerright)
for (i, speed) in enumerate(speeds)
    Plots.plot!(aoa_plot, speed_histories[i], aoa_histories[i];
        label="$(speed) m/s", linewidth=1.2)
end
Plots.savefig(aoa_plot, joinpath(output_directory, "pazy_LCO_wingtip_aoa_histories.png"))

# Plot half peak-to-peak AoA over the final window as a compact LCO amplitude
# curve. This includes any steady AoA offset only through the oscillation range.
amplitude_plot = Plots.plot(speeds, response_amplitudes;
    marker=:circle, linewidth=2, xlabel="Airspeed [m/s]",
    ylabel="Wingtip AoA amplitude [deg]",
    title="Pazy wingtip LCO amplitude (final $(lco_window_duration) s)",
    label="Half peak-to-peak")
Plots.savefig(amplitude_plot,
    joinpath(output_directory, "pazy_LCO_wingtip_aoa_amplitude.png"))

amplitude_table = hcat(speeds, response_amplitudes)
writedlm(joinpath(output_directory, "pazy_LCO_amplitude_summary.csv"),
    vcat(permutedims(["airspeed_mps", "wingtip_aoa_amplitude_deg"]), amplitude_table), ',')
println("\nLCO sweep outputs saved in: $output_directory")
