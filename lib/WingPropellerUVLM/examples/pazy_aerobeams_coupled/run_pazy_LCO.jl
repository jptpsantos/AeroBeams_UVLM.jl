# Pazy LCO analyses. Direct execution uses independent equilibrium-initialized
# linear upsweep/downsweep runs; the staircase function remains available below.
include(joinpath(@__DIR__, "PazyWingUVLMCoupling.jl"))
include(joinpath(@__DIR__, "PazyWingUVLMVisualization.jl"))
using DelimitedFiles
using Statistics: mean

"""Run one increasing-speed trajectory and plot tip AoA and imposed airspeed.

Defaults: root incidence 3 degrees, 45:1:57 m/s, 20 s per plateau (260 s total).
The final-window half peak-to-peak response is an amplitude estimate; inspect
settling and repeat with longer holds before interpreting it as an LCO branch.
`wingtip_aoa_deg` denotes root incidence + elastic tip twist, without induced
flow or local motion corrections to aerodynamic effective incidence.
"""
function run_pazy_LCO(;
    initial_speed=45.0, final_speed=57.0, speed_step=1.0,
    duration_per_speed=20.0, lco_window_duration=5.0,
    density=1.225, angle_of_attack=deg2rad(3.0), sideslip=0.0,
    chordwise_panels=5, spanwise_panels=16, time_step_chords=0.25,
    wake_length_chords=3.0,
    settling_time=1.0, perturbation_amplitude=0.0, perturbation_duration=0.1,
    coupling_scheme=:strong, coupling_maximum_iterations=20,
    coupling_relaxation=0.5, coupling_geometry_tolerance=1e-4,
    coupling_load_tolerance=1e-4, progress_frequency=100,
    output_directory=joinpath(@__DIR__, "output", "LCO_sweep"))

    @assert all(isfinite, (initial_speed, final_speed, speed_step,
        duration_per_speed, lco_window_duration, time_step_chords, wake_length_chords))
    @assert initial_speed > 0 && final_speed >= initial_speed && speed_step > 0
    @assert 0 < lco_window_duration <= duration_per_speed
    @assert 0 <= settling_time < duration_per_speed && perturbation_duration > 0
    @assert time_step_chords > 0 && wake_length_chords > 0
    speeds = collect(Float64(initial_speed):Float64(speed_step):Float64(final_speed))
    if !isapprox(last(speeds), final_speed; atol=1e-10, rtol=0)
        push!(speeds, final_speed)
    end

    # Use one fixed dt, sized for the highest speed. Shorten it slightly so
    # every hold contains an integer number of steps and no transition drifts.
    _, _, chord, _ = AeroBeams.geometrical_properties_Pazy()
    nominal_dt = time_step_chords * chord / maximum(speeds)
    steps_per_speed = ceil(Int, duration_per_speed / nominal_dt)
    @assert steps_per_speed >= 2 "Each plateau must contain at least two steps"
    dt = duration_per_speed / steps_per_speed
    @assert lco_window_duration >= 2dt "The amplitude window must span at least two timesteps"
    aligned_time_step_chords = dt * maximum(speeds) / chord
    # At a plateau endpoint the OLD speed completes that interval; the new
    # speed is applied on the first integration step of the following plateau.
    plateau_index(t) = clamp(cld(round(Int, t / dt), steps_per_speed), 1, length(speeds))
    airspeed_profile(t) = speeds[plateau_index(t)]
    duration = length(speeds) * duration_per_speed
    # Fixed dt means the lowest speed has the shortest retained wake length.
    maximum_wake_rows = ceil(Int, wake_length_chords * chord / (minimum(speeds) * dt))

    println("\nContinuous Pazy LCO upsweep: $(first(speeds)) -> $(last(speeds)) m/s")
    println("$(length(speeds)) plateaus, $duration_per_speed s each; total $duration s")
    println("One startup only; no state reset or repeated pulse at speed changes.")
    mkpath(output_directory)
    result = run_pazy_wing_uvlm(;
        airspeed=maximum(speeds), airspeed_profile, density, angle_of_attack, sideslip,
        initial_airspeed_fraction=1.0, airspeed_ramp_duration=0.0,
        duration, settling_time, perturbation_amplitude, perturbation_duration,
        chordwise_panels, spanwise_panels, symmetric_wing=iszero(sideslip),
        time_step_chords=aligned_time_step_chords, maximum_wake_rows,
        core_radius=1e-3, wake_shedding_fraction=0.1, aerodynamic_interaction=true,
        save_uvlm_history=false, uvlm_save_frequency=1,
        newton_maximum_iterations=20, newton_absolute_tolerance=1e-5,
        newton_relative_tolerance=1e-5, newton_display_iterations=false,
        newton_always_update_jacobian=true,
        coupling_scheme, coupling_maximum_iterations, coupling_relaxation,
        coupling_geometry_tolerance, coupling_load_tolerance,
        coupling_display_iterations=false, animation_frames=2,
        save_animation_history=false, save_structural_history=false,
        save_aerodynamic_load_history=false, progress_frequency)

    wingtip_aoa_degrees = rad2deg(angle_of_attack) .+ result.tip_twist_degrees
    response = hcat(result.time, result.airspeed_history,
        result.tip_out_of_plane, result.tip_twist_degrees, wingtip_aoa_degrees)
    response_header = ["time_s", "airspeed_mps", "tip_out_of_plane_m",
        "tip_twist_deg", "wingtip_aoa_deg"]
    writedlm(joinpath(output_directory, "pazy_LCO_continuous_history.csv"),
        vcat(permutedims(response_header), response), ',')

    response_amplitudes = Float64[]
    speed_histories = Vector{Float64}[]
    aoa_histories = Vector{Float64}[]
    summary = zeros(length(speeds), 8)
    for (i, speed) in enumerate(speeds)
        first_step = (i - 1) * steps_per_speed + 1
        last_step = i * steps_per_speed + 1
        indices = first_step:last_step
        start_time, end_time = result.time[first_step], result.time[last_step]
        # Exclude the previous plateau's endpoint from amplitude statistics.
        window_indices = [j for j in first_step+1:last_step
            if result.time[j] >= end_time - lco_window_duration]
        length(window_indices) >= 2 || throw(ArgumentError(
            "The final-window amplitude requires at least two samples; increase lco_window_duration"))
        tail = wingtip_aoa_degrees[window_indices]
        low, high = extrema(tail)
        amplitude = (high - low) / 2
        push!(response_amplitudes, amplitude)
        push!(speed_histories, result.time[indices] .- start_time)
        push!(aoa_histories, wingtip_aoa_degrees[indices])
        summary[i, :] .= (speed, start_time, end_time, result.time[first(window_indices)],
            mean(tail), low, high, amplitude)

        # Both global and local time are saved; the first row is the inherited
        # boundary state, including its previously applied speed.
        plateau_response = hcat(result.time[indices] .- start_time, response[indices, :])
        writedlm(joinpath(output_directory, "pazy_LCO_V$(speed).csv"),
            vcat(permutedims(vcat(["plateau_time_s"], response_header)), plateau_response), ',')
        println("U = $speed m/s: final-window wingtip AoA amplitude = $amplitude deg")
    end
    writedlm(joinpath(output_directory, "pazy_LCO_amplitude_summary.csv"),
        vcat(permutedims(["airspeed_mps", "plateau_start_s", "plateau_end_s",
            "window_start_s", "wingtip_aoa_mean_deg", "wingtip_aoa_min_deg",
            "wingtip_aoa_max_deg", "wingtip_aoa_amplitude_deg"]), summary), ',')

    aoa_plot = Plots.plot(result.time, wingtip_aoa_degrees;
        xlabel="Time [s]", ylabel="Wingtip AoA [deg]", color=:dodgerblue,
        linewidth=0.8, legend=false, yguidefontcolor=:dodgerblue,
        ytickfontcolor=:dodgerblue,
        title="Root incidence = $(round(rad2deg(angle_of_attack); digits=3)) deg, increasing velocity",
        titlefontsize=11,
        right_margin=8Plots.mm)
    speed_axis = Plots.twinx(aoa_plot)
    # Draw the command steps at exact plateau boundaries. The sampled speed
    # history uses the integration endpoint convention described above.
    command_time = collect(0:length(speeds)) .* duration_per_speed
    command_speed = vcat(speeds, last(speeds))
    Plots.plot!(speed_axis, command_time, command_speed;
        seriestype=:steppost, ylabel="Airspeed [m/s]", color=:orangered,
        linewidth=1.4, legend=false, grid=false,
        yguidefontcolor=:orangered, ytickfontcolor=:orangered)
    Plots.savefig(aoa_plot, joinpath(output_directory, "pazy_LCO_wingtip_aoa_histories.png"))

    amplitude_plot = Plots.plot(speeds, response_amplitudes;
        marker=:circle, linewidth=2, xlabel="Airspeed [m/s]",
        ylabel="Wingtip AoA amplitude [deg]",
        title="Final-window response (last $(lco_window_duration) s per plateau)",
        label="Half peak-to-peak")
    Plots.savefig(amplitude_plot,
        joinpath(output_directory, "pazy_LCO_wingtip_aoa_amplitude.png"))
    println("\nLCO sweep outputs saved in: $output_directory")

    return (; speeds, response_amplitudes, speed_histories, aoa_histories,
        wingtip_aoa_degrees, result, summary, output_directory)
end

include(joinpath(@__DIR__, "PazyLinearLCOSweep.jl"))

# Direct execution runs independent linear upsweep and downsweep cases at 5 deg,
# each initialized from its own coupled static equilibrium and developed wake.
# Set PAZY_LCO_RUN=false before including to load without running.
# The earlier staircase sweep remains available as run_pazy_LCO().
# Run plot_pazy_LCO.jl to replot strain versus speed from saved branch CSVs.
if get(ENV, "PAZY_LCO_RUN", "true") == "true"
    lco_results = run_pazy_linear_LCO(;
        initial_speed=35.0, peak_speed=52.0,
        upsweep_duration=20.0, downsweep_duration=10.0,
        angle_of_attack=deg2rad(5.0))
end
