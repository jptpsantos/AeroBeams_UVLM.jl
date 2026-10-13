# Run Pazy wing cases, estimate damping, and plot the stability boundary.
let root = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))
    root in LOAD_PATH || pushfirst!(LOAD_PATH, root)
end
using Plots, DelimitedFiles, ColorSchemes
include(joinpath(@__DIR__, "..", "..", "studies", "src", "MovingBlockDamping.jl"))
using .MovingBlockDamping: moving_block_metrics
include(joinpath(@__DIR__, "PazyStabilitySweepTools.jl"))

# Run every combination of these pitch angles and airspeeds.
pitch_angles_deg = [0.0, 1.0, 3.0, 5.0, 7.0]
airspeeds_mps = collect(30.0:10.0:90.0)
run_simulations = true                 # false: read the saved response.csv files
output_directory = joinpath(@__DIR__, "output", "stability_sweep")
if run_simulations
    include(joinpath(@__DIR__, "PazyWingUVLMCoupling.jl"))
end

# Flight and simulation settings, using the same names as run_pazy_wing_uvlm.jl.
density = 1.225
sideslip = 0.0
initial_airspeed_fraction = 1.0       # UVLM requires a positive starting speed
airspeed_ramp_duration = 0.0
duration = 2.0
settling_time = 0.2
perturbation_amplitude = 10.0
perturbation_duration = 0.1

chordwise_panels = 10
spanwise_panels = 20
symmetric_wing = true
match_wake_panel_chord = false       # true: nominal U*dt = chord/chordwise_panels
time_step_chords = 0.25
if match_wake_panel_chord
    time_step_chords = 1 / chordwise_panels
end
wake_length_chords = 5.0
maximum_wake_rows = ceil(Int, wake_length_chords / time_step_chords)
core_radius = 1e-3
wake_shedding_fraction = 0.1
aerodynamic_interaction = true

save_uvlm_history = false
uvlm_save_frequency = 1
newton_maximum_iterations = 30
newton_absolute_tolerance = 1e-5
newton_relative_tolerance = 1e-5
newton_display_iterations = false
newton_always_update_jacobian = true
coupling_scheme = :strong
coupling_maximum_iterations = 30
coupling_relaxation = 0.5
coupling_geometry_tolerance = 1e-4
coupling_load_tolerance = 1e-4
coupling_display_iterations = false
animation_frames = 2
animation_time_step = 3.0
save_animation_history = false       # Skip frame copies during the stability sweep
progress_frequency = 100

# Moving-block damping uses the full wingtip history and the original defaults.
response_to_fit = :twist             # :twist or :bending
boundary_minimum_fit_r_squared = 0.90 # only for accepting a boundary estimate
growth_deadband = 0.05                # Treat near-zero growth as unresolved
midpoint_bracket_angles_deg = [4.0, 6.0] # Saturating responses: use speed-bracket midpoints

@assert response_to_fit in (:twist, :bending)
@assert all(diff(airspeeds_mps) .> 0) "List speeds in increasing order."
mkpath(output_directory)
println("Running $(length(pitch_angles_deg)*length(airspeeds_mps)) angle/speed cases.")

summary_rows = Any[]
boundary_rows = Any[]
bracket_rows = Any[]
growth = fill(NaN, length(pitch_angles_deg), length(airspeeds_mps))
raw_growth = fill(NaN, size(growth))
case_status = fill("failed", size(growth))

# Run each case and estimate its damping.
for (i, angle) in enumerate(pitch_angles_deg), (j, speed) in enumerate(airspeeds_mps)
    label = "aoa$(sweep_number(angle))_V$(sweep_number(speed))"
    directory = joinpath(output_directory, label)
    mkpath(directory)
    println("\nRunning $label")
    summary_row = [angle, speed, NaN, NaN, NaN, NaN, "failed"]

    try
        if run_simulations
            result = run_pazy_wing_uvlm(;
                airspeed=speed, density, angle_of_attack=deg2rad(angle), sideslip,
                initial_airspeed_fraction, airspeed_ramp_duration,
                duration, settling_time, perturbation_amplitude, perturbation_duration,
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
            response = hcat(result.time, result.airspeed_history,
                result.tip_bending_displacement, result.tip_twist_degrees)
            write_sweep_csv(joinpath(directory, "response.csv"),
                ["time_s", "airspeed_mps", "tip_bending_m", "tip_twist_deg"], response)
        else
            response = readdlm(joinpath(directory, "response.csv"), ',', Float64; skipstart=1)
        end
        column = response_to_fit == :twist ? 4 : 3
        all(u -> isapprox(u, speed; rtol=1e-8), response[:, 2]) ||
            error("Airspeed is not constant in the response")
        fit = moving_block_metrics(response[:, 1], response[:, column])
        s = fit.summary
        raw_growth[i, j] = s.growth_rate_per_s

        status = "unresolved"
        if s.available && s.fit_r_squared >= boundary_minimum_fit_r_squared &&
           abs(s.growth_rate_per_s) > growth_deadband
            growth[i, j] = s.growth_rate_per_s
            status = growth[i, j] < 0 ? "decaying" : "growing"
        end
        summary_row = [angle, speed, s.growth_rate_per_s, s.damping_ratio_percent,
            s.frequency_hz, s.fit_r_squared, status]
        println("$label: $status, damping = $(s.damping_ratio_percent) %")

        response_plot = plot(response[:, 1], response[:, column];
            label=string(response_to_fit), xlabel="Time [s]",
            ylabel=response_to_fit == :twist ? "Twist [deg]" : "Bending [m]", title=label)
        if s.available
            vspan!(response_plot, [s.analyzed_start_s, s.analyzed_end_s];
                alpha=0.15, label="Analyzed segment")
        end
        diagnostics = fit.diagnostics
        fit_plot = plot(diagnostics.block_start_s, diagnostics.log_amplitude;
            label="Moving blocks", xlabel="Block start [s]", ylabel="Log amplitude")
        plot!(fit_plot, diagnostics.block_start_s, diagnostics.fitted_log_amplitude;
            label="Linear fit")
        savefig(plot(response_plot, fit_plot; layout=(2, 1), size=(800, 650)),
            joinpath(directory, "moving_block.png"))
    catch err
        err isa InterruptException && rethrow()
        # A solver failure is not evidence of physical instability.
        growth[i, j] = NaN
        raw_growth[i, j] = NaN
        summary_row = [angle, speed, NaN, NaN, NaN, NaN, "failed"]
        write(joinpath(directory, "failure.txt"), sprint(showerror, err))
        @warn "Case $label failed; continuing" exception=err
    end

    case_status[i, j] = summary_row[7]
    push!(summary_rows, summary_row)
    write_sweep_csv(joinpath(output_directory, "damping_summary.csv"),
        ["angle_deg", "speed_mps", "growth_per_s", "damping_percent", "frequency_hz", "fit_r_squared", "status"],
        permutedims(hcat(summary_rows...)))
end

# Provisional hump boundaries use a positive-growth unresolved fit and the
# neighboring confirmed decaying case. Keep the strict boundaries separately.
estimated_boundary_rows = Any[]
for (i, angle) in enumerate(pitch_angles_deg)
    indices = eachindex(airspeeds_mps)
    first_positive = findfirst(j -> isfinite(raw_growth[i, j]) &&
        raw_growth[i, j] > growth_deadband, indices)
    isnothing(first_positive) && continue
    onset_left = findlast(j -> j < first_positive && case_status[i, j] == "decaying", indices)
    offset_right = findfirst(j -> j > first_positive && case_status[i, j] == "decaying", indices)
    offset_left = isnothing(offset_right) ? nothing :
        findlast(j -> j < offset_right && raw_growth[i, j] > growth_deadband, indices)
    for (kind, left, right) in (("onset", onset_left, first_positive),
                               ("offset", offset_left, offset_right))
        (isnothing(left) || isnothing(right)) && continue
        all(s -> s in ("decaying", "growing", "unresolved"), case_status[i, left:right]) || continue
        u1, u2 = airspeeds_mps[left], airspeeds_mps[right]
        g1, g2 = raw_growth[i, left], raw_growth[i, right]
        speed = u1 - g1*(u2-u1)/(g2-g1)
        quality = any(==("unresolved"), case_status[i, left:right]) ?
            "provisional_unresolved" : "resolved"
        if angle in midpoint_bracket_angles_deg && quality == "provisional_unresolved"
            # The fitted slope mixes initial growth and amplitude saturation.
            # A midpoint expresses the speed bracket without using that slope.
            speed = (u1 + u2)/2
            quality = "bracket_midpoint_unresolved"
        end
        push!(estimated_boundary_rows, [angle, speed, kind, quality, u1, u2, g1, g2])
    end
end
estimated_data = isempty(estimated_boundary_rows) ? Matrix{Any}(undef, 0, 8) :
    permutedims(hcat(estimated_boundary_rows...))
write_sweep_csv(joinpath(output_directory, "uvlm_boundary_updated.csv"),
    ["angle_deg", "speed_mps", "boundary", "quality", "lower_speed_mps", "upper_speed_mps",
     "lower_growth_per_s", "upper_growth_per_s"], estimated_data)

# Interpolate sign changes to estimate onset and restabilization speeds.
for (i, angle) in enumerate(pitch_angles_deg)
    for crossing in sweep_crossings(airspeeds_mps, growth[i, :])
        push!(boundary_rows, [angle, crossing.speed, crossing.kind])
    end
    # If an unresolved case lies between opposite signs, report the speed
    # interval without interpolating a precise crossing through missing data.
    valid = findall(isfinite, growth[i, :])
    for (left, right) in zip(valid[1:end-1], valid[2:end])
        right > left + 1 || continue
        all(==("unresolved"), case_status[i, left+1:right-1]) || continue
        growth[i, left] * growth[i, right] < 0 || continue
        kind = growth[i, left] < 0 ? "onset" : "offset"
        push!(bracket_rows, [angle, airspeeds_mps[left], airspeeds_mps[right], kind])
    end
end
boundary_data = isempty(boundary_rows) ? Matrix{Any}(undef, 0, 3) : permutedims(hcat(boundary_rows...))
write_sweep_csv(joinpath(output_directory, "uvlm_boundary.csv"),
    ["angle_deg", "speed_mps", "boundary"], boundary_data)
bracket_data = isempty(bracket_rows) ? Matrix{Any}(undef, 0, 4) : permutedims(hcat(bracket_rows...))
write_sweep_csv(joinpath(output_directory, "uvlm_boundary_brackets.csv"),
    ["angle_deg", "lower_speed_mps", "upper_speed_mps", "boundary"], bracket_data)

# Edit and rerun only this block to change the plot after the sweep.
# Use the AeroBeams reference colors for tests, UM/NAST, and SHARPy.
aerobeams_color = get(ColorSchemes.darkrainbow, 0.05)
aerobeams_styles = (:solid, :solid)
experimental_colors = (:red, :green)
literature_colors = (:gold, :brown, :magenta)
uvlm_color = :red
uvlm_fill_alpha = 0.25
comparison_labels = (
    uvlm_region="AeroBeams/UVLM time domain (estimated)",
    aerobeams_onset="AeroBeams eigenvalue onset", aerobeams_offset="",
    test_onset_up="Test onset up", test_offset_up="Test offset up",
    test_onset_down="Test onset down", test_offset_down="Test offset down",
    umnast_loss="UM/NAST (Exp. loss)", umnast_panel="UM/NAST (Panel coeffs.)",
    sharpy="Sharpy",
)
ts = 10
fs = 16
lfs = 9
lw = 2
ms = 10
ms2 = 4
msw = 0
comparison_plot = plot_sweep_comparison(estimated_boundary_rows, pitch_angles_deg;
    bracket_rows, aerobeams_color, aerobeams_styles, literature_colors,
    experimental_colors, uvlm_color, uvlm_fill_alpha, labels=comparison_labels,
    line_width=lw, marker_size=ms, literature_marker_size=ms2,
    marker_stroke_width=msw)
plot!(comparison_plot;
    xlabel="Airspeed [m/s]", ylabel="Root pitch angle [deg]",
    xlims=[30, 90], ylims=[0, 7.25],
    xticks=collect(30:10:90), yticks=collect(0:1:7),
    legend=:topright, tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
savefig(comparison_plot, joinpath(output_directory, "pazy_stability_boundary.png"))
savefig(comparison_plot, joinpath(output_directory, "pazy_stability_boundary.svg"))
display(comparison_plot)
println("Results saved to $output_directory")
