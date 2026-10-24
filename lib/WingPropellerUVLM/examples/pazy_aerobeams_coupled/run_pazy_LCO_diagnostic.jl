# Constant-speed pulse and numerical-sensitivity checks, separate from the ramps.
if !isdefined(@__MODULE__, :run_pazy_wing_uvlm_static)
    include(joinpath(@__DIR__, "PazyWingUVLMStatic.jl"))
end
using DelimitedFiles, LinearAlgebra, Statistics, Printf
import Plots

"""Post-pulse windows; detrending removes slow drift, not the oscillatory signal."""
function pazy_diagnostic_windows(time, strain; pulse_end, window_duration)
    rows = Vector{Float64}[]
    start = pulse_end
    while start + window_duration <= last(time) + 10eps(last(time))
        stop = start + window_duration
        indices = findall(t -> start <= t < stop, time)
        length(indices) >= 4 || throw(ArgumentError("Response windows need at least four samples"))
        t = time[indices] .- mean(time[indices])
        values = strain[indices]
        residual = values .- mean(values)
        residual .-= t .* (dot(t, residual) / dot(t, t))
        low, high = extrema(values)
        amplitude = (maximum(residual) - minimum(residual)) / 2
        rms = sqrt(mean(abs2, residual))
        push!(rows, [start, stop, mean(values), low, high, (high-low)/2, amplitude, rms])
        start = stop
    end
    isempty(rows) && throw(ArgumentError("Duration must contain a complete post-pulse response window"))
    return permutedims(hcat(rows...))
end

"""Compare fixed-speed forced/unforced responses, half dt, and double wake.

Each distinct speed/dt/wake combination has its own converged static equilibrium
and developed wake. Baseline and unforced runs COPY the same equilibrium. Mesh,
coupling tolerances, and physical tip pulse are identical across refinements.
The timestep and wake row rules reproduce run_pazy_linear_LCO's default sweeps
(52 m/s reference speed, 35 m/s minimum speed, 20 s time alignment).

Window amplitudes, RMS, drift, and tail growth are diagnostics, not proof that an
LCO has settled. Increase `duration` if the window amplitudes are still changing.
"""
function run_pazy_LCO_diagnostic(;
    airspeeds=(45.0,), duration=5.0, density=1.225, angle_of_attack=deg2rad(5.0),
    chordwise_panels=10, spanwise_panels=20,
    reference_speed=52.0, minimum_sweep_speed=35.0, time_step_chords=0.25,
    wake_length_chords=5.0, pulse_amplitude=0.1, pulse_start=0.1,
    pulse_duration=0.04, window_duration=0.25,
    gauge_y=0.0, gauge_z=2.25e-3/2,
    coupling_maximum_iterations=20, coupling_relaxation=0.5,
    coupling_geometry_tolerance=1e-4, coupling_load_tolerance=1e-4,
    progress_frequency=100, equilibrium_options=(;),
    output_directory=joinpath(@__DIR__, "output", "LCO_constant_speed_diagnostic"))

    @assert !isempty(airspeeds) && all(u -> isfinite(u) && u > 0, airspeeds)
    @assert all(isfinite, (duration, reference_speed, minimum_sweep_speed,
        time_step_chords, wake_length_chords, pulse_amplitude, pulse_start,
        pulse_duration, window_duration, gauge_y, gauge_z))
    @assert duration > 0 && reference_speed >= maximum(airspeeds)
    @assert minimum_sweep_speed > 0 && time_step_chords > 0 && wake_length_chords > 0
    @assert 0 <= pulse_start && pulse_duration > 0 && window_duration > 0
    @assert pulse_start + pulse_duration + 2window_duration < duration
    @assert !iszero(pulse_amplitude)
    _, _, chord, _ = AeroBeams.geometrical_properties_Pazy()
    nominal_dt = time_step_chords * chord / reference_speed
    baseline_dt = 20.0 / ceil(Int, 20.0 / nominal_dt)
    baseline_wake_rows = ceil(Int, wake_length_chords * chord / (minimum_sweep_speed * baseline_dt))
    @assert pulse_duration >= 8baseline_dt "Resolve the pulse with at least eight baseline steps"
    @assert window_duration >= 8baseline_dt
    static_settings = merge((;
        maximum_aerodynamic_time=0.20, minimum_aerodynamic_time=0.05,
        aerodynamic_window_time=0.01, aerodynamic_force_tolerance=1e-4,
        aerodynamic_circulation_tolerance=1e-4, maximum_static_iterations=100,
        static_relaxation=0.3, geometry_tolerance=1e-5, load_tolerance=1e-4,
        consecutive_equilibrium_iterations=2, newton_maximum_iterations=50,
        newton_absolute_tolerance=1e-8, newton_relative_tolerance=1e-8,
        newton_display_iterations=false, progress_frequency), equilibrium_options)
    variants = ((name="baseline", dt_factor=1.0, wake_factor=1.0, force=pulse_amplitude),
        (name="unforced", dt_factor=1.0, wake_factor=1.0, force=0.0),
        (name="half_dt", dt_factor=0.5, wake_factor=1.0, force=pulse_amplitude),
        (name="double_wake", dt_factor=1.0, wake_factor=2.0, force=pulse_amplitude))
    mkpath(output_directory)
    settings_path = joinpath(output_directory, "diagnostic_settings.txt")
    open(settings_path, "w") do io
        println(io, "Constant-speed diagnostic; root incidence [deg]: ", rad2deg(angle_of_attack))
        println(io, "Speeds [m/s]: $airspeeds; requested duration [s]: $duration")
        println(io, "Baseline dt [s]: $baseline_dt; baseline nominal wake [chords at minimum sweep speed]: $wake_length_chords")
        println(io, "Physical pulse: F_A1 = amplitude*sin(2pi*phase)*sinpi(phase)^2, 0<phase<1")
        println(io, "Pulse coefficient [N]: $pulse_amplitude; start [s]: $pulse_start; duration [s]: $pulse_duration")
        println(io, "Unforced is zero pulse; same copied equilibrium as baseline.")
        println(io, "Window duration [s]: $window_duration; gauge y/z [m]: $gauge_y / $gauge_z")
        println(io, "Mesh: $chordwise_panels x $spanwise_panels; gravity zero; first-element midpoint strain.")
        println(io, "Strong coupling: max iterations=$coupling_maximum_iterations; relaxation=$coupling_relaxation; geometry tolerance=$coupling_geometry_tolerance; load tolerance=$coupling_load_tolerance")
        println(io, "Static controls: ", static_settings)
    end
    summary_header = ["airspeed_mps", "variant", "dt_s", "wake_rows", "actual_wake_chords",
        "duration_s", "pulse_coefficient_N", "static_strain_ue", "last_window_amplitude_ue",
        "last_window_rms_ue", "tail_log_rms_slope_per_s", "max_geometry_residual",
        "max_load_residual", "max_coupling_iterations", "elapsed_wall_s"]
    summary_rows = Vector{Any}[]
    cases = []
    for speed in airspeeds
        equilibria = Dict{Tuple{Float64,Int},Any}()
        speed_tag = replace(string(Float64(speed)), "." => "p")
        for variant in variants
            started = time_ns()
            dt = baseline_dt * variant.dt_factor
            # Refine the already rounded baseline capacity to preserve its
            # physical wake age exactly, including for noninteger chord ratios.
            wake_rows = ceil(Int, baseline_wake_rows * variant.wake_factor / variant.dt_factor)
            case_name = "U$(speed_tag)_$(variant.name)"
            case_directory = joinpath(output_directory, case_name)
            mkpath(case_directory)
            println("\nDIAGNOSTIC $case_name: dt=$dt, wake rows=$wake_rows, pulse=$(variant.force) N")
            flush(stdout)
            key = (dt, wake_rows)
            equilibrium = get!(equilibria, key) do
                run_pazy_wing_uvlm_static(; static_settings..., airspeed=speed,
                    density, angle_of_attack, chordwise_panels, spanwise_panels,
                    symmetric_wing=true, maximum_wake_rows=wake_rows,
                    time_step_chords=dt * speed / chord, core_radius=1e-3,
                    wake_shedding_fraction=0.1, aerodynamic_interaction=true)
            end
            writedlm(joinpath(case_directory, "initialization.csv"),
                vcat(permutedims(equilibrium.history_columns), equilibrium.history), ',')
            open(settings_path, "a") do io
                println(io, "$case_name initialization: $(equilibrium.termination_reason); converged=$(equilibrium.converged); r_geometry=$(equilibrium.geometry_residual); r_load=$(equilibrium.load_residual)")
            end
            equilibrium.converged || error("Cannot start $case_name: $(equilibrium.termination_reason)")
            result = run_pazy_wing_uvlm(; airspeed=speed,
                initial_equilibrium=equilibrium, density, angle_of_attack, sideslip=0.0,
                initial_airspeed_fraction=1.0, airspeed_ramp_duration=0.0, duration,
                settling_time=pulse_start, perturbation_amplitude=variant.force,
                perturbation_duration=pulse_duration, chordwise_panels, spanwise_panels,
                symmetric_wing=true, time_step_chords=dt * speed / chord,
                maximum_wake_rows=wake_rows, core_radius=1e-3,
                wake_shedding_fraction=0.1, aerodynamic_interaction=true,
                save_uvlm_history=false, uvlm_save_frequency=1,
                newton_maximum_iterations=20, newton_absolute_tolerance=1e-5,
                newton_relative_tolerance=1e-5, newton_display_iterations=false,
                newton_always_update_jacobian=true,
                coupling_scheme=:strong, coupling_maximum_iterations,
                coupling_relaxation, coupling_geometry_tolerance, coupling_load_tolerance,
                coupling_display_iterations=false, animation_frames=2,
                save_animation_history=false, save_structural_history=false,
                save_aerodynamic_load_history=false, progress_frequency)
            strain = 1e6 .* (result.root_axial_strain .+
                gauge_z .* result.root_curvature_2 .- gauge_y .* result.root_curvature_3)
            phase = (result.time .- pulse_start) ./ pulse_duration
            force = [0 < p < 1 ? variant.force * sin(2pi*p) * sinpi(p)^2 : 0.0 for p in phase]
            history_header = ["time_s", "airspeed_mps", "root_longitudinal_strain_ue",
                "tip_out_of_plane_m", "tip_twist_deg", "tip_force_A1_N",
                "root_axial_strain", "root_curvature_2_per_m", "root_curvature_3_per_m",
                "coupling_geometry_residual", "coupling_load_residual", "coupling_iterations"]
            history = hcat(result.time, result.airspeed_history, strain,
                result.tip_out_of_plane, result.tip_twist_degrees, force,
                result.root_axial_strain, result.root_curvature_2, result.root_curvature_3,
                result.coupling_geometry_residual, result.coupling_load_residual, result.coupling_iterations)
            writedlm(joinpath(case_directory, "history.csv"),
                vcat(permutedims(history_header), history), ',')
            windows = pazy_diagnostic_windows(result.time, strain;
                pulse_end=pulse_start+pulse_duration, window_duration)
            window_header = ["window_start_s", "window_end_s", "mean_strain_ue",
                "minimum_strain_ue", "maximum_strain_ue", "raw_half_range_ue",
                "detrended_half_range_ue", "detrended_rms_ue"]
            writedlm(joinpath(case_directory, "windows.csv"),
                vcat(permutedims(window_header), windows), ',')
            tail = max(1, size(windows,1)-3):size(windows,1)
            centers = (windows[tail,1] + windows[tail,2]) / 2
            centered_time = centers .- mean(centers)
            log_rms = log.(max.(windows[tail,8], eps(Float64)))
            tail_slope = dot(centered_time, log_rms .- mean(log_rms)) / dot(centered_time, centered_time)
            elapsed = (time_ns()-started)/1e9
            row = Any[speed, variant.name, dt, wake_rows, wake_rows*speed*dt/chord,
                last(result.time), variant.force, first(strain), windows[end,7], windows[end,8],
                tail_slope, maximum(result.coupling_geometry_residual),
                maximum(result.coupling_load_residual), maximum(result.coupling_iterations), elapsed]
            push!(summary_rows, row)
            writedlm(joinpath(output_directory, "diagnostic_summary.csv"),
                vcat(permutedims(summary_header), permutedims(hcat(summary_rows...))), ',')
            push!(cases, (; name=case_name, speed, variant=variant.name, time=result.time,
                strain, force, windows, dt, wake_rows, initial_strain=first(strain), tail_slope,
                case_directory))
            open(settings_path, "a") do io
                println(io, "$case_name completed; actual duration=$(last(result.time)) s; wall time=$elapsed s; tail log RMS slope=$tail_slope /s")
            end
            @printf("Completed %s: last-window amplitude %.6g microstrain; tail log RMS slope %.6g /s; wall %.1f s\n",
                case_name, windows[end,7], tail_slope, elapsed)
            flush(stdout)
            # Save comparisons after each completed case, even during a long suite.
            pazy_plot_LCO_diagnostic(cases, output_directory)
        end
    end
    return (; cases, summary=summary_rows, output_directory)
end

function pazy_plot_LCO_diagnostic(cases, output_directory)
    traces = Plots.plot(; xlabel="Time [s]", ylabel="Strain change from equilibrium [microstrain]",
        legend=:outerright, size=(1100, 500))
    amplitudes = Plots.plot(; xlabel="Window center time [s]",
        ylabel="Detrended half range [microstrain]", legend=:outerright, size=(1100, 500))
    for case in cases
        label = "$(case.speed) m/s, $(case.variant)"
        Plots.plot!(traces, case.time, case.strain .- case.initial_strain; label, linewidth=0.6)
        centers = (case.windows[:,1] + case.windows[:,2]) / 2
        Plots.plot!(amplitudes, centers, case.windows[:,7]; label, marker=:circle, linewidth=1.2)
    end
    Plots.savefig(traces, joinpath(output_directory, "strain_histories.png"))
    Plots.savefig(amplitudes, joinpath(output_directory, "window_amplitudes.png"))
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    # Small dense solves can suffer from thread overhead. This affects this process only.
    BLAS.set_num_threads(1)
    run_pazy_LCO_diagnostic()
end
