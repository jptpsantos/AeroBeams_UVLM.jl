# Included by run_pazy_LCO.jl after the coupling and plotting helpers.
include(joinpath(@__DIR__, "plot_pazy_LCO.jl"))
if !isdefined(@__MODULE__, :run_pazy_wing_uvlm_static)
    include(joinpath(@__DIR__, "PazyWingUVLMStatic.jl"))
end

"""One linear ramp, on its own local time axis starting at zero."""
function pazy_linear_airspeed(t; start_speed, end_speed, duration)
    return start_speed + (end_speed - start_speed) * clamp(t / duration, 0.0, 1.0)
end

"""Time for a linear ramp at a prescribed positive speed-change magnitude [m/s^2]."""
function pazy_lco_ramp_duration(start_speed, end_speed, rate)
    all(isfinite, (start_speed, end_speed, rate)) &&
        start_speed > 0 && end_speed > 0 && start_speed != end_speed && rate > 0 ||
        throw(ArgumentError("Ramp endpoints must be positive, finite and distinct; rate must be positive and finite"))
    return abs(end_speed-start_speed)/rate
end

"""Align simple fractional ramp durations without reducing dt to rounding noise."""
function pazy_lco_aligned_timestep(upsweep_duration, downsweep_duration, nominal_dt;
    fixed_time_step=nothing)
    if !isnothing(fixed_time_step)
        dt = Float64(fixed_time_step)
        isfinite(dt) && 0 < dt <= nominal_dt * (1+1e-12) || throw(ArgumentError(
            "fixed_time_step must be positive and no larger than the nominal timestep"))
        all(duration -> isapprox(round(Int,duration/dt)*dt,duration; atol=1e-8,rtol=0),
            (upsweep_duration,downsweep_duration)) || throw(ArgumentError(
            "fixed_time_step must align both ramp endpoints"))
        return dt
    end
    durations = rationalize.(Float64[upsweep_duration, downsweep_duration]; tol=1e-8)
    unit_duration = Float64(gcd(durations...))
    dt = unit_duration / ceil(Int, unit_duration / nominal_dt)
    dt >= nominal_dt / 2 || throw(ArgumentError(
        "Ramp durations require excessive timestep reduction to align their endpoints; " *
        "choose durations with a simpler ratio"))
    return dt
end

pazy_root_strain_ue(result, y, z) = 1e6 .* (result.root_axial_strain .+
    z .* result.root_curvature_2 .- y .* result.root_curvature_3)

"""Run independent upsweep and downsweep simulations at fixed root incidence.

Defaults: 35 -> 52 m/s in 40 s and 52 -> 35 m/s in 20 s, at 5 degrees.
With durations omitted, `upsweep_rate=0.425` and `downsweep_rate=0.85` [m/s^2]
set positive rate magnitudes. Durations follow each branch's speed range, so
changing endpoints preserves dU/dt = +0.425 and -0.85 m/s^2, respectively.
An explicitly supplied duration overrides the rate for that branch, preserving
duration-based calls used for numerical diagnostics and earlier comparisons.
`downsweep_initial_speed` defaults to `peak_speed`; set it independently to
start the downsweep higher without extending the upsweep. The shared timestep
is sized for the highest speed of either branch.
Each branch starts from its own convergence-checked coupled static equilibrium
and developed wake. Starting structural velocities are zero. No upsweep state
is supplied to the downsweep, and each branch has its own clock starting at zero.
`equilibrium_options` overrides static-solver convergence controls; an
unconverged initialization stops that branch before dynamic integration.
Gauge offsets are local beam coordinates (+y forward, +z upper surface).
`fixed_time_step` optionally matches a comparison run's dt exactly.
`upsweep_pulse_amplitude` and `downsweep_pulse_amplitude` are coefficients [N]
of the smooth out-of-plane tip force along structural A1. Zero disables the
pulse in the corresponding branch.
Each branch's `*_pulse_start` and `*_pulse_duration` are local branch times
[s]. The bipolar waveform has zero net impulse and a peak magnitude about
0.6495 times the coefficient. A nonzero pulse must fit inside its branch
and span at least eight timesteps. Both coefficients default to zero.
UVLM mesh, timestep, symmetry, interaction, core radius, and shedding controls
are shared by static initialization and dynamics. `maximum_wake_rows=nothing`
sizes wake capacity from `wake_length_chords`; an integer overrides that sizing.
`run_branches` selects :both (default), :upsweep, or :downsweep for simulation.
Skipped branches have no static solve or dynamic integration and return nothing.
Timestep and wake sizing retain both configured ramps for consistent comparisons.
`plot_branches` defaults to `run_branches` and independently selects saved
figures from the computed branches. No previous CSV is used in automatic plots.
By default, save to `output_base_directory/AoA<degrees>/run_001`, then run_002,
etc. An explicit `output_directory` must be empty; existing results are never
overwritten. The returned `output_directory` identifies this run's actual folder.
"""
function run_pazy_linear_LCO(; initial_speed=35.0, peak_speed=52.0,
    downsweep_initial_speed=peak_speed,
    upsweep_rate=0.425, downsweep_rate=0.85,
    upsweep_duration=nothing, downsweep_duration=nothing,
    density=1.225, angle_of_attack=deg2rad(5.0),
    chordwise_panels=10, spanwise_panels=20, time_step_chords=0.25,
    fixed_time_step=nothing, maximum_wake_rows=nothing,
    core_radius=1e-3, wake_shedding_fraction=0.1,
    symmetric_wing=true, aerodynamic_interaction=true,
    wake_length_chords=5.0, gauge_y=0.0, gauge_z=2.25e-3/2,
    upsweep_pulse_amplitude=0.0, upsweep_pulse_start=0.1,
    upsweep_pulse_duration=0.04,
    downsweep_pulse_amplitude=0.0, downsweep_pulse_start=0.1,
    downsweep_pulse_duration=0.04,
    coupling_scheme=:strong, coupling_maximum_iterations=20,
    coupling_relaxation=0.5, coupling_geometry_tolerance=1e-4,
    coupling_load_tolerance=1e-4, progress_frequency=100,
    equilibrium_options=(;), run_branches=:both, plot_branches=run_branches,
    output_base_directory=joinpath(@__DIR__, "output", "LCO_independent_sweeps"),
    output_directory=nothing)

    upsweep_duration = isnothing(upsweep_duration) ?
        pazy_lco_ramp_duration(initial_speed, peak_speed, upsweep_rate) : upsweep_duration
    downsweep_duration = isnothing(downsweep_duration) ?
        pazy_lco_ramp_duration(downsweep_initial_speed, initial_speed, downsweep_rate) : downsweep_duration
    @assert all(isfinite, (initial_speed, peak_speed, downsweep_initial_speed, upsweep_duration,
        downsweep_duration, time_step_chords, wake_length_chords, gauge_y, gauge_z,
        downsweep_pulse_amplitude, downsweep_pulse_start, downsweep_pulse_duration,
        upsweep_pulse_amplitude, upsweep_pulse_start, upsweep_pulse_duration,
        angle_of_attack))
    @assert 0 < initial_speed < peak_speed
    @assert downsweep_initial_speed > initial_speed
    @assert upsweep_duration > 0 && downsweep_duration > 0
    @assert time_step_chords > 0 && wake_length_chords > 0
    chordwise_panels isa Integer && chordwise_panels > 0 &&
        spanwise_panels isa Integer && spanwise_panels > 0 ||
        throw(ArgumentError("UVLM panel counts must be positive integers"))
    isnothing(maximum_wake_rows) ||
        (maximum_wake_rows isa Integer && maximum_wake_rows > 0) ||
        throw(ArgumentError("maximum_wake_rows must be nothing or a positive integer"))
    UVLM.create_UVLMSolver(; coreRadius=core_radius,
        wakeSheddingFraction=wake_shedding_fraction, interaction=aerodynamic_interaction)
    run_names = pazy_lco_branch_names(run_branches; option_name="run_branches")
    selected_names = pazy_lco_branch_names(plot_branches)
    any(name -> name in selected_names, run_names) || throw(ArgumentError(
        "plot_branches must include at least one branch selected by run_branches"))
    _, _, chord, _ = AeroBeams.geometrical_properties_Pazy()
    maximum_speed = max(peak_speed, downsweep_initial_speed)
    nominal_dt = time_step_chords * chord / maximum_speed
    dt = pazy_lco_aligned_timestep(upsweep_duration, downsweep_duration, nominal_dt; fixed_time_step)
    pulses = Dict(
        "upsweep" => (; amplitude=upsweep_pulse_amplitude,
            start=upsweep_pulse_start, duration=upsweep_pulse_duration),
        "downsweep" => (; amplitude=downsweep_pulse_amplitude,
            start=downsweep_pulse_start, duration=downsweep_pulse_duration))
    for (name, duration) in (("upsweep", upsweep_duration), ("downsweep", downsweep_duration))
        pulse = pulses[name]
        name in run_names && !iszero(pulse.amplitude) || continue
        0 <= pulse.start && pulse.duration > 0 && pulse.start + pulse.duration <= duration ||
            throw(ArgumentError("The tip pulse must fit inside the $name"))
        pulse.duration >= 8dt * (1-1e-12) || throw(ArgumentError(
            "Resolve the $name tip pulse with at least eight timesteps"))
    end
    aligned_time_step_chords = dt * maximum_speed / chord
    maximum_wake_rows = isnothing(maximum_wake_rows) ?
        ceil(Int, wake_length_chords * chord / (initial_speed * dt)) : Int(maximum_wake_rows)
    static_settings = merge((;
        maximum_aerodynamic_time=0.20, minimum_aerodynamic_time=0.05,
        aerodynamic_window_time=0.01, aerodynamic_force_tolerance=1e-4,
        aerodynamic_circulation_tolerance=1e-4, maximum_static_iterations=100,
        static_relaxation=0.3, geometry_tolerance=1e-5, load_tolerance=1e-4,
        consecutive_equilibrium_iterations=2, newton_maximum_iterations=50,
        newton_absolute_tolerance=1e-8, newton_relative_tolerance=1e-8,
        newton_display_iterations=false, progress_frequency), equilibrium_options)
    output_directory = pazy_lco_prepare_output_directory(output_directory;
        base_directory=output_base_directory, angle_of_attack_deg=rad2deg(angle_of_attack))
    println("New LCO run at $(rad2deg(angle_of_attack)) degrees: $output_directory")
    settings_path = joinpath(output_directory, "pazy_linear_LCO_settings.txt")
    open(settings_path, "w") do io
        println(io, "Independent simulations; separate coupled static equilibrium and developed wake for each branch.")
        println(io, "Run branches: $run_branches; unselected branches have no static or dynamic solve")
        println(io, "Upsweep: $initial_speed -> $peak_speed m/s in $upsweep_duration s")
        println(io, "Downsweep: $downsweep_initial_speed -> $initial_speed m/s in $downsweep_duration s")
        println(io, "Upsweep rate [m/s^2]: ", (peak_speed-initial_speed)/upsweep_duration)
        println(io, "Downsweep rate [m/s^2]: ", (initial_speed-downsweep_initial_speed)/downsweep_duration)
        println(io, "Timestep reference airspeed [m/s]: $maximum_speed")
        println(io, "Each branch starts at t=0 with zero structural velocities; no preparation hold or state transfer between branches.")
        println(io, "Upsweep tip pulse coefficient [N]: $upsweep_pulse_amplitude; start [local s]: $upsweep_pulse_start; duration [s]: $upsweep_pulse_duration")
        println(io, "Downsweep tip pulse coefficient [N]: $downsweep_pulse_amplitude; start [local s]: $downsweep_pulse_start; duration [s]: $downsweep_pulse_duration")
        println(io, "Tip force along structural A1: coefficient*sin(2pi*phase)*sinpi(phase)^2 for 0<phase<1; zero otherwise. Zero net impulse; peak magnitude approximately 0.6495*abs(coefficient).")
        println(io, "Root incidence [deg]: ", rad2deg(angle_of_attack))
        println(io, "dt [s]: $dt; gauge y / z [m]: $gauge_y / $gauge_z")
        println(io, "UVLM mesh: $chordwise_panels x $spanwise_panels; symmetric wing: $symmetric_wing; aerodynamic interaction: $aerodynamic_interaction")
        println(io, "UVLM core radius [m]: $core_radius; wake shedding fraction: $wake_shedding_fraction")
        println(io, "Wake capacity [rows]: $maximum_wake_rows; retained age at capacity [s]: ", maximum_wake_rows*dt)
        println(io, "Nominal retained wake length [chords] at initial/maximum speed: ",
            initial_speed*maximum_wake_rows*dt/chord, " / ", maximum_speed*maximum_wake_rows*dt/chord)
        println(io, "Plot branches: $plot_branches; automatic plots use only branches computed in this run")
        println(io, "First-element midpoint gauge approximation; gravity is zero in the present coupling model.")
        println(io, "Static initialization controls: ", static_settings)
    end
    header = ["time_s", "airspeed_mps", "root_longitudinal_strain_ue",
        "root_axial_strain", "root_curvature_2_per_m", "root_curvature_3_per_m",
        "tip_out_of_plane_m", "tip_twist_deg", "wingtip_aoa_deg", "tip_force_A1_N"]
    branches = []
    for (name, start_speed, end_speed, duration) in (
        ("upsweep", initial_speed, peak_speed, upsweep_duration),
        ("downsweep", downsweep_initial_speed, initial_speed, downsweep_duration))
        name in run_names || continue
        println("\nInitializing independent $name at $start_speed m/s")
        equilibrium = run_pazy_wing_uvlm_static(;
            static_settings..., airspeed=start_speed, density, angle_of_attack,
            chordwise_panels, spanwise_panels, symmetric_wing, maximum_wake_rows,
            time_step_chords=dt * start_speed / chord,
            core_radius, wake_shedding_fraction, aerodynamic_interaction)
        writedlm(joinpath(output_directory, "pazy_linear_LCO_$(name)_initialization.csv"),
            vcat(permutedims(equilibrium.history_columns), equilibrium.history), ',')
        open(settings_path, "a") do io
            println(io, "$name initialization: ", equilibrium.termination_reason,
                "; converged=", equilibrium.converged, "; r_geometry=", equilibrium.geometry_residual,
                "; r_load=", equilibrium.load_residual)
        end
        equilibrium.converged || error("Cannot start $name: $(equilibrium.termination_reason)")
        airspeed_profile(t) = pazy_linear_airspeed(t; start_speed, end_speed, duration)
        pulse = pulses[name]
        pulse_amplitude = pulse.amplitude
        pulse_start = iszero(pulse_amplitude) ? 0.0 : pulse.start
        pulse_duration = iszero(pulse_amplitude) ? 0.1 : pulse.duration
        println("Running independent $name: $start_speed -> $end_speed m/s in $duration s")
        result = run_pazy_wing_uvlm(;
            airspeed=maximum_speed, airspeed_profile, initial_equilibrium=equilibrium,
            density, angle_of_attack, sideslip=0.0,
            initial_airspeed_fraction=1.0, airspeed_ramp_duration=0.0, duration,
            settling_time=pulse_start, perturbation_amplitude=pulse_amplitude,
            perturbation_duration=pulse_duration,
            chordwise_panels, spanwise_panels, symmetric_wing,
            time_step_chords=aligned_time_step_chords, maximum_wake_rows,
            core_radius, wake_shedding_fraction, aerodynamic_interaction,
            save_uvlm_history=false, uvlm_save_frequency=1,
            newton_maximum_iterations=20, newton_absolute_tolerance=1e-5,
            newton_relative_tolerance=1e-5, newton_display_iterations=false,
            newton_always_update_jacobian=true,
            coupling_scheme, coupling_maximum_iterations, coupling_relaxation,
            coupling_geometry_tolerance, coupling_load_tolerance,
            coupling_display_iterations=false, animation_frames=2,
            save_animation_history=false, save_structural_history=false,
            save_aerodynamic_load_history=false, progress_frequency)
        root_strain_ue = pazy_root_strain_ue(result, gauge_y, gauge_z)
        wingtip_aoa_degrees = rad2deg(angle_of_attack) .+ result.tip_twist_degrees
        tip_force_A1 = map(result.time) do t
            phase = (t-pulse_start)/pulse_duration
            0 < phase < 1 ? pulse_amplitude*sin(2pi*phase)*sinpi(phase)^2 : 0.0
        end
        response = hcat(result.time, result.airspeed_history, root_strain_ue,
            result.root_axial_strain, result.root_curvature_2, result.root_curvature_3,
            result.tip_out_of_plane, result.tip_twist_degrees, wingtip_aoa_degrees,
            tip_force_A1)
        writedlm(joinpath(output_directory, "pazy_linear_LCO_$(name).csv"),
            vcat(permutedims(header), response), ',')
        push!(branches, (; name, equilibrium, result, root_strain_ue, wingtip_aoa_degrees,
            tip_force_A1))
    end
    by_name = Dict(branch.name => branch for branch in branches)
    upsweep = get(by_name, "upsweep", nothing)
    downsweep = get(by_name, "downsweep", nothing)
    plotted = filter(branch -> branch.name in selected_names, branches)
    strain_plot = pazy_LCO_strain_plot(map(plotted) do branch
        (; airspeed=branch.result.airspeed_history, strain=branch.root_strain_ue,
            label=branch.name == "upsweep" ? "Upsweep" : "Downsweep",
            color=branch.name == "upsweep" ? :blue : :red)
    end)
    suffix = length(plotted) == 2 ? "" : "_$(only(plotted).name)"
    Plots.savefig(strain_plot, joinpath(output_directory, "pazy_linear_LCO_root_strain_vs_airspeed$suffix.png"))
    traces = map(plotted) do branch
        trace = Plots.plot(branch.result.time, branch.root_strain_ue;
            xlabel="Branch time [s]", ylabel="Root strain [microstrain]",
            title=branch.name, titlefontsize=11,
            color=branch.name == "upsweep" ? :blue : :red, linewidth=0.8,
            legend=false, right_margin=8Plots.mm)
        Plots.plot!(Plots.twinx(trace), branch.result.time, branch.result.airspeed_history;
            ylabel="Airspeed [m/s]", legend=false, color=:orangered, linewidth=1.2, grid=false)
        trace
    end
    trace_plot = Plots.plot(traces...; layout=(length(plotted),1), size=(800,325 * length(plotted)))
    Plots.savefig(trace_plot, joinpath(output_directory, "pazy_linear_LCO_root_strain_history$suffix.png"))
    println("Independent LCO sweep outputs saved in: $output_directory")
    return (; upsweep, downsweep, dt, output_directory, run_branches, plot_branches)
end
