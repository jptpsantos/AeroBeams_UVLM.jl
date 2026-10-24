# Included by run_pazy_LCO.jl after the coupling and plotting helpers.
if !isdefined(@__MODULE__, :run_pazy_wing_uvlm_static)
    include(joinpath(@__DIR__, "PazyWingUVLMStatic.jl"))
end

"""One linear ramp, on its own local time axis starting at zero."""
function pazy_linear_airspeed(t; start_speed, end_speed, duration)
    return start_speed + (end_speed - start_speed) * clamp(t / duration, 0.0, 1.0)
end

pazy_root_strain_ue(result, y, z) = 1e6 .* (result.root_axial_strain .+
    z .* result.root_curvature_2 .- y .* result.root_curvature_3)

"""Run independent upsweep and downsweep simulations at fixed root incidence.

Defaults: 35 -> 52 m/s in 40 s and 52 -> 35 m/s in 20 s, at 5 degrees.
Each branch starts from its own convergence-checked coupled static equilibrium
and developed wake. Starting structural velocities are zero. No upsweep state
is supplied to the downsweep, and each branch has its own clock starting at zero.
`equilibrium_options` overrides static-solver convergence controls; an
unconverged initialization stops that branch before dynamic integration.
Gauge offsets are local beam coordinates (+y forward, +z upper surface).
"""
function run_pazy_linear_LCO(; initial_speed=35.0, peak_speed=52.0,
    upsweep_duration=40.0, downsweep_duration=20.0,
    density=1.225, angle_of_attack=deg2rad(5.0),
    chordwise_panels=10, spanwise_panels=20, time_step_chords=0.25,
    wake_length_chords=5.0, gauge_y=0.0, gauge_z=2.25e-3/2,
    coupling_scheme=:strong, coupling_maximum_iterations=20,
    coupling_relaxation=0.5, coupling_geometry_tolerance=1e-4,
    coupling_load_tolerance=1e-4, progress_frequency=100,
    equilibrium_options=(;),
    output_directory=joinpath(@__DIR__, "output", "LCO_independent_sweeps"))

    @assert all(isfinite, (initial_speed, peak_speed, upsweep_duration,
        downsweep_duration, time_step_chords, wake_length_chords, gauge_y, gauge_z))
    @assert 0 < initial_speed < peak_speed
    @assert upsweep_duration > 0 && downsweep_duration > 0
    @assert time_step_chords > 0 && wake_length_chords > 0
    _, _, chord, _ = AeroBeams.geometrical_properties_Pazy()
    nominal_dt = time_step_chords * chord / peak_speed
    unit_duration = AeroBeams.gcd_floats_approx(Float64[upsweep_duration, downsweep_duration])
    @assert unit_duration > 0
    dt = unit_duration / ceil(Int, unit_duration / nominal_dt)
    aligned_time_step_chords = dt * peak_speed / chord
    maximum_wake_rows = ceil(Int, wake_length_chords * chord / (initial_speed * dt))
    static_settings = merge((;
        maximum_aerodynamic_time=0.20, minimum_aerodynamic_time=0.05,
        aerodynamic_window_time=0.01, aerodynamic_force_tolerance=1e-4,
        aerodynamic_circulation_tolerance=1e-4, maximum_static_iterations=100,
        static_relaxation=0.3, geometry_tolerance=1e-5, load_tolerance=1e-4,
        consecutive_equilibrium_iterations=2, newton_maximum_iterations=50,
        newton_absolute_tolerance=1e-8, newton_relative_tolerance=1e-8,
        newton_display_iterations=false, progress_frequency), equilibrium_options)
    mkpath(output_directory)
    settings_path = joinpath(output_directory, "pazy_linear_LCO_settings.txt")
    open(settings_path, "w") do io
        println(io, "Independent simulations; separate coupled static equilibrium and developed wake for each branch.")
        println(io, "Upsweep: $initial_speed -> $peak_speed m/s in $upsweep_duration s")
        println(io, "Downsweep: $peak_speed -> $initial_speed m/s in $downsweep_duration s")
        println(io, "Each branch starts at t=0 with zero structural velocities; no preparation hold or state transfer between branches.")
        println(io, "Root incidence [deg]: ", rad2deg(angle_of_attack))
        println(io, "dt [s]: $dt; gauge y / z [m]: $gauge_y / $gauge_z")
        println(io, "First-element midpoint gauge approximation; gravity is zero in the present coupling model.")
        println(io, "Static initialization controls: ", static_settings)
    end
    header = ["time_s", "airspeed_mps", "root_longitudinal_strain_ue",
        "root_axial_strain", "root_curvature_2_per_m", "root_curvature_3_per_m",
        "tip_out_of_plane_m", "tip_twist_deg", "wingtip_aoa_deg"]
    branches = []
    for (name, start_speed, end_speed, duration) in (
        ("upsweep", initial_speed, peak_speed, upsweep_duration),
        ("downsweep", peak_speed, initial_speed, downsweep_duration))
        println("\nInitializing independent $name at $start_speed m/s")
        equilibrium = run_pazy_wing_uvlm_static(;
            static_settings..., airspeed=start_speed, density, angle_of_attack,
            chordwise_panels, spanwise_panels, symmetric_wing=true, maximum_wake_rows,
            time_step_chords=dt * start_speed / chord,
            core_radius=1e-3, wake_shedding_fraction=0.1, aerodynamic_interaction=true)
        writedlm(joinpath(output_directory, "pazy_linear_LCO_$(name)_initialization.csv"),
            vcat(permutedims(equilibrium.history_columns), equilibrium.history), ',')
        open(settings_path, "a") do io
            println(io, "$name initialization: ", equilibrium.termination_reason,
                "; converged=", equilibrium.converged, "; r_geometry=", equilibrium.geometry_residual,
                "; r_load=", equilibrium.load_residual)
        end
        equilibrium.converged || error("Cannot start $name: $(equilibrium.termination_reason)")
        airspeed_profile(t) = pazy_linear_airspeed(t; start_speed, end_speed, duration)
        println("Running independent $name: $start_speed -> $end_speed m/s in $duration s")
        result = run_pazy_wing_uvlm(;
            airspeed=peak_speed, airspeed_profile, initial_equilibrium=equilibrium,
            density, angle_of_attack, sideslip=0.0,
            initial_airspeed_fraction=1.0, airspeed_ramp_duration=0.0, duration,
            settling_time=0.0, perturbation_amplitude=0.0, perturbation_duration=0.1,
            chordwise_panels, spanwise_panels, symmetric_wing=true,
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
        root_strain_ue = pazy_root_strain_ue(result, gauge_y, gauge_z)
        wingtip_aoa_degrees = rad2deg(angle_of_attack) .+ result.tip_twist_degrees
        response = hcat(result.time, result.airspeed_history, root_strain_ue,
            result.root_axial_strain, result.root_curvature_2, result.root_curvature_3,
            result.tip_out_of_plane, result.tip_twist_degrees, wingtip_aoa_degrees)
        writedlm(joinpath(output_directory, "pazy_linear_LCO_$(name).csv"),
            vcat(permutedims(header), response), ',')
        push!(branches, (; name, equilibrium, result, root_strain_ue, wingtip_aoa_degrees))
    end
    upsweep, downsweep = branches
    strain_plot = Plots.plot(upsweep.result.airspeed_history, upsweep.root_strain_ue;
        color=:blue, linewidth=0.8, label="Upsweep ($upsweep_duration s)",
        xlabel="Airspeed [m/s]", ylabel="Root longitudinal strain [microstrain]",
        title="Root incidence = $(round(rad2deg(angle_of_attack); digits=3)) deg", titlefontsize=11)
    Plots.plot!(strain_plot, downsweep.result.airspeed_history, downsweep.root_strain_ue;
        color=:red, linewidth=0.8, label="Downsweep ($downsweep_duration s)")
    Plots.savefig(strain_plot, joinpath(output_directory, "pazy_linear_LCO_root_strain_vs_airspeed.png"))
    traces = map(branches) do branch
        trace = Plots.plot(branch.result.time, branch.root_strain_ue;
            xlabel="Branch time [s]", ylabel="Root strain [microstrain]",
            title=branch.name, titlefontsize=11, color=:blue, linewidth=0.8,
            legend=false, right_margin=8Plots.mm)
        Plots.plot!(Plots.twinx(trace), branch.result.time, branch.result.airspeed_history;
            ylabel="Airspeed [m/s]", legend=false, color=:orangered, linewidth=1.2, grid=false)
        trace
    end
    trace_plot = Plots.plot(traces...; layout=(2,1), size=(800,650))
    Plots.savefig(trace_plot, joinpath(output_directory, "pazy_linear_LCO_root_strain_history.png"))
    println("Independent LCO sweep outputs saved in: $output_directory")
    return (; upsweep, downsweep, dt, output_directory)
end
