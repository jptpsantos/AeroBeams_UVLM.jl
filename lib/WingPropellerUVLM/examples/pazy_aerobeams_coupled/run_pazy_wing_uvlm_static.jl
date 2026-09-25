# Single-condition STATIC Pazy analysis. Run this file from Julia or the IDE.
# Follows AeroBeams' PazyWingPitchRange.jl structural procedure: SteadyProblem,
# nonlinear Newton balance, and the same wingtip OOP/IP displacement and twist.
# UVLM replaces strip aerodynamics; no linearized UVLM or dynamic beam settling.
include(joinpath(@__DIR__, "PazyWingUVLMStatic.jl"))
using DelimitedFiles

# Flight condition. Incidence is imposed through the freestream, once only.
airspeed = 55.0
density = 1.225
angle_of_attack = deg2rad(2.0)

# UVLM controls. The structural mesh remains the original 15-element Pazy mesh.
chordwise_panels = 4
spanwise_panels = 15
symmetric_wing = true
wake_length_chords = 10.0
maximum_wake_rows = ceil(Int, wake_length_chords*chordwise_panels)
core_radius = 1e-3                    # metres
wake_shedding_fraction = 0.1
aerodynamic_interaction = true
# dt = chord/(chordwise_panels*airspeed), as in the dynamic run.

# Time marching at EACH fixed trial shape. A fresh free wake is developed each
# time, avoiding fictitious aerodynamic forces from static shape updates.
maximum_aerodynamic_time = 0.20       # seconds per trial shape, NOT total runtime
minimum_aerodynamic_time = 0.05       # also requires at least 2 retained wake ages
aerodynamic_window_time = 0.01        # force/circulation peak-to-peak check window
aerodynamic_force_tolerance = 1e-4
aerodynamic_circulation_tolerance = 1e-4

# Static outer iteration: relax the applied loads, solve the nonlinear beam,
# settle UVLM on that deformed shape, and check the unrelaxed load mismatch.
maximum_static_iterations = 100
static_relaxation = 0.3
geometry_tolerance = 1e-5            # maximum grid change / chord
load_tolerance = 1e-4               # same force/moment scaling as dynamic coupling
consecutive_equilibrium_iterations = 2

# Newton-Raphson controls (structural STATIC equilibrium, no Newmark step).
newton_maximum_iterations = 50
newton_absolute_tolerance = 1e-8
newton_relative_tolerance = 1e-8
newton_display_iterations = false
progress_frequency = 100             # UVLM wake steps; every static iteration prints

# Dedicated output folder: existing dynamic results are not overwritten.
output_directory = joinpath(@__DIR__, "output", "static_uvlm")
save_convergence_plot = true

static_result = run_pazy_wing_uvlm_static(;
    airspeed, density, angle_of_attack, chordwise_panels, spanwise_panels,
    symmetric_wing, maximum_wake_rows, core_radius, wake_shedding_fraction,
    aerodynamic_interaction, maximum_aerodynamic_time, minimum_aerodynamic_time,
    aerodynamic_window_time, aerodynamic_force_tolerance, aerodynamic_circulation_tolerance,
    maximum_static_iterations, static_relaxation, geometry_tolerance, load_tolerance,
    consecutive_equilibrium_iterations, newton_maximum_iterations, newton_absolute_tolerance,
    newton_relative_tolerance, newton_display_iterations, progress_frequency)

mkpath(output_directory)
writedlm(joinpath(output_directory,"static_convergence.csv"),
    vcat(permutedims(static_result.history_columns), static_result.history), ',')
# Always write the convergence flag: a time/iteration limit is not equilibrium.
writedlm(joinpath(output_directory,"static_summary.csv"), [
    "converged" "U_m_s" "root_incidence_deg" "tip_OOP_m" "tip_IP_m" "tip_twist_deg";
    static_result.converged airspeed rad2deg(angle_of_attack) static_result.tip_out_of_plane static_result.tip_in_plane static_result.tip_twist_degrees
], ',')
open(joinpath(output_directory,"static_status.txt"), "w") do io
    println(io, static_result.termination_reason)
    println(io, "Time in the convergence history is cumulative aerodynamic settling time, not physical wing motion.")
end
if save_convergence_plot
    static_plotter = AeroBeams.Plots
    t = static_result.history[:,2]
    p1 = static_plotter.plot(t, static_result.history[:,7], ylabel="Tip OOP [m]", label=false)
    p2 = static_plotter.plot(t, static_result.history[:,9], ylabel="Tip twist [deg]", label=false,
        xlabel="Cumulative aerodynamic settling time [s]")
    figure = static_plotter.plot(p1,p2,layout=(2,1),size=(900,650),
        plot_title=static_result.converged ? "Static equilibrium converged" : "NOT converged - last iterates")
    static_plotter.savefig(figure, joinpath(output_directory,"static_convergence.png"))
end
println("Tip OOP [m]: ",static_result.tip_out_of_plane)
println("Tip IP [m]: ",static_result.tip_in_plane)
println("Tip twist [deg]: ",static_result.tip_twist_degrees)
println("Converged: ",static_result.converged,"; results in ",output_directory)
static_result.converged || @warn static_result.termination_reason
