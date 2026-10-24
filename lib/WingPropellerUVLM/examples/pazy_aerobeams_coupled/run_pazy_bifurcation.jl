# Independent constant-speed runs; edit this configuration before a full study.
if !isdefined(@__MODULE__, :PazyBifurcation)
    include(joinpath(@__DIR__, "PazyBifurcation.jl"))
end
using .PazyBifurcation

if abspath(PROGRAM_FILE) == @__FILE__
    bifurcation = (
        enabled=true,
        U_min=35.0, U_max=52.0, dU=1.0,
        simulation_time=10.0, transient_time=8.0,
        angle_of_attack=deg2rad(5.0),
        initialization=:equilibrium,
        pulse_amplitude=0.1, pulse_start=0.1, pulse_duration=0.04,
        responses=(:tip_out_of_plane_percent, :tip_in_plane_m, :tip_twist_deg, :root_strain_ue),
        check_convergence=true,
        save_time_histories=false, plot_results=true, plot_amplitude=true,
        zoom_limits=(40.0, 50.0),
    )
    bifurcation_results = run_pazy_bifurcation(; bifurcation...)
end
