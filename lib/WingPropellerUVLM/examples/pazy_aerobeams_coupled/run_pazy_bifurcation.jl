# Independent constant-speed runs; edit this configuration before a full study.
if !isdefined(@__MODULE__, :PazyBifurcation)
    include(joinpath(@__DIR__, "PazyBifurcation.jl"))
end
using .PazyBifurcation

#if abspath(PROGRAM_FILE) == @__FILE__
    bifurcation = (
        enabled=true,
        U_min=40.0, U_max=50.0, dU=0.25,
        simulation_time=10.0, transient_time=8.0,
        angle_of_attack=deg2rad(5.0),
        initialization=:equilibrium,
        pulse_amplitude=0.0, pulse_start=0.1, pulse_duration=0.04,
        responses=(:tip_out_of_plane_percent, :root_strain_ue,
            :tip_twist_deg, :tip_geometric_incidence_deg),
        check_convergence=true,
        save_time_histories=false, save_local_extrema=true,
        plot_results=true, plot_amplitude=true, plot_paper=true,
        # Axis limits and tick spacing from the reference figure (panels a–d).
        paper_plot_options=(xlims=(40.0, 50.0), xticks=40:5:50,
            ylims=((15.0, 35.0), (-5000.0, 0.0), (-15.0, 15.0), (-10.0, 20.0)),
            yticks=(15:5:35, -5000:1000:0, -15:5:15, -10:5:20),
            formats=("png", "pdf"), save_panels=true, linewidth = 2.0, markersize=4.0),
        zoom_limits=(40.0, 50.0),
    )
    bifurcation_results = run_pazy_bifurcation(; bifurcation...)
#end
