# Independent constant-speed runs; edit this configuration before a full study.
if !isdefined(@__MODULE__, :PazyBifurcation)
    include(joinpath(@__DIR__, "PazyBifurcation.jl"))
end
if !isdefined(PazyBifurcation, :BIFURCATION_API_VERSION) || PazyBifurcation.BIFURCATION_API_VERSION < 2
    error("This Julia session has an older PazyBifurcation module loaded. " *
        "Restart the Julia REPL, or run this script in a fresh terminal process, " *
        "to load save_local_extrema, plot_paper and paper_plot_options support.")
end
using .PazyBifurcation

#if abspath(PROGRAM_FILE) == @__FILE__
    bifurcation = (
        enabled=true,
        U_min=40.0, U_max=50.0, dU=0.25,
        simulation_time=5.0, transient_time=4.0,
        angle_of_attack=deg2rad(5.0),
        initialization=:equilibrium,
        pulse_amplitude=1.0, pulse_start=0.1, pulse_duration=0.04,
        responses=(:tip_out_of_plane_percent, :root_strain_ue,
            :tip_twist_deg, :tip_geometric_incidence_deg),
        check_convergence=true,
        save_time_histories=false, save_local_extrema=true,
        plot_results=true, plot_amplitude=true, plot_paper=true,
        # Reference style, with your selected 40–50 m/s view (panels a–d).
        paper_plot_options=(xlims=(40.0, 50.0), xticks=40:5:50,
            ylims=((15.0, 35.0), (-5000.0, 0.0), (-15.0, 15.0), (-10.0, 20.0)),
            yticks=(15:5:35, -5000:1000:0, -15:5:15, -10:5:20),
            formats=("png", "pdf"), save_panels=true,
            fontfamily="Times", captionfontfamily="Times Bold",
            guidefontsize=12, tickfontsize=11, legendfontsize=10, captionfontsize=11,
            linewidth=2.0, markersize=1.35, markerstrokewidth=0.3,
            figure_size=(1100,860), panel_size=(550,420), dpi=300, gridalpha=0.18),
        zoom_limits=(40.0, 50.0),
    )
    if "--check-config" in ARGS
        # Keyword dispatch is checked, but enabled=false constructs no models or outputs.
        bifurcation_results = PazyBifurcation.run_pazy_bifurcation(; bifurcation..., enabled=false)
        println("Pazy bifurcation runner keywords accepted; no simulations were run.")
    else
        bifurcation_results = PazyBifurcation.run_pazy_bifurcation(; bifurcation...)
    end
#end
