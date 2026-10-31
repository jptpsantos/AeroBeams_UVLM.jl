# A separate entry point; no existing simulation or plotting files are modified.
if !isdefined(@__MODULE__,:PazySteppedBifurcation)
    include(joinpath(@__DIR__,"PazySteppedBifurcation.jl"))
end
using .PazySteppedBifurcation

function pazy_stepped_configuration()
    return (
        U_min=40.0,U_max=50.0,dU=1.0,
        plateau_time=10.0,transient_time=8.0, # increase if convergence is unverified
        sweep_order=:up_down,              # also :down_up, :upsweep or :downsweep
        angle_of_attack=deg2rad(5.0),
        chordwise_panels=10,spanwise_panels=20,
        time_step_chords=0.25,fixed_time_step=nothing,wake_length_chords=5.0,
        pulse_amplitude=0.1,pulse_start=0.1,pulse_duration=0.04, # first plateau ONLY
        extrema_method=:quadratic,          # derivative-zero peak refinement; or :sampled
        check_convergence=true,on_unconverged=:warn, # :error saves then stops
        save_time_histories=false,save_boundary_states=false,
        plot_results=true,
        plot_options=(linewidth=2.8,markersize=3.4,formats=("png","pdf")),
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    all(a -> a in ("--check-config","--help"),ARGS) || error("Use --check-config or --help, or no arguments to run")
    if "--help" in ARGS
        println("Run: julia --project=. run_pazy_stepped_bifurcation.jl")
        println("--check-config validates the schedule without simulations or output files.")
        println("Edit pazy_stepped_configuration() in this NEW file to set the study.")
    else
        configuration=pazy_stepped_configuration()
        if "--check-config" in ARGS
            checked=run_pazy_stepped_bifurcation(;configuration...,dry_run=true)
            p=checked.plan
            println("Configuration accepted. $(length(p.plateaus)) plateaus, $(p.total_time) s total, dt=$(p.dt) s.")
            println("Speed order [m/s]: ",[x.airspeed_mps for x in p.plateaus])
            println("No simulations were run and no output files were created.")
        else
            stepped_bifurcation_results=run_pazy_stepped_bifurcation(;configuration...)
            println("Stepped bifurcation outputs: ",stepped_bifurcation_results.output_directory)
        end
    end
end
