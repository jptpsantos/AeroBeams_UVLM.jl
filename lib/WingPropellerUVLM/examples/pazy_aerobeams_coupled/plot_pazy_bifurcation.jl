# CSV-only plotting: this file does not load or run either aeroelastic solver.
if !isdefined(@__MODULE__, :TimeDomainBifurcation)
    include(joinpath(@__DIR__, "..", "bifurcation", "TimeDomainBifurcation.jl"))
end
using .TimeDomainBifurcation: plot_bifurcation
if !isdefined(@__MODULE__, :PazyBifurcationPlots)
    include(joinpath(@__DIR__, "PazyBifurcationPlots.jl"))
end
using .PazyBifurcationPlots: plot_pazy_bifurcation

if abspath(PROGRAM_FILE) == @__FILE__
    isempty(ARGS) && error("Usage: julia --project=. plot_pazy_bifurcation.jl path/to/bifurcation_summary.csv [other_summary.csv ...]")
    # Incidence comes from each run's settings.txt; accepts later 3°/7° overlays.
    paths = plot_pazy_bifurcation(abspath.(ARGS))
    foreach(println, paths)
end
