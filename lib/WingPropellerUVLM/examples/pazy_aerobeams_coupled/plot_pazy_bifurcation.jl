# CSV-only plotting: this file does not load or run either aeroelastic solver.
include(joinpath(@__DIR__, "..", "bifurcation", "TimeDomainBifurcation.jl"))
using .TimeDomainBifurcation: plot_bifurcation

if abspath(PROGRAM_FILE) == @__FILE__
    isempty(ARGS) && error("Usage: julia --project=. plot_pazy_bifurcation.jl path/to/bifurcation_summary.csv")
    plot_bifurcation(abspath(ARGS[1]))
end
