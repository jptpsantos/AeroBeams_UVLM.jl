ENV["GKSwstype"] = "100"
include(joinpath(@__DIR__, "..", "lib", "WingPropellerUVLM", "examples",
    "pazy_aerobeams_coupled", "PazyBifurcationPlots.jl"))
using .PazyBifurcationPlots
summary = joinpath(@__DIR__, "..", "lib", "WingPropellerUVLM", "examples",
    "pazy_aerobeams_coupled", "output", "bifurcation", "20261024_165844546", "bifurcation_summary.csv")
paths = plot_pazy_bifurcation(summary; xlims=(40.,50.), xticks=40:5:50,
    linewidth=0.55, markersize=1.35, save_panels=false,
    output_directory=joinpath(@__DIR__, "pazy_reference_style_controls_final"))
foreach(println, paths)
