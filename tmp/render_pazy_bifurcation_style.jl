# Replot saved results only; no physical simulations are launched.
include(joinpath(@__DIR__, "..", "lib", "WingPropellerUVLM", "examples",
    "pazy_aerobeams_coupled", "PazyBifurcationPlots.jl"))
using .PazyBifurcationPlots
synthetic_summary = joinpath(@__DIR__, "pazy_bifurcation_style_validation_590936394985200",
    "synthetic", "bifurcation_summary.csv")
plot_pazy_bifurcation(synthetic_summary; xlims=(40.,50.), save_panels=false)
existing_summary = joinpath(@__DIR__, "..", "lib", "WingPropellerUVLM", "examples",
    "pazy_aerobeams_coupled", "output", "bifurcation", "20261024_165844546", "bifurcation_summary.csv")
paths = plot_pazy_bifurcation(existing_summary;
    output_directory=joinpath(@__DIR__, "pazy_bifurcation_reference_5deg_final"),
    filename="pazy_bifurcation_reference_5deg", save_panels=false)
foreach(println, paths)
