# Replot the completed 35–55 m/s sweep using saved CSVs only.
# Run from the repository root: julia --project=. tmp/replot_pazy_bifurcation_limits_20261103.jl
ENV["GKSwstype"] = "100"
example = joinpath(@__DIR__, "..", "lib", "WingPropellerUVLM", "examples", "pazy_aerobeams_coupled")
include(joinpath(example, "PazyBifurcationPlots.jl"))
using .PazyBifurcationPlots

summary = joinpath(example, "output", "bifurcation", "20261031_112128589", "bifurcation_summary.csv")
options = (
    xlims=(35.0, 55.0), xticks=35:5:55,
    ylims=((10.0, 35.0), (-5000.0, 0.0), (-15.0, 15.0), (-10.0, 20.0)),
    yticks=(10:5:35, -5000:1000:0, -15:5:15, -10:5:20),
    formats=("png", "pdf"), save_panels=true,
    fontfamily="Times", captionfontfamily="Times Bold",
    guidefontsize=12, tickfontsize=11, legendfontsize=10, captionfontsize=11,
    linewidth=3.0, markersize=2.0, markerstrokewidth=0.3,
    figure_size=(1100,860), panel_size=(550,420), dpi=300, gridalpha=0.18)

# Check that both the window envelopes and local extrema fit in every panel.
angle = PazyBifurcationPlots.saved_incidence(summary)
for path in (summary, joinpath(dirname(summary), "bifurcation_local_extrema.csv"))
    data, columns = PazyBifurcationPlots.read_table(path)
    @assert extrema(Float64.(data[:, columns["airspeed_mps"]])) == options.xlims
    for (i, response) in enumerate(PazyBifurcationPlots.PANEL_RESPONSES)
        group, offset = PazyBifurcationPlots.response_group(data, columns, response, angle)
        fields = path == summary ? ("minimum", "maximum", "static_equilibrium") : ("value",)
        for field in fields
            values = Float64.(group[:, columns[field]]) .+ offset
            @assert all(v -> options.ylims[i][1] <= v <= options.ylims[i][2], values)
        end
    end
end

paths = plot_pazy_bifurcation(summary; options...)
@assert length(paths) == 10 && all(p -> isfile(p) && filesize(p) > 1000, paths)
plots = PazyBifurcationPlots.Plots
figure = plots.current()
for (i, subplot) in enumerate((1, 2, 5, 6))
    @assert plots.xlims(figure[subplot]) == options.xlims
    @assert plots.ylims(figure[subplot]) == options.ylims[i]
end
@assert !isdefined(Main, :AeroBeams)
println("Verified: all saved extrema fit in the four panels; no solvers loaded.")
foreach(println, paths)
