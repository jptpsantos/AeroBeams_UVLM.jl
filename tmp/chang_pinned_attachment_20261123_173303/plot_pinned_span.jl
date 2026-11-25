using Printf

root = normpath(joinpath(@__DIR__,"..",".."))
convergence = joinpath(root,"lib","WingPropellerUVLM","examples",
    "chang_linear_aeroelastic","studies","convergence")
include(joinpath(convergence,"plot_chang_convergence_sensitivity.jl"))

originalPath,pinnedPath,figurePath = ARGS[1:3]
fullSweep = ARGS[4:end] == ["--full"]
original = Saved.load_current_results(originalPath;metrics=WING_METRICS)
pinned = Saved.load_current_results(pinnedPath;metrics=WING_METRICS)
panelCounts = pinned.settings["sweeps"]["wing_span"]

for key in ("physical","structural","integration","coupling","response","damping",
    "nominal","core_mode","damping_method_sha256","propeller_moment_projection")
    original.settings[key] == pinned.settings[key] || error("Study settings differ for $key")
end

originalSweep = filter(r -> r.family == "wing_span",original.rows)
originalRows = filter(r -> r.x in panelCounts,originalSweep)
pinnedRows = filter(r -> r.family == "wing_span",pinned.rows)
for rows in (originalRows,pinnedRows)
    [r.x for r in rows] == panelCounts || error("Expected spanwise panel counts $(panelCounts)")
    all(r -> r.valid && all(isfinite,r.values),rows) || error("Incomplete or invalid damping fits")
end

for (before,after) in zip(originalRows,pinnedRows)
    @printf("%2d panels | original bending % .8f, torsion % .8f | pinned bending % .8f, torsion % .8f\n",
        after.x,before.values[1:2]...,after.values[1:2]...)
end

for i in 2:length(panelCounts)
    before,after = originalRows[i],pinnedRows[i]
    originalDifference = abs.(before.values .- originalRows[i-1].values)
    @printf("%2d to %2d | original %.8f, %.8f | pinned %.8f, %.8f | pinned percent %.5f, %.5f\n",
        panelCounts[i-1],panelCounts[i],originalDifference[1:2]...,after.previous[1:2]...,
        after.successive_percent[1:2]...)
end

p = sensitivity_panel(pinnedRows,"wing_span","Spanwise panels",:percent,WING_SIGNALS;showLegend=false)
plot!(p;ylabel="Difference (%)")
if panelCounts == [20,40,80]
    plot!(p;xlims=(38,82),xticks=[40,60,80])
else
    plot!(p;xlims=(0,last(panelCounts)+3))
end
if fullSweep
    all(r -> r.valid && all(isfinite,r.values),originalSweep) ||
        error("Incomplete or invalid damping fits in the original sweep")
    originalPlot = sensitivity_panel(originalSweep,"wing_span","Spanwise panels",:percent,
        WING_SIGNALS;showLegend=false)
    for panel in (originalPlot,p)
        plot!(panel;ylabel="Difference (%)",xlims=(0,83),xticks=0:10:80,
            margin=6Plots.mm)
    end
    p = plot(originalPlot,p;layout=(1,2),size=(1200,400),dpi=300)
end
savefig(p,figurePath)
