using Plots
using DelimitedFiles

# Same definitions as plot_chang_convergence_aeroelastic.jl:
# successive (%) = 100 * abs(q_i - q_previous) / abs(q_i)
# reference  (%) = 100 * abs(q_i - q_finest) / abs(q_finest)
# Smaller azimuthal steps are finer; all other parameters increase with refinement.
function propeller_percentage_data(dataDir)
    studies = (
        (file="results_sensitivity_span.txt", caption="(a) Spanwise discretization",
         xlabel="Spanwise panels", xlims=(0, 40), xticks=0:10:40, reverse=false),
        (file="results_sensitivity_chord.txt", caption="(b) Chordwise discretization",
         xlabel="Chordwise panels", xlims=(0, 40), xticks=0:10:40, reverse=false),
        (file="results_sensitivity_step_2.txt", caption="(c) Azimuthal step",
         xlabel="Azimuthal step (deg.)", xlims=(0, 20), xticks=0:5:20, reverse=true),
        (file="results_sensitivity_wake.txt", caption="(d) Wake cutoff",
         xlabel="Wake Cutoff (rev.)", xlims=(1, 8), xticks=1:8, reverse=false),
    )
    return map(studies) do study
        values, header = readdlm(joinpath(dataDir, study.file), '\t', Float64, '\n'; header=true)
        names = strip.(String.(vec(header)))
        columns = [findfirst(==(name), names) for name in ("CT_mean", "CQ_mean")]
        any(isnothing, columns) && error("Missing coefficient column in $(study.file)")
        size(values, 1) >= 2 || error("Need two refinements in $(study.file)")
        order = sortperm(values[:, 1]; rev=study.reverse)
        x = values[order, 1]
        q = values[order, Int.(columns)]
        all(isfinite, x) && all(isfinite, q) || error("Non-finite data in $(study.file)")
        allunique(x) || error("Repeated refinement in $(study.file)")
        successive = fill(NaN, size(q))
        reference = fill(NaN, size(q))
        for j in axes(q, 2), i in axes(q, 1)
            if i > 1 && !iszero(q[i, j])
                successive[i, j] = 100 * abs(q[i, j] - q[i-1, j]) / abs(q[i, j])
            end
            if !iszero(q[end, j])
                reference[i, j] = 100 * abs(q[i, j] - q[end, j]) / abs(q[end, j])
            end
        end
        (; study..., x, q, successive, reference)
    end
end

function plot_propeller_percentages(data; comparison=:successive)
    comparison in (:successive, :reference) || error("Use :successive or :reference")
    # A single linear percentage range for all eight axes.
    percentages = [v for d in data for v in getproperty(d, comparison) if isfinite(v)]
    upperLimit = max(1.08maximum(percentages), 1e-3)
    groups = Plots.Plot[]
    for d in data
        y = getproperty(d, comparison)
        axesPlots = Plots.Plot[]
        for j in 1:2
            # Omit the first successive case, as in the aeroelastic plots.
            indices = comparison == :successive ? (2:length(d.x)) : eachindex(d.x)
            # Preserve increasing x axes in the original figure, including step.
            indices = sort(collect(indices); by=i -> d.x[i])
            p = plot(d.x[indices], y[indices, j];
                color=j == 1 ? :black : "#EF4444",
                linestyle=j == 1 ? :solid : :dash, linewidth=1.1,
                label="", legend=false, grid=false, markershape=:none,
                framestyle=:axes, mirror=false, tick_direction=:in,
                fontfamily="Times", guidefontsize=9, tickfontsize=8,
                foreground_color_axis=:black, foreground_color_border=:black,
                foreground_color_text=:black, background_color=:white,
                ylabel=j == 1 ? "Thrust diff. (%)" : "Torque diff. (%)",
                xlabel=j == 2 ? d.xlabel : "", xlims=d.xlims, xticks=d.xticks,
                ylims=(0, upperLimit), yscale=:identity, widen=false,
                left_margin=6Plots.mm, right_margin=4Plots.mm,
                top_margin=2Plots.mm, bottom_margin=j == 2 ? 3Plots.mm : 0Plots.mm,
            )
            push!(axesPlots, p)
        end
        # A small caption row keeps the original four groups of stacked plots.
        caption = plot(; framestyle=:none, legend=false, grid=false,
                       xlims=(0, 1), ylims=(0, 1), margin=0Plots.mm)
        annotate!(caption, 0.5, 0.45, text(d.caption, 10, "Times", :center))
        groupLayout = grid(3, 1; heights=[0.42, 0.45, 0.13])
        push!(groups, plot(axesPlots..., caption; layout=groupLayout))
    end
    return plot(groups...; layout=(2, 2), size=(1000, 700), dpi=300,
                background_color=:white)
end

function save_propeller_percentages(dataDir=@__DIR__; outputDir=dataDir)
    data = propeller_percentage_data(dataDir)
    mkpath(outputDir)
    figures = map((:successive, :reference)) do comparison
        figure = plot_propeller_percentages(data; comparison)
        stem = "sensitivity_propeller_$(comparison)_percent"
        for extension in (".png", ".pdf", ".svg")
            savefig(figure, joinpath(outputDir, stem * extension))
        end
        figure
    end
    # Keep the previous output name pointing to the corrected successive plot.
    for extension in (".png", ".pdf")
        savefig(figures[1], joinpath(outputDir, "sensitivity_propeller_percent" * extension))
    end
    # Export the plotted values so both definitions can be checked directly.
    open(joinpath(outputDir, "sensitivity_propeller_percentage_differences.csv"), "w") do io
        println(io, "study,parameter,CT_mean,CQ_mean,successive_CT_percent,successive_CQ_percent,reference_CT_percent,reference_CQ_percent")
        for d in data, i in eachindex(d.x)
            println(io, join((d.file, d.x[i], d.q[i, 1], d.q[i, 2],
                d.successive[i, 1], d.successive[i, 2],
                d.reference[i, 1], d.reference[i, 2]), ','))
        end
    end
    return figures
end

dataDir = isempty(ARGS) ? (@__DIR__) : abspath(ARGS[1])
outputDir = length(ARGS) < 2 ? dataDir : abspath(ARGS[2])
plot_successive, plot_reference = save_propeller_percentages(dataDir; outputDir)
display(plot_successive)
display(plot_reference)
