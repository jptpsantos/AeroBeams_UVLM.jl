using Plots
using DelimitedFiles

# Propeller-only study: one panel per refinement parameter, with CT and CQ
# together. No wing results are available in these input files.
const SENSITIVITY_STUDIES = (
    (file="results_sensitivity_chord.txt", title="(a) Chordwise refinement",
     xlabel="Chordwise panels", reverse=false),
    (file="results_sensitivity_span.txt", title="(b) Spanwise refinement",
     xlabel="Spanwise panels", reverse=false),
    (file="results_sensitivity_step_2.txt", title="(c) Azimuthal-step refinement",
     xlabel="Azimuthal step (deg.)", reverse=true),
    (file="results_sensitivity_wake.txt", title="(d) Wake-length refinement",
     xlabel="Wake length (propeller revolutions)", reverse=false),
)

function refinement_changes(filepath; reverse=false)
    values, header = readdlm(filepath, '\t', Float64, '\n'; header=true)
    columns = strip.(String.(vec(header)))
    ct_column = findfirst(==("CT_mean"), columns)
    cq_column = findfirst(==("CQ_mean"), columns)
    isnothing(ct_column) && error("Missing CT_mean column in $filepath")
    isnothing(cq_column) && error("Missing CQ_mean column in $filepath")
    size(values, 1) >= 2 || error("At least two refinements are needed: $filepath")

    # Always compare in the physical coarse-to-fine direction. Smaller
    # azimuthal steps are finer; larger panel counts/wake lengths are finer.
    order = sortperm(values[:, 1]; rev=reverse)
    x = values[order, 1]
    all(isfinite, x) || error("Non-finite refinement parameter in $filepath")
    allunique(x) || error("Duplicate refinement parameters in $filepath")

    function percent_change(column)
        q = values[order, column]
        all(isfinite, q) || error("Non-finite coefficient in $filepath")
        previous = q[1:end-1]
        any(iszero, previous) && error("Percentage change from zero is undefined: $filepath")
        # Magnitude of the change relative to the previous, coarser case:
        # 100 * abs(q_fine - q_coarse) / abs(q_coarse).
        # Recompute from means; stored Change columns may use another order.
        return 100 .* abs.(diff(q)) ./ abs.(previous)
    end

    # The first case has no predecessor: omit it instead of assigning 0%.
    return (x=x[2:end], ct=percent_change(ct_column), cq=percent_change(cq_column))
end

function plot_sensitivity(data_dir=@__DIR__; output_dir=data_dir, yscale=:log10)
    changes = [refinement_changes(joinpath(data_dir, study.file); reverse=study.reverse)
               for study in SENSITIVITY_STUDIES]

    # Share the same limits and scale across all four panels. Logarithmic
    # percentages keep small wake-length changes visible alongside mesh changes.
    # Use yscale=:identity for a shared linear axis instead.
    percentages = vcat([vcat(change.ct, change.cq) for change in changes]...)
    low, peak = extrema(percentages)
    if yscale == :log10
        low > 0 || error("Zero changes require a linear axis: use yscale=:identity")
        lower_power = floor(Int, log10(low))
        upper_power = max(lower_power + 1, ceil(Int, log10(peak)))
        limits = (10.0^lower_power, 10.0^upper_power)
    elseif yscale == :identity
        limits = (0.0, peak > 0 ? 1.10 * peak : 1.0)
    else
        error("Use yscale=:log10 or yscale=:identity")
    end
    panels = Plots.Plot[]
    for (study, change) in zip(SENSITIVITY_STUDIES, changes)
        p = Plots.plot(
            change.x, change.ct;
            label="Propeller CT", color=:black, linestyle=:solid,
            markershape=:circle, markersize=4, linewidth=2,
            title=study.title, titlefontsize=12,
            xlabel=study.xlabel,
            ylabel="Absolute relative change (%)",
            ylims=limits, yscale=yscale, xflip=study.reverse,
            tickfontsize=10, labelfontsize=11, legendfontsize=10,
            legend=:bottomleft, foreground_color_legend=nothing,
            grid=false, framestyle=:box, margin=5Plots.mm,
        )
        Plots.plot!(p, change.x, change.cq;
            label="Propeller CQ", color=:steelblue, linestyle=:dash,
            markershape=:rect, markersize=4, markercolor=:white,
            markerstrokecolor=:steelblue, linewidth=2,
        )
        push!(panels, p)
    end

    # Lines only connect sampled cases; spline smoothing can create artificial
    # extrema and obscure whether the actual refinements have converged.
    figure = Plots.plot(panels...; layout=(2, 2), link=:y,
                        size=(1100, 760), dpi=300)
    mkpath(output_dir)
    savefig(figure, joinpath(output_dir, "sensitivity_propeller_percent.png"))
    savefig(figure, joinpath(output_dir, "sensitivity_propeller_percent.pdf"))
    return figure
end

# Run directly, or include this script in the Julia IDE.
# Optional command-line arguments: input directory, output directory.
data_dir = isempty(ARGS) ? (@__DIR__) : abspath(ARGS[1])
output_dir = length(ARGS) < 2 ? data_dir : abspath(ARGS[2])
plot_final = plot_sensitivity(data_dir; output_dir=output_dir)
display(plot_final)
