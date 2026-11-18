using Plots
using DelimitedFiles

# Same successive-difference definition as the aeroelastic convergence plots:
# difference_i (%) = 100 * abs(q_i - q_previous) / abs(q_i).
# Calculate differences at the sampled refinements BEFORE interpolation.
function propeller_percentage_data(dataDir)
    studies = (
        (name="span", file="results_sensitivity_span.txt",
         xlabel="Spanwise panels", xlims=(0, 40), xticks=0:10:40, reverse=false),
        (name="chord", file="results_sensitivity_chord.txt",
         xlabel="Chordwise panels", xlims=(0, 40), xticks=0:10:40, reverse=false),
        (name="step", file="results_sensitivity_step_2.txt",
         xlabel="Azimuthal step (deg.)", xlims=(0, 20), xticks=0:5:20, reverse=true),
        (name="wake", file="results_sensitivity_wake.txt",
         xlabel="Wake cutoff (rev.)", xlims=(1, 8), xticks=1:8, reverse=false),
    )
    return map(studies) do study
        values, header = readdlm(joinpath(dataDir, study.file), '\t', Float64, '\n'; header=true)
        names = strip.(String.(vec(header)))
        columns = [findfirst(==(name), names) for name in ("CT_mean", "CQ_mean")]
        any(isnothing, columns) && error("Missing coefficient column in $(study.file)")
        size(values, 1) >= 3 || error("Need three refinements in $(study.file)")
        # Smaller azimuthal steps are finer; panel counts and wake lengths increase.
        order = sortperm(values[:, 1]; rev=study.reverse)
        x = values[order, 1]
        q = values[order, Int.(columns)]
        all(isfinite, x) && all(isfinite, q) || error("Non-finite data in $(study.file)")
        allunique(x) || error("Repeated refinement in $(study.file)")
        any(iszero, q[2:end, :]) && error("Percentage difference from zero is undefined")
        differences = 100 .* abs.(diff(q; dims=1)) ./ abs.(q[2:end, :])
        # Omit the coarsest case, which has no previous refinement. Then restore
        # increasing x for the original plotting convention, including step.
        refinedX = x[2:end]
        plotOrder = sortperm(refinedX)
        (; study..., x=refinedX[plotOrder], differences=differences[plotOrder, :])
    end
end

# Shape-preserving cubic Hermite interpolation (PCHIP). This interpolates the
# supplied samples exactly, retains turning points, and avoids negative
# overshoot when interpolating nonnegative percentage differences.
function interpolate_percentages(x, y, xFine)
    n = length(x)
    n == length(y) && n >= 2 || error("Need matching interpolation samples")
    h = diff(x)
    all(>(0), h) || error("Interpolation coordinates must increase")
    first(xFine) >= first(x) && last(xFine) <= last(x) ||
        error("Interpolate only inside the sampled range")
    slopes = diff(y) ./ h
    derivatives = zeros(Float64, n)
    if n == 2
        derivatives .= slopes[1]
    else
        for i in 2:n-1
            if slopes[i-1] * slopes[i] > 0
                w1 = 2h[i] + h[i-1]
                w2 = h[i] + 2h[i-1]
                derivatives[i] = (w1 + w2) / (w1 / slopes[i-1] + w2 / slopes[i])
            end
        end
        function endpoint(h1, h2, slope1, slope2)
            d = ((2h1 + h2) * slope1 - h1 * slope2) / (h1 + h2)
            sign(d) != sign(slope1) && return 0.0
            if sign(slope1) != sign(slope2) && abs(d) > 3abs(slope1)
                return 3slope1
            end
            return d
        end
        derivatives[1] = endpoint(h[1], h[2], slopes[1], slopes[2])
        derivatives[end] = endpoint(h[end], h[end-1], slopes[end], slopes[end-1])
    end
    return map(xFine) do xi
        i = clamp(searchsortedlast(x, xi), 1, n-1)
        t = (xi - x[i]) / h[i]
        (2t^3 - 3t^2 + 1) * y[i] + (t^3 - 2t^2 + t) * h[i] * derivatives[i] +
        (-2t^3 + 3t^2) * y[i+1] + (t^3 - t^2) * h[i] * derivatives[i+1]
    end
end

function save_propeller_percentages(dataDir=@__DIR__; outputDir=dataDir)
    data = propeller_percentage_data(dataDir)
    mkpath(outputDir)
    figures = Plots.Plot[]
    for d in data
        xFine = range(first(d.x), last(d.x); length=600)
        ctFine = interpolate_percentages(d.x, d.differences[:, 1], xFine)
        cqFine = interpolate_percentages(d.x, d.differences[:, 2], xFine)
        # Fit each figure's linear axis to BOTH coefficient differences, with
        # headroom and round ticks so the wake changes remain visible.
        peak = max(maximum(ctFine), maximum(cqFine), 1e-3)
        rawStep = 1.08peak / 6
        magnitude = 10.0^floor(log10(rawStep))
        niceSteps = [1.0, 2.0, 2.5, 5.0, 10.0]
        tickStep = magnitude * niceSteps[argmin(abs.(niceSteps .- rawStep / magnitude))]
        upperLimit = tickStep * ceil(1.08peak / tickStep)
        p = plot(xFine, ctFine;
            label="Propeller thrust (CT)", color=:black, linestyle=:solid, linewidth=3,
            fontfamily="sans-serif", tickfontsize=10, legendfontsize=10,
            labelfontsize=12, margin=4Plots.mm, grid=false,
            legend=d.name == "span" ? :topright : false,
            foreground_color_legend=nothing,
            background_color_legend=nothing,
            framestyle=:axes, mirror=false, tick_direction=:in,
            xlabel=d.xlabel, ylabel="Difference between refinements (%)",
            xlims=d.xlims, xticks=d.xticks,
            yscale=:identity, ylims=(0, upperLimit), yticks=0:tickStep:upperLimit,
            size=(600, 400), dpi=300,
        )
        plot!(p, xFine, cqFine;
            label="Propeller torque (CQ)", color=:red, linestyle=:dash, linewidth=3)
        for extension in (".pdf", ".png")
            savefig(p, joinpath(outputDir, "sensitivity_$(d.name)" * extension))
        end
        push!(figures, p)
    end
    open(joinpath(outputDir, "sensitivity_percentage_differences.csv"), "w") do io
        println(io, "study,parameter,successive_CT_percent,successive_CQ_percent")
        for d in data, i in eachindex(d.x)
            println(io, join((d.name, d.x[i], d.differences[i, 1], d.differences[i, 2]), ','))
        end
    end
    return figures
end

dataDir = isempty(ARGS) ? (@__DIR__) : abspath(ARGS[1])
outputDir = length(ARGS) < 2 ? dataDir : abspath(ARGS[2])
figures = save_propeller_percentages(dataDir; outputDir)
foreach(display, figures)
