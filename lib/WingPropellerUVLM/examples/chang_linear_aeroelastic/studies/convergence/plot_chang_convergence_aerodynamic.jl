# Plot saved results without loading a study entry or rerunning simulations.
# Run with --project=lib/WingPropellerUVLM/studies; an optional first argument
# selects a different results directory.
module ChangConvergenceAerodynamicPlots
using DelimitedFiles, Plots, Printf, TOML

const DEFAULT_DIRECTORY = normpath(joinpath(@__DIR__, "..", "..", "output", "convergence_aerodynamic"))
const FAMILY_LABELS = Dict(
    "wing_span" => ("Wing span mesh", "Spanwise panels"),
    "wing_chord" => ("Wing chord mesh", "Chordwise panels"),
    "prop_radial" => ("Propeller radial mesh", "Radial panels per blade"),
    "prop_chord" => ("Propeller chord mesh", "Chordwise panels per blade"),
    "wake_revolutions" => ("Wake length", "Retained wake (rev)"),
    "azimuth_deg" => ("Time step", "Azimuth increment (deg)"),
    "core" => ("Vortex core", "Core value (saved core law)"),
)
const COLORS = (:steelblue, :chocolate3)

number(x) = x isa Number ? Float64(x) : parse(Float64, string(x))
truth(x) = lowercase(string(x)) == "true"

function read_table(path)
    raw, header = readdlm(path, ','; header=true)
    columns = Dict(string(name) => i for (i, name) in enumerate(vec(header)))
    return [Dict(name => raw[i, column] for (name, column) in columns) for i in axes(raw, 1)]
end

function load_results(directory)
    saved = read_table(joinpath(directory, "convergence.csv"))
    cases = read_table(joinpath(directory, "case_matrix.csv"))
    matrix = Dict((string(c["family"]), Int(number(c["level"]))) => c for c in cases)
    rows = map(saved) do r
        family = string(r["family"])
        level = Int(number(r["level"]))
        c = matrix[(family, level)]
        string(r["label"]) == string(c["label"]) || error("Case label mismatch")
        string(r["status"]) == "completed" && truth(r["valid"]) ||
            error("Cannot plot an incomplete or invalid case: $(r["label"])")
        values = (number(r["CT"]), number(r["CQ"]))
        limits = (number(r["limit_CT"]), number(r["limit_CQ"]))
        all(isfinite, values) && all(x -> isfinite(x) && x > 0, limits) ||
            error("Invalid coefficients or tolerances: $(r["label"])")
        (; family, level, x=number(c[family]), values, limits,
            previous=(number(r["previous_change_CT"]), number(r["previous_change_CQ"])),
            accepted=truth(r["accepted"]), stable=truth(r["reference_stable"]),
            complete=truth(r["family_complete"]))
    end
    isempty(rows) && error("No saved convergence results in $directory")
    settings = TOML.parsefile(joinpath(directory, "study.toml"))["settings"]
    return rows, settings
end

function panel_options(family, group)
    title, xlabel = get(FAMILY_LABELS, family, (family, family))
    ticks = sort(unique(r.x for r in group))
    padding = max(0.04(last(ticks) - first(ticks)), 0.01abs(first(ticks)), 1e-8)
    return (; title, xlabel, xticks=(ticks, [@sprintf("%g", x) for x in ticks]),
        xlims=(first(ticks)-padding, last(ticks)+padding),
        xflip=family in ("azimuth_deg", "core"),
        titlefontsize=11, guidefontsize=10, tickfontsize=9,
        gridalpha=0.16, framestyle=:box, left_margin=7Plots.mm,
        bottom_margin=7Plots.mm, top_margin=4Plots.mm, right_margin=4Plots.mm)
end

function save_panels(panels, title, path)
    figure = plot(panels...; layout=(cld(length(panels), 2), 2),
        size=(1300, 340cld(length(panels), 2)), plot_title=title,
        plot_titlefontsize=14, background_color=:white, dpi=200)
    paths = [path * extension for extension in (".png", ".svg")]
    foreach(p -> savefig(figure, p), paths)
    return paths
end

function coefficient_panels(rows, metric)
    scale, ylabel = metric == 1 ? (1e3, "1000 C_T") : (1e4, "10000 C_Q")
    panels = Any[]
    for family in unique(r.family for r in rows)
        group = sort(filter(r -> r.family == family, rows); by=r -> r.level)
        values = [r.values[metric] * scale for r in group]
        # A minimum 2% full axis span keeps roundoff from appearing as a trend.
        lo, hi = extrema(values)
        padding = max(0.10(hi - lo), 0.01maximum(abs, values), 1e-8)
        p = plot([r.x for r in group], values; panel_options(family, group)...,
            ylabel, ylims=(lo-padding, hi+padding), color=COLORS[metric],
            marker=:circle, markersize=5, linewidth=2, legend=false,
            yformatter=x -> @sprintf("%.4g", x))
        push!(panels, p)
    end
    return panels
end

function error_panels(rows)
    panels = Any[]
    for family in unique(r.family for r in rows)
        group = sort(filter(r -> r.family == family, rows); by=r -> r.level)
        refined = filter(r -> r.level > 1, group)
        isempty(refined) && error("$family has no successive-level comparison")
        ratios = [[r.previous[i] / r.limits[i] for r in refined] for i in 1:2]
        all(x -> isfinite(x) && x >= 0, vcat(ratios...)) || error("Invalid errors in $family")
        # Exact zeros, if present, are displayed at machine epsilon on log axes.
        visible = [max.(v, eps(Float64)) for v in ratios]
        lo = 10.0^floor(log10(min(minimum(vcat(visible...)), 0.1))) / 1.5
        hi = 10.0^ceil(log10(max(maximum(vcat(visible...)), 10.0))) * 1.5
        stable = all(r -> r.stable && r.complete, group)
        title = first(FAMILY_LABELS[family]) * (stable ? " | stable" : " | unresolved")
        p = plot(; panel_options(family, group)..., title, ylabel="Successive error / tolerance",
            yscale=:log10, ylims=(lo, hi), legend=:bottomleft, legendfontsize=8)
        hspan!(p, [lo, 1.0]; color=:seagreen, alpha=0.07, label="")
        hline!(p, [1.0]; color=:gray35, linestyle=:dash, linewidth=1.3, label="Limit")
        for (i, label) in enumerate(("C_T", "C_Q"))
            plot!(p, [r.x for r in refined], visible[i]; label, color=COLORS[i],
                marker=i == 1 ? :circle : :diamond, markersize=5, linewidth=2)
        end
        push!(panels, p)
    end
    return panels
end

function plot_saved_convergence(directory=DEFAULT_DIRECTORY)
    directory = abspath(directory)
    rows, settings = load_results(directory)
    physical = settings["physical"]
    context = @sprintf("U = %g m/s, alpha = %g deg", physical["speed_mps"], physical["alpha_deg"])
    out = joinpath(directory, "plots")
    mkpath(out)
    paths = String[]
    for (metric, name, title) in ((1, "convergence_CT", "Mean thrust coefficient"),
            (2, "convergence_CQ", "Mean torque coefficient"))
        append!(paths, save_panels(coefficient_panels(rows, metric),
            "$title | $context", joinpath(out, name)))
    end
    append!(paths, save_panels(error_panels(rows),
        "Convergence: successive-level errors | tolerance = 1",
        joinpath(out, "convergence_errors")))
    open(joinpath(out, "plot_notes.md"), "w") do io
        println(io, "# Saved Chang aerodynamic convergence results\n")
        println(io, "Figures read convergence.csv, case_matrix.csv and study.toml; no simulations are run.")
        println(io, "Axes use saved case values, which can differ from the current entry script.")
        println(io, "Operating point: $context; RPM = $(physical["rpm"]).")
        println(io, "Core law: $(settings["core_mode"]); nominal value = $(settings["nominal"]["core"]).")
        println(io, "Simulated revolutions: $(settings["simulated_revolutions"]); averaged: $(settings["averaged_revolutions"]).\n")
        println(io, "CT/CQ figures show signed means. Their axes retain a minimum span to avoid magnifying roundoff.")
        println(io, "Successive error is max(mean difference, phase-waveform RMS difference), divided by the saved family tolerance.")
        println(io, "Tolerance is max(absolute tolerance, relative tolerance times the absolute finest-reference mean).")
        println(io, "The first level has no successive error. Exact zero errors are displayed at machine epsilon.")
        println(io, "Stable/unresolved follows the saved family_complete and reference_stable fields. It includes validity checks, not just plotted mean changes.\n")
        println(io, "| Family | Saved values, in refinement order | Accepted values |")
        println(io, "|---|---|---|")
        for family in unique(r.family for r in rows)
            group = sort(filter(r -> r.family == family, rows); by=r -> r.level)
            levels = join([@sprintf("%g", r.x) for r in group], ", ")
            accepted = join([@sprintf("%g", r.x) for r in group if r.accepted], ", ")
            println(io, "| $family | $levels | $(isempty(accepted) ? "none" : accepted) |")
        end
        if !physical["interaction"]
            println(io, "\nWing-propeller aerodynamic interaction is disabled. Flat CT/CQ under wing refinement does not establish convergence of wing loads.")
        end
        selection = joinpath(directory, "aerodynamic_selection.toml")
        if isfile(selection)
            println(io, "\nSaved selection status:\n\n```toml\n", read(selection, String), "```")
        end
    end
    foreach(println, paths)
    return paths
end
end

if abspath(PROGRAM_FILE) == @__FILE__
    ChangConvergenceAerodynamicPlots.plot_saved_convergence(isempty(ARGS) ?
        ChangConvergenceAerodynamicPlots.DEFAULT_DIRECTORY : only(ARGS))
end
