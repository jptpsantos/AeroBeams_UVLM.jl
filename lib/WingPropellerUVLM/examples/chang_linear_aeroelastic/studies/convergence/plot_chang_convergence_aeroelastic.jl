# Plot archived damping results without loading a study or rerunning simulations.
# julia --project=lib/WingPropellerUVLM/studies path/to/this_script.jl [results_directory]
module ChangConvergenceAeroelasticPlots
using DelimitedFiles, Plots, Printf, TOML

const DEFAULT_DIRECTORY = normpath(joinpath(@__DIR__, "..", "..", "output", "convergence_aeroelastic_explicit"))
const CASE_NAMES = ("First case", "Second case")
const FAMILIES = ("prop_radial", "prop_chord")
const CONTROLS = ("wing_span", "wing_chord", "prop_radial", "prop_chord",
    "wake_revolutions", "core", "azimuth_deg")
const METRICS = ("pitch_lambda_per_s", "yaw_lambda_per_s", "pitch_frequency_hz", "yaw_frequency_hz")
const COLORS = (:steelblue, :chocolate3)
const MARKERS = (:circle, :diamond)

number(x) = x isa Number ? Float64(x) : parse(Float64, string(x))
truth(x) = lowercase(string(x)) == "true"

function read_table(path)
    raw, header = readdlm(path, ','; header=true)
    columns = Dict(string(name) => i for (i, name) in enumerate(vec(header)))
    return [Dict(name => raw[i, column] for (name, column) in columns) for i in axes(raw, 1)]
end

# Retain fingerprints and all settings that affect the response or estimator.
function comparison_settings(document)
    comparable = deepcopy(document)
    for key in ("families", "sweeps", "output_directory", "dry_run", "make_plots")
        pop!(comparable["settings"], key, nothing)
    end
    return comparable
end

same_numbers(a, b) = all(isapprox.(a, b; rtol=1e-10, atol=1e-12))

# A relative error is undefined at zero growth/decay rate.
relative_error_percent(error, value) = iszero(value) ? NaN : 100error / abs(value)
successive_errors_percent(row) = ntuple(i -> relative_error_percent(row.previous[i], row.values[i]), 2)
relative_limits_percent(row) = ntuple(i -> relative_error_percent(row.limits[i], row.values[i]), 2)
reference_errors_percent(row, reference) =
    ntuple(i -> relative_error_percent(row.errors[i], reference.values[i]), 2)

"""Combine the two archived blade sweeps using full controls as the case identity.

`cases` contains each physical/numerical configuration once. `rows` contains its
membership in each sweep, so the common baseline can be used in both panels.
Saved errors and level numbers are discarded and recomputed on the merged grid.
"""
function load_results(directory=DEFAULT_DIRECTORY; case_names=CASE_NAMES)
    isempty(case_names) && error("At least one archived case is required")
    cases = NamedTuple[]
    memberships = Dict{Tuple{String,Float64},Int}()
    identities = Dict{Tuple,Int}()
    comparison = nothing
    settings = nothing
    input_count = 0
    for case_name in case_names
        archive = joinpath(directory, case_name)
        document = TOML.parsefile(joinpath(archive, "study.toml"))
        comparable = comparison_settings(document)
        if isnothing(comparison)
            comparison = comparable
            settings = document["settings"]
        else
            isequal(comparison, comparable) || error("Incompatible saved study settings or fingerprints: $case_name")
        end
        matrix = Dict((string(c["family"]), Int(number(c["level"]))) => c
            for c in read_table(joinpath(archive, "case_matrix.csv")))
        saved = filter(r -> string(r["family"]) in FAMILIES,
            read_table(joinpath(archive, "convergence.csv")))
        isempty(saved) && error("No blade convergence results in $archive")
        for r in saved
            input_count += 1
            family = string(r["family"])
            c = matrix[(family, Int(number(r["level"])))]
            label = string(r["label"])
            string(c["label"]) == label || error("Case label mismatch: $case_name/$label")
            string(r["status"]) == "completed" && truth(r["valid"]) ||
                error("Cannot compare an incomplete or invalid case: $case_name/$label")
            controls = Tuple(number(c[k]) for k in CONTROLS)
            all(x -> isfinite(x) && x > 0, controls) || error("Invalid controls: $case_name/$label")
            all(isinteger, controls[1:4]) || error("Noninteger panel counts: $case_name/$label")
            for (k, value) in zip(CONTROLS, controls)
                k == family && continue
                value == number(settings["nominal"][k]) ||
                    error("Non-swept control $k differs from nominal: $case_name/$label")
            end
            string(c["core_mode"]) == settings["core_mode"] || error("Core law mismatch: $case_name/$label")
            isapprox(number(c["dt_s"]), controls[7] / (6settings["physical"]["rpm"]); rtol=1e-10) ||
                error("Time step mismatch: $case_name/$label")
            values = Tuple(number(r[k]) for k in METRICS)
            limits = Tuple(number(r["limit_$k"]) for k in METRICS)
            all(isfinite, values) && all(x -> isfinite(x) && x > 0, limits) ||
                error("Invalid damping results or tolerances: $case_name/$label")
            all(>(0), values[3:4]) || error("Invalid fitted frequency: $case_name/$label")
            d = settings["damping"]
            expected_limits = (d["lambda_tolerance_per_s"], d["lambda_tolerance_per_s"],
                d["frequency_tolerance_hz"], d["frequency_tolerance_hz"])
            same_numbers(limits, expected_limits) || error("Tolerance mismatch: $case_name/$label")
            identity = (controls..., string(c["core_mode"]))
            source = "$case_name/$label"
            index = get(identities, identity, 0)
            if index == 0
                index = length(cases) + 1
                identities[identity] = index
                push!(cases, (; id=@sprintf("case_%03d", index), controls, values, limits,
                    sources=String[source], families=String[family]))
            else
                canonical = cases[index]
                same_numbers(values, canonical.values) && same_numbers(limits, canonical.limits) ||
                    error("Conflicting duplicate results: $source and $(first(canonical.sources))")
                push!(canonical.sources, source)
                family in canonical.families || push!(canonical.families, family)
            end
            member = (family, number(c[family]))
            get(memberships, member, index) == index ||
                error("Different configurations at the same $family panel count")
            memberships[member] = index
        end
    end
    rows = NamedTuple[]
    for family in FAMILIES
        levels = sort([x for (f, x) in keys(memberships) if f == family])
        length(levels) >= 3 || error("$family needs at least three distinct panel counts")
        reference = cases[memberships[(family, last(levels))]]
        for (i, x) in enumerate(levels)
            c = cases[memberships[(family, x)]]
            previous = i == 1 ? ntuple(_ -> NaN, 4) :
                abs.(c.values .- cases[memberships[(family, levels[i-1])]].values)
            errors = abs.(c.values .- reference.values)
            push!(rows, (; family, x, case_id=c.id, values=c.values, limits=c.limits,
                previous, errors, reference_panels=last(levels)))
        end
    end
    return (; rows, cases, settings, input_count)
end

function panel_options(family, group, nominal)
    radial = family == "prop_radial"
    ticks = [r.x for r in group]
    padding = 0.04(last(ticks) - first(ticks))
    fixed_panels = nominal[radial ? "prop_chord" : "prop_radial"]
    return (; title=radial ? "Blade radial panels (chord = $fixed_panels)" : "Blade chord panels (radial = $fixed_panels)",
        xlabel=radial ? "Radial panels per blade" : "Chordwise panels per blade",
        xticks=(ticks, [@sprintf("%g", x) for x in ticks]),
        xlims=(first(ticks)-padding, last(ticks)+padding),
        titlefontsize=11, guidefontsize=10, tickfontsize=9,
        legendfontsize=9, gridalpha=0.16, framestyle=:box,
        left_margin=8Plots.mm, bottom_margin=7Plots.mm,
        top_margin=4Plots.mm, right_margin=4Plots.mm)
end

function damping_panels(rows, nominal)
    panels = Any[]
    lo, hi = extrema([r.values[i] for r in rows for i in 1:2])
    padding = max(0.08(hi-lo), 0.01max(abs(lo), abs(hi)))
    error_values = [max(value, eps(Float64)) for r in rows for value in successive_errors_percent(r) if isfinite(value)]
    tolerance_values = [value for r in rows for value in relative_limits_percent(r) if isfinite(value)]
    isempty(error_values) && error("No finite relative damping errors to plot")
    error_lo = min(minimum(error_values), minimum(tolerance_values)) / 1.5
    error_hi = max(maximum(error_values), maximum(tolerance_values)) * 1.5
    error_ticks = filter(x -> error_lo <= x <= error_hi,
        [m * 10.0^e for e in floor(Int, log10(error_lo)):ceil(Int, log10(error_hi)) for m in (1, 2, 5)])
    for family in FAMILIES
        group = filter(r -> r.family == family, rows)
        p = plot(; panel_options(family, group, nominal)..., ylabel="Damping rate, λ (s⁻¹)",
            ylims=(lo-padding, hi+padding), legend=:bottomright)
        for (i, name) in enumerate(("Pitch", "Yaw"))
            plot!(p, [r.x for r in group], [r.values[i] for r in group];
                label=name, color=COLORS[i], marker=MARKERS[i],
                markersize=i == 1 ? 5 : 3.5, markerstrokecolor=COLORS[i],
                linestyle=i == 1 ? :solid : :dash, linewidth=2)
        end
        push!(panels, p)
    end
    for family in FAMILIES
        group = filter(r -> r.family == family, rows)
        refined = group[2:end]
        p = plot(; panel_options(family, group, nominal)...,
            ylabel="Successive relative damping error (%)",
            yscale=:log10, ylims=(error_lo, error_hi), legend=:topright,
            yticks=(error_ticks, [@sprintf("%g", x) for x in error_ticks]))
        for (i, name) in enumerate(("Pitch", "Yaw"))
            plot!(p, [r.x for r in refined], [max(successive_errors_percent(r)[i], eps(Float64)) for r in refined];
                label=name, color=COLORS[i], marker=MARKERS[i],
                markersize=i == 1 ? 5 : 3.5, markerstrokecolor=COLORS[i],
                linestyle=i == 1 ? :solid : :dash, linewidth=2)
        end
        for (i, name) in enumerate(("Pitch", "Yaw"))
            plot!(p, [r.x for r in refined], [relative_limits_percent(r)[i] for r in refined];
                label="$name tolerance (%)", color=COLORS[i], linestyle=:dot,
                linewidth=1.3, alpha=0.7)
        end
        push!(panels, p)
    end
    return panels
end

"""Show successive and finest-reference damping differences on linear percent axes."""
function damping_difference_panels(rows, nominal)
    panels = Any[]
    for comparison in (:successive, :reference)
        curves = map(FAMILIES) do family
            group = sort(filter(r -> r.family == family, rows); by=r -> r.x)
            reference = last(group)
            points = comparison == :successive ? group[2:end] : group
            errors = [comparison == :successive ? successive_errors_percent(r) :
                reference_errors_percent(r, reference) for r in points]
            limits = [relative_limits_percent(comparison == :successive ? r : reference)
                for r in points]
            (; family, group, reference, points, errors, limits)
        end
        finite_values = [value for curve in curves for pair in (curve.errors..., curve.limits...)
            for value in pair if isfinite(value)]
        ymax = isempty(finite_values) ? 1.0 : max(1.12maximum(finite_values), 1e-3)
        for curve in curves
            (; family, group, reference, points, errors, limits) = curve
            options = panel_options(family, group, nominal)
            subtitle = comparison == :successive ? "Successive cases" :
                @sprintf("Finest reference: %g panels", reference.x)
            options = merge(options, (; title=options.title * "\n" * subtitle))
            ylabel = comparison == :successive ? "Successive damping difference (%)" :
                "Difference from finest damping (%)"
            p = plot(; options..., ylabel, ylims=(0.0, ymax), legend=:topright)
            for (i, name) in enumerate(("Pitch", "Yaw"))
                plot!(p, [r.x for r in points], [pair[i] for pair in errors];
                    label=name, color=COLORS[i], marker=MARKERS[i],
                    markersize=i == 1 ? 5 : 3.5, markerstrokecolor=COLORS[i],
                    linestyle=i == 1 ? :solid : :dash, linewidth=2)
                plot!(p, [r.x for r in points], [pair[i] for pair in limits];
                    label="$name tolerance (%)", color=COLORS[i], linestyle=:dot,
                    linewidth=1.3, alpha=0.7)
            end
            push!(panels, p)
        end
    end
    return panels
end

function write_tables(out, results)
    open(joinpath(out, "unique_cases.csv"), "w") do io
        writedlm(io, permutedims(["case_id", CONTROLS..., METRICS..., "families", "sources"]), ',')
        for c in results.cases
            writedlm(io, permutedims(Any[c.id, c.controls..., c.values...,
                join(c.families, "; "), join(c.sources, "; ")]), ',')
        end
    end
    open(joinpath(out, "damping_convergence.csv"), "w") do io
        headers = ["family", "panels", "case_id", "reference_panels", METRICS...,
            ("previous_change_$k" for k in METRICS)...,
            ("error_$k" for k in METRICS)...,
            "previous_relative_pitch_percent", "previous_relative_yaw_percent",
            "reference_relative_pitch_percent", "reference_relative_yaw_percent",
            "limit_pitch_lambda_per_s", "limit_yaw_lambda_per_s"]
        writedlm(io, permutedims(headers), ',')
        for r in results.rows
            group = filter(a -> a.family == r.family, results.rows)
            relative_previous = successive_errors_percent(r)
            reference = last(group)
            relative_reference = reference_errors_percent(r, reference)
            writedlm(io, permutedims(Any[r.family, r.x, r.case_id, r.reference_panels,
                r.values..., r.previous..., r.errors..., relative_previous...,
                relative_reference..., r.limits[1:2]...]), ',')
        end
    end
end

function write_notes(out, results, case_names)
    s = results.settings
    open(joinpath(out, "plot_notes.md"), "w") do io
        println(io, "# Combined Chang blade-mesh damping convergence\n")
        println(io, "Sources: $(join(case_names, " + ")). Only archived convergence.csv, case_matrix.csv and study.toml are read; no simulations run.")
        println(io, "Source fingerprints and response, structural, physical, integration, coupling and damping settings match.")
        println(io, "Operating point: U = $(s["physical"]["speed_mps"]) m/s; alpha = $(s["physical"]["alpha_deg"]) deg; RPM = $(s["physical"]["rpm"]).")
        println(io, "Fixed wing mesh: $(s["nominal"]["wing_span"]) span x $(s["nominal"]["wing_chord"]) chord; azimuth increment = $(s["nominal"]["azimuth_deg"]) deg; wake = $(s["nominal"]["wake_revolutions"]) rev.")
        println(io, "Core law = $(s["core_mode"]); core value = $(s["nominal"]["core"]). Wing/propeller interaction = $(s["physical"]["interaction"]).\n")
        println(io, "$(results.input_count) input blade-family rows become $(length(results.rows)) distinct sweep points and $(length(results.cases)) unique configurations.")
        println(io, "unique_cases.csv stores each configuration once with all source labels. Configurations shared by both curves use one case_id.")
        println(io, "Duplicate identity uses all seven numerical controls and the core law; conflicting duplicate results cause an error. Saved level labels are not identities.\n")
        println(io, "The plotted damping metric is the saved moving-block lambda (1/s); negative lambda means decay. All selected fits are completed and valid.")
        println(io, "Plotted successive error (%) = 100 * abs(lambda(N_i) - lambda(N_(i-1))) / abs(lambda(N_i)) after sorting the combined panel counts. The first point has no previous error.")
        println(io, "Dotted pitch/yaw tolerance curves (%) = 100 * saved absolute lambda tolerance / abs(current lambda). The saved acceptance criterion is unchanged.")
        println(io, "Zero errors are displayed at machine epsilon on the log axis. Relative errors at zero lambda are undefined and saved as NaN.")
        println(io, "The CSV also contains errors against each merged finest-panel reference, recalculated rather than copied from the archives.")
        println(io, "Successive relative error (%) = 100 * successive absolute error / abs(current lambda); reference relative error uses abs(finest lambda).\n")
        println(io, "damping_percentage_differences.png/.svg shows both comparisons on linear percentage axes: successive cases in the top row and each case versus its family's finest reference in the bottom row. Radial/chord sweeps occupy the left/right columns.")
        println(io, "The finest reference has exactly zero reference difference, displayed at zero on this linear axis. Dotted tolerance curves use the same denominator as the corresponding difference; no simulation acceptance criteria are changed.\n")
        println(io, "| Sweep | Combined panels | First count whose remaining successive pitch/yaw changes meet tolerance | First count whose remaining pitch/yaw reference errors meet tolerance |")
        println(io, "|---|---|---|---|")
        for family in FAMILIES
            group = filter(r -> r.family == family, results.rows)
            successive = findfirst(i -> i > 1 && all(r -> all(r.previous[1:2] .<= r.limits[1:2]), group[i:end]), eachindex(group))
            reference = findfirst(i -> all(r -> all(r.errors[1:2] .<= r.limits[1:2]), group[i:end]), eachindex(group))
            format_count(i) = isnothing(i) ? "none" : @sprintf("%g", group[i].x)
            println(io, "| $family | $(join([@sprintf("%g", r.x) for r in group], ", ")) | $(format_count(successive)) | $(format_count(reference)) |")
        end
        println(io, "\nThese are separate one-parameter mesh sensitivity sweeps at fixed time step, not a test with both blade directions refined together. The finest mesh is a numerical reference, not an exact solution.")
        println(io, "The root directory's separate azimuth convergence run is excluded. These figures do not establish time-step convergence or a flutter boundary.")
    end
end

function plot_saved_convergence(directory=DEFAULT_DIRECTORY; case_names=CASE_NAMES)
    directory = abspath(directory)
    results = load_results(directory; case_names)
    out = joinpath(directory, "plots")
    mkpath(out)
    panels = damping_panels(results.rows, results.settings["nominal"])
    p = results.settings["physical"]
    title = @sprintf("%s | U = %g m/s, alpha = %g deg", join(case_names, " + "), p["speed_mps"], p["alpha_deg"])
    figure = plot(panels...; layout=(2, 2), size=(1400, 900), plot_title=title,
        plot_titlefontsize=14, background_color=:white, dpi=180)
    paths = [joinpath(out, "damping_convergence" * extension) for extension in (".png", ".svg")]
    foreach(path -> savefig(figure, path), paths)
    difference_figure = plot(damping_difference_panels(results.rows, results.settings["nominal"])...;
        layout=(2, 2), size=(1400, 950),
        plot_title="Damping convergence: percentage differences | " * title,
        plot_titlefontsize=13, background_color=:white, dpi=180)
    difference_paths = [joinpath(out, "damping_percentage_differences" * extension)
        for extension in (".png", ".svg")]
    foreach(path -> savefig(difference_figure, path), difference_paths)
    append!(paths, difference_paths)
    write_tables(out, results)
    write_notes(out, results, case_names)
    println("Merged $(results.input_count) archived rows into $(length(results.cases)) unique configurations.")
    foreach(println, paths)
    return paths
end
end

if abspath(PROGRAM_FILE) == @__FILE__
    ChangConvergenceAeroelasticPlots.plot_saved_convergence(isempty(ARGS) ?
        ChangConvergenceAeroelasticPlots.DEFAULT_DIRECTORY : only(ARGS))
end
