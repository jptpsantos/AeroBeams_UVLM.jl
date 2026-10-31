# Standalone CSV postprocessor. Existing solvers, runners and saved files are read only.
module PazyInterpolatedBifurcation

using DelimitedFiles, Statistics, Plots
include(joinpath(@__DIR__, "PazyBifurcationPlots.jl"))
import .PazyBifurcationPlots as Paper
export plot_interpolated_bifurcation

# Match dev/sweptPazy/S0/PazyWingPitchRangeURangeLCOPlotGenerator.jl exactly.
# Keep the full incidence list so a 5-degree case keeps its color when plotted alone.
const AEROBEAMS_INCIDENCES = [1.0, 3.0, 5.0, 7.0, 10.0]
const AEROBEAMS_COLORS = cgrad(:rainbow, length(AEROBEAMS_INCIDENCES); categorical=true)
const DATA_ROOT = joinpath(@__DIR__, "output", "bifurcation")
const OUTPUT_ROOT = joinpath(@__DIR__, "output", "bifurcation_interpolated")

Plots.gr_text(x, y, s::Paper.NativeTeXLabel) = Plots.GR.textext(x, y, String(s))
Plots.gr_inqtext(x, y, s::Paper.NativeTeXLabel) = Plots.GR.inqtextext(x, y, String(s))

function incidence_color(angle)
    index = findfirst(a -> isapprox(a, angle; atol=1e-8), AEROBEAMS_INCIDENCES)
    # Other incidences use the same palette, continuously between its reference entries.
    isnothing(index) || return AEROBEAMS_COLORS[index]
    return AEROBEAMS_COLORS[clamp((angle-1)/9, 0.0, 1.0)]
end

"""Shape-preserving cubic Hermite interpolation (PCHIP), without extrapolation."""
function pchip(x, y, query)
    length(x) == length(y) > 0 || throw(ArgumentError("Empty or mismatched interpolation data"))
    all(isfinite, x) && all(isfinite, y) && all(isfinite, query) ||
        throw(ArgumentError("PCHIP requires finite values"))
    all(diff(x) .> 0) || throw(ArgumentError("Interpolation speeds must be strictly increasing"))
    all(q -> first(x) <= q <= last(x), query) || throw(ArgumentError("Extrapolation is disabled"))
    length(x) == 1 && return fill(Float64(first(y)), length(query))
    h = diff(x)
    delta = diff(y) ./ h
    slope = zeros(length(x))
    if length(x) == 2
        slope .= first(delta)
    else
        for k in 2:length(x)-1
            if delta[k-1]*delta[k] > 0
                w1, w2 = 2h[k]+h[k-1], h[k]+2h[k-1]
                slope[k] = (w1+w2)/(w1/delta[k-1]+w2/delta[k])
            end
        end
        function endpoint(h1, h2, d1, d2)
            d = ((2h1+h2)*d1-h1*d2)/(h1+h2)
            sign(d) != sign(d1) && return 0.0
            sign(d1) != sign(d2) && abs(d) > 3abs(d1) && return 3d1
            return d
        end
        slope[1] = endpoint(h[1], h[2], delta[1], delta[2])
        slope[end] = endpoint(h[end], h[end-1], delta[end], delta[end-1])
    end
    return map(query) do q
        k = min(searchsortedlast(x, q), length(x)-1)
        t = (q-x[k])/h[k]
        (2t^3-3t^2+1)*y[k] + (t^3-2t^2+t)*h[k]*slope[k] +
            (-2t^3+3t^2)*y[k+1] + (t^3-t^2)*h[k]*slope[k+1]
    end
end

# Split at NaNs and missing speeds. A user-supplied max_gap overrides the inferred
# 1.5 * median speed step (use Inf for an intentionally irregular grid).
function valid_segments(x, valid; max_gap=nothing)
    gap = isnothing(max_gap) ? (length(x) > 1 ? 1.5median(diff(x)) : Inf) : max_gap
    gap > 0 || throw(ArgumentError("max_gap must be positive"))
    segments = UnitRange{Int}[]
    start = 0
    for i in eachindex(x)
        if start > 0 && (!valid[i] || x[i]-x[i-1] > gap)
            push!(segments, start:i-1)
            start = 0
        end
        valid[i] && start == 0 && (start = i)
    end
    start > 0 && push!(segments, start:length(x))
    return segments
end

function interpolate_envelope(speed, lower, upper, steady; points_per_interval=32, max_gap=nothing)
    length(speed) == length(lower) == length(upper) == length(steady) ||
        throw(ArgumentError("Envelope vectors must have matching lengths"))
    all(isfinite, speed) && all(diff(speed) .> 0) ||
        throw(ArgumentError("Summary requires unique, finite, increasing speeds"))
    points_per_interval isa Integer && points_per_interval >= 2 ||
        throw(ArgumentError("points_per_interval must be an integer >= 2"))
    valid = isfinite.(lower) .& isfinite.(upper)
    all(lower[valid] .<= upper[valid]) || throw(ArgumentError("Minimum exceeds maximum"))
    curves = NamedTuple[]
    for indices in valid_segments(speed, valid; max_gap)
        x = speed[indices]
        query = length(x) == 1 ? copy(x) : vcat([
            collect(range(x[k], x[k+1]; length=points_per_interval+1))[1:end-1]
            for k in 1:length(x)-1]..., last(x))
        center = (lower[indices] .+ upper[indices]) ./ 2
        amplitude = (upper[indices] .- lower[indices]) ./ 2
        mid = pchip(x, center, query)
        amp = max.(0.0, pchip(x, amplitude, query))
        static = fill(NaN, length(query))
        for local_indices in valid_segments(x, isfinite.(steady[indices]); max_gap)
            sx = x[local_indices]
            select = findall(q -> first(sx) <= q <= last(sx), query)
            static[select] = pchip(sx, steady[indices][local_indices], query[select])
        end
        push!(curves, (; speed=query, minimum=mid.-amp, maximum=mid.+amp, static))
    end
    return curves
end

function load_case(path, angle; points_per_interval, max_gap)
    data, columns = Paper.read_table(path)
    all(haskey(columns, k) for k in ("airspeed_mps", "response", "minimum", "maximum", "static_equilibrium")) ||
        throw(ArgumentError("Required summary columns missing in $path"))
    peak_path = joinpath(dirname(path), "bifurcation_local_extrema.csv")
    peaks, peak_columns = isfile(peak_path) ? Paper.read_table(peak_path) : (nothing, nothing)
    panels = map(Paper.PANEL_RESPONSES) do response
        group, offset = Paper.response_group(data, columns, response, angle)
        size(group, 1) > 0 || throw(ArgumentError("Response $response missing in $path"))
        speed = Float64.(group[:,columns["airspeed_mps"]])
        lower = Float64.(group[:,columns["minimum"]]) .+ offset
        upper = Float64.(group[:,columns["maximum"]]) .+ offset
        steady = Float64.(group[:,columns["static_equilibrium"]]) .+ offset
        status = haskey(columns, "convergence_status") ? String.(group[:,columns["convergence_status"]]) :
            fill("not_checked", length(speed))
        verified = [s in ("stationary", "converged") for s in status]
        curves = interpolate_envelope(speed, lower, upper, steady; points_per_interval, max_gap)
        isempty(curves) && throw(ArgumentError("No finite envelope for $response in $path"))
        peak_speed, peak_value, peak_kind = Float64[], Float64[], String[]
        if !isnothing(peaks)
            pg, po = Paper.response_group(peaks, peak_columns, response, angle)
            peak_speed = Float64.(pg[:,peak_columns["airspeed_mps"]])
            peak_value = Float64.(pg[:,peak_columns["value"]]) .+ po
            peak_kind = String.(pg[:,peak_columns["kind"]])
        end
        (; response, speed, lower, upper, steady, verified, status, curves, peak_speed, peak_value, peak_kind)
    end
    return (; path, angle, color=incidence_color(angle), panels, has_local=!isnothing(peaks),
        run=basename(dirname(path)))
end

function saved_summaries()
    isdir(DATA_ROOT) || throw(ArgumentError("No saved bifurcation directory: $DATA_ROOT"))
    paths = [joinpath(d, "bifurcation_summary.csv") for d in readdir(DATA_ROOT; join=true) if isdir(d)]
    return sort(filter(isfile, paths))
end

"""
    plot_interpolated_bifurcation(summary_paths; linewidth=2.8, markersize=3.4, ...)

Read saved CSVs only. Export a four-panel figure, individual panels (PNG/PDF),
and interpolated_curves.csv into a NEW directory. Multiple paths overlay cases;
different simulation settings should normally be plotted separately.

Colors match AeroBeams' five-incidence categorical rainbow palette at 1/3/5/7/10
degrees. Solid lines interpolate the final-window envelope through every saved
speed; triangles retain raw local extrema (or window extrema for legacy CSVs).
PCHIP interpolates envelope center and nonnegative amplitude to prevent crossing.
This is an envelope interpolation, not continuation of periodic solutions: it
cannot resolve individual branches in a multi-valued cloud or identify a Hopf
point. Unverified envelope samples remain visible as black crosses.

No extrapolation or interpolation across missing/nonfinite samples is performed.
The fourth response is GEOMETRIC incidence, reconstructed from root incidence +
twist if absent; it is not effective aerodynamic angle of attack.

Examples from the repository root:
  julia --project=. lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/plot_pazy_bifurcation_interpolated.jl
  julia --project=. lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/plot_pazy_bifurcation_interpolated.jl --all
  # Or include this file and call with one/multiple summary paths and keywords.
"""
function plot_interpolated_bifurcation(summary_paths;
    angle_of_attack_deg=nothing, output_directory=nothing, linewidth=2.8,
    steady_linewidth=1.8, markersize=3.4, markerstrokewidth=0.65,
    points_per_interval=32, max_gap=nothing, xlims=nothing,
    ylims=(nothing,nothing,nothing,nothing), xticks=:auto,
    formats=("png", "pdf"), save_panels=true, show_local_extrema=true,
    figure_size=(1300,1040), panel_size=(650,460), dpi=300)
    paths = summary_paths isa AbstractString ? [abspath(summary_paths)] : abspath.(collect(summary_paths))
    isempty(paths) && throw(ArgumentError("Supply at least one summary CSV"))
    all(isfile, paths) || throw(ArgumentError("One or more summaries do not exist"))
    angles = isnothing(angle_of_attack_deg) ? Paper.saved_incidence.(paths) :
        angle_of_attack_deg isa Real ? [Float64(angle_of_attack_deg)] : Float64.(collect(angle_of_attack_deg))
    length(angles) == length(paths) && all(isfinite, angles) || throw(ArgumentError("Supply one finite incidence per CSV"))
    length(ylims) == 4 || throw(ArgumentError("ylims needs four limits or nothing entries"))
    !isempty(formats) && all(f -> f in ("png", "pdf", "svg"), formats) || throw(ArgumentError("Use png, pdf or svg formats"))
    all(v -> isfinite(v) && v > 0, (linewidth, steady_linewidth, markersize, dpi)) ||
        throw(ArgumentError("Linewidths, markersize and dpi must be positive"))
    isfinite(markerstrokewidth) && markerstrokewidth >= 0 || throw(ArgumentError("Invalid marker stroke width"))
    cases = [load_case(path, angle; points_per_interval, max_gap) for (path,angle) in zip(paths,angles)]
    directory = if isnothing(output_directory)
        mkpath(OUTPUT_ROOT)
        mktempdir(OUTPUT_ROOT; prefix="plots_", cleanup=false)
    else
        dest = abspath(output_directory)
        ispath(dest) && (!isdir(dest) || !isempty(readdir(dest))) &&
            throw(ArgumentError("Output directory must be new or empty: $dest"))
        mkpath(dest)
        dest
    end
    previous_backend = Plots.backend_name()
    try
        Plots.gr()
        return render_cases(cases; directory, linewidth, steady_linewidth, markersize, markerstrokewidth,
            xlims, ylims, xticks, formats, save_panels, show_local_extrema, figure_size, panel_size, dpi)
    finally
        Plots.backend(previous_backend)
    end
end

function render_cases(cases; directory, linewidth, steady_linewidth, markersize, markerstrokewidth,
    xlims, ylims, xticks, formats, save_panels, show_local_extrema, figure_size, panel_size, dpi)
    common = (; fontfamily="Times", guidefontsize=13, tickfontsize=11, legendfontsize=9,
        framestyle=:box, grid=true, gridalpha=0.18, foreground_color_grid=:gray,
        tick_direction=:in, background_color=:white, background_color_legend=:white,
        left_margin=6Plots.mm, right_margin=4Plots.mm, bottom_margin=4Plots.mm, top_margin=3Plots.mm)
    panels, captions, saved = Any[], Any[], String[]
    for (i, response) in enumerate(Paper.PANEL_RESPONSES)
        p = plot(; common..., xlabel=Paper.NativeTeXLabel("Airspeed, U_{∞} [m/s]"),
            ylabel=Paper.PANEL_LABELS[i], legend=i == 1 ? :topleft : false, xticks)
        for (j, case) in enumerate(cases)
            d = case.panels[i]
            # Symbol size is in points. Draw raw peaks first, then the curves.
            if show_local_extrema && case.has_local
                for (kind, marker) in (("minimum", :dtriangle), ("maximum", :utriangle))
                    chosen = d.peak_kind .== kind
                    scatter!(p, d.peak_speed[chosen], d.peak_value[chosen]; marker, markersize,
                        markercolor=:white, markerstrokecolor=case.color, markerstrokewidth,
                        markeralpha=0.6, label=false)
                end
            end
            for curve in d.curves
                plot!(p, curve.speed, curve.static; color=case.color, linestyle=:dash,
                    linewidth=steady_linewidth, label=false)
                plot!(p, curve.speed, curve.minimum; color=case.color, linewidth, label=false)
                plot!(p, curve.speed, curve.maximum; color=case.color, linewidth, label=false)
            end
            # Always expose the actual interpolation knots, including stationary cases.
            for (values, marker) in ((d.lower, :dtriangle), (d.upper, :utriangle))
                scatter!(p, d.speed, values; marker, markersize, markercolor=:white,
                    markerstrokecolor=case.color, markerstrokewidth, label=false)
                bad = .!d.verified
                any(bad) && scatter!(p, d.speed[bad], values[bad]; marker=:xcross,
                    markersize=markersize+0.7, markercolor=:black, markerstrokecolor=:black,
                    markerstrokewidth=1.1, label=false)
            end
            if i == 1
                angle_label = raw"\alpha_{r} = " * string(round(case.angle; digits=2)) * "°"
                duplicate_angle = count(c -> isapprox(c.angle,case.angle; atol=1e-8), cases) > 1
                case_label = duplicate_angle ? angle_label * " (" * case.run * ")" : angle_label
                plot!(p, [NaN], [NaN]; color=case.color, linewidth, label=Paper.NativeTeXLabel(case_label))
            end
        end
        if i == 1
            for (label, marker) in (("Minima", :dtriangle), ("Maxima", :utriangle))
                scatter!(p, [NaN], [NaN]; marker, markersize, markercolor=:white,
                    markerstrokecolor=:black, markerstrokewidth, label)
            end
            plot!(p, [NaN], [NaN]; color=:black, linewidth, label="Interpolated envelope")
            plot!(p, [NaN], [NaN]; color=:black, linestyle=:dash, linewidth=steady_linewidth,
                label="Static equilibrium")
            any(c -> any(.!c.panels[i].verified), cases) && scatter!(p, [NaN], [NaN];
                marker=:xcross, markersize=markersize+0.7, markercolor=:black,
                markerstrokewidth=1.1, label="Unverified envelope sample")
        end
        # Automatic axes show the complete data; optional limits can match the reference.
        isnothing(xlims) || Plots.xlims!(p, xlims)
        isnothing(ylims[i]) || Plots.ylims!(p, ylims[i])
        push!(panels, p)
        cap = plot(; xlims=(0,1), ylims=(0,1), axis=false, ticks=false,
            grid=false, framestyle=:none, margin=0Plots.mm, legend=false)
        annotate!(cap, 0.5, 0.5, text(Paper.PANEL_CAPTIONS[i], 12, "Times Bold"))
        push!(captions, cap)
        if save_panels
            plot!(p; size=panel_size, dpi)
            for format in formats
                path = joinpath(directory, "pazy_interpolated_$(response).$(format)")
                savefig(p, path)
                push!(saved, path)
            end
        end
    end
    footer = plot(; xlims=(0,1), ylims=(0,1), axis=false, ticks=false,
        grid=false, framestyle=:none, margin=0Plots.mm, legend=false)
    note = "Solid curves: interpolated final-window envelopes; triangles: recorded extrema."
    any(c -> !c.has_local, cases) && (note *= " Legacy CSV: window extrema only.")
    annotate!(footer, 0.5, 0.68, text(note, 10, "Times"))
    annotate!(footer, 0.5, 0.25, text("Fourth panel: geometric incidence (root incidence + twist).", 10, "Times"))
    title = length(cases) == 1 ? "Pazy bifurcation | saved run " * only(cases).run : "Pazy bifurcation | saved cases"
    figure = plot(panels[1], panels[2], captions[1], captions[2], panels[3], panels[4],
        captions[3], captions[4], footer;
        layout=@layout([grid(4,2, heights=[0.44,0.06,0.44,0.06]); footer{0.07h}]),
        size=figure_size, dpi, plot_title=title, plot_titlefontsize=14, fontfamily="Times")
    for format in formats
        path = joinpath(directory, "pazy_bifurcation_interpolated.$format")
        savefig(figure, path)
        push!(saved, path)
    end
    csv_path = joinpath(directory, "interpolated_curves.csv")
    open(csv_path, "w") do io
        println(io, "case,angle_of_attack_deg,response,airspeed_mps,minimum,maximum,static_equilibrium")
        for (j,case) in enumerate(cases), panel in case.panels, curve in panel.curves, k in eachindex(curve.speed)
            println(io, join((j, case.angle, panel.response, curve.speed[k], curve.minimum[k],
                curve.maximum[k], curve.static[k]), ','))
        end
    end
    push!(saved, csv_path)
    open(joinpath(directory, "plot_notes.txt"), "w") do io
        println(io, "PCHIP interpolation of final-window envelope center and half-range; no extrapolation.")
        println(io, "Interpolation does not identify stability, individual periodic branches, or Hopf points.")
        println(io, "Black crosses: samples without verified final-window convergence.")
        println(io, "Palette: cgrad(:rainbow, 5; categorical=true), incidences [1,3,5,7,10] deg.")
        println(io, "Linewidth=$linewidth; steady_linewidth=$steady_linewidth; markersize=$markersize pt.")
        for (j,case) in enumerate(cases)
            println(io, "\nCase $j: ", case.path, "\nRoot incidence [deg]: ", case.angle,
                "\nColor: ", case.color, "\nLocal extrema available: ", case.has_local)
            for panel in case.panels
                println(io, panel.response, ": ", length(panel.speed), " speeds, ",
                    count(!, panel.verified), " unverified; statuses: ", join(unique(panel.status), ", "))
            end
        end
    end
    return (; output_directory=directory, paths=saved, cases)
end

function main(args)
    if "--help" in args
        println("Usage: julia --project=. plot_pazy_bifurcation_interpolated.jl [--all | summary.csv ...]")
        println("No arguments: latest saved run. --all: separate figures for every saved run.")
        println("Multiple explicit CSV paths: overlay. All outputs go to new directories.")
        return
    end
    ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")
    if args == ["--all"]
        paths = saved_summaries()
        isempty(paths) && error("No saved summaries found")
        mkpath(OUTPUT_ROOT)
        root = mktempdir(OUTPUT_ROOT; prefix="all_cases_", cleanup=false)
        for path in paths
            result = plot_interpolated_bifurcation(path; output_directory=joinpath(root,basename(dirname(path))))
            println(result.output_directory)
        end
    else
        any(startswith("--"), args) && error("Unknown option; use --help")
        paths = isempty(args) ? saved_summaries() : args
        isempty(paths) && error("No saved summaries found")
        result = plot_interpolated_bifurcation(isempty(args) ? last(paths) : paths)
        foreach(println, result.paths)
    end
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    PazyInterpolatedBifurcation.main(ARGS)
end
