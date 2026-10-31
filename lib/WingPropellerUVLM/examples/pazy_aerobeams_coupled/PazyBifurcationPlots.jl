module PazyBifurcationPlots

using DelimitedFiles
export plot_pazy_bifurcation

const REFERENCE_AXES = (
    xlims=(30.0, 60.0), xticks=30:5:60,
    ylims=((15.0, 35.0), (-5000.0, 0.0), (-15.0, 15.0), (-10.0, 20.0)),
    yticks=(15:5:35, -5000:1000:0, -15:5:15, -10:5:20))

const PANEL_RESPONSES = (:tip_out_of_plane_percent, :root_strain_ue,
    :tip_twist_deg, :tip_geometric_incidence_deg)
# A dedicated label type selects GR's built-in TeX for Latin subscripts too.
# Plain strings containing only u_3 would otherwise be drawn as literal text.
struct NativeTeXLabel <: AbstractString
    value::String
end
Base.String(s::NativeTeXLabel) = s.value
Base.string(s::NativeTeXLabel) = s
Base.ncodeunits(s::NativeTeXLabel) = ncodeunits(s.value)
Base.codeunit(::Type{NativeTeXLabel}) = UInt8
Base.codeunit(s::NativeTeXLabel, i::Integer) = codeunit(s.value, i)
Base.isvalid(s::NativeTeXLabel, i::Integer) = isvalid(s.value, i)
Base.iterate(s::NativeTeXLabel, state...) = iterate(s.value, state...)
Base.getindex(s::NativeTeXLabel, i::Int) = s.value[i]

const PANEL_LABELS = (NativeTeXLabel("Tip OOP displacement, u_{3} [% semispan]"),
    raw"Root longitudinal strain, \epsilon_{11} [\mu]", raw"Tip twist, \theta [deg]",
    raw"Tip geometric incidence, \alpha_{g} [deg]")
const PANEL_CAPTIONS = ("(a) Tip OOP displacement.",
    "(b) Root longitudinal strain.", "(c) Tip twist angle.",
    "(d) Tip geometric incidence.")

function read_table(path)
    data, header = readdlm(path, ','; header=true)
    return data, Dict(String(n) => i for (i,n) in enumerate(vec(header)))
end

function saved_incidence(summary_path)
    path = joinpath(dirname(summary_path), "settings.txt")
    isfile(path) || throw(ArgumentError("Supply angle_of_attack_deg when settings.txt is unavailable"))
    found = match(r"\bangle_of_attack\s*=\s*([+\-0-9.eE]+)", read(path, String))
    isnothing(found) && throw(ArgumentError("Root incidence missing in settings.txt; supply angle_of_attack_deg"))
    return rad2deg(parse(Float64, found.captures[1]))
end

function response_group(data, columns, name, angle)
    selected = data[:,columns["response"]] .== String(name)
    offset = 0.0
    # Legacy CSVs can reconstruct geometric incidence exactly from tip twist.
    if !any(selected) && name == :tip_geometric_incidence_deg
        selected = data[:,columns["response"]] .== "tip_twist_deg"
        offset = angle
    end
    group = data[selected, :]
    if size(group, 1) > 0
        group = group[sortperm(Float64.(group[:,columns["airspeed_mps"]])), :]
    end
    return group, offset
end

curve_color(angle) = isapprox(angle, 3; atol=1e-8) ? "#0072BD" :
    isapprox(angle, 7; atol=1e-8) ? "#EDB120" : "#D95319"

"""Paper-style four-panel Pazy figure, using saved CSV data only.

Accept one summary path, or a vector of paths for overlays. Root incidences are
read from settings.txt unless angle_of_attack_deg is supplied (scalar/vector).
Every local minimum/maximum is read from bifurcation_local_extrema.csv. Legacy
summaries without that file use explicitly labelled final-window envelopes.
The fourth panel is geometric incidence (root incidence + elastic twist).
Defaults match the supplied figure: x=30:5:60 m/s and panel y ranges
(15,35), (-5000,0), (-15,15), (-10,20), with tick steps 5,1000,5,5.
Use xlims=nothing and four nothing ylims entries for automatic limits;
xticks=:auto and four :auto yticks entries restore automatic ticks.
The Hopf legend key is shown by default; show_hopf_key=false hides it if unused.
Hopf diamonds require independently identified points supplied as
`hopf_points=((angle_of_attack_deg=5.0, airspeed_mps=U,
responses=(tip_out_of_plane_percent=..., root_strain_ue=...,
tip_twist_deg=..., tip_geometric_incidence_deg=...)),)`.
No bifurcation location is inferred from an amplitude or a convergence flag.
Appearance controls are accepted directly here and in the runner's
paper_plot_options: linewidth (steady curves), markersize, markerstrokewidth,
fontfamily, guidefontsize, tickfontsize, legendfontsize, captionfontsize,
captionfontfamily, figure_size, panel_size, dpi, and gridalpha.
Latin and Greek subscripts use GR's built-in TeX, without external LaTeX.
"""
function plot_pazy_bifurcation(summary_paths;
    angle_of_attack_deg=nothing, output_directory=nothing,
    xlims=REFERENCE_AXES.xlims, xticks=REFERENCE_AXES.xticks,
    ylims=REFERENCE_AXES.ylims, yticks=REFERENCE_AXES.yticks,
    hopf_points=(), formats=("png", "pdf"), filename="pazy_bifurcation",
    save_panels=true, show_hopf_key=true,
    fontfamily="Times", guidefontsize=12, tickfontsize=11, legendfontsize=10,
    captionfontsize=11, captionfontfamily="Times Bold",
    linewidth=0.55, markersize=1.35, markerstrokewidth=0.4,
    figure_size=(1100,860), panel_size=(550,420), dpi=300, gridalpha=0.18)
    paths = summary_paths isa AbstractString ? [abspath(summary_paths)] : abspath.(summary_paths)
    isempty(paths) && throw(ArgumentError("Supply at least one summary CSV"))
    angles = isnothing(angle_of_attack_deg) ? saved_incidence.(paths) :
        angle_of_attack_deg isa Real ? [Float64(angle_of_attack_deg)] : collect(Float64, angle_of_attack_deg)
    length(angles) == length(paths) && all(isfinite, angles) ||
        throw(ArgumentError("Supply one finite root incidence per summary CSV"))
    length(ylims) == 4 || throw(ArgumentError("ylims must contain four limits or nothing entries"))
    length(yticks) == 4 || throw(ArgumentError("yticks must contain four tick specifications"))
    all(f -> f in ("png", "pdf", "svg"), formats) || throw(ArgumentError("Use png, pdf or svg formats"))
    all(v -> v isa Integer && v > 0,
        (guidefontsize, tickfontsize, legendfontsize, captionfontsize)) ||
        throw(ArgumentError("Font sizes must be positive integer point sizes"))
    all(v -> v isa Real && isfinite(v) && v > 0, (markersize, dpi)) ||
        throw(ArgumentError("markersize and dpi must be finite and positive"))
    all(v -> v isa Real && isfinite(v) && v >= 0, (linewidth, markerstrokewidth)) ||
        throw(ArgumentError("Line and marker stroke widths must be finite and nonnegative"))
    all(s -> length(s) == 2 && all(v -> v isa Integer && v > 0, s), (figure_size, panel_size)) ||
        throw(ArgumentError("figure_size and panel_size must be positive integer (width,height) tuples"))
    isfinite(gridalpha) && 0 <= gridalpha <= 1 || throw(ArgumentError("gridalpha must lie in [0,1]"))
    style = (; fontfamily, guidefontsize, tickfontsize, legendfontsize, captionfontsize, captionfontfamily,
        linewidth, markersize, markerstrokewidth, figure_size, panel_size, dpi, gridalpha)
    directory = isnothing(output_directory) ? dirname(first(paths)) : abspath(output_directory)
    datasets = map(zip(paths, angles)) do (path, angle)
        data, columns = read_table(path)
        for name in PANEL_RESPONSES
            group, _ = response_group(data, columns, name, angle)
            size(group, 1) > 0 || throw(ArgumentError("Response $name missing in $path"))
        end
        peak_path = joinpath(dirname(path), "bifurcation_local_extrema.csv")
        peaks, peak_columns = isfile(peak_path) ? read_table(peak_path) : (nothing, nothing)
        isnothing(peaks) && @warn "Local-extrema CSV unavailable; plotting final-window envelope only" path
        (; data, columns, peaks, peak_columns, angle, color=curve_color(angle))
    end
    if !isdefined(@__MODULE__, :Plots)
        @eval import Plots
    end
    return Base.invokelatest(plot_with_native_labels, datasets; directory, xlims, xticks, ylims, yticks,
        hopf_points, formats, filename, save_panels, show_hopf_key, style)
end

const NATIVE_LABELS_READY = Ref(false)

function plot_with_native_labels(datasets; kwargs...)
    previous_backend = Plots.backend_name()
    try
        Plots.gr()
        if !NATIVE_LABELS_READY[]
            @eval begin
                Plots.gr_text(x, y, s::NativeTeXLabel) = Plots.GR.textext(x, y, s.value)
                Plots.gr_inqtext(x, y, s::NativeTeXLabel) = Plots.GR.inqtextext(x, y, s.value)
            end
            NATIVE_LABELS_READY[] = true
        end
        return Base.invokelatest(plot_loaded, datasets; kwargs...)
    finally
        Plots.backend(previous_backend)
    end
end

function reference_legend!(p, datasets; all_local, show_hopf_key, style)
    conventions = [(all_local ? "Minima" : "Window minima", :dtriangle),
        (all_local ? "Maxima" : "Window maxima", :utriangle), ("Steady", :none)]
    show_hopf_key && push!(conventions, ("Hopf point", :diamond))
    # Draw the two independent columns directly. Empty native legend entries
    # trigger GR bitmap warnings, and otherwise force the columns to share rows.
    xmin, xmax = Plots.xlims(p)
    ymin, ymax = Plots.ylims(p)
    dx, dy = xmax-xmin, ymax-ymin
    left, top = xmin+0.015dx, ymax-0.005dy
    width = (all_local ? 0.44 : 0.60)*dx*style.legendfontsize/10
    row_height = 0.055dy*style.legendfontsize/10
    rows = max(length(conventions), length(datasets))
    right, bottom = left+width, top-(rows+0.3)*row_height
    box = Plots.Shape([left,right,right,left], [bottom,bottom,top,top])
    Plots.plot!(p, box; fillcolor=:white, fillalpha=1.0, linecolor=:black,
        linewidth=0.5, label=false)
    for i in 1:rows
        y = top-(i-0.4)*row_height
        if i <= length(conventions)
            label, marker = conventions[i]
            if marker == :none
                Plots.plot!(p, left .+ [0.006dx,0.042dx], [y,y];
                    color=:black, linestyle=:dash, linewidth=0.6, label=false)
            else
                Plots.scatter!(p, [left+0.025dx], [y]; marker, markersize=1.6,
                    markercolor=:white, markerstrokecolor=:black,
                    markerstrokewidth=0.4, label=false)
            end
            Plots.annotate!(p, left+0.055dx, y,
                Plots.text(label, style.legendfontsize, style.fontfamily, :left))
        end
        if i <= length(datasets)
            d = datasets[i]
            line_start = left+0.58width
            Plots.plot!(p, [line_start,line_start+0.045dx], [y,y];
                color=d.color, linewidth=1.6, label=false)
            angle_text = isinteger(d.angle) ? string(Int(d.angle)) : string(round(d.angle; digits=3))
            Plots.annotate!(p, line_start+0.05dx, y,
                Plots.text(raw"\alpha_{r} = " * angle_text * "°",
                    style.legendfontsize, style.fontfamily, :left))
        end
    end
end

function plot_loaded(datasets; directory, xlims, xticks, ylims, yticks,
    hopf_points, formats, filename, save_panels, show_hopf_key, style)
    mkpath(directory)
    all_speeds = vcat([Float64.(d.data[:,d.columns["airspeed_mps"]]) for d in datasets]...)
    limits = isnothing(xlims) ? extrema(all_speeds) : xlims
    if limits[1] == limits[2]
        limits = (limits[1]-0.5, limits[2]+0.5)
    end
    common = (; style.fontfamily, style.guidefontsize, style.tickfontsize,
        style.legendfontsize, framestyle=:axes, grid=true, style.gridalpha,
        foreground_color_grid=:gray, background_color=:white,
        foreground_color_legend=:gray, background_color_legend=:white,
        tick_direction=:in, left_margin=5Plots.mm, right_margin=3Plots.mm,
        bottom_margin=3Plots.mm, top_margin=2Plots.mm)
    panels, captions = [], []
    saved = String[]
    all_local = all(d -> !isnothing(d.peaks), datasets)
    for (i, name) in enumerate(PANEL_RESPONSES)
        p = Plots.plot(; common..., xlabel=NativeTeXLabel("Airspeed, U_{∞}  [m∕s]"),
            ylabel=PANEL_LABELS[i], xlims=limits, xticks, yticks=yticks[i],
            legend=false, widen=false)
        for d in datasets
            group, offset = response_group(d.data, d.columns, name, d.angle)
            speed = Float64.(group[:,d.columns["airspeed_mps"]])
            static = Float64.(group[:,d.columns["static_equilibrium"]]) .+ offset
            Plots.plot!(p, speed, static; color=d.color, linestyle=:dash,
                linewidth=style.linewidth, label=false)
            if isnothing(d.peaks)
                for (key, marker) in (("minimum", :dtriangle), ("maximum", :utriangle))
                    Plots.scatter!(p, speed, Float64.(group[:,d.columns[key]]) .+ offset;
                        marker, markersize=style.markersize, markercolor=:white,
                        markerstrokecolor=d.color, markerstrokewidth=style.markerstrokewidth, label=false)
                end
            else
                peaks, peak_offset = response_group(d.peaks, d.peak_columns, name, d.angle)
                for (kind, marker) in (("minimum", :dtriangle), ("maximum", :utriangle))
                    selected = peaks[:,d.peak_columns["kind"]] .== kind
                    Plots.scatter!(p, Float64.(peaks[selected,d.peak_columns["airspeed_mps"]]),
                        Float64.(peaks[selected,d.peak_columns["value"]]) .+ peak_offset;
                        marker, markersize=style.markersize, markercolor=:white,
                        markerstrokecolor=d.color, markerstrokewidth=style.markerstrokewidth, label=false)
                end
            end
            for point in hopf_points
                isapprox(point.angle_of_attack_deg, d.angle; atol=1e-8) || continue
                hasproperty(point.responses, name) || continue
                Plots.scatter!(p, [point.airspeed_mps], [getproperty(point.responses, name)];
                    marker=:diamond, markersize=3, markercolor=d.color,
                    markerstrokecolor=:black, markerstrokewidth=0.6, label=false)
            end
        end
        isnothing(ylims[i]) || Plots.ylims!(p, ylims[i])
        if i == 1
            reference_legend!(p, datasets; all_local, show_hopf_key=show_hopf_key || !isempty(hopf_points), style)
        end
        push!(panels, p)
        caption = Plots.plot(; xlims=(0,1), ylims=(0,1), axis=false, ticks=false,
            grid=false, framestyle=:none, margin=0Plots.mm, legend=false)
        Plots.annotate!(caption, 0.5, 0.5,
            Plots.text(PANEL_CAPTIONS[i], style.captionfontsize, style.captionfontfamily))
        push!(captions, caption)
        if save_panels
            Plots.plot!(p; size=style.panel_size, style.dpi)
            for format in formats
                path = joinpath(directory, "$(filename)_$(name).$(format)")
                Plots.savefig(p, path)
                push!(saved, path)
            end
        end
    end
    figure = Plots.plot(panels[1], panels[2], captions[1], captions[2],
        panels[3], panels[4], captions[3], captions[4];
        layout=Plots.grid(4,2; heights=[0.44,0.06,0.44,0.06]), size=style.figure_size, style.dpi)
    for format in formats
        path = joinpath(directory, "$(filename).$(format)")
        Plots.savefig(figure, path)
        push!(saved, path)
    end
    return saved
end

end # module
