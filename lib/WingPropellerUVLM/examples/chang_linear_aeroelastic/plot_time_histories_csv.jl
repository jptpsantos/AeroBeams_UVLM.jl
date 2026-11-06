# Plot saved CSV time histories without loading or running the aeroelastic model.
# From the repository root:
#   julia --project=. lib/WingPropellerUVLM/examples/chang_linear_aeroelastic/plot_time_histories_csv.jl
#   julia --project=. lib/WingPropellerUVLM/examples/chang_linear_aeroelastic/plot_time_histories_csv.jl "input.csv" "output.pdf"
# Edit the settings below, or pass CSV/output paths on the command line.
module CSVTimeHistoryPlots

# Use GR without opening a window. Set this to true to display after saving.
const DISPLAY_PLOT = false
DISPLAY_PLOT || get!(ENV, "GKSwstype", "100")

using DelimitedFiles
using Plots

# ------------------------- INPUT AND OUTPUT -------------------------
# The default input is relative to this script; command-line paths are relative
# to the current working directory. Absolute paths work in both cases.
# For an edited Windows path, use raw"C:\Users\Acer\Documents\history.csv" or /.
const INPUT_CSV = joinpath(@__DIR__, "output", "isolated_wing", "V90",
    "chang_isolated_wing_uvlm_history.csv")
const OUTPUT_PDF = nothing  # nothing: same directory and basename as input CSV.
const TIME_COLUMN = "time_s"
const TIME_SCALE = 1.0      # Plotted time = CSV time * TIME_SCALE + TIME_OFFSET.
const TIME_OFFSET = 0.0
const SPAN_M = 7.5         # Used only by the default displacement conversion.

# ------------------------- FIGURE SETTINGS --------------------------
const FIGURE_LAYOUT = nothing  # nothing: one panel per row; e.g. (1, 2).
const FIGURE_OPTIONS = (
    size = (1000, 1250),        # Width, height in Plots.jl layout units.
    dpi = 300,                  # Affects raster elements; PDF curves are vector.
    background_color = :white, # :transparent for a transparent background.
    foreground_color = :black,
    plot_title = "",            # Overall figure title.
    plot_titlefontsize = 20,
    link = :x,                  # :x, :y, :both, or :none.
)

# ---------------------- SHARED AXES AND STYLE -----------------------
# Add any Plots.jl keyword here. A panel's `options` override shared settings.
# Use :auto for automatic axis limits/ticks. No data are discarded by limits.
const AXES_OPTIONS = (
    fontfamily = "sans-serif", # e.g. "Times Roman", "Computer Modern".
    guidefontsize = 32,         # Axis label size, in points.
    tickfontsize = 28,
    titlefontsize = 30,         # Individual panel title size.
    legendfontsize = 14,
    xlabel = "Time (s)",
    xlims = :auto,              # :auto uses the full CSV time range; or (0, 2).
    ylims = :auto,              # Or e.g. (-0.3, 0.3), overridden per panel below.
    xticks = :auto,             # Or 0:0.5:5; or ([0, 1], ["start", "end"]).
    yticks = :auto,
    xscale = :identity,         # :identity, :log10, :ln, etc.
    yscale = :identity,
    xflip = false,
    yflip = false,
    xrotation = 0,
    yrotation = 0,
    framestyle = :axes,         # :axes, :box, :origin, :none, etc.
    tick_direction = :in,       # :in, :out, or :none.
    grid = false,
    gridalpha = 0.2,
    gridlinewidth = 0.6,
    gridstyle = :solid,
    minorgrid = false,
    minorticks = 0,
    legend = false,             # false, true, :topright, :bottomleft, etc.
    legend_background_color = :white,
    linewidth = 3.5,
    linestyle = :solid,         # :solid, :dash, :dot, :dashdot.
    linealpha = 1.0,
    color = :black,             # Named color or a hex string, e.g. "#2864A0".
    markershape = :none,        # :none, :circle, :square, :diamond, etc.
    markersize = 3,
    markerstrokewidth = 0.5,
    left_margin = 22Plots.mm,
    right_margin = 8Plots.mm,
    top_margin = 6Plots.mm,
    bottom_margin = 12Plots.mm,
)

# --------------------------- PANEL SETTINGS -------------------------
# Copy/remove entries to plot any numeric columns from another CSV.
# Each panel displays y = CSV[column] * scale + offset.
# Default Chang convention: negative displacement / span * 100.
# To plot displacement in metres, set scale=1.0 and ylabel="Tip displacement (m)".
# Twist values in this CSV are ALREADY degrees: no rad2deg conversion is needed.
# All panel-specific Plots.jl keywords belong in `options`, including overrides
# of shared fonts, margins, line styles, axis scales, ticks, limits, and legends.
const PANELS = [
    (
        column = "tip_displacement_m",
        scale = -100.0 / SPAN_M,
        offset = 0.0,
        options = (
            ylabel = "Tip displacement (%)",
            xlabel = "",        # Only label the bottom panel's time axis.
            title = "",
            label = "Tip displacement",
            ylims = (-5.0, 7.0),       # Example: (-5.0, 7.0).
            yticks = -5:2:7,      # Example: -4:2:6.
        ),
    ),
    (
        column = "tip_twist_deg",
        scale = 1.0,
        offset = 0.0,
        options = (
            ylabel = "Tip twist angle (deg)",
            title = "",
            label = "Tip twist",
            ylims = (-0.4, 0.3),       # Example: (-0.35, 0.35).
            yticks = -0.4:0.1:0.3,      # Example: -0.3:0.1:0.3.
        ),
    ),
]

"""Read a headered, comma-separated CSV, retaining mixed diagnostic columns."""
function read_history_csv(path)
    isfile(path) || throw(ArgumentError("CSV file not found: $path"))
    data, header = readdlm(path, ','; header = true)
    names = [strip(replace(string(name), '\ufeff' => "")) for name in vec(header)]
    length(unique(names)) == length(names) ||
        throw(ArgumentError("CSV contains duplicate column names: $path"))
    size(data, 1) >= 2 || throw(ArgumentError("CSV needs at least two data rows: $path"))
    size(data, 2) == length(names) ||
        throw(ArgumentError("CSV data width does not match its header: $path"))
    return data, Dict(name => i for (i, name) in enumerate(names))
end

function numeric_column(data, columns, name)
    haskey(columns, name) || throw(ArgumentError(
        "Column '$name' not found. Available columns: $(join(sort(collect(keys(columns))), ", "))"))
    return map(axes(data, 1)) do row
        value = tryparse(Float64, strip(string(data[row, columns[name]])))
        isnothing(value) && throw(ArgumentError(
            "Column '$name' has a nonnumeric value at CSV row $(row + 1)"))
        isfinite(value) || throw(ArgumentError(
            "Column '$name' has a nonfinite value at CSV row $(row + 1)"))
        value
    end
end

"""
    plot_time_histories(csv_path=INPUT_CSV; output_path=OUTPUT_PDF, ...)

Read selected numeric columns and save one vector PDF. Returns `(figure, path)`.
Optional keywords: `time_column`, `time_scale`, `time_offset`, `panels`,
`figure_options`, `axes_options`, `layout`, and `display_plot`. Paths are resolved
against the working directory. With no output path, use the CSV basename + .pdf.
Shared and per-panel plotting options are NamedTuples of Plots.jl attributes.
"""
function plot_time_histories(csv_path::AbstractString = INPUT_CSV;
    output_path = OUTPUT_PDF,
    time_column = TIME_COLUMN,
    time_scale = TIME_SCALE,
    time_offset = TIME_OFFSET,
    panels = PANELS,
    figure_options = FIGURE_OPTIONS,
    axes_options = AXES_OPTIONS,
    layout = FIGURE_LAYOUT,
    display_plot = DISPLAY_PLOT,
)
    isempty(panels) && throw(ArgumentError("Configure at least one panel"))
    source = abspath(csv_path)
    target = isnothing(output_path) ? splitext(source)[1] * ".pdf" : abspath(output_path)
    lowercase(splitext(target)[2]) == ".pdf" ||
        throw(ArgumentError("Output path must have a .pdf extension: $target"))
    (source == target || (isfile(source) && isfile(target) && samefile(source, target))) &&
        throw(ArgumentError("Input and output must be different files"))

    data, columns = read_history_csv(source)
    time = numeric_column(data, columns, time_column) .* time_scale .+ time_offset
    all(isfinite, time) && all(>(0), diff(time)) || throw(ArgumentError(
        "Time values must be finite and strictly increasing after scaling"))

    # Validate all requested data before creating a figure or writing output.
    values = map(panels) do panel
        y = numeric_column(data, columns, panel.column) .* panel.scale .+ panel.offset
        all(isfinite, y) || throw(ArgumentError(
            "Scaling column '$(panel.column)' produced nonfinite values"))
        y
    end

    gr()
    subplots = map(panels, values) do panel, y
        options = merge(axes_options, panel.options)
        # Automatic x limits mean exactly the recorded time range.
        if get(options, :xlims, :auto) === :auto
            options = merge(options, (; xlims = extrema(time)))
        end
        plot(time, y; options...)
    end
    final_layout = isnothing(layout) ? (length(panels), 1) : layout
    options = merge(figure_options, (; layout = final_layout))
    figure = plot(subplots...; options...)
    mkpath(dirname(target))
    savefig(figure, target)
    display_plot && display(figure)
    println("PDF written to $target ($(length(time)) samples, $(length(panels)) panels)")
    return (; figure, path = target)
end

function main(args = ARGS)
    if args in (["--help"], ["-h"])
        println("Usage: julia --project=. plot_time_histories_csv.jl [input.csv] [output.pdf]")
        println("No arguments: use editable settings in this script.")
        println("One argument: plot that CSV; OUTPUT_PDF=nothing saves beside it.")
        return nothing
    end
    length(args) <= 2 || throw(ArgumentError("Expected [input.csv] [output.pdf]; use --help"))
    input = isempty(args) ? INPUT_CSV : args[1]
    output = length(args) == 2 ? args[2] : OUTPUT_PDF
    return plot_time_histories(input; output_path = output)
end

end # module


    CSVTimeHistoryPlots.main()

