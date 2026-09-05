# Plotting and report generation for the legacy one-stage convergence study.

"""
Create publication-ready plots and a Markdown report from a completed Chang
UVLM convergence sweep.

Run from the repository root with

    julia --project=lib/WingPropellerUVLM `
      lib/WingPropellerUVLM/examples/chang_linear_aeroelastic/studies/convergence/plot_chang_uvlm_convergence_report.jl `
      lib/WingPropellerUVLM/examples/chang_linear_aeroelastic/output/uvlm_convergence_YYYYMMDD_HHMMSS

Set `CHANG_REPORT_PRE_CORE_CORRECTION=true` when documenting a sweep generated
before the corrected finite-core kernel was introduced. The report then marks
the finite-core family as historical rather than current-model evidence.
"""

ENV["GKSwstype"] = "100"

using Dates
using DelimitedFiles
using Printf
using Plots

const EXAMPLE_DIR = normpath(joinpath(@__DIR__, "..", ".."))

const FAMILY_ORDER = [
    "wing_mesh",
    "propeller_mesh",
    "wake_length",
    "finite_core",
    "time_step",
]
const FAMILY_TITLES = Dict(
    "wing_mesh" => "Wing mesh",
    "propeller_mesh" => "Propeller mesh",
    "wake_length" => "Wake length",
    "finite_core" => "Finite core",
    "time_step" => "Time step",
)
const PITCH_COLOR = RGB(0.10, 0.36, 0.67)
const YAW_COLOR = RGB(0.85, 0.33, 0.10)
const QUALITY_THRESHOLD = 0.5
const LAMBDA_TOLERANCE = 0.01

as_string(value) = strip(string(value))

function as_float(value)
    value isa Number && return Float64(value)
    text = as_string(value)
    isempty(text) && return NaN
    return something(tryparse(Float64, text), NaN)
end

function as_bool(value)
    text = lowercase(as_string(value))
    return text in ("true", "1", "yes", "on")
end

function environment_bool(name, default)
    value = lowercase(strip(get(ENV, name, string(default))))
    value in ("true", "1", "yes", "on") && return true
    value in ("false", "0", "no", "off") && return false
    error("$name must be true/false, yes/no, on/off, or 1/0")
end

function read_summary(path)
    data, header = readdlm(path, ',', Any, '\n'; header = true)
    names = Tuple(Symbol.(String.(vec(header))))
    return [NamedTuple{names}(Tuple(data[index, :])) for index in axes(data, 1)]
end

function family_rows(rows, family)
    selected = filter(row -> as_string(row.family) == family, rows)
    return sort(selected; by = row -> as_float(row.level))
end

function display_label(row)
    label = replace(as_string(row.label), '_' => ' ')
    label = replace(label, r"(?<=\d)p(?=\d)" => ".")
    label = replace(label, " ds" => " Δs", " deg" => "°")
    return label
end

function channel_quality(row, channel)
    key = channel == :pitch ? :pitch_fit_r_squared : :yaw_fit_r_squared
    return as_string(row.status) != "failed" && as_float(getproperty(row, key)) >= QUALITY_THRESHOLD
end

function mark_low_quality!(panel, x, values, group, channel, color; show_label = false)
    indices = [
        index for index in eachindex(group)
        if isfinite(values[index]) && !channel_quality(group[index], channel)
    ]
    isempty(indices) && return panel
    scatter!(
        panel,
        x[indices],
        values[indices];
        marker = :xcross,
        markersize = 7,
        markerstrokewidth = 2,
        color,
        label = show_label ? "fit R² < 0.5" : "",
    )
    return panel
end

function metric_figure(
    rows,
    pitch_key,
    yaw_key,
    ylabel,
    title,
    path;
    zero_line = false,
    quality_markers = true,
)
    panels = Any[]
    for (family_index, family) in enumerate(FAMILY_ORDER)
        group = family_rows(rows, family)
        x = collect(eachindex(group))
        pitch = [as_float(getproperty(row, pitch_key)) for row in group]
        yaw = [as_float(getproperty(row, yaw_key)) for row in group]
        labels = display_label.(group)

        panel = plot(
            x,
            pitch;
            marker = :circle,
            markersize = 5,
            linewidth = 2,
            color = PITCH_COLOR,
            label = family_index == 1 ? "pitch" : "",
            xticks = (x, labels),
            xrotation = 24,
            ylabel,
            title = FAMILY_TITLES[family],
            framestyle = :box,
            gridalpha = 0.22,
            legend = family_index == 1 ? :best : false,
            bottom_margin = 5Plots.mm,
        )
        plot!(
            panel,
            x,
            yaw;
            marker = :diamond,
            markersize = 5,
            linewidth = 2,
            color = YAW_COLOR,
            label = family_index == 1 ? "yaw" : "",
        )
        zero_line && hline!(panel, [0.0]; color = :black, linestyle = :dash, label = "")
        if quality_markers
            mark_low_quality!(
                panel,
                x,
                pitch,
                group,
                :pitch,
                PITCH_COLOR;
                show_label = family_index == 1,
            )
            mark_low_quality!(panel, x, yaw, group, :yaw, YAW_COLOR)
        end
        push!(panels, panel)
    end

    figure = plot(
        panels...;
        layout = (3, 2),
        size = (1500, 1250),
        plot_title = title,
        plot_titlefontsize = 16,
        left_margin = 4Plots.mm,
    )
    savefig(figure, path)
    return path
end

function quality_figure(rows, path)
    panels = Any[]
    for (family_index, family) in enumerate(FAMILY_ORDER)
        group = family_rows(rows, family)
        x = collect(eachindex(group))
        pitch = [as_float(row.pitch_fit_r_squared) for row in group]
        yaw = [as_float(row.yaw_fit_r_squared) for row in group]
        labels = display_label.(group)
        panel = plot(
            x,
            pitch;
            marker = :circle,
            linewidth = 2,
            color = PITCH_COLOR,
            label = family_index == 1 ? "pitch" : "",
            xticks = (x, labels),
            xrotation = 24,
            ylabel = "Fit R²",
            ylims = (0.0, 1.0),
            title = FAMILY_TITLES[family],
            framestyle = :box,
            gridalpha = 0.22,
            legend = family_index == 1 ? :best : false,
            bottom_margin = 5Plots.mm,
        )
        plot!(
            panel,
            x,
            yaw;
            marker = :diamond,
            linewidth = 2,
            color = YAW_COLOR,
            label = family_index == 1 ? "yaw" : "",
        )
        hline!(
            panel,
            [QUALITY_THRESHOLD];
            color = :black,
            linestyle = :dash,
            linewidth = 1.5,
            label = family_index == 1 ? "acceptance threshold" : "",
        )
        push!(panels, panel)
    end
    figure = plot(
        panels...;
        layout = (3, 2),
        size = (1500, 1250),
        plot_title = "Moving-block damping-fit quality",
        plot_titlefontsize = 16,
        left_margin = 4Plots.mm,
    )
    savefig(figure, path)
    return path
end

function reference_difference_figure(rows, path)
    panels = Any[]
    for (family_index, family) in enumerate(FAMILY_ORDER)
        group = family_rows(rows, family)
        x = collect(eachindex(group))
        labels = display_label.(group)
        pitch = [as_float(row.pitch_moving_block_lambda_per_s) for row in group]
        yaw = [as_float(row.yaw_moving_block_lambda_per_s) for row in group]
        pitch_reference = pitch[end]
        yaw_reference = yaw[end]
        pitch_difference = abs.(pitch .- pitch_reference)
        yaw_difference = abs.(yaw .- yaw_reference)

        panel = plot(
            x,
            pitch_difference;
            marker = :circle,
            linewidth = 2,
            color = PITCH_COLOR,
            label = family_index == 1 ? "pitch" : "",
            xticks = (x, labels),
            xrotation = 24,
            ylabel = "|Δλ| to finest (1/s)",
            title = FAMILY_TITLES[family],
            framestyle = :box,
            gridalpha = 0.22,
            legend = family_index == 1 ? :best : false,
            bottom_margin = 5Plots.mm,
        )
        plot!(
            panel,
            x,
            yaw_difference;
            marker = :diamond,
            linewidth = 2,
            color = YAW_COLOR,
            label = family_index == 1 ? "yaw" : "",
        )
        hline!(
            panel,
            [LAMBDA_TOLERANCE];
            color = :black,
            linestyle = :dash,
            linewidth = 1.5,
            label = family_index == 1 ? "0.01 1/s tolerance" : "",
        )
        push!(panels, panel)
    end
    figure = plot(
        panels...;
        layout = (3, 2),
        size = (1500, 1250),
        plot_title = "Raw decay-rate difference from each family's finest case",
        plot_titlefontsize = 16,
        left_margin = 4Plots.mm,
    )
    savefig(figure, path)
    return path
end

function runtime_figure(rows, path)
    executed = filter(row -> !as_bool(row.reused) && as_float(row.elapsed_s) > 0.0, rows)
    labels = ["$(FAMILY_TITLES[as_string(row.family)]): $(display_label(row))" for row in executed]
    hours = [as_float(row.elapsed_s) / 3600 for row in executed]
    colors = [as_string(row.status) == "failed" ? RGB(0.70, 0.15, 0.15) : RGB(0.25, 0.52, 0.38) for row in executed]
    positions = collect(eachindex(executed))
    figure = scatter(
        reverse(hours),
        positions;
        color = reverse(colors),
        markersize = 9,
        markerstrokecolor = :black,
        markerstrokewidth = 0.8,
        yticks = (positions, reverse(labels)),
        legend = false,
        xlabel = "Elapsed wall time (h, logarithmic scale)",
        xscale = :log10,
        framestyle = :box,
        gridalpha = 0.22,
        size = (1450, 900),
        left_margin = 42Plots.mm,
        title = "Executed-case computational cost",
    )
    savefig(figure, path)
    return path
end

format_number(value; digits = 4) = isfinite(value) ? @sprintf("%.*f", digits, value) : "--"

function maximum_quality(row)
    return min(as_float(row.pitch_fit_r_squared), as_float(row.yaw_fit_r_squared))
end

function last_step_change(group, key)
    length(group) < 2 && return NaN
    current = as_float(getproperty(group[end], key))
    previous = as_float(getproperty(group[end - 1], key))
    return isfinite(current) && isfinite(previous) ? abs(current - previous) : NaN
end

function write_report(path, rows, source_directory; pre_core_correction = false)
    completed_count = count(row -> as_string(row.status) == "completed", rows)
    indeterminate_count = count(row -> as_string(row.status) == "indeterminate", rows)
    failed_count = count(row -> as_string(row.status) == "failed", rows)

    open(path, "w") do stream
        println(stream, "# Chang UVLM convergence-analysis report")
        println(stream)
        println(stream, "Generated: $(Dates.format(now(), "yyyy-mm-dd HH:MM"))")
        println(stream)
        println(stream, "Source dataset: [`uvlm_convergence_summary.csv`](../uvlm_convergence_summary.csv)")
        println(stream)
        if pre_core_correction
            println(stream, "> [!WARNING]")
            println(stream, "> This dataset predates the corrected finite-core Biot--Savart kernel and the 3-D segment-length definition. The finite-core family is retained for provenance only and must not be used to select the current production core factor.")
            println(stream)
        end

        println(stream, "## Executive assessment")
        println(stream)
        println(stream, "The sweep contains **$(length(rows)) cases**: **$completed_count completed**, **$indeterminate_count indeterminate**, and **$failed_count failed** under the configured moving-block acceptance rules. Only estimates with both pitch and yaw fit quality of at least R² = 0.5 are suitable for the formal damping-convergence decision.")
        println(stream)
        println(stream, "The numerical frequencies are relatively stable for the wing mesh, propeller mesh, and wake-length families. The damping estimates are not sufficiently reliable to demonstrate convergence because nearly all moving-block fits have low R². The time-step family shows a persistent change in both frequency and decay rate, while the historical finite-core family is strongly sensitive to core size.")
        println(stream)
        println(stream, "## Plot conventions")
        println(stream)
        println(stream, "Blue circles denote pitch and orange diamonds denote yaw. Crosses superposed on a point indicate that its channel has R² < 0.5. Lines show raw numerical estimates and do not imply that an indeterminate family has converged.")
        println(stream)

        for (heading, filename, description) in (
            ("Moving-block decay rate", "01_decay_rate_convergence.png", "The directly fitted moving-block slope. Negative values indicate decay."),
            ("Derived damping ratio", "02_damping_ratio_convergence.png", "Damping ratio derived from the decay rate and fitted frequency."),
            ("Dominant response frequency", "03_frequency_convergence.png", "Dominant pitch and yaw frequencies extracted during moving-block processing."),
            ("Fit quality", "04_fit_quality_convergence.png", "The dashed line is the configured acceptance threshold, R² = 0.5."),
            ("Difference from the finest tested level", "05_lambda_difference_to_finest.png", "Raw decay-rate differences relative to the final configured case in each family. This is diagnostic when the reference fit itself is indeterminate."),
            ("Computational cost", "06_runtime_convergence.png", "Wall time for cases that were actually executed; reused nominal cases are omitted."),
        )
            println(stream, "## $heading")
            println(stream)
            println(stream, description)
            println(stream)
            println(stream, "![$heading]($filename)")
            println(stream)
        end

        println(stream, "## Numerical results")
        println(stream)
        for family in FAMILY_ORDER
            group = family_rows(rows, family)
            println(stream, "### $(FAMILY_TITLES[family])")
            println(stream)
            println(stream, "| Case | Status | Pitch λ (1/s) | Yaw λ (1/s) | Pitch ζ (%) | Yaw ζ (%) | Pitch f (Hz) | Yaw f (Hz) | Pitch R² | Yaw R² |")
            println(stream, "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|")
            for row in group
                println(
                    stream,
                    "| $(display_label(row)) | $(as_string(row.status)) | " *
                    "$(format_number(as_float(row.pitch_moving_block_lambda_per_s))) | " *
                    "$(format_number(as_float(row.yaw_moving_block_lambda_per_s))) | " *
                    "$(format_number(as_float(row.pitch_damping_percent))) | " *
                    "$(format_number(as_float(row.yaw_damping_percent))) | " *
                    "$(format_number(as_float(row.pitch_frequency_hz))) | " *
                    "$(format_number(as_float(row.yaw_frequency_hz))) | " *
                    "$(format_number(as_float(row.pitch_fit_r_squared); digits = 3)) | " *
                    "$(format_number(as_float(row.yaw_fit_r_squared); digits = 3)) |",
                )
            end
            println(stream)
            pitch_change = last_step_change(group, :pitch_moving_block_lambda_per_s)
            yaw_change = last_step_change(group, :yaw_moving_block_lambda_per_s)
            println(
                stream,
                "Final-step raw changes: pitch " *
                "$(format_number(pitch_change; digits = 5)) 1/s; yaw " *
                "$(format_number(yaw_change; digits = 5)) 1/s.",
            )
            println(stream)
        end

        println(stream, "## Interpretation by parameter")
        println(stream)
        println(stream, "- **Wing mesh:** the 30-to-40-panel decay-rate changes are below 0.01 1/s, and the frequencies are nearly unchanged. Nevertheless, the 20×3 case failed and the remaining fits have R² ≈ 0.11, so damping convergence is not established.")
        println(stream, "- **Propeller mesh:** the 10×10-to-12×12 decay-rate changes remain approximately 0.016 1/s, above the selected tolerance. The low fit quality prevents a formal conclusion.")
        println(stream, "- **Wake length:** the raw two-to-three-revolution change is below 0.005 1/s, suggesting a plateau, but every fit is below the quality threshold.")
        if pre_core_correction
            println(stream, "- **Finite core:** the large trend was obtained with the former defective kernel. It demonstrates sensitivity of the old calculation only; all finite-core cases must be rerun with the corrected kernel.")
        else
            println(stream, "- **Finite core:** the tested values show strong sensitivity and do not demonstrate convergence before the finest level.")
        end
        println(stream, "- **Time step:** both decay rate and frequency continue changing between 5° and 2.5° azimuth steps. A smaller step and a longer fit interval are required.")
        println(stream)
        println(stream, "## Recommended next analysis")
        println(stream)
        println(stream, "1. Rerun the sweep after the finite-core correction; do not reuse the histories in this directory.")
        println(stream, "2. Extend the response duration so that the moving-block envelope contains enough blocks and produces acceptable fit quality.")
        println(stream, "3. Repeat at least the 30×5 and 40×7 wing meshes, 10×10 and 12×12 propeller meshes, two and three wake revolutions, and azimuth steps of 5°, 2.5°, and a finer level.")
        println(stream, "4. Treat finite-core convergence separately from mesh convergence because the core radius scales with the local three-dimensional vortex-segment length.")
        println(stream)
        println(stream, "The source logs and histories remain in `$(source_directory)`.")
    end
    return path
end

function main()
    isempty(ARGS) && error(
        "Pass the convergence output directory as the first argument; " *
        "see the usage example at the top of this file.",
    )
    length(ARGS) <= 2 || error("Expected at most an input and report directory")
    source_directory = normpath(ARGS[1])
    summary_path = joinpath(source_directory, "uvlm_convergence_summary.csv")
    isfile(summary_path) || error("Convergence summary not found: $summary_path")

    report_directory = length(ARGS) >= 2 ? normpath(ARGS[2]) :
        joinpath(source_directory, "report_with_plots")
    mkpath(report_directory)

    rows = read_summary(summary_path)
    metric_figure(
        rows,
        :pitch_moving_block_lambda_per_s,
        :yaw_moving_block_lambda_per_s,
        "Moving-block λ (1/s)",
        "Chang UVLM moving-block decay-rate convergence",
        joinpath(report_directory, "01_decay_rate_convergence.png");
        zero_line = true,
    )
    metric_figure(
        rows,
        :pitch_damping_percent,
        :yaw_damping_percent,
        "Damping ratio ζ (%)",
        "Chang UVLM derived damping-ratio convergence",
        joinpath(report_directory, "02_damping_ratio_convergence.png"),
    )
    metric_figure(
        rows,
        :pitch_frequency_hz,
        :yaw_frequency_hz,
        "Dominant frequency (Hz)",
        "Chang UVLM dominant-frequency convergence",
        joinpath(report_directory, "03_frequency_convergence.png"),
    )
    quality_figure(rows, joinpath(report_directory, "04_fit_quality_convergence.png"))
    reference_difference_figure(
        rows,
        joinpath(report_directory, "05_lambda_difference_to_finest.png"),
    )
    runtime_figure(rows, joinpath(report_directory, "06_runtime_convergence.png"))

    pre_core_correction = environment_bool("CHANG_REPORT_PRE_CORE_CORRECTION", false)
    report_path = write_report(
        joinpath(report_directory, "chang_uvlm_convergence_report.md"),
        rows,
        source_directory;
        pre_core_correction,
    )
    println("Convergence report: $report_path")
    return report_path
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
