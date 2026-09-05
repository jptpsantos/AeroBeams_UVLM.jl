# Aerodynamic-first UVLM convergence sweep for the rigid Chang wing--propeller
# system.  This driver compares periodic wing CL and propeller CT before the
# substantially more expensive aeroelastic damping sweep is attempted.

using Dates
using DelimitedFiles
using Printf
using Statistics
using Plots

const AERO_SWEEP_DIR = @__DIR__
const AERO_EXAMPLE_DIR = normpath(joinpath(AERO_SWEEP_DIR, "..", ".."))
const AERO_CASE_DRIVER = joinpath(AERO_SWEEP_DIR, "run_chang_coupled_aerodynamic_analysis.jl")
const AERO_PROJECT_DIR = normpath(joinpath(AERO_EXAMPLE_DIR, "..", ".."))

const AERO_NOMINAL_WING_SPAN = 30
const AERO_NOMINAL_WING_CHORD = 10
const AERO_NOMINAL_PROP_RADIAL = 10
const AERO_NOMINAL_PROP_CHORD = 10
const AERO_NOMINAL_WAKE_REVOLUTIONS = 2.0
const AERO_NOMINAL_CORE_FACTOR = 0.25
const AERO_NOMINAL_AZIMUTH_DEG = 5.0

# Only one coordinate is refined in each family.  Unlike the aeroelastic
# sweep, the wing is rigid here, so spanwise wing refinement is purely
# aerodynamic and cannot change a structural natural frequency.
const AERO_WING_SPAN_LEVELS = [20, 30, 40]
const AERO_WING_CHORD_LEVELS = [10, 20, 30]
const AERO_PROP_RADIAL_LEVELS = [10, 15, 20]
const AERO_PROP_CHORD_LEVELS = [10, 15, 20]
const AERO_WAKE_LEVELS = [1.0, 2.0, 3.0]
# This is a regularization-sensitivity family. Decreasing the core is not, by
# itself, proof of increasing physical accuracy.
const AERO_CORE_LEVELS = [0.1, 0.05, 0.025, 0.01]
# Retain three levels for a genuine time-step trend.  The 5-degree case is the
# nominal setting; 10 degrees supplies a coarse point and 2.5 degrees is the
# refined reference.
const AERO_AZIMUTH_LEVELS = [5.0, 2.5, 1.0]

aero_sweep_float(name, default) = parse(Float64, get(ENV, name, string(default)))
aero_sweep_int(name, default) = parse(Int, get(ENV, name, string(default)))

function aero_sweep_bool(name, default)
    value = lowercase(strip(get(ENV, name, string(default))))
    value in ("1", "true", "yes", "on") && return true
    value in ("0", "false", "no", "off") && return false
    error("$name must be true/false, yes/no, on/off, or 1/0")
end

function number_token(value)
    return replace(@sprintf("%.3f", value), "." => "p", "-" => "m")
end

Base.@kwdef struct AerodynamicSweepCase
    family::Symbol
    level::Int
    label::String
    wing_span::Int = AERO_NOMINAL_WING_SPAN
    wing_chord::Int = AERO_NOMINAL_WING_CHORD
    prop_radial::Int = AERO_NOMINAL_PROP_RADIAL
    prop_chord::Int = AERO_NOMINAL_PROP_CHORD
    wake_revolutions::Float64 = AERO_NOMINAL_WAKE_REVOLUTIONS
    core_factor::Float64 = AERO_NOMINAL_CORE_FACTOR
    azimuth_deg::Float64 = AERO_NOMINAL_AZIMUTH_DEG
end

case_key(case) = (
    case.wing_span,
    case.wing_chord,
    case.prop_radial,
    case.prop_chord,
    case.wake_revolutions,
    case.core_factor,
    case.azimuth_deg,
)

function requested_families()
    raw = lowercase(strip(get(
        ENV,
        "CHANG_AERO_SWEEP_FAMILIES",
        "wing_span,wing_chord,prop_radial,prop_chord,wake_length,finite_core,time_step",
    )))
    families = unique(Symbol.(filter(!isempty, strip.(split(raw, ',')))))
    allowed = Set((
        :wing_span,
        :wing_chord,
        :prop_radial,
        :prop_chord,
        :wake_length,
        :finite_core,
        :time_step,
    ))
    invalid = filter(family -> !(family in allowed), families)
    isempty(invalid) || error("Unknown aerodynamic sweep families: $(join(invalid, ", "))")
    isempty(families) && error("At least one aerodynamic sweep family is required")
    return families
end

function build_cases(families)
    cases = AerodynamicSweepCase[]
    if :wing_span in families
        for (level, value) in enumerate(AERO_WING_SPAN_LEVELS)
            push!(cases, AerodynamicSweepCase(
                family = :wing_span,
                level = level,
                label = "span_$(value)",
                wing_span = value,
            ))
        end
    end
    if :wing_chord in families
        for (level, value) in enumerate(AERO_WING_CHORD_LEVELS)
            push!(cases, AerodynamicSweepCase(
                family = :wing_chord,
                level = level,
                label = "chord_$(value)",
                wing_chord = value,
            ))
        end
    end
    if :prop_radial in families
        for (level, value) in enumerate(AERO_PROP_RADIAL_LEVELS)
            push!(cases, AerodynamicSweepCase(
                family = :prop_radial,
                level = level,
                label = "radial_$(value)",
                prop_radial = value,
            ))
        end
    end
    if :prop_chord in families
        for (level, value) in enumerate(AERO_PROP_CHORD_LEVELS)
            push!(cases, AerodynamicSweepCase(
                family = :prop_chord,
                level = level,
                label = "chord_$(value)",
                prop_chord = value,
            ))
        end
    end
    if :wake_length in families
        for (level, value) in enumerate(AERO_WAKE_LEVELS)
            push!(cases, AerodynamicSweepCase(
                family = :wake_length,
                level = level,
                label = "wake_$(number_token(value))_rev",
                wake_revolutions = value,
            ))
        end
    end
    if :finite_core in families
        for (level, value) in enumerate(AERO_CORE_LEVELS)
            push!(cases, AerodynamicSweepCase(
                family = :finite_core,
                level = level,
                label = "core_$(number_token(value))_ds",
                core_factor = value,
            ))
        end
    end
    if :time_step in families
        for (level, value) in enumerate(AERO_AZIMUTH_LEVELS)
            push!(cases, AerodynamicSweepCase(
                family = :time_step,
                level = level,
                label = "azimuth_$(number_token(value))_deg",
                azimuth_deg = value,
            ))
        end
    end
    return cases
end

function aerodynamic_output_matches_case(case_directory, case; warn_on_mismatch = true)
    summary_path = joinpath(case_directory, "coupled_aerodynamic_summary.txt")
    isfile(summary_path) || return false
    lines = readlines(summary_path)
    line_with(prefix) = begin
        index = findfirst(line -> startswith(line, prefix), lines)
        isnothing(index) ? nothing : lines[index]
    end
    wing_line = line_with("Wing mesh:")
    propeller_line = line_with("Propeller mesh:")
    azimuth_line = line_with("Azimuth step:")
    wake_line = line_with("Retained wake:")
    core_line = line_with("Corrected finite core:")
    interaction_line = line_with("Interaction enabled:")
    speed_line = line_with("Flow speed:")
    aoa_line = line_with("Angle of attack:")
    rpm_line = line_with("RPM:")
    any(isnothing, (
        wing_line,
        propeller_line,
        azimuth_line,
        wake_line,
        core_line,
        interaction_line,
        speed_line,
        aoa_line,
        rpm_line,
    )) &&
        return false

    wing_match = match(r"Wing mesh: (\d+) span x (\d+) chord panels", wing_line)
    propeller_match = match(
        r"Propeller mesh: (\d+) radial x (\d+) chord panels per blade",
        propeller_line,
    )
    azimuth_match = match(r"Azimuth step: ([^ ]+) deg", azimuth_line)
    wake_match = match(r"Retained wake: ([^ ]+) revolutions", wake_line)
    core_match = match(r"Corrected finite core: max\(([^ ]+) ds, ([^ ]+) c\)", core_line)
    interaction_match = match(r"Interaction enabled: (true|false)", interaction_line)
    speed_match = match(r"Flow speed: ([^ ]+) m/s", speed_line)
    aoa_match = match(r"Angle of attack: ([^ ]+) deg", aoa_line)
    rpm_match = match(r"RPM: ([^ ]+)", rpm_line)
    any(isnothing, (
        wing_match,
        propeller_match,
        azimuth_match,
        wake_match,
        core_match,
        interaction_match,
        speed_match,
        aoa_match,
        rpm_match,
    )) &&
        return false

    expected_speed = aero_sweep_float("CHANG_AERO_SPEED_MPS", 80.0)
    expected_aoa = aero_sweep_float("CHANG_AERO_AOA_DEG", 0.0)
    reference_rpm = aero_sweep_float("CHANG_AERO_REFERENCE_RPM", 1217.6962)
    reference_speed = aero_sweep_float("CHANG_AERO_REFERENCE_SPEED_MPS", 65.0)
    expected_rpm = reference_rpm * expected_speed / reference_speed

    matches =
        parse(Int, wing_match.captures[1]) == case.wing_span &&
        parse(Int, wing_match.captures[2]) == case.wing_chord &&
        parse(Int, propeller_match.captures[1]) == case.prop_radial &&
        parse(Int, propeller_match.captures[2]) == case.prop_chord &&
        isapprox(parse(Float64, wake_match.captures[1]), case.wake_revolutions) &&
        isapprox(parse(Float64, core_match.captures[1]), case.core_factor) &&
        isapprox(parse(Float64, core_match.captures[2]), 0.0; atol = eps(Float64)) &&
        isapprox(parse(Float64, azimuth_match.captures[1]), case.azimuth_deg) &&
        parse(Bool, interaction_match.captures[1]) &&
        isapprox(parse(Float64, speed_match.captures[1]), expected_speed) &&
        isapprox(parse(Float64, aoa_match.captures[1]), expected_aoa) &&
        isapprox(parse(Float64, rpm_match.captures[1]), expected_rpm; rtol = 1e-8)
    if !matches && warn_on_mismatch
        @warn "Existing aerodynamic output does not match its requested case and will not be reused" case_directory case
    end
    return matches
end

function read_periodic_result(case_directory, case, averaged_revolutions)
    aerodynamic_output_matches_case(case_directory, case; warn_on_mismatch = false) ||
        error("Aerodynamic output metadata does not match case $(case.family)/$(case.label)")
    history_path = joinpath(case_directory, "coupled_aerodynamic_history.csv")
    revolution_path = joinpath(case_directory, "coupled_aerodynamic_revolutions.csv")
    raw_history, history_header = readdlm(history_path, ',', header = true)
    raw_revolutions, revolution_header = readdlm(revolution_path, ',', header = true)
    history_headers = String.(vec(history_header))
    revolution_headers = String.(vec(revolution_header))
    history_column(name) = Float64.(raw_history[:, something(findfirst(==(name), history_headers))])
    revolution_column(name) = Float64.(raw_revolutions[:, something(findfirst(==(name), revolution_headers))])

    cl = history_column("wing_CL")
    ct = history_column("propeller_CT")
    steps_per_revolution = round(Int, 360 / case.azimuth_deg)
    average_count = averaged_revolutions * steps_per_revolution
    length(cl) >= average_count || error("History is shorter than the requested average")
    average_indices = (length(cl) - average_count + 1):length(cl)
    revolution_cl = revolution_column("mean_wing_CL")
    revolution_ct = revolution_column("mean_propeller_CT")
    return (;
        mean_cl = mean(cl[average_indices]),
        mean_ct = mean(ct[average_indices]),
        std_cl = std(cl[average_indices]),
        std_ct = std(ct[average_indices]),
        drift_cl = revolution_cl[end] - revolution_cl[end - 1],
        drift_ct = revolution_ct[end] - revolution_ct[end - 1],
        history_path,
        revolution_path,
    )
end

function child_environment(case, case_directory; simulated_revolutions, averaged_revolutions)
    environment = copy(ENV)
    environment["CHANG_AERO_OUTPUT_DIR"] = case_directory
    environment["CHANG_AERO_WING_SPAN_PANELS"] = string(case.wing_span)
    environment["CHANG_AERO_WING_CHORD_PANELS"] = string(case.wing_chord)
    environment["CHANG_AERO_PROP_RADIAL_PANELS"] = string(case.prop_radial)
    environment["CHANG_AERO_PROP_CHORD_PANELS"] = string(case.prop_chord)
    environment["CHANG_AERO_RETAINED_WAKE_REVOLUTIONS"] = string(case.wake_revolutions)
    environment["CHANG_AERO_FCORE_SEGMENT_FACTOR"] = string(case.core_factor)
    environment["CHANG_AERO_FCORE_CHORD_FACTOR"] = "0.0"
    environment["CHANG_AERO_AZIMUTH_STEP_DEG"] = string(case.azimuth_deg)
    environment["CHANG_AERO_SIMULATED_REVOLUTIONS"] = string(simulated_revolutions)
    environment["CHANG_AERO_AVERAGED_REVOLUTIONS"] = string(averaged_revolutions)
    environment["CHANG_AERO_INTERACTION"] = "true"
    return environment
end

"""Run one child while teeing its combined stdout/stderr to its log and terminal."""
function run_logged_child(command, environment, log_path; live_log)
    if !live_log
        return open(log_path, "w") do stream
            process = run(
                pipeline(setenv(command, environment), stdout = stream, stderr = stream);
                wait = false,
            )
            wait(process)
            success(process)
        end
    end

    return open(log_path, "w") do stream
        process = run(
            pipeline(setenv(command, environment), stdout = stream, stderr = stream);
            wait = false,
        )
        byte_offset = 0
        last_heartbeat = time()
        function show_new_log_text()
            flush(stream)
            isfile(log_path) || return
            open(log_path, "r") do reader
                seek(reader, byte_offset)
                new_text = read(reader, String)
                byte_offset = position(reader)
                if !isempty(new_text)
                    print(new_text)
                    flush(stdout)
                end
            end
        end
        while Base.process_running(process)
            sleep(2.0)
            show_new_log_text()
            if time() - last_heartbeat >= 60.0
                @printf("[aero sweep] child still running (%.1f min)\n", (time() - last_heartbeat) / 60)
                flush(stdout)
                last_heartbeat = time()
            end
        end
        wait(process)
        show_new_log_text()
        return success(process)
    end
end

function run_case(
    case,
    output_directory;
    simulated_revolutions,
    averaged_revolutions,
    reuse_existing,
    live_log,
)
    output_label = "$(case.family)_L$(case.level)_$(case.label)"
    case_directory = joinpath(output_directory, output_label)
    history_path = joinpath(case_directory, "coupled_aerodynamic_history.csv")
    revolution_path = joinpath(case_directory, "coupled_aerodynamic_revolutions.csv")
    log_path = joinpath(output_directory, output_label * ".log")
    reusable_files = isfile(history_path) && isfile(revolution_path)
    reused = reuse_existing && reusable_files &&
        aerodynamic_output_matches_case(case_directory, case)
    elapsed_s = 0.0

    if !reused
        mkpath(case_directory)
        environment = child_environment(
            case,
            case_directory;
            simulated_revolutions,
            averaged_revolutions,
        )
        command = `$(Base.julia_cmd()) --startup-file=no --project=$(AERO_PROJECT_DIR) $(AERO_CASE_DRIVER)`
        println("\n[aero sweep] $(case.family) L$(case.level) $(case.label) -> $log_path")
        start_time = time()
        succeeded = run_logged_child(command, environment, log_path; live_log)
        elapsed_s = time() - start_time
        if !succeeded || !isfile(history_path) || !isfile(revolution_path)
            return (;
                case,
                status = "failed",
                reason = succeeded ? "result files missing" : "child process failed",
                reused = false,
                elapsed_s,
                mean_cl = NaN,
                mean_ct = NaN,
                std_cl = NaN,
                std_ct = NaN,
                drift_cl = NaN,
                drift_ct = NaN,
                history_path,
                revolution_path,
                log_path,
            )
        end
    else
        println("\n[aero sweep] Reusing $output_label")
    end

    periodic = read_periodic_result(case_directory, case, averaged_revolutions)
    return merge((;
        case,
        status = "completed",
        reason = "periodic coefficient history available",
        reused,
        elapsed_s,
        log_path,
    ), periodic)
end

finite_or_blank(value) = value isa Real && isfinite(value) ? @sprintf("%.12g", value) : ""

function expected_final_level(family)
    levels = Dict(
        :wing_span => AERO_WING_SPAN_LEVELS,
        :wing_chord => AERO_WING_CHORD_LEVELS,
        :prop_radial => AERO_PROP_RADIAL_LEVELS,
        :prop_chord => AERO_PROP_CHORD_LEVELS,
        :wake_length => AERO_WAKE_LEVELS,
        :finite_core => AERO_CORE_LEVELS,
        :time_step => AERO_AZIMUTH_LEVELS,
    )
    haskey(levels, family) || error("Unknown aerodynamic family: $family")
    return length(levels[family])
end

function annotate_results(results; cl_absolute_tolerance, ct_absolute_tolerance, relative_tolerance)
    annotated = NamedTuple[]
    for family in unique(result.case.family for result in results)
        group = sort(filter(result -> result.case.family == family, results); by = r -> r.case.level)
        reference = last(group)
        # During a long run, the last *completed* case is not necessarily the
        # configured finest reference. Keep convergence flags false until the
        # actual last level of this family is available.
        valid_reference = reference.status == "completed" &&
            reference.case.level == expected_final_level(family)
        cl_limit = valid_reference ? max(cl_absolute_tolerance, relative_tolerance * abs(reference.mean_cl)) : NaN
        ct_limit = valid_reference ? max(ct_absolute_tolerance, relative_tolerance * abs(reference.mean_ct)) : NaN
        for result in group
            cl_error = result.status == "completed" && valid_reference ?
                abs(result.mean_cl - reference.mean_cl) : NaN
            ct_error = result.status == "completed" && valid_reference ?
                abs(result.mean_ct - reference.mean_ct) : NaN
            push!(annotated, merge(result, (;
                reference_label = reference.case.label,
                cl_reference_error = cl_error,
                ct_reference_error = ct_error,
                cl_tolerance = cl_limit,
                ct_tolerance = ct_limit,
                within_tolerance = isfinite(cl_error) && isfinite(ct_error) &&
                    cl_error <= cl_limit && ct_error <= ct_limit,
            )))
        end
    end
    return annotated
end

function write_summary(path, results)
    headers = [
        "family", "level", "label", "status", "reason", "reused",
        "wing_span_panels", "wing_chord_panels", "prop_radial_panels",
        "prop_chord_panels", "wake_revolutions", "core_segment_factor",
        "azimuth_step_deg", "mean_wing_CL", "mean_propeller_CT",
        "std_wing_CL", "std_propeller_CT", "last_revolution_CL_change",
        "last_revolution_CT_change", "reference_label", "CL_reference_error",
        "CT_reference_error", "CL_tolerance", "CT_tolerance", "within_tolerance",
        "elapsed_s", "history_path", "log_path",
    ]
    open(path, "w") do stream
        println(stream, join(headers, ','))
        for result in results
            case = result.case
            values = [
                case.family, case.level, case.label, result.status, result.reason,
                result.reused, case.wing_span, case.wing_chord, case.prop_radial,
                case.prop_chord, case.wake_revolutions, case.core_factor,
                case.azimuth_deg, finite_or_blank(result.mean_cl),
                finite_or_blank(result.mean_ct), finite_or_blank(result.std_cl),
                finite_or_blank(result.std_ct), finite_or_blank(result.drift_cl),
                finite_or_blank(result.drift_ct), result.reference_label,
                finite_or_blank(result.cl_reference_error),
                finite_or_blank(result.ct_reference_error),
                finite_or_blank(result.cl_tolerance), finite_or_blank(result.ct_tolerance),
                result.within_tolerance, finite_or_blank(result.elapsed_s),
                result.history_path, result.log_path,
            ]
            println(stream, join(string.(values), ','))
        end
    end
    return path
end

function write_report(path, results; simulated_revolutions, averaged_revolutions)
    open(path, "w") do stream
        println(stream, "# Chang rigid coupled aerodynamic convergence")
        println(stream)
        println(stream, "Generated: $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
        println(stream, "Interaction: enabled (wing and propeller are in one UVLM solve).")
        println(stream, "Finite-core kernel: corrected regularized segment, epsilon = factor times 3-D segment length.")
        println(stream, "Simulated/averaged revolutions: $simulated_revolutions/$averaged_revolutions.")
        println(stream)
        for family in unique(result.case.family for result in results)
            println(stream, "## $(replace(string(family), '_' => ' '))")
            println(stream)
            println(stream, "| Level | Case | Wing CL | Propeller CT | Delta CL/rev | Delta CT/rev | Converged |")
            println(stream, "|---:|---|---:|---:|---:|---:|:---:|")
            group = sort(filter(result -> result.case.family == family, results); by = r -> r.case.level)
            for result in group
                println(
                    stream,
                    "| $(result.case.level) | $(result.case.label) | " *
                    "$(finite_or_blank(result.mean_cl)) | $(finite_or_blank(result.mean_ct)) | " *
                    "$(finite_or_blank(result.drift_cl)) | $(finite_or_blank(result.drift_ct)) | " *
                    "$(result.within_tolerance ? "yes" : "no") |",
                )
            end
            println(stream)
        end
        println(stream, "The finite-core table is a sensitivity test. Its smallest core is not automatically the physically preferred value.")
        println(stream, "Only settings with small revolution-to-revolution drift should be used to judge spatial convergence.")
    end
    return path
end

function plot_results(path, results)
    families = unique(result.case.family for result in results)
    panels = Any[]
    for family in families
        group = sort(filter(result -> result.case.family == family, results); by = r -> r.case.level)
        x = collect(eachindex(group))
        labels = [result.case.label for result in group]
        cl = [result.mean_cl for result in group]
        ct = [result.mean_ct for result in group]
        panel = plot(
            x,
            cl;
            marker = :circle,
            linewidth = 2,
            label = "wing CL",
            xticks = (x, labels),
            xrotation = 25,
            framestyle = :box,
            gridalpha = 0.25,
            title = replace(string(family), '_' => ' '),
        )
        plot!(panel, x, ct; marker = :diamond, linewidth = 2, label = "propeller CT")
        push!(panels, panel)
    end
    columns = min(2, length(panels))
    rows = ceil(Int, length(panels) / columns)
    figure = plot(
        panels...;
        layout = (rows, columns),
        size = (760 * columns, 430 * rows),
        plot_title = "Rigid coupled aerodynamic convergence",
    )
    savefig(figure, path)
    return path
end

function main()
    families = requested_families()
    cases = build_cases(families)
    simulated_revolutions = aero_sweep_int("CHANG_AERO_SWEEP_SIMULATED_REVOLUTIONS", 4)
    averaged_revolutions = aero_sweep_int("CHANG_AERO_SWEEP_AVERAGED_REVOLUTIONS", 1)
    reuse_existing = aero_sweep_bool("CHANG_AERO_SWEEP_REUSE_EXISTING", true)
    live_log = aero_sweep_bool("CHANG_AERO_SWEEP_LIVE_LOG", true)
    dry_run = aero_sweep_bool("CHANG_AERO_SWEEP_DRY_RUN", false)
    cl_absolute_tolerance = aero_sweep_float("CHANG_AERO_SWEEP_CL_ABS_TOL", 1e-4)
    # CT is O(1e-3) in this windmilling/interference case, so the CL floor of
    # 1e-4 would permit an unacceptably large relative thrust error.  Use a
    # smaller CT floor; the common relative tolerance still scales at larger CT.
    ct_absolute_tolerance = aero_sweep_float("CHANG_AERO_SWEEP_CT_ABS_TOL", 1e-6)
    relative_tolerance = aero_sweep_float("CHANG_AERO_SWEEP_REL_TOL", 0.01)
    run_stamp = Dates.format(now(), "yyyymmdd_HHMMSS")
    output_directory = abspath(get(
        ENV,
        "CHANG_AERO_SWEEP_OUTPUT_DIR",
        joinpath(AERO_EXAMPLE_DIR, "output", "coupled_aerodynamic_sweep_$run_stamp"),
    ))

    println("Chang aerodynamic-first UVLM convergence sweep")
    println("Generated $(length(cases)) entries ($(length(unique(case_key.(cases)))) unique cases)")
    for case in cases
        @printf(
            "  %-12s L%d %-20s wing=%dx%d prop=%dx%d wake=%.3g rev core=%.3g ds dpsi=%.3g deg\n",
            string(case.family), case.level, case.label, case.wing_span,
            case.wing_chord, case.prop_radial, case.prop_chord,
            case.wake_revolutions, case.core_factor, case.azimuth_deg,
        )
    end
    dry_run && return cases

    mkpath(output_directory)
    summary_path = joinpath(output_directory, "coupled_aerodynamic_convergence.csv")
    report_path = joinpath(output_directory, "coupled_aerodynamic_convergence.md")
    plot_path = joinpath(output_directory, "coupled_aerodynamic_convergence.png")
    cache = Dict{Tuple,NamedTuple}()
    results = NamedTuple[]
    for case in cases
        key = case_key(case)
        if haskey(cache, key)
            source = cache[key]
            result = merge(source, (; case, reused = true, elapsed_s = 0.0))
            println("\n[aero sweep] Reusing $(source.case.family)/$(source.case.label) for $(case.family)/$(case.label)")
        else
            result = run_case(
                case,
                output_directory;
                simulated_revolutions,
                averaged_revolutions,
                reuse_existing,
                live_log,
            )
            cache[key] = result
        end
        push!(results, result)
        @printf(
            "[aero sweep] %-9s CL=%+.8f CT=%+.8f drift=%+.2e/%+.2e elapsed=%.1f s\n",
            uppercase(result.status), result.mean_cl, result.mean_ct,
            result.drift_cl, result.drift_ct, result.elapsed_s,
        )
        partial = annotate_results(
            results;
            cl_absolute_tolerance,
            ct_absolute_tolerance,
            relative_tolerance,
        )
        write_summary(summary_path, partial)
        write_report(report_path, partial; simulated_revolutions, averaged_revolutions)
        try
            plot_results(plot_path, partial)
        catch exception
            @warn "Could not refresh the partial aerodynamic convergence plot" exception
        end
    end

    results = annotate_results(
        results;
        cl_absolute_tolerance,
        ct_absolute_tolerance,
        relative_tolerance,
    )
    write_summary(summary_path, results)
    write_report(report_path, results; simulated_revolutions, averaged_revolutions)
    plot_results(plot_path, results)
    println("\nSummary: $summary_path")
    println("Report:  $report_path")
    println("Plot:    $plot_path")
    return results
end

#if abspath(PROGRAM_FILE) == @__FILE__
    main()
#end
