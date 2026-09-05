# Aerodynamic-first UVLM convergence sweep for the rigid Chang wing--propeller
# system.  This driver compares periodic wing CL and propeller CT before the
# substantially more expensive aeroelastic damping sweep is attempted.

using Dates
using DelimitedFiles
using Printf
using Statistics
include(joinpath(@__DIR__, "chang_aerodynamic_study.jl"))
using .ChangAerodynamicStudy

const AERO_SWEEP_DIR = @__DIR__
const AERO_EXAMPLE_DIR = normpath(joinpath(AERO_SWEEP_DIR, "..", ".."))
const AERO_CASE_DRIVER = joinpath(AERO_SWEEP_DIR, "run_chang_coupled_aerodynamic_analysis.jl")
const AERO_PROJECT_DIR = normpath(joinpath(AERO_EXAMPLE_DIR, "..", ".."))

aero_sweep_float(name, default) = aero_float(name, default)
aero_sweep_int(name, default) = aero_int(name, default)
aero_sweep_bool(name, default) = aero_bool(name, default)

# Chord mode holds a physical radius fixed as span/radial grids change.
# Segment mode intentionally couples regularization to those mesh sizes.
const AERO_CORE_MODE = Symbol(lowercase(get(ENV, "CHANG_AERO_SWEEP_CORE_MODE", "chord")))
AERO_CORE_MODE in (:chord, :segment) || error("Core mode must be chord or segment")
const AERO_NOMINAL_WING_SPAN = aero_sweep_int("CHANG_AERO_WING_SPAN_PANELS", 30)
const AERO_NOMINAL_WING_CHORD = aero_sweep_int("CHANG_AERO_WING_CHORD_PANELS", 10)
const AERO_NOMINAL_PROP_RADIAL = aero_sweep_int("CHANG_AERO_PROP_RADIAL_PANELS", 10)
const AERO_NOMINAL_PROP_CHORD = aero_sweep_int("CHANG_AERO_PROP_CHORD_PANELS", 10)
const AERO_NOMINAL_WAKE_REVOLUTIONS = aero_sweep_float("CHANG_AERO_RETAINED_WAKE_REVOLUTIONS", 2.0)
const AERO_NOMINAL_CORE_FACTOR = AERO_CORE_MODE == :segment ?
    aero_sweep_float("CHANG_AERO_FCORE_SEGMENT_FACTOR", 0.25) : 0.0
const AERO_NOMINAL_CHORD_CORE_FACTOR = AERO_CORE_MODE == :chord ?
    aero_sweep_float("CHANG_AERO_FCORE_CHORD_FACTOR", 0.01) : 0.0
const AERO_NOMINAL_AZIMUTH_DEG = aero_sweep_float("CHANG_AERO_AZIMUTH_STEP_DEG", 5.0)

function sweep_levels(name, defaults, type = Float64; decreasing = false)
    values = parse.(type, strip.(split(get(ENV, name, join(defaults, ',')), ',')))
    length(values) >= 3 || error("$name requires at least three levels")
    all(value -> isfinite(value) && value > 0, values) || error("$name must be finite and positive")
    all(delta -> decreasing ? delta < 0 : delta > 0, diff(values)) ||
        error("$name must be strictly ordered from coarse to fine")
    return values
end

const AERO_WING_SPAN_LEVELS = sweep_levels("CHANG_AERO_SWEEP_WING_SPAN_LEVELS", [20,30,40], Int)
const AERO_WING_CHORD_LEVELS = sweep_levels("CHANG_AERO_SWEEP_WING_CHORD_LEVELS", [10,20,30], Int)
const AERO_PROP_RADIAL_LEVELS = sweep_levels("CHANG_AERO_SWEEP_PROP_RADIAL_LEVELS", [10,15,20], Int)
const AERO_PROP_CHORD_LEVELS = sweep_levels("CHANG_AERO_SWEEP_PROP_CHORD_LEVELS", [10,15,20], Int)
const AERO_WAKE_LEVELS = sweep_levels("CHANG_AERO_SWEEP_WAKE_LEVELS", [1,2,3])
const AERO_CORE_LEVELS = sweep_levels("CHANG_AERO_SWEEP_CORE_LEVELS",
    AERO_CORE_MODE == :chord ? [0.04,0.02,0.01,0.005] : [0.5,0.25,0.125,0.0625]; decreasing = true)
const AERO_AZIMUTH_LEVELS = sweep_levels("CHANG_AERO_SWEEP_AZIMUTH_LEVELS", [5,2.5,1]; decreasing = true)

function number_token(value)
    return replace(string(value), "." => "p", "-" => "m")
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
    chord_core_factor::Float64 = AERO_NOMINAL_CHORD_CORE_FACTOR
    azimuth_deg::Float64 = AERO_NOMINAL_AZIMUTH_DEG
end

case_key(case) = (
    case.wing_span,
    case.wing_chord,
    case.prop_radial,
    case.prop_chord,
    case.wake_revolutions,
    case.core_factor,
    case.chord_core_factor,
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
                label = "core_$(number_token(value))_$(AERO_CORE_MODE)",
                core_factor = AERO_CORE_MODE == :segment ? value : 0.0,
                chord_core_factor = AERO_CORE_MODE == :chord ? value : 0.0,
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

function case_options(case, case_directory; simulated_revolutions, averaged_revolutions)
    options = options_from_environment(child_environment(case, case_directory;
        simulated_revolutions, averaged_revolutions))
    validate_options(options)
    return options
end

function aerodynamic_output_matches_case(case_directory, case;
    simulated_revolutions = aero_sweep_int("CHANG_AERO_SWEEP_SIMULATED_REVOLUTIONS", 8),
    averaged_revolutions = aero_sweep_int("CHANG_AERO_SWEEP_AVERAGED_REVOLUTIONS", 2),
    warn_on_mismatch = true)
    options = case_options(case, case_directory; simulated_revolutions, averaged_revolutions)
    return output_matches_options(case_directory, options; warn_on_mismatch)
end

function read_periodic_result(case_directory, case, averaged_revolutions;
    simulated_revolutions = aero_sweep_int("CHANG_AERO_SWEEP_SIMULATED_REVOLUTIONS", 8))
    aerodynamic_output_matches_case(case_directory, case; simulated_revolutions,
        averaged_revolutions, warn_on_mismatch = false) ||
        error("Aerodynamic output metadata does not match the requested case")
    options = case_options(case, case_directory; simulated_revolutions, averaged_revolutions)
    history_path = joinpath(case_directory, "coupled_aerodynamic_history.csv")
    revolution_path = joinpath(case_directory, "coupled_aerodynamic_revolutions.csv")
    raw, header = readdlm(history_path, ',', header = true)
    headers = String.(vec(header))
    column(name) = Float64.(raw[:, something(findfirst(==(name), headers))])
    steps = validate_options(options)
    count = simulated_revolutions * steps
    size(raw, 1) == count || error("Incomplete coefficient history")
    dt = deg2rad(case.azimuth_deg) /
        (options.reference_rpm * 2pi / 60 * options.flow_speed_mps / options.reference_speed_mps)
    all(isapprox.(column("time_s"), collect(1:count) .* dt; rtol = 1e-10, atol = 1e-12)) ||
        error("Coefficient history has an inconsistent time grid")
    phases = mod.(collect(1:count) .* case.azimuth_deg, 360)
    phase_delta = mod.(column("azimuth_deg") .- phases .+ 180, 360) .- 180
    all(abs.(phase_delta) .<= 1e-8) || error("Coefficient history has an inconsistent azimuth grid")
    cl = periodic_metrics(column("wing_CL"), steps, averaged_revolutions)
    ct = periodic_metrics(column("propeller_CT"), steps, averaged_revolutions)
    cq = periodic_metrics(column("propeller_CQ"), steps, averaged_revolutions)
    return (; mean_cl = cl.mean, mean_ct = ct.mean, mean_cq = cq.mean,
        std_cl = cl.std, std_ct = ct.std, drift_cl = cl.drift, drift_ct = ct.drift,
        periodic_cl = cl.periodic_rms, periodic_ct = ct.periodic_rms, periodic_cq = cq.periodic_rms,
        phase_cl = cl.phase, phase_ct = ct.phase, phase_cq = cq.phase,
        history_path, revolution_path)
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
    environment["CHANG_AERO_FCORE_CHORD_FACTOR"] = string(case.chord_core_factor)
    environment["CHANG_AERO_PLOT_RESULTS"] = "false"
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
        aerodynamic_output_matches_case(case_directory, case; simulated_revolutions, averaged_revolutions)
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
        succeeded = try
            run_logged_child(command, environment, log_path; live_log)
        catch exception
            exception isa InterruptException && rethrow()
            @warn "Could not run aerodynamic child" exception
            false
        end
        elapsed_s = time() - start_time
        if !succeeded || !isfile(history_path) || !isfile(revolution_path)
            return failed_result(case, history_path, revolution_path, log_path, elapsed_s,
                succeeded ? "result files missing" : "child process failed")
        end
    else
        println("\n[aero sweep] Reusing $output_label")
    end

    periodic = try
        read_periodic_result(case_directory, case, averaged_revolutions; simulated_revolutions)
    catch exception
        exception isa InterruptException && rethrow()
        @warn "Invalid aerodynamic result" case_directory exception
        return failed_result(case, history_path, revolution_path, log_path, elapsed_s,
            sprint(showerror, exception))
    end
    return merge((; case, status = "completed", reason = "coefficient history verified",
        reused, elapsed_s, log_path), periodic)
end

function failed_result(case, history_path, revolution_path, log_path, elapsed_s, reason)
    return (; case, status = "failed", reason, reused = false, elapsed_s,
        mean_cl = NaN, mean_ct = NaN, mean_cq = NaN, std_cl = NaN, std_ct = NaN,
        drift_cl = NaN, drift_ct = NaN, periodic_cl = NaN, periodic_ct = NaN, periodic_cq = NaN,
        phase_cl = Float64[], phase_ct = Float64[], phase_cq = Float64[],
        history_path, revolution_path, log_path)
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

function annotate_results(results; cl_absolute_tolerance, ct_absolute_tolerance, relative_tolerance,
    cq_absolute_tolerance = aero_sweep_float("CHANG_AERO_SWEEP_CQ_ABS_TOL", 1e-6),
    periodic_cl_tolerance = aero_sweep_float("CHANG_AERO_SWEEP_PERIODIC_CL_TOL", 1e-4),
    periodic_ct_tolerance = aero_sweep_float("CHANG_AERO_SWEEP_PERIODIC_CT_TOL", 1e-6),
    periodic_cq_tolerance = aero_sweep_float("CHANG_AERO_SWEEP_PERIODIC_CQ_TOL", 1e-6))
    tolerances = (cl_absolute_tolerance, ct_absolute_tolerance, cq_absolute_tolerance,
        relative_tolerance, periodic_cl_tolerance, periodic_ct_tolerance, periodic_cq_tolerance)
    all(value -> isfinite(value) && value >= 0, tolerances) || error("Tolerances must be finite and nonnegative")
    periodic(result) = result.status == "completed" &&
        result.periodic_cl <= periodic_cl_tolerance &&
        result.periodic_ct <= periodic_ct_tolerance && result.periodic_cq <= periodic_cq_tolerance
    annotated = NamedTuple[]
    for family in unique(result.case.family for result in results)
        group = sort(filter(result -> result.case.family == family, results); by = r -> r.case.level)
        reference = last(group)
        expected = expected_final_level(family)
        complete = [r.case.level for r in group] == collect(1:expected) &&
            all(r -> r.status == "completed", group)
        # Keep diagnostic errors visible even when periodicity fails; the
        # acceptance path below still requires periodic candidate/reference loads.
        valid_reference = complete
        cl_limit = max(cl_absolute_tolerance, relative_tolerance * abs(reference.mean_cl))
        ct_limit = max(ct_absolute_tolerance, relative_tolerance * abs(reference.mean_ct))
        cq_limit = max(cq_absolute_tolerance, relative_tolerance * abs(reference.mean_cq))
        # A finest case cannot prove convergence by agreeing with itself.
        agrees(a, b) = periodic(a) && periodic(b) &&
            abs(a.mean_cl - b.mean_cl) <= cl_limit &&
            abs(a.mean_ct - b.mean_ct) <= ct_limit &&
            abs(a.mean_cq - b.mean_cq) <= cq_limit &&
            phase_rms_error(a.phase_cl, b.phase_cl) <= cl_limit &&
            phase_rms_error(a.phase_ct, b.phase_ct) <= ct_limit &&
            phase_rms_error(a.phase_cq, b.phase_cq) <= cq_limit
        reference_stable = valid_reference && agrees(group[end - 1], reference)
        for (index, result) in enumerate(group)
            errors = valid_reference && result.status == "completed" ?
                (abs(result.mean_cl - reference.mean_cl), abs(result.mean_ct - reference.mean_ct),
                 abs(result.mean_cq - reference.mean_cq),
                 phase_rms_error(result.phase_cl, reference.phase_cl),
                 phase_rms_error(result.phase_ct, reference.phase_ct),
                 phase_rms_error(result.phase_cq, reference.phase_cq)) : ntuple(_ -> NaN, 6)
            push!(annotated, merge(result, (;
                expected_family_levels = expected, reference_label = reference.case.label,
                cl_reference_error = errors[1], ct_reference_error = errors[2], cq_reference_error = errors[3],
                cl_phase_error = errors[4], ct_phase_error = errors[5], cq_phase_error = errors[6],
                cl_tolerance = cl_limit, ct_tolerance = ct_limit, cq_tolerance = cq_limit,
                periodic_cl_tolerance, periodic_ct_tolerance, periodic_cq_tolerance,
                periodic_converged = periodic(result), reference_stable,
                within_tolerance = reference_stable && all(r -> agrees(r, reference), group[index:end]),
            )))
        end
    end
    return annotated
end

function write_summary(path, results)
    headers = [
        "family", "level", "label", "status", "reason", "reused",
        "wing_span_panels", "wing_chord_panels", "prop_radial_panels", "prop_chord_panels",
        "wake_revolutions", "core_segment_factor", "core_chord_factor", "azimuth_step_deg",
        "mean_wing_CL", "mean_propeller_CT", "mean_propeller_CQ", "std_wing_CL", "std_propeller_CT",
        "last_revolution_CL_change", "last_revolution_CT_change",
        "periodic_CL_RMS", "periodic_CT_RMS", "periodic_CQ_RMS", "periodic_converged",
        "periodic_CL_tolerance", "periodic_CT_tolerance", "periodic_CQ_tolerance",
        "reference_label", "CL_reference_error", "CT_reference_error", "CQ_reference_error",
        "CL_phase_RMS_error", "CT_phase_RMS_error", "CQ_phase_RMS_error",
        "CL_tolerance", "CT_tolerance", "CQ_tolerance", "reference_stable",
        "within_tolerance", "expected_family_levels", "elapsed_s", "history_path", "log_path",
    ]
    open(path, "w") do stream
        writedlm(stream, permutedims(headers), ',')
        for result in results
            case = result.case
            values = [
                case.family, case.level, case.label, result.status, result.reason, result.reused,
                case.wing_span, case.wing_chord, case.prop_radial, case.prop_chord,
                case.wake_revolutions, case.core_factor, case.chord_core_factor, case.azimuth_deg,
                result.mean_cl, result.mean_ct, result.mean_cq, result.std_cl, result.std_ct,
                result.drift_cl, result.drift_ct, result.periodic_cl, result.periodic_ct,
                result.periodic_cq, result.periodic_converged,
                result.periodic_cl_tolerance, result.periodic_ct_tolerance, result.periodic_cq_tolerance,
                result.reference_label,
                result.cl_reference_error, result.ct_reference_error, result.cq_reference_error,
                result.cl_phase_error, result.ct_phase_error, result.cq_phase_error,
                result.cl_tolerance, result.ct_tolerance, result.cq_tolerance, result.reference_stable,
                result.within_tolerance, result.expected_family_levels, result.elapsed_s,
                result.history_path, result.log_path,
            ]
            writedlm(stream, permutedims(values), ',')
        end
    end
    return path
end

function write_report(path, results; simulated_revolutions, averaged_revolutions)
    open(path, "w") do stream
        println(stream, "# Chang rigid coupled aerodynamic convergence\n")
        println(stream, "Generated: $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
        println(stream, "Interaction enabled. Core mode: $AERO_CORE_MODE; radius = max(segment factor * ds, chord factor * c).")
        println(stream, "Here ds is the 3-D span/radial bound-edge length, and c is the full local chord.")
        println(stream, "Simulated/averaged revolutions: $simulated_revolutions/$averaged_revolutions.\n")
        println(stream, "Acceptance requires periodic CL/CT/CQ waveforms, agreement of the final two levels, and agreement of all finer levels with the reference.")
        println(stream, "Waveform RMS errors use both azimuth grids; standard deviation is a physical fluctuation measure, not a statistical confidence interval.\n")
        for family in unique(result.case.family for result in results)
            println(stream, "## $(replace(string(family), '_' => ' '))\n")
            println(stream, "| Level | Case | Status | CL | CT | CQ | Periodic | Accepted |")
            println(stream, "|---:|---|---|---:|---:|---:|:---:|:---:|")
            group = sort(filter(result -> result.case.family == family, results); by = r -> r.case.level)
            for result in group
                println(stream, "| $(result.case.level) | $(result.case.label) | $(result.status) | " *
                    "$(finite_or_blank(result.mean_cl)) | $(finite_or_blank(result.mean_ct)) | " *
                    "$(finite_or_blank(result.mean_cq)) | $(result.periodic_converged) | $(result.within_tolerance) |")
            end
            println(stream)
        end
        println(stream, "Finite core is a regularization-sensitivity study. A small radius is not proof of physical accuracy.")
        println(stream, "If periodicity fails, increase simulated revolutions before interpreting grid differences.")
        println(stream, "If the final pair disagrees, extend the levels. Confirm the selected settings together in a combined refinement study.")
    end
    return path
end

function plot_results(path, results)
    families = unique(result.case.family for result in results)
    panels = Any[]
    ratio(error, tolerance) = tolerance > 0 ? error / tolerance : iszero(error) ? 0.0 : NaN
    for family in families
        group = sort(filter(result -> result.case.family == family, results); by = r -> r.case.level)
        x = collect(eachindex(group))
        labels = [begin
            c = result.case
            family == :time_step ? "$(c.azimuth_deg) deg" :
            family == :wake_length ? "$(c.wake_revolutions) rev" :
            family == :finite_core ? (c.chord_core_factor > 0 ? "$(c.chord_core_factor) c" : "$(c.core_factor) ds") :
            string(getfield(c, family))
        end for result in group]
        cl = [ratio(max(r.cl_reference_error, r.cl_phase_error), r.cl_tolerance) for r in group]
        ct = [ratio(max(r.ct_reference_error, r.ct_phase_error), r.ct_tolerance) for r in group]
        cq = [ratio(max(r.cq_reference_error, r.cq_phase_error), r.cq_tolerance) for r in group]
        panel = plot(x, cl; marker = :circle, linewidth = 2, label = "CL",
            ylabel = "Error / tolerance", xticks = (x, labels), xrotation = 0,
            framestyle = :box, gridalpha = 0.25,
            title = replace(string(family), '_' => ' '))
        plot!(panel, x, ct; marker = :diamond, linewidth = 2, label = "CT")
        plot!(panel, x, cq; marker = :square, linewidth = 2, label = "CQ")
        hline!(panel, [1.0]; color = :black, linestyle = :dash, label = "tolerance")
        # A zero reference error is only a comparison, not a convergence result.
        accepted = findall(r -> r.within_tolerance, group)
        scatter!(panel, x[accepted], zeros(length(accepted)); marker = :star5,
            color = :green, markersize = 8, label = "accepted")
        push!(panels, panel)
    end
    columns = min(2, length(panels))
    rows = ceil(Int, length(panels) / columns)
    figure = plot(panels...; layout = (rows, columns),
        size = (760 * columns, 430 * rows),
        plot_title = "Aerodynamic convergence: mean and waveform errors",
        plot_titlefontsize = 12, bottom_margin = 5Plots.mm)
    savefig(figure, path)
    return path
end

function main()
    families = requested_families()
    cases = build_cases(families)
    simulated_revolutions = aero_sweep_int("CHANG_AERO_SWEEP_SIMULATED_REVOLUTIONS", 8)
    averaged_revolutions = aero_sweep_int("CHANG_AERO_SWEEP_AVERAGED_REVOLUTIONS", 2)
    reuse_existing = aero_sweep_bool("CHANG_AERO_SWEEP_REUSE_EXISTING", true)
    live_log = aero_sweep_bool("CHANG_AERO_SWEEP_LIVE_LOG", true)
    dry_run = aero_sweep_bool("CHANG_AERO_SWEEP_DRY_RUN", false)
    cl_absolute_tolerance = aero_sweep_float("CHANG_AERO_SWEEP_CL_ABS_TOL", 1e-4)
    # CT is O(1e-3) in this windmilling/interference case, so the CL floor of
    # 1e-4 would permit an unacceptably large relative thrust error.  Use a
    # smaller CT floor; the common relative tolerance still scales at larger CT.
    ct_absolute_tolerance = aero_sweep_float("CHANG_AERO_SWEEP_CT_ABS_TOL", 1e-6)
    relative_tolerance = aero_sweep_float("CHANG_AERO_SWEEP_REL_TOL", 0.01)
    # Validate every case before starting an expensive batch, including dry runs.
    for case in cases
        case_options(case, ""; simulated_revolutions, averaged_revolutions)
    end
    annotate_results(NamedTuple[]; cl_absolute_tolerance, ct_absolute_tolerance, relative_tolerance)
    make_plots = aero_sweep_bool("CHANG_AERO_SWEEP_PLOT_RESULTS", true)
    if make_plots && !dry_run
        @eval using Plots
    end
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
            "  %-12s L%d %-20s wing=%dx%d prop=%dx%d wake=%.3g rev core=max(%.3g ds, %.3g c) dpsi=%.3g deg\n",
            string(case.family), case.level, case.label, case.wing_span,
            case.wing_chord, case.prop_radial, case.prop_chord,
            case.wake_revolutions, case.core_factor, case.chord_core_factor, case.azimuth_deg,
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
            make_plots && Base.invokelatest(plot_results, plot_path, partial)
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
    make_plots && Base.invokelatest(plot_results, plot_path, results)
    println("\nSummary: $summary_path")
    println("Report:  $report_path")
    println("Plot:    $plot_path")
    return results
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
