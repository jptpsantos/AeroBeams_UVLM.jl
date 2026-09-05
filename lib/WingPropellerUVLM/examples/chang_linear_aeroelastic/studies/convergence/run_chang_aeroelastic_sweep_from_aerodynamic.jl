# Gated aeroelastic damping validation driven by the preceding rigid
# aerodynamic CL/CT convergence sweep.
#
# The script reads `coupled_aerodynamic_convergence.csv`, verifies that the
# required aerodynamic families are complete, selects the first level that
# remains within the CL/CT tolerances, and builds a small aeroelastic case
# matrix.  The default `combined` stage compares the selected aerodynamic
# setup with a combined refined setup.  Use `attribution` only if that first
# comparison shows a material damping difference.

using Dates
using DelimitedFiles
using Printf

include(joinpath(@__DIR__, "chang_aeroelastic_convergence.jl"))
include(joinpath(@__DIR__, "chang_aerodynamic_study.jl"))
using .ChangAerodynamicStudy: metadata_options

const REQUIRED_AERODYNAMIC_FAMILIES = (
    :wing_span,
    :wing_chord,
    :prop_radial,
    :prop_chord,
    :wake_length,
    :finite_core,
    :time_step,
)

function latest_aerodynamic_summary()
    output_root = joinpath(CONVERGENCE_EXAMPLE_DIR, "output")
    isdir(output_root) || error("Aerodynamic output directory does not exist: $output_root")
    candidates = String[]
    for directory in readdir(output_root; join = true)
        startswith(basename(directory), "coupled_aerodynamic_sweep_") || continue
        summary = joinpath(directory, "coupled_aerodynamic_convergence.csv")
        isfile(summary) && push!(candidates, summary)
    end
    isempty(candidates) && error(
        "No coupled_aerodynamic_convergence.csv was found. Finish the aerodynamic sweep first.",
    )
    sort!(candidates; by = path -> stat(path).mtime)
    return last(candidates)
end

function read_aerodynamic_table(path)
    raw, header_matrix = readdlm(path, ',', header = true)
    headers = String.(vec(header_matrix))
    function index(name)
        column = findfirst(==(name), headers)
        isnothing(column) && error("Column '$name' is missing from $path")
        return column
    end
    rows = NamedTuple[]
    for row in axes(raw, 1)
        text(name) = string(raw[row, index(name)])
        number(name) = parse(Float64, text(name))
        integer(name) = round(Int, number(name))
        boolean(name) = lowercase(text(name)) in ("true", "1", "yes")
        push!(rows, (;
            family = Symbol(text("family")),
            level = integer("level"),
            label = text("label"),
            status = text("status"),
            wing_span = integer("wing_span_panels"),
            wing_chord = integer("wing_chord_panels"),
            prop_radial = integer("prop_radial_panels"),
            prop_chord = integer("prop_chord_panels"),
            wake_revolutions = number("wake_revolutions"),
            core_factor = number("core_segment_factor"),
            chord_core_factor = number("core_chord_factor"),
            expected_family_levels = integer("expected_family_levels"),
            periodic_converged = boolean("periodic_converged"),
            periodic_cl = number("periodic_CL_RMS"),
            periodic_ct = number("periodic_CT_RMS"),
            reference_stable = boolean("reference_stable"),
            azimuth_deg = number("azimuth_step_deg"),
            mean_cl = number("mean_wing_CL"),
            mean_ct = number("mean_propeller_CT"),
            drift_cl = number("last_revolution_CL_change"),
            drift_ct = number("last_revolution_CT_change"),
            within_tolerance = boolean("within_tolerance"),
            history_path = text("history_path"),
        ))
    end
    return rows
end

function verify_aerodynamic_table_provenance(rows)
    options_by_directory = Dict{String,Any}()
    operating_points = NamedTuple[]
    for row in rows
        row.status == "completed" || continue
        directory = dirname(row.history_path)
        options = get!(options_by_directory, directory) do
            metadata_options(directory)
        end
        # Check every row even when several families share one cached history.
        (options.wing_span_panels, options.wing_chord_panels,
            options.propeller_radial_panels, options.propeller_chord_panels,
            options.retained_wake_revolutions, options.finite_core_segment_factor,
            options.finite_core_chord_factor, options.azimuth_step_deg) ==
        (row.wing_span, row.wing_chord, row.prop_radial, row.prop_chord,
            row.wake_revolutions, row.core_factor, row.chord_core_factor, row.azimuth_deg) ||
            error("Aerodynamic metadata differs from row $(row.family)/$(row.label)")
        options.interaction_on || error("The gated aerodynamic study requires interaction")
        # These controls are fixed in the production Chang adapter.
        (options.air_density_kgpm3, options.sideslip_deg,
            options.collective_pitch_offset_deg, options.wake_relaxation) == (1.225, 0.0, 0.0, 0.1) ||
            error("Density, sideslip, collective pitch or shedding fraction differ from the production adapter; align them before the aeroelastic study")
        push!(operating_points, (;
            speed_mps = options.flow_speed_mps,
            angle_of_attack_deg = options.angle_of_attack_deg,
            rpm = options.reference_rpm * options.flow_speed_mps / options.reference_speed_mps))
    end
    isempty(operating_points) && error("The aerodynamic table has no completed result")
    all(==(first(operating_points)), operating_points) ||
        error("Aerodynamic cases do not share one operating point")
    return first(operating_points)
end

function completed_family(rows, family)
    group = sort(filter(row -> row.family == family, rows); by = row -> row.level)
    isempty(group) && error("Aerodynamic family '$family' is missing")
    expected = first(group).expected_family_levels
    expected >= 3 && all(row -> row.expected_family_levels == expected, group) &&
        [row.level for row in group] == collect(1:expected) ||
        error("Aerodynamic family '$family' is incomplete or has duplicate levels")
    all(row -> row.status == "completed", group) || error(
        "Aerodynamic family '$family' contains a failed or incomplete case",
    )
    return group
end

function selected_row(
    group;
    periodic_cl_tolerance,
    periodic_ct_tolerance,
    allow_finest_only,
)
    acceptable(row) = row.within_tolerance && row.periodic_converged && row.reference_stable &&
        row.periodic_cl <= periodic_cl_tolerance && row.periodic_ct <= periodic_ct_tolerance &&
        abs(row.drift_cl) <= periodic_cl_tolerance &&
        abs(row.drift_ct) <= periodic_ct_tolerance
    for index in eachindex(group)
        all(acceptable, group[index:end]) || continue
        if index == length(group) && !allow_finest_only
            error(
                "$(group[1].family) converged only to its own finest reference. " *
                "Add a finer aerodynamic level or set " *
                "CHANG_AE_ALLOW_FINEST_ONLY=true after reviewing the result.",
            )
        end
        return group[index]
    end
    error("No periodically converged aerodynamic level was found for $(group[1].family)")
end

function select_aerodynamic_configuration(
    rows;
    periodic_cl_tolerance,
    periodic_ct_tolerance,
    allow_finest_only,
)
    groups = Dict(
        family => completed_family(rows, family)
        for family in REQUIRED_AERODYNAMIC_FAMILIES
    )
    selected = Dict(
        family => selected_row(
            groups[family];
            periodic_cl_tolerance,
            periodic_ct_tolerance,
            allow_finest_only,
        )
        for family in REQUIRED_AERODYNAMIC_FAMILIES
    )

    # The production aeroelastic model is still collocated spanwise: changing
    # wing span panels also changes the structural beam. Structural modal
    # convergence established 30 elements, so require the rigid aerodynamic
    # 30-panel point to pass, then keep 30 in every aeroelastic case.
    wing_30_index = findfirst(row -> row.wing_span == 30, groups[:wing_span])
    isnothing(wing_30_index) && error("The wing-span family must contain 30 panels")
    wing_30 = groups[:wing_span][wing_30_index]
    wing_30.within_tolerance && wing_30.periodic_converged && wing_30.reference_stable || error(
        "Thirty span panels did not pass the rigid aerodynamic tolerance. " *
        "A work-conjugate noncollocated aero/structure mapping is required before refinement.",
    )
    max(abs(wing_30.drift_cl), wing_30.periodic_cl) <= periodic_cl_tolerance || error(
        "The 30-panel wing CL is not periodic enough for aeroelastic validation",
    )
    max(abs(wing_30.drift_ct), wing_30.periodic_ct) <= periodic_ct_tolerance || error(
        "The 30-panel wing CT is not periodic enough for aeroelastic validation",
    )

    configuration = (;
        wing_span = 30,
        wing_chord = selected[:wing_chord].wing_chord,
        prop_radial = selected[:prop_radial].prop_radial,
        prop_chord = selected[:prop_chord].prop_chord,
        wake_revolutions = selected[:wake_length].wake_revolutions,
        core_factor = selected[:finite_core].core_factor,
        chord_core_factor = selected[:finite_core].chord_core_factor,
        azimuth_deg = selected[:time_step].azimuth_deg,
    )
    finest = (;
        wing_span = 30,
        wing_chord = last(groups[:wing_chord]).wing_chord,
        prop_radial = last(groups[:prop_radial]).prop_radial,
        prop_chord = last(groups[:prop_chord]).prop_chord,
        wake_revolutions = last(groups[:wake_length]).wake_revolutions,
        # Core is a model sensitivity, not a refinement direction. Keep the
        # selected plateau value in the combined numerical-refinement case.
        core_factor = configuration.core_factor,
        chord_core_factor = configuration.chord_core_factor,
        azimuth_deg = last(groups[:time_step]).azimuth_deg,
    )
    return (; configuration, finest, groups, selected)
end

function convergence_case(family, level, label, configuration)
    return validate_case(UVLMConvergenceCase(
        family = family,
        level = level,
        label = label,
        wing_span_panels = configuration.wing_span,
        wing_chord_panels = configuration.wing_chord,
        prop_radial_panels = configuration.prop_radial,
        prop_chord_panels = configuration.prop_chord,
        wake_revolutions = configuration.wake_revolutions,
        fcore_segment_factor = configuration.core_factor,
        fcore_chord_factor = configuration.chord_core_factor,
        azimuth_step_deg = configuration.azimuth_deg,
    ))
end

function build_gated_aeroelastic_cases(selection, stage)
    baseline = selection.configuration
    finest = selection.finest
    cases = UVLMConvergenceCase[
        convergence_case(:combined_validation, 1, "aero_selected", baseline),
        convergence_case(:combined_validation, 2, "combined_refined", finest),
    ]
    stage == :combined && return cases
    stage == :attribution || error("CHANG_AE_VALIDATION_STAGE must be combined or attribution")

    controls = (
        (:wing_chord, :wing_chord),
        (:prop_radial, :prop_radial),
        (:prop_chord, :prop_chord),
        (:wake_length, :wake_revolutions),
        (:time_step, :azimuth_deg),
    )
    for (family, field) in controls
        selected_value = getfield(baseline, field)
        finer_value = getfield(finest, field)
        selected_value == finer_value && continue
        refined = merge(baseline, NamedTuple{(field,)}((finer_value,)))
        push!(cases, convergence_case(family, 1, "selected", baseline))
        push!(cases, convergence_case(family, 2, "refined", refined))
    end

    core_group = selection.groups[:finite_core]
    core_index = findfirst(row -> row.core_factor == baseline.core_factor && row.chord_core_factor == baseline.chord_core_factor, core_group)
    if !isnothing(core_index) && core_index < length(core_group)
        core_refined = merge(
            baseline,
            (core_factor = core_group[core_index + 1].core_factor,
             chord_core_factor = core_group[core_index + 1].chord_core_factor,),
        )
        push!(cases, convergence_case(:finite_core_sensitivity, 1, "selected", baseline))
        push!(cases, convergence_case(:finite_core_sensitivity, 2, "smaller_core", core_refined))
    end
    return cases
end

function write_configuration(path, source_path, selection, cases, stage, analysis)
    open(path, "w") do stream
        println(stream, "# Gated Chang aeroelastic convergence configuration")
        println(stream)
        println(stream, "Aerodynamic source: `$source_path`")
        println(stream, "Stage: `$stage`")
        println(stream, "Structural/spanwise wing mesh: 30 elements/panels (fixed)")
        println(stream, "Flow speed: $(analysis.speed_mps) m/s")
        println(stream, "Angle of attack: $(analysis.angle_of_attack_deg) deg")
        println(stream, "Propeller speed: $(analysis.rpm) RPM")
        println(stream, "Wing--propeller aerodynamic interaction: $(analysis.interaction_on)")
        println(stream, "Generalized-alpha rho_inf: $(analysis.ga_rho_inf)")
        println(stream, "Coupling relaxation: $(analysis.coupling_relaxation)")
        println(stream, "Simulation/fit interval: 0--$(analysis.end_time_s) s / " *
            "$(analysis.fit_start_s)--$(analysis.fit_end_s) s")
        println(stream, "Moving block: $(analysis.block_duration_s) s, " *
            "overlap=$(analysis.block_overlap), band=$(analysis.frequency_min_hz)--" *
            "$(analysis.frequency_max_hz) Hz, Hann=$(analysis.hann_window)")
        println(stream)
        println(stream, "Selected aerodynamic configuration:")
        println(stream, "```text")
        for name in propertynames(selection.configuration)
            println(stream, "$name = $(getfield(selection.configuration, name))")
        end
        println(stream, "```")
        println(stream)
        println(stream, "Aeroelastic cases:")
        println(stream)
        println(stream, "| Family | Level | Label | Wing | Propeller | Wake (rev) | Core/ds | Core/c | dpsi (deg) |")
        println(stream, "|---|---:|---|---:|---:|---:|---:|---:|---:|")
        for case in cases
            println(
                stream,
                "| $(case.family) | $(case.level) | $(case.label) | " *
                "$(case.wing_span_panels)x$(case.wing_chord_panels) | " *
                "$(case.prop_radial_panels)x$(case.prop_chord_panels) | " *
                "$(case.wake_revolutions) | $(case.fcore_segment_factor) | " *
                "$(case.fcore_chord_factor) | $(case.azimuth_step_deg) |",
            )
        end
    end
    return path
end

function main_aeroelastic_from_aerodynamic()
    # Do not evaluate latest_aerodynamic_summary() when the caller supplied an
    # explicit table: Julia evaluates ordinary function arguments eagerly.
    aerodynamic_summary = if haskey(ENV, "CHANG_AEROELASTIC_AERO_SUMMARY")
        abspath(ENV["CHANG_AEROELASTIC_AERO_SUMMARY"])
    else
        abspath(latest_aerodynamic_summary())
    end
    rows = read_aerodynamic_table(aerodynamic_summary)
    operating_point = verify_aerodynamic_table_provenance(rows)
    lowercase(get(ENV, "CHANG_CONVERGENCE_FORCE_MODEL", "imperial")) == "imperial" ||
        error("The rigid aerodynamic validation uses Imperial forces; the gated study must use the same force model")
    allow_finest_only = convergence_bool("CHANG_AE_ALLOW_FINEST_ONLY", false)
    selection = select_aerodynamic_configuration(
        rows;
        periodic_cl_tolerance = convergence_float("CHANG_AE_PERIODIC_CL_TOL", 1e-4),
        periodic_ct_tolerance = convergence_float("CHANG_AE_PERIODIC_CT_TOL", 1e-6),
        allow_finest_only,
    )
    stage = Symbol(lowercase(get(ENV, "CHANG_AE_VALIDATION_STAGE", "combined")))
    cases = build_gated_aeroelastic_cases(selection, stage)

    # Freeze all non-UVLM controls during the numerical comparison.  In
    # particular rho_inf=1 avoids adding high-frequency algorithmic damping
    # while the physical aeroelastic damping is the measured quantity.
    end_time_s = convergence_float("CHANG_AE_END_TIME_S", 6.0)
    fit_start_s = convergence_float("CHANG_AE_FIT_START_S", 0.8)
    fit_end_s = convergence_float("CHANG_AE_FIT_END_S", end_time_s)
    analysis = (;
        speed_mps = operating_point.speed_mps,
        angle_of_attack_deg = operating_point.angle_of_attack_deg,
        rpm = operating_point.rpm,
        interaction_on = convergence_bool("CHANG_AE_INTERACTION", true),
        ga_rho_inf = convergence_float("CHANG_AE_GA_RHO_INF", 1.0),
        coupling_relaxation = convergence_float("CHANG_AE_COUPLING_RELAXATION", 1.0),
        end_time_s,
        fit_start_s,
        fit_end_s,
        block_duration_s = convergence_float("CHANG_AE_MOVING_BLOCK_DURATION_S", 1.5),
        block_overlap = convergence_float("CHANG_AE_MOVING_BLOCK_OVERLAP", 0.9),
        frequency_min_hz = convergence_float("CHANG_AE_FREQUENCY_MIN_HZ", 3.0),
        frequency_max_hz = convergence_float("CHANG_AE_FREQUENCY_MAX_HZ", 6.0),
        hann_window = convergence_bool("CHANG_AE_HANN_WINDOW", true),
    )
    0 <= analysis.fit_start_s < analysis.fit_end_s <= analysis.end_time_s || error(
        "The damping-fit interval must lie inside the simulation interval",
    )
    0 <= analysis.ga_rho_inf <= 1 || error("CHANG_AE_GA_RHO_INF must be in [0,1]")
    0 < analysis.coupling_relaxation <= 1 || error(
        "CHANG_AE_COUPLING_RELAXATION must be in (0,1]",
    )
    analysis.block_duration_s > 0 || error("Moving-block duration must be positive")
    0 <= analysis.block_overlap < 1 || error("Moving-block overlap must be in [0,1)")
    0 <= analysis.frequency_min_hz < analysis.frequency_max_hz || error(
        "Moving-block frequency bounds must satisfy 0 <= min < max",
    )

    run_stamp = Dates.format(now(), "yyyymmdd_HHMMSS")
    output_directory = abspath(get(
        ENV,
        "CHANG_AE_OUTPUT_DIR",
        joinpath(CONVERGENCE_EXAMPLE_DIR, "output", "gated_aeroelastic_$run_stamp"),
    ))
    configuration_path = joinpath(output_directory, "aeroelastic_case_matrix.md")
    dry_run = convergence_bool("CHANG_AE_DRY_RUN", true)
    mkpath(output_directory)
    write_configuration(
        configuration_path,
        aerodynamic_summary,
        selection,
        cases,
        stage,
        analysis,
    )

    println("Gated aeroelastic convergence stage: $stage")
    println("Aerodynamic source: $aerodynamic_summary")
    println("Case matrix: $configuration_path")
    print_case_matrix(cases, analysis.speed_mps, analysis.rpm, analysis.speed_mps, 1.15)
    dry_run && return cases

    lambda_tolerance = convergence_float("CHANG_AE_LAMBDA_TOLERANCE_PER_S", 0.01)
    results = NamedTuple[]
    cache = Dict{Tuple,NamedTuple}()
    for case in cases
        key = case_key(case)
        result = if haskey(cache, key)
            merge(cache[key], (; case, reused = true, elapsed_s = 0.0))
        else
            evaluated = run_case(
                case;
                output_directory,
                speed_mps = analysis.speed_mps,
                trim_rpm = analysis.rpm,
                trim_speed_mps = analysis.speed_mps,
                radius_m = 1.15,
                end_time_s = analysis.end_time_s,
                fit_start_s = analysis.fit_start_s,
                fit_end_s = analysis.fit_end_s,
                minimum_peaks = convergence_int("CHANG_AE_MIN_PEAKS", 8),
                minimum_fit_r_squared = convergence_float("CHANG_AE_MIN_FIT_R2", 0.8),
                moving_block_initial_size = 512,
                moving_block_size_ratio_lower = 0.25,
                moving_block_size_ratio_upper = 0.50,
                moving_block_peak_from_start = 1,
                moving_block_peak_from_end = 0,
                hard_angle_deg = convergence_float("CHANG_AE_ABORT_ANGLE_DEG", 15.0),
                moving_block_duration_s = analysis.block_duration_s,
                moving_block_overlap = analysis.block_overlap,
                moving_block_frequency_min_hz = analysis.frequency_min_hz,
                moving_block_frequency_max_hz = analysis.frequency_max_hz,
                moving_block_apply_hann_window = analysis.hann_window,
                angle_of_attack_deg = analysis.angle_of_attack_deg,
                sideslip_deg = 0.0,
                interaction_on = analysis.interaction_on,
                ga_rho_inf = analysis.ga_rho_inf,
                coupling_relaxation = analysis.coupling_relaxation,
            )
            cache[key] = evaluated
            evaluated
        end
        push!(results, result)
        print_case_result(result)
        partial = annotated_results(results, lambda_tolerance)
        write_summary_csv(joinpath(output_directory, "aeroelastic_convergence.csv"), partial)
        try
            plot_convergence(joinpath(output_directory, "aeroelastic_convergence.png"), partial)
        catch exception
            @warn "Could not refresh partial aeroelastic plot" exception
        end
    end

    results = annotated_results(results, lambda_tolerance)
    summary_path = joinpath(output_directory, "aeroelastic_convergence.csv")
    plot_path = joinpath(output_directory, "aeroelastic_convergence.png")
    write_summary_csv(summary_path, results)
    plot_convergence(plot_path, results)
    println("Aeroelastic summary: $summary_path")
    println("Aeroelastic plot:    $plot_path")
    return results
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_aeroelastic_from_aerodynamic()
end
