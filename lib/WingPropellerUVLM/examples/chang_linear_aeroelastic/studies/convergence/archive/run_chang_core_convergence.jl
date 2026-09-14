# Run from the repository root with --project=lib/WingPropellerUVLM.
module ChangCoreConvergence
using Dates, DelimitedFiles, TOML
module CaseModel
include(joinpath(@__DIR__, "..", "..", "..", "src", "ChangAeroelastic.jl"))
end
module Rigid
include(joinpath(@__DIR__, "run_chang_coupled_aerodynamic_analysis.jl"))
end
const Integrity = Rigid.ChangAerodynamicStudy
include("core_convergence_case.jl")

const MESH_FIELDS = (:wing_span, :wing_chord, :prop_radial, :prop_chord)
panel_count(x) = x isa Integer && !(x isa Bool) && x > 0
mesh_controls(mesh) = NamedTuple{MESH_FIELDS}(Tuple(Int(getproperty(mesh,k)) for k in MESH_FIELDS))
effective_core_mesh(mesh,s) = merge(mesh_controls(mesh),
    (; wing_span=s.mode == :aeroelastic ? Int(s.nominal_mesh.wing_span) : Int(mesh.wing_span)))

function validate_mesh(mesh, name)
    all(k -> hasproperty(mesh,k) && panel_count(getproperty(mesh,k)), MESH_FIELDS) ||
        error("$name must specify positive integer wing_span, wing_chord, prop_radial and prop_chord counts")
end

function validate_study(s, config)
    s.mode in (:aerodynamic, :aeroelastic) || error("mode must be :aerodynamic or :aeroelastic")
    for (name, values, decreasing) in (
        ("radii_m", s.radii_m, true),
        ("azimuth_levels_deg", s.azimuth_levels_deg, true),
        ("wake_levels_revolutions", s.wake_levels_revolutions, false))
        length(values) >= 3 || error("$name requires at least three ordered levels")
        all(x -> isfinite(x) && x > 0, values) || error("$name must be finite and positive")
        all(x -> decreasing ? x < 0 : x > 0, diff(values)) || error("$name must be strictly ordered")
    end
    validate_mesh(s.nominal_mesh,"nominal_mesh")
    for k in MESH_FIELDS
        hasproperty(s.panel_levels,k) || error("Missing panel_levels.$k")
        levels = getproperty(s.panel_levels,k)
        length(levels) >= 3 && all(panel_count,levels) && all(>(0),diff(levels)) ||
            error("panel_levels.$k requires at least three strictly increasing positive integer counts")
    end
    length(s.core_meshes) >= 3 || error("core_meshes requires at least three meshes")
    for mesh in s.core_meshes
        validate_mesh(mesh,"core_meshes entry")
    end
    effective = effective_core_mesh.(s.core_meshes,Ref(s))
    for i in 2:length(effective)
        changes = [getproperty(effective[i],k)-getproperty(effective[i-1],k) for k in MESH_FIELDS]
        all(>=(0),changes) && any(>(0),changes) ||
            error("Each effective core mesh must refine at least one direction without coarsening others; aeroelastic wing span is fixed")
    end
    isfinite(s.nominal_radius_m) && s.nominal_radius_m > 0 || error("Nominal radius must be positive")
    isfinite(s.nominal_wake_revolutions) && s.nominal_wake_revolutions > 0 || error("Wake age must be positive")
    allowed = (:core_mesh, :wing_span, :wing_chord, :prop_radial, :prop_chord, :time_step, :wake_length)
    !isempty(s.families) && all(f -> f in allowed, s.families) || error("Unknown or empty study families")
    length(unique(s.families)) == length(s.families) || error("Duplicate study families")
    length(config.propeller.attachment_eta) == 1 || error("This study supports one propeller")
    config.simulation.near_field_force_model == :imperial ||
        error("The rigid reference uses :imperial loads; select that model for this comparison")
    all(x -> isfinite(x) && x >= 0, (s.coefficient_absolute_tolerances...,
        s.coefficient_relative_tolerance, s.periodic_tolerances...,
        s.lambda_tolerance_per_s, s.frequency_tolerance_hz)) || error("Invalid comparison tolerance")
    if s.mode == :aeroelastic
        start = fit_start(s, config)
        0 <= start < config.simulation.end_time_s || error("Damping fit must start inside the response")
        0 < s.block_duration_s < config.simulation.end_time_s - start || error("FFT window must fit the response interval")
        0 <= s.block_overlap < 1 || error("Overlap must be in [0,1)")
        0 <= s.minimum_fit_r_squared <= 1 || error("Fit R2 threshold must be in [0,1]")
        s.minimum_peaks >= 3 || error("At least three fit points are required")
        lo, hi = s.frequency_band_hz
        isfinite(lo) && isfinite(hi) && 0 <= lo < hi || error("Invalid modal frequency band")
    end
    return s
end

rotor_rpm(config) = config.propeller.rotation_rpm *
    config.simulation.freestream_speed_mps / config.propeller.trim_speed_mps
fit_start(s, config) = isnothing(s.fit_start_s) ?
    max(0.75, config.excitation.impulse_start_s + config.excitation.impulse_duration_s + 60 / rotor_rpm(config)) : s.fit_start_s

function build_cases(s, config)
    validate_study(s, config)
    base = merge(mesh_controls(s.nominal_mesh), (; radius_m = Float64(s.nominal_radius_m),
        azimuth_deg = config.simulation.azimuth_step_deg, wake_revolutions = s.nominal_wake_revolutions))
    cases = NamedTuple[]
    add(family, level, controls) = push!(cases, merge(base, controls,
        (; family, level, label = "$(family)_L$(level)")))
    for family in s.families
        if family == :core_mesh
            for (mesh, counts) in enumerate(s.core_meshes), (level, radius) in enumerate(s.radii_m)
                add(Symbol("core_mesh_$mesh"), level, merge(effective_core_mesh(counts,s),
                    (; radius_m = Float64(radius))))
            end
            # Reuse the same simulations to quantify mesh convergence at each
            # fixed radius, as well as core sensitivity on each fixed mesh.
            core_cases = copy(cases)
            for (radius_level, radius) in enumerate(s.radii_m), mesh_level in eachindex(s.core_meshes)
                source = only(filter(c -> c.family == Symbol("core_mesh_$mesh_level") &&
                    c.radius_m == radius, core_cases))
                add(Symbol("mesh_at_core_$radius_level"), mesh_level, source)
            end
        elseif family == :time_step
            for (level, value) in enumerate(s.azimuth_levels_deg)
                add(family, level, (; azimuth_deg = value))
            end
        elseif family == :wake_length
            for (level, value) in enumerate(s.wake_levels_revolutions)
                add(family, level, (; wake_revolutions = value))
            end
        else
            family == :wing_span && s.mode == :aeroelastic && continue
            for (level, count) in enumerate(getproperty(s.panel_levels,family))
                add(family, level, NamedTuple{(family,)}((Int(count),)))
            end
        end
    end
    isempty(cases) && error("No cases remain for the selected mode")
    for c in cases
        steps = 360 / c.azimuth_deg
        isfinite(steps) && steps >= 4 && isapprox(steps, round(steps); atol=1e-10, rtol=0) ||
            error("Every azimuth step must divide 360 degrees and be at most 90 degrees")
        Integrity.validate_options(rigid_options(c, s, config))
        s.mode == :aeroelastic && s.frequency_band_hz[2] >= 0.5 / time_step(c, config) &&
            error("Modal frequency band reaches the Nyquist frequency")
    end
    return cases
end

case_key(c) = (c.radius_m, c.wing_span, c.wing_chord, c.prop_radial, c.prop_chord,
    c.azimuth_deg, c.wake_revolutions)
time_step(c, config) = c.azimuth_deg / (6 * rotor_rpm(config))
wake_rows(c) = ceil(Int, c.wake_revolutions * 360 / c.azimuth_deg)

function rigid_options(c, s, config)
    w, p, sim, aero = config.wing, config.propeller, config.simulation, config.aerodynamic
    return Integrity.ChangCoupledAerodynamicOptions(
        flow_speed_mps=sim.freestream_speed_mps, air_density_kgpm3=sim.air_density_kgpm3,
        angle_of_attack_deg=sim.angle_of_attack_deg, sideslip_deg=sim.sideslip_deg,
        wing_span_m=w.span_m, wing_root_chord_m=w.root_chord_m, wing_tip_chord_m=w.tip_chord_m,
        wing_span_panels=c.wing_span, wing_chord_panels=c.wing_chord,
        wing_symmetric=w.symmetric,
        elastic_axis_fraction=aero.elastic_axis_fraction, propeller_radius_m=p.radius_m,
        propeller_chord_m=p.chord_m, propeller_blades=p.blades,
        propeller_radial_panels=c.prop_radial, propeller_chord_panels=c.prop_chord,
        propeller_attachment_eta=only(p.attachment_eta), pylon_length_m=config.structural.pylon_length_m,
        reference_rpm=p.rotation_rpm, reference_speed_mps=p.trim_speed_mps,
        collective_pitch_offset_deg=p.collective_pitch_offset_deg,
        azimuth_step_deg=c.azimuth_deg, simulated_revolutions=s.simulated_revolutions,
        averaged_revolutions=s.averaged_revolutions, retained_wake_revolutions=c.wake_revolutions,
        core_radius_m=c.radius_m, finite_core_segment_factor=0, finite_core_chord_factor=0,
        interaction_on=sim.interaction_on)
end

function rigid_run(c, s, config, directory)
    options = rigid_options(c, s, config)
    mkpath(directory)
    log = joinpath(directory, "run.log")
    reused = s.reuse_existing && Integrity.output_matches_options(directory, options; warn_on_mismatch=false)
    if !reused
        open(log, "w") do io
            redirect_stdout(io) do
                redirect_stderr(io) do
                    result = Rigid.simulate_coupled_aerodynamics(options)
                    Rigid.write_results(result, directory; plot_results=false)
                end
            end
        end
    end
    # Validate hashes and options again even for freshly written histories.
    Integrity.output_matches_options(directory, options; warn_on_mismatch=false) || error("Result integrity check failed")
    history = joinpath(directory, "coupled_aerodynamic_history.csv")
    raw, header = readdlm(history, ',', header=true)
    names = String.(vec(header))
    column(name) = Float64.(raw[:, something(findfirst(==(name), names))])
    steps = Integrity.validate_options(options)
    length(column("time_s")) == steps * s.simulated_revolutions || error("Incomplete history")
    all(isapprox.(column("time_s"), (1:size(raw,1)) .* time_step(c, config); rtol=1e-10, atol=1e-12)) || error("Invalid time grid")
    expected = mod.((1:size(raw,1)) .* c.azimuth_deg, 360)
    all(abs.(mod.(column("azimuth_deg") .- expected .+ 180, 360) .- 180) .< 1e-8) || error("Invalid phase grid")
    metrics = map(name -> Integrity.periodic_metrics(column(name), steps, s.averaged_revolutions),
        ("wing_CL", "propeller_CT", "propeller_CQ"))
    return (; status="completed", reason="verified periodic coefficient history", reused,
        values=map(m -> m.mean, metrics), phases=map(m -> m.phase, metrics),
        periodic=map(m -> m.periodic_rms, metrics),
        valid=all(i -> metrics[i].periodic_rms <= s.periodic_tolerances[i], 1:3), history, log)
end

function load_elastic_backend()
    if !isdefined(@__MODULE__, :Elastic)
        @eval module Elastic
            include(joinpath(@__DIR__, "chang_aeroelastic_convergence.jl"))
        end
    end
end

function elastic_run(c, s, config, directory)
    mkpath(directory)
    case = Elastic.UVLMConvergenceCase(family=c.family, level=c.level, label=c.label,
        wing_span_panels=c.wing_span, wing_chord_panels=c.wing_chord,
        prop_radial_panels=c.prop_radial, prop_chord_panels=c.prop_chord,
        wake_revolutions=c.wake_revolutions, core_radius_m=c.radius_m,
        fcore_segment_factor=0, fcore_chord_factor=0, azimuth_step_deg=c.azimuth_deg)
    e = config.excitation
    # The child reads the unchanged case file. Override the legacy driver's
    # excitation defaults explicitly, and discard unrelated session overrides.
    environment = Pair{String,Union{Nothing,String}}[
        name => nothing for name in keys(ENV) if startswith(name, "CHANG_")]
    append!(environment, [
        "CHANG_CONVERGENCE_TRIM_REVOLUTIONS" => string(e.trim_revolutions),
        "CHANG_CONVERGENCE_TRIM_AVERAGE_REVOLUTIONS" => string(e.trim_average_revolutions),
        "CHANG_CONVERGENCE_IMPULSE_MAGNITUDE" => string(e.impulse_magnitude_nm),
        "CHANG_CONVERGENCE_IMPULSE_START_S" => string(e.impulse_start_s),
        "CHANG_CONVERGENCE_IMPULSE_DURATION_S" => string(e.impulse_duration_s),
        "CHANG_CONVERGENCE_FORCE_MODEL" => string(config.simulation.near_field_force_model),
        "CHANG_CONVERGENCE_PROP_MOMENT_PROJECTION" => string(config.simulation.propeller_moment_projection)])
    result = withenv(environment...) do
        Elastic.run_case(case; output_directory=directory,
            speed_mps=config.simulation.freestream_speed_mps, trim_rpm=config.propeller.rotation_rpm,
            trim_speed_mps=config.propeller.trim_speed_mps, radius_m=config.propeller.radius_m,
            end_time_s=config.simulation.end_time_s, fit_start_s=fit_start(s, config),
            fit_end_s=config.simulation.end_time_s, minimum_peaks=s.minimum_peaks,
            minimum_fit_r_squared=s.minimum_fit_r_squared, moving_block_initial_size=512,
            moving_block_size_ratio_lower=0.25, moving_block_size_ratio_upper=0.50,
            moving_block_peak_from_start=1, moving_block_peak_from_end=0,
            hard_angle_deg=config.integration.propeller_angle_limit_deg,
            moving_block_duration_s=s.block_duration_s, moving_block_overlap=s.block_overlap,
            moving_block_frequency_min_hz=s.frequency_band_hz[1],
            moving_block_frequency_max_hz=s.frequency_band_hz[2], moving_block_apply_hann_window=true,
            angle_of_attack_deg=config.simulation.angle_of_attack_deg,
            sideslip_deg=config.simulation.sideslip_deg, interaction_on=config.simulation.interaction_on,
            wing_symmetric=config.wing.symmetric,
            ga_rho_inf=config.integration.rho_inf, coupling_relaxation=config.coupling.relaxation)
    end
    write_fit_diagnostics(directory,result)
    return (; status=result.status, reason=result.reason, reused=false,
        values=(result.pitch.moving_block_lambda_per_s, result.yaw.moving_block_lambda_per_s,
                result.pitch.frequency_hz, result.yaw.frequency_hz), phases=(), periodic=(),
        valid=result.status == "completed" && result.pitch.valid_for_convergence && result.yaw.valid_for_convergence,
        history=result.history_path, log=result.log_path)
end

function write_fit_diagnostics(directory,result)
    open(joinpath(directory,"damping_fit.toml"),"w") do io
        TOML.print(io, Dict("status"=>result.status,"reason"=>result.reason,
            "pitch"=>Dict(string(k)=>v for (k,v) in pairs(result.pitch)),
            "yaw"=>Dict(string(k)=>v for (k,v) in pairs(result.yaw))); sorted=true)
    end
end

function compare_results(results, cases, s)
    rows = NamedTuple[]
    for family in unique(c.family for c in cases)
        group = sort(filter(r -> r.case.family == family, results); by=r -> r.case.level)
        isempty(group) && continue
        expected = count(c -> c.family == family, cases)
        reference = last(group)
        limits = s.mode == :aerodynamic ? ntuple(i -> max(s.coefficient_absolute_tolerances[i],
            s.coefficient_relative_tolerance * abs(reference.values[i])), 3) :
            (s.lambda_tolerance_per_s, s.lambda_tolerance_per_s,
             s.frequency_tolerance_hz, s.frequency_tolerance_hz)
        errors(a, b) = map(eachindex(limits)) do i
            mean_error = abs(a.values[i] - b.values[i])
            s.mode == :aerodynamic ? max(mean_error, Integrity.phase_rms_error(a.phases[i], b.phases[i])) : mean_error
        end
        agrees(a,b) = a.valid && b.valid && all(errors(a,b) .<= collect(limits))
        complete = length(group) == expected && [r.case.level for r in group] == collect(1:expected) &&
            all(r -> r.status == "completed", group)
        stable = complete && length(group) >= 2 && agrees(group[end-1], reference)
        for (i, row) in enumerate(group)
            err = row.valid && reference.valid ? errors(row, reference) : fill(NaN, length(limits))
            previous = i > 1 && row.valid && group[i-1].valid ? errors(row,group[i-1]) : fill(NaN,length(limits))
            accepted = stable && all(r -> agrees(r,reference), group[i:end])
            push!(rows, merge(row, (; reference=reference.case.label, complete, reference_stable=stable,
                accepted, errors=err, previous_errors=previous, limits)))
        end
    end
    return rows
end

metric_names(s) = s.mode == :aerodynamic ? ("CL", "CT", "CQ") :
    ("pitch_lambda_per_s", "yaw_lambda_per_s", "pitch_frequency_hz", "yaw_frequency_hz")

function write_tables(directory, rows, cases, s, config)
    open(joinpath(directory,"case_matrix.csv"), "w") do io
        writedlm(io, permutedims(["family","level","label","core_radius_m","wing_span","wing_chord",
            "prop_radial","prop_chord","azimuth_deg","dt_s","wake_revolutions","wake_rows"]), ',')
        for c in cases
            writedlm(io, permutedims(Any[c.family,c.level,c.label,c.radius_m,c.wing_span,c.wing_chord,
                c.prop_radial,c.prop_chord,c.azimuth_deg,time_step(c,config),c.wake_revolutions,wake_rows(c)]), ',')
        end
    end
    names = metric_names(s)
    open(joinpath(directory,"convergence.csv"), "w") do io
        headers = vcat(["family","level","label","core_radius_m","status","reason","reused","valid",
            "family_complete","reference_stable","accepted","reference"], collect(names),
            ["error_$n" for n in names], ["previous_change_$n" for n in names],
            ["limit_$n" for n in names], s.mode == :aerodynamic ? ["periodic_$n" for n in names] : String[],
            ["elapsed_s","history_path","log_path"])
        writedlm(io, permutedims(headers), ',')
        for r in rows
            values = Any[r.case.family,r.case.level,r.case.label,r.case.radius_m,r.status,r.reason,r.reused,
                r.valid,r.complete,r.reference_stable,r.accepted,r.reference,r.values...,r.errors...,
                r.previous_errors...,r.limits...,r.periodic...,r.elapsed_s,r.history,r.log]
            writedlm(io, permutedims(values), ',')
        end
    end
    open(joinpath(directory,"report.md"), "w") do io
        println(io,"# Chang fixed-core convergence\n\nMode: $(s.mode). Speed: $(config.simulation.freestream_speed_mps) m/s; RPM: $(rotor_rpm(config)); interaction: $(config.simulation.interaction_on).")
        println(io,"\nWing image symmetry across y=0: $(config.wing.symmetric); propeller blade symmetry: false.")
        println(io,"\n$(length(cases)) family entries; $(length(unique(case_key.(cases)))) unique cases. Radii are in metres. Both wakes retain the listed age; rows scale with time step.")
        println(io,"\nThe last listed level is a comparison reference, not an exact solution. Acceptance requires a complete family, valid candidate/reference signals, agreement of the final pair, and agreement of all finer levels. Core-family acceptance indicates a sensitivity plateau, not a calibrated physical radius. The mesh_at_core families reuse those simulations to check mesh convergence at each fixed radius.")
        println(io,"\nAeroelastic mode fixes wing span/structural resolution. Positive lambda indicates growth. Inspect the raw responses and tracked frequencies for mode switching and beating; a fitted growth rate alone does not establish a flutter boundary. Rigid periodic coefficients do not establish converged damping.")
        println(io,"\n| Family | Level | Radius (m) | Status | Valid | Accepted |\n|---|---:|---:|---|---|---|")
        for r in rows
            println(io,"| $(r.case.family) | $(r.case.level) | $(r.case.radius_m) | $(r.status) | $(r.valid) | $(r.accepted) |")
        end
        println(io,"\nSee case_matrix.csv, convergence.csv and each case's history/log. No retrimming is performed: RPM, collective and physical conditions are held fixed. Use matched-condition retrimming separately if required.")
    end
end

function plot_results(directory, rows, s)
    names = metric_names(s)
    panels = Any[]
    core_families = filter(f -> startswith(string(f), "core_mesh_"), unique(r.case.family for r in rows))
    for (i, name) in enumerate(names)
        radii = sort(s.radii_m)
        p = plot(; xlabel="Core radius (m)", ylabel=name, xscale=:log10, legend=:best,
            xticks=(radii,string.(radii)), left_margin=8 * Plots.mm)
        for family in core_families
            group = sort(filter(r -> r.case.family == family && r.status == "completed", rows); by=r -> r.case.radius_m)
            isempty(group) && continue
            plot!(p, [r.case.radius_m for r in group], [r.values[i] for r in group]; marker=:circle, label=string(family))
        end
        push!(panels,p)
    end
    savefig(plot(panels...; layout=(length(panels),1), size=(950,300length(panels))), joinpath(directory,"core_sensitivity.png"))
    panels = Any[]
    for family in unique(r.case.family for r in rows)
        group = sort(filter(r -> r.case.family == family, rows); by=r -> r.case.level)
        valid_comparison = any(r -> any(isfinite,r.errors),group)
        title = string(family) * (valid_comparison ? "" : "\nNo valid comparison")
        p = plot(; title, xlabel="Refinement level (see case_matrix.csv)", ylabel="Error / tolerance",
            left_margin=8 * Plots.mm, xticks=[r.case.level for r in group])
        for (i, name) in enumerate(names)
            ratio(r) = r.limits[i] > 0 ? r.errors[i]/r.limits[i] : (r.errors[i] == 0 ? 0.0 : NaN)
            plot!(p, [r.case.level for r in group], ratio.(group); marker=:circle, label=name)
        end
        hline!(p,[1.0]; color=:black, linestyle=:dash, label="tolerance")
        push!(panels,p)
    end
    savefig(plot(panels...; layout=(cld(length(panels),2),2), size=(1200,320cld(length(panels),2))),
        joinpath(directory,"convergence_errors.png"))
end

function run_study(s, config)
    cases = build_cases(s, config)
    directory = isnothing(s.output_directory) ? joinpath(@__DIR__, "..", "..", "..","output",
        "core_convergence_$(s.mode)_$(Dates.format(now(),"yyyymmdd_HHMMSS"))") : s.output_directory
    directory = abspath(directory)
    mkpath(directory)
    fingerprint = Integrity.source_fingerprint()
    open(joinpath(directory,"study_snapshot.txt"), "w") do io
        println(io,"Source SHA256: $fingerprint\nJulia: $VERSION\nSettings:\n$(repr(s))\nResolved case:\n$(repr(config))")
    end
    write_tables(directory, NamedTuple[], cases, s, config)
    println("$(s.mode): $(length(cases)) entries, $(length(unique(case_key.(cases)))) unique simulations")
    println("Case matrix: $(joinpath(directory,"case_matrix.csv"))")
    s.mode == :aeroelastic && println("Wing span/structural mesh held fixed; inspect fit band and response duration.")
    s.dry_run && return (; directory, cases, results=NamedTuple[])
    s.mode == :aeroelastic && load_elastic_backend()
    if s.make_plots
        @eval using Plots
    end
    results = NamedTuple[]
    cache = Dict{Tuple,NamedTuple}()
    for (index,c) in enumerate(cases)
        Integrity.source_fingerprint() == fingerprint || error("Source or case changed during study; restart into a new directory")
        key = case_key(c)
        start = time()
        println("[$index/$(length(cases))] $(c.label): rc=$(c.radius_m) m, wing=$(c.wing_span)x$(c.wing_chord), blade=$(c.prop_radial)x$(c.prop_chord), wake=$(wake_rows(c)) rows")
        result = if haskey(cache,key)
            merge(cache[key], (; reused=true))
        else
            folder = joinpath(directory,c.label)
            r = try
                s.mode == :aerodynamic ? rigid_run(c,s,config,folder) : Base.invokelatest(elastic_run,c,s,config,folder)
            catch exception
                exception isa InterruptException && rethrow()
                @warn "Case failed; continuing study" label=c.label exception
                n = length(metric_names(s))
                (; status="failed",reason=sprint(showerror,exception),reused=false,values=ntuple(_->NaN,n),
                    phases=s.mode == :aerodynamic ? (Float64[],Float64[],Float64[]) : (),
                    periodic=s.mode == :aerodynamic ? (NaN,NaN,NaN) : (), valid=false,
                    history="",log=joinpath(folder,"run.log"))
            end
            cache[key] = r
            r
        end
        Integrity.source_fingerprint() == fingerprint || error("Source or case changed during a simulation; result not accepted")
        push!(results,merge(result,(; case=c,elapsed_s=time()-start)))
        rows = compare_results(results,cases,s)
        write_tables(directory,rows,cases,s,config)
        println("  $(result.status); signal valid=$(result.valid); values=$(result.values)")
    end
    rows = compare_results(results,cases,s)
    s.make_plots && Base.invokelatest(plot_results,directory,rows,s)
    println("Results: $directory")
    return (; directory,cases,results=rows)
end

function main()
    # Resolve from the editable case, without unrelated interactive ENV overrides.
    config = CaseModel.ChangAeroelastic.load_chang_configuration(; env=Dict{String,String}())
    return run_study(core_convergence_defaults(config),config)
end
end

if abspath(PROGRAM_FILE) == @__FILE__
    ChangCoreConvergence.main()
end
