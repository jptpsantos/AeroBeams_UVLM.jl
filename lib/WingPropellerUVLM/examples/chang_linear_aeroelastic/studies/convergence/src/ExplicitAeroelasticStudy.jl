# Entry-level explicit-input workflow. Keep it outside the shared numerical
# sources so an aerodynamic study already in progress retains its fingerprint.
mutable struct ExplicitAeroelasticStudy{S} <: WingPropellerUVLMStudies.ConvergenceStudy
    settings::S
    results::Union{Nothing,NamedTuple}
end

function create_study(settings::NamedTuple; kwargs...)
    study = WingPropellerUVLMStudies.create_AeroelasticConvergenceStudy(settings; kwargs...)
    validate_run_mode(WingPropellerUVLMStudies.legacy_settings(study))
    return isnothing(study.settings.selectionFile) ?
           ExplicitAeroelasticStudy(study.settings,nothing) : study
end

function validate_run_mode(s)
    mode = get(s,:run_mode,:full)
    mode in (:preview,:full) || throw(ArgumentError("runMode must be :preview or :full"))
    mode == :preview && !isnothing(s.selection_file) &&
        throw(ArgumentError(
            "Preview uses explicit coarse inputs; use runMode=:full with a verified selection"))
    return nothing
end

explicit_driver_fingerprint() = bytes2hex(Convergence.sha256(read(@__FILE__)))
const LOADED_EXPLICIT_DRIVER_SHA256 = explicit_driver_fingerprint()
struct ExplicitAeroelasticSettings{S}
    movingBlock::S
    driverSha256::String
end
Base.getproperty(s::ExplicitAeroelasticSettings,k::Symbol) =
    k in (:movingBlock,:driverSha256) ?
    getfield(s,k) : getproperty(getfield(s,:movingBlock),k)
Base.propertynames(s::ExplicitAeroelasticSettings) = (propertynames(s.movingBlock)...,:driverSha256)
Convergence.data(s::ExplicitAeroelasticSettings) = merge(Convergence.data(s.movingBlock),
    Dict("input_mode" => "explicit","aerodynamic_verified" => false,
        "explicit_driver_sha256" => s.driverSha256))

function explicit_selection(s)
    for k in (:physical,:core_mode,:nominal)
        hasproperty(s,k) || throw(ArgumentError("Explicit aeroelastic inputs require $k"))
    end
    s.nominal isa NamedTuple || throw(ArgumentError("Explicit mode requires absolute nominal controls"))
    for k in Convergence.FIELDS
        hasproperty(s.nominal,k) && !isnothing(getproperty(s.nominal,k)) ||
            throw(ArgumentError("Explicit mode requires an absolute nominal value for $k"))
    end
    Convergence.validate_controls(s.nominal)
    reference = (; physical=s.physical,core_mode=s.core_mode,selected=Convergence.controls(s.nominal),
        input_mode="explicit",aerodynamic_verified=false)
    return merge(operating_selection(reference,s),(; input_mode="explicit",aerodynamic_verified=false))
end

function Convergence.validate_elastic(s::ExplicitAeroelasticSettings,selection)
    p = selection.physical
    p.wake_shedding_fraction == 0.1 || error("Production aeroelastic solver requires wake_shedding_fraction=0.1")
    nominal = Convergence.elastic_nominal(selection.selected,s)
    # An explicit span sweep rebuilds both meshes; no aerodynamic verification
    # or fixed structural-mesh claim is attached to these cases.
    cases = Convergence.build_cases(nominal,Convergence.elastic_sweeps(nominal,s),s.families)
    d = s.damping
    options = WingPropellerUVLMStudies.moving_block_options(d)
    MovingBlockDamping.moving_block_metrics(Float64[],Float64[]; options...)
    start = Convergence.fit_start(s,p)
    0 <= start < s.response.end_time_s || error("Damping fit must start inside the response interval")
    isnothing(options.block_duration_s) || 0 < options.block_duration_s < s.response.end_time_s - start ||
        error("Damping window must fit inside the post-excitation interval")
    all(x -> isfinite(x) && x >= 0,(d.lambda_tolerance_per_s,d.frequency_tolerance_hz)) ||
        error("Invalid damping tolerances")
    for c in cases
        d.frequency_band_hz[2] < 3p.rpm / c.azimuth_deg || error("Damping frequency band reaches Nyquist")
        defaults = Convergence.model_defaults(p,selection.core_mode,c,s,s.output_directory)
        Convergence.Model.ChangAeroelastic.load_chang_configuration(defaults; env=Dict{String,String}())
    end
    return cases
end

function Convergence.damping_result(history,s::ExplicitAeroelasticSettings,p)
    explicit_driver_fingerprint() == s.driverSha256 || error("Explicit aeroelastic driver changed during the study")
    fits = Convergence.damping_result(history,s.movingBlock,p)
    # A preview remains a solver check even if a user extends its time window.
    return get(s.movingBlock.settings,:run_mode,:full) == :preview ?
           merge(fits,(; valid=false,analysis_mode="preview")) : fits
end

"""Count physical steps and aerodynamic panels; these are not wall-time predictions."""
function runtime_plan(s,selection)
    nominal = Convergence.elastic_nominal(selection.selected,s)
    cases = Convergence.build_cases(nominal,Convergence.elastic_sweeps(nominal,s),s.families)
    seen = Set{Tuple}()
    return map(cases) do c
        dt_s = c.azimuth_deg / (6selection.physical.rpm)
        steps = length(0.0:dt_s:s.response.end_time_s) - 1
        bound_panels = c.wing_span * c.wing_chord + selection.physical.blades * c.prop_radial * c.prop_chord
        wake_rows = ceil(Int,c.wake_revolutions * 360 / c.azimuth_deg)
        wake_panels_capacity = wake_rows * (c.wing_span + selection.physical.blades * c.prop_radial)
        reused = Convergence.key(c) in seen
        push!(seen,Convergence.key(c))
        (; label=c.label,family=c.family,level=c.level,reused,dt_s,steps,
            bound_panels,wake_rows,wake_panels_capacity)
    end
end

function write_runtime_plan(s,selection)
    plan = runtime_plan(s,selection)
    directory = abspath(s.output_directory)
    mkpath(directory)
    open(joinpath(directory,"runtime_plan.csv"),"w") do io
        Convergence.writedlm(io,permutedims(string.(collect(keys(first(plan))))),',')
        for row in plan
            Convergence.writedlm(io,permutedims(collect(values(row))),',')
        end
    end
    unique_cases = filter(row -> !row.reused,plan)
    println("Run mode: $(get(s,:run_mode,:full)); requested duration: $(s.response.end_time_s) s")
    println(
        "Planned work: $(length(unique_cases)) unique cases, $(sum(row.steps for row in unique_cases)) total physical steps"
    )
    println(
        "Bound panels per case: $(minimum(row.bound_panels for row in plan)) to $(maximum(row.bound_panels for row in plan)); wake rows: $(maximum(row.wake_rows for row in plan)) maximum"
    )
    if get(s,:run_mode,:full) == :preview
        println(
            "PREVIEW ONLY: coarse mesh/short wake; damping convergence is disabled. An indeterminate damping fit is expected."
        )
    end
    return plan
end

function run_explicit_aeroelastic(s)
    validate_run_mode(s)
    WingPropellerUVLMStudies.check_loaded_source()
    explicit_driver_fingerprint() == LOADED_EXPLICIT_DRIVER_SHA256 ||
        error("Explicit aeroelastic driver changed after loading; reload the entry file")
    selection = explicit_selection(s)
    println("Explicit aeroelastic inputs: aerodynamic convergence is unverified")
    println("Aeroelastic operating point: $(selection.physical.speed_mps) m/s, $(selection.physical.rpm) RPM")
    println("Aeroelastic coupling: $(s.coupling.scheme); structural integrator: $(s.integration.time_integrator)")
    s.make_plots && !s.dry_run && (@eval WingPropellerUVLMStudies using Plots)
    Convergence.load_model()
    wrapped = ExplicitAeroelasticSettings(MovingBlockStudy(s),explicit_driver_fingerprint())
    # Validate before writing the plan or launching any expensive case.
    Base.invokelatest(Convergence.validate_elastic,wrapped,selection)
    plan = write_runtime_plan(s,selection)
    result = Base.invokelatest(Convergence.run_aeroelastic_loaded,wrapped,selection)
    report = "# Aeroelastic convergence\n\nInput mode: explicit. Aerodynamic convergence is **unverified**; no aerodynamic selection was loaded or certified.\n\n"
    report *= "Operating point: $(selection.physical.speed_mps) m/s, $(selection.physical.rpm) RPM.\n\n"
    report *= "Aeroelastic coupling: `$(s.coupling.scheme)`. Structural integrator: `$(s.integration.time_integrator)`.\n\n"
    if get(s,:run_mode,:full) == :preview
        report *= "**Preview only:** this coarse, short run checks solver execution and coupling. Damping convergence is disabled; an indeterminate fit is expected. Use the full mode for damping sensitivity.\n\n"
    end
    report *= "Requested duration: $(s.response.end_time_s) s. Planned physical steps across unique cases: $(sum(row.steps for row in plan if !row.reused)). See runtime_plan.csv for time steps, bound panels and allocated wake capacity; these counts do not predict wall time.\n\n"
    report *= "Nominal numerical controls:\n\n```toml\n$(Convergence.toml_text(result.nominal))```\n\n"
    if :wing_span in s.families
        report *= "The wing span-panel sweep changes both the aerodynamic and structural meshes.\n\n"
    end
    report *= "Positive lambda means growth; negative lambda means decay. Compare pitch/yaw growth rates and frequencies, and inspect damping_fit.toml and raw histories. Damping convergence does not establish aerodynamic convergence or a flutter boundary.\n"
    write(joinpath(result.directory,"report.md"),report)
    return merge(result,(; aerodynamic_reference=nothing))
end

function WingPropellerUVLMStudies.solve!(study::ExplicitAeroelasticStudy)
    study.results = run_explicit_aeroelastic(WingPropellerUVLMStudies.legacy_settings(study))
    return study
end
