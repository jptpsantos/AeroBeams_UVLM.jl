# Staged studies use the production time marcher and the optional threaded
# wake commit. Keep this orchestration outside the shared solver fingerprint.
using LinearAlgebra: BLAS
using Printf: @sprintf

module StagedPitchFit
include(joinpath(@__DIR__,"..","..","..","analyze_pitch_phase.jl"))
end

struct StagedDeadline <: Exception end

function staged_plan(s)
    nominal = Convergence.controls(s.nominal)
    sweeps = Convergence.elastic_sweeps(nominal,s)
    cases = Convergence.build_cases(nominal,sweeps,s.families)
    pilot = merge(nominal,(; family=:pilot,level=0,label="pilot",stage=:pilot))
    screen = NamedTuple[]
    followup = NamedTuple[]
    for family in s.families
        group = filter(c -> c.family == family,cases)
        atNominal = findfirst(c -> Convergence.key(c) == Convergence.key(pilot),group)
        isnothing(atNominal) && error("Include the nominal $family value in its staged sweep")
        atNominal < length(group) || error("Include a finer $family value in its staged sweep")
        if family == :azimuth_deg && atNominal > 1
            push!(screen,merge(group[atNominal-1],(; stage=:coarser)))
        end
        push!(screen,merge(group[atNominal+1],(; stage=:screen)))
        append!(followup,(merge(c,(; stage=:followup)) for c in group[atNominal+2:end]))
    end
    # Start the time-step comparison early. Give every other control its
    # screening run before the more expensive fine-step calculation.
    sort!(screen; by=c -> c.stage == :coarser ? 0 : c.family == :azimuth_deg ? 3 :
        c.family == :wake_revolutions ? 2 : 1)
    return pilot,screen,followup
end

function staged_work(c,p,s)
    dt = c.azimuth_deg / (6p.rpm)
    steps = length(0.0:dt:s.response.end_time_s) - 1
    rows = ceil(Int,c.wake_revolutions * 360 / c.azimuth_deg)
    bound = c.wing_span * c.wing_chord + p.blades * c.prop_radial * c.prop_chord
    wake = rows * (c.wing_span + p.blades * c.prop_radial)
    return (; dt_s=dt,steps,bound_panels=bound,wake_rows=rows,wake_panels_capacity=wake)
end

function staged_fit(history,s,p)
    signal = StagedPitchFit.read_pitch_history(history)
    stop = min(last(signal.time),s.response.end_time_s)
    firstFit = Convergence.fit_start(s,p)
    # All runs use identical physical windows to within one discrete step.
    fitStarts = unique(max.(firstFit,s.response.end_time_s .- (2.4,1.95,1.55)))
    length(fitStarts) >= 2 && last(fitStarts) < stop ||
        error("The staged response needs two distinct free-response fit windows")
    fits = map(fitStarts) do start
        m = StagedPitchFit.identify_pitch_mode(signal.time,signal.wing,signal.pylon;
            fit_start=start,fit_end=stop)
        return (; start_s=start,end_s=stop,lambda_per_s=m.growth_rate_per_s,frequency_hz=m.frequency_hz,
            residuals=m.normalized_residuals,reliable=m.reliable,reasons=m.reasons)
    end
    main = first(fits)
    lambdaSpread = maximum(abs(m.lambda_per_s - main.lambda_per_s) for m in fits)
    frequencySpread = maximum(abs(m.frequency_hz - main.frequency_hz) for m in fits)
    lo,hi = s.damping.frequency_band_hz
    valid = all(m -> m.reliable && lo <= m.frequency_hz <= hi,fits) &&
        lambdaSpread <= s.damping.lambda_tolerance_per_s / 2 &&
        frequencySpread <= s.damping.frequency_tolerance_hz / 2
    return (; valid,lambda_per_s=main.lambda_per_s,frequency_hz=main.frequency_hz,
        lambda_window_shift_per_s=lambdaSpread,frequency_window_shift_hz=frequencySpread,windows=collect(fits))
end

function staged_agrees(a,b,s)
    a.status == "completed" && b.status == "completed" || return false
    return abs(a.lambda_per_s - b.lambda_per_s) <= s.damping.lambda_tolerance_per_s &&
        abs(a.frequency_hz - b.frequency_hz) <= s.damping.frequency_tolerance_hz &&
        sign(a.lambda_per_s) == sign(b.lambda_per_s)
end

function staged_march(config,deadline)
    C = Convergence.Model.ChangAeroelastic
    F = Convergence.Model.ChangRuntimeOptimizations
    model = C.build_chang_model(config)
    config = model.config
    p = model.parameters
    workspace = F.build_fast_workspace(model; save_surface_history=false)
    wake = F.fast_wake_context(workspace; interaction_on=config.simulation.interaction_on)
    excitation = C.chang_excitation_options(config,p)
    integration = C.chang_integration_options(config,model.structural,p)
    C.report_chang_solver_options(excitation,integration)
    animation = C.chang_animation_options(workspace.system,workspace.iwake,config.output)
    aerodynamic_load = (snapshot,state,step) -> begin
        time() < deadline || throw(StagedDeadline())
        C.aero_load_for_state!(model,workspace,snapshot,state,step)
    end
    solution = C.solve_chang_aeroelastic!(workspace.system,model.structural;
        time=p.t,time_steps=p.dt,freestream_history=workspace.fs_vec,aerodynamic_load,
        wake,excitation,integration,animation)
    solution.last_step == length(p.dt) && all(solution.coupling_converged) ||
        error("Aeroelastic case stopped before completing its accepted time grid")
    C.write_chang_results(solution;
        wing=config.wing,propeller=config.propeller,simulation=config.simulation,
        time=p.t,time_steps=p.dt,wing_node_count=p.nnodes,dofs_per_node=p.ndof,
        density=p.ref.rho,visualization=config.output,output_directory=config.output.directory,
        output_label=config.output.label)
    return nothing
end

staged_history_digest(path) = open(io -> bytes2hex(Convergence.sha256(io)),path)

function staged_cache_record(path,history,id)
    isfile(path) && isfile(history) || return nothing
    try
        record = Convergence.TOML.parsefile(path)
        record["identity"] == id && record["history_sha256"] == staged_history_digest(history) ||
            return nothing
        elapsed = record["solver_elapsed_s"]
        isfinite(elapsed) && elapsed > 0 || return nothing
        return elapsed
    catch exception
        exception isa InterruptException && rethrow()
        return nothing
    end
end

function write_staged_table(path,rows,columns)
    open(path * ".tmp","w") do io
        Convergence.writedlm(io,permutedims(string.(collect(columns))),',')
        for row in rows
            Convergence.writedlm(io,permutedims(Any[getproperty(row,k) for k in columns]),',')
        end
    end
    mv(path * ".tmp",path; force=true)
end

function staged_case(c,s,selection,directory,fingerprint,deadline,timings)
    p = selection.physical
    defaults = Convergence.model_defaults(p,selection.core_mode,c,s,directory)
    C = Convergence.Model.ChangAeroelastic
    config = C.load_chang_configuration(defaults; env=Dict{String,String}())
    identity = Convergence.digest((; fingerprint,backend="exact_threaded_wake",
        config=Base.structdiff(config,(; output=nothing))))
    folder = joinpath(directory,"cases",identity[1:20])
    config = merge(config,(; output=merge(config.output,(; directory=folder))))
    history = joinpath(folder,"response_history.csv")
    log = joinpath(folder,"run.log")
    cachePath = joinpath(folder,"completed.toml")
    elapsed = staged_cache_record(cachePath,history,identity)
    reused = !isnothing(elapsed)
    work = staged_work(c,p,s)
    # Estimate queue fit from the nearest finished workload on this machine.
    # The step-level deadline remains authoritative if the estimate is low.
    workload = work.steps * (work.bound_panels + work.wake_panels_capacity)^2
    if !reused
        remaining = deadline - time()
        estimate = isempty(timings) ? 0.0 : begin
            nearest = argmin(t -> abs(log2(t.workload / workload)),timings)
            1.2 * nearest.elapsed * (workload / nearest.workload)
        end
        remaining > 0 && estimate <= remaining ||
            return (; status="deferred",reason="Estimated workload exceeds the remaining wall budget",
                lambda_per_s=NaN,frequency_hz=NaN,lambda_window_shift_per_s=NaN,
                elapsed_s=0.0,reused=false,history="",log="")
        mkpath(folder)
        isfile(cachePath) && rm(cachePath)
        Convergence.save_toml(joinpath(folder,"resolved_configuration.toml"),config)
        println("  $(work.steps) steps, $(work.wake_rows) wake rows, log: $log")
        flush(stdout)
        start = time()
        try
            open(log,"w") do io
                redirect_stdout(io) do
                    redirect_stderr(io) do
                        staged_march(config,deadline)
                    end
                end
            end
        catch exception
            exception isa StagedDeadline || rethrow()
            return (; status="deferred",reason="Wall budget reached during this case",
                lambda_per_s=NaN,frequency_hz=NaN,lambda_window_shift_per_s=NaN,
                elapsed_s=time() - start,reused=false,history="",log)
        end
        elapsed = time() - start
        Convergence.SOURCE() == fingerprint || error("Solver source changed during the case")
        Convergence.save_toml(cachePath,(; identity,solver_elapsed_s=elapsed,
            history_sha256=staged_history_digest(history)))
    end
    push!(timings,(; workload,elapsed))
    fit = staged_fit(history,s,p)
    Convergence.save_toml(joinpath(folder,"damping_fit.toml"),fit)
    return (; status=fit.valid ? "completed" : "indeterminate",
        reason=fit.valid ? "" : "Mode fit or fit-window agreement failed",
        fit.lambda_per_s,fit.frequency_hz,fit.lambda_window_shift_per_s,
        elapsed_s=reused ? 0.0 : elapsed,reused,history,log)
end

function staged_report(directory,selection,s,status,pairs)
    open(joinpath(directory,"report.md"),"w") do io
        println(io,"# Chang staged damping sensitivity\n")
        println(io,"$(selection.physical.speed_mps) m/s, $(selection.physical.rpm) RPM, $(s.response.end_time_s) s response.")
        println(io,"Aerodynamic interaction: $(selection.physical.interaction). Julia threads: $(Threads.nthreads()).")
        println(io,"Status: $status\n")
        println(io,"The fit tracks a common mode in wing-tip twist and relative pylon pitch. Positive growth rate indicates growth.")
        println(io,"Treat the stability sign as unresolved if the rate magnitude falls within the numerical tolerance of $(s.damping.lambda_tolerance_per_s) s⁻¹.")
        println(io,"Wing span scales both the aerodynamic and structural meshes. The wake cutoff applies to both wakes.\n")
        println(io,"| Family | Comparison | Result |")
        println(io,"|---|---|---|")
        for (family,a,b,result) in pairs
            println(io,"| $family | $a → $b | $result |")
        end
        println(io,"\nSee \`progress.csv\` for signed growth rates, fit-window shifts, runtimes, and case histories.")
        println(io,"Each completed history can be reused when its resolved inputs, solver source and file hash match.")
    end
end

function run_staged_aeroelastic(s,selection)
    s.wall_hours > 0 || throw(ArgumentError("wallHours must be positive"))
    deadline = time() + 3600s.wall_hours
    directory = abspath(s.output_directory)
    mkpath(directory)
    lockPath = joinpath(directory,"running.lock")
    isdir(lockPath) && error("Staged run directory is locked: $lockPath")
    mkdir(lockPath)
    try
        return run_staged_locked(s,selection,directory,deadline)
    finally
        rm(lockPath)
    end
end

function run_staged_locked(s,selection,directory,deadline)
    pilot,screen,followup = staged_plan(s)
    cases = [pilot; screen; followup]
    p = selection.physical
    Convergence.load_model()
    wrapped = ExplicitAeroelasticSettings(MovingBlockStudy(s),explicit_driver_fingerprint())
    Base.invokelatest(Convergence.validate_elastic,wrapped,selection)
    plan = [merge(c,staged_work(c,p,s)) for c in cases]
    write_staged_table(joinpath(directory,"runtime_plan.csv"),plan,keys(first(plan)))
    println("Staged Chang study: $(p.speed_mps) m/s, $(p.rpm) RPM, $(Threads.nthreads()) Julia threads")
    println("Results: $directory")
    println("$(1+length(screen)) pilot/screen cases and $(length(followup)) conditional refinements")
    rows = NamedTuple[]
    result = (; directory,cases,results=rows,selection,nominal=Convergence.controls(pilot))
    s.dry_run && return result

    modelScope = Base.invokelatest(getproperty,Convergence,:Model)
    isdefined(modelScope,:ChangRuntimeOptimizations) ||
        Base.include(modelScope,joinpath(Convergence.EXAMPLE,"src","ChangRuntimeOptimizations.jl"))
    fingerprint = Convergence.SOURCE()
    Convergence.save_toml(joinpath(directory,"study.toml"),(; source_sha256=fingerprint,
        settings=s,selection,threads=Threads.nthreads()))
    previousThreads = BLAS.get_num_threads()
    BLAS.set_num_threads(1)
    try
        Base.invokelatest(execute_staged_plan!,rows,pilot,screen,followup,s,selection,directory,
            fingerprint,deadline)
    finally
        BLAS.set_num_threads(previousThreads)
    end
    return result
end

function execute_staged_plan!(rows,pilot,screen,followup,s,selection,directory,fingerprint,deadline)
    timings = NamedTuple[]
    comparisons = Tuple[]
    function run_case(c)
        Convergence.SOURCE() == fingerprint || error("Solver source changed during the study")
        println("$(c.stage): $(c.label)")
        r = try
            staged_case(c,s,selection,directory,fingerprint,deadline,timings)
        catch exception
            exception isa InterruptException && rethrow()
            status = exception isa StagedDeadline ? "deferred" : "failed"
            reason = exception isa StagedDeadline ? "Wall budget reached during this case" :
                sprint(showerror,exception)
            (; status,reason,lambda_per_s=NaN,frequency_hz=NaN,lambda_window_shift_per_s=NaN,
                elapsed_s=0.0,reused=false,history="",log="")
        end
        row = merge(c,r)
        push!(rows,row)
        write_staged_table(joinpath(directory,"progress.csv"),rows,keys(row))
        println("  $(r.status), growth rate = $(r.lambda_per_s) s^-1 $(r.reason)")
        return row
    end

    staged_report(directory,selection,s,"Running pilot",comparisons)
    baseline = run_case(pilot)
    if baseline.status != "completed"
        staged_report(directory,selection,s,"Pilot $(baseline.status). The remaining cases require an identifiable pilot mode.",comparisons)
        return
    end
    selected = Convergence.controls(pilot)
    lastCases = Dict{Symbol,NamedTuple}()
    compared = Dict{Symbol,Bool}()
    for c in screen
        fine = run_case(c)
        matches = staged_agrees(baseline,fine,s)
        if c.stage == :coarser
            push!(comparisons,(c.family,getproperty(c,c.family),getproperty(baseline,c.family),
                matches ? "within tolerance" : fine.status == "completed" ? "sensitive" : fine.status))
            staged_report(directory,selection,s,"Screening",comparisons)
            continue
        end
        compared[c.family] = matches
        lastCases[c.family] = fine
        if !matches && fine.status == "completed"
            selected = merge(selected,NamedTuple{(c.family,)}((getproperty(c,c.family),)))
        end
        push!(comparisons,(c.family,getproperty(baseline,c.family),getproperty(c,c.family),
            matches ? "within tolerance" : fine.status == "completed" ? "sensitive" : fine.status))
        staged_report(directory,selection,s,"Screening",comparisons)
    end
    # Test affordable sensitive refinements before the combined run. The 1.25°
    # case is deliberately last because shortening the step also adds wake rows.
    for c in sort(followup; by=c -> c.family == :azimuth_deg)
        coarse = lastCases[c.family]
        coarse.status == "completed" && !compared[c.family] || continue
        fine = run_case(c)
        matches = staged_agrees(coarse,fine,s)
        compared[c.family] = matches
        if fine.status == "completed"
            lastCases[c.family] = fine
            if !matches
                selected = merge(selected,NamedTuple{(c.family,)}((getproperty(c,c.family),)))
            end
        end
        push!(comparisons,(c.family,getproperty(coarse,c.family),getproperty(c,c.family),
            matches ? "within tolerance" : fine.status == "completed" ? "sensitive" : fine.status))
    end
    if any(r -> r.status != "completed",values(lastCases))
        staged_report(directory,selection,s,"Some finer comparisons are unresolved. See completed and deferred cases in progress.csv.",comparisons)
        return
    end
    candidate = Convergence.key(selected) == Convergence.key(pilot) ? baseline :
        run_case(merge(selected,(; family=:combined,level=0,label="combined_candidate",stage=:combined)))
    cross = merge(selected,(; azimuth_deg=selected.azimuth_deg / 2,
        wake_revolutions=selected.wake_revolutions + 1))
    verified = false
    if candidate.status == "completed"
        fine = run_case(merge(cross,(; family=:combined,level=1,label="combined_step_wake",stage=:combined)))
        verified = staged_agrees(candidate,fine,s)
        push!(comparisons,(:combined,"$(selected.azimuth_deg)° / $(selected.wake_revolutions) rev",
            "$(cross.azimuth_deg)° / $(cross.wake_revolutions) rev",
            verified ? "within tolerance" : fine.status == "completed" ? "sensitive" : fine.status))
    end
    status = if verified && all(values(compared))
        "Combined refinement agrees within the requested tolerances."
    else
        "Screening complete. Damping convergence is unresolved at the requested tolerances."
    end
    staged_report(directory,selection,s,status,comparisons)
end
