module PazyFlutterSweep

let project = normpath(joinpath(@__DIR__,".."))
    project in LOAD_PATH || pushfirst!(LOAD_PATH,project)
end

using DelimitedFiles
include(joinpath(@__DIR__,"PazySteppedBifurcation.jl"))
import .PazySteppedBifurcation as Stepped
include(joinpath(@__DIR__,"PazyFlutterAnalysis.jl"))
export run_pazy_flutter_sweep

const EQUILIBRIUM_COLUMNS = [:angle_deg,:speed_mps,:variant,:dt_s,:wake_rows,:wake_chords,
    :status,:reason,:growth_per_s,:frequency_hz,:mode_ratio_real,:mode_ratio_imag,
    :trim_bending_m,:trim_twist_deg,:trim_strain_ue,:background_bending_m,:background_twist_deg]
const BOUNDARY_COLUMNS = [:angle_deg,:quantity,:path,:transition,:lower_speed_mps,
    :upper_speed_mps,:width_mps,:status,:estimate_mps]
const LCO_COLUMNS = [:angle_deg,:path,:plateau_index,:direction,:speed_mps,:status,
    :bending_amplitude_m,:twist_amplitude_deg,:frequency_hz,:bending_floor_m,:twist_floor_deg]
const LCO_RESPONSES = (:tip_out_of_plane_m,:tip_twist_deg,:root_strain_ue)

function lco_paths(slice,step)
    lo,hi = slice.lco_limits
    seed = slice.lco_seed
    grid(a,b) = Stepped.TD.bifurcation_speeds(a,b,step)
    return (
        (; name=:seed_down,direction=:downsweep,seeded=true,speeds=grid(lo,seed)),
        (; name=:seed_up,direction=:upsweep,seeded=true,speeds=grid(seed,hi)),
        (; name=:entry_up,direction=:upsweep,seeded=false,speeds=grid(lo,hi)),
        (; name=:entry_down,direction=:downsweep,seeded=false,speeds=grid(lo,hi)))
end

function number_tag(x)
    return replace(replace(string(Float64(x)),r"\.0$" => ""),"." => "p","-" => "m")
end

function validate_settings(slices,o,stage)
    stage in (:pilot,:equilibrium,:lco,:all) || throw(ArgumentError("Invalid flutter sweep stage"))
    !isempty(slices) && allunique(s.angle_deg for s in slices) ||
        throw(ArgumentError("Supply distinct root incidences"))
    for s in slices
        lo,hi = s.lco_limits
        all(isfinite,(s.angle_deg,s.pilot_speed,lo,hi,s.lco_seed)) &&
            0 < lo < s.lco_seed < hi && s.pilot_speed > 0 &&
            length(s.speeds) >= 3 && all(u -> isfinite(u) && u > 0,s.speeds) &&
            all(>(0),diff(collect(s.speeds))) ||
            throw(ArgumentError("Invalid angle, airspeed grid, or LCO seed/extents"))
    end
    all(isfinite,(o.density,o.core_radius,o.wake_shedding_fraction,o.coarse_dt,o.fine_dt,
        o.wake_length_chords,o.wake_factor,o.speed_tolerance,o.lco_step,o.duration,
        o.pulse_start,o.pulse_duration,o.lco_pulse,o.plateau_time)) &&
        o.density > 0 && o.core_radius > 0 && 0 <= o.wake_shedding_fraction <= 1 &&
        all(n -> n isa Integer && n > 0,(o.chordwise_panels,o.spanwise_panels)) &&
        0 < o.fine_dt <= o.coarse_dt/2 && o.wake_length_chords > 0 && o.wake_factor > 1 &&
        o.maximum_refinements isa Integer && o.maximum_refinements >= 0 &&
        0 < o.speed_tolerance <= o.lco_step && length(o.pulse_forces) == 2 &&
        all(isfinite,o.pulse_forces) &&
        0 < o.pulse_forces[1] < o.pulse_forces[2] && o.pulse_start >= 0 &&
        o.pulse_duration >= 8o.coarse_dt && o.lco_pulse > 0 ||
        throw(ArgumentError("Invalid discretization, refinement or pulse controls"))
    g = o.growth_options
    length(g.windows) == 2 && o.pulse_start+o.pulse_duration <= g.windows[1][1] &&
        g.windows[1][2] <= g.windows[2][1] &&
        all(w -> w[2]-w[1] >= 2.5g.block_time,g.windows) &&
        g.block_time >= 8o.coarse_dt && o.duration >= last(last(g.windows)) &&
        0 < g.frequency_band[1] < g.frequency_band[2] < 0.5/o.coarse_dt ||
        throw(ArgumentError("Provide two post-pulse fit windows, each with enough blocks"))
    o.lco_options.tail_time >= 8o.fine_dt && o.lco_options.minimum_cycles >= 2 &&
        o.plateau_time-2o.lco_options.tail_time >= o.pulse_start+o.pulse_duration ||
        throw(ArgumentError("Two LCO tail windows must fit after the pulse"))
    isapprox(round(o.plateau_time/o.fine_dt)*o.fine_dt,o.plateau_time;atol=1e-10,rtol=1e-12) ||
        throw(ArgumentError("The fine timestep must divide plateau_time"))
    all(k -> k in Stepped.MODEL_KEYS,keys(o.model_options)) ||
        throw(ArgumentError("model_options accepts only structural Pazy settings"))
    for (given,defaults) in ((o.dynamic_options,Stepped.DYNAMIC_DEFAULTS),
        (o.equilibrium_options,Stepped.STATIC_DEFAULTS))
        all(k -> k in keys(defaults),keys(given)) || throw(ArgumentError("Unknown solver control"))
    end
end

"""
    run_pazy_flutter_sweep(; slices, stage=:pilot, ...)

Estimate small-disturbance hump-mode growth and follow finite-amplitude LCOs
in separate experiments. `stage=:pilot` runs the timestep and wake diagnostic,
`:equilibrium` runs the fixed-speed sweep, `:lco` runs the continuation paths,
and `:all` runs both sweeps. Each run gets a fresh output directory.
"""
function run_pazy_flutter_sweep(;
    slices=((angle_deg=5.0,pilot_speed=45.0,speeds=(40.,43.,44.,45.,49.,50.,50.5,52.),
        lco_seed=45.0,lco_limits=(35.0,53.0)),),
    stage=:pilot,dry_run=false,output_directory=nothing,
    model_options=(; mesh=:refined,rootTorsionCorr=false),
    density=1.225,chordwise_panels=10,spanwise_panels=20,symmetric_wing=true,
    core_radius=1e-3,wake_shedding_fraction=0.1,
    coarse_dt=0.0004,fine_dt=0.0001,wake_length_chords=6.5,wake_factor=2.0,
    trim_atol=(2e-4,0.02,2.0),trim_rtol=0.02,
    duration=2.25,pulse_start=0.10,pulse_duration=0.04,pulse_forces=(0.025,0.05),
    speed_tolerance=0.5,maximum_refinements=12,
    lco_step=2.0,lco_pulse=0.1,plateau_time=10.0,refine_lco=true,include_entry_paths=true,
    growth_options=(;),lco_options=(;),dynamic_options=(;),equilibrium_options=(;),
    save_histories=true,plot_results=true)
    o = (; model_options,density,chordwise_panels,spanwise_panels,symmetric_wing,
        core_radius,wake_shedding_fraction,coarse_dt,fine_dt,wake_length_chords,wake_factor,
        trim_atol,trim_rtol,duration,pulse_start,pulse_duration,pulse_forces,speed_tolerance,
        maximum_refinements,lco_step,lco_pulse,plateau_time,refine_lco,include_entry_paths,
        growth_options=merge(GROWTH_DEFAULTS,growth_options),
        lco_options=merge(LCO_DEFAULTS,lco_options),dynamic_options,equilibrium_options,
        save_histories,plot_results)
    validate_settings(slices,o,stage)
    if dry_run
        return (; stage,slices,options=o,
            paths=[filter(p -> p.seeded || include_entry_paths,lco_paths(s,lco_step)) for s in slices])
    end
    Stepped.load_backend()
    return Base.invokelatest(run_loaded,slices,o,stage,output_directory)
end

function evaluate_speed(B,angle,speed,variant,dt,wakeChords,chord,o,root)
    tag = "aoa$(number_tag(angle))_V$(number_tag(speed))_$variant"
    directory = mkpath(joinpath(root,tag))
    wakeRows = ceil(Int,wakeChords*chord/(speed*dt))
    common = (; o.model_options...,airspeed=speed,angle_of_attack=deg2rad(angle),
        o.density,o.chordwise_panels,o.spanwise_panels,o.symmetric_wing,o.core_radius,
        o.wake_shedding_fraction,aerodynamic_interaction=true,
        maximum_wake_rows=wakeRows,time_step_chords=dt*speed/chord)
    println("Fixed speed $tag: dt=$dt s, wake=$wakeRows rows")
    flush(stdout)
    row = (; angle_deg=angle,speed_mps=speed,variant,dt_s=dt,wake_rows=wakeRows,
        wake_chords=wakeRows*speed*dt/chord)
    try
        staticControls = merge(Stepped.STATIC_DEFAULTS,o.equilibrium_options)
        # Static settling requires at least two retained-wake ages.
        if !haskey(o.equilibrium_options,:maximum_aerodynamic_time)
            staticControls = merge(staticControls,(;
                maximum_aerodynamic_time=max(staticControls.maximum_aerodynamic_time,3wakeRows*dt)))
        end
        equilibrium = B.run_pazy_wing_uvlm_static(; staticControls...,common...)
        writedlm(joinpath(directory,"initialization.csv"),
            vcat(permutedims(equilibrium.history_columns),equilibrium.history),',')
        equilibrium.converged || error("Static initialization: $(equilibrium.termination_reason)")
        dynamicControls = merge(Stepped.DYNAMIC_DEFAULTS,o.dynamic_options)
        responses = Matrix{Float64}[]
        time = Float64[]
        trim = (NaN,NaN,NaN)
        # Diagnose a failed perturbation before spending a full run on its control.
        for force in (o.pulse_forces...,0.0)
            println("$tag: integrating pulse coefficient $force N for $(o.duration) s")
            flush(stdout)
            solution = B.run_pazy_wing_uvlm(; dynamicControls...,common...,
                initial_equilibrium=equilibrium,duration=o.duration,sideslip=0.0,
                initial_airspeed_fraction=1.0,airspeed_ramp_duration=0.0,
                settling_time=o.pulse_start,perturbation_amplitude=force,
                perturbation_duration=o.pulse_duration,save_uvlm_history=false,uvlm_save_frequency=1,
                animation_frames=2,save_animation_history=false,save_structural_history=false,
                save_aerodynamic_load_history=false)
            if isempty(time)
                time = solution.time
                trim = (first(solution.tip_out_of_plane),first(solution.tip_twist_degrees),
                    1e6 * (first(solution.root_axial_strain)+0.001125*first(solution.root_curvature_2)))
            end
            response = hcat(solution.tip_out_of_plane,solution.tip_twist_degrees)
            all(isfinite,response) || error("Nonfinite dynamic response")
            push!(responses,response)
            if o.save_histories
                writedlm(joinpath(directory,"pulse_$(number_tag(force)).csv"),
                    vcat(["time_s" "tip_bending_m" "tip_twist_deg"],hcat(time,response)),',')
            end
        end
        fit = classify_equilibrium(time,responses[3],responses[1],responses[2],
            o.pulse_forces,o.growth_options)
        fitRows = NamedTuple[]
        for j in axes(fit.fits,2), i in axes(fit.fits,1)
            f = fit.fits[i,j]
            isnothing(f) && continue
            push!(fitRows,(; pulse_N=o.pulse_forces[i],start_s=o.growth_options.windows[j][1],
                end_s=o.growth_options.windows[j][2],frequency_hz=f.frequency,
                bending_frequency_hz=f.bending_frequency,bending_growth=f.growth[1],twist_growth=f.growth[2],
                bending_r_squared=f.r_squared[1],twist_r_squared=f.r_squared[2],
                bending_amplitude=f.amplitudes[1],twist_amplitude=f.amplitudes[2],
                bending_background=f.background[1],twist_background=f.background[2],
                same_mode=f.same_mode,above_noise=f.above_noise,observable_growth=f.observable_growth))
        end
        isempty(fitRows) || Stepped.write_rows(joinpath(directory,"fits.csv"),fitRows,collect(keys(first(fitRows))))
        bg = ntuple(c -> maximum((f.background[c] for f in fit.fits if !isnothing(f));init=0.0),2)
        return (; row...,fit.status,fit.reason,growth_per_s=fit.rate,frequency_hz=fit.frequency,
            mode_ratio_real=real(fit.signature),mode_ratio_imag=imag(fit.signature),
            trim_bending_m=trim[1],trim_twist_deg=trim[2],trim_strain_ue=trim[3],
            background_bending_m=bg[1],background_twist_deg=bg[2])
    catch err
        err isa InterruptException && rethrow()
        open(joinpath(directory,"failure.txt"),"w") do io
            showerror(io,err,catch_backtrace())
        end
        @warn "Fixed-speed case failed" tag exception=err
        return (; row...,status=:failed,reason="see case failure.txt",growth_per_s=NaN,
            frequency_hz=NaN,mode_ratio_real=NaN,mode_ratio_imag=NaN,
            trim_bending_m=NaN,trim_twist_deg=NaN,trim_strain_ue=NaN,
            background_bending_m=NaN,background_twist_deg=NaN)
    end
end

function same_mode(a,b,o)
    return abs(a.frequency_hz-b.frequency_hz) <= o.frequency_tolerance &&
        mode_close(complex(a.mode_ratio_real,a.mode_ratio_imag),
            complex(b.mode_ratio_real,b.mode_ratio_imag),o.mode_rtol)
end

function numerical_agreement(a,b,o)
    a.status == b.status && a.status in (:decay,:growth) &&
        same_mode(a,b,o.growth_options) && rate_close(a.growth_per_s,b.growth_per_s,o.growth_options) ||
        return false
    fields = (:trim_bending_m,:trim_twist_deg,:trim_strain_ue)
    return all(i -> begin
        x,y = getproperty(a,fields[i]),getproperty(b,fields[i])
        abs(x-y) <= max(o.trim_atol[i],o.trim_rtol*max(abs(x),abs(y)))
    end,eachindex(fields))
end

function equilibrium_brackets(cases,o)
    rows = sort!(collect(values(cases));by=r -> r.speed_mps)
    return [(; lower=rows[b.first],upper=rows[b.last],b.has_gap,
             same_mode=same_mode(rows[b.first],rows[b.last],o.growth_options),
             speeds=Float64[r.speed_mps for r in rows[b.first:b.last]])
            for b in observed_transitions(rows,(:decay,:growth))]
end

function refine_equilibrium!(cases,evaluate,o)
    for _ in 1:o.maximum_refinements
        brackets = filter(b -> last(b.speeds)-first(b.speeds) > o.speed_tolerance,
            equilibrium_brackets(cases,o))
        isempty(brackets) && break
        gaps = [(a,b) for bracket in brackets for (a,b) in zip(bracket.speeds[1:end-1],bracket.speeds[2:end])]
        i = argmax([b-a for (a,b) in gaps])
        a,b = gaps[i]
        b-a >= o.speed_tolerance/2 || break
        speed = (a+b)/2
        cases[speed] = evaluate(speed,:fine,o.fine_dt,o.wake_length_chords)
    end
end

function save_equilibrium_boundaries!(boundaries,cases,evaluate,o,root)
    refine_equilibrium!(cases,evaluate,o)
    for bracket in equilibrium_brackets(cases,o)
        a,b = bracket.lower,bracket.upper
        width = b.speed_mps-a.speed_mps
        status = :mode_unresolved
        if bracket.same_mode
            checkA = evaluate(a.speed_mps,:edge_dt_half,o.fine_dt/2,o.wake_length_chords)
            checkB = evaluate(b.speed_mps,:edge_dt_half,o.fine_dt/2,o.wake_length_chords)
            numerical = numerical_agreement(a,checkA,o) && numerical_agreement(b,checkB,o)
            if numerical
                wakeA = evaluate(a.speed_mps,:long_wake,o.fine_dt,o.wake_factor*o.wake_length_chords)
                wakeB = evaluate(b.speed_mps,:long_wake,o.fine_dt,o.wake_factor*o.wake_length_chords)
                numerical = numerical_agreement(a,wakeA,o) && numerical_agreement(b,wakeB,o)
            end
            status = !numerical ? :resolution_unverified : bracket.has_gap ? :unresolved_between :
                width > o.speed_tolerance ? :coarse_bracket : :resolved
        end
        estimate = status == :resolved ? a.speed_mps -
            a.growth_per_s*width/(b.growth_per_s-a.growth_per_s) : NaN
        push!(boundaries,(; a.angle_deg,quantity=:equilibrium,path=:increasing_speed,
            transition=a.status == :decay ? :onset : :restabilization,
            lower_speed_mps=a.speed_mps,upper_speed_mps=b.speed_mps,width_mps=width,status,
            estimate_mps=estimate))
        Stepped.write_rows(joinpath(root,"boundaries.csv"),boundaries,BOUNDARY_COLUMNS)
    end
end

struct InitialLCOStateUnverified <: Exception
    status::Symbol
end
Base.showerror(io::IO,e::InitialLCOStateUnverified) =
    print(io,"Initial plateau was classified as ",e.status)

function lco_brackets(rows,slice,path,tolerance)
    return map(observed_transitions(rows,(:steady,:lco))) do b
        a,c = rows[b.first],rows[b.last]
        lo,hi = minmax(a.speed_mps,c.speed_mps)
        width = hi-lo
        (; angle_deg=slice.angle_deg,quantity=:lco,path=path.name,
            transition=a.status == :steady ? :entry : :extinction,
            lower_speed_mps=lo,upper_speed_mps=hi,width_mps=width,
            status=b.has_gap ? :unresolved_between : width > tolerance ? :coarse_bracket : :resolved,
            estimate_mps=NaN)
    end
end

function run_lco_path(slice,path,phase,o,floor,chord,root)
    name = "aoa$(number_tag(slice.angle_deg))_$(path.name)_$phase"
    directory = joinpath(root,name)
    rows = NamedTuple[]
    function observe(plateau,time,values)
        response = hcat(values.tip_out_of_plane_m,values.tip_twist_deg)
        measured = classify_lco(time,response,floor,o.lco_options)
        push!(rows,(; angle_deg=slice.angle_deg,path=path.name,
            plateau_index=plateau.index,direction=plateau.direction,speed_mps=plateau.airspeed_mps,
            measured.status,bending_amplitude_m=measured.amplitudes[1],
            twist_amplitude_deg=measured.amplitudes[2],frequency_hz=measured.frequency,
            bending_floor_m=floor[1],twist_floor_deg=floor[2]))
        Stepped.write_rows(joinpath(directory,"lco_classification.csv"),rows,LCO_COLUMNS)
        expected = path.seeded ? :lco : :steady
        plateau.index == 1 && measured.status != expected &&
            throw(InitialLCOStateUnverified(measured.status))
    end
    try
        # Retain at least the requested nominal wake length over all paths.
        wakeRows = ceil(Int,o.wake_length_chords*chord/(first(slice.lco_limits)*o.fine_dt))
        staticControls = !haskey(o.equilibrium_options,:maximum_aerodynamic_time) ? merge(o.equilibrium_options,(;
            maximum_aerodynamic_time=max(Stepped.STATIC_DEFAULTS.maximum_aerodynamic_time,
                3wakeRows*o.fine_dt))) : o.equilibrium_options
        Stepped.run_pazy_stepped_bifurcation(; airspeeds=path.speeds,sweep_order=path.direction,
            angle_of_attack=deg2rad(slice.angle_deg),o.density,o.chordwise_panels,o.spanwise_panels,
            o.symmetric_wing,o.core_radius,o.wake_shedding_fraction,o.model_options,
            fixed_time_step=o.fine_dt,time_step_chords=max(0.25,o.fine_dt*last(path.speeds)/chord),
            maximum_wake_rows=wakeRows,plateau_time=o.plateau_time,
            transient_time=o.plateau_time-2o.lco_options.tail_time,
            pulse_amplitude=o.lco_pulse,pulse_start=o.pulse_start,pulse_duration=o.pulse_duration,
            responses=LCO_RESPONSES,save_time_histories=o.save_histories,
            on_unconverged=:continue,plot_results=false,
            o.dynamic_options,equilibrium_options=staticControls,plateau_callback=observe,
            output_directory=directory)
    catch err
        err isa InterruptException && rethrow()
        expected = path.seeded ? :lco : :steady
        @warn "LCO continuation stopped" name expected exception=err
        mkpath(directory)
        if !isfile(joinpath(directory,"failure.txt"))
            open(joinpath(directory,"failure.txt"),"w") do io
                showerror(io,err,catch_backtrace())
            end
        end
    end
    expected = path.seeded ? :lco : :steady
    if isempty(rows) || first(rows).status != expected
        return NamedTuple[]
    end
    return lco_brackets(rows,slice,path,o.speed_tolerance)
end

function run_lco_paths!(boundaries,slice,o,floor,chord,root)
    for path in lco_paths(slice,o.lco_step)
        path.seeded || o.include_entry_paths || continue
        coarse = run_lco_path(slice,path,:coarse,o,floor,chord,root)
        append!(boundaries,coarse)
        Stepped.write_rows(joinpath(root,"boundaries.csv"),boundaries,BOUNDARY_COLUMNS)
        candidates = filter(b -> b.width_mps > o.speed_tolerance,coarse)
        if o.refine_lco && !isempty(candidates)
            speedGrid = copy(path.speeds)
            for b in candidates
                n = ceil(Int,b.width_mps/o.speed_tolerance)
                append!(speedGrid,range(b.lower_speed_mps,b.upper_speed_mps;length=n+1))
            end
            finePath = merge(path,(; speeds=sort!(unique!(speedGrid))))
            # Keep coarse observations as a separate path to expose rate sensitivity.
            fine = run_lco_path(slice,finePath,:fine,o,floor,chord,root)
            append!(boundaries,[merge(b,(; path=Symbol(b.path,:_refined))) for b in fine])
            Stepped.write_rows(joinpath(root,"boundaries.csv"),boundaries,BOUNDARY_COLUMNS)
        end
    end
end

function run_loaded(slices,o,stage,outputDirectory)
    B = Stepped.PazyBifurcation
    _,chord = Stepped.pazy_span_and_chord(B,o.model_options)
    if isnothing(outputDirectory)
        parent = mkpath(joinpath(@__DIR__,"output","flutter_sweep"))
        root = mktempdir(parent;prefix="run_",cleanup=false)
    else
        root = abspath(outputDirectory)
        # mkdir reserves the destination exclusively.
        mkpath(dirname(root))
        mkdir(root)
    end
    println("Flutter sweep outputs: $root")
    open(joinpath(root,"settings.txt"),"w") do io
        println(io,"Julia version: ",VERSION,"; threads: ",Threads.nthreads(),
            "; BLAS threads: ",BLAS.get_num_threads())
        println(io,"Stage: $stage")
        println(io,"Slices: $slices")
        println(io,"Settings: $o")
        println(io,"Trim controls: ",merge(Stepped.STATIC_DEFAULTS,o.equilibrium_options))
        println(io,"Dynamic controls: ",merge(Stepped.DYNAMIC_DEFAULTS,o.dynamic_options))
    end
    equilibriumRows,boundaries,pilots = NamedTuple[],NamedTuple[],NamedTuple[]
    for (name,rows,columns) in (("equilibrium_cases.csv",equilibriumRows,EQUILIBRIUM_COLUMNS),
        ("boundaries.csv",boundaries,BOUNDARY_COLUMNS))
        Stepped.write_rows(joinpath(root,name),rows,columns)
    end
    for (i,slice) in enumerate(slices)
        println("Incidence $i/$(length(slices)): $(slice.angle_deg) deg ($stage)")
        flush(stdout)
        cache = Dict{Tuple{Float64,Symbol},NamedTuple}()
        function evaluate(speed,variant,dt,wake)
            return get!(cache,(speed,variant)) do
                row = evaluate_speed(B,slice.angle_deg,speed,variant,dt,wake,chord,o,root)
                push!(equilibriumRows,row)
                Stepped.write_rows(joinpath(root,"equilibrium_cases.csv"),equilibriumRows,EQUILIBRIUM_COLUMNS)
                println("U=$speed, $variant: $(row.status), growth=$(row.growth_per_s) s^-1")
                flush(stdout)
                row
            end
        end
        if stage == :pilot
            speed = slice.pilot_speed
            evaluate(speed,:coarse,o.coarse_dt,o.wake_length_chords)
            medium = evaluate(speed,:medium,2o.fine_dt,o.wake_length_chords)
            fine = evaluate(speed,:fine,o.fine_dt,o.wake_length_chords)
            wake = evaluate(speed,:long_wake,o.fine_dt,o.wake_factor*o.wake_length_chords)
            dtPassed = numerical_agreement(medium,fine,o)
            wakePassed = numerical_agreement(fine,wake,o)
            push!(pilots,(; slice.angle_deg,speed_mps=speed,dt_passed=dtPassed,wake_passed=wakePassed,
                status=dtPassed && wakePassed ? :passed : :unverified))
            Stepped.write_rows(joinpath(root,"pilot.csv"),pilots,collect(keys(first(pilots))))
        end
        if stage in (:equilibrium,:all)
            cases = Dict{Float64,NamedTuple}()
            for speed in slice.speeds
                cases[speed] = evaluate(speed,:fine,o.fine_dt,o.wake_length_chords)
            end
            save_equilibrium_boundaries!(boundaries,cases,evaluate,o,root)
        end
        if stage in (:lco,:all)
            floor = o.lco_options.amplitude_floor
            seed = get(cache,(slice.lco_seed,:fine),nothing)
            if !isnothing(seed) && all(isfinite,(seed.background_bending_m,seed.background_twist_deg))
                bg = (seed.background_bending_m,seed.background_twist_deg)
                floor = max.(floor,o.lco_options.noise_multiplier .* bg)
            end
            run_lco_paths!(boundaries,slice,o,floor,chord,root)
        end
    end
    plotted = filter(r -> r.status in (:resolved,:coarse_bracket,:unresolved_between),boundaries)
    if o.plot_results && !isempty(plotted)
        @eval import Plots
        Base.invokelatest(write_boundary_plot,plotted,root)
    end
    return (; output_directory=root,pilots,equilibrium_cases=equilibriumRows,boundaries)
end

function write_boundary_plot(rows,root)
    p = Plots.plot(; xlabel="Airspeed [m/s]",ylabel="Root incidence [deg]",
        legend=:outerright,size=(1000,560))
    groups = unique([(r.quantity,r.path,r.transition) for r in rows])
    for (i,group) in enumerate(groups)
        selected = filter(r -> (r.quantity,r.path,r.transition) == group,rows)
        label = join(replace.(String.(group),"_"=>" ")," · ")
        for (j,row) in enumerate(selected)
            x = [row.lower_speed_mps,row.upper_speed_mps]
            y = fill(row.angle_deg,2)
            Plots.plot!(p,x,y; color=i,marker=:vline,
                linestyle=row.status == :resolved ? :solid : :dot,label=j == 1 ? label : false)
            if isfinite(row.estimate_mps)
                Plots.scatter!(p,[row.estimate_mps],[row.angle_deg]; color=i,marker=:circle,label=false)
            end
        end
    end
    Plots.savefig(p,joinpath(root,"pazy_flutter_boundaries.png"))
    return nothing
end

end # module
