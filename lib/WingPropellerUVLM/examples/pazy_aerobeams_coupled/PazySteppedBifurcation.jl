module PazySteppedBifurcation

using DelimitedFiles, LinearAlgebra, Dates
include(joinpath(@__DIR__, "..", "bifurcation", "TimeDomainBifurcation.jl"))
import .TimeDomainBifurcation as TD
export run_pazy_stepped_bifurcation, stepped_plan, plateau_at_step, quadratic_extrema

const DEFAULT_RESPONSES = (:tip_out_of_plane_percent, :root_strain_ue,
    :tip_twist_deg, :tip_geometric_incidence_deg)
const DEFINITIONS = (
    tip_out_of_plane_m=(unit="m", label="Tip out-of-plane displacement", atol=1e-6),
    tip_out_of_plane_percent=(unit="% semispan", label="Tip out-of-plane displacement", atol=1e-4),
    tip_in_plane_m=(unit="m", label="Tip in-plane displacement (+downstream)", atol=1e-6),
    tip_twist_deg=(unit="deg", label="Tip twist", atol=1e-3),
    root_strain_ue=(unit="microstrain", label="Near-root longitudinal strain", atol=0.1),
    tip_geometric_incidence_deg=(unit="deg", label="Root incidence + elastic tip twist", atol=1e-3))
const DYNAMIC_DEFAULTS = (; newton_maximum_iterations=20, newton_absolute_tolerance=1e-5,
    newton_relative_tolerance=1e-5, newton_display_iterations=false,
    newton_always_update_jacobian=true, coupling_scheme=:strong,
    coupling_maximum_iterations=20, coupling_relaxation=0.5,
    coupling_geometry_tolerance=1e-4, coupling_load_tolerance=1e-4,
    coupling_display_iterations=false, progress_frequency=1000)
const STATIC_DEFAULTS = (; maximum_aerodynamic_time=0.20, minimum_aerodynamic_time=0.05,
    aerodynamic_window_time=0.01, aerodynamic_force_tolerance=1e-4,
    aerodynamic_circulation_tolerance=1e-4, maximum_static_iterations=100,
    static_relaxation=0.3, geometry_tolerance=1e-5, load_tolerance=1e-4,
    consecutive_equilibrium_iterations=2, newton_maximum_iterations=50,
    newton_absolute_tolerance=1e-8, newton_relative_tolerance=1e-8,
    newton_display_iterations=false, progress_frequency=1000)

"""Plan two sequences of fixed-speed plateaus. The reversal speed is held twice,
once in each direction, matching the supplied propeller procedure. `chord` is
supplied by the physical model; no Pazy dimensions are hard-coded here.
"""
function stepped_plan(; U_min=40.0, U_max=50.0, dU=1.0, plateau_time=10.0,
    transient_time=8.0, sweep_order=:up_down, chord, time_step_chords=0.25,
    fixed_time_step=nothing)
    grid = TD.bifurcation_speeds(U_min, U_max, dU)
    sweep_order in (:up_down, :down_up, :upsweep, :downsweep) ||
        throw(ArgumentError("sweep_order must be :up_down, :down_up, :upsweep or :downsweep"))
    all(isfinite, (plateau_time, transient_time, chord, time_step_chords)) &&
        plateau_time > 0 && 0 <= transient_time < plateau_time && chord > 0 && time_step_chords > 0 ||
        throw(ArgumentError("Require positive duration/chord/timestep and 0 <= transient_time < plateau_time"))
    nominal_dt = time_step_chords*chord/maximum(grid)
    steps = ceil(Int, plateau_time/nominal_dt)
    dt = plateau_time/steps
    if !isnothing(fixed_time_step)
        dt = Float64(fixed_time_step)
        isfinite(dt) && 0 < dt <= nominal_dt*(1+1e-12) ||
            throw(ArgumentError("fixed_time_step must be positive and <= time_step_chords*chord/U_max"))
        steps = round(Int, plateau_time/dt)
        isapprox(steps*dt, plateau_time; atol=1e-10, rtol=1e-12) ||
            throw(ArgumentError("fixed_time_step must divide plateau_time"))
    end
    # Exclude the shared boundary sample: it was integrated at the preceding speed.
    observation_start = max(transient_time, dt)
    local_time = collect(0:steps).*dt
    count(t -> t >= observation_start, local_time) >= 3 ||
        throw(ArgumentError("At least three samples must remain after transient removal"))
    directions = sweep_order == :up_down ? (:upsweep,:downsweep) :
        sweep_order == :down_up ? (:downsweep,:upsweep) : (sweep_order,)
    plateaus = NamedTuple[]
    for direction in directions, speed in (direction == :upsweep ? grid : reverse(grid))
        k = length(plateaus)+1
        push!(plateaus, (; index=k, direction, airspeed_mps=speed,
            start_step=(k-1)*steps, end_step=k*steps,
            start_time_s=(k-1)*steps*dt, end_time_s=k*steps*dt))
    end
    return (; grid, plateaus, dt, steps_per_plateau=steps, local_time,
        plateau_time=steps*dt, transient_time, observation_start, sweep_order,
        total_steps=length(plateaus)*steps, total_time=length(plateaus)*steps*dt,
        nominal_dt, chord)
end

"""The step ending at a plateau boundary belongs to the preceding plateau."""
function plateau_at_step(plan, step::Integer)
    0 <= step <= plan.total_steps || throw(ArgumentError("Step outside the planned integration"))
    return clamp(cld(step,plan.steps_per_plateau), 1, length(plan.plateaus))
end

"""Refine observed extrema by fitting three neighboring response samples and
solving dq/dt=0 for the local quadratic. This approximates the propeller's
derivative-zero procedure for scalar outputs without an ODE dense interpolant.
The same observable is differentiated (including strain/twist), not a surrogate
nodal velocity. Fits never cross the observation window or a speed change.
Flat peaks retain the sampled convention used by TimeDomainBifurcation.
"""
function quadratic_extrema(time, q; transient_time)
    raw = TD.response_local_extrema(time,q; transient_time)
    start = searchsortedfirst(time,transient_time)
    return map(raw) do peak
        i = searchsortedfirst(time,peak.time_s)
        hl, hr = time[i]-time[i-1], time[i+1]-time[i]
        if i <= start || q[i] == q[i-1] || q[i] == q[i+1]
            return peak
        end
        dl, dr = (q[i]-q[i-1])/hl, (q[i+1]-q[i])/hr
        a = (dr-dl)/(hl+hr)
        b = dl+a*hl
        iszero(a) && return peak
        offset = -b/(2a)
        -hl <= offset <= hr || return peak
        (peak.kind == "maximum" ? a < 0 : a > 0) || return peak
        return (; kind=peak.kind, time_s=time[i]+offset, value=q[i]+b*offset+a*offset^2)
    end
end

const SUMMARY_COLUMNS = vcat(Symbol.(TD.SUMMARY_COLUMNS),
    [:plateau_index,:direction,:global_window_start_s,:global_window_end_s,
     :sample_minimum,:sample_maximum,:extrema_method])
const PEAK_COLUMNS = [:airspeed_mps,:response,:unit,:kind,:time_s,:value,
    :plateau_index,:direction,:global_time_s,:extrema_method]
const BOUNDARY_COLUMNS = [:plateau_index,:direction,:airspeed_mps,:global_time_s,
    :accepted_step,:aerodynamic_step,:aerodynamic_time_s,:retained_wake_rows,
    :circulation_norm,:structural_state_norm,:maximum_section_velocity_norm,
    :structural_problem_id,:aerodynamic_problem_id]

function write_rows(path, rows, columns)
    open(path,"w") do io
        writedlm(io,permutedims(String.(columns)),',')
        for row in rows
            writedlm(io,permutedims(Any[getproperty(row,k) for k in columns]),',')
        end
    end
end

function analyze_plateau(plan, plateau, values, initial_static;
    responses, angle_of_attack, extrema_method, check_convergence, rtol, convergence_cycles,
    response_definitions)
    rows, peaks = NamedTuple[], NamedTuple[]
    for response in responses
        q = values[response]
        definition = getproperty(response_definitions,response)
        stats = TD.response_statistics(plan.local_time,q; transient_time=plan.observation_start,
            check_convergence, atol=definition.atol, rtol, convergence_cycles)
        local_peaks = extrema_method == :quadratic ?
            quadratic_extrema(plan.local_time,q; transient_time=plan.observation_start) :
            TD.response_local_extrema(plan.local_time,q; transient_time=plan.observation_start)
        low, high = stats.minimum, stats.maximum
        for peak in local_peaks
            low, high = min(low,peak.value), max(high,peak.value)
            push!(peaks, (; airspeed_mps=plateau.airspeed_mps,response=String(response),
                unit=definition.unit,peak...,plateau_index=plateau.index,
                direction=String(plateau.direction),global_time_s=plateau.start_time_s+peak.time_s,
                extrema_method=String(extrema_method)))
        end
        # Only the initially solved static equilibrium is known. No response mean
        # is substituted for an independently calculated static reference.
        static = plateau.airspeed_mps == first(plan.plateaus).airspeed_mps ?
            getproperty(initial_static,response) : NaN
        push!(rows,(; airspeed_mps=plateau.airspeed_mps,response=String(response),unit=definition.unit,
            stats...,minimum=low,maximum=high,amplitude=(high-low)/2,
            static_equilibrium=static,label=definition.label,plateau_index=plateau.index,
            direction=String(plateau.direction),
            global_window_start_s=plateau.start_time_s+stats.window_start_s,
            global_window_end_s=plateau.start_time_s+stats.window_end_s,
            sample_minimum=stats.minimum,sample_maximum=stats.maximum,
            extrema_method=String(extrema_method)))
    end
    return (; rows,peaks)
end

function load_backend()
    isdefined(@__MODULE__,:PazyBifurcation) ||
        Base.include(@__MODULE__,joinpath(@__DIR__,"PazyBifurcation.jl"))
end

"""
    run_pazy_stepped_bifurcation(; U_min=40., U_max=50., dU=1., ...)

One static initialization followed by ONE coupled integration with fixed-speed
plateaus, preserving beam motion, circulation and wake through all speed changes.
Default :up_down includes a separate hold at the reversal speed in each direction.
A single initial tip pulse seeds the response; later plateaus have no new pulse.
`plateau_time=10`, `transient_time=8` selects the final two seconds of each hold.
These are observation settings, not a guarantee of settled LCOs. Increase them
when convergence remains unverified; `on_unconverged=:error` stops AFTER saving
the completed plateau. No adaptive extension or unstable-orbit tracking is done.

`:quadratic` refines local peaks using zeros of a fitted response derivative;
`:sampled` keeps raw sampled extrema. Statistics retain sampled extrema in
sample_minimum/sample_maximum; convergence is assessed on the raw samples.

`dry_run=true` validates controls and returns the schedule without constructing
models or outputs. Including this file runs no simulations. Existing files and
nonempty output directories are never overwritten. Detailed solver histories
are disabled, but the existing solver still retains scalar arrays for the whole
integration. Use `keep_solver_result=true` only for diagnostics.
"""
function run_pazy_stepped_bifurcation(; enabled=true, dry_run=false,
    U_min=40.0,U_max=50.0,dU=1.0,plateau_time=10.0,transient_time=8.0,sweep_order=:up_down,
    angle_of_attack=deg2rad(5.0),density=1.225,chordwise_panels=10,spanwise_panels=20,
    time_step_chords=0.25,fixed_time_step=nothing,wake_length_chords=5.0,maximum_wake_rows=nothing,
    core_radius=1e-3,wake_shedding_fraction=0.1,aerodynamic_interaction=true,
    symmetric_wing=true,pulse_amplitude=0.1,pulse_start=0.1,pulse_duration=0.04,
    gauge_y=0.0,gauge_z=2.25e-3/2,responses=DEFAULT_RESPONSES,
    extrema_method=:quadratic,check_convergence=true,rtol=0.02,convergence_cycles=3,
    response_definitions=DEFINITIONS,on_unconverged=:warn,
    save_time_histories=false,save_boundary_states=false,keep_solver_result=false,
    plot_results=true,plot_options=(;),dynamic_options=(;),equilibrium_options=(;),
    output_directory=nothing)
    enabled || return nothing
    all(isfinite,(angle_of_attack,density,time_step_chords,wake_length_chords,
        core_radius,wake_shedding_fraction,pulse_amplitude,pulse_start,pulse_duration,gauge_y,gauge_z,rtol)) &&
        density > 0 && wake_length_chords > 0 && rtol >= 0 || throw(ArgumentError("Invalid finite model controls"))
    all(n -> n isa Integer && n > 0,(chordwise_panels,spanwise_panels)) || throw(ArgumentError("Invalid panel counts"))
    isnothing(maximum_wake_rows) || (maximum_wake_rows isa Integer && maximum_wake_rows > 0) ||
        throw(ArgumentError("maximum_wake_rows must be a positive integer or nothing"))
    extrema_method in (:quadratic,:sampled) || throw(ArgumentError("Use :quadratic or :sampled extrema"))
    on_unconverged in (:warn,:error,:continue) || throw(ArgumentError("Use :warn, :error or :continue for on_unconverged"))
    convergence_cycles isa Integer && convergence_cycles >= 2 || throw(ArgumentError("convergence_cycles must be >= 2"))
    selected = Tuple(Symbol.(responses))
    !isempty(selected) && allunique(selected) && all(r -> r in keys(DEFINITIONS),selected) ||
        throw(ArgumentError("Select distinct available Pazy responses"))
    !plot_results || all(r -> r in selected,(:tip_out_of_plane_percent,:root_strain_ue,:tip_twist_deg)) ||
        throw(ArgumentError("Four-panel plots require tip_out_of_plane_percent, root_strain_ue and tip_twist_deg; otherwise set plot_results=false"))
    for response in selected
        hasproperty(response_definitions,response) || throw(ArgumentError("Missing response definition: $response"))
        d = getproperty(response_definitions,response)
        isfinite(d.atol) && d.atol >= 0 || throw(ArgumentError("Invalid absolute convergence tolerance"))
    end
    for (given,defaults) in ((dynamic_options,DYNAMIC_DEFAULTS),(equilibrium_options,STATIC_DEFAULTS))
        all(k -> k in keys(defaults),keys(given)) || throw(ArgumentError("Only numerical solver controls may be overridden"))
    end
    pulse_start >= 0 && pulse_duration > 0 || throw(ArgumentError("Invalid pulse timing"))
    iszero(pulse_amplitude) || pulse_start+pulse_duration <= transient_time ||
        throw(ArgumentError("The initial pulse must finish before the first observation window"))
    if !isnothing(output_directory)
        path = abspath(output_directory)
        ispath(path) && (!isdir(path) || !isempty(readdir(path))) &&
            throw(ArgumentError("Output directory must be new or empty: $path"))
    end
    load_backend()
    return Base.invokelatest(run_loaded; U_min,U_max,dU,plateau_time,transient_time,sweep_order,
        angle_of_attack,density,chordwise_panels,spanwise_panels,time_step_chords,fixed_time_step,
        wake_length_chords,maximum_wake_rows,core_radius,wake_shedding_fraction,aerodynamic_interaction,
        symmetric_wing,pulse_amplitude,pulse_start,pulse_duration,gauge_y,gauge_z,responses=selected,
        extrema_method,check_convergence,rtol,convergence_cycles,response_definitions,on_unconverged,
        save_time_histories,save_boundary_states,keep_solver_result,plot_results,plot_options,
        dynamic_options=merge(DYNAMIC_DEFAULTS,dynamic_options),
        equilibrium_options=merge(STATIC_DEFAULTS,equilibrium_options),output_directory,dry_run)
end

function run_loaded(; U_min,U_max,dU,plateau_time,transient_time,sweep_order,
    angle_of_attack,density,chordwise_panels,spanwise_panels,time_step_chords,fixed_time_step,
    wake_length_chords,maximum_wake_rows,core_radius,wake_shedding_fraction,aerodynamic_interaction,
    symmetric_wing,pulse_amplitude,pulse_start,pulse_duration,gauge_y,gauge_z,responses,
    extrema_method,check_convergence,rtol,convergence_cycles,response_definitions,on_unconverged,
    save_time_histories,save_boundary_states,keep_solver_result,plot_results,plot_options,
    dynamic_options,equilibrium_options,output_directory,dry_run)
    B = PazyBifurcation
    _,span,chord,_ = B.AeroBeams.geometrical_properties_Pazy()
    plan = stepped_plan(; U_min,U_max,dU,plateau_time,transient_time,sweep_order,
        chord,time_step_chords,fixed_time_step)
    !iszero(pulse_amplitude) && pulse_duration < 8plan.dt*(1-1e-12) &&
        throw(ArgumentError("Resolve the initial pulse with at least eight timesteps"))
    wake_rows = isnothing(maximum_wake_rows) ? ceil(Int,wake_length_chords*chord/(U_min*plan.dt)) : maximum_wake_rows
    B.UVLM.create_UVLMSolver(;coreRadius=core_radius,wakeSheddingFraction=wake_shedding_fraction,
        interaction=aerodynamic_interaction)
    dry_run && return (;plan,maximum_wake_rows=wake_rows)
    directory = if isnothing(output_directory)
        root=joinpath(@__DIR__,"output","bifurcation_stepped")
        mkpath(root)
        mktempdir(root;prefix="run_",cleanup=false)
    else
        mkpath(abspath(output_directory))
        abspath(output_directory)
    end
    options=(;angle_of_attack,density,chordwise_panels,spanwise_panels,symmetric_wing,
        core_radius,wake_shedding_fraction,aerodynamic_interaction,pulse_amplitude,pulse_start,pulse_duration,
        gauge_y,gauge_z,responses,extrema_method,check_convergence,rtol,convergence_cycles,
        response_definitions,on_unconverged,save_time_histories,dynamic_options,equilibrium_options)
    for direction in unique(p.direction for p in plan.plateaus)
        mkpath(joinpath(directory,String(direction)))
    end
    for dest in vcat([directory],[joinpath(directory,String(d)) for d in unique(p.direction for p in plan.plateaus)])
        open(joinpath(dest,"settings.txt"),"w") do io
            println(io,"Stepped fixed-speed sweep; one coupled integration; no state resets.")
            println(io,"Model/analysis controls: ",options)
            println(io,"Plateau speeds [m/s]: ",[p.airspeed_mps for p in plan.plateaus])
            println(io,"dt=$(plan.dt); steps per plateau=$(plan.steps_per_plateau); total time=$(plan.total_time)")
            println(io,"Wake capacity [rows]: $wake_rows; transient_time=$transient_time")
            println(io,"Peak speed is held once per direction. Only one initial pulse is applied.")
            println(io,"Convergence checks diagnose late-window stationarity, not orbital stability.")
            println(io,"static_equilibrium is NaN except at the initially solved airspeed.")
        end
    end
    schedule_rows=[(;plateau_index=p.index,direction=String(p.direction),p.airspeed_mps,
        p.start_step,p.end_step,p.start_time_s,p.end_time_s) for p in plan.plateaus]
    write_rows(joinpath(directory,"plateau_schedule.csv"),schedule_rows,collect(keys(first(schedule_rows))))
    rows,peaks,boundaries,boundary_states=NamedTuple[],NamedTuple[],NamedTuple[],NamedTuple[]
    branch_paths=Dict(String(d)=>joinpath(directory,String(d),"bifurcation_summary.csv")
        for d in unique(p.direction for p in plan.plateaus))
    function save_completed()
        write_rows(joinpath(directory,"bifurcation_summary.csv"),rows,SUMMARY_COLUMNS)
        write_rows(joinpath(directory,"bifurcation_local_extrema.csv"),peaks,PEAK_COLUMNS)
        write_rows(joinpath(directory,"plateau_boundaries.csv"),boundaries,BOUNDARY_COLUMNS)
        for (direction,path) in branch_paths
            write_rows(path,filter(r -> r.direction == direction,rows),SUMMARY_COLUMNS)
            write_rows(joinpath(dirname(path),"bifurcation_local_extrema.csv"),
                filter(r -> r.direction == direction,peaks),PEAK_COLUMNS)
        end
    end
    save_completed()
    try
        start_speed=first(plan.plateaus).airspeed_mps
        common=(;density,angle_of_attack,chordwise_panels,spanwise_panels,symmetric_wing,
            maximum_wake_rows=wake_rows,core_radius,wake_shedding_fraction,aerodynamic_interaction)
        equilibrium=B.run_pazy_wing_uvlm_static(;equilibrium_options...,common...,
            airspeed=start_speed,time_step_chords=plan.dt*start_speed/chord)
        writedlm(joinpath(directory,"initialization.csv"),
            vcat(permutedims(equilibrium.history_columns),equilibrium.history),',')
        equilibrium.converged || error("Static initialization failed: $(equilibrium.termination_reason)")
        initial_static=B.model_responses(equilibrium.structural_problem.model,angle_of_attack,span,gauge_y,gauge_z)
        values=Dict(r=>fill(NaN,plan.steps_per_plateau+1) for r in responses)
        for r in responses
            values[r][1]=getproperty(initial_static,r)
        end
        snapshots(structure,aerodynamic,t)=(;time=t,structural_state=copy(structure.x),
            aerodynamic_snapshot=B.UVLM.snapshot_uvlm(aerodynamic.state.system))
        function observe(i,t,structure,aerodynamic)
            step=i-1
            k=plateau_at_step(plan,step)
            plateau=plan.plateaus[k]
            local_step=step-plateau.start_step
            speed=aerodynamic.operatingPoint.airspeed
            speed == plateau.airspeed_mps || error("Plateau speed mismatch at step $step")
            response=B.model_responses(structure.model,angle_of_attack,span,gauge_y,gauge_z)
            for r in responses
                values[r][local_step+1]=getproperty(response,r)
            end
            if step == plateau.end_step
                analyzed=analyze_plateau(plan,plateau,values,initial_static;
                    responses,angle_of_attack,extrema_method,check_convergence,rtol,convergence_cycles,response_definitions)
                append!(rows,analyzed.rows)
                append!(peaks,analyzed.peaks)
                push!(boundaries,(;plateau_index=k,direction=String(plateau.direction),airspeed_mps=speed,
                    global_time_s=t,accepted_step=step,aerodynamic_step=aerodynamic.state.stepIndex,
                    aerodynamic_time_s=aerodynamic.state.timeNow,retained_wake_rows=aerodynamic.state.system.nwake[1],
                    circulation_norm=norm(aerodynamic.state.system.Γ),structural_state_norm=norm(structure.x),
                    maximum_section_velocity_norm=maximum(norm(vcat(e.states.V,e.states.Ω)) for e in structure.model.elements),
                    structural_problem_id=objectid(structure),aerodynamic_problem_id=objectid(aerodynamic)))
                save_boundary_states && push!(boundary_states,snapshots(structure,aerodynamic,t))
                prefix=joinpath(directory,"plateau_$(lpad(k,4,'0'))")
                write_rows(prefix*"_extrema.csv",analyzed.rows,SUMMARY_COLUMNS)
                write_rows(prefix*"_local_extrema.csv",analyzed.peaks,PEAK_COLUMNS)
                if save_time_histories
                    open(prefix*"_history.csv","w") do io
                        writedlm(io,permutedims(vcat(["global_time_s","plateau_time_s","airspeed_mps"],String.(collect(responses)))),',')
                        # Exclude the shared pre-change boundary from the constant-speed history.
                        writedlm(io,hcat(plateau.start_time_s.+plan.local_time[2:end],plan.local_time[2:end],
                            fill(speed,plan.steps_per_plateau),(values[r][2:end] for r in responses)...),',')
                    end
                end
                save_completed()
                bad=filter(r -> !r.converged,analyzed.rows)
                println("Saved plateau $k/$(length(plan.plateaus)): $(plateau.direction), U=$speed m/s; ",
                    join([r.response*"="*r.convergence_status for r in analyzed.rows],", "))
                flush(stdout)
                if !isempty(bad) && on_unconverged != :continue
                    message="Late-window convergence unverified at plateau $k, U=$speed: "*join(getproperty.(bad,:response),", ")
                    on_unconverged == :error && error(message*". Completed plateau data were saved.")
                    @warn message
                end
                for r in responses
                    carried=last(values[r])
                    fill!(values[r],NaN)
                    values[r][1]=carried
                end
            end
        end
        profile(t)=plan.plateaus[plateau_at_step(plan,round(Int,t/plan.dt))].airspeed_mps
        result=B.run_pazy_wing_uvlm(;dynamic_options...,common...,
            airspeed=U_max,airspeed_profile=profile,initial_equilibrium=equilibrium,
            time_step_chords=plan.dt*U_max/chord,duration=plan.total_time,sideslip=0.0,
            initial_airspeed_fraction=1.0,airspeed_ramp_duration=0.0,
            settling_time=pulse_start,perturbation_amplitude=pulse_amplitude,perturbation_duration=pulse_duration,
            accepted_step_callback=observe,save_uvlm_history=false,uvlm_save_frequency=1,
            animation_frames=2,save_animation_history=false,save_structural_history=false,
            save_aerodynamic_load_history=false)
        length(boundaries) == length(plan.plateaus) || error("Integration did not complete every plateau")
        open(joinpath(directory,"settings.txt"),"a") do io
            println(io,"Completed all plateaus. Global final time=$(last(result.time)); accepted steps=$(length(result.time)-1).")
        end
        plot_paths=String[]
        if plot_results
            isdefined(@__MODULE__,:PazySteppedBifurcationPlots) ||
                Base.include(@__MODULE__,joinpath(@__DIR__,"PazySteppedBifurcationPlots.jl"))
            plot_paths=Base.invokelatest(PazySteppedBifurcationPlots.plot_stepped_bifurcation,
                directory;plot_options...)
        end
        return (;plan,dt=plan.dt,maximum_wake_rows=wake_rows,rows,local_extrema=peaks,
            boundaries,boundary_states,summary_path=joinpath(directory,"bifurcation_summary.csv"),
            branch_summary_paths=branch_paths,output_directory=directory,plot_paths,
            solver_result=keep_solver_result ? result : nothing,
            initial_equilibrium=keep_solver_result ? equilibrium : nothing)
    catch exception
        open(joinpath(directory,"failure.txt"),"w") do io
            println(io,"Completed plateaus: ",length(boundaries)," / ",length(plan.plateaus))
            showerror(io,exception,catch_backtrace())
        end
        rethrow()
    end
end

end # module
