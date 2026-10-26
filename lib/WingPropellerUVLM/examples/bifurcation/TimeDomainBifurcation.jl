module TimeDomainBifurcation

using DelimitedFiles, Statistics, Distributed, Dates
export bifurcation_speeds, response_statistics, response_local_extrema, run_bifurcation, plot_bifurcation

"""Positive inclusive grid; append U_max if dU does not divide the interval."""
function bifurcation_speeds(U_min, U_max, dU)
    all(isfinite, (U_min, U_max, dU)) && 0 < U_min <= U_max && dU > 0 ||
        throw(ArgumentError("Require finite 0 < U_min <= U_max and dU > 0"))
    speeds = collect(Float64(U_min):Float64(dU):Float64(U_max))
    isapprox(last(speeds), U_max; atol=10eps(Float64(U_max)), rtol=0) ||
        push!(speeds, Float64(U_max))
    return speeds
end

# Include flat-topped peaks once, at the last sample of the plateau.
function turning_points(q)
    maxima, minima = Int[], Int[]
    previous = 0
    for i in 1:length(q)-1
        direction = sign(q[i+1] - q[i])
        iszero(direction) && continue
        previous > 0 && direction < 0 && push!(maxima, i)
        previous < 0 && direction > 0 && push!(minima, i)
        previous = direction
    end
    return maxima, minima
end

"""Raw extrema after transient removal; conservative cycle-extrema stationarity check.

`atol` is in response units. `rtol` scales with oscillation amplitude, not the
static offset. This diagnostic does not establish asymptotic stability. Signals
with too few complete cycles remain `insufficient_cycles`, unless their entire
window range is within 2atol of a stationary response.
"""
function response_statistics(time, q; transient_time, check_convergence=true,
    atol=1e-6, rtol=0.02, convergence_cycles=3)
    length(time) == length(q) >= 2 || throw(ArgumentError("Need matching histories with at least two samples"))
    all(isfinite, time) && all(isfinite, q) && all(diff(time) .> 0) ||
        throw(ArgumentError("Histories must be finite and time strictly increasing"))
    isfinite(transient_time) && first(time) <= transient_time < last(time) ||
        throw(ArgumentError("transient_time must lie before the final sample"))
    isfinite(atol) && atol >= 0 && isfinite(rtol) && rtol >= 0 &&
        convergence_cycles isa Integer && convergence_cycles >= 2 ||
        throw(ArgumentError("Convergence tolerances must be nonnegative; cycles >= 2"))
    start = searchsortedfirst(time, transient_time)
    length(time)-start+1 >= 2 || throw(ArgumentError("Observation window needs at least two samples"))
    tail = @view q[start:end]
    low, high = extrema(tail)
    amplitude = (high-low)/2
    status = "not_checked"
    if check_convergence
        if high-low <= 2atol
            status = "stationary"
        else
            peaks, _ = turning_points(tail)
            if length(peaks) >= convergence_cycles+1
                peaks = peaks[end-convergence_cycles:end]
                cycles = [extrema(@view tail[a:b]) for (a,b) in zip(peaks[1:end-1], peaks[2:end])]
                lows, highs = first.(cycles), last.(cycles)
                tolerance = atol + rtol * mean((highs .- lows) ./ 2)
                status = max(maximum(lows)-minimum(lows), maximum(highs)-minimum(highs)) <= tolerance ?
                    "converged" : "unconverged"
            else
                status = "insufficient_cycles"
            end
        end
    end
    return (; minimum=low, maximum=high, mean=mean(tail), amplitude,
        convergence_status=status, converged=status in ("stationary", "converged"),
        window_start_s=time[start], window_end_s=last(time), sample_count=length(tail))
end

"""Every interior turning point in the final window, with its original sample time.

Flat peaks are included once, at their last sample. Window endpoints and
monotone/stationary windows do not create artificial local extrema.
"""
function response_local_extrema(time, q; transient_time)
    response_statistics(time, q; transient_time, check_convergence=false)
    start = searchsortedfirst(time, transient_time)
    tail = @view q[start:end]
    maxima, minima = turning_points(tail)
    points = [(; kind, time_s=time[start+i-1], value=tail[i])
        for (kind, indices) in (("minimum", minima), ("maximum", maxima)) for i in indices]
    return sort!(points; by=p -> p.time_s)
end

const SUMMARY_COLUMNS = ["airspeed_mps", "response", "unit", "minimum", "maximum",
    "mean", "amplitude", "convergence_status", "converged", "window_start_s",
    "window_end_s", "sample_count", "static_equilibrium", "label"]

function write_summary(path, rows)
    open(path, "w") do io
        writedlm(io, permutedims(SUMMARY_COLUMNS), ',')
        for r in rows
            writedlm(io, permutedims(Any[r.airspeed_mps, r.response, r.unit,
                r.minimum, r.maximum, r.mean, r.amplitude, r.convergence_status,
                r.converged, r.window_start_s, r.window_end_s, r.sample_count,
                r.static_equilibrium, r.label]), ',')
        end
    end
end

const LOCAL_EXTREMA_COLUMNS = ["airspeed_mps", "response", "unit", "kind", "time_s", "value"]

function write_local_extrema(path, rows)
    open(path, "w") do io
        writedlm(io, permutedims(LOCAL_EXTREMA_COLUMNS), ',')
        for r in rows
            writedlm(io, permutedims(Any[r.airspeed_mps, r.response, r.unit,
                r.kind, r.time_s, r.value]), ',')
        end
    end
end

function analyze_speed(job, simulate, settings, directory)
    index, speed = job
    # simulate must construct new structural/aerodynamic/wake states on EVERY call.
    result = simulate(speed)
    time = result.time
    isapprox(first(time), 0; atol=1e-12, rtol=0) &&
        isapprox(last(time), settings.simulation_time; atol=1e-10, rtol=1e-10) ||
        throw(ArgumentError("Each independent simulation must cover [0, simulation_time]"))
    local_extrema = Any[]
    rows = map(settings.responses) do name
        hasproperty(result.responses, name) || throw(ArgumentError("Response $name unavailable"))
        definition = getproperty(settings.definitions, name)
        statistics = response_statistics(time, getproperty(result.responses, name);
            settings.transient_time, settings.check_convergence, settings.rtol,
            atol=definition.atol, settings.convergence_cycles)
        if settings.save_local_extrema
            for point in response_local_extrema(time, getproperty(result.responses, name);
                settings.transient_time)
                push!(local_extrema, (; airspeed_mps=speed, response=String(name),
                    unit=definition.unit, point...))
            end
        end
        static_value = hasproperty(result, :static_equilibrium) ?
            get(result.static_equilibrium, name, NaN) : NaN
        (; airspeed_mps=speed, response=String(name), unit=definition.unit,
            statistics..., static_equilibrium=static_value, label=definition.label)
    end
    # Unique case files also preserve completed cases if another solve fails.
    prefix = joinpath(directory, "case_$(lpad(index, 4, '0'))")
    write_summary(prefix * "_extrema.csv", rows)
    settings.save_local_extrema && write_local_extrema(prefix * "_local_extrema.csv", local_extrema)
    if settings.save_time_histories
        open(prefix * "_history.csv", "w") do io
            writedlm(io, permutedims(vcat(["time_s"], String.(settings.responses))), ',')
            writedlm(io, hcat(time, (getproperty(result.responses, n) for n in settings.responses)...), ',')
        end
    end
    return (; rows, local_extrema) # Release full histories after extracting peaks.
end

"""Run fresh fixed-speed simulations through `simulate(U)`.

The callback returns `(time, responses, static_equilibrium)`; responses is a
NamedTuple of selected scalar vectors, and static_equilibrium is optional.
For process parallelism, load the adapter/module on `worker_ids` first. Workers
must see the same absolute output directory. No worker creates plots.
"""
function run_bifurcation(simulate; U_min, U_max, dU, definitions,
    responses=collect(keys(definitions)), simulation_time=10.0, transient_time=8.0,
    check_convergence=true, rtol=0.02, convergence_cycles=3,
    save_time_histories=false, save_local_extrema=true, plot_results=true, plot_amplitude=true,
    zoom_limits=nothing, worker_ids=Int[], metadata=(;),
    output_directory=joinpath(pwd(), "output", "bifurcation_" * Dates.format(now(), "yyyymmdd_HHMMSSsss")))
    speeds = bifurcation_speeds(U_min, U_max, dU)
    all(isfinite, (simulation_time, transient_time, rtol)) &&
        0 <= transient_time < simulation_time && rtol >= 0 ||
        throw(ArgumentError("Require 0 <= transient_time < simulation_time and rtol >= 0"))
    convergence_cycles isa Integer && convergence_cycles >= 2 ||
        throw(ArgumentError("convergence_cycles must be an integer >= 2"))
    selected = collect(Symbol, responses)
    !isempty(selected) && allunique(selected) || throw(ArgumentError("Select distinct response names"))
    for name in selected
        hasproperty(definitions, name) || throw(ArgumentError("Unknown response $name"))
        occursin(r"^[A-Za-z][A-Za-z0-9_]*$", String(name)) || throw(ArgumentError("Invalid response name"))
        d = getproperty(definitions, name)
        isfinite(d.atol) && d.atol >= 0 || throw(ArgumentError("Invalid absolute tolerance for $name"))
    end
    isnothing(zoom_limits) || (length(zoom_limits) == 2 && all(isfinite, zoom_limits) &&
        zoom_limits[1] < zoom_limits[2]) || throw(ArgumentError("zoom_limits must be (U_low, U_high)"))
    all(id -> id != myid() && id in workers(), worker_ids) ||
        throw(ArgumentError("worker_ids must name active remote workers"))
    directory = abspath(output_directory)
    mkpath(directory)
    isempty(readdir(directory)) || throw(ArgumentError("Use an empty output directory to preserve run provenance: $directory"))
    settings = (; responses=selected, definitions, simulation_time, transient_time,
        check_convergence, rtol, convergence_cycles, save_time_histories, save_local_extrema)
    open(joinpath(directory, "settings.txt"), "w") do io
        println(io, "Independent constant-speed simulations; no state transfer between speeds.")
        println(io, "Speeds [m/s]: ", speeds)
        println(io, "Analysis: ", settings)
        println(io, "Worker IDs: ", worker_ids)
        println(io, "Model settings: ", metadata)
    end
    jobs = collect(enumerate(speeds))
    rows = Any[]
    local_extrema = Any[]
    summary_path = joinpath(directory, "bifurcation_summary.csv")
    local_extrema_path = save_local_extrema ? joinpath(directory, "bifurcation_local_extrema.csv") : nothing
    operation = job -> analyze_speed(job, simulate, settings, directory)
    if isempty(worker_ids)
        for job in jobs
            case = operation(job)
            append!(rows, case.rows)
            append!(local_extrema, case.local_extrema)
            write_summary(summary_path, rows)
            save_local_extrema && write_local_extrema(local_extrema_path, local_extrema)
        end
    else
        for case in pmap(operation, WorkerPool(worker_ids), jobs)
            append!(rows, case.rows)
            append!(local_extrema, case.local_extrema)
        end
        write_summary(summary_path, rows)
        save_local_extrema && write_local_extrema(local_extrema_path, local_extrema)
    end
    if plot_results
        plot_bifurcation(summary_path; output_directory=directory, plot_amplitude, zoom_limits)
    end
    return (; speeds, rows, summary_path, local_extrema_path, output_directory=directory)
end

"""Regenerate figures from CSV alone; means are never substituted for static equilibria."""
function plot_bifurcation(summary_path; output_directory=dirname(summary_path),
    plot_amplitude=true, zoom_limits=nothing)
    # Load plotting only for requested figures; serial and worker solves need no GR state.
    if !isdefined(@__MODULE__, :Plots)
        @eval import Plots
    end
    # invokelatest handles the first dynamic import's world age.
    return Base.invokelatest(plot_bifurcation_loaded, summary_path;
        output_directory, plot_amplitude, zoom_limits)
end

function plot_bifurcation_loaded(summary_path; output_directory, plot_amplitude, zoom_limits)
    data, header = readdlm(summary_path, ','; header=true)
    columns = Dict(String(n) => i for (i,n) in enumerate(vec(header)))
    paths = String[]
    mkpath(output_directory)
    for name in unique(data[:,columns["response"]])
        group = data[data[:,columns["response"]] .== name, :]
        group = group[sortperm(Float64.(group[:,columns["airspeed_mps"]])), :]
        values(key) = Float64.(group[:,columns[key]])
        speed = values("airspeed_mps")
        ylabel = "$(group[1,columns["label"]]) [$(group[1,columns["unit"]])]"
        p = Plots.plot(speed, values("minimum"); marker=:dtriangle, label="Minimum",
            xlabel="Freestream velocity [m/s]", ylabel, linewidth=1.2,
            size=(800,500), left_margin=5Plots.mm, bottom_margin=4Plots.mm,
            legend=:outertopright)
        Plots.plot!(p, speed, values("maximum"); marker=:utriangle, label="Maximum", linewidth=1.2)
        static = values("static_equilibrium")
        any(isfinite, static) && Plots.plot!(p, speed, static;
            label="Static equilibrium", linestyle=:dash, color=:gray)
        flagged = [string(s) != "true" for s in group[:,columns["converged"]]]
        if any(flagged)
            Plots.scatter!(p, vcat(speed[flagged], speed[flagged]),
                vcat(values("minimum")[flagged], values("maximum")[flagged]);
                marker=:xcross, color=:red, label="Not convergence-verified")
        end
        path = joinpath(output_directory, "bifurcation_$(name).png")
        Plots.savefig(p, path); push!(paths, path)
        if !isnothing(zoom_limits)
            Plots.xlims!(p, zoom_limits)
            path = joinpath(output_directory, "bifurcation_$(name)_zoom.png")
            Plots.savefig(p, path); push!(paths, path)
        end
        if plot_amplitude
            p = Plots.plot(speed, values("amplitude"); marker=:circle,
                xlabel="Freestream velocity [m/s]",
                ylabel="Oscillation amplitude [$(group[1,columns["unit"]])]",
                title=String(group[1,columns["label"]]), titlefontsize=10,
                label="Final-window half peak-to-peak", size=(800,500), legend=:outertopright)
            any(flagged) && Plots.scatter!(p, speed[flagged], values("amplitude")[flagged];
                marker=:xcross, color=:red, label="Not convergence-verified")
            path = joinpath(output_directory, "bifurcation_$(name)_amplitude.png")
            Plots.savefig(p, path); push!(paths, path)
        end
    end
    return paths
end

end # module
