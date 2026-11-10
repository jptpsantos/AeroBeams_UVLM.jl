"""Fit a perturbation growth exponent from nonoverlapping RMS windows.

Use a *difference* between phase-matched perturbed and unperturbed trajectories.
A periodically forced rotor response on its own is not a flutter envelope.
The returned R² and number of windows help identify an inadequate fit interval.
"""
function perturbation_growth_rate(time,signal;start_time,window_duration)
    length(time) == length(signal) || throw(DimensionMismatch("Signal/time lengths differ"))
    window_duration > 0 || throw(ArgumentError("Positive RMS window duration required"))
    centers,amplitudes = Float64[],Float64[]
    left = Float64(start_time)
    while left+window_duration <= last(time)
        indices = findall(t -> left <= t < left+window_duration,time)
        if length(indices) >= 3
            amplitude = sqrt(sum(abs2,signal[indices])/length(indices))
            if isfinite(amplitude) && amplitude > 100eps(Float64)
                push!(centers,left+window_duration/2)
                push!(amplitudes,amplitude)
            end
        end
        left += window_duration
    end
    length(centers) >= 4 || throw(ArgumentError("At least four resolved RMS windows required"))
    logamp = log.(amplitudes)
    A = hcat(ones(length(centers)),centers)
    fit = A\logamp
    variance = sum(abs2,logamp.-sum(logamp)/length(logamp))
    r_squared = variance > 0 ? 1-sum(abs2,logamp.-A*fit)/variance : NaN
    return (;growth_rate_per_s=fit[2],r_squared,window_centers=centers,rms_amplitudes=amplitudes)
end

"""Run phase-matched baseline/perturbed cases independently at each airspeed.

Each speed gets its own mean trim and developed wake. This is a small-signal
growth study, not branch continuation or a Floquet eigenanalysis. The complete
case results remain available for independent signal and convergence checks.
"""
function run_chang_aerobeams_speed_sweep(config,speeds;
    fit_start_s=config.solver.impulse_start_s+config.solver.impulse_duration_s+0.2,
    window_duration_s=0.1)
    config.initialization.mode == :mean_equilibrium || throw(ArgumentError(
        "Use independently converged mean equilibria for this speed sweep"))
    rows = NamedTuple[]
    runs = Any[]
    for speed in speeds
        speed > 0 || throw(ArgumentError("Positive airspeeds required"))
        reference = merge(config.reference,(simulation=merge(config.reference.simulation,
            (freestream_speed_mps=Float64(speed),)),))
        localconfig = merge(config,(reference=reference,output=merge(config.output,
            (write_results=false,plot_results=false,)),))
        baseline_config = merge(localconfig,(solver=merge(localconfig.solver,
            (impulse_magnitude_nm=0.0,)),))
        baseline = run_chang_aerobeams(baseline_config)
        perturbed = run_chang_aerobeams(localconfig)
        isapprox(baseline.phase_time_offset,perturbed.phase_time_offset;atol=1e-12,rtol=0) ||
            error("Baseline/perturbed initialization phases differ")
        perturbation = perturbed.results[:rotor_angles].-baseline.results[:rotor_angles]
        for rotor in eachindex(perturbed.rotors), coordinate in 1:2
            fit = perturbation_growth_rate(perturbed.time,perturbation[2rotor-2+coordinate,:];
                start_time=fit_start_s,window_duration=window_duration_s)
            push!(rows,(;speed_mps=Float64(speed),rotor,coordinate=coordinate == 1 ? "pitch" : "yaw",
                growth_rate_per_s=fit.growth_rate_per_s,r_squared=fit.r_squared))
        end
        push!(runs,(;baseline,perturbed))
    end
    if config.output.write_results
        mkpath(config.output.directory)
        path = joinpath(config.output.directory,config.output.label*"_speed_sweep.csv")
        open(path,"w") do io
            println(io,"speed_mps,rotor,coordinate,growth_rate_per_s,r_squared")
            for row in rows
                println(io,join(values(row),','))
            end
        end
    end
    return (;rows,runs)
end
