module MovingBlockDamping
using FFTW,Statistics,Plots
using LinearAlgebra: mul!

export moving_block_damping,moving_block_metrics

function unavailable(reason; peak_count=0,block_size=0)
    return (;
        summary=(; available=false,valid_for_convergence=false,reason,
            peak_count,block_size,block_count=0,moving_block_lambda_per_s=NaN,
            growth_rate_per_s=NaN,frequency_hz=NaN,omega_radps=NaN,
            frequency_resolution_hz=NaN,damping_ratio_percent=NaN,
            fit_r_squared=NaN,envelope_ratio=NaN,analyzed_start_s=NaN,
            analyzed_end_s=NaN),
        diagnostics=(; time_s=Float64[],signal=Float64[],block_start_s=Float64[],
            log_amplitude=Float64[],fitted_log_amplitude=Float64[],frequency_hz=Float64[]))
end

"""
    moving_block_metrics(t, x; ...)

Moving-block FFT/log-amplitude regression following `moving_block_2.jl`.
Default blocks occupy 25--50% of the peak-trimmed record, shift one sample,
remove their own mean, and use the dominant non-DC FFT component. No Hann
window is applied. The frequency band, if supplied, checks the measured
dominant frequency; it never restricts the search to spectral leakage in-band.

`block_duration_s` and `block_overlap` optionally select a fixed physical
window and fractional stride for time-step studies. Frequency is an FFT bin
estimate; `frequency_resolution_hz` reports the spacing, not an error bound.
"""
function moving_block_metrics(t::AbstractVector,x::AbstractVector;
    block_size::Int=512,size_ratio_lb::Real=0.25,size_ratio_ub::Real=0.50,
    peak_from_start::Int=1,peak_from_end::Int=0,
    minimum_peaks::Int=2,minimum_fit_r_squared::Real=0.0,
    fit_start_s::Real=-Inf,fit_end_s::Real=Inf,
    block_duration_s=nothing,block_overlap=nothing,frequency_band_hz=nothing)

    length(t) == length(x) || throw(DimensionMismatch("Time and signal lengths differ"))
    block_size >= 2 || throw(ArgumentError("block_size must be at least 2"))
    0 < size_ratio_lb <= size_ratio_ub < 1 || throw(ArgumentError("Require 0 < size_ratio_lb <= size_ratio_ub < 1"))
    peak_from_start >= 1 && peak_from_end >= 0 && minimum_peaks >= 2 ||
        throw(ArgumentError("Invalid peak selection"))
    0 <= minimum_fit_r_squared <= 1 || throw(ArgumentError("Fit R² threshold must be in [0,1]"))
    fit_start_s < fit_end_s || throw(ArgumentError("Fit start must precede fit end"))
    isnothing(block_duration_s) || (isfinite(block_duration_s) && block_duration_s > 0) ||
        throw(ArgumentError("Block duration must be finite and positive"))
    isnothing(block_overlap) || 0 <= block_overlap < 1 || throw(ArgumentError("Block overlap must be in [0,1)"))
    isnothing(frequency_band_hz) ||
        (length(frequency_band_hz) == 2 &&
         0 <= frequency_band_hz[1] < frequency_band_hz[2]) || throw(ArgumentError("Invalid frequency band"))
    length(t) >= 3 || return unavailable("Insufficient samples")
    all(isfinite,t) && all(isfinite,x) || throw(ArgumentError("Time and signal must be finite"))
    intervals = diff(t)
    all(>(0),intervals) || throw(ArgumentError("Time must be strictly increasing"))
    dt = mean(intervals)
    timeTolerance = 100eps(Float64) * max(1,maximum(abs,t))
    all(h -> isapprox(h,dt; rtol=1e-7,atol=timeTolerance),intervals) ||
        throw(ArgumentError("Moving-block FFT requires uniformly sampled time"))

    keep = findall(i -> fit_start_s <= t[i] <= fit_end_s,eachindex(t))
    length(keep) >= 3 || return unavailable("Insufficient samples in fitting interval")
    time = Float64.(t[keep])
    signal = Float64.(x[keep])
    # Preserve the reference's positive-peak threshold and peak-based cut.
    peakMin = maximum(abs,signal) / 100
    peaks = [i for i = 2:length(signal)-1 if signal[i] > signal[i-1] &&
     signal[i] >= signal[i+1] && signal[i] >= peakMin]
    count = length(peaks)
    count >= max(minimum_peaks,peak_from_start + peak_from_end + 1) ||
        return unavailable("Insufficient positive peaks"; peak_count=count)
    cut = peaks[peak_from_start]:peaks[end-peak_from_end]
    time = time[cut]
    signal = signal[cut]
    n = length(time)
    dt = mean(diff(time))
    if isnothing(block_duration_s)
        while block_size / n < size_ratio_lb
            block_size *= 2
        end
        while block_size / n > size_ratio_ub
            block_size = div(block_size,2)
        end
        size_ratio_lb <= block_size / n <= size_ratio_ub ||
            return unavailable("No doubled/halved block size meets the requested ratios"; peak_count=count,block_size)
    else
        block_size = round(Int,block_duration_s / dt)
    end
    2 <= block_size < n || return unavailable(
        "Block must contain at least two samples and leave multiple windows";
        peak_count=count,
        block_size
    )
    stride = isnothing(block_overlap) ? 1 : max(1,round(Int,block_size * (1 - block_overlap)))
    starts = collect(1:stride:n-block_size+1)
    last(starts) == n - block_size + 1 || push!(starts,n - block_size + 1)
    length(starts) >= 2 || return unavailable("Insufficient moving blocks"; peak_count=count,block_size)

    # Reuse one real FFT plan and input/output buffers for all windows.
    buffer = zeros(block_size)
    transform = plan_rfft(buffer; flags=FFTW.ESTIMATE)
    spectrum = transform * buffer
    log_amplitude = Vector{Float64}(undef,length(starts))
    frequency = similar(log_amplitude)
    spacing = 1 / (block_size * dt)
    for (k,start) in enumerate(starts)
        copyto!(buffer,1,signal,start,block_size)
        buffer .-= mean(buffer)
        mul!(spectrum,transform,buffer)
        # Double positive-frequency bins except Nyquist for an even FFT size.
        amplitudes = abs.(spectrum) ./ block_size
        lastDoubled = iseven(block_size) ? length(amplitudes) - 1 : length(amplitudes)
        amplitudes[2:lastDoubled] .*= 2
        index = argmax(@view amplitudes[2:end]) + 1
        log_amplitude[k] = log(max(amplitudes[index],eps(Float64)))
        frequency[k] = (index - 1) * spacing
    end
    blockTime = time[starts]
    centeredTime = blockTime .- mean(blockTime)
    centeredLog = log_amplitude .- mean(log_amplitude)
    lambda = sum(centeredTime .* centeredLog) / sum(abs2,centeredTime)
    fitted = mean(log_amplitude) .+ lambda .* centeredTime
    total = sum(abs2,centeredLog)
    r_squared = total > eps(Float64) ? 1 - sum(abs2,log_amplitude .- fitted) / total : 1.0
    f = median(frequency)
    omega = 2pi * f
    inBand = isnothing(frequency_band_hz) || all(v -> frequency_band_hz[1] <= v <= frequency_band_hz[2],frequency)
    valid = isfinite(lambda) && r_squared >= minimum_fit_r_squared && inBand
    reason =
        !inBand ? "Dominant frequency leaves the requested band" :
        r_squared < minimum_fit_r_squared ? "Moving-block log-amplitude fit below R² threshold" :
        "Moving-block fit completed"
    summary = (; available=true,valid_for_convergence=valid,reason,peak_count=count,
        block_size,block_count=length(starts),moving_block_lambda_per_s=lambda,
        growth_rate_per_s=lambda,frequency_hz=f,omega_radps=omega,
        frequency_resolution_hz=spacing,damping_ratio_percent=-100lambda / hypot(lambda,omega),
        fit_r_squared=r_squared,envelope_ratio=exp(lambda * (blockTime[end] - blockTime[1])),
        analyzed_start_s=first(time),analyzed_end_s=last(time))
    return (;
        summary,
        diagnostics=(; time_s=time,signal,block_start_s=blockTime,
            log_amplitude,fitted_log_amplitude=fitted,frequency_hz=frequency)
    )
end

"""
    moving_block_damping(t, x; kwargs...)

Applies the moving block approach to estimate damping and dominant frequency
from a given time history.

Returns:
- `lambda_approx`: exponential growth/decay rate, in 1/s
- `frequency_approx`: dominant frequency, in Hz
- `omega_approx`: dominant angular frequency, in rad/s
"""
function moving_block_damping(t::AbstractVector,x::AbstractVector;
    block_size::Int=512,
    size_ratio_lb::Float64=0.25,
    size_ratio_ub::Float64=0.50,
    peak_from_start::Int=1,
    peak_from_end::Int=0)

    if length(t) != length(x)
        error("Time vector and signal vector must have the same length.")
    end

    peakMin = maximum(abs.(x)) / 100.0

    # Local maxima equivalent to Peaks.argmaxima, without an extra dependency.
    locsAll = [i for i = 2:length(x)-1 if x[i] > x[i-1] && x[i] >= x[i+1]]
    locs = filter(i -> x[i] >= peakMin,locsAll)

    if length(locs) < (peak_from_start + peak_from_end + 1)
        error("Not enough peaks found in the signal. Try adjusting peak_from_start/end.")
    end

    startIdx = locs[peak_from_start]
    endIdx = locs[end-peak_from_end]

    tTest = t[startIdx:endIdx]
    xTest = x[startIdx:endIdx]

    sampleSize = length(tTest)
    sizeRatio = block_size / sampleSize

    while sizeRatio < size_ratio_lb
        block_size *= 2
        sizeRatio = block_size / sampleSize
    end

    while sizeRatio > size_ratio_ub
        block_size = div(block_size,2)
        sizeRatio = block_size / sampleSize
    end

    nBlocks = sampleSize - block_size + 1

    if nBlocks < 1
        error("Block size is larger than the analyzed signal segment.")
    end

    dt = mean(diff(tTest))
    fs = 1.0 / dt

    lnXbar = zeros(Float64,nBlocks)
    freqBlock = zeros(Float64,nBlocks)

    for i = 1:nBlocks
        x_i = xTest[i:i+block_size-1]

        # Remove DC offset so the 0 Hz bin does not dominate.
        x_i = x_i .- mean(x_i)

        X_i = fft(x_i)
        Xbar_i = abs.(X_i ./ block_size)

        halfIdx = div(block_size,2) + 1
        Xbar_i = Xbar_i[1:halfIdx]
        Xbar_i[2:end-1] .*= 2.0

        freqs = (0:halfIdx-1) .* fs ./ block_size

        # Ignore DC bin when finding dominant oscillation.
        peakIdx = argmax(Xbar_i[2:end]) + 1

        lnXbar[i] = log(max(Xbar_i[peakIdx],eps()))
        freqBlock[i] = freqs[peakIdx]
    end

    tau = tTest[1:nBlocks]

    # Degree-one least squares, equivalent to Polynomials.fit(tau, lnXbar, 1).
    coeff = hcat(ones(nBlocks),tau) \ lnXbar
    p_fit(time) = coeff[1] + coeff[2] * time

    lambda_approx = coeff[2]
    frequency_approx = median(freqBlock)
    omega_approx = 2.0 * pi * frequency_approx

    pTest = plot(tTest,xTest,
        linestyle=:dash,
        linewidth=1,
        label="Analyzed Segment",
        title="Cut Signal",
        xlabel="Time (s)")

    pFitPlot = plot(tau,lnXbar,
        linewidth=1,
        label="lnXbar")

    plot!(pFitPlot,tau,p_fit.(tau),
        linestyle=:dash,
        linewidth=1,
        label="Linear Fit",
        xlabel="Block start time (s)",
        ylabel="Log of moving block function (-)",
        title="Moving Block Fit")

    display(plot(pTest,pFitPlot,layout=(2,1),size=(800,600)))

    return lambda_approx,frequency_approx,omega_approx
end
end
