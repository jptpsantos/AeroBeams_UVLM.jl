module MovingBlockDamping
using FFTW, Statistics
using LinearAlgebra: mul!

export moving_block_damping, moving_block_metrics

function unavailable(reason; peak_count=0, block_size=0)
    return (;summary=(;available=false, valid_for_convergence=false, reason,
        peak_count, block_size, block_count=0, moving_block_lambda_per_s=NaN,
        growth_rate_per_s=NaN, frequency_hz=NaN, omega_radps=NaN,
        frequency_resolution_hz=NaN, damping_ratio_percent=NaN,
        fit_r_squared=NaN, envelope_ratio=NaN, analyzed_start_s=NaN,
        analyzed_end_s=NaN),
        diagnostics=(;time_s=Float64[], signal=Float64[], block_start_s=Float64[],
            log_amplitude=Float64[], fitted_log_amplitude=Float64[], frequency_hz=Float64[]))
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
function moving_block_metrics(t::AbstractVector, x::AbstractVector;
    block_size::Int=512, size_ratio_lb::Real=0.25, size_ratio_ub::Real=0.50,
    peak_from_start::Int=1, peak_from_end::Int=0,
    minimum_peaks::Int=2, minimum_fit_r_squared::Real=0.0,
    fit_start_s::Real=-Inf, fit_end_s::Real=Inf,
    block_duration_s=nothing, block_overlap=nothing, frequency_band_hz=nothing)

    length(t)==length(x) || throw(DimensionMismatch("Time and signal lengths differ"))
    block_size>=2 || throw(ArgumentError("block_size must be at least 2"))
    0<size_ratio_lb<=size_ratio_ub<1 || throw(ArgumentError("Require 0 < size_ratio_lb <= size_ratio_ub < 1"))
    peak_from_start>=1 && peak_from_end>=0 && minimum_peaks>=2 ||
        throw(ArgumentError("Invalid peak selection"))
    0<=minimum_fit_r_squared<=1 || throw(ArgumentError("Fit R² threshold must be in [0,1]"))
    fit_start_s<fit_end_s || throw(ArgumentError("Fit start must precede fit end"))
    isnothing(block_duration_s) || (isfinite(block_duration_s) && block_duration_s>0) ||
        throw(ArgumentError("Block duration must be finite and positive"))
    isnothing(block_overlap) || 0<=block_overlap<1 || throw(ArgumentError("Block overlap must be in [0,1)"))
    isnothing(frequency_band_hz) || (length(frequency_band_hz)==2 &&
        0<=frequency_band_hz[1]<frequency_band_hz[2]) || throw(ArgumentError("Invalid frequency band"))
    length(t)>=3 || return unavailable("Insufficient samples")
    all(isfinite,t) && all(isfinite,x) || throw(ArgumentError("Time and signal must be finite"))
    intervals=diff(t)
    all(>(0),intervals) || throw(ArgumentError("Time must be strictly increasing"))
    dt=mean(intervals)
    time_tolerance=100eps(Float64)*max(1,maximum(abs,t))
    all(h -> isapprox(h,dt;rtol=1e-7,atol=time_tolerance),intervals) ||
        throw(ArgumentError("Moving-block FFT requires uniformly sampled time"))

    keep=findall(i -> fit_start_s<=t[i]<=fit_end_s,eachindex(t))
    length(keep)>=3 || return unavailable("Insufficient samples in fitting interval")
    time=Float64.(t[keep]); signal=Float64.(x[keep])
    # Preserve the reference's positive-peak threshold and peak-based cut.
    peak_min=maximum(abs,signal)/100
    peaks=[i for i in 2:length(signal)-1 if signal[i]>signal[i-1] &&
        signal[i]>=signal[i+1] && signal[i]>=peak_min]
    count=length(peaks)
    count>=max(minimum_peaks,peak_from_start+peak_from_end+1) ||
        return unavailable("Insufficient positive peaks";peak_count=count)
    cut=peaks[peak_from_start]:peaks[end-peak_from_end]
    time=time[cut]; signal=signal[cut]; n=length(time)
    dt=mean(diff(time))
    if isnothing(block_duration_s)
        while block_size/n<size_ratio_lb
            block_size*=2
        end
        while block_size/n>size_ratio_ub
            block_size=div(block_size,2)
        end
        size_ratio_lb<=block_size/n<=size_ratio_ub ||
            return unavailable("No doubled/halved block size meets the requested ratios";peak_count=count,block_size)
    else
        block_size=round(Int,block_duration_s/dt)
    end
    2<=block_size<n || return unavailable("Block must contain at least two samples and leave multiple windows";peak_count=count,block_size)
    stride=isnothing(block_overlap) ? 1 : max(1,round(Int,block_size*(1-block_overlap)))
    starts=collect(1:stride:n-block_size+1)
    last(starts)==n-block_size+1 || push!(starts,n-block_size+1)
    length(starts)>=2 || return unavailable("Insufficient moving blocks";peak_count=count,block_size)

    # Reuse one real FFT plan and input/output buffers for all windows.
    buffer=zeros(block_size)
    transform=plan_rfft(buffer;flags=FFTW.ESTIMATE)
    spectrum=transform*buffer
    log_amplitude=Vector{Float64}(undef,length(starts))
    frequency=similar(log_amplitude)
    spacing=1/(block_size*dt)
    for (k,start) in enumerate(starts)
        copyto!(buffer,1,signal,start,block_size)
        buffer .-= mean(buffer)
        mul!(spectrum,transform,buffer)
        # Double positive-frequency bins except Nyquist for an even FFT size.
        amplitudes=abs.(spectrum)./block_size
        last_doubled=iseven(block_size) ? length(amplitudes)-1 : length(amplitudes)
        amplitudes[2:last_doubled].*=2
        index=argmax(@view amplitudes[2:end])+1
        log_amplitude[k]=log(max(amplitudes[index],eps(Float64)))
        frequency[k]=(index-1)*spacing
    end
    block_time=time[starts]
    centered_time=block_time .- mean(block_time)
    centered_log=log_amplitude .- mean(log_amplitude)
    lambda=sum(centered_time.*centered_log)/sum(abs2,centered_time)
    fitted=mean(log_amplitude) .+ lambda.*centered_time
    total=sum(abs2,centered_log)
    r_squared=total>eps(Float64) ? 1-sum(abs2,log_amplitude.-fitted)/total : 1.0
    f=median(frequency); omega=2pi*f
    in_band=isnothing(frequency_band_hz) || all(v -> frequency_band_hz[1]<=v<=frequency_band_hz[2],frequency)
    valid=isfinite(lambda) && r_squared>=minimum_fit_r_squared && in_band
    reason=!in_band ? "Dominant frequency leaves the requested band" :
        r_squared<minimum_fit_r_squared ? "Moving-block log-amplitude fit below R² threshold" : "Moving-block fit completed"
    summary=(;available=true,valid_for_convergence=valid,reason,peak_count=count,
        block_size,block_count=length(starts),moving_block_lambda_per_s=lambda,
        growth_rate_per_s=lambda,frequency_hz=f,omega_radps=omega,
        frequency_resolution_hz=spacing,damping_ratio_percent=-100lambda/hypot(lambda,omega),
        fit_r_squared=r_squared,envelope_ratio=exp(lambda*(block_time[end]-block_time[1])),
        analyzed_start_s=first(time),analyzed_end_s=last(time))
    return (;summary,diagnostics=(;time_s=time,signal,block_start_s=block_time,
        log_amplitude,fitted_log_amplitude=fitted,frequency_hz=frequency))
end

"""Return `(lambda, frequency_hz, omega_radps)`, as in moving_block_2.jl."""
function moving_block_damping(t::AbstractVector,x::AbstractVector;kwargs...)
    result=moving_block_metrics(t,x;kwargs...).summary
    result.available || throw(ArgumentError(result.reason))
    return result.moving_block_lambda_per_s,result.frequency_hz,result.omega_radps
end
end
