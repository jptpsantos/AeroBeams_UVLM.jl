# Included by PazyFlutterSweep.
using FFTW,LinearAlgebra,Statistics

const GROWTH_DEFAULTS = (;
    windows=((0.20,1.20),(1.20,2.20)),block_time=0.24,
    frequency_band=(22.0,43.0),frequency_tolerance=3.0,
    rate_atol=0.15,rate_rtol=0.25,minimum_rate=0.05,minimum_log_change=0.06,
    minimum_r_squared=0.85,minimum_snr=5.0,scale_rtol=0.30,mode_rtol=0.45,
    noise_floor=(1e-8,1e-5))

const LCO_DEFAULTS = (; tail_time=1.0,minimum_cycles=20,
    amplitude_rtol=0.05,repeat_rtol=0.12,noise_multiplier=5.0,
    amplitude_floor=(1e-6,1e-3))

function detrend_signal(t,q)
    centeredTime = t .- mean(t)
    signal = q .- mean(q)
    signal .-= centeredTime .* (dot(centeredTime,signal) / sum(abs2,centeredTime))
    return signal
end

# Hann-weighted spectrum. Interpolate the peak so reporting is not quantized
# by the length of the fitting window.
function signal_frequency(t,q,band)
    n = length(t)
    dt = (last(t)-first(t)) / (n-1)
    taper = [sinpi((i-1)/(n-1))^2 for i in 1:n]
    spectrum = abs.(rfft(taper .* detrend_signal(t,q)))
    spacing = 1 / (n*dt)
    bins = [i for i in 2:length(spectrum) if band[1] <= (i-1)*spacing <= band[2]]
    isempty(bins) && return NaN
    i = bins[argmax(spectrum[bins])]
    if i < length(spectrum)
        a,b,c = spectrum[i-1],spectrum[i],spectrum[i+1]
        curvature = a - 2b + c
        if curvature < 0
            return (i-1 + clamp((a-c)/(2curvature),-0.5,0.5))*spacing
        end
    end
    return (i-1)*spacing
end

function log_envelope_fit(t,a)
    ln = log.(max.(a,floatmin(Float64)))
    centeredTime = t .- mean(t)
    centeredLog = ln .- mean(ln)
    rate = dot(centeredTime,centeredLog) / sum(abs2,centeredTime)
    total = sum(abs2,centeredLog)
    rSquared = total > eps() ? 1 - sum(abs2,centeredLog .- rate .* centeredTime)/total : 0.0
    return rate,rSquared
end

function modal_window(t,forced,unforced,window,options)
    firstIndex = searchsortedfirst(t,window[1])
    lastIndex = searchsortedfirst(t,window[2])-1
    indices = firstIndex:lastIndex
    timeWindow = t[indices]
    response = forced[indices,:] .- unforced[indices,:]
    blockSize = round(Int,options.block_time/(t[2]-t[1]))
    if length(indices) < blockSize+2 || length(indices) < 3
        return nothing
    end
    frequency = signal_frequency(timeWindow,view(response,:,2),options.frequency_band)
    isfinite(frequency) || return nothing
    bendingFrequency = signal_frequency(timeWindow,view(response,:,1),options.frequency_band)
    starts = collect(1:max(1,div(blockSize,2)):length(indices)-blockSize+1)
    endpoint = length(indices)-blockSize+1
    last(starts) == endpoint || push!(starts,endpoint)
    length(starts) >= 4 || return nothing
    τ = (0:blockSize-1) .* (t[2]-t[1])
    design = hcat(ones(blockSize),τ .- mean(τ),cos.(2π*frequency .* τ),sin.(2π*frequency .* τ))
    # Square-root Hann weights reduce leakage from other modes into the fitted
    # oscillation, especially when tip bending is dominated by a lower mode.
    taper = [sinpi((i-1)/(blockSize-1)) for i in 1:blockSize]
    fit = qr(taper .* design)
    coefficients = Matrix{ComplexF64}(undef,length(starts),2)
    background = zeros(2)
    for (i,start) in enumerate(starts)
        block = start:start+blockSize-1
        pulseMode = fit \ (taper .* response[block,:])
        baseMode = fit \ (taper .* unforced[indices[block],:])
        coefficients[i,:] = complex.(pulseMode[3,:],pulseMode[4,:])
        background .= max.(background,hypot.(baseMode[3,:],baseMode[4,:]))
    end
    amplitudes = abs.(coefficients)
    rates = ntuple(c -> log_envelope_fit(timeWindow[starts],amplitudes[:,c]),2)
    signatures = coefficients[:,1] ./ coefficients[:,2]
    signature = mean(signatures)
    sameMode = isfinite(bendingFrequency) && abs(frequency-bendingFrequency) <= options.frequency_tolerance &&
               all(s -> mode_close(s,signature,options.mode_rtol),signatures)
    signalAboveNoise = all(c -> minimum(amplitudes[:,c]) >= options.minimum_snr *
        max(options.noise_floor[c],background[c]),1:2)
    observableGrowth = all(r -> r[2] >= options.minimum_r_squared &&
        abs(r[1])* (timeWindow[last(starts)]-timeWindow[first(starts)]) >= options.minimum_log_change,rates)
    return (; frequency,bending_frequency=bendingFrequency,signature,
        growth=(rates[1][1],rates[2][1]),r_squared=(rates[1][2],rates[2][2]),
        amplitudes=Tuple(median(amplitudes;dims=1)[1,:]),background=Tuple(background),
        same_mode=sameMode,above_noise=signalAboveNoise,observable_growth=observableGrowth)
end

mode_close(a,b,rtol) = isfinite(a) && isfinite(b) && abs(a-b) <= rtol*max(abs(a),abs(b))
rate_close(a,b,o) = abs(a-b) <= max(o.rate_atol,o.rate_rtol*max(abs(a),abs(b)))

function classify_equilibrium(t,unforced,pulseLow,pulseHigh,forces,o)
    fits = [modal_window(t,pulse,unforced,window,o)
            for pulse in (pulseLow,pulseHigh),window in o.windows]
    any(isnothing,fits) && return (; status=:unresolved,reason="insufficient fitting data",
        rate=NaN,frequency=NaN,signature=complex(NaN),fits)
    rates = [g for f in fits for g in f.growth]
    rate = median(rates)
    frequency = median([f.frequency for f in fits])
    signature = mean(f.signature for f in fits)
    repeatable = all(j -> all(c -> mode_close(fits[1,j].amplitudes[c]/forces[1],
        fits[2,j].amplitudes[c]/forces[2],o.scale_rtol),1:2),axes(fits,2))
    reason = if !all(f -> f.above_noise,fits)
        "modal response below signal-to-noise threshold"
    elseif !all(f -> f.same_mode && abs(f.frequency-frequency) <= o.frequency_tolerance &&
        mode_close(f.signature,signature,o.mode_rtol),fits)
        "frequency or bending/twist mode signature changed"
    elseif !repeatable
        "modal amplitudes did not scale with the pulse"
    elseif !all(f -> f.observable_growth,fits) || !all(g -> abs(g) > o.minimum_rate,rates)
        "growth sign too weak for these windows"
    elseif !all(g -> sign(g) == sign(rate) && rate_close(g,rate,o),rates)
        "growth rates disagree across channels, windows or pulses"
    else
        "accepted"
    end
    status = reason == "accepted" ? (rate > 0 ? :growth : :decay) : :unresolved
    return (; status,reason,rate,frequency,signature,fits)
end

# Return opposite observed states in traversal order. An unresolved observation
# widens the interval. A failed case breaks continuity of the evidence.
function observed_transitions(rows,states)
    brackets = NamedTuple[]
    previous = nothing
    for (i,row) in enumerate(rows)
        if row.status == :failed
            previous = nothing
        elseif row.status in states
            if !isnothing(previous) && rows[previous].status != row.status
                push!(brackets,(; first=previous,last=i,has_gap=i > previous+1))
            end
            previous = i
        end
    end
    return brackets
end

function lco_window(t,q,band,repeatTolerance)
    signal = detrend_signal(t,q)
    amplitude = (maximum(signal)-minimum(signal))/2
    frequency = signal_frequency(t,q,band)
    isfinite(frequency) || return (; amplitude,frequency,repeat_error=Inf,center=mean(q))
    dt = (last(t)-first(t))/(length(t)-1)
    scale = max(sqrt(mean(abs2,signal)),eps())
    repeatError = Inf
    repeatFrequency = frequency
    # The largest spectral peak may be a harmonic of a repeating waveform.
    for periods in 1:4
        candidate = frequency/periods
        candidate >= band[1] || break
        lag = round(Int,1/(dt*candidate))
        1 <= lag < length(t)/2 || break
        repeatError = sqrt(mean(abs2,view(signal,1:length(t)-lag) .-
            view(signal,1+lag:length(t))))/scale
        repeatFrequency = candidate
        repeatError <= repeatTolerance && break
    end
    return (; amplitude,frequency=repeatFrequency,repeat_error=repeatError,center=mean(q))
end

function classify_lco(t,response,floor,o)
    # Two disjoint tail windows, both after the plateau transient.
    tail = last(t)-2o.tail_time
    if first(t) > tail
        return (; status=:unresolved,amplitudes=(NaN,NaN),frequency=NaN)
    end
    edges = searchsortedfirst.(Ref(t),(tail,tail+o.tail_time))
    windows = (edges[1]:edges[2]-1,edges[2]:length(t))
    # Search from the minimum resolvable cycle frequency to well below Nyquist.
    band = (o.minimum_cycles/o.tail_time,0.25/(t[2]-t[1]))
    measures = [lco_window(t[ix],view(response,ix,c),band,o.repeat_rtol) for ix in windows,c in 1:2]
    amplitudes = ntuple(c -> measures[2,c].amplitude,2)
    frequency = mean(m.frequency for m in measures)
    drift = ntuple(c -> abs(measures[1,c].center-measures[2,c].center),2)
    quiet = all(c -> maximum(abs,response[first(windows[1]):end,c] .-
        mean(view(response,first(windows[1]):length(t),c))) <= floor[c],1:2)
    periodic = isfinite(frequency) && all(m -> m.repeat_error <= o.repeat_rtol &&
        abs(m.frequency-frequency) <= max(1/o.tail_time,0.05frequency),measures)
    stationary = all(c -> minimum(m.amplitude for m in measures[:,c]) > floor[c] &&
        mode_close(measures[1,c].amplitude,amplitudes[c],o.amplitude_rtol) &&
        drift[c] <= max(floor[c],o.amplitude_rtol*amplitudes[c]),1:2)
    status = quiet ? :steady : periodic && stationary ? :lco : :unresolved
    return (; status,amplitudes,frequency=status == :lco ? frequency : NaN)
end
