# Shared moving-block damping analysis, extracted from the archived driver.
using DelimitedFiles,Statistics,FFTW
function read_response_history(history_path)
    # Only the columns needed by the convergence metric are loaded. The full
    # production history remains on disk for later inspection.
    raw,header_matrix = readdlm(history_path,',',header=true)
    headers = String.(vec(header_matrix))
    column(name) = begin
        index = findfirst(==(name),headers)
        isnothing(index) && error("Column '$name' is missing from $history_path")
        Float64.(raw[:,index])
    end
    return (;
        time=column("time_s"),
        pitch=column("propeller_1_pitch_deg"),
        yaw=column("propeller_1_yaw_deg")
    )
end

"""
    positive_peak_indices(signal)

Returns indices of positive-direction local maxima

# Arguments
- `signal`: sampled response vector
"""
function positive_peak_indices(signal)
    # This local-maximum operation keeps the sweep independent of Peaks.jl.
    # Positive maxima define the beginning and end of the response segment.
    peak_indices = Int[]
    for index = 2:(length(signal)-1)
        if signal[index] > signal[index-1] && signal[index] >= signal[index+1]
            push!(peak_indices,index)
        end
    end
    return peak_indices
end

"""
    unavailable_metrics(peak_count=0; block_size=0,block_count=0)

Creates a damping-result record for a signal that could not be evaluated

# Arguments
- `peak_count`: number of valid peaks found before the analysis failed

# Keyword arguments
- `block_size`: final moving-block size, when available
- `block_count`: number of overlapping blocks, when available
"""
function unavailable_metrics(peak_count=0; block_size=0,block_count=0)
    # Return one consistent record shape even when damping cannot be measured.
    # This keeps the CSV writer and comparison logic simple and preserves the
    # reason through the surrounding case status.
    return (;
        available=false,
        valid_for_convergence=false,
        peak_count,
        moving_block_lambda_per_s=NaN,
        growth_rate_per_s=NaN,
        damping_ratio_percent=NaN,
        frequency_hz=NaN,
        fit_r_squared=NaN,
        envelope_ratio=NaN,
        block_size,
        block_count,
        analyzed_start_s=NaN,
        analyzed_end_s=NaN
    )
end

"""
    moving_block_metrics(time,signal; kwargs...)

Apply the locally implemented moving-block FFT damping method to one response channel.

# Arguments
- `time`: sampled time vector in s
- `signal`: sampled pitch or yaw response vector

# Keyword arguments
- `fit_start_s`: beginning of the permitted damping-fit interval in s
- `fit_end_s`: end of the permitted damping-fit interval in s
- `minimum_peaks`: minimum number of retained positive peaks
- `minimum_fit_r_squared`: minimum accepted linear-fit quality
- `initial_block_size`: initial power-of-two FFT block length in samples
- `size_ratio_lower`: lower block-to-record size ratio
- `size_ratio_upper`: upper block-to-record size ratio
- `peak_from_start`: first positive peak retained from the beginning
- `peak_from_end`: number of positive peaks discarded from the end
- `block_duration_s`: optional target FFT-block duration; when supplied it
  replaces the sample-ratio rule and preserves frequency resolution across dt
- `block_overlap`: fractional overlap between consecutive FFT blocks
- `frequency_min_hz`, `frequency_max_hz`: tracked modal-frequency band
- `apply_hann_window`: remove the block mean and apply a Hann window
"""
function moving_block_metrics(
    time,
    signal;
    fit_start_s,
    fit_end_s,
    minimum_peaks,
    minimum_fit_r_squared,
    initial_block_size,
    size_ratio_lower,
    size_ratio_upper,
    peak_from_start,
    peak_from_end,
    block_duration_s=nothing,
    block_overlap=nothing,
    frequency_min_hz=0.0,
    frequency_max_hz=Inf,
    apply_hann_window=false
)
    length(time) == length(signal) || error("Time and signal vectors must have equal length")

    # Step 1: discard startup/trim data outside the user-selected fit window.
    fitIndices = findall(index -> fit_start_s <= time[index] <= fit_end_s,eachindex(time))
    length(fitIndices) >= 3 || return unavailable_metrics()
    fitTime = time[fitIndices]
    fitSignal = signal[fitIndices]

    # Step 2: peaks below one percent of the maximum absolute response are
    # ignored as numerical noise.
    peakMinimum = maximum(abs,fitSignal) / 100.0
    allPeakIndices = positive_peak_indices(fitSignal)
    peak_indices = filter(index -> fitSignal[index] >= peakMinimum,allPeakIndices)
    requiredPeaks = max(minimum_peaks,peak_from_start + peak_from_end + 1)
    length(peak_indices) >= requiredPeaks || return unavailable_metrics(length(peak_indices))

    # Step 3: cut from the requested peak near the start to the requested peak
    # counted back from the end. Defaults retain first through last valid peak.
    startIndex = peak_indices[peak_from_start]
    endIndex = peak_indices[end-peak_from_end]
    analyzedTime = fitTime[startIndex:endIndex]
    analyzedSignal = fitSignal[startIndex:endIndex]
    sampleSize = length(analyzedTime)

    # Step 4: choose the FFT block. The original method uses a fraction of the
    # sample record. Corrected studies can instead request a physical duration,
    # which preserves spectral resolution when dt changes.
    sampleIntervals = diff(analyzedTime)
    isempty(sampleIntervals) && return unavailable_metrics(length(peak_indices))
    sampleDt = median(sampleIntervals)
    sampleDt > 0 || return unavailable_metrics(length(peak_indices))
    if isnothing(block_duration_s)
        block_size = initial_block_size
        sizeRatio = block_size / sampleSize
        while sizeRatio < size_ratio_lower
            block_size *= 2
            sizeRatio = block_size / sampleSize
        end
        while sizeRatio > size_ratio_upper
            block_size = div(block_size,2)
            sizeRatio = block_size / sampleSize
        end
    else
        block_duration_s > 0 || error("Moving-block duration must be positive")
        targetSamples = max(2.0,block_duration_s / sampleDt)
        # FFTW supports arbitrary lengths. Power-of-two rounding would change
        # the physical fit window when time steps are not refined by factors
        # of two, contaminating a time-step convergence comparison.
        block_size = round(Int,targetSamples)
    end
    block_size >= 1 || return unavailable_metrics(length(peak_indices))
    block_size <= sampleSize || return unavailable_metrics(
        length(peak_indices);
        block_size
    )

    0.0 <= frequency_min_hz < frequency_max_hz || error(
        "Moving-block frequency bounds must satisfy 0 <= min < max"
    )
    stride =
        isnothing(block_overlap) ? 1 : begin
            0.0 <= block_overlap < 1.0 || error("Moving-block overlap must be in [0,1)")
            max(1,round(Int,block_size * (1.0 - block_overlap)))
        end
    blockStarts = collect(1:stride:(sampleSize-block_size+1))
    lastStart = sampleSize - block_size + 1
    last(blockStarts) == lastStart || push!(blockStarts,lastStart)
    block_count = length(blockStarts)
    logBlockAmplitude = Vector{Float64}(undef,block_count)
    trackedFrequency = Vector{Float64}(undef,block_count)
    hann_window = block_size == 1 ? ones(1) :
                  0.5 .- 0.5 .* cos.(2pi .* (0:(block_size-1)) ./ (block_size - 1))
    frequencyBins = collect(0:div(block_size,2)) ./ (block_size * sampleDt)
    frequencyIndices = findall(
        frequency -> frequency_min_hz <= frequency <= frequency_max_hz,
        frequencyBins
    )
    isempty(frequencyIndices) && return unavailable_metrics(
        length(peak_indices);
        block_size,
        block_count
    )
    for (blockIndex,blockStart) in enumerate(blockStarts)
        block = collect(@view analyzedSignal[blockStart:blockStart+block_size-1])
        if apply_hann_window
            block .-= mean(block)
            block .*= hann_window
        end
        spectrum = abs.(fft(block)) ./ block_size
        halfIndex = div(block_size,2) + 1
        singleSided = spectrum[1:halfIndex]
        length(singleSided) > 2 && (singleSided[2:end-1] .*= 2.0)
        localIndex = argmax(@view singleSided[frequencyIndices])
        spectralIndex = frequencyIndices[localIndex]
        maximumAmplitude = singleSided[spectralIndex]
        maximumAmplitude > 0 || return unavailable_metrics(
            length(peak_indices);
            block_size,
            block_count
        )
        logBlockAmplitude[blockIndex] = log(maximumAmplitude)
        trackedFrequency[blockIndex] = frequencyBins[spectralIndex]
    end

    # Step 6: least-squares fit
    #
    #     log(Xbar) = intercept + lambda * block_start_time
    #
    # lambda < 0 indicates decay; lambda > 0 indicates growth.
    blockTime = analyzedTime[blockStarts]
    meanTime = sum(blockTime) / block_count
    meanLogAmplitude = sum(logBlockAmplitude) / block_count
    centeredTime = blockTime .- meanTime
    denominator = sum(abs2,centeredTime)
    denominator > eps(Float64) || return unavailable_metrics(
        length(peak_indices);
        block_size,
        block_count
    )

    movingBlockLambda = sum(
        centeredTime .* (logBlockAmplitude .- meanLogAmplitude)
    ) /
                        denominator
    fitted_log_amplitude = meanLogAmplitude .+ movingBlockLambda .* centeredTime
    residualSumSquares = sum(abs2,logBlockAmplitude .- fitted_log_amplitude)
    totalSumSquares = sum(abs2,logBlockAmplitude .- meanLogAmplitude)
    fit_r_squared = totalSumSquares > eps(Float64) ?
                    1.0 - residualSumSquares / totalSumSquares : 1.0

    # Step 7: use the frequency of the same in-band spectral component whose
    # amplitude was fitted. This prevents the damping ratio from mixing one
    # mode's decay rate with a different mode's peak spacing.
    frequency_hz = median(trackedFrequency)
    dampedOmega = 2pi * frequency_hz
    damping_ratio_percent = -100.0 * movingBlockLambda /
                            hypot(movingBlockLambda,dampedOmega)
    envelope_ratio = exp(movingBlockLambda * (blockTime[end] - blockTime[1]))

    return (;
        available=true,
        valid_for_convergence=fit_r_squared >= minimum_fit_r_squared,
        peak_count=length(peak_indices),
        moving_block_lambda_per_s=movingBlockLambda,
        growth_rate_per_s=movingBlockLambda,
        damping_ratio_percent,
        frequency_hz,
        fit_r_squared,
        envelope_ratio,
        block_size,
        block_count,
        analyzed_start_s=analyzedTime[1],
        analyzed_end_s=analyzedTime[end]
    )
end

"""
    completed_requested_window(time,requested_end_time)

Checks whether a child history reached its requested end time

# Arguments
- `time`: time vector written by the child simulation
- `requested_end_time`: requested simulation end time in s
"""
function completed_requested_window(time,requested_end_time)
    isempty(time) && return false
    length(time) == 1 && return time[end] >= requested_end_time - eps(requested_end_time)
    finalStep = time[end] - time[end-1]
    return time[end] >= requested_end_time - 1.01 * finalStep
end

