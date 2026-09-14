# Shared moving-block damping analysis, extracted from the archived driver.
using DelimitedFiles, Statistics, FFTW
function read_response_history(history_path)
    # Only the columns needed by the convergence metric are loaded. The full
    # production history remains on disk for later inspection.
    raw, header_matrix = readdlm(history_path, ',', header = true)
    headers = String.(vec(header_matrix))
    column(name) = begin
        index = findfirst(==(name), headers)
        isnothing(index) && error("Column '$name' is missing from $history_path")
        Float64.(raw[:, index])
    end
    return (;
        time = column("time_s"),
        pitch = column("propeller_1_pitch_deg"),
        yaw = column("propeller_1_yaw_deg"),
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
    for index in 2:(length(signal) - 1)
        if signal[index] > signal[index - 1] && signal[index] >= signal[index + 1]
            push!(peak_indices, index)
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
function unavailable_metrics(peak_count = 0; block_size = 0, block_count = 0)
    # Return one consistent record shape even when damping cannot be measured.
    # This keeps the CSV writer and comparison logic simple and preserves the
    # reason through the surrounding case status.
    return (;
        available = false,
        valid_for_convergence = false,
        peak_count,
        moving_block_lambda_per_s = NaN,
        growth_rate_per_s = NaN,
        damping_ratio_percent = NaN,
        frequency_hz = NaN,
        fit_r_squared = NaN,
        envelope_ratio = NaN,
        block_size,
        block_count,
        analyzed_start_s = NaN,
        analyzed_end_s = NaN,
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
    block_duration_s = nothing,
    block_overlap = nothing,
    frequency_min_hz = 0.0,
    frequency_max_hz = Inf,
    apply_hann_window = false,
)
    length(time) == length(signal) || error("Time and signal vectors must have equal length")

    # Step 1: discard startup/trim data outside the user-selected fit window.
    fit_indices = findall(index -> fit_start_s <= time[index] <= fit_end_s, eachindex(time))
    length(fit_indices) >= 3 || return unavailable_metrics()
    fit_time = time[fit_indices]
    fit_signal = signal[fit_indices]

    # Step 2: peaks below one percent of the maximum absolute response are
    # ignored as numerical noise.
    peak_minimum = maximum(abs, fit_signal) / 100.0
    all_peak_indices = positive_peak_indices(fit_signal)
    peak_indices = filter(index -> fit_signal[index] >= peak_minimum, all_peak_indices)
    required_peaks = max(minimum_peaks, peak_from_start + peak_from_end + 1)
    length(peak_indices) >= required_peaks || return unavailable_metrics(length(peak_indices))

    # Step 3: cut from the requested peak near the start to the requested peak
    # counted back from the end. Defaults retain first through last valid peak.
    start_index = peak_indices[peak_from_start]
    end_index = peak_indices[end - peak_from_end]
    analyzed_time = fit_time[start_index:end_index]
    analyzed_signal = fit_signal[start_index:end_index]
    sample_size = length(analyzed_time)

    # Step 4: choose the FFT block. The original method uses a fraction of the
    # sample record. Corrected studies can instead request a physical duration,
    # which preserves spectral resolution when dt changes.
    sample_intervals = diff(analyzed_time)
    isempty(sample_intervals) && return unavailable_metrics(length(peak_indices))
    sample_dt = median(sample_intervals)
    sample_dt > 0 || return unavailable_metrics(length(peak_indices))
    if isnothing(block_duration_s)
        block_size = initial_block_size
        size_ratio = block_size / sample_size
        while size_ratio < size_ratio_lower
            block_size *= 2
            size_ratio = block_size / sample_size
        end
        while size_ratio > size_ratio_upper
            block_size = div(block_size, 2)
            size_ratio = block_size / sample_size
        end
    else
        block_duration_s > 0 || error("Moving-block duration must be positive")
        target_samples = max(2.0, block_duration_s / sample_dt)
        # FFTW supports arbitrary lengths. Power-of-two rounding would change
        # the physical fit window when time steps are not refined by factors
        # of two, contaminating a time-step convergence comparison.
        block_size = round(Int, target_samples)
    end
    block_size >= 1 || return unavailable_metrics(length(peak_indices))
    block_size <= sample_size || return unavailable_metrics(
        length(peak_indices);
        block_size,
    )

    0.0 <= frequency_min_hz < frequency_max_hz || error(
        "Moving-block frequency bounds must satisfy 0 <= min < max",
    )
    stride = isnothing(block_overlap) ? 1 : begin
        0.0 <= block_overlap < 1.0 || error("Moving-block overlap must be in [0,1)")
        max(1, round(Int, block_size * (1.0 - block_overlap)))
    end
    block_starts = collect(1:stride:(sample_size - block_size + 1))
    last_start = sample_size - block_size + 1
    last(block_starts) == last_start || push!(block_starts, last_start)
    block_count = length(block_starts)
    log_block_amplitude = Vector{Float64}(undef, block_count)
    tracked_frequency = Vector{Float64}(undef, block_count)
    hann_window = block_size == 1 ? ones(1) :
        0.5 .- 0.5 .* cos.(2pi .* (0:(block_size - 1)) ./ (block_size - 1))
    frequency_bins = collect(0:div(block_size, 2)) ./ (block_size * sample_dt)
    frequency_indices = findall(
        frequency -> frequency_min_hz <= frequency <= frequency_max_hz,
        frequency_bins,
    )
    isempty(frequency_indices) && return unavailable_metrics(
        length(peak_indices);
        block_size,
        block_count,
    )
    for (block_index, block_start) in enumerate(block_starts)
        block = collect(@view analyzed_signal[block_start:block_start+block_size-1])
        if apply_hann_window
            block .-= mean(block)
            block .*= hann_window
        end
        spectrum = abs.(fft(block)) ./ block_size
        half_index = div(block_size, 2) + 1
        single_sided = spectrum[1:half_index]
        length(single_sided) > 2 && (single_sided[2:end-1] .*= 2.0)
        local_index = argmax(@view single_sided[frequency_indices])
        spectral_index = frequency_indices[local_index]
        maximum_amplitude = single_sided[spectral_index]
        maximum_amplitude > 0 || return unavailable_metrics(
            length(peak_indices);
            block_size,
            block_count,
        )
        log_block_amplitude[block_index] = log(maximum_amplitude)
        tracked_frequency[block_index] = frequency_bins[spectral_index]
    end

    # Step 6: least-squares fit
    #
    #     log(Xbar) = intercept + lambda * block_start_time
    #
    # lambda < 0 indicates decay; lambda > 0 indicates growth.
    block_time = analyzed_time[block_starts]
    mean_time = sum(block_time) / block_count
    mean_log_amplitude = sum(log_block_amplitude) / block_count
    centered_time = block_time .- mean_time
    denominator = sum(abs2, centered_time)
    denominator > eps(Float64) || return unavailable_metrics(
        length(peak_indices);
        block_size,
        block_count,
    )

    moving_block_lambda = sum(
        centered_time .* (log_block_amplitude .- mean_log_amplitude),
    ) /
        denominator
    fitted_log_amplitude = mean_log_amplitude .+ moving_block_lambda .* centered_time
    residual_sum_squares = sum(abs2, log_block_amplitude .- fitted_log_amplitude)
    total_sum_squares = sum(abs2, log_block_amplitude .- mean_log_amplitude)
    fit_r_squared = total_sum_squares > eps(Float64) ?
        1.0 - residual_sum_squares / total_sum_squares : 1.0

    # Step 7: use the frequency of the same in-band spectral component whose
    # amplitude was fitted. This prevents the damping ratio from mixing one
    # mode's decay rate with a different mode's peak spacing.
    frequency_hz = median(tracked_frequency)
    damped_omega = 2pi * frequency_hz
    damping_ratio_percent = -100.0 * moving_block_lambda /
        hypot(moving_block_lambda, damped_omega)
    envelope_ratio = exp(moving_block_lambda * (block_time[end] - block_time[1]))

    return (;
        available = true,
        valid_for_convergence = fit_r_squared >= minimum_fit_r_squared,
        peak_count = length(peak_indices),
        moving_block_lambda_per_s = moving_block_lambda,
        growth_rate_per_s = moving_block_lambda,
        damping_ratio_percent,
        frequency_hz,
        fit_r_squared,
        envelope_ratio,
        block_size,
        block_count,
        analyzed_start_s = analyzed_time[1],
        analyzed_end_s = analyzed_time[end],
    )
end

"""
    completed_requested_window(time,requested_end_time)

Checks whether a child history reached its requested end time

# Arguments
- `time`: time vector written by the child simulation
- `requested_end_time`: requested simulation end time in s
"""
function completed_requested_window(time, requested_end_time)
    isempty(time) && return false
    length(time) == 1 && return time[end] >= requested_end_time - eps(requested_end_time)
    final_step = time[end] - time[end - 1]
    return time[end] >= requested_end_time - 1.01 * final_step
end

