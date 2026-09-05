# Reusable cases, damping metrics, process orchestration, and reporting for the
# legacy one-stage aeroelastic convergence study.

"""
Run a one-parameter-at-a-time UVLM convergence study for the Chang
wing--propeller aeroelastic model.

The study varies:

  1. wing aerodynamic/structural mesh,
  2. propeller aerodynamic mesh,
  3. retained wake length,
  4. vortex finite-core factor, and
  5. azimuthal time step.

Each case runs the production entry point in a fresh Julia process, so all
structural, aerodynamic, and coupling routines are preserved.
The post-impulse pitch and yaw responses are evaluated with the moving-block
FFT method implemented below. For every overlapping signal block, the maximum
single-sided FFT magnitude is retained and fitted as

    log(Xbar(t)) = a + lambda * t.

The moving-block slope `lambda` (1/s) is the primary convergence quantity. It
is also converted to a conventional dimensionless damping ratio using

    zeta = -lambda / sqrt(lambda^2 + (2*pi*f_d)^2).

Within each parameter family, the last (finest) configured level is the
reference. A level is considered converged when the pitch and yaw
moving-block slopes both differ from that reference by no more than
`MOVING_BLOCK_LAMBDA_TOLERANCE_PER_S`.

Run from the repository root with

    julia --project=lib/WingPropellerUVLM lib/WingPropellerUVLM/examples/chang_linear_aeroelastic/studies/convergence/run_chang_uvlm_convergence_sweep.jl

Useful PowerShell examples:

    # Print the generated case matrix without running simulations.
    \$env:CHANG_CONVERGENCE_DRY_RUN = "true"

    # Run only the wake and time-step studies.
    \$env:CHANG_CONVERGENCE_FAMILIES = "wake_length,time_step"

    # Use a longer damping-fit window.
    \$env:CHANG_CONVERGENCE_END_TIME_S = "5.0"

The edit-friendly parameter tables are immediately below.
"""

using Dates
using DelimitedFiles
using FFTW
using Printf
using Statistics
using Plots

# HOW THIS DRIVER REACHES THE PRODUCTION SOLVER
# ---------------------------------------------
# This file is the parent/orchestrator; it does not assemble or time-march the
# aeroelastic model itself. For every unique case it:
#
#   1. creates a UVLMConvergenceCase below;
#   2. converts that case into production CHANG_* environment overrides;
#   3. launches run_chang_linear_aeroelastic.jl in a fresh process;
#   4. reads the resulting pitch/yaw history CSV;
#   5. calculates moving-block damping; and
#   6. compares every refinement level with the last level of its family.
#
# A fresh process avoids constant/type redefinition problems between cases and
# prevents one UVLM wake state from leaking into the next simulation.

# ### Study definition

# These nominal values define the retained one-stage convergence study and are
# intentionally independent of the lighter interactive defaults in
# chang_case.jl. Wake row counts are generated from wake age and azimuth step.
const NOMINAL_WING_SPAN_PANELS = 30
const NOMINAL_WING_CHORD_PANELS = 5
const NOMINAL_PROP_RADIAL_PANELS = 10
const NOMINAL_PROP_CHORD_PANELS = 10
const NOMINAL_WAKE_REVOLUTIONS = 1.0
const NOMINAL_FCORE_SEGMENT_FACTOR = 0.25
const NOMINAL_AZIMUTH_STEP_DEG = 5.0

# The study is one-parameter-at-a-time (OAT): only the values listed in the
# active family change; every other parameter remains at the NOMINAL_* value.
# Levels must be ordered from coarse to fine. The last entry in each table is
# used as that family's numerical reference. Therefore, do not put the nominal
# value last unless it is genuinely the finest value you want as reference.

# Wing mesh is listed as span x chord. In the current Chang implementation the
# spanwise wing count also sets the structural beam mesh, so this family is a
# coupled aero-structural discretization check rather than a pure UVLM check.
const WING_MESH_LEVELS = [
    (label = "coarse_20x10", span = 20, chord = 10),
    (label = "nominal_30x10", span = 30, chord = 10),
    (label = "fine_40x10", span = 40, chord = 10),
]

# Propeller mesh is radial x chord for each blade.
const PROPELLER_MESH_LEVELS = [
    (label = "coarse_10x10", radial = 10, chord = 10),
    (label = "nominal_15x15", radial = 15, chord = 15),
    (label = "fine_20x20", radial = 20, chord = 20),
]

# Wake length is specified as physical rotor revolutions, not as row count.
# wake_rows() converts it to rows for the selected azimuth step.
const WAKE_LENGTH_LEVELS_REVOLUTIONS = [1.0, 2.0, 3.0]

# The production core law is epsilon = factor * local 3-D vortex-segment length. Smaller
# factors are ordered later and treated as the refinement direction here.
const FINITE_CORE_LEVELS = [0.25, 0.1, 0.05]

# Because the time step is derived from rotor azimuth, smaller angles are
# finer. Wake rows are rescaled for every entry to retain the same wake age.
const AZIMUTH_STEP_LEVELS_DEG = [5.0, 2.5]

# Moving-block settings used by the local implementation. The FFT block must be a
# power of two. It is doubled/halved until it occupies 25--50% of the selected
# response segment. The lambda tolerance is an absolute difference in 1/s;
# this remains meaningful close to flutter, where relative errors do not.
const MOVING_BLOCK_INITIAL_SIZE = 512
const MOVING_BLOCK_SIZE_RATIO_LOWER = 0.25
const MOVING_BLOCK_SIZE_RATIO_UPPER = 0.50
const MOVING_BLOCK_PEAK_FROM_START = 1
const MOVING_BLOCK_PEAK_FROM_END = 0
const MOVING_BLOCK_LAMBDA_TOLERANCE_PER_S = 0.01

# ### Runtime controls and case generation

# Directory containing this study, the example root, the package project, and
# the production runner used by each isolated child process.
const CONVERGENCE_STUDY_DIR = @__DIR__
const CONVERGENCE_EXAMPLE_DIR = normpath(joinpath(CONVERGENCE_STUDY_DIR, "..", ".."))
const CONVERGENCE_PACKAGE_DIR = normpath(joinpath(CONVERGENCE_EXAMPLE_DIR, "..", ".."))
const CONVERGENCE_CASE_RUNNER = joinpath(CONVERGENCE_EXAMPLE_DIR, "run_chang_linear_aeroelastic.jl")

# These helpers make environment variables optional: when a variable is not
# defined, the provided Julia default is used. This permits both editing the
# constants above and overriding run-level controls from PowerShell.
"""
    convergence_float(name,default)

Reads a floating-point environment variable, or returns the default value when absent

# Arguments
- `name`: environment-variable name
- `default`: value used when the variable is absent
"""
convergence_float(name, default) = parse(Float64, get(ENV, name, string(default)))

"""
    convergence_int(name,default)

Reads an integer environment variable, or returns the default value when absent

# Arguments
- `name`: environment-variable name
- `default`: value used when the variable is absent
"""
convergence_int(name, default) = parse(Int, get(ENV, name, string(default)))

"""
    convergence_bool(name,default)

Reads a Boolean environment variable. Accepted values are `true/false`, `yes/no`, `on/off`, and `1/0`, without case sensitivity

# Arguments
- `name`: environment-variable name
- `default`: value used when the variable is absent
"""
function convergence_bool(name, default)
    value = lowercase(strip(get(ENV, name, string(default))))
    value in ("1", "true", "yes", "on") && return true
    value in ("0", "false", "no", "off") && return false
    error("$name must be true/false, yes/no, on/off, or 1/0")
end

"""
    number_token(value; digits=3)

Converts a number to a deterministic file-name token. Decimal points become `p` and minus signs become `m`

# Arguments
- `value`: number to convert

# Keyword arguments
- `digits`: number of decimal digits in the token
"""
function number_token(value; digits = 3)
    # Decimal points and minus signs are converted so labels are safe in file
    # names, for example 2.5 -> "2p500" and -1.0 -> "m1p000".
    text = @sprintf("%.*f", digits, value)
    return replace(text, "-" => "m", "." => "p")
end

"""
    UVLMConvergenceCase composite type

Defines the numerical settings and result identifiers of one convergence case

# Fields
- `family::Symbol`: parameter family used to group convergence results
- `level::Int`: coarse-to-fine position inside the family
- `label::String`: readable case and file label
- `wing_span_panels::Int`: wing spanwise panel count
- `wing_chord_panels::Int`: wing chordwise panel count
- `prop_radial_panels::Int`: propeller radial panel count per blade
- `prop_chord_panels::Int`: propeller chordwise panel count per blade
- `wake_revolutions::Float64`: retained wake age in rotor revolutions
- `fcore_segment_factor::Float64`: finite-core radius divided by 3-D segment length
- `azimuth_step_deg::Float64`: rotor azimuth increment per time step
"""
Base.@kwdef struct UVLMConvergenceCase
    # family/level/label identify where the result belongs in the convergence
    # plot. All remaining fields fully identify the numerical UVLM setup.
    family::Symbol
    level::Int
    label::String
    wing_span_panels::Int = NOMINAL_WING_SPAN_PANELS
    wing_chord_panels::Int = NOMINAL_WING_CHORD_PANELS
    prop_radial_panels::Int = NOMINAL_PROP_RADIAL_PANELS
    prop_chord_panels::Int = NOMINAL_PROP_CHORD_PANELS
    wake_revolutions::Float64 = NOMINAL_WAKE_REVOLUTIONS
    fcore_segment_factor::Float64 = NOMINAL_FCORE_SEGMENT_FACTOR
    fcore_chord_factor::Float64 = 0.0
    azimuth_step_deg::Float64 = NOMINAL_AZIMUTH_STEP_DEG
end

"""
    validate_case(case)

Checks that all mesh counts and continuous UVLM controls are valid and returns the original case

# Arguments
- `case::UVLMConvergenceCase`: case to validate
"""
function validate_case(case::UVLMConvergenceCase)
    case.wing_span_panels > 0 || error("Wing spanwise panels must be positive")
    case.wing_chord_panels > 0 || error("Wing chordwise panels must be positive")
    case.prop_radial_panels > 0 || error("Propeller radial panels must be positive")
    case.prop_chord_panels > 0 || error("Propeller chordwise panels must be positive")
    case.wake_revolutions > 0 || error("Wake length must be positive")
    all(x -> isfinite(x) && x >= 0, (case.fcore_segment_factor, case.fcore_chord_factor)) ||
        error("Finite-core factors must be finite and nonnegative")
    max(case.fcore_segment_factor, case.fcore_chord_factor) > 0 || error("Free wakes require a positive core")
    case.azimuth_step_deg > 0 || error("Azimuth step must be positive")
    return case
end

"""
    wake_rows(case)

Converts retained wake age from rotor revolutions to shed wake rows

# Arguments
- `case::UVLMConvergenceCase`: case containing wake age and azimuth step
"""
function wake_rows(case::UVLMConvergenceCase)
    # One wake row is shed per time step. ceil prevents the retained wake from
    # being shorter than the requested physical age.
    #
    # rows = (360 deg/revolution * requested revolutions) / deg per step
    return max(1, ceil(Int, 360.0 * case.wake_revolutions / case.azimuth_step_deg))
end

"""
    physical_time_step(case,speed_mps,trim_rpm,trim_speed_mps,radius_m)

Computes the physical time step associated with the case's azimuth increment while preserving the Chang advance ratio

# Arguments
- `case::UVLMConvergenceCase`: case containing the azimuth step
- `speed_mps`: case flow speed in m/s
- `trim_rpm`: propeller speed at the trim operating point in RPM
- `trim_speed_mps`: flow speed at the trim operating point in m/s
- `radius_m`: propeller radius in m
"""
function physical_time_step(case, speed_mps, trim_rpm, trim_speed_mps, radius_m)
    # The Chang model preserves the advance ratio when flow speed changes.
    # First recover that reference advance ratio, then find the case rotation
    # rate and convert an azimuth increment into seconds. This value is written
    # to the summary for interpretation; the production driver independently
    # calculates the same dt for its actual integration.
    trim_omega = trim_rpm * 2pi / 60.0
    advance_ratio_mu = trim_speed_mps / (trim_omega * radius_m)
    case_omega = speed_mps / (advance_ratio_mu * radius_m)
    return deg2rad(case.azimuth_step_deg) / case_omega
end

"""
    case_key(case)

Returns the numerical settings that uniquely identify a simulation

# Arguments
- `case::UVLMConvergenceCase`: case to identify
"""
function case_key(case::UVLMConvergenceCase)
    # family and label are deliberately excluded. The nominal case appears in
    # several families but has the same numerical inputs, so main() can run it
    # once and reuse its history and damping result.
    return (
        case.wing_span_panels,
        case.wing_chord_panels,
        case.prop_radial_panels,
        case.prop_chord_panels,
        case.wake_revolutions,
        case.fcore_segment_factor,
        case.fcore_chord_factor,
        case.azimuth_step_deg,
    )
end

"""
    requested_families()

Parses the parameter families selected through `CHANG_CONVERGENCE_FAMILIES`
"""
function requested_families()
    # Example: CHANG_CONVERGENCE_FAMILIES="wake_length,time_step". The alias
    # "mesh" expands to both wing_mesh and propeller_mesh.
    raw = lowercase(strip(get(
        ENV,
        "CHANG_CONVERGENCE_FAMILIES",
        "wing_mesh,propeller_mesh,wake_length,finite_core,time_step",
    )))
    tokens = filter(!isempty, strip.(split(raw, ',')))
    "mesh" in tokens && append!(tokens, ["wing_mesh", "propeller_mesh"])
    families = unique(Symbol.(filter(!=("mesh"), tokens)))
    allowed = Set((:wing_mesh, :propeller_mesh, :wake_length, :finite_core, :time_step))
    invalid = filter(family -> !(family in allowed), families)
    isempty(invalid) || error(
        "Unknown CHANG_CONVERGENCE_FAMILIES value(s): $(join(string.(invalid), ", "))",
    )
    isempty(families) && error("At least one convergence family is required")
    return families
end

"""
    build_cases(families)

Expands selected families into ordered, one-parameter-at-a-time convergence cases

# Arguments
- `families`: vector of requested parameter-family symbols
"""
function build_cases(families)
    # Each block changes only its own parameter and relies on the struct
    # defaults for every other setting. This is what enforces the OAT design.
    cases = UVLMConvergenceCase[]

    if :wing_mesh in families
        for (level, mesh) in enumerate(WING_MESH_LEVELS)
            push!(cases, validate_case(UVLMConvergenceCase(
                family = :wing_mesh,
                level = level,
                label = mesh.label,
                wing_span_panels = mesh.span,
                wing_chord_panels = mesh.chord,
            )))
        end
    end

    if :propeller_mesh in families
        for (level, mesh) in enumerate(PROPELLER_MESH_LEVELS)
            push!(cases, validate_case(UVLMConvergenceCase(
                family = :propeller_mesh,
                level = level,
                label = mesh.label,
                prop_radial_panels = mesh.radial,
                prop_chord_panels = mesh.chord,
            )))
        end
    end

    if :wake_length in families
        for (level, revolutions) in enumerate(WAKE_LENGTH_LEVELS_REVOLUTIONS)
            push!(cases, validate_case(UVLMConvergenceCase(
                family = :wake_length,
                level = level,
                label = "wake_$(number_token(revolutions))_rev",
                wake_revolutions = revolutions,
            )))
        end
    end

    if :finite_core in families
        for (level, factor) in enumerate(FINITE_CORE_LEVELS)
            push!(cases, validate_case(UVLMConvergenceCase(
                family = :finite_core,
                level = level,
                label = "core_$(number_token(factor))_ds",
                fcore_segment_factor = factor,
            )))
        end
    end

    if :time_step in families
        for (level, step_deg) in enumerate(AZIMUTH_STEP_LEVELS_DEG)
            push!(cases, validate_case(UVLMConvergenceCase(
                family = :time_step,
                level = level,
                label = "azimuth_$(number_token(step_deg))_deg",
                azimuth_step_deg = step_deg,
            )))
        end
    end

    return cases
end

# ### Damping extraction

"""
    read_response_history(history_path)

Reads time, propeller pitch, and propeller yaw from a production history CSV

# Arguments
- `history_path`: path to the production history CSV
"""
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
        block_size = 2^round(Int, log2(target_samples))
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

# ### Child-process execution

"""
    child_environment(case,output_directory,output_label; kwargs...)

Creates the environment dictionary consumed by the production Chang driver.

# Arguments
- `case::UVLMConvergenceCase`: numerical settings for the child simulation
- `output_directory`: directory for child output files
- `output_label`: common prefix for the child output files

# Keyword arguments
- `speed_mps`: flow speed in m/s
- `trim_rpm`: propeller RPM at the trim reference speed
- `trim_speed_mps`: trim reference speed in m/s
- `end_time_s`: requested simulation end time in s
- `hard_angle_deg`: propeller-angle safety-abort limit in degrees
"""
function child_environment(
    case,
    output_directory,
    output_label;
    speed_mps,
    trim_rpm,
    trim_speed_mps,
    end_time_s,
    hard_angle_deg,
    angle_of_attack_deg = nothing,
    sideslip_deg = nothing,
    interaction_on = nothing,
    ga_rho_inf = nothing,
    coupling_relaxation = nothing,
)
    # Direct CHANG_* overrides keep the production case and every sweep on one
    # configuration path. copy(ENV) preserves Julia/package settings.
    environment = copy(ENV)
    rows = wake_rows(case)

    # Physical operating point and duration.
    environment["CHANG_SPEED_MPS"] = string(speed_mps)
    environment["CHANG_AOA_DEG"] = isnothing(angle_of_attack_deg) ?
        get(ENV, "CHANG_CONVERGENCE_AOA_DEG", "0.0") : string(angle_of_attack_deg)
    !isnothing(sideslip_deg) && (environment["CHANG_SIDESLIP_DEG"] = string(sideslip_deg))
    environment["CHANG_ROTATION_RPM"] = string(trim_rpm)
    environment["CHANG_TRIM_SPEED_MPS"] = string(trim_speed_mps)
    environment["CHANG_END_TIME_S"] = string(end_time_s)

    # Aerodynamic grids and azimuth-based time step.
    environment["CHANG_WING_SPAN_PANELS"] = string(case.wing_span_panels)
    environment["CHANG_WING_CHORD_PANELS"] = string(case.wing_chord_panels)
    environment["CHANG_PROP_RADIAL_PANELS"] = string(case.prop_radial_panels)
    environment["CHANG_PROP_CHORD_PANELS"] = string(case.prop_chord_panels)
    environment["CHANG_AZIMUTH_STEP_DEG"] = string(case.azimuth_step_deg)

    # Wake age and Biot--Savart finite-core regularization. Wing and propeller
    # retain the same physical wake age even when the azimuth step changes.
    environment["CHANG_WAKE_ROWS_WING"] = string(rows)
    environment["CHANG_WAKE_ROWS_PROPELLER"] = string(rows)
    environment["CHANG_FCORE_SEGMENT_FACTOR"] = string(case.fcore_segment_factor)
    environment["CHANG_FCORE_CHORD_FACTOR"] = string(case.fcore_chord_factor)

    # Aerodynamic-model switches used directly by chang_case.jl.
    environment["CHANG_INTERACTION"] = isnothing(interaction_on) ? get(
        ENV,
        "CHANG_CONVERGENCE_INTERACTION",
        "false",
    ) : string(interaction_on)
    environment["CHANG_NEAR_FIELD_FORCE_MODEL"] = get(
        ENV,
        "CHANG_CONVERGENCE_FORCE_MODEL",
        "imperial",
    )
    environment["CHANG_PROP_MOMENT_PROJECTION"] = get(
        ENV,
        "CHANG_CONVERGENCE_PROP_MOMENT_PROJECTION",
        "exact_virtual_work",
    )
    environment["CHANG_PLOT_RESULTS"] = "false"
    environment["CHANG_ANIMATE_WAKE"] = "false"

    # A gated aeroelastic study supplies these explicitly so an unrelated
    # interactive-session override cannot change the structural algorithm or
    # partitioned-coupling relaxation between aerodynamic and aeroelastic
    # stages.  Legacy callers retain their existing environment/default path.
    !isnothing(ga_rho_inf) &&
        (environment["CHANG_GA_RHO_INF"] = string(ga_rho_inf))
    !isnothing(coupling_relaxation) &&
        (environment["CHANG_COUPLING_RELAXATION"] = string(coupling_relaxation))

    # Keep every case's history and console output separate. The hard angle is
    # a safety stop, not a convergence criterion.
    environment["CHANG_OUTPUT_DIR"] = output_directory
    environment["CHANG_OUTPUT_LABEL"] = output_label
    environment["CHANG_PROP_ANGLE_ABORT_DEG"] = string(hard_angle_deg)

    # Establish the periodic aerodynamic trim before applying the pitch
    # impulse whose free response is used for damping estimation.
    environment["CHANG_TRIM_REVOLUTIONS"] = get(
        ENV,
        "CHANG_CONVERGENCE_TRIM_REVOLUTIONS",
        "10.0",
    )
    environment["CHANG_TRIM_AVERAGE_REVOLUTIONS"] = get(
        ENV,
        "CHANG_CONVERGENCE_TRIM_AVERAGE_REVOLUTIONS",
        "1.0",
    )
    environment["CHANG_IMPULSE_MAGNITUDE"] = get(
        ENV,
        "CHANG_CONVERGENCE_IMPULSE_MAGNITUDE",
        "1000.0",
    )

    # An unrelated interactive-session override must not silently change the
    # convergence study's excitation timing. A convergence-specific override
    # remains available when an absolute start time is desired.
    pop!(environment, "CHANG_IMPULSE_START_S", nothing)
    if haskey(ENV, "CHANG_CONVERGENCE_IMPULSE_START_S")
        environment["CHANG_IMPULSE_START_S"] = ENV["CHANG_CONVERGENCE_IMPULSE_START_S"]
    end
    if haskey(ENV, "CHANG_CONVERGENCE_IMPULSE_DURATION_S")
        environment["CHANG_IMPULSE_DURATION_S"] = ENV["CHANG_CONVERGENCE_IMPULSE_DURATION_S"]
    end
    return environment
end

"""
    run_case(case; kwargs...)

Executes one coupled Chang simulation in a fresh Julia process and post-processes its response

# Arguments
- `case::UVLMConvergenceCase`: numerical settings for the child simulation

# Keyword arguments
- `output_directory`: directory for case histories and logs
- `speed_mps`: flow speed in m/s
- `trim_rpm`: propeller RPM at the trim reference speed
- `trim_speed_mps`: trim reference speed in m/s
- `radius_m`: propeller radius in m
- `end_time_s`: requested simulation end time in s
- `fit_start_s`: beginning of the moving-block fit window in s
- `fit_end_s`: end of the moving-block fit window in s
- `minimum_peaks`: minimum peaks required by the moving-block analysis
- `minimum_fit_r_squared`: minimum accepted moving-block fit quality
- `moving_block_initial_size`: initial FFT block size in samples
- `moving_block_size_ratio_lower`: lower block-to-record size ratio
- `moving_block_size_ratio_upper`: upper block-to-record size ratio
- `moving_block_peak_from_start`: first peak retained from the start
- `moving_block_peak_from_end`: peaks discarded from the end
- `hard_angle_deg`: propeller-angle safety-abort limit in degrees
- `moving_block_duration_s`: optional physical FFT-window duration in s
- `moving_block_overlap`: optional fractional overlap between FFT windows
- `moving_block_frequency_min_hz`, `moving_block_frequency_max_hz`: modal band
- `moving_block_apply_hann_window`: remove the mean and apply a Hann window
- `angle_of_attack_deg`: optional explicit angle of attack for the child case
- `interaction_on`: optional explicit wing--propeller interaction switch
- `ga_rho_inf`: optional explicit generalized-alpha spectral radius
- `coupling_relaxation`: optional explicit partitioned-coupling relaxation
"""
function run_case(
    case;
    output_directory,
    speed_mps,
    trim_rpm,
    trim_speed_mps,
    radius_m,
    end_time_s,
    fit_start_s,
    fit_end_s,
    minimum_peaks,
    minimum_fit_r_squared,
    moving_block_initial_size,
    moving_block_size_ratio_lower,
    moving_block_size_ratio_upper,
    moving_block_peak_from_start,
    moving_block_peak_from_end,
    hard_angle_deg,
    moving_block_duration_s = nothing,
    moving_block_overlap = nothing,
    moving_block_frequency_min_hz = 0.0,
    moving_block_frequency_max_hz = Inf,
    moving_block_apply_hann_window = false,
    angle_of_attack_deg = nothing,
    sideslip_deg = nothing,
    interaction_on = nothing,
    ga_rho_inf = nothing,
    coupling_relaxation = nothing,
)
    # The label is deterministic and carries family/level information into all
    # files produced by the child simulation.
    output_label = "$(case.family)_L$(case.level)_$(case.label)"
    history_path = joinpath(output_directory, output_label * "_history.csv")
    log_path = joinpath(output_directory, output_label * ".log")
    # Execute the same production entry point used by a standalone case.
    command = `$(Base.julia_cmd()) --startup-file=no --project=$(CONVERGENCE_PACKAGE_DIR) $(CONVERGENCE_CASE_RUNNER)`
    environment = child_environment(
        case,
        output_directory,
        output_label;
        speed_mps,
        trim_rpm,
        trim_speed_mps,
        end_time_s,
        hard_angle_deg,
        angle_of_attack_deg,
        sideslip_deg,
        interaction_on,
        ga_rho_inf,
        coupling_relaxation,
    )

    println(
        "\n[convergence] $(case.family) level $(case.level): $(case.label) " *
        "-> $log_path",
    )
    start_time = time()
    # Redirect verbose coupled-solver output to a per-case log. wait(process)
    # makes this driver sequential by default and limits simultaneous memory.
    process_succeeded = open(log_path, "w") do stream
        process = run(
            pipeline(setenv(command, environment), stdout = stream, stderr = stream);
            wait = false,
        )
        wait(process)
        return success(process)
    end
    elapsed_s = time() - start_time

    # Information common to successful and failed cases. time_step_s is
    # reconstructed here only for the report; the child calculated its own dt.
    base_result = (;
        case,
        reused = false,
        reused_from = "",
        elapsed_s,
        time_step_s = physical_time_step(
            case,
            speed_mps,
            trim_rpm,
            trim_speed_mps,
            radius_m,
        ),
        history_path,
        log_path,
    )

    if !process_succeeded || !isfile(history_path)
        reason = process_succeeded ? "history file was not produced" : "child process failed"
        return merge(base_result, (;
            status = "failed",
            reason,
            completed = false,
            max_pitch_deg = NaN,
            max_yaw_deg = NaN,
            pitch = unavailable_metrics(),
            yaw = unavailable_metrics(),
        ))
    end

    # The child has finished. Read its response, check basic integrity, then
    # evaluate pitch and yaw independently with identical moving-block rules.
    response = read_response_history(history_path)
    all_finite = all(isfinite, response.time) && all(isfinite, response.pitch) &&
        all(isfinite, response.yaw)
    completed = completed_requested_window(response.time, end_time_s)
    max_pitch = isempty(response.pitch) ? NaN : maximum(abs, response.pitch)
    max_yaw = isempty(response.yaw) ? NaN : maximum(abs, response.yaw)
    pitch = moving_block_metrics(
        response.time,
        response.pitch;
        fit_start_s,
        fit_end_s,
        minimum_peaks,
        minimum_fit_r_squared,
        initial_block_size = moving_block_initial_size,
        size_ratio_lower = moving_block_size_ratio_lower,
        size_ratio_upper = moving_block_size_ratio_upper,
        peak_from_start = moving_block_peak_from_start,
        peak_from_end = moving_block_peak_from_end,
        block_duration_s = moving_block_duration_s,
        block_overlap = moving_block_overlap,
        frequency_min_hz = moving_block_frequency_min_hz,
        frequency_max_hz = moving_block_frequency_max_hz,
        apply_hann_window = moving_block_apply_hann_window,
    )
    yaw = moving_block_metrics(
        response.time,
        response.yaw;
        fit_start_s,
        fit_end_s,
        minimum_peaks,
        minimum_fit_r_squared,
        initial_block_size = moving_block_initial_size,
        size_ratio_lower = moving_block_size_ratio_lower,
        size_ratio_upper = moving_block_size_ratio_upper,
        peak_from_start = moving_block_peak_from_start,
        peak_from_end = moving_block_peak_from_end,
        block_duration_s = moving_block_duration_s,
        block_overlap = moving_block_overlap,
        frequency_min_hz = moving_block_frequency_min_hz,
        frequency_max_hz = moving_block_frequency_max_hz,
        apply_hann_window = moving_block_apply_hann_window,
    )

    # Status precedence distinguishes solver failure, an incomplete/unsafe
    # response, and a completed simulation whose damping fit is inconclusive.
    if !all_finite
        status, reason = "failed", "non-finite response"
    elseif !completed
        status, reason = "failed", "simulation stopped before the requested end time"
    elseif max(max_pitch, max_yaw) >= hard_angle_deg
        status, reason = "failed", "hard propeller-angle limit reached"
    elseif !pitch.available || !yaw.available
        status, reason = "indeterminate", "moving-block analysis is unavailable for at least one channel"
    elseif !pitch.valid_for_convergence || !yaw.valid_for_convergence
        status, reason = "indeterminate", "at least one moving-block fit is below the R-squared threshold"
    else
        status, reason = "completed", "pitch and yaw moving-block damping fits are available"
    end

    return merge(base_result, (;
        status,
        reason,
        completed,
        max_pitch_deg = max_pitch,
        max_yaw_deg = max_yaw,
        pitch,
        yaw,
    ))
end

# ### Convergence comparison and output

"""
    finite_or_blank(value)

Formats a finite real number for CSV and report output, or returns a blank string

# Arguments
- `value`: value to format
"""
finite_or_blank(value) = value isa Real && isfinite(value) ? @sprintf("%.10g", value) : ""

"""
    csv_field(value)

Escapes a value as one CSV field when it contains delimiters

# Arguments
- `value`: value to convert to CSV text
"""
function csv_field(value)
    text = string(value)
    if occursin(',', text) || occursin('"', text) || occursin('\n', text)
        return "\"" * replace(text, "\"" => "\"\"") * "\""
    end
    return text
end

"""
    channel_error(value,reference)

Returns the absolute scalar error, or `NaN` when either input is non-finite

# Arguments
- `value`: value to assess
- `reference`: comparison value
"""
function channel_error(value, reference)
    return isfinite(value) && isfinite(reference) ? abs(value - reference) : NaN
end

"""
    annotated_results(results,lambda_tolerance)

Appends reference errors and convergence flags to the case results

# Arguments
- `results`: evaluated case-result records
- `lambda_tolerance`: maximum accepted pitch or yaw lambda error in 1/s
"""
function annotated_results(results, lambda_tolerance)
    # Add convergence errors after simulations have been evaluated. The final
    # configured level in each family is the reference. Both pitch and yaw
    # must be valid, and the larger of their absolute lambda errors controls
    # acceptance so one channel cannot hide a poorly converged other channel.
    annotated = NamedTuple[]
    families = unique(result.case.family for result in results)
    for family in families
        group = sort(filter(result -> result.case.family == family, results); by = r -> r.case.level)
        reference = last(group)
        reference_is_valid = reference.pitch.valid_for_convergence &&
            reference.yaw.valid_for_convergence && reference.status == "completed"
        # Error to the finest reference answers "how far from the best tested
        # result?" Change from previous answers "did the last refinement still
        # move the answer?" Both quantities are written to the CSV.
        previous = nothing

        for result in group
            pitch_reference_error = reference_is_valid ? channel_error(
                result.pitch.moving_block_lambda_per_s,
                reference.pitch.moving_block_lambda_per_s,
            ) : NaN
            yaw_reference_error = reference_is_valid ? channel_error(
                result.yaw.moving_block_lambda_per_s,
                reference.yaw.moving_block_lambda_per_s,
            ) : NaN
            maximum_reference_error = result.pitch.valid_for_convergence &&
                result.yaw.valid_for_convergence && reference_is_valid ?
                max(pitch_reference_error, yaw_reference_error) : NaN

            if isnothing(previous)
                maximum_previous_change = NaN
            else
                pitch_change = channel_error(
                    result.pitch.moving_block_lambda_per_s,
                    previous.pitch.moving_block_lambda_per_s,
                )
                yaw_change = channel_error(
                    result.yaw.moving_block_lambda_per_s,
                    previous.yaw.moving_block_lambda_per_s,
                )
                maximum_previous_change = result.pitch.valid_for_convergence &&
                    result.yaw.valid_for_convergence &&
                    previous.pitch.valid_for_convergence &&
                    previous.yaw.valid_for_convergence ? max(pitch_change, yaw_change) : NaN
            end

            push!(annotated, merge(result, (;
                reference_label = reference.case.label,
                pitch_lambda_reference_error_per_s = pitch_reference_error,
                yaw_lambda_reference_error_per_s = yaw_reference_error,
                maximum_lambda_reference_error_per_s = maximum_reference_error,
                maximum_lambda_previous_change_per_s = maximum_previous_change,
                within_lambda_tolerance = isfinite(maximum_reference_error) &&
                    maximum_reference_error <= lambda_tolerance,
            )))
            previous = result
        end
    end
    return annotated
end

"""
    write_summary_csv(path,results)

Writes the machine-readable convergence table

# Arguments
- `path`: destination CSV path
- `results`: annotated case-result records
"""
function write_summary_csv(path, results)
    # This is the machine-readable output. It contains numerical inputs,
    # moving-block diagnostics, reference errors, status, and paths back to the
    # raw history/log so every convergence point remains traceable.
    headers = [
        "family", "level", "label", "status", "reason", "reused", "reused_from",
        "wing_span_panels", "wing_chord_panels", "prop_radial_panels",
        "prop_chord_panels", "wake_revolutions", "wake_rows",
        "finite_core_segment_factor", "finite_core_chord_factor", "azimuth_step_deg", "time_step_s",
        "pitch_moving_block_lambda_per_s", "yaw_moving_block_lambda_per_s",
        "pitch_damping_percent", "yaw_damping_percent",
        "pitch_frequency_hz", "yaw_frequency_hz", "pitch_fit_r_squared",
        "yaw_fit_r_squared", "pitch_peak_count", "yaw_peak_count",
        "pitch_block_size", "yaw_block_size", "pitch_block_count", "yaw_block_count",
        "pitch_analyzed_start_s", "pitch_analyzed_end_s",
        "yaw_analyzed_start_s", "yaw_analyzed_end_s",
        "max_pitch_deg", "max_yaw_deg", "reference_label",
        "pitch_lambda_reference_error_per_s", "yaw_lambda_reference_error_per_s",
        "maximum_lambda_reference_error_per_s",
        "maximum_lambda_previous_change_per_s", "within_lambda_tolerance",
        "elapsed_s", "history_path", "log_path",
    ]
    open(path, "w") do stream
        println(stream, join(headers, ','))
        for result in results
            case = result.case
            values = [
                case.family,
                case.level,
                case.label,
                result.status,
                result.reason,
                result.reused,
                result.reused_from,
                case.wing_span_panels,
                case.wing_chord_panels,
                case.prop_radial_panels,
                case.prop_chord_panels,
                finite_or_blank(case.wake_revolutions),
                wake_rows(case),
                finite_or_blank(case.fcore_segment_factor),
                finite_or_blank(case.fcore_chord_factor),
                finite_or_blank(case.azimuth_step_deg),
                finite_or_blank(result.time_step_s),
                finite_or_blank(result.pitch.moving_block_lambda_per_s),
                finite_or_blank(result.yaw.moving_block_lambda_per_s),
                finite_or_blank(result.pitch.damping_ratio_percent),
                finite_or_blank(result.yaw.damping_ratio_percent),
                finite_or_blank(result.pitch.frequency_hz),
                finite_or_blank(result.yaw.frequency_hz),
                finite_or_blank(result.pitch.fit_r_squared),
                finite_or_blank(result.yaw.fit_r_squared),
                result.pitch.peak_count,
                result.yaw.peak_count,
                result.pitch.block_size,
                result.yaw.block_size,
                result.pitch.block_count,
                result.yaw.block_count,
                finite_or_blank(result.pitch.analyzed_start_s),
                finite_or_blank(result.pitch.analyzed_end_s),
                finite_or_blank(result.yaw.analyzed_start_s),
                finite_or_blank(result.yaw.analyzed_end_s),
                finite_or_blank(result.max_pitch_deg),
                finite_or_blank(result.max_yaw_deg),
                result.reference_label,
                finite_or_blank(result.pitch_lambda_reference_error_per_s),
                finite_or_blank(result.yaw_lambda_reference_error_per_s),
                finite_or_blank(result.maximum_lambda_reference_error_per_s),
                finite_or_blank(result.maximum_lambda_previous_change_per_s),
                result.within_lambda_tolerance,
                finite_or_blank(result.elapsed_s),
                result.history_path,
                result.log_path,
            ]
            println(stream, join(csv_field.(values), ','))
        end
    end
    return path
end

"""
    first_converged_level(group,tolerance)

Returns the first level for which that result and every finer result remain within tolerance

# Arguments
- `group`: ordered results belonging to one parameter family
- `tolerance`: maximum accepted absolute lambda error in 1/s
"""
function first_converged_level(group, tolerance)
    # A point is reported as the first converged level only when it and every
    # subsequent/finer level remain within tolerance. This avoids declaring a
    # fortuitous intermediate crossing as convergence.
    for index in eachindex(group)
        remaining = group[index:end]
        if all(
            result -> result.status == "completed" &&
                isfinite(result.maximum_lambda_reference_error_per_s) &&
                result.maximum_lambda_reference_error_per_s <= tolerance,
            remaining,
        )
            return index
        end
    end
    return nothing
end

"""
    write_text_report(path,results; kwargs...)

Writes the human-readable convergence report

# Arguments
- `path`: destination report path
- `results`: annotated case-result records

# Keyword arguments
- `lambda_tolerance`: maximum accepted pitch or yaw lambda error in 1/s
- `speed_mps`: study flow speed in m/s
- `end_time_s`: simulation duration in s
- `fit_start_s`: beginning of the moving-block fit window in s
- `fit_end_s`: end of the moving-block fit window in s
- `minimum_fit_r_squared`: minimum accepted moving-block fit quality
- `moving_block_initial_size`: initial FFT block size in samples
- `moving_block_size_ratio_lower`: lower block-to-record size ratio
- `moving_block_size_ratio_upper`: upper block-to-record size ratio
- `moving_block_peak_from_start`: first peak retained from the start
- `moving_block_peak_from_end`: peaks discarded from the end
"""
function write_text_report(
    path,
    results;
    lambda_tolerance,
    speed_mps,
    end_time_s,
    fit_start_s,
    fit_end_s,
    minimum_fit_r_squared,
    moving_block_initial_size,
    moving_block_size_ratio_lower,
    moving_block_size_ratio_upper,
    moving_block_peak_from_start,
    moving_block_peak_from_end,
)
    # Human-readable companion to the CSV. It records the method and controls
    # alongside the numerical assessment so the result can be interpreted
    # without reopening this source file.
    families = unique(result.case.family for result in results)
    open(path, "w") do stream
        println(stream, "Chang UVLM damping-convergence study")
        println(stream, "====================================")
        println(stream, "Generated: $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
        println(stream, "Flow speed: $speed_mps m/s")
        println(stream, "Simulation duration: $end_time_s s")
        println(stream, "Damping fit window: $fit_start_s to $fit_end_s s")
        println(stream, "Minimum fit R-squared: $minimum_fit_r_squared")
        println(stream, "Initial moving-block size: $moving_block_initial_size samples")
        println(
            stream,
            "Accepted block/sample ratio: $moving_block_size_ratio_lower to " *
            "$moving_block_size_ratio_upper",
        )
        println(
            stream,
            "Peak cut: start=$moving_block_peak_from_start, end=$moving_block_peak_from_end",
        )
        println(stream, "Moving-block lambda tolerance: $lambda_tolerance 1/s")
        println(stream)
        println(stream, "Method")
        println(stream, "------")
        println(stream, "The signal is cut between selected positive peaks.")
        println(stream, "For each overlapping FFT block, the maximum single-sided magnitude Xbar")
        println(stream, "is retained. A linear fit of log(Xbar) against block start time gives")
        println(stream, "the moving-block slope lambda in 1/s, which is the convergence metric.")
        println(stream, "lambda < 0 denotes decay and lambda > 0 denotes growth.")
        println(stream, "The derived modal ratio is zeta=-lambda/sqrt(lambda^2+(2*pi*f_d)^2),")
        println(stream, "so positive zeta denotes stable damping.")
        println(stream, "The last configured level of each family is its comparison reference.")
        println(stream, "Time-step cases retain a constant wake age in rotor revolutions by")
        println(stream, "changing the maximum number of wake rows with the azimuth step.")
        println(stream)

        for family in families
            group = sort(filter(result -> result.case.family == family, results); by = r -> r.case.level)
            println(stream, "Family: $family")
            println(stream, repeat('-', 8 + length(string(family))))
            println(
                stream,
                "label | pitch lambda (1/s) | yaw lambda (1/s) | pitch zeta (%) | " *
                "yaw zeta (%) | max lambda error (1/s) | status",
            )
            for result in group
                println(
                    stream,
                    "$(result.case.label) | " *
                    "$(finite_or_blank(result.pitch.moving_block_lambda_per_s)) | " *
                    "$(finite_or_blank(result.yaw.moving_block_lambda_per_s)) | " *
                    "$(finite_or_blank(result.pitch.damping_ratio_percent)) | " *
                    "$(finite_or_blank(result.yaw.damping_ratio_percent)) | " *
                    "$(finite_or_blank(result.maximum_lambda_reference_error_per_s)) | " *
                    result.status,
                )
            end

            first_index = first_converged_level(group, lambda_tolerance)
            if isnothing(first_index)
                println(stream, "Assessment: convergence could not be assessed from the valid cases.")
            elseif first_index == length(group)
                println(
                    stream,
                    "Assessment: convergence was not demonstrated before the finest tested level; " *
                    "only the reference itself is within tolerance.",
                )
            else
                println(
                    stream,
                    "Assessment: $(group[first_index].case.label) is the first level after which " *
                    "all finer levels remain within tolerance.",
                )
            end
            println(stream)
        end

        println(stream, "Interpretation notes")
        println(stream, "--------------------")
        println(stream, "- The wing spanwise panel count is also the collocated structural beam")
        println(stream, "  discretization in this implementation. The wing-mesh family therefore")
        println(stream, "  measures coupled aero-structural discretization sensitivity, not a purely")
        println(stream, "  aerodynamic grid change.")
        println(stream, "- The finite core is proportional to local 3-D vortex-segment length. That family is a")
        println(stream, "  regularization-sensitivity study; the smallest core is not automatically")
        println(stream, "  the most physically correct core.")
        println(stream, "- A poor moving-block linear-fit R-squared, too few peaks, or an")
        println(stream, "  insufficient analyzed segment makes a case indeterminate.")
        println(stream, "  Increase the end time or adjust the fit start before interpreting it.")
    end
    return path
end

"""
    plot_convergence(path,results)

Creates and saves the pitch and yaw moving-block convergence figure

# Arguments
- `path`: destination PNG path
- `results`: annotated case-result records
"""
function plot_convergence(path, results)
    # One panel per parameter family. Lambda, rather than derived zeta, is
    # plotted because lambda is the direct moving-block output and the actual
    # quantity used by the convergence criterion.
    families = unique(result.case.family for result in results)
    panels = Any[]
    for family in families
        group = sort(filter(result -> result.case.family == family, results); by = r -> r.case.level)
        x = collect(eachindex(group))
        labels = [result.case.label for result in group]
        pitch = [result.pitch.valid_for_convergence ?
            result.pitch.moving_block_lambda_per_s : NaN for result in group]
        yaw = [result.yaw.valid_for_convergence ?
            result.yaw.moving_block_lambda_per_s : NaN for result in group]
        panel = plot(
            x,
            pitch;
            marker = :circle,
            linewidth = 2,
            label = "pitch",
            xticks = (x, labels),
            xrotation = 25,
            ylabel = "Moving-block lambda (1/s)",
            title = replace(string(family), '_' => ' '),
            framestyle = :box,
            gridalpha = 0.25,
        )
        plot!(panel, x, yaw; marker = :diamond, linewidth = 2, label = "yaw")
        push!(panels, panel)
    end

    columns = min(2, length(panels))
    rows = ceil(Int, length(panels) / columns)
    figure = plot(
        panels...;
        layout = (rows, columns),
        size = (750 * columns, 430 * rows),
        plot_title = "Chang UVLM damping convergence",
    )
    savefig(figure, path)
    return path
end

"""
    print_case_matrix(cases,speed_mps,trim_rpm,trim_speed_mps,radius_m)

Prints the complete planned convergence-study matrix

# Arguments
- `cases`: generated convergence cases
- `speed_mps`: study flow speed in m/s
- `trim_rpm`: propeller RPM at the trim reference speed
- `trim_speed_mps`: trim reference speed in m/s
- `radius_m`: propeller radius in m
"""
function print_case_matrix(cases, speed_mps, trim_rpm, trim_speed_mps, radius_m)
    println("Generated $(length(cases)) family entries ($(length(unique(case_key.(cases)))) unique simulations):")
    for case in cases
        dt = physical_time_step(case, speed_mps, trim_rpm, trim_speed_mps, radius_m)
        println(
            @sprintf(
                "  %-15s L%d %-22s wing=%dx%d prop=%dx%d wake=%g rev/%d rows core=max(%g ds, %g c) dpsi=%g deg dt=%.7g s",
                string(case.family),
                case.level,
                case.label,
                case.wing_span_panels,
                case.wing_chord_panels,
                case.prop_radial_panels,
                case.prop_chord_panels,
                case.wake_revolutions,
                wake_rows(case),
                case.fcore_segment_factor,
                case.fcore_chord_factor,
                case.azimuth_step_deg,
                dt,
            ),
        )
    end
end

"""
    print_case_result(result)

Prints one compact terminal summary of lambda, damping ratio, frequency, and fit quality

# Arguments
- `result`: evaluated case-result record
"""
function print_case_result(result)
    println(
        @sprintf(
            "[convergence] %-13s pitch/yaw lambda=%+.6f/%+.6f 1/s, elapsed=%.1f s",
            uppercase(result.status),
            result.pitch.moving_block_lambda_per_s,
            result.yaw.moving_block_lambda_per_s,
            result.elapsed_s,
        ),
    )
    println(
        @sprintf(
            "              pitch/yaw zeta=%+.5f/%+.5f%%, f=%.4f/%.4f Hz, R2=%.3f/%.3f; %s",
            result.pitch.damping_ratio_percent,
            result.yaw.damping_ratio_percent,
            result.pitch.frequency_hz,
            result.yaw.frequency_hz,
            result.pitch.fit_r_squared,
            result.yaw.fit_r_squared,
            result.reason,
        ),
    )
end

"""
    main()

Runs the complete Chang UVLM damping-convergence workflow
"""
function main()
    # A. Expand the requested families into the complete ordered case list.
    families = requested_families()
    cases = build_cases(families)

    # B. Read controls that are common to all cases. These values describe the
    # operating point and damping analysis, whereas the case objects describe
    # the UVLM parameter being refined. The convergence study uses one flow
    # speed per invocation; CHANG_CONVERGENCE_SPEED_MPS selects that speed.
    speed_mps = convergence_float("CHANG_CONVERGENCE_SPEED_MPS", 80.0)
    trim_rpm = convergence_float("CHANG_CONVERGENCE_TRIM_RPM", 1217.6962)
    trim_speed_mps = convergence_float("CHANG_CONVERGENCE_TRIM_SPEED_MPS", 65.0)
    # Radius is fixed by the Chang case; it cancels from the
    # constant-advance-ratio time-step relation but is retained explicitly in
    # physical_time_step for clarity.
    radius_m = 1.15
    end_time_s = convergence_float("CHANG_CONVERGENCE_END_TIME_S", 2.0)
    fit_start_s = convergence_float("CHANG_CONVERGENCE_FIT_START_S", 0.7)
    fit_end_s = convergence_float("CHANG_CONVERGENCE_FIT_END_S", end_time_s)
    minimum_peaks = convergence_int("CHANG_CONVERGENCE_MIN_PEAKS", 4)
    minimum_fit_r_squared = convergence_float("CHANG_CONVERGENCE_MIN_FIT_R2", 0.5)
    moving_block_initial_size = convergence_int(
        "CHANG_CONVERGENCE_MOVING_BLOCK_SIZE",
        MOVING_BLOCK_INITIAL_SIZE,
    )
    moving_block_size_ratio_lower = convergence_float(
        "CHANG_CONVERGENCE_MOVING_BLOCK_RATIO_LOWER",
        MOVING_BLOCK_SIZE_RATIO_LOWER,
    )
    moving_block_size_ratio_upper = convergence_float(
        "CHANG_CONVERGENCE_MOVING_BLOCK_RATIO_UPPER",
        MOVING_BLOCK_SIZE_RATIO_UPPER,
    )
    moving_block_peak_from_start = convergence_int(
        "CHANG_CONVERGENCE_MOVING_BLOCK_PEAK_FROM_START",
        MOVING_BLOCK_PEAK_FROM_START,
    )
    moving_block_peak_from_end = convergence_int(
        "CHANG_CONVERGENCE_MOVING_BLOCK_PEAK_FROM_END",
        MOVING_BLOCK_PEAK_FROM_END,
    )
    hard_angle_deg = convergence_float("CHANG_CONVERGENCE_ABORT_ANGLE_DEG", 15.0)
    lambda_tolerance = convergence_float(
        "CHANG_CONVERGENCE_LAMBDA_TOLERANCE_PER_S",
        MOVING_BLOCK_LAMBDA_TOLERANCE_PER_S,
    )
    dry_run = convergence_bool("CHANG_CONVERGENCE_DRY_RUN", false)
    run_stamp = Dates.format(now(), "yyyymmdd_HHMMSS")
    output_directory = abspath(get(
        ENV,
        "CHANG_CONVERGENCE_OUTPUT_DIR",
        joinpath(CONVERGENCE_EXAMPLE_DIR, "output", "uvlm_convergence_$run_stamp"),
    ))

    # C. Reject invalid inputs before starting any expensive child simulation.
    speed_mps > 0 || error("CHANG_CONVERGENCE_SPEED_MPS must be positive")
    trim_rpm > 0 || error("CHANG_CONVERGENCE_TRIM_RPM must be positive")
    trim_speed_mps > 0 || error("CHANG_CONVERGENCE_TRIM_SPEED_MPS must be positive")
    0 <= fit_start_s < fit_end_s <= end_time_s || error(
        "The damping-fit window must satisfy 0 <= start < end <= simulation end time",
    )
    minimum_peaks >= 3 || error("CHANG_CONVERGENCE_MIN_PEAKS must be at least 3")
    0 <= minimum_fit_r_squared <= 1 || error("CHANG_CONVERGENCE_MIN_FIT_R2 must be in [0, 1]")
    ispow2(moving_block_initial_size) || error(
        "CHANG_CONVERGENCE_MOVING_BLOCK_SIZE must be a positive power of two",
    )
    0 < moving_block_size_ratio_lower <= moving_block_size_ratio_upper <= 1 || error(
        "Moving-block ratios must satisfy 0 < lower <= upper <= 1",
    )
    moving_block_peak_from_start >= 1 || error(
        "CHANG_CONVERGENCE_MOVING_BLOCK_PEAK_FROM_START must be at least 1",
    )
    moving_block_peak_from_end >= 0 || error(
        "CHANG_CONVERGENCE_MOVING_BLOCK_PEAK_FROM_END must be nonnegative",
    )
    hard_angle_deg > 0 || error("CHANG_CONVERGENCE_ABORT_ANGLE_DEG must be positive")
    lambda_tolerance >= 0 || error("Moving-block lambda tolerance must be nonnegative")

    println("Chang UVLM damping-convergence study")
    println("Flow speed: $speed_mps m/s; fit window: $fit_start_s:$fit_end_s s")
    println(
        "Moving block: initial size=$moving_block_initial_size, ratio=" *
        "$moving_block_size_ratio_lower:$moving_block_size_ratio_upper",
    )
    println("Lambda convergence tolerance: $lambda_tolerance 1/s")
    print_case_matrix(cases, speed_mps, trim_rpm, trim_speed_mps, radius_m)

    # Dry-run mode intentionally stops here. It is useful for checking exactly
    # which cases, wake rows, and physical time steps would be used.
    dry_run && return cases

    # D. A timestamped default directory prevents a new study from overwriting
    # earlier histories. Every child writes into this common study directory
    # with a unique family/level label.
    mkpath(output_directory)
    println("Results directory: $output_directory")
    summary_path = joinpath(output_directory, "uvlm_convergence_summary.csv")
    report_path = joinpath(output_directory, "uvlm_convergence_report.txt")
    figure_path = joinpath(output_directory, "uvlm_convergence_damping.png")

    # E. Run cases sequentially. Identical nominal configurations occur in
    # several families; case_key() lets the cache reuse their response instead
    # of spending time on duplicate simulations.
    cache = Dict{Tuple,NamedTuple}()
    results = NamedTuple[]
    for case in cases
        key = case_key(case)
        if haskey(cache, key)
            source_result = cache[key]
            println(
                "\n[convergence] Reusing $(source_result.case.family)/" *
                "$(source_result.case.label) for $(case.family)/$(case.label)",
            )
            result = merge(source_result, (;
                case,
                reused = true,
                reused_from = "$(source_result.case.family)/$(source_result.case.label)",
                elapsed_s = 0.0,
            ))
        else
            result = run_case(
                case;
                output_directory,
                speed_mps,
                trim_rpm,
                trim_speed_mps,
                radius_m,
                end_time_s,
                fit_start_s,
                fit_end_s,
                minimum_peaks,
                minimum_fit_r_squared,
                moving_block_initial_size,
                moving_block_size_ratio_lower,
                moving_block_size_ratio_upper,
                moving_block_peak_from_start,
                moving_block_peak_from_end,
                hard_angle_deg,
            )
            cache[key] = result
        end
        push!(results, result)
        print_case_result(result)

        # Rewrite a partial CSV after every case. If a long study is stopped,
        # the completed results and their file paths are still available.
        write_summary_csv(summary_path, annotated_results(results, lambda_tolerance))
    end

    # F. Once every family has its true finest reference, recompute annotations
    # and generate the final tabular, narrative, and graphical products.
    results = annotated_results(results, lambda_tolerance)
    write_summary_csv(summary_path, results)
    write_text_report(
        report_path,
        results;
        lambda_tolerance,
        speed_mps,
        end_time_s,
        fit_start_s,
        fit_end_s,
        minimum_fit_r_squared,
        moving_block_initial_size,
        moving_block_size_ratio_lower,
        moving_block_size_ratio_upper,
        moving_block_peak_from_start,
        moving_block_peak_from_end,
    )
    try
        plot_convergence(figure_path, results)
    catch exception
        @warn "Could not generate convergence figure" exception
    end

    println("\nConvergence summary: $summary_path")
    println("Convergence report:  $report_path")
    isfile(figure_path) && println("Convergence figure:  $figure_path")
    return results
end
