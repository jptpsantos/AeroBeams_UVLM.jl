# Postprocess saved root-strain histories; no aeroelastic solver is loaded.
# From the repository root:
# julia --project=. lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/plot_pazy_LCO_spectrogram.jl
# Optional command-line arguments: INPUT_DIRECTORY [upsweep|downsweep|both]
module PazyLCOSpectrogram

using DelimitedFiles
using FFTW: rfft
using LinearAlgebra: dot
using Statistics: mean
import Plots
include(joinpath(@__DIR__, "PazyLCOOutput.jl"))
using .PazyLCOOutput: pazy_lco_latest_output_directory

export pazy_strain_stft, plot_pazy_LCO_spectrogram

"""One-sided, Hann-windowed STFT of signed strain samples in microstrain.

Only complete windows are used. Each window is detrended before applying the
Hann taper. Non-DC/non-Nyquist amplitudes are `2abs(FFT)/sum(window)`; DC and
Nyquist are not doubled. Thus an on-bin sinusoid of peak amplitude A microstrain
has amplitude A, up to negative-frequency leakage. This is an amplitude
spectrum, not a power spectral density. `amplitude_db = 20log10(A/reference)`.

The frame airspeed is interpolated at the window's center time. Samples must
have uniform, increasing times. Zero padding refines the displayed frequency
grid; the window duration still determines the actual frequency resolution.
"""
function pazy_strain_stft(time, airspeed, strain_ue;
    window_seconds=1.0, overlap=0.95, nfft=nothing,
    frequency_limits=(0.0, 100.0), detrend=:linear,
    reference_strain_ue=1.0, db_floor=-80.0)

    length(time) == length(airspeed) == length(strain_ue) ||
        throw(ArgumentError("Time, airspeed and strain must have equal lengths"))
    length(time) >= 4 || throw(ArgumentError("At least four samples are required"))
    all(isfinite, time) && all(isfinite, airspeed) && all(isfinite, strain_ue) ||
        throw(ArgumentError("Time, airspeed and strain must be finite"))
    isfinite(window_seconds) && window_seconds > 0 ||
        throw(ArgumentError("window_seconds must be positive and finite"))
    isfinite(overlap) && 0 <= overlap < 1 ||
        throw(ArgumentError("overlap must satisfy 0 <= overlap < 1"))
    isfinite(reference_strain_ue) && reference_strain_ue > 0 ||
        throw(ArgumentError("reference_strain_ue must be positive and finite"))
    isfinite(db_floor) || throw(ArgumentError("db_floor must be finite"))
    detrend in (:linear, :mean, :none) ||
        throw(ArgumentError("detrend must be :linear, :mean, or :none"))

    steps = diff(time)
    dt = mean(steps)
    all(>(0), steps) || throw(ArgumentError("Time must be strictly increasing"))
    all(step -> isapprox(step, dt; rtol=1e-6, atol=1e-12), steps) ||
        throw(ArgumentError("STFT requires uniformly sampled time; resample first"))
    fmin, fmax = frequency_limits
    isfinite(fmin) && isfinite(fmax) && 0 <= fmin < fmax <= 0.5/dt ||
        throw(ArgumentError("frequency_limits must lie between 0 and Nyquist ($(0.5/dt) Hz)"))

    window_samples = round(Int, window_seconds/dt)
    4 <= window_samples <= length(time) ||
        throw(ArgumentError("The window must contain at least four samples and fit in the record"))
    hop_samples = max(1, round(Int, (1-overlap)*window_samples))
    if isnothing(nfft)
        nfft = nextpow(2, 4window_samples)
    end
    nfft isa Integer && nfft >= window_samples ||
        throw(ArgumentError("nfft must be an integer at least as large as the window"))

    # Symmetric Hann window; sum normalization preserves sinusoidal peak amplitude.
    sample_offsets = collect(0:window_samples-1) .- (window_samples-1)/2
    taper = 0.5 .- 0.5cos.(2pi .* (0:window_samples-1) ./ (window_samples-1))
    taper_sum = sum(taper)
    slope_denominator = dot(sample_offsets, sample_offsets)
    all_frequencies = collect(0:fld(nfft, 2)) ./ (nfft*dt)
    frequency_bins = findall(f -> fmin <= f <= fmax, all_frequencies)
    isempty(frequency_bins) && throw(ArgumentError("No FFT bins in frequency_limits"))
    starts = 1:hop_samples:(length(time)-window_samples+1)
    frame_time = Vector{Float64}(undef, length(starts))
    frame_speed = similar(frame_time)
    amplitude_ue = Matrix{Float64}(undef, length(frequency_bins), length(starts))
    buffer = zeros(nfft)

    for (column, start) in enumerate(starts)
        stop = start+window_samples-1
        segment = Float64.(strain_ue[start:stop])
        if detrend != :none
            segment .-= mean(segment)
            if detrend == :linear
                slope = dot(sample_offsets, segment)/slope_denominator
                segment .-= slope .* sample_offsets
            end
        end
        fill!(buffer, 0)
        buffer[1:window_samples] .= segment .* taper
        transform = rfft(buffer)
        for (row, bin) in enumerate(frequency_bins)
            is_endpoint = bin == 1 || (iseven(nfft) && bin == length(transform))
            amplitude_ue[row, column] = (is_endpoint ? 1 : 2)*abs(transform[bin])/taper_sum
        end
        center = (start+stop)/2
        left, right = floor(Int, center), ceil(Int, center)
        frame_time[column] = (time[left]+time[right])/2
        frame_speed[column] = (airspeed[left]+airspeed[right])/2
    end
    # Floor only the dB representation, preserving the original amplitude array.
    amplitude_db = max.(db_floor,
        20 .* log10.(amplitude_ue ./ reference_strain_ue))
    return (; frequency_hz=all_frequencies[frequency_bins], time_s=frame_time,
        airspeed_mps=frame_speed, amplitude_ue, amplitude_db,
        dt, window_samples, hop_samples, nfft,
        window_seconds=window_samples*dt, overlap=1-hop_samples/window_samples,
        frequency_bin_hz=1/(nfft*dt), detrend, reference_strain_ue, db_floor)
end

function read_strain_history(path)
    isfile(path) || throw(ArgumentError("Saved branch CSV not found: $path"))
    data, header = readdlm(path, ',', Float64; header=true)
    names = vec(header)
    function column(name; required=true)
        index = findfirst(==(name), names)
        if isnothing(index)
            required && throw(ArgumentError("Missing column '$name' in $path"))
            return nothing
        end
        return data[:, index]
    end
    return (; time=column("time_s"), airspeed=column("airspeed_mps"),
        strain_ue=column("root_longitudinal_strain_ue"),
        force=column("tip_force_A1_N"; required=false))
end

function saved_incidence(input_directory)
    path = joinpath(input_directory, "pazy_linear_LCO_settings.txt")
    isfile(path) || return nothing
    matched = match(r"Root incidence \[deg\]:\s*([-+0-9.eE]+)", read(path, String))
    return isnothing(matched) ? nothing : parse(Float64, matched.captures[1])
end

"""Plot saved Pazy root-strain STFTs against airspeed, as in paper Fig. 4.

Defaults read the newest 5-degree `output/LCO_downsweep_pulse/AoA5/run_NNN`
upsweep, falling back to legacy 5-degree data. Set `angle_of_attack_deg=3.0`
or `7.0` to select another angle, or supply an exact `input_directory`.
Use `plot_branches=:both` or `:downsweep` for the separate pulsed branch.
Incidence is read from the saved settings, rather than the current run script.
PNG, SVG and processing settings are saved under `spectrograms/`. Optional
`save_spectra=true` also writes a long-format CSV with unclipped amplitudes.
`color_limits` clips colors only; `db_floor` limits the saved dB representation.

Including this file only loads functions. No histories are concatenated,
no absolute value is taken before the FFT, and no per-frame normalization is
applied. Frequency content and its airspeed evolution come from the saved data.
"""
function plot_pazy_LCO_spectrogram(;
    angle_of_attack_deg=5.0, input_directory=nothing, output_directory=nothing,
    plot_branches=:upsweep, window_seconds=1.0, overlap=0.95, nfft=nothing,
    frequency_limits=(0.0, 100.0), detrend=:linear,
    reference_strain_ue=1.0, db_floor=-80.0,
    color_limits=(-40.0, 60.0), airspeed_limits=:auto,
    mark_pulse=true, save_spectra=false)

    branch_names = plot_branches == :both ? (:upsweep, :downsweep) :
        plot_branches in (:upsweep, :downsweep) ? (plot_branches,) :
        throw(ArgumentError("plot_branches must be :upsweep, :downsweep, or :both"))
    all(isfinite, color_limits) && color_limits[1] < color_limits[2] ||
        throw(ArgumentError("color_limits must be a finite increasing pair"))
    if isnothing(input_directory)
        input_directory = pazy_lco_latest_output_directory(
            joinpath(@__DIR__, "output", "LCO_downsweep_pulse"), angle_of_attack_deg)
    end
    output_directory = isnothing(output_directory) ?
        joinpath(input_directory, "spectrograms") : output_directory
    incidence = saved_incidence(input_directory)
    results = Dict{Symbol, Any}()
    # Load all requested branches before writing any files, so missing data is explicit.
    histories = Dict(name => read_strain_history(joinpath(input_directory,
        "pazy_linear_LCO_$(name).csv")) for name in branch_names)

    for name in branch_names
        history = histories[name]
        speed_steps = diff(history.airspeed)
        all(>(0), speed_steps) || all(<(0), speed_steps) ||
            throw(ArgumentError("Airspeed must be strictly monotonic within each branch"))
        spectrum = pazy_strain_stft(history.time, history.airspeed, history.strain_ue;
            window_seconds, overlap, nfft, frequency_limits, detrend,
            reference_strain_ue, db_floor)
        length(spectrum.time_s) >= 2 ||
            throw(ArgumentError("At least two complete STFT windows are needed to plot"))
        # Sort frames only after computing the STFT in chronological order. This
        # keeps descending-speed frames aligned with their amplitudes and times.
        order = sortperm(spectrum.airspeed_mps)
        speed = spectrum.airspeed_mps[order]
        limits = airspeed_limits == :auto ? extrema(history.airspeed) : airspeed_limits
        angle_label = isnothing(incidence) ? "" : ", α = $(incidence)°"
        title = "AeroBeams + UVLM$(angle_label) ($(name))"
        figure = Plots.heatmap(speed, spectrum.frequency_hz, spectrum.amplitude_db[:, order];
            xlabel="Airspeed [m/s]", ylabel="Frequency [Hz]", title,
            color=:jet, clims=color_limits, colorbar=true,
            colorbar_title="Strain amplitude [dB re $(reference_strain_ue) με]",
            xlims=limits, ylims=frequency_limits, legend=false, grid=false,
            fontfamily="Computer Modern", guidefontsize=13, tickfontsize=11,
            titlefontsize=13, size=(900, 600), dpi=250,
            margin=6Plots.mm, right_margin=12Plots.mm)
        pulse_indices = isnothing(history.force) ? Int[] : findall(!iszero, history.force)
        if mark_pulse && !isempty(pulse_indices)
            pulse_speed = history.airspeed[first(pulse_indices)]
            Plots.vline!(figure, [pulse_speed]; color=:white, linestyle=:dash,
                linewidth=1, label=false)
            label_height = frequency_limits[2]-0.06*(frequency_limits[2]-frequency_limits[1])
            Plots.annotate!(figure, pulse_speed, label_height,
                Plots.text("tip pulse", 10, :white, :right))
        end

        mkpath(output_directory)
        stem = joinpath(output_directory, "pazy_linear_LCO_root_strain_spectrogram_$(name)")
        png_path, svg_path, settings_path = stem*".png", stem*".svg", stem*"_settings.txt"
        Plots.savefig(figure, png_path)
        Plots.savefig(figure, svg_path)
        spectra_path = save_spectra ? stem*".csv" : nothing
        if save_spectra
            open(spectra_path, "w") do io
                println(io, "time_s,airspeed_mps,frequency_hz,amplitude_ue,amplitude_db")
                for column in eachindex(spectrum.time_s), row in eachindex(spectrum.frequency_hz)
                    println(io, join((spectrum.time_s[column], spectrum.airspeed_mps[column],
                        spectrum.frequency_hz[row], spectrum.amplitude_ue[row, column],
                        spectrum.amplitude_db[row, column]), ','))
                end
            end
        end
        open(settings_path, "w") do io
            println(io, "Source: ", abspath(joinpath(input_directory, "pazy_linear_LCO_$(name).csv")))
            println(io, "Branch: $name; root incidence [deg]: $incidence")
            println(io, "Signal: signed root_longitudinal_strain_ue; near-root first-element midpoint gauge")
            println(io, "Sampling step [s]: $(spectrum.dt); sample count: $(length(history.time))")
            println(io, "Window: symmetric Hann; requested [s]: $window_seconds; actual N*dt [s]: $(spectrum.window_seconds)")
            println(io, "Window samples: $(spectrum.window_samples); hop samples: $(spectrum.hop_samples); actual overlap: $(spectrum.overlap)")
            println(io, "FFT length: $(spectrum.nfft); padded bin spacing [Hz]: $(spectrum.frequency_bin_hz)")
            println(io, "Nominal frequency resolution scale 1/(N*dt) [Hz]: $(1/spectrum.window_seconds)")
            println(io, "Window-center time range [s]: $(extrema(spectrum.time_s)); frame count: $(length(spectrum.time_s))")
            println(io, "Window-center speed range [m/s]: $(extrema(spectrum.airspeed_mps))")
            println(io, "Airspeed change across a window [m/s], maximum: ",
                maximum(abs(history.airspeed[i+spectrum.window_samples-1]-history.airspeed[i])
                    for i in 1:spectrum.hop_samples:length(history.time)-spectrum.window_samples+1))
            println(io, "Per-window detrending: $detrend; no global static solution is subtracted")
            println(io, "One-sided peak amplitude: 2*abs(rfft(window*signal))/sum(window); DC/Nyquist factor=1")
            println(io, "dB = 20*log10(amplitude_ue/$reference_strain_ue); dB floor: $db_floor")
            println(io, "Color limits [dB]: $color_limits; frequency limits [Hz]: $frequency_limits")
            println(io, "Only complete windows; no time padding/extrapolation. Descending frames sorted for display only.")
            println(io, "STFT window, overlap, detrending and dB reference are processing choices; the paper does not specify them.")
            if !isempty(pulse_indices)
                println(io, "Applied nonzero tip-force samples: time range [s] ",
                    extrema(history.time[pulse_indices]), "; speed range [m/s] ",
                    extrema(history.airspeed[pulse_indices]))
                println(io, "The pulse and its transient are retained; its start is marked when mark_pulse=true.")
            end
        end
        results[name] = (; figure, spectrum, png_path, svg_path, settings_path, spectra_path)
        println("$name spectrogram: $png_path")
        println("  $(length(spectrum.time_s)) frames; window $(spectrum.window_seconds) s; bin spacing $(spectrum.frequency_bin_hz) Hz")
    end
    return results
end

end # module

using .PazyLCOSpectrogram: pazy_strain_stft, plot_pazy_LCO_spectrogram

if abspath(PROGRAM_FILE) == @__FILE__
    # Edit these processing/plotting choices when running from the IDE.
    angle_of_attack_deg = 5.0             # selects the latest saved run for this angle
    input_directory = isempty(ARGS) ? nothing : ARGS[1]
    plot_branches = length(ARGS) < 2 ? :upsweep : Symbol(ARGS[2])
    spectrogram_results = plot_pazy_LCO_spectrogram(;
        angle_of_attack_deg, input_directory, plot_branches,
        window_seconds=1.0, overlap=0.95,
        frequency_limits=(0.0, 100.0), color_limits=(-40.0, 60.0),
        detrend=:linear, reference_strain_ue=1.0,
        airspeed_limits=:auto, save_spectra=false)
end
