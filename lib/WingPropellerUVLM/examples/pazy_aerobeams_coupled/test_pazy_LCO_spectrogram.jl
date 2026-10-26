# Signal-processing checks only; no aeroelastic simulation.
# julia --project=. lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/test_pazy_LCO_spectrogram.jl
using Test
using DelimitedFiles
include(joinpath(@__DIR__, "plot_pazy_LCO_spectrogram.jl"))

@testset "Calibrated strain STFT" begin
    time = collect(0:0.001:1.999)
    speed = 35 .+ 0.425time
    strain = 1000 .+ 7time .+ 100sin.(2pi*28time) .+
        40cos.(2pi*56time) .+ 20sin.(2pi*84time)
    spectrum = pazy_strain_stft(time, speed, strain;
        overlap=0.75, nfft=1000, frequency_limits=(0.0, 100.0))
    @test spectrum.window_samples == 1000
    @test spectrum.hop_samples == 250
    @test spectrum.time_s ≈ [0.4995, 0.7495, 0.9995, 1.2495, 1.4995]
    @test spectrum.airspeed_mps ≈ 35 .+ 0.425spectrum.time_s
    for (frequency, amplitude) in ((28.0, 100.0), (56.0, 40.0), (84.0, 20.0))
        bin = findfirst(==(frequency), spectrum.frequency_hz)
        @test all(isapprox.(spectrum.amplitude_ue[bin, :], amplitude; rtol=1e-3))
        @test all(isapprox.(spectrum.amplitude_db[bin, :], 20log10(amplitude); atol=0.01))
    end
    padded = pazy_strain_stft(time, speed, strain; overlap=0.75, nfft=4000)
    @test padded.amplitude_ue[1:4:end, :] ≈ spectrum.amplitude_ue
    @test padded.frequency_bin_hz ≈ 0.25
    referenced = pazy_strain_stft(time, speed, strain;
        overlap=0.75, nfft=1000, reference_strain_ue=10.0)
    bin28 = findfirst(==(28.0), spectrum.frequency_hz)
    @test referenced.amplitude_db[bin28, :] ≈ spectrum.amplitude_db[bin28, :] .- 20

    # A changing equilibrium contributes no oscillatory content after detrending.
    ramp = pazy_strain_stft(time, speed, 1000 .+ 7time; nfft=1000)
    @test maximum(ramp.amplitude_ue) < 1e-10
    @test all(==(-80.0), ramp.amplitude_db)
    @test all(isfinite, ramp.amplitude_db)
    down = pazy_strain_stft(time, 52 .- 0.85time, reverse(strain);
        overlap=0.75, nfft=1000)
    @test down.airspeed_mps ≈ 52 .- 0.85down.time_s
    @test down.amplitude_ue ≈ reverse(spectrum.amplitude_ue; dims=2) atol=1e-10

    # DC and the real-FFT Nyquist bin must not receive the one-sided factor of 2.
    endpoint_time = collect(0:511)./512
    endpoint_signal = 3 .+ 2 .* (-1.0).^(0:511)
    endpoints = pazy_strain_stft(endpoint_time, ones(512), endpoint_signal;
        nfft=512, frequency_limits=(0.0, 256.0), detrend=:none)
    @test endpoints.amplitude_ue[1, 1] ≈ 3.0 atol=1e-10
    @test endpoints.amplitude_ue[end, 1] ≈ 2.0 atol=1e-10

    @test_throws ArgumentError pazy_strain_stft(time[1:end-1], speed, strain)
    @test_throws ArgumentError pazy_strain_stft(time, speed, strain; overlap=1.0)
    @test_throws ArgumentError pazy_strain_stft(time, speed, strain; nfft=500)
    @test_throws ArgumentError pazy_strain_stft(time, speed, strain; window_seconds=3.0)
    @test_throws ArgumentError pazy_strain_stft(time, speed, strain; frequency_limits=(0, 501))
    uneven_time = copy(time)
    uneven_time[10] += 0.0001
    @test_throws ArgumentError pazy_strain_stft(uneven_time, speed, strain)
end

@testset "Saved branch postprocessing" begin
    mktempdir() do directory
        time = collect(0:0.002:1.998)
        for (name, speed) in (("upsweep", 35 .+ 0.425time), ("downsweep", 52 .- 0.85time))
            force = Float64.((0.2 .< time) .& (time .< 0.3))
            # Reordered columns ensure names, rather than column positions, are used.
            open(joinpath(directory, "pazy_linear_LCO_$(name).csv"), "w") do io
                println(io, "root_longitudinal_strain_ue,airspeed_mps,tip_force_A1_N,time_s")
                writedlm(io, hcat(100sin.(2pi*28time), speed, force, time), ',')
            end
        end
        write(joinpath(directory, "pazy_linear_LCO_settings.txt"), "Root incidence [deg]: 5.0\n")
        results = plot_pazy_LCO_spectrogram(; input_directory=directory,
            plot_branches=:both, overlap=0.5, nfft=500, save_spectra=true)
        for name in (:upsweep, :downsweep)
            result = results[name]
            @test filesize(result.png_path) > 1000
            @test filesize(result.svg_path) > 1000
            @test occursin("root incidence [deg]: 5.0", read(result.settings_path, String))
            @test occursin("Applied nonzero tip-force", read(result.settings_path, String))
            csv = readdlm(result.spectra_path, ',', Float64; skipstart=1)
            @test size(csv, 1) == length(result.spectrum.amplitude_ue)
            @test csv[1, 1:3] ≈ [result.spectrum.time_s[1],
                result.spectrum.airspeed_mps[1], result.spectrum.frequency_hz[1]]
            # Plot x-coordinates increase for both branches; the STFT stays chronological.
            series = first(result.figure.series_list)
            @test issorted(series[:x])
            @test series[:z].surf == result.spectrum.amplitude_db[:, sortperm(result.spectrum.airspeed_mps)]
        end
        @test !isfile(joinpath(directory, "pazy_linear_LCO_root_strain_vs_airspeed.png"))
    end
end
