# Standalone numerical checks for pitch-mode postprocessing; no simulation runs.
using Test
include(joinpath(@__DIR__,"..","analyze_pitch_phase.jl"))

@testset "Common-frequency pitch amplitude and phase" begin
    time = collect(0.0:0.002:5.0)
    τ = time .- 3.95
    for growth in (0.07,-0.12)
        wing = 0.03 .+ 0.2 .* exp.(growth .* τ) .* cos.(2π*4.3 .* τ .- deg2rad(30))
        propeller = -0.4 .+ 4.0 .* exp.(growth .* τ) .* cos.(2π*4.3 .* τ .+ deg2rad(20))
        mode = identify_pitch_mode(time,wing,propeller)
        @test mode.reliable
        @test mode.frequency_hz ≈ 4.3 atol=2e-5
        @test mode.growth_rate_per_s ≈ growth atol=1e-4
        @test mode.wing_amplitude_deg ≈ 0.2 atol=1e-5
        @test mode.propeller_amplitude_deg ≈ 4.0 atol=1e-4
        @test mode.wing_to_propeller_amplitude_ratio ≈ 0.05 atol=1e-5
        @test mode.wing_minus_propeller_phase_deg ≈ -50.0 atol=0.01
    end
    wing = 0.2 .* cos.(2π*4.3 .* τ .+ deg2rad(179))
    propeller = 4.0 .* cos.(2π*4.3 .* τ .- deg2rad(179))
    @test identify_pitch_mode(time,wing,propeller).wing_minus_propeller_phase_deg ≈ -2.0 atol=0.01

    weak = identify_pitch_mode(time,wing .* 1e-5,propeller)
    @test !weak.reliable
    @test isnan(weak.wing_minus_propeller_phase_deg)
    @test isnan(weak.wing_to_propeller_amplitude_ratio)
    mixed = identify_pitch_mode(time,wing .+ 0.2 .* sin.(2π*7.0 .* τ),propeller)
    @test !mixed.reliable
    @test isnan(mixed.wing_minus_propeller_phase_deg)
    @test maximum(mixed.normalized_residuals) > 0.1
    @test_throws ErrorException identify_pitch_mode(time,wing[1:end-1],propeller)
    @test_throws ErrorException identify_pitch_mode(reverse(time),wing,propeller)
    @test_throws ErrorException identify_pitch_mode(time,zero.(wing),propeller)
end

if "--results" in ARGS
    run_pitch_reviewer_analysis()
end
