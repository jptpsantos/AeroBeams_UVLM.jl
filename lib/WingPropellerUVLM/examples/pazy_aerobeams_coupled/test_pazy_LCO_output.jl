# Filesystem and driver-routing regression; no physical aeroelastic simulation.
# julia --project=. lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/test_pazy_LCO_output.jl
using Test
include(joinpath(@__DIR__, "PazyLCOOutput.jl"))
using .PazyLCOOutput

@testset "Angle-separated, non-overwriting LCO directories" begin
    mktempdir() do base
        @test basename(pazy_lco_angle_directory(base, 3)) == "AoA3"
        @test basename(pazy_lco_angle_directory(base, rad2deg(deg2rad(3.0)))) == "AoA3"
        @test basename(pazy_lco_angle_directory(base, 7)) == "AoA7"
        @test basename(pazy_lco_angle_directory(base, -3.5)) == "AoAm3p5"
        @test basename(pazy_lco_angle_directory(base, -0.0)) == "AoA0"
        @test_throws ArgumentError pazy_lco_angle_directory(base, NaN)

        # Legacy 5-degree data can be found but cannot satisfy a 3/7-degree request.
        settings = joinpath(base, "pazy_linear_LCO_settings.txt")
        write(settings, "Root incidence [deg]: 5.0\n")
        @test pazy_lco_latest_output_directory(base, 5) == base
        @test_throws ArgumentError pazy_lco_latest_output_directory(base, 3)
        @test_throws ArgumentError pazy_lco_latest_output_directory(base, 7)
        legacy5 = pazy_lco_angle_directory(base, 5)
        mkpath(legacy5)
        write(joinpath(legacy5, "pazy_linear_LCO_settings.txt"), read(settings))
        @test pazy_lco_latest_output_directory(base, 5) == legacy5

        first3 = pazy_lco_new_run_directory(base, 3)
        write(joinpath(first3, "history.csv"), "preserve this run")
        second3 = pazy_lco_new_run_directory(base, 3)
        first7 = pazy_lco_new_run_directory(base, 7)
        @test first3 == joinpath(base, "AoA3", "run_001")
        @test second3 == joinpath(base, "AoA3", "run_002")
        @test first7 == joinpath(base, "AoA7", "run_001")
        @test pazy_lco_latest_output_directory(base, 3) == second3
        @test pazy_lco_latest_output_directory(base, 7) == first7
        @test read(joinpath(first3, "history.csv"), String) == "preserve this run"
        @test read(settings, String) == "Root incidence [deg]: 5.0\n"
        @test_throws ArgumentError pazy_lco_prepare_output_directory(first3;
            base_directory=base, angle_of_attack_deg=3)
        @test pazy_lco_prepare_output_directory(second3;
            base_directory=base, angle_of_attack_deg=3) == second3
        explicit = joinpath(base, "custom_empty")
        @test pazy_lco_prepare_output_directory(explicit;
            base_directory=base, angle_of_attack_deg=3) == explicit
        @test isdir(explicit)

        # Run ordering remains numeric above 999, and new runs never fill old gaps.
        mkpath(joinpath(base, "AoA3", "run_999"))
        mkpath(joinpath(base, "AoA3", "run_1000"))
        @test pazy_lco_latest_output_directory(base, 3) == joinpath(base, "AoA3", "run_1000")
        @test pazy_lco_new_run_directory(base, 3) == joinpath(base, "AoA3", "run_1001")
    end
end

# These minimal stand-ins exercise the real CSV/settings/plot export and both
# entry points without loading a structural model or running a physical sweep.
module LCOSweepOutputHarness
using DelimitedFiles
module UVLM
create_UVLMSolver(; kwargs...) = nothing
end
module AeroBeams
geometrical_properties_Pazy() = (nothing, nothing, 0.1, nothing)
end
const requested_angles = Float64[]
function run_pazy_wing_uvlm_static(; angle_of_attack, kwargs...)
    push!(requested_angles, angle_of_attack)
    return (; history_columns=["test_only"], history=ones(1, 1),
        termination_reason="mock output regression", converged=true,
        geometry_residual=0.0, load_residual=0.0)
end
function run_pazy_wing_uvlm(; angle_of_attack, duration, airspeed_profile, kwargs...)
    push!(requested_angles, angle_of_attack)
    time = [0.0, duration/2, duration]
    return (; time, airspeed_history=airspeed_profile.(time), root_axial_strain=zeros(3),
        root_curvature_2=ones(3), root_curvature_3=zeros(3),
        tip_out_of_plane=zeros(3), tip_twist_degrees=zeros(3))
end
include(joinpath(@__DIR__, "PazyLinearLCOSweep.jl"))
end

@testset "Ramp durations preserve the 5-degree speed slopes" begin
    ramp_duration = LCOSweepOutputHarness.pazy_lco_ramp_duration
    airspeed = LCOSweepOutputHarness.pazy_linear_airspeed
    @test ramp_duration(35.0, 52.0, 0.425) == 40.0
    @test ramp_duration(52.0, 35.0, 0.85) == 20.0
    for (low, high) in ((40.0, 61.0), (32.0, 50.0), (35.5, 52.8))
        up_time = ramp_duration(low, high, 0.425)
        down_time = ramp_duration(high, low, 0.85)
        @test up_time ≈ 2down_time
        t = collect(range(0, up_time; length=11))
        speeds = airspeed.(t; start_speed=low, end_speed=high, duration=up_time)
        @test first(speeds) ≈ low
        @test last(speeds) ≈ high
        @test all(isapprox.(diff(speeds)./diff(t), 0.425; atol=1e-12))
        t = collect(range(0, down_time; length=11))
        speeds = airspeed.(t; start_speed=high, end_speed=low, duration=down_time)
        @test all(isapprox.(diff(speeds)./diff(t), -0.85; atol=1e-12))
        # Fractional durations must still align to the shared integration step.
        nominal_dt = 0.25*0.1/high
        dt = LCOSweepOutputHarness.pazy_lco_aligned_timestep(up_time, down_time, nominal_dt)
        @test nominal_dt/2 <= dt <= nominal_dt
        @test round(Int, up_time/dt)*dt ≈ up_time
        @test round(Int, down_time/dt)*dt ≈ down_time
    end
    for rate in (0.0, -0.425, NaN, Inf)
        @test_throws ArgumentError ramp_duration(35.0, 52.0, rate)
    end
    @test_throws ArgumentError ramp_duration(35.0, 35.0, 0.425)
    @test_throws ArgumentError ramp_duration(-35.0, 52.0, 0.425)

    mktempdir() do base
        # Direct API calls also preserve the rate, including a separate descending start.
        result = LCOSweepOutputHarness.run_pazy_linear_LCO(; initial_speed=40.0,
            peak_speed=61.0, downsweep_initial_speed=63.0, output_base_directory=base)
        @test last(result.upsweep.result.time) ≈ 21/0.425
        @test last(result.downsweep.result.time) ≈ 23/0.85
        settings = read(joinpath(result.output_directory, "pazy_linear_LCO_settings.txt"), String)
        @test occursin("Upsweep rate [m/s^2]: 0.425", settings)
        @test occursin("Downsweep rate [m/s^2]: -0.85", settings)
        # Short duration-based diagnostics retain their existing behavior.
        explicit = LCOSweepOutputHarness.run_pazy_linear_LCO(; initial_speed=35.0,
            peak_speed=52.0, upsweep_duration=0.006, downsweep_duration=0.003,
            output_base_directory=base)
        @test last(explicit.upsweep.result.time) ≈ 0.006
        @test last(explicit.downsweep.result.time) ≈ 0.003
    end
end

@testset "Both strain entry points route 3/7 degrees and preserve prior data" begin
    mktempdir() do base
        legacy = joinpath(base, "LCO_downsweep_pulse")
        mkpath(joinpath(legacy, "AoA5"))
        for directory in (legacy, joinpath(legacy, "AoA5"))
            write(joinpath(directory, "pazy_linear_LCO_settings.txt"), "Root incidence [deg]: 5.0\n")
            write(joinpath(directory, "pazy_linear_LCO_upsweep.csv"), "old 5-degree upsweep")
            write(joinpath(directory, "pazy_linear_LCO_downsweep.csv"), "old 5-degree downsweep")
        end
        original = Dict(joinpath(root, name) => read(joinpath(root, name))
            for (root, _, names) in walkdir(legacy) for name in names)
        previous_run = get(ENV, "PAZY_LCO_DOWNSWEEP_PULSE_RUN", nothing)
        ENV["PAZY_LCO_DOWNSWEEP_PULSE_RUN"] = "true"
        try
            for (run_number, driver) in enumerate(("run_pazy_LCO_strains.jl", "run_pazy_LCO_downsweep_pulse.jl"))
                path = joinpath(@__DIR__, driver)
                for degrees in (3.0, 7.0)
                    low, high = degrees == 3.0 ? (40.0, 61.0) : (32.0, 50.0)
                    # Keep the real entry-point configuration/control flow, routing
                    # its output into the temporary test tree and selecting the angle.
                    source = replace(read(path, String),
                        "include(joinpath(@__DIR__, \"PazyLinearLCOSweep.jl\"))" => "",
                        "angle_of_attack_deg = 5.0" => "angle_of_attack_deg = $degrees",
                        "sweep_initial_speed = 35.0" => "sweep_initial_speed = $low",
                        "sweep_peak_speed = 52.0" => "sweep_peak_speed = $high",
                        "joinpath(@__DIR__, \"output\", output_name)" => repr(legacy))
                    calls_before = length(LCOSweepOutputHarness.requested_angles)
                    Base.include_string(LCOSweepOutputHarness, source, path)
                    result = LCOSweepOutputHarness.pulse_lco_results
                    expected = joinpath(legacy, "AoA$(Int(degrees))", "run_" * lpad(string(run_number), 3, '0'))
                    @test result.output_directory == expected
                    angles = LCOSweepOutputHarness.requested_angles[calls_before+1:end]
                    @test length(angles) == 4 # initialization + dynamics in each branch
                    @test all(angle -> angle ≈ deg2rad(degrees), angles)
                    settings = read(joinpath(expected, "pazy_linear_LCO_settings.txt"), String)
                    @test occursin("Root incidence [deg]: $(rad2deg(deg2rad(degrees)))", settings)
                    @test isfile(joinpath(expected, "pazy_linear_LCO_upsweep.csv"))
                    @test isfile(joinpath(expected, "pazy_linear_LCO_downsweep.csv"))
                    @test isfile(joinpath(expected, "pazy_linear_LCO_root_strain_vs_airspeed.png"))
                    @test isfile(joinpath(expected, "pazy_linear_LCO_root_strain_history.png"))
                    for (branch, expected_rate) in ((result.upsweep, 0.425), (result.downsweep, -0.85))
                        @test all(isapprox.(diff(branch.result.airspeed_history)./diff(branch.result.time),
                            expected_rate; atol=1e-12))
                    end
                    @test LCOSweepOutputHarness.pulse_start ≈ (high-49.0)/0.85

                    # Plot-only must resolve this angle's existing run, not allocate one.
                    plot_source = replace(source, "plot_only = false" => "plot_only = true")
                    Base.include_string(LCOSweepOutputHarness, plot_source, path)
                    @test LCOSweepOutputHarness.output_directory == expected
                    @test length(LCOSweepOutputHarness.requested_angles) == calls_before+4
                end
            end
            calls_before = length(LCOSweepOutputHarness.requested_angles)
            @test_throws ArgumentError LCOSweepOutputHarness.run_pazy_linear_LCO(;
                angle_of_attack=deg2rad(7.0), output_directory=legacy)
            @test length(LCOSweepOutputHarness.requested_angles) == calls_before
            for (path, contents) in original
                @test read(path) == contents
            end
        finally
            isnothing(previous_run) ? delete!(ENV, "PAZY_LCO_DOWNSWEEP_PULSE_RUN") :
                (ENV["PAZY_LCO_DOWNSWEEP_PULSE_RUN"] = previous_run)
        end
    end
end
