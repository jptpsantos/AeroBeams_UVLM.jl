using Test, DelimitedFiles
include(joinpath(@__DIR__, "run_pazy_LCO_diagnostic.jl"))
BLAS.set_num_threads(1)

@testset "Constant-speed pulse diagnostic" begin
    mktempdir() do directory
        checked = run_pazy_LCO_diagnostic(; airspeeds=(20.0,), reference_speed=22.0,
            minimum_sweep_speed=20.0, duration=0.08, angle_of_attack=0.01,
            chordwise_panels=2, spanwise_panels=4,
            time_step_chords=0.25, wake_length_chords=2.0,
            pulse_start=0.005, pulse_duration=0.02, window_duration=0.02,
            equilibrium_options=(; maximum_static_iterations=80, static_relaxation=0.5,
                aerodynamic_force_tolerance=5e-4, aerodynamic_circulation_tolerance=5e-4,
                geometry_tolerance=2e-6, load_tolerance=5e-4),
            progress_frequency=10000, output_directory=directory)
        baseline, unforced, half_dt, double_wake = checked.cases
        @test length(checked.cases) == 4
        @test baseline.initial_strain == unforced.initial_strain
        @test baseline.dt == unforced.dt == double_wake.dt
        @test half_dt.dt == baseline.dt / 2
        @test half_dt.wake_rows * half_dt.dt ≈ baseline.wake_rows * baseline.dt
        @test double_wake.wake_rows ≈ 2baseline.wake_rows atol=1
        @test any(!iszero, baseline.force)
        @test all(iszero, unforced.force)
        for case in checked.cases
            data, header = readdlm(joinpath(case.case_directory, "history.csv"), ',', Float64; header=true)
            column(name) = data[:, findfirst(==(name), vec(header))]
            @test all(==(20.0), column("airspeed_mps"))
            @test column("root_longitudinal_strain_ue") ≈
                1e6 .* (column("root_axial_strain") .+ 0.001125 .* column("root_curvature_2_per_m"))
            @test all(column("coupling_geometry_residual") .<= 1e-4)
            @test all(column("coupling_load_residual") .<= 1e-4)
            @test all(isfinite, case.windows)
            @test isfinite(case.tail_slope)
            @test all(t -> 0.005 < t < 0.025,
                case.time[.!iszero.(case.force)])
        end
        @test filesize(joinpath(directory, "diagnostic_summary.csv")) > 100
        @test filesize(joinpath(directory, "strain_histories.png")) > 1000
        @test filesize(joinpath(directory, "window_amplitudes.png")) > 1000
    end
end
