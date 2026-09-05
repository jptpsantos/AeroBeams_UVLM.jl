module AerodynamicConvergenceTests
using Test

# Keep the library suite independent of a user's sweep environment.
withenv((name => nothing for name in keys(ENV) if startswith(name, "CHANG_AERO_"))...) do
    include(joinpath(@__DIR__, "..", "examples", "chang_linear_aeroelastic",
        "studies", "convergence", "run_chang_coupled_aerodynamic_sweep.jl"))
end

withenv((name => nothing for name in keys(ENV) if startswith(name, "CHANG_AERO_"))...) do
    @testset "Aerodynamic convergence integrity" begin
        defaults = ChangCoupledAerodynamicOptions()
        @test ChangAerodynamicStudy.options_dict(options_from_environment(Dict())) ==
            ChangAerodynamicStudy.options_dict(defaults)
        @test validate_options(defaults) == 72
        @test_throws ErrorException validate_options(ChangCoupledAerodynamicOptions(
            finite_core_segment_factor = 0, finite_core_chord_factor = 0))
        @test_throws ErrorException validate_options(ChangCoupledAerodynamicOptions(flow_speed_mps = Inf))
        @test_throws ErrorException validate_options(ChangCoupledAerodynamicOptions(azimuth_step_deg = 7))
        @test_throws ErrorException validate_options(ChangCoupledAerodynamicOptions(simulated_revolutions = 4))
        @test_throws ErrorException validate_options(ChangCoupledAerodynamicOptions(wake_relaxation = 1.1))
        @test_throws ErrorException sweep_levels("UNUSED_SWEEP_LEVELS", [1,1,2])
        @test number_token(0.0001) != number_token(0.0002)

        n = 72
        wave = sin.(2pi .* (1:n) ./ n)
        periodic = periodic_metrics(repeat(wave, 8), n, 2)
        @test periodic.periodic_rms == 0
        @test periodic.phase ≈ wave
        changing = periodic_metrics(vcat((i .* wave for i in 1:8)...), n, 2)
        @test abs(changing.drift) < 1e-14
        @test changing.periodic_rms > 0.7
        @test_throws ErrorException periodic_metrics(ones(73), n, 1)
        @test_throws ErrorException periodic_metrics(fill(NaN, 8n), n, 2)
        @test phase_rms_error(wave, wave) < 1e-14
        fine = sin.(2pi .* (1:360) ./ 360)
        @test phase_rms_error(wave, fine) < 0.001
        # Fine-grid oscillations hidden on the coarse grid must still be seen.
        aliased = sin.(2pi .* 12 .* (1:120) ./ 120)
        @test phase_rms_error(zeros(12), aliased) > 0.6

        cases = build_cases(requested_families())
        @test length(cases) == 22
        @test length(unique(case_key.(cases))) == 16
        @test all(c -> c.core_factor == 0, cases)
        @test all(c -> c.family == :finite_core || c.chord_core_factor == 0.01, cases)

        group = build_cases([:wing_span])
        record(case) = (; case, status = "completed", reason = "verified, \"quoted\"\nrow",
            reused = false, elapsed_s = 1.0, mean_cl = 0.3, mean_ct = 0.001, mean_cq = 0.0,
            std_cl = 0.0, std_ct = 0.0, drift_cl = 0.0, drift_ct = 0.0,
            periodic_cl = 0.0, periodic_ct = 0.0, periodic_cq = 0.0,
            phase_cl = fill(0.3, n), phase_ct = fill(0.001, n), phase_cq = zeros(n),
            history_path = "history,with,commas.csv", revolution_path = "revolutions.csv", log_path = "case.log")
        rows = record.(group)
        annotate(rows) = annotate_results(rows; cl_absolute_tolerance = 1e-4,
            ct_absolute_tolerance = 1e-6, relative_tolerance = 0.01)
        @test all(r -> r.within_tolerance, annotate(rows))
        @test !any(r -> r.within_tolerance, annotate(rows[1:2]))
        unsteady = copy(rows)
        unsteady[end] = merge(rows[end], (; periodic_ct = 2e-6))
        @test !any(r -> r.within_tolerance, annotate(unsteady))
        unstable_reference = copy(rows)
        unstable_reference[end] = merge(rows[end], (; mean_cl = 0.32))
        @test !any(r -> r.within_tolerance, annotate(unstable_reference))
        different_wave = copy(rows)
        different_wave[end] = merge(rows[end], (; phase_cl = 0.3 .+ 0.01wave))
        @test !any(r -> r.within_tolerance, annotate(different_wave))
        torque_change = copy(rows)
        torque_change[end] = merge(rows[end], (; mean_cq = 1e-4))
        @test !any(r -> r.within_tolerance, annotate(torque_change))

        mktempdir() do directory
            table = joinpath(directory, "summary.csv")
            write_summary(table, annotate(rows))
            raw, header = readdlm(table, ',', header = true)
            @test size(raw, 2) == length(header)
            @test raw[1, 5] == rows[1].reason

            case = group[1]
            options = case_options(case, directory; simulated_revolutions = 8, averaged_revolutions = 2)
            count = 8n
            dt = deg2rad(case.azimuth_deg) / (options.reference_rpm * 2pi / 60 *
                options.flow_speed_mps / options.reference_speed_mps)
            history_path = joinpath(directory, "coupled_aerodynamic_history.csv")
            open(history_path, "w") do io
                println(io, "time_s,azimuth_deg,wing_CL,propeller_CT,propeller_CQ")
                writedlm(io, hcat(collect(1:count) .* dt, mod.(collect(1:count) .* case.azimuth_deg, 360),
                    fill(0.3, count), fill(0.001, count), zeros(count)), ',')
            end
            write(joinpath(directory, "coupled_aerodynamic_revolutions.csv"), "revolution,mean_wing_CL,mean_propeller_CT\n")
            write(joinpath(directory, "coupled_aerodynamic_summary.txt"), "synthetic test data\n")
            write_metadata(directory, options)
            @test output_matches_options(directory, options; warn_on_mismatch = false)
            loaded = read_periodic_result(directory, case, 2; simulated_revolutions = 8)
            @test loaded.mean_cl ≈ 0.3
            @test loaded.periodic_cl == 0
            for override in ("CHANG_AERO_COLLECTIVE_OFFSET_DEG" => "1",
                             "CHANG_AERO_DENSITY_KGPM3" => "1.1",
                             "CHANG_AERO_BETA_DEG" => "1",
                             "CHANG_AERO_WAKE_RELAXATION" => "0.2",
                             "CHANG_AERO_SIMULATED_REVOLUTIONS" => "9",
                             "CHANG_AERO_AVERAGED_REVOLUTIONS" => "1")
                environment = child_environment(case, directory; simulated_revolutions = 8, averaged_revolutions = 2)
                environment[first(override)] = last(override)
                @test !output_matches_options(directory, options_from_environment(environment); warn_on_mismatch = false)
            end
            write(history_path, "truncated\n")
            @test !output_matches_options(directory, options; warn_on_mismatch = false)
            @test_throws ErrorException write_metadata(directory, options; expected_fingerprint = "stale")
        end
    end
end
end
