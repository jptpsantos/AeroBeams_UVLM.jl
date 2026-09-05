# Focused checks of the aerodynamic-to-aeroelastic handoff; no time marching.
using Test
include(joinpath(@__DIR__, "..", "studies", "convergence",
    "run_chang_aeroelastic_sweep_from_aerodynamic.jl"))

@testset "Aerodynamic-to-aeroelastic gate" begin
    row(level; expected = 4, chord = 0.04 / 2^(level - 1)) = (;
        family = :finite_core, level, expected_family_levels = expected,
        status = "completed", label = "core_$level", wing_span = 30, wing_chord = 10,
        prop_radial = 10, prop_chord = 10, wake_revolutions = 2.0,
        core_factor = 0.0, chord_core_factor = chord, azimuth_deg = 5.0,
        within_tolerance = true, periodic_converged = true, reference_stable = true,
        periodic_cl = 0.0, periodic_ct = 0.0, drift_cl = 0.0, drift_ct = 0.0)
    rows = [row(i) for i in 1:4]
    @test length(completed_family(rows, :finite_core)) == 4
    @test_throws ErrorException completed_family(rows[1:3], :finite_core)
    @test_throws ErrorException completed_family([rows[1], rows[2], rows[2], rows[4]], :finite_core)
    @test selected_row(rows; periodic_cl_tolerance = 1e-4, periodic_ct_tolerance = 1e-6,
        allow_finest_only = false) == rows[1]
    unsteady = [merge(r, (; periodic_ct = 2e-6)) for r in rows]
    @test_throws ErrorException selected_row(unsteady; periodic_cl_tolerance = 1e-4,
        periodic_ct_tolerance = 1e-6, allow_finest_only = false)
    config = (; wing_span = 30, wing_chord = 10, prop_radial = 10, prop_chord = 10,
        wake_revolutions = 2.0, core_factor = 0.0, chord_core_factor = 0.01, azimuth_deg = 5.0)
    case = convergence_case(:combined_validation, 1, "test", config)
    @test case.fcore_segment_factor == 0
    @test case.fcore_chord_factor == 0.01
    other = convergence_case(:combined_validation, 1, "test", merge(config, (; chord_core_factor = 0.02)))
    @test case_key(case) != case_key(other)
    environment = child_environment(case, "unused", "unused"; speed_mps = 65.0,
        trim_rpm = 1217.6962, trim_speed_mps = 65.0, end_time_s = 6.0, hard_angle_deg = 15.0,
        sideslip_deg = 0.0)
    @test environment["CHANG_FCORE_SEGMENT_FACTOR"] == "0.0"
    @test environment["CHANG_FCORE_CHORD_FACTOR"] == "0.01"
    @test environment["CHANG_SIDESLIP_DEG"] == "0.0"
    print_case_matrix([case], 65.0, 1217.6962, 65.0, 1.15)
end
