# Configuration isolation and interleaved aerodynamic trials in one Julia session.
# Run with --project=lib/WingPropellerUVLM from the repository root.
using Test
using LinearAlgebra
using WingPropellerUVLM: snapshot_uvlm, legacy_imperial_segment_forces!

include(joinpath(@__DIR__, "..", "src", "ChangAeroelastic.jl"))
using .ChangAeroelastic
using .ChangAeroelastic: chang_aerodynamic_options, chang_excitation_options,
    chang_integration_options

BLAS.set_num_threads(1)

@testset "One resolved configuration per case" begin
    defaults = chang_case_defaults()
    environment = Dict(
        "CHANG_WING_SPAN_PANELS" => "4", "CHANG_WING_CHORD_PANELS" => "2",
        "CHANG_PROP_RADIAL_PANELS" => "2", "CHANG_PROP_CHORD_PANELS" => "2",
        "CHANG_FREESTREAM_SPEED_MPS" => "80.0", "CHANG_SPEED_MPS" => "85.0",
        "CHANG_END_TIME_S" => "0.003", "CHANG_PLOT_RESULTS" => "false",
        "CHANG_ANIMATE_WAKE" => "false", "CHANG_INTERACTION" => "true",
        "CHANG_FCORE_SEGMENT_FACTOR" => "0.02", "CHANG_FCORE_CHORD_FACTOR" => "0.03",
        "CHANG_WAKE_ROWS_PROPELLER" => "8", "CHANG_TRIM_REVOLUTIONS" => "0.5",
        "CHANG_TRIM_AVERAGE_REVOLUTIONS" => "0.25",
        "CHANG_GA_RHO_INF" => "0.8", "CHANG_COUPLING_TOL_U" => "2e-6",
        "CHANG_COUPLING_MAX_ITER" => "12", "CHANG_ANIMATION_FPS" => "20",
    )
    first_config = load_chang_configuration(defaults; env = environment)
    integer_defaults = merge(defaults, (
        simulation = merge(defaults.simulation, (; freestream_speed_mps = 80)),
    ))
    @test load_chang_configuration(integer_defaults;
        env = Dict("CHANG_SPEED_MPS" => "82.5")).simulation.freestream_speed_mps == 82.5
    @test first_config.simulation.freestream_speed_mps == 85.0
    alias_environment = copy(environment)
    delete!(alias_environment, "CHANG_SPEED_MPS")
    @test load_chang_configuration(defaults; env = alias_environment).simulation.freestream_speed_mps == 80.0
    @test first_config.wake.maximum_rows_wing == defaults.wake.wing_rows_per_chord_panel * 2
    @test first_config.wake.maximum_rows_propeller == 8
    @test first_config.output.animation_fps == 20
    @test first_config.structural.damping_reference_frequency_hz == defaults.structural.pitch_frequency_hz

    # Editing defaults or the source environment after loading cannot change a case.
    defaults.propeller.attachment_eta[1] = 0.4
    environment["CHANG_FCORE_CHORD_FACTOR"] = "0.9"
    @test first_config.propeller.attachment_eta == [0.83]
    withenv("CHANG_FCORE_CHORD_FACTOR" => "0.8", "CHANG_GA_RHO_INF" => "0.1") do
        @test chang_aerodynamic_options(first_config).finite_core(2.0, 0.1) == 0.06
    end

    first_model = build_chang_model(first_config)
    first_workspace = build_chang_workspace(first_model)
    integration = chang_integration_options(first_model.config, first_model.structural, first_model.parameters)
    excitation = chang_excitation_options(first_model.config, first_model.parameters)
    @test integration.generalized_alpha.rho_inf == 0.8
    @test integration.coupling.state_tolerance == 2e-6
    @test integration.coupling.maximum_iterations == 12
    @test excitation.start_time ≈ 0.5 * 2pi / abs(first_model.parameters.Ω)

    # Distinct geometry, rotor count, force model, and wake limits in the same module.
    second_defaults = chang_case_defaults()
    second_defaults = merge(second_defaults, (
        propeller = merge(second_defaults.propeller, (; attachment_eta = [0.4, 0.8])),
    ))
    second_environment = merge(environment, Dict(
        "CHANG_SPEED_MPS" => "70.0", "CHANG_WING_SPAN_PANELS" => "5",
        "CHANG_WING_CHORD_PANELS" => "3", "CHANG_INTERACTION" => "false",
        "CHANG_NEAR_FIELD_FORCE_MODEL" => "legacy_imperial_segments",
        "CHANG_FCORE_CHORD_FACTOR" => "0.04", "CHANG_WAKE_ROWS_WING" => "6",
        "CHANG_WAKE_ROWS_PROPELLER" => "7", "CHANG_IMPULSE_START_S" => "0.1",
        "CHANG_PROP_MOMENT_PROJECTION" => "fixed_aero_axes",
    ))
    second_config = load_chang_configuration(second_defaults; env = second_environment)
    second_model = build_chang_model(second_config)
    second_workspace = build_chang_workspace(second_model)
    @test second_model.near_field_force_function === legacy_imperial_segment_forces!
    @test second_config.excitation.impulse_start_s == 0.1
    @test size(first_model.structural.M) == (26, 26)
    @test size(second_model.structural.M) == (34, 34)
    @test first_workspace.nwake == [20, 8, 8, 8, 8]
    @test second_workspace.nwake == [6, 7, 7, 7, 7, 7, 7, 7, 7]

    first_state = fill(1e-5, first_model.structural.ndof_free)
    second_state = fill(-2e-5, second_model.structural.ndof_free)
    first_snapshot = snapshot_uvlm(first_workspace.system)
    second_snapshot = snapshot_uvlm(second_workspace.system)
    before = aero_load_for_state!(first_model, first_workspace, first_snapshot, first_state, 1)
    other = aero_load_for_state!(second_model, second_workspace, second_snapshot, second_state, 1)
    after = aero_load_for_state!(first_model, first_workspace, first_snapshot, first_state, 1)
    @test all(isfinite, before) && norm(before) > 0
    @test all(isfinite, other) && norm(other) > 0
    @test before == after
    @test first_config.propeller.attachment_eta == [0.83]

    # Reusing a model still allocates independent mutable aerodynamic storage.
    another_workspace = build_chang_workspace(first_model)
    @test another_workspace.system !== first_workspace.system
    first_workspace.iwake[1] = 1
    @test another_workspace.iwake[1] == 0
    first_config.propeller.attachment_eta[1] = 0.2
    @test first_model.config.propeller.attachment_eta == [0.83]

    @test_throws ErrorException load_chang_configuration(; env = Dict("CHANG_GA_RHO_INF" => "bad"))
    @test_throws ErrorException load_chang_configuration(; env = Dict("CHANG_INTERACTION" => "maybe"))
    @test_throws ErrorException load_chang_configuration(; env = Dict("CHANG_SPEED_MPS" => "0"))
end
