# Run from the repository root with --project=lib/WingPropellerUVLM.
using Test
using LinearAlgebra
using WingPropellerUVLM
include(joinpath(@__DIR__, "..", "src", "ChangAeroelastic.jl"))
using .ChangAeroelastic

BLAS.set_num_threads(1)

@testset "Chang isolated wing" begin
    defaults = chang_isolated_wing_case_defaults()
    original = chang_case_defaults()
    @test !isempty(original.propeller.attachment_eta)
    @test isempty(defaults.propeller.attachment_eta)
    @test isempty(defaults.simulation.impulse_propeller_indices)

    env = Dict(
        "CHANG_WING_SPAN_PANELS" => "2",
        "CHANG_WING_CHORD_PANELS" => "1",
        "CHANG_TIME_STEP_S" => "0.001",
        "CHANG_END_TIME_S" => "0.012",
        "CHANG_WAKE_ROWS_WING" => "32",
        "CHANG_IMPULSE_START_S" => "0.003",
        "CHANG_TRIM_AVERAGE_DURATION_S" => "0.002",
        "CHANG_IMPULSE_DURATION_S" => "0.004",
        "CHANG_IMPULSE_MAGNITUDE" => "0.1",
        "CHANG_PLOT_RESULTS" => "false",
        "CHANG_ANIMATE_WAKE" => "false",
        "CHANG_COUPLING_MAX_ITER" => "50",
    )
    config = load_chang_isolated_wing_configuration(defaults; env)
    model = build_chang_model(config)
    structural = model.structural
    @test model.parameters.Npropellers == 0
    @test model.parameters.Ω == 0
    @test structural.ndof_free == 6config.wing.spanwise_panels
    @test structural.ndof_prop_free == 0
    wing_free = 7:size(structural.Ms_W, 1)
    @test structural.M ≈ structural.Ms_W[wing_free, wing_free]
    @test structural.K ≈ structural.Ks_W[wing_free, wing_free]
    @test iszero(structural.C)
    @test iszero(structural.Gs_W)
    @test all(model.parameters.dt .== 0.001)

    # Rotor operating point, pylon properties, and blade mesh cannot affect a
    # bare wing's matrices, time grid, aerodynamic geometry, or trial load.
    changed = merge(defaults, (;
        propeller = merge(defaults.propeller, (;
            rotation_rpm = 2original.propeller.rotation_rpm,
            radial_panels = 1, chordwise_panels = 1,
        )),
        structural = merge(defaults.structural, (;
            pylon_length_m = 4.0, blade_mass_kg = 100.0,
            pylon_mass_per_length_kgpm = 100.0,
        )),
    ))
    other = build_chang_model(load_chang_isolated_wing_configuration(changed; env))
    @test other.structural.M == structural.M
    @test other.structural.K == structural.K
    @test other.structural.C == structural.C
    @test other.parameters.t == model.parameters.t
    @test other.parameters.dt == model.parameters.dt
    workspaces = [build_chang_workspace(m) for m in (model, other)]
    loads = map((model, other), workspaces) do m, w
        @test w.nsurf == 1
        @test isempty(w.nodal_forces_prop)
        @test isempty(w.prop_surface_indices)
        @test w.nwake == [32]
        @test size(w.system.surfaces[1]) == (1, 2)
        copy_surfaces_to_previous!(w.system, 1)
        return aero_load_for_state!(m, w, snapshot_uvlm(w.system),
            zeros(m.structural.ndof_free), 1)
    end
    @test all(isfinite, loads[1])
    @test loads[1] == loads[2]
    @test !iszero(loads[1])

    for integrator in (:newmark_beta, :generalized_alpha),
            scheme in (:loose_explicit, :implicit_predictor_corrector)
        mktempdir() do directory
            active_env = merge(env, Dict(
                "CHANG_TIME_INTEGRATOR" => string(integrator),
                "CHANG_COUPLING_SCHEME" => string(scheme),
                "CHANG_OUTPUT_DIR" => directory,
            ))
            run = run_chang(load_chang_isolated_wing_configuration(defaults; env = active_env))
            solution = run.solution
            @test run.results.validation_passed
            @test solution.last_step == length(run.model.parameters.dt)
            @test run.workspace.iwake == [solution.last_step]
            @test all(solution.coupling_converged)
            @test all(q -> all(isfinite, q), solution.displacement_history)
            @test all(iszero, solution.displacement_history[1:4])
            @test any(q -> !iszero(q), solution.displacement_history[5:end])
            @test isempty(run.results.propeller_pitch)
            @test isempty(run.results.propeller_yaw)
            @test !occursin("propeller_", first(readlines(run.results.history_path)))
            @test occursin("number_of_propellers = 0", read(run.results.summary_path, String))
            @test occursin("isolated_wing = true", read(run.results.summary_path, String))
        end
    end

    # Both alternative wing-moment directions select the expected free-node DOF.
    for dof in (4, 5, 6)
        variant = merge(defaults, (; excitation = merge(defaults.excitation, (;
            wing_impulse_dof = dof, wing_impulse_eta = 0.5,
        ))))
        active = load_chang_isolated_wing_configuration(variant; env)
        @test ChangAeroelastic.chang_excitation_options(active, model.parameters).dof_indices == [dof]
    end
    for dt in ("0", "-0.001", "NaN", "Inf", "1")
        @test_throws ErrorException load_chang_isolated_wing_configuration(defaults;
            env = merge(env, Dict("CHANG_TIME_STEP_S" => dt)))
    end
    @test_throws ErrorException load_chang_isolated_wing_configuration(original; env)
    limit = ChangAeroelastic.chang_state_limit_status(zeros(structural.ndof_free),
        structural.ndof_wing_free; state_norm_limit = 1.0, propeller_angle_limit_deg = Inf)
    @test !limit.exceeded
    @test limit.maximum_propeller_angle_deg == 0.0
end
