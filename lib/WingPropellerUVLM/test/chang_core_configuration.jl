module ChangCoreConfigurationTests
using Test
using StaticArrays
import WingPropellerUVLM as UVLM

const EXAMPLE = joinpath(@__DIR__, "..", "examples", "chang_linear_aeroelastic")
include(joinpath(EXAMPLE, "src", "ChangAeroelastic.jl"))
using .ChangAeroelastic
include(joinpath(EXAMPLE, "src", "chang_propeller_trim.jl"))

@testset "Chang fixed-core configuration and propagation" begin
    defaults = chang_case_defaults()
    empty_env = Dict{String,String}()
    config = load_chang_configuration(defaults; env = empty_env)
    @test config.aerodynamic.core_radius_m == defaults.aerodynamic.core_radius_m
    @test config.wing.symmetric == defaults.wing.symmetric
    @test_throws ErrorException load_chang_configuration(defaults;
        env=Dict("CHANG_WING_SYMMETRIC"=>"invalid"))

    # Fixed radius takes precedence even when both old factors are zero.
    environment = Dict(
        "CHANG_CORE_RADIUS_M" => "0.002",
        "CHANG_FCORE_SEGMENT_FACTOR" => "0",
        "CHANG_FCORE_CHORD_FACTOR" => "0",
        "CHANG_WING_SPAN_PANELS" => "4", "CHANG_WING_CHORD_PANELS" => "2",
        "CHANG_PROP_RADIAL_PANELS" => "2", "CHANG_PROP_CHORD_PANELS" => "2",
        "CHANG_WAKE_ROWS_WING" => "2", "CHANG_WAKE_ROWS_PROPELLER" => "2",
        "CHANG_END_TIME_S" => "0.001", "CHANG_PLOT_RESULTS" => "false",
        "CHANG_ANIMATE_WAKE" => "false",
    )
    fixed = load_chang_configuration(defaults; env = environment)
    radius = 0.002
    law = ChangAeroelastic.chang_aerodynamic_options(fixed).finite_core
    @test law(1.8, 0.375) == radius
    @test law(0.197, 0.01) == radius
    @test (@inferred law(0.197, 0.01)) == radius
    withenv("CHANG_CORE_RADIUS_M" => "0.9") do
        @test ChangAeroelastic.chang_aerodynamic_options(fixed).finite_core(2.0, 1.0) == radius
    end

    for value in ("0", "-0.1", "NaN", "Inf", "bad")
        @test_throws ErrorException load_chang_configuration(defaults;
            env = merge(environment, Dict("CHANG_CORE_RADIUS_M" => value)))
    end
    @test_throws ErrorException load_chang_configuration(defaults;
        env = merge(environment, Dict("CHANG_CORE_RADIUS_M" => "nothing")))
    scaled = load_chang_configuration(defaults; env = merge(environment, Dict(
        "CHANG_CORE_RADIUS_M" => "nothing", "CHANG_FCORE_SEGMENT_FACTOR" => "0.1",
        "CHANG_FCORE_CHORD_FACTOR" => "0.02",
    )))
    scaled_law = ChangAeroelastic.chang_aerodynamic_options(scaled).finite_core
    @test isnothing(scaled.aerodynamic.core_radius_m)
    @test scaled_law(2.0, 0.1) == 0.04
    @test scaled_law(0.2, 1.0) == 0.1
    fixed_defaults = merge(defaults, (; aerodynamic = merge(defaults.aerodynamic,
        (; core_radius_m = 0.003))))
    @test load_chang_configuration(fixed_defaults; env = empty_env).aerodynamic.core_radius_m == 0.003

    # Exercise the actual Chang initialization and moving-geometry adapter.
    for ns in (2, 4)
        refined = load_chang_configuration(defaults; env = merge(environment, Dict(
            "CHANG_PROP_RADIAL_PANELS" => string(ns),
            "CHANG_WING_SPAN_PANELS" => string(2ns),
            "CHANG_WING_SYMMETRIC" => string(ns == 2),
        )))
        model = build_chang_model(refined)
        workspace = build_chang_workspace(model)
        @test model.parameters.symmetric_wing == (ns == 2)
        @test !model.parameters.mirror_wing
        blade_surface_count = defaults.propeller.blades *
            length(refined.propeller.attachment_eta)
        @test workspace.system.symmetric == vcat(ns == 2, fill(false, blade_surface_count))
        @test size(workspace.system.surfaces[1]) == (2,2ns)
        cores_match() = all(s -> all(p -> p.core_size == radius, s), workspace.system.surfaces)
        @test cores_match()

        # Image induction must match a separately reflected vortex with the
        # same circulation, for both singular and regularized kernels.
        panel = workspace.system.surfaces[1][1,end]
        reflected = UVLM.reflect(panel)
        rc = SVector(.8,.4,.7)
        for finite_core in (false,true)
            velocity(p; symmetric=false) = first(UVLM.ring_induced_velocity(rc,p;
                symmetric,finite_core))
            @test velocity(panel;symmetric=true) ≈ velocity(panel) + velocity(reflected)
            plane_velocity = first(UVLM.ring_induced_velocity(SVector(.8,0.,.7),panel;
                symmetric=true,finite_core))
            @test abs(plane_velocity[2]) < 1e-12
            wake_panel = UVLM.WakePanel(panel.rtl,panel.rtr,panel.rbl,panel.rbr,radius,1.)
            wake_image = UVLM.WakePanel(reflected.rtl,reflected.rtr,reflected.rbl,reflected.rbr,radius,1.)
            @test velocity(wake_panel;symmetric=true) ≈ velocity(wake_panel) + velocity(wake_image)
        end
        q = fill(1e-5, model.structural.ndof_free)
        update_aero_geometry_for_state!(model, workspace, q, model.parameters.dt[1])
        @test workspace.system.symmetric == vcat(ns == 2, fill(false, blade_surface_count))
        @test cores_match()

        # Shed from every component, then convect a wake panel.
        for surface in workspace.system.surfaces
            nspan = size(surface, 2)
            locations = vcat([surface[end, 1].rbl], [surface[end, j].rbr for j in 1:nspan])
            wake = Matrix{UVLM.WakePanel{Float64}}(undef, 1, nspan)
            velocities = fill(SVector(1.0, 0.0, 0.0), 2, nspan + 1)
            UVLM.shed_wake!(wake, locations, velocities, 0.01, surface, ones(length(surface)), 0)
            @test all(p -> p.core_size == radius, wake)
            convected = UVLM.translate_wake(wake[1, 1], velocities[:, 1:2], 0.01)
            @test convected.core_size == radius
        end
    end

    trim = ChangWindmillingTrimOptions(core_radius_m = radius,
        vortex_core_span_fraction = 0, vortex_core_chord_fraction = 0)
    @test validate_chang_windmilling_options(trim) == 72
    @test chang_trim_finite_core(trim)(0.197, 0.01) == radius
    @test (@inferred chang_trim_finite_core(trim)(0.197, 0.01)) == radius
    @test occursin("fixed", chang_trim_core_description(trim))
    for bad_radius in (0.0, -0.1, NaN, Inf)
        @test_throws ErrorException validate_chang_windmilling_options(
            ChangWindmillingTrimOptions(core_radius_m = bad_radius))
    end
    @test_throws ErrorException validate_chang_windmilling_options(
        ChangWindmillingTrimOptions(vortex_core_span_fraction = 0, vortex_core_chord_fraction = 0))
    @test chang_trim_finite_core(ChangWindmillingTrimOptions())(0.197, 0.1) == 0.05
end

@testset "Direct propeller attachment to a single wing node" begin
    defaults = chang_case_defaults()
    env = Dict(
        "CHANG_WING_SPAN_PANELS" => "8",
        "CHANG_WING_CHORD_PANELS" => "3",
        "CHANG_PROP_RADIAL_PANELS" => "4",
        "CHANG_PROP_CHORD_PANELS" => "2",
        "CHANG_END_TIME_S" => "0.001",
        "CHANG_PLOT_RESULTS" => "false",
        "CHANG_ANIMATE_WAKE" => "false",
    )
    config = load_chang_configuration(defaults; env = env)
    model = build_chang_model(config)
    structural = ChangAeroelastic.assemble_structural_model(model.parameters)

    node = model.parameters.prop_attach_nodes[1]
    target = collect((6 * (node - 1) + 1):(6 * node))
    all_dofs = collect(1:size(structural.attachment_operators[1], 2))
    inactive = setdiff(all_dofs, target)

    @test maximum(abs.(structural.attachment_operators[1][:, inactive])) == 0.0
    @test maximum(abs.(structural.attachment_operators[1][:, target])) > 0.0

    workspace = build_chang_workspace(model)
    @test workspace.attach_node_y[1] == model.parameters.span_nodes[node]
    state = zeros(structural.ndof_free)
    free_node_start = 6 * (node - 2)
    state[free_node_start + 1] = 0.03
    state[free_node_start + 2] = -0.02
    state[free_node_start + 3] = 0.01
    update_aero_geometry_for_state!(model, workspace, state, 0.0)
    expected_node_position = SVector(
        workspace.ea_x_aero[1] - 0.02,
        model.parameters.span_nodes[node] + 0.03,
        -0.01,
    )
    @test workspace.T_pivot_A_current[1] ≈ expected_node_position
    @test workspace.T_hub_A_current[1] == workspace.T_pivot_A_current[1]
    @test workspace.T_load_A_current[1] == workspace.T_pivot_A_current[1]
end
end
