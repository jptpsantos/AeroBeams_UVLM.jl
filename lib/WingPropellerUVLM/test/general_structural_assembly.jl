module GeneralStructuralAssemblyTests

using Test
using LinearAlgebra
using StaticArrays
using WingPropellerUVLM

const GENERAL_EXAMPLE = joinpath(
    @__DIR__,
    "..",
    "examples",
    "wing_propeller_aeroelastic",
)
include(joinpath(GENERAL_EXAMPLE, "src", "WingPropellerAeroelastic.jl"))
using .WingPropellerAeroelastic

@testset "Bohnisch active consistent mass embedding" begin
    length = 5.7 / 20
    rho = 18.0
    rhoJ = 4.5
    offset = 0.22 * 0.625
    full = consistent_distributed_element_mass(length, rho, rhoJ, offset)
    active = [3, 5, 4, 9, 11, 10]
    reference = rho * length / 420 .* [
        156 22length 0 54 -13length 0;
        22length 4length^2 0 13length -3length^2 0;
        0 0 0 0 0 0;
        54 13length 0 156 -22length 0;
        -13length -3length^2 0 -22length 4length^2 0;
        0 0 0 0 0 0
    ] + length / 420 .* [
        0 0 0 0 0 0;
        0 0 0 0 0 0;
        0 0 140rhoJ 0 0 70rhoJ;
        0 0 0 0 0 0;
        0 0 0 0 0 0;
        0 0 70rhoJ 0 0 140rhoJ
    ] + rho * length * offset .* [
        0 0 7/20 0 0 3/20;
        0 0 length/20 0 0 length/30;
        7/20 length/20 0 3/20 -length/30 0;
        0 0 3/20 0 0 7/20;
        0 0 -length/30 0 0 -length/20;
        3/20 length/30 0 7/20 -length/20 0
    ]
    @test full[active, active] == reference
    @test full[3, 4] != 0.0
    @test full == full'
    @test minimum(eigvals(Symmetric(full))) > 0.0
end

@testset "Direct-node structural attachment and speed laws" begin
    case = bohnisch_case(inactive_stiffness_multiplier = 1.0e3)
    second_propeller = merge(case.propellers[1], (;
        attachment_node = 8,
        radius = 0.7,
        blades = 3,
        blade_chord = 0.16,
        radial_panels = 2,
        chordwise_panels = 2,
        blade_twists = case.propellers[1].blade_twists[[1, 3, end]],
        spin_inertia = 1.7,
        speed_model = :constant_advance_ratio,
        advance_ratio = 1.8,
    ))
    propellers = [case.propellers[1], second_propeller]
    structural = assemble_wing_propeller_structure(
        case.wing,
        propellers,
        case.operating_condition,
    )
    @test structural.attachment_nodes == [21, 8]
    @test case.propellers[1].speed_model == :fixed_rpm
    @test case.propellers[1].rpm == 2500.0
    @test structural.angular_speeds[1] == 2500 * 2pi / 60
    @test structural.angular_speeds[2] ≈
        pi * case.operating_condition.freestream_speed / (1.8 * 0.7)

    pitched = bohnisch_case(
        propeller_radial_panels = 4,
        propeller_rpm = 2100.0,
        propeller_beta75_deg = 12.0,
    )
    @test only(pitched.propellers).speed_model == :fixed_rpm
    @test only(pitched.propellers).rpm == 2100.0
    @test only(pitched.propellers).beta75_deg == 12.0
    @test only(pitched.propellers).blade_twists[4] ≈ deg2rad(12.0)
    @test propeller_angular_speed(
        only(pitched.propellers),
        pitched.operating_condition,
    ) ≈ 2100 * 2pi / 60

    for (index, node) in enumerate(structural.attachment_nodes)
        target = collect(wing_node_dofs(node))
        inactive = setdiff(1:size(structural.attachment_operators[index], 2), target)
        @test all(iszero, structural.attachment_operators[index][:, inactive])
        @test structural.attachment_operators[index][:, target] == Matrix{Float64}(I, 6, 6)
    end

    loads = zeros(size(structural.M_wing, 1))
    wrench = collect(1.0:6.0)
    add_direct_node_wrench!(loads, wrench, 8)
    @test loads[wing_node_dofs(8)] == wrench
    @test all(iszero, loads[wing_node_dofs(7)])
    @test all(iszero, loads[wing_node_dofs(9)])

    force = SVector(1.0, 2.0, 3.0)
    moment = SVector(4.0, 5.0, 6.0)
    point = SVector(0.5, 1.0, -0.2)
    mapped = colocated_hub_wrench(force, moment, point, point)
    @test mapped.force == force
    @test mapped.moment == moment
    @test_throws ArgumentError colocated_hub_wrench(
        force,
        moment,
        point,
        point + SVector(1.0e-3, 0.0, 0.0),
    )

    free_state = zeros(structural.ndof_free)
    node = 8
    prescribed = Dict(1 => 0.2, 2 => 0.1, 3 => 0.05)
    for (local_dof, value) in prescribed
        global_dof = first(wing_node_dofs(node)) + local_dof - 1
        reduced_dof = findfirst(==(global_dof), structural.free_dofs)
        free_state[reduced_dof] = value
    end
    kinematics = full_beam_aerodynamic_kinematics(structural, free_state)
    @test kinematics.u_x[node] == 0.1
    @test kinematics.u_y[node] == 0.2
    @test kinematics.u_z[node] == -0.05
    hubs = direct_propeller_hub_positions(case.wing, propellers, structural, free_state)
    geometry = case.wing.geometry
    fraction = case.wing.span_nodes[node] / last(case.wing.span_nodes)
    chord = geometry.root_chord +
        (geometry.tip_chord - geometry.root_chord) * fraction
    xle = geometry.xle_root + (geometry.xle_tip - geometry.xle_root) * fraction
    @test hubs[2] == SVector(
        xle + geometry.elastic_axis_fraction * chord + 0.1,
        case.wing.span_nodes[node] + 0.2,
        -0.05,
    )

    aerodynamic_case = merge(case, (; propellers))
    uvlm = build_uvlm_system(aerodynamic_case; time_steps = [1.0e-3])
    @test uvlm.blade_counts_prop == [4, 3]
    @test size(uvlm.nodal_forces_prop[2][1]) == (3, 3)
    @test uvlm.attach_node_y == case.wing.span_nodes[[21, 8]]
    @test uvlm.T_pivot_global_init[2] == direct_propeller_hub_positions(
        case.wing,
        propellers,
        structural,
        zeros(structural.ndof_free),
    )[2]
end

@testset "Chang structural regression through general interface" begin
    case = chang_case()
    generalized = build_structural_model(case)
    legacy = WingPropellerAeroelastic.ChangAeroelastic.assemble_structural_model(
        case.legacy_parameters,
    )
    @test generalized.mass_model == :lumped_nodal
    @test generalized.M == legacy.M
    @test generalized.C == legacy.C
    @test generalized.K == legacy.K
    @test structural_natural_frequencies(generalized) ==
        structural_natural_frequencies(legacy)

    remeshed = chang_case(;
        element_count = 10,
        wing_spanwise_panels = 7,
        wing_chordwise_panels = 3,
        propeller_radial_panels = 3,
        propeller_chordwise_panels = 2,
        propeller_attachment_eta = [0.4, 0.8],
    )
    @test getproperty.(remeshed.propellers, :attachment_node) == [5, 9]
    @test remeshed.aerodynamic.spanwise_panels == 7
    @test remeshed.aerodynamic.chordwise_panels == 3
    @test all(propeller -> propeller.radial_panels == 3, remeshed.propellers)
    @test all(propeller -> propeller.chordwise_panels == 2, remeshed.propellers)
end

@testset "Xu reference structural frequencies" begin
    case = xu_case()
    wing = assemble_wing_propeller_structure(
        case.wing,
        NamedTuple[],
        case.operating_condition,
    )
    coupled = build_structural_model(case)
    reference_wing = [2.88, 16.68, 18.06, 46.59, 50.57]
    reference_coupled = [2.85, 5.48, 7.00, 17.75, 19.48, 49.25, 51.37]
    @test structural_natural_frequencies(wing; count = 5) ≈ reference_wing atol = 0.015
    @test structural_natural_frequencies(coupled; count = 7) ≈
        reference_coupled atol = 0.04
    @test case.wing.span_nodes[case.propellers[1].attachment_node] ≈ 1.6

    remeshed = xu_case(
        element_count = 12,
        wing_spanwise_panels = 9,
        wing_chordwise_panels = 3,
        propeller_radial_panels = 4,
        propeller_chordwise_panels = 2,
        propeller_attachment_eta = 0.4,
    )
    expected_node = attachment_node_from_span_fraction(12, 0.4)
    @test remeshed.propellers[1].attachment_node == expected_node
    @test remeshed.wing.span_nodes[expected_node] == 0.4 * 5.7
    @test remeshed.aerodynamic.spanwise_panels == 9
    @test remeshed.aerodynamic.chordwise_panels == 3
    @test remeshed.propellers[1].radial_panels == 4
    @test remeshed.propellers[1].chordwise_panels == 2
    @test length(remeshed.propellers[1].blade_twists) == 5

    @test attachment_node_from_span_fraction(20, 0.42) == 9
    @test attachment_node_from_span_fraction(20, 0.83) == 18
    @test_throws ArgumentError attachment_node_from_span_fraction(0, 0.5)
    @test_throws ArgumentError attachment_node_from_span_fraction(10, 1.1)

    bohnisch_remeshed = bohnisch_case(
        element_count = 12,
        wing_spanwise_panels = 8,
        propeller_attachment_eta = 0.5,
    )
    @test bohnisch_remeshed.propellers[1].attachment_node == 7
    @test bohnisch_remeshed.wing.span_nodes[7] ≈ 2.85
    @test bohnisch_remeshed.aerodynamic.spanwise_panels == 8
end

@testset "Common integrators and coupling schemes" begin
    base = bohnisch_case()
    for integrator in (:newmark_beta, :generalized_alpha)
        for coupling in (:loose_explicit, :implicit_predictor_corrector)
            case = merge(base, (;
                solver_options = (;
                    time_integrator = integrator,
                    coupling_scheme = coupling,
                ),
            ))
            structural = build_structural_model(case)
            result = run_time_domain_analysis(
                case;
                aerodynamic_load = (state, step, time) -> zeros(structural.ndof_free),
                time = [0.0, 1.0e-3],
                time_steps = [1.0e-3],
            )
            @test length(result.displacement_history) == 2
            @test all(iszero, result.displacement_history[end])
        end
    end

    configured = merge(base, (;
        solver_options = (;
            time_integrator = :generalized_alpha,
            generalized_alpha_rho_infinity = 0.7,
            coupling_scheme = :loose_explicit,
            coupling_maximum_iterations = 4,
            coupling_state_tolerance = 2.0e-5,
            coupling_load_tolerance = 3.0e-2,
            coupling_equilibrium_tolerance = 4.0e-10,
            coupling_coupled_equilibrium_tolerance = 5.0e-5,
            coupling_relaxation = 0.6,
        ),
    ))
    configured_structural = build_structural_model(configured)
    configured_result = run_time_domain_analysis(
        configured;
        aerodynamic_load = (state, step, time) -> zeros(configured_structural.ndof_free),
        time = [0.0, 1.0e-3],
        time_steps = [1.0e-3],
    )
    @test configured_result.integration_parameters.rho_infinity == 0.7
    @test configured_result.coupling_options.maximum_iterations == 4
    @test configured_result.coupling_options.state_tolerance == 2.0e-5
    @test configured_result.coupling_options.load_tolerance == 3.0e-2
    @test configured_result.coupling_options.equilibrium_tolerance == 4.0e-10
    @test configured_result.coupling_options.coupled_equilibrium_tolerance == 5.0e-5
    @test configured_result.coupling_options.relaxation == 0.6
end

@testset "Common full UVLM time-domain path" begin
    base = bohnisch_case()
    propeller = merge(only(base.propellers), (;
        blades = 2,
        radial_panels = 1,
        chordwise_panels = 1,
        blade_twists = only(base.propellers).blade_twists[[1, end]],
    ))
    case = merge(base, (;
        propellers = [propeller],
        aerodynamic = (;
            spanwise_panels = 2,
            chordwise_panels = 1,
            mirror = false,
            symmetric = false,
        ),
        solver_options = (;
            time_integrator = :newmark_beta,
            coupling_scheme = :loose_explicit,
        ),
    ))
    result = run_uvlm_time_domain_analysis(
        case;
        time = [0.0, 1.0e-3],
        time_steps = [1.0e-3],
        maximum_wake_rows_wing = 1,
        maximum_wake_rows_propeller = 1,
    )
    @test all(isfinite, result.displacement_history[end])
    @test result.uvlm.iwake == result.uvlm.nwake
    accepted_hub = only(result.aerodynamic.last_kinematics[].hubs)
    direct_hub = only(direct_propeller_hub_positions(
        case.wing,
        case.propellers,
        result.structural,
        result.displacement_history[end],
    ))
    @test accepted_hub == direct_hub
end

end
