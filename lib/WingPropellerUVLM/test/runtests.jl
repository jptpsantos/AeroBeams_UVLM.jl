using Test
using WingPropellerUVLM
using StaticArrays
using LinearAlgebra

include(joinpath(
    @__DIR__,
    "..",
    "examples",
    "chang_linear_aeroelastic",
    "chang_structural_matrices.jl",
))

@testset "WingPropellerUVLM" begin
    @testset "Propeller grid" begin
        twists = zeros(5)
        grids = generate_propeller_blades_grid(1.0,0.1,4,2,twists,3)
        @test length(grids) == 3
        @test all(size(grid) == (3,3,5) for grid in grids)
        @test grids[1][:,1,1] ≈ zeros(3)
    end

    @testset "Reusable interpolation and deformed wing grid" begin
        @test linear_interpolate_1d([0.0, 1.0], [2.0, 4.0], [-1.0, 0.5, 2.0]) ==
            [2.0, 3.0, 4.0]
        @test_throws ArgumentError linear_interpolate_1d([0.0, 0.0], [1.0, 2.0], [0.0])

        grid = generate_panel_grid_and_interpolate(
            1.0,
            [1.0, 1.0],
            [0.0, 0.0],
            1,
            1,
            zeros(2),
            zeros(2),
            zeros(2),
            zeros(2),
            zeros(2),
            zeros(2);
            elastic_axis_fraction = 0.30,
        )
        @test grid[:, 1, 1] ≈ @SVector [0.0, 0.0, 0.0]
        @test grid[:, 2, 2] ≈ @SVector [1.0, 1.0, 0.0]
    end

    @testset "Generalized-alpha integration" begin
        parameters = generalized_alpha_parameters(0.7)
        @test parameters.rho_inf == 0.7
        @test parameters.beta > 0
        @test_throws ArgumentError generalized_alpha_parameters(1.1)

        correction = generalized_alpha_corrector(
            ones(1, 1),
            zeros(1, 1),
            ones(1, 1),
            zeros(1),
            zeros(1),
            zeros(1),
            ones(1),
            zeros(1),
            zeros(1),
            zeros(1),
            0.01,
            parameters,
        )
        @test all(isfinite, correction.displacement)
        @test correction.equilibrium_residual <= eps(Float64)

        options = PartitionedCouplingOptions(
            maximum_iterations = 4,
            state_tolerance = 1.0e-8,
            load_tolerance = 1.0e-8,
            equilibrium_tolerance = 1.0e-12,
            relaxation = 1.0,
        )
        step = partitioned_generalized_alpha_step(
            ones(1, 1),
            zeros(1, 1),
            ones(1, 1),
            zeros(1),
            zeros(1),
            zeros(1),
            zeros(1),
            zeros(1),
            zeros(1),
            0.01,
            parameters,
            _ -> zeros(1);
            options,
            require_load_convergence = false,
        )
        @test step.converged
        @test step.iterations == 1

        last_aerodynamic_state = Ref(NaN)
        coupled_options = PartitionedCouplingOptions(
            maximum_iterations = 20,
            state_tolerance = 1.0e-10,
            load_tolerance = 1.0e-10,
            equilibrium_tolerance = 1.0e-12,
            coupled_equilibrium_tolerance = 1.0e-10,
            relaxation = 1.0,
        )
        coupled_step = partitioned_generalized_alpha_step(
            ones(1, 1),
            zeros(1, 1),
            2.0 .* ones(1, 1),
            zeros(1),
            zeros(1),
            zeros(1),
            zeros(1),
            ones(1),
            zeros(1),
            0.05,
            parameters,
            state -> begin
                last_aerodynamic_state[] = only(state)
                return 0.2 .* state
            end;
            options = coupled_options,
            state_scale = [0.1],
            load_scale = [1.0],
        )
        @test coupled_step.converged
        @test only(coupled_step.displacement) == last_aerodynamic_state[]
        @test coupled_step.trial_load ≈ 0.2 .* coupled_step.displacement
        @test coupled_step.coupled_equilibrium_residual <=
            coupled_options.coupled_equilibrium_tolerance

        failed_step = partitioned_generalized_alpha_step(
            ones(1, 1),
            zeros(1, 1),
            ones(1, 1),
            zeros(1),
            zeros(1),
            zeros(1),
            zeros(1),
            ones(1),
            zeros(1),
            0.05,
            parameters,
            _ -> zeros(1);
            options = PartitionedCouplingOptions(
                maximum_iterations = 1,
                state_tolerance = 1.0e-14,
                load_tolerance = 1.0e-14,
                equilibrium_tolerance = 1.0e-12,
                coupled_equilibrium_tolerance = 1.0e-14,
                relaxation = 1.0,
            ),
            require_load_convergence = false,
        )
        @test !failed_step.converged
        @test_throws DimensionMismatch partitioned_generalized_alpha_step(
            ones(1, 1), zeros(1, 1), ones(1, 1), zeros(1), zeros(1), zeros(1),
            zeros(1), zeros(1), zeros(1), 0.05, parameters, _ -> zeros(1);
            state_scale = ones(2),
        )
    end

    @testset "Chang point-mass parallel-axis correction" begin
        mass = 2.0
        center_of_mass_inertia = Matrix(Diagonal([100.0, 110.0, 120.0]))
        offset = [0.4, -0.2, 0.1]
        skew_offset = chang_tilde(offset)
        corrected = chang_point_mass_matrix(mass, center_of_mass_inertia, offset)
        legacy = chang_point_mass_matrix(
            mass,
            center_of_mass_inertia,
            offset;
            inertia_reference = :beam_axis,
        )

        @test corrected[4:6, 4:6] ≈
            center_of_mass_inertia - mass .* (skew_offset * skew_offset)
        @test legacy[4:6, 4:6] == center_of_mass_inertia
        @test corrected ≈ corrected'

        translational_velocity = [1.2, -0.5, 0.7]
        angular_velocity = [0.1, -0.2, 0.3]
        generalized_velocity = [translational_velocity; angular_velocity]
        center_of_mass_velocity = translational_velocity - skew_offset * angular_velocity
        expected_kinetic_energy = 0.5 * mass * dot(center_of_mass_velocity, center_of_mass_velocity) +
            0.5 * dot(angular_velocity, center_of_mass_inertia * angular_velocity)
        @test 0.5 * dot(generalized_velocity, corrected * generalized_velocity) ≈
            expected_kinetic_energy
        @test_throws ArgumentError chang_point_mass_matrix(
            mass,
            center_of_mass_inertia,
            offset;
            inertia_reference = :unknown,
        )
    end

    @testset "Smooth pulse excitation" begin
        @test smooth_hann_pulse(0.0; start_time = 1.0, duration = 2.0) == 0.0
        @test smooth_hann_pulse(2.0; start_time = 1.0, duration = 2.0) ≈ 1.0
        pulse_load = smooth_hann_pulse_load(
            2.0,
            4,
            [2, 4];
            magnitude = 3.0,
            start_time = 1.0,
            duration = 2.0,
        )
        @test pulse_load == [0.0, 3.0, 0.0, 3.0]
    end

    @testset "Wake-row commit" begin
        active = [0,1,2]
        maximum = [2,1,3]
        @test commit_wake_rows!(active,maximum) == [1,1,3]
    end

    @testset "Imperial College force equations" begin
        @test WingPropellerUVLM._to_imperial_circulation(2.0) == -2.0
        gamma = [2.0 3.0; 5.0 7.0]
        @test WingPropellerUVLM._imperial_spanwise_circulation_jump(gamma,1,1) == -2.0
        @test WingPropellerUVLM._imperial_spanwise_circulation_jump(gamma,2,1) == -3.0
        @test WingPropellerUVLM._imperial_chordwise_circulation_jump(gamma,1,1) == 2.0
        @test WingPropellerUVLM._imperial_chordwise_circulation_jump(gamma,1,2) == 1.0
        @test WingPropellerUVLM._imperial_chordwise_circulation_jump(gamma,1,3) == -3.0

        velocity = @SVector [10.0,0.0,0.0]
        r1 = @SVector [0.0,0.0,0.0]
        r2 = @SVector [0.0,2.0,0.0]
        @test WingPropellerUVLM.imperial_segment_force(velocity,r1,r2,-3.0,1.2) ≈
            @SVector [0.0,0.0,-72.0]

        grid,_ = wing_to_grid(
            [0.0,0.0],[0.0,1.0],[0.0,0.0],ones(2),zeros(2),zeros(2),1,1,
        )
        _,_,surface = grid_to_surface_panels(grid)
        panel = surface[1,1]
        @test WingPropellerUVLM.imperial_panel_area(panel) ≈ 0.75
        @test WingPropellerUVLM.imperial_unsteady_panel_force(panel,4.0,1.2) ≈
            @SVector [0.0,0.0,-3.6]
    end

    @testset "Imperial force-to-node transfer" begin
        span = fill(@SVector([0.0,0.0,0.0]),2,1)
        chord = fill(@SVector([0.0,0.0,0.0]),1,2)
        unsteady = fill(@SVector([0.0,0.0,0.0]),1,1)
        span[1,1] = @SVector [2.0,0.0,0.0]
        chord[1,1] = @SVector [0.0,4.0,0.0]
        chord[1,2] = @SVector [0.0,0.0,6.0]
        unsteady[1,1] = @SVector [8.0,10.0,12.0]

        nodal = imperial_nodal_forces(span,chord,unsteady)
        @test nodal[1,1] ≈ span[1,1]/2 + chord[1,1]/2 + unsteady[1,1]/4
        @test nodal[1,2] ≈ span[1,1]/2 + chord[1,2]/2 + unsteady[1,1]/4
        @test nodal[2,1] ≈ chord[1,1]/2
        @test nodal[2,2] ≈ chord[1,2]/2
        @test sum(nodal) ≈ sum(span) + sum(chord) + unsteady[1,1]/2
    end

    @testset "Imperial vortex-node positions" begin
        grid,_ = wing_to_grid(
            [0.0,0.0],[0.0,2.0],[0.0,0.0],ones(2),zeros(2),zeros(2),2,2,
        )
        _,_,surface = grid_to_surface_panels(grid)
        positions = imperial_nodal_positions(surface)

        @test size(positions) == (size(surface) .+ 1)
        @test positions[1,1] == WingPropellerUVLM.top_left(surface[1,1])
        @test positions[1,end] == WingPropellerUVLM.top_right(surface[1,end])
        @test positions[end,1] == WingPropellerUVLM.bottom_left(surface[end,1])
        @test positions[end,end] == WingPropellerUVLM.bottom_right(surface[end,end])
    end

    @testset "Imperial system-load consistency" begin
        span = 2.0
        chord = 1.0
        speed = 20.0
        density = 1.225
        grid,_ = wing_to_grid(
            [0.0,0.0],[-span/2,span/2],[0.0,0.0],fill(chord,2),
            zeros(2),zeros(2),4,2,
        )
        reference = Reference(span*chord,chord,span,zeros(3),speed,density)
        freestream = Freestream(speed,deg2rad(5.0),0.0,zeros(3))
        system = steady_analysis([grid],reference,freestream;derivatives=false)

        coefficient,_ = body_forces(system)
        dimensional = sum(imperial_nodal_forces(system)[1])
        dynamic_pressure_area = density*speed^2*reference.S/2
        @test dimensional ≈ dynamic_pressure_area*coefficient
        @test coefficient[3] > 0

        # Density is owned by Reference. Circulation and nondimensional
        # coefficients are unchanged when rho changes, while every dimensional
        # Imperial load scales linearly with it. The Trefftz path must use the
        # same reference density so that its coefficient also remains invariant.
        denser_reference = Reference(
            span*chord, chord, span, zeros(3), speed, 2*density,
        )
        denser_system = steady_analysis(
            [grid], denser_reference, freestream; derivatives=false,
        )
        denser_coefficient,_ = body_forces(denser_system)
        denser_dimensional = sum(imperial_nodal_forces(denser_system)[1])
        @test denser_coefficient ≈ coefficient
        @test denser_dimensional ≈ 2*dimensional
        @test far_field_drag(denser_system) ≈ far_field_drag(system)

        @test_throws ArgumentError steady_analysis(
            [grid],reference,freestream;derivatives=true,
        )
    end

    @testset "UVLM snapshot" begin
        xle = [0.0,0.0]
        yle = [0.0,1.0]
        zle = [0.0,0.0]
        chord = [1.0,1.0]
        theta = zeros(2)
        phi = zeros(2)
        grid,_ = wing_to_grid(xle,yle,zle,chord,theta,phi,1,1)
        system = System([grid]; nw=[1])

        _,ratios,surface = grid_to_surface_panels(grid)
        system.ratios[1] = ratios
        system.surfaces[1] .= surface
        system.previous_surfaces[1] .= surface
        system.reference[] = Reference(1.0,1.0,1.0,zeros(3),1.0)
        system.freestream[] = Freestream(1.0,0.0,0.0,zeros(3))
        system.Γ .= 2.0

        snapshot = snapshot_uvlm(system)
        system.Γ .= 3.0
        system.surfaces[1][1] = translate(system.surfaces[1][1],[1.0;0.0;0.0])
        restore_uvlm!(system,snapshot)

        @test system.Γ == snapshot.Γ
        @test system.surfaces == snapshot.surfaces

        activeWakeRows = [0]
        repeatedPoints = repeated_trailing_edge_points(system.surfaces)
        advance_uvlm_trial!(system,snapshot,system.freestream[],0.01;
            repeatedPoints=repeatedPoints,
            activeWakeRows=activeWakeRows,
        )
        ΓFirst = copy(system.Γ)
        wakesFirst = deepcopy(system.wakes)

        advance_uvlm_trial!(system,snapshot,system.freestream[],0.01;
            repeatedPoints=repeatedPoints,
            activeWakeRows=activeWakeRows,
        )
        @test system.Γ ≈ ΓFirst
        @test system.wakes == wakesFirst

        # A mature physical step must be identical whether it is propagated in
        # one call (legacy/default API) or split into a circulation/load trial
        # followed by one accepted wake commit. Start with one initialized wake
        # row so the comparison covers both convection and new-row shedding.
        activeWakeRows .= 1
        acceptedSnapshot = snapshot_uvlm(system)
        advance_uvlm_trial!(
            system,
            acceptedSnapshot,
            system.freestream[],
            0.01;
            repeatedPoints=repeatedPoints,
            activeWakeRows=activeWakeRows,
        )
        ΓFull = copy(system.Γ)
        dΓdtFull = copy(system.dΓdt)
        propertiesFull = deepcopy(system.properties)
        chordForcesFull = deepcopy(system.chord_seg_forces)
        spanForcesFull = deepcopy(system.span_seg_forces)
        unsteadyForcesFull = deepcopy(system.unsteady_forces)
        wakeVelocitiesFull = deepcopy(system.V)
        wakesFull = deepcopy(system.wakes)

        advance_uvlm_trial!(
            system,
            acceptedSnapshot,
            system.freestream[],
            0.01;
            repeatedPoints=repeatedPoints,
            activeWakeRows=activeWakeRows,
            advanceWake=false,
        )
        @test system.Γ == ΓFull
        @test system.dΓdt == dΓdtFull
        @test system.properties == propertiesFull
        @test system.chord_seg_forces == chordForcesFull
        @test system.span_seg_forces == spanForcesFull
        @test system.unsteady_forces == unsteadyForcesFull
        @test system.wakes == acceptedSnapshot.wakes

        advance_wake!(
            system,
            system.freestream[],
            0.01;
            repeated_points=repeatedPoints,
            nwake=activeWakeRows,
        )
        @test system.V == wakeVelocitiesFull
        @test system.wakes == wakesFull
    end
end
