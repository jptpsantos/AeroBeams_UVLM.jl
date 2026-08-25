using Test
using WingPropellerUVLM
using StaticArrays

@testset "WingPropellerUVLM" begin
    @testset "Propeller grid" begin
        twists = zeros(5)
        grids = generate_propeller_blades_grid(1.0,0.1,4,2,twists,3)
        @test length(grids) == 3
        @test all(size(grid) == (3,3,5) for grid in grids)
        @test grids[1][:,1,1] ≈ zeros(3)
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
    end
end
