using Test
using WingPropellerUVLM

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
    end
end
