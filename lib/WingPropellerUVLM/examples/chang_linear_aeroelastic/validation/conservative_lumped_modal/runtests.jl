using Test
using LinearAlgebra

include(joinpath(@__DIR__, "ConservativeLumpedModal.jl"))
using .ConservativeLumpedModal

@testset "conservative lumped spatial-inertia remap" begin
    reference_total = reduce(+, source_spatial_inertia_blocks())
    source_model = assemble_wing_model(16; mesh = :source, remap = :source_lumps)
    @test source_model.blocks[2:end] == source_spatial_inertia_blocks()
    for ne in (16, 20, 30, 40, 80)
        wing = assemble_wing_model(ne)
        audit = remap_audit(wing)
        properties = remapped_nodal_properties(ne)
        target_total = reduce(+, wing.blocks)
        @test isapprox(target_total, reference_total; rtol = 1e-12, atol = 1e-12)
        @test audit.mass_minimum_eigenvalue > 0
        @test audit.block_minimum_eigenvalue > 0
        @test audit.stiffness_minimum_eigenvalue > 0
        @test all(properties.mass[2:end] .> 0)
        @test all(isfinite, properties.cg_x)
        @test all(
            minimum(eigvals(Symmetric([
                properties.Ixx[i] properties.Ixy[i] properties.Ixz[i];
                properties.Ixy[i] properties.Iyy[i] properties.Iyz[i];
                properties.Ixz[i] properties.Iyz[i] properties.Izz[i]
            ]))) > 0 for i in 2:length(properties.mass)
        )
        @test all(isfinite(mode.frequency_hz) for mode in modal_analysis(
            wing; number_of_modes = 8, classify_wing = true,
        ))
    end
end

@testset "20-element failure reproduction" begin
    failing = assemble_wing_model(20; remap = :signed_scaled)
    audit = remap_audit(failing)
    @test audit.minimum_cg_inertia_eigenvalue < 0
    @test audit.mass_minimum_eigenvalue < 0
    @test audit.maximum_offset > 1.0
end

@testset "coupled system" begin
    wing = assemble_wing_model(20)
    coupled = assemble_coupled_model(wing)
    @test minimum(eigvals(Symmetric(coupled.M))) > 0
    @test minimum(eigvals(Symmetric(coupled.K))) > 0
    modes = modal_analysis(coupled; number_of_modes = 10)
    @test length(modes) == 10
    @test issorted([mode.frequency_hz for mode in modes])
end
