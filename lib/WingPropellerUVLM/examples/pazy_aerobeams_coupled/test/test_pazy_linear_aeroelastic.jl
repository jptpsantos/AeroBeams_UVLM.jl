using Test
using LinearAlgebra
include(joinpath(@__DIR__,"..","PazyWingUVLMStatic.jl"))
include(joinpath(@__DIR__,"..","PazyLinearAeroelastic.jl"))

@testset "Coupled eigenvalue continuation utilities" begin
    values = reshape(ComplexF64[-1 + 2im,1 + 2im,-1 + 2im],3,1)
    crossings = PazyLinearAeroelastic.growth_rate_crossings([10.0,20.0,30.0],values)
    @test getproperty.(crossings,:transition) == [:onset,:restabilization]
    @test getproperty.(crossings,:estimate_mps) ≈ [15.0,25.0]
    @test_throws DimensionMismatch PazyLinearAeroelastic.growth_rate_crossings(
        [10.0,20.0],values)

    previous_values = ComplexF64[-1 + 2im,-1 - 2im]
    current_values = reverse(previous_values)
    previous_vectors = Matrix{ComplexF64}(I,2,2)
    current_vectors = previous_vectors[:,[2,1]]
    @test PazyLinearAeroelastic.match_eigenmodes(previous_values,previous_vectors,
        current_values,current_vectors) == [2,1]
end

@testset "Pazy coupled frozen-wake eigenanalysis" begin
    options = (;airspeed=15.0,density=1e-8,angle_of_attack=deg2rad(0.6),
        chordwise_panels=2,spanwise_panels=3,symmetric_wing=true,maximum_wake_rows=4,
        core_radius=1e-3,wake_shedding_fraction=0.1,aerodynamic_interaction=true,
        maximum_aerodynamic_time=0.2,minimum_aerodynamic_time=0.0,
        aerodynamic_window_time=0.01,aerodynamic_force_tolerance=5e-4,
        aerodynamic_circulation_tolerance=5e-4,maximum_static_iterations=10,
        static_relaxation=0.5,geometry_tolerance=1e-5,load_tolerance=5e-4,
        consecutive_equilibrium_iterations=2,newton_maximum_iterations=50,
        newton_absolute_tolerance=1e-8,newton_relative_tolerance=1e-8,
        newton_display_iterations=false,progress_frequency=1000,
        mesh=:nominal,withWingtipRod=false)
    reference = run_pazy_wing_uvlm_static(;options...)
    @test reference.converged
    speeds = [15.0,16.0]
    sweep = PazyLinearAeroelastic.sweep_pazy_linear_flutter(speeds,
        speed -> run_pazy_wing_uvlm_static(;merge(options,(;airspeed=speed))...);
        number_of_modes=1,finite_difference_step=1e-6)
    result = first(sweep.cases).result
    @test result.analysis == :coupled_frozen_wake
    @test size(result.coupled_state_matrix,1) ==
          size(result.uvlm.A,1) + 2 * length(result.structural_eigenvalues)
    @test all(isfinite,result.coupled_state_matrix)
    @test all(isfinite,result.discrete_eigenvalues)
    @test all(isfinite,result.growth_rates_per_s)
    for structural_eigenvalue in result.structural_eigenvalues
        distance = minimum(abs.(result.continuous_eigenvalues .- structural_eigenvalue))
        @test distance <= 1e-3 * max(abs(structural_eigenvalue),1.0)
    end
    @test getproperty.(sweep.cases,:speed_mps) == speeds
    @test size(sweep.tracked_continuous_eigenvalues,1) == length(speeds)
    mktempdir() do directory
        path = PazyLinearAeroelastic.save_pazy_coupled_modes(
            joinpath(directory,"modes.csv"),result)
        @test isfile(path)
        @test occursin("frequency_hz",read(path,String))
        sweepPath = PazyLinearAeroelastic.save_pazy_flutter_sweep(
            joinpath(directory,"sweep.csv"),sweep)
        @test isfile(sweepPath)
        @test length(readlines(sweepPath)) > 2
        boundaryPath = PazyLinearAeroelastic.save_pazy_flutter_boundaries(
            joinpath(directory,"boundaries.csv"),sweep)
        @test isfile(boundaryPath)
        @test startswith(read(boundaryPath,String),"mode,transition,")
    end
end