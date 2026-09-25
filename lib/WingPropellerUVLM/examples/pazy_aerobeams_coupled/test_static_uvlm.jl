using Test
include(joinpath(@__DIR__, "PazyWingUVLMStatic.jl"))

# Small meshes test the procedure, not spatial convergence of the benchmark.
options = (; airspeed=15.0, density=1.225, angle_of_attack=deg2rad(3.0),
    chordwise_panels=2, spanwise_panels=3, symmetric_wing=true, maximum_wake_rows=6,
    core_radius=1e-3, wake_shedding_fraction=0.1, aerodynamic_interaction=true,
    maximum_aerodynamic_time=0.3, minimum_aerodynamic_time=0.02,
    aerodynamic_window_time=0.01, aerodynamic_force_tolerance=1e-4,
    aerodynamic_circulation_tolerance=1e-4, maximum_static_iterations=60,
    static_relaxation=0.5, geometry_tolerance=1e-5, load_tolerance=1e-4,
    consecutive_equilibrium_iterations=2, newton_maximum_iterations=50,
    newton_absolute_tolerance=1e-8, newton_relative_tolerance=1e-8,
    newton_display_iterations=false, progress_frequency=1000)

@testset "Static Pazy time-marched UVLM" begin
    @testset "Zero incidence" begin
        result = run_pazy_wing_uvlm_static(; options..., angle_of_attack=0.0)
        @test result.converged
        @test result.structural_problem isa AeroBeams.SteadyProblem
        @test norm(result.tip_displacement_A) < 1e-10
        @test abs(result.tip_twist_degrees) < 1e-10
        @test result.dt ≈ 0.0989/(2*15)
        @test size(result.history,2) == length(result.history_columns)
        @test all(diff(result.history[:,2]) .> 0)
    end
    @testset "Nonzero static deflection" begin
        result = run_pazy_wing_uvlm_static(; options...)
        @test result.converged
        # Positive UVLM lift is upward for the Pazy model (-A-x).
        @test result.tip_out_of_plane > 1e-5
        @test all(isfinite, result.tip_displacement_A)
        @test isfinite(result.tip_twist_degrees)
        @test result.geometry_residual <= options.geometry_tolerance
        @test result.load_residual <= options.load_tolerance
        @test all(iszero(v) for v in result.aerodynamic_problem.state.system.Vcp[1])
        @test all(result.history[:,14] .== 1)
        @test all(bc.currentValue ≈ result.applied_nodal_loads[:,n]
            for (n,bc) in enumerate(result.structural_problem.model.BCs[2:end]))
        @test length(result.structural_problem.nodalStatesOverσ) == 1
        # An independent geometric reconstruction: final wake and static shape agree.
        ne, _, chord, spar = AeroBeams.geometrical_properties_Pazy()
        eta = AeroBeams.nodal_positions_Pazy()
        W = zeros(options.spanwise_panels+1,ne+1)
        for j in axes(W,1)
            station = (j-1)/options.spanwise_panels
            left = min(searchsortedlast(eta,station),ne)
            f = (station-eta[left])/(eta[left+1]-eta[left])
            W[j,left], W[j,left+1] = 1-f, f
        end
        grid, _, _ = wing_geometry(result.structural_problem.model,chord,spar,
            options.chordwise_panels,W)
        @test grid ≈ result.aerodynamic_problem.state.system.grids[1]
    end
    @testset "Limits must not report false equilibrium" begin
        wake_limit = run_pazy_wing_uvlm_static(; options..., maximum_aerodynamic_time=0.004)
        @test !wake_limit.converged
        @test occursin("Wake did not settle",wake_limit.termination_reason)
        iteration_limit = run_pazy_wing_uvlm_static(; options..., maximum_static_iterations=1)
        @test !iteration_limit.converged
        @test occursin("Maximum static iterations",iteration_limit.termination_reason)
    end
end
