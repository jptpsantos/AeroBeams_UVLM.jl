const PAZY_EXAMPLE_DIR = normpath(joinpath(@__DIR__,".."))

using Test
include(joinpath(PAZY_EXAMPLE_DIR,"PazyWingUVLMStatic.jl"))

# Small meshes test the procedure, not spatial convergence of the benchmark.
options = (; airspeed=40.0,density=1.225,angle_of_attack=deg2rad(5.0),
    chordwise_panels=2,spanwise_panels=4,symmetric_wing=true,maximum_wake_rows=8,
    core_radius=1e-4,wake_shedding_fraction=0.1,aerodynamic_interaction=true,
    maximum_aerodynamic_time=0.5,minimum_aerodynamic_time=0.02,
    aerodynamic_window_time=0.01,aerodynamic_force_tolerance=1e-5,
    aerodynamic_circulation_tolerance=1e-5,maximum_static_iterations=100,
    static_relaxation=0.5,geometry_tolerance=1e-5,load_tolerance=1e-5,
    consecutive_equilibrium_iterations=2,newton_maximum_iterations=50,
    newton_absolute_tolerance=1e-8,newton_relative_tolerance=1e-8,
    newton_display_iterations=false,progress_frequency=10^9)
fast = (; relaxation_mode=:anderson,static_relaxation=1.0,geometry_cap_chords=2.0,wake_reuse=true)

# Equilibria reached by different paths must coincide within the convergence tolerances
same_equilibrium(a,b) = isapprox(a.tip_out_of_plane,b.tip_out_of_plane; rtol=2e-3) &&
    isapprox(a.tip_twist_degrees,b.tip_twist_degrees; rtol=2e-3)

@testset "Accelerated static Pazy UVLM solver" begin
    @testset "Anderson step on a linear fixed-point map" begin
        # G(x) = M x + b with eigenvalues of both signs: plain relaxation with any single
        # factor is slow, the depth-2 Anderson mixing finds the solution of a 2D map exactly
        M = [0.06 0.0; 0.0 -0.78]
        b = [1.0,2.0]
        xStar = (I - M) \ b
        W = ones(2,1)
        x = zeros(2,1); xPrevious = rPrevious = nothing
        dX = Matrix{Float64}[]; dR = Matrix{Float64}[]
        for _ = 1:4
            r = reshape(M * vec(x) .+ b .- vec(x),2,1)
            step = anderson_step!(dX,dR,x,xPrevious,r,rPrevious,1.0,2,W)
            xPrevious,rPrevious = copy(x),r
            x = x .+ step
        end
        @test vec(x) ≈ xStar rtol=1e-6
    end
    @testset "Anderson step with a degenerate history" begin
        # Zero loads (zero incidence) give identical zero residuals: no singular least-squares problem
        x = zeros(6,3); r = zeros(6,3)
        dX = Matrix{Float64}[]; dR = Matrix{Float64}[]
        step = anderson_step!(dX,dR,x,nothing,r,nothing,1.0,4,ones(6,3))
        step = anderson_step!(dX,dR,x,copy(x),r,copy(r),1.0,4,ones(6,3))
        @test all(iszero,step)
        r2 = fill(1e-12,6,3)
        step = anderson_step!(dX,dR,x,copy(x),r2,copy(r),1.0,4,ones(6,3))   # identical x, different r
        @test all(isfinite,step)
    end
    @testset "Zero incidence with the accelerated options" begin
        result = run_pazy_wing_uvlm_static(; options...,fast...,angle_of_attack=0.0)
        @test result.converged
        @test abs(result.tip_out_of_plane) < 1e-10
    end
    reference = run_pazy_wing_uvlm_static(; options...)
    @test reference.converged
    @test reference.tip_out_of_plane > 1e-3
    @testset "Anderson mixing, wake reuse and geometry cap" begin
        result = run_pazy_wing_uvlm_static(; options...,fast...)
        @test result.converged
        @test same_equilibrium(result,reference)
        @test size(result.history,1) < size(reference.history,1)
        @test result.total_aerodynamic_steps < reference.total_aerodynamic_steps
        @test result.geometry_residual <= options.geometry_tolerance
        @test result.load_residual <= options.load_tolerance
        # A very small cap forces many halved load steps and must reach the same equilibrium
        capped = run_pazy_wing_uvlm_static(; options...,fast...,geometry_cap_chords=0.02)
        @test capped.converged
        @test same_equilibrium(capped,reference)
    end
    @testset "Warm start across speeds" begin
        slower = run_pazy_wing_uvlm_static(; options...,fast...,airspeed=30.0)
        U = options.airspeed
        warm = run_pazy_wing_uvlm_static(; options...,fast...,
            warm_start=(; slower.warm_start...,aero=nothing),
            warm_start_target_loads=slower.warm_start.loads .* (U / 30.0)^2)
        @test warm.converged
        @test same_equilibrium(warm,reference)
        @test warm.total_aerodynamic_steps <= run_pazy_wing_uvlm_static(; options...,fast...).total_aerodynamic_steps
    end
    @testset "Warm start across incidences with the same wake" begin
        Δ = deg2rad(0.1)
        central = run_pazy_wing_uvlm_static(; options...,fast...)
        for sign in (-1,1)
            angle = options.angle_of_attack + sign * Δ
            cold = run_pazy_wing_uvlm_static(; options...,angle_of_attack=angle)
            warm = run_pazy_wing_uvlm_static(; options...,fast...,angle_of_attack=angle,
                warm_start=central.warm_start)
            @test warm.converged
            @test same_equilibrium(warm,cold)
            @test warm.total_aerodynamic_steps < cold.total_aerodynamic_steps
            @test warm.aerodynamic_problem === central.aerodynamic_problem  # the wake was reused
            central = warm
        end
    end
    @testset "Options are validated" begin
        @test_throws AssertionError run_pazy_wing_uvlm_static(; options...,relaxation_mode=:aitken)
        @test_throws AssertionError run_pazy_wing_uvlm_static(; options...,geometry_cap_chords=0.0)
    end
end

include(joinpath(PAZY_EXAMPLE_DIR,"PazyWingUVLMStaticSweep.jl"))
sweepOptions = (; density=1.225,chordwise_panels=2,spanwise_panels=3,
    symmetric_wing=true,maximum_wake_rows=6,core_radius=1e-3,
    wake_shedding_fraction=0.1,aerodynamic_interaction=true,
    aerodynamic_force_tolerance=1e-5,aerodynamic_circulation_tolerance=1e-5,
    maximum_static_iterations=65,static_relaxation=0.5,geometry_tolerance=1e-6,
    load_tolerance=1e-5,consecutive_equilibrium_iterations=2,
    newton_maximum_iterations=50,newton_absolute_tolerance=1e-9,
    newton_relative_tolerance=1e-9,newton_display_iterations=false,progress_frequency=10000)

@testset "Static sweep with warm start and accelerated iteration" begin
    directory = mktempdir(joinpath(REPOSITORY_ROOT,"tmp"); cleanup=true)
    common = (; airspeeds=[30.0,40.0],angles_of_attack=[deg2rad(3.0)],sweep_angles=[0.0],
        slope_modes=(:local_ratio,:equilibrium,:fixed_shape),alpha_step=deg2rad(0.1),
        maximum_aerodynamic_chords=120.0,save_plots=false)
    reference = run_pazy_wing_uvlm_static_sweep(; common...,static_options=sweepOptions,
        output_directory=joinpath(directory,"reference"))
    fastOptions = (; sweepOptions...,relaxation_mode=:anderson,static_relaxation=1.0,
        geometry_cap_chords=2.0,wake_reuse=true)
    fast = run_pazy_wing_uvlm_static_sweep(; common...,static_options=fastOptions,warm_start=true,
        output_directory=joinpath(directory,"fast"))
    @test length(fast.cases) == length(reference.cases) == 2
    for (a,b) in zip(fast.cases,reference.cases)
        @test a.central_converged && a.equilibrium_slope_converged && a.fixed_shape_slope_converged
        @test a.section.cn ≈ b.section.cn rtol=2e-3
        @test a.cn_alpha_local ≈ b.cn_alpha_local rtol=2e-3
        @test a.cn_alpha_equilibrium ≈ b.cn_alpha_equilibrium rtol=2e-2
        @test a.cn_alpha_fixed_shape ≈ b.cn_alpha_fixed_shape rtol=2e-2
        @test a.tip_out_of_plane ≈ b.tip_out_of_plane rtol=2e-3
    end
    @test isfile(joinpath(directory,"fast","aoa3deg","spanwise_cn_alpha.csv"))
end

@testset "Sweep warm start is restricted to the same incidence and comparable speed" begin
    directory = mktempdir(joinpath(REPOSITORY_ROOT,"tmp"); cleanup=true)
    fastOptions = (; sweepOptions...,relaxation_mode=:anderson,static_relaxation=1.0,
        geometry_cap_chords=2.0,wake_reuse=true)
    sweep = run_pazy_wing_uvlm_static_sweep(; airspeeds=[0.01,30.0,40.0],
        angles_of_attack=deg2rad.([3.0,5.0]),sweep_angles=[0.0],static_options=fastOptions,
        maximum_aerodynamic_chords=120.0,save_plots=false,warm_start=true,
        output_directory=directory)
    firstTipOOP = [case.diagnostics[1].history[1,7] for case in sweep.cases]
    @test all(case -> case.central_converged,sweep.cases)
    # Cases: (3°: 0.01, 30, 40), (5°: 0.01, 30, 40). Cold starts: first speed of each incidence, and 30 m/s
    # after 0.01 m/s (dynamic pressure ratio 9e6). Warm: 40 m/s after 30 m/s at the same incidence.
    @test firstTipOOP[[1,2,4,5]] == zeros(4)
    @test firstTipOOP[3] > 0 && firstTipOOP[6] > 0
end
