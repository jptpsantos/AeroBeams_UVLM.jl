using Test
include(joinpath(@__DIR__,"..","PazyWingUVLMStaticSweep.jl"))

# Small meshes check load conservation, sign, slope definitions and wake
# convergence gates; they do not establish spatial convergence of Pazy data.
options = (; density=1.225,chordwise_panels=2,spanwise_panels=3,
    symmetric_wing=true,maximum_wake_rows=6,core_radius=1e-3,
    wake_shedding_fraction=0.1,aerodynamic_interaction=true,
    aerodynamic_force_tolerance=1e-5,aerodynamic_circulation_tolerance=1e-5,
    maximum_static_iterations=65,static_relaxation=0.5,geometry_tolerance=1e-6,
    load_tolerance=1e-5,consecutive_equilibrium_iterations=2,
    newton_maximum_iterations=50,newton_absolute_tolerance=1e-9,
    newton_relative_tolerance=1e-9,newton_display_iterations=false,progress_frequency=10000)

# Prescribed streamwise sections with independently known incidence and force.
function local_incidence_problem(betaNodes; alpha=deg2rad(3.0),ratios=[4.0,5.0,6.0])
    ns,nc = length(betaNodes) - 1,2
    span,chord,U,rho = 0.6,0.1,15.0,1.225
    grid = zeros(3,nc + 1,ns + 1)
    for j = 1:ns+1,i = 1:nc+1
        distance = ((i - 1) / nc - 0.5) * chord
        grid[:,i,j] = [distance * cos(betaNodes[j]),(j - 1) * span / ns,distance * sin(betaNodes[j])]
    end
    zeroForce = UVLM.SVector(0.0,0.0,0.0)
    spanForces = fill(zeroForce,nc + 1,ns)
    for j = 1:ns
        beta = (betaNodes[j] + betaNodes[j+1]) / 2
        normal = UVLM.SVector(-sin(beta),0.0,cos(beta))
        spanForces[1,j] = (rho * U^2 / 2) * (chord * span / ns) * ratios[j] * (alpha - beta) * normal
    end
    _,_,surface = UVLM.grid_to_surface_panels(grid)
    system = (; grids=[grid],surfaces=[surface],span_seg_forces=[spanForces],
        chord_seg_forces=[fill(zeroForce,nc,ns + 1)],unsteady_forces=[fill(zeroForce,nc,ns)],
        freestream=Ref(UVLM.Freestream(U,alpha,0.0,zeros(3))))
    reference = UVLM.Reference(span * chord,chord,span,zeros(3),U,rho)
    return (; state=(; system),model=(; reference)),span
end

@testset "Static spanwise cn_alpha" begin
    @testset "Local incidence and secant definition" begin
        betaNodes = deg2rad.([-3.0,-1.0,1.0,2.0])
        problem,span = local_incidence_problem(betaNodes)
        section = pazy_spanwise_cn(problem,span,0.0)
        @test section.alpha_local ≈ deg2rad.([6.0,4.0,2.0,1.0]) atol=1e-14
        betaStrips = (betaNodes[1:end-1] + betaNodes[2:end]) / 2
        stripCn = [4.0,5.0,6.0] .* (deg2rad(3.0) .- betaStrips)
        expectedCn = [stripCn[1] * cos(betaStrips[1] - betaNodes[1]),
            (stripCn[1] * cos(betaStrips[1] - betaNodes[2]) + stripCn[2] * cos(betaStrips[2] - betaNodes[2])) / 2,
            (stripCn[2] * cos(betaStrips[2] - betaNodes[3]) + stripCn[3] * cos(betaStrips[3] - betaNodes[3])) / 2,
            stripCn[3] * cos(betaStrips[3] - betaNodes[4])]
        @test section.cn ≈ expectedCn atol=1e-14
        @test section.cn_alpha_local ≈ expectedCn ./ deg2rad.([6.0,4.0,2.0,1.0]) atol=1e-12
        @test all(section.local_alpha_valid)
        flat,span = local_incidence_problem(zeros(4); alpha=deg2rad(-2.0))
        @test pazy_spanwise_cn(flat,span,0.0).alpha_local ≈ fill(deg2rad(-2.0),4)
        uniform,span = local_incidence_problem(zeros(4); ratios=fill(5.0,3))
        @test pazy_spanwise_cn(uniform,span,0.0).cn_alpha_local ≈ fill(5.0,4)
        for alpha in (0.0,1e-9)
            zeroIncidence,span = local_incidence_problem(zeros(4); alpha)
            section = pazy_spanwise_cn(zeroIncidence,span,0.0)
            @test !any(section.local_alpha_valid)
            @test all(isnan,section.cn_alpha_local)
        end
    end
    @testset "Pitching moment about physical chord references" begin
        problem,span = local_incidence_problem(zeros(4); ratios=fill(5.0,3))
        section = pazy_spanwise_cn(problem,span,0.0; spar_fraction=0.4)
        # With two chordwise panels, the first bound vortex is at x/c=0.125.
        @test section.cm_quarter_chord ≈ 0.125 .* section.cn
        @test section.cm_spar ≈ 0.275 .* section.cn
        # Move those same forces aft to the trailing-edge vortex vertices.
        forces = problem.state.system.span_seg_forces[1]
        forces[end,:] .= forces[1,:]
        forces[1,:] .= Ref(zero(eltype(forces)))
        aft = pazy_spanwise_cn(problem,span,0.0; spar_fraction=0.4)
        @test aft.cn == section.cn
        @test aft.cm_quarter_chord ≈ -0.75 .* aft.cn
        @test aft.cm_spar ≈ -0.6 .* aft.cn
    end
    @testset "Nodal force redistribution and tributary areas" begin
        problem,span = local_incidence_problem(zeros(4))
        system = problem.state.system
        # Flat sections make all three force components independently
        # checkable, including shared segments and residual unsteady loads.
        fill!(system.span_seg_forces[1],UVLM.SVector(2.0,4.0,6.0))
        fill!(system.chord_seg_forces[1],UVLM.SVector(3.0,5.0,7.0))
        fill!(system.unsteady_forces[1],UVLM.SVector(4.0,6.0,8.0))
        section = pazy_spanwise_cn(problem,span,0.0)
        expectedForces = [13.0 20.0 20.0 13.0; 22.0 34.0 34.0 22.0; 31.0 48.0 48.0 31.0]
        @test section.forces_UVLM ≈ expectedForces
        @test section.eta == [0.0,1 / 3,2 / 3,1.0]
        @test section.strip_area_m2 ≈ [0.01,0.02,0.02,0.01]
        @test sum(section.strip_area_m2) ≈ problem.model.reference.S
        q = problem.model.reference.rho * problem.model.reference.V^2 / 2
        @test section.cn ≈ expectedForces[3,:] ./ (q .* section.strip_area_m2)
        @test sum(q .* section.strip_area_m2 .* section.cn) ≈ sum(expectedForces[3,:])
        @test section.positions_UVLM[2,:] ≈ [0.0,0.2,0.4,0.6]
        # A single spanwise panel has two half-width endpoint sections.
        single,span = local_incidence_problem(zeros(2); ratios=[5.0])
        section = pazy_spanwise_cn(single,span,0.0)
        @test section.eta == [0.0,1.0]
        @test section.strip_area_m2 ≈ fill(single.model.reference.S / 2,2)
        @test section.cn_alpha_local ≈ [5.0,5.0]
        # The production 20-panel case must include all 21 stations.
        fine,span = local_incidence_problem(zeros(21); ratios=fill(5.0,20))
        section = pazy_spanwise_cn(fine,span,0.0)
        @test section.eta ≈ collect(0:20) ./ 20
        @test section.cn_alpha_local ≈ fill(5.0,21)
    end
    @testset "Local incidence changes in central differences" begin
        minus = (; cn=[0.1,0.2,0.3],alpha_local=[0.01,0.02,0.03],
            cm_quarter_chord=zeros(3),cm_spar=zeros(3))
        plus = (; minus...,cn=[0.15,0.24,0.35],alpha_local=[0.02,0.01,0.03 + 1e-10])
        slopes = pazy_static_derivatives((minus,plus),0.01)
        @test slopes.cn ≈ [2.5,2.0,2.5]
        @test slopes.local_incidence.valid == [true,true,false]
        @test slopes.local_incidence.cn_alpha[1:2] ≈ [5.0,-4.0]
        @test isnan(slopes.local_incidence.cn_alpha[3])
    end
    @testset "Validation" begin
        @test_throws ArgumentError run_pazy_wing_uvlm_static_sweep(;
            airspeeds=[0.0],angles_of_attack=[0.0],static_options=options)
        @test_throws ArgumentError run_pazy_wing_uvlm_static_sweep(;
            airspeeds=[15.0],angles_of_attack=[0.0],alpha_step=0.0,
            slope_modes=(:equilibrium,),static_options=options)
        @test_throws ArgumentError run_pazy_wing_uvlm_static_sweep(;
            airspeeds=[15.0],angles_of_attack=[0.0],slope_modes=(:bad,),static_options=options)
        @test_throws ArgumentError run_pazy_wing_uvlm_static_sweep(;
            airspeeds=[15.0],angles_of_attack=deg2rad.([1e6]),static_options=options)
        @test_throws ArgumentError run_pazy_wing_uvlm_static_sweep(;
            airspeeds=[15.0],angles_of_attack=deg2rad.([89.9]),
            slope_modes=(:fixed_shape,),static_options=options)
        @test_throws AssertionError run_pazy_wing_uvlm_static_sweep(;
            airspeeds=[15.0],angles_of_attack=[0.0],sweep_angles=deg2rad.([5.0]),static_options=options)
    end
    @testset "Conservation on a deformed swept wing" begin
        equilibrium = run_pazy_wing_uvlm_static(; options...,airspeed=15.0,
            angle_of_attack=deg2rad(2.0),sweepAngle=deg2rad(20.0),
            maximum_aerodynamic_time=0.4,minimum_aerodynamic_time=0.02,aerodynamic_window_time=0.02)
        @test equilibrium.converged
        interface = pazy_static_interface(equilibrium)
        problem = equilibrium.aerodynamic_problem
        section = pazy_spanwise_cn(problem,interface.span,deg2rad(20.0))
        @test interface.grid ≈ problem.state.system.grids[1]
        @test all(>(0),section.cn)
        @test vec(sum(section.forces_UVLM; dims=2)) ≈ sum(UVLM.dimensional_loads(problem).forces[1])
        @test sum(section.strip_area_m2) ≈ problem.model.reference.S
        @test all(norm(section.normals_UVLM[:,j]) ≈ 1 for j in eachindex(section.cn))
        @test section.eta == [0.0,1 / 3,2 / 3,1.0]
        @test section.strip_area_m2[[1,end]] ≈ fill(problem.model.reference.S / 6,2)
        @test section.y_reference_m ≈ section.s_m .* cosd(20)
        @test all(section.local_alpha_valid)
        @test section.cn_alpha_local ≈ section.cn ./ section.alpha_local
    end
    @testset "Single equilibrium without perturbations" begin
        directory = mktempdir(joinpath(REPOSITORY_ROOT,"tmp"); cleanup=!("--plots" in ARGS))
        result = run_pazy_wing_uvlm_static_sweep(; airspeeds=[15.0],angles_of_attack=deg2rad.([3.0]),
            sweep_angles=deg2rad.([10.0]),static_options=options,alpha_step=0.0,
            maximum_aerodynamic_chords=80.0,output_directory=directory,save_plots="--plots" in ARGS)
        case = only(result.cases)
        @test result.slope_modes == (:local_ratio,)
        @test case.central_converged
        @test only(case.diagnostics).solve == :central
        @test case.angle_of_attack == deg2rad(3.0)
        @test isnan(case.alpha_step)
        @test all(case.cn_alpha_local_valid)
        @test all(>(0),case.cn_alpha_local)
        @test case.cn_alpha_local ≈ case.section.cn ./ case.section.alpha_local
        @test maximum(case.section.alpha_local) - minimum(case.section.alpha_local) > 1e-6
        @test all(isnan,case.cn_alpha_equilibrium)
        @test all(isnan,case.cn_alpha_fixed_shape)
        angleDirectory = joinpath(directory,"aoa3deg")
        table = readdlm(joinpath(angleDirectory,"spanwise_cn_alpha.csv"),',',Any)
        @test size(table) == (5,23)
        @test table[2,6] == 0.0 && table[end,6] == 1.0
        settings = read(joinpath(angleDirectory,"settings.txt"),String)
        @test occursin("half-strip averages",settings)
        @test table[1,20:23] == ["alpha_local_deg","alpha_local_rad","cn_alpha_local_per_rad","cn_alpha_local_valid"]
        moments = readdlm(joinpath(angleDirectory,"spanwise_cm.csv"),',',Any)
        @test isequal(moments[:,1:8],table[:,1:8])
        @test moments[2:end,9] == case.section.cm_quarter_chord
        @test moments[2:end,10] == case.section.cm_spar
        @test all(isnan,moments[2:end,11:14])
        compact = readdlm(joinpath(angleDirectory,"spanwise_cn_alpha_clean.csv"),',',Any)
        @test size(compact) == (5,12)
        @test compact[2:end,8] == case.cn_alpha_local
        @test all(isempty,compact[2:end,[5,9,10]])
        @test compact[2:end,11:12] == moments[2:end,9:10]
        @test size(readdlm(joinpath(angleDirectory,"spanwise_cn_alpha_perturbations.csv"),',',Any),1) == 1
        @test size(readdlm(joinpath(angleDirectory,"status.csv"),',',Any),1) == 2
        println("Local-incidence equilibrium test output: ",directory)
    end
    @testset "Speed scaling and perturbation refinement" begin
        directory = mktempdir(joinpath(REPOSITORY_ROOT,"tmp"); cleanup=!("--plots" in ARGS))
        result = run_pazy_wing_uvlm_static_sweep(; airspeeds=[10.0,20.0],angles_of_attack=[0.0],
            sweep_angles=deg2rad.([0.0,20.0]),static_options=options,slope_modes=(:equilibrium,:fixed_shape),
            maximum_aerodynamic_chords=80.0,output_directory=directory,save_plots="--plots" in ARGS)
        @test length(result.cases) == 4
        for case in result.cases
            @test case.central_converged
            @test case.equilibrium_slope_converged
            @test case.fixed_shape_slope_converged
            @test maximum(abs,case.section.cn) < 1e-10
            @test all(>(0),case.cn_alpha_equilibrium)
            @test all(>(0),case.cn_alpha_fixed_shape)
            @test all(<(2pi),case.cn_alpha_fixed_shape)
            @test maximum(abs,case.section.cm_spar) < 1e-10
            @test all(isfinite,case.cm_quarter_chord_alpha_fixed_shape)
            @test all(isfinite,case.cm_spar_alpha_fixed_shape)
            fixedLocal = case.cn_alpha_local_perturbations.fixed_shape
            @test all(fixedLocal.valid)
            @test fixedLocal.delta_alpha_rad ≈ fill(2 * case.alpha_step,length(case.section.eta))
            @test fixedLocal.cn_alpha ≈ case.cn_alpha_fixed_shape
            @test all(case.cn_alpha_local_perturbations.equilibrium.valid)
        end
        for i in (1,3)
            # At zero incidence the central wing stays flat. Convective wake
            # times remove U-dependence of this incompressible fixed-shape slope.
            @test result.cases[i].cn_alpha_fixed_shape ≈ result.cases[i+1].cn_alpha_fixed_shape rtol=1e-8
        end
        @test !isapprox(result.cases[1].cn_alpha_equilibrium,result.cases[2].cn_alpha_equilibrium; rtol=1e-3)
        finer = run_pazy_wing_uvlm_static_sweep(; airspeeds=[10.0],angles_of_attack=[0.0],
            static_options=options,slope_modes=:fixed_shape,alpha_step=deg2rad(0.125),
            maximum_aerodynamic_chords=80.0)
        @test finer.slope_modes == (:fixed_shape,)
        @test finer.cases[1].fixed_shape_slope_converged
        @test finer.cases[1].cn_alpha_fixed_shape ≈ result.cases[1].cn_alpha_fixed_shape rtol=0.01
        @test finer.cases[1].cm_spar_alpha_fixed_shape ≈ result.cases[1].cm_spar_alpha_fixed_shape rtol=0.01
        @test all(isnan,finer.cases[1].cn_alpha_equilibrium)
        finerEquilibrium = run_pazy_wing_uvlm_static_sweep(; airspeeds=[10.0],angles_of_attack=[0.0],
            static_options=options,slope_modes=(:equilibrium,),alpha_step=deg2rad(0.125),maximum_aerodynamic_chords=80.0)
        @test finerEquilibrium.cases[1].equilibrium_slope_converged
        @test finerEquilibrium.cases[1].cn_alpha_equilibrium ≈ result.cases[1].cn_alpha_equilibrium rtol=0.01
        @test all(isnan,finerEquilibrium.cases[1].cn_alpha_fixed_shape)
        angleDirectory = joinpath(directory,"aoa0deg")
        table = readdlm(joinpath(angleDirectory,"spanwise_cn_alpha.csv"),',',Any)
        @test size(table) == (17,23)
        @test "cn_alpha_fixed_shape_per_rad" in table[1,:]
        @test isfile(joinpath(angleDirectory,"settings.txt"))
        @test isfile(joinpath(angleDirectory,"case_1_central_convergence.csv"))
        @test size(readdlm(joinpath(angleDirectory,"status.csv"),',',Any),1) == 21
        perturbations = readdlm(joinpath(angleDirectory,"spanwise_cn_alpha_perturbations.csv"),',',Any)
        @test size(perturbations) == (33,13)
        @test perturbations[1,9:end] == ["mode","delta_cn","delta_alpha_local_rad",
            "cn_alpha_local_perturbation_per_rad","valid"]
        @test perturbations[2:3,9] == ["equilibrium","fixed_shape"]
        @test Float64.(perturbations[2:3,12]) ≈ [
            result.cases[1].cn_alpha_local_perturbations.equilibrium.cn_alpha[1],
            result.cases[1].cn_alpha_fixed_shape[1]]
        compact = readdlm(joinpath(angleDirectory,"spanwise_cn_alpha_clean.csv"),',',Any)
        @test size(compact) == (17,12)
        @test all(isempty,compact[2:end,8])
        @test compact[2:5,9] == result.cases[1].cn_alpha_local_perturbations.fixed_shape.cn_alpha
        @test compact[2:5,10] == result.cases[1].cn_alpha_local_perturbations.fixed_shape.delta_alpha_rad
        if "--plots" in ARGS
            @test count(file -> startswith(file,"spanwise_cm_") && endswith(file,".png"),readdir(angleDirectory)) == 2
        end
        println("Spanwise sweep test output: ",directory)
    end
    @testset "Small-incidence local ratio and pitching moment" begin
        directory = mktempdir(joinpath(REPOSITORY_ROOT,"tmp"))
        alpha = deg2rad(1e-6)
        result = run_pazy_wing_uvlm_static_sweep(; airspeeds=[60.0],
            angles_of_attack=[alpha,100alpha],local_alpha_tolerance=1e-10,
            static_options=(; options...,load_residual_force_floor=1e-9),
            maximum_aerodynamic_chords=80.0,output_directory=directory,save_plots=false)
        proxy,check = result.cases
        @test all(case -> case.central_converged && all(case.cn_alpha_local_valid),result.cases)
        @test only(proxy.diagnostics).solve == :central
        @test proxy.cn_alpha_local ≈ check.cn_alpha_local rtol=0.01
        @test 100 .* proxy.section.cm_spar ≈ check.section.cm_spar rtol=0.01
        momentTable = readdlm(joinpath(directory,pazy_static_angle_folder(alpha),"spanwise_cm.csv"),',',Any)
        @test momentTable[2,3] ≈ 1e-6
        @test momentTable[2:end,10] == proxy.section.cm_spar
    end
    @testset "Convergence gates" begin
        directory = mktempdir(joinpath(REPOSITORY_ROOT,"tmp"))
        failed = @test_logs (:warn,r"Unconverged static/wake solves") run_pazy_wing_uvlm_static_sweep(;
            airspeeds=[15.0],angles_of_attack=[0.0],static_options=options,maximum_aerodynamic_chords=0.5,
            output_directory=directory,save_plots=false)
        @test !failed.cases[1].central_converged
        @test all(isnan,failed.cases[1].cn_alpha_equilibrium)
        @test all(isnan,failed.cases[1].cn_alpha_fixed_shape)
        @test all(isnan,failed.cases[1].cn_alpha_local)
        @test !any(failed.cases[1].cn_alpha_local_valid)
        moments = readdlm(joinpath(directory,"aoa0deg","spanwise_cm.csv"),',',Any)
        @test all(isnan,moments[2:end,9:14])
        compact = readdlm(joinpath(directory,"aoa0deg","spanwise_cn_alpha_clean.csv"),',',Any)
        @test all(isempty,compact[2:end,5:12])
        # Zero-incidence wakes settle exactly; nonzero-incidence
        # perturbations fail an impossible tolerance.
        failedPerturbation = @test_logs (:warn,r"Unconverged static/wake solves") run_pazy_wing_uvlm_static_sweep(;
            airspeeds=[15.0],angles_of_attack=[0.0],slope_modes=(:fixed_shape,),
            static_options=(; options...,aerodynamic_force_tolerance=1e-30,
                aerodynamic_circulation_tolerance=1e-30),maximum_aerodynamic_chords=15.0,
            minimum_aerodynamic_chords=6.0,output_directory=directory,save_plots=false)
        @test failedPerturbation.cases[1].central_converged
        @test !failedPerturbation.cases[1].fixed_shape_slope_converged
        @test isnothing(failedPerturbation.cases[1].cn_alpha_local_perturbations.fixed_shape)
        @test all(isnan,failedPerturbation.cases[1].cm_spar_alpha_fixed_shape)
        perturbations = readdlm(joinpath(directory,"aoa0deg","spanwise_cn_alpha_perturbations.csv"),',',Any)
        @test all(isnan,perturbations[2:end,10:12])
        @test all(value -> string(value) == "false",perturbations[2:end,13])
        @test all(isnan,failedPerturbation.cases[1].cn_alpha_fixed_shape)
    end
end
