using Test, AeroBeams, LinearAlgebra, SparseArrays, DelimitedFiles
import AeroBeams as AB

function efficiency_structure(steps; distributed=false)
    beam = create_Beam(name="efficiency regression",length=1.0,nElements=2,
        S=[isotropic_stiffness_matrix(EA=1e6,EIy=20.0,EIz=20.0,GJ=10.0)],
        I=[inertia_matrix(ρA=1.0,ρIy=0.01,ρIz=0.01)])
    distributed && add_loads_to_beam!(beam;loadTypes=["f_A_of_x1t"],
        loadFuns=[(x,t)->[1+x+2t,0.0,0.0]])
    root = create_BC(name="clamp",beam=beam,node=1,
        types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=zeros(6))
    model = create_Model(name="efficiency regression",beams=[beam],BCs=[root])
    problem = create_DynamicProblem(model=model,timeVector=collect(0:steps) .* 0.01,
        trackingTimeSteps=false,displayProgress=false)
    AB.precompute_distributed_loads!(problem)
    return problem
end

@testset "Distributed load storage, interpolation and Newton checkpoints" begin
    short = efficiency_structure(2)
    long = efficiency_structure(256)
    @test all(isempty(getfield(e,f)) for e in long.model.elements for f in AB.DISTRIBUTED_LOAD_TABLES)
    @test sum(sizeof(getfield(e,f)) for e in short.model.elements for f in AB.DISTRIBUTED_LOAD_TABLES) ==
        sum(sizeof(getfield(e,f)) for e in long.model.elements for f in AB.DISTRIBUTED_LOAD_TABLES) == 0
    active = efficiency_structure(2;distributed=true)
    for e in active.model.elements, index in eachindex(active.timeVector)
        t = active.timeVector[index]
        base, ell = 1+e.x1_n1+2t,e.Δℓ
        @test e.f_A[1,1,index] ≈ ell*(base/2+ell/6)
        @test e.f_A[2,1,index] ≈ ell*(base/2+ell/3)
        @test iszero(e.f_A[:,2:3,index])
        @test all(isempty(getfield(e,f)) for f in AB.DISTRIBUTED_LOAD_TABLES if f != :f_A)
    end
    AB.update_time_variables!(active,2)
    for t in (active.timeBeginTimeStep,0.005,active.timeEndTimeStep), factor in (0.0,0.3,1.0)
        active.timeNow = t
        array = active.model.elements[1].f_A
        reference = AB.interpolate_distributed_loads(active,factor .* array)
        result = AB.interpolate_distributed_loads(active,array,factor)
        @test all(isapprox(a,b;atol=1e-14) for (a,b) in zip(reference,result))
        saved = copy(array)
        result[1][1] += 10
        @test array == saved
    end
    snapshot = AB.checkpoint_model(active.model)
    @test snapshot.elements[1].f_A === active.model.elements[1].f_A
    @test snapshot.R_A_ofTime === active.model.R_A_ofTime
    @test snapshot.elements[1] === snapshot.beams[1].elements[1]
    @test snapshot.BCs[1].beam === snapshot.beams[1]
    @test snapshot.elements[1].states.u !== active.model.elements[1].states.u
    @test snapshot.BCs[1].currentValue !== active.model.BCs[1].currentValue
    active.model.elements[1].states.u[1] = 2.0
    active.model.BCs[1].currentValue[1] = 3.0
    @test iszero(snapshot.elements[1].states.u[1])
    @test iszero(snapshot.BCs[1].currentValue[1])
end

@testset "Cached sparse Newton increments" begin
    J = sparse([4.0 1 0;1 5 1;0 1 6])
    rhs = [1.0,-2.0,3.0]
    step = similar(rhs)
    cache = AB.NewtonLinearCache()
    AB.solve_newton_increment!(step,J,rhs,cache,false)
    @test norm(J*step+rhs) < 1e-12
    factor = cache.factorization
    rhs .*= 2
    AB.solve_newton_increment!(step,J,rhs,cache,true)
    @test cache.factorization === factor
    @test norm(J*step+rhs) < 1e-12
    J[2,2] += 1
    AB.solve_newton_increment!(step,J,rhs,cache,false)
    @test cache.factorization === factor
    @test norm(J*step+rhs) < 1e-12
    J[1,3] = 0.5 # Change the pattern: symbolic reuse must be invalidated.
    AB.solve_newton_increment!(step,J,rhs,cache,false)
    @test norm(J*step+rhs) < 1e-12
end

@testset "Saved state ownership and fixed quadrature" begin
    problem = efficiency_structure(2)
    AB.save_time_step_data!(problem,0.0)
    problem.x[1] = 2.0
    AB.save_time_step_data!(problem,0.01)
    @test problem.xOverTime[1][1] == 0.0
    @test problem.xOverTime[2][1] == 2.0
    @test problem.xOverTime[1] !== problem.xOverTime[2]
    steady = create_SteadyProblem(model=problem.model)
    AB.save_load_factor_data!(steady,0.5,steady.x)
    steady.x[1] = 3.0
    @test steady.xOverσ[1][1] == 0.0
    for n in 2:7
        @test gauss_legendre(x->[1.0,x,x^2],-1.0,1.0;n) ≈ [2.0,0.0,2/3] atol=2e-15
    end
    A = [1.0 2 3;4 5 6;7 8 10]
    b = [0.4,-0.2,0.8]
    @test AB.mul3(A,A',2A,b) ≈ hcat(A*b,A'*b,2A*b)
end

@testset "Existing analytical spring case and modified Newton" begin
    scope = Module(:EfficiencySpringCase)
    Base.include(scope,joinpath(@__DIR__,"examples","axialDoublyAttachedSpringCantilever.jl"))
    @test scope.u1[end] ≈ scope.uTipAnalytical rtol=1e-10
    @test norm(scope.spring.Fs) ≈ scope.FsAnalytical rtol=1e-10
    p = scope.problem
    p.systemSolver.alwaysUpdateJacobian = false
    AB.assemble_system_arrays!(p)
    J = copy(p.jacobian)
    p.skipJacobianUpdate = true
    AB.assemble_system_arrays!(p)
    @test p.jacobian == J
    p.skipJacobianUpdate = false
    AB.solve_NewtonRaphson!(p)
    @test p.systemSolver.convergedFinalSolution
    @test p.model.elements[end].nodalStates.u_n2[1] ≈ scope.uTipAnalytical rtol=1e-9
end

@testset "Strip aerodynamics with AD and tip quadrature" begin
    _,_,chord,spar = geometrical_properties_Pazy()
    model,_ = create_Pazy(aeroSolver=QuasiSteady(),derivationMethod=AD(),
        airfoil=deepcopy(flatPlate),upright=true,θ=3π/180,airspeed=30.0,g=0,
        tipMass=0.01,ηtipMass=[0.0,-chord*(1-spar)-0.04,0.0],smallAngles=true)
    p = create_SteadyProblem(model=model,systemSolver=create_NewtonRaphson(displayStatus=false))
    solve!(p)
    @test p.systemSolver.convergedFinalSolution
    tip = p.nodalStatesOverσ[end][end]
    R,_ = rotation_tensor_WM(tip.p_n2_b)
    outputs = (tip_OOP=tip.u_n2_b[3],tip_IP=-tip.u_n2_b[2],tip_twist=asind((R*[0.0,1.0,0.0])[3]))
    for (name,value) in pairs(outputs)
        reference = readdlm(joinpath(@__DIR__,"newTestDataGenerators","PazyWingPitchRange","$name.txt"))[1,31]
        @test value ≈ reference atol=1e-4 rtol=1e-4
    end
end
