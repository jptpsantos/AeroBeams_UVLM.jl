using Test,AeroBeams,LinearAlgebra,SparseArrays
import AeroBeams as AB

function efficiency_structure(steps; distributed=false,tracking=false)
    beam = create_Beam(name="efficiency regression",length=1.0,nElements=2,
        S=[isotropic_stiffness_matrix(EA=1e6,EIy=20.0,EIz=20.0,GJ=10.0)],
        I=[inertia_matrix(ρA=1.0,ρIy=0.01,ρIz=0.01)])
    if distributed
        add_distributed_loads_to_beam!(beam=beam,
            loads=Dict(:f_A_of_x1t => (x,t) -> [1+x+2t,0.0,0.0]))
    end
    root = create_BC(name="clamp",beam=beam,node=1,
        conditions=Dict(type=>0.0 for type in (:u1A,:u2A,:u3A,:p1A,:p2A,:p3A)))
    model = create_Model(name="efficiency regression",beams=[beam],BCs=[root])
    problem = create_DynamicProblem(model=model,timeVector=collect(0:steps) .* 0.01,
        systemSolver=create_NewtonRaphson(),trackingTimeSteps=tracking,displayProgress=false)
    AB.precompute_distributed_loads!(problem)
    return problem
end

@testset "Distributed load storage and analytical nodal resultants" begin
    tables = (:f_A,:m_A,:f_b,:m_b,:ff_A,:mf_A,:ff_b,:mf_b)
    for steps in (2,256)
        problem = efficiency_structure(steps)
        @test all(isempty(getfield(e.distLoads,f)) for
            e in problem.model.assemblyData.elements for f in tables)
    end
    active = efficiency_structure(2; distributed=true)
    for element in active.model.assemblyData.elements,index in eachindex(active.time.timeVector)
        t = active.time.timeVector[index]
        base,ell = 1+element.geometry.x1_n1+2t,element.geometry.Δℓ
        table = element.distLoads.f_A
        @test table[1,1,index] ≈ ell*(base/2+ell/6)
        @test table[2,1,index] ≈ ell*(base/2+ell/3)
        @test iszero(table[:,2:3,index])
        @test all(isempty(getfield(element.distLoads,f)) for f in tables if f != :f_A)
    end
end

@testset "Native sparse factorization reuse" begin
    problem = efficiency_structure(2)
    J = sparse([4.0 1 0;1 5 1;0 1 6])
    rhs = [1.0,-2.0,3.0]
    factor = AB.factorize_jacobian!(problem,J)
    @test norm(J*(factor\rhs)-rhs) < 1e-12
    J[2,2] += 1
    nextFactor = AB.factorize_jacobian!(problem,J)
    @test nextFactor === factor
    @test norm(J*(nextFactor\rhs)-rhs) < 1e-12
    J[1,3] = 0.5
    changedFactor = AB.factorize_jacobian!(problem,J)
    @test norm(J*(changedFactor\rhs)-rhs) < 1e-12
end

@testset "Saved state ownership and quadrature" begin
    problem = efficiency_structure(2; tracking=true)
    element = problem.model.assemblyData.elements[1]
    AB.save_time_step_data!(problem)
    element.state.midpoint.u[1] = 2.0
    AB.update_time_state!(problem,2)
    AB.save_time_step_data!(problem)
    history = problem.solutionHistory.states
    @test history[1][1].midpoint.u[1] == 0.0
    @test history[2][1].midpoint.u[1] == 2.0
    @test history[1][1].midpoint.u !== history[2][1].midpoint.u
    for n in 2:7
        @test [gauss_legendre(x->x^degree,-1.0,1.0; n) for degree in 0:2] ≈ [2.0,0.0,2/3] atol=2e-15
    end
    A = [1.0 2 3;4 5 6;7 8 10]
    b = [0.4,-0.2,0.8]
    @test AB.mul3(A,A',2A,b) ≈ hcat(A*b,A'*b,2A*b)
end
