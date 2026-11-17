@testset "Invariance under rigid-body motions and equilibrium of the nodal loads of a doubly-attached spring, for random nodal states and several locations of the spring" begin
    include_example("examples/numeric/doublyAttachedSpringRigidBodyInvariance.jl")
    @test all(ϵKinematics .< 1e-12)
    @test all(ϵLoads .< 1e-12)
    @test all(ϵEquilibrium .< 1e-12)
end

@testset "Analytical Jacobian and inertia matrices of hinged beams, with and without aerodynamic surface, against finite differences" begin
    include_example("examples/numeric/hingedBeamJacobianInertia.jl")
    @test all(jacobianMismatches .== 0)
    @test all(inertiaErrors .< 1e-8)
    @test all(aeroJacobianMismatches .== 0)
end
