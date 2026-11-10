# Pure mechanical checks; no aerodynamic package loads or output files.
# Run: julia --compiled-modules=existing --project=. lib/WingPropellerUVLM/examples/chang_aerobeams_coupled/test_tip_rotor.jl
using Test, LinearAlgebra, ForwardDiff
include(joinpath(@__DIR__,"src","ReducedRotor.jl"))
include(joinpath(@__DIR__,"src","TipRotor.jl"))

reference = (;
    parameters=(;Npropellers=1,Ix_prop=2.6,Ω=120.0,mP_prop=100.0),
    config=(;propeller=(;blades=4),structural=(;blade_mass_kg=1.44)),
)
tip = TipRotor(reference,1,9)
z = zeros(3)
empty_eta = Float64[]
E3 = Matrix{Float64}(I,3,3)
mechanics(rotor,R,v,w,a,alpha;kwargs...) = rotor_mechanics(
    rotor,R,v,w,a,alpha,empty_eta,empty_eta,empty_eta;kwargs...)

@testset "Physical rotor mass and inertia ownership" begin
    @test tip.mass == 4*1.44
    @test tip.mass != reference.parameters.mP_prop
    @test tip.matrix_J ≈ Diagonal([2.6,1.3,1.3])
    @test tip.h == 2.6*120
    @test tip.attachment_node == 9
    @test tip.reference_index == 1
    custom = TipRotor(reference,1,7;rotor_mass=8.2,spin_inertia=1.8,
        inertia_tensor=Diagonal([3.0,4.0,5.0]))
    @test custom.mass == 8.2
    @test custom.h == 1.8*120
    @test custom.matrix_J ≈ Diagonal([3.,4.,5.])
    @test_throws ArgumentError TipRotor(reference,1,9;inertia_tensor=Diagonal([1.,1.,4.]))
    @test_throws ArgumentError TipRotor(reference,1,9;rotor_mass=-1)
    @test_throws DimensionMismatch rotor_mechanics(tip,E3,z,z,z,z,[0.],empty_eta,empty_eta)
    @test eigmin(Symmetric(rotor_linear_matrices(tip).M)) > 0
end

@testset "Full physical inertia and spin tangent" begin
    function response_from_state(acc,vel)
        r = mechanics(tip,E3,vel[1:3],vel[4:6],acc[1:3],acc[4:6])
        @test isempty(r.eta_residual)
        r.base_wrench
    end
    Ma = ForwardDiff.jacobian(a->response_from_state(a,zeros(6)),zeros(6))
    Gv = ForwardDiff.jacobian(v->response_from_state(zeros(6),v),zeros(6))
    linear = rotor_linear_matrices(tip)
    @test Ma ≈ linear.M atol=1e-12
    @test Gv ≈ linear.G atol=1e-12
    @test linear.G+linear.G' ≈ zeros(6,6) atol=1e-12
    @test linear.G[5:6,5:6] ≈ [0. -tip.h;tip.h 0.]
    @test rotor_spin_data(tip).H ≈ [-tip.h,0.,0.]
    reversed = TipRotor(reference,1,9;spin_scale=-1)
    @test rotor_linear_matrices(reversed).G ≈ -linear.G
    axial = mechanics(tip,E3,z,[0.7,0.,0.],z,z)
    @test axial.base_wrench ≈ zeros(6) atol=1e-12
end

@testset "Finite tip angular momentum, energy and objectivity" begin
    # A measured rotor/nacelle inertia need not be axisymmetric in carrier axes.
    rotor = TipRotor(reference,1,9;rotor_mass=7.0,spin_inertia=2.6,
        inertia_tensor=[3.0 0.1 0.2;0.1 4.0 -0.1;0.2 -0.1 5.0])
    ex,ey = [1.,0.,0.],[0.,1.,0.]
    function trajectory(t)
        f,fd,fdd = 0.2t+0.15t^2,0.2+0.3t,0.3
        g,gd,gdd = -0.4t+0.1t^2,-0.4+0.2t,0.2
        R1 = exp(_rr_skew(ex)*f)
        R = R1*exp(_rr_skew(ey)*g)
        w = fd*ex+gd*R1*ey
        alpha = fdd*ex+gdd*R1*ey+gd*fd*cross(ex,R1*ey)
        v,a = [cos(t),2t,0.9t^2],[-sin(t),2.,1.8t]
        H = R*rotor.matrix_J*R'*w+R*[-rotor.h,0.,0.]
        E = rotor.mass*dot(v,v)/2+dot(w,R*rotor.matrix_J*R'*w)/2
        (;R,w,alpha,v,a,H,E)
    end
    t,dt = 0.7,1e-5
    q,qp,qm = trajectory(t),trajectory(t+dt),trajectory(t-dt)
    r = mechanics(rotor,q.R,q.v,q.w,q.a,q.alpha)
    @test r.base_wrench[4:6] ≈ (qp.H-qm.H)/(2dt) rtol=1e-9 atol=1e-8
    @test r.base_wrench[1:3] ≈ rotor.mass*(qp.v-qm.v)/(2dt) rtol=1e-9 atol=1e-8
    power = dot(q.v,r.base_wrench[1:3])+dot(q.w,r.base_wrench[4:6])
    @test power ≈ (qp.E-qm.E)/(2dt) rtol=1e-9 atol=1e-8
    @test abs(r.gyro_power) < 1e-12
    Q = exp(_rr_skew([0.3,-0.5,0.7]))
    rotated = mechanics(rotor,Q*q.R,Q*q.v,Q*q.w,Q*q.a,Q*q.alpha)
    @test rotated.base_wrench ≈ vcat(Q*r.base_wrench[1:3],Q*r.base_wrench[4:6]) atol=1e-10
    boosted = mechanics(rotor,q.R,q.v+[4.,2.,-7.],q.w,q.a,q.alpha)
    @test boosted.base_wrench ≈ r.base_wrench atol=1e-12
    gravity = [0.,0.,-9.81]
    grav = mechanics(rotor,q.R,q.v,q.w,q.a,q.alpha;gravity)
    @test grav.base_wrench-r.base_wrench ≈ vcat(-rotor.mass*gravity,zeros(3)) atol=1e-12
end
