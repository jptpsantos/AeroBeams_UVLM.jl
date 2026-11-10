# Standalone mechanical checks; no aerodynamic time stepping or files written.
# Run: julia --project=. lib/WingPropellerUVLM/examples/chang_aerobeams_coupled/test_reduced_rotor.jl
using Test, LinearAlgebra, ForwardDiff
push!(LOAD_PATH,normpath(joinpath(@__DIR__,"..","..")))
include(joinpath(@__DIR__,"..","chang_linear_aeroelastic","src","ChangAeroelastic.jl"))
using .ChangAeroelastic
include(joinpath(@__DIR__,"src","ReducedRotor.jl"))

config = load_chang_configuration(;env=Dict{String,String}())
reference = build_chang_model(config)
rotor = ReducedRotor(reference,1)
E3 = Matrix{Float64}(I,3,3)
z2,z3 = zeros(2),zeros(3)

function mechanics_from_quasi(rotor,u,acc;eta=zeros(2))
    rotor_mechanics(rotor,E3,u[1:3],u[4:6],acc[1:3],acc[4:6],
        eta,u[7:8],acc[7:8])
end

@testset "Reference Chang component blocks" begin
    s,p = reference.structural,reference.parameters
    node,rel = s.attach_dofs_all[1],1:2
    RUS = [0.0 1 0; 1 0 0; 0 0 -1]
    Q = zeros(8,8)
    Q[1:3,1:3] .= RUS
    Q[4:6,4:6] .= RUS
    Q[7:8,7:8] .= Matrix{Float64}(I,2,2)
    Mref = [s.Gs_W[node,node] s.Fs_W[node,rel];
            s.Bs_P[rel,node] s.Ms_P[rel,rel]]
    CP = s.Cs_P[rel,rel]
    Gref = [s.Cs_W[node,node] s.Hs_W[node,rel];
            s.Ds_P[rel,node] (CP-CP')/2]
    linear = rotor_linear_matrices(rotor)
    @test linear.M ≈ Q*Mref*Q' atol=1e-12
    @test linear.G ≈ Q*Gref*Q' atol=1e-12
    @test linear.G+linear.G' ≈ zeros(8,8) atol=1e-13
    Mfd = ForwardDiff.jacobian(zeros(8)) do acc
        m = mechanics_from_quasi(rotor,zeros(8),acc)
        vcat(m.base_wrench,m.eta_residual)
    end
    Cfd = ForwardDiff.jacobian(zeros(8)) do u
        m = mechanics_from_quasi(rotor,u,zeros(8))
        vcat(m.base_wrench,m.eta_residual)
    end
    @test Mfd ≈ linear.M atol=1e-11
    @test Cfd ≈ linear.G+linear.D atol=1e-11
    @test rotor_linear_matrices(ReducedRotor(reference,1;spin_scale=-1)).G ≈ -linear.G
    for eta in ([0.,0.],[0.12,-0.09],[-0.08,0.13])
        M = rotor_mass_data(rotor,eta).M
        @test isapprox(M,M';atol=1e-13)
        @test eigmin(Symmetric(M)) >= -1e-10
        @test eigmin(Symmetric(M[7:8,7:8])) > 0
    end
end

@testset "Finite carrier objectivity, Galilean invariance and power" begin
    rot_y(a) = [cos(a) 0 sin(a); 0 1 0; -sin(a) 0 cos(a)]
    rot_z(a) = [cos(a) -sin(a) 0; sin(a) cos(a) 0; 0 0 1]
    R,Q = rot_z(0.41)*rot_y(-0.28),rot_y(0.77)*rot_z(-0.19)
    v,w,a,alpha = [2.,-1.,0.4],[0.3,-0.7,0.6],[-0.4,0.8,0.2],[0.1,0.2,-0.5]
    eta,edot,eddot = [0.07,-0.09],[0.4,-0.2],[-0.5,0.3]
    r = rotor_mechanics(rotor,R,v,w,a,alpha,eta,edot,eddot)
    rotated = rotor_mechanics(rotor,Q*R,Q*v,Q*w,Q*a,Q*alpha,eta,edot,eddot)
    @test rotated.base_wrench ≈ vcat(Q*r.base_wrench[1:3],Q*r.base_wrench[4:6]) atol=1e-10
    @test rotated.eta_residual ≈ r.eta_residual atol=1e-10
    boost = rotor_mechanics(rotor,R,v+[5.,-3.,8.],w,a,alpha,eta,edot,eddot)
    @test boost.base_wrench ≈ r.base_wrench atol=1e-10
    @test boost.eta_residual ≈ r.eta_residual atol=1e-10
    @test abs(r.gyro_power) < 1e-11
    vb,wb = R'*v,R'*w
    u = vcat(vb,wb,edot)
    udot = vcat(R'*a-cross(wb,vb),R'*alpha,eddot)
    md = rotor_mass_data(rotor,eta)
    Mdot = md.derivatives[1]*edot[1]+md.derivatives[2]*edot[2]
    kinetic_rate = dot(u,md.M*udot)+dot(u,Mdot*u)/2
    power = dot(v,r.base_wrench[1:3])+dot(w,r.base_wrench[4:6])+
        dot(edot,r.eta_residual)-r.damping_power-dot(edot,rotor.K*eta)
    @test power ≈ kinetic_rate atol=1e-10
end

@testset "Dry two-coordinate gyro frequencies" begin
    M,K = rotor.M0[7:8,7:8],rotor.K
    G = rotor_linear_matrices(rotor).G[7:8,7:8]
    A = [zeros(2,2) Matrix{Float64}(I,2,2); -M\K -M\G]
    computed = sort(imag.(filter(x->imag(x)>0,eigvals(A))))
    j1,j2,k1,k2 = M[1,1],M[2,2],K[1,1],K[2,2]
    b = k1*j2+k2*j1+rotor.h^2
    disc = sqrt(b^2-4j1*j2*k1*k2)
    analytic = sqrt.([(b-disc)/(2j1*j2),(b+disc)/(2j1*j2)])
    @test computed ≈ analytic rtol=1e-12
    @test maximum(abs.(real.(eigvals(A)))) < 1e-10
end
