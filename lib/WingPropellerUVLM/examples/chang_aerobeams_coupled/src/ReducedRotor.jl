# The coordinates below retain the Chang reduced flexible-support model. This
# is NOT a geometrically exact pylon beam: its two relative coordinates and
# intrinsic modal inertia are a small-relative-deformation approximation. The
# carrier may nevertheless undergo arbitrary finite translation and rotation.

"""Chang two-coordinate support with an objective moving wing attachment.

`M0` uses carrier coordinates U=[downstream, span, up], with generalized
velocities [v_U; omega_U; theta_dot; psi_dot]. In particular, +psi is the
Chang structural yaw: the relative rotor orientation is Ry(theta)*Rz(-psi).
The prescribed blade spin is Rx(-Omega*t).

The ordinary kinetic model preserves every reference Chang mass block at
eta=0. Its moving centre of mass is reconstructed linearly in eta, while its
intrinsic (centre-of-mass) reduced/modal inertia is frozen at eta=0. This is a
consistent finite-carrier extension of the reduced support, not a replacement
for a nonlinear flexible pylon. The spin angular momentum uses the exact
relative orientation and includes base/base, base/relative and relative/relative
gyroscopic coupling.
"""
struct ReducedRotor{T<:Real}
    M0::Matrix{T}
    K::Matrix{T}
    D::Matrix{T}
    mass::T
    c0::Vector{T}
    cm_modes::Matrix{T}
    h::T
    attachment_node::Int
    reference_index::Int
end

_rr_skew(a) = [zero(a[1]) -a[3] a[2]; a[3] zero(a[1]) -a[1]; -a[2] a[1] zero(a[1])]

"""Build one reduced component from an existing ChangModel without mutating it.

The source structural axes S=[span, downstream, down] map to U through
R_U_S=[0 1 0; 1 0 0; 0 0 -1]. `spin_scale=-1` reverses the structural
spin momentum only; an aerodynamic spin reversal must also be requested by
the caller. Rayleigh damping on the reduced mount is included by default.
"""
function ReducedRotor(reference,ip::Integer; spin_scale::Real=1.0,
    include_stiffness_damping::Bool=true)
    p,s = reference.parameters,reference.structural
    1 <= ip <= p.Npropellers || throw(ArgumentError("invalid propeller index $ip"))
    nodes = s.attach_dofs_all[ip]
    rel = (2ip-1):(2ip)
    # Individual attachment inertia comes from physical parameters rather than
    # Gs_W[nodes,nodes], which can contain another rotor at a coincident node.
    A = zeros(6,6)
    A[1,1] = A[2,2] = A[3,3] = p.mP_prop
    A[3,4] = A[4,3] = p.SαP_prop
    A[1,6] = A[6,1] = -p.SγP_prop
    A[4,4],A[6,6] = p.IαP_prop,p.IγP_prop
    B = [0.0 0.0 p.SθP_prop p.IθαP_prop 0.0 0.0;
        -p.SψP_prop 0.0 0.0 0.0 0.0 p.IψγP_prop]
    M_S = [A B'; B Matrix(s.Ms_P[rel,rel])]
    R_U_S = [0.0 1.0 0.0; 1.0 0.0 0.0; 0.0 0.0 -1.0]
    Q = zeros(8,8)
    Q[1:3,1:3] .= R_U_S
    Q[4:6,4:6] .= R_U_S
    Q[7:8,7:8] .= Matrix{Float64}(I,2,2)
    M0 = Q * M_S * Q'
    K = Matrix(s.Ks_P[rel,rel])
    Cref = Matrix(s.Cs_P[rel,rel])
    D = (Cref + Cref') / 2
    if include_stiffness_damping
        D += s.stiffness_damping_coefficient * K
    end
    mass = Float64(p.mP_prop)
    mass > 0 || throw(ArgumentError("reduced support mass must be positive"))
    # M_vomega=-skew(first moment), M_veta=mass*cm_modes.
    C = M0[1:3,4:6]
    c0 = [C[2,3],C[3,1],C[1,2]] / mass
    cm_modes = M0[1:3,7:8] / mass
    return ReducedRotor(M0,K,D,mass,c0,cm_modes,
        Float64(spin_scale * p.Ix_prop * p.Ω),Int(p.prop_attach_nodes[ip]),Int(ip))
end

"""Mass matrix and its two shape derivatives in the moving carrier frame."""
function rotor_mass_data(rotor::ReducedRotor,eta::AbstractVector)
    length(eta) == 2 || throw(DimensionMismatch("eta must contain pitch and yaw"))
    c = rotor.c0 + rotor.cm_modes * eta
    T = promote_type(eltype(rotor.M0),eltype(eta))
    E3 = Matrix{T}(I,3,3)
    L = hcat(E3,-_rr_skew(c),rotor.cm_modes)
    L0 = hcat(E3,-_rr_skew(rotor.c0),rotor.cm_modes)
    M = rotor.M0 + rotor.mass * (L' * L - L0' * L0)
    derivatives = ntuple(2) do i
        Li = hcat(zeros(T,3,3),-_rr_skew(rotor.cm_modes[:,i]),zeros(T,3,2))
        rotor.mass * (Li' * L + L' * Li)
    end
    return (; M,derivatives,c)
end

"""Exact relative-orientation spin data for Ry(theta)*Rz(-psi), shaft -U-x."""
function rotor_spin_data(rotor::ReducedRotor,eta::AbstractVector)
    theta,psi = eta
    ct,st,cp,sp = cos(theta),sin(theta),cos(psi),sin(psi)
    h = rotor.h
    # H=-h*Ry(theta)*Rz(-psi)*ex. Heta differentiates its carrier components.
    H = [-h * ct * cp,h * sp,h * st * cp]
    Heta = [h*st*cp h*ct*sp;
        zero(h * theta) h*cp;
        h*ct*cp -h*st*sp]
    T = promote_type(eltype(H),eltype(Heta))
    G = zeros(T,8,8)
    G[4:6,4:6] .= -_rr_skew(H)
    G[4:6,7:8] .= Heta
    G[7:8,4:6] .= -Heta'
    # Spin one-form in relative coordinates is [h*sin(psi),0]. Its curl
    # gives the Chang +h/-h pitch-yaw convention at eta=0.
    G[7,8],G[8,7] = h * cp,-h * cp
    return (; H,Heta,G)
end

"""Return reference 8-coordinate mass, gyro, damping and stiffness matrices."""
function rotor_linear_matrices(rotor::ReducedRotor)
    z = zeros(eltype(rotor.M0),2)
    G = rotor_spin_data(rotor,z).G
    D,K = zeros(eltype(rotor.M0),8,8),zeros(eltype(rotor.M0),8,8)
    D[7:8,7:8] .= rotor.D
    K[7:8,7:8] .= rotor.K
    return (; M=copy(rotor.M0),G,D,K)
end

"""Required inertial base wrench and reduced support equation.

Inputs R, v, omega, a, alpha are the base deformation rotation (initially I)
and absolute velocities/accelerations in inertial aerodynamic U axes. eta,
eta_dot and eta_ddot retain the Chang pitch/yaw convention. Returned
`base_wrench` is [force; physical moment] in inertial U axes on the LEFT side
of equilibrium; therefore the component applies its negative to the wing.
`eta_residual` includes inertia, all gyro terms, mount stiffness, physical
and reference Rayleigh damping, and gravity, but excludes applied UVLM loads.

Internal mount springs/dampers appear only in the relative equation in these
unified coordinates. Do not also apply the mount moment or complete hub
aerodynamic wrench a second time as a wing attachment reaction.

The invariant reduced Lagrangian is T=1/2*u'*M(eta)*u with
u=[R'*v; R'*omega; eta_dot]. The base equation uses SE(3) momentum transport;
the relative equation includes d/dt(dT/deta_dot)-dT/deta. Hence moving-base
inertia is retained even when finite wing rotations are present.
"""
function rotor_mechanics(rotor::ReducedRotor,R::AbstractMatrix,
    v::AbstractVector,omega::AbstractVector,a::AbstractVector,alpha::AbstractVector,
    eta::AbstractVector,eta_dot::AbstractVector,eta_ddot::AbstractVector;
    gravity::AbstractVector=zeros(3))
    vb,wb = R' * v,R' * omega
    u = vcat(vb,wb,eta_dot)
    udot = vcat(R' * a - cross(wb,vb),R' * alpha,eta_ddot)
    md = rotor_mass_data(rotor,eta)
    Mdot = md.derivatives[1] * eta_dot[1] + md.derivatives[2] * eta_dot[2]
    momentum = md.M * u
    momentumDot = md.M * udot + Mdot * u
    p,pi = momentum[1:3],momentum[4:6]
    f = momentumDot[1:3] + cross(wb,p)
    moment = momentumDot[4:6] + cross(wb,pi) + cross(vb,p)
    kineticShapeGradient = [dot(u,md.derivatives[i] * u) / 2 for i = 1:2]
    relative = momentumDot[7:8] - kineticShapeGradient +
               rotor.K * eta + rotor.D * eta_dot
    spin = rotor_spin_data(rotor,eta)
    gyro = spin.G * u
    moment += gyro[4:6]
    relative += gyro[7:8]
    gb = R' * gravity
    f -= rotor.mass * gb
    moment -= cross(rotor.mass * md.c,gb)
    relative -= rotor.mass * rotor.cm_modes' * gb
    base_wrench = vcat(R * f,R * moment)
    return (; base_wrench,eta_residual=relative,
        kinetic_energy=dot(u,md.M * u) / 2,
        mount_potential=dot(eta,rotor.K * eta) / 2,
        damping_power=dot(eta_dot,rotor.D * eta_dot),
        gyro_power=dot(u,gyro),body_velocity=u,body_momentum=momentum,
        mass_matrix=md.M,gyro_matrix=spin.G,cm_body=md.c)
end
