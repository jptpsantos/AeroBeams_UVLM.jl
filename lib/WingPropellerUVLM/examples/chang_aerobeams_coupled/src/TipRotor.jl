"""Rigid axisymmetric-mean spinning rotor at a nonlinear pylon tip.

The rotor centre of mass and hub coincide with `attachment_node`. `matrix_J`
is the ordinary inertia tensor about that centre, expressed in the initially
aerodynamic U axes [downstream, span, up]. `h=spin_inertia*Omega` is signed;
the initial spin momentum vector is [-h,0,0], consistent with Rx(-Omega*t).

All rotor/nacelle mass and ordinary inertia belong to this nodal component.
The pylon beam supplies its own distributed mass. Do not add this component's
mass to the pylon section inertia or add a duplicate AeroBeams PointInertia.
There are no independent relative pitch/yaw coordinates in this variation:
the rotor carrier follows the finite pylon-tip orientation.
"""
struct TipRotor{T<:Real}
    mass::T
    matrix_J::Matrix{T}
    spin_inertia::T
    h::T
    attachment_node::Int
    reference_index::Int
end

"""Construct the physical rotor-only tip component from a Chang reference.

Default mass is blades*blade_mass_kg, excluding all reduced-pylon mass. The
default spin inertia is the reference calibrated Ix_prop. With no measured
transverse inertia, Jtransverse=Ix_prop/2 assumes an axisymmetric thin disk;
this is an explicit inertia assumption, not an additional Chang datum.

`rotor_mass` and `inertia_tensor` can describe a measured rotor/nacelle assembly
about the hub/CM. `spin_inertia` remains the independently spinning rotor's
polar inertia, even when stationary nacelle inertia is included in matrix_J.
`axial_inertia` and `transverse_inertia` offer diagonal alternatives to the
full `inertia_tensor`. `transverse_inertia` may be one scalar or two values.
"""
function TipRotor(reference,ip::Integer,node::Integer;
    transverse_inertia=nothing,axial_inertia=nothing,rotor_mass=nothing,
    spin_inertia=nothing,inertia_tensor=nothing,spin_scale::Real=1.0)
    p,c = reference.parameters,reference.config
    1 <= ip <= p.Npropellers || throw(ArgumentError("invalid propeller index $ip"))
    node >= 1 || throw(ArgumentError("attachment node must be positive"))
    m = Float64(isnothing(rotor_mass) ?
        c.propeller.blades*c.structural.blade_mass_kg : rotor_mass)
    Js = Float64(isnothing(spin_inertia) ? p.Ix_prop : spin_inertia)
    isfinite(m) && m > 0 || throw(ArgumentError("rotor/nacelle mass must be positive and finite"))
    isfinite(Js) && Js > 0 || throw(ArgumentError("spin inertia must be positive and finite"))
    if isnothing(inertia_tensor)
        Ja = isnothing(axial_inertia) ? Js : Float64(axial_inertia)
        Jt = isnothing(transverse_inertia) ? (Js/2,Js/2) :
            transverse_inertia isa Real ? (transverse_inertia,transverse_inertia) : transverse_inertia
        length(Jt) == 2 || throw(DimensionMismatch("provide one transverse inertia or a pair"))
        J = Matrix(Diagonal(Float64[Ja,Jt[1],Jt[2]]))
    else
        isnothing(axial_inertia) && isnothing(transverse_inertia) ||
            throw(ArgumentError("use inertia_tensor or diagonal inertia options, not both"))
        size(inertia_tensor) == (3,3) || throw(DimensionMismatch("inertia tensor must be 3 by 3"))
        J = Matrix{Float64}(inertia_tensor)
    end
    all(isfinite,J) || throw(ArgumentError("inertia tensor must be finite"))
    isapprox(J,J';rtol=1e-12,atol=0.0) || throw(ArgumentError("inertia tensor must be symmetric"))
    principal = eigvals(Symmetric(J))
    minimum(principal) > 0 || throw(ArgumentError("inertia tensor must be positive definite"))
    maximum(principal) <= sum(principal)/2+1e-12*maximum(principal) ||
        throw(ArgumentError("principal inertias violate the rigid-body triangle inequality"))
    h = Float64(spin_scale*Js*p.Ω)
    isfinite(h) || throw(ArgumentError("spin angular momentum must be finite"))
    return TipRotor(m,J,Js,h,Int(node),Int(ip))
end

function _tip_rotor_empty_coordinates(eta,eta_dot=eta,eta_ddot=eta)
    isempty(eta) && isempty(eta_dot) && isempty(eta_ddot) ||
        throw(DimensionMismatch("TipRotor follows its pylon tip and has no relative coordinates"))
    return nothing
end

"""Six-coordinate ordinary inertia at the physical hub/CM, in carrier axes."""
function rotor_mass_data(rotor::TipRotor,eta::AbstractVector=Float64[])
    _tip_rotor_empty_coordinates(eta)
    T = eltype(rotor.matrix_J)
    M = zeros(T,6,6)
    M[1:3,1:3] .= rotor.mass*Matrix{T}(I,3,3)
    M[4:6,4:6] .= rotor.matrix_J
    return (; M,derivatives=(),c=zeros(T,3))
end

"""Spin momentum and its skew gyro matrix in the pylon-tip carrier axes."""
function rotor_spin_data(rotor::TipRotor,eta::AbstractVector=Float64[])
    _tip_rotor_empty_coordinates(eta)
    H = [-rotor.h,zero(rotor.h),zero(rotor.h)]
    G = zeros(eltype(H),6,6)
    G[4:6,4:6] .= -_rr_skew(H)
    return (; H,Heta=zeros(eltype(H),3,0),G)
end

function rotor_linear_matrices(rotor::TipRotor)
    M = rotor_mass_data(rotor).M
    G = rotor_spin_data(rotor).G
    return (; M,G,D=zeros(eltype(M),6,6),K=zeros(eltype(M),6,6))
end

"""Required physical tip wrench for prescribed constant rotor spin.

`R` is the finite tip deformation rotation, initially I in U coordinates;
v, omega, a and alpha are absolute inertial-U velocities and accelerations.
The component's ordinary inertia and spin momentum are rotated by R. The
returned `base_wrench` is on the left side of equilibrium, so its negative
is the rotor reaction applied to the pylon. eta arguments must be empty.

The rotational residual is J_world*alpha + omega x (J_world*omega+H_spin).
This includes ordinary rigid-body inertia and prescribed-spin gyroscopic
transport once each. Gravity acts at the coincident hub and CM.
"""
function rotor_mechanics(rotor::TipRotor,R::AbstractMatrix,
    v::AbstractVector,omega::AbstractVector,a::AbstractVector,alpha::AbstractVector,
    eta::AbstractVector,eta_dot::AbstractVector,eta_ddot::AbstractVector;
    gravity::AbstractVector=zeros(3))
    _tip_rotor_empty_coordinates(eta,eta_dot,eta_ddot)
    wb = R'*omega
    Jworld = R*rotor.matrix_J*R'
    spin = rotor_spin_data(rotor)
    Hworld = R*spin.H
    gyro_moment = cross(omega,Hworld)
    force = rotor.mass*(a-gravity)
    moment = Jworld*alpha+cross(omega,Jworld*omega)+gyro_moment
    base_wrench = vcat(force,moment)
    vb = R'*v
    body_velocity = vcat(vb,wb)
    mass_matrix = rotor_mass_data(rotor).M
    return (; base_wrench,eta_residual=Vector{eltype(base_wrench)}(),
        kinetic_energy=rotor.mass*dot(v,v)/2+dot(wb,rotor.matrix_J*wb)/2,
        mount_potential=zero(eltype(base_wrench)),
        damping_power=zero(eltype(base_wrench)),gyro_power=dot(omega,gyro_moment),
        body_velocity,body_momentum=mass_matrix*body_velocity,
        mass_matrix,gyro_matrix=spin.G,cm_body=zeros(eltype(base_wrench),3),
        spin_momentum=Hworld)
end
