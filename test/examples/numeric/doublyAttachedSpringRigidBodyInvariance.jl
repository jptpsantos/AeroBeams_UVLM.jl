# Invariance and equilibrium of a doubly-attached spring
# A doubly-attached spring connects two nodes, at positions `r1` and `r2`, through rigid arms to a point located at a fraction `α` of the segment joining them. The spring's deformations and nodal loads are defined with respect to the mean rotation of its two ends. Therefore, for arbitrary (large) displacements and rotations of the nodes:
#
# - the spring's deformations must be invariant to rigid-body motions of the system (a rotation `R` and a translation `c`, applied to both nodes), whereas its nodal loads must rotate with `R`;
# - the nodal loads must be in equilibrium, i.e., their resultant force, and their moment about the origin, must vanish.
#
# No beams are involved in this example. The spring's kinematics and nodal loads are evaluated directly, for random states of the nodes and several locations of the spring.

# Given the positions of the nodes in the undeformed configuration (`r1` and `r2`), their displacements (`u1` and `u2`) and rotation parameters (`p1` and `p2`), and the spring's force (`Fs`) and moment (`Ms`), the function below returns the spring's deformations (`Δu` and `Δp` of the first node relative to the second) and the loads it exerts on each node (forces `F1` and `F2`, and moments `M1` and `M2`), all resolved in basis A. It relies on the internal functions of the package that compute the spring's kinematics and nodal loads.

using AeroBeams, LinearAlgebra, Random

# Kinematics and nodal loads of a doubly-attached spring
function spring_response(; r1,u1,p1,r2,u2,p2,α,Fs,Ms)
    # Arms from the nodes to the spring's location, resolved in basis A in the undeformed configuration
    d1,d2 = α*(r2-r1),(α-1)*(r2-r1)
    # Compute the spring's kinematics and nodal loads
    Δu,Δp,Rm,a1,a2 = AeroBeams.double_attachment_spring_kinematics(u1,p1,u2,p2,d1,d2)
    F1,M1,F2,M2 = AeroBeams.double_attachment_spring_nodal_loads(Fs,Ms,Rm,a1,a2)
    return (Δu=Δu,Δp=Δp,F1=F1,M1=M1,F2=F2,M2=M2)
end

# Displacement and rotation parameters of a node after a rigid-body motion (rotation R and translation c) of the system, given its displacement u, rotation parameters p and undeformed position r
function rigid_motion(; u,p,r,R,c)
    return c+R*(r+u)-r, rotation_parameters_WM(R*rotation_tensor_WM(p))
end

# Range of spring positions along the segment connecting the nodes
αRange = [0,0.3,1]

# For each location of the spring and each case, the nodes are placed at random positions, with random (large) displacements and rotations, and the spring has a random force and moment. The nodes then undergo a random rigid-body motion, and the residuals are: the difference between the spring's deformations after and before the motion, and between its nodal loads after the motion and those before it rotated by `R` (both must vanish), and the equilibrium residual of its nodal loads before the motion (the resultant force and the moment about the origin must vanish).
ϵKinematics,ϵLoads,ϵEquilibrium = (zeros(length(αRange)) for _ in 1:3)
for (i,α) in enumerate(αRange)
    # Random undeformed positions of the nodes, displacements and rotations of the nodes, and force and moment of the spring
    r1,r2 = randn(3),randn(3)
    u1,p1,u2,p2 = randn(3),2randn(3),randn(3),2randn(3)
    Fs,Ms = randn(3),randn(3)
    # Random rigid-body motion of the system
    R,c = rotation_tensor_E321(randn(3)),randn(3)
    # Response of the spring before and after the rigid-body motion
    resp = spring_response(r1=r1,u1=u1,p1=p1,r2=r2,u2=u2,p2=p2,α=α,Fs=Fs,Ms=Ms)
    u1R,p1R = rigid_motion(u=u1,p=p1,r=r1,R=R,c=c)
    u2R,p2R = rigid_motion(u=u2,p=p2,r=r2,R=R,c=c)
    respR = spring_response(r1=r1,u1=u1R,p1=p1R,r2=r2,u2=u2R,p2=p2R,α=α,Fs=Fs,Ms=Ms)
    ϵKinematics[i] = maximum(abs, vcat(respR.Δu-resp.Δu, respR.Δp-resp.Δp))
    ϵLoads[i] = maximum(abs, vcat(respR.F1-R*resp.F1, respR.M1-R*resp.M1, respR.F2-R*resp.F2, respR.M2-R*resp.M2))
    ϵEquilibrium[i] = maximum(abs, vcat(resp.F1+resp.F2, cross(r1+u1,resp.F1)+resp.M1+cross(r2+u2,resp.F2)+resp.M2))
end

# The maximum residuals over the random cases are shown for each location of the spring. All of them must be at the level of the round-off error.
table = pretty_table(String,hcat(αRange,maximum(ϵKinematics,dims=2),maximum(ϵLoads,dims=2),maximum(ϵEquilibrium,dims=2)),column_labels=["α","Kinematics","Nodal loads","Equilibrium"],formatters=[(v,i,j) -> round(v,sigdigits=3)])
print(table)

println("Finished doublyAttachedSpringRigidBodyInvariance.jl")