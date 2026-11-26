using AeroBeams, LinearAlgebra

# Beam/assembly properties
L = 1
h = 0.5
EIy, EIz, GJ = 2, 3, 1
∞ = 1e12
nElem = 20

# Location of the springs as a fraction of the segment joining the tips of the cantilevers
α = 0.3

# Tip load applied to the first cantilever (in the reference assembly)
P = [0.2,-0.3,0.6]

# Rotations of the assembly about its clamp, given by the Euler angles (3-2-1 sequence): reference (none) and rotated by 45 deg about axes 3 and 2
anglesRange = [[0,0,0], [π/4,π/4,0]]

# Spring properties: translational (ku) or rotational (kp) stiffness, acting along (or about) the first axis of the basis given by the rotation parameters p. The first spring has a nonlinear force function, with a cubic term
ku = [5,8,3,0,0,0]
kp = [0,0,0,2,0.5,1.5]
pSprings = [[0,0,0],[0.3,-1.2,0.4],[π/2,0,0],[0.2,0.1,-0.3],[π/2,0,0],[0,-π/2,0]]
nonlinearForce = (K,x) -> K*x+10*K*x.^3
nSprings = length(ku)

# Initialize outputs
problem = Array{SteadyProblem}(undef,length(anglesRange))
springs = Array{Vector{Spring}}(undef,length(anglesRange))
R = Array{Matrix{Float64}}(undef,length(anglesRange))
tipDisp = Array{Vector{Float64}}(undef,length(anglesRange))
clampLoads = Array{Vector{Float64}}(undef,length(anglesRange))

# Loop over the rotations of the assembly and solve the corresponding problems
for (i,angles) in enumerate(anglesRange)
    # Rotation tensor of the assembly, and the rotated tip load
    R[i] = rotation_tensor_E321(angles)
    RP = R[i]*P
    # Cantilevers, connected through a rigid bar at their roots (along x2)
    S = [isotropic_stiffness_matrix(∞=∞,EIy=EIy,EIz=EIz,GJ=GJ)]
    cantilever1 = create_Beam(name="cantilever1",length=L,nElements=nElem,S=S)
    rootBar = create_Beam(name="rootBar",length=h,nElements=1,S=[isotropic_stiffness_matrix(∞=∞)],p0=transform_to_WM([π/2,0,0],from=:E321))
    cantilever2 = create_Beam(name="cantilever2",length=L,nElements=nElem,S=S)
    # Beam connections
    connections = [
        (rootBar => cantilever1, 1 => 1),
        (cantilever2 => rootBar, 1 => 2),
    ]
    # Springs between the tips of the cantilevers (created anew, as they store their states)
    springs[i] = [create_Spring(ku=ku[j],kp=kp[j],forceFunction=(j==1 ? nonlinearForce : (K,x)->K*x)) for j in 1:nSprings]
    for j in 1:nSprings
        add_spring_to_beams!(beams=[cantilever1,cantilever2],spring=springs[i][j],nodesLocalIDs=[nElem+1,nElem+1],p=pSprings[j],locationFraction=α)
    end
    # BCs: clamp at the root of the first cantilever, rotated by R, and tip load
    pR = rotation_parameters_WM(R[i])
    clamp = create_BC(name="clamp",beam=cantilever1,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>pR[1],:p2A=>pR[2],:p3A=>pR[3]))
    tipLoad = create_BC(name="tipLoad",beam=cantilever1,node=nElem+1,conditions=Dict(:F1A=>RP[1],:F2A=>RP[2],:F3A=>RP[3]))
    # Model
    model = create_Model(name="nonlinearDoublyAttachedSprings",beams=[cantilever1,rootBar,cantilever2],BCs=[clamp,tipLoad],connections=connections)
    # Create and solve problem
    problem[i] = create_SteadyProblem(model=model,systemSolver=create_NewtonRaphson(absTol=1e-12,relTol=1e-12))
    solve!(problem[i])
    # Get outputs
    tipDisp[i] = problem[i].loadstepHistory.states[end][nElem].nodal.u2
    specialNodes = problem[i].model.assemblyData.specialNodes
    clampNode = specialNodes[findfirst(n -> n.assemblyData.globalID == cantilever1.assemblyData.nodesGlobalIDs[1], specialNodes)]
    clampLoads[i] = vcat(clampNode.state.F,clampNode.state.M)
end

# Undeformed position of the tip of the first cantilever
rTip = [L,0,0]

# Errors of the rotated assembly w.r.t. the reference one: the springs' deformations must be invariant, and their loads must rotate with the assembly; the tip displacement must be the rigidly rotated one (up to the discretization error, as the interpolation of the rotation parameters is not exactly objective)
Rrot = R[2]
ϵSprings = maximum(maximum(abs, vcat(sRot.state.Δu-sRef.state.Δu, sRot.state.Δp-sRef.state.Δp, Rrot'*sRot.state.Fs-sRef.state.Fs, Rrot'*sRot.state.Ms-sRef.state.Ms)) for (sRef,sRot) in zip(springs[1],springs[2]))
ϵTip = maximum(abs, (Rrot-AeroBeams.I3)*rTip + Rrot*tipDisp[1] - tipDisp[2])
println("Max. error of the springs' deformations and loads = $ϵSprings")
println("Max. error of the tip displacement = $ϵTip")

# Equilibrium of the clamp's reactions (force, and moment about the origin) with the tip load, for each assembly
ϵEquilibrium = zeros(length(anglesRange))
for i in eachindex(anglesRange)
    RP = R[i]*P
    ΣF = clampLoads[i][1:3] + RP
    ΣM = clampLoads[i][4:6] + cross(rTip+tipDisp[i],RP)
    ϵEquilibrium[i] = max(norm(ΣF),norm(ΣM))
end
println("Max. error of the equilibrium of the clamp's reactions = $(maximum(ϵEquilibrium))")

## Plots
# Set paths
relPath = "/test/outputs/figures/staticStructural/nonlinearDoublyAttachedSprings"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Deformed shapes of the reference and rotated assemblies
deformationRef = plot_steady_deformation(problem[1],view=(30,30),plotUndeformed=true,save=true,savePath=string(relPath,"/nonlinearDoublyAttachedSprings_deformation_ref.pdf"))
display(deformationRef)
deformationRot = plot_steady_deformation(problem[2],view=(30,30),plotUndeformed=true,save=true,savePath=string(relPath,"/nonlinearDoublyAttachedSprings_deformation_rot.pdf"))
display(deformationRot)

println("Finished nonlinearDoublyAttachedSprings.jl")