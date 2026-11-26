using AeroBeams

# Beam/assembly properties
L = 1
h = 0.5
EI, GJ = 1, 1
∞ = 1e12
nElem = 50

# First cantilever
cantilever1 = create_Beam(name="cantilever1",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(∞=∞,EIy=EI,EIz=EI,GJ=GJ)])

# Rigid bar connecting the roots of the cantilevers (along x2)
rootBar = create_Beam(name="rootBar",length=h,nElements=1,S=[isotropic_stiffness_matrix(∞=∞)],p0=transform_to_WM([π/2,0,0],from=:E321))

# Second cantilever
cantilever2 = create_Beam(name="cantilever2",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(∞=∞,EIy=EI,EIz=EI,GJ=GJ)])

# Beam connections
connections = [
    (rootBar => cantilever1, 1 => 1),
    (cantilever2 => rootBar, 1 => 2),
]

# Spring stiffness, applied load (small, to keep the response linear) and range of spring positions along the connecting segment
k = 1e0
P = 1e-3
αRange = vcat(0:0.25:1)

# BCs
clamp1 = create_BC(name="clamp1",beam=rootBar,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))
clamp2 = create_BC(name="clamp2",beam=rootBar,node=2,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))
tipLoad = create_BC(name="tipLoad",beam=cantilever1,node=nElem+1,conditions=Dict(:F3A=>P))

# Normalized compliances and analytical solution (The cantilevers are compliant in bending, with `cb=L³/(3EIy)`, corrected by the discretization factor '(1-1/(4nElem^2))', and in torsion, with `ct=L/GJ`. For a spring force `F=kΔ`, where `Δ` is the transverse relative displacement of the spring's ends, the tips are displaced by `w1=(P-F)cb` and `w2=F*cb`, and twisted by `θ1=-αFh*ct` and `θ2=-(1-α)Fh*ct`, due to the moments `αhF` and `(1-α)hF` exerted by the spring on them. The relative displacement of the spring's ends is `Δ=w1-w2+αhθ1+(1-α)hθ2`, which yields a closed-form expression for `Δ`, and hence for all the other quantities. To be dimensionless, the displacements are normalized by `P*cb` and the twists by `P*h*ct`, so that the solution is of order unity whatever the load.)
cb = L^3/(3EI)*(1-1/(4nElem^2))
ct = L/GJ
function analytical_solution(α)
    Δ = 1/(1+2k*cb+k*h^2*ct*(α^2+(1-α)^2))
    F = k*cb*Δ
    return Δ,1-F,F,-α*F,-(1-α)*F
end

# Initialize outuputs
problem = Array{SteadyProblem}(undef,length(αRange))
springDisp = Array{Float64}(undef,length(αRange))
tipsDisp = Array{Vector{Float64}}(undef,length(αRange))
tipsTwist = Array{Vector{Float64}}(undef,length(αRange))

# Loop over the range of spring positions and solve the corresponding problems
for (i,α) in enumerate(αRange)
    # Remove all previously added springs
    remove_all_springs_from_beams!(beams=[cantilever1,cantilever2])
    # Create a new spring for the current position
    spring = create_Spring(ku=k)
    # Set spring between the tips of the cantilevers
    add_spring_to_beams!(beams=[cantilever1,cantilever2],spring=spring,nodesLocalIDs=[nElem+1,nElem+1],p=[0,-π/2,0],locationFraction=α)
    # Model
    model = create_Model(beams=[cantilever1,rootBar,cantilever2],BCs=[clamp1,clamp2,tipLoad],connections=connections)
    # Create and solve problem
    problem[i] = create_SteadyProblem(model=model)
    solve!(problem[i])
    # Get outputs
    tip1 = problem[i].loadstepHistory.states[end][nElem].nodal
    tip2 = problem[i].loadstepHistory.states[end][end].nodal
    springDisp[i] = spring.state.Δu[3]/(P*cb)
    tipsDisp[i] = [tip1.u2[3], tip2.u2[3]]/(P*cb)
    tipsTwist[i] = [tip1.p2[1], tip2.p2[1]]/(P*h*ct)
end

# Errors relative to the analytical solution
ϵΔ,ϵw1,ϵw2,ϵθ1,ϵθ2 = (zeros(length(αRange)) for _ in eachindex(αRange))
for (i,α) in enumerate(αRange)
    Δa,w1a,w2a,θ1a,θ2a = analytical_solution(α)
    ϵΔ[i] = springDisp[i]/Δa - 1
    ϵw1[i] = tipsDisp[i][1]/w1a - 1
    ϵw2[i] = tipsDisp[i][2]/w2a - 1
    ϵθ1[i] = (tipsTwist[i][1]-θ1a)/(k*cb*Δa)
    ϵθ2[i] = (tipsTwist[i][2]-θ2a)/(k*cb*Δa)
end
ϵAnalytical = maximum(abs, vcat(ϵΔ,ϵw1,ϵw2,ϵθ1,ϵθ2))
println("Max. error w.r.t. the analytical solution = $ϵAnalytical")

## Plots
# Set paths
relPath = "/test/outputs/figures/staticStructural/transverseDoublyAttachedSpringCantilevers"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Deformed shape for the last location of the spring (α=1)
deformationPlot = plot_steady_deformation(problem[end],plotUndeformed=true,save=true,savePath=string(relPath,"/transverseDoublyAttachedSpringCantilevers_deformation.pdf"))
display(deformationPlot)

# Plot configurations
lw = 2
gr()

# Analytical solutions as continuous functions of the location of the spring
αc = range(0,1,length=101)
sols = analytical_solution.(αc)
Δc,w1c,w2c,θ1c,θ2c = (getindex.(sols,i) for i in eachindex(sols))

# Transverse displacements of the spring and of the tips
plt_u3 = plot(αc,Δc,c=:black,lw=lw,label="Spring (analytical)",xlabel="\$\\alpha\$",ylabel="Normalized transverse displacement",xlims=[0,1])
plot!(αc,w1c,c=:blue,lw=lw,label="Tip 1 (analytical)")
plot!(αc,w2c,c=:red,lw=lw,label="Tip 2 (analytical)")
scatter!(αRange,springDisp,c=:black,label="Spring (AeroBeams)")
scatter!(αRange,[tipsDisp[i][1] for i in eachindex(αRange)],c=:blue,label="Tip 1 (AeroBeams)")
scatter!(αRange,[tipsDisp[i][2] for i in eachindex(αRange)],c=:red,label="Tip 2 (AeroBeams)")
display(plt_u3)
savefig(string(absPath,"/transverseDoublyAttachedSpringCantilevers_u3.pdf"))

# Tip twists
plt_θ = plot(αc,θ1c,c=:blue,lw=lw,label="Tip 1 (analytical)",xlabel="\$\\alpha\$",ylabel="Normalized tip twist",xlims=[0,1])
plot!(αc,θ2c,c=:red,lw=lw,label="Tip 2 (analytical)")
scatter!(αRange,[tipsTwist[i][1] for i in eachindex(αRange)],c=:blue,label="Tip 1 (AeroBeams)")
scatter!(αRange,[tipsTwist[i][2] for i in eachindex(αRange)],c=:red,label="Tip 2 (AeroBeams)")
display(plt_θ)
savefig(string(absPath,"/transverseDoublyAttachedSpringCantilevers_theta.pdf"))

println("Finished transverseDoublyAttachedSpringCantilevers.jl")