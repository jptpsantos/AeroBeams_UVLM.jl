using AeroBeams

# Beam
L = 0.5
b,h = 5e-2, 1e-3
A = b*h
Iy, Iz = b*h^3/12, b^3*h/12
E = 70e9
ρ = 2.7e3
μ = 1e-4
∞ = 1e12
nElements = 20
beam = create_Beam(name="beam",length=L,nElements=nElements,S=[isotropic_stiffness_matrix(EA=E*A,EIy=E*Iy,EIz=E*Iz,∞=∞)],I=[inertia_matrix(ρA=ρ*A,ρIy=ρ*Iy,ρIz=ρ*Iz)],μ=μ)

# Cable stiffness range, as a multiple of the beam's tip bending stiffness
kBeam = 3*E*Iy/L^3
γRange = [0, 1, 4]
kCable = kBeam*γRange

# Cable force function: the cable lies above the tip, so it is stretched, and exerts an upwards force on the tip. Since the spring's force function receives the stiffness matrix (K, which has a single non-zero entry, along x3, as the spring's local x1 axis is rotated to x3) and the tip displacement vector (x), both resolved in basis A, the tip displacement is clipped to its negative part to make the cable tension-only
cableForce = (K,x) -> K*[x[1],x[2],min(x[3],0)]

# Boundary conditions (sinusoidal tip force with frequency as a multiple of the first natural frequency and root clamp)
F = 1e0
ω1 = (1.875^2) * sqrt(E*Iy/(ρ*A*L^4))
ω = 0.5 * ω1
force = create_BC(name="force",beam=beam,node=nElements+1,conditions=Dict(:F3A=>t->F*sin(ω*t)))
clamp = create_BC(name="clamp",beam=beam,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))

# Static tip displacement of the cantilever under the force amplitude alone
u3Static = F/kBeam

# Time variables
T = 2π/ω
cycles = 40
tf = cycles * T
nSteps = 200
Δt = T/nSteps

# Integration method
integrationMethod = GeneralizedAlpha(ρ∞=0.9)

# Initialize outputs
problem = Array{DynamicProblem}(undef, length(kCable))
t = Array{Vector{Float64}}(undef, length(kCable))
u3_tipNorm = Array{Vector{Float64}}(undef, length(kCable))
F3_rootNorm = Array{Vector{Float64}}(undef, length(kCable))
M2_rootNorm = Array{Vector{Float64}}(undef, length(kCable))

# Loop over cable stiffness range
for (i,k) in enumerate(kCable)
    println("Solving for cable/beam stiffness ratio γ = $(round(γRange[i], digits=2))")
    # Remove any existing cable before adding a new one
    remove_all_springs_from_beams!(beams=[beam])
    # Add the cable to the tip of the beam, in the perpendicular direction
    add_spring_to_beam!(beam=beam,spring=create_Spring(ku=k,forceFunction=cableForce),nodeID=nElements+1,p=[0,-π/2,0])
    # Model
    model = create_Model(name="tipSineLoadedCantileverWithCable",beams=[beam],BCs=[clamp,force])
    # Create and solve the problem
    problem[i] = create_DynamicProblem(model=model,finalTime=tf,Δt=Δt,integrationMethod=integrationMethod)
    solve!(problem[i])
    # Outputs
    sol = problem[i].solutionHistory
    t[i] = sol.t
    states = sol.states
    u3_tipNorm[i] = [states[j][end].nodal.u2[3] for j in eachindex(t[i])] / u3Static
    F3_rootNorm[i] = [states[j][1].nodal.F1[3] for j in eachindex(t[i])] / F
    M2_rootNorm[i] = [states[j][1].nodal.M1[2] for j in eachindex(t[i])] / (F*L)
end

## Plots
# Set paths
relPath = "/test/outputs/figures/dynamicStructural/tipSineLoadedCantileverWithCable"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
anim = plot_dynamic_deformations(problem,plotFrequency=div(nSteps,20),fps=24,plotLimits=([0,L],[-L/2,L/2],[-L/2,L/2]),normalizeTimeBy=T,normalizeByDim=(L,"\$L\$"),maxPlotTime=10T,legendEntries=["\$\\gamma = $(γRange[i])\$" for i in eachindex(γRange)],legendPos=:topleft,save=true,savePath=string(relPath,"/tipSineLoadedCantileverWithCable_deformation.gif"),displayProgress=true)
display(anim)

# Plot configurations
colors = cgrad(:rainbow, length(kCable), categorical=true)
labels = ["\$\\gamma = $(γRange[i])\$" for i in eachindex(γRange)]
ts = 10
fs = 16
lfs = 12
lw = 2
gr()

# Normalized tip displacement
plt1 = plot(xlabel="\$t/T\$", ylabel="Tip \$u_3\$ / static deflection", legendfont=font(lfs), tickfont=font(ts), guidefont=font(fs))
for (i,k) in enumerate(kCable)
    plot!(t[i]/T,u3_tipNorm[i], c=colors[i], lw=lw, label=labels[i])
end
display(plt1)
savefig(string(absPath,"/tipSineLoadedCantileverWithCable_u3.pdf"))

# Normalized root F3
plt2 = plot(xlabel="\$t/T\$", ylabel="Normalized root \$F_3^*\$", legendfont=font(lfs), tickfont=font(ts), guidefont=font(fs))
for (i,k) in enumerate(kCable)
    plot!(t[i]/T,F3_rootNorm[i], c=colors[i], lw=lw, label=labels[i])
end
display(plt2)
savefig(string(absPath,"/tipSineLoadedCantileverWithCable_F3.pdf"))

# Normalized root M2
plt3 = plot(xlabel="\$t/T\$", ylabel="Normalized root \$M_2^*\$", legendfont=font(lfs), tickfont=font(ts), guidefont=font(fs))
for (i,k) in enumerate(kCable)
    plot!(t[i]/T,M2_rootNorm[i], c=colors[i], lw=lw, label=labels[i])
end
display(plt3)
savefig(string(absPath,"/tipSineLoadedCantileverWithCable_M2.pdf"))

# Tip displacement PSD (Welch estimate: the spectrogram's power is averaged over its time segments, and the DC component is discarded for the logarithmic frequency axis)
plt4 = plot(xlabel="\$\\omega/\\omega_1\$", ylabel="Normalized tip \$u_3\$ PSD [dB]", xscale=:log10, xlims=[1e-2,5], ylims=[-75,Inf], xticks=([1e-2,5e-2,1e-1,2e-1,5e-1,1,2,5],["\$0.01\$","\$0.05\$","\$0.1\$","\$0.2\$","\$0.5\$","\$1\$","\$2\$","\$5\$"]), legendfont=font(lfs), tickfont=font(ts), guidefont=font(fs), legend=:topleft)
for (i, k) in enumerate(kCable)
    u3tipdev = u3_tipNorm[i] .- average(u3_tipNorm[i])
    f, _, psd = get_FFT_and_PSD(t[i], u3tipdev)
    plot!(2π*f/ω1, 10*log10.(psd), lw=lw, c=colors[i], label=labels[i])
end
display(plt4)
savefig(string(absPath,"/tipSineLoadedCantileverWithCable_psd.pdf"))

println("Finished tipSineLoadedCantileverWithCable.jl")