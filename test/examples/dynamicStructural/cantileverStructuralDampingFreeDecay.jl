using AeroBeams

# Free decay of a cantilever with stiffness-proportional (Kelvin-Voigt) structural damping, released from the first bending mode shape
L = 1.0
EI = 1.0
ρA = 1.0
μs = 1e-2
nElem = 50
δ = 1e-3

# First bending mode of the cantilever (Euler-Bernoulli)
βL = 1.87510407
β = βL/L
σβ = (cosh(βL) + cos(βL)) / (sinh(βL) + sin(βL))
φ(x) = 1/2*(cosh(β*x) - cos(β*x) - σβ*(sinh(β*x) - sin(β*x)))
dφ(x) = 1/2*β*(sinh(β*x) + sin(β*x) - σβ*(cosh(β*x) - cos(β*x)))

# Analytical frequency, undamped period, damping ratio and tip displacement of the first mode
ω1 = βL^2 * sqrt(EI/ρA) / L^2
T1 = 2π/ω1
ζ1 = μs * ω1 / 2
ωd = ω1*sqrt(1-ζ1^2)
u3Tip_analytic(t) = δ * φ(L) * exp.(-ζ1 * ω1 * t) .* (cos.(ωd*t) .+ ζ1/sqrt(1-ζ1^2) * sin.(ωd*t))

# Initial displacement and rotation parameters (small rotation) along the mode
u0(x1) = [0; 0; δ*φ(x1)]
p0(x1) = [0; 4*tan(-δ*dφ(x1)/4); 0]

# Beam
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(∞=1e6,EIy=EI)],I=[inertia_matrix(ρA=ρA)],μ=μs)
add_initial_conditions_to_beam!(beam=beam, conditions=Dict(:u0_of_x1 => u0, :p0_of_x1 => p0))

# BCs: clamped-free
clamp = create_BC(name="clamp",beam=beam,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))

# Model
cantileverStructuralDampingFreeDecay = create_Model(name="cantileverStructuralDampingFreeDecay",beams=[beam],BCs=[clamp])

# Integration method
integrationMethod = :newmark

# Time variables
N = 10
tf = N*T1
Δt = T1/100

# Create and solve the dynamic problem
problem = create_DynamicProblem(model=cantileverStructuralDampingFreeDecay,finalTime=tf,Δt=Δt,integrationMethod=integrationMethod)
solve!(problem)

# Outputs
sol = problem.solutionHistory
t = sol.t
tNorm = t/T1
u3Tip = [sol.states[i][nElem].nodal.u2[3] for i in eachindex(t)]

# Peak amplitude of the tip deflection within a time window
peak_amplitude(t, u, t0) = maximum(abs, u[(t .>= t0) .& (t .< t0 + T1)])

# Ratios of the amplitudes after N-1 periods, and their analytical values
ratio = peak_amplitude(t, u3Tip, (N-1)*T1) / peak_amplitude(t, u3Tip, 0.0)
ratioExpected = exp(-ζ1 * ω1 * (N-1)*T1)

## Plots
# Set paths
relPath = "/test/outputs/figures/dynamicStructural/cantileverStructuralDampingFreeDecay"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
anim = plot_dynamic_deformation(problem,scale=1/δ,plotLimits=([0,L],[-L/2,L/2],[-1,1]),plotFrequency=10,displayProgress=true,save=true,savePath=string(relPath,"/cantileverStructuralDampingFreeDecay_deformation.gif"))
display(anim)

# Plot configurations
lw = 2
ms = 3
msw = 0
gr()

# Normalized tip displacement
step = 5
plt1 = plot(xlabel="Undamped period \$t/T\$", ylabel="Tip \$u_3/\\delta\$")
plot!(tNorm,u3Tip/δ, c=:black, lw=lw, label="Numerical")
scatter!(tNorm[1:step:end],u3Tip_analytic.(t[1:step:end])/δ, c=:blue, ms=ms, msw=msw, label="Analytical")
display(plt1)
savefig(string(absPath,"/cantileverStructuralDampingFreeDecay_tipDisp.pdf"))

println("Finished cantileverStructuralDampingFreeDecay.jl")