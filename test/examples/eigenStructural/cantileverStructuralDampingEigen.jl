using AeroBeams

# Beam with stiffness-proportional (Kelvin-Voigt) structural damping
L = 1.0
EI = 1.0
ρA = 1.0
μs = 1e-3
nElem = 20
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(∞=1e6,EIy=EI)],I=[inertia_matrix(ρA=ρA)],μ=μs)

# BCs: clamped-free
clamp = create_BC(name="clamp",beam=beam,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))

# Model
cantileverStructuralDampingEigen = create_Model(name="cantileverStructuralDampingEigen",beams=[beam],BCs=[clamp])

# Eigenproblem
problem = create_EigenProblem(model=cantileverStructuralDampingEigen,systemSolver=create_LinearSolver(),eigenmethod=:krylovkit,nev=12)
solve!(problem)

# Oscillatory eigenvalues of the first bending modes
λ = problem.eigensolution.λ
λ = λ[imag.(λ) .> 1e-6]
λ = λ[sortperm(imag.(λ))][1:4]

# Analytical undamped frequencies of the Euler-Bernoulli cantilever
βL = [1.87510407, 4.69409113, 7.85475744, 10.99554073]
ωAnalytical = (βL/L).^2 * sqrt(EI/ρA)

# Undamped frequencies (modulus of the eigenvalues) and damping ratios of the discrete model
ωNum = abs.(λ)
ζNum = -real.(λ) ./ ωNum

# Rayleigh damping: ζ = μs*ω/2 for each mode, with ω the undamped frequency
ζExpected = μs * ωNum / 2

# Relative error between numerical and expected damping ratios
ζRelError = abs.(ζNum .- ζExpected) ./ ζExpected

# Compare numerical and expected damping ratios
println("Relative error between numerical and expected damping ratios:")
for i in 1:length(ζNum)
    println("Mode $i: $(round(ζRelError[i] * 100, sigdigits=2)) %")
end

println("Finished cantileverStructuralDampingEigen.jl")