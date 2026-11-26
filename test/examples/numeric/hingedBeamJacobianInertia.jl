using AeroBeams, Random, LinearAlgebra

# Hinged beams (hinge at a node, with one or two hinged rotational DoFs), with and without structural damping, and with and without an aerodynamic surface along the beam. The analytical Jacobian and inertia matrices are checked against finite differences of the residual, with the states perturbed around the solution

# Beam with a hinge, clamped at one end and loaded at the other
function hinged_beam_problem(; hingedNodes, hingedNodesDoF, μ, method, ppSolver=nothing)
    S = isotropic_stiffness_matrix(∞=1e4,EIy=10.0,EIz=20.0,GJ=5.0)
    Ii = inertia_matrix(ρA=1.0,ρIy=1e-2,ρIz=1e-2,ρIyz=1e-3,e2=0.01,e3=-0.02)
    surfaceKwargs = isnothing(ppSolver) ? (;) : (; surface=create_AeroSurface(ppSolver=ppSolver,airfoil=create_Airfoil(name="FlatPlate"),c=0.2,referenceAxis=0.3))
    beam = create_Beam(name="beam",length=1.0,nElements=4,S=[S],I=[Ii],μ=μ,k=[0.3,0.2,0.4],p0=[0.1,-0.2,0.05],hingedNodes=hingedNodes,hingedNodesDoF=hingedNodesDoF;surfaceKwargs...)
    clamp = create_BC(name="clamp",beam=beam,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))
    tip = create_BC(name="tip",beam=beam,node=5,conditions=Dict(:F3A=>t->0.5*sin(20t),:M2A=>t->0.2*cos(15t)))
    model = create_Model(name="hingedBeam",beams=[beam],BCs=[clamp,tip],v_A=isnothing(ppSolver) ? [0,0,0] : [0,20,0])
    problem = create_DynamicProblem(model=model,finalTime=0.01,Δt=0.005,integrationMethod=method,displayProgress=false,systemSolver=create_NewtonRaphson(displayStatus=false))
    solve!(problem)
    return problem
end

# Perturbs the states (reproducibly) and checks the analytical Jacobian against finite differences
function jacobian_mismatches(problem)
    Random.seed!(3)
    problem.systemArrays.x .+= 5e-2 .* randn(length(problem.systemArrays.x)) .* max.(abs.(problem.systemArrays.x), 1e-2)
    AeroBeams.f!(problem)
    problem.systemSolver.state.skipJacobianUpdate = false
    J = copy(AeroBeams.j!(problem))
    labels = AeroBeams.build_jacobian_labels(problem)
    return AeroBeams.compare_jacobians(problem, J, labels; h0=1e-7, rtol=1e-5, atol=1e-6, tag="")
end

# Inertia matrix, which is the derivative of the Jacobian w.r.t. the coefficient c of the states in their rates (the rates themselves are kept fixed)
function inertia_error(problem)
    elements = problem.model.assemblyData.elements
    x0 = copy(problem.systemArrays.x)
    sa = AeroBeams.initialize_system_arrays!(problem.model,problem.systemSolver,zeros(0); isEigenProblem=true)
    sa.x .= x0
    setfield!(problem, :systemArrays, sa)
    AeroBeams.f!(problem)
    AeroBeams.i!(problem)
    M = copy(problem.systemArrays.inertia)
    c0 = elements[1].stateCache.equivStatesRates.c
    h = 1e-5 * c0
    Jat(c) = (for e in elements; e.stateCache.equivStatesRates.c = c; end; AeroBeams.j!(problem); copy(problem.systemArrays.jacobian))
    Mfd = (Jat(c0+h) - Jat(c0-h)) / (2h)
    for e in elements; e.stateCache.equivStatesRates.c = c0; end
    return maximum(abs, Matrix(M - Mfd)) / maximum(abs, Matrix(M))
end

hingeCases = [
    (hingedNodes=[3], hingedNodesDoF=[[false,true,false]]),
    (hingedNodes=[2], hingedNodesDoF=[[true,true,false]]),
]
jacobianMismatches = Int[]
inertiaErrors = Float64[]
for hinge in hingeCases, μ in (0.0, 1e-3), method in (:newmark, :bdf2, :gen_alpha)
    problem = hinged_beam_problem(; hinge..., μ=μ, method=method)
    push!(jacobianMismatches, length(jacobian_mismatches(problem)))
    push!(inertiaErrors, inertia_error(problem))
end

# With an aerodynamic surface along the beam (the elements to the right of the hinge carry aerodynamic loads in the hinged moment equations), with airfoil parameters updated with the Mach number
aeroHingeCases = [
    (hingedNodes=[3], hingedNodesDoF=[[false,true,false]]),
    (hingedNodes=[2], hingedNodesDoF=[[true,true,false]]),
    (hingedNodes=[3], hingedNodesDoF=[[true,false,true]]),
]
aeroJacobianMismatches = Int[]
for hinge in aeroHingeCases, ppSolver in (QuasiSteady(), Indicial()), μ in (0.0, 1e-3), method in (:newmark, :bdf2, :gen_alpha)
    problem = hinged_beam_problem(; hinge..., μ=μ, method=method, ppSolver=ppSolver)
    push!(aeroJacobianMismatches, length(jacobian_mismatches(problem)))
end

println("Finished hingedBeamJacobianInertia.jl")