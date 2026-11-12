# Synthetic structural history for rendering QA, not an aeroelastic simulation.
using Test
const ROOT = normpath(joinpath(@__DIR__,"..","..",".."))
const EXAMPLE = joinpath(ROOT,"lib","WingPropellerUVLM","examples","pazy_aerobeams_coupled")
include(joinpath(EXAMPLE,"PazyWingUVLMCoupling.jl"))
include(joinpath(EXAMPLE,"PazyWingUVLMVisualization.jl"))
ne,span,chord,spar,nodes = AeroBeams.geometrical_properties_Pazy()
beam = AeroBeams.create_Beam(length=span,nElements=ne,normalizedNodalPositions=nodes,
    p0=AeroBeams.transform_to_WM([-pi/2,0.0,0.0]; from=:E231),
    S=AeroBeams.stiffness_matrices_Pazy(withSkin=true,sweepStructuralCorrections=false,GAy=1e16,GAz=1e16,Λ=0.0),
    I=AeroBeams.inertia_matrices_Pazy(withSkin=true))
clamp = AeroBeams.create_BC(beam=beam,node=1,
    conditions=Dict(type=>0.0 for type in (:u1A,:u2A,:u3A,:p1A,:p2A,:p3A)))
model = AeroBeams.create_Model(beams=[beam],BCs=[clamp])
times = collect(0:4) .* 0.001
problem = AeroBeams.create_DynamicProblem(model=model,timeVector=times,displayProgress=false)
prepare_structural_problem!(problem)
for (index,t) in enumerate(times)
    for (e,element) in enumerate(model.assemblyData.elements)
        element.state.nodal.u1[1] = -0.02*(index-1)/4 * nodes[e]^2
        element.state.nodal.u2[1] = -0.02*(index-1)/4 * nodes[e+1]^2
    end
    AeroBeams.save_time_step_data!(problem)
end
saved = deepcopy(problem.solutionHistory.states)
result = (; structural_problem=problem,time=times,dt=0.001,animation_stride=3,
    aerodynamic_nodal_load_history=zeros(6,ne+1,length(times)))
folder = joinpath(@__DIR__,"render_qa")
mkpath(folder)
path = save_structural_animation(result; output_path=joinpath(folder,"pazy_synthetic.gif"),fps=2)
Plots.savefig(joinpath(folder,"pazy_synthetic_final.png"))
@test isfile(path)
@test all(problem.solutionHistory.states[k][e].nodal.u2 == saved[k][e].nodal.u2
    for k in eachindex(saved) for e in eachindex(saved[k]))
@test beam.attachments.surface === nothing
@test model.BCs[1] === clamp
println("Synthetic animation and final frame verified; source model unchanged")
