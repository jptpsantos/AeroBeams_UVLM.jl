using LinearAlgebra: BLAS
using Dates
import WingPropellerUVLM

BLAS.set_num_threads(1)
Threads.nthreads() == 10 || error("Start this study with --threads=10")
package = joinpath(@__DIR__,"WingPropellerUVLM")
normpath(dirname(dirname(pathof(WingPropellerUVLM)))) == normpath(package) ||
    error("Start Julia with the isolated studies project")
example = joinpath(package,"examples","chang_linear_aeroelastic")
include(joinpath(example,"studies","convergence","chang_convergence_aeroelastic.jl"))
const C = ChangConvergenceAeroelastic
const Study = C.Convergence

workspace = normpath(joinpath(@__DIR__,"..",".."))
root = joinpath(workspace,"lib","WingPropellerUVLM","examples","chang_linear_aeroelastic","output")
fullSweep = "--full" in ARGS
outputArgs = filter(!=("--full"),ARGS)
panelCounts = fullSweep ? [5,10,20,40,60,80] : [20,40,80]
prefix = fullSweep ? "convergence_pinned_full_span_10threads_" : "convergence_pinned_span_10threads_"
directory = isempty(outputArgs) ? mktempdir(root;prefix,cleanup=false) : abspath(only(outputArgs))
isdir(directory) || mkpath(directory)
isempty(readdir(directory)) || error("Use an empty output directory for this campaign")

settings = merge(C.aeroelastic_study(; dampingTarget=:wing),(;
    attachment_reference_panels=first(panelCounts),
    output_directory=directory,
    make_plots=false,
    families=[:wing_span],
    sweeps=(; wing_span=panelCounts)
))
selection = C.explicit_selection(settings)
@assert selection.physical.speed_mps == 65.0
@assert selection.selected.azimuth_deg == 10.0
@assert selection.selected.wake_revolutions == 1.0

Study.load_model()
const Chang = Study.Model.ChangAeroelastic

function verify_pinned_model(c,settings,selection)
    defaults = Study.model_defaults(selection.physical,selection.core_mode,c,settings,
        joinpath(settings.output_directory,c.label))
    config = Chang.load_chang_configuration(defaults;env=Dict{String,String}())
    model = Chang.build_chang_model(config)
    parameters = model.parameters
    nodes = parameters.span_nodes
    node = only(parameters.prop_attach_nodes)
    pivotY = selection.physical.attachment_eta * selection.physical.wing_span_m
    @assert length(nodes) == c.wing_span + 1
    @assert all(>(0),diff(nodes))
    @assert nodes[node] == pivotY
    workspace = Chang.build_chang_workspace(model)
    @assert isapprox(workspace.grid_wing_init[2,1,:],nodes;rtol=0,atol=1e-12)
    q = zeros(model.structural.ndof_wing_free + parameters.ndof_P)
    Chang.update_aero_geometry_for_state!(model,workspace,q,0.0)
    @assert isapprox(workspace.attach_node_y[1],pivotY;rtol=0,atol=1e-12)
    @assert isapprox(workspace.T_pivot_A_current[1][2],pivotY;rtol=0,atol=1e-12)
    @assert isapprox(workspace.T_hub_A_current[1][2],pivotY;rtol=0,atol=1e-12)
    return nodes
end

cases = Study.build_cases(settings.nominal,settings.sweeps,settings.families)
nodes = [verify_pinned_model(c,settings,selection) for c in cases]
@assert [c.wing_span for c in cases] == panelCounts
for i in eachindex(panelCounts), j in i+1:length(panelCounts)
    refinement,remainder = divrem(panelCounts[j],panelCounts[i])
    if iszero(remainder)
        @assert isapprox(nodes[i],nodes[j][1:refinement:end];rtol=0,atol=1e-12)
    end
end
println("Verified beam and UVLM attachment at y = 6.225 m on panel counts $(panelCounts)")
println("RESULTS_DIRECTORY=",directory)
flush(stdout)
result = C.run_explicit_aeroelastic(settings)
println("FINISHED=",result.directory)
