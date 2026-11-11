# Bare Chang wing: edit chang_isolated_wing_case.jl, then run this script.
if Base.active_project() != normpath(joinpath(@__DIR__,"..","..","Project.toml"))
    import Pkg
    Pkg.activate(normpath(joinpath(@__DIR__,"..","..")))
end

if !isdefined(@__MODULE__,:ChangAeroelastic) ||
   !isdefined(ChangAeroelastic,:load_chang_isolated_wing_configuration)
    include(joinpath(@__DIR__,"src","ChangAeroelastic.jl"))
end
using .ChangAeroelastic

function run_chang_isolated_wing_example()
    return run_chang(load_chang_isolated_wing_configuration())
end

# Also runs when included from the IDE; the case file is reread on every run.
chang_isolated_wing_run = run_chang_isolated_wing_example()
results = chang_isolated_wing_run.results
