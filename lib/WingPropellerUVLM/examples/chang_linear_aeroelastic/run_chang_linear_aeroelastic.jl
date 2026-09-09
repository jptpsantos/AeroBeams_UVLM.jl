# Chang wing–propeller aeroelastic response.
# Edit chang_case.jl for the case, numerical controls, and output.
# CHANG_* environment variables override those defaults when config is loaded.

import Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..", "..")))

if !isdefined(@__MODULE__, :ChangAeroelastic)
    include(joinpath(@__DIR__, "src", "ChangAeroelastic.jl"))
else
    # Pick up edits to the case when rerunning in the same Julia/IDE session.
    Base.include(ChangAeroelastic, joinpath(@__DIR__, "chang_case.jl"))
end
using .ChangAeroelastic

config = load_chang_configuration()
chang_run = run_chang(config)

# Inspect chang_run.model, chang_run.workspace, or chang_run.solution as needed.
results = chang_run.results
