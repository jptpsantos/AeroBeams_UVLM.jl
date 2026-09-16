# Chang wing–propeller aeroelastic response.
# Edit chang_case.jl for the case, numerical controls, and output.
# CHANG_* environment variables override those defaults when config is loaded.

import Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..", "..")))

if !isdefined(@__MODULE__, :ChangAeroelastic) ||
        !isdefined(ChangAeroelastic, :load_chang_configuration)
    include(joinpath(@__DIR__, "src", "ChangAeroelastic.jl"))
end
using .ChangAeroelastic

function run_chang_example()
    config = load_chang_configuration()
    return run_chang(config)
end

# Inspect chang_run.model, chang_run.workspace, or chang_run.solution as needed.
#if abspath(PROGRAM_FILE) == @__FILE__
    chang_run = run_chang_example()
    results = chang_run.results
#end
