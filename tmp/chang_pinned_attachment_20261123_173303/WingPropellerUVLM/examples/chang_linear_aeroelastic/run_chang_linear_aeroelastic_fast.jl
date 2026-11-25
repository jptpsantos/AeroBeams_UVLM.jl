# Optional optimized runner; all existing source and case files stay intact.
# From the repository root:
# julia --threads=4 --project=lib/WingPropellerUVLM lib/WingPropellerUVLM/examples/chang_linear_aeroelastic/run_chang_linear_aeroelastic_fast.jl

# Select the package environment without Pkg.activate's global usage-log write.
let project = normpath(joinpath(@__DIR__,"..",".."))
    if Base.active_project() != joinpath(project,"Project.toml")
        pushfirst!(LOAD_PATH,project)
    end
end

if !isdefined(@__MODULE__,:ChangAeroelastic)
    include(joinpath(@__DIR__,"src","ChangAeroelastic.jl"))
end
if !isdefined(@__MODULE__,:ChangRuntimeOptimizations)
    include(joinpath(@__DIR__,"src","ChangRuntimeOptimizations.jl"))
end

function run_chang_fast_example()
    config = ChangAeroelastic.load_chang_configuration()
    saveHistory = ChangAeroelastic.environment_flag("CHANG_SAVE_SURFACE_HISTORY",false)
    threaded = ChangAeroelastic.environment_flag("CHANG_THREADED_WAKE",true)
    return ChangRuntimeOptimizations.run_chang_fast(config;
        save_surface_history=saveHistory,threaded_wake=threaded)
end

#if abspath(PROGRAM_FILE) == @__FILE__
    chang_run = run_chang_fast_example()
    results = chang_run.results
#end
