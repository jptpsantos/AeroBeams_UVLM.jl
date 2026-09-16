"""Build the existing Chang benchmark from explicit camelCase or legacy settings.
No CHANG_* environment overrides are applied by this constructor.
"""
function create_ChangModel(settings::NamedTuple)
    Convergence.load_model()
    moduleRef=Convergence.Model.ChangAeroelastic
    config=Base.invokelatest(moduleRef.load_chang_configuration,legacy_case_settings(settings);env=Dict{String,String}())
    return Base.invokelatest(moduleRef.build_chang_model,config)
end
legacy_case_settings(settings::NamedTuple)=translate_settings(settings,
    key -> key==:azimuthStepDeg ? :azimuth_step_deg : legacy_key(key))
mutable struct ChangDynamicProblem{M}
    model::M
    results::Union{Nothing,NamedTuple}
end
create_ChangDynamicProblem(;model)=ChangDynamicProblem(model,nothing)
function solve!(p::ChangDynamicProblem)
    check_loaded_source()
    p.results===nothing || throw(ArgumentError("Create a fresh ChangDynamicProblem for another run"))
    m=Convergence.Model.ChangAeroelastic
    Base.invokelatest(m.load_chang_visualization,p.model.config.output)
    run=Base.invokelatest(m.run_chang_case,p.model.config;model=p.model)
    completed=run.solution.last_step==length(p.model.parameters.dt)
    solverConverged=completed && all(run.solution.coupling_converged)
    p.results=merge(run,(;completed,solverConverged))
    return p
end
