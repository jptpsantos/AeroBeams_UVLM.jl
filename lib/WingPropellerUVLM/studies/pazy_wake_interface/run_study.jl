# Run from any directory with Julia 1.13. The source and project copies are temporary.
const package_root = normpath(joinpath(@__DIR__,"..",".."))
const repo_root = normpath(joinpath(package_root,"..",".."))
const original_project = joinpath(package_root,"examples")

function write_experiment_environment(scratch)
    pkg = joinpath(scratch,"WingPropellerUVLM")
    env = joinpath(scratch,"environment")
    mkpath(pkg)
    mkpath(env)
    cp(joinpath(package_root,"src"),joinpath(pkg,"src"))
    cp(joinpath(package_root,"Project.toml"),joinpath(pkg,"Project.toml"))
    backend = joinpath(pkg,"src","backend")
    cp(joinpath(@__DIR__,"wake_interface_correction.jl"),joinpath(backend,"wake_interface_correction.jl"))
    analyses = joinpath(backend,"analyses.jl")
    source = read(analyses,String)
    # Both the complete and cached matrix paths run before this hook. At this
    # point the original trailing-edge columns have been rebuilt in either path.
    anchor = "    # calculate RHS\n    if derivatives\n        # Derivative RHS follows"
    source = replace(source,"\r\n" => "\n")
    length(findall(anchor,source)) == 1 || error("Cannot identify the UVLM wake RHS")
    correction = "    match_wake_interface_core!(AIC,currentSurfaces,wakes,nwake,surface_id,wake_finite_core,\n" *
                 "        symmetric,trailing_vortices,interaction_id,interaction)\n\n"
    write(analyses,"include(\"wake_interface_correction.jl\")\n\n" * replace(source,anchor => correction * anchor))

    toml_path(path) = "\"" * replace(path,'\\' => '/') * "\""
    for filename in ("Project.toml","Manifest.toml")
        contents = read(joinpath(original_project,filename),String)
        patched = replace(contents,"path = \"../../..\"" => "path = " * toml_path(repo_root),
            "path = \"..\"" => "path = " * toml_path(pkg))
        write(joinpath(env,filename),patched)
    end
    return env
end

function run_study(args)
    parent = mkpath(joinpath(repo_root,"tmp"))
    mktempdir(parent) do scratch
        env = write_experiment_environment(scratch)
        runner = joinpath(@__DIR__,"analyze.jl")
        run(`$(Base.julia_cmd()) --startup-file=no --compiled-modules=existing --project=$env --threads=2 $runner $scratch $args`)
    end
end

run_study(ARGS)
