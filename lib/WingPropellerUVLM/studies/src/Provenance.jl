# Included inside ChangAerodynamicStudy (SHA and TOML are already imported).
# Deliberately conservative: each numerical identity covers indirect study and
# Chang dependencies. Existing schema-1 seals are never rewritten.
function numerical_source_fingerprint()
    package=normpath(joinpath(@__DIR__,"..",".."))
    paths=String[]
    for directory in (joinpath(package,"src"),@__DIR__,
        joinpath(package,"examples","chang_linear_aeroelastic","src"))
        for (root,_,files) in walkdir(directory), file in files
            endswith(file,".jl") && push!(paths,joinpath(root,file))
        end
    end
    for path in (joinpath(package,"Project.toml"),joinpath(package,"Manifest.toml"),
        joinpath(package,"studies","Project.toml"),joinpath(package,"studies","Manifest.toml"))
        isfile(path) && push!(paths,path)
    end
    io=IOBuffer(); print(io,"uvlm-migration-v2",'\0',VERSION,'\0')
    for path in sort(paths)
        print(io,relpath(path,package),'\0'); write(io,read(path)); write(io,UInt8(0))
    end
    return bytes2hex(sha256(take!(io)))
end

