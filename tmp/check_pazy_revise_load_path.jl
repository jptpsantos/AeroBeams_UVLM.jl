# Isolated reproduction probe. No source files or package environments change.
using Revise
repo = normpath(joinpath(@__DIR__,".."))
uvlm = joinpath(repo,"lib","WingPropellerUVLM")
for package_root in (repo,uvlm)
    package_root in LOAD_PATH || pushfirst!(LOAD_PATH,package_root)
end
import PrettyTables, Interpolations
println("Active project: ",Base.active_project())
println("LOAD_PATH: ",LOAD_PATH)
for (pkg,mod) in Base.loaded_modules
    if pkg.name in ("Revise","PrettyTables","Adapt")
        println(pkg.name,": ",Base.pkgversion(mod)," ",pathof(mod))
    end
end
Revise.revise()
Revise.errors()
PrettyTables.pretty_table([1 2;3 4])
println("Completed isolated Revise/load-path probe.")
