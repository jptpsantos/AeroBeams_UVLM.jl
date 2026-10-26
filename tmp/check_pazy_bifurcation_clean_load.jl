# Read-only startup check; including the guarded entry point runs no study.
using LinearAlgebra
BLAS.set_num_threads(1)
include(joinpath(@__DIR__, "..", "lib", "WingPropellerUVLM", "examples",
    "pazy_aerobeams_coupled", "run_pazy_bifurcation.jl"))
println("Bifurcation entry point loaded without executing the study.")
println("Julia version: ", VERSION)
println("Active project: ", Base.active_project())
for (pkg,mod) in Base.loaded_modules
    if pkg.name in ("Revise", "PrettyTables", "Adapt")
        println(pkg.name, " loaded, version=", Base.pkgversion(mod), "; source=", pathof(mod))
    end
end
@assert !any(pkg -> pkg.name == "Revise",keys(Base.loaded_modules))
println("PASS: clean process loaded the bifurcation module without Revise.")
