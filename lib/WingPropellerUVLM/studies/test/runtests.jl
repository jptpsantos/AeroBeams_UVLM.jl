# This project declares the study, plotting, FFT and test dependencies.
empty!(LOAD_PATH)
append!(LOAD_PATH,["@","@stdlib"])
include(joinpath(@__DIR__,"..","..","test","migration_api.jl"))
include(joinpath(@__DIR__,"..","..","test","runtests.jl"))
for file in ("aerodynamic_convergence.jl","core_convergence.jl",
    "independent_convergence.jl","moving_block_damping.jl")
    include(joinpath(@__DIR__,"..","..","test",file))
end
include("migration_smoke.jl")
