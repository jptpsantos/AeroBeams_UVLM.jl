push!(ARGS,"--dynamic","--animations")
include(joinpath(@__DIR__,"..","lib","WingPropellerUVLM","examples",
    "pazy_aerobeams_coupled","test","runtests.jl"))
scope = Module(:PazyStaticRegression)
Core.eval(scope,:(include(path::AbstractString) = Base.include($scope,path)))
Base.include(scope,joinpath(@__DIR__,"..","lib","WingPropellerUVLM","examples",
    "pazy_aerobeams_coupled","test","test_static_uvlm.jl"))
