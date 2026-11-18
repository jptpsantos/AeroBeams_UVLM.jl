using LinearAlgebra
BLAS.set_num_threads(1)
push!(ARGS,"--plots")
include(joinpath(@__DIR__,"..","lib","WingPropellerUVLM","examples",
    "pazy_aerobeams_coupled","test","test_static_spanwise_cn_alpha.jl"))
