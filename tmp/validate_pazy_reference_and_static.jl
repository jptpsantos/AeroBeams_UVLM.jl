using LinearAlgebra
BLAS.set_num_threads(1)
push!(ARGS,"--plots")
for file in ("test_pazy_stability_sweep.jl","test_static_uvlm.jl")
    scope = Module(gensym(:PazyFinalRegression))
    Core.eval(scope,:(include(path::AbstractString) = Base.include($scope,path)))
    Base.include(scope,joinpath(@__DIR__,"..","lib","WingPropellerUVLM","examples",
        "pazy_aerobeams_coupled","test",file))
end
