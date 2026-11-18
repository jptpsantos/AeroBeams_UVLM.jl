using LinearAlgebra
BLAS.set_num_threads(1)
for file in ("test_pazy_model_options.jl","test_static_uvlm.jl")
    scope = Module(gensym(:StaticSpanwiseRegression))
    Core.eval(scope,:(include(path::AbstractString) = Base.include($scope,path)))
    Base.include(scope,joinpath(@__DIR__,"..","lib","WingPropellerUVLM","examples",
        "pazy_aerobeams_coupled","test",file))
end
