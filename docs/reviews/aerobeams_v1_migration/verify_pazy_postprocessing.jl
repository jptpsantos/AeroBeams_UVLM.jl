using Test
const ROOT = normpath(joinpath(@__DIR__,"..","..",".."))
const TEST_DIRECTORY = joinpath(ROOT,"lib","WingPropellerUVLM","examples",
    "pazy_aerobeams_coupled","test")
@testset "Pazy FFTW-dependent postprocessing" begin
    for file in ("test_pazy_stability_sweep.jl","test_pazy_LCO_spectrogram.jl")
        scope = Module(gensym(:PazyPostprocessing))
        Core.eval(scope,:(include(path::AbstractString) = Base.include($scope,path)))
        Base.include(scope,joinpath(TEST_DIRECTORY,file))
    end
end
