const PAZY_EXAMPLE_DIR = normpath(joinpath(@__DIR__,".."))

using Test
include(joinpath(PAZY_EXAMPLE_DIR,"PazyWingUVLMCoupling.jl"))
include(joinpath(PAZY_EXAMPLE_DIR,"PazyWingUVLMVisualization.jl"))

@testset "AeroBeams animation output paths" begin
    root = abspath(REPOSITORY_ROOT)
    @test path_for_aerobeams(joinpath(root,"output","preview.gif")) == "/output/preview.gif"
    @test_throws ArgumentError path_for_aerobeams(joinpath(dirname(root),"preview.gif"))
    @test_throws ArgumentError path_for_aerobeams(joinpath(root * "_other","preview.gif"))
    @test_throws ArgumentError path_for_aerobeams(joinpath(root,"output","..","..","preview.gif"))
end
