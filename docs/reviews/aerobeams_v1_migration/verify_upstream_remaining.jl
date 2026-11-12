# Run suites after the upstream Helios comparison, which stops runtests.jl early.
using Test,LinearAlgebra
BLAS.set_num_threads(1)
const ROOT = normpath(joinpath(@__DIR__,"..","..",".."))
@testset "Remaining unmodified upstream suites" begin
    for file in ("aerodynamicTests.jl","trimAeroelasticTests.jl","steadyAeroelasticTests.jl",
        "dynamicStructuralTests.jl","eigenStructuralTests.jl","staticStructuralTests.jl",
        "trimStructuralTests.jl")
        @testset "$file" begin
            scope = Module(gensym(:UpstreamSuite))
            Core.eval(scope,:(include(path::AbstractString) = Base.include($scope,path)))
            Core.eval(scope,quote
                using Test,DelimitedFiles
                runExamplePlots = false
                is_ci = false
                SELFatol = 1e-6
                SELFrtol = 1e-4
            end)
            Base.include(scope,joinpath(ROOT,"test","testData.jl"))
            Base.include(scope,joinpath(ROOT,"test",file))
        end
    end
end
