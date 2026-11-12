using Test,DelimitedFiles,Random,LinearAlgebra
BLAS.set_num_threads(1)
Random.seed!(1234)
const ROOT = normpath(joinpath(@__DIR__,"..","..",".."))
runExamplePlots = false
SELFatol = 1e-6
SELFrtol = 1e-4
include(joinpath(ROOT,"test","testData.jl"))
include_example("examples/eigenAeroelastic/heliosFlutterPLambdaRange.jl")
@testset "Isolated upstream Helios reference comparison" begin
    for (name,actual) in (("freqs",freqsOsc),("damps",dampsOsc))
        reference = readdlm(joinpath(ROOT,"test","newTestDataGenerators","heliosFlutterPLambdaRange",name*".txt"),Float64)
        values = test_data_matrix(actual)
        println(name,": actual shape=",size(values),", reference shape=",size(reference))
        if size(values) == size(reference)
            println("NaN mask differences: ",count(isnan.(values) .!= isnan.(reference)))
            mask = isfinite.(values) .& isfinite.(reference)
            println("Maximum finite absolute difference: ",maximum(abs.(values[mask] .- reference[mask])))
        end
        @test matches_test_data(actual,"heliosFlutterPLambdaRange",name*".txt")
    end
end
