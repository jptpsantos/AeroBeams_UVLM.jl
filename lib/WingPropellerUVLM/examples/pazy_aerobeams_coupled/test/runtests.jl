using Test

# Keep standalone test options available, and restore ARGS when included from an IDE.
originalArgs = copy(ARGS)
fullSuite = "--full" in ARGS
testFiles = [
    "test_coupling.jl",
    "verify_pazy_coupling.jl",
    "test_visualization.jl",
    "test_pazy_stability_sweep.jl",
    "test_pazy_LCO_spectrogram.jl",
    "test_pazy_LCO_output.jl",
    "test_pazy_bifurcation.jl",
    "test_pazy_bifurcation_plots.jl",
    "test_pazy_stepped_bifurcation.jl",
]
if fullSuite
    "--dynamic" in ARGS || push!(ARGS,"--dynamic")
    append!(testFiles,[
        "test_static_uvlm.jl",
        "test_pazy_LCO.jl",
        "test_pazy_continuous_LCO.jl",
        "test_pazy_LCO_diagnostic.jl",
    ])
end

try
    @testset "Pazy UVLM" begin
        for file in testFiles
            @testset "$file" begin
                println("Running ",file)
                # Each standalone script defines its own model helpers and constants.
                scope = Module(gensym(:PazyTest))
                Core.eval(scope,:(include(path::AbstractString) = Base.include($scope,path)))
                Base.include(scope,joinpath(@__DIR__,file))
            end
        end
    end
finally
    empty!(ARGS)
    append!(ARGS,originalArgs)
end
