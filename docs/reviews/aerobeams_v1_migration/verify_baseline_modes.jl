using Test,LinearAlgebra
const ROOT = normpath(joinpath(@__DIR__,"..","..",".."))
const EXAMPLE = joinpath(ROOT,"lib","WingPropellerUVLM","examples","chang_aerobeams_coupled")
include(joinpath(EXAMPLE,"src","ChangAeroBeams.jl"))
const CA = ChangAeroBeams
BLAS.set_num_threads(1)

# The committed pre-migration log records the independent AeroBeams 0.8.1
# result, rather than recomputing the expected frequencies with v1.
bytes = read(joinpath(EXAMPLE,"validation","modal_comparison.log"))
baseline = bytes[1:2] == [0xff,0xfe] ? transcode(String,reinterpret(UInt16,bytes[3:end])) : String(bytes)
matches = collect(eachmatch(r"(?s)Ne=(\d+): Chang=\[.*?\], AeroBeams=\[([^\]]+)\]",baseline))
@testset "AeroBeams 0.8.1 versus v1 dry-wing frequencies" begin
    @test length(matches) == 3
    for match in matches
        elements = parse(Int,match.captures[1])
        expected = parse.(Float64,strip.(split(match.captures[2],',')))
        config = CA.chang_aerobeams_defaults()
        defaults = merge(config.reference,(
            wing=merge(config.reference.wing,(spanwise_panels=elements,chordwise_panels=2)),
            simulation=merge(config.reference.simulation,(end_time_s=0.002,))))
        reference = CA.Chang.build_chang_model(CA.Chang.load_chang_configuration(defaults;
            env=Dict{String,String}()))
        wing = CA.build_wing(reference)
        problem = CA.AB.create_EigenProblem(model=deepcopy(wing.problem.model),
            systemSolver=CA.AB.create_LinearSolver(),eigenmethod=:fulldecomp)
        CA.AB.solve!(problem)
        frequencies = sort(imag.(filter(λ -> imag(λ)>2pi*0.01,problem.eigensolution.λ)))[1:6] ./ (2pi)
        println("Ne=$elements: baseline=$expected, v1=$frequencies")
        println("Maximum relative frequency difference: ",maximum(abs.(frequencies ./ expected .- 1)))
        @test frequencies ≈ expected rtol=1e-5 atol=1e-7
    end
end
