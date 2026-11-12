using Test,DelimitedFiles
const ROOT = normpath(joinpath(@__DIR__,"..","..",".."))
const ANALYSIS = joinpath(ROOT,"lib","WingPropellerUVLM","examples","chang_linear_aeroelastic")
include(joinpath(ANALYSIS,"analyze_pitch_phase.jl"))
result = run_pitch_reviewer_analysis(; output_directory=joinpath(@__DIR__,"pitch_analysis_qa"))
reference = readdlm(joinpath(ANALYSIS,"output","multi_prop","pitch_phase_analysis",
    "reviewer_response","pitch_mode_results.csv"),','; header=true)
actual = readdlm(joinpath(result.output_directory,"pitch_mode_results.csv"),','; header=true)
@testset "Saved-history pitch analysis preserved" begin
    @test actual[2] == reference[2]
    @test actual[1][:,1] == reference[1][:,1]
    @test isapprox(Float64.(actual[1][:,2:13]),Float64.(reference[1][:,2:13]); atol=1e-8,rtol=1e-8)
    @test actual[1][:,14:end] == reference[1][:,14:end]
end
