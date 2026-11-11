using LinearAlgebra,Test,DelimitedFiles
include(joinpath(@__DIR__,"..","src","ChangAeroBeams.jl"))
using .ChangAeroBeams
const CA = ChangAeroBeams
BLAS.set_num_threads(1)
rows = Any[]
@testset "Chang structural data and modal comparison" begin
    for elements in (8,16,32)
        c = chang_aerobeams_defaults()
        d = merge(
            c.reference,
            (
                wing=merge(c.reference.wing,(spanwise_panels=elements,chordwise_panels=2)),
                simulation=merge(c.reference.simulation,(end_time_s=0.002,)))
        )
        ref = CA.Chang.build_chang_model(CA.Chang.load_chang_configuration(d; env=Dict{String,String}()))
        wing = build_wing(ref)
        eigenproblem = CA.AB.create_EigenProblem(model=deepcopy(wing.problem.model),
            getLinearSolution=true,nModes=6,frequencyFilterLimits=[0.01,Inf])
        CA.AB.solve!(eigenproblem)
        af = eigenproblem.frequenciesOscillatory ./ (2pi)
        cf = CA.Chang.chang_wing_modal_frequencies(ref.structural; number_of_modes=6)
        @test length(af) == length(cf) == 6
        @test wing.metadata.spatial_block_conservation_error <= 1e-12
        @test wing.metadata.mass_kg ≈ sum(ref.parameters.m_node_vec) rtol = 1e-12
        for mode = 1:6
            push!(rows,[elements,mode,cf[mode],af[mode],(af[mode] / cf[mode] - 1) * 100])
        end
        println("Ne=$elements: Chang=$cf, AeroBeams=$af")
    end
end
directory = joinpath(@__DIR__,"..","output")
mkpath(directory)
open(joinpath(directory,"wing_modal_comparison.csv"),"w") do io
    println(io,"elements,mode,chang_hz,aerobeams_hz,difference_percent")
    for row in rows
        writedlm(io,permutedims(row),',')
    end
end
