include(joinpath(@__DIR__,"..","lib","WingPropellerUVLM","examples",
    "pazy_aerobeams_coupled","run_pazy_s0_steady_pitch_range.jl"))
LinearAlgebra.BLAS.set_num_threads(1)
directory = joinpath(@__DIR__,"pazy_s0_steady_resume")
failures = run_s0_steady_pitch_range(outputDirectory=directory)
plot_s0_steady_pitch_range(directory)
iszero(failures) || error("$failures UVLM conditions did not converge")
