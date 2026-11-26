example = joinpath(@__DIR__,"..","lib","WingPropellerUVLM","examples","pazy_aerobeams_coupled")
include(joinpath(example,"run_pazy_s0_steady_pitch_range.jl"))
directory = joinpath(example,"output","s0_fig2_comparison")
plot_s0_steady_pitch_range(joinpath(example,"output","s0_steady_pitch_range");
    mesh=:nominal,modelLabel="S0 beam-rod model",
    figureReference=joinpath(directory,"figure2_aerobeams_samples.csv"),saveDirectory=directory)
