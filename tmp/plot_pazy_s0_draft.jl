include(joinpath(@__DIR__,"..","lib","WingPropellerUVLM","examples",
    "pazy_aerobeams_coupled","run_pazy_s0_steady_pitch_range.jl"))
plot_s0_steady_pitch_range(joinpath(@__DIR__,"pazy_s0_steady_draft"))
