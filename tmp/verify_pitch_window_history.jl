include(joinpath(pwd(),"lib","WingPropellerUVLM","examples","chang_linear_aeroelastic","plot_pitch_amplitude_phase_history.jl"))
using Test
const PH = PitchAmplitudePhaseHistory
pair = PH.read_pitch_pair(joinpath(pwd(),"tmp","hub_pitch_test","synthetic_history.txt"))
result = PH.pitch_window_history(pair)
expected = 0.1*cis(deg2rad(10.0)) + 2.0*cis(deg2rad(30.0))
@test pair.useHub
@test maximum(abs.(result.ratios .- abs(expected)/0.12)) < 0.005
@test maximum(abs.(result.phases .- rad2deg(angle(expected)))) < 0.02
@test maximum(result.residuals) < 0.01
println("Known amplitude ratio and phase recovered from every moving window; automatic attachment-pitch composition passed.")
PH.plot_pitch_amplitude_phase_history(joinpath(pwd(),"lib","WingPropellerUVLM","examples","chang_linear_aeroelastic","output","multi_prop","single_propeller_history.txt"))
