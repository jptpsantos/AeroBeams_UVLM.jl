include(joinpath(pwd(),"lib","WingPropellerUVLM","examples","chang_linear_aeroelastic","analyze_pitch_phase.jl"))
using Test
x = collect(range(-1.0,1.0; length=1600))
w = 0.2 .+ 0.12 .* exp.(0.08 .* x) .* cos.(2π * 4.2 .* x .- 0.3)
p = -0.4 .+ 4.5 .* exp.(0.08 .* x) .* cos.(2π * 4.2 .* x .+ 0.4)
m = fit_mode(x,w,p)
@test abs(m.frequency-4.2) < 0.005
@test abs(m.growth-0.08) < 0.01
@test abs(m.delta-rad2deg(0.7)) < 0.2
@test abs(m.amplitude[2]/m.amplitude[1]-37.5) < 0.1
println("Synthetic recovery passed: ", (m.frequency,m.growth,m.delta))
path = joinpath(pwd(),"lib","WingPropellerUVLM","examples","chang_linear_aeroelastic","output","multi_prop","single_propeller_history.txt")
t, w, p = read_history(path)
for (a,b) in [(2.0,3.4),(3.5,4.9)]
    inds = findall(z -> a <= z <= b,t)[1:4:end]
    m = fit_mode(t[inds] .- (a+b)/2,w[inds],p[inds])
    println("window=",(a,b)," frequency=",m.frequency," growth=",m.growth," phase_deg=",m.delta," amplitude_ratio=",m.amplitude[2]/m.amplitude[1]," residuals=",m.residuals)
end
