using DelimitedFiles, LinearAlgebra, Statistics, Printf
path = joinpath(pwd(),"lib","WingPropellerUVLM","examples","chang_linear_aeroelastic","output","multi_prop","single_propeller_history.txt")
x = readdlm(path,'\t',Float64;skipstart=5)
rows = findall(t -> 3.0 <= t <= 4.9,x[:,1])
t = x[rows,1] .- 3.95
E = exp.(0.028 .* t)
A = hcat(ones(length(t)),E .* cos.(2π*3.9496 .* t),E .* sin.(2π*3.9496 .* t))
y = x[rows,[5,2,3,4]]
c = A \ y
phase = [atan(-c[3,i],c[2,i])*180/π for i in 1:4]
for i in 1:4
    amp=hypot(c[2,i],c[3,i]); delta=mod(phase[i]-phase[1]+180,360)-180
    err=norm(y[:,i]-A*c[:,i])/norm(y[:,i].-mean(y[:,i]))
    @printf("%s amplitude=%.8f phase=%.6f residual=%.6f mean=%.8f\n",["wing_twist_deg","pylon_pitch_deg","pylon_yaw_deg","tip_displacement_m"][i],amp,delta,err,c[1,i])
end
