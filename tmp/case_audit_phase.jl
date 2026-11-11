using DelimitedFiles,LinearAlgebra,Statistics,Printf
path=joinpath(pwd(),"lib","WingPropellerUVLM","examples","chang_linear_aeroelastic","output","multi_prop","single_propeller_history.txt")
d = readdlm(path,'\t',Float64;skipstart=5)
sel = (3 .<= d[:,1]) .& (d[:,1] .<= 4.9)
t = d[sel,1] .- 3.95
f = 3.9496; sigma = 0.028
H = hcat(ones(length(t)),exp.(sigma.*t).*cos.(2π*f.*t),exp.(sigma.*t).*sin.(2π*f.*t))
D = d[sel,2:5]
coef = H \ D
for i in 1:4
 amp = hypot(coef[2,i],coef[3,i]); phase=rad2deg(atan(-coef[3,i],coef[2,i])); res=norm(D[:,i]-H*coef[:,i])/norm(D[:,i].-mean(D[:,i]))
 @printf("i=%d amp=%.9f phase=%.6f residual=%.6f\n",i,amp,phase,res)
end
