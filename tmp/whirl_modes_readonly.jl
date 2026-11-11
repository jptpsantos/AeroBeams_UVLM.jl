using DelimitedFiles,LinearAlgebra,Statistics,Printf
p = joinpath(pwd(),"lib","WingPropellerUVLM","examples","chang_linear_aeroelastic","output","multi_prop","single_propeller_history.txt")
a = readdlm(p,'\t',Float64;skipstart=5)
i = findall(x->3 <= x <= 4.9,a[:,1]); t = a[i,1] .- 3.95
f=3.9496;g=0.028
X=hcat(ones(length(t)),exp.(g.*t).*cos.(2π*f.*t),exp.(g.*t).*sin.(2π*f.*t))
y=a[i,2:5];c=X\y; phases=[atan(-c[3,k],c[2,k]) for k in 1:4]
for k in 1:4
 amp=hypot(c[2,k],c[3,k]); d=mod(rad2deg(phases[k]-phases[1])+180,360)-180
 @printf("col %d amplitude %.8g relative_phase_deg %.6f residual %.6f\n",k+1,amp,d,norm(y[:,k]-X*c[:,k])/norm(y[:,k].-mean(y[:,k])))
end
