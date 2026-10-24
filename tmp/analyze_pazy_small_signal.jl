using DelimitedFiles, LinearAlgebra, Statistics, Printf
import Plots
const BASE = normpath(joinpath(@__DIR__,"..","lib","WingPropellerUVLM","examples","pazy_aerobeams_coupled"))
const OUT = joinpath(BASE,"output","LCO_small_signal_audit")
readcase(name) = first(readdlm(joinpath(OUT,name*".csv"),',',Float64;header=true))
tags = [tag for tag in ("default","tight","refined") if isfile(joinpath(OUT,tag*"_larger.csv"))]
figures = []
for (size,force) in (("small",1e-5),("larger",.01))
    figure = Plots.plot(;xlabel="Time [s]",ylabel="Strain / coefficient [microstrain/N]",
        title="Pulse coefficient = $force N; unforced trajectory subtracted",
        titlefontsize=10,legend=:outerright,guidefontsize=10,margin=6Plots.mm)
    Plots.vspan!(figure,[.01,.05];color=:gray,alpha=.1,label="Pulse interval")
    for (tag,color) in zip(("default","tight","refined"),(:orange,:blue,:green))
        tag in tags || continue
        quiet,forced = readcase(tag*"_quiet"),readcase(tag*"_"*size)
        @assert quiet[:,1] == forced[:,1]
        delta = (forced[:,2] .- quiet[:,2]) ./ force
        Plots.plot!(figure,forced[:,1],delta;label=tag,color,linewidth=1.3)
    end
    push!(figures,figure)
end
Plots.savefig(Plots.plot(figures...;layout=(2,1),size=(1000,800),dpi=180),
    joinpath(OUT,"normalized_pulse_response.png"))
open(joinpath(OUT,"linearity_checks.txt"),"w") do io
    for tag in tags
        quiet,small,large = readcase(tag*"_quiet"),readcase(tag*"_small"),readcase(tag*"_larger")
        a = (small[:,2] .- quiet[:,2]) ./ 1e-5
        b = (large[:,2] .- quiet[:,2]) ./ .01
        idx = findall(t -> .06 <= t <= .49,quiet[:,1])
        @printf(io,"%s: normalized small/large relative L2 mismatch %.8g percent; max quiet drift %.8g ue\n",
            tag,100norm(a[idx]-b[idx])/norm(b[idx]),maximum(abs.(quiet[:,2] .- quiet[1,2])))
        @printf(io,"  Mean coupling iterations quiet/small/large: %.6f / %.6f / %.6f\n",
            mean(quiet[2:end,5]),mean(small[2:end,5]),mean(large[2:end,5]))
    end
end
println(read(joinpath(OUT,"linearity_checks.txt"),String))
snapshot = joinpath(OUT,"matched_upsweep_snapshot.csv")
if isfile(snapshot)
    include(joinpath(BASE,"plot_pazy_LCO.jl"))
    data,header = readdlm(snapshot,',',Float64;header=true)
    speed = data[:,findfirst(==("airspeed_mps"),vec(header))]
    strain = data[:,findfirst(==("root_longitudinal_strain_ue"),vec(header))]
    figure = pazy_LCO_strain_plot([(;airspeed=speed,strain,
        label="Unforced upsweep (latest pulse run)",color=:blue)])
    Plots.savefig(figure,joinpath(OUT,"latest_unforced_upsweep.png"))
end
