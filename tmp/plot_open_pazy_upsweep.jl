using DelimitedFiles
import Plots

directory = normpath(joinpath(@__DIR__, "..", "lib", "WingPropellerUVLM", "examples",
    "pazy_aerobeams_coupled", "output", "LCO_downsweep_pulse"))
data, header = readdlm(joinpath(directory, "pazy_linear_LCO_upsweep.csv"), ',', Float64; header=true)
speed = data[:,findfirst(==("airspeed_mps"),vec(header))]
strain = data[:,findfirst(==("root_longitudinal_strain_ue"),vec(header))]
figure = Plots.plot(speed, strain; color=:blue, linewidth=0.5, linealpha=0.65,
    label="Upsweep", xlabel="Airspeed [m/s]", ylabel="Root longitudinal strain [με]",
    fontfamily="Computer Modern", guidefontsize=16, tickfontsize=12, legendfontsize=12,
    legend=:bottomleft, framestyle=:axes, grid=true, gridalpha=0.2,
    size=(800,550), dpi=300, margin=5Plots.mm)
path = joinpath(directory,"pazy_linear_LCO_root_strain_vs_airspeed_upsweep.png")
Plots.savefig(figure,path)
println(path)
