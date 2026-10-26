# Use the saved 40 s upsweep and 20 s downsweep. No simulations are executed.
base = normpath(joinpath(@__DIR__,"..","lib","WingPropellerUVLM","examples","pazy_aerobeams_coupled"))
path = joinpath(base,"plot_pazy_LCO.jl")
source = read(path,String)
# The user's plot script ends with an unconditional default plotting call.
# Load only its definitions to avoid regenerating unrelated case directories.
definitions = first(split(source,"\n#if abspath(PROGRAM_FILE)";limit=2))
Base.include_string(Main,definitions,path)
directory = joinpath(base,"output","LCO_downsweep_pulse")
for (branch,duration) in (("upsweep",40.),("downsweep",20.))
    data,header = readdlm(joinpath(directory,"pazy_linear_LCO_$branch.csv"),',',Float64;header=true)
    time = data[:,findfirst(==("time_s"),vec(header))]
    @assert isapprox(last(time),duration;atol=1e-8,rtol=0)
    @assert all(isfinite,data)
end
for selection in (:both,:upsweep,:downsweep)
    figure = plot_pazy_LCO(;output_directory=directory,plot_branches=selection,linealpha=1.0)
    @assert length(figure.series_list) == (selection == :both ? 2 : 1)
    @assert all(series -> series[:linealpha] == 1.0,figure.series_list)
end
println("Verified fully opaque lines for the latest saved 40 s / 20 s responses.")
