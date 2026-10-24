# Plot saved histories that match the ramp durations in the latest settings.
# No simulation is run. Older branch responses are reported and omitted.
using DelimitedFiles, SHA
import Plots

function plot_pazy_LCO_latest(;
    source_directory=joinpath(@__DIR__,"output","LCO_downsweep_pulse"),
    output_directory=joinpath(@__DIR__,"output","LCO_latest_configured"))
    settings = read(joinpath(source_directory,"pazy_linear_LCO_settings.txt"),String)
    branches = []
    status = String[]
    missing = String[]
    for (name,heading,color) in (("upsweep","Upsweep",:blue),("downsweep","Downsweep",:red))
        pattern = Regex(heading * ": ([0-9.eE+-]+) -> ([0-9.eE+-]+) m/s in ([0-9.eE+-]+) s")
        ramp = match(pattern,settings)
        isnothing(ramp) && error("Missing $heading ramp in settings")
        start_speed,end_speed,duration = parse.(Float64,ramp.captures)
        path = joinpath(source_directory,"pazy_linear_LCO_$name.csv")
        if !isfile(path)
            push!(missing,"$heading ($(duration) s): response not saved")
            push!(status,"$heading: no saved response at $path")
            continue
        end
        bytes = read(path)
        data,header = readdlm(IOBuffer(bytes),',',Float64;header=true)
        column(key) = data[:,something(findfirst(==(key),vec(header)))]
        time,speed,strain = column("time_s"),column("airspeed_mps"),column("root_longitudinal_strain_ue")
        metadata = "$heading source: $(abspath(path))\nSHA256: $(bytes2hex(sha256(bytes)))\n" *
            "Saved time: $(first(time)) -> $(last(time)) s; speed: $(first(speed)) -> $(last(speed)) m/s"
        if !(isapprox(first(time),0;atol=1e-8,rtol=0) &&
             isapprox(last(time),duration;atol=1e-8,rtol=0) &&
             isapprox(first(speed),start_speed;atol=1e-8,rtol=0) &&
             isapprox(last(speed),end_speed;atol=1e-8,rtol=0))
            push!(missing,"$heading ($(duration) s): no matching saved response")
            push!(status,metadata * "\nOMITTED: saved endpoints/duration differ from latest settings.")
            continue
        end
        all(isfinite,strain) && all(diff(time) .> 0) || error("Invalid history in $path")
        push!(branches,(;speed,strain,label="$heading ($(duration) s)",color))
        push!(status,metadata * "\nINCLUDED: ramp endpoints and duration match settings.")
    end
    isempty(branches) && error("No saved histories match the latest configured ramps")
    title = isempty(missing) ? "Latest configured case" : join(missing,"\n")
    figure = Plots.plot(;xlabel="Airspeed [m/s]",ylabel="Root longitudinal strain [microstrain]",
        fontfamily="Computer Modern",guidefontsize=16,tickfontsize=12,
        legendfontsize=11,legend=:bottomleft,framestyle=:axes,grid=true,gridalpha=.2,
        size=(800,550),dpi=300,margin=5Plots.mm,title,titlefontsize=10)
    for branch in branches
        Plots.plot!(figure,branch.speed,branch.strain;label=branch.label,
            color=branch.color,linewidth=.5,linealpha=.65)
    end
    mkpath(output_directory)
    path = joinpath(output_directory,"pazy_LCO_latest_root_strain_vs_airspeed.png")
    Plots.savefig(figure,path)
    write(joinpath(output_directory,"settings_snapshot.txt"),settings)
    write(joinpath(output_directory,"plot_status.txt"),join(status,"\n\n") *
        "\n\nEndpoint/duration checks prevent use of older ramps; they do not establish mesh or model provenance when those values are absent from a CSV.\n")
    println(join(status,"\n\n"))
    println("Saved plot: ",path)
    return figure
end

if abspath(PROGRAM_FILE) == @__FILE__
    plot_pazy_LCO_latest()
end
