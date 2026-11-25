# Plot the whirl-like mode represented by a saved single-propeller history.
# Propeller coordinates are relative pylon angles, not absolute hub angles.
include(joinpath(@__DIR__,"analyze_pitch_phase.jl"))

function plot_whirl_mode(path,start_time,end_time)
    lines = readlines(path)
    header = findfirst(s -> startswith(strip(s),"time_s\t"),lines)
    isnothing(header) && error("Missing tab-separated time_s header")
    names = split(strip(lines[header]),'\t')
    needed = ["time_s","tip_twist_deg","propeller_1_pitch_deg","propeller_1_yaw_deg","tip_displacement_m"]
    indexes = [findfirst(==(name),names) for name in needed]
    any(isnothing,indexes) && error("Missing required history columns")
    all_data = readdlm(IOBuffer(join(lines[header+1:end],"\n")),'\t',Float64)
    data = all_data[:,indexes]
    all(isfinite,data) || error("Non-finite data")
    t = data[:,1]
    all(diff(t) .> 0) || error("Time must be strictly increasing")
    first(t) <= start_time < end_time <= last(t) || error("Window outside history")
    inds = findall(x -> start_time <= x <= end_time,t)
    length(inds) >= 100 || error("Analysis window has fewer than 100 samples")
    reduced = inds[1:max(1,cld(length(inds),2500)):end]
    midpoint = (start_time + end_time)/2
    mode = fit_mode(t[reduced] .- midpoint,data[reduced,2],data[reduced,3])
    _, coefficients, prediction = fit_at(t[inds] .- midpoint,data[inds,2:5],mode.frequency,mode.growth)
    amplitudes = [hypot(coefficients[2,i],coefficients[3,i]) for i in 1:4]
    phases = [atan(-coefficients[3,i],coefficients[2,i]) for i in 1:4]
    relative = [mod(rad2deg(phases[i]-phases[1])+180,360)-180 for i in 1:4]
    residuals = [norm(data[inds,i+1]-prediction[:,i])/norm(data[inds,i+1] .- mean(data[inds,i+1])) for i in 1:4]
    labels = ["Wing-tip twist", "Relative pylon pitch", "Relative pylon yaw", "Tip disp. (+down)"]
    colors = ["#0072B2", "#D55E00", "#009E73", "#CC79A7"]
    displays = [@sprintf("%.3f°",amplitudes[1]),@sprintf("%.3f°",amplitudes[2]),
        @sprintf("%.3f°",amplitudes[3]),@sprintf("%.2f mm",1000amplitudes[4])]
    period = 1/mode.frequency
    elapsed = range(0,2period,length=900)
    waves = plot(;xlabel="Time relative to wing-twist peak (s)", ylabel="Normalized oscillation",
        title="Fitted oscillations (mean and envelope removed)", gridalpha=0.15,
        ylims=(-1.15,1.15),legend=:bottomright)
    for i in 1:4
        plot!(waves,elapsed,cos.(2π*mode.frequency .* elapsed .+ deg2rad(relative[i]));
            label="$(labels[i]) ($(displays[i]))", color=colors[i],linewidth=2,
            linestyle=i == 1 ? :solid : i == 2 ? :dash : i == 3 ? :dot : :dashdot)
    end
    cycle = range(0,period,length=500)
    px = amplitudes[2] .* cos.(2π*mode.frequency .* cycle .+ deg2rad(relative[2]))
    py = amplitudes[3] .* cos.(2π*mode.frequency .* cycle .+ deg2rad(relative[3]))
    orbit = plot(px,py; xlabel="Relative pylon pitch (deg)",ylabel="Relative pylon yaw (deg)",
        title="Pylon pitch–yaw orbit",label="",color="#333333",linewidth=2,
        aspect_ratio=:equal,gridalpha=0.15)
    scatter!(orbit,[px[1]],[py[1]];label="Wing-twist peak",color=colors[1],markersize=5)
    k=60
    plot!(orbit,px[k:k+3],py[k:k+3];label="",color=colors[2],linewidth=3,arrow=:closed)
    phase_plot = plot(; xlims=(-190,190),ylims=(0.4,4.6),
        yticks=([1,2,3,4],reverse(labels)),xlabel="Phase relative to tip twist (deg)",
        title="Relative phase of fitted mode",gridalpha=0.15,legend=false)
    vline!(phase_plot,[0];color=:gray,linestyle=:dot,label="")
    for i in 1:4
        y=5-i
        plot!(phase_plot,[0,relative[i]],[y,y];color=colors[i],linewidth=2,label="")
        scatter!(phase_plot,[relative[i]],[y];color=colors[i],markersize=5,label="")
        annotate!(phase_plot,relative[i],y+0.25,text(@sprintf("%+.1f°",relative[i]),9,:center))
    end
    figure=plot(waves,orbit,phase_plot;layout=@layout([a{0.50h}; b c]),size=(1400,950),
        margin=6Plots.mm,plot_title=@sprintf("Whirl-like response | %.3f Hz | %.1f–%.1f s",mode.frequency,start_time,end_time))
    output=joinpath(dirname(path),"pitch_phase_analysis")
    mkpath(output)
    for extension in ("png","pdf")
        target=joinpath(output,"whirl_mode."*extension)
        savefig(figure,target)
        println("Figure written to $target")
    end
    report=joinpath(output,"whirl_mode_summary.txt")
    open(report,"w") do io
        println(io,"Source: $(abspath(path))")
        @printf(io,"Window %.3f–%.3f s; amplitude reference time %.3f s\n",start_time,end_time,midpoint)
        @printf(io,"Frequency %.5f Hz; common growth rate %.5f 1/s\n",mode.frequency,mode.growth)
        println(io,"Relative phases use wing-tip twist as 0°. Positive phase leads; angles follow saved structural sign conventions.")
        for i in 1:4
            @printf(io,"%-26s amplitude %s, phase %+.3f deg, normalized residual %.5f\n",labels[i],displays[i],relative[i],residuals[i])
        end
        @printf(io,"Pylon yaw phase minus pylon pitch phase: %+.3f deg\n",mod(relative[3]-relative[2]+180,360)-180)
        println(io,"Pitch and yaw are pylon angles relative to the local wing; wing vertical displacement follows the structural down-positive convention.")
        println(io,"The attachment-node rotation and aerodynamic loads are not in this file. This mode shape does not by itself establish which aerodynamic contribution drives the instability.")
    end
    println(read(report,String))
end

if abspath(PROGRAM_FILE) == @__FILE__
    default_path=joinpath(@__DIR__,"output","multi_prop","single_propeller_history.txt")
    path=length(ARGS)>=1 ? ARGS[1] : default_path
    start_time=length(ARGS)>=2 ? parse(Float64,ARGS[2]) : 3.0
    end_time=length(ARGS)>=3 ? parse(Float64,ARGS[3]) : 4.9
    plot_whirl_mode(path,start_time,end_time)
end
