# Analyze wing and hub pitch from a Chang history with attachment-node twist.
# Usage: julia --project=lib/WingPropellerUVLM plot_hub_pitch_mode.jl [history.txt or .csv] [start_s] [end_s]
# Hub pitch is the combined intermediate y-rotation coordinate: attachment wing
# twist + relative pylon pitch. Global hub elevation additionally depends on
# the remaining finite rotations; the coordinate is a small-angle approximation
# to global pitch, not a general Euler-angle extraction.
include(joinpath(@__DIR__,"analyze_pitch_phase.jl"))

function read_hub_history(path)
    lines = readlines(path)
    header = findfirst(line -> startswith(strip(line),"time_s\t") || startswith(strip(line),"time_s,"),lines)
    isnothing(header) && error("Could not find a time_s header in $path")
    delim = occursin('\t',lines[header]) ? '\t' : ','
    columns = split(strip(lines[header]),delim)
    required = ["time_s","tip_twist_deg","propeller_1_attachment_wing_twist_deg","propeller_1_pitch_deg","propeller_1_yaw_deg"]
    indices = [findfirst(==(name),columns) for name in required]
    missing = required[findall(isnothing,indices)]
    isempty(missing) || error("Missing columns $(join(missing,", ")). Export a matching history with attachment twist first.")
    data = readdlm(IOBuffer(join(lines[header+1:end],"\n")),delim,Float64)[:,indices]
    all(isfinite,data) || error("History contains non-finite values")
    all(diff(data[:,1]) .> 0) || error("Time must be strictly increasing")
    return data
end

function hub_pitch_mode(data,start_time,end_time)
    t = data[:,1]
    first(t) <= start_time < end_time <= last(t) || error("Window outside history")
    selected = findall(x -> start_time <= x <= end_time,t)
    length(selected) >= 100 || error("Analysis window has fewer than 100 samples")
    reduced = selected[1:max(1,cld(length(selected),2500)):end]
    midpoint = (start_time+end_time)/2
    hub = data[:,3] .+ data[:,4]
    mode = fit_mode(t[reduced] .- midpoint,data[reduced,2],hub[reduced])
    signals = hcat(data[selected,2],data[selected,3],hub[selected],data[selected,4],data[selected,5])
    _, coefficients, prediction = fit_at(t[selected] .- midpoint,signals,mode.frequency,mode.growth)
    amplitudes = [hypot(coefficients[2,i],coefficients[3,i]) for i in 1:5]
    phase = [atan(-coefficients[3,i],coefficients[2,i]) for i in 1:5]
    relative = [mod(rad2deg(phase[i]-phase[1])+180,360)-180 for i in 1:5]
    residual = [norm(signals[:,i]-prediction[:,i])/norm(signals[:,i] .- mean(signals[:,i])) for i in 1:5]
    return (; midpoint,frequency=mode.frequency,growth=mode.growth,amplitudes,relative,residual)
end

function plot_hub_pitch_mode(path,start_time,end_time)
    mode = hub_pitch_mode(read_hub_history(path),start_time,end_time)
    labels = ["Wing-tip pitch", "Attachment wing pitch", "Combined hub pitch", "Relative pylon pitch", "Relative pylon yaw"]
    colors = ["#0072B2","#009E73","#D55E00","#CC79A7","#777777"]
    elapsed = range(0,2/mode.frequency,length=800)
    waves = plot(; xlabel="Time relative to wing-tip pitch peak (s)",ylabel="Normalized oscillation",
        title="Fitted pitch oscillations (mean and envelope removed)",ylim=(-1.15,1.15),
        legend=:bottomright,gridalpha=0.15)
    for i in 1:3
        plot!(waves,elapsed,cos.(2π*mode.frequency .* elapsed .+ deg2rad(mode.relative[i]));
            label=@sprintf("%s (%.3f°)",labels[i],mode.amplitudes[i]),color=colors[i],
            linewidth=2,linestyle=i==1 ? :solid : i==2 ? :dot : :dash)
    end
    amps = bar([1,2,3],mode.amplitudes[1:3];label="",color=colors[1:3],
        xticks=([1,2,3],["Wing tip","Attachment","Combined hub"]),
        title=@sprintf("Amplitude at %.2f s",mode.midpoint),ylabel="Amplitude (deg)",
        ylim=(0,maximum(mode.amplitudes[1:3])*1.25),gridalpha=0.15)
    for i in 1:3
        annotate!(amps,i,mode.amplitudes[i]+maximum(mode.amplitudes[1:3])*0.05,
            text(@sprintf("%.3f°",mode.amplitudes[i]),9,:center))
    end
    phase_limit=max(30.0,ceil(maximum(abs.(mode.relative[1:3]))/10)*10+10)
    ph=plot(;xlims=(-phase_limit,phase_limit),ylim=(0.5,3.5),
        yticks=([1,2,3],["Combined hub","Attachment wing","Wing tip"]),
        xlabel="Phase relative to wing-tip pitch (deg)",title="Relative phase",legend=false,gridalpha=0.15)
    vline!(ph,[0.0];label="",color=:gray,linestyle=:dot)
    for i in 1:3
        y=4-i
        plot!(ph,[0.0,mode.relative[i]],[y,y];label="",color=colors[i],linewidth=2)
        scatter!(ph,[mode.relative[i]],[y];label="",color=colors[i],markersize=5)
        annotate!(ph,mode.relative[i],y+0.18,text(@sprintf("%+.2f°",mode.relative[i]),9,:center))
    end
    figure=plot(waves,amps,ph;layout=@layout([a{0.50h}; b c]),size=(1400,950),margin=6Plots.mm,
        plot_title=@sprintf("Wing and hub pitch | %.3f Hz | %.1f–%.1f s",mode.frequency,start_time,end_time))
    output=joinpath(dirname(path),"pitch_phase_analysis")
    mkpath(output)
    for ext in ("png","pdf")
        target=joinpath(output,"hub_pitch_amplitude_phase."*ext)
        savefig(figure,target)
        println("Figure written to $target")
    end
    report=joinpath(output,"hub_pitch_summary.txt")
    open(report,"w") do io
        println(io,"Source: $(abspath(path))")
        @printf(io,"Window %.3f–%.3f s; reference time %.3f s\n",start_time,end_time,mode.midpoint)
        @printf(io,"Frequency %.5f Hz; growth rate %.5f 1/s\n",mode.frequency,mode.growth)
        println(io,"Positive phase leads the wing-tip pitch. Angles follow the structural sign convention.")
        for i in 1:5
            @printf(io,"%-22s amplitude %.6f deg; relative phase %+.3f deg; normalized residual %.5f\n",labels[i],mode.amplitudes[i],mode.relative[i],mode.residual[i])
        end
        @printf(io,"Combined hub / tip amplitude ratio: %.5f\n",mode.amplitudes[3]/mode.amplitudes[1])
        @printf(io,"Combined hub / attachment amplitude ratio: %.5f\n",mode.amplitudes[3]/mode.amplitudes[2])
        println(io,"Combined hub pitch is attachment wing twist + relative pylon pitch; this is the composed intermediate y-rotation coordinate.")
        println(io,"Global hub elevation depends on the other finite rotations; this coordinate approximates global pitch for small rotations.")
    end
    println(read(report,String))
end

if abspath(PROGRAM_FILE) == @__FILE__
    default_path=joinpath(@__DIR__,"output","multi_prop","single_propeller_history.txt")
    path=length(ARGS)>=1 ? ARGS[1] : default_path
    start_time=length(ARGS)>=2 ? parse(Float64,ARGS[2]) : 3.0
    end_time=length(ARGS)>=3 ? parse(Float64,ARGS[3]) : 4.9
    plot_hub_pitch_mode(path,start_time,end_time)
end
