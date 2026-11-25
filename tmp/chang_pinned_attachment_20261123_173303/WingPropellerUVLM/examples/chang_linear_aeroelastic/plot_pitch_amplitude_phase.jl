# Visual comparison of the quantities in the saved single-propeller history.
# Usage: julia --project=lib/WingPropellerUVLM plot_pitch_amplitude_phase.jl [history.txt] [start_s] [end_s]
# Propeller pitch in this file is pylon rotation relative to the local wing.
include(joinpath(@__DIR__,"analyze_pitch_phase.jl"))

function plot_pitch_amplitude_phase(path,start_time,end_time)
    t, wing, pylon = read_history(path)
    first(t) <= start_time < end_time <= last(t) || error("Window outside history")
    selected = findall(x -> start_time <= x <= end_time,t)
    length(selected) >= 100 || error("Analysis window has fewer than 100 samples")
    stride = max(1,cld(length(selected),2500))
    reduced = selected[1:stride:end]
    midpoint = (start_time+end_time)/2
    mode = fit_mode(t[reduced] .- midpoint,wing[reduced],pylon[reduced])
    _, coefficients, _ = fit_at(t[selected] .- midpoint,
        hcat(wing[selected],pylon[selected]),mode.frequency,mode.growth)
    amplitudes = [hypot(coefficients[2,i],coefficients[3,i]) for i in 1:2]
    phase = mod(rad2deg(atan(-coefficients[3,2],coefficients[2,2]) -
        atan(-coefficients[3,1],coefficients[2,1]))+180,360)-180
    lead = phase / (360 * mode.frequency)
    direction = phase >= 0 ? "lead" : "lag"
    phase_limit = max(30.0,ceil(abs(phase)/10)*10+10)
    time = range(-0.13/mode.frequency * 4,0.13/mode.frequency * 4,length=800)
    wing_signal = cos.(2π * mode.frequency .* time)
    pylon_signal = cos.(2π * mode.frequency .* time .+ deg2rad(phase))
    blue, orange = "#0072B2", "#D55E00"
    waveform = plot(time,wing_signal; label="Wing-tip twist", color=blue, linewidth=2.6,
        xlabel="Time relative to wing-tip peak (s)", ylabel="Normalized oscillation",
        title="Fitted oscillations (offset and envelope removed)", legend=:bottomright,
        ylims=(-1.28,1.32), gridalpha=0.18)
    plot!(waveform,time,pylon_signal; label="Relative pylon pitch", color=orange,
        linewidth=2.3, linestyle=:dash)
    vline!(waveform,[0.0]; label="",color=blue,linestyle=:dot,linewidth=1.1)
    vline!(waveform,[-lead]; label="",color=orange,linestyle=:dot,linewidth=1.1)
    annotate!(waveform,0.04,1.16,text(@sprintf("Phase difference = %+.2f°  (%.2f ms %s)",phase,1000abs(lead),direction),9,:left))
    amplitude_plot = bar([1,2],amplitudes; color=[blue,orange],label="",bar_width=0.52,
        xticks=([1,2],["Wing-tip\ntwist","Relative pylon\npitch"]),
        ylabel="Amplitude (deg)", title=@sprintf("Amplitude at %.2f s",midpoint),
        ylim=(0,maximum(amplitudes)*1.24), gridalpha=0.18)
    for i in 1:2
        annotate!(amplitude_plot,i,amplitudes[i]+maximum(amplitudes)*0.055,
            text(@sprintf("%.3f°",amplitudes[i]),10,:center))
    end
    phase_plot = plot([-phase_limit,phase_limit],[0.0,0.0]; label="",color=:gray,linewidth=1.5,
        xlims=(-phase_limit,phase_limit),ylims=(-0.6,0.8),yticks=false,
        xlabel="Pylon phase − wing-tip phase (deg)", title="Phase difference",
        grid=false, framestyle=:box)
    scatter!(phase_plot,[0.0],[0.0]; label="Wing tip: 0°",color=blue,markersize=6)
    scatter!(phase_plot,[phase],[0.0]; label=@sprintf("Relative pylon: %+.2f°",phase),
        color=orange,markersize=6,legend=:topright)
    plot!(phase_plot,[0,phase],[0,0]; label="",color=orange,linewidth=3)
    annotate!(phase_plot,phase,0.19,text(@sprintf("%+.2f°",phase),10,:center))
    layout = @layout [a{0.54h}; b c]
    figure = plot(waveform,amplitude_plot,phase_plot; layout=layout,size=(1200,850),
        margin=6Plots.mm, plot_title=@sprintf("Wing-tip and relative pylon pitch | %.3f Hz | %.1f–%.1f s",mode.frequency,start_time,end_time))
    output = joinpath(dirname(path),"pitch_phase_analysis")
    mkpath(output)
    for extension in ("png","pdf")
        target = joinpath(output,"pitch_amplitude_phase." * extension)
        savefig(figure,target)
        println("Figure written to $target")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    default_path = joinpath(@__DIR__,"output","multi_prop","single_propeller_history.txt")
    path = length(ARGS) >= 1 ? ARGS[1] : default_path
    start_time = length(ARGS) >= 2 ? parse(Float64,ARGS[2]) : 3.0
    end_time = length(ARGS) >= 3 ? parse(Float64,ARGS[3]) : 4.9
    plot_pitch_amplitude_phase(path,start_time,end_time)
end
