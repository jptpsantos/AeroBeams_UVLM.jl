# Moving-window pitch amplitude and phase comparison from a recorded history.
# From the repository root:
# julia --project=. lib/WingPropellerUVLM/examples/chang_linear_aeroelastic/plot_pitch_amplitude_phase_history.jl [history.txt or .csv]
module PitchAmplitudePhaseHistory

include(joinpath(@__DIR__,"analyze_pitch_phase.jl"))

function read_pitch_pair(path; wing_location=:tip,propeller_pitch=:auto)
    lines = readlines(path)
    headerIndex = findfirst(line -> startswith(strip(line),"time_s\t") ||
                                   startswith(strip(line),"time_s,"),lines)
    isnothing(headerIndex) && error("No time_s header found in $path")
    delimiter = occursin('\t',lines[headerIndex]) ? '\t' : ','
    names = split(strip(lines[headerIndex]),delimiter)
    data = readdlm(IOBuffer(join(lines[headerIndex+1:end],"\n")),delimiter)
    function column(name)
        index = findfirst(==(name),names)
        isnothing(index) && error("Missing column $name")
        return [value isa Number ? Float64(value) : parse(Float64,value) for value in data[:,index]]
    end
    attachmentName = "propeller_1_attachment_wing_twist_deg"
    hasAttachment = attachmentName in names
    wing_location in (:tip,:attachment) || error("wing_location must be :tip or :attachment")
    propeller_pitch in (:auto,:relative,:combined_hub) || error("Invalid propeller_pitch option")
    useHub = propeller_pitch == :combined_hub || (propeller_pitch == :auto && hasAttachment)
    (wing_location == :attachment || useHub) && !hasAttachment &&
        error("This history does not contain attachment-node twist")
    time = column("time_s")
    wing = column(wing_location == :tip ? "tip_twist_deg" : attachmentName)
    relativePitch = column("propeller_1_pitch_deg")
    propeller = copy(relativePitch)
    if useHub
        propeller .+= column(attachmentName)
    end
    all(isfinite,time) && all(isfinite,wing) && all(isfinite,propeller) || error("Non-finite pitch data")
    all(diff(time) .> 0) || error("Time must be strictly increasing")
    wingLabel = wing_location == :tip ? "Wing-tip pitch" : "Attachment wing pitch"
    propellerLabel = useHub ? "Combined hub pitch" : "Hub pitch relative to wing"
    definition = useHub ?
        "Propeller pitch = attachment wing twist + relative pylon pitch (combined intermediate pitch coordinate)." :
        "Propeller pitch is the recorded hub/pylon pitch rotation relative to the local wing."
    yaw = "propeller_1_yaw_deg" in names ? column("propeller_1_yaw_deg") : nothing
    return (; time,wing,propeller,relativePitch,yaw,wingLabel,propellerLabel,definition,useHub)
end

function pitch_window_history(pair; fit_start=3.0,fit_end=4.9,
    response_start=0.25,cycles_per_window=3.0,center_step=0.025)
    time, wing, propeller = pair.time,pair.wing,pair.propeller
    first(time) <= fit_start < fit_end <= last(time) || error("Late-time fit window lies outside history")
    selected = findall(t -> fit_start <= t <= fit_end,time)
    length(selected) >= 100 || error("Too few samples in late-time window")
    reduced = selected[1:max(1,cld(length(selected),2500)):end]
    referenceTime = (fit_start+fit_end)/2
    mode = fit_mode(time[reduced] .- referenceTime,wing[reduced],propeller[reduced])
    windowLength = cycles_per_window/mode.frequency
    firstCenter = max(first(time),response_start) + windowLength/2
    lastCenter = last(time) - windowLength/2
    firstCenter < lastCenter || error("History is shorter than the moving window")
    centers = collect(firstCenter:center_step:lastCenter)
    amplitudes = zeros(length(centers),2)
    phases = zeros(length(centers))
    residuals = zeros(length(centers),2)
    for (index,center) in enumerate(centers)
        firstIndex = searchsortedfirst(time,center-windowLength/2)
        lastIndex = searchsortedlast(time,center+windowLength/2)
        values = hcat(wing[firstIndex:lastIndex],propeller[firstIndex:lastIndex])
        _, coefficients, prediction = fit_at(time[firstIndex:lastIndex] .- center,
            values,mode.frequency,mode.growth)
        for signal = 1:2
            amplitudes[index,signal] = hypot(coefficients[2,signal],coefficients[3,signal])
            residuals[index,signal] = norm(values[:,signal]-prediction[:,signal]) /
                                    norm(values[:,signal] .- mean(values[:,signal]))
        end
        phiWing = atan(-coefficients[3,1],coefficients[2,1])
        phiPropeller = atan(-coefficients[3,2],coefficients[2,2])
        phases[index] = mod(rad2deg(phiPropeller-phiWing)+180,360)-180
    end
    ratios = amplitudes[:,2] ./ amplitudes[:,1]
    return (; centers,amplitudes,phases,ratios,residuals,windowLength,referenceTime,
        frequency=mode.frequency,growth=mode.growth)
end

function history_panel(time,values; ylabel,ylims,title="",xlabel="")
    plot(time,values; color=:black,linewidth=1.5,label="",xlabel,ylabel,ylims,title,
        xlims=(0,5),xticks=0:1:5,grid=false,framestyle=:axes,tick_direction=:in,
        fontfamily="Times Roman",guidefontsize=12,tickfontsize=10,titlefontsize=11,
        left_margin=5Plots.mm,right_margin=4Plots.mm,top_margin=2Plots.mm,
        bottom_margin=2Plots.mm)
end

function write_reviewer_interpretation(path,pair,result,outputDirectory,fit_start,fit_end)
    selected = findall(t -> fit_start <= t <= fit_end,pair.time)
    centeredTime = pair.time[selected] .- result.referenceTime
    _, coefficients, prediction = fit_at(centeredTime,
        hcat(pair.wing[selected],pair.propeller[selected]),result.frequency,result.growth)
    amplitudes = [hypot(coefficients[2,i],coefficients[3,i]) for i in 1:2]
    phases = [atan(-coefficients[3,i],coefficients[2,i]) for i in 1:2]
    phaseDifference = mod(rad2deg(phases[2]-phases[1])+180,360)-180
    amplitudeRatio = amplitudes[2]/amplitudes[1]
    dampingRatio = -result.growth/hypot(result.growth,2π*result.frequency)
    growthFactor = exp(result.growth*(fit_end-fit_start))
    yawDifference = nothing
    yawRatio = nothing
    if !isnothing(pair.yaw)
        _, whirlCoefficients, _ = fit_at(centeredTime,
            hcat(pair.relativePitch[selected],pair.yaw[selected]),result.frequency,result.growth)
        pitchPhase = atan(-whirlCoefficients[3,1],whirlCoefficients[2,1])
        yawPhase = atan(-whirlCoefficients[3,2],whirlCoefficients[2,2])
        yawDifference = mod(rad2deg(yawPhase-pitchPhase)+180,360)-180
        yawRatio = hypot(whirlCoefficients[2,2],whirlCoefficients[3,2]) /
                   hypot(whirlCoefficients[2,1],whirlCoefficients[3,1])
    end
    reportPath = joinpath(outputDirectory,"reviewer_pitch_mode_interpretation.md")
    open(reportPath,"w") do io
        println(io,"# Pitch relationship and unstable-mode interpretation\n")
        println(io,"Source: $(abspath(path))\n")
        println(io,"$(pair.wingLabel) is compared with $(lowercase(pair.propellerLabel)). $(pair.definition)\n")
        println(io,"## Quantitative result\n")
        @printf(io,"The %.2f–%.2f s history is fitted by q(t) = q_mean + A exp[σ(t−t_ref)] cos[2πf(t−t_ref)+φ]. Amplitudes are single-sided and refer to t_ref = %.2f s.\n\n",fit_start,fit_end,result.referenceTime)
        @printf(io,"| Quantity | Result |\n|---|---:|\n| Frequency | %.5f Hz |\n| Fitted growth rate σ | %+.5f s⁻¹ |\n| Equivalent fitted damping ratio | %+.4f%% |\n| Wing pitch amplitude | %.6f° |\n| Propeller pitch amplitude | %.6f° |\n| Propeller/wing amplitude ratio | %.5f |\n| Propeller phase minus wing phase | %+.3f° |\n| Propeller time lead | %.3f ms |\n| Envelope factor across fit window | %.5f |\n\n",
            result.frequency,result.growth,100dampingRatio,amplitudes[1],amplitudes[2],
            amplitudeRatio,phaseDifference,1000phaseDifference/(360result.frequency),growthFactor)
        if !isnothing(yawDifference)
            @printf(io,"The relative pylon yaw/pitch amplitude ratio is %.5f and yaw leads pitch by %+.3f° in the saved structural coordinates. Thus the relative pylon response is very close to equal-amplitude quadrature, producing a nearly circular whirl orbit.\n\n",yawRatio,yawDifference)
        end
        println(io,"## Mode interpretation supported by the history\n")
        println(io,"The growing wing and propeller pitch oscillations maintain nearly constant amplitude ratio and relative phase after the early transient. This is evidence that the late response is governed by one dominant coupled aeroelastic mode. The small pitch phase difference means that these two recorded pitch coordinates are nearly in phase. The pylon pitch–yaw quadrature identifies the whirl component, while the wing deforms at the same frequency. This combination characterizes a wing–pylon whirl mode.\n")
        println(io,"The positive fitted growth rate describes weak growth over this window. It is a response-based modal estimate, not an independently established flutter boundary. The motion alone does not quantify separate wing/propeller aerodynamic work or prove that cross-aerodynamic interaction destabilizes the mode.\n")
        println(io,"## Suggested manuscript paragraph\n")
        if !pair.useHub
            println(io,"> Here wing pitch denotes wing-tip twist, while propeller pitch denotes the hub/pylon pitch rotation relative to the attachment wing, consistent with the recorded model coordinate.\n")
        end
        @printf(io,"> To clarify the mode associated with the instability, the wing-pitch and propeller-pitch histories were analysed using a common exponentially varying sinusoid. The dominant response has a frequency of %.3f Hz and a fitted growth rate of %+.3f s⁻¹. At t = %.2f s, the wing-pitch and propeller-pitch amplitudes are %.3f° and %.3f°, respectively. Their amplitude ratio is %.2f, and the propeller-pitch coordinate leads wing pitch by %.2f°. The nearly constant amplitude ratio and phase difference indicate a dominant coupled mode after the initial transient.",
            result.frequency,result.growth,result.referenceTime,amplitudes[1],amplitudes[2],amplitudeRatio,phaseDifference)
        if !isnothing(yawDifference)
            @printf(io," In addition, relative pylon pitch and yaw have an amplitude ratio of %.3f and a phase separation of %.2f°, indicating a nearly circular whirl component coupled to wing deformation.",yawRatio,abs(yawDifference))
        end
        println(io,"\n")
        println(io,"## Suggested response to the reviewer\n")
        @printf(io,"> Thank you for this suggestion. We have added an analysis of the time history using wing-tip pitch as the wing reference. The new figure shows the wing and propeller pitch histories with their amplitude envelopes, their amplitude ratio, and the phase difference. The late-time propeller-to-wing pitch amplitude ratio is %.2f and propeller pitch leads wing-tip pitch by %.2f°. Both motions occur at %.3f Hz. We also discuss the near-quadrature pitch–yaw relationship, which identifies the whirl component of the coupled wing–pylon mode. The coordinate definitions and phase convention have been stated explicitly.\n\n",amplitudeRatio,phaseDifference,result.frequency)
        if !pair.useHub
            println(io,"Coordinate statement for the caption and response: propeller pitch is the hub/pylon pitch rotation relative to the attachment wing, as stored in the original response file; wing pitch is wing-tip twist. The reported phase difference is between these defined model coordinates. No absolute hub orientation has been inferred.\n")
        end
        println(io,"Suggested caption: Wing–propeller pitch relationship in the aeroelastic response: (a) wing pitch and (b) propeller pitch, with dashed local amplitude envelopes; (c) propeller-to-wing pitch amplitude ratio; and (d) phase difference, defined as propeller phase minus wing phase. Positive phase denotes propeller lead. Envelopes and phase are estimated in three-period moving windows after the initial forcing interval.\n")
        println(io,"## Figures\n")
        println(io,"- [Pitch histories, amplitude ratio and phase difference](reviewer_pitch_mode.pdf)\n- [Whirl orbit and modal phase relationships](whirl_mode.pdf)\n")
        println(io,"The response text is ready for insertion alongside these figures. No manuscript document has been edited. The analysis characterizes the growing mode; it does not determine individual aerodynamic work contributions or a forward/backward whirl classification.")
    end
    write(joinpath(outputDirectory,"whirl_interpretation.md"),read(reportPath,String))
    return reportPath,coefficients
end

function plot_pitch_amplitude_phase_history(path; output_directory=nothing,
    wing_location=:tip,propeller_pitch=:auto,fit_start=3.0,fit_end=4.9,
    response_start=0.25,cycles_per_window=3.0)
    pair = read_pitch_pair(path; wing_location,propeller_pitch)
    result = pitch_window_history(pair; fit_start,fit_end,response_start,cycles_per_window)
    outputDirectory = isnothing(output_directory) ? joinpath(dirname(path),"pitch_phase_analysis") : output_directory
    mkpath(outputDirectory)
    # Both quantities share the same time windows; positive phase means that
    # propeller pitch leads wing pitch under the cosine phase convention.
    ratioLimits = (floor(minimum(result.ratios)-0.5),ceil(maximum(result.ratios)+0.5))
    phaseLimits = (floor(minimum(result.phases)-0.5),ceil(maximum(result.phases)+0.5))
    ratioPanel = history_panel(result.centers,result.ratios;
        ylabel="Amplitude ratio (A_p / A_w)",ylims=ratioLimits,title="(a)")
    phasePanel = history_panel(result.centers,result.phases;
        ylabel="Phase difference (deg)",ylims=phaseLimits,title="(b)",xlabel="Time (s)")
    figure = plot(ratioPanel,phasePanel;layout=(2,1),size=(850,620),dpi=300,
        background_color=:white,foreground_color=:black)
    wingPanel = history_panel(result.centers,result.amplitudes[:,1];
        ylabel="Wing amplitude (deg)",ylims=(0,1.08maximum(result.amplitudes[:,1])),
        title="(a) "*pair.wingLabel)
    propellerPanel = history_panel(result.centers,result.amplitudes[:,2];
        ylabel="Propeller amplitude (deg)",ylims=(0,1.08maximum(result.amplitudes[:,2])),
        title="(b) "*pair.propellerLabel,xlabel="Time (s)")
    amplitudesFigure = plot(wingPanel,propellerPanel;layout=(2,1),size=(850,620),dpi=300,
        background_color=:white,foreground_color=:black)
    combinedFigure = plot(
        history_panel(result.centers,result.amplitudes[:,1];ylabel="Wing amplitude (deg)",
            ylims=(0,1.08maximum(result.amplitudes[:,1])),title="(a) "*pair.wingLabel),
        history_panel(result.centers,result.amplitudes[:,2];ylabel="Propeller amplitude (deg)",
            ylims=(0,1.08maximum(result.amplitudes[:,2])),title="(b) "*pair.propellerLabel),
        history_panel(result.centers,result.ratios;ylabel="Amplitude ratio (A_p / A_w)",
            ylims=ratioLimits,title="(c)"),
        history_panel(result.centers,result.phases;ylabel="Phase difference (deg)",
            ylims=phaseLimits,title="(d)",xlabel="Time (s)");
        layout=(4,1),size=(850,1080),dpi=300,background_color=:white,foreground_color=:black)
    interpretationPath,lateCoefficients = write_reviewer_interpretation(
        path,pair,result,outputDirectory,fit_start,fit_end)
    pitchWingPanel = history_panel(pair.time,pair.wing;
        ylabel="Wing pitch (deg)",ylims=(-1.15maximum(abs,pair.wing),1.15maximum(abs,pair.wing)),
        title="(a) "*pair.wingLabel)
    pitchPropellerPanel = history_panel(pair.time,pair.propeller;
        ylabel="Propeller pitch (deg)",ylims=(-1.15maximum(abs,pair.propeller),1.15maximum(abs,pair.propeller)),
        title="(b) "*pair.propellerLabel)
    for (panel,signal) in ((pitchWingPanel,1),(pitchPropellerPanel,2))
        for sign in (-1,1)
            plot!(panel,result.centers,lateCoefficients[1,signal] .+ sign .* result.amplitudes[:,signal];
                color=:gray,linestyle=:dash,linewidth=1.3,label="")
        end
    end
    reviewerFigure = plot(pitchWingPanel,pitchPropellerPanel,
        history_panel(result.centers,result.ratios;ylabel="Amplitude ratio (A_p / A_w)",
            ylims=ratioLimits,title="(c)"),
        history_panel(result.centers,result.phases;ylabel="Phase difference (deg)",
            ylims=phaseLimits,title="(d)",xlabel="Time (s)");
        layout=(4,1),size=(850,1080),dpi=300,background_color=:white,foreground_color=:black)
    for extension in ("pdf","png")
        savefig(figure,joinpath(outputDirectory,"pitch_amplitude_ratio_phase_history."*extension))
        savefig(amplitudesFigure,joinpath(outputDirectory,"pitch_amplitude_histories."*extension))
        savefig(combinedFigure,joinpath(outputDirectory,"pitch_amplitude_phase_history."*extension))
        savefig(reviewerFigure,joinpath(outputDirectory,"reviewer_pitch_mode."*extension))
    end
    open(joinpath(outputDirectory,"pitch_amplitude_phase_history.csv"),"w") do io
        println(io,"window_center_s,wing_amplitude_deg,propeller_amplitude_deg,propeller_to_wing_amplitude_ratio,propeller_minus_wing_phase_deg,wing_normalized_fit_residual,propeller_normalized_fit_residual")
        for index in eachindex(result.centers)
            println(io,join((result.centers[index],result.amplitudes[index,1],result.amplitudes[index,2],
                result.ratios[index],result.phases[index],result.residuals[index,1],result.residuals[index,2]),','))
        end
    end
    late = result.centers .>= fit_start
    reportPath = joinpath(outputDirectory,"pitch_amplitude_phase_history_caption.md")
    open(reportPath,"w") do io
        println(io,"Source: $(abspath(path))\n")
        println(io,"Wing coordinate: $(pair.wingLabel). $(pair.definition)\n")
        @printf(io,"Common fitted frequency: %.5f Hz. Growth rate: %.5f 1/s.\n\n",result.frequency,result.growth)
        @printf(io,"Moving-window duration: %.5f s (%.1f periods). Window centers span %.3f–%.3f s; only complete windows after %.3f s are used.\n\n",
            result.windowLength,cycles_per_window,first(result.centers),last(result.centers),response_start)
        @printf(io,"For centers at t >= %.2f s: mean amplitude ratio %.5f; mean phase difference %+.4f deg; ranges %.5f–%.5f and %+.4f–%+.4f deg.\n\n",
            fit_start,mean(result.ratios[late]),mean(result.phases[late]),extrema(result.ratios[late])...,extrema(result.phases[late])...)
        println(io,"Caption (four-panel figure): Time-dependent comparison of wing pitch and propeller pitch extracted from the response history: (a) wing-pitch amplitude, (b) propeller-pitch amplitude, (c) pitch amplitude ratio A_p/A_w and (d) phase difference Δφ = φ_p − φ_w. Amplitude is the single-sided oscillation envelope, not peak-to-peak amplitude. A positive phase difference means that propeller pitch leads wing pitch. Each estimate uses a three-period local sinusoidal fit at the common late-time frequency with an exponential envelope. The initial forcing interval and incomplete windows at the endpoints are excluded.\n")
        println(io,"Amplitude and phase are derived from the recorded signals; curves between window centers are linearly joined. The local fits also estimate a mean offset. The CSV includes normalized fit residuals so mixed-mode transients can be assessed.\n")
        if !pair.useHub
            println(io,"For this saved file, propeller pitch denotes the hub/pylon rotation relative to the local wing, and wing pitch denotes wing-tip twist. The plotted phase difference is between those model coordinates.\n")
        end
        println(io,"No airspeed label is assigned because the selected history does not establish its operating-condition metadata.")
    end
    println(read(reportPath,String))
    println("Figures and numerical data written to $(abspath(outputDirectory))")
    println("Reviewer interpretation written to $(abspath(interpretationPath))")
    return (; pair,result,outputDirectory)
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    path = isempty(ARGS) ? joinpath(@__DIR__,"output","multi_prop","single_propeller_history.txt") : ARGS[1]
    output = length(ARGS) >= 2 ? ARGS[2] : nothing
    PitchAmplitudePhaseHistory.plot_pitch_amplitude_phase_history(path; output_directory=output)
end
