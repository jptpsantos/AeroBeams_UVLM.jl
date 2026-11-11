# Analyze the wing-tip twist and relative pylon pitch in a saved Chang history.
# Usage: julia --project=lib/WingPropellerUVLM analyze_pitch_phase.jl [history.txt] [start_s] [end_s]
# The common exponentially growing/decaying sinusoid is estimated by variable
# projection. Amplitudes refer to the middle of the chosen window.
using DelimitedFiles, LinearAlgebra, Statistics, Printf
get!(ENV,"GKSwstype","100")
using Plots
using ColorSchemes

function read_pitch_history(path)
    lines = readlines(path)
    index = findfirst(line -> startswith(strip(line),"time_s\t") || startswith(strip(line),"time_s,"),lines)
    isnothing(index) && error("Could not find the time_s header in $path")
    delimiter = occursin('\t',lines[index]) ? '\t' : ','
    columns = split(strip(lines[index]),delimiter)
    data = readdlm(IOBuffer(join(lines[index+1:end],"\n")),delimiter)
    column(name) = begin
        columnIndex = findfirst(==(name),columns)
        isnothing(columnIndex) && error("Missing column $name in $path")
        return Float64.(data[:,columnIndex])
    end
    t, wing, pylon = column("time_s"),column("tip_twist_deg"),column("propeller_1_pitch_deg")
    all(isfinite,t) && all(isfinite,wing) && all(isfinite,pylon) || error("History contains nonfinite pitch values")
    all(diff(t) .> 0) || error("Time must be strictly increasing")
    attachmentName = "propeller_1_attachment_wing_twist_deg"
    attachment = attachmentName in columns ? column(attachmentName) : nothing
    yaw = "propeller_1_yaw_deg" in columns ? column("propeller_1_yaw_deg") : nothing
    return (; time=t,wing,pylon,attachment,yaw)
end

function read_history(path)
    history = read_pitch_history(path)
    return history.time,history.wing,history.pylon
end

function fit_at(t, signals, frequency, growth)
    envelope = exp.(growth .* t)
    angle = 2π * frequency .* t
    design = hcat(ones(length(t)),envelope .* cos.(angle),envelope .* sin.(angle))
    coefficients = design \ signals
    prediction = design * coefficients
    score = sum(abs2,signals - prediction) / sum(abs2,signals .- mean(signals,dims=1))
    return score, coefficients, prediction
end

function fit_mode(t, wing, pylon)
    # Normalize the two signals for equal influence on the shared parameters.
    values = hcat(wing,pylon)
    scales = vec(std(values,dims=1))
    all(scales .> 1e-12) || error("A pitch signal has negligible variance")
    scaled = values ./ scales'
    best = (Inf,0.0,0.0)
    for frequency in 1.0:0.1:20.0, growth in -2.0:0.5:2.0
        score, = fit_at(t,scaled,frequency,growth)
        if score < best[1]
            best = (score,frequency,growth)
        end
    end
    # Refine both nonlinear parameters on successively finer local grids.
    for (df,dg) in ((0.02,0.1),(0.004,0.02),(0.0008,0.004))
        _, frequency0, growth0 = best
        for frequency in max(0.01,frequency0-5df):df:frequency0+5df,
            growth in growth0-5dg:dg:growth0+5dg
            score, = fit_at(t,scaled,frequency,growth)
            if score < best[1]
                best = (score,frequency,growth)
            end
        end
    end
    _, frequency, growth = best
    _, coefficients, prediction = fit_at(t,values,frequency,growth)
    residuals = [norm(values[:,i] - prediction[:,i]) / norm(values[:,i] .- mean(values[:,i])) for i in 1:2]
    amplitude = [hypot(coefficients[2,i],coefficients[3,i]) for i in 1:2]
    phase = [atan(-coefficients[3,i],coefficients[2,i]) for i in 1:2]
    delta = mod(rad2deg(phase[2]-phase[1])+180,360)-180
    return (; frequency,growth,coefficients,prediction,residuals,amplitude,delta)
end

function analyze(path, start_time, end_time)
    t, wing, pylon = read_history(path)
    first(t) <= start_time < end_time <= last(t) || error("Window must lie within $(first(t))–$(last(t)) s")
    selected = findall(x -> start_time <= x <= end_time,t)
    length(selected) >= 100 || error("Analysis window has fewer than 100 samples")
    # Limit samples for the search; the full data are used for the final fit.
    stride = max(1,cld(length(selected),2500))
    reduced = selected[1:stride:end]
    midpoint = (start_time+end_time)/2
    mode = fit_mode(t[reduced] .- midpoint,wing[reduced],pylon[reduced])
    values = hcat(wing[selected],pylon[selected])
    _, coefficients, prediction = fit_at(t[selected] .- midpoint,values,mode.frequency,mode.growth)
    amplitudes = [hypot(coefficients[2,i],coefficients[3,i]) for i in 1:2]
    delta = mod(rad2deg(atan(-coefficients[3,2],coefficients[2,2])-atan(-coefficients[3,1],coefficients[2,1]))+180,360)-180
    rms = [norm(values[:,i]-prediction[:,i]) / norm(values[:,i] .- mean(values[:,i])) for i in 1:2]
    output = joinpath(dirname(path),"pitch_phase_analysis")
    mkpath(output)
    report = joinpath(output,"pitch_phase_summary.txt")
    open(report,"w") do io
        println(io,"Source: $(abspath(path))")
        @printf(io,"Window: %.6f to %.6f s; amplitude reference time: %.6f s\n",start_time,end_time,midpoint)
        println(io,"Fit: offset + exp(growth*(t-reference)) * (a*cos(2*pi*f*(t-reference)) + b*sin(2*pi*f*(t-reference)))")
        @printf(io,"Frequency: %.5f Hz\nGrowth rate: %.5f 1/s\n",mode.frequency,mode.growth)
        @printf(io,"Wing-tip twist amplitude: %.6f deg\nRelative pylon pitch amplitude: %.6f deg\n",amplitudes...)
        @printf(io,"Relative pylon / wing-tip amplitude: %.5f\n",amplitudes[2]/amplitudes[1])
        @printf(io,"Relative pylon phase minus wing-tip phase: %+.3f deg\n",delta)
        @printf(io,"Normalized RMS fit residual, wing tip: %.5f; relative pylon: %.5f\n",rms...)
        println(io,"Phase convention: cos(2*pi*f*(t-reference)+phase); positive difference means relative pylon leads wing tip.")
        println(io,"The exported propeller pitch is pylon pitch relative to the local wing, not absolute hub pitch.")
        println(io,"The local wing twist at eta=0.83 is absent, so absolute hub pitch cannot be reconstructed from this file.")
        println(io,"Amplitudes and phase describe the fitted oscillation, not peak-to-peak values; check residuals and window sensitivity.")
        println(io,"The frequency search initially covers 1–20 Hz and growth rates -2 to 2 1/s, then refines locally.")
    end
    println(read(report,String))
    wing_plot = plot(t[selected],wing[selected]; label="Wing-tip twist", xlabel="Time (s)", ylabel="Angle (deg)", linewidth=1.3)
    plot!(wing_plot,t[selected],prediction[:,1]; label="Single-mode fit", linewidth=1.5, linestyle=:dash)
    pylon_plot = plot(t[selected],pylon[selected]; label="Relative pylon pitch", xlabel="Time (s)", ylabel="Angle (deg)", linewidth=1.3)
    plot!(pylon_plot,t[selected],prediction[:,2]; label="Single-mode fit", linewidth=1.5, linestyle=:dash)
    savefig(plot(wing_plot,pylon_plot; layout=(2,1), size=(1000,700),margin=5Plots.mm),joinpath(output,"pitch_phase_fit.png"))
    return report
end

"""
    identify_pitch_mode(time, wing, propeller; fit_start=3.0, fit_end=4.9, ...)

Fit two synchronized rotations (degrees) to one exponentially varying sinusoid.
The reported convention is R = A_wing/A_propeller and Δφ = φ_wing−φ_propeller.
Positive Δφ means wing lead. A failed quality check returns phase=NaN.
"""
function identify_pitch_mode(time,wing,propeller;
    fit_start=3.0,fit_end=4.9,minimum_amplitude_deg=1e-4,
    minimum_amplitude_ratio=1e-4,maximum_residual=0.10,
    maximum_split_phase_deg=5.0,maximum_split_ratio_change=0.05)
    length(time) == length(wing) == length(propeller) || error("Time and pitch vectors must have identical lengths")
    all(isfinite,time) && all(isfinite,wing) && all(isfinite,propeller) || error("Nonfinite pitch data")
    all(diff(time) .> 0) || error("Time must be strictly increasing")
    first(time) <= fit_start < fit_end <= last(time) || error("Fit interval lies outside the history")
    indices = findall(t -> fit_start <= t <= fit_end,time)
    length(indices) >= 100 || error("At least 100 synchronized samples are required")
    referenceTime = (fit_start+fit_end)/2
    fitTime = time[indices] .- referenceTime
    values = hcat(wing[indices],propeller[indices])
    scales = vec(std(values,dims=1))
    all(scales .> 1e-12) || error("A pitch signal has negligible variance; its phase is unidentifiable")
    reduced = 1:max(1,cld(length(indices),2500)):length(indices)
    mode = fit_mode(fitTime[reduced],values[reduced,1],values[reduced,2])
    frequency, growth = mode.frequency,mode.growth
    # Refine the existing variable-projection search without an optimizer dependency.
    scaled = values[reduced,:] ./ scales'
    bestScore, = fit_at(fitTime[reduced],scaled,frequency,growth)
    for (df,dg) in ((0.00016,0.0008),(0.000032,0.00016),(0.0000064,0.000032))
        frequency0, growth0 = frequency,growth
        for trialFrequency in frequency0-5df:df:frequency0+5df,
            trialGrowth in growth0-5dg:dg:growth0+5dg
            score, = fit_at(fitTime[reduced],scaled,trialFrequency,trialGrowth)
            if score < bestScore
                bestScore, frequency, growth = score,trialFrequency,trialGrowth
            end
        end
    end
    _, coefficients, prediction = fit_at(fitTime,values,frequency,growth)
    amplitude = [hypot(coefficients[2,i],coefficients[3,i]) for i in 1:2]
    phases = [atan(-coefficients[3,i],coefficients[2,i]) for i in 1:2]
    wrap_phase(x) = mod(rad2deg(x)+180,360)-180
    phaseDifference = wrap_phase(phases[1]-phases[2])
    ratio = amplitude[1]/amplitude[2]
    residuals = [norm(values[:,i]-prediction[:,i])/norm(values[:,i] .- mean(values[:,i])) for i in 1:2]
    # Independently fit the frequency of each channel. Raw local-maxima counts
    # are unreliable when tiny high-frequency ripples create extra wing peaks.
    independentFrequencies = fill(NaN,2)
    for i in 1:2
        bestIndependentScore = Inf
        for trialFrequency in 1.0:0.02:20.0
            score, = fit_at(fitTime[reduced],values[reduced,i:i],trialFrequency,growth)
            if score < bestIndependentScore
                bestIndependentScore = score
                independentFrequencies[i] = trialFrequency
            end
        end
    end
    splitPhases, splitRatios = zeros(2),zeros(2)
    middle = div(length(fitTime),2)
    for (i,part) in enumerate((1:middle,middle+1:length(fitTime)))
        _, localCoefficients, _ = fit_at(fitTime[part],values[part,:],frequency,growth)
        splitPhases[i] = wrap_phase(atan(-localCoefficients[3,1],localCoefficients[2,1]) -
                                          atan(-localCoefficients[3,2],localCoefficients[2,2]))
        splitRatios[i] = hypot(localCoefficients[2,1],localCoefficients[3,1]) /
                         hypot(localCoefficients[2,2],localCoefficients[3,2])
    end
    splitPhaseDifference = abs(mod(splitPhases[1]-splitPhases[2]+180,360)-180)
    splitRatioChange = abs(splitRatios[1]-splitRatios[2])/ratio
    reasons = String[]
    minimum(amplitude) >= minimum_amplitude_deg || push!(reasons,"pitch amplitude below $(minimum_amplitude_deg) deg")
    minimum(amplitude)/maximum(amplitude) >= minimum_amplitude_ratio || push!(reasons,"one pitch is too weak relative to the other")
    maximum(residuals) <= maximum_residual || push!(reasons,"single-mode fit residual exceeds $(maximum_residual)")
    (fit_end-fit_start)*frequency >= 5 || push!(reasons,"fewer than five oscillation periods")
    frequencyTolerance = max(0.02frequency,0.1/(fit_end-fit_start))
    all(isfinite,independentFrequencies) && maximum(abs.(independentFrequencies .- frequency)) <= frequencyTolerance ||
        push!(reasons,"independent signal frequencies disagree or are unavailable")
    splitPhaseDifference <= maximum_split_phase_deg || push!(reasons,"phase is sensitive to the fit interval")
    splitRatioChange <= maximum_split_ratio_change || push!(reasons,"amplitude ratio is sensitive to the fit interval")
    maximum(diff(time[indices]))*frequency <= 0.1 || push!(reasons,"fewer than ten samples per period")
    reliable = isempty(reasons)
    return (; reliable,reasons,frequency_hz=frequency,growth_rate_per_s=growth,
        damping_ratio_percent=-100growth/hypot(growth,2π*frequency),
        reference_time_s=referenceTime,fit_start_s=fit_start,fit_end_s=fit_end,
        wing_amplitude_deg=amplitude[1],propeller_amplitude_deg=amplitude[2],
        wing_to_propeller_amplitude_ratio=reliable ? ratio : NaN,
        wing_minus_propeller_phase_deg=reliable ? phaseDifference : NaN,
        normalized_residuals=residuals,independent_frequencies_hz=independentFrequencies,
        split_phase_difference_deg=splitPhaseDifference,split_ratio_relative_change=splitRatioChange,
        time=time[indices],signals=values,coefficients,prediction)
end

function run_pitch_reviewer_analysis(; output_directory=joinpath(@__DIR__,"output","multi_prop","pitch_phase_analysis","reviewer_response"))
    cases = [
        # Airspeed supplied by the author; the original TXT has no speed metadata.
        (label="V = 83 m/s",path=joinpath(@__DIR__,"output","multi_prop","single_propeller_history.txt"),summary=nothing)
    ]
    histories = [read_pitch_history(case.path) for case in cases]
    modes = [identify_pitch_mode(history.time,history.wing,history.pylon) for history in histories]
    mkpath(output_directory)
    panels = Plots.Plot[]
    for (index,(case,mode)) in enumerate(zip(cases,modes))
        if !isnothing(case.summary)
            summaryText = read(case.summary,String)
            endMatch = match(r"(?m)^integrated_end_time_s = (.+)",summaryText)
            stepMatch = match(r"(?m)^integrated_steps = (.+)",summaryText)
            !isnothing(endMatch) && isapprox(parse(Float64,endMatch[1]),last(histories[index].time);atol=1e-9) || error("Summary/history end times do not match")
            !isnothing(stepMatch) && parse(Int,stepMatch[1])+1 == length(histories[index].time) || error("Summary/history sample counts do not match")
        end
        if mode.reliable
            # Same reference interval for the phase and amplitude plots.
            visible = findall(t -> 3.65 <= t <= 4.25,mode.time)
            time = mode.time[visible]
            rotations = (mode.signals[visible,:] .- mode.coefficients[1,:]') ./
                exp.(mode.growth_rate_per_s .* (time .- mode.reference_time_s))
            wingColor = get(ColorSchemes.colorschemes[:darkrainbow],1.0)
            for normalized in (true,false)
                plotTime = normalized ? time : histories[index].time
                values = normalized ? rotations ./ [mode.wing_amplitude_deg mode.propeller_amplitude_deg] :
                                      hcat(histories[index].wing,histories[index].pylon)
                panel = plot(plotTime,values[:,1];label="Wing-tip twist",color=wingColor,
                    linestyle=:dash,linewidth=2.5,xlabel="Time (s)",
                    ylabel=normalized ? "Normalized angle" : "Angle (deg)",
                    xlims=normalized ? (3.65,4.25) : (0.0,5.0),
                    xticks=normalized ? collect(range(3.65,4.25;length=7)) : collect(0.0:1.0:5.0),
                    ylims=normalized ? (-1.5,1.5) : (-7.5,7.5),
                    yticks=normalized ? collect(-1.5:0.5:1.5) : collect(-7.5:2.5:7.5),
                    grid=false,title="",legend=:top,legend_column=2,
                    legend_foreground_color=:transparent,legend_background_color=:transparent,
                    fontfamily="DejaVu Sans",guidefontsize=18,tickfontsize=16,legendfontsize=16,
                    framestyle=:axes,tick_direction=:in,margin=6Plots.mm,size=(1200,540),dpi=300)
                plot!(panel,plotTime,values[:,2];label="Pylon pitch",color=:black,linewidth=2.5)
                if normalized
                    annotate!(panel,3.95,-1.25,
                        text(@sprintf("Δφ = %+.2f°",mode.wing_minus_propeller_phase_deg),"DejaVu Sans",16,:center))
                end
                push!(panels,panel)
            end
            println("Wing curve color, darkrainbow red: ",wingColor)
        else
            push!(panels,plot(;title="",annotations=(0.5,0.5,"Unreliable phase: "*join(mode.reasons,"; ")),
                xlims=(0,1),ylims=(0,1),axis=false))
            push!(panels,plot(;axis=false,legend=false))
        end
    end
    for (figure,basename) in zip(panels,("pitch_mode_comparison_83mps","pitch_amplitude_83mps")), extension in ("pdf","png")
        # Write a fresh file before replacing a figure open in a Windows preview.
        temporaryPath = tempname(output_directory)*"."*extension
        savefig(figure,temporaryPath)
        mv(temporaryPath,joinpath(output_directory,basename*"."*extension);force=true)
    end
    open(joinpath(output_directory,"pitch_mode_results.csv"),"w") do io
        println(io,"case,frequency_hz,growth_rate_per_s,wing_amplitude_deg,propeller_amplitude_deg,wing_to_propeller_amplitude_ratio,wing_minus_propeller_phase_deg,wing_nrmse,propeller_nrmse,wing_independent_frequency_hz,propeller_independent_frequency_hz,split_phase_difference_deg,split_ratio_relative_change,reliable,quality_notes")
        for (case,mode) in zip(cases,modes)
            println(io,join((case.label,mode.frequency_hz,mode.growth_rate_per_s,mode.wing_amplitude_deg,
                mode.propeller_amplitude_deg,mode.wing_to_propeller_amplitude_ratio,mode.wing_minus_propeller_phase_deg,
                mode.normalized_residuals...,mode.independent_frequencies_hz...,mode.split_phase_difference_deg,
                mode.split_ratio_relative_change,mode.reliable,join(mode.reasons,"; ")),','))
        end
    end
    for (case,mode) in zip(cases,modes)
        println(case.label,"\n",mode.reliable ? @sprintf("f=%.6f Hz, growth=%+.6f 1/s, R=%.7f, phase(wing-propeller)=%+.4f deg; residuals=%.5f, %.5f",
            mode.frequency_hz,mode.growth_rate_per_s,mode.wing_to_propeller_amplitude_ratio,
            mode.wing_minus_propeller_phase_deg,mode.normalized_residuals...) : join(mode.reasons,"; "))
    end
    return (; cases,histories,modes,output_directory)
end

if abspath(PROGRAM_FILE) == @__FILE__
    if !isempty(ARGS) && ARGS[1] == "--reviewer"
        run_pitch_reviewer_analysis()
    else
        default_path = joinpath(@__DIR__,"output","multi_prop","single_propeller_history.txt")
        path = length(ARGS) >= 1 ? ARGS[1] : default_path
        start_time = length(ARGS) >= 2 ? parse(Float64,ARGS[2]) : 3.0
        end_time = length(ARGS) >= 3 ? parse(Float64,ARGS[3]) : 4.9
        analyze(path,start_time,end_time)
    end
end
