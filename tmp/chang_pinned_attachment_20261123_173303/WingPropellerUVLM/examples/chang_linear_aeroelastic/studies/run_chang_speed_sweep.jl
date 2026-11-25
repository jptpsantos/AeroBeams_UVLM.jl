# Chang airspeed stability sweep.

"""
Run consecutive Chang aeroelastic cases and stop at the first divergent case.

The default sweep is 80:1:84 m/s. Each case runs in a fresh Julia process and
writes its full console output to a per-speed log. A case is divergent when it
does not reach the requested end time, exceeds the hard propeller-angle limit,
or has a fitted pitch/yaw envelope growth rate and ratio above the configured
thresholds.
"""

using DelimitedFiles
using Printf
using Plots

const SWEEP_STUDY_DIR = @__DIR__
const SWEEP_EXAMPLE_DIR = normpath(joinpath(SWEEP_STUDY_DIR,".."))
const SWEEP_PACKAGE_DIR = normpath(joinpath(SWEEP_EXAMPLE_DIR,"..",".."))
const SWEEP_CASE_RUNNER = joinpath(SWEEP_EXAMPLE_DIR,"run_chang_linear_aeroelastic.jl")

sweep_float(name,default) = parse(Float64,get(ENV,name,string(default)))

function speed_values(firstSpeed,lastSpeed,speed_step)
    speed_step > 0 || error("CHANG_SWEEP_SPEED_STEP_MPS must be positive")
    lastSpeed >= firstSpeed || error(
        "CHANG_SWEEP_SPEED_STOP_MPS must not be less than CHANG_SWEEP_SPEED_START_MPS"
    )
    count = floor(Int,(lastSpeed - firstSpeed) / speed_step + 1.0e-10)
    values = [firstSpeed + index * speed_step for index = 0:count]
    if isempty(values) || values[end] < lastSpeed - 1.0e-10
        push!(values,lastSpeed)
    end
    return values
end

function local_absolute_peaks(time,signal; fit_start_s)
    indices = findall(value -> value >= fit_start_s,time)
    length(indices) >= 3 || return Float64[],Float64[]

    baseline = sum(signal[indices]) / length(indices)
    amplitude = abs.(signal .- baseline)
    peak_indices = Int[]
    for index = first(indices)+1:last(indices)-1
        if amplitude[index] > amplitude[index-1] &&
           amplitude[index] >= amplitude[index+1] &&
           amplitude[index] > 100 * eps(Float64)
            push!(peak_indices,index)
        end
    end
    return time[peak_indices],amplitude[peak_indices]
end

function oscillation_metrics(time,signal; fit_start_s,minimum_peaks=4)
    peakTime,peak_amplitude = local_absolute_peaks(time,signal; fit_start_s)
    if length(peakTime) < minimum_peaks
        return (;
            available=false,
            peak_count=length(peakTime),
            growth_rate_per_s=NaN,
            fitted_envelope_ratio=NaN,
            frequency_hz=NaN,
            fit_r_squared=NaN
        )
    end

    log_amplitude = log.(peak_amplitude)
    meanTime = sum(peakTime) / length(peakTime)
    meanLogAmplitude = sum(log_amplitude) / length(log_amplitude)
    centeredTime = peakTime .- meanTime
    denominator = sum(abs2,centeredTime)
    denominator > eps(Float64) || return (;
        available=false,
        peak_count=length(peakTime),
        growth_rate_per_s=NaN,
        fitted_envelope_ratio=NaN,
        frequency_hz=NaN,
        fit_r_squared=NaN
    )

    growth_rate = sum(centeredTime .* (log_amplitude .- meanLogAmplitude)) /
                  denominator
    fitted_log_amplitude = meanLogAmplitude .+ growth_rate .* centeredTime
    residualSumSquares = sum(abs2,log_amplitude .- fitted_log_amplitude)
    totalSumSquares = sum(abs2,log_amplitude .- meanLogAmplitude)
    fit_r_squared = totalSumSquares > eps(Float64) ?
                    1.0 - residualSumSquares / totalSumSquares : 1.0
    envelope_ratio = exp(growth_rate * (peakTime[end] - peakTime[1]))
    # Consecutive maxima of |signal| are normally half a cycle apart.
    meanHalfPeriod = sum(diff(peakTime)) / (length(peakTime) - 1)
    frequency_hz = 1.0 / (2.0 * meanHalfPeriod)

    return (;
        available=true,
        peak_count=length(peakTime),
        growth_rate_per_s=growth_rate,
        fitted_envelope_ratio=envelope_ratio,
        frequency_hz,
        fit_r_squared
    )
end

function read_response_history(history_path)
    raw,header_matrix = readdlm(history_path,',',header=true)
    headers = String.(vec(header_matrix))
    column(name) = begin
        index = findfirst(==(name),headers)
        isnothing(index) && error("Column '$name' is missing from $history_path")
        Float64.(raw[:,index])
    end
    return (;
        time=column("time_s"),
        tip_displacement=column("tip_displacement_m"),
        tip_twist=column("tip_twist_deg"),
        pitch=column("propeller_1_pitch_deg"),
        yaw=column("propeller_1_yaw_deg")
    )
end

function completed_requested_window(time,requested_end_time)
    isempty(time) && return false
    length(time) == 1 && return time[end] >= requested_end_time - eps(requested_end_time)
    finalStep = time[end] - time[end-1]
    return time[end] >= requested_end_time - 1.01 * finalStep
end

finite_or_blank(value) = isfinite(value) ? @sprintf("%.8g",value) : ""

function csv_field(value)
    text = string(value)
    if occursin(',',text) || occursin('"',text) || occursin('\n',text)
        return "\"" * replace(text,"\"" => "\"\"") * "\""
    end
    return text
end

function write_sweep_summary(path,results)
    headers = [
        "speed_mps",
        "status",
        "reason",
        "completed",
        "max_pitch_deg",
        "max_yaw_deg",
        "pitch_growth_rate_per_s",
        "yaw_growth_rate_per_s",
        "pitch_envelope_ratio",
        "yaw_envelope_ratio",
        "pitch_frequency_hz",
        "yaw_frequency_hz",
        "pitch_fit_r_squared",
        "yaw_fit_r_squared",
        "history_path",
        "log_path"
    ]
    open(path,"w") do stream
        println(stream,join(headers,','))
        for result in results
            values = [
                result.speed_mps,
                result.status,
                result.reason,
                result.completed,
                finite_or_blank(result.max_pitch_deg),
                finite_or_blank(result.max_yaw_deg),
                finite_or_blank(result.pitch.growth_rate_per_s),
                finite_or_blank(result.yaw.growth_rate_per_s),
                finite_or_blank(result.pitch.fitted_envelope_ratio),
                finite_or_blank(result.yaw.fitted_envelope_ratio),
                finite_or_blank(result.pitch.frequency_hz),
                finite_or_blank(result.yaw.frequency_hz),
                finite_or_blank(result.pitch.fit_r_squared),
                finite_or_blank(result.yaw.fit_r_squared),
                result.history_path,
                result.log_path
            ]
            println(stream,join(csv_field.(values),','))
        end
    end
    return path
end

"""Plot the wing and propeller histories from every completed speed case."""
function plot_sweep_time_histories(results,output_directory)
    availableResults = filter(result -> isfile(result.history_path),results)
    isempty(availableResults) && return nothing

    common_style = (
        linewidth=2.0,
        xlabel="Time (s)",
        framestyle=:box,
        gridalpha=0.25
    )
    tipDisplacementPlot = plot(
        ; ylabel="Tip displacement (mm)",title="Wing-tip displacement",
        common_style...
    )
    tipTwistPlot = plot(
        ; ylabel="Tip twist (deg)",title="Wing-tip twist",
        common_style...
    )
    pitchPlot = plot(
        ; ylabel="Pitch (deg)",title="Propeller pitch",
        common_style...
    )
    yawPlot = plot(
        ; ylabel="Yaw (deg)",title="Propeller yaw",
        common_style...
    )

    for (colorIndex,result) in enumerate(availableResults)
        response = read_response_history(result.history_path)
        label = @sprintf("%.0f m/s",result.speed_mps)
        plot!(
            tipDisplacementPlot,
            response.time,
            1.0e3 .* response.tip_displacement;
            color=colorIndex,
            label
        )
        plot!(
            tipTwistPlot,
            response.time,
            response.tip_twist;
            color=colorIndex,
            label
        )
        plot!(
            pitchPlot,
            response.time,
            response.pitch;
            color=colorIndex,
            label
        )
        plot!(
            yawPlot,
            response.time,
            response.yaw;
            color=colorIndex,
            label
        )
    end

    figure = plot(
        tipDisplacementPlot,
        tipTwistPlot,
        pitchPlot,
        yawPlot;
        layout=(2,2),
        size=(1400,900),
        plot_title="Chang speed-sweep time histories"
    )
    output_path = joinpath(output_directory,"speed_sweep_time_histories.png")
    savefig(figure,output_path)
    display(figure)
    println("[sweep] Time-history figure: $output_path")
    return output_path
end

function log_tail(path; line_count=30)
    isfile(path) || return String[]
    lines = readlines(path)
    return lines[max(1,length(lines) - line_count + 1):end]
end

function speed_token(speed)
    return replace(@sprintf("%.3f",speed),"-" => "m","." => "p")
end

function child_environment(speed,output_directory,output_label,endTime,hardAngle)
    environment = copy(ENV)
    environment["CHANG_SPEED_MPS"] = string(speed)
    environment["CHANG_END_TIME_S"] = string(endTime)
    environment["CHANG_WING_SPAN_PANELS"] = get(
        ENV,"CHANG_SWEEP_WING_SPAN_PANELS","30"
    )
    environment["CHANG_WING_CHORD_PANELS"] = get(
        ENV,"CHANG_SWEEP_WING_CHORD_PANELS","5"
    )
    environment["CHANG_PROP_RADIAL_PANELS"] = get(
        ENV,"CHANG_SWEEP_PROP_RADIAL_PANELS","10"
    )
    environment["CHANG_PROP_CHORD_PANELS"] = get(
        ENV,"CHANG_SWEEP_PROP_CHORD_PANELS","10"
    )
    environment["CHANG_ROTATION_RPM"] = get(
        ENV,"CHANG_SWEEP_ROTATION_RPM","1217.6962"
    )
    environment["CHANG_TRIM_SPEED_MPS"] = get(
        ENV,"CHANG_SWEEP_TRIM_SPEED_MPS","65.0"
    )
    environment["CHANG_AOA_DEG"] = get(ENV,"CHANG_SWEEP_AOA_DEG","0.0")
    environment["CHANG_AZIMUTH_STEP_DEG"] = get(
        ENV,"CHANG_SWEEP_AZIMUTH_STEP_DEG","2.5"
    )
    environment["CHANG_INTERACTION"] = get(
        ENV,"CHANG_SWEEP_INTERACTION","false"
    )
    environment["CHANG_NEAR_FIELD_FORCE_MODEL"] = get(
        ENV,"CHANG_SWEEP_BASE_FORCE_MODEL","legacy_imperial_segments"
    )
    environment["CHANG_PROP_MOMENT_PROJECTION"] = get(
        ENV,"CHANG_SWEEP_PROP_MOMENT_PROJECTION","exact_virtual_work"
    )
    environment["CHANG_WAKE_ROWS_WING"] = get(
        ENV,"CHANG_SWEEP_WAKE_ROWS_WING","50"
    )
    environment["CHANG_WAKE_ROWS_PROPELLER"] = get(
        ENV,"CHANG_SWEEP_WAKE_ROWS_PROPELLER","72"
    )
    environment["CHANG_HUB_LOAD_ARM_FACTOR"] = get(
        ENV,"CHANG_SWEEP_HUB_LOAD_ARM_FACTOR","0.5"
    )
    # Preserve this study's factor defaults unless a fixed radius is requested.
    environment["CHANG_CORE_RADIUS_M"] = get(ENV,"CHANG_SWEEP_CORE_RADIUS_M","nothing")
    environment["CHANG_FCORE_SEGMENT_FACTOR"] = get(
        ENV,"CHANG_SWEEP_FCORE_SEGMENT_FACTOR","0.025"
    )
    environment["CHANG_FCORE_CHORD_FACTOR"] = get(
        ENV,"CHANG_SWEEP_FCORE_CHORD_FACTOR","0.01"
    )
    environment["CHANG_PLOT_RESULTS"] = "false"
    environment["CHANG_ANIMATE_WAKE"] = "false"
    environment["CHANG_OUTPUT_DIR"] = output_directory
    environment["CHANG_OUTPUT_LABEL"] = output_label
    environment["CHANG_PROP_ANGLE_ABORT_DEG"] = string(hardAngle)
    return environment
end

function run_speed_case(
    speed;
    output_directory,
    requested_end_time,
    fit_start_s,
    growth_threshold,
    envelope_ratio_threshold,
    minimum_fit_r_squared,
    hard_angle_deg
)
    token = speed_token(speed)
    output_label = "chang_speed_$(token)"
    history_path = joinpath(output_directory,output_label * "_history.csv")
    log_path = joinpath(output_directory,output_label * ".log")
    command = `$(Base.julia_cmd()) --project=$(SWEEP_PACKAGE_DIR) $(SWEEP_CASE_RUNNER)`
    environment = child_environment(
        speed,
        output_directory,
        output_label,
        requested_end_time,
        hard_angle_deg
    )

    println("\n[sweep] Running $(speed) m/s; full output -> $log_path")
    processSucceeded = open(log_path,"w") do stream
        process = run(
            pipeline(setenv(command,environment),stdout=stream,stderr=stream);
            wait=false
        )
        wait(process)
        return success(process)
    end

    unavailable_metrics = (;
        available=false,
        peak_count=0,
        growth_rate_per_s=NaN,
        fitted_envelope_ratio=NaN,
        frequency_hz=NaN,
        fit_r_squared=NaN
    )
    if !processSucceeded || !isfile(history_path)
        reason = processSucceeded ? "history file was not produced" : "child process failed"
        return (;
            speed_mps=speed,
            status="failed",
            reason,
            completed=false,
            max_pitch_deg=NaN,
            max_yaw_deg=NaN,
            pitch=unavailable_metrics,
            yaw=unavailable_metrics,
            history_path,
            log_path
        )
    end

    response = read_response_history(history_path)
    all_finite = all(isfinite,response.time) && all(isfinite,response.pitch) &&
                 all(isfinite,response.yaw)
    # The production grid is 0:dT:Tend, so a completed run may end by less
    # than one dT before Tend when Tend is not an integer multiple of dT.
    completed = completed_requested_window(response.time,requested_end_time)
    maxPitch = isempty(response.pitch) ? NaN : maximum(abs,response.pitch)
    maxYaw = isempty(response.yaw) ? NaN : maximum(abs,response.yaw)
    pitch = oscillation_metrics(response.time,response.pitch; fit_start_s)
    yaw = oscillation_metrics(response.time,response.yaw; fit_start_s)

    hardDivergence = !all_finite || !completed ||
                     max(maxPitch,maxYaw) >= hard_angle_deg
    pitchGrowth =
        pitch.available &&
        pitch.growth_rate_per_s >= growth_threshold &&
        pitch.fitted_envelope_ratio >= envelope_ratio_threshold &&
        pitch.fit_r_squared >= minimum_fit_r_squared
    yawGrowth =
        yaw.available &&
        yaw.growth_rate_per_s >= growth_threshold &&
        yaw.fitted_envelope_ratio >= envelope_ratio_threshold &&
        yaw.fit_r_squared >= minimum_fit_r_squared

    if hardDivergence
        status = "divergent"
        reason =
            !all_finite ? "non-finite response" :
            !completed ? "simulation stopped before requested end time" :
            "hard propeller-angle limit reached"
    elseif pitchGrowth || yawGrowth
        status = "divergent"
        growingChannels = join(
            [name for (name,growing) in (("pitch",pitchGrowth),("yaw",yawGrowth)) if growing],
            " and "
        )
        reason = "$growingChannels envelope is growing"
    elseif !pitch.available && !yaw.available
        status = "indeterminate"
        reason = "too few post-impulse peaks for a growth fit"
    else
        status = "stable"
        reason = "no channel met the divergence thresholds"
    end

    return (;
        speed_mps=speed,
        status,
        reason,
        completed,
        max_pitch_deg=maxPitch,
        max_yaw_deg=maxYaw,
        pitch,
        yaw,
        history_path,
        log_path
    )
end

function print_case_result(result)
    println(
        @sprintf(
            "[sweep] %6.2f m/s  %-13s max |pitch/yaw| = %8.4f / %8.4f deg",
            result.speed_mps,
            uppercase(result.status),
            result.max_pitch_deg,
            result.max_yaw_deg,
        )
    )
    for (name,metrics) in (("pitch",result.pitch),("yaw",result.yaw))
        if metrics.available
            println(
                @sprintf(
                    "        %-5s sigma=%+.5f /s, ratio=%.4f, f=%.4f Hz, R2=%.3f, peaks=%d",
                    name,
                    metrics.growth_rate_per_s,
                    metrics.fitted_envelope_ratio,
                    metrics.frequency_hz,
                    metrics.fit_r_squared,
                    metrics.peak_count,
                )
            )
        else
            println("        $name: insufficient peaks ($(metrics.peak_count))")
        end
    end
    println("        reason: $(result.reason)")
end

function main()
    firstSpeed = sweep_float("CHANG_SWEEP_SPEED_START_MPS",80.0)
    lastSpeed = sweep_float("CHANG_SWEEP_SPEED_STOP_MPS",84.0)
    speed_step = sweep_float("CHANG_SWEEP_SPEED_STEP_MPS",1.0)
    requested_end_time = sweep_float("CHANG_SWEEP_END_TIME_S",3.0)
    fit_start_s = sweep_float("CHANG_SWEEP_FIT_START_S",0.7)
    growth_threshold = sweep_float("CHANG_SWEEP_DIVERGENCE_RATE_PER_S",0.1)
    envelope_ratio_threshold = sweep_float("CHANG_SWEEP_ENVELOPE_RATIO",1.05)
    minimum_fit_r_squared = sweep_float("CHANG_SWEEP_MIN_FIT_R2",0.5)
    hard_angle_deg = sweep_float("CHANG_SWEEP_ABORT_ANGLE_DEG",15.0)
    output_directory = abspath(get(
        ENV,
        "CHANG_SWEEP_OUTPUT_DIR",
        joinpath(SWEEP_EXAMPLE_DIR,"output","speed_sweep_80_84")
    ))

    requested_end_time > fit_start_s || error(
        "CHANG_SWEEP_END_TIME_S must be greater than CHANG_SWEEP_FIT_START_S"
    )
    envelope_ratio_threshold > 1 || error("CHANG_SWEEP_ENVELOPE_RATIO must exceed 1")
    0 <= minimum_fit_r_squared <= 1 || error("CHANG_SWEEP_MIN_FIT_R2 must be in [0, 1]")
    hard_angle_deg > 0 || error("CHANG_SWEEP_ABORT_ANGLE_DEG must be positive")
    mkpath(output_directory)
    summary_path = joinpath(output_directory,"speed_sweep_summary.csv")
    speeds = speed_values(firstSpeed,lastSpeed,speed_step)

    println("Chang speed sweep: $(first(speeds)):$(speed_step):$(last(speeds)) m/s")
    println(
        "Divergence: sigma >= $growth_threshold /s, fitted ratio >= " *
        "$envelope_ratio_threshold, R2 >= $minimum_fit_r_squared, or " *
        "|propeller angle| >= $hard_angle_deg deg"
    )
    println("Results directory: $output_directory")

    results = NamedTuple[]
    for speed in speeds
        result = run_speed_case(
            speed;
            output_directory,
            requested_end_time,
            fit_start_s,
            growth_threshold,
            envelope_ratio_threshold,
            minimum_fit_r_squared,
            hard_angle_deg
        )
        push!(results,result)
        write_sweep_summary(summary_path,results)
        print_case_result(result)

        if result.status in ("divergent","failed")
            println("\n[sweep] STOPPED at $(result.speed_mps) m/s ($(result.status)).")
            println("[sweep] Last lines of the case output:")
            foreach(line -> println("  " * line),log_tail(result.log_path))
            break
        end
    end

    println("\n[sweep] Summary: $summary_path")
    plot_sweep_time_histories(results,output_directory)
    return results
end

# Start when launched as a script. From an included REPL session, call main().
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
