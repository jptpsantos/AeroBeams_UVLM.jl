# Shared moving-block postprocessing. The conservative numerical fingerprint and
# a separate estimator fingerprint identify the code used for each study.
struct MovingBlockStudy{S}
    settings::S
    estimator_sha256::String
end
Base.getproperty(s::MovingBlockStudy,key::Symbol) = key in (:settings,:estimator_sha256) ?
    getfield(s,key) : getproperty(getfield(s,:settings),key)
Base.propertynames(s::MovingBlockStudy) = (propertynames(s.settings)...,:estimator_sha256)
function moving_block_fingerprint()
    paths=(joinpath(@__DIR__,"MovingBlockDamping.jl"),@__FILE__)
    return bytes2hex(Convergence.sha256(join(read.(paths,String),'\0')))
end
MovingBlockStudy(s) = MovingBlockStudy(s,moving_block_fingerprint())
Convergence.data(s::MovingBlockStudy) = Convergence.data(merge(s.settings,
    (;damping_method="moving_block_2",damping_method_sha256=s.estimator_sha256)))

function moving_block_options(d)
    mode=get(d,:block_size_mode,:record_fraction)
    mode in (:record_fraction,:duration) || throw(ArgumentError("block_size_mode must be :record_fraction or :duration"))
    return (;block_size=get(d,:block_size,512),size_ratio_lb=get(d,:size_ratio_lb,.25),
        size_ratio_ub=get(d,:size_ratio_ub,.5),peak_from_start=get(d,:peak_from_start,1),
        peak_from_end=get(d,:peak_from_end,0),minimum_peaks=d.minimum_peaks,
        minimum_fit_r_squared=d.minimum_fit_r_squared,
        block_duration_s=mode==:duration ? d.block_duration_s : nothing,
        block_overlap=mode==:duration ? d.block_overlap : nothing,
        frequency_band_hz=d.frequency_band_hz)
end

function Convergence.validate_elastic(s::MovingBlockStudy,selection)
    d=s.damping
    options=moving_block_options(d)
    # Validate FFT settings before any simulation (an empty signal is sufficient).
    MovingBlockDamping.moving_block_metrics(Float64[],Float64[];options...)
    duration=isnothing(options.block_duration_s) ?
        (s.response.end_time_s-Convergence.fit_start(s,selection.physical))*options.size_ratio_lb : options.block_duration_s
    # Reuse the original physical/model validation with a representative window.
    # Its legacy overlap and window checks must not constrain record-fraction mode.
    legacy=merge(s.settings,(;damping=merge(d,(;block_duration_s=duration,
        block_overlap=something(options.block_overlap,0.0)))))
    return Convergence.validate_elastic(legacy,selection)
end

function write_moving_block_diagnostics(directory,label,result)
    d=result.diagnostics
    open(joinpath(directory,"$(label)_moving_block.csv"),"w") do io
        println(io,"block_start_s,log_amplitude,fitted_log_amplitude,frequency_hz")
        Convergence.writedlm(io,hcat(d.block_start_s,d.log_amplitude,
            d.fitted_log_amplitude,d.frequency_hz),',')
    end
    return nothing
end
function plot_moving_block(directory,label,result)
    d=result.diagnostics
    p1=Plots.plot(d.time_s,d.signal;label="Analyzed segment",xlabel="Time (s)",
        ylabel="Angle (deg)",title="$label response")
    p2=Plots.plot(d.block_start_s,d.log_amplitude;label="Moving-block log amplitude",
        xlabel="Block start time (s)",ylabel="Log amplitude")
    Plots.plot!(p2,d.block_start_s,d.fitted_log_amplitude;label="Linear fit",linestyle=:dash)
    Plots.savefig(Plots.plot(p1,p2;layout=(2,1),size=(900,650)),
        joinpath(directory,"$(label)_moving_block.png"))
end

function Convergence.damping_result(history,s::MovingBlockStudy,p)
    moving_block_fingerprint()==s.estimator_sha256 || error("Moving-block source changed during the study")
    response=Convergence.Damping.read_response_history(history)
    Convergence.Damping.completed_requested_window(response.time,s.response.end_time_s) ||
        error("Response stopped before the requested end time")
    options=moving_block_options(s.damping)
    metric(signal)=MovingBlockDamping.moving_block_metrics(response.time,signal;
        fit_start_s=Convergence.fit_start(s,p),fit_end_s=s.response.end_time_s,options...)
    pitch=metric(response.pitch); yaw=metric(response.yaw)
    directory=dirname(history)
    for (label,result) in (("pitch",pitch),("yaw",yaw))
        write_moving_block_diagnostics(directory,label,result)
        if s.make_plots && result.summary.available
            @eval using Plots
            Base.invokelatest(plot_moving_block,directory,label,result)
        end
    end
    amplitude_ok=max(maximum(abs,response.pitch),maximum(abs,response.yaw))<s.integration.propeller_angle_limit_deg
    return (;pitch=pitch.summary,yaw=yaw.summary,
        valid=pitch.summary.valid_for_convergence && yaw.summary.valid_for_convergence && amplitude_ok,
        method="moving_block_2",method_sha256=s.estimator_sha256)
end

# The previously documented nested entry now uses the same operating-point and
# moving-block path as main(). Specialize on study settings without overwriting
# the generic backend implementation.
Convergence.run_aeroelastic(s::NamedTuple) = run_aeroelastic(s)
