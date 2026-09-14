function load_model()
    if !isdefined(@__MODULE__,:Model)
        @eval module Model
            include(joinpath(@__DIR__,"..","..","..","src","ChangAeroelastic.jl"))
        end
    end
    return nothing
end
function model_defaults(p,mode,c,s,directory)
    f=core_factors(mode,c.core)
    return (;
        wing=(root_chord_m=p.wing_root_chord_m,tip_chord_m=p.wing_tip_chord_m,span_m=p.wing_span_m,
            symmetric=p.wing_symmetric,spanwise_panels=c.wing_span,chordwise_panels=c.wing_chord),
        propeller=(radius_m=p.propeller_radius_m,chord_m=p.propeller_chord_m,blades=p.blades,
            radial_panels=c.prop_radial,chordwise_panels=c.prop_chord,rotation_rpm=p.rpm,
            trim_speed_mps=p.speed_mps,attachment_eta=[p.attachment_eta],collective_pitch_offset_deg=p.collective_offset_deg),
        simulation=(air_density_kgpm3=p.density_kgpm3,freestream_speed_mps=p.speed_mps,
            angle_of_attack_deg=p.alpha_deg,sideslip_deg=p.beta_deg,azimuth_step_deg=c.azimuth_deg,
            end_time_s=s.response.end_time_s,interaction_on=p.interaction,near_field_force_model=:imperial,
            propeller_moment_projection=s.propeller_moment_projection,impulse_propeller_indices=[1]),
        aerodynamic=(core_radius_m=f.radius,segment_core_factor=f.segment,chord_core_factor=f.chord,
            elastic_axis_fraction=p.elastic_axis_fraction,hub_load_arm_factor=s.hub_load_arm_factor),
        wake=(maximum_rows_wing=ceil(Int,c.wake_revolutions*360/c.azimuth_deg),wing_rows_per_chord_panel=10,
            maximum_rows_propeller=ceil(Int,c.wake_revolutions*360/c.azimuth_deg)),
        excitation=(trim_revolutions=s.response.trim_revolutions,trim_average_revolutions=s.response.trim_average_revolutions,
            impulse_start_s=s.response.impulse_start_s,impulse_duration_s=s.response.impulse_duration_s,
            impulse_magnitude_nm=s.response.impulse_magnitude_nm),
        integration=s.integration,coupling=s.coupling,
        structural=merge(s.structural,(;pylon_length_m=p.pylon_length_m)),
        output=(directory=abspath(directory),label="response",plot_results=false,plot_time_limit_s=s.response.end_time_s,
            animate_wake=false,animation_stride=5,animation_fps=15,
            animation_axis_limits=((-3.,5.),(0.,8.),(-4.,4.)),animation_tick_spacing_m=1.))
end
function fit_start(s,p)
    pulse_start=isnothing(s.response.impulse_start_s) ? s.response.trim_revolutions*60/p.rpm : s.response.impulse_start_s
    automatic=max(.75,pulse_start+s.response.impulse_duration_s+60/p.rpm)
    return something(s.damping.fit_start_s,automatic)
end
function elastic_sweeps(selected,s)
    pairs=map(s.families) do k
        k in FIELDS || error("Unknown family $k")
        hasproperty(s.sweeps,k) || error("Missing aeroelastic sweep $k")
        value=getproperty(s.sweeps,k); base=getproperty(selected,k)
        if isnothing(value)
            value=k in MESH ? [base,base+max(1,cld(base,4)),base+max(2,cld(base,2))] :
                k==:wake_revolutions ? [base,base+1,base+2] : [base,base/2,base/4]
        end
        k=>value
    end
    return (;pairs...)
end
function validate_elastic(s,selection)
    # The production solver currently fixes this shedding fraction at 0.1.
    selection.physical.wake_shedding_fraction==.1 || error("Production aeroelastic solver requires wake_shedding_fraction=0.1; rerun the aerodynamic stage with that value")
    cases=build_cases(selection.selected,elastic_sweeps(selection.selected,s),s.families;elastic=true)
    d=s.damping; start=fit_start(s,selection.physical)
    0<=start<s.response.end_time_s || error("Damping fit must start inside the response interval")
    0<d.block_duration_s<s.response.end_time_s-start || error("Damping window must fit inside the post-excitation interval")
    0<=d.block_overlap<1 && d.minimum_peaks>=3 && 0<=d.minimum_fit_r_squared<=1 || error("Invalid damping fit controls")
    all(x->isfinite(x)&&x>=0,(d.lambda_tolerance_per_s,d.frequency_tolerance_hz)) || error("Invalid damping tolerances")
    lo,hi=d.frequency_band_hz
    isfinite(lo) && isfinite(hi) && 0<=lo<hi || error("Invalid damping frequency band")
    for c in cases
        hi<3selection.physical.rpm/c.azimuth_deg || error("Damping frequency band reaches Nyquist")
        # Full explicit configuration validation before any costly simulation.
        defaults=model_defaults(selection.physical,selection.core_mode,c,s,s.output_directory)
        Model.ChangAeroelastic.load_chang_configuration(defaults;env=Dict{String,String}())
    end
    return cases
end
function damping_result(history,s,p)
    response=Damping.read_response_history(history)
    all(isfinite,response.time) && all(isfinite,response.pitch) && all(isfinite,response.yaw) || error("Nonfinite aeroelastic response")
    all(>(0),diff(response.time)) || error("Nonmonotone response times")
    Damping.completed_requested_window(response.time,s.response.end_time_s) || error("Response stopped before the requested end time")
    d=s.damping
    metric(signal)=Damping.moving_block_metrics(response.time,signal;
        fit_start_s=fit_start(s,p),fit_end_s=s.response.end_time_s,
        minimum_peaks=d.minimum_peaks,minimum_fit_r_squared=d.minimum_fit_r_squared,
        initial_block_size=512,size_ratio_lower=.25,size_ratio_upper=.5,
        peak_from_start=1,peak_from_end=0,block_duration_s=d.block_duration_s,block_overlap=d.block_overlap,
        frequency_min_hz=d.frequency_band_hz[1],frequency_max_hz=d.frequency_band_hz[2],apply_hann_window=true)
    pitch=metric(response.pitch); yaw=metric(response.yaw)
    amplitude_ok=max(maximum(abs,response.pitch),maximum(abs,response.yaw))<s.integration.propeller_angle_limit_deg
    valid=pitch.valid_for_convergence && yaw.valid_for_convergence && amplitude_ok
    return (;pitch,yaw,valid)
end
function aeroelastic_case(c,s,selection,directory)
    mkpath(directory)
    defaults=model_defaults(selection.physical,selection.core_mode,c,s,directory)
    config=Model.ChangAeroelastic.load_chang_configuration(defaults;env=Dict{String,String}())
    save_toml(joinpath(directory,"resolved_configuration.toml"),config)
    log=joinpath(directory,"run.log")
    open(log,"w") do io
        redirect_stdout(io) do
            redirect_stderr(io) do
                Model.ChangAeroelastic.run_chang(config)
            end
        end
    end
    history=joinpath(directory,"response_history.csv")
    fits=damping_result(history,s,selection.physical)
    save_toml(joinpath(directory,"damping_fit.toml"),fits)
    return (;status=fits.valid ? "completed" : "indeterminate",reason=fits.valid ? "Pitch/yaw damping fits verified" : "Damping fit unavailable, poor quality, or amplitude limit exceeded",
        values=(fits.pitch.moving_block_lambda_per_s,fits.yaw.moving_block_lambda_per_s,
            fits.pitch.frequency_hz,fits.yaw.frequency_hz),phases=(),periodic=(),valid=fits.valid,history,log,reused=false)
end
function run_aeroelastic(s)
    selection=load_selection(s.selection_file)
    load_model()
    return Base.invokelatest(run_aeroelastic_loaded,s,selection)
end
function run_aeroelastic_loaded(s,selection)
    cases=validate_elastic(s,selection)
    directory=abspath(s.output_directory); mkpath(directory)
    save_toml(joinpath(directory,"aerodynamic_input.toml"),selection)
    results=execute_cases(cases,s,selection.physical,selection.core_mode;elastic=true,
        runner=(c,folder)->aeroelastic_case(c,s,selection,folder))
    rows=annotate(results,cases,s;elastic=true)
    if !s.dry_run && s.make_plots
        @eval using Plots
        Base.invokelatest(plot_results,directory,rows;elastic=true)
    end
    write(joinpath(directory,"report.md"),"# Aeroelastic convergence\n\nAerodynamic selection: $(s.selection_file)\n\nWing span/structural elements: $(selection.selected.wing_span). Positive lambda means growth; negative lambda means decay. Compare both pitch/yaw growth rates and tracked frequencies, and inspect damping_fit.toml and raw histories. This is a damping sensitivity study at one operating point, not a flutter boundary calculation.\n")
    println("Results: $directory")
    return (;directory,cases,results=rows,selection)
end
