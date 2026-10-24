# Independent ramps with optional tip-force pulses in either or both branches.
# The original unforced run_pazy_LCO.jl remains a separate entry point.
include(joinpath(@__DIR__, "PazyLinearLCOSweep.jl"))

# Set PAZY_LCO_DOWNSWEEP_PULSE_RUN=false before including to load without running.
if get(ENV, "PAZY_LCO_DOWNSWEEP_PULSE_RUN", "true") == "true"
    # Select the analyses to execute independently from plotting.
    run_branches = :both                 # :both, :upsweep, or :downsweep
    plot_only = false                    # true: regenerate strain-versus-speed plot only
    plot_branches = run_branches         # or explicitly :both, :upsweep, :downsweep
    pulse_branches = :downsweep          # :both, :upsweep, :downsweep, or :none
    # Keep each pulse selection separate from existing results.
    output_name = pulse_branches == :downsweep ? "LCO_downsweep_pulse" :
        pulse_branches == :both ? "LCO_both_pulses" :
        pulse_branches == :upsweep ? "LCO_upsweep_pulse" : "LCO_no_pulse"
    output_directory = joinpath(@__DIR__, "output", output_name)

    # Shared UVLM settings for each branch's static initialization and dynamics.
    uvlm_options = (;
        chordwise_panels=10,
        spanwise_panels=20,
        time_step_chords=0.25,           # maximum U*dt/chord, at the highest speed
        fixed_time_step=nothing,        # optional dt [s]; must align both ramp endpoints
        wake_length_chords=5.0,          # sizes wake capacity at the minimum speed
        maximum_wake_rows=nothing,      # integer overrides wake_length_chords sizing
        core_radius=1e-3,               # m
        wake_shedding_fraction=0.1,      # trailing-edge offset fraction of U*dt
        symmetric_wing=true,
        aerodynamic_interaction=true)

    sweep_initial_speed = 35.0           # m/s, lower endpoint of both branches
    sweep_peak_speed = 52.0              # m/s, upper endpoint of both branches
    upsweep_duration = 40.0              # s, independent upsweep
    downsweep_duration = 20.0            # s, independent downsweep with pulse

    # Excite the wing when the descending speed reaches this value.
    pulse_speed = 49.0                   # m/s, diagnostic trigger
    pulse_coefficient = 1.0              # N; waveform peak is about 0.6495 times this
    pulse_duration = 0.1                # s; smooth bipolar pulse, zero net impulse
    # Separate ascending trigger; 45 m/s is a diagnostic choice, not a boundary.
    upsweep_pulse_speed = 45.0           # m/s
    upsweep_pulse_coefficient = 1.0      # N
    upsweep_pulse_duration = 0.1        # s
    pulse_branches in (:both, :upsweep, :downsweep, :none) ||
        throw(ArgumentError("pulse_branches must be :both, :upsweep, :downsweep, or :none"))
    if plot_only
        pulse_lco_plot = plot_pazy_LCO(; output_directory, plot_branches)
    else
        run_names = pazy_lco_branch_names(run_branches; option_name="run_branches")
        pulse_names = pulse_branches == :none ? () :
            pazy_lco_branch_names(pulse_branches; option_name="pulse_branches")
        pulse_start = 0.0
        upsweep_pulse_start = 0.0
        if "downsweep" in run_names && "downsweep" in pulse_names && !iszero(pulse_coefficient)
            @assert sweep_initial_speed < pulse_speed <= sweep_peak_speed
            pulse_start = (sweep_peak_speed-pulse_speed) /
                (sweep_peak_speed-sweep_initial_speed) * downsweep_duration
        end
        if "upsweep" in run_names && "upsweep" in pulse_names && !iszero(upsweep_pulse_coefficient)
            @assert sweep_initial_speed <= upsweep_pulse_speed < sweep_peak_speed
            upsweep_pulse_start = (upsweep_pulse_speed-sweep_initial_speed) /
                (sweep_peak_speed-sweep_initial_speed) * upsweep_duration
        end
        # Alternatively, set pulse_start=0.1 to excite near the branch beginning.
        pulse_lco_results = run_pazy_linear_LCO(;
            initial_speed=sweep_initial_speed, peak_speed=sweep_peak_speed,
            downsweep_initial_speed=sweep_peak_speed,
            upsweep_duration, downsweep_duration,
            angle_of_attack=deg2rad(5.0),
            upsweep_pulse_amplitude="upsweep" in pulse_names ? upsweep_pulse_coefficient : 0.0,
            upsweep_pulse_start, upsweep_pulse_duration,
            downsweep_pulse_amplitude="downsweep" in pulse_names ? pulse_coefficient : 0.0,
            downsweep_pulse_start=pulse_start,
            downsweep_pulse_duration=pulse_duration,
            uvlm_options..., run_branches, plot_branches, output_directory)
    end
end
