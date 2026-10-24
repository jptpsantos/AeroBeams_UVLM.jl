# One continuous upsweep/downsweep, preserving beam motion and wake at reversal.
include(joinpath(@__DIR__, "PazyContinuousLCOSweep.jl"))

# Set PAZY_CONTINUOUS_LCO_RUN=false before including to load without running.
# Results are separate from run_pazy_LCO.jl in output/LCO_continuous_sweeps.
if get(ENV, "PAZY_CONTINUOUS_LCO_RUN", "true") == "true"
    sweep_initial_speed = 35.0           # m/s
    sweep_peak_speed = 52.0              # m/s, reversal speed
    downsweep_deceleration = 0.85        # m/s^2, magnitude
    # Match dt and wake capacity of the completed independent 60 m/s experiment.
    # Set nothing to size dt automatically from the current ramp parameters.
    comparison_time_step = 0.00041207376120325537
    continuous_lco_results = run_pazy_continuous_LCO(;
        initial_speed=sweep_initial_speed, peak_speed=sweep_peak_speed,
        upsweep_duration=40.0,
        downsweep_duration=(sweep_peak_speed-sweep_initial_speed)/downsweep_deceleration,
        fixed_time_step=comparison_time_step,
        angle_of_attack=deg2rad(5.0))
end
