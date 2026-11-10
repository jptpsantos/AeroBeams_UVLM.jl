include(joinpath(@__DIR__,"src","ChangAeroBeams.jl"))
using .ChangAeroBeams

if abspath(PROGRAM_FILE) == @__FILE__
    case = include(joinpath(@__DIR__,"chang_aerobeams_case.jl"))
    case = merge(case,(
        reference=merge(case.reference,(simulation=merge(case.reference.simulation,(end_time_s=5.0,)),)),
        solver=merge(case.solver,(impulse_start_s=1.0,)),))
    speed_sweep = run_chang_aerobeams_speed_sweep(case,50.0:10.0:90.0;
        fit_start_s=2.0,window_duration_s=0.2)
end
