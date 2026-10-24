using Test, DelimitedFiles
include(joinpath(@__DIR__, "..", "WingPropellerBifurcation.jl"))
using .WingPropellerBifurcation
const WPA = WingPropellerBifurcation.WingPropellerAeroelastic

@testset "Wing-propeller bifurcation adapter" begin
    case = WPA.xu_case(; element_count=2, wing_spanwise_panels=2,
        wing_chordwise_panels=1, propeller_radial_panels=2, propeller_chordwise_panels=1)
    base_speed = case.operating_condition.freestream_speed
    options = (; simulation_time=0.0003, time_step_s=0.0001, azimuth_step_deg=5.,
        solver_options=(maximum_wake_rows_wing=3, maximum_wake_rows_propeller=3,))
    first_run = WingPropellerBifurcation.simulate_wing_propeller_speed(60.,case,options)
    second_run = WingPropellerBifurcation.simulate_wing_propeller_speed(62.,case,options)
    @test second_run.angular_speeds ≈ first_run.angular_speeds .* (62/60)
    @test case.operating_condition.freestream_speed == base_speed
    @test first_run.time[1] == second_run.time[1] == 0
    @test last(first_run.time) ≈ options.simulation_time
    @test first_run.responses.tip_out_of_plane_m !== second_run.responses.tip_out_of_plane_m
    @test all(iszero, first(v) for v in values(first_run.responses))
    @test all(isfinite, first_run.responses.propeller_1_pitch_rad)
    @test all(isfinite, first_run.responses.propeller_1_yaw_rad)
    mktempdir() do directory
        sweep = run_wing_propeller_bifurcation(case;
            U_min=60., U_max=62., dU=2., simulation_time=0.0003, transient_time=0.0001,
            time_step_s=options.time_step_s, solver_options=options.solver_options,
            plot_results=false, save_time_histories=false, output_directory=directory)
        @test length(sweep.rows) == 10
        @test all(r -> isfinite(r.minimum) && isfinite(r.maximum), sweep.rows)
        @test all(r -> isnan(r.static_equilibrium), sweep.rows)
        @test !isfile(joinpath(directory,"case_0001_history.csv"))
        data = readdlm(sweep.summary_path, ','; skipstart=1)
        @test size(data,1) == 10
    end
end
