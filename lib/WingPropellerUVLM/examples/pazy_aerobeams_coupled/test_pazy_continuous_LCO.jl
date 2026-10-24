using Test, DelimitedFiles
include(joinpath(@__DIR__, "PazyContinuousLCOSweep.jl"))
BLAS.set_num_threads(1)

@testset "Continuous LCO reversal preserves the evolving state" begin
    @test pazy_continuous_airspeed(0.;initial_speed=20.,peak_speed=22.,
        upsweep_duration=0.006,downsweep_duration=0.003) == 20.
    @test pazy_continuous_airspeed(0.006;initial_speed=20.,peak_speed=22.,
        upsweep_duration=0.006,downsweep_duration=0.003) == 22.
    @test pazy_continuous_airspeed(0.009;initial_speed=20.,peak_speed=22.,
        upsweep_duration=0.006,downsweep_duration=0.003) ≈ 20.
    @test_throws ArgumentError pazy_lco_aligned_timestep(0.006,0.003,0.001;
        fixed_time_step=0.0007)
    @test_throws ArgumentError pazy_lco_aligned_timestep(0.006,0.003,0.001;
        fixed_time_step=0.002)
    options = (;initial_speed=20.,peak_speed=22.,upsweep_duration=0.006,
        downsweep_duration=0.003,angle_of_attack=0.01,chordwise_panels=2,
        spanwise_panels=4,wake_length_chords=2.,fixed_time_step=0.0005,
        progress_frequency=10000,
        equilibrium_options=(;maximum_aerodynamic_time=0.08,minimum_aerodynamic_time=0.02,
            aerodynamic_window_time=0.006,aerodynamic_force_tolerance=5e-4,
            aerodynamic_circulation_tolerance=5e-4,maximum_static_iterations=80,
            static_relaxation=0.5,geometry_tolerance=2e-6,load_tolerance=5e-4))
    mktempdir() do directory
        continuous = run_pazy_continuous_LCO(;options...,
            output_directory=joinpath(directory,"continuous"))
        reference = run_pazy_linear_LCO(;options...,
            output_directory=joinpath(directory,"independent"))
        result,up,down = continuous.result,continuous.upsweep,continuous.downsweep
        turn,nextstep = continuous.turnaround,continuous.first_downsweep_step
        @test result.dt == continuous.dt == reference.dt == 0.0005
        @test result.time[1] == 0.0
        @test last(result.time) ≈ 0.009
        @test last(up.time) == first(down.time) ≈ 0.006
        @test first(down.branch_time) == 0.0
        @test last(down.branch_time) ≈ 0.003
        @test last(up.airspeed_history) == first(down.airspeed_history) == 22.0
        @test last(up.root_strain_ue) == first(down.root_strain_ue)
        @test all(diff(up.airspeed_history) .> 0)
        @test all(diff(down.airspeed_history) .< 0)
        @test length(up.time)+length(down.time)-1 == length(result.time)
        @test up.root_strain_ue ≈ reference.upsweep.root_strain_ue rtol=1e-10 atol=1e-10
        @test up.tip_out_of_plane ≈ reference.upsweep.result.tip_out_of_plane rtol=1e-10 atol=1e-12
        @test abs(first(down.root_strain_ue)-first(reference.downsweep.root_strain_ue)) > 1e-4
        @test turn.maximum_velocity > 1e-10
        @test turn.aerodynamic_time == turn.time
        @test turn.step_index == turn.index
        @test nextstep.step_index == turn.step_index+1
        @test nextstep.aerodynamic_time == nextstep.time > turn.time
        @test nextstep.nwake == min.(turn.aerodynamic_snapshot.nwake .+ 1,
            continuous.maximum_wake_rows)
        @test result.aerodynamic_problem.state.stepIndex == length(result.time)
        @test result.aerodynamic_problem.state.system.nwake[1] ==
            min(result.initial_equilibrium_state.aerodynamic_snapshot.nwake[1]+
                length(result.time)-1,continuous.maximum_wake_rows)
        @test maximum(result.coupling_geometry_residual) <= 1e-4
        @test maximum(result.coupling_load_residual) <= 1e-4
        for (branch,name) in ((up,"upsweep"),(down,"downsweep"))
            data = readdlm(joinpath(continuous.output_directory,"pazy_linear_LCO_$name.csv"),
                ',',Float64;skipstart=1)
            @test data[:,1] == branch.time
            @test data[:,2] == branch.airspeed_history
            @test data[:,3] == branch.root_strain_ue
            @test data[:,end] == branch.branch_time
        end
        data = readdlm(joinpath(continuous.output_directory,"pazy_linear_LCO_continuous_history.csv"),
            ',',Float64;skipstart=1)
        @test data[:,1] == result.time
        @test data[:,3] ≈ pazy_root_strain_ue(result,0.,0.001125) atol=1e-10
        @test !isfile(joinpath(continuous.output_directory,"pazy_linear_LCO_downsweep_initialization.csv"))
        for name in ("pazy_linear_LCO_root_strain_vs_airspeed.png",
                     "pazy_linear_LCO_root_strain_history.png","pazy_linear_LCO_turnaround.csv")
            @test filesize(joinpath(continuous.output_directory,name)) > 100
        end
        plot_pazy_LCO(; output_directory=continuous.output_directory)
        if "--keep-output" in ARGS
            preview = mktempdir(joinpath(REPOSITORY_ROOT,"tmp");
                prefix="pazy_continuous_validation_",cleanup=false)
            for name in readdir(continuous.output_directory)
                cp(joinpath(continuous.output_directory,name),joinpath(preview,name))
            end
            println("Short verification outputs (not physical LCO results): $preview")
        end
    end
end
