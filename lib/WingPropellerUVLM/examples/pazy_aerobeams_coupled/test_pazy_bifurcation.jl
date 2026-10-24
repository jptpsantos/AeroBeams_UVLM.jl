using Test, DelimitedFiles, Distributed
include(joinpath(@__DIR__, "PazyBifurcation.jl"))
using .PazyBifurcation
using .PazyBifurcation.TimeDomainBifurcation

@testset "Final-window statistics" begin
    t = collect(0:0.001:10)
    q = 3 .+ 2sin.(2pi .* 5 .* t)
    q[t .< 8] .= 1e4 # This startup must not affect any reported statistic.
    s = response_statistics(t, q; transient_time=8, atol=1e-8)
    @test s.minimum ≈ 1
    @test s.maximum ≈ 5
    @test s.mean ≈ 3 atol=1e-12
    @test s.amplitude ≈ 2
    @test s.window_start_s == 8
    @test s.sample_count == 2001
    @test s.converged && s.convergence_status == "converged"
    @test response_statistics(t, fill(12., length(t)); transient_time=8).convergence_status == "stationary"
    growth = 1e5 .+ exp.(t) .* sin.(2pi .* 5 .* t)
    @test response_statistics(t, growth; transient_time=8).convergence_status == "unconverged"
    @test response_statistics(t, t; transient_time=8).convergence_status == "insufficient_cycles"
    @test response_statistics(t, q; transient_time=8, check_convergence=false).convergence_status == "not_checked"
    @test !response_statistics(t, q; transient_time=8, check_convergence=false).converged
    @test_throws ArgumentError response_statistics(t, q; transient_time=10)
    @test_throws ArgumentError response_statistics(t, fill(NaN,length(t)); transient_time=8)
    @test_throws ArgumentError response_statistics([0.,1.,1.], [0.,1.,0.]; transient_time=0)
    @test bifurcation_speeds(20.,22.,1.) == [20.,21.,22.]
    @test bifurcation_speeds(20.,22.,1.5) == [20.,21.5,22.]
    @test bifurcation_speeds(20.,20.,1.) == [20.]
    @test_throws ArgumentError bifurcation_speeds(0,10,1)
    @test_throws ArgumentError bifurcation_speeds(10,20,0)
end

@testset "Reusable driver, CSV, plots and independent callback invocations" begin
    mktempdir() do directory
        calls = Float64[]
        # A pitch output also verifies that extraction is not hardcoded to Pazy.
        definitions = (propeller_pitch=(unit="deg", label="Propeller pitch", atol=1e-6),)
        function simulate(speed)
            push!(calls, speed)
            t = collect(0:0.001:1)
            q = speed .+ sin.(20pi .* t)
            q[t .< 0.5] .= 1e6
            (; time=t, responses=(propeller_pitch=q,), static_equilibrium=(propeller_pitch=speed-0.5,))
        end
        sweep = run_bifurcation(simulate; U_min=20., U_max=22., dU=2., definitions,
            simulation_time=1., transient_time=0.5, save_time_histories=true,
            zoom_limits=(20.5,21.5), output_directory=directory)
        @test calls == [20.,22.]
        @test length(sweep.rows) == 2
        @test all(r -> r.minimum ≈ r.airspeed_mps-1 && r.maximum ≈ r.airspeed_mps+1, sweep.rows)
        @test all(r -> r.static_equilibrium == r.airspeed_mps-0.5, sweep.rows)
        data, header = readdlm(sweep.summary_path, ','; header=true)
        @test size(data) == (2,14)
        @test vec(header)[1:4] == ["airspeed_mps","response","unit","minimum"]
        @test isfile(joinpath(directory,"case_0001_history.csv"))
        @test isfile(joinpath(directory,"case_0002_extrema.csv"))
        paths = plot_bifurcation(sweep.summary_path; zoom_limits=(20.5,21.5))
        @test length(paths) == 3
        @test all(p -> filesize(p)>1000, paths)
        @test_throws ArgumentError run_bifurcation(simulate; U_min=20.,U_max=20.,dU=1.,
            definitions, simulation_time=1., transient_time=0.5, output_directory=directory)
        if "--keep-output" in ARGS
            preview = joinpath(@__DIR__,"..","..","..","..","tmp","bifurcation_plot_validation")
            mkpath(preview)
            foreach(p -> cp(p,joinpath(preview,basename(p)); force=true), paths)
            println("Synthetic plotting validation (not physical results): $preview")
        end
    end
    @test isnothing(run_pazy_bifurcation(; enabled=false, U_min=-1))
    @test_throws ArgumentError run_pazy_bifurcation(; responses=(:propeller_pitch,))
    @test_throws ArgumentError run_pazy_bifurcation(; transient_time=0.12) # pulse ends at 0.14 s
    @test_throws ArgumentError run_pazy_bifurcation(; dynamic_options=(airspeed=42.,))
end

@testset "Two-speed physical smoke test (not settled LCO results)" begin
    equilibrium_options = (; maximum_aerodynamic_time=0.08, minimum_aerodynamic_time=0.02,
        aerodynamic_window_time=0.006, aerodynamic_force_tolerance=5e-4,
        aerodynamic_circulation_tolerance=5e-4, maximum_static_iterations=80,
        static_relaxation=0.5, geometry_tolerance=2e-6, load_tolerance=5e-4,
        progress_frequency=1000)
    short = (; U_min=20., U_max=22., dU=2., simulation_time=0.006, transient_time=0.003,
        angle_of_attack=0.01, chordwise_panels=2, spanwise_panels=4,
        wake_length_chords=2., pulse_amplitude=0., pulse_start=0.,
        responses=(:tip_out_of_plane_m,:tip_in_plane_m,:tip_twist_deg,:root_strain_ue),
        equilibrium_options, dynamic_options=(progress_frequency=1000,), plot_results=false)
    mktempdir() do directory
        sweep = run_pazy_bifurcation(; short..., save_time_histories=true,
            output_directory=joinpath(directory,"serial"))
        @test length(sweep.rows) == 8
        @test sort(unique(getproperty.(sweep.rows,:airspeed_mps))) == [20.,22.]
        for i in 1:2
            history = readdlm(joinpath(sweep.output_directory,"case_$(lpad(i,4,'0'))_history.csv"), ','; skipstart=1)
            @test first(history[:,1]) == 0
            @test last(history[:,1]) ≈ short.simulation_time
            tail = history[:,1] .>= short.transient_time
            for (j,row) in enumerate(sweep.rows[4(i-1)+1:4i])
                @test row.minimum == minimum(history[tail,j+1])
                @test row.maximum == maximum(history[tail,j+1])
                @test row.window_start_s >= short.transient_time
                @test isfinite(row.static_equilibrium)
            end
        end
        # Repeating U=20 after U=22 must reproduce its own initialization and history.
        repeat = run_pazy_bifurcation(; short..., U_max=20., save_time_histories=true,
            output_directory=joinpath(directory,"repeat"))
        a = readdlm(joinpath(sweep.output_directory,"case_0001_history.csv"), ','; skipstart=1)
        b = readdlm(joinpath(repeat.output_directory,"case_0001_history.csv"), ','; skipstart=1)
        @test a ≈ b rtol=1e-10 atol=1e-12
        @test !isfile(joinpath(directory,"repeat","case_0002_history.csv"))
        # Exercise separate-process solves and compare compact summaries to serial.
        if "--parallel" in ARGS
            ids = addprocs(2; exeflags=`--startup-file=no --compiled-modules=existing --project=$(dirname(Base.active_project()))`)
            try
                parallel = run_pazy_bifurcation(; short..., worker_ids=ids,
                    output_directory=joinpath(directory,"parallel"))
                @test getproperty.(parallel.rows,:minimum) ≈ getproperty.(sweep.rows,:minimum)
                @test getproperty.(parallel.rows,:maximum) ≈ getproperty.(sweep.rows,:maximum)
                @test !isfile(joinpath(parallel.output_directory,"case_0001_history.csv"))
            finally
                rmprocs(ids)
            end
        end
    end
end
