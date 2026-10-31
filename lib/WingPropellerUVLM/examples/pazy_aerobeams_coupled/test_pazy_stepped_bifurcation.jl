using Test, DelimitedFiles, LinearAlgebra
include(joinpath(@__DIR__,"PazySteppedBifurcation.jl"))
using .PazySteppedBifurcation
const P=PazySteppedBifurcation

@testset "Fixed-speed schedule and boundary ownership" begin
    plan=stepped_plan(;U_min=20.,U_max=22.,dU=1.,plateau_time=0.01,transient_time=0.005,
        chord=0.1,time_step_chords=0.25,fixed_time_step=0.0005)
    @test getproperty.(plan.plateaus,:airspeed_mps) == [20.,21.,22.,22.,21.,20.]
    @test getproperty.(plan.plateaus,:direction) == [:upsweep,:upsweep,:upsweep,:downsweep,:downsweep,:downsweep]
    @test plan.total_steps == 120
    @test plan.total_time ≈ 0.06
    @test plateau_at_step(plan,0) == 1
    @test plateau_at_step(plan,20) == 1
    @test plateau_at_step(plan,21) == 2
    @test plateau_at_step(plan,60) == 3
    @test plateau_at_step(plan,61) == 4
    @test plateau_at_step(plan,120) == 6
    @test_throws ArgumentError plateau_at_step(plan,121)
    for p in plan.plateaus
        @test all(plateau_at_step(plan,k) == p.index for k in p.start_step+1:p.end_step)
    end
    for (order,speeds) in ((:down_up,[22.,20.,20.,22.]),(:upsweep,[20.,22.]),(:downsweep,[22.,20.]))
        p=stepped_plan(;U_min=20.,U_max=22.,dU=2.,plateau_time=0.01,transient_time=0.005,chord=0.1,sweep_order=order)
        @test getproperty.(p.plateaus,:airspeed_mps) == speeds
    end
    p=stepped_plan(;U_min=20.,U_max=22.5,dU=1.,plateau_time=0.01,transient_time=0.,chord=0.1)
    @test p.grid == [20.,21.,22.,22.5]
    @test p.observation_start == p.dt
    @test p.dt <= p.nominal_dt
    @test_throws ArgumentError stepped_plan(;chord=0.1,fixed_time_step=0.0007)
    @test_throws ArgumentError stepped_plan(;chord=0.1,fixed_time_step=1.)
    @test_throws ArgumentError stepped_plan(;chord=0.1,transient_time=10.)
    @test_throws ArgumentError stepped_plan(;chord=0.1,sweep_order=:other)
end

@testset "Derivative-zero peak refinement and transient exclusion" begin
    t=collect(0.:0.1:1.)
    q=3 .- 2 .* (t .- 0.53).^2
    peak=only(quadratic_extrema(t,q;transient_time=0.))
    @test peak.kind == "maximum"
    @test peak.time_s ≈ 0.53
    @test peak.value ≈ 3.
    low=only(quadratic_extrema(t,-q;transient_time=0.))
    @test low.kind == "minimum"
    @test low.time_s ≈ 0.53
    @test low.value ≈ -3.
    uneven=[0.,0.12,0.23,0.4,0.57,0.76,1.]
    peak=only(quadratic_extrema(uneven,3 .- 2 .* (uneven .- 0.53).^2;transient_time=0.))
    @test peak.time_s ≈ 0.53
    @test peak.value ≈ 3.
    @test isempty(quadratic_extrema(t,fill(2.,length(t));transient_time=0.))
    @test isempty(quadratic_extrema(t,t;transient_time=0.))
    @test isempty(quadratic_extrema(t,q;transient_time=0.6))
    flat=quadratic_extrema(collect(0.:5.),[0.,1.,2.,2.,1.,0.];transient_time=0.)
    @test only(flat).time_s == 3.
    @test only(flat).value == 2.
end

@testset "Per-plateau statistics and refined envelope" begin
    plan=stepped_plan(;U_min=20.,U_max=22.,dU=2.,plateau_time=2.,transient_time=1.,
        chord=0.1,time_step_chords=0.25)
    response=:tip_out_of_plane_percent
    q=20 .+ sin.(8pi.*plan.local_time)
    q[plan.local_time .< 1.] .+= 100. # This startup must never enter the envelope.
    static=(;tip_out_of_plane_percent=20.)
    first=P.analyze_plateau(plan,plan.plateaus[1],Dict(response=>q),static;
        responses=(response,),angle_of_attack=0.,extrema_method=:quadratic,
        check_convergence=true,rtol=0.02,convergence_cycles=3,response_definitions=P.DEFINITIONS)
    row=only(first.rows)
    @test row.minimum ≈ 19. atol=1e-5
    @test row.maximum ≈ 21. atol=1e-5
    @test row.minimum <= row.sample_minimum
    @test row.maximum >= row.sample_maximum
    @test row.converged
    @test all(p.time_s >= 1. for p in first.peaks)
    @test row.window_start_s >= 1.
    @test row.static_equilibrium == 20.
    second=P.analyze_plateau(plan,plan.plateaus[2],Dict(response=>q),static;
        responses=(response,),angle_of_attack=0.,extrema_method=:sampled,
        check_convergence=true,rtol=0.02,convergence_cycles=3,response_definitions=P.DEFINITIONS)
    row=only(second.rows)
    @test isnan(row.static_equilibrium)
    @test row.minimum == row.sample_minimum
    @test row.global_window_start_s >= 3.
    @test all(p.global_time_s >= 3. for p in second.peaks)
    @test isnothing(run_pazy_stepped_bifurcation(;enabled=false,U_min=-1.))
    @test_throws ArgumentError run_pazy_stepped_bifurcation(;dynamic_options=(airspeed=40.,))
    @test_throws ArgumentError run_pazy_stepped_bifurcation(;responses=(:propeller_pitch,))
    @test_throws ArgumentError run_pazy_stepped_bifurcation(;pulse_start=8.)
    @test !isdefined(P,:PazyBifurcation) # Pure checks do not load physics.
end

if "--physical" in ARGS
    @testset "One physical integration preserves state across plateaus" begin
        ENV["GKSwstype"]="100"
        BLAS.set_num_threads(1)
        root=normpath(joinpath(@__DIR__,"..","..","..","..","tmp"))
        mkpath(root)
        directory=mktempdir(root;prefix="pazy_stepped_smoke_",cleanup=false)
        equilibrium_options=(;maximum_aerodynamic_time=0.08,minimum_aerodynamic_time=0.02,
            aerodynamic_window_time=0.006,aerodynamic_force_tolerance=5e-4,
            aerodynamic_circulation_tolerance=5e-4,maximum_static_iterations=80,
            static_relaxation=0.5,geometry_tolerance=2e-6,load_tolerance=5e-4,
            progress_frequency=10000)
        options=(;U_min=20.,U_max=22.,dU=2.,plateau_time=0.006,transient_time=0.0045,
            fixed_time_step=0.0005,angle_of_attack=0.01,chordwise_panels=2,spanwise_panels=4,
            wake_length_chords=2.,pulse_amplitude=0.01,pulse_start=0.,pulse_duration=0.004,
            on_unconverged=:continue,save_time_histories=true,save_boundary_states=true,
            keep_solver_result=true,plot_results=false,dynamic_options=(progress_frequency=10000,),equilibrium_options)
        checked=run_pazy_stepped_bifurcation(;options...,dry_run=true,output_directory=joinpath(directory,"dry_run"))
        @test !ispath(joinpath(directory,"dry_run"))
        @test length(checked.plan.plateaus) == 4
        sweep=run_pazy_stepped_bifurcation(;options...,output_directory=joinpath(directory,"sweep"))
        reference=run_pazy_stepped_bifurcation(;options...,U_max=20.,sweep_order=:upsweep,
            output_directory=joinpath(directory,"reference"))
        n=sweep.plan.steps_per_plateau
        @test sweep.solver_result.tip_out_of_plane[1:n+1] ≈ reference.solver_result.tip_out_of_plane rtol=1e-10 atol=1e-12
        @test getproperty.(sweep.boundaries,:accepted_step) == [12,24,36,48]
        @test getproperty.(sweep.boundaries,:aerodynamic_step) == [13,25,37,49]
        @test getproperty.(sweep.boundaries,:global_time_s) ≈ [0.006,0.012,0.018,0.024]
        @test getproperty.(sweep.boundaries,:aerodynamic_time_s) ≈ getproperty.(sweep.boundaries,:global_time_s)
        @test length(unique(getproperty.(sweep.boundaries,:structural_problem_id))) == 1
        @test length(unique(getproperty.(sweep.boundaries,:aerodynamic_problem_id))) == 1
        @test first(sweep.boundaries).maximum_section_velocity_norm > 1e-10
        @test length(sweep.boundary_states) == 4
        @test sweep.boundary_states[1].structural_state ≈ reference.boundary_states[1].structural_state rtol=1e-10 atol=1e-11
        @test sweep.solver_result.aerodynamic_problem.state.stepIndex == 49
        @test length(sweep.rows) == 16
        @test count(r -> r.direction == "upsweep",sweep.rows) == 8
        @test count(r -> r.direction == "downsweep",sweep.rows) == 8
        @test !isfile(joinpath(sweep.output_directory,"failure.txt"))
        @test_throws ArgumentError run_pazy_stepped_bifurcation(;options...,output_directory=sweep.output_directory)
        for plateau in sweep.plan.plateaus
            prefix=joinpath(sweep.output_directory,"plateau_$(lpad(plateau.index,4,'0'))")
            history=readdlm(prefix*"_history.csv",',',Float64;skipstart=1)
            @test all(history[:,3] .== plateau.airspeed_mps)
            @test size(history,1) == n
            @test first(history[:,1]) > plateau.start_time_s
            @test last(history[:,1]) ≈ plateau.end_time_s
            @test all(diff(history[:,1]) .> 0)
            rows=filter(r -> r.plateau_index == plateau.index,sweep.rows)
            selected=history[:,2] .>= sweep.plan.observation_start
            for (j,row) in enumerate(rows)
                @test row.sample_minimum == minimum(history[selected,j+3])
                @test row.sample_maximum == maximum(history[selected,j+3])
                @test row.window_start_s >= sweep.plan.observation_start
                @test row.window_end_s ≈ options.plateau_time
            end
        end
        # With :error the completed plateau must survive before the run stops.
        @test_throws ErrorException run_pazy_stepped_bifurcation(;options...,check_convergence=false,
            on_unconverged=:error,output_directory=joinpath(directory,"stopped"))
        stopped=readdlm(joinpath(directory,"stopped","bifurcation_summary.csv"),',';skipstart=1)
        @test size(stopped,1) == 4
        @test isfile(joinpath(directory,"stopped","failure.txt"))
        include(joinpath(@__DIR__,"PazySteppedBifurcationPlots.jl"))
        plot_paths=PazySteppedBifurcationPlots.plot_stepped_bifurcation(sweep.output_directory;
            save_panels=false,formats=("png","pdf"),figure_size=(1100,860),dpi=150)
        @test length(plot_paths) == 3
        @test all(path -> isfile(path) && filesize(path) > 1000,plot_paths)
        @test all(isfile,values(sweep.branch_summary_paths))
        println("Short physical plumbing validation (not settled LCO results): ",directory)
        println("Preview: ",first(plot_paths))
    end
end
