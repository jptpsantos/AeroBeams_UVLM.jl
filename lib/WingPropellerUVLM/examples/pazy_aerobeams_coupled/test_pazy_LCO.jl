# Short physical regression; all generated CSVs and plots use a temporary folder.
# julia --startup-file=no --compiled-modules=existing --project=. <this file>
using Test
previous_run = get(ENV, "PAZY_LCO_RUN", nothing)
ENV["PAZY_LCO_RUN"] = "false"
try
    include(joinpath(@__DIR__, "run_pazy_LCO.jl"))
finally
    isnothing(previous_run) ? delete!(ENV, "PAZY_LCO_RUN") : (ENV["PAZY_LCO_RUN"] = previous_run)
end

@testset "Independent linear sweeps from coupled equilibria" begin
    _, _, chord, _ = AeroBeams.geometrical_properties_Pazy()
    nominal_dt = 0.25 * chord / 60.0
    fractional_duration = (60.0-35.0)/0.85
    aligned_dt = pazy_lco_aligned_timestep(40.0, fractional_duration, nominal_dt)
    @test nominal_dt / 2 <= aligned_dt <= nominal_dt
    @test round(Int, 40.0/aligned_dt) * aligned_dt ≈ 40.0
    @test round(Int, fractional_duration/aligned_dt) * aligned_dt ≈ fractional_duration
    @test pazy_lco_aligned_timestep(40.0, 20.0, nominal_dt) ==
        20.0 / ceil(Int, 20.0/nominal_dt)
    @test_throws ArgumentError pazy_lco_aligned_timestep(40.0, sqrt(2.0), nominal_dt)
    @test pazy_linear_airspeed(0.; start_speed=35., end_speed=52., duration=40.) == 35.
    @test pazy_linear_airspeed(20.; start_speed=35., end_speed=52., duration=40.) == 43.5
    @test pazy_linear_airspeed(40.; start_speed=35., end_speed=52., duration=40.) == 52.
    @test pazy_linear_airspeed(0.; start_speed=52., end_speed=35., duration=20.) == 52.
    @test pazy_linear_airspeed(10.; start_speed=52., end_speed=35., duration=20.) == 43.5
    @test pazy_linear_airspeed(20.; start_speed=52., end_speed=35., duration=20.) == 35.
    @test pazy_linear_airspeed(1.; start_speed=35., end_speed=52., duration=40.) - 35 ≈ 0.425
    @test pazy_linear_airspeed(1.; start_speed=52., end_speed=35., duration=20.) - 52 ≈ -0.85
    equilibrium_options = (; maximum_aerodynamic_time=0.08, minimum_aerodynamic_time=0.02,
        aerodynamic_window_time=0.006, aerodynamic_force_tolerance=5e-4,
        aerodynamic_circulation_tolerance=5e-4, maximum_static_iterations=80,
        static_relaxation=0.5, geometry_tolerance=2e-6, load_tolerance=5e-4)
    case_options = (; initial_speed=20.0, peak_speed=22.0,
        upsweep_duration=0.006, downsweep_duration=0.003, angle_of_attack=0.01,
        chordwise_panels=2, spanwise_panels=4, wake_length_chords=2.0,
        progress_frequency=1000)
    mktempdir() do directory
        checked = run_pazy_linear_LCO(; case_options..., equilibrium_options,
            output_directory=directory)
        up, down = checked.upsweep, checked.downsweep
        @test up.result.time[1] == down.result.time[1] == 0.0
        @test last(up.result.time) ≈ 0.006
        @test last(down.result.time) ≈ 0.003
        @test up.result.airspeed_history[1] == down.result.airspeed_history[end] == 20.0
        @test down.result.airspeed_history[1] == up.result.airspeed_history[end] == 22.0
        @test all(diff(up.result.airspeed_history) .> 0)
        @test all(diff(down.result.airspeed_history) .< 0)
        @test up.result.structural_problem.model !== down.result.structural_problem.model
        @test up.equilibrium.structural_problem.model !== down.equilibrium.structural_problem.model
        @test !isapprox(down.result.initial_equilibrium_state.structural_state,
            up.result.structural_problem.x; rtol=1e-5, atol=1e-9)
        for branch in (up, down)
            result, equilibrium = branch.result, branch.equilibrium
            initial = result.initial_equilibrium_state
            @test equilibrium.converged
            @test initial.structural_state == equilibrium.structural_problem.x
            @test initial.maximum_velocity <= 1e-10
            @test result.structural_problem.model !== equilibrium.structural_problem.model
            @test initial.aerodynamic_snapshot.Γ == equilibrium.aerodynamic_problem.state.system.Γ
            @test initial.aerodynamic_snapshot.nwake == equilibrium.aerodynamic_problem.state.system.nwake
            @test all(initial.aerodynamic_snapshot.nwake .> 0)
            @test UVLM.top_left(initial.aerodynamic_snapshot.wakes[1][1,1]) ==
                UVLM.top_left(equilibrium.aerodynamic_problem.state.system.wakes[1][1,1])
            @test equilibrium.dt ≈ result.dt
            @test result.aerodynamic_problem.state.timeNow ≈ last(result.time)
            @test result.aerodynamic_problem.state.stepIndex == length(result.time)
            @test all(isapprox.(diff(result.time), result.dt; rtol=1e-12))
            root = result.structural_problem.model.elements[1].compStates
            fiber_tangent = [1.,0.,0.] + root.γ + cross(root.κ, [0.,0.,2.25e-3/2])
            @test last(branch.root_strain_ue) ≈ 1e6 * (fiber_tangent[1] - 1) atol=1e-8
            @test branch.root_strain_ue[1] < 0
            @test result.tip_twist_degrees[1] ≈ wingtip_twist_degrees(equilibrium.structural_problem.model)
            @test result.tip_out_of_plane[1] ≈ -equilibrium.structural_problem.model.elements[end].nodalStates.u_n2[1]
            @test result.root_axial_strain[1] ≈ equilibrium.structural_problem.model.elements[1].compStates.γ[1]
            history = readdlm(joinpath(directory, "pazy_linear_LCO_$(branch.name).csv"), ',', Float64; skipstart=1)
            @test history[:,1] == result.time
            @test history[:,2] == result.airspeed_history
            @test history[:,3] == branch.root_strain_ue
            @test isfile(joinpath(directory, "pazy_linear_LCO_$(branch.name)_initialization.csv"))
        end
        for name in ("pazy_linear_LCO_root_strain_vs_airspeed.png", "pazy_linear_LCO_root_strain_history.png")
            @test filesize(joinpath(directory, name)) > 1000
        end
        @test !isfile(joinpath(directory, "pazy_linear_LCO_history.csv"))
        if "--keep-output" in ARGS
            preview_directory = mktempdir(; cleanup=false)
            for name in readdir(directory)
                cp(joinpath(directory, name), joinpath(preview_directory, name))
            end
            println("Short independent-sweep validation outputs (not LCO results): $preview_directory")
        end
    end
    mktempdir() do directory
        higher = run_pazy_linear_LCO(; case_options..., downsweep_initial_speed=30.0,
            equilibrium_options, output_directory=directory)
        up, down = higher.upsweep, higher.downsweep
        @test first(up.result.airspeed_history) == 20.0
        @test last(up.result.airspeed_history) == 22.0
        @test first(down.result.airspeed_history) == 30.0
        @test last(down.result.airspeed_history) == 20.0
        @test down.equilibrium.airspeed == 30.0
        @test down.result.dt == up.result.dt == higher.dt
        @test last(up.result.time) ≈ 0.006
        @test last(down.result.time) ≈ 0.003
        _, _, chord, _ = AeroBeams.geometrical_properties_Pazy()
        @test higher.dt * 30.0 / chord <= 0.25
        @test down.equilibrium.dt ≈ higher.dt
        @test all(diff(down.result.airspeed_history) .< 0)
        settings = read(joinpath(directory, "pazy_linear_LCO_settings.txt"), String)
        @test occursin("Downsweep: 30.0 -> 20.0", settings)
        @test occursin("Timestep reference airspeed [m/s]: 30.0", settings)
    end
    mktempdir() do directory
        pulse_options = merge(case_options, (; fixed_time_step=0.00025,
            maximum_wake_rows=44, core_radius=0.0008, wake_shedding_fraction=0.2,
            aerodynamic_interaction=false, symmetric_wing=false))
        quiet = run_pazy_linear_LCO(; pulse_options..., equilibrium_options,
            output_directory=joinpath(directory, "quiet"))
        forced = run_pazy_linear_LCO(; pulse_options..., equilibrium_options,
            downsweep_pulse_amplitude=0.1, downsweep_pulse_start=0.00025,
            downsweep_pulse_duration=0.002, plot_branches=:downsweep,
            output_directory=joinpath(directory, "forced"))
        @test all(iszero, quiet.downsweep.tip_force_A1)
        @test all(iszero, forced.upsweep.tip_force_A1)
        @test isapprox(forced.upsweep.result.tip_out_of_plane,
            quiet.upsweep.result.tip_out_of_plane; rtol=1e-10, atol=1e-12)
        @test isapprox(forced.upsweep.root_strain_ue,
            quiet.upsweep.root_strain_ue; rtol=1e-10, atol=1e-9)
        down = forced.downsweep
        @test minimum(down.tip_force_A1) < 0 < maximum(down.tip_force_A1)
        @test isapprox(sum(down.tip_force_A1)*forced.dt, 0; atol=1e-15)
        @test all(iszero, down.tip_force_A1[down.result.time .<= 0.00025])
        @test all(iszero, down.tip_force_A1[down.result.time .>= 0.00225])
        @test down.result.airspeed_history == quiet.downsweep.result.airspeed_history
        @test down.root_strain_ue[1] == quiet.downsweep.root_strain_ue[1]
        @test maximum(abs.(down.result.tip_out_of_plane .-
            quiet.downsweep.result.tip_out_of_plane)) > 1e-10
        data, header = readdlm(joinpath(directory,"forced",
            "pazy_linear_LCO_downsweep.csv"), ',', Float64; header=true)
        force_column = findfirst(==("tip_force_A1_N"), vec(header))
        @test !isnothing(force_column)
        @test data[:,force_column] == down.tip_force_A1
        settings = read(joinpath(directory,"forced","pazy_linear_LCO_settings.txt"),String)
        @test occursin("Downsweep tip pulse coefficient [N]: 0.1",settings)
        @test occursin("Upsweep tip pulse coefficient [N]: 0.0",settings)
        for branch in (forced.upsweep, forced.downsweep)
            dynamic = branch.result.aerodynamic_problem
            @test dynamic.aeroSolver == branch.equilibrium.aerodynamic_problem.aeroSolver
            @test dynamic.aeroSolver.coreRadius == 0.0008
            @test dynamic.aeroSolver.wakeSheddingFraction == 0.2
            @test !dynamic.aeroSolver.interaction
            @test dynamic.model.symmetric == [false]
            @test size(dynamic.state.system.wakes[1],1) == 44
        end
        @test isfile(joinpath(directory,"forced","pazy_linear_LCO_upsweep.csv"))
        @test isfile(joinpath(directory,"forced","pazy_linear_LCO_downsweep.csv"))
        @test isfile(joinpath(directory,"forced","pazy_linear_LCO_root_strain_history_downsweep.png"))
        @test !isfile(joinpath(directory,"forced","pazy_linear_LCO_root_strain_vs_airspeed.png"))
        combined_path = joinpath(directory,"quiet","pazy_linear_LCO_root_strain_vs_airspeed.png")
        combined_bytes = read(combined_path)
        up_figure = plot_pazy_LCO(; output_directory=joinpath(directory,"quiet"),
            plot_branches=:upsweep)
        @test only(up_figure.series_list)[:label] == "Upsweep"
        @test only(up_figure.series_list)[:x] == quiet.upsweep.result.airspeed_history
        @test only(up_figure.series_list)[:y] == quiet.upsweep.root_strain_ue
        down_figure = plot_pazy_LCO(; output_directory=joinpath(directory,"quiet"),
            plot_branches=:downsweep)
        @test only(down_figure.series_list)[:label] == "Downsweep"
        @test only(down_figure.series_list)[:y] == quiet.downsweep.root_strain_ue
        @test read(combined_path) == combined_bytes
        @test isfile(joinpath(directory,"quiet","pazy_linear_LCO_root_strain_vs_airspeed_upsweep.png"))
        @test isfile(joinpath(directory,"quiet","pazy_linear_LCO_root_strain_vs_airspeed_downsweep.png"))
        both_forced = run_pazy_linear_LCO(; pulse_options..., equilibrium_options,
            upsweep_pulse_amplitude=0.1, upsweep_pulse_start=0.0005,
            upsweep_pulse_duration=0.002,
            downsweep_pulse_amplitude=0.1, downsweep_pulse_start=0.00025,
            downsweep_pulse_duration=0.002,
            output_directory=joinpath(directory,"both_forced"))
        up = both_forced.upsweep
        @test minimum(up.tip_force_A1) < 0 < maximum(up.tip_force_A1)
        @test isapprox(sum(up.tip_force_A1)*both_forced.dt, 0; atol=1e-15)
        @test all(iszero, up.tip_force_A1[up.result.time .<= 0.0005])
        @test all(iszero, up.tip_force_A1[up.result.time .>= 0.0025])
        @test up.result.airspeed_history == quiet.upsweep.result.airspeed_history
        @test up.root_strain_ue[1] == quiet.upsweep.root_strain_ue[1]
        @test maximum(abs.(up.result.tip_out_of_plane .-
            quiet.upsweep.result.tip_out_of_plane)) > 1e-10
        @test isapprox(both_forced.downsweep.root_strain_ue, forced.downsweep.root_strain_ue;
            rtol=1e-10, atol=1e-9)
        @test both_forced.downsweep.tip_force_A1 == forced.downsweep.tip_force_A1
        up_data, up_header = readdlm(joinpath(directory,"both_forced",
            "pazy_linear_LCO_upsweep.csv"), ',', Float64; header=true)
        @test up_data[:,findfirst(==("tip_force_A1_N"),vec(up_header))] == up.tip_force_A1
        @test occursin("Upsweep tip pulse coefficient [N]: 0.1; start [local s]: 0.0005",
            read(joinpath(directory,"both_forced","pazy_linear_LCO_settings.txt"),String))
        up_only = run_pazy_linear_LCO(; pulse_options..., equilibrium_options,
            run_branches=:upsweep, downsweep_pulse_amplitude=0.1,
            downsweep_pulse_start=1000.0,
            output_directory=joinpath(directory,"up_only"))
        @test isnothing(up_only.downsweep)
        @test up_only.run_branches == up_only.plot_branches == :upsweep
        @test up_only.dt == quiet.dt
        @test isapprox(up_only.upsweep.root_strain_ue, quiet.upsweep.root_strain_ue;
            rtol=1e-10, atol=1e-9)
        @test all(iszero, up_only.upsweep.tip_force_A1)
        @test !isfile(joinpath(directory,"up_only","pazy_linear_LCO_downsweep_initialization.csv"))
        @test !isfile(joinpath(directory,"up_only","pazy_linear_LCO_downsweep.csv"))
        @test isfile(joinpath(directory,"up_only","pazy_linear_LCO_root_strain_vs_airspeed_upsweep.png"))
        down_only = run_pazy_linear_LCO(; pulse_options..., equilibrium_options,
            run_branches=:downsweep, downsweep_pulse_amplitude=0.1,
            upsweep_pulse_amplitude=0.1, upsweep_pulse_start=1000.0,
            downsweep_pulse_start=0.00025, downsweep_pulse_duration=0.002,
            output_directory=joinpath(directory,"down_only"))
        @test isnothing(down_only.upsweep)
        @test down_only.run_branches == down_only.plot_branches == :downsweep
        @test down_only.dt == forced.dt
        @test isapprox(down_only.downsweep.root_strain_ue, forced.downsweep.root_strain_ue;
            rtol=1e-10, atol=1e-9)
        @test down_only.downsweep.tip_force_A1 == forced.downsweep.tip_force_A1
        @test !isfile(joinpath(directory,"down_only","pazy_linear_LCO_upsweep_initialization.csv"))
        @test !isfile(joinpath(directory,"down_only","pazy_linear_LCO_upsweep.csv"))
        @test isfile(joinpath(directory,"down_only","pazy_linear_LCO_root_strain_vs_airspeed_downsweep.png"))
        @test_throws ArgumentError run_pazy_linear_LCO(; pulse_options...,
            downsweep_pulse_amplitude=0.1, downsweep_pulse_start=0.002,
            downsweep_pulse_duration=0.002)
        @test_throws ArgumentError run_pazy_linear_LCO(; pulse_options...,
            downsweep_pulse_amplitude=0.1, downsweep_pulse_start=0.0,
            downsweep_pulse_duration=0.001)
        @test_throws ArgumentError run_pazy_linear_LCO(; pulse_options...,
            upsweep_pulse_amplitude=0.1, upsweep_pulse_start=0.005,
            upsweep_pulse_duration=0.002)
        @test_throws ArgumentError run_pazy_linear_LCO(; pulse_options...,
            upsweep_pulse_amplitude=0.1, upsweep_pulse_start=0.0,
            upsweep_pulse_duration=0.001)
    end
    @test_throws AssertionError run_pazy_linear_LCO(; case_options..., downsweep_initial_speed=19.0)
    @test_throws ArgumentError run_pazy_linear_LCO(; case_options..., core_radius=0.0)
    @test_throws ArgumentError run_pazy_linear_LCO(; case_options..., wake_shedding_fraction=1.1)
    @test_throws ArgumentError run_pazy_linear_LCO(; case_options..., maximum_wake_rows=0)
    @test_throws ArgumentError run_pazy_linear_LCO(; case_options..., chordwise_panels=0)
    @test_throws ArgumentError run_pazy_linear_LCO(; case_options..., plot_branches=:invalid)
    @test_throws ArgumentError run_pazy_linear_LCO(; case_options..., run_branches=:invalid)
    @test_throws ArgumentError run_pazy_linear_LCO(; case_options...,
        run_branches=:downsweep, plot_branches=:upsweep)
    @test_throws ArgumentError plot_pazy_LCO(; plot_branches=:invalid)
    mktempdir() do directory
        @test_throws ErrorException run_pazy_linear_LCO(; case_options...,
            equilibrium_options=merge(equilibrium_options, (;maximum_static_iterations=1)),
            output_directory=directory)
        @test isfile(joinpath(directory, "pazy_linear_LCO_upsweep_initialization.csv"))
        @test !isfile(joinpath(directory, "pazy_linear_LCO_upsweep.csv"))
        @test !isfile(joinpath(directory, "pazy_linear_LCO_downsweep.csv"))
    end
end
@testset "Continuous Pazy speed sweep" begin
    mktempdir() do directory
        sweep = run_pazy_LCO(; initial_speed=30.0, final_speed=33.0, speed_step=2.0,
            duration_per_speed=0.006, lco_window_duration=0.002,
            chordwise_panels=3, spanwise_panels=6, wake_length_chords=5.0,
            settling_time=0.0008, perturbation_amplitude=0.05,
            perturbation_duration=0.002, output_directory=directory)
        result = sweep.result
        n = (length(result.time) - 1) ÷ 3
        @test sweep.speeds == [30.0, 32.0, 33.0] # includes non-divisible endpoint
        @test result.airspeed_history[1:n+1] == fill(30.0, n+1)
        @test result.airspeed_history[n+2:2n+1] == fill(32.0, n)
        @test result.airspeed_history[2n+2:end] == fill(33.0, n)
        @test all(isapprox.(diff(result.time), result.dt; rtol=1e-12))
        @test result.time[end] ≈ 0.018
        @test sweep.summary[:, 3] ≈ [0.006, 0.012, 0.018]
        @test result.aerodynamic_problem.state.stepIndex == length(result.time)
        @test result.aerodynamic_problem.operatingPoint.airspeed == 33.0
        @test result.aerodynamic_problem.state.system.reference[].V == 33.0
        @test result.aerodynamic_problem.state.system.nwake[1] ==
            min(length(result.time)-1, size(result.aerodynamic_problem.state.system.wakes[1], 1))
        @test result.aerodynamic_problem.state.system.nwake[1] > n
        @test isempty(result.structural_problem.nodalStatesOverTime)
        @test sweep.aoa_histories[1][end] == sweep.aoa_histories[2][1]
        @test sweep.aoa_histories[2][end] == sweep.aoa_histories[3][1]

        # Independently run the first plateau from rest with the legacy constant-
        # speed API and the same dt. Its response must match the sweep's prefix.
        _, _, chord, _ = AeroBeams.geometrical_properties_Pazy()
        base_options = (; density=1.225, angle_of_attack=deg2rad(3.0), sideslip=0.0,
            initial_airspeed_fraction=1.0, airspeed_ramp_duration=0.0,
            duration=0.006, settling_time=0.0008,
            chordwise_panels=3, spanwise_panels=6, symmetric_wing=true,
            maximum_wake_rows=size(result.aerodynamic_problem.state.system.wakes[1], 1),
            core_radius=1e-3, wake_shedding_fraction=0.1, aerodynamic_interaction=true,
            save_uvlm_history=false, uvlm_save_frequency=1,
            newton_maximum_iterations=20, newton_absolute_tolerance=1e-5,
            newton_relative_tolerance=1e-5, newton_display_iterations=false,
            newton_always_update_jacobian=true,
            perturbation_amplitude=0.05, perturbation_duration=0.002,
            animation_frames=2, progress_frequency=1000,
            time_step_chords=result.dt * 30 / chord,
            save_animation_history=false, save_structural_history=false,
            save_aerodynamic_load_history=false, coupling_scheme=:strong,
            coupling_maximum_iterations=20, coupling_relaxation=0.5,
            coupling_geometry_tolerance=1e-4, coupling_load_tolerance=1e-4)
        constant = run_pazy_wing_uvlm(; base_options..., airspeed=30.0)
        @test constant.tip_twist_degrees ≈ result.tip_twist_degrees[1:n+1] rtol=1e-9 atol=1e-11
        @test constant.tip_out_of_plane ≈ result.tip_out_of_plane[1:n+1] rtol=1e-9 atol=1e-11
        prescribed = run_pazy_wing_uvlm(; base_options..., airspeed=30.0,
            airspeed_profile=t -> 30.0)
        @test prescribed.structural_problem.x ≈ constant.structural_problem.x rtol=1e-10 atol=1e-12
        @test prescribed.tip_twist_degrees ≈ constant.tip_twist_degrees rtol=1e-10 atol=1e-12

        history = readdlm(joinpath(directory, "pazy_LCO_continuous_history.csv"), ',', Float64; skipstart=1)
        @test history[:, 1] == result.time
        @test history[:, 2] == result.airspeed_history
        @test history[:, 5] ≈ rad2deg(3pi/180) .+ result.tip_twist_degrees
        for speed in sweep.speeds
            @test isfile(joinpath(directory, "pazy_LCO_V$(speed).csv"))
        end
        for name in ("pazy_LCO_wingtip_aoa_histories.png", "pazy_LCO_wingtip_aoa_amplitude.png")
            @test filesize(joinpath(directory, name)) > 1000
        end
        if "--keep-output" in ARGS
            preview_directory = mktempdir(; cleanup=false)
            for name in readdir(directory)
                cp(joinpath(directory, name), joinpath(preview_directory, name))
            end
            println("Short validation outputs (not a settled LCO): $preview_directory")
        end
        @test_throws AssertionError run_pazy_wing_uvlm(; base_options..., airspeed=30.0,
            airspeed_profile=t -> NaN)
    end
end
