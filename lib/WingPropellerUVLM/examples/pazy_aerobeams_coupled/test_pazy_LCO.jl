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
