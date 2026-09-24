using Test
include(joinpath(@__DIR__, "PazyWingUVLMCoupling.jl"))

@testset "Pazy clamped-root geometry and UVLM motion" begin
    # Minimal beam outputs, with a deliberately small recovery error at the clamp.
    R0 = [0.0 0.0 -1.0; 0.0 1.0 0.0; 1.0 0.0 0.0]
    root_error = [0.0, 0.0, 1e-14]
    node_states = (u_n1=root_error, p_n1=zeros(3),
        u_n2=zeros(3), p_n2=zeros(3), p_n2_b=zeros(3))
    element = (r_n1=zeros(3), r_n2=[0.0, 0.0, 0.55],
        R0_n1=R0, R0_n2=R0, nodalStates=node_states)
    model = (elements=[element],)
    weights = hcat(1 .- collect(0:15)./15, collect(0:15)./15)
    grid, positions = wing_geometry(model, 0.0989, 0.44096, 4, weights)
    @test all(iszero, grid[2, :, 1])
    @test positions[:, 1] == zeros(3)
    @test root_error[3] == 1e-14 # The transfer must not mutate beam outputs.
    @test wingtip_twist_degrees(model) == 0.0

    zero_loads = zeros(6, 2)
    @test interface_load_residual(zero_loads, zero_loads, 0.1) == 0.0
    changed_loads = copy(zero_loads)
    changed_loads[1, 2] = 1.0
    @test interface_load_residual(changed_loads, zero_loads, 0.1) == 1.0

    for symmetric in (true, false)
        aero_model = UVLM.create_UVLMModel(surfaces=[grid], symmetric=symmetric,
            reference=UVLM.Reference(0.55*0.0989, 0.0989, 0.55, zeros(3), 40.0, 1.225))
        p = UVLM.create_UVLMDynamicProblem(model=aero_model,
            operatingPoint=UVLM.create_OperatingPoint(airspeed=40.0, angleOfAttack=deg2rad(3.0)),
            timeVector=collect(0:2).*(0.0989/(4*40)), maximumWakeRows=2)
        UVLM.begin_time_step!(p)
        loads = UVLM.evaluate_trial!(p; surfaces=[grid])
        @test all(f -> all(isfinite, f), loads.forces[1])
        UVLM.commit_time_step!(p)

        # A moving trial must give nonzero surface velocity without advancing
        # the accepted wake. Only commit_time_step! advances the wake.
        moving_grid = copy(grid)
        moving_grid[3, :, :] .+= reshape(collect(0:15).*1e-7, 1, 16)
        UVLM.begin_time_step!(p)
        loads = UVLM.evaluate_trial!(p; surfaces=[moving_grid])
        @test all(f -> all(isfinite, f), loads.forces[1])
        @test any(v -> norm(v) > 0, p.workspace.system.Vcp[1])
        @test p.state.system.nwake == [1]
        UVLM.evaluate_trial!(p; surfaces=[moving_grid])
        @test p.state.system.nwake == [1]
        UVLM.commit_time_step!(p)
        @test p.state.system.nwake == [2]
    end
end

# Optional, slower integration check: ramp, tip pulse, and subsequent response.
# This compressed startup is a numerical test, not a flutter benchmark.
if "--dynamic" in ARGS
    @testset "Pazy coupled startup and pulse" begin
        result = run_pazy_wing_uvlm(;
            airspeed=40.0, density=1.225, angle_of_attack=deg2rad(3.0), sideslip=0.0,
            initial_airspeed_fraction=0.05, airspeed_ramp_duration=0.01,
            duration=0.06, settling_time=0.04,
            chordwise_panels=4, spanwise_panels=15, symmetric_wing=true,
            maximum_wake_rows=80, core_radius=1e-3, wake_shedding_fraction=0.1,
            aerodynamic_interaction=true, save_uvlm_history=true, uvlm_save_frequency=1,
            newton_maximum_iterations=50, newton_absolute_tolerance=1e-8,
            newton_relative_tolerance=1e-8, newton_display_iterations=false,
            newton_always_update_jacobian=false, perturbation_amplitude=0.05,
            perturbation_duration=0.01, animation_frames=150, progress_frequency=100)
        @test result.aerodynamic_problem.results.completed
        @test all(isfinite, result.tip_out_of_plane)
        @test all(isfinite, result.tip_twist_degrees)
        @test length(result.tip_twist_degrees) == length(result.time)
        @test maximum(abs, result.tip_out_of_plane) > 1e-6
        @test size(result.aerodynamic_nodal_load_history) ==
            (6, 16, length(result.time))
        @test all(isfinite, result.aerodynamic_nodal_load_history)
        @test maximum(abs, result.aerodynamic_nodal_load_history) > 0
        @test all(result.coupling_iterations[2:end] .>= 2)
        @test all(result.coupling_geometry_residual[2:end] .<= 1e-5)
        @test all(result.coupling_load_residual[2:end] .<= 1e-3)
        @test length(result.structural_problem.savedTimeVector) == length(result.time)
        @test all(f -> all(isfinite, f), result.aerodynamic_problem.results.forceOverTime)
    end
end
