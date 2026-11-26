using Test
using LinearAlgebra
include(joinpath(@__DIR__,"..","PazyWingUVLMStatic.jl"))
include(joinpath(@__DIR__,"..","PazyLinearUVLM.jl"))

@testset "Pazy fixed-wake UVLM state-space block" begin
    options = (;airspeed=15.0,density=1.225,angle_of_attack=deg2rad(0.6),
        chordwise_panels=2,spanwise_panels=3,symmetric_wing=true,maximum_wake_rows=4,
        core_radius=1e-3,wake_shedding_fraction=0.1,aerodynamic_interaction=true,
        maximum_aerodynamic_time=0.2,minimum_aerodynamic_time=0.0,
        aerodynamic_window_time=0.01,aerodynamic_force_tolerance=5e-4,
        aerodynamic_circulation_tolerance=5e-4,maximum_static_iterations=10,
        static_relaxation=0.5,geometry_tolerance=1e-5,load_tolerance=5e-4,
        consecutive_equilibrium_iterations=2,newton_maximum_iterations=50,
        newton_absolute_tolerance=1e-8,newton_relative_tolerance=1e-8,
        newton_display_iterations=false,progress_frequency=1000,
        mesh=:nominal,withWingtipRod=false)
    reference = run_pazy_wing_uvlm_static(;options...)
    @test reference.converged
    aerodynamic = reference.aerodynamic_problem
    system = aerodynamic.state.system
    grid = system.grids[1]
    displacement = zeros(size(grid))
    displacement[1,:,:] .= 1.0
    velocity = zeros(size(grid))
    velocity[3,:,:] .= 1.0
    original_gamma = copy(system.Γ)
    original_wake = [[UVLM.circulation_strength(panel) for panel in wake]
        for wake in system.wakes]

    model = PazyLinearUVLM.linearize_uvlm_step(aerodynamic;
        displacement_directions=[[displacement]],
        velocity_directions=[[velocity]],
        finite_difference_step=1e-6)
    state_count = length(model.state_reference)
    force_count = length(model.force_reference)
    @test size(model.A) == (state_count,state_count)
    @test size(model.B) == (state_count,2)
    @test size(model.C) == (force_count,state_count)
    @test size(model.D) == (force_count,2)
    @test all(isfinite,model.A) && all(isfinite,model.B)
    @test all(isfinite,model.C) && all(isfinite,model.D)
    @test model.relative_state_drift <= 1e-3
    @test model.relative_force_drift <= 1e-3
    @test system.Γ == original_gamma
    @test [[UVLM.circulation_strength(panel) for panel in wake]
        for wake in system.wakes] == original_wake

    bound_count = length(system.Γ)
    perturbed_state = copy(model.state_reference)
    perturbed_state[bound_count+1:end] .+= collect(1:length(perturbed_state)-bound_count) ./ 10
    original_wake_geometry = deepcopy(system.wakes)
    workspace = deepcopy(system)
    zero_grid = [zeros(size(grid))]
    next_state,_ = PazyLinearUVLM.evaluate_step!(workspace,system,aerodynamic,
        model.time_step,perturbed_state,zero_grid,zero_grid)
    wake_geometry(wakes) = [[(panel.rtl,panel.rtr,panel.rbl,panel.rbr)
        for panel in wake] for wake in wakes]
    actual_geometry = wake_geometry(workspace.wakes)
    expected_geometry = wake_geometry(original_wake_geometry)
    geometry_error = 0.0
    for surface in eachindex(actual_geometry),panel in eachindex(actual_geometry[surface]),
        corner in 1:4
        geometry_error = max(geometry_error,
            norm(actual_geometry[surface][panel][corner] -
                 expected_geometry[surface][panel][corner]))
    end
    @test geometry_error <= 1e-14
    wake_state = next_state[bound_count+1:end]
    expected_wake = Float64[]
    panels = LinearIndices(size(system.surfaces[1]))
    for span_panel in axes(system.wakes[1],2)
        push!(expected_wake,workspace.Γ[panels[end,span_panel]])
    end
    for row in 1:system.nwake[1]-1,span_panel in axes(system.wakes[1],2)
        previous_index = bound_count + (row - 1) * size(system.wakes[1],2) + span_panel
        push!(expected_wake,perturbed_state[previous_index])
    end
    @test wake_state == expected_wake

    @test_throws DimensionMismatch PazyLinearUVLM.linearize_uvlm_step(aerodynamic;
        displacement_directions=[[displacement]],velocity_directions=[])
end