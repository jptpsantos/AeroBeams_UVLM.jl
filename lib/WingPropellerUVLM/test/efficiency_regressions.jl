using Test, WingPropellerUVLM, LinearAlgebra, StaticArrays
import WingPropellerUVLM as EU

function efficiency_grid(; shift=0.0, warp=0.0)
    grid = zeros(3,4,5)
    for j in 1:5,i in 1:4
        grid[:,i,j] .= ((i-1)/3,(j-1)/4+shift,warp*(i-1)*(j-1)/12)
    end
    return grid
end

@testset "Full AIC already contains trailing-edge contributions" begin
    for warp in (0.0,0.04), symmetric in (false,true), shared_id in (false,true),
        trailing in (false,true), interaction in (false,true)
        grids = [efficiency_grid(;warp),efficiency_grid(;shift=1.0,warp)]
        panels = [grid_to_surface_panels(grid;fcore=(c,ds)->0.001)[3] for grid in grids]
        shedding = [[EU.bottom_left(panel[end,j])+SVector(0.05,0,0.01)
            for j in axes(panel,2)] for panel in panels]
        for i in eachindex(panels)
            push!(shedding[i],EU.bottom_right(panels[i][end,end])+SVector(0.05,0,0.01))
        end
        N = sum(length,panels)
        A = zeros(N,N)
        kwargs = (;symmetric=fill(symmetric,2),surface_id=shared_id ? [1,1] : [1,2],
            wake_shedding_locations=shedding,trailing_vortices=fill(trailing,2),
            interaction_id=[1,2],interaction)
        EU.influence_coefficients!(A,panels;kwargs...)
        updated = copy(A)
        EU.update_trailing_edge_coefficients!(updated,panels;kwargs...)
        @test updated ≈ A atol=1e-12 rtol=0
        if !interaction
            @test iszero(A[1:length(panels[1]),length(panels[1])+1:end])
        end
    end
end

function efficiency_aero()
    grid = efficiency_grid(;warp=0.02)
    model = create_UVLMModel(surfaces=[grid],reference=Reference(1.0,1.0,1.0,zeros(3),10.0))
    point = create_OperatingPoint(airspeed=10.0,angleOfAttack=0.04)
    create_UVLMDynamicProblem(model=model,operatingPoint=point,
        timeVector=[0.0,0.01,0.02,0.03],maximumWakeRows=2,trackingTimeSteps=false)
end

@testset "Fixed geometry updates only trailing-edge AIC columns" begin
    p = efficiency_aero()
    begin_time_step!(p)
    evaluate_trial!(p)
    commit_time_step!(p)
    for moved in (false,true)
        full = deepcopy(p.state.system)
        cached = deepcopy(full)
        if moved
            grid = copy(p.model.grids[1])
            grid[3,:,:] .+= 0.001
            for sys in (full,cached)
                EU.apply_trial_geometry!(sys,[grid],p.aeroSolver)
            end
        end
        for (sys,recompute) in ((full,true),(cached,false))
            EU.propagate_system!(sys,sys.freestream[],0.01;
                additional_velocity=nothing,repeated_points=nothing,nwake=sys.nwake,
                eta=0.1,calculate_influence_matrix=recompute,near_field_analysis=true,
                derivatives=false,advance_wake=false)
        end
        @test cached.AIC ≈ full.AIC atol=1e-12 rtol=0
        @test cached.Γ ≈ full.Γ atol=1e-12 rtol=0
        @test cached.chord_seg_forces ≈ full.chord_seg_forces
        @test cached.span_seg_forces ≈ full.span_seg_forces
    end
end

function efficiency_forces(sys; buffers=false, unsteady=true)
    kwargs = buffers ? (;
        chord_force_buffers=sys.chord_seg_forces,span_force_buffers=sys.span_seg_forces,
        unsteady_force_buffers=sys.unsteady_forces,span_velocity_buffers=sys.span_force_velocities) : NamedTuple()
    near_field_forces!(sys.properties,sys.surfaces,sys.wakes,sys.reference[],sys.freestream[],sys.Γ;
        dΓdt=unsteady ? sys.dΓdt : nothing,additional_velocity=nothing,Vh=sys.Vh,Vv=sys.Vv,
        symmetric=sys.symmetric,nwake=sys.nwake,surface_id=sys.surface_id,
        wake_finite_core=sys.wake_finite_core,wake_shedding_locations=sys.wake_shedding_locations,
        trailing_vortices=sys.trailing_vortices,xhat=sys.xhat[],kwargs...)
end

@testset "Reusable near-field buffers clear every trial" begin
    p = efficiency_aero()
    begin_time_step!(p)
    evaluate_trial!(p)
    sys = p.workspace.system
    for unsteady in (true,false,true,false)
        expected = deepcopy(efficiency_forces(sys;unsteady))
        for group in (sys.chord_seg_forces,sys.span_seg_forces,sys.unsteady_forces,sys.span_force_velocities), array in group
            fill!(array,SVector(8.0,9.0,10.0))
        end
        actual = efficiency_forces(sys;buffers=true,unsteady)
        @test actual[2] === sys.chord_seg_forces
        @test actual[3] === sys.span_seg_forces
        @test actual[4] === sys.unsteady_forces
        for i in 2:4
            @test actual[i] ≈ expected[i]
        end
        for (a,b) in zip(actual[1][1],expected[1][1]), field in fieldnames(typeof(a))
            @test getfield(a,field) ≈ getfield(b,field)
        end
        @test all(iszero,sys.span_seg_forces[1][end,:])
        !unsteady && @test all(iszero,sys.unsteady_forces[1])
    end
end

@testset "Trial buffers, rollback and failed commit isolation" begin
    p = efficiency_aero()
    accepted = p.state.system
    trial = p.workspace.system
    begin_time_step!(p)
    evaluate_trial!(p)
    first_gamma = copy(trial.Γ)
    @test p.workspace.system === trial
    moved = deepcopy(p.model.grids)
    moved[1][3,:,:] .+= 0.002
    evaluate_trial!(p;surfaces=moved)
    evaluate_trial!(p)
    @test p.workspace.system === trial
    @test trial.Γ ≈ first_gamma
    @test p.state.system === accepted
    @test iszero(accepted.Γ)
    @test all(trial_array !== accepted_array for
        field in (:grids,:surfaces,:wakes,:chord_seg_forces,:span_force_velocities)
        for (trial_array,accepted_array) in zip(getfield(trial,field),getfield(accepted,field)))
    trial.nwake[1] = 3 # Invalid capacity exercises the separate commit buffer.
    @test_throws ArgumentError commit_time_step!(p)
    @test p.state.system === accepted && iszero(accepted.Γ)
    @test p.workspace.stepOpen && p.workspace.trialReady
    trial.nwake .= accepted.nwake
    commit_time_step!(p)
    @test p.state.system !== trial
    @test p.state.system.Γ ≈ first_gamma
    @test p.state.system.nwake == [1]
    begin_time_step!(p)
    bad = deepcopy(p.model.grids)
    bad[1][1] = NaN
    @test_throws ArgumentError evaluate_trial!(p;surfaces=bad)
    before = copy(p.state.system.Γ)
    rollback_time_step!(p)
    @test p.state.timeNow == 0.01
    @test p.state.system.Γ == before
    @test p.workspace.system === trial
    @test solve!(p).results.completed
end

@testset "Cached linear effective systems" begin
    M = [2.0 0;0 3]
    C = [0.2 0.3;-0.3 0.4]
    K = [12.0 -1;-1 15]
    u,v,a = [0.01,-0.02],[0.03,0.01],[-0.01,0.02]
    f,fn = [0.2,-0.1],[0.1,0.1]
    dt = 0.01
    for method in (:newmark_beta,:generalized_alpha)
        params = structural_integration_parameters(method)
        cache = EU.structural_effective_system(method,M,C,K,dt,params)
        expected = loose_explicit_aeroelastic_step(method,M,C,K,u,v,a,f,fn,zeros(2),zeros(2),dt,params)
        actual = loose_explicit_aeroelastic_step(method,M,C,K,u,v,a,f,fn,zeros(2),zeros(2),dt,params;effective_system=cache)
        @test actual.displacement ≈ expected.displacement
        @test actual.velocity ≈ expected.velocity
        @test actual.acceleration ≈ expected.acceleration
        @test actual.converged == expected.converged
    end
end
