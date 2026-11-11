module ChangRuntimeOptimizationTests

using Test,StaticArrays,LinearAlgebra
import WingPropellerUVLM as W
const EXAMPLE = joinpath(@__DIR__,"..","examples","chang_linear_aeroelastic")
include(joinpath(EXAMPLE,"src","ChangAeroelastic.jl"))
include(joinpath(EXAMPLE,"src","ChangRuntimeOptimizations.jl"))
const C = ChangAeroelastic
const F = ChangRuntimeOptimizations

function wake_fixture()
    panels = map(((0.0,1.0,0.0,2),(1.0,2.0,0.0,3),(0.3,0.9,-0.7,2))) do spec
        lo,hi,shift,ns = spec
        grid = zeros(3,3,ns + 1)
        for j = 1:ns+1,i = 1:3
            y = lo + (hi - lo) * (j - 1) / ns
            x = shift + (i - 1) / 2
            grid[:,i,j] .= (x,y,0.02y + 0.01x * y)
        end
        W.grid_to_surface_panels(grid; fcore=(c,ds) -> 0.004)[3]
    end
    system = W.System(collect(panels); nw=[4,4,3])
    system.surfaces .= panels
    system.previous_surfaces .= panels
    system.reference[] = W.Reference(2.0,1.0,2.0,zeros(3),12.0)
    system.freestream[] = W.Freestream(12.0,0.04,0.01,SVector(0.01,-0.02,0.03))
    system.Γ .= range(0.2,0.8; length=length(system.Γ))
    for isurf in eachindex(system.surfaces)
        surface = system.surfaces[isurf]
        locations = vcat([W.bottom_left(surface[end,j]) for j in axes(surface,2)],
            [W.bottom_right(surface[end,end])])
        system.wake_shedding_locations[isurf] .= locations
        system.Vte[isurf] .= Ref(SVector(0.002,-0.001,0.003))
        for j in axes(system.wakes[isurf],2),i in axes(system.wakes[isurf],1)
            a,b = (i - 1) * SVector(0.09,0.004,0.015),i * SVector(0.09,0.004,0.015)
            left,right = W.bottom_left(surface[end,j]),W.bottom_right(surface[end,j])
            system.wakes[isurf][i,j] = W.WakePanel(
                left + a,right + a,left + b,right + b,0.004,0.4 + 0.01i + 0.02j)
        end
        fill!(system.V[isurf],SVector(-9.0,-8.0,-7.0))
    end
    return system
end

function reference_velocities!(system; kwargs...)
    W.get_wake_velocities!(system.V,system.surfaces,system.wakes,
        system.reference[],system.freestream[],system.Γ,
        kwargs[:additional_velocity],system.Vte,system.symmetric,
        kwargs[:repeated_points],kwargs[:nwake],system.surface_id,
        system.wake_finite_core,system.wake_shedding_locations,
        system.trailing_vortices,system.xhat[];
        interaction_id=kwargs[:interaction_id],interaction=kwargs[:interaction])
end

@testset "Serial and threaded wake vertices preserve the backend result" begin
    fixture = wake_fixture()
    repeated = W.repeated_trailing_edge_points(fixture.surfaces)
    @test !isempty(repeated) # Exercise shared-edge vertex exclusion.
    additional = r -> SVector(0.01r[2],0.02r[3],-0.01r[1])
    for rows in ([0,0,0],[1,1,1],[2,2,1],[4,4,3],[4,1,3],[1,4,3]),interaction in (false,true),
        images in (false,true),trailing in (false,true),shared_ids in (false,true)

        input = deepcopy(fixture)
        input.nwake .= rows
        input.symmetric .= [images,false,false]
        input.trailing_vortices .= trailing
        input.surface_id .= shared_ids ? [1,1,1] : [1,2,3]
        input.wake_finite_core .= [false,true,true]
        options = (; additional_velocity=additional,repeated_points=repeated,
            nwake=rows,interaction_id=[1,1,2],interaction)
        expected = deepcopy(input)
        reference_velocities!(expected; options...)
        @test all(group -> all(v -> all(isfinite,v),group),expected.V)
        for threaded in (false,true)
            actual = deepcopy(input)
            F.fast_wake_velocities!(actual,actual.freestream[]; options...,threaded)
            @test actual.V == expected.V
            @test actual.wakes == input.wakes
            @test actual.Γ == input.Γ
        end
        expected = deepcopy(input)
        W.advance_wake!(expected,expected.freestream[],0.001; options...)
        actual = deepcopy(input)
        F.advance_fast_wake!(actual,actual.freestream[],0.001; options...)
        @test actual.V == expected.V
        @test actual.wakes == expected.wakes
        @test actual.nwake == expected.nwake
        @test actual.Γ == expected.Γ
    end
    input = deepcopy(fixture)
    input.nwake .= [4,4,3]
    expected = deepcopy(input)
    W.advance_wake!(expected,expected.freestream[],0.001)
    for _ = 1:4
        actual = deepcopy(input)
        fetch(Threads.@spawn F.advance_fast_wake!(actual,actual.freestream[],0.001))
        @test actual.V == expected.V
        @test actual.wakes == expected.wakes
    end
    @test_throws ArgumentError F.advance_fast_wake!(input,input.freestream[],0.0)
    @test_throws DimensionMismatch F.advance_fast_wake!(input,input.freestream[],0.01; nwake=[1])
    @test_throws ArgumentError F.advance_fast_wake!(input,input.freestream[],0.01; nwake=[5,4,3])
end

@testset "A shorter shared-edge wake is evaluated at the receiving coordinate" begin
    input = wake_fixture()
    input.nwake .= [4,1,3]
    input.Γ .= 0.0
    input.surface_id .= 1:3
    input.wake_finite_core .= true
    input.trailing_vortices .= false
    # Isolate induction from the short sending wake. The other wake panels
    # retain their geometry but carry zero circulation.
    for isurf in (1,3),I in CartesianIndices(input.wakes[isurf])
        panel = input.wakes[isurf][I]
        input.wakes[isurf][I] = W.WakePanel(panel.rtl,panel.rtr,panel.rbl,panel.rbr,
            W.get_core_size(panel),0.0)
    end
    I = CartesianIndex(4,3)
    rc = W.top_right(input.wakes[1][4,end])
    expected =
        W.freestream_velocity(input.freestream[]) +
        W.rotational_velocity(rc,input.freestream[],input.reference[]) +
        W.induced_velocity(rc,input.wakes[2]; nc=1,finite_core=true,
            symmetric=false,trailing_vortices=false,xhat=input.xhat[])
    options = (; additional_velocity=nothing,
        repeated_points=W.repeated_trailing_edge_points(input.surfaces),
        nwake=input.nwake,interaction_id=input.surface_id,interaction=true)
    reference = deepcopy(input)
    reference_velocities!(reference; options...)
    @test reference.V[1][I] ≈ expected rtol = 1e-14
    for threaded in (false,true)
        actual = deepcopy(input)
        F.fast_wake_velocities!(actual,actual.freestream[]; options...,threaded)
        @test actual.V[1][I] ≈ expected rtol = 1e-14
        @test actual.V == reference.V
    end
end

function small_configuration(directory; integrator=:newmark_beta,coupling=:loose_explicit)
    defaults = C.chang_case_defaults()
    defaults = merge(
        defaults,
        (;
            excitation=merge(
                defaults.excitation,
                (; impulse_start_s=0.003,
                    impulse_duration_s=0.002,impulse_magnitude_nm=1.0)
            ),
            integration=merge(defaults.integration,(; time_integrator=integrator)),
            coupling=merge(defaults.coupling,(; scheme=coupling))
        )
    )
    return C.load_chang_configuration(
        defaults;
        env=Dict(
            "CHANG_WING_SPAN_PANELS" => "4","CHANG_WING_CHORD_PANELS" => "2",
            "CHANG_PROP_RADIAL_PANELS" => "2","CHANG_PROP_CHORD_PANELS" => "2",
            "CHANG_WAKE_ROWS_WING" => "3","CHANG_WAKE_ROWS_PROPELLER" => "4",
            "CHANG_END_TIME_S" => "0.01","CHANG_AZIMUTH_STEP_DEG" => "5",
            "CHANG_PLOT_RESULTS" => "false","CHANG_ANIMATE_WAKE" => "false",
            "CHANG_OUTPUT_DIR" => directory
        )
    )
end

function compare_runs(reference,actual)
    for field in (:displacement_history,:velocity_history,:acceleration_history,
        :trim_load,:coupling_iterations,:coupling_state_residual,
        :coupling_load_residual,:coupling_equilibrium_residual,:coupling_converged,:last_step)
        @test isequal(getproperty(reference.solution,field),getproperty(actual.solution,field))
    end
    @test actual.solution.last_step == length(actual.model.parameters.dt)
    @test any(q -> !iszero(q),actual.solution.displacement_history)
    @test actual.workspace.system.wakes == reference.workspace.system.wakes
    @test actual.workspace.system.Γ == reference.workspace.system.Γ
end

@testset "Full Chang solver, history opt-in and fresh output directories" begin
    for integrator in (:newmark_beta,:generalized_alpha),
        coupling in (:loose_explicit,:implicit_predictor_corrector)

        mktempdir() do directory
            config = small_configuration(directory; integrator,coupling)
            reference = C.run_chang(config)
            baselineFiles = Dict(
                name => read(joinpath(directory,name))
                for name in readdir(directory) if isfile(joinpath(directory,name))
            )
            actual = F.run_chang_fast(config; write_results=false)
            compare_runs(reference,actual)
            @test isempty(actual.workspace.save)
            @test isempty(actual.workspace.surface_history)
            @test Base.summarysize(actual.workspace.surface_history) <
                  Base.summarysize(reference.workspace.surface_history)
            saved = F.run_chang_fast(config; save_surface_history=true,write_results=false)
            compare_runs(reference,saved)
            @test saved.workspace.surface_history == reference.workspace.surface_history
            if integrator == :newmark_beta && coupling == :loose_explicit
                exported = F.run_chang_fast(config)
                compare_runs(reference,exported)
                @test isdir(exported.config.output.directory)
                @test dirname(exported.config.output.directory) == joinpath(directory,"fast_runs")
                @test config.output.directory == directory
                @test isfile(joinpath(exported.config.output.directory,
                    config.output.label * "_history.csv"))
                @test !isempty(readdir(exported.config.output.directory))
            end
            for (name,bytes) in baselineFiles
                @test read(joinpath(directory,name)) == bytes
            end
        end
    end
end

@testset "Per-propeller wake capacities" begin
    mktempdir() do directory
        config = small_configuration(directory)
        config = merge(config,(; propeller=merge(config.propeller,
            (; attachment_eta=[0.4,0.8]))))
        model = C.build_chang_model(config)
        options = merge(model.aerodynamic_options,(; maximum_wake_rows_propeller=[0,2]))
        workspace = C.initialize_chang_uvlm(model.parameters,options)
        blades = config.propeller.blades
        @test workspace.nwake == [3; zeros(Int,blades); fill(2,blades)]
        @test size.(workspace.system.wakes,1) == workspace.nwake
        @test_throws ArgumentError C.initialize_chang_uvlm(model.parameters,
            merge(options,(; maximum_wake_rows_propeller=[-1,2])))
        @test_throws DimensionMismatch C.initialize_chang_uvlm(model.parameters,
            merge(options,(; maximum_wake_rows_propeller=[2])))
    end
end

@testset "Limit stops retain the final accepted surface geometry" begin
    mktempdir() do directory
        config = small_configuration(directory)
        config = merge(config,(; integration=merge(config.integration,
            (; state_norm_limit=eps(Float64)))))
        reference = C.run_chang(config)
        actual = F.run_chang_fast(config; save_surface_history=true,write_results=false)
        last_step = reference.solution.last_step
        @test 0 < last_step < length(reference.model.parameters.dt)
        @test !reference.results.validation_passed
        @test actual.solution.last_step == last_step
        @test actual.solution.displacement_history[1:last_step+1] ==
              reference.solution.displacement_history[1:last_step+1]
        @test actual.workspace.system.wakes == reference.workspace.system.wakes
        for run in (reference,actual)
            @test isassigned(run.workspace.surface_history,last_step)
            if isassigned(run.workspace.surface_history,last_step)
                @test run.workspace.surface_history[last_step] == run.workspace.system.surfaces
            end
        end
    end
end

@testset "Animation frame storage is independent of full geometry history" begin
    mktempdir() do directory
        model = C.build_chang_model(small_configuration(directory))
        workspace = F.build_fast_workspace(model)
        wake = F.fast_wake_context(workspace; interaction_on=model.config.simulation.interaction_on)
        function record(surfaces,wakes,rows,times,time)
            push!(surfaces,deepcopy(workspace.system.surfaces))
            push!(wakes,deepcopy(workspace.system.wakes))
            push!(rows,copy(workspace.iwake))
            push!(times,time)
        end
        solution = C.solve_chang_aeroelastic!(workspace.system,model.structural;
            time=model.parameters.t,time_steps=model.parameters.dt,
            freestream_history=workspace.fs_vec,wake,
            aerodynamic_load=(snapshot,state,step) -> C.aero_load_for_state!(
                model,workspace,snapshot,state,step),
            excitation=C.chang_excitation_options(model.config,model.parameters),
            integration=C.chang_integration_options(model.config,model.structural,model.parameters),
            animation=(; enabled=true,stride=3,record_frame=record))
        @test isempty(workspace.surface_history)
        @test length(solution.animation_time_history) > 2
        @test length(solution.animation_surface_history) == length(solution.animation_time_history)
        @test last(solution.animation_time_history) == model.parameters.t[end]
    end
end

@info "Chang runtime optimization tests completed" threads = Threads.nthreads()
end # module
