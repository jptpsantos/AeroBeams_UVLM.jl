# Allocation measurements of the implemented routines; preserves the original
# measurements.toml and cpu_profile.txt. Run from the root project.
include(joinpath(@__DIR__,"measure_allocations.jl"))

function measure_implementation()
    empty!(MEASUREMENTS)
    metadata = Dict{String,Any}("julia_version"=>string(VERSION),
        "julia_threads"=>Threads.nthreads(),"blas_threads"=>BLAS.get_num_threads(),
        "date"=>"2026-10-19","repetitions"=>7,
        "baseline_commit"=>"0b02c08730d8072ea46383766863bc888e761017",
        "scope"=>"warmed isolated operations; no full LCO sweep")
    sizes = Dict{String,Any}[]
    problem = make_structure(256)
    for steps in (256,4096)
        p = steps == 256 ? problem : make_structure(steps)
        push!(sizes,Dict("steps"=>steps,"model_bytes"=>Base.summarysize(p.model),
            "distributed_load_bytes"=>sum(sizeof(getfield(e,f)) for e in p.model.elements for f in AB.DISTRIBUTED_LOAD_TABLES)))
        measure("AeroBeams model deepcopy, $steps steps",()->deepcopy(p.model))
        measure("AeroBeams Newton checkpoint, $steps steps",()->AB.checkpoint_model(p.model))
        AB.add_loads_to_beam!(p.model.beams[1];loadTypes=["f_A_of_x1t"],
            loadFuns=[(x,t)->[sin(x)+t,0.0,0.0]])
        AB.precompute_distributed_loads!(p)
        measure("AeroBeams active-load Newton checkpoint, $steps steps",()->AB.checkpoint_model(p.model))
        measure("AeroBeams active-load deepcopy, $steps steps",()->deepcopy(p.model))
    end
    # Assembly fixture is unloaded, like the original measurement.
    problem = make_structure(256)
    AB.update_time_variables!(problem,2)
    AB.get_equivalent_states_rates!(problem)
    measure("AeroBeams element rotation variables",()->AB.element_rotation_variables!(problem,problem.model.elements[1]))
    measure("AeroBeams system assembly, 15 elements",()->AB.assemble_system_arrays!(problem))
    for steps in (256,4096)
        array = reshape(sin.(1.0:(6*(steps+1))),2,3,steps+1)
        measure("Scaled distributed-load interpolation, $steps steps",
            ()->AB.interpolate_distributed_loads(problem,array,0.5))
    end
    append!(problem.model.R_A_ofTime,[copy(problem.model.R_A) for _ in 1:4096])
    measure("AeroBeams Newton checkpoint plus 4096 orientation samples",()->AB.checkpoint_model(problem.model))
    measure("AeroBeams save complete structural step",()->AB.save_time_step_data!(problem,problem.timeNow))
    metadata["saved_x_aliases_active_x"] = problem.xOverTime[end] === problem.x

    grid = rectangular_grid()
    nn = length(problem.model.elements)+1
    weights = zeros(21,nn)
    beam_nodes = AB.nodal_positions_Pazy()
    for j in 1:21
        eta = (j-1)/20
        left = min(searchsortedlast(beam_nodes,eta),nn-1)
        fraction = (eta-beam_nodes[left])/(beam_nodes[left+1]-beam_nodes[left])
        weights[j,left] = 1-fraction
        weights[j,left+1] = fraction
    end
    map = PazySpanMap(weights)
    offsets = reshape(sin.(1.0:(3*11*nn)),3,11,nn) .* 0.02
    forces = [UVLM.SVector(sin(k+j),cos(k-2j),sin(k-3j)) for k in 1:11,j in 1:21]
    loads = (;forces=[forces])
    out = zeros(6,nn)
    measure("Pazy beam_loads! reusable",()->beam_loads!(out,loads,offsets,map))
    geometry_grid = similar(grid)
    geometry_offsets = similar(offsets)
    geometry_work = (positions=zeros(3,nn),beam_grid=zeros(3,11,nn))
    measure("Pazy wing_geometry! reusable",()->wing_geometry!(geometry_grid,geometry_offsets,
        problem.model,0.0989,0.44096,10,map,geometry_work))

    model = UVLM.create_UVLMModel(surfaces=[grid],
        reference=UVLM.Reference(0.549843728*0.0989,0.0989,0.549843728,zeros(3),30.0,1.225),
        symmetric=true,interactionGroups=[1])
    aero = UVLM.create_UVLMDynamicProblem(model=model,
        operatingPoint=UVLM.create_OperatingPoint(airspeed=30.0,density=1.225,angleOfAttack=0.01,sideslip=0.0),
        timeVector=collect(0:22) .* 0.0005,maximumWakeRows=20,trackingTimeSteps=false)
    for _ in 1:2
        UVLM.begin_time_step!(aero)
        UVLM.evaluate_trial!(aero;surfaces=[grid])
        UVLM.commit_time_step!(aero)
    end
    UVLM.begin_time_step!(aero)
    measure("UVLM full trial, 200 bound panels, 2 active wake rows",()->UVLM.evaluate_trial!(aero;surfaces=[grid]))
    sys = aero.workspace.system
    metadata["uvlm_system_bytes"] = Base.summarysize(sys)
    measure("UVLM copy_system! into existing trial",()->UVLM.copy_system!(sys,aero.state.system))
    UVLM.evaluate_trial!(aero;surfaces=[grid])
    measure("UVLM apply_trial_geometry!",()->UVLM.apply_trial_geometry!(sys,[grid],aero.aeroSolver))
    measure("UVLM near-field forces with System buffers",()->UVLM.near_field_forces!(sys.properties,
        sys.surfaces,sys.wakes,sys.reference[],sys.freestream[],sys.Γ;
        dΓdt=sys.dΓdt,additional_velocity=nothing,Vh=sys.Vh,Vv=sys.Vv,
        symmetric=sys.symmetric,nwake=sys.nwake,surface_id=sys.surface_id,
        wake_finite_core=sys.wake_finite_core,wake_shedding_locations=sys.wake_shedding_locations,
        trailing_vortices=sys.trailing_vortices,xhat=sys.xhat[],interaction=true,
        chord_force_buffers=sys.chord_seg_forces,span_force_buffers=sys.span_seg_forces,
        unsteady_force_buffers=sys.unsteady_forces,span_velocity_buffers=sys.span_force_velocities))
    UVLM.commit_time_step!(aero)
    while aero.state.stepIndex < 21
        UVLM.begin_time_step!(aero)
        UVLM.evaluate_trial!(aero;surfaces=[grid])
        UVLM.commit_time_step!(aero)
    end
    UVLM.begin_time_step!(aero)
    measure("UVLM full trial, 200 bound panels, 20 active wake rows",()->UVLM.evaluate_trial!(aero;surfaces=[grid]))
    full = aero.workspace.system
    measure("UVLM wake velocity evaluation, 20 active rows",()->UVLM.get_wake_velocities!(full.V,
        full.surfaces,full.wakes,full.reference[],full.freestream[],full.Γ,nothing,full.Vte,full.symmetric,
        UVLM.repeated_trailing_edge_points(full.surfaces),full.nwake,full.surface_id,full.wake_finite_core,
        full.wake_shedding_locations,full.trailing_vortices,full.xhat[]))
    Profile.clear()
    Profile.@profile repeat_trial!(aero,grid,24)
    open(joinpath(@__DIR__,"cpu_profile_after.txt"),"w") do io
        Profile.print(IOContext(io,:displaysize=>(1000,220));format=:flat,sortedby=:count,mincount=4)
    end
    open(joinpath(@__DIR__,"measurements_after.toml"),"w") do io
        TOML.print(io,Dict("metadata"=>metadata,"model_sizes"=>sizes,"operations"=>MEASUREMENTS))
    end
end

abspath(PROGRAM_FILE) == (@__FILE__) && measure_implementation()
