# Bounded performance diagnostics; this does not run the LCO speed sweep.
# Run from the repository root:
#   julia --startup-file=no --project=. docs/reviews/efficiency_2026-10-19/measure_allocations.jl
using LinearAlgebra, Statistics, TOML, Profile
# Load the root project before the example adds its local package environments.
using AeroBeams

const ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(ROOT, "lib", "WingPropellerUVLM", "examples",
    "pazy_aerobeams_coupled", "PazyWingUVLMCoupling.jl"))
const AB = AeroBeams
const MEASUREMENTS = Dict{String,Any}[]

function measure(label, f; repetitions=7)
    f() # Compile the operation before recording it.
    GC.gc()
    elapsed = Float64[]
    allocated = Int[]
    for _ in 1:repetitions
        sample = @timed f()
        push!(elapsed, sample.time)
        push!(allocated, sample.bytes)
        sample = nothing
    end
    entry = Dict{String,Any}("operation" => label,
        "median_seconds" => median(elapsed),
        "minimum_seconds" => minimum(elapsed),
        "median_bytes" => median(allocated))
    push!(MEASUREMENTS, entry)
    println(entry)
    return entry
end

function make_structure(steps)
    ne, span, _, _ = AB.geometrical_properties_Pazy()
    beam = AB.create_Beam(name="performance review", length=span, nElements=ne,
        normalizedNodalPositions=AB.nodal_positions_Pazy(),
        S=AB.stiffness_matrices_Pazy(withSkin=true,
            sweepStructuralCorrections=false, GAy=1e16, GAz=1e16, Λ=0.0),
        I=AB.inertia_matrices_Pazy(withSkin=true),
        rotationParametrization="E231", p0=[-pi/2, 0.0, 0.0], aeroSurface=nothing)
    root = AB.create_BC(name="root", beam=beam, node=1,
        types=["u1A", "u2A", "u3A", "p1A", "p2A", "p3A"], values=zeros(6))
    model = AB.create_Model(name="performance review", beams=[beam], BCs=[root],
        gravityVector=zeros(3), v_A=t -> zeros(3))
    problem = AB.create_DynamicProblem(model=model,
        timeVector=collect(0:steps) .* 0.0001,
        trackingTimeSteps=false, displayProgress=false)
    AB.precompute_distributed_loads!(problem)
    return problem
end

function rectangular_grid(nc=10, ns=20)
    grid = zeros(3, nc+1, ns+1)
    for j in 1:ns+1, i in 1:nc+1
        grid[1,i,j] = (i-1)/nc * 0.0989
        grid[2,i,j] = (j-1)/ns * 0.549843728
    end
    return grid
end

# Review-only prototype. The stencils are prepared once from the exact weights
# used by the existing routine, preserving its force and moment transfer.
function stencil_beam_loads!(out, loads, offsets, weights, stencils)
    fill!(out, 0)
    forces = loads.forces[1]
    @inbounds for j in axes(forces,2), k in axes(forces,1)
        f = forces[k,j]
        force = UVLM.SVector(-f[3], -f[1], f[2])
        for node in stencils[j]
            weight = weights[j,node]
            offset = UVLM.SVector(offsets[1,k,node], offsets[2,k,node], offsets[3,k,node])
            moment = cross(offset,force)
            for d in 1:3
                out[d,node] += weight * force[d]
                out[d+3,node] += weight * moment[d]
            end
        end
    end
    return out
end

function interpolate_then_scale(problem, array, factor)
    first_node, second_node = AB.interpolate_distributed_loads(problem,array)
    first_node .*= factor
    second_node .*= factor
    return first_node, second_node
end

function repeat_trial!(aero, grid, count)
    for _ in 1:count
        UVLM.evaluate_trial!(aero;surfaces=[grid])
    end
    return nothing
end

function main()
    metadata = Dict{String,Any}("julia_version" => string(VERSION),
        "julia_threads" => Threads.nthreads(), "blas_threads" => BLAS.get_num_threads(),
        "date" => "2026-10-19", "repetitions" => 7,
        "scope" => "isolated operations; no full coupled simulation")
    println(metadata)
    distributed_fields = (:f_A,:m_A,:f_b,:m_b,:ff_A,:mf_A,:ff_b,:mf_b)
    sizes = Dict{String,Any}[]
    short_problem = nothing
    for steps in (256,4096)
        problem = make_structure(steps)
        short_problem = steps == 256 ? problem : short_problem
        load_bytes = sum(sizeof(getproperty(e,field)) for e in problem.model.elements
            for field in distributed_fields)
        all_zero = all(all(iszero,getproperty(e,field)) for e in problem.model.elements
            for field in distributed_fields)
        entry = Dict("steps" => steps, "model_bytes" => Base.summarysize(problem.model),
            "distributed_load_bytes" => load_bytes, "all_load_arrays_zero" => all_zero)
        push!(sizes,entry)
        println(entry)
        measure("AeroBeams model deepcopy, $steps steps", () -> deepcopy(problem.model))
    end

    problem = short_problem
    AB.update_time_variables!(problem,2)
    AB.get_equivalent_states_rates!(problem)
    measure("AeroBeams element rotation variables", () -> AB.element_rotation_variables!(problem,problem.model.elements[1]))
    measure("AeroBeams system assembly, 15 elements", () -> AB.assemble_system_arrays!(problem))
    for steps in (256,4096)
        array = reshape(sin.(1.0:(6*(steps+1))),2,3,steps+1)
        reference = AB.interpolate_distributed_loads(problem,0.5*array)
        candidate = interpolate_then_scale(problem,array,0.5)
        all(isapprox(a,b) for (a,b) in zip(reference,candidate)) || error("Scaled interpolation differs")
        measure("Scale whole load history then interpolate, $steps steps",
            () -> AB.interpolate_distributed_loads(problem,0.5*array))
        measure("Interpolate then scale six values, $steps steps",
            () -> interpolate_then_scale(problem,array,0.5))
    end
    base_rotation = copy(problem.model.R_A)
    append!(problem.model.R_A_ofTime,[copy(base_rotation) for _ in 1:4096])
    measure("AeroBeams model deepcopy, 256 steps plus 4096 orientation samples", () -> deepcopy(problem.model))

    grid = rectangular_grid()
    model = UVLM.create_UVLMModel(surfaces=[grid],
        reference=UVLM.Reference(0.549843728*0.0989,0.0989,0.549843728,zeros(3),30.0,1.225),
        symmetric=true, interactionGroups=[1])
    aero = UVLM.create_UVLMDynamicProblem(model=model,
        operatingPoint=UVLM.create_OperatingPoint(airspeed=30.0,density=1.225,angleOfAttack=0.01,sideslip=0.0),
        timeVector=collect(0:22) .* 0.0005, maximumWakeRows=20, trackingTimeSteps=false)
    for _ in 1:2
        UVLM.begin_time_step!(aero)
        UVLM.evaluate_trial!(aero;surfaces=[grid])
        UVLM.commit_time_step!(aero)
    end
    UVLM.begin_time_step!(aero)
    measure("UVLM full trial, 200 bound panels, 2 active wake rows", () -> UVLM.evaluate_trial!(aero;surfaces=[grid]))
    sys = aero.workspace.system
    metadata["uvlm_system_bytes"] = Base.summarysize(sys)
    measure("UVLM System deepcopy, wake capacity 20 rows", () -> deepcopy(sys))
    snapshot = UVLM.snapshot_uvlm(sys)
    measure("UVLM snapshot allocation", () -> UVLM.snapshot_uvlm(sys))
    measure("UVLM snapshot restore in place", () -> UVLM.restore_uvlm!(sys,snapshot))
    measure("UVLM apply_trial_geometry!", () -> UVLM.apply_trial_geometry!(sys,[grid],aero.aeroSolver))
    measure("UVLM existing update_surface_panels!", () -> UVLM.update_surface_panels!(sys.surfaces[1],grid;
        ratios=sys.ratios[1],fcore=(c,ds)->aero.aeroSolver.coreRadius))

    keywords = (; symmetric=sys.symmetric, wake_shedding_locations=sys.wake_shedding_locations,
        surface_id=sys.surface_id, trailing_vortices=sys.trailing_vortices, xhat=sys.xhat[],
        interaction_id=model.interactionGroups, interaction=true)
    trailing_edge_checks = Dict{String,Any}[]
    for symmetry in (false,true), finite_core in (false,true)
        local_keywords = (; symmetric=symmetry, finite_core,
            wake_shedding_locations=sys.wake_shedding_locations[1],
            trailing_vortices=sys.trailing_vortices[1], xhat=sys.xhat[])
        original = similar(sys.AIC)
        UVLM.influence_coefficients!(original,sys.surfaces[1],sys.surfaces[1];local_keywords...)
        updated = copy(original)
        UVLM.update_trailing_edge_coefficients!(updated,sys.surfaces[1],sys.surfaces[1];local_keywords...)
        difference = maximum(abs,updated-original)
        push!(trailing_edge_checks,Dict("symmetric"=>symmetry,"finite_core"=>finite_core,
            "maximum_difference"=>difference))
        println((check="full AIC versus full AIC plus TE update",symmetry,finite_core,difference))
        difference <= 1e-12 || error("Trailing-edge update was not redundant for this configuration")
    end
    metadata["trailing_edge_checks"] = trailing_edge_checks
    measure("UVLM full AIC assembly", () -> UVLM.influence_coefficients!(sys.AIC,sys.surfaces;keywords...))
    measure("UVLM trailing-edge coefficient update", () -> UVLM.update_trailing_edge_coefficients!(sys.AIC,sys.surfaces;keywords...))
    original_blas_threads = BLAS.get_num_threads()
    try
        for count in unique([original_blas_threads,1])
            BLAS.set_num_threads(count)
            measure("UVLM dense LU and solve, BLAS $count threads", () -> UVLM.circulation!(sys.Γ,sys.AIC,sys.w))
        end
    finally
        BLAS.set_num_threads(original_blas_threads)
    end
    measure("UVLM near-field forces", () -> UVLM.near_field_forces!(sys.properties,
        sys.surfaces,sys.wakes,sys.reference[],sys.freestream[],sys.Γ;
        dΓdt=sys.dΓdt,additional_velocity=nothing,Vh=sys.Vh,Vv=sys.Vv,
        symmetric=sys.symmetric,nwake=sys.nwake,surface_id=sys.surface_id,
        wake_finite_core=sys.wake_finite_core,wake_shedding_locations=sys.wake_shedding_locations,
        trailing_vortices=sys.trailing_vortices,xhat=sys.xhat[],interaction=true))

    measure("AeroBeams save complete structural step", () -> AB.save_time_step_data!(problem,problem.timeNow))
    metadata["saved_x_aliases_active_x"] = problem.xOverTime[end] === problem.x
    println((saved_x_aliases_active_x=metadata["saved_x_aliases_active_x"],))
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
    stencils = [findall(!iszero,weights[j,:]) for j in axes(weights,1)]
    offsets = reshape(sin.(1.0:(3*11*nn)),3,11,nn) .* 0.02
    forces = [UVLM.SVector(sin(k+j),cos(k-2j),sin(k-3j)) for k in 1:11,j in 1:21]
    loads = (; forces=[forces])
    expected = beam_loads(loads,offsets,weights)
    out = similar(expected)
    stencil_beam_loads!(out,loads,offsets,weights,stencils)
    metadata["stencil_load_transfer_max_difference"] = maximum(abs,out-expected)
    isapprox(out,expected;rtol=1e-13,atol=1e-13) || error("Load transfer prototype differs")
    measure("Pazy beam_loads original", () -> beam_loads(loads,offsets,weights))
    measure("Pazy beam_loads stencil prototype", () -> stencil_beam_loads!(out,loads,offsets,weights,stencils))
    measure("Pazy wing_geometry original", () -> wing_geometry(problem.model,0.0989,0.44096,10,weights))

    # Fill the wake by actual accepted aerodynamic steps, then repeat trials
    # from one accepted state. Convection remains a separate physical update.
    UVLM.evaluate_trial!(aero;surfaces=[grid])
    UVLM.commit_time_step!(aero)
    while aero.state.stepIndex < 21
        UVLM.begin_time_step!(aero)
        UVLM.evaluate_trial!(aero;surfaces=[grid])
        UVLM.commit_time_step!(aero)
    end
    UVLM.begin_time_step!(aero)
    measure("UVLM full trial, 200 bound panels, 20 active wake rows",
        () -> UVLM.evaluate_trial!(aero;surfaces=[grid]))
    full_sys = aero.workspace.system
    measure("UVLM wake velocity evaluation, 20 active rows",
        () -> UVLM.get_wake_velocities!(full_sys.V,full_sys.surfaces,full_sys.wakes,
            full_sys.reference[],full_sys.freestream[],full_sys.Γ,nothing,full_sys.Vte,full_sys.symmetric,
            UVLM.repeated_trailing_edge_points(full_sys.surfaces),full_sys.nwake,
            full_sys.surface_id,full_sys.wake_finite_core,full_sys.wake_shedding_locations,
            full_sys.trailing_vortices,full_sys.xhat[]))
    repeat_trial!(aero,grid,1)
    Profile.clear()
    Profile.@profile repeat_trial!(aero,grid,24)
    open(joinpath(@__DIR__,"cpu_profile.txt"),"w") do io
        Profile.print(IOContext(io,:displaysize=>(1000,220));
            format=:flat,sortedby=:count,mincount=4)
    end

    open(joinpath(@__DIR__,"measurements.toml"),"w") do io
        TOML.print(io,Dict("metadata"=>metadata,"model_sizes"=>sizes,"operations"=>MEASUREMENTS))
    end
end

main()
