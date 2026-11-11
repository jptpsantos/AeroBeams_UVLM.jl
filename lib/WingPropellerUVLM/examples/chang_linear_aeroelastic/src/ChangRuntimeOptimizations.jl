module ChangRuntimeOptimizations

import WingPropellerUVLM as W
import ..ChangAeroelastic as C
using Dates

export run_chang_fast,build_fast_workspace,fast_wake_context,advance_fast_wake!

"""Wake controls owned by this optional runner, with a specialized commit hook."""
struct FastChangWake{R,S,H}
    surface_count::Int
    repeated_points::R
    maximum_rows::Vector{Int}
    active_rows::Vector{Int}
    interaction_ids::Vector{Int}
    interaction_on::Bool
    saved_steps::S
    surface_history::H
    threaded::Bool
end

function build_fast_workspace(model; save_surface_history::Bool=false)
    workspace = C.build_chang_workspace(model)
    save_surface_history && return workspace
    return merge(workspace,(;
        save=1:0,
        surface_history=empty(workspace.surface_history)
    ))
end

function fast_wake_context(workspace; interaction_on::Bool,threaded::Bool=true)
    return FastChangWake(
        workspace.nsurf,workspace.repeated_points,workspace.nwake,
        workspace.iwake,workspace.surface_interaction_id,interaction_on,
        workspace.save,workspace.surface_history,threaded
    )
end

# Each call owns its accumulation variables. Sending panels and circulation
# are read-only; no reduction order changes when receiving vertices run in parallel.
function wake_vertex_velocity(isurf,I,data)
    (; surfaces,wakes,ref,fs,circulation,additional_velocity,symmetric,
        repeated_points,nwake,surface_id,wake_finite_core,shedding_locations,
        trailing_vortices,xhat,interaction_id,interaction,offsets) = data
    nw,ns = nwake[isurf],size(surfaces[isurf],2)
    if I[1] <= nw && I[2] <= ns
        rc = W.top_left(wakes[isurf][I])
    elseif I[1] == nw + 1 && I[2] <= ns
        rc = W.bottom_left(wakes[isurf][I[1]-1,I[2]])
    elseif I[1] <= nw && I[2] == ns + 1
        rc = W.top_right(wakes[isurf][I[1],I[2]-1])
    else
        rc = W.bottom_right(wakes[isurf][I[1]-1,I[2]-1])
    end

    velocity = W.freestream_velocity(fs)
    velocity += W.rotational_velocity(rc,fs,ref)
    if !isnothing(additional_velocity)
        velocity += additional_velocity(rc)
    end
    for jsurf in eachindex(surfaces)
        if !interaction && interaction_id[isurf] != interaction_id[jsurf]
            continue
        end
        sameSurface = isurf == jsurf
        js = I[2]
        if !sameSurface && haskey(repeated_points,(isurf,I[2]))
            duplicates = repeated_points[(isurf,I[2])]
            index = findfirst(point -> point[1] == jsurf,duplicates)
            if !isnothing(index)
                sameSurface = true
                js = duplicates[index][2]
            end
        end
        gamma = view(circulation,offsets[jsurf]+1:offsets[jsurf+1])
        velocity += W.induced_velocity(rc,surfaces[jsurf],gamma;
            finite_core=surface_id[isurf] != surface_id[jsurf],
            wake_shedding_locations=shedding_locations[jsurf],
            symmetric=symmetric[jsurf],trailing_vortices=false,xhat)
        if sameSurface && I[1] <= nwake[jsurf] + 1
            velocity += W.induced_velocity(CartesianIndex(I[1],js),wakes[jsurf];
                finite_core=wake_finite_core[jsurf] || surface_id[isurf] != surface_id[jsurf],
                symmetric=symmetric[jsurf],nc=nwake[jsurf],
                trailing_vortices=trailing_vortices[jsurf],xhat)
        else
            velocity += W.induced_velocity(rc,wakes[jsurf];
                finite_core=wake_finite_core[jsurf] || surface_id[isurf] != surface_id[jsurf],
                symmetric=symmetric[jsurf],nc=nwake[jsurf],
                trailing_vortices=trailing_vortices[jsurf],xhat)
        end
    end
    return velocity
end

function fill_wake_surface_velocities!(output,isurf,data,Vte,threaded)
    ns = size(data.surfaces[isurf],2)
    for j = 1:ns+1
        rc = data.shedding_locations[isurf][j]
        velocity = W.freestream_velocity(data.fs)
        velocity += W.rotational_velocity(rc,data.fs,data.ref)
        if !isnothing(data.additional_velocity)
            velocity += data.additional_velocity(rc)
        end
        if !isnothing(Vte)
            velocity += Vte[isurf][j]
        end
        output[1,j] = velocity
    end
    vertices = CartesianIndices((2:data.nwake[isurf]+1,1:ns+1))
    if threaded && Threads.nthreads() > 1 && length(vertices) > 1
        # Dynamic scheduling also works when invoked from a worker task.
        Threads.@threads :dynamic for index = 1:length(vertices)
            I = vertices[index]
            output[I] = wake_vertex_velocity(isurf,I,data)
        end
    else
        for I in vertices
            output[I] = wake_vertex_velocity(isurf,I,data)
        end
    end
    return output
end

"""Evaluate independent wake vertices; callbacks must support concurrent reads."""
function fast_wake_velocities!(system,fs;
    additional_velocity=nothing,
    repeated_points=W.repeated_trailing_edge_points(system.surfaces),
    nwake=system.nwake,interaction_id=system.surface_id,
    interaction::Bool=true,threaded::Bool=true)

    offsets = cumsum(vcat(0,length.(system.surfaces)))
    data = (;
        surfaces=system.surfaces,wakes=system.wakes,ref=system.reference[],fs,
        circulation=system.Γ,additional_velocity,symmetric=system.symmetric,
        repeated_points,nwake,surface_id=system.surface_id,
        wake_finite_core=system.wake_finite_core,
        shedding_locations=system.wake_shedding_locations,
        trailing_vortices=system.trailing_vortices,xhat=system.xhat[],
        interaction_id,interaction,offsets
    )
    for isurf in eachindex(system.surfaces)
        fill_wake_surface_velocities!(system.V[isurf],isurf,data,system.Vte,threaded)
    end
    return system.V
end

"""Commit one wake step with the original shedding and stretching equations."""
function advance_fast_wake!(system,fs,dt;
    additional_velocity=nothing,
    repeated_points=W.repeated_trailing_edge_points(system.surfaces),
    nwake=system.nwake,interaction_id=system.surface_id,
    interaction::Bool=true,threaded::Bool=true)

    dt > 0 || throw(ArgumentError("dt must be positive"))
    length(nwake) == length(system.surfaces) || throw(DimensionMismatch(
        "nwake must contain one active-row count per surface"))
    all((0 .<= nwake) .& (nwake .<= size.(system.wakes,1))) ||
        throw(ArgumentError("active wake rows must lie within allocated capacity"))
    system.freestream[] = fs
    system.nwake .= nwake
    fast_wake_velocities!(system,fs;
        additional_velocity,repeated_points,nwake,interaction_id,interaction,threaded)
    W.shed_wake!(system.wakes,system.wake_shedding_locations,system.V,
        dt,system.surfaces,system.Γ,nwake)
    return system
end

# This adds a method for our own context type. Existing package methods and
# the original runner keep their original behavior, even in the same session.
function C.commit_chang_wake!(system,freestream,time_step,wake::FastChangWake)
    if wake.threaded && Threads.nthreads() > 1
        advance_fast_wake!(system,freestream,time_step;
            repeated_points=wake.repeated_points,nwake=wake.active_rows,
            interaction_id=wake.interaction_ids,interaction=wake.interaction_on)
    else
        W.advance_wake!(system,freestream,time_step;
            repeated_points=wake.repeated_points,nwake=wake.active_rows,
            interaction_id=wake.interaction_ids,interaction=wake.interaction_on)
    end
    for isurf in eachindex(wake.active_rows)
        wake.active_rows[isurf] < wake.maximum_rows[isurf] || continue
        wake.active_rows[isurf] += 1
    end
    return nothing
end

"""
    run_chang_fast(config=ChangAeroelastic.load_chang_configuration();
        threaded_wake=true, save_surface_history=false, write_results=true)

Run the existing Chang solver with optional parallel wake convection and
optional full surface history. Animation uses the existing independent frame
storage. Result files go into a fresh directory under `output/fast_runs`, or
under the configured output base. Existing outputs are never overwritten.
Start Julia with multiple threads, e.g. `julia --threads=4`, for parallel work.
"""
function run_chang_fast(config=C.load_chang_configuration();
    threaded_wake::Bool=true,save_surface_history::Bool=false,
    write_results::Bool=true)

    if write_results || config.output.animate_wake
        C.load_chang_visualization(config.output)
    end
    return Base.invokelatest(run_fast_case,config;
        threaded_wake,save_surface_history,write_results)
end

function run_fast_case(config; threaded_wake,save_surface_history,write_results)
    if write_results
        base = joinpath(config.output.directory,"fast_runs")
        mkpath(base)
        directory = mktempdir(base; prefix=Dates.format(now(),"yyyymmdd_HHMMSS") * "_",
            cleanup=false)
        config = merge(config,(; output=merge(config.output,(; directory))))
    end
    model = C.build_chang_model(config)
    config = model.config
    parameters,structural = model.parameters,model.structural
    C.report_and_validate_structural_model(structural)
    println("Aeroelastic near-field force model: $(config.simulation.near_field_force_model)")
    println("Propeller moment projection: $(config.simulation.propeller_moment_projection)")
    println("Matrices after BCs. Total DOFs (free): $(structural.ndof_free)")
    println("Full surface geometry history: $save_surface_history")
    println("Threaded wake requested: $threaded_wake; Julia threads: $(Threads.nthreads())")
    if threaded_wake && Threads.nthreads() == 1
        println("Start Julia with --threads=4 to enable parallel wake convection.")
    end
    workspace = build_fast_workspace(model; save_surface_history)
    wake = fast_wake_context(workspace;
        interaction_on=config.simulation.interaction_on,threaded=threaded_wake)
    animation = C.chang_animation_options(workspace.system,workspace.iwake,config.output)
    excitation = C.chang_excitation_options(config,parameters)
    integration = C.chang_integration_options(config,structural,parameters)
    C.report_chang_solver_options(excitation,integration)
    aerodynamic_load = (snapshot,state,step) ->
        C.aero_load_for_state!(model,workspace,snapshot,state,step)
    solution = C.solve_chang_aeroelastic!(workspace.system,structural;
        time=parameters.t,time_steps=parameters.dt,
        freestream_history=workspace.fs_vec,aerodynamic_load,
        wake,excitation,integration,animation)
    results = if write_results
        println("Writing fast-run results to: $(config.output.directory)")
        C.write_chang_results(solution;
            wing=config.wing,propeller=config.propeller,simulation=config.simulation,
            time=parameters.t,time_steps=parameters.dt,
            wing_node_count=parameters.nnodes,dofs_per_node=parameters.ndof,
            density=parameters.ref.rho,visualization=config.output,
            output_directory=config.output.directory,output_label=config.output.label)
    else
        nothing
    end
    return (; config,model,workspace,solution,results)
end

end # module
