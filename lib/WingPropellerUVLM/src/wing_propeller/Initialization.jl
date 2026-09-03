# ==============================================================================
# UVLM system initialization for the Bohnisch wing-propeller time-domain model.
# ==============================================================================

using StaticArrays

function initialize_bohnisch_uvlm_system(;
    xle, yle, zle, chord_geo, theta_geo, phi_geo,
    ns_wing, nc_wing, mirror_wing, spacing_s_wing, spacing_c_wing,
    R_prop, c_prop, ns_prop, nc_prop, blade_twists_prop, Nb_prop,
    Npropellers, span_nodes, prop_attach_nodes, chord, xle_distribution,
    propeller_span_positions = nothing,
    ref, symmetric_wing, fs, dt, nnodes,
    prop_pivot_offset_from_ea_A,
    hub_center_prop_A,
    fcore,
    elastic_axis_fraction::Real = 0.35,
    maximum_wake_rows_wing::Integer = 36,
    maximum_wake_rows_propeller::Integer = 36,
    verbose::Bool = false)

    0.0 <= elastic_axis_fraction <= 1.0 ||
        throw(ArgumentError("elastic_axis_fraction must be between 0 and 1"))
    maximum_wake_rows_wing >= 0 ||
        throw(ArgumentError("maximum_wake_rows_wing must be nonnegative"))
    maximum_wake_rows_propeller >= 0 ||
        throw(ArgumentError("maximum_wake_rows_propeller must be nonnegative"))

    verbose && println("Initializing global UVLM system...")

    #fcore = (c, ds) -> 0.1*ds

    grid_wing_init, ratio_wing = wing_to_grid(xle, yle, zle, chord_geo, theta_geo, phi_geo,
                                              ns_wing, nc_wing;
                                              mirror=mirror_wing,
                                              spacing_s=spacing_s_wing,
                                              spacing_c=spacing_c_wing)
    grids_prop_ref_single = generate_propeller_blades_grid(R_prop, c_prop, ns_prop, nc_prop,
                                                           blade_twists_prop, Nb_prop)
    grids_prop_ref = [deepcopy(grids_prop_ref_single) for _ in 1:Npropellers]
    grids_prop_initial_global = [deepcopy(grids_prop_ref[ip]) for ip in 1:Npropellers]

    if isnothing(propeller_span_positions)
        attach_node_y = [span_nodes[prop_attach_nodes[ip]] for ip in 1:Npropellers]
        attach_node_chord = [chord[prop_attach_nodes[ip]] for ip in 1:Npropellers]
        attach_node_xle = [xle_distribution[prop_attach_nodes[ip]] for ip in 1:Npropellers]
    else
        length(propeller_span_positions) == Npropellers || throw(DimensionMismatch(
            "propeller_span_positions must contain Npropellers entries",
        ))
        attach_node_y = collect(Float64, propeller_span_positions)
        all(first(span_nodes) .<= attach_node_y .<= last(span_nodes)) ||
            throw(ArgumentError("Every propeller span position must lie on the wing"))
        attach_node_chord = linear_interpolate_1d(
            span_nodes,
            chord,
            attach_node_y,
        )
        attach_node_xle = linear_interpolate_1d(
            span_nodes,
            xle_distribution,
            attach_node_y,
        )
    end
    ea_x_aero = attach_node_xle .+ attach_node_chord .* elastic_axis_fraction

    hub_offset_from_ea_A = prop_pivot_offset_from_ea_A + hub_center_prop_A
    T_pivot_global_init = [SVector(ea_x_aero[ip], attach_node_y[ip], 0.0) + prop_pivot_offset_from_ea_A for ip in 1:Npropellers]
    initialize_propeller_grids!(grids_prop_initial_global, grids_prop_ref,
                                T_pivot_global_init, hub_center_prop_A, Nb_prop)

    _, returned_ratio_wing, surface_wing_init = grid_to_surface_panels(grid_wing_init; ratios=ratio_wing, fcore=fcore)
    surfaces_prop_init = Vector{typeof(surface_wing_init)}()
    prop_surface_indices = [Vector{Int}(undef, Nb_prop) for _ in 1:Npropellers]
    for ip in 1:Npropellers
        for k in 1:Nb_prop
            push!(surfaces_prop_init, grid_to_surface_panels(grids_prop_initial_global[ip][k]; fcore=fcore)[3])
            prop_surface_indices[ip][k] = 1 + length(surfaces_prop_init)
        end
    end

    surfaces = [surface_wing_init, surfaces_prop_init...]
    nsurf = length(surfaces)

    surface_interaction_id = Vector{Int}(undef, nsurf)
    surface_interaction_id[1] = 1
    for ip in 1:Npropellers
        for k in 1:Nb_prop
            surface_interaction_id[prop_surface_indices[ip][k]] = ip + 1
        end
    end

    nwake_wing = Int(maximum_wake_rows_wing)
    nwake_prop = fill(Int(maximum_wake_rows_propeller), Npropellers * Nb_prop)
    nwake = vcat(nwake_wing, nwake_prop)
    system = System(surfaces; nw=nwake)
    verbose && println("UVLM system created.")

    system.reference[] = ref
    system.symmetric .= vcat(symmetric_wing, fill(false, Npropellers*Nb_prop))
    system.surface_id .= 1:nsurf
    system.wake_finite_core .= fill(true, nsurf)
    system.trailing_vortices .= false
    if isdefined(system, :ratios) && length(system.ratios) == nsurf
        system.ratios[1] = returned_ratio_wing
    end
    system.previous_surfaces[:] = surfaces[:]
    system.surfaces[:] = surfaces[:]

    initial_circulation = zero(system.Γ)
    initial_wakes = [Matrix{WakePanel{Float64}}(undef, 0, size(surfaces[i], 2)) for i = 1:nsurf]
    repeated_points = repeated_trailing_edge_points(surfaces)
    iwake = [min(size(initial_wakes[isurf], 1), nwake[isurf]) for isurf = 1:nsurf]
    fs_vec = isa(fs, Freestream) ? fill(fs, length(dt)) : fs
    system.freestream[] = fs_vec[1]
    system.Γ .= initial_circulation
    for isurf = 1:nsurf
        for I in CartesianIndices(initial_wakes[isurf])
            system.wakes[isurf][I] = initial_wakes[isurf][I]
        end
    end

    save = 1:length(dt)
    TF = eltype(system)
    surface_history = Vector{Vector{Matrix{SurfacePanel{TF}}}}(undef, length(save))
    nodal_forces_wing = Matrix{SVector{3,TF}}(undef, nc_wing+1, ns_wing+1)
    nodal_moments_wing = Matrix{SVector{3,TF}}(undef, nc_wing+1, ns_wing+1)
    EA_nodes_wing = Matrix{SVector{3,TF}}(undef, 1, nnodes)
    nodal_forces_prop = [[Matrix{SVector{3,TF}}(undef, nc_prop+1, ns_prop+1) for _ in 1:Nb_prop] for _ in 1:Npropellers]
    grids_prop_current = deepcopy(grids_prop_initial_global)
    T_pivot_A_current = Vector{SVector{3,Float64}}(undef, Npropellers)

    return (
        fcore = fcore,
        grid_wing_init = grid_wing_init,
        ratio_wing = ratio_wing,
        grids_prop_ref_single = grids_prop_ref_single,
        grids_prop_ref = grids_prop_ref,
        grids_prop_initial_global = grids_prop_initial_global,
        attach_node_y = attach_node_y,
        attach_node_chord = attach_node_chord,
        attach_node_xle = attach_node_xle,
        ea_x_aero = ea_x_aero,
        elastic_axis_fraction = elastic_axis_fraction,
        prop_pivot_offset_from_ea_A = prop_pivot_offset_from_ea_A,
        hub_center_prop_A = hub_center_prop_A,
        hub_offset_from_ea_A = hub_offset_from_ea_A,
        T_pivot_global_init = T_pivot_global_init,
        returned_ratio_wing = returned_ratio_wing,
        surface_wing_init = surface_wing_init,
        surfaces_prop_init = surfaces_prop_init,
        prop_surface_indices = prop_surface_indices,
        surfaces = surfaces,
        nsurf = nsurf,
        surface_interaction_id = surface_interaction_id,
        nwake_wing = nwake_wing,
        nwake_prop = nwake_prop,
        nwake = nwake,
        system = system,
        initial_circulation = initial_circulation,
        initial_wakes = initial_wakes,
        repeated_points = repeated_points,
        iwake = iwake,
        fs_vec = fs_vec,
        save = save,
        TF = TF,
        surface_history = surface_history,
        nodal_forces_wing = nodal_forces_wing,
        nodal_moments_wing = nodal_moments_wing,
        EA_nodes_wing = EA_nodes_wing,
        nodal_forces_prop = nodal_forces_prop,
        grids_prop_current = grids_prop_current,
        T_pivot_A_current = T_pivot_A_current,
    )
end
