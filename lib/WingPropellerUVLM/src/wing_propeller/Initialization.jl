# ==============================================================================
# Configuration-neutral UVLM system initialization for a wing with directly
# attached propellers. The aerodynamic hub station is always obtained from the
# selected structural node.
# ==============================================================================

using StaticArrays

function _per_propeller_scalar(value, count::Integer, name::AbstractString)
    if value isa Number
        return fill(value, count)
    end
    length(value) == count || throw(DimensionMismatch(
        "$name must be a scalar or contain one entry per propeller",
    ))
    return collect(value)
end

function _per_propeller_twists(value, count::Integer)
    if value isa AbstractVector{<:Real}
        return [collect(Float64, value) for _ in 1:count]
    end
    length(value) == count || throw(DimensionMismatch(
        "blade_twists_prop must be a shared vector or one vector per propeller",
    ))
    return [collect(Float64, twists) for twists in value]
end

function initialize_wing_propeller_uvlm_system(;
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
    maximum_wake_rows_propeller = 36,
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
    length(prop_attach_nodes) == Npropellers || throw(DimensionMismatch(
        "prop_attach_nodes must contain Npropellers entries",
    ))
    all(1 .<= prop_attach_nodes .<= length(span_nodes)) || throw(ArgumentError(
        "Every propeller attachment node must lie on the wing mesh",
    ))

    radii_prop = _per_propeller_scalar(R_prop, Npropellers, "R_prop")
    chords_prop = _per_propeller_scalar(c_prop, Npropellers, "c_prop")
    radial_panels_prop = Int.(_per_propeller_scalar(ns_prop, Npropellers, "ns_prop"))
    chordwise_panels_prop = Int.(_per_propeller_scalar(nc_prop, Npropellers, "nc_prop"))
    blade_counts_prop = Int.(_per_propeller_scalar(Nb_prop, Npropellers, "Nb_prop"))
    twists_prop = _per_propeller_twists(blade_twists_prop, Npropellers)
    maximum_wake_rows_props = Int.(_per_propeller_scalar(
        maximum_wake_rows_propeller,
        Npropellers,
        "maximum_wake_rows_propeller",
    ))
    all(radial_panels_prop .> 0) || throw(ArgumentError("ns_prop must be positive"))
    all(chordwise_panels_prop .> 0) || throw(ArgumentError("nc_prop must be positive"))
    all(blade_counts_prop .> 0) || throw(ArgumentError("Nb_prop must be positive"))
    all(maximum_wake_rows_props .>= 0) || throw(ArgumentError(
        "maximum_wake_rows_propeller must be nonnegative",
    ))

    # Each propeller may have independent radius, blade count, discretization,
    # and twist while using the same unchanged UVLM panel backend.
    grids_prop_ref = [
        generate_propeller_blades_grid(
            radii_prop[ip],
            chords_prop[ip],
            radial_panels_prop[ip],
            chordwise_panels_prop[ip],
            twists_prop[ip],
            blade_counts_prop[ip],
        ) for ip in 1:Npropellers
    ]
    grids_prop_ref_single = first(grids_prop_ref)
    grids_prop_initial_global = [deepcopy(grids_prop_ref[ip]) for ip in 1:Npropellers]

    attach_node_y = [span_nodes[prop_attach_nodes[ip]] for ip in 1:Npropellers]
    attach_node_chord = [chord[prop_attach_nodes[ip]] for ip in 1:Npropellers]
    attach_node_xle = [xle_distribution[prop_attach_nodes[ip]] for ip in 1:Npropellers]
    if !isnothing(propeller_span_positions)
        length(propeller_span_positions) == Npropellers || throw(DimensionMismatch(
            "propeller_span_positions must contain Npropellers entries",
        ))
        all(isapprox.(propeller_span_positions, attach_node_y; atol = 1e-12, rtol = 0.0)) ||
            throw(ArgumentError(
                "propeller_span_positions may only confirm, not override, attachment-node stations",
            ))
    end
    ea_x_aero = attach_node_xle .+ attach_node_chord .* elastic_axis_fraction

    hub_offset_from_ea_A = prop_pivot_offset_from_ea_A + hub_center_prop_A
    T_pivot_global_init = [SVector(ea_x_aero[ip], attach_node_y[ip], 0.0) + prop_pivot_offset_from_ea_A for ip in 1:Npropellers]
    initialize_propeller_grids!(grids_prop_initial_global, grids_prop_ref,
                                T_pivot_global_init, hub_center_prop_A, blade_counts_prop)

    _, returned_ratio_wing, surface_wing_init = grid_to_surface_panels(grid_wing_init; ratios=ratio_wing, fcore=fcore)
    surfaces_prop_init = Vector{typeof(surface_wing_init)}()
    prop_surface_indices = [Vector{Int}(undef, blade_counts_prop[ip]) for ip in 1:Npropellers]
    for ip in 1:Npropellers
        for k in 1:blade_counts_prop[ip]
            push!(surfaces_prop_init, grid_to_surface_panels(grids_prop_initial_global[ip][k]; fcore=fcore)[3])
            prop_surface_indices[ip][k] = 1 + length(surfaces_prop_init)
        end
    end

    surfaces = [surface_wing_init, surfaces_prop_init...]
    nsurf = length(surfaces)

    surface_interaction_id = Vector{Int}(undef, nsurf)
    surface_interaction_id[1] = 1
    for ip in 1:Npropellers
        for k in 1:blade_counts_prop[ip]
            surface_interaction_id[prop_surface_indices[ip][k]] = ip + 1
        end
    end

    nwake_wing = Int(maximum_wake_rows_wing)
    nwake_prop = reduce(vcat, [
        fill(maximum_wake_rows_props[ip], blade_counts_prop[ip])
        for ip in 1:Npropellers
    ]; init = Int[])
    nwake = vcat(nwake_wing, nwake_prop)
    system = System(surfaces; nw=nwake)
    verbose && println("UVLM system created.")

    system.reference[] = ref
    system.symmetric .= vcat(symmetric_wing, fill(false, sum(blade_counts_prop)))
    verbose && println("Wing image symmetry across y=0: $symmetric_wing; propeller blade symmetry: false")
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
    nodal_forces_prop = [[
        Matrix{SVector{3,TF}}(
            undef,
            chordwise_panels_prop[ip] + 1,
            radial_panels_prop[ip] + 1,
        ) for _ in 1:blade_counts_prop[ip]
    ] for ip in 1:Npropellers]
    grids_prop_current = deepcopy(grids_prop_initial_global)
    T_pivot_A_current = Vector{SVector{3,Float64}}(undef, Npropellers)

    return (
        fcore = fcore,
        grid_wing_init = grid_wing_init,
        ratio_wing = ratio_wing,
        grids_prop_ref_single = grids_prop_ref_single,
        grids_prop_ref = grids_prop_ref,
        radii_prop = radii_prop,
        chords_prop = chords_prop,
        radial_panels_prop = radial_panels_prop,
        chordwise_panels_prop = chordwise_panels_prop,
        blade_counts_prop = blade_counts_prop,
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

"""Compatibility alias for the former configuration-specific entry point."""
initialize_bohnisch_uvlm_system(; kwargs...) =
    initialize_wing_propeller_uvlm_system(; kwargs...)
