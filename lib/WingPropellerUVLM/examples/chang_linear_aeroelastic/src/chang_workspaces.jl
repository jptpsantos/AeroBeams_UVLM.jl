# Shared workspaces for the Chang geometry/load adapter.
#
# Include after structural, aerodynamic_options, and uvlm have been created.
# These are aliases to the initialized arrays, not copies. Keeping the names
# also lets the standalone validation scripts inspect the accepted UVLM state.

free_dofs = structural.free_dofs
ndof_wing_free = structural.ndof_wing_free
FCORE = aerodynamic_options.finite_core
elastic_axis_fraction = aerodynamic_options.elastic_axis_fraction
prop_pivot_offset_from_ea_A = aerodynamic_options.propeller_pivot_offset_A
hub_center_prop_A = aerodynamic_options.physical_hub_center_A
hub_center_load_A = aerodynamic_options.load_center_A

# Geometry, surface groups, and wake storage.
(;
    ratio_wing, grids_prop_ref, attach_node_y, ea_x_aero, prop_surface_indices,
    nsurf, surface_interaction_id, nwake, system, repeated_points, iwake,
    fs_vec, save, TF, surface_history,
    nodal_forces_wing, nodal_moments_wing, EA_nodes_wing,
    nodal_forces_prop, grids_prop_current, T_pivot_A_current,
) = uvlm

# Current hub locations and instantaneous axes, updated for every state trial.
T_hub_A_current = Vector{SVector{3,Float64}}(undef, Npropellers)
T_load_A_current = Vector{SVector{3,Float64}}(undef, Npropellers)
pitch_axis_A_current = Vector{SVector{3,Float64}}(undef, Npropellers)
yaw_axis_A_current = Vector{SVector{3,Float64}}(undef, Npropellers)
