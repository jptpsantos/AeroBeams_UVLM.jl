# ==============================================================================
# Kinematics helpers for the Bohnisch wing-propeller coupled system.
# ==============================================================================

using LinearAlgebra
using StaticArrays
using FLOWMath

function copy_surfaces_to_previous!(system, nsurf::Int)
    for i = 1:nsurf
        system.previous_surfaces[i] .= system.surfaces[i]
    end
    return nothing
end

function wing_kinematics_from_free_state(q_W_free, ndof::Int, nnodes::Int)
    h_down_S = vcat(0.0, q_W_free[1:ndof:end])
    slope_x_S = vcat(0.0, q_W_free[2:ndof:end])
    torsion_y_S = vcat(0.0, q_W_free[3:ndof:end])

    u_x_A = zeros(nnodes)
    u_y_A = zeros(nnodes)
    u_z_A = -h_down_S
    theta_x_A = -slope_x_S
    theta_y_A = torsion_y_S
    theta_z_A = zeros(nnodes)

    return (
        h_down_S = h_down_S,
        slope_x_S = slope_x_S,
        torsion_y_S = torsion_y_S,
        u_x_A = u_x_A,
        u_y_A = u_y_A,
        u_z_A = u_z_A,
        theta_x_A = theta_x_A,
        theta_y_A = theta_y_A,
        theta_z_A = theta_z_A,
    )
end

function initialize_propeller_grids!(grids_prop_initial_global, grids_prop_ref,
    T_pivot_global_init, hub_center_prop_A, Nb_prop::Int)

    R_whirl_init = RotationMatrix(0.0, 2) * RotationMatrix(0.0, 3)
    R_spin_init = I(3)

    for ip in eachindex(grids_prop_ref)
        for k in 1:Nb_prop
            grid_ref_k = grids_prop_ref[ip][k]
            grid_global_k = grids_prop_initial_global[ip][k]
            for i in 1:size(grid_ref_k, 2), j in 1:size(grid_ref_k, 3)
                p_ref_A = SVector{3}(grid_ref_k[:, i, j])
                p_spun_A = R_spin_init * p_ref_A
                p_pivot_rel_A = p_spun_A + hub_center_prop_A
                p_whirled_A = R_whirl_init * p_pivot_rel_A
                grid_global_k[:, i, j] = T_pivot_global_init[ip] + p_whirled_A
            end
        end
    end

    return grids_prop_initial_global
end

function update_propeller_grids!(grids_prop_current, T_pivot_A_current, grids_prop_ref,
    q_P_free, prop_attach_nodes, ea_x_aero, attach_node_y,
    u_z_A, theta_x_A, theta_y_A, theta_z_A,
    prop_pivot_offset_from_ea_A, hub_center_prop_A,
    Omega_prop, t_next, Npropellers::Int, Nb_prop::Int)

    for ip in 1:Npropellers
        pos = prop_attach_nodes[ip]
        prop_pitch_A = q_P_free[2*(ip-1)+1]
        prop_yaw_A = -q_P_free[2*(ip-1)+2]

        R_A_node = RotationMatrix(theta_z_A[pos], 3) * RotationMatrix(theta_x_A[pos], 1) * RotationMatrix(theta_y_A[pos], 2)
        T_elastic_axis_A = SVector(ea_x_aero[ip], attach_node_y[ip], u_z_A[pos])
        T_pivot_A = T_elastic_axis_A + R_A_node * prop_pivot_offset_from_ea_A
        T_pivot_A_current[ip] = T_pivot_A

        R_whirl_A = RotationMatrix(prop_pitch_A, 2) * RotationMatrix(prop_yaw_A, 3)
        R_spin_A = RotationMatrix(-Omega_prop * t_next, 1)

        for k in 1:Nb_prop
            grid_ref_blade = grids_prop_ref[ip][k]
            current_grid_blade = grids_prop_current[ip][k]
            for i in 1:size(grid_ref_blade, 2), j in 1:size(grid_ref_blade, 3)
                p_ref_A = SVector{3}(grid_ref_blade[:, i, j])
                p_spun_A = R_spin_A * p_ref_A
                p_pivot_rel_A = p_spun_A + hub_center_prop_A
                p_whirled_A = R_whirl_A * p_pivot_rel_A
                p_deformed_A = R_A_node * p_whirled_A
                current_grid_blade[:, i, j] = T_pivot_A + p_deformed_A
            end
        end
    end

    return grids_prop_current
end

function update_system_surfaces!(system, grid_wing, grids_prop_current,
    ratio_wing, Npropellers::Int, Nb_prop::Int;
    fcore = (c, ds) -> 1e-3)

    _, _, current_surface_wing = grid_to_surface_panels(grid_wing; ratios=ratio_wing, fcore=fcore)
    current_surfaces_prop = Vector{typeof(current_surface_wing)}()

    for ip in 1:Npropellers
        for k in 1:Nb_prop
            push!(current_surfaces_prop, grid_to_surface_panels(grids_prop_current[ip][k]; fcore=fcore)[3])
        end
    end

    system.surfaces[1] = current_surface_wing
    for k in 1:length(current_surfaces_prop)
        system.surfaces[k+1] = current_surfaces_prop[k]
    end

    return current_surface_wing, current_surfaces_prop
end

