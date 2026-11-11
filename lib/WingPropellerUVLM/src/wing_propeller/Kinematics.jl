# ==============================================================================
# Kinematics helpers for the general wing-propeller coupled system.
# ==============================================================================

using LinearAlgebra
using StaticArrays
using FLOWMath

"""Store the accepted geometry used by UVLM to calculate surface velocities."""
function copy_surfaces_to_previous!(system,nsurf::Int)
    for i = 1:nsurf
        system.previous_surfaces[i] .= system.surfaces[i]
    end
    return nothing
end

"""Map the reduced Bohnisch wing coordinates to aerodynamic translations and rotations."""
function wing_kinematics_from_free_state(q_W_free,ndof::Int,nnodes::Int)
    h_down_S = vcat(0.0,q_W_free[1:ndof:end])
    slope_x_S = vcat(0.0,q_W_free[2:ndof:end])
    torsion_y_S = vcat(0.0,q_W_free[3:ndof:end])

    u_x_A = zeros(nnodes)
    u_y_A = zeros(nnodes)
    u_z_A = -h_down_S
    theta_x_A = -slope_x_S
    theta_y_A = torsion_y_S
    theta_z_A = zeros(nnodes)

    return (
        h_down_S=h_down_S,
        slope_x_S=slope_x_S,
        torsion_y_S=torsion_y_S,
        u_x_A=u_x_A,
        u_y_A=u_y_A,
        u_z_A=u_z_A,
        theta_x_A=theta_x_A,
        theta_y_A=theta_y_A,
        theta_z_A=theta_z_A
    )
end

# Position vectors use r; p denotes AeroBeams rotation parameters.
"""Place the undeformed propeller grids at their initial global pivots."""
function initialize_propeller_grids!(grids_prop_initial_global,grids_prop_ref,
    T_pivot_global_init,hub_center_prop_A,blade_counts_prop)

    R_whirl_init = RotationMatrix(0.0,2) * RotationMatrix(0.0,3)
    R_spin_init = I(3)

    for ip in eachindex(grids_prop_ref)
        for k in eachindex(grids_prop_ref[ip])
            gridRefK = grids_prop_ref[ip][k]
            gridGlobalK = grids_prop_initial_global[ip][k]
            for i = 1:size(gridRefK,2),j = 1:size(gridRefK,3)
                r_ref_A = SVector{3}(gridRefK[:,i,j])
                r_spun_A = R_spin_init * r_ref_A
                r_pivot_rel_A = r_spun_A + hub_center_prop_A
                r_whirled_A = R_whirl_init * r_pivot_rel_A
                gridGlobalK[:,i,j] = T_pivot_global_init[ip] + r_whirled_A
            end
        end
    end

    return grids_prop_initial_global
end

"""Apply attachment deformation, pylon whirl, and prescribed spin at `t_next`."""
function update_propeller_grids!(grids_prop_current,T_pivot_A_current,grids_prop_ref,
    q_P_free,prop_attach_nodes,ea_x_aero,attach_node_y,
    u_z_A,theta_x_A,theta_y_A,theta_z_A,
    prop_pivot_offset_from_ea_A,hub_center_prop_A,
    Omega_prop,t_next,Npropellers::Int,blade_counts_prop)

    for ip = 1:Npropellers
        pos = prop_attach_nodes[ip]
        prop_pitch_A = q_P_free[2*(ip-1)+1]
        prop_yaw_A = -q_P_free[2*(ip-1)+2]

        # Wing motion carries the propeller pivot and its local rotation axes.
        R_A_node =
            RotationMatrix(theta_z_A[pos],3) * RotationMatrix(theta_x_A[pos],1) * RotationMatrix(theta_y_A[pos],2)
        T_elastic_axis_A = SVector(ea_x_aero[ip],attach_node_y[ip],u_z_A[pos])
        T_pivot_A = T_elastic_axis_A + R_A_node * prop_pivot_offset_from_ea_A
        T_pivot_A_current[ip] = T_pivot_A

        # Rotor azimuth is evaluated from physical time, so repeated coupling
        # trials at the same t[n+1] do not advance it more than once.
        R_whirl_A = RotationMatrix(prop_pitch_A,2) * RotationMatrix(prop_yaw_A,3)
        omega = Omega_prop isa Number ? Omega_prop : Omega_prop[ip]
        R_spin_A = RotationMatrix(-omega * t_next,1)

        for k in eachindex(grids_prop_ref[ip])
            gridRefBlade = grids_prop_ref[ip][k]
            currentGridBlade = grids_prop_current[ip][k]
            for i = 1:size(gridRefBlade,2),j = 1:size(gridRefBlade,3)
                r_ref_A = SVector{3}(gridRefBlade[:,i,j])
                r_spun_A = R_spin_A * r_ref_A
                r_pivot_rel_A = r_spun_A + hub_center_prop_A
                r_whirled_A = R_whirl_A * r_pivot_rel_A
                r_deformed_A = R_A_node * r_whirled_A
                currentGridBlade[:,i,j] = T_pivot_A + r_deformed_A
            end
        end
    end

    return grids_prop_current
end

"""Convert deformed grids to panels and replace only the System surface geometry."""
function update_system_surfaces!(system,gridWing,grids_prop_current,
    ratio_wing,Npropellers::Int,blade_counts_prop;
    fcore=(c,ds) -> 1e-3)

    _,_,currentSurfaceWing = grid_to_surface_panels(gridWing; ratios=ratio_wing,fcore=fcore)
    currentSurfacesProp = Vector{typeof(currentSurfaceWing)}()

    for ip = 1:Npropellers
        for k in eachindex(grids_prop_current[ip])
            push!(currentSurfacesProp,grid_to_surface_panels(grids_prop_current[ip][k]; fcore=fcore)[3])
        end
    end

    system.surfaces[1] = currentSurfaceWing
    for k = 1:length(currentSurfacesProp)
        system.surfaces[k+1] = currentSurfacesProp[k]
    end

    return currentSurfaceWing,currentSurfacesProp
end
