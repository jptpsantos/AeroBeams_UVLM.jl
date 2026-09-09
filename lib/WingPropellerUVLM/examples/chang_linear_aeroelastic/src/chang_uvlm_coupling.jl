# Chang linear-structure to UVLM adapter.
#
# These functions intentionally remain with the example because they encode
# its structural DOF order and the Z-down structural/Z-up aerodynamic sign
# conversion. Generic grid deformation, UVLM stepping, segment-to-nodal loads,
# and generalized-alpha integration live in WingPropellerUVLM itself.

"""Initialize the case UVLM system from `chang_aerodynamic_options`."""
function initialize_chang_uvlm(options)
    return initialize_chang_uvlm(;
        finite_core = options.finite_core,
        elastic_axis_fraction = options.elastic_axis_fraction,
        propeller_pivot_offset_A = options.propeller_pivot_offset_A,
        physical_hub_center_A = options.physical_hub_center_A,
        maximum_wake_rows_wing = options.maximum_wake_rows_wing,
        maximum_wake_rows_propeller = options.maximum_wake_rows_propeller,
    )
end

"""Collect the wake controls while sharing the initialized mutable storage."""
function chang_wake_context(uvlm; interaction_on::Bool)
    return (;
        surface_count = uvlm.nsurf,
        repeated_points = uvlm.repeated_points,
        maximum_rows = uvlm.nwake,
        active_rows = uvlm.iwake,
        interaction_ids = uvlm.surface_interaction_id,
        interaction_on,
        saved_steps = uvlm.save,
        surface_history = uvlm.surface_history,
    )
end

"""
    initialize_chang_uvlm(; kwargs...)

Create the coupled wing--propeller UVLM state from the geometry and operating
point loaded by `chang_model_parameters.jl`. The entry-point file supplies only
the numerical regularization, reference locations, and retained-wake limits.
"""
function initialize_chang_uvlm(;
    finite_core,
    elastic_axis_fraction,
    propeller_pivot_offset_A,
    physical_hub_center_A,
    maximum_wake_rows_wing,
    maximum_wake_rows_propeller,
)
    return initialize_bohnisch_uvlm_system(
        xle = xle,
        yle = yle,
        zle = zle,
        chord_geo = chord_geo,
        theta_geo = theta_geo,
        phi_geo = phi_geo,
        ns_wing = ns_wing,
        nc_wing = nc_wing,
        mirror_wing = mirror_wing,
        spacing_s_wing = spacing_s_wing,
        spacing_c_wing = spacing_c_wing,
        R_prop = R_prop,
        c_prop = c_prop,
        ns_prop = ns_prop,
        nc_prop = nc_prop,
        blade_twists_prop = blade_twists_prop,
        Nb_prop = Nb_prop,
        Npropellers = Npropellers,
        span_nodes = span_nodes,
        prop_attach_nodes = prop_attach_nodes,
        propeller_span_positions = propeller_span_positions,
        chord = chord,
        xle_distribution = xle_distribution,
        ref = ref,
        symmetric_wing = symmetric_wing,
        fs = fs,
        dt = dt,
        nnodes = nnodes,
        prop_pivot_offset_from_ea_A = propeller_pivot_offset_A,
        hub_center_prop_A = physical_hub_center_A,
        fcore = finite_core,
        elastic_axis_fraction = elastic_axis_fraction,
        maximum_wake_rows_wing = maximum_wake_rows_wing,
        maximum_wake_rows_propeller = maximum_wake_rows_propeller,
        verbose = true,
    )
end

"""
    update_aero_geometry_for_state!(system, q_free, time_np1)

Map one structural displacement guess to the UVLM geometry at `time_np1`.

The Chang structural basis is `[span, chord, down]`; the UVLM aerodynamic
basis is `[chord, span, up]`. The function reconstructs the clamped root,
applies all six wing DOFs, adds each propeller's pitch/yaw modal motion and
prescribed spin, rebuilds the panel surfaces, and returns the aerodynamic-frame
wing kinematics needed to form moment arms during load transfer.

This function changes geometry only. It does not solve circulation or advance
the wake.
"""

function update_aero_geometry_for_state!(system, q_free::AbstractVector, time_np1::Real)
    # Split the global free-state vector according to the assembly convention:
    # all wing free DOFs first, then [pitch, yaw] for each propeller.
    q_wing_free = q_free[1:ndof_wing_free]
    q_propeller_free = q_free[(ndof_wing_free + 1):end]

    # Reinsert the clamped root DOFs and extract each component along the beam.
    displacement_span_structural = vcat(0.0, q_wing_free[1:ndof:end])
    displacement_chord_structural = vcat(0.0, q_wing_free[2:ndof:end])
    displacement_vertical_structural = vcat(0.0, q_wing_free[3:ndof:end])
    rotation_span_structural = vcat(0.0, q_wing_free[4:ndof:end])
    rotation_chord_structural = vcat(0.0, q_wing_free[5:ndof:end])
    rotation_vertical_structural = vcat(0.0, q_wing_free[6:ndof:end])

    # Proper structural-to-aerodynamic basis mapping:
    # (span, chord, down) -> (y, x, -z), with the same mapping for rotations.
    u_x_aero = displacement_chord_structural
    u_y_aero = displacement_span_structural
    u_z_aero = -displacement_vertical_structural
    theta_x_aero = rotation_chord_structural
    theta_y_aero = rotation_span_structural
    theta_z_aero = -rotation_vertical_structural

    # Interpolate nodal translations/rotations spanwise and move every wing
    # lattice point about the elastic axis.
    grid_wing = generate_panel_grid_and_interpolate(
        span_length,
        chord,
        xle_distribution,
        ns_wing,
        nc_wing,
        u_x_aero,
        u_y_aero,
        u_z_aero,
        theta_x_aero,
        theta_y_aero,
        theta_z_aero;
        elastic_axis_fraction = elastic_axis_fraction,
    )

    # Propeller motion is composed in this order: prescribed rotor spin,
    # pylon pitch/yaw, then the interpolated attachment rotation/translation.
    for propeller_index in 1:Npropellers
        left_node, right_node = prop_attachment_node_pairs[propeller_index]
        left_weight, right_weight = prop_attachment_weights[propeller_index]
        pitch_aero = q_propeller_free[2 * (propeller_index - 1) + 1]
        yaw_aero = -q_propeller_free[2 * (propeller_index - 1) + 2]

        u_x_attachment = left_weight * u_x_aero[left_node] +
            right_weight * u_x_aero[right_node]
        u_y_attachment = left_weight * u_y_aero[left_node] +
            right_weight * u_y_aero[right_node]
        u_z_attachment = left_weight * u_z_aero[left_node] +
            right_weight * u_z_aero[right_node]
        theta_x_attachment = left_weight * theta_x_aero[left_node] +
            right_weight * theta_x_aero[right_node]
        theta_y_attachment = left_weight * theta_y_aero[left_node] +
            right_weight * theta_y_aero[right_node]
        theta_z_attachment = left_weight * theta_z_aero[left_node] +
            right_weight * theta_z_aero[right_node]

        pivot = SVector(
            ea_x_aero[propeller_index] + u_x_attachment,
            attach_node_y[propeller_index] + u_y_attachment,
            u_z_attachment,
        )
        T_pivot_A_current[propeller_index] = pivot

        wing_rotation = RotationMatrix(theta_z_attachment, 3) *
            RotationMatrix(theta_x_attachment, 1) *
            RotationMatrix(theta_y_attachment, 2)
        pitch_rotation = RotationMatrix(pitch_aero, 2)
        yaw_rotation = RotationMatrix(yaw_aero, 3)
        whirl_rotation = pitch_rotation * yaw_rotation
        spin_rotation = RotationMatrix(-Ω * time_np1, 1)

        # Work-conjugate axes for R_whirl = R_y(pitch) R_z(-yaw). The pitch
        # axis is the wing-rotated y axis. The structural yaw coordinate is
        # opposite aerodynamic z and its instantaneous axis is also carried
        # by the preceding pitch rotation.
        pitch_axis_A_current[propeller_index] = wing_rotation *
            SVector(0.0, 1.0, 0.0)
        yaw_axis_A_current[propeller_index] = -wing_rotation * pitch_rotation *
            SVector(0.0, 0.0, 1.0)

        T_hub_A_current[propeller_index] = pivot + wing_rotation *
            (hub_center_prop_A + (whirl_rotation * hub_center_load_A - hub_center_load_A))
        T_load_A_current[propeller_index] = pivot + wing_rotation *
            (whirl_rotation * hub_center_load_A)

        for blade_index in 1:Nb_prop
            reference_grid = grids_prop_ref[propeller_index][blade_index]
            current_grid = grids_prop_current[propeller_index][blade_index]
            for chord_index in axes(reference_grid, 2), radial_index in axes(reference_grid, 3)
                reference_point = SVector{3}(reference_grid[:, chord_index, radial_index])
                spun_point = spin_rotation * reference_point
                static_point = spun_point + hub_center_prop_A
                modal_point = spun_point + hub_center_load_A
                whirled_point = static_point + (whirl_rotation * modal_point - modal_point)
                current_grid[:, chord_index, radial_index] =
                    pivot + wing_rotation * whirled_point
            end
        end
    end

    # Replace the panel geometry in the existing System. Keeping the System
    # object itself preserves all preallocated circulation, force, and wake
    # storage used by `propagate_system!`.
    _, _, wing_surface = grid_to_surface_panels(
        grid_wing;
        ratios = ratio_wing,
        fcore = FCORE,
    )
    propeller_surfaces = Vector{typeof(wing_surface)}()
    for propeller_index in 1:Npropellers, blade_index in 1:Nb_prop
        push!(
            propeller_surfaces,
            grid_to_surface_panels(
                grids_prop_current[propeller_index][blade_index];
                fcore = FCORE,
            )[3],
        )
    end

    system.surfaces[1] = wing_surface
    for surface_index in eachindex(propeller_surfaces)
        system.surfaces[surface_index + 1] = propeller_surfaces[surface_index]
    end

    return (;
        u_x_A = u_x_aero,
        u_y_A = u_y_aero,
        u_z_A = u_z_aero,
        theta_x_A = theta_x_aero,
        theta_y_A = theta_y_aero,
        theta_z_A = theta_z_aero,
        propeller_pitch_axes_A = pitch_axis_A_current,
        propeller_yaw_axes_A = yaw_axis_A_current,
    )
end

"""
    chang_wing_generalized_moment(moment_A, theta_x_A, theta_z_A)

Project a spatial moment in the aerodynamic basis onto the wing's structural
Euler coordinates `[span, chord, down]`. Geometry uses
`R = Rz(theta_z_A) Rx(theta_x_A) Ry(theta_y_A)`, with structural down rotation
equal to `-theta_z_A`. The work-conjugate axes are `Rz Rx e_y`, `Rz e_x`, and
`-e_z`; they do not depend on the innermost angle `theta_y_A`.
"""
function chang_wing_generalized_moment(moment_A, theta_x_A, theta_z_A)
    sx, cx = sincos(theta_x_A)
    sz, cz = sincos(theta_z_A)
    return SVector(
        -sz * cx * moment_A[1] + cz * cx * moment_A[2] + sx * moment_A[3],
        cz * moment_A[1] + sz * moment_A[2],
        -moment_A[3],
    )
end

"""
    assemble_structural_aero_load!(system, kinematics; step=0, print_loads=false)

Transfer the dimensional UVLM vertex forces to the Chang free-DOF
ordering. Wing forces are summed chordwise and moments are formed about the
deformed elastic axis and projected onto the instantaneous wing rotation axes.
Blade forces are reduced to propeller hub/pivot wrenches. Their wing moments
are projected at the interpolated attachment angles before distribution over
the two attachment nodes; modal moments use the selected pitch/yaw projection.

Returns one generalized-load vector ordered exactly like the reduced
structural state used by `M`, `C`, and `K`.
"""
function assemble_structural_aero_load!(system, kinematics;
    step::Int = 0, print_loads::Bool = false)

    zero_vector = SVector{3,TF}(0.0, 0.0, 0.0)
    fill!(nodal_forces_wing, zero_vector)
    fill!(nodal_moments_wing, zero_vector)

    # These arrays are dimensional forces and their matching vortex-lattice
    # vertex positions. Using the paired positions preserves moment arms.
    surface_forces_aero = imperial_nodal_forces(system)
    surface_positions_aero = imperial_nodal_positions(system)
    nodal_forces_wing .= surface_forces_aero[1]

    # Build the deformed elastic-axis line used as the wing moment reference.
    for node_index in 1:nnodes
        elastic_axis_x = xle_distribution[node_index] +
            chord[node_index] * elastic_axis_fraction + kinematics.u_x_A[node_index]
        EA_nodes_wing[node_index] = SVector(
            elastic_axis_x,
            span_nodes[node_index] + kinematics.u_y_A[node_index],
            kinematics.u_z_A[node_index],
        )
    end
    for chord_index in 1:(nc_wing + 1), span_index in 1:(ns_wing + 1)
        nodal_moments_wing[chord_index, span_index] = cross(
            surface_positions_aero[1][chord_index, span_index] - EA_nodes_wing[span_index],
            nodal_forces_wing[chord_index, span_index],
        )
    end

    wing_loads = zeros(NDOF)
    force_x = sum(getindex.(nodal_forces_wing, 1), dims = 1)
    force_y = sum(getindex.(nodal_forces_wing, 2), dims = 1)
    force_z = sum(getindex.(nodal_forces_wing, 3), dims = 1)
    moment_x = sum(getindex.(nodal_moments_wing, 1), dims = 1)
    moment_y = sum(getindex.(nodal_moments_wing, 2), dims = 1)
    moment_z = sum(getindex.(nodal_moments_wing, 3), dims = 1)

    for node_index in 1:nnodes
        offset = ndof * (node_index - 1)
        generalized_moment = chang_wing_generalized_moment(
            SVector(moment_x[node_index], moment_y[node_index], moment_z[node_index]),
            kinematics.theta_x_A[node_index],
            kinematics.theta_z_A[node_index],
        )
        wing_loads[offset + 1] = force_y[node_index]
        wing_loads[offset + 2] = force_x[node_index]
        wing_loads[offset + 3] = -force_z[node_index]
        wing_loads[offset + 4] = generalized_moment[1]
        wing_loads[offset + 5] = generalized_moment[2]
        wing_loads[offset + 6] = generalized_moment[3]
    end

    propeller_loads = zeros(ndof_P)
    # Reduce all blade vertex forces of each propeller to one resultant force
    # and moment, then apply them to both the pylon modal and wing attachment
    # coordinates without discarding the appropriate lever arms.
    for propeller_index in 1:Npropellers
        total_force = zero_vector
        total_moment_about_hub = zero_vector
        pivot = T_pivot_A_current[propeller_index]
        hub = T_hub_A_current[propeller_index]
        modal_load_point = T_load_A_current[propeller_index]

        for blade_index in 1:Nb_prop
            surface_index = prop_surface_indices[propeller_index][blade_index]
            blade_forces = nodal_forces_prop[propeller_index][blade_index]
            blade_forces .= surface_forces_aero[surface_index]
            blade_positions = surface_positions_aero[surface_index]

            for chord_index in axes(blade_forces, 1), radial_index in axes(blade_forces, 2)
                force = SVector{3,TF}(blade_forces[chord_index, radial_index])
                position = blade_positions[chord_index, radial_index]
                total_force += force
                total_moment_about_hub += cross(position - hub, force)
            end
        end

        modal_moment = total_moment_about_hub + cross(modal_load_point - pivot, total_force)
        wing_moment = total_moment_about_hub + cross(hub - pivot, total_force)

        if AEROELASTIC_PROPELLER_MOMENT_PROJECTION == :exact_virtual_work
            pitch_axis = kinematics.propeller_pitch_axes_A[propeller_index]
            yaw_axis = kinematics.propeller_yaw_axes_A[propeller_index]
            propeller_loads[2 * (propeller_index - 1) + 1] =
                dot(modal_moment, pitch_axis)
            propeller_loads[2 * (propeller_index - 1) + 2] =
                dot(modal_moment, yaw_axis)
        else
            # Original small-angle projection onto fixed aerodynamic axes.
            propeller_loads[2 * (propeller_index - 1) + 1] = modal_moment[2]
            propeller_loads[2 * (propeller_index - 1) + 2] = -modal_moment[3]
        end

        # The geometry interpolates Euler angles before rotating the propeller.
        # Apply that same chain rule: project at the attachment, then multiply
        # by the nodal weights. Individual node axes would be inconsistent when
        # the neighboring rotations differ.
        left_node, right_node = prop_attachment_node_pairs[propeller_index]
        left_weight, right_weight = prop_attachment_weights[propeller_index]
        theta_x_attachment = left_weight * kinematics.theta_x_A[left_node] +
            right_weight * kinematics.theta_x_A[right_node]
        theta_z_attachment = left_weight * kinematics.theta_z_A[left_node] +
            right_weight * kinematics.theta_z_A[right_node]
        generalized_wing_moment = chang_wing_generalized_moment(
            wing_moment,
            theta_x_attachment,
            theta_z_attachment,
        )
        structural_wrench = (
            total_force[2],
            total_force[1],
            -total_force[3],
            generalized_wing_moment[1],
            generalized_wing_moment[2],
            generalized_wing_moment[3],
        )
        for (node, weight) in (
            (left_node, left_weight),
            (right_node, right_weight),
        )
            node_offset = ndof * (node - 1)
            for local_dof in 1:ndof
                wing_loads[node_offset + local_dof] +=
                    weight * structural_wrench[local_dof]
            end
        end

        if print_loads
            println("\n--- Step $step Propeller P$(propeller_index) (partitioned GA) ---")
            println("  total_force_A            = $total_force")
            println("  total_moment_about_hub_A = $total_moment_about_hub")
            println("  total_moment_modal_A     = $modal_moment")
            println("  total_moment_wing_A      = $wing_moment")
            println(
                "  propeller moment projection = " *
                "$AEROELASTIC_PROPELLER_MOMENT_PROJECTION",
            )
        end
    end

    free_wing_loads = wing_loads[free_dofs[1:ndof_wing_free]]
    return vcat(free_wing_loads, propeller_loads)
end

"""
    aero_load_for_state!(system, snapshot, state, step; print_loads=false)

Evaluate the complete aerodynamic generalized-load operator for one structural
fixed-point guess at `t[step+1]`.

Every call restores the beginning-of-step snapshot, updates geometry, solves
one UVLM circulation/load trial, and transfers its loads. Wake convection and
shedding are deliberately disabled here. Therefore repeated calls during the
same coupled step are deterministic trials from `t[step]`, not successive wake
steps. On return, `system` contains the circulation/load trial corresponding to
`state`; the run driver advances its wake only after convergence.
"""
function aero_load_for_state!(system, snapshot, state::AbstractVector, step::Int;
    print_loads::Bool = false)

    # Transaction rollback: discard the previous coupling iterate's wake,
    # circulation, geometry, and surface-history changes.
    restore_uvlm!(system, snapshot)
    # Geometry is evaluated at the end of the physical step, including the
    # prescribed rotor azimuth at that time.
    kinematics = update_aero_geometry_for_state!(system, state, t[step + 1])
    # Solve the unsteady aerodynamic trial through circulation, gamma-dot, and
    # selected near-field segment loads. The expensive free-wake convection is a
    # physical-time update, so it is deferred until this structural state has
    # converged and is accepted by the run driver.
    propagate_system!(
        system,
        fs_vec[step],
        dt[step];
        additional_velocity = nothing,
        repeated_points = repeated_points,
        nwake = iwake,
        eta = 0.1,
        calculate_influence_matrix = true,
        near_field_analysis = true,
        derivatives = false,
        near_field_force_function = AEROELASTIC_NEAR_FIELD_FORCE_FUNCTION,
        interaction_id = surface_interaction_id,
        interaction = INTERACTION_ON,
        advance_wake = false,
    )
    # Convert the aerodynamic trial into a load vector for the structural
    # corrector. The caller decides whether to subtract the trim baseline.
    return assemble_structural_aero_load!(
        system,
        kinematics;
        step = step,
        print_loads = print_loads,
    )
end
