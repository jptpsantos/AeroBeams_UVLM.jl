# Configuration-independent adapter between the full six-DOF beam and the
# unchanged UVLM backend. Wing surface motion/loads may be interpolated along
# beam elements; each propeller remains attached to exactly one beam node.

function _general_wing_geometry(case, stations)
    geometry = case.wing.geometry
    span_nodes = case.wing.span_nodes
    span = last(span_nodes) - first(span_nodes)
    fraction = (stations .- first(span_nodes)) ./ span
    chord = geometry.root_chord .+
        (geometry.tip_chord - geometry.root_chord) .* fraction
    xle = geometry.xle_root .+
        (geometry.xle_tip - geometry.xle_root) .* fraction
    return chord, xle
end

function _beam_interpolation_support(span_nodes, position)
    position <= first(span_nodes) && return ((1, 1.0),)
    position >= last(span_nodes) && return ((length(span_nodes), 1.0),)
    right = searchsortedfirst(span_nodes, position)
    isapprox(position, span_nodes[right]; atol = 1.0e-12, rtol = 0.0) &&
        return ((right, 1.0),)
    left = right - 1
    right_weight = (position - span_nodes[left]) /
        (span_nodes[right] - span_nodes[left])
    return ((left, 1 - right_weight), (right, right_weight))
end

"""Work-conjugate projection for the implemented wing Euler-angle sequence."""
function _generalized_wing_moment(moment, theta_x, theta_z)
    sine_x, cosine_x = sincos(theta_x)
    sine_z, cosine_z = sincos(theta_z)
    return SVector(
        -sine_z * cosine_x * moment[1] +
            cosine_z * cosine_x * moment[2] + sine_x * moment[3],
        cosine_z * moment[1] + sine_z * moment[2],
        -moment[3],
    )
end

function _update_general_uvlm_geometry!(case, structural, uvlm, state, time, finite_core)
    wing_kinematics = full_beam_aerodynamic_kinematics(structural, state)
    span_nodes = collect(Float64, case.wing.span_nodes)
    aero_span = collect(range(
        first(span_nodes),
        last(span_nodes);
        length = case.aerodynamic.spanwise_panels + 1,
    ))
    chord_aero, xle_aero = _general_wing_geometry(case, aero_span)
    interpolate(values) = linear_interpolate_1d(span_nodes, values, aero_span)
    u_x = interpolate(wing_kinematics.u_x)
    u_y = interpolate(wing_kinematics.u_y)
    u_z = interpolate(wing_kinematics.u_z)
    theta_x = interpolate(wing_kinematics.theta_x)
    theta_y = interpolate(wing_kinematics.theta_y)
    theta_z = interpolate(wing_kinematics.theta_z)
    span = last(span_nodes) - first(span_nodes)
    grid_wing = generate_panel_grid_and_interpolate(
        span,
        chord_aero,
        xle_aero,
        case.aerodynamic.spanwise_panels,
        case.aerodynamic.chordwise_panels,
        u_x,
        u_y,
        u_z,
        theta_x,
        theta_y,
        theta_z;
        elastic_axis_fraction = case.wing.geometry.elastic_axis_fraction,
    )

    hubs = direct_propeller_hub_positions(
        case.wing,
        case.propellers,
        structural,
        state,
    )
    pitch_axes = Vector{SVector{3,Float64}}(undef, length(case.propellers))
    yaw_axes = similar(pitch_axes)
    for (index, propeller) in enumerate(case.propellers)
        node = propeller.attachment_node
        node_rotation = RotationMatrix(wing_kinematics.theta_z[node], 3) *
            RotationMatrix(wing_kinematics.theta_x[node], 1) *
            RotationMatrix(wing_kinematics.theta_y[node], 2)
        pitch = wing_kinematics.propeller_pitch[index]
        yaw = wing_kinematics.propeller_yaw[index]
        pitch_rotation = RotationMatrix(pitch, 2)
        whirl_rotation = pitch_rotation * RotationMatrix(yaw, 3)
        spin_rotation = RotationMatrix(-structural.angular_speeds[index] * time, 1)
        pitch_axes[index] = node_rotation * SVector(0.0, 1.0, 0.0)
        yaw_axes[index] = -node_rotation * pitch_rotation * SVector(0.0, 0.0, 1.0)
        uvlm.T_pivot_A_current[index] = hubs[index]
        for blade in eachindex(uvlm.grids_prop_ref[index])
            reference_grid = uvlm.grids_prop_ref[index][blade]
            current_grid = uvlm.grids_prop_current[index][blade]
            for chord_index in axes(reference_grid, 2), radial_index in axes(reference_grid, 3)
                point = SVector{3}(reference_grid[:, chord_index, radial_index])
                current_grid[:, chord_index, radial_index] = hubs[index] +
                    node_rotation * whirl_rotation * spin_rotation * point
            end
        end
    end

    _, _, wing_surface = grid_to_surface_panels(
        grid_wing;
        ratios = uvlm.ratio_wing,
        fcore = finite_core,
    )
    uvlm.system.surfaces[1] = wing_surface
    surface_offset = 1
    for index in eachindex(case.propellers)
        for blade in eachindex(uvlm.grids_prop_current[index])
            surface_offset += 1
            uvlm.system.surfaces[surface_offset] = grid_to_surface_panels(
                uvlm.grids_prop_current[index][blade];
                fcore = finite_core,
            )[3]
        end
    end
    return (;
        wing = wing_kinematics,
        aero_span,
        chord_aero,
        xle_aero,
        u_x,
        u_y,
        u_z,
        theta_x,
        theta_y,
        theta_z,
        hubs,
        pitch_axes,
        yaw_axes,
    )
end

function _general_uvlm_structural_load(case, structural, uvlm, kinematics)
    forces = imperial_nodal_forces(uvlm.system)
    positions = imperial_nodal_positions(uvlm.system)
    wing_dof_count = WING_DOFS_PER_NODE * structural.node_count
    global_load = zeros(size(structural.M_global, 1))
    elastic_axis_fraction = case.wing.geometry.elastic_axis_fraction

    # Wing surface-to-beam interpolation remains independent of the direct
    # propeller attachment rule. Each aerodynamic station distributes its
    # work-conjugate wrench to the enclosing beam element.
    for station in eachindex(kinematics.aero_span)
        elastic_axis = SVector(
            kinematics.xle_aero[station] +
                elastic_axis_fraction * kinematics.chord_aero[station] +
                kinematics.u_x[station],
            kinematics.aero_span[station] + kinematics.u_y[station],
            kinematics.u_z[station],
        )
        total_force = SVector(0.0, 0.0, 0.0)
        total_moment = SVector(0.0, 0.0, 0.0)
        for chord_index in axes(forces[1], 1)
            force = forces[1][chord_index, station]
            total_force += force
            total_moment += cross(
                positions[1][chord_index, station] - elastic_axis,
                force,
            )
        end
        generalized_moment = _generalized_wing_moment(
            total_moment,
            kinematics.theta_x[station],
            kinematics.theta_z[station],
        )
        wrench = (
            total_force[2],
            total_force[1],
            -total_force[3],
            generalized_moment[1],
            generalized_moment[2],
            generalized_moment[3],
        )
        for (node, weight) in _beam_interpolation_support(
            case.wing.span_nodes,
            kinematics.aero_span[station],
        )
            add_direct_node_wrench!(global_load, weight .* wrench, node)
        end
    end

    for (index, propeller) in enumerate(case.propellers)
        total_force = SVector(0.0, 0.0, 0.0)
        total_moment = SVector(0.0, 0.0, 0.0)
        hub = kinematics.hubs[index]
        for surface_index in uvlm.prop_surface_indices[index]
            for vertex in eachindex(forces[surface_index])
                force = forces[surface_index][vertex]
                total_force += force
                total_moment += cross(positions[surface_index][vertex] - hub, force)
            end
        end
        # Hub and attachment are the same physical point. Keep the UVLM hub
        # moment itself and introduce no extra hub-to-node cross product.
        wrench_at_node = colocated_hub_wrench(total_force, total_moment, hub, hub)
        node = propeller.attachment_node
        generalized_moment = _generalized_wing_moment(
            wrench_at_node.moment,
            kinematics.wing.theta_x[node],
            kinematics.wing.theta_z[node],
        )
        structural_wrench = (
            total_force[2],
            total_force[1],
            -total_force[3],
            generalized_moment[1],
            generalized_moment[2],
            generalized_moment[3],
        )
        add_direct_node_wrench!(global_load, structural_wrench, node)
        propeller_start = wing_dof_count + 2 * (index - 1)
        global_load[propeller_start + 1] = dot(total_moment, kinematics.pitch_axes[index])
        global_load[propeller_start + 2] = dot(total_moment, kinematics.yaw_axes[index])
    end
    return global_load[structural.free_dofs]
end

"""
Create transactional UVLM callbacks for the common time-domain solver. Every
implicit trial is restored to the same beginning-of-step wake, and the wake is
advanced exactly once after the accepted structural state.
"""
function make_uvlm_aerodynamic_coupling(
    case,
    structural,
    uvlm,
    time_steps;
    interaction::Bool = true,
    finite_core = uvlm.fcore,
    near_field_force_function = near_field_forces!,
)
    length(time_steps) > 0 || throw(ArgumentError("time_steps must not be empty"))
    snapshot = Ref{Any}(nothing)
    last_kinematics = Ref{Any}(nothing)

    begin_step = function (step, state, time)
        copy_surfaces_to_previous!(uvlm.system, uvlm.nsurf)
        snapshot[] = snapshot_uvlm(uvlm.system)
        return nothing
    end
    load = function (state, step, time)
        isnothing(snapshot[]) && throw(ArgumentError(
            "begin_step must be called before evaluating an aerodynamic load",
        ))
        restore_uvlm!(uvlm.system, snapshot[])
        kinematics = _update_general_uvlm_geometry!(
            case,
            structural,
            uvlm,
            state,
            time,
            finite_core,
        )
        last_kinematics[] = kinematics
        propagate_system!(
            uvlm.system,
            uvlm.fs_vec[step],
            time_steps[step];
            additional_velocity = nothing,
            repeated_points = uvlm.repeated_points,
            nwake = uvlm.iwake,
            eta = 0.1,
            calculate_influence_matrix = true,
            near_field_analysis = true,
            derivatives = false,
            near_field_force_function,
            interaction_id = uvlm.surface_interaction_id,
            interaction,
            advance_wake = false,
        )
        return _general_uvlm_structural_load(case, structural, uvlm, kinematics)
    end
    commit_step = function (step, state, time, time_step)
        load(state, step, time)
        advance_wake!(
            uvlm.system,
            uvlm.fs_vec[step],
            time_step;
            additional_velocity = nothing,
            repeated_points = uvlm.repeated_points,
            nwake = uvlm.iwake,
            interaction_id = uvlm.surface_interaction_id,
            interaction,
        )
        commit_wake_rows!(uvlm.iwake, uvlm.nwake)
        return nothing
    end
    return (; load, begin_step, commit_step, last_kinematics)
end
