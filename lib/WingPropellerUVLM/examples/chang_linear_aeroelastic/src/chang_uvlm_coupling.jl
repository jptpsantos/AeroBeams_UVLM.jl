# Chang linear-structure to UVLM adapter.
#
# These functions intentionally remain with the example because they encode
# its structural DOF order and the Z-down structural/Z-up aerodynamic sign
# conversion. Generic grid deformation, UVLM stepping, segment-to-nodal loads,
# and generalized-alpha integration live in WingPropellerUVLM itself.

"""Initialize the case UVLM system from `chang_aerodynamic_options`."""
function initialize_chang_uvlm(parameters,options)
    return initialize_chang_uvlm(parameters;
        finite_core=options.finite_core,
        elastic_axis_fraction=options.elastic_axis_fraction,
        propeller_pivot_offset_A=options.propeller_pivot_offset_A,
        physical_hub_center_A=options.physical_hub_center_A,
        maximum_wake_rows_wing=options.maximum_wake_rows_wing,
        maximum_wake_rows_propeller=options.maximum_wake_rows_propeller
    )
end

"""Collect the wake controls while sharing the initialized mutable storage."""
function chang_wake_context(uvlm; interaction_on::Bool)
    return (;
        surface_count=uvlm.nsurf,
        repeated_points=uvlm.repeated_points,
        maximum_rows=uvlm.nwake,
        active_rows=uvlm.iwake,
        interaction_ids=uvlm.surface_interaction_id,
        interaction_on,
        saved_steps=uvlm.save,
        surface_history=uvlm.surface_history
    )
end

"""
    initialize_chang_uvlm(parameters; kwargs...)

Create the coupled wing--propeller UVLM state from the geometry and operating
point in `parameters`. The options supply
the numerical regularization, reference locations, and retained-wake limits.
"""
function initialize_chang_uvlm(parameters;
    finite_core,
    elastic_axis_fraction,
    propeller_pivot_offset_A,
    physical_hub_center_A,
    maximum_wake_rows_wing,
    maximum_wake_rows_propeller
)
    return initialize_wing_propeller_uvlm_system(
        xle=parameters.xle,
        yle=parameters.yle,
        zle=parameters.zle,
        chord_geo=parameters.chord_geo,
        theta_geo=parameters.theta_geo,
        phi_geo=parameters.phi_geo,
        ns_wing=parameters.ns_wing,
        nc_wing=parameters.nc_wing,
        mirror_wing=parameters.mirror_wing,
        spacing_s_wing=parameters.spacing_s_wing,
        spacing_c_wing=parameters.spacing_c_wing,
        R_prop=parameters.R_prop,
        c_prop=parameters.c_prop,
        ns_prop=parameters.ns_prop,
        nc_prop=parameters.nc_prop,
        blade_twists_prop=parameters.blade_twists_prop,
        Nb_prop=parameters.Nb_prop,
        Npropellers=parameters.Npropellers,
        span_nodes=parameters.span_nodes,
        prop_attach_nodes=parameters.prop_attach_nodes,
        chord=parameters.chord,
        xle_distribution=parameters.xle_distribution,
        ref=parameters.ref,
        symmetric_wing=parameters.symmetric_wing,
        fs=parameters.fs,
        dt=parameters.dt,
        nnodes=parameters.nnodes,
        prop_pivot_offset_from_ea_A=propeller_pivot_offset_A,
        hub_center_prop_A=physical_hub_center_A,
        fcore=finite_core,
        elastic_axis_fraction=elastic_axis_fraction,
        maximum_wake_rows_wing=maximum_wake_rows_wing,
        maximum_wake_rows_propeller=maximum_wake_rows_propeller,
        verbose=true
    )
end

"""
    update_aero_geometry_for_state!(model, workspace, q_free, time_np1)

Map one structural displacement guess to the UVLM geometry at `time_np1`.

The Chang structural basis is `[span, chord, down]`; the UVLM aerodynamic
basis is `[chord, span, up]`. The function reconstructs the clamped root,
applies all six wing DOFs, adds each propeller's pitch/yaw modal motion and
prescribed spin, rebuilds the panel surfaces, and returns the aerodynamic-frame
wing kinematics needed to form moment arms during load transfer.

This function changes geometry only. It does not solve circulation or advance
the wake.
"""

function update_aero_geometry_for_state!(model,workspace,q_free::AbstractVector,time_np1::Real)
    (; ndof,span_length,chord,xle_distribution,ns_wing,nc_wing,
        Npropellers,prop_attach_nodes,Ω,Nb_prop) = model.parameters
    ndof_wing_free = model.structural.ndof_wing_free
    (; system,ea_x_aero,attach_node_y,T_pivot_A_current,pitch_axis_A_current,
        yaw_axis_A_current,T_hub_A_current,T_load_A_current,grids_prop_ref,
        grids_prop_current,ratio_wing) = workspace
    (; elastic_axis_fraction) = model.aerodynamic_options
    FCORE = model.aerodynamic_options.finite_core
    hub_center_prop_A = model.aerodynamic_options.physical_hub_center_A
    hub_center_load_A = model.aerodynamic_options.load_center_A

    # Split the global free-state vector according to the assembly convention:
    # all wing free DOFs first, then [pitch, yaw] for each propeller.
    qWingFree = q_free[1:ndof_wing_free]
    qPropellerFree = q_free[(ndof_wing_free+1):end]

    # Reinsert the clamped root DOFs and extract each component along the beam.
    displacementSpanStructural = vcat(0.0,qWingFree[1:ndof:end])
    displacementChordStructural = vcat(0.0,qWingFree[2:ndof:end])
    displacementVerticalStructural = vcat(0.0,qWingFree[3:ndof:end])
    rotationSpanStructural = vcat(0.0,qWingFree[4:ndof:end])
    rotationChordStructural = vcat(0.0,qWingFree[5:ndof:end])
    rotationVerticalStructural = vcat(0.0,qWingFree[6:ndof:end])

    # Proper structural-to-aerodynamic basis mapping:
    # (span, chord, down) -> (y, x, -z), with the same mapping for rotations.
    uXAero = displacementChordStructural
    uYAero = displacementSpanStructural
    uZAero = -displacementVerticalStructural
    thetaXAero = rotationChordStructural
    thetaYAero = rotationSpanStructural
    thetaZAero = -rotationVerticalStructural

    # Interpolate nodal translations/rotations spanwise and move every wing
    # lattice point about the elastic axis.
    gridWing = generate_panel_grid_and_interpolate(
        span_length,
        chord,
        xle_distribution,
        ns_wing,
        nc_wing,
        uXAero,
        uYAero,
        uZAero,
        thetaXAero,
        thetaYAero,
        thetaZAero;
        elastic_axis_fraction=elastic_axis_fraction
    )

    # Propeller motion is composed in this order: prescribed rotor spin,
    # pylon pitch/yaw, then the selected attachment-node rotation/translation.
    for propeller_index = 1:Npropellers
        attachNode = prop_attach_nodes[propeller_index]
        pitchAero = qPropellerFree[2*(propeller_index-1)+1]
        yawAero = -qPropellerFree[2*(propeller_index-1)+2]

        # Direct attachment: the aerodynamic hub is colocated with the selected
        # structural node. No interpolation weights are used to split the propeller
        # motion between adjacent beam nodes.
        uXAttachment = uXAero[attachNode]
        uYAttachment = uYAero[attachNode]
        uZAttachment = uZAero[attachNode]
        thetaXAttachment = thetaXAero[attachNode]
        thetaYAttachment = thetaYAero[attachNode]
        thetaZAttachment = thetaZAero[attachNode]

        pivot = SVector(
            ea_x_aero[propeller_index] + uXAttachment,
            attach_node_y[propeller_index] + uYAttachment,
            uZAttachment
        )
        T_pivot_A_current[propeller_index] = pivot

        wingRotation =
            RotationMatrix(thetaZAttachment,3) *
            RotationMatrix(thetaXAttachment,1) *
            RotationMatrix(thetaYAttachment,2)
        pitchRotation = RotationMatrix(pitchAero,2)
        yawRotation = RotationMatrix(yawAero,3)
        whirlRotation = pitchRotation * yawRotation
        spinRotation = RotationMatrix(-Ω * time_np1,1)

        # Work-conjugate axes for R_whirl = R_y(pitch) R_z(-yaw). The pitch
        # axis is the wing-rotated y axis. The structural yaw coordinate is
        # opposite aerodynamic z and its instantaneous axis is also carried
        # by the preceding pitch rotation.
        pitch_axis_A_current[propeller_index] = wingRotation *
                                                SVector(0.0,1.0,0.0)
        yaw_axis_A_current[propeller_index] = -wingRotation * pitchRotation *
                                              SVector(0.0,0.0,1.0)

        # The pylon pivot is the selected elastic-axis node. Keep the physical
        # hub and assumed-mode load point separate so the rotor disk retains its
        # full static pylon offset and its pitch/yaw motion is work-conjugate to
        # the structural modal coordinates.
        T_hub_A_current[propeller_index] =
            pivot + wingRotation *
                    (hub_center_prop_A + (whirlRotation * hub_center_load_A - hub_center_load_A))
        T_load_A_current[propeller_index] = pivot + wingRotation *
                                                    (whirlRotation * hub_center_load_A)

        for bladeIndex = 1:Nb_prop
            referenceGrid = grids_prop_ref[propeller_index][bladeIndex]
            currentGrid = grids_prop_current[propeller_index][bladeIndex]
            for chordIndex in axes(referenceGrid,2),radialIndex in axes(referenceGrid,3)
                referencePoint = SVector{3}(referenceGrid[:,chordIndex,radialIndex])
                spunPoint = spinRotation * referencePoint
                staticPoint = spunPoint + hub_center_prop_A
                modalPoint = spunPoint + hub_center_load_A
                whirledPoint = staticPoint + (whirlRotation * modalPoint - modalPoint)
                currentGrid[:,chordIndex,radialIndex] =
                    pivot + wingRotation * whirledPoint
            end
        end
    end

    # Replace the panel geometry in the existing System. Keeping the System
    # object itself preserves all preallocated circulation, force, and wake
    # storage used by `propagate_system!`.
    _,_,wingSurface = grid_to_surface_panels(
        gridWing;
        ratios=ratio_wing,
        fcore=FCORE
    )
    propeller_surfaces = Vector{typeof(wingSurface)}()
    for propeller_index = 1:Npropellers,bladeIndex = 1:Nb_prop
        push!(
            propeller_surfaces,
            grid_to_surface_panels(
                grids_prop_current[propeller_index][bladeIndex];
                fcore=FCORE
            )[3]
        )
    end

    system.surfaces[1] = wingSurface
    for surface_index in eachindex(propeller_surfaces)
        system.surfaces[surface_index+1] = propeller_surfaces[surface_index]
    end

    return (;
        u_x_A=uXAero,
        u_y_A=uYAero,
        u_z_A=uZAero,
        theta_x_A=thetaXAero,
        theta_y_A=thetaYAero,
        theta_z_A=thetaZAero,
        propeller_pitch_axes_A=pitch_axis_A_current,
        propeller_yaw_axes_A=yaw_axis_A_current
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
function chang_wing_generalized_moment(moment_A,theta_x_A,theta_z_A)
    sx,cx = sincos(theta_x_A)
    sz,cz = sincos(theta_z_A)
    return SVector(
        -sz * cx * moment_A[1] + cz * cx * moment_A[2] + sx * moment_A[3],
        cz * moment_A[1] + sz * moment_A[2],
        -moment_A[3]
    )
end

"""
    assemble_structural_aero_load!(model, workspace, kinematics; step=0, print_loads=false)

Transfer the dimensional UVLM vertex forces to the Chang free-DOF
ordering. Wing forces are summed chordwise and moments are formed about the
deformed elastic axis and projected onto the instantaneous wing rotation axes.
Blade forces are reduced to propeller hub wrenches and transferred directly to
the same structural attachment node. Modal moments use the selected pitch/yaw
projection.

Returns one generalized-load vector ordered exactly like the reduced
structural state used by `M`, `C`, and `K`.
"""
function assemble_structural_aero_load!(model,workspace,kinematics;
    step::Int=0,print_loads::Bool=false)

    (; nnodes,ndof,NDOF,ndof_P,chord,xle_distribution,span_nodes,
        nc_wing,ns_wing,Npropellers,Nb_prop,prop_attach_nodes) = model.parameters
    (; free_dofs,ndof_wing_free) = model.structural
    (; system,TF,nodal_forces_wing,nodal_moments_wing,EA_nodes_wing,
        nodal_forces_prop,prop_surface_indices,T_pivot_A_current,
        T_hub_A_current,T_load_A_current) = workspace
    elastic_axis_fraction = model.aerodynamic_options.elastic_axis_fraction
    projection = model.config.simulation.propeller_moment_projection

    zeroVector = SVector{3,TF}(0.0,0.0,0.0)
    fill!(nodal_forces_wing,zeroVector)
    fill!(nodal_moments_wing,zeroVector)

    # These arrays are dimensional forces and their matching vortex-lattice
    # vertex positions. Using the paired positions preserves moment arms.
    surfaceForcesAero = imperial_nodal_forces(system)
    surfacePositionsAero = imperial_nodal_positions(system)
    nodal_forces_wing .= surfaceForcesAero[1]

    # Build the deformed elastic-axis line used as the wing moment reference.
    for nodeIndex = 1:nnodes
        elastic_axis_x =
            xle_distribution[nodeIndex] +
            chord[nodeIndex] * elastic_axis_fraction + kinematics.u_x_A[nodeIndex]
        EA_nodes_wing[nodeIndex] = SVector(
            elastic_axis_x,
            span_nodes[nodeIndex] + kinematics.u_y_A[nodeIndex],
            kinematics.u_z_A[nodeIndex]
        )
    end
    for chordIndex = 1:(nc_wing+1),spanIndex = 1:(ns_wing+1)
        nodal_moments_wing[chordIndex,spanIndex] = cross(
            surfacePositionsAero[1][chordIndex,spanIndex] - EA_nodes_wing[spanIndex],
            nodal_forces_wing[chordIndex,spanIndex]
        )
    end

    wing_loads = zeros(NDOF)
    force_x = sum(getindex.(nodal_forces_wing,1),dims=1)
    force_y = sum(getindex.(nodal_forces_wing,2),dims=1)
    force_z = sum(getindex.(nodal_forces_wing,3),dims=1)
    moment_x = sum(getindex.(nodal_moments_wing,1),dims=1)
    moment_y = sum(getindex.(nodal_moments_wing,2),dims=1)
    moment_z = sum(getindex.(nodal_moments_wing,3),dims=1)

    for nodeIndex = 1:nnodes
        offset = ndof * (nodeIndex - 1)
        generalizedMoment = chang_wing_generalized_moment(
            SVector(moment_x[nodeIndex],moment_y[nodeIndex],moment_z[nodeIndex]),
            kinematics.theta_x_A[nodeIndex],
            kinematics.theta_z_A[nodeIndex]
        )
        wing_loads[offset+1] = force_y[nodeIndex]
        wing_loads[offset+2] = force_x[nodeIndex]
        wing_loads[offset+3] = -force_z[nodeIndex]
        wing_loads[offset+4] = generalizedMoment[1]
        wing_loads[offset+5] = generalizedMoment[2]
        wing_loads[offset+6] = generalizedMoment[3]
    end

    propeller_loads = zeros(ndof_P)
    # Reduce all blade vertex forces of each propeller to one resultant force
    # and moment, then apply them to both the pylon modal and wing attachment
    # coordinates without discarding the appropriate lever arms.
    for propeller_index = 1:Npropellers
        attachNode = prop_attach_nodes[propeller_index]
        totalForce = zeroVector
        totalMomentAboutHub = zeroVector
        hub = T_hub_A_current[propeller_index]

        for bladeIndex = 1:Nb_prop
            surface_index = prop_surface_indices[propeller_index][bladeIndex]
            bladeForces = nodal_forces_prop[propeller_index][bladeIndex]
            bladeForces .= surfaceForcesAero[surface_index]
            bladePositions = surfacePositionsAero[surface_index]

            for chordIndex in axes(bladeForces,1),radialIndex in axes(bladeForces,2)
                force = SVector{3,TF}(bladeForces[chordIndex,radialIndex])
                position = bladePositions[chordIndex,radialIndex]
                totalForce += force
                totalMomentAboutHub += cross(position - hub,force)
            end
        end

        pivot = T_pivot_A_current[propeller_index]
        modalLoadPoint = T_load_A_current[propeller_index]
        # Pitch/yaw generalized loads use the assumed-mode load arm. The wing
        # node receives the same physical wrench translated from hub to pivot.
        modalMoment = totalMomentAboutHub +
                      cross(modalLoadPoint - pivot,totalForce)
        wing_moment = totalMomentAboutHub + cross(hub - pivot,totalForce)

        if projection == :exact_virtual_work
            pitchAxis = kinematics.propeller_pitch_axes_A[propeller_index]
            yawAxis = kinematics.propeller_yaw_axes_A[propeller_index]
            propeller_loads[2*(propeller_index-1)+1] =
                dot(modalMoment,pitchAxis)
            propeller_loads[2*(propeller_index-1)+2] =
                dot(modalMoment,yawAxis)
        else
            # Original small-angle projection onto fixed aerodynamic axes.
            propeller_loads[2*(propeller_index-1)+1] = modalMoment[2]
            propeller_loads[2*(propeller_index-1)+2] = -modalMoment[3]
        end

        # The propeller aerodynamic load is transferred directly to the single
        # structural attachment node without introducing an artificial
        # neighboring-node lever arm.
        thetaXAttachment = kinematics.theta_x_A[attachNode]
        thetaZAttachment = kinematics.theta_z_A[attachNode]
        generalizedWingMoment = chang_wing_generalized_moment(
            wing_moment,
            thetaXAttachment,
            thetaZAttachment
        )
        structuralWrench = (
            totalForce[2],
            totalForce[1],
            -totalForce[3],
            generalizedWingMoment[1],
            generalizedWingMoment[2],
            generalizedWingMoment[3]
        )
        add_direct_node_wrench!(wing_loads,structuralWrench,attachNode;
            dofs_per_node=ndof)

        if print_loads
            println("\n--- Step $step Propeller P$(propeller_index) aerodynamic loads ---")
            println("  total_force_A            = $totalForce")
            println("  total_moment_about_hub_A = $totalMomentAboutHub")
            println("  total_moment_modal_A     = $modalMoment")
            println("  total_moment_wing_A      = $wing_moment")
            println(
                "  propeller moment projection = " *
                "$projection"
            )
        end
    end

    freeWingLoads = wing_loads[free_dofs[1:ndof_wing_free]]
    return vcat(freeWingLoads,propeller_loads)
end

"""
    aero_load_for_state!(model, workspace, snapshot, state, step; print_loads=false)

Evaluate the complete aerodynamic generalized-load operator for one structural
fixed-point guess at `t[step+1]`.

Every call restores the beginning-of-step snapshot, updates geometry, solves
one UVLM circulation/load trial, and transfers its loads. Wake convection and
shedding are deliberately disabled here. Therefore repeated calls during the
same coupled step are deterministic trials from `t[step]`, not successive wake
steps. On return, `system` contains the circulation/load trial corresponding to
`state`; the run driver advances its wake only after convergence.
"""
function aero_load_for_state!(model,workspace,snapshot,state::AbstractVector,step::Int;
    print_loads::Bool=false)

    (; t,dt) = model.parameters
    (; system,fs_vec,repeated_points,iwake,surface_interaction_id) = workspace

    # Transaction rollback: discard the previous coupling iterate's wake,
    # circulation, geometry, and surface-history changes.
    restore_uvlm!(system,snapshot)
    # Geometry is evaluated at the end of the physical step, including the
    # prescribed rotor azimuth at that time.
    kinematics = update_aero_geometry_for_state!(model,workspace,state,t[step+1])
    # Solve the unsteady aerodynamic trial through circulation, gamma-dot, and
    # selected near-field segment loads. The expensive free-wake convection is a
    # physical-time update, so it is deferred until this structural state has
    # converged and is accepted by the run driver.
    propagate_system!(
        system,
        fs_vec[step],
        dt[step];
        additional_velocity=nothing,
        repeated_points=repeated_points,
        nwake=iwake,
        eta=0.1,
        calculate_influence_matrix=true,
        near_field_analysis=true,
        derivatives=false,
        near_field_force_function=model.near_field_force_function,
        interaction_id=surface_interaction_id,
        interaction=model.config.simulation.interaction_on,
        advance_wake=false
    )
    # Convert the aerodynamic trial into a load vector for the structural
    # corrector. The caller decides whether to subtract the trim baseline.
    return assemble_structural_aero_load!(
        model,workspace,
        kinematics;
        step=step,
        print_loads=print_loads
    )
end
