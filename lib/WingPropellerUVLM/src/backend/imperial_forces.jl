"""
    imperial_panel_area(panel)

Return the quadrilateral area using the two-diagonal average used by Imperial
College London's C++ UVLM implementation.
"""
function imperial_panel_area(panel::SurfacePanel)
    vertices = (top_left(panel), bottom_left(panel), bottom_right(panel), top_right(panel))
    sides = ntuple(i -> norm(vertices[mod1(i + 1, 4)] - vertices[i]), 4)

    triangle_area(a, b, c) = begin
        s = (a + b + c) / 2
        sqrt(max(zero(s), s * (s - a) * (s - b) * (s - c)))
    end

    diagonal_02 = norm(vertices[3] - vertices[1])
    diagonal_13 = norm(vertices[4] - vertices[2])
    area_02 = triangle_area(sides[1], sides[2], diagonal_02) +
              triangle_area(sides[3], sides[4], diagonal_02)
    area_13 = triangle_area(sides[2], sides[3], diagonal_13) +
              triangle_area(sides[1], sides[4], diagonal_13)
    return (area_02 + area_13) / 2
end

"""
    imperial_segment_force(velocity, r1, r2, delta_gamma, rho)

Evaluate the dimensional Joukovski force on one vortex segment. Segment
orientation is from `r1` to `r2`, matching `f = rho*delta_gamma*v×(r2-r1)` in
Imperial's `calculate_static_forces_unsteady` routine.
"""
imperial_segment_force(velocity, r1, r2, delta_gamma, rho) =
    rho * delta_gamma * cross(velocity, r2 - r1)

"""
    imperial_unsteady_panel_force(panel, gamma_dot, rho)

Evaluate the dimensional non-circulatory panel force
`-rho*area*normal*gamma_dot` used by Imperial's `calculate_dynamic_forces`.
"""
imperial_unsteady_panel_force(panel, gamma_dot, rho) =
    -rho * imperial_panel_area(panel) * normal(panel) * gamma_dot

# VortexLattice.jl traverses a positive-circulation ring in the opposite
# direction to Imperial's zeta vertex ordering. Convert at this boundary and
# keep all force equations below in the Imperial convention.
_to_imperial_circulation(gamma) = -gamma

_imperial_spanwise_circulation_jump(gamma, i, j) =
    i == 1 ? -gamma[i, j] : gamma[i - 1, j] - gamma[i, j]

function _imperial_chordwise_circulation_jump(gamma, i, j)
    ns = size(gamma, 2)
    if j == 1
        return gamma[i, 1]
    elseif j == ns + 1
        return -gamma[i, ns]
    end
    return gamma[i, j] - gamma[i, j - 1]
end

function _imperial_base_velocity(rc, ref, fs, additional_velocity, segment_motion)
    velocity = freestream_velocity(fs) + rotational_velocity(rc, fs, ref)
    if !isnothing(additional_velocity)
        velocity += additional_velocity(rc)
    end
    if !isnothing(segment_motion)
        velocity += segment_motion
    end
    return velocity
end

function _imperial_induced_velocity(rc, receiving_index, segment_kind,
    surfaces, wakes, Γ, isurf; symmetric, nwake, surface_id,
    wake_finite_core, wake_shedding_locations, trailing_vortices, xhat,
    interaction_id, interaction)

    induced = zero(rc)
    jΓ = 0
    for jsurf in eachindex(surfaces)
        sending = surfaces[jsurf]
        number_of_panels = length(sending)
        same_group = interaction || (interaction_id[isurf] == interaction_id[jsurf])
        if !same_group
            jΓ += number_of_panels
            continue
        end

        wake_panels = nwake[jsurf] > 0
        shedding_locations = isnothing(wake_shedding_locations) ?
            nothing : wake_shedding_locations[jsurf]
        circulation = view(Γ, jΓ + 1:jΓ + number_of_panels)

        skip_top = ()
        skip_bottom = ()
        skip_left = ()
        skip_right = ()
        if isurf == jsurf
            i, j = Tuple(receiving_index)
            if segment_kind === :span
                skip_top = (CartesianIndex(i, j),)
                skip_bottom = (CartesianIndex(i - 1, j),)
            elseif segment_kind === :chord
                if j <= size(sending, 2)
                    skip_left = (CartesianIndex(i, j),)
                end
                if j > 1
                    skip_right = (CartesianIndex(i, j - 1),)
                end
            end
        end

        induced += induced_velocity(rc, sending, circulation;
            finite_core = surface_id[isurf] != surface_id[jsurf],
            wake_shedding_locations = shedding_locations,
            symmetric = symmetric[jsurf],
            trailing_vortices = trailing_vortices[jsurf] && !wake_panels,
            xhat, skip_top, skip_bottom, skip_left, skip_right)

        if wake_panels
            induced += induced_velocity(rc, wakes[jsurf];
                finite_core = wake_finite_core[jsurf] ||
                    (surface_id[isurf] != surface_id[jsurf]),
                symmetric = symmetric[jsurf], nc = nwake[jsurf],
                trailing_vortices = trailing_vortices[jsurf], xhat)
        end
        jΓ += number_of_panels
    end
    return induced
end

function _imperial_panel_properties!(props, isurf, panels, Γjulia,
    span_forces, chord_forces, unsteady_forces, span_velocities, ref)

    nc, ns = size(panels)
    dynamic_pressure_area = (ref.rho * ref.V^2 / 2) * ref.S
    zero_force = zero(eltype(span_forces))

    for j in 1:ns, i in 1:nc
        span_force = span_forces[i, j]

        left_weight = j == 1 ? one(eltype(Γjulia)) : one(eltype(Γjulia)) / 2
        right_weight = j == ns ? one(eltype(Γjulia)) : one(eltype(Γjulia)) / 2
        left_force = left_weight * chord_forces[i, j]
        right_force = right_weight * chord_forces[i, j + 1]

        # Imperial transfers the unsteady force to all four corners except at
        # the trailing edge, where the downstream pair is deliberately skipped.
        if i == nc
            span_force += unsteady_forces[i, j] / 2
        else
            left_force += unsteady_forces[i, j] / 2
            right_force += unsteady_forces[i, j] / 2
        end

        props[isurf][i, j] = PanelProperties(
            Γjulia[i, j] / ref.V,
            span_velocities[i, j] / ref.V,
            span_force / dynamic_pressure_area,
            left_force / dynamic_pressure_area,
            right_force / dynamic_pressure_area,
            zero_force,
        )
    end
    return props
end

"""
    near_field_forces!(properties, surfaces, wakes, reference, freestream, Γ; ...)

Calculate dimensional segment and panel forces following Imperial College
London's C++ UVLM force implementation. The spanwise and chordwise circulation
jumps, Joukovski cross-product orientation, unsteady pressure-rate force, and
trailing-edge treatment mirror `include/postproc.h`.
"""
function near_field_forces!(props, surfaces, wakes, ref, fs, Γ;
    dΓdt, additional_velocity, Vh, Vv, symmetric, nwake, surface_id,
    wake_finite_core, wake_shedding_locations, trailing_vortices, xhat,
    interaction_id = surface_id, interaction::Bool = true)

    nsurf = length(surfaces)
    TF = eltype(Γ)
    chord_seg_forces = Vector{Matrix{SVector{3, TF}}}(undef, nsurf)
    span_seg_forces = Vector{Matrix{SVector{3, TF}}}(undef, nsurf)
    unsteady_forces = Vector{Matrix{SVector{3, TF}}}(undef, nsurf)

    iΓ = 0
    for isurf in eachindex(surfaces)
        panels = surfaces[isurf]
        nc, ns = size(panels)
        linear_indices = LinearIndices(panels)
        Γjulia = reshape(view(Γ, iΓ + 1:iΓ + length(panels)), nc, ns)
        Γimperial = _to_imperial_circulation.(Γjulia)

        chord_forces = fill(zero(SVector{3, TF}), nc, ns + 1)
        span_forces = fill(zero(SVector{3, TF}), nc + 1, ns)
        panel_unsteady = fill(zero(SVector{3, TF}), nc, ns)
        span_velocities = fill(zero(SVector{3, TF}), nc, ns)

        # Spanwise segments. As in Imperial's loop, the final chordwise row of
        # spanwise segments is allocated but remains zero.
        for j in 1:ns, i in 1:nc
            panel = panels[i, j]
            r1 = top_left(panel)
            r2 = top_right(panel)
            rc = (r1 + r2) / 2
            motion = isnothing(Vh) ? nothing : Vh[isurf][i, j]
            velocity = _imperial_base_velocity(rc, ref, fs, additional_velocity, motion)
            velocity += _imperial_induced_velocity(rc, CartesianIndex(i, j), :span,
                surfaces, wakes, Γ, isurf; symmetric, nwake, surface_id,
                wake_finite_core, wake_shedding_locations, trailing_vortices, xhat,
                interaction_id, interaction)

            delta_gamma = _imperial_spanwise_circulation_jump(Γimperial, i, j)
            span_forces[i, j] = imperial_segment_force(
                velocity, r1, r2, delta_gamma, ref.rho)
            span_velocities[i, j] = velocity
        end

        # Chordwise segments, including both outer spanwise boundaries.
        for j in 1:ns + 1, i in 1:nc
            if j <= ns
                panel = panels[i, j]
                r1 = top_left(panel)
                r2 = bottom_left(panel)
            else
                panel = panels[i, ns]
                r1 = top_right(panel)
                r2 = bottom_right(panel)
            end
            rc = (r1 + r2) / 2
            motion = isnothing(Vv) ? nothing : Vv[isurf][i, j]
            velocity = _imperial_base_velocity(rc, ref, fs, additional_velocity, motion)
            velocity += _imperial_induced_velocity(rc, CartesianIndex(i, j), :chord,
                surfaces, wakes, Γ, isurf; symmetric, nwake, surface_id,
                wake_finite_core, wake_shedding_locations, trailing_vortices, xhat,
                interaction_id, interaction)

            delta_gamma = _imperial_chordwise_circulation_jump(Γimperial, i, j)
            chord_forces[i, j] = imperial_segment_force(
                velocity, r1, r2, delta_gamma, ref.rho)
        end

        if !isnothing(dΓdt)
            for j in 1:ns, i in 1:nc
                gamma_dot = _to_imperial_circulation(
                    dΓdt[iΓ + linear_indices[i, j]])
                panel_unsteady[i, j] = imperial_unsteady_panel_force(
                    panels[i, j], gamma_dot, ref.rho)
            end
        end

        _imperial_panel_properties!(props, isurf, panels, Γjulia,
            span_forces, chord_forces, panel_unsteady, span_velocities, ref)
        chord_seg_forces[isurf] = chord_forces
        span_seg_forces[isurf] = span_forces
        unsteady_forces[isurf] = panel_unsteady
        iΓ += length(panels)
    end

    return props, chord_seg_forces, span_seg_forces, unsteady_forces
end

"""
    imperial_nodal_forces(span_forces, chord_forces, unsteady_forces)

Transfer dimensional segment forces to vortex vertices using Imperial's
half-to-each-endpoint rule. Unsteady panel forces contribute one quarter to
each corner, except that the downstream pair of trailing-edge vertices is
skipped exactly as in the C++ implementation.
"""
function imperial_nodal_forces(span_forces::AbstractMatrix,
    chord_forces::AbstractMatrix, unsteady_forces::AbstractMatrix)

    nc, ns = size(unsteady_forces)
    @assert size(span_forces) == (nc + 1, ns)
    @assert size(chord_forces) == (nc, ns + 1)
    nodal = fill(zero(eltype(span_forces)), nc + 1, ns + 1)

    for j in 1:ns + 1, i in 1:nc + 1
        if j > 1
            nodal[i, j] += span_forces[i, j - 1] / 2
        end
        if j <= ns
            nodal[i, j] += span_forces[i, j] / 2
        end
        if i > 1
            nodal[i, j] += chord_forces[i - 1, j] / 2
        end
        if i <= nc
            nodal[i, j] += chord_forces[i, j] / 2
        end
    end

    for j in 1:ns, i in 1:nc
        nodal[i, j] += unsteady_forces[i, j] / 4
        nodal[i, j + 1] += unsteady_forces[i, j] / 4
        if i < nc
            nodal[i + 1, j] += unsteady_forces[i, j] / 4
            nodal[i + 1, j + 1] += unsteady_forces[i, j] / 4
        end
    end
    return nodal
end

"""
    imperial_nodal_forces(system)

Return one `(nc+1, ns+1)` matrix of dimensional body-frame nodal forces per
surface. A near-field analysis must have been performed first.
"""
function imperial_nodal_forces(system::System)
    @assert system.near_field_analysis[] "Near field analysis required"
    return [imperial_nodal_forces(system.span_seg_forces[i],
        system.chord_seg_forces[i], system.unsteady_forces[i])
        for i in eachindex(system.surfaces)]
end
