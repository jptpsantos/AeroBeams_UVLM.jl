# Corrected legacy segment-force formulation for the wing-propeller
# aeroelastic model. Only the block derived from Imperial College London's C++
# UVLM implementation is kept in this file.

function _legacy_imperial_base_velocity(
    location,
    reference,
    freestream,
    additional_velocity,
    surface_motion,
)
    velocity = freestream_velocity(freestream)
    velocity += rotational_velocity(location, freestream, reference)
    if !isnothing(additional_velocity)
        velocity += additional_velocity(location)
    end
    if !isnothing(surface_motion)
        velocity += surface_motion
    end
    return velocity
end

"""
    _legacy_imperial_induced_velocity(location, panel_index, segment_kind,
        surfaces, wakes, circulation, receiving_surface; ...)

Return the velocity induced at a legacy force segment. Induction from the
receiving bound surface is retained, while only the coincident filament(s) are
excluded to avoid evaluating a vortex segment on itself. Wake induction is
retained for both segment types.
"""
function _legacy_imperial_induced_velocity(
    location,
    panel_index,
    segment_kind::Symbol,
    surfaces,
    wakes,
    surface_circulations,
    receiving_surface;
    symmetric,
    nwake,
    wake_finite_core,
    wake_is_active,
    shedding_locations,
    trailing_vortices,
    xhat,
    same_interaction_group,
    different_surface_id,
)
    segment_kind in (:spanwise, :chordwise) ||
        throw(ArgumentError("segment_kind must be :spanwise or :chordwise"))

    velocity = zero(location)
    for sending_surface in eachindex(surfaces)
        if !same_interaction_group[receiving_surface, sending_surface]
            continue
        end

        sending = surfaces[sending_surface]
        sending_circulation = surface_circulations[sending_surface]
        ids_are_different =
            different_surface_id[receiving_surface, sending_surface]

        if receiving_surface == sending_surface
            if segment_kind === :spanwise
                velocity += induced_velocity(
                    panel_index,
                    sending,
                    sending_circulation;
                    finite_core = ids_are_different,
                    wake_shedding_locations = shedding_locations[sending_surface],
                    symmetric = symmetric[sending_surface],
                    trailing_vortices = trailing_vortices[sending_surface] &&
                        !wake_is_active[sending_surface],
                    xhat = xhat,
                )
            else
                chordwise_index, spanwise_segment_index = Tuple(panel_index)
                skip_left = spanwise_segment_index <= size(sending, 2) ?
                    (CartesianIndex(chordwise_index, spanwise_segment_index),) : ()
                skip_right = spanwise_segment_index > 1 ?
                    (CartesianIndex(chordwise_index, spanwise_segment_index - 1),) : ()
                velocity += induced_velocity(
                    location,
                    sending,
                    sending_circulation;
                    finite_core = ids_are_different,
                    wake_shedding_locations = shedding_locations[sending_surface],
                    symmetric = symmetric[sending_surface],
                    trailing_vortices = trailing_vortices[sending_surface] &&
                        !wake_is_active[sending_surface],
                    xhat = xhat,
                    skip_left,
                    skip_right,
                )
            end
        else
            velocity += induced_velocity(location, sending, sending_circulation;
                finite_core = ids_are_different,
                wake_shedding_locations = shedding_locations[sending_surface],
                symmetric = symmetric[sending_surface],
                trailing_vortices = trailing_vortices[sending_surface] &&
                    !wake_is_active[sending_surface],
                xhat = xhat)
        end

        if wake_is_active[sending_surface]
            velocity += induced_velocity(location, wakes[sending_surface];
                finite_core = wake_finite_core[sending_surface] ||
                    ids_are_different,
                symmetric = symmetric[sending_surface],
                nc = nwake[sending_surface],
                trailing_vortices = trailing_vortices[sending_surface],
                xhat = xhat)
        end
    end

    return velocity
end

"""
    legacy_imperial_segment_forces!(properties, surfaces, wakes, reference,
        freestream, circulation; ...)

Calculate the dimensional chordwise, spanwise, and unsteady panel loads using
only the Imperial College-inspired segment-force block from the original
wing-propeller UVLM implementation.

This compatibility formulation intentionally retains its original choices:

- the trailing-edge spanwise segment force is zero;
- only coincident bound filaments are excluded from segment velocities;
- unsteady force uses the product of the spanwise and chordwise segment lengths.

`properties` is returned unchanged. Aeroelastic callers must consume the three
dimensional segment-force arrays returned after it.
"""
function legacy_imperial_segment_forces!(
    properties,
    surfaces,
    wakes,
    reference,
    freestream,
    circulation;
    dΓdt,
    additional_velocity,
    Vh,
    Vv,
    symmetric,
    nwake,
    surface_id,
    wake_finite_core,
    wake_shedding_locations,
    trailing_vortices,
    xhat,
    interaction_id = surface_id,
    interaction::Bool = true,
    threaded::Bool = true,
    chord_force_buffers = nothing,
    span_force_buffers = nothing,
    unsteady_force_buffers = nothing,
)
    number_of_surfaces = length(surfaces)
    scalar_type = eltype(circulation)
    zero_force = zero(SVector{3, scalar_type})

    chord_segment_forces = isnothing(chord_force_buffers) ?
        Vector{Matrix{SVector{3, scalar_type}}}(undef, number_of_surfaces) :
        chord_force_buffers
    span_segment_forces = isnothing(span_force_buffers) ?
        Vector{Matrix{SVector{3, scalar_type}}}(undef, number_of_surfaces) :
        span_force_buffers
    unsteady_forces = isnothing(unsteady_force_buffers) ?
        Vector{Matrix{SVector{3, scalar_type}}}(undef, number_of_surfaces) :
        unsteady_force_buffers
    length(chord_segment_forces) == number_of_surfaces ||
        throw(DimensionMismatch("one chord-force buffer is required per surface"))
    length(span_segment_forces) == number_of_surfaces ||
        throw(DimensionMismatch("one span-force buffer is required per surface"))
    length(unsteady_forces) == number_of_surfaces ||
        throw(DimensionMismatch("one unsteady-force buffer is required per surface"))

    panel_counts = length.(surfaces)
    circulation_ends = cumsum(panel_counts)
    circulation_starts = circulation_ends .- panel_counts .+ 1
    surface_circulations = [
        view(circulation, circulation_starts[index]:circulation_ends[index])
        for index in eachindex(surfaces)
    ]
    wake_is_active = nwake .> 0
    shedding_locations = isnothing(wake_shedding_locations) ?
        fill(nothing, number_of_surfaces) : wake_shedding_locations
    same_interaction_group = [
        interaction || interaction_id[receiving] == interaction_id[sending]
        for receiving in eachindex(surfaces), sending in eachindex(surfaces)
    ]
    different_surface_id = [
        surface_id[receiving] != surface_id[sending]
        for receiving in eachindex(surfaces), sending in eachindex(surfaces)
    ]

    for receiving_surface in eachindex(surfaces)
        receiving = surfaces[receiving_surface]
        number_chordwise, number_spanwise = size(receiving)
        circulation_offset = circulation_starts[receiving_surface] - 1

        expected_chord_size = (number_chordwise, number_spanwise + 1)
        expected_span_size = (number_chordwise + 1, number_spanwise)
        expected_unsteady_size = (number_chordwise, number_spanwise)

        current_chord = isnothing(chord_force_buffers) ?
            fill(zero_force, expected_chord_size) :
            chord_segment_forces[receiving_surface]
        current_span = isnothing(span_force_buffers) ?
            fill(zero_force, expected_span_size) :
            span_segment_forces[receiving_surface]
        current_unsteady = isnothing(unsteady_force_buffers) ?
            fill(zero_force, expected_unsteady_size) :
            unsteady_forces[receiving_surface]
        size(current_chord) == expected_chord_size ||
            throw(DimensionMismatch("chord-force buffer has the wrong size"))
        size(current_span) == expected_span_size ||
            throw(DimensionMismatch("span-force buffer has the wrong size"))
        size(current_unsteady) == expected_unsteady_size ||
            throw(DimensionMismatch("unsteady-force buffer has the wrong size"))
        fill!(current_chord, zero_force)
        fill!(current_span, zero_force)
        fill!(current_unsteady, zero_force)

        function calculate_panel_forces!(linear_index)
            panel_index = CartesianIndices(receiving)[linear_index]
            chordwise_index, spanwise_index = Tuple(panel_index)
            panel = receiving[panel_index]
            global_index = circulation_offset + linear_index

            span_location = top_center(panel)
            span_motion = isnothing(Vh) ? nothing :
                Vh[receiving_surface][panel_index]
            span_velocity = _legacy_imperial_base_velocity(
                span_location,
                reference,
                freestream,
                additional_velocity,
                span_motion,
            )
            span_velocity += _legacy_imperial_induced_velocity(
                span_location,
                panel_index,
                :spanwise,
                surfaces,
                wakes,
                surface_circulations,
                receiving_surface;
                symmetric,
                nwake,
                wake_finite_core,
                wake_is_active,
                shedding_locations,
                trailing_vortices,
                xhat,
                same_interaction_group,
                different_surface_id,
            )
            spanwise_jump = chordwise_index == 1 ?
                circulation[global_index] :
                circulation[global_index] - circulation[global_index - 1]
            current_span[chordwise_index, spanwise_index] =
                reference.rho * spanwise_jump *
                cross(span_velocity, top_vector(panel))

            left_location = left_center(panel)
            left_motion = isnothing(Vv) ? nothing :
                Vv[receiving_surface][chordwise_index, spanwise_index]
            left_velocity = _legacy_imperial_base_velocity(
                left_location,
                reference,
                freestream,
                additional_velocity,
                left_motion,
            )
            left_velocity += _legacy_imperial_induced_velocity(
                left_location,
                panel_index,
                :chordwise,
                surfaces,
                wakes,
                surface_circulations,
                receiving_surface;
                symmetric,
                nwake,
                wake_finite_core,
                wake_is_active,
                shedding_locations,
                trailing_vortices,
                xhat,
                same_interaction_group,
                different_surface_id,
            )
            chordwise_jump = spanwise_index == 1 ?
                circulation[global_index] :
                circulation[global_index] -
                    circulation[global_index - number_chordwise]
            current_chord[chordwise_index, spanwise_index] =
                reference.rho * chordwise_jump *
                cross(left_velocity, left_vector(panel))

            if spanwise_index == number_spanwise
                right_location = right_center(panel)
                right_motion = isnothing(Vv) ? nothing :
                    Vv[receiving_surface][chordwise_index, spanwise_index + 1]
                right_velocity = _legacy_imperial_base_velocity(
                    right_location,
                    reference,
                    freestream,
                    additional_velocity,
                    right_motion,
                )
                right_velocity += _legacy_imperial_induced_velocity(
                    right_location,
                    CartesianIndex(chordwise_index, spanwise_index + 1),
                    :chordwise,
                    surfaces,
                    wakes,
                    surface_circulations,
                    receiving_surface;
                    symmetric,
                    nwake,
                    wake_finite_core,
                    wake_is_active,
                    shedding_locations,
                    trailing_vortices,
                    xhat,
                    same_interaction_group,
                    different_surface_id,
                )
                current_chord[chordwise_index, spanwise_index + 1] =
                    reference.rho * circulation[global_index] *
                    cross(right_velocity, right_vector(panel))
            end

            if !isnothing(dΓdt)
                approximate_area =
                    norm(top_vector(panel)) * norm(left_vector(panel))
                current_unsteady[chordwise_index, spanwise_index] =
                    reference.rho * approximate_area * normal(panel) *
                    dΓdt[global_index]
            end
            return nothing
        end

        if threaded && Threads.nthreads() > 1
            Threads.@threads :static for linear_index in eachindex(receiving)
                calculate_panel_forces!(linear_index)
            end
        else
            for linear_index in eachindex(receiving)
                calculate_panel_forces!(linear_index)
            end
        end

        chord_segment_forces[receiving_surface] = current_chord
        span_segment_forces[receiving_surface] = current_span
        unsteady_forces[receiving_surface] = current_unsteady
    end

    return properties, chord_segment_forces, span_segment_forces, unsteady_forces
end

"""
    legacy_near_field_forces!(args...; kwargs...)

Backward-compatible name for [`legacy_imperial_segment_forces!`](@ref). The
former non-Imperial panel-property reconstruction has been removed.
"""
legacy_near_field_forces!(args...; kwargs...) =
    legacy_imperial_segment_forces!(args...; kwargs...)
