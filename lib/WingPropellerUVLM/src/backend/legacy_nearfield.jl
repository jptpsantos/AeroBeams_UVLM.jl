# Corrected legacy segment-force formulation for the wing-propeller
# aeroelastic model. Only the block derived from Imperial College London's C++
# UVLM implementation is kept in this file.

function _legacy_imperial_base_velocity(
    location,
    reference,
    freestream,
    additional_velocity,
    surface_motion
)
    velocity = freestream_velocity(freestream)
    velocity += rotational_velocity(location,freestream,reference)
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
    surfaceCirculations,
    receiving_surface;
    symmetric,
    nwake,
    wake_finite_core,
    wake_is_active,
    shedding_locations,
    trailing_vortices,
    xhat,
    same_interaction_group,
    different_surface_id
)
    segment_kind in (:spanwise,:chordwise) ||
        throw(ArgumentError("segment_kind must be :spanwise or :chordwise"))

    velocity = zero(location)
    for sendingSurface in eachindex(surfaces)
        if !same_interaction_group[receiving_surface,sendingSurface]
            continue
        end

        sending = surfaces[sendingSurface]
        sendingCirculation = surfaceCirculations[sendingSurface]
        idsAreDifferent =
            different_surface_id[receiving_surface,sendingSurface]

        if receiving_surface == sendingSurface
            if segment_kind === :spanwise
                velocity += induced_velocity(
                    panel_index,
                    sending,
                    sendingCirculation;
                    finite_core=idsAreDifferent,
                    wake_shedding_locations=shedding_locations[sendingSurface],
                    symmetric=symmetric[sendingSurface],
                    trailing_vortices=trailing_vortices[sendingSurface] &&
                                      !wake_is_active[sendingSurface],
                    xhat=xhat
                )
            else
                chordwiseIndex,spanwiseSegmentIndex = Tuple(panel_index)
                skip_left =
                    spanwiseSegmentIndex <= size(sending,2) ?
                    (CartesianIndex(chordwiseIndex,spanwiseSegmentIndex),) : ()
                skip_right = spanwiseSegmentIndex > 1 ?
                             (CartesianIndex(chordwiseIndex,spanwiseSegmentIndex - 1),) : ()
                velocity += induced_velocity(
                    location,
                    sending,
                    sendingCirculation;
                    finite_core=idsAreDifferent,
                    wake_shedding_locations=shedding_locations[sendingSurface],
                    symmetric=symmetric[sendingSurface],
                    trailing_vortices=trailing_vortices[sendingSurface] &&
                                      !wake_is_active[sendingSurface],
                    xhat=xhat,
                    skip_left,
                    skip_right
                )
            end
        else
            velocity += induced_velocity(location,sending,sendingCirculation;
                finite_core=idsAreDifferent,
                wake_shedding_locations=shedding_locations[sendingSurface],
                symmetric=symmetric[sendingSurface],
                trailing_vortices=trailing_vortices[sendingSurface] &&
                                  !wake_is_active[sendingSurface],
                xhat=xhat)
        end

        if wake_is_active[sendingSurface]
            velocity += induced_velocity(location,wakes[sendingSurface];
                finite_core=wake_finite_core[sendingSurface] ||
                            idsAreDifferent,
                symmetric=symmetric[sendingSurface],
                nc=nwake[sendingSurface],
                trailing_vortices=trailing_vortices[sendingSurface],
                xhat=xhat)
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
    interaction_id=surface_id,
    interaction::Bool=true,
    threaded::Bool=true,
    chord_force_buffers=nothing,
    span_force_buffers=nothing,
    unsteady_force_buffers=nothing
)
    nSurfaces = length(surfaces)
    scalarType = eltype(circulation)
    zero_force = zero(SVector{3,scalarType})

    chordSegmentForces =
        isnothing(chord_force_buffers) ?
        Vector{Matrix{SVector{3,scalarType}}}(undef,nSurfaces) :
        chord_force_buffers
    spanSegmentForces =
        isnothing(span_force_buffers) ?
        Vector{Matrix{SVector{3,scalarType}}}(undef,nSurfaces) :
        span_force_buffers
    unsteady_forces =
        isnothing(unsteady_force_buffers) ?
        Vector{Matrix{SVector{3,scalarType}}}(undef,nSurfaces) :
        unsteady_force_buffers
    length(chordSegmentForces) == nSurfaces ||
        throw(DimensionMismatch("one chord-force buffer is required per surface"))
    length(spanSegmentForces) == nSurfaces ||
        throw(DimensionMismatch("one span-force buffer is required per surface"))
    length(unsteady_forces) == nSurfaces ||
        throw(DimensionMismatch("one unsteady-force buffer is required per surface"))

    panelCounts = length.(surfaces)
    circulationEnds = cumsum(panelCounts)
    circulationStarts = circulationEnds .- panelCounts .+ 1
    surfaceCirculations = [
        view(circulation,circulationStarts[index]:circulationEnds[index])
        for index in eachindex(surfaces)
    ]
    wake_is_active = nwake .> 0
    shedding_locations = isnothing(wake_shedding_locations) ?
                         fill(nothing,nSurfaces) : wake_shedding_locations
    same_interaction_group = [
        interaction || interaction_id[receiving] == interaction_id[sending]
        for receiving in eachindex(surfaces),sending in eachindex(surfaces)
    ]
    different_surface_id = [
        surface_id[receiving] != surface_id[sending]
        for receiving in eachindex(surfaces),sending in eachindex(surfaces)
    ]

    for receiving_surface in eachindex(surfaces)
        receiving = surfaces[receiving_surface]
        numberChordwise,numberSpanwise = size(receiving)
        circulationOffset = circulationStarts[receiving_surface] - 1

        expectedChordSize = (numberChordwise,numberSpanwise + 1)
        expectedSpanSize = (numberChordwise + 1,numberSpanwise)
        expectedUnsteadySize = (numberChordwise,numberSpanwise)

        currentChord =
            isnothing(chord_force_buffers) ?
            fill(zero_force,expectedChordSize) :
            chordSegmentForces[receiving_surface]
        currentSpan =
            isnothing(span_force_buffers) ?
            fill(zero_force,expectedSpanSize) :
            spanSegmentForces[receiving_surface]
        currentUnsteady =
            isnothing(unsteady_force_buffers) ?
            fill(zero_force,expectedUnsteadySize) :
            unsteady_forces[receiving_surface]
        size(currentChord) == expectedChordSize ||
            throw(DimensionMismatch("chord-force buffer has the wrong size"))
        size(currentSpan) == expectedSpanSize ||
            throw(DimensionMismatch("span-force buffer has the wrong size"))
        size(currentUnsteady) == expectedUnsteadySize ||
            throw(DimensionMismatch("unsteady-force buffer has the wrong size"))
        fill!(currentChord,zero_force)
        fill!(currentSpan,zero_force)
        fill!(currentUnsteady,zero_force)

        function calculate_panel_forces!(linearIndex)
            panel_index = CartesianIndices(receiving)[linearIndex]
            chordwiseIndex,spanwiseIndex = Tuple(panel_index)
            panel = receiving[panel_index]
            globalIndex = circulationOffset + linearIndex

            spanLocation = top_center(panel)
            spanMotion = isnothing(Vh) ? nothing :
                         Vh[receiving_surface][panel_index]
            spanVelocity = _legacy_imperial_base_velocity(
                spanLocation,
                reference,
                freestream,
                additional_velocity,
                spanMotion
            )
            spanVelocity += _legacy_imperial_induced_velocity(
                spanLocation,
                panel_index,
                :spanwise,
                surfaces,
                wakes,
                surfaceCirculations,
                receiving_surface;
                symmetric,
                nwake,
                wake_finite_core,
                wake_is_active,
                shedding_locations,
                trailing_vortices,
                xhat,
                same_interaction_group,
                different_surface_id
            )
            spanwiseJump =
                chordwiseIndex == 1 ?
                circulation[globalIndex] :
                circulation[globalIndex] - circulation[globalIndex-1]
            currentSpan[chordwiseIndex,spanwiseIndex] =
                reference.rho * spanwiseJump *
                cross(spanVelocity,top_vector(panel))

            leftLocation = left_center(panel)
            leftMotion = isnothing(Vv) ? nothing :
                         Vv[receiving_surface][chordwiseIndex,spanwiseIndex]
            leftVelocity = _legacy_imperial_base_velocity(
                leftLocation,
                reference,
                freestream,
                additional_velocity,
                leftMotion
            )
            leftVelocity += _legacy_imperial_induced_velocity(
                leftLocation,
                panel_index,
                :chordwise,
                surfaces,
                wakes,
                surfaceCirculations,
                receiving_surface;
                symmetric,
                nwake,
                wake_finite_core,
                wake_is_active,
                shedding_locations,
                trailing_vortices,
                xhat,
                same_interaction_group,
                different_surface_id
            )
            chordwiseJump =
                spanwiseIndex == 1 ?
                circulation[globalIndex] :
                circulation[globalIndex] -
                circulation[globalIndex-numberChordwise]
            currentChord[chordwiseIndex,spanwiseIndex] =
                reference.rho * chordwiseJump *
                cross(leftVelocity,left_vector(panel))

            if spanwiseIndex == numberSpanwise
                rightLocation = right_center(panel)
                rightMotion = isnothing(Vv) ? nothing :
                              Vv[receiving_surface][chordwiseIndex,spanwiseIndex+1]
                rightVelocity = _legacy_imperial_base_velocity(
                    rightLocation,
                    reference,
                    freestream,
                    additional_velocity,
                    rightMotion
                )
                rightVelocity += _legacy_imperial_induced_velocity(
                    rightLocation,
                    CartesianIndex(chordwiseIndex,spanwiseIndex + 1),
                    :chordwise,
                    surfaces,
                    wakes,
                    surfaceCirculations,
                    receiving_surface;
                    symmetric,
                    nwake,
                    wake_finite_core,
                    wake_is_active,
                    shedding_locations,
                    trailing_vortices,
                    xhat,
                    same_interaction_group,
                    different_surface_id
                )
                currentChord[chordwiseIndex,spanwiseIndex+1] =
                    reference.rho * circulation[globalIndex] *
                    cross(rightVelocity,right_vector(panel))
            end

            if !isnothing(dΓdt)
                approximateArea =
                    norm(top_vector(panel)) * norm(left_vector(panel))
                currentUnsteady[chordwiseIndex,spanwiseIndex] =
                    reference.rho * approximateArea * normal(panel) *
                    dΓdt[globalIndex]
            end
            return nothing
        end

        if threaded && Threads.nthreads() > 1
            Threads.@threads :static for linearIndex in eachindex(receiving)
                calculate_panel_forces!(linearIndex)
            end
        else
            for linearIndex in eachindex(receiving)
                calculate_panel_forces!(linearIndex)
            end
        end

        chordSegmentForces[receiving_surface] = currentChord
        spanSegmentForces[receiving_surface] = currentSpan
        unsteady_forces[receiving_surface] = currentUnsteady
    end

    return properties,chordSegmentForces,spanSegmentForces,unsteady_forces
end

"""
    legacy_near_field_forces!(args...; kwargs...)

Backward-compatible name for [`legacy_imperial_segment_forces!`](@ref). The
former non-Imperial panel-property reconstruction has been removed.
"""
legacy_near_field_forces!(args...; kwargs...) =
    legacy_imperial_segment_forces!(args...; kwargs...)
