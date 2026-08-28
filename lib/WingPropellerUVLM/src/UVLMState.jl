"""
    UVLMSnapshot

Deep copy of the history variables required to repeat a UVLM time-step trial
without advancing circulation, geometry, or wake state more than once.
"""
struct UVLMSnapshot{TW,TΓ,TDΓ,TΓT,TS,TWS,TF}
    wakes::TW
    Γ::TΓ
    dΓ::TDΓ
    dΓdt::TΓT
    surfaces::TS
    previousSurfaces::TS
    wakeSheddingLocations::TWS
    nwake::Vector{Int}
    freestream::TF
end

"""
    snapshot_uvlm(system::System)

Capture the mutable aerodynamic history at the beginning of a physical time
step. Restore this snapshot before every rejected or repeated coupling trial.
"""
function snapshot_uvlm(system::System)
    return UVLMSnapshot(
        deepcopy(system.wakes),
        copy(system.Γ),
        tuple((copy(value) for value in system.dΓ)...),
        copy(system.dΓdt),
        deepcopy(system.surfaces),
        deepcopy(system.previous_surfaces),
        deepcopy(system.wake_shedding_locations),
        copy(system.nwake),
        deepcopy(system.freestream[]),
    )
end

function restore_nested_arrays!(destination,source)
    @assert length(destination) == length(source)
    for i in eachindex(destination)
        @assert size(destination[i]) == size(source[i])
        copyto!(destination[i],source[i])
    end
    return destination
end

"""
    restore_uvlm!(system::System,snapshot::UVLMSnapshot)

Restore a UVLM system in place. Array identities are retained so other objects
may safely hold references to the system storage.
"""
function restore_uvlm!(system::System,snapshot::UVLMSnapshot)
    restore_nested_arrays!(system.wakes,snapshot.wakes)
    restore_nested_arrays!(system.surfaces,snapshot.surfaces)
    restore_nested_arrays!(system.previous_surfaces,snapshot.previousSurfaces)

    system.Γ .= snapshot.Γ
    system.dΓdt .= snapshot.dΓdt
    for i in eachindex(system.dΓ)
        system.dΓ[i] .= snapshot.dΓ[i]
    end

    @assert length(system.wake_shedding_locations) ==
        length(snapshot.wakeSheddingLocations)
    for i in eachindex(system.wake_shedding_locations)
        system.wake_shedding_locations[i] .=
            snapshot.wakeSheddingLocations[i]
    end

    system.nwake .= snapshot.nwake
    system.freestream[] = snapshot.freestream
    return system
end

"""
    advance_uvlm_trial!(system,snapshot,freestream,Δt; kwargs...)

Restore `snapshot` and propagate one trial from the accepted state. Set
`advanceWake=false` to leave the wake unchanged during a partitioned-coupling
trial; after convergence, call `advance_wake!` once using the accepted state.
The caller separately owns and commits the active wake-row counter.
"""
function advance_uvlm_trial!(system::System,snapshot::UVLMSnapshot,
    freestream::Freestream,Δt::Real;
    additionalVelocity=nothing,
    repeatedPoints=repeated_trailing_edge_points(system.surfaces),
    activeWakeRows=system.nwake,
    η::Real=0.1,
    calculateInfluenceMatrix::Bool=true,
    nearFieldAnalysis::Bool=true,
    nearFieldForceFunction=near_field_forces!,
    derivatives::Bool=false,
    interactionID=system.surface_id,
    interaction::Bool=true,
    advanceWake::Bool=true)

    @assert Δt > 0
    @assert length(activeWakeRows) == length(system.surfaces)
    restore_uvlm!(system,snapshot)

    propagate_system!(
        system,
        freestream,
        Δt;
        additional_velocity=additionalVelocity,
        repeated_points=repeatedPoints,
        nwake=activeWakeRows,
        eta=η,
        calculate_influence_matrix=calculateInfluenceMatrix,
        near_field_analysis=nearFieldAnalysis,
        near_field_force_function=nearFieldForceFunction,
        derivatives=derivatives,
        interaction_id=interactionID,
        interaction=interaction,
        advance_wake=advanceWake,
    )

    return system
end

"""
    commit_wake_rows!(activeWakeRows,maximumWakeRows)

Advance each active wake length once, without exceeding its allocated maximum.
"""
function commit_wake_rows!(activeWakeRows::Vector{Int},
    maximumWakeRows::AbstractVector{<:Integer})

    @assert length(activeWakeRows) == length(maximumWakeRows)
    for i in eachindex(activeWakeRows)
        @assert 0 <= activeWakeRows[i] <= maximumWakeRows[i]
        activeWakeRows[i] = min(activeWakeRows[i]+1,maximumWakeRows[i])
    end
    return activeWakeRows
end
