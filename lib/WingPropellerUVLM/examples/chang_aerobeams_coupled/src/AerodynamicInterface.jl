# Additive nonlinear-beam/UVLM interface. All positions, spatial rotations,
# forces and moments in this file are in the right-handed UVLM global frame U.
const CHANG_DIRECTORS_TO_U = [0.0 1.0 0.0; 1.0 0.0 0.0; 0.0 0.0 -1.0]

function _ca_validate_poses(poses, span_nodes)
    size(poses.positions) == (3, length(span_nodes)) ||
        throw(DimensionMismatch("poses.positions must be 3 by structural node count"))
    length(poses.rotations) == length(span_nodes) ||
        throw(DimensionMismatch("One director rotation is required per structural node"))
    length(span_nodes) >= 2 && all(>(0), diff(span_nodes)) ||
        throw(ArgumentError("Structural span stations must be strictly increasing"))
    all(isfinite, poses.positions) || throw(ArgumentError("Nonfinite wing positions"))
    for rotation in poses.rotations
        size(rotation) == (3, 3) && all(isfinite, rotation) &&
            isapprox(rotation' * rotation, I; atol=1e-8, rtol=1e-8) &&
            isapprox(det(rotation), 1.0; atol=1e-8, rtol=1e-8) ||
            throw(ArgumentError("Wing directors must be proper 3 by 3 rotations"))
    end
    return nothing
end

function _ca_span_weights(stations, query)
    tolerance = 1e-10 * max(1.0, abs(last(stations) - first(stations)))
    first(stations)-tolerance <= query <= last(stations)+tolerance ||
        throw(ArgumentError("Aerodynamic span station is outside the structural mesh"))
    q = clamp(query, first(stations), last(stations))
    left = clamp(searchsortedlast(stations, q), 1, length(stations)-1)
    fraction = (q-stations[left]) / (stations[left+1]-stations[left])
    return left, left+1, 1.0-fraction, fraction
end

function _ca_interp(stations, values, query)
    left, right, wl, wr = _ca_span_weights(stations, query)
    return wl*values[left] + wr*values[right]
end

function _ca_reference_grids(reference)
    p = reference.parameters
    wing, _ = UVLM.wing_to_grid(p.xle, p.yle, p.zle, p.chord_geo,
        p.theta_geo, p.phi_geo, p.ns_wing, p.nc_wing;
        mirror=p.mirror_wing, spacing_s=p.spacing_s_wing, spacing_c=p.spacing_c_wing)
    blades = UVLM.generate_propeller_blades_grid(p.R_prop, p.c_prop,
        p.ns_prop, p.nc_prop, p.blade_twists_prop, p.Nb_prop)
    return wing, blades
end

function _ca_attachment_nodes(reference, span_nodes)
    p = reference.parameters
    nodes = Int[]
    for reference_node in p.prop_attach_nodes
        station = p.span_nodes[reference_node]
        node = findfirst(x -> isapprox(x, station; atol=1e-10, rtol=1e-10), span_nodes)
        isnothing(node) && throw(ArgumentError(
            "The structural mesh must contain each physical propeller attachment station ($station m)"))
        push!(nodes, node)
    end
    return nodes
end

function _ca_wing_geometry(reference, poses, structural_span_nodes)
    p = reference.parameters
    _ca_validate_poses(poses, structural_span_nodes)
    wing_reference, blade_reference = _ca_reference_grids(reference)
    count = length(structural_span_nodes)
    nchord, nspan = size(wing_reference, 2), size(wing_reference, 3)
    weights = zeros(count, nspan)
    chord_fractions = zeros(nchord, nspan)
    chord_nodes = [_ca_interp(p.span_nodes, p.chord, y) for y in structural_span_nodes]
    wing_grid = zeros(3, nchord, nspan)
    for j in 1:nspan
        station = wing_reference[2, 1, j]
        left, right, wl, wr = _ca_span_weights(structural_span_nodes, station)
        weights[left,j], weights[right,j] = wl, wr
        chord_here = _ca_interp(p.span_nodes, p.chord, station)
        xle_here = _ca_interp(p.span_nodes, p.xle_distribution, station)
        for k in 1:nchord
            fraction = (wing_reference[1,k,j] - xle_here) / chord_here
            chord_fractions[k,j] = fraction
            for (node, weight) in ((left,wl), (right,wr))
                offset = (fraction-reference.aerodynamic_options.elastic_axis_fraction) *
                    chord_nodes[node] * poses.rotations[node][:,2]
                wing_grid[:,k,j] .+= weight .* (poses.positions[:,node] + offset)
            end
        end
    end
    # grid_to_surface_panels moves each bound row to the quarter-panel point,
    # but the final row remains the physical trailing edge.
    vortex_fractions = copy(chord_fractions)
    for k in 1:nchord-1
        vortex_fractions[k,:] .= 0.75 .* chord_fractions[k,:] .+ 0.25 .* chord_fractions[k+1,:]
    end
    geometry = (;weights,chord_nodes,chord_fractions,vortex_fractions,
        structural_span_nodes=collect(structural_span_nodes))
    return wing_grid, blade_reference, geometry
end

"""
    coupled_grids(reference, poses, eta, time; structural_span_nodes=reference.parameters.span_nodes)

Return `(grids, geometry)` for a nonlinear wing and the Chang relative pitch/yaw
coordinates `eta = [pitch1,yaw1,pitch2,yaw2,...]`. `poses.positions` contains
elastic-axis positions in U; each `poses.rotations[i]` maps the initial
`[span,chord,down]` directors to U. Physical rotor hubs retain the full pylon
length while their relative translation retains Chang's assumed-mode load arm.

Wing aerodynamic stations may differ from structural stations. Their vertex
positions are weighted sums of actual structural-node chord-point positions;
the load operator below is the exact transpose of this motion map. Each rotor
attachment must be an actual structural node. Geometry carries the *vortex*
quarter-panel chord fractions used for the paired UVLM vertex forces.
"""
function coupled_grids(reference, poses, eta::AbstractVector, time::Real;
    structural_span_nodes=reference.parameters.span_nodes)
    p = reference.parameters
    length(eta) == 2p.Npropellers || throw(DimensionMismatch("Two coordinates per rotor are required"))
    all(isfinite, eta) && isfinite(time) || throw(ArgumentError("Nonfinite rotor coordinates or time"))
    wing_grid,blade_reference,wing_geometry = _ca_wing_geometry(reference,poses,structural_span_nodes)
    grids = Array{Float64,3}[wing_grid]
    attachment_nodes = _ca_attachment_nodes(reference, structural_span_nodes)
    hubs, pivots, centers = Vector{Float64}[], Vector{Float64}[], Vector{Float64}[]
    pitch_axes, yaw_axes = Vector{Float64}[], Vector{Float64}[]
    propeller_surfaces = Vector{Int}[]
    physical_offset = collect(reference.aerodynamic_options.physical_hub_center_A)
    modal_offset = collect(reference.aerodynamic_options.load_center_A)
    spin = UVLM.RotationMatrix(-p.Ω*time, 1)
    for rotor in 1:p.Npropellers
        node = attachment_nodes[rotor]
        deformation = poses.rotations[node] * CHANG_DIRECTORS_TO_U'
        pitch = UVLM.RotationMatrix(eta[2rotor-1], 2)
        yaw = UVLM.RotationMatrix(-eta[2rotor], 3)
        whirl = pitch*yaw
        pivot = collect(poses.positions[:,node])
        center = pivot + deformation*(physical_offset-modal_offset)
        hub = center + deformation*whirl*modal_offset
        push!(pivots,pivot); push!(hubs,hub); push!(centers,center)
        push!(pitch_axes,collect(deformation*[0.0,1.0,0.0]))
        push!(yaw_axes,collect(-deformation*pitch*[0.0,0.0,1.0]))
        surface_indices = Int[]
        for blade in blade_reference
            grid = similar(blade, Float64)
            for index in CartesianIndices((size(blade,2),size(blade,3)))
                k,j = Tuple(index)
                grid[:,k,j] .= center + deformation*whirl*(modal_offset+spin*blade[:,k,j])
            end
            push!(grids,grid); push!(surface_indices,length(grids))
        end
        push!(propeller_surfaces,surface_indices)
    end
    geometry = (; wing_geometry...,
        attachment_nodes, hubs, pivots, centers, pitch_axes, yaw_axes,
        propeller_surfaces)
    return grids, geometry
end

"""
    aerodynamic_loads(reference, poses, eta, time, vertex_forces, vertex_positions; ...)

Reduce paired dimensional UVLM vortex-vertex forces/application points to
physical wing-node wrenches in U and exact relative-rotor generalized forces.
Wing rows 1:3 are spatial forces; rows 4:6 are spatial moments about the current
elastic axis. They are not moments conjugate to AeroBeams rotation parameters.
The structural adapter performs that final parameterization-specific conversion.
"""
function _ca_wing_loads(reference, poses, forces, positions, geometry)
    wing_loads = zeros(6,size(poses.positions,2))
    size(forces) == size(positions) == size(geometry.vortex_fractions) ||
        throw(DimensionMismatch("Wing loads must have the vortex-grid topology"))
    for index in CartesianIndices(forces)
        k,j = Tuple(index)
        force = forces[index]
        for node in axes(geometry.weights,1)
            weight = geometry.weights[node,j]
            iszero(weight) && continue
            offset = (geometry.vortex_fractions[k,j] -
                reference.aerodynamic_options.elastic_axis_fraction) *
                geometry.chord_nodes[node] * poses.rotations[node][:,2]
            wing_loads[1:3,node] .+= weight .* force
            wing_loads[4:6,node] .+= weight .* cross(offset,force)
        end
    end
    return wing_loads
end

function aerodynamic_loads(reference, poses, eta, time, vertex_forces, vertex_positions;
    geometry=nothing, structural_span_nodes=reference.parameters.span_nodes)
    if isnothing(geometry)
        _, geometry = coupled_grids(reference,poses,eta,time;structural_span_nodes)
    end
    length(vertex_forces) == length(vertex_positions) == 1 +
        reference.parameters.Npropellers*reference.parameters.Nb_prop ||
        throw(DimensionMismatch("Expected one wing surface and every rotor blade"))
    wing_loads = _ca_wing_loads(reference,poses,vertex_forces[1],vertex_positions[1],geometry)
    rotor_loads = zeros(length(eta))
    hub_forces = zeros(3,reference.parameters.Npropellers)
    hub_moments = zeros(3,reference.parameters.Npropellers)
    for rotor in eachindex(geometry.attachment_nodes)
        force, moment = zeros(3), zeros(3)
        hub, pivot = geometry.hubs[rotor], geometry.pivots[rotor]
        for surface in geometry.propeller_surfaces[rotor]
            size(vertex_forces[surface]) == size(vertex_positions[surface]) ||
                throw(DimensionMismatch("Rotor force and position topologies differ"))
            for index in eachindex(vertex_forces[surface])
                f = vertex_forces[surface][index]
                force .+= f
                moment .+= cross(vertex_positions[surface][index]-hub,f)
            end
        end
        node = geometry.attachment_nodes[rotor]
        wing_loads[1:3,node] .+= force
        wing_loads[4:6,node] .+= moment + cross(hub-pivot,force)
        # Relative coordinates rotate about the effective Chang modal center.
        modal_moment = moment + cross(hub-geometry.centers[rotor],force)
        rotor_loads[2rotor-1] = dot(geometry.pitch_axes[rotor],modal_moment)
        rotor_loads[2rotor] = dot(geometry.yaw_axes[rotor],modal_moment)
        hub_forces[:,rotor] .= force
        hub_moments[:,rotor] .= moment
    end
    return (;wing_loads,rotor_loads,hub_forces,hub_moments)
end

coupled_aerodynamic_loads(reference,poses,eta,time,loads;kwargs...) =
    aerodynamic_loads(reference,poses,eta,time,loads.forces,loads.positions;kwargs...)

"""
    build_aerodynamic(reference, poses, eta, time_vector; ...)

Create the multi-surface `UVLMDynamicProblem` at the initial physical time.
Call `coupled_grids` and pass its grids to `UVLM.evaluate_trial!` during FSI.
Unique physical surface IDs are retained separately from aerodynamic interaction
groups in all three UVLM transaction buffers, without modifying the backend.
"""
function build_aerodynamic(reference, poses, eta, time_vector;
    structural_span_nodes=reference.parameters.span_nodes, kwargs...)
    times = collect(Float64,time_vector)
    length(times) >= 2 || throw(ArgumentError("At least two aerodynamic times are required"))
    grids,_ = coupled_grids(reference,poses,eta,first(times);structural_span_nodes)
    return _ca_build_aerodynamic(reference,grids,times;kwargs...)
end

function _ca_build_aerodynamic(reference, grids, time_vector; solver=nothing,
    maximum_wake_rows=nothing, tracking_time_steps=true, tracking_frequency=1,
    save_geometry=false)
    p = reference.parameters
    reference.config.simulation.near_field_force_model == :imperial ||
        throw(ArgumentError("The dynamic aerodynamic adapter requires the corrected :imperial load model"))
    times = collect(Float64,time_vector)
    length(times) >= 2 || throw(ArgumentError("At least two aerodynamic times are required"))
    groups = vcat(1,repeat(collect(2:p.Npropellers+1);inner=p.Nb_prop))
    symmetric = vcat(p.symmetric_wing,fill(false,p.Npropellers*p.Nb_prop))
    model = UVLM.create_UVLMModel(;surfaces=grids,reference=p.ref,
        name="Chang nonlinear AeroBeams wing with 2DOF propellers",
        symmetric=false,interactionGroups=groups)
    # A vector Bool broadcast in the existing convenience constructor becomes a
    # BitVector, but UVLMModel stores Vector{Bool}. Allocate through its scalar
    # path, then set the wing-only symmetry locally without changing core code.
    model.symmetric .= symmetric
    if isnothing(solver)
        core = reference.config.aerodynamic.core_radius_m
        isnothing(core) && throw(ArgumentError(
            "The UVLM dynamic facade requires a fixed core radius; pass solver explicitly"))
        solver = UVLM.create_UVLMSolver(;coreRadius=core,
            interaction=reference.config.simulation.interaction_on)
    end
    rows = isnothing(maximum_wake_rows) ? vcat(
        reference.aerodynamic_options.maximum_wake_rows_wing,
        fill(reference.aerodynamic_options.maximum_wake_rows_propeller,p.Npropellers*p.Nb_prop)) :
        maximum_wake_rows
    point = UVLM.create_OperatingPoint(;airspeed=p.Vinf,density=p.ref.rho,
        angleOfAttack=p.fs.alpha,sideslip=p.fs.beta,rotationRPM=p.Ω*60/(2pi))
    problem = UVLM.create_UVLMDynamicProblem(;model,aeroSolver=solver,
        operatingPoint=point,timeVector=times,maximumWakeRows=rows,
        trackingTimeSteps=tracking_time_steps,trackingFrequency=tracking_frequency,
        saveGeometry=save_geometry)
    # The public dynamic facade currently uses group IDs as physical IDs.
    # Correct it locally: finite-core decisions use unique physical IDs, while
    # evaluate/commit still receive model.interactionGroups for interference.
    for system in (problem.state.system,problem.workspace.system,problem.workspace.commitSystem)
        system.surface_id .= eachindex(system.surfaces)
    end
    return problem
end

"""
    pylon_coupled_grids(reference, wing_poses, hub_poses, time; ...)

Wing and rigid-blade grids for the beam-pylon variation. Each column of
`hub_poses.positions` is the *physical* pylon-tip/rotor-hub position in U.
`hub_poses.rotations[i]` maps the reference rotor aerodynamic axes to global U;
it is I for an initially aligned rotor. Convert actual pylon section directors
with `R_tip * R_tip_initial'` before calling. There are no relative two-DOF
coordinates, artificial modal centers, or assumed-mode half-length load arms.
"""
function pylon_coupled_grids(reference, wing_poses, hub_poses, time::Real;
    structural_span_nodes=reference.parameters.span_nodes)
    p = reference.parameters
    size(hub_poses.positions) == (3,p.Npropellers) &&
        length(hub_poses.rotations) == p.Npropellers ||
        throw(DimensionMismatch("One physical hub position and rotation per propeller are required"))
    isfinite(time) && all(isfinite,hub_poses.positions) ||
        throw(ArgumentError("Nonfinite pylon hub positions or time"))
    for rotation in hub_poses.rotations
        size(rotation) == (3,3) && all(isfinite,rotation) &&
            isapprox(rotation'*rotation,I;atol=1e-8,rtol=1e-8) &&
            isapprox(det(rotation),1.0;atol=1e-8,rtol=1e-8) ||
            throw(ArgumentError("Hub reference-to-global rotations must be proper rotations"))
    end
    wing_grid,blade_reference,wing_geometry = _ca_wing_geometry(reference,wing_poses,structural_span_nodes)
    grids = Array{Float64,3}[wing_grid]
    spin = UVLM.RotationMatrix(-p.Ω*time,1)
    hubs = [collect(hub_poses.positions[:,i]) for i in 1:p.Npropellers]
    propeller_surfaces = Vector{Int}[]
    for rotor in 1:p.Npropellers
        indices = Int[]
        for blade in blade_reference
            grid = similar(blade,Float64)
            for index in CartesianIndices((size(blade,2),size(blade,3)))
                k,j = Tuple(index)
                grid[:,k,j] .= hubs[rotor] + hub_poses.rotations[rotor]*spin*blade[:,k,j]
            end
            push!(grids,grid); push!(indices,length(grids))
        end
        push!(propeller_surfaces,indices)
    end
    return grids,(;wing_geometry...,hubs,propeller_surfaces,support=:beam_pylon)
end

"""
    pylon_aerodynamic_loads(reference, wing_poses, hub_poses, time, forces, positions; ...)

Return wing-only aerodynamic wrenches plus `pylon_loads::6×Nrotor`, the physical
tip force and moment about each hub. Apply the hub wrench once to the pylon tip;
its reaction reaches the wing through the assembled nonlinear beam joint. No
direct rotor wrench is added to the wing and `rotor_loads` is empty.
"""
function pylon_aerodynamic_loads(reference, wing_poses, hub_poses, time,
    vertex_forces, vertex_positions;
    geometry=nothing, structural_span_nodes=reference.parameters.span_nodes)
    if isnothing(geometry)
        _,geometry = pylon_coupled_grids(reference,wing_poses,hub_poses,time;structural_span_nodes)
    end
    p = reference.parameters
    length(vertex_forces) == length(vertex_positions) == 1+p.Npropellers*p.Nb_prop ||
        throw(DimensionMismatch("Expected one wing surface and every rotor blade"))
    wing_loads = _ca_wing_loads(reference,wing_poses,vertex_forces[1],vertex_positions[1],geometry)
    pylon_loads = zeros(6,p.Npropellers)
    for rotor in 1:p.Npropellers
        for surface in geometry.propeller_surfaces[rotor]
            size(vertex_forces[surface]) == size(vertex_positions[surface]) ||
                throw(DimensionMismatch("Rotor force and position topologies differ"))
            for index in eachindex(vertex_forces[surface])
                force = vertex_forces[surface][index]
                pylon_loads[1:3,rotor] .+= force
                pylon_loads[4:6,rotor] .+=
                    cross(vertex_positions[surface][index]-geometry.hubs[rotor],force)
            end
        end
    end
    return (;wing_loads,pylon_loads,rotor_loads=Float64[],
        hub_forces=copy(pylon_loads[1:3,:]),hub_moments=copy(pylon_loads[4:6,:]))
end

pylon_aerodynamic_loads(reference,wing_poses,hub_poses,time,loads;kwargs...) =
    pylon_aerodynamic_loads(reference,wing_poses,hub_poses,time,loads.forces,loads.positions;kwargs...)

"""Build the UVLM transaction problem using actual beam-pylon hub poses."""
function pylon_build_aerodynamic(reference, wing_poses, hub_poses, time_vector;
    structural_span_nodes=reference.parameters.span_nodes,kwargs...)
    times = collect(Float64,time_vector)
    length(times) >= 2 || throw(ArgumentError("At least two aerodynamic times are required"))
    grids,_ = pylon_coupled_grids(reference,wing_poses,hub_poses,first(times);structural_span_nodes)
    return _ca_build_aerodynamic(reference,grids,times;kwargs...)
end
