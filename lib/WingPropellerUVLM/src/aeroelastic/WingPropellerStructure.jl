# Configuration-independent structural assembly for a six-DOF 3D beam with
# directly attached two-DOF pitch/yaw propeller-nacelle systems.

const WING_DOFS_PER_NODE = 6

"""
    attachment_node_from_span_fraction(element_count, span_fraction)

Return the nearest node on a uniform `element_count`-element beam for a
physical propeller station expressed as a span fraction between zero and one.
The returned index follows Julia's one-based node numbering.
"""
function attachment_node_from_span_fraction(
    element_count::Integer,
    span_fraction::Real,
)
    element_count > 0 || throw(ArgumentError("element_count must be positive"))
    isfinite(span_fraction) && 0 <= span_fraction <= 1 || throw(ArgumentError(
        "span_fraction must be finite and between zero and one",
    ))
    return clamp(round(Int, span_fraction * element_count) + 1, 1, element_count + 1)
end

"""Return the six global wing coordinates belonging to `node`."""
function wing_node_dofs(node::Integer; dofs_per_node::Integer = WING_DOFS_PER_NODE)
    node > 0 || throw(ArgumentError("node must be positive"))
    dofs_per_node > 0 || throw(ArgumentError("dofs_per_node must be positive"))
    first_dof = dofs_per_node * (node - 1) + 1
    return first_dof:(first_dof + dofs_per_node - 1)
end

"""
Return a force/moment pair for a hub colocated with its structural attachment.
The moment is already expressed about that shared point, so no additional
`r x F` term is introduced.
"""
function colocated_hub_wrench(force, moment_about_hub, hub_position, attachment_position;
    tolerance::Real = 1.0e-12)
    isapprox(hub_position, attachment_position; atol = tolerance, rtol = 0.0) ||
        throw(ArgumentError("aerodynamic hub and structural attachment must be colocated"))
    return (; force, moment = moment_about_hub)
end

"""Add one six-component structural wrench only to the selected wing node."""
function add_direct_node_wrench!(wing_loads, wrench, node::Integer;
    dofs_per_node::Integer = WING_DOFS_PER_NODE)
    length(wrench) == dofs_per_node || throw(DimensionMismatch(
        "wrench length must equal dofs_per_node",
    ))
    node_dofs = wing_node_dofs(node; dofs_per_node)
    last(node_dofs) <= length(wing_loads) || throw(BoundsError(wing_loads, node_dofs))
    wing_loads[node_dofs] .+= wrench
    return wing_loads
end

"""
    full_beam_aerodynamic_kinematics(structural, free_state)

Reconstruct all wing-node coordinates and map the implemented structural basis
`(span, chord, down)` into the aerodynamic basis `(chord, span, up)`. Appended
propeller pitch/yaw coordinates are returned separately.
"""
function full_beam_aerodynamic_kinematics(structural, free_state)
    length(free_state) == structural.ndof_free || throw(DimensionMismatch(
        "free_state length must match the reduced structural system",
    ))
    global_state = zeros(promote_type(Float64, eltype(free_state)), size(structural.M_global, 1))
    global_state[structural.free_dofs] .= free_state
    wing_state = reshape(
        global_state[1:(WING_DOFS_PER_NODE * structural.node_count)],
        WING_DOFS_PER_NODE,
        structural.node_count,
    )
    propeller_state = global_state[(WING_DOFS_PER_NODE * structural.node_count + 1):end]
    return (;
        u_x = collect(wing_state[2, :]),
        u_y = collect(wing_state[1, :]),
        u_z = collect(-wing_state[3, :]),
        theta_x = collect(wing_state[5, :]),
        theta_y = collect(wing_state[4, :]),
        theta_z = collect(-wing_state[6, :]),
        propeller_pitch = collect(propeller_state[1:2:end]),
        propeller_yaw = collect(-propeller_state[2:2:end]),
    )
end

"""Return aerodynamic hub positions taken only from their attachment nodes."""
function direct_propeller_hub_positions(
    wing,
    propellers,
    structural,
    free_state,
)
    kinematics = full_beam_aerodynamic_kinematics(structural, free_state)
    span_nodes = wing.span_nodes
    geometry = wing.geometry
    span = last(span_nodes) - first(span_nodes)
    span > 0 || throw(ArgumentError("wing span must be positive"))
    fraction = (span_nodes .- first(span_nodes)) ./ span
    chord = geometry.root_chord .+
        (geometry.tip_chord - geometry.root_chord) .* fraction
    xle = geometry.xle_root .+ (geometry.xle_tip - geometry.xle_root) .* fraction
    elastic_axis_x = xle .+ geometry.elastic_axis_fraction .* chord
    return [
        begin
            node = Int(propeller.attachment_node)
            SVector(
                elastic_axis_x[node] + kinematics.u_x[node],
                span_nodes[node] + kinematics.u_y[node],
                kinematics.u_z[node],
            )
        end for propeller in propellers
    ]
end

_physical_field(data, name::Symbol, default) =
    hasproperty(data, name) ? getproperty(data, name) : default

function _linear_property_value(stations, values, position)
    length(stations) == length(values) || throw(DimensionMismatch(
        "spanwise property stations and values must have equal lengths",
    ))
    issorted(stations) || throw(ArgumentError("spanwise property stations must be sorted"))
    position <= first(stations) && return first(values)
    position >= last(stations) && return last(values)
    right = searchsortedfirst(stations, position)
    left = right - 1
    fraction = (position - stations[left]) / (stations[right] - stations[left])
    return (1 - fraction) * values[left] + fraction * values[right]
end

"""
    spanwise_property_value(property, position, span_nodes, element)

Evaluate a constant or spanwise structural property. A vector may contain one
value per element or one per node. A named tuple may use `(positions, values)`
in metres or `(eta, values)` in normalized span coordinates.
"""
function spanwise_property_value(property, position, span_nodes, element::Integer)
    property isa Number && return property
    if property isa AbstractVector
        element_count = length(span_nodes) - 1
        length(property) == element_count && return property[element]
        length(property) == length(span_nodes) && return _linear_property_value(
            span_nodes,
            property,
            position,
        )
        throw(DimensionMismatch(
            "a property vector must contain one value per element or per node",
        ))
    end
    if property isa NamedTuple
        values = _physical_field(property, :values, nothing)
        isnothing(values) && throw(ArgumentError(
            "a spanwise property named tuple must contain `values`",
        ))
        if hasproperty(property, :positions)
            return _linear_property_value(property.positions, values, position)
        elseif hasproperty(property, :eta)
            span = last(span_nodes) - first(span_nodes)
            span > 0 || throw(ArgumentError("wing span must be positive"))
            eta = (position - first(span_nodes)) / span
            return _linear_property_value(property.eta, values, eta)
        end
        throw(ArgumentError(
            "a spanwise property named tuple must contain `positions` or `eta`",
        ))
    end
    throw(ArgumentError("unsupported spanwise property representation"))
end

function _property_breakpoints(property, left, right, span_nodes)
    property isa NamedTuple || return Float64[]
    positions = if hasproperty(property, :positions)
        property.positions
    elseif hasproperty(property, :eta)
        first(span_nodes) .+ property.eta .* (last(span_nodes) - first(span_nodes))
    else
        return Float64[]
    end
    return Float64[position for position in positions if left < position < right]
end

function beam_strain_displacement_matrix(length::Real, xi::Real)
    length > 0 || throw(ArgumentError("beam element length must be positive"))
    B = zeros(4, 12)
    # Axial strain from spanwise translation.
    B[1, 1] = -1 / length
    B[1, 7] = 1 / length
    # In-plane curvature from chord translation and vertical-axis rotation.
    B[2, 2] = 6xi / length^2
    B[2, 6] = (3xi - 1) / length
    B[2, 8] = -6xi / length^2
    B[2, 12] = (3xi + 1) / length
    # Out-of-plane curvature from downward translation and chord-axis rotation.
    B[3, 3] = 6xi / length^2
    B[3, 5] = (3xi - 1) / length
    B[3, 9] = -6xi / length^2
    B[3, 11] = (3xi + 1) / length
    # Torsion from rotation about the span axis.
    B[4, 4] = -1 / length
    B[4, 10] = 1 / length
    return B
end

function beam_constitutive_matrix(stiffness, position, span_nodes, element)
    EA = spanwise_property_value(stiffness.EA, position, span_nodes, element)
    EI_in_plane = spanwise_property_value(
        stiffness.EI_in_plane,
        position,
        span_nodes,
        element,
    )
    EI_out_of_plane = spanwise_property_value(
        stiffness.EI_out_of_plane,
        position,
        span_nodes,
        element,
    )
    EI_coupling = spanwise_property_value(
        _physical_field(stiffness, :EI_coupling, 0.0),
        position,
        span_nodes,
        element,
    )
    GJ = spanwise_property_value(stiffness.GJ, position, span_nodes, element)
    return [EA 0.0 0.0 0.0;
            0.0 EI_in_plane EI_coupling 0.0;
            0.0 EI_coupling EI_out_of_plane 0.0;
            0.0 0.0 0.0 GJ]
end

"""Integrate the general 3D beam stiffness over one possibly nonuniform element."""
function general_beam_element_stiffness(stiffness, left, right, span_nodes, element)
    right > left || throw(ArgumentError("beam element length must be positive"))
    element_length = right - left
    element_center = (left + right) / 2
    fields = (:EA, :EI_in_plane, :EI_out_of_plane, :EI_coupling, :GJ)
    breaks = Float64[]
    for field in fields
        property = _physical_field(stiffness, field, 0.0)
        append!(breaks, _property_breakpoints(property, left, right, span_nodes))
    end
    edges = unique(sort([Float64(left); breaks; Float64(right)]))
    gauss_points = (-sqrt(3 / 5), 0.0, sqrt(3 / 5))
    gauss_weights = (5 / 9, 8 / 9, 5 / 9)
    matrix = zeros(12, 12)
    for interval in 1:(length(edges) - 1)
        interval_left = edges[interval]
        interval_right = edges[interval + 1]
        interval_length = interval_right - interval_left
        interval_center = (interval_left + interval_right) / 2
        for (point, weight) in zip(gauss_points, gauss_weights)
            position = interval_center + point * interval_length / 2
            xi = 2 * (position - element_center) / element_length
            B = beam_strain_displacement_matrix(element_length, xi)
            C = beam_constitutive_matrix(stiffness, position, span_nodes, element)
            matrix .+= B' * C * B * (weight * interval_length / 2)
        end
    end
    return 0.5 .* (matrix .+ matrix')
end

function _skew_matrix(vector)
    length(vector) == 3 || throw(DimensionMismatch("offset must have three components"))
    return [0.0 -vector[3] vector[2];
            vector[3] 0.0 -vector[1];
            -vector[2] vector[1] 0.0]
end

"""Return a six-by-six nodal rigid-body spatial inertia about the beam axis."""
function nodal_spatial_inertia(mass, inertia_at_cg, offset)
    mass >= 0 || throw(ArgumentError("nodal mass must be nonnegative"))
    size(inertia_at_cg) == (3, 3) || throw(DimensionMismatch(
        "nodal inertia tensor must be three-by-three",
    ))
    skew = _skew_matrix(offset)
    inertia_at_axis = Matrix(inertia_at_cg) .- mass .* (skew * skew)
    identity3 = Matrix{promote_type(Float64, eltype(inertia_at_cg))}(I, 3, 3)
    return [mass .* identity3 -mass .* skew;
            mass .* skew inertia_at_axis]
end

"""
    consistent_distributed_element_mass(length, mass_per_length,
        torsional_inertia_per_length, cg_offset_chord)

Return the full 12-by-12 consistent beam mass. The active Bohnisch subspace is
`[w_down, theta_chord, theta_span]` at each node. Its translational, torsional,
and bending-torsion coupling blocks reproduce the supplied reduced element
matrix exactly. Axial and in-plane consistent mass terms keep the full 3D beam
mass positive definite while their stiffness can be raised independently.
"""
function consistent_distributed_element_mass(
    length::Real,
    mass_per_length::Real,
    torsional_inertia_per_length::Real,
    cg_offset_chord::Real,
)
    length > 0 || throw(ArgumentError("beam element length must be positive"))
    mass_per_length > 0 || throw(ArgumentError("mass_per_length must be positive"))
    torsional_inertia_per_length > 0 || throw(ArgumentError(
        "torsional_inertia_per_length must be positive",
    ))

    matrix = zeros(12, 12)

    # Distributed axial translational mass.
    axial = mass_per_length * length / 6 .* [2.0 1.0; 1.0 2.0]
    matrix[[1, 7], [1, 7]] .+= axial

    # Distributed transverse mass with Euler-Bernoulli interpolation.
    bending = mass_per_length * length / 420 .* [
        156.0 22length 54.0 -13length;
        22length 4length^2 13length -3length^2;
        54.0 13length 156.0 -22length;
        -13length -3length^2 -22length 4length^2
    ]
    matrix[[2, 6, 8, 12], [2, 6, 8, 12]] .+= bending
    matrix[[3, 5, 9, 11], [3, 5, 9, 11]] .+= bending

    # Distributed sectional torsional inertia.
    torsion = torsional_inertia_per_length * length / 6 .* [2.0 1.0; 1.0 2.0]
    matrix[[4, 10], [4, 10]] .+= torsion

    # Bending-torsion inertial coupling caused by the signed offset between the
    # sectional center of mass and the beam reference/elastic axis.
    coupling = mass_per_length * length * cg_offset_chord .* [
        0.0 0.0 7 / 20 0.0 0.0 3 / 20;
        0.0 0.0 length / 20 0.0 0.0 length / 30;
        7 / 20 length / 20 0.0 3 / 20 -length / 30 0.0;
        0.0 0.0 3 / 20 0.0 0.0 7 / 20;
        0.0 0.0 -length / 30 0.0 0.0 -length / 20;
        3 / 20 length / 30 0.0 7 / 20 -length / 20 0.0
    ]
    active = [3, 5, 4, 9, 11, 10]
    matrix[active, active] .+= coupling
    return 0.5 .* (matrix .+ matrix')
end

function _assemble_lumped_wing_mass!(matrix, inertia, node_count)
    blocks = _physical_field(inertia, :spatial_inertia_blocks, nothing)
    if !isnothing(blocks)
        length(blocks) == node_count || throw(DimensionMismatch(
            "spatial_inertia_blocks must contain one block per wing node",
        ))
        for node in 1:node_count
            block = Matrix{Float64}(blocks[node])
            size(block) == (6, 6) || throw(DimensionMismatch(
                "each spatial inertia block must be six-by-six",
            ))
            matrix[wing_node_dofs(node), wing_node_dofs(node)] .+=
                0.5 .* (block .+ block')
        end
        return matrix
    end

    masses = inertia.nodal_mass
    offsets = inertia.cg_offset
    tensors = inertia.inertia_at_cg
    length(masses) == node_count || throw(DimensionMismatch(
        "nodal_mass must contain one value per wing node",
    ))
    length(offsets) == node_count || throw(DimensionMismatch(
        "cg_offset must contain one vector per wing node",
    ))
    length(tensors) == node_count || throw(DimensionMismatch(
        "inertia_at_cg must contain one tensor per wing node",
    ))
    for node in 1:node_count
        block = nodal_spatial_inertia(masses[node], tensors[node], offsets[node])
        matrix[wing_node_dofs(node), wing_node_dofs(node)] .+= block
    end
    return matrix
end

function _assemble_consistent_wing_mass!(matrix, inertia, span_nodes)
    for element in 1:(length(span_nodes) - 1)
        left, right = span_nodes[element], span_nodes[element + 1]
        center = (left + right) / 2
        mass_per_length = spanwise_property_value(
            inertia.mass_per_length,
            center,
            span_nodes,
            element,
        )
        torsional_inertia = spanwise_property_value(
            inertia.torsional_inertia_per_length,
            center,
            span_nodes,
            element,
        )
        cg_offset = spanwise_property_value(
            _physical_field(inertia, :cg_offset_chord, 0.0),
            center,
            span_nodes,
            element,
        )
        element_mass = consistent_distributed_element_mass(
            right - left,
            mass_per_length,
            torsional_inertia,
            cg_offset,
        )
        element_dofs = first(wing_node_dofs(element)):last(wing_node_dofs(element + 1))
        matrix[element_dofs, element_dofs] .+= element_mass
    end
    return matrix
end

"""Resolve one propeller's angular speed from its physical operating law."""
function propeller_angular_speed(propeller, operating_condition)
    model = propeller.speed_model
    if model == :fixed_omega
        omega = propeller.omega_rad_s
    elseif model == :fixed_rpm
        omega = propeller.rpm * 2pi / 60
    elseif model == :constant_advance_ratio
        advance_ratio = propeller.advance_ratio
        advance_ratio != 0 || throw(ArgumentError("advance_ratio must be nonzero"))
        propeller.radius > 0 || throw(ArgumentError("propeller radius must be positive"))
        omega = pi * operating_condition.freestream_speed /
            (advance_ratio * propeller.radius)
    else
        throw(ArgumentError("unknown propeller speed model: $model"))
    end
    isfinite(omega) || throw(ArgumentError("propeller angular speed must be finite"))
    return omega
end

function _propeller_damping_ratio(propeller, name)
    common = _physical_field(propeller, :damping_ratio, 0.0)
    return _physical_field(propeller, name, common)
end

function _constraint_dofs(wing, wing_dof_count)
    boundary = _physical_field(
        wing,
        :boundary_conditions,
        (; constrained_nodes = [1], constrained_local_dofs = collect(1:6)),
    )
    nodes = _physical_field(boundary, :constrained_nodes, [1])
    local_dofs = _physical_field(boundary, :constrained_local_dofs, collect(1:6))
    node_count = wing_dof_count ÷ WING_DOFS_PER_NODE
    all(1 .<= nodes .<= node_count) || throw(ArgumentError(
        "a constrained wing node lies outside the mesh",
    ))
    all(1 .<= local_dofs .<= WING_DOFS_PER_NODE) || throw(ArgumentError(
        "a constrained local wing DOF lies outside 1:6",
    ))
    return sort!(unique!(Int[
        WING_DOFS_PER_NODE * (node - 1) + local_dof
        for node in nodes for local_dof in local_dofs
    ]))
end

"""
    assemble_wing_propeller_structure(wing, propellers, operating_condition)

Assemble a common full-beam `M`, `C`, and `K` from physical inputs. `wing` uses
the general six-DOF beam stiffness and selects either `:lumped_nodal` or
`:consistent_distributed` inertia. Every propeller is attached directly to its
configured wing node and may have independent properties and speed law.
"""
function assemble_wing_propeller_structure(wing, propellers, operating_condition)
    span_nodes = collect(Float64, wing.span_nodes)
    issorted(span_nodes) || throw(ArgumentError("wing span_nodes must be sorted"))
    length(span_nodes) >= 2 || throw(ArgumentError("wing requires at least two nodes"))
    any(diff(span_nodes) .<= 0) && throw(ArgumentError("wing elements must have positive length"))
    node_count = length(span_nodes)
    element_count = node_count - 1
    wing_dof_count = WING_DOFS_PER_NODE * node_count

    wing_stiffness = zeros(wing_dof_count, wing_dof_count)
    for element in 1:element_count
        element_matrix = general_beam_element_stiffness(
            wing.stiffness,
            span_nodes[element],
            span_nodes[element + 1],
            span_nodes,
            element,
        )
        element_dofs = first(wing_node_dofs(element)):last(wing_node_dofs(element + 1))
        wing_stiffness[element_dofs, element_dofs] .+= element_matrix
    end

    wing_mass = zeros(wing_dof_count, wing_dof_count)
    mass_model = wing.mass_model
    if mass_model in (:lumped, :lumped_nodal)
        _assemble_lumped_wing_mass!(wing_mass, wing.inertia, node_count)
        mass_model = :lumped_nodal
    elseif mass_model == :consistent_distributed
        _assemble_consistent_wing_mass!(wing_mass, wing.inertia, span_nodes)
    else
        throw(ArgumentError("unknown wing mass model: $(wing.mass_model)"))
    end

    wing_damping = zeros(wing_dof_count, wing_dof_count)
    damping = _physical_field(wing, :damping, (;))
    mass_coefficient = _physical_field(damping, :mass_coefficient, 0.0)
    stiffness_coefficient = _physical_field(damping, :stiffness_coefficient, 0.0)
    wing_damping .+= mass_coefficient .* wing_mass .+
        stiffness_coefficient .* wing_stiffness

    propeller_count = length(propellers)
    propeller_dof_count = 2 * propeller_count
    propeller_mass = zeros(propeller_dof_count, propeller_dof_count)
    propeller_damping = zeros(propeller_dof_count, propeller_dof_count)
    propeller_stiffness = zeros(propeller_dof_count, propeller_dof_count)
    propeller_to_wing_mass = zeros(propeller_dof_count, wing_dof_count)
    propeller_to_wing_damping = zeros(propeller_dof_count, wing_dof_count)
    wing_to_propeller_mass = zeros(wing_dof_count, propeller_dof_count)
    attached_wing_mass = zeros(wing_dof_count, wing_dof_count)
    wing_to_propeller_damping = zeros(wing_dof_count, propeller_dof_count)
    attachment_operators = Vector{Matrix{Float64}}(undef, propeller_count)
    attachment_nodes = Vector{Int}(undef, propeller_count)
    angular_speeds = Vector{Float64}(undef, propeller_count)

    for (index, propeller) in enumerate(propellers)
        node = Int(propeller.attachment_node)
        1 <= node <= node_count || throw(ArgumentError(
            "propeller $index attachment node $node lies outside the wing mesh",
        ))
        attachment_nodes[index] = node
        node_dofs = wing_node_dofs(node)
        propeller_dofs = (2index - 1):(2index)

        # Concentrated propeller/nacelle inertia is applied directly to the
        # selected structural wing node; no neighboring-node weights exist.
        attachment = zeros(6, wing_dof_count)
        attachment[:, node_dofs] .= Matrix{Float64}(I, 6, 6)
        attachment_operators[index] = attachment

        pitch_inertia = propeller.pitch_inertia
        yaw_inertia = propeller.yaw_inertia
        pitch_stiffness = propeller.pitch_stiffness
        yaw_stiffness = propeller.yaw_stiffness
        pitch_inertia > 0 && yaw_inertia > 0 || throw(ArgumentError(
            "propeller pitch/yaw inertias must be positive",
        ))
        pitch_stiffness > 0 && yaw_stiffness > 0 || throw(ArgumentError(
            "propeller pitch/yaw stiffnesses must be positive",
        ))
        omega = propeller_angular_speed(propeller, operating_condition)
        angular_speeds[index] = omega
        pitch_ratio = _propeller_damping_ratio(propeller, :pitch_damping_ratio)
        yaw_ratio = _propeller_damping_ratio(propeller, :yaw_damping_ratio)

        local_propeller_mass = [pitch_inertia 0.0; 0.0 yaw_inertia]
        local_propeller_stiffness = [pitch_stiffness 0.0; 0.0 yaw_stiffness]
        local_propeller_damping = [
            2pitch_ratio * sqrt(pitch_inertia * pitch_stiffness) propeller.spin_inertia * omega;
            -propeller.spin_inertia * omega 2yaw_ratio * sqrt(yaw_inertia * yaw_stiffness)
        ]
        propeller_mass[propeller_dofs, propeller_dofs] .= local_propeller_mass
        propeller_damping[propeller_dofs, propeller_dofs] .= local_propeller_damping
        propeller_stiffness[propeller_dofs, propeller_dofs] .= local_propeller_stiffness

        pitch_first_moment = _physical_field(propeller, :pitch_first_moment, 0.0)
        yaw_first_moment = _physical_field(propeller, :yaw_first_moment, 0.0)
        pitch_cross_inertia = _physical_field(propeller, :pitch_cross_inertia, 0.0)
        yaw_cross_inertia = _physical_field(propeller, :yaw_cross_inertia, 0.0)
        wing_pitch_first_moment = _physical_field(
            propeller,
            :wing_pitch_first_moment,
            pitch_first_moment,
        )
        wing_yaw_first_moment = _physical_field(
            propeller,
            :wing_yaw_first_moment,
            yaw_first_moment,
        )
        wing_pitch_inertia = _physical_field(propeller, :wing_pitch_inertia, 0.0)
        wing_yaw_inertia = _physical_field(propeller, :wing_yaw_inertia, 0.0)

        local_cross_mass = [
            0.0 0.0 pitch_first_moment pitch_cross_inertia 0.0 0.0;
            -yaw_first_moment 0.0 0.0 0.0 0.0 yaw_cross_inertia
        ]
        propeller_to_wing_mass[propeller_dofs, :] .+= local_cross_mass * attachment
        wing_to_propeller_mass[:, propeller_dofs] .+= attachment' * local_cross_mass'

        local_attached_mass = zeros(6, 6)
        local_attached_mass[1, 1] = propeller.mass
        local_attached_mass[2, 2] = propeller.mass
        local_attached_mass[3, 3] = propeller.mass
        local_attached_mass[3, 4] = wing_pitch_first_moment
        local_attached_mass[4, 3] = wing_pitch_first_moment
        local_attached_mass[4, 4] = wing_pitch_inertia
        local_attached_mass[1, 6] = -wing_yaw_first_moment
        local_attached_mass[6, 1] = -wing_yaw_first_moment
        local_attached_mass[6, 6] = wing_yaw_inertia
        attached_wing_mass .+= attachment' * local_attached_mass * attachment

        # Propeller gyroscopic pitch-yaw coupling proportional to Ix*Omega.
        local_wing_to_propeller_gyro = [
            0.0 0.0;
            0.0 0.0;
            0.0 0.0;
            0.0 propeller.spin_inertia * omega;
            0.0 0.0;
            -propeller.spin_inertia * omega 0.0
        ]
        wing_to_propeller_damping[:, propeller_dofs] .+=
            attachment' * local_wing_to_propeller_gyro
        propeller_to_wing_damping[propeller_dofs, :] .+=
            -local_wing_to_propeller_gyro' * attachment
        local_wing_gyro = zeros(6, 6)
        local_wing_gyro[4, 6] = propeller.spin_inertia * omega
        local_wing_gyro[6, 4] = -propeller.spin_inertia * omega
        wing_damping .+= attachment' * local_wing_gyro * attachment
    end

    total_dof_count = wing_dof_count + propeller_dof_count
    M_global = [wing_mass + attached_wing_mass wing_to_propeller_mass;
                propeller_to_wing_mass propeller_mass]
    C_global = [wing_damping wing_to_propeller_damping;
                propeller_to_wing_damping propeller_damping]
    K_global = [wing_stiffness zeros(wing_dof_count, propeller_dof_count);
                zeros(propeller_dof_count, wing_dof_count) propeller_stiffness]
    # Optional system-level Rayleigh terms preserve configurations whose
    # validated damping definition applies to both wing and mount stiffness.
    C_global .+= _physical_field(damping, :global_mass_coefficient, 0.0) .* M_global
    C_global .+= _physical_field(damping, :global_stiffness_coefficient, 0.0) .* K_global

    constrained_dofs = _constraint_dofs(wing, wing_dof_count)
    free_dofs = setdiff(collect(1:total_dof_count), constrained_dofs)
    M = M_global[free_dofs, free_dofs]
    C = C_global[free_dofs, free_dofs]
    K = K_global[free_dofs, free_dofs]
    wing_free_dofs = count(<=(wing_dof_count), free_dofs)

    return (;
        M,
        C,
        K,
        M_global,
        C_global,
        K_global,
        M_wing = wing_mass,
        C_wing = wing_damping,
        K_wing = wing_stiffness,
        M_attached = attached_wing_mass,
        M_propeller = propeller_mass,
        C_propeller = propeller_damping,
        K_propeller = propeller_stiffness,
        propeller_to_wing_mass,
        wing_to_propeller_mass,
        propeller_to_wing_damping,
        wing_to_propeller_damping,
        attachment_operators,
        attachment_nodes,
        angular_speeds,
        constrained_dofs,
        free_dofs,
        ndof_free = length(free_dofs),
        ndof_wing_free = wing_free_dofs,
        ndof_propeller_free = propeller_dof_count,
        dofs_per_node = WING_DOFS_PER_NODE,
        node_count,
        element_count,
        mass_model,
    )
end

"""Return the lowest positive undamped natural frequencies in hertz."""
function structural_natural_frequencies(structural; count::Integer = 10)
    count > 0 || throw(ArgumentError("count must be positive"))
    solution = eigen(Symmetric(structural.K), Symmetric(structural.M))
    values = real.(solution.values)
    positive = sort(values[isfinite.(values) .& (values .> 0.0)])
    return sqrt.(positive[1:min(count, length(positive))]) ./ (2pi)
end
