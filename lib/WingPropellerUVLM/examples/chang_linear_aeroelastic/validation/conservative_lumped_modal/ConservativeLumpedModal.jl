module ConservativeLumpedModal

using LinearAlgebra

export SOURCE_NODES_M,
       REFERENCE_WING_MODES,
       assemble_wing_model,
       assemble_coupled_model,
       modal_analysis,
       remap_audit,
       remapped_nodal_properties,
       source_spatial_inertia_blocks

const FT_TO_M = 0.3048
const LBF_TO_N = 4.4482216152605
const SLUG_TO_KG = 14.5939029372064
const SLUG_FT2_TO_KGM2 = SLUG_TO_KG * FT_TO_M^2
const LBF_FT2_TO_NM2 = LBF_TO_N * FT_TO_M^2

# The spreadsheet stations close at 7.5 m to its printed precision.  Scaling
# by the printed end station removes only round-off and makes every remap cover
# exactly the same [0, span] interval.
const SPAN_M = 7.5
const SOURCE_NODES_FT = [
    0.0, 0.54675, 1.0936, 2.18722, 3.28084, 4.51115, 5.74147,
    6.97178, 8.2021, 10.2526, 12.30315, 14.353675, 16.4042,
    18.4547, 20.50525, 22.55577, 24.6063,
]
const SOURCE_NODES_M = (SOURCE_NODES_FT .* FT_TO_M) .* (
    SPAN_M / (SOURCE_NODES_FT[end] * FT_TO_M)
)

# Each row belongs to the flexible element ending at SOURCE_NODES_M[i + 1].
# The reference root node has no independent inertia record.
const MASS_SLUG = [
    0.4276211, 0.4139885, 0.736678, 0.4139886, 0.9436633,
    0.7336013, 0.3851399, 0.9932094, 0.9932177, 0.9932237,
    0.9931985, 0.971823, 0.9687671, 0.9109462, 0.9367904,
    0.3350064,
]
const IXX_SLUG_FT2 = [
    0.6055, 0.6281, 1.92, 0.6281, 2.3536, 1.909, 0.5645, 2.4813,
    2.4813, 2.4814, 2.4813, 2.3998, 2.3877, 2.2748, 2.3653, 0.459,
]
const IYY_SLUG_FT2 = [
    0.0969, 0.1028, 0.1269, 0.1028, 0.2878, 0.1266, 0.0967,
    0.3565, 0.3656, 0.3656, 0.3656, 0.3577, 0.3561, 0.3028,
    0.3572, 0.0797,
]
const IZZ_SLUG_FT2 = [
    0.6475, 0.6777, 1.9703, 0.6777, 2.3194, 1.9593, 0.6118,
    2.388, 2.388, 2.388, 2.3879, 2.3118, 2.3009, 2.2058,
    2.2509, 0.4844,
]
const IXY_SLUG_FT2 = [
    0.000483499, 0.000545311, -0.003364902, 0.000546007,
    -0.003767303, -0.003368369, 0.000500823, -0.004768277,
    -0.004781691, -0.004764693, -0.00476285, -0.004980494,
    -0.004983901, -0.005318305, -0.005590199, 0.000500334,
]
const IXZ_SLUG_FT2 = [
    0.001076907, -1.83141e-7, 0.001721633, -4.2213e-7,
    0.004030091, -0.001715619, 0.000169958, 0.002555006,
    0.002553573, 0.002555579, 0.002557012, 0.002683975,
    0.00262101, 0.003720733, 0.004553367, -0.001762152,
]
const IYZ_SLUG_FT2 = [
    -0.000821912, 6.13447e-6, 0.026378458, 7.88634e-6,
    0.061766804, -0.026536753, -7.79545e-6, 0.031499987,
    0.031497299, 0.03149642, 0.031500523, 0.030281192,
    0.03089678, 0.024007273, 0.035641812, 0.001067331,
]

# Literal workbook coordinates: X spanwise, Y chordwise, Z vertical.
const CG_X_FT = [
    0.027173898, 0.037438118, 0.021032451, 0.037438576,
    0.02464125, 0.021088361, 0.040243862, 0.031203696,
    0.03121296, 0.03120211, 0.031192549, 0.031900145,
    0.031983194, 0.034021895, 0.03582014, 0.039275863,
]
const CG_Y_FT = [
    -0.66071232, -0.668605558, -0.919885764, -0.668609685,
    -0.864794967, -0.921385197, -0.671753688, -0.857220095,
    -0.857206007, -0.857200098, -0.857224387, -0.863614669,
    -0.864554528, -0.874965167, -0.869545335, -0.664042589,
]
const CG_Z_FT = [
    -0.245439677, -0.17398455, 0.016601162, -0.014913878,
    0.121183138, -0.199929292, 0.137702621, -0.081678808,
    -0.059307771, -0.03693769, -0.014573273, 0.01200473,
    0.03224826, 0.09014999, 0.048423281, 0.248940617,
]

const ETA_STIFFNESS = [
    0.04443, 0.04444, 0.08888, 0.08889, 0.13332, 0.13333,
    0.18332, 0.18333, 0.23332, 0.23333, 0.28332, 0.28333,
    0.33333, 0.33334, 0.41666, 0.41667, 0.49999, 0.5,
    0.58332, 0.58333, 0.66666, 0.66667, 0.74999, 0.75,
    0.83333, 0.83334, 0.91667, 1.0,
]
const EIYY_LBF_FT2 = [
    32836890.0, 32836890.0, 32836890.0, 32836890.0, 32836890.0,
    20120070.0, 20120070.0, 20120070.0, 20120070.0,
    16181100.0, 16181100.0, 16181100.0, 16181100.0,
    19251360.0, 19251360.0, 19251360.0, 19251360.0,
    18286380.0, 18286380.0, 18286380.0, 18286380.0,
    17183860.0, 17183860.0, 17183860.0, 17183860.0,
    11198940.0, 11198940.0, 11198940.0,
]
const EIZZ_LBF_FT2 = [
    325381340.0, 325381340.0, 325381340.0, 325381340.0, 325381340.0,
    208531310.0, 208531310.0, 208531310.0, 208531310.0,
    163131180.0, 163131180.0, 163131180.0, 163131180.0,
    194913580.0, 194913580.0, 194913580.0, 194913580.0,
    181528760.0, 181528760.0, 181528760.0, 181528760.0,
    163800810.0, 163800810.0, 163800810.0, 163800810.0,
    128273870.0, 128273870.0, 128273870.0,
]
const EIZY_LBF_FT2 = [
    4.7103163e6, 4.7103163e6, 4.7103163e6, 4.7103163e6,
    4.7103163e6, 2.9685812e6, 2.9685812e6, 2.9685812e6,
    2.9685812e6, 2.5060481e6, 2.5060481e6, 2.5060481e6,
    2.5060481e6, 3.3875434e6, 3.3875434e6, 3.3875434e6,
    3.3875434e6, 3.1960560e6, 3.1960560e6, 3.1960560e6,
    3.1960560e6, 3.3050893e6, 3.3050893e6, 3.3050893e6,
    3.3050893e6, 2.4234529e6, 2.4234529e6, 2.4234529e6,
]
const GJ_LBF_FT2 = [
    3.05e7, 3.05e7, 3.05e7, 3.05e7, 3.05e7,
    1.84e7, 1.84e7, 1.84e7, 1.84e7,
    1.40e7, 1.40e7, 1.40e7, 1.40e7,
    1.62e7, 1.62e7, 1.62e7, 1.62e7,
    1.52e7, 1.52e7, 1.52e7, 1.52e7,
    1.30e7, 1.30e7, 1.30e7, 1.30e7,
    9.82e6, 9.82e6, 9.82e6,
]
const EA_LBF = [
    168148740.0, 168148740.0, 168148740.0, 168148740.0, 168148740.0,
    111286380.0, 111286380.0, 111286380.0, 111286380.0,
    85691907.0, 85691907.0, 85691907.0, 85691907.0,
    102776990.0, 102776990.0, 102776990.0, 102776990.0,
    96393539.0, 96393539.0, 96393539.0, 96393539.0,
    88126997.0, 88126997.0, 88126697.0, 88126997.0,
    100925580.0, 100925580.0, 100925580.0,
]

const REFERENCE_WING_MODES = [
    (label = "OOP1", family = :out_of_plane, family_order = 1, frequency_hz = 6.63),
    (label = "IP1", family = :in_plane, family_order = 1, frequency_hz = 20.93),
    (label = "OOP2", family = :out_of_plane, family_order = 2, frequency_hz = 36.66),
    (label = "T1", family = :torsion, family_order = 1, frequency_hz = 43.50),
    (label = "OOP3", family = :out_of_plane, family_order = 3, frequency_hz = 93.31),
]

skew3(v) = [0.0 -v[3] v[2]; v[3] 0.0 -v[1]; -v[2] v[1] 0.0]

function spatial_inertia_block(mass, inertia_at_cg, offset)
    S = skew3(offset)
    I3 = Matrix{Float64}(I, 3, 3)
    return [mass .* I3 -mass .* S; mass .* S inertia_at_cg .- mass .* S * S]
end

"""Return the 16 spreadsheet spatial-inertia blocks about the beam axis."""
function source_spatial_inertia_blocks()
    blocks = Matrix{Float64}[]
    for i in eachindex(MASS_SLUG)
        mass = MASS_SLUG[i] * SLUG_TO_KG
        Jcg = SLUG_FT2_TO_KGM2 .* [
            IXX_SLUG_FT2[i] IXY_SLUG_FT2[i] IXZ_SLUG_FT2[i];
            IXY_SLUG_FT2[i] IYY_SLUG_FT2[i] IYZ_SLUG_FT2[i];
            IXZ_SLUG_FT2[i] IYZ_SLUG_FT2[i] IZZ_SLUG_FT2[i]
        ]
        offset = FT_TO_M .* [CG_X_FT[i], CG_Y_FT[i], CG_Z_FT[i]]
        push!(blocks, spatial_inertia_block(mass, Jcg, offset))
    end
    return blocks
end

"""
    integrated_linear_density_remap(target_nodes)

Interpret the spreadsheet blocks as nodal reference values. Divide every
nodal block by its source tributary length, linearly interpolate the complete
6-by-6 spatial-inertia density between reference nodes, and integrate that
density exactly over every target element. The integrated target property is
lumped at the element's outer node.

The source tributary length equals the integral of its piecewise-linear hat
function. Consequently, this interpolation exactly preserves the sum of all
source blocks without a posteriori component scaling. All interpolation and
integration weights are nonnegative, so symmetry and positive semidefiniteness
are preserved as well. The fixed root node intentionally remains massless;
all free nodes have a physical lumped inertia.
"""
function integrated_linear_density_remap(target_nodes::AbstractVector{<:Real})
    issorted(target_nodes) || throw(ArgumentError("target_nodes must be sorted"))
    isapprox(target_nodes[1], 0.0; atol = 1e-12) ||
        throw(ArgumentError("target mesh must begin at zero"))
    isapprox(target_nodes[end], SPAN_M; atol = 1e-10) ||
        throw(ArgumentError("target mesh must end at SPAN_M"))

    source_blocks = vcat([zeros(6, 6)], source_spatial_inertia_blocks())
    source_lengths = nodal_tributary_lengths(SOURCE_NODES_M)
    source_density = [
        source_blocks[i] ./ source_lengths[i] for i in eachindex(source_blocks)
    ]
    target_blocks = [zeros(6, 6) for _ in eachindex(target_nodes)]
    for target_element in 1:(length(target_nodes) - 1)
        left = target_nodes[target_element]
        right = target_nodes[target_element + 1]
        block = target_blocks[target_element + 1]
        for source_element in 1:(length(SOURCE_NODES_M) - 1)
            source_left = SOURCE_NODES_M[source_element]
            source_right = SOURCE_NODES_M[source_element + 1]
            overlap_left = max(left, source_left)
            overlap_right = min(right, source_right)
            overlap_right <= overlap_left && continue

            interval = source_right - source_left
            t_left = (overlap_left - source_left) / interval
            t_right = (overlap_right - source_left) / interval
            # Exact integrals of (1-t) and t over the overlap interval.
            left_weight = interval * (
                (t_right - t_right^2 / 2) - (t_left - t_left^2 / 2)
            )
            right_weight = interval * (t_right^2 - t_left^2) / 2
            block .+= left_weight .* source_density[source_element]
            block .+= right_weight .* source_density[source_element + 1]
        end
        block .= 0.5 .* (block .+ block')
    end
    return target_blocks
end

"""
    control_volume_spatial_remap(target_nodes)

Integrate the same linearly interpolated spatial-inertia density over nodal
control volumes bounded by adjacent element midpoints.  This locates each lump
at the center of its tributary span and removes the half-element tipward bias
of outer-node element lumping.  The root control-volume block is transferred
to the first free node because the reference root record is zero and the root
DOFs are eliminated.  The full free-node spatial-inertia total is therefore
conserved with no massless free node.
"""
function control_volume_spatial_remap(target_nodes::AbstractVector{<:Real})
    issorted(target_nodes) || throw(ArgumentError("target_nodes must be sorted"))
    isapprox(target_nodes[1], 0.0; atol = 1e-12) ||
        throw(ArgumentError("target mesh must begin at zero"))
    isapprox(target_nodes[end], SPAN_M; atol = 1e-10) ||
        throw(ArgumentError("target mesh must end at SPAN_M"))

    source_blocks = vcat([zeros(6, 6)], source_spatial_inertia_blocks())
    source_lengths = nodal_tributary_lengths(SOURCE_NODES_M)
    source_density = [
        source_blocks[i] ./ source_lengths[i] for i in eachindex(source_blocks)
    ]
    control_edges = [
        target_nodes[1];
        (target_nodes[1:end-1] .+ target_nodes[2:end]) ./ 2;
        target_nodes[end]
    ]
    target_blocks = [zeros(6, 6) for _ in eachindex(target_nodes)]

    for target in eachindex(target_nodes)
        left = control_edges[target]
        right = control_edges[target + 1]
        block = target_blocks[target]
        for source_element in 1:(length(SOURCE_NODES_M) - 1)
            source_left = SOURCE_NODES_M[source_element]
            source_right = SOURCE_NODES_M[source_element + 1]
            overlap_left = max(left, source_left)
            overlap_right = min(right, source_right)
            overlap_right <= overlap_left && continue
            interval = source_right - source_left
            t_left = (overlap_left - source_left) / interval
            t_right = (overlap_right - source_left) / interval
            left_weight = interval * (
                (t_right - t_right^2 / 2) - (t_left - t_left^2 / 2)
            )
            right_weight = interval * (t_right^2 - t_left^2) / 2
            block .+= left_weight .* source_density[source_element]
            block .+= right_weight .* source_density[source_element + 1]
        end
        block .= 0.5 .* (block .+ block')
    end

    target_blocks[2] .+= target_blocks[1]
    target_blocks[1] .= 0.0
    return target_blocks
end

"""
    sampled_basis_normalized_remap(target_nodes)

Diagnostic conservative variant of the user's former nodal interpolation. First form the
same linearly interpolated spatial-inertia density used by the old code and
sample it at every target node.  The target tributary length supplies the
lumped quadrature weight.  Instead of globally rescaling each signed tensor
component, normalize the contribution of each nonnegative source hat function
before combining the complete source spatial-inertia blocks.

It reproduces the source blocks when the source mesh is requested, but it is
not suitable for a coarse target mesh: a narrow source hat may contain no
target sampling point. The production-quality remap below instead integrates
the interpolated density exactly.
"""
function sampled_basis_normalized_remap(target_nodes::AbstractVector{<:Real})
    issorted(target_nodes) || throw(ArgumentError("target_nodes must be sorted"))
    isapprox(target_nodes[1], 0.0; atol = 1e-12) ||
        throw(ArgumentError("target mesh must begin at zero"))
    isapprox(target_nodes[end], SPAN_M; atol = 1e-10) ||
        throw(ArgumentError("target mesh must end at SPAN_M"))

    source_blocks = vcat([zeros(6, 6)], source_spatial_inertia_blocks())
    source_lengths = nodal_tributary_lengths(SOURCE_NODES_M)
    target_lengths = nodal_tributary_lengths(target_nodes)
    number_of_targets = length(target_nodes)
    number_of_sources = length(SOURCE_NODES_M)
    transfer = zeros(number_of_targets, number_of_sources)

    for target in eachindex(target_nodes)
        x = target_nodes[target]
        if x <= SOURCE_NODES_M[1]
            transfer[target, 1] = target_lengths[target] / source_lengths[1]
        elseif x >= SOURCE_NODES_M[end]
            transfer[target, end] = target_lengths[target] / source_lengths[end]
        else
            right = searchsortedfirst(SOURCE_NODES_M, x)
            left = right - 1
            alpha = (x - SOURCE_NODES_M[left]) /
                (SOURCE_NODES_M[right] - SOURCE_NODES_M[left])
            transfer[target, left] =
                target_lengths[target] * (1 - alpha) / source_lengths[left]
            transfer[target, right] =
                target_lengths[target] * alpha / source_lengths[right]
        end
    end

    # Normalize each source basis separately.  Every complete source block is
    # therefore conserved without dividing by a small signed component total.
    for source in 2:number_of_sources
        contribution_sum = sum(transfer[:, source])
        contribution_sum > eps(Float64) || error(
            "Target mesh is too coarse to sample source inertia station $source",
        )
        transfer[:, source] ./= contribution_sum
    end

    target_blocks = Matrix{Float64}[]
    for target in 1:number_of_targets
        block = zeros(6, 6)
        for source in 2:number_of_sources
            block .+= transfer[target, source] .* source_blocks[source]
        end
        push!(target_blocks, 0.5 .* (block .+ block'))
    end
    return target_blocks
end

# Public name used by the retained validation scripts.
# Nodal control volumes match a lumped nodal model without the outer-node bias.
conservative_spatial_remap(target_nodes) =
    control_volume_spatial_remap(target_nodes)

"""
    element_overlap_spatial_remap(target_nodes)

Alternative interpretation retained for sensitivity analysis. It assumes each
of the 16 inertia records is the integrated property of the preceding
spanwise flexible beam element, constructs a piecewise-constant density in
that element, and transfers it by exact interval overlap. The workbook
connectivity shows that the records actually belong to separate rigid offset
elements, so this is not the preferred reference interpretation.
"""
function element_overlap_spatial_remap(target_nodes)
    source_blocks = source_spatial_inertia_blocks()
    target_blocks = [zeros(6, 6) for _ in eachindex(target_nodes)]
    for target_element in 1:(length(target_nodes) - 1)
        left = target_nodes[target_element]
        right = target_nodes[target_element + 1]
        for source_element in eachindex(source_blocks)
            source_left = SOURCE_NODES_M[source_element]
            source_right = SOURCE_NODES_M[source_element + 1]
            overlap = max(0.0, min(right, source_right) - max(left, source_left))
            overlap == 0.0 && continue
            target_blocks[target_element + 1] .+=
                overlap / (source_right - source_left) .* source_blocks[source_element]
        end
        target_blocks[target_element + 1] .= 0.5 .* (
            target_blocks[target_element + 1] .+
            target_blocks[target_element + 1]'
        )
    end
    return target_blocks
end

function properties_from_spatial_block(block)
    mass = tr(block[1:3, 1:3]) / 3
    mass <= 100eps(Float64) && return (
        mass = 0.0,
        cg = zeros(3),
        inertia_at_cg = zeros(3, 3),
    )
    S = -block[1:3, 4:6] ./ mass
    S = 0.5 .* (S .- S')
    cg = [S[3, 2], S[1, 3], S[2, 1]]
    inertia_at_cg = block[4:6, 4:6] .+ mass .* S * S
    inertia_at_cg = 0.5 .* (inertia_at_cg .+ inertia_at_cg')
    return (; mass, cg, inertia_at_cg)
end

"""
    remapped_nodal_properties(Ne; mesh=:uniform)

Return the conventional mass, CG-offset, and center-of-mass-inertia arrays for
an arbitrary lumped target mesh. These arrays have the same physical meaning
as `m_node_vec`, `cg_*_node_vec`, and `I**_node_vec` in the active Chang case.
"""
function remapped_nodal_properties(number_of_elements; mesh = :uniform)
    nodes = target_nodes(number_of_elements; mesh)
    blocks = conservative_spatial_remap(nodes)
    properties = properties_from_spatial_block.(blocks)
    return (;
        nodes,
        blocks,
        mass = [p.mass for p in properties],
        cg_x = [p.cg[1] for p in properties],
        cg_y = [p.cg[2] for p in properties],
        cg_z = [p.cg[3] for p in properties],
        Ixx = [p.inertia_at_cg[1, 1] for p in properties],
        Iyy = [p.inertia_at_cg[2, 2] for p in properties],
        Izz = [p.inertia_at_cg[3, 3] for p in properties],
        Ixy = [p.inertia_at_cg[1, 2] for p in properties],
        Ixz = [p.inertia_at_cg[1, 3] for p in properties],
        Iyz = [p.inertia_at_cg[2, 3] for p in properties],
    )
end

function linear_interpolate(x, y, xq)
    xq <= x[1] && return y[1]
    xq >= x[end] && return y[end]
    right = searchsortedfirst(x, xq)
    left = right - 1
    alpha = (xq - x[left]) / (x[right] - x[left])
    return (1 - alpha) * y[left] + alpha * y[right]
end

function interpolate_vector(x, y, xq)
    return [linear_interpolate(x, y, value) for value in xq]
end

function nodal_tributary_lengths(nodes)
    lengths = zeros(length(nodes))
    element_lengths = diff(nodes)
    lengths[1] = element_lengths[1] / 2
    lengths[end] = element_lengths[end] / 2
    for i in 2:(length(nodes) - 1)
        lengths[i] = (element_lengths[i - 1] + element_lengths[i]) / 2
    end
    return lengths
end

function source_nodal_components()
    mass = [0.0; MASS_SLUG .* SLUG_TO_KG]
    cgx = [CG_X_FT[1]; CG_X_FT] .* FT_TO_M
    cgy = [CG_Y_FT[1]; CG_Y_FT] .* FT_TO_M
    cgz = [CG_Z_FT[1]; CG_Z_FT] .* FT_TO_M
    ixx = [0.0; IXX_SLUG_FT2 .* SLUG_FT2_TO_KGM2]
    iyy = [0.0; IYY_SLUG_FT2 .* SLUG_FT2_TO_KGM2]
    izz = [0.0; IZZ_SLUG_FT2 .* SLUG_FT2_TO_KGM2]
    ixy = [0.0; IXY_SLUG_FT2 .* SLUG_FT2_TO_KGM2]
    ixz = [0.0; IXZ_SLUG_FT2 .* SLUG_FT2_TO_KGM2]
    iyz = [0.0; IYZ_SLUG_FT2 .* SLUG_FT2_TO_KGM2]
    return (; mass, cgx, cgy, cgz, ixx, iyy, izz, ixy, ixz, iyz)
end

function scalar_density_remap(values, target_nodes; signed_rescale)
    source_lengths = nodal_tributary_lengths(SOURCE_NODES_M)
    target_lengths = nodal_tributary_lengths(target_nodes)
    density = values ./ source_lengths
    remapped = interpolate_vector(SOURCE_NODES_M, density, target_nodes) .* target_lengths
    if signed_rescale
        total = sum(remapped)
        tolerance = 100eps(Float64) * max(sum(abs, remapped), abs(sum(values)), 1.0)
        abs(total) > tolerance || error("near-zero remapped signed total")
        remapped .*= sum(values) / total
    end
    return remapped
end

function componentwise_remap(target_nodes; signed_rescale)
    p = source_nodal_components()
    m = scalar_density_remap(p.mass, target_nodes; signed_rescale)
    mx = scalar_density_remap(p.mass .* p.cgx, target_nodes; signed_rescale)
    my = scalar_density_remap(p.mass .* p.cgy, target_nodes; signed_rescale)
    mz = scalar_density_remap(p.mass .* p.cgz, target_nodes; signed_rescale)
    fallback_x = interpolate_vector(SOURCE_NODES_M, p.cgx, target_nodes)
    fallback_y = interpolate_vector(SOURCE_NODES_M, p.cgy, target_nodes)
    fallback_z = interpolate_vector(SOURCE_NODES_M, p.cgz, target_nodes)
    cgx = [m[i] > eps() ? mx[i] / m[i] : fallback_x[i] for i in eachindex(m)]
    cgy = [m[i] > eps() ? my[i] / m[i] : fallback_y[i] for i in eachindex(m)]
    cgz = [m[i] > eps() ? mz[i] / m[i] : fallback_z[i] for i in eachindex(m)]

    ixx_axis = p.ixx .+ p.mass .* (p.cgy .^ 2 .+ p.cgz .^ 2)
    iyy_axis = p.iyy .+ p.mass .* (p.cgx .^ 2 .+ p.cgz .^ 2)
    izz_axis = p.izz .+ p.mass .* (p.cgx .^ 2 .+ p.cgy .^ 2)
    ixy_axis = p.ixy .- p.mass .* p.cgx .* p.cgy
    ixz_axis = p.ixz .- p.mass .* p.cgx .* p.cgz
    iyz_axis = p.iyz .- p.mass .* p.cgy .* p.cgz
    axis_components = map(
        values -> scalar_density_remap(values, target_nodes; signed_rescale),
        (ixx_axis, iyy_axis, izz_axis, ixy_axis, ixz_axis, iyz_axis),
    )
    ixx_a, iyy_a, izz_a, ixy_a, ixz_a, iyz_a = axis_components

    blocks = Matrix{Float64}[]
    minimum_cg_inertia_eigenvalue = Inf
    maximum_offset = 0.0
    for i in eachindex(m)
        Jcg = [
            ixx_a[i] - m[i] * (cgy[i]^2 + cgz[i]^2)  ixy_a[i] + m[i] * cgx[i] * cgy[i]  ixz_a[i] + m[i] * cgx[i] * cgz[i];
            ixy_a[i] + m[i] * cgx[i] * cgy[i]  iyy_a[i] - m[i] * (cgx[i]^2 + cgz[i]^2)  iyz_a[i] + m[i] * cgy[i] * cgz[i];
            ixz_a[i] + m[i] * cgx[i] * cgz[i]  iyz_a[i] + m[i] * cgy[i] * cgz[i]  izz_a[i] - m[i] * (cgx[i]^2 + cgy[i]^2)
        ]
        Jcg = 0.5 .* (Jcg .+ Jcg')
        minimum_cg_inertia_eigenvalue = min(
            minimum_cg_inertia_eigenvalue,
            minimum(eigvals(Symmetric(Jcg))),
        )
        maximum_offset = max(maximum_offset, norm([cgx[i], cgy[i], cgz[i]]))
        push!(blocks, spatial_inertia_block(m[i], Jcg, [cgx[i], cgy[i], cgz[i]]))
    end
    return (; blocks, minimum_cg_inertia_eigenvalue, maximum_offset)
end

function strain_displacement_matrix(L, xi)
    B = zeros(4, 12)
    B[1, 1] = -1 / L
    B[1, 7] = 1 / L
    B[2, 2] = 6xi / L^2
    B[2, 6] = (3xi - 1) / L
    B[2, 8] = -6xi / L^2
    B[2, 12] = (3xi + 1) / L
    B[3, 3] = 6xi / L^2
    B[3, 5] = (3xi - 1) / L
    B[3, 9] = -6xi / L^2
    B[3, 11] = (3xi + 1) / L
    B[4, 4] = -1 / L
    B[4, 10] = 1 / L
    return B
end

function beam_element_stiffness(L, C)
    K = zeros(12, 12)
    Csym = 0.5 .* (C .+ C')
    for xi in (-inv(sqrt(3.0)), inv(sqrt(3.0)))
        B = strain_displacement_matrix(L, xi)
        K .+= B' * Csym * B * (L / 2)
    end
    return 0.5 .* (K .+ K')
end

function element_constitutive(eta)
    conversion_moment = LBF_FT2_TO_NM2
    EIy = conversion_moment * linear_interpolate(ETA_STIFFNESS, EIYY_LBF_FT2, eta)
    EIz = conversion_moment * linear_interpolate(ETA_STIFFNESS, EIZZ_LBF_FT2, eta)
    EIzy = conversion_moment * linear_interpolate(ETA_STIFFNESS, EIZY_LBF_FT2, eta)
    GJ = conversion_moment * linear_interpolate(ETA_STIFFNESS, GJ_LBF_FT2, eta)
    EA = LBF_TO_N * linear_interpolate(ETA_STIFFNESS, EA_LBF, eta)
    return [EA 0.0 0.0 0.0; 0.0 EIz EIzy 0.0; 0.0 EIzy EIy 0.0; 0.0 0.0 0.0 GJ]
end

"""
Integrate an element whose constitutive matrix follows the complete workbook
distribution. The integration interval is split at every tabulated stiffness
breakpoint, and a three-point Gauss rule is used within each smooth piece.
This avoids the mesh-alignment oscillation caused by one midpoint sample.
"""
function distributed_beam_element_stiffness(x_left, x_right)
    L = x_right - x_left
    element_center = (x_left + x_right) / 2
    stiffness_stations = ETA_STIFFNESS .* SPAN_M
    internal_breaks = filter(x -> x_left < x < x_right, stiffness_stations)
    integration_intervals = unique(sort([x_left; internal_breaks; x_right]))
    gauss_points = (-sqrt(3 / 5), 0.0, sqrt(3 / 5))
    gauss_weights = (5 / 9, 8 / 9, 5 / 9)
    K = zeros(12, 12)

    for interval in 1:(length(integration_intervals) - 1)
        sub_left = integration_intervals[interval]
        sub_right = integration_intervals[interval + 1]
        sub_length = sub_right - sub_left
        sub_center = (sub_left + sub_right) / 2
        for (gauss, weight) in zip(gauss_points, gauss_weights)
            x = sub_center + gauss * sub_length / 2
            xi = 2 * (x - element_center) / L
            B = strain_displacement_matrix(L, xi)
            C = element_constitutive(x / SPAN_M)
            K .+= B' * C * B * (weight * sub_length / 2)
        end
    end
    return 0.5 .* (K .+ K')
end

function target_nodes(number_of_elements; mesh)
    number_of_elements > 0 || throw(ArgumentError("number_of_elements must be positive"))
    if mesh == :uniform
        return collect(range(0.0, SPAN_M, length = number_of_elements + 1))
    elseif mesh == :source
        number_of_elements == 16 || throw(ArgumentError("source mesh has exactly 16 elements"))
        return copy(SOURCE_NODES_M)
    end
    throw(ArgumentError("mesh must be :uniform or :source"))
end

"""
    assemble_wing_model(Ne; mesh=:uniform, remap=:conservative_spatial)

Assemble the fixed-root beam before root elimination.  Available remaps are:
`:conservative_spatial` (recommended), `:source_lumps` (the exact 16-record
reference, available only on the source mesh), `:legacy_linear` (the formerly
working componentwise interpolation), and `:signed_scaled` (the current
failing componentwise signed-total scaling, retained only for diagnosis).
"""
function assemble_wing_model(
    number_of_elements::Integer;
    mesh::Symbol = :uniform,
    remap::Symbol = :conservative_spatial,
    stiffness_model::Symbol = :distributed_integrated,
)
    nodes = target_nodes(number_of_elements; mesh)
    blocks_and_diagnostic = if remap == :conservative_spatial
        (blocks = conservative_spatial_remap(nodes),
         minimum_cg_inertia_eigenvalue = NaN,
         maximum_offset = NaN)
    elseif remap == :integrated_linear_density
        (blocks = integrated_linear_density_remap(nodes),
         minimum_cg_inertia_eigenvalue = NaN,
         maximum_offset = NaN)
    elseif remap == :control_volume
        (blocks = control_volume_spatial_remap(nodes),
         minimum_cg_inertia_eigenvalue = NaN,
         maximum_offset = NaN)
    elseif remap == :element_overlap
        (blocks = element_overlap_spatial_remap(nodes),
         minimum_cg_inertia_eigenvalue = NaN,
         maximum_offset = NaN)
    elseif remap == :source_lumps
        mesh == :source || throw(ArgumentError(
            ":source_lumps requires the 16-element source mesh",
        ))
        (blocks = vcat([zeros(6, 6)], source_spatial_inertia_blocks()),
         minimum_cg_inertia_eigenvalue = NaN,
         maximum_offset = NaN)
    elseif remap == :legacy_linear
        componentwise_remap(nodes; signed_rescale = false)
    elseif remap == :signed_scaled
        componentwise_remap(nodes; signed_rescale = true)
    else
        throw(ArgumentError("unknown remap $remap"))
    end
    blocks = blocks_and_diagnostic.blocks

    ndof = 6
    total_dofs = ndof * length(nodes)
    M = zeros(total_dofs, total_dofs)
    K = zeros(total_dofs, total_dofs)
    for node in eachindex(nodes)
        dofs = (6node - 5):(6node)
        M[dofs, dofs] .+= blocks[node]
    end
    for element in 1:number_of_elements
        L = nodes[element + 1] - nodes[element]
        eta = ((nodes[element + 1] + nodes[element]) / 2) / SPAN_M
        dofs = (6element - 5):(6element + 6)
        if stiffness_model == :distributed_integrated
            K[dofs, dofs] .+= distributed_beam_element_stiffness(
                nodes[element], nodes[element + 1],
            )
        elseif stiffness_model == :midpoint
            C = element_constitutive(eta)
            K[dofs, dofs] .+= beam_element_stiffness(L, C)
        else
            throw(ArgumentError(
                "stiffness_model must be :distributed_integrated or :midpoint",
            ))
        end
    end
    M .= 0.5 .* (M .+ M')
    K .= 0.5 .* (K .+ K')
    free_dofs = 7:total_dofs
    return (;
        number_of_elements,
        mesh,
        remap,
        stiffness_model,
        nodes,
        blocks,
        M_full = M,
        K_full = K,
        free_dofs,
        M = M[free_dofs, free_dofs],
        K = K[free_dofs, free_dofs],
        minimum_cg_inertia_eigenvalue = blocks_and_diagnostic.minimum_cg_inertia_eigenvalue,
        maximum_offset = blocks_and_diagnostic.maximum_offset,
    )
end

"""Add the current rigid-pylon pitch/yaw coordinates at exact eta = 0.83."""
function assemble_coupled_model(wing; attachment_eta = 0.83)
    nodes = wing.nodes
    attachment = attachment_eta * SPAN_M
    right = clamp(searchsortedfirst(nodes, attachment), 2, length(nodes))
    left = right - 1
    alpha = (attachment - nodes[left]) / (nodes[right] - nodes[left])
    A = zeros(6, size(wing.M_full, 1))
    A[:, (6left - 5):(6left)] .= (1 - alpha) .* Matrix{Float64}(I, 6, 6)
    A[:, (6right - 5):(6right)] .= alpha .* Matrix{Float64}(I, 6, 6)

    f_pitch = 7.97 * 2pi
    f_yaw = 7.97 * 2pi
    K_pitch = 19220.0
    K_yaw = 18916.0
    In_pitch = K_pitch / f_pitch^2
    In_yaw = K_yaw / f_yaw^2
    pylon_length = 5.6 * FT_TO_M
    pylon_mass = (0.0506 * SLUG_TO_KG / FT_TO_M) * pylon_length
    rotor_mass = 4 * 1.44
    total_propeller_mass = pylon_mass + rotor_mass
    S_root = -(rotor_mass * pylon_length + pylon_mass * pylon_length / 2)
    S_modal = -(rotor_mass * pylon_length / 2 + pylon_mass * pylon_length / 6)
    I_root = rotor_mass * pylon_length^2 + pylon_mass * pylon_length^2 / 3
    I_cross = rotor_mass * pylon_length^2 / 2 + pylon_mass * pylon_length^2 / 8

    B = [0.0 0.0 S_modal I_cross 0.0 0.0;
         -S_modal 0.0 0.0 0.0 0.0 I_cross]
    G = zeros(6, 6)
    G[1, 1] = total_propeller_mass
    G[2, 2] = total_propeller_mass
    G[3, 3] = total_propeller_mass
    G[3, 4] = S_root
    G[4, 3] = S_root
    G[4, 4] = I_root
    G[1, 6] = -S_root
    G[6, 1] = -S_root
    G[6, 6] = I_root

    wing_dofs = size(wing.M_full, 1)
    Mfull = zeros(wing_dofs + 2, wing_dofs + 2)
    Kfull = zeros(wing_dofs + 2, wing_dofs + 2)
    Mfull[1:wing_dofs, 1:wing_dofs] .= wing.M_full .+ A' * G * A
    Mfull[1:wing_dofs, (wing_dofs + 1):end] .= A' * B'
    Mfull[(wing_dofs + 1):end, 1:wing_dofs] .= B * A
    Mfull[(wing_dofs + 1):end, (wing_dofs + 1):end] .= [In_pitch 0.0; 0.0 In_yaw]
    Kfull[1:wing_dofs, 1:wing_dofs] .= wing.K_full
    Kfull[(wing_dofs + 1):end, (wing_dofs + 1):end] .= [K_pitch 0.0; 0.0 K_yaw]
    free_dofs = [collect(7:wing_dofs); wing_dofs + 1; wing_dofs + 2]
    return (;
        wing...,
        attachment_eta,
        attachment_nodes = (left, right),
        attachment_weights = (1 - alpha, alpha),
        M_full_coupled = Mfull,
        K_full_coupled = Kfull,
        M = Mfull[free_dofs, free_dofs],
        K = Kfull[free_dofs, free_dofs],
        free_dofs_coupled = free_dofs,
    )
end

function modal_analysis(model; number_of_modes = 12, classify_wing = false)
    solution = eigen(Symmetric(model.K), Symmetric(model.M))
    positive = findall(value -> isfinite(value) && value > 1e-9, real.(solution.values))
    positive = positive[sortperm(real.(solution.values[positive]))]
    positive = positive[1:min(number_of_modes, length(positive))]
    frequencies = [sqrt(real(solution.values[i])) / (2pi) for i in positive]
    classify_wing || return [
        (index = order, frequency_hz = frequencies[order]) for order in eachindex(frequencies)
    ]

    K = Matrix(Symmetric(model.K))
    wing_dofs = 6 * model.number_of_elements
    family_dofs = (
        out_of_plane = reduce(vcat, (start .+ [3, 5] for start in 0:6:(wing_dofs - 6))),
        in_plane = reduce(vcat, (start .+ [2, 6] for start in 0:6:(wing_dofs - 6))),
        torsion = collect(4:6:wing_dofs),
        axial = collect(1:6:wing_dofs),
    )
    counts = Dict(name => 0 for name in keys(family_dofs))
    modes = NamedTuple[]
    for (order, eigen_index) in enumerate(positive)
        shape = solution.vectors[:, eigen_index]
        energies = map(family_dofs) do indices
            projected = zeros(eltype(shape), length(shape))
            projected[indices] .= shape[indices]
            max(real(dot(projected, K * projected)), 0.0)
        end
        names = collect(keys(energies))
        values_for_families = collect(values(energies))
        family = names[argmax(values_for_families)]
        counts[family] += 1
        push!(modes, (;
            index = order,
            frequency_hz = frequencies[order],
            family,
            family_order = counts[family],
            family_energies = energies,
        ))
    end
    return modes
end

function remap_audit(wing)
    source_total = reduce(+, source_spatial_inertia_blocks())
    target_total = reduce(+, wing.blocks)
    relative_spatial_inertia_error = norm(target_total - source_total, Inf) /
        max(norm(source_total, Inf), 1.0)
    block_minimum_eigenvalue = minimum(
        minimum(eigvals(Symmetric(block))) for block in wing.blocks[2:end]
    )
    mass_minimum_eigenvalue = minimum(eigvals(Symmetric(wing.M)))
    stiffness_minimum_eigenvalue = minimum(eigvals(Symmetric(wing.K)))
    mass_kg = sum(tr(block[1:3, 1:3]) / 3 for block in wing.blocks)
    return (;
        mass_kg,
        relative_spatial_inertia_error,
        block_minimum_eigenvalue,
        mass_minimum_eigenvalue,
        stiffness_minimum_eigenvalue,
        minimum_cg_inertia_eigenvalue = wing.minimum_cg_inertia_eigenvalue,
        maximum_offset = wing.maximum_offset,
    )
end

end # module
