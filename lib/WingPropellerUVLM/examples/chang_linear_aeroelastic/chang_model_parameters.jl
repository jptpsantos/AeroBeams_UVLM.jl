# Chang-specific geometry, structural distributions, propeller properties,
# and simulation environment. Reusable algorithms belong in WingPropellerUVLM.

## ------------------ 3.1. Wing Parameters ------------------
Cr = WING_CONFIG.root_chord_m
Ct = WING_CONFIG.tip_chord_m
b = WING_CONFIG.span_m
span_length = b
taper_ratio = Ct / Cr
S_wing = b * (Cr + Ct) / 2.0
xle_tip = 0.50 * (Cr - Ct)
xle = [0.0, xle_tip]
yle = [0.0, b]
zle = [0.0, 0.0]
chord_geo = [Cr, Ct]
theta_geo = [0, 0]
phi_geo = [0, 0]
spacing_s_wing = Uniform()
spacing_c_wing = Uniform()
mirror_wing = false
symmetric_wing = false

ns_wing = WING_CONFIG.spanwise_panels
nc_wing = WING_CONFIG.chordwise_panels

Ne = ns_wing
le = span_length / Ne
nnodes = Ne + 1
ndof = 6
NDOF = nnodes * ndof

X_nodes_new = collect(range(0, span_length, length=nnodes))
X_cen_new = [(i - 0.5) * le for i in 1:Ne]
eta_cen_new = X_cen_new ./ span_length

chord = collect(range(Cr, Ct, length=nnodes))
xle_distribution = collect(range(xle[1], xle[2], length=nnodes))

span_nodes = X_nodes_new
L_node_new = fill(le, nnodes)
L_node_new[1] = le / 2.0
L_node_new[end] = le / 2.0

# --- 3.2 Nodal Mass Data Processing ---
chang_ft_to_m = 0.3048
chang_slug_to_kg = 14.5939029372064
chang_slug_ft2_to_kgm2 = chang_slug_to_kg * chang_ft_to_m^2
chang_lbf_to_n = 4.4482216152605
chang_lbf_ft2_to_nm2 = chang_lbf_to_n * chang_ft_to_m^2

X_nodes_ft = [0.0, 0.54675, 1.0936, 2.18722, 3.28084, 4.51115, 5.74147, 6.97178, 8.2021, 10.2526, 12.30315, 14.353675, 16.4042, 18.4547, 20.50525, 22.55577, 24.6063]
X_nodes_orig_unscaled = X_nodes_ft .* chang_ft_to_m
# The workbook tip coordinate differs from 7.5 m only by printed precision.
# Scaling all source stations by the same factor makes the remap cover exactly
# the active [0, span_length] interval and also keeps the data usable if the
# reference span is changed parametrically.
X_nodes_orig = X_nodes_orig_unscaled .* (span_length / X_nodes_orig_unscaled[end])
L_orig_elements = diff(X_nodes_orig)

L_node_orig = zeros(17)
L_node_orig[1] = L_orig_elements[1] / 2.0
for i in 2:16
    L_node_orig[i] = (L_orig_elements[i-1] + L_orig_elements[i]) / 2.0
end
L_node_orig[17] = L_orig_elements[16] / 2.0

Mass_slug = [0.4276211, 0.4139885, 0.736678, 0.4139886, 0.9436633, 0.7336013, 0.3851399, 0.9932094, 0.9932177, 0.9932237, 0.9931985, 0.971823, 0.9687671, 0.9109462, 0.9367904, 0.3350064]
Ixx_slug  = [0.6055, 0.6281, 1.92, 0.6281, 2.3536, 1.909, 0.5645, 2.4813, 2.4813, 2.4814, 2.4813, 2.3998, 2.3877, 2.2748, 2.3653, 0.459]
# Keep the inertia components in the literal Model Properties.xlsx column
# order: X is spanwise, Y is chordwise, and Z is vertical.  The later
# structural-to-aerodynamic adapter performs the required basis conversion;
# applying another Y/Z permutation here is incorrect.
Iyy_slug  = [0.0969, 0.1028, 0.1269, 0.1028, 0.2878, 0.1266, 0.0967, 0.3565, 0.3656, 0.3656, 0.3656, 0.3577, 0.3561, 0.3028, 0.3572, 0.0797]
Izz_slug  = [0.6475, 0.6777, 1.9703, 0.6777, 2.3194, 1.9593, 0.6118, 2.388, 2.388, 2.388, 2.3879, 2.3118, 2.3009, 2.2058, 2.2509, 0.4844]
Ixy_slug  = [0.000483499, 0.000545311, -0.003364902, 0.000546007, -0.003767303, -0.003368369, 0.000500823, -0.004768277, -0.004781691, -0.004764693, -0.00476285, -0.004980494, -0.004983901, -0.005318305, -0.005590199, 0.000500334]
Ixz_slug  = [0.001076907, -1.83141E-07, 0.001721633, -4.2213E-07, 0.004030091, -0.001715619, 0.000169958, 0.002555006, 0.002553573, 0.002555579, 0.002557012, 0.002683975, 0.00262101, 0.003720733, 0.004553367, -0.001762152]
Iyz_slug  = [-0.000821912, 6.13447E-06, 0.026378458, 7.88634E-06, 0.061766804, -0.026536753, -7.79545E-06, 0.031499987, 0.031497299, 0.03149642, 0.031500523, 0.030281192, 0.03089678, 0.024007273, 0.035641812, 0.001067331]

cg_x_ft = [0.027173898, 0.037438118, 0.021032451, 0.037438576, 0.02464125, 0.021088361, 0.040243862, 0.031203696, 0.03121296, 0.03120211, 0.031192549, 0.031900145, 0.031983194, 0.034021895, 0.03582014, 0.039275863]
cg_y_ft = [-0.66071232, -0.668605558, -0.919885764, -0.668609685, -0.864794967, -0.921385197, -0.671753688, -0.857220095, -0.857206007, -0.857200098, -0.857224387, -0.863614669, -0.864554528, -0.874965167, -0.869545335, -0.664042589]
cg_z_ft = [-0.245439677, -0.17398455, 0.016601162, -0.014913878, 0.121183138, -0.199929292, 0.137702621, -0.081678808, -0.059307771, -0.03693769, -0.014573273, 0.01200473, 0.03224826, 0.09014999, 0.048423281, 0.248940617]

M_orig_nodes = [0.0; Mass_slug .* chang_slug_to_kg]
Ixx_orig_nodes = [0.0; Ixx_slug .* chang_slug_ft2_to_kgm2]
Iyy_orig_nodes = [0.0; Iyy_slug .* chang_slug_ft2_to_kgm2]
Izz_orig_nodes = [0.0; Izz_slug .* chang_slug_ft2_to_kgm2]
Ixy_orig_nodes = [0.0; Ixy_slug .* chang_slug_ft2_to_kgm2]
Ixz_orig_nodes = [0.0; Ixz_slug .* chang_slug_ft2_to_kgm2]
Iyz_orig_nodes = [0.0; Iyz_slug .* chang_slug_ft2_to_kgm2]

cg_x_orig_nodes = [cg_x_ft[1] * chang_ft_to_m; cg_x_ft .* chang_ft_to_m]
cg_y_orig_nodes = [cg_y_ft[1] * chang_ft_to_m; cg_y_ft .* chang_ft_to_m]
cg_z_orig_nodes = [cg_z_ft[1] * chang_ft_to_m; cg_z_ft .* chang_ft_to_m]

chang_skew3(v) = [
     0.0  -v[3]   v[2];
     v[3]  0.0   -v[1];
    -v[2]  v[1]   0.0
]

"""Build a complete nodal spatial inertia about the beam reference axis."""
function chang_spatial_inertia_block(mass, inertia_at_cg, offset)
    skew_offset = chang_skew3(offset)
    identity3 = Matrix{Float64}(I, 3, 3)
    rotational_inertia = inertia_at_cg .- mass .* skew_offset * skew_offset
    return [mass .* identity3       -mass .* skew_offset;
            mass .* skew_offset      rotational_inertia]
end

"""
    chang_control_volume_spatial_remap(source_nodes, source_blocks, target_nodes)

Convert the source nodal spatial inertias to a piecewise-linear matrix density
and integrate it exactly over target nodal control volumes.  Complete 6-by-6
blocks are transferred with nonnegative weights, so mass, first moments, all
products of inertia, symmetry, and positive semidefiniteness remain coupled.
The root control-volume contribution is moved to the first free node because
the root DOFs are constrained and the reference root inertia is zero.
"""
function chang_control_volume_spatial_remap(
    source_nodes::AbstractVector{<:Real},
    source_blocks::AbstractVector{<:AbstractMatrix},
    target_nodes::AbstractVector{<:Real},
)
    length(source_nodes) == length(source_blocks) || throw(DimensionMismatch(
        "source_nodes and source_blocks must have the same length",
    ))
    issorted(source_nodes) || throw(ArgumentError("source_nodes must be sorted"))
    issorted(target_nodes) || throw(ArgumentError("target_nodes must be sorted"))
    isapprox(target_nodes[1], source_nodes[1]; atol = 1e-12, rtol = 0.0) ||
        throw(ArgumentError("source and target meshes must have the same root"))
    isapprox(target_nodes[end], source_nodes[end]; atol = 1e-10, rtol = 0.0) ||
        throw(ArgumentError("source and target meshes must have the same tip"))

    source_lengths = zeros(length(source_nodes))
    source_element_lengths = diff(source_nodes)
    source_lengths[1] = source_element_lengths[1] / 2
    source_lengths[end] = source_element_lengths[end] / 2
    for node in 2:(length(source_nodes) - 1)
        source_lengths[node] =
            (source_element_lengths[node - 1] + source_element_lengths[node]) / 2
    end
    source_density = [
        Matrix{Float64}(source_blocks[node]) ./ source_lengths[node]
        for node in eachindex(source_nodes)
    ]

    control_edges = [
        target_nodes[1];
        (target_nodes[1:end-1] .+ target_nodes[2:end]) ./ 2;
        target_nodes[end]
    ]
    target_blocks = [zeros(6, 6) for _ in eachindex(target_nodes)]

    for target_node in eachindex(target_nodes)
        target_left = control_edges[target_node]
        target_right = control_edges[target_node + 1]
        target_block = target_blocks[target_node]

        for source_element in 1:(length(source_nodes) - 1)
            source_left = source_nodes[source_element]
            source_right = source_nodes[source_element + 1]
            overlap_left = max(target_left, source_left)
            overlap_right = min(target_right, source_right)
            overlap_right <= overlap_left && continue

            interval = source_right - source_left
            t_left = (overlap_left - source_left) / interval
            t_right = (overlap_right - source_left) / interval
            left_weight = interval * (
                (t_right - t_right^2 / 2) - (t_left - t_left^2 / 2)
            )
            right_weight = interval * (t_right^2 - t_left^2) / 2
            target_block .+= left_weight .* source_density[source_element]
            target_block .+= right_weight .* source_density[source_element + 1]
        end
        target_block .= 0.5 .* (target_block .+ target_block')
    end

    # The clamped root block would be discarded during boundary-condition
    # elimination. Move it to the first free station to retain the complete
    # reconstructed wing inertia in the dynamic model.
    target_blocks[2] .+= target_blocks[1]
    fill!(target_blocks[1], 0.0)
    return target_blocks
end

source_spatial_inertia_blocks = Matrix{Float64}[]
for node in eachindex(X_nodes_orig)
    inertia_at_cg = [
        Ixx_orig_nodes[node] Ixy_orig_nodes[node] Ixz_orig_nodes[node];
        Ixy_orig_nodes[node] Iyy_orig_nodes[node] Iyz_orig_nodes[node];
        Ixz_orig_nodes[node] Iyz_orig_nodes[node] Izz_orig_nodes[node]
    ]
    offset = [
        cg_x_orig_nodes[node],
        cg_y_orig_nodes[node],
        cg_z_orig_nodes[node],
    ]
    push!(
        source_spatial_inertia_blocks,
        chang_spatial_inertia_block(M_orig_nodes[node], inertia_at_cg, offset),
    )
end

spatial_inertia_node_blocks = chang_control_volume_spatial_remap(
    X_nodes_orig,
    source_spatial_inertia_blocks,
    X_nodes_new,
)

function chang_properties_from_spatial_inertia(block)
    mass = tr(block[1:3, 1:3]) / 3
    mass <= 100eps(Float64) && return (
        mass = 0.0,
        cg = zeros(3),
        inertia_at_cg = zeros(3, 3),
    )
    skew_offset = -block[1:3, 4:6] ./ mass
    skew_offset = 0.5 .* (skew_offset .- skew_offset')
    offset = [skew_offset[3, 2], skew_offset[1, 3], skew_offset[2, 1]]
    inertia_at_cg = block[4:6, 4:6] .+ mass .* skew_offset * skew_offset
    inertia_at_cg = 0.5 .* (inertia_at_cg .+ inertia_at_cg')
    return (; mass, cg = offset, inertia_at_cg)
end

remapped_inertia_properties =
    chang_properties_from_spatial_inertia.(spatial_inertia_node_blocks)
m_node_vec = [property.mass for property in remapped_inertia_properties]
cg_x_node_vec = [property.cg[1] for property in remapped_inertia_properties]
cg_y_node_vec = [property.cg[2] for property in remapped_inertia_properties]
cg_z_node_vec = [property.cg[3] for property in remapped_inertia_properties]
Ixx_node_vec = [property.inertia_at_cg[1, 1] for property in remapped_inertia_properties]
Iyy_node_vec = [property.inertia_at_cg[2, 2] for property in remapped_inertia_properties]
Izz_node_vec = [property.inertia_at_cg[3, 3] for property in remapped_inertia_properties]
Ixy_node_vec = [property.inertia_at_cg[1, 2] for property in remapped_inertia_properties]
Ixz_node_vec = [property.inertia_at_cg[1, 3] for property in remapped_inertia_properties]
Iyz_node_vec = [property.inertia_at_cg[2, 3] for property in remapped_inertia_properties]

# Retain conventional first-moment and beam-axis inertia arrays for diagnostics
# and for the explicit legacy-inertia comparison.
mass_moment_x_node = m_node_vec .* cg_x_node_vec
mass_moment_y_node = m_node_vec .* cg_y_node_vec
mass_moment_z_node = m_node_vec .* cg_z_node_vec
Ixx_beam_orig = [block[4, 4] for block in source_spatial_inertia_blocks]
Iyy_beam_orig = [block[5, 5] for block in source_spatial_inertia_blocks]
Izz_beam_orig = [block[6, 6] for block in source_spatial_inertia_blocks]
Ixy_beam_orig = [block[4, 5] for block in source_spatial_inertia_blocks]
Ixz_beam_orig = [block[4, 6] for block in source_spatial_inertia_blocks]
Iyz_beam_orig = [block[5, 6] for block in source_spatial_inertia_blocks]
Ixx_beam_node = [block[4, 4] for block in spatial_inertia_node_blocks]
Iyy_beam_node = [block[5, 5] for block in spatial_inertia_node_blocks]
Izz_beam_node = [block[6, 6] for block in spatial_inertia_node_blocks]
Ixy_beam_node = [block[4, 5] for block in spatial_inertia_node_blocks]
Ixz_beam_node = [block[4, 6] for block in spatial_inertia_node_blocks]
Iyz_beam_node = [block[5, 6] for block in spatial_inertia_node_blocks]

source_spatial_inertia_total = reduce(+, source_spatial_inertia_blocks)
remapped_spatial_inertia_total = reduce(+, spatial_inertia_node_blocks)
wing_spatial_inertia_relative_error = norm(
    remapped_spatial_inertia_total - source_spatial_inertia_total,
    Inf,
) / max(norm(source_spatial_inertia_total, Inf), 1.0)
wing_spatial_inertia_relative_error <= 1e-12 || error(
    "Wing spatial-inertia remap is not conservative: " *
    "relative error=$wing_spatial_inertia_relative_error",
)
all(m_node_vec[2:end] .> 0.0) || error(
    "Every free structural wing node must have a positive remapped mass",
)
wing_spatial_block_minimum_eigenvalue = minimum(
    minimum(eigvals(Symmetric(block)))
    for block in spatial_inertia_node_blocks[2:end]
)
wing_spatial_block_minimum_eigenvalue > 0.0 || error(
    "A remapped free-node spatial inertia is not positive definite: " *
    "minimum eigenvalue=$wing_spatial_block_minimum_eigenvalue",
)
println(
    "Wing inertia remap: control-volume spatial blocks, Ne=$Ne, " *
    "mass=$(sum(m_node_vec)) kg, relative error=$wing_spatial_inertia_relative_error, " *
    "minimum block eigenvalue=$wing_spatial_block_minimum_eigenvalue",
)

# --- 3.3 Full 28-Element Stiffness Data ---
eta_stiff = [0.04443, 0.04444, 0.08888, 0.08889, 0.13332, 0.13333, 0.18332, 0.18333, 0.23332, 0.23333, 0.28332, 0.28333, 0.33333, 0.33334, 0.41666, 0.41667, 0.49999, 0.5, 0.58332, 0.58333, 0.66666, 0.66667, 0.74999, 0.75, 0.83333, 0.83334, 0.91667, 1.0]

EIyy_lbf = [32836890.0, 32836890.0, 32836890.0, 32836890.0, 32836890.0, 20120070.0, 20120070.0, 20120070.0, 20120070.0, 16181100.0, 16181100.0, 16181100.0, 16181100.0, 19251360.0, 19251360.0, 19251360.0, 19251360.0, 18286380.0, 18286380.0, 18286380.0, 18286380.0, 17183860.0, 17183860.0, 17183860.0, 17183860.0, 11198940.0, 11198940.0, 11198940.0]
EIzz_lbf = [325381340.0, 325381340.0, 325381340.0, 325381340.0, 325381340.0, 208531310.0, 208531310.0, 208531310.0, 208531310.0, 163131180.0, 163131180.0, 163131180.0, 163131180.0, 194913580.0, 194913580.0, 194913580.0, 194913580.0, 181528760.0, 181528760.0, 181528760.0, 181528760.0, 163800810.0, 163800810.0, 163800810.0, 163800810.0, 128273870.0, 128273870.0, 128273870.0]
EIzy_lbf = [4.7103163E+06, 4.7103163E+06, 4.7103163E+06, 4.7103163E+06, 4.7103163E+06, 2.9685812E+06, 2.9685812E+06, 2.9685812E+06, 2.9685812E+06, 2.5060481E+06, 2.5060481E+06, 2.5060481E+06, 2.5060481E+06, 3.3875434E+06, 3.3875434E+06, 3.3875434E+06, 3.3875434E+06, 3.1960560E+06, 3.1960560E+06, 3.1960560E+06, 3.1960560E+06, 3.3050893E+06, 3.3050893E+06, 3.3050893E+06, 3.3050893E+06, 2.4234529E+06, 2.4234529E+06, 2.4234529E+06]
GJ_lbf   = [3.05E+07, 3.05E+07, 3.05E+07, 3.05E+07, 3.05E+07, 1.84E+07, 1.84E+07, 1.84E+07, 1.84E+07, 1.40E+07, 1.40E+07, 1.40E+07, 1.40E+07, 1.62E+07, 1.62E+07, 1.62E+07, 1.62E+07, 1.52E+07, 1.52E+07, 1.52E+07, 1.52E+07, 1.30E+07, 1.30E+07, 1.30E+07, 1.30E+07, 9.82E+06, 9.82E+06, 9.82E+06]
EA_lbf   = [168148740.0, 168148740.0, 168148740.0, 168148740.0, 168148740.0, 111286380.0, 111286380.0, 111286380.0, 111286380.0, 85691907.0, 85691907.0, 85691907.0, 85691907.0, 102776990.0, 102776990.0, 102776990.0, 102776990.0, 96393539.0, 96393539.0, 96393539.0, 96393539.0, 88126997.0, 88126997.0, 88126697.0, 88126997.0, 100925580.0, 100925580.0, 100925580.0]

EIy_si = EIyy_lbf .* chang_lbf_ft2_to_nm2
EIz_si = EIzz_lbf .* chang_lbf_ft2_to_nm2
EIzy_si = EIzy_lbf .* chang_lbf_ft2_to_nm2
GJ_si  = GJ_lbf .* chang_lbf_ft2_to_nm2
EA_si  = EA_lbf .* chang_lbf_to_n

EIy_vec = linear_interpolate_1d(eta_stiff, EIy_si, eta_cen_new)
EIz_vec = linear_interpolate_1d(eta_stiff, EIz_si, eta_cen_new)
EIzy_vec = linear_interpolate_1d(eta_stiff, EIzy_si, eta_cen_new)
GJ_vec  = linear_interpolate_1d(eta_stiff, GJ_si, eta_cen_new)
EA_vec  = linear_interpolate_1d(eta_stiff, EA_si, eta_cen_new)

# The element-center vectors remain available for diagnostics and legacy
# comparisons. The active structural assembly receives the complete
# distribution below and integrates B'C(x)B through all stiffness breakpoints.
wing_stiffness_distribution = (
    eta = eta_stiff,
    EIy = EIy_si,
    EIz = EIz_si,
    EIzy = EIzy_si,
    GJ = GJ_si,
    EA = EA_si,
)

## ------------------ 3.2. Propeller Parameters ------------------
Vinf = SIMULATION_CONFIG.freestream_speed_mps  # increase to find the flutter speed
R_prop = PROPELLER_CONFIG.radius_m
Ω_wind = PROPELLER_CONFIG.rotation_rpm
V_trim = PROPELLER_CONFIG.trim_speed_mps
μ_prop = V_trim / (Ω_wind * 2 * pi / 60)/(R_prop)
J = μ_prop * pi
Ω = Vinf / R_prop / μ_prop
Nb_prop = PROPELLER_CONFIG.blades
c_prop = PROPELLER_CONFIG.chord_m

fθ_prop = 7.97 * 2 * pi
fψ_prop = 7.97 * 2 * pi
Kθ_prop = 19220.0
Kψ_prop = 18916.0

f_twist = 12.73 * 2 * pi
K_twist = 16835.0

ξ_prop = 0.0
# Global stiffness-proportional Rayleigh damping. For C = βK, the modal
# damping is ξ(ω) = βω/2. Calibrate β to 1% at the nominal 7.97 Hz pylon
# pitch frequency; lower-frequency modes receive proportionally less damping.
stiffness_damping_ratio = 0.01
stiffness_damping_reference_omega = fθ_prop
η0_prop = 0.0 / R_prop

ns_prop = PROPELLER_CONFIG.radial_panels
nc_prop = PROPELLER_CONFIG.chordwise_panels

# The Chang data already define the full blade-angle distribution. Keep this
# offset at zero unless a deliberate collective-pitch variation is required.
collective_pitch_offset_deg = 0.0
nodal_radii, _, twists_at_nodes = get_nodal_properties_chang(ns_prop)
blade_twists_prop = 1.0 .* (
    (twists_at_nodes .- 90 .+ collective_pitch_offset_deg) .|> deg2rad
)

slug_to_kg = chang_slug_to_kg
ft_to_m = chang_ft_to_m

L_pylon = 5.6 * ft_to_m
m_pylon = (0.0506 * slug_to_kg / ft_to_m) * L_pylon
m_blade = 1.44
m_rotor = Nb_prop * m_blade
mP_prop = m_pylon + m_rotor

cg_rotor = -L_pylon
cg_pylon = -0.5 * L_pylon

S_root = -(m_rotor * L_pylon + m_pylon * L_pylon / 2.0)
S_modal = -(m_rotor * L_pylon / 2.0 + m_pylon * L_pylon / 6.0)
I_root = m_rotor * L_pylon^2 + m_pylon * L_pylon^2 / 3.0
I_cross = m_rotor * L_pylon^2 / 2.0 + m_pylon * L_pylon^2 / 8.0

Inθ_prop = Kθ_prop / (fθ_prop^2)
Inψ_prop = Kψ_prop / (fψ_prop^2)
Ix_prop  = K_twist / (f_twist^2)     # positive axial inertia; Cgyro carries the sign convention

d_EA_to_pivot = 0.0

SθP_prop = S_modal
SψP_prop = S_modal
SαP_prop = S_root + mP_prop * d_EA_to_pivot
SγP_prop = S_root + mP_prop * d_EA_to_pivot
IθαP_prop = I_cross + d_EA_to_pivot * SθP_prop
IψγP_prop = I_cross + d_EA_to_pivot * SψP_prop
IαP_prop = I_root + 2.0 * d_EA_to_pivot * SθP_prop + mP_prop * d_EA_to_pivot^2
IγP_prop = I_root + 2.0 * d_EA_to_pivot * SψP_prop + mP_prop * d_EA_to_pivot^2

println("Flexible-pylon modal structural definition:")
println("  L_pylon = $(round(L_pylon, digits=6)) m,  mP = $(round(mP_prop, digits=6)) kg")
println("  S_root = $(round(S_root, digits=6)), S_modal = $(round(S_modal, digits=6))")
println("  I_root = $(round(I_root, digits=6)), I_cross = $(round(I_cross, digits=6))")
println("  Ix = $(round(Ix_prop, digits=6))  (positive axial inertia; Cgyro sign follows spin/yaw convention)")

propeller_eta = PROPELLER_CONFIG.attachment_eta
propeller_span_positions = propeller_eta .* span_length
Npropellers = length(propeller_span_positions)
prop_attach_nodes = [clamp(round(Int, Lp / le) + 1, 1, nnodes) for Lp in propeller_span_positions]
ndof_P = 2 * Npropellers

# Work-conjugate attachment interpolation. The nearest-node indices are kept
# for compatibility with UVLM allocation utilities, while the structural
# matrices, propeller kinematics, and load transfer use these exact brackets.
prop_attachment_node_pairs = Tuple{Int,Int}[]
prop_attachment_weights = Tuple{Float64,Float64}[]
for position in propeller_span_positions
    right_node = clamp(searchsortedfirst(span_nodes, position), 2, nnodes)
    left_node = right_node - 1
    fraction = (position - span_nodes[left_node]) /
        (span_nodes[right_node] - span_nodes[left_node])
    push!(prop_attachment_node_pairs, (left_node, right_node))
    push!(prop_attachment_weights, (1 - fraction, fraction))
end

for ip in 1:Npropellers
    left_node, right_node = prop_attachment_node_pairs[ip]
    left_weight, right_weight = prop_attachment_weights[ip]
    println(
        "Propeller P$(ip) at η=$(propeller_eta[ip]), y=$(propeller_span_positions[ip]) m; " *
        "structural nodes=($left_node,$right_node), " *
        "weights=($(round(left_weight, digits=6)),$(round(right_weight, digits=6)))",
    )
end
structural_pivot_prop = SVector(0.0, 0.0, 0.0)

## ------------------ 3.4. Simulation Environment & Time ------------------
azimutal_step_deg = SIMULATION_CONFIG.azimuth_step_deg
t_end = parse(Float64, get(ENV, "CHANG_END_TIME_S", string(SIMULATION_CONFIG.end_time_s)))
t_end > 0 || error("CHANG_END_TIME_S must be positive")
t_step = (azimutal_step_deg * pi / 180) / Ω
t = range(0.0, t_end, step=t_step)
dt = fill(t_step, length(t) - 1)
println("Timestep based on propeller rotation: dt = $t_step s")

INTERACTION_ON = SIMULATION_CONFIG.interaction_on
INTERACTION_MATRIX = nothing
println("Aerodynamic Interaction: $INTERACTION_ON")

alpha_deg = SIMULATION_CONFIG.angle_of_attack_deg
alpha = alpha_deg * pi / 180
beta = SIMULATION_CONFIG.sideslip_deg * pi / 180
Omega_fs = [0.0, 0.0, 0.0]
fs = Freestream(Vinf, alpha, beta, Omega_fs)

Sref = S_wing
cref_wing = (2/3) * Cr * (1 + taper_ratio + taper_ratio^2) / (1 + taper_ratio)
cref = cref_wing
bref = b
rref = [0.3 * Cr, 0.0, 0.0]
ref = Reference(Sref, cref, bref, rref, Vinf, AIR_DENSITY)
