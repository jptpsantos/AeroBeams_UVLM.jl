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
X_nodes_ft = [0.0, 0.54675, 1.0936, 2.18722, 3.28084, 4.51115, 5.74147, 6.97178, 8.2021, 10.2526, 12.30315, 14.353675, 16.4042, 18.4547, 20.50525, 22.55577, 24.6063]
X_nodes_orig = X_nodes_ft .* 0.3048
L_orig_elements = diff(X_nodes_orig)

L_node_orig = zeros(17)
L_node_orig[1] = L_orig_elements[1] / 2.0
for i in 2:16
    L_node_orig[i] = (L_orig_elements[i-1] + L_orig_elements[i]) / 2.0
end
L_node_orig[17] = L_orig_elements[16] / 2.0

Mass_slug = [0.4276211, 0.4139885, 0.736678, 0.4139886, 0.9436633, 0.7336013, 0.3851399, 0.9932094, 0.9932177, 0.9932237, 0.9931985, 0.971823, 0.9687671, 0.9109462, 0.9367904, 0.3350064]
Ixx_slug  = [0.6055, 0.6281, 1.92, 0.6281, 2.3536, 1.909, 0.5645, 2.4813, 2.4813, 2.4814, 2.4813, 2.3998, 2.3877, 2.2748, 2.3653, 0.459]
# >>> FIX 3: Iyy/Izz, Ixy/Ixz, cg_y/cg_z assigned as in the validated reference
#            (main_vgf_chang.jl). The [0.0969] array is Izz, [0.6475] is Iyy, etc.
Izz_slug  = [0.0969, 0.1028, 0.1269, 0.1028, 0.2878, 0.1266, 0.0967, 0.3565, 0.3656, 0.3656, 0.3656, 0.3577, 0.3561, 0.3028, 0.3572, 0.0797]
Iyy_slug  = [0.6475, 0.6777, 1.9703, 0.6777, 2.3194, 1.9593, 0.6118, 2.388, 2.388, 2.388, 2.3879, 2.3118, 2.3009, 2.2058, 2.2509, 0.4844]
Ixz_slug  = [0.000483499, 0.000545311, -0.003364902, 0.000546007, -0.003767303, -0.003368369, 0.000500823, -0.004768277, -0.004781691, -0.004764693, -0.00476285, -0.004980494, -0.004983901, -0.005318305, -0.005590199, 0.000500334]
Ixy_slug  = [0.001076907, -1.83141E-07, 0.001721633, -4.2213E-07, 0.004030091, -0.001715619, 0.000169958, 0.002555006, 0.002553573, 0.002555579, 0.002557012, 0.002683975, 0.00262101, 0.003720733, 0.004553367, -0.001762152]
Iyz_slug  = [-0.000821912, 6.13447E-06, 0.026378458, 7.88634E-06, 0.061766804, -0.026536753, -7.79545E-06, 0.031499987, 0.031497299, 0.03149642, 0.031500523, 0.030281192, 0.03089678, 0.024007273, 0.035641812, 0.001067331]

cg_x_ft = [0.027173898, 0.037438118, 0.021032451, 0.037438576, 0.02464125, 0.021088361, 0.040243862, 0.031203696, 0.03121296, 0.03120211, 0.031192549, 0.031900145, 0.031983194, 0.034021895, 0.03582014, 0.039275863]
# >>> FIX 3: cg_z is the large (~-0.66 ft) offset; cg_y is the small (~-0.25 ft) one.
cg_z_ft = [-0.66071232, -0.668605558, -0.919885764, -0.668609685, -0.864794967, -0.921385197, -0.671753688, -0.857220095, -0.857206007, -0.857200098, -0.857224387, -0.863614669, -0.864554528, -0.874965167, -0.869545335, -0.664042589]
cg_y_ft = [-0.245439677, -0.17398455, 0.016601162, -0.014913878, 0.121183138, -0.199929292, 0.137702621, -0.081678808, -0.059307771, -0.03693769, -0.014573273, 0.01200473, 0.03224826, 0.09014999, 0.048423281, 0.248940617]

M_orig_nodes = [0.0; Mass_slug .* 14.5939]
Ixx_orig_nodes = [0.0; Ixx_slug .* 1.355818]
Iyy_orig_nodes = [0.0; Iyy_slug .* 1.355818]
Izz_orig_nodes = [0.0; Izz_slug .* 1.355818]
Ixy_orig_nodes = [0.0; Ixy_slug .* 1.355818]
Ixz_orig_nodes = [0.0; Ixz_slug .* 1.355818]
Iyz_orig_nodes = [0.0; Iyz_slug .* 1.355818]

cg_x_orig_nodes = [cg_x_ft[1]*0.3048; cg_x_ft .* 0.3048]
cg_y_orig_nodes = [cg_y_ft[1]*0.3048; cg_y_ft .* 0.3048]
cg_z_orig_nodes = [cg_z_ft[1]*0.3048; cg_z_ft .* 0.3048]

rho_M_orig = M_orig_nodes ./ L_node_orig
rho_Ixx_orig = Ixx_orig_nodes ./ L_node_orig
rho_Iyy_orig = Iyy_orig_nodes ./ L_node_orig
rho_Izz_orig = Izz_orig_nodes ./ L_node_orig
rho_Ixy_orig = Ixy_orig_nodes ./ L_node_orig
rho_Ixz_orig = Ixz_orig_nodes ./ L_node_orig
rho_Iyz_orig = Iyz_orig_nodes ./ L_node_orig

rho_M_new = linear_interpolate_1d(X_nodes_orig, rho_M_orig, X_nodes_new)
rho_Ixx_new = linear_interpolate_1d(X_nodes_orig, rho_Ixx_orig, X_nodes_new)
rho_Iyy_new = linear_interpolate_1d(X_nodes_orig, rho_Iyy_orig, X_nodes_new)
rho_Izz_new = linear_interpolate_1d(X_nodes_orig, rho_Izz_orig, X_nodes_new)
rho_Ixy_new = linear_interpolate_1d(X_nodes_orig, rho_Ixy_orig, X_nodes_new)
rho_Ixz_new = linear_interpolate_1d(X_nodes_orig, rho_Ixz_orig, X_nodes_new)
rho_Iyz_new = linear_interpolate_1d(X_nodes_orig, rho_Iyz_orig, X_nodes_new)

m_node_vec = rho_M_new .* L_node_new
Ixx_node_vec = rho_Ixx_new .* L_node_new
Iyy_node_vec = rho_Iyy_new .* L_node_new
Izz_node_vec = rho_Izz_new .* L_node_new
Ixy_node_vec = rho_Ixy_new .* L_node_new
Ixz_node_vec = rho_Ixz_new .* L_node_new
Iyz_node_vec = rho_Iyz_new .* L_node_new

cg_x_node_vec = linear_interpolate_1d(X_nodes_orig, cg_x_orig_nodes, X_nodes_new)
cg_y_node_vec = linear_interpolate_1d(X_nodes_orig, cg_y_orig_nodes, X_nodes_new)
cg_z_node_vec = linear_interpolate_1d(X_nodes_orig, cg_z_orig_nodes, X_nodes_new)

# --- 3.3 Full 28-Element Stiffness Data ---
eta_stiff = [0.04443, 0.04444, 0.08888, 0.08889, 0.13332, 0.13333, 0.18332, 0.18333, 0.23332, 0.23333, 0.28332, 0.28333, 0.33333, 0.33334, 0.41666, 0.41667, 0.49999, 0.5, 0.58332, 0.58333, 0.66666, 0.66667, 0.74999, 0.75, 0.83333, 0.83334, 0.91667, 1.0]

EIyy_lbf = [3.28E+07, 3.28E+07, 3.28E+07, 3.28E+07, 3.28E+07, 2.01E+07, 2.01E+07, 2.01E+07, 2.01E+07, 1.62E+07, 1.62E+07, 1.62E+07, 1.62E+07, 1.93E+07, 1.93E+07, 1.93E+07, 1.93E+07, 1.83E+07, 1.83E+07, 1.83E+07, 1.83E+07, 1.72E+07, 1.72E+07, 1.72E+07, 1.72E+07, 1.12E+07, 1.12E+07, 1.12E+07]
EIzz_lbf = [3.25E+08, 3.25E+08, 3.25E+08, 3.25E+08, 3.25E+08, 2.09E+08, 2.09E+08, 2.09E+08, 2.09E+08, 1.63E+08, 1.63E+08, 1.63E+08, 1.63E+08, 1.95E+08, 1.95E+08, 1.95E+08, 1.95E+08, 1.82E+08, 1.82E+08, 1.82E+08, 1.82E+08, 1.64E+08, 1.64E+08, 1.64E+08, 1.64E+08, 1.28E+08, 1.28E+08, 1.28E+08]
EIzy_lbf = [4.7103163E+06, 4.7103163E+06, 4.7103163E+06, 4.7103163E+06, 4.7103163E+06, 2.9685812E+06, 2.9685812E+06, 2.9685812E+06, 2.9685812E+06, 2.5060481E+06, 2.5060481E+06, 2.5060481E+06, 2.5060481E+06, 3.3875434E+06, 3.3875434E+06, 3.3875434E+06, 3.3875434E+06, 3.1960560E+06, 3.1960560E+06, 3.1960560E+06, 3.1960560E+06, 3.3050893E+06, 3.3050893E+06, 3.3050893E+06, 3.3050893E+06, 2.4234529E+06, 2.4234529E+06, 2.4234529E+06]
GJ_lbf   = [3.05E+07, 3.05E+07, 3.05E+07, 3.05E+07, 3.05E+07, 1.84E+07, 1.84E+07, 1.84E+07, 1.84E+07, 1.40E+07, 1.40E+07, 1.40E+07, 1.40E+07, 1.62E+07, 1.62E+07, 1.62E+07, 1.62E+07, 1.52E+07, 1.52E+07, 1.52E+07, 1.52E+07, 1.30E+07, 1.30E+07, 1.30E+07, 1.30E+07, 9.82E+06, 9.82E+06, 9.82E+06]
EA_lbf   = [1.68E+08, 1.68E+08, 1.68E+08, 1.68E+08, 1.68E+08, 1.11E+08, 1.11E+08, 1.11E+08, 1.11E+08, 8.57E+07, 8.57E+07, 8.57E+07, 8.57E+07, 1.03E+08, 1.03E+08, 1.03E+08, 1.03E+08, 9.64E+07, 9.64E+07, 9.64E+07, 9.64E+07, 8.81E+07, 8.81E+07, 8.81E+07, 8.81E+07, 1.01E+08, 1.01E+08, 1.01E+08]

EIy_si = EIyy_lbf .* 0.413138
EIz_si = EIzz_lbf .* 0.413138
EIzy_si = EIzy_lbf .* 0.413138
GJ_si  = GJ_lbf .* 0.413138
EA_si  = EA_lbf .* 4.44822

EIy_vec = linear_interpolate_1d(eta_stiff, EIy_si, eta_cen_new)
EIz_vec = linear_interpolate_1d(eta_stiff, EIz_si, eta_cen_new)
EIzy_vec = linear_interpolate_1d(eta_stiff, EIzy_si, eta_cen_new)
GJ_vec  = linear_interpolate_1d(eta_stiff, GJ_si, eta_cen_new)
EA_vec  = linear_interpolate_1d(eta_stiff, EA_si, eta_cen_new)

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
η0_prop = 0.0 / R_prop

ns_prop = PROPELLER_CONFIG.radial_panels
nc_prop = PROPELLER_CONFIG.chordwise_panels

β75 = 0.0
nodal_radii, _, twists_at_nodes = get_nodal_properties_chang(ns_prop)
blade_twists_prop = 1.0 .* ((twists_at_nodes .- 90 .+ β75) .|> deg2rad)

slug_to_kg = 14.5939
ft_to_m = 0.3048

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

for ip in 1:Npropellers
    println("Propeller P$(ip) at η=$(propeller_eta[ip]) attached at node $(prop_attach_nodes[ip]) (y=$(round(span_nodes[prop_attach_nodes[ip]], digits=4)) m)")
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
ref = Reference(Sref, cref, bref, rref, Vinf, RHO)
