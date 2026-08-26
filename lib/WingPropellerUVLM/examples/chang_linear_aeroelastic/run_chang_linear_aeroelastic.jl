# ==============================================================================
# v9 UVLM time-domain model with frequency-domain structural definitions
#   (PARTITIONED GENERALIZED-ALPHA build — five corrections tagged "# >>> FIX")
#
#   FIX 1: Generalized-alpha partitioned predictor-corrector with iterative aerodynamic load correction (α_nb = 0.05).
#   FIX 2: propeller modal motion/load arm = -0.5*L_pylon
#          (reference eθ), while the physical hub remains at -L_pylon.
#   FIX 3: nodal inertia/CG axis assignment matched to the validated reference
#          (Iyy<->Izz, Ixy<->Ixz, cg_y<->cg_z).
#   FIX 4: INTERACTION_ON = false (isolated prop aero, as in the reference).
#   FIX 5: full structural deformation in the wing grid and full wing-node
#          rotation onto the propeller (Rz*Rx*Ry), not torsion-only.
#
# Structural nodal DOF order:
#   [u_span, v_chord, w_vertical_down, θ_span/torsion, θ_chord, θ_vertical]
# ==============================================================================

import Pkg

# Make the copied input directly runnable from VS Code or a terminal without
# installing WingPropellerUVLM in the user's global environment.
Pkg.activate(normpath(joinpath(@__DIR__, "..", "..")))

using LinearAlgebra
using StaticArrays
using DelimitedFiles
using Statistics

include(joinpath(@__DIR__, "chang_case.jl"))

const PLOT_RESULTS = environment_flag("CHANG_PLOT_RESULTS", true)
PLOT_RESULTS && (@eval using Plots)

using WingPropellerUVLM:
    Uniform,
    Freestream,
    Reference,
    RotationMatrix,
    SurfacePanel,
    initialize_bohnisch_uvlm_system,
    get_nodal_properties_chang,
    grid_to_surface_panels,
    copy_surfaces_to_previous!,
    propagate_system!,
    snapshot_uvlm,
    restore_uvlm!,
    imperial_nodal_forces,
    imperial_nodal_positions,
    generate_panel_grid_and_interpolate,
    linear_interpolate_1d,
    generalized_alpha_parameters,
    PartitionedCouplingOptions,
    partitioned_generalized_alpha_step,
    smooth_hann_pulse_load

# ==============================================================================
# 1. INCLUDES
# ==============================================================================
PLOT_RESULTS && include(joinpath(@__DIR__, "chang_plotting.jl"))
include(joinpath(@__DIR__, "chang_postprocessing.jl"))

RHO = 1.225 # Air density (kg/m^3)
const VERIFY_OUTPUT_DIR = joinpath(@__DIR__, "output")
const VERIFY_LABEL = "chang_linear_imperial_uvlm"
mkpath(VERIFY_OUTPUT_DIR)

include(joinpath(@__DIR__, "chang_structural_matrices.jl"))

# ==============================================================================
# 2. CHANG MODEL PARAMETERS
# ==============================================================================
include(joinpath(@__DIR__, "chang_model_parameters.jl"))


# ==============================================================================
# 3. STRUCTURAL MATRICES
# ==============================================================================
println("Assembling Chang structural matrices (Z-DOWN)...")

structural = assemble_chang_structural_matrices(
    Ne = Ne,
    le = le,
    ndof = ndof,
    NDOF = NDOF,
    nnodes = nnodes,
    EIy_vec = EIy_vec,
    EIz_vec = EIz_vec,
    EIzy_vec = EIzy_vec,
    GJ_vec = GJ_vec,
    EA_vec = EA_vec,
    m_node_vec = m_node_vec,
    Ixx_node_vec = Ixx_node_vec,
    Iyy_node_vec = Iyy_node_vec,
    Izz_node_vec = Izz_node_vec,
    Ixy_node_vec = Ixy_node_vec,
    Ixz_node_vec = Ixz_node_vec,
    Iyz_node_vec = Iyz_node_vec,
    cg_x_node_vec = cg_x_node_vec,
    cg_y_node_vec = cg_y_node_vec,
    cg_z_node_vec = cg_z_node_vec,
    ndof_P = ndof_P,
    Npropellers = Npropellers,
    prop_attach_nodes = prop_attach_nodes,
    Inθ_prop = Inθ_prop,
    Inψ_prop = Inψ_prop,
    Kθ_prop = Kθ_prop,
    Kψ_prop = Kψ_prop,
    ξ_prop = ξ_prop,
    Ix_prop = Ix_prop,
    Ω = Ω,
    mP_prop = mP_prop,
    SθP_prop = SθP_prop,
    SψP_prop = SψP_prop,
    SαP_prop = SαP_prop,
    SγP_prop = SγP_prop,
    IθαP_prop = IθαP_prop,
    IψγP_prop = IψγP_prop,
    IαP_prop = IαP_prop,
    IγP_prop = IγP_prop,
)

Ks_W = structural.Ks_W
Ms_W = structural.Ms_W
Cs_W = structural.Cs_W
Ms_P = structural.Ms_P
Cs_P = structural.Cs_P
Ks_P = structural.Ks_P
Bs_P = structural.Bs_P
Ds_P = structural.Ds_P
Fs_W = structural.Fs_W
Gs_W = structural.Gs_W
Hs_W = structural.Hs_W
attach_dofs_all = structural.attach_dofs_all
M_global = structural.M_global
C_global = structural.C_global
K_global = structural.K_global
free_dofs = structural.free_dofs
M = structural.M
C = structural.C
K = structural.K
ndof_free = structural.ndof_free
ndof_wing_free = structural.ndof_wing_free
ndof_prop_free = structural.ndof_prop_free
println("Matrices after BCs. Total DOFs (free): $ndof_free")

U = Vector{Vector{Float64}}(undef, length(t))
Ud = similar(U)
Udd = similar(U)
U0 = zeros(ndof_free)
U0d = zeros(ndof_free)
U0dd = isempty(M) ? zeros(ndof_free) : M \ (zeros(ndof_free) - C*U0d - K*U0)
U[1] = U0
Ud[1] = U0d
Udd[1] = U0dd
# ==============================================================================
# 4. UVLM INITIALIZATION
# ==============================================================================
println("Initializing global UVLM system...")
FCORE = (c, Δs) -> 0.5 * Δs

# Chang reference locations: the wing elastic axis is at 30% chord, the
# physical propeller hub is one pylon length ahead of the attachment, and the
# two pylon modal coordinates act at half that length.
elastic_axis_fraction = 0.30
prop_pivot_offset_from_ea_A = SVector(0.0, 0.0, 0.0)
hub_center_prop_A = SVector(-L_pylon, 0.0, 0.0)
hub_center_load_A = SVector(-0.5 * L_pylon, 0.0, 0.0)

uvlm = initialize_bohnisch_uvlm_system(
    xle=xle, yle=yle, zle=zle,
    chord_geo=chord_geo, theta_geo=theta_geo, phi_geo=phi_geo,
    ns_wing=ns_wing, nc_wing=nc_wing,
    mirror_wing=mirror_wing,
    spacing_s_wing=spacing_s_wing, spacing_c_wing=spacing_c_wing,
    R_prop=R_prop, c_prop=c_prop, ns_prop=ns_prop, nc_prop=nc_prop,
    blade_twists_prop=blade_twists_prop, Nb_prop=Nb_prop,
    Npropellers=Npropellers,
    span_nodes=span_nodes, prop_attach_nodes=prop_attach_nodes,
    chord=chord, xle_distribution=xle_distribution,
    ref=ref, symmetric_wing=symmetric_wing, fs=fs, dt=dt, nnodes=nnodes,
    prop_pivot_offset_from_ea_A=prop_pivot_offset_from_ea_A,
    hub_center_prop_A=hub_center_prop_A,
    fcore=FCORE,
    elastic_axis_fraction=elastic_axis_fraction,
    maximum_wake_rows_wing=10 * nc_wing,
    maximum_wake_rows_propeller=72,
    verbose=true,
)

ratio_wing = uvlm.ratio_wing
grids_prop_ref = uvlm.grids_prop_ref
grids_prop_initial_global = uvlm.grids_prop_initial_global
attach_node_y = uvlm.attach_node_y
ea_x_aero = uvlm.ea_x_aero
prop_surface_indices = uvlm.prop_surface_indices
surfaces = uvlm.surfaces
nsurf = uvlm.nsurf
surface_interaction_id = uvlm.surface_interaction_id
nwake = uvlm.nwake
system = uvlm.system
repeated_points = uvlm.repeated_points
iwake = uvlm.iwake
fs_vec = uvlm.fs_vec
save = uvlm.save
TF = uvlm.TF
surface_history = uvlm.surface_history
nodal_forces_wing = uvlm.nodal_forces_wing
nodal_moments_wing = uvlm.nodal_moments_wing
EA_nodes_wing = uvlm.EA_nodes_wing
nodal_forces_prop = uvlm.nodal_forces_prop
grids_prop_current = uvlm.grids_prop_current
T_pivot_A_current = uvlm.T_pivot_A_current
T_hub_A_current   = Vector{SVector{3,Float64}}(undef, Npropellers)
T_load_A_current  = Vector{SVector{3,Float64}}(undef, Npropellers)
include(joinpath(@__DIR__, "chang_uvlm_coupling.jl"))

# ==============================================================================
# 5. SIMULATION
# ==============================================================================
println("Starting coupled aeroelastic simulation...")

pitch_dof_indices = [ndof_wing_free + 2*(ip - 1) + 1 for ip in 1:Npropellers]
impulse_magnitude = parse(Float64, get(ENV, "CHANG_IMPULSE_MAGNITUDE", "1000.0"))
impulse_start_time = parse(Float64, get(ENV, "CHANG_IMPULSE_START_S", "0.2"))
impulse_duration = parse(Float64, get(ENV, "CHANG_IMPULSE_DURATION_S", "0.15"))
impulse_start_time >= 0 || error("CHANG_IMPULSE_START_S must be nonnegative")
impulse_duration > 0 || error("CHANG_IMPULSE_DURATION_S must be positive")
impulse_triggered_msg = false

F0_struct = zeros(ndof_free)
have_F0 = false
printed_steady = false
U_ABORT = 1.0e3
N_LAST = length(dt)

const GA_RHO_INF = parse(Float64, get(ENV, "CHANG_GA_RHO_INF", "0.7"))
const GA_PARAMS = generalized_alpha_parameters(GA_RHO_INF)
const COUPLING_MAX_ITER = parse(Int, get(ENV, "CHANG_COUPLING_MAX_ITER", "10"))
const COUPLING_TOL_U = parse(Float64, get(ENV, "CHANG_COUPLING_TOL_U", "1.0e-5"))
const COUPLING_TOL_F = parse(Float64, get(ENV, "CHANG_COUPLING_TOL_F", "1.0e-2"))
const COUPLING_TOL_EQ = parse(Float64, get(ENV, "CHANG_COUPLING_TOL_EQ", "1.0e-10"))
const COUPLING_RELAXATION = parse(Float64, get(ENV, "CHANG_COUPLING_RELAXATION", "0.5"))
const COUPLING_OPTIONS = PartitionedCouplingOptions(
    maximum_iterations = COUPLING_MAX_ITER,
    state_tolerance = COUPLING_TOL_U,
    load_tolerance = COUPLING_TOL_F,
    equilibrium_tolerance = COUPLING_TOL_EQ,
    relaxation = COUPLING_RELAXATION,
)
println("Partitioned generalized-alpha: rho_inf=$(GA_PARAMS.rho_inf), alpha_m=$(GA_PARAMS.alpha_m), alpha_f=$(GA_PARAMS.alpha_f), gamma=$(GA_PARAMS.gamma), beta=$(GA_PARAMS.beta)")
println("Coupling correction: max_iter=$COUPLING_MAX_ITER, tol_u=$COUPLING_TOL_U, tol_f=$COUPLING_TOL_F, tol_eq=$COUPLING_TOL_EQ, relaxation=$COUPLING_RELAXATION")

coupling_iterations = fill(0, length(dt))
coupling_disp_residual = fill(NaN, length(dt))
coupling_load_residual = fill(NaN, length(dt))
coupling_equilibrium_residual = fill(NaN, length(dt))
coupling_converged = fill(false, length(dt))
F_struct_n = zeros(ndof_free)
F_pert_n = zeros(ndof_free)

for it = 1:length(dt)
    copy_surfaces_to_previous!(system, nsurf)
    snap = snapshot_uvlm(system)
    dt_i = dt[it]

    f_ext_n = smooth_hann_pulse_load(
        t[it],
        ndof_free,
        pitch_dof_indices;
        magnitude = impulse_magnitude,
        start_time = impulse_start_time,
        duration = impulse_duration,
    )
    f_ext_np1 = smooth_hann_pulse_load(
        t[it + 1],
        ndof_free,
        pitch_dof_indices;
        magnitude = impulse_magnitude,
        start_time = impulse_start_time,
        duration = impulse_duration,
    )
    if !impulse_triggered_msg && (maximum(abs.(f_ext_n)) > 0.0 || maximum(abs.(f_ext_np1)) > 0.0)
        println("<<<<<<<<<<< Applying Smooth Pitch Impulse >>>>>>>>>>>")
        global impulse_triggered_msg = true
    end

    correction = partitioned_generalized_alpha_step(
        M,
        C,
        K,
        U[it],
        Ud[it],
        Udd[it],
        F_pert_n,
        f_ext_np1,
        f_ext_n,
        dt_i,
        GA_PARAMS,
        state_guess -> begin
            aerodynamic_load = aero_load_for_state!(
                system,
                snap,
                state_guess,
                it;
                print_loads = false,
            )
            return have_F0 ? aerodynamic_load .- F0_struct : zeros(ndof_free)
        end;
        options = COUPLING_OPTIONS,
        require_load_convergence = have_F0,
    )

    U_corr = correction.displacement
    Ud_corr = correction.velocity
    Udd_corr = correction.acceleration
    F_pert_guess = correction.trial_load
    state_res = correction.state_residual
    force_res = correction.load_residual
    equilibrium_res = correction.equilibrium_residual
    converged = correction.converged
    iter_count = correction.iterations

    final_print = (it <= 5 || it % 500 == 0)
    U_final = copy(U_corr)
    Ud_final = copy(Ud_corr)
    Udd_final = copy(Udd_corr)

    # Recompute once at the accepted state and leave this call's wake and
    # circulation as the committed state for the next time step. The trial
    # calls above were restored from `snap` and therefore do not advance the
    # aerodynamic history multiple times.
    F_struct_final = aero_load_for_state!(system, snap, U_final, it; print_loads = final_print)
    F_pert_final = have_F0 ? (F_struct_final .- F0_struct) : zeros(ndof_free)
    commit_force_res = norm(F_pert_final .- F_pert_guess) / max(norm(F_pert_final), 1.0)
    force_res = max(force_res, commit_force_res)

    if !converged && state_res <= 10.0 * COUPLING_TOL_U &&
        equilibrium_res <= 10.0 * COUPLING_TOL_EQ &&
        (!have_F0 || commit_force_res <= 10.0 * COUPLING_TOL_F)
        converged = true
    elseif commit_force_res > 10.0 * COUPLING_TOL_F && have_F0
        converged = false
    end

    U[it+1] = U_final
    Ud[it+1] = Ud_final
    Udd[it+1] = Udd_final

    trim_ready = t[it] >= impulse_start_time - 2 * dt_i
    if !have_F0 && trim_ready
        F0_struct .= F_struct_final
        global have_F0 = true
        F_pert_final .= 0.0
        if !printed_steady
            println("\n=== TRIM generalized load captured at t=$(round(t[it], digits=4)) s ===")
            for ip in 1:Npropellers
                println("  P$(ip): pitch F0 = $(round(F0_struct[ndof_wing_free+2*(ip-1)+1], digits=3)) N*m,  yaw F0 = $(round(F0_struct[ndof_wing_free+2*(ip-1)+2], digits=3)) N*m")
            end
            println("  |F0_wing| = $(round(norm(F0_struct[1:ndof_wing_free]), digits=2)),  |F0_prop| = $(round(norm(F0_struct[ndof_wing_free+1:end]), digits=2))")
            global printed_steady = true
        end
    end

    F_struct_n .= F_struct_final
    F_pert_n .= F_pert_final
    coupling_iterations[it] = iter_count
    coupling_disp_residual[it] = state_res
    coupling_load_residual[it] = force_res
    coupling_equilibrium_residual[it] = equilibrium_res
    coupling_converged[it] = converged

    for isurf in 1:nsurf
        if iwake[isurf] < nwake[isurf]
            iwake[isurf] += 1
        end
    end

    status = converged ? "converged" : "not converged"
    println("Step $it/$N_LAST (t=$(round(t[it+1], digits=4)) s) partitioned correction: $status in $iter_count iterations, state_res=$(round(state_res, sigdigits=4)), load_res=$(round(force_res, sigdigits=4)), equilibrium_res=$(round(equilibrium_res, sigdigits=4))")

    nrm = maximum(abs, U[it+1])
    if any(isnan, U[it+1]) || nrm > U_ABORT
        rng = max(1, it-400):it
        env = [abs(U[k][pitch_dof_indices[1]]) for k in rng]
        tt = [t[k] for k in rng]
        gi = findall(>(1e-12), env)
        if length(gi) > 5
            growth_sigma = (log(env[gi[end]]) - log(env[gi[1]])) / (tt[gi[end]] - tt[gi[1]])
            println("\n>>> ABORT at t=$(round(t[it+1], digits=4)) s, |U|=$(round(nrm, sigdigits=4)); growth sigma(P1 pitch) approx $(round(growth_sigma, digits=4)) /s")
        else
            println("\n>>> ABORT at t=$(round(t[it+1], digits=4)) s, |U|=$(round(nrm, sigdigits=4))")
        end
        global N_LAST = it
        break
    end

    if it in save
        surface_history[it] = [copy(s) for s in system.surfaces]
    end
end
println("Simulation finished.")
# ==============================================================================
# 6. VALIDATION OUTPUT
# ==============================================================================
plot_time_limit_s = parse(Float64, get(ENV, "CHANG_PLOT_END_TIME_S", "5.0"))
results = write_chang_results(
    displacement_history = U,
    time = t,
    time_steps = dt,
    last_step = N_LAST,
    wing_node_count = nnodes,
    degrees_of_freedom_per_node = ndof,
    number_of_propellers = Npropellers,
    number_of_blades = Nb_prop,
    propeller_eta,
    span_length = b,
    density = RHO,
    freestream_speed = Vinf,
    interaction_on = INTERACTION_ON,
    requested_end_time = t_end,
    coupling_iterations,
    coupling_state_residual = coupling_disp_residual,
    coupling_load_residual,
    coupling_equilibrium_residual,
    coupling_converged,
    output_directory = VERIFY_OUTPUT_DIR,
    output_label = VERIFY_LABEL,
    plot_results = PLOT_RESULTS,
    plot_time_limit = plot_time_limit_s,
)
