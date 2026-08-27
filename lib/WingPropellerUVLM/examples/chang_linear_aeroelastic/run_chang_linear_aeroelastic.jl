# Chang linear wing-propeller aeroelastic simulation.
# Structural DOFs per node: [span, chord, down, torsion, chord rotation, yaw].
# Each physical step solves the partitioned structural/UVLM iteration first,
# then advances the accepted wake exactly once.

import Pkg

# Use the local WingPropellerUVLM package.
Pkg.activate(normpath(joinpath(@__DIR__, "..", "..")))

using LinearAlgebra
using StaticArrays
using DelimitedFiles
using Statistics

include(joinpath(@__DIR__, "chang_case.jl"))

# User controls.
const PLOT_RESULTS = true
const ANIMATE_WAKE = true
const PLOTS_REQUIRED = PLOT_RESULTS || ANIMATE_WAKE
PLOTS_REQUIRED && (@eval using Plots)

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
    advance_wake!,
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
# 1. SUPPORT FILES AND OUTPUT CONTROLS
# ==============================================================================
PLOT_RESULTS && include(joinpath(@__DIR__, "chang_plotting.jl"))
ANIMATE_WAKE && include(joinpath(@__DIR__, "chang_animation.jl"))
include(joinpath(@__DIR__, "chang_postprocessing.jl"))

const AIR_DENSITY = 1.225 # kg/m^3; `ref.rho` is used by the UVLM.
const VERIFY_OUTPUT_DIR = normpath(get(
    ENV,
    "CHANG_OUTPUT_DIR",
    joinpath(@__DIR__, "output"),
))
const VERIFY_LABEL = get(ENV, "CHANG_OUTPUT_LABEL", "chang_linear_imperial_uvlm")
const WAKE_ANIMATION_STRIDE = parse(Int, get(ENV, "CHANG_ANIMATION_STRIDE", "5"))
const WAKE_ANIMATION_FPS = parse(Int, get(ENV, "CHANG_ANIMATION_FPS", "15"))
WAKE_ANIMATION_STRIDE > 0 || error("CHANG_ANIMATION_STRIDE must be positive")
WAKE_ANIMATION_FPS > 0 || error("CHANG_ANIMATION_FPS must be positive")
mkpath(VERIFY_OUTPUT_DIR)

# ==============================================================================
# 2. PHYSICAL AND STRUCTURAL MODEL
# ==============================================================================
include(joinpath(@__DIR__, "chang_model_parameters.jl"))
include(joinpath(@__DIR__, "chang_structural_model.jl"))

println("Assembling Chang structural matrices (Z-DOWN)...")

structural = assemble_chang_structural_model()
structural_diagnostics = chang_structural_diagnostics(structural)
println("Wing-only modal frequencies (Hz): $(round.(structural_diagnostics.wing_modal_frequencies_hz, digits=4))")
println(
    "Structural checks: min eig(M)=$(structural_diagnostics.minimum_mass_eigenvalue), " *
    "min eig(K)=$(structural_diagnostics.minimum_stiffness_eigenvalue), " *
    "symmetry(M/K)=($(structural_diagnostics.mass_symmetry_error), " *
    "$(structural_diagnostics.stiffness_symmetry_error))",
)

# Expose the matrices needed by the driver and coupling adapter.
(;
    Ks_W, Ms_W, Cs_W, Ms_P, Cs_P, Ks_P, Bs_P, Ds_P,
    Fs_W, Gs_W, Hs_W, attach_dofs_all,
    M_global, C_global, K_global, free_dofs,
    M, C, K, ndof_free, ndof_wing_free, ndof_prop_free,
) = structural
println("Matrices after BCs. Total DOFs (free): $ndof_free")

# Structural displacement, velocity, and acceleration histories.
U = Vector{Vector{Float64}}(undef, length(t))
Ud = similar(U)
Udd = similar(U)
U0 = zeros(ndof_free)
U0d = zeros(ndof_free)
# Initial acceleration from structural equilibrium.
U0dd = isempty(M) ? zeros(ndof_free) : M \ (zeros(ndof_free) - C*U0d - K*U0)
U[1] = U0
Ud[1] = U0d
Udd[1] = U0dd

# ==============================================================================
# 3. UVLM MODEL
# ==============================================================================
println("Initializing global UVLM system...")
# Vortex-core radius proportional to local segment width.
FCORE = (c, Δs) -> 0.5 * Δs

# Wing elastic axis, hub position, and modal load point.
elastic_axis_fraction = 0.30
prop_pivot_offset_from_ea_A = SVector(0.0, 0.0, 0.0)
hub_center_prop_A = SVector(-L_pylon, 0.0, 0.0)
# Test it equal to L_Pylon
hub_center_load_A = SVector(-0.5 * L_pylon, 0.0, 0.0)

# Build one UVLM system containing the wing, all blades, and their wakes.
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

# Expose aerodynamic state and reusable coupling workspaces.
(;
    ratio_wing, grids_prop_ref, grids_prop_initial_global,
    attach_node_y, ea_x_aero, prop_surface_indices, surfaces, nsurf,
    surface_interaction_id, nwake, system, repeated_points, iwake,
    fs_vec, save, TF, surface_history, nodal_forces_wing,
    nodal_moments_wing, EA_nodes_wing, nodal_forces_prop,
    grids_prop_current, T_pivot_A_current,
) = uvlm
T_hub_A_current   = Vector{SVector{3,Float64}}(undef, Npropellers)
T_load_A_current  = Vector{SVector{3,Float64}}(undef, Npropellers)

# Chang-specific geometry update and aerodynamic load transfer.
include(joinpath(@__DIR__, "chang_uvlm_coupling.jl"))

# Accepted states retained for optional wake animation.
animation_surface_history = Vector{typeof(system.surfaces)}()
animation_wake_history = Vector{typeof(system.wakes)}()
animation_active_wake_rows_history = Vector{Vector{Int}}()
animation_time_history = Float64[]
if ANIMATE_WAKE
    record_chang_animation_frame!(
        animation_surface_history,
        animation_wake_history,
        animation_active_wake_rows_history,
        animation_time_history,
        system,
        iwake,
        t[1],
    )
end

# ==============================================================================
# 4. TIME-INTEGRATION AND COUPLING CONTROLS
# ==============================================================================
# Select the propellers that receive the pitch impulse; for example `[1, 2]`.
impulse_propeller_indices = SIMULATION_CONFIG.impulse_propeller_indices
pitch_dof_indices = [
    ndof_wing_free + 2*(ip - 1) + 1
    for ip in impulse_propeller_indices
]
impulse_magnitude = parse(Float64, get(ENV, "CHANG_IMPULSE_MAGNITUDE", "1000.0"))
propeller_revolution_period = 2π / abs(Ω)
trim_revolutions = parse(Float64, get(ENV, "CHANG_TRIM_REVOLUTIONS", "10.0"))
trim_average_revolutions = parse(Float64, get(ENV, "CHANG_TRIM_AVERAGE_REVOLUTIONS", "1.0"))
default_impulse_start_time = trim_revolutions * propeller_revolution_period
impulse_start_time = parse(Float64, get(
    ENV,
    "CHANG_IMPULSE_START_S",
    string(default_impulse_start_time),
))
impulse_duration = parse(Float64, get(ENV, "CHANG_IMPULSE_DURATION_S", "0.15"))
impulse_start_time >= 0 || error("CHANG_IMPULSE_START_S must be nonnegative")
impulse_duration > 0 || error("CHANG_IMPULSE_DURATION_S must be positive")
trim_revolutions > 0 || error("CHANG_TRIM_REVOLUTIONS must be positive")
trim_average_revolutions > 0 || error("CHANG_TRIM_AVERAGE_REVOLUTIONS must be positive")
impulse_triggered_msg = false

# Mean periodic load used as the perturbation baseline after wake startup.
F0_struct = zeros(ndof_free)
F0_accumulator = zeros(ndof_free)
F0_sample_count = 0
have_F0 = false
printed_steady = false
U_ABORT = 1.0e3
N_LAST = length(dt)
trim_average_start_time = max(
    0.0,
    impulse_start_time - trim_average_revolutions * propeller_revolution_period,
)

const GA_RHO_INF = parse(Float64, get(ENV, "CHANG_GA_RHO_INF", "0.7"))
# Generalized-alpha parameters and fixed-point tolerances.
const GA_PARAMS = generalized_alpha_parameters(GA_RHO_INF)
const COUPLING_MAX_ITER = parse(Int, get(ENV, "CHANG_COUPLING_MAX_ITER", "10"))
const COUPLING_TOL_U = parse(Float64, get(ENV, "CHANG_COUPLING_TOL_U", "1.0e-5"))
const COUPLING_TOL_F = parse(Float64, get(ENV, "CHANG_COUPLING_TOL_F", "1.0e-2"))
const COUPLING_TOL_EQ = parse(Float64, get(ENV, "CHANG_COUPLING_TOL_EQ", "1.0e-10"))
const COUPLING_TOL_COUPLED_EQ = parse(Float64, get(
    ENV,
    "CHANG_COUPLING_TOL_COUPLED_EQ",
    "1.0e-4",
))
const COUPLING_RELAXATION = parse(Float64, get(ENV, "CHANG_COUPLING_RELAXATION", "0.5"))
const COUPLING_OPTIONS = PartitionedCouplingOptions(
    maximum_iterations = COUPLING_MAX_ITER,
    state_tolerance = COUPLING_TOL_U,
    load_tolerance = COUPLING_TOL_F,
    equilibrium_tolerance = COUPLING_TOL_EQ,
    coupled_equilibrium_tolerance = COUPLING_TOL_COUPLED_EQ,
    relaxation = COUPLING_RELAXATION,
)

# Scales for dimensionless displacement and load residuals.
const COUPLING_STATE_SCALE = ones(ndof_free)
const COUPLING_LOAD_SCALE = ones(ndof_free)
const REFERENCE_FORCE_SCALE = max(0.5 * ref.rho * Vinf^2 * Sref, 1.0)
const REFERENCE_MOMENT_SCALE = max(REFERENCE_FORCE_SCALE * cref, 1.0)
for inode in 1:(nnodes - 1)
    node_start = ndof * (inode - 1)
    COUPLING_STATE_SCALE[node_start .+ (1:3)] .= max(cref, eps(Float64))
    COUPLING_STATE_SCALE[node_start .+ (4:6)] .= 1.0
    COUPLING_LOAD_SCALE[node_start .+ (1:3)] .= REFERENCE_FORCE_SCALE
    COUPLING_LOAD_SCALE[node_start .+ (4:6)] .= REFERENCE_MOMENT_SCALE
end
propeller_moment_scale = max(
    0.5 * ref.rho * Vinf^2 * π * R_prop^2 * R_prop,
    1.0,
)
COUPLING_LOAD_SCALE[(ndof_wing_free + 1):end] .= propeller_moment_scale
println(
    "Partitioned generalized-alpha: rho_inf=$(GA_PARAMS.rho_inf), " *
    "alpha_m=$(GA_PARAMS.alpha_m), alpha_f=$(GA_PARAMS.alpha_f), " *
    "gamma=$(GA_PARAMS.gamma), beta=$(GA_PARAMS.beta)",
)
println(
    "Coupling correction: max_iter=$COUPLING_MAX_ITER, tol_u=$COUPLING_TOL_U, " *
    "tol_f=$COUPLING_TOL_F, tol_linear=$COUPLING_TOL_EQ, " *
    "tol_coupled=$COUPLING_TOL_COUPLED_EQ, relaxation=$COUPLING_RELAXATION",
)
println(
    "Trim baseline: average $trim_average_revolutions revolution(s), " *
    "from t=$(round(trim_average_start_time, digits=4)) s " *
    "to t=$(round(impulse_start_time, digits=4)) s",
)

coupling_iterations = fill(0, length(dt))
coupling_disp_residual = fill(NaN, length(dt))
coupling_load_residual = fill(NaN, length(dt))
coupling_equilibrium_residual = fill(NaN, length(dt))
coupling_converged = fill(false, length(dt))
# Accepted aerodynamic perturbation load at t[n].
F_pert_n = zeros(ndof_free)

# ==============================================================================
# 5. COUPLED TIME MARCHING
# ==============================================================================
println("Starting coupled aeroelastic simulation...")

for it = 1:length(dt)
    # A. Freeze and snapshot the accepted aerodynamic state at t[n].
    copy_surfaces_to_previous!(system, nsurf)
    snap = snapshot_uvlm(system)
    dt_i = dt[it]

    # B. Evaluate the external load at both integration endpoints.
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

    # C. Iterate structure and UVLM loads with the wake held at t[n].
    # The callback restores `snap`, updates geometry, solves circulation/loads,
    # and transfers the Imperial loads to structural DOFs.
    last_full_aerodynamic_load = zeros(ndof_free)
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
            last_full_aerodynamic_load .= aero_load_for_state!(
                system,
                snap,
                state_guess,
                it;
                print_loads = false,
            )
            return have_F0 ? last_full_aerodynamic_load .- F0_struct : zeros(ndof_free)
        end;
        options = COUPLING_OPTIONS,
        require_load_convergence = have_F0,
        state_scale = COUPLING_STATE_SCALE,
        load_scale = COUPLING_LOAD_SCALE,
    )

    # Read the final structural/aerodynamic fixed-point pair.
    U_corr = correction.displacement
    Ud_corr = correction.velocity
    Udd_corr = correction.acceleration
    F_pert_guess = correction.trial_load
    state_res = correction.state_residual
    force_res = correction.load_residual
    linear_equilibrium_res = correction.equilibrium_residual
    equilibrium_res = correction.coupled_equilibrium_residual
    converged = correction.converged
    iter_count = correction.iterations

    if !converged
        # Reject the physical step and restore its initial aerodynamic state.
        restore_uvlm!(system, snap)
        error(
            "Partitioned coupling failed at step $it (t=$(t[it + 1]) s) after " *
            "$iter_count iterations: state_res=$state_res, load_res=$force_res, " *
            "coupled_equilibrium_res=$equilibrium_res, " *
            "linear_equilibrium_res=$linear_equilibrium_res. The time step was not committed.",
        )
    end

    # D. Store the converged structural state and aerodynamic load.
    F_struct_final = copy(last_full_aerodynamic_load)
    F_pert_final = copy(F_pert_guess)

    U[it+1] = copy(U_corr)
    Ud[it+1] = copy(Ud_corr)
    Udd[it+1] = copy(Udd_corr)

    # E. Accumulate and activate the mean periodic trim load.
    time_np1 = t[it + 1]
    trim_sample = !have_F0 &&
        time_np1 >= trim_average_start_time - eps(time_np1) &&
        time_np1 <= impulse_start_time + eps(time_np1)
    if trim_sample
        F0_accumulator .+= F_struct_final
        global F0_sample_count += 1
    end
    if !have_F0 && time_np1 >= impulse_start_time - eps(time_np1)
        F0_sample_count > 0 || error("No aerodynamic samples were available for trim averaging")
        F0_struct .= F0_accumulator ./ F0_sample_count
        global have_F0 = true
        F_pert_final .= 0.0
        if !printed_steady
            println(
                "\n=== Mean trim generalized load captured at " *
                "t=$(round(time_np1, digits=4)) s from $F0_sample_count samples ===",
            )
            for ip in 1:Npropellers
                pitch_index = ndof_wing_free + 2*(ip - 1) + 1
                yaw_index = pitch_index + 1
                println(
                    "  P$ip: pitch F0 = $(round(F0_struct[pitch_index], digits=3)) N*m, " *
                    "yaw F0 = $(round(F0_struct[yaw_index], digits=3)) N*m",
                )
            end
            println(
                "  |F0_wing| = $(round(norm(F0_struct[1:ndof_wing_free]), digits=2)), " *
                "|F0_prop| = $(round(norm(F0_struct[ndof_wing_free+1:end]), digits=2))",
            )
            global printed_steady = true
        end
    end

    # Save the accepted load and convergence diagnostics.
    F_pert_n .= F_pert_final
    coupling_iterations[it] = iter_count
    coupling_disp_residual[it] = state_res
    coupling_load_residual[it] = force_res
    coupling_equilibrium_residual[it] = equilibrium_res
    coupling_converged[it] = converged

    # F. Convect and shed the accepted wake once. Do not also advance it inside
    # `aero_load_for_state!`.
    advance_wake!(
        system,
        fs_vec[it],
        dt_i;
        additional_velocity = nothing,
        repeated_points = repeated_points,
        nwake = iwake,
        interaction_id = surface_interaction_id,
        interaction = INTERACTION_ON,
    )

    # Activate the new wake row for the next physical step.
    for isurf in 1:nsurf
        if iwake[isurf] < nwake[isurf]
            iwake[isurf] += 1
        end
    end

    # G. Save optional animation data and report convergence.
    if ANIMATE_WAKE &&
        (it == 1 || it % WAKE_ANIMATION_STRIDE == 0 || it == length(dt))
        record_chang_animation_frame!(
            animation_surface_history,
            animation_wake_history,
            animation_active_wake_rows_history,
            animation_time_history,
            system,
            iwake,
            t[it + 1],
        )
    end

    println(
        "Step $it/$N_LAST (t=$(round(t[it+1], digits=4)) s): " *
        "$iter_count iterations, state_res=$(round(state_res, sigdigits=4)), " *
        "load_res=$(round(force_res, sigdigits=4)), " *
        "coupled_eq_res=$(round(equilibrium_res, sigdigits=4)), " *
        "linear_eq_res=$(round(linear_equilibrium_res, sigdigits=4))",
    )

    # Abort before a divergent solution contaminates later history.
    nrm = maximum(abs, U[it+1])
    if any(isnan, U[it+1]) || nrm > U_ABORT
        rng = max(1, it-400):it
        env = [abs(U[k][pitch_dof_indices[1]]) for k in rng]
        tt = [t[k] for k in rng]
        gi = findall(>(1e-12), env)
        if length(gi) > 5
            growth_sigma = (log(env[gi[end]]) - log(env[gi[1]])) / (tt[gi[end]] - tt[gi[1]])
            println(
                "\n>>> ABORT at t=$(round(t[it+1], digits=4)) s, " *
                "|U|=$(round(nrm, sigdigits=4)); growth sigma(P1 pitch) " *
                "approx $(round(growth_sigma, digits=4)) /s",
            )
        else
            println("\n>>> ABORT at t=$(round(t[it+1], digits=4)) s, |U|=$(round(nrm, sigdigits=4))")
        end
        global N_LAST = it
        break
    end

    # Store accepted surfaces only.
    if it in save
        surface_history[it] = [copy(s) for s in system.surfaces]
    end
end
println("Simulation finished.")

# Ensure the animation ends at the final accepted state.
if ANIMATE_WAKE && animation_time_history[end] != t[N_LAST + 1]
    record_chang_animation_frame!(
        animation_surface_history,
        animation_wake_history,
        animation_active_wake_rows_history,
        animation_time_history,
        system,
        iwake,
        t[N_LAST + 1],
    )
end

# ==============================================================================
# 6. VALIDATION OUTPUT
# ==============================================================================
# Write accepted histories, diagnostics, plots, and optional animation.
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
    propeller_eta = propeller_eta,
    span_length = b,
    density = ref.rho,
    freestream_speed = Vinf,
    interaction_on = INTERACTION_ON,
    requested_end_time = t_end,
    coupling_iterations = coupling_iterations,
    coupling_state_residual = coupling_disp_residual,
    coupling_load_residual = coupling_load_residual,
    coupling_equilibrium_residual = coupling_equilibrium_residual,
    coupling_converged = coupling_converged,
    output_directory = VERIFY_OUTPUT_DIR,
    output_label = VERIFY_LABEL,
    plot_results = PLOT_RESULTS,
    plot_time_limit = plot_time_limit_s,
)

if ANIMATE_WAKE
    wake_animation_path = joinpath(
        VERIFY_OUTPUT_DIR,
        VERIFY_LABEL * "_wing_wake.gif",
    )
    animate_chang_wing_wake(
        animation_surface_history,
        animation_wake_history,
        animation_active_wake_rows_history,
        animation_time_history;
        output_path = wake_animation_path,
        fps = WAKE_ANIMATION_FPS,
        axis_limits = ((-3.0, 5.0), (0.0, 8.0), (-4.0, 4.0)),
        tick_spacing = 1.0,
    )
    results = merge(results, (; wake_animation_path))
end
