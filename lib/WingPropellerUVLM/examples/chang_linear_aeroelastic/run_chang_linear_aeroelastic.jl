# Chang wing–propeller aeroelastic response.
# Edit chang_case.jl for the physical case, numerical controls, and output.
# CHANG_* environment variables override those defaults.
# The five steps below prepare the models, solve the response, and save it.

import Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..", "..")))

using LinearAlgebra
using StaticArrays
using WingPropellerUVLM:
    Uniform, Freestream, Reference, RotationMatrix,
    initialize_bohnisch_uvlm_system, get_nodal_properties_chang,
    grid_to_surface_panels, generate_panel_grid_and_interpolate, linear_interpolate_1d,
    copy_surfaces_to_previous!, propagate_system!, advance_wake!, snapshot_uvlm, restore_uvlm!,
    near_field_forces!, legacy_imperial_segment_forces!, imperial_nodal_forces, imperial_nodal_positions,
    generalized_alpha_parameters, PartitionedCouplingOptions,
    partitioned_generalized_alpha_step, smooth_hann_pulse_load

# 1. Read the case and choose the output.
include(joinpath(@__DIR__, "chang_case.jl"))
include(joinpath(@__DIR__, "src", "chang_simulation.jl"))
include(joinpath(@__DIR__, "src", "chang_uvlm_coupling.jl"))
include(joinpath(@__DIR__, "src", "chang_postprocessing.jl"))

const AIR_DENSITY = SIMULATION_DEFAULTS.air_density_kgpm3
const AEROELASTIC_NEAR_FIELD_FORCE_MODEL = SIMULATION_CONFIG.near_field_force_model
const AEROELASTIC_PROPELLER_MOMENT_PROJECTION = SIMULATION_CONFIG.propeller_moment_projection
const AEROELASTIC_NEAR_FIELD_FORCE_FUNCTION =
    AEROELASTIC_NEAR_FIELD_FORCE_MODEL == :legacy_imperial_segments ?
        legacy_imperial_segment_forces! : near_field_forces!

println("Aeroelastic near-field force model: $AEROELASTIC_NEAR_FIELD_FORCE_MODEL")
println("Propeller moment projection: $AEROELASTIC_PROPELLER_MOMENT_PROJECTION")

visualization_options = chang_visualization_options()
output_directory, output_label = chang_output_paths(@__DIR__, AEROELASTIC_NEAR_FIELD_FORCE_MODEL)
if visualization_options.plot_results || visualization_options.animate_wake
    using Plots
end
if visualization_options.plot_results
    include(joinpath(@__DIR__, "src", "chang_plotting.jl"))
end
if visualization_options.animate_wake
    include(joinpath(@__DIR__, "src", "chang_animation.jl"))
end

# 2. Build and check the structural model.
# Parameters define geometry, physical properties, rotor speed Ω, and time
# points t and step sizes dt in seconds. Assembly applies the clamped root.
# M, C, K are the free-DOF mass, damping/gyroscopic, and stiffness matrices.
include(joinpath(@__DIR__, "src", "chang_model_parameters.jl"))
include(joinpath(@__DIR__, "src", "chang_structural_model.jl"))
println("Assembling Chang structural matrices (Z-DOWN)...")
structural = assemble_chang_structural_model(inertia_reference = :center_of_mass)
structural_diagnostics = report_and_validate_structural_model(structural)
println("Matrices after BCs. Total DOFs (free): $(structural.ndof_free)")

# 3. Build the aerodynamic model and wake.
# The core law uses full local chord c and the 3-D span/radial edge length Δs.
aerodynamic_options = chang_aerodynamic_options(nc_wing, L_pylon)
println(
    "Finite-core radius: max($(aerodynamic_options.segment_core_factor) Δs, " *
    "$(aerodynamic_options.chord_core_factor) c)",
)
uvlm = initialize_chang_uvlm(aerodynamic_options)

# Bind the geometry/load workspaces used by the adapter and validation scripts.
include(joinpath(@__DIR__, "src", "chang_workspaces.jl"))
wake_context = chang_wake_context(uvlm; interaction_on = INTERACTION_ON)
animation_options = chang_animation_options(system, iwake, visualization_options)

# 4. Configure the excitation and solve the coupled response.
# A pre-impulse mean load is subtracted to obtain the perturbation response.
excitation_options = chang_excitation_options(
    Ω, Npropellers, SIMULATION_CONFIG.impulse_propeller_indices,
)
integration_options = chang_integration_options(
    structural;
    reference = ref, freestream_speed = Vinf,
    reference_area = Sref, reference_chord = cref, propeller_radius = R_prop,
    wing_node_count = nnodes, dofs_per_node = ndof,
)
report_chang_solver_options(excitation_options, integration_options)

aerodynamic_load = (snapshot, state, step) ->
    aero_load_for_state!(system, snapshot, state, step; print_loads = false)

# Each step converges loads and motion with a fixed wake, then advances it once.
solution = solve_chang_aeroelastic!(
    system, structural;
    time = t,
    time_steps = dt,
    freestream_history = fs_vec,
    aerodynamic_load,
    wake = wake_context,
    excitation = excitation_options,
    integration = integration_options,
    animation = animation_options,
)

# 5. Save the history and diagnostics; create the requested plot/animation.
# solution contains the full accepted state histories; results contains the
# extracted wing/propeller responses and output paths.
results = write_chang_results(
    solution;
    wing = WING_CONFIG, propeller = PROPELLER_CONFIG, simulation = SIMULATION_CONFIG,
    time = t, time_steps = dt, wing_node_count = nnodes, dofs_per_node = ndof,
    density = ref.rho, visualization = visualization_options,
    output_directory, output_label,
)
