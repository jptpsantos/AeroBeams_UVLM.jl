# # Linear aeroelastic response of the Chang wing--propeller model
#
# This example couples a linear beam/pylon model to a free-wake UVLM model.
# Each physical step converges the structural state and aerodynamic load while
# holding the wake fixed; the accepted wake is then advanced exactly once.
# A mean pre-impulse aerodynamic load is removed so the saved motion represents
# the perturbation response to a smooth propeller-pitch impulse.
#
# The structural basis is `[span, chord, down]`, with nodal rotations in the
# same order. Detailed coordinate and load-transfer conventions are documented
# in the example README and `src/chang_uvlm_coupling.jl`.

import Pkg

# Use the WingPropellerUVLM project that owns this example.
Pkg.activate(normpath(joinpath(@__DIR__, "..", "..")))

using LinearAlgebra
using StaticArrays
using WingPropellerUVLM:
    Uniform,
    Freestream,
    Reference,
    RotationMatrix,
    initialize_bohnisch_uvlm_system,
    get_nodal_properties_chang,
    grid_to_surface_panels,
    copy_surfaces_to_previous!,
    propagate_system!,
    near_field_forces!,
    legacy_imperial_segment_forces!,
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

# ### Problem setup

## Case definition
# Edit chang_case.jl for routine geometry, mesh, and operating-point changes.
# CHANG_* environment variables provide temporary overrides for automated runs.
include(joinpath(@__DIR__, "chang_case.jl"))
include(joinpath(@__DIR__, "src", "chang_simulation.jl"))

const AIR_DENSITY = 1.225 # kg/m^3; model parameters use this to create `ref`.

# Choose the dimensional near-field loads and propeller moment projection.
const AEROELASTIC_NEAR_FIELD_FORCE_MODEL =
    SIMULATION_CONFIG.near_field_force_model
const AEROELASTIC_NEAR_FIELD_FORCE_FUNCTION =
    AEROELASTIC_NEAR_FIELD_FORCE_MODEL == :legacy_imperial_segments ?
        legacy_imperial_segment_forces! : near_field_forces!
const AEROELASTIC_PROPELLER_MOMENT_PROJECTION =
    SIMULATION_CONFIG.propeller_moment_projection

println("Aeroelastic near-field force model: $AEROELASTIC_NEAR_FIELD_FORCE_MODEL")
println("Propeller moment projection: $AEROELASTIC_PROPELLER_MOMENT_PROJECTION")

## Output and visualization
# Plotting is useful for an interactive run. Animation is opt-in because its
# retained wake geometry and GIF can be large; studies disable both.
visualization_options = chang_visualization_options()

(visualization_options.plot_results || visualization_options.animate_wake) &&
    (@eval using Plots)
visualization_options.plot_results &&
    include(joinpath(@__DIR__, "src", "chang_plotting.jl"))
visualization_options.animate_wake &&
    include(joinpath(@__DIR__, "src", "chang_animation.jl"))
include(joinpath(@__DIR__, "src", "chang_postprocessing.jl"))

output_directory = normpath(get(
    ENV,
    "CHANG_OUTPUT_DIR",
    joinpath(@__DIR__, "output"),
))
output_label = get(
    ENV,
    "CHANG_OUTPUT_LABEL",
    "chang_linear_$(AEROELASTIC_NEAR_FIELD_FORCE_MODEL)_uvlm",
)
mkpath(output_directory)

## Structural model
# Physical mass, inertia, stiffness, damping, and the time grid are defined in
# chang_model_parameters.jl. The assembly applies the clamped-root boundary
# condition and uses center-of-mass spatial inertia with parallel-axis terms.
include(joinpath(@__DIR__, "src", "chang_model_parameters.jl"))
include(joinpath(@__DIR__, "src", "chang_structural_model.jl"))

println("Assembling Chang structural matrices (Z-DOWN)...")
structural = assemble_chang_structural_model(inertia_reference = :center_of_mass)
structural_diagnostics = report_and_validate_structural_model(structural)
println("Matrices after BCs. Total DOFs (free): $(structural.ndof_free)")

# The example-specific aerodynamic adapter uses these reduced-system indices.
free_dofs = structural.free_dofs
ndof_wing_free = structural.ndof_wing_free

## Aerodynamic model
# The finite-core law uses the local 3-D vortex-segment length together with an
# optional chord-based floor. Both factors are collected in the aerodynamic
# options below and can be overridden for dedicated sensitivity studies.
aerodynamic_options = chang_aerodynamic_options(nc_wing, L_pylon)
FCORE = aerodynamic_options.finite_core
println(
    "Finite-core radius: max($(aerodynamic_options.segment_core_factor) Δs, " *
    "$(aerodynamic_options.chord_core_factor) c)",
)

# Aerodynamic locations are expressed in frame A. The physical hub and the
# load-reduction point are distinct so their moment arms remain explicit.
elastic_axis_fraction = aerodynamic_options.elastic_axis_fraction
prop_pivot_offset_from_ea_A = aerodynamic_options.propeller_pivot_offset_A
hub_center_prop_A = aerodynamic_options.physical_hub_center_A
hub_center_load_A = aerodynamic_options.load_center_A

# The adapter contains the detailed structural/aerodynamic coordinate mapping.
include(joinpath(@__DIR__, "src", "chang_uvlm_coupling.jl"))
uvlm = initialize_chang_uvlm(
    finite_core = FCORE,
    elastic_axis_fraction = elastic_axis_fraction,
    propeller_pivot_offset_A = prop_pivot_offset_from_ea_A,
    physical_hub_center_A = hub_center_prop_A,
    maximum_wake_rows_wing = aerodynamic_options.maximum_wake_rows_wing,
    maximum_wake_rows_propeller = aerodynamic_options.maximum_wake_rows_propeller,
)

# These workspaces are shared with the example-specific geometry/load adapter.
(;
    ratio_wing,
    grids_prop_ref,
    attach_node_y,
    ea_x_aero,
    prop_surface_indices,
    nsurf,
    surface_interaction_id,
    nwake,
    system,
    repeated_points,
    iwake,
    fs_vec,
    save,
    TF,
    surface_history,
    nodal_forces_wing,
    nodal_moments_wing,
    EA_nodes_wing,
    nodal_forces_prop,
    grids_prop_current,
    T_pivot_A_current,
) = uvlm
T_hub_A_current = Vector{SVector{3,Float64}}(undef, Npropellers)
T_load_A_current = Vector{SVector{3,Float64}}(undef, Npropellers)
pitch_axis_A_current = Vector{SVector{3,Float64}}(undef, Npropellers)
yaw_axis_A_current = Vector{SVector{3,Float64}}(undef, Npropellers)

## Excitation and coupling solver
excitation_options = chang_excitation_options(
    Ω,
    Npropellers,
    SIMULATION_CONFIG.impulse_propeller_indices,
)
integration_options = chang_integration_options(
    structural;
    reference = ref,
    freestream_speed = Vinf,
    reference_area = Sref,
    reference_chord = cref,
    propeller_radius = R_prop,
    wing_node_count = nnodes,
    dofs_per_node = ndof,
)
report_chang_solver_options(excitation_options, integration_options)

wake_context = (;
    surface_count = nsurf,
    repeated_points,
    maximum_rows = nwake,
    active_rows = iwake,
    interaction_ids = surface_interaction_id,
    interaction_on = INTERACTION_ON,
    saved_steps = save,
    surface_history,
)

record_animation_frame = if visualization_options.animate_wake
    (surface_frames, wake_frames, active_rows_frames, frame_times, frame_time) ->
        record_chang_animation_frame!(
            surface_frames,
            wake_frames,
            active_rows_frames,
            frame_times,
            system,
            iwake,
            frame_time,
        )
else
    (arguments...) -> nothing
end
animation_options = (;
    enabled = visualization_options.animate_wake,
    stride = visualization_options.animation_stride,
    record_frame = record_animation_frame,
)

# ### Problem solution
# The detailed loop lives in chang_simulation.jl. Its key transaction is:
# snapshot -> fixed-wake coupling trials -> accept or restore -> advance wake.
aerodynamic_load = (aerodynamic_snapshot, state, step) -> aero_load_for_state!(
    system,
    aerodynamic_snapshot,
    state,
    step;
    print_loads = false,
)
solution = solve_chang_aeroelastic!(
    system,
    structural;
    time = t,
    time_steps = dt,
    freestream_history = fs_vec,
    aerodynamic_load,
    wake = wake_context,
    excitation = excitation_options,
    integration = integration_options,
    animation = animation_options,
)

# ### Post-processing
# Export the accepted response and coupling diagnostics. Plotting and animation
# are optional views of the same accepted states.
results = write_chang_results(
    displacement_history = solution.displacement_history,
    time = t,
    time_steps = dt,
    last_step = solution.last_step,
    wing_node_count = nnodes,
    degrees_of_freedom_per_node = ndof,
    number_of_propellers = Npropellers,
    number_of_blades = Nb_prop,
    propeller_eta = propeller_eta,
    span_length = b,
    density = ref.rho,
    freestream_speed = Vinf,
    interaction_on = INTERACTION_ON,
    near_field_force_model = AEROELASTIC_NEAR_FIELD_FORCE_MODEL,
    propeller_moment_projection = AEROELASTIC_PROPELLER_MOMENT_PROJECTION,
    requested_end_time = t_end,
    coupling_iterations = solution.coupling_iterations,
    coupling_state_residual = solution.coupling_state_residual,
    coupling_load_residual = solution.coupling_load_residual,
    coupling_equilibrium_residual = solution.coupling_equilibrium_residual,
    coupling_converged = solution.coupling_converged,
    output_directory = output_directory,
    output_label = output_label,
    plot_results = visualization_options.plot_results,
    plot_time_limit = visualization_options.plot_time_limit_s,
)

if visualization_options.animate_wake
    wake_animation_path = joinpath(output_directory, output_label * "_wing_wake.gif")
    animate_chang_wing_wake(
        solution.animation_surface_history,
        solution.animation_wake_history,
        solution.animation_active_wake_rows_history,
        solution.animation_time_history;
        output_path = wake_animation_path,
        fps = visualization_options.animation_fps,
        axis_limits = ((-3.0, 5.0), (0.0, 8.0), (-4.0, 4.0)),
        tick_spacing = 1.0,
    )
    results = merge(results, (; wake_animation_path))
end
