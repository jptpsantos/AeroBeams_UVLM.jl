module ChangAeroelastic

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

export chang_case_defaults, load_chang_configuration, ChangModel,
       build_chang_model, build_chang_workspace, run_chang,
       update_aero_geometry_for_state!, assemble_structural_aero_load!, aero_load_for_state!

# Includes define functions only. A case is constructed explicitly below.
include(joinpath(@__DIR__, "..", "chang_case.jl"))
include("chang_configuration.jl")
include("chang_model_parameters.jl")
include("chang_structural_model.jl")
include("chang_uvlm_coupling.jl")
include("chang_workspaces.jl")
include("chang_simulation.jl")
include("chang_postprocessing.jl")

"""One case's resolved settings, derived model data, and structural matrices.

Mutable aerodynamic state is allocated separately by build_chang_workspace.
Treat model data as read-only after construction; build a new model to change a case.
"""
struct ChangModel{C,P,S,A,F}
    config::C
    parameters::P
    structural::S
    aerodynamic_options::A
    near_field_force_function::F
end

"""Build a model with its own copy of the resolved configuration."""
function build_chang_model(config)
    config = deepcopy(config)
    validate_configuration(config)
    parameters = build_chang_model_parameters(config)
    structural = assemble_chang_structural_model(parameters)
    aerodynamic_options = chang_aerodynamic_options(config)
    force_function = config.simulation.near_field_force_model == :legacy_imperial_segments ?
        legacy_imperial_segment_forces! : near_field_forces!
    return ChangModel(config, parameters, structural, aerodynamic_options, force_function)
end

# Plotting dependencies are loaded only for runs requesting visual output.
function load_chang_visualization(output)
    if output.plot_results || output.animate_wake
        isdefined(@__MODULE__, :Plots) || (@eval using Plots)
    end
    if output.plot_results && !isdefined(@__MODULE__, :plot_chang_time_histories)
        include(joinpath(@__DIR__, "chang_plotting.jl"))
    end
    if output.animate_wake && !isdefined(@__MODULE__, :record_chang_animation_frame!)
        include(joinpath(@__DIR__, "chang_animation.jl"))
    end
    return nothing
end

"""
    run_chang(config=load_chang_configuration())

Build a fresh model/workspace, solve the response, and save the requested output.
Return (; config, model, workspace, solution, results) for inspection or audits.
The supplied configuration is copied so separate runs own their settings.
"""
function run_chang(config = load_chang_configuration())
    load_chang_visualization(config.output)
    # Optional plotting methods may have been defined during this call.
    return Base.invokelatest(run_chang_case, config)
end

function run_chang_case(config)
    # 1. Build and check the physical and aerodynamic models.
    model = build_chang_model(config)
    config = model.config
    parameters, structural = model.parameters, model.structural
    println("Aeroelastic near-field force model: $(config.simulation.near_field_force_model)")
    println("Propeller moment projection: $(config.simulation.propeller_moment_projection)")
    report_and_validate_structural_model(structural)
    println("Matrices after BCs. Total DOFs (free): $(structural.ndof_free)")
    println(
        "Finite-core radius: max($(config.aerodynamic.segment_core_factor) Δs, " *
        "$(config.aerodynamic.chord_core_factor) c)",
    )
    workspace = build_chang_workspace(model)
    wake = chang_wake_context(workspace; interaction_on = config.simulation.interaction_on)
    animation = chang_animation_options(workspace.system, workspace.iwake, config.output)
    excitation = chang_excitation_options(config, parameters)
    integration = chang_integration_options(config, structural, parameters)
    report_chang_solver_options(excitation, integration)
    output_directory, output_label = chang_output_paths(config.output)

    # 2. Converge each structural/aerodynamic step, then advance its wake once.
    aerodynamic_load = (snapshot, state, step) ->
        aero_load_for_state!(model, workspace, snapshot, state, step)
    solution = solve_chang_aeroelastic!(
        workspace.system, structural;
        time = parameters.t,
        time_steps = parameters.dt,
        freestream_history = workspace.fs_vec,
        aerodynamic_load, wake, excitation, integration, animation,
    )

    # 3. Extract the accepted response and write the selected output.
    results = write_chang_results(
        solution;
        wing = config.wing, propeller = config.propeller, simulation = config.simulation,
        time = parameters.t, time_steps = parameters.dt,
        wing_node_count = parameters.nnodes, dofs_per_node = parameters.ndof,
        density = parameters.ref.rho, visualization = config.output,
        output_directory, output_label,
    )
    return (; config, model, workspace, solution, results)
end

end # module
