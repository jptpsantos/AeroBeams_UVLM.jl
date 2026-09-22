module WingPropellerAeroelastic

using LinearAlgebra
using StaticArrays
using WingPropellerUVLM

export chang_case, xu_case, bohnisch_case
export build_structural_model, build_uvlm_system, make_uvlm_aerodynamic_coupling
export run_time_domain_analysis, run_uvlm_time_domain_analysis

include(joinpath(@__DIR__, "..", "cases", "chang_case.jl"))
include(joinpath(@__DIR__, "..", "cases", "xu_case.jl"))
include(joinpath(@__DIR__, "..", "cases", "bohnisch_case.jl"))
include(joinpath(@__DIR__, "general_uvlm_coupling.jl"))

"Return an optional solver setting while preserving compact legacy case files."
solver_option(options, name::Symbol, default) =
    hasproperty(options, name) ? getproperty(options, name) : default

"""Build the common structural matrices for any physical case definition."""
build_structural_model(case) = assemble_wing_propeller_structure(
    case.wing,
    case.propellers,
    case.operating_condition,
)

"""
    build_uvlm_system(case; time_steps, ...)

Construct the unchanged UVLM backend from the same physical case used by the
structural solver. Every propeller origin is obtained from `attachment_node`;
there is no independent aerodynamic span station or hub-to-node offset.
"""
function build_uvlm_system(
    case;
    time_steps = [1.0e-3],
    finite_core = (chord, segment) -> max(1.0e-3chord, 1.0e-3segment),
    maximum_wake_rows_wing::Integer = 36,
    maximum_wake_rows_propeller = 36,
    verbose::Bool = false,
)
    isempty(time_steps) && throw(ArgumentError("time_steps must not be empty"))
    wing = case.wing
    geometry = wing.geometry
    aero = case.aerodynamic
    propellers = case.propellers
    operating = case.operating_condition
    span_nodes = collect(Float64, wing.span_nodes)
    span = last(span_nodes) - first(span_nodes)
    span_fraction = (span_nodes .- first(span_nodes)) ./ span
    chord = geometry.root_chord .+
        (geometry.tip_chord - geometry.root_chord) .* span_fraction
    xle_distribution = geometry.xle_root .+
        (geometry.xle_tip - geometry.xle_root) .* span_fraction
    area = span * (geometry.root_chord + geometry.tip_chord) / 2
    reference_chord = (geometry.root_chord + geometry.tip_chord) / 2
    reference = Reference(
        area,
        reference_chord,
        span,
        SVector(
            geometry.xle_root + geometry.elastic_axis_fraction * geometry.root_chord,
            first(span_nodes),
            0.0,
        ),
        operating.freestream_speed,
        operating.air_density,
    )
    freestream = Freestream(
        operating.freestream_speed,
        operating.angle_of_attack,
        operating.sideslip,
        SVector(0.0, 0.0, 0.0),
    )
    return initialize_wing_propeller_uvlm_system(
        xle = [geometry.xle_root, geometry.xle_tip],
        yle = [first(span_nodes), last(span_nodes)],
        zle = [0.0, 0.0],
        chord_geo = [geometry.root_chord, geometry.tip_chord],
        theta_geo = [0.0, 0.0],
        phi_geo = [0.0, 0.0],
        ns_wing = aero.spanwise_panels,
        nc_wing = aero.chordwise_panels,
        mirror_wing = hasproperty(aero, :mirror) ? aero.mirror : false,
        spacing_s_wing = Uniform(),
        spacing_c_wing = Uniform(),
        R_prop = [propeller.radius for propeller in propellers],
        c_prop = [propeller.blade_chord for propeller in propellers],
        ns_prop = [propeller.radial_panels for propeller in propellers],
        nc_prop = [propeller.chordwise_panels for propeller in propellers],
        blade_twists_prop = [propeller.blade_twists for propeller in propellers],
        Nb_prop = [propeller.blades for propeller in propellers],
        Npropellers = length(propellers),
        span_nodes = span_nodes,
        prop_attach_nodes = [propeller.attachment_node for propeller in propellers],
        chord = chord,
        xle_distribution = xle_distribution,
        ref = reference,
        symmetric_wing = aero.symmetric,
        fs = freestream,
        dt = collect(Float64, time_steps),
        nnodes = length(span_nodes),
        prop_pivot_offset_from_ea_A = SVector(0.0, 0.0, 0.0),
        hub_center_prop_A = SVector(0.0, 0.0, 0.0),
        fcore = finite_core,
        elastic_axis_fraction = geometry.elastic_axis_fraction,
        maximum_wake_rows_wing = maximum_wake_rows_wing,
        maximum_wake_rows_propeller = maximum_wake_rows_propeller,
        verbose = verbose,
    )
end

"""
    run_time_domain_analysis(case; aerodynamic_load, time, time_steps, ...)

Run the common structural time integration/coupling machinery for a case. The
caller supplies the aerodynamic load callback and wake transaction because the
validated UVLM backend owns that mutable state. Configuration names never enter
the integration path.
"""
function run_time_domain_analysis(
    case;
    aerodynamic_load,
    time,
    time_steps,
    initial_displacement = nothing,
    initial_velocity = nothing,
    external_load = nothing,
    begin_step = (step, state, time) -> nothing,
    commit_step = (step, state, time, time_step) -> nothing,
)
    structural = build_structural_model(case)
    state_count = length(time)
    length(time_steps) == state_count - 1 || throw(DimensionMismatch(
        "time_steps must contain one fewer entry than time",
    ))
    displacement = isnothing(initial_displacement) ? zeros(structural.ndof_free) :
        collect(Float64, initial_displacement)
    velocity = isnothing(initial_velocity) ? zeros(structural.ndof_free) :
        collect(Float64, initial_velocity)
    acceleration = structural.M \ (
        zeros(structural.ndof_free) .- structural.C * velocity .-
        structural.K * displacement
    )
    displacement_history = [copy(displacement)]
    velocity_history = [copy(velocity)]
    acceleration_history = [copy(acceleration)]
    previous_load = zeros(structural.ndof_free)
    solver_options = case.solver_options
    parameters = structural_integration_parameters(
        solver_options.time_integrator;
        newmark_alpha = solver_option(solver_options, :newmark_alpha, 0.05),
        generalized_alpha_rho_infinity = solver_option(
            solver_options,
            :generalized_alpha_rho_infinity,
            1.0,
        ),
    )
    coupling_scheme = validate_aeroelastic_coupling_scheme(
        solver_options.coupling_scheme,
    )
    coupling = PartitionedCouplingOptions(
        maximum_iterations = solver_option(
            solver_options,
            :coupling_maximum_iterations,
            10,
        ),
        state_tolerance = solver_option(
            solver_options,
            :coupling_state_tolerance,
            1.0e-5,
        ),
        load_tolerance = solver_option(
            solver_options,
            :coupling_load_tolerance,
            1.0e-2,
        ),
        equilibrium_tolerance = solver_option(
            solver_options,
            :coupling_equilibrium_tolerance,
            1.0e-10,
        ),
        coupled_equilibrium_tolerance = solver_option(
            solver_options,
            :coupling_coupled_equilibrium_tolerance,
            1.0e-5,
        ),
        relaxation = solver_option(solver_options, :coupling_relaxation, 0.5),
    )

    for step in eachindex(time_steps)
        begin_step(step, displacement, time[step])
        load_n = isnothing(external_load) ? zeros(structural.ndof_free) :
            external_load(time[step], displacement)
        load_np1 = isnothing(external_load) ? zeros(structural.ndof_free) :
            external_load(time[step + 1], displacement)
        load_at_state = state -> aerodynamic_load(state, step, time[step + 1])
        solution = if coupling_scheme == :loose_explicit
            current_load = load_at_state(displacement)
            loose_explicit_aeroelastic_step(
                solver_options.time_integrator,
                structural.M,
                structural.C,
                structural.K,
                displacement,
                velocity,
                acceleration,
                current_load,
                previous_load,
                load_np1,
                load_n,
                time_steps[step],
                parameters;
                options = coupling,
            )
        else
            partitioned_aeroelastic_step(
                solver_options.time_integrator,
                structural.M,
                structural.C,
                structural.K,
                displacement,
                velocity,
                acceleration,
                previous_load,
                load_np1,
                load_n,
                time_steps[step],
                parameters,
                load_at_state;
                options = coupling,
            )
        end
        displacement = solution.displacement
        velocity = solution.velocity
        acceleration = solution.acceleration
        previous_load = solution.trial_load
        commit_step(step, displacement, time[step + 1], time_steps[step])
        push!(displacement_history, copy(displacement))
        push!(velocity_history, copy(velocity))
        push!(acceleration_history, copy(acceleration))
    end
    return (;
        structural,
        displacement_history,
        velocity_history,
        acceleration_history,
        integration_parameters = parameters,
        coupling_options = coupling,
    )
end

"""Run the common structural/coupling algorithm with the general UVLM adapter."""
function run_uvlm_time_domain_analysis(
    case;
    time,
    time_steps,
    interaction::Bool = true,
    finite_core = (chord, segment) -> max(1.0e-3chord, 1.0e-3segment),
    maximum_wake_rows_wing::Integer = 36,
    maximum_wake_rows_propeller = 36,
    kwargs...,
)
    structural = build_structural_model(case)
    uvlm = build_uvlm_system(
        case;
        time_steps,
        finite_core,
        maximum_wake_rows_wing,
        maximum_wake_rows_propeller,
    )
    aerodynamic = make_uvlm_aerodynamic_coupling(
        case,
        structural,
        uvlm,
        time_steps;
        interaction,
        finite_core,
    )
    results = run_time_domain_analysis(
        case;
        aerodynamic_load = aerodynamic.load,
        begin_step = aerodynamic.begin_step,
        commit_step = aerodynamic.commit_step,
        time,
        time_steps,
        kwargs...,
    )
    return merge(results, (; uvlm, aerodynamic))
end

end
