const CHANG_EXAMPLE_DIRECTORY = normpath(joinpath(@__DIR__, "..", "..", "chang_linear_aeroelastic"))

if !isdefined(@__MODULE__, :ChangAeroelastic)
    include(joinpath(CHANG_EXAMPLE_DIRECTORY, "src", "ChangAeroelastic.jl"))
end

"""Map the validated Chang tables into the common physical case interface."""
function chang_case(
    defaults = ChangAeroelastic.chang_case_defaults();
    env = Dict{String,String}(),
    element_count::Integer = defaults.wing.spanwise_panels,
    wing_spanwise_panels::Integer = element_count,
    wing_chordwise_panels::Integer = defaults.wing.chordwise_panels,
    propeller_radial_panels::Integer = defaults.propeller.radial_panels,
    propeller_chordwise_panels::Integer = defaults.propeller.chordwise_panels,
    propeller_attachment_eta = defaults.propeller.attachment_eta,
    propeller_advance_ratio = nothing,
)
    element_count > 0 || throw(ArgumentError("element_count must be positive"))
    wing_spanwise_panels > 0 || throw(ArgumentError(
        "wing_spanwise_panels must be positive",
    ))
    wing_chordwise_panels > 0 || throw(ArgumentError(
        "wing_chordwise_panels must be positive",
    ))
    propeller_radial_panels > 0 || throw(ArgumentError(
        "propeller_radial_panels must be positive",
    ))
    propeller_chordwise_panels > 0 || throw(ArgumentError(
        "propeller_chordwise_panels must be positive",
    ))
    attachment_eta = propeller_attachment_eta isa Real ?
        [Float64(propeller_attachment_eta)] : collect(Float64, propeller_attachment_eta)
    all(eta -> isfinite(eta) && 0 <= eta <= 1, attachment_eta) ||
        throw(ArgumentError("every propeller attachment eta must lie in [0, 1]"))
    configured_defaults = merge(defaults, (;
        wing = merge(defaults.wing, (;
            spanwise_panels = Int(element_count),
            chordwise_panels = Int(wing_chordwise_panels),
        )),
        propeller = merge(defaults.propeller, (;
            radial_panels = Int(propeller_radial_panels),
            chordwise_panels = Int(propeller_chordwise_panels),
            attachment_eta,
        )),
    ))
    config = ChangAeroelastic.load_chang_configuration(configured_defaults; env)
    parameters = ChangAeroelastic.build_chang_model_parameters(config)
    distribution = parameters.wing_stiffness_distribution
    wing = (;
        span_nodes = parameters.span_nodes,
        stiffness = (;
            EA = (; eta = distribution.eta, values = distribution.EA),
            EI_in_plane = (; eta = distribution.eta, values = distribution.EIz),
            EI_out_of_plane = (; eta = distribution.eta, values = distribution.EIy),
            EI_coupling = (; eta = distribution.eta, values = distribution.EIzy),
            GJ = (; eta = distribution.eta, values = distribution.GJ),
        ),
        mass_model = :lumped_nodal,
        inertia = (; spatial_inertia_blocks = parameters.spatial_inertia_node_blocks),
        damping = (;
            global_stiffness_coefficient = 2parameters.stiffness_damping_ratio /
                parameters.stiffness_damping_reference_omega,
        ),
        boundary_conditions = (;
            constrained_nodes = [1],
            constrained_local_dofs = collect(1:6),
        ),
        geometry = (;
            root_chord = config.wing.root_chord_m,
            tip_chord = config.wing.tip_chord_m,
            xle_root = parameters.xle[1],
            xle_tip = parameters.xle[2],
            elastic_axis_fraction = config.aerodynamic.elastic_axis_fraction,
        ),
    )

    reference_advance_ratio = pi * parameters.Vinf / (parameters.Ω * parameters.R_prop)
    advance_ratio = something(propeller_advance_ratio, reference_advance_ratio)
    isfinite(advance_ratio) && advance_ratio != 0 || throw(ArgumentError(
        "propeller_advance_ratio must be finite and nonzero",
    ))
    propellers = [
        (;
            attachment_node = parameters.prop_attach_nodes[index],
            radius = parameters.R_prop,
            blades = parameters.Nb_prop,
            blade_chord = parameters.c_prop,
            radial_panels = parameters.ns_prop,
            chordwise_panels = parameters.nc_prop,
            blade_twists = copy(parameters.blade_twists_prop),
            mass = parameters.mP_prop,
            spin_inertia = parameters.Ix_prop,
            pitch_inertia = parameters.Inθ_prop,
            yaw_inertia = parameters.Inψ_prop,
            pitch_stiffness = parameters.Kθ_prop,
            yaw_stiffness = parameters.Kψ_prop,
            damping_ratio = parameters.ξ_prop,
            pitch_first_moment = parameters.SθP_prop,
            yaw_first_moment = parameters.SψP_prop,
            pitch_cross_inertia = parameters.IθαP_prop,
            yaw_cross_inertia = parameters.IψγP_prop,
            wing_pitch_first_moment = parameters.SαP_prop,
            wing_yaw_first_moment = parameters.SγP_prop,
            wing_pitch_inertia = parameters.IαP_prop,
            wing_yaw_inertia = parameters.IγP_prop,
            speed_model = :constant_advance_ratio,
            advance_ratio,
        )
        for index in 1:parameters.Npropellers
    ]
    operating_condition = (;
        freestream_speed = parameters.Vinf,
        air_density = parameters.ref.rho,
        angle_of_attack = parameters.alpha,
        sideslip = config.simulation.sideslip_deg * pi / 180,
    )
    aerodynamic = (;
        spanwise_panels = wing_spanwise_panels,
        chordwise_panels = wing_chordwise_panels,
        mirror = parameters.mirror_wing,
        symmetric = parameters.symmetric_wing,
    )
    solver_options = (;
        time_integrator = config.integration.time_integrator,
        coupling_scheme = config.coupling.scheme,
    )
    return (; wing, propellers, operating_condition, aerodynamic, solver_options,
        legacy_config = config, legacy_parameters = parameters)
end
