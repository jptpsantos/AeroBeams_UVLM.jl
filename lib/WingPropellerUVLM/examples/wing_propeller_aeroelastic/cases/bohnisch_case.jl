"""Return the supplied Bohnisch reference as physical inputs to the common solver."""
function bohnisch_case(;
    inactive_stiffness_multiplier::Real = 1.0e3,
    element_count::Integer = 20,
    wing_spanwise_panels::Integer = element_count,
    wing_chordwise_panels::Integer = 5,
    propeller_radial_panels::Integer = 5,
    propeller_chordwise_panels::Integer = 5,
    propeller_attachment_eta::Real = 1.0,
    propeller_rpm::Real = 2500.0,
    propeller_beta75_deg::Real = 0.0,
)
    inactive_stiffness_multiplier >= 1 || throw(ArgumentError(
        "inactive_stiffness_multiplier must be at least one",
    ))

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
    isfinite(propeller_rpm) && propeller_rpm > 0 || throw(ArgumentError(
        "propeller_rpm must be finite and positive",
    ))
    isfinite(propeller_beta75_deg) || throw(ArgumentError(
        "propeller_beta75_deg must be finite",
    ))

    root_chord = 1.25
    tip_chord = 1.25
    span = 5.7
    span_nodes = collect(range(0.0, span; length = element_count + 1))
    semichord = root_chord / 2
    elastic_axis = 2 * (0.35 - 0.5)
    center_of_mass = 2 * (0.46 - 0.5)
    cg_offset_chord = (center_of_mass - elastic_axis) * semichord
    EI = 2.09e6
    GJ = 2.7e5

    # The reduced reference contains out-of-plane bending and torsion. The full
    # beam's in-plane bending and axial coordinates are suppressed with a
    # configurable finite stiffness ratio, not an extreme sentinel value.
    wing = (;
        span_nodes,
        stiffness = (;
            EA = inactive_stiffness_multiplier * EI / span^2,
            EI_in_plane = inactive_stiffness_multiplier * EI,
            EI_out_of_plane = EI,
            EI_coupling = 0.0,
            GJ,
        ),
        mass_model = :consistent_distributed,
        inertia = (;
            mass_per_length = 18.0,
            torsional_inertia_per_length = 4.5,
            cg_offset_chord,
        ),
        damping = (; stiffness_coefficient = 0.001),
        boundary_conditions = (;
            constrained_nodes = [1],
            constrained_local_dofs = collect(1:6),
        ),
        geometry = (;
            root_chord,
            tip_chord,
            xle_root = 0.0,
            xle_tip = 0.0,
            elastic_axis_fraction = 0.35,
        ),
    )

    radius = 0.9
    spin_inertia = 2.38
    transverse_inertia = 4.85
    pitch_stiffness = 65000.0
    yaw_stiffness = 65000.0
    rotor_mass = 8.8
    motor_mass = 25.0
    rotor_cg = -0.4
    motor_cg = -0.3
    total_mass = rotor_mass + motor_mass
    pitch_first_moment = rotor_mass * rotor_cg + motor_mass * motor_cg
    propeller_to_ea_offset = -0.2
    ltheta = propeller_to_ea_offset + elastic_axis * semichord
    pitch_cross_inertia = transverse_inertia +
        (ltheta - elastic_axis * semichord) * pitch_first_moment
    wing_pitch_first_moment = pitch_first_moment +
        (ltheta - elastic_axis * semichord) * total_mass
    wing_pitch_inertia = transverse_inertia +
        (ltheta - elastic_axis * semichord)^2 * total_mass +
        2 * (ltheta - elastic_axis * semichord) * pitch_first_moment

    radial_panels = propeller_radial_panels
    _, _, blade_twists_degrees = get_nodal_properties(radial_panels)
    # The reference twist distribution is zero at 0.75R. Adding beta75 as a
    # collective offset therefore sets the requested blade angle at 0.75R.
    blade_twists_degrees .+= propeller_beta75_deg
    propellers = [(;
        attachment_node = attachment_node_from_span_fraction(
            element_count,
            propeller_attachment_eta,
        ),
        radius,
        blades = 4,
        blade_chord = get_chord_over_R(0.75) * radius,
        radial_panels,
        chordwise_panels = propeller_chordwise_panels,
        blade_twists = deg2rad.(blade_twists_degrees),
        beta75_deg = Float64(propeller_beta75_deg),
        mass = total_mass,
        spin_inertia,
        pitch_inertia = transverse_inertia,
        yaw_inertia = transverse_inertia,
        pitch_stiffness,
        yaw_stiffness,
        damping_ratio = 0.0,
        pitch_first_moment,
        yaw_first_moment = 0.0,
        pitch_cross_inertia,
        yaw_cross_inertia = 0.0,
        wing_pitch_first_moment,
        wing_yaw_first_moment = 0.0,
        wing_pitch_inertia,
        wing_yaw_inertia = 0.0,
        speed_model = :fixed_rpm,
        rpm = Float64(propeller_rpm),
    )]

    operating_condition = (;
        freestream_speed = 130.0,
        air_density = 1.225,
        angle_of_attack = 0.0,
        sideslip = 0.0,
    )
    aerodynamic = (;
        spanwise_panels = wing_spanwise_panels,
        chordwise_panels = wing_chordwise_panels,
        mirror = false,
        symmetric = false,
    )
    solver_options = (;
        time_integrator = :newmark_beta,
        coupling_scheme = :implicit_predictor_corrector,
    )
    return (; wing, propellers, operating_condition, aerodynamic, solver_options)
end
