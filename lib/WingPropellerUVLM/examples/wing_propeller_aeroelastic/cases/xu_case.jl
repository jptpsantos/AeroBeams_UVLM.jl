"""
Return the baseline physical model from V. Q. Liu Xu, *Propeller-Wing
Whirl Flutter* (TU Delft, 2020), Tables 9.1, 9.3, and 9.4.
"""
function xu_case(;
    inactive_stiffness_multiplier::Real = 1.0e4,
    element_count::Integer = 20,
    wing_spanwise_panels::Integer = element_count,
    wing_chordwise_panels::Integer = 5,
    propeller_radial_panels::Integer = 5,
    propeller_chordwise_panels::Integer = 5,
    propeller_attachment_eta::Real = 1.6 / (11.4 / 2),
    propeller_advance_ratio::Real = 1.96,
    propeller_beta75_deg::Real = 35.0,
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
    isfinite(propeller_advance_ratio) && propeller_advance_ratio != 0 ||
        throw(ArgumentError("propeller_advance_ratio must be finite and nonzero"))
    isfinite(propeller_beta75_deg) || throw(ArgumentError(
        "propeller_beta75_deg must be finite",
    ))

    semispan = 11.4 / 2
    attachment_node = attachment_node_from_span_fraction(
        element_count,
        propeller_attachment_eta,
    )
    attachment_position = propeller_attachment_eta * semispan
    elements_inboard = attachment_node - 1
    elements_outboard = element_count - elements_inboard
    # Calculate the attachment index from eta and Ne, then anchor that node at
    # the requested physical station. The remaining elements are uniform on
    # each side, retaining exact aerodynamic/structural hub colocation.
    span_nodes = if elements_inboard == 0
        collect(range(0.0, semispan; length = element_count + 1))
    elseif elements_outboard == 0
        collect(range(0.0, semispan; length = element_count + 1))
    else
        vcat(
            collect(range(0.0, attachment_position; length = elements_inboard + 1)),
            collect(range(
                attachment_position,
                semispan;
                length = elements_outboard + 1,
            ))[2:end],
        )
    end
    root_chord = 1.25
    tip_chord = 0.8
    chord_at_nodes = root_chord .+
        (tip_chord - root_chord) .* (span_nodes ./ semispan)
    mass_per_length = 25.0
    EI = 7.0e5
    GJ = 2.0e5

    wing = (;
        span_nodes,
        stiffness = (;
            EA = inactive_stiffness_multiplier * EI / semispan^2,
            EI_in_plane = inactive_stiffness_multiplier * EI,
            EI_out_of_plane = EI,
            EI_coupling = 0.0,
            GJ,
        ),
        mass_model = :consistent_distributed,
        inertia = (;
            mass_per_length,
            # Xu specifies a sectional radius of gyration of 25% local chord.
            torsional_inertia_per_length =
                mass_per_length .* (0.25 .* chord_at_nodes) .^ 2,
            cg_offset_chord = 0.0,
        ),
        damping = (; stiffness_coefficient = 0.0),
        boundary_conditions = (;
            constrained_nodes = [1],
            constrained_local_dofs = collect(1:6),
        ),
        geometry = (;
            root_chord,
            tip_chord,
            xle_root = 0.0,
            xle_tip = (root_chord - tip_chord) / 2,
            elastic_axis_fraction = 0.5,
        ),
    )

    pitch_yaw_inertia = 36.65
    mounting_frequency = 7.0
    mounting_stiffness = pitch_yaw_inertia * (2pi * mounting_frequency)^2
    first_moment = -39.38
    radial_panels = propeller_radial_panels
    pitch_at_75_percent_radius = deg2rad(propeller_beta75_deg)
    radial_fraction = collect(range(0.0, 1.0; length = radial_panels + 1))
    # Xu specifies a fixed-pitch propeller and its blade angle at 0.75R. A
    # constant geometric pitch therefore gives tan(beta) * r = constant.
    blade_twists = [
        radius_fraction == 0 ? pi / 2 : atan(
            0.75 * tan(pitch_at_75_percent_radius) / radius_fraction,
        ) for radius_fraction in radial_fraction
    ]
    propellers = [(;
        attachment_node,
        radius = 1.524 / 2,
        blades = 3,
        blade_chord = 0.094,
        radial_panels,
        chordwise_panels = propeller_chordwise_panels,
        blade_twists,
        beta75_deg = Float64(propeller_beta75_deg),
        mass = 43.0,
        spin_inertia = 1.55,
        pitch_inertia = pitch_yaw_inertia,
        yaw_inertia = pitch_yaw_inertia,
        pitch_stiffness = mounting_stiffness,
        yaw_stiffness = mounting_stiffness,
        damping_ratio = 0.005,
        pitch_first_moment = first_moment,
        yaw_first_moment = first_moment,
        pitch_cross_inertia = pitch_yaw_inertia,
        yaw_cross_inertia = pitch_yaw_inertia,
        wing_pitch_first_moment = first_moment,
        wing_yaw_first_moment = first_moment,
        wing_pitch_inertia = pitch_yaw_inertia,
        wing_yaw_inertia = pitch_yaw_inertia,
        speed_model = :constant_advance_ratio,
        advance_ratio = Float64(propeller_advance_ratio),
    )]

    operating_condition = (;
        freestream_speed = 77.0,
        air_density = 0.962870,
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
