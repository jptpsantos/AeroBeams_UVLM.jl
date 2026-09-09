# Chang wing–propeller case: edit the run settings in this file.
# Lengths are in metres, time in seconds, and angles in degrees.
# CHANG_* environment variables override these defaults when set (e.g. by sweeps).
# nothing selects an automatic value where explained below.

# Return a fresh set of editable defaults for one Chang case.
function chang_case_defaults()
    # 1. Wing geometry and mesh
    wing = (
        root_chord_m = 1.8,
        tip_chord_m = 1.8,
        span_m = 7.5,              # Modelled span, from the clamped root to the tip.

        spanwise_panels = 20,      # Also sets the number of structural beam elements.
        chordwise_panels = 5,
    )

    # 2. Propeller geometry, mesh, and installation
    propeller = (
        radius_m = 1.15,
        chord_m = 0.197,           # Constant blade chord.
        blades = 4,

        radial_panels = 5,         # Panels along each blade.
        chordwise_panels = 5,      # Panels across each blade chord.

        # RPM is specified at trim_speed_mps. The solver scales RPM with airspeed
        # to keep the advance ratio fixed.
        rotation_rpm = 1212.0,#1217.6962,
        trim_speed_mps = 65.0,

        # One entry per propeller: 0 = wing root, 1 = wing tip.
        attachment_eta = [0.83],
        collective_pitch_offset_deg = 0.0, # Added to the Chang blade-angle distribution.
    )

    # 3. Flow, time stepping, and coupling
    simulation = (
        air_density_kgpm3 = 1.225,
        freestream_speed_mps = 85.0,
        angle_of_attack_deg = 3.0,
        sideslip_deg = 0.0,

        azimuth_step_deg = 5.0,    # Rotor angle advanced per time step.
        end_time_s = 5.0,          # Total requested simulation duration.

        # true: include wing–propeller and propeller–propeller aerodynamic influence.
        # false: isolate those groups; blades within each propeller still interact.
        # Structural wing–propeller coupling remains active in both modes.
        interaction_on = false,

        # Aerodynamic loads: :imperial (corrected) or :legacy_imperial_segments.
        near_field_force_model = :imperial,

        # Propeller pitch/yaw moments: :exact_virtual_work or :fixed_aero_axes.
        propeller_moment_projection = :exact_virtual_work,

        impulse_propeller_indices = [1],   # Propellers receiving the pitch impulse.
    )

    # 4. Aerodynamic finite core and load geometry
    aerodynamic = (
        # Core radius = max(segment_core_factor * Δs, chord_core_factor * c).
        # Δs is the local span/radial edge length; c is the full local chord.
        segment_core_factor = 1e-6,
        chord_core_factor = 1e-6,

        elastic_axis_fraction = 0.30,      # Chord fraction measured from the leading edge.
        hub_load_arm_factor = 0.5,         # Modal load arm divided by pylon length.
    )

    # 5. Retained wake length
    wake = (
        # One wake row is shed per time step: retained wake age ≈ rows * Δt.
        # Wing: use an integer for a fixed row count, or nothing for automatic sizing.
        maximum_rows_wing = nothing,
        wing_rows_per_chord_panel = 10,    # Automatic count = this * active chordwise panels.
        maximum_rows_propeller = 72,
    )

    # 6. Trim baseline and pitch impulse
    excitation = (
        trim_revolutions = 10.0,           # Revolutions before the automatic impulse start.
        trim_average_revolutions = 1.0,    # Revolutions used for the mean baseline load.
        impulse_start_s = nothing,        # nothing: start after trim_revolutions.
        impulse_duration_s = 0.15,
        impulse_magnitude_nm = 1200.0,     # Peak pitch moment per selected propeller.
    )

    # 7. Time integration and stop limits
    integration = (
        # Generalized-alpha high-frequency spectral radius, between 0 and 1.
        # 1 = no algorithmic damping; smaller values increase numerical damping.
        rho_inf = 1.0,
        # The actual time step comes from simulation.azimuth_step_deg and RPM.
        state_norm_limit = 1.0e3,          # Maximum absolute state component (mixed m/rad).
        propeller_angle_limit_deg = Inf,  # Inf disables the propeller-angle stop limit.
    )

    # 8. Iterations coupling the aerodynamic loads and structural motion
    coupling = (
        maximum_iterations = 10,
        state_tolerance = 1.0e-5,               # Scaled change in displacement.
        load_tolerance = 1.0e-2,                # Scaled change in aerodynamic load.
        equilibrium_tolerance = 1.0e-10,        # Structural linear-solve residual.
        coupled_equilibrium_tolerance = 1.0e-4, # Equilibrium with the updated aero load.
        relaxation = 1.0,                      # 1 = full update; smaller values under-relax.
    )

    # 9. Files, plots, and wake animation
    output = (
        directory = "output",             # Relative to this example, or an absolute path.
        label = nothing,                   # nothing: chang_linear_<force_model>_uvlm.
        plot_results = true,
        plot_time_limit_s = 5.0,
        animate_wake = false,
        animation_stride = 5,              # Record every N accepted steps.
        animation_fps = 15,
        animation_axis_limits = ((-3.0, 5.0), (0.0, 8.0), (-4.0, 4.0)),
        animation_tick_spacing_m = 1.0,
    )

    # 10. Structural damping and pylon/rotor properties
    # The tabulated Chang wing mass/stiffness and blade distributions remain reference
    # data in the model helpers. The scalar controls for this run are collected here.
    structural = (
        propeller_damping_ratio = 0.0,
        stiffness_damping_ratio = 0.0,     # Global stiffness-proportional Rayleigh damping.
        damping_reference_frequency_hz = nothing, # nothing: use the pitch frequency below.

        pitch_frequency_hz = 7.97,
        yaw_frequency_hz = 7.97,
        pitch_stiffness_nm_per_rad = 19220.0,
        yaw_stiffness_nm_per_rad = 18916.0,
        twist_frequency_hz = 12.73,        # Frequency/stiffness pair defines rotor axial inertia.
        twist_stiffness_nm_per_rad = 16835.0,

        pylon_length_m = 5.6 * 0.3048,     # Original Chang value: 5.6 ft.
        pylon_mass_per_length_kgpm = 0.0506 * 14.5939029372064 / 0.3048, # 0.0506 slug/ft.
        blade_mass_kg = 1.44,
    )

    return (; wing, propeller, simulation, aerodynamic, wake, excitation, integration, coupling, output, structural)
end
