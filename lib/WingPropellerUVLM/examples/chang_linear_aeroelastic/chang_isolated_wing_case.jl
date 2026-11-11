# Chang isolated-wing variant. Edit this file, then run run_chang_isolated_wing.jl.
# Units: metres, seconds, degrees, and N m. chang_case.jl is left unchanged.
# Shared numerical defaults are inherited from that case; the settings below
# configure the bare wing. No propeller or pylon is assembled or discretized.
# Existing CHANG_* environment overrides still apply; CHANG_TIME_STEP_S controls dt.

function chang_isolated_wing_case_defaults(base)
    wing = (
        root_chord_m=1.8,
        tip_chord_m=1.8,
        span_m=7.5,
        symmetric=false,
        spanwise_panels=20,# Also the number of structural beam elements.
        chordwise_panels=10
    )

    simulation = merge(
        base.simulation,
        (;
            air_density_kgpm3=1.225,
            freestream_speed_mps=50.0,
            angle_of_attack_deg=3.0,
            sideslip_deg=0.0,
            time_step_s=0.00025,# Independent of RPM and propeller azimuth.
            end_time_s=5.0,
            near_field_force_model=:imperial,
            interaction_on=false,
            impulse_propeller_indices=Int[]
        )
    )

    aerodynamic = merge(base.aerodynamic,(;
        core_radius_m=0.002758,# Same fixed radius as the original case.
        elastic_axis_fraction=0.30
    ))
    wake = merge(
        base.wake,
        (;
            maximum_rows_wing=nothing,
            wing_rows_per_chord_panel=10,
            maximum_rows_propeller=0,
            retained_revolutions_propeller=0.0
        )
    )

    # As in the original solver, hold structural perturbations at zero while
    # developing the aerodynamic wake, average the baseline load, then subtract
    # it from subsequent loads. This does not compute a deflected static trim.
    excitation = merge(
        base.excitation,
        (;
            impulse_start_s=0.5,
            trim_average_duration_s=0.1,
            impulse_duration_s=0.15,
            impulse_magnitude_nm=-100000.0,
            wing_impulse_eta=1.0,# Select the nearest free beam node; 1 = tip.
            wing_impulse_dof=5      # Moment: 4 = torsion, 5 = bending, 6 = in-plane.
        )
    )

    integration = merge(
        base.integration,
        (;
            time_integrator=:newmark_beta,
            newmark_alpha=0.05,
            rho_inf=1.0,
            state_norm_limit=1.0e3
        )
    )
    coupling = merge(base.coupling,(;
        scheme=:loose_explicit,# Also supports :implicit_predictor_corrector.
    ))
    structural = merge(base.structural,(;
        stiffness_damping_ratio=0.0,
        damping_reference_frequency_hz=7.97 # Explicit wing damping reference.
    ))
    output = merge(
        base.output,
        (;
            directory="output/isolated_wing",
            label="chang_isolated_wing_uvlm",
            plot_results=true,
            plot_time_limit_s=5.0,
            animate_wake=false
        )
    )

    # An empty attachment list removes every propeller/pylon mass, stiffness,
    # gyroscopic term, modal DOF, blade surface, and propeller wake.
    propeller = merge(base.propeller,(; attachment_eta=Float64[]))
    return (; wing,propeller,simulation,aerodynamic,wake,excitation,
        integration,coupling,output,structural)
end
