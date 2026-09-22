# General wing-propeller aeroelastic response.
# Edit the grouped inputs below; the reference case files do not need to be
# changed for ordinary discretization or run-control studies.

import Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..", "..")))

include(joinpath(@__DIR__, "src", "WingPropellerAeroelastic.jl"))
using .WingPropellerAeroelastic
using WingPropellerUVLM

function wing_propeller_case_defaults()
    # 1. Reference configuration
    configuration = :xu # :xu, :bohnisch, or :chang

    # 2. Wing structural and aerodynamic discretization
    wing = (;
        structural_elements = 20,
        spanwise_panels = 20,
        chordwise_panels = 5,
        # `nothing` retains the Xu/Bohnisch reference stiffness multiplier.
        inactive_stiffness_multiplier = nothing,
    )

    # 3. Propeller discretization, installation, and operating parameter
    propeller = (;
        radial_panels = 5,
        chordwise_panels = 5,
        # Root=0 and tip=1. `nothing` retains the reference station.
        # Direct attachment node = round(eta*structural_elements)+1.
        attachment_eta = nothing,
        # Bohnisch: fixed RPM and variable beta75.
        # Chang/Xu: constant advance ratio. `nothing` retains case values.
        rpm = nothing,
        advance_ratio = nothing,
        beta75_deg = nothing,
    )

    # 4. Flow condition and simulation time
    simulation = (;
        # `nothing` retains the corresponding reference-case value.
        freestream_speed_mps = nothing,
        air_density_kgpm3 = nothing,
        angle_of_attack_deg = nothing,
        sideslip_deg = nothing,
        # Unless time_step_s is provided, dt follows the fastest rotor.
        azimuth_step_deg = 5.0,
        time_step_s = nothing,
        end_time_s = 0.001, # Short smoke run; increase for response studies.
        interaction_on = true,
    )

    # 5. Aerodynamic finite core
    aerodynamic = (;
        # A numeric value selects a fixed radius in metres. `nothing` uses
        # max(segment_core_factor*segment_length, chord_core_factor*chord).
        core_radius_m = nothing,
        segment_core_factor = 0.025,
        chord_core_factor = 0.025,
    )

    # 6. Retained wake length
    wake = (;
        retained_revolutions = 0.5,
        maximum_rows_wing = nothing,
        maximum_rows_propeller = nothing,
    )

    # 7. Structural time integration
    integration = (;
        time_integrator = :newmark_beta, # or :generalized_alpha
        newmark_alpha = 0.05,
        generalized_alpha_rho_infinity = 1.0,
    )

    # 8. Partitioned aeroelastic coupling
    coupling = (;
        scheme = :loose_explicit, # or :implicit_predictor_corrector
        maximum_iterations = 10,
        state_tolerance = 1.0e-5,
        load_tolerance = 1.0e-2,
        equilibrium_tolerance = 1.0e-10,
        coupled_equilibrium_tolerance = 1.0e-5,
        relaxation = 0.5,
    )

    return (;
        configuration,
        wing,
        propeller,
        simulation,
        aerodynamic,
        wake,
        integration,
        coupling,
    )
end

function run_wing_propeller_case(config = wing_propeller_case_defaults())
    mesh = (;
        element_count = config.wing.structural_elements,
        wing_spanwise_panels = config.wing.spanwise_panels,
        wing_chordwise_panels = config.wing.chordwise_panels,
        propeller_radial_panels = config.propeller.radial_panels,
        propeller_chordwise_panels = config.propeller.chordwise_panels,
    )

    # The case defines its physical speed law. The runner supplies only its
    # physical parameter: J for Chang/Xu, or RPM and beta75 for Bohnisch.
    case = if config.configuration == :chang
        eta = something(config.propeller.attachment_eta, [0.42, 0.83])
        chang_case(;
            mesh...,
            propeller_attachment_eta = eta isa Real ? [eta] : eta,
            propeller_advance_ratio = config.propeller.advance_ratio,
        )
    elseif config.configuration == :xu
        eta = something(config.propeller.attachment_eta, 1.6 / 5.7)
        xu_case(;
            mesh...,
            propeller_attachment_eta = eta isa Real ? eta : only(eta),
            inactive_stiffness_multiplier = something(
                config.wing.inactive_stiffness_multiplier,
                1.0e4,
            ),
            propeller_advance_ratio = something(config.propeller.advance_ratio, 1.96),
            propeller_beta75_deg = something(config.propeller.beta75_deg, 35.0),
        )
    elseif config.configuration == :bohnisch
        eta = something(config.propeller.attachment_eta, 1.0)
        bohnisch_case(;
            mesh...,
            propeller_attachment_eta = eta isa Real ? eta : only(eta),
            inactive_stiffness_multiplier = something(
                config.wing.inactive_stiffness_multiplier,
                1.0e3,
            ),
            propeller_rpm = something(config.propeller.rpm, 2500.0),
            propeller_beta75_deg = something(config.propeller.beta75_deg, 0.0),
        )
    else
        throw(ArgumentError("configuration must be :xu, :bohnisch, or :chang"))
    end

    simulation = config.simulation
    operating_condition = merge(case.operating_condition, (;
        freestream_speed = something(
            simulation.freestream_speed_mps,
            case.operating_condition.freestream_speed,
        ),
        air_density = something(
            simulation.air_density_kgpm3,
            case.operating_condition.air_density,
        ),
        angle_of_attack = isnothing(simulation.angle_of_attack_deg) ?
            case.operating_condition.angle_of_attack : deg2rad(simulation.angle_of_attack_deg),
        sideslip = isnothing(simulation.sideslip_deg) ?
            case.operating_condition.sideslip : deg2rad(simulation.sideslip_deg),
    ))
    solver_options = (;
        time_integrator = config.integration.time_integrator,
        newmark_alpha = config.integration.newmark_alpha,
        generalized_alpha_rho_infinity =
            config.integration.generalized_alpha_rho_infinity,
        coupling_scheme = config.coupling.scheme,
        coupling_maximum_iterations = config.coupling.maximum_iterations,
        coupling_state_tolerance = config.coupling.state_tolerance,
        coupling_load_tolerance = config.coupling.load_tolerance,
        coupling_equilibrium_tolerance = config.coupling.equilibrium_tolerance,
        coupling_coupled_equilibrium_tolerance =
            config.coupling.coupled_equilibrium_tolerance,
        coupling_relaxation = config.coupling.relaxation,
    )
    case = merge(case, (; operating_condition, solver_options))

    structural = build_structural_model(case)
    azimuth_step = Float64(simulation.azimuth_step_deg)
    end_time = Float64(simulation.end_time_s)
    isfinite(azimuth_step) && azimuth_step > 0 || throw(ArgumentError(
        "simulation.azimuth_step_deg must be finite and positive",
    ))
    time_step = isnothing(simulation.time_step_s) ?
        deg2rad(azimuth_step) / maximum(abs, structural.angular_speeds) :
        Float64(simulation.time_step_s)
    all(isfinite, (end_time, time_step)) && end_time > 0 && time_step > 0 ||
        throw(ArgumentError(
            "simulation end time and time step must be positive",
        ))
    time = collect(0:floor(Int, end_time / time_step)) .* time_step
    if last(time) < end_time - 10eps(end_time)
        push!(time, end_time)
    else
        time[end] = end_time
    end
    time_steps = diff(time)

    core_radius = config.aerodynamic.core_radius_m
    finite_core = isnothing(core_radius) ?
        (chord, segment) -> max(
            config.aerodynamic.segment_core_factor * segment,
            config.aerodynamic.chord_core_factor * chord,
        ) :
        (chord, segment) -> Float64(core_radius)

    default_wake_rows = max(
        1,
        ceil(Int, config.wake.retained_revolutions * 360 / azimuth_step),
    )
    maximum_wake_rows_wing = something(
        config.wake.maximum_rows_wing,
        default_wake_rows,
    )
    maximum_wake_rows_propeller = something(
        config.wake.maximum_rows_propeller,
        default_wake_rows,
    )

    println("Running $(config.configuration) with direct attachment nodes ",
        getproperty.(case.propellers, :attachment_node))
    results = run_uvlm_time_domain_analysis(
        case;
        time,
        time_steps,
        interaction = simulation.interaction_on,
        finite_core,
        maximum_wake_rows_wing,
        maximum_wake_rows_propeller,
    )
    return (; config, case, time, time_steps, results)
end

# =============================================================================
# RUN THE SELECTED CASE
# =============================================================================
# julia --project=lib/WingPropellerUVLM lib/WingPropellerUVLM/examples/wing_propeller_aeroelastic/run_case.jl
config = wing_propeller_case_defaults()
case_run = run_wing_propeller_case(config)
results = case_run.results
