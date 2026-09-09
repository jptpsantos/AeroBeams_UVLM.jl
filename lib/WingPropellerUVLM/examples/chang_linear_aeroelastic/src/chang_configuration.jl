# Resolve all case settings once. Simulation and model functions consume the
# returned configuration and never consult ENV or chang_case.jl themselves.

function environment_value(name::AbstractString, default, aliases::AbstractString...; env = ENV)
    candidates = (name, aliases...)
    key = findfirst(candidate -> haskey(env, candidate), candidates)
    return isnothing(key) ? string(default) : env[candidates[key]]
end

function environment_flag(name, default::Bool, aliases::AbstractString...; env = ENV)
    raw = lowercase(strip(environment_value(name, default, aliases...; env)))
    raw in ("1", "true", "yes", "on") && return true
    raw in ("0", "false", "no", "off") && return false
    error("$name must be true/false, yes/no, on/off, or 1/0")
end

function environment_number(::Type{T}, name, default, aliases::AbstractString...; env = ENV) where {T<:Real}
    raw = environment_value(name, default, aliases...; env)
    value = tryparse(T, raw)
    isnothing(value) && error("$name must be a valid $T value, received '$raw'")
    return value
end

environment_symbol(name, default, aliases::AbstractString...; env = ENV) =
    Symbol(lowercase(strip(environment_value(name, default, aliases...; env))))

"""Apply typed overrides; numeric parsing stays stable when a default is written as an integer."""
function environment_overrides(defaults::NamedTuple, bindings::NamedTuple, env)
    overrides = map(keys(bindings)) do field
        T, names = first(bindings[field]), Base.tail(bindings[field])
        default = defaults[field]
        value = if T == Bool
            environment_flag(names[1], default, names[2:end]...; env)
        elseif T <: Real
            environment_number(T, names[1], default, names[2:end]...; env)
        elseif T == Symbol
            environment_symbol(names[1], default, names[2:end]...; env)
        else
            environment_value(names[1], default, names[2:end]...; env)
        end
        return value
    end
    return merge(defaults, NamedTuple{keys(bindings)}(overrides))
end

"""Rotor angular speed at the case airspeed, preserving the specified advance ratio."""
function chang_rotor_angular_speed(config)
    prop, sim = config.propeller, config.simulation
    inflow_ratio = prop.trim_speed_mps / (prop.rotation_rpm * 2 * pi / 60) / prop.radius_m
    return sim.freestream_speed_mps / prop.radius_m / inflow_ratio
end

"""
    load_chang_configuration(defaults=chang_case_defaults(); env=ENV)

Return one complete nested NamedTuple containing the active case. Resolve all
CHANG_* overrides and automatic values here; `env=Dict{String,String}()` uses
only the supplied defaults. Each call owns its vector settings. To change a
case programmatically, merge changes into defaults before calling this loader.
"""
function load_chang_configuration(defaults = chang_case_defaults(); env = ENV)
    defaults = deepcopy(defaults)
    wing = environment_overrides(defaults.wing, (
        spanwise_panels = (Int, "CHANG_WING_SPAN_PANELS"),
        chordwise_panels = (Int, "CHANG_WING_CHORD_PANELS"),
    ), env)
    propeller = environment_overrides(defaults.propeller, (
        radial_panels = (Int, "CHANG_PROP_RADIAL_PANELS"),
        chordwise_panels = (Int, "CHANG_PROP_CHORD_PANELS"),
        rotation_rpm = (Float64, "CHANG_ROTATION_RPM"),
        trim_speed_mps = (Float64, "CHANG_TRIM_SPEED_MPS"),
    ), env)
    simulation = environment_overrides(defaults.simulation, (
        freestream_speed_mps = (Float64, "CHANG_SPEED_MPS", "CHANG_FREESTREAM_SPEED_MPS"),
        angle_of_attack_deg = (Float64, "CHANG_AOA_DEG"),
        sideslip_deg = (Float64, "CHANG_SIDESLIP_DEG"),
        azimuth_step_deg = (Float64, "CHANG_AZIMUTH_STEP_DEG"),
        end_time_s = (Float64, "CHANG_END_TIME_S"),
        interaction_on = (Bool, "CHANG_INTERACTION"),
        near_field_force_model = (Symbol, "CHANG_NEAR_FIELD_FORCE_MODEL"),
        propeller_moment_projection = (Symbol, "CHANG_PROP_MOMENT_PROJECTION"),
    ), env)
    aerodynamic = environment_overrides(defaults.aerodynamic, (
        segment_core_factor = (Float64, "CHANG_FCORE_SEGMENT_FACTOR"),
        chord_core_factor = (Float64, "CHANG_FCORE_CHORD_FACTOR"),
        hub_load_arm_factor = (Float64, "CHANG_HUB_LOAD_ARM_FACTOR"),
    ), env)
    wing_rows = isnothing(defaults.wake.maximum_rows_wing) ?
        defaults.wake.wing_rows_per_chord_panel * wing.chordwise_panels :
        defaults.wake.maximum_rows_wing
    wake = environment_overrides(merge(defaults.wake, (; maximum_rows_wing = wing_rows)), (
        maximum_rows_wing = (Int, "CHANG_WAKE_ROWS_WING"),
        maximum_rows_propeller = (Int, "CHANG_WAKE_ROWS_PROPELLER"),
    ), env)
    excitation = environment_overrides(defaults.excitation, (
        trim_revolutions = (Float64, "CHANG_TRIM_REVOLUTIONS"),
        trim_average_revolutions = (Float64, "CHANG_TRIM_AVERAGE_REVOLUTIONS"),
        impulse_duration_s = (Float64, "CHANG_IMPULSE_DURATION_S"),
        impulse_magnitude_nm = (Float64, "CHANG_IMPULSE_MAGNITUDE"),
    ), env)
    revolution_period = 2pi / abs(chang_rotor_angular_speed((; propeller, simulation)))
    start = isnothing(excitation.impulse_start_s) ?
        excitation.trim_revolutions * revolution_period : excitation.impulse_start_s
    excitation = merge(excitation, (;
        impulse_start_s = environment_number(Float64, "CHANG_IMPULSE_START_S", start; env),
    ))
    integration = environment_overrides(defaults.integration, (
        rho_inf = (Float64, "CHANG_GA_RHO_INF"),
        state_norm_limit = (Float64, "CHANG_STATE_ABORT_NORM"),
        propeller_angle_limit_deg = (Float64, "CHANG_PROP_ANGLE_ABORT_DEG"),
    ), env)
    coupling = environment_overrides(defaults.coupling, (
        maximum_iterations = (Int, "CHANG_COUPLING_MAX_ITER"),
        state_tolerance = (Float64, "CHANG_COUPLING_TOL_U"),
        load_tolerance = (Float64, "CHANG_COUPLING_TOL_F"),
        equilibrium_tolerance = (Float64, "CHANG_COUPLING_TOL_EQ"),
        coupled_equilibrium_tolerance = (Float64, "CHANG_COUPLING_TOL_COUPLED_EQ"),
        relaxation = (Float64, "CHANG_COUPLING_RELAXATION"),
    ), env)
    directory = isabspath(defaults.output.directory) ? defaults.output.directory :
        joinpath(@__DIR__, "..", defaults.output.directory)
    label = isnothing(defaults.output.label) ?
        "chang_linear_$(simulation.near_field_force_model)_uvlm" : defaults.output.label
    output = environment_overrides(merge(defaults.output, (; directory, label)), (
        directory = (String, "CHANG_OUTPUT_DIR"),
        label = (String, "CHANG_OUTPUT_LABEL"),
        plot_results = (Bool, "CHANG_PLOT_RESULTS"),
        plot_time_limit_s = (Float64, "CHANG_PLOT_END_TIME_S"),
        animate_wake = (Bool, "CHANG_ANIMATE_WAKE"),
        animation_stride = (Int, "CHANG_ANIMATION_STRIDE"),
        animation_fps = (Int, "CHANG_ANIMATION_FPS"),
    ), env)
    output = merge(output, (; directory = abspath(output.directory)))
    reference_frequency = isnothing(defaults.structural.damping_reference_frequency_hz) ?
        defaults.structural.pitch_frequency_hz : defaults.structural.damping_reference_frequency_hz
    structural = merge(defaults.structural, (; damping_reference_frequency_hz = reference_frequency))
    config = (; wing, propeller, simulation, aerodynamic, wake, excitation,
              integration, coupling, output, structural)
    validate_configuration(config)
    return config
end

# Validate geometry, mesh, model, and excitation settings.
function validate_configuration(config)
    wing, prop, sim = config.wing, config.propeller, config.simulation
    wing.root_chord_m > 0 || error("Wing root chord must be positive")
    wing.tip_chord_m > 0 || error("Wing tip chord must be positive")
    wing.span_m > 0 || error("Wing span must be positive")
    wing.spanwise_panels > 0 || error("Wing spanwise panel count must be positive")
    wing.chordwise_panels > 0 || error("Wing chordwise panel count must be positive")

    prop.radius_m > 0 || error("Propeller radius must be positive")
    prop.chord_m > 0 || error("Propeller chord must be positive")
    prop.blades > 0 || error("Number of propeller blades must be positive")
    prop.radial_panels > 0 || error("Propeller radial panel count must be positive")
    prop.chordwise_panels > 0 || error("Propeller chordwise panel count must be positive")
    prop.rotation_rpm > 0 || error("Propeller rotation speed must be positive")
    prop.trim_speed_mps > 0 || error("Propeller trim reference speed must be positive")
    all(0.0 .<= prop.attachment_eta .<= 1.0) || error("Propeller attachment eta values must be between 0 and 1")
    length(unique(prop.attachment_eta)) == length(prop.attachment_eta) ||
        error("Propeller attachment eta values must be unique")

    sim.freestream_speed_mps > 0 || error("A positive airspeed is required for the fixed-advance-ratio time grid")
    sim.azimuth_step_deg > 0 || error("Azimuth step must be positive")
    sim.end_time_s > 0 || error("Simulation end time must be positive")
    sim.near_field_force_model in (:imperial, :legacy_imperial_segments) ||
        error("Near-field force model must be :imperial or :legacy_imperial_segments")
    sim.propeller_moment_projection in (:fixed_aero_axes, :exact_virtual_work) ||
        error(
            "Propeller moment projection must be :fixed_aero_axes or " *
            ":exact_virtual_work",
        )
    isempty(sim.impulse_propeller_indices) &&
        error("At least one impulse propeller index is required")
    all(
        (1 .<= sim.impulse_propeller_indices) .&
        (sim.impulse_propeller_indices .<= length(prop.attachment_eta)),
    ) ||
        error("Impulse propeller indices must be between 1 and $(length(prop.attachment_eta))")
    length(unique(sim.impulse_propeller_indices)) ==
        length(sim.impulse_propeller_indices) ||
        error("Impulse propeller indices must be unique")
    for (label, value) in (
        ("Air density", sim.air_density_kgpm3),
        ("Trim revolutions", config.excitation.trim_revolutions),
        ("Trim averaging revolutions", config.excitation.trim_average_revolutions),
        ("Impulse duration", config.excitation.impulse_duration_s),
        ("Animation stride", config.output.animation_stride),
        ("Animation FPS", config.output.animation_fps),
    )
        isfinite(value) && value > 0 || error("$label must be finite and positive")
    end
    aero = config.aerodynamic
    all(x -> isfinite(x) && x >= 0, (aero.segment_core_factor, aero.chord_core_factor)) ||
        error("Finite-core factors must be finite and nonnegative")
    max(aero.segment_core_factor, aero.chord_core_factor) > 0 ||
        error("A positive finite core is required for free-wake self induction")
    isfinite(aero.hub_load_arm_factor) && aero.hub_load_arm_factor >= 0 ||
        error("Hub load-arm factor must be finite and nonnegative")
    config.wake.maximum_rows_wing >= 0 || error("Wing wake rows must be nonnegative")
    config.wake.maximum_rows_propeller >= 0 || error("Propeller wake rows must be nonnegative")
    isfinite(config.excitation.impulse_start_s) && config.excitation.impulse_start_s >= 0 ||
        error("Impulse start time must be finite and nonnegative")
    isfinite(config.excitation.impulse_magnitude_nm) || error("Impulse magnitude must be finite")
    config.integration.state_norm_limit > 0 || error("State limit must be positive")
    config.integration.propeller_angle_limit_deg > 0 || error("Propeller-angle limit must be positive")
    generalized_alpha_parameters(config.integration.rho_inf)
    PartitionedCouplingOptions(; config.coupling...)
    return nothing
end
