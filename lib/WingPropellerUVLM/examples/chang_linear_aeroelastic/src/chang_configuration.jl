# Configuration types, environment overrides, and validation for the Chang case.
# Included by chang_case.jl after its editable *_DEFAULTS groups are defined.

function environment_value(name::AbstractString, default, aliases::AbstractString...)
    key = findfirst(candidate -> haskey(ENV, candidate), (name, aliases...))
    return isnothing(key) ? string(default) : ENV[(name, aliases...)[key]]
end

function environment_flag(name::AbstractString, default::Bool, aliases::AbstractString...)
    raw_value = lowercase(strip(environment_value(name, default, aliases...)))
    raw_value in ("1", "true", "yes", "on") && return true
    raw_value in ("0", "false", "no", "off") && return false
    error("$name must be true/false, yes/no, on/off, or 1/0")
end

function environment_number(::Type{T}, name, default, aliases::AbstractString...) where {T<:Real}
    raw_value = environment_value(name, default, aliases...)
    value = tryparse(T, raw_value)
    isnothing(value) && error("$name must be a valid $T value, received '$raw_value'")
    return value
end

environment_symbol(name, default, aliases::AbstractString...) =
    Symbol(lowercase(strip(environment_value(name, default, aliases...))))

# Defaults come from chang_case.jl, so each editable value has a single source.
# Explicit constructor keywords continue to take precedence over these defaults.
Base.@kwdef struct WingConfig
    root_chord_m::Float64 = WING_DEFAULTS.root_chord_m
    tip_chord_m::Float64 = WING_DEFAULTS.tip_chord_m
    span_m::Float64 = WING_DEFAULTS.span_m
    spanwise_panels::Int = environment_number(
        Int, "CHANG_WING_SPAN_PANELS", WING_DEFAULTS.spanwise_panels,
    )
    chordwise_panels::Int = environment_number(
        Int, "CHANG_WING_CHORD_PANELS", WING_DEFAULTS.chordwise_panels,
    )
end

Base.@kwdef struct PropellerConfig
    radius_m::Float64 = PROPELLER_DEFAULTS.radius_m
    chord_m::Float64 = PROPELLER_DEFAULTS.chord_m
    blades::Int = PROPELLER_DEFAULTS.blades
    radial_panels::Int = environment_number(
        Int, "CHANG_PROP_RADIAL_PANELS", PROPELLER_DEFAULTS.radial_panels,
    )
    chordwise_panels::Int = environment_number(
        Int, "CHANG_PROP_CHORD_PANELS", PROPELLER_DEFAULTS.chordwise_panels,
    )
    rotation_rpm::Float64 = environment_number(
        Float64, "CHANG_ROTATION_RPM", PROPELLER_DEFAULTS.rotation_rpm,
    )
    trim_speed_mps::Float64 = environment_number(
        Float64, "CHANG_TRIM_SPEED_MPS", PROPELLER_DEFAULTS.trim_speed_mps,
    )
    attachment_eta::Vector{Float64} = copy(PROPELLER_DEFAULTS.attachment_eta)
end

Base.@kwdef struct SimulationConfig
    freestream_speed_mps::Float64 = environment_number(
        Float64, "CHANG_SPEED_MPS", SIMULATION_DEFAULTS.freestream_speed_mps,
        "CHANG_FREESTREAM_SPEED_MPS",
    )
    angle_of_attack_deg::Float64 = environment_number(
        Float64, "CHANG_AOA_DEG", SIMULATION_DEFAULTS.angle_of_attack_deg,
    )
    sideslip_deg::Float64 = environment_number(
        Float64, "CHANG_SIDESLIP_DEG", SIMULATION_DEFAULTS.sideslip_deg,
    )
    azimuth_step_deg::Float64 = environment_number(
        Float64, "CHANG_AZIMUTH_STEP_DEG", SIMULATION_DEFAULTS.azimuth_step_deg,
    )
    end_time_s::Float64 = environment_number(
        Float64, "CHANG_END_TIME_S", SIMULATION_DEFAULTS.end_time_s,
    )
    interaction_on::Bool = environment_flag(
        "CHANG_INTERACTION", SIMULATION_DEFAULTS.interaction_on,
    )
    near_field_force_model::Symbol = environment_symbol(
        "CHANG_NEAR_FIELD_FORCE_MODEL", SIMULATION_DEFAULTS.near_field_force_model,
    )
    propeller_moment_projection::Symbol = environment_symbol(
        "CHANG_PROP_MOMENT_PROJECTION", SIMULATION_DEFAULTS.propeller_moment_projection,
    )
    impulse_propeller_indices::Vector{Int} = copy(SIMULATION_DEFAULTS.impulse_propeller_indices)
end

# Validate geometry, mesh, model, and excitation settings.
function validate_configuration(wing::WingConfig, prop::PropellerConfig, sim::SimulationConfig)
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

    sim.freestream_speed_mps >= 0 || error("Freestream speed cannot be negative")
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
    return nothing
end
