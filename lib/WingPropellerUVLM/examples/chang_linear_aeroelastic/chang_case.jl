# User-controlled configuration for the Chang UVLM case.
# Physical quantities use SI units unless the variable name says otherwise.

function environment_flag(name::AbstractString, default::Bool)
    raw_value = lowercase(strip(get(ENV, name, string(default))))
    raw_value in ("1", "true", "yes", "on") && return true
    raw_value in ("0", "false", "no", "off") && return false
    error("$name must be true/false, yes/no, on/off, or 1/0")
end

Base.@kwdef struct WingConfig
    root_chord_m::Float64 = 1.8
    tip_chord_m::Float64 = 1.8
    span_m::Float64 = 7.5
    spanwise_panels::Int = 30
    chordwise_panels::Int = 5
end

Base.@kwdef struct PropellerConfig
    radius_m::Float64 = 1.15
    chord_m::Float64 = 0.197
    blades::Int = 4
    radial_panels::Int = 10
    chordwise_panels::Int = 5
    rotation_rpm::Float64 = 1207.96
    trim_speed_mps::Float64 = 65.0
    attachment_eta::Vector{Float64} = [0.83]
end

Base.@kwdef struct SimulationConfig
    freestream_speed_mps::Float64 = 85.0
    angle_of_attack_deg::Float64 = 0.0
    sideslip_deg::Float64 = 0.0
    azimuth_step_deg::Float64 = 5.0
    end_time_s::Float64 = 2
    interaction_on::Bool = false
    near_field_force_model::Symbol = :legacy_imperial_segments
    impulse_propeller_indices::Vector{Int} = [1]
end

const WING_CONFIG = WingConfig()
const PROPELLER_CONFIG = PropellerConfig()
const SIMULATION_CONFIG = SimulationConfig(
    impulse_propeller_indices = [1], # Use [1, 2] to excite both propellers.
)

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
    all(0.0 .<= prop.attachment_eta .<= 1.0) || error("Propeller attachment eta values must be between 0 and 1")
    length(unique(prop.attachment_eta)) == length(prop.attachment_eta) ||
        error("Propeller attachment eta values must be unique")

    sim.freestream_speed_mps >= 0 || error("Freestream speed cannot be negative")
    sim.azimuth_step_deg > 0 || error("Azimuth step must be positive")
    sim.end_time_s > 0 || error("Simulation end time must be positive")
    sim.near_field_force_model in (:imperial, :legacy_imperial_segments) ||
        error("Near-field force model must be :imperial or :legacy_imperial_segments")
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

validate_configuration(WING_CONFIG, PROPELLER_CONFIG, SIMULATION_CONFIG)
