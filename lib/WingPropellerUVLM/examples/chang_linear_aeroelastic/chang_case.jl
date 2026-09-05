# User-controlled configuration for the Chang UVLM case.
# Physical quantities use SI units unless the variable name says otherwise.

# Read one environment override while retaining the editable defaults below.
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

# Wing geometry and aerodynamic mesh:
# - root_chord_m and tip_chord_m define the chord distribution.
# - span_m is the full modeled span.
# - spanwise_panels and chordwise_panels define the wing UVLM grid.
Base.@kwdef struct WingConfig
    root_chord_m::Float64 = 1.8
    tip_chord_m::Float64 = 1.8
    span_m::Float64 = 7.5
    spanwise_panels::Int = environment_number(
        Int, "CHANG_WING_SPAN_PANELS", 20,
    )
    chordwise_panels::Int = environment_number(
        Int, "CHANG_WING_CHORD_PANELS", 5,
    )
end

# Propeller geometry, mesh, rotation, and wing attachment:
# - radius_m, chord_m, and blades define the rotor geometry.
# - radial_panels and chordwise_panels define each blade UVLM grid.
# - rotation_rpm is specified at trim_speed_mps and preserves its advance ratio.
# - attachment_eta gives each propeller span location from root (0) to tip (1).
Base.@kwdef struct PropellerConfig
    radius_m::Float64 = 1.15
    chord_m::Float64 = 0.197
    blades::Int = 4
    radial_panels::Int = environment_number(
        Int, "CHANG_PROP_RADIAL_PANELS", 5,
    )
    chordwise_panels::Int = environment_number(
        Int, "CHANG_PROP_CHORD_PANELS", 5,
    )
    rotation_rpm::Float64 = environment_number(
        Float64, "CHANG_ROTATION_RPM", 1217.6962,
    )
    trim_speed_mps::Float64 = environment_number(
        Float64, "CHANG_TRIM_SPEED_MPS", 65.0,
    )
    attachment_eta::Vector{Float64} = [0.83]
end

# Simulation controls:
# - Flow speed and angles define the incoming air.
# - azimuth_step_deg sets the aerodynamic time step.
# - end_time_s sets the simulation duration.
# - interaction_on enables interaction between aerodynamic surface groups.
# - near_field_force_model selects :imperial (corrected direct model) or
#   :legacy_imperial_segments (original-compatible model).
# - propeller_moment_projection selects :exact_virtual_work (instantaneous
#   axes) or :fixed_aero_axes (original small-angle axes).
# - impulse_propeller_indices selects which propellers are excited.
Base.@kwdef struct SimulationConfig
    freestream_speed_mps::Float64 = environment_number(
        Float64, "CHANG_SPEED_MPS", 80.0, "CHANG_FREESTREAM_SPEED_MPS",
    )
    angle_of_attack_deg::Float64 = environment_number(
        Float64, "CHANG_AOA_DEG", 3.0,
    )
    sideslip_deg::Float64 = environment_number(Float64, "CHANG_SIDESLIP_DEG", 0.0)
    azimuth_step_deg::Float64 = environment_number(
        Float64, "CHANG_AZIMUTH_STEP_DEG", 5.0,
    )
    end_time_s::Float64 = environment_number(
        Float64, "CHANG_END_TIME_S", 5.0,
    )
    interaction_on::Bool = environment_flag(
        "CHANG_INTERACTION", false,
    )
    near_field_force_model::Symbol = environment_symbol(
        "CHANG_NEAR_FIELD_FORCE_MODEL", :imperial,
    )
    propeller_moment_projection::Symbol = environment_symbol(
        "CHANG_PROP_MOMENT_PROJECTION", :exact_virtual_work,
    )
    impulse_propeller_indices::Vector{Int} = [1]
end

const WING_CONFIG = WingConfig()
const PROPELLER_CONFIG = PropellerConfig()

# Active case used by run_chang_linear_aeroelastic.jl. Edit the defaults above
# for routine work; the environment overrides are intended for automated sweeps.
const SIMULATION_CONFIG = SimulationConfig()

# Check geometry, mesh, model selections, and excitation indices before running.
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

# Stop immediately if any case input is invalid.
validate_configuration(WING_CONFIG, PROPELLER_CONFIG, SIMULATION_CONFIG)
