# Runtime-controlled configuration used only by the Chang validation harness.

diagnostic_float(name, default) = parse(Float64, get(ENV, name, string(default)))
diagnostic_int(name, default) = parse(Int, get(ENV, name, string(default)))

function diagnostic_bool(name, default)
    value = lowercase(strip(get(ENV, name, string(default))))
    value in ("1", "true", "yes", "on") && return true
    value in ("0", "false", "no", "off") && return false
    error("$name must be true/false, yes/no, on/off, or 1/0")
end

diagnostic_symbol(name, default) = Symbol(lowercase(strip(get(ENV, name, string(default)))))

Base.@kwdef struct WingConfig
    root_chord_m::Float64 = 1.8
    tip_chord_m::Float64 = 1.8
    span_m::Float64 = 7.5
    spanwise_panels::Int = diagnostic_int("CHANG_TEST_WING_SPAN_PANELS", 30)
    chordwise_panels::Int = diagnostic_int("CHANG_TEST_WING_CHORD_PANELS", 5)
end

Base.@kwdef struct PropellerConfig
    radius_m::Float64 = 1.15
    chord_m::Float64 = 0.197
    blades::Int = 4
    radial_panels::Int = diagnostic_int("CHANG_TEST_PROP_RADIAL_PANELS", 5)
    chordwise_panels::Int = diagnostic_int("CHANG_TEST_PROP_CHORD_PANELS", 5)
    rotation_rpm::Float64 = diagnostic_float("CHANG_TEST_TRIM_RPM", 1207.96)
    trim_speed_mps::Float64 = diagnostic_float("CHANG_TEST_TRIM_SPEED_MPS", 65.0)
    attachment_eta::Vector{Float64} = [0.83]
end

Base.@kwdef struct SimulationConfig
    freestream_speed_mps::Float64 = diagnostic_float("CHANG_TEST_SPEED_MPS", 80.0)
    angle_of_attack_deg::Float64 = diagnostic_float("CHANG_TEST_AOA_DEG", 0.0)
    sideslip_deg::Float64 = 0.0
    azimuth_step_deg::Float64 = diagnostic_float("CHANG_TEST_AZIMUTH_STEP_DEG", 5.0)
    end_time_s::Float64 = diagnostic_float("CHANG_TEST_END_TIME_S", 2.0)
    interaction_on::Bool = diagnostic_bool("CHANG_TEST_INTERACTION", false)
    near_field_force_model::Symbol = diagnostic_symbol(
        "CHANG_TEST_BASE_FORCE_MODEL",
        :imperial,
    )
    impulse_propeller_indices::Vector{Int} = [1]
end

const WING_CONFIG = WingConfig()
const PROPELLER_CONFIG = PropellerConfig()
const SIMULATION_CONFIG = SimulationConfig()

function validate_case_configuration(
    wing::WingConfig,
    prop::PropellerConfig,
    sim::SimulationConfig,
)
    wing.root_chord_m > 0 || error("Wing root chord must be positive")
    wing.tip_chord_m > 0 || error("Wing tip chord must be positive")
    wing.span_m > 0 || error("Wing span must be positive")
    wing.spanwise_panels > 0 || error("Wing spanwise panel count must be positive")
    wing.chordwise_panels > 0 || error("Wing chordwise panel count must be positive")
    prop.radius_m > 0 || error("Propeller radius must be positive")
    prop.chord_m > 0 || error("Propeller chord must be positive")
    prop.blades > 0 || error("Propeller blade count must be positive")
    prop.radial_panels > 0 || error("Propeller radial panel count must be positive")
    prop.chordwise_panels > 0 || error("Propeller chordwise panel count must be positive")
    prop.rotation_rpm > 0 || error("Propeller rotation speed must be positive")
    prop.trim_speed_mps > 0 || error("Trim speed must be positive")
    all(0.0 .<= prop.attachment_eta .<= 1.0) ||
        error("Propeller attachment eta values must be between 0 and 1")
    sim.freestream_speed_mps >= 0 || error("Freestream speed cannot be negative")
    sim.azimuth_step_deg > 0 || error("Azimuth step must be positive")
    sim.end_time_s > 0 || error("Simulation end time must be positive")
    sim.near_field_force_model in (:imperial, :legacy_imperial_segments) ||
        error("Near-field force model must be :imperial or :legacy_imperial_segments")
    return nothing
end

validate_case_configuration(WING_CONFIG, PROPELLER_CONFIG, SIMULATION_CONFIG)

const TEST_WAKE_ROWS_WING = diagnostic_int("CHANG_TEST_WAKE_ROWS_WING", 50)
const TEST_WAKE_ROWS_PROPELLER = diagnostic_int("CHANG_TEST_WAKE_ROWS_PROP", 72)
const TEST_HUB_LOAD_ARM_FACTOR = diagnostic_float("CHANG_TEST_HUB_LOAD_ARM_FACTOR", 0.5)
const TEST_FCORE_SEGMENT_FACTOR = diagnostic_float("CHANG_TEST_FCORE_SEGMENT_FACTOR", 0.5)
const TEST_FCORE_CHORD_FACTOR = diagnostic_float("CHANG_TEST_FCORE_CHORD_FACTOR", 0.0)
TEST_WAKE_ROWS_WING >= 0 || error("CHANG_TEST_WAKE_ROWS_WING must be nonnegative")
TEST_WAKE_ROWS_PROPELLER >= 0 || error("CHANG_TEST_WAKE_ROWS_PROP must be nonnegative")
TEST_HUB_LOAD_ARM_FACTOR >= 0 || error("CHANG_TEST_HUB_LOAD_ARM_FACTOR must be nonnegative")
TEST_FCORE_SEGMENT_FACTOR >= 0 || error("CHANG_TEST_FCORE_SEGMENT_FACTOR must be nonnegative")
TEST_FCORE_CHORD_FACTOR >= 0 || error("CHANG_TEST_FCORE_CHORD_FACTOR must be nonnegative")
