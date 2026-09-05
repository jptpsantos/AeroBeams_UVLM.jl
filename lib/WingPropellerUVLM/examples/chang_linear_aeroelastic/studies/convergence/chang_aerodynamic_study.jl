# Shared configuration and result integrity for rigid aerodynamic studies.
module ChangAerodynamicStudy
using Statistics
using SHA
using TOML
export ChangCoupledAerodynamicOptions, options_from_environment, validate_options,
    aero_float, aero_int, aero_bool, write_metadata, output_matches_options,
    periodic_metrics, phase_rms_error, metadata_options

aero_float(name, default; environment = ENV) = parse(Float64, get(environment, name, string(default)))
aero_int(name, default; environment = ENV) = parse(Int, get(environment, name, string(default)))

function aero_bool(name, default; environment = ENV)
    value = lowercase(strip(get(environment, name, string(default))))
    value in ("1", "true", "yes", "on") && return true
    value in ("0", "false", "no", "off") && return false
    error("$name must be true/false, yes/no, on/off, or 1/0")
end

"""Physical and numerical controls for one rigid coupled aerodynamic run."""
Base.@kwdef struct ChangCoupledAerodynamicOptions
    flow_speed_mps::Float64 = 65.0
    air_density_kgpm3::Float64 = 1.225
    angle_of_attack_deg::Float64 = 3.0
    sideslip_deg::Float64 = 0.0

    wing_span_m::Float64 = 7.5
    wing_root_chord_m::Float64 = 1.8
    wing_tip_chord_m::Float64 = 1.8
    wing_span_panels::Int = 30
    wing_chord_panels::Int = 10
    elastic_axis_fraction::Float64 = 0.30

    propeller_radius_m::Float64 = 1.15
    propeller_chord_m::Float64 = 0.197
    propeller_blades::Int = 4
    propeller_radial_panels::Int = 10
    propeller_chord_panels::Int = 10
    propeller_attachment_eta::Float64 = 0.83
    pylon_length_m::Float64 = 5.6 * 0.3048
    reference_rpm::Float64 = 1217.6962
    reference_speed_mps::Float64 = 65.0
    collective_pitch_offset_deg::Float64 = 0.0

    azimuth_step_deg::Float64 = 5.0
    simulated_revolutions::Int = 8
    averaged_revolutions::Int = 2
    retained_wake_revolutions::Float64 = 2.0
    wake_relaxation::Float64 = 0.1
    finite_core_segment_factor::Float64 = 0.0
    finite_core_chord_factor::Float64 = 0.01
    interaction_on::Bool = true
end

function options_from_environment(environment = ENV)
    defaults = ChangCoupledAerodynamicOptions()
    fields = (
        :flow_speed_mps => "SPEED_MPS",
        :air_density_kgpm3 => "DENSITY_KGPM3",
        :angle_of_attack_deg => "AOA_DEG",
        :sideslip_deg => "BETA_DEG",
        :wing_span_panels => "WING_SPAN_PANELS",
        :wing_chord_panels => "WING_CHORD_PANELS",
        :propeller_radial_panels => "PROP_RADIAL_PANELS",
        :propeller_chord_panels => "PROP_CHORD_PANELS",
        :reference_rpm => "REFERENCE_RPM",
        :reference_speed_mps => "REFERENCE_SPEED_MPS",
        :collective_pitch_offset_deg => "COLLECTIVE_OFFSET_DEG",
        :azimuth_step_deg => "AZIMUTH_STEP_DEG",
        :simulated_revolutions => "SIMULATED_REVOLUTIONS",
        :averaged_revolutions => "AVERAGED_REVOLUTIONS",
        :retained_wake_revolutions => "RETAINED_WAKE_REVOLUTIONS",
        :wake_relaxation => "WAKE_RELAXATION",
        :finite_core_segment_factor => "FCORE_SEGMENT_FACTOR",
        :finite_core_chord_factor => "FCORE_CHORD_FACTOR",
        :interaction_on => "INTERACTION",
    )
    values = map(fields) do (field, suffix)
        default = getfield(defaults, field)
        name = "CHANG_AERO_" * suffix
        value = default isa Bool ? aero_bool(name, default; environment) :
            parse(typeof(default), get(environment, name, string(default)))
        field => value
    end
    return ChangCoupledAerodynamicOptions(; values...)
end

function validate_options(options)
    all(name -> isfinite(getfield(options, name)), fieldnames(typeof(options))) ||
        error("All aerodynamic controls must be finite")
    all(value -> value > 0, (options.wing_span_m, options.wing_root_chord_m,
        options.wing_tip_chord_m, options.propeller_radius_m, options.propeller_chord_m,
        options.propeller_blades, options.pylon_length_m)) || error("Geometry must be positive")
    0 <= options.elastic_axis_fraction <= 1 || error("Elastic axis must be in [0,1]")
    max(options.finite_core_segment_factor, options.finite_core_chord_factor) > 0 ||
        error("A positive finite core is required for free-wake self induction")
    options.flow_speed_mps > 0 || error("Flow speed must be positive")
    options.air_density_kgpm3 > 0 || error("Air density must be positive")
    options.wing_span_panels > 0 || error("Wing spanwise panel count must be positive")
    options.wing_chord_panels > 0 || error("Wing chordwise panel count must be positive")
    options.propeller_radial_panels > 0 || error("Propeller radial panel count must be positive")
    options.propeller_chord_panels > 0 || error("Propeller chordwise panel count must be positive")
    options.reference_rpm > 0 || error("Reference RPM must be positive")
    options.reference_speed_mps > 0 || error("Reference speed must be positive")
    0 < options.azimuth_step_deg <= 90 || error("Azimuth step must be in (0,90]")
    options.simulated_revolutions >= 2 || error("Simulate at least two revolutions")
    1 <= options.averaged_revolutions < options.simulated_revolutions ||
        error("Averaged revolutions must be positive and smaller than the total")
    options.retained_wake_revolutions > 0 || error("Retained wake age must be positive")
    0 <= options.wake_relaxation <= 1 || error("Wake shedding fraction must be in [0,1]")
    options.finite_core_segment_factor >= 0 || error("Segment core factor must be nonnegative")
    options.finite_core_chord_factor >= 0 || error("Chord core factor must be nonnegative")
    0 <= options.propeller_attachment_eta <= 1 || error("Propeller eta must be in [0,1]")

    steps_per_revolution = round(Int, 360 / options.azimuth_step_deg)
    isapprox(
        steps_per_revolution * options.azimuth_step_deg,
        360.0;
        atol = 100eps(Float64),
        rtol = 0.0,
    ) || error("Azimuth step must divide 360 degrees exactly")
    minimum_revolutions = ceil(Int, options.retained_wake_revolutions) +
        max(options.averaged_revolutions, 2) + 1
    options.simulated_revolutions >= minimum_revolutions || error(
        "Simulate at least $minimum_revolutions revolutions: fill the wake before the averaging and periodicity windows")
    return steps_per_revolution
end

options_dict(options) = Dict(string(name) => getfield(options, name)
    for name in fieldnames(typeof(options)))

"""Fingerprint the numerical source, including uncommitted solver edits."""
function source_fingerprint()
    package = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))
    paths = [@__FILE__, joinpath(@__DIR__, "run_chang_coupled_aerodynamic_analysis.jl")]
    for (root, _, files) in walkdir(joinpath(package, "src"))
        append!(paths, [joinpath(root, file) for file in files if endswith(file, ".jl")])
    end
    manifest = joinpath(package, "Manifest.toml")
    isfile(manifest) && push!(paths, manifest)
    buffer = IOBuffer()
    print(buffer, VERSION, '\0')
    for path in sort(paths)
        print(buffer, relpath(path, package), '\0')
        write(buffer, read(path))
        write(buffer, UInt8(0))
    end
    return bytes2hex(sha256(take!(buffer)))
end

const RESULT_FILES = ("coupled_aerodynamic_history.csv", "coupled_aerodynamic_revolutions.csv",
    "coupled_aerodynamic_summary.txt")

function write_metadata(directory, options; expected_fingerprint = source_fingerprint())
    source_fingerprint() == expected_fingerprint || error("Numerical source changed during the aerodynamic run")
    metadata = Dict("schema" => 2, "source_sha256" => expected_fingerprint,
        "options" => options_dict(options),
        "files" => Dict(file => bytes2hex(sha256(read(joinpath(directory, file))))
            for file in RESULT_FILES))
    # The completion marker is published only after all result files exist.
    path = joinpath(directory, "coupled_aerodynamic_metadata.toml")
    open(path * ".tmp", "w") do io
        TOML.print(io, metadata; sorted = true)
    end
    mv(path * ".tmp", path; force = true)
    return path
end

function metadata_options(directory)
    metadata = TOML.parsefile(joinpath(directory, "coupled_aerodynamic_metadata.toml"))
    get(metadata, "schema", 0) == 2 || error("Unsupported aerodynamic metadata; recompute the case")
    metadata["source_sha256"] == source_fingerprint() || error("Aerodynamic numerical source changed")
    for file in RESULT_FILES
        bytes2hex(sha256(read(joinpath(directory, file)))) == metadata["files"][file] ||
            error("Aerodynamic result file changed or is incomplete: $file")
    end
    values = metadata["options"]
    options = ChangCoupledAerodynamicOptions(;
        (name => values[string(name)] for name in fieldnames(ChangCoupledAerodynamicOptions))...)
    validate_options(options)
    return options
end

function output_matches_options(directory, options; warn_on_mismatch = true)
    try
        options_dict(metadata_options(directory)) == options_dict(options) ||
            error("Aerodynamic controls differ from the saved run")
        return true
    catch exception
        exception isa InterruptException && rethrow()
        warn_on_mismatch && @warn "Aerodynamic output cannot be reused" directory exception
        return false
    end
end

"""Compare complete, phase-aligned revolutions after the wake has filled.

The worst RMS change in the final two (or averaging-count) cycle pairs is
retained. Equal means alone cannot hide a changing oscillation amplitude.
"""
function periodic_metrics(values, steps_per_revolution, averaged_revolutions)
    all(isfinite, values) || error("Coefficient history contains nonfinite values")
    length(values) % steps_per_revolution == 0 || error("Incomplete final revolution")
    cycles = reshape(values, steps_per_revolution, :)
    pairs = max(2, averaged_revolutions)
    size(cycles, 2) > pairs || error("Insufficient complete revolutions for periodicity")
    tail = @view cycles[:, (end - averaged_revolutions + 1):end]
    changes = [sqrt(mean(abs2, cycles[:, i] - cycles[:, i - 1]))
        for i in (size(cycles, 2) - pairs + 1):size(cycles, 2)]
    return (; mean = mean(tail), std = std(vec(tail); corrected = false),
        drift = mean(cycles[:, end]) - mean(cycles[:, end - 1]),
        periodic_rms = maximum(changes), phase = vec(mean(tail; dims = 2)))
end

"""RMS waveform difference on the union of both periodic azimuth grids.

Linear interpolation preserves phase zero at the final sample. Using both
grids catches fine-grid oscillations that alias to zero on the coarse grid.
"""
function phase_rms_error(a, b)
    (isempty(a) || isempty(b)) && return NaN
    sample(v, phase) = begin
        index = phase * length(v)
        lower = floor(Int, index)
        fraction = index - lower
        (1 - fraction) * v[mod1(lower, length(v))] + fraction * v[mod1(lower + 1, length(v))]
    end
    phases = sort!(unique!(vcat(collect(1:length(a)) ./ length(a),
        collect(1:length(b)) ./ length(b))))
    return sqrt(mean(phase -> abs2(sample(a, phase) - sample(b, phase)), phases))
end

end
