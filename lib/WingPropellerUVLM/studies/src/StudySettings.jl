# Public settings use camelCase; the backend and existing file schemas retain
# their established names/units. This is the single translation boundary.
const SETTING_NAMES=Dict(
    :speed_mps=>:airspeed, :end_time_s=>:finalTime,
    :trim_revolutions=>:baselineRevolutions,
    :trim_average_revolutions=>:baselineAverageRevolutions,
    :spanwise_panels=>:nSpanwisePanels, :chordwise_panels=>:nChordwisePanels,
    :radial_panels=>:nRadialPanels, :core_radius_m=>:coreRadius,
    :wing_span=>:wingSpanPanels, :wing_chord=>:wingChordPanels,
    :prop_radial=>:propellerRadialPanels, :prop_chord=>:propellerChordPanels,
    :azimuth_deg=>:azimuthStepDeg)
const LEGACY_NAMES=Dict(v=>k for (k,v) in SETTING_NAMES)
function canonical_key(k::Symbol)
    haskey(SETTING_NAMES,k) && return SETTING_NAMES[k]
    parts=split(string(k),'_')
    return Symbol(first(parts)*join(uppercasefirst.(parts[2:end])))
end
function legacy_key(k::Symbol)
    haskey(LEGACY_NAMES,k) && return LEGACY_NAMES[k]
    k in (:CT,:CQ) && return k
    separated=replace(string(k),r"([A-Z])([A-Z][a-z])"=>s"\1_\2")
    return Symbol(lowercase(replace(separated,r"([a-z0-9])([A-Z])"=>s"\1_\2")))
end
function translate_settings(x::NamedTuple,keyfn)
    pairsOut=Pair{Symbol,Any}[]
    seen=Set{Symbol}()
    for (key,value) in pairs(x)
        mapped=keyfn(key)
        mapped in seen && throw(ArgumentError("Duplicate aliases for $mapped"))
        push!(seen,mapped)
        converted=key==:families ? keyfn.(value) : translate_settings(value,keyfn)
        push!(pairsOut,mapped=>converted)
    end
    return (;pairsOut...)
end
translate_settings(x,keyfn)=deepcopy(x)
canonical_settings(x::NamedTuple)=translate_settings(x,canonical_key)
legacy_settings(x::NamedTuple)=translate_settings(x,legacy_key)
function merge_settings(base::NamedTuple,overrides::NamedTuple)
    unknown=setdiff(keys(overrides),keys(base))
    isempty(unknown) || throw(ArgumentError("Unknown settings: $unknown"))
    return (; (k=>(hasproperty(overrides,k) ?
        (v isa NamedTuple && overrides[k] isa NamedTuple ? merge_settings(v,overrides[k]) : deepcopy(overrides[k])) : deepcopy(v))
        for (k,v) in pairs(base))...)
end

abstract type ConvergenceStudy end
mutable struct AerodynamicConvergenceStudy{S} <: ConvergenceStudy
    settings::S
    results::Union{Nothing,NamedTuple}
end
mutable struct AeroelasticConvergenceStudy{S} <: ConvergenceStudy
    settings::S
    results::Union{Nothing,NamedTuple}
end
study_settings(s::ConvergenceStudy)=deepcopy(s.settings)
legacy_settings(s::ConvergenceStudy)=legacy_settings(s.settings)

function create_AerodynamicConvergenceStudy(settings::NamedTuple;kwargs...)
    resolved=merge_settings(canonical_settings(settings),canonical_settings((;kwargs...)))
    Convergence.validate_aerodynamic(legacy_settings(resolved))
    return AerodynamicConvergenceStudy(resolved,nothing)
end
create_AerodynamicConvergenceStudy(;kwargs...)=create_AerodynamicConvergenceStudy((;kwargs...))
function create_AeroelasticConvergenceStudy(settings::NamedTuple;kwargs...)
    resolved=merge_settings(canonical_settings(settings),canonical_settings((;kwargs...)))
    legacy=legacy_settings(resolved)
    hasproperty(legacy,:selection_file) || throw(ArgumentError("selectionFile is required"))
    MovingBlockDamping.moving_block_metrics(Float64[],Float64[];moving_block_options(legacy.damping)...)
    return AeroelasticConvergenceStudy(resolved,nothing)
end
create_AeroelasticConvergenceStudy(;kwargs...)=create_AeroelasticConvergenceStudy((;kwargs...))

function check_loaded_source()
    Convergence.SOURCE()==LOADED_SOURCE || error("Numerical source changed after this module was loaded; start a fresh Julia process")
end
function solve!(s::AerodynamicConvergenceStudy)
    check_loaded_source()
    settings=legacy_settings(s)
    settings.make_plots && (@eval using Plots)
    s.results=Base.invokelatest(Convergence.run_aerodynamic,settings)
    return s
end
function solve!(s::AeroelasticConvergenceStudy)
    check_loaded_source()
    settings=legacy_settings(s)
    settings.make_plots && (@eval using Plots)
    s.results=run_aeroelastic(settings)
    return s
end

"""Read an archived selection without asserting current numerical verification.
The returned record cannot be used as a verified selection by solve!.
"""
function inspect_selection(path)
    payload=Convergence.TOML.parsefile(path)
    return (;verifiedForCurrentSource=false,record=payload)
end
