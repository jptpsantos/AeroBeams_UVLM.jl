# Run from any directory with --project=<repo>/lib/WingPropellerUVLM.
# These are diagnostic counterexamples, not assertions of physical convergence.
using Test, TOML
const ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const STUDY = joinpath(ROOT, "lib", "WingPropellerUVLM", "examples",
    "chang_linear_aeroelastic", "studies", "convergence")
module DampingAudit
include(joinpath(Main.STUDY, "src", "damping_metrics.jl"))
end

t = collect(0.0:0.001:6.0)
function metric(signal)
    DampingAudit.moving_block_metrics(t, signal;
        fit_start_s=0.75, fit_end_s=6.0, minimum_peaks=8,
        minimum_fit_r_squared=0.8, initial_block_size=512,
        size_ratio_lower=0.25, size_ratio_upper=0.5,
        peak_from_start=1, peak_from_end=0, block_duration_s=1.0,
        block_overlap=0.9, frequency_min_hz=3.0, frequency_max_hz=6.0,
        apply_hann_window=true)
end
records = Dict{String,Any}[]
for frequency in (4.1, 4.4, 8.0)
    m = metric(exp.(-0.15 .* t) .* sin.(2pi * frequency .* t))
    record = Dict("actual_frequency_hz"=>frequency,
        "reported_frequency_hz"=>m.frequency_hz,
        "reported_lambda_per_s"=>m.moving_block_lambda_per_s,
        "valid"=>m.valid_for_convergence, "r_squared"=>m.fit_r_squared)
    push!(records, record)
    println(record)
end
signal = exp.(-0.15 .* t) .* sin.(2pi * 4.1 .* t)
offset_result = (;unshifted_valid=metric(signal).valid_for_convergence,
    negative_offset_valid=metric(signal .- 2.0).valid_for_convergence)
println(offset_result)

# Check project isolation without changing the caller's global environment.
isolated = copy(LOAD_PATH)
dependency_results = try
    empty!(LOAD_PATH); append!(LOAD_PATH, ["@", "@stdlib"])
    Dict(name => !isnothing(Base.find_package(name)) for name in ("FFTW", "Plots"))
finally
    empty!(LOAD_PATH); append!(LOAD_PATH, isolated)
end
println(dependency_results)
open(joinpath(@__DIR__, "diagnostics.toml"), "w") do io
    TOML.print(io, Dict("julia_version"=>string(VERSION), "frequencies"=>records,
        "offset"=>Dict(string(k)=>v for (k,v) in pairs(offset_result)),
        "available_in_isolated_project"=>dependency_results))
end
