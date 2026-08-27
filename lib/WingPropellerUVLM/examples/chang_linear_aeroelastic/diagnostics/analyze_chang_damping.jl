using DelimitedFiles
using Statistics

length(ARGS) >= 1 || error("Usage: julia analyze_chang_damping.jl HISTORY.csv [fit_start_s]")
history_path = abspath(ARGS[1])
fit_start = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 0.7

raw, header = readdlm(history_path, ',', header = true)
names = vec(String.(header))
time_index = findfirst(==("time_s"), names)
pitch_index = findfirst(==("propeller_1_pitch_deg"), names)
isnothing(time_index) && error("time_s column not found")
isnothing(pitch_index) && error("propeller_1_pitch_deg column not found")

time = Float64.(raw[:, time_index])
pitch = Float64.(raw[:, pitch_index])
peak_indices = Int[]
for i in 2:(length(pitch) - 1)
    if time[i] >= fit_start && pitch[i] > pitch[i - 1] && pitch[i] >= pitch[i + 1]
        push!(peak_indices, i)
    end
end
length(peak_indices) >= 3 || error("At least three positive peaks are required after $fit_start s")

peak_times = time[peak_indices]
peak_values = abs.(pitch[peak_indices])
valid = peak_values .> eps(Float64)
peak_times = peak_times[valid]
peak_values = peak_values[valid]

design = hcat(ones(length(peak_times)), peak_times)
coefficients = design \ log.(peak_values)
growth_rate = coefficients[2]
frequency = 1 / mean(diff(peak_times))
fitted = design * coefficients
ss_residual = sum(abs2, log.(peak_values) .- fitted)
ss_total = sum(abs2, log.(peak_values) .- mean(log.(peak_values)))
r_squared = ss_total > 0 ? 1 - ss_residual / ss_total : 1.0

println("history = $history_path")
println("fit_start_s = $fit_start")
println("number_of_positive_peaks = $(length(peak_times))")
println("first_peak_deg = $(first(peak_values))")
println("last_peak_deg = $(last(peak_values))")
println("frequency_hz = $frequency")
println("growth_rate_per_s = $growth_rate")
println("log_envelope_r_squared = $r_squared")
