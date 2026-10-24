using DelimitedFiles, LinearAlgebra, Statistics, Printf
BLAS.set_num_threads(1)
directory = joinpath(@__DIR__, "..", "lib", "WingPropellerUVLM", "examples", "pazy_aerobeams_coupled", "output", "LCO_independent_sweeps")
for name in ("upsweep", "downsweep")
    data, header = readdlm(joinpath(directory, "pazy_linear_LCO_$name.csv"), ',', Float64; header=true)
    col(key) = data[:,findfirst(==(key), vec(header))]
    t, u, strain = col("time_s"), col("airspeed_mps"), col("root_longitudinal_strain_ue")
    reconstruction = 1e6 .* (col("root_axial_strain") .+ 0.001125 .* col("root_curvature_2_per_m"))
    @assert all(isfinite, data)
    @assert all(diff(t) .> 0)
    @assert all((name == "upsweep" ? (>) : (<)).(diff(u), 0))
    @assert maximum(abs.(strain .- reconstruction)) < 1e-8
    @printf("%s: samples=%d, speed=%.6f->%.6f, duration=%.6f, dt=%.12g, reconstruction error=%.3g\n", name, length(t), first(u), last(u), last(t), t[2]-t[1], maximum(abs.(strain .- reconstruction)))
    best = (amplitude=-Inf, speed=NaN)
    for center in 35.5:0.5:51.5
        indices = findall(abs.(u .- center) .< 0.25)
        x = t[indices] .- mean(t[indices])
        design = hcat(ones(length(x)), x, x.^2)
        residual = strain[indices] .- design*(design\strain[indices])
        amplitude = (maximum(residual)-minimum(residual))/2
        amplitude > best.amplitude && (best = (; amplitude, speed=center))
        center in (38, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51) &&
            @printf("  U=%.1f, t=%.4f, detrended half range=%.6g ue, RMS=%.6g ue\n", center, mean(t[indices]), amplitude, sqrt(mean(abs2,residual)))
    end
    println("Peak bin: ", best)
end
# Sample the exact current downsweep pulse and speed formulas without loading
# or launching the simulation. These checks do not validate LCO physics.
dt = 10 / ceil(Int, 10 / (0.25 * 0.0989 / 52))
time = collect(0:round(Int,10/dt)) .* dt
speed = 52 .- 17 .* clamp.(time ./ 10, 0, 1)
start = (52-49)/17*10
phase = (time .- start) ./ 0.1
force = [0 < p < 1 ? sin(2pi*p)*sinpi(p)^2 : 0.0 for p in phase]
@assert speed[1] == 52 && speed[end] == 35
@assert all(diff(speed) .< 0)
@assert isapprox(52-1.7*start,49; atol=1e-12)
@assert all(iszero, force[time .<= start])
@assert all(iszero, force[time .>= start+0.1])
@assert abs(sum(force)*dt) < 1e-8
@printf("Current pulse: start=%.9f s, end=%.9f s, U=49->%.6f m/s, dt=%.12g, nonzero samples=%d, min/max=%.8f/%.8f N, net impulse=%.3g Ns\n",start,start+0.1,52-1.7*(start+0.1),dt,count(!iszero,force),minimum(force),maximum(force),sum(force)*dt)
for duration in (0.1, 0.04)
    tau = collect(range(0, duration; length=10001))
    waveform = sin.(2pi .* tau ./ duration) .* sinpi.(tau ./ duration).^2
    for frequency in (28., 30., 31.)
        projection = abs(sum(waveform .* exp.(-2pi*im*frequency .* tau)) * (tau[2]-tau[1]))
        @printf("Pulse spectral magnitude per 1 N coefficient: T=%.2f s, f=%.1f Hz, |integral F exp(-iwt) dt|=%.9g Ns\n", duration, frequency, projection)
        duration == 0.1 && frequency == 30 && @assert projection < 1e-12
    end
end
