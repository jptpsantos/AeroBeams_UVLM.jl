using DelimitedFiles, LinearAlgebra, Statistics, Printf
BLAS.set_num_threads(1)
directory = joinpath(@__DIR__, "..", "lib", "WingPropellerUVLM", "examples", "pazy_aerobeams_coupled", "output", "LCO_downsweep_pulse")
for name in ("upsweep", "downsweep")
    data, header = readdlm(joinpath(directory, "pazy_linear_LCO_$name.csv"), ',', Float64; header=true)
    col(key) = data[:,findfirst(==(key), vec(header))]
    t, u, strain, force = col("time_s"), col("airspeed_mps"), col("root_longitudinal_strain_ue"), col("tip_force_A1_N")
    @assert all(isfinite,data)
    @printf("%s: U %.1f -> %.1f; T %.6f s; samples %d; force min/max %.6g/%.6g N\n", name, first(u), last(u), last(t), length(t), minimum(force), maximum(force))
    best = (amplitude=-Inf,speed=NaN)
    for center in 35.5:0.5:51.5
        idx = findall(abs.(u .- center) .< 0.25)
        x = t[idx] .- mean(t[idx])
        design = hcat(ones(length(x)),x,x.^2)
        residual = strain[idx] .- design*(design\strain[idx])
        amplitude = (maximum(residual)-minimum(residual))/2
        amplitude > best.amplitude && (best=(;amplitude,speed=center))
        center in (40,43,45,47,48,49,50,51) && @printf("  U %.1f: detrended half range %.6g ue\n",center,amplitude)
    end
    println("Peak bin: ", best)
end
