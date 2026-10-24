using DelimitedFiles, LinearAlgebra, Statistics, Printf
directory = joinpath(@__DIR__,"..","lib","WingPropellerUVLM","examples","pazy_aerobeams_coupled","output","LCO_independent_sweeps")
rows = Vector{Any}[]
for name in ("upsweep","downsweep")
    data,header = readdlm(joinpath(directory,"pazy_linear_LCO_$name.csv"),',',Float64;header=true)
    col(key) = data[:,findfirst(==(key),vec(header))]
    t,u,y = col("time_s"),col("airspeed_mps"),col("root_longitudinal_strain_ue")
    @assert all(isfinite,data)
    @assert all((name=="upsweep" ? (>) : (<)).(diff(u),0))
    @assert maximum(abs.(y .- 1e6 .* (col("root_axial_strain") .+ 0.001125 .* col("root_curvature_2_per_m")))) < 1e-8
    @printf("%s: %d samples, dt %.12g, speed %.1f -> %.1f, duration %.8f\n",name,length(t),t[2]-t[1],first(u),last(u),last(t))
    best = (amplitude=-Inf,speed=NaN,time=NaN)
    for center in 35.5:0.5:59.5
        indices = findall(abs.(u .- center) .< 0.25)
        length(indices)>=4 || continue
        x=t[indices] .- mean(t[indices])
        design=hcat(ones(length(x)),x,x.^2)
        residual=y[indices] .- design*(design\y[indices])
        amplitude=(maximum(residual)-minimum(residual))/2
        push!(rows,Any[name,center,mean(t[indices]),amplitude,sqrt(mean(abs2,residual))])
        amplitude>best.amplitude && (best=(;amplitude,speed=center,time=mean(t[indices])))
        if amplitude>1 || center in (40,43,45,47,49,51,52,55,59)
            @printf("  %.1f m/s: detrended amplitude %.9g microstrain at %.4f s\n",center,amplitude,mean(t[indices]))
        end
    end
    println("Peak bin: ",best)
end
writedlm(joinpath(directory,"latest_sweep_window_analysis.csv"),
    vcat(permutedims(["branch","bin_center_speed_mps","bin_center_time_s","detrended_half_range_ue","detrended_RMS_ue"]),permutedims(hcat(rows...))),',')
