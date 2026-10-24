using DelimitedFiles, LinearAlgebra, Statistics, Printf, SHA
import Plots
BLAS.set_num_threads(1)
base = joinpath(@__DIR__, "..", "lib", "WingPropellerUVLM", "examples", "pazy_aerobeams_coupled")
out = joinpath(base,"output","LCO_upsweep_investigation")
mkpath(out)
traces = []
details = String[]
for (folder, branch, label, color, style) in (
    ("LCO_independent_sweeps", "upsweep", "Original unforced upsweep", :blue, :solid),
    ("LCO_downsweep_pulse", "upsweep", "Unforced upsweep in pulse run", :darkcyan, :dash),
    ("LCO_downsweep_pulse", "downsweep", "Pulsed downsweep", :red, :solid))
    path = joinpath(base,"output",folder,"pazy_linear_LCO_$branch.csv")
    bytes = read(path)
    data, header = readdlm(IOBuffer(bytes),',',Float64;header=true)
    col(name) = data[:,something(findfirst(==(name),vec(header)))]
    time, speed, strain = col("time_s"), col("airspeed_mps"), col("root_longitudinal_strain_ue")
    @assert all(isfinite,data) && all(diff(time) .> 0)
    @assert all((branch == "upsweep" ? (>) : (<)).(diff(speed),0))
    @assert maximum(abs.(strain .- 1e6 .* (col("root_axial_strain") .+ 0.001125 .* col("root_curvature_2_per_m")))) < 1e-8
    info = IOBuffer()
    println(info,label,"\nSource: ",abspath(path),"\nSHA256: ",bytes2hex(sha256(bytes)))
    @printf(info,"N=%d; U=%.3f -> %.3f m/s; T=%.6f s; dt=%.12g s; rate=%.6f m/s^2\n", length(time),first(speed),last(speed),last(time),time[2]-time[1],(last(speed)-first(speed))/(last(time)-first(time)))
    @printf(info,"Initial strain=%.6f ue; initial tip displacement=%.9f m; initial tip twist=%.6f deg\n",first(strain),first(col("tip_out_of_plane_m")),first(col("tip_twist_deg")))
    force_idx = findfirst(==("tip_force_A1_N"),vec(header))
    if !isnothing(force_idx)
        force = data[:,force_idx]
        active = findall(!iszero,force)
        @printf(info,"Force min/max=%.9f/%.9f N; nonzero samples=%d\n",minimum(force),maximum(force),length(active))
        isempty(active) || @printf(info,"Sampled pulse interval: t=%.9f -> %.9f s, U=%.6f -> %.6f m/s\n",time[first(active)],time[last(active)],speed[first(active)],speed[last(active)])
    else
        println(info,"No force column in this older CSV; original runner is unforced.")
    end
    rows = []
    for center in 35.5:0.5:51.5
        idx = findall(abs.(speed .- center) .< 0.25)
        x = time[idx] .- mean(time[idx])
        design = hcat(ones(length(idx)),x,x.^2)
        residual = strain[idx] .- design*(design\strain[idx])
        half_range = (maximum(residual)-minimum(residual))/2
        rms = sqrt(mean(abs2,residual))
        push!(rows,[center,mean(time[idx]),half_range,rms])
        center in (43,45,47,48,49,50,51) && @printf(info,"U=%.1f: detrended half range=%.9g ue; RMS=%.9g ue\n",center,half_range,rms)
    end
    metrics = reduce(vcat,permutedims.(rows))
    peak = argmax(metrics[:,3])
    @printf(info,"Largest bin half range=%.6f ue at U=%.1f m/s, t=%.6f s\n",metrics[peak,3],metrics[peak,1],metrics[peak,2])
    writedlm(joinpath(out,"$(folder)_$(branch)_windows.csv"),vcat(
        ["speed_bin_mps" "mean_time_s" "detrended_half_range_ue" "detrended_RMS_ue"],metrics),',')
    push!(traces,(;metrics,label,color,style))
    push!(details,String(take!(info)))
end
figure = Plots.plot(;xlabel="Airspeed [m/s]",ylabel="Detrended half range [microstrain]",
    yscale=:log10,yticks=10.0 .^ (-5:3),legend=:topleft,size=(950,600),dpi=200,
    title="Saved responses: 0.5 m/s windows, quadratic trend removed",
    titlefontsize=11,margin=6Plots.mm,top_margin=10Plots.mm)
for trace in traces
    Plots.plot!(figure,trace.metrics[:,1],trace.metrics[:,3];label=trace.label,
        color=trace.color,linestyle=trace.style,linewidth=2)
end
Plots.savefig(figure,joinpath(out,"oscillation_amplitude_comparison.png"))
write(joinpath(out,"measurements.txt"),join(details,"\n") *
    "\nWindow measures describe transient residuals, including trend error/numerical noise at tiny amplitudes; they are not settled LCO amplitudes or eigenvalue estimates.\n")
println(join(details,"\n"))
println("Outputs: ",abspath(out))
