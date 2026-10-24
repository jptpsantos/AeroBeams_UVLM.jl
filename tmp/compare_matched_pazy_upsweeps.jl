using DelimitedFiles, LinearAlgebra, Statistics, Printf, SHA
const BASE = normpath(joinpath(@__DIR__,"..","lib","WingPropellerUVLM","examples","pazy_aerobeams_coupled","output"))
const OUT = joinpath(BASE,"LCO_small_signal_audit")
original_path = joinpath(BASE,"LCO_independent_sweeps","pazy_linear_LCO_upsweep.csv")
current_path = joinpath(BASE,"LCO_downsweep_pulse","pazy_linear_LCO_upsweep.csv")
original_bytes,current_bytes = read(original_path),read(current_path)
a,ha = readdlm(IOBuffer(original_bytes),',',Float64;header=true)
b,hb = readdlm(IOBuffer(current_bytes),',',Float64;header=true)
@assert size(a,1) == size(b,1)
open(joinpath(OUT,"matched_upsweep_comparison.txt"),"w") do io
    println(io,"Original SHA256: ",bytes2hex(sha256(original_bytes)))
    println(io,"Current SHA256: ",bytes2hex(sha256(current_bytes)))
    println(io,"Rows: ",size(a,1))
    for key in vec(ha)
        x,y = a[:,findfirst(==(key),vec(ha))],b[:,findfirst(==(key),vec(hb))]
        @printf(io,"%s: exactly_equal=%s; max_abs_difference=%.15g; RMS_difference=%.15g\n",key,x==y,maximum(abs.(x-y)),sqrt(mean(abs2,x-y)))
    end
    force = b[:,findfirst(==("tip_force_A1_N"),vec(hb))]
    println(io,"Current upsweep force identically zero: ",all(iszero,force))
end
write(joinpath(OUT,"matched_upsweep_snapshot.csv"),current_bytes)
println(read(joinpath(OUT,"matched_upsweep_comparison.txt"),String))
