using AeroBeams

# Set comparison case
case = 1

# Airspeed range [m/s]
URange = vcat(10:10:60)

# Sweep angle
Λ = case == 1 ? 
    deg2rad(0) : 
    deg2rad(20)

# Pitch angle [rad]
θ = deg2rad(5)

# Span domain
s = 0:0.01:1

# Load reference data (from AePW4 meetings)
clalapha_Lambda0_aoa5_NastranTechnion = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazy/clalapha_Lambda0_aoa5_NastranTechnion.txt")
clalapha_Lambda20_aoa5_NastranTechnion = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazy/clalapha_Lambda20_aoa5_NastranTechnion.txt")
refData = case == 1 ? 
          clalapha_Lambda0_aoa5_NastranTechnion : clalapha_Lambda20_aoa5_NastranTechnion
ylim = case == 1 ? 
       7.5 : 
       6

## Plots       
# Set paths
relPath = "/dev/Pazy/General/outputs/tipLossComparisonVLM_Nastran"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot normal force distribution at the specified condition
colors = cgrad(:rainbow, length(URange), categorical=true)
plt = plot(xlims=[0,1], ylims=[0,ylim], xlabel="Normalized span", ylabel="\$C_{n_\\alpha}\$ [1/rad]", title=string("\$\\Lambda=",round(Int,rad2deg(Λ)),"^\\circ\$",", \$\\alpha_r=",round(Int,rad2deg(θ)),"^\\circ\$"), yticks=0:1:ylim, tickfont=font(10), guidefont=font(16), legendfontsize=11, legend=:bottomleft)
plot!([NaN], [NaN], c=:black, lw=2, label="AeroBeams - VLM")
scatter!([NaN], [NaN], c=:black, ms=4, msw=0, label="Nastran (Technion)")
for (i,U) in enumerate(URange)
    W = tip_loss_function_Pazy(:VLM_tan_fixed,θ,U,Λ)
    cnαVLMdef = 2π * W.(s)
    plot!(s, cnαVLMdef, c=colors[i], lw=2, label="\$U=$U\$ m/s")
    j = Int(U/10)
    scatter!(refData[2*j-1,:],refData[2*j,:], c=colors[i], ms=4, msw=0, label=false)
end
# if case != 1
#     plot!(legend=false)
# end
display(plt)
savefig(string(absPath,"/cnalpha_compare_VLM_Nastran_case",case,".pdf"))

println("Finished tipLossComparisonVLM_Nastran.jl")