using AeroBeams, Plots, DelimitedFiles, ColorSchemes

# Set case
case = 2

# Airspeed range
URange = vcat(10:10:60)

# Sweep angle
Λ = case == 1 ? 0*π/180 : 20*π/180

# Pitch angle
θ = 5*π/180

# Domain
s = 0:0.01:1

# Load reference data (from AePW4 meetings)
clalapha_Lambda0_aoa5_NastranTechnion = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazy/clalapha_Lambda0_aoa5_NastranTechnion.txt")
clalapha_Lambda20_aoa5_NastranTechnion = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazy/clalapha_Lambda20_aoa5_NastranTechnion.txt")

ylim = 0
refData = 0
if Λ == 0 && θ == 5*π/180
    refData = clalapha_Lambda0_aoa5_NastranTechnion
    ylim = 7.5
elseif Λ == 20*π/180 && θ == 5*π/180
    refData = clalapha_Lambda20_aoa5_NastranTechnion
    ylim = 6
end

# Set paths
relPath = "/dev/Pazy/General/tipLossComparisonVLM_Nastran"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot normal force distribution at the specified condition
colors = cgrad(:rainbow, length(URange), categorical=true)
plt = plot(xlims=[0,1], ylims=[0,ylim], xlabel="Normalized span", ylabel="\$C_{n_\\alpha}\$ [1/rad]", title=string("\$\\Lambda=",round(Int,Λ*180/pi),"^\\circ\$",", \$\\alpha_r=",round(Int,θ*180/pi),"^\\circ\$"), yticks=0:1:ylim, tickfont=font(10), guidefont=font(16), legendfontsize=11, legend=:bottomleft)
plot!([NaN], [NaN], c=:black, lw=2, label="Present VLM")
scatter!([NaN], [NaN], c=:black, ms=4, msw=0, label="Nastran (Technion)")
for (i,U) in enumerate(URange)
    β = sqrt(1-(U/343)^2)
    W = tip_loss_function_Pazy("VLM-def",θ,U,Λ)
    cnαVLMdef = 2π/β * W.(s)
    plot!(s, cnαVLMdef, c=colors[i], lw=2, label="\$U=$U\$ m/s")
    j = Int(U/10)
    scatter!(refData[2*j-1,:],refData[2*j,:], c=colors[i], ms=4, msw=0, label=false)
end
# if case != 1
#     plot!(legend=false)
# end
display(plt)
savefig(string(absPath,"/cnalpha_compare_VLM_Nastran_case",case,".pdf"))
