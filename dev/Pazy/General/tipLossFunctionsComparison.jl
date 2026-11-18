using AeroBeams, Plots, DelimitedFiles, ColorSchemes

# Sweep angle
Λ = 0*π/180

# Pitch angle
θ = 0*π/180

# Airspeed
U = 100

# Compressibility factor
M = U/343
β = sqrt(1-M^2)

# Domain
s = 0:0.01:1

# No tip loss
WConst = s -> 1

# Exponential loss
WExp = tip_loss_function_Pazy("Exponential",θ,U,Λ)

# VLM - deformed shape
WVLMdef = tip_loss_function_Pazy("VLM-def",θ,U,Λ)

# Load reference data (from AePW4 meetings)
clalapha_Lambda0_aoa5_NastranTechnion = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazy/clalapha_Lambda0_aoa5_NastranTechnion.txt")
clalapha_Lambda20_aoa5_NastranTechnion = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazy/clalapha_Lambda20_aoa5_NastranTechnion.txt")

# Set paths
relPath = "/dev/Pazy/General/tipLossFunctionsComparison"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot normal force distribution at the specified condition
colors = cgrad(:rainbow, 3, categorical=true)
plt = plot(xlims=[0,1], ylims=[0,1.25], xlabel="Normalized span", ylabel="\$ W(s) = C_{n_\\alpha}/(2\\pi/\\beta)\$", title=string("\$\\Lambda=",round(Int,Λ*180/pi),"^\\circ\$",", \$\\alpha_r=",round(Int,θ*180/pi),"^\\circ\$",", \$U=",U," \$ m/s"), tickfont=font(10), guidefont=font(16), legendfontsize=12)
plot!(s, WConst.(s), c=colors[1], lw=2, ls=:solid, label="No loss")
plot!(s, WExp.(s), c=colors[2], lw=2, ls=:dash, label="Exponential loss")
plot!(s, WVLMdef.(s), c=colors[3], lw=2, ls=:dashdot, label="VLM - deformed")
if Λ == 0 && θ == 5*π/180
    i = Int(U/10)
    scatter!(clalapha_Lambda0_aoa5_NastranTechnion[2*i-1,:],clalapha_Lambda0_aoa5_NastranTechnion[2*i,:]/(2π/β), c=:black, ms=4, msw=0, label="Nastran (Technion)")
elseif Λ == 20*π/180 && θ == 5*π/180
    i = Int(U/10)
    scatter!(clalapha_Lambda20_aoa5_NastranTechnion[2*i-1,:],clalapha_Lambda20_aoa5_NastranTechnion[2*i,:]/(2π/β), c=:black, ms=4, msw=0, label="Nastran (Technion)")
end
display(plt)
savefig(string(absPath,"/cnalpha.pdf"))
