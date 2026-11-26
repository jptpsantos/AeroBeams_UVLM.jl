using AeroBeams

# Sweep angle
Λ = deg2rad(20)

# Pitch angle
θ = deg2rad(5)

# Airspeed [m/s]
U = 50

# Span domain
s = 0:0.01:1

# No tip loss
WConst = s -> 1

# VLM distribution
W_VLM = tip_loss_function_Pazy(:VLM_tan_fixed,θ,U,Λ)

# Load reference data (from AePW4 meetings)
clalapha_Lambda0_aoa5_NastranTechnion = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazy/clalapha_Lambda0_aoa5_NastranTechnion.txt")
clalapha_Lambda20_aoa5_NastranTechnion = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazy/clalapha_Lambda20_aoa5_NastranTechnion.txt")

# ID string
IDstr = string("Lambda",round(Int,rad2deg(Λ)),"_aoa",round(Int,rad2deg(θ)),"_U",U)

## Plots
# Set paths
relPath = "/dev/Pazy/General/outputs/tipLossFunctionsComparison"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot normalized normal force distribution at the specified condition
colors = cgrad(:rainbow, 2, categorical=true)
plt = plot(xlims=[0,1], ylims=[0,1.25], xlabel="Normalized span", ylabel="\$ W(s) = C_{n_\\alpha}/2\\pi\$", title=string("\$\\Lambda=",round(Int,rad2deg(Λ)),"^\\circ\$",", \$\\alpha_r=",round(Int,rad2deg(θ)),"^\\circ\$",", \$U=",U," \$ m/s"), tickfont=font(10), guidefont=font(16), legendfontsize=12)
plot!(s, WConst.(s), c=colors[1], lw=2, ls=:solid, label="No loss")
plot!(s, W_VLM.(s), c=colors[2], lw=2, ls=:dashdot, label="VLM")
if Λ ≈ deg2rad(0) && θ ≈ deg2rad(5)
    i = Int(U/10)
    scatter!(clalapha_Lambda0_aoa5_NastranTechnion[2*i-1,:],clalapha_Lambda0_aoa5_NastranTechnion[2*i,:]/2π, c=:black, ms=4, msw=0, label="Nastran (Technion)")
elseif Λ ≈ deg2rad(20) && θ ≈ deg2rad(5)
    i = Int(U/10)
    scatter!(clalapha_Lambda20_aoa5_NastranTechnion[2*i-1,:],clalapha_Lambda20_aoa5_NastranTechnion[2*i,:]/2π, c=:black, ms=4, msw=0, label="Nastran (Technion)")
end
display(plt)
savefig(string(absPath,"/cnalpha_",IDstr,".pdf"))

println("Finished tipLossFunctionsComparison.jl")