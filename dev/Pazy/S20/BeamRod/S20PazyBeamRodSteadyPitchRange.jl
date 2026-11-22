using AeroBeams, JLD2

# Sweep angle [rad]
Λ = deg2rad(20)

# Root pitch angle offset (assumed) [rad]
θoffset = deg2rad(0.5)

# Root pitch angle range [rad]
θRange = deg2rad.([0,1,3,5,7,10])

# Airspeed range [m/s]
URange = vcat(1e-3,1:1:70)

# Aerodynamic solver
aeroSolver = Indicial()

# Mesh
mesh = :refined

# Number of elements of the wingtip rod beam
nElementsWingtipRod = 12

# Geometric properties
nElem,L,chord,referenceAxis,nodalPositions = geometrical_properties_Pazy(Λ=Λ, mesh=mesh)

# Tip mass configuration
tipMassConfig = :TE

# Strain gauge coordinates on the cross-section (the spar cs is 60 x 2.25 mm)
x2SG = -50e-3/2      # TE strains
x3SG = 2.25e-3/2     # On top of spar

# Strain gauges' number, corresponding element numbers, and spanwise positions
gaugeNumber,SGelem,gaugeSpanPct = mesh == :refined ? 
                                 (7,12,40) : 
                                 (10,9,60)

# Initialize outputs
problem = Array{SteadyProblem}(undef,length(θRange),length(URange))
tipOOP = Array{Float64}(undef,length(θRange),length(URange))
tipAoA = Array{Float64}(undef,length(θRange),length(URange))
epsSG = Array{Float64}(undef,length(θRange),length(URange))

# Sweep pitch angle
for (i,θ) in enumerate(θRange)
    # Model
    model,_ = create_Pazy(aeroSolver=aeroSolver,sweepAngle=Λ,pitchAngle=θ+θoffset,tipMassConfig=tipMassConfig,mesh=mesh,beamWingtipRod=true,nElementsWingtipRod=nElementsWingtipRod)
    # Initialize guess solution
    x0 = zeros(0)
    # Sweep airspeed
    for (j,U) in enumerate(URange)
        # Display progress
        println("Solving for θ = $(round(rad2deg(θ),digits=1)) deg, U = $U m/s")
        # Update airspeed
        update_v_A!(model, v_A=[0,U,0])
        # Create and solve problem
        problem[i,j] = create_SteadyProblem(model=model,x0=x0)
        solve!(problem[i,j])
        # Final states
        state = problem[i,j].loadstepHistory.states[end]
        aeroVars = problem[i,j].loadstepHistory.aeroVariables[end]
        # Outputs
        tipOOP[i,j] = -state[nElem].nodal.u2[1]
        tipAoA[i,j] = rad2deg(aeroVars[nElem].unsteadyKin.αₑ)
        ϵ11SG = state[SGelem].midpoint.γ[1]
        κ2SG = state[SGelem].midpoint.κ[2]
        κ3SG = state[SGelem].midpoint.κ[3]
        epsSG[i,j] = (ϵ11SG - κ2SG*x3SG - κ3SG*x2SG)
        # Update initial guess for next airspeed
        x0 = problem[i,j].systemArrays.x
    end
end

# Load static results derived from raw data (https://doi.org/10.5281/zenodo.16912763)
s = load(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS20/static_results.jld2")
# Rows: airspeed, tip OOP displacement [mm], std of tip OOP displacement [mm], std of airspeed [m/s]
tipdisp_aoa0_ref = permutedims(hcat(s["V"][1:6], s["tip"][1:6,3], s["tip_std"][1:6,3], s["V_std"][1:6]))
tipdisp_aoa1_ref = permutedims(hcat(s["V"][7:12], s["tip"][7:12,3], s["tip_std"][7:12,3], s["V_std"][7:12]))
tipdisp_aoa3_ref = permutedims(hcat(s["V"][13:18], s["tip"][13:18,3], s["tip_std"][13:18,3], s["V_std"][13:18]))
tipdisp_aoa5_ref = permutedims(hcat(s["V"][19:24], s["tip"][19:24,3], s["tip_std"][19:24,3], s["V_std"][19:24]))
tipdisp_aoa7_ref = permutedims(hcat(s["V"][25:30], s["tip"][25:30,3], s["tip_std"][25:30,3], s["V_std"][25:30]))
tipdisp_aoa10_ref = permutedims(hcat(s["V"][31:end], s["tip"][31:end,3], s["tip_std"][31:end,3], s["V_std"][31:end]))

tipAoA_aoa0_ref = permutedims(hcat(s["V"][1:6], 0 .+ s["twist"][1:6]))
tipAoA_aoa1_ref = permutedims(hcat(s["V"][7:12], 1 .+ s["twist"][7:12]))
tipAoA_aoa3_ref = permutedims(hcat(s["V"][13:18], 3 .+ s["twist"][13:18]))
tipAoA_aoa5_ref = permutedims(hcat(s["V"][19:24], 5 .+ s["twist"][19:24]))
tipAoA_aoa7_ref = permutedims(hcat(s["V"][25:30], 7 .+ s["twist"][25:30]))
tipAoA_aoa10_ref = permutedims(hcat(s["V"][31:end], 10 .+ s["twist"][31:end]))

# Strains at the TE gauge of the mesh-matching span position (strain_mean_TE: gauges 1..15 from root to tip)
# Rows: airspeed, strain [µε], std of strain [µε], std of airspeed [m/s]
epsSG_aoa0_ref = permutedims(hcat(s["V"][1:6], s["strain_mean_TE"][1:6,gaugeNumber], s["strain_std_TE"][1:6,gaugeNumber], s["V_std"][1:6]))
epsSG_aoa1_ref = permutedims(hcat(s["V"][7:12], s["strain_mean_TE"][7:12,gaugeNumber], s["strain_std_TE"][7:12,gaugeNumber], s["V_std"][7:12]))
epsSG_aoa3_ref = permutedims(hcat(s["V"][13:18], s["strain_mean_TE"][13:18,gaugeNumber], s["strain_std_TE"][13:18,gaugeNumber], s["V_std"][13:18]))
epsSG_aoa5_ref = permutedims(hcat(s["V"][19:24], s["strain_mean_TE"][19:24,gaugeNumber], s["strain_std_TE"][19:24,gaugeNumber], s["V_std"][19:24]))
epsSG_aoa7_ref = permutedims(hcat(s["V"][25:30], s["strain_mean_TE"][25:30,gaugeNumber], s["strain_std_TE"][25:30,gaugeNumber], s["V_std"][25:30]))
epsSG_aoa10_ref = permutedims(hcat(s["V"][31:end], s["strain_mean_TE"][31:end,gaugeNumber], s["strain_std_TE"][31:end,gaugeNumber], s["V_std"][31:end]))

# Load tip angle of attack data (digitized from https://doi.org/10.2514/1.J066240, seems more accurate than derived data)
steady_tipaoa_aoa0_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS20/steady_tipaoa_aoa0.txt")
steady_tipaoa_aoa1_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS20/steady_tipaoa_aoa1.txt")
steady_tipaoa_aoa3_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS20/steady_tipaoa_aoa3.txt")
steady_tipaoa_aoa5_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS20/steady_tipaoa_aoa5.txt")
steady_tipaoa_aoa7_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS20/steady_tipaoa_aoa7.txt")
steady_tipaoa_aoa10_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS20/steady_tipaoa_aoa10.txt")

## Plots
# Set paths
relPath = "/dev/Pazy/S20/BeamRod/outputs/S20PazyBeamRodSteadyPitchRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
colors = cgrad(:rainbow, length(θRange), categorical=true)
ts = 10
fs = 16
lfs = 12
lw = 2
ms = 4
msw = 0
gr()

# Tip OOP displacement vs. airspeed
plt_tipOOP = plot(xlabel="Airspeed [m/s]", ylabel="Tip OOP disp. [% semispan]", xlims=[0,71], xticks=0:10:70, tickfont=font(ts), guidefont=font(fs), legend=:topleft, legendfontsize=lfs)
plot!([NaN], [NaN], lw=lw, c=:black, label="AeroBeams")
scatter!([NaN], [NaN], ms=ms, msw=msw, c=:black, label="Exp.")
for (i,θ) in enumerate(θRange)
    plot!(URange, tipOOP[i,:]/L*100, c=colors[i], lw=lw, ls=:solid, label="\$\\theta = $(round(Int,rad2deg(θ))) ^\\circ\$")
    if θ ≈ 0
        scatter!(tipdisp_aoa0_ref[1,:], tipdisp_aoa0_ref[2,:]*1e-3*100/L, xerror=tipdisp_aoa0_ref[4,:], yerror=tipdisp_aoa0_ref[3,:]*1e-3*100/L, ms=ms, msw=msw, c=colors[i], label=false)
    elseif θ ≈ deg2rad(1)
        scatter!(tipdisp_aoa1_ref[1,:], tipdisp_aoa1_ref[2,:]*1e-3*100/L, xerror=tipdisp_aoa1_ref[4,:], yerror=tipdisp_aoa1_ref[3,:]*1e-3*100/L, ms=ms, msw=msw, c=colors[i], label=false)
    elseif θ ≈ deg2rad(3)
        scatter!(tipdisp_aoa3_ref[1,:], tipdisp_aoa3_ref[2,:]*1e-3*100/L, xerror=tipdisp_aoa3_ref[4,:], yerror=tipdisp_aoa3_ref[3,:]*1e-3*100/L, ms=ms, msw=msw, c=colors[i], label=false)
    elseif θ ≈ deg2rad(5)
        scatter!(tipdisp_aoa5_ref[1,:], tipdisp_aoa5_ref[2,:]*1e-3*100/L, xerror=tipdisp_aoa5_ref[4,:], yerror=tipdisp_aoa5_ref[3,:]*1e-3*100/L, ms=ms, msw=msw, c=colors[i], label=false)
    elseif θ ≈ deg2rad(7)
        scatter!(tipdisp_aoa7_ref[1,:], tipdisp_aoa7_ref[2,:]*1e-3*100/L, xerror=tipdisp_aoa7_ref[4,:], yerror=tipdisp_aoa7_ref[3,:]*1e-3*100/L, ms=ms, msw=msw, c=colors[i], label=false)
    elseif θ ≈ deg2rad(10)
        scatter!(tipdisp_aoa10_ref[1,:], tipdisp_aoa10_ref[2,:]*1e-3*100/L, xerror=tipdisp_aoa10_ref[4,:], yerror=tipdisp_aoa10_ref[3,:]*1e-3*100/L, ms=ms, msw=msw, c=colors[i], label=false)
    end
end
display(plt_tipOOP)
savefig(string(absPath,"/S20PazyBeamRodSteadyPitchRange_tipOOP.pdf"))

# Tip AoA vs. airspeed
plt_tipAOA = plot(xlabel="Airspeed [m/s]", ylabel="Tip angle of attack [deg]", xlims=[0,71], xticks=0:10:70, ylims=[-1,10], yticks=0:2:10, tickfont=font(ts), guidefont=font(fs))
for (i,θ) in enumerate(θRange)
    plot!(URange, tipAoA[i,:], c=colors[i], lw=lw, label=false)
    if θ ≈ 0
        scatter!(steady_tipaoa_aoa0_ref[1,:], steady_tipaoa_aoa0_ref[2,:], ms=ms, msw=msw, c=colors[i], label=false)
    elseif θ ≈ deg2rad(1)
        scatter!(steady_tipaoa_aoa1_ref[1,:], steady_tipaoa_aoa1_ref[2,:], ms=ms, msw=msw, c=colors[i], label=false)
    elseif θ ≈ deg2rad(3)
        scatter!(steady_tipaoa_aoa3_ref[1,:], steady_tipaoa_aoa3_ref[2,:], ms=ms, msw=msw, c=colors[i], label=false)
    elseif θ ≈ deg2rad(5)
        scatter!(steady_tipaoa_aoa5_ref[1,:], steady_tipaoa_aoa5_ref[2,:], ms=ms, msw=msw, c=colors[i], label=false)
    elseif θ ≈ deg2rad(7)
        scatter!(steady_tipaoa_aoa7_ref[1,:], steady_tipaoa_aoa7_ref[2,:], ms=ms, msw=msw, c=colors[i], label=false)
    elseif θ ≈ deg2rad(10)
        scatter!(steady_tipaoa_aoa10_ref[1,:], steady_tipaoa_aoa10_ref[2,:], ms=ms, msw=msw, c=colors[i], label=false)
    end
end
display(plt_tipAOA)
savefig(string(absPath,"/S20PazyBeamRodSteadyPitchRange_tipAoA.pdf"))

# Strains at the TE gauge located at gaugeSpanPct span vs. airspeed
plt_eps = plot(xlabel="Airspeed [m/s]", ylabel="Strains @ $gaugeSpanPct% span (\$\\mu\\varepsilon\$)", xlims=[0,71], xticks=0:10:70, tickfont=font(ts), guidefont=font(fs), legend=:topleft, legendfontsize=lfs)
for (i,θ) in enumerate(θRange)
    plot!(URange, epsSG[i,:]*1e6, c=colors[i], lw=lw, ls=:solid, label=false)
    if θ ≈ deg2rad(0)
        scatter!(epsSG_aoa0_ref[1,:], epsSG_aoa0_ref[2,:], xerror=epsSG_aoa0_ref[4,:], yerror=epsSG_aoa0_ref[3,:], ms=ms, msw=msw, c=colors[i], label=false)
    elseif θ ≈ deg2rad(1)
        scatter!(epsSG_aoa1_ref[1,:], epsSG_aoa1_ref[2,:], xerror=epsSG_aoa1_ref[4,:], yerror=epsSG_aoa1_ref[3,:], ms=ms, msw=msw, c=colors[i], label=false)
    elseif θ ≈ deg2rad(3)
        scatter!(epsSG_aoa3_ref[1,:], epsSG_aoa3_ref[2,:], xerror=epsSG_aoa3_ref[4,:], yerror=epsSG_aoa3_ref[3,:], ms=ms, msw=msw, c=colors[i], label=false)
    elseif θ ≈ deg2rad(5)
        scatter!(epsSG_aoa5_ref[1,:], epsSG_aoa5_ref[2,:], xerror=epsSG_aoa5_ref[4,:], yerror=epsSG_aoa5_ref[3,:], ms=ms, msw=msw, c=colors[i], label=false)
    elseif θ ≈ deg2rad(7)
        scatter!(epsSG_aoa7_ref[1,:], epsSG_aoa7_ref[2,:], xerror=epsSG_aoa7_ref[4,:], yerror=epsSG_aoa7_ref[3,:], ms=ms, msw=msw, c=colors[i], label=false)
    elseif θ ≈ deg2rad(10)
        scatter!(epsSG_aoa10_ref[1,:], epsSG_aoa10_ref[2,:], xerror=epsSG_aoa10_ref[4,:], yerror=epsSG_aoa10_ref[3,:], ms=ms, msw=msw, c=colors[i], label=false)
    end
end
display(plt_eps)
savefig(string(absPath,"/S20PazyBeamRodSteadyPitchRange_eps.pdf"))

println("Finished S20PazyBeamRodSteadyPitchRange.jl")