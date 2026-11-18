using AeroBeams

# Sweep angle [rad]
Λ = deg2rad(10)

# Tip mass configurations
tipMassConfigs = [:LE, :TE]

# Flag for ad hoc correction of torsional inertia of wingtip rod
sweptTipMassInertiaCorrection = true

# Mesh
mesh = :nominal

# Root pitch angle range [rad]
θRange = deg2rad.([0,1,3,5,7,8.5,10])

# Airspeed range [m/s]
URange = collect(30:1:90)

# Tip correction function type
tipLossType = "VLM-def"

# Aerodynamic solver
aeroSolver = Indicial()

# Airfoil section
airfoil = deepcopy(NACA0018)

# Geometric properties
nElem,L,chord,normSparPos,nodalPositions = geometrical_properties_Pazy(Λ=Λ,mesh=mesh)

# Number of modes to track
nModes = 3

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :krylovkit
nev = 40

# Initialize outputs
eigensols = Array{EigenSolution}(undef,length(tipMassConfigs),length(θRange),length(URange))
tipOOP = Array{Float64}(undef,length(tipMassConfigs),length(θRange),length(URange))
freqsOsc = Array{Matrix{Float64}}(undef,length(tipMassConfigs),length(θRange))
dampsOsc = Array{Matrix{Float64}}(undef,length(tipMassConfigs),length(θRange))
flutterOnsetSpeed = fill(Inf,length(tipMassConfigs),length(θRange))
flutterOnsetFreq = fill(NaN,length(tipMassConfigs),length(θRange))
flutterOnsetTipOOP = fill(NaN,length(tipMassConfigs),length(θRange))

# Sweep tip mass configurations
for (c,config) in enumerate(tipMassConfigs)
    # Sweep root pitch angle
    for (i,θ) in enumerate(θRange)
        # Initialize guess solution
        x0 = zeros(0)
        # Sweep airspeed
        for (j,U) in enumerate(URange)
            # Display progress
            println("Solving for tip mass at $config, θ = $(round(rad2deg(θ),digits=1)) deg, U = $U m/s")
            # Model
            model,_ = create_Pazy(aeroSolver=aeroSolver,mesh=mesh,airfoil=airfoil,upright=true,hasTipCorrection=true,tipLossType=tipLossType,sweepAngle=Λ,pitchAngle=θ,airspeed=U,sweptTipMassConfig=config,sweptTipMassInertiaCorrection=sweptTipMassInertiaCorrection)
            # Create and solve problem
            problem = create_EigenProblem(model=model,eigenmethod=eigenmethod,nev=nev,x0=x0)
            solve!(problem)
            # Eigensolution
            eigensols[c,i,j] = problem.eigensolution
            # Tip OOP deflection
            tipOOP[c,i,j] = -problem.loadstepHistory.states[end][nElem].nodal.u2[1]
            # Set initial guess solution for next iteration
            x0 = problem.systemArrays.x
        end
        # Apply mode tracking
        freqs,damps,_ = mode_tracking(eigensols[c,i,:])
        # Get frequencies and dampings of oscillatory modes
        freqsOsc[c,i],dampsOsc[c,i],_ = filter_tracked_modes(freqs,damps; fmin=Hz2rads(3.75), fmax=Hz2rads(50), maxTracks=nModes)
        # Flutter onset speed, frequency and tip OOP (lowest onset)
        onsetSpeeds,onsetFreqs,_,onsetTipOOPs,_ = find_flutter(URange,freqsOsc[c,i],dampsOsc[c,i],tipOOP[c,i,:])
        if !isempty(onsetSpeeds)
            flutterOnsetSpeed[c,i] = onsetSpeeds[1]
            flutterOnsetFreq[c,i] = onsetFreqs[1]
            flutterOnsetTipOOP[c,i] = onsetTipOOPs[1]
        end
    end
end

# Load reference data (experiments from Technion - https://doi.org/10.5281/zenodo.16354530)
Uf_vs_aoa_LE = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS10/Uf_vs_aoa_LE.txt")
Uf_vs_aoa_TE = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS10/Uf_vs_aoa_TE.txt")
Uf_vs_tipdisp_LE = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS10/Uf_vs_tipdisp_LE.txt")
Uf_vs_tipdisp_TE = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS10/Uf_vs_tipdisp_TE.txt")
omegaf_vs_aoa_LE = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS10/omegaf_vs_aoa_LE.txt")
omegaf_vs_aoa_TE = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS10/omegaf_vs_aoa_TE.txt")
omegaf_vs_tipdisp_LE = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS10/omegaf_vs_tipdisp_LE.txt")
omegaf_vs_tipdisp_TE = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS10/omegaf_vs_tipdisp_TE.txt")

## Plots
# Set paths
relPath = "/dev/Pazy/S10/InertiaRod/outputs/S10PazyFlutterOnsetPitchRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
configColors = cgrad(:rainbow, length(tipMassConfigs), categorical=true)
ts = 10
fs = 10
lfs = 8
lw = 2
ms = 4
msw = 1
gr()

# Flutter onset speed vs. AoA
plt1 = plot(ylabel="Flutter speed [m/s]", xlims=[0,10], ylims=[0,85], tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, legend=:bottomleft)
plot!([NaN], [NaN], lw=lw, c=:black, label="AeroBeams")
plot!([NaN], [NaN], lw=lw, ls=:dash, c=:black, marker=:circle, ms=ms, msw=msw, mc=:white, msc=:black, label="Exp.")
plot!([NaN], [NaN], lw=lw, c=configColors[1], label="LE weight")
plot!([NaN], [NaN], lw=lw, c=configColors[2], label="TE weight")
for (c,config) in enumerate(tipMassConfigs)
    plot!(rad2deg.(θRange), flutterOnsetSpeed[c,:], c=configColors[c], lw=lw, label=false)
    if config == :LE
        plot!(Uf_vs_aoa_LE[1,:], Uf_vs_aoa_LE[2,:], lw=lw, ls=:dash, c=configColors[c], marker=:circle, ms=ms, msw=msw, mc=:white, msc=configColors[c], label=false)
    elseif config == :TE
        plot!(Uf_vs_aoa_TE[1,:], Uf_vs_aoa_TE[2,:], lw=lw, ls=:dash, c=configColors[c], marker=:circle, ms=ms, msw=msw, mc=:white, msc=configColors[c], label=false)
    end
end

# Flutter onset frequency vs. AoA
plt2 = plot(xlabel="Root pitch angle [deg]", ylabel="Flutter frequency [Hz]", xlims=[0,10], ylims=[0,35], tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
plot!([NaN], [NaN], lw=lw, c=:black, label=false)
plot!([NaN], [NaN], lw=lw, ls=:dash, c=:black, marker=:circle, ms=ms, msw=msw, mc=:white, msc=:black, label=false)
plot!([NaN], [NaN], lw=lw, c=configColors[1], label=false)
plot!([NaN], [NaN], lw=lw, c=configColors[2], label=false)
for (c,config) in enumerate(tipMassConfigs)
    plot!(rad2deg.(θRange), rads2Hz.(flutterOnsetFreq[c,:]), c=configColors[c], lw=lw, label=false)
    if config == :LE
        plot!(omegaf_vs_aoa_LE[1,:], omegaf_vs_aoa_LE[2,:], lw=lw, ls=:dash, c=configColors[c], marker=:circle, ms=ms, msw=msw, mc=:white, msc=configColors[c], label=false)
    elseif config == :TE
        plot!(omegaf_vs_aoa_TE[1,:], omegaf_vs_aoa_TE[2,:], lw=lw, ls=:dash, c=configColors[c], marker=:circle, ms=ms, msw=msw, mc=:white, msc=configColors[c], label=false)
    end
end

# Flutter onset speed vs. normalized tip displacement
plt3 = plot(ylabel="Flutter speed [m/s]", xlims=[0,40], ylims=[0,85], tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, legend=:bottomleft)
plot!([NaN], [NaN], lw=lw, c=:black, label="AeroBeams")
plot!([NaN], [NaN], lw=lw, ls=:dash, c=:black, marker=:circle, ms=ms, msw=msw, mc=:white, msc=:black, label="Exp.")
plot!([NaN], [NaN], lw=lw, c=configColors[1], label="LE weight")
plot!([NaN], [NaN], lw=lw, c=configColors[2], label="TE weight")
for (c,config) in enumerate(tipMassConfigs)
    plot!(flutterOnsetTipOOP[c,:]/L*100, flutterOnsetSpeed[c,:], c=configColors[c], lw=lw, label=false)
    if config == :LE
        plot!(Uf_vs_tipdisp_LE[1,:]*100, Uf_vs_tipdisp_LE[2,:], lw=lw, ls=:dash, c=configColors[c], marker=:circle, ms=ms, msw=msw, mc=:white, msc=configColors[c], label=false)
    elseif config == :TE
        plot!(Uf_vs_tipdisp_TE[1,:]*100, Uf_vs_tipdisp_TE[2,:], lw=lw, ls=:dash, c=configColors[c], marker=:circle, ms=ms, msw=msw, mc=:white, msc=configColors[c], label=false)
    end
end

# Flutter onset frequency vs. tip displacement
plt4 = plot(xlabel="Normalized tip displacement [% semispan]", ylabel="Flutter frequency [Hz]", xlims=[0,40], ylims=[0,35], tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
plot!([NaN], [NaN], lw=lw, c=:black, label=false)
plot!([NaN], [NaN], lw=lw, ls=:dash, c=:black, marker=:circle, ms=ms, msw=msw, mc=:white, msc=:black, label=false)
plot!([NaN], [NaN], lw=lw, c=configColors[1], label=false)
plot!([NaN], [NaN], lw=lw, c=configColors[2], label=false)
for (c,config) in enumerate(tipMassConfigs)
    plot!(flutterOnsetTipOOP[c,:]/L*100, rads2Hz.(flutterOnsetFreq[c,:]), c=configColors[c], lw=lw, label=false)
    if config == :LE
        plot!(omegaf_vs_tipdisp_LE[1,:]*100, omegaf_vs_tipdisp_LE[2,:], lw=lw, ls=:dash, c=configColors[c], marker=:circle, ms=ms, msw=msw, mc=:white, msc=configColors[c], label=false)
    elseif config == :TE
        plot!(omegaf_vs_tipdisp_TE[1,:]*100, omegaf_vs_tipdisp_TE[2,:], lw=lw, ls=:dash, c=configColors[c], marker=:circle, ms=ms, msw=msw, mc=:white, msc=configColors[c], label=false)
    end
end

# Combine plots and save
plt_aoa = plot(plt1, plt2, layout=(2,1))
display(plt_aoa)
savefig(plt_aoa,string(absPath,"/S10PazyFlutterOnsetPitchRange_Uf_omegaf_vs_aoa.pdf"))

plt_disp = plot(plt3, plt4, layout=(2,1))
display(plt_disp)
savefig(plt_disp,string(absPath,"/S10PazyFlutterOnsetPitchRange_Uf_omegaf_vs_tipdisp.pdf"))

println("Finished S10PazyFlutterOnsetPitchRange.jl")