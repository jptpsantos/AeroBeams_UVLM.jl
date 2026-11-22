using AeroBeams, LinearInterpolations, JLD2

# Reference data IDs (available: 3020, 3021, 3022, 3023, 3031, 3033, 3034, 3035, 3036, 3037, 3038, 3201, 3202, 3203)
refIDs = [3023]

# Angle of sweep [rad]
Λ = deg2rad(20)

# Flag to smooth airspeed signal and corresponding moving average window
smoothAirspeedSignal = true
smoothingWindow = Int(2e3)

# Time integration method
integrationMethod = GeneralizedAlpha(ρ∞=0.9)

# Mesh: :nominal or :refined
mesh = :refined

# Aerodynamic solver
aeroSolver = Indicial()

# Structural damping coefficients for the wing
μWing = 0

# Time variables [s]
Δt = 1e-4
trackingFrequency = clamp(round(Int,1e-3/Δt), 1, 10)

# Set system solver options
systemSolver = aeroSolver isa DynamicStallSolver ?
               create_NewtonRaphson(absTol=1e-7,displayStatus=false) :
               create_NLsolver(method=:newton,linesearch=false,displayStatus=false)

# Root strain gauge coordinates on the cross-section (the spar cs is 60 x 2.25 mm)
x2SG_LE = 50e-3/2    # x2-position: on LE
x2SG_TE = -50e-3/2   # x2-position: on TE
x3SG = 2.25e-3/2     # x3-position: on top

# Strain gauge used in the comparisons, and the element whose midpoint (where strains are computed) is closest to it: g7 (≈40% of the semispan) for the :refined mesh, g10 (≈60%) for the :nominal mesh
gaugeName,SGelem,gaugeSpanPct = mesh == :refined ? 
                                ("g7",12,40) : 
                                ("g10",9,60)

# Initialize outputs
dynamicProblem = Array{DynamicProblem}(undef,length(refIDs))
Uprofile = Array{Function}(undef,length(refIDs))
t = Array{Vector{Float64}}(undef,length(refIDs))
tipAoA = Array{Vector{Float64}}(undef,length(refIDs))
tipOOP = Array{Vector{Float64}}(undef,length(refIDs))
tipTwist = Array{Vector{Float64}}(undef,length(refIDs))
rootEpsLE = Array{Vector{Float64}}(undef,length(refIDs))
rootEpsTE = Array{Vector{Float64}}(undef,length(refIDs))
SGEpsLE = Array{Vector{Float64}}(undef,length(refIDs))
SGEpsTE = Array{Vector{Float64}}(undef,length(refIDs))

# Sweep IDs
for (i,refID) in enumerate(refIDs)
    # Load reference data (https://zenodo.org/communities/technion_ae_lab/records?q=&l=list&p=1&s=10)
    data = load(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS20/ID_"*string(refID)*".jld2")
    θ = data["aoa"]
    tipMassConfig = contains(data["tip_weight_config"], "LE") ? :LE : :TE
    t_exp = data["t"]
    U_exp = data["V"]
    # Display configuration
    println("ID $refID configuration: Λ=$(round(Int, rad2deg(Λ))) deg, θ=$(θ) deg, $(data["tip_weight_config"]),")
    println("Runtime of $(round(t_exp[end]-t_exp[1],digits=1)) s, airspeed ranging from $(round(minimum(U_exp),digits=1)) to $(round(maximum(U_exp),digits=1)) m/s")
    # Geometric properties
    nElem,L,chord,referenceAxis,nodalPositions = geometrical_properties_Pazy(Λ=Λ,mesh=mesh)
    # Smooth airspeed signal, if applicable
    U_exp = smoothAirspeedSignal ? moving_average(U_exp,smoothingWindow) : U_exp
    # Airspeed profile (interpolate linearly between experimental time steps, constant outside)
    Uprofile[i] = t -> LinearInterpolations.interpolate(t_exp, U_exp, clamp(t,t_exp[1],t_exp[end]))
    # Find flutter speed of configuration
    URangeFlutter = collect(0.5:0.5:100)
    nModes = 3
    eigenModel,_ = create_Pazy(aeroSolver=aeroSolver,mesh=mesh,sweepAngle=Λ,pitchAngle=deg2rad(θ),tipMassConfig=tipMassConfig,tipRodInertiaCorrection=true,μWing=μWing)
    eigensols = Array{EigenSolution}(undef,length(URangeFlutter))
    x0 = zeros(0)
    for (j,U) in enumerate(URangeFlutter)
        update_v_A!(eigenModel; v_A=[0,U,0])
        problem = create_EigenProblem(model=eigenModel,eigenmethod=:fulldecomp)
        solve!(problem)
        eigensols[j] = problem.eigensolution
        x0 = problem.systemArrays.x
    end
    freqs,damps,_ = mode_tracking(eigensols)
    freqsOsc,dampsOsc,_ = filter_tracked_modes(freqs,damps; fmin=Hz2rads(3.75), fmax=Hz2rads(50), maxTracks=nModes)
    onsetSpeeds,onsetFreqs,onsetModes,offsetSpeeds,offsetFreqs,offsetModes = find_flutter(URangeFlutter,freqsOsc,dampsOsc)
    if !isempty(onsetSpeeds)
        flutterOnsetSpeed = onsetSpeeds[1]
        flutterOnsetFreq = onsetFreqs[1]
        println("Flutter onset: $(round(flutterOnsetSpeed,digits=2)) m/s, $(round(rads2Hz(flutterOnsetFreq),digits=1)) Hz")
    else
        println("Flutter onset: not found")
    end
    if !isempty(offsetSpeeds)
        flutterOffsetSpeed = offsetSpeeds[1]
        flutterOffsetFreq = offsetFreqs[1]
        println("Flutter offset: $(round(flutterOffsetSpeed,digits=2)) m/s, $(round(rads2Hz(flutterOffsetFreq),digits=1)) Hz")
    else
        println("Flutter offset: not found")
    end
    # Plot flutter results
    modeColors = cgrad(:rainbow, nModes, categorical=true)
    ts = 10
    fs = 14
    lw = 2
    gr()
    plt_Vf = plot(ylabel="Frequency [Hz]", xlims=extrema(URangeFlutter), ylims=[0,50], xticks=vcat(0:10:100), tickfont=font(ts), guidefont=font(fs))
    for mode in axes(freqsOsc,1)
        plot!(URangeFlutter, rads2Hz.(freqsOsc[mode,:]), c=modeColors[mode], lw=lw, label=false)
    end
    plt_Vg = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", xlims=extrema(URangeFlutter), ylims=[-0.2,0.1], xticks=vcat(0:10:100), tickfont=font(ts), guidefont=font(fs))
    plot!(URangeFlutter, zeros(length(URangeFlutter)), c=:black, lw=lw, ls=:dash, label=false)
    for mode in axes(dampsOsc,1)
        plot!(URangeFlutter, dampsOsc[mode,:]./freqsOsc[mode,:], c=modeColors[mode], lw=lw, label=false)
    end
    plt_Vgf = plot(plt_Vf,plt_Vg, layout=(2,1))
    display(plt_Vgf)
    # Model for steady problem
    steadyModel,_ = create_Pazy(aeroSolver=aeroSolver,mesh=mesh,sweepAngle=Λ,pitchAngle=deg2rad(θ),airspeed=U_exp[1],tipMassConfig=tipMassConfig,tipRodInertiaCorrection=true,μWing=μWing)
    # Create and solve steady problem for initial airspeed
    steadyProblem = create_SteadyProblem(model=steadyModel)
    solve!(steadyProblem)
    # Model for dynamic problem
    dynamicModel,_ = create_Pazy(aeroSolver=aeroSolver,mesh=mesh,sweepAngle=Λ,pitchAngle=deg2rad(θ),airspeed=Uprofile[i],tipLossFunctionIsAirspeedDependent=true,tipMassConfig=tipMassConfig,tipRodInertiaCorrection=true,μWing=μWing)
    # Create and solve dynamic problem
    dynamicProblem[i] = create_DynamicProblem(model=dynamicModel,integrationMethod=integrationMethod,initialTime=t_exp[1],finalTime=t_exp[end],Δt=Δt,systemSolver=systemSolver,x0=steadyProblem.systemArrays.x,trackingFrequency=trackingFrequency)
    solve!(dynamicProblem[i])
    # Unpack numerical solution
    sol = dynamicProblem[i].solutionHistory
    t[i] = sol.t
    tipAoA[i] = [sol.aeroVariables[j][nElem].unsteadyKin.αₑ for j in eachindex(t[i])]
    tipOOP[i] = -[sol.states[j][nElem].nodal.u2[1] for j in eachindex(t[i])]
    tipTwist[i] = [twist_angle(sol.states[j][nElem].nodal.p2_b) for j in eachindex(t[i])]
    ϵ11Root = [sol.states[j][1].midpoint.γ[1] for j in eachindex(t[i])]
    κ2Root = [sol.states[j][1].midpoint.κ[2] for j in eachindex(t[i])]
    κ3Root = [sol.states[j][1].midpoint.κ[3] for j in eachindex(t[i])]
    rootEpsLE[i] = (ϵ11Root .- κ2Root*x3SG .- κ3Root*x2SG_LE)
    rootEpsTE[i] = (ϵ11Root .- κ2Root*x3SG .- κ3Root*x2SG_TE)
    ϵ11SG = [sol.states[j][SGelem].midpoint.γ[1] for j in eachindex(t[i])]
    κ2SG = [sol.states[j][SGelem].midpoint.κ[2] for j in eachindex(t[i])]
    κ3SG = [sol.states[j][SGelem].midpoint.κ[3] for j in eachindex(t[i])]
    SGEpsLE[i] = (ϵ11SG .- κ2SG*x3SG .- κ3SG*x2SG_LE)
    SGEpsTE[i] = (ϵ11SG .- κ2SG*x3SG .- κ3SG*x2SG_TE)
end

## Plots
# Set paths
relPath = "/dev/Pazy/S20/InertiaRod/outputs/S20MatchingAirspeedSweep"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
showTitle = true
showAirspeedInPlots = false
interactivePlots = false
PSDplot = true
savePlots = true
colors = cgrad(:rainbow, 2, categorical=true)
epsColors = cgrad(:rainbow, 3, categorical=true)
ts = 10
fs = 16
lfs = 8
lw = 2
ms = 6
msw = 1
α = 0.5
windowFFT = 2^14
lsExp = :solid
if interactivePlots
    plotlyjs()
else
    gr()
end 

# Sweep IDs
for (i,refID) in enumerate(refIDs)
    # Case ID
    idString = "ID_" * string(refID) * "_" * string(mesh) * "mesh_" * aeroSolver.name * "_μWing" * string(μWing) * "_ρGenα" * string(ρGenα)
    # Load reference data (https://zenodo.org/communities/technion_ae_lab/records?q=&l=list&p=1&s=10)
    data = load(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS20/ID_"*string(refID)*".jld2")
    Λ = deg2rad(10)
    θ = data["aoa"]
    tipMassConfig = contains(data["tip_weight_config"], "LE") ? :LE : :TE
    t_exp = data["t"]
    t_240Hz_exp = data["t_cam"]
    U_exp = data["V"]
    root_eps_LE_exp = data["strain_LE_root"]
    root_eps_TE_exp = data["strain_TE_root"]
    gauge_eps_LE_exp = data["strain_LE_"*gaugeName]
    gauge_eps_TE_exp = data["strain_TE_"*gaugeName]
    tip_twist_exp = data["twist"]
    tip_disp_exp = (data["tip_disp"][:,3])
    # Geometric properties
    nElem,L,_ = geometrical_properties_Pazy(Λ=Λ,mesh=mesh)
    # Title
    title = showTitle ? 
            "ID $refID, Λ=$(round(Int, rad2deg(Λ))), θ=$(θ), $(data["tip_weight_config"])" : 
            ""
    # Airspeed
    plt_U = plot(xlabel="Time [s]", ylabel="Airspeed [m/s]", title=title, xlims=extrema(t_exp), tickfont=font(ts), guidefont=font(fs))
    if smoothAirspeedSignal
        plot!(t_exp, U_exp, lw=lw, c=:blue, alpha=α, label=false)
    end
    plot!(t[i], Uprofile[i].(t[i]), lw=lw, c=:black, label=false)
    display(plt_U)
    savePlots && savefig(string(absPath,"/S20MatchingAirspeedSweep_",idString,"_U.pdf"))
    # Tip OOP displacement
    plt_tipOOP = plot(xlabel="Time [s]", ylabel="Tip OOP disp. [% semispan]", title=title, xlims=extrema(t_exp), tickfont=font(ts), guidefont=font(fs), legend=:best, legendfontsize=lfs)
    plot!(t[i], tipOOP[i]/L*100, lw=lw, c=epsColors[1], alpha=α, label="AeroBeams")
    plot!(t_240Hz_exp, tip_disp_exp*1e-3*100/L, lw=lw, ls=lsExp, c=epsColors[2], alpha=α, label="Revivo & Raveh (2025)")
    if showAirspeedInPlots
        plt_tipOOP_twin = twinx(plt_tipOOP)
        plot!(plt_tipOOP_twin, ylabel="Airspeed [m/s]", tickfont=font(ts), guidefont=font(fs), xlims=extrema(t_exp))
        plot!(plt_tipOOP_twin, t[i], Uprofile[i].(t[i]), lw=lw, c=:black, label=false)
    end
    display(plt_tipOOP)
    savePlots && savefig(string(absPath,"/S20MatchingAirspeedSweep_",idString,"_tipOOP.pdf"))
    # Tip twist
    plt_tipTwist = plot(xlabel="Time [s]", ylabel="Tip twist [deg]", title=title, xlims=extrema(t_exp), tickfont=font(ts), guidefont=font(fs))
    plot!(t[i], tipTwist[i], lw=lw, c=epsColors[1], alpha=α, label="AeroBeams")
    plot!(t_240Hz_exp, tip_twist_exp, lw=lw, ls=lsExp, c=epsColors[2], alpha=α, label="Revivo & Raveh (2025)")
    if showAirspeedInPlots
        plt_tipTwist_twin = twinx(plt_tipTwist)
        plot!(plt_tipTwist_twin, ylabel="Airspeed [m/s]", tickfont=font(ts), guidefont=font(fs), xlims=extrema(t_exp))
        plot!(plt_tipTwist_twin, t[i], Uprofile[i].(t[i]), lw=lw, c=:black, label=false)
    end
    display(plt_tipTwist)
    savePlots && savefig(string(absPath,"/S20MatchingAirspeedSweep_",idString,"_tipTwist.pdf"))
    # Tip AoA
    plt_tipAOA = plot(xlabel="Time [s]", ylabel="Tip angle of attack [deg]", title=title, xlims=extrema(t_exp), tickfont=font(ts), guidefont=font(fs))
    plot!(t[i], rad2deg.(tipAoA[i]), lw=lw, c=colors[1], label=false)
    plot!([NaN], [NaN], lw=lw, ls=lsExp, c=colors[2], label=false)
    if showAirspeedInPlots
        plt_tipAOA_twin = twinx(plt_tipAOA)
        plot!(plt_tipAOA_twin, ylabel="Airspeed [m/s]", tickfont=font(ts), guidefont=font(fs), xlims=extrema(t_exp))
        plot!(plt_tipAOA_twin, t[i], Uprofile[i].(t[i]), lw=lw, c=:black, label=false)
    end
    display(plt_tipAOA)
    savePlots && savefig(string(absPath,"/S20MatchingAirspeedSweep_",idString,"_tipAoA.pdf"))
    # Root LE strains
    plt_rootEpsLE = plot(xlabel="Time [s]", ylabel="Near root LE strains (\$\\mu\$)", title=title, xlims=extrema(t_exp), tickfont=font(ts), guidefont=font(fs), legend=:best, legendfontsize=lfs)
    plot!(t[i], rootEpsLE[i]*1e6, lw=lw, c=epsColors[1], alpha=α, label="AeroBeams")
    plot!(t_exp, root_eps_LE_exp, lw=lw, ls=lsExp, c=epsColors[2], alpha=α, label="Revivo & Raveh (2025)")
    if showAirspeedInPlots
        plt_rootEpsLE_twin = twinx(plt_rootEpsLE)
        plot!(plt_rootEpsLE_twin, ylabel="Airspeed [m/s]", tickfont=font(ts), guidefont=font(fs), xlims=extrema(t_exp))
        plot!(plt_rootEpsLE_twin, t[i], Uprofile[i].(t[i]), lw=lw, c=:black, label=false)
    end
    display(plt_rootEpsLE)
    savePlots && savefig(string(absPath,"/S20MatchingAirspeedSweep_",idString,"_rootEpsLE.pdf"))
    # Root TE strains
    plt_rootEpsTE = plot(xlabel="Time [s]", ylabel="Near root TE strains (\$\\mu\$)", title=title, xlims=extrema(t_exp), tickfont=font(ts), guidefont=font(fs), legend=:best, legendfontsize=lfs)
    plot!(t[i], rootEpsTE[i]*1e6, lw=lw, c=epsColors[1], alpha=α, label=false)
    plot!(t_exp, root_eps_TE_exp, lw=lw, ls=lsExp, c=epsColors[2], alpha=α, label=false)
    if showAirspeedInPlots
        plt_rootEpsTE_twin = twinx(plt_rootEpsTE)
        plot!(plt_rootEpsTE_twin, ylabel="Airspeed [m/s]", tickfont=font(ts), guidefont=font(fs), xlims=extrema(t_exp))
        plot!(plt_rootEpsTE_twin, t[i], Uprofile[i].(t[i]), lw=lw, c=:black, label=false)
    end
    display(plt_rootEpsTE)
    savePlots && savefig(string(absPath,"/S20MatchingAirspeedSweep_",idString,"_rootEpsTE.pdf"))
    # Selected strain gauge LE strains
    plt_SGEpsLE = plot(xlabel="Time [s]", ylabel="$(uppercasefirst(gaugeName)) LE strains (\$\\mu\$)", title=title, xlims=extrema(t_exp), tickfont=font(ts), guidefont=font(fs), legend=:best, legendfontsize=lfs)
    plot!(t[i], SGEpsLE[i]*1e6, lw=lw, c=epsColors[1], alpha=α, label="AeroBeams")
    plot!(t_exp, gauge_eps_LE_exp, lw=lw, ls=lsExp, c=epsColors[2], alpha=α, label="Revivo & Raveh (2025)")
    if showAirspeedInPlots
        plt_SGEpsLE_twin = twinx(plt_SGEpsLE)
        plot!(plt_SGEpsLE_twin, ylabel="Airspeed [m/s]", tickfont=font(ts), guidefont=font(fs), xlims=extrema(t_exp))
        plot!(plt_SGEpsLE_twin, t[i], Uprofile[i].(t[i]), lw=lw, c=:black, label=false)
    end
    display(plt_SGEpsLE)
    savePlots && savefig(string(absPath,"/S20MatchingAirspeedSweep_",idString,"_",gaugeName,"EpsLE.pdf"))
    # Selected strain gauge TE strains
    plt_SGEpsTE = plot(xlabel="Time [s]", ylabel="$(uppercasefirst(gaugeName)) TE strains (\$\\mu\$)", title=title, xlims=extrema(t_exp), tickfont=font(ts), guidefont=font(fs), legend=:best, legendfontsize=lfs)
    plot!(t[i], SGEpsTE[i]*1e6, lw=lw, c=epsColors[1], alpha=α, label=false)
    plot!(t_exp, gauge_eps_TE_exp, lw=lw, ls=lsExp, c=epsColors[2], alpha=α, label=false)
    if showAirspeedInPlots
        plt_SGEpsTE_twin = twinx(plt_SGEpsTE)
        plot!(plt_SGEpsTE_twin, ylabel="Airspeed [m/s]", tickfont=font(ts), guidefont=font(fs), xlims=extrema(t_exp))
        plot!(plt_SGEpsTE_twin, t[i], Uprofile[i].(t[i]), lw=lw, c=:black, label=false)
    end
    display(plt_SGEpsTE)
    savePlots && savefig(string(absPath,"/S20MatchingAirspeedSweep_",idString,"_",gaugeName,"EpsTE.pdf"))
    # Root TE strains spectrogram
    f_rootEps, t_rootEps, S_rootEps = spectrogram(t[i], rootEpsTE[i]*1e6, window=windowFFT)
    tt = title * ", root TE"
    plt_rootEps_spectrum = heatmap(t_rootEps, f_rootEps, S_rootEps; ylims=[0,100], title=tt, xlabel="Time [s]", ylabel="Frequency [Hz]", c=:rainbow, colorbar_title="Power [dB]", aspect_ratio=:auto, tickfont=font(ts), guidefont=font(fs))
    display(plt_rootEps_spectrum)
    savePlots && savefig(string(absPath,"/S20MatchingAirspeedSweep_",idString,"_rootEpsTESpectrogram.pdf"))
    # Experimental root LE strains spectrogram
    f_eps, t_eps, S_eps = spectrogram(t_exp, root_eps_TE_exp, window=2^12)
    tt = title * ", root TE Exp."
    plt_rootEps_spectrum_exp = heatmap(t_eps, f_eps, S_eps; ylims=[0,100], title=tt, xlabel="Time [s]", ylabel="Frequency [Hz]", c=:rainbow, colorbar_title="Power [dB]", aspect_ratio=:auto, tickfont=font(ts), guidefont=font(fs))
    display(plt_rootEps_spectrum_exp)
    savePlots && savefig(string(absPath,"/S20MatchingAirspeedSweep_",idString,"_rootEpsTESpectrogram_exp.pdf"))
    # Selected gauge TE strains spectrogram
    f_gEps, t_gEps, S_gEps = spectrogram(t[i], SGEpsTE[i]*1e6, window=windowFFT)
    tt = title * ", " * gaugeName
    plt_gEps_spectrum = heatmap(t_gEps, f_gEps, S_gEps; ylims=[0,100], title=tt, xlabel="Time [s]", ylabel="Frequency [Hz]", c=:rainbow, colorbar_title="Power [dB]", aspect_ratio=:auto, tickfont=font(ts), guidefont=font(fs))
    display(plt_gEps_spectrum)
    savePlots && savefig(string(absPath,"/S20MatchingAirspeedSweep_",idString,"_",gaugeName,"EpsTESpectrogram.pdf"))
    # Experimental gauge TE strains spectrogram
    f_gEps, t_gEps, S_gEps = spectrogram(t_exp, gauge_eps_TE_exp, window=2^12)
    tt = title * ", " * gaugeName * " Exp."
    plt_gEps_spectrum_exp = heatmap(t_gEps, f_gEps, S_gEps; ylims=[0,100], title=tt, xlabel="Time [s]", ylabel="Frequency [Hz]", c=:rainbow, colorbar_title="Power [dB]", aspect_ratio=:auto, tickfont=font(ts), guidefont=font(fs))
    display(plt_gEps_spectrum_exp)
    savePlots && savefig(string(absPath,"/S20MatchingAirspeedSweep_",idString,"_",gaugeName,"EpsTESpectrogram_exp.pdf"))
    # PSD of the selected gauge TE strains
    if PSDplot
        minTime = 50
        maxTime = t[i][end]
        timeRangeInd = findfirst(t[i] .>= minTime) : findlast(t[i] .<= maxTime)
        timeRangeIndExp = findfirst(t_exp .>= minTime) : findlast(t_exp .<= maxTime)
        f_sim, _, psd_sim = get_FFT_and_PSD(t[i][timeRangeInd], SGEpsTE[i][timeRangeInd]*1e6)
        f_exp, _, psd_exp = get_FFT_and_PSD(t_exp[timeRangeIndExp], gauge_eps_TE_exp[timeRangeIndExp])
        plt_gEps_psd = plot(title=tt, xlabel="Frequency [Hz]", ylabel="Strain PSD [µε^2/Hz] in dB", xlims=[0,100], ylims=[-40,80], tickfont=font(ts), guidefont=font(fs))
        plot!(f_sim, 10*log10.(psd_sim); lw=lw, c=epsColors[1], label="AeroBeams")
        plot!(f_exp, 10*log10.(psd_exp); lw=lw, c=epsColors[2], label="Exp.")
        display(plt_gEps_psd)
        savePlots && savefig(string(absPath,"/S20MatchingAirspeedSweep_",idString,"_",gaugeName,"EpsTEPSD.pdf"))
    end
end

println("Finished S20MatchingAirspeedSweep.jl")