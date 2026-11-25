using AeroBeams

# Root pitch angle range [rad]
θRange = deg2rad.([0,3,5,7])

# Airspeed range [m/s]
URange = vcat(0:1:100)

# Mesh
mesh = :refined

# Flag to model the wingtip rod as a beam (true) or with its inertia only (false)
beamWingtipRod = true

# Number of elements of the wingtip rod beam (used only if beamWingtipRod is true)
nElementsWingtipRod = 16

# Structural damping coefficients for the wing and wingtip rod
μWing = 5e-6
μRod = 5e-5

# Aerodynamic solver
aeroSolver = Indicial()

# Number of modes to track
nModes = 3

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = aeroSolver isa DynamicStallSolver ? :fulldecomp : :krylovkit
nev = 65

# Geometric properties
nElem,L,chord,referenceAxis,nodalPositions = geometrical_properties_Pazy(mesh=mesh)

# Initialize outputs
eigensols = Array{EigenSolution}(undef,length(θRange),length(URange))
tipOOP = Array{Float64}(undef,length(θRange),length(URange))
freqsOsc = Array{Matrix{Float64}}(undef,length(θRange))
dampsOsc = Array{Matrix{Float64}}(undef,length(θRange))

# Sweep pitch angle
for (i,θ) in enumerate(θRange)
    # Initialize guess solution
    x0 = zeros(0)
    # Sweep airspeed
    for (j,U) in enumerate(URange)
        # Display progress
        println("Solving for θ = $(round(rad2deg(θ),digits=1)) deg, U = $U m/s")
        # Model
        model,_ = create_Pazy(aeroSolver=aeroSolver,mesh=mesh,pitchAngle=θ,airspeed=U,beamWingtipRod=beamWingtipRod,nElementsWingtipRod=nElementsWingtipRod,μWing=μWing,μRod=μRod)
        # Create and solve problem
        problem = create_EigenProblem(model=model,eigenmethod=eigenmethod,nev=nev,x0=x0)
        solve!(problem)
        # Eigensolution
        eigensols[i,j] = problem.eigensolution
        # Tip OOP deflection
        tipOOP[i,j] = -problem.loadstepHistory.states[end][nElem].nodal.u2[1]
        # Set initial guess solution for next iteration
        x0 = problem.systemArrays.x
    end
    # Apply mode tracking
    freqs,damps,_ = mode_tracking(eigensols[i,:])
    # Get frequencies and dampings of oscillatory modes
    freqsOsc[i],dampsOsc[i],_ = filter_tracked_modes(freqs,damps; fmin=Hz2rads(1e-1), fmax=Hz2rads(50), maxTracks=nModes)
end

## Plots
# Set paths
relPath = "/dev/Pazy/S0/outputs/S0FlutterPitchRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Name of the wingtip rod model, for the outputs
rodModel = beamWingtipRod ? "BeamRod" : "InertiaRod"

# Plot configurations
modeColors = cgrad(:rainbow, nModes, categorical=true)
modeLabels = ["OOP1", "OOP2", "T1"]
ts = 10
fs = 14
lfs = 8
lw = 2
ms = 4
msw = 1
gr()

# V-g-f
plt_Vf = plot(ylabel="Frequency [Hz]", xlims=[0,URange[end]], ylims=[0,45], xticks=vcat(0:10:URange[end]), tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, legend=(0.125,0.5), right_margin=5mm)
for (i,θ) in enumerate(θRange)
    plot!(plt_Vf, [NaN], [NaN], c=:black, alpha=0.1+0.9*i/length(θRange), lw=lw, label="\$\\theta=$(round(rad2deg(θ),digits=1))^\\circ\$")
end
plt_Vg = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", xlims=[0,URange[end]], ylims=[-0.1,0.05], xticks=vcat(0:10:URange[end]), tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, legend=(0.125,0.97), right_margin=5mm)
plot!(plt_Vg, URange, zeros(length(URange)), c=:black, lw=lw, ls=:dash, label=false)
for (i,θ) in enumerate(θRange)
    for mode in axes(freqsOsc[i],1)
        plot!(plt_Vf, URange, rads2Hz.(freqsOsc[i][mode,:]), c=modeColors[mode], lw=lw, alpha=0.1+0.9*i/length(θRange), label=false)
    end
    for mode in axes(dampsOsc[i],1)
        plot!(plt_Vg, URange, dampsOsc[i][mode,:]./freqsOsc[i][mode,:], c=modeColors[mode], lw=lw, alpha=0.1+0.9*i/length(θRange), label=i==length(θRange) ? modeLabels[mode] : false)
    end
end
plt_Vgf = plot(plt_Vf,plt_Vg, layout=(2,1))
display(plt_Vgf)
savefig(plt_Vgf,string(absPath,"/S0$(rodModel)FlutterPitchRange_Vgf.pdf"))

println("Finished S0FlutterPitchRange.jl")