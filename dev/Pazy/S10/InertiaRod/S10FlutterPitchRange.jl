using AeroBeams

# Sweep angle [rad]
Λ = deg2rad(10)

# Tip mass configurations
tipMassConfigs = [:LE, :TE]

# Root pitch angle range [rad]
θRange = deg2rad.([0,3,7,10])

# Airspeed range [m/s]
URange = vcat(0:1:90)

# Mesh
mesh = :nominal

# Aerodynamic solver
aeroSolver = Indicial()

# Flag for ad hoc correction of torsional inertia of wingtip rod
tipRodInertiaCorrection = true

# Number of modes to track
nModes = 3

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = aeroSolver isa DynamicStallSolver ? :fulldecomp : :krylovkit
nev = 50

# Geometric properties
nElem,L,chord,referenceAxis,nodalPositions = geometrical_properties_Pazy(Λ=Λ,mesh=mesh)

# Initialize outputs
eigensols = Array{EigenSolution}(undef,length(tipMassConfigs),length(θRange),length(URange))
tipOOP = Array{Float64}(undef,length(tipMassConfigs),length(θRange),length(URange))
freqsOsc = Array{Matrix{Float64}}(undef,length(tipMassConfigs),length(θRange))
dampsOsc = Array{Matrix{Float64}}(undef,length(tipMassConfigs),length(θRange))

# Sweep tip weight configurations
for (c,config) in enumerate(tipMassConfigs)
    # Sweep pitch angle
    for (i,θ) in enumerate(θRange)
        # Initialize guess solution
        x0 = zeros(0)
        # Sweep airspeed
        for (j,U) in enumerate(URange)
            # Display progress
            println("Solving for tip mass at $config, θ = $(round(rad2deg(θ),digits=1)) deg, U = $U m/s")
            # Model
            model,_ = create_Pazy(aeroSolver=aeroSolver,mesh=mesh,sweepAngle=Λ,pitchAngle=θ,airspeed=U,tipMassConfig=config,tipRodInertiaCorrection=tipRodInertiaCorrection)
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
    end
end

## Plots
# Set paths
relPath = "/dev/Pazy/S10/InertiaRod/outputs/S10FlutterPitchRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
modeColors = cgrad(:rainbow, nModes, categorical=true)
ts = 10
fs = 14
lfs = 8
lw = 2
ms = 4
msw = 1
gr()

# V-g-f
for (c,config) in enumerate(tipMassConfigs)
    plt_Vf = plot(ylabel="Frequency [Hz]", xlims=[0,URange[end]], ylims=[0,40], xticks=vcat(0:10:URange[end]), tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
    plot!([NaN], [NaN], lw=lw, c=:black, label="AeroBeams")
    for (i,θ) in enumerate(θRange)
        for mode in axes(freqsOsc[c,i],1)
            plot!(URange, rads2Hz.(freqsOsc[c,i][mode,:]), c=modeColors[mode], lw=lw, alpha=0.1+0.9*i/length(θRange), label=false)
        end
    end
    plt_Vg = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", xlims=[0,URange[end]], ylims=[-0.2,0.1], xticks=vcat(0:10:URange[end]), tickfont=font(ts), guidefont=font(fs), legend=:topleft)
    plot!(URange, zeros(length(URange)), c=:black, lw=lw, ls=:dash, label=false)
    for (i,θ) in enumerate(θRange)
        for mode in axes(dampsOsc[c,i],1)
            plot!(URange, dampsOsc[c,i][mode,:]./freqsOsc[c,i][mode,:], c=modeColors[mode], lw=lw, alpha=0.1+0.9*i/length(θRange), label=false)
        end
    end
    plt_Vgf = plot(plt_Vf,plt_Vg, layout=(2,1))
    display(plt_Vgf)
    savefig(plt_Vgf,string(absPath,"/S10FlutterPitchRange_Vgf_",config,".pdf"))
end

println("Finished S10FlutterPitchRange.jl")