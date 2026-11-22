using AeroBeams

# Sweep angle [rad]
Λ = deg2rad(20)

# Tip mass configurations
tipMassConfigs = [:baseline, :LE, :TE]

# Mesh
mesh = :refined

# Number of elements of the wingtip rod beam
nElementsWingtipRod = 12

# Number of modes to track
nModes = 4

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :krylovkit
nev = mesh == :nominal ? 45 : 60

# Initialize outputs
problem = Array{EigenProblem}(undef,length(tipMassConfigs))
eigensols = Array{EigenSolution}(undef,length(tipMassConfigs))
freqsOsc = Array{Matrix{Float64}}(undef,length(tipMassConfigs))

# Loop tip mass configurations
for (i, tipMassConfig) in enumerate(tipMassConfigs)
    display("Solving for tipMassConfig = $tipMassConfig")
    # Model
    model,_ = create_Pazy(aeroSolver=QuasiSteady(),sweepAngle=Λ,mesh=mesh,tipMassConfig=tipMassConfig,beamWingtipRod=true,nElementsWingtipRod=nElementsWingtipRod)
    # Create and solve problem
    problem[i] = create_EigenProblem(model=model,eigenmethod=eigenmethod,nev=nev)
    solve!(problem[i])
    # Eigensolution
    eigensols[i] = problem[i].eigensolution
    # Modal frequencies
    freqs = imag.(eigensols[i].λ)
    # Oscillatory frequencies in Hz
    freqsOsc[i],_ = filter_tracked_modes(rads2Hz.(freqs); fmin=eps(), maxTracks=nModes)
end

# Load GVT reference data (derived from https://doi.org/10.2514/1.J066240)
GVTS20CL = [4.5, 29.5, 65, 86.1]
GVTS20LE = [4.2, 28.4, 34.0, 84.4]
GVTS20TE = [4.1, 25.1, 38.1, 83.8]

# Compute relative differences [%]
ϵ_gS20CL = 100*(freqsOsc[1] ./ GVTS20CL .- 1)
ϵ_gS20LE = 100*(freqsOsc[2] ./ GVTS20LE .- 1)
ϵ_gS20TE = 100*(freqsOsc[3] ./ GVTS20TE .- 1)

# Mode labels
modeLabels = ["OOP1", "OOP2", "T1", "OOP3"]
modeLabels = ["OOP1", "OOP2", "T1", "OOP3"]

# Pretty tables for modal frequency comparison
dataS20CL = hcat(
    [modeLabels[i] for i in 1:nModes], 
    freqsOsc[1],
    GVTS20CL,
    ϵ_gS20CL,
    )
dataS20LE = hcat(
    [modeLabels[i] for i in 1:nModes], 
    freqsOsc[2],
    GVTS20LE,
    ϵ_gS20LE,
    )
dataS20TE = hcat(
    [modeLabels[i] for i in 1:nModes], 
    freqsOsc[3],
    GVTS20TE,
    ϵ_gS20TE,
)

tableS20CL = pretty_table(
    String,
    dataS20CL;
    column_labels    = ["Mode", 
                        "AeroBeams [Hz]", 
                        "GVT [Hz]",
                        "Δ GVT [%]"],
    title            = "Modal Frequency Comparison - S20 baseline (clean)",
    formatters       = [
    (v, i, j) -> (j in (2, 3) && v isa Number) ? round(v, digits=2) : v,
    (v, i, j) -> (j in (4) && v isa Number) ? string(round(v, digits=2), " %") : v,
    ],
)
tableS20LE = pretty_table(
    String,
    dataS20LE;
    column_labels    = ["Mode", 
                        "AeroBeams [Hz]", 
                        "GVT [Hz]",
                        "Δ GVT [%]"],
    title            = "Modal Frequency Comparison - S20 LE",
    formatters       = [
    (v, i, j) -> (j in (2, 3) && v isa Number) ? round(v, digits=2) : v,
    (v, i, j) -> (j in (4) && v isa Number) ? string(round(v, digits=2), " %") : v,
    ],
)
tableS20TE = pretty_table(
    String,
    dataS20TE;
    column_labels    = ["Mode", 
                        "AeroBeams [Hz]", 
                        "GVT [Hz]",
                        "Δ GVT [%]"],
    title            = "Modal Frequency Comparison - S20 TE",
    formatters       = [
    (v, i, j) -> (j in (2, 3) && v isa Number) ? round(v, digits=2) : v,
    (v, i, j) -> (j in (4) && v isa Number) ? string(round(v, digits=2), " %") : v,
    ],
)
print(tableS20CL)
print(tableS20LE)
print(tableS20TE)

## Plots
# Set paths
relPath = "/dev/Pazy/S20/BeamRod/outputs/S20BeamRodModalGVT"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot mode shapes
for (i,str) in enumerate(string.(tipMassConfigs))
    plt_swept_modeShapes = plot_mode_shapes(problem[i],nModes=nModes,view=(30,30),modeLabels=modeLabels,save=true,savePath=string(relPath,"/S20BeamRodModalGVT_modeShapes_$(str).pdf"))
    display(plt_swept_modeShapes)
end

println("Finished S20BeamRodModalGVT.jl")