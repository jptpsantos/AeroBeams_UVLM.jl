using AeroBeams

# Meshes
meshes = [:nominal, :refined]

# Number of elements of the wingtip rod beam
nElementsWingtipRod = 16

# Flag to model the wingtip rod as a beam (true) or with its inertia only (false)
beamWingtipRod = true

# Number of elements of the wingtip rod beam (used only if beamWingtipRod is true)
nElementsWingtipRod = 16

# Number of modes to track
nModes = 6

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :krylovkit
nev = 60

# Initialize outputs
problem = Array{EigenProblem}(undef,length(meshes))
eigensols = Array{EigenSolution}(undef,length(meshes))
freqsOsc = Array{Matrix{Float64}}(undef,length(meshes))

# Loop tip mass configurations
for (i, mesh) in enumerate(meshes)
    display("Solving for mesh = $mesh")
    # Model
    model,_ = create_PazyDelft(aeroSolver=QuasiSteady(),mesh=mesh,beamWingtipRod=beamWingtipRod,nElementsWingtipRod=nElementsWingtipRod,tipRodInertiaCorrection=!beamWingtipRod)
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

# Load GVT reference data (from Dussler's IFASD 2024 paper: Effect of Gust Models on the Response Prediction of a Very Flexible Wind Tunnel Wing Model)
GVT = [3.21, 22.47, 29.49, 64.99, missing, 119.70]

# Compute relative differences [%]
ϵ_nom = 100*(freqsOsc[1] ./ GVT .- 1)
ϵ_ref = 100*(freqsOsc[2] ./ GVT .- 1)

# Mode labels
modeLabels = ["OOP1", "OOP2", "T1", "OOP3", "IP1", "T2"]

# Pretty tables for modal frequency comparison
dataNom = hcat(
    [modeLabels[i] for i in 1:nModes], 
    freqsOsc[1],
    GVT,
    ϵ_nom,
    )
dataRef = hcat(
    [modeLabels[i] for i in 1:nModes], 
    freqsOsc[2],
    GVT,
    ϵ_ref,
    )

tableNom = pretty_table(
    String,
    dataNom;
    column_labels    = ["Mode", 
                        "AeroBeams [Hz]", 
                        "GVT [Hz]",
                        "Δ GVT [%]"],
    title            = "Modal Frequency Comparison - Nominal mesh",
    formatters       = [
    (v, i, j) -> (j in (2, 3) && v isa Number) ? round(v, digits=2) : v,
    (v, i, j) -> (j in (4) && v isa Number) ? string(round(v, digits=2), " %") : v,
    ],
)
tableRef = pretty_table(
    String,
    dataRef;
    column_labels    = ["Mode", 
                        "AeroBeams [Hz]", 
                        "GVT [Hz]",
                        "Δ GVT [%]"],
    title            = "Modal Frequency Comparison - Refined mesh",
    formatters       = [
    (v, i, j) -> (j in (2, 3) && v isa Number) ? round(v, digits=2) : v,
    (v, i, j) -> (j in (4) && v isa Number) ? string(round(v, digits=2), " %") : v,
    ],
)
print(tableNom)
print(tableRef)

## Plots
# Set paths
relPath = "/dev/PazyDelft/outputs/PazyDelftModalGVT"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Name of the wingtip rod model, for the outputs
rodModel = beamWingtipRod ? "BeamRod" : "InertiaRod"

# Plot mode shapes
for (i,str) in enumerate(string.(meshes))
    plt_swept_modeShapes = plot_mode_shapes(problem[i],nModes=nModes,view=(30,30),modeLabels=modeLabels,save=true,savePath=string(relPath,"/PazyDelft$(rodModel)ModalGVT_modeShapes_$(str)mesh.pdf"))
    display(plt_swept_modeShapes)
end

println("Finished PazyDelftModalGVT.jl")