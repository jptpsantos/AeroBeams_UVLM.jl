using AeroBeams

# Meshes
meshes = [:nominal, :refined]

# Number of modes to track
nModes = 4

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
    model,_ = create_Pazy(aeroSolver=QuasiSteady(),mesh=mesh)
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

# Load GVT reference data 
GVT = [4.26, 28.50, 42.00, 81.50]

# Compute relative differences [%]
ϵ_nom = 100*(freqsOsc[1] ./ GVT .- 1)
ϵ_ref = 100*(freqsOsc[2] ./ GVT .- 1)

# Mode labels
modeLabels = ["OOP1", "OOP2", "T1", "OOP3"]

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
relPath = "/dev/Pazy/S0/InertiaRod/outputs/PazyWingModalGVT"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot mode shapes
for (i,str) in enumerate(string.(meshes))
    plt_swept_modeShapes = plot_mode_shapes(problem[i],nModes=nModes,view=(30,30),modeLabels=modeLabels,save=true,savePath=string(relPath,"/PazyWingModalGVT_modeShapes_$(str)mesh.pdf"))
    display(plt_swept_modeShapes)
end

println("Finished PazyWingModalGVT.jl")