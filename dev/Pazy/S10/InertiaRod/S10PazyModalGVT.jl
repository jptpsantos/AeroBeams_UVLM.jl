using AeroBeams

# Sweep angle [rad]
Λ = deg2rad(10)

# Tip mass configurations
tipMassConfigs = [:baseline, :LE, :TE]

# Mesh
mesh = :refined

# Flag for ad hoc correction of torsional inertia of wingtip rod
tipRodInertiaCorrection = true

# Number of modes to track
nModes = 8

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
    model,_ = create_Pazy(aeroSolver=QuasiSteady(),sweepAngle=Λ,mesh=mesh,tipMassConfig=tipMassConfig,tipRodInertiaCorrection=tipRodInertiaCorrection)
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
GVTS10CL = [4.5, 30.2, 59.0, 87.1, missing, missing, missing, missing]
GVTS10LE = [4.3, 28.9, 37.5, 86.2, missing, missing, 170.2, missing]
GVTS10TE = [4.3, 27.4, 38.7, 85.2, missing, missing, 171.4, missing]

# Compute relative differences [%]
ϵ_gS10CL = 100*(freqsOsc[1] ./ GVTS10CL .- 1)
ϵ_gS10LE = 100*(freqsOsc[2] ./ GVTS10LE .- 1)
ϵ_gS10TE = 100*(freqsOsc[3] ./ GVTS10TE .- 1)

# Mode labels
modeLabelsBaseline = ["OOP1", "OOP2", "T1", "OOP3", "IP1", "OOP4", "T2", "OOP5"]
modeLabelsLETE = ["OOP1", "OOP2", "T1", "OOP3", "IP1", "T2", "OOP4", "OOP5"]

# Pretty tables for modal frequency comparison
dataS10CL = hcat(
    [modeLabelsBaseline[i] for i in 1:nModes], 
    freqsOsc[1],
    GVTS10CL,
    ϵ_gS10CL,
    )
dataS10LE = hcat(
    [modeLabelsLETE[i] for i in 1:nModes], 
    freqsOsc[2],
    GVTS10LE,
    ϵ_gS10LE,
    )
dataS10TE = hcat(
    [modeLabelsLETE[i] for i in 1:nModes], 
    freqsOsc[3],
    GVTS10TE,
    ϵ_gS10TE,
)

tableS10CL = pretty_table(
    String,
    dataS10CL;
    column_labels    = ["Mode", 
                        "AeroBeams [Hz]", 
                        "GVT [Hz]",
                        "Δ GVT [%]"],
    title            = "Modal Frequency Comparison - S10 baseline (clean)",
    formatters       = [
    (v, i, j) -> (j in (2, 3) && v isa Number) ? round(v, digits=2) : v,
    (v, i, j) -> (j in (4) && v isa Number) ? string(round(v, digits=2), " %") : v,
    ],
)
tableS10LE = pretty_table(
    String,
    dataS10LE;
    column_labels    = ["Mode", 
                        "AeroBeams [Hz]", 
                        "GVT [Hz]",
                        "Δ GVT [%]"],
    title            = "Modal Frequency Comparison - S10 LE",
    formatters       = [
    (v, i, j) -> (j in (2, 3) && v isa Number) ? round(v, digits=2) : v,
    (v, i, j) -> (j in (4) && v isa Number) ? string(round(v, digits=2), " %") : v,
    ],
)
tableS10TE = pretty_table(
    String,
    dataS10TE;
    column_labels    = ["Mode", 
                        "AeroBeams [Hz]", 
                        "GVT [Hz]",
                        "Δ GVT [%]"],
    title            = "Modal Frequency Comparison - S10 TE",
    formatters       = [
    (v, i, j) -> (j in (2, 3) && v isa Number) ? round(v, digits=2) : v,
    (v, i, j) -> (j in (4) && v isa Number) ? string(round(v, digits=2), " %") : v,
    ],
)
print(tableS10CL)
print(tableS10LE)
print(tableS10TE)

## Plots
# Set paths
relPath = "/dev/Pazy/S10/InertiaRod/outputs/S10PazyModalGVT"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot mode shapes
for (i,str) in enumerate(string.(tipMassConfigs))
    modeLabels = i == 1 ? modeLabelsBaseline : modeLabelsLETE
    plt_swept_modeShapes = plot_mode_shapes(problem[i],nModes=nModes,view=(30,30),modeLabels=modeLabels,save=true,savePath=string(relPath,"/S10PazyModalGVT_modeShapes_$(str).pdf"))
    display(plt_swept_modeShapes)
end

println("Finished S10PazyModalGVT.jl")