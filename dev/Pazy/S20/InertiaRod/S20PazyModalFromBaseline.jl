# Results should match S20PazyModal.jl when matchGVT is false (builds the TE and LE configurations from the baseline)
using AeroBeams

# Sweep angle [rad]
Λ = deg2rad(20)

# Flag to model the GVT data
matchGVT = false

# Flag for ad hoc correction of torsional inertia of wingtip rod
sweptTipMassInertiaCorrection = matchGVT

# Gravity
applyGravity = matchGVT
g = applyGravity ? standard_atmosphere(0).g : 0.0

# Mesh
mesh = :refined

# Tip mass configurations: baseline model + point masses removed from the LE/TE beam mass models
tipMasses, ηtipMasses, tipMassesInertias, tipMassesElementIDs = get_sweptPazy_tip_inertias_from_baseline(; Λ=Λ, mesh=mesh, configs=(:baseline, :LE, :TE), withFPSMass=true, withAdapterMass=true, withCalibrationMass=true, sweptTipMassInertiaCorrection=sweptTipMassInertiaCorrection)

# Number of modes to track
nModes = 7

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :fulldecomp
nev = 50

# Initialize outputs
problem = Array{EigenProblem}(undef,length(tipMasses))
eigensols = Array{EigenSolution}(undef,length(tipMasses))
freqsOsc = Array{Matrix{Float64}}(undef,length(tipMasses))

# Loop tip mass configurations
for (i,m,η,I,elementIDs) in zip(1:length(tipMasses), tipMasses, ηtipMasses, tipMassesInertias, tipMassesElementIDs)
    display("Solving for configuration = $i")
    # Model
    model,_ = create_Pazy(sweepAngle=Λ,mesh=mesh,upright=true,gravityVector=[0,0,-g],sweptTipMassConfig=:baseline,tipMasses=m,ηtipMasses=η,tipMassesInertias=I,tipMassesElementIDs=elementIDs)
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

# Load GVT reference data (https://doi.org/10.2514/1.J066240)
GVTS20CL = [4.5, 29.5, 65, 86.1, missing, missing, missing]
GVTS20LE = [4.2, 28.4, 34.0, 84.4, missing, missing, missing]
GVTS20TE = [4.1, 25.1, 38.1, 83.8, missing, missing, missing]

# Load NASTRAN reference data (derived from https://doi.org/10.2514/1.J066240)
nastranS20CLF = [4.53, 29.28, 57.05, 84.80, 103.38, 161.62, 183.58] # clean, flexible wingtip rod
nastranS20CLS = [4.53, 29.29, 57.48, 84.95, 103.44, 163.34, 190.02] # clean, stiffened wingtip rod
nastranS20LEF = [4.16, 28.18, 34.45, 83.84, 95.02, 131.09, 161.57] # LE, flexible wingtip rod
nastranS20LES = [4.17, 28.43, 36.40, 84.27, 94.91, 163.25, 168.45] # LE, stiffened wingtip rod
nastranS20TEF = [4.14, 25.53, 38.32, 80.95, 85.38, 125.56, 168.75] # TE, flexible wingtip rod
nastranS20TES = [4.14, 26.00, 40.00, 82.50, 94.34, 155.14, 177.14] # TE, stiffened wingtip rod

# Compute relative differences [%]
ϵ_gS20CL = 100*(freqsOsc[1] ./ GVTS20CL .- 1)
ϵ_gS20LE = 100*(freqsOsc[2] ./ GVTS20LE .- 1)
ϵ_gS20TE = 100*(freqsOsc[3] ./ GVTS20TE .- 1)
ϵ_nS20CLS = 100*(freqsOsc[1] ./ nastranS20CLS .- 1)
ϵ_nS20CLF = 100*(freqsOsc[1] ./ nastranS20CLF .- 1)
ϵ_nS20LEF = 100*(freqsOsc[2] ./ nastranS20LEF .- 1)
ϵ_nS20LES = 100*(freqsOsc[2] ./ nastranS20LES .- 1)
ϵ_nS20TEF = 100*(freqsOsc[3] ./ nastranS20TEF .- 1)
ϵ_nS20TES = 100*(freqsOsc[3] ./ nastranS20TES .- 1)

# Mode labels
modeLabels = ["OOP1", "OOP2", "T1", "OOP3", "IP1", "T2", "OOP4"]

# Pretty tables for modal frequency comparison
dataS20CL = hcat(
    [modeLabels[i] for i in eachindex(modeLabels)], 
    freqsOsc[1],
    GVTS20CL,
    nastranS20CLF,
    nastranS20CLS,
    ϵ_gS20CL,
    ϵ_nS20CLF,
    ϵ_nS20CLS,
    )
dataS20LE = hcat(
    [modeLabels[i] for i in eachindex(modeLabels)], 
    freqsOsc[2],
    GVTS20LE,
    nastranS20LEF,
    nastranS20LES,
    ϵ_gS20LE,
    ϵ_nS20LEF,
    ϵ_nS20LES,
    )
dataS20TE = hcat(
    [modeLabels[i] for i in eachindex(modeLabels)], 
    freqsOsc[3],
    GVTS20TE,
    nastranS20TEF,
    nastranS20TES,
    ϵ_gS20TE,
    ϵ_nS20TEF,
    ϵ_nS20TES,
)

tableS20CL = pretty_table(
    String,
    dataS20CL;
    column_labels    = ["Mode", 
                        "AeroBeams [Hz]", 
                        "GVT [Hz]",
                        "NASTRAN (flex) [Hz]", 
                        "NASTRAN (stiff) [Hz]", 
                        "Δ GVT [%]",
                        "Δ NASTRAN (flex) [%]",
                        "Δ NASTRAN (stiff) [%]"],
    title            = "Modal Frequency Comparison - S20 baseline (clean)",
    formatters       = [
    (v, i, j) -> (j in (2, 3, 4, 5) && v isa Number) ? round(v, digits=2) : v,
    (v, i, j) -> (j in (6, 7, 8) && v isa Number) ? string(round(v, digits=2), " %") : v,
    ],
)
tableS20LE = pretty_table(
    String,
    dataS20LE;
    column_labels    = ["Mode", 
                        "AeroBeams [Hz]", 
                        "GVT [Hz]",
                        "NASTRAN (flex) [Hz]", 
                        "NASTRAN (stiff) [Hz]",
                        "Δ GVT [%]",
                        "Δ NASTRAN (flex) [%]",
                        "Δ NASTRAN (stiff) [%]"],
    title            = "Modal Frequency Comparison - S20 LE",
    formatters       = [
    (v, i, j) -> (j in (2, 3, 4, 5) && v isa Number) ? round(v, digits=2) : v,
    (v, i, j) -> (j in (6, 7, 8) && v isa Number) ? string(round(v, digits=2), " %") : v,
    ],
)
tableS20TE = pretty_table(
    String,
    dataS20TE;
    column_labels    = ["Mode", 
                        "AeroBeams [Hz]", 
                        "GVT [Hz]",
                        "NASTRAN (flex) [Hz]",
                        "NASTRAN (stiff) [Hz]",
                        "Δ GVT [%]",
                        "Δ NASTRAN (flex) [%]",
                        "Δ NASTRAN (stiff) [%]"],
    title            = "Modal Frequency Comparison - S20 TE",
    formatters       = [
    (v, i, j) -> (j in (2, 3, 4, 5) && v isa Number) ? round(v, digits=2) : v,
    (v, i, j) -> (j in (6, 7, 8) && v isa Number) ? string(round(v, digits=2), " %") : v,
    ],
)
print(tableS20CL)
print(tableS20LE)
print(tableS20TE)

## Plots
# Set paths
relPath = "/dev/Pazy/S20/InertiaRod/outputs/S20PazyModalFromBaseline"
absPath = string(pwd(), relPath)
mkpath(absPath)

# Plot mode shapes
for i in eachindex(tipMasses)
    plt_swept_modeShapes = plot_mode_shapes(problem[i],nModes=nModes,view=(30,30),modeLabels=modeLabels,save=true,savePath=string(relPath,"/S20PazyModal2_modeShapes_$(i).pdf"))
    display(plt_swept_modeShapes)
end

println("Finished S20PazyModalFromBaseline.jl")