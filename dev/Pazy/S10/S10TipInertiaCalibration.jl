# This script is used to calibrate the T1 mode by adjusting the factors applied when applyTipInertiaCorrection=true, correcting Iyy[end], in inertia_matrices_Pazy for the S10 Pazy model
using AeroBeams

# Sweep angle [rad]
Λ = deg2rad(10)

# Tip mass configurations
tipMassConfigs = [:LE, :TE]

# Mesh
mesh = :refined

# Number of modes to track
nModes = 7

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
    model,_ = create_Pazy(aeroSolver=QuasiSteady(),sweepAngle=Λ,mesh=mesh,tipMassConfig=tipMassConfig,tipRodInertiaCorrection=true)
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

# Load NASTRAN reference data (derived from https://doi.org/10.2514/1.J066240)
nastranS10LEF = [4.35, 28.84, 37.83, 85.32, 103.27, 167.41, 182.49] # LE, flexible wingtip rod
nastranS10TEF = [4.31, 27.05, 38.71, 84.81, 104.28, 170.09, 194.99] # TE, flexible wingtip rod

# Compute relative differences [%]
ϵ_gS10LE = 100*(freqsOsc[1] ./ nastranS10LEF .- 1)
ϵ_gS10TE = 100*(freqsOsc[2] ./ nastranS10TEF .- 1)

# Mode labels
modeLabels = ["OOP1", "OOP2", "T1", "OOP3", "IP1", "T2", "OOP4"]

# Pretty tables for modal frequency comparison
dataS10LE = hcat(
    [modeLabels[i] for i in 1:nModes], 
    freqsOsc[1],
    nastranS10LEF,
    ϵ_gS10LE,
    )
dataS10TE = hcat(
    [modeLabels[i] for i in 1:nModes], 
    freqsOsc[2],
    nastranS10TEF,
    ϵ_gS10TE,
)

tableS10LE = pretty_table(
    String,
    dataS10LE;
    column_labels    = ["Mode", 
                        "AeroBeams [Hz]", 
                        "Nastran [Hz]",
                        "Δ Nastran [%]"],
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
                        "Nastran [Hz]",
                        "Δ Nastran [%]"],
    title            = "Modal Frequency Comparison - S10 TE",
    formatters       = [
    (v, i, j) -> (j in (2, 3) && v isa Number) ? round(v, digits=2) : v,
    (v, i, j) -> (j in (4) && v isa Number) ? string(round(v, digits=2), " %") : v,
    ],
)
print(tableS10LE)
print(tableS10TE)

println("Finished S10TipInertiaCalibration.jl")