# This script is used to calibrate the T1 mode by adjusting the factors applied when applyTipInertiaCorrection=true, correcting Iyy[end], in inertia_matrices_Pazy for the S20 Pazy model
using AeroBeams

# Sweep angle [rad]
Λ = deg2rad(20)

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
nastranS20LEF = [4.16, 28.18, 34.45, 83.84, 95.02, 131.09, 161.57] # LE, flexible wingtip rod
nastranS20TEF = [4.14, 25.53, 38.32, 80.95, 85.38, 125.56, 168.75] # TE, flexible wingtip rod

# Compute relative differences [%]
ϵ_gS20LE = 100*(freqsOsc[1] ./ nastranS20LEF .- 1)
ϵ_gS20TE = 100*(freqsOsc[2] ./ nastranS20TEF .- 1)

# Mode labels
modeLabels = ["OOP1", "OOP2", "T1", "OOP3", "IP1", "T2", "OOP4"]

# Pretty tables for modal frequency comparison
dataS20LE = hcat(
    [modeLabels[i] for i in 1:nModes], 
    freqsOsc[1],
    nastranS20LEF,
    ϵ_gS20LE,
    )
dataS20TE = hcat(
    [modeLabels[i] for i in 1:nModes], 
    freqsOsc[2],
    nastranS20TEF,
    ϵ_gS20TE,
)

tableS20LE = pretty_table(
    String,
    dataS20LE;
    column_labels    = ["Mode", 
                        "AeroBeams [Hz]", 
                        "Nastran [Hz]",
                        "Δ Nastran [%]"],
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
                        "Nastran [Hz]",
                        "Δ Nastran [%]"],
    title            = "Modal Frequency Comparison - S20 TE",
    formatters       = [
    (v, i, j) -> (j in (2, 3) && v isa Number) ? round(v, digits=2) : v,
    (v, i, j) -> (j in (4) && v isa Number) ? string(round(v, digits=2), " %") : v,
    ],
)
print(tableS20LE)
print(tableS20TE)

println("Finished S20TipInertiaCalibration.jl")