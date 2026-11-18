using AeroBeams

# Pitch angle range [rad]
θRange = unique([vcat(0:0.5:1)...,vcat(1:1:7)...])

# Airspeed range [m/s]
URange = vcat(30:1:120)

# Tip correction function types
tipLossTypes = ["Exponential", "VLM-def"]

# Aerodynamic solver
aeroSolver = Indicial()

# Airfoil section
mod = modified_database(
    "FlatPlate";
    attached = (
        cm₀ = 0e-3,
        x_ac = 0.23,
    ),
)
airfoil = create_Airfoil(database=mod)
# airfoil = deepcopy(FlatPlate)

# Mesh: :nominal or :refined
mesh = :nominal

# Flag to use old data matrices from AePW3
useAePW3Matrices = false

# Fixed geometrical and discretization properties
nElem,L,chord,normSparPos,_ = geometrical_properties_Pazy(; mesh=mesh, useAePW3Matrices=useAePW3Matrices)

# Number of modes to track
nModes = 5

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :krylovkit
nev = mesh == :nominal ? 45 : 65

# Initialize outputs
eigensols = Array{EigenSolution}(undef,length(tipLossTypes),length(θRange),length(URange))
tip_OOP = Array{Float64}(undef,length(tipLossTypes),length(θRange),length(URange))
freqsOsc = Array{Matrix{Float64}}(undef,length(tipLossTypes),length(θRange))
dampsOsc = Array{Matrix{Float64}}(undef,length(tipLossTypes),length(θRange))
flutterOnsetSpeeds = Array{Vector{Float64}}(undef,length(tipLossTypes),length(θRange))
flutterOnsetFreqs = Array{Vector{Float64}}(undef,length(tipLossTypes),length(θRange))
flutterOnsetModes = Array{Vector{Int}}(undef,length(tipLossTypes),length(θRange))
flutterOnsetDisp = Array{Vector{Float64}}(undef,length(tipLossTypes),length(θRange))
flutterOffsetSpeeds = Array{Vector{Float64}}(undef,length(tipLossTypes),length(θRange))
flutterOffsetFreqs = Array{Vector{Float64}}(undef,length(tipLossTypes),length(θRange))
flutterOffsetModes = Array{Vector{Int}}(undef,length(tipLossTypes),length(θRange))
flutterOffsetDisp = Array{Vector{Float64}}(undef,length(tipLossTypes),length(θRange))

# Sweep tip loss types
for (k,tipLossType) in enumerate(tipLossTypes)
    # Sweep root angle
    for (i,θ) in enumerate(θRange)
        # Initialize guess solution
        x0 = zeros(0)
        # Sweep airspeed
        for (j,U) in enumerate(URange)
            # Display progress
            println("Solving for tip loss type = $tipLossType, θ = $(round(θ, digits=2)) deg, U = $U m/s")
            # Model
            PazyWingFlutterPitchRangeTipLossRange,_ = create_Pazy(aeroSolver=aeroSolver,airfoil=airfoil,hasTipCorrection=true,tipLossType=tipLossType,upright=true,pitchAngle=deg2rad(θ),airspeed=U,mesh=mesh,useAePW3Matrices=useAePW3Matrices)
            # Create and solve problem
            problem = create_EigenProblem(model=PazyWingFlutterPitchRangeTipLossRange,eigenmethod=eigenmethod,nev=nev,x0=x0)
            solve!(problem)
            # Tip OOP displacement at midchord
            states = problem.loadstepHistory.states
            tipTwist = twist_angle(states[end][nElem].nodal.p2_b, outputUnit=:deg)
            tip_OOP[k,i,j] = -(states[end][nElem].nodal.u2[1] - chord*(1/2-normSparPos)*sind(tipTwist))/L*100
            # Eigensolution
            eigensols[k,i,j] = problem.eigensolution
            # Set initial guess solution for next iteration
            x0 = problem.systemArrays.x
        end
        # Apply mode tracking
        freqs,damps,_ = mode_tracking(eigensols[k,i,:])
        rads2Hz!(freqs)
        # Get frequencies and dampings of oscillatory modes
        freqsOsc[k,i],dampsOsc[k,i],_ = filter_tracked_modes(freqs,damps; fmin=1e-1, maxTracks=nModes)
        # Flutter onset and offset speeds, frequencies and tip OOP displacements
        flutterOnsetSpeeds[k,i],flutterOnsetFreqs[k,i],flutterOnsetModes[k,i],flutterOnsetDisp[k,i],flutterOffsetSpeeds[k,i],flutterOffsetFreqs[k,i],flutterOffsetModes[k,i],flutterOffsetDisp[k,i] = find_flutter(URange,freqsOsc[k,i],dampsOsc[k,i],tip_OOP[k,i,:])
    end
end

# Load reference data
rootPitchVelUp = [3, 5, 7]
rootPitchVelDown = [3, 5]
flutterOnsetVelUp = [49, 43, 38]
flutterOffsetVelUp = [58, 51, 46]
flutterOnsetVelDown = [55, 48]
flutterOffsetVelDown = [40, 36]
flutterOnsetDispUp = [23, 24.5, 26]
flutterOnsetDispDown = [28.5, 31.5]

## Plots
# Set paths
relPath = "/dev/Pazy/S0/InertiaRod/outputs/PazyWingFlutterPitchRangeTipLossRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
ts = 10
fs = 16
lfs = 9
lw = 2
ms = 10
ms2 = 4
msw = 0
gr()

# Colors and line styles per tip loss type
tipLossColors = cgrad([:red,:blue], length(tipLossTypes), categorical=true)
tipLossStyles = [:dash,:dot]

# Extract flutter onset and offset data for the hump mode (NaN where the hump mode is not found)
humpMode = 3
UfOn = fill(NaN, length(θRange), length(tipLossTypes))
UfOff = fill(NaN, length(θRange), length(tipLossTypes))
ufOn = fill(NaN, length(θRange), length(tipLossTypes))
ufOff = fill(NaN, length(θRange), length(tipLossTypes))
for k in eachindex(tipLossTypes)
    for i in eachindex(θRange)
        for m in eachindex(flutterOnsetModes[k,i])
            if flutterOnsetModes[k,i][m] == humpMode
                UfOn[i,k] = flutterOnsetSpeeds[k,i][m]
                ufOn[i,k] = flutterOnsetDisp[k,i][m]
                break
            end
        end
        for m in eachindex(flutterOffsetModes[k,i])
            if flutterOffsetModes[k,i][m] == humpMode
                UfOff[i,k] = flutterOffsetSpeeds[k,i][m]
                ufOff[i,k] = flutterOffsetDisp[k,i][m]
                break
            end
        end
    end
end

# Angles at which the hump mode flutter region is complete (onset and offset found), per tip loss type
validHump = [.!isnan.(UfOn[:,k]) .& .!isnan.(UfOff[:,k]) .& .!isnan.(ufOn[:,k]) .& .!isnan.(ufOff[:,k]) for k in eachindex(tipLossTypes)]

# Flutter onset and offset speeds vs root pitch angle
plt1 = plot(xlabel="Airspeed [m/s]", ylabel="Root pitch angle [deg]", xlims=[30,90], ylims=[0,7.25], xticks=collect(30:10:90), yticks=collect(0:1:7), legend=:topright, tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
for k in eachindex(tipLossTypes)
    valid = validHump[k]
    if any(valid)
        plot!(Shape(vcat(UfOn[valid,k],reverse(UfOff[valid,k])),vcat(θRange[valid],reverse(θRange[valid]))), fillcolor = plot_color(tipLossColors[k], 0.25), lw=lw, label="Flutter region - $(tipLossTypes[k])")
    end
end
scatter!(flutterOnsetVelUp, rootPitchVelUp, shape=:rtriangle, mc=:red, ms=ms, msw=msw, label="Test onset up")
scatter!(flutterOffsetVelUp, rootPitchVelUp, shape=:rtriangle, mc=:green, ms=ms, msw=msw, label="Test offset up")
scatter!(flutterOnsetVelDown, rootPitchVelDown, shape=:ltriangle, mc=:red, ms=ms, msw=msw, label="Test onset down")
scatter!(flutterOffsetVelDown, rootPitchVelDown, shape=:ltriangle, mc=:green, ms=ms, msw=msw, label="Test offset down")
display(plt1)
savefig(string(absPath,"/PazyWingFlutterPitchRangeTipLossRange_flutterBoundaryPitch.pdf"))

# Flutter onset and offset speeds vs tip OOP displacement for varying root pitch angle
θ2plot = [0.5,1,2,3,5,7]
indθ2plot = findall(x -> any(isapprox(x, θ; atol=1e-10) for θ in θ2plot), θRange)
θcolors = cgrad([:blue, :red], length(θ2plot), categorical=true)
xθ = [ 5,  9, 12, 14, 16, 17.5]
yθ = [59, 55, 48, 44, 38, 34]
plt2 = plot(xlabel="Tip OOP displacement [% semispan]", ylabel="Airspeed [m/s]", xlims=[0,32], ylims=[30,90], xticks=collect(0:5:30), yticks=collect(30:10:90), legend=:bottomleft, tickfont=font(ts), guidefont=font(fs), legendfontsize=10)
for k in eachindex(tipLossTypes)
    valid = validHump[k]
    if any(valid)
        plot!(Shape(vcat(ufOn[valid,k],reverse(ufOff[valid,k])),vcat(UfOn[valid,k],reverse(UfOff[valid,k]))), fillcolor = plot_color(tipLossColors[k], 0.25), lw=lw, label="Flutter region - $(tipLossTypes[k])")
    end
end
for (n,ind) in enumerate(indθ2plot)
    θ = θ2plot[n]
    for k in eachindex(tipLossTypes)
        plot!(tip_OOP[k,ind,:], URange, c=θcolors[n], ls=tipLossStyles[k], lw=lw, label=false)
    end
    annotate!([xθ[n]],[yθ[n]], text("\$$(round(θ2plot[n],digits=1)) ^\\circ\$", 10, :bottom, θcolors[n]))
end
scatter!(flutterOnsetDispUp, flutterOnsetVelUp, shape=:rtriangle, mc=:red, ms=ms, msw=msw, label=false)
scatter!(flutterOnsetDispDown, flutterOnsetVelDown, shape=:ltriangle, mc=:red, ms=ms, msw=msw, label=false)
display(plt2)
savefig(string(absPath,"/PazyWingFlutterPitchRangeTipLossRange_flutterBoundaryDisp.pdf"))

println("Finished PazyWingFlutterPitchRangeTipLossRange.jl")