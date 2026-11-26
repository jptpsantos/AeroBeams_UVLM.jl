using AeroBeams

# Calibration of the wingtip rod of the Pazy wing with 10 degrees of sweep modeled as a beam (create_Pazy with beamWingtipRod=true) with the natural frequencies and mode shapes of the Nastran models with the flexible tip rod (SOL 103, clamped at the root, in vacuo) of the LE and TE configurations. The calibrated parameters are:
# - the limits of the effectively rigid parts of the rod, glued to the wing along the chord, for the bending in the plane of the wing and out of it (`tied` in `wingtip_rod_Pazy`): the effective clamps of the cantilevers of the rod beyond the chord
# - the factors of the springs around the hinge of the breaking beam (`hinge.springFactors`), in bending in the plane of the wing and out of it, for each configuration, which multiply the stiffness EI/L of the breaking beam. The one in torsion is not calibrated: the FPS weights have no inertia about the axis of the rod, so the modes do not depend on it
# The modes of the beam model are paired with those of Nastran by their shapes (modal assurance criterion, as in S10Modal.jl), and the parameters minimize the relative errors of the frequencies of the rod modes and of the wing modes affected by the rod, with a weak prior at the starting values (Levenberg-Marquardt with finite difference Jacobian). The results are the values in `wingtip_rod_Pazy`

# Tip mass configurations, meshes and sweep angle [rad]
configs = [:LE, :TE]
meshes = [:nominal, :refined]
Λ = deg2rad(10)

# Number of elements of the wingtip rod beam (without the extension to the FPS weights), number of modes of the beam model to pair, and eigensolver
nElementsWingtipRod = 16
nModes = 20
nev = 80

# Nastran modes in the calibration for each configuration (the first torsion mode of the wing, the modes of the rod in the plane of the wing, in which the FPS weights move against it, and the modes of the wing and of the rod in out-of-plane bending and torsion that are affected by them) and their weights (the first mode of the rod in the plane is the one that matters most to the flutter)
calibratedModes = Dict(:LE => [3,4,5,6,9,11,12,13], :TE => [3,4,5,7,9,12])
weights = Dict((:LE,4) => 6.0, (:TE,4) => 6.0)

# Starting values (the limits of the rigid parts are x-coordinates of the Nastran frame, in mm, x = 106.4 to 173.37 for both planes; the springs are the nominal ones, with factors of 1, whose logarithms are used), their bounds, and the standard deviations of the weak prior and of the errors of the frequencies [mm, log of the factors, %]
parameterNames = ["inPlane start [mm]", "inPlane end [mm]", "outOfPlane start [mm]", "outOfPlane end [mm]", "ln kIP LE", "ln kOOP LE", "ln kIP TE", "ln kOOP TE"]
p0 = [106.4, 173.37, 106.4, 173.37, 0.0, 0.0, 0.0, 0.0]
xSpar = 0.138871796846351 # x-coordinate of the tip node of the spar [m]
xLERod, xTERod = 0.05289, 0.23278 # ends of the rod [m]
lowerBounds = [1000*xLERod+8, 1000*xSpar+8, 1000*xLERod+8, 1000*xSpar+8, log(0.1), log(0.1), log(0.1), log(0.1)]
upperBounds = [1000*xSpar-8, 1000*xTERod-8, 1000*xSpar-8, 1000*xTERod-8, log(10), log(10), log(10), log(10)]
priorStd = [15.0, 15.0, 15.0, 15.0, 0.7, 0.7, 0.7, 0.7]
frequencyStd = 0.5

# Load reference data: natural frequencies of the Nastran models and their mode shapes: translations (Nastran basic frame: x aft along the chord, y along the span, unswept, z vertical) of the leading and trailing edge points of the sections of the spar, along the span from the root, and of the grids of the tip rod, along the chord from the leading edge of the section at the rod (see extractModeShapes.py of BeamROM)
referenceDir = joinpath(pkgdir(AeroBeams),"test","referenceData","Pazy")
modalFreqsNastran = Dict(config => vec(readdlm(joinpath(referenceDir,"S10_$(config)_modal_frequencies_flextip.csv"),',')[:,3]) for config in configs)
wingShapesNastran = Dict(config => readdlm(joinpath(referenceDir,"S10_$(config)_modal_shapes_flextip_wing.csv"),',',Float64,skipstart=1) for config in configs)
rodShapesNastran = Dict(config => readdlm(joinpath(referenceDir,"S10_$(config)_modal_shapes_flextip_rod.csv"),',',Float64,skipstart=1) for config in configs)
nNastranModes = length(modalFreqsNastran[configs[1]])

# Chord of the wing and normalized position of the reference axis
_,_,chord,referenceAxis,_ = geometrical_properties_Pazy(Λ=Λ)

# Linear interpolation of the data `ys` at the increasing points `xs`, evaluated at `x` (constant beyond the ends)
function interpolate_data(x,xs,ys)
    x <= xs[1] && return ys[1]
    x >= xs[end] && return ys[end]
    k = searchsortedlast(xs,x)
    return ys[k] + (x-xs[k])/(xs[k+1]-xs[k]) * (ys[k+1]-ys[k])
end

# Modal assurance criterion between two mode shape vectors
mac(a,b) = dot(a,b)^2 / (dot(a,a)*dot(b,b))

# Pair the modes of the beam model with those of Nastran by their shapes: the pair of largest MAC among the modes not paired yet is taken, until all the modes of Nastran are paired (the modes of the beam that are not paired have a zero)
function pair_modes(MAC)
    pairedNastran = zeros(Int,size(MAC,1))
    availableBeam = trues(size(MAC,1))
    availableNastran = trues(size(MAC,2))
    for _ in 1:min(size(MAC)...)
        _,idx = findmax(MAC .* (availableBeam * availableNastran'))
        pairedNastran[idx[1]] = idx[2]
        availableBeam[idx[1]] = false
        availableNastran[idx[2]] = false
    end
    return pairedNastran
end

# Data of the rod for the parameters `p`, from the default data of `wingtip_rod_Pazy`
function rod_data(config,p)
    rod = AeroBeams.wingtip_rod_Pazy(Λ,config)
    tied = (inPlane = (p[1]/1000 - xSpar, p[2]/1000 - xSpar), outOfPlane = (p[3]/1000 - xSpar, p[4]/1000 - xSpar))
    kIP,kOOP = config == :LE ? exp.(p[5:6]) : exp.(p[7:8])
    return merge(rod, (tied = tied, hinge = merge(rod.hinge, (springFactors = (torsion = 1.0, inPlane = kIP, outOfPlane = kOOP),))))
end

# Modes of the beam model of the configuration and mesh with the rod `rod`, paired with those of Nastran: frequency [Hz] and MAC of the beam mode paired with each Nastran mode (NaN and zero when it is not paired)
function paired_modes(config,mesh,rod)
    model,nElemWing,_ = create_Pazy(aeroSolver=QuasiSteady(),sweepAngle=Λ,tipMassConfig=config,mesh=mesh,beamWingtipRod=true,nElementsWingtipRod=nElementsWingtipRod,wingtipRodData=rod,gravityVector=zeros(3))
    problem = create_EigenProblem(model=model,eigenmethod=:krylovkit,nev=nev)
    solve!(problem)
    modeShapes = sort(filter(m -> m.frequency > 0, problem.eigensolution.modeShapes), by=m -> m.frequency)[1:nModes]
    # Sampling points: midpoints of the elements of the wing and of the rod (relative to the clamp, in basis A), without the two elements at the hinge side, with the adapter and the FPS weights, beyond the end of the rod of the Nastran data
    midpoints = [(model.assemblyData.r_n[n[1]] + model.assemblyData.r_n[n[2]])/2 for n in model.assemblyData.elementNodes]
    wingElements = 1:nElemWing
    rodElements = (nElemWing + 2 + (config == :LE ? 2 : 0)):(model.assemblyData.nElementsTotal - (config == :TE ? 2 : 0))
    # Translations (in basis A) of the leading and trailing edge points of the wing and of the rod at the sampling points, for the beam and for the Nastran modes
    beam_vector(m) = vcat([vcat([m.state[e].midpoint.u + cross(m.state[e].midpoint.p,[0,η,0]) for η in (referenceAxis*chord, (referenceAxis-1)*chord)]...) for e in wingElements]..., [m.state[e].midpoint.u for e in rodElements]...)
    function nastran_vector(mode)
        wing = wingShapesNastran[config][wingShapesNastran[config][:,1] .== mode, :]
        rodShape = rodShapesNastran[config][rodShapesNastran[config][:,1] .== mode, :]
        toBasisA(u) = [u[2], -u[1], u[3]]
        v = Float64[]
        for e in wingElements, columns in (3:5, 6:8)
            append!(v, toBasisA([interpolate_data(midpoints[e][1],wing[:,2],wing[:,c]) for c in columns]))
        end
        for e in rodElements
            x = referenceAxis*chord - midpoints[e][2] - tan(Λ)*midpoints[e][1]
            append!(v, toBasisA([interpolate_data(x,rodShape[:,2],rodShape[:,c]) for c in 3:5]))
        end
        return v
    end
    beamVectors = [beam_vector(m) for m in modeShapes]
    nastranVectors = [nastran_vector(mode) for mode in 1:nNastranModes]
    MAC = [mac(a,b) for a in beamVectors, b in nastranVectors]
    pairing = pair_modes(MAC)
    freqs = fill(NaN,nNastranModes)
    macs = zeros(nNastranModes)
    for (k,mode) in enumerate(pairing)
        mode == 0 && continue
        freqs[mode] = rads2Hz(modeShapes[k].frequency)
        macs[mode] = MAC[k,mode]
    end
    return freqs,macs
end

# Residuals of the calibration for the parameters `p`: weighted relative errors [%] of the frequencies of the calibrated modes, divided by their standard deviation, and the weak prior
function residuals(p)
    r = Float64[]
    for config in configs, mesh in meshes
        freqs,_ = paired_modes(config,mesh,rod_data(config,p))
        for mode in calibratedModes[config]
            error = isnan(freqs[mode]) ? 20.0 : 100*log(freqs[mode]/modalFreqsNastran[config][mode])
            push!(r, get(weights,(config,mode),1.0) * error / frequencyStd)
        end
    end
    return vcat(r, (p .- p0) ./ priorStd)
end

# Levenberg-Marquardt, with the finite difference Jacobian and the parameters within their bounds
function levenberg_marquardt(p; maxIterations=30, steps=[0.2,0.2,0.2,0.2,0.01,0.01,0.01,0.01])
    λ = 1.0
    r = residuals(p)
    cost = sum(abs2,r)
    println("Start: cost = ", round(cost,digits=2))
    for iteration in 1:maxIterations
        previousCost = cost
        J = hcat([(residuals(p + steps[j]*(1:length(p) .== j)) - r) / steps[j] for j in eachindex(p)]...)
        H = J'*J
        g = J'*r
        improved = false
        for _ in 1:8
            δ = -(H + λ*Diagonal([H[i,i] + 1e-8 for i in axes(H,1)])) \ g
            pNew = clamp.(p + δ, lowerBounds, upperBounds)
            rNew = residuals(pNew)
            if sum(abs2,rNew) < cost
                p, r, cost = pNew, rNew, sum(abs2,rNew)
                λ = max(λ/3, 1e-6)
                improved = true
                break
            end
            λ *= 4
        end
        println("Iteration ", iteration, ": cost = ", round(cost,digits=3), ", parameters = ", round.(p,digits=2))
        (improved && previousCost - cost > 1e-5*previousCost) || break
    end
    return p
end

# Calibrate
pOpt = levenberg_marquardt(copy(p0))

# Frequencies of the paired modes with the calibrated parameters, compared with those of Nastran
results = [(config, mesh, paired_modes(config,mesh,rod_data(config,pOpt))...) for config in configs for mesh in meshes]
data = vcat([hcat(fill(string(config," ",mesh),nNastranModes), 1:nNastranModes, freqs, modalFreqsNastran[config], 100 .* (freqs ./ modalFreqsNastran[config] .- 1), macs) for (config,mesh,freqs,macs) in results]...)
table = pretty_table(
    String,
    data;
    column_labels = ["Configuration and mesh", "Nastran mode", "AeroBeams [Hz]", "Nastran [Hz]", "Δ Nastran [%]", "MAC"],
    title         = "Modal Frequency Comparison with the Calibrated Wingtip Rod",
    subtitle      = "Modes used in the calibration: " * join(["$config $(calibratedModes[config])" for config in configs], ", "),
    formatters    = [
        (v, i, j) -> (j in (3, 4) && v isa Number) ? round(v, digits=2) : v,
        (v, i, j) -> (j == 5 && v isa Number) ? string(round(v, digits=2), " %") : v,
        (v, i, j) -> (j == 6 && v isa Number) ? round(v, digits=3) : v,
    ],
)
print(table)

# Calibrated data, in the format of `wingtip_rod_Pazy` (the limits are x-coordinates of the Nastran frame [m])
println("tied = (inPlane = (", round(pOpt[1]/1000,digits=4), " - xSpar, ", round(pOpt[2]/1000,digits=4), " - xSpar), outOfPlane = (", round(pOpt[3]/1000,digits=4), " - xSpar, ", round(pOpt[4]/1000,digits=4), " - xSpar))")
println("LE springFactors = (torsion = 1.0, inPlane = ", round(exp(pOpt[5]),digits=2), ", outOfPlane = ", round(exp(pOpt[6]),digits=2), ")")
println("TE springFactors = (torsion = 1.0, inPlane = ", round(exp(pOpt[7]),digits=2), ", outOfPlane = ", round(exp(pOpt[8]),digits=2), ")")

println("Finished S10RodCalibration.jl")