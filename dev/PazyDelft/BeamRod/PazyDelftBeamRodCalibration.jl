using AeroBeams

# Calibration of the wingtip rod of the TU Delft Pazy wing modeled as a beam (create_PazyDelft with beamWingtipRod=true) with the natural frequencies and mode shapes of the Nastran model with the flexible tip rod (SOL 103, clamped at the root, in vacuo). The calibrated parameters are the limits of the effectively rigid parts of the rod, glued to the wing along the chord, for the bending in the plane of the wing and out of it (`tied` in `wingtip_rod_PazyDelft`): the effective clamps of the cantilevers of the rod beyond the chord
# The modes of the beam model are paired with those of Nastran by their shapes (modal assurance criterion, as in PazyDelftBeamRodModal.jl), and the parameters minimize the relative errors of the frequencies of the modes, with a weak prior at the starting values (Levenberg-Marquardt with finite difference Jacobian). The results are the values in `wingtip_rod_PazyDelft`

# Meshes
meshes = [:nominal, :refined]

# Number of elements of the wingtip rod beam, number of modes of the beam model to pair, and eigensolver
nElementsWingtipRod = 12
nModes = 20
nev = 80

# Nastran modes in the calibration (all of them) and their weights (the first mode of the rod in the plane of the wing, mode 8, in which the rod bends in the plane against the wing, is the one that matters most to the flutter)
calibratedModes = 1:15
weights = Dict(8 => 6.0)

# Starting values (the limits of the rigid parts are x-coordinates of the Nastran frame, in mm, for both planes, from x = 0 to the chord of 98.9 mm), their bounds, and the standard deviations of the weak prior and of the errors of the frequencies [mm, %]
chord = AeroBeams.chord_PazyDelft
p0 = [0, chord*1e3, 0, chord*1e3]
xSpar = chord*AeroBeams.shear_center_PazyDelft # x-coordinate of the tip node of the spar [m]
xLERod, xTERod = -0.0841, 0.1841 # ends of the rod [m]
lowerBounds = [1000*xLERod+8, 1000*xSpar+8, 1000*xLERod+8, 1000*xSpar+8]
upperBounds = [1000*xSpar-8, 1000*xTERod-8, 1000*xSpar-8, 1000*xTERod-8]
priorStd = [15.0, 15.0, 15.0, 15.0]
frequencyStd = 0.5

# Load reference data: natural frequencies of the Nastran model and their mode shapes: translations (Nastran basic frame: x aft along the chord, y along the span, unswept, z vertical) of the leading and trailing edge points of the sections of the spar, along the span from the root, and of the grids of the tip rod, along the chord from the leading edge of the section at the rod (see extractModeShapes.py of BeamROM)
referenceDir = joinpath(pkgdir(AeroBeams),"test","referenceData","PazyDelft")
modalFreqsNastran = vec(readdlm(joinpath(referenceDir,"pazy_delft_modal_frequencies_flextip.csv"),',')[:,3])
wingShapesNastran = readdlm(joinpath(referenceDir,"pazy_delft_modal_shapes_flextip_wing.csv"),',',Float64,skipstart=1)
rodShapesNastran = readdlm(joinpath(referenceDir,"pazy_delft_modal_shapes_flextip_rod.csv"),',',Float64,skipstart=1)
nNastranModes = length(modalFreqsNastran)

# Normalized position of the reference axis (shear center)
referenceAxis = AeroBeams.shear_center_PazyDelft

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

# Data of the rod for the parameters `p`, from the default data of `wingtip_rod_PazyDelft`
function rod_data(p)
    rod = AeroBeams.wingtip_rod_PazyDelft()
    tied = (inPlane = (p[1]/1000 - xSpar, p[2]/1000 - xSpar), outOfPlane = (p[3]/1000 - xSpar, p[4]/1000 - xSpar))
    return merge(rod, (tied = tied,))
end

# Modes of the beam model with the rod `rod` and the mesh `mesh`, paired with those of Nastran: frequency [Hz] and MAC of the beam mode paired with each Nastran mode (NaN and zero when it is not paired)
function paired_modes(mesh,rod)
    model,nElemWing,_ = create_PazyDelft(aeroSolver=QuasiSteady(),mesh=mesh,beamWingtipRod=true,nElementsWingtipRod=nElementsWingtipRod,wingtipRodData=rod,gravityVector=zeros(3))
    problem = create_EigenProblem(model=model,eigenmethod=:krylovkit,nev=nev)
    solve!(problem)
    modeShapes = sort(filter(m -> m.frequency > 0, problem.eigensolution.modeShapes), by=m -> m.frequency)[1:nModes]
    # Sampling points: midpoints of the elements of the wing and of the rod (relative to the clamp, in basis A), beyond the rigid link
    midpoints = [(model.assemblyData.r_n[n[1]] + model.assemblyData.r_n[n[2]])/2 for n in model.assemblyData.elementNodes]
    wingElements = 1:nElemWing
    rodElements = (nElemWing+2):model.assemblyData.nElementsTotal
    # Translations (in basis A) of the leading and trailing edge points of the wing and of the rod at the sampling points, for the beam and for the Nastran modes
    beam_vector(m) = vcat([vcat([m.state[e].midpoint.u + cross(m.state[e].midpoint.p,[0,η,0]) for η in (referenceAxis*chord, (referenceAxis-1)*chord)]...) for e in wingElements]..., [m.state[e].midpoint.u for e in rodElements]...)
    function nastran_vector(mode)
        wing = wingShapesNastran[wingShapesNastran[:,1] .== mode, :]
        rodShape = rodShapesNastran[rodShapesNastran[:,1] .== mode, :]
        toBasisA(u) = [u[2], -u[1], u[3]]
        v = Float64[]
        for e in wingElements, columns in (3:5, 6:8)
            append!(v, toBasisA([interpolate_data(midpoints[e][1],wing[:,2],wing[:,c]) for c in columns]))
        end
        for e in rodElements
            x = referenceAxis*chord - midpoints[e][2] # chordwise position of the element, from the leading edge of the section at its span position
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
    for mesh in meshes
        freqs,_ = paired_modes(mesh,rod_data(p))
        for mode in calibratedModes
            error = isnan(freqs[mode]) ? 20.0 : 100*log(freqs[mode]/modalFreqsNastran[mode])
            push!(r, get(weights,mode,1.0) * error / frequencyStd)
        end
    end
    return vcat(r, (p .- p0) ./ priorStd)
end

# Levenberg-Marquardt, with the finite difference Jacobian and the parameters within their bounds
function levenberg_marquardt(p; maxIterations=30, steps=[0.2,0.2,0.2,0.2])
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
results = [(mesh, paired_modes(mesh,rod_data(pOpt))...) for mesh in meshes]
data = vcat([hcat(fill(string(mesh),nNastranModes), 1:nNastranModes, freqs, modalFreqsNastran, 100 .* (freqs ./ modalFreqsNastran .- 1), macs) for (mesh,freqs,macs) in results]...)
table = pretty_table(
    String,
    data;
    column_labels = ["Mesh", "Nastran mode", "AeroBeams [Hz]", "Nastran [Hz]", "Δ Nastran [%]", "MAC"],
    title         = "Modal Frequency Comparison with the Calibrated Wingtip Rod",
    subtitle      = "Modes used in the calibration: " * join(string.(calibratedModes), ", "),
    formatters    = [
        (v, i, j) -> (j in (3, 4) && v isa Number) ? round(v, digits=2) : v,
        (v, i, j) -> (j == 5 && v isa Number) ? string(round(v, digits=2), " %") : v,
        (v, i, j) -> (j == 6 && v isa Number) ? round(v, digits=3) : v,
    ],
)
print(table)

# Calibrated data, in the format of `wingtip_rod_PazyDelft` (the limits are x-coordinates of the Nastran frame [m])
println("tied = (inPlane = (", round(pOpt[1]/1000,digits=4), " - xSpar, ", round(pOpt[2]/1000,digits=4), " - xSpar), outOfPlane = (", round(pOpt[3]/1000,digits=4), " - xSpar, ", round(pOpt[4]/1000,digits=4), " - xSpar))")

println("Finished PazyDelftBeamRodCalibration.jl")
