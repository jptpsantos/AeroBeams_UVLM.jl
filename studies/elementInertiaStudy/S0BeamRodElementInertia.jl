using AeroBeams

# Study: nodal inertia data (half of each node's lumped mass to each adjacent element, as in `inertia_matrices_from_Nastran_data`) vs inertia data condensed directly to the midpoints of the elements (BeamROM `condense_to_elements`, see generate_element_inertia.jl), for the straight Pazy wing with the tip rod modeled as a beam. The modes are compared with the Nastran model with the flexible tip rod, as in S0BeamRodModal.jl

meshes = [:nominal, :refined]
nElementsWingtipRod = 12
nModes = 15
eigenmethod = :krylovkit
nev = 80
rodModeRatio = 0.07
studyDir = @__DIR__

_,_,chord,referenceAxis,_ = geometrical_properties_Pazy()
referenceDir = joinpath(pkgdir(AeroBeams),"test","referenceData","Pazy")
modalFreqsNastran = vec(readdlm(joinpath(referenceDir,"S0_modal_frequencies_flextip.csv"),',')[:,3])
wingShapesNastran = readdlm(joinpath(referenceDir,"S0_modal_shapes_flextip_wing.csv"),',',Float64,skipstart=1)
rodShapesNastran = readdlm(joinpath(referenceDir,"S0_modal_shapes_flextip_rod.csv"),',',Float64,skipstart=1)

function interpolate_data(x,xs,ys)
    x <= xs[1] && return ys[1]
    x >= xs[end] && return ys[end]
    k = searchsortedlast(xs,x)
    return ys[k] + (x-xs[k])/(xs[k+1]-xs[k]) * (ys[k+1]-ys[k])
end

function sampling_elements(model,nElemWing)
    midpoints = [(model.assemblyData.r_n[n[1]] + model.assemblyData.r_n[n[2]])/2 for n in model.assemblyData.elementNodes]
    return midpoints, 1:nElemWing, (nElemWing+2):model.assemblyData.nElementsTotal
end

function beam_mode_vector(modeShape,wingElements,rodElements)
    toNastranBasis(w) = [w[3], w[2], -w[1]]
    v = Float64[]
    for e in wingElements
        u = modeShape.state[e].midpoint.u
        p = modeShape.state[e].midpoint.p
        for η in (referenceAxis*chord, (referenceAxis-1)*chord)
            append!(v, toNastranBasis(u + cross(p,[0,η,0])))
        end
    end
    for e in rodElements
        append!(v, toNastranBasis(modeShape.state[e].midpoint.u))
    end
    return v
end

function nastran_mode_vector(mode,midpoints,wingElements,rodElements)
    wing = wingShapesNastran[wingShapesNastran[:,1] .== mode, :]
    rod = rodShapesNastran[rodShapesNastran[:,1] .== mode, :]
    toBasisA(u) = [u[2], -u[1], u[3]]
    v = Float64[]
    for e in wingElements, columns in (3:5, 6:8)
        append!(v, toBasisA([interpolate_data(midpoints[e][3],wing[:,2],wing[:,c]) for c in columns]))
    end
    for e in rodElements
        x = referenceAxis*chord - midpoints[e][2]
        append!(v, toBasisA([interpolate_data(x,rod[:,2],rod[:,c]) for c in 3:5]))
    end
    return v
end

mac(a,b) = dot(a,b)^2 / (dot(a,a)*dot(b,b))

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

# Solves the modal problem of a mesh and returns the frequencies, the paired Nastran frequencies and the MACs
function solve_case(mesh)
    model,nElemWing,_ = create_Pazy(aeroSolver=QuasiSteady(),mesh=mesh,beamWingtipRod=true,nElementsWingtipRod=nElementsWingtipRod,gravityVector=zeros(3))
    problem = create_EigenProblem(model=model,eigenmethod=eigenmethod,nev=nev)
    solve!(problem)
    modeShapes = sort(filter(m -> m.frequency > 0, problem.eigensolution.modeShapes), by=m -> m.frequency)[1:nModes]
    freqs = [rads2Hz(m.frequency) for m in modeShapes]
    midpoints,wingElements,rodElements = sampling_elements(model,nElemWing)
    beamVectors = [beam_mode_vector(m,wingElements,rodElements) for m in modeShapes]
    nastranVectors = [nastran_mode_vector(mode,midpoints,wingElements,rodElements) for mode in 1:nModes]
    MAC = [mac(a,b) for a in beamVectors, b in nastranVectors]
    paired = pair_modes(MAC)
    return freqs, modalFreqsNastran[paired], [MAC[k,paired[k]] for k in 1:nModes], paired
end

# Nodal inertia data (the current implementation)
nodal = Dict(mesh => solve_case(mesh) for mesh in meshes)

# Inertia matrices as a blend of the nodal data (the nodes' masses split between the adjacent elements, as in the current implementation: the lumped mass matrix of the beam) and the data at the midpoints of the elements (offsets from the midpoints, in the Nastran frame: the consistent mass matrix): α I_nodal + (1-α) I_element, with α set through `lumpedWeight`. The sectional inertia matrices are additive in the mass, so the blend has the same total mass
lumpedWeight = Ref(1.0)
AeroBeams.eval(quote
    function inertia_matrices_Pazy(Λ; mesh::Symbol=:nominal, kwargs...)
        nElem,L,_,_,nodalPositions = geometrical_properties_Pazy(mesh=mesh,Λ=Λ)
        # Nodal data
        file = Pazy_data_file(:inertia, Λ, mesh; withWingtipRod=false)
        d, h = readdlm(file, ','; header=true)
        c(name) = Float64.(d[:, findfirst(==(name), vec(h))])
        I = inertia_matrices_from_Nastran_data(Λ=Λ,nElem=nElem,L=L,nodalPositions=nodalPositions,m=c("m"),m_x=c("dx"),m_y=c("dy"),m_z=c("dz"),Ixx=c("Ixx"),Iyy=c("Iyy"),Izz=c("Izz"),Ixy=c("Ixy"),Ixz=c("Ixz"),Iyz=c("Iyz"))
        # Element data
        data, header = readdlm(joinpath($studyDir,"S$(round(Int,rad2deg(Λ)))_norod_inertia_matrices_$(mesh)mesh_elements.csv"), ','; header=true)
        column(name) = Float64.(data[:, findfirst(==(name), vec(header))])
        @assert size(data,1) == nElem
        R = [ sin(Λ) cos(Λ) 0; -cos(Λ) sin(Λ) 0; 0 0 1]
        m,dx,dy,dz = column("m"),column("dx"),column("dy"),column("dz")
        Ixx,Iyy,Izz,Ixy,Ixz,Iyz = column("Ixx"),column("Iyy"),column("Izz"),column("Ixy"),column("Ixz"),column("Iyz")
        α = $lumpedWeight[]
        return [α*I[k] + (1-α)*sectional_inertia_elem(m[k], R*[dx[k],dy[k],dz[k]], R*[Ixx[k] Ixy[k] Ixz[k]; Ixy[k] Iyy[k] Iyz[k]; Ixz[k] Iyz[k] Izz[k]]*R', L*(nodalPositions[k+1]-nodalPositions[k])) for k in 1:nElem]
    end
end)
sets = [("element (consistent)", 0.0), ("blend α=1/4", 0.25), ("blend α=1/2", 0.5), ("blend α=3/4", 0.75), ("nodal (lumped)", 1.0)]
results = Dict{Tuple{String,Symbol},Any}()
for (name,α) in sets, mesh in meshes
    lumpedWeight[] = α
    results[(name,mesh)] = solve_case(mesh)
end

# Report: mean and maximum |Δ| to Nastran over the modes, and mean MAC
avg(x) = sum(x)/length(x)
r(x,d) = round(x,digits=d)
rows = []
for mesh in meshes, (name,_) in sets
    f,fN,M,_ = results[(name,mesh)]
    push!(rows, [mesh, name, r(100*avg(abs.(f./fN.-1)),3), r(100*maximum(abs.(f./fN.-1)),3), r(avg(M),4)])
end
print(pretty_table(String, reduce(vcat,permutedims.(rows)); column_labels=["Mesh","Inertia data","mean |Δ| [%]","max |Δ| [%]","mean MAC"], title="Modes 1-$nModes of the S0 beam-rod wing vs Nastran (flextip)"))
for mesh in meshes
    deltas = [r.(100*(results[(name,mesh)][1]./results[(name,mesh)][2].-1),2) for (name,_) in sets]
    print(pretty_table(String, hcat(1:nModes, deltas...); column_labels=["Mode", [string(name," Δ [%]") for (name,_) in sets]...], title="Mesh :$mesh"))
end
println("Finished S0BeamRodElementInertia.jl")
