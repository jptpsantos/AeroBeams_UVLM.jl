# Vibration modes of the Pazy wing (any sweep angle and tip mass configuration) in the undeformed state, without gravity, from the eigenproblem of AeroBeams without airflow, and the data needed to compute their aerodynamic forces and the flutter in their basis: the generalized masses (from the strain energy of the mode, M = 2U/ω²), the bending displacement of the reference axis along the lift direction (h) and its twist, nose up positive (θ), at the nodes of the wing
using AeroBeams, LinearAlgebra

function pazy_modal_data(; sweepAngle=0.0, tipMassConfig=:baseline, mesh=:nominal, beamWingtipRod=false, tipRodInertiaCorrection=true, nModes=12, nev=40)

    # Options of the model: the rod is modeled in one of the two ways (lumped in the tip node, with or without the ad hoc correction of its inertia, or as a beam)
    straight = isapprox(sweepAngle, 0; atol=1e-8)
    kwargs = straight ?
             (;) :
             (; tipMassConfig=tipMassConfig, beamWingtipRod=beamWingtipRod, tipRodInertiaCorrection=!beamWingtipRod && tipRodInertiaCorrection)

    # Model without airflow (the quasi-steady solver has no aerodynamic states)
    model,_ = create_Pazy(aeroSolver=QuasiSteady(), mesh=mesh, gravityVector=zeros(3), sweepAngle=sweepAngle; kwargs...)
    problem = create_EigenProblem(model=model, eigenmethod=:krylovkit, nev=nev, x0=zeros(0))
    solve!(problem)
    assemblyData = problem.model.assemblyData
    elements = assemblyData.elements
    eigensolution = problem.eigensolution

    # Wing elements (the first ones, as the wing is the first beam), and the arclength of their nodes along the reference axis
    nElem,L,chord,referenceAxis,_ = geometrical_properties_Pazy(Λ=sweepAngle, mesh=mesh)
    s = [[elements[e].geometry.x1_n1 for e in 1:nElem]; elements[nElem].geometry.x1_n2]

    # Oscillatory modes in ascending order of frequency
    ω = imag.(eigensolution.λ)
    selected = findall(>(1e-3), ω)
    selected = selected[sortperm(ω[selected])][1:nModes]

    # Indices of the displacement-like states, to get real snapshots of the mode shapes
    idx = reduce(vcat, vcat(e.assemblyIndices.DOF_u, e.assemblyIndices.DOF_p, e.assemblyIndices.DOF_F, e.assemblyIndices.DOF_M) for e in elements)
    modes = map(selected) do m
        vector = AeroBeams.real_eigenvector(@view(eigensolution.V[:,m]), idx)
        states = [AeroBeams.ModalState() for _ in elements]
        for (e,element) in enumerate(elements)
            AeroBeams.element_modal_states!(element, states[e], vector, assemblyData.forceScaling)
        end
        # Strain energy of the mode (all elements, including those of the rod), which gives the generalized stiffness and mass
        energy = sum(0.5*element.geometry.Δℓ*(dot(state.midpoint.γ, state.midpoint.F) + dot(state.midpoint.κ, state.midpoint.M)) for (element,state) in zip(elements,states))
        K = 2energy
        h = [[states[e].nodal.u1_b[3] for e in 1:nElem]; states[nElem].nodal.u2_b[3]]
        θ = [[states[e].nodal.p1_b[1] for e in 1:nElem]; states[nElem].nodal.p2_b[1]]
        (; frequency=ω[m]/2π, K=K, M=K/ω[m]^2, h=h, θ=θ)
    end

    return (; modes, s, chord, referenceAxis, length=L, sweepAngle, tipMassConfig, beamWingtipRod)

end
