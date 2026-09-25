# Coupling regression tests. No simulation outputs are written.
# Run: julia --project=lib/WingPropellerUVLM lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/verify_pazy_coupling.jl
# Tests intentionally assert the desired physical properties, so defects give FAIL,
# not an apparently successful test run. This file does not change either solver.
using Test
include(joinpath(@__DIR__, "PazyWingUVLMCoupling.jl"))

function audit_geometry()
    ne, span, chord, spar = AeroBeams.geometrical_properties_Pazy()
    eta = AeroBeams.nodal_positions_Pazy()
    R0 = [0.0 0.0 -1.0; 0.0 1.0 0.0; 1.0 0.0 0.0]
    elements = [(r_n1=[0.0, 0.0, span*eta[e]],
        r_n2=[0.0, 0.0, span*eta[e+1]], R0_n1=R0, R0_n2=R0,
        nodalStates=(u_n1=zeros(3), u_n2=zeros(3), p_n1=zeros(3),
            p_n2=zeros(3), p_n2_b=zeros(3))) for e in 1:ne]
    weights = zeros(16, ne+1)
    for j in 1:16
        station = (j-1)/15
        left = min(searchsortedlast(eta, station), ne)
        fraction = (station-eta[left])/(eta[left+1]-eta[left])
        weights[j,left] = 1-fraction
        weights[j,left+1] = fraction
    end
    return (elements=elements,), chord, spar, weights
end

function audit_vertices(model, chord, spar, weights; nc=4)
    grid, positions, offsets = wing_geometry(model, chord, spar, nc, weights)
    _, _, panels = UVLM.grid_to_surface_panels(grid)
    return UVLM.imperial_nodal_positions(panels), positions, offsets
end

function virtual_work_check(weights)
    model, chord, spar, _ = audit_geometry()
    points, positions, offsets = audit_vertices(model, chord, spar, weights)
    forces = [zeros(3) for k in axes(points,1), j in axes(points,2)]
    forces[2,8] = A_TO_UVLM * [1.0, 0.0, 0.0] # 1 N normal force.
    mapped = beam_loads((forces=[forces], positions=[points]), offsets, weights)

    # Perturb only node 8's rotation about the chord, holding positions fixed.
    # Its chord points must not move, so this virtual motion does no aero work.
    # At p=0, the WM parameter increment equals the infinitesimal rotation.
    epsilon = 1e-6
    model.elements[7].nodalStates.p_n2[2] = epsilon
    plus, _ = audit_vertices(model, chord, spar, weights)
    model.elements[7].nodalStates.p_n2[2] = -epsilon
    minus, _ = audit_vertices(model, chord, spar, weights)
    aerodynamic_work = sum(dot(forces[i], (plus[i]-minus[i])/(2epsilon))
        for i in eachindex(forces))
    structural_work = mapped[5,8]

    force_aero = sum(A_TO_UVLM' * f for f in forces)
    moment_aero = sum(cross(A_TO_UVLM' * points[i], A_TO_UVLM' * forces[i])
        for i in eachindex(forces))
    force_beam = vec(sum(mapped[1:3,:]; dims=2))
    moment_beam = sum(cross(positions[:,n], mapped[1:3,n]) + mapped[4:6,n]
        for n in axes(positions,2))
    return (; aerodynamic_work, structural_work, force_aero, force_beam,
        moment_aero, moment_beam)
end

function deformed_work_check(nc)
    model, chord, spar, weights = audit_geometry()
    for (e, element) in enumerate(model.elements)
        element.nodalStates.u_n2 .= [0.03sin(e/9), 0.01cos(e/8), -0.04sin(e/7)]
        element.nodalStates.p_n2 .= [0.4sin(e/9), -0.3cos(e/8), 0.2sin(e/7)]
    end
    points, positions, offsets = audit_vertices(model, chord, spar, weights; nc)
    forces = [[sin(k+j), cos(2k+j), sin(k-2j)]
        for k in axes(points,1), j in axes(points,2)]
    loads = beam_loads((forces=[forces], positions=[points]), offsets, weights)

    # Independent finite difference of the full geometry, about finite rotations.
    # Physical moments are conjugate to delta-theta, not directly to WM delta-p.
    epsilon = 1e-6
    plus_model, minus_model = deepcopy(model), deepcopy(model)
    structural_work = 0.0
    for (e, element) in enumerate(model.elements)
        du = [0.01cos(e), 0.02sin(e), 0.01cos(2e)]
        dp = [0.2sin(e), 0.1cos(e), -0.3sin(2e)]
        plus_model.elements[e].nodalStates.u_n2 .+= epsilon .* du
        minus_model.elements[e].nodalStates.u_n2 .-= epsilon .* du
        plus_model.elements[e].nodalStates.p_n2 .+= epsilon .* dp
        minus_model.elements[e].nodalStates.p_n2 .-= epsilon .* dp
        R = first(AeroBeams.rotation_tensor_WM(element.nodalStates.p_n2))
        Rp = first(AeroBeams.rotation_tensor_WM(plus_model.elements[e].nodalStates.p_n2))
        Rm = first(AeroBeams.rotation_tensor_WM(minus_model.elements[e].nodalStates.p_n2))
        spin = ((Rp-Rm)/(2epsilon)) * R'
        dtheta = [spin[3,2]-spin[2,3], spin[1,3]-spin[3,1], spin[2,1]-spin[1,2]]/2
        structural_work += dot(loads[1:3,e+1], du) + dot(loads[4:6,e+1], dtheta)
    end
    plus, _, _ = audit_vertices(plus_model, chord, spar, weights; nc)
    minus, _, _ = audit_vertices(minus_model, chord, spar, weights; nc)
    aerodynamic_work = sum(dot(forces[i], (plus[i]-minus[i])/(2epsilon))
        for i in eachindex(forces))
    @test isapprox(structural_work, aerodynamic_work; atol=1e-9, rtol=1e-7)
    @test vec(sum(loads[1:3,:]; dims=2)) ≈ sum(A_TO_UVLM' * f for f in forces)
    moment_aero = sum(cross(A_TO_UVLM' * points[i], A_TO_UVLM' * forces[i])
        for i in eachindex(forces))
    moment_beam = sum(cross(positions[:,n], loads[1:3,n]) + loads[4:6,n]
        for n in axes(positions,2))
    @test moment_aero ≈ moment_beam

    # Offsets must reconstruct every vortex vertex, including the trailing edge.
    # Test a converged geometry, the FSI predictor, and a relaxed trial geometry.
    grid, _, _ = wing_geometry(model, chord, spar, nc, weights)
    next_grid, next_positions, next_offsets = wing_geometry(plus_model, chord, spar, nc, weights)
    for fraction in (0.0, 2.0, 0.3)
        trial_grid = (1-fraction).*grid .+ fraction.*next_grid
        trial_positions = (1-fraction).*positions .+ fraction.*next_positions
        trial_offsets = (1-fraction).*offsets .+ fraction.*next_offsets
        _, _, panels = UVLM.grid_to_surface_panels(trial_grid)
        trial_points = UVLM.imperial_nodal_positions(panels)
        reconstructed = [A_TO_UVLM * sum(weights[j,n] .* (trial_positions[:,n] + trial_offsets[:,k,n])
            for n in axes(weights,2)) for k in axes(points,1), j in axes(points,2)]
        @test all(isapprox.(reconstructed, trial_points; atol=1e-12))
    end
end

function retry_check()
    # Use the same numeric load BCs and active-model update as the Pazy loop.
    loads = zeros(6,5)
    loads[3,end] = 30.0
    beam = AeroBeams.create_Beam(length=1.0, nElements=4,
        S=[Matrix(Diagonal([1e5, 1e5, 1e5, 100.0, 10.0, 10.0]))])
    clamp = AeroBeams.create_BC(beam=beam, node=1,
        types=["u1A", "u2A", "u3A", "p1A", "p2A", "p3A"], values=zeros(6))
    bcs = [clamp]
    for node in 1:5
        push!(bcs, AeroBeams.create_BC(beam=beam, node=node,
            types=["F1A", "F2A", "F3A", "M1A", "M2A", "M3A"], values=zeros(6)))
    end
    model = AeroBeams.create_Model(beams=[beam], BCs=bcs)
    solver = AeroBeams.create_NewtonRaphson(maximumIterations=4,
        absoluteTolerance=1e-8, relativeTolerance=1e-8)
    problem = AeroBeams.create_SteadyProblem(model=model, systemSolver=solver)
    apply_pazy_nodal_loads!(problem, loads, 0.0)
    # This calls the same Newton retry implementation as solve_time_step!.
    AeroBeams.solve_NewtonRaphson!(problem)
    replaced = problem.model !== model
    converged = solver.convergedFinalSolution
    loads[3,end] = 31.0
    apply_pazy_nodal_loads!(problem, loads, 0.0)
    solver_value = problem.model.BCs[end].currentValue[3]
    special_node_value = problem.model.specialNodes[end].BCs[end].currentValue[3]
    previous_solution = copy(problem.x)
    AeroBeams.solve_NewtonRaphson!(problem)
    return (; replaced, converged, live_value=loads[3,end], solver_value,
        special_node_value, second_converged=solver.convergedFinalSolution,
        solution_change=norm(problem.x-previous_solution))
end

@testset "Pazy coupling audit (failures identify implementation defects)" begin
    @testset "AeroBeams and UVLM axis convention" begin
        model, chord, spar, weights = audit_geometry()
        grid, _, _ = wing_geometry(model, chord, spar, 4, weights)
        root_le_A = A_TO_UVLM' * grid[:, 1, 1]
        root_te_A = A_TO_UVLM' * grid[:, end, 1]

        # AeroBeams' Pazy airfoil points forward along +A-y. Therefore the
        # leading edge is ahead of the spar and the trailing edge is behind it.
        @test root_le_A ≈ [0.0, spar * chord, 0.0]
        @test root_te_A ≈ [0.0, -(1-spar) * chord, 0.0]
        @test grid[1, 1, 1] < grid[1, end, 1] # UVLM orders LE -> TE downstream.
        @test isapprox(det(A_TO_UVLM), 1.0; atol=1e-14)

        alpha = deg2rad(3.0)
        freestream_A = A_TO_UVLM' * [cos(alpha), 0.0, sin(alpha)]
        @test freestream_A ≈ [-sin(alpha), -cos(alpha), 0.0]
    end
    @testset "Coincident-grid control" begin
        check = virtual_work_check(Matrix{Float64}(I,16,16))
        @test check.force_aero ≈ check.force_beam
        @test check.moment_aero ≈ check.moment_beam
        @test isapprox(check.structural_work, check.aerodynamic_work; atol=1e-12)
    end
    @testset "Finite-rotation virtual work and vortex offsets" begin
        for nc in (1, 4, 7)
            deformed_work_check(nc)
        end
    end
    @testset "Current noncoincident-grid virtual work" begin
        _, _, _, weights = audit_geometry()
        check = virtual_work_check(weights)
        println("Virtual work per unit rotation: aero=", check.aerodynamic_work,
            ", structure=", check.structural_work, " N*m")
        @test check.force_aero ≈ check.force_beam
        @test check.moment_aero ≈ check.moment_beam
        @test isapprox(check.structural_work, check.aerodynamic_work; atol=1e-12)
    end
    @testset "Newton retry preserves externally updated load BCs" begin
        check = retry_check()
        println("Newton recovery: ", check)
        @test check.replaced # Confirm this exercise actually took the retry path.
        @test check.converged
        @test check.solver_value == check.live_value
        @test check.special_node_value == check.live_value
        @test check.second_converged
        @test check.solution_change > 1e-6
    end
end
