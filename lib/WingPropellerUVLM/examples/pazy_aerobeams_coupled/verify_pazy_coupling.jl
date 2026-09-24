# Read-only numerical audit of the current coupling. No simulation outputs are written.
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

function audit_vertices(model, chord, spar, weights)
    grid, positions = wing_geometry(model, chord, spar, 4, weights)
    _, _, panels = UVLM.grid_to_surface_panels(grid)
    return UVLM.imperial_nodal_positions(panels), positions
end

function virtual_work_check(weights)
    model, chord, spar, _ = audit_geometry()
    points, positions = audit_vertices(model, chord, spar, weights)
    forces = [zeros(3) for k in axes(points,1), j in axes(points,2)]
    forces[2,8] = A_TO_UVLM * [1.0, 0.0, 0.0] # 1 N normal force.
    mapped = beam_loads((forces=[forces], positions=[points]), positions, weights)

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

function retry_check()
    # Use the same captured-array load BC pattern as run_pazy_wing_uvlm.
    load = [30.0]
    beam = AeroBeams.create_Beam(length=1.0, nElements=4,
        S=[Matrix(Diagonal([1e5, 1e5, 1e5, 100.0, 10.0, 10.0]))])
    clamp = AeroBeams.create_BC(beam=beam, node=1,
        types=["u1A", "u2A", "u3A", "p1A", "p2A", "p3A"], values=zeros(6))
    force = AeroBeams.create_BC(beam=beam, node=5, types=["F3A"],
        values=[t -> load[1]])
    model = AeroBeams.create_Model(beams=[beam], BCs=[clamp, force])
    solver = AeroBeams.create_NewtonRaphson(maximumIterations=4,
        absoluteTolerance=1e-8, relativeTolerance=1e-8)
    problem = AeroBeams.create_SteadyProblem(model=model, systemSolver=solver)
    # This calls the same Newton retry implementation as solve_time_step!.
    AeroBeams.solve_NewtonRaphson!(problem)
    replaced = problem.model !== model
    load[1] = 31.0
    AeroBeams.update_BC_data!(model.BCs[end], 0.0)
    AeroBeams.update_BC_data!(problem.model.BCs[end], 0.0)
    return (; replaced, converged=solver.convergedFinalSolution,
        live_value=model.BCs[end].currentValue[3],
        solver_value=problem.model.BCs[end].currentValue[3])
end

@testset "Pazy coupling audit (failures identify implementation defects)" begin
    @testset "Coincident-grid control" begin
        check = virtual_work_check(Matrix{Float64}(I,16,16))
        @test check.force_aero ≈ check.force_beam
        @test check.moment_aero ≈ check.moment_beam
        @test isapprox(check.structural_work, check.aerodynamic_work; atol=1e-12)
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
    end
end
