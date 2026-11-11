using Test
using LinearAlgebra
using Random

# Run in isolation from the structural solver: these tests exercise geometry,
# virtual work, aerodynamic transactions and parity with the existing Chang
# geometry using frozen forces. They do not march a structural simulation.
module AerodynamicInterfaceHarness
using LinearAlgebra
const REPOSITORY_ROOT = normpath(joinpath(@__DIR__,"..","..","..","..",".."))
const UVLM_ROOT = normpath(joinpath(@__DIR__,"..","..",".."))
for root in (REPOSITORY_ROOT,UVLM_ROOT)
    root in LOAD_PATH || pushfirst!(LOAD_PATH,root)
end

import WingPropellerUVLM as UVLM
include(joinpath(@__DIR__,"..","..","chang_linear_aeroelastic","src","ChangAeroelastic.jl"))
const Chang = ChangAeroelastic
include(joinpath(@__DIR__,"..","src","AerodynamicInterface.jl"))
end

const CAI = AerodynamicInterfaceHarness

function interface_reference()
    defaults = CAI.Chang.chang_case_defaults()
    defaults = merge(
        defaults,
        (;
            wing=merge(defaults.wing,(; spanwise_panels=4,chordwise_panels=2)),
            propeller=merge(defaults.propeller,(; attachment_eta=[0.5],radial_panels=3,chordwise_panels=2)),
            simulation=merge(defaults.simulation,(; end_time_s=0.02,interaction_on=true)),
            output=merge(defaults.output,(; plot_results=false,animate_wake=false)))
    )
    return CAI.Chang.build_chang_model(CAI.Chang.load_chang_configuration(defaults; env=Dict{String,String}()))
end

function reference_poses(reference)
    p = reference.parameters
    positions = hcat(
        (
            [p.xle_distribution[i] + reference.aerodynamic_options.elastic_axis_fraction * p.chord[i],
                p.span_nodes[i],0.0] for i = 1:p.nnodes
        )...
    )
    return (; positions,rotations=[copy(CAI.CHANG_DIRECTORS_TO_U) for _ = 1:p.nnodes])
end

interface_skew(v) = [0.0 -v[3] v[2]; v[3] 0.0 -v[1]; -v[2] v[1] 0.0]

function vortex_positions(grids)
    return [CAI.UVLM.imperial_nodal_positions(CAI.UVLM.grid_to_surface_panels(grid)[3]) for grid in grids]
end

@testset "Chang nonlinear aerodynamic interface" begin
    reference = interface_reference()
    p = reference.parameters
    poses = reference_poses(reference)
    eta = deg2rad.([4.0,3.0])
    time = 0.0123

    @testset "Existing Chang rigid geometry parity" begin
        workspace = CAI.Chang.build_chang_workspace(reference)
        state = zeros(reference.structural.ndof_free)
        state[reference.structural.ndof_wing_free+1:end] .= eta
        CAI.Chang.update_aero_geometry_for_state!(reference,workspace,state,time)
        grids,geometry = CAI.coupled_grids(reference,poses,eta,time)
        @test vortex_positions(grids)[1] ≈
              CAI.UVLM.imperial_nodal_positions(workspace.system.surfaces[1]) atol = 1e-12
        for blade = 1:p.Nb_prop
            @test grids[1+blade] ≈ workspace.grids_prop_current[1][blade] atol = 1e-12
        end
        @test norm(geometry.hubs[1] - geometry.pivots[1]) >
              norm(geometry.hubs[1] - geometry.centers[1])
    end

    @testset "Finite rotations, quarter-vortex work and wrench balance" begin
        rng = MersenneTwister(92031)
        finite_poses = (;
            positions=poses.positions + 0.04randn(rng,size(poses.positions)),
            rotations=[exp(interface_skew(0.2randn(rng,3))) * r for r in poses.rotations])
        grids,geometry = CAI.coupled_grids(reference,finite_poses,eta,time)
        points = vortex_positions(grids)
        forces = [map(_ -> randn(rng,3),surface) for surface in points]
        loads = CAI.aerodynamic_loads(reference,finite_poses,eta,time,forces,points; geometry)
        @test geometry.vortex_fractions[1,1] != geometry.chord_fractions[1,1]
        @test geometry.vortex_fractions[end,:] == geometry.chord_fractions[end,:]
        translation = randn(rng,size(poses.positions))
        rotation = randn(rng,size(poses.positions))
        delta_eta = randn(rng,length(eta))
        epsilon = 1e-7
        perturbed = []
        for sign in (1.0,-1.0)
            trial = (; positions=finite_poses.positions + sign * epsilon * translation,
                rotations=[
                    exp(sign * epsilon * interface_skew(rotation[:,i])) * finite_poses.rotations[i]
                    for i in eachindex(finite_poses.rotations)
                ])
            moved,_ = CAI.coupled_grids(reference,trial,eta + sign * epsilon * delta_eta,time)
            push!(perturbed,vortex_positions(moved))
        end
        aerodynamic_work = sum(
            dot(forces[s][i],(perturbed[1][s][i] - perturbed[2][s][i]) / (2epsilon))
            for s in eachindex(forces) for i in eachindex(forces[s])
        )
        structural_work =
            sum(loads.wing_loads[1:3,:] .* translation) +
            sum(loads.wing_loads[4:6,:] .* rotation) + dot(loads.rotor_loads,delta_eta)
        @test aerodynamic_work ≈ structural_work rtol = 2e-7 atol = 2e-7
        force_total = sum((f for surface in forces for f in surface))
        moment_total = sum((cross(points[s][i],forces[s][i])
                            for s in eachindex(forces) for i in eachindex(forces[s])))
        wing_moment = sum((
            loads.wing_loads[4:6,i] +
            cross(finite_poses.positions[:,i],loads.wing_loads[1:3,i]) for i = 1:p.nnodes
        ))
        @test vec(sum(loads.wing_loads[1:3,:]; dims=2)) ≈ force_total atol = 1e-12
        @test wing_moment ≈ moment_total atol = 1e-11
    end

    @testset "Independent aerodynamic and structural span meshes" begin
        span = collect(range(0,p.span_length; length=9))
        positions = hcat(
            (
                [_ca_x,station,0.0] for (station,_ca_x) in
                zip(
                    span,
                    [
                        CAI._ca_interp(p.span_nodes,p.xle_distribution,y) +
                        reference.aerodynamic_options.elastic_axis_fraction * CAI._ca_interp(p.span_nodes,p.chord,y)
                        for y in span
                    ]
                )
            )...
        )
        fine = (; positions,rotations=[copy(CAI.CHANG_DIRECTORS_TO_U) for _ in span])
        coarse_grid,_ = CAI.coupled_grids(reference,poses,eta,time)
        fine_grid,_ = CAI.coupled_grids(reference,fine,eta,time; structural_span_nodes=span)
        @test all(coarse_grid[i] ≈ fine_grid[i] for i in eachindex(coarse_grid))
        # Coarser structural stations require genuine two-node interpolation.
        # Check its transpose under finite, nonuniform director rotations.
        rng = MersenneTwister(184)
        nodes = [1,3,5]
        span = p.span_nodes[nodes]
        mixed = (; positions=poses.positions[:,nodes] + 0.03randn(rng,3,3),
            rotations=[exp(interface_skew(0.25randn(rng,3))) * CAI.CHANG_DIRECTORS_TO_U for _ in nodes])
        grid,geometry = CAI.coupled_grids(reference,mixed,eta,time; structural_span_nodes=span)
        points = vortex_positions(grid)
        forces = [map(_ -> randn(rng,3),surface) for surface in points]
        load = CAI.aerodynamic_loads(reference,mixed,eta,time,forces,points; geometry)
        @test count(>(0),geometry.weights[:,2]) == 2
        dr,dw,dq = randn(rng,3,3),randn(rng,3,3),randn(rng,2)
        epsilon = 1e-7
        moved = []
        for sign in (1.0,-1.0)
            trial = (; positions=mixed.positions + sign * epsilon * dr,
                rotations=[exp(sign * epsilon * interface_skew(dw[:,i])) * mixed.rotations[i] for i = 1:3])
            g,_ = CAI.coupled_grids(reference,trial,eta + sign * epsilon * dq,time; structural_span_nodes=span)
            push!(moved,vortex_positions(g))
        end
        aero_work = sum(
            dot(forces[s][i],(moved[1][s][i] - moved[2][s][i]) / (2epsilon))
            for s in eachindex(forces) for i in eachindex(forces[s])
        )
        struct_work =
            sum(load.wing_loads[1:3,:] .* dr) +
            sum(load.wing_loads[4:6,:] .* dw) + dot(load.rotor_loads,dq)
        @test aero_work ≈ struct_work rtol = 2e-7 atol = 2e-7
    end

    @testset "Unique surface IDs and deterministic aerodynamic rollback" begin
        problem = CAI.build_aerodynamic(reference,poses,zeros(2),[0.0,0.0001,0.0002];
            maximum_wake_rows=3,tracking_time_steps=false)
        expected_ids = collect(1:length(problem.model.grids))
        @test problem.model.interactionGroups == vcat(1,fill(2,p.Nb_prop))
        for sys in (problem.state.system,problem.workspace.system,problem.workspace.commitSystem)
            @test sys.surface_id == expected_ids
        end
        grids,_ = CAI.coupled_grids(reference,poses,eta,0.0001)
        accepted = deepcopy(problem.state.system.Γ)
        CAI.UVLM.begin_time_step!(problem)
        first_load = CAI.UVLM.evaluate_trial!(problem; surfaces=grids)
        first_gamma = copy(problem.workspace.system.Γ)
        second_load = CAI.UVLM.evaluate_trial!(problem; surfaces=grids)
        @test problem.workspace.system.Γ == first_gamma
        @test first_load.forces == second_load.forces
        @test problem.state.system.Γ == accepted
        @test problem.workspace.system.surface_id == expected_ids
        CAI.UVLM.rollback_time_step!(problem)
        @test problem.state.timeNow == 0.0
        @test all(iszero,problem.state.system.nwake)
    end
end

@testset "Physical beam-pylon tip virtual work" begin
    reference = interface_reference()
    poses = reference_poses(reference)
    p = reference.parameters
    pivot = poses.positions[:,p.prop_attach_nodes[1]]
    hubs = (positions=reshape(pivot + [-p.L_pylon,0.0,0.0],3,1),
        rotations=[exp(interface_skew([0.12,-0.08,0.06]))])
    time = 0.0123
    grids,geometry = CAI.pylon_coupled_grids(reference,poses,hubs,time)
    points = vortex_positions(grids)
    rng = MersenneTwister(817)
    forces = [map(_ -> randn(rng,3),surface) for surface in points]
    mapped = CAI.pylon_aerodynamic_loads(reference,poses,hubs,time,forces,points; geometry)
    dr,dw = randn(rng,3),randn(rng,3)
    epsilon = 1e-7
    moved = []
    for sign in (1.0,-1.0)
        trial = (positions=hubs.positions .+ sign * epsilon .* reshape(dr,3,1),
            rotations=[exp(sign * epsilon * interface_skew(dw)) * hubs.rotations[1]])
        g,_ = CAI.pylon_coupled_grids(reference,poses,trial,time)
        push!(moved,vortex_positions(g))
    end
    aero_work = sum(
        dot(forces[s][i],(moved[1][s][i] - moved[2][s][i]) / (2epsilon))
        for s in eachindex(forces) for i in eachindex(forces[s])
    )
    structural_work = dot(mapped.pylon_loads[1:3,1],dr) + dot(mapped.pylon_loads[4:6,1],dw)
    @test aero_work ≈ structural_work rtol = 2e-7 atol = 2e-7
    @test isempty(mapped.rotor_loads)
    rotor_only = deepcopy(forces)
    fill!(rotor_only[1],[0.0,0.0,0.0])
    mapped_only = CAI.pylon_aerodynamic_loads(reference,poses,hubs,time,rotor_only,points; geometry)
    @test iszero(mapped_only.wing_loads)
    @test norm(mapped_only.pylon_loads) > 0
end
