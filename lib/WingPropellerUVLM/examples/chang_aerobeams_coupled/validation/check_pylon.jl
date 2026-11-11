using Test,LinearAlgebra
include(joinpath(@__DIR__,"..","src","ChangAeroBeams.jl"))
using .ChangAeroBeams
const CA = ChangAeroBeams
c = chang_aerobeams_defaults()
ref = merge(
    c.reference,
    (
        wing=merge(c.reference.wing,(spanwise_panels=4,chordwise_panels=2)),
        propeller=merge(c.reference.propeller,(attachment_eta=[0.5],radial_panels=2,chordwise_panels=2)),
        simulation=merge(c.reference.simulation,(end_time_s=0.002,azimuth_step_deg=5.0)),
        wake=merge(c.reference.wake,(maximum_rows_wing=8,maximum_rows_propeller=8)))
)
c = merge(
    c,
    (reference=ref,solver=merge(c.solver,(support=:beam_pylon,impulse_magnitude_nm=0.0)),
        pylon=merge(c.pylon,(elements=2,)),initialization=merge(c.initialization,(mode=:unloaded,)),
        output=merge(c.output,(write_results=false,plot_results=false)))
)

@testset "Nonlinear beam pylon and physical rotor" begin
    cp = build_chang_aerobeams(c)
    @test isempty(cp.eta)
    @test cp.rotors[1] isa TipRotor
    @test cp.wing.node_count == cp.wing.wing_node_count + 2
    hub = CA.pylon_hub_poses(cp.wing.problem,cp.wing.pylon_tip_nodes)
    pivot = CA.wing_node_poses(cp.wing.problem).positions[:,cp.reference.parameters.prop_attach_nodes[1]]
    @test hub.positions[:,1] ≈ pivot + [-cp.reference.parameters.L_pylon,0.0,0.0] atol = 1e-8
    @test hub.rotations[1] ≈ I atol = 1e-10
    solve_chang_aerobeams!(cp)
    @test cp.completed
    @test maximum(cp.results[:geometry_residual]) <= c.solver.geometry_tolerance
    @test maximum(cp.results[:load_residual]) <= c.solver.load_tolerance
    @test all(isfinite,cp.results[:node_positions])
    @test any(abs.(cp.results[:node_positions][3,:,:]) .> 0)
    @test cp.aerodynamic.state.stepIndex == length(cp.time)
end
