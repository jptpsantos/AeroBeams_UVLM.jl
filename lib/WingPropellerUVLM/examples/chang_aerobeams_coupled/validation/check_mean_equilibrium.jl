using Test,LinearAlgebra
include(joinpath(@__DIR__,"..","src","ChangAeroBeams.jl"))
using .ChangAeroBeams
c = chang_aerobeams_defaults()
ref = merge(
    c.reference,
    (
        wing=merge(c.reference.wing,(spanwise_panels=4,chordwise_panels=2)),
        propeller=merge(c.reference.propeller,(attachment_eta=[0.5],radial_panels=2,chordwise_panels=2)),
        simulation=merge(c.reference.simulation,(end_time_s=0.01,azimuth_step_deg=5.0)),
        wake=merge(c.reference.wake,(maximum_rows_wing=12,maximum_rows_propeller=12)))
)
c = merge(
    c,
    (reference=ref,
        initialization=merge(
            c.initialization,
            (relaxation=1.0,maximum_iterations=12,
                minimum_revolutions=3,maximum_revolutions=10,periodic_load_tolerance=2e-3)
        ),
        solver=merge(c.solver,(impulse_magnitude_nm=0.0,)),
        output=merge(c.output,(label="verification_mean_equilibrium",write_results=true,plot_results=true)))
)
@testset "Mean-loaded nonlinear equilibrium and developed-wake response" begin
    cp = run_chang_aerobeams(c)
    @test cp.completed
    @test haskey(cp.results,:mean_equilibrium)
    @test cp.phase_time_offset > 0
    @test abs(cp.results[:node_positions][3,cp.wing.wing_node_count,1]) > 1e-8
    @test cp.results[:mean_equilibrium].periodic_load_error <= c.initialization.periodic_load_tolerance
    @test maximum(cp.results[:load_residual]) <= c.solver.load_tolerance
    @test isfile(cp.results[:paths].history)
end
