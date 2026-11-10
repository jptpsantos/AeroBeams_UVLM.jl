using Test, LinearAlgebra, DelimitedFiles
include(joinpath(@__DIR__,"..","src","ChangAeroBeams.jl"))
using .ChangAeroBeams
BLAS.set_num_threads(1)
ENV["GKSwstype"] = "100"
summaries = Any[]

@testset "Perturbed Chang responses in both support formulations" begin
    for support in (:reduced,:beam_pylon)
        c = chang_aerobeams_defaults()
        reference = merge(c.reference,(
            wing=merge(c.reference.wing,(spanwise_panels=8,chordwise_panels=3,)),
            propeller=merge(c.reference.propeller,(
                attachment_eta=[0.5],radial_panels=3,chordwise_panels=2,)),
            simulation=merge(c.reference.simulation,(
                end_time_s=0.25,azimuth_step_deg=5.0,interaction_on=true,)),
            wake=merge(c.reference.wake,(maximum_rows_wing=36,maximum_rows_propeller=72,)),))
        c = merge(c,(;reference,
            solver=merge(c.solver,(;support,impulse_start_s=0.05,impulse_duration_s=0.02,
                impulse_magnitude_nm=2.0,maximum_coupling_iterations=40,progress_frequency=50)),
            # The long free-wake check is an unloaded-start transient. Its
            # rigid-shape periodic trim is not assumed to exist/converge.
            # Mean equilibrium is checked independently in runtests.jl.
            initialization=merge(c.initialization,(mode=:unloaded,)),
            output=merge(c.output,(label="verification_"*String(support),)),))
        cp = run_chang_aerobeams(c)
        r = cp.results
        @test cp.completed
        @test maximum(r[:geometry_residual]) <= c.solver.geometry_tolerance
        @test maximum(r[:velocity_residual]) <= c.solver.velocity_tolerance
        @test maximum(r[:load_residual]) <= c.solver.load_tolerance
        @test maximum(r[:mechanical_residual]) <= c.solver.coupled_equilibrium_tolerance
        @test all(isfinite,r[:node_positions])
        @test length(cp.time) > 300
        tip = cp.wing.wing_node_count
        push!(summaries,[String(support),last(cp.time),length(cp.time)-1,
            r[:node_positions][3,tip,1],minimum(r[:node_positions][3,tip,:]),
            maximum(r[:node_positions][3,tip,:]),
            maximum(abs,rad2deg.(r[:rotor_angles][1,:])),
            maximum(abs,rad2deg.(r[:rotor_angles][2,:])),
            maximum(r[:coupling_iterations]),maximum(r[:load_residual]),
            maximum(r[:mechanical_residual])])
    end
end
directory = joinpath(@__DIR__,"..","output")
open(joinpath(directory,"perturbed_response_summary.csv"),"w") do io
    println(io,"support,duration_s,steps,initial_tip_z_m,min_tip_z_m,max_tip_z_m,max_abs_pitch_deg,max_abs_yaw_deg,max_fsi_iterations,max_load_residual,max_mechanical_residual")
    for row in summaries
        writedlm(io,permutedims(row),',')
    end
end
