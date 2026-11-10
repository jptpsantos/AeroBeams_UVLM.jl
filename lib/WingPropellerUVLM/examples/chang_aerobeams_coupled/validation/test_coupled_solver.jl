using Test, LinearAlgebra
include(joinpath(@__DIR__,"..","src","ChangAeroBeams.jl"))
using .ChangAeroBeams
const CA = ChangAeroBeams

function small_config(;duration=0.002,impulse=0.0,initial_eta=nothing,interaction=true)
    c = chang_aerobeams_defaults()
    ref = merge(c.reference,(
        wing=merge(c.reference.wing,(spanwise_panels=4,chordwise_panels=2,)),
        propeller=merge(c.reference.propeller,(
            attachment_eta=[0.5],radial_panels=2,chordwise_panels=2,)),
        simulation=merge(c.reference.simulation,(
            end_time_s=duration,azimuth_step_deg=5.0,interaction_on=interaction,)),
        wake=merge(c.reference.wake,(maximum_rows_wing=8,maximum_rows_propeller=8,)),
    ))
    return merge(c,(reference=ref,
        initialization=merge(c.initialization,(mode=:unloaded,)),
        solver=merge(c.solver,(impulse_magnitude_nm=impulse,initial_eta=initial_eta,
            maximum_coupling_iterations=40,coupling_relaxation=0.7,progress_frequency=1000,)),
        output=merge(c.output,(write_results=false,plot_results=false,)),))
end

@testset "Augmented beam/rotor Jacobian" begin
    cp = build_chang_aerobeams(small_config())
    p = cp.wing.problem
    CA.AB.update_time_variables!(p,2)
    CA.AB.get_equivalent_states_rates!(p)
    z = vcat(p.x,cp.eta)
    loads = zeros(6,cp.wing.node_count)
    rl = zeros(length(cp.eta))
    args = (cp.node_history,cp.eta,cp.eta_velocity,cp.eta_acceleration,p.Δt)
    r,J = CA.coupled_mechanical_arrays!(cp,z,loads,rl,args...)
    @test norm(r,Inf) < 1e-8
    # Independent central difference of the COMPLETE residual, not just the
    # component helper used when assembling the augmented Jacobian.
    direction = zeros(length(z))
    sn = CA._special_node(p,cp.rotors[1].attachment_node)
    direction[sn.DOF_uF] .= [0.02,-0.01,0.03]
    direction[sn.DOF_pM] .= [0.01,0.02,-0.01]
    direction[end-1:end] .= [0.015,-0.025]
    h = 1e-6
    rp = CA.coupled_mechanical_arrays!(cp,z.+h.*direction,loads,rl,args...;jacobian=false)
    rm = CA.coupled_mechanical_arrays!(cp,z.-h.*direction,loads,rl,args...;jacobian=false)
    fd = (rp.-rm)./(2h)
    @test norm(J*direction-fd,Inf)/max(norm(fd,Inf),1.0) < 2e-5
end

@testset "Interacting nonlinear Chang response" begin
    cp = run_chang_aerobeams(small_config())
    @test cp.completed
    @test all(isfinite,cp.results[:node_positions])
    @test maximum(cp.results[:geometry_residual]) <= cp.config.solver.geometry_tolerance
    @test maximum(cp.results[:load_residual]) <= cp.config.solver.load_tolerance
    @test maximum(cp.results[:velocity_residual]) <= cp.config.solver.velocity_tolerance
    @test maximum(cp.results[:mechanical_residual]) <= cp.config.solver.coupled_equilibrium_tolerance
    @test cp.aerodynamic.state.stepIndex == length(cp.time)
    @test any(abs.(cp.results[:eta]) .> 0)
    @test any(abs.(cp.results[:node_positions][3,:,:]) .> 0)
    @test all(cp.aerodynamic.state.system.nwake .== length(cp.time)-1)
    @test !cp.wing.problem.getLinearSolution
end

@testset "Growth-rate fit on a known perturbation" begin
    time = collect(0:0.001:3.0)
    signal = exp.(-0.7.*time).*sin.(2pi*20 .* time)
    fit = perturbation_growth_rate(time,signal;start_time=0.2,window_duration=0.1)
    @test fit.growth_rate_per_s ≈ -0.7 atol=2e-3
    @test fit.r_squared > 0.999
end

@testset "Failed step leaves accepted clocks and state" begin
    c = small_config()
    c = merge(c,(solver=merge(c.solver,(maximum_coupling_iterations=1,)),))
    cp = build_chang_aerobeams(c)
    accepted_x = copy(cp.wing.problem.x)
    accepted_clock = cp.wing.problem.timeNow
    @test_throws ErrorException solve_chang_aerobeams!(cp)
    @test cp.wing.problem.x == accepted_x
    @test cp.wing.problem.timeNow == accepted_clock
    @test cp.aerodynamic.state.stepIndex == 1
    @test all(iszero,cp.aerodynamic.state.system.nwake)
    @test !cp.aerodynamic.workspace.stepOpen
end
