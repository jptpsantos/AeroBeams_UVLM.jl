# Bounded dynamic and virtual-work checks. No convergence or flutter claim.
using Test, LinearAlgebra, TOML
pushfirst!(DEPOT_PATH, mktempdir())
const ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const STUDY = joinpath(ROOT, "lib", "WingPropellerUVLM", "examples",
    "chang_linear_aeroelastic", "studies", "convergence")
Base.include(ex -> ex == :(ChangConvergenceAerodynamic.main()) ? nothing : ex,
    @__MODULE__, joinpath(STUDY, "chang_convergence_aerodynamic.jl"))
include(joinpath(STUDY, "chang_convergence_aeroelastic.jl"))
const E = ChangConvergenceAeroelastic
const C = E.Convergence
C.load_model()
const M = C.Model.ChangAeroelastic
a = ChangConvergenceAerodynamic.aerodynamic_study()
e = E.aeroelastic_study()
c = C.case(:core, 1, (;wing_span=4, wing_chord=1, prop_radial=2, prop_chord=1,
    wake_revolutions=1.0, core=0.02, azimuth_deg=30.0))
e = merge(e, (;speed_mps=85.0, make_plots=false,
    response=merge(e.response, (;end_time_s=0.12, trim_revolutions=1.0,
        trim_average_revolutions=0.5, impulse_duration_s=0.015,
        impulse_magnitude_nm=0.01))))
reference = (;physical=a.physical, selected=C.controls(c), core_mode=:fixed)
operating = E.operating_selection(reference, e)
defaults = C.model_defaults(operating.physical, :fixed, c, e, joinpath(@__DIR__, "smoke_output"))
config = M.load_chang_configuration(defaults; env=Dict{String,String}())
run = M.run_chang(config)
@testset "Bounded coupled speed-override run" begin
    @test run.model.parameters.Vinf == 85.0
    @test run.model.parameters.Ω ≈ (1212.0 * 85/65) * 2pi/60
    @test run.solution.last_step == length(run.model.parameters.dt)
    @test all(run.solution.coupling_converged)
    @test any(>(1), run.solution.coupling_iterations)
    @test any(q -> norm(q)>0, run.solution.displacement_history)
    @test all(q -> all(isfinite,q), run.solution.displacement_history)
    @test run.model.structural.C ≈ -run.model.structural.C'
end

interacting_config=merge(config,(;
    simulation=merge(config.simulation,(;interaction_on=true)),
    output=merge(config.output,(;directory=joinpath(@__DIR__,"smoke_interacting_output")))))
interacting_run=M.run_chang(interacting_config)
@testset "Bounded interacting wing-propeller run" begin
    @test interacting_run.solution.last_step==length(interacting_run.model.parameters.dt)
    @test all(interacting_run.solution.coupling_converged)
    @test all(q -> all(isfinite,q),interacting_run.solution.displacement_history)
    @test !isapprox(last(interacting_run.solution.displacement_history),last(run.solution.displacement_history))
end

# Check actual vortex positions and the same frozen loads used in the adapter.
# Nonuniform wing rotations exercise interpolation at the propeller attachment.
q = copy(last(run.solution.displacement_history))
for i in eachindex(q)
    q[i] = 0.015sin(i)
end
model, workspace = run.model, run.workspace
time = last(model.parameters.t)
kinematics = M.update_aero_geometry_for_state!(model,workspace,q,time)
load = M.assemble_structural_aero_load!(model,workspace,kinematics)
forces = M.imperial_nodal_forces(workspace.system)
errors = Float64[]
for dof in eachindex(q)
    h=1e-7
    plus=copy(q); minus=copy(q)
    plus[dof]+=h; minus[dof]-=h
    M.update_aero_geometry_for_state!(model,workspace,plus,time)
    xp=M.imperial_nodal_positions(workspace.system)
    M.update_aero_geometry_for_state!(model,workspace,minus,time)
    xm=M.imperial_nodal_positions(workspace.system)
    work=sum(dot(forces[s][v],(xp[s][v]-xm[s][v])/(2h))
        for s in eachindex(forces) for v in eachindex(forces[s]))
    push!(errors,abs(work-load[dof])/max(1.0,abs(work),abs(load[dof])))
end
M.update_aero_geometry_for_state!(model,workspace,q,time)
@testset "All free DOFs: frozen-load virtual work" begin
    @test maximum(errors)<1e-6
end
panel_bytes=sizeof(eltype(workspace.system.surfaces[1]))
summary=Dict("julia_version"=>string(VERSION), "steps"=>run.solution.last_step,
    "source_sha256"=>C.SOURCE(),
    "interacting_steps"=>interacting_run.solution.last_step,
    "max_coupling_iterations"=>maximum(run.solution.coupling_iterations),
    "maximum_virtual_work_relative_error"=>maximum(errors),
    "surface_panel_bytes"=>panel_bytes,
    "retained_surface_history_bytes"=>Base.summarysize(workspace.surface_history))
println(summary)
open(joinpath(@__DIR__,"smoke_summary.toml"),"w") do io
    TOML.print(io,summary)
end
