module MigrationSmoke
using Test, LinearAlgebra, TOML
import WingPropellerUVLM as U
import WingPropellerUVLMStudies as S
const EXAMPLE=normpath(joinpath(@__DIR__,"..","..","examples","chang_linear_aeroelastic"))
include(joinpath(EXAMPLE,"studies","convergence","chang_convergence_aerodynamic.jl"))
include(joinpath(EXAMPLE,"studies","convergence","chang_convergence_aeroelastic.jl"))
include(joinpath(EXAMPLE,"..","prescribed_wing.jl"))

@testset "Executable examples and Chang facade parity" begin
    @test PrescribedWingExample.main().results.completed
    C=S.Convergence
    a=ChangConvergenceAerodynamic.aerodynamic_study()
    e=ChangConvergenceAeroelastic.aeroelastic_study()
    controls=C.case(:core,1,(;wing_span=4,wing_chord=1,prop_radial=2,prop_chord=1,
        wake_revolutions=1.,core=.02,azimuth_deg=30.))
    e=merge(e,(;speed_mps=85.,make_plots=false,
        response=merge(e.response,(;end_time_s=.12,trim_revolutions=1.,
            trim_average_revolutions=.5,impulse_duration_s=.015,impulse_magnitude_nm=.01))))
    reference=(;physical=a.physical,selected=C.controls(controls),core_mode=:fixed)
    operating=S.operating_selection(reference,e)
    mktempdir() do folder
        settings=C.model_defaults(operating.physical,:fixed,controls,e,joinpath(folder,"new"))
        model=S.create_ChangModel(S.canonical_settings(settings))
        problem=S.create_ChangDynamicProblem(;model)
        S.solve!(problem)
        @test problem.results.completed && problem.results.solverConverged
        @test problem.results.model===model
        @test model.parameters.Vinf==85.
        @test model.parameters.Ω ≈ 1212*85/65*2pi/60
        @test any(q->norm(q)>0,problem.results.solution.displacement_history)
        oldSettings=merge(settings,(;output=merge(settings.output,(;directory=joinpath(folder,"old")))))
        M=C.Model.ChangAeroelastic
        config=Base.invokelatest(M.load_chang_configuration,oldSettings;env=Dict{String,String}())
        old=Base.invokelatest(M.run_chang,config)
        @test old.solution.displacement_history==problem.results.solution.displacement_history
        @test old.solution.velocity_history==problem.results.solution.velocity_history
        @test old.solution.coupling_iterations==problem.results.solution.coupling_iterations
        @test_throws ArgumentError S.solve!(problem)
        # A physical limit stopping the run must fail before damping acceptance.
        stopped=merge(e,(;integration=merge(e.integration,(;state_norm_limit=1e-20))))
        failedDirectory=joinpath(folder,"stopped")
        @test_throws ErrorException Base.invokelatest(C.aeroelastic_case,
            controls,stopped,operating,failedDirectory)
        status=TOML.parsefile(joinpath(failedDirectory,"analysis_status.toml"))
        @test !status["completed"] && !status["solverConverged"]
        @test status["acceptedSteps"] < status["expectedSteps"]
    end
end
end
