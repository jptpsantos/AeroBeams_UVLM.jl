using Test
@testset "Canonical convergence settings and dry run" begin
    import WingPropellerUVLMStudies as Studies
    entries=joinpath(@__DIR__,"..","examples","chang_linear_aeroelastic","studies","convergence")
    scope=Module(:MigratedEntryTests)
    Base.include(scope,joinpath(entries,"chang_convergence_aerodynamic.jl"))
    Base.include(scope,joinpath(entries,"chang_convergence_aeroelastic.jl"))
    a=scope.ChangConvergenceAerodynamic.create_study()
    e=scope.ChangConvergenceAeroelastic.create_study()
    @test a.results===nothing && e.results===nothing
    @test Studies.legacy_settings(a)==scope.ChangConvergenceAerodynamic.aerodynamic_study()
    @test Studies.canonical_settings(Studies.legacy_settings(a))==Studies.study_settings(a)
    @test Studies.canonical_settings(Studies.legacy_settings(e))==Studies.study_settings(e)
    Base.include(scope,joinpath(entries,"..","..","chang_case.jl"))
    @test Studies.legacy_case_settings(Studies.canonical_settings(scope.chang_case_defaults()))==scope.chang_case_defaults()
    @test_throws ArgumentError Studies.canonical_settings((;speed_mps=65.,airspeed=85.))
    @test_throws ArgumentError Studies.create_AerodynamicConvergenceStudy(a.settings;typo=true)
    changed=Studies.create_AerodynamicConvergenceStudy(a.settings;physical=(;airspeed=70.))
    @test changed.settings.physical.airspeed==70.
    @test changed.settings.physical.rpm==a.settings.physical.rpm
    mktempdir() do folder
        dry=Studies.create_AerodynamicConvergenceStudy(a.settings;dryRun=true,makePlots=false,outputDirectory=folder)
        Studies.solve!(dry)
        @test isfile(joinpath(folder,"case_matrix.csv"))
        @test !isfile(joinpath(folder,"aerodynamic_selection.toml"))
        @test dry.results!==nothing
    end
end
