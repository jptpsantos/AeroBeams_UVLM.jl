# Bounded real-solver checks. These coarse/short runs are NOT convergence evidence.
using Test
const ROOT=joinpath(@__DIR__,"..","studies","convergence")
include(joinpath(ROOT,"chang_convergence_aerodynamic.jl"))
include(joinpath(ROOT,"chang_convergence_aeroelastic.jl"))
const C=ChangConvergenceAerodynamic.Convergence
const OUTPUT=normpath(joinpath(@__DIR__,"..","output","independent_convergence_smoke"))
a=ChangConvergenceAerodynamic.aerodynamic_study()
small=merge(a,(;
    nominal=(wing_span=2,wing_chord=1,prop_radial=2,prop_chord=1,wake_revolutions=1.,core=.01,azimuth_deg=90.),
    sweeps=(core=[.02,.01,.005],),families=[:core],
    simulated_revolutions=5,averaged_revolutions=1,dry_run=false,make_plots=true,
    output_directory=joinpath(OUTPUT,"aerodynamic")))
result=C.run_aerodynamic(small)
@testset "Independent aerodynamic solver smoke" begin
    @test length(result.results)==3
    @test all(r->r.status=="completed",result.results)
    @test all(r->all(isfinite,r.values),result.results)
    @test isnothing(result.selection_file) # Partial sweeps cannot authorize stage 2.
    @test_throws ErrorException C.load_selection(joinpath(result.directory,"aerodynamic_selection.toml"))
    @test isfile(joinpath(result.directory,"convergence.png"))
    c=first(result.cases)
    @test C.aerodynamic_case(c,small,joinpath(result.directory,c.label)).reused
end

C.load_model()
e=ChangConvergenceAeroelastic.aeroelastic_study()
short=merge(e,(;response=merge(e.response,(;end_time_s=.02)),make_plots=false))
selection=(;physical=small.physical,core_mode=small.core_mode,selected=small.nominal)
c=C.case(:core,1,small.nominal)
# Make fallback to the interactive case fail visibly, even during time marching.
@eval C.Model.ChangAeroelastic chang_case_defaults()=error("Unexpected interactive-case dependency")
elastic_result=C.aeroelastic_case(c,short,selection,joinpath(OUTPUT,"aeroelastic"))
@testset "Independent production solver smoke" begin
    @test isfile(elastic_result.history)
    @test isfile(joinpath(OUTPUT,"aeroelastic","resolved_configuration.toml"))
    @test isfile(joinpath(OUTPUT,"aeroelastic","damping_fit.toml"))
    @test elastic_result.status=="indeterminate" # Too short for a damping estimate.
    @test !elastic_result.valid
end
