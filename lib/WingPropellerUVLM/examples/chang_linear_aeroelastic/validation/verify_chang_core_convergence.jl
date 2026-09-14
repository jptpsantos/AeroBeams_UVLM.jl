# Bounded end-to-end verification; these coarse grids are NOT production results.
using Test
include(joinpath(@__DIR__, "..", "studies", "convergence", "archive", "run_chang_core_convergence.jl"))
const C = ChangCoreConvergence
config = C.CaseModel.ChangAeroelastic.load_chang_configuration(; env=Dict{String,String}())
coarse = merge(config, (;
    wing=merge(config.wing,(; spanwise_panels=2,chordwise_panels=1)),
    propeller=merge(config.propeller,(; radial_panels=2,chordwise_panels=1)),
    simulation=merge(config.simulation,(; azimuth_step_deg=90.))))
settings = merge(C.core_convergence_defaults(coarse), (;
    dry_run=false, make_plots=true, reuse_existing=true,
    nominal_mesh=(wing_span=2,wing_chord=1,prop_radial=2,prop_chord=1),
    core_meshes=[(wing_span=2,wing_chord=1,prop_radial=2,prop_chord=1),
        (wing_span=4,wing_chord=2,prop_radial=4,prop_chord=2),
        (wing_span=6,wing_chord=3,prop_radial=6,prop_chord=3)],
    radii_m=[.02,.01,.005], nominal_radius_m=.01,
    families=[:core_mesh,:time_step], azimuth_levels_deg=[90.,45.,30.],
    nominal_wake_revolutions=1., simulated_revolutions=5, averaged_revolutions=1,
    output_directory=normpath(joinpath(@__DIR__,"..","output","core_convergence_smoke"))))
result = C.run_study(settings,coarse)
@testset "Fixed-core rigid solver and report smoke" begin
    @test length(result.results) == 21
    @test all(r -> r.status == "completed",result.results)
    @test all(r -> all(isfinite,r.values),result.results)
    @test count(r -> r.reused,result.results) >= 1
    @test all(file -> isfile(joinpath(result.directory,file)),
        ("core_sensitivity.png","convergence_errors.png","case_matrix.csv","convergence.csv","report.md"))
    c = first(result.cases)
    reused = C.rigid_run(c,settings,coarse,joinpath(result.directory,c.label))
    @test reused.reused
end

# Exercise the production-child interface without running a costly full damping study.
C.load_elastic_backend()
const E = C.Elastic
@testset "Fixed-core aeroelastic child configuration" begin
    c = E.UVLMConvergenceCase(family=:core,level=1,label="fixed",
        core_radius_m=.002, fcore_segment_factor=0.,fcore_chord_factor=0.)
    E.validate_case(c)
    environment = E.child_environment(c,result.directory,"smoke";
        speed_mps=85.,trim_rpm=1212.,trim_speed_mps=65.,end_time_s=5.,hard_angle_deg=Inf)
    @test environment["CHANG_CORE_RADIUS_M"] == "0.002"
    @test environment["CHANG_FCORE_SEGMENT_FACTOR"] == "0.0"
    @test environment["CHANG_WAKE_ROWS_WING"] == environment["CHANG_WAKE_ROWS_PROPELLER"]
    @test E.case_key(c) != E.case_key(E.UVLMConvergenceCase(family=:core,level=1,label="fixed",
        core_radius_m=.001,fcore_segment_factor=0.,fcore_chord_factor=0.))
    dry = merge(C.core_convergence_defaults(config), (;
        mode=:aeroelastic,make_plots=false,dry_run=true,
        output_directory=joinpath(result.directory,"aeroelastic_preview")))
    preview = C.run_study(dry,config)
    @test isempty(preview.results)
    @test !isempty(preview.cases)
    for dt in (.007,.011)
        time = collect(0.:dt:6.)
        signal = exp.(-.1time) .* sin.(2pi*4.2time)
        metrics = E.moving_block_metrics(time,signal; fit_start_s=.5,fit_end_s=6.,
            minimum_peaks=8,minimum_fit_r_squared=.8,initial_block_size=512,
            size_ratio_lower=.25,size_ratio_upper=.5,peak_from_start=1,peak_from_end=0,
            block_duration_s=1.,block_overlap=.9,frequency_min_hz=3.,frequency_max_hz=6.,
            apply_hann_window=true)
        @test metrics.valid_for_convergence
        @test isapprox(metrics.block_size*dt,1.;atol=dt/2+1e-12)
        @test isapprox(metrics.moving_block_lambda_per_s,-.1;atol=.01)
    end
end
