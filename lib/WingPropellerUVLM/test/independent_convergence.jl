module IndependentConvergenceTests
using Test,DelimitedFiles,TOML
const ROOT = joinpath(@__DIR__,"..","examples","chang_linear_aeroelastic","studies","convergence")
# Load definitions even if the interactive entry's execution guard is disabled.
Base.include(ex -> ex == :(ChangConvergenceAerodynamic.main()) ? nothing : ex,
    @__MODULE__,joinpath(ROOT,"chang_convergence_aerodynamic.jl"))
Base.include(ex -> ex == :(ChangConvergenceAeroelastic.main()) ? nothing : ex,
    @__MODULE__,joinpath(ROOT,"chang_convergence_aeroelastic.jl"))
const C = ChangConvergenceAerodynamic.Convergence
const E = ChangConvergenceAeroelastic

function fixture()
    s = ChangConvergenceAerodynamic.aerodynamic_study()
    return merge(
        s,
        (; dry_run=true,make_plots=false,families=collect(C.FIELDS),
            physical=merge(s.physical,(; speed_mps=65.0,wake_shedding_fraction=0.1)),core_mode=:fixed,
            nominal=(
                wing_span=2,
                wing_chord=1,
                prop_radial=2,
                prop_chord=1,
                wake_revolutions=1.0,
                core=0.01,
                azimuth_deg=90.0
            ),
            sweeps=(wing_span=[2,3,4],wing_chord=[1,2,3],prop_radial=[2,3,4],prop_chord=[1,2,3],
                wake_revolutions=[1.0,2.0,3.0],core=[0.02,0.01,0.005],azimuth_deg=[90.0,45.0,30.0]),
            simulated_revolutions=6,averaged_revolutions=1)
    )
end

function elastic_fixture()
    # Planning tests own their controls rather than the editable entry preset.
    return merge(
        E.aeroelastic_study(runMode=:full),
        (; dry_run=true,make_plots=false,
            speed_mps=nothing,nominal=nothing,families=collect(C.FIELDS[2:end]),
            sweeps=(wing_chord=nothing,prop_radial=nothing,prop_chord=nothing,
                wake_revolutions=nothing,core=nothing,azimuth_deg=nothing))
    )
end

@testset "Independent two-stage convergence" begin
    s = fixture()
    cases = C.validate_aerodynamic(s)
    @test length(cases) == 21
    @test length(unique(C.key.(cases))) == 15
    @test C.metric_names(false) == ("CT","CQ")
    @test C.core_factors(:fixed,0.01) == (radius=0.01,chord=0.0,segment=0.0)
    @test C.core_factors(:chord,0.01) == (radius=nothing,chord=0.01,segment=0.0)
    @test C.core_factors(:segment,0.01) == (radius=nothing,chord=0.0,segment=0.01)
    @test_throws ErrorException C.core_factors(:unknown,0.1)
    @test_throws ErrorException C.build_cases(s.nominal,merge(s.sweeps,(; azimuth_deg=[90.0,45.0])),s.families)
    @test_throws ErrorException C.build_cases(s.nominal,merge(s.sweeps,(; wing_span=[2,2,4])),s.families)
    wave = zeros(12)
    record(c) = (; case=c,status="completed",reason="fixture",reused=false,elapsed_s=0.0,
        values=(0.02,0.003),phases=(copy(wave),copy(wave)),periodic=(0.0,0.0),valid=true,history="",log="")
    results = record.(cases)
    rows = C.annotate(results,cases,s)
    @test all(r -> r.accepted,rows)
    @test !any(r -> r.accepted,C.annotate(results[1:2],cases,s))
    changed = copy(results)
    changed[3] = merge(changed[3],(; values=(0.2,0.003)))
    @test !any(r -> r.accepted,filter(r -> r.case.family == :wing_span,C.annotate(changed,cases,s)))
    changed = copy(results)
    changed[3] = merge(changed[3],(; phases=(copy(wave),fill(0.01,12))))
    @test !any(r -> r.accepted,filter(r -> r.case.family == :wing_span,C.annotate(changed,cases,s)))
    selected = C.selected_controls(rows,s)
    @test selected.wing_span == 2
    @test selected.core == 0.02
    @test length(C.confirmation_cases(selected,s)) == 2
    @test_throws ErrorException C.selected_controls(rows,merge(s,(; families=[:core])))
    @test_throws ErrorException C.load_selection(joinpath(ROOT,"missing_selection.toml"))

    # Test handoff with synthetic, hash-verified histories, not asserted flags.
    mktempdir() do directory
        write_fixture(c) = begin
            folder = joinpath(directory,c.label)
            mkpath(folder)
            options = C.aerodynamic_options(s.physical,s.core_mode,c,s)
            steps = C.Integrity.validate_options(options)
            n = steps * s.simulated_revolutions
            dt = c.azimuth_deg / (6s.physical.rpm)
            open(joinpath(folder,"coupled_aerodynamic_history.csv"),"w") do io
                println(io,"time_s,azimuth_deg,wing_CL,propeller_CT,propeller_CQ")
                # Deliberately unsettled CL must not veto a CT/CQ selection.
                writedlm(
                    io,
                    hcat((1:n) .* dt,mod.((1:n) .* c.azimuth_deg,360),collect(1:n),fill(0.02,n),fill(0.003,n)),
                    ','
                )
            end
            for file in C.Integrity.RESULT_FILES[2:3]
                write(joinpath(folder,file),"fixture\n")
            end
            C.Integrity.write_metadata(folder,options)
            merge(C.read_aerodynamic(folder,options,s),(; case=c,reused=false,elapsed_s=0.0))
        end
        evidence = write_fixture.(vcat(cases,C.confirmation_cases(selected,s)))
        path = joinpath(directory,"selection.toml")
        C.publish_selection(path,s,selected,evidence)
        loaded = C.load_selection(path)
        @test loaded.selected == selected
        @test C.data(loaded.physical) == C.data(s.physical)
        @test !haskey(TOML.parsefile(path)["study"],"chang_case")
        original = read(path)
        elastic = C.run_aeroelastic(
            merge(elastic_fixture(),(;
                selection_file=path,output_directory=joinpath(directory,"elastic")))
        )
        @test length(elastic.cases) == 18 && isempty(elastic.results)
        @test elastic.nominal == loaded.selected
        @test elastic.selection == loaded
        @test read(path) == original
        @test isfile(joinpath(elastic.directory,"aerodynamic_input.toml"))
        overridden = E.run_aeroelastic(
            merge(
                elastic_fixture(),
                (;
                    speed_mps=85.0,selection_file=path,output_directory=joinpath(directory,"elastic_U85"))
            )
        )
        @test read(path) == original
        @test overridden.selection.selected == loaded.selected
        @test overridden.selection.physical.speed_mps == 85.0
        @test overridden.selection.physical.rpm ≈ s.physical.rpm * 85 / 65
        @test overridden.aerodynamic_reference == loaded
        recorded = TOML.parsefile(joinpath(overridden.directory,"aerodynamic_input.toml"))
        @test recorded["physical"]["speed_mps"] == 85.0
        @test recorded["aerodynamic_reference"]["physical"]["speed_mps"] == 65.0
        @test !haskey(recorded,"selection_sha256")
        matrix,header = readdlm(joinpath(overridden.directory,"case_matrix.csv"),',',header=true)
        dt_column = findfirst(==("dt_s"),vec(header))
        @test Float64.(matrix[:,dt_column]) ≈
              [c.azimuth_deg / (6overridden.selection.physical.rpm) for c in overridden.cases]
        @test occursin("85.0 m/s",read(joinpath(overridden.directory,"report.md"),String))
        @test TOML.parsefile(joinpath(overridden.directory,"study.toml"))["settings"]["damping_method"] ==
              "moving_block_2"
        custom = E.run_aeroelastic(
            merge(
                elastic_fixture(),
                (;
                    nominal=(wing_chord=3,prop_radial=6,core=nothing),
                    families=[:prop_radial],sweeps=(prop_radial=[2,4,8],),
                    selection_file=path,output_directory=joinpath(directory,"elastic_custom"))
            )
        )
        @test isempty(custom.results) && length(custom.cases) == 3
        @test custom.nominal == merge(loaded.selected,(; wing_chord=3,prop_radial=6))
        @test [c.prop_radial for c in custom.cases] == [2,4,8]
        @test all(c -> c.wing_chord == 3 && c.wing_span == loaded.selected.wing_span,custom.cases)
        @test custom.selection == loaded && custom.aerodynamic_reference == loaded
        @test C.load_selection(path) == loaded && read(path) == original
        custom_matrix,customHeader = readdlm(joinpath(custom.directory,"case_matrix.csv"),',',header=true)
        radial_column = findfirst(==("prop_radial"),vec(customHeader))
        chord_column = findfirst(==("wing_chord"),vec(customHeader))
        @test Int.(custom_matrix[:,radial_column]) == [2,4,8]
        @test Int.(custom_matrix[:,chord_column]) == [3,3,3]
        report = read(joinpath(custom.directory,"report.md"),String)
        @test occursin("Nominal numerical controls (after overrides):",report)
        @test occursin("wing_chord = 3",report)
        nested = E.Convergence.run_aeroelastic(
            merge(
                elastic_fixture(),
                (;
                    speed_mps=85.0,selection_file=path,output_directory=joinpath(directory,"nested_U85"))
            )
        )
        @test nested.selection.physical.speed_mps == 85.0
        payload = TOML.parsefile(path)
        payload["selected"]["core"] *= 2
        C.save_toml(path,payload)
        @test_throws ErrorException C.load_selection(path)
        C.publish_selection(path,s,selected,evidence)
        write(first(evidence).history,"truncated")
        @test_throws ErrorException C.load_selection(path)
    end
    C.load_model()
end

# The nominal baseline and the varied levels are independent, as in stage 1.
@testset "Aeroelastic nominal controls and sweep families" begin
    selected = fixture().nominal
    e = elastic_fixture()
    without_nominal = (; (k => v for (k,v) in pairs(e) if k != :nominal)...)
    @test C.elastic_nominal(selected,without_nominal) == selected
    @test C.elastic_nominal(selected,e) == selected
    @test C.elastic_nominal(selected,merge(e,(; nominal=(;)))) == selected
    inherited = (; (k => nothing for k in C.FIELDS)...)
    @test C.elastic_nominal(selected,merge(e,(; nominal=inherited))) == selected
    @test C.elastic_nominal(selected,merge(e,(; nominal=(wing_span=selected.wing_span,)))) == selected

    overrides = (wing_chord=3,prop_radial=6,azimuth_deg=45.0,core=nothing)
    customized = merge(
        e,
        (; nominal=overrides,families=[:prop_radial,:azimuth_deg],
            sweeps=(prop_radial=[2,4,8],azimuth_deg=[90.0,45.0,30.0]))
    )
    nominal = C.elastic_nominal(selected,customized)
    @test nominal == merge(selected,(; wing_chord=3,prop_radial=6,azimuth_deg=45.0))
    selection = (; physical=fixture().physical,core_mode=:fixed,selected)
    plan = C.validate_elastic(customized,selection)
    @test length(plan) == 6
    @test [c.prop_radial for c in plan if c.family == :prop_radial] == [2,4,8]
    @test [c.azimuth_deg for c in plan if c.family == :azimuth_deg] == [90.0,45.0,30.0]
    # Every unvaried direction keeps the resolved baseline, including overrides.
    for c in plan,k in C.FIELDS
        k == c.family && continue
        @test getproperty(c,k) == getproperty(nominal,k)
    end
    @test selection.selected == selected

    automatic = C.elastic_sweeps(nominal,merge(e,(; families=[:prop_radial,:core,:azimuth_deg])))
    @test automatic.prop_radial == [6,8,9]
    @test automatic.core == [selected.core,selected.core / 2,selected.core / 4]
    @test automatic.azimuth_deg == [45.0,22.5,11.25]
    auto_plan = C.validate_elastic(merge(e,(; nominal=overrides,
            families=[:prop_radial],sweeps=(prop_radial=nothing,))),selection)
    @test [c.prop_radial for c in auto_plan] == [6,8,9]
    @test all(c -> c.wing_chord == 3 && c.azimuth_deg == 45.0,auto_plan)

    for invalid in ((typo=2,),[1,2,3],(wing_span=selected.wing_span + 1,),
        (wing_chord=0,),(prop_radial=2.5,),(azimuth_deg=7.0,))
        @test_throws ErrorException C.elastic_nominal(selected,merge(e,(; nominal=invalid)))
    end
    @test_throws ErrorException C.validate_elastic(merge(customized,(;
            nominal=(wing_span=selected.wing_span + 1,))),selection)
    @test_throws ErrorException C.validate_elastic(merge(e,(;
            families=[:wing_span],sweeps=(wing_span=[2,3,4],))),selection)
    @test_throws ErrorException C.build_cases(nominal,(prop_radial=[2,4],),[:prop_radial]; elastic=true)
    @test_throws ErrorException C.build_cases(nominal,(prop_radial=[2,2,4],),[:prop_radial]; elastic=true)
end

# Execute after lazy module loading, so the explicit model API is in scope.
@testset "Independent aeroelastic configuration" begin
    s = fixture()
    e = elastic_fixture()
    c = C.case(:core,1,s.nominal)
    defaults = C.model_defaults(s.physical,:chord,c,e,mktempdir())
    # Failing default loader proves the adapter does not fall back to chang_case.
    @eval C.Model.ChangAeroelastic chang_case_defaults() = error("Interactive case must not be read")
    config = C.Model.ChangAeroelastic.load_chang_configuration(defaults; env=Dict{String,String}())
    @test config.wing.spanwise_panels == 2
    @test config.wing.symmetric == s.physical.wing_symmetric
    @test isnothing(config.aerodynamic.core_radius_m)
    @test config.aerodynamic.chord_core_factor == 0.01
    @test config.aerodynamic.hub_load_arm_factor == 0.5
    @test config.propeller.rotation_rpm == s.physical.rpm
    @test config.wake.maximum_rows_wing == 4
    @test config.wake.maximum_rows_propeller == 4
    @test config.simulation.freestream_speed_mps == s.physical.speed_mps
    selection = (; physical=s.physical,core_mode=:fixed,selected=s.nominal)
    # Own the no-override fixture; the editable entry may select another speed.
    @test E.operating_selection(selection,merge(e,(; speed_mps=nothing))) == selection
    fast = E.operating_selection(selection,merge(e,(; speed_mps=85.0)))
    fast_defaults = C.model_defaults(fast.physical,:fixed,c,e,mktempdir())
    fast_config = C.Model.ChangAeroelastic.load_chang_configuration(fast_defaults; env=Dict{String,String}())
    @test fast_config.simulation.freestream_speed_mps == 85.0
    @test fast_config.propeller.trim_speed_mps == 85.0
    @test fast_config.propeller.rotation_rpm ≈ s.physical.rpm * 85 / 65
    @test fast.physical.rpm / fast.physical.speed_mps ≈ s.physical.rpm / s.physical.speed_mps
    @test C.fit_start(e,fast.physical) ≈ max(
        0.75,
        e.response.trim_revolutions * 60 / fast.physical.rpm +
        e.response.impulse_duration_s + 60 / fast.physical.rpm
    )
    for speed in (0.0,-1.0,Inf,NaN,true,"85")
        @test_throws ArgumentError E.operating_selection(selection,merge(e,(; speed_mps=speed)))
    end
    plan = C.validate_elastic(e,selection)
    @test length(plan) == 18
    @test all(c -> c.wing_span == 2,plan)
    @test_throws ErrorException C.validate_elastic(
        merge(e,(; families=[:wing_span],sweeps=(wing_span=[2,3,4],))),
        selection
    )
    @test_throws ErrorException C.validate_elastic(
        e,
        merge(selection,(; physical=merge(s.physical,(; wake_shedding_fraction=0.2))))
    )
end

@testset "Explicit aeroelastic inputs without aerodynamic certification" begin
    s = fixture()
    base = E.create_study()
    @test base isa E.ExplicitAeroelasticStudy
    @test isnothing(base.settings.selectionFile)
    @test all(!isnothing,values(base.settings.nominal))
    driver = E.explicit_driver_fingerprint()
    backend = C.SOURCE()
    mktempdir() do directory
        # The incomplete selection is deliberately unusable; explicit mode
        # must neither read it nor change its certification status.
        path = joinpath(directory,"incomplete_selection.toml")
        C.save_toml(path,(; schema=1,status="incomplete"))
        original = read(path)
        settings = merge(
            E.aeroelastic_study(runMode=:full),
            (; selection_file=nothing,
                physical=s.physical,core_mode=s.core_mode,nominal=s.nominal,speed_mps=nothing,
                families=[:wing_span,:prop_radial],
                sweeps=(wing_span=[2,3,4],prop_radial=[1,2,3]),
                output_directory=joinpath(directory,"explicit"),dry_run=true,make_plots=false)
        )
        study = E.create_study(settings)
        @test E.WingPropellerUVLMStudies.solve!(study) === study
        result = study.results
        @test isempty(result.results) && length(result.cases) == 6
        @test result.selection.aerodynamic_verified === false
        @test result.selection.input_mode == "explicit"
        @test isnothing(result.aerodynamic_reference)
        @test result.nominal == s.nominal
        @test [c.wing_span for c in result.cases if c.family == :wing_span] == [2,3,4]
        @test all(c -> c.wing_span == 2,filter(c -> c.family == :prop_radial,result.cases))
        recorded = TOML.parsefile(joinpath(result.directory,"aerodynamic_input.toml"))
        @test recorded["aerodynamic_verified"] === false
        @test recorded["input_mode"] == "explicit"
        @test !haskey(recorded,"selection_sha256")
        saved = TOML.parsefile(joinpath(result.directory,"study.toml"))["settings"]
        @test saved["aerodynamic_verified"] === false
        @test saved["explicit_driver_sha256"] == driver
        @test saved["damping_method"] == "moving_block_2"
        report = read(joinpath(result.directory,"report.md"),String)
        @test occursin("Aerodynamic convergence is **unverified**",report)
        @test occursin("both the aerodynamic and structural meshes",report)
        @test read(path) == original
        @test_throws ErrorException E.run_aeroelastic(merge(settings,(; selection_file=path)))

        # Confirm that the span sweep reaches the production structural model.
        for c in filter(c -> c.family == :wing_span,result.cases)
            defaults = C.model_defaults(s.physical,s.core_mode,c,settings,directory)
            config = C.Model.ChangAeroelastic.load_chang_configuration(defaults; env=Dict{String,String}())
            model = C.Model.ChangAeroelastic.build_chang_model(config)
            @test model.parameters.Ne == c.wing_span
        end
        fast = E.explicit_selection(merge(settings,(; speed_mps=85.0)))
        @test fast.physical.rpm ≈ s.physical.rpm * 85 / 65
        @test fast.aerodynamic_verified === false
        @test_throws ArgumentError E.explicit_selection(merge(settings,(; nominal=nothing)))
        @test_throws ArgumentError E.explicit_selection(merge(settings,(;
            nominal=merge(s.nominal,(; core=nothing)))))
        @test_throws ErrorException E.run_aeroelastic(
            merge(settings,(;
                damping=merge(settings.damping,(; lambda_tolerance_per_s=-1.0))))
        )
        wrapped = E.ExplicitAeroelasticSettings(E.MovingBlockStudy(settings),driver)
        @test_throws ErrorException C.damping_result(
            "unused",
            E.ExplicitAeroelasticSettings(wrapped.movingBlock,"stale"),
            s.physical
        )
    end
    @test C.SOURCE() == backend
end

@testset "Moving-block case diagnostics" begin
    mktempdir() do directory
        path = joinpath(directory,"response_history.csv")
        t = collect(0.0:0.002:6.0)
        wave = exp.(-0.15t) .* sin.(2pi * 4.1t)
        open(path,"w") do io
            println(io,"time_s,propeller_1_pitch_deg,propeller_1_yaw_deg")
            writedlm(io,hcat(t,wave,0.5wave),',')
        end
        base = elastic_fixture()
        e = merge(
            base,
            (; response=merge(base.response,(; end_time_s=last(t))),
                damping=merge(
                    base.damping,
                    (; fit_start_s=0.75,frequency_band_hz=(3.0,6.0),
                        minimum_peaks=8,minimum_fit_r_squared=0.8)
                ))
        )
        wrapped = E.MovingBlockStudy(e)
        result = E.Convergence.damping_result(path,wrapped,fixture().physical)
        @test result.method == "moving_block_2"
        @test result.valid
        @test result.pitch.moving_block_lambda_per_s ≈ -0.15 atol = 0.01
        @test result.yaw.moving_block_lambda_per_s ≈ result.pitch.moving_block_lambda_per_s atol = 1e-10
        @test result.pitch.block_count > 100
        @test isfile(joinpath(directory,"pitch_moving_block.csv"))
        @test isfile(joinpath(directory,"yaw_moving_block.csv"))
        @test result.method_sha256 == wrapped.estimator_sha256
        explicit = E.ExplicitAeroelasticSettings(wrapped,E.explicit_driver_fingerprint())
        @test E.Convergence.damping_result(path,explicit,fixture().physical).pitch == result.pitch
    end
end
end
