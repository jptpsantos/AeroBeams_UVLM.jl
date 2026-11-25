# julia --project=lib/WingPropellerUVLM/studies path/to/aeroelastic_damping_target.jl [--solver]
module AeroelasticDampingTargetTests
using Test,DelimitedFiles,TOML
include(joinpath(@__DIR__,"..","examples","chang_linear_aeroelastic","studies",
    "convergence","chang_convergence_aeroelastic.jl"))
const E = ChangConvergenceAeroelastic
const S = E.WingPropellerUVLMStudies
const C = E.Convergence

function settings(target)
    base = E.aeroelastic_study(; dampingTarget=target)
    return merge(base,(; make_plots=false,
        response=merge(base.response,(; end_time_s=6.0)),
        damping=merge(base.damping,(; fit_start_s=0.75,frequency_band_hz=(3.0,6.0),
            minimum_peaks=8,minimum_fit_r_squared=0.8))))
end

@testset "Damping target settings and backwards compatibility" begin
    wing = E.create_study()
    propeller = E.create_study(; dampingTarget=:propeller)
    @test wing.settings.damping.target == :wing
    @test propeller.settings.damping.target == :propeller
    @test wing.settings.outputDirectory != propeller.settings.outputDirectory
    @test wing.settings.families == propeller.settings.families
    @test_throws ArgumentError E.create_study(; dampingTarget=:unknown)
    @test_throws ArgumentError E.cli_main(["--damping-target=unknown","--dry-run"])
    @test_throws ArgumentError C.Damping.response_channels("wing")
    base = settings(:wing)
    oldDamping = (; (k => v for (k,v) in pairs(base.damping) if k != :target)...)
    old = merge(base,(; damping=oldDamping))
    restored = S.create_AeroelasticConvergenceStudy(old)
    @test restored.settings.damping.target == :propeller
    overridden = S.create_AeroelasticConvergenceStudy(old; damping=(target=:wing,))
    @test overridden.settings.damping.target == :wing
    @test S.legacy_settings(overridden).damping.target == :wing
    @test C.metric_names(true) == ("pitch_lambda_per_s","yaw_lambda_per_s",
        "pitch_frequency_hz","yaw_frequency_hz")
    @test C.metric_names(true; damping_target=:wing) == ("wing_bending_lambda_per_s",
        "wing_torsion_lambda_per_s","wing_bending_frequency_hz","wing_torsion_frequency_hz")
end

@testset "Wing and propeller fits select different physical signals" begin
    mktempdir() do directory
        t = collect(0.0:0.002:6.0)
        lambdas = (-0.15,-0.25,-0.35,-0.45)
        scales = (0.2,0.1,20.0,0.5)
        signals = ntuple(i -> scales[i] .* exp.(lambdas[i] .* t) .* sin.(2pi * 4.1 .* t),4)
        path = joinpath(directory,"response_history.csv")
        open(path,"w") do io
            println(io,"time_s,propeller_1_pitch_deg,propeller_1_yaw_deg,tip_displacement_m,tip_twist_deg")
            writedlm(io,hcat(t,signals...),',')
        end
        for (target,indices) in ((:propeller,(1,2)),(:wing,(3,4)))
            s = settings(target)
            wrapped = E.MovingBlockStudy(s)
            result = C.damping_result(path,wrapped,s.physical)
            channels = C.Damping.response_channels(target)
            response = C.Damping.read_response_history(path; target)
            @test result.target == target
            @test result.valid
            @test response.time == t
            for (i,name) in enumerate(channels.names)
                @test getproperty(response,name) == signals[indices[i]]
                fit = getproperty(result,name)
                @test isapprox(fit.moving_block_lambda_per_s,lambdas[indices[i]]; atol=0.01)
                @test getproperty(result.signals,name).column == channels.columns[i]
                @test getproperty(result.signals,name).unit == channels.units[i]
                @test isfile(joinpath(directory,"$(channels.prefixes[i])_moving_block.csv"))
            end
            # The explicit entry delegates to the same selected estimator.
            explicit = E.ExplicitAeroelasticSettings(wrapped,E.explicit_driver_fingerprint())
            fromEntry = C.damping_result(path,explicit,s.physical)
            @test getproperty(fromEntry,first(channels.names)) == getproperty(result,first(channels.names))
            preview = E.ExplicitAeroelasticSettings(E.MovingBlockStudy(merge(s,(; run_mode=:preview))),
                E.explicit_driver_fingerprint())
            @test !C.damping_result(path,preview,s.physical).valid
            # The legacy backend also selects the requested component.
            legacy = C.damping_result(path,s,s.physical)
            @test legacy.target == target
            @test isapprox(getproperty(legacy,first(channels.names)).moving_block_lambda_per_s,
                lambdas[first(indices)]; atol=0.02)
            C.save_toml(joinpath(directory,"$(target)_fit.toml"),result)
            saved = TOML.parsefile(joinpath(directory,"$(target)_fit.toml"))
            @test saved["target"] == string(target)
            @test saved["signals"][string(first(channels.names))]["unit"] == first(channels.units)
            @test C.data(wrapped)["damping_target"] == string(target)
        end
        wing = settings(:wing)
        limited = merge(wing,(; integration=merge(wing.integration,(; propeller_angle_limit_deg=0.3))))
        response = C.Damping.read_response_history(path; target=:wing)
        # A 20-m synthetic bending signal must not be compared with an angular limit.
        @test C.damping_amplitude_ok(path,limited,response)
        @test C.damping_result(path,E.MovingBlockStudy(limited),limited.physical).valid
        unsafe = merge(limited,(; integration=merge(limited.integration,(; propeller_angle_limit_deg=0.01))))
        @test !C.damping_result(path,E.MovingBlockStudy(unsafe),unsafe.physical).valid
        wingOnly = joinpath(directory,"wing_only.csv")
        open(wingOnly,"w") do io
            println(io,"time_s,tip_displacement_m,tip_twist_deg")
            writedlm(io,hcat(t,signals[3],signals[4]),',')
        end
        @test C.damping_result(wingOnly,E.MovingBlockStudy(wing),wing.physical).valid
        @test_throws ErrorException C.Damping.read_response_history(wingOnly)
        @test_throws ErrorException C.damping_result(wingOnly,E.MovingBlockStudy(limited),limited.physical)
    end
end

@testset "Selected target reaches checkpoint CSVs" begin
    mktempdir() do directory
        base = settings(:wing)
        nominal = merge(base.nominal,(; wing_span=2,wing_chord=1,prop_radial=2,prop_chord=1))
        cases = C.build_cases(nominal,(; wing_span=[2,3,4]),[:wing_span])
        for target in (:wing,:propeller)
            s = merge(base,(; damping=merge(base.damping,(; target)),dry_run=false,
                output_directory=joinpath(directory,string(target))))
            runner(c,folder) = (; status="completed",reason="Synthetic fit",values=(-0.2,-0.3,4.0,4.1),
                phases=(),periodic=(),valid=true,history="fixture.csv",log="",reused=false)
            rows = C.execute_cases(cases,E.MovingBlockStudy(s),s.physical,s.core_mode; elastic=true,runner)
            @test length(rows) == 3
            _,header = readdlm(joinpath(s.output_directory,"convergence.csv"),','; header=true)
            names = string.(vec(header))
            metrics = C.metric_names(true; damping_target=target)
            @test all(name -> name in names,metrics)
            @test all(name -> "previous_change_$name" in names,metrics)
            @test all(name -> "limit_$name" in names,metrics)
            @test ("pitch_lambda_per_s" in names) == (target == :propeller)
            @test ("wing_bending_lambda_per_s" in names) == (target == :wing)
            saved = TOML.parsefile(joinpath(s.output_directory,"study.toml"))["settings"]
            @test saved["damping"]["target"] == string(target)
            @test saved["damping_target"] == string(target)
        end
    end
end

if "--solver" in ARGS
    @testset "Coupled solver wing-response smoke case" begin
        mktempdir() do directory
            base = E.aeroelastic_study(; runMode=:preview,dampingTarget=:wing)
            s = merge(base,(; output_directory=directory,make_plots=false))
            selection = E.explicit_selection(s)
            C.load_model()
            wrapped = E.ExplicitAeroelasticSettings(E.MovingBlockStudy(s),E.explicit_driver_fingerprint())
            cases = Base.invokelatest(C.validate_elastic,wrapped,selection)
            result = Base.invokelatest(C.aeroelastic_case,first(cases),wrapped,selection,directory)
            @test result.status == "indeterminate" # Preview is too short to certify damping.
            @test !result.valid
            saved = TOML.parsefile(joinpath(directory,"damping_fit.toml"))
            @test saved["target"] == "wing"
            @test haskey(saved,"bending") && haskey(saved,"torsion")
            @test !haskey(saved,"pitch") && !haskey(saved,"yaw")
            @test saved["signals"]["bending"]["unit"] == "m"
            response = C.Damping.read_response_history(result.history; target=:wing)
            @test C.Damping.completed_requested_window(response.time,s.response.end_time_s)
            @test maximum(abs,response.bending) > 0
            @test maximum(abs,response.torsion) > 0
            @test isfile(joinpath(directory,"wing_bending_moving_block.csv"))
            @test isfile(joinpath(directory,"wing_torsion_moving_block.csv"))
        end
    end
end
end
