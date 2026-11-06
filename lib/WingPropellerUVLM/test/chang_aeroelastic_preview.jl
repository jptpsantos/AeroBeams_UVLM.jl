module ChangAeroelasticPreviewTests
using Test, TOML, DelimitedFiles

const ENTRY=joinpath(@__DIR__,"..","examples","chang_linear_aeroelastic",
    "studies","convergence","chang_convergence_aeroelastic.jl")
# Including the interactive entry must only define the study API.
include(ENTRY)
const E=ChangConvergenceAeroelastic
const C=E.Convergence

@testset "Chang aeroelastic preview planning" begin
    preview=E.create_study(runMode=:preview)
    full=E.create_study()
    @test preview.results===nothing && full.results===nothing
    @test preview isa E.ExplicitAeroelasticStudy
    @test preview.settings.runMode==:preview
    @test preview.settings.response.finalTime==0.15
    @test !preview.settings.makePlots
    @test preview.settings.coupling==full.settings.coupling
    @test preview.settings.damping.minimumPeaks==full.settings.damping.minimumPeaks
    @test preview.settings.outputDirectory==full.settings.outputDirectory*"_preview"
    @test full.settings.runMode==:full
    @test full.settings.coupling.scheme==:loose_explicit
    @test full.settings.response.finalTime==4.0
    @test_throws ArgumentError E.create_study(runMode=:unknown)
    @test_throws ArgumentError E.cli_main(["--unknown"])

    s=E.aeroelastic_study(runMode=:preview)
    @test s.run_mode==:preview
    pulse_end=s.response.trim_revolutions*60/s.physical.rpm+s.response.impulse_duration_s
    @test 0<pulse_end<s.damping.fit_start_s<s.response.end_time_s
    @test_throws ArgumentError E.create_study(s;runMode=:unknown)
    @test_throws ArgumentError E.create_study(s;selectionFile="unused_selection.toml")
    @test_throws ArgumentError E.run_aeroelastic(merge(s,(;selection_file="unused_selection.toml")))
    plan=E.runtime_plan(s,E.explicit_selection(s))
    @test length(plan)==3
    @test [r.family for r in plan]==fill(:prop_radial,3)
    @test [r.level for r in plan]==[1,2,3]
    @test all(r->!r.reused,plan)
    @test all(r->r.steps==110,plan)
    @test all(r->isapprox(r.dt_s,10/(6s.physical.rpm)),plan)
    @test [r.bound_panels for r in plan]==[16,20,24]
    @test [r.wake_rows for r in plan]==[18,18,18]
    @test [r.wake_panels_capacity for r in plan]==[216,288,360]
    sf=E.aeroelastic_study(runMode=:full)
    fullplan=E.runtime_plan(sf,E.explicit_selection(sf))
    @test length(fullplan)==15
    @test count(r->!r.reused,fullplan)==11
    timeplan=filter(r->r.family==:azimuth_deg,fullplan)
    @test length(timeplan)==3
    @test all(<(0),diff([r.dt_s for r in timeplan]))

    mktempdir() do directory
        full_directory=joinpath(directory,"convergence_aeroelastic")
        mkpath(full_directory)
        sentinel=joinpath(full_directory,"existing_result.txt")
        write(sentinel,"preserve existing full study\n")
        preview_directory=full_directory*"_preview"
        dry=E.create_study(s;dryRun=true,outputDirectory=preview_directory)
        @test E.WingPropellerUVLMStudies.solve!(dry)===dry
        result=dry.results
        @test result.directory==abspath(preview_directory)
        @test length(result.cases)==3
        @test isempty(result.results)
        @test all(c->!isdir(joinpath(result.directory,c.label)),result.cases)
        @test read(sentinel,String)=="preserve existing full study\n"
        @test readdir(full_directory)==["existing_result.txt"]
        @test isfile(joinpath(result.directory,"runtime_plan.csv"))
        raw,header=readdlm(joinpath(result.directory,"runtime_plan.csv"),',',header=true)
        @test size(raw,1)==3
        steps_column=findfirst(==("steps"),vec(header))
        panels_column=findfirst(==("bound_panels"),vec(header))
        @test Int.(raw[:,steps_column])==[110,110,110]
        @test Int.(raw[:,panels_column])==[16,20,24]
        rows,headers=readdlm(joinpath(result.directory,"convergence.csv"),',',header=true)
        @test "accepted" in vec(headers)
        @test isempty(rows)
        saved=TOML.parsefile(joinpath(result.directory,"study.toml"))["settings"]
        @test saved["run_mode"]=="preview"
        @test saved["dry_run"]===true
        @test occursin("Preview only",read(joinpath(result.directory,"report.md"),String))
    end
end

@testset "Preview cannot certify damping even with a long signal" begin
    mktempdir() do directory
        history=joinpath(directory,"response_history.csv")
        t=collect(0.:.002:6.)
        wave=exp.(-.15t).*sin.(2pi*4.1t)
        open(history,"w") do io
            println(io,"time_s,propeller_1_pitch_deg,propeller_1_yaw_deg")
            writedlm(io,hcat(t,wave,.5wave),',')
        end
        base=E.aeroelastic_study(runMode=:full)
        settings=merge(base,(;response=merge(base.response,(;end_time_s=last(t))),
            damping=merge(base.damping,(;fit_start_s=.75))))
        explicit(s)=E.ExplicitAeroelasticSettings(E.MovingBlockStudy(s),E.explicit_driver_fingerprint())
        full=C.damping_result(history,explicit(settings),settings.physical)
        preview=C.damping_result(history,explicit(merge(settings,(;run_mode=:preview))),settings.physical)
        @test full.valid
        @test full.pitch.moving_block_lambda_per_s≈-.15 atol=.01
        @test !preview.valid
        @test preview.analysis_mode=="preview"
        @test preview.pitch==full.pitch
        @test preview.yaw==full.yaw
    end
end
end
