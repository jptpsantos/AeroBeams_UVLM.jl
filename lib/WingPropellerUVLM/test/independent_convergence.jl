module IndependentConvergenceTests
using Test, DelimitedFiles, TOML
const ROOT=joinpath(@__DIR__,"..","examples","chang_linear_aeroelastic","studies","convergence")
# Load definitions even if the interactive entry's execution guard is disabled.
Base.include(ex -> ex == :(ChangConvergenceAerodynamic.main()) ? nothing : ex,
    @__MODULE__,joinpath(ROOT,"chang_convergence_aerodynamic.jl"))
include(joinpath(ROOT,"chang_convergence_aeroelastic.jl"))
const C=ChangConvergenceAerodynamic.Convergence
const E=ChangConvergenceAeroelastic

function fixture()
    s=ChangConvergenceAerodynamic.aerodynamic_study()
    return merge(s,(;dry_run=true,make_plots=false,
        nominal=(wing_span=2,wing_chord=1,prop_radial=2,prop_chord=1,wake_revolutions=1.,core=.01,azimuth_deg=90.),
        sweeps=(wing_span=[2,3,4],wing_chord=[1,2,3],prop_radial=[2,3,4],prop_chord=[1,2,3],
            wake_revolutions=[1.,2.,3.],core=[.02,.01,.005],azimuth_deg=[90.,45.,30.]),
        simulated_revolutions=6,averaged_revolutions=1))
end

@testset "Independent two-stage convergence" begin
    s=fixture(); cases=C.validate_aerodynamic(s)
    @test length(cases)==21
    @test length(unique(C.key.(cases)))==15
    @test C.metric_names(false)==("CT","CQ")
    @test C.core_factors(:fixed,.01)==(radius=.01,chord=0.,segment=0.)
    @test C.core_factors(:chord,.01)==(radius=nothing,chord=.01,segment=0.)
    @test C.core_factors(:segment,.01)==(radius=nothing,chord=0.,segment=.01)
    @test_throws ErrorException C.core_factors(:unknown,.1)
    @test_throws ErrorException C.build_cases(s.nominal,merge(s.sweeps,(;azimuth_deg=[90.,45.])),s.families)
    @test_throws ErrorException C.build_cases(s.nominal,merge(s.sweeps,(;wing_span=[2,2,4])),s.families)
    wave=zeros(12)
    record(c)= (;case=c,status="completed",reason="fixture",reused=false,elapsed_s=0.,
        values=(.02,.003),phases=(copy(wave),copy(wave)),periodic=(0.,0.),valid=true,history="",log="")
    results=record.(cases)
    rows=C.annotate(results,cases,s)
    @test all(r->r.accepted,rows)
    @test !any(r->r.accepted,C.annotate(results[1:2],cases,s))
    changed=copy(results); changed[3]=merge(changed[3],(;values=(.2,.003)))
    @test !any(r->r.accepted,filter(r->r.case.family==:wing_span,C.annotate(changed,cases,s)))
    changed=copy(results); changed[3]=merge(changed[3],(;phases=(copy(wave),fill(.01,12))))
    @test !any(r->r.accepted,filter(r->r.case.family==:wing_span,C.annotate(changed,cases,s)))
    selected=C.selected_controls(rows,s)
    @test selected.wing_span==2
    @test selected.core==.02
    @test length(C.confirmation_cases(selected,s))==2
    @test_throws ErrorException C.selected_controls(rows,merge(s,(;families=[:core])))
    @test_throws ErrorException C.load_selection(joinpath(ROOT,"missing_selection.toml"))

    # Test handoff with synthetic, hash-verified histories, not asserted flags.
    mktempdir() do directory
        write_fixture(c)=begin
            folder=joinpath(directory,c.label); mkpath(folder)
            options=C.aerodynamic_options(s.physical,s.core_mode,c,s)
            steps=C.Integrity.validate_options(options); n=steps*s.simulated_revolutions
            dt=c.azimuth_deg/(6s.physical.rpm)
            open(joinpath(folder,"coupled_aerodynamic_history.csv"),"w") do io
                println(io,"time_s,azimuth_deg,wing_CL,propeller_CT,propeller_CQ")
                # Deliberately unsettled CL must not veto a CT/CQ selection.
                writedlm(io,hcat((1:n).*dt,mod.((1:n).*c.azimuth_deg,360),collect(1:n),fill(.02,n),fill(.003,n)),',')
            end
            for file in C.Integrity.RESULT_FILES[2:3]; write(joinpath(folder,file),"fixture\n"); end
            C.Integrity.write_metadata(folder,options)
            merge(C.read_aerodynamic(folder,options,s),(;case=c,reused=false,elapsed_s=0.))
        end
        evidence=write_fixture.(vcat(cases,C.confirmation_cases(selected,s)))
        path=joinpath(directory,"selection.toml")
        C.publish_selection(path,s,selected,evidence)
        loaded=C.load_selection(path)
        @test loaded.selected==selected
        @test C.data(loaded.physical)==C.data(s.physical)
        @test !haskey(TOML.parsefile(path)["study"],"chang_case")
        elastic=C.run_aeroelastic(merge(E.aeroelastic_study(),(;dry_run=true,make_plots=false,
            selection_file=path,output_directory=joinpath(directory,"elastic"))))
        @test length(elastic.cases)==18 && isempty(elastic.results)
        @test isfile(joinpath(elastic.directory,"aerodynamic_input.toml"))
        original=read(path)
        overridden=E.run_aeroelastic(merge(E.aeroelastic_study(),(;dry_run=true,make_plots=false,
            speed_mps=85.0,selection_file=path,output_directory=joinpath(directory,"elastic_U85"))))
        @test read(path)==original
        @test overridden.selection.selected==loaded.selected
        @test overridden.selection.physical.speed_mps==85.0
        @test overridden.selection.physical.rpm≈s.physical.rpm*85/65
        @test overridden.aerodynamic_reference==loaded
        recorded=TOML.parsefile(joinpath(overridden.directory,"aerodynamic_input.toml"))
        @test recorded["physical"]["speed_mps"]==85.0
        @test recorded["aerodynamic_reference"]["physical"]["speed_mps"]==65.0
        @test !haskey(recorded,"selection_sha256")
        matrix,header=readdlm(joinpath(overridden.directory,"case_matrix.csv"),',',header=true)
        dt_column=findfirst(==("dt_s"),vec(header))
        @test Float64.(matrix[:,dt_column])≈[c.azimuth_deg/(6overridden.selection.physical.rpm) for c in overridden.cases]
        @test occursin("85.0 m/s",read(joinpath(overridden.directory,"report.md"),String))
        @test TOML.parsefile(joinpath(overridden.directory,"study.toml"))["settings"]["damping_method"]=="moving_block_2"
        nested=E.Convergence.run_aeroelastic(merge(E.aeroelastic_study(),(;dry_run=true,make_plots=false,
            speed_mps=85.0,selection_file=path,output_directory=joinpath(directory,"nested_U85"))))
        @test nested.selection.physical.speed_mps==85.0
        payload=TOML.parsefile(path); payload["selected"]["core"]*=2
        C.save_toml(path,payload)
        @test_throws ErrorException C.load_selection(path)
        C.publish_selection(path,s,selected,evidence)
        write(first(evidence).history,"truncated")
        @test_throws ErrorException C.load_selection(path)
    end
    C.load_model()
end

# Execute after lazy module loading, so the explicit model API is in scope.
@testset "Independent aeroelastic configuration" begin
    s=fixture(); e=E.aeroelastic_study(); c=C.case(:core,1,s.nominal)
    defaults=C.model_defaults(s.physical,:chord,c,e,mktempdir())
    # Failing default loader proves the adapter does not fall back to chang_case.
    @eval C.Model.ChangAeroelastic chang_case_defaults()=error("Interactive case must not be read")
    config=C.Model.ChangAeroelastic.load_chang_configuration(defaults;env=Dict{String,String}())
    @test config.wing.spanwise_panels==2
    @test config.wing.symmetric==s.physical.wing_symmetric
    @test isnothing(config.aerodynamic.core_radius_m)
    @test config.aerodynamic.chord_core_factor==.01
    @test config.propeller.rotation_rpm==s.physical.rpm
    @test config.wake.maximum_rows_wing==4
    @test config.wake.maximum_rows_propeller==4
    @test config.simulation.freestream_speed_mps==s.physical.speed_mps
    selection=(;physical=s.physical,core_mode=:fixed,selected=s.nominal)
    @test E.operating_selection(selection,e)==selection
    fast=E.operating_selection(selection,merge(e,(;speed_mps=85.0)))
    fast_defaults=C.model_defaults(fast.physical,:fixed,c,e,mktempdir())
    fast_config=C.Model.ChangAeroelastic.load_chang_configuration(fast_defaults;env=Dict{String,String}())
    @test fast_config.simulation.freestream_speed_mps==85.0
    @test fast_config.propeller.trim_speed_mps==85.0
    @test fast_config.propeller.rotation_rpm≈s.physical.rpm*85/65
    @test fast.physical.rpm/fast.physical.speed_mps≈s.physical.rpm/s.physical.speed_mps
    @test C.fit_start(e,fast.physical)≈max(.75,e.response.trim_revolutions*60/fast.physical.rpm+
        e.response.impulse_duration_s+60/fast.physical.rpm)
    for speed in (0.0,-1.0,Inf,NaN,true,"85")
        @test_throws ArgumentError E.operating_selection(selection,merge(e,(;speed_mps=speed)))
    end
    plan=C.validate_elastic(e,selection)
    @test length(plan)==18
    @test all(c->c.wing_span==2,plan)
    @test_throws ErrorException C.validate_elastic(merge(e,(;families=[:wing_span],sweeps=(wing_span=[2,3,4],))),selection)
    @test_throws ErrorException C.validate_elastic(e,merge(selection,(;physical=merge(s.physical,(;wake_shedding_fraction=.2)))))
end

@testset "Moving-block case diagnostics" begin
    mktempdir() do directory
        path=joinpath(directory,"response_history.csv")
        t=collect(0.:.002:6.)
        wave=exp.(-.15t).*sin.(2pi*4.1t)
        open(path,"w") do io
            println(io,"time_s,propeller_1_pitch_deg,propeller_1_yaw_deg")
            writedlm(io,hcat(t,wave,.5wave),',')
        end
        e=merge(E.aeroelastic_study(),(;make_plots=false))
        wrapped=E.MovingBlockStudy(e)
        result=E.Convergence.damping_result(path,wrapped,fixture().physical)
        @test result.method=="moving_block_2"
        @test result.valid
        @test result.pitch.moving_block_lambda_per_s≈-.15 atol=.01
        @test result.yaw.moving_block_lambda_per_s≈result.pitch.moving_block_lambda_per_s atol=1e-10
        @test result.pitch.block_count>100
        @test isfile(joinpath(directory,"pitch_moving_block.csv"))
        @test isfile(joinpath(directory,"yaw_moving_block.csv"))
        @test result.method_sha256==wrapped.estimator_sha256
    end
end
end
