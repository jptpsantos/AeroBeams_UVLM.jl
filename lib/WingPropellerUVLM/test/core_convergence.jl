module CoreConvergenceTests
using Test, DelimitedFiles
include(joinpath(@__DIR__, "..", "examples", "chang_linear_aeroelastic",
    "studies", "convergence", "archive", "run_chang_core_convergence.jl"))
const C = ChangCoreConvergence
const I = C.Integrity

@testset "Fixed-core convergence study" begin
    config = C.CaseModel.ChangAeroelastic.load_chang_configuration(; env=Dict{String,String}())
    # Study presets are user-editable; numerical tests own their level lists.
    s = merge(C.core_convergence_defaults(config),(;
        mode=:aerodynamic,radii_m=[.01,.003,.001],nominal_radius_m=.001,
        nominal_mesh=(wing_span=20,wing_chord=5,prop_radial=5,prop_chord=5),
        panel_levels=(wing_span=[20,40,60],wing_chord=[5,10,15],
            prop_radial=[5,10,15],prop_chord=[5,10,15]),
        core_meshes=[(wing_span=20,wing_chord=5,prop_radial=5,prop_chord=5),
            (wing_span=40,wing_chord=10,prop_radial=10,prop_chord=10),
            (wing_span=60,wing_chord=15,prop_radial=15,prop_chord=15)],
        families=[:core_mesh,:wing_span,:wing_chord,:prop_radial,:prop_chord,:time_step,:wake_length],
        azimuth_levels_deg=[5.,2.5,1.25],wake_levels_revolutions=[1.,2.,3.],
        nominal_wake_revolutions=2.,simulated_revolutions=8,averaged_revolutions=2,
        coefficient_absolute_tolerances=(1e-4,1e-6,1e-6),coefficient_relative_tolerance=.01,
        periodic_tolerances=(1e-4,1e-6,1e-6),lambda_tolerance_per_s=.01,frequency_tolerance_hz=.1))
    cases = C.build_cases(s, config)
    @test s.nominal_radius_m in s.radii_m
    @test length(filter(c -> startswith(string(c.family), "core_mesh_"), cases)) ==
        length(s.radii_m) * length(s.core_meshes)
    @test length(unique(C.case_key.(cases))) < length(cases)
    for (i,radius) in enumerate(s.radii_m)
        mesh_cases = filter(c -> c.family == Symbol("mesh_at_core_$i"),cases)
        @test length(mesh_cases) == length(s.core_meshes)
        @test all(c -> c.radius_m == radius,mesh_cases)
        @test [c.wing_span for c in mesh_cases] == [m.wing_span for m in s.core_meshes]
    end
    c = first(cases)
    @test C.case_key(c) != C.case_key(merge(c, (; radius_m=c.radius_m/2)))
    options = C.rigid_options(c,s,config)
    @test options.core_radius_m == c.radius_m
    @test options.flow_speed_mps == config.simulation.freestream_speed_mps
    @test options.interaction_on == config.simulation.interaction_on
    @test options.wing_symmetric == config.wing.symmetric
    opposite = merge(config,(;wing=merge(config.wing,(;symmetric=!config.wing.symmetric))))
    @test C.rigid_options(c,s,opposite).wing_symmetric == !options.wing_symmetric
    @test I.options_from_environment(Dict("CHANG_AERO_WING_SYMMETRIC"=>"true")).wing_symmetric
    @test !I.options_from_environment(Dict("CHANG_AERO_WING_SYMMETRIC"=>"false")).wing_symmetric
    @test options.propeller_attachment_eta == only(config.propeller.attachment_eta)
    @test options.collective_pitch_offset_deg == config.propeller.collective_pitch_offset_deg
    timecases = filter(c -> c.family == :time_step, cases)
    @test all(c -> isapprox(C.time_step(c,config)*C.wake_rows(c),
        s.nominal_wake_revolutions*60/C.rotor_rpm(config)), timecases)
    @test all(c -> c.radius_m == s.nominal_radius_m, timecases)
    ae = merge(s, (; mode=:aeroelastic))
    @test all(c -> c.wing_span == s.nominal_mesh.wing_span, C.build_cases(ae,config))
    for change in ((; radii_m=[.01,.01,.001]), (; radii_m=[.01,0.,-.001]),
        (; azimuth_levels_deg=[7.,2.5,1.25]), (; simulated_revolutions=2),
        (; panel_levels=merge(s.panel_levels,(; wing_span=[20,30.5,40]))),
        (; panel_levels=merge(s.panel_levels,(; wing_chord=[5,5,10]))),
        (; nominal_mesh=merge(s.nominal_mesh,(; prop_radial=0))),
        (; core_meshes=reverse(s.core_meshes)), (; families=[:unknown]))
        @test_throws ErrorException C.build_cases(merge(s,change),config)
    end

    # Nonproportional counts and unequal family lengths must remain exact.
    explicit = merge(s,(;
        nominal_mesh=(wing_span=24,wing_chord=6,prop_radial=7,prop_chord=3),
        panel_levels=(wing_span=[18,27,43,55],wing_chord=[3,7,11],
            prop_radial=[6,9,14],prop_chord=[2,4,7]),
        core_meshes=[(wing_span=18,wing_chord=3,prop_radial=6,prop_chord=2),
            (wing_span=27,wing_chord=7,prop_radial=9,prop_chord=4),
            (wing_span=43,wing_chord=11,prop_radial=14,prop_chord=7)]))
    explicit_cases = C.build_cases(explicit,config)
    for direction in C.MESH_FIELDS
        group = filter(c -> c.family == direction,explicit_cases)
        @test [getproperty(c,direction) for c in group] == getproperty(explicit.panel_levels,direction)
        @test all(c -> all(k -> k == direction || getproperty(c,k) ==
            getproperty(explicit.nominal_mesh,k),C.MESH_FIELDS),group)
    end
    @test all(c -> C.mesh_controls(c) == explicit.nominal_mesh,
        filter(c -> c.family in (:time_step,:wake_length),explicit_cases))
    for (i,mesh) in enumerate(explicit.core_meshes)
        @test all(c -> C.mesh_controls(c) == mesh,
            filter(c -> c.family == Symbol("core_mesh_$i"),explicit_cases))
    end
    changed_case = merge(config,(; wing=merge(config.wing,(;spanwise_panels=99,chordwise_panels=17)),
        propeller=merge(config.propeller,(;radial_panels=23,chordwise_panels=19))))
    @test C.case_key.(C.build_cases(explicit,changed_case)) == C.case_key.(explicit_cases)
    @test all(c -> c.wing_span == 24,C.build_cases(merge(explicit,(;mode=:aeroelastic)),changed_case))
    span_only = [merge(explicit.nominal_mesh,(;wing_span=n)) for n in (18,27,43)]
    @test_throws ErrorException C.build_cases(merge(explicit,(;mode=:aeroelastic,core_meshes=span_only)),config)

    # A finest case cannot pass by agreeing with itself. Missing cases, an
    # unsettled final pair, or a changed waveform must prevent acceptance.
    group = [merge(c,(; family=:core_mesh_1,level=i,label="test_$i")) for i in 1:3]
    record(c) = (; case=c,status="completed",reason="verified, \"quoted\"\nrow",reused=false,
        valid=true,values=(.3,.001,.0001),phases=(fill(.3,12),fill(.001,12),fill(.0001,12)),
        periodic=(0.,0.,0.),elapsed_s=1.,history="history,with,commas.csv",log="run.log")
    rows = record.(group)
    annotate(r) = C.compare_results(r,group,s)
    @test all(r -> r.accepted, annotate(rows))
    @test !any(r -> r.accepted, annotate(rows[1:2]))
    for change in ((; valid=false), (; status="failed",valid=false),
        (; values=(.32,.001,.0001)),
        (; phases=(.3 .+ .01sin.(2pi*(1:12)/12),fill(.001,12),fill(.0001,12))))
        altered = copy(rows)
        altered[end] = merge(rows[end],change)
        @test !any(r -> r.accepted,annotate(altered))
    end
    altered = copy(rows)
    altered[1] = merge(rows[1],(; values=(.32,.001,.0001)))
    @test [r.accepted for r in annotate(altered)] == [false,true,true]
    elastic_rows = [merge(r,(; values=(-.1,-.12,4.,4.2),phases=(),periodic=())) for r in rows]
    @test all(r -> r.accepted,C.compare_results(elastic_rows,group,ae))
    elastic_rows[end] = merge(elastic_rows[end],(; values=(-.1,-.12,4.3,4.2)))
    @test !any(r -> r.accepted,C.compare_results(elastic_rows,group,ae))

    mktempdir() do directory
        C.write_tables(directory,annotate(rows),group,s,config)
        raw,header = readdlm(joinpath(directory,"convergence.csv"),',',header=true)
        @test size(raw,2) == length(header)
        @test raw[1,6] == rows[1].reason
        C.write_fit_diagnostics(directory,(; status="indeterminate",reason="fixture",
            pitch=(; fit_r_squared=NaN,valid_for_convergence=false),
            yaw=(; fit_r_squared=.99,valid_for_convergence=true)))
        @test isfile(joinpath(directory,"damping_fit.toml"))
        for file in I.RESULT_FILES
            write(joinpath(directory,file),"synthetic metadata fixture\n")
        end
        I.write_metadata(directory,options)
        @test I.metadata_options(directory).core_radius_m == c.radius_m
        @test I.output_matches_options(directory,options;warn_on_mismatch=false)
        @test !I.output_matches_options(directory,C.rigid_options(c,s,opposite);warn_on_mismatch=false)
        different = C.rigid_options(merge(c,(; radius_m=c.radius_m/2)),s,config)
        @test !I.output_matches_options(directory,different;warn_on_mismatch=false)
        write(joinpath(directory,first(I.RESULT_FILES)),"truncated")
        @test !I.output_matches_options(directory,options;warn_on_mismatch=false)
    end
end
end
