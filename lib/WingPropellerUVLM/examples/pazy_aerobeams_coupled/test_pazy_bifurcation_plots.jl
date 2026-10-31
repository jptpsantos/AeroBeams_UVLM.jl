# Synthetic signals validate extraction and presentation, not physical branches.
using Test, DelimitedFiles
include(joinpath(@__DIR__, "..", "bifurcation", "TimeDomainBifurcation.jl"))
include(joinpath(@__DIR__, "PazyBifurcationPlots.jl"))
using .TimeDomainBifurcation, .PazyBifurcationPlots

@testset "Local extrema preserve final-window samples" begin
    t = collect(0.0:8.0)
    q = [8., 0., 1., 3., 3., 1., 0., 2., 4.]
    points = response_local_extrema(t, q; transient_time=2.)
    @test points == [(kind="maximum", time_s=4., value=3.),
                     (kind="minimum", time_s=6., value=0.)]
    @test isempty(response_local_extrema(t, fill(1.,9); transient_time=2.))
    @test isempty(response_local_extrema(t, t; transient_time=2.))
    @test_throws ArgumentError response_local_extrema(t, q; transient_time=8.)
    @test_throws ArgumentError response_local_extrema(t, q[1:end-1]; transient_time=2.)
end

@testset "CSV-only paper-style Pazy figures" begin
    directory = "--keep-output" in ARGS ?
        joinpath(@__DIR__, "..", "..", "..", "..", "tmp", "pazy_bifurcation_style_validation_" * string(time_ns())) : mktempdir()
    definitions = (
        tip_out_of_plane_percent=(unit="% semispan", label="Tip OOP displacement", atol=1e-4),
        root_strain_ue=(unit="microstrain", label="Near-root longitudinal strain", atol=0.1),
        tip_twist_deg=(unit="deg", label="Tip twist", atol=1e-3),
        tip_geometric_incidence_deg=(unit="deg", label="Tip geometric incidence", atol=1e-3))
    function simulate(speed)
        t = collect(0:0.002:2)
        waveform = sin.(16pi .* t) .+ 0.25sin.(26pi .* t)
        amplitude = max(0., 1 - ((speed-45)/5)^2)
        static = (; tip_out_of_plane_percent=21+0.6(speed-40),
            root_strain_ue=-1800-65(speed-40), tip_twist_deg=1+0.04(speed-40),
            tip_geometric_incidence_deg=6+0.04(speed-40))
        scales = (3., 600., 5., 5.)
        responses = NamedTuple{keys(static)}(Tuple(value .+ scale*amplitude .* waveform
            for (value,scale) in zip(values(static),scales)))
        (; time=t, responses, static_equilibrium=static)
    end
    try
        sweep = run_bifurcation(simulate; U_min=40., U_max=50., dU=1., definitions,
            simulation_time=2., transient_time=1., plot_results=false,
            metadata=(angle_of_attack=deg2rad(5.), note="SYNTHETIC plotting validation only"),
            output_directory=joinpath(directory,"synthetic"))
        @test isfile(sweep.local_extrema_path)
        peaks, header = readdlm(sweep.local_extrema_path, ','; header=true)
        @test size(peaks,2) == 6
        @test size(peaks,1) > 11*4*2 # More than a min/max envelope.
        @test all(Float64.(peaks[:,5]) .>= 1.)
        @test Set(peaks[:,4]) == Set(["minimum","maximum"])
        case_peaks = readdlm(joinpath(sweep.output_directory,"case_0002_local_extrema.csv"), ','; skipstart=1)
        reference = simulate(41.)
        for r in eachrow(case_peaks)
            index = findfirst(==(Float64(r[5])), reference.time)
            @test !isnothing(index)
            @test r[6] == getproperty(reference.responses, Symbol(r[2]))[index]
        end
        # These keywords are also passed directly by the editable runner.
        paths = plot_pazy_bifurcation(sweep.summary_path; xlims=(40.,50.),
            linewidth=0.55, markersize=1.35, markerstrokewidth=0.4,
            fontfamily="Times", guidefontsize=12, tickfontsize=11,
            legendfontsize=10, captionfontsize=11, captionfontfamily="Times Bold",
            figure_size=(1100,860), panel_size=(550,420), dpi=300, gridalpha=0.18)
        @test length(paths) == 10
        @test all(p -> isfile(p) && filesize(p)>1000, paths)
        @test PazyBifurcationPlots.saved_incidence(sweep.summary_path) ≈ 5.
        # Old summaries lack peaks and the fourth response. Reconstruct only the
        # exactly known geometric incidence; explicitly label envelope markers.
        legacy = joinpath(directory,"legacy")
        mkpath(legacy)
        rows = filter(r -> r.response != "tip_geometric_incidence_deg", sweep.rows)
        legacy_summary = joinpath(legacy,"bifurcation_summary.csv")
        TimeDomainBifurcation.write_summary(legacy_summary, rows)
        legacy_paths = @test_logs (:warn, r"Local-extrema CSV unavailable") plot_pazy_bifurcation(
            legacy_summary; angle_of_attack_deg=5., formats=("png",), save_panels=false)
        @test length(legacy_paths) == 1 && filesize(first(legacy_paths))>1000
        # Header-only peak CSVs are valid for strictly stationary signals.
        stationary = run_bifurcation(U -> (; time=[0.,1.,2.],
                responses=(; (n => fill(1.,3) for n in keys(definitions))...),
                static_equilibrium=(; (n => 1. for n in keys(definitions))...));
            U_min=40., U_max=40., dU=1., definitions, simulation_time=2., transient_time=1.,
            plot_results=false, output_directory=joinpath(directory,"stationary"))
        @test length(readlines(stationary.local_extrema_path)) == 1
        stationary_paths = plot_pazy_bifurcation(stationary.summary_path;
            angle_of_attack_deg=5., formats=("png",), save_panels=false)
        @test length(stationary_paths) == 1 && filesize(first(stationary_paths))>1000
        @test !isdefined(Main, :AeroBeams) # Plotting does not load a solver.
        println("Synthetic figure validation: ", sweep.output_directory)
    finally
        "--keep-output" in ARGS || rm(directory; recursive=true)
    end
end
