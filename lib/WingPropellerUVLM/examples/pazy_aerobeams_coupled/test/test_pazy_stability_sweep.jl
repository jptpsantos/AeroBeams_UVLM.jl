const PAZY_EXAMPLE_DIR = normpath(joinpath(@__DIR__,".."))

using Test,Statistics
include(joinpath(PAZY_EXAMPLE_DIR,"PazyStabilitySweepTools.jl"))
include(joinpath(PAZY_EXAMPLE_DIR,"..","..","studies","src","MovingBlockDamping.jl"))

@testset "Simple sweep" begin
    speeds = [40.0,45.0,50.0,55.0,60.0]
    growth = [-1.0,-0.2,0.4,0.2,-0.3]
    c = sweep_crossings(speeds,growth)
    @test length(c) == 2
    @test c[1].kind == "onset"
    @test c[1].speed ≈ 46 + 2 / 3
    @test c[2].kind == "offset"
    @test c[2].speed ≈ 57
    growth[3] = NaN
    @test length(sweep_crossings(speeds,growth)) == 1
    growth[4] = NaN
    @test isempty(sweep_crossings(speeds,growth))
    @test isempty(sweep_crossings([40.0,50.0],[0.0,1.0]))

    angles,onset,offset = sweep_aerobeams_boundary()
    @test angles == [0.0,0.5,1.0,2.0,3.0,4.0,5.0,6.0,7.0]
    @test all(isfinite,onset) && all(isfinite,offset)
    i = findfirst(==(3.0),angles)
    @test 40 < onset[i] < offset[i] < 70
    @test sweep_number(3) == "3"
    @test sweep_number(2.5) == "2p5"

    t = collect(0.0:0.002:3.0)
    for lambda in (-0.3,0.3)
        x = 2 .+ 0.01 .* exp.(lambda .* t) .* sin.(2pi * 25 .* t)
        fit = MovingBlockDamping.moving_block_metrics(
            t,x .- mean(x); block_duration_s=0.2,minimum_fit_r_squared=0.9,minimum_peaks=6)
        @test fit.summary.valid_for_convergence
        @test fit.summary.growth_rate_per_s ≈ lambda atol = 0.01
        @test sign(fit.summary.damping_ratio_percent) == -sign(lambda)
    end
    mktempdir() do directory
        path = joinpath(directory,"response.csv")
        write_sweep_csv(path,["t","x"],[0.0 1.0; 1.0 2.0])
        @test readdlm(path,',',Float64; skipstart=1) == [0.0 1.0; 1.0 2.0]
    end
end

# Optional: render the native comparison with synthetic UVLM boundaries.
if "--plots" in ARGS
    @testset "Native stability comparison plot" begin
        mktempdir() do directory
            labels = (; uvlm_region="Synthetic UVLM",aerobeams_onset="AeroBeams onset",
                aerobeams_offset="AeroBeams offset",test_onset_up="Test onset up",
                test_offset_up="Test offset up",test_onset_down="Test onset down",
                test_offset_down="Test offset down",umnast_loss="UM/NAST",
                umnast_panel="UM/NAST panels",sharpy="SHARPy")
            boundaries = [[3.0,45.0,"onset"],[3.0,55.0,"offset"],
                [5.0,40.0,"onset"],[5.0,50.0,"offset"]]
            figure = plot_sweep_comparison(boundaries,[3.0,5.0];
                aerobeams_color=:blue,literature_colors=(:gold,:brown,:magenta),labels)
            path = joinpath(directory,"pazy_stability_boundary.png")
            savefig(figure,path)
            @test filesize(path) > 1000
        end
    end
end
