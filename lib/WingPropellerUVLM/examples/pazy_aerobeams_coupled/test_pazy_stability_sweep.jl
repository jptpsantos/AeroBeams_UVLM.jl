using Test, Statistics
include(joinpath(@__DIR__,"PazyStabilitySweepTools.jl"))
include(joinpath(@__DIR__,"..","..","test","moving_block_damping.jl"))

@testset "Simple sweep" begin
    speeds = [40.,45.,50.,55.,60.]
    growth = [-1.,-0.2,0.4,0.2,-0.3]
    c = sweep_crossings(speeds,growth)
    @test length(c) == 2
    @test c[1].kind == "onset"
    @test c[1].speed ≈ 46+2/3
    @test c[2].kind == "offset"
    @test c[2].speed ≈ 57
    growth[3] = NaN
    @test length(sweep_crossings(speeds,growth)) == 1
    growth[4] = NaN
    @test isempty(sweep_crossings(speeds,growth))
    @test isempty(sweep_crossings([40.,50.],[0.,1.]))

    angles,onset,offset = sweep_aerobeams_boundary()
    @test length(angles) == 17
    @test all(isfinite,onset) && all(isfinite,offset)
    i = findfirst(==(3.0),angles)
    @test 40 < onset[i] < offset[i] < 70
    @test sweep_number(3) == "3"
    @test sweep_number(2.5) == "2p5"

    t = collect(0.:0.002:3.)
    for lambda in (-0.3,0.3)
        x = 2 .+ 0.01 .* exp.(lambda .* t) .* sin.(2pi*25 .* t)
        fit = MovingBlockDampingTests.MovingBlockDamping.moving_block_metrics(
            t,x .- mean(x);block_duration_s=0.2,minimum_fit_r_squared=0.9,minimum_peaks=6)
        @test fit.summary.valid_for_convergence
        @test fit.summary.growth_rate_per_s ≈ lambda atol=0.01
        @test sign(fit.summary.damping_ratio_percent) == -sign(lambda)
    end
    mktempdir() do directory
        path = joinpath(directory,"response.csv")
        write_sweep_csv(path,["t","x"],[0. 1.;1. 2.])
        @test readdlm(path,',',Float64;skipstart=1) == [0. 1.;1. 2.]
    end
end

# Optional: exercise the actual postprocessing loop using synthetic data.
if "--plots" in ARGS
    driver = joinpath(@__DIR__,"run_pazy_stability_sweep.jl")
    source = replace(read(driver,String),
        "run_simulations = true" => "run_simulations = false",
        "sweep_result = run_stability_sweep()" => "")
    include_string(Main,source,driver)
    global pitch_angles_deg = [3.0]
    global airspeeds_mps = [40.0,50.0,60.0]
    global output_directory = mktempdir(;cleanup=false)
    for (u,lambda) in zip(airspeeds_mps,[-0.3,0.3,-0.3])
        folder = joinpath(output_directory,"aoa3_V$(sweep_number(u))")
        mkpath(folder)
        t = collect(0.0:0.002:3.0)
        x = exp.(lambda .* t) .* sin.(2pi*25 .* t)
        write_sweep_csv(joinpath(folder,"response.csv"),["t","U","bending","twist"],
            hcat(t,fill(u,length(t)),0.001 .* x,0.01 .* x))
    end
    checked = run_stability_sweep()
    @test length(checked.boundary_rows) == 2
    @test all(isfinite,checked.growth)
    @test isfile(joinpath(output_directory,"pazy_stability_boundary.png"))
    println("Synthetic test plots (not physical results): $output_directory")
end
