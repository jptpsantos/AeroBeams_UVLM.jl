ENV["GKSwstype"] = "100"
using Test, DelimitedFiles
include(joinpath(@__DIR__, "..", "..", "lib", "WingPropellerUVLM", "examples",
    "pazy_aerobeams_coupled", "plot_pazy_bifurcation_interpolated.jl"))
const P = PazyInterpolatedBifurcation

@testset "Shape-preserving interpolation" begin
    x = [0., 0.3, 1.5, 2., 4.]
    for y in ([1., 2., 2., 3., 5.], [0., 1., 0., 0., 2.], fill(3.,5), 2 .* x .+ 1)
        @test P.pchip(x, y, x) ≈ y
        for i in 1:length(x)-1
            q = collect(range(x[i],x[i+1]; length=101))
            v = P.pchip(x,y,q)
            @test all(min(y[i],y[i+1])-1e-12 .<= v .<= max(y[i],y[i+1])+1e-12)
        end
    end
    q = collect(range(0.,4.; length=101))
    @test P.pchip(x, 2 .* x .+ 1, q) ≈ 2 .* q .+ 1
    @test P.pchip([2.], [5.], [2.]) == [5.]
    @test P.pchip([0.,2.], [1.,5.], [0.,1.,2.]) == [1.,3.,5.]
    @test_throws ArgumentError P.pchip(x, x, [-0.1])
    @test_throws ArgumentError P.pchip([0.,0.], [1.,2.], [0.])
    @test_throws ArgumentError P.pchip(x, [NaN,1.,2.,3.,4.], x)
end

@testset "Ordered envelopes, exact knots and gaps" begin
    speed = [40.,41.,42.,43.,44.]
    lower = [2.,3.,1.,1.,4.]
    upper = [2.,8.,1.1,2.,4.]
    steady = [2.,3.,NaN,1.2,4.]
    curves = P.interpolate_envelope(speed,lower,upper,steady)
    @test length(curves) == 1
    c = only(curves)
    @test all(c.minimum .<= c.maximum)
    for (k,u) in enumerate(speed)
        i = findfirst(==(u),c.speed)
        @test c.minimum[i] ≈ lower[k]
        @test c.maximum[i] ≈ upper[k]
    end
    @test all(isnan, c.static[(c.speed .> 41) .& (c.speed .< 43)])
    @test extrema(c.speed) == extrema(speed)
    with_nan = P.interpolate_envelope(speed,[2.,3.,NaN,1.,4.],upper,steady)
    @test length(with_nan) == 2
    @test last(with_nan[1].speed) == 41.
    @test first(with_nan[2].speed) == 43.
    @test P.valid_segments([1.,2.,3.,8.,9.], trues(5)) == [1:3,4:5]
    @test_throws ArgumentError P.interpolate_envelope(speed,upper,lower,steady)
    @test_throws ArgumentError P.interpolate_envelope(speed,lower,upper,steady; points_per_interval=1)
end

@testset "Saved cases and AeroBeams palette" begin
    paths = P.saved_summaries()
    @test length(paths) == 4
    for path in paths
        angle = P.Paper.saved_incidence(path)
        case = P.load_case(path,angle; points_per_interval=32,max_gap=nothing)
        @test case.angle ≈ 5.
        @test case.color == P.AEROBEAMS_COLORS[3]
        @test length(case.panels) == 4
        for d in case.panels, c in d.curves
            @test all(isfinite,c.minimum) && all(isfinite,c.maximum)
            @test all(c.minimum .<= c.maximum)
            @test first(c.speed) >= first(d.speed)
            @test last(c.speed) <= last(d.speed)
            for k in eachindex(d.speed)
                i = findfirst(==(d.speed[k]),c.speed)
                isnothing(i) && continue
                @test isapprox(c.minimum[i],d.lower[k]; atol=1e-10)
                @test isapprox(c.maximum[i],d.upper[k]; atol=1e-10)
            end
        end
        twist, geometric = case.panels[3],case.panels[4]
        @test geometric.lower ≈ twist.lower .+ angle
        @test geometric.upper ≈ twist.upper .+ angle
    end
    @test !isdefined(Main,:AeroBeams)
    println("Validated palette for 5 degrees: ", P.incidence_color(5.))
end

P.main(["--all"])
