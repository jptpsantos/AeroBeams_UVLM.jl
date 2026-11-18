using Test

script = joinpath(@__DIR__, "plot_sensitivity_percent.jl")
definitions = first(split(read(script, String), "\ndataDir = isempty(ARGS)"; limit=2))
include_string(Main, definitions, script)
data = propeller_percentage_data(ARGS[1])

@testset "Percentage interpolation" begin
    for d in data, j in 1:2
        y = d.differences[:, j]
        @test interpolate_percentages(d.x, y, d.x) ≈ y
        for i in 1:length(d.x)-1
            fine = range(d.x[i], d.x[i+1]; length=101)
            interpolated = interpolate_percentages(d.x, y, fine)
            lower, upper = extrema(y[i:i+1])
            @test all(v -> lower - 1e-12 <= v <= upper + 1e-12, interpolated)
            @test minimum(interpolated) >= -1e-12
        end
    end
    chord = only(filter(d -> d.name == "chord", data))
    step = only(filter(d -> d.name == "step", data))
    @test chord.x[1] == 5
    @test chord.differences[1, 1] ≈ 31.519236236689512
    @test step.x[1] == 1
    @test step.x[end] == 15
    @test step.differences[1, 1] ≈ 1.504440144680282
end
