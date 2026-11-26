@testset "Trim analysis (reaction loads check) of a beam loaded at the middle" begin
    include_example("examples/trimStructural/freeBeamTrim.jl")
    @test matches_test_data(F3, "trimStructural/freeBeamTrim", "F3.txt")
    @test matches_test_data(M2, "trimStructural/freeBeamTrim", "M2.txt")
end

@testset "Trim analysis (reaction loads check) of a simply-supported beam loaded at the middle" begin
    include_example("examples/trimStructural/midDispBeamTrim.jl")
    @test matches_test_data(F3, "trimStructural/midDispBeamTrim", "F3.txt")
    @test matches_test_data(M2, "trimStructural/midDispBeamTrim", "M2.txt")
end

@testset "Trim analysis (reaction loads check) of a right-angled frame" begin
    include_example("examples/trimStructural/rightAngledFrameTrim.jl")
    @test matches_test_data(balanceHorizontalForce, "trimStructural/rightAngledFrameTrim", "balanceHorizontalForce.txt")
    @test matches_test_data(balanceVerticalForce, "trimStructural/rightAngledFrameTrim", "balanceVerticalForce.txt")
end

@testset "Trim analysis of a cantilever with tip force" begin
    include_example("examples/trimStructural/tipLoadedCantileverTrim.jl")
    @test matches_test_data(u3, "trimStructural/tipLoadedCantileverTrim", "u3.txt")
    @test matches_test_data(F3, "trimStructural/tipLoadedCantileverTrim", "F3.txt")
    @test matches_test_data(M2, "trimStructural/tipLoadedCantileverTrim", "M2.txt")
end