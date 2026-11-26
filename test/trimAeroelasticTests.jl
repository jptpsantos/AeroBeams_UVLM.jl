@testset "Trim analysis of the Blended-Wing-Body flying wing in free flight" begin
    include_example("examples/trimAeroelastic/BWBtrim.jl")
    @test matches_test_data(trimAoA, "trimAeroelastic/BWBtrim", "trimAoA.txt")
    @test matches_test_data(trimThrust, "trimAeroelastic/BWBtrim", "trimThrust.txt")
    @test matches_test_data(trimδ, "trimAeroelastic/BWBtrim", "trimDelta.txt")
end

@testset "Trim analysis the conventional HALE aircraft in free flight (considering aerodynamics from stabilizers and thrust)" begin
    include_example("examples/trimAeroelastic/cHALEfullTrim.jl")
    @test matches_test_data(trimAoA, "trimAeroelastic/cHALEfullTrim", "trimAoA.txt")
    @test matches_test_data(trimThrust, "trimAeroelastic/cHALEfullTrim", "trimThrust.txt")
    @test matches_test_data(trimδ, "trimAeroelastic/cHALEfullTrim", "trimDelta.txt")
end

@testset "Trim analysis the conventional HALE aircraft in free flight at rigid and flexible configurations (neglecting aerodynamics from stabilizers and thrust)" begin
    include_example("examples/trimAeroelastic/cHALEtrim.jl")
    @test matches_test_data(trimAoA, "trimAeroelastic/cHALEtrim", "trimAoA.txt")
end

@testset "Trim analysis of the Helios flying-wing" begin
    include_example("examples/trimAeroelastic/heliosTrim.jl")
    @test matches_test_data(trimAoA, "trimAeroelastic/heliosTrim", "trimAoA.txt")
    @test matches_test_data(trimThrust, "trimAeroelastic/heliosTrim", "trimThrust.txt")
    @test matches_test_data(trimδ, "trimAeroelastic/heliosTrim", "trimDelta.txt")
end