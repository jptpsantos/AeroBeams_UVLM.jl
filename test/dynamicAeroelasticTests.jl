@testset "Dynamic analysis of the Blended-Wing-Body vehicle undergoing a checked pitch maneuver" begin
    include_example("examples/dynamicAeroelastic/BWBcheckedPitchManeuver.jl")
    @test matches_test_data(rootAoA, "dynamicAeroelastic/BWBcheckedPitchManeuver", "rootAoA.txt")
    @test matches_test_data(Δu3, "dynamicAeroelastic/BWBcheckedPitchManeuver", "Deltau3.txt")
end

!is_ci && @testset "Dynamic analysis of the conventional HALE aircraft undergoing a checked pitch maneuver" begin
    include_example("examples/dynamicAeroelastic/cHALECheckedPitchManeuver.jl")
    @test matches_test_data(wingAoA, "dynamicAeroelastic/cHALECheckedPitchManeuver", "wingAoA.txt")
    @test matches_test_data(Δu3, "dynamicAeroelastic/cHALECheckedPitchManeuver", "Deltau3.txt")
end

@testset "Dynamic analysis of the conventional HALE aircraft undergoing a coordinated turn maneuver" begin
    include_example("examples/dynamicAeroelastic/cHALECheckedRollManeuver.jl")
    @test matches_test_data(wingAoA, "dynamicAeroelastic/cHALECheckedRollManeuver", "wingAoA.txt")
    @test matches_test_data(Δu3, "dynamicAeroelastic/cHALECheckedRollManeuver", "Deltau3.txt")
end

!is_ci && @testset "Dynamic analysis of the baseline Healy FFWT wing under a series of one-minus-cosine gusts" begin
    include_example("examples/dynamicAeroelastic/HealyBaselineFFWTOMCGustFloating.jl")
    @test matches_test_data(M2root, "dynamicAeroelastic/HealyBaselineFFWTOMCGustFloating", "M2root.txt")
    @test matches_test_data(ϕ, "dynamicAeroelastic/HealyBaselineFFWTOMCGustFloating", "fold.txt")
end

@testset "Dynamic analysis of the baseline Healy FFWT wing for a specific one-minus-cosine gust" begin
    include_example("examples/dynamicAeroelastic/HealyBaselineFFWTOMCGustFloating2.jl")
    @test matches_test_data(M2root, "dynamicAeroelastic/HealyBaselineFFWTOMCGustFloating2", "M2root.txt")
    @test matches_test_data(ϕ, "dynamicAeroelastic/HealyBaselineFFWTOMCGustFloating2", "fold.txt")
end

@testset "Dynamic analysis of the Helios flying-wing undergoing a checked pitch maneuver" begin
    include_example("examples/dynamicAeroelastic/heliosCheckedPitchManeuver.jl")
    @test matches_test_data(rootAoA, "dynamicAeroelastic/heliosCheckedPitchManeuver", "rootAoA.txt")
    @test matches_test_data(Δu3, "dynamicAeroelastic/heliosCheckedPitchManeuver", "Deltau3.txt")
end

!is_ci && @testset "Dynamic analysis of the Pazy wing encountering a continuous, 1-dimensional gust defined over time" begin
    include_example("examples/dynamicAeroelastic/PazyWingContinuous1DGust.jl")
    @test matches_test_data(tipAoA, "dynamicAeroelastic/PazyWingContinuous1DGust", "tipAoA.txt")
    @test matches_test_data(tipOOP, "dynamicAeroelastic/PazyWingContinuous1DGust", "tipOOP.txt")
    @test matches_test_data(tqSpan_cn, "dynamicAeroelastic/PazyWingContinuous1DGust", "tqSpan_cn.txt")
    @test matches_test_data(tqSpan_cm, "dynamicAeroelastic/PazyWingContinuous1DGust", "tqSpan_cm.txt")
    @test matches_test_data(tqSpan_ct, "dynamicAeroelastic/PazyWingContinuous1DGust", "tqSpan_ct.txt")
end

!is_ci && @testset "Dynamic analysis of the Pazy wing encountering a continuous, 1-dimensional gust defined over space" begin
    include_example("examples/dynamicAeroelastic/PazyWingContinuous1DSpaceGust.jl")
    @test matches_test_data(tipAoA, "dynamicAeroelastic/PazyWingContinuous1DSpaceGust", "tipAoA.txt")
    @test matches_test_data(tipOOP, "dynamicAeroelastic/PazyWingContinuous1DSpaceGust", "tipOOP.txt")
    @test matches_test_data(tqSpan_cn, "dynamicAeroelastic/PazyWingContinuous1DSpaceGust", "tqSpan_cn.txt")
    @test matches_test_data(tqSpan_cm, "dynamicAeroelastic/PazyWingContinuous1DSpaceGust", "tqSpan_cm.txt")
    @test matches_test_data(tqSpan_ct, "dynamicAeroelastic/PazyWingContinuous1DSpaceGust", "tqSpan_ct.txt")
end

!is_ci && @testset "Dynamic analysis of the Pazy wing encountering a continuous, 2-dimensional gust" begin
    include_example("examples/dynamicAeroelastic/PazyWingContinuous2DSpaceGust.jl")
    @test matches_test_data(tipAoA, "dynamicAeroelastic/PazyWingContinuous2DSpaceGust", "tipAoA.txt")
    @test matches_test_data(tipOOP, "dynamicAeroelastic/PazyWingContinuous2DSpaceGust", "tipOOP.txt")
    @test matches_test_data(tqSpan_cn, "dynamicAeroelastic/PazyWingContinuous2DSpaceGust", "tqSpan_cn.txt")
    @test matches_test_data(tqSpan_cm, "dynamicAeroelastic/PazyWingContinuous2DSpaceGust", "tqSpan_cm.txt")
    @test matches_test_data(tqSpan_ct, "dynamicAeroelastic/PazyWingContinuous2DSpaceGust", "tqSpan_ct.txt")
end

!is_ci && @testset "Dynamic analysis of the Pazy wing encountering a DARPA gust" begin
    include_example("examples/dynamicAeroelastic/PazyWingDARPAGust.jl")
    @test matches_test_data(tipAoA, "dynamicAeroelastic/PazyWingDARPAGust", "tipAoA.txt")
    @test matches_test_data(tipOOP, "dynamicAeroelastic/PazyWingDARPAGust", "tipOOP.txt")
    @test matches_test_data(tqSpan_cn, "dynamicAeroelastic/PazyWingDARPAGust", "tqSpan_cn.txt")
    @test matches_test_data(tqSpan_cm, "dynamicAeroelastic/PazyWingDARPAGust", "tqSpan_cm.txt")
    @test matches_test_data(tqSpan_ct, "dynamicAeroelastic/PazyWingDARPAGust", "tqSpan_ct.txt")
end

!is_ci && @testset "Dynamic analysis of the Pazy wing encountering a one-minus-cosine gust" begin
    include_example("examples/dynamicAeroelastic/PazyWingOMCGust.jl")
    @test matches_test_data(tipAoA, "dynamicAeroelastic/PazyWingOMCGust", "tipAoA.txt")
    @test matches_test_data(tipOOP, "dynamicAeroelastic/PazyWingOMCGust", "tipOOP.txt")
    @test matches_test_data(tqSpan_cn, "dynamicAeroelastic/PazyWingOMCGust", "tqSpan_cn.txt")
    @test matches_test_data(tqSpan_cm, "dynamicAeroelastic/PazyWingOMCGust", "tqSpan_cm.txt")
    @test matches_test_data(tqSpan_ct, "dynamicAeroelastic/PazyWingOMCGust", "tqSpan_ct.txt")
end

# Reduce CI time
!is_ci && @testset "Dynamic analysis of the Pazy wing with a tip impulse force" begin
    include_example("examples/dynamicAeroelastic/PazyWingTipImpulse.jl")
    @test matches_test_data(tipAoA, "dynamicAeroelastic/PazyWingTipImpulse", "tipAoA.txt")
    @test matches_test_data(tipOOP, "dynamicAeroelastic/PazyWingTipImpulse", "tipOOP.txt")
    @test matches_test_data(tqSpan_cn, "dynamicAeroelastic/PazyWingTipImpulse", "tqSpan_cn.txt")
    @test matches_test_data(tqSpan_cm, "dynamicAeroelastic/PazyWingTipImpulse", "tqSpan_cm.txt")
    @test matches_test_data(tqSpan_ct, "dynamicAeroelastic/PazyWingTipImpulse", "tqSpan_ct.txt")
end

!is_ci && @testset "Dynamic analysis of the baseline Healy locked FFWT wing under a series of one-minus-cosine gusts" begin
    include_example("examples/dynamicAeroelastic/HealyBaselineFFWTlockedOMCGustFloating.jl")
    @test matches_test_data(M2root, "dynamicAeroelastic/HealyBaselineFFWTlockedOMCGustFloating", "M2root.txt")
end

!is_ci && @testset "Dynamic analysis of the baseline Healy free and locked FFWT wings under one-minus-cosine gusts" begin
    include_example("examples/dynamicAeroelastic/HealyBaselineFFWTlockedVfreeOMCGustFloating.jl")
    @test matches_test_data(ΔM2peaks, "dynamicAeroelastic/HealyBaselineFFWTlockedVfreeOMCGustFloating", "dM2peaks.txt")
end

!is_ci && @testset "Dynamic analysis of the Healy sideslip FFWT wing under a one-minus-cosine gust" begin
    include_example("examples/dynamicAeroelastic/HealySideslipFFWTOMCGustFloating.jl")
    @test matches_test_data(tipOOP, "dynamicAeroelastic/HealySideslipFFWTOMCGustFloating", "tipOOP.txt")
    @test matches_test_data(M2_root, "dynamicAeroelastic/HealySideslipFFWTOMCGustFloating", "M2_root.txt")
    @test matches_test_data(root_αₑ, "dynamicAeroelastic/HealySideslipFFWTOMCGustFloating", "root_alphae.txt")
    @test matches_test_data(root_cn, "dynamicAeroelastic/HealySideslipFFWTOMCGustFloating", "root_cn.txt")
end

!is_ci && @testset "Dynamic analysis of the Pazy FFWT wing under a one-minus-cosine gust" begin
    include_example("examples/dynamicAeroelastic/PazyFFWTOMCGustFloating.jl")
    @test matches_test_data(ϕ, "dynamicAeroelastic/PazyFFWTOMCGustFloating", "fold.txt")
    @test matches_test_data(tipOOP, "dynamicAeroelastic/PazyFFWTOMCGustFloating", "tipOOP.txt")
    @test matches_test_data(M2_root, "dynamicAeroelastic/PazyFFWTOMCGustFloating", "M2_root.txt")
    @test matches_test_data(tqSpan_cn, "dynamicAeroelastic/PazyFFWTOMCGustFloating", "tqSpan_cn.txt")
end

!is_ci && @testset "Dynamic analysis of a rotating wing using the dynamic stall model" begin
    include_example("examples/dynamicAeroelastic/rotorDStest.jl")
    @test matches_test_data(cl, "dynamicAeroelastic/rotorDStest", "cl.txt")
    @test matches_test_data(cm, "dynamicAeroelastic/rotorDStest", "cm.txt")
    @test matches_test_data(cdrag, "dynamicAeroelastic/rotorDStest", "cdrag.txt")
end

@testset "Dynamic analysis of a pitching wing with dynamic stall models" begin
    include_example("examples/dynamicAeroelastic/wingDStest.jl")
    @test matches_test_data(cl, "dynamicAeroelastic/wingDStest", "cl.txt")
    @test matches_test_data(cm, "dynamicAeroelastic/wingDStest", "cm.txt")
    @test matches_test_data(cdrag, "dynamicAeroelastic/wingDStest", "cdrag.txt")
end