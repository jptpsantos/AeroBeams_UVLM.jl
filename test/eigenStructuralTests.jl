@testset "Modal analysis of the axial vibration of a beam under clamped-clamped boundary conditions" begin
    include_example("examples/eigenStructural/beamAxialVibrationCC.jl")
    @test matches_test_data(freqs, "eigenStructural/beamAxialVibrationCC", "freqs.txt")
end

@testset "Modal analysis of the axial vibration of a beam under clamped-free boundary conditions" begin
    include_example("examples/eigenStructural/beamAxialVibrationCF.jl")
    @test matches_test_data(freqs, "eigenStructural/beamAxialVibrationCF", "freqs.txt")
end

@testset "Modal analysis of the axial vibration of a beam under free-free boundary conditions" begin
    include_example("examples/eigenStructural/beamAxialVibrationFF.jl")
    @test matches_test_data(freqs, "eigenStructural/beamAxialVibrationFF", "freqs.txt")
end

@testset "Modal analysis of the bending vibration of a beam under clamped-clamped boundary conditions" begin
    include_example("examples/eigenStructural/beamBendingVibrationCC.jl")
    @test matches_test_data(freqs, "eigenStructural/beamBendingVibrationCC", "freqs.txt")
end

@testset "Modal analysis of the bending vibration of a beam under clamped-free boundary conditions" begin
    include_example("examples/eigenStructural/beamBendingVibrationCF.jl")
    @test matches_test_data(freqs, "eigenStructural/beamBendingVibrationCF", "freqs.txt")
end

@testset "Modal analysis of the bending vibration of a beam under clamped-pinned boundary conditions" begin
    include_example("examples/eigenStructural/beamBendingVibrationCP.jl")
    @test matches_test_data(freqs, "eigenStructural/beamBendingVibrationCP", "freqs.txt")
end

@testset "Modal analysis of the bending vibration of a beam under clamped-sliding boundary conditions" begin
    include_example("examples/eigenStructural/beamBendingVibrationCS.jl")
    @test matches_test_data(freqs, "eigenStructural/beamBendingVibrationCS", "freqs.txt")
end

@testset "Modal analysis of the bending vibration of a beam under free-free boundary conditions" begin
    include_example("examples/eigenStructural/beamBendingVibrationFF.jl")
    @test matches_test_data(freqs, "eigenStructural/beamBendingVibrationFF", "freqs.txt")
end

@testset "Modal analysis of the bending vibration of a beam under pinned-pinned boundary conditions" begin
    include_example("examples/eigenStructural/beamBendingVibrationPP.jl")
    @test matches_test_data(freqs, "eigenStructural/beamBendingVibrationPP", "freqs.txt")
end

@testset "Modal analysis of the torsional vibration of a beam under clamped-clamped boundary conditions" begin
    include_example("examples/eigenStructural/beamTorsionalVibrationCC.jl")
    @test matches_test_data(freqs, "eigenStructural/beamTorsionalVibrationCC", "freqs.txt")
end

@testset "Modal analysis of the torsional vibration of a beam under clamped-free boundary conditions" begin
    include_example("examples/eigenStructural/beamTorsionalVibrationCF.jl")
    @test matches_test_data(freqs, "eigenStructural/beamTorsionalVibrationCF", "freqs.txt")
end

@testset "Modal analysis of the torsional vibration of a beam under free-free boundary conditions" begin
    include_example("examples/eigenStructural/beamTorsionalVibrationFF.jl")
    @test matches_test_data(freqs, "eigenStructural/beamTorsionalVibrationFF", "freqs.txt")
end

@testset "Modal analysis of a cantilever with stiffness-proportional structural damping" begin
    include_example("examples/eigenStructural/cantileverStructuralDampingEigen.jl")
    # Damping ratios follow ζ = μs*ω/2 (Rayleigh damping)
    @test all(abs.(ζNum ./ ζExpected .- 1) .< 1e-4)
    @test matches_test_data(ωNum, "eigenStructural/cantileverStructuralDampingEigen", "freqs.txt")
    @test matches_test_data(ζNum, "eigenStructural/cantileverStructuralDampingEigen", "damps.txt")
end

@testset "Modal analysis of a cantilever beam with a tip axial inertia" begin
    include_example("examples/eigenStructural/cantileverWithTipAxialMassEigen.jl")
    @test matches_test_data(freqsNorm, "eigenStructural/cantileverWithTipAxialMassEigen", "freqsNorm.txt")
end

@testset "Modal analysis of a cantilever beam with a tip axial spring" begin
    include_example("examples/eigenStructural/cantileverWithTipAxialSpringEigen.jl")
    @test matches_test_data(freqsNorm, "eigenStructural/cantileverWithTipAxialSpringEigen", "freqsNorm.txt")
end

@testset "Modal analysis of a cantilever beam with a tip spring in bending" begin
    include_example("examples/eigenStructural/cantileverWithTipSpringEigen.jl")
    @test matches_test_data(freqsNorm, "eigenStructural/cantileverWithTipSpringEigen", "freqsNorm.txt")
end

@testset "Modal analysis of a cantilever beam with a tip torsional inertia" begin
    include_example("examples/eigenStructural/cantileverWithTipTorsionalInertiaEigen.jl")
    @test matches_test_data(freqsNorm, "eigenStructural/cantileverWithTipTorsionalInertiaEigen", "freqsNorm.txt")
end

@testset "Modal analysis of a cantilever beam with a tip torsional spring" begin
    include_example("examples/eigenStructural/cantileverWithTipTorsionalSpringEigen.jl")
    @test matches_test_data(freqsNorm, "eigenStructural/cantileverWithTipTorsionalSpringEigen", "freqsNorm.txt")
end

@testset "Modal analysis a beam clamped at one end, simply-supported at the other and with a tip inertia" begin
    include_example("examples/eigenStructural/clampedSSBeamWIthTipInertiaEigen.jl")
    @test matches_test_data(freqs, "eigenStructural/clampedSSBeamWIthTipInertiaEigen", "freqs.txt")
end

@testset "Modal analysis of the baseline Healy free FFWT wing without gravity" begin
    include_example("examples/eigenStructural/HealyBaselineFFWTModalFree.jl")
    @test matches_test_data(freqs, "eigenStructural/HealyBaselineFFWTModalFree", "freqs.txt")
end

@testset "Modal analysis of the baseline Healy locked FFWT wing" begin
    include_example("examples/eigenStructural/HealyLCOFFWTLockedModal.jl")
    @test matches_test_data(freqs, "eigenStructural/HealyLCOFFWTLockedModal", "freqs.txt")
end

@testset "Modal analysis of the baseline Healy FFWT wing with wingtip removed" begin
    include_example("examples/eigenStructural/HealyLCOFFWTRemovedModal.jl")
    @test matches_test_data(freqs, "eigenStructural/HealyLCOFFWTRemovedModal", "freqs.txt")
end

@testset "Modal analyses of the Pazy wing in horizontal and vertical positions" begin
    include_example("examples/eigenStructural/PazyWingModal.jl")
    @test matches_test_data(freqs, "eigenStructural/PazyWingModal", "freqs.txt")
end

@testset "Modal analysis of the Pazy wing with the wingtip rod modeled as a beam" begin
    include_example("examples/eigenStructural/PazyWingBeamRodModal.jl")
    @test matches_test_data(freqs, "eigenStructural/PazyWingBeamRodModal", "freqs.txt")
end

@testset "Modal analysis of the Pazy wing without the wingtip rod" begin
    include_example("examples/eigenStructural/PazyWingNoRodModal.jl")
    @test matches_test_data(freqs, "eigenStructural/PazyWingNoRodModal", "freqs.txt")
end

@testset "Modal analysis of the Pazy wing with 10 degrees of sweep and the wingtip rod modeled as a beam" begin
    include_example("examples/eigenStructural/PazyS10BeamRodModal.jl")
    @test matches_test_data(freqs, "eigenStructural/PazyS10BeamRodModal", "freqs.txt")
end

@testset "Modal analysis of a beam pinned at one end and transversely springed at the other" begin
    include_example("examples/eigenStructural/pinnedSpringedBeamEigen.jl")
    @test matches_test_data(freqs, "eigenStructural/pinnedSpringedBeamEigen", "freqs.txt")
end

@testset "Modal analysis of the sixteen-meter-wing" begin
    include_example("examples/eigenStructural/SMWModal.jl")
    @test matches_test_data(freqs, "eigenStructural/SMWModal", "freqs.txt")
end

@testset "Modal analysis of a straight rotor under varying angular velocities" begin
    include_example("examples/eigenStructural/straightRotor.jl")
    @test matches_test_data(freqs, "eigenStructural/straightRotor", "freqs.txt")
end

@testset "Modal analysis of the undeformed swept Pazy wing" begin
    include_example("examples/eigenStructural/sweptPazyModal.jl")
    @test matches_test_data(freqs, "eigenStructural/sweptPazyModal", "freqs.txt")
end

@testset "Modal analysis of the swept Pazy wing with varying tip displacement" begin
    include_example("examples/eigenStructural/sweptPazyModalTipDispRange.jl")
    @test matches_test_data(tipOOP, "eigenStructural/sweptPazyModalTipDispRange", "tipOOP.txt")
    @test matches_test_data(tip_twist, "eigenStructural/sweptPazyModalTipDispRange", "tip_twist.txt")
    @test matches_test_data(freqsOsc, "eigenStructural/sweptPazyModalTipDispRange", "freqs.txt")
end

@testset "Modal analysis of a swept-tip rotor under varying angular velocities" begin
    include_example("examples/eigenStructural/sweptTipRotor.jl")
    @test matches_test_data(freqs, "eigenStructural/sweptTipRotor", "freqs.txt")
end

@testset "Modal analysis of a cantilevered tapered beam" begin
    include_example("examples/eigenStructural/taperedBeamEigen.jl")
    @test matches_test_data(freqs, "eigenStructural/taperedBeamEigen", "freqs.txt")
end

@testset "Modal analysis of the Tang&Dowell wing at varying pitch angles" begin
    include_example("examples/eigenStructural/TDWingPitchRange.jl")
    @test matches_test_data(tip_u3, "eigenStructural/TDWingPitchRange", "tip_u3.txt")
    @test matches_test_data(tip_u2, "eigenStructural/TDWingPitchRange", "tip_u2.txt")
    @test matches_test_data(tip_twist, "eigenStructural/TDWingPitchRange", "tip_twist.txt")
    @test matches_test_data(freqs, "eigenStructural/TDWingPitchRange", "freqs.txt")
end

@testset "Modal analysis of a 2-story frame" begin
    include_example("examples/eigenStructural/twoStoryFrame.jl")
    @test matches_test_data(freqs, "eigenStructural/twoStoryFrame", "freqs.txt")
end