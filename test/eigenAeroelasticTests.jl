# Seed the random initial vector of the Krylov eigensolver for reproducible modes
using Random

!is_ci && @testset "Flutter analysis of the Blended-Wing-Body flying wing" begin
    Random.seed!(1234)
    include_example("examples/eigenAeroelastic/BWBflutter.jl")
    @test matches_test_data(freqsOsc, "eigenAeroelastic/BWBflutter", "freqs.txt")
    @test matches_test_data(dampsOsc, "eigenAeroelastic/BWBflutter", "damps.txt")
end

!is_ci && @testset "Flutter analysis of the matched wing of the Blended-Wing-Body flying wing" begin
    Random.seed!(1234)
    include_example("examples/eigenAeroelastic/BWBmatchedWingFlutter.jl")
    @test matches_test_data(freqsOsc, "eigenAeroelastic/BWBmatchedWingFlutter", "freqsOsc.txt")
    @test matches_test_data(dampsOsc, "eigenAeroelastic/BWBmatchedWingFlutter", "dampsOsc.txt")
    @test matches_test_data(trimAoA, "eigenAeroelastic/BWBmatchedWingFlutter", "trimAoA.txt")
    @test matches_test_data(trimThrust, "eigenAeroelastic/BWBmatchedWingFlutter", "trimThrust.txt")
    @test matches_test_data(trimδ, "eigenAeroelastic/BWBmatchedWingFlutter", "trimDelta.txt")
    @test matches_test_data(flutterSpeeds, "eigenAeroelastic/BWBmatchedWingFlutter", "flutterSpeeds.txt")
end

!is_ci && @testset "Flutter analysis of the conventional HALE aircraft in free flight with structural stiffness as the varying parameter" begin
    Random.seed!(1234)
    include_example("examples/eigenAeroelastic/cHALELambdaRange.jl")
    @test matches_test_data(freqsOsc, "eigenAeroelastic/cHALELambdaRange", "freqs.txt")
    @test matches_test_data(dampsOsc, "eigenAeroelastic/cHALELambdaRange", "damps.txt")
end

!is_ci && @testset "Flutter analysis of the conventional HALE aircraft in free flight with airspeed and structural stiffness as the varying parameters" begin
    Random.seed!(1234)
    include_example("examples/eigenAeroelastic/cHALELURange.jl")
    @test matches_test_data(freqsOsc, "eigenAeroelastic/cHALELURange", "freqs.txt")
    @test matches_test_data(dampsOsc, "eigenAeroelastic/cHALELURange", "damps.txt")
end

@testset "Flutter analysis of the conventional HALE aircraft in free flight with airspeed as the varying parameter" begin
    Random.seed!(1234)
    include_example("examples/eigenAeroelastic/cHALEURange.jl")
    @test matches_test_data(freqsOsc, "eigenAeroelastic/cHALEURange", "freqs.txt")
    @test matches_test_data(dampsOsc, "eigenAeroelastic/cHALEURange", "damps.txt")
end

@testset "Flutter analysis of the Goland wing" begin
    Random.seed!(1234)
    include_example("examples/eigenAeroelastic/GolandWingFlutter.jl")
    @test matches_test_data(freqsOsc, "eigenAeroelastic/GolandWingFlutter", "freqs.txt")
    @test matches_test_data(dampsOsc, "eigenAeroelastic/GolandWingFlutter", "damps.txt")
end

@testset "Flutter analysis of the baseline Healy free FFWT wing with tip loss, root pitch and airspeed as the varying parameters" begin
    Random.seed!(1234)
    include_example("examples/eigenAeroelastic/HealyBaselineFFWTfreeFlutterAoARangeURange.jl")
    @test matches_test_data(freqsOsc, "eigenAeroelastic/HealyBaselineFFWTfreeFlutterAoARangeURange", "freqs.txt")
    @test matches_test_data(dampsOsc, "eigenAeroelastic/HealyBaselineFFWTfreeFlutterAoARangeURange", "damps.txt")
end

!is_ci && @testset "Flutter analysis of the baseline Healy free FFWT wing with flare angle and airspeed as the varying parameters" begin
    Random.seed!(1234)
    include_example("examples/eigenAeroelastic/HealyBaselineFFWTfreeFlutterFlareRangeURange.jl")
    @test matches_test_data(freqsOsc, "eigenAeroelastic/HealyBaselineFFWTfreeFlutterFlareRangeURange", "freqs.txt")
    @test matches_test_data(dampsOsc, "eigenAeroelastic/HealyBaselineFFWTfreeFlutterFlareRangeURange", "damps.txt")
end

!is_ci && @testset "Flutter analysis of the baseline Healy locked FFWT wing with airspeed as the varying parameters" begin
    Random.seed!(1234)
    include_example("examples/eigenAeroelastic/HealyBaselineFFWTlockedFlutterURange.jl")
    @test matches_test_data(freqsOsc, "eigenAeroelastic/HealyBaselineFFWTlockedFlutterURange", "freqs.txt")
    @test matches_test_data(dampsOsc, "eigenAeroelastic/HealyBaselineFFWTlockedFlutterURange", "damps.txt")
end

!is_ci && @testset "Flutter analysis of the Healy LCO FFWT wing with root pitch angle as the varying parameter" begin
    Random.seed!(1234)
    include_example("examples/eigenAeroelastic/HealyLCOFFWTfreeFlutterAoARange.jl")
    @test matches_test_data(freqsOsc, "eigenAeroelastic/HealyLCOFFWTfreeFlutterAoARange", "freqs.txt")
    @test matches_test_data(dampsOsc, "eigenAeroelastic/HealyLCOFFWTfreeFlutterAoARange", "damps.txt")
end

!is_ci && @testset "Flutter analysis of the Healy FFWT wing with sideslip angle as the varying parameter" begin
    Random.seed!(1234)
    include_example("examples/eigenAeroelastic/HealySideslipFFWTflutterSideslipRange.jl")
    @test matches_test_data(freqsOsc, "eigenAeroelastic/HealySideslipFFWTflutterSideslipRange", "freqs.txt")
    @test matches_test_data(dampsOsc, "eigenAeroelastic/HealySideslipFFWTflutterSideslipRange", "damps.txt")
end

!is_ci && @testset "Flutter analysis of the Healy sideslip FFWT wing with airspeed as the varying parameter" begin
    Random.seed!(1234)
    include_example("examples/eigenAeroelastic/HealySideslipFFWTflutterURange.jl")
    @test matches_test_data(freqsOsc, "eigenAeroelastic/HealySideslipFFWTflutterURange", "freqs.txt")
    @test matches_test_data(dampsOsc, "eigenAeroelastic/HealySideslipFFWTflutterURange", "damps.txt")
end

!is_ci && @testset "Flutter analysis of the Helios flying-wing in free flight with payload and structural stiffness as the varying parameters" begin
    Random.seed!(1234)
    include_example("examples/eigenAeroelastic/heliosFlutterPLambdaRange.jl")
    @test matches_test_data(freqsOsc, "eigenAeroelastic/heliosFlutterPLambdaRange", "freqs.txt")
    @test matches_test_data(dampsOsc, "eigenAeroelastic/heliosFlutterPLambdaRange", "damps.txt")
end

!is_ci && @testset "Flutter analysis of the Helios flying-wing in free flight with payload as the varying parameter" begin
    Random.seed!(1234)
    include_example("examples/eigenAeroelastic/heliosFlutterPRange.jl")
    @test matches_test_data(freqsOsc, "eigenAeroelastic/heliosFlutterPRange", "freqs.txt")
    @test matches_test_data(dampsOsc, "eigenAeroelastic/heliosFlutterPRange", "damps.txt")
end

@testset "Flutter analysis of the Helios flying-wing in free flight with airspeed as the varying parameter" begin
    Random.seed!(1234)
    include_example("examples/eigenAeroelastic/heliosFlutterURange.jl")
    @test matches_test_data(freqsOsc, "eigenAeroelastic/heliosFlutterURange", "freqs.txt")
    @test matches_test_data(dampsOsc, "eigenAeroelastic/heliosFlutterURange", "damps.txt")
end

@testset "Flutter analysis of the wing of the Helios flying-wing" begin
    Random.seed!(1234)
    include_example("examples/eigenAeroelastic/heliosWingFlutter.jl")
    @test matches_test_data(freqsOsc, "eigenAeroelastic/heliosWingFlutter", "freqs.txt")
    @test matches_test_data(dampsOsc, "eigenAeroelastic/heliosWingFlutter", "damps.txt")
end

!is_ci && @testset "Flutter analysis of the Pazy FFWT wing with airspeed as the varying parameter" begin
    Random.seed!(1234)
    include_example("examples/eigenAeroelastic/PazyFFWTflutterURange.jl")
    @test matches_test_data(freqsOsc, "eigenAeroelastic/PazyFFWTflutterURange", "freqs.txt")
    @test matches_test_data(dampsOsc, "eigenAeroelastic/PazyFFWTflutterURange", "damps.txt")
    @test matches_test_data(ϕHinge, "eigenAeroelastic/PazyFFWTflutterURange", "phiHinge.txt")
end

@testset "Flutter analysis of the Pazy wing" begin
    Random.seed!(1234)
    include_example("examples/eigenAeroelastic/PazyWingFlutter.jl")
    @test matches_test_data(flutterOnsetSpeeds, "eigenAeroelastic/PazyWingFlutter", "flutterOnsetSpeeds.txt")
    @test matches_test_data(flutterOnsetFreqs, "eigenAeroelastic/PazyWingFlutter", "flutterOnsetFreqs.txt")
    @test matches_test_data(flutterOnsetDisp, "eigenAeroelastic/PazyWingFlutter", "flutterOnsetDisp.txt")
    @test matches_test_data(flutterOffsetSpeeds, "eigenAeroelastic/PazyWingFlutter", "flutterOffsetSpeeds.txt")
    @test matches_test_data(flutterOffsetFreqs, "eigenAeroelastic/PazyWingFlutter", "flutterOffsetFreqs.txt")
    @test matches_test_data(flutterOffsetDisp, "eigenAeroelastic/PazyWingFlutter", "flutterOffsetDisp.txt")
end

!is_ci && @testset "Flutter and divergence analysis of the Pazy wing" begin
    Random.seed!(1234)
    include_example("examples/eigenAeroelastic/PazyWingFlutterAndDivergence.jl")
    @test matches_test_data(flutterOnsetSpeeds, "eigenAeroelastic/PazyWingFlutterAndDivergence", "flutterOnsetSpeeds.txt")
    @test matches_test_data(flutterOnsetFreqs, "eigenAeroelastic/PazyWingFlutterAndDivergence", "flutterOnsetFreqs.txt")
end

!is_ci && @testset "Flutter analysis of the Pazy wing with varying root pitch angle" begin
    Random.seed!(1234)
    include_example("examples/eigenAeroelastic/PazyWingFlutterPitchRange.jl")
    @test matches_test_data(freqsOsc, "eigenAeroelastic/PazyWingFlutterPitchRange", "freqs.txt")
    @test matches_test_data(dampsOsc, "eigenAeroelastic/PazyWingFlutterPitchRange", "damps.txt")
end

@testset "Flutter analysis of the Pazy wing with the wingtip rod modeled as a beam" begin
    Random.seed!(1234)
    include_example("examples/eigenAeroelastic/PazyWingBeamRodFlutter.jl")
    @test matches_test_data(flutterOnsetSpeeds, "eigenAeroelastic/PazyWingBeamRodFlutter", "flutterOnsetSpeeds.txt")
    @test matches_test_data(flutterOnsetFreqs, "eigenAeroelastic/PazyWingBeamRodFlutter", "flutterOnsetFreqs.txt")
    @test matches_test_data(flutterOffsetSpeeds, "eigenAeroelastic/PazyWingBeamRodFlutter", "flutterOffsetSpeeds.txt")
    @test matches_test_data(flutterOffsetFreqs, "eigenAeroelastic/PazyWingBeamRodFlutter", "flutterOffsetFreqs.txt")
end

@testset "Flutter analysis of the Pazy wing without the wingtip rod" begin
    Random.seed!(1234)
    include_example("examples/eigenAeroelastic/PazyWingNoRodFlutter.jl")
    # Discard the events of the fifth tracked mode (undamped, at ~125 Hz), which are intermittent: its damping at low airspeeds is numerical noise
    keepOnset = flutterOnsetModes .!= 5
    keepOffset = flutterOffsetModes .!= 5
    @test matches_test_data(flutterOnsetSpeeds[keepOnset], "eigenAeroelastic/PazyWingNoRodFlutter", "flutterOnsetSpeeds.txt")
    @test matches_test_data(flutterOnsetFreqs[keepOnset], "eigenAeroelastic/PazyWingNoRodFlutter", "flutterOnsetFreqs.txt")
    @test matches_test_data(flutterOffsetSpeeds[keepOffset], "eigenAeroelastic/PazyWingNoRodFlutter", "flutterOffsetSpeeds.txt")
    @test matches_test_data(flutterOffsetFreqs[keepOffset], "eigenAeroelastic/PazyWingNoRodFlutter", "flutterOffsetFreqs.txt")
end

!is_ci && @testset "Flutter analysis of the Pazy wing with the wingtip rod modeled as a beam with varying root pitch angle" begin
    Random.seed!(1234)
    include_example("examples/eigenAeroelastic/PazyWingBeamRodFlutterPitchRange.jl")
    @test matches_test_data(freqsOsc, "eigenAeroelastic/PazyWingBeamRodFlutterPitchRange", "freqs.txt")
    @test matches_test_data(dampsOsc, "eigenAeroelastic/PazyWingBeamRodFlutterPitchRange", "damps.txt")
end

!is_ci && @testset "Flutter analysis of the Pazy wing without the wingtip rod with varying root pitch angle" begin
    Random.seed!(1234)
    include_example("examples/eigenAeroelastic/PazyWingNoRodFlutterPitchRange.jl")
    @test matches_test_data(freqsOsc, "eigenAeroelastic/PazyWingNoRodFlutterPitchRange", "freqs.txt")
    @test matches_test_data(dampsOsc, "eigenAeroelastic/PazyWingNoRodFlutterPitchRange", "damps.txt")
end

@testset "Flutter analysis of the Pazy wing with varying tip mass positions" begin
    Random.seed!(1234)
    include_example("examples/eigenAeroelastic/PazyWingFlutterTipMassRange.jl")
    @test matches_test_data(freqsOsc, "eigenAeroelastic/PazyWingFlutterTipMassRange", "freqs.txt")
    @test matches_test_data(dampsOsc, "eigenAeroelastic/PazyWingFlutterTipMassRange", "damps.txt")
end

@testset "Frequency analysis of the Pazy wing" begin
    Random.seed!(1234)
    include_example("examples/eigenAeroelastic/PazyWingFreqsEvolution.jl")
    @test matches_test_data(freqsOsc, "eigenAeroelastic/PazyWingFreqsEvolution", "freqs.txt")
end

!is_ci && @testset "Flutter analysis of the sixteen-meter-wing" begin
    Random.seed!(1234)
    include_example("examples/eigenAeroelastic/SMWFlutter.jl")
    @test matches_test_data(x1_def[end]/L, "eigenAeroelastic/SMWFlutter", "x1_def.txt")
    @test matches_test_data(x3_def[end]/L, "eigenAeroelastic/SMWFlutter", "x3_def.txt")
    @test matches_test_data(α_of_x1[end], "eigenAeroelastic/SMWFlutter", "alpha.txt")
    @test matches_test_data(freqsOsc, "eigenAeroelastic/SMWFlutter", "freqs.txt")
    @test matches_test_data(dampsOsc, "eigenAeroelastic/SMWFlutter", "damps.txt")
end

!is_ci && @testset "Flutter boundary analysis of the sixteen-meter-wing as a function of the pitch angle" begin
    Random.seed!(1234)
    include_example("examples/eigenAeroelastic/SMWFlutterPitchRange.jl")
    @test matches_test_data(freqsOsc, "eigenAeroelastic/SMWFlutterPitchRange", "freqs.txt")
    @test matches_test_data(dampsOsc, "eigenAeroelastic/SMWFlutterPitchRange", "damps.txt")
end
    
!is_ci && @testset "Flutter boundary analysis of the sixteen-meter-wing as a function of the wing curvature with tip load as the varying parameter" begin
    Random.seed!(1234)
    include_example("examples/eigenAeroelastic/SMWFlutterPrecurvatureRange.jl")
    @test matches_test_data(flutterSpeed, "eigenAeroelastic/SMWFlutterPrecurvatureRange", "flutterSpeed.txt")
end

!is_ci && @testset "Flutter boundary analysis of the sixteen-meter-wing as a function of the wing curvature with root angle as the varying parameter" begin
    Random.seed!(1234)
    include_example("examples/eigenAeroelastic/SMWFlutterPrecurvatureRange2.jl")
    @test matches_test_data(freqsOsc, "eigenAeroelastic/SMWFlutterPrecurvatureRange2", "freqs.txt")
    @test matches_test_data(dampsOsc, "eigenAeroelastic/SMWFlutterPrecurvatureRange2", "damps.txt")
end

!is_ci && @testset "Flutter boundary analysis of the sixteen-meter-wing as a function of the bending-torsion coupling factor" begin
    Random.seed!(1234)
    include_example("examples/eigenAeroelastic/SMWFlutterStructuralCouplingRange.jl")
    @test matches_test_data(flutterSpeed, "eigenAeroelastic/SMWFlutterStructuralCouplingRange", "flutterSpeed.txt")
end

!is_ci && @testset "Flutter boundary analysis of the sixteen-meter-wing as a function of the tip displacement" begin
    Random.seed!(1234)
    include_example("examples/eigenAeroelastic/SMWFlutterTipDispRange.jl")
    @test matches_test_data(flutterSpeed, "eigenAeroelastic/SMWFlutterTipDispRange", "flutterSpeed.txt")
    @test matches_test_data(flutterFreq, "eigenAeroelastic/SMWFlutterTipDispRange", "flutterFreq.txt")
end

!is_ci && @testset "Linear flutter analysis of the sixteen-meter-wing" begin
    Random.seed!(1234)
    include_example("examples/eigenAeroelastic/SMWLinearFlutter.jl")
    @test matches_test_data(flutterSpeed, "eigenAeroelastic/SMWLinearFlutter", "flutterSpeed.txt")
    @test matches_test_data(flutterFreq, "eigenAeroelastic/SMWLinearFlutter", "flutterFreq.txt")
end

!is_ci && @testset "Coupled bending-torsion divergence analysis of a straight wing" begin
    Random.seed!(1234)
    include_example("examples/eigenAeroelastic/straightWingBendingTorsionalDivergence.jl")
    @test matches_test_data(dampsNonOsc, "eigenAeroelastic/straightWingBendingTorsionalDivergence", "dampsNonOsc.txt")
    @test matches_test_data(UD, "eigenAeroelastic/straightWingBendingTorsionalDivergence", "UD.txt")
end

!is_ci && @testset "Torsional divergence analysis of a straight wing" begin
    Random.seed!(1234)
    include_example("examples/eigenAeroelastic/straightWingTorsionalDivergence.jl")
    @test matches_test_data(dampsNonOsc, "eigenAeroelastic/straightWingTorsionalDivergence", "damps.txt")
end

!is_ci && @testset "Flutter analysis of the swept Pazy wing" begin
    Random.seed!(1234)
    include_example("examples/eigenAeroelastic/sweptPazyFlutterPitchRange.jl")
    @test matches_test_data(freqsOsc, "eigenAeroelastic/sweptPazyFlutterPitchRange", "freqs.txt")
    @test matches_test_data(dampsOsc, "eigenAeroelastic/sweptPazyFlutterPitchRange", "damps.txt")
end

@testset "Flutter analysis of the undeformed swept Pazy wing" begin
    Random.seed!(1234)
    include_example("examples/eigenAeroelastic/sweptPazyUndeformedFlutterSweepRange.jl")
    @test matches_test_data(freqsOsc, "eigenAeroelastic/sweptPazyUndeformedFlutterSweepRange", "freqs.txt")
    @test matches_test_data(dampsOsc, "eigenAeroelastic/sweptPazyUndeformedFlutterSweepRange", "damps.txt")
end

@testset "Bending divergence analyses of swept wing" begin
    Random.seed!(1234)
    include_example("examples/eigenAeroelastic/sweptWingBendingDivergenceSweepRange.jl")
    @test matches_test_data(UD, "eigenAeroelastic/sweptWingBendingDivergenceSweepRange", "UD.txt")
end

@testset "Coupled bending-torsion divergence analyses of swept wings" begin
    Random.seed!(1234)
    include_example("examples/eigenAeroelastic/sweptWingBendingTorsionalDivergenceSweepRange.jl")
    @test matches_test_data(UD, "eigenAeroelastic/sweptWingBendingTorsionalDivergenceSweepRange", "UD.txt")
end

@testset "Torsional divergence analyses of swept wings" begin
    Random.seed!(1234)
    include_example("examples/eigenAeroelastic/sweptWingTorsionalDivergenceSweepRange.jl")
    @test matches_test_data(UD, "eigenAeroelastic/sweptWingTorsionalDivergenceSweepRange", "UD.txt")
end

@testset "Bending divergence analysis of a swept wing" begin
    Random.seed!(1234)
    include_example("examples/eigenAeroelastic/sweptWingBendingDivergence.jl")
    @test matches_test_data(dampsNonOsc, "eigenAeroelastic/sweptWingBendingDivergence", "dampsNonOsc.txt")
    @test matches_test_data(UD, "eigenAeroelastic/sweptWingBendingDivergence", "UD.txt")
end

@testset "Coupled bending-torsion divergence analysis of a swept wing" begin
    Random.seed!(1234)
    include_example("examples/eigenAeroelastic/sweptWingBendingTorsionalDivergence.jl")
    @test matches_test_data(dampsNonOsc, "eigenAeroelastic/sweptWingBendingTorsionalDivergence", "dampsNonOsc.txt")
    @test matches_test_data(UD, "eigenAeroelastic/sweptWingBendingTorsionalDivergence", "UD.txt")
end

@testset "Torsional divergence analysis of a swept wing" begin
    Random.seed!(1234)
    include_example("examples/eigenAeroelastic/sweptWingTorsionalDivergence.jl")
    @test matches_test_data(dampsNonOsc, "eigenAeroelastic/sweptWingTorsionalDivergence", "dampsNonOsc.txt")
    @test matches_test_data(UD, "eigenAeroelastic/sweptWingTorsionalDivergence", "UD.txt")
end

@testset "Steady aeroelastic analysis of the Tang&Dowell wing at varying airspeed" begin
    Random.seed!(1234)
    include_example("examples/eigenAeroelastic/TDWingAirspeedRange.jl")
    @test matches_test_data(tip_u3, "eigenAeroelastic/TDWingAirspeedRange", "tip_u3.txt")
    @test matches_test_data(tip_twist, "eigenAeroelastic/TDWingAirspeedRange", "tip_twist.txt")
end

@testset "Flutter and divergence analysis of a typical section" begin
    Random.seed!(1234)
    include_example("examples/eigenAeroelastic/typicalSectionFlutterAndDivergence.jl")
    @test matches_test_data(freqsOsc, "eigenAeroelastic/typicalSectionFlutterAndDivergence", "freqsOsc.txt")
    @test matches_test_data(dampsOsc, "eigenAeroelastic/typicalSectionFlutterAndDivergence", "dampsOsc.txt")
    @test matches_test_data(freqsNonOsc, "eigenAeroelastic/typicalSectionFlutterAndDivergence", "freqsNonOsc.txt")
    @test matches_test_data(dampsNonOsc, "eigenAeroelastic/typicalSectionFlutterAndDivergence", "dampsNonOsc.txt")
end

@testset "Torsional divergence analysis of a typical section" begin
    Random.seed!(1234)
    include_example("examples/eigenAeroelastic/typicalSectionPitchDivergence.jl")
    @test matches_test_data(dampsNonOsc, "eigenAeroelastic/typicalSectionPitchDivergence", "damps.txt")
end