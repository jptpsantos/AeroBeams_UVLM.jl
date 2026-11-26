using Random

# Flutter analysis of the Blended-Wing-Body flying wing
Random.seed!(1234)
include("../examples/eigenAeroelastic/BWBflutter.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenAeroelastic/BWBflutter"))
writedlm("test/newTestDataGenerators/eigenAeroelastic/BWBflutter/freqs.txt", test_data_matrix(freqsOsc))
writedlm("test/newTestDataGenerators/eigenAeroelastic/BWBflutter/damps.txt", test_data_matrix(dampsOsc))

# Flutter analysis of the matched wing of the Blended-Wing-Body flying wing
Random.seed!(1234)
include("../examples/eigenAeroelastic/BWBmatchedWingFlutter.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenAeroelastic/BWBmatchedWingFlutter"))
writedlm("test/newTestDataGenerators/eigenAeroelastic/BWBmatchedWingFlutter/freqsOsc.txt", test_data_matrix(freqsOsc))
writedlm("test/newTestDataGenerators/eigenAeroelastic/BWBmatchedWingFlutter/dampsOsc.txt", test_data_matrix(dampsOsc))
writedlm("test/newTestDataGenerators/eigenAeroelastic/BWBmatchedWingFlutter/trimAoA.txt", test_data_matrix(trimAoA))
writedlm("test/newTestDataGenerators/eigenAeroelastic/BWBmatchedWingFlutter/trimThrust.txt", test_data_matrix(trimThrust))
writedlm("test/newTestDataGenerators/eigenAeroelastic/BWBmatchedWingFlutter/trimDelta.txt", test_data_matrix(trimδ))
writedlm("test/newTestDataGenerators/eigenAeroelastic/BWBmatchedWingFlutter/flutterSpeeds.txt", test_data_matrix(flutterSpeeds))

# Flutter analysis of the conventional HALE aircraft in free flight with structural stiffness as the varying parameter
Random.seed!(1234)
include("../examples/eigenAeroelastic/cHALELambdaRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenAeroelastic/cHALELambdaRange"))
writedlm("test/newTestDataGenerators/eigenAeroelastic/cHALELambdaRange/freqs.txt", test_data_matrix(freqsOsc))
writedlm("test/newTestDataGenerators/eigenAeroelastic/cHALELambdaRange/damps.txt", test_data_matrix(dampsOsc))

# Flutter analysis of the conventional HALE aircraft in free flight with airspeed and structural stiffness as the varying parameters
Random.seed!(1234)
include("../examples/eigenAeroelastic/cHALELURange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenAeroelastic/cHALELURange"))
writedlm(string("test/newTestDataGenerators/eigenAeroelastic/cHALELURange/freqs.txt"), test_data_matrix(freqsOsc))
writedlm(string("test/newTestDataGenerators/eigenAeroelastic/cHALELURange/damps.txt"), test_data_matrix(dampsOsc))

# Flutter analysis of the conventional HALE aircraft in free flight with airspeed as the varying parameter
Random.seed!(1234)
include("../examples/eigenAeroelastic/cHALEURange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenAeroelastic/cHALEURange"))
writedlm("test/newTestDataGenerators/eigenAeroelastic/cHALEURange/freqs.txt", test_data_matrix(freqsOsc))
writedlm("test/newTestDataGenerators/eigenAeroelastic/cHALEURange/damps.txt", test_data_matrix(dampsOsc))

# Flutter analysis of the Goland wing
Random.seed!(1234)
include("../examples/eigenAeroelastic/GolandWingFlutter.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenAeroelastic/GolandWingFlutter"))
writedlm("test/newTestDataGenerators/eigenAeroelastic/GolandWingFlutter/freqs.txt", test_data_matrix(freqsOsc))
writedlm("test/newTestDataGenerators/eigenAeroelastic/GolandWingFlutter/damps.txt", test_data_matrix(dampsOsc))

# Flutter analysis of the baseline Healy free FFWT wing with tip loss, root pitch and airspeed as the varying parameters
Random.seed!(1234)
include("../examples/eigenAeroelastic/HealyBaselineFFWTfreeFlutterAoARangeURange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenAeroelastic/HealyBaselineFFWTfreeFlutterAoARangeURange"))
writedlm(string("test/newTestDataGenerators/eigenAeroelastic/HealyBaselineFFWTfreeFlutterAoARangeURange/freqs.txt"), test_data_matrix(freqsOsc))
writedlm(string("test/newTestDataGenerators/eigenAeroelastic/HealyBaselineFFWTfreeFlutterAoARangeURange/damps.txt"), test_data_matrix(dampsOsc))

# Flutter analysis of the baseline Healy free FFWT wing with flare angle and airspeed as the varying parameters
Random.seed!(1234)
include("../examples/eigenAeroelastic/HealyBaselineFFWTfreeFlutterFlareRangeURange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenAeroelastic/HealyBaselineFFWTfreeFlutterFlareRangeURange"))
writedlm(string("test/newTestDataGenerators/eigenAeroelastic/HealyBaselineFFWTfreeFlutterFlareRangeURange/freqs.txt"), test_data_matrix(freqsOsc))
writedlm(string("test/newTestDataGenerators/eigenAeroelastic/HealyBaselineFFWTfreeFlutterFlareRangeURange/damps.txt"), test_data_matrix(dampsOsc))

# Flutter analysis of the baseline Healy locked FFWT wing with airspeed as the varying parameters
Random.seed!(1234)
include("../examples/eigenAeroelastic/HealyBaselineFFWTlockedFlutterURange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenAeroelastic/HealyBaselineFFWTlockedFlutterURange"))
writedlm("test/newTestDataGenerators/eigenAeroelastic/HealyBaselineFFWTlockedFlutterURange/freqs.txt", test_data_matrix(freqsOsc))
writedlm("test/newTestDataGenerators/eigenAeroelastic/HealyBaselineFFWTlockedFlutterURange/damps.txt", test_data_matrix(dampsOsc))

# Flutter analysis of the Healy LCO FFWT wing with root pitch angle as the varying parameter
Random.seed!(1234)
include("../examples/eigenAeroelastic/HealyLCOFFWTfreeFlutterAoARange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenAeroelastic/HealyLCOFFWTfreeFlutterAoARange"))
writedlm(string("test/newTestDataGenerators/eigenAeroelastic/HealyLCOFFWTfreeFlutterAoARange/freqs.txt"), test_data_matrix(freqsOsc))
writedlm(string("test/newTestDataGenerators/eigenAeroelastic/HealyLCOFFWTfreeFlutterAoARange/damps.txt"), test_data_matrix(dampsOsc))

# Flutter analysis of the Healy sideslip FFWT wing with sideslip angle as the varying parameter
Random.seed!(1234)
include("../examples/eigenAeroelastic/HealySideslipFFWTflutterSideslipRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenAeroelastic/HealySideslipFFWTflutterSideslipRange"))
writedlm("test/newTestDataGenerators/eigenAeroelastic/HealySideslipFFWTflutterSideslipRange/freqs.txt", test_data_matrix(freqsOsc))
writedlm("test/newTestDataGenerators/eigenAeroelastic/HealySideslipFFWTflutterSideslipRange/damps.txt", test_data_matrix(dampsOsc))

# Flutter analysis of the Healy sideslip FFWT wing with airspeed as the varying parameter
Random.seed!(1234)
include("../examples/eigenAeroelastic/HealySideslipFFWTflutterURange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenAeroelastic/HealySideslipFFWTflutterURange"))
writedlm("test/newTestDataGenerators/eigenAeroelastic/HealySideslipFFWTflutterURange/freqs.txt", test_data_matrix(freqsOsc))
writedlm("test/newTestDataGenerators/eigenAeroelastic/HealySideslipFFWTflutterURange/damps.txt", test_data_matrix(dampsOsc))

# Flutter analysis of the Helios flying-wing in free flight with payload and structural stiffness as the varying parameters
Random.seed!(1234)
include("../examples/eigenAeroelastic/heliosFlutterPLambdaRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenAeroelastic/heliosFlutterPLambdaRange"))
writedlm(string("test/newTestDataGenerators/eigenAeroelastic/heliosFlutterPLambdaRange/freqs.txt"), test_data_matrix(freqsOsc))
writedlm(string("test/newTestDataGenerators/eigenAeroelastic/heliosFlutterPLambdaRange/damps.txt"), test_data_matrix(dampsOsc))

# Flutter analysis of the Helios flying-wing in free flight with payload as the varying parameter
Random.seed!(1234)
include("../examples/eigenAeroelastic/heliosFlutterPRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenAeroelastic/heliosFlutterPRange"))
writedlm("test/newTestDataGenerators/eigenAeroelastic/heliosFlutterPRange/freqs.txt", test_data_matrix(freqsOsc))
writedlm("test/newTestDataGenerators/eigenAeroelastic/heliosFlutterPRange/damps.txt", test_data_matrix(dampsOsc))

# Flutter analysis of the Helios flying-wing in free flight with airspeed as the varying parameter
Random.seed!(1234)
include("../examples/eigenAeroelastic/heliosFlutterURange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenAeroelastic/heliosFlutterURange"))
writedlm("test/newTestDataGenerators/eigenAeroelastic/heliosFlutterURange/freqs.txt", test_data_matrix(freqsOsc))
writedlm("test/newTestDataGenerators/eigenAeroelastic/heliosFlutterURange/damps.txt", test_data_matrix(dampsOsc))

# Flutter analysis of the wing of the Helios flying-wing
Random.seed!(1234)
include("../examples/eigenAeroelastic/heliosWingFlutter.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenAeroelastic/heliosWingFlutter"))
writedlm("test/newTestDataGenerators/eigenAeroelastic/heliosWingFlutter/freqs.txt", test_data_matrix(freqsOsc))
writedlm("test/newTestDataGenerators/eigenAeroelastic/heliosWingFlutter/damps.txt", test_data_matrix(dampsOsc))

# Flutter analysis of the Pazy FFWT wing with airspeed as the varying parameter
Random.seed!(1234)
include("../examples/eigenAeroelastic/PazyFFWTflutterURange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenAeroelastic/PazyFFWTflutterURange"))
writedlm("test/newTestDataGenerators/eigenAeroelastic/PazyFFWTflutterURange/freqs.txt", test_data_matrix(freqsOsc))
writedlm("test/newTestDataGenerators/eigenAeroelastic/PazyFFWTflutterURange/damps.txt", test_data_matrix(dampsOsc))
writedlm("test/newTestDataGenerators/eigenAeroelastic/PazyFFWTflutterURange/phiHinge.txt", test_data_matrix(ϕHinge))

# Flutter analysis of the Pazy wing
Random.seed!(1234)
include("../examples/eigenAeroelastic/PazyWingFlutter.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenAeroelastic/PazyWingFlutter"))
writedlm("test/newTestDataGenerators/eigenAeroelastic/PazyWingFlutter/flutterOnsetSpeeds.txt", test_data_matrix(flutterOnsetSpeeds))
writedlm("test/newTestDataGenerators/eigenAeroelastic/PazyWingFlutter/flutterOnsetFreqs.txt", test_data_matrix(flutterOnsetFreqs))
writedlm("test/newTestDataGenerators/eigenAeroelastic/PazyWingFlutter/flutterOnsetDisp.txt", test_data_matrix(flutterOnsetDisp))
writedlm("test/newTestDataGenerators/eigenAeroelastic/PazyWingFlutter/flutterOffsetSpeeds.txt", test_data_matrix(flutterOffsetSpeeds))
writedlm("test/newTestDataGenerators/eigenAeroelastic/PazyWingFlutter/flutterOffsetFreqs.txt", test_data_matrix(flutterOffsetFreqs))
writedlm("test/newTestDataGenerators/eigenAeroelastic/PazyWingFlutter/flutterOffsetDisp.txt", test_data_matrix(flutterOffsetDisp))

# Flutter and divergence analysis of the Pazy wing
Random.seed!(1234)
include("../examples/eigenAeroelastic/PazyWingFlutterAndDivergence.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenAeroelastic/PazyWingFlutterAndDivergence"))
writedlm("test/newTestDataGenerators/eigenAeroelastic/PazyWingFlutterAndDivergence/flutterOnsetSpeeds.txt", test_data_matrix(flutterOnsetSpeeds))
writedlm("test/newTestDataGenerators/eigenAeroelastic/PazyWingFlutterAndDivergence/flutterOnsetFreqs.txt", test_data_matrix(flutterOnsetFreqs))

# Flutter analysis of the Pazy wing with varying root pitch angle
Random.seed!(1234)
include("../examples/eigenAeroelastic/PazyWingFlutterPitchRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenAeroelastic/PazyWingFlutterPitchRange"))
writedlm(string("test/newTestDataGenerators/eigenAeroelastic/PazyWingFlutterPitchRange/freqs.txt"), test_data_matrix(freqsOsc))
writedlm(string("test/newTestDataGenerators/eigenAeroelastic/PazyWingFlutterPitchRange/damps.txt"), test_data_matrix(dampsOsc))

# Flutter analysis of the Pazy wing with the wingtip rod modeled as a beam
Random.seed!(1234)
include("../examples/eigenAeroelastic/PazyWingBeamRodFlutter.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenAeroelastic/PazyWingBeamRodFlutter"))
writedlm("test/newTestDataGenerators/eigenAeroelastic/PazyWingBeamRodFlutter/flutterOnsetSpeeds.txt", test_data_matrix(flutterOnsetSpeeds))
writedlm("test/newTestDataGenerators/eigenAeroelastic/PazyWingBeamRodFlutter/flutterOnsetFreqs.txt", test_data_matrix(flutterOnsetFreqs))
writedlm("test/newTestDataGenerators/eigenAeroelastic/PazyWingBeamRodFlutter/flutterOffsetSpeeds.txt", test_data_matrix(flutterOffsetSpeeds))
writedlm("test/newTestDataGenerators/eigenAeroelastic/PazyWingBeamRodFlutter/flutterOffsetFreqs.txt", test_data_matrix(flutterOffsetFreqs))

# Flutter analysis of the Pazy wing without the wingtip rod
Random.seed!(1234)
include("../examples/eigenAeroelastic/PazyWingNoRodFlutter.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenAeroelastic/PazyWingNoRodFlutter"))
keepOnset = flutterOnsetModes .!= 5 # discard the events of the fifth tracked mode, which are intermittent
keepOffset = flutterOffsetModes .!= 5
writedlm("test/newTestDataGenerators/eigenAeroelastic/PazyWingNoRodFlutter/flutterOnsetSpeeds.txt", test_data_matrix(flutterOnsetSpeeds[keepOnset]))
writedlm("test/newTestDataGenerators/eigenAeroelastic/PazyWingNoRodFlutter/flutterOnsetFreqs.txt", test_data_matrix(flutterOnsetFreqs[keepOnset]))
writedlm("test/newTestDataGenerators/eigenAeroelastic/PazyWingNoRodFlutter/flutterOffsetSpeeds.txt", test_data_matrix(flutterOffsetSpeeds[keepOffset]))
writedlm("test/newTestDataGenerators/eigenAeroelastic/PazyWingNoRodFlutter/flutterOffsetFreqs.txt", test_data_matrix(flutterOffsetFreqs[keepOffset]))

# Flutter analysis of the Pazy wing with the wingtip rod modeled as a beam with varying root pitch angle
Random.seed!(1234)
include("../examples/eigenAeroelastic/PazyWingBeamRodFlutterPitchRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenAeroelastic/PazyWingBeamRodFlutterPitchRange"))
writedlm(string("test/newTestDataGenerators/eigenAeroelastic/PazyWingBeamRodFlutterPitchRange/freqs.txt"), test_data_matrix(freqsOsc))
writedlm(string("test/newTestDataGenerators/eigenAeroelastic/PazyWingBeamRodFlutterPitchRange/damps.txt"), test_data_matrix(dampsOsc))

# Flutter analysis of the Pazy wing without the wingtip rod with varying root pitch angle
Random.seed!(1234)
include("../examples/eigenAeroelastic/PazyWingNoRodFlutterPitchRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenAeroelastic/PazyWingNoRodFlutterPitchRange"))
writedlm(string("test/newTestDataGenerators/eigenAeroelastic/PazyWingNoRodFlutterPitchRange/freqs.txt"), test_data_matrix(freqsOsc))
writedlm(string("test/newTestDataGenerators/eigenAeroelastic/PazyWingNoRodFlutterPitchRange/damps.txt"), test_data_matrix(dampsOsc))

# Flutter analysis of the Pazy wing with varying tip mass positions
Random.seed!(1234)
include("../examples/eigenAeroelastic/PazyWingFlutterTipMassRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenAeroelastic/PazyWingFlutterTipMassRange"))
writedlm(string("test/newTestDataGenerators/eigenAeroelastic/PazyWingFlutterTipMassRange/freqs.txt"), test_data_matrix(freqsOsc))
writedlm(string("test/newTestDataGenerators/eigenAeroelastic/PazyWingFlutterTipMassRange/damps.txt"), test_data_matrix(dampsOsc))

# Frequency analysis of the Pazy wing
Random.seed!(1234)
include("../examples/eigenAeroelastic/PazyWingFreqsEvolution.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenAeroelastic/PazyWingFreqsEvolution"))
writedlm("test/newTestDataGenerators/eigenAeroelastic/PazyWingFreqsEvolution/freqs.txt", test_data_matrix(freqsOsc))

# Flutter analysis of the sixteen-meter-wing
Random.seed!(1234)
include("../examples/eigenAeroelastic/SMWFlutter.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenAeroelastic/SMWFlutter"))
writedlm("test/newTestDataGenerators/eigenAeroelastic/SMWFlutter/x1_def.txt", test_data_matrix(x1_def[end]/L))
writedlm("test/newTestDataGenerators/eigenAeroelastic/SMWFlutter/x3_def.txt", test_data_matrix(x3_def[end]/L))
writedlm("test/newTestDataGenerators/eigenAeroelastic/SMWFlutter/alpha.txt", test_data_matrix(α_of_x1[end]))
writedlm(string("test/newTestDataGenerators/eigenAeroelastic/SMWFlutter/freqs.txt"), test_data_matrix(freqsOsc))
writedlm(string("test/newTestDataGenerators/eigenAeroelastic/SMWFlutter/damps.txt"), test_data_matrix(dampsOsc))

# Flutter boundary analysis of the sixteen-meter-wing as a function of the pitch angle
Random.seed!(1234)
include("../examples/eigenAeroelastic/SMWFlutterPitchRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenAeroelastic/SMWFlutterPitchRange"))
writedlm(string("test/newTestDataGenerators/eigenAeroelastic/SMWFlutterPitchRange/freqs.txt"),test_data_matrix(freqsOsc))
writedlm(string("test/newTestDataGenerators/eigenAeroelastic/SMWFlutterPitchRange/damps.txt"),test_data_matrix(dampsOsc))

# Flutter boundary analysis of the sixteen-meter-wing as a function of the wing curvature with tip load as the varying parameter
Random.seed!(1234)
include("../examples/eigenAeroelastic/SMWFlutterPrecurvatureRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenAeroelastic/SMWFlutterPrecurvatureRange"))
writedlm(string("test/newTestDataGenerators/eigenAeroelastic/SMWFlutterPrecurvatureRange/flutterSpeed.txt"), test_data_matrix(flutterSpeed))

# Flutter boundary analysis of the sixteen-meter-wing as a function of the wing curvature with root angle as the varying parameter
Random.seed!(1234)
include("../examples/eigenAeroelastic/SMWFlutterPrecurvatureRange2.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenAeroelastic/SMWFlutterPrecurvatureRange2"))
writedlm(string("test/newTestDataGenerators/eigenAeroelastic/SMWFlutterPrecurvatureRange2/freqs.txt"), test_data_matrix(freqsOsc))
writedlm(string("test/newTestDataGenerators/eigenAeroelastic/SMWFlutterPrecurvatureRange2/damps.txt"), test_data_matrix(dampsOsc))

# Flutter boundary analysis of the sixteen-meter-wing as a function of the bending-torsion coupling factor
Random.seed!(1234)
include("../examples/eigenAeroelastic/SMWFlutterStructuralCouplingRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenAeroelastic/SMWFlutterStructuralCouplingRange"))
writedlm("test/newTestDataGenerators/eigenAeroelastic/SMWFlutterStructuralCouplingRange/flutterSpeed.txt", test_data_matrix(flutterSpeed))

# Flutter boundary analysis of the sixteen-meter-wing as a function of the tip displacement
Random.seed!(1234)
include("../examples/eigenAeroelastic/SMWFlutterTipDispRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenAeroelastic/SMWFlutterTipDispRange"))
writedlm("test/newTestDataGenerators/eigenAeroelastic/SMWFlutterTipDispRange/flutterSpeed.txt", test_data_matrix(flutterSpeed))
writedlm("test/newTestDataGenerators/eigenAeroelastic/SMWFlutterTipDispRange/flutterFreq.txt", test_data_matrix(flutterFreq))

# Linear flutter analysis of the sixteen-meter-wing
Random.seed!(1234)
include("../examples/eigenAeroelastic/SMWLinearFlutter.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenAeroelastic/SMWLinearFlutter"))
writedlm("test/newTestDataGenerators/eigenAeroelastic/SMWLinearFlutter/flutterSpeed.txt", test_data_matrix(flutterSpeed))
writedlm("test/newTestDataGenerators/eigenAeroelastic/SMWLinearFlutter/flutterFreq.txt", test_data_matrix(flutterFreq))

# Torsional divergence analysis of a straight wing
Random.seed!(1234)
include("../examples/eigenAeroelastic/straightWingTorsionalDivergence.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenAeroelastic/straightWingTorsionalDivergence"))
writedlm("test/newTestDataGenerators/eigenAeroelastic/straightWingTorsionalDivergence/damps.txt", test_data_matrix(dampsNonOsc))

# Coupled bending-torsion divergence analysis of a straight wing
Random.seed!(1234)
include("../examples/eigenAeroelastic/straightWingBendingTorsionalDivergence.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenAeroelastic/straightWingBendingTorsionalDivergence"))
writedlm("test/newTestDataGenerators/eigenAeroelastic/straightWingBendingTorsionalDivergence/dampsNonOsc.txt", test_data_matrix(dampsNonOsc))
writedlm("test/newTestDataGenerators/eigenAeroelastic/straightWingBendingTorsionalDivergence/UD.txt", test_data_matrix(UD))

# Bending divergence analysis of a swept wing
Random.seed!(1234)
include("../examples/eigenAeroelastic/sweptWingBendingDivergence.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenAeroelastic/sweptWingBendingDivergence"))
writedlm("test/newTestDataGenerators/eigenAeroelastic/sweptWingBendingDivergence/dampsNonOsc.txt", test_data_matrix(dampsNonOsc))
writedlm("test/newTestDataGenerators/eigenAeroelastic/sweptWingBendingDivergence/UD.txt", test_data_matrix(UD))

# Coupled bending-torsion divergence analysis of a swept wing
Random.seed!(1234)
include("../examples/eigenAeroelastic/sweptWingBendingTorsionalDivergence.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenAeroelastic/sweptWingBendingTorsionalDivergence"))
writedlm("test/newTestDataGenerators/eigenAeroelastic/sweptWingBendingTorsionalDivergence/dampsNonOsc.txt", test_data_matrix(dampsNonOsc))
writedlm("test/newTestDataGenerators/eigenAeroelastic/sweptWingBendingTorsionalDivergence/UD.txt", test_data_matrix(UD))

# Torsional divergence analysis of a swept wing
Random.seed!(1234)
include("../examples/eigenAeroelastic/sweptWingTorsionalDivergence.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenAeroelastic/sweptWingTorsionalDivergence"))
writedlm("test/newTestDataGenerators/eigenAeroelastic/sweptWingTorsionalDivergence/dampsNonOsc.txt", test_data_matrix(dampsNonOsc))
writedlm("test/newTestDataGenerators/eigenAeroelastic/sweptWingTorsionalDivergence/UD.txt", test_data_matrix(UD))

# Flutter analysis of the swept Pazy wing
Random.seed!(1234)
include("../examples/eigenAeroelastic/sweptPazyFlutterPitchRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenAeroelastic/sweptPazyFlutterPitchRange"))
writedlm(string("test/newTestDataGenerators/eigenAeroelastic/sweptPazyFlutterPitchRange/freqs.txt"), test_data_matrix(freqsOsc))
writedlm(string("test/newTestDataGenerators/eigenAeroelastic/sweptPazyFlutterPitchRange/damps.txt"), test_data_matrix(dampsOsc))

# Flutter analysis of the undeformed swept Pazy wing
Random.seed!(1234)
include("../examples/eigenAeroelastic/sweptPazyUndeformedFlutterSweepRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenAeroelastic/sweptPazyUndeformedFlutterSweepRange"))
writedlm(string("test/newTestDataGenerators/eigenAeroelastic/sweptPazyUndeformedFlutterSweepRange/freqs.txt"), test_data_matrix(freqsOsc))
writedlm(string("test/newTestDataGenerators/eigenAeroelastic/sweptPazyUndeformedFlutterSweepRange/damps.txt"), test_data_matrix(dampsOsc))

# Bending divergence analyses of swept wings
Random.seed!(1234)
include("../examples/eigenAeroelastic/sweptWingBendingDivergenceSweepRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenAeroelastic/sweptWingBendingDivergenceSweepRange"))
writedlm(string("test/newTestDataGenerators/eigenAeroelastic/sweptWingBendingDivergenceSweepRange/UD.txt"), test_data_matrix(UD))

# Coupled bending-torsion divergence analyses of swept wings
Random.seed!(1234)
include("../examples/eigenAeroelastic/sweptWingBendingTorsionalDivergenceSweepRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenAeroelastic/sweptWingBendingTorsionalDivergenceSweepRange"))
writedlm(string("test/newTestDataGenerators/eigenAeroelastic/sweptWingBendingTorsionalDivergenceSweepRange/UD.txt"), test_data_matrix(UD))

# Torsional divergence analyses of swept wings
Random.seed!(1234)
include("../examples/eigenAeroelastic/sweptWingTorsionalDivergenceSweepRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenAeroelastic/sweptWingTorsionalDivergenceSweepRange"))
writedlm(string("test/newTestDataGenerators/eigenAeroelastic/sweptWingTorsionalDivergenceSweepRange/UD.txt"), test_data_matrix(UD))

# Eigen aeroelastic analysis of the Tang&Dowell wing at varying airspeed
Random.seed!(1234)
include("../examples/eigenAeroelastic/TDWingAirspeedRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenAeroelastic/TDWingAirspeedRange"))
writedlm("test/newTestDataGenerators/eigenAeroelastic/TDWingAirspeedRange/tip_u3.txt", test_data_matrix(tip_u3))
writedlm("test/newTestDataGenerators/eigenAeroelastic/TDWingAirspeedRange/tip_twist.txt", test_data_matrix(tip_twist))

# Pitch divergence analysis of a typical section
Random.seed!(1234)
include("../examples/eigenAeroelastic/typicalSectionPitchDivergence.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenAeroelastic/typicalSectionPitchDivergence"))
writedlm("test/newTestDataGenerators/eigenAeroelastic/typicalSectionPitchDivergence/damps.txt", test_data_matrix(dampsNonOsc))

# Flutter and divergence analysis of a typical section
Random.seed!(1234)
include("../examples/eigenAeroelastic/typicalSectionFlutterAndDivergence.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenAeroelastic/typicalSectionFlutterAndDivergence"))
writedlm("test/newTestDataGenerators/eigenAeroelastic/typicalSectionFlutterAndDivergence/freqsOsc.txt", test_data_matrix(freqsOsc))
writedlm("test/newTestDataGenerators/eigenAeroelastic/typicalSectionFlutterAndDivergence/dampsOsc.txt", test_data_matrix(dampsOsc))
writedlm("test/newTestDataGenerators/eigenAeroelastic/typicalSectionFlutterAndDivergence/freqsNonOsc.txt", test_data_matrix(freqsNonOsc))
writedlm("test/newTestDataGenerators/eigenAeroelastic/typicalSectionFlutterAndDivergence/dampsNonOsc.txt", test_data_matrix(dampsNonOsc))