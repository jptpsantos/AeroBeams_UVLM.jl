# Modal analysis of the axial vibration of a beam under clamped-clamped boundary conditions
include("../examples/eigenStructural/beamAxialVibrationCC.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenStructural/beamAxialVibrationCC"))
writedlm("test/newTestDataGenerators/eigenStructural/beamAxialVibrationCC/freqs.txt", test_data_matrix(freqs))
writedlm("test/newTestDataGenerators/eigenStructural/beamAxialVibrationCC/u1_modeShapes.txt", test_data_matrix(u1_modeShapes))

# Modal analysis of the axial vibration of a beam under clamped-free boundary conditions
include("../examples/eigenStructural/beamAxialVibrationCF.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenStructural/beamAxialVibrationCF"))
writedlm("test/newTestDataGenerators/eigenStructural/beamAxialVibrationCF/freqs.txt", test_data_matrix(freqs))
writedlm("test/newTestDataGenerators/eigenStructural/beamAxialVibrationCF/u1_modeShapes.txt", test_data_matrix(u1_modeShapes))

# Modal analysis of the axial vibration of a beam under free-free boundary conditions
include("../examples/eigenStructural/beamAxialVibrationFF.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenStructural/beamAxialVibrationFF"))
writedlm("test/newTestDataGenerators/eigenStructural/beamAxialVibrationFF/freqs.txt", test_data_matrix(freqs))
writedlm("test/newTestDataGenerators/eigenStructural/beamAxialVibrationFF/u1_modeShapes.txt", test_data_matrix(u1_modeShapes))

# Modal analysis of the bending vibration of a beam under clamped-clamped boundary conditions
include("../examples/eigenStructural/beamBendingVibrationCC.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenStructural/beamBendingVibrationCC"))
writedlm("test/newTestDataGenerators/eigenStructural/beamBendingVibrationCC/freqs.txt", test_data_matrix(freqs))
writedlm("test/newTestDataGenerators/eigenStructural/beamBendingVibrationCC/u3_modeShapes.txt", test_data_matrix(u3_modeShapes))

# Modal analysis of the bending vibration of a beam under clamped-free boundary conditions
include("../examples/eigenStructural/beamBendingVibrationCF.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenStructural/beamBendingVibrationCF"))
writedlm("test/newTestDataGenerators/eigenStructural/beamBendingVibrationCF/freqs.txt", test_data_matrix(freqs))
writedlm("test/newTestDataGenerators/eigenStructural/beamBendingVibrationCF/u3_modeShapes.txt", test_data_matrix(u3_modeShapes))

# Modal analysis of the bending vibration of a beam under clamped-pinned boundary conditions
include("../examples/eigenStructural/beamBendingVibrationCP.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenStructural/beamBendingVibrationCP"))
writedlm("test/newTestDataGenerators/eigenStructural/beamBendingVibrationCP/freqs.txt", test_data_matrix(freqs))
writedlm("test/newTestDataGenerators/eigenStructural/beamBendingVibrationCP/u3_modeShapes.txt", test_data_matrix(u3_modeShapes))

# Modal analysis of the bending vibration of a beam under clamped-sliding boundary conditions
include("../examples/eigenStructural/beamBendingVibrationCS.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenStructural/beamBendingVibrationCS"))
writedlm("test/newTestDataGenerators/eigenStructural/beamBendingVibrationCS/freqs.txt", test_data_matrix(freqs))
writedlm("test/newTestDataGenerators/eigenStructural/beamBendingVibrationCS/u3_modeShapes.txt", test_data_matrix(u3_modeShapes))

# Modal analysis of the bending vibration of a beam under free-free boundary conditions
include("../examples/eigenStructural/beamBendingVibrationFF.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenStructural/beamBendingVibrationFF"))
writedlm("test/newTestDataGenerators/eigenStructural/beamBendingVibrationFF/freqs.txt", test_data_matrix(freqs))
writedlm("test/newTestDataGenerators/eigenStructural/beamBendingVibrationFF/u3_modeShapes.txt", test_data_matrix(u3_modeShapes))

# Modal analysis of the bending vibration of a beam under pinned-pinned boundary conditions
include("../examples/eigenStructural/beamBendingVibrationPP.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenStructural/beamBendingVibrationPP"))
writedlm("test/newTestDataGenerators/eigenStructural/beamBendingVibrationPP/freqs.txt", test_data_matrix(freqs))
writedlm("test/newTestDataGenerators/eigenStructural/beamBendingVibrationPP/u3_modeShapes.txt", test_data_matrix(u3_modeShapes))

# Modal analysis of the torsional vibration of a beam under clamped-clamped boundary conditions
include("../examples/eigenStructural/beamTorsionalVibrationCC.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenStructural/beamTorsionalVibrationCC"))
writedlm("test/newTestDataGenerators/eigenStructural/beamTorsionalVibrationCC/freqs.txt", test_data_matrix(freqs))
writedlm("test/newTestDataGenerators/eigenStructural/beamTorsionalVibrationCC/p1_modeShapes.txt", test_data_matrix(p1_modeShapes))

# Modal analysis of the torsional vibration of a beam under clamped-free boundary conditions
include("../examples/eigenStructural/beamTorsionalVibrationCF.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenStructural/beamTorsionalVibrationCF"))
writedlm("test/newTestDataGenerators/eigenStructural/beamTorsionalVibrationCF/freqs.txt", test_data_matrix(freqs))
writedlm("test/newTestDataGenerators/eigenStructural/beamTorsionalVibrationCF/p1_modeShapes.txt", test_data_matrix(p1_modeShapes))

# Modal analysis of the torsional vibration of a beam under free-free boundary conditions
include("../examples/eigenStructural/beamTorsionalVibrationFF.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenStructural/beamTorsionalVibrationFF"))
writedlm("test/newTestDataGenerators/eigenStructural/beamTorsionalVibrationFF/freqs.txt", test_data_matrix(freqs))
writedlm("test/newTestDataGenerators/eigenStructural/beamTorsionalVibrationFF/p1_modeShapes.txt", test_data_matrix(p1_modeShapes))

# Modal analysis of a cantilever with stiffness-proportional structural damping
include("../examples/eigenStructural/cantileverStructuralDampingEigen.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenStructural/cantileverStructuralDampingEigen"))
writedlm("test/newTestDataGenerators/eigenStructural/cantileverStructuralDampingEigen/freqs.txt", test_data_matrix(ωNum))
writedlm("test/newTestDataGenerators/eigenStructural/cantileverStructuralDampingEigen/damps.txt", test_data_matrix(ζNum))

# Modal analysis of a cantilever beam with a tip axial inertia
include("../examples/eigenStructural/cantileverWithTipAxialMassEigen.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenStructural/cantileverWithTipAxialMassEigen"))
writedlm("test/newTestDataGenerators/eigenStructural/cantileverWithTipAxialMassEigen/freqsNorm.txt", test_data_matrix(freqsNorm))
writedlm("test/newTestDataGenerators/eigenStructural/cantileverWithTipAxialMassEigen/u1_modeShapes.txt", test_data_matrix(u1_modeShapes))

# Modal analysis of a cantilever beam with a tip axial spring
include("../examples/eigenStructural/cantileverWithTipAxialSpringEigen.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenStructural/cantileverWithTipAxialSpringEigen"))
writedlm("test/newTestDataGenerators/eigenStructural/cantileverWithTipAxialSpringEigen/freqsNorm.txt", test_data_matrix(freqsNorm))
writedlm("test/newTestDataGenerators/eigenStructural/cantileverWithTipAxialSpringEigen/u1_modeShapes.txt", test_data_matrix(u1_modeShapes))

# Modal analysis of a cantilever beam with a tip spring in bending
include("../examples/eigenStructural/cantileverWithTipSpringEigen.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenStructural/cantileverWithTipSpringEigen"))
writedlm("test/newTestDataGenerators/eigenStructural/cantileverWithTipSpringEigen/freqsNorm.txt", test_data_matrix(freqsNorm))
writedlm("test/newTestDataGenerators/eigenStructural/cantileverWithTipSpringEigen/u3_modeShapes.txt", test_data_matrix(u3_modeShapes))

# Modal analysis of a cantilever beam with a tip torsional inertia
include("../examples/eigenStructural/cantileverWithTipTorsionalInertiaEigen.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenStructural/cantileverWithTipTorsionalInertiaEigen"))
writedlm("test/newTestDataGenerators/eigenStructural/cantileverWithTipTorsionalInertiaEigen/freqsNorm.txt", test_data_matrix(freqsNorm))
writedlm("test/newTestDataGenerators/eigenStructural/cantileverWithTipTorsionalInertiaEigen/p1_modeShapes.txt", test_data_matrix(p1_modeShapes))

# Modal analysis of a cantilever beam with a tip torsional spring
include("../examples/eigenStructural/cantileverWithTipTorsionalSpringEigen.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenStructural/cantileverWithTipTorsionalSpringEigen"))
writedlm("test/newTestDataGenerators/eigenStructural/cantileverWithTipTorsionalSpringEigen/freqsNorm.txt", test_data_matrix(freqsNorm))
writedlm("test/newTestDataGenerators/eigenStructural/cantileverWithTipTorsionalSpringEigen/p1_modeShapes.txt", test_data_matrix(p1_modeShapes))

# Modal analysis a beam clamped at one end, simply-supported at the other and with a tip inertia
include("../examples/eigenStructural/clampedSSBeamWIthTipInertiaEigen.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenStructural/clampedSSBeamWIthTipInertiaEigen"))
writedlm("test/newTestDataGenerators/eigenStructural/clampedSSBeamWIthTipInertiaEigen/freqs.txt", test_data_matrix(freqs))
writedlm("test/newTestDataGenerators/eigenStructural/clampedSSBeamWIthTipInertiaEigen/u3_modeShapes.txt", test_data_matrix(u3_modeShapes))

# Modal analysis of the baseline Healy free FFWT wing without gravity
include("../examples/eigenStructural/HealyBaselineFFWTModalFree.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenStructural/HealyBaselineFFWTModalFree"))
writedlm("test/newTestDataGenerators/eigenStructural/HealyBaselineFFWTModalFree/freqs.txt", test_data_matrix(freqs))

# Modal analysis of the baseline Healy locked FFWT wing
include("../examples/eigenStructural/HealyLCOFFWTLockedModal.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenStructural/HealyLCOFFWTLockedModal"))
writedlm("test/newTestDataGenerators/eigenStructural/HealyLCOFFWTLockedModal/freqs.txt", test_data_matrix(freqs))

# Modal analysis of the baseline Healy FFWT wing with wingtip removed
include("../examples/eigenStructural/HealyLCOFFWTRemovedModal.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenStructural/HealyLCOFFWTRemovedModal"))
writedlm("test/newTestDataGenerators/eigenStructural/HealyLCOFFWTRemovedModal/freqs.txt", test_data_matrix(freqs))

# Modal analyses of the Pazy wing in horizontal and vertical positions
include("../examples/eigenStructural/PazyWingModal.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenStructural/PazyWingModal"))
writedlm("test/newTestDataGenerators/eigenStructural/PazyWingModal/freqs.txt", test_data_matrix(freqs))

# Modal analysis of the Pazy wing with the wingtip rod modeled as a beam
include("../examples/eigenStructural/PazyWingBeamRodModal.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenStructural/PazyWingBeamRodModal"))
writedlm("test/newTestDataGenerators/eigenStructural/PazyWingBeamRodModal/freqs.txt", test_data_matrix(freqs))

# Modal analysis of the Pazy wing without the wingtip rod
include("../examples/eigenStructural/PazyWingNoRodModal.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenStructural/PazyWingNoRodModal"))
writedlm("test/newTestDataGenerators/eigenStructural/PazyWingNoRodModal/freqs.txt", test_data_matrix(freqs))

# Modal analysis of the Pazy wing with 10 degrees of sweep and the wingtip rod modeled as a beam
include("../examples/eigenStructural/PazyS10BeamRodModal.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenStructural/PazyS10BeamRodModal"))
writedlm("test/newTestDataGenerators/eigenStructural/PazyS10BeamRodModal/freqs.txt", test_data_matrix(freqs))

# Modal analysis of a beam pinned at one end and transversely springed at the other
include("../examples/eigenStructural/pinnedSpringedBeamEigen.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenStructural/pinnedSpringedBeamEigen"))
writedlm("test/newTestDataGenerators/eigenStructural/pinnedSpringedBeamEigen/freqs.txt", test_data_matrix(freqs))

# Modal analysis of the sixteen-meter-wing
include("../examples/eigenStructural/SMWModal.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenStructural/SMWModal"))
writedlm("test/newTestDataGenerators/eigenStructural/SMWModal/freqs.txt", test_data_matrix(freqs))

# Modal analysis of a straight rotor under varying angular velocities
include("../examples/eigenStructural/straightRotor.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenStructural/straightRotor"))
writedlm("test/newTestDataGenerators/eigenStructural/straightRotor/freqs.txt", test_data_matrix(freqs))

# Modal analysis of the undeformed swept Pazy wing
include("../examples/eigenStructural/sweptPazyModal.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenStructural/sweptPazyModal"))
writedlm("test/newTestDataGenerators/eigenStructural/sweptPazyModal/freqs.txt", test_data_matrix(freqs))

# Modal analysis of the swept Pazy wing with varying tip displacement
include("../examples/eigenStructural/sweptPazyModalTipDispRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenStructural/sweptPazyModalTipDispRange"))
writedlm("test/newTestDataGenerators/eigenStructural/sweptPazyModalTipDispRange/tipOOP.txt", test_data_matrix(tipOOP))
writedlm("test/newTestDataGenerators/eigenStructural/sweptPazyModalTipDispRange/tip_twist.txt", test_data_matrix(tip_twist))
writedlm("test/newTestDataGenerators/eigenStructural/sweptPazyModalTipDispRange/freqs.txt", test_data_matrix(freqsOsc))

# Modal analysis of a swept-tip rotor under varying angular velocities
include("../examples/eigenStructural/sweptTipRotor.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenStructural/sweptTipRotor"))
writedlm("test/newTestDataGenerators/eigenStructural/sweptTipRotor/freqs.txt", test_data_matrix(freqs))

# Modal analysis of a cantilevered tapered beam
include("../examples/eigenStructural/taperedBeamEigen.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenStructural/taperedBeamEigen"))
writedlm("test/newTestDataGenerators/eigenStructural/taperedBeamEigen/freqs.txt", test_data_matrix(freqs))

# Modal analysis of the Tang&Dowell wing at varying pitch angles
include("../examples/eigenStructural/TDWingPitchRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenStructural/TDWingPitchRange"))
writedlm("test/newTestDataGenerators/eigenStructural/TDWingPitchRange/tip_u3.txt", test_data_matrix(tip_u3))
writedlm("test/newTestDataGenerators/eigenStructural/TDWingPitchRange/tip_u2.txt", test_data_matrix(tip_u2))
writedlm("test/newTestDataGenerators/eigenStructural/TDWingPitchRange/tip_twist.txt", test_data_matrix(tip_twist))
writedlm("test/newTestDataGenerators/eigenStructural/TDWingPitchRange/freqs.txt", test_data_matrix(freqs))

# Modal analysis of a 2-story frame
include("../examples/eigenStructural/twoStoryFrame.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/eigenStructural/twoStoryFrame"))
writedlm("test/newTestDataGenerators/eigenStructural/twoStoryFrame/freqs.txt", test_data_matrix(freqs))