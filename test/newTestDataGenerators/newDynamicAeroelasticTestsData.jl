# Dynamic analysis of the Blended-Wing-Body vehicle undergoing a checked pitch maneuver
include("../examples/dynamicAeroelastic/BWBcheckedPitchManeuver.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/dynamicAeroelastic/BWBcheckedPitchManeuver"))
writedlm("test/newTestDataGenerators/dynamicAeroelastic/BWBcheckedPitchManeuver/rootAoA.txt", test_data_matrix(rootAoA))
writedlm("test/newTestDataGenerators/dynamicAeroelastic/BWBcheckedPitchManeuver/Deltau3.txt", test_data_matrix(Δu3))

# Dynamic analysis of the conventional HALE aircraft undergoing a checked pitch maneuver
include("../examples/dynamicAeroelastic/cHALECheckedPitchManeuver.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/dynamicAeroelastic/cHALECheckedPitchManeuver"))
writedlm("test/newTestDataGenerators/dynamicAeroelastic/cHALECheckedPitchManeuver/wingAoA.txt", test_data_matrix(wingAoA))
writedlm("test/newTestDataGenerators/dynamicAeroelastic/cHALECheckedPitchManeuver/Deltau3.txt", test_data_matrix(Δu3))

# Dynamic analysis of the conventional HALE aircraft undergoing a coordinated turn maneuver
include("../examples/dynamicAeroelastic/cHALECheckedRollManeuver.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/dynamicAeroelastic/cHALECheckedRollManeuver"))
writedlm("test/newTestDataGenerators/dynamicAeroelastic/cHALECheckedRollManeuver/wingAoA.txt", test_data_matrix(wingAoA))
writedlm("test/newTestDataGenerators/dynamicAeroelastic/cHALECheckedRollManeuver/Deltau3.txt", test_data_matrix(Δu3))

# Dynamic analysis of the baseline Healy FFWT wing under a series of one-minus-cosine gusts
include("../examples/dynamicAeroelastic/HealyBaselineFFWTOMCGustFloating.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/dynamicAeroelastic/HealyBaselineFFWTOMCGustFloating"))
writedlm(string("test/newTestDataGenerators/dynamicAeroelastic/HealyBaselineFFWTOMCGustFloating/M2root.txt"),test_data_matrix(M2root))
writedlm(string("test/newTestDataGenerators/dynamicAeroelastic/HealyBaselineFFWTOMCGustFloating/fold.txt"),test_data_matrix(ϕ))

# Dynamic analysis of the baseline Healy FFWT wing for a specific one-minus-cosine gust
include("../examples/dynamicAeroelastic/HealyBaselineFFWTOMCGustFloating2.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/dynamicAeroelastic/HealyBaselineFFWTOMCGustFloating2"))
writedlm(string("test/newTestDataGenerators/dynamicAeroelastic/HealyBaselineFFWTOMCGustFloating2/M2root.txt"),test_data_matrix(M2root))
writedlm(string("test/newTestDataGenerators/dynamicAeroelastic/HealyBaselineFFWTOMCGustFloating2/fold.txt"),test_data_matrix(ϕ))

# Dynamic analysis of the Helios flying-wing undergoing a checked pitch maneuver
include("../examples/dynamicAeroelastic/heliosCheckedPitchManeuver.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/dynamicAeroelastic/heliosCheckedPitchManeuver"))
writedlm("test/newTestDataGenerators/dynamicAeroelastic/heliosCheckedPitchManeuver/rootAoA.txt", test_data_matrix(rootAoA))
writedlm("test/newTestDataGenerators/dynamicAeroelastic/heliosCheckedPitchManeuver/Deltau3.txt", test_data_matrix(Δu3))

# Dynamic analysis of the Pazy wing encountering a continuous, 1-dimensional gust defined over time
include("../examples/dynamicAeroelastic/PazyWingContinuous1DGust.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/dynamicAeroelastic/PazyWingContinuous1DGust"))
writedlm("test/newTestDataGenerators/dynamicAeroelastic/PazyWingContinuous1DGust/tipAoA.txt", test_data_matrix(tipAoA))
writedlm("test/newTestDataGenerators/dynamicAeroelastic/PazyWingContinuous1DGust/tipOOP.txt", test_data_matrix(tipOOP))
writedlm("test/newTestDataGenerators/dynamicAeroelastic/PazyWingContinuous1DGust/tqSpan_cn.txt", test_data_matrix(tqSpan_cn))
writedlm("test/newTestDataGenerators/dynamicAeroelastic/PazyWingContinuous1DGust/tqSpan_cm.txt", test_data_matrix(tqSpan_cm))
writedlm("test/newTestDataGenerators/dynamicAeroelastic/PazyWingContinuous1DGust/tqSpan_ct.txt", test_data_matrix(tqSpan_ct))

# Dynamic analysis of the Pazy wing encountering a continuous, 1-dimensional gust defined over space
include("../examples/dynamicAeroelastic/PazyWingContinuous1DSpaceGust.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/dynamicAeroelastic/PazyWingContinuous1DSpaceGust"))
writedlm("test/newTestDataGenerators/dynamicAeroelastic/PazyWingContinuous1DSpaceGust/tipAoA.txt", test_data_matrix(tipAoA))
writedlm("test/newTestDataGenerators/dynamicAeroelastic/PazyWingContinuous1DSpaceGust/tipOOP.txt", test_data_matrix(tipOOP))
writedlm("test/newTestDataGenerators/dynamicAeroelastic/PazyWingContinuous1DSpaceGust/tqSpan_cn.txt", test_data_matrix(tqSpan_cn))
writedlm("test/newTestDataGenerators/dynamicAeroelastic/PazyWingContinuous1DSpaceGust/tqSpan_cm.txt", test_data_matrix(tqSpan_cm))
writedlm("test/newTestDataGenerators/dynamicAeroelastic/PazyWingContinuous1DSpaceGust/tqSpan_ct.txt", test_data_matrix(tqSpan_ct))

# Dynamic analysis of the Pazy wing encountering a continuous, 2-dimensional gust
include("../examples/dynamicAeroelastic/PazyWingContinuous2DSpaceGust.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/dynamicAeroelastic/PazyWingContinuous2DSpaceGust"))
writedlm("test/newTestDataGenerators/dynamicAeroelastic/PazyWingContinuous2DSpaceGust/tipAoA.txt", test_data_matrix(tipAoA))
writedlm("test/newTestDataGenerators/dynamicAeroelastic/PazyWingContinuous2DSpaceGust/tipOOP.txt", test_data_matrix(tipOOP))
writedlm("test/newTestDataGenerators/dynamicAeroelastic/PazyWingContinuous2DSpaceGust/tqSpan_cn.txt", test_data_matrix(tqSpan_cn))
writedlm("test/newTestDataGenerators/dynamicAeroelastic/PazyWingContinuous2DSpaceGust/tqSpan_cm.txt", test_data_matrix(tqSpan_cm))
writedlm("test/newTestDataGenerators/dynamicAeroelastic/PazyWingContinuous2DSpaceGust/tqSpan_ct.txt", test_data_matrix(tqSpan_ct))

# Dynamic analysis of the Pazy wing encountering a DARPA gust
include("../examples/dynamicAeroelastic/PazyWingDARPAGust.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/dynamicAeroelastic/PazyWingDARPAGust"))
writedlm("test/newTestDataGenerators/dynamicAeroelastic/PazyWingDARPAGust/tipAoA.txt", test_data_matrix(tipAoA))
writedlm("test/newTestDataGenerators/dynamicAeroelastic/PazyWingDARPAGust/tipOOP.txt", test_data_matrix(tipOOP))
writedlm("test/newTestDataGenerators/dynamicAeroelastic/PazyWingDARPAGust/tqSpan_cn.txt", test_data_matrix(tqSpan_cn))
writedlm("test/newTestDataGenerators/dynamicAeroelastic/PazyWingDARPAGust/tqSpan_cm.txt", test_data_matrix(tqSpan_cm))
writedlm("test/newTestDataGenerators/dynamicAeroelastic/PazyWingDARPAGust/tqSpan_ct.txt", test_data_matrix(tqSpan_ct))

# Dynamic analysis of the Pazy wing encountering a one-minus-cosine gust
include("../examples/dynamicAeroelastic/PazyWingOMCGust.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/dynamicAeroelastic/PazyWingOMCGust"))
writedlm("test/newTestDataGenerators/dynamicAeroelastic/PazyWingOMCGust/tipAoA.txt", test_data_matrix(tipAoA))
writedlm("test/newTestDataGenerators/dynamicAeroelastic/PazyWingOMCGust/tipOOP.txt", test_data_matrix(tipOOP))
writedlm("test/newTestDataGenerators/dynamicAeroelastic/PazyWingOMCGust/tqSpan_cn.txt", test_data_matrix(tqSpan_cn))
writedlm("test/newTestDataGenerators/dynamicAeroelastic/PazyWingOMCGust/tqSpan_cm.txt", test_data_matrix(tqSpan_cm))
writedlm("test/newTestDataGenerators/dynamicAeroelastic/PazyWingOMCGust/tqSpan_ct.txt", test_data_matrix(tqSpan_ct))

# Dynamic analysis of the Pazy wing with a tip impulse force
include("../examples/dynamicAeroelastic/PazyWingTipImpulse.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/dynamicAeroelastic/PazyWingTipImpulse"))
writedlm("test/newTestDataGenerators/dynamicAeroelastic/PazyWingTipImpulse/tipAoA.txt", test_data_matrix(tipAoA))
writedlm("test/newTestDataGenerators/dynamicAeroelastic/PazyWingTipImpulse/tipOOP.txt", test_data_matrix(tipOOP))
writedlm("test/newTestDataGenerators/dynamicAeroelastic/PazyWingTipImpulse/tqSpan_cn.txt", test_data_matrix(tqSpan_cn))
writedlm("test/newTestDataGenerators/dynamicAeroelastic/PazyWingTipImpulse/tqSpan_cm.txt", test_data_matrix(tqSpan_cm))
writedlm("test/newTestDataGenerators/dynamicAeroelastic/PazyWingTipImpulse/tqSpan_ct.txt", test_data_matrix(tqSpan_ct))

# Dynamic analysis of the baseline Healy locked FFWT wing under a series of one-minus-cosine gusts
include("../examples/dynamicAeroelastic/HealyBaselineFFWTlockedOMCGustFloating.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/dynamicAeroelastic/HealyBaselineFFWTlockedOMCGustFloating"))
writedlm(string("test/newTestDataGenerators/dynamicAeroelastic/HealyBaselineFFWTlockedOMCGustFloating/M2root.txt"),test_data_matrix(M2root))

# Dynamic analysis of the baseline Healy free and locked FFWT wings under one-minus-cosine gusts
include("../examples/dynamicAeroelastic/HealyBaselineFFWTlockedVfreeOMCGustFloating.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/dynamicAeroelastic/HealyBaselineFFWTlockedVfreeOMCGustFloating"))
writedlm(string("test/newTestDataGenerators/dynamicAeroelastic/HealyBaselineFFWTlockedVfreeOMCGustFloating/dM2peaks.txt"),test_data_matrix(ΔM2peaks))

# Dynamic analysis of the Healy sideslip FFWT wing under a one-minus-cosine gust
include("../examples/dynamicAeroelastic/HealySideslipFFWTOMCGustFloating.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/dynamicAeroelastic/HealySideslipFFWTOMCGustFloating"))
writedlm("test/newTestDataGenerators/dynamicAeroelastic/HealySideslipFFWTOMCGustFloating/tipOOP.txt", test_data_matrix(tipOOP))
writedlm("test/newTestDataGenerators/dynamicAeroelastic/HealySideslipFFWTOMCGustFloating/M2_root.txt", test_data_matrix(M2_root))
writedlm("test/newTestDataGenerators/dynamicAeroelastic/HealySideslipFFWTOMCGustFloating/root_alphae.txt", test_data_matrix(root_αₑ))
writedlm("test/newTestDataGenerators/dynamicAeroelastic/HealySideslipFFWTOMCGustFloating/root_cn.txt", test_data_matrix(root_cn))

# Dynamic analysis of the Pazy FFWT wing under a one-minus-cosine gust
include("../examples/dynamicAeroelastic/PazyFFWTOMCGustFloating.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/dynamicAeroelastic/PazyFFWTOMCGustFloating"))
writedlm("test/newTestDataGenerators/dynamicAeroelastic/PazyFFWTOMCGustFloating/fold.txt", test_data_matrix(ϕ))
writedlm("test/newTestDataGenerators/dynamicAeroelastic/PazyFFWTOMCGustFloating/tipOOP.txt", test_data_matrix(tipOOP))
writedlm("test/newTestDataGenerators/dynamicAeroelastic/PazyFFWTOMCGustFloating/M2_root.txt", test_data_matrix(M2_root))
writedlm("test/newTestDataGenerators/dynamicAeroelastic/PazyFFWTOMCGustFloating/tqSpan_cn.txt", test_data_matrix(tqSpan_cn))

# Dynamic analysis of a rotating wing using the dynamic stall model
include("../examples/dynamicAeroelastic/rotorDStest.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/dynamicAeroelastic/rotorDStest"))
writedlm("test/newTestDataGenerators/dynamicAeroelastic/rotorDStest/cl.txt", test_data_matrix(cl))
writedlm("test/newTestDataGenerators/dynamicAeroelastic/rotorDStest/cm.txt", test_data_matrix(cm))
writedlm("test/newTestDataGenerators/dynamicAeroelastic/rotorDStest/cdrag.txt", test_data_matrix(cdrag))

# Dynamic analysis of a pitching wing with dynamic stall models
include("../examples/dynamicAeroelastic/wingDStest.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/dynamicAeroelastic/wingDStest"))
writedlm("test/newTestDataGenerators/dynamicAeroelastic/wingDStest/cl.txt", test_data_matrix(cl))
writedlm("test/newTestDataGenerators/dynamicAeroelastic/wingDStest/cm.txt", test_data_matrix(cm))
writedlm("test/newTestDataGenerators/dynamicAeroelastic/wingDStest/cdrag.txt", test_data_matrix(cdrag))