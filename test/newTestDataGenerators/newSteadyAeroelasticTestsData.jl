# Steady aeroelastic analysis of the clamped conventional HALE
include("../examples/steadyAeroelastic/cHALEclampedSteady.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/steadyAeroelastic/cHALEclampedSteady"))
writedlm("test/newTestDataGenerators/steadyAeroelastic/cHALEclampedSteady/x1_def.txt", test_data_matrix(x1_def))
writedlm("test/newTestDataGenerators/steadyAeroelastic/cHALEclampedSteady/x3_def.txt", test_data_matrix(x3_def))

# Steady analysis of the baseline Healy free FFWT wing with varying root pitch and airspeed
include("../examples/steadyAeroelastic/HealyBaselineFFWTsteadyAoARangeURangeCoast.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/steadyAeroelastic/HealyBaselineFFWTsteadyAoARangeURangeCoast"))
writedlm("test/newTestDataGenerators/steadyAeroelastic/HealyBaselineFFWTsteadyAoARangeURangeCoast/phiHinge.txt", test_data_matrix(ϕHinge))
writedlm("test/newTestDataGenerators/steadyAeroelastic/HealyBaselineFFWTsteadyAoARangeURangeCoast/u3Hinge.txt", test_data_matrix(u3Hinge))
writedlm("test/newTestDataGenerators/steadyAeroelastic/HealyBaselineFFWTsteadyAoARangeURangeCoast/M2root.txt", test_data_matrix(M2root))

# Steady analysis of the Healy sideslip FFWT wing with varying flare angle and root pitch angle
include("../examples/steadyAeroelastic/HealySideslipFFWTsteadyFlareRangeAoARangeCoast.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/steadyAeroelastic/HealySideslipFFWTsteadyFlareRangeAoARangeCoast"))
writedlm("test/newTestDataGenerators/steadyAeroelastic/HealySideslipFFWTsteadyFlareRangeAoARangeCoast/phiHinge.txt", test_data_matrix(ϕHinge))

# Steady analysis of the Healy sideslip FFWT wing with varying flare angle, root pitch angle and airspeed
include("../examples/steadyAeroelastic/HealySideslipFFWTsteadyFlareRangeURangeAoARangeCoast.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/steadyAeroelastic/HealySideslipFFWTsteadyFlareRangeURangeAoARangeCoast"))
writedlm("test/newTestDataGenerators/steadyAeroelastic/HealySideslipFFWTsteadyFlareRangeURangeAoARangeCoast/phiHinge.txt", test_data_matrix(ϕHinge))

# Steady analysis of the Healy sideslip FFWT wing with varying wingtip twist, root pitch angle and sideslip angle
include("../examples/steadyAeroelastic/HealySideslipFFWTsteadyTwistRangeAoARangeSideslipRangeCoast.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/steadyAeroelastic/HealySideslipFFWTsteadyTwistRangeAoARangeSideslipRangeCoast"))
writedlm("test/newTestDataGenerators/steadyAeroelastic/HealySideslipFFWTsteadyTwistRangeAoARangeSideslipRangeCoast/phiHinge.txt", test_data_matrix(ϕHinge))

# Steady analysis of the flat-plate wing with a coasting FFWT
include("../examples/steadyAeroelastic/FlatPlateFFWTsteadyCoast.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/steadyAeroelastic/FlatPlateFFWTsteadyCoast"))
writedlm("test/newTestDataGenerators/steadyAeroelastic/FlatPlateFFWTsteadyCoast/phiHinge.txt", test_data_matrix(ϕHinge))

# Steady analysis of the flat-plate left wing with a coasting FFWT
include("../examples/steadyAeroelastic/FlatPlateFFWTsteadyCoastLeft.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/steadyAeroelastic/FlatPlateFFWTsteadyCoastLeft"))
writedlm("test/newTestDataGenerators/steadyAeroelastic/FlatPlateFFWTsteadyCoastLeft/phiHinge.txt", test_data_matrix(ϕHinge))

# Steady analysis of the Pazy wing with a coasting FFWT
include("../examples/steadyAeroelastic/PazyFFWTsteadyCoast.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/steadyAeroelastic/PazyFFWTsteadyCoast"))
writedlm("test/newTestDataGenerators/steadyAeroelastic/PazyFFWTsteadyCoast/phiHinge.txt", test_data_matrix(ϕHinge))

# Steady analysis of the Pazy wing with a FFWT at a fixed fold angle
include("../examples/steadyAeroelastic/PazyFFWTsteadyFixedFold.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/steadyAeroelastic/PazyFFWTsteadyFixedFold"))
writedlm("test/newTestDataGenerators/steadyAeroelastic/PazyFFWTsteadyFixedFold/phiHinge.txt", test_data_matrix(ϕHinge))
writedlm("test/newTestDataGenerators/steadyAeroelastic/PazyFFWTsteadyFixedFold/hingeMoment.txt", test_data_matrix(hingeMoment))

# Steady analysis of the Pazy wing with a coasting FFWT, at varying airspeed and root pitch angle
include("../examples/steadyAeroelastic/PazyFFWTsteadyURangeAoARangeCoast.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/steadyAeroelastic/PazyFFWTsteadyURangeAoARangeCoast"))
writedlm("test/newTestDataGenerators/steadyAeroelastic/PazyFFWTsteadyURangeAoARangeCoast/phiHinge.txt", test_data_matrix(ϕHinge))

# Steady analysis of the Pazy wing with a coasting FFWT, at varying airspeed
include("../examples/steadyAeroelastic/PazyFFWTsteadyURangeCoast.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/steadyAeroelastic/PazyFFWTsteadyURangeCoast"))
writedlm("test/newTestDataGenerators/steadyAeroelastic/PazyFFWTsteadyURangeCoast/phiHinge.txt", test_data_matrix(ϕHinge))

# Steady analysis of the Pazy wing with a FFWT at a fixed fold angle, at varying airspeed
include("../examples/steadyAeroelastic/PazyFFWTsteadyURangeFixedFold.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/steadyAeroelastic/PazyFFWTsteadyURangeFixedFold"))
writedlm("test/newTestDataGenerators/steadyAeroelastic/PazyFFWTsteadyURangeFixedFold/phiHinge.txt", test_data_matrix(ϕHinge))
writedlm("test/newTestDataGenerators/steadyAeroelastic/PazyFFWTsteadyURangeFixedFold/hingeMoment.txt", test_data_matrix(hingeMoment))

# Steady analysis of the Pazy wing with varying root pitch angle
include("../examples/steadyAeroelastic/PazyWingPitchRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/steadyAeroelastic/PazyWingPitchRange"))
writedlm("test/newTestDataGenerators/steadyAeroelastic/PazyWingPitchRange/tip_AoA.txt", test_data_matrix(tip_AoA))
writedlm("test/newTestDataGenerators/steadyAeroelastic/PazyWingPitchRange/tip_OOP.txt", test_data_matrix(tip_OOP))
writedlm("test/newTestDataGenerators/steadyAeroelastic/PazyWingPitchRange/tip_IP.txt", test_data_matrix(tip_IP))
writedlm("test/newTestDataGenerators/steadyAeroelastic/PazyWingPitchRange/tip_twist.txt", test_data_matrix(tip_twist))

# Steady aeroelastic analysis of the sixteen-meter-wing
include("../examples/steadyAeroelastic/SMWSteady.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/steadyAeroelastic/SMWSteady"))
writedlm("test/newTestDataGenerators/steadyAeroelastic/SMWSteady/tip_u3.txt", test_data_matrix(tip_u3))
writedlm("test/newTestDataGenerators/steadyAeroelastic/SMWSteady/tip_twist.txt", test_data_matrix(tip_twist))

# Steady analysis of the Healy LCO FFWT wing with varying airspeed and root pitch angle
include("../examples/steadyAeroelastic/HealyLCOFFWTsteadyURangeAoARangeCoast.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/steadyAeroelastic/HealyLCOFFWTsteadyURangeAoARangeCoast"))
writedlm("test/newTestDataGenerators/steadyAeroelastic/HealyLCOFFWTsteadyURangeAoARangeCoast/phiHinge.txt", test_data_matrix(ϕHinge))

# Creating the conventional HALE model
include("../examples/steadyAeroelastic/cHALEmodel.jl")

# Steady analysis of the swept Pazy wing with varying root pitch angle
include("../examples/steadyAeroelastic/sweptPazySteadyPitchRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/steadyAeroelastic/sweptPazySteadyPitchRange"))
writedlm("test/newTestDataGenerators/steadyAeroelastic/sweptPazySteadyPitchRange/tipOOP.txt", test_data_matrix(tipOOP))
writedlm("test/newTestDataGenerators/steadyAeroelastic/sweptPazySteadyPitchRange/tipTwist.txt", test_data_matrix(tipTwist))
writedlm("test/newTestDataGenerators/steadyAeroelastic/sweptPazySteadyPitchRange/tipAoA.txt", test_data_matrix(tipAoA))

# Steady analysis of the swept Pazy wing with varying sweep and root pitch angles
include("../examples/steadyAeroelastic/sweptPazySteadySweepRangePitchRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/steadyAeroelastic/sweptPazySteadySweepRangePitchRange"))
writedlm("test/newTestDataGenerators/steadyAeroelastic/sweptPazySteadySweepRangePitchRange/tipOOP.txt", test_data_matrix(tipOOP))
writedlm("test/newTestDataGenerators/steadyAeroelastic/sweptPazySteadySweepRangePitchRange/tipTwist.txt", test_data_matrix(tipTwist))

# Steady lift distribution of swept wings
include("../examples/steadyAeroelastic/sweptWingsLiftDistribution.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/steadyAeroelastic/sweptWingsLiftDistribution"))
writedlm("test/newTestDataGenerators/steadyAeroelastic/sweptWingsLiftDistribution/cn.txt", test_data_matrix(cn))