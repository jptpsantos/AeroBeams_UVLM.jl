# Dynamic analysis of an airfoils with harmonic flap deflection profile
include("../examples/aerodynamic/flapOscillation.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/aerodynamic/flapOscillation"))
writedlm("test/newTestDataGenerators/aerodynamic/flapOscillation/cn.txt", test_data_matrix(cn))
writedlm("test/newTestDataGenerators/aerodynamic/flapOscillation/cm.txt", test_data_matrix(cm))

# Dynamic analysis of two airfoils with linked harmonic flap deflection profiles
include("../examples/aerodynamic/flapOscillationLinked.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/aerodynamic/flapOscillationLinked"))
writedlm("test/newTestDataGenerators/aerodynamic/flapOscillationLinked/cnMaster.txt", test_data_matrix(cnMaster))
writedlm("test/newTestDataGenerators/aerodynamic/flapOscillationLinked/cmMaster.txt", test_data_matrix(cmMaster))
writedlm("test/newTestDataGenerators/aerodynamic/flapOscillationLinked/cnSlave.txt", test_data_matrix(cnSlave))
writedlm("test/newTestDataGenerators/aerodynamic/flapOscillationLinked/cmSlave.txt", test_data_matrix(cmSlave))

# One-minus-cosine gust response of an airfoil section at several pitch angles
include("../examples/aerodynamic/OMCgustTests.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/aerodynamic/OMCgustTests"))
writedlm("test/newTestDataGenerators/aerodynamic/OMCgustTests/dcl.txt", test_data_matrix(Δcl))

# One-minus-cosine gust response of an airfoil section at several pitch angles
include("../examples/aerodynamic/OMCgustTests2.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/aerodynamic/OMCgustTests2"))
writedlm("test/newTestDataGenerators/aerodynamic/OMCgustTests2/dcl.txt", test_data_matrix(Δcl))

# Dynamic analysis of a harmonically pitching airfoil, using the dynamic stall model
include("../examples/aerodynamic/pitchingAirfoilDSModelTest.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/aerodynamic/pitchingAirfoilDSModelTest"))
writedlm("test/newTestDataGenerators/aerodynamic/pitchingAirfoilDSModelTest/"*frameString*"_cn.txt", test_data_matrix(cn))
writedlm("test/newTestDataGenerators/aerodynamic/pitchingAirfoilDSModelTest/"*frameString*"_cm.txt", test_data_matrix(cm))
writedlm("test/newTestDataGenerators/aerodynamic/pitchingAirfoilDSModelTest/"*frameString*"_ct.txt", test_data_matrix(ct))

# Dynamic analysis of a harmonically plunging airfoil, using the dynamic stall model
include("../examples/aerodynamic/plungingAirfoilDSModelTest.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/aerodynamic/plungingAirfoilDSModelTest"))
writedlm("test/newTestDataGenerators/aerodynamic/plungingAirfoilDSModelTest/cn.txt", test_data_matrix(cn))
writedlm("test/newTestDataGenerators/aerodynamic/plungingAirfoilDSModelTest/cm.txt", test_data_matrix(cm))

# Sharp-edged gust response of an airfoil section at several pitch angles
include("../examples/aerodynamic/SEgustTests.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/aerodynamic/SEgustTests"))
writedlm("test/newTestDataGenerators/aerodynamic/SEgustTests/dcl.txt", test_data_matrix(Δcl))

# Dynamic analysis of an airfoil section sinusoidally surging (facing a time-varying freestream)
include("../examples/aerodynamic/timeVaryingFreestream.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/aerodynamic/timeVaryingFreestream"))
writedlm("test/newTestDataGenerators/aerodynamic/timeVaryingFreestream/cn.txt", test_data_matrix(cn[end]))
writedlm("test/newTestDataGenerators/aerodynamic/timeVaryingFreestream/cm.txt", test_data_matrix(cm[end]))

# Dynamic analysis of an airfoil section sinusoidally pitching and surging (facing a time-varying freestream)
include("../examples/aerodynamic/timeVaryingFreestreamAndPitch.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/aerodynamic/timeVaryingFreestreamAndPitch"))
writedlm("test/newTestDataGenerators/aerodynamic/timeVaryingFreestreamAndPitch/cn.txt", test_data_matrix(cn[end]))
writedlm("test/newTestDataGenerators/aerodynamic/timeVaryingFreestreamAndPitch/cm.txt", test_data_matrix(cm[end]))