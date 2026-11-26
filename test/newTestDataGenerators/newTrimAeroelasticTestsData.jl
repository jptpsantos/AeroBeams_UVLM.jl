# Trim analysis of the Blended-Wing-Body flying wing in free flight
include("../examples/trimAeroelastic/BWBtrim.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/trimAeroelastic/BWBtrim"))
writedlm("test/newTestDataGenerators/trimAeroelastic/BWBtrim/trimAoA.txt", test_data_matrix(trimAoA))
writedlm("test/newTestDataGenerators/trimAeroelastic/BWBtrim/trimThrust.txt", test_data_matrix(trimThrust))
writedlm("test/newTestDataGenerators/trimAeroelastic/BWBtrim/trimDelta.txt", test_data_matrix(trimδ))

# Trim analysis of the conventional HALE aircraft in free flight (considering aerodynamics from stabilizers and thrust)
include("../examples/trimAeroelastic/cHALEfullTrim.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/trimAeroelastic/cHALEfullTrim"))
writedlm("test/newTestDataGenerators/trimAeroelastic/cHALEfullTrim/trimAoA.txt", test_data_matrix(trimAoA))
writedlm("test/newTestDataGenerators/trimAeroelastic/cHALEfullTrim/trimThrust.txt", test_data_matrix(trimThrust))
writedlm("test/newTestDataGenerators/trimAeroelastic/cHALEfullTrim/trimDelta.txt", test_data_matrix(trimδ))

# Trim analysis of the conventional HALE aircraft in free flight at rigid and flexible configurations (neglecting aerodynamics from stabilizers and thrust)
include("../examples/trimAeroelastic/cHALEtrim.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/trimAeroelastic/cHALEtrim"))
writedlm("test/newTestDataGenerators/trimAeroelastic/cHALEtrim/trimAoA.txt", test_data_matrix(trimAoA))

# Trim analysis of the Helios flying-wing
include("../examples/trimAeroelastic/heliosTrim.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/trimAeroelastic/heliosTrim"))
writedlm("test/newTestDataGenerators/trimAeroelastic/heliosTrim/trimAoA.txt", test_data_matrix(trimAoA))
writedlm("test/newTestDataGenerators/trimAeroelastic/heliosTrim/trimThrust.txt", test_data_matrix(trimThrust))
writedlm("test/newTestDataGenerators/trimAeroelastic/heliosTrim/trimDelta.txt", test_data_matrix(trimδ))