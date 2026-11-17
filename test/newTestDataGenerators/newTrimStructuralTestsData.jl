# Trim analysis (reaction loads check) of a beam loaded at the middle
include("../examples/trimStructural/freeBeamTrim.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/trimStructural/freeBeamTrim"))
writedlm("test/newTestDataGenerators/trimStructural/freeBeamTrim/F3.txt", test_data_matrix(F3))
writedlm("test/newTestDataGenerators/trimStructural/freeBeamTrim/M2.txt", test_data_matrix(M2))

# Trim analysis (reaction loads check) of a simply-supported beam loaded at the middle
include("../examples/trimStructural/midDispBeamTrim.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/trimStructural/midDispBeamTrim"))
writedlm("test/newTestDataGenerators/trimStructural/midDispBeamTrim/F3.txt", test_data_matrix(F3))
writedlm("test/newTestDataGenerators/trimStructural/midDispBeamTrim/M2.txt", test_data_matrix(M2))

# Trim analysis (reaction loads check) of a right-angled frame
include("../examples/trimStructural/rightAngledFrameTrim.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/trimStructural/rightAngledFrameTrim"))
writedlm("test/newTestDataGenerators/trimStructural/rightAngledFrameTrim/balanceHorizontalForce.txt", test_data_matrix(balanceHorizontalForce))
writedlm("test/newTestDataGenerators/trimStructural/rightAngledFrameTrim/balanceVerticalForce.txt", test_data_matrix(balanceVerticalForce))

# Trim analysis of a cantilever with tip force
include("../examples/trimStructural/tipLoadedCantileverTrim.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/trimStructural/tipLoadedCantileverTrim"))
writedlm("test/newTestDataGenerators/trimStructural/tipLoadedCantileverTrim/u3.txt", test_data_matrix(u3))
writedlm("test/newTestDataGenerators/trimStructural/tipLoadedCantileverTrim/F3.txt", test_data_matrix(F3))
writedlm("test/newTestDataGenerators/trimStructural/tipLoadedCantileverTrim/M2.txt", test_data_matrix(M2))