if !isdefined(@__MODULE__,:ChangAeroBeams)
    include(joinpath(@__DIR__,"src","ChangAeroBeams.jl"))
end
using .ChangAeroBeams

if abspath(PROGRAM_FILE) == @__FILE__
    configuration = include(joinpath(@__DIR__,"chang_aerobeams_case.jl"))
    chang_aerobeams = run_chang_aerobeams(configuration)
end
