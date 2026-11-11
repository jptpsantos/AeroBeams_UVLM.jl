if !isdefined(@__MODULE__,:ChangAeroBeams)
    include(joinpath(@__DIR__,"src","ChangAeroBeams.jl"))
end
using .ChangAeroBeams

if abspath(PROGRAM_FILE) == @__FILE__
    base = include(joinpath(@__DIR__,"chang_aerobeams_case.jl"))
    configuration = merge(
        base,
        (
            solver=merge(base.solver,(support=:beam_pylon,)),
            output=merge(base.output,(label="chang_aerobeams_beam_pylon",)))
    )
    chang_beam_pylon = run_chang_aerobeams(configuration)
end
