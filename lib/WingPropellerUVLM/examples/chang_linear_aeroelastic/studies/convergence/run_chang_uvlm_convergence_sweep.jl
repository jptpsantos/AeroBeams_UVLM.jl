# Legacy one-stage Chang aeroelastic convergence entry point.

include(joinpath(@__DIR__, "chang_aeroelastic_convergence.jl"))

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
