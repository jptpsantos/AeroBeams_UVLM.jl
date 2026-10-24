# Plot the continuous test's saved strain histories without running a simulation.
include(joinpath(@__DIR__, "plot_pazy_LCO.jl"))

plot_pazy_continuous_LCO(; kwargs...) = plot_pazy_LCO(; sweep_mode=:continuous, kwargs...)

if abspath(PROGRAM_FILE) == @__FILE__
    plot_pazy_continuous_LCO()
end
