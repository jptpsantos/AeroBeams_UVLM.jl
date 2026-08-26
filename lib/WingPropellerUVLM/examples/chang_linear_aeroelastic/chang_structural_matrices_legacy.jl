# Legacy compatibility entry point for the pre-parallel-axis Chang matrices.
#
# This file intentionally preserves the former interpretation in which the
# tabulated nodal inertia tensor was inserted directly at the beam reference
# axis. New analyses should call `assemble_chang_structural_matrices`, whose
# default treats the tabulated tensor as a center-of-mass inertia.

isdefined(@__MODULE__, :assemble_chang_structural_matrices) ||
    include(joinpath(@__DIR__, "chang_structural_matrices.jl"))

function assemble_chang_structural_matrices_legacy(; kwargs...)
    legacy_kwargs = merge((; kwargs...), (; inertia_reference = :beam_axis))
    return assemble_chang_structural_matrices(; legacy_kwargs...)
end
