import Pkg

Pkg.activate(normpath(joinpath(@__DIR__, "..", "..")))

using LinearAlgebra
using StaticArrays
using WingPropellerUVLM:
    Freestream,
    Reference,
    Uniform,
    get_nodal_properties_chang,
    linear_interpolate_1d

const AIR_DENSITY = 1.225

include(joinpath(@__DIR__, "chang_case.jl"))
include(joinpath(@__DIR__, "chang_model_parameters.jl"))
include(joinpath(@__DIR__, "chang_structural_model.jl"))

structural = assemble_chang_structural_model()
legacy_structural = assemble_chang_structural_model(inertia_reference = :beam_axis)
diagnostics = chang_structural_diagnostics(structural)
legacy_diagnostics = chang_structural_diagnostics(legacy_structural)

corrected_modes = chang_wing_modal_analysis(structural)
legacy_modes = chang_wing_modal_analysis(legacy_structural)

# Isolated-wing values from Chang et al., reproduced in Table 2 of Santos,
# Marques, and Riso, "Wing-Propeller Aeroelastic Stability Under Aerodynamic
# Interactions Using the Unsteady Vortex-Lattice Method," VFS Forum 82, 2026.
reference_modes = [
    (label = "OOP1", family = :out_of_plane, family_order = 1, frequency_hz = 6.63),
    (label = "IP1",  family = :in_plane,     family_order = 1, frequency_hz = 20.93),
    (label = "OOP2", family = :out_of_plane, family_order = 2, frequency_hz = 36.66),
    (label = "T1",   family = :torsion,      family_order = 1, frequency_hz = 43.50),
    (label = "OOP3", family = :out_of_plane, family_order = 3, frequency_hz = 93.31),
]

function matching_frequency(modes, reference)
    index = findfirst(mode ->
        mode.family == reference.family && mode.family_order == reference.family_order,
        modes,
    )
    isnothing(index) && error("Could not identify reference mode $(reference.label)")
    return modes[index].frequency_hz
end

comparison = map(reference_modes) do reference
    corrected_hz = matching_frequency(corrected_modes, reference)
    legacy_hz = matching_frequency(legacy_modes, reference)
    return (;
        reference.label,
        reference_hz = reference.frequency_hz,
        corrected_hz,
        corrected_error_percent = 100 * (corrected_hz / reference.frequency_hz - 1),
        legacy_hz,
        legacy_error_percent = 100 * (legacy_hz / reference.frequency_hz - 1),
    )
end

println("\nChang structural-model verification")
println("  inertia reference: $(structural.inertia_reference)")
println("  computed wing modes (Hz): $(round.(diagnostics.wing_modal_frequencies_hz, digits=4))")
println("  legacy wing modes (Hz): $(round.(legacy_diagnostics.wing_modal_frequencies_hz, digits=4))")
println("  corrected mode classification: ", [
    "$(mode.family)$(mode.family_order)" for mode in corrected_modes
])
println("  legacy mode classification: ", [
    "$(mode.family)$(mode.family_order)" for mode in legacy_modes
])
println("\n  Mode      Reference    Corrected   Error (%)     Legacy   Error (%)")
for row in comparison
    println(
        "  ", rpad(row.label, 8),
        lpad(round(row.reference_hz, digits=4), 10),
        lpad(round(row.corrected_hz, digits=4), 13),
        lpad(round(row.corrected_error_percent, digits=3), 12),
        lpad(round(row.legacy_hz, digits=4), 11),
        lpad(round(row.legacy_error_percent, digits=3), 12),
    )
end
println("  min eig(M): $(diagnostics.minimum_mass_eigenvalue)")
println("  min eig(K): $(diagnostics.minimum_stiffness_eigenvalue)")
println("  M symmetry error: $(diagnostics.mass_symmetry_error)")
println("  K symmetry error: $(diagnostics.stiffness_symmetry_error)")
println("  norm of symmetric part of C: $(diagnostics.damping_symmetric_part_norm)")

@assert structural.inertia_reference == :center_of_mass
@assert diagnostics.minimum_mass_eigenvalue > 0.0
@assert diagnostics.minimum_stiffness_eigenvalue > 0.0
@assert diagnostics.mass_symmetry_error <= 1.0e-12
@assert diagnostics.stiffness_symmetry_error <= 1.0e-12
# The present 20-element, lumped-mass discretization is accepted within 8%.
# OOP3 is the limiting mode; the first four reference modes remain within 4%.
@assert maximum(abs(row.corrected_error_percent) for row in comparison) <= 8.0
