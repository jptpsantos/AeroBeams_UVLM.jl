# Verify the Chang structural discretization against published modal data.

import Pkg

const EXAMPLE_DIR = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(normpath(joinpath(EXAMPLE_DIR, "..", "..")))

using LinearAlgebra
using StaticArrays
using WingPropellerUVLM:
    Freestream,
    Reference,
    Uniform,
    get_nodal_properties_chang,
    linear_interpolate_1d

const AIR_DENSITY = 1.225

# The published-frequency tolerance below was established for 30 beam elements.
# Fix that discretization here so this check is independent of the interactive
# case default or inherited shell settings.
ENV["CHANG_WING_SPAN_PANELS"] = "30"

include(joinpath(EXAMPLE_DIR, "chang_case.jl"))
include(joinpath(EXAMPLE_DIR, "src", "chang_model_parameters.jl"))
include(joinpath(EXAMPLE_DIR, "src", "chang_structural_model.jl"))

structural = assemble_chang_structural_model()
legacy_structural = assemble_chang_structural_model(inertia_reference = :beam_axis)
diagnostics = chang_structural_diagnostics(structural)
legacy_diagnostics = chang_structural_diagnostics(legacy_structural)

corrected_modes = chang_wing_modal_analysis(structural)
legacy_modes = chang_wing_modal_analysis(legacy_structural)
coupled_solution = eigen(Symmetric(structural.K), Symmetric(structural.M))
coupled_positive = findall(value -> isfinite(value) && value > 0.0, coupled_solution.values)
coupled_positive = coupled_positive[sortperm(coupled_solution.values[coupled_positive])]
coupled_frequency_count = min(10, length(coupled_positive))
coupled_frequencies_hz = [
    sqrt(coupled_solution.values[index]) / (2π)
    for index in coupled_positive[1:coupled_frequency_count]
]

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
println("  wing mass model: $(structural.mass_model)")
println("  wing stiffness model: $(structural.stiffness_model)")
println("  computed wing modes (Hz): $(round.(diagnostics.wing_modal_frequencies_hz, digits=4))")
println("  coupled modes (Hz): $(round.(coupled_frequencies_hz, digits=4))")
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
println("  workbook wing mass (kg): $(sum(M_orig_nodes))")
println("  remapped wing mass (kg): $(sum(m_node_vec))")

conserved_quantity_pairs = (
    (m_node_vec, M_orig_nodes),
    (mass_moment_x_node, M_orig_nodes .* cg_x_orig_nodes),
    (mass_moment_y_node, M_orig_nodes .* cg_y_orig_nodes),
    (mass_moment_z_node, M_orig_nodes .* cg_z_orig_nodes),
    (Ixx_beam_node, Ixx_beam_orig),
    (Iyy_beam_node, Iyy_beam_orig),
    (Izz_beam_node, Izz_beam_orig),
    (Ixy_beam_node, Ixy_beam_orig),
    (Ixz_beam_node, Ixz_beam_orig),
    (Iyz_beam_node, Iyz_beam_orig),
)
maximum_conservation_error = maximum(
    abs(sum(remapped) - sum(workbook)) / max(abs(sum(workbook)), 1.0)
    for (remapped, workbook) in conserved_quantity_pairs
)
println("  maximum conserved-total error: $maximum_conservation_error")

# A maps the complete wing-coordinate vector to the exact attachment motion.
# Its transpose is therefore the only load map that preserves virtual work.
attachment_operator = structural.attachment_operators[1]
attachment_position_from_weights = sum(
    weight * span_nodes[node]
    for (node, weight) in zip(
        prop_attachment_node_pairs[1],
        prop_attachment_weights[1],
    )
)
trial_wing_state = collect(range(-0.2, 0.3; length = NDOF))
trial_attachment_wrench = collect(range(1.0, 2.0; length = ndof))
attachment_virtual_work_error = abs(
    dot(attachment_operator * trial_wing_state, trial_attachment_wrench) -
    dot(trial_wing_state, attachment_operator' * trial_attachment_wrench),
)
println("  attachment position from weights (m): $attachment_position_from_weights")
println("  attachment virtual-work error: $attachment_virtual_work_error")

@assert structural.inertia_reference == :center_of_mass
@assert structural.mass_model == :control_volume_spatial_blocks
@assert structural.stiffness_model == :distributed_integrated
@assert diagnostics.minimum_mass_eigenvalue > 0.0
@assert diagnostics.minimum_stiffness_eigenvalue > 0.0
@assert diagnostics.mass_symmetry_error <= 1.0e-12
@assert diagnostics.stiffness_symmetry_error <= 1.0e-12
@assert all(
    isapprox(sum(remapped), sum(workbook); rtol = 1.0e-12, atol = 1.0e-12)
    for (remapped, workbook) in conserved_quantity_pairs
)
@assert isapprox(
    attachment_position_from_weights,
    propeller_span_positions[1];
    atol = 1.0e-12,
    rtol = 0.0,
)
@assert attachment_virtual_work_error <= 1.0e-12
# The active 30-element model must retain the workbook-axis agreement found in
# the reference modal comparison.  A 2% limit catches the former Y/Z swap,
# whose OOP3 error exceeded 7%.
@assert maximum(abs(row.corrected_error_percent) for row in comparison) <= 2.0
