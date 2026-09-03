# This file is injected immediately after the production
# chang_model_parameters.jl by run_conservative_diagnostic_case.jl.
# It replaces only the nodal wing mass/CG/inertia arrays; distributed stiffness,
# pylon/propeller properties, aerodynamics, and time integration remain those of
# the production driver.

include(joinpath(@__DIR__, "ConservativeLumpedModal.jl"))
using .ConservativeLumpedModal

experimental_properties = remapped_nodal_properties(Ne; mesh = :uniform)
isapprox(experimental_properties.nodes, X_nodes_new; rtol = 0.0, atol = 1e-10) ||
    error("Experimental inertia-remap nodes do not match the active beam mesh")

global m_node_vec = experimental_properties.mass
global cg_x_node_vec = experimental_properties.cg_x
global cg_y_node_vec = experimental_properties.cg_y
global cg_z_node_vec = experimental_properties.cg_z
global Ixx_node_vec = experimental_properties.Ixx
global Iyy_node_vec = experimental_properties.Iyy
global Izz_node_vec = experimental_properties.Izz
global Ixy_node_vec = experimental_properties.Ixy
global Ixz_node_vec = experimental_properties.Ixz
global Iyz_node_vec = experimental_properties.Iyz

experimental_total_block = reduce(+, experimental_properties.blocks)
experimental_source_block = reduce(+, source_spatial_inertia_blocks())
experimental_conservation_error = norm(
    experimental_total_block - experimental_source_block,
    Inf,
) / max(norm(experimental_source_block, Inf), 1.0)

@assert all(m_node_vec[2:end] .> 0.0)
@assert experimental_conservation_error < 1e-12
println(
    "[conservative-remap] applied to Ne=$Ne: mass=$(sum(m_node_vec)) kg, " *
    "spatial-inertia error=$experimental_conservation_error, " *
    "all $Ne free nodal lumps positive",
)
