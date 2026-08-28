module WingPropellerUVLM

using CCBlade
using DelimitedFiles
using FLOWMath
using Interpolations
using LinearAlgebra
using StaticArrays
using VSPGeom
using WriteVTK

# VortexLattice main functions
include("backend/nonlinear.jl")
include("backend/rotors.jl")
include("backend/panel.jl")
include("backend/wake.jl")
include("backend/geometry.jl")
include("backend/vspgeom.jl")
include("backend/reference.jl")
include("backend/freestream.jl")
include("backend/induced.jl")
include("backend/circulation.jl")
include("backend/system.jl")
include("backend/analyses.jl")
include("backend/legacy_nearfield.jl")
include("backend/nearfield.jl")
include("backend/nearfield_postprocessing.jl")
include("backend/farfield.jl")
include("backend/stability.jl")
include("backend/visualization.jl")

# Geometry and Kinematics helpers.
include("wing_propeller/BladeGeometry.jl")
include("wing_propeller/GridUtilities.jl")
include("wing_propeller/Kinematics.jl")
include("wing_propeller/Initialization.jl")
include("UVLMState.jl")
include("aeroelastic/GeneralizedAlpha.jl")
include("aeroelastic/Excitations.jl")

export SectionProperties, grid_to_sections, nonlinear_analysis!
export generate_rotor
export SurfacePanel, WakePanel, TrefftzPanel, Wake
export reflect, set_normal
export AbstractSpacing, Uniform, Sine, Cosine
export grid_to_surface_panels, wing_to_grid
export lifting_line_geometry, lifting_line_geometry!
export translate, translate!, rotate, rotate!
export repeated_trailing_edge_points
export read_degengeom, import_vsp
export Reference
export AbstractFrame, Body, Stability, Wind
export Freestream, trajectory_to_freestream
export System, PanelProperties, get_surface_properties
export steady_analysis, steady_analysis!
export unsteady_analysis, unsteady_analysis!, propagate_system!, advance_wake!
export spanwise_force_coefficients
export body_forces, body_forces_history
export lifting_line_coefficients, lifting_line_coefficients!
export near_field_forces!, legacy_near_field_forces!, legacy_imperial_segment_forces!
export imperial_nodal_forces, imperial_nodal_positions
export far_field_drag
export body_derivatives, stability_derivatives
export write_vtk

export span_position_to_node_index
export span_positions_to_node_indices
export propeller_attachment_nodes_from_eta
export linear_interpolate_1d
export generate_panel_grid_and_interpolate
export generate_aero_panel_grid_and_interpolate
export generate_propeller_blades_grid
export get_chord_over_R, get_twist_deg, get_twist_deg_chang
export get_twist_deg_interp
export get_nodal_properties, get_nodal_properties_chang
export RotationMatrix
export copy_surfaces_to_previous!
export wing_kinematics_from_free_state
export initialize_propeller_grids!, update_propeller_grids!
export update_system_surfaces!
export initialize_bohnisch_uvlm_system

export UVLMSnapshot
export snapshot_uvlm, restore_uvlm!
export advance_uvlm_trial!, commit_wake_rows!

export generalized_alpha_parameters, generalized_alpha_kinematics
export generalized_alpha_corrector, generalized_alpha_equilibrium_residual
export PartitionedCouplingOptions, partitioned_generalized_alpha_step
export smooth_hann_pulse, smooth_hann_pulse_load

end # module
