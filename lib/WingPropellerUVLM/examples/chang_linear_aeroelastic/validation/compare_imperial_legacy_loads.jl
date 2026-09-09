# Compare the production Imperial load implementation with the retained legacy
# implementation at exactly the same final aeroelastic state.

const EXAMPLE_DIR = normpath(joinpath(@__DIR__, ".."))
get!(ENV, "CHANG_PLOT_RESULTS", "false")
get!(ENV, "CHANG_ANIMATE_WAKE", "false")
get!(ENV, "CHANG_END_TIME_S", "0.01")
get!(ENV, "CHANG_WING_SPAN_PANELS", "4")
get!(ENV, "CHANG_WING_CHORD_PANELS", "2")
get!(ENV, "CHANG_PROP_RADIAL_PANELS", "2")
get!(ENV, "CHANG_PROP_CHORD_PANELS", "2")
get!(ENV, "CHANG_NEAR_FIELD_FORCE_MODEL", "imperial")
get!(ENV, "CHANG_OUTPUT_DIR", joinpath(EXAMPLE_DIR, "output", "validation", "force_comparison"))
include(joinpath(EXAMPLE_DIR, "run_chang_linear_aeroelastic.jl"))

using LinearAlgebra
using StaticArrays
using WingPropellerUVLM: imperial_nodal_forces, imperial_nodal_positions

# Inspect the objects returned by the run; the adapter receives them explicitly.
(; model, workspace, solution) = chang_run
structural = model.structural
(; ndof_wing_free) = structural
(; t) = model.parameters
(; system, surface_interaction_id) = workspace
import WingPropellerUVLM

final_state = solution.displacement_history[solution.last_step + 1]
kinematics = update_aero_geometry_for_state!(
    model, workspace,
    final_state,
    t[solution.last_step + 1],
)
imperial_load = assemble_structural_aero_load!(model, workspace, kinematics)
imperial_chord = deepcopy(system.chord_seg_forces)
imperial_span = deepcopy(system.span_seg_forces)
imperial_unsteady = deepcopy(system.unsteady_forces)

legacy_properties = deepcopy(system.properties)
_, legacy_chord, legacy_span, legacy_unsteady =
    WingPropellerUVLM.legacy_near_field_forces!(
    legacy_properties,
    system.surfaces,
    system.wakes,
    system.reference[],
    system.freestream[],
    system.Γ;
    dΓdt = system.dΓdt,
    additional_velocity = nothing,
    Vh = system.Vh,
    Vv = system.Vv,
    symmetric = system.symmetric,
    nwake = system.nwake,
    surface_id = system.surface_id,
    wake_finite_core = system.wake_finite_core,
    wake_shedding_locations = system.wake_shedding_locations,
    trailing_vortices = system.trailing_vortices,
    xhat = system.xhat[],
    interaction_id = surface_interaction_id,
    interaction = model.config.simulation.interaction_on,
)

for surface_index in eachindex(system.surfaces)
    system.chord_seg_forces[surface_index] .= legacy_chord[surface_index]
    system.span_seg_forces[surface_index] .= legacy_span[surface_index]
    system.unsteady_forces[surface_index] .= legacy_unsteady[surface_index]
end
legacy_load = assemble_structural_aero_load!(model, workspace, kinematics)

function component_load(chord, span, unsteady, selected::Symbol)
    zero_force = SVector{3,Float64}(0.0, 0.0, 0.0)
    for surface_index in eachindex(system.surfaces)
        fill!(system.chord_seg_forces[surface_index], zero_force)
        fill!(system.span_seg_forces[surface_index], zero_force)
        fill!(system.unsteady_forces[surface_index], zero_force)
        selected == :chord &&
            (system.chord_seg_forces[surface_index] .= chord[surface_index])
        selected == :span &&
            (system.span_seg_forces[surface_index] .= span[surface_index])
        selected == :unsteady &&
            (system.unsteady_forces[surface_index] .= unsteady[surface_index])
    end
    return assemble_structural_aero_load!(model, workspace, kinematics)
end

imperial_components = Dict(
    component => component_load(
        imperial_chord, imperial_span, imperial_unsteady, component,
    ) for component in (:chord, :span, :unsteady)
)
legacy_components = Dict(
    component => component_load(
        legacy_chord, legacy_span, legacy_unsteady, component,
    ) for component in (:chord, :span, :unsteady)
)

function scaled_difference(a, b)
    return norm(a - b) / max(norm(a), norm(b), eps(Float64))
end

wing_indices = 1:ndof_wing_free
propeller_indices = (ndof_wing_free + 1):length(imperial_load)

println("IMPERIAL_LEGACY_LOAD_COMPARISON")
println("time_s=", t[solution.last_step + 1])
println("all_finite=", all(isfinite, imperial_load) && all(isfinite, legacy_load))
println("relative_full_load_difference=", scaled_difference(imperial_load, legacy_load))
println("relative_wing_load_difference=", scaled_difference(
    imperial_load[wing_indices], legacy_load[wing_indices]))
println("relative_propeller_load_difference=", scaled_difference(
    imperial_load[propeller_indices], legacy_load[propeller_indices]))
println("imperial_propeller_load=", imperial_load[propeller_indices])
println("legacy_propeller_load=", legacy_load[propeller_indices])
for component in (:chord, :span, :unsteady)
    imperial_component = imperial_components[component]
    legacy_component = legacy_components[component]
    println(component, "_relative_full_difference=",
        scaled_difference(imperial_component, legacy_component))
    println(component, "_imperial_propeller_load=",
        imperial_component[propeller_indices])
    println(component, "_legacy_propeller_load=",
        legacy_component[propeller_indices])
end

println("imperial_trailing_span_row_norms=", [
    norm(imperial_span[index][end, :]) for index in eachindex(imperial_span)
])
println("legacy_trailing_span_row_norms=", [
    norm(legacy_span[index][end, :]) for index in eachindex(legacy_span)
])
