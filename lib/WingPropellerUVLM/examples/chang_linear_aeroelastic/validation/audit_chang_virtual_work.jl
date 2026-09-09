# Run a short case, then verify aerodynamic load transfer by virtual work.

const EXAMPLE_DIR = normpath(joinpath(@__DIR__, ".."))
get!(ENV, "CHANG_PLOT_RESULTS", "false")
get!(ENV, "CHANG_ANIMATE_WAKE", "false")
get!(ENV, "CHANG_END_TIME_S", "0.01")
get!(ENV, "CHANG_WING_SPAN_PANELS", "4")
get!(ENV, "CHANG_WING_CHORD_PANELS", "2")
get!(ENV, "CHANG_PROP_RADIAL_PANELS", "2")
get!(ENV, "CHANG_PROP_CHORD_PANELS", "2")
get!(ENV, "CHANG_OUTPUT_DIR", joinpath(EXAMPLE_DIR, "output", "validation", "virtual_work"))
include(joinpath(EXAMPLE_DIR, "run_chang_linear_aeroelastic.jl"))

using LinearAlgebra
using StaticArrays
using WingPropellerUVLM: imperial_nodal_forces, imperial_nodal_positions

# Inspect the objects returned by the run; the adapter receives them explicitly.
(; model, workspace, solution) = chang_run
structural = model.structural
(; ndof_wing_free) = structural
(;
   span_length, nnodes, ndof, chord, xle_distribution,
   span_nodes, Npropellers, prop_attachment_node_pairs, prop_attachment_weights, t
) = model.parameters
(; system, surface_interaction_id, T_pivot_A_current, prop_surface_indices) = workspace
elastic_axis_fraction = model.aerodynamic_options.elastic_axis_fraction

base_state = copy(solution.displacement_history[solution.last_step + 1])
base_time = t[solution.last_step + 1]
if haskey(ENV, "CHANG_AUDIT_PROP_PITCH_DEG")
    base_state[ndof_wing_free + 1] = deg2rad(parse(
        Float64,
        ENV["CHANG_AUDIT_PROP_PITCH_DEG"],
    ))
end
if haskey(ENV, "CHANG_AUDIT_PROP_YAW_DEG")
    base_state[ndof_wing_free + 2] = deg2rad(parse(
        Float64,
        ENV["CHANG_AUDIT_PROP_YAW_DEG"],
    ))
end
base_kinematics = update_aero_geometry_for_state!(model, workspace, base_state, base_time)
generalized_load = assemble_structural_aero_load!(model, workspace, base_kinematics)
nodal_forces = imperial_nodal_forces(system)

function vortex_positions_at(state)
    update_aero_geometry_for_state!(model, workspace, state, base_time)
    return imperial_nodal_positions(system)
end

function virtual_work_load(index; state = base_state, step = 1.0e-7)
    plus_state = copy(state)
    minus_state = copy(state)
    plus_state[index] += step
    minus_state[index] -= step
    plus_positions = vortex_positions_at(plus_state)
    minus_positions = vortex_positions_at(minus_state)

    work = 0.0
    for surface_index in eachindex(nodal_forces)
        for vertex_index in eachindex(nodal_forces[surface_index])
            derivative = (
                plus_positions[surface_index][vertex_index] -
                minus_positions[surface_index][vertex_index]
            ) / (2 * step)
            work += dot(nodal_forces[surface_index][vertex_index], derivative)
        end
    end
    return work
end

left_attachment_node, right_attachment_node = prop_attachment_node_pairs[1]
indices = Tuple{String,Int}[]
for node in 2:nnodes
    prefix = node == left_attachment_node ? "left_attachment" :
        node == right_attachment_node ? "right_attachment" : "wing_node_$node"
    attachment_start = ndof * (node - 2)
    append!(indices, [
        ("$(prefix)_span_translation", attachment_start + 1),
        ("$(prefix)_chord_translation", attachment_start + 2),
        ("$(prefix)_vertical_translation", attachment_start + 3),
        ("$(prefix)_span_rotation", attachment_start + 4),
        ("$(prefix)_chord_rotation", attachment_start + 5),
        ("$(prefix)_vertical_rotation", attachment_start + 6),
    ])
end
for propeller_index in 1:Npropellers
    prefix = propeller_index == 1 ? "propeller" : "propeller_$propeller_index"
    offset = ndof_wing_free + 2 * (propeller_index - 1)
    append!(indices, [("$(prefix)_pitch", offset + 1), ("$(prefix)_yaw", offset + 2)])
end

function audit_state_virtual_work(state; label)
    kinematics = update_aero_geometry_for_state!(model, workspace, state, base_time)
    mapped_load = assemble_structural_aero_load!(model, workspace, kinematics)
    errors = Dict{String,Float64}()
    maximum_checked_error = 0.0
    println("\nstate=$label")
    println("quantity,generalized_load,virtual_work_load,scaled_error")
    try
        for (name, index) in indices
            work_load = virtual_work_load(index; state)
            assembled_load = mapped_load[index]
            scale = max(abs(assembled_load), abs(work_load), 1.0)
            error = abs(assembled_load - work_load) / scale
            errors[name] = error
            println("$name,$assembled_load,$work_load,$error")
            # Wing projection is exact for either propeller modal option.
            # The retained fixed-axis modal option is intentionally approximate.
            if index <= ndof_wing_free ||
                model.config.simulation.propeller_moment_projection == :exact_virtual_work
                maximum_checked_error = max(maximum_checked_error, error)
                @assert error <= 1.0e-6 "$label: virtual-work mismatch for $name ($error)"
            end
        end
    finally
        update_aero_geometry_for_state!(model, workspace, base_state, base_time)
    end
    println("maximum_checked_virtual_work_error = $maximum_checked_error")
    return (; errors, maximum_checked_error)
end

println("\nChang aerodynamic load-transfer virtual-work audit")
println("propeller_moment_projection = $(model.config.simulation.propeller_moment_projection)")
println("attachment_nodes = $(prop_attachment_node_pairs[1])")
println("attachment_weights = $(prop_attachment_weights[1])")
println(
    "audit_propeller_angles_deg = " *
    "($(rad2deg(base_state[ndof_wing_free + 1])), " *
    "$(rad2deg(base_state[ndof_wing_free + 2])))",
)
reference_audit = audit_state_virtual_work(base_state; label = "reference")
errors = reference_audit.errors
maximum_error = reference_audit.maximum_checked_error

# Frozen forces isolate load transfer from the circulation solver. Exercise
# both a uniform rotation and the interpolation chain rule at attachments.
for pattern in (:uniform, :nonuniform)
    trial_state = copy(base_state)
    for node in 2:nnodes
        offset = ndof * (node - 2)
        eta = span_nodes[node] / span_length
        trial_state[offset+1:offset+3] .= (0.001eta, -0.002eta^2, 0.003eta)
        trial_state[offset+4:offset+6] .= pattern == :uniform ? fill(deg2rad(1.0), 3) :
            deg2rad(5.0) .* [sin(pi*eta/2), -0.7eta^2, 0.6sin(pi*eta)]
    end
    for propeller_index in 1:Npropellers
        offset = ndof_wing_free + 2 * (propeller_index - 1)
        trial_state[offset+1:offset+2] .= deg2rad.([4.0, 3.0])
    end
    audit = audit_state_virtual_work(trial_state; label = string(pattern))
    global maximum_error = max(maximum_error, audit.maximum_checked_error)
end
println("maximum_scaled_virtual_work_error = $maximum_error")
