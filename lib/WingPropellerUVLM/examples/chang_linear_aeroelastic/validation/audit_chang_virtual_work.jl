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
base_kinematics = update_aero_geometry_for_state!(system, base_state, base_time)
generalized_load = assemble_structural_aero_load!(system, base_kinematics)
nodal_forces = imperial_nodal_forces(system)

function vortex_positions_at(state)
    update_aero_geometry_for_state!(system, state, base_time)
    return imperial_nodal_positions(system)
end

function virtual_work_load(index; step = 1.0e-7)
    plus_state = copy(base_state)
    minus_state = copy(base_state)
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
for (side, node) in (
    ("left", left_attachment_node),
    ("right", right_attachment_node),
)
    attachment_start = ndof * (node - 2)
    append!(indices, [
        ("$(side)_attachment_chord_translation", attachment_start + 2),
        ("$(side)_attachment_vertical_translation", attachment_start + 3),
        ("$(side)_attachment_span_rotation", attachment_start + 4),
        ("$(side)_attachment_chord_rotation", attachment_start + 5),
        ("$(side)_attachment_vertical_rotation", attachment_start + 6),
    ])
end
append!(indices, [
    ("propeller_pitch", ndof_wing_free + 1),
    ("propeller_yaw", ndof_wing_free + 2),
])

println("\nChang aerodynamic load-transfer virtual-work audit")
println("propeller_moment_projection = $AEROELASTIC_PROPELLER_MOMENT_PROJECTION")
println("attachment_nodes = $(prop_attachment_node_pairs[1])")
println("attachment_weights = $(prop_attachment_weights[1])")
println(
    "audit_propeller_angles_deg = " *
    "($(rad2deg(base_state[ndof_wing_free + 1])), " *
    "$(rad2deg(base_state[ndof_wing_free + 2])))",
)
println("quantity,generalized_load,virtual_work_load,scaled_error")
maximum_error = 0.0
errors = Dict{String,Float64}()
for (name, index) in indices
    work_load = virtual_work_load(index)
    assembled_load = generalized_load[index]
    scale = max(abs(assembled_load), abs(work_load), 1.0)
    error = abs(assembled_load - work_load) / scale
    errors[name] = error
    global maximum_error = max(maximum_error, error)
    println("$name,$assembled_load,$work_load,$error")
end
update_aero_geometry_for_state!(system, base_state, base_time)
println("maximum_scaled_virtual_work_error = $maximum_error")
if AEROELASTIC_PROPELLER_MOMENT_PROJECTION == :exact_virtual_work
    @assert errors["propeller_pitch"] <= 1.0e-6
    @assert errors["propeller_yaw"] <= 1.0e-6
end
