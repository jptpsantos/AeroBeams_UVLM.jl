# Run a short diagnostic case, then verify aerodynamic load transfer by virtual work.

include(joinpath(@__DIR__, "run_chang_diagnostic_case.jl"))

base_state = copy(U[N_LAST + 1])
base_time = t[N_LAST + 1]
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

attachment_node = prop_attach_nodes[1]
attachment_start = ndof * (attachment_node - 2)
indices = [
    ("attachment_chord_translation", attachment_start + 2),
    ("attachment_vertical_translation", attachment_start + 3),
    ("attachment_span_rotation", attachment_start + 4),
    ("attachment_chord_rotation", attachment_start + 5),
    ("attachment_vertical_rotation", attachment_start + 6),
    ("propeller_pitch", ndof_wing_free + 1),
    ("propeller_yaw", ndof_wing_free + 2),
]

println("\nChang aerodynamic load-transfer virtual-work audit")
println("quantity,generalized_load,virtual_work_load,scaled_error")
maximum_error = 0.0
for (name, index) in indices
    work_load = virtual_work_load(index)
    assembled_load = generalized_load[index]
    scale = max(abs(assembled_load), abs(work_load), 1.0)
    error = abs(assembled_load - work_load) / scale
    global maximum_error = max(maximum_error, error)
    println("$name,$assembled_load,$work_load,$error")
end
update_aero_geometry_for_state!(system, base_state, base_time)
println("maximum_scaled_virtual_work_error = $maximum_error")
