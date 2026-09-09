# Regression audit of the corrected production map against geometry derivatives.
# Run from the repository root:
# julia --startup-file=no --compiled-modules=existing --project=lib/WingPropellerUVLM <this file>
#
# The existing audit provides a small UVLM case and independent geometry FD.
# Explicit settings reproduce the original review, regardless of editable defaults.
const INVESTIGATION_DIR = joinpath(@__DIR__, "..", "output", "wing_moment_investigation")
mkpath(INVESTIGATION_DIR)
pushfirst!(DEPOT_PATH, abspath(joinpath(INVESTIGATION_DIR, "depot")))
for (key, value) in (
    "CHANG_SPEED_MPS" => "84.0", "CHANG_INTERACTION" => "true",
    "CHANG_FCORE_SEGMENT_FACTOR" => "0.025", "CHANG_FCORE_CHORD_FACTOR" => "0.025",
    "CHANG_PROP_MOMENT_PROJECTION" => "exact_virtual_work",
    "CHANG_AUDIT_PROP_PITCH_DEG" => "4.0", "CHANG_AUDIT_PROP_YAW_DEG" => "3.0",
    "CHANG_PLOT_RESULTS" => "false", "CHANG_ANIMATE_WAKE" => "false",
    "CHANG_END_TIME_S" => "0.01", "CHANG_WING_SPAN_PANELS" => "4",
    "CHANG_WING_CHORD_PANELS" => "2", "CHANG_PROP_RADIAL_PANELS" => "2",
    "CHANG_PROP_CHORD_PANELS" => "2", "CHANG_OUTPUT_DIR" => abspath(INVESTIGATION_DIR),
)
    ENV[key] = value
end
include(joinpath(@__DIR__, "audit_chang_virtual_work.jl"))
using Test

# For R = Rz(theta_z) Rx(theta_x) Ry(theta_y):
# a_span = Rz Rx e_y, a_chord = Rz e_x, a_down = -e_z.
# theta_z is the negative of the structural down-rotation coordinate.
function analytic_wing_moments(moment, theta_x, theta_z)
    sx, cx = sincos(theta_x)
    sz, cz = sincos(theta_z)
    return SVector(
        -sz * cx * moment[1] + cz * cx * moment[2] + sx * moment[3],
        cz * moment[1] + sz * moment[2],
        -moment[3],
    )
end

# The correction uses spatial moments, not finite-difference derivatives.
# Keep separate wing/propeller parts to identify the source of the discrepancy.
function projection_correction(kinematics)
    positions = imperial_nodal_positions(system)
    forces = imperial_nodal_forces(system)
    delta_wing, delta_prop = zeros(structural.ndof_free), zeros(structural.ndof_free)
    old_wing, old_prop = zeros(structural.ndof_free), zeros(structural.ndof_free)

    function accumulate!(delta, old, node, moment, theta_x, theta_z, weight)
        node == 1 && return
        offset = ndof * (node - 2)
        fixed = SVector(moment[2], moment[1], -moment[3])
        projected = analytic_wing_moments(moment, theta_x, theta_z)
        old[offset+4:offset+6] .+= weight .* fixed
        delta[offset+4:offset+6] .+= weight .* (projected - fixed)
    end

    for node in 1:nnodes
        ea = SVector(
            xle_distribution[node] + chord[node] * elastic_axis_fraction + kinematics.u_x_A[node],
            span_nodes[node] + kinematics.u_y_A[node], kinematics.u_z_A[node])
        moment = sum(cross(positions[1][j,node] - ea, forces[1][j,node])
            for j in axes(forces[1], 1))
        accumulate!(delta_wing, old_wing, node, moment,
            kinematics.theta_x_A[node], kinematics.theta_z_A[node], 1.0)
    end
    for p in 1:Npropellers
        pivot = T_pivot_A_current[p]
        moment = sum(cross(positions[s][j] - pivot, forces[s][j])
            for s in prop_surface_indices[p] for j in eachindex(forces[s]))
        left, right = prop_attachment_node_pairs[p]
        wl, wr = prop_attachment_weights[p]
        tx = wl * kinematics.theta_x_A[left] + wr * kinematics.theta_x_A[right]
        tz = wl * kinematics.theta_z_A[left] + wr * kinematics.theta_z_A[right]
        for (node, weight) in ((left, wl), (right, wr))
            accumulate!(delta_prop, old_prop, node, moment, tx, tz, weight)
        end
    end
    return (; delta_wing, delta_prop, old_wing, old_prop)
end

const INVESTIGATION_REFERENCE_STATE = copy(base_state)
function investigation_state(angle; pattern=:uniform, propeller_motion=true)
    state = copy(INVESTIGATION_REFERENCE_STATE)
    state[1:ndof_wing_free] .= 0.0
    propeller_motion || (state[ndof_wing_free+1:end] .= 0.0)
    for node in 2:nnodes
        offset = ndof * (node - 2)
        if pattern == :nonuniform
            # Different rotations at neighboring nodes exercise the attachment
            # interpolation chain rule. Root translations/rotations stay zero.
            eta = span_nodes[node] / span_length
            state[offset+1:offset+3] .= (0.001eta, -0.002eta^2, 0.003sin(pi*eta/2))
            state[offset+4:offset+6] .= deg2rad(angle) .* (
                sin(pi*eta/2), -0.7eta^2, 0.6sin(pi*eta))
        elseif pattern == :uniform
            state[offset+4:offset+6] .= deg2rad(angle)
        else
            component = Dict(:span => 4, :chord => 5, :vertical => 6)[pattern]
            state[offset+component] = deg2rad(angle)
        end
    end
    return state
end

function evaluate_projection(state)
    kinematics = update_aero_geometry_for_state!(model, workspace, state, base_time)
    mapped = assemble_structural_aero_load!(model, workspace, kinematics)
    parts = projection_correction(kinematics)
    # Reconstruct the former fixed-axis result for historical comparison.
    # The production result itself is tested against independent geometry FD.
    fixed = copy(mapped)
    for node in 2:nnodes, component in 4:6
        index = ndof * (node - 2) + component
        fixed[index] = parts.old_wing[index] + parts.old_prop[index]
    end
    return (; mapped, fixed, parts)
end

function investigation_finite_difference(state, step)
    # virtual_work_load uses the existing audit's base_state and fixed nodal forces.
    base_state .= state
    return [virtual_work_load(index; step) for index in eachindex(state)]
end

function dof_name(index)
    index > ndof_wing_free && return "propeller_$(index - ndof_wing_free)"
    node = div(index - 1, ndof) + 2
    component = ("span_translation", "chord_translation", "down_translation",
        "span_rotation", "chord_rotation", "down_rotation")[mod1(index, ndof)]
    return "node_$(node)_$component"
end

scenarios = [(label="uniform_$angle", angle, pattern=:uniform, propeller_motion=true)
    for angle in (0.0, 0.001, 0.01, 0.1, 1.0, 5.0)]
append!(scenarios, [(label="pure_$pattern", angle=1.0, pattern, propeller_motion=true)
    for pattern in (:span, :chord, :vertical)])
append!(scenarios, [
    (label="nonuniform", angle=5.0, pattern=:nonuniform, propeller_motion=true),
    (label="zero_propeller_angles", angle=1.0, pattern=:uniform, propeller_motion=false)])
fd_steps = (1e-5, 1e-6, 1e-7)
summary_rows = NamedTuple[]
@testset "Production wing projection against actual geometry" begin
    open(joinpath(INVESTIGATION_DIR, "components.csv"), "w") do stream
        println(stream, "scenario,fd_step,dof,fixed_axis,production,finite_difference,fixed_axis_scaled_error,production_scaled_error,wing_correction,attachment_correction")
        for scenario in scenarios
            state = investigation_state(scenario.angle; pattern=scenario.pattern,
                propeller_motion=scenario.propeller_motion)
            evaluation = evaluate_projection(state)
            for step in fd_steps
                fd_load = investigation_finite_difference(state, step)
                scales = max.(abs.(evaluation.mapped), abs.(fd_load), 1.0)
                fixed_scales = max.(abs.(evaluation.fixed), abs.(fd_load), 1.0)
                fixed_error = abs.(evaluation.fixed - fd_load) ./ fixed_scales
                production_error = abs.(evaluation.mapped - fd_load) ./ scales
                # Independent geometry derivatives verify ALL free coordinates.
                @test maximum(production_error) < 1e-6
                @test evaluation.mapped ≈ evaluation.fixed +
                    evaluation.parts.delta_wing + evaluation.parts.delta_prop
                if scenario.angle == 0 || scenario.pattern == :span
                    @test maximum(fixed_error) < 1e-6
                end
                worst = argmax(fixed_error)
                push!(summary_rows, (; scenario=scenario.label, step,
                    worst_dof=dof_name(worst), fixed_error=maximum(fixed_error),
                    production_error=maximum(production_error), fixed=evaluation.fixed[worst],
                    production=evaluation.mapped[worst], finite_difference=fd_load[worst],
                    wing_delta=evaluation.parts.delta_wing[worst],
                    attachment_delta=evaluation.parts.delta_prop[worst]))
                for index in eachindex(state)
                    println(stream, join((scenario.label, step, dof_name(index),
                        evaluation.fixed[index], evaluation.mapped[index], fd_load[index],
                        fixed_error[index], production_error[index],
                        evaluation.parts.delta_wing[index], evaluation.parts.delta_prop[index]), ","))
                end
            end
        end
    end
end
open(joinpath(INVESTIGATION_DIR, "summary.csv"), "w") do stream
    println(stream, join(string.(keys(first(summary_rows))), ","))
    for row in summary_rows
        println(stream, join(values(row), ","))
        row.step == 1e-6 && println("PROJECTION_RESULT: ", row)
    end
end

# A nonzero preload makes the restored axis derivative first-order in q.
# This check freezes aerodynamic forces and uses zero wing AND propeller angles.
# It therefore isolates the load-map contribution to the tangent, not the
# complete aeroelastic stiffness or flutter damping.
@testset "Preload contribution near zero state" begin
    angle_step = 1e-5
    state_plus = investigation_state(rad2deg(angle_step); propeller_motion=false)
    state_minus = investigation_state(-rad2deg(angle_step); propeller_motion=false)
    plus, minus = evaluate_projection(state_plus), evaluate_projection(state_minus)
    tangent = ((plus.mapped - plus.fixed) - (minus.mapped - minus.fixed)) / (2angle_step)
    @test maximum(abs, tangent) > 1.0
    finer_plus = evaluate_projection(investigation_state(rad2deg(angle_step/2); propeller_motion=false))
    finer_minus = evaluate_projection(investigation_state(-rad2deg(angle_step/2); propeller_motion=false))
    finer_tangent = ((finer_plus.mapped-finer_plus.fixed) -
        (finer_minus.mapped-finer_minus.fixed)) / angle_step
    @test isapprox(tangent, finer_tangent; rtol=1e-7, atol=1e-6)
    open(joinpath(INVESTIGATION_DIR, "restored_preload_tangent.csv"), "w") do stream
        println(stream, "dof,restored_directional_tangent_Nm_per_rad")
        for i in eachindex(tangent)
            println(stream, "$(dof_name(i)),$(tangent[i])")
        end
    end
    println("RESTORED_PRELOAD_TANGENT_MAX: ", maximum(abs, tangent), " N*m/rad")
end
base_state .= INVESTIGATION_REFERENCE_STATE
update_aero_geometry_for_state!(model, workspace, base_state, base_time)
println("Investigation outputs: ", abspath(INVESTIGATION_DIR))
