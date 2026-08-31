function finite_maximum(values)
    finite_values = filter(isfinite, collect(values))
    return isempty(finite_values) ? NaN : maximum(finite_values)
end

"""Write the Chang response history, validation summary, and optional plot."""
function write_chang_results(;
    displacement_history,
    time,
    time_steps,
    last_step::Integer,
    wing_node_count::Integer,
    degrees_of_freedom_per_node::Integer,
    number_of_propellers::Integer,
    number_of_blades::Integer,
    propeller_eta,
    span_length::Real,
    density::Real,
    freestream_speed::Real,
    interaction_on::Bool,
    near_field_force_model::Symbol,
    propeller_moment_projection::Symbol,
    requested_end_time::Real,
    coupling_iterations,
    coupling_state_residual,
    coupling_load_residual,
    coupling_equilibrium_residual,
    coupling_converged,
    output_directory::AbstractString,
    output_label::AbstractString,
    plot_results::Bool,
    plot_time_limit::Real,
)
    println("Extracting validation results...")

    active_step_count = min(last_step, length(time_steps))
    number_of_states = min(active_step_count + 1, length(time))
    active_steps = 1:max(active_step_count, 1)
    time_history = collect(time[1:number_of_states])

    wing_dofs = (wing_node_count - 1) * degrees_of_freedom_per_node
    tip_displacement_index = wing_dofs - degrees_of_freedom_per_node + 3
    tip_twist_index = wing_dofs - degrees_of_freedom_per_node + 4
    tip_displacement = [
        displacement_history[index][tip_displacement_index]
        for index in 1:number_of_states
    ]
    tip_twist = [
        displacement_history[index][tip_twist_index]
        for index in 1:number_of_states
    ]
    propeller_pitch = [
        [
            displacement_history[index][wing_dofs + 2 * (propeller_index - 1) + 1]
            for index in 1:number_of_states
        ]
        for propeller_index in 1:number_of_propellers
    ]
    propeller_yaw = [
        [
            displacement_history[index][wing_dofs + 2 * (propeller_index - 1) + 2]
            for index in 1:number_of_states
        ]
        for propeller_index in 1:number_of_propellers
    ]

    history_path = joinpath(output_directory, output_label * "_history.csv")
    open(history_path, "w") do stream
        propeller_headers = String[]
        for propeller_index in 1:number_of_propellers
            push!(propeller_headers, "propeller_$(propeller_index)_pitch_deg")
            push!(propeller_headers, "propeller_$(propeller_index)_yaw_deg")
        end
        println(stream, join([
            "time_s",
            "tip_displacement_m",
            "tip_twist_deg",
            propeller_headers...,
            "coupling_iterations",
            "coupling_state_residual",
            "coupling_load_residual",
            "coupling_equilibrium_residual",
            "coupling_converged",
        ], ","))

        for state_index in 1:number_of_states
            step_index = clamp(state_index - 1, 1, max(active_step_count, 1))
            propeller_values = Float64[]
            for propeller_index in 1:number_of_propellers
                push!(propeller_values, rad2deg(propeller_pitch[propeller_index][state_index]))
                push!(propeller_values, rad2deg(propeller_yaw[propeller_index][state_index]))
            end
            values = Any[
                time_history[state_index],
                tip_displacement[state_index],
                rad2deg(tip_twist[state_index]),
                propeller_values...,
                state_index == 1 ? 0 : coupling_iterations[step_index],
                state_index == 1 ? 0.0 : coupling_state_residual[step_index],
                state_index == 1 ? 0.0 : coupling_load_residual[step_index],
                state_index == 1 ? 0.0 : coupling_equilibrium_residual[step_index],
                state_index == 1 ? true : coupling_converged[step_index],
            ]
            println(stream, join(values, ","))
        end
    end

    history_is_finite = all(
        all(isfinite, displacement_history[index]) for index in 1:number_of_states
    )
    completed_all_steps = active_step_count == length(time_steps)
    all_steps_converged = active_step_count > 0 &&
        all(coupling_converged[1:active_step_count])
    validation_passed = history_is_finite && completed_all_steps && all_steps_converged

    summary_path = joinpath(output_directory, output_label * "_summary.txt")
    open(summary_path, "w") do stream
        println(stream, "Chang linear aeroelastic / UVLM validation")
        println(stream, "source_case = run_chang_partitioned_ga_equilibrium_residual.jl")
        println(stream, "air_density_kg_m3 = $density")
        println(stream, "freestream_speed_m_s = $freestream_speed")
        println(stream, "number_of_propellers = $number_of_propellers")
        println(stream, "number_of_blades_per_propeller = $number_of_blades")
        println(stream, "interaction_on = $interaction_on")
        println(stream, "near_field_force_model = $near_field_force_model")
        println(stream, "propeller_moment_projection = $propeller_moment_projection")
        println(stream, "requested_end_time_s = $requested_end_time")
        println(stream, "integrated_end_time_s = $(time_history[end])")
        println(stream, "integrated_steps = $active_step_count")
        println(stream, "completed_all_steps = $completed_all_steps")
        println(stream, "history_is_finite = $history_is_finite")
        println(stream, "all_coupling_steps_converged = $all_steps_converged")
        println(stream, "validation_passed = $validation_passed")
        println(stream, "max_coupling_iterations = $(maximum(coupling_iterations[active_steps]))")
        println(stream, "max_coupling_state_residual = $(finite_maximum(coupling_state_residual[active_steps]))")
        println(stream, "max_coupling_load_residual = $(finite_maximum(coupling_load_residual[active_steps]))")
        println(stream, "max_coupling_equilibrium_residual = $(finite_maximum(coupling_equilibrium_residual[active_steps]))")
        println(stream, "max_abs_tip_displacement_m = $(maximum(abs, tip_displacement))")
        println(stream, "max_abs_tip_twist_deg = $(maximum(abs, rad2deg.(tip_twist)))")
        for propeller_index in 1:number_of_propellers
            println(stream, "max_abs_propeller_$(propeller_index)_pitch_deg = $(maximum(abs, rad2deg.(propeller_pitch[propeller_index])))")
            println(stream, "max_abs_propeller_$(propeller_index)_yaw_deg = $(maximum(abs, rad2deg.(propeller_yaw[propeller_index])))")
        end
    end

    println("Validation history written to $history_path")
    println("Validation summary written to $summary_path")
    println("VALIDATION_PASSED=$validation_passed")

    plot_path = nothing
    if plot_results
        plot_path = joinpath(output_directory, output_label * "_history.png")
        plot_chang_time_histories(
            time_history,
            tip_displacement,
            tip_twist,
            propeller_pitch,
            propeller_yaw,
            propeller_eta;
            span_length,
            time_limit_s = plot_time_limit,
            output_path = plot_path,
        )
    end

    return (;
        validation_passed,
        history_path,
        summary_path,
        plot_path,
        time = time_history,
        tip_displacement,
        tip_twist,
        propeller_pitch,
        propeller_yaw,
    )
end
