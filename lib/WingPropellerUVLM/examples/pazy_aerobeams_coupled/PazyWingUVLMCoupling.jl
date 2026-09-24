using LinearAlgebra

# Load the two local packages regardless of the environment selected in the IDE.
const REPOSITORY_ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))
const UVLM_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
for package_root in (REPOSITORY_ROOT, UVLM_ROOT)
    package_root in LOAD_PATH || pushfirst!(LOAD_PATH, package_root)
end
import AeroBeams
import WingPropellerUVLM as UVLM

# Convert AeroBeams axes to UVLM axes: (x, y, z) -> (y, z, x).
const A_TO_UVLM = [0.0 1.0 0.0; 0.0 0.0 1.0; 1.0 0.0 0.0]

function run_pazy_wing_uvlm(; airspeed, density, angle_of_attack, sideslip,
    initial_airspeed_fraction, airspeed_ramp_duration,
    duration, settling_time, chordwise_panels, spanwise_panels, symmetric_wing,
    maximum_wake_rows, core_radius, wake_shedding_fraction, aerodynamic_interaction,
    save_uvlm_history, uvlm_save_frequency,
    newton_maximum_iterations, newton_absolute_tolerance,
    newton_relative_tolerance, newton_display_iterations,
    newton_always_update_jacobian, perturbation_amplitude,
    perturbation_duration, animation_frames, progress_frequency,
    animation_time_step=nothing,
    coupling_maximum_iterations=20, coupling_relaxation=0.3,
    coupling_geometry_tolerance=1e-5, coupling_load_tolerance=1e-3,
    coupling_display_iterations=false)

    @assert airspeed > 0 && density > 0 && duration > 0
    @assert chordwise_panels > 0 && spanwise_panels > 0 && maximum_wake_rows > 0
    @assert 0 <= settling_time < duration && perturbation_duration > 0
    @assert 0 < initial_airspeed_fraction <= 1
    @assert 0 <= airspeed_ramp_duration <= settling_time
    @assert animation_frames >= 2
    @assert isnothing(animation_time_step) || animation_time_step > 0
    @assert progress_frequency >= 1
    @assert coupling_maximum_iterations >= 2
    @assert 0 < coupling_relaxation <= 1
    @assert coupling_geometry_tolerance > 0 && coupling_load_tolerance > 0
    @assert !symmetric_wing || iszero(sideslip) "Nonzero sideslip requires symmetric_wing=false"
    # The UVLM constructor checks the radius, shedding fraction and save frequency.
    aero_solver = UVLM.create_UVLMSolver(
        coreRadius=core_radius, wakeSheddingFraction=wake_shedding_fraction,
        interaction=aerodynamic_interaction)

    # 1. Create the nonlinear Pazy beam, clamped at its root.
    n_elements, span, chord, spar_fraction = AeroBeams.geometrical_properties_Pazy()
    beam_nodes = AeroBeams.nodal_positions_Pazy()
    stiffness = AeroBeams.stiffness_matrices_Pazy(
        withSkin=true, sweepStructuralCorrections=false,
        GAy=1e16, GAz=1e16, Λ=0.0)
    inertia = AeroBeams.inertia_matrices_Pazy(withSkin=true)
    beam = AeroBeams.create_Beam(
        name="Pazy wing", length=span, nElements=n_elements,
        normalizedNodalPositions=beam_nodes, S=stiffness, I=inertia,
        rotationParametrization="E231", p0=[-pi/2; 0.0; 0.0],
        aeroSurface=nothing) # All aerodynamic forces come from UVLM.
    root_clamp = AeroBeams.create_BC(
        name="root clamp", beam=beam, node=1,
        types=["u1A", "u2A", "u3A", "p1A", "p2A", "p3A"], values=zeros(6))

    # Each column contains [Fx, Fy, Fz, Mx, My, Mz] at one beam node.
    # AeroBeams reads this array through the load functions at every step.
    nodal_loads = zeros(6, n_elements + 1)
    boundary_conditions = [root_clamp]
    for node in 1:n_elements+1
        load_functions = [t -> nodal_loads[k, node] for k in 1:6]
        push!(boundary_conditions, AeroBeams.create_BC(
            name="UVLM loads $node", beam=beam, node=node,
            types=["F1A", "F2A", "F3A", "M1A", "M2A", "M3A"],
            values=load_functions))
    end
    model = AeroBeams.create_Model(
        name="Pazy with UVLM", beams=[beam], BCs=boundary_conditions,
        gravityVector=zeros(3), v_A=t -> zeros(3),
        units=AeroBeams.create_UnitsSystem(frequency="Hz"))

    # 2. Choose dt so U*dt = chord / number of chordwise panels.
    dt = chord / (chordwise_panels * airspeed)
    n_steps = round(Int, duration / dt)
    @assert n_steps >= 1 "Duration must include at least one timestep"
    time = collect(0:n_steps) .* dt
    newton = AeroBeams.create_NewtonRaphson(
        maximumIterations=newton_maximum_iterations,
        absoluteTolerance=newton_absolute_tolerance,
        relativeTolerance=newton_relative_tolerance,
        displayStatus=newton_display_iterations,
        alwaysUpdateJacobian=newton_always_update_jacobian,
        minConvRateJacUpdate=1.2)
    structure = AeroBeams.create_DynamicProblem(
        model=model, timeVector=time, systemSolver=newton,
        trackingTimeSteps=true, trackingFrequency=1,
        displayProgress=true, saveInitialSolution=true)
    AeroBeams.precompute_distributed_loads!(structure)
    AeroBeams.solve_initial_dynamic!(structure)
    AeroBeams.save_time_step_data!(structure, 0.0)

    # 3. Place uniform UVLM span stations between the original beam nodes.
    # Each row of weights has just two nonzero entries, which sum to one.
    weights = zeros(spanwise_panels + 1, n_elements + 1)
    for j in 1:spanwise_panels+1
        eta = (j - 1) / spanwise_panels
        left = min(searchsortedlast(beam_nodes, eta), n_elements)
        fraction = (eta - beam_nodes[left]) / (beam_nodes[left+1] - beam_nodes[left])
        weights[j, left] = 1 - fraction
        weights[j, left+1] = fraction
    end
    grid, positions = wing_geometry(model, chord, spar_fraction, chordwise_panels, weights)
    previous_grid, previous_positions = copy(grid), copy(positions)
    initial_airspeed = airspeed * initial_airspeed_fraction
    reference = UVLM.Reference(span*chord, chord, span, zeros(3), airspeed, density)
    aerodynamic_model = UVLM.create_UVLMModel(
        surfaces=[grid], reference=reference, name="Pazy UVLM",
        symmetric=symmetric_wing, interactionGroups=[1])
    operating_point = UVLM.create_OperatingPoint(
        airspeed=initial_airspeed, density=density, angleOfAttack=angle_of_attack,
        sideslip=sideslip)
    aerodynamic = UVLM.create_UVLMDynamicProblem(
        model=aerodynamic_model, aeroSolver=aero_solver,
        operatingPoint=operating_point, timeVector=time,
        maximumWakeRows=maximum_wake_rows, trackingTimeSteps=save_uvlm_history,
        trackingFrequency=uvlm_save_frequency)

    # Save the tip motion at every step and uniformly spaced wing/wake frames.
    # This is the same stride used by AeroBeams' structural animation, so frame
    # k in both GIFs represents the same physical time and deformation.
    tip_out_of_plane = zeros(length(time))
    tip_twist_degrees = zeros(length(time))
    airspeed_history = zeros(length(time))
    aerodynamic_nodal_load_history = zeros(6, n_elements + 1, length(time))
    coupling_iterations_history = zeros(Int, length(time))
    coupling_geometry_residual_history = zeros(length(time))
    coupling_load_residual_history = zeros(length(time))
    airspeed_history[1] = initial_airspeed
    surface_history, wake_history, animation_time = Any[], Any[], Float64[]
    save_wing_frame!(aerodynamic, surface_history, wake_history, animation_time)
    animation_stride = isnothing(animation_time_step) ?
        max(1, cld(n_steps + 1, animation_frames)) :
        max(1, round(Int, animation_time_step / dt))
    actual_animation_time_step = animation_stride * dt
    animation_steps = Set(0:animation_stride:n_steps)
    println("UVLM mesh: $chordwise_panels x $spanwise_panels panels")
    println("dt = $dt s; final time = $(last(time)) s; steps = $n_steps")
    println("Wake capacity: $maximum_wake_rows rows, nominal length = " *
        "$(maximum_wake_rows/chordwise_panels) chords ($(maximum_wake_rows*airspeed*dt) m)")
    println("Retained wake age at capacity = $(maximum_wake_rows*dt) s")
    println("Airspeed ramp: $initial_airspeed -> $airspeed m/s over " *
        "$airspeed_ramp_duration s; tip pulse starts at $settling_time s")
    println("Strong coupling: maximum $coupling_maximum_iterations iterations, " *
        "relaxation = $coupling_relaxation")
    println("Animation sampling interval = $actual_animation_time_step s")

    # 4. Strong partitioned coupling loop. At each physical time step, UVLM
    # and AeroBeams are iterated to a common interface geometry and load. UVLM
    # trials always restart from the accepted wake at time n. The wake is
    # advanced exactly once, after the coupled iteration has converged.
    for i in 2:length(time)
        AeroBeams.update_time_variables!(structure, i)
        AeroBeams.update_basis_A_orientation!(structure)

        # Cubic smoothstep ramp: its slope is zero at the beginning and end.
        ramp_fraction = airspeed_ramp_duration == 0 ? 1.0 :
            clamp(time[i] / airspeed_ramp_duration, 0.0, 1.0)
        smooth_ramp = ramp_fraction^2 * (3 - 2*ramp_fraction)
        current_airspeed = airspeed * (initial_airspeed_fraction +
            (1 - initial_airspeed_fraction) * smooth_ramp)
        airspeed_history[i] = current_airspeed

        # UVLM uses this speed in both its freestream and dimensional loads.
        aerodynamic.operatingPoint = UVLM.create_OperatingPoint(
            airspeed=current_airspeed, density=density,
            angleOfAttack=angle_of_attack, sideslip=sideslip)
        aerodynamic.state.system.reference[] = UVLM.Reference(
            span*chord, chord, span, zeros(3), current_airspeed, density)

        # The extrapolated geometry is the first n+1 interface guess. Reusing
        # the accepted grid would incorrectly impose zero wing velocity on the
        # first UVLM trial.
        grid_guess = 2 .* grid .- previous_grid
        positions_guess = 2 .* positions .- previous_positions

        # The external pulse is constant throughout all coupling iterations at
        # this physical time. It is not part of the aerodynamic load residual.
        phase = (time[i] - settling_time) / perturbation_duration
        tip_force = 0 < phase < 1 ?
            perturbation_amplitude * sin(2pi*phase) * sinpi(phase)^2 : 0.0

        # The equivalent rates contain the accepted structural state at time n
        # and must be formed only once. Every inner solve is for the same n+1.
        AeroBeams.get_equivalent_states_rates!(structure)
        UVLM.begin_time_step!(aerodynamic)

        previous_aerodynamic_loads = nothing
        accepted_aerodynamic_loads = nothing
        candidate_grid, candidate_positions = grid_guess, positions_guess
        geometry_residual = Inf
        load_residual = Inf
        coupling_iterations = 0
        coupling_converged = false

        try
            for coupling_iteration in 1:coupling_maximum_iterations
                coupling_iterations = coupling_iteration

                # Aerodynamic trial at the current relaxed interface guess.
                loads = UVLM.evaluate_trial!(aerodynamic; surfaces=[grid_guess])
                aerodynamic_loads = beam_loads(loads, positions_guess, weights)

                # Apply aerodynamic loads and the prescribed pulse, then solve
                # the full nonlinear AeroBeams step with its Newton iterations.
                nodal_loads .= aerodynamic_loads
                nodal_loads[1, end] += tip_force
                for bc in model.BCs
                    AeroBeams.update_BC_data!(bc, time[i])
                end
                AeroBeams.solve_time_step!(structure)
                structure.systemSolver.convergedFinalSolution ||
                    error("AeroBeams failed at t=$(time[i]), coupling iteration $coupling_iteration")

                candidate_grid, candidate_positions = wing_geometry(
                    model, chord, spar_fraction, chordwise_panels, weights)
                geometry_residual = maximum(abs, candidate_grid .- grid_guess) / chord
                load_residual = isnothing(previous_aerodynamic_loads) ? Inf :
                    interface_load_residual(
                        aerodynamic_loads, previous_aerodynamic_loads, chord)

                coupling_display_iterations && println(
                    "  FSI $coupling_iteration: r_geometry=$(geometry_residual), " *
                    "r_load=$(load_residual)")

                if geometry_residual <= coupling_geometry_tolerance &&
                    load_residual <= coupling_load_tolerance
                    # Make the uncommitted UVLM trial exactly match the accepted
                    # structural geometry. This still does not advance the wake.
                    final_loads = UVLM.evaluate_trial!(
                        aerodynamic; surfaces=[candidate_grid])
                    final_aerodynamic_loads = beam_loads(
                        final_loads, candidate_positions, weights)
                    final_load_residual = interface_load_residual(
                        final_aerodynamic_loads, aerodynamic_loads, chord)
                    if final_load_residual <= coupling_load_tolerance
                        accepted_aerodynamic_loads = final_aerodynamic_loads
                        load_residual = max(load_residual, final_load_residual)
                        coupling_converged = true
                        break
                    end
                    load_residual = max(load_residual, final_load_residual)
                end

                previous_aerodynamic_loads = copy(aerodynamic_loads)
                grid_guess .= (1 - coupling_relaxation) .* grid_guess .+
                    coupling_relaxation .* candidate_grid
                positions_guess .= (1 - coupling_relaxation) .* positions_guess .+
                    coupling_relaxation .* candidate_positions
            end

            coupling_converged || error(
                "Strong coupling failed at t=$(time[i]) after " *
                "$coupling_iterations iterations: r_geometry=$geometry_residual, " *
                "r_load=$load_residual")

            # Leave the load BCs consistent with the accepted aerodynamic trial.
            nodal_loads .= accepted_aerodynamic_loads
            aerodynamic_nodal_load_history[:, :, i] .= accepted_aerodynamic_loads
            nodal_loads[1, end] += tip_force

            # Accept the structural geometry and advance the wake exactly once.
            previous_grid, previous_positions = grid, positions
            grid, positions = candidate_grid, candidate_positions
            UVLM.commit_time_step!(aerodynamic)
        catch
            UVLM.rollback_time_step!(aerodynamic)
            rethrow()
        end

        coupling_iterations_history[i] = coupling_iterations
        coupling_geometry_residual_history[i] = geometry_residual
        coupling_load_residual_history[i] = load_residual

        AeroBeams.save_time_step_data!(structure, time[i])
        tip_out_of_plane[i] = -model.elements[end].nodalStates.u_n2[1]
        tip_twist_degrees[i] = wingtip_twist_degrees(model)
        if (i-1) in animation_steps
            save_wing_frame!(aerodynamic, surface_history, wake_history, animation_time)
        end
        if (i-1) % progress_frequency == 0 || i == length(time)
            percent = round(100 * (i-1) / n_steps; digits=1)
            println("Step $(i-1)/$n_steps ($percent%), " *
                "t=$(round(time[i]; digits=4)) s, " *
                "U=$(round(current_airspeed; digits=3)) m/s, " *
                "FSI=$coupling_iterations, " *
                "r_x=$(round(geometry_residual; sigdigits=3)), " *
                "r_F=$(round(load_residual; sigdigits=3))")
        end
    end

    return (structural_problem=structure, aerodynamic_problem=aerodynamic,
        time=time, dt=dt, airspeed_history=airspeed_history,
        tip_out_of_plane=tip_out_of_plane,
        tip_bending_displacement=tip_out_of_plane,
        tip_twist_degrees=tip_twist_degrees,
        aerodynamic_nodal_load_history=aerodynamic_nodal_load_history,
        coupling_iterations=coupling_iterations_history,
        coupling_geometry_residual=coupling_geometry_residual_history,
        coupling_load_residual=coupling_load_residual_history,
        aerodynamic_surface_history=surface_history, wake_history=wake_history,
        animation_time=animation_time, animation_frames=animation_frames,
        animation_stride=animation_stride,
        animation_time_step=actual_animation_time_step)
end

# Relative interface-load change. Forces and moments are normalized separately;
# moments use force*chord as their minimum physical scale.
function interface_load_residual(current, previous, chord)
    force_scale = max(
        maximum(abs, current[1:3, :]),
        maximum(abs, previous[1:3, :]),
        1.0,
    )
    moment_scale = max(
        maximum(abs, current[4:6, :]),
        maximum(abs, previous[4:6, :]),
        force_scale * chord,
    )
    force_residual = maximum(abs, current[1:3, :] .- previous[1:3, :]) /
        force_scale
    moment_residual = maximum(abs, current[4:6, :] .- previous[4:6, :]) /
        moment_scale
    return max(force_residual, moment_residual)
end

# Same wingtip-twist definition used by the AeroBeams Pazy examples: rotate
# the local chord direction and measure its angle out of the wing plane.
function wingtip_twist_degrees(model)
    tip_rotation_local = model.elements[end].nodalStates.p_n2_b
    rotation, _ = AeroBeams.rotation_tensor_WM(tip_rotation_local)
    chord_direction = rotation * [0.0; 1.0; 0.0]
    return asind(clamp(chord_direction[3], -1.0, 1.0))
end

# Beam deformation -> aerodynamic grid. The chord rotates with each beam node.
function wing_geometry(model, chord, spar_fraction, chordwise_panels, weights)
    n_nodes = length(model.elements) + 1
    positions = zeros(3, n_nodes)
    beam_grid = zeros(3, chordwise_panels + 1, n_nodes)
    for node in 1:n_nodes
        if node == 1
            element = model.elements[1]
            # This root is fully clamped. Enforce its prescribed values exactly:
            # recovered beam outputs can have small residuals, which move the
            # root off y=0 and make its mirrored UVLM vortex nearly coincident.
            positions[:, node] = element.r_n1
            rotation = zeros(3)
            reference_rotation = element.R0_n1
        else
            element = model.elements[node-1]
            positions[:, node] = element.r_n2 + element.nodalStates.u_n2
            rotation = element.nodalStates.p_n2
            reference_rotation = element.R0_n2
        end
        R, _ = AeroBeams.rotation_tensor_WM(rotation)
        for k in 1:chordwise_panels+1
            offset = reference_rotation * [0.0; 1.0; 0.0] *
                (((k-1) / chordwise_panels - spar_fraction) * chord)
            beam_grid[:, k, node] = A_TO_UVLM * (positions[:, node] + R * offset)
        end
    end
    grid = zeros(3, chordwise_panels + 1, size(weights, 1))
    for j in axes(weights, 1), node in axes(weights, 2)
        grid[:, :, j] .+= weights[j, node] .* beam_grid[:, :, node]
    end
    return grid, positions
end

# Aerodynamic forces -> beam nodal forces and moments (moment = arm x force).
function beam_loads(loads, positions, weights)
    nodal_loads = zeros(6, size(positions, 2))
    forces, points = loads.forces[1], loads.positions[1]
    for j in axes(forces, 2), k in axes(forces, 1)
        force = A_TO_UVLM' * forces[k, j]
        point = A_TO_UVLM' * points[k, j]
        for node in axes(weights, 2)
            weight = weights[j, node]
            nodal_loads[1:3, node] .+= weight .* force
            nodal_loads[4:6, node] .+= weight .* cross(point - positions[:, node], force)
        end
    end
    return nodal_loads
end

# Copy the current geometry, because the solver changes it at the next step.
function save_wing_frame!(aerodynamic, surfaces, wakes, times)
    system = aerodynamic.state.system
    push!(surfaces, deepcopy(system.surfaces))
    push!(wakes, [deepcopy(system.wakes[1][1:system.nwake[1], :])])
    push!(times, aerodynamic.state.timeNow)
end
