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
    perturbation_duration, animation_frames, progress_frequency)

    @assert airspeed > 0 && density > 0 && duration > 0
    @assert chordwise_panels > 0 && spanwise_panels > 0 && maximum_wake_rows > 0
    @assert 0 <= settling_time < duration && perturbation_duration > 0
    @assert 0 < initial_airspeed_fraction <= 1
    @assert 0 <= airspeed_ramp_duration <= settling_time
    @assert animation_frames >= 2
    @assert progress_frequency >= 1
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

    # Save the tip motion at every step and a few wing/wake animation frames.
    # Half of the frames are assigned to the initial wake growth (up to the
    # retained-row limit); the remainder follows the rest of the simulation.
    # A uniform frame stride would miss this short startup for long runs.
    tip_out_of_plane = zeros(length(time))
    tip_twist_degrees = zeros(length(time))
    airspeed_history = zeros(length(time))
    aerodynamic_nodal_load_history = zeros(6, n_elements + 1, length(time))
    airspeed_history[1] = initial_airspeed
    surface_history, wake_history, animation_time = Any[], Any[], Float64[]
    save_wing_frame!(aerodynamic, surface_history, wake_history, animation_time)
    number_of_animation_frames = min(animation_frames, n_steps + 1)
    wake_fill_step = min(maximum_wake_rows, n_steps)
    if wake_fill_step == n_steps
        animation_steps = round.(Int,
            range(0, n_steps; length=number_of_animation_frames))
    else
        early_frame_count = min(wake_fill_step + 1,
            cld(number_of_animation_frames, 2))
        late_frame_count = min(number_of_animation_frames - early_frame_count,
            n_steps - wake_fill_step)
        early_frame_count = number_of_animation_frames - late_frame_count
        early_steps = early_frame_count == 1 ? [0] : round.(Int,
            range(0, wake_fill_step; length=early_frame_count))
        late_steps = late_frame_count == 1 ? [n_steps] : round.(Int,
            range(wake_fill_step + 1, n_steps; length=late_frame_count))
        animation_steps = vcat(early_steps, late_steps)
    end
    animation_steps = Set(animation_steps)
    println("UVLM mesh: $chordwise_panels x $spanwise_panels panels")
    println("dt = $dt s; final time = $(last(time)) s; steps = $n_steps")
    println("Wake capacity: $maximum_wake_rows rows, nominal length = " *
        "$(maximum_wake_rows/chordwise_panels) chords ($(maximum_wake_rows*airspeed*dt) m)")
    println("Retained wake age at capacity = $(maximum_wake_rows*dt) s")
    println("Airspeed ramp: $initial_airspeed -> $airspeed m/s over " *
        "$airspeed_ramp_duration s; tip pulse starts at $settling_time s")

    # 4. Coupling loop: aerodynamic loads -> beam motion -> new wing and wake.
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

        # Extrapolate the wing's last motion to the new time. Reusing grid here
        # would tell UVLM that the wing velocity is zero at every structural step.
        predicted_grid = 2 .* grid .- previous_grid
        predicted_positions = 2 .* positions .- previous_positions
        UVLM.begin_time_step!(aerodynamic)
        loads = UVLM.evaluate_trial!(aerodynamic; surfaces=[predicted_grid])
        nodal_loads .= beam_loads(loads, predicted_positions, weights)
        aerodynamic_nodal_load_history[:, :, i] .= nodal_loads

        # After the settling interval, apply a short, smooth tip force [N].
        phase = (time[i] - settling_time) / perturbation_duration
        if 0 < phase < 1
            nodal_loads[1, end] += perturbation_amplitude * sin(2pi*phase) * sinpi(phase)^2
        end

        # Advance the nonlinear structure with these nodal forces and moments.
        for bc in model.BCs
            AeroBeams.update_BC_data!(bc, time[i])
        end
        AeroBeams.get_equivalent_states_rates!(structure)
        AeroBeams.solve_time_step!(structure)
        structure.systemSolver.convergedFinalSolution || error("AeroBeams failed at t=$(time[i])")

        # Update the aerodynamic geometry, recompute circulation, shed one wake row.
        previous_grid, previous_positions = grid, positions
        grid, positions = wing_geometry(model, chord, spar_fraction, chordwise_panels, weights)
        UVLM.evaluate_trial!(aerodynamic; surfaces=[grid])
        UVLM.commit_time_step!(aerodynamic)

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
                "U=$(round(current_airspeed; digits=3)) m/s, NR=converged")
        end
    end

    return (structural_problem=structure, aerodynamic_problem=aerodynamic,
        time=time, dt=dt, airspeed_history=airspeed_history,
        tip_out_of_plane=tip_out_of_plane,
        tip_bending_displacement=tip_out_of_plane,
        tip_twist_degrees=tip_twist_degrees,
        aerodynamic_nodal_load_history=aerodynamic_nodal_load_history,
        aerodynamic_surface_history=surface_history, wake_history=wake_history,
        animation_time=animation_time, animation_frames=animation_frames)
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
