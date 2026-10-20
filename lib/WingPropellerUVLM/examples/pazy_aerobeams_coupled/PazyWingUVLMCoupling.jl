using LinearAlgebra

# Load the two local packages regardless of the environment selected in the IDE.
const REPOSITORY_ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))
const UVLM_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
for package_root in (REPOSITORY_ROOT, UVLM_ROOT)
    package_root in LOAD_PATH || pushfirst!(LOAD_PATH, package_root)
end
import AeroBeams
import WingPropellerUVLM as UVLM

# Convert AeroBeams axes to the UVLM right-handed body axes:
#   UVLM +x (downstream) = AeroBeams -y
#   UVLM +y (right span) = AeroBeams +z
#   UVLM +z (up)         = AeroBeams -x
const A_TO_UVLM = [0.0 -1.0 0.0; 0.0 0.0 1.0; -1.0 0.0 0.0]

"""Run the Pazy coupling, retaining scalar responses at every accepted step.

`save_structural_history` and `save_aerodynamic_load_history` default to true.
Disable them for response-only sweeps; their result arrays are then empty.
`save_animation_history` controls the separate lattice/wake frame capture.
Structural animations require structural history; load glyphs also require
nodal load history. Wake-only animations need only lattice/wake frames.

Supply `airspeed_profile(t)` for a continuous simulation with prescribed speed
changes. It replaces the startup ramp (use fraction 1 and ramp duration 0).
`airspeed` still sets the fixed timestep; use the maximum prescribed speed to
keep U*dt/chord at or below `time_step_chords`. No state is reset at speed changes.

`initial_equilibrium` may be a converged `run_pazy_wing_uvlm_static` result.
Its operating point, mesh, wake capacity and timestep must match the starting
condition. The beam state and developed UVLM state are copied into a fresh run;
the supplied equilibrium is preserved. No unloaded startup solve is performed.
"""
function run_pazy_wing_uvlm(; airspeed, density, angle_of_attack, sideslip,
    initial_airspeed_fraction, airspeed_ramp_duration,
    duration, settling_time, chordwise_panels, spanwise_panels, symmetric_wing,
    maximum_wake_rows, core_radius, wake_shedding_fraction, aerodynamic_interaction,
    save_uvlm_history, uvlm_save_frequency,
    newton_maximum_iterations, newton_absolute_tolerance,
    newton_relative_tolerance, newton_display_iterations,
    newton_always_update_jacobian, perturbation_amplitude,
    perturbation_duration, animation_frames, progress_frequency,
    time_step_chords=1 / chordwise_panels,
    airspeed_profile=nothing,
    initial_equilibrium=nothing,
    animation_time_step=nothing,
    save_animation_history=true,
    save_structural_history=true,
    save_aerodynamic_load_history=true,
    coupling_scheme=:strong,
    coupling_maximum_iterations=20, coupling_relaxation=0.3,
    coupling_geometry_tolerance=1e-5, coupling_load_tolerance=1e-3,
    coupling_display_iterations=false)

    @assert airspeed > 0 && density > 0 && duration > 0
    @assert chordwise_panels > 0 && spanwise_panels > 0 && maximum_wake_rows > 0
    @assert isfinite(time_step_chords) && time_step_chords > 0
    @assert 0 <= settling_time < duration && perturbation_duration > 0
    @assert 0 < initial_airspeed_fraction <= 1
    @assert 0 <= airspeed_ramp_duration <= settling_time
    @assert isnothing(airspeed_profile) ||
        (initial_airspeed_fraction == 1 && iszero(airspeed_ramp_duration)) "A custom airspeed profile replaces the startup ramp"
    @assert animation_frames >= 2
    @assert isnothing(animation_time_step) || animation_time_step > 0
    @assert progress_frequency >= 1
    @assert coupling_scheme in (:loose, :strong) "Choose :loose or :strong coupling"
    @assert coupling_scheme == :loose || coupling_maximum_iterations >= 2
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
    # Copy these values into the active model's BCs before every structural solve.
    # Numeric BCs remain valid if Newton restores a deep copy of the model.
    nodal_loads = zeros(6, n_elements + 1)
    boundary_conditions = [root_clamp]
    for node in 1:n_elements+1
        push!(boundary_conditions, AeroBeams.create_BC(
            name="UVLM loads $node", beam=beam, node=node,
            types=["F1A", "F2A", "F3A", "M1A", "M2A", "M3A"],
            values=zeros(6)))
    end
    model = AeroBeams.create_Model(
        name="Pazy with UVLM", beams=[beam], BCs=boundary_conditions,
        gravityVector=zeros(3), v_A=t -> zeros(3),
        units=AeroBeams.create_UnitsSystem(frequency="Hz"))
    if !isnothing(initial_equilibrium)
        @assert initial_equilibrium.converged "Initial coupled static equilibrium must have converged"
        model = deepcopy(initial_equilibrium.structural_problem.model)
        @assert length(model.elements) == n_elements
    end

    # 2. Choose the physical timestep independently of the aerodynamic mesh.
    # The nondimensional step is U*dt/chord = time_step_chords.
    dt = time_step_chords * chord / airspeed
    n_steps = round(Int, duration / dt)
    @assert n_steps >= 1 "Duration must include at least one timestep"
    time = collect(0:n_steps) .* dt
    airspeed_history = if isnothing(airspeed_profile)
        map(time) do t
            ramp_fraction = airspeed_ramp_duration == 0 ? 1.0 :
                clamp(t / airspeed_ramp_duration, 0.0, 1.0)
            # Preserve the original initial sample when the ramp duration is zero.
            t == 0 && return airspeed * initial_airspeed_fraction
            smooth_ramp = ramp_fraction^2 * (3 - 2*ramp_fraction)
            airspeed * (initial_airspeed_fraction +
                (1 - initial_airspeed_fraction) * smooth_ramp)
        end
    else
        Float64[airspeed_profile(t) for t in time]
    end
    @assert all(u -> isfinite(u) && u > 0, airspeed_history) "Prescribed airspeeds must be finite and positive"
    if !isnothing(initial_equilibrium)
        point = initial_equilibrium.aerodynamic_problem.operatingPoint
        @assert isapprox(initial_equilibrium.dt, dt; rtol=1e-10, atol=0) "Equilibrium and dynamic timesteps must match to preserve wake spacing"
        @assert isapprox(point.airspeed, first(airspeed_history); rtol=1e-12)
        @assert point.density == density && point.angleOfAttack == angle_of_attack && point.sideslip == sideslip "Equilibrium operating point must match the initial condition"
        @assert initial_equilibrium.aerodynamic_problem.model.symmetric == [symmetric_wing]
    end
    newton = AeroBeams.create_NewtonRaphson(
        maximumIterations=newton_maximum_iterations,
        absoluteTolerance=newton_absolute_tolerance,
        relativeTolerance=newton_relative_tolerance,
        displayStatus=newton_display_iterations,
        alwaysUpdateJacobian=newton_always_update_jacobian,
        minConvRateJacUpdate=1.2)
    structure = AeroBeams.create_DynamicProblem(
        model=model, timeVector=time, systemSolver=newton,
        x0=isnothing(initial_equilibrium) ? Float64[] : copy(initial_equilibrium.structural_problem.x),
        trackingTimeSteps=save_structural_history, trackingFrequency=1,
        displayProgress=true, saveInitialSolution=true)
    AeroBeams.precompute_distributed_loads!(structure)
    isnothing(initial_equilibrium) && AeroBeams.solve_initial_dynamic!(structure)
    save_structural_history && AeroBeams.save_time_step_data!(structure, 0.0)

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
    span_map = PazySpanMap(weights)
    geometry_workspace = (positions=zeros(3,n_elements+1),
        beam_grid=zeros(3,chordwise_panels+1,n_elements+1))
    grid = zeros(3,chordwise_panels+1,spanwise_panels+1)
    vortex_offsets = zeros(3,chordwise_panels+1,n_elements+1)
    wing_geometry!(grid,vortex_offsets,structure.model,chord,spar_fraction,
        chordwise_panels,span_map,geometry_workspace)
    previous_grid, previous_offsets = copy(grid), copy(vortex_offsets)
    candidate_grid, candidate_offsets = similar(grid), similar(vortex_offsets)
    grid_guess, offsets_guess = similar(grid), similar(vortex_offsets)
    aerodynamic_loads = similar(nodal_loads)
    final_aerodynamic_loads = similar(nodal_loads)
    previous_aerodynamic_loads = similar(nodal_loads)
    initial_airspeed = first(airspeed_history)
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
    initial_equilibrium_state = nothing
    if !isnothing(initial_equilibrium)
        source = initial_equilibrium.aerodynamic_problem.state.system
        @assert length(source.grids) == 1 && size(source.grids[1]) == size(grid) "Equilibrium UVLM mesh must match"
        @assert size(source.wakes[1]) == size(aerodynamic.state.system.wakes[1]) "Equilibrium wake capacity must match"
        @assert maximum_absolute_difference(source.grids[1], grid) / chord <= 1e-10 "Equilibrium structural and aerodynamic geometry must match"
        @assert initial_equilibrium.aerodynamic_problem.aeroSolver == aero_solver "Equilibrium UVLM solver settings must match"
        UVLM.copy_system!(aerodynamic.state.system, source)
        # Fresh problem clocks remain at zero; only physical states are seeded.
        initial_equilibrium_state = (
            structural_state=copy(structure.x),
            aerodynamic_snapshot=UVLM.snapshot_uvlm(aerodynamic.state.system),
            applied_nodal_loads=copy(initial_equilibrium.applied_nodal_loads),
            maximum_velocity=maximum(norm(vcat(e.states.V, e.states.Ω)) for e in structure.model.elements))
    end

    # Save the tip motion at every step and uniformly spaced wing/wake frames.
    # This is the same stride used by AeroBeams' structural animation, so frame
    # k in both GIFs represents the same physical time and deformation.
    tip_out_of_plane = zeros(length(time))
    tip_twist_degrees = zeros(length(time))
    if !isnothing(initial_equilibrium)
        tip_out_of_plane[1] = -structure.model.elements[end].nodalStates.u_n2[1]
        tip_twist_degrees[1] = wingtip_twist_degrees(structure.model)
    end
    # Near-root section values (first element midpoint), retained even when
    # detailed structural histories are disabled. The driver reconstructs
    # surface gauges from these local beam section quantities.
    root_axial_strain = zeros(length(time))
    root_curvature_2 = zeros(length(time))
    root_curvature_3 = zeros(length(time))
    root_section = structure.model.elements[1].compStates
    root_axial_strain[1] = root_section.γ[1]
    root_curvature_2[1] = root_section.κ[2]
    root_curvature_3[1] = root_section.κ[3]
    aerodynamic_nodal_load_history = zeros(6,n_elements+1,save_aerodynamic_load_history ? length(time) : 0)
    if !isnothing(initial_equilibrium) && save_aerodynamic_load_history
        aerodynamic_nodal_load_history[:, :, 1] .= initial_equilibrium.aerodynamic_nodal_loads
    end
    coupling_iterations_history = zeros(Int, length(time))
    coupling_geometry_residual_history = zeros(length(time))
    coupling_load_residual_history = zeros(length(time))
    surface_history, wake_history, animation_time = Any[], Any[], Float64[]
    if save_animation_history
        save_wing_frame!(aerodynamic, surface_history, wake_history, animation_time)
    end
    animation_stride = isnothing(animation_time_step) ?
        max(1, cld(n_steps + 1, animation_frames)) :
        max(1, round(Int, animation_time_step / dt))
    actual_animation_time_step = animation_stride * dt
    animation_steps = save_animation_history ? Set(0:animation_stride:n_steps) : Set{Int}()
    println("UVLM mesh: $chordwise_panels x $spanwise_panels panels")
    println("dt = $dt s (U*dt/chord = $time_step_chords); " *
        "final time = $(last(time)) s; steps = $n_steps")
    nominal_wake_length_chords = maximum_wake_rows * airspeed * dt / chord
    println("Wake capacity: $maximum_wake_rows rows, nominal length = " *
        "$nominal_wake_length_chords chords ($(maximum_wake_rows*airspeed*dt) m)")
    println("Retained wake age at capacity = $(maximum_wake_rows*dt) s")
    if isnothing(airspeed_profile)
        println("Airspeed ramp: $initial_airspeed -> $airspeed m/s over " *
            "$airspeed_ramp_duration s; tip pulse starts at $settling_time s")
    else
        println("Prescribed airspeed: $(minimum(airspeed_history)) to " *
            "$(maximum(airspeed_history)) m/s; fixed dt; tip pulse starts at $settling_time s")
    end
    if coupling_scheme == :strong
        println("Coupling scheme: strong partitioned, maximum " *
            "$coupling_maximum_iterations iterations, relaxation = $coupling_relaxation")
    else
        println("Coupling scheme: loose (one aerodynamic/structural pass per timestep)")
    end
    if save_animation_history
        println("Animation sampling interval = $actual_animation_time_step s")
    else
        println("Animation frame capture disabled")
    end

    # 4. Partitioned coupling loop. Strong coupling iterates UVLM and
    # AeroBeams to a common interface state. Loose coupling evaluates UVLM once
    # on the extrapolated n+1 geometry and advances after one structural solve.
    # In both schemes, the wake advances exactly once per physical timestep.
    for i in 2:length(time)
        AeroBeams.update_time_variables!(structure, i)
        AeroBeams.update_basis_A_orientation!(structure;update_R_A_array=save_structural_history)

        current_airspeed = airspeed_history[i]

        # UVLM uses this speed in both its freestream and dimensional loads.
        aerodynamic.operatingPoint = UVLM.create_OperatingPoint(
            airspeed=current_airspeed, density=density,
            angleOfAttack=angle_of_attack, sideslip=sideslip)
        aerodynamic.state.system.reference[] = UVLM.Reference(
            span*chord, chord, span, zeros(3), current_airspeed, density)

        # The extrapolated geometry is the first n+1 interface guess. Reusing
        # the accepted grid would incorrectly impose zero wing velocity on the
        # first UVLM trial.
        grid_guess .= 2 .* grid .- previous_grid
        offsets_guess .= 2 .* vortex_offsets .- previous_offsets

        # The external pulse is constant throughout all coupling iterations at
        # this physical time. It is not part of the aerodynamic load residual.
        phase = (time[i] - settling_time) / perturbation_duration
        tip_force = 0 < phase < 1 ?
            perturbation_amplitude * sin(2pi*phase) * sinpi(phase)^2 : 0.0

        # The equivalent rates contain the accepted structural state at time n
        # and must be formed only once. Every inner solve is for the same n+1.
        AeroBeams.get_equivalent_states_rates!(structure)
        UVLM.begin_time_step!(aerodynamic)

        have_previous_loads = false
        accepted_aerodynamic_loads = nothing
        geometry_residual = Inf
        load_residual = Inf
        coupling_iterations = 0
        coupling_converged = false

        try
            if coupling_scheme == :loose
                coupling_iterations = 1

                # Explicit staggered pass: evaluate UVLM on the predicted n+1
                # geometry, then apply those loads in one AeroBeams solve. The
                # resulting interface mismatch is accepted without iteration.
                loads = UVLM.evaluate_trial!(aerodynamic; surfaces=[grid_guess])
                accepted_aerodynamic_loads = beam_loads!(
                    aerodynamic_loads,loads,offsets_guess,span_map)
                nodal_loads .= accepted_aerodynamic_loads
                nodal_loads[1, end] += tip_force
                apply_pazy_nodal_loads!(structure, nodal_loads, time[i])
                AeroBeams.solve_time_step!(structure)
                structure.systemSolver.convergedFinalSolution ||
                    error("AeroBeams failed at t=$(time[i]), loose coupling pass")

                wing_geometry!(candidate_grid,candidate_offsets,structure.model,chord,
                    spar_fraction,chordwise_panels,span_map,geometry_workspace)
                geometry_residual = maximum_absolute_difference(candidate_grid,grid_guess) / chord
                # A load convergence residual requires at least two FSI passes.
                load_residual = NaN
                coupling_converged = true

                coupling_display_iterations && println(
                    "  Loose FSI pass: r_geometry=$(geometry_residual), r_load=not applicable")
            else
                for coupling_iteration in 1:coupling_maximum_iterations
                    coupling_iterations = coupling_iteration

                    # Aerodynamic trial at the current relaxed interface guess.
                    loads = UVLM.evaluate_trial!(aerodynamic; surfaces=[grid_guess])
                    beam_loads!(aerodynamic_loads,loads,offsets_guess,span_map)

                    # Apply aerodynamic loads and the prescribed pulse, then solve
                    # the full nonlinear AeroBeams step with its Newton iterations.
                    nodal_loads .= aerodynamic_loads
                    nodal_loads[1, end] += tip_force
                    apply_pazy_nodal_loads!(structure, nodal_loads, time[i])
                    AeroBeams.solve_time_step!(structure)
                    structure.systemSolver.convergedFinalSolution ||
                        error("AeroBeams failed at t=$(time[i]), coupling iteration $coupling_iteration")

                    # Newton may replace structure.model during a retry. Always
                    # read the active model, not the original construction reference.
                    wing_geometry!(candidate_grid,candidate_offsets,structure.model,chord,
                        spar_fraction,chordwise_panels,span_map,geometry_workspace)
                    geometry_residual = maximum_absolute_difference(candidate_grid,grid_guess) / chord
                    load_residual = !have_previous_loads ? Inf :
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
                        beam_loads!(final_aerodynamic_loads,final_loads,candidate_offsets,span_map)
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

                    previous_aerodynamic_loads .= aerodynamic_loads
                    have_previous_loads = true
                    grid_guess .= (1 - coupling_relaxation) .* grid_guess .+
                        coupling_relaxation .* candidate_grid
                    # Keep the moment arms consistent with the relaxed trial grid.
                    offsets_guess .= (1 - coupling_relaxation) .* offsets_guess .+
                        coupling_relaxation .* candidate_offsets
                end

                coupling_converged || error(
                    "Strong coupling failed at t=$(time[i]) after " *
                    "$coupling_iterations iterations: r_geometry=$geometry_residual, " *
                    "r_load=$load_residual")
            end

            # Leave the load BCs consistent with the accepted aerodynamic trial.
            nodal_loads .= accepted_aerodynamic_loads
            if save_aerodynamic_load_history
                aerodynamic_nodal_load_history[:, :, i] .= accepted_aerodynamic_loads
            end
            nodal_loads[1, end] += tip_force
            apply_pazy_nodal_loads!(structure, nodal_loads, time[i])

            # Accept the structural geometry and advance the wake exactly once.
            previous_grid, grid, candidate_grid = grid, candidate_grid, previous_grid
            previous_offsets, vortex_offsets, candidate_offsets = vortex_offsets, candidate_offsets, previous_offsets
            UVLM.commit_time_step!(aerodynamic)
        catch
            UVLM.rollback_time_step!(aerodynamic)
            rethrow()
        end

        coupling_iterations_history[i] = coupling_iterations
        coupling_geometry_residual_history[i] = geometry_residual
        coupling_load_residual_history[i] = load_residual

        save_structural_history && AeroBeams.save_time_step_data!(structure, time[i])
        tip_out_of_plane[i] = -structure.model.elements[end].nodalStates.u_n2[1]
        tip_twist_degrees[i] = wingtip_twist_degrees(structure.model)
        root_section = structure.model.elements[1].compStates
        root_axial_strain[i] = root_section.γ[1]
        root_curvature_2[i] = root_section.κ[2]
        root_curvature_3[i] = root_section.κ[3]
        if save_animation_history && (i-1) in animation_steps
            save_wing_frame!(aerodynamic, surface_history, wake_history, animation_time)
        end
        if (i-1) % progress_frequency == 0 || i == length(time)
            percent = round(100 * (i-1) / n_steps; digits=1)
            load_residual_text = isfinite(load_residual) ?
                string(round(load_residual; sigdigits=3)) : "n/a"
            println("Step $(i-1)/$n_steps ($percent%), " *
                "t=$(round(time[i]; digits=4)) s, " *
                "U=$(round(current_airspeed; digits=3)) m/s, " *
                "FSI=$coupling_iterations, " *
                "r_x=$(round(geometry_residual; sigdigits=3)), " *
                "r_F=$load_residual_text")
        end
    end

    return (structural_problem=structure, aerodynamic_problem=aerodynamic,
        initial_equilibrium_state=initial_equilibrium_state,
        coupling_scheme=coupling_scheme,
        time=time, dt=dt, time_step_chords=time_step_chords,
        airspeed_history=airspeed_history,
        tip_out_of_plane=tip_out_of_plane,
        tip_bending_displacement=tip_out_of_plane,
        tip_twist_degrees=tip_twist_degrees,
        root_axial_strain=root_axial_strain,
        root_curvature_2=root_curvature_2, root_curvature_3=root_curvature_3,
        aerodynamic_nodal_load_history=aerodynamic_nodal_load_history,
        save_structural_history=save_structural_history,
        save_aerodynamic_load_history=save_aerodynamic_load_history,
        coupling_iterations=coupling_iterations_history,
        coupling_geometry_residual=coupling_geometry_residual_history,
        coupling_load_residual=coupling_load_residual_history,
        aerodynamic_surface_history=surface_history, wake_history=wake_history,
        animation_time=animation_time, animation_frames=animation_frames,
        animation_stride=animation_stride,
        animation_time_step=actual_animation_time_step)
end

# This example's BC order is: root clamp, then one six-component load per node.
# Updating numeric values on the ACTIVE model avoids copied closure/load arrays
# after a Newton retry. update_BC_data! also refreshes AeroBeams' cached values.
function apply_pazy_nodal_loads!(problem, nodal_loads, time)
    @assert size(nodal_loads, 1) == 6
    @assert length(problem.model.BCs) == size(nodal_loads, 2) + 1
    for node in axes(nodal_loads, 2)
        bc = problem.model.BCs[node+1]
        @views bc.values .= nodal_loads[:, node]
    end
    for bc in problem.model.BCs
        AeroBeams.update_BC_data!(bc, time)
    end
    return nothing
end

# Relative interface-load change. Forces and moments are normalized separately;
# moments use force*chord as their minimum physical scale.
function interface_load_residual(current, previous, chord)
    size(current) == size(previous) || throw(DimensionMismatch("Interface load shapes differ"))
    force_scale, moment_scale = 1.0, 0.0
    force_change, moment_change = 0.0, 0.0
    for node in axes(current,2), d in 1:3
        force_scale = max(force_scale,abs(current[d,node]),abs(previous[d,node]))
        moment_scale = max(moment_scale,abs(current[d+3,node]),abs(previous[d+3,node]))
        force_change = max(force_change,abs(current[d,node]-previous[d,node]))
        moment_change = max(moment_change,abs(current[d+3,node]-previous[d+3,node]))
    end
    return max(force_change/force_scale,moment_change/max(moment_scale,force_scale*chord))
end

maximum_absolute_difference(a,b) = maximum(abs(a[i]-b[i]) for i in eachindex(a,b))

"""Precomputed nonzero span support; linear stations have at most two nodes."""
struct PazySpanMap
    nodes::Vector{Vector{Int}}
    weights::Vector{Vector{Float64}}
    node_count::Int
end
function PazySpanMap(weights::AbstractMatrix)
    nodes = [findall(!iszero,view(weights,j,:)) for j in axes(weights,1)]
    values = [Float64[weights[j,node] for node in nodes[j]] for j in axes(weights,1)]
    return PazySpanMap(nodes,values,size(weights,2))
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
# Also return the beam-node-to-vortex-vertex offsets in AeroBeams axes, so the
# reverse load transfer uses the same motion mapping (preserving virtual work).
function wing_geometry(model, chord, spar_fraction, chordwise_panels, weights)
    n_nodes = length(model.elements) + 1
    workspace = (positions=zeros(3,n_nodes),beam_grid=zeros(3,chordwise_panels+1,n_nodes))
    offsets = zeros(3,chordwise_panels+1,n_nodes)
    grid = zeros(3,chordwise_panels+1,size(weights,1))
    wing_geometry!(grid,offsets,model,chord,spar_fraction,chordwise_panels,PazySpanMap(weights),workspace)
    return grid,workspace.positions,offsets
end

function wing_geometry!(grid,vortex_offsets,model,chord,spar_fraction,chordwise_panels,
    map::PazySpanMap,workspace)
    n_nodes = length(model.elements)+1
    positions, beam_grid = workspace.positions, workspace.beam_grid
    size(grid) == (3,chordwise_panels+1,length(map.nodes)) &&
        size(vortex_offsets) == size(beam_grid) == (3,chordwise_panels+1,n_nodes) &&
        size(positions) == (3,n_nodes) && map.node_count == n_nodes ||
        throw(DimensionMismatch("Pazy geometry workspace does not match mesh"))
    for node in 1:n_nodes
        if node == 1
            element = model.elements[1]
            # This root is fully clamped. Enforce its prescribed values exactly:
            # recovered beam outputs can have small residuals, which move the
            # root off y=0 and make its mirrored UVLM vortex nearly coincident.
            positions[:, node] = element.r_n1
            rotation = UVLM.SVector(0.0,0.0,0.0)
            reference_rotation = element.R0_n1
        else
            element = model.elements[node-1]
            for d in 1:3
                positions[d,node] = element.r_n2[d]+element.nodalStates.u_n2[d]
            end
            rotation = element.nodalStates.p_n2
            reference_rotation = element.R0_n2
        end
        R, _ = AeroBeams.rotation_tensor_WM(rotation)
        # AeroBeams +local-y points forward, from the spar toward the leading
        # edge. UVLM chord fractions instead increase from leading to trailing.
        chord_direction = UVLM.SMatrix{3,3}(R) *
            UVLM.SVector(reference_rotation[1,2],reference_rotation[2,2],reference_rotation[3,2])
        for k in 1:chordwise_panels+1
            chord_fraction = (k-1) / chordwise_panels
            offset = chord_direction *
                ((spar_fraction - chord_fraction) * chord)
            beam_grid[1,k,node] = -(positions[2,node]+offset[2])
            beam_grid[2,k,node] = positions[3,node]+offset[3]
            beam_grid[3,k,node] = -(positions[1,node]+offset[1])
            # UVLM places each vortex row 1/4 panel aft of its grid row,
            # except the last row, which stays at the physical trailing edge.
            vortex_fraction = k <= chordwise_panels ?
                (k-1+0.25) / chordwise_panels : 1.0
            distance = (spar_fraction-vortex_fraction)*chord
            for d in 1:3
                vortex_offsets[d,k,node] = chord_direction[d]*distance
            end
        end
    end
    fill!(grid,0)
    for j in eachindex(map.nodes), entry in eachindex(map.nodes[j])
        node, weight = map.nodes[j][entry], map.weights[j][entry]
        for k in axes(grid,2), d in 1:3
            grid[d,k,j] += weight * beam_grid[d,k,node]
        end
    end
    return grid, positions, vortex_offsets
end

# Reverse of wing_geometry: x_a = sum(weight * (node_position + node_offset)).
# F_node = sum(weight * F_a); M_node = sum(weight * node_offset x F_a).
# Using the interpolated aero point minus the beam position instead would add
# spurious spanwise lever arms and would not preserve virtual work.
function beam_loads(loads, vortex_offsets, weights)
    nodal_loads = zeros(6, size(weights, 2))
    return beam_loads!(nodal_loads,loads,vortex_offsets,PazySpanMap(weights))
end

function beam_loads!(nodal_loads,loads,vortex_offsets,map::PazySpanMap)
    forces = loads.forces[1]
    @assert size(vortex_offsets) == (3,size(forces,1),map.node_count)
    @assert size(forces,2) == length(map.nodes)
    @assert size(nodal_loads) == (6,map.node_count)
    fill!(nodal_loads,0)
    for j in axes(forces, 2), k in axes(forces, 1)
        f = forces[k,j]
        force = UVLM.SVector(-f[3],-f[1],f[2])
        for entry in eachindex(map.nodes[j])
            node, weight = map.nodes[j][entry], map.weights[j][entry]
            offset = UVLM.SVector(vortex_offsets[1,k,node],vortex_offsets[2,k,node],vortex_offsets[3,k,node])
            moment = cross(offset,force)
            for d in 1:3
                nodal_loads[d,node] += weight*force[d]
                nodal_loads[d+3,node] += weight*moment[d]
            end
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
