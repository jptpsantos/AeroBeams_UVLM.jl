# Nonlinear STATIC beam balance with time-marched UVLM aerodynamics.
# Reuse the verified coordinate, geometry, and load-transfer functions.
if !isdefined(@__MODULE__, :run_pazy_wing_uvlm)
    include(joinpath(@__DIR__, "PazyWingUVLMCoupling.jl"))
end

# Hold the wing fixed and let its free wake develop. Each new static trial
# starts with an empty wake: a shape iteration is NOT a physical wing velocity.
function settle_static_uvlm(grid, offsets, weights, reference, operating_point, solver;
    symmetric_wing, maximum_wake_rows, dt, maximum_time, minimum_time,
    window_time, force_tolerance, circulation_tolerance, progress_frequency)

    steps = max(1, floor(Int, maximum_time / dt))
    window_steps = max(2, ceil(Int, window_time / dt) + 1)
    # Do not accept the initial wake transient as equilibrium.
    minimum_steps = max(ceil(Int, minimum_time / dt), 2maximum_wake_rows)
    model = UVLM.create_UVLMModel(surfaces=[grid], reference=reference,
        symmetric=symmetric_wing, interactionGroups=[1])
    aero = UVLM.create_UVLMDynamicProblem(model=model, aeroSolver=solver,
        operatingPoint=operating_point, timeVector=collect(0:steps).*dt,
        maximumWakeRows=maximum_wake_rows, trackingTimeSteps=false)
    force_window = Matrix{Float64}[]
    circulation_window = Vector{Float64}[]
    nodal_loads = zeros(6, size(weights,2))
    force_residual = circulation_residual = Inf
    converged = false
    for step in 1:steps
        UVLM.begin_time_step!(aero)
        try
            loads = UVLM.evaluate_trial!(aero; surfaces=[grid])
            nodal_loads = beam_loads(loads, offsets, weights)
            UVLM.commit_time_step!(aero)
        catch
            UVLM.rollback_time_step!(aero)
            rethrow()
        end
        push!(force_window, copy(nodal_loads))
        push!(circulation_window, copy(aero.state.system.Γ))
        if length(force_window) > window_steps
            popfirst!(force_window)
            popfirst!(circulation_window)
        end
        if length(force_window) == window_steps
            # Peak-to-peak change across the WHOLE window, not just two steps.
            fmax = reduce((a,b) -> max.(a,b), force_window)
            fmin = reduce((a,b) -> min.(a,b), force_window)
            force_residual = interface_load_residual(fmax, fmin, reference.c)
            gmax = reduce((a,b) -> max.(a,b), circulation_window)
            gmin = reduce((a,b) -> min.(a,b), circulation_window)
            circulation_scale = max(maximum(abs,gmax), maximum(abs,gmin),
                1e-12*reference.V*reference.c)
            circulation_residual = maximum(abs,gmax-gmin) / circulation_scale
            converged = step >= minimum_steps && force_residual <= force_tolerance &&
                circulation_residual <= circulation_tolerance
        end
        if step % progress_frequency == 0 || converged || step == steps
            println("  Wake step $step/$steps, t=$(round(step*dt; digits=5)) s: " *
                "r_force=$(round(force_residual; sigdigits=3)), " *
                "r_Gamma=$(round(circulation_residual; sigdigits=3))")
        end
        converged && break
    end
    return (; aerodynamic_problem=aero, nodal_loads, converged,
        time=aero.state.timeNow, force_residual, circulation_residual)
end

function run_pazy_wing_uvlm_static(; airspeed, density, angle_of_attack,
    chordwise_panels, spanwise_panels, symmetric_wing, maximum_wake_rows,
    core_radius, wake_shedding_fraction, aerodynamic_interaction,
    maximum_aerodynamic_time, minimum_aerodynamic_time, aerodynamic_window_time,
    aerodynamic_force_tolerance, aerodynamic_circulation_tolerance,
    maximum_static_iterations, static_relaxation, geometry_tolerance,
    load_tolerance, consecutive_equilibrium_iterations,
    newton_maximum_iterations, newton_absolute_tolerance,
    newton_relative_tolerance, newton_display_iterations, progress_frequency)

    @assert airspeed > 0 && density > 0
    @assert chordwise_panels >= 1 && spanwise_panels >= 1 && maximum_wake_rows >= 1
    @assert maximum_aerodynamic_time > 0 && minimum_aerodynamic_time >= 0
    @assert aerodynamic_window_time > 0
    @assert aerodynamic_force_tolerance > 0 && aerodynamic_circulation_tolerance > 0
    @assert maximum_static_iterations >= 1 && 0 < static_relaxation <= 1
    @assert geometry_tolerance > 0 && load_tolerance > 0
    @assert consecutive_equilibrium_iterations >= 1 && progress_frequency >= 1

    # Same skin-on Pazy structural properties as the dynamic example and
    # AeroBeams' Pazy model. Gravity is zero, as in PazyWingPitchRange.jl.
    # Incidence is applied ONCE, through UVLM's freestream angle, not root p0.
    ne, span, chord, spar = AeroBeams.geometrical_properties_Pazy()
    beam_nodes = AeroBeams.nodal_positions_Pazy()
    beam = AeroBeams.create_Beam(name="Pazy static UVLM", length=span, nElements=ne,
        normalizedNodalPositions=beam_nodes,
        S=AeroBeams.stiffness_matrices_Pazy(withSkin=true,
            sweepStructuralCorrections=false, GAy=1e16, GAz=1e16, Λ=0.0),
        I=AeroBeams.inertia_matrices_Pazy(withSkin=true),
        rotationParametrization="E231", p0=[-pi/2;0.0;0.0], aeroSurface=nothing)
    clamp = AeroBeams.create_BC(name="root clamp", beam=beam, node=1,
        types=["u1A","u2A","u3A","p1A","p2A","p3A"], values=zeros(6))
    bcs = [clamp]
    for node in 1:ne+1
        push!(bcs, AeroBeams.create_BC(name="UVLM loads $node", beam=beam, node=node,
            types=["F1A","F2A","F3A","M1A","M2A","M3A"], values=zeros(6)))
    end
    model = AeroBeams.create_Model(beams=[beam], BCs=bcs,
        gravityVector=zeros(3), v_A=t -> zeros(3))
    newton = AeroBeams.create_NewtonRaphson(maximumIterations=newton_maximum_iterations,
        absoluteTolerance=newton_absolute_tolerance, relativeTolerance=newton_relative_tolerance,
        displayStatus=newton_display_iterations, alwaysUpdateJacobian=true,
        trackingLoadSteps=false)
    # SteadyProblem: nonlinear internal-force balance, with NO structural inertia.
    structure = AeroBeams.create_SteadyProblem(model=model, systemSolver=newton)
    AeroBeams.solve!(structure)
    structure.systemSolver.convergedFinalSolution || error("Initial static solve failed")

    weights = zeros(spanwise_panels+1, ne+1)
    for j in 1:spanwise_panels+1
        eta = (j-1)/spanwise_panels
        left = min(searchsortedlast(beam_nodes,eta), ne)
        fraction = (eta-beam_nodes[left])/(beam_nodes[left+1]-beam_nodes[left])
        weights[j,left] = 1-fraction
        weights[j,left+1] = fraction
    end
    dt = chord / (chordwise_panels*airspeed)
    @assert maximum_aerodynamic_time >= dt "Aerodynamic time must include at least one step"
    reference = UVLM.Reference(span*chord,chord,span,zeros(3),airspeed,density)
    point = UVLM.create_OperatingPoint(airspeed=airspeed, density=density,
        angleOfAttack=angle_of_attack, sideslip=0.0)
    solver = UVLM.create_UVLMSolver(coreRadius=core_radius,
        wakeSheddingFraction=wake_shedding_fraction, interaction=aerodynamic_interaction)

    applied_loads = zeros(6,ne+1)
    geometry_residual = Inf
    settled_iterations = 0
    total_aerodynamic_time = 0.0
    history_rows = Vector{Float64}[]
    history_columns = ["iteration", "cumulative_aero_time_s", "wake_time_s",
        "tip_u1_A_m", "tip_u2_A_m", "tip_u3_A_m", "tip_OOP_m", "tip_IP_m",
        "tip_twist_deg", "geometry_residual", "load_residual",
        "wake_force_residual", "wake_circulation_residual", "wake_converged"]
    wake = nothing
    converged = false
    reason = "Maximum static iterations reached; equilibrium NOT established"
    println("STATIC Pazy / UVLM: U=$airspeed m/s, dt=$dt s, wake=$maximum_wake_rows rows")
    println("Each trial holds the wing fixed while UVLM settles, then solves nonlinear static balance.")

    for iteration in 1:maximum_static_iterations
        println("Static iteration $iteration/$maximum_static_iterations")
        grid, _, offsets = wing_geometry(structure.model, chord, spar, chordwise_panels, weights)
        wake = settle_static_uvlm(grid, offsets, weights, reference, point, solver;
            symmetric_wing, maximum_wake_rows, dt,
            maximum_time=maximum_aerodynamic_time, minimum_time=minimum_aerodynamic_time,
            window_time=aerodynamic_window_time, force_tolerance=aerodynamic_force_tolerance,
            circulation_tolerance=aerodynamic_circulation_tolerance, progress_frequency)
        total_aerodynamic_time += wake.time
        load_residual = interface_load_residual(wake.nodal_loads, applied_loads, chord)
        tip = structure.model.elements[end].nodalStates
        twist = wingtip_twist_degrees(structure.model)
        push!(history_rows, [iteration, total_aerodynamic_time, wake.time,
            tip.u_n2..., tip.u_n2_b[3], -tip.u_n2_b[2], twist,
            geometry_residual, load_residual, wake.force_residual,
            wake.circulation_residual, Float64(wake.converged)])
        println("  Tip OOP=$(round(tip.u_n2_b[3]; digits=7)) m, " *
            "twist=$(round(twist; digits=5)) deg; r_geometry=$geometry_residual, r_load=$load_residual")

        if !wake.converged
            reason = "Wake did not settle within maximum_aerodynamic_time; equilibrium NOT established"
            break
        end
        # Check the UNRELAXED aero/structural load mismatch, so a small relaxation
        # factor cannot masquerade as equilibrium. The final wake and structure
        # have exactly the same geometry: on convergence we do not move it again.
        settled_iterations = geometry_residual <= geometry_tolerance &&
            load_residual <= load_tolerance ? settled_iterations+1 : 0
        if settled_iterations >= consecutive_equilibrium_iterations
            converged = true
            reason = "Static aeroelastic equilibrium reached"
            break
        end
        iteration == maximum_static_iterations && break
        applied_loads .= (1-static_relaxation).*applied_loads .+
            static_relaxation.*wake.nodal_loads
        apply_pazy_nodal_loads!(structure, applied_loads, 0.0)
        AeroBeams.solve!(structure)
        structure.systemSolver.convergedFinalSolution ||
            error("Static Newton solve failed at iteration $iteration; no equilibrium claimed")
        next_grid, _, _ = wing_geometry(structure.model, chord, spar, chordwise_panels, weights)
        geometry_residual = maximum(abs,next_grid-grid)/chord
    end

    # Preserve a standard AeroBeams steady solution for later inspection/plotting.
    AeroBeams.save_load_factor_data!(structure, structure.σ, structure.x)
    tip = structure.model.elements[end].nodalStates
    println(reason)
    return (; converged, termination_reason=reason, structural_problem=structure,
        aerodynamic_problem=wake.aerodynamic_problem, dt, airspeed, angle_of_attack,
        history=permutedims(hcat(history_rows...)), history_columns,
        tip_displacement_A=copy(tip.u_n2), tip_out_of_plane=tip.u_n2_b[3],
        tip_in_plane=-tip.u_n2_b[2], tip_twist_degrees=wingtip_twist_degrees(structure.model),
        applied_nodal_loads=copy(applied_loads), aerodynamic_nodal_loads=wake.nodal_loads,
        geometry_residual, load_residual=interface_load_residual(wake.nodal_loads,applied_loads,chord),
        total_aerodynamic_time)
end
