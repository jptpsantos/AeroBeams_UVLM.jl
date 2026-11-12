# Nonlinear STATIC beam balance with time-marched UVLM aerodynamics.
# Reuse the verified coordinate, geometry, and load-transfer functions.
if !isdefined(@__MODULE__,:run_pazy_wing_uvlm)
    include(joinpath(@__DIR__,"PazyWingUVLMCoupling.jl"))
end

# Hold the wing fixed and let its free wake develop. Each new static trial
# starts with an empty wake: a shape iteration is NOT a physical wing velocity.
function settle_static_uvlm(grid,offsets,weights,reference,operating_point,solver;
    symmetric_wing,maximum_wake_rows,dt,maximum_time,minimum_time,
    window_time,force_tolerance,circulation_tolerance,progress_frequency)

    steps = max(1,floor(Int,maximum_time / dt))
    windowSteps = max(2,ceil(Int,window_time / dt) + 1)
    # Do not accept the initial wake transient as equilibrium.
    minimumSteps = max(ceil(Int,minimum_time / dt),2maximum_wake_rows)
    model = UVLM.create_UVLMModel(surfaces=[grid],reference=reference,
        symmetric=symmetric_wing,interactionGroups=[1])
    aero = UVLM.create_UVLMDynamicProblem(model=model,aeroSolver=solver,
        operatingPoint=operating_point,timeVector=collect(0:steps) .* dt,
        maximumWakeRows=maximum_wake_rows,trackingTimeSteps=false)
    forceWindow = Matrix{Float64}[]
    circulationWindow = Vector{Float64}[]
    nodal_loads = zeros(6,size(weights,2))
    spanMap = PazySpanMap(weights)
    force_residual = circulation_residual = Inf
    converged = false
    for step = 1:steps
        UVLM.begin_time_step!(aero)
        try
            loads = UVLM.evaluate_trial!(aero; surfaces=[grid])
            beam_loads!(nodal_loads,loads,offsets,spanMap)
            UVLM.commit_time_step!(aero)
        catch
            UVLM.rollback_time_step!(aero)
            rethrow()
        end
        push!(forceWindow,copy(nodal_loads))
        push!(circulationWindow,copy(aero.state.system.Γ))
        if length(forceWindow) > windowSteps
            popfirst!(forceWindow)
            popfirst!(circulationWindow)
        end
        if length(forceWindow) == windowSteps
            # Peak-to-peak change across the WHOLE window, not just two steps.
            fmax = reduce((a,b) -> max.(a,b),forceWindow)
            fmin = reduce((a,b) -> min.(a,b),forceWindow)
            force_residual = interface_load_residual(fmax,fmin,reference.c)
            gmax = reduce((a,b) -> max.(a,b),circulationWindow)
            gmin = reduce((a,b) -> min.(a,b),circulationWindow)
            circulationScale = max(maximum(abs,gmax),maximum(abs,gmin),
                1e-12 * reference.V * reference.c)
            circulation_residual = maximum(abs,gmax - gmin) / circulationScale
            converged =
                step >= minimumSteps && force_residual <= force_tolerance &&
                circulation_residual <= circulation_tolerance
        end
        if step % progress_frequency == 0 || converged || step == steps
            println(
                "  Wake step $step/$steps, t=$(round(step*dt; digits=5)) s: " *
                "r_force=$(round(force_residual; sigdigits=3)), " *
                "r_Gamma=$(round(circulation_residual; sigdigits=3))"
            )
        end
        converged && break
    end
    return (; aerodynamic_problem=aero,nodal_loads,converged,
        time=aero.state.timeNow,force_residual,circulation_residual)
end

function run_pazy_wing_uvlm_static(; airspeed,density,angle_of_attack,
    chordwise_panels,spanwise_panels,symmetric_wing,maximum_wake_rows,
    core_radius,wake_shedding_fraction,aerodynamic_interaction,
    maximum_aerodynamic_time,minimum_aerodynamic_time,aerodynamic_window_time,
    aerodynamic_force_tolerance,aerodynamic_circulation_tolerance,
    maximum_static_iterations,static_relaxation,geometry_tolerance,
    load_tolerance,consecutive_equilibrium_iterations,
    newton_maximum_iterations,newton_absolute_tolerance,
    newton_relative_tolerance,newton_display_iterations,progress_frequency,
    time_step_chords=1 / chordwise_panels)

    @assert airspeed > 0 && density > 0
    @assert chordwise_panels >= 1 && spanwise_panels >= 1 && maximum_wake_rows >= 1
    @assert maximum_aerodynamic_time > 0 && minimum_aerodynamic_time >= 0
    @assert aerodynamic_window_time > 0
    @assert aerodynamic_force_tolerance > 0 && aerodynamic_circulation_tolerance > 0
    @assert maximum_static_iterations >= 1 && 0 < static_relaxation <= 1
    @assert geometry_tolerance > 0 && load_tolerance > 0
    @assert consecutive_equilibrium_iterations >= 1 && progress_frequency >= 1
    @assert isfinite(time_step_chords) && time_step_chords > 0

    # Same skin-on Pazy structural properties as the dynamic example and
    # AeroBeams' Pazy model. Gravity is zero, as in PazyWingPitchRange.jl.
    # Incidence is applied ONCE, through UVLM's freestream angle, not root p0.
    ne,span,chord,spar = AeroBeams.geometrical_properties_Pazy()
    beamNodes = last(AeroBeams.geometrical_properties_Pazy())
    beam = AeroBeams.create_Beam(name="Pazy static UVLM",
        length=span,
        nElements=ne,
        normalizedNodalPositions=beamNodes,
        S=AeroBeams.stiffness_matrices_Pazy(withSkin=true,
            sweepStructuralCorrections=false,GAy=1e16,GAz=1e16,Λ=0.0),
        I=AeroBeams.inertia_matrices_Pazy(withSkin=true),
        p0=AeroBeams.transform_to_WM([-pi / 2; 0.0; 0.0]; from=:E231),
        surface=nothing)
    clamp = AeroBeams.create_BC(name="root clamp",
        beam=beam,
        node=1,
        conditions=Dict(Symbol.(["u1A","u2A","u3A","p1A","p2A","p3A"]) .=> zeros(6)))
    bcs = [clamp]
    for node = 1:ne+1
        push!(
            bcs,
            AeroBeams.create_BC(name="UVLM loads $node",
        beam=beam,
        node=node,
        conditions=Dict(Symbol.(["F1A","F2A","F3A","M1A","M2A","M3A"]) .=> zeros(6)))
        )
    end
    model = AeroBeams.create_Model(beams=[beam],BCs=bcs,
        gravityVector=zeros(3),v_A=t -> zeros(3))
    newton = AeroBeams.create_NewtonRaphson(maxIter=newton_maximum_iterations,
        absTol=newton_absolute_tolerance,
        relTol=newton_relative_tolerance,
        displayStatus=newton_display_iterations,
        method=:Newton,
        trackLoadSteps=false)
    # SteadyProblem: nonlinear internal-force balance, with NO structural inertia.
    structure = AeroBeams.create_SteadyProblem(model=model,systemSolver=newton)
    AeroBeams.solve!(structure)
    structure.systemSolver.loadStepState.convergedFinalSolution || error("Initial static solve failed")

    weights = zeros(spanwise_panels + 1,ne + 1)
    for j = 1:spanwise_panels+1
        eta = (j - 1) / spanwise_panels
        left = min(searchsortedlast(beamNodes,eta),ne)
        fraction = (eta - beamNodes[left]) / (beamNodes[left+1] - beamNodes[left])
        weights[j,left] = 1 - fraction
        weights[j,left+1] = fraction
    end
    dt = time_step_chords * chord / airspeed
    @assert maximum_aerodynamic_time >= dt "Aerodynamic time must include at least one step"
    reference = UVLM.Reference(span * chord,chord,span,zeros(3),airspeed,density)
    point = UVLM.create_OperatingPoint(airspeed=airspeed,density=density,
        angleOfAttack=angle_of_attack,sideslip=0.0)
    solver = UVLM.create_UVLMSolver(coreRadius=core_radius,
        wakeSheddingFraction=wake_shedding_fraction,interaction=aerodynamic_interaction)

    appliedLoads = zeros(6,ne + 1)
    geometry_residual = Inf
    settledIterations = 0
    total_aerodynamic_time = 0.0
    historyRows = Vector{Float64}[]
    history_columns = ["iteration","cumulative_aero_time_s","wake_time_s",
        "tip_u1_A_m","tip_u2_A_m","tip_u3_A_m","tip_OOP_m","tip_IP_m",
        "tip_twist_deg","geometry_residual","load_residual",
        "wake_force_residual","wake_circulation_residual","wake_converged"]
    wake = nothing
    converged = false
    reason = "Maximum static iterations reached; equilibrium NOT established"
    println("STATIC Pazy / UVLM: U=$airspeed m/s, dt=$dt s, wake=$maximum_wake_rows rows")
    println("Each trial holds the wing fixed while UVLM settles, then solves nonlinear static balance.")

    for iteration = 1:maximum_static_iterations
        println("Static iteration $iteration/$maximum_static_iterations")
        grid,_,offsets = wing_geometry(structure.model,chord,spar,chordwise_panels,weights)
        wake = settle_static_uvlm(grid,offsets,weights,reference,point,solver;
            symmetric_wing,maximum_wake_rows,dt,
            maximum_time=maximum_aerodynamic_time,minimum_time=minimum_aerodynamic_time,
            window_time=aerodynamic_window_time,force_tolerance=aerodynamic_force_tolerance,
            circulation_tolerance=aerodynamic_circulation_tolerance,progress_frequency)
        total_aerodynamic_time += wake.time
        load_residual = interface_load_residual(wake.nodal_loads,appliedLoads,chord)
        tip = structure.model.assemblyData.elements[end].state.nodal
        twist = wingtip_twist_degrees(structure.model)
        push!(
            historyRows,
            [iteration,total_aerodynamic_time,wake.time,
                tip.u2...,tip.u2_b[3],-tip.u2_b[2],twist,
                geometry_residual,load_residual,wake.force_residual,
                wake.circulation_residual,Float64(wake.converged)]
        )
        println(
            "  Tip OOP=$(round(tip.u2_b[3]; digits=7)) m, " *
            "twist=$(round(twist; digits=5)) deg; r_geometry=$geometry_residual, r_load=$load_residual"
        )

        if !wake.converged
            reason = "Wake did not settle within maximum_aerodynamic_time; equilibrium NOT established"
            break
        end
        # Check the UNRELAXED aero/structural load mismatch, so a small relaxation
        # factor cannot masquerade as equilibrium. The final wake and structure
        # have exactly the same geometry: on convergence we do not move it again.
        settledIterations =
            geometry_residual <= geometry_tolerance &&
            load_residual <= load_tolerance ? settledIterations + 1 : 0
        if settledIterations >= consecutive_equilibrium_iterations
            converged = true
            reason = "Static aeroelastic equilibrium reached"
            break
        end
        iteration == maximum_static_iterations && break
        appliedLoads .= (1 - static_relaxation) .* appliedLoads .+
                        static_relaxation .* wake.nodal_loads
        apply_pazy_nodal_loads!(structure,appliedLoads,0.0)
        AeroBeams.solve!(structure)
        structure.systemSolver.loadStepState.convergedFinalSolution ||
            error("Static Newton solve failed at iteration $iteration; no equilibrium claimed")
        nextGrid,_,_ = wing_geometry(structure.model,chord,spar,chordwise_panels,weights)
        geometry_residual = maximum(abs,nextGrid - grid) / chord
    end

    # Preserve a standard AeroBeams steady solution for later inspection/plotting.
    AeroBeams.save_loadstep_data!(structure)
    tip = structure.model.assemblyData.elements[end].state.nodal
    println(reason)
    return (; converged,termination_reason=reason,structural_problem=structure,
        aerodynamic_problem=wake.aerodynamic_problem,dt,airspeed,angle_of_attack,
        history=permutedims(hcat(historyRows...)),history_columns,
        tip_displacement_A=copy(tip.u2),tip_out_of_plane=tip.u2_b[3],
        tip_in_plane=-tip.u2_b[2],tip_twist_degrees=wingtip_twist_degrees(structure.model),
        applied_nodal_loads=copy(appliedLoads),aerodynamic_nodal_loads=wake.nodal_loads,
        geometry_residual,load_residual=interface_load_residual(wake.nodal_loads,appliedLoads,chord),
        total_aerodynamic_time)
end
