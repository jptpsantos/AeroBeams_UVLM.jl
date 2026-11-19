# Nonlinear STATIC beam balance with time-marched UVLM aerodynamics.
# Reuse the verified coordinate, geometry, and load-transfer functions.
if !isdefined(@__MODULE__,:run_pazy_wing_uvlm)
    include(joinpath(@__DIR__,"PazyWingUVLMCoupling.jl"))
end

# Elementwise maximum and minimum over a window of arrays, written in place
function window_extrema!(xmax,xmin,window)
    copyto!(xmax,window[1])
    copyto!(xmin,window[1])
    for k = 2:length(window)
        x = window[k]
        for i in eachindex(xmax,xmin,x)
            xmax[i] = max(xmax[i],x[i])
            xmin[i] = min(xmin[i],x[i])
        end
    end
    return nothing
end

# Hold the wing fixed and let its free wake develop. By default each static trial
# starts with an empty wake: a shape iteration is NOT a physical wing velocity.
# With `aero` given (wake reuse), the accepted UVLM state is moved onto the new
# shape WITHOUT imposing a surface velocity and its wake and circulation are the
# initial condition: the settled wake of a nearby shape needs far fewer steps than
# a fresh one. The acceptance test (window of force and circulation changes) is
# the same; `reuse_minimum_steps` replaces the two-wake-ages minimum (default: the
# window length). The operating point may differ from the one of `aero` (e.g. a
# perturbed incidence) but not the time step.
function settle_static_uvlm(grid,offsets,weights,reference,operating_point,solver;
    symmetric_wing,maximum_wake_rows,dt,maximum_time,minimum_time,
    window_time,force_tolerance,circulation_tolerance,progress_frequency,
    aero=nothing,reuse_minimum_steps=nothing)

    steps = max(1,floor(Int,maximum_time / dt))
    windowSteps = max(2,ceil(Int,window_time / dt) + 1)
    # Do not accept the initial wake transient as equilibrium.
    minimumSteps = max(ceil(Int,minimum_time / dt),2maximum_wake_rows)
    if isnothing(aero)
        model = UVLM.create_UVLMModel(surfaces=[grid],reference=reference,
            symmetric=symmetric_wing,interactionGroups=[1])
        aero = UVLM.create_UVLMDynamicProblem(model=model,aeroSolver=solver,
            operatingPoint=operating_point,timeVector=collect(0:steps) .* dt,
            maximumWakeRows=maximum_wake_rows,trackingTimeSteps=false)
    else
        aero.operatingPoint = operating_point
        for system in (aero.state.system,aero.workspace.system,aero.workspace.commitSystem)
            system.freestream[] = UVLM.freestream(operating_point)
        end
        UVLM.apply_trial_geometry!(aero.state.system,[grid],aero.aeroSolver)
        UVLM.copy_surfaces_to_previous!(aero.state.system,length(aero.state.system.surfaces))
        timeVector = aero.timeVector
        while length(timeVector) < aero.state.stepIndex + steps + 1
            push!(timeVector,timeVector[end] + dt)
        end
        steps = min(steps,length(timeVector) - aero.state.stepIndex)
        minimumSteps = something(reuse_minimum_steps,windowSteps)
    end
    nodal_loads = zeros(6,size(weights,2))
    # Circular buffers holding the last windowSteps loads and circulations, and
    # buffers for their extrema, all reused at every step
    forceWindow = [similar(nodal_loads) for _ = 1:windowSteps]
    fmax = similar(nodal_loads)
    fmin = similar(nodal_loads)
    circulationWindow = Vector{Float64}[]
    gmax = gmin = Float64[]
    spanMap = PazySpanMap(weights)
    force_residual = circulation_residual = Inf
    converged = false
    performedSteps = 0
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
        performedSteps = step
        if step == 1
            # Allocate the circulation buffers once its size is known
            Γ = aero.state.system.Γ
            circulationWindow = [similar(Γ) for _ = 1:windowSteps]
            gmax = similar(Γ)
            gmin = similar(Γ)
        end
        slot = mod1(step,windowSteps)
        copyto!(forceWindow[slot],nodal_loads)
        copyto!(circulationWindow[slot],aero.state.system.Γ)
        if step >= windowSteps
            # Peak-to-peak change across the WHOLE window, not just two steps.
            window_extrema!(fmax,fmin,forceWindow)
            force_residual = interface_load_residual(fmax,fmin,reference.c)
            window_extrema!(gmax,gmin,circulationWindow)
            circulationScale = max(maximum(abs,gmax),maximum(abs,gmin),
                1e-12 * reference.V * reference.c)
            circulation_residual = mapreduce((a,b) -> abs(a - b),max,gmax,gmin) / circulationScale
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
        time=performedSteps * dt,steps=performedSteps,force_residual,circulation_residual)
end

# Apply the loads x + step to the structure and solve. If the Newton solve fails,
# or the wing geometry (a function returning the UVLM grid) moves by more than
# `geometry_cap` (same length unit as the grid) from the one before the call, the
# structural state is restored and the step halved, up to `maximum_halvings` times.
# Returns the applied loads, the step fraction used and the number of halvings.
function solve_structure_safely!(structure,x,step,geometry,geometry_cap; maximum_halvings=6)
    stateGood = copy(structure.systemArrays.x)
    geometryBefore = isfinite(geometry_cap) ? geometry() : nothing
    fraction = 1.0
    for attempt = 0:maximum_halvings
        loads = x .+ fraction .* step
        apply_pazy_nodal_loads!(structure,loads,0.0)
        converged = try
            AeroBeams.solve!(structure)
            structure.systemSolver.loadStepState.convergedFinalSolution
        catch
            false
        end
        if converged && attempt < maximum_halvings && !isnothing(geometryBefore) &&
            maximum(abs,geometry() - geometryBefore) > geometry_cap
            converged = false
        end
        converged && return loads,fraction,attempt
        structure.systemArrays.x .= stateGood
        fraction /= 2
    end
    error("Static Newton solve failed even after halving the load step $maximum_halvings times; no equilibrium claimed")
end

# Anderson-accelerated update of the applied loads x for the fixed-point map
# G(x) = UVLM loads on the shape that the structure takes under x. With
# r = G(x) - x, the plain (relaxed) step is β r. The depth-m history of
# increments (dX, dR) is used in a least-squares problem on the residual, with
# moments scaled by the chord, so that the combined step also removes the slowest
# modes (the map has modes of both signs, which a single relaxation factor cannot
# damp simultaneously). Returns the increment to add to x.
function anderson_step!(dX,dR,x,xPrevious,r,rPrevious,β,depth,loadWeights)
    step = β .* r
    if !isnothing(xPrevious)
        push!(dX,x .- xPrevious)
        push!(dR,r .- rPrevious)
        length(dX) > depth && (popfirst!(dX); popfirst!(dR))
        A = hcat([vec(loadWeights .* d) for d in dR]...)
        # Rank-revealing solve: the history can be degenerate (e.g. zero loads at zero incidence
        # give zero residual increments), which makes the normal equations singular
        γ = pinv(A; rtol=1e-8) * vec(loadWeights .* r)
        for (k,g) in enumerate(γ)
            step .-= g .* (dX[k] .+ β .* dR[k])
        end
    end
    return step
end

# Options of the acceleration (defaults reproduce the original fixed-relaxation,
# fresh-wake, cold-start algorithm):
#  - relaxation_mode=:fixed: loads <- (1-ω) loads + ω G(loads), ω = static_relaxation.
#    :anderson: Anderson mixing of depth anderson_depth, static_relaxation is its
#    mixing parameter β (1.0 recommended; no tuning to the case is needed).
#  - geometry_cap_chords: maximum change of the UVLM grid per load update, in chords.
#    The load step is halved (and the Anderson history dropped) until it complies or
#    the structural solve converges. Avoids huge jumps at high loads, from which the
#    wake of a very deformed wing may never settle (recommended: 2).
#  - wake_reuse: keep wake and circulation between static iterations (see above).
#  - warm_start: `warm_start` field of a previous result (same structural options
#    and spanwise/chordwise discretization); the structure is reused and moved to
#    `warm_start_target_loads` (default its own loads, e.g. scaled by the dynamic
#    pressure ratio). Its `aero` (wake) is used only if wake_reuse and the time step
#    is the same; set it to `nothing` otherwise.
function run_pazy_wing_uvlm_static(; airspeed,density,angle_of_attack,
    chordwise_panels,spanwise_panels,symmetric_wing,maximum_wake_rows,
    core_radius,wake_shedding_fraction,aerodynamic_interaction,
    maximum_aerodynamic_time,minimum_aerodynamic_time,aerodynamic_window_time,
    aerodynamic_force_tolerance,aerodynamic_circulation_tolerance,
    maximum_static_iterations,static_relaxation,geometry_tolerance,
    load_tolerance,consecutive_equilibrium_iterations,
    newton_maximum_iterations,newton_absolute_tolerance,
    newton_relative_tolerance,newton_display_iterations,progress_frequency,
    mesh=:nominal,sweepAngle=0.0,useAePW3Matrices=false,
    withWingtipRod=true,rootTorsionCorr=true,beamWingtipRod=false,
    nElementsWingtipRod=12,sweptTipMassConfig=:baseline,
    sweptTipMassInertiaCorrection=false,
    time_step_chords=1 / chordwise_panels,load_residual_force_floor=1.0,
    relaxation_mode=:anderson,anderson_depth=4,geometry_cap_chords=2,
    wake_reuse=true,warm_start=nothing,warm_start_target_loads=nothing)

    @assert airspeed > 0 && density > 0
    @assert chordwise_panels >= 1 && spanwise_panels >= 1 && maximum_wake_rows >= 1
    @assert maximum_aerodynamic_time > 0 && minimum_aerodynamic_time >= 0
    @assert aerodynamic_window_time > 0
    @assert aerodynamic_force_tolerance > 0 && aerodynamic_circulation_tolerance > 0
    @assert maximum_static_iterations >= 1 && 0 < static_relaxation <= 1
    @assert geometry_tolerance > 0 && load_tolerance > 0
    @assert isfinite(load_residual_force_floor) && load_residual_force_floor > 0
    @assert consecutive_equilibrium_iterations >= 1 && progress_frequency >= 1
    @assert isfinite(time_step_chords) && time_step_chords > 0
    @assert relaxation_mode in (:fixed,:anderson) "relaxation_mode must be :fixed or :anderson"
    @assert anderson_depth >= 1 && geometry_cap_chords > 0

    # Share the native Pazy structure with the dynamic solver. Incidence is
    # applied once through UVLM's freestream; gravity and root pitch are zero.
    pazy = isnothing(warm_start) ? create_pazy_uvlm_structure(; mesh,sweepAngle,useAePW3Matrices, withWingtipRod,rootTorsionCorr,beamWingtipRod,nElementsWingtipRod, sweptTipMassConfig,sweptTipMassInertiaCorrection) : warm_start.pazy
    (; model,span,chord,beamNodes) = pazy
    ne,spar = pazy.n_elements,pazy.sparFraction
    projectedSpan = span * cos(sweepAngle)
    structure = if isnothing(warm_start)
        # newton = AeroBeams.create_NewtonRaphson(maxIter=newton_maximum_iterations,
        #     absTol=newton_absolute_tolerance,
        #     relTol=newton_relative_tolerance,
        #     displayStatus=newton_display_iterations,
        #     method=:Newton,
        #     trackLoadSteps=false)
        newton = AeroBeams.create_NLsolver(method=:newton,
            maxIter=newton_maximum_iterations,
            ftol=newton_absolute_tolerance,
            displayStatus=newton_display_iterations,trackLoadSteps=false)
        # SteadyProblem: nonlinear internal-force balance, with NO structural inertia.
        problem = AeroBeams.create_SteadyProblem(model=model,systemSolver=newton)
        AeroBeams.solve!(problem)
        problem.systemSolver.loadStepState.convergedFinalSolution || error("Initial static solve failed")
        problem
    else
        warm_start.structure
    end

    weights = zeros(spanwise_panels + 1,ne + 1)
    for j = 1:spanwise_panels+1
        eta = (j - 1) / spanwise_panels
        left = min(searchsortedlast(beamNodes,eta),ne)
        fraction = (eta - beamNodes[left]) / (beamNodes[left+1] - beamNodes[left])
        weights[j,left] = 1 - fraction
        weights[j,left+1] = fraction
    end
    currentGrid() = wing_geometry(structure.model,chord,spar,chordwise_panels,weights; sweepAngle)[1]
    dt = time_step_chords * chord / airspeed
    @assert maximum_aerodynamic_time >= dt "Aerodynamic time must include at least one step"
    reference = UVLM.Reference(projectedSpan * chord,chord,projectedSpan,zeros(3),airspeed,density)
    point = UVLM.create_OperatingPoint(airspeed=airspeed,density=density,
        angleOfAttack=angle_of_attack,sideslip=0.0)
    solver = UVLM.create_UVLMSolver(coreRadius=core_radius,
        wakeSheddingFraction=wake_shedding_fraction,interaction=aerodynamic_interaction)

    appliedLoads = zeros(6,ne + 1)
    if !isnothing(warm_start)
        # The structure sits at warm_start.loads: move it to the target loads.
        target = isnothing(warm_start_target_loads) ? warm_start.loads : warm_start_target_loads
        appliedLoads .= warm_start.loads
        loads,_,_ = solve_structure_safely!(structure,warm_start.loads,target .- warm_start.loads,
            currentGrid,geometry_cap_chords * chord)
        appliedLoads .= loads
    end
    # Loads scaled for the least-squares problem of the Anderson mixing: forces as
    # they are, moments divided by the chord
    loadWeights = ones(6,ne + 1)
    loadWeights[4:6,:] .= 1 / chord
    andersonX = Matrix{Float64}[]
    andersonR = Matrix{Float64}[]
    xPrevious = rPrevious = nothing
    geometry_residual = Inf
    settledIterations = 0
    total_aerodynamic_time = 0.0
    total_aerodynamic_steps = 0
    historyRows = Vector{Float64}[]
    history_columns = ["iteration","cumulative_aero_time_s","wake_time_s",
        "tip_u1_A_m","tip_u2_A_m","tip_u3_A_m","tip_OOP_m","tip_IP_m",
        "tip_twist_deg","geometry_residual","load_residual",
        "wake_force_residual","wake_circulation_residual","wake_converged"]
    wake = nothing
    aero = wake_reuse && !isnothing(warm_start) ? warm_start.aero : nothing
    converged = false
    reason = "Maximum static iterations reached; equilibrium NOT established"
    println("STATIC Pazy / UVLM: U=$airspeed m/s, dt=$dt s, wake=$maximum_wake_rows rows")
    println("Each trial holds the wing fixed while UVLM settles, then solves nonlinear static balance.")

    for iteration = 1:maximum_static_iterations
        println("Static iteration $iteration/$maximum_static_iterations")
        grid,_,offsets = wing_geometry(structure.model,chord,spar,chordwise_panels,weights; sweepAngle)
        wake = settle_static_uvlm(grid,offsets,weights,reference,point,solver;
            symmetric_wing,maximum_wake_rows,dt,
            maximum_time=maximum_aerodynamic_time,minimum_time=minimum_aerodynamic_time,
            window_time=aerodynamic_window_time,force_tolerance=aerodynamic_force_tolerance,
            circulation_tolerance=aerodynamic_circulation_tolerance,progress_frequency,
            aero=wake_reuse ? aero : nothing)
        aero = wake.aerodynamic_problem
        total_aerodynamic_time += wake.time
        total_aerodynamic_steps += wake.steps
        load_residual = interface_load_residual(wake.nodal_loads,appliedLoads,chord;
            force_floor=load_residual_force_floor)
        tip = structure.model.beams[1].elements[end].state.nodal
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
        x = copy(appliedLoads)
        r = wake.nodal_loads .- x
        step = if relaxation_mode == :fixed
            static_relaxation .* r
        else
            anderson_step!(andersonX,andersonR,x,xPrevious,r,rPrevious,static_relaxation,
                anderson_depth,loadWeights)
        end
        xPrevious,rPrevious = x,r
        appliedLoads,_,halvings = solve_structure_safely!(structure,x,step,currentGrid,
            geometry_cap_chords * chord)
        if halvings > 0
            println("  Load step halved $halvings time(s) (failed structural solve or geometry change above the cap)")
            empty!(andersonX); empty!(andersonR)
            xPrevious = rPrevious = nothing
        end
        nextGrid = currentGrid()
        geometry_residual = maximum(abs,nextGrid - grid) / chord
    end

    # Preserve a standard AeroBeams steady solution for later inspection/plotting.
    AeroBeams.save_loadstep_data!(structure)
    tip = structure.model.beams[1].elements[end].state.nodal
    println(reason)
    return (; converged,termination_reason=reason,structural_problem=structure,
        pazy_parameters=pazy.parameters,
        aerodynamic_problem=wake.aerodynamic_problem,dt,airspeed,angle_of_attack,
        history=permutedims(hcat(historyRows...)),history_columns,
        tip_displacement_A=copy(tip.u2),tip_out_of_plane=tip.u2_b[3],
        tip_in_plane=-tip.u2_b[2],tip_twist_degrees=wingtip_twist_degrees(structure.model),
        applied_nodal_loads=copy(appliedLoads),aerodynamic_nodal_loads=wake.nodal_loads,
        geometry_residual,load_residual=interface_load_residual(wake.nodal_loads,appliedLoads,chord;
            force_floor=load_residual_force_floor),
        total_aerodynamic_time,total_aerodynamic_steps,
        warm_start=(; pazy,structure,loads=copy(appliedLoads),aero=wake.aerodynamic_problem))
end
