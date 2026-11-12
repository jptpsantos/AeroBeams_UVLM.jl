"""Editable defaults for the additive nonlinear Chang verification example.

Unlike the legacy perturbation solver, mean aerodynamic loads are retained.
The default outboard installation and enabled interference exercise structural
and aerodynamic wing/rotor coupling. Original Chang defaults are not edited.
"""
function chang_aerobeams_defaults()
    base = Chang.chang_case_defaults()
    reference = merge(
        base,
        (
            wing=merge(base.wing,(symmetric=false,)),
            propeller=merge(base.propeller,(attachment_eta=[0.65],)),
            simulation=merge(
                base.simulation,
                (
                    freestream_speed_mps=65.0,end_time_s=0.25,interaction_on=true,
                    near_field_force_model=:imperial,propeller_moment_projection=:exact_virtual_work)
            ),
            output=merge(base.output,(plot_results=false,animate_wake=false))
        )
    )
    solver = (support=:reduced,time_step_s=nothing,speed_policy=:fixed_rpm,
        maximum_newton_iterations=25,newton_tolerance=1e-8,
        maximum_coupling_iterations=30,coupling_relaxation=0.6,
        geometry_tolerance=1e-7,velocity_tolerance=1e-6,load_tolerance=1e-5,
        coupled_equilibrium_tolerance=1e-6,
        maximum_relative_angle=0.35,progress_frequency=100,
        impulse_start_s=0.15,impulse_duration_s=0.02,impulse_magnitude_nm=2.0,
        impulse_rotor=1,initial_eta=nothing,initial_eta_velocity=nothing)
    initialization = (mode=:mean_equilibrium,maximum_iterations=15,
        relaxation=0.6,tolerance=1e-5,minimum_revolutions=4,
        maximum_revolutions=20,periodic_load_tolerance=1e-3)
    output = (directory=joinpath(@__DIR__,"..","output"),label="chang_aerobeams",
        write_results=true,plot_results=true,save_structural_history=false)
    pylon = (elements=4,EI_pitch=nothing,EI_yaw=nothing,EA=1e10,GA=1e10,GJ=1e8,
        rotary_inertia_per_length=1e-3,rotor_transverse_inertia=nothing,
        rotor_inertia_tensor=nothing,rotor_mass=nothing,spin_inertia=nothing)
    return (; reference,solver,initialization,output,pylon)
end

mutable struct ChangAeroBeamsProblem
    reference::Any
    wing::Any
    rotors::Vector
    aerodynamic::Any
    time::Vector{Float64}
    eta::Vector{Float64}
    eta_velocity::Vector{Float64}
    eta_acceleration::Vector{Float64}
    node_history::Any
    config::Any
    phase_time_offset::Float64
    initialized::Bool
    completed::Bool
    results::Dict{Symbol,Any}
end

function _node_history(problem)
    q = wing_node_parameters(problem)
    z = zeros(size(q.u))
    return (; u=copy(q.u),p=copy(q.p),udot=copy(z),pdot=copy(z),
        uddot=copy(z),pddot=copy(z))
end

"""Build new structural/aerodynamic storage without running a response."""
function build_chang_aerobeams(config=chang_aerobeams_defaults())
    config = deepcopy(config)
    defaults = config.reference
    settings = config.solver
    0 < settings.coupling_relaxation <= 1 || throw(ArgumentError("Coupling relaxation must lie in (0,1]"))
    settings.maximum_newton_iterations >= 2 && settings.maximum_coupling_iterations >= 1 ||
        throw(ArgumentError("Invalid iteration limits"))
    all(
        x -> isfinite(x) && x > 0,
        (settings.newton_tolerance,settings.geometry_tolerance,
            settings.velocity_tolerance,settings.load_tolerance,settings.coupled_equilibrium_tolerance,
            settings.impulse_duration_s,settings.maximum_relative_angle)
    ) ||
        throw(ArgumentError("Tolerances, pulse duration and reduced angle limit must be positive"))
    config.solver.speed_policy in (:fixed_rpm,:fixed_advance_ratio) ||
        throw(ArgumentError("speed_policy must be :fixed_rpm or :fixed_advance_ratio"))
    if config.solver.speed_policy == :fixed_rpm
        defaults = merge(
            defaults,
            (propeller=merge(defaults.propeller,
                (trim_speed_mps=defaults.simulation.freestream_speed_mps,)),)
        )
    end
    resolved = Chang.load_chang_configuration(defaults; env=Dict{String,String}())
    reference = Chang.build_chang_model(resolved)
    p = reference.parameters
    p.Npropellers > 0 || throw(ArgumentError("This example requires at least one propeller"))
    1 <= settings.impulse_rotor <= p.Npropellers || throw(ArgumentError("Invalid impulse rotor"))
    reference.config.simulation.near_field_force_model == :imperial ||
        throw(ArgumentError("The nonlinear interface uses the corrected :imperial forces"))
    dt = isnothing(config.solver.time_step_s) ? first(p.dt) : config.solver.time_step_s
    isfinite(dt) && dt > 0 || throw(ArgumentError("Positive finite timestep required"))
    count = floor(Int,resolved.simulation.end_time_s / dt + 1e-10)
    count >= 1 || throw(ArgumentError("Duration must cover at least one timestep"))
    time = collect(0:count) .* dt
    resolved.structural.stiffness_damping_ratio == 0 || throw(
        ArgumentError(
            "Global Chang Rayleigh damping is not a nonlinear beam material law; set stiffness_damping_ratio=0")
    )
    config.solver.support in (:reduced,:beam_pylon) || throw(ArgumentError("Unknown support formulation"))
    beam_pylon = config.solver.support == :beam_pylon
    if beam_pylon
        resolved.structural.propeller_damping_ratio == 0 || throw(
            ArgumentError(
                "Beam pylon needs its own physical damping model; reduced support damping must be zero")
        )
        q = config.pylon
        wing = build_wing_with_beam_pylons(reference; time_vector=time,
            pylon_elements=q.elements,EI_pitch=q.EI_pitch,EI_yaw=q.EI_yaw,
            EA=q.EA,GA=q.GA,GJ=q.GJ,rotary_inertia_per_length=q.rotary_inertia_per_length)
    else
        wing = build_wing(reference; time_vector=time)
    end
    nrelative = beam_pylon ? 0 : 2p.Npropellers
    eta = isnothing(config.solver.initial_eta) ? zeros(nrelative) :
          collect(Float64,config.solver.initial_eta)
    velocity =
        isnothing(config.solver.initial_eta_velocity) ? zeros(length(eta)) :
        collect(Float64,config.solver.initial_eta_velocity)
    length(eta) == length(velocity) == nrelative ||
        throw(DimensionMismatch("Incorrect number of relative support coordinates"))
    all(isfinite,eta) && all(isfinite,velocity) || throw(ArgumentError("Nonfinite initial support state"))
    (isempty(eta) || maximum(abs,eta) <= settings.maximum_relative_angle) ||
        throw(ArgumentError("Initial relative angles exceed the reduced model's configured limit"))
    if beam_pylon
        q = config.pylon
        rotors = [
            TipRotor(reference,i,wing.pylon_tip_nodes[i];
                transverse_inertia=q.rotor_transverse_inertia,inertia_tensor=q.rotor_inertia_tensor,
                rotor_mass=q.rotor_mass,spin_inertia=q.spin_inertia)
            for i = 1:p.Npropellers
        ]
    else
        rotors = [ReducedRotor(reference,i) for i = 1:p.Npropellers]
    end
    poses = wing_node_poses(wing.problem; nodes=1:wing.wing_node_count)
    aerodynamic =
        beam_pylon ?
        pylon_build_aerodynamic(reference,poses,
            pylon_hub_poses(wing.problem,wing.pylon_tip_nodes),time; tracking_time_steps=false) :
        build_aerodynamic(reference,poses,eta,time; tracking_time_steps=false)
    return ChangAeroBeamsProblem(reference,wing,rotors,aerodynamic,time,eta,
        velocity,zeros(length(eta)),_node_history(wing.problem),config,
        0.0,false,false,Dict{Symbol,Any}())
end

_beam_pylon(coupled) = coupled.config.solver.support == :beam_pylon
_aero_wing_poses(coupled) = wing_node_poses(coupled.wing.problem; nodes=1:coupled.wing.wing_node_count)

function _aero_geometry(coupled,eta,time)
    poses = _aero_wing_poses(coupled)
    if _beam_pylon(coupled)
        hubs = pylon_hub_poses(coupled.wing.problem,coupled.wing.pylon_tip_nodes)
        return pylon_coupled_grids(coupled.reference,poses,hubs,time)
    end
    return coupled_grids(coupled.reference,poses,eta,time)
end

function _aero_load_map(coupled,eta,time,loads,geometry)
    poses = _aero_wing_poses(coupled)
    if _beam_pylon(coupled)
        hubs = pylon_hub_poses(coupled.wing.problem,coupled.wing.pylon_tip_nodes)
        mapped = pylon_aerodynamic_loads(coupled.reference,poses,hubs,time,loads; geometry)
        nodal = zeros(6,coupled.wing.node_count)
        nodal[:,1:coupled.wing.wing_node_count] .= mapped.wing_loads
        for (ip,node) in enumerate(coupled.wing.pylon_tip_nodes)
            nodal[:,node] .+= mapped.pylon_loads[:,ip]
        end
        return merge(mapped,(wing_loads=nodal,))
    end
    return coupled_aerodynamic_loads(coupled.reference,poses,eta,time,loads; geometry)
end

function _new_aerodynamic(coupled,times; eta=coupled.eta)
    poses = _aero_wing_poses(coupled)
    if _beam_pylon(coupled)
        return pylon_build_aerodynamic(coupled.reference,poses,
            pylon_hub_poses(coupled.wing.problem,coupled.wing.pylon_tip_nodes),times;
            tracking_time_steps=false)
    end
    return build_aerodynamic(coupled.reference,poses,eta,times; tracking_time_steps=false)
end

function _special_node(problem,node)
    name = "Chang wing load node $node"
    index = findfirst(s -> any(b -> b.definition.name == name,s.BCs),problem.model.assemblyData.specialNodes)
    isnothing(index) && error("Missing explicit nodal load map for node $node")
    return problem.model.assemblyData.specialNodes[index]
end

function _insert_wrench_residual!(residual,problem,node,wrench)
    s = _special_node(problem,node)
    scale = problem.model.assemblyData.forceScaling
    # v1 assembles one three-component force/moment balance per shared node.
    residual[s.assemblyData.eqs_Fu] .+= wrench[1:3] ./ scale
    residual[s.assemblyData.eqs_Fp] .+= wrench[4:6] ./ scale
    return residual
end

_newmark_rates(q,qn,vn,an,dt) = (
    velocity=2 .* (q .- qn) ./ dt .- vn,
    acceleration=4 .* (q .- qn) ./ dt^2 .- 4 .* vn ./ dt .- an)

function _component_columns(coupled)
    problem = coupled.wing.problem
    columns = Int[]
    for rotor in coupled.rotors
        sn = _special_node(problem,rotor.attachment_node)
        append!(columns,sn.assemblyData.DOF_uF[.!sn.flags.uIsPrescribed])
        append!(columns,sn.assemblyData.DOF_pM[.!sn.flags.pIsPrescribed])
    end
    append!(columns,length(problem.systemArrays.x) .+ (1:length(coupled.eta)))
    return unique(columns)
end

function _component_residual(coupled,z,history,eta_n,velocity_n,acceleration_n,dt;
    static=false,rotor_loads=zeros(length(coupled.eta)))
    problem = coupled.wing.problem
    nx = length(problem.systemArrays.x)
    problem.systemArrays.x .= view(z,1:nx)
    eta = view(z,nx+1:length(z))
    rates =
        static ? (velocity=zeros(length(eta)),acceleration=zeros(length(eta))) :
        _newmark_rates(eta,eta_n,velocity_n,acceleration_n,dt)
    correction = zeros(length(z))
    for rotor in coupled.rotors
        node = rotor.attachment_node
        if static
            poses = wing_node_poses(problem)
            R =
                rotor isa TipRotor ? AB.rotation_tensor_WM(wing_node_parameters(problem).p[:,node]) :
                poses.rotations[node] * CHANG_DIRECTORS_TO_U'
            motion = (velocity=zeros(3),omega=zeros(3),acceleration=zeros(3),alpha=zeros(3))
        else
            motion = wing_node_motion(problem,node,history,dt)
            R = rotor isa TipRotor ? AB.rotation_tensor_WM(motion.p) :
                motion.rotation * CHANG_DIRECTORS_TO_U'
        end
        indices = rotor isa TipRotor ? (1:0) : (2rotor.reference_index-1):(2rotor.reference_index)
        mechanics = rotor_mechanics(rotor,R,motion.velocity,motion.omega,
            motion.acceleration,motion.alpha,eta[indices],
            rates.velocity[indices],rates.acceleration[indices])
        _insert_wrench_residual!(correction,problem,node,mechanics.base_wrench)
        correction[nx.+indices] .= (mechanics.eta_residual .- rotor_loads[indices]) ./
                                   problem.model.assemblyData.forceScaling
    end
    return correction
end

"""Assemble the native beam residual plus moving-base rotor contributions.

Aerodynamic loads remain fixed during a mechanical Newton solve. The native
beam Jacobian is analytic; only the small attachment/support block is centrally
differenced, including finite-rotation spin and acceleration derivatives.
"""
function coupled_mechanical_arrays!(coupled,z,wing_loads,rotor_loads,history,
    eta_n,velocity_n,acceleration_n,dt; static=false,jacobian=true)
    problem = coupled.wing.problem
    nx = length(problem.systemArrays.x)
    problem.systemArrays.x .= z[1:nx]
    problem.systemSolver.state.skipJacobianUpdate = !jacobian
    problem.systemSolver.loadStepState.σ = 1.0
    set_wing_loads!(problem,wing_loads,problem.time.state.timeNow)
    AB.fj!(problem)
    correction = _component_residual(coupled,z,history,eta_n,velocity_n,
        acceleration_n,dt; static,rotor_loads)
    residual = vcat(copy(problem.systemArrays.residual),zeros(length(coupled.eta))) + correction
    jacobian || return residual
    J = zeros(length(z),length(z))
    J[1:nx,1:nx] .= Matrix(problem.systemArrays.jacobian)
    for column in _component_columns(coupled)
        h = cbrt(eps(Float64)) * max(abs(z[column]),1e-3)
        zp,zm = copy(z),copy(z)
        zp[column] += h
        zm[column] -= h
        plus = _component_residual(coupled,zp,history,eta_n,velocity_n,
            acceleration_n,dt; static,rotor_loads)
        minus = _component_residual(coupled,zm,history,eta_n,velocity_n,
            acceleration_n,dt; static,rotor_loads)
        J[:,column] .+= (plus .- minus) ./ (2h)
    end
    problem.systemArrays.x .= z[1:nx]
    return residual,J
end

function _solve_mechanics!(coupled,z,wing_loads,rotor_loads,history,
    eta_n,velocity_n,acceleration_n,dt; static=false)
    tolerance = coupled.config.solver.newton_tolerance
    firstNorm = NaN
    for iteration = 1:coupled.config.solver.maximum_newton_iterations
        r,J = coupled_mechanical_arrays!(coupled,z,wing_loads,rotor_loads,
            history,eta_n,velocity_n,acceleration_n,dt; static)
        residualNorm = norm(r,Inf)
        iteration == 1 && (firstNorm = residualNorm)
        if residualNorm <= tolerance
            coupled.wing.problem.systemArrays.x .= z[1:length(coupled.wing.problem.systemArrays.x)]
            AB.update_states!(coupled.wing.problem)
            return (; z,iterations=iteration,residual=residualNorm)
        end
        all(isfinite,r) || error("Nonfinite augmented mechanical residual")
        increment = -(J \ r)
        all(isfinite,increment) || error("Nonfinite mechanical Newton increment")
        accepted = false
        step = 1.0
        for trial = 1:12
            candidate = z .+ step .* increment
            rc = coupled_mechanical_arrays!(coupled,candidate,wing_loads,rotor_loads,
                history,eta_n,velocity_n,acceleration_n,dt; static,jacobian=false)
            if norm(rc,Inf) < residualNorm || norm(rc,Inf) <= tolerance
                z = candidate
                accepted = true
                break
            end
            step /= 2
        end
        accepted || error("Augmented Newton line search failed: residual=$residualNorm")
    end
    error("Augmented wing/rotor Newton did not converge (initial residual=$firstNorm)")
end

function _loads_error(first,second,chord)
    a = vcat(vec(first.wing_loads[1:3,:]),vec(first.wing_loads[4:6,:]) ./ chord,
        first.rotor_loads ./ chord)
    b = vcat(vec(second.wing_loads[1:3,:]),vec(second.wing_loads[4:6,:]) ./ chord,
        second.rotor_loads ./ chord)
    return norm(a .- b,Inf) / max(norm(a,Inf),norm(b,Inf),1.0)
end

function _interface_velocity(coupled,z,history,eta_n,velocity_n,acceleration_n,dt)
    nx = length(coupled.wing.problem.systemArrays.x)
    coupled.wing.problem.systemArrays.x .= z[1:nx]
    values = Float64[]
    p = coupled.reference.parameters
    for rotor in coupled.rotors
        motion = wing_node_motion(coupled.wing.problem,rotor.attachment_node,history,dt)
        append!(values,motion.velocity ./ max(p.Vinf,1.0))
        append!(values,motion.omega ./ max(abs(p.Ω),1.0))
    end
    if !isempty(coupled.eta)
        rates = _newmark_rates(z[nx+1:end],eta_n,velocity_n,acceleration_n,dt)
        append!(values,rates.velocity ./ max(abs(p.Ω),1.0))
    end
    return values
end

function _evaluate_aero(coupled,z,time)
    problem = coupled.wing.problem
    nx = length(problem.systemArrays.x)
    problem.systemArrays.x .= z[1:nx]
    eta = z[nx+1:end]
    grids,geometry = _aero_geometry(coupled,eta,time)
    loads = UVLM.evaluate_trial!(coupled.aerodynamic; surfaces=grids)
    mapped = _aero_load_map(coupled,eta,time,loads,geometry)
    return mapped,grids
end

function _zero_element_rates!(problem)
    for e in problem.model.assemblyData.elements
        for name in fieldnames(typeof(e.state.rates))
            value = getfield(e.state.rates,name)
            value isa AbstractArray && fill!(value,0.0)
        end
        for name in (:udotEquiv,:pdotEquiv,:VdotEquiv,:ΩdotEquiv,:χdotEquiv)
            fill!(getfield(e.stateCache.equivStatesRates,name),0.0)
        end
    end
    return problem
end

function _settle_fixed_shape(coupled,poses,eta)
    p = coupled.reference.parameters
    settings = coupled.config.initialization
    dt = coupled.time[2] - coupled.time[1]
    nperiod = max(2,round(Int,2pi / (abs(p.Ω) * dt)))
    # Mean windows below use complete discrete revolutions. For exact periodic
    # comparisons use an azimuth-based step dividing 360 degrees (the default).
    isapprox(nperiod * dt * abs(p.Ω),2pi; atol=1e-8,rtol=1e-8) ||
        throw(ArgumentError("Mean initialization needs an integer number of timesteps per rotor revolution"))
    steps = nperiod * settings.maximum_revolutions
    times = collect(0:steps) .* dt
    aero = _new_aerodynamic(coupled,times; eta)
    oldMean = nothing
    oldEndLoads = nothing
    oldGamma = nothing
    cycleWing = zeros(6,coupled.wing.node_count)
    cycleEta = zeros(length(eta))
    periodic_error = Inf
    meanError = endpointError = circulationError = Inf
    for k = 2:length(times)
        UVLM.begin_time_step!(aero)
        grids,geometry = _aero_geometry(coupled,eta,times[k])
        loads = UVLM.evaluate_trial!(aero; surfaces=grids)
        mapped = _aero_load_map(coupled,eta,times[k],loads,geometry)
        cycleWing .+= mapped.wing_loads ./ nperiod
        cycleEta .+= mapped.rotor_loads ./ nperiod
        UVLM.commit_time_step!(aero)
        if mod(k - 1,nperiod) == 0
            meanLoads = (wing_loads=copy(cycleWing),rotor_loads=copy(cycleEta))
            revolution = (k - 1) ÷ nperiod
            if !isnothing(oldMean)
                meanError = _loads_error(meanLoads,oldMean,p.cref)
                endpointError = _loads_error(mapped,oldEndLoads,p.cref)
                circulationError = norm(aero.state.system.Γ .- oldGamma,Inf) / max(norm(aero.state.system.Γ,Inf),1.0)
                periodic_error = max(meanError,endpointError,circulationError)
                if revolution >= settings.minimum_revolutions &&
                   periodic_error <= settings.periodic_load_tolerance
                    return (; loads=meanLoads,aerodynamic=aero,periodic_error,revolutions=revolution)
                end
            end
            oldMean = meanLoads
            oldEndLoads = deepcopy(mapped)
            oldGamma = copy(aero.state.system.Γ)
            fill!(cycleWing,0.0)
            fill!(cycleEta,0.0)
        end
    end
    error(
        "Fixed-shape rotor trim did not converge after $(settings.maximum_revolutions) revolutions: " *
        "mean=$meanError, endpoint=$endpointError, circulation=$circulationError. " *
        "Refine the aerodynamic setup or use :unloaded for a transient; no periodic trim was accepted."
    )
end

"""Solve mean-loaded nonlinear equilibrium and retain a developed rotor wake.

This initializes mean equilibrium, not a fully coupled periodic structural
orbit. Subsequent unperturbed dynamic settling remains necessary for damping
identification. No aerodynamic baseline is subtracted during time marching.
"""
function initialize_mean_equilibrium!(coupled)
    problem = coupled.wing.problem
    settings = coupled.config.initialization
    oldDt = problem.time.state.Δt
    problem.time.state.Δt = Inf
    _zero_element_rates!(problem)
    nx = length(problem.systemArrays.x)
    z = vcat(problem.systemArrays.x,coupled.eta)
    previousLoads = nothing
    records = NamedTuple[]
    converged = false
    finalSettle = nothing
    try
        for iteration = 1:settings.maximum_iterations
            problem.systemArrays.x .= z[1:nx]
            poses = wing_node_poses(problem)
            settled = _settle_fixed_shape(coupled,poses,z[nx+1:end])
            solved = _solve_mechanics!(coupled,copy(z),settled.loads.wing_loads,
                settled.loads.rotor_loads,coupled.node_history,coupled.eta,
                coupled.eta_velocity,coupled.eta_acceleration,Inf; static=true)
            change = norm(solved.z .- z,Inf) / max(norm(solved.z,Inf),1.0)
            load_error =
                isnothing(previousLoads) ? Inf :
                _loads_error(settled.loads,previousLoads,coupled.reference.parameters.cref)
            push!(
                records,
                (; iteration,state_change=change,load_error,
                    periodic_load_error=settled.periodic_error,revolutions=settled.revolutions)
            )
            println("Mean equilibrium $iteration: state=$change, load=$load_error")
            flush(stdout)
            if change <= settings.tolerance && load_error <= settings.tolerance
                z = solved.z
                converged = true
                break
            end
            z .+= settings.relaxation .* (solved.z .- z)
            previousLoads = settled.loads
        end
        converged || error("Mean wing/propeller equilibrium did not converge")
        problem.systemArrays.x .= z[1:nx]
        coupled.eta .= z[nx+1:end]
        (isempty(coupled.eta) || maximum(abs,coupled.eta) <= coupled.config.solver.maximum_relative_angle) ||
            error("Mean equilibrium exceeds the reduced support angle limit")
        fill!(coupled.eta_velocity,0.0)
        fill!(coupled.eta_acceleration,0.0)
        AB.update_states!(problem)
        finalSettle = _settle_fixed_shape(coupled,wing_node_poses(problem),coupled.eta)
    finally
        problem.time.state.Δt = oldDt
    end
    coupled.node_history = _node_history(problem)
    coupled.phase_time_offset = finalSettle.aerodynamic.state.timeNow
    coupled.aerodynamic = _new_aerodynamic(coupled,coupled.time .+ coupled.phase_time_offset)
    UVLM.copy_system!(coupled.aerodynamic.state.system,finalSettle.aerodynamic.state.system)
    coupled.results[:initialization] = records
    coupled.results[:mean_equilibrium] = (state=copy(z),loads=finalSettle.loads,
        periodic_load_error=finalSettle.periodic_error,phase_time_offset=coupled.phase_time_offset)
    coupled.initialized = true
    return coupled
end

function _record!(coupled,index,mapped,iterations,geometryError,load_error,mechanicalError;
    velocity_error=0.0)
    r = coupled.results
    poses = wing_node_poses(coupled.wing.problem)
    r[:node_positions][:,:,index] .= poses.positions
    r[:eta][:,index] .= coupled.eta
    r[:eta_velocity][:,index] .= coupled.eta_velocity
    r[:eta_acceleration][:,index] .= coupled.eta_acceleration
    tip = coupled.wing.wing_node_count
    deformation = poses.rotations[tip] * CHANG_DIRECTORS_TO_U'
    chord = deformation * [1.0,0.0,0.0]
    r[:tip_twist_rad][index] = atan(chord[3],chord[1])
    r[:coupling_iterations][index] = iterations
    r[:geometry_residual][index] = geometryError
    r[:load_residual][index] = load_error
    r[:mechanical_residual][index] = mechanicalError
    r[:velocity_residual][index] = velocity_error
    if _beam_pylon(coupled)
        hubs = pylon_hub_poses(coupled.wing.problem,coupled.wing.pylon_tip_nodes)
        for ip in eachindex(coupled.rotors)
            axis = hubs.rotations[ip] * [1.0,0.0,0.0]
            r[:rotor_angles][2ip-1,index] = atan(-axis[3],axis[1])
            r[:rotor_angles][2ip,index] = atan(-axis[2],axis[1])
        end
    else
        r[:rotor_angles][:,index] .= coupled.eta
    end
    if !isnothing(mapped)
        r[:hub_forces][:,:,index] .= mapped.hub_forces
        r[:hub_moments][:,:,index] .= mapped.hub_moments
    end
    if coupled.config.output.save_structural_history
        AB.save_time_step_data!(coupled.wing.problem)
    end
    return nothing
end

"""Strongly partitioned nonlinear wing/2DOF-rotor/UVLM integration."""
function solve_chang_aerobeams!(coupled)
    coupled.completed && throw(ArgumentError("Build a fresh problem for another response"))
    coupled.aerodynamic.state.stepIndex == 1 || throw(ArgumentError(
        "A partial response cannot be restarted at t=0; build a fresh problem"))
    if !coupled.initialized
        mode = coupled.config.initialization.mode
        mode in (:mean_equilibrium,:unloaded) || throw(ArgumentError("Unknown initialization mode"))
        mode == :mean_equilibrium && initialize_mean_equilibrium!(coupled)
        coupled.initialized = true
    end
    problem = coupled.wing.problem
    nx = length(problem.systemArrays.x)
    count,nodes,np = length(coupled.time),coupled.wing.node_count,length(coupled.rotors)
    merge!(
        coupled.results,
        Dict{Symbol,Any}(
            :time => copy(coupled.time),:node_positions => zeros(3,nodes,count),
            :eta => zeros(length(coupled.eta),count),:eta_velocity => zeros(length(coupled.eta),count),
            :eta_acceleration => zeros(length(coupled.eta),count),:tip_twist_rad => zeros(count),
            :rotor_angles => zeros(2np,count),
            :hub_forces => zeros(3,np,count),:hub_moments => zeros(3,np,count),
            :coupling_iterations => zeros(Int,count),:geometry_residual => zeros(count),
            :load_residual => zeros(count),:mechanical_residual => zeros(count),
            :velocity_residual => zeros(count))
    )
    initialLoads = nothing
    if haskey(coupled.results,:mean_equilibrium)
        _,geometry = _aero_geometry(coupled,coupled.eta,coupled.phase_time_offset)
        initialLoads = _aero_load_map(coupled,coupled.eta,coupled.phase_time_offset,
            UVLM.dimensional_loads(coupled.aerodynamic),geometry)
    end
    _record!(coupled,1,initialLoads,0,0.0,0.0,0.0)
    for index = 2:count
        dt = coupled.time[index] - coupled.time[index-1]
        stateSaved = structural_checkpoint(problem)
        eta_n,velocity_n,acceleration_n =
            copy(coupled.eta),copy(coupled.eta_velocity),copy(coupled.eta_acceleration)
        history = coupled.node_history
        AB.update_time_state!(problem,index)
        AB.update_equivalent_states_rates!(problem)
        z = vcat(problem.systemArrays.x,eta_n .+ dt .* velocity_n .+ 0.5dt^2 .* acceleration_n)
        # Predict interface nodal displacement/rotation without replacing the
        # native beam states or their accepted equivalent-rate history.
        for node = 2:nodes
            sn = _special_node(problem,node)
            z[sn.assemblyData.DOF_uF] .= history.u[:,node] .+ dt .* history.udot[:,node] .+ 0.5dt^2 .* history.uddot[:,node]
            z[sn.assemblyData.DOF_pM] .= history.p[:,node] .+ dt .* history.pdot[:,node] .+ 0.5dt^2 .* history.pddot[:,node]
        end
        time = coupled.phase_time_offset + coupled.time[index]
        UVLM.begin_time_step!(coupled.aerodynamic)
        oldLoads = nothing
        converged = false
        acceptedLoads = nothing
        geometryError = velocity_error = load_error = mechanicalError = Inf
        iterations = 0
        try
            for iteration = 1:coupled.config.solver.maximum_coupling_iterations
                iterations = iteration
                mapped,grids = _evaluate_aero(coupled,z,time)
                rotor_loads = copy(mapped.rotor_loads)
                wing_loads = copy(mapped.wing_loads)
                settings = coupled.config.solver
                phase = (coupled.time[index] - settings.impulse_start_s) / settings.impulse_duration_s
                if 0 < phase < 1
                    pulse = settings.impulse_magnitude_nm * sinpi(phase)^2
                    if _beam_pylon(coupled)
                        ip = settings.impulse_rotor
                        hubs = pylon_hub_poses(problem,coupled.wing.pylon_tip_nodes)
                        wing_loads[4:6,coupled.wing.pylon_tip_nodes[ip]] .+=
                            pulse .* hubs.rotations[ip][:,2]
                    else
                        rotor_loads[2settings.impulse_rotor-1] += pulse
                    end
                end
                solved = _solve_mechanics!(coupled,copy(z),wing_loads,rotor_loads,
                    history,eta_n,velocity_n,acceleration_n,dt)
                problem.systemArrays.x .= solved.z[1:nx]
                candidate_grids,_ = _aero_geometry(coupled,solved.z[nx+1:end],time)
                geometryError =
                    maximum(maximum(abs.(candidate_grids[i] .- grids[i])) for i in eachindex(grids)) /
                    coupled.reference.parameters.cref
                velocity_error = norm(
                    _interface_velocity(coupled,solved.z,history,eta_n,velocity_n,acceleration_n,dt) .-
                    _interface_velocity(coupled,z,history,eta_n,velocity_n,acceleration_n,dt),Inf)
                load_error =
                    isnothing(oldLoads) ? Inf : _loads_error(mapped,oldLoads,coupled.reference.parameters.cref)
                mechanicalError = solved.residual
                if geometryError <= settings.geometry_tolerance && velocity_error <= settings.velocity_tolerance &&
                   load_error <= settings.load_tolerance
                    finalLoads,_ = _evaluate_aero(coupled,solved.z,time)
                    finalError = _loads_error(finalLoads,mapped,coupled.reference.parameters.cref)
                    finalWing = wing_loads .+ (finalLoads.wing_loads .- mapped.wing_loads)
                    finalRotor = rotor_loads .+ (finalLoads.rotor_loads .- mapped.rotor_loads)
                    finalR = coupled_mechanical_arrays!(coupled,solved.z,finalWing,finalRotor,
                        history,eta_n,velocity_n,acceleration_n,dt; jacobian=false)
                    finalMechanical = norm(finalR,Inf)
                    if finalError <= settings.load_tolerance &&
                       finalMechanical <= settings.coupled_equilibrium_tolerance
                        z = solved.z
                        acceptedLoads = finalLoads
                        load_error = max(load_error,finalError)
                        mechanicalError = finalMechanical
                        converged = true
                        break
                    end
                end
                oldLoads = mapped
                z .+= settings.coupling_relaxation .* (solved.z .- z)
            end
            converged || error("FSI failed at t=$(coupled.time[index]): geometry=$geometryError, load=$load_error")
            (isempty(coupled.eta) || maximum(abs,z[nx+1:end]) <= coupled.config.solver.maximum_relative_angle) ||
                error("Relative support angles exceed the reduced model's configured limit")
            problem.systemArrays.x .= z[1:nx]
            AB.update_states!(problem)
            coupled.eta .= z[nx+1:end]
            erates = _newmark_rates(coupled.eta,eta_n,velocity_n,acceleration_n,dt)
            coupled.eta_velocity .= erates.velocity
            coupled.eta_acceleration .= erates.acceleration
            q = wing_node_parameters(problem)
            ur = _newmark_rates(q.u,history.u,history.udot,history.uddot,dt)
            pr = _newmark_rates(q.p,history.p,history.pdot,history.pddot,dt)
            coupled.node_history = (; u=copy(q.u),p=copy(q.p),udot=ur.velocity,pdot=pr.velocity,
                uddot=ur.acceleration,pddot=pr.acceleration)
            # Commit the accepted aerodynamic trial/wakes exactly once.
            UVLM.commit_time_step!(coupled.aerodynamic)
        catch
            UVLM.rollback_time_step!(coupled.aerodynamic)
            restore_structural_checkpoint!(problem,stateSaved)
            coupled.eta .= eta_n
            coupled.eta_velocity .= velocity_n
            coupled.eta_acceleration .= acceleration_n
            coupled.node_history = history
            rethrow()
        end
        _record!(coupled,index,acceptedLoads,iterations,geometryError,load_error,mechanicalError; velocity_error)
        if mod(index - 1,coupled.config.solver.progress_frequency) == 0 || index == count
            println(
                "Chang AeroBeams: t=$(coupled.time[index]), FSI=$iterations, geometry=$geometryError, load=$load_error"
            )
            flush(stdout)
        end
    end
    coupled.completed = true
    return coupled
end

function save_chang_aerobeams(coupled)
    coupled.completed || throw(ArgumentError("Only a completed analysis can be saved"))
    output,r = coupled.config.output,coupled.results
    mkpath(output.directory)
    prefix = joinpath(output.directory,output.label)
    columns = ["time_s","tip_x_m","tip_y_m","tip_z_m","tip_twist_deg",
        "fsi_iterations","geometry_residual","velocity_residual","load_residual","mechanical_residual"]
    tip = coupled.wing.wing_node_count
    data = hcat(r[:time],permutedims(r[:node_positions][:,tip,:]),
        rad2deg.(r[:tip_twist_rad]),r[:coupling_iterations],r[:geometry_residual],r[:velocity_residual],
        r[:load_residual],r[:mechanical_residual])
    for rotor in eachindex(coupled.rotors)
        append!(
            columns,
            ["rotor$(rotor)_pitch_deg","rotor$(rotor)_yaw_deg",
                "rotor$(rotor)_force_x_N","rotor$(rotor)_force_y_N","rotor$(rotor)_force_z_N"]
        )
        data = hcat(data,rad2deg.(r[:rotor_angles][2rotor-1,:]),rad2deg.(r[:rotor_angles][2rotor,:]),
            permutedims(r[:hub_forces][:,rotor,:]))
    end
    open(prefix * "_history.csv","w") do io
        println(io,join(columns,','))
        writedlm(io,data,',')
    end
    open(prefix * "_settings.txt","w") do io
        println(io,"Nonlinear AeroBeams Chang wing + support=$(coupled.config.solver.support) + interacting UVLM")
        println(io,"Completed: ",coupled.completed)
        println(
            io,
            "Reference configuration (legacy integration, coupling, excitation and output controls are not used by this driver): ",
            coupled.reference.config
        )
        println(
            io,
            "Actual integration: AeroBeams first-order trapezoidal; support/attachment Newmark gamma=0.5, beta=0.25"
        )
        println(io,"Actual coupling: strong partitioned UVLM with augmented nonlinear mechanical Newton")
        println(io,"Applied perturbation uses solver.impulse_* settings below")
        println(io,"Coupling settings: ",coupled.config.solver)
        println(io,"Initialization: ",coupled.config.initialization)
        println(io,"Wing representation: ",coupled.wing.metadata)
        println(io,"Rotor phase time offset: ",coupled.phase_time_offset)
        println(io,"Final relative time: ",last(coupled.time))
        println(io,"Timestep: ",coupled.time[2] - coupled.time[1])
        println(io,"Maximum geometry residual: ",maximum(r[:geometry_residual]))
        println(io,"Maximum interface velocity residual: ",maximum(r[:velocity_residual]))
        println(io,"Maximum load residual: ",maximum(r[:load_residual]))
        println(io,"Maximum coupled mechanical residual: ",maximum(r[:mechanical_residual]))
        println(io,"Mean aerodynamic loads retained. Support: ",
            _beam_pylon(coupled) ? "nonlinear pylon beam; rotor inertia assumptions recorded above" :
            "small-relative-deformation reduced model on finite-rotation wing")
        if haskey(r,:mean_equilibrium)
            println(io,"Mean equilibrium iterations: ",length(r[:initialization]))
            println(io,"Initial tip z [m]: ",r[:node_positions][3,tip,1])
            println(io,"Periodic aerodynamic initialization residual: ",r[:mean_equilibrium].periodic_load_error)
        end
    end
    if haskey(r,:initialization)
        open(prefix * "_initialization.csv","w") do io
            println(io,"iteration,state_change,load_error,periodic_load_error,revolutions")
            for row in r[:initialization]
                println(io,join(values(row),','))
            end
        end
    end
    if output.plot_results
        P = AB.Plots
        wingplot =
            P.plot(r[:time],r[:node_positions][3,tip,:],xlabel="Time [s]",ylabel="Wing tip z [m]",label=false)
        twistplot =
            P.plot(r[:time],rad2deg.(r[:tip_twist_rad]),xlabel="Time [s]",ylabel="Wing tip twist [deg]",label=false)
        pitchplot = P.plot(
            r[:time],
            permutedims(rad2deg.(r[:rotor_angles][1:2:end,:])),
            xlabel="Time [s]",
            ylabel="Rotor pitch [deg]",
            label=false
        )
        yawplot = P.plot(
            r[:time],
            permutedims(rad2deg.(r[:rotor_angles][2:2:end,:])),
            xlabel="Time [s]",
            ylabel="Rotor yaw [deg]",
            label=false
        )
        P.savefig(
            P.plot(wingplot,twistplot,pitchplot,yawplot; layout=(2,2),size=(1000,700)),
            prefix * "_history.png"
        )
    end
    return (; history=prefix * "_history.csv",settings=prefix * "_settings.txt")
end

function run_chang_aerobeams(config=chang_aerobeams_defaults())
    coupled = build_chang_aerobeams(config)
    solve_chang_aerobeams!(coupled)
    config.output.write_results && (coupled.results[:paths] = save_chang_aerobeams(coupled))
    return coupled
end
