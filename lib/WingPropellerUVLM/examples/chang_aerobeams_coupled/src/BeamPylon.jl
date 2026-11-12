"""Build nonlinear pylon beams connected to the physical Chang wing nodes.

The default EI=K_relative*L matches a cantilever's tip rotational compliance
under a pure moment. It does not assert equivalence of the distributed pylon
modes or effective inertia to Chang's calibrated two-coordinate support.
Provide measured sectional properties for a physical comparison.
"""
function build_wing_with_beam_pylons(reference; time_vector,pylon_elements=4,
    EI_pitch=nothing,EI_yaw=nothing,EA=1e10,GA=1e10,GJ=1e8,
    rotary_inertia_per_length=1e-3,kwargs...)
    pylon_elements >= 1 || throw(ArgumentError("At least one pylon element required"))
    base = build_wing(reference; time_vector,kwargs...)
    wing = _wing_beam(base.problem)
    p = reference.parameters
    L = p.L_pylon
    L > 0 || throw(ArgumentError("Positive physical pylon length required"))
    pitch = isnothing(EI_pitch) ? p.Kθ_prop * L : EI_pitch
    yaw = isnothing(EI_yaw) ? p.Kψ_prop * L : EI_yaw
    all(x -> isfinite(x) && x > 0,(pitch,yaw,EA,GA,GJ,rotary_inertia_per_length)) ||
        throw(ArgumentError("Pylon sectional properties must be finite and positive"))
    # Local1=-U-x runs from wing to upstream rotor; local2=U-y,
    # local3=-U-z. Local2 bending is physical rotor pitch.
    section = Matrix(Diagonal([EA,GA,GA,GJ,pitch,yaw]))
    mu = reference.config.structural.pylon_mass_per_length_kgpm
    mu > 0 || throw(ArgumentError("Positive distributed pylon mass required"))
    inertia = Matrix(Diagonal([mu,mu,mu,rotary_inertia_per_length,
        rotary_inertia_per_length,rotary_inertia_per_length]))
    beams = AB.Beam[wing]
    connections = Tuple{Pair{AB.Beam,AB.Beam},Pair{Int,Int}}[]
    for ip = 1:p.Npropellers
        beam = AB.create_Beam(name="Chang nonlinear pylon $ip",
            length=L,
            nElements=pylon_elements,
            S=[section],
            I=[inertia],
            p0=[0.0,4.0,0.0], # Exact pi rotation about U-y, in WM parameters.
            surface=nothing)
        push!(beams,beam)
        # v1 attaches the current (later) beam to an already assembled beam.
        push!(connections,(beam => wing,1 => p.prop_attach_nodes[ip]))
    end
    # Resolve model global IDs before constructing explicitly mapped nodal BCs.
    preliminary = AB.create_Model(name="Chang nonlinear wing and beam pylons",
        beams=beams,connections=connections,BCs=AB.BC[],gravityVector=zeros(3),v_A=t -> zeros(3))
    clamp = AB.create_BC(name="Chang wing root clamp",
        beam=wing,
        node=1,
        conditions=Dict(Symbol.(["u1A","u2A","u3A","p1A","p2A","p3A"]) .=> zeros(6)))
    bcs = AB.BC[clamp]
    for node = 1:preliminary.assemblyData.nNodesTotal
        owner = findfirst(b -> node in b.assemblyData.nodesGlobalIDs,beams)
        beam = beams[owner]
        localnode = findfirst(==(node),beam.assemblyData.nodesGlobalIDs)
        push!(
            bcs,
            AB.create_BC(name=_wing_load_bc_name(node),
                beam=beam,
                node=localnode,
                conditions=Dict(Symbol.(["F1A","F2A","F3A","M1A","M2A","M3A"]) .=> zeros(6)))
        )
    end
    model = AB.create_Model(name="Chang nonlinear wing and beam pylons",
        beams=beams,connections=connections,BCs=bcs,gravityVector=zeros(3),v_A=t -> zeros(3),
        units=AB.create_UnitsSystem(frequency="Hz"))
    problem = AB.create_DynamicProblem(model=model,
        timeVector=Float64.(time_vector),
        systemSolver=AB.create_NewtonRaphson(trackLoadSteps=false),
        trackingTimeSteps=false,
        displayProgress=false)
    prepare_structural_problem!(problem)
    AB.solve_initial_dynamic!(problem)
    tips = [beam.elements[end].assemblyData.nodesGlobalIDs[2] for beam in model.beams[2:end]]
    metadata = merge(
        base.metadata,
        (pylon_elements,EI_pitch=pitch,EI_yaw=yaw,EA,GA,GJ,
            pylon_mass_per_length=mu,pylon_rotary_inertia_per_length=rotary_inertia_per_length,
            pylon_calibration=:tip_moment_rotational_compliance,
            pylon_note="Reduced support mass/stiffness removed; physical rotor inertia supplied once at tip.")
    )
    return merge(
        base,
        (; problem,node_count=model.assemblyData.nNodesTotal,wing_node_count=base.node_count,
            pylon_tip_nodes=tips,metadata)
    )
end

"""Physical rotor hub positions and incremental pylon-tip rotations in U."""
function pylon_hub_poses(problem,tipNodes)
    parameters = wing_node_parameters(problem)
    poses = wing_node_poses(problem; nodes=tipNodes)
    rotations = [Matrix(AB.rotation_tensor_WM(parameters.p[:,node])) for node in tipNodes]
    return (; positions=poses.positions,rotations)
end
