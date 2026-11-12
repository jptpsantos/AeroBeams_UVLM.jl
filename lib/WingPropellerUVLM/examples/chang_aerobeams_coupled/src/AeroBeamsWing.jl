# A physical geometrically nonlinear AeroBeams wing. All public poses, spatial
# velocities and wrenches use the UVLM frame U = (chord aft, span, up).
# The beam's right-handed section frame is b = (span, chord aft, down).
const CHANG_WING_NAME = "Chang nonlinear wing"
const CHANG_WING_R0 = [0.0 1.0 0.0; 1.0 0.0 0.0; 0.0 0.0 -1.0]
# The reference is a pi rotation about (1,1,0)/sqrt(2). Specify its WM
# parameters directly: v1's Euler-to-quaternion conversion can round the
# quaternion vector norm above one at this exact half-turn.
const CHANG_WING_P0 = [2sqrt(2.0),2sqrt(2.0),0.0]

_wing_load_bc_name(node) = "Chang wing load node $node"

function _wing_control_lengths(nodes)
    length(nodes) >= 2 || throw(ArgumentError("At least two wing nodes are required"))
    all(>(0),diff(nodes)) || throw(ArgumentError("Wing nodes must increase strictly"))
    return [diff(nodes)[1] / 2;
        (nodes[3:end] .- nodes[1:end-2]) ./ 2;
        diff(nodes)[end] / 2]
end

"""Integrate the complete Chang spatial-block density over beam elements.

Chang's conservative nodal blocks are converted to a piecewise-linear density
using nodal control lengths. The integral preserves the sum of all 6-by-6 local
spatial blocks, including all first moments and products of inertia. This is a
distributed reconstruction of the same data; it does not reproduce Chang's
nodal-lumped mass matrix. Its modal discretization must be compared explicitly.
"""
function _wing_element_inertias(source_nodes,source_blocks,target_nodes)
    length(source_nodes) == length(source_blocks) || throw(DimensionMismatch(
        "One spatial inertia block is required at every source node"))
    controlLengths = _wing_control_lengths(source_nodes)
    density = [Matrix{Float64}(block) ./ width
               for (block,width) in zip(source_blocks,controlLengths)]
    result = Matrix{Float64}[]
    for element = 1:(length(target_nodes)-1)
        left,right = target_nodes[element],target_nodes[element+1]
        integral = zeros(6,6)
        for segment = 1:(length(source_nodes)-1)
            a = max(left,source_nodes[segment])
            b = min(right,source_nodes[segment+1])
            b <= a && continue
            width = source_nodes[segment+1] - source_nodes[segment]
            ta = (a - source_nodes[segment]) / width
            tb = (b - source_nodes[segment]) / width
            wl = width * ((tb - tb^2 / 2) - (ta - ta^2 / 2))
            wr = width * (tb^2 - ta^2) / 2
            integral .+= wl .* density[segment] .+ wr .* density[segment+1]
        end
        inertia = integral ./ (right - left)
        inertia = 0.5 .* (inertia .+ inertia')
        isposdef(Symmetric(inertia)) || error(
            "Reconstructed AeroBeams wing sectional inertia $element is not positive definite")
        push!(result,inertia)
    end
    integrated = sum(result[e] .* (target_nodes[e+1] - target_nodes[e])
                     for e in eachindex(result))
    original = sum(source_blocks)
    relative_error = norm(integrated - original,Inf) / max(norm(original,Inf),1.0)
    relative_error <= 1e-12 || error("Wing spatial-block conservation failed: $relative_error")
    return result,relative_error
end

# AeroBeams has constant compliance per element. Average compliance across every
# tabulated material interval, splitting at the sharp Chang transition stations.
# This preserves the integrated local force/strain law for constant resultants.
function _wing_element_stiffness(parameters,nodes; shear_stiffness=1e12)
    isfinite(shear_stiffness) && shear_stiffness > 0 ||
        throw(ArgumentError("shear_stiffness must be finite and positive"))
    distribution = parameters.wing_stiffness_distribution
    stations = distribution.eta .* parameters.span_length
    gauss_x = (-sqrt(3 / 5),0.0,sqrt(3 / 5))
    gaussW = (5 / 9,8 / 9,5 / 9)
    result = Matrix{Float64}[]
    for element = 1:(length(nodes)-1)
        left,right = nodes[element],nodes[element+1]
        edges = sort!(unique([left; filter(x -> left < x < right,stations); right]))
        compliance = zeros(6,6)
        for interval = 1:(length(edges)-1)
            width = edges[interval+1] - edges[interval]
            center = (edges[interval+1] + edges[interval]) / 2
            for (xi,weight) in zip(gauss_x,gaussW)
                eta = (center + width * xi / 2) / parameters.span_length
                D = Chang.chang_constitutive_matrix_at_eta(distribution,eta)
                S = zeros(6,6)
                S[1,1] = D[1,1]
                S[2,2] = S[3,3] = shear_stiffness
                S[4,4] = D[4,4]
                # With b=(span,chord,down), kappa_2=-w_down'' and
                # kappa_3=v_chord''. Transform Chang's [v'',w''] law.
                S[5,5],S[6,6] = D[3,3],D[2,2]
                S[5,6] = S[6,5] = -D[2,3]
                compliance .+= (weight * width / 2) .* inv(Symmetric(S))
            end
        end
        push!(result,Matrix(inv(Symmetric(compliance ./ (right - left)))))
    end
    return result
end

"""Build and initialize an unloaded nonlinear AeroBeams Chang wing.

The returned `problem` contains the actual AeroBeams mixed beam states and
nonlinear residual, not a reduced matrix surrogate. Rotor inertia belongs to
the separate hub subsystem and is intentionally absent here. `time_vector`
must be shared with the coupled driver. Large, finite shear stiffness provides
the Euler-Bernoulli limit; it is an explicit approximation to test in refinement.
"""
function build_wing(reference; time_vector=Float64.(reference.parameters.t),
    n_elements=reference.parameters.Ne,shear_stiffness=1e12,
    newton_maximum_iterations=35,newton_absolute_tolerance=1e-9,
    newton_relative_tolerance=1e-9,newton_display_iterations=false,
    tracking_time_steps=false)
    parameters = reference.parameters
    n_elements isa Integer && n_elements >= 1 || throw(ArgumentError(
        "n_elements must be a positive integer"))
    times = Float64.(time_vector)
    length(times) >= 2 && all(isfinite,times) && all(>(0),diff(times)) ||
        throw(ArgumentError("Provide at least two finite increasing time samples"))
    nodes =
        n_elements == parameters.Ne ? Float64.(parameters.span_nodes) :
        collect(range(0.0,parameters.span_length; length=n_elements + 1))
    elastic_axis_fraction = reference.config.aerodynamic.elastic_axis_fraction
    ea_x = parameters.xle_distribution .+ elastic_axis_fraction .* parameters.chord
    maximum(abs.(ea_x .- first(ea_x))) <= 1e-10 || throw(ArgumentError(
        "This Chang nonlinear-wing adapter currently requires a constant-x elastic axis"))
    stiffness = _wing_element_stiffness(parameters,nodes; shear_stiffness)
    inertia,conservationError = _wing_element_inertias(
        Float64.(parameters.span_nodes),parameters.spatial_inertia_node_blocks,nodes)
    beam = AB.create_Beam(name=CHANG_WING_NAME,
        length=parameters.span_length,
        nElements=n_elements,
        normalizedNodalPositions=nodes ./ parameters.span_length,
        S=stiffness,
        I=inertia,
        p0=copy(CHANG_WING_P0),
        firstNodeOffset=[first(ea_x),0.0,0.0],
        surface=nothing)
    isapprox(beam.geometry.R0,CHANG_WING_R0; atol=1e-13,rtol=0) ||
        error("Unexpected Chang AeroBeams reference frame")
    clamp = AB.create_BC(name="Chang wing root clamp",
        beam=beam,
        node=1,
        conditions=Dict(Symbol.(["u1A","u2A","u3A","p1A","p2A","p3A"]) .=> zeros(6)))
    bcs = [clamp]
    for node in eachindex(nodes)
        push!(
            bcs,
            AB.create_BC(name=_wing_load_bc_name(node),
                beam=beam,
                node=node,
                conditions=Dict(Symbol.(["F1A","F2A","F3A","M1A","M2A","M3A"]) .=> zeros(6)))
        )
    end
    model = AB.create_Model(name="Chang wing with UVLM and hub subsystem",
        beams=[beam],BCs=bcs,gravityVector=zeros(3),v_A=t -> zeros(3),
        units=AB.create_UnitsSystem(frequency="Hz"))
    newton = AB.create_NewtonRaphson(maxIter=newton_maximum_iterations,
        absTol=newton_absolute_tolerance,
        relTol=newton_relative_tolerance,
        displayStatus=newton_display_iterations,
        method=:Newton,
        trackLoadSteps=false)
    problem = AB.create_DynamicProblem(model=model,
        timeVector=times,
        systemSolver=newton,
        trackingTimeSteps=tracking_time_steps,
        displayProgress=false)
    prepare_structural_problem!(problem)
    AB.solve_initial_dynamic!(problem)
    # Check initialized structural arrays independently of the startup solver flag.
    all(isfinite,problem.systemArrays.x) || error("Nonfinite unloaded Chang wing initialization")
    poses = wing_node_poses(problem)
    maximum(abs.(poses.positions[2,:] .- nodes)) <= 1e-10 ||
        error("Chang AeroBeams span stations do not match the aerodynamic frame")
    mass = sum(inertia[e][1,1] * diff(nodes)[e] for e in eachindex(inertia))
    metadata = (; mass_kg=mass,spatial_block_conservation_error=conservationError,
        inertia_representation=:conservative_distributed_spatial_density,
        stiffness_representation=:element_mean_compliance,
        shear_stiffness,section_frame="span, chord aft, down",
        global_frame="chord aft, span, up",structural_damping=:none,
        comparison_note="Distributed AeroBeams inertia differs from Chang nodal lumping; compare modes and refine.")
    return (; problem,node_count=length(nodes),wing_node_count=length(nodes),span_nodes=nodes,
        load_bc_indices=collect(2:(length(nodes)+1)),beam_name=CHANG_WING_NAME,
        reference_rotation=copy(CHANG_WING_R0),metadata)
end

function _wing_beam(problem)
    index = findfirst(beam -> beam.name == CHANG_WING_NAME,problem.model.beams)
    isnothing(index) && error("The active AeroBeams model has no Chang nonlinear wing")
    return problem.model.beams[index]
end

function _wing_special_node(problem,node)
    1 <= node <= problem.model.assemblyData.nNodesTotal || throw(BoundsError(1:problem.model.assemblyData.nNodesTotal,node))
    index = findfirst(s -> s.assemblyData.globalID == node,problem.model.assemblyData.specialNodes)
    isnothing(index) && error("Wing node $node has no explicit special-node coordinates")
    return problem.model.assemblyData.specialNodes[index]
end

"""Exact nodal displacements and WM parameters from the active mixed unknowns."""
function wing_node_parameters(problem)
    count = problem.model.assemblyData.nNodesTotal
    u,p = zeros(3,count),zeros(3,count)
    for node = 1:count
        special = _wing_special_node(problem,node)
        # This adapter prescribes an exactly fixed root. Its unknowns are
        # reaction loads and must never be interpreted as displacements.
        node == 1 && continue
        u[:,node] .= problem.systemArrays.x[special.assemblyData.DOF_uF]
        p[:,node] .= problem.systemArrays.x[special.assemblyData.DOF_pM]
    end
    return (; u,p)
end

function set_wing_node_parameters!(problem,u,p)
    count = problem.model.assemblyData.nNodesTotal
    size(u) == size(p) == (3,count) || throw(DimensionMismatch(
        "Wing u and p must have size (3,$count)"))
    for node = 2:count
        special = _wing_special_node(problem,node)
        problem.systemArrays.x[special.assemblyData.DOF_uF] .= u[:,node]
        problem.systemArrays.x[special.assemblyData.DOF_pM] .= p[:,node]
    end
    return problem
end

"""Reference-axis positions and full section rotations b -> U, from problem.systemArrays.x."""
function wing_node_poses(problem; nodes=1:problem.model.assemblyData.nNodesTotal)
    parameters = wing_node_parameters(problem)
    count = length(nodes)
    positions = zeros(3,count)
    rotations = Matrix{Float64}[]
    for (j,node) in enumerate(nodes)
        index = findfirst(e -> node in e.assemblyData.nodesGlobalIDs,problem.model.assemblyData.elements)
        element = problem.model.assemblyData.elements[index]
        firstSide = element.assemblyData.nodesGlobalIDs[1] == node
        referencePosition = firstSide ? element.assemblyData.r_n1 : element.assemblyData.r_n2
        reference_rotation = firstSide ? element.geometry.R0_n1 : element.geometry.R0_n2
        positions[:,j] .= referencePosition .+ parameters.u[:,node]
        R = AB.rotation_tensor_WM(parameters.p[:,node])
        push!(rotations,Matrix(R * reference_rotation))
    end
    return (; positions,rotations)
end

"""Special-node columns and physical force/moment balance rows for augmentation."""
function wing_node_dof_map(problem,node)
    special = _wing_special_node(problem,node)
    force_rows = collect(special.assemblyData.eqs_Fu)
    moment_rows = collect(special.assemblyData.eqs_Fp)
    return (; special_node=special,prescribed=node == 1,
        displacement_columns=node == 1 ? Int[] : collect(special.assemblyData.DOF_uF),
        rotation_columns=node == 1 ? Int[] : collect(special.assemblyData.DOF_pM),
        columns=node == 1 ? Int[] : vcat(special.assemblyData.DOF_uF,special.assemblyData.DOF_pM),
        force_rows,moment_rows,
        equilibrium_rows=vcat(force_rows,moment_rows),
        force_scaling=problem.model.assemblyData.forceScaling)
end

"""Apply numeric spatial nodal wrenches [Fx,Fy,Fz,Mx,My,Mz] in frame U."""
function set_wing_loads!(problem,loads,time=problem.time.state.timeNow)
    count = problem.model.assemblyData.nNodesTotal
    size(loads) == (6,count) || throw(DimensionMismatch(
        "Wing wrenches must have size (6,$count)"))
    all(isfinite,loads) || throw(ArgumentError("Nonfinite wing wrench"))
    for node = 1:count
        index = findfirst(bc -> bc.definition.name == _wing_load_bc_name(node),problem.model.BCs)
        isnothing(index) && error("Missing numeric load BC at wing node $node")
        bc = problem.model.BCs[index]
        set_nodal_wrench!(bc,view(loads,:,node),time)
    end
    return problem
end

_wing_vee_skew(A) = [A[3,2] - A[2,3],A[1,3] - A[3,1],
    A[2,1] - A[1,2]] ./ 2

"""Spatial angular velocity/acceleration of an exact WM rotation.

`pdot` and `pddot` are parameter rates, not angular vectors. Differentiating
R(p) gives omega=vee(Rdot*R') and alpha=vee(Rddot*R'+Rdot*Rdot').
"""
function wm_spatial_angular_kinematics(p,pdot,pddot)
    R = AB.rotation_tensor_WM(p)
    rotationPath = t -> vec(AB.rotation_tensor_WM(
        p .+ t .* pdot .+ (t^2 / 2) .* pddot))
    Rdot = reshape(AB.ForwardDiff.derivative(rotationPath,0.0),3,3)
    Rddot = reshape(AB.ForwardDiff.derivative(
            t -> AB.ForwardDiff.derivative(rotationPath,t),0.0),3,3)
    omega = _wing_vee_skew(Rdot * R')
    alpha = _wing_vee_skew(Rddot * R' + Rdot * Rdot')
    return (; omega,alpha)
end

"""Node motion from accepted nodal parameter history and Newmark (.5,.25).

`history` contains six 3-by-N arrays u,p,udot,pdot,uddot,pddot at t_n.
No accepted history is changed by this trial evaluation. In particular a second
FSI or Newton trial at t_(n+1) must use the same history and dt.
"""
function wing_node_motion(problem,node,history,dt)
    isfinite(dt) && dt > 0 || throw(ArgumentError("dt must be finite and positive"))
    parameters = wing_node_parameters(problem)
    poses = wing_node_poses(problem)
    u,p = parameters.u[:,node],parameters.p[:,node]
    udot = 2 / dt .* (u .- history.u[:,node]) .- history.udot[:,node]
    pdot = 2 / dt .* (p .- history.p[:,node]) .- history.pdot[:,node]
    uddot = 2 / dt .* (udot .- history.udot[:,node]) .- history.uddot[:,node]
    pddot = 2 / dt .* (pdot .- history.pdot[:,node]) .- history.pddot[:,node]
    angular = wm_spatial_angular_kinematics(p,pdot,pddot)
    return (; position=poses.positions[:,node],rotation=poses.rotations[node],
        velocity=udot,omega=angular.omega,acceleration=uddot,alpha=angular.alpha,
        u,p,udot,pdot,uddot,pddot)
end

"""Initial zero-rate history for an unloaded stationary wing."""
function initial_wing_node_history(problem)
    parameters = wing_node_parameters(problem)
    return (; parameters.u,parameters.p,udot=zero(parameters.u),pdot=zero(parameters.p),
        uddot=zero(parameters.u),pddot=zero(parameters.p))
end
