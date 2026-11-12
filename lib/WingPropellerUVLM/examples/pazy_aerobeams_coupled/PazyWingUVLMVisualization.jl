# Included after PazyWingUVLMCoupling.jl; uses its AeroBeams and UVLM imports.
const Plots = AeroBeams.Plots

const DEFAULT_OUTPUT_DIRECTORY = joinpath(@__DIR__,"output")

function path_for_aerobeams(absolutePath::String)
    root = abspath(REPOSITORY_ROOT)
    path = abspath(absolutePath)
    rootKey = replace(root,'\\' => '/')
    pathKey = replace(path,'\\' => '/')
    if Sys.iswindows()
        rootKey,pathKey = lowercase(rootKey),lowercase(pathKey)
    end
    (pathKey == rootKey || startswith(pathKey,rootKey * "/")) ||
        throw(ArgumentError("AeroBeams animations must be saved inside the repository"))
    relativePath = relpath(path,root)
    (isabspath(relativePath) || first(splitpath(relativePath)) == "..") &&
        throw(ArgumentError("AeroBeams animations must be saved inside the repository"))
    return "/" * replace(relativePath,'\\' => '/')
end

function maximum_aerodynamic_force(result)
    history = result.aerodynamic_nodal_load_history
    maximumForce = 0.0
    for step in axes(history,3),node in axes(history,2)
        force = history[1:3,node,step]
        maximumForce = max(maximumForce,sqrt(sum(abs2,force)))
    end
    return maximumForce
end

function force_history_index(result,time)
    return clamp(round(Int,time / result.dt) + 1,1,length(result.time))
end

function beam_node_positions(result,time_index)
    problem = result.structural_problem
    isempty(problem.solutionHistory.states) && throw(ArgumentError(
        "Enable save_structural_history when generating structural animations"))
    elements = problem.model.assemblyData.elements
    states = problem.solutionHistory.states[time_index]
    positions = zeros(3,length(elements) + 1)
    positions[:,1] = elements[1].assemblyData.r_n1
    for (elementIndex,element) in enumerate(elements)
        positions[:,elementIndex+1] =
            element.assemblyData.r_n2 + states[element.assemblyData.globalID].nodal.u2
    end
    return positions
end

"""Create the usual AeroBeams nonlinear structural-deformation GIF."""
function save_structural_animation(result;
    output_path::String=joinpath(DEFAULT_OUTPUT_DIRECTORY,"pazy_structure.gif"),
    fps::Int=30,deformation_scale::Real=1.0,backend::Symbol=:gr,camera=(45,45),
    surface_alpha::Real=0.5,elastic_axis_linewidth::Real=2.5,
    show_clamped_root::Bool=true,show_force_vectors::Bool=false,
    show_moment_vectors::Bool=false,force_vector_scale::Real=1.0)

    problem = result.structural_problem
    isempty(problem.solutionHistory.states) && throw(ArgumentError(
        "Enable save_structural_history when generating structural animations"))
    if (show_force_vectors || show_moment_vectors) && isempty(result.aerodynamic_nodal_load_history)
        throw(ArgumentError("Enable save_aerodynamic_load_history to draw load vectors"))
    end
    backend in (:gr,:plotlyjs) || throw(ArgumentError(
        "AeroBeams v1 structural animations support :gr and :plotlyjs"))
    mkpath(dirname(output_path))

    # Build a display-only model with a surface at construction. The solved
    # UVLM structure and its typed aerodynamic/state containers stay untouched.
    originalBeam = problem.model.beams[1]
    _,span,chord,sparFraction = AeroBeams.geometrical_properties_Pazy()
    surface = AeroBeams.create_AeroSurface(airfoil=deepcopy(AeroBeams.NACA0018),
        c=chord,normSparPos=sparFraction,ppSolver=AeroBeams.QuasiSteady())
    beam = AeroBeams.create_Beam(name="Pazy UVLM visualization",length=span,
        nElements=originalBeam.discretization.nElements,
        normalizedNodalPositions=originalBeam.discretization.normalizedNodalPositions,
        S=originalBeam.sectional.S,I=originalBeam.sectional.I,
        p0=originalBeam.geometry.p0,surface=surface)
    bcs = AeroBeams.BC[]
    if show_clamped_root
        push!(bcs,AeroBeams.create_BC(beam=beam,node=1,
            conditions=Dict(type=>0.0 for type in (:u1A,:u2A,:u3A,:p1A,:p2A,:p3A))))
    end
    if show_force_vectors || show_moment_vectors
        components = Int[]
        show_force_vectors && push!(components,1)
        show_moment_vectors && push!(components,6)
        types = (:F1A,:F2A,:F3A,:M1A,:M2A,:M3A)
        for node in axes(result.aerodynamic_nodal_load_history,2)
            conditions = Dict(types[d] => let d=d,node=node
                t -> result.aerodynamic_nodal_load_history[d,node,force_history_index(result,t)]
            end for d in components)
            push!(bcs,AeroBeams.create_BC(beam=beam,node=node,conditions=conditions))
        end
    end
    model = AeroBeams.create_Model(beams=[beam],BCs=bcs)
    history = problem.solutionHistory
    savedCount = history.index - 1
    savedCount >= 2 || throw(ArgumentError("Structural animation requires two saved states"))
    frameIndices = unique(vcat(collect(1:result.animation_stride:savedCount),savedCount))
    displayProblem = AeroBeams.create_DynamicProblem(model=model,
        timeVector=history.t[frameIndices],trackingTimeSteps=true,displayProgress=false)
    # Rendering reads selected accepted states, including the final frame.
    # No structural/aerodynamic solve or mutation of the source history occurs.
    for field in (:t,:states,:hacStates,:aeroVariables,:bliAeroVariables,:bloAeroVariables)
        setfield!(displayProblem.solutionHistory,field,getfield(history,field)[frameIndices])
    end
    displayProblem.solutionHistory.index = length(frameIndices) + 1
    model.kinematicsBasisA = deepcopy(problem.model.kinematicsBasisA)
    timeIndices = 1 .+ (frameIndices .- 1) .* problem.trackingFrequency
    for field in (:u_AVec,:v_AVec,:ω_AVec,:vdot_AVec,:ωdot_AVec,:R_AVec)
        setfield!(model.kinematicsBasisA,field,
            getfield(problem.model.kinematicsBasisA,field)[timeIndices])
    end
    AeroBeams.precompute_BCs_history!(displayProblem)
    AeroBeams.plot_dynamic_deformation(displayProblem; backend=backend,follow=:basisA,
        view=camera,plotFrequency=1,plotUndeformed=false,
        plotBCs=show_clamped_root || show_force_vectors || show_moment_vectors,
        plotDistLoads=false,plotAeroSurf=true,surfα=Float64(surface_alpha),
        lw=elastic_axis_linewidth,colorDef=:blue,showTimeStamp=true,showScale=true,
        scale=deformation_scale,loadsSizeScaler=force_vector_scale,
        plotLimits=([-span/2,span/2],[-span/2,span/2],[0.0,span]),fps=fps,
        save=true,savePath=path_for_aerobeams(output_path),displayProgress=true)
    println("Structural animation written to $(abspath(output_path))")
    return abspath(output_path)
end

function wake_vertex_positions(wake::AbstractMatrix)
    wake_rows,spanwise_panels = size(wake)
    wake_rows > 0 || return Matrix{Vector{Float64}}(undef,0,0)
    positions = Matrix{typeof(wake[1,1].rtl)}(
        undef,wake_rows + 1,spanwise_panels + 1
    )
    for row = 1:wake_rows,span = 1:spanwise_panels
        positions[row,span] = wake[row,span].rtl
    end
    for row = 1:wake_rows
        positions[row,end] = wake[row,end].rtr
    end
    for span = 1:spanwise_panels
        positions[end,span] = wake[end,span].rbl
    end
    positions[end,end] = wake[end,end].rbr
    return positions
end

function lattice_lines(positions::AbstractMatrix)
    x,y,z = Float64[],Float64[],Float64[]
    function add_line(points)
        for point in points
            push!(x,point[1])
            push!(y,point[2])
            push!(z,point[3])
        end
        push!(x,NaN)
        push!(y,NaN)
        push!(z,NaN)
    end
    for row in axes(positions,1)
        add_line(positions[row,:])
    end
    for column in axes(positions,2)
        add_line(positions[:,column])
    end
    return x,y,z
end

function animation_limits(result; extra_padding=0.0)
    minimumPoint = fill(Inf,3)
    maximumPoint = fill(-Inf,3)
    function include!(positions)
        for point in positions,component = 1:3
            minimumPoint[component] = min(minimumPoint[component],point[component])
            maximumPoint[component] = max(maximumPoint[component],point[component])
        end
    end
    # Include every saved frame once, then reuse these limits throughout the
    # GIF. This prevents the axes and apparent camera position from moving as
    # the initially empty wake grows.
    for frame in eachindex(result.animation_time)
        surfacePositions = UVLM.imperial_nodal_positions(
            result.aerodynamic_surface_history[frame][1])
        include!(surfacePositions)
        wake = result.wake_history[frame][1]
        if size(wake,1) > 0
            include!(wake_vertex_positions(wake))
        end
    end

    # A fixed cube keeps the same limits and metres-per-unit scale on all three
    # axes for every animation frame.
    span = maximum(maximumPoint - minimumPoint)
    padding = max(0.04 * span,0.005) + extra_padding
    halfWidth = span / 2 + padding
    center = (minimumPoint + maximumPoint) / 2
    return ntuple(3) do component
        (center[component] - halfWidth,center[component] + halfWidth)
    end
end

"""Create a GIF of the deformed UVLM lattice and its free wake."""
function save_wake_animation(result;
    output_path::String=joinpath(DEFAULT_OUTPUT_DIRECTORY,"pazy_uvlm_wake.gif"),
    fps::Int=30,
    show_force_vectors::Bool=false,
    force_vector_scale::Real=1.0,
    camera=(45,30))

    isempty(result.animation_time) && throw(ArgumentError(
        "Enable save_animation_history when generating wake animations"))
    structuralTimes = result.structural_problem.solutionHistory.t
    if !isempty(structuralTimes)
        structural_frame_times = structuralTimes[1:result.animation_stride:end]
        length(structural_frame_times) == length(result.animation_time) ||
            error("Structural and wake animations contain different frame counts")
        all(isapprox.(structural_frame_times,result.animation_time)) ||
            error("Structural and wake animation times are not synchronized")
    end
    if show_force_vectors && (isempty(structuralTimes) || isempty(result.aerodynamic_nodal_load_history))
        throw(ArgumentError("Enable save_structural_history and save_aerodynamic_load_history to draw load vectors"))
    end
    mkpath(dirname(output_path))
    Plots.gr()
    animation = Plots.Animation()
    maximumForce = show_force_vectors ? maximum_aerodynamic_force(result) : 0.0
    _,span,_,_ = AeroBeams.geometrical_properties_Pazy()
    arrowMargin = show_force_vectors && maximumForce > 0 ?
                  force_vector_scale * span / 10 : 0.0
    limits = animation_limits(result; extra_padding=arrowMargin)
    for frame in eachindex(result.animation_time)
        # Plot in conventional aerodynamic axes. The UVLM data already use
        # x=downstream, y=span and z=up; converting back to AeroBeams axes here
        # would put the span on the plot's vertical axis and make the wake look
        # like a wall.
        surfacePositions = UVLM.imperial_nodal_positions(
            result.aerodynamic_surface_history[frame][1])
        wake = result.wake_history[frame][1]
        sx,sy,sz = lattice_lines(surfacePositions)
        roundedTime = round(result.animation_time[frame]; sigdigits=2)
        figure = Plots.plot3d(
            sx,sy,sz;
            color=:blue,
            linewidth=2,
            label=false,
            xlims=limits[1],
            ylims=limits[2],
            zlims=limits[3],
            xlabel="Downstream [m]",
            ylabel="Span [m]",
            zlabel="Up [m]",
            title="Time = $roundedTime s",
            titlefontsize=14,
            guidefontsize=10,
            tickfontsize=8,
            camera=camera,
            aspect_ratio=:equal,
            size=(900,700),
            dpi=120,
            grid=true,
            axis=true,
            legend=false,
            widen=false
        )
        if size(wake,1) > 0
            wakePositions = wake_vertex_positions(wake)
            wx,wy,wz = lattice_lines(wakePositions)
            Plots.plot3d!(
                figure,wx,wy,wz;
                color=:gray45,
                linewidth=0.8,
                label=false
            )
        end
        if show_force_vectors && maximumForce > 0
            time_index = force_history_index(result,result.animation_time[frame])
            # Beam positions and transferred forces are stored in AeroBeams A
            # axes. Rotate both into the same aerodynamic view axes as the wake.
            positions = A_TO_UVLM * beam_node_positions(result,time_index)
            forces = A_TO_UVLM *
                     result.aerodynamic_nodal_load_history[1:3,:,time_index]
            vectors = (force_vector_scale * span / (10 * maximumForce)) .* forces
            origins = positions .- vectors
            Plots.quiver!(
                figure,
                vec(origins[1,:]),vec(origins[2,:]),vec(origins[3,:]);
                quiver=(vec(vectors[1,:]),vec(vectors[2,:]),vec(vectors[3,:])),
                color=:green,
                linewidth=2,
                quiverhead=0.5,
                label=false
            )
        end
        Plots.frame(animation,figure)
    end
    Plots.gif(animation,output_path; fps=fps,show_msg=false)
    println("UVLM wake animation written to $(abspath(output_path))")
    return abspath(output_path)
end

function save_animations(result; output_directory::String=DEFAULT_OUTPUT_DIRECTORY,
    case_label::String="pazy_response",
    case_suffix::String="",
    structure_backend::Symbol=:gr,
    structure_camera=(45,45),structure_surface_alpha::Real=0.5,
    structure_axis_linewidth::Real=2.5,show_clamped_root::Bool=true,
    fps::Int=30,show_force_vectors::Bool=false,
    show_moment_vectors::Bool=false,
    force_vector_scale::Real=1.0,wake_camera=(45,30))

    suffix = isempty(case_suffix) ? "" : "_$(case_suffix)"
    structural = save_structural_animation(
        result;
        output_path=joinpath(output_directory,"$(case_label)_structure$(suffix).gif"),
        fps=fps,
        backend=structure_backend,
        camera=structure_camera,
        surface_alpha=structure_surface_alpha,
        elastic_axis_linewidth=structure_axis_linewidth,
        show_clamped_root=show_clamped_root,
        show_force_vectors=show_force_vectors,
        show_moment_vectors=show_moment_vectors,
        force_vector_scale=force_vector_scale
    )
    wake = save_wake_animation(
        result;
        output_path=joinpath(output_directory,"$(case_label)_uvlm_wake$(suffix).gif"),
        fps=fps,
        show_force_vectors=show_force_vectors,
        force_vector_scale=force_vector_scale,
        camera=wake_camera
    )
    return (; structural,wake)
end

"""Plot and save the wingtip bending-displacement and twist time histories."""
function save_tip_time_histories(result;
    output_path::String=joinpath(DEFAULT_OUTPUT_DIRECTORY,"pazy_tip_time_histories.png"))

    mkpath(dirname(output_path))
    bendingPlot = Plots.plot(
        result.time,
        result.tip_bending_displacement;
        color=:black,
        linewidth=4,
        xlabel="Time (s)",
        ylabel="Tip displacement (m)",
        tickfontsize=20,
        labelfontsize=24,
        title="",
        label=false,
        grid=false,
        left_margin=22Plots.mm,
        bottom_margin=12Plots.mm,
        top_margin=5Plots.mm,
        right_margin=5Plots.mm
    )
    twistPlot = Plots.plot(
        result.time,
        result.tip_twist_degrees;
        color=:black,
        linewidth=4,
        xlabel="Time (s)",
        ylabel="Tip twist (deg)",
        tickfontsize=20,
        labelfontsize=24,
        title="",
        label=false,
        grid=false,
        left_margin=22Plots.mm,
        bottom_margin=12Plots.mm,
        top_margin=5Plots.mm,
        right_margin=5Plots.mm
    )

    figure = Plots.plot(bendingPlot,twistPlot; layout=(2,1),size=(1000,900))
    Plots.savefig(figure,output_path)
    println("Wingtip time histories written to $(abspath(output_path))")
    return abspath(output_path)
end
