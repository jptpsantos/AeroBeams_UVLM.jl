# Included after PazyWingUVLMCoupling.jl; uses its AeroBeams and UVLM imports.
const Plots = AeroBeams.Plots

const DEFAULT_OUTPUT_DIRECTORY = joinpath(@__DIR__, "output")

function path_for_aerobeams(absolute_path::String)
    root = abspath(REPOSITORY_ROOT)
    path = abspath(absolute_path)
    startswith(lowercase(path), lowercase(root)) ||
        throw(ArgumentError("AeroBeams animations must be saved inside the repository"))
    return "/" * replace(relpath(path, root), '\\' => '/')
end

function maximum_aerodynamic_force(result)
    history = result.aerodynamic_nodal_load_history
    maximum_force = 0.0
    for step in axes(history, 3), node in axes(history, 2)
        force = history[1:3, node, step]
        maximum_force = max(maximum_force, sqrt(sum(abs2, force)))
    end
    return maximum_force
end

function force_history_index(result, time)
    return clamp(round(Int, time / result.dt) + 1, 1, length(result.time))
end

function beam_node_positions(result, time_index)
    problem = result.structural_problem
    elements = problem.model.elements
    states = problem.nodalStatesOverTime[time_index]
    positions = zeros(3, length(elements) + 1)
    positions[:, 1] = elements[1].r_n1
    for (element_index, element) in enumerate(elements)
        positions[:, element_index + 1] =
            element.r_n2 + states[element.globalID].u_n2
    end
    return positions
end

"""Create the usual AeroBeams nonlinear structural-deformation GIF."""
function save_structural_animation(result;
    output_path::String = joinpath(DEFAULT_OUTPUT_DIRECTORY, "pazy_structure.gif"),
    fps::Int = 30,
    deformation_scale::Real = 1.0,
    show_force_vectors::Bool = true,
    force_vector_scale::Real = 1.0)

    mkpath(dirname(output_path))
    stride = result.animation_stride
    problem = result.structural_problem

    # The structural model has no AeroBeams aerodynamic loads, because UVLM
    # supplies them. Temporarily attach only the standard Pazy NACA 0018
    # geometry so AeroBeams can draw its usual airfoil surfaces in the GIF.
    _, span, chord, spar_fraction = AeroBeams.geometrical_properties_Pazy()
    beam = problem.model.beams[1]
    original_surface = beam.aeroSurface
    original_element_aero = [element.aero for element in problem.model.elements]
    force_bcs = problem.model.BCs[2:end]
    original_force_bcs = [(
        types=bc.types, values=bc.values, toBeTrimmed=bc.toBeTrimmed,
        Fmax=bc.Fmax, Mmax=bc.Mmax,
    ) for bc in force_bcs]
    display_surface = AeroBeams.create_AeroSurface(
        airfoil=deepcopy(AeroBeams.NACA0018), c=chord, normSparPos=spar_fraction)
    try
        beam.aeroSurface = display_surface
        for element in problem.model.elements
            element.aero = AeroBeams.AeroProperties(
                display_surface, beam.rotationParametrization, element.R0,
                element.x1, element.x1_norm,
                element.x1_n1_norm, element.x1_n2_norm)
        end
        if show_force_vectors
            maximum_force = max(maximum_aerodynamic_force(result), eps(Float64))
            for (node, bc) in enumerate(force_bcs)
                bc.types = ["F1A", "F2A", "F3A"]
                bc.values = [let component=component, node=node
                    time -> result.aerodynamic_nodal_load_history[
                        component, node, force_history_index(result, time)]
                end for component in 1:3]
                bc.toBeTrimmed = falses(3)
                # One reference value makes arrow lengths comparable between
                # every node and every frame.
                bc.Fmax = maximum_force
            end
        end
        # Use the earlier oblique view. The default AeroBeams camera sees this
        # particular wing nearly edge-on, hiding most of its deformation.
        AeroBeams.plot_dynamic_deformation(
            problem;
            refBasis = "A",
            plotFrequency = stride,
            plotUndeformed = true,
            plotBCs = show_force_vectors,
            plotDistLoads = false,
            plotAeroSurf = true,
            surfα = 0.55,
            view = (35, 20),
            scale = deformation_scale,
            loadsSizeScaler = force_vector_scale,
            # Equal 0.64 m ranges avoid geometric distortion between axes.
            plotLimits = ([-0.32, 0.32], [-0.32, 0.32], [-0.04, 0.60]),
            fps = fps,
            save = true,
            savePath = path_for_aerobeams(output_path),
            displayProgress = true,
        )
    finally
        beam.aeroSurface = original_surface
        for (element, original_aero) in zip(problem.model.elements, original_element_aero)
            element.aero = original_aero
        end
        for (bc, original) in zip(force_bcs, original_force_bcs)
            bc.types = original.types
            bc.values = original.values
            bc.toBeTrimmed = original.toBeTrimmed
            bc.Fmax = original.Fmax
            bc.Mmax = original.Mmax
            AeroBeams.update_BC_data!(bc, last(result.time))
        end
    end
    println("Structural animation written to $(abspath(output_path))")
    return abspath(output_path)
end

function wake_vertex_positions(wake::AbstractMatrix)
    wake_rows, spanwise_panels = size(wake)
    wake_rows > 0 || return Matrix{Vector{Float64}}(undef, 0, 0)
    positions = Matrix{typeof(wake[1, 1].rtl)}(
        undef, wake_rows + 1, spanwise_panels + 1,
    )
    for row in 1:wake_rows, span in 1:spanwise_panels
        positions[row, span] = wake[row, span].rtl
    end
    for row in 1:wake_rows
        positions[row, end] = wake[row, end].rtr
    end
    for span in 1:spanwise_panels
        positions[end, span] = wake[end, span].rbl
    end
    positions[end, end] = wake[end, end].rbr
    return positions
end

function lattice_lines(positions::AbstractMatrix)
    x, y, z = Float64[], Float64[], Float64[]
    function add_line(points)
        for point in points
            push!(x, point[1]); push!(y, point[2]); push!(z, point[3])
        end
        push!(x, NaN); push!(y, NaN); push!(z, NaN)
    end
    for row in axes(positions, 1)
        add_line(positions[row, :])
    end
    for column in axes(positions, 2)
        add_line(positions[:, column])
    end
    return x, y, z
end

function animation_limits(result, frame; extra_padding=0.0)
    minimum_point = fill(Inf, 3)
    maximum_point = fill(-Inf, 3)
    function include!(positions)
        for point in positions, component in 1:3
            minimum_point[component] = min(minimum_point[component], point[component])
            maximum_point[component] = max(maximum_point[component], point[component])
        end
    end
    surface_positions = UVLM.imperial_nodal_positions(
        result.aerodynamic_surface_history[frame][1])
    include!(surface_positions)
    wake = result.wake_history[frame][1]
    if size(wake, 1) > 0
        include!(wake_vertex_positions(wake))
    end

    # Fit the cube to the wake that exists in this frame. This shows the wing
    # clearly while the wake starts shedding, then zooms out as the wake grows.
    # A cube keeps the same metres-per-unit scale on all three axes.
    span = maximum(maximum_point - minimum_point)
    padding = max(0.04 * span, 0.005) + extra_padding
    half_width = span / 2 + padding
    center = (minimum_point + maximum_point) / 2
    return ntuple(3) do component
        (center[component] - half_width, center[component] + half_width)
    end
end

"""Create a GIF of the deformed UVLM lattice and its free wake."""
function save_wake_animation(result;
    output_path::String = joinpath(DEFAULT_OUTPUT_DIRECTORY, "pazy_uvlm_wake.gif"),
    fps::Int = 30,
    show_force_vectors::Bool = true,
    force_vector_scale::Real = 1.0,
    camera = (45, 30))

    mkpath(dirname(output_path))
    Plots.gr()
    animation = Plots.Animation()
    structural_frame_times = result.structural_problem.savedTimeVector[
        1:result.animation_stride:end]
    length(structural_frame_times) == length(result.animation_time) ||
        error("Structural and wake animations contain different frame counts")
    all(isapprox.(structural_frame_times, result.animation_time)) ||
        error("Structural and wake animation times are not synchronized")
    maximum_force = maximum_aerodynamic_force(result)
    _, span, _, _ = AeroBeams.geometrical_properties_Pazy()
    arrow_margin = show_force_vectors && maximum_force > 0 ?
        force_vector_scale * span / 10 : 0.0
    for frame in eachindex(result.animation_time)
        limits = animation_limits(result, frame; extra_padding=arrow_margin)
        # Plot in conventional aerodynamic axes. The UVLM data already use
        # x=downstream, y=span and z=up; converting back to AeroBeams axes here
        # would put the span on the plot's vertical axis and make the wake look
        # like a wall.
        surface_positions = UVLM.imperial_nodal_positions(
            result.aerodynamic_surface_history[frame][1])
        wake = result.wake_history[frame][1]
        sx, sy, sz = lattice_lines(surface_positions)
        rounded_time = round(result.animation_time[frame]; sigdigits=2)
        figure = Plots.plot3d(
            sx, sy, sz;
            color = :blue,
            linewidth = 2,
            label = false,
            xlims = limits[1],
            ylims = limits[2],
            zlims = limits[3],
            xlabel = "Downstream [m]",
            ylabel = "Span [m]",
            zlabel = "Up [m]",
            title = "Time = $rounded_time s",
            titlefontsize = 14,
            guidefontsize = 10,
            tickfontsize = 8,
            camera = camera,
            aspect_ratio = :equal,
            size = (900, 700),
            dpi = 120,
            grid = true,
            axis = true,
            legend = false,
            widen = false,
        )
        if size(wake, 1) > 0
            wake_positions = wake_vertex_positions(wake)
            wx, wy, wz = lattice_lines(wake_positions)
            Plots.plot3d!(
                figure, wx, wy, wz;
                color = :gray45,
                linewidth = 0.8,
                label = false,
            )
        end
        if show_force_vectors && maximum_force > 0
            time_index = force_history_index(result, result.animation_time[frame])
            # Beam positions and transferred forces are stored in AeroBeams A
            # axes. Rotate both into the same aerodynamic view axes as the wake.
            positions = A_TO_UVLM * beam_node_positions(result, time_index)
            forces = A_TO_UVLM *
                result.aerodynamic_nodal_load_history[1:3, :, time_index]
            vectors = (force_vector_scale * span / (10 * maximum_force)) .* forces
            origins = positions .- vectors
            Plots.quiver!(
                figure,
                vec(origins[1, :]), vec(origins[2, :]), vec(origins[3, :]);
                quiver = (vec(vectors[1, :]), vec(vectors[2, :]), vec(vectors[3, :])),
                color = :green,
                linewidth = 2,
                quiverhead = 0.5,
                label = false,
            )
        end
        Plots.frame(animation, figure)
    end
    Plots.gif(animation, output_path; fps = fps, show_msg = false)
    println("UVLM wake animation written to $(abspath(output_path))")
    return abspath(output_path)
end

function save_animations(result; output_directory::String = DEFAULT_OUTPUT_DIRECTORY,
    fps::Int = 30, show_force_vectors::Bool = true,
    force_vector_scale::Real = 1.0, wake_camera = (45, 30))

    structural = save_structural_animation(
        result;
        output_path = joinpath(output_directory, "pazy_structure.gif"),
        fps = fps,
        show_force_vectors = show_force_vectors,
        force_vector_scale = force_vector_scale,
    )
    wake = save_wake_animation(
        result;
        output_path = joinpath(output_directory, "pazy_uvlm_wake.gif"),
        fps = fps,
        show_force_vectors = show_force_vectors,
        force_vector_scale = force_vector_scale,
        camera = wake_camera,
    )
    return (; structural, wake)
end

"""Plot and save the wingtip bending-displacement and twist time histories."""
function save_tip_time_histories(result;
    output_path::String = joinpath(DEFAULT_OUTPUT_DIRECTORY, "pazy_tip_time_histories.png"))

    mkpath(dirname(output_path))
    bending_plot = Plots.plot(
        result.time,
        result.tip_bending_displacement;
        color = :blue,
        linewidth = 2,
        xlabel = "Time [s]",
        ylabel = "Bending displacement [m]",
        title = "Pazy wingtip bending",
        label = false,
        grid = true,
    )
    twist_plot = Plots.plot(
        result.time,
        result.tip_twist_degrees;
        color = :red,
        linewidth = 2,
        xlabel = "Time [s]",
        ylabel = "Twist [deg]",
        title = "Pazy wingtip twist",
        label = false,
        grid = true,
    )
    figure = Plots.plot(bending_plot, twist_plot; layout = (2, 1), size = (900, 700))
    Plots.savefig(figure, output_path)
    println("Wingtip time histories written to $(abspath(output_path))")
    return abspath(output_path)
end
