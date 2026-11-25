# Wake-history capture and rendering helpers for the Chang example.

"""
    chang_wake_vertex_positions(wake)

Reconstruct the `(nwake + 1, nspan + 1)` wake-vertex lattice from a matrix of
`WakePanel`s. The wake matrix passed to this function must contain active rows
only, ordered from the newest row at the lifting-surface trailing edge to the
oldest downstream row.
"""
function chang_wake_vertex_positions(wake::AbstractMatrix)
    wake_rows,spanwise_panels = size(wake)
    wake_rows > 0 || throw(ArgumentError("wake must contain at least one active row"))
    spanwise_panels > 0 || throw(ArgumentError("wake must contain spanwise panels"))

    pointType = typeof(wake[1,1].rtl)
    positions = Matrix{pointType}(undef,wake_rows + 1,spanwise_panels + 1)

    for row = 1:wake_rows,span = 1:spanwise_panels
        positions[row,span] = wake[row,span].rtl
    end
    for row = 1:wake_rows
        positions[row,spanwise_panels+1] = wake[row,spanwise_panels].rtr
    end
    for span = 1:spanwise_panels
        positions[wake_rows+1,span] = wake[wake_rows,span].rbl
    end
    positions[wake_rows+1,spanwise_panels+1] =
        wake[wake_rows,spanwise_panels].rbr
    return positions
end

"""Return NaN-separated x/y/z arrays for plotting all lines of a lattice."""
function chang_lattice_lines(positions::AbstractMatrix)
    x = Float64[]
    y = Float64[]
    z = Float64[]

    function append_point!(point)
        push!(x,point[1])
        push!(y,point[2])
        push!(z,point[3])
    end
    function end_polyline!()
        push!(x,NaN)
        push!(y,NaN)
        push!(z,NaN)
    end

    for row in axes(positions,1)
        for column in axes(positions,2)
            append_point!(positions[row,column])
        end
        end_polyline!()
    end
    for column in axes(positions,2)
        for row in axes(positions,1)
            append_point!(positions[row,column])
        end
        end_polyline!()
    end
    return x,y,z
end

function chang_animation_limits(surface_history,wake_history;
    padding_fraction::Real=0.05)

    0 <= padding_fraction < 1 || throw(ArgumentError(
        "padding_fraction must lie in [0, 1)"
    ))
    minimumCoordinates = fill(Inf,3)
    maximumCoordinates = fill(-Inf,3)

    function include_positions!(positions)
        for point in positions
            for component = 1:3
                minimumCoordinates[component] = min(
                    minimumCoordinates[component],
                    point[component]
                )
                maximumCoordinates[component] = max(
                    maximumCoordinates[component],
                    point[component]
                )
            end
        end
    end

    for frameSurfaces in surface_history
        for surface in frameSurfaces
            include_positions!(imperial_nodal_positions(surface))
        end
    end
    for frameWakes in wake_history
        for wake in frameWakes
            size(wake,1) == 0 && continue
            include_positions!(chang_wake_vertex_positions(wake))
        end
    end
    all(isfinite,minimumCoordinates) || error("No geometry was supplied")

    # Use the same numerical range on all three axes. Plots otherwise expands
    # each axis independently to fill the plotting box, which makes a 2.3 m
    # propeller disk appear comparable to the 7.5 m modeled wing span.
    coordinateSpans = maximumCoordinates - minimumCoordinates
    commonSpan = maximum(coordinateSpans)
    commonSpan > 0 || (commonSpan = 1.0)
    paddedHalfSpan = 0.5 * commonSpan * (1 + 2 * padding_fraction)

    return ntuple(3) do component
        center = 0.5 * (
            minimumCoordinates[component] + maximumCoordinates[component]
        )
        (center - paddedHalfSpan,center + paddedHalfSpan)
    end
end

"""Return every requested tick inside a closed axis interval."""
function chang_animation_ticks(limits,spacing::Real)
    spacing > 0 || throw(ArgumentError("tick spacing must be positive"))
    firstTick = ceil(limits[1] / spacing) * spacing
    lastTick = floor(limits[2] / spacing) * spacing
    firstTick <= lastTick || throw(ArgumentError(
        "tick spacing does not place a tick inside the axis limits"
    ))
    return collect(firstTick:spacing:lastTick)
end

"""
    record_chang_animation_frame!(
        surface_history,
        wake_history,
        active_wake_rows_history,
        time_history,
        system,
        active_wake_rows,
        time,
    )

Store one accepted UVLM state for later animation. Only active wake rows are
copied, which avoids reading uninitialized preallocated wake panels and keeps
the animation history substantially smaller than a full solver snapshot.
"""
function record_chang_animation_frame!(
    surface_history,
    wake_history,
    active_wake_rows_history,
    time_history,
    system,
    active_wake_rows,
    time::Real
)
    length(active_wake_rows) == length(system.wakes) || throw(DimensionMismatch(
        "One active wake-row count is required per aerodynamic surface"
    ))
    for surface_index in eachindex(system.wakes)
        0 <= active_wake_rows[surface_index] <= size(system.wakes[surface_index],1) ||
            throw(ArgumentError("Active wake-row count is outside allocated storage"))
    end

    push!(surface_history,[copy(surface) for surface in system.surfaces])
    push!(
        wake_history,
        [
            copy(view(system.wakes[surface_index],
                1:active_wake_rows[surface_index],:))
            for surface_index in eachindex(system.wakes)
        ]
    )
    push!(active_wake_rows_history,copy(active_wake_rows))
    push!(time_history,Float64(time))
    return nothing
end

"""
    animate_chang_wing_wake(
        surface_history,
        wake_history,
        active_wake_rows_history,
        time_history;
        output_path,
        fps=15,
        camera=(35, 25),
        axis_limits=nothing,
        tick_spacing=1,
    )

Create a GIF showing the accepted deformed wing and propeller lattices together
with their active free wakes. Every entry in `wake_history` must already be
cropped to its active row count; uninitialized allocated wake storage must not
be passed to this function.

The animation is intended for coupled-solver diagnostics: its frames must be
recorded only after a physical time step converges. Intermediate partitioned
coupling trials should never be included.
"""
function animate_chang_wing_wake(
    surface_history,
    wake_history,
    active_wake_rows_history,
    time_history;
    output_path::AbstractString,
    fps::Integer=15,
    camera=(35,25),
    figure_size=(900,900),
    axis_limits=nothing,
    tick_spacing::Real=1.0
)
    nFrames = length(time_history)
    nFrames > 0 || throw(ArgumentError("At least one frame is required"))
    length(surface_history) == nFrames || throw(DimensionMismatch(
        "surface_history and time_history must have the same length"
    ))
    length(wake_history) == nFrames || throw(DimensionMismatch(
        "wake_history and time_history must have the same length"
    ))
    length(active_wake_rows_history) == nFrames || throw(DimensionMismatch(
        "active_wake_rows_history and time_history must have the same length"
    ))
    fps > 0 || throw(ArgumentError("fps must be positive"))
    endswith(lowercase(output_path),".gif") || throw(ArgumentError(
        "output_path must have a .gif extension"
    ))

    if isnothing(axis_limits)
        xlimits,ylimits,zlimits = chang_animation_limits(
            surface_history,
            wake_history
        )
    else
        length(axis_limits) == 3 || throw(DimensionMismatch(
            "axis_limits must contain X, Y, and Z intervals"
        ))
        for limits in axis_limits
            length(limits) == 2 || throw(DimensionMismatch(
                "Each axis interval must contain a lower and upper limit"
            ))
            all(isfinite,limits) || throw(ArgumentError(
                "Axis limits must be finite"
            ))
            limits[1] < limits[2] || throw(ArgumentError(
                "Each lower axis limit must be smaller than its upper limit"
            ))
        end
        xlimits,ylimits,zlimits = axis_limits
    end
    xtickValues = chang_animation_ticks(xlimits,tick_spacing)
    ytickValues = chang_animation_ticks(ylimits,tick_spacing)
    ztickValues = chang_animation_ticks(zlimits,tick_spacing)
    animation = Plots.Animation()

    for frameIndex = 1:nFrames
        frameSurfaces = surface_history[frameIndex]
        frameWakes = wake_history[frameIndex]
        active_rows = active_wake_rows_history[frameIndex]
        length(frameSurfaces) == length(frameWakes) || throw(DimensionMismatch(
            "Each animation frame requires one wake matrix per surface"
        ))
        length(active_rows) == length(frameWakes) || throw(DimensionMismatch(
            "Each animation frame requires one wake-row count per surface"
        ))

        figure = Plots.plot3d(
            xlabel="X (m)",
            ylabel="Y (m)",
            zlabel="Z (m)",
            xlims=xlimits,
            ylims=ylimits,
            zlims=zlimits,
            xticks=xtickValues,
            yticks=ytickValues,
            zticks=ztickValues,
            aspect_ratio=:equal,
            camera=camera,
            size=figure_size,
            # Keep enough canvas around the 3-D box for the first and last tick
            # labels on every axis; zero margins can clip those endpoint labels.
            margin=6Plots.mm,
            legend=false,
            grid=true,
            foreground_color_grid=:gray80,
            background_color=:white,
            guidefontsize=10,
            tickfontsize=8
        )

        for surface_index in eachindex(frameSurfaces)
            surfacePositions = imperial_nodal_positions(frameSurfaces[surface_index])
            x,y,z = chang_lattice_lines(surfacePositions)
            Plots.plot3d!(
                figure,
                x,
                y,
                z;
                color=:black,
                linewidth=surface_index == 1 ? 2.2 : 1.4,
                alpha=0.95
            )

            wake = frameWakes[surface_index]
            size(wake,1) == 0 && continue
            size(wake,1) == active_rows[surface_index] || throw(DimensionMismatch(
                "Stored wake rows do not match active_wake_rows_history"
            ))
            wakePositions = chang_wake_vertex_positions(wake)
            xw,yw,zw = chang_lattice_lines(wakePositions)
            Plots.plot3d!(
                figure,
                xw,
                yw,
                zw;
                color=surface_index == 1 ? :darkblue : :darkred,
                linewidth=surface_index == 1 ? 1.1 : 0.8,
                alpha=0.62
            )
        end
        Plots.frame(animation,figure)
    end

    mkpath(dirname(output_path))
    Plots.gif(animation,output_path; fps=fps,show_msg=false)
    println("Wing-response/wake animation written to $output_path")
    return output_path
end
