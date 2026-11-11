# ==============================================================================
# Aerodynamic grid utilities for the Bohnisch time-domain driver.
# ==============================================================================

using StaticArrays
using FLOWMath

"""
    linear_interpolate_1d(x_reference, y_reference, x_target)

Linearly interpolate a one-dimensional data set and hold the first or last
value outside the reference interval. Reference coordinates must be strictly
increasing.
"""
function linear_interpolate_1d(x_reference::AbstractVector{<:Real},
    y_reference::AbstractVector{<:Real},x_target::AbstractVector{<:Real})

    length(x_reference) == length(y_reference) ||
        throw(DimensionMismatch("Reference coordinates and values must have equal lengths"))
    isempty(x_reference) && throw(ArgumentError("Reference data must not be empty"))
    all(diff(x_reference) .> 0) ||
        throw(ArgumentError("Reference coordinates must be strictly increasing"))

    outputType = promote_type(Float64,eltype(x_reference),eltype(y_reference),eltype(x_target))
    interpolated = Vector{outputType}(undef,length(x_target))
    for (targetIndex,targetCoordinate) in pairs(x_target)
        if targetCoordinate <= first(x_reference)
            interpolated[targetIndex] = first(y_reference)
        elseif targetCoordinate >= last(x_reference)
            interpolated[targetIndex] = last(y_reference)
        else
            upperIndex = searchsortedfirst(x_reference,targetCoordinate)
            lowerIndex = upperIndex - 1
            fraction =
                (targetCoordinate - x_reference[lowerIndex]) /
                (x_reference[upperIndex] - x_reference[lowerIndex])
            interpolated[targetIndex] = (1 - fraction) * y_reference[lowerIndex] +
                                        fraction * y_reference[upperIndex]
        end
    end
    return interpolated
end

function span_position_to_node_index(span_nodes::AbstractVector{<:Real},y_position::Real)
    isempty(span_nodes) && error("span_nodes must not be empty")
    yMin,yMax = extrema(span_nodes)
    yMin <= y_position <= yMax || error("Propeller span position $y_position is outside the wing span [$yMin, $yMax]")

    return argmin(abs.(span_nodes .- y_position))
end

function span_positions_to_node_indices(propeller_span_positions::AbstractVector{<:Real},
    span_nodes::AbstractVector{<:Real})

    return [span_position_to_node_index(span_nodes,y_position) for y_position in propeller_span_positions]
end

function propeller_attachment_nodes_from_eta(propeller_eta::AbstractVector{<:Real},
    span_length::Real,span_nodes::AbstractVector{<:Real})

    propeller_span_positions = collect(propeller_eta) .* span_length
    prop_attach_nodes = span_positions_to_node_indices(propeller_span_positions,span_nodes)
    return propeller_span_positions,prop_attach_nodes
end

function generate_panel_grid_and_interpolate(span_length::Float64,chordDistribution::Vector{Float64},
    xle_distribution::Vector{Float64},ns::Int,nc::Int,
    u_x::Vector{Float64},u_y::Vector{Float64},u_z::Vector{Float64},
    theta_x::Vector{Float64},theta_y::Vector{Float64},theta_z::Vector{Float64};
    elastic_axis_fraction::Real=0.35)

    0.0 <= elastic_axis_fraction <= 1.0 ||
        throw(ArgumentError("elastic_axis_fraction must be between 0 and 1"))

    panelGrid = zeros(Float64,3,nc + 1,ns + 1)
    spanPositions = range(0,span_length,length=ns + 1)

    for (i,y) in enumerate(spanPositions)
        local_chord = chordDistribution[i]
        localXle = xle_distribution[i]
        chordPositions = range(0,local_chord,length=nc + 1)

        eaLocal = local_chord * elastic_axis_fraction
        ea_global_x = localXle + eaLocal + u_x[i]
        r_EA = SVector(ea_global_x,y + u_y[i],u_z[i])

        R = RotationMatrix(theta_z[i],3) * RotationMatrix(theta_x[i],1) * RotationMatrix(theta_y[i],2)

        for (j,xLocal) in enumerate(chordPositions)
            xRelative = xLocal - eaLocal
            pLocal = SVector(xRelative,0.0,0.0)
            panelGrid[:,j,i] = r_EA + R * pLocal
        end
    end
    return panelGrid
end

function generate_aero_panel_grid_and_interpolate(span_length::Float64,chordDistribution::Vector{Float64},
    xle_distribution::Vector{Float64},ns::Int,nc::Int,
    u_x::Vector{Float64},u_y::Vector{Float64},u_z::Vector{Float64},
    theta_x::Vector{Float64},theta_y::Vector{Float64},theta_z::Vector{Float64};
    elastic_axis_fraction::Real=0.35)

    0.0 <= elastic_axis_fraction <= 1.0 ||
        throw(ArgumentError("elastic_axis_fraction must be between 0 and 1"))

    panelGrid = zeros(Float64,3,nc + 1,ns + 1)
    spanPositions = range(0,span_length,length=ns + 1)

    for (i,y) in enumerate(spanPositions)
        local_chord = chordDistribution[i]
        localXle = xle_distribution[i]
        chordPositions = range(0.25 * local_chord / nc,local_chord + 0.25 * local_chord / nc,length=nc + 1)

        eaLocal = local_chord * elastic_axis_fraction
        ea_global_x = localXle + eaLocal + u_x[i]
        r_EA = SVector(ea_global_x,y + u_y[i],u_z[i])

        R = RotationMatrix(theta_z[i],3) * RotationMatrix(theta_x[i],1) * RotationMatrix(theta_y[i],2)

        for (j,xLocal) in enumerate(chordPositions)
            xRelative = xLocal - eaLocal
            pLocal = SVector(xRelative,0.0,0.0)
            panelGrid[:,j,i] = r_EA + R * pLocal
        end
    end
    return panelGrid
end

function generate_propeller_blades_grid(totalRadius::Float64,chordLength::Float64,
    ns::Int,nc::Int,blade_twists::Vector{Float64},numBlades::Int)

    length(blade_twists) == ns + 1 || error("blade_twists must have length ns+1")
    totalRadius > 0 || error("Total radius must be positive")

    baseGrid = zeros(3,nc + 1,ns + 1)
    radialPositions = range(0.0,totalRadius,length=ns + 1)
    chordPositions = range(0.0,chordLength,length=nc + 1)

    for (i,r) in enumerate(radialPositions)
        theta = blade_twists[i]
        for (j,xChord) in enumerate(chordPositions)
            x = xChord * cos(theta)
            z = -xChord * sin(theta)
            baseGrid[:,j,i] = [x,r,z]
        end
    end

    bladeGrids = [similar(baseGrid) for _ = 1:numBlades]
    angularStep = 2pi / numBlades
    for (bladeIdx,grid) in enumerate(bladeGrids)
        phi = angularStep * (bladeIdx - 1)
        cphi,sphi = cos(phi),sin(phi)
        for i = 1:size(baseGrid,3)
            for j = 1:size(baseGrid,2)
                x,y,z = baseGrid[:,j,i]
                grid[:,j,i] = [x,y * cphi - z * sphi,y * sphi + z * cphi]
            end
        end
    end
    return bladeGrids
end

