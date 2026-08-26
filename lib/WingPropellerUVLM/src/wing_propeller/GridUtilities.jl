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
    y_reference::AbstractVector{<:Real}, x_target::AbstractVector{<:Real})

    length(x_reference) == length(y_reference) ||
        throw(DimensionMismatch("Reference coordinates and values must have equal lengths"))
    isempty(x_reference) && throw(ArgumentError("Reference data must not be empty"))
    all(diff(x_reference) .> 0) ||
        throw(ArgumentError("Reference coordinates must be strictly increasing"))

    output_type = promote_type(Float64, eltype(x_reference), eltype(y_reference), eltype(x_target))
    interpolated = Vector{output_type}(undef, length(x_target))
    for (target_index, target_coordinate) in pairs(x_target)
        if target_coordinate <= first(x_reference)
            interpolated[target_index] = first(y_reference)
        elseif target_coordinate >= last(x_reference)
            interpolated[target_index] = last(y_reference)
        else
            upper_index = searchsortedfirst(x_reference, target_coordinate)
            lower_index = upper_index - 1
            fraction = (target_coordinate - x_reference[lower_index]) /
                (x_reference[upper_index] - x_reference[lower_index])
            interpolated[target_index] = (1 - fraction) * y_reference[lower_index] +
                fraction * y_reference[upper_index]
        end
    end
    return interpolated
end

function span_position_to_node_index(span_nodes::AbstractVector{<:Real}, y_position::Real)
    isempty(span_nodes) && error("span_nodes must not be empty")
    y_min, y_max = extrema(span_nodes)
    y_min <= y_position <= y_max || error("Propeller span position $y_position is outside the wing span [$y_min, $y_max]")

    return argmin(abs.(span_nodes .- y_position))
end

function span_positions_to_node_indices(propeller_span_positions::AbstractVector{<:Real},
    span_nodes::AbstractVector{<:Real})

    return [span_position_to_node_index(span_nodes, y_position) for y_position in propeller_span_positions]
end

function propeller_attachment_nodes_from_eta(propeller_eta::AbstractVector{<:Real},
    span_length::Real, span_nodes::AbstractVector{<:Real})

    propeller_span_positions = collect(propeller_eta) .* span_length
    prop_attach_nodes = span_positions_to_node_indices(propeller_span_positions, span_nodes)
    return propeller_span_positions, prop_attach_nodes
end

function generate_panel_grid_and_interpolate(span_length::Float64, chord_distribution::Vector{Float64},
    xle_distribution::Vector{Float64}, ns::Int, nc::Int,
    u_x::Vector{Float64}, u_y::Vector{Float64}, u_z::Vector{Float64},
    theta_x::Vector{Float64}, theta_y::Vector{Float64}, theta_z::Vector{Float64};
    elastic_axis_fraction::Real = 0.35)

    0.0 <= elastic_axis_fraction <= 1.0 ||
        throw(ArgumentError("elastic_axis_fraction must be between 0 and 1"))

    panel_grid = zeros(Float64, 3, nc+1, ns+1)
    span_positions = range(0, span_length, length=ns+1)

    for (i, y) in enumerate(span_positions)
        local_chord = chord_distribution[i]
        local_xle = xle_distribution[i]
        chord_positions = range(0, local_chord, length=nc+1)

        ea_local = local_chord * elastic_axis_fraction
        ea_global_x = local_xle + ea_local + u_x[i]
        r_EA = SVector(ea_global_x, y + u_y[i], u_z[i])

        R = RotationMatrix(theta_z[i], 3) * RotationMatrix(theta_x[i], 1) * RotationMatrix(theta_y[i], 2)

        for (j, x_local) in enumerate(chord_positions)
            x_relative = x_local - ea_local
            p_local = SVector(x_relative, 0.0, 0.0)
            panel_grid[:, j, i] = r_EA + R * p_local
        end
    end
    return panel_grid
end

function generate_aero_panel_grid_and_interpolate(span_length::Float64, chord_distribution::Vector{Float64},
    xle_distribution::Vector{Float64}, ns::Int, nc::Int,
    u_x::Vector{Float64}, u_y::Vector{Float64}, u_z::Vector{Float64},
    theta_x::Vector{Float64}, theta_y::Vector{Float64}, theta_z::Vector{Float64};
    elastic_axis_fraction::Real = 0.35)

    0.0 <= elastic_axis_fraction <= 1.0 ||
        throw(ArgumentError("elastic_axis_fraction must be between 0 and 1"))

    panel_grid = zeros(Float64, 3, nc+1, ns+1)
    span_positions = range(0, span_length, length=ns+1)

    for (i, y) in enumerate(span_positions)
        local_chord = chord_distribution[i]
        local_xle = xle_distribution[i]
        chord_positions = range(0.25*local_chord/nc, local_chord + 0.25*local_chord/nc, length=nc+1)

        ea_local = local_chord * elastic_axis_fraction
        ea_global_x = local_xle + ea_local + u_x[i]
        r_EA = SVector(ea_global_x, y + u_y[i], u_z[i])

        R = RotationMatrix(theta_z[i], 3) * RotationMatrix(theta_x[i], 1) * RotationMatrix(theta_y[i], 2)

        for (j, x_local) in enumerate(chord_positions)
            x_relative = x_local - ea_local
            p_local = SVector(x_relative, 0.0, 0.0)
            panel_grid[:, j, i] = r_EA + R * p_local
        end
    end
    return panel_grid
end

function generate_propeller_blades_grid(total_radius::Float64, chord_length::Float64,
    ns::Int, nc::Int, blade_twists::Vector{Float64}, num_blades::Int)

    length(blade_twists) == ns+1 || error("blade_twists must have length ns+1")
    total_radius > 0 || error("Total radius must be positive")

    base_grid = zeros(3, nc+1, ns+1)
    radial_positions = range(0.0, total_radius, length=ns+1)
    chord_positions = range(0.0, chord_length, length=nc+1)

    for (i, r) in enumerate(radial_positions)
        theta = blade_twists[i]
        for (j, x_chord) in enumerate(chord_positions)
            x = x_chord * cos(theta)
            z = -x_chord * sin(theta)
            base_grid[:, j, i] = [x, r, z]
        end
    end

    blade_grids = [similar(base_grid) for _ in 1:num_blades]
    angular_step = 2pi / num_blades
    for (blade_idx, grid) in enumerate(blade_grids)
        phi = angular_step * (blade_idx - 1)
        cphi, sphi = cos(phi), sin(phi)
        for i in 1:size(base_grid, 3)
            for j in 1:size(base_grid, 2)
                x, y, z = base_grid[:, j, i]
                grid[:, j, i] = [x, y*cphi - z*sphi, y*sphi + z*cphi]
            end
        end
    end
    return blade_grids
end

