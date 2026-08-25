# ==============================================================================
# Aerodynamic grid utilities for the Bohnisch time-domain driver.
# ==============================================================================

using StaticArrays
using FLOWMath

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
    theta_x::Vector{Float64}, theta_y::Vector{Float64}, theta_z::Vector{Float64})

    panel_grid = zeros(Float64, 3, nc+1, ns+1)
    span_positions = range(0, span_length, length=ns+1)

    for (i, y) in enumerate(span_positions)
        local_chord = chord_distribution[i]
        local_xle = xle_distribution[i]
        chord_positions = range(0, local_chord, length=nc+1)

        ea_local = local_chord * 0.35
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
    theta_x::Vector{Float64}, theta_y::Vector{Float64}, theta_z::Vector{Float64})

    panel_grid = zeros(Float64, 3, nc+1, ns+1)
    span_positions = range(0, span_length, length=ns+1)

    for (i, y) in enumerate(span_positions)
        local_chord = chord_distribution[i]
        local_xle = xle_distribution[i]
        chord_positions = range(0.25*local_chord/nc, local_chord + 0.25*local_chord/nc, length=nc+1)

        ea_local = local_chord * 0.35
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


