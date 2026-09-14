# --- Functions to define blade geometry ---

"""
    get_chord_over_R(non_dim_radius::Float64) -> Float64

Calculates the chord/R ratio at a given non-dimensional radius.
Based on the provided graph, chord/R is constant at 0.3.

# Arguments
- `non_dim_radius`: Non-dimensional radius (0 at root, 1 at tip).

# Returns
- Chord/R ratio.
"""
function get_chord_over_R(non_dim_radius::Float64)
    if !(0.0 <= non_dim_radius <= 1.0)
        error("Non-dimensional radius must be between 0.0 and 1.0. Received: ", non_dim_radius)
    end
    return 0.3
end

"""
    get_twist_deg(non_dim_radius::Float64) -> Float64

Calculates the twist angle in degrees at a given non-dimensional radius.
The twist distribution is piecewise linear:
- T(r_nd) = 25 - 42.5 * r_nd for 0 <= r_nd <= 0.2
- T(r_nd) = -30 * r_nd + 22.5 for 0.2 < r_nd <= 1
This ensures T(0.75) = 0 degrees.

# Arguments
- `non_dim_radius`: Non-dimensional radius (0 at root, 1 at tip).

# Returns
- Twist angle in degrees.
"""
function get_twist_deg(non_dim_radius::Float64)
    if !(0.0 <= non_dim_radius <= 1.0)
        error("Non-dimensional radius must be between 0.0 and 1.0. Received: ", non_dim_radius)
    end

    r_nd = non_dim_radius
    local twist::Float64 # Ensure type stability for twist

    if r_nd <= 0.2
        twist = 25.0 - 42.5 * r_nd
    else # 0.2 < r_nd <= 1.0
        twist = -30.0 * r_nd + 22.5
    end
    return twist
end

function get_twist_deg_chang(non_dim_radius::Float64)
    if !(0.0 <= non_dim_radius <= 1.0)
        error("Non-dimensional radius must be between 0.0 and 1.0. Received: ", non_dim_radius)
    end

    # Linear twist from 44 deg (root) to 0 deg (tip)
    # Equation: twist = -44.0 * r_nd + 44.0
    # This can be factored to: 44.0 * (1.0 - r_nd)
    twist = 44.0 * (1.0 - non_dim_radius)
    
    return twist
end


# --- Function to get properties at panel nodes ---

"""
    get_nodal_properties(num_panels::Int) -> Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}

Calculates the non-dimensional radius, chord/R, and twist (in degrees) 
at each node of the spanwise panels. For N panels, there are N+1 nodes.

# Arguments
- `num_panels`: The number of spanwise panels.

# Returns
- `nodal_radii`: A vector of non-dimensional radii at each node.
- `nodal_chords_over_R`: A vector of chord/R ratios at each node (will be constant).
- `nodal_twists_deg`: A vector of twist angles (in degrees) at each node.
"""
function get_nodal_properties(num_panels::Int)
    if num_panels <= 0
        error("Number of panels must be positive. Received: ", num_panels)
    end

    N = num_panels
    num_nodes = N + 1
    
    nodal_radii = zeros(Float64, num_nodes)
    nodal_chords_over_R = zeros(Float64, num_nodes)
    nodal_twists_deg = zeros(Float64, num_nodes)

    for j in 0:N # Loop from node 0 to node N (N+1 nodes in total)
        # Calculate non-dimensional radius of the current node
        # Nodes are at r_nd = 0/N, 1/N, 2/N, ..., N/N
        r_node = Float64(j) / N 
        
        # Julia arrays are 1-indexed, so node j (0 to N) maps to array index j+1
        array_index = j + 1
        
        nodal_radii[array_index] = r_node
        nodal_chords_over_R[array_index] = get_chord_over_R(r_node)
        nodal_twists_deg[array_index] = get_twist_deg(r_node)
    end

    return nodal_radii, nodal_chords_over_R, nodal_twists_deg
end

using Interpolations

# --- 1. Define the Data ---

# Model Properties.xlsx, "Propeller blade": B2:B18 contains radial nodes
# in feet; E21:E36 assigns twist to the 16 segments between those nodes.
# Locate each sample at its segment midpoint, normalized by the tip radius.
const CHANG_RADIAL_NODES_FT = [
    0.0, 0.235813, 0.471625, 0.707437, 0.943250, 1.179063,
    1.414875, 1.650687, 1.886500, 2.122313, 2.358125, 2.593937,
    2.829750, 3.065563, 3.301375, 3.537187, 3.773000,
]
const R_DATA = (CHANG_RADIAL_NODES_FT[1:end-1] .+ CHANG_RADIAL_NODES_FT[2:end]) ./
    (2 * CHANG_RADIAL_NODES_FT[end])

# Preserve the tabulated angles converted to degrees; only their locations
# change. Do not replace them with rounded nominal root/tip angles.
const TWIST_DATA = [
    62.59324988,
    59.84464529,
    57.09603497,
    54.34743037,
    51.59882578,
    48.85022118,
    46.10161659,
    43.35300627,
    40.60440167,
    37.85579708,
    35.10719249,
    32.35858216,
    29.60997757,
    26.86137297,
    24.11276838,
    21.36416379
]

# Linearly interpolate between segment centers and extrapolate over the
# half-segments at the root and tip to obtain the UVLM nodal twist angles.
const twist_interpolator = linear_interpolation(R_DATA, TWIST_DATA, extrapolation_bc=Line())


# --- 2. New Geometry Function ---

"""
    get_twist_deg_interp(non_dim_radius::Float64) -> Float64

Calculate the blade angle in degrees from the 16 segment-midpoint samples.
The root and tip values are linearly extrapolated from the nearest samples.
"""
function get_twist_deg_interp(non_dim_radius::Float64)
    # Sanity check for bounds (with small tolerance for floating point errors)
    if non_dim_radius < -1e-6 || non_dim_radius > 1.0 + 1e-6
        error("Non-dimensional radius must be between 0.0 and 1.0. Received: ", non_dim_radius)
    end

    # The interpolator object can be called like a function
    return twist_interpolator(non_dim_radius)
end


# --- 3. New Nodal Properties Function ---

"""
    get_nodal_properties_chang(num_panels::Int) -> Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}

Generate nodal properties from the Chang segment-midpoint twist distribution.
"""
function get_nodal_properties_chang(num_panels::Int)
    if num_panels <= 0
        error("Number of panels must be positive. Received: ", num_panels)
    end

    N = num_panels
    num_nodes = N + 1
    
    nodal_radii = zeros(Float64, num_nodes)
    nodal_chords_over_R = zeros(Float64, num_nodes)
    nodal_twists_deg = zeros(Float64, num_nodes)

    for j in 0:N 
        r_node = Float64(j) / N 
        array_index = j + 1
        
        nodal_radii[array_index] = r_node
        nodal_chords_over_R[array_index] = get_chord_over_R(r_node) # Uses your existing constant chord function
        nodal_twists_deg[array_index] = get_twist_deg_interp(r_node) # Uses new interpolation
    end

    return nodal_radii, nodal_chords_over_R, nodal_twists_deg
end

## --- Example Usage for 15 panels ---
#num_example_panels = 20
#radii_at_nodes, chords_R_at_nodes, twists_at_nodes = get_nodal_properties(num_example_panels)
#
#println("--- Nodal Properties for ", num_example_panels, " Panels (", num_example_panels + 1, " nodes) ---")
#
#println("\nNodal Non-dimensional Radii (r_nd):")
#display(radii_at_nodes)
#
## Chord/R is constant, but shown for completeness
## println("\nNodal Chords/R:")
## display(chords_R_at_nodes) 
#
#println("\nNodal Twists (degrees):")
#display(twists_at_nodes)

# You can also print them in a more formatted way if desired:
# println("\nDetailed Nodal Values:")
# for i in 1:(num_example_panels + 1)
#     println("Node ", i-1, # 0-indexed node number
#             ": r_nd=", round(radii_at_nodes[i], digits=4), 
#             ", chord/R=", round(chords_R_at_nodes[i], digits=2), 
#             ", twist=", round(twists_at_nodes[i], digits=3), " deg")
# end
