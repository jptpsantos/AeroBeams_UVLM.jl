# Structural matrix assembly for the Chang wing-propeller time-domain model.
#
# Wing structural DOF order per node:
#   [u_span, v_chord, w_vertical_down, theta_span, theta_chord, theta_vertical]
#
# Propeller DOF order:
#   [pitch_about_span, yaw_about_vertical]

using LinearAlgebra

function chang_tilde(v)
    return [ 0.0  -v[3]   v[2]
             v[3]   0.0  -v[1]
            -v[2]   v[1]   0.0]
end

"""
    chang_point_mass_matrix(mass, inertia, offset; inertia_reference=:center_of_mass)

Return the consistent 6-by-6 rigid point-mass matrix at a beam node. `offset`
points from the beam reference axis to the mass center. By default, `inertia`
is interpreted about the mass center and is shifted to the beam node with the
parallel-axis theorem. Set `inertia_reference=:beam_axis` only to reproduce
the former implementation, where the supplied inertia was used directly.
"""
function chang_point_mass_matrix(
    mass::Real,
    inertia::AbstractMatrix,
    offset::AbstractVector;
    inertia_reference::Symbol = :center_of_mass,
)
    size(inertia) == (3, 3) || throw(DimensionMismatch("inertia must be 3-by-3"))
    length(offset) == 3 || throw(DimensionMismatch("offset must have three components"))
    mass >= 0 || throw(ArgumentError("mass must be nonnegative"))
    inertia_reference in (:center_of_mass, :beam_axis) || throw(ArgumentError(
        "inertia_reference must be :center_of_mass or :beam_axis",
    ))

    skew_offset = chang_tilde(offset)
    rotational_inertia = inertia_reference == :center_of_mass ?
        Matrix(inertia) .- mass .* (skew_offset * skew_offset) : Matrix(inertia)
    identity3 = Matrix{promote_type(Float64, eltype(inertia))}(I, 3, 3)
    return [mass .* identity3       -mass .* skew_offset;
            mass .* skew_offset     rotational_inertia]
end

function chang_beam_strain_displacement_matrix(L::Real, xi::Real)
    B = zeros(4, 12)
    B[1, 1]  = -1 / L
    B[1, 7]  =  1 / L

    B[2, 2]  =  (6xi) / L^2
    B[2, 6]  =  (3xi - 1) / L
    B[2, 8]  = -(6xi) / L^2
    B[2, 12] =  (3xi + 1) / L

    B[3, 3]  =  (6xi) / L^2
    B[3, 5]  =  (3xi - 1) / L
    B[3, 9]  = -(6xi) / L^2
    B[3, 11] =  (3xi + 1) / L

    B[4, 4]  = -1 / L
    B[4, 10] =  1 / L
    return B
end

function chang_beam_element_stiffness_matrix(L::Real, C::AbstractMatrix)
    # B'CB is quadratic in the parent coordinate, so the two-point
    # Gauss-Legendre rule integrates this element matrix exactly for constant C.
    gauss_points = (-inv(sqrt(3.0)), inv(sqrt(3.0)))
    K = zeros(12, 12)
    Csym = 0.5 .* (Matrix{Float64}(C) .+ Matrix{Float64}(C)')

    for xi in gauss_points
        B = chang_beam_strain_displacement_matrix(L, xi)
        K .+= B' * Csym * B * (L / 2)
    end

    return 0.5 .* (K .+ K')
end

function chang_linear_interpolate_scalar(nodes, values, query)
    length(nodes) == length(values) || throw(DimensionMismatch(
        "Interpolation nodes and values must have the same length",
    ))
    query <= nodes[1] && return values[1]
    query >= nodes[end] && return values[end]
    right = searchsortedfirst(nodes, query)
    left = right - 1
    fraction = (query - nodes[left]) / (nodes[right] - nodes[left])
    return (1 - fraction) * values[left] + fraction * values[right]
end

function chang_constitutive_matrix_at_eta(stiffness_distribution, eta)
    stations = stiffness_distribution.eta
    return [
        chang_linear_interpolate_scalar(stations, stiffness_distribution.EA, eta) 0.0 0.0 0.0;
        0.0 chang_linear_interpolate_scalar(stations, stiffness_distribution.EIz, eta) chang_linear_interpolate_scalar(stations, stiffness_distribution.EIzy, eta) 0.0;
        0.0 chang_linear_interpolate_scalar(stations, stiffness_distribution.EIzy, eta) chang_linear_interpolate_scalar(stations, stiffness_distribution.EIy, eta) 0.0;
        0.0 0.0 0.0 chang_linear_interpolate_scalar(stations, stiffness_distribution.GJ, eta)
    ]
end

"""
    chang_distributed_beam_element_stiffness_matrix(
        x_left, x_right, span_length, stiffness_distribution,
    )

Integrate B'C(x)B while splitting the element at every tabulated stiffness
station. This removes the target-mesh alignment error caused by assigning one
element-center property to an element that crosses a sharp stiffness change.
"""
function chang_distributed_beam_element_stiffness_matrix(
    x_left::Real,
    x_right::Real,
    span_length::Real,
    stiffness_distribution,
)
    x_right > x_left || throw(ArgumentError("Element length must be positive"))
    span_length > 0 || throw(ArgumentError("Wing span must be positive"))
    element_length = x_right - x_left
    element_center = (x_left + x_right) / 2
    stiffness_stations = stiffness_distribution.eta .* span_length
    internal_breaks = filter(x -> x_left < x < x_right, stiffness_stations)
    integration_edges = unique(sort([x_left; internal_breaks; x_right]))
    gauss_points = (-sqrt(3 / 5), 0.0, sqrt(3 / 5))
    gauss_weights = (5 / 9, 8 / 9, 5 / 9)
    K = zeros(12, 12)

    for interval in 1:(length(integration_edges) - 1)
        interval_left = integration_edges[interval]
        interval_right = integration_edges[interval + 1]
        interval_length = interval_right - interval_left
        interval_center = (interval_left + interval_right) / 2

        for (gauss_point, gauss_weight) in zip(gauss_points, gauss_weights)
            x = interval_center + gauss_point * interval_length / 2
            xi = 2 * (x - element_center) / element_length
            B = chang_beam_strain_displacement_matrix(element_length, xi)
            C = chang_constitutive_matrix_at_eta(
                stiffness_distribution,
                x / span_length,
            )
            K .+= B' * C * B * (gauss_weight * interval_length / 2)
        end
    end
    return 0.5 .* (K .+ K')
end

function assemble_chang_structural_matrices(;
    Ne, le, ndof, NDOF, nnodes,
    EIy_vec, EIz_vec, EIzy_vec, GJ_vec, EA_vec,
    span_nodes = nothing,
    span_length = nothing,
    stiffness_distribution = nothing,
    m_node_vec, Ixx_node_vec, Iyy_node_vec, Izz_node_vec,
    Ixy_node_vec, Ixz_node_vec, Iyz_node_vec,
    cg_x_node_vec, cg_y_node_vec, cg_z_node_vec,
    spatial_inertia_blocks = nothing,
    ndof_P, Npropellers, prop_attach_nodes,
    prop_attachment_node_pairs = nothing,
    prop_attachment_weights = nothing,
    Inθ_prop, Inψ_prop, Kθ_prop, Kψ_prop, ξ_prop,
    stiffness_damping_ratio, stiffness_damping_reference_omega,
    Ix_prop, Ω,
    mP_prop, SθP_prop, SψP_prop, SαP_prop, SγP_prop,
    IθαP_prop, IψγP_prop, IαP_prop, IγP_prop,
    inertia_reference::Symbol = :center_of_mass)

    Ks_W = zeros(NDOF, NDOF)
    Ms_W = zeros(NDOF, NDOF)
    Cs_W = zeros(NDOF, NDOF)

    use_distributed_stiffness = !isnothing(stiffness_distribution)
    if use_distributed_stiffness
        isnothing(span_nodes) && error(
            "span_nodes are required for distributed stiffness integration",
        )
        isnothing(span_length) && error(
            "span_length is required for distributed stiffness integration",
        )
        length(span_nodes) == nnodes || throw(DimensionMismatch(
            "span_nodes must contain nnodes entries",
        ))
    end

    for i = 1:Ne
        ks_W = if use_distributed_stiffness
            chang_distributed_beam_element_stiffness_matrix(
                span_nodes[i],
                span_nodes[i + 1],
                span_length,
                stiffness_distribution,
            )
        else
            D = [EA_vec[i]   0.0         0.0         0.0;
                 0.0         EIz_vec[i]  EIzy_vec[i] 0.0;
                 0.0         EIzy_vec[i] EIy_vec[i]  0.0;
                 0.0         0.0         0.0         GJ_vec[i]]
            chang_beam_element_stiffness_matrix(le, D)
        end
        idx_start = ndof * (i - 1) + 1
        idx_end = ndof * i + ndof
        Ks_W[idx_start:idx_end, idx_start:idx_end] += ks_W
    end

    if !isnothing(spatial_inertia_blocks)
        length(spatial_inertia_blocks) == nnodes || throw(DimensionMismatch(
            "spatial_inertia_blocks must contain nnodes entries",
        ))
    end
    for n = 1:nnodes
        ms_W = if isnothing(spatial_inertia_blocks)
            m_val = m_node_vec[n]
            eta = [cg_x_node_vec[n]; cg_y_node_vec[n]; cg_z_node_vec[n]]
            inertia_matrix = [Ixx_node_vec[n] Ixy_node_vec[n] Ixz_node_vec[n];
                              Ixy_node_vec[n] Iyy_node_vec[n] Iyz_node_vec[n];
                              Ixz_node_vec[n] Iyz_node_vec[n] Izz_node_vec[n]]
            chang_point_mass_matrix(
                m_val,
                inertia_matrix,
                eta;
                inertia_reference,
            )
        else
            block = Matrix{Float64}(spatial_inertia_blocks[n])
            size(block) == (6, 6) || throw(DimensionMismatch(
                "Every spatial inertia block must be 6-by-6",
            ))
            0.5 .* (block .+ block')
        end
        Ms_W[6n-5:6n, 6n-5:6n] += ms_W
    end

    total_prop_dofs = ndof_P
    Ms_p_local = [Inθ_prop 0.0; 0.0 Inψ_prop]
    Ks_p_local = [Kθ_prop 0.0; 0.0 Kψ_prop]
    omega_local = sqrt.(diag(Ks_p_local ./ Ms_p_local))
    Cvisc = [2 * ξ_prop * Inθ_prop * omega_local[1] 0.0;
             0.0 2 * ξ_prop * Inψ_prop * omega_local[2]]
    Cgyro = [0.0 Ix_prop * Ω; -Ix_prop * Ω 0.0]
    Cs_p_local = Cvisc + Cgyro

    Ms_P = zeros(total_prop_dofs, total_prop_dofs)
    Cs_P = zeros(total_prop_dofs, total_prop_dofs)
    Ks_P = zeros(total_prop_dofs, total_prop_dofs)
    Bs_P = zeros(total_prop_dofs, NDOF)
    Ds_P = zeros(total_prop_dofs, NDOF)
    Fs_W = zeros(NDOF, total_prop_dofs)
    Gs_W = zeros(NDOF, NDOF)
    Hs_W = zeros(NDOF, total_prop_dofs)

    attach_dofs_all = Vector{UnitRange{Int}}(undef, Npropellers)
    attachment_operators = Vector{Matrix{Float64}}(undef, Npropellers)
    use_interpolated_attachment =
        !isnothing(prop_attachment_node_pairs) || !isnothing(prop_attachment_weights)
    if use_interpolated_attachment
        isnothing(prop_attachment_node_pairs) && error(
            "prop_attachment_node_pairs and prop_attachment_weights must be provided together",
        )
        isnothing(prop_attachment_weights) && error(
            "prop_attachment_node_pairs and prop_attachment_weights must be provided together",
        )
        length(prop_attachment_node_pairs) == Npropellers || throw(DimensionMismatch(
            "prop_attachment_node_pairs must contain Npropellers entries",
        ))
        length(prop_attachment_weights) == Npropellers || throw(DimensionMismatch(
            "prop_attachment_weights must contain Npropellers entries",
        ))
    end

    for ip in 1:Npropellers
        pos = prop_attach_nodes[ip]
        prop_dofs = (2 * (ip - 1) + 1):(2 * ip)
        node_dofs = (ndof * (pos - 1) + 1):(ndof * pos)
        attach_dofs_all[ip] = node_dofs

        attachment_operator = zeros(6, NDOF)
        if use_interpolated_attachment
            left_node, right_node = prop_attachment_node_pairs[ip]
            left_weight, right_weight = prop_attachment_weights[ip]
            1 <= left_node <= nnodes || throw(ArgumentError(
                "Left attachment node $left_node is outside the wing mesh",
            ))
            1 <= right_node <= nnodes || throw(ArgumentError(
                "Right attachment node $right_node is outside the wing mesh",
            ))
            isapprox(left_weight + right_weight, 1.0; atol = 1e-12) ||
                throw(ArgumentError("Attachment weights must sum to one"))
            left_dofs = (ndof * (left_node - 1) + 1):(ndof * left_node)
            right_dofs = (ndof * (right_node - 1) + 1):(ndof * right_node)
            attachment_operator[:, left_dofs] .+=
                left_weight .* Matrix{Float64}(I, 6, 6)
            attachment_operator[:, right_dofs] .+=
                right_weight .* Matrix{Float64}(I, 6, 6)
        else
            attachment_operator[:, node_dofs] .= Matrix{Float64}(I, 6, 6)
        end
        attachment_operators[ip] = attachment_operator

        Ms_P[prop_dofs, prop_dofs] .= Ms_p_local
        Cs_P[prop_dofs, prop_dofs] .= Cs_p_local
        Ks_P[prop_dofs, prop_dofs] .= Ks_p_local

        Blocal = [0.0       0.0  SθP_prop  IθαP_prop  0.0  0.0;
                 -SψP_prop  0.0  0.0       0.0        0.0  IψγP_prop]
        Bs_P[prop_dofs, :] .+= Blocal * attachment_operator
        Fs_W[:, prop_dofs] .+= attachment_operator' * Blocal'

        Glocal = zeros(6, 6)
        Glocal[1, 1] = mP_prop
        Glocal[2, 2] = mP_prop
        Glocal[3, 3] = mP_prop
        Glocal[3, 4] = SαP_prop
        Glocal[4, 3] = SαP_prop
        Glocal[4, 4] = IαP_prop
        Glocal[1, 6] = -SγP_prop
        Glocal[6, 1] = -SγP_prop
        Glocal[6, 6] = IγP_prop
        Gs_W .+= attachment_operator' * Glocal * attachment_operator

        Hlocal = [0.0 0.0;
                  0.0 0.0;
                  0.0 0.0;
                  0.0 Ix_prop * Ω;
                  0.0 0.0;
                 -Ix_prop * Ω 0.0]
        Hs_W[:, prop_dofs] .+= attachment_operator' * Hlocal
        Ds_P[prop_dofs, :] .+= -Hlocal' * attachment_operator

        local_wing_gyro = [0 0 0 0 0 0;
                           0 0 0 0 0 0;
                           0 0 0 0 0 0;
                           0 0 0 0 0 Ix_prop * Ω;
                           0 0 0 0 0 0;
                           0 0 0 -Ix_prop * Ω 0 0]
        Cs_W .+= attachment_operator' * local_wing_gyro * attachment_operator
    end

    ndof_total_full = NDOF + ndof_P
    M_global = [Ms_W + Gs_W  Fs_W; Bs_P Ms_P]
    C_global = [Cs_W Hs_W; Ds_P Cs_P]
    K_global = [Ks_W zeros(NDOF, total_prop_dofs); zeros(total_prop_dofs, NDOF) Ks_P]

    free_dofs = (ndof + 1):ndof_total_full
    M = M_global[free_dofs, free_dofs]
    K = K_global[free_dofs, free_dofs]
    0.0 <= stiffness_damping_ratio < 1.0 || throw(ArgumentError(
        "stiffness_damping_ratio must lie in [0, 1)",
    ))
    stiffness_damping_reference_omega > 0.0 || throw(ArgumentError(
        "stiffness_damping_reference_omega must be positive",
    ))
    stiffness_damping_coefficient =
        2.0 * stiffness_damping_ratio / stiffness_damping_reference_omega
    C_rayleigh = stiffness_damping_coefficient .* K
    C = C_global[free_dofs, free_dofs] + C_rayleigh

    return (
        Ks_W = Ks_W,
        Ms_W = Ms_W,
        Cs_W = Cs_W,
        Ms_P = Ms_P,
        Cs_P = Cs_P,
        Ks_P = Ks_P,
        Bs_P = Bs_P,
        Ds_P = Ds_P,
        Fs_W = Fs_W,
        Gs_W = Gs_W,
        Hs_W = Hs_W,
        attach_dofs_all = attach_dofs_all,
        attachment_operators = attachment_operators,
        M_global = M_global,
        C_global = C_global,
        K_global = K_global,
        free_dofs = free_dofs,
        M = M,
        C = C,
        K = K,
        stiffness_damping_ratio,
        stiffness_damping_reference_omega,
        stiffness_damping_coefficient,
        ndof_free = length(free_dofs),
        ndof_wing_free = NDOF - ndof,
        ndof_prop_free = ndof_P,
        inertia_reference,
        mass_model = isnothing(spatial_inertia_blocks) ?
            :component_arrays : :control_volume_spatial_blocks,
        stiffness_model = use_distributed_stiffness ?
            :distributed_integrated : :element_center,
    )
end
