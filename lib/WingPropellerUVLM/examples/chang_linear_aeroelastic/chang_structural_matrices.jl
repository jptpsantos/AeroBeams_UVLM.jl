# ==============================================================================
# Structural matrix assembly for the Chang wing-propeller time-domain model.
#
# Wing structural DOF order per node:
#   [u_span, v_chord, w_vertical_down, theta_span, theta_chord, theta_vertical]
#
# Propeller DOF order:
#   [pitch_about_span, yaw_about_vertical]
# ==============================================================================

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

function chang_beam_element_stiffness_matrix(L::Float64, C::Matrix{Float64})
    # B'CB is quadratic in the parent coordinate, so the two-point
    # Gauss-Legendre rule integrates this element matrix exactly.
    gauss_points = (-inv(sqrt(3.0)), inv(sqrt(3.0)))
    gauss_weights = (1.0, 1.0)
    K = zeros(12, 12)
    Csym = 0.5 .* (Matrix{Float64}(C) .+ Matrix{Float64}(C)')

    for (xi, w) in zip(gauss_points, gauss_weights)
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

        K += B' * Csym * B * (w * L / 2)
    end

    return 0.5 .* (K .+ K')
end

function assemble_chang_structural_matrices(;
    Ne, le, ndof, NDOF, nnodes,
    EIy_vec, EIz_vec, EIzy_vec, GJ_vec, EA_vec,
    m_node_vec, Ixx_node_vec, Iyy_node_vec, Izz_node_vec,
    Ixy_node_vec, Ixz_node_vec, Iyz_node_vec,
    cg_x_node_vec, cg_y_node_vec, cg_z_node_vec,
    ndof_P, Npropellers, prop_attach_nodes,
    Inθ_prop, Inψ_prop, Kθ_prop, Kψ_prop, ξ_prop,
    stiffness_damping_ratio, stiffness_damping_reference_omega,
    Ix_prop, Ω,
    mP_prop, SθP_prop, SψP_prop, SαP_prop, SγP_prop,
    IθαP_prop, IψγP_prop, IαP_prop, IγP_prop,
    inertia_reference::Symbol = :center_of_mass)

    Ks_W = zeros(NDOF, NDOF)
    Ms_W = zeros(NDOF, NDOF)
    Cs_W = zeros(NDOF, NDOF)

    for i = 1:Ne
        D = [EA_vec[i]   0.0         0.0         0.0;
             0.0         EIz_vec[i]  EIzy_vec[i] 0.0;
             0.0         EIzy_vec[i] EIy_vec[i]  0.0;
             0.0         0.0         0.0         GJ_vec[i]]

        ks_W = chang_beam_element_stiffness_matrix(le, D)
        idx_start = ndof * (i - 1) + 1
        idx_end = ndof * i + ndof
        Ks_W[idx_start:idx_end, idx_start:idx_end] += ks_W
    end

    for n = 1:nnodes
        m_val = m_node_vec[n]
        eta = [cg_x_node_vec[n]; cg_y_node_vec[n]; cg_z_node_vec[n]]
        inertia_matrix = [Ixx_node_vec[n] Ixy_node_vec[n] Ixz_node_vec[n];
                          Ixy_node_vec[n] Iyy_node_vec[n] Iyz_node_vec[n];
                          Ixz_node_vec[n] Iyz_node_vec[n] Izz_node_vec[n]]
        ms_W = chang_point_mass_matrix(
            m_val,
            inertia_matrix,
            eta;
            inertia_reference,
        )
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

    for ip in 1:Npropellers
        pos = prop_attach_nodes[ip]
        prop_dofs = (2 * (ip - 1) + 1):(2 * ip)
        node_dofs = (ndof * (pos - 1) + 1):(ndof * pos)
        attach_dofs_all[ip] = node_dofs

        Ms_P[prop_dofs, prop_dofs] .= Ms_p_local
        Cs_P[prop_dofs, prop_dofs] .= Cs_p_local
        Ks_P[prop_dofs, prop_dofs] .= Ks_p_local

        Blocal = [0.0       0.0  SθP_prop  IθαP_prop  0.0  0.0;
                 -SψP_prop  0.0  0.0       0.0        0.0  IψγP_prop]
        Bs_P[prop_dofs, node_dofs] .= Blocal
        Fs_W[node_dofs, prop_dofs] .= Blocal'

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
        Gs_W[node_dofs, node_dofs] .+= Glocal

        Hlocal = [0.0 0.0;
                  0.0 0.0;
                  0.0 0.0;
                  0.0 Ix_prop * Ω;
                  0.0 0.0;
                 -Ix_prop * Ω 0.0]
        Hs_W[node_dofs, prop_dofs] .= Hlocal
        Ds_P[prop_dofs, node_dofs] .= -Hlocal'

        Cs_W[node_dofs, node_dofs] += [0 0 0 0 0 0;
                                       0 0 0 0 0 0;
                                       0 0 0 0 0 0;
                                       0 0 0 0 0 Ix_prop * Ω;
                                       0 0 0 0 0 0;
                                       0 0 0 -Ix_prop * Ω 0 0]
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
    )
end
