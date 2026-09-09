# Case-level construction and diagnostics for the Chang structural model.
# The numerical element/mass routines remain in `chang_structural_matrices.jl`.

isdefined(@__MODULE__, :assemble_chang_structural_matrices) ||
    include(joinpath(@__DIR__, "chang_structural_matrices.jl"))

function assemble_chang_structural_model(parameters;
    inertia_reference::Symbol = :center_of_mass,
)
    matrices = assemble_chang_structural_matrices(
        Ne = parameters.Ne,
        le = parameters.le,
        ndof = parameters.ndof,
        NDOF = parameters.NDOF,
        nnodes = parameters.nnodes,
        EIy_vec = parameters.EIy_vec,
        EIz_vec = parameters.EIz_vec,
        EIzy_vec = parameters.EIzy_vec,
        GJ_vec = parameters.GJ_vec,
        EA_vec = parameters.EA_vec,
        span_nodes = parameters.span_nodes,
        span_length = parameters.span_length,
        stiffness_distribution = parameters.wing_stiffness_distribution,
        m_node_vec = parameters.m_node_vec,
        Ixx_node_vec = parameters.Ixx_node_vec,
        Iyy_node_vec = parameters.Iyy_node_vec,
        Izz_node_vec = parameters.Izz_node_vec,
        Ixy_node_vec = parameters.Ixy_node_vec,
        Ixz_node_vec = parameters.Ixz_node_vec,
        Iyz_node_vec = parameters.Iyz_node_vec,
        cg_x_node_vec = parameters.cg_x_node_vec,
        cg_y_node_vec = parameters.cg_y_node_vec,
        cg_z_node_vec = parameters.cg_z_node_vec,
        spatial_inertia_blocks = inertia_reference == :center_of_mass ?
            parameters.spatial_inertia_node_blocks : nothing,
        ndof_P = parameters.ndof_P,
        Npropellers = parameters.Npropellers,
        prop_attach_nodes = parameters.prop_attach_nodes,
        prop_attachment_node_pairs = parameters.prop_attachment_node_pairs,
        prop_attachment_weights = parameters.prop_attachment_weights,
        Inθ_prop = parameters.Inθ_prop,
        Inψ_prop = parameters.Inψ_prop,
        Kθ_prop = parameters.Kθ_prop,
        Kψ_prop = parameters.Kψ_prop,
        ξ_prop = parameters.ξ_prop,
        stiffness_damping_ratio = parameters.stiffness_damping_ratio,
        stiffness_damping_reference_omega = parameters.stiffness_damping_reference_omega,
        Ix_prop = parameters.Ix_prop,
        Ω = parameters.Ω,
        mP_prop = parameters.mP_prop,
        SθP_prop = parameters.SθP_prop,
        SψP_prop = parameters.SψP_prop,
        SαP_prop = parameters.SαP_prop,
        SγP_prop = parameters.SγP_prop,
        IθαP_prop = parameters.IθαP_prop,
        IψγP_prop = parameters.IψγP_prop,
        IαP_prop = parameters.IαP_prop,
        IγP_prop = parameters.IγP_prop,
        inertia_reference = inertia_reference,
    )
    return merge(matrices, (; dofs_per_node = parameters.ndof))
end

"""
    chang_wing_modal_frequencies(structural; number_of_modes=6)

Return the lowest fixed-root, wing-only undamped natural frequencies in hertz.
The propeller/pylon coordinates are intentionally excluded so these values can
be compared with the published isolated-wing modal targets.
"""
function chang_wing_modal_frequencies(structural; number_of_modes::Int = 6)
    return [
        mode.frequency_hz for mode in chang_wing_modal_analysis(
            structural;
            number_of_modes,
        )
    ]
end

"""
    chang_wing_modal_analysis(structural; number_of_modes=6)

Solve the fixed-root, wing-only eigenproblem and classify every mode from the
dominant strain-energy family. The families follow the structural DOF order:
out-of-plane bending uses `(w_vertical, theta_chord)`, in-plane bending uses
`(v_chord, theta_vertical)`, torsion uses `theta_span`, and axial deformation
uses `u_span`.
"""
function chang_wing_modal_analysis(structural; number_of_modes::Int = 6)
    number_of_modes > 0 || throw(ArgumentError("number_of_modes must be positive"))
    ndof = structural.dofs_per_node
    wing_free_dofs = (ndof + 1):size(structural.Ms_W, 1)
    wing_mass = Symmetric(structural.Ms_W[wing_free_dofs, wing_free_dofs])
    wing_stiffness = Symmetric(structural.Ks_W[wing_free_dofs, wing_free_dofs])
    solution = eigen(wing_stiffness, wing_mass)
    positive_indices = findall(>(0.0), real.(solution.values))
    positive_indices = positive_indices[sortperm(real.(solution.values[positive_indices]))]
    count = min(number_of_modes, length(positive_indices))

    family_dofs = (
        out_of_plane = reduce(vcat, (node_start .+ [3, 5] for node_start in 0:ndof:(length(wing_free_dofs) - ndof))),
        in_plane = reduce(vcat, (node_start .+ [2, 6] for node_start in 0:ndof:(length(wing_free_dofs) - ndof))),
        torsion = collect(4:ndof:length(wing_free_dofs)),
        axial = collect(1:ndof:length(wing_free_dofs)),
    )
    stiffness_matrix = Matrix(wing_stiffness)
    family_counts = Dict(name => 0 for name in keys(family_dofs))
    modes = NamedTuple[]

    for eigen_index in positive_indices[1:count]
        shape = solution.vectors[:, eigen_index]
        family_energies = map(family_dofs) do indices
            projected_shape = zeros(eltype(shape), length(shape))
            projected_shape[indices] .= shape[indices]
            return max(real(dot(projected_shape, stiffness_matrix * projected_shape)), 0.0)
        end
        family_names = collect(keys(family_energies))
        family_values = collect(values(family_energies))
        family = family_names[argmax(family_values)]
        family_counts[family] += 1
        push!(modes, (;
            frequency_hz = sqrt(real(solution.values[eigen_index])) / (2π),
            family,
            family_order = family_counts[family],
            family_energies,
        ))
    end

    return modes
end

function chang_structural_diagnostics(structural)
    mass_symmetry_error = norm(structural.M - structural.M', Inf) /
        max(norm(structural.M, Inf), 1.0)
    stiffness_symmetry_error = norm(structural.K - structural.K', Inf) /
        max(norm(structural.K, Inf), 1.0)
    minimum_mass_eigenvalue = minimum(eigvals(Symmetric(structural.M)))
    minimum_stiffness_eigenvalue = minimum(eigvals(Symmetric(structural.K)))
    damping_symmetric_part = 0.5 .* (structural.C .+ structural.C')
    return (;
        mass_symmetry_error,
        stiffness_symmetry_error,
        minimum_mass_eigenvalue,
        minimum_stiffness_eigenvalue,
        damping_symmetric_part_norm = norm(damping_symmetric_part, Inf),
        wing_modal_frequencies_hz = chang_wing_modal_frequencies(structural),
    )
end
