import Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..", "..", "..")))

using LinearAlgebra
using Printf
using WingPropellerUVLM
include(joinpath(@__DIR__, "..", "src", "WingPropellerAeroelastic.jl"))
using .WingPropellerAeroelastic

function reference_bohnisch_matrices(case)
    span_nodes = case.wing.span_nodes
    element_count = length(span_nodes) - 1
    node_count = length(span_nodes)
    element_length = span_nodes[2] - span_nodes[1]
    reduced_dofs_per_node = 3
    wing_dof_count = reduced_dofs_per_node * node_count
    wing_mass = zeros(wing_dof_count, wing_dof_count)
    wing_stiffness = zeros(wing_dof_count, wing_dof_count)
    rho = case.wing.inertia.mass_per_length
    rhoJ = case.wing.inertia.torsional_inertia_per_length
    offset = case.wing.inertia.cg_offset_chord
    EI = case.wing.stiffness.EI_out_of_plane
    GJ = case.wing.stiffness.GJ

    for element in 1:element_count
        length = element_length
        element_mass = rho * length / 420 .* [
            156 22length 0 54 -13length 0;
            22length 4length^2 0 13length -3length^2 0;
            0 0 0 0 0 0;
            54 13length 0 156 -22length 0;
            -13length -3length^2 0 -22length 4length^2 0;
            0 0 0 0 0 0
        ] + length / 420 .* [
            0 0 0 0 0 0;
            0 0 0 0 0 0;
            0 0 140rhoJ 0 0 70rhoJ;
            0 0 0 0 0 0;
            0 0 0 0 0 0;
            0 0 70rhoJ 0 0 140rhoJ
        ] + rho * length * offset .* [
            0 0 7/20 0 0 3/20;
            0 0 length/20 0 0 length/30;
            7/20 length/20 0 3/20 -length/30 0;
            0 0 3/20 0 0 7/20;
            0 0 -length/30 0 0 -length/20;
            3/20 length/30 0 7/20 -length/20 0
        ]
        element_stiffness = EI / length^3 .* [
            12 6length 0 -12 6length 0;
            6length 4length^2 0 -6length 2length^2 0;
            0 0 0 0 0 0;
            -12 -6length 0 12 -6length 0;
            6length 2length^2 0 -6length 4length^2 0;
            0 0 0 0 0 0
        ] + GJ / length .* [
            0 0 0 0 0 0;
            0 0 0 0 0 0;
            0 0 1 0 0 -1;
            0 0 0 0 0 0;
            0 0 0 0 0 0;
            0 0 -1 0 0 1
        ]
        dofs = (reduced_dofs_per_node * (element - 1) + 1):(
            reduced_dofs_per_node * (element + 1)
        )
        wing_mass[dofs, dofs] .+= element_mass
        wing_stiffness[dofs, dofs] .+= element_stiffness
    end
    wing_damping = 0.001 .* wing_stiffness

    propeller = only(case.propellers)
    omega = propeller_angular_speed(propeller, case.operating_condition)
    propeller_mass = [propeller.pitch_inertia 0.0; 0.0 propeller.yaw_inertia]
    propeller_stiffness = [propeller.pitch_stiffness 0.0; 0.0 propeller.yaw_stiffness]
    propeller_damping = [
        0.0 propeller.spin_inertia * omega;
        -propeller.spin_inertia * omega 0.0
    ]
    B = zeros(2, wing_dof_count)
    F = zeros(wing_dof_count, 2)
    G = zeros(wing_dof_count, wing_dof_count)
    H = zeros(wing_dof_count, 2)
    D = zeros(2, wing_dof_count)
    node_dofs = (reduced_dofs_per_node * (propeller.attachment_node - 1) + 1):(
        reduced_dofs_per_node * propeller.attachment_node
    )
    B[:, node_dofs] .= [
        propeller.pitch_first_moment 0.0 propeller.pitch_cross_inertia;
        0.0 0.0 0.0
    ]
    F[node_dofs, :] .= B[:, node_dofs]'
    G[node_dofs, node_dofs] .+= [
        propeller.mass 0.0 propeller.wing_pitch_first_moment;
        0.0 0.0 0.0;
        propeller.wing_pitch_first_moment 0.0 propeller.wing_pitch_inertia
    ]
    H[node_dofs, :] .= [
        0.0 0.0;
        0.0 0.0;
        0.0 propeller.spin_inertia * omega
    ]
    D[:, node_dofs] .= [
        0.0 0.0 0.0;
        0.0 0.0 -propeller.spin_inertia * omega
    ]

    M_global = [wing_mass + G F; B propeller_mass]
    C_global = [wing_damping H; D propeller_damping]
    K_global = [wing_stiffness zeros(wing_dof_count, 2);
                zeros(2, wing_dof_count) propeller_stiffness]
    free_dofs = (reduced_dofs_per_node + 1):size(M_global, 1)
    return (;
        M_global,
        C_global,
        K_global,
        M = M_global[free_dofs, free_dofs],
        C = C_global[free_dofs, free_dofs],
        K = K_global[free_dofs, free_dofs],
    )
end

function state_space_modes(structural)
    n = size(structural.M, 1)
    matrix = [zeros(n, n) Matrix{Float64}(I, n, n);
              -(structural.M \ structural.K) -(structural.M \ structural.C)]
    solution = eigen(matrix)
    indices = findall(value -> imag(value) > 1.0e-6, solution.values)
    frequencies = imag.(solution.values[indices]) ./ (2pi)
    order = sortperm(frequencies)
    return frequencies[order], solution.vectors[:, indices[order]]
end

println("\n=== Chang regression ===")
chang = chang_case()
chang_general = build_structural_model(chang)
chang_legacy = WingPropellerAeroelastic.ChangAeroelastic.assemble_structural_model(
    chang.legacy_parameters,
)
for name in (:M, :C, :K)
    general_matrix = getproperty(chang_general, name)
    legacy_matrix = getproperty(chang_legacy, name)
    relative_error = norm(general_matrix - legacy_matrix) / max(norm(legacy_matrix), 1.0)
    @printf("%s relative error: %.3e\n", name, relative_error)
    @assert relative_error <= 1.0e-14
end
println("coupled frequencies [Hz]: ", structural_natural_frequencies(chang_general))

println("\n=== Xu reference validation ===")
xu = xu_case()
xu_wing = assemble_wing_propeller_structure(
    xu.wing,
    NamedTuple[],
    xu.operating_condition,
)
xu_coupled = build_structural_model(xu)
xu_wing_frequencies = structural_natural_frequencies(xu_wing; count = 5)
xu_coupled_frequencies = structural_natural_frequencies(xu_coupled; count = 7)
println("wing frequencies [Hz]: ", xu_wing_frequencies)
println("coupled frequencies at zero RPM [Hz]: ", xu_coupled_frequencies)
@assert maximum(abs.(xu_wing_frequencies .- [2.88, 16.68, 18.06, 46.59, 50.57])) <= 0.015
@assert maximum(abs.(xu_coupled_frequencies .- [2.85, 5.48, 7.00, 17.75, 19.48, 49.25, 51.37])) <= 0.04

println("\n=== Bohnisch active-subspace validation ===")
bohnisch = bohnisch_case()
bohnisch_general = build_structural_model(bohnisch)
bohnisch_reference = reference_bohnisch_matrices(bohnisch)
active_full_dofs = reduce(vcat, (
    collect(wing_node_dofs(node))[[3, 5, 4]]
    for node in 1:bohnisch_general.node_count
))
append!(active_full_dofs, (6bohnisch_general.node_count + 1):(
    size(bohnisch_general.M_global, 1)
))
for name in (:M_global, :C_global, :K_global)
    general_matrix = getproperty(bohnisch_general, name)[active_full_dofs, active_full_dofs]
    reference_matrix = getproperty(bohnisch_reference, name)
    maximum_error = maximum(abs.(general_matrix - reference_matrix))
    relative_error = norm(general_matrix - reference_matrix) /
        max(norm(reference_matrix), 1.0)
    @printf("%s active maximum error: %.3e (relative %.3e)\n",
        name, maximum_error, relative_error)
    @assert relative_error <= 1.0e-12
end

println("\n=== Bohnisch inactive-stiffness convergence ===")
multipliers = (1.0e2, 1.0e3, 1.0e4, 1.0e5)
rows = NamedTuple[]
for multiplier in multipliers
    case = bohnisch_case(inactive_stiffness_multiplier = multiplier)
    structural = build_structural_model(case)
    frequencies, vectors = state_space_modes(structural)
    n = structural.ndof_free
    active_reduced_dofs = reduce(vcat, (
        6(node - 2) .+ [3, 5, 4]
        for node in 2:structural.node_count
    ))
    append!(active_reduced_dofs, (structural.ndof_wing_free + 1):n)
    participation = [
        norm(vectors[active_reduced_dofs, mode]) / norm(vectors[1:n, mode])
        for mode in axes(vectors, 2)
    ]
    selected = [
        mode for mode in eachindex(frequencies)
        if frequencies[mode] > 1.0 && participation[mode] > 0.7
    ]
    active_frequencies = frequencies[selected[1:5]]
    push!(rows, (;
        multiplier,
        condition_K = cond(structural.K),
        active_frequencies,
    ))
    @printf("%7.0e  cond(K)=%10.3e  active Hz=%s\n",
        multiplier, cond(structural.K), string(round.(active_frequencies; digits = 6)))
end
reference_frequencies = rows[end].active_frequencies
for row in rows
    relative_change = maximum(abs.(row.active_frequencies ./ reference_frequencies .- 1))
    @printf("ratio=%7.0e  max change from 1e5=%.3e\n", row.multiplier, relative_change)
end
selected_row = rows[2] # 1e3
@assert maximum(abs.(selected_row.active_frequencies ./ reference_frequencies .- 1)) < 1.0e-4
println("Selected inactive stiffness multiplier: 1e3")
