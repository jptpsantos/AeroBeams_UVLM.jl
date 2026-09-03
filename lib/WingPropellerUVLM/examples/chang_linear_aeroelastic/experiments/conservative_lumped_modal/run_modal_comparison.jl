using LinearAlgebra
using Dates
using Printf

include(joinpath(@__DIR__, "ConservativeLumpedModal.jl"))
using .ConservativeLumpedModal

const OUTPUT_DIR = joinpath(@__DIR__, "output")
mkpath(OUTPUT_DIR)

function matching_frequency(modes, reference)
    index = findfirst(
        mode -> mode.family == reference.family &&
                mode.family_order == reference.family_order,
        modes,
    )
    return isnothing(index) ? NaN : modes[index].frequency_hz
end

cases = [
    (label = "source_16", ne = 16, mesh = :source, remap = :source_lumps),
    (label = "uniform_16", ne = 16, mesh = :uniform, remap = :conservative_spatial),
    (label = "uniform_20", ne = 20, mesh = :uniform, remap = :conservative_spatial),
    (label = "uniform_24", ne = 24, mesh = :uniform, remap = :conservative_spatial),
    (label = "uniform_30", ne = 30, mesh = :uniform, remap = :conservative_spatial),
    (label = "uniform_32", ne = 32, mesh = :uniform, remap = :conservative_spatial),
    (label = "uniform_36", ne = 36, mesh = :uniform, remap = :conservative_spatial),
    (label = "uniform_40", ne = 40, mesh = :uniform, remap = :conservative_spatial),
    (label = "uniform_48", ne = 48, mesh = :uniform, remap = :conservative_spatial),
    (label = "uniform_60", ne = 60, mesh = :uniform, remap = :conservative_spatial),
    (label = "uniform_80", ne = 80, mesh = :uniform, remap = :conservative_spatial),
    (label = "uniform_120", ne = 120, mesh = :uniform, remap = :conservative_spatial),
    (label = "uniform_160", ne = 160, mesh = :uniform, remap = :conservative_spatial),
    (label = "uniform_240", ne = 240, mesh = :uniform, remap = :conservative_spatial),
    (label = "uniform_320", ne = 320, mesh = :uniform, remap = :conservative_spatial),
]

results = NamedTuple[]
for case in cases
    println("[modal] $(case.label): assembling conservative lumped model")
    wing = assemble_wing_model(case.ne; mesh = case.mesh, remap = case.remap)
    audit = remap_audit(wing)
    wing_modes = modal_analysis(wing; number_of_modes = 16, classify_wing = true)
    coupled = assemble_coupled_model(wing)
    coupled_modes = modal_analysis(coupled; number_of_modes = 12)
    push!(results, (; case..., wing, audit, wing_modes, coupled, coupled_modes))
    @printf(
        "        min eig(M)=%+.6e, remap error=%.3e, wing f1=%.5f Hz, coupled f1=%.5f Hz\n",
        audit.mass_minimum_eigenvalue,
        audit.relative_spatial_inertia_error,
        wing_modes[1].frequency_hz,
        coupled_modes[1].frequency_hz,
    )
end

# Compare the interpretation requested by the user with the former code and
# with the alternative (but workbook-inconsistent) spanwise-element reading.
remap_comparisons = NamedTuple[]
for method in (:legacy_linear, :conservative_spatial, :element_overlap)
    model = assemble_wing_model(30; remap = method)
    audit = remap_audit(model)
    modes = modal_analysis(model; number_of_modes = 16, classify_wing = true)
    push!(remap_comparisons, (; method, model, audit, modes))
end

# Reproduce both the formerly working interpolation and the current failing
# signed-total rescaling at the requested 20-element mesh.
diagnostic_methods = [:legacy_linear, :signed_scaled]
diagnostics = NamedTuple[]
for method in diagnostic_methods
    println("[diagnostic] uniform_20 with $method")
    model = assemble_wing_model(20; remap = method)
    audit = remap_audit(model)
    modal_status = "not_run"
    first_frequency = NaN
    try
        modes = modal_analysis(model; number_of_modes = 6, classify_wing = true)
        modal_status = "success"
        first_frequency = modes[1].frequency_hz
    catch error
        modal_status = string(nameof(typeof(error)))
    end
    push!(diagnostics, (; method, model, audit, modal_status, first_frequency))
    @printf(
        "             min eig(M)=%+.6e, min eig(Jcg)=%+.6e, max |r|=%.6f m, modal=%s\n",
        audit.mass_minimum_eigenvalue,
        audit.minimum_cg_inertia_eigenvalue,
        audit.maximum_offset,
        modal_status,
    )
end

frequency_path = joinpath(OUTPUT_DIR, "modal_frequency_comparison.csv")
open(frequency_path, "w") do io
    println(io, "case,mesh,elements,system,mode_index,family,family_order,frequency_hz")
    for result in results
        for mode in result.wing_modes
            println(io, join((
                result.label, result.mesh, result.ne, "wing", mode.index,
                mode.family, mode.family_order, mode.frequency_hz,
            ), ','))
        end
        for mode in result.coupled_modes
            println(io, join((
                result.label, result.mesh, result.ne, "coupled", mode.index,
                "ordered", mode.index, mode.frequency_hz,
            ), ','))
        end
    end
end

audit_path = joinpath(OUTPUT_DIR, "mass_matrix_audit.csv")
open(audit_path, "w") do io
    println(io, "case,mesh,elements,mass_kg,relative_spatial_inertia_error,min_block_eigenvalue,min_mass_eigenvalue,min_stiffness_eigenvalue")
    for result in results
        a = result.audit
        println(io, join((
            result.label, result.mesh, result.ne, a.mass_kg,
            a.relative_spatial_inertia_error, a.block_minimum_eigenvalue,
            a.mass_minimum_eigenvalue, a.stiffness_minimum_eigenvalue,
        ), ','))
    end
end

diagnostic_path = joinpath(OUTPUT_DIR, "current_resampler_diagnostics.csv")
open(diagnostic_path, "w") do io
    println(io, "method,elements,mass_kg,relative_spatial_inertia_error,min_Jcg_eigenvalue,max_cg_offset_m,min_mass_eigenvalue,modal_status,first_frequency_hz")
    for result in diagnostics
        a = result.audit
        println(io, join((
            result.method, 20, a.mass_kg, a.relative_spatial_inertia_error,
            a.minimum_cg_inertia_eigenvalue, a.maximum_offset,
            a.mass_minimum_eigenvalue, result.modal_status,
            result.first_frequency,
        ), ','))
    end
end

uniform_results = filter(result -> result.mesh == :uniform, results)
finest_uniform = uniform_results[end]
convergence_path = joinpath(OUTPUT_DIR, "structural_mesh_convergence.csv")
open(convergence_path, "w") do io
    println(io, "system,elements,mode,frequency_hz,error_vs_320_percent,successive_change_percent,published_error_percent")
    previous_wing = nothing
    previous_coupled = nothing
    reference_wing = [
        matching_frequency(finest_uniform.wing_modes, ref) for ref in REFERENCE_WING_MODES
    ]
    reference_coupled = [mode.frequency_hz for mode in finest_uniform.coupled_modes[1:10]]
    for result in uniform_results
        current_wing = [
            matching_frequency(result.wing_modes, ref) for ref in REFERENCE_WING_MODES
        ]
        current_coupled = [mode.frequency_hz for mode in result.coupled_modes[1:10]]
        for i in eachindex(REFERENCE_WING_MODES)
            successive = isnothing(previous_wing) ? NaN :
                100 * (current_wing[i] / previous_wing[i] - 1)
            published = 100 * (
                current_wing[i] / REFERENCE_WING_MODES[i].frequency_hz - 1
            )
            println(io, join((
                "wing", result.ne, REFERENCE_WING_MODES[i].label,
                current_wing[i], 100 * (current_wing[i] / reference_wing[i] - 1),
                successive, published,
            ), ','))
        end
        for i in eachindex(current_coupled)
            successive = isnothing(previous_coupled) ? NaN :
                100 * (current_coupled[i] / previous_coupled[i] - 1)
            println(io, join((
                "coupled", result.ne, "f$i", current_coupled[i],
                100 * (current_coupled[i] / reference_coupled[i] - 1),
                successive, NaN,
            ), ','))
        end
        previous_wing = current_wing
        previous_coupled = current_coupled
    end
end

remap_comparison_path = joinpath(OUTPUT_DIR, "remap_method_comparison_30.csv")
open(remap_comparison_path, "w") do io
    println(io, "method,mass_kg,spatial_inertia_error,min_mass_eigenvalue,mode,frequency_hz,published_error_percent")
    for result in remap_comparisons
        for reference in REFERENCE_WING_MODES
            frequency = matching_frequency(result.modes, reference)
            error_percent = 100 * (frequency / reference.frequency_hz - 1)
            println(io, join((
                result.method, result.audit.mass_kg,
                result.audit.relative_spatial_inertia_error,
                result.audit.mass_minimum_eigenvalue, reference.label,
                frequency, error_percent,
            ), ','))
        end
    end
end

report_path = joinpath(OUTPUT_DIR, "modal_analysis_report.md")
validation_summary_path = joinpath(
    OUTPUT_DIR,
    "aeroelastic_20_trim_impulse",
    "conservative_lumped_20_trim_impulse_summary.txt",
)
validation_values = Dict{String,String}()
if isfile(validation_summary_path)
    for line in readlines(validation_summary_path)
        occursin(" = ", line) || continue
        key, value = split(line, " = "; limit = 2)
        validation_values[key] = value
    end
end
timestamp = Dates.format(now(), dateformat"yyyy-mm-dd HH:MM:SS")
open(report_path, "w") do io
    println(io, "# Chang wing-propeller lumped-inertia remesh audit")
    println(io)
    println(io, "Generated: $timestamp")
    println(io)
    println(io, "## Implemented interpretation")
    println(io)
    println(io, "The spreadsheet stiffness is treated as a distributed spanwise constitutive field. The 16 mass/CG/inertia records, plus the zero root record, define nodal reference values. Each complete 6-by-6 nodal spatial-inertia block is divided by its source tributary length, the resulting density is linearly interpolated, and that density is integrated exactly over every target element. Each target integral is lumped at its outer node. No free structural node is massless.")
    println(io)
    println(io, "The integral of each source-node hat function equals its tributary length. Therefore the transfer uses nonnegative weights and conserves the sum of the complete spatial-inertia matrix exactly. It conserves mass, first moments, and beam-axis inertia simultaneously, without separately rescaling signed components.")
    println(io)
    println(io, "Workbook connectivity is important: inertia elements 52-67 connect each non-root flexible wing station to a short offset node at the same span coordinate. They are rigid mass elements, not the flexible spanwise beam elements 20-35. Thus source_16 is the literal concentrated-inertia reference. The uniform cases intentionally construct an equivalent continuous distribution from those nodal reference records.")
    println(io)
    println(io, "## Mass-matrix verification")
    println(io)
    println(io, "| Case | Elements | Mass (kg) | Spatial-inertia error | min eig(M) | min nodal-block eig |")
    println(io, "|---|---:|---:|---:|---:|---:|")
    for result in results
        a = result.audit
        @printf(io, "| %s | %d | %.8f | %.3e | %+.6e | %+.6e |\n",
            result.label, result.ne, a.mass_kg,
            a.relative_spatial_inertia_error, a.mass_minimum_eigenvalue,
            a.block_minimum_eigenvalue)
    end
    println(io)
    println(io, "## Thirty-element remap interpretation sensitivity")
    println(io)
    println(io, "| Method | Mass (kg) | Spatial-inertia error | min eig(M) | OOP1 | IP1 | OOP2 | T1 | OOP3 |")
    println(io, "|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    for result in remap_comparisons
        frequencies = [matching_frequency(result.modes, ref) for ref in REFERENCE_WING_MODES]
        @printf(io, "| %s | %.8f | %.3e | %+.6e | %.5f | %.5f | %.5f | %.5f | %.5f |\n",
            result.method, result.audit.mass_kg,
            result.audit.relative_spatial_inertia_error,
            result.audit.mass_minimum_eigenvalue, frequencies...)
    end
    println(io)
    println(io, "legacy_linear is the former point-sampled nodal-density interpolation. conservative_spatial is the corrected exact-integration version. element_overlap assumes the records belong to the preceding spanwise flexible elements; it is included only to quantify that alternative interpretation, which conflicts with the workbook connectivity.")
    println(io)
    println(io, "## Isolated-wing natural frequencies")
    println(io)
    println(io, "| Case | OOP1 | IP1 | OOP2 | T1 | OOP3 |")
    println(io, "|---|---:|---:|---:|---:|---:|")
    for result in results
        values = [matching_frequency(result.wing_modes, ref) for ref in REFERENCE_WING_MODES]
        @printf(io, "| %s | %.5f | %.5f | %.5f | %.5f | %.5f |\n", result.label, values...)
    end
    @printf(io, "| Published reference | %.5f | %.5f | %.5f | %.5f | %.5f |\n",
        (ref.frequency_hz for ref in REFERENCE_WING_MODES)...)
    println(io)
    finest = results[end]
    finest_wing = [matching_frequency(finest.wing_modes, ref) for ref in REFERENCE_WING_MODES]
    finest_errors = [
        100 * (finest_wing[i] / REFERENCE_WING_MODES[i].frequency_hz - 1)
        for i in eachindex(REFERENCE_WING_MODES)
    ]
    source_wing = [matching_frequency(results[1].wing_modes, ref) for ref in REFERENCE_WING_MODES]
    source_errors = [
        100 * (source_wing[i] / REFERENCE_WING_MODES[i].frequency_hz - 1)
        for i in eachindex(REFERENCE_WING_MODES)
    ]
    uniform_20 = only(filter(result -> result.label == "uniform_20", results))
    uniform_20_wing = [
        matching_frequency(uniform_20.wing_modes, ref) for ref in REFERENCE_WING_MODES
    ]
    uniform_20_errors = [
        100 * (uniform_20_wing[i] / REFERENCE_WING_MODES[i].frequency_hz - 1)
        for i in eachindex(REFERENCE_WING_MODES)
    ]
    println(io, "The exact 16-record/source-station model differs from the published values by " * join(
        [@sprintf("%s %+.3f%%", REFERENCE_WING_MODES[i].label, source_errors[i]) for i in eachindex(source_errors)],
        ", ",
    ) * ".")
    println(io)
    println(io, "The valid 20-element uniform remesh is not yet frequency-converged: its differences from the published values are " * join(
        [@sprintf("%s %+.3f%%", REFERENCE_WING_MODES[i].label, uniform_20_errors[i]) for i in eachindex(uniform_20_errors)],
        ", ",
    ) * ".")
    println(io)
    println(io, "At the finest uniform mesh ($(finest.ne) elements), the differences from the published isolated-wing values are " * join(
        [@sprintf("%s %+.3f%%", REFERENCE_WING_MODES[i].label, finest_errors[i]) for i in eachindex(finest_errors)],
        ", ",
    ) * ".")
    println(io)
    println(io, "## Coupled wing-pylon-propeller natural frequencies")
    println(io)
    println(io, "The table gives the first 10 undamped, non-gyroscopic structural frequencies. The pitch/yaw attachment is kept at the exact span fraction 0.83 by work-conjugate interpolation between adjacent beam nodes.")
    println(io, "This intentionally removes attachment-node snapping from the modal remesh comparison. The unmodified production driver still selects the nearest structural node; in the 20-element end-to-end test it attached at y = 6.375 m instead of the exact 6.225 m position, so that separate source of mesh dependence remains outside the inertia-remap correction.")
    println(io)
    println(io, "| Case | f1 | f2 | f3 | f4 | f5 | f6 | f7 | f8 | f9 | f10 |")
    println(io, "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for result in results
        values = [mode.frequency_hz for mode in result.coupled_modes[1:10]]
        @printf(io, "| %s | %.5f | %.5f | %.5f | %.5f | %.5f | %.5f | %.5f | %.5f | %.5f | %.5f |\n",
            result.label, values...)
    end
    println(io)
    source_coupled = [mode.frequency_hz for mode in results[1].coupled_modes[1:10]]
    finest_coupled = [mode.frequency_hz for mode in finest.coupled_modes[1:10]]
    coupled_errors = 100 .* (finest_coupled ./ source_coupled .- 1)
    println(io, "Relative to the exact 16-record/source-station model, the finest uniform-mesh changes in the first ten coupled frequencies are: " * join(
        [@sprintf("f%d %+.3f%%", i, coupled_errors[i]) for i in eachindex(coupled_errors)],
        ", ",
    ) * ".")
    println(io)
    println(io, "## Structural mesh convergence relative to the finest uniform mesh")
    println(io)
    println(io, "| Elements | Maximum wing-mode error | Maximum coupled-mode error |")
    println(io, "|---:|---:|---:|")
    reference_wing = [matching_frequency(finest.wing_modes, ref) for ref in REFERENCE_WING_MODES]
    reference_coupled = [mode.frequency_hz for mode in finest.coupled_modes[1:10]]
    for result in uniform_results
        current_wing = [matching_frequency(result.wing_modes, ref) for ref in REFERENCE_WING_MODES]
        current_coupled = [mode.frequency_hz for mode in result.coupled_modes[1:10]]
        wing_error = maximum(abs.(100 .* (current_wing ./ reference_wing .- 1)))
        coupled_error = maximum(abs.(100 .* (current_coupled ./ reference_coupled .- 1)))
        @printf(io, "| %d | %.4f%% | %.4f%% |\n", result.ne, wing_error, coupled_error)
    end
    println(io)
    println(io, "## Why the present 20-element approach fails")
    println(io)
    println(io, "| Remap | min eig(M) | min eig(recovered J_CG) | max recovered CG offset (m) | Modal solve |")
    println(io, "|---|---:|---:|---:|---|")
    for result in diagnostics
        a = result.audit
        @printf(io, "| %s | %+.6e | %+.6e | %.6f | %s |\n",
            result.method, a.mass_minimum_eigenvalue,
            a.minimum_cg_inertia_eigenvalue, a.maximum_offset,
            result.modal_status)
    end
    println(io)
    println(io, "The signed-scaled method globally rescales each first moment and product of inertia by its own ratio of source total to interpolated total. Several of these totals are small because positive and negative entries cancel. At 20 elements, the independent scale factors amplify interpolation error, produce unrealistic CG offsets, and make recovered center-of-mass inertia tensors and the global mass matrix indefinite. The failure therefore occurs before Generalized-alpha time integration; it is not caused by the GA parameters.")
    println(io)
    println(io, "The old unscaled componentwise interpolation can remain positive definite, but it does not conserve the complete spatial inertia. The conservative spatial-block remap satisfies both requirements.")
    if !isempty(validation_values)
        println(io)
        println(io, "## End-to-end 20-element aeroelastic validation")
        println(io)
        println(io, "A separate in-memory launcher injected the conservative arrays into the actual production UVLM/Generalized-alpha driver. The short trim-plus-impulse test is a solver smoke test, not a damping-convergence result.")
        println(io)
        validation_fields = (
            ("integrated steps", "integrated_steps"),
            ("finite history", "history_is_finite"),
            ("all coupling steps converged", "all_coupling_steps_converged"),
            ("maximum coupling iterations", "max_coupling_iterations"),
            ("maximum state residual", "max_coupling_state_residual"),
            ("maximum load residual", "max_coupling_load_residual"),
            ("maximum coupled-equilibrium residual", "max_coupling_equilibrium_residual"),
            ("validation passed", "validation_passed"),
        )
        for (label, key) in validation_fields
            println(io, "- $label: " * validation_values[key])
        end
    end
end

# Assertions turn this script into an executable regression test.
@assert all(result.audit.mass_minimum_eigenvalue > 0 for result in results)
@assert all(result.audit.block_minimum_eigenvalue > 0 for result in results)
@assert all(result.audit.relative_spatial_inertia_error < 1e-12 for result in results)
signed = only(filter(result -> result.method == :signed_scaled, diagnostics))
@assert signed.audit.mass_minimum_eigenvalue < 0

println("[modal] wrote $frequency_path")
println("[modal] wrote $audit_path")
println("[modal] wrote $diagnostic_path")
println("[modal] wrote $convergence_path")
println("[modal] wrote $remap_comparison_path")
println("[modal] wrote $report_path")
