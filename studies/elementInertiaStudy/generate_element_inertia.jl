# Generates the equivalent beam mass model of the S0 Pazy wing without the tip rod condensed to the midpoints of the elements (BeamROM, condense_to_elements), and compares it with the nodal one (condense_to_beam), for the study of the element-wise inertia distribution
using BeamROM, LinearAlgebra, Printf
include(joinpath(pkgdir(BeamROM), "nastran", "workflow.jl"))

reference_axis = 0.4653
outdir = @__DIR__

for mesh in [:nominal, :refined]
    bdf = pazy_no_rod_bdf()
    nodal = pazy_beam_mass_model(bdf; nodes=mesh, reference_axis, sweep=0, verbose=false)
    elem = pazy_beam_mass_model(bdf; nodes=mesh, reference_axis, sweep=0, verbose=false, to_elements=true)
    @printf("%s: nodes %d, elements %d, mass nodal %.10f, element %.10f\n", mesh, length(nodal.m), length(elem.m), sum(nodal.m), sum(elem.m))
    # Element midpoints (y relative to the first node, as in the nodal files)
    xmid = [(elem.xref[k] + elem.xref[k+1]) / 2 for k in 1:length(elem.m)]
    xmid0 = [[x[1], x[2] - elem.xref[1][2], x[3]] for x in xmid]
    write_beam_mass_csv(joinpath(outdir, "S0_norod_inertia_matrices_$(mesh)mesh_elements.csv"), xmid0, elem.m, elem.offsets, elem.tensors)
    # Check against the stored nodal data
    ref = joinpath(pkgdir(BeamROM), "nastran", "Pazy_Technion", "S0_norod", "S0_norod_RA0p4653_inertia_matrices_$(mesh)mesh.csv")
    stored = readlines(ref)[2:end]
    mstored = [parse(Float64, split(l, ",")[5]) for l in stored]
    @printf("   nodal vs stored file: max |Δm| = %.2e\n", maximum(abs.(mstored - nodal.m)))
    # Total first moment and inertia about the root, nodal vs element
    cgy(m, off, pos) = sum(m[k] * (pos[k][2] + off[k][2]) for k in eachindex(m)) / sum(m)
    @printf("   span cg nodal %.6f, element %.6f\n", cgy(nodal.m, nodal.offsets, nodal.xref), cgy(elem.m, elem.offsets, xmid))
end
