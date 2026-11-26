# Generalized masses of the Nastran modes of the S0 Pazy wing without the tip rod (S0_norod_modal.op2, dumped by dump_norod_modes.py), with different mass representations, independently of the stiffness of the beam model:
#   M_true:  lumped GFEM masses with the displacements and rotations of their grids in the Nastran mode
#   M_fe:    lumped GFEM masses moving with the rigid body motion of the spar sections (fitted at the beam nodes) linearly interpolated along the span, i.e. the exact reduction of the GFEM to the kinematics of the beam
#   M_nodal: beam nodes' mass models (condense_to_beam), half of each node's mass to each adjacent element (as inertia_matrices_from_Nastran_data), with the kinematics of the element at its midpoint (V and Ω the averages of those of its nodes, as in the beam element)
#   M_elem:  mass model condensed to the midpoints of the elements (condense_to_elements), with the same kinematics
# Run with the BeamROM project: julia --project=<BeamROM> modal_mass_study.jl
using BeamROM, LinearAlgebra, DelimitedFiles, Printf
include(joinpath(pkgdir(BeamROM), "nastran", "workflow.jl"))

reference_axis = 0.4653
dir = @__DIR__
nmodes = 10

# Nastran data
G = readdlm(joinpath(dir, "S0_norod_modal_grids.csv"), ','; skipstart=1)
disp = Dict((Int(r[1]), Int(r[2])) => Float64.(r[3:8]) for r in eachrow(G))
S = readdlm(joinpath(dir, "S0_norod_modal_sections.csv"), ','; skipstart=1)
sections = Dict((Int(r[1]), round(r[2], digits=7)) => Float64.(r[3:8]) for r in eachrow(S))
sectionYs = sort(unique(round.(S[:, 2], digits=7)))

bdf = pazy_no_rod_bdf()
gfem = read_gfem(bdf)
lumped = lump_by_grid(lump_gfem(gfem))
grids = gfem.grids

# Kinetic energy (twice) of a rigid body of mass m with the center of mass at η from the point with the translation V and the rotation Ω, and the inertia tensor J about the center of mass
body_M(m, η, J, V, Ω) = m * norm(V + cross(Ω, η))^2 + dot(Ω, J * Ω)

function M_true(mode)
    sum(lumped) do lm
        d = disp[(mode, lm.grid)]
        body_M(lm.m, lm.p - grids[lm.grid], lm.J, d[1:3], d[4:6])
    end
end

# Section motion at the node k of the beam (the row of spar grids at its y)
function node_motion(mode, x)
    y = sectionYs[argmin(abs.(sectionYs .- x[2]))]
    abs(y - x[2]) < 1e-5 || error("no section of the spar at y = $(x[2])")
    return sections[(mode, y)]
end

function M_fe(mode, xref)
    mo = [node_motion(mode, x) for x in xref]
    ys = [x[2] for x in xref]
    sum(lumped) do lm
        y = grids[lm.grid][2]
        k = clamp(searchsortedlast(ys, y), 1, length(ys) - 1)
        s = clamp((y - ys[k]) / (ys[k+1] - ys[k]), 0, 1)
        f = (1 - s) * mo[k] + s * mo[k+1]
        ref = [xref[k][1], y, xref[k][3]]
        body_M(lm.m, lm.p - ref, lm.J, f[1:3], f[4:6])
    end
end

# Elements' kinematics at the midpoints (V, Ω averages of the nodal motions), and the midpoints
function element_kinematics(mode, xref)
    mo = [node_motion(mode, x) for x in xref]
    return [((mo[k][1:3] + mo[k+1][1:3]) / 2, (mo[k][4:6] + mo[k+1][4:6]) / 2) for k in 1:length(xref)-1], [(xref[k] + xref[k+1]) / 2 for k in 1:length(xref)-1]
end

function M_nodal(mode, nodal)
    xref = nodal.xref
    kin, xmid = element_kinematics(mode, xref)
    n = length(xmid)
    total = 0.0
    for e in 1:n
        V, Ω = kin[e]
        for (k, w) in ((e, e == 1 ? 1.0 : 0.5), (e + 1, e + 1 == n + 1 ? 1.0 : 0.5))
            total += w * body_M(nodal.m[k], xref[k] + nodal.offsets[k] - xmid[e], nodal.tensors[k], V, Ω) # the weights of the mass (a fraction w of the node's mass, and of its inertia)
        end
    end
    return total
end

function M_elem(mode, elem)
    kin, xmid = element_kinematics(mode, elem.xref)
    sum(body_M(elem.m[e], xmid[e] + elem.offsets[e] - xmid[e], elem.tensors[e], kin[e]...) for e in eachindex(xmid))
end

# Reference with the sections at all the rows of spar grids as beam nodes (no interpolation error of the beam mesh, only the rigid section hypothesis)
rowAxis = [[spar_reference_axis(gfem, :nominal; reference_axis, sweep=0)[1][1], y, spar_reference_axis(gfem, :nominal; reference_axis, sweep=0)[1][3]] for y in sectionYs]

for mesh in [:nominal, :refined]
    nodal = pazy_beam_mass_model(bdf; nodes=mesh, reference_axis, sweep=0, verbose=false)
    elem = pazy_beam_mass_model(bdf; nodes=mesh, reference_axis, sweep=0, verbose=false, to_elements=true)
    println("\nMesh :$mesh (generalized masses M = 2T; Nastran modes mass-normalized: M_true should be 1 for the wing alone)")
    @printf("%4s %9s | %9s %9s | %9s %9s | %8s %8s\n", "mode", "M_true", "rows/true", "fe/true", "nodal/true", "elem/true", "nodal/fe", "elem/fe")
    for mode in 1:nmodes
        Mt, Mf, Mn, Me = M_true(mode), M_fe(mode, nodal.xref), M_nodal(mode, nodal), M_elem(mode, elem)
        Mr = M_fe(mode, rowAxis)
        @printf("%4d %9.5f | %9.4f %9.4f | %9.4f %9.4f | %8.4f %8.4f\n", mode, Mt, Mr / Mt, Mf / Mt, Mn / Mt, Me / Mt, Mn / Mf, Me / Mf)
    end
end
