# Standalone kernel verification and diagnostic examples; no production run.
# The final NaN prints deliberately exercise unsupported singular inputs.
# The geometry diagnostic uses the 20x5 wing / 5x5 blade case inspected on
# 2026-09-10; the blade angles are converted to the grid's radian convention.
# Upstream comparison formula: VortexLattice.jl commit
# b89ecd089b3fac5737378afeecc73a47ed1006d7, src/induced.jl.
# See FINITE_CORE_SELECTION.md for interpretation and the run command.
using WingPropellerUVLM, StaticArrays, LinearAlgebra, Test, Printf
const segment = WingPropellerUVLM.bound_induced_velocity
const trailing = WingPropellerUVLM.trailing_induced_velocity
function upstream_segment(r1, r2, rc)
    a, b, q = norm(r1), norm(r2), dot(r1, r2)
    cross(r1, r2) / (a^2*b^2-q^2+rc^2*(a^2+b^2-2a*b)) *
        ((a^2-q)/sqrt(a^2+rc^2)+(b^2-q)/sqrt(b^2+rc^2)) / (4pi)
end
function midpoint_reference(r1, r2, rc; n=100000)
    line = r1-r2
    result = zero(r1)
    for i in 1:n
        r = r1-((i-0.5)/n)*line
        result += cross(line,r)/(dot(r,r)+rc^2)^(3/2)
    end
    result/(4pi*n)
end
@testset "Independent core audit" begin
    axis = normalize(SVector(1.0,2.0,-1.0))
    rot = 2axis*axis'-I
    worst_error = 0.0
    for (r1,r2) in ((SVector(0.5,0.0,0.0),SVector(-0.5,0.0,0.0)),
        (SVector(0.0,0.0,0.0),SVector(-1.0,0.0,0.0)),
        (SVector(0.5,0.01,0.02),SVector(-0.5,0.01,0.02)),
        (SVector(0.02,0.04,0.01),SVector(-0.98,0.04,0.01)),
        (SVector(2.0,0.02,0.05),SVector(1.0,0.02,0.05)),
        (SVector(0.3,0.7,-0.4),SVector(-0.2,0.4,0.9)))
        for rc in (0.01,0.1,0.3)
            v = segment(r1,r2,true,rc)
            ref = midpoint_reference(r1,r2,rc)
            @test all(isfinite,v)
            @test isapprox(v,ref;rtol=2e-7,atol=1e-12)
            @test isapprox(segment(rot*r1,rot*r2,true,rc),rot*v;rtol=1e-10,atol=1e-12)
            @test isapprox(segment(r2,r1,true,rc),-v;rtol=1e-12,atol=1e-12)
            @test isapprox(segment(100r1,100r2,true,100rc),v/100;rtol=1e-10,atol=1e-12)
            worst_error = max(worst_error,norm(v-ref)/max(norm(ref),1e-12))
        end
    end
    @printf("Maximum segment quadrature relative error: %.3e\n",worst_error)
    for direction in (SVector(1.0,0.0,0.0),axis)
        for p in (SVector(0.1,0.03,-0.02),SVector(-0.1,0.3,0.04))
            for rc in (0.01,0.1,0.3)
                @test isapprox(trailing(p,direction,true,rc),segment(p,p-1e5direction,true,rc);rtol=1e-8)
            end
        end
    end
    # Independent circular-line integral at ring center: Gamma R^2/[2(R^2+rc^2)^(3/2)].
    for rc in (0.01,0.1,0.3)
        previous_error = Inf
        for n in (16,64,256)
            velocity = zero(SVector{3,Float64})
            for i in 1:n
                a,b = 2pi*(i-1)/n,2pi*i/n
                velocity += segment(-SVector(cos(a),sin(a),0.0),-SVector(cos(b),sin(b),0.0),true,rc)
            end
            reference = SVector(0.0,0.0,1/(2*(1+rc^2)^(3/2)))
            err = norm(velocity-reference)/norm(reference)
            @test err < previous_error
            previous_error = err
        end
        @test previous_error < 6e-5
        @printf("Ring center at rc=%.3f: 256-segment relative error %.3e\n",rc,previous_error)
    end
end
println("Bisector comparison: unit segment, Gamma=1, radius=0.1")
for h in (0.01,0.001,0.0001)
    r1,r2 = SVector(0.5,h,0.0),SVector(-0.5,h,0.0)
    @printf("h=%.4g, local speed=%.8g, upstream speed=%.8g\n",h,norm(segment(r1,r2,true,0.1)),norm(upstream_segment(r1,r2,0.1)))
end
for (label,r1,r2,fc,rc) in (("zero-core midpoint",SVector(0.5,0.0,0.0),SVector(-0.5,0.0,0.0),true,0.0),
    ("inviscid midpoint",SVector(0.5,0.0,0.0),SVector(-0.5,0.0,0.0),false,0.0),
    ("inviscid endpoint",SVector(0.0,0.0,0.0),SVector(-1.0,0.0,0.0),false,0.0))
    println(label," => ",segment(r1,r2,fc,rc))
end
_,_,twists = get_nodal_properties_chang(5)
grid = generate_propeller_blades_grid(1.15,0.197,5,5,deg2rad.(twists .- 90.0),4)[1]
_,_,panels = grid_to_surface_panels(grid;fcore=(c,ds)->0.1ds)
cores = getproperty.(panels,:core_size)
@printf("Blade core range: %.9f .. %.9f m\n",minimum(cores),maximum(cores))
worst_core_jump = 0.0
worst_residual = 0.0
for j in axes(panels,2), i in 1:size(panels,1)-1
    a,b = panels[i,j],panels[i+1,j]
    @assert a.rbl == b.rtl && a.rbr == b.rtr
    rstart,rend = a.rbl,a.rbr
    line = rend-rstart
    normal = normalize(cross(line,SVector(1.0,0.0,0.0)))
    p = (rstart+rend)/2+0.001normal
    va = segment(p-rstart,p-rend,true,a.core_size)
    vb = segment(p-rend,p-rstart,true,b.core_size)
    global worst_core_jump = max(worst_core_jump,abs(a.core_size-b.core_size)/min(a.core_size,b.core_size))
    global worst_residual = max(worst_residual,norm(va+vb)/max(norm(va),norm(vb)))
end
@printf("Chord-adjacent shared-edge core jump: %.5f%%; unit equal-circulation cancellation residual at 1mm: %.5f%%\n",100worst_core_jump,100worst_residual)


