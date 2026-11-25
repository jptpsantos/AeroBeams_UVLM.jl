# Frequency-domain unsteady vortex-lattice method (VLM) for the generalized aerodynamic forces (GAFs) of the Pazy wing in small harmonic motion about its undeformed (planar, zero lift) state, and p-k flutter solution in the basis of the vibration modes of AeroBeams. Axes: x downstream, y spanwise (projected, root at y=0, which is a symmetry plane with an image wing), z up (lift direction). The leading edge of the root section is at x=0, and the reference axis (spar) is at x = xs0 + y tan(Λ)
using LinearAlgebra

## Vortex lattice kernels (all vortices and points lie in the plane z=0)
# Normal (z) velocity induced at (px,py) by a unit-strength straight vortex segment from (ax,ay) to (bx,by), with a small core to avoid singularities (relative to the segment length)
@inline function segmentInducedVelocity(px, py, ax, ay, bx, by, ε²)
    r1x = px - ax; r1y = py - ay
    r2x = px - bx; r2y = py - by
    r0x = bx - ax; r0y = by - ay
    cz = r1x * r2y - r1y * r2x
    l1 = sqrt(r1x^2 + r1y^2); l2 = sqrt(r2x^2 + r2y^2)
    (l1 < 1e-12 || l2 < 1e-12) && return 0.0
    d = cz^2 + ε² * (r0x^2 + r0y^2)
    d < 1e-30 && return 0.0
    return (r0x * (r1x / l1 - r2x / l2) + r0y * (r1y / l1 - r2y / l2)) * cz / (4π * d)
end

# Normal velocity at (px,py) induced by a unit-strength vortex ring with corners P (front at y1, front at y2, rear at y2, rear at y1; the front segment is directed along +y, so that a positive circulation gives positive lift), and by its image across the plane y=0
@inline function ringInducedVelocity(px, py, P, ε²)
    w = 0.0
    for q in 1:4
        a = P[q]; b = P[mod1(q+1, 4)]
        w += segmentInducedVelocity(px, py, a[1], a[2], b[1], b[2], ε²)
    end
    Q = ((P[2][1], -P[2][2]), (P[1][1], -P[1][2]), (P[4][1], -P[4][2]), (P[3][1], -P[3][2]))
    for q in 1:4
        a = Q[q]; b = Q[mod1(q+1, 4)]
        w += segmentInducedVelocity(px, py, a[1], a[2], b[1], b[2], ε²)
    end
    return w
end

## Vortex lattice: M chordwise and N spanwise panels of equal size, bound vortices at the quarter-chord and collocation points at the three-quarter-chord of each panel, and a planar wake of Nw rows of rings (convected with the freestream, one panel length per row) shed from the trailing edge
struct VortexLattice
    c::Float64          # streamwise chord
    Lp::Float64         # projected semispan
    Λ::Float64          # sweep angle (aft positive)
    xs0::Float64        # chordwise position of the reference axis at the root, from the leading edge
    M::Int
    N::Int
    Δx::Float64
    Δy::Float64
    xcp::Matrix{Float64}    # collocation points
    yc::Vector{Float64}     # centers of the strips
    xbound::Matrix{Float64} # midpoints of the bound vortices
    xcenter::Matrix{Float64}# centers of the panels
    Ab::Matrix{Float64}     # influence of the bound rings on the collocation points
    Aw::Array{Float64,3}    # influence of the wake rings (Nc × N × Nw)
    Nw::Int
end

function create_vortex_lattice(; c, Lp, Λ=0.0, xs0, M=8, N=30, wakeChords=20.0, core=1e-6)
    Δx = c/M; Δy = Lp/N; tanΛ = tan(Λ)
    Nw = ceil(Int, wakeChords*M)
    xLE(y) = y*tanΛ
    yc = [(j-0.5)*Δy for j in 1:N]
    xcp = [xLE(yc[j]) + (i-1)*Δx + 0.75Δx for i in 1:M, j in 1:N]
    xbound = [xLE(yc[j]) + (i-1)*Δx + 0.25Δx for i in 1:M, j in 1:N]
    xcenter = [xLE(yc[j]) + (i-1)*Δx + 0.5Δx for i in 1:M, j in 1:N]
    ε² = (core/Δx)^2
    Nc = M*N
    idx(i,j) = (j-1)*M + i
    Ab = zeros(Nc, Nc)
    for jj in 1:N
        y1 = (jj-1)*Δy; y2 = jj*Δy
        for ii in 1:M
            xf(y) = xLE(y) + (ii-1)*Δx + 0.25Δx
            xr(y) = xLE(y) + ii*Δx + 0.25Δx
            P = ((xf(y1),y1), (xf(y2),y2), (xr(y2),y2), (xr(y1),y1))
            for j in 1:N, i in 1:M
                Ab[idx(i,j), idx(ii,jj)] = ringInducedVelocity(xcp[i,j], yc[j], P, ε²)
            end
        end
    end
    Aw = zeros(Nc, N, Nw)
    for m in 1:Nw, jj in 1:N
        y1 = (jj-1)*Δy; y2 = jj*Δy
        xf(y) = xLE(y) + c + 0.25Δx + (m-1)*Δx
        xr(y) = xLE(y) + c + 0.25Δx + m*Δx
        P = ((xf(y1),y1), (xf(y2),y2), (xr(y2),y2), (xr(y1),y1))
        for j in 1:N, i in 1:M
            Aw[idx(i,j), jj, m] = ringInducedVelocity(xcp[i,j], yc[j], P, ε²)
        end
    end
    return VortexLattice(c, Lp, Λ, xs0, M, N, Δx, Δy, xcp, yc, xbound, xcenter, Ab, Aw, Nw)
end

# Influence matrix at reduced frequency k = ωb/U (b = c/2): the wake row m carries the circulation of the trailing-edge row m time steps before (phase exp(-i m ω Δx/U))
function influence_matrix(L::VortexLattice, k)
    M, N = L.M, L.N
    φ = (k/(L.c/2))*L.Δx
    A = complex.(L.Ab)
    wake = zeros(ComplexF64, M*N, N)
    for m in 1:L.Nw
        @views wake .+= cis(-φ*m) .* L.Aw[:, :, m]
    end
    for j in 1:N
        A[:, (j-1)*M + M] .+= wake[:, j]
    end
    return A
end

# Steady lift slope of the planform, CL/(2π α), for a uniform incidence (the area-weighted mean of the normalized lift slope)
function steady_lift_slope(L::VortexLattice; k=1e-4)
    Γ = influence_matrix(L, k) \ fill(-1.0 + 0im, L.M*L.N)
    ΓTE = [real(Γ[(j-1)*L.M + L.M]) for j in 1:L.N]
    return 2*sum(ΓTE)*L.Δy/(L.c*L.Lp)/(2π)
end

## Vibration modes: bending displacement h (along the lift direction) and twist θ (nose up positive) of the reference axis as functions of its arclength s
function cubic_spline_with_derivative(x, y)
    n = length(x); h = diff(x)
    A = zeros(n, n); r = zeros(n)
    A[1,1] = 1; A[n,n] = 1
    for i in 2:n-1
        A[i,i-1] = h[i-1]; A[i,i] = 2(h[i-1]+h[i]); A[i,i+1] = h[i]
        r[i] = 3*((y[i+1]-y[i])/h[i] - (y[i]-y[i-1])/h[i-1])
    end
    cc = A \ r
    b = [(y[i+1]-y[i])/h[i] - h[i]*(2cc[i]+cc[i+1])/3 for i in 1:n-1]
    d = [(cc[i+1]-cc[i])/(3h[i]) for i in 1:n-1]
    f = s -> (i = clamp(searchsortedlast(x, s), 1, n-1); Δ = s - x[i]; y[i] + b[i]*Δ + cc[i]*Δ^2 + d[i]*Δ^3)
    df = s -> (i = clamp(searchsortedlast(x, s), 1, n-1); Δ = s - x[i]; b[i] + 2cc[i]*Δ + 3d[i]*Δ^2)
    return f, df
end

struct ModeField
    h::Function
    dh::Function
    θ::Function
    dθ::Function
end
function mode_fields(modes, s)
    map(modes) do m
        h, dh = cubic_spline_with_derivative(s, m.h)
        θ, dθ = cubic_spline_with_derivative(s, m.θ)
        ModeField(h, dh, θ, dθ)
    end
end

# Vertical displacement z and its streamwise slope ∂z/∂x of a mode at the point (x,y) of the planform: the point is displaced as the perpendicular to the (swept) reference axis through it (at arclength s), by h(s) + θ(s) d, where d is its distance to the axis toward the leading edge
function mode_displacement(L::VortexLattice, mode::ModeField, x, y, axisLength)
    sinΛ, cosΛ = sin(L.Λ), cos(L.Λ)
    s = clamp((x - L.xs0)*sinΛ + y*cosΛ, 0.0, axisLength)
    d = -(x - L.xs0)*cosΛ + y*sinΛ
    return mode.h(s) + mode.θ(s)*d, sinΛ*(mode.dh(s) + mode.dθ(s)*d) - mode.θ(s)*cosΛ
end

# Generalized aerodynamic forces per ρU² at reduced frequency k: Q[s,r] is the generalized force in the mode s due to a unit harmonic motion of the mode r. The normal velocity at the collocation points is iωz + U ∂z/∂x, and the force on each panel is ρ[U (Γ - Γ_prev) Δy + ∂Γ/∂t Δx Δy] (Kutta-Joukowski and unsteady term), applied at the bound vortex and at the center of the panel, respectively
function generalized_forces(L::VortexLattice, modes::Vector{ModeField}, axisLength, k)
    M, N = L.M, L.N
    ω = k/(L.c/2)
    R = length(modes)
    ν = zeros(ComplexF64, M*N, R)
    zBound = zeros(R, M, N); zCenter = zeros(R, M, N)
    for r in 1:R, j in 1:N, i in 1:M
        z, dzdx = mode_displacement(L, modes[r], L.xcp[i,j], L.yc[j], axisLength)
        ν[(j-1)*M + i, r] = im*ω*z + dzdx
        zBound[r,i,j] = mode_displacement(L, modes[r], L.xbound[i,j], L.yc[j], axisLength)[1]
        zCenter[r,i,j] = mode_displacement(L, modes[r], L.xcenter[i,j], L.yc[j], axisLength)[1]
    end
    Γ = influence_matrix(L, k) \ ν
    Q = zeros(ComplexF64, R, R)
    for r in 1:R, j in 1:N, i in 1:M
        g = Γ[(j-1)*M + i, r]
        gPrev = i > 1 ? Γ[(j-1)*M + i - 1, r] : zero(g)
        Fb = (g - gPrev)*L.Δy
        Fc = im*ω*g*L.Δx*L.Δy
        for s in 1:R
            Q[s,r] += Fb*zBound[s,i,j] + Fc*zCenter[s,i,j]
        end
    end
    return Q
end

## Strip theory (Theodorsen's function from Jones' approximation of Wagner's function, as in the Indicial solver of AeroBeams) of a straight wing, with the tip loss function W at the strips, to compare with the vortex lattice. With finiteWingLag, the amplitudes of the exponentials of the lag are scaled by κ = (1-(1-ΣA)/W) / ΣA (not less than zero)
function strip_generalized_forces(modes, s, c, Lp, xs, k, W; xac=0.25, β=1.0, A=(0.165, 0.335), bC=(0.0455, 0.3), finiteWingLag=false)
    N = length(W); b = c/2; ω = k/b; Δy = Lp/N; a = (xs - c/2)/b
    yc = [(j-0.5)*Δy for j in 1:N]
    R = length(modes)
    h = [cubic_spline_with_derivative(s, m.h)[1] for m in modes]
    θ = [cubic_spline_with_derivative(s, m.θ)[1] for m in modes]
    Q = zeros(ComplexF64, R, R)
    for j in 1:N
        ΣA = sum(A)
        κ = finiteWingLag ? max((1 - (1-ΣA)/W[j])/ΣA, 0.0) : 1.0
        C = 1 - κ*sum(A[n]*im*k/(im*k + bC[n]) for n in 1:2)
        for r in 1:R
            hr, θr = h[r](yc[j]), θ[r](yc[j])
            hd = -im*ω*hr; hdd = ω^2*hr; θd = im*ω*θr; θdd = -ω^2*θr    # Theodorsen's convention: displacement h positive downward
            circulatory = hd + θr + b*(0.5 - a)*θd
            Lc = 2π*b*(C/β)*W[j]*circulatory
            Ln = π*b^2*(hdd + θd - b*a*θdd)
            Mc = Lc*(xs - xac*c)
            Mn = π*b^2*(b*a*hdd - b*(0.5 - a)*θd - b^2*(1/8 + a^2)*θdd)
            for q in 1:R
                Q[q,r] += Δy*((Ln + Lc)*h[q](yc[j]) + (Mn + Mc)*θ[q](yc[j]))
            end
        end
    end
    return Q
end

## p-k flutter solution: roots p = σ + iω of [p² M + K - ρU² Q(k)] η = 0, with k = Im(p) b/U, for each airspeed, following each mode (the roots at an airspeed start from those of the previous one)
function pk_flutter(M, K, Qfun, b, ρ, URange)
    R = length(M)
    p = [im*sqrt(K[r]/M[r]) for r in 1:R]
    vecs = [begin v = zeros(ComplexF64, R); v[r] = 1.0; v end for r in 1:R]
    roots = zeros(ComplexF64, length(URange), R)
    for (iU,U) in enumerate(URange), r in 1:R
        pr = p[r]; vr = vecs[r]
        for _ in 1:80
            k = max(imag(pr), 1e-3)*b/U
            F = eigen(Diagonal(1 ./ M)*(Diagonal(K) - ρ*U^2*Qfun(k)))
            candidates = [im*sqrt(μ) for μ in F.values]
            scores = [abs(dot(vr, F.vectors[:,j]/norm(F.vectors[:,j]))) - 0.02*abs(candidates[j]-pr)/max(abs(pr),1.0) - (imag(candidates[j]) < 0 ? 10 : 0) for j in eachindex(candidates)]
            j = argmax(scores)
            pNew = candidates[j]; vecs[r] = F.vectors[:,j]/norm(F.vectors[:,j])
            converged = abs(pNew - pr) < 1e-9*max(abs(pr), 1.0)
            pr = converged ? pNew : 0.5*pr + 0.5*pNew
            converged && break
        end
        p[r] = pr
        roots[iU,r] = pr
    end
    return roots
end

# Airspeeds at which the real part of a branch crosses zero upward (flutter onset) and downward (offset), by linear interpolation
function stability_crossings(URange, branch; upward=true)
    crossings = Float64[]
    for i in 1:length(URange)-1
        σ1, σ2 = real(branch[i]), real(branch[i+1])
        if upward ? (σ1 < 0 && σ2 >= 0) : (σ1 >= 0 && σ2 < 0)
            push!(crossings, URange[i] + (URange[i+1]-URange[i])*σ1/(σ1-σ2))
        end
    end
    return crossings
end
