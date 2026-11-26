# =============================================================================
# Numerical / conceptual verification of section 1.5.3.2 ("Free-flight stability analyses") of theory/main.tex
# =============================================================================
# Companion of verify_chapter1.jl, verify_C.jl and verify_2_*.jl (same conventions: check(...) kinds pass / typo / note).
#
# What is checked (equation labels of main.tex)
#   A   continuum: the global balance laws of a beam (eq:rigid_body_motion, corrected) follow from the mixed variational statement
#       eq:hamilton_lagrange6 with a rigid-body virtual motion (δu = δa + δθ~ (r+u), δψ = δθ). Checked with manufactured smooth
#       fields that satisfy the strong form of eq:hamilton_lagrange6 exactly (loads f, m defined as the residuals), integrated by
#       Gauss quadrature, with time and space derivatives by forward-mode AD. The version AS PRINTED before the red \fix{} (no
#       v_A~ ∫P term, no ∫(r+u)~ f moment) must fail.
#   B   finite element model of AeroBeams: in a dynamic simulation of a FREE assembly (two beams joined at an interior node, CG
#       offsets, point inertia, initially curved beam, basis A with v_A ≠ 0 and ω_A ≠ 0, load pulse), the momenta computed from the
#       states and their rates satisfy the discrete global balance laws at every time step: exactly for the linear momentum, exactly
#       for the angular momentum with straight elements, and up to O(Δℓ²) for initially curved elements. That is, the balance
#       laws are implied by the FE equations of a free assembly (they are NOT additional equations).
#   C   free-flight eigenanalysis WITHOUT attachment springs (code): the rigid-body translations relative to A span the null space of
#       the Jacobian, the shifted eigenproblem gives the same eigenvalues as a dense generalized eigensolver, and the attachment
#       springs distort the trim and the flight dynamic modes.
#
# Result categories (last column of the report)
#   ok                         agreement
#   ERRO NO TEXTO CONFIRMADO   the formula evaluated exactly AS PRINTED in main.tex (before the red fix) disagrees
#   NOTA                       statement/approximation that is not exact (see report)
#
# Run:   julia --project=<AeroBeams root> theory/verification/verify_1_5_3.jl      (~2 min)
# =============================================================================

using LinearAlgebra, ForwardDiff, QuadGK, Random, Printf, SparseArrays
using AeroBeams

const TOL = 1e-9

# ----------------------------------------------------------------------------- notation
tilde(v) = [0 -v[3] v[2]; v[3] 0 -v[1]; -v[2] v[1] 0]
const a1 = [1.0, 0.0, 0.0]
function R_exp(p)                      # rotation vector -> rotation matrix
    φ = sqrt(sum(abs2, p) + 1e-300)
    N = tilde(p / φ)
    Matrix{eltype(p)}(I, 3, 3) + sin(φ) * N + (1 - cos(φ)) * N * N
end

# ----------------------------------------------------------------------------- bookkeeping
struct Res
    sec::String
    label::String
    desc::String
    err::Float64
    kind::Symbol
    ok::Bool
end
const RESULTS = Res[]
function check(sec, label, desc, a, b; kind=:pass, tol=TOL)
    err = maximum(abs.(a .- b)) / max(1.0, maximum(abs.(a)), maximum(abs.(b)))
    ok = kind == :pass ? err < tol : err > 1e-6
    push!(RESULTS, Res(sec, label, desc, err, kind, ok))
    return ok
end
function checkb(sec, label, desc, cond::Bool; kind=:pass, err=NaN)
    ok = kind == :pass ? cond : !cond
    push!(RESULTS, Res(sec, label, desc, isnan(err) ? (cond ? 0.0 : 1.0) : err, kind, ok))
    return ok
end
function mutant(sec, label, desc, truth, wrong)
    err = maximum(abs.(truth .- wrong)) / max(1.0, maximum(abs.(truth)), maximum(abs.(wrong)))
    push!(RESULTS, Res(sec, label, "[mutation] " * desc, err, :pass, err > 1e-6))
end

# =============================================================================
# A  Continuum: global balance laws from the variational statement
# =============================================================================
function section_continuum()
    s = "1.5.3-A"
    rng = MersenneTwister(1532)
    ℓ = 1.7
    gx, gw = QuadGK.gauss(Float64, 48)
    xq = ℓ / 2 .* (gx .+ 1); wq = ℓ / 2 .* gw
    t0 = 0.37

    # Undeformed geometry: Rz(x) = Rz0 exp(x k~) (constant curvature), r' = Rz a1 (arclength parametrization), r(0) = r0
    Rz0 = R_exp(0.8 * randn(rng, 3)); k = 0.6 * randn(rng, 3); r0 = randn(rng, 3)
    Rz(x) = Rz0 * R_exp(x * k)
    rpos(x) = r0 + (x / 2) * sum(w * (Rz(x / 2 * (z + 1)) * a1) for (z, w) in zip(QuadGK.gauss(Float64, 24)...))

    # Motion of basis A (components in A): any values (only enter algebraically at time t)
    ω = 0.7 * randn(rng, 3); vA = 2.0 * randn(rng, 3)

    # Smooth random fields (x, t)
    C = [randn(rng, 3) for _ in 1:20]
    u(x, t) = 0.3 * (C[1] * sin(1.1x + 2t) + C[2] * x^2 * cos(t) + C[3] * t * x)
    p(x, t) = 0.5 * (C[4] * sin(1.3x + t) + C[5] * x * t + C[6])
    free(x) = sin(π * x / ℓ)                                  # vanishes at both ends (free ends)
    Fs(x, t; fr=false) = (fr ? free(x) : 1.0) * (C[7] * cos(x + t) + C[8] * x + C[9])
    Ms(x, t; fr=false) = (fr ? free(x) : 1.0) * (C[10] * sin(2x - t) + C[11])
    Ps(x, t) = C[12] * cos(0.7x + 1.5t) + C[13] * x
    Hs(x, t) = C[14] * sin(0.5x - t) + C[15] * x * t

    RR(x, t) = R_exp(p(x, t)) * Rz(x)
    X(x, t) = rpos(x) + u(x, t)                               # position of B relative to A
    udot(x, t) = ForwardDiff.derivative(τ -> u(x, τ), t)
    Xp(x, t) = Rz(x) * a1 + ForwardDiff.derivative(ξ -> u(ξ, t), x)
    # kinematic relations enforced: γ* from the compatibility (δF) equation, V* from the velocity (δP) equation
    γs(x, t) = RR(x, t)' * Xp(x, t) - a1
    Vs(x, t) = RR(x, t)' * (vA + tilde(ω) * X(x, t) + udot(x, t))
    # quantities in A
    PA(x, t) = RR(x, t) * Ps(x, t)
    HA(x, t) = RR(x, t) * Hs(x, t)
    FA(x, t; fr=false) = RR(x, t) * Fs(x, t; fr=fr)
    MA(x, t; fr=false) = RR(x, t) * Ms(x, t; fr=fr)
    dt(fun, x, t) = ForwardDiff.derivative(τ -> fun(x, τ), t)
    dx(fun, x, t) = ForwardDiff.derivative(ξ -> fun(ξ, t), x)
    # manufactured loads: strong form of eq:hamilton_lagrange6 (equations of δu and δψ) satisfied exactly
    f(x, t; fr=false) = dt(PA, x, t) + tilde(ω) * PA(x, t) - dx((ξ, τ) -> FA(ξ, τ; fr=fr), x, t)
    m(x, t; fr=false) = dt(HA, x, t) + tilde(ω) * HA(x, t) + RR(x, t) * (tilde(Vs(x, t)) * Ps(x, t) - tilde(a1 + γs(x, t)) * Fs(x, t; fr=fr)) - dx((ξ, τ) -> MA(ξ, τ; fr=fr), x, t)

    ∫(g) = sum(w * g(x) for (x, w) in zip(xq, wq))

    for fr in (false, true)
        tag = fr ? " (free ends)" : ""
        L = ∫(x -> PA(x, t0)); Ldot = ∫(x -> dt(PA, x, t0))
        HAt(x, t) = HA(x, t) + tilde(X(x, t)) * PA(x, t)
        Hs_A = ∫(x -> HAt(x, t0)); Hs_Adot = ∫(x -> dt(HAt, x, t0))
        Fext = ∫(x -> f(x, t0; fr=fr)) + FA(ℓ, t0; fr=fr) - FA(0.0, t0; fr=fr)
        Mext = ∫(x -> m(x, t0; fr=fr) + tilde(X(x, t0)) * f(x, t0; fr=fr)) +
               (MA(ℓ, t0; fr=fr) + tilde(X(ℓ, t0)) * FA(ℓ, t0; fr=fr)) - (MA(0.0, t0; fr=fr) + tilde(X(0.0, t0)) * FA(0.0, t0; fr=fr))
        check(s, "linear momentum" * tag, "(d/dt + ω~)∫P dx1 = ∫f dx1 + [F̂]_0^ℓ   (δa: δu = δa, δψ = 0)", Ldot + tilde(ω) * L, Fext)
        check(s, "angular momentum" * tag, "(d/dt + ω~)∫(H + (r+u)~P) dx1 + v_A~ ∫P dx1 = ∫(m + (r+u)~f) dx1 + [M̂ + (r+u)~F̂]_0^ℓ   (δθ)", Hs_Adot + tilde(ω) * Hs_A + tilde(vA) * L, Mext)
        if fr
            # text AS PRINTED before the fix: (d/dt + ω~)∫(H + (r+u)~P) = ∫ m
            check(s, "eq:rigid_body_motion", "AS PRINTED: (d/dt + ω~)∫(H + (r+u)~P) dx1 = ∫m dx1  (no v_A~∫P, no moment of f)", Hs_Adot + tilde(ω) * Hs_A, ∫(x -> m(x, t0; fr=fr)); kind=:typo)
            mutant(s, "no v_A~∫P", "angular balance without the term v_A~ ∫P must disagree", Hs_Adot + tilde(ω) * Hs_A + tilde(vA) * L, Hs_Adot + tilde(ω) * Hs_A)
            mutant(s, "no (r+u)~f", "angular balance without the moment of the distributed forces must disagree", Mext, ∫(x -> m(x, t0; fr=fr)))
        end
    end

    # Pointwise steps of the derivation (random data)
    for i in 1:5
        R = R_exp(randn(rng, 3)); R0 = R_exp(randn(rng, 3)); γ = 0.1 * randn(rng, 3); Fst = randn(rng, 3)
        Xprime = R * R0 * (a1 + γ)
        δθ = randn(rng, 3)
        # δu' = δθ~ X'  ->  δu'^T (R R0 F*) - δθ^T R R0 (a1~ + γ~) F* = 0
        check(s, "cancellation of F", "δu'ᵀ(RR₀F*) + δψᵀ[-RR₀(ã₁+γ̃)F*] = 0 for δu' = δθ~X', δψ = δθ, X' = RR₀(a₁+γ*)  (eq. of δF)", dot(tilde(δθ) * Xprime, R * R0 * Fst) - dot(δθ, R * R0 * tilde(a1 + γ) * Fst), 0.0)
        Xv, ud, P, Pd, w, v = randn(rng, 3), randn(rng, 3), randn(rng, 3), randn(rng, 3), randn(rng, 3), randn(rng, 3)
        V = v + tilde(w) * Xv + ud                              # R R0 V* (eq. of δP), with v_b = v_A + ω~ r
        check(s, "identity", "X~(Ṗ + ω~P) + V~P = (X~P)˙ + ω~X~P + v_A~P  with V = v_A + ω~X + u̇, Ẋ = u̇", tilde(Xv) * (Pd + tilde(w) * P) + tilde(V) * P, tilde(ud) * P + tilde(Xv) * Pd + tilde(w) * tilde(Xv) * P + tilde(v) * P)
    end
    return nothing
end

# =============================================================================
# B  Finite element model: discrete balance laws implied by the FE equations of a free assembly
# =============================================================================
const τp = 0.15
const F0 = [3.0, -2.0, 4.0]
const M0 = [0.5, 1.0, -0.8]
pulse(t) = t <= τp ? sin(π * t / τp)^2 : 0.0

function free_assembly(; n2=6, curved=true, v_A=[2.0, -1.0, 0.5], ω_A=[0.3, -0.2, 0.4], p_A0=[0.1, -0.2, 0.15])
    S = isotropic_stiffness_matrix(EA=2e4, GAy=1e4, GAz=1e4, GJ=2e2, EIy=3e2, EIz=5e2)
    Ib = inertia_matrix(ρA=1.0, ρIy=0.02, ρIz=0.05, e2=0.03, e3=-0.02)
    b1 = create_Beam(name="b1", length=2.0, nElements=10, S=[S], I=[Ib], p0=transform_to_WM([0.2, 0.1, 0.3], from=:E321))
    b2 = create_Beam(name="b2", length=1.5, nElements=n2, S=[S], I=[Ib], k=(curved ? [0.05, 0.3, -0.2] : zeros(3)), p0=transform_to_WM([1.0, -0.3, 0.2], from=:E321))
    add_point_inertia_to_beam!(beam=b2, elementID=n2 - 1, pointInertia=create_PointInertia(mass=0.4, η=[0.05, 0.1, -0.03], Ixx=0.01, Iyy=0.02, Izz=0.015))
    bc = create_BC(name="pulse", beam=b2, node=n2 + 1, conditions=Dict(:F1A => t -> F0[1] * pulse(t), :F2A => t -> F0[2] * pulse(t), :F3A => t -> F0[3] * pulse(t), :M1A => t -> M0[1] * pulse(t), :M2A => t -> M0[2] * pulse(t), :M3A => t -> M0[3] * pulse(t)))
    model = create_Model(name="freeAssembly", beams=[b1, b2], connections=[(b2 => b1, 1 => 6)], BCs=[bc], v_A=v_A, ω_A=ω_A, p_A0=p_A0)
    return model, b2
end

# Discrete momenta (basis A) and their rates from the states and states' rates of the elements
function discrete_momenta(model, st)
    L = zeros(3); HA = zeros(3); Ldot = zeros(3); HAdot = zeros(3)
    for (e, element) in enumerate(model.assemblyData.elements)
        Δℓ, R0, r = element.geometry.Δℓ, element.geometry.R0, element.assemblyData.r
        mid, rt = st[e].midpoint, st[e].rates
        u, p, P, H = Vector(mid.u), Vector(mid.p), Vector(mid.P), Vector(mid.H)
        udot, pdot, Pdot, Hdot = Vector(rt.udot), Vector(rt.pdot), Vector(rt.Pdot), Vector(rt.Hdot)
        R = rotation_tensor_WM(p)
        Rdot = ForwardDiff.derivative(x -> rotation_tensor_WM(p + x * pdot), 0.0)
        PA, HA_ = R * R0 * P, R * R0 * H
        PAdot, HAdot_ = Rdot * R0 * P + R * R0 * Pdot, Rdot * R0 * H + R * R0 * Hdot
        X = r + u
        L += Δℓ * PA; HA += Δℓ * (HA_ + tilde(X) * PA)
        Ldot += Δℓ * PAdot; HAdot += Δℓ * (HAdot_ + tilde(udot) * PA + tilde(X) * PAdot)
    end
    return L, HA, Ldot, HAdot
end

# Maximum residuals of the global balance laws over the time steps (the initial one is skipped: its rates are extrapolated)
function balance_residuals(; tf=0.4, Δt=1e-3, kw...)
    model, b2 = free_assembly(; kw...)
    problem = create_DynamicProblem(model=model, finalTime=tf, Δt=Δt, systemSolver=create_NewtonRaphson(absTol=1e-12, relTol=1e-6, maxIter=50), displayProgress=false)
    solve!(problem)
    sol = problem.solutionHistory
    kin = model.kinematicsBasisA
    tipElem = b2.elements[end]
    out = Dict(k => 0.0 for k in (:L, :H, :Htext, :HnoV, :HnoXF, :scaleL, :scaleH))
    for n in 2:length(sol.t)
        st = sol.states[n]
        L, HA, Ldot, HAdot = discrete_momenta(model, st)
        RA = kin.R_AVec[n]; ω = RA' * kin.ω_AVec[n]; vA = RA' * kin.v_AVec[n]
        F = F0 * pulse(sol.t[n]); M = M0 * pulse(sol.t[n])
        Xtip = tipElem.assemblyData.r_n2 + Vector(st[tipElem.assemblyData.globalID].nodal.u2)
        lhsH = HAdot + tilde(ω) * HA
        out[:L] = max(out[:L], norm(Ldot + tilde(ω) * L - F))
        out[:H] = max(out[:H], norm(lhsH + tilde(vA) * L - (M + tilde(Xtip) * F)))
        out[:Htext] = max(out[:Htext], norm(lhsH - M))
        out[:HnoV] = max(out[:HnoV], norm(lhsH - (M + tilde(Xtip) * F)))
        out[:HnoXF] = max(out[:HnoXF], norm(lhsH + tilde(vA) * L - M))
        out[:scaleL] = max(out[:scaleL], norm(Ldot), norm(tilde(ω) * L), norm(F))
        out[:scaleH] = max(out[:scaleH], norm(HAdot), norm(tilde(ω) * HA), norm(tilde(vA) * L))
    end
    return out, model, b2
end

function section_fem()
    s = "1.5.3-B"
    rs, _, _ = balance_residuals(curved=false)
    check(s, "linear (straight)", "Σ_e Δℓ_e[(R_eR₀_eP*_e)˙ + ω~R_eR₀_eP*_e] = F̂_ext   at every time step (free assembly, v_A, ω_A ≠ 0)", rs[:L] / rs[:scaleL], 0.0)
    check(s, "angular (straight)", "Σ_e Δℓ_e[(d/dt + ω~)(H + X~P)_e + v_A~P_e] = M̂_ext + X~F̂_ext,  X_e = r_e + u_e", rs[:H] / rs[:scaleH], 0.0; tol=1e-8)
    check(s, "text (FE)", "AS PRINTED (no v_A~L, no moment of the nodal force) applied to the FE solution", rs[:Htext] / rs[:scaleH], 0.0; kind=:typo)
    checkb(s, "mutation v_A~L", "[mutation] without v_A~L the FE solution violates the balance", rs[:HnoV] / rs[:scaleH] > 1e-3; err=rs[:HnoV] / rs[:scaleH])
    checkb(s, "mutation X~F", "[mutation] without the moment X~F̂ of the nodal force the FE solution violates the balance", rs[:HnoXF] / rs[:scaleH] > 1e-3; err=rs[:HnoXF] / rs[:scaleH])
    # initially curved elements: linear momentum still exact, angular momentum up to O(Δℓ²)
    rc6, model6, b26 = balance_residuals(curved=true, n2=6)
    rc12, _, _ = balance_residuals(curved=true, n2=12)
    check(s, "linear (curved)", "linear momentum balance, initially curved beam: still exact (internal forces cancel element by element)", rc6[:L] / rc6[:scaleL], 0.0)
    check(s, "angular (curved)", "angular momentum balance, initially curved beam (6 elements): NOT exact, since r_e (arc midpoint) ≠ (r_n1+r_n2)/2 and r_n2 - r_n1 ≠ Δℓ R₀_e a₁", rc6[:H] / rc6[:scaleH], 0.0; kind=:note)
    ratio = rc6[:H] / rc12[:H]
    checkb(s, "O(Δℓ²)", @sprintf("curved: error(6 elem)/error(12 elem) = %.2f ≈ 4 (second order in Δℓ)", ratio), 3.5 < ratio < 4.5; err=abs(ratio - 4) / 4)
    dmid = maximum(norm(e.assemblyData.r - (e.assemblyData.r_n1 + e.assemblyData.r_n2) / 2) for e in b26.elements)
    dchord = maximum(norm(e.assemblyData.r_n2 - e.assemblyData.r_n1 - e.geometry.Δℓ * e.geometry.R0[:, 1]) for e in b26.elements)
    checkb(s, "curved geometry", @sprintf("curved elements: |r_e - (r_n1+r_n2)/2| = %.1e (≈ kΔℓ²/8), |r_n2 - r_n1 - Δℓ R₀_e a₁| = %.1e (≈ k²Δℓ³/24)", dmid, dchord), dmid > dchord; kind=:pass, err=dmid)
    return nothing
end

# =============================================================================
# C  Free-flight eigenanalysis without attachment springs (code)
# =============================================================================
function section_eigen()
    s = "1.5.3-C"
    U = 30.0
    kw = (aeroSolver=Indicial(), airspeed=U, nElemWing=20, nElemTailBoom=5, nElemHorzStabilizer=4, stiffnessFactor=1.0, stabilizersAero=true, includeVS=true, wingCd0=1e-2, stabsCd0=1e-2)
    out = Dict{Bool,Any}()
    for springs in (false, true)
        NR = create_NewtonRaphson(trimSolverMethod=:LM, displayStatus=false)
        trimModel, _ = create_cHALE(; kw..., δElevIsTrimVariable=true, thrustIsTrimVariable=true, withAttachmentSprings=springs, μ=1e-1)
        trimProblem = create_TrimProblem(model=trimModel, systemSolver=NR)
        solve!(trimProblem)
        T, δ = trimProblem.trimSolution.trimLoads[1], trimProblem.trimSolution.trimδ[1]
        eigenModel, _ = create_cHALE(; kw..., δElev=δ, thrust=T, withAttachmentSprings=springs, μ=1e-1)
        eig = create_EigenProblem(model=eigenModel, eigenmethod=:arnoldi, nev=100, refTrimProblem=trimProblem)
        solve!(eig)
        out[springs] = (T=T, δ=δ, eig=eig, model=eigenModel)
    end
    eig, model = out[false].eig, out[false].model
    J, M = eig.systemArrays.jacobian, eig.systemArrays.inertia
    N = AeroBeams.neutral_translations(model, J)
    checkb(s, "null space", "free flight: the 3 rigid-body translations relative to A are in the null space of the Jacobian (‖J N‖ = 0)", size(N, 2) == 3 && norm(J * N) <= 1e-12 * opnorm(J, 1); err=norm(J * N))
    Nspring = AeroBeams.neutral_translations(out[true].model, out[true].eig.systemArrays.jacobian)
    checkb(s, "springs", "with attachment springs, the translations are not neutral", size(Nspring, 2) == 0)
    # shifted (sparse) eigensolution vs dense generalized eigensolution
    λd = filter(isfinite, eigen(Matrix(J), -Matrix(M)).values)
    λs = eig.eigensolution.λ
    near0 = count(x -> abs(x) < 1e-5, λd)
    checkb(s, "neutral modes", "dense generalized eigenproblem: ≥ 3 null eigenvalues (translations; heading)", near0 >= 3; err=Float64(near0))
    osc = sort(filter(x -> imag(x) > 1e-3 && abs(x) < 15, λd), by=abs)
    err = maximum(minimum(abs.(λs .- x)) / abs(x) for x in osc)
    check(s, "shift-invert", @sprintf("flight dynamic and elastic modes |λ| < 15 (%d pairs) from the shifted sparse eigensolver = dense generalized eigenvalues", length(osc)), err, 0.0; tol=1e-6)
    # attachment springs: trim and modes
    check(s, "springs: trim", @sprintf("trim thrust with springs (μ = 0.1) %.3f N vs. free flight %.3f N", out[true].T, out[false].T), out[true].T, out[false].T; kind=:note)
    lowosc(e, k) = sort(filter(x -> imag(x) > 1e-3, e.eigensolution.λ), by=abs)[k]
    fmt(x) = @sprintf("%.4f%+.4fi", real(x), imag(x))
    check(s, "springs: modes", "two lowest oscillatory modes with springs " * fmt(lowosc(out[true].eig, 1)) * ", " * fmt(lowosc(out[true].eig, 2)) * " vs. free flight " * fmt(lowosc(out[false].eig, 1)) * ", " * fmt(lowosc(out[false].eig, 2)) * " (spurious mode / distorted phugoid)", [reim(lowosc(out[true].eig, 1))...], [reim(lowosc(out[false].eig, 1))...]; kind=:note)
    return nothing
end

# =============================================================================
println("="^100)
println("Verification of section 1.5.3.2 of theory/main.tex   -   tolerance = $TOL")
println("="^100)
section_continuum()
section_fem()
section_eigen()

# ----------------------------------------------------------------------------- report
println()
let lastsec = ""
    for r in RESULTS
        mark = r.kind == :pass ? (r.ok ? "ok" : "FALHOU") :
               r.kind == :typo ? (r.ok ? "ERRO NO TEXTO CONFIRMADO" : "!! esperava-se discrepancia e nao houve") :
                                 (r.ok ? "NOTA (ver relatorio)" : "!! esperava-se diferenca e nao houve")
        if r.sec != lastsec
            println("\n------------------------- ", r.sec, " -------------------------")
            lastsec = r.sec
        end
        @printf("%-8s %-22s %-150s %9.2e  %s\n", r.sec, first(r.label, 22), first(r.desc, 150), r.err, mark)
    end
end
nok = count(r -> r.kind == :pass && r.ok, RESULTS)
nbad = count(r -> r.kind == :pass && !r.ok, RESULTS)
ntypo = count(r -> r.kind == :typo && r.ok, RESULTS)
nnote = count(r -> r.kind == :note && r.ok, RESULTS)
nmiss = count(r -> r.kind != :pass && !r.ok, RESULTS)
println("\nSUMMARY: ", nok, " checks passed, ", nbad, " FAILED, ", ntypo, " errors-in-text confirmed (formula as printed disagrees), ", nnote, " notes, ", nmiss, " expected discrepancies that did not occur")
