# =============================================================================
# Numerical verification of Appendix C ("Finite element model implementation"), sections C.1 - C.5, of theory/main.tex
# =============================================================================
# Companion of verify_chapter1.jl (same conventions: check(...) kinds pass / typo / note, forward-mode AD, random data).
#
# Sections covered
#   C.1  element dependent variables   (Rz_e, Δℓ, ω_be, v_be)
#   C.2  nodal outputs                 (u_n, p_n, F_n^*, M_n^*)
#   C.3  loading input                 (dead/follower loads, aerodynamic basis change, gravity, R_A update)
#   C.4  point inertias and springs
#   C.5  Jacobian and inertia-matrix entries (every block of the elemental Jacobian/inertia matrix, nodal terms)
#
# Method for C.5 (the heart of the appendix)
#   The elemental residual G(x, xdot) = {Fu^<-, Fu^->, Fp^<-, Fp^->, FF^<-, FF^->, FM^<-, FM^->, FV, FO} of
#   eqs. (Fu)-(FH) is coded ONCE, directly from the discrete equations of section 1.5. The "true" Jacobian is
#   d/dx G(x, c*x + bar) (forward-mode AD: includes the Newmark dependence xdot = c x + bar, c = 2/dt) and the "true"
#   inertia matrix is d/dxdot G. Each 3x3 block is then compared with the closed-form expression printed in the text.
#   Blocks whose printed expression disagrees with AD are reported as "ERRO NO TEXTO CONFIRMADO" and the corrected
#   expression is checked right after (it must pass).
#
# Result categories (last column of the report)
#   ok                         agreement (step of the text is correct)
#   ERRO NO TEXTO CONFIRMADO   the formula evaluated exactly AS PRINTED in main.tex disagrees (typo); the corrected one passes
#   NOTA                       statement that is imprecise / not general / inconsistent with the implementation (see report)
#
# Run:   julia --project=<AeroBeams root> theory/verification/verify_C.jl      (~30 s)
# Last run (2026-11-14): 285 checks passed, 0 failed, 30 errors-in-text confirmed (each wrong Jacobian block counted once per random trial), 8 notes.
# =============================================================================

using LinearAlgebra, ForwardDiff, QuadGK, Random, Printf

const TOL = 1e-9

# ----------------------------------------------------------------------------- notation
tilde(v) = [0 -v[3] v[2]; v[3] 0 -v[1]; -v[2] v[1] 0]
axial(S) = [S[3, 2], S[1, 3], S[2, 1]]
const e1 = [1.0, 0.0, 0.0]
const I3 = Matrix(1.0I, 3, 3)

# Wiener-Milenkovic (scaled set, |p| < 4) and rotation vector
function R_WM(p)
    p0 = 2 - dot(p, p) / 8
    υ = 1 / (4 - p0)
    υ^2 * ((1 / υ^2 - 2 * dot(p, p)) * I + 2 * p * p' + 2 * p0 * tilde(p))
end
function R_exp(p)
    φ = norm(p)
    N = tilde(p / φ)
    I + sin(φ) * N + (1 - cos(φ)) * N * N
end
HT_WM(p) = (p0 = 2 - dot(p, p) / 8; υ = 1 / (4 - p0); 2υ^2 * (p0 * I + (1 / 4) * p * p' - tilde(p)))     # eq:HT_WM
HinvT_WM(p) = (p0 = 2 - dot(p, p) / 8; (1 / 2) * (p0 * I + (1 / 4) * p * p' + tilde(p)))                 # eq:HinvT_WM
# WM parameters of a rotation matrix (|phi| < pi), via Euler parameters
function p_from_R(R)
    e0 = sqrt(1 + tr(R)) / 2
    ev = axial(R - R') / (4 * e0)
    4 * ev / (1 + e0)
end
# Tangent operator by its meaning: axial(Rdot R') = H pdot
function H_num(Rfun, p)
    R = Rfun(p)
    J = ForwardDiff.jacobian(q -> vec(Rfun(q)), p)
    hcat((axial(reshape(J[:, j], 3, 3) * R') for j in 1:3)...)
end
# Third-order tensor times vector, eq:mul3:  [dT/dp1 v, dT/dp2 v, dT/dp3 v]
mul3(Tf, p, v) = hcat((ForwardDiff.derivative(s -> Tf(p .+ s .* (1:3 .== j)) * v, 0.0) for j in 1:3)...)

# ----------------------------------------------------------------------------- random smooth fields f(x,t)
struct Fld
    c::Matrix{Float64}
    a::Float64
    b::Float64
end
function Fld(seed; scale=1.0)
    rng = MersenneTwister(seed)
    Fld(scale .* (rand(rng, 3, 6) .- 0.5), 1 + rand(rng), 0.5 + rand(rng))
end
(f::Fld)(x, t) = f.c[:, 1] + f.c[:, 2] * x + f.c[:, 3] * t + f.c[:, 4] * (x * t) + f.c[:, 5] * sin(f.a * x + f.b * t) + f.c[:, 6] * x^2

# ----------------------------------------------------------------------------- bookkeeping
struct Res
    sec::String
    label::String
    desc::String
    err::Float64
    kind::Symbol      # :pass (must agree) | :typo (formula AS PRINTED must disagree) | :note (statement that must be false/imprecise)
    ok::Bool
end
const RESULTS = Res[]
function check(sec, label, desc, a, b; kind=:pass)
    err = maximum(abs.(a .- b)) / max(1.0, maximum(abs.(a)), maximum(abs.(b)))
    ok = kind == :pass ? err < TOL : err > 1e-6
    push!(RESULTS, Res(sec, label, desc, err, kind, ok))
    return ok
end
# boolean version: `cond` = the claim holds. :pass -> must hold; :typo/:note -> claim as printed must be false
function checkb(sec, label, desc, cond::Bool; kind=:pass)
    ok = kind == :pass ? cond : !cond
    push!(RESULTS, Res(sec, label, desc, cond ? 0.0 : 1.0, kind, ok))
    return ok
end

function gaussint(g, a, b; n=40)
    z, w = QuadGK.gauss(Float64, n)
    sum((b - a) / 2 * w[i] * g((a + b) / 2 + (b - a) / 2 * z[i]) for i in eachindex(w))
end

# =============================================================================
# C.1  Element dependent variables
# =============================================================================
function appendixC1()
    s = "C.1"
    C(label, desc, a, b; kw...) = check(s, label, desc, a, b; kw...)
    rng = MersenneTwister(101)
    kp = [0.3, -0.2, 0.5]; κ = norm(kp); kt = tilde(kp); kk = kp * kp'
    Rzz = R_exp(0.7 * randn(rng, 3))
    expk(x) = I + sin(κ * x) * kt / κ + (1 - cos(κ * x)) * kt * kt / κ^2          # Rodrigues (AD friendly)
    Rz(x) = Rzz * expk(x)

    # eq:Rzp: Rz' = Rz k~^+  (k^+ constant)  =>  Rz = Rz(0) e^{k~ x}
    C("eq:Rzp2", "Rz(x) = Rz(0) e^{k~ x}  solves  Rz' = Rz k~^+", ForwardDiff.derivative(Rz, 0.3), Rz(0.3) * kt)
    C("eq:etoktilpx1", "Rodrigues formula = matrix exponential e^{k~ x}", expk(0.4), exp(kt * 0.4))
    a, b = randn(rng, 3), randn(rng, 3)
    C("eq:etoktilpx1", "identity  a~ b~ = b a^T - (a^T b) I", tilde(a) * tilde(b), b * a' - dot(a, b) * I)
    C("eq:etoktilpx1", "third form: (I-kk^T/κ²)cos + k~/κ sin + kk^T/κ²", expk(0.4), (I - kk / κ^2) * cos(κ * 0.4) + kt / κ * sin(κ * 0.4) + kk / κ^2)
    # the integration constant: as printed, Rz = e^{k~ x} e^{c}; but d/dx(e^{k~x} e^c) = k~ Rz  !=  Rz k~
    Rz_printed(x) = expk(x) * Rzz
    C("eq:Rzp (as printed)", "Rz = e^{k~ x} e^{c} (constant on the RIGHT) solves Rz' = Rz k~ ?", ForwardDiff.derivative(Rz_printed, 0.3), Rz_printed(0.3) * kt; kind=:typo)
    C("eq:Rzp / eq:Rzp2", "NOTE: e^{k~x}e^{c} != e^{c}e^{k~x}: the constant must stand on the LEFT (matrices do not commute)", Rz_printed(0.3), Rz(0.3); kind=:note)

    # eq:Rze
    xoz, Δl = 0.7, 0.3
    xom = xoz + Δl / 2
    C("eq:Rze", "Rze = Rz(0) e^{k~ (xoz + Δℓ/2)}", Rz(xom), Rzz * exp(kt * xom))
    C("note after eq:rmr0", "Rz(xoz+u) = Rz(xoz) e^{k~ u}  (k^+ constant)", Rz(xoz + 0.2), Rz(xoz) * exp(kt * 0.2))
    C("eq:Rze", "Rze is a rotation tensor", Rz(xom) * Rz(xom)', I3)

    # eq:rmr0 : r(x)-r(0) = Rz(0) ε(x) a1 ; integral of Rz a1 by quadrature (independent of the closed form)
    eps_closed(k, x) = (κk = norm(k); kkk = k * k'; (I - kkk / κk^2) * sin(κk * x) / κk + tilde(k) / κk^2 * (1 - cos(κk * x)) + kkk / κk^2 * x)
    rq(x0, x1) = gaussint(σ -> Rz(σ) * e1, x0, x1)
    C("eq:rmr0", "r(x)-r(0) = Rz(0) ε(x) a1  (quadrature of Rz a1 vs closed form)", rq(0.0, 0.9), Rzz * eps_closed(kp, 0.9) * e1)
    C("eq:rmr0 / note", "r(xoz+Δℓ)-r(xoz) = Rz(xoz) ε(Δℓ) a1", rq(xoz, xoz + Δl), Rz(xoz) * eps_closed(kp, Δl) * e1)
    # eq:rmr02 (as corrected in the text: element form, (κΔℓ)²)
    Lsq = dot(rq(xoz, xoz + Δl), rq(xoz, xoz + Δl))
    f_sin(k, D) = (κk = norm(k); (k[1]^2 * (κk * D)^2 + 4 * (κk^2 - k[1]^2) * sin(κk * D / 2)^2) / κk^4)
    f_cos(k, D) = (κk = norm(k); (k[1]^2 * (κk * D)^2 + 2 * (κk^2 - k[1]^2) * (1 - cos(κk * D))) / κk^4)
    C("eq:rmr02", "Δℓ² = |r(xoz+Δℓ)-r(xoz)|² = [k1²(κΔℓ)² + 4(κ²-k1²) sin²(κΔℓ/2)]/κ⁴", [Lsq], [f_sin(kp, Δl)])
    C("eq:rmr02", "... = [k1²(κΔℓ)² + 2(κ²-k1²)(1-cos κΔℓ)]/κ⁴  (last line)", [Lsq], [f_cos(kp, Δl)])
    C("eq:rmr02", "(ε^T ε)_11 = same expression (Rz(xoz) orthogonal)", [(eps_closed(kp, Δl)' * eps_closed(kp, Δl))[1, 1]], [f_sin(kp, Δl)])
    # Newton solution of eq:rmr02 starting from the initial guess of eq:Deltal_k230
    chord = sqrt(Lsq)
    D = chord
    for it in 1:30
        g(d) = f_sin(kp, d) - chord^2
        D -= g(D) / ForwardDiff.derivative(g, D)
    end
    C("eq:Deltal_k230 (initial guess)", "Newton on eq:rmr02 from Δℓ0 = |Δr| converges to the true Δℓ", [D], [Δl])
    checkb(s, "eq:Deltal_k230", "initial guess |Δr| <= Δℓ (chord shorter than arc)", chord <= Δl)
    # special cases
    k0 = [0.0, 0.3, 0.5]; κ0 = norm(k0)
    ch0 = norm(gaussint(σ -> Rzz * (I + sin(κ0 * σ) * tilde(k0) / κ0 + (1 - cos(κ0 * σ)) * tilde(k0)^2 / κ0^2) * e1, 0.0, Δl))   # untwisted: root chord = element chord
    C("eq:Deltal_k10", "k1 = 0: Δℓ = (2/κ) asin(κ |Δr| / 2)", [2 / κ0 * asin(κ0 / 2 * ch0)], [Δl])
    k1 = [0.4, 0.0, 0.0]; κ1 = norm(k1)
    ch1 = norm(gaussint(σ -> Rzz * (I + sin(κ1 * σ) * tilde(k1) / κ1 + (1 - cos(κ1 * σ)) * tilde(k1)^2 / κ1^2) * e1, 0.0, Δl))
    C("eq:Deltal_k230", "k2 = k3 = 0 (straight, twisted): Δℓ = |Δr|", [ch1], [Δl])
    C("eq:rmr02", "k2=k3=0 reduces eq:rmr02 to Δℓ²", [f_sin(k1, Δl)], [Δl^2])

    # "constant over the element" approximation (note in sec. 1.5.1): error of replacing Rz by Rz(xom) is O(Δℓ²)
    err(h) = norm(gaussint(σ -> Rz(σ), xom - h / 2, xom + h / 2) / h - Rz(xom))
    ratio = err(0.4) / err(0.2)
    checkb(s, "eq:vars_e2 (note)", "mean of Rz over the element - Rz(midpoint) = O(Δℓ²) (error ratio ≈ 4 when Δℓ is halved; got $(round(ratio, digits=2)))", abs(ratio - 4) < 0.2)

    # ω_be and v_be (eq:vbe) against differentiating the position of the element midpoint
    ψA = Fld(41; scale=1.2); RA(t) = R_exp(ψA(0.0, t)); rA(t) = Fld(51; scale=1.0)(0.0, t)
    t0 = 0.61
    re = Rz(xom) * eps_closed(kp, 0.5) * e1 + [0.9, 0.2, -0.1]          # any constant vector (components in A)
    rbo(t) = rA(t) + RA(t) * re
    vbo = ForwardDiff.derivative(rbo, t0)
    RAdot = ForwardDiff.derivative(RA, t0); RAn = RA(t0)
    ωtilAo, ωtilA = RAdot * RAn', RAn' * RAdot
    vAo = ForwardDiff.derivative(rA, t0)
    C("eq:vbe", "v_be = R_A^T v_b^o (independent: derivative of r_A^o + R_A r_e)", RAn' * vbo, RAn' * vAo + ωtilA * re)
    C("eq:vbe", "v_be = R_A^T (v_A^o + ω~_A^o R_A r_e)  (second line)", RAn' * vbo, RAn' * (vAo + ωtilAo * RAn * re))
    Rb(t) = RA(t) * Rzz * expk(xom)
    ωtilb = ForwardDiff.derivative(Rb, t0) * Rb(t0)'
    C("eq:vars_e2", "ω_be = ω_b = ω_A = R_A^T ω_A^o: angular velocity of b (R_A Rz) equals that of A (Rz constant in time)", RAn' * axial(ωtilb), axial(ωtilA))
    return nothing
end

# =============================================================================
# kinematic machinery shared by C.2 / C.3: random fields u(x,t), p(x,t), R_0(x), r(x), frame A
# =============================================================================
function kin_setup()
    θ0 = Fld(11; scale=1.6)
    R0(x) = R_exp(θ0(x, 0.0))
    gx, gw = QuadGK.gauss(Float64, 48)
    r0 = [0.3, -0.2, 0.5]
    rr(x) = r0 + (x / 2) * reduce(+, (gw[i] * R0((x / 2) * (gx[i] + 1))[:, 1] for i in eachindex(gx)))
    ufld, pfld = Fld(21; scale=0.4), Fld(31; scale=0.9)
    ψA = Fld(41; scale=1.2); RA(t) = R_exp(ψA(0.0, t)); rA(t) = Fld(51; scale=1.0)(0.0, t)
    Dx(f, x, t) = ForwardDiff.derivative(σ -> f(σ, t), x)
    Dt(f, x, t) = ForwardDiff.derivative(σ -> f(x, σ), t)
    R(x, t) = R_WM(pfld(x, t))
    RRz(x, t) = R(x, t) * R0(x)
    Rpl(x, t) = R0(x)' * R(x, t) * R0(x)
    upl(x, t) = R0(x)' * ufld(x, t)
    ωAtil(t) = RA(t)' * ForwardDiff.derivative(RA, t)
    ωA(t) = axial(ωAtil(t))
    vAo(t) = ForwardDiff.derivative(rA, t)
    vbA(x, t) = RA(t)' * vAo(t) + ωAtil(t) * rr(x)
    Vs(x, t) = Rpl(x, t)' * (R0(x)' * vbA(x, t) + tilde(R0(x)' * ωA(t)) * upl(x, t) + Dt(upl, x, t))                # eq:V_star
    Os(x, t) = axial(Rpl(x, t)' * Dt(Rpl, x, t)) + Rpl(x, t)' * (R0(x)' * ωA(t))                                   # eq:Omega_star_rotparp
    γp(x, t) = R0(x)' * R(x, t)' * (R0(x)[:, 1] + Dx(ufld, x, t)) - e1                                              # eq:gamp_A
    κp(x, t) = R0(x)' * axial(R(x, t)' * Dx(R, x, t))                                                              # eq:kappa_plus_A
    Kt(x, t) = RRz(x, t)' * Dx(RRz, x, t)                                                                          # eq:K_def (B^* basis)
    return (; R0, rr, ufld, pfld, RA, rA, Dx, Dt, R, RRz, Rpl, upl, ωAtil, ωA, Vs, Os, γp, κp, Kt)
end

# =============================================================================
# C.2  Nodal outputs
# =============================================================================
function appendixC2()
    s = "C.2"
    C(label, desc, a, b; kw...) = check(s, label, desc, a, b; kw...)
    K = kin_setup()
    (; R0, ufld, pfld, Dx, Dt, R, RRz, ωAtil, Vs, Os, γp, κp, Kt) = K
    x0, t0 = 0.37, 0.61

    # --- generalized strain-displacement relations in basis A (used by eq:un, eq:pn)
    up = Dx(ufld, x0, t0)
    C("eq:un", "u' = R Rz (a1 + γ) - Rz a1   (u' by AD of the field; γ^+ from eq:gamp_A)", up, RRz(x0, t0) * (e1 + γp(x0, t0)) - R0(x0) * e1)
    pprime = Dx(pfld, x0, t0)
    C("eq:pn", "p' = H^{-T} Rz κ   (p' by AD of the field; κ^+ from eq:kappa_plus_A)", pprime, HinvT_WM(pfld(x0, t0)) * R0(x0) * κp(x0, t0))
    C("eq:pn", "H^T(p) from the closed form (eq:HT_WM) = tangent operator by definition", HT_WM(pfld(x0, t0)), H_num(R_WM, pfld(x0, t0))')
    # first-order finite differences: error O(h²)
    for (fld, name, label) in ((ufld, "u", "eq:un"), (pfld, "p", "eq:pn"))
        pr = Dx(fld, x0, t0)
        errh(h) = norm(fld(x0 - h / 2, t0) - (fld(x0, t0) - h / 2 * pr))
        ratio = errh(0.2) / errh(0.1)
        checkb(s, label, "first-node value $(name)_n = $(name)_e - Δℓ/2 $(name)'_e is O(Δℓ²) accurate (error ratio ≈ 4; got $(round(ratio, digits=2)))", abs(ratio - 4) < 0.5)
    end
    # compatibility equations <=> equality of the nodal values obtained from both elements of an interior node (eq:F_middle)
    rng = MersenneTwister(202)
    Δ1, Δ2 = 0.31, 0.22
    ue, Qe, Qn = randn(rng, 3), randn(rng, 3), randn(rng, 3)
    un_ = ue + Δ1 / 2 * Qe + Δ2 / 2 * Qn                 # imposes FF^->_e + FF^<-_{e+1} = 0  (eq:FF with upper/lower signs)
    FFr_e = -ue - Δ1 / 2 * Qe
    FFl_n = un_ - Δ2 / 2 * Qn
    C("eq:F_middle / eq:un", "FF^->_e + FF^<-_{e+1} = 0  <=>  u_n from element e (+) equals u_n from element e+1 (-)", FFr_e + FFl_n, zeros(3))
    C("eq:un", "u_n (element e, second node) = u_n (element e+1, first node)", ue + Δ1 / 2 * Qe, un_ - Δ2 / 2 * Qn)

    # --- strong form behind eqs. (Fsn)-(Msn2): the A-form used by the discrete equations (eq:Fu, eq:Fp) vs the B^*-form of the text
    Fs, Ms_, Ps, Hs = Fld(81), Fld(82), Fld(83), Fld(84)
    ωA = K.ωA
    fA(x, t) = Dt((xx, tt) -> RRz(xx, tt) * Ps(xx, tt), x, t) + tilde(ωA(t)) * RRz(x, t) * Ps(x, t) - Dx((xx, tt) -> RRz(xx, tt) * Fs(xx, tt), x, t)
    function mA(x, t; variant=:ok)
        γt = tilde(γp(x, t))
        sg = variant == :flip ? -1.0 : 1.0
        Dt((xx, tt) -> RRz(xx, tt) * Hs(xx, tt), x, t) + tilde(ωA(t)) * RRz(x, t) * Hs(x, t) +
        RRz(x, t) * (tilde(Vs(x, t)) * Ps(x, t) - sg * (tilde(e1) + γt) * Fs(x, t)) - Dx((xx, tt) -> RRz(xx, tt) * Ms_(xx, tt), x, t)
    end
    function fstar(x, t; noK=false)                       # F' = Pdot + Ω~ P - K~ F - f*  =>  f* = Pdot + Ω~ P - F' - K~ F
        Kf = noK ? zeros(3, 3) : Kt(x, t)
        Dt(Ps, x, t) + tilde(Os(x, t)) * Ps(x, t) - Dx(Fs, x, t) - Kf * Fs(x, t)
    end
    function mstar(x, t; noK=false)
        Kf = noK ? zeros(3, 3) : Kt(x, t)
        Dt(Hs, x, t) + tilde(Os(x, t)) * Hs(x, t) + tilde(Vs(x, t)) * Ps(x, t) - (tilde(e1) + tilde(γp(x, t))) * Fs(x, t) - Dx(Ms_, x, t) - Kf * Ms_(x, t)
    end
    C("sec:nodal_out / eq:Fsn", "A-form (Fu): f = (R Rz P)˙ + ω~ R Rz P - (R Rz F)'   equals   R Rz f*,  f* = Ṗ* + Ω~* P* - F*' - K~* F*", fA(x0, t0), RRz(x0, t0) * fstar(x0, t0))
    C("sec:nodal_out / eq:Msn", "A-form (Fp): m = (R Rz H)˙ + ω~ R Rz H + R Rz(V~P-(a~1+γ~)F) - (R Rz M)'   equals   R Rz m*, m* = Ḣ* + Ω~*H* + V~*P* - (a~1+γ~)F* - M*' - K~* M*", mA(x0, t0), RRz(x0, t0) * mstar(x0, t0))
    C("eq:Fsn", "K~* = k~^+ + κ~^+  (K~*_e in the text), k~^+ = R_0^T R_0'", Kt(x0, t0), R0(x0)' * Dx((xx, tt) -> R0(xx), x0, t0) + tilde(κp(x0, t0)))
    # mutation tests (the check must be able to fail)
    C("MUTAÇÃO", "dropping K~* F* in F' breaks the identity", fA(x0, t0), RRz(x0, t0) * fstar(x0, t0; noK=true); kind=:typo)
    C("MUTAÇÃO", "flipping the sign of (a~1+γ~)F in M' breaks the identity", mA(x0, t0; variant=:flip), RRz(x0, t0) * mstar(x0, t0); kind=:typo)

    # --- nodal values F_n^*, M_n^* (eqs. Fsn2/Msn2): accuracy of the formula against the exact field values
    xm = 0.37
    function nodal_errors(h; noK=false)
        Fe, Me, Pe, He = Fs(xm, t0), Ms_(xm, t0), Ps(xm, t0), Hs(xm, t0)
        RRze = RRz(xm, t0)
        Ωe, Ve, γe = Os(xm, t0), Vs(xm, t0), γp(xm, t0)
        Ke = noK ? zeros(3, 3) : Kt(xm, t0)
        Pdot, Hdot = Dt(Ps, xm, t0), Dt(Hs, xm, t0)
        Fprime = Pdot + tilde(Ωe) * Pe - Ke * Fe
        Mprime = Hdot + tilde(Ωe) * He + tilde(Ve) * Pe - (tilde(e1) + tilde(γe)) * Fe - Ke * Me
        # nodal resultants of the distributed loads (A components), eq:fmp_mmp, weights φ<- = 1-ζ and φ-> = ζ
        xz(ζ) = xm + h * (ζ - 1 / 2)
        fres(φ) = h * gaussint(ζ -> φ(ζ) * fA(xz(ζ), t0), 0.0, 1.0; n=8)
        mres(φ) = h * gaussint(ζ -> φ(ζ) * mA(xz(ζ), t0), 0.0, 1.0; n=8)
        fL, fRr, mL, mRr = fres(ζ -> 1 - ζ), fres(ζ -> ζ), mres(ζ -> 1 - ζ), mres(ζ -> ζ)
        F1 = Fe - (h / 2 * Fprime - RRze' * fL);  F2 = Fe + (h / 2 * Fprime - RRze' * fRr)
        M1 = Me - (h / 2 * Mprime - RRze' * mL);  M2 = Me + (h / 2 * Mprime - RRze' * mRr)
        eF = max(norm(F1 - Fs(xz(0.0), t0)), norm(F2 - Fs(xz(1.0), t0)))
        eM = max(norm(M1 - Ms_(xz(0.0), t0)), norm(M2 - Ms_(xz(1.0), t0)))
        return eF, eM
    end
    eF1, eM1 = nodal_errors(0.2); eF2, eM2 = nodal_errors(0.1)
    checkb(s, "eq:Fsn2", "F_n^* (formula with f_e^<-> nodal resultants, K~* term) is O(Δℓ²) accurate: error ratio $(round(eF1 / eF2, digits=2)) ≈ 4", abs(eF1 / eF2 - 4) < 0.7)
    checkb(s, "eq:Msn2", "M_n^* is O(Δℓ²) accurate: error ratio $(round(eM1 / eM2, digits=2)) ≈ 4", abs(eM1 / eM2 - 4) < 0.7)
    nF1, nM1 = nodal_errors(0.2; noK=true); nF2, nM2 = nodal_errors(0.1; noK=true)
    checkb(s, "eq:Fsn2 (K~ term)", "without the K~* F* term F_n^* is only O(Δℓ) (ratio $(round(nF1 / nF2, digits=2)) ≈ 2): the term is needed (basis B^* rotates between midpoint and node)", abs(nF1 / nF2 - 2) < 0.5)
    checkb(s, "eq:Msn2 (K~ term)", "same for M_n^* (ratio $(round(nM1 / nM2, digits=2)) ≈ 2)", abs(nM1 / nM2 - 2) < 0.5)

    # --- boundary nodes: relation between F_n^* and the hatted (applied) loads of eq:F_boundaries (algebra with random element data)
    Δℓ = 0.31
    Re, Rze = R_WM(0.5 * randn(rng, 3)), R_exp(randn(rng, 3))
    RRze = Re * Rze
    Fe, Me, Pe, He, Vs_, Om = ntuple(_ -> randn(rng, 3), 6)
    γe = 0.1 * randn(rng, 3); dotP, dotH, ωb = ntuple(_ -> randn(rng, 3), 3)
    Kte = tilde(0.3 * randn(rng, 3))
    f = randn(rng, 3)                                                  # constant distributed load, resolved in A
    fL = fR = Δℓ / 2 * f
    XP = dotP + tilde(ωb) * RRze * Pe                                  # (R Rz P)˙ + ω~ R Rz P  (dotP stands for the first term)
    Fu_left = -RRze * Fe - fL + Δℓ / 2 * XP                            # eq:Fu, first node
    Fu_right = RRze * Fe - fR + Δℓ / 2 * XP
    # F_n^* with the K~-free formula (same basis as the element midpoint): A components  R Rz F_n^*
    Xstar = RRze' * XP                                                 # Pdot + Ω~ P in B^* (by definition of the data)
    F1 = Fe - (Δℓ / 2 * Xstar - RRze' * fL);  F2 = Fe + (Δℓ / 2 * Xstar - RRze' * fR)
    C("sec:nodal_out (BC)", "last node: R Rz F_n^* = +F^_{N+1}  (from Fu^-> - F^_{N+1} = 0)", RRze * F2, Fu_right)
    C("sec:nodal_out (BC)", "first node: R Rz F_n^* = -F^_1  (from Fu^<- - F^_1 = 0): opposite sign", RRze * F1, -Fu_left)
    C("sec:nodal_out (BC)", "NOTE: 'F_n^* = F^_n^*' holds with + sign only at the last node; at the first node F_n^* = -(R Rz)^T F^_1", RRze * F1, Fu_left; kind=:note)
    return nothing
end

# =============================================================================
# C.3  Loading input
# =============================================================================
function appendixC3()
    s = "C.3"
    C(label, desc, a, b; kw...) = check(s, label, desc, a, b; kw...)
    rng = MersenneTwister(303)
    R = R_WM(1.2 * randn(rng, 3)); Rz = R_exp(randn(rng, 3))
    f0, f0p = randn(rng, 3), randn(rng, 3)
    # eq:fd .. eq:ffp : a vector attached to the (deformed) beam
    C("eq:ffs", "follower (defined in A): f_f = R Rz R_z^T f_f0 = R f_f0", R * Rz * (Rz' * f0), R * f0)
    C("eq:ffp", "follower (defined in b^+, B^* and b^+ aligned at t=0): f_f = R Rz f_f0^+", R * Rz * f0p, R * (Rz * f0p))
    C("eq:fdp", "dead in b^+: A components Rz f_d0^+ (b^+ does not move)", Rz * f0p, Rz * f0p)
    # initial coincidence of the descriptions at the undeformed state R = I
    C("eq:ffs/eq:ffp", "R = I: follower in A and in b^+ coincide when f_f0 = Rz f_f0^+", I3 * f0, I3 * Rz * (Rz' * f0))

    # eq:fa : aerodynamic basis; Rw (A components) = Rz Rw^+ Rz^T
    Rwp = R_exp(randn(rng, 3))                      # Rw^+ (components in b^+)
    Rw = Rz * Rwp * Rz'
    fsq = randn(rng, 3)                             # f^□ components in the deformed aerodynamic basis W
    # W_i = 𝐑 w_i^+ : A components of w_j^+ are the columns of Rz Rw^+ ; applying the deformation R to them:
    Wcols = R * (Rz * Rwp)
    C("eq:fa", "f_a = R Rw Rz f^□  equals  W f^□ (columns of W = deformed aerodynamic basis in A)", R * Rw * Rz * fsq, Wcols * fsq)
    C("eq:fa", "R Rw Rz = R Rz Rw^+  (consequence of Rw = Rz Rw^+ Rz^T)", R * Rw * Rz, R * Rz * Rwp)
    vstar = (R * Rz)' * (Wcols * fsq)               # B^* components of the same vector
    C("eq:vstar_to_vsquare", "v^* = Rw^+ v^□   (b^+ components of Rw)", vstar, Rwp * fsq)
    C("eq:vstar_to_vsquare / eq:Us (Ch.2)", "NOTE: with Rw in A components, v^□ = Rw^T v^* would be WRONG (needs Rw^{+T}): differs unless Rz commutes", Rw' * vstar, Rwp' * vstar; kind=:note)

    # gravity: direct sum over the point masses of a section
    nP = 30
    ρdS = 0.5 .+ rand(rng, nP); ξ = [[0.0, 0.3 * (rand(rng) - 0.5), 0.1 * (rand(rng) - 0.5)] for _ in 1:nP]
    μ = sum(ρdS); η = sum(ρdS[i] * ξ[i] for i in 1:nP) / μ
    g = 9.81; i3 = [0.0, 0.0, 1.0]
    RA = R_exp(randn(rng, 3)); Re = R_WM(0.8 * randn(rng, 3)); Rze = Rz; RRze = Re * Rze
    gA = RA' * (-g * i3)                                                   # gravity acceleration in A components
    fg_direct = sum(ρdS[i] * gA for i in 1:nP)
    mg_direct = sum(ρdS[i] * tilde(RRze * ξ[i]) * gA for i in 1:nP)        # Σ m_i (R Rz ξ_i) × g
    fgo = -μ * g * i3
    C("eq:fg_final", "f_g = R_A^T f_g^o = -μ g R_A^T i3  (direct sum of point masses)", fg_direct, RA' * fgo)
    etao = RA * RRze * tilde(η) * RRze' * RA'
    C("eq:etatiloe", "η~^o = R_A (R Rz) η~^* (R Rz)^T R_A^T  (transformation rule)", etao, tilde(RA * RRze * η))
    mgo = etao * fgo ./ 1.0
    C("eq:mgo", "m_g^o = η~^o f_g^o  (moment of the weight about the reference line; with μ η = Σ m_i ξ_i)", RA' * tilde(RA * RRze * η) * fgo, mg_direct)
    C("eq:mg_final", "m_g = (R Rz) η~^* (R Rz)^T f_g   (direct sum of Σ m_i (R Rz ξ_i) × g_A)", RRze * tilde(η) * RRze' * fg_direct, mg_direct)
    C("eq:fmg_nodal", "nodal resultants of a constant distributed load: ∫ φ f Δℓ dζ = Δℓ/2 f", [0.37 * gaussint(ζ -> 1 - ζ, 0.0, 1.0)], [0.37 / 2])
    C("eq:fgo", "-1/Δℓ ∫_0^{Δℓ} μ g i3 dx1 = -μ g i3", [-(1 / 0.37) * gaussint(σ -> 2.0 * g, 0.0, 0.37)], [-2.0 * g])

    # R_A update (eqs. RAdot, RApartial, RA): trapezoid vs Cayley transform; accuracy and orthogonality
    ωo(t) = [sin(1.3t) + 0.5, cos(0.7t), 0.4 + 0.3 * sin(2t)]
    Ωt(t) = tilde(ωo(t))
    RA0 = R_exp([0.3, -0.5, 0.2])
    function rk4(T, N)
        Rm = copy(RA0); h = T / N; t = 0.0
        for _ in 1:N
            k1 = Ωt(t) * Rm; k2 = Ωt(t + h / 2) * (Rm + h / 2 * k1); k3 = Ωt(t + h / 2) * (Rm + h / 2 * k2); k4 = Ωt(t + h) * (Rm + h * k3)
            Rm = Rm + h / 6 * (k1 + 2k2 + 2k3 + k4); t += h
        end
        Rm
    end
    T = 1.6
    Rex = rk4(T, 20000)
    function march(scheme, N)
        Rm = copy(RA0); Δt = T / N; a = 2 / Δt; t = 0.0
        for _ in 1:N
            Ωn, Ωn1, Ωm = Ωt(t), Ωt(t + Δt), Ωt(t + Δt / 2)
            Rm = scheme == :cayley_end ? (a * I - Ωn1) \ ((a * I + Ωn1) * Rm) :          # previous version of eq:RA / old code (ω~ at the end of the step)
                 scheme == :trapezoid ? (a * I - Ωn1) \ ((a * I + Ωn) * Rm) :            # plain trapezoidal rule (ω~(t) and ω~(t+Δt) at the two ends)
                                        (a * I - Ωm) \ ((a * I + Ωm) * Rm)               # Cayley with ω at mid-step: eq:RA / eq:RAdot_mid as now written, and as coded in Problem.jl
            t += Δt
        end
        Rm
    end
    errs = Dict(sch => [norm(march(sch, N) - Rex) for N in (16, 32, 64)] for sch in (:cayley_end, :trapezoid, :cayley_mid))
    orth = Dict(sch => norm(march(sch, 32) * march(sch, 32)' - I) for sch in (:cayley_end, :trapezoid, :cayley_mid))
    order(sch) = log2(errs[sch][2] / errs[sch][3])
    checkb(s, "eq:RA (old form)", "Cayley with ω~(t+Δt) on both sides: orthogonal to round-off (|R R^T - I| = $(@sprintf("%.1e", orth[:cayley_end])))", orth[:cayley_end] < 1e-12)
    checkb(s, "eq:RA (old form)", "... but only FIRST-order accurate for time-varying ω (observed order $(round(order(:cayley_end), digits=2)))", abs(order(:cayley_end) - 1) < 0.2)
    checkb(s, "eq:RAdot", "NOTE: plain trapezoid rule (ω~(t) and ω~(t+Δt)): 2nd order ($(round(order(:trapezoid), digits=2))) but NOT orthogonal (|R R^T - I| = $(@sprintf("%.1e", orth[:trapezoid]))) -> why eq:RAdot_mid uses the mid-step ω~", orth[:trapezoid] < 1e-12; kind=:note)
    checkb(s, "eq:RA / eq:RAdot_mid", "Cayley with ω~ at mid-step is orthogonal AND second order (observed $(round(order(:cayley_mid), digits=2)))", orth[:cayley_mid] < 1e-12 && abs(order(:cayley_mid) - 2) < 0.3)
    # constant ω: the three schemes coincide
    Ωc = tilde([0.3, -0.2, 0.9]); a = 2 / 0.1
    C("eq:RA", "for constant ω the printed form and the trapezoid form coincide", (a * I - Ωc) \ ((a * I + Ωc) * RA0), (a * I - Ωc) \ ((a * I + Ωc) * RA0))

    # steady-state remark at the end of C.3
    ω3 = 0.7
    RAt(t) = R_exp([0.0, 0.0, ω3 * t]) * RA0                                # rotation about the inertial vertical at constant rate
    dRAT_i3 = ForwardDiff.derivative(t -> RAt(t)' * i3, 0.4)
    C("sec:loads (steady-state)", "ω_A^o = ω3 i3 (constant)  =>  R_A^T i3 constant (gravity loads steady in A)", dRAT_i3, zeros(3))
    ωA = RAt(0.4)' * (ω3 * i3)
    C("sec:loads (steady-state)", "NOTE: ω_A = R_A^T ω_A^o = ω3 R_A^T i3, which is {0,0,ω3}^T only if R_A^T i3 = i3", ωA, [0.0, 0.0, ω3]; kind=:note)
    C("sec:loads (steady-state)", "general condition: d/dt(R_A^T i3) = -ω~_A R_A^T i3 = 0  <=>  ω_A ∥ R_A^T i3", -tilde(ωA) * (RAt(0.4)' * i3), zeros(3))
    return nothing
end

# =============================================================================
# C.4  Point inertias and springs
# =============================================================================
function appendixC4()
    s = "C.4"
    C(label, desc, a, b; kw...) = check(s, label, desc, a, b; kw...)
    rng = MersenneTwister(404)
    # rigid cluster of point masses -> sectional inertia matrix by direct summation (eqs. sec:mass, sec:mom_inertia_*, M)
    nP = 8
    m = 0.5 .+ rand(rng, nP); ξ = [0.3 * randn(rng, 3) for _ in 1:nP]
    m0 = sum(m); η0 = sum(m[i] * ξ[i] for i in 1:nP) / m0
    ϱ = sum(m[i] * (dot(ξ[i] - η0, ξ[i] - η0) * I - (ξ[i] - η0) * (ξ[i] - η0)') for i in 1:nP)   # inertia tensor about the cluster's own CM
    Mdirect = sum(m[i] * [I3 tilde(ξ[i])'; tilde(ξ[i]) tilde(ξ[i]) * tilde(ξ[i])'] for i in 1:nP)
    Δℓ = 0.4
    M0 = 1 / Δℓ * [m0 * I3 m0 * tilde(η0)'; m0 * tilde(η0) ϱ + m0 * tilde(η0) * tilde(η0)']
    C("eq:M_PM", "M_0 = 1/Δℓ [m I, m η~^T; m η~, ϱ + m η~ η~^T]  equals the direct sum over the cluster (÷Δℓ)", M0, Mdirect / Δℓ)
    C("eq:M_PM", "parallel-axis theorem: ι = ϱ + m η~ η~^T = ϱ + m(|η|² I - η η^T)", ϱ + m0 * tilde(η0) * tilde(η0)', ϱ + m0 * (dot(η0, η0) * I - η0 * η0'))
    M_e0 = let A = randn(rng, 6, 6); A * A' + 6I end
    C("eq:M_PM", "Mass_e = Mass_e0 + Mass_0 is symmetric if both are", (M_e0 + M0) - (M_e0 + M0)', zeros(6, 6))

    # springs (eqs. F_spring, M_spring and the Jacobians of sec:JM)
    Ku = let A = randn(rng, 3, 3); A * A' + I end
    uh, ur = randn(rng, 3), randn(rng, 3)
    Vspring(Δu) = 0.5 * dot(Δu, Ku * Δu)
    Fnode_correct = -ForwardDiff.gradient(Vspring, uh - ur)                 # restoring force on the node
    Fnode_printed = Ku * (uh - ur)                                          # eq:F_spring as printed
    C("eq:F_spring", "force on the node = -∂V/∂u = -K_u (û_n - u_r) (restoring)", Fnode_correct, -Ku * (uh - ur))
    C("eq:F_spring (as printed)", "printed F^_n = +K_u(û_n - u_r) is NOT the restoring force", Fnode_printed, Fnode_correct; kind=:typo)
    Jprinted = ForwardDiff.jacobian(x -> -(Ku * (x - ur)), uh)             # ∂(-F^_n)/∂û_n with F^ as printed
    Jcorrect = ForwardDiff.jacobian(x -> -(-(Ku * (x - ur))), uh)
    C("sec:JM springs", "printed ∂(-F^_n)/∂û_n = K_u is consistent with F^ = -K_u(û_n-u_r), not with the printed F^ = +K_u(...)", Jprinted, Ku; kind=:typo)
    C("sec:JM springs", "with F^_n = -K_u(û_n - u_r):  ∂(-F^_n)/∂û_n = K_u", Jcorrect, Ku)
    C("sec:JM springs (double attach.)", "∂(-F^_n)/∂u_r = -K_u (missing in the text: contribution of the other end, needed when u_r is an unknown)", ForwardDiff.jacobian(x -> Ku * (uh - x), ur), -Ku)

    # rotational spring: p_h = rotation parameters of R(p_n) R(p_r)^T (code: rotation_between_WM(p_r, p_n))
    pr, pn = [0.5, -0.4, 0.3], [0.2, 0.7, -0.5]
    ph(pn_, pr_) = p_from_R(R_WM(pn_) * R_WM(pr_)')
    C("eq:M_spring", "p_h composition: R(p_h) = R(p_n) R(p_r)^T", R_WM(ph(pn, pr)), R_WM(pn) * R_WM(pr)')
    checkb(s, "eq:M_spring", "p_h != p_n - p_r in general (as remarked in the text)", norm(ph(pn, pr) - (pn - pr)) > 1e-3)
    dph_dpn = ForwardDiff.jacobian(x -> ph(x, pr), pn)
    C("sec:JM springs", "NOTE: ∂p_h/∂p_n = I only if p_r = 0; for p_r != 0, ∂(-M^_n)/∂p_n = K_p ∂p_h/∂p_n != K_p", dph_dpn, I3; kind=:note)
    C("sec:JM springs", "∂p_h/∂p_n = I when p_r = 0 (then p_h = p_n)", ForwardDiff.jacobian(x -> ph(x, zeros(3)), pn), I3)
    C("sec:JM springs (double attach.)", "NOTE: ∂p_h/∂p_r != 0: the rotational spring also couples to the other node's rotation (not in the text)", ForwardDiff.jacobian(x -> ph(pn, x), pr), zeros(3, 3); kind=:note)
    # WM parameters vs rotation vector as the spring deformation measure
    ψ = π / 2
    pWM = 4 * tan(ψ / 4)
    checkb(s, "eq:M_spring", "NOTE: |p_h| = 4 tan(ψ/4) vs ψ: for ψ = 90° the spring moment K_p p_h overestimates K_p ψ by $(round(100 * (pWM / ψ - 1), digits=1))% (p_h is the WM parameter vector, not the rotation vector)", abs(pWM / ψ - 1) < 1e-3; kind=:note)
    return nothing
end

# =============================================================================
# C.5  Jacobian and inertia matrix entries: every block against forward-mode AD of the discrete residual
# =============================================================================
const ROWS = (:Fu1, :Fu2, :Fp1, :Fp2, :FF1, :FF2, :FM1, :FM2, :FV, :FO)
const COLS = (:u, :p, :F, :M, :V, :Om)
const RATES = (:ud, :pd, :Vd, :Omd)
const ROWDESC = Dict(:Fu1 => "Fu^<-", :Fu2 => "Fu^->", :Fp1 => "Fp^<-", :Fp2 => "Fp^->", :FF1 => "FF^<-", :FF2 => "FF^->", :FM1 => "FM^<-", :FM2 => "FM^->", :FV => "FV", :FO => "FO")

function appendixC5()
    s = "C.5"
    for trial in 1:2
        rng = MersenneTwister(500 + trial)
        Δℓ = 0.2 + 0.3 * rand(rng)
        Rze = R_exp(randn(rng, 3)); ωb, vb = randn(rng, 3), randn(rng, 3)
        # sectional matrices: compliance (SPD) and an inertia matrix of a point-mass cluster
        Sm = let A = randn(rng, 6, 6); A * A' + 6I end
        Cm = inv(Sm)
        nP = 8; mk = 0.5 .+ rand(rng, nP); ξk = [0.3 * randn(rng, 3) for _ in 1:nP]
        Mm = sum(mk[i] * [I3 tilde(ξk[i])'; tilde(ξk[i]) tilde(ξk[i]) * tilde(ξk[i])'] for i in 1:nP)
        M11, M12, M21, M22 = Mm[1:3, 1:3], Mm[1:3, 4:6], Mm[4:6, 1:3], Mm[4:6, 4:6]
        C11, C12, C21, C22 = Cm[1:3, 1:3], Cm[1:3, 4:6], Cm[4:6, 1:3], Cm[4:6, 4:6]
        μ = Mm[1, 1]; ηs = axial(Mm[4:6, 1:3]) / μ                      # η* from M21 = μ η~*
        Δt = 0.05; c = 2 / Δt
        # distributed loads (nodal resultants): follower in A, follower in b^+, dead, gravity
        f0AL, f0AR, m0AL, m0AR, f0BL, f0BR, m0BL, m0BR, fdL, fdR, mdL, mdR, fg = ntuple(_ -> randn(rng, 3), 13)
        loads(p) = begin
            local R, mgrav
            R = R_WM(p)
            mgrav = Δℓ / 2 * R * Rze * tilde(ηs) * Rze' * R' * fg                                  # eq:mg_final times Δℓ/2
            (R * f0AL + R * Rze * f0BL + fdL + Δℓ / 2 * fg, R * f0AR + R * Rze * f0BR + fdR + Δℓ / 2 * fg,
             R * m0AL + R * Rze * m0BL + mdL + mgrav, R * m0AR + R * Rze * m0BR + mdR + mgrav)
        end
        bu, bp, bV, bO = ntuple(_ -> randn(rng, 3), 4)
        # ---- elemental residual (eqs. Fu ... FH), coded from the discrete equations
        function elres(x, xd)
            local u, p, F, M, V, Om, ud, pd, Vd, Omd, R, RRz, HT, HiT, Rdot, γ, κ, P, H, Pd, Hd, f1, f2, m1, m2, XP, XH, Q, Ks     # (avoid capturing the outer variables)
            u, p, F, M, V, Om = x[1:3], x[4:6], x[7:9], x[10:12], x[13:15], x[16:18]
            ud, pd, Vd, Omd = xd[1:3], xd[4:6], xd[7:9], xd[10:12]
            R = R_WM(p); RRz = R * Rze; HT = HT_WM(p); HiT = HinvT_WM(p)
            Rdot = R * tilde(HT * pd)
            γ = (Cm * vcat(F, M))[1:3]; κ = (Cm * vcat(F, M))[4:6]
            P = (Mm * vcat(V, Om))[1:3]; H = (Mm * vcat(V, Om))[4:6]
            Pd = (Mm * vcat(Vd, Omd))[1:3]; Hd = (Mm * vcat(Vd, Omd))[4:6]
            f1, f2, m1, m2 = loads(p)
            XP = Rdot * Rze * P + RRz * Pd + tilde(ωb) * RRz * P
            XH = Rdot * Rze * H + RRz * Hd + tilde(ωb) * RRz * H + RRz * (tilde(V) * P - (tilde(e1) + tilde(γ)) * F)
            Q = RRz * (e1 + γ) - Rze * e1
            Ks = HiT * Rze * κ
            vcat(-RRz * F - f1 + Δℓ / 2 * XP, RRz * F - f2 + Δℓ / 2 * XP,
                 -RRz * M - m1 + Δℓ / 2 * XH, RRz * M - m2 + Δℓ / 2 * XH,
                 u - Δℓ / 2 * Q, -u - Δℓ / 2 * Q,
                 p - Δℓ / 2 * Ks, -p - Δℓ / 2 * Ks,
                 RRz * V - vb - tilde(ωb) * u - ud,
                 RRz * Om - ωb - R * HT * pd)
        end
        x0 = vcat(0.3 * randn(rng, 3), 0.8 * randn(rng, 3), randn(rng, 3), randn(rng, 3), randn(rng, 3), randn(rng, 3))
        rates(x) = vcat(c * x[1:3] + bu, c * x[4:6] + bp, c * x[13:15] + bV, c * x[16:18] + bO)      # Newmark (eq:Newmark_ydot3): ydot = c y + bar
        G(x) = elres(x, rates(x))
        J_AD = ForwardDiff.jacobian(G, x0)                                  # 30 x 18 total derivative
        xd0 = rates(x0)
        M_AD = ForwardDiff.jacobian(xd -> elres(x0, xd), xd0)               # 30 x 12 inertia matrix
        # consistency of the residual coding itself
        u, p, F, M, V, Om = x0[1:3], x0[4:6], x0[7:9], x0[10:12], x0[13:15], x0[16:18]
        ud, pd, Vd, Omd = xd0[1:3], xd0[4:6], xd0[7:9], xd0[10:12]
        R = R_WM(p); RRz = R * Rze; HT = HT_WM(p); HiT = HinvT_WM(p)
        γ = (Cm * vcat(F, M))[1:3]; κ = (Cm * vcat(F, M))[4:6]
        P = (Mm * vcat(V, Om))[1:3]; H = (Mm * vcat(V, Om))[4:6]
        Pd = (Mm * vcat(Vd, Omd))[1:3]; Hd = (Mm * vcat(Vd, Omd))[4:6]
        Rdot = R * tilde(HT * pd)
        check(s, "setup (t$trial)", "Rdot = R tilde(H^T pdot) equals d/dt R(p(t))", ForwardDiff.derivative(σ -> R_WM(p + σ * pd), 0.0), Rdot)
        check(s, "setup (t$trial)", "H^{-T} (closed form of C.6) = (H^T)^{-1}", HiT, inv(HT))
        check(s, "setup (t$trial)", "R^T Rdot = tilde(H^T pdot)  (H^T pdot = body-frame angular velocity)", R' * Rdot, tilde(HT * pd))

        # ---- derivatives of rotation quantities (third-order tensor times vector, eq:mul3)
        dR(v) = mul3(R_WM, p, v)                                                   # ∂R/∂p · v
        dHT(v) = mul3(HT_WM, p, v)                                                 # ∂H^T/∂p · v
        dHiT(v) = mul3(HinvT_WM, p, v)                                             # ∂H^{-T}/∂p · v
        dRdot_fixed(v) = mul3(q -> R_WM(q) * tilde(HT_WM(q) * pd), p, v)           # ∂Rdot/∂p with pdot held fixed (eq:dRdot_dpext of C.6: ∂pdot/∂p = 0)
        dRdot_tot(v) = dRdot_fixed(v) + c * dR(v)                                  # total derivative with pdot = c p + bar (Newmark)
        # derivatives of the loads (eqs. dfmp_dmmp, dmgp)
        df1 = dR(f0AL) + dR(Rze * f0BL); df2 = dR(f0AR) + dR(Rze * f0BR)
        dmgrav = Δℓ / 2 * (dR(Rze * tilde(ηs) * (R * Rze)' * fg) + (R * Rze * tilde(ηs) * Rze') * mul3(q -> R_WM(q)', p, fg))      # eq:dmgp
        dm1 = dR(m0AL) + dR(Rze * m0BL) + dmgrav; dm2 = dR(m0AR) + dR(Rze * m0BR) + dmgrav
        check(s, "eq:dmgp (t$trial)", "∂m_g^<->/∂p of eq:dmgp equals AD of m_g(p) = Δℓ/2 R Rz η~* Rz^T R^T f_g", dmgrav, ForwardDiff.jacobian(q -> Δℓ / 2 * R_WM(q) * Rze * tilde(ηs) * Rze' * R_WM(q)' * fg, p))
        check(s, "eq:dfmp_dmmp (t$trial)", "∂f^<->/∂p = ∂R/∂p (resultant of the follower load in the undeformed state), dead/gravity forces independent of p", hcat(df1, df2), hcat(ForwardDiff.jacobian(q -> loads(q)[1], p), ForwardDiff.jacobian(q -> loads(q)[2], p)))
        check(s, "eq:dfmp_dmmp (t$trial)", "∂m^<->/∂p (follower + gravity) vs AD of the moment resultants", hcat(dm1, dm2), hcat(ForwardDiff.jacobian(q -> loads(q)[3], p), ForwardDiff.jacobian(q -> loads(q)[4], p)))

        # ---- block formulas printed in the text. v = :printed | :corrected
        Ĩ = tilde(e1)
        Mdot = R * tilde(HT * pd) * R'   # Rdot R^T
        function textJ(v)
            J = Dict{Tuple{Symbol,Symbol},Matrix{Float64}}((r, cl) => zeros(3, 3) for r in ROWS for cl in COLS)
            for (rowF, sgn, df) in ((:Fu1, -1.0, df1), (:Fu2, 1.0, df2))                                     # sgn = ∓ (upper sign for <-)
                dRd = v == :printed ? dRdot_fixed : dRdot_tot
                J[(rowF, :p)] = sgn * dR(Rze * F) - df + Δℓ / 2 * (dRd(Rze * P) + tilde(ωb) * dR(Rze * P) + dR(Rze * Pd))     # eq:dFump_dp
                J[(rowF, :F)] = sgn * RRz                                                                      # eq:dFump_dF
                J[(rowF, :V)] = Δℓ / 2 * (Rdot * R' + tilde(ωb) + c * I3) * RRz * M11                           # eq:dFump_dV
                J[(rowF, :Om)] = Δℓ / 2 * (Rdot * R' + tilde(ωb) + c * I3) * RRz * M12                          # eq:dFump_dOm
            end
            for (rowF, sgn, dm) in ((:Fp1, -1.0, dm1), (:Fp2, 1.0, dm2))
                dRd = v == :printed ? dRdot_fixed : dRdot_tot
                J[(rowF, :p)] = sgn * dR(Rze * M) - dm + Δℓ / 2 * (dRd(Rze * H) + tilde(ωb) * dR(Rze * H) + dR(Rze * (Hd + tilde(V) * P - (Ĩ + tilde(γ)) * F)))   # eq:dFpmp_dp
                J[(rowF, :F)] = v == :printed ? Δℓ / 2 * RRz * ((Ĩ + tilde(γ)) - tilde(F) * C11) :
                                                Δℓ / 2 * RRz * (-(Ĩ + tilde(γ)) + tilde(F) * C11)             # eq:dFpmp_dF
                J[(rowF, :M)] = v == :printed ? sgn * RRz - Δℓ / 2 * RRz * tilde(F) * C12 :
                                                sgn * RRz + Δℓ / 2 * RRz * tilde(F) * C12                      # eq:dFpmp_dM
                J[(rowF, :V)] = Δℓ / 2 * ((Rdot * R' + tilde(ωb) + c * I3) * RRz * M21 + RRz * (tilde(V) * M11 - tilde(P)))   # eq:dFpmp_dV
                J[(rowF, :Om)] = Δℓ / 2 * ((Rdot * R' + tilde(ωb) + c * I3) * RRz * M22 + RRz * tilde(V) * M12)               # eq:dFpmp_dOm
            end
            for (rowF, sgn) in ((:FF1, -1.0), (:FF2, 1.0))                                                    # FF = ±u ...: +u for <- (sgn = -1)
                J[(rowF, :u)] = -sgn * I3                                                                       # eq:dFFmp_du
                J[(rowF, :p)] = -Δℓ / 2 * dR(Rze * (e1 + γ))                                                    # eq:dFFmp_dp
                J[(rowF, :F)] = -Δℓ / 2 * RRz * C11                                                             # eq:dFFmp_dF
                J[(rowF, :M)] = -Δℓ / 2 * RRz * C12                                                             # eq:dFFmp_dM
            end
            for (rowF, sgn) in ((:FM1, -1.0), (:FM2, 1.0))
                J[(rowF, :p)] = -sgn * I3 - Δℓ / 2 * dHiT(Rze * κ)                                              # eq:dFMmp_dp
                J[(rowF, :F)] = v == :printed ? Δℓ / 2 * HiT * Rze * C21 : -Δℓ / 2 * HiT * Rze * C21            # eq:dFMmp_dF
                J[(rowF, :M)] = v == :printed ? Δℓ / 2 * HiT * Rze * C22 : -Δℓ / 2 * HiT * Rze * C22            # eq:dFMmp_dM
            end
            J[(:FV, :u)] = -(tilde(ωb) + c * I3)                                                                 # eq:dFP_du
            J[(:FV, :p)] = dR(Rze * V)                                                                           # eq:dFP_dp
            J[(:FV, :V)] = RRz                                                                                   # eq:dFP_dV
            J[(:FO, :p)] = dR(Rze * Om - HT * pd) - R * (dHT(pd) + c * HT)                                       # eq:dFO_dp
            J[(:FO, :Om)] = RRz                                                                                  # eq:dFO_dOm
            return J
        end
        Jp, Jc = textJ(:printed), textJ(:corrected)
        rowrng(r) = (3 * (findfirst(==(r), ROWS) - 1) + 1):(3 * findfirst(==(r), ROWS))
        colrng(cl) = (3 * (findfirst(==(cl), COLS) - 1) + 1):(3 * findfirst(==(cl), COLS))
        for r in ROWS, cl in COLS
            blkAD = J_AD[rowrng(r), colrng(cl)]
            lab = "∂$(ROWDESC[r])/∂$(cl) (t$trial)"
            if norm(Jp[(r, cl)] - Jc[(r, cl)]) < 1e-12
                check(s, lab, "block of the elemental Jacobian vs AD of the discrete residual (Newmark included)", blkAD, Jp[(r, cl)])
            else
                check(s, lab, "AS PRINTED in the text", blkAD, Jp[(r, cl)]; kind=:typo)
                check(s, lab, "corrected expression", blkAD, Jc[(r, cl)])
            end
        end
        # ---- inertia matrix (derivatives with respect to the rates): eqs. dFump_dpdot ... dFaero_daerodot
        Mt = Dict{Tuple{Symbol,Symbol},Matrix{Float64}}((r, cl) => zeros(3, 3) for r in ROWS for cl in RATES)
        for (rowF) in (:Fu1, :Fu2)
            Mt[(rowF, :pd)] = Δℓ / 2 * dR(Rze * P)                      # eq:dFump_dpdot
            Mt[(rowF, :Vd)] = Δℓ / 2 * RRz * M11                        # eq:dFump_dVdot
            Mt[(rowF, :Omd)] = Δℓ / 2 * RRz * M12                       # eq:dFump_dOmdot
        end
        for (rowF) in (:Fp1, :Fp2)
            Mt[(rowF, :pd)] = Δℓ / 2 * dR(Rze * H)                      # eq:dFpmp_dpdot
            Mt[(rowF, :Vd)] = Δℓ / 2 * RRz * M21                        # eq:dFpmp_dVdot
            Mt[(rowF, :Omd)] = Δℓ / 2 * RRz * M22                       # eq:dFpmp_dOmdot
        end
        Mt[(:FV, :ud)] = -I3                                            # eq:dFP_dudot
        Mt[(:FO, :pd)] = -R * HT                                        # eq:dFO_dpdot
        ratecol(cl) = (3 * (findfirst(==(cl), RATES) - 1) + 1):(3 * findfirst(==(cl), RATES))
        for r in ROWS, cl in RATES
            check(s, "∂$(ROWDESC[r])/∂$(cl) (inertia, t$trial)", "block of the inertia matrix vs AD of the discrete residual", M_AD[rowrng(r), ratecol(cl)], Mt[(r, cl)])
        end
        # ---- nodal Jacobian terms (eq:FMhatn, eq:dfMhatn)
        Fn0, pn = randn(rng, 3), 0.9 * randn(rng, 3)
        check(s, "eq:dfMhatn (t$trial)", "follower nodal load F^_n = R(p_n) F^_n0:  ∂(-F^_n)/∂p_n = -∂R_n/∂p_n F^_n0", ForwardDiff.jacobian(q -> -R_WM(q) * Fn0, pn), -mul3(R_WM, pn, Fn0))
        check(s, "eq:mul3 (t$trial)", "(∂Λ/∂a) b := [∂Λ/∂a1 b, ∂Λ/∂a2 b, ∂Λ/∂a3 b] equals ∂(Λ b)/∂a", ForwardDiff.jacobian(q -> R_WM(q) * Fn0, pn), mul3(R_WM, pn, Fn0))
    end

    # ---- aerodynamic equations (structure only): F_a = χdot - A_a χ - B_a ; Newmark makes ∂/∂V a TOTAL derivative
    rng = MersenneTwister(555)
    nχ = 5; c = 2 / 0.05
    A0 = randn(rng, nχ, nχ); Aw = [randn(rng, nχ, nχ) for _ in 1:6]; Bw1, Bw2 = randn(rng, nχ, 6), randn(rng, nχ, 6)
    Aa(w) = A0 + sum(w[i] * Aw[i] for i in 1:6)                       # A_a depends on (V*, Ω*) only
    Ba(w, wd) = Bw1 * w + Bw2 * wd + 0.1 * sin.(Bw1 * w)             # B_a depends on (V*, Ω*) and on their rates (non-circulatory terms)
    Fa(χ, χd, w, wd) = χd - Aa(w) * χ - Ba(w, wd)
    χ0, w0, bχ, bw = randn(rng, nχ), randn(rng, 6), randn(rng, nχ), randn(rng, 6)
    tot(z) = Fa(z[1:nχ], c * z[1:nχ] + bχ, z[nχ+1:end], c * z[nχ+1:end] + bw)
    Jtot = ForwardDiff.jacobian(tot, vcat(χ0, w0))
    wd0 = c * w0 + bw
    dA_dw_χ = hcat((Aw[i] * χ0 for i in 1:6)...)
    dB_dw = ForwardDiff.jacobian(w -> Ba(w, wd0), w0); dB_dwd = ForwardDiff.jacobian(wd -> Ba(w0, wd), wd0)
    check(s, "eq:dFaero_daero", "∂F_a/∂χ = (2/Δt) I - A_a", Jtot[:, 1:nχ], c * I - Aa(w0))
    check(s, "eq:dFaero_dV", "AS PRINTED: ∂F_a/∂V = -∂A_a/∂V χ - ∂B_a/∂V (partial derivative only)", Jtot[:, nχ+1:end], -dA_dw_χ - dB_dw; kind=:typo)
    check(s, "eq:dFaero_dV", "total derivative: ... + (2/Δt) ∂/∂Vdot [-B_a] (Newmark; coded as F_χ_Vdot*c)", Jtot[:, nχ+1:end], -dA_dw_χ - dB_dw - c * dB_dwd)
    M_aero = ForwardDiff.jacobian(wd -> Fa(χ0, c * χ0 + bχ, w0, wd), wd0)
    check(s, "eq:dFaero_dVdot", "∂F_a/∂Vdot = -∂A_a/∂Vdot χ - ∂B_a/∂Vdot (A_a indep. of rates)", M_aero, -dB_dwd)
    check(s, "eq:dFaero_daerodot", "∂F_a/∂χdot = I", ForwardDiff.jacobian(xd -> Fa(χ0, xd, w0, wd0), c * χ0 + bχ), Matrix(1.0I, nχ, nχ))
    return nothing
end

# =============================================================================
println("="^100)
println("Verification of Appendix C (sections C.1 - C.5) of theory/main.tex   -   tolerance = $TOL")
println("="^100)
appendixC1()
appendixC2()
appendixC3()
appendixC4()
appendixC5()

# ----------------------------------------------------------------------------- report
println()
let lastsec = ""
    for r in RESULTS
        mark = r.kind == :pass ? (r.ok ? "ok" : "FALHOU") :
               r.kind == :typo ? (r.ok ? "ERRO NO TEXTO CONFIRMADO" : "!! esperava-se discrepancia e nao houve") :
                                 (r.ok ? "NOTA (afirmacao imprecisa/nao geral, ver relatorio)" : "!! esperava-se diferenca e nao houve")
        if r.sec != lastsec
            println("\n------------------------- ", r.sec, " -------------------------")
            lastsec = r.sec
        end
        @printf("%-5s %-34s %-110s %9.2e  %s\n", r.sec, first(r.label, 34), first(r.desc, 110), r.err, mark)
    end
end
nok = count(r -> r.kind == :pass && r.ok, RESULTS)
nbad = count(r -> r.kind == :pass && !r.ok, RESULTS)
ntypo = count(r -> r.kind == :typo && r.ok, RESULTS)
nnote = count(r -> r.kind == :note && r.ok, RESULTS)
nmiss = count(r -> r.kind != :pass && !r.ok, RESULTS)
println("\nSUMMARY: ", nok, " checks passed, ", nbad, " FAILED, ", ntypo, " errors-in-text confirmed (formula as printed disagrees), ", nnote, " notes, ", nmiss, " expected discrepancies that did not occur")
