# =============================================================================
# Numerical verification of Chapter 1 ("Beam model"), sections 1.1 to 1.4.3, of theory/main.tex
# =============================================================================
# Method
# ------
# Every identity of the text (including each intermediate line of each derivation) is evaluated
# numerically on random smooth fields: u(x1,t), p(x1,t), R_0(x1), r(x1), the motion of frame A,
# and a random cross-sectional mass distribution. All derivatives (d/dx1, d/dt) and variations
# (delta = d/d(epsilon)) are obtained with forward-mode automatic differentiation (ForwardDiff.jl),
# so they are exact to machine precision (same idea as jax.jvp / jax.jacfwd). Left- and right-hand
# sides are compared; a step "passes" when the relative difference is below TOL.
#
# Result categories (last column of the report)
#   ok                         the two sides agree (the step of the text is correct)
#   ERRO NO TEXTO CONFIRMADO   a formula evaluated exactly AS WRITTEN in main.tex disagrees (typo); the
#                              corrected expression is checked right after and passes
#   NOTA                       two quantities that the text could be read as equating are in fact different
#                              (supports a conceptual remark; not an error of the algebra)
#
# Independent (non-algebraic) tests are included so that the check does not merely re-encode the text:
# the velocity of every material point is obtained by differentiating its position; strains must vanish
# for a rigid-body motion of the whole beam; P^*, H^* are compared with int(rho v_P) and int(rho xi x v_P);
# the variations of Appendix A are compared with d/d(epsilon) of the perturbed strains/velocities.
# (A mutation test with 7 deliberately wrong formulas -- sign, transpose -- was caught in all cases.)
#
# Everything is repeated for two different rotation parametrizations (Wiener-Milenkovic, the one of
# the code, and the rotation vector) because Chapter 1 should not depend on the parametrization.
#
# Run:   julia --project=<AeroBeams root> theory/verification/verify_chapter1.jl      (~30 s)
# Symbols: Rp = R^+, up = u^+, pp = p^+, Rz = R_0, Rb = R_b^o = R_A R_0, "o" = inertial basis (circle),
#          D(f) = derivative, primes on b^+-components are derivatives OF THE COMPONENTS (e.g. u'^+ = d(u^+)/dx1)
# Last run (2026-11-14): 247 checks passed, 0 failed, 5 errors-in-text confirmed, 6 notes.
# =============================================================================

using LinearAlgebra, ForwardDiff, QuadGK, Random, Printf

const TOL = 1e-9

# ----------------------------------------------------------------------------- notation
tilde(v) = [0 -v[3] v[2]; v[3] 0 -v[1]; -v[2] v[1] 0]      # eq:skew-sym_def
axial(S) = [S[3, 2], S[1, 3], S[2, 1]]                      # inverse of tilde
const e1 = [1.0, 0.0, 0.0]                                  # \aone

# ----------------------------------------------------------------------------- rotation parametrizations
# Wiener-Milenkovic, scaled set (|p|<4): eq:R_WM / eq:T_explicit with p0 = 2 - p'p/8
function R_WM(p)
    p0 = 2 - dot(p, p) / 8
    υ = 1 / (4 - p0)
    υ^2 * ((1 / υ^2 - 2 * dot(p, p)) * I + 2 * p * p' + 2 * p0 * tilde(p))
end
# Rotation vector (exponential map)
function R_exp(p)
    φ = norm(p)
    N = tilde(p / φ)
    I + sin(φ) * N + (1 - cos(φ)) * N * N
end
# Tangent operator defined by its meaning: axial(Rdot R') = H pdot  (columns obtained by AD)
function H_num(Rfun, p)
    R = Rfun(p)
    J = ForwardDiff.jacobian(q -> vec(Rfun(q)), p)
    hcat((axial(reshape(J[:, j], 3, 3) * R') for j in 1:3)...)
end

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
    param::String
    sec::String
    label::String
    desc::String
    err::Float64
    kind::Symbol      # :pass (must agree) | :typo (formula AS WRITTEN must disagree) | :note (two quantities that must differ)
    ok::Bool
end
const RESULTS = Res[]
function check(param, sec, label, desc, a, b; kind=:pass)
    err = maximum(abs.(a .- b)) / max(1.0, maximum(abs.(a)), maximum(abs.(b)))
    ok = kind == :pass ? err < TOL : err > 1e-6
    push!(RESULTS, Res(param, sec, label, desc, err, kind, ok))
    return ok
end

# =============================================================================
function run_checks(Rfun, pname)
    C(sec, label, desc, a, b; kw...) = check(pname, sec, label, desc, a, b; kw...)
    Dx(f, x, t) = ForwardDiff.derivative(s -> f(s, t), x)      # d/dx1
    Dt(f, x, t) = ForwardDiff.derivative(s -> f(x, s), t)      # d/dt

    # ------------------------------------------------------------- fields
    θ0 = Fld(11; scale=1.6)                      # R_0(x1) = exp(θ0(x1))  (undeformed beam basis b^+ w.r.t. A)
    R0(x) = R_exp(θ0(x, 0.0))
    gx, gw = QuadGK.gauss(Float64, 48)           # arclength-consistent reference line: r' = b_1 = R_0 e_1
    r0 = [0.3, -0.2, 0.5]
    rr(x) = r0 + (x / 2) * reduce(+, (gw[i] * R0((x / 2) * (gx[i] + 1))[:, 1] for i in eachindex(gx)))
    ufld = Fld(21; scale=0.4)                    # u(x1,t)  (components in A)
    pfld = Fld(31; scale=0.9)                    # p(x1,t)  (components in A)
    ψA = Fld(41; scale=1.2)                      # frame A rotation vector (function of t)
    RA(t) = R_exp(ψA(0.0, t))
    rA(t) = Fld(51; scale=1.0)(0.0, t)           # r_A^o(t)

    x0, t0 = 0.37, 0.61
    u(x, t) = ufld(x, t)
    p(x, t) = pfld(x, t)
    R(x, t) = Rfun(p(x, t))                      # R (components in A)
    Rp(x, t) = R0(x)' * R(x, t) * R0(x)          # R^+ (transformation rule, 2nd-order tensor)
    up(x, t) = R0(x)' * u(x, t)                  # u^+
    pp(x, t) = R0(x)' * p(x, t)                  # p^+  (eq:rotpar: p = R_0 p^+)
    rp(x) = R0(x)' * rr(x)                       # r^+
    Rb(x, t) = RA(t) * R0(x)                     # R_b^o = R_A R_0  (eq:RboT); also R_b^+
    rbo(x, t) = rA(t) + RA(t) * rr(x)            # r_b^o = r_A^o + r^o, r^o = R_A r  (eq:rbo)
    Ro(x, t) = Rb(x, t) * Rp(x, t) * Rb(x, t)'   # R^o (transformation rule)
    RoRb(x, t) = Ro(x, t) * Rb(x, t)             # R^o R_b^o
    println("[$pname] |p| at test point = ", round(norm(p(x0, t0)), digits=3), " (WM validity: |p|<4)")

    # numerical constants at the test point
    R0n, Rn, Rpn, upn, ppn = R0(x0), R(x0, t0), Rp(x0, t0), up(x0, t0), pp(x0, t0)
    RAn, Rbn = RA(t0), Rb(x0, t0)

    # ------------------------------------------------------------- 1.2.1  displacement field
    s = "1.2.1"
    ξs = [0.0, 0.31, -0.22]                      # ξ^* = {0,x2,x3}: components of ξ in b^+ (= those in B^*)
    Rhat_truth(x, t, ξ) = rA(t) + RA(t) * (rr(x) + u(x, t) + R(x, t) * R0(x) * ξ)   # eq:R_hat (w=0) mapped to the inertial frame
    C(s, "eq:RboT", "R_b^o = R_0^o R_A^o = R_A R_0", Rbn, (RAn * R0n * RAn') * RAn)
    C(s, "eq:Rbp", "R_b^+ = R_b^oT R_b^o R_b^o = R_b^o", Rbn, Rbn' * Rbn * Rbn)
    C(s, "transf.", "R^+ = R(p^+)  (equivariance of the parametrization, p^+ = R_0^T p)", Rpn, Rfun(R0n' * p(x0, t0)))
    C(s, "eq:RboRp", "R_b^o R^+ = R^o R_b^o", Rbn * Rpn, Ro(x0, t0) * Rbn)
    C(s, "eq:R_hat_0", "R_hat^o = r_b^o + R_b^+ u^+ + (R_b^o R^+) xi^*   [partial]", Rhat_truth(x0, t0, ξs), rbo(x0, t0) + Rbn * upn + (Rbn * Rpn) * ξs)
    C(s, "eq:R_hat_0", "R_hat^o = r_b^o + R_b^+ u^+ + (R^o R_b^o) xi^*   [final]", Rhat_truth(x0, t0, ξs), rbo(x0, t0) + Rbn * upn + (Ro(x0, t0) * Rbn) * ξs)
    C(s, "setup", "r' = b_1 = R_0 e_1 (arclength parametrization of the reference line)", Dx((x, t) -> rr(x), x0, t0), R0n[:, 1])

    # ------------------------------------------------------------- 1.2.2  velocity field
    s = "1.2.2"
    vP_truth = Dt((x, t) -> Rhat_truth(x, t, ξs), x0, t0)                       # d/dt of the position of P
    Rbdot, Rodot = Dt(Rb, x0, t0), Dt(Ro, x0, t0)
    RoRbdot = Dt(RoRb, x0, t0)
    ωtil_o = Rbdot * Rbn'                                                       # eq:omegatil_def (first form)
    ωtil_o2 = Dt(Rb, x0, t0) * Rbn'                                             # (R_b^+ = R_b^o)
    Ωtil_o = RoRbdot * (RoRb(x0, t0))'                                          # eq:Omegatil_def (first form)
    Ωtil_o2 = Dt((x, t) -> Rb(x, t) * Rp(x, t), x0, t0) * (Rbn * Rpn)'          # eq:Omegatil_def (second form)
    C(s, "eq:omegatil_def", "skew-symmetry of Rbdot R_b^T", ωtil_o + ωtil_o', 0 * ωtil_o)
    C(s, "eq:Omegatil_def", "skew-symmetry of (R^o R_b^o)dot (R^o R_b^o)^T", Ωtil_o + Ωtil_o', 0 * Ωtil_o)
    C(s, "eq:Omegatil_def", "(R^o R_b^o)dot(...)^T = (R_b^+ R^+)dot(...)^T", Ωtil_o, Ωtil_o2)
    ξo = RoRb(x0, t0) * ξs                                                      # xi^o = (R^o R_b^o) xi^*
    uo = Rbn * upn                                                              # u^o = R_b^+ u^+
    updot = Dt(up, x0, t0)                                                      # \dot{u}^+ = d(u^+)/dt
    vbo_ad = Dt(rbo, x0, t0)                                                    # \dot r_b^o
    # chain eq:v_P_o, line by line
    L1 = vbo_ad + Dt((x, t) -> Rb(x, t) * up(x, t), x0, t0) + RoRbdot * ξs
    L2 = vbo_ad + Rbdot * upn + Rbn * updot + RoRbdot * RoRb(x0, t0)' * (RoRb(x0, t0) * ξs)
    L3 = vbo_ad + Rbdot * Rbn' * Rbn * upn + Rbn * updot + Ωtil_o * (RoRb(x0, t0) * ξs)
    L4 = vbo_ad + ωtil_o * uo + Rbn * updot + Ωtil_o * ξo
    C(s, "eq:v_P_o", "line 1: vP^o = d/dt[R_hat^o]  (AD of the position) = \\dot r_b + (R_b u^+)dot + (R^o R_b)dot xi^*", vP_truth, L1)
    C(s, "eq:v_P_o", "line 2: product rule", vP_truth, L2)
    C(s, "eq:v_P_o", "line 3: Rbpdot Rbp^T Rbp u^+ and Omega^o (R^o R_b) xi^*", vP_truth, L3)
    C(s, "eq:v_P_o", "line 4: v_b^o + omega~^o u^o + R_b^+ udot^+ + Omega~^o xi^o", vP_truth, L4)

    # relations between velocities of b and of A
    vAo = Dt((x, t) -> rA(t), x0, t0)
    RAdot = ForwardDiff.derivative(RA, t0)
    ωtilAo = RAdot * RAn'                                                       # omega~_A^o := \dot R_A R_A^T
    ωtilA = RAn' * RAdot                                                        # omega~_A   := R_A^T \dot R_A
    ro = RAn * rr(x0)
    C(s, "eq:Rbodot", "Rbdot = RAdot R_0  (R_0 constant in time)", Rbdot, RAdot * R0n)
    C(s, "eq:omtilo_omtilA", "omega~ = R_A^T omega~_b^o R_A = R_A^T Rbdot R_b^T R_A", RAn' * ωtil_o * RAn, RAn' * Rbdot * Rbn' * RAn)
    C(s, "eq:omtilo_omtilA", "... = R_A^T (RAdot R_0)(R_0^T R_A^T) R_A", RAn' * ωtil_o * RAn, RAn' * (RAdot * R0n) * (R0n' * RAn') * RAn)
    C(s, "eq:omtilo_omtilA", "... = R_A^T RAdot = omega~_A", RAn' * ωtil_o * RAn, ωtilA)
    C(s, "eq:omegatil_def", "omega~_b^o = RAdot R_A^T = omega~_A^o (free vector: omega_b = omega_A)", ωtil_o, ωtilAo)
    C(s, "eq:vbo", "v_b^o = d/dt[r_A^o + R_A r] = v_A^o + omega~_A^o r^o", vbo_ad, vAo + ωtilAo * ro)
    vb = RAn' * vbo_ad                                                          # components of v_b in A
    C(s, "eq:vb", "v_b = R_A^T v_b^o = R_A^T (v_A^o + omega~_A^o R_A r)", vb, RAn' * (vAo + ωtilAo * RAn * rr(x0)))
    C(s, "eq:vb", "(equivalent physical form) v_b = v_A + omega_A x r  (A components)", vb, RAn' * vAo + ωtilA * rr(x0))
    # V^o, V^*
    Vo = Dt((x, t) -> rbo(x, t) + Rb(x, t) * up(x, t), x0, t0)                  # velocity of point B (origin of B^*)
    C(s, "eq:V_o", "V^o = v_b^o + omega~^o u^o + R_b^+ udot^+  (velocity of point B)", Vo, vbo_ad + ωtil_o * uo + Rbn * updot)
    ωp = R0n' * axial(ωtilA)                                                    # omega^+ = R_0^T omega
    ωtilp = tilde(ωp)
    vbp = R0n' * vb                                                             # v_b^+ = R_0^T v_b
    C(s, "transf.", "omega~^+ = R_b^T omega~^o R_b = R_0^T omega~_A R_0", ωtilp, Rbn' * ωtil_o * Rbn)
    Vs1 = (Rbn * Rpn)' * Vo
    Vs2 = Rpn' * Rbn' * (vbo_ad + ωtil_o * uo + Rbn * updot)
    Vs3 = Rpn' * (vbp + ωtilp * upn + updot)
    C(s, "eq:V_star_initial", "V^* = (R_b^+ R^+)^T V^o", Vs1, Vs2)
    C(s, "eq:V_star", "V^* = R^+T (v_b^+ + omega~^+ u^+ + udot^+)", Vs1, Vs3)
    # Omega~^*
    Ωtil_s1 = (Rbn * Rpn)' * Ωtil_o * (Rbn * Rpn)                               # eq:Omegatil_initial
    Rpdot = Dt(Rp, x0, t0)
    Ωtil_s2 = Rpn' * Rbn' * (Rbdot * Rpn + Rbn * Rpdot)
    Ωtil_s3 = Rpn' * ωtilp * Rpn + Rpn' * Rpdot
    Ωtil_s4 = -Rpdot' * Rpn + Rpn' * ωtilp * Rpn
    C(s, "eq:Omegatil_star_partial", "R_b^T Rbdot = omega~^+ (R_b^+T Rbpdot)", ωtilp, Rbn' * Rbdot)
    C(s, "eq:Omegatil_initial", "Omega~^* = (R_b R^+)^T (R_b R^+)dot ... (2 forms)", Ωtil_s1, Ωtil_s2)
    C(s, "eq:Omegatil_star_partial", "Omega~^* = R^+T omega~^+ R^+ + R^+T Rpdot", Ωtil_s1, Ωtil_s3)
    C(s, "eq:Omegatil_star", "skew-symmetry of R^+T Rpdot", Rpn' * Rpdot + (Rpn' * Rpdot)', 0 * Rpn)
    C(s, "eq:Omegatil_star", "Omega~^* = -Rpdot^T R^+ + R^+T omega~^+ R^+", Ωtil_s1, Ωtil_s4)
    # Omega^* in terms of rotation parameters
    Hpn = H_num(Rfun, ppn)                                                      # H^+ = H(p^+)
    Hn = H_num(Rfun, p(x0, t0))                                                 # H in A
    ppdot = Dt(pp, x0, t0)                                                      # \dot p^+ = d(p^+)/dt
    pdot = Dt(p, x0, t0)                                                        # \dot p (A)
    Oms = axial(Ωtil_s1)
    C(s, "transf.", "H^+ = H(p^+) = R_0^T H R_0", Hpn, R0n' * Hn * R0n)
    C(s, "eq:Omega_star_rotparp", "axial(R^+T Rpdot) = H^+T pdot^+", axial(Rpn' * Rpdot), Hpn' * ppdot)
    C(s, "eq:Omega_star_rotparp", "Omega^* = H^+T pdot^+ + R^+T omega^+", Oms, Hpn' * ppdot + Rpn' * ωp)
    C(s, "eq:Omega_star_rotpar", "Omega^* = R_0^T (H^T pdot + R^T omega)", Oms, R0n' * (Hn' * pdot + Rn' * axial(ωtilA)))
    C(s, "eq:Omega_def", "Omega^o = axial((R^o R_b^o)dot(..)^T)  =>  Omega^* = (R_b R^+)^T Omega^o", Oms, (Rbn * Rpn)' * axial(Ωtil_o))
    # v_P^*
    vPs_def = (Rbn * Rpn)' * vP_truth                                           # definition: v_P^* = (R_b^+ R^+)^T v_P^o
    C(s, "eq:v_P_star", "v_P^* = V^* + Omega~^* xi^*", vPs_def, Vs1 + Ωtil_s1 * ξs)
    asWritten = Vs1 + (Rbn * Rpn)' * Ωtil_o * (Rbn * Rpn) * (Rbn * Rpn)' * ξs   # line 806 of main.tex: ends with (R_b R^+)^T xi^*
    corrected = Vs1 + (Rbn * Rpn)' * Ωtil_o * (Rbn * Rpn) * (Rbn * Rpn)' * ξo   # ... should end with (R_b R^+)^T xi^o
    C(s, "eq:v_P_star", "AS WRITTEN (3rd line): ...(R_b R^+)(R_b R^+)^T xi^*  [typo: xi^*]", vPs_def, asWritten; kind=:typo)
    C(s, "eq:v_P_star", "corrected: ...(R_b R^+)(R_b R^+)^T xi^o", vPs_def, corrected)

    # ------------------------------------------------------------- 1.3  strain measures
    s = "1.3.1"
    R0p = Dx((x, t) -> R0(x), x0, t0)                                           # R_0'
    ktil = R0p * R0n'                                                           # eq:k_def:  k~ = R_0' R_0^T
    k = axial(ktil)
    ktilp = R0n' * ktil * R0n
    C(s, "eq:ktilp", "k~^+ = R_0^T k~ R_0 = R_0^T R_0'", ktilp, R0n' * R0p)
    C(s, "eq:ktilp", "skew-symmetry of R_0^T R_0'", ktilp + ktilp', 0 * ktilp)
    kp = R0n' * k
    C(s, "transf.", "k^+ = R_0^T k  and  k~^+ = tilde(k^+)", ktilp, tilde(kp))
    rprime = Dx((x, t) -> rr(x), x0, t0)
    rpprime = Dx((x, t) -> rp(x), x0, t0)                                       # r'^+ = d(r^+)/dx1
    C(s, "eq:Rztrprime", "R_0^T r' = R_0^T (R_0 r^+)' = r'^+ + k~^+ r^+ (line 2: R_0^T R_0' r^+ + r'^+)", R0n' * rprime, R0n' * R0p * rp(x0) + rpprime)
    C(s, "eq:Rztrprime", "R_0^T r' = r'^+ + k~^+ r^+", R0n' * rprime, rpprime + ktilp * rp(x0))
    C(s, "eq:Rztrprime_simple", "R_0^T r' = R_0^T b_1 = b_1^+ = e_1", R0n' * rprime, e1)
    upprime = Dx(up, x0, t0)                                                    # u'^+ = d(u^+)/dx1
    uprime = Dx(u, x0, t0)
    C(s, "eq:Rztruprime", "R_0^T (r+u)' = (r'^+ + u'^+) + k~^+ (r^+ + u^+)", R0n' * (rprime + uprime), (rpprime + upprime) + ktilp * (rp(x0) + upn))
    C(s, "eq:Rztruprime", "... = e_1 + u'^+ + k~^+ u^+", R0n' * (rprime + uprime), e1 + upprime + ktilp * upn)
    C(s, "transf.", "NOTE: (u')^+ = R_0^T u' is NOT u'^+ = d(u^+)/dx1 (the text uses the latter)", R0n' * uprime, upprime; kind=:note)
    # gamma
    γ = Rn' * (rprime + uprime) - rprime                                        # eq:gamma_def
    γp = R0n' * γ
    C(s, "eq:gamp_A", "gamma^+ = R_0^T R^T (R_0 e_1 + u') - e_1", γp, R0n' * Rn' * (R0n * e1 + uprime) - e1)
    C(s, "eq:gamp_A", "... = R_0^T R^T R_0 R_0^T (r'+u') - e_1", γp, R0n' * Rn' * R0n * R0n' * (rprime + uprime) - e1)
    C(s, "eq:gamp_b", "gamma^+ = R^+T (b_1^+ + u'^+ + k~^+ u^+) - b_1^+", γp, Rpn' * (e1 + upprime + ktilp * upn) - e1)
    C(s, "eq:gamp_b_rearranged", "e_1 + u'^+ + k~^+ u^+ = R^+ (e_1 + gamma^+)", e1 + upprime + ktilp * upn, Rpn * (e1 + γp))
    s = "1.3.2"
    Rx = Dx(R, x0, t0)                                                          # R' (components in A)
    RRz(x, t) = R(x, t) * R0(x)
    Ktil = Dx(RRz, x0, t0) * (RRz(x0, t0))'                                     # eq:K_def: K~ = (R R_0)'(R R_0)^T
    K = axial(Ktil)
    Ktil2 = Dx((x, t) -> R0(x) * Rp(x, t), x0, t0) * (R0n * Rpn)'               # second form: (R_0 R^+)'(R_0 R^+)^T
    C(s, "eq:K_def", "(R R_0)'(R R_0)^T = (R_0 R^+)'(R_0 R^+)^T", Ktil, Ktil2)
    C(s, "eq:K_def", "skew-symmetry of K~", Ktil + Ktil', 0 * Ktil)
    κ = Rn' * K - k                                                             # eq:kappa_def
    κp = R0n' * κ
    Ks = (R0n * Rpn)' * K                                                       # K^* = (R R_0)^T K
    C(s, "eq:kappa_plus", "kappa^+ = R_0^T R^T K - R_0^T k", κp, R0n' * Rn' * K - R0n' * k)
    C(s, "eq:kappa_plus", "kappa^+ = (R R_0)^T K - k^+ = K^* - k^+", κp, Ks - kp)
    Rpx = Dx(Rp, x0, t0)                                                        # R'^+ = d(R^+)/dx1
    Ktils1 = (R0n * Rpn)' * Ktil * (R0n * Rpn)
    C(s, "eq:Ktils1", "K~^* = (R_0 R^+)^T (R_0 R^+)'", Ktils1, (R0n * Rpn)' * (R0p * Rpn + R0n * Rpx))
    C(s, "eq:Ktils1", "K~^* = R^+T k~^+ R^+ + R^+T R'^+", Ktils1, Rpn' * ktilp * Rpn + Rpn' * Rpx)
    C(s, "eq:Ktils", "K~^* = -R'^+T R^+ + R^+T k~^+ R^+", Ktils1, -Rpx' * Rpn + Rpn' * ktilp * Rpn)
    C(s, "transf.", "K^* = axial(K~^*)", Ks, axial(Ktils1))
    ppx = Dx(pp, x0, t0)                                                        # p'^+ = d(p^+)/dx1
    px = Dx(p, x0, t0)                                                          # p' (A)
    C(s, "eq:Ks_rotpar", "axial(R^+T R'^+) = H^+T p'^+", axial(Rpn' * Rpx), Hpn' * ppx)
    C(s, "eq:Ks_rotpar", "K^* = H^+T p'^+ + R^+T k^+", Ks, Hpn' * ppx + Rpn' * kp)
    C(s, "eq:kappa_plus2", "kappa^+ = H^+T p'^+ + R^+T k^+ - k^+", κp, Hpn' * ppx + Rpn' * kp - kp)
    κtil_p = Ktils1 - ktilp
    C(s, "eq:kappatil_plus", "kappa~^+ = K~^* - k~^+ = R_0^T R^T R' R_0", κtil_p, R0n' * Rn' * Rx * R0n)
    C(s, "eq:kappatil_plus", "kappa~^+ = tilde(kappa^+)", κtil_p, tilde(κp))
    C(s, "eq:kappa_plus_A", "kappa^+ = R_0^T axial(R^T R')", κp, R0n' * axial(Rn' * Rx))
    C(s, "eq:kappa_plus_A", "kappa^+ = R_0^T H^T p'  (H in A)", κp, R0n' * Hn' * px)
    C(s, "transf.", "R_0 p'^+ = p' + R_0 p~^+ k^+   (relation between p' and p'^+)", R0n * ppx, px + R0n * tilde(ppn) * kp)
    # --- independent (physical) test: a rigid-body motion of the whole beam, u = (Q-I) r + c, R = Q, must give zero strains
    pq = [0.3, -0.5, 0.4]; Q = Rfun(pq); cc = [0.2, 0.1, -0.3]
    up_rb(x) = R0(x)' * ((Q - I) * rr(x) + cc)
    ppq(x) = R0(x)' * pq
    Rp_rb = R0n' * Q * R0n
    γp_rb = Rp_rb' * (e1 + Dx((x, t) -> up_rb(x), x0, t0) + ktilp * up_rb(x0)) - e1
    κp_rb = H_num(Rfun, ppq(x0))' * Dx((x, t) -> ppq(x), x0, t0) + Rp_rb' * kp - kp
    C(s, "eq:gamp_b", "rigid-body motion of the whole beam gives gamma^+ = 0 (objectivity of the force-strains)", γp_rb, zeros(3))
    C(s, "eq:kappa_plus2", "rigid-body motion of the whole beam gives kappa^+ = 0 (needs R^T - I = -H^T p~)", κp_rb, zeros(3))
    # --- conceptual notes about the text
    C(s, "eq:gamma_def", "NOTE: |(r+u)'| = ds/dx1 (= 1 + stretch) is not 1, so (r+u)' is NOT a unit tangent vector", [norm(rprime + uprime)], [1.0]; kind=:note)

    # ------------------------------------------------------------- 1.4.1  strain energy (trivial, quadratic U)
    s = "1.4.1"
    rng = MersenneTwister(7)
    Sm = let A = randn(rng, 6, 6); A * A' + 6I end                              # random SPD 6x6 sectional stiffness
    εv = vcat(γp, κp)
    U(ε) = 0.5 * dot(ε, Sm * ε)
    gradU = ForwardDiff.gradient(U, εv)
    C(s, "eq:force_res/mom_res", "{F^*;M^*} = (dU/d{gamma^+;kappa^+})^T = S^* {gamma^+;kappa^+}", gradU, Sm * εv)
    C(s, "eq:constitutive", "{gamma^+;kappa^+} = C^* {F^*;M^*},  C^* = S^{*-1}", εv, inv(Sm) * gradU)
    δε = randn(rng, 6)
    C(s, "eq:delU", "dU = dgamma^T F + dkappa^T M  (chain rule)", ForwardDiff.derivative(h -> U(εv + h * δε), 0.0), dot(δε[1:3], gradU[1:3]) + dot(δε[4:6], gradU[4:6]))
    C(s, "Sec. 1.3.1/1.3.2", "NOTE: the scalar product gamma^T F + kappa^T M is 2U (linear case), not U: conjugacy holds for the variations (eq:delU)", [dot(εv, gradU)], [U(εv)]; kind=:note)

    # ------------------------------------------------------------- 1.4.2  kinetic energy
    s = "1.4.2"
    nP = 40
    x2 = 0.6 .* (rand(rng, nP) .- 0.5); x3 = 0.2 .* (rand(rng, nP) .- 0.5); ρdS = 0.5 .+ rand(rng, nP)   # point-mass quadrature of the section
    ξi = [[0.0, x2[i], x3[i]] for i in 1:nP]
    V, Ω = randn(rng, 3), randn(rng, 3)
    μ = sum(ρdS)                                                                # eq:sec:mass
    η = sum(ρdS[i] * ξi[i] for i in 1:nP) / μ                                   # eq:sec:mom_inertia_1
    ι = sum(ρdS[i] * tilde(ξi[i]) * tilde(ξi[i])' for i in 1:nP)                # eq:sec:mom_inertia_2
    C(s, "eq:sec:mom_inertia_2", "iota^* symmetric", ι, ι')
    C(s, "eq:sec:mom_inertia_2", "iota^* = int rho (|xi|^2 I - xi xi^T) dS", ι, sum(ρdS[i] * (dot(ξi[i], ξi[i]) * I - ξi[i] * ξi[i]') for i in 1:nP))
    K_direct = 0.5 * sum(ρdS[i] * dot(V + tilde(Ω) * ξi[i], V + tilde(Ω) * ξi[i]) for i in 1:nP)
    K_l2 = 0.5 * sum(ρdS[i] * (dot(V, V) + dot(V, tilde(Ω) * ξi[i]) + dot(tilde(Ω) * ξi[i], V) + dot(tilde(Ω) * ξi[i], tilde(Ω) * ξi[i])) for i in 1:nP)
    K_l3 = 0.5 * sum(ρdS[i] * dot(V, V) + 2ρdS[i] * dot(V, tilde(Ω) * ξi[i]) + ρdS[i] * dot(tilde(ξi[i])' * Ω, tilde(ξi[i])' * Ω) for i in 1:nP)
    K_l4 = 0.5 * sum(ρdS[i] * dot(V, V) + 2ρdS[i] * Ω' * tilde(V)' * ξi[i] + Ω' * (ρdS[i] * tilde(ξi[i]) * tilde(ξi[i])') * Ω for i in 1:nP)
    K_l5 = 0.5 * sum(ρdS[i] * dot(V, V) - 2 * Ω' * tilde(V) * (ρdS[i] * ξi[i]) + Ω' * (ρdS[i] * tilde(ξi[i]) * tilde(ξi[i])') * Ω for i in 1:nP)
    K_simpl = 0.5 * (μ * dot(V, V) - 2μ * Ω' * tilde(V) * η + Ω' * ι * Ω)
    C(s, "eq:kinetic_energy_expanded", "line 2 (expanded square)", [K_direct], [K_l2])
    C(s, "eq:kinetic_energy_expanded", "line 3 ((Omega~ xi)^T(Omega~ xi) = (xi~^T Omega)^T(xi~^T Omega))", [K_direct], [K_l3])
    C(s, "eq:kinetic_energy_expanded", "line 4 (Omega^T V~^T xi, Omega^T rho xi~ xi~^T Omega)", [K_direct], [K_l4])
    C(s, "eq:kinetic_energy_expanded", "line 5 (-2 Omega^T V~ rho xi)", [K_direct], [K_l5])
    C(s, "eq:kinetic_energy_simplified", "K = 1/2 (mu V^T V - 2 mu Omega^T V~ eta + Omega^T iota Omega)", [K_direct], [K_simpl])
    Kfun(w) = 0.5 * sum(ρdS[i] * dot(w[1:3] + tilde(w[4:6]) * ξi[i], w[1:3] + tilde(w[4:6]) * ξi[i]) for i in 1:nP)
    g = ForwardDiff.gradient(Kfun, vcat(V, Ω))                                   # (dK/dV*, dK/dOmega*) from the DIRECT sum
    Ps, Hs = g[1:3], g[4:6]
    C(s, "eq:lin_mom_expanded", "P^* = dK/dV^{*T} = mu (V^* + eta~^T Omega^*)", Ps, μ * (V + tilde(η)' * Ω))
    C(s, "eq:ang_mom_expanded", "H^* = dK/dOmega^{*T} = mu eta~ V^* + iota Omega^*", Hs, μ * tilde(η) * V + ι * Ω)
    Mm = [μ * I(3) μ * tilde(η)'; μ * tilde(η) ι]
    C(s, "eq:mom_vel / eq:M", "{P^*;H^*} = M^* {V^*;Omega^*},  M^* = [mu I, mu eta~^T; mu eta~, iota]", vcat(Ps, Hs), Mm * vcat(V, Ω))
    C(s, "eq:M", "M^* symmetric and positive semi-definite", [minimum(eigvals(Symmetric(Mm))) >= -1e-12 ? 0.0 : 1.0], [0.0])
    C(s, "eq:M", "M^* symmetric", Mm, Mm')
    # independent (physical) meaning: P^* = int rho v_P dS  and  H^* = int rho xi x v_P dS  (angular momentum about B)
    vPi = [V + tilde(Ω) * ξi[i] for i in 1:nP]
    C(s, "eq:lin_mom", "P^* = int rho v_P^* dS (linear momentum per unit length)", Ps, sum(ρdS[i] * vPi[i] for i in 1:nP))
    C(s, "eq:ang_mom", "H^* = int rho xi~ v_P^* dS (angular momentum per unit length about B)", Hs, sum(ρdS[i] * tilde(ξi[i]) * vPi[i] for i in 1:nP))
    C(s, "Sec. 1.4.2 (text)", "iota^*: with xi_1 = 0 only torsion, two bending-rotary and ONE product of inertia (1,2),(1,3) = 0; no translation coupling", ι[1, 2:3], zeros(2))
    C(s, "Sec. 1.4.2 (text)", "iota^*(2,3) = -int rho x2 x3 dS (coupling of the two bending rotations)", [ι[2, 3]], [-sum(ρdS .* x2 .* x3)])
    # also with the actual V^*, Omega^* of the kinematic section, and the formula of eq:v_P_star
    Kreal = 0.5 * sum(ρdS[i] * dot(Vs1 + Ωtil_s1 * ξi[i], Vs1 + Ωtil_s1 * ξi[i]) for i in 1:nP)
    Kinertial = 0.5 * sum(ρdS[i] * dot(Dt((x, t) -> Rhat_truth(x, t, ξi[i]), x0, t0), Dt((x, t) -> Rhat_truth(x, t, ξi[i]), x0, t0)) for i in 1:nP)
    C(s, "eq:kinetic_energy", "K from V^*, Omega^* of the kinematics = K from the inertial velocity of every material point", [Kreal], [Kinertial])

    # ------------------------------------------------------------- 1.4.3  virtual work and virtual action
    s = "1.4.3"
    ηu = Fld(61; scale=0.5); ηp = Fld(71; scale=0.5)                             # arbitrary variations of u and p (A components)
    δu = ηu(x0, t0)
    δp = ηp(x0, t0)
    δR = ForwardDiff.derivative(ε -> Rfun(p(x0, t0) + ε * δp), 0.0)              # delta R (A components)
    δRp = ForwardDiff.derivative(ε -> Rfun(R0n' * (p(x0, t0) + ε * δp)), 0.0)    # delta R^+
    δψtil_s = Rpn' * δRp                                                         # eq:del_psitils:  R^+T delta R^+
    δψs = axial(δψtil_s)                                                         # eq:del_psis
    C(s, "eq:del_psitils", "R^+T dR^+ is skew-symmetric", δψtil_s + δψtil_s', 0 * δψtil_s)
    C(s, "eq:del_psitils", "R^+T dR^+ = -(dR^+T R^+)", δψtil_s, -(δRp' * Rpn))
    C(s, "eq:del_Rp", "delta R^+ = R^+ delta psi~^*", δRp, Rpn * δψtil_s)
    C(s, "eq:del_RpT", "delta R^+T = -delta psi~^* R^+T", δRp', -δψtil_s * Rpn')
    δuplus = ForwardDiff.derivative(ε -> R0n' * (u(x0, t0) + ε * δu), 0.0)       # delta u^+
    δqp = δuplus                                                                 # eq:delqp
    δqs = Rpn' * δqp                                                             # eq:delqs
    δq = R0n * δqp                                                               # eq:delq
    C(s, "eq:deltaq_expand", "delta q = R_0 delta(R_0^T u) = R_0 R_0^T delta u = delta u  (delta R_0 = 0)", δq, δu)
    C(s, "eq:delq", "delta q = R_0 R^+ delta q^*", δq, R0n * Rpn * δqs)
    C(s, "eq:deltau", "delta u^T = delta q^{*T} (R_0 R^+)^T", δu, R0n * Rpn * δqs)
    δψ_spatial = axial(δR * Rn')                                                 # virtual rotation (A components) = axial(dR R^T)
    C(s, "eq:delpsiT_f_m", "delta psi (A) = R_0 R^+ delta psi^*  [spatial virtual rotation = pushed-forward body one]", δψ_spatial, R0n * Rpn * δψs)
    C(s, "eq:delpsiT_f_m", "R_0 R^+ = R R_0 (composition)", R0n * Rpn, Rn * R0n)
    fA, mA = randn(rng, 3), randn(rng, 3)
    fs, ms_ = (R0n * Rpn)' * fA, (R0n * Rpn)' * mA                               # f^* = (R_0R^+)^T f, m^* = ...
    C(s, "eq:delpsiT_f_m", "f = (R_0 R^+) f^*,  m = (R_0 R^+) m^*", fA, R0n * Rpn * fs)
    C(s, "eq:deltaWs", "dW = du^T f + dpsi^T m  =  dq^{*T} f^* + dpsi^{*T} m^*", [dot(δu, fA) + dot(δψ_spatial, mA)], [dot(δqs, fs) + dot(δψs, ms_)])

    # ------------------------------------------------------------- extra: App. A pointwise variations (they turn eq:delU, eq:delK into eq:delU_final, eq:delK_final)
    s = "A.1/A.2"
    uε(x, t, ε) = u(x, t) + ε * ηu(x, t)
    pε(x, t, ε) = p(x, t) + ε * ηp(x, t)
    Rpε(x, t, ε) = Rfun(R0(x)' * pε(x, t, ε))
    upε(x, t, ε) = R0(x)' * uε(x, t, ε)
    δRp_f(x, t) = ForwardDiff.derivative(ε -> Rpε(x, t, ε), 0.0)
    δψs_f(x, t) = axial(Rp(x, t)' * δRp_f(x, t))            # virtual rotation delta psi^*(x,t)   (eq:del_psis)
    δqs_f(x, t) = Rp(x, t)' * (R0(x)' * ηu(x, t))           # virtual displacement delta q^*(x,t) (eq:delqs)
    ktilp_f(x) = R0(x)' * ForwardDiff.derivative(R0, x)
    γp_f(x, t, ε) = Rpε(x, t, ε)' * (e1 + ForwardDiff.derivative(sx -> upε(sx, t, ε), x) + ktilp_f(x) * upε(x, t, ε)) - e1   # eq:gamp_b
    function κp_f(x, t, ε)                                   # eq:kappa_plus_A
        Rε(sx) = Rfun(pε(sx, t, ε))
        R0(x)' * axial(Rε(x)' * ForwardDiff.derivative(Rε, x))
    end
    Ωs_f(x, t, ε) = axial(Rpε(x, t, ε)' * ForwardDiff.derivative(st -> Rpε(x, st, ε), t)) + Rpε(x, t, ε)' * ωp                     # eq:Omega_star_rotparp
    Vs_f(x, t, ε) = Rpε(x, t, ε)' * (vbp + ωtilp * upε(x, t, ε) + ForwardDiff.derivative(st -> upε(x, st, ε), t))              # eq:V_star
    Dε(f) = ForwardDiff.derivative(ε -> f(ε), 0.0)
    δγp = Dε(ε -> γp_f(x0, t0, ε)); δκp = Dε(ε -> κp_f(x0, t0, ε)); δΩs = Dε(ε -> Ωs_f(x0, t0, ε)); δVs = Dε(ε -> Vs_f(x0, t0, ε))
    C(s, "setup", "perturbed fields reproduce the unperturbed ones at epsilon = 0 (gamma^+, kappa^+, Omega^*, V^*)", vcat(γp_f(x0, t0, 0.0), κp_f(x0, t0, 0.0), Ωs_f(x0, t0, 0.0), Vs_f(x0, t0, 0.0)), vcat(γp, κp, Oms, Vs1))
    C(s, "eq:del_kappa_p", "delta kappa^+ = delta psi'^* + K~^* delta psi^*", δκp, Dx(δψs_f, x0, t0) + Ktils1 * δψs_f(x0, t0))
    C(s, "eq:delgamp", "delta gamma^+ = delta q'^* + K~^* delta q^* + (a~_1 + gamma~^+) delta psi^*", δγp, Dx(δqs_f, x0, t0) + Ktils1 * δqs_f(x0, t0) + (tilde(e1) + tilde(γp)) * δψs_f(x0, t0))
    C(s, "eq:del_OmsT", "delta Omega^* = delta psi_dot^* + Omega~^* delta psi^*  (<=> eq:del_OmsT)", δΩs, Dt(δψs_f, x0, t0) + Ωtil_s1 * δψs_f(x0, t0))
    C(s, "eq:delVsT", "delta V^* = delta q_dot^* + Omega~^* delta q^* + V~^* delta psi^*  (<=> eq:delVsT)", δVs, Dt(δqs_f, x0, t0) + Ωtil_s1 * δqs_f(x0, t0) + tilde(Vs1) * δψs_f(x0, t0))
    return nothing
end

# =============================================================================
# Extra: facts about the rotation parametrization used by Chapter 1 (properties of H) and the closed
# forms of Chapter 5 for the Wiener-Milenkovic set (outside the requested scope, cheap to check)
function parametrization_facts(Rfun, pname; wm=false)
    C(sec, label, desc, a, b; kw...) = check(pname, sec, label, desc, a, b; kw...)
    s = "param."
    rng = MersenneTwister(99)
    for trial in 1:3
        p = 1.1 * (rand(rng, 3) .- 0.5) .* 2.5
        R = Rfun(p); H = H_num(Rfun, p)
        C(s, "def. of H", "R orthogonal, det = +1 (trial $trial)", R * R', Matrix(1.0I, 3, 3))
        C(s, "(1.2.2)", "H^T = R^T H  (so that axial(R^T Rdot) = H^T pdot) (trial $trial)", H', R' * H)
        C(s, "(App. A, ptil H = R - I)", "p~ H = R - I  and  R^T - I = -H^T p~ (trial $trial)", tilde(p) * H, R - I)
        C(s, "(App. A)", "R^T - I = -H^T p~ (trial $trial)", R' - I, -H' * tilde(p))
        Q = R_exp([0.4, -0.7, 0.2])
        C(s, "(1.2.2)", "equivariance: R(Q p) = Q R(p) Q^T and H(Q p) = Q H(p) Q^T (trial $trial)", vcat(Rfun(Q * p), H_num(Rfun, Q * p)), vcat(Q * R * Q', Q * H * Q'))
        if wm
            p0 = 2 - dot(p, p) / 8; υ = 1 / (4 - p0)
            HT_doc = 2υ^2 * (p0 * I + (1 / 4) * p * p' - tilde(p))               # eq:HT_WM
            HinvT_doc = (1 / 2) * (p0 * I + (1 / 4) * p * p' + tilde(p))         # eq:HinvT_WM
            C("5.x", "eq:HT_WM", "H^T = 2 upsilon^2 (p0 I + pp^T/4 - p~)  (trial $trial)", H', HT_doc)
            C("5.x", "eq:HinvT_WM", "H^{-T} = 1/2 (p0 I + pp^T/4 + p~)  (trial $trial)", inv(H)', HinvT_doc)
            R_doc = υ^2 * ((1 / υ^2 - 2 * dot(p, p)) * I + 2 * p * p' + 2 * p0 * tilde(p))
            C("5.x", "eq:R_WM", "R = upsilon^2[(1/upsilon^2 - 2|p|^2) I + 2pp^T + 2 p0 p~]  (trial $trial)", R, R_doc)
            # eq:p0 AS WRITTEN: p0 = 2 - |p|/8 (instead of |p|^2/8)
            p0w = 2 - norm(p) / 8; υw = 1 / (4 - p0w)
            Rw = υw^2 * ((1 / υw^2 - 2 * dot(p, p)) * I + 2 * p * p' + 2 * p0w * tilde(p))
            C("5.x", "eq:p0", "AS WRITTEN: p0 = 2 - |p|/8 gives an orthogonal R? (trial $trial)", Rw * Rw', Matrix(1.0I, 3, 3); kind=:typo)
        end
    end
end

# =============================================================================
println("="^100)
println("Verification of Chapter 1 (sections 1.1 - 1.4.3) of theory/main.tex   -   tolerance = $TOL")
println("="^100)
run_checks(R_WM, "WM ")
parametrization_facts(R_WM, "WM "; wm=true)
run_checks(R_exp, "exp")
parametrization_facts(R_exp, "exp")

# =============================================================================
# Section 1.5 (finite element model): algebra of the discrete weak form, element geometry,
# solver identities, hinge constraint
# =============================================================================
function gaussint(g, a, b)                         # ∫_a^b g(s) ds, 40-point Gauss (exact for smooth g)
    z, w = QuadGK.gauss(Float64, 40)
    sum((b - a) / 2 * w[i] * g((a + b) / 2 + (b - a) / 2 * z[i]) for i in eachindex(w))
end

function section15_checks()
    C(sec, label, desc, a, b; kw...) = check("FE  ", sec, label, desc, a, b; kw...)
    rng = MersenneTwister(15)
    v3() = randn(rng, 3)
    M3() = randn(rng, 3, 3) + 3I
    gz, gw = QuadGK.gauss(Float64, 3)
    ζg, wg = (gz .+ 1) ./ 2, gw ./ 2
    quad01(g) = sum(wg[i] * g(ζg[i]) for i in eachindex(wg))   # ∫_0^1 g dζ (exact up to degree 5)

    # ------------------------------------------------------------- 1.5.1 weak form: eq:motion_discretized -> eq:motion_discretized3
    for e in 1:4
        Δℓ = 0.1 + rand(rng)
        Re, Rze, He = M3(), M3(), M3()
        HinvT = inv(He)'
        Pse, Fse, Mse, Vse, Ose, Hse = ntuple(_ -> v3(), 6)
        ue, pe, vbe, ombe, udote, pdote, κe, γe = ntuple(_ -> v3(), 8)
        a1 = [1.0, 0.0, 0.0]
        fa, fb, ma, mb, dotPb, dotHb = ntuple(_ -> v3(), 6)            # f(ζ)=fa+fb ζ, m(ζ)=ma+mb ζ; time derivatives as symbols
        du0, du1, dp0, dp1, dF0, dF1, dM0, dM1 = ntuple(_ -> v3(), 8)  # nodal variations at n=e and n+1
        dP, dH = v3(), v3()                                            # elemental variations (constant over the element)
        A = Re * Rze
        AF, AM, Pb, Hb = A * Fse, A * Mse, A * Pse, A * Hse
        bP = dotPb + tilde(ombe) * Pb
        bH = dotHb + tilde(ombe) * Hb + A * (tilde(Vse) * Pse - tilde(a1 + γe) * Fse)
        Qk = A * (a1 + γe) - Rze * a1
        Ks = HinvT * Rze * κe
        f(ζ) = fa + fb * ζ
        m(ζ) = ma + mb * ζ
        integrand(ζ) = (dot(-du0 + du1, AF) + dot(-dp0 + dp1, AM) - dot(-dF0 + dF1, ue) - dot(-dM0 + dM1, pe)) / Δℓ +
                       dot((1 - ζ) * du0 + ζ * du1, bP - f(ζ)) + dot((1 - ζ) * dp0 + ζ * dp1, bH - m(ζ)) +
                       - dot((1 - ζ) * dF0 + ζ * dF1, Qk) - dot((1 - ζ) * dM0 + ζ * dM1, Ks) +
                       dot(dP, A * Vse - vbe - tilde(ombe) * ue - udote) + dot(dH, A * Ose - ombe - Re * He' * pdote)
        LHS = Δℓ * quad01(integrand)
        vint(g) = Δℓ * quad01(g)
        f_left, f_right = vint(ζ -> (1 - ζ) * f(ζ)), vint(ζ -> ζ * f(ζ))
        m_left, m_right = vint(ζ -> (1 - ζ) * m(ζ)), vint(ζ -> ζ * m(ζ))
        FUm = -AF - f_left + Δℓ / 2 * bP
        FUp = AF - f_right + Δℓ / 2 * bP
        FPm = -AM - m_left + Δℓ / 2 * bH
        FPp = AM - m_right + Δℓ / 2 * bH
        FFm, FFp = ue - Δℓ / 2 * Qk, -ue - Δℓ / 2 * Qk
        FMm, FMp = pe - Δℓ / 2 * Ks, -pe - Δℓ / 2 * Ks
        FV = Δℓ * (A * Vse - vbe - tilde(ombe) * ue - udote)
        FO = Δℓ * (A * Ose - ombe - Re * He' * pdote)
        RHS = dot(du0, FUm) + dot(du1, FUp) + dot(dp0, FPm) + dot(dp1, FPp) + dot(dF0, FFm) + dot(dF1, FFp) +
              dot(dM0, FMm) + dot(dM1, FMp) + dot(dP, FV) + dot(dH, FO)
        C("1.5.1", "eq:motion_discretized3 e=$e", "element weak form = coefficient form (Fum..FO), random data", [LHS], [RHS])
    end
    C("1.5.1", "eq:integrals_zeta", "∫_0^1 (1-ζ) Δℓ dζ = Δℓ/2", [0.7 * quad01(ζ -> 1 - ζ)], [0.35])

    # ------------------------------------------------------------- 1.5.1 element geometry (App. C)
    kp = [0.3, -0.2, 0.5]                                  # k1 = initial twist, (k2,k3) = initial curvature
    R00 = R_exp(v3())
    κ = norm(kp)
    kt, kk = tilde(kp), kp * kp'
    Rz(x) = R00 * (I + sin(κ * x) * kt / κ + (1 - cos(κ * x)) * kt * kt / κ^2)   # Rodrigues (AD-friendly)
    rq(x) = R00 * gaussint(s -> exp(kt * s) * [1.0, 0.0, 0.0], 0.0, x)     # r(x) - r(0)
    C("1.5.1", "eq:Rzp", "Rz' = Rz k~^+ (constant k^+)", ForwardDiff.derivative(Rz, 0.3), Rz(0.3) * kt)
    C("1.5.1", "eq:etoktilpx1", "Rodrigues formula = exp map", I + sin(κ * 0.4) * kt / κ + (1 - cos(κ * 0.4)) * kt * kt / κ^2, exp(kt * 0.4))
    C("1.5.1", "eq:etoktilpx1", "third form of Rodrigues (kk^T/κ²)", I + sin(κ * 0.4) * kt / κ + (1 - cos(κ * 0.4)) * kt * kt / κ^2,
      (I - kk / κ^2) * cos(κ * 0.4) + kt / κ * sin(κ * 0.4) + kk / κ^2)
    Lr = 0.9
    eps_closed(x) = (I - kk / κ^2) * sin(κ * x) / κ + kt / κ^2 * (1 - cos(κ * x)) + kk / κ^2 * x
    C("1.5.1", "eq:rmr0", "r(x)-r(0) = Rzz ε(x) e1 (closed form)", rq(Lr), R00 * eps_closed(Lr) * [1.0, 0.0, 0.0])
    Lsq = dot(rq(Lr), rq(Lr))
    C("1.5.1", "eq:rmr02 as written", "|r(L)-r(0)|² = [k1² κL + 4(κ²-k1²)sin²(κL/2)]/κ⁴", [Lsq],
      [(kp[1]^2 * (κ * Lr) + 4 * (κ^2 - kp[1]^2) * sin(κ * Lr / 2)^2) / κ^4]; kind=:typo)
    C("1.5.1", "eq:rmr02 corrected", "|r(L)-r(0)|² = [k1² (κL)² + 4(κ²-k1²)sin²(κL/2)]/κ⁴", [Lsq],
      [(kp[1]^2 * (κ * Lr)^2 + 4 * (κ^2 - kp[1]^2) * sin(κ * Lr / 2)^2) / κ^4])
    xoz, Δl = 0.7, 0.3
    chord_elem = norm(rq(xoz + Δl) - rq(xoz))              # chord between the element's own nodes
    chord_text = norm(rq(xoz + Δl))                        # what eq:rmr02 uses: r(xoz+Δℓ) - r(0), i.e. rz = r(0)
    C("1.5.1", "eq:rmr02 element form", "|r(xoz+Δℓ)-r(xoz)|² = [k1²(κΔℓ)² + 4(κ²-k1²)sin²(κΔℓ/2)]/κ⁴ (no xoz)", [chord_elem^2],
      [(kp[1]^2 * (κ * Δl)^2 + 4 * (κ^2 - kp[1]^2) * sin(κ * Δl / 2)^2) / κ^4])
    C("1.5.1", "eq:Deltal_k10 / rz as written", "NOTE: rz=r(0) (root) instead of r(xoz): chord from the root ≠ element chord for e>1", [chord_text], [chord_elem]; kind=:note)
    C("1.5.1", "eq:vars_e2 (note)", "NOTE: Rz is NOT constant over the element (midpoint value is an approximation)", [norm(Rz(xoz + Δl) - Rz(xoz))], [0.0]; kind=:note)

    # ------------------------------------------------------------- 1.5.2 Newmark identities
    β, α, Δt = 0.3, 0.5, 0.07
    y0, yd0, ydd0, yddn1 = v3(), v3(), v3(), v3()
    y1 = y0 + yd0 * Δt + ((1 - 2β) * ydd0 + 2β * yddn1) * Δt^2 / 2
    yd1 = yd0 + ((1 - α) * ydd0 + α * yddn1) * Δt
    C("1.5.2", "eq:Newmark_ydot2", "ẏ(t+dt) from the two Newmark relations (general α, β)", yd1,
      (1 - α / β) * yd0 + (1 - α / (2β)) * Δt * ydd0 + α / (β * Δt) * (y1 - y0))
    α2 = 2β
    yd1b = yd0 + ((1 - α2) * ydd0 + α2 * yddn1) * Δt
    C("1.5.2", "eq:Newmark_ydot3", "α = 2β: ẏ(t+dt) = (2/dt) y(t+dt) - [ẏ(t) + (2/dt) y(t)]", yd1b,
      (2 / Δt) * y1 - (yd0 + (2 / Δt) * y0))

    # ------------------------------------------------------------- 1.5.2 / 1.5.3 linear-algebra identities
    Jw, Rw = randn(rng, 3, 5), v3()
    C("1.5.2", "eq:sys_trim_sol", "Δ = -J⁺R for the wide trim Jacobian (J⁺ = Jᵀ(JJᵀ)⁻¹, Moore-Penrose)", -Jw' * inv(Jw * Jw') * Rw, -pinv(Jw) * Rw)
    nn = normalize(v3())
    Cfun(ph) = ph - dot(ph, nn) * nn                       # eq:hinge_constraint with ψ = p_h·n
    Jc = ForwardDiff.jacobian(Cfun, v3())
    C("1.5.3", "eq:hinge_constraint", "rank ∂C/∂p_h = 2 (one free component)", [Float64(rank(Jc))], [2.0])
    C("1.5.3", "eq:hinge_constraint", "the free direction of ∂C/∂p_h is the hinge axis n", Jc * nn, zeros(3))
    ψh = 1.1
    C("1.5.3", "eq:hinge_constraint (as written)", "WM: p_h = ψ n gives the rotation by ψ about n", R_WM(ψh * nn), R_exp(ψh * nn); kind=:typo)
    C("1.5.3", "eq:hinge_constraint (corrected)", "WM: p_h = g(ψ) n, g(ψ)=4tan(ψ/4), gives the rotation by ψ about n", R_WM(4 * tan(ψh / 4) * nn), R_exp(ψh * nn))
    return nothing
end
section15_checks()

# =============================================================================
# Appendix C (rotation parametrization derivatives, Wiener-Milenkovic): eq:pscaled_and_pext ... eq:dT_dp1
# =============================================================================
function appendixC_rotation_checks()
    C(sec, label, desc, a, b; kw...) = check("AppC", sec, label, desc, a, b; kw...)
    rng = MersenneTwister(33)
    Th_m(p) = begin                                   # eq:T_explicit (p0 = 2 - |p|²/8, eq:p0)
        p0 = 2 - dot(p, p) / 8
        [p0^2+p[1]^2-p[2]^2-p[3]^2  2*(p[1]*p[2]-p0*p[3])  2*(p[1]*p[3]+p0*p[2]);
         2*(p[1]*p[2]+p0*p[3])  p0^2-p[1]^2+p[2]^2-p[3]^2  2*(p[2]*p[3]-p0*p[1]);
         2*(p[1]*p[3]-p0*p[2])  2*(p[2]*p[3]+p0*p[1])  p0^2-p[1]^2-p[2]^2+p[3]^2]
    end
    for trial in 1:3
        pext = 9 * randn(rng, 3) / 3                  # |p_ext| > 4: rescaling is active
        qv = round(norm(pext) / 8)                    # eq:q (the rounding is locally constant)
        lam(pe) = 1 - 8 * qv / norm(pe)               # eq:scaling_fac
        pscl(pe) = lam(pe) * pe                       # eq:pscaled_and_pext
        λv = lam(pext)
        C("C.rot", "eq:pscaled_and_pext (t$trial)", "|p̊| < 4 after rescaling", [norm(pscl(pext)) < 4 ? 0.0 : 1.0], [0.0])
        C("C.rot", "eq:dlambda_dpext (t$trial)", "∂λ/∂p_ext = (1-λ)/|p_ext|² p_ext", ForwardDiff.gradient(lam, pext), (1 - λv) / dot(pext, pext) * pext)
        C("C.rot", "eq:dpscaled_dpext (t$trial)", "∂p̊/∂p_ext = (∂λ/∂p_ext) p_extᵀ + λ I", ForwardDiff.jacobian(pscl, pext),
          ForwardDiff.gradient(lam, pext) * pext' + λv * I)
        ps = pscl(pext)
        ups2(p) = (1 / (2 + dot(p, p) / 8))^2         # υ² with υ = (4-p0)⁻¹, p0 = 2 - |p|²/8
        C("C.rot", "eq:dups_dpscaled (t$trial)", "∂υ²/∂p̊ = -υ³ p̊/2", ForwardDiff.gradient(ups2, ps),
          -(1 / (2 + dot(ps, ps) / 8))^3 * ps / 2)
        C("C.rot", "eq:R_WM / eq:T_explicit (t$trial)", "R = υ² Θ with Θ from eq:T_explicit", (1 / (2 + dot(ps, ps) / 8))^2 * Th_m(ps), R_WM(ps))
        p1, p2, p3 = ps
        p0 = 2 - dot(ps, ps) / 8
        dT1_text = 2 * [p1*(1 - p0/4)  p2+p1*p3/4  p3-p1*p2/4;
                        p2-p1*p3/4  -p1*(1+p0/4)  -(p0-p1^2/4);
                        p3+p1*p2/4  p0-p1^2/4  -p1*(1+p0/4)]
        C("C.rot", "eq:dT_dp1 (t$trial)", "∂Θ/∂p̊₁ matrix of the text", ForwardDiff.derivative(s -> Th_m([s, p2, p3]), p1), dT1_text)
    end
    return nothing
end
appendixC_rotation_checks()

# ----------------------------------------------------------------------------- report
println()
let lastp = ""
    for r in RESULTS
        mark = r.kind == :pass ? (r.ok ? "ok" : "FALHOU") :
               r.kind == :typo ? (r.ok ? "ERRO NO TEXTO CONFIRMADO" : "!! esperava-se discrepancia e nao houve") :
                                 (r.ok ? "NOTA (grandezas diferentes, ver relatorio)" : "!! esperava-se diferenca e nao houve")
        if r.param != lastp
            println("\n------------------------- parametrization: ", r.param, " -------------------------")
            lastp = r.param
        end
        @printf("%-5s %-30s %-90s %9.2e  %s\n", r.sec, r.label, first(r.desc, 90), r.err, mark)
    end
end
nok = count(r -> r.kind == :pass && r.ok, RESULTS)
nbad = count(r -> r.kind == :pass && !r.ok, RESULTS)
ntypo = count(r -> r.kind == :typo && r.ok, RESULTS)
nnote = count(r -> r.kind == :note && r.ok, RESULTS)
nmiss = count(r -> r.kind != :pass && !r.ok, RESULTS)
println("\nSUMMARY: ", nok, " checks passed, ", nbad, " FAILED, ", ntypo, " errors-in-text confirmed (formula as written disagrees), ", nnote, " notes (quantities that differ), ", nmiss, " expected discrepancies that did not occur")
