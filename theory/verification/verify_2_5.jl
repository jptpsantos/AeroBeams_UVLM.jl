# =============================================================================
# Numerical / conceptual verification of section 2.5 ("Kinematical description in terms of beam variables") of theory/main.tex
# =============================================================================
# Companion of verify_chapter1.jl and verify_C.jl (same conventions: check(...) kinds pass / typo / note, random data, forward-mode AD).
#
# What is checked (equation labels of main.tex)
#   G   geometry of the aerodynamic basis: Rw^+ from (Λ, φ) [order of the rotations, sign of Λ and φ, explicit matrix], Rw = Rz Rw^+ Rz^T
#   B   chain of bases: eq:vstar_to_vsquare, eq:Us, eq:Ugs (every line of the chain, incl. the macros Ro, Rwo, Rzo, RAo) and the
#       version of the text BEFORE the red \fix{} (that one must fail)
#   W   rigid-body kinematics of a chord point P: eq:Usp, eq:Uspdot, eq:Usp3, eq:Uspdot3 and the meaning of "U-dot" (rate of the components
#       in the moving basis W, not components of the inertial acceleration)
#   E   equivalences with thin-airfoil variables: eq:U_equiv, eq:thetadot_equv, eq:a_equiv, eq:d_equiv, eq:w_pp_equiv/eq:v_pp_equiv,
#       eq:vg_equiv, eq:wg_equiv (+ eqs. w_34pp_QS, w_12pp_QS of section 2.1, which are the s_P = 3/4, 1/2 cases of eq:Usp3)
#   S   sweep: "chord becomes c secΛ", simple sweep theory (independence principle) vs. the loads of eqs. f_a / m_a / cfC_cfNC
#   C   the SAME formulas evaluated with the actual code of AeroBeams (geometry of the element, aero_steady_kinematics!,
#       aero_unsteady_kinematics!, local_gust_velocity!), so that text and code are compared with each other
#   M   mutation tests: wrong variants of the key formulas must be detected by the checks above
#
# Result categories (last column of the report)
#   ok                         agreement
#   ERRO NO TEXTO CONFIRMADO   the formula evaluated exactly AS PRINTED in main.tex (before the red fix) disagrees
#   NOTA                       statement that is imprecise / not general / inconsistent with the implementation (see report)
#
# Run:   julia --project=<AeroBeams root> theory/verification/verify_2_5.jl      (~20 s)
# =============================================================================

using LinearAlgebra, ForwardDiff, Random, Printf
using AeroBeams, StaticArrays

const TOL = 1e-9
const TEXFILE = joinpath(@__DIR__, "..", "main.tex")

# ----------------------------------------------------------------------------- notation
tilde(v) = [0 -v[3] v[2]; v[3] 0 -v[1]; -v[2] v[1] 0]
axial(S) = [S[3, 2], S[1, 3], S[2, 1]]
const e1 = [1.0, 0.0, 0.0]
const e2 = [0.0, 1.0, 0.0]
const e3 = [0.0, 0.0, 1.0]
const I3 = Matrix(1.0I, 3, 3)
Rxm(a) = [1 0 0; 0 cos(a) -sin(a); 0 sin(a) cos(a)]
Rym(a) = [cos(a) 0 sin(a); 0 1 0; -sin(a) 0 cos(a)]
Rzm(a) = [cos(a) -sin(a) 0; sin(a) cos(a) 0; 0 0 1]
function R_exp(p)                      # rotation vector -> rotation matrix (AD friendly away from p = 0)
    φ = norm(p)
    N = tilde(p / φ)
    I3 + sin(φ) * N + (1 - cos(φ)) * N * N
end
# Euler 3-2-1 copy of AeroBeams.rotation_tensor_E321 (src/Utilities.jl), the one used for Rw^+ = rotation_tensor_E321([Λ; 0; φ])
function R_E321(p)
    yaw, roll, pitch = p
    sy, cy = sincos(yaw); sr, cr = sincos(roll); sp, cp = sincos(pitch)
    [cr*cy cy*sp*sr-cp*sy sp*sy+cp*cy*sr;
     cr*sy cy*cp+sp*sr*sy cp*sr*sy-cy*sp;
     -sr cr*sp cr*cp]
end

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
function checkb(sec, label, desc, cond::Bool; kind=:pass)
    ok = kind == :pass ? cond : !cond
    push!(RESULTS, Res(sec, label, desc, cond ? 0.0 : 1.0, kind, ok))
    return ok
end
# mutation test: the WRONG variant must disagree with the truth (shows that the check above has power)
function mutant(sec, label, desc, truth, wrong)
    err = maximum(abs.(truth .- wrong)) / max(1.0, maximum(abs.(truth)), maximum(abs.(wrong)))
    push!(RESULTS, Res(sec, label, "[mutation] " * desc, err, :pass, err > 1e-6))
end

# =============================================================================
# G  Geometry of the aerodynamic basis
# =============================================================================
function section_geometry()
    s = "2.5-G"
    rng = MersenneTwister(2501)
    Λ, φ = deg2rad(27.0), deg2rad(-6.0)
    cΛ, sΛ, cφ, sφ = cos(Λ), sin(Λ), cos(φ), sin(φ)

    # Rw^+ in the code: rotation_tensor_E321([Λ; 0; φ]); explicit matrix = Rz(Λ) Rx(φ)
    Rwp = R_E321([Λ, 0, φ])
    explicit = [cΛ -sΛ*cφ sΛ*sφ; sΛ cΛ*cφ -cΛ*sφ; 0 sφ cφ]
    check(s, "Rw^+ (code)", "rotation_tensor_E321([Λ;0;φ]) = Rz(Λ)·Rx(φ)  (twist about b1 FIRST, then sweep about b3; fixed axes of b^+)", Rwp, Rzm(Λ) * Rxm(φ))
    check(s, "Rw^+ explicit", "= [cΛ -sΛcφ sΛsφ; sΛ cΛcφ -cΛsφ; 0 sφ cφ]", Rwp, explicit)
    # other possible reading of "rotate about b3 and b1, respectively": sweep first, then twist about the FIXED b1: Rx(φ) Rz(Λ)
    check(s, "order of rotations", "text read literally (Λ about b3, then φ about the fixed b1) = Rx(φ)Rz(Λ) ≠ code; code = twist first, sweep after (as in the side-view figure 'before being swept')", Rwp, Rxm(φ) * Rzm(Λ); kind=:note)
    # equivalently: sweep first and twist about the swept span axis w1
    check(s, "order of rotations", "code = Rz(Λ) then twist about the NEW axis w1 (intrinsic z-x' sequence)", Rwp, Rzm(Λ) * Rxm(φ))

    # columns of Rw^+ are the w_i in b^+ components; b1 in w components = Rw^+T e1
    b1w = Rwp' * e1
    φ0 = 0.0
    b1w0 = (Rzm(Λ) * Rxm(φ0))' * e1
    checkb(s, "sign of Λ", "Λ > 0 is AFT sweep: b1 (beam axis) has a component −w2 (towards the TE, since b2 points to the LE)", b1w0[2] < 0 && b1w0[1] > 0)
    # positive twist φ is nose-up: flow along -b2 gives α = +φ  (U_w = -Rw^+T V with V = U e2 for Λ = 0)
    U0 = 10.0
    Uw = -(Rzm(0.0) * Rxm(φ))' * (U0 * e2)
    check(s, "sign of φ", "φ > 0 (Λ=0): α = atan(-U3/U2) = φ  (nose-up = wash-in)", atan(-Uw[3] / Uw[2]), φ)

    # ---- Rw = Rz Rw^+ Rz^T (second-order tensor change of basis b^+ -> A)
    Rz0 = R_exp(1.1 * randn(rng, 3))                       # components in A of the basis b^+ (columns b_i^A)
    bA = [Rz0[:, i] for i in 1:3]
    wb = [Rwp[:, i] for i in 1:3]                          # w_i resolved in b^+ (what Rw^+ means)
    wA = [Rz0 * wb[i] for i in 1:3]                        # w_i resolved in A
    Bm = hcat(bA...); Wm = hcat(wA...)
    Rw_tensor_A = Wm * inv(Bm)                             # tensor that maps b_i -> w_i, resolved in A
    check(s, "Rw = Rz Rw^+ Rz^T", "tensor with w_i = Rw b_i, components in A, equals Rz Rw^+ Rz^T (as printed)", Rw_tensor_A, Rz0 * Rwp * Rz0')
    mutant(s, "Rw = Rz Rw^+ Rz^T", "wrong order Rz^T Rw^+ Rz is not the tensor in A", Rw_tensor_A, Rz0' * Rwp * Rz0)
    check(s, "RwR0 (code)", "code: RwR0 = R0·Rw^+  equals  Rw·Rz  (matrix of the w_i in A components)", Rz0 * Rwp, (Rz0 * Rwp * Rz0') * Rz0)

    # ---- deformed bases: B_i = R b_i (text), W_i = R w_i; Rw^+ is also the components of the tensor B* -> W resolved in B*
    pR = 0.9 * randn(rng, 3)
    R = R_exp(pR)
    Rp = Rz0' * R * Rz0                                    # R^+ (R resolved in b^+)
    Bstar = [R * bA[i] for i in 1:3]; Wsq = [R * wA[i] for i in 1:3]
    Tb = hcat(Wsq...) * inv(hcat(Bstar...))                # tensor B_i -> W_i resolved in A
    BstarM = hcat(Bstar...)
    check(s, "Rw^+ in B*", "the tensor B*->W has the SAME components Rw^+ in basis B* as Rw in b^+ (justifies eq:vstar_to_vsquare in one line)", BstarM' * Tb * BstarM, Rwp)
    check(s, "R B_i = R^+ ...", "R Rz = Rz R^+ (R^+ = Rz^T R Rz)", R * Rz0, Rz0 * Rp)
    return nothing
end

# =============================================================================
# B  chain of bases: eq:vstar_to_vsquare, eq:Us, eq:Ugs
# =============================================================================
function section_bases()
    s = "2.5-B"
    rng = MersenneTwister(2502)
    Λ, φ = deg2rad(-17.0), deg2rad(5.0)
    for trial in 1:3
        Rwp = R_E321([Λ, 0, φ])
        Rz0 = R_exp(1.2 * randn(rng, 3))
        R = R_exp(0.8 * randn(rng, 3))
        RA = R_exp(1.3 * randn(rng, 3))
        Rw = Rz0 * Rwp * Rz0'
        Rp = Rz0' * R * Rz0
        vA = randn(rng, 3)                                 # a vector resolved in A
        # truth: components along the basis vectors W_i = R w_i, w_i = Rw b_i, b_i = Rz e_i (all resolved in A)
        W = [R * (Rz0 * Rwp * e) for e in (e1, e2, e3)]
        vsq = [dot(W[i], vA) for i in 1:3]
        vplus = Rz0' * vA                                  # components in b^+
        vstar = [dot((R * Rz0 * e), vA) for e in (e1, e2, e3)]    # components in B*
        t = " (t$trial)"
        check(s, "B* comps" * t, "v^* = (R Rz)^T v^A = (R^+)^T v^+", vstar, Rp' * vplus)
        check(s, "eq:vstar_to_vsquare" * t, "line 1 (red fix): v^□ = (R^+ Rw^+)^T v^+", vsq, (Rp * Rwp)' * vplus)
        check(s, "eq:vstar_to_vsquare" * t, "line 2 (red fix): v^□ = (Rw^+)^T (R^+)^T v^+", vsq, Rwp' * Rp' * vplus)
        check(s, "eq:vstar_to_vsquare" * t, "line 3 (red fix): v^□ = (Rw^+)^T v^*", vsq, Rwp' * vstar)
        check(s, "eq:vstar_to_vsquare" * t, "version BEFORE the fix: v^□ = (R Rw)^T v^+ = Rw^T v^*  (mixes A and b^+ components)", vsq, (R * Rw)' * vplus; kind=:typo)
        check(s, "eq:vstar_to_vsquare" * t, "version BEFORE the fix, last line: v^□ = Rw^T v^*", vsq, Rw' * vstar; kind=:typo)
        check(s, "(R Rw Rz)^T" * t, "v^□ = (R Rw Rz)^T v^A  (A components; the chain used in eq:Ugs)", vsq, (R * Rw * Rz0)' * vA)
        # ---- eq:Us:  U^□ = -(Rw^+)^T V^*  with V^* the inertial sectional velocity in B*
        VA = randn(rng, 3)                                 # inertial velocity of point b, components in A
        Vstar = (R * Rz0)' * VA
        Uphys = [dot(W[i], -VA) for i in 1:3]              # relative wind = -velocity, components along W_i
        check(s, "eq:Us" * t, "U^□ = -(Rw^+)^T V^*  (fixed version)", Uphys, -Rwp' * Vstar)
        check(s, "eq:Us" * t, "version BEFORE the fix: U^□ = -Rw^T V^*", Uphys, -Rw' * Vstar; kind=:typo)
        # ---- eq:Ugs: all lines, with the meaning of the macros  (\Ro = R^∘ etc. are the A->I components of each tensor)
        Ugo = randn(rng, 3)
        Ro = RA * R * RA'; Rwo = RA * Rw * RA'; Rzo = RA * Rz0 * RA'; RAo = RA
        truth = vsq_from_inertial(Ugo, RA, Rz0, R, Rwp)    # step by step: I -> A -> b^+ -> B* -> W
        check(s, "eq:Ugs" * t, "line 1: (Ro Rwo Rzo RAo)^T Ugo   [Ro = RA R RA^T, ...]", truth, (Ro * Rwo * Rzo * RAo)' * Ugo)
        check(s, "eq:Ugs" * t, "line 2: [(RAo R RAo^T)(RAo Rw RAo^T)(RAo Rz RAo^T) RAo]^T Ugo", truth, ((RAo * R * RAo') * (RAo * Rw * RAo') * (RAo * Rz0 * RAo') * RAo)' * Ugo)
        check(s, "eq:Ugs" * t, "line 3: (RA R Rw Rz)^T Ugo   (uses RA = RAo)", truth, (RA * R * Rw * Rz0)' * Ugo)
        check(s, "eq:Ugs" * t, "line 4: (R Rw Rz)^T RA^T Ugo   (macro \\RRwRzT)", truth, (R * Rw * Rz0)' * RA' * Ugo)
        check(s, "RA = RAo" * t, "components in A of the tensor I->A: RA^T RAo RA = RAo  (so RA = RAo, used in line 3; stated only in C.3)", RA' * RAo * RA, RAo)
        check(s, "eq:Ugs vs code" * t, "code (local_gust_velocity!): (R·RwR0)' · R_AT · Ugo,  RwR0 = R0 Rw^+,  R_AT = RA^T", truth, (R * (Rz0 * Rwp))' * (RA') * Ugo)
        mutant(s, "eq:Ugs" * t, "using Rw^+ in place of Rw in (R Rw Rz)^T RA^T Ugo", truth, (R * Rwp * Rz0)' * RA' * Ugo)
        mutant(s, "eq:Ugs" * t, "omitting RA^T", truth, (R * Rw * Rz0)' * Ugo)
    end
    return nothing
end
# I -> A -> b^+ -> B* -> W  (vector components, step by step)
function vsq_from_inertial(Ugo, RA, Rz0, R, Rwp)
    vA = RA' * Ugo
    vplus = Rz0' * vA
    Rp = Rz0' * R * Rz0
    vstar = Rp' * vplus
    return Rwp' * vstar
end

# =============================================================================
# W  wind at a chord point P of a moving rigid section (eq:Usp, eq:Uspdot, eq:Usp3, eq:Uspdot3)
# =============================================================================
function section_wind()
    s = "2.5-W"
    rng = MersenneTwister(2503)
    c, sb, sP = 0.37, 0.41, 0.73                           # chord, s (attachment point), s_P (generalized point), both from the LE (L)
    M0 = R_exp(1.0 * randn(rng, 3))
    a1, a2, a3 = randn(rng, 3), 0.5 * randn(rng, 3), 0.6 * randn(rng, 3)
    r0, r1, r2, r3 = randn(rng, 3), randn(rng, 3), 0.5 * randn(rng, 3), 0.4 * randn(rng, 3)
    M(t) = M0 * R_exp(a1 * t + a2 * t^2 + a3 * sin(1.3t + 0.4))       # columns = W_i (resolved in I)
    rb(t) = r0 + r1 * t + r2 * t^2 + r3 * sin(0.9t)                    # inertial position of the attachment point b
    xL(t) = rb(t) + M(t) * (sb * c * e2)                               # LE: +w2 side, at distance s c from b
    xT(t) = xL(t) - M(t) * (c * e2)                                    # TE: chord c towards -w2
    xP(t) = xL(t) - M(t) * (sP * c * e2)                               # P at s_P c from L towards T
    d(f) = t -> ForwardDiff.derivative(f, t)
    Vsq(t) = M(t)' * d(rb)(t)                                          # V^□ (inertial velocity of b, components in W)
    Osq(t) = axial(M(t)' * d(M)(t))                                    # Ω^□  (M^T Mdot = tilde(Ω^□))
    vP(t) = M(t)' * d(xP)(t)
    t0 = 0.31
    # geometry of the points
    ξP = c * [0.0, sb - sP, 0.0]
    check(s, "ξ_P", "|P-L| = s_P c,  |T-L| = c,  P on the chordline, ξ_P^□ = M^T (x_P - x_b) = c{0, s - s_P, 0}", [norm(xP(t0) - xL(t0)), norm(xT(t0) - xL(t0)), norm(cross(xT(t0) - xL(t0), xP(t0) - xL(t0)))], [sP * c, c, 0.0])
    check(s, "ξ_P", "ξ_P^□ = c{0, s - s_P, 0}  (eq:Usp3 text, position of P with respect to b, in W)", M(t0)' * (xP(t0) - rb(t0)), ξP)
    U = -Vsq(t0); Ω = Osq(t0)
    UP = -vP(t0)                                           # relative wind at P, components in W (the truth)
    check(s, "eq:Usp", "U_P = U - Ω~ ξ_P  (ξ_P = c{0,s-s_P,0}, all in W)", UP, U - tilde(Ω) * ξP)
    check(s, "eq:Usp3", "U_P3 = U_3 - Ω_1 c (s - s_P)", UP[3], U[3] - Ω[1] * c * (sb - sP))
    mutant(s, "eq:Usp", "wrong sign U + Ω~ ξ_P", UP, U + tilde(Ω) * ξP)
    mutant(s, "eq:Usp", "ξ_P with (s_P - s)", UP, U - tilde(Ω) * (-ξP))
    # rates: dU_P/dt (components) = Udot - Ωdot~ ξ_P   (ξ_P^□ is constant in W)
    Udot = ForwardDiff.derivative(t -> -Vsq(t), t0)
    Ωdot = ForwardDiff.derivative(Osq, t0)
    UPdot = ForwardDiff.derivative(t -> -vP(t), t0)
    check(s, "eq:Uspdot", "d/dt of the components of U_P = Udot - tilde(Ωdot) ξ_P", UPdot, Udot - tilde(Ωdot) * ξP)
    check(s, "eq:Uspdot3", "Udot_P3 = Udot_3 - Ωdot_1 c (s - s_P)", UPdot[3], Udot[3] - Ωdot[1] * c * (sb - sP))
    # meaning of the dot: rate of the COMPONENTS in the moving basis W  vs  components of the inertial acceleration
    aP_sq = M(t0)' * ForwardDiff.derivative(t -> ForwardDiff.derivative(xP, t), t0)       # components of the inertial acceleration of P
    check(s, "meaning of Udot", "components of -(inertial acceleration of P) = Udot_P + Ω~ U_P  (transport term)", -aP_sq, UPdot + tilde(Ω) * UP)
    check(s, "meaning of Udot", "Udot_P (text/code) is NOT -(acceleration of P) in W: differs by Ω~ U_P", UPdot, -aP_sq; kind=:note)
    # thin-airfoil argument: pitching airfoil in a steady uniform wind (inertially constant, nothing accelerates): needs d/dt of the COMPONENT
    U∞, θ0, θ1 = 12.0, 0.1, 0.8
    Mp(t) = Rxm(θ0 + θ1 * t)                                           # pitching about w1 = ι1 (nose-up positive)
    rp(t) = [0.0, U∞ * t, 0.0]                                         # airfoil moves towards its LE at speed U∞ (wind along -ι2)
    Vp(t) = Mp(t)' * ForwardDiff.derivative(rp, t)
    Op(t) = axial(Mp(t)' * ForwardDiff.derivative(Mp, t))
    Up(t) = -Vp(t)
    tt = 0.2; θ = θ0 + θ1 * tt
    Upd = ForwardDiff.derivative(Up, tt)
    check(s, "thin airfoil", "w^QS = U3 = U sinθ  (eq:w_pp_QS with hdot = 0) and  θdot = Ω_1^□", [Up(tt)[3], Op(tt)[1]], [U∞ * sin(θ), θ1])
    check(s, "thin airfoil", "d/dt(w^QS) = U∞ cosθ θdot  (what the apparent-mass loads need) = Udot_3 (rate of the component)", Upd[3], U∞ * cos(θ) * θ1)
    check(s, "thin airfoil", "inertial acceleration of the (steady, uniform) wind is zero: Udot + Ω~ U = 0, yet Udot_3 ≠ 0", Upd + tilde(Op(tt)) * Up(tt), zeros(3))
    checkb(s, "thin airfoil", "so Udot must be the rate of the components (Udot_3 = $(round(Upd[3], digits=3)) ≠ 0), not the components of the acceleration (= 0)", abs(Upd[3]) > 1e-3)
    return nothing
end

# =============================================================================
# E  equivalences with thin-airfoil variables
# =============================================================================
function section_equivalences()
    s = "2.5-E"
    rng = MersenneTwister(2504)
    c = 0.23; b = c / 2
    sb = 0.37; sf = 0.8
    a = 2 * sb - 1; dflap = 2 * sf - 1
    # geometry: a = semichord-normalized position of b aft of midchord; coordinates along w2 (LE = +w2) from b
    xLE = sb * c; xTE = xLE - c; xmid = (xLE + xTE) / 2; xh = xLE - sf * c            # hinge at s_f from LE
    check(s, "eq:a_equiv", "a = (aft distance of b from midchord)/b = 2s - 1", (xmid - 0.0) / b, a)
    check(s, "eq:d_equiv", "d = (aft distance of the hinge from midchord)/b = 2 s_f - 1", (xmid - xh) / b, dflap)
    # w_{3/4}, w_{1/2} (section 2.1) as the s_P = 3/4, 1/2 cases of eq:Usp3, with U3 ≡ U sinθ + hdot cosθ, thetadot ≡ Ω1
    U, θ, hdot, θdot = 14.0, 0.21, -1.3, 0.9
    U3 = U * sin(θ) + hdot * cos(θ); Ω1 = θdot
    check(s, "eq:w_34pp_QS", "U_P3 (s_P=3/4) = U3 - Ω1 c (s - 3/4)  =  U sinθ + hdot cosθ + b(1/2 - a) θdot", U3 - Ω1 * c * (sb - 3 / 4), U * sin(θ) + hdot * cos(θ) + b * (1 / 2 - a) * θdot)
    check(s, "eq:w_12pp_QS", "U_P3 (s_P=1/2) = U3 - Ω1 c (s - 1/2)  =  U sinθ + hdot cosθ - a b θdot", U3 - Ω1 * c * (sb - 1 / 2), U * sin(θ) + hdot * cos(θ) - a * b * θdot)
    # eq:w_pp_equiv / eq:v_pp_equiv and U ≡ sqrt(U2^2+U3^2)
    v = U * cos(θ) - hdot * sin(θ); w = U * sin(θ) + hdot * cos(θ)                  # eqs. v_pp_QS, w_pp_QS
    check(s, "eq:alpha_pp_QS", "α^QS = atan(w/v) = atan(-U3/U2) with v ≡ -U2, w ≡ U3  (code: α = atan(-Uₙ/Uₜ))", atan(w / v), atan(-w / (-v)))
    check(s, "eq:U_equiv", "sqrt(v_pp^2 + w_pp^2) = sqrt(U^2 + hdot^2): equals U only if hdot = 0 (the three equivalences U≡U_ip, v≡-U2, w≡U3 force hdot≡0, θ≡α^QS)", hypot(v, w), U; kind=:note)
    check(s, "eq:U_equiv", "identity hypot(v,w)^2 = U^2 + hdot^2", hypot(v, w)^2, U^2 + hdot^2)
    # code: Uᵢ = sqrt(Uₜ² + Uₙ²), α = atan(-Uₙ/Uₜ), pitch rate = Ω1^□;  Uₙ = U_ip sin α  for the code's α (Uₜ<0)
    Ut, Un = -9.0, 1.7
    α = atan(-Un / Ut); Ui = hypot(Ut, Un)
    check(s, "α vs Uₙ", "with v = -U2 > 0:  U3 = U_ip sin α  and  -U2 = U_ip cos α  (so U sinθ + hdot cosθ with hdot=0, θ=α)", [Un, -Ut], [Ui * sin(α), Ui * cos(α)])
    # gust: relative air velocity = air - body = U + Ug  (components in W);  v_g ≡ -Ug2,  w_g ≡ Ug3
    Ug = randn(rng, 3); Vb = randn(rng, 3)
    Urel = Ug - Vb                                         # relative wind including the gust
    Uq = -Vb
    check(s, "eq:wg_equiv", "total normalwash = U3 + Ug3 = w_pp^QS + w_g^QS,  w_g ≡ Ug3", Urel[3], Uq[3] + Ug[3])
    check(s, "eq:vg_equiv", "total tangential velocity -(U2 + Ug2) = v_pp^QS + v_g^QS,  v_g ≡ -Ug2  (used in eq:alpha_E; code: Uₜ + UₜGust)", -Urel[2], -Uq[2] + (-Ug[2]))
    # label/content swap of eq:w_pp_equiv / eq:v_pp_equiv (cosmetic): read main.tex
    if isfile(TEXFILE)
        L = readlines(TEXFILE)
        iw = findfirst(l -> occursin("\\label{eq:w_pp_equiv}", l), L)
        iv = findfirst(l -> occursin("\\label{eq:v_pp_equiv}", l), L)
        swapped = !isnothing(iw) && !isnothing(iv) && occursin("v_{pp}^{QS}", L[iw-1]) && occursin("w_{pp}^{QS}", L[iv-1])
        checkb(s, "labels", "labels eq:w_pp_equiv / eq:v_pp_equiv were swapped in the original text (cosmetic); after the audit they name the w_pp^QS / v_pp^QS equations, respectively", !swapped)
    end
    return nothing
end

# =============================================================================
# S  sweep: chord c secΛ and simple sweep theory
# =============================================================================
function section_sweep()
    s = "2.5-S"
    c, Λ, Δℓ, ρ, U, cnα, cmα, α = 0.19, deg2rad(32.0), 0.11, 1.225, 21.0, 2π, -0.4, 0.013
    cΛ = cos(Λ)
    # top view (figure aero_section_top): beam axis b1 at angle Λ from w1 (towards -w2); strip = parallelogram with edges
    # along b1 (length Δℓ) and along w2 (streamwise chord cs).  Chord NORMAL to the beam axis = c  =>  cs = c secΛ
    b1 = [cΛ, -sin(Λ)]                                     # (w1,w2) components of b1  (Rw^+T e1, φ=0)
    b2 = [sin(Λ), cΛ]                                      # in-plane normal to b1 (towards the LE)
    cs = c / cΛ
    edge_chord = cs * [0.0, 1.0]                           # streamwise chord vector (w2 direction)
    check(s, "c secΛ", "streamwise chord cs such that the chord normal to b1 is c:  b2^T ... : (cs w2)·b2 = c  =>  cs = c secΛ", dot(edge_chord, b2), c)
    area = abs(Δℓ * b1[1] * edge_chord[2] - Δℓ * b1[2] * edge_chord[1])      # |Δℓ b1 × cs w2|
    check(s, "c secΛ", "strip area = Δℓ·c  (parallelogram: Δℓ along b1 times chord c normal to b1) = (c secΛ)(Δℓ cosΛ)", area, Δℓ * c)
    checkb(s, "c secΛ", "claim as printed: 'upon being swept (rotated by Λ about b3) the chord becomes c secΛ' -- a rigid rotation keeps the chord c; c secΛ is an assumption on the strip (parallelogram)", abs(norm(R_E321([Λ, 0, 0.0]) * (c * e2)) - c / cΛ) < 1e-14; kind=:note)
    # forces per element: simple sweep theory (independence principle): section normal to the axis, chord c, speed U cosΛ, AoA αn = α/cosΛ
    ncode = 0.5 * ρ * U^2 * (c / cΛ) * cΛ^2 * Δℓ * cnα * α                               # f_a,n with k = cos²Λ on the circulatory coefficient
    nsweep = 0.5 * ρ * (U * cΛ)^2 * c * Δℓ * cnα * (α / cΛ)
    check(s, "eq:f_a", "½ρU² (c secΛ) Δℓ · cos²Λ cnα α  =  ½ρ(U cosΛ)² c Δℓ cnα (α/cosΛ)  (independence principle, c = chord normal to the beam axis)", ncode, nsweep)
    mutant(s, "eq:f_a", "without the cos²Λ of eq:cfC_cfNC", nsweep, 0.5 * ρ * U^2 * (c / cΛ) * Δℓ * cnα * α)
    mcode = 0.5 * ρ * U^2 * (c / cΛ)^2 * Δℓ * (cΛ^2 * cmα * α) * (cΛ)                     # eq:m_a: (c secΛ)² cm [b1^T... = cosΛ]
    msweep = 0.5 * ρ * (U * cΛ)^2 * c^2 * Δℓ * cmα * (α / cΛ)
    check(s, "eq:m_a", "moment about the beam axis: ½ρU² (c secΛ)² Δℓ cos²Λ cmα α · cosΛ = ½ρ(U cosΛ)² c² Δℓ cmα (α/cosΛ)", mcode, msweep)
    # flow-aligned basis: beam yawed by -Λ w.r.t. A (as in test/examples/.../sweptWingBendingDivergence.jl, p0 = E321(-Λ,0,θ)): Rz Rw^+ = I
    Rz0 = R_E321([-Λ, 0, 0.0]); Rwp = R_E321([Λ, 0, 0.0])
    check(s, "w = A", "beam yawed -Λ in A and Rw^+ = Rz(Λ):  Rz Rw^+ = I  (w basis = flow-aligned basis A; b1 = (cosΛ,-sinΛ,0) swept aft)", Rz0 * Rwp, I3)
    Vstar = Rz0' * [0.0, U, 0.0]                           # A moves with +U a2 (wind along -a2), R = I
    Uw = -Rwp' * Vstar
    check(s, "U_ip", "then U^□ = (0,-U,0): spanwise component U_1 = 0 and U_ip = U (NOT U cosΛ): the 2-D section is the streamwise one", Uw, [0.0, -U, 0.0])
    # time scale remark: Θ = U_ip/b with b = c/2 (code) vs. independence principle U cosΛ/(c/2)  (modelling convention, see report)
    check(s, "Θ = U_ip/b", "reduced-time scale in the code/section 2.1 is U_ip/(c/2); the normal-section scale would be U_ip cosΛ/(c/2) = U_ip/(b secΛ)  (differs by cosΛ = $(round(cΛ, digits=3)))", U / (c / 2), U * cΛ / (c / 2); kind=:note)
    return nothing
end

# =============================================================================
# C  the same formulas against the code of AeroBeams
# =============================================================================
function section_code()
    s = "2.5-C"
    rng = MersenneTwister(2506)
    Λ, φ, sb, c = deg2rad(23.0), deg2rad(3.5), 0.37, 0.21
    surf = create_AeroSurface(ppSolver=Indicial(), airfoil=deepcopy(FlatPlate), c=c, Λ=Λ, φ=φ, normSparPos=sb, updateAirfoilParameters=false)
    p0 = [0.3, -0.5, 0.7]
    pA0 = [0.4, -0.2, 0.6]
    beam = create_Beam(name="b", length=1.0, nElements=2, S=[isotropic_stiffness_matrix(∞=1e9, EIy=10)], I=[inertia_matrix(ρA=1)], surface=surf, p0=p0)
    gust = create_SharpEdgedGust(initialTime=0.0, duration=10.0, convectiveVelocity=1.3, verticalVelocity=2.1, p=[0.2, 0.1, -0.3])
    clampBC = create_BC(name="clamp", beam=beam, node=1, conditions=Dict(:u1A => 0, :u2A => 0, :u3A => 0, :p1A => 0, :p2A => 0, :p3A => 0))
    model = create_Model(name="m", beams=[beam], BCs=[clampBC], p_A0=pA0, gust=gust)
    el = model.assemblyData.elements[1]
    g = el.aero.geometry
    cache = el.aero.stateCache
    Rwp = R_E321([Λ, 0, φ])
    R0 = AeroBeams.rotation_tensor_WM(p0)                  # straight beam: R0 (= Rz of the text) constant
    # the WM matrix of my notes vs the code's
    check(s, "geometry", "code Rw^+ (get_aero_geometric_properties) = Rz(Λ)Rx(φ)", Matrix(g.Rw⁺), Rwp)
    check(s, "geometry", "code Rw^+T = transpose", Matrix(g.Rw⁺T), Rwp')
    check(s, "geometry", "code RwR0 = R0 Rw^+  (= Rw Rz of the text, with R0 = Rz of the element)", Matrix(g.RwR0), R0 * Rwp)
    check(s, "geometry", "code RwR0 = (Rz Rw^+ Rz^T) Rz  with Rw = Rz Rw^+ Rz^T", Matrix(g.RwR0), (R0 * Rwp * R0') * R0)
    check(s, "geometry", "code a_h = 2s - 1, b = c/2, csecΛ = c/cosΛ, cosΛ² = cos²Λ", [g.aₕ, g.b, g.csecΛ, g.cosΛ²], [2sb - 1, c / 2, c / cos(Λ), cos(Λ)^2])
    # ---- steady kinematics: eq:Us, eq:U_equiv, eq:thetadot_equv, eq:Usp3 (s_P = 1/2, 3/4)
    V = randn(rng, 3); Ω = randn(rng, 3)
    AeroBeams.aero_steady_kinematics!(el, cache, V, Ω)
    k = cache.attachedVars.steadyKin
    Utext = -Rwp' * V                                      # eq:Us (V = V^*)
    Ωsq = Rwp' * Ω                                         # Ω^□ = (Rw^+)^T Ω^*  (same rule as eq:vstar_to_vsquare; not written in the text)
    check(s, "eq:Us", "code U^□ = -(Rw^+)^T V^*:  (Uₛ,Uₜ,Uₙ)", [k.Uₛ, k.Uₜ, k.Uₙ], Utext)
    check(s, "eq:U_equiv", "code Uᵢ = sqrt(Uₜ²+Uₙ²) = sqrt(U2²+U3²)", k.Uᵢ, hypot(Utext[2], Utext[3]))
    check(s, "eq:thetadot_equv", "code Ωₐ = dot(Rw⁺T[1,:], Ω) = Ω_1^□", k.Ωₐ, Ωsq[1])
    check(s, "eq:Usp3", "code UₙMid = U3 - Ω1 c (s - 1/2),  UₙTQC = U3 - Ω1 c (s - 3/4)", [k.UₙMid, k.UₙTQC], [Utext[3] - Ωsq[1] * c * (sb - 1 / 2), Utext[3] - Ωsq[1] * c * (sb - 3 / 4)])
    check(s, "eq:Usp", "full vector version: U_P = U - Ω~^□ ξ_P, normal component, s_P = 3/4", (Utext - tilde(Ωsq) * (c * [0.0, sb - 3 / 4, 0.0]))[3], k.UₙTQC)
    check(s, "α", "code α = atan(-Uₙ/Uₜ) = atan(w_pp/v_pp)  with v ≡ -U2, w ≡ U3", k.α, atan(Utext[3] / (-Utext[2])))
    # ---- rates: Udot = -Rw^+T Vdot (Vdot = d/dt of the components V^*), checked against AD of the scalar signals
    V0, V1, V2 = randn(rng, 3), randn(rng, 3), randn(rng, 3)
    O0, O1, O2 = randn(rng, 3), randn(rng, 3), randn(rng, 3)
    Vt(t) = V0 + V1 * t + V2 * sin(1.1t); Ot(t) = O0 + O1 * t + O2 * cos(0.7t)
    t0 = 0.4
    Vd = ForwardDiff.derivative(Vt, t0); Od = ForwardDiff.derivative(Ot, t0)
    AeroBeams.aero_steady_kinematics!(el, cache, Vt(t0), Ot(t0))
    AeroBeams.aero_unsteady_kinematics!(el, cache, Vd, Od)
    u = cache.attachedVars.unsteadyKin
    Ut(t) = -Rwp' * Vt(t); Ωt(t) = (Rwp' * Ot(t))[1]
    wtqc(t) = Ut(t)[3] - Ωt(t) * c * (sb - 3 / 4); wmid(t) = Ut(t)[3] - Ωt(t) * c * (sb - 1 / 2)
    Uip(t) = hypot(Ut(t)[2], Ut(t)[3]); αt(t) = atan(-Ut(t)[3] / Ut(t)[2])
    check(s, "eq:Uspdot3", "code UₙdotMid, UₙdotTQC = d/dt of the scalar signals w_{1/2}^QS, w_{3/4}^QS", [u.UₙdotMid, u.UₙdotTQC], [ForwardDiff.derivative(wmid, t0), ForwardDiff.derivative(wtqc, t0)])
    check(s, "eq:Uspdot", "code Udot = -Rw^+T Vdot,  Ωₐdot = Ω̇_1^□", [u.Udot..., u.Ωₐdot], [(-Rwp' * Vd)..., (Rwp' * Od)[1]])
    check(s, "Uᵢdot, αdot", "code Uᵢdot, αdot = d/dt sqrt(U2²+U3²), d/dt atan(-U3/U2)", [u.Uᵢdot, u.αdot], [ForwardDiff.derivative(Uip, t0), ForwardDiff.derivative(αt, t0)])
    # ---- gust: eq:Ugs vs local_gust_velocity! (p = WM parameters of R at the element midpoint)
    p = 0.7 * randn(rng, 3)
    tnow = 1.7
    AeroBeams.local_gust_velocity!(tnow, model, el, cache, p)
    k = cache.attachedVars.steadyKin
    R = AeroBeams.rotation_tensor_WM(p); RA = AeroBeams.rotation_tensor_WM(pA0)
    Ugo = Vector(gust.UGustInertial(tnow))
    Rw = R0 * Rwp * R0'
    Ugtext = (R * Rw * R0)' * RA' * Ugo                    # eq:Ugs (line 4)
    check(s, "eq:Ugs", "code (UₜGust, UₙGust) = (U_g2, U_g3) of (R Rw Rz)^T RA^T Ugo   [vg ≡ -Ug2, wg ≡ Ug3]", [k.UₜGust, k.UₙGust], Ugtext[2:3])
    check(s, "R_AT", "model.kinematicsBasisA.state.R_AT = RA^T (maps I components to A components; its docstring 'from A to I' refers to the tensor of the opposite direction)", Matrix(model.kinematicsBasisA.state.R_AT), RA')
    return nothing
end

# =============================================================================
println("="^100)
println("Verification of section 2.5 of theory/main.tex   -   tolerance = $TOL")
println("="^100)
section_geometry()
section_bases()
section_wind()
section_equivalences()
section_sweep()
section_code()

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
        @printf("%-7s %-24s %-150s %9.2e  %s\n", r.sec, first(r.label, 24), first(r.desc, 150), r.err, mark)
    end
end
nok = count(r -> r.kind == :pass && r.ok, RESULTS)
nbad = count(r -> r.kind == :pass && !r.ok, RESULTS)
ntypo = count(r -> r.kind == :typo && r.ok, RESULTS)
nnote = count(r -> r.kind == :note && r.ok, RESULTS)
nmiss = count(r -> r.kind != :pass && !r.ok, RESULTS)
println("\nSUMMARY: ", nok, " checks passed, ", nbad, " FAILED, ", ntypo, " errors-in-text confirmed (formula as printed disagrees), ", nnote, " notes, ", nmiss, " expected discrepancies that did not occur")
