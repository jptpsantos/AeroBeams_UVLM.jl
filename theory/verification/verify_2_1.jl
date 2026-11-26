# =============================================================================
# Conceptual + numerical verification of section 2.1 ("Attached-flow aerodynamics") of theory/main.tex
# =============================================================================
# Companion of verify_chapter1.jl / verify_C.jl (same conventions: check(...) kinds pass / typo / note).
#
# Sections covered
#   2.1      indicial (Duhamel) machinery
#   2.1.1    pitch-plunge-induced loads (quasi-steady normalwash, circulatory states, compressible and incompressible inertial loads)
#   2.1.2    finite-state induced-inflow (Peters)
#   2.1.3    trailing-edge flap loads (Theodorsen flap constants)
#   2.1.4    gust-induced loads
#   2.1.5    composition of loads (alpha^E, c_n^C, c_t^C, c_m^C, induced/parasite drag, c_l, c_d)
#   2.1.6    ad hoc flap formulation
#
# Method
#   (1) every equation is re-implemented exactly AS PRINTED in main.tex (line numbers refer to the file at the time of writing);
#   (2) independent references: explicit rigid-body kinematics (positions differentiated by AD), thin-airfoil integrals
#       (Birnbaum/Glauert) for the 3/4-chord normalwash and flap constants, Theodorsen's function C(k) (Bessel functions by
#       quadrature, validated by the Wronskian), piston theory integrals, exact solutions of the exponential-kernel
#       convolutions, Theodorsen's non-circulatory loads (Theodorsen 1935/Fung, h positive downward);
#   (3) comparison with src/Aerodynamics.jl, src/AeroSolver.jl, src/AeroProperties.jl (the script also greps the code lines
#       it replicates, so that it fails loudly if the code is edited).
#
# Result categories (last column of the report)
#   ok                         agreement (step/statement of the text is correct)
#   ERRO NO TEXTO CONFIRMADO   the formula evaluated exactly AS PRINTED in main.tex disagrees (typo); the corrected one passes
#   NOTA                       statement that is imprecise / not general / inconsistent with the implementation (see report)
#   For "checkv" rows the "err" column is the value of the metric, not a relative error.
#
# Run:   julia --project=<AeroBeams root> theory/verification/verify_2_1.jl      (~20 s)
# Last run (2026-11-14): 81 checks passed, 0 failed, 10 errors-in-text confirmed (T_nα counted once per indicial set), 9 notes.
# After the audit the corrections were applied to main.tex as red \fix{} (the script evaluates the formulas AS ORIGINALLY PRINTED, transcribed here).
# =============================================================================

using LinearAlgebra, ForwardDiff, QuadGK, Random, Printf

const TOL = 1e-9

# ----------------------------------------------------------------------------- real code (optional)
const HAVE_AB = try
    @eval using AeroBeams
    true
catch err
    @warn "AeroBeams could not be loaded: checks against the real implementation are skipped" exception = err
    false
end
const SRC = joinpath(@__DIR__, "..", "..", "src")
const CODE_AERO = read(joinpath(SRC, "Aerodynamics.jl"), String)
const CODE_SOLV = read(joinpath(SRC, "AeroSolver.jl"), String)
const CODE_PROP = read(joinpath(SRC, "AeroProperties.jl"), String)
squash(s) = replace(s, r"\s+" => " ")
code_has(code, pattern) = occursin(squash(pattern), squash(code))

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
function check(sec, label, desc, a, b; kind=:pass, tol=TOL)
    err = maximum(abs.(a .- b)) / max(1.0, maximum(abs.(a)), maximum(abs.(b)))
    ok = kind == :pass ? err < tol : err > max(1e-6, tol)
    push!(RESULTS, Res(sec, label, desc, err, kind, ok))
    return ok
end
# boolean version: `cond` = the claim holds. :pass -> must hold; :typo/:note -> claim as printed must be false
function checkb(sec, label, desc, cond::Bool; kind=:pass)
    ok = kind == :pass ? cond : !cond
    push!(RESULTS, Res(sec, label, desc, cond ? 0.0 : 1.0, kind, ok))
    return ok
end
const DIAG = String[]          # numbers quoted in the report
diag(args...) = push!(DIAG, string(args...))
const UNVERIFIED = String[]    # statements that cannot be checked numerically

# ----------------------------------------------------------------------------- tools
function rk4(f, x0, t0, t1, n)
    h = (t1 - t0) / n
    x = copy(x0)
    t = t0
    for _ in 1:n
        k1 = f(t, x); k2 = f(t + h / 2, x + h / 2 * k1); k3 = f(t + h / 2, x + h / 2 * k2); k4 = f(t + h, x + h * k3)
        x = x + h / 6 * (k1 + 2k2 + 2k3 + k4)
        t += h
    end
    return x
end
d(f, t) = ForwardDiff.derivative(f, t)

# Bessel functions by quadrature (no SpecialFunctions in the project); validated by the Wronskian J1 Y0 - J0 Y1 = 2/(pi x)
Jn(n, x) = quadgk(τ -> cos(n * τ - x * sin(τ)), 0, π; rtol=1e-13)[1] / π
Yn(n, x) = quadgk(τ -> sin(x * sin(τ) - n * τ), 0, π; rtol=1e-13)[1] / π -
           quadgk(t -> (exp(n * t) + (-1)^n * exp(-n * t)) * exp(-x * sinh(t)), 0, asinh(700 / x); rtol=1e-12)[1] / π
# Theodorsen's function C(k) = H1/(H1 + i H0), H = J - iY, input e^{ikt}, t in semichords (b/U)
Cth(k) = (H1 = Jn(1, k) - im * Yn(1, k); H0 = Jn(0, k) - im * Yn(0, k); H1 / (H1 + im * H0))

# Indicial parameters (src/AeroSolver.jl)
const WAGNER = ([0.165, 0.335], [0.0455, 0.3])
const BEDDOES = ([0.3, 0.7], [0.14, 0.53])
const JOSE = ([0.3493, 0.6507], [0.0984, 0.7759])
const AI = [1.5, -0.5, 1.0]; const bI = [0.25, 0.1, 5.0]

# Peters' finite-state arrays (copy of AeroSolver.jl inflow_arrays)
function inflow_arrays(N::Int)
    D = zeros(N, N); D[1, 2] = -1 / 2; D[end, end-1] = 1 / (2N)
    for n = 2:N-1
        D[n, n+1] = -1 / (2n); D[n, n-1] = 1 / (2n)
    end
    b = zeros(N)
    for n = 1:N-1
        b[n] = (-1)^(n - 1) * factorial(N + n - 1) / factorial(N - n - 1) / (factorial(n))^2
    end
    b[N] = (-1)^(N + 1)
    c = [2 / n for n in 1:N]
    dd = zeros(N); dd[1] = 1 / 2
    A = D + dd * b' + c * dd' + 1 / 2 * c * b'
    return A, inv(A), b, c
end

# =============================================================================
# 2.1  indicial (Duhamel) machinery   [l. 1888-1908, 1927-1970]
# =============================================================================
function sec_indicial()
    s = "2.1-indicial"
    C(l, dsc, a, b; kw...) = check(s, l, dsc, a, b; kw...)
    B(l, dsc, c; kw...) = checkb(s, l, dsc, c; kw...)
    A, bb = WAGNER
    Θ = 0.8                                           # constant beta^2 U/b
    Φ(t) = 1 - sum(A .* exp.(-bb .* Θ .* t))
    dΦ(t) = sum(A .* bb .* Θ .* exp.(-bb .* Θ .* t))
    α(t) = 0.1 + 0.3 * sin(1.3t) + 0.05t
    t = 3.7

    # eq:alpha_qs_conv
    C("eq:alpha_qs_conv", "α(t) = α(0)u(t) + ∫α̇ u(t-σ)dσ (trivial: Duhamel for the quasi-steady value)", α(t), α(0) + quadgk(σ -> d(α, σ), 0, t)[1])
    # eq:alpha_E_conv, 1st line = 2nd line ?
    l1 = α(0) * Φ(t) + quadgk(σ -> d(α, σ) * Φ(t - σ), 0, t)[1]
    l2_printed = α(t) * Φ(0) - quadgk(σ -> dΦ(t - σ) * α(σ), 0, t)[1]        # "- ∫ ∂Φ(t-σ)/∂t α(σ)dσ"  (l. 1902)
    l2_fixed = α(t) * Φ(0) + quadgk(σ -> dΦ(t - σ) * α(σ), 0, t)[1]
    C("eq:alpha_E_conv 2a linha", "α(0)Φ(t)+∫α̇Φ(t-σ)dσ = α(t)Φ(0) - ∫ [∂Φ(t-σ)/∂t] α dσ  AS PRINTED (∂/∂t of the argument)", l1, l2_printed; kind=:typo)
    C("eq:alpha_E_conv corrigida", "... = α(t)Φ(0) + ∫ Φ'(t-σ) α(σ)dσ  [= -∫ ∂Φ(t-σ)/∂σ α dσ, integration by parts]", l1, l2_fixed)
    diag(@sprintf("eq:alpha_E_conv: lhs=%.6f  printed(-)=%.6f  corrected(+)=%.6f   (Φ(0)=%.3f, Wagner)", l1, l2_printed, l2_fixed, Φ(0)))

    # eq:cnaw_pp_E_recurrence + eq:state_space_indicial_pp: E = QS - ΣX_i, Ẋ_i = -b_iΘ X_i + A_i Q̇, X_i(0)=A_i Q(0)
    Q(t) = 0.7 + 0.4 * sin(0.9t) + 0.1t                                       # (c_nα w)^QS
    Qd(t) = d(Q, t)
    conv = Q(0) * Φ(t) + quadgk(σ -> Qd(σ) * Φ(t - σ), 0, t)[1]               # eq:cnaw_pp_E_conv
    X = rk4((τ, x) -> -bb .* Θ .* x .+ A .* Qd(τ), A .* Q(0), 0.0, t, 4000)
    C("eq:cnaw_pp_E_recurrence", "(c_nα w)^E = QS - ΣX_i, with Ẋ_i = -b_iΘX_i + A_iQ̇ (X_i(0)=A_iQ(0)) equals the convolution of eq:cnaw_pp_E_conv", Q(t) - sum(X), conv; tol=1e-9)
    B("X_i(0) = A_i Q(0)", "initial condition that makes E(0)=Q(0)(1-ΣA)=Q(0)Φ(0); if instead X(0)=0 the Duhamel initial term is lost", abs((Q(0) - sum(A .* Q(0))) - Q(0) * Φ(0)) < 1e-14)

    # eq:recursive_formula: X^n = X^{n-1} e^{-bΘΔt} + A ΔQ  -> order of accuracy
    texact = 3.6
    Xex = rk4((τ, x) -> -bb .* Θ .* x .+ A .* Qd(τ), A .* Q(0), 0.0, texact, 20000)
    function recur(Δt; mid=false)
        n = round(Int, texact / Δt)
        X = A .* Q(0)
        for k in 1:n
            ΔQ = Q(k * Δt) - Q((k - 1) * Δt)
            X = X .* exp.(-bb .* Θ .* Δt) .+ A .* ΔQ .* (mid ? exp.(-bb .* Θ .* Δt / 2) : 1.0)
        end
        return X
    end
    e1(Δt; kw...) = norm(recur(Δt; kw...) - Xex)
    p_print = log2(e1(0.2) / e1(0.1)); p_print2 = log2(e1(0.1) / e1(0.05))
    p_mid = log2(e1(0.2; mid=true) / e1(0.1; mid=true)); p_mid2 = log2(e1(0.1; mid=true) / e1(0.05; mid=true))
    B("eq:recursive_formula ordem 1", "printed recurrence (increment A_iΔQ without e^{-b_iΘΔt/2}) is FIRST-order accurate in Δt (observed order ≈ 1)", abs(p_print - 1) < 0.15 && abs(p_print2 - 1) < 0.15; kind=:note)
    B("recorrência c/ ponto médio", "with the factor e^{-b_iΘΔt/2} on the increment (Leishman) the recurrence is second-order (observed ≈ 2)", abs(p_mid - 2) < 0.3 && abs(p_mid2 - 2) < 0.3)
    diag(@sprintf("recurrence orders: printed %.2f/%.2f ; with midpoint factor %.2f/%.2f", p_print, p_print2, p_mid, p_mid2))

    # eq:pre_state_space_pp_indicial algebra and the Euler limit
    Δt = 0.1; x0 = [0.3, 0.2]; ΔQ = 0.07
    lhs = (x0 .* exp.(-bb .* Θ .* Δt) .+ A .* ΔQ .- x0) ./ Δt
    rhs = (-x0 .* (1 .- exp.(-bb .* Θ .* Δt)) .+ A .* ΔQ) ./ Δt
    C("eq:pre_state_space", "(X^n - X^{n-1})/Δt = [-X^{n-1}(1-e^{-bΘΔt}) + AΔQ]/Δt (algebra)", lhs, rhs; tol=1e-13)
    push!(UNVERIFIED, "eq:state_space_indicial_pp with time-varying U, M (Θ=Θ(t), c_nα=c_nα(M(t))): the convolution with Φ(t-σ) of constant Θ is replaced by the state equation with Θ(t) (Beddoes-Leishman / Jose et al. practice); this is an approximation, not an identity, when U or M vary.")

    # sums of A_i: Φ(0)
    C("Φ_p(0) Wagner", "Σ A_i^C = 0.5 => Φ_p^C(0) = 0.5 (Wagner's value at t=0)", 1 - sum(WAGNER[1]), 0.5)
    B("Beddoes/Jose: ΣA=1", "for the compressible sets Σ A_i^C = 1, i.e. Φ_p^C(0)=0 (the initial response is purely non-circulatory)", abs(sum(BEDDOES[1]) - 1) < 1e-12 && abs(sum(JOSE[1]) - 1) < 1e-12)
    B("[code] parametros", "src/AeroSolver.jl has Wagner (0.165,0.335;0.0455,0.3) [text: l. 1943] and Küssner (0.5,0.5;0.13,1.0) [l. 2166]",
        code_has(CODE_SOLV, ":Wagner  => (SVector(0.165, 0.335), SVector(0.0455, 0.3))") && code_has(CODE_SOLV, ":Kussner       => (SVector(0.5, 0.5), SVector(0.13, 1.0))"))
    B("l. 1906 'rate of change'", "text: 'Φ defines the RATE OF CHANGE of an output in response to a step change in an input'. Φ is the step RESPONSE itself (Φ' is the impulse response)", false; kind=:note)
end

# =============================================================================
# 2.1.1  kinematics: quasi-steady normalwash   [l. 1920-1925, 1989-1996, 2084-2088 + section 2.5 equivalences]
# =============================================================================
function sec_kinematics()
    s = "2.1-cinemat"
    C(l, dsc, a, b; kw...) = check(s, l, dsc, a, b; kw...)
    B(l, dsc, c; kw...) = checkb(s, l, dsc, c; kw...)
    # explicit 2-D rigid-body kinematics. Inertial axes: X aft, Z up. The airfoil flies toward -X at speed U in still air; its pitch axis
    # moves DOWN by h(t) (h>0 down, Theodorsen's convention). Body axes: x_c aft along the chord, z_c normal (up); pitch θ nose-up:
    # x_c = (cosθ, -sinθ), z_c = (sinθ, cosθ).
    U0 = 23.0; b = 0.45; a = -0.37
    h(t) = 0.1 * t^2 + 0.05 * sin(2t); θ(t) = 0.2 + 0.35 * sin(1.1t)            # large pitch: no small-angle assumption
    xc(t) = [cos(θ(t)), -sin(θ(t))]; zc(t) = [sin(θ(t)), cos(θ(t))]
    Xb(t) = [-U0 * t, -h(t)]
    XP(t, ξ) = Xb(t) + ξ * xc(t)                                           # point at distance ξ aft of the pitch axis, on the chordline
    t = 0.8
    rel_wind(ξ) = -ForwardDiff.derivative(τ -> XP(τ, ξ), t)                  # air velocity relative to the point (air at rest)
    hd = d(h, t); θd = d(θ, t)
    w_pp = U0 * sin(θ(t)) + hd * cos(θ(t)); v_pp = U0 * cos(θ(t)) - hd * sin(θ(t))
    # eq:v_pp_QS, eq:w_pp_QS: components of the relative wind AT THE PITCH AXIS in the body axes
    C("eq:v_pp_QS", "v^QS_pp = U cosθ - ḣ sinθ = chordwise component of the relative wind at the pitch axis (ḣ>0 downward)", dot(rel_wind(0.0), xc(t)), v_pp; tol=1e-12)
    C("eq:w_pp_QS", "w^QS_pp = U sinθ + ḣ cosθ = normal component of the relative wind at the pitch axis", dot(rel_wind(0.0), zc(t)), w_pp; tol=1e-12)
    # eq:w_34pp_QS: 3/4-chord is (1/2-a)b aft of the pitch axis ; eq:w_12pp_QS: 1/2-chord is -ab aft
    C("eq:w_34pp_QS", "w^QS_{3/4} = U sinθ + ḣ cosθ + b(1/2-a)θ̇  (exact, no small angles)", dot(rel_wind((1 / 2 - a) * b), zc(t)), w_pp + b * (1 / 2 - a) * θd; tol=1e-12)
    C("eq:w_12pp_QS", "w^QS_{1/2} = U sinθ + ḣ cosθ - abθ̇", dot(rel_wind(-a * b), zc(t)), w_pp - a * b * θd; tol=1e-12)
    # if ḣ were positive UP the printed formulas would be wrong => the text must state the sign of ḣ
    hdup = -hd
    C("sinal de ḣ", "with ḣ POSITIVE UP the printed w_{3/4} = U sinθ + ḣ cosθ + ... would NOT hold (it holds for ḣ positive DOWN = +W3 direction)", dot(rel_wind((1 / 2 - a) * b), zc(t)), U0 * sin(θ(t)) + hdup * cos(θ(t)) + b * (1 / 2 - a) * θd; kind=:note, tol=1e-6)
    # eq:alpha_pp_QS
    C("eq:alpha_pp_QS", "α^QS = arctan(w/v) = angle between chordline and relative wind at the pitch axis", atan(dot(rel_wind(0.0), zc(t)), dot(rel_wind(0.0), xc(t))), atan(w_pp, v_pp); tol=1e-12)
    # section 2.5 equivalences: U_P3 = U3 - Ω1 c (s - sP), a = 2s-1 ; U3 ≡ w_pp, Ω1 ≡ θ̇
    c = 2b; sP_ = (a + 1) / 2
    for (sP, nm, w) in ((3 / 4, "3/4", w_pp + b * (1 / 2 - a) * θd), (1 / 2, "1/2", w_pp - a * b * θd))
        C("U_P3 ≡ w_$nm", "U_P3 = U_3 - Ω_1 c (s - s_P) (eq:Usp3) with s=(a+1)/2, U_3≡w_pp, Ω_1≡θ̇ gives w^QS at s_P=$nm", w_pp - θd * c * (sP_ - sP), w; tol=1e-12)
    end
    push!(UNVERIFIED, "The sign convention of ḣ (positive in the +W3 direction, i.e. DOWNWARD for a wing seen from above) is not stated in the text; it is the one that makes eq:w_pp_equiv (w_pp ≡ U_3) consistent.")
end

# =============================================================================
# 2.1.1  inertial loads (incompressible and compressible)   [l. 1972-2088]
# =============================================================================
function sec_inertial()
    s = "2.1-inercial"
    C(l, dsc, a, b; kw...) = check(s, l, dsc, a, b; kw...)
    B(l, dsc, c; kw...) = checkb(s, l, dsc, c; kw...)
    rng = MersenneTwister(21)
    # ---- incompressible: eq:cnI_pp_inc, eq:cmI_pp_inc  vs Theodorsen (1935), non-circulatory part (h DOWN, nose-up θ, a aft of midchord):
    #   L_up = πρb²(ḧ + Uθ̇ - abθ̈),   M_nc (nose-up, about the axis) = πρb²[ abḧ - Ub(1/2-a)θ̇ - b²(1/8+a²)θ̈ ]
    for trial in 1:3
        a = 2rand(rng) - 1; b = 0.2 + rand(rng); U = 10 + 30rand(rng); hdd, θd, θdd = randn(rng, 3)
        wdot = hdd + U * θd - a * b * θdd                                   # d/dt of w_{1/2} (small angles, U const)
        cn_text = π * b / U^2 * wdot
        cm_text = -π * b / (2U^2) * (-a * wdot + b * θdd / 8 + U * θd / 2)
        L = π * b^2 * (hdd + U * θd - a * b * θdd); M = π * b^2 * (a * b * hdd - U * b * (1 / 2 - a) * θd - b^2 * (1 / 8 + a^2) * θdd)
        trial == 1 && C("eq:cnI_pp_inc", "c_n^I = (πb/U²) ẇ_{1/2} = non-circulatory lift of Theodorsen / (ρU²b)", cn_text, L / (U^2 * b); tol=1e-12)
        trial == 1 && C("eq:cmI_pp_inc", "c_m^I = -(πb/2U²)[-aẇ_{1/2} + bθ̈/8 + Uθ̇/2] = Theodorsen's non-circulatory moment about the axis / (2ρU²b²)", cm_text, M / (2U^2 * b^2); tol=1e-12)
        trial > 1 && C("eq:cnI/cmI", "idem (random trial $trial)", [cn_text, cm_text], [L / (U^2 * b), M / (2U^2 * b^2)]; tol=1e-12)
    end
    # ---- compressible: initial (t=0+) values from piston theory: Δp = 2ρā w(x), w(x) = Uα + θ̇ (x - ab), x aft of midchord in [-b,b]
    a = 0.31; b = 0.4; U = 120.0; ā = 340.0; M = U / ā; α = 0.07; θd = 0.9
    w(x) = U * α + θd * (x - a * b)
    cn_pist = quadgk(x -> 2 * ā * w(x), -b, b)[1] / (U^2 * b)                            # ρ cancels: Δp/(½ρU²·2b)
    cm_pist = -quadgk(x -> 2 * ā * w(x) * (x - a * b), -b, b)[1] / (2U^2 * b^2)         # nose-up moment: force aft of the axis is nose-down
    C("piston c_n(0+)", "c_n(0+) = 4α/M - 4abθ̇/(UM)  (coefficients of eq:cnI_comp_conv)", cn_pist, 4α / M - 4 * a * b * θd / (U * M); tol=1e-12)
    C("piston c_m(0+)", "c_m(0+) = 2aα/M - 2b(a²+1/3)θ̇/(UM)  (coefficients of eq:cmI_comp_conv: NEGATIVE sign of the θ̇ term)", cm_pist, 2a * α / M - 2 * b * (a^2 + 1 / 3) * θd / (U * M); tol=1e-12)
    # M term: step ΔU at fixed α: Δp = 2ρā α ΔU ; normalized by ½ρU²·2b
    ΔM = 0.01; ΔU = ΔM * ā
    C("piston termo Ṁ (cn)", "step ΔM at fixed α: Δc_n(0+) = +4αΔM/M² (sign + of eq:cnI_comp_conv second term)", 4 * ā * α * ΔU / U^2, 4α * ΔM / M^2; tol=1e-12)
    # M→0 consistency of the last term of eq:cmI_comp: -πbθ̇/(4Uβ) -> incompressible Uθ̇/2 term of eq:cmI_pp_inc
    C("limite M→0 (termo -πbθ̇/4Uβ)", "last term of eq:cmI_comp at β=1 equals the  -(πb/2U²)(Uθ̇/2) term of eq:cmI_pp_inc", -π * b * θd / (4U), -π * b / (2U^2) * (U * θd / 2); tol=1e-14)

    # ---- convolution (printed) vs state-space recast (printed): initial states steady, X_j(0) = τ_j u_j(0)
    Ma0 = 0.25; τ0 = 0.013                                                             # time scales (s)
    Tnα, TnM, Tnθ̇, Tmα, TmM, Tmθ̇ = 0.011, 0.012, 0.013, 0.014, 0.0145, 0.016
    αt(t) = 0.05 + 0.04 * sin(30t) + 0.1t; Mt(t) = Ma0 + 0.03 * sin(20t) + 0.2t; θdt(t) = 0.3 + 0.5 * sin(25t); Ut(t) = 80.0 + 5 * sin(10t)
    a = 0.2; b = 0.3; t = 0.07
    # exponential kernels  ∫ u̇(σ) e^{-(t-σ)/τ} dσ  (conv) and  u - X/τ  (recast)
    conv(u, τ) = quadgk(σ -> d(u, σ) * exp(-(t - σ) / τ), 0, t; rtol=1e-12)[1]
    Xs(u, τ) = τ * u(0) * exp(-t / τ) + quadgk(σ -> u(σ) * exp(-(t - σ) / τ), 0, t; rtol=1e-12)[1]   # exact solution of Ẋ=-X/τ+u, X(0)=τu(0)
    rec(u, τ) = u(t) - Xs(u, τ) / τ
    # c_n  (eq:cnI_comp_conv -> eq:cnI_comp)
    cn_conv = 4 / Mt(t) * conv(αt, Tnα) + 4 * αt(t) / Mt(t)^2 * conv(Mt, TnM) + (-4 * a * b / (Ut(t) * Mt(t))) * conv(θdt, Tnθ̇)
    cn_rec = 4 / Mt(t) * rec(αt, Tnα) + 4 * αt(t) / Mt(t)^2 * rec(Mt, TnM) - 4 * a * b / (Ut(t) * Mt(t)) * rec(θdt, Tnθ̇)
    C("eq:cnI_comp_conv→cnI_comp", "the 3 terms of c_n^I (printed convolution form) = recast with states (steady initial states; the α(0)Φ, M(0)Φ, θ̇(0)Φ terms absent in the printed conv. are exactly what X(0)=τu(0) accounts for)", cn_conv, cn_rec; tol=1e-9)
    # c_m, term by term
    convm(u, T) = AI[1] * conv(u, bI[1] * T) + AI[2] * conv(u, bI[2] * T)              # ∫u̇ Φ_m(t-σ)dσ, Φ_m = A1 e^{-t/(b1 T)} + A2 e^{-t/(b2 T)}
    recm(u, T) = AI[1] * rec(u, bI[1] * T) + AI[2] * rec(u, bI[2] * T)                 # (u - A1X4/(b1T) - A2X5/(b2T)) since A1+A2=1
    cm_conv_α = -(-2a / Mt(t)) * convm(αt, Tmα); cm_rec_α = 2a / Mt(t) * recm(αt, Tmα)
    cm_conv_M = -(-2a * αt(t) / Mt(t)^2) * convm(Mt, TmM); cm_rec_M = 2a * αt(t) / Mt(t)^2 * recm(Mt, TmM)
    cm_conv_θ = -(-2b * (a^2 + 1 / 3) / (Ut(t) * Mt(t))) * conv(θdt, Tmθ̇)               # AS PRINTED (l. 1981): "-∫ (-2b(a²+1/3)/(UM)) θ̈ Φ"
    cm_rec_θ = -2b * (a^2 + 1 / 3) / (Ut(t) * Mt(t)) * rec(θdt, Tmθ̇)                  # eq:cmI_comp (l. 2045) and src/Aerodynamics.jl
    C("eq:cmI_comp_conv termo α", "-∫(-2a/M) α̇ Φ_mα = +(2a/M)(α - A1X4/(b1T) - A2X5/(b2T))", cm_conv_α, cm_rec_α; tol=1e-9)
    C("eq:cmI_comp_conv termo Ṁ", "-∫(-2aα/M²) Ṁ Φ_mM = +(2aα/M²)(M - A1X6/(b1T) - A2X7/(b2T))", cm_conv_M, cm_rec_M; tol=1e-9)
    C("eq:cmI_comp_conv termo θ̈", "-∫(-2b(a²+1/3)/(UM)) θ̈ Φ_mθ̇ [AS PRINTED, = +…] vs recast -(2b(a²+1/3)/(UM))(θ̇ - X8/T_mθ̇): OPPOSITE sign", cm_conv_θ, cm_rec_θ; kind=:typo)
    C("eq:cmI_comp_conv corrigida", "with the integral term written  -∫ (2b(a²+1/3)/(UM)) θ̈ Φ_mθ̇  (same sign as the c_n term and as the piston-theory value) conv = recast", -(2b * (a^2 + 1 / 3) / (Ut(t) * Mt(t))) * conv(θdt, Tmθ̇), cm_rec_θ; tol=1e-9)
    # code consistency
    B("[code] cnI compressivel", "src: cnI = 4/Ma*(α-χ1/Tnα) + -4aₕ b/(Uᵢ Ma)*(Ωₐ-χ2/Tnθ̇) + 4α/Ma²*(Ma-χ6/TnM)", code_has(CODE_AERO, "4/Ma * (α - χ1/Tnα) + -4*aₕ*b/(Uᵢ*Ma) * (Ωₐ - χ2/Tnθ̇) + 4*α/Ma^2 * (Ma - χ6/TnM)"))
    B("[code] cmI compressivel", "src: cmI = 2aₕ/Ma(α-…) + -2b(aₕ²+1/3)/(Uᵢ Ma)(Ωₐ-χ5/Tmθ̇) + 2aₕα/Ma²(Ma-…) - πbΩₐ/(4Uᵢβₚ)  (negative sign of the θ̇ term)", code_has(CODE_AERO, "+ -2*b*(aₕ^2+1/3)/(Uᵢ*Ma) * (Ωₐ - χ5/Tmθ̇)") && code_has(CODE_AERO, "- π*b*Ωₐ/(4*Uᵢ*βₚ)"))
    B("[code] cnI/cmI incompressivel", "src: cnI = πbUₙdotMid/Uᵢ² ; cmI = -πb/(2Uᵢ²)(-aₕ UₙdotMid + bΩₐdot/8 + UᵢΩₐ/2)", code_has(CODE_AERO, "π*b*UₙdotMid/Uᵢ^2") && code_has(CODE_AERO, "-π*b/(2*Uᵢ^2) * (-aₕ * UₙdotMid + b * Ωₐdot / 8 + Uᵢ * Ωₐ / 2)"))
    # state equations and their numbering
    B("eq:xdot5_cmI_comp", "text: Ẋ_8^I = -X_8^I/T_{nθ̇} + θ̇ ; but X_8 enters c_m through (θ̇ - X_8/T_{mθ̇}) (eq:cmI_comp) and the code uses 1/Tmθ̇", false; kind=:typo)
    B("[code] ordem dos estados", "src: B[inertial] = (α, Ωₐ, α, α, Ωₐ, Ma, Ma, Ma): ORDER (nα, nθ̇, mα, mα, mθ̇, nM, mM, mM) differs from the text X1..X8 (nα, nM, nθ̇, mα, mα, mM, mM, mθ̇)", !(code_has(CODE_AERO, "SVector(α, Ωₐ, α, α, Ωₐ, Ma, Ma, Ma)")); kind=:note)

    # ---- time scales eqs. (T_n_alpha ... T_m_thetadot)
    c = 0.6; āa = 340.0
    # initial-slope matching with piston theory:  d/dt[c_n^I + c_n^C]_{0+} = -(4/M)(1-M) ā/c   (α step; c_nα = 2π/β, Φ = 1 - ΣA e^{-bβ²Ut/b})
    for (nm, (Ak, bk)) in (("Wagner", WAGNER), ("Beddoes", BEDDOES), ("Jose", JOSE))
        M = 0.4; β = sqrt(1 - M^2); Σ = sum(Ak .* bk); U = M * āa; bsc = c / 2
        T_print = (c / āa) / ((1 - M) + π * M^2 * β^2 * Σ)                          # eq:T_n_alpha AS PRINTED (β²), k=1
        T_code = (c / āa) / ((1 - M) + π * β * M^2 * Σ)                             # src/Aerodynamics.jl (β), = Leishman K_α
        slope(T) = -(4 / M) / T + (2π / β) * Σ * β^2 * U / bsc                      # total initial slope: non-circulatory + circulatory
        target = -(4 / M) * (1 - M) * āa / c
        C("T_nα ($nm) printed", "initial slope of c_n^I+c_n^C with T_nα AS PRINTED (πM²β²ΣAb) vs piston-theory slope -(4/M)(1-M)ā/c", slope(T_print), target; kind=:typo, tol=1e-6)
        C("T_nα ($nm) código", "idem with πβM²ΣAb (code, Leishman & Beddoes K_α): matches", slope(T_code), target; tol=1e-12)
    end
    B("[code] T_nα", "src: Tnα = Tᵢ / (1 - Ma + π * βₚ * Ma^2 * ΣγACbC)  (β, not β²)", code_has(CODE_AERO, "Tnα = Tᵢ / (1 - Ma + π * βₚ * Ma^2 * ΣγACbC)"))
    B("[code] T_nθ̇, T_mα, T_mθ̇", "src formulas = eq:T_n_thetadot (2πM²β), eq:T_m_alpha, eq:T_m_thetadot (7/[15(1-M)+3πM²β A3 b3]); k=0.75 (l. 2037)",
        code_has(CODE_AERO, "Tnθ̇ = Tᵢ / (1 - Ma + 2π * βₚ * Ma^2 * ΣγACbC)") && code_has(CODE_AERO, "Tmα = Tᵢ * (AI[1]*bI[2] + AI[2]*bI[1]) / (bI[1] * bI[2] * (1-Ma))") &&
        code_has(CODE_AERO, "Tmθ̇ = 7*Tᵢ / (15*(1-Ma) + 3π * βₚ * Ma^2 * AI[3]*bI[3])") && code_has(CODE_AERO, "λ = 0.75"))
    B("[code] T_nM", "src: TnM = Tᵢ / (1 - Ma + π / βₚ * Ma^2 * ΣγACbC) = eq:T_n_M (β^{-1}); TmM = Tmα", code_has(CODE_AERO, "TnM = Tᵢ / (1 - Ma + π / βₚ * Ma^2 * ΣγACbC)") && code_has(CODE_AERO, "TmM = Tmα"))
    # T_mα: initial relative slope of Φ_mα
    M = 0.4; Tmα_ = (c / āa) * (AI[1] * bI[2] + AI[2] * bI[1]) / (bI[1] * bI[2] * (1 - M))
    C("T_mα", "Φ_mα'(0)/Φ_mα(0) = -(A1/b1+A2/b2)/T_mα = -(1-M)ā/c (same piston-theory relative decay rate as for c_n)", -(AI[1] / bI[1] + AI[2] / bI[2]) / Tmα_, -(1 - M) * āa / c; tol=1e-12)
    B("eq:T_m_M (símbolo)", "text l. 2031 labels the 5th time scale 'T_{n_M}' a 2nd time; it is T_{m_M} (used in eq:Phi_m_M_I and eq:cmI_comp)", false; kind=:typo)
    push!(UNVERIFIED, "eq:T_n_M (T_{n_M} with πM²β^{-1}ΣAb) and the M-rate terms (4α/M²)Ṁ, (2aα/M²)Ṁ: they come from Jose et al. (2006); only the sign/initial value of the M-step response was checked here (piston theory).")
    push!(UNVERIFIED, "Target of the initial-slope matching, -(4/M)(1-M)ā/c, is Leishman & Beddoes' (1989) piston-theory result (propagation of the TE disturbance at relative speed ā-U); not re-derived here.")
    push!(UNVERIFIED, "In the compressible inertial loads U, M, α^QS multiplying the convolutions (e.g. (4α/M²) before ∫Ṁ Φ) are frozen at the current time (approximation of the recast).")
end

# =============================================================================
# 2.1.2  Peters finite-state inflow  [l. 2093-2111]  + Wagner (Jones) approximation
# =============================================================================
function sec_inflow()
    s = "2.1-inflow"
    C(l, dsc, a, b; kw...) = check(s, l, dsc, a, b; kw...)
    B(l, dsc, c; kw...) = checkb(s, l, dsc, c; kw...)
    ks = [0.01, 0.05, 0.1, 0.3, 0.5, 1.0, 2.0, 3.0, 5.0]
    Cex = [Cth(k) for k in ks]
    # Bessel validation
    x = 0.7
    C("Bessel (Wronskiano)", "J1 Y0 - J0 Y1 = 2/(πx): validates the quadrature Bessel functions used for Theodorsen's C(k)", Jn(1, x) * Yn(0, x) - Jn(0, x) * Yn(1, x), 2 / (π * x); tol=1e-10)
    C("C(0.1) conhecido", "Theodorsen C(0.1) = 0.8319 - 0.1723i (tables)", Cth(0.1), 0.8319 - 0.1723im; tol=2e-4)
    # Peters: A λ̇ + λ = ẇ c (U/b=1), w^E = w - ½ bᵀλ  =>  w^E/w^QS = C_P(ik) = 1 - ½ ik bᵀ (ik A + I)^{-1} c  (text eqs. 2097-2108)
    for N in (4, 6, 8)
        A, Ainv, bv, cv = inflow_arrays(N)
        CP(k) = 1 - 0.5 * im * k * dot(bv, (im * k * A + I) \ cv)
        CPwrong(k) = 1 + 0.5 * im * k * dot(bv, (im * k * A + I) \ cv)
        err = maximum(abs(CP(k) - Cex[i]) for (i, k) in enumerate(ks))
        errw = maximum(abs(CPwrong(k) - Cex[i]) for (i, k) in enumerate(ks))
        C("Peters N=$N vs Theodorsen", "w^E/w^QS = 1 - ½bᵀλ̂ with Aλ̇ + (U/b)λ = ẇ^QS c reproduces Theodorsen's C(k) (k∈[0.01,5]); max|ΔC| = $(round(err, sigdigits=2)) [tol 0.04]", err, 0.0; tol=0.04)
        diag(@sprintf("Peters N=%d: max|C_P - C_Theodorsen| = %.4f ; C_P(k→0)=%.6f ; C_P(k→∞)=%.5f ; mutation (w^E = w^QS + λ_0): max err %.3f", N, err, real(CP(1e-6)), real(CP(1e3)), errw))
        N == 6 && C("Peters N=6 limites", "C_P(0)=1 and C_P(∞)=1/2 (lift deficiency limits)", [real(CP(1e-6)), real(CP(1e4))], [1.0, 0.5]; tol=5e-3)
    end
    A6, A6inv, b6, c6 = inflow_arrays(6)
    C("[code] inflow_arrays", "copy of inflow_arrays equals AeroBeams.inflow_arrays (N=6)", HAVE_AB ? vcat(vec(A6), b6, c6) : zeros(1), HAVE_AB ? vcat(vec(AeroBeams.inflow_arrays(6)[1]), AeroBeams.inflow_arrays(6)[3], AeroBeams.inflow_arrays(6)[4]) : zeros(1); tol=1e-12)
    # Wagner/Jones exponential approximation vs Theodorsen in the frequency domain: C_J(ik) = 1 - ΣA_i ik/(ik+b_i)
    CJ(k) = 1 - sum(WAGNER[1][i] * im * k / (im * k + WAGNER[2][i]) for i in 1:2)
    errJ = maximum(abs(CJ(k) - Cex[i]) for (i, k) in enumerate(ks))
    C("Jones vs Theodorsen", "Wagner approx (0.165,0.335;0.0455,0.3) vs exact C(k): max|ΔC| = $(round(errJ, sigdigits=2)) [tol 0.02]", errJ, 0.0; tol=0.02)
    diag(@sprintf("Jones (2-term Wagner) max|C_J - C_Theodorsen| = %.4f", errJ))
    # text vs code on the time scale of the inflow equation
    B("[code] Θ no Peters", "src: Ac .= -Θ*AₚInvVec, A[...] = -Θ*AₚInv with Θ = Uᵢ βₚ²/b ; text eq:lambda_statespace has (U/b)λ (no β²)", !(code_has(CODE_AERO, "A[circulatoryPitchPlungeStatesRange,circulatoryPitchPlungeStatesRange] .= -Θ * AₚInv") && code_has(CODE_AERO, "Θ = Uᵢ/b*βₚ²")); kind=:note)
    B("[code] w^E Peters", "src: wₑ = UₙTQC - 1/2 dot(bₚ, λ); B = UₙdotTQC * AₚInv cₚ  (λ̇ = -ΘA⁻¹λ + ẇ A⁻¹ c): consistent with eq:w_34_eff_peters, eq:lambda_0, eq:lambda_statespace", code_has(CODE_AERO, "return UₙTQC - 1/2*dot(ppSolver.bₚ, view(χ,circulatoryPitchPlungeStatesRange))") && code_has(CODE_AERO, "B[circulatoryPitchPlungeStatesRange] .= UₙdotTQC * AₚInvcₚ"))
end

# =============================================================================
# 2.1.3  flap loads  [l. 2116-2152] , 2.1.6 ad hoc  [l. 2249-2275]
# =============================================================================
function sec_flap()
    s = "2.1-flap"
    C(l, dsc, a, b; kw...) = check(s, l, dsc, a, b; kw...)
    B(l, dsc, c; kw...) = checkb(s, l, dsc, c; kw...)
    # thin-airfoil (Birnbaum/Glauert) equivalent normalwash: w_{3/4}^QS = (1/π) ∫_{-1}^{1} w(x) sqrt((1+x)/(1-x)) dx, x = cosφ-like substitution
    # x = cosφ: sqrt((1+x)/(1-x)) dx → (1+cosφ)dφ (smooth)
    a = -0.4
    I_pitch = quadgk(φ -> (cos(φ) - a) * (1 + cos(φ)), 0, π)[1] / π
    C("w_{3/4} (arfagem)", "(1/π)∫(x-a)√((1+x)/(1-x))dx = 1/2 - a: the pitch-rate part b(1/2-a)θ̇ of eq:w_34pp_QS is the thin-airfoil equivalent normalwash", I_pitch, 1 / 2 - a; tol=1e-12)
    for dh in (-0.3, 0.4, 0.8)
        φd = acos(dh)
        T10 = sqrt(1 - dh^2) + acos(dh); T11 = (2 - dh) * sqrt(1 - dh^2) + (1 - 2dh) * acos(dh)
        I0 = quadgk(φ -> 1 + cos(φ), 0, φd)[1]                                       # ∫_d^1 √((1+x)/(1-x)) dx
        I1 = quadgk(φ -> (cos(φ) - dh) * (1 + cos(φ)), 0, φd)[1]                     # ∫_d^1 (x-d) √((1+x)/(1-x)) dx
        C("T10 (d=$dh)", "T10 = ∫_d^1 √((1+x)/(1-x))dx: w_{3/4,f} = (T10/π)Uδ (eq:w_34f_QS, 1st term)", I0, T10; tol=1e-12)
        C("T11 (d=$dh)", "T11/2 = ∫_d^1 (x-d)√((1+x)/(1-x))dx: w_{3/4,f} ∋ (b T11/2π) δ̇ (eq:w_34f_QS, 2nd term)", I1, T11 / 2; tol=1e-12)
        if HAVE_AB && dh == 0.4
            Th = AeroBeams.theodorsen_flap_constants(0.3, 0.7)
            C("[code] Th[10],Th[11]", "AeroBeams.theodorsen_flap_constants: T10, T11 equal Theodorsen's closed forms", [Th[10], Th[11]], [T10, T11]; tol=1e-12)
            C("[code] w_f", "src: wf = Th[17] Uᵢ δ + b Th[18] δ̇ with Th[17]=T10/π, Th[18]=T11/(2π) = eq:w_34f_QS", [Th[17], Th[18]], [T10 / π, T11 / (2π)]; tol=1e-12)
        end
    end
    # quasi-steady flap lift slope: c_n = c_nα T10 δ/π = 2 T10 δ = 2(π - θ_h + sinθ_h) δ (classical thin-airfoil flap effectiveness)
    dh = 0.4; θh = acos(-dh)
    C("dcl/dδ", "2π·(T10/π) = 2T10 = 2(π - θ_h + sinθ_h), θ_h = hinge angle from the LE (classical thin-airfoil flap effectiveness)", 2 * (sqrt(1 - dh^2) + acos(dh)), 2 * (π - θh + sin(θh)); tol=1e-12)
    # eq:cmI_f dimension: printed has an overall b/(2U²) in front of [(T4+T10)U²δ + (...)Ubδ̇ + (...)b²δ̈]: not dimensionless
    Th = [0.3, -0.2, 0.1, 0.5, 0.2, -0.2, 0.15, 0.1, 0.05, 2.0, 1.3, 0.9, 0.2, 0, 0, 0, 0, 0]     # arbitrary numbers: only the structure matters
    T1, T4, T7, T8, T10, T11 = Th[1], Th[4], Th[7], Th[8], Th[10], Th[11]
    aa, dd = -0.3, 0.2; εm = 0.96; δ, δd, δdd = 0.1, 0.5, 2.0
    cm_f(b, U, printed) = -εm * (printed ? b : 1.0) / (2U^2) * ((T4 + T10) * U^2 * δ + (T1 - T8 - T4 * (dd - aa) + T11 / 2) * U * b * δd + (T7 + T1 * (dd - aa)) * b^2 * δdd)
    # dimensional analysis: scale lengths by λ with time fixed: b→λb, U→λU ; δ̇, δ̈ unchanged (rates in 1/T). c_m must be invariant.
    λ = 3.0
    C("eq:cmI_f (dimensões)", "printed c_m^f (prefactor b/2U²): NOT invariant under a change of length unit (b→λb, U→λU): factor λ", cm_f(0.4λ, 20λ, true), cm_f(0.4, 20.0, true); kind=:typo)
    C("eq:cmI_f corrigida", "c_m^f = -ε_{f_m}/(2U²)[…] (no b in front, as in src/Aerodynamics.jl) is dimensionless", cm_f(0.4λ, 20λ, false), cm_f(0.4, 20.0, false); tol=1e-12)
    B("[code] cmI flap", "src: cmI -= ϵₘ/(2*Uᵢ^2) * (Uᵢ^2*Th[14]*δ + Uᵢ*b*Th[15]*δdot + b^2*Th[16]*δddot)  (no b in front)", code_has(CODE_AERO, "cmI -= ϵₘ/(2*Uᵢ^2) * (Uᵢ^2*Th[14]*δ + Uᵢ*b*Th[15]*δdot + b^2*Th[16]*δddot)"))
    B("[code] cnI flap", "src: cnI -= ϵₙ*b/Uᵢ^2 * (Uᵢ*Th[4]*δdot + b*Th[1]*δddot) = eq:cnI_f", code_has(CODE_AERO, "cnI -= ϵₙ*b/Uᵢ^2 * (Uᵢ*Th[4]*δdot + b*Th[1]*δddot)"))
    # efficiency also multiplies the circulatory flap normalwash in the code
    B("[code] ε_n na parte circulatória", "src: wₑ = wₑp + ϵₙ*wₑf + wₑg: ε_{f_n} also scales the CIRCULATORY flap normalwash; eq:alpha_E (l. 2191) and the text of l. 2152 ('ad hoc efficiency factors') present it only in the inertial coefficients", !(code_has(CODE_AERO, "wₑ = wₑp + ϵₙ*wₑf + wₑg")); kind=:note)
    # eq:w_34f_E and state equations vs code
    B("[code] flap states", "src: wₑf = wf - Σχ_f/cnα ; Ẋ = -bΘX + A ċ_nαw_f (B = cnαwfdot*ACf) = eq:w_34f_E, eq:state_space_indicial_flap", code_has(CODE_AERO, "wₑf = wf - sum(view(χ, flapStatesRange)) / cnα") && code_has(CODE_AERO, "B[flapStatesRange] .= cnαwfdot * ACf"))
    # ad hoc / TableLookup version (2.1.6)
    B("[code] TableLookup: c_t", "src: ctC -= cdδ*abs(δ)/cos(α): flap drag increment ∝ |δ| (even), whereas eq:ctfC_adhoc is linear, c_{t_f}^C = c_{tδ}δ_f", !(code_has(CODE_AERO, "ctC -= cdδ * abs(δ) / cos(α)")); kind=:note)
    B("[code] TableLookup: c_n→c_m,c_t", "src adds cnδ δ INSIDE cnC (cnC += cnδ*δ) BEFORE cmC = cnC*(referenceAxis-x_ac) and ctC = η cnC tan(…): the flap normal force also produces moment about the axis (arm s-x_ac) and tangential force; eq:cm_attached_alternative/eq:ct_attached_alternative show c_m^C (from the flapless c_n^C) + c_{m_f}^C and c_t^C + c_{t_f}^C", !(code_has(CODE_AERO, "cnC += cnδ * δ")); kind=:note)
    push!(UNVERIFIED, "Theodorsen flap constants T1..T9, T12, T13 and the non-circulatory flap loads (eq:cnI_f, eq:cmI_f): the closed forms in src/AeroProperties.jl coincide with Theodorsen's (1935) as remembered, and T10, T11 were verified by thin-airfoil integrals; the apparent-mass flap loads were not re-derived.")
end

# =============================================================================
# 2.1.4  gust-induced loads   [l. 2157-2184]
# =============================================================================
function sec_gust()
    s = "2.1-gust"
    C(l, dsc, a, b; kw...) = check(s, l, dsc, a, b; kw...)
    B(l, dsc, c; kw...) = checkb(s, l, dsc, c; kw...)
    for (nm, Ag, bg) in (("Küssner", [0.5, 0.5], [0.13, 1.0]), ("Berci-Righi", [0.3694, 0.0550, 0.2654, 0.1829, 0.0861, 0.0412], [0.3733, 0.0179, 0.1096, 1.6003, 10.428, 170.93]))
        C("ΣA^G = 1 ($nm)", "Φ_g^C(0) = 1 - ΣA_i^G = 0 (used in l. 2173 to drop the first term of eq:w_34g_conv)", sum(Ag), 1.0; tol=1e-4)
    end
    Ag, bg = [0.5, 0.5], [0.13, 1.0]; Θ = 0.9
    Φg(t) = 1 - sum(Ag .* exp.(-bg .* Θ .* t)); dΦg(t) = sum(Ag .* bg .* Θ .* exp.(-bg .* Θ .* t))
    wqs(t) = 0.8 + 0.5 * sin(1.7t) + 0.1t
    t = 2.9
    # eq:w_34g_conv as printed vs eq:w_34g_E (state form)
    conv_printed = wqs(t) * Φg(0) - quadgk(σ -> dΦg(t - σ) * wqs(σ), 0, t)[1]             # l. 2159 ( - ∫ ∂Φ(t-σ)/∂t w dσ )
    conv_fixed = wqs(t) * Φg(0) + quadgk(σ -> dΦg(t - σ) * wqs(σ), 0, t)[1]
    X = [quadgk(σ -> exp(-bg[i] * Θ * (t - σ)) * wqs(σ), 0, t)[1] for i in 1:2]             # eq:gust_states
    wE_states = Θ * sum(Ag .* bg .* X)                                                     # eq:w_34g_E
    C("eq:w_34g_conv printed", "w^E = w^QS Φ_g(0) - ∫ ∂Φ_g(t-σ)/∂t w^QS dσ AS PRINTED vs eq:w_34g_E (states): opposite sign", conv_printed, wE_states; kind=:typo)
    C("eq:w_34g_conv corrigida", "w^E = w^QS Φ_g(0) + ∫ Φ_g'(t-σ) w^QS dσ = βU/b... = Θ ΣA_i b_i X_i (eq:w_34g_E): agree", conv_fixed, wE_states; tol=1e-10)
    # state equation (eq:state_space_gust)
    Xode = rk4((τ, x) -> -bg .* Θ .* x .+ wqs(τ), zeros(2), 0.0, t, 6000)
    C("eq:state_space_gust", "Ẋ_i^G = -b_iΘX_i + w^QS (X(0)=0) integrates to eq:gust_states", Xode, X; tol=1e-9)
    # step gust: w^E(t)/w0 must rise from 0 to 1 (Küssner)
    Xstep(t) = (1 .- exp.(-bg .* Θ .* t)) ./ (bg .* Θ)
    C("rampa do gust (degrau)", "unit step gust: w^E(0)=0, w^E(∞)=1 (Küssner-type response)", [Θ * sum(Ag .* bg .* Xstep(0.0)), Θ * sum(Ag .* bg .* Xstep(1e4))], [0.0, 1.0]; tol=1e-9)
    # code
    B("[code] gust", "src: wₑg = Θ*dot(AGbG, χ_g) [AGbG = A∘b] and B_g = UₙGust, Ag = -Θ*bG = eq:w_34g_E, eq:state_space_gust", code_has(CODE_AERO, "Θ * dot(AGbG, view(χ, gustStatesRange))") && code_has(CODE_AERO, "B[gustStatesRange] .= UₙGust") && code_has(CODE_AERO, "Ag .= - Θ * bG"))
    B("[code] velocidade relativa com gust", "src: α_E = atan(-wₑ/(Uₜ+UₜGust)) with Uₜ = U_2, v_g^QS ≡ -U_g2, v^QS ≡ -U_2: = eq:alpha_E; but U_i (dynamic pressure, Mach, Re, Θ) uses only -V (gust excluded: 'quasi-stationary gust' assumption, l. 2184)", code_has(CODE_AERO, "αₑ = atan( - wₑ / (Uₜ + UₜGust) )"))
    push!(UNVERIFIED, "The Küssner function approximation (0.5,0.5;0.13,1.0) and the Berci-Righi coefficients are taken from the literature; only ΣA=1 and the algebra of the state-space form were verified.")
end

# =============================================================================
# 2.1.5  composition of loads   [l. 2189-2244]
# =============================================================================
function sec_composition()
    s = "2.1-composicao"
    C(l, dsc, a, b; kw...) = check(s, l, dsc, a, b; kw...)
    B(l, dsc, c; kw...) = checkb(s, l, dsc, c; kw...)
    rng = MersenneTwister(5)
    # eq:alpha_E vs code
    we, v, vg = 0.7, 22.0, 1.3
    Ut, Ug = -v, -vg                                                                       # U_2 = -v ; U_g2 = -v_g
    C("eq:alpha_E ↔ código", "arctan(w^E/(v+v_g)) = atan(-wₑ/(Uₜ+UₜGust)) with v ≡ -U_2, v_g ≡ -U_g2", atan(we / (v + vg)), atan(-we / (Ut + Ug)); tol=1e-14)
    # force decomposition: c_l, c_d from c_n, c_t (c_t positive toward the LEADING edge), α = angle between chord and relative wind
    cn, ct, α = 0.9, 0.04, 0.17
    # relative wind direction (air velocity) in chord axes: (cosα, sinα): flow comes from slightly below for α>0
    ŵ = [cos(α), sin(α)]; l̂ = [-sin(α), cos(α)]
    F = [-ct, cn]                                                                          # F = N ẑ + T_aft x̂, T_aft = -c_t (c_t positive toward LE)
    C("eq:cl_attached", "c_l = c_n cosα + c_t sinα = force ⟂ relative wind (c_t positive toward the LE)", dot(F, l̂), cn * cos(α) + ct * sin(α); tol=1e-14)
    C("eq:cd_attached", "c_d = c_n sinα - c_t cosα = force ∥ relative wind", dot(F, ŵ), cn * sin(α) - ct * cos(α); tol=1e-14)
    # eq:ct0: pure parasite drag  (c_n=0)
    cd0 = 0.011
    C("eq:ct0", "c_t0 = -c_d0/cosα gives c_d = c_d0 when c_n = 0", -(-cd0 / cos(α)) * cos(α), cd0; tol=1e-14)
    # eq:ctLI: c_d = c_n²/(π AR) from c_t^LI = -c_n²/(π AR cosα)
    AR = 8.0
    C("eq:ctLI", "c_t^{LI} = -c_n²/(π AR cosα) gives c_d = c_n²/(π AR) (c_di = c_l²/(π AR))", -(-cn^2 / (π * AR * cos(α))) * cos(α), cn^2 / (π * AR); tol=1e-14)
    # eq:cmC_attached
    for x_ac in (0.25, 0.259)
        aₕ = -0.3; sp = (aₕ + 1) / 2; c = 0.6; Fn = 0.8
        # explicit moment about the attachment point b of a normal force at the AC; b2 points to the LE; s, x_ac measured from the LE
        r = [0.0, (sp - x_ac) * c, 0.0]                                                    # AC - b (LE direction = +b2)
        m = cross(r, [0.0, 0.0, Fn])[1] / (Fn * c)                                         # nose-up (+b1) moment / (F c)
        printed = x_ac + aₕ / 2
        if x_ac == 0.25
            C("eq:cmC_attached (x_ac=1/4)", "c_m^C = c_n^C (x_ac + a/2) coincides with the exact lever arm (s - x_ac) only for x_ac = 1/4", printed, m; tol=1e-12)
        else
            C("eq:cmC_attached (x_ac=$x_ac)", "c_m^C = c_n^C (x_ac + a/2) AS PRINTED vs exact lever arm (s - x_ac) = a/2 + 1/2 - x_ac, x_ac measured from the LE (as in src/Airfoil.jl x_ac=0.259 for some airfoils)", printed, m; kind=:typo)
            C("eq:cmC_attached corrigida", "c_m^C = c_n^C (a/2 + 1/2 - x_ac)", aₕ / 2 + 1 / 2 - x_ac, m; tol=1e-12)
        end
    end
    B("[code] cmC", "src: cmC = cnC * (referenceAxis - x_ac), x_ac = fraction of the chord from the LE (Airfoil.jl)", code_has(CODE_AERO, "cmC = cnC * (referenceAxis - x_ac)"))
    B("[code] ct", "src: ctC = η*cnC*tan(αₑ+b*Ωₐ/(2Uᵢ)) ; ct₀ = -cd₀/cos(α) ; induced: ctC -= cnC²/(π AR)/cos(α)", code_has(CODE_AERO, "η * cnC * tan(αₑ+b*Ωₐ/(2*Uᵢ))") && code_has(CODE_AERO, "-cd₀/cos(α)") && code_has(CODE_AERO, "ctC -= cnC^2 / (π*AR) / cos(α)"))
    B("[code] cn circulatório", "src: cnC = cnα sin(αₑ)cos(αₑ) (smallAngles=false) = eq:cnC_attached", code_has(CODE_AERO, "cnα * sin(αₑ) * cos(αₑ)"))
    push!(UNVERIFIED, "eq:ctC_attached: c_t^C = η c_n^C tan(α^E + bθ̇/(2U)) (Peters 2007 eq. 41 / Leishman eq. 8.173): the argument α^E + bθ̇/2U (leading-edge suction angle) could not be re-derived independently here.")
    push!(UNVERIFIED, "The use of α^QS (not α^E) to resolve c_n, c_t into c_l, c_d and in c_t0, c_t^{LI} is conceptually right (geometric angle between chord and wind) and agrees with the code.")
end

# =============================================================================
# derivative of c_nα (feeds the circulatory input of the states)
# =============================================================================
function sec_cnalpha()
    s = "2.1-cnα"
    C(l, dsc, a, b; kw...) = check(s, l, dsc, a, b; kw...)
    B(l, dsc, c; kw...) = checkb(s, l, dsc, c; kw...)
    c0 = 2π; ā = 340.0
    cn(M) = c0 / sqrt(1 - M^2)                                                             # eq:cnalpha
    U(t) = 100.0 + 30t; M(t) = U(t) / ā
    true_rate = d(t -> cn(M(t)), 0.0)
    Ui = U(0.0); Uid = 30.0; Ma = M(0.0); βp = sqrt(1 - Ma^2)
    βpdot = -Ui * Uid / ((Ui / Ma)^2 * βp)                                                 # as in src/Aerodynamics.jl cnα_rate
    C("β̇ (código)", "βₚdot = -Uᵢ U̇ᵢ/((Uᵢ/Ma)² βₚ) = d/dt sqrt(1-M²)", βpdot, d(t -> sqrt(1 - M(t)^2), 0.0); tol=1e-12)
    code_old = cn(Ma) * βpdot / βp^2                                                       # BEFORE the audit: cnαdot = cnα * βₚdot / βₚ²
    code_rate = - cn(Ma) * βpdot / βp                                                      # src (corrected): cnαdot = - cnα * βₚdot / βₚ
    C("ċ_nα (código antigo)", "BEFORE the audit src cnα_rate had cnαdot = cnα βₚdot/βₚ² vs d/dt[c_nα0/β(M(t))] = -c_nα β̇/β: opposite sign and an extra 1/β", code_old, true_rate; kind=:note)
    C("ċ_nα corrigida (código)", "src cnα_rate (corrected) = -cnα βₚdot/βₚ = d/dt[c_nα0/β(M(t))]", code_rate, true_rate; tol=1e-12)
    B("[code] cnα_rate", "src: cnαdot = - cnα * βₚdot / βₚ (corrected)", code_has(CODE_AERO, "cnαdot = - cnα * βₚdot / βₚ"))
    diag(@sprintf("cnα rate: true=%.5f  old code=%.5f (ratio %.4f = -1/β = %.4f)  corrected code=%.5f", true_rate, code_old, code_old / true_rate, -1 / βp, code_rate))
end

# =============================================================================
println("="^100)
println("Verification of section 2.1 (Attached-flow aerodynamics) of theory/main.tex   -   AeroBeams loaded: $HAVE_AB")
println("="^100)
sec_indicial()
sec_kinematics()
sec_inertial()
sec_inflow()
sec_flap()
sec_gust()
sec_composition()
sec_cnalpha()

println()
let lastsec = ""
    for r in RESULTS
        mark = r.kind == :pass ? (r.ok ? "ok" : "FALHOU") :
               r.kind == :typo ? (r.ok ? "ERRO NO TEXTO CONFIRMADO" : "!! esperava-se discrepancia e nao houve") :
                                 (r.ok ? "NOTA (afirmacao imprecisa/nao geral/inconsistente com o code, ver relatorio)" : "!! esperava-se diferenca e nao houve")
        if r.sec != lastsec
            println("\n------------------------- ", r.sec, " -------------------------")
            lastsec = r.sec
        end
        @printf("%-14s %-30s %-170s %10.3e  %s\n", r.sec, first(r.label, 30), first(r.desc, 170), r.err, mark)
    end
end
println("\n--- DIAGNOSTICS (numbers quoted in the report) ---")
foreach(println, DIAG)
println("\n--- NOT VERIFIABLE NUMERICALLY ---")
foreach(x -> println(" * ", x), UNVERIFIED)
nok = count(r -> r.kind == :pass && r.ok, RESULTS)
nbad = count(r -> r.kind == :pass && !r.ok, RESULTS)
ntypo = count(r -> r.kind == :typo && r.ok, RESULTS)
nnote = count(r -> r.kind == :note && r.ok, RESULTS)
nmiss = count(r -> r.kind != :pass && !r.ok, RESULTS)
println("\nSUMMARY: ", nok, " checks passed, ", nbad, " FAILED, ", ntypo, " errors-in-text confirmed (formula as printed disagrees), ", nnote, " notes, ", nmiss, " expected discrepancies that did not occur")
