# =============================================================================
# Conceptual + numerical verification of section 2.4 ("Gust field profiles") of theory/main.tex
# =============================================================================
# Companion of verify_chapter1.jl / verify_C.jl (same conventions: check(...) kinds pass / typo / note).
#
# Method
#   (1) The formulas are re-implemented exactly AS PRINTED in main.tex (the script also reads main.tex and asserts that the
#       printed formula is there. After the audit the \\fix{} corrections were applied to main.tex, so the "[texto]" rows now confirm that the corrected wording is present).
#   (2) They are compared with the implementation of src/Gust.jl (replicas of the code recipes written here, and, when AeroBeams
#       can be loaded, the real functions create_*Gust / stochastic_gust_velocity_from_white_noise).
#   (3) Independent numerical tests: quadrature of the PSDs (normalization, isotropy relation), Parseval for the sum of
#       sinusoids, exact expectation of the white-noise coloring (a real(ifft(fft(x).*H)) filter acts as the Hermitian filter
#       (H_k + H_{N-k})/2, so the generated PSD is known in closed form) validated by Monte Carlo (fixed seeds),
#       root-bending-moment argument for the DARPA profile.
#
#   Result categories (last column of the report)
#     ok                         agreement (statement/step of the text is correct)
#     ERRO NO TEXTO CONFIRMADO   the formula evaluated AS PRINTED disagrees (typo); the corrected one passes
#     NOTA                       statement imprecise / not general / inconsistent with the implementation
#   For the "checkv" rows the "err" column is the value of the metric (see the description), not a relative error.
#
# Run:   julia --project=<AeroBeams root> theory/verification/verify_2_4.jl      (~15 s)
# Last run (2026-11-14, after the code corrections): see SUMMARY line; before: 43 checks passed, 0 failed, 2 errors-in-text confirmed (ΔΩ), 12 notes (text imprecisions / code discrepancies).
# =============================================================================

using LinearAlgebra, QuadGK, Random, Printf, FFTW, ForwardDiff

const TOL = 1e-9

# ----------------------------------------------------------------------------- real code (optional)
const CODE_GUST = read(joinpath(@__DIR__, "..", "..", "src", "Gust.jl"), String)
const HAVE_AB = try
    @eval using AeroBeams
    true
catch err
    @warn "AeroBeams could not be loaded: checks against the real implementation are skipped" exception = err
    false
end

# ----------------------------------------------------------------------------- printed text of the section
const TXT = let L = readlines(joinpath(@__DIR__, "..", "main.tex"))
    i0 = findfirst(l -> occursin("\\section{Gust field profiles}", l), L)
    i1 = findfirst(l -> occursin("\\section{Kinematical description in terms of beam variables}", l), L)
    join(L[i0:i1-1], "\n")
end
printed(p) = occursin(p, TXT)

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
# metric version: :pass -> value < thr ; :typo/:note -> value > thr (the claim "value < thr" is false). err column = value
function checkv(sec, label, desc, value, thr; kind=:pass)
    ok = kind == :pass ? value < thr : value > thr
    push!(RESULTS, Res(sec, label, desc, Float64(value), kind, ok))
    return ok
end
const DIAG = String[]          # numbers quoted in the report
diag(args...) = push!(DIAG, string(args...))
const UNVERIFIED = String[]    # statements that cannot be checked numerically

# ----------------------------------------------------------------------------- PSDs (as printed / as coded)
const aVK = 1.339
Ψu(Ω; σ=1.0, L=762.0, a=aVK) = σ^2 * 2L / π / (1 + (a * L * Ω)^2)^(5 / 6)
Ψw(Ω; σ=1.0, L=762.0, a=aVK) = σ^2 * L / π * (1 + 8 / 3 * (a * L * Ω)^2) / (1 + (a * L * Ω)^2)^(11 / 6)
Dstd_u(Ω; σ=1.0, L=762.0) = σ^2 * 2L / π / (1 + (L * Ω)^2)                  # Dryden longitudinal (MIL-F-8785C)
Dcode_u(Ω; σ=1.0, L=762.0) = σ^2 * 2L / π / (1 + (L * Ω)^2)                 # as coded in src/Gust.jl (exponent 1 after the audit; it was 2)
D_w(Ω; σ=1.0, L=762.0) = σ^2 * L / π * (1 + 3 * (L * Ω)^2) / (1 + (L * Ω)^2)^2

# =============================================================================
# (a) von Kármán PSDs
# =============================================================================
function sec_psd()
    s = "2.4-PSD"
    C(l, d, a, b; kw...) = check(s, l, d, a, b; kw...)
    B(l, d, c; kw...) = checkb(s, l, d, c; kw...)
    B("[texto] eq:vkpsd_vw", "printed Ψu and Ψv=Ψw (2L/π, 5/6, L/π, 8/3, 11/6, 1.339)",
        printed("\\frac{2L}{\\pi}\\frac{1}{[1+(1.339L\\Omega)^2]^{5/6}}") &&
        printed("\\frac{L}{\\pi}\\frac{1+\\frac{8}{3}(1.339L\\Omega)^2}{[1+(1.339L\\Omega)^2]^{11/6}}"))
    Iu = quadgk(Ω -> Ψu(Ω), 0, Inf; rtol=1e-12)[1]
    Iw = quadgk(Ω -> Ψw(Ω), 0, Inf; rtol=1e-12)[1]
    aex = 2 / π * quadgk(x -> (1 + x^2)^(-5 / 6), 0, Inf; rtol=1e-12)[1]
    C("a=1.339", "constant 1.339 = (2/π)∫₀^∞(1+x²)^(-5/6)dx (exact 1.33902...) [rounding only]", aex, aVK; tol=1e-4)
    C("∫₀^∞Ψu=σ²", "one-sided normalization of Ψu (σ=1, L=762 m)", Iu, 1.0; tol=1e-4)
    C("∫₀^∞Ψw=σ²", "one-sided normalization of Ψv=Ψw", Iw, 1.0; tol=1e-4)
    C("∫_{-∞}^{∞}Ψu=2σ²", "=> the printed Ψ are ONE-SIDED PSDs (per rad/m, Ω≥0)", quadgk(Ω -> Ψu(Ω), -Inf, Inf; rtol=1e-12)[1], 2.0; tol=1e-4)
    B("[texto] convenção", "text now states that Ψ is one-sided (∫₀^∞Ψ dΩ = σ²) [added as \\fix{} after the audit]", printed(r"one-sided|single-sided|unilateral"))
    Ωs = [1e-4, 5e-4, 1e-3, 4e-3, 2e-2]
    dΨu(Ω) = ForwardDiff.derivative(x -> Ψu(x), Ω)
    C("isotropia vK", "Ψw = ½(Ψu − Ω dΨu/dΩ) (isotropic turbulence, same L): 8/3 and 11/6 are consistent with 2L/π and 5/6",
        [Ψw(Ω) for Ω in Ωs], [0.5 * (Ψu(Ω) - Ω * dΨu(Ω)) for Ω in Ωs]; tol=1e-12)
    C("L=2500 ft", "2500 ft = 762 m", 2500 * 0.3048, 762.0; tol=1e-12)
    # Dryden spectrum implemented in src/Gust.jl (not in the text)
    dDu(Ω) = ForwardDiff.derivative(x -> Dstd_u(x), Ω)
    C("isotropia Dryden", "Dryden: Φw = ½(Φu − ΩΦu') holds with the standard Φu (denominator exponent 1)", D_w.(Ωs), 0.5 .* (Dstd_u.(Ωs) .- Ωs .* dDu.(Ωs)); tol=1e-12)
    IDc = quadgk(Ω -> Dcode_u(Ω), 0, Inf; rtol=1e-12)[1]
    C("∫Φu Dryden (padrão)", "standard Dryden Φu integrates to σ²", quadgk(Ω -> Dstd_u(Ω), 0, Inf; rtol=1e-12)[1], 1.0; tol=1e-8)
    C("∫Φw Dryden", "Dryden Φv=Φw integrates to σ²", quadgk(Ω -> D_w(Ω), 0, Inf; rtol=1e-12)[1], 1.0; tol=1e-8)
    C("∫Φu Dryden (código)", "[CÓDIGO corrigido] src/Gust.jl Dryden Φu now has denominator exponent 1 (it was 2, which integrated to σ²/2): integrates to σ²", IDc, 1.0; tol=1e-8)
    checkb(s, "[código] Dryden Φu", "src/Gust.jl has Φu = 2Lσ²/π/(1+(LΩ)²) in both Dryden generators (no ^2)", count("σ^2*2*L/π * 1 / (1+(L*Ω)^2)\n", CODE_GUST) == 1 && count("σ^2*2*L/π * 1 / (1+(LΩ)^2)\n", CODE_GUST) == 1)
    diag("∫₀^∞Ψu=", Iu, " ∫₀^∞Ψw=", Iw, " a_exato=", aex, " ∫Dryden-u(código)=", IDc)
    push!(UNVERIFIED, "MIL-STD-1797A: 'low altitude' = below 1000 ft (not 2000 ft) and σ_w = 0.1 u20 refers to the VERTICAL component; for h<1000 ft L_w=h (not 2500 ft). From memory, not checked against the standard.")
end

# =============================================================================
# (d) discrete profiles
# =============================================================================
function sec_discrete()
    s = "2.4-discreto"
    C(l, d, a, b; kw...) = check(s, l, d, a, b; kw...)
    B(l, d, c; kw...) = checkb(s, l, d, c; kw...)
    B("[texto] eq:wgse/wg1c/DARPA", "printed: w0[h(t)-h(t-2bH/U)], (1/2)w0[1-cos(πUt/(bH))], -(1/2)w0[...]cos(2πs̄)",
        printed("w_0\\lb h(t)-h \\lp t-\\frac{2bH}{U}\\rp\\rb") && printed("\\cos{\\lp\\frac{\\pi Ut}{bH}\\rp}") && printed("-\\frac{1}{2} w_0 \\lb1-\\cos"))
    B("[texto] H", "text now says H is the gust gradient distance in semichords (2H = length) [\\fix{} after the audit]", printed("gust gradient distance"))
    w0, b, H, U = 1.7, 0.05, 12.0, 25.0
    dur = 2b * H / U
    h(x) = x ≥ 0 ? 1.0 : 0.0
    w1c_text(t) = 0.5 * w0 * (1 - cos(π * U * t / (b * H))) * (h(t) - h(t - 2b * H / U))
    w1c_code(t) = 0 < t < dur ? 0.5 * w0 * (1 - cos(2π * t / dur)) : 0.0       # create_OneMinusCosineGust (initialTime=0)
    ts = range(0.013, dur - 0.011, length=17)
    C("eq:wg1c vs código", "one-minus-cosine: text (πUt/(bH)) = code (2πt/duration) with duration = 2bH/U", w1c_text.(ts), w1c_code.(ts); tol=1e-12)
    C("1-cos: pico", "maximum w0 at t = bH/U", w1c_text(b * H / U), w0; tol=1e-12)
    dw(t) = ForwardDiff.derivative(x -> 0.5 * w0 * (1 - cos(π * U * x / (b * H))), t)
    C("1-cos: continuidade C¹", "w(0)=w(2bH/U)=0 and dw/dt=0 at both ends (smooth entry/exit)", [w1c_text(0.0), w1c_text(dur - 1e-12), dw(0.0), dw(dur)], zeros(4); tol=1e-9)
    C("duração = 2bH/U", "gust length = U·duration = 2bH = 2H semichords", U * dur / b, 2H; tol=1e-12)
    B("H = comprimento", "claim 'H is the length of the gust in semichords' (it is the gradient distance = HALF the length)", isapprox(U * dur / b, H); kind=:note)
    # sharp-edged
    C("eq:wgse", "sharp-edged: w0 for 0≤t<2bH/U (top-hat of length 2H semichords), 0 elsewhere",
        [w0 * (h(t) - h(t - dur)) for t in (0.3dur, 0.9dur, 1.1dur)], [w0, w0, 0.0]; tol=1e-12)
    # DARPA: spatial factor
    zeronet = quadgk(x -> cos(2π * x), -0.5, 0.5)[1]
    C("DARPA: força líquida nula", "∫_{-1/2}^{1/2} cos(2πs̄) ds̄ = 0 (no net lift: pure bending excitation)", zeronet, 0.0; tol=1e-12)
    Mroot(k) = quadgk(x -> x * (-cos(2π * k * x)), 0, 0.5)[1]          # root bending moment of ONE semispan, nondimensional, text sign
    Ms = [Mroot(k) for k in 1:10]
    C("DARPA: momento na raiz", "M_root = ∫₀^{1/2} s̄ w ds̄ for w=-cos(2π s̄) is +1/(2π²)", Ms[1], 1 / (2π^2); tol=1e-12)
    checkb(s, "DARPA: k=1 é o máximo", "among w=-cos(2πk s̄) (zero net lift) k=1 maximizes |M_root| (text: 'sinusoidal ... maximal bending moment at the root')", argmax(abs.(Ms)) == 1)
    C("DARPA vs limite", "|M_root| = 8/π² (81%) of the bang-bang bound 1/16 for |w|≤1 with zero net lift", Ms[1] / (1 / 16), 8 / π^2; tol=1e-12)
    # sign text vs code
    ℓ, wd = 30.0, 6.0
    wdar_text(t, sb) = -0.5 * w0 * (1 - cos(π * U * t / (b * H))) * cos(2π * sb)
    wdar_code(x1, x2) = -0.5 * w0 * cos(2π * x1 / wd) * (1 - cos(2π * x2 / ℓ))      # create_DiscreteSpaceGust(:DARPA), c0 = 0 (minus sign added after the audit)
    pts = [(0.31, 0.2), (0.0, 0.5), (-0.4, 0.77), (0.1, 0.5)]
    tt(l) = l * 2b * H / U                                                          # Ut/(2bH) = x2/ℓ
    C("DARPA: |texto| = |código|", "same modulus at corresponding points (t↔x2/ℓ = Ut/(2bH), s̄↔x1/wd)", [abs(wdar_text(tt(l), sb)) for (sb, l) in pts], [abs(wdar_code(sb * wd, l * ℓ)) for (sb, l) in pts]; tol=1e-12)
    C("DARPA: sinal texto vs código", "text and code now agree: minus sign (negative/minimum at the center s̄=0, positive/maximum at the tips s̄=±1/2: also the sign of figures/DARPA_gust.pdf)",
        [wdar_text(tt(l), sb) for (sb, l) in pts], [wdar_code(sb * wd, l * ℓ) for (sb, l) in pts]; tol=1e-12)
    # real code
    if HAVE_AB
        g = AeroBeams.create_DiscreteSpaceGust(type=:DARPA, gustLength=ℓ, gustWidth=wd, verticalVelocity=w0, c0=zeros(3))
        rng = MersenneTwister(24)
        xs = [[wd * (rand(rng) - 0.5), ℓ * rand(rng), 0.0] for _ in 1:8]
        C("código real: DARPA", "AeroBeams.create_DiscreteSpaceGust(:DARPA) = replica (negative at the center)", [g.UGustInertial(x, 0.0)[3] for x in xs], [wdar_code(x[1], x[2]) for x in xs]; tol=1e-12)
        g1 = AeroBeams.create_OneMinusCosineGust(initialTime=0.0, duration=dur, verticalVelocity=w0)
        C("código real: 1-cos", "AeroBeams.create_OneMinusCosineGust = eq:wg1c with duration=2bH/U", [g1.UGustInertial(t)[3] for t in ts], w1c_text.(ts); tol=1e-12)
        g2 = AeroBeams.create_SharpEdgedGust(initialTime=0.0, duration=dur, verticalVelocity=w0)
        C("código real: sharp-edged", "AeroBeams.create_SharpEdgedGust = eq:wgse", [g2.UGustInertial(t)[3] for t in (0.3dur, 0.9dur, 1.1dur)], [w0, w0, 0.0]; tol=1e-12)
    end
    # convective velocity (code): profile defined w.r.t. c0 but the window follows c(t) = c0 + Vc t
    ℓ2, Vc, tn = 10.0, 2.0, 2.0
    front = Vc * tn
    if HAVE_AB
        gc = AeroBeams.create_DiscreteSpaceGust(type=:OneMinusCosine, gustLength=ℓ2, gustWidth=4.0, convectiveVelocity=Vc, verticalVelocity=1.0)
        wfront = gc.UGustInertial([0.0, front + 1e-9, 0.0], tn)[3]
        wmid = gc.UGustInertial([0.0, front + ℓ2 / 2, 0.0], tn)[3]
    else
        l(x2) = (x2 - front) / ℓ2          # replica of `l(x,t) = (x[2]-c(t))/gustLength` (corrected: the profile follows the moving window)
        wfront = 0.5 * (1 - cos(2π * l(front + 1e-9)))
        wmid = 0.5 * (1 - cos(2π * l(front + ℓ2 / 2)))
    end
    checkb(s, "convectiveVelocity≠0", "[CÓDIGO corrigido] DiscreteSpaceGust(1-cos) moving window: w=0 at the gust front and w=w0 at mid-gust (code gives $(round(wfront, digits=3)) and $(round(wmid, digits=3)); before the audit 0.905 and 0.095)",
        abs(wfront) < 1e-6 && abs(wmid - 1) < 1e-6)
    diag("convective 1-cos: w(front)=", wfront, " w(mid)=", wmid, "  |  DARPA M_root(k=1..4)=", Ms[1:4])
end

# =============================================================================
# (b) sum of sinusoids
# =============================================================================
sin_grid(; T, ωmax=200π, ωmin=0.0) = (dω = min(2π / T, (2π)^2 / ωmax); (collect(ωmin:dω:ωmax), dω))    # as in create_Continuous1DGust
function sin_var(Ψ; U, T, L=762.0, σ=1.0, ωmax=200π, mode=:code)
    ω, dω = sin_grid(; T, ωmax)
    wt = ones(length(ω))
    mode == :skipDC && (wt[1] = 0.0)
    mode == :halfDC && (wt[1] = 0.5)
    return sum(wt .* Ψ.(ω ./ U; σ, L)) * dω / U
end
function sec_sinusoids()
    s = "2.4-senoides"
    C(l, d, a, b; kw...) = check(s, l, d, a, b; kw...)
    B(l, d, c; kw...) = checkb(s, l, d, c; kw...)
    B("[texto] ΔΩ", "text now has ΔΩ = 2π/(UT) [\\fix{} after the audit]", printed("\\Delta\\Omega = \\fix{2\\pi/(UT)}"))
    B("[texto] φ_i", "text now uses independent phases φ_i^(u), φ_i^(v), φ_i^(w) [\\fix{} after the audit]", printed("\\varphi_i^{(u)}") && printed("\\varphi_i^{(w)}"))
    # dimensional analysis: (exponent of m, exponent of s)
    dΨ = (3, -2)                 # (m/s)²/(rad/m)
    dΩt = (1 - 0, -1 - 1)        # 2πU/T: (m/s)/s
    dΩc = (-1, 0)                # 2π/(U T): 1/((m/s) s)
    dvar = (2, -2)               # (m/s)²
    B("ΔΩ dimensional (texto)", "Ψ·ΔΩ has units (m/s)² with ΔΩ = 2πU/T (as printed): m/s² ≠ rad/m", (dΨ[1] + dΩt[1], dΨ[2] + dΩt[2]) == dvar; kind=:typo)
    B("ΔΩ dimensional (corrigido)", "Ψ·ΔΩ has units (m/s)² with ΔΩ = 2π/(UT)", (dΨ[1] + dΩc[1], dΨ[2] + dΩc[2]) == dvar)
    # numerical consequence: grid of spacing Δω=2π/T (spatial ΔΩ_c) but amplitudes weighted with the printed ΔΩ
    U, T, L, σ = 25.0, 3000.0, 762.0, 1.0
    dΩ_c = 2π / (U * T); dΩ_t = 2π * U / T
    Ω = collect(0:dΩ_c:8.0)
    Vc = sum(Ψw.(Ω)) * dΩ_c; Vt = sum(Ψw.(Ω)) * dΩ_c * (dΩ_t / dΩ_c)
    C("RMS com ΔΩ impresso", "weights 2Ψ·ΔΩ with ΔΩ=2πU/T on the frequency grid of spacing 2π/(UT) multiply the RMS by U", sqrt(Vt / Vc), U; tol=1e-12)
    checkb(s, "RMS = σ?", "RMS of the generated gust equals σ using the printed ΔΩ (it is U·σ)", abs(sqrt(Vt) - σ) < 0.05σ; kind=:typo)
    # Ω_i U t equals ω_i t
    ωi = 0.37; C("Ω_i U t", "phase Ω_i U t = ω_i t (code uses ω t)", (ωi / U) * U * 3.1, ωi * 3.1; tol=1e-14)
    # Parseval: time average of the realized signal over a common period equals sum of amplitude² / 2 (+ DC)
    ωmax, Tp = 20π, 50.0
    ω, dω = sin_grid(; T=Tp, ωmax)
    Uref = 25.0; dΩ = dω / Uref
    amp = sqrt.(2 .* dΩ .* Ψw.(ω ./ Uref))
    ψ = 2π .* rand(MersenneTwister(7), length(ω))
    Tcommon = 2π / dω
    Ns = 8192
    tt = (0:Ns-1) .* (Tcommon / Ns)
    sig(t) = sum(amp .* cos.(ω .* t .+ ψ))
    msq = sum(abs2, sig.(tt)) / Ns
    pars = sum(amp[2:end] .^ 2) / 2 + amp[1]^2 * cos(ψ[1])^2
    C("Parseval", "time-mean-square of Σ√(2ΨΔΩ)cos(ωt+φ) over a common period = Σ amp²/2 (+ DC term amp₀² cos²φ₀)", msq, pars; tol=1e-9)
    # expected variance (ensemble) = Σ_i Ψ(Ω_i)ΔΩ (i = 0 included with full weight) vs σ²
    for (T, U) in ((100.0, 25.0), (1000.0, 25.0), (10.0, 25.0))
        vu, vw = sin_var(Ψu; U, T), sin_var(Ψw; U, T)
        du, dw_ = sin_var(Ψu; U, T, mode=:skipDC), sin_var(Ψw; U, T, mode=:skipDC)
        hu, hw = sin_var(Ψu; U, T, mode=:halfDC), sin_var(Ψw; U, T, mode=:halfDC)
        diag(@sprintf("sinusoids U=%g T=%g: Var/σ² (code) u=%.3f w=%.3f | skipDC u=%.3f w=%.3f | halfDC u=%.3f w=%.3f   [ΔΩL=%.2f]", U, T, vu, vw, du, dw_, hu, hw, L * min(2π / T, 4π^2 / (200π)) / U))
    end
    vw100 = sin_var(Ψw; U=25.0, T=100.0); vu100 = sin_var(Ψu; U=25.0, T=100.0)
    hw100 = sin_var(Ψw; U=25.0, T=100.0, mode=:halfDC); hu100 = sin_var(Ψu; U=25.0, T=100.0, mode=:halfDC)
    C("Var/σ² (w, T=100 s, código corrigido)", "[CÓDIGO corrigido] Ω=0 term with half weight (src/Gust.jl ampWeight): w variance within 6% of σ² (U=25 m/s, T=100 s, L=762 m); with full weight it was 1.26", hw100, 1.0; tol=0.06)
    C("Var/σ² (u, T=100 s, código corrigido)", "same for u: still ≈12% above σ² because ΔΩ·1.339L is not ≪ 1 (the caveat now stated in the text)", hu100, 1.0; kind=:note, tol=0.05)
    C("Var/σ² (T=1000 s)", "T=1000 s (ΔΩ·1.339 L ≈ 0.3): variance within 5% of σ² (w), half-weight DC", sin_var(Ψw; U=25.0, T=1000.0, mode=:halfDC), 1.0; tol=0.05)
    checkb(s, "[código] DC meio peso", "src/Gust.jl: ampWeight[1] = 1/sqrt(2) when ωmin == 0 (Ω=0 term given half weight)", occursin("ωmin == 0 && (ampWeight[1] = 1/sqrt(2))", CODE_GUST))
    # same phase for u and w => correlated components
    Ωg = collect(0:2π/(25.0 * 5000.0):8.0); dΩg = Ωg[2] - Ωg[1]
    ρsame = sum(sqrt.(Ψu.(Ωg) .* Ψw.(Ωg))) * dΩg / sqrt(sum(Ψu.(Ωg)) * dΩg * sum(Ψw.(Ωg)) * dΩg)
    checkv(s, "correlação u–w (mesma fase)", "text uses ONE phase φ_i for u, v, w: correlation coefficient u–w (isotropic turbulence along a line: 0; code uses independent phases)", ρsame, 0.05; kind=:note)
    diag("rho_uw(same phase)=", ρsame)
    push!(UNVERIFIED, "frozen turbulence (Taylor): the time-defined gust is uniform over the whole wing and x = U_ref t; not tested numerically.")
end

# =============================================================================
# (c) white-noise coloring
# =============================================================================
# effective (Hermitian) filter of real(ifft(fft(x) .* H)) for x real:  (H_k + H_{N-k})/2
effH(H) = (N = length(H); [(H[k+1] + H[mod(N - k, N)+1]) / 2 for k in 0:N-1])
# one-sided PSD (per unit of the sampled coordinate frequency, angular) predicted for y = real(ifft(fft(x).*H)) rescaled to RMS σ
function predicted_S(H, δx, σ)
    N = length(H); He = effH(H); s2 = σ^2 * N / sum(abs2, He)
    k = 1:N÷2
    return 2π .* k ./ (N * δx), (δx / π) .* abs2.(He[k.+1]) .* s2          # Ω_k (true bin), S_k [per rad]
end
# Monte Carlo periodogram of rescaled realizations
function mc_S(gen, N, δx, M)
    S = zeros(N ÷ 2)
    for m in 1:M
        y = gen(m); S .+= (δx / (π * N)) .* abs2.(fft(y))[2:N÷2+1]
    end
    return S ./ M
end
logbin(x, n=25) = [mean_(x[i:min(i + n - 1, end)]) for i in 1:n:length(x)]
mean_(v) = sum(v) / length(v)

function sec_whitenoise()
    s = "2.4-ruidoBranco"
    C(l, d, a, b; kw...) = check(s, l, d, a, b; kw...)
    B(l, d, c; kw...) = checkb(s, l, d, c; kw...)
    B("[texto] eq:uvk_whitenoise", "text now states DFT frequencies Ω_k=2πk/ℓ mirrored (Hermitian spectrum) and that re-scaling is required [\\fix{} after the audit]", printed("\\mathbb{R}\\lp\\mathcal{F}^{-1}") && printed("Hermitian") && printed("must be re-scaled"))
    # ---------------- space version (create_Continuous1DSpaceGust): N nodes over ℓ
    N, ℓ, L, σ = 1001, 20000.0, 762.0, 1.0
    δx = ℓ / (N - 1)
    Ωhat(k) = 2π * min(k, N - k) / (N * δx)                                  # Hermitian frequencies of the DFT bins
    Hherm = [sqrt(Ψw(Ωhat(k))) for k in 0:N-1]
    Htext = sqrt.(Ψw.(collect(LinRange(0, 2π * N / ℓ, N))))                  # literal reading of "0 ≤ Ω ≤ 2πN/ℓ"
    Hcode = sqrt.(Ψw.((2π / ℓ) .* min.(0:N-1, N .- (0:N-1))))                  # src/Gust.jl (corrected): LΩ = 2πL/ℓ * min(k, N-k)
    Ωk, Sh = predicted_S(Hherm, δx, σ)
    _, St = predicted_S(Htext, δx, σ)
    _, Sc = predicted_S(Hcode, δx, σ)
    Sdef = Ψw.(Ωk)
    band = findall(k -> 0.5 <= L * Ωk[k] <= 100, eachindex(Ωk))              # LΩ ∈ [0.5, 100]
    ks = [round(Int, x * ℓ / (2π * L)) for x in (1, 10, 100)]                # LΩ = 1, 10, 100
    spread(S) = (r = S[band] ./ Sdef[band]; maximum(r) / minimum(r))
    diag(@sprintf("space N=%d ℓ=%g L=%g: S_gen/Ψ_def at LΩ=1,10,100:  Hermitian %s | text-literal %s | code %s", N, ℓ, L,
        string(round.((Sh./Sdef)[ks], digits=3)), string(round.((St./Sdef)[ks], digits=3)), string(round.((Sc./Sdef)[ks], digits=3))))
    diag(@sprintf("spread max/min of S_gen/Ψ_def over LΩ∈[0.5,100]: Hermitian %.3f | text-literal %.3f | code %.3f", spread(Sh), spread(St), spread(Sc)))
    checkv(s, "forma (Hermitiano)", "Hermitian filter √Ψ(|Ω_k|): generated PSD / definition constant (max/min over LΩ∈[0.5,100])", spread(Sh), 1.05)
    checkv(s, "forma (texto literal)", "text read literally (Ψ at LinRange(0,2πN/ℓ), no mirroring): PSD distorted (max/min of generated/definition)", spread(St), 1.25; kind=:note)
    checkv(s, "forma (código)", "[CÓDIGO corrigido] Ω_k = 2π min(k,N-k)/ℓ (DFT spacing, mirrored; it was LinRange(0,(N-1)/2)·2π/ℓ): generated/definition constant (max/min)", spread(Sc), 1.05)
    # unscaled variance
    Hs = Hherm
    varraw = sum(abs2, Hs) / N
    C("variância sem escala", "unscaled Hermitian coloring: Var = (1/N)Σ|H_k|² = σ² Δx/π (1−tail), Δx=ℓ/N → rescaling is mandatory, not optional", varraw, σ^2 * δx / π; tol=0.05)
    diag(@sprintf("unscaled variance / σ² = %.3f ; σ²Δx/π = %.3f", varraw / σ^2, δx / π))
    # R(·) is redundant with a Hermitian filter, necessary (and harmful) otherwise
    rng = MersenneTwister(5); x = randn(rng, N)
    imH = maximum(abs.(imag(ifft(fft(x) .* Hherm)))) / maximum(abs.(real(ifft(fft(x) .* Hherm))))
    imC = maximum(abs.(imag(ifft(fft(x) .* Hcode)))) / maximum(abs.(real(ifft(fft(x) .* Hcode))))
    checkv(s, "R(·) c/ filtro Hermitiano", "with a Hermitian filter ifft(...) is already real (imag/real ~ 1e-16): R(·) redundant", imH, 1e-10)
    checkv(s, "R(·) c/ filtro do código", "[CÓDIGO corrigido] with the mirrored filter imag/real ~ 1e-16 (real(·) is redundant)", imC, 1e-10)
    # validate the closed-form effective filter by Monte Carlo (replica of the code recipe, rescale per realization like the code)
    Mmc = 300
    gen_code(m) = (r = MersenneTwister(1000 + m); y = real(ifft(fft(randn(r, N)) .* Hcode)); y .* (σ / sqrt(sum(abs2, y) / N)))
    gen_herm(m) = (r = MersenneTwister(5000 + m); y = real(ifft(fft(randn(r, N)) .* Hherm)); y .* (σ / sqrt(sum(abs2, y) / N)))
    Smc_c = mc_S(gen_code, N, δx, Mmc); Smc_h = mc_S(gen_herm, N, δx, Mmc)
    C("MC vs forma fechada (código)", "Monte Carlo PSD (300 realizations, 25-bin averages) = (δx/π)|(H_k+H_{N-k})/2|² σ²N/Σ|·|² for the code filter", logbin(Smc_c) ./ logbin(Sc), ones(length(logbin(Sc))); tol=0.08)
    C("MC vs forma fechada (Hermit.)", "same for the Hermitian filter", logbin(Smc_h) ./ logbin(Sh), ones(length(logbin(Sh))); tol=0.08)
    if HAVE_AB
        Wreal(m) = (g = AeroBeams.create_Continuous1DSpaceGust(spectrum=:vK, gustLength=ℓ, N=N, L=L, σ=σ, seed=100 + m); X = collect(LinRange(0, ℓ, N)); [g.W([0.0, x, 0.0]) for x in X])
        Smc_r = mc_S(Wreal, N, δx, Mmc)
        C("código real vs replica", "create_Continuous1DSpaceGust: MC PSD of W = closed form of the replica (Hcode)", logbin(Smc_r) ./ logbin(Sc), ones(length(logbin(Sc))); tol=0.08)
        rr = logbin(Smc_r) ./ logbin(Sdef)
        ib(x) = findfirst(i -> L * Ωk[1+25(i-1)] > x, eachindex(rr))
        diag(@sprintf("REAL create_Continuous1DSpaceGust: PSD(W)/Ψw_def (25-bin averages): first=%.3f, LΩ>3: %.3f, LΩ>10: %.3f, LΩ>50: %.3f, last=%.3f", rr[1], rr[ib(3)], rr[ib(10)], rr[ib(50)], rr[end]))
    end

    # ---------------- time version (stochastic_gust_velocity_from_white_noise, vK): Δt = π/ωmax, N samples
    U, Tgen, ωmax = 25.0, 1000.0, 200π
    Δt = π / ωmax; Nt = length(0:Δt:Tgen-Δt); fmax = 1 / (2Δt)
    a = aVK; TL = a * L / U; Cc = σ * sqrt(TL / (a * π))
    fc = fftfreq(Nt, 1 / Δt)                                                    # frequency vector of the code (signed DFT frequencies, corrected)
    Gvw = @. Cc * ((1 + 2 * sqrt(2 / 3) * im * 2π * fc * TL) / (1 + im * 2π * fc * TL)^(11 / 6))
    Φw_code = abs2.(Gvw)
    Ht = sqrt.(Φw_code)
    He = effH(Ht); s2 = σ^2 * Nt / sum(abs2, He)
    kk = 1:Nt÷2; ftrue = kk ./ (Nt * Δt)
    Sgen = 2Δt .* abs2.(He[kk.+1]) .* s2                                        # per Hz
    Sdeft = (2π / U) .* Ψw.(2π .* ftrue ./ U; L, σ)                              # S(f) = Ψ(Ω) dΩ/df, Ω = 2πf/U
    Hh = [sqrt((2π / U) * Ψw(2π * (min(k, Nt - k) / (Nt * Δt)) / U)) for k in 0:Nt-1]   # corrected filter: Hermitian, argument 2πf/U
    Heh = effH(Hh); s2h = σ^2 * Nt / sum(abs2, Heh); Sgenh = 2Δt .* abs2.(Heh[kk.+1]) .* s2h
    bandt = findall(f -> 0.002 <= f <= 20, ftrue)
    spreadt(S) = (r = S[bandt] ./ Sdeft[bandt]; maximum(r) / minimum(r))
    fidx = [findmin(abs.(ftrue .- f))[2] for f in (0.003, 0.03, 0.3, 3.0)]
    diag(@sprintf("time version U=%g L=%g: S_gen/S_def(f=0.003,0.03,0.3,3 Hz): code %s ; Hermitian-corrected %s", U, L, string(round.((Sgen./Sdeft)[fidx], digits=3)), string(round.((Sgenh./Sdeft)[fidx], digits=3))))
    diag(@sprintf("time version: spread max/min over f∈[0.002,20] Hz: code %.3g | corrected %.3g ; corner f_def=%.4f Hz (1/(2πT_L)) vs filter corner 1/T_L=%.4f Hz → ratio 2π", spreadt(Sgen), spreadt(Sgenh), 1 / (2π * TL), 1 / TL))
    checkv(s, "tempo: forma (corrigido)", "corrected time recipe (Hermitian filter, argument ω=2πf in the von Kármán form): generated/definition constant over f∈[0.002,20] Hz", spreadt(Sgenh), 1.05)
    checkv(s, "tempo: forma (código)", "[CÓDIGO corrigido] stochastic_gust_velocity_from_white_noise now uses signed DFT frequencies and 2πf (it used f in Hz on LinRange(0,fmax,N)): generated/definition constant (max/min)", spreadt(Sgen), 1.05)
    if HAVE_AB
        tvec = collect(0:Δt:Tgen-Δt)
        Random.seed!(3)
        _, _, _, f_r, Φu_r, Φv_r, Φw_r = AeroBeams.stochastic_gust_velocity_from_white_noise(:vK, tvec, U, L, σ)
        C("código real: Φw", "replica of Φw=|G_vw|² (non-negative frequencies) equals the function output", Φw_r, Φw_code[1:Nt÷2+1]; tol=1e-12)
        Mt = 40
        function genreal(m)
            Random.seed!(900 + m)
            _, _, W, _, _, _, _ = AeroBeams.stochastic_gust_velocity_from_white_noise(:vK, tvec, U, L, σ); W
        end
        Smc_t = zeros(Nt ÷ 2)
        for m in 1:Mt
            Smc_t .+= 2Δt / Nt .* abs2.(fft(genreal(m)))[2:Nt÷2+1]
        end
        Smc_t ./= Mt
        C("MC tempo vs forma fechada", "REAL stochastic_gust_velocity_from_white_noise: MC PSD (40 realizations, 200-bin averages) = closed form (code filter)", logbin(Smc_t, 200) ./ logbin(Sgen, 200), ones(length(logbin(Sgen, 200))); tol=0.15)  # lowest bin: the per-realization RMS rescale makes it noisy (low frequencies dominate the variance)
        rr = logbin(Smc_t, 200) ./ logbin(Sdeft, 200)
        diag(@sprintf("REAL time function: PSD(W)/S_def (200-bin averages): min=%.3g max=%.3g", minimum(rr), maximum(rr)))
    end
end

# =============================================================================
# (e) coherence with the use of the gust in the aerodynamics (rotation chain)
# =============================================================================
function sec_chain()
    s = "2.4-cadeia"
    rng = MersenneTwister(11)
    rot() = (M = Matrix(LinearAlgebra.qr(randn(rng, 3, 3)).Q); det(M) < 0 ? M * Diagonal([1.0, 1.0, -1.0]) : M)
    R, Rz, Rwp, RA = rot(), rot(), rot(), rot(); Ug = randn(rng, 3)
    Rw = Rz * Rwp * Rz'
    text = (RA * R * Rw * Rz)' * Ug                  # eq:Ugs last line: (R_A R R_w R_z)^T U_g^o
    code = (R * (Rz * Rwp))' * RA' * Ug              # Aerodynamics.jl: (R*RwR0)' * R_AT * UGustInertial, RwR0 = R0*Rw⁺
    check(s, "eq:Ugs vs código", "(R_A R R_w R_z)ᵀ U_g^o = (R R_z Rw⁺)ᵀ R_Aᵀ U_g^o (code), with R_w = R_z Rw⁺ R_zᵀ", text, code; tol=1e-12)
    # inertial axes used by the code: î₂ = longitudinal (convective velocity, space coordinate x[2]), î₃ vertical, î₁ lateral (DARPA width x[1])
    checkb(s, "eixos (código)", "[CÓDIGO corrigido] flight direction is î₂: the 2nd inertial component carries the longitudinal PSD Φu (and the convective velocity), the 1st the lateral Φv (src/Gust.jl)", occursin("sum( sqrt2dΩΦu .* cos.(ω*t.+ψv) ) .+ convectiveVelocity", CODE_GUST) && occursin("sum( sqrt2dΩΦv .* cos.(ω*t.+ψu) )", CODE_GUST))
end

# =============================================================================
println("="^100)
println("Verification of section 2.4 (Gust field profiles) of theory/main.tex   -   AeroBeams loaded: $HAVE_AB")
println("="^100)
sec_psd()
sec_discrete()
sec_sinusoids()
sec_whitenoise()
sec_chain()

println()
let lastsec = ""
    for r in RESULTS
        mark = r.kind == :pass ? (r.ok ? "ok" : "FALHOU") :
               r.kind == :typo ? (r.ok ? "ERRO NO TEXTO CONFIRMADO" : "!! esperava-se discrepancia e nao houve") :
                                 (r.ok ? "NOTA (afirmacao imprecisa/nao geral/inconsistente com o codigo, ver relatorio)" : "!! esperava-se diferenca e nao houve")
        if r.sec != lastsec
            println("\n------------------------- ", r.sec, " -------------------------")
            lastsec = r.sec
        end
        @printf("%-13s %-30s %-150s %10.3e  %s\n", r.sec, first(r.label, 30), first(r.desc, 150), r.err, mark)
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
