# Aerodynamic center, lift level and tangent lift slope of the Pazy wing

> Update 2026-11-21: the `:VLM` tip loss of the package is now the tangent lift slope (`dev/Pazy/General/generatePazyTipLossFunctionData.jl`), the best single choice for the three wings with the updated `NACA0018Timmer` (MAE 1.28 vs. 1.81 m/s with the secant). The secant table is kept here as `tipLossSecantVLM.csv` (`TLMODE=sec`); `TLMODE=pkg` (default) uses the package. The tables below were computed before this change, with the secant `:VLM`.

Study of 2026-11-21 on what the calibrated aerodynamic center of `NACA0018Pazy` (`x_ac` as a function of the Mach number) is absorbing in the flutter boundaries of the Pazy wings S0, S10 and S20 (beam rod models, `Indicial()` with the default `finiteWingLag=true`, `"VLM"` tip loss), and which of its effects can be replaced by data-based or physical corrections. Nothing in `src` was changed: the two prototypes (spanwise shift of the aerodynamic center) are in `aeroCenterShift.patch`, to be applied to a copy of the package.

Errors are the AeroBeams flutter speed minus the experimental one, in m/s, over 26 points: S0 hump, lower edge (test onset with the airspeed increasing, θ = 3, 5, 7°) and upper edge (test onset with the airspeed decreasing, θ = 3, 5°), and the first flutter onset of S10 LE/TE and S20 LE/TE vs. root pitch angle. MAE is the mean absolute error over the 26 points.

## Files

- `runVariant.jl`: runs one of the beam rod flutter scripts (or the S0 static one) with a variant of the aerodynamic model and prints the errors relative to the experiments.
- `variants.jl`: the variants, selected with environment variables (airfoil, finite-wing lag, tip loss scale, Nastran-level tip loss, secant/tangent lift slope, shift of the aerodynamic center).
- `generateStudyTables.jl`: writes `aeroCenterShiftVLM.csv` (Δx = x_cp − 1/4 = −cm_c4/cn at the VLM equilibria) and `tipLossTangentVLM.csv` (tangent lift slope, from the fixed-shape perturbations of the VLM) from `test/referenceData/PazyVLM`, with the anchors and rules of `generateSweptPazyTipLossFunctionDataDeformed.jl`.
- `aeroCenterShift.patch`: prototype of the spanwise shift of the aerodynamic center (`git apply` or `patch -p1` on a copy of the package, then `--project=<copy>`; the variants with `XS=1` need it).
- `neededAerodynamicCenter.png`: uniform `x_ac` that makes each flutter point match the experiment, vs. the `x_ac(U)` of `NACA0018Pazy` and of the data of Timmer.

Example: `AF=timmerfull TLMODE=swap julia -t 1 --project=. dev/Pazy/General/aeroCenterStudy/runVariant.jl s10 timmerSwap` (about one minute per script; `s0`, `s10`, `s20`, `s0static`).

## Findings

**1. The lift level is not the cause.** With the VLM table regenerated with `core_radius = 1e-6`, the integral of the lift slope along the span of the VLM relative to the Nastran model of Technion (θ = 5°) is 0.995 to 0.964 for Λ = 0° and 1.007 to 0.984 for Λ = 20°, from 10 to 60 m/s (the largest local difference is at the root of S0 at 60 m/s: 5.96 vs 6.78). Multiplying the tip loss by the Nastran-to-VLM ratio moves the S0 hump onset by −0.1 to −0.3 m/s. A uniform increase of the lift also lowers the upper edge of the hump (×1.10: −2.7 m/s at θ = 3°) and narrows the region, whereas lowering `x_ac` widens it and leaves the upper edge almost unchanged: the two are not equivalent.

**2. The cm of the VLM is a shift of the center of pressure, not a tip loss.** The cm about the quarter-chord of a 2D flat plate is zero, so there is nothing to normalize it with; the well-posed quantity is Δx = x_cp − 1/4 = −cm_c4/cn (`spanwise_cm.csv` and `spanwise_cn_alpha.csv`). For S0 (θ = 3 to 7°, U = 30 to 70 m/s) x_cp is 0.245–0.250 up to two thirds of the span, 0.240 at 80%, 0.21–0.22 at 93%, and 0.21 to 0.16 (θ = 3°) at the tip station, down to negative values at θ = 7° and 70 m/s, where cn vanishes; the load-weighted average is 0.237–0.245. For Λ = 20° the root is aft (0.26–0.265). The secant (−cm/cn) and the incremental (fixed-shape slopes) definitions agree except at the tip station. Adding Δx(s) to the airfoil `x_ac` in the moment of the circulatory load (`aero_load_at_ζ`, prototype in the patch) lowers the S0 hump onset by 0.5–0.6 m/s, i.e. as a uniform shift of about −0.008 for S0 and −0.012 for the LE configurations of S10 and S20, and lowers the TE speeds by 1.4 (S10) and 1.2 (S20) m/s.

**3. Uniform `x_ac` (S0, S10 and S20 with the VLM tip loss).** Lowering it widens the S0 hump (lower edge −1.8 m/s on average from 0.25 to 0.22; upper edge by ±0.5 m/s), lowers the LE flutter speeds (−2.8 m/s per −0.03 at θ = 0 in S10) and the hard mode, and also lowers the TE speeds, which are already below the experiments. The `x_ac` that matches each experimental point (`neededAerodynamicCenter.png`) is 0.19–0.26 (mostly 0.20–0.23) for the S0 hump and the LE configurations (38–56 m/s), but 0.27–0.29 (S10) and 0.34–0.36 (S20) for the TE ones (73–78 m/s): a uniform `x_ac` cannot fit both, and the `x_ac(U)` of `NACA0018Pazy` (0.21 at 0 m/s to 0.25 at 100 m/s) interpolates between the two groups. The trend of the data of Timmer (`NACA0018Timmer`) is the opposite one (x_ac from 0.258 to 0.214 as the Reynolds number goes from 0.15e6 to 1e6).

**4. What `x_ac` ≈ 0.22 is compensating.** Flutter of the undeformed wings with the unsteady VLM in the basis of the modes of AeroBeams (`../frequencyDomainUVLMFlutter`, same beam rod structures; onset in m/s at θ = 0; experiments: S10 LE 56.0, S10 TE 78.1, S20 LE 52.3, S20 TE 76.7 at 3°):

| | UVLM | strip theory, x_ac = 0.25 | experiment |
|---|---|---|---|
| S0 | 69.9 | 71.4 | |
| S10 LE | 54.4 | 60.1 | 56.0 |
| S10 TE | 76.3 | 76.4 | 78.1 |
| S20 LE | 50.5 | 54.6 | 52.3 |
| S20 TE | 75.1 | 72.5 | 76.7 |

The potential-flow 3D model is 1.6–1.8 m/s below the four experiments, while the strip theory is above it for LE and below it for TE: the same pattern that the calibration of `x_ac(U)` follows, without any viscous effect. In the generalized aerodynamic forces of the S0 modes (strip theory with the tip loss and the finite-wing lag vs. UVLM, k = 0.1 to 0.3) the torsion–torsion term agrees within 1% (real part, which is the one that `x_ac` scales: `x_ac` = 0.22 would make it 14% larger), the aerodynamic damping of the bending modes is 11–15% higher in the strip theory and the lift on the second bending mode due to the twist 7–8% lower: nonlocal 3D effects of the mode shapes that no local tip loss function reproduces. The least-squares fit of those forces with free lag coefficients, a scale of the tip loss and `x_ac` gave `x_ac` = 0.233. So the calibrated `x_ac(U)` mostly compensates for the 3D unsteady effects of the modes, and cannot be both right for the LE (high k) and the TE (low k) modes.

**5. Secant vs. tangent lift slope.** `"VLM"` is the secant, cn/α_local, of the static aeroelastic solution. For the swept wings at high airspeed it differs a lot from the tangent (Δcn/Δα_local, fixed shape), which is the slope seen by the perturbations of the eigenproblem and is close to the undeformed distribution: Λ = 20°, θ = 3°, U = 70 m/s, root to tip, 0.60 0.72 0.83 0.91 0.93 0.90 0.81 0.58 0.48 (secant) vs. 0.75 0.77 0.79 0.80 0.79 0.77 0.71 0.53 0.44 (tangent; undeformed wing: 0.75 0.77 0.79 0.80 0.79 0.77 0.71 0.53 0.45). Using the secant for the steady solution and the tangent for the eigenproblem (`TLMODE=swap`: re-linearization at the same equilibrium; equivalent to a spanwise zero-lift angle offset for the static preload) raises the TE flutter speeds by 2.6 m/s (S20) and 0.7 m/s (S10) and the LE ones by at most 0.5 m/s, without changing the S0 results; keeping the undeformed tip loss (`TLTH0=1`) or the tangent for statics too (`TLMODE=tan`) gives similar results. The jump of the secant between θ = 0 (anchor) and θ ≥ 1° is what lowers the S20 TE flutter speed from 72.7 m/s at θ = 0 to 70.0 m/s at θ = 1° and 69.7 m/s at θ = 3° (72.6 and 72.3 m/s with the undeformed tip loss).

**6. Calibration-free alternative.** The airfoil of the data of Timmer (`NACA0018Timmer`, Reynolds-dependent cnα and `x_ac`) with the tangent lift slope in the eigenproblem has a lower MAE than the calibrated `NACA0018Pazy`, with the opposite error distribution (better TE, worse S0 and LE):

| model | S0 lower edge θ=3/5/7 | S0 upper edge θ=3/5 | S10 LE | S10 TE | S20 LE | S20 TE | MAE |
|---|---|---|---|---|---|---|---|
| x_ac = 0.25 (flat plate), secant | +1.9/+1.5/+2.7 | +0.9/−0.8 | +1.7 | −2.7 | +1.7 | −6.6 | 2.64 |
| NACA0018Pazy, secant | −0.3/−0.2/+0.9 | 0.0/−1.2 | +0.2 | −3.7 | +0.2 | −7.6 | 2.32 |
| NACA0018Timmer, secant | +1.8/+1.2/+1.9 | +1.6/−0.4 | +1.8 | −1.3 | +1.5 | −5.3 | 2.13 |
| NACA0018Timmer + tangent | +1.7/+1.1/+1.9 | +1.6/−0.4 | +2.0 | −0.4 | +2.0 | −2.3 | 1.65 |
| NACA0018Timmer + tangent + Δx(s) | +1.2/+0.7/+1.5 | +1.0/−0.7 | +1.2 | −1.9 | +1.3 | −3.9 | 1.82 |
| NACA0018Timmer + tangent, x_ac − 0.01 | +1.0/+0.6/+1.4 | +1.2/−0.5 | +1.4 | −1.3 | +1.5 | −3.1 | 1.63 |
| NACA0018Pazy + tangent | −0.3/−0.2/+0.9 | 0.0/−1.1 | +0.4 | −2.9 | +0.7 | −5.0 | 1.85 |
| x_ac = 0.25, Nastran-level tip loss | +1.8/+1.4/+2.4 | +1.1/−0.5 | +1.3 | −2.3 | +1.1 | −6.5 | 2.32 |
| x_ac = 0.25, tip loss ×1.10 | −0.2/−0.4/+0.9 | −1.8/−3.0 | +0.1 | −4.8 | +0.2 | −8.6 | 3.04 |

For the lower edge of the hump, the independent 3D potential-flow solution of Sharpy (`flutterBoundary_UVsPitch_Lambda0_Sharpy.txt`) is +1.3/+0.8/+1.3 m/s above the experiments, i.e. the residual of the data-based models is of the order of that of the other potential-flow methods. Static S0 (tip OOP rms in % of the semispan at θ = 3/5/7°; tip twist error at θ = 5/7°): x_ac = 0.25: 1.8/1.7/2.6, −10/−30%; NACA0018Pazy: 1.9/1.7/2.0, +5/−17%; NACA0018Timmer: 2.1/2.1/2.9, −6/−27%; NACA0018Timmer + Δx(s): 2.0/2.1/2.7, −1/−23%.

**7. Finite-wing lag.** All of the above has `finiteWingLag=true` (the default now). Without it the S0 torsion flutter at θ ≈ 0 disappears below 100 m/s (UVLM, Sharpy and ZAERO: 64.5–70 m/s), the TE flutter frequencies are 2–4 Hz low (14–15 vs. 16–19 Hz), and the TE speeds are 1.4 (S20) to 3.2 (S10) m/s higher; e.g., `NACA0018Pazy` without it has MAE 2.18 (+1.9/+0.9/+1.7 for the S0 lower edge).

## Caveats

- The prototypes (`XS=1`) and the swap of the lift slope are not in `src`; `solve_swap!` re-linearizes the eigenproblem at the steady solution of the secant model, which is exact for the Jacobian but leaves the residual of the tangent model nonzero. A model with a zero-lift angle offset (spanwise, from the VLM) would be the equivalent for any analysis.
- The tangent table replaces 89 entries (S0 tip station at θ ≥ 5° and U ≥ 70 m/s, and mid-span at U ≥ 100 m/s) by those of the undeformed wing, and Δx is clamped to [−0.15, 0.05]; none of them is in the airspeed range of the hump.
- The `x_ac` of the data of Timmer is the secant between 0 and 7° of the lift and moment coefficients; the flutter perturbations are small, so the slopes at 2–4° could differ.
- The experimental flutter onsets have an uncertainty of the order of 1–2 m/s (hysteresis of the S0 hump, discrete tunnel speeds); differences below that should not be over-interpreted.
