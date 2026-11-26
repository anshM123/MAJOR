# TAIL_CERT: the large-impact-parameter tail L ≥ L_tail of full AdS positivity at Δ_g = 1000 (running document)

*26 Nov 2026. Research agent. Unreviewed, not committed. No existing file edited; nothing under `prl/` touched.*
*Code: `ads4/code/tailcert/` (outputs in `cache/`, git-ignored). Compute: ≤ 2 processes, `OMP_NUM_THREADS=1`.*

Grades: **[certified]** Arb balls / exact rationals with every remainder bounded; **[P]** proved on paper here; **[DV]** sampled
numerics; **[author]** claim not reproduced.

---

## R. Item 3 and the author's Lemma R (coordinator request, done first)

### R.1 Signed low-frequency residual, measured [DV] (`tc_lownu_dv.py`, `cache/tc_lownu_dv.txt`)
G_r^low(β) := (2/Λ)∫₀^{ν₁} r(ν)P_{−½+iν}(cosh β)dν, r = E − ψ₂(ν/Λ)/k₀, E = 2Re[χR₋](ν−i), repaired weight (θ = 1 to 10⁻²⁰⁰ here).
This is the low-ν part of (1−c)³G^r(c) (FULLPOS3 §4.1(c)), the term whose *absolute* bound grows like β² in FULLPOS5 §3.
P from the two-branch formula (checked against mpmath `legenp` to 14 digits at β = 30).

| β | F_cert | signed G_r^low/F_cert, ν₁ = 0.31 | ν₁ = 1.0 | absolute bound (2/Λ)∫₀^{ν₁}\|r\|·P_{−½}/F_cert, ν₁ = 0.31 |
|---|---|---|---|---|
| 30 | 1.176·10⁻¹² (LRINF cell [30, 30.5]) | **+0.0021** | +0.0021 | 0.011 |
| 64 | 3.61·10⁻²⁰ (fcert_tail) | **+0.0031** | +0.0031 | 0.032 |
| 100 | 3.52·10⁻²⁸ | **+0.0036** | (see cache) | 0.077 |
| 200 | 3.40·10⁻⁵⁰ | **+0.0040** | (see cache) | 0.306 |

The signed piece is **positive** and saturates at a few ×10⁻³ of the floor, while the absolute bound grows like β² (0.31 at β = 200),
exactly as FULLPOS5 §3 found. So the β² growth is an artefact of taking |r|.

### R.2 Lemma R (ANGULAR_ROUTE §81) reproduced in Arb [certified, given EPSW's certified inputs; two θ-corrections [P]]
(`tc_lemmaR.py`, `cache/tc_lemmaR.txt`)
* **Inputs.** From `epswcert/cache/ew_cert_1000_re_rep.json` (EPSW_CERT, certified): r(0) = 9.6925·10⁻¹⁴, ε_W ≥ ∫\|r\| = 2.68136·10⁻⁶,
  B_EPSW ≥ \|r′(0⁺)\| + ∫\|r″\| = 1.89971·10⁻³. For r̃ = r_θ − r(0)sech(πν):
  * E_in = ε_W + r(0)/2 + δE_θ ≤ **2.68236·10⁻⁶** (author: < 2.683·10⁻⁶ ✓);
  * B_in = B_EPSW + π²r(0) + δB_θ ≤ **1.89971·10⁻³** (author: < 1.901·10⁻³ ✓).
  * Here ∫sech(πν) = ½, ∫\|(sech πν)″\| ≤ π² (sech″ = sech − 2sech³), and δE_θ, δB_θ ≤ 10⁻⁹ is **[P]**: on [0, Λ+30], 1 − θ ≤ erfc(7)
    and \|θ′\|, \|θ″\| ≤ e⁻⁴⁸, against ∫\|E\|, sup\|E′\| ≤ B_EPSW + sup\|ψ′\|; beyond Λ+30 all RBF centres are ≤ Λ+30, so θ′E′ is ≤ e^{−48}.
    These are e^{−40}-size, far below the 10⁻⁹ allowance, but the sup bounds were not ball-evaluated.
* **Mehler step, re-derived [P].** W(t) = ∫r̃ cos νt, \|W\| ≤ min(E_in, B_in/t²) (two integrations by parts), ∫₀^∞W = (π/2)r̃(0) = 0.
  With k(t) = (cosh β − cosh t)^{−½} on [0, β) and K₀ = k(0):
  I_{r̃} = (√2/π)∫₀^∞W(t)[k(t)1_{t<β} − K₀]dt, and
  **\|I_{r̃}\| ≤ (√2/π)K₀[B/β + 8B log2/β² + E·ε₁(β)]**, ε₁ = (sinh(β/2) − β/2)/(2(cosh β − cosh(β/2))) ≈ e^{−β/2}/2.
  * [0, β/2]: √(1+x) ≤ 1 + x/2 gives k − K₀ ≤ K₀(cosh t − 1)/(2(cosh β − cosh t)).
  * [β/2, β]: \|W\| ≤ 4B/β², and ∫₀^β(k/K₀ − 1)dt ≤ ∫₀¹((1−x)^{−½} − 1)dx/x = 2 log 2 exactly (substitute x = (cosh t−1)/(cosh β−1);
    dt/dx ≤ 1/x). This makes the author's 8B log2/β² term exact, for every β.
  * [β, ∞): K₀∫_β^∞B/t² = K₀B/β.
* **Result [certified arithmetic].**
  * β ∈ [30, 64], against LRINF's certified cells: (2/Λ)\|I_{r̃}\| ≤ **0.0249·F_cert** (worst cell [30, 30.5]); author: < 0.027F_*.
  * β ≥ 64, against `fcert_tail`: the ratio 2Λ[B + 8B log2/β + Eβε₁]/((1−2.1·10⁻⁷)A₁ℓ) decreases in β, so ≤ its β = 64 value
    **0.01443** (author: < 0.015); limit 0.01328. A₁ = 2ψ₂′(0) − ψ₃′(0) = 327.149, ℓ = 64L(64) = 0.87401.
  * The extracted piece r(0)sech(πν) contributes +(2/Λ)r(0)/(2cosh(β/2)) > 0 (and, at finite mass, r(0)/Λ·w³𝒱_b((1−z)^{−5/2}) ≥ 0,
    all coefficients nonnegative — consistent with the author).
* **Consequence.** In the FULLPOS3/5 h-expansion route the m = 0 realised-residual term G^r is ≥ −0.025F_cert for every β ≥ 30, uniformly:
  the β² obstruction of FULLPOS5 §3.2 is removed **for the m = 0 term**. What Lemma R does not cover (author agrees): the finite-mass
  (ladder) transfer of r̃ and the high-frequency residual pieces in the local terms h^mP_m(ϑ)G^r, which in FULLPOS5 grow like β.

### R.3 How Lemma R sits in the contour (tail) route
In the contour route (§§1–3 below) the realised f is never split into ψ + r: the main term
𝓜(L) = √w[f(0)/Λ + (2/Λ)∫₀^{3/4}W e^{−Ly}dy] uses the realised f exactly (W ≥ 0 certified by CONTOUR, f(0) = r(0) > 0 by EPSW), and the
finite-mass transfer is the factor (1 + D) with D the ladder deviation of a pure power, applied to the whole f. So the contour route
needs neither Lemma R nor any bound on r: the "higher-order finite-mass transfer of r̃" that Lemma R leaves open is, in this route,
the D-bound of §2, which acts on f as a whole.

---
(main sections follow; in progress)
