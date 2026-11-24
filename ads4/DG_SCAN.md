# Δ_g dependence of Results D and Q: bound values as functions of Λ, positivity margins vs Δ_g, and what can be stated over a range

*22 Nov 2026, Δ_g-scan agent. Unrefereed, not committed. No existing document or code file was edited. New code:
`ads4/code/dgscan/bounds.py`, `tol.py`, `analyze.py`. Caches: `ads4/code/dgscan/cache/` (git-ignored through its own `.gitignore`).
Read-only imports: `omegau/nf_common` (NF02), `resultq/functionals_q3.json`, `spinfunc/ee1000_common` (the `hyb` near-gap states of
`q3_analysis`/`w2_chain`). Read-only parsed caches: `spinfunc/cache/y2_chain100_jhi*.txt`, `w2_chain_final_xc0.txt`. Nothing in
`spinfunc/lj_*` was imported or run (the `T_B` formula of `lj_contour.run` is re-typed in `tol.py`). Compute: ≈ 15 CPU-min,
≤ 2 processes, `ulimit -v 3 GB`, `OMP_NUM_THREADS=1`.*

Grades: **[P]** proved; **[D]** derived; **[D-sketch]**; **[DV]** derived and checked numerically; **[NE]** numerical evidence;
**[M]** model; **[A]** assumption; **[I]** open. `Λ := √((Δ_g−2)(Δ_g−1)) = Δ_g − 3/2 − 1/(8Δ_g) + O(Δ_g⁻²)`.

---

> **Update (GAP_SCALING.md, unreviewed): extending positivity to all Δ_g ≥ 1000 by scaling is PARTIALLY viable. The certified layers carry to every larger gap; the beyond-leading-Regge corrections carry over only under tested scaling hypotheses.**
> * **Weights:** the flat weights ψ_k are exactly Λ-independent in gap units. The k=2 de-shift weight is ψ₂ + Λ⁻²·(explicit symbol) for p ≥ 0.1 [DV]. The leading-Regge action is exactly ∫ψP at every Λ.
> * **Certified layers:**
>   * The curvature is g(J, m/Λ)/Λ² [DV].
>   * The **Lemma C bound is monotone decreasing in Λ [P]**, so "certified flat floor − Lemma C" only improves with the gap; it is positive near the twist edge up to J = 13 / 17 / 23 / 35 / 61 at Δ_g = 1000 / 2000 / 5000 / 2·10⁴ / 10⁵.
>   * **Lemma LR∞ is certified uniformly for 998.5 ≤ Λ ≤ 10⁶** (454 Arb Λ-balls, remainder ratio ≤ 0.9988), so the Theorem L′ floor and J\*(m) hold for every gap in that range [certified].
> * **Binding state at Δ_g = 2000 [NE±]:** own-mass row Δ = 2005 gives λ/N_s = +0.20041, budget 6.3·10⁻⁴, L = 320 (it was 78 at 1000). The truncation scales like Λ⁻² and the window like Λ⁻⁴; the graviton correction profile converges (|Γ| ≤ 41 c₂).
> * **What does not carry over automatically:**
>   * Sign-free budgets at fixed β do not improve (pointwise corrections do not shrink).
>   * EDGE-R's +− transfer factor grows like (Δ_g/J)², so at fixed spin it gets harder.
>   * The graviton edge integral is uncomputed at every gap.
>   * There is a seam strip between the 1000 and new twist edges with no 1000 counterpart (J up to ≈ 344 / 424 / 552 at Δ_g = 2000 / 3000 / 5000).
>   * The structural premises (AG_U, EDGE-R, P1′, W) carry over as premises.
> * **Plan (3–5 days):**
>   1. Write up Lemma C monotonicity and uniform LR∞.
>   2. Own-mass K₂ rows at Δ_g = 2000, 5000 for every error-budgeted spin.
>   3. Direct n1/n3 nodes at 2000 for J = 12, 16, 21, 27, 29 (external machine).
>   4. Parameterise Δ_g in the large-spin code and re-run the P2/P4W-χ costs at the chart-B images.
>
>   Target: error-controlled positivity at Δ_g = 1000, 2000 and 5000, plus a transfer to all Δ_g ≥ 1000 under named, tested scaling hypotheses.

## 0. Summary

1. **The bound values are explicit functions of `Λ` [DV].** With the functionals fixed (NF02 for `Ω^EE`; `Q3-full-delta0.2` for `Ω^Q_up`)
   the light side gives closed forms (exact for polynomial weights; checked against quadrature to `10⁻³⁰` and against every quoted digit
   at 100/300/1000):
   * **cubic:** `√(t₄²+α_T²) ≤ 240√(16L_A(Λ))/Λ⁴ = 1.3194·10⁴ · √(ln(Λ/8.7415) + 6.384/Λ + O(ln Λ/Λ²)) / Λ⁴`
     (relative error of the bracketed form `4.5·10⁻⁵` at 300, `3.6·10⁻⁶` at 1000; without the `6.384/Λ` term `−0.30 %` / `−0.067 %`);
   * **quartic upper end:** `2L_AdS(Λ) = 17.2407·ln(2Λ/3) − 19.2782 + 10.314/Λ + (860.3·ln(2Λ/3) − 2508.6)/Λ² + O(Λ⁻³ ln Λ)`.
     **The drift `B_eff = −19.220 (300) → −19.265 (1000)` is the finite-Λ tail of the exact AdS pole** `1/(p² + 9/(4Λ²))`, converging to
     the flat constant `−19.2782` (= the stored `B` of `Q3-full-delta0.2`). Of the difference `0.0450`: `0.0242` from the `1/Λ` term
     (`p²` coefficient of `ψ₂^Q`), `0.0199` from the `Λ⁻² ln Λ` terms (its large `p³` coefficient, −191), `0.0009` higher order.
   * **lower end:** `−O(Λ⁻²)` (no explicit constant exists; unchanged).
   * **error terms:** `c₁`: `∓19941/(Λ²√(16L_A))` per unit `c₁` (D; `0.216 % / 0.0167 %`); `|ε_{c₁}| = 2.032|c₁t₄|Λ²` (Q) `≤ 487.7√(16L_A)/Λ² |c₁|`
     with D's `t₄` (`0.567 / 0.0586`). H17 band `±212.5/Λ²` on the bound (`±0.24 % / ±0.021 %`). EE `δk(0) = −68.08/Λ²·(1 + O(Λ⁻¹))`
     (bound factor `1 + 34.0/Λ²`).
2. **Positivity vs Δ_g.** Only 100 (fails at model level), 300 and 1000 (hold, [NE + M]) have been computed. Using proxies that reuse
   the computed weights and corrections (§2.2, validated on the 100/300/1000 verdicts):
   * every correction-driven margin **improves monotonically with Δ_g**. The one margin that shrinks is the hybrid baseline, which
     approaches the certified flat value 0.200 monotonically from above (`0.2256 → 0.2043 → 0.2009 → 0.2002` at 100/300/1000/5000)
     and so never drops below it in the model. At fixed `(J, m/Λ)` a beyond-Regge correction costs
     `∝ Λ⁻ᵖ` with `p = 2` (bulk), `→ 4` (small-ν window), empirically `p = 2.83` (median; 1.4–3.2) between 300 and 1000; the
     large-spin window tolerance `T_B` grows as `Λ^{1.20}` (`Λ^{1.0}` at fixed `J/Λ`) [D for the weights, DV; M for the corrections];
   * **threshold `Δ_g*`: not determined. Model bracket `140 ≲ Δ_g* ≲ 300`** [M]. Transferring the adopted 300 corrections at fixed
     `(J, m/Λ)` (spins `J ≤ 32`) gives `Δ_g* ≈ 139 / 172 / 205` for `p = 2 / 2.83 / 4`; requiring the sign-free large-spin window bound
     (refined remainder ≈ 250 c₂, as measured at 300) gives `Δ_g ≳ 275`; with the unrefined window size (450 c₂), `≳ 440`.
     The most dangerous state is always near-gap large spin (`J = 32` in the anchored set, then the `H_LJ` region), never the binding
     `J = 5` state.
3. **Recommended statement (§3).** State the implication for every Δ_g with the explicit Λ-functions; state positivity as computed at
   {300, 1000} [NE + M]; extend to "all `Δ_g ≥ 300`" only as a **conditional corollary** under a named scaling hypothesis `H_mono`
   [D-sketch + M], plus `H_LJ(Δ_g)` whose window part is [I] at 300 and becomes sign-free-closable (refined form) above `Δ_g ≈ 275`, with
   margin ≥ 2 above `≈ 490`. Do not claim anything below 300.

---

## 1. Bound values as explicit functions of Λ (functionals fixed)

### 1.1 Formulas [D]; constants [DV]

* `L(ψ; a) := ∫₀¹ψ(p)dp/(p² + a²)`, `a = 3/(2Λ)` (the exact graviton pole `16πG/(ν² + 9/4)`, H7). For a polynomial
  `ψ = Σ c_n pⁿ` (`c₀ = 0`) it is exact term by term: `I₀ = atan(1/a)/a`, `I₁ = ½ln(1 + a⁻²)`, `I_n = 1/(n−1) − a²I_{n−2}`.
* **Expansion [D]** (`ψ = c₁p + c₂p² + c₃p³ + R`, `R = O(p⁴)`):
  `L = c₁ ln(1/a) + [c₂ + c₃/2 + ∫R/p²] − (π/2)c₂a + a²[c₁/2 + c₂ − c₃ ln(1/a) − ∫R/p⁴] + O(a³)`.
* **Result D** (`16L_A = 16L(ψ₂^{NF}; a)`, `k(0) = 1`):
  `16L_A = A ln(2Λ/3) + B_∞ + b₁/Λ + (b₂L ln(2Λ/3) + b₂)/Λ² + O(Λ⁻³ ln Λ)` with
  `A = 3022.3377`, `B_∞ = −5327.2096` (= `B_AdS + A ln(3/2)`, `B_AdS = −6552.662`), `b₁ = 19293.4`, `b₂L = −13549`, `b₂ = −14705`.
  Equivalently `16L_A = A ln(Λ/Λ₀^D) + b₁/Λ + …`, `Λ₀^D = (3/2)e^{−B_∞/A} = 8.7415`.
* **Result Q** (`2L_AdS = 2L(ψ₂^Q; a)`): `A_Q = 17.240713`, `B_∞ = −19.278176`, `b₁ = 10.3137`, `b₂L = 860.35`, `b₂ = −2508.6`;
  `Λ₀^Q = 4.5889`.

### 1.2 Table (`dgscan/bounds.py` → `cache/bounds.txt`)

| Δ_g | Λ | `16L_A` | **cubic bound** (`k = 1`) | `δk(0)` factor | H17 band | `c₁` (D), per unit | `2L_AdS` | `B_eff` | `ε_{c₁}` (Q, with D's `t₄`) | ABJM `2ζ(3)/2L` |
|---|---|---|---|---|---|---|---|---|---|---|
| 100 | 98.499 | 7508.51 | 2.2094·10⁻⁴ | ×1.0035 | ±2.3 % | ∓2.37 % | 53.107 | −19.0378 | 4.36 \|c₁\| | 4.5 % |
| 150 | 148.50 | 8687.21 | 4.5999·10⁻⁵ | ×1.0015 | ±0.98 % | ∓0.97 % | 60.087 | −19.1361 | 2.06 \|c₁\| | 4.0 % |
| 200 | 198.50 | 9533.02 | 1.5094·10⁻⁵ | ×1.0009 | ±0.54 % | ∓0.52 % | 65.046 | −19.1801 | 1.21 \|c₁\| | 3.7 % |
| 250 | 248.50 | 10193.14 | 6.3542·10⁻⁶ | ×1.0006 | ±0.35 % | ∓0.32 % | 68.895 | −19.2045 | 0.80 \|c₁\| | 3.5 % |
| **300** | 298.50 | **10734.61** | **3.1321·10⁻⁶** | ×1.0004 | ±0.24 % | ∓0.216 % | **72.041** | **−19.2197** | **0.567 \|c₁\|** | 3.3 % |
| 500 | 498.50 | 12259.23 | 4.3031·10⁻⁷ | ×1.0001 | ±0.086 % | ∓0.073 % | 80.855 | −19.2473 | 0.217 \|c₁\| | 3.0 % |
| **1000** | 998.50 | **14339.59** | **2.8913·10⁻⁸** | ×1.00003 | ±0.021 % | ∓0.0167 % | **92.814** | **−19.2647** | **0.0586 \|c₁\|** | 2.6 % |
| 2000 | 1998.5 | 16427.19 | 1.9283·10⁻⁹ | — | ±0.0053 % | ∓0.0039 % | 104.769 | −19.2721 | 0.0157 \|c₁\| | 2.3 % |
| 5000 | 4998.5 | 19192.12 | 5.3262·10⁻¹¹ | — | ±0.00085 % | ∓0.0006 % | 120.571 | −19.2759 | 0.0027 \|c₁\| | 2.0 % |

* **Reproduction [DV].** Every stated digit at 100/300/1000 is reproduced: `16L_A = 7508.51/10734.61/14339.59`, bounds
  `2.2094·10⁻⁴/3.1321·10⁻⁶/2.8913·10⁻⁸`, `2L_AdS = 53.107/72.041/92.814`, `B_eff = −19.0378/−19.2197/−19.2647`, `c₁` `0.22 %/0.017 %`,
  `ε_{c₁}` `0.57/0.059`, `∫ψ₂^{NF} = 10.386`, `∫ψ₂^Q = 1.016`. One inconsistency in the record: RESULT_D §3's table quotes
  "≈ ∓1.7 %" per unit `c₁` at Δ_g = 100; the same linearised formula that gives the 300/1000 entries gives **∓2.37 %** (status-only
  column; nothing is claimed at 100).
* **Above the table, the rows other than 100/300/1000 are values of the implication only.** They are bounds at Δ_g only if per-state
  positivity holds there (§2).

### 1.3 Asymptotic forms and checks

| form | error at 300 | at 1000 | grade |
|---|---|---|---|
| cubic `240√(A ln(Λ/8.7415))/Λ⁴` | −0.30 % | −0.067 % | [DV] |
| cubic `240√(A ln(Λ/8.7415) + 19293/Λ)/Λ⁴` | +4.5·10⁻⁵ | +3.6·10⁻⁶ | [DV] |
| cubic in Δ_g: `240√(A ln(Δ_g/8.7415))(1 + 6/Δ_g)/Δ_g⁴` | −0.25 % | −0.054 % | [DV] |
| quartic `17.2407 ln(Λ/4.5889)` | −0.058 | −0.013 | [DV] |
| quartic `+ 10.314/Λ` | −0.024 | −0.0031 | [DV] |
| quartic `+ (860.3 ln(2Λ/3) − 2508.6)/Λ²` | −0.0009 | −0.00002 | [DV] |

* **Leading behaviour.** `√(t₄²+α_T²) ≲ 1.32·10⁴·√(ln Λ − 2.168)/Λ⁴`, i.e. **`p = 4` with a `√ln` enhancement**, and the next correction
  is relative `+3.19/(Λ ln(Λ/8.74))` (`+0.30 %` at 300). Quartic: `17.24 ln Λ − 26.27 + 10.3/Λ`, i.e. **`O(1/Λ)` with coefficient
  `−(3π/2)c₂ = 10.31`**, then `O(ln Λ/Λ²)` with the large coefficient 860 — so the `1/Λ` form alone is good to 3 decimals only for
  Δ_g ≳ 1000.
* **Error terms as functions of Λ [D].**
  * `c₁` (D): relative shift of the bound per unit `c₁` = `½·16Λ²t₄*∫ψ₂^{NF}/(16L_A) = 1920∫ψ₂/(Λ²√(16L_A)) = 19941/(Λ²√(16L_A))`
    `∼ 363/(Λ²√ln(Λ/8.74))`.
  * `c₁` (Q): `|ε_{c₁}| = 2∫ψ₂^Q·|c₁t₄|Λ² = 2.032|c₁t₄|Λ²`; with D's `t₄`: `487.7√(16L_A)/Λ²·|c₁| ∼ 2.68·10⁴√ln(Λ/8.74)/Λ²`.
  * H17 (model, `h3_H17_model.py`'s heuristic): `|δk_H17| = 424.97/Λ²` exactly; bound factor `(1 ∓ δk)^{−½}`, i.e. `±212.5/Λ²`. [M]
  * EE `δk(0)`: the first-order symbol `X − ψ₂ = Λ⁻²[½ψ₂″ − (23/2)(ψ₂/p)′]` (H7_I14 §23.2) integrates by parts to
    **`Λ²δk(0) = −168∫ψ₂p⁴ = −68.08`** [D], against the exact `−69.9 / −68.3 / −68.1` (100/300/1000) [DV to 3 % / 0.3 % / 0.02 %].
  * `κ = 0.1099` and `c_L` are flat numbers plus `O(Λ⁻²)` [D-sketch + M]; the lower end `−O(Λ⁻²)` has no computed constant.

---

## 2. Per-state positivity as a function of Δ_g

### 2.1 What is known and what fails at 100

| Δ_g | status | source |
|---|---|---|
| 100 | computed set (even `J ≤ 10`, odd 5, 7, 9): `n = 1`, min `+0.2286`, not robust (odd `J = 9` extra row `+0.03`, MC `0.005 ± 0.044`); **near-gap spins fail in both transfer models** (`J ≥ 16` with rows held: `−0.19 … −3.30`; `J = 18, 24, 28` with `Γ ∝ m`: `−0.07, −0.62, −0.57`; at `1.2–1.3Λ`) [M] | EE_DG100 §2.5 |
| 300 | holds on the computed set, min `+0.2038`; most fragile `J = 32` (`+0.783`); `H_LJ(300)` window [I] for `J ≥ 34` | RESULT_D §5.3 |
| 1000 | holds on the computed set, min `+0.2009`; `H_LJ(1000)` [D-sketch + M] | RESULT_D §5.2 |

**Mechanism of the failure at 100** (EE_DG100 §2.5, EE_NEARGAP §0.3, re-quantified here). It is the **diagonal beyond-Regge graviton
correction concentrated at small `ν` (`ν ≲ 25`)** at near-gap large spin:
* its absolute size in `c₂` units is roughly `J`-independent (`|Γ| ~ 10²–10³`; the fold window is a `J`- and `Δ`-independent fibre
  integral, HLJ §1 [D]), while `N_s ∝ b̂⁻³` with `b̂ = (2J+1)Λ/m`, so its weight per state grows like `J³`;
* a correction supported at fixed `ν < ν_c` has per-state weight `∝ ν_c²/(Λ²m²) ∝ Λ⁻⁴` at fixed `(J, m/Λ)`; at Δ_g = 100 the "window"
  `ν < 25` is a quarter of the whole `p`-range;
* `K₂` (the row engine) is valid only for `J²/m ≲ 1.5`, i.e. `J ≲ 12` near the gap at 100, so the `J ≥ 16` states there have no
  own-mass data at all; the failure is a transfer-model statement [M], not a computed one.

The uniform-window failure threshold `C*_25` (a uniform correction on `ν < 25` that flips the state; `tol.py`) makes this concrete:

| `|C*_25|` (c₂) | J = 8 | 16 | 24 | 32 | 48 |
|---|---|---|---|---|---|
| Δ_g = 100 | 228 | 59 | 26 | 16 | — |
| 200 | 2.5·10³ | 603 | 220 | 115 | 51 |
| 300 | 1.1·10⁴ | 2.5·10³ | 867 | 422 | 167 |
| 1000 | 1.2·10⁶ | 2.4·10⁵ | 7.4·10⁴ | 3.3·10⁴ | 1.1·10⁴ |

Measured small-ν graviton values are `|Γ| ≈ 430–480` (300, `J = 24, 28`, fold fits) and `≤ 350` (1000, `J = 64`, `ν = 13–33`).

### 2.2 Proxies used, and their validation

No new crystal data were computed. Each proxy combines exact per-state weights (any Δ_g, from `ee1000_common.hyb`) with corrections
measured at 300/1000.

* **Baseline `λ₀`** = the hybrid near-gap state (flat NF02 + leading-Regge curvature), min over helicity, `J ≤ 0.345(Δ_g−1)`
  (`m_gap ≤ 1.3Λ`), `m/Λ ∈ {m_gap, 1.02, …, 1.3}` [M]. *Validation:* reproduces the band/hybrid minima `0.2256 (100)`, `0.2043 (300, 1.0166Λ)`,
  `0.2009 (1000)` of EE_DG100, REVIEW_ROUTEA §4, EE_DG1000_LARGEJ exactly.
* **Weights.** `C*_ν<c` (as `ee1000_smallnu.py`), the bulk uniform tolerance (SMOOTH_BOUND), and the `H_LJ` window tolerance `T_B(c)`
  (HLJ §2, without `W_X`). *Validation:* `T_B(c=3/5)` at 300 `J = 24, 34, 64, 100`: `524/290, 390/211, 287/154, 298/158` vs HLJ
  `522/291, 392/214, 286/153, 297/159`; at 1000 `J = 48, 100, 300`: `2497, 1415, 1579` vs `2456, 1409, 1582`; `C*_50(1000, J = 48) =
  −2784` vs `−2800`; bulk `28.5 / 4.8 / 1.1` vs SMOOTH_BOUND `29 / 5 / 1.1` [DV].
* **P1 (correction transfer).** `corr_J(Δ_g) = corr_J(300)·(298.5/Λ)ᵖ` at the `m/Λ` of each 300 minimum, with `corr_J(300) = λ₀ −`
  (adopted 300 value, RESULT_D §5.3). `p = 1` (`Γ ∝ m`, EE_DG100's second transfer), `p = 2` (`Γ(x)` held: per-state weight `∝ m⁻²`),
  `p = 4` (`Γ(ν)` held at small `ν`). **Empirical `p`** from the same code at both gaps on own-mass rows (convention A, `t*_J`,
  `w2_chain_final_xc0.txt`): **median 2.83, range 1.4–3.2** over `J = 8–24` [NE]. *Validation at 100:* P1 reproduces the model
  failure at `J = 16, 24, 28` (`p = 2`: `−0.09, −0.38, −0.70`; `p = 3`: `−2.3, −3.1, −4.1`) bracketing the 100 chain (`−0.19, −3.30, −2.98`
  held; `+0.53, −0.62, −0.57` with `Γ ∝ m`); it misses the 100 chain's `J = 18, 20` failures, because the adopted 300 values there are
  almost uncorrected (direct + edge data) while the 100 chain used K₂ rows. *Validation at 1000:* predicted `0.68 / 0.75 / 0.95 / 0.97`
  (`J = 12, 14, 24, 28`) vs 1000 chain `0.68 / 0.71 / 0.91 / 0.93`.
* **P2 (sign-free large-spin window).** `min_{J ≥ √(1.5Λ)} T_B(c=3)` against a Δ-independent window correction: the refined remainder
  `≈ 250 c₂` (after removing the leading-Regge-shaped `D_h` class; measured 150–250 at 300, `(364, 64)`, `(400, 100)`) or the full
  measured size `≈ 450 c₂`. Assumption: the window correction at fixed `(β, ν)` is Δ_g-independent [M], supported by the fibre
  structure (HLJ §1, [D]) and by the 300/1000 data.

### 2.3 Margins vs Δ_g (`dgscan/analyze.py` → `cache/analyze.txt` §H)

| Δ_g | Λ | baseline min `λ₀` [M] | P1: min margin over the anchored `J ≤ 32`, `p = 2 / 2.83 / 4` [M] | window proxy `min_{J≤32}|C*_25|/450` | P2: `min T_B3/250` (refined) , `/450` | Q: `s_fail/s_d` |
|---|---|---|---|---|---|---|
| 100 | 98.5 | +0.2256 | **−0.87 / −3.64 / −15.9** (J = 32) | 0.035 | 0.32 , 0.18 | 2.76 |
| 150 | 148.5 | +0.2125 | +0.17 / **−0.46 / −2.29** | 0.11 | 0.50 , 0.28 | 3.16 |
| 200 | 198.5 | +0.2079 | +0.21 / +0.21 / **−0.05** | 0.25 | 0.69 , 0.39 | 3.32 |
| 250 | 248.5 | +0.2056 | +0.20 / +0.20 / +0.20 | 0.54 | 0.89 , 0.49 | 3.40 |
| **300** | 298.5 | +0.2043 | +0.20 (computed: **+0.2038**) | 0.94 | **1.1** , 0.62 | 3.45 |
| 500 | 498.5 | +0.2021 | +0.20 | 5.6 | 2.1 , 1.1 | 3.52 |
| **1000** | 998.5 | +0.2009 | +0.20 (computed: **+0.2009**) | 73 | 4.8 , 2.7 | 3.55 |
| 2000 | 1998.5 | +0.2004 | +0.20 | 1.1·10³ | 11 , 6.3 | 3.56 |
| 5000 | 4998.5 | +0.2002 | +0.20 | 4.5·10⁴ | 35 , 20 | 3.56 |

Reading:
* A P1 entry of `+0.20` means the binding `J = 5` twist-edge state (uncorrected to `10⁻³`) is the minimum again.
* The window proxy is pessimistic (it treats the whole measured size as uniformly negative on `ν < 25`; at 300 the measured data at
  `ν ≥ 8.3` enter with their signs and give `+0.783` at `J = 32`), but it is the sign-free statement.
* **The binding state never moves the threshold.** It sits at `0.2256 … 0.2002`, decreasing monotonically to the certified flat `0.200`.

### 2.4 The threshold Δ_g*

* **Computed:** none between 100 and 300. **Not determined.**
* **Model bracket [M]: `140 ≲ Δ_g* ≲ 300`.**
  * P1 on the anchored spins `J ≤ 32`: `Δ_g* ≈ 139 (p = 2), 172 (p = 2.83), 179 (p = 3), 205 (p = 4)`, always driven by `J = 32` near
    the gap (`m ≈ 1.15Λ`). With `p = 1` (the most favourable transfer) the anchored spins would pass even at 100, which
    contradicts both 100 transfer models and is excluded by `p_emp ≥ 1.4`.
  * The spins `32 < J ≤ 0.345Δ_g` (`J_cap = 51 / 68 / 85` at 150 / 200 / 250) have no anchor below 300. There the only argument is
    `H_LJ`: the refined sign-free window bound closes only for `min T_B3 ≥ 250`, i.e. **`Δ_g ≳ 275`** (`≳ 440` with the full 450).
    At 300 itself it is marginal (`×1.1`), which is exactly RESULT_D's "[I] in the window for `J ≥ 34`".
  * So: if the `D_h`-shape argument (HLJ §5.2, [D-sketch]) carries over to lower Δ_g, the model threshold is `≈ 140–205`
    (central `≈ 170`); if only sign-free bounds are used, the existing argument structure works from `≈ 275` up, i.e. **300 is the
    lowest gap at which the current method closes, and only marginally**.
* **Why 300 and 1000 (answer to the referee).** 300 is the lowest gap where the computed near-gap spins pass with margin
  (`≥ +0.78`) *and* the large-spin window tolerance reaches the measured window corrections; 1000 is where `H_LJ` becomes
  [D-sketch + M] with the sign-free tolerance `×4.8` (refined) above them. Both are computed; nothing in between is needed beyond the
  scaling of §2.5.

### 2.5 Does positivity, once it holds, persist for Δ_g ≥ Δ_g*? [D-sketch + M]

Per state, `λ/N_s = λ₀(J, m/Λ; Λ) + Σ_k ∫ w_k(p; J, m/Λ, Λ) Γ_k dp`, `w_k = ψ_k env_k/(m² N_s)`.
1. **Baseline.** `λ₀` decreases monotonically on the grid to the certified flat minimum `0.200` (NF02 certificate). Its minimum is
   the `J = 5` `+−` twist-edge state at all Δ_g; the excess over 0.200 is `0.0256 / 0.0043 / 0.0009 / 0.00016` at 100/300/1000/5000,
   fitted by `≈ 0.73/Λ + 1.8·10²/Λ²`. [M + NE]
2. **Fixed `(J, m/Λ)` [D for the weights, DV].** `N_s` depends only on `b̂ = (2J+1)/(m/Λ)` and `u = (Λ/m)²` (up to `O(Λ⁻²)`
   curvature), and `w ∝ ψ(p)/m²`. So a correction fixed in `x = ν/m` costs `∝ Λ⁻²` (measured exponent of the bulk tolerance
   `1.95–2.00` at `J = 5`; `→ 1.8` at `J = 16–32` by 2000), and one fixed in `ν` at small `ν` costs `∝ Λ⁻⁴` (`C*_25` exponent
   `3.85 → 4.09` at `J = 5`, `2.7 → 3.9` at `J = 32`). Even the most adverse transfer considered (`Γ ∝ m`) costs `∝ Λ⁻¹`. **Every one
   of these is monotone decreasing in Λ**; the empirical `p = 2.83` sits between them.
3. **New near-gap states as Δ_g grows.** `J_cap ≈ 0.345Δ_g` grows, while `K₂` covers only `J ≲ √(1.5Λ)`; the extra states are the
   `H_LJ` region. Along fixed `J/Λ` (fixed `β`), `ν_B = c·m/(2J+1)` is Λ-independent, `N_s ∝ b̂⁻³ ∝ Λ⁻³`, the window weight `∝ Λ⁻¹`:
   **`T_B ∝ Λ`** [D]; measured `T_B3(J = 0.1Λ) = 157 / 432 / 1404 / 2974` at 100/300/1000/2000, and `min_J T_B3 ≈ 0.30 Λ^{1.20}` [DV].
   With a Δ-independent window correction [M], the window margin grows `∝ Λ^{1.0–1.2}`.
4. **States between the 300 and Δ_g′ twist edges.** At fixed `J`, `m_gap/Λ = √(1 + 2J/(Δ_g−1)·…)` decreases with Δ_g, so a state at
   Δ_g′ > 300 with `m_gap(Δ_g′) ≤ m/Λ < m_gap(300)` has no 300 counterpart. Its baseline is computed (≥ 0.2), its correction is covered
   only by continuity in `m/Λ` (`H_cont`, [NE] at 300/1000).
5. **Large `b̂`, `m ≥ 1.3Λ` (`H_R`).** Measured corrections `−0.04 … −0.49` (300, `J = 12–28`) vs `≤ 0.005` (1000, `J ≤ 32`): a factor
   ≥ 98 for `Λ` ratio 3.35, i.e. `p ≳ 3.8` [NE].

**Conclusion.** Under the scaling hypothesis
> **`H_mono`**: at fixed `(J, m/Λ)` the beyond-leading-Regge correction in `c₂` units grows at most like `m` (relative correction
> `≲ J²/m`), and at fixed `(β, ν)` in the Regge window it is bounded uniformly in Δ_g,

positivity at Δ_g₀ implies positivity at every Δ_g ≥ Δ_g₀ for the states with a Δ_g₀ counterpart and for the `H_LJ` window (whose
tolerance grows), modulo `H_cont` for the new near-edge states. This is a **[D-sketch]** given `H_mono` + `H_cont`; `H_mono` itself is
**[M]** (consistent with all 300/1000 data: `p_emp ≥ 1.4 > 1`). It is not a proof: `Γ` is measured on rows only where `K₂` or direct
data exist, and the Δ_g-independence of the window correction is a structural argument, not a computation.

### 2.6 Which hypotheses get easier or harder with Δ_g

| hypothesis | trend with Δ_g | quantitative handle |
|---|---|---|
| A1 (gap), K1, K1′, `H_K`, Adm, H7, L0, Lemma A + D_LR, (F)+(P) structure | Δ_g-independent | — |
| **A2** (`c_T ≫ Δ_g²`) | **harder** (needs larger `c_T`) | physical requirement, not a computation |
| H17, `δk(0)`, `c₁` | easier | `∝ Λ⁻²` (§1.3) |
| `H_M` (`K₂` validity `J²/m ≲ 1.5`) | per state easier; **coverage harder** | `J_K2 ≈ √(1.5Λ)` vs `J_cap ≈ 0.345Δ_g`: the K₂-covered fraction of near-gap spins falls like `Λ^{−1/2}` (13/34 at 100, 22/103 at 300, 39/344 at 1000) |
| `H_LJ(Δ_g)` | **easier** | `T_B ∝ Λ^{1.0–1.2}`; refined sign-free form closes from `Δ_g ≈ 275`, `×2` at ≈ 490, `×4.8` at 1000 |
| `H_edge` | easier | edge-amplitude margin `×1.6` (300, `J = 34`) → `×3–19` (1000) |
| `H_R` | easier | corrections `p ≳ 3.8` (300 → 1000) |
| `H_OO` | easier | identity residuals fixed in `c₂` (`≤ 0.05 / 0.12`) → relative `∝ Λ⁻²` |
| `H_RS` (Q binding) | easier | `Δs_d(H_RS) = 0` to all orders; truncation `a⁴` term `0.015 / 0.0015` (300/1000) `∝ Δ⁻²`; RD far region `e^{−O(Δ)}` |
| De-shift existence (Q) | easier | sc3 moment ratio `1.9·10³/Λ²` per order: 0.20 / 0.086 / 0.048 / 0.021 / 0.002 at 100 / 150 / 200 / 300 / 1000 |
| `H_fwd^Q` (300 J-drift) | easier | drift tied to `J²/Δ` |
| `H_smear^Q`, `H_LJ^Q`, HYB | easier or neutral | `O(Λ⁻²)` and `J²/m` effects |
| `H_EFT`(a) (term-wise (P) of heavy-exchange double traces) | **possibly harder** [I] | toy error grows with `Δ_h` (`10⁻⁹` at 20 → `1.5·10⁻⁶` at 100; untested at `Δ_h ≥ 300`), still `≥ 10⁵×` below sensitivity at 100 |
| `κ ≥ 0`, `c_L > 0` | neutral | flat values + `O(Λ⁻²)` |

### 2.7 Result Q's positivity vs Δ_g

* Q's binding pair (`++`/OO `J = 0` at the gap) is controlled by `s_d`, measured `+0.509 (300)`, `+0.508 (1000)`. From the per-order
  engine contributions (BELOW_FOLD2 §4.2: `+1.66/−3.95/+2.52/+0.27/+0.016` at 300, `+1.66/−4.21/+2.96/+0.09/+0.002` at 1000), with a
  `1/Δ` model per order, **`s_d(Δ) ≈ 0.500 + 1440/Δ² + O(Δ⁻³)`** (the `1/Δ` coefficients cancel to `−0.1`) [NE at 300/1000; M as a law].
  Against `s_fail = 1.78`: margin `×2.8 / 3.2 / 3.3 / 3.45 / 3.55` at 100 / 150 / 200 / 300 / 1000.
* So Q's binding state is Δ_g-stable and does not set a threshold. What limits Q below 300 is not computed: the de-shift moment
  control (ratio `≥ 0.05` per order below Δ_g ≈ 200), the `J²/Δ` drift of the `+−` channels (already `−0.07` at `J = 7`, 300), `H_LJ^Q`
  at larger `J²/m`, and the import of D's `t₄` (needed only to bound `ε_{c₁}`; carrying `c₁t₄` explicitly avoids it).
* **Q's statement can therefore be extended over the same range as D's at no extra cost in grade**, and plausibly lower, but nothing
  below 300 is claimed.

---

## 3. Proposed statement

> **Theorem (conditional; any gap).** Let `𝒯` satisfy A1–A2 with gap `Δ_g` and `Λ = √((Δ_g−2)(Δ_g−1))`, and assume the analytic inputs
> of UNIFIED_THEOREM §2. If `Ω^EE`, `Ω^Q_up`, `Ω^Q_low` act positively on every twist-allowed heavy block **at this Δ_g**, then
> * `√(t₄²+α_T²) ≤ 240√(16L_A(Λ))/Λ⁴·(1 ∓ 19941 c₁/(Λ²√(16L_A)) ± 212.5/Λ²)`, with
>   `16L_A(Λ) = 3022.34·ln(Λ/8.7415) + 19293/Λ + O(ln Λ/Λ²)` (exact closed form in `bounds.py`);
> * `−O(Λ⁻²) ≤ g̃₄Λ⁶/(8πG) ≤ 17.2407·ln(2Λ/3) − 19.2782 + 10.314/Λ + (860.3 ln(2Λ/3) − 2508.6)/Λ² + O(Λ⁻³ ln Λ) + ε_{c₁}`,
>   `|ε_{c₁}| = 2.032|c₁t₄|Λ²`.
>
> **Positivity** is computed [NE + M] at `Δ_g = 300` and `1000`, giving `3.13·10⁻⁶ / 2.89·10⁻⁸` and `72.04 / 92.81`, under the
> hypotheses of UNIFIED_THEOREM §3 (including `H_LJ(300)`, [I] in its Regge window).
>
> **Corollary (conditional; range).** Assume in addition `H_mono` and `H_cont` (§2.5) [M]. Then positivity holds, and the bounds above
> apply, **for every `Δ_g ≥ 300`** [D-sketch + M]. For `300 ≤ Δ_g ≲ 490` this inherits the `H_LJ(300)` window caveat in weakening form
> (refined sign-free window margin `×1.1 → ×2`); for `Δ_g ≳ 490` the sign-free refined window bound has margin `≥ 2` [M].
>
> **Not claimed:** any `Δ_g < 300`. At `Δ_g = 100` positivity fails at model level (near-gap `J ≥ 16`); the model threshold lies in
> `140 ≲ Δ_g* ≲ 300` [M] and was not computed.

Wording for a PRL abstract/first paragraph: *"…for a twist gap `Δ_g` (`Λ ≈ Δ_g − 3/2`), `√(t₄²+α_T²) ≤ 1.32·10⁴√(ln(Λ/8.74) + 6.4/Λ)/Λ⁴`
and `g̃₄Λ⁶/(8πG) ≤ 17.24 ln(Λ/4.59) + 10.3/Λ`, established at `Δ_g = 300, 1000` and expected for all `Δ_g ≥ 300`."*
(Errors at 300: `+4.5·10⁻⁵` relative and `−0.024`, the latter on the unsafe side by 0.03 %; quote the exact `2L_AdS` (72.04 / 92.81)
where digits matter. Note that `Λ⁴` must not be replaced by `Δ_g⁴`: that alone moves the cubic bound by `−2.2 %` at 300.) If the referee asks about intermediate or larger gaps: the table of §1.2 is the implication, the table of
§2.3 shows that every computed margin improves monotonically with Δ_g.

**Grades, compactly.**

| item | grade |
|---|---|
| closed forms of `16L_A(Λ)`, `2L_AdS(Λ)`, bound values on the grid, reproduction of all quoted digits | [DV] |
| asymptotic expansions (`1/Λ`, `ln Λ/Λ²` coefficients), origin of the `B_eff` drift | [D] + [DV] |
| `c₁`, `ε_{c₁}` as functions of Λ | [D] (linearised) |
| H17 band `±212.5/Λ²` | [M] (heuristic carried over) |
| `δk(0) = −68.08/Λ²` | [D] (first-order symbol) + [DV] |
| baseline margin monotone, `→ 0.200` | [M] on the grid; limit [certified flat] |
| weight scalings `Λ⁻²`, `Λ⁻⁴`, `T_B ∝ Λ^{1.0–1.2}` | [D] + [DV] |
| correction exponent `p_emp = 2.83` (1.4–3.2) | [NE] (one convention, 300 vs 1000) |
| threshold bracket `140 ≲ Δ_g* ≲ 300` | [M] |
| persistence for Δ_g ≥ 300 | [D-sketch] given `H_mono` + `H_cont` ([M]) + `H_LJ(Δ_g)` |
| Q binding `s_d(Δ) ≈ 0.50 + 1440/Δ²` | [NE] at 300/1000 + [M] extrapolation |

---

## 4. What remains

1. **One intermediate gap computed for real**, preferably `Δ_g = 200` (P1 central threshold region): own-mass rows (`Δ ≈ 200–260`) for
   even `J = 12–32` and direct `n1/n3` nodes at `x ≤ 0.2`. `K₂` fails there for `J²/m ≳ 1.5` (`J ≳ 17`), so direct sub-rows are
   needed. This would turn the bracket into a computed threshold side.
2. **`H_mono` beyond the rows:** a measurement of `Γ` at fixed `(J, ν)` across three Δ at `J ≥ 32` (rows exist at 336/500/1600 for
   `J = 24, 28`; their small-`x` K₂ values are known to be unreliable, EE_DG300_SMALLNU).
3. **Window corrections at 300 for `J ≥ 34`** (the `H_LJ(300)` [I]); then 300 becomes the lowest gap with a closed window.
4. **The `D_h`-shape argument at lower Δ_g**; if it holds, the threshold moves to `≈ 170`.
5. Q below 300: de-shift weights with smaller high moments (the `1.9·10³/Λ²` ratio), `s_d` at 150–200 from the engine (cheap: the
   engine runs at any Δ), the `+−` `J²/Δ` drift.
6. `H_EFT`(a) at `Δ_h ≥ 300` (the only item that may get harder with Δ_g).

## 5. Reproduce

```
cd ads4/code/dgscan
ulimit -v 3000000; export OMP_NUM_THREADS=1
python3 bounds.py                                   # §1: ~1 s  -> cache/bounds.txt, bounds.json
python3 tol.py 100,150,200,250,300,500              # §2: ~3 min -> cache/tol_<Dg>.json
python3 tol.py 1000,2000,5000                       #     ~9 min (second process)
python3 analyze.py                                  # §2.3-2.7: ~6 min -> cache/analyze.txt (tables A-H), analyze.json
```
`analyze.py` parses (read-only) `spinfunc/cache/y2_chain100_jhi.txt`, `y2_chain100_jhi_scale.txt`, `w2_chain_final_xc0.txt`, and
hard-codes the adopted 300 values of RESULT_D_STATEMENT §5.3 and the per-order `s_d` of BELOW_FOLD2 §4.2 (both cited in the code).
