# Result Q, V₄″ test: the u-channel contact term with the actual de-shift weights, the O(G) twist-1 corner, and strip analyticity of the C4 weight

*21 Nov 2026, V₄″ agent. Unrefereed, not committed. No existing document or code file was edited; no `bj_*` file was touched or run. New
code: `ads4/code/resultq/v4d_core.py`, `v4d_wu.py`, `v4d_strip.py`, `v4d_og.py`; outputs `resultq/cache/v4d_wu.txt`, `v4d_strip.txt`,
`v4d_og.txt` (git-ignored). Read-only inputs: `functionals_q4_deshift.json` (stored Stage-4 weights), `functionals_q7_xprec.json` (Stage-7 XP
weights, 64 digits), `v4c_core.py`, `v4c_wu.py` (only for cross-checks). Compute: ≈ 10 CPU-minutes, ≤ 2 cores, < 1 GB.*

Grades: **[P]** proved (incl. exact algebra); **[D]** derived; **[D-sketch]**; **[DV]** derived and numerically verified; **[NE]** numerical
evidence; **[M]** model; **[I]** open.

---

## 0. Answer in brief

1. **The de-shift does not cure the u-channel contact term. With the actual weights the "cure" has no number: the de-shifted contact term
   diverges as ε → 0 in the timelike-tube prescription, for both weight sets, both pieces and both tube orientations [DV].**
   * **Exact inputs [P/DV].** The contact profile has a closed form, `U = (π²/64)T₄₂ε⁴(1+c)²(1−X²)/X⁶`, `X = −(1−c)/2 + 2iε + ε²`. It agrees
     with V₄′'s Cauchy–FFT residue to 3·10⁻⁶. The per-ν pairing is evaluated exactly through Fourier modes, at any ν. Its small-ε form is
     `G_ν(ε) = (π²/16)T₄₂ε^{−1/4}[a(ν)J(5/2+iν)ε^{iν/2} + b(ν)J(5/2−iν)ε^{−iν/2}]`, with closed-form `a`, `b`, `J`; it holds to 3·10⁻⁶ at
     ε = 10⁻⁶.
   * **Normalisation.** The weight multiplies the heavy-normalised kernel `Ψ_ν/N_com(ν)`, with `N_com = 2¹¹π²|Γ(11/4 + iν/2)|⁴ ~ ν⁹e^{−πν}`.
     Matched to `N̂₄`, this is `w(ν) = χ(ν)N_com(0)/N_com(ν)` on the V₄′ kernel.
   * **Three obstructions, each sufficient:**
     * **(i) Large ν [DV].** The dominant branch decays only like `ν^{9/2}e^{−πν/4}`. The heavy normalisation grows like `e^{πν}`, and the
       weights decay like `e^{−(π/2)ν}`. So `|w·G_ν|` grows like `e^{πν/4}`: `10⁻⁸` (ν ≤ 5) → `10` (ν = 20) → `2·10¹⁰` (ν = 40) → `8·10⁴⁵`
       (ν = 120, XP). The ν-sum is set by where the weight is truncated:

       | ν_c | 10 | 20 | 41 | 90 |
       |---|---|---|---|---|
       | `|G_ds|` (XP value, ε = 10⁻²), V₄′ kernel units | `3.5·10⁻⁸` | `0.065` | `3·10⁸` | `5·10²⁷` |

       For comparison, the raw per-ν value is `0.014`. The integral does not converge with the ν-integral inside, nor (at fixed ε) outside.
     * **(ii) The ν = 0 junction [P + DV].** Exactly, `b(ν) = a(−ν)(1 + i sinh πν)/cosh πν`. So the two branches of the ν > 0 integral do
       not continue into each other: `a₀ − b₀ = −π`. Every contour shift leaves the endpoint ν = 0, and
       `G_ds → (π²/16)T₄₂ε^{−1/4}·2πiχ₂J(5/2)/L³`, `L = ½log(1/ε)`. That is **ε^{−1/4}/log³(1/ε) → ∞**. It holds even with a smooth
       Gaussian taper that removes (i): the ratio to the asymptote is 0.46 / 0.64 / 0.78 / 0.98 at L = 10 / 20 / 40 / 60.
     * **(iii) Truncated sums.** For every truncation ν_c = 2…20 and ε = 10⁻²…10⁻²⁴ the effective exponent is `−0.15 → −0.23`. The raw value
       is `−0.31`; the claimed value was `+1.19`.
   * **Formal pole structure [NE + D].** V₄′'s contour argument is right about the poles, and only about them:
     * the weight has poles at ±i/2 and ±5i/2, and zeros at ±i and ±2i, all forced by the Abel equation (item 3);
     * these cancel exactly against the zeros and poles of `a`, `b`;
     * the only singularity left in the strip is the edge pole at ±3i, which would give ε^{5/4}log ε.
     * Obstructions (i) and (ii) are not about poles, and they decide.
   * **Prescription dependence.**
     * Tube `x₁ → x₁ ± iεe`: both diverge, with identical magnitudes.
     * Cone-preserving / light-ray-parameter prescription `t → t + iε`: exactly 0, per ν and de-shifted alike [P].
     * That prescription is the one consistent with `L[T]|Ω⟩ = 0`. The W_u Wick term contains `⟨Ω|L[T₃]L[T₁]|Ω⟩ = 0` [D-sketch].
   * **Size of the remnant.** Under the same truncated weight it is `R = 0.2–0.4` (ε = 10⁻²), `8–23` (ε = 10⁻¹²) times an O(1) regular light
     family (ν_c = 6), and `10³–10⁵` (ν_c = 20). An O(1) family is `O(C_T)` times the O(G) content. **Leverage on the bound: unbounded**
     (divergent) unless the light-ray prescription is adopted, in which case it is exactly 0.
2. **O(G), u-channel twist-1 (T) exchange, spinning [D + DV; the Einstein-specific cancellation is I].**
   * **Homogeneity lemma [D].** With light-ray polarisations every (13) leading-twist structure of a twist-τ exchange scales as `ε^{−1+τ/2}`
     in the contact region, and as `Δφ^{−2+τ}` pointwise. The corner integrand is therefore `Δφ^{τ−3/2}`. For **τ = 1 it is absolutely
     convergent per ν, margin `r^{1/2}`, for every ⟨TTT⟩ structure.**
   * **Exact light-ray structures `D₁^mD₃^n X^{−1/2}`.** These are the singular factors of the twist-1 family, which contains T. They confirm
     the lemma (§2), and the raw per-ν pairing tends to a finite limit with an `ε^{1/4}` approach.
   * **De-shifted:** **the power-counting claim `ε^{1.69}` is not confirmed.** The same obstructions (i) and (ii) apply, and the per-ν
     pointwise pairing decays only like a power of ν. In heavy normalisation the ν-sum therefore diverges with ν_c (§2).
     * It is finite only per ν (ν outside).
     * The ε-approach is `ε^{1/4}` (up to logs), not `ε^{7/4}`.
   * **Whether the four O(G) light-family residues cancel with spin as in the scalar proxy: not computed [I].** This needs the Einstein
     ⟨TTT⟩ u-channel T block and the `[TT]₀` anomalous dimensions in this gauge.
3. **Strip analyticity of the C4 weight (hypothesis 5) [NE, 80-digit Paley–Wiener on the XP weights].**
   * **The products `χR_j` are analytic in `0 < |Im ν| < 3`, with the edge pole at ±3i.**
     * `e^{3t}ĝ₃(t) = 2Φ(t)` holds for t ≥ 8 (value: 2.0; derivative: `−t²`).
     * `ĝ₁` decays at a rate going to 3 (2.90 at t = 30).
   * **`χ` itself is not analytic in the strip.** `χ̂(t) ∝ e^{−t/2}` exactly (rate 0.5000 for t = 10…30). So χ has simple poles at ±i/2,
     placed on the common zeros of all four `R_j`, and zeros at ±i, ±2i (the poles of `R₃`).
   * **Grade.**
     * Hypothesis 5 as worded ("χ analytic in |Im ν| ≤ 3") is false.
     * The form the heavy side uses, "`χR_j` analytic in the strip up to edge poles", is supported [NE].
     * The light side needs neither: it fails on (i) and (ii).
4. **Verdict.**
   * **V₄″ as stated is false** [DV]. The total de-shifted kernel with the ν-integral inside does not converge, not even as a distribution
     tested on the MFT contact term. The de-shift gains only `log³` at the light-cone-edge corner, not three powers.
   * **The upper bound survives only via a per-ν (ν-outside) light side with the light-ray prescription.** Call this V₄‴:
     * **O(1):** exactly 0 in that prescription. The regular families cancel pointwise [P, V₄′]; W_u ≡ 0 [P].
     * **O(G):** absolutely convergent per ν, margin `r^{1/2}` [D].
     * **Plus two new hypotheses:**
       * **H_LT** (definitional): the position-space realisation of `L[T₁]`, `L[T₃]` uses the light-ray (cone-preserving) regularisation at
         the collinear point. [D-sketch]
       * **H_νout:** the per-ν light side of tree ⟨TTTT⟩, combined over families at fixed ν, is polynomially bounded in heavy
         normalisation, so that `∫χ` converges. True for polynomial / CHMRS-type light sides; [I] for spinning exchanges.
   * **Grade of V₄‴: [D-sketch]**, with its O(1) part [P given H_LT].
   * The swap that makes the de-shifted functional a valid sum rule mixes "ν inside" (heavy side) with "ν outside" (light side). That is the
     review's circularity (§3.5), now realised concretely at the corner.

---

## 1. W_u with the actual χ (`v4d_core.py`, `v4d_wu.py` → `cache/v4d_wu.txt`)

### 1.1 Exact ingredients

**Contact profile [P].** In the timelike tube `x₁ → x₁ + i·sgn·εe`, with `x₁₃ = −t₁z₁ − t₃z₃ + iεe`, `e·z = −1`, `e² = −1`:
* `X₁₃ = 2t₁t₃z₁·z₃ + 2iε(t₁+t₃) + ε²`.
* The numerator `N = (z₁·z₃)X₁₃ − 2(z₁·x₁₃)(z₃·x₁₃) = ε²(1+c)`, **independent of t**. It vanishes for ε = 0.
* `Res f₄ = −(1−X²)/(4X)` (from `f₄ ≃ (v′²−1)/(4u′)`, `v′ = −1/X₁₃` at `t = ½`).
* So the edge residue of `t₁³t₃³f₄W_u` is

  `U(Δφ) = (π²/64)T₄₂ε⁴(1+c)²(1−X²)/X⁶`, with `X = −(1−c)/2 + 2i·sgn·ε + ε²`, `T₄₂ = −1/64`.

  It agrees with `v4c_wu.U_eps` to 3·10⁻⁶ (FFT precision) at ε = 10⁻³, 10⁻⁵ and Δφ = 0.3√ε … 2.

**Exact pairing.** `∫dφ₁dφ₃ C g(φ₁−φ₃) = (2π)²Σ_n f̂_n²ĝ_n`, where `f̂_n = (−1)ⁿΓ(2s+1)/(Γ(s+n+1)Γ(s−n+1))` are the exact Fourier
coefficients of `(2−2cos φ)^s`, `s = −δ/2`.
* There is no oscillatory quadrature, so any ν is reachable.
* Check: it reproduces V₄′'s `I_t(0)`, `I_t(1)` to 10⁻³ and 6·10⁻³ (V₄′'s trapezoid error).

**Corner expansion [P + DV].**
* `K_ν(Δφ) = a + bΔφ^{−2iν} + O(Δφ²)`, with
  * `a = 2πΓ(1−2δ)/Γ(1−δ)²`,
  * `b = 2^{1+2iν}[B(½,1−δ) + B(δ−½,1−δ)]/2`.
* Against `Kfun`: 1.5·10⁻⁵.
* Leading-ε pairing:
  * `G_ν = (π²/16)T₄₂ε^{−1/4}[c₊ε^{iν/2} + c₋ε^{−iν/2}]`, with `c₊ = aJ(5/2+iν)`, `c₋ = bJ(5/2−iν)`;
  * `J(σ) = 2^{σ+1}c^{α−6}Γ(α)Γ(6−α)/Γ(6)`, `α = (σ+1)/2`, `c = 2e^{−i·sgn·π/2}`.
* Against the exact pairing:

  | ε | ν = 1 | ν = 3 | ν = 8 | ν = 14 |
  |---|---|---|---|---|
  | 10⁻⁴ | `4·10⁻⁴` | | | |
  | 10⁻⁶ | `3·10⁻⁶` | `1.3·10⁻⁶` | `7·10⁻⁵` | `10⁻³` |

**Exact identity [P].** `b(ν) = a(−ν)(1 + i sinh πν)/cosh πν`, by the duplication and reflection formulas (checked to 12 digits).

**Heavy normalisation [D].**
* Stage 4: `h_j = N_com R_j`, with `N_com = 2πK_A^{ref}K_B^{ref} = 2¹¹π²|Γ(11/4 + iν/2)|⁴` (from `v3_regge_kernel`'s `K(c,λ)` with
  `c_A = 4−δ`, `λ_A = 1−δ`; `c_B = 3+δ`, `λ_B = δ`).
* **V₄′'s gauge is the CHMRS frame of v3 with `x = e`.** The dual side gives `(−z₀·e)^{−(1−δ)} = 1`.
* For the celestial side, with `z₁ ∝ ŷ+n` and `z₃ ∝ ŷ−n`, one finds `A₁₃^{(2+δ)/2}(A₁₀A₃₀)^{−δ/2} = 4|y|²(−z₀·ŷ)^{−δ}|sin θ|^{−δ}`. That
  is exactly v3's `B` structure, with all bases positive, so no ν-phase enters.
* So the heavy-normalised kernel is `C_ν/N_com(ν)` up to a ν-independent (Faddeev–Popov) constant.
* Any ν-only factor dropped consistently in both cancels.

### 1.2 Sizes along the real ν axis (XP value weight; tube +)

| ν | 1 | 3 | 5 | 10 | 15 | 20 | 40 | 80 | 120 |
|---|---|---|---|---|---|---|---|---|---|
| `|c₊|` | 0.13 | 0.15 | 0.12 | 0.029 | 3.0·10⁻³ | 2.0·10⁻⁴ | 6.5·10⁻¹⁰ | 3.3·10⁻²² | 4.6·10⁻³⁵ |
| `N(0)/N(ν)` | 1.24 | 6.4 | 112 | 2.9·10⁶ | 5.8·10¹¹ | 3.1·10¹⁷ | 1.2·10⁴² | 9·10⁹³ | 9·10¹⁴⁶ |
| `|χ|` | 8.9·10⁻⁸ | 3.7·10⁻⁸ | 1.1·10⁻⁹ | 5.5·10⁻¹¹ | 5.4·10⁻¹² | 1.6·10⁻¹³ | 2.6·10⁻²³ | 5.7·10⁻⁴⁷ | 1.9·10⁻⁶⁷ |
| `|χc₊N₀/N|` | 1.4·10⁻⁸ | 3.4·10⁻⁸ | 1.6·10⁻⁸ | 4.7·10⁻⁶ | 9.5·10⁻³ | 9.7 | 2.1·10¹⁰ | 1.7·10²⁶ | 7.8·10⁴⁵ |

* **Weight decay rates** `−dlog|χ|/dν`:
  * XP: `1.04 / 1.21 / 1.35 / 1.43 / 1.46` on `[20,30] / [30,40] / [50,60] / [80,90] / [100,110]`.
  * Stored: `1.18 / 1.46` on `[20,30] / [30,40]`.
  * Asymptotically `π/2` with a power prefactor.
  * Convergence would need `> 3π/4 = 2.36`.
* **Tube −.** `c₊ ↔ c₋` exchange magnitudes, so one branch always decays only like `e^{−πν/4}`.
* **Near ν = 0.** `χ/ν² → 3.24·10⁻⁷` (value) and `−8.40·10⁻⁸` (derivative), for stored and XP alike. So χ = O(ν²) and the ν = 0 principal
  value term vanishes. The junction term (ii) is what remains.

### 1.3 The de-shifted contact term `G_ds(ε; ν_c) = Σ_{ν_k ≤ ν_c} dν w_k G_{ν_k}(ε)`

Units are those of the V₄′ kernel. References: the raw per-ν `|G₀| = 0.014 / 0.062 / 0.125` (ε = 10⁻² / 10⁻⁴ / 10⁻⁵); `|I_t(0)| = 0.015`.

| weight (tube +) | ν_c | ε = 10⁻² | 10⁻⁴ | 10⁻⁶ | 10⁻⁸ | 10⁻¹² |
|---|---|---|---|---|---|---|
| stored value | 10 | 1.2·10⁻⁸ | 2.3·10⁻⁸ | 5.2·10⁻⁸ | 1.3·10⁻⁷ | 8.6·10⁻⁷ |
| | 20 | 7.6·10⁻³ | 0.014 | 0.030 | 0.072 | 0.49 |
| | 30 | 149 | 263 | 572 | 1.4·10³ | 9.4·10³ |
| XP value | 10 | 3.5·10⁻⁸ | 6.2·10⁻⁸ | 1.4·10⁻⁷ | 3.3·10⁻⁷ | 2.2·10⁻⁶ |
| | 20 | 0.065 | 0.12 | 0.26 | 0.63 | 4.3 |
| | 41 | 3.2·10⁸ | 5.6·10⁸ | 1.2·10⁹ | 3.0·10⁹ | 2.0·10¹⁰ |
| | 90 | 5.5·10²⁷ | 9.6·10²⁷ | 2.1·10²⁸ | 5.1·10²⁸ | 3.5·10²⁹ |
| XP deriv | 10 | 1.4·10⁻⁹ | 2.6·10⁻⁹ | 5.6·10⁻⁹ | 1.4·10⁻⁸ | 9.3·10⁻⁸ |
| | 41 | 1.3·10⁷ | 2.3·10⁷ | 5.0·10⁷ | 1.2·10⁸ | 8.3·10⁸ |

* **Effective ε-exponent** (all four weights, ν_c = 6, 10, 20, ε = 10⁻² … 10⁻²⁴): **−0.15, −0.20, −0.22, −0.22, −0.23** in 4-decade
  windows.
* **Tube −.** Identical to 3 digits for ν_c ≥ 10.

**Gaussian taper `w → w·e^{−(ν/Λ)²}`** (XP value; columns are ε):

| Λ | 10⁻² | 10⁻⁶ | 10⁻²⁰ | 10⁻⁴⁰ | 10⁻⁸⁰ |
|---|---|---|---|---|---|
| 3 | 2.6·10⁻¹¹ | 2.6·10⁻¹¹ | 4.1·10⁻⁹ | 6.6·10⁻⁵ | 4.1·10⁵ |
| 20 | 6.5·10²⁹ | 6.8·10³⁰ | 2.0·10³⁴ | | |

* Fixed Λ diverges as ε^{−1/4}/log³: this is (ii), and it matches the closed-form asymptote.
* Fixed ε diverges as Λ → ∞: this is (i).

**Relative to an O(1) regular light family under the same weight.** `R = |G_ds| / Σdν|w||I_t|`:

| ν_c | weight | ε = 10⁻² | 10⁻⁶ | 10⁻¹² |
|---|---|---|---|---|
| 6 | stored value | 0.22 | 0.92 | 15 |
| 6 | XP value | 0.36 | 1.4 | 23 |
| 6 | XP deriv | 0.19 | 0.48 | 8.3 |
| 10 | (all three) | 4–6 | 17–24 | 280–400 |
| 20 | (all three) | ≈ 1.6·10³ | ≈ 6·10³ | ≈ 10⁵ |

**Individual regular families are not integrable either.**
* `|I_t(ν)|N₀/N` = `7·10⁻³` (ν = 1), `7` (8), `3·10⁸` (20), `10¹⁵` (30).
* Only the pointwise sum `I_s + I_t ≡ 0` is harmless.

### 1.4 Where V₄′'s contour argument stands

**Lower strip (`c₊`).**
* `a(ν)` has poles at ν = 0, −i, −2i, −3i and zeros at −i/2, −3i/2, −5i/2.
* χ (§3) has poles at ±i/2, ±5i/2 and zeros at ±i, ±2i. All interior singularities of `χa` cancel.
* Leftover: the double pole at −3i, giving ε^{5/4}log ε.

**Upper strip (`c₋`).**
* `b` has zeros at i/2, 5i/2 and poles at i, 2i, 3i. Same cancellation.

**So the pole bookkeeping of V₄′ is right.** The two non-formal requirements fail:
1. A contour shift of a half-line integral leaves the segment `[0, ∓3i]`. Its contribution cancels between the branches only if
   `b(iy) = a(−iy)`, and **it does not**: the factor is `(1 + i sinh πν)/cosh πν`. Result: ε^{−1/4}/L³.
2. The integrand must decay along `Re ν → ∞` in the strip. On the real axis it grows like `e^{πν/4}`.

V₄′'s test (`δ → δ+3` on the A₁₃ power at ν = 0) saw neither.

---

## 2. O(G): u-channel twist-1 exchange at the corner (`v4d_og.py` → `cache/v4d_og.txt`)

**Lemma [D].** At the regulated corner the invariants `r = z₁·z₃`, `p = z₁·x₁₃`, `q = z₃·x₁₃` and `X = x₁₃²` are all O(ε). Every contraction
of the exchanged operator with `z₁`, `z₃`, `x₁₃` (all ≃ ∝ z) is O(1).
* A leading-twist (13) term of a twist-τ exchange is `X^{(τ−6)/2}` times a polynomial of z₁- and z₃-degree 2 each.
* Degree and dimension counting give **`ε^{−1+τ/2}` for every structure**, and `Δφ^{−2+τ}` pointwise.
  * Identity: ε⁻¹ (pointwise 0, because `N ≡ 0`).
  * Twist 1: ε^{−1/2}.
* Corner integrand: `Δφ^{τ−3/2}`. **Twist 1 is absolutely convergent per ν, margin `r^{1/2}`.**

**Exact light-ray structures.** These are the (13) legs `D₁^mD₃^nX^{−1/2}`, m, n ≤ 2, the singular factors of the twist-1 family, with closed
forms in `v4d_og.py`. Results: OG_RESULTS

**Not done [I].** The Einstein ⟨TTT⟩ u-channel T block with its `[TT]₀` partners, and hence whether the four O(G) light-family residues
cancel with spin as in the scalar Mellin proxy. With the per-ν (ν-outside) reading, that cancellation is the content of Lemma A at O(G).

---

## 3. Strip analyticity of the C4 weight (`v4d_strip.py` → `cache/v4d_strip.txt`)

**Method: Paley–Wiener on the full line.**
* `g_j(ν) = χ(|ν|)R_j(ν)`, with `g_j(−ν) = χR_{−j}(ν)`.
* `ĝ_j(t) = C[χA_j] − S[χB_j]`.
* The Abel equation is `W = ½Σ_{j=1,3}[e^{jt}ĝ_j(t) + e^{−jt}ĝ_j(−t)] = Φ`.
* XP weights at 80-digit arithmetic. `|W/Φ − 1| ≤ 3·10⁻¹⁵` up to t = 30, as Stage 7 found.

| t | 4 | 8 | 12 | 16 | 20 | 24 | 30 |
|---|---|---|---|---|---|---|---|
| rate of `χ̂` (value) | 0.90 | 0.499 | 0.5000 | 0.5000 | 0.5000 | 0.5000 | 0.5000 |
| `e^{3t}|ĝ₃|` (value) | 1.999 | 2.000 | 2.000 | 2.000 | 2.000 | 2.000 | 2.000 |
| rate of `ĝ₁` (value) | 4.7 | 2.64 | 2.82 | 2.86 | 2.88 | 2.89 | 2.90 |
| `e^{3t}|ĝ₃|` (deriv) | 16.0 | 64.0 | 144 | 256 | 400 | 576 | 900 |
| rate of `ĝ₁` (deriv) | 2.2 | 3.7 | 2.44 | 2.69 | 2.77 | 2.81 | 2.85 |

The derivative row of `e^{3t}|ĝ₃|` equals `t²`, i.e. `−2Φ`.

**Reading [NE].**
* **`g₃` is analytic in `0 < Im ν < 3`**, with the edge pole carrying exactly `2Φ` (simple pole for Φ = 1, triple for `−t²/2`).
* **`g₁` is analytic to width → 3**: its rate is still rising at t = 30, consistent with a higher-order edge pole.
* **`χ̂` decays at exactly ½**, so χ has simple poles at ±i/2. These are the common zeros `(2ν ∓ i)` of all four `R_j`.
* Since `R₃` has poles at i and 2i, analyticity of `g₃` forces χ(±i) = χ(±2i) = 0.

**Grade of hypothesis 5.**
* As worded ("continuum de-shift analytic in |Im ν| ≤ 3"), it is **false for χ** (poles at ±i/2).
* It is **true in the form the heavy side uses**, "χR_j analytic in the strip up to edge poles" [NE]. This is a statement about the grid
  solution, which Stage 7 realises to t = 30.
* Hypothesis 5 should be reworded accordingly.
* It does not rescue the light side (§1.4).

---

## 4. Verdict and changes to the hypothesis list

| claim | grade |
|---|---|
| closed form of the W_u edge profile; `N = ε²(1+c)` | [P] |
| `b(ν) = a(−ν)(1 + i sinh πν)/cosh πν`; `a₀ − b₀ = −π` | [P] |
| de-shifted W_u contact term with actual χ (stored, XP; value, deriv; tube ±) diverges: ε^{−0.23} effective; ε^{−1/4}/log³ asymptotically; truncation-dominated amplitude `e^{πν_c/4}` | **[DV]** |
| light-ray (cone-preserving) prescription: W_u ≡ 0, raw and de-shifted | [P] (algebra) |
| light-ray prescription is the one consistent with `L[T]\|Ω⟩ = 0` (H_LT) | [D-sketch] |
| V₄′'s pole bookkeeping (ε^{5/4}log) | confirmed as formal [NE + D]; not realised (junction, large ν) |
| twist-τ homogeneity lemma; twist-1 corner absolutely convergent per ν, margin `r^{1/2}` | [D] + [DV] (exact light-ray structures) |
| de-shifted O(G) `ε^{1.69}` | **not confirmed**: ν-sum diverges with ν_c; ε-approach `ε^{1/4}` |
| Einstein spinning O(G) four-residue cancellation | [I] |
| hypothesis 5: `χR_j` analytic in the strip up to edge poles | [NE]; χ itself has poles at ±i/2 |

**V₄″ → false as stated.**

**Replacement V₄‴ (per ν, light-ray prescription).** For every ν the swapped form of `Ψ_{4,ν}` equals the commutator form at the corner
`u′ → 0`, `r → 0`, `ρ → 1`.
* O(1): exactly, given H_LT [P given H_LT].
* O(G): with absolute convergence, margin `r^{1/2}` [D].
* Q's `Ĉ₄^{ds}` light side is `∫χ(ν)·(per-ν light side)dν`, which needs H_νout.
* **Grade: [D-sketch].**

**Hypothesis list.**
* **Item 3 (V₄)** is replaced by V₄‴ + **H_LT** + **H_νout**.
* **Item 5** is reworded to "`χR_j` analytic in `|Im ν| ≤ 3` up to edge poles" [NE].
* The validity argument for the de-shifted functionals (item 2) now visibly uses "ν inside" on the heavy side and "ν outside" on the light
  side. Their compatibility is the review's §3.5 circularity, and it stays [I].

**Effect on the bound.**
* With H_LT the O(1) light side is exactly 0 and nothing changes numerically.
* Without H_LT the MFT u-channel contact term is a divergent, prescription-dependent O(1)·C_T² term, of unbounded leverage. The de-shift does
  not remove it.

**Next (cheapest).**
1. Derive H_LT from the CLS-D definition: show that the f₄-weighted light transform of the W_u term is the ray-parameter limit.
2. The spinning O(G) Lemma A at fixed ν (the Einstein u-channel T block + `[TT]₀` with `γ_J`, `δp_J`) in this gauge.

---

## 5. Files (`ads4/code/resultq/`)

| file | content | output | time |
|---|---|---|---|
| `v4d_core.py` | closed-form `U`, exact Fourier pairing, `a`, `b`, `J`, `N_com`, weight loaders | — | — |
| `v4d_wu.py` | §1: checks, per-ν sizes, de-shifted sums (4 weights × tube ±), ε-exponents, taper, junction, ratio to `I_t` | `cache/v4d_wu.txt` | 30 s |
| `v4d_strip.py` | §3: Paley–Wiener test (80 digits) | `cache/v4d_strip.txt` | 15 s |
| `v4d_og.py` | §2: twist-1 structures, raw and de-shifted pairings | `cache/v4d_og.txt` | ≈ 5 min |
