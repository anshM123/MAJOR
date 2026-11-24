# De-shift closure: a fully fixed Ω^Q_up (REFEREE_REPORT §M6)

*23 Nov 2026. Unreviewed, not committed. No existing file was edited.*

* *Code (new):* `ads4/code/deshift/`:
  * `kern.py`, `ab_solve.py`, `ab_check.py`, `lr_scan.py`;
  * `moments.py`, `c5_coeffs.py`, `c5_binding.py`, `d_abel.py`, `export_data.py`.
* *Caches:* `ads4/code/deshift/cache/` (git-ignored).
* *Data:* `ads4/prl/data/` (README there).
* *Read-only inputs:*
  * `resultq/functionals_q{3,4_deshift,5_deshift,7_xprec}.json`;
  * `robust/functionals_noforward.json`;
  * `spinfunc/cache/u2_hclosed.pkl`, `u3_chi_EE_{300,1000}.pkl`;
  * read-only imports of `resultq/stage{5_positivity,7_s5check,12_forward_regge}.py` and `repro/q_binding.py`, used for the
    positivity scan and the binding-state numbers.
* *Compute:* about 1.5 CPU-hours, ≤ 2 processes, ≤ 1.2 GB.

Grades:
* **[P]** proved; **[D]** derived; **[D-sketch]**; **[DV]** derived and checked numerically;
* **[NE]** numerical evidence; **[M]** model; **[A]** assumption; **[I]** open.

"Certified" is not claimed anywhere: the bounds below are computed in 384-bit arithmetic with fine quadratures, but they are not
interval-verified.

---

> **Update (NU_TAIL_CHECK.md, unreviewed): the exponential growth above a state's mass affects NO computed number, but the published grid-weight DEFINITION at Δ_g = 1000 is defective, and the fix is free.**
> * **No effect on computed numbers [P, by inspection]:** every computed per-state number (crystal/D′, K₂ rows, direct n1/n3, GK edge nodes, windows, EDGE-R, hybrid floors) uses ν ≤ Λ ≤ m. The light-side constants are unaffected, and so is everything built from the continuum weights or the exact min-norm X: X and X^Q decay like e^{−(π/2)(ν−Λ)} beyond Λ, and tails beyond m are ≤ 4.4·10⁻⁷ (D) and ≤ 1.6·10⁻⁴ (Q) of a state's action.
> * **Problem (c):** the stored grids `prl/data/d_abel_X_1000.txt` and `prl/data/k2/xq_abel_X_1000.txt` run to Λ+312 and carry a non-decaying double-precision floor (|X| ≈ 10⁻¹⁵–10⁻¹⁶) beyond Λ. Read literally with θ(Λ+100, 10), that floor times the exact image growth (Γ ≈ 10¹⁶ at Λ+100) gives tails up to 650× (D) and 1300× (Q) a near-gap state's action. The sign is set by round-off, so (Pos) for that grid-defined object is untested at m − Λ ≲ 60–100. HK_CLOSE §5 and §6 claims are false for that object, though true for the computed values. Δ_g = 300 grids are unaffected (they end at Λ+29).
> * **Fix (changes no number):** define the functional by the continuum weights; or read the grids as zero beyond Λ+10; or use the cutoff (Λ+60, 6) for the k = 2 pieces (same leakage, tail 5·10⁻⁹ / 4.6·10⁻⁷).
> * **Growth numbers [DV]:** exact scalar Arb series. Γ(m) = 1.8·10³ and Γ(Λ+100) = 1.6·10¹⁶ at Δ_g = 1000, J = 0; the rate per unit ν is 0.13–0.19 at x = 1, rising to π/2 only as x → ∞.

## 0. Verdict

| M6 item | result | grade |
|---|---|---|
| C5 `j′ = 1` coefficients "[I]" | **Stale label.** The coefficients were computed in RESULT_Q_STAGE5 §1 [P + DV]. I re-checked every structural property the bound uses (§1.1). **I did both**: (a) the coefficients exist; (b) a bound shows the result is **insensitive to any admissible C5 weight**. The upper end 72.04 / 92.81 is exactly independent of it [D]. The LR per-state values are exactly independent [D]. Beyond LR, the binding state would need `|Γ⁽⁵⁾| ≳ 3·10⁴ c₂` (300) per harmonic to move, against `~10 c₂` for the C4 derivative piece. | [P + DV] coefficients (Stage 5); insensitivity [D] + [M] (scale of Γ) |
| k = 4, 5 de-shift "[NE] to t = 30" | **Closed for the functional actually used.** I built explicit continuum (analytic, super-Gaussian cut-off) weights that satisfy the de-shift identity **uniformly for all `t ≥ 0`**: `|W − Φ| ≤ 1.2·10⁻¹³` (C4 value), `4.2·10⁻¹¹` (C4 derivative), `2.4·10⁻⁷` (C5). Hence `|image − target| ≤ ε_W P_{−1/2}(cosh β)` for **all** β. The C4 weights equal the stored XP weights to `10⁻¹⁶` on their grid. | [DV] |
| C4 existence (exact DS) | Exact DS ⇔ a third-order difference equation with rational coefficients (step 2i), plus prescribed singularities [D]. The exact solution is not constructed analytically. **The theorem no longer needs it**: define the functional by the explicit weight; ε-DS is [DV]. Exact existence: [NE] (geometric convergence of the approximants) + [D-sketch] (compactness). | ε-DS [DV]; exact DS [NE + D-sketch]; [D] gap stated in §2.6 |
| "C5 poor beyond t ≈ 15" | **Fixed** (the new C5 weight is good for all t). With the worst-sign realised images, the HYB/LR region (β up to 45) is re-scanned. Stage 7's dip `+0.9645` (Λ = 100, β = 30) disappears. Margins equal the exact-continuum ones to `≤ 6·10⁻¹¹`. **No binding value moves.** | [DV] (images) + [M] (LR model) |
| moments "disagree by up to 5×" | **No disagreement.** The two documents measure different objects. V4F's `2.15, 82.6, 1.36·10⁴, 3.55·10⁶` are the absolute moments of the **EE3 XP** weights (the theorem's). FH's `1.49, 33, 3.5·10³, 6.6·10⁵` are those of the **sc3 Stage-4** scalar-proxy weights. Both reproduce to 3–4 digits. RESULT_Q_STATEMENT's "sc3" attribution of V4F's numbers is the error. **The growing signed high moments (`m₄ = −1.4·10⁻²`, …) belong to the Stage-4 grid weights only**. For the theorem's XP weights `|m_n − δ| ≤ 2.5·10⁻¹⁶` (n ≤ 9). For the continuum weights they are bounded by `∫|E|μ^{2n}` (§3). | [DV] |
| Abel weights never reproduced | **Reproduced independently.** D's k = 2 weights agree at 300 and 1000 (own discretisation, FFT transform, CGLS): `X` to `≤ 1.5·10⁻⁶` (ν < 2) and `≤ 3·10⁻¹¹` (ν > 25) of max ψ₂, `dk₀` to `10⁻¹⁰`, V3 exact. Q's C4 weights are reproduced by an entirely different method (shifted form, analytic basis) to `10⁻¹⁶`. Q's C5 XP weight is **not** close to any exact solution; the new one is a different gauge, with the same LR image. | [DV] |
| export | ψ₂, ψ₃ (NF02), Q3-full-δ0.2, the kernels, the XP grid weights, the continuum weights and D's Abel weights are in `prl/data/`, with a README. | — |

**Bound values and margins: unchanged.**
* `2L_AdS = 72.04 / 92.81`, the cubic bounds, the binding FF/N `+0.1427 / +0.1428`, and D's binding states.
* The only per-state numbers that change are Stage 7's realised-image LR minima, which go **up**: `+0.9645 → +0.9866` (Λ = 100, β ≥ 6),
  and `+0.9960 → +0.9977` (Λ = 300, β ≥ 12).

---

## 1. C5 spinning `j′ = 1` coefficients

### 1.1 Status
* **Computed.** RESULT_Q_STAGE5 §1.3 gives the exact `R_j` (`j = 0, ±2`) of the raw spinning `Ψ_{5,ν}[o_H q]`. The grade is [P] (sympy)
  + [DV]:
  * rank 0 reproduces q1 to `10⁻¹⁵`;
  * the rank-2 scalar reproduces the Stage-4 proxy to `3·10⁻¹⁶`;
  * an independent finite-difference q1-engine check agrees to `1.3·10⁻⁶`.
* **Where the "[I]" comes from.** It goes back to RESULT_Q_STAGE4 §3 (written before Stage 5). V4FINAL row 5 ("C5 spinning [I] as before") and
  UT/RQS hypothesis 7 inherited it. What Stage 5 actually left open is a **canonical gauge** for the 3-shift (§5 there). That is not needed
  (§1.3).
* **Checks done here** (`c5_coeffs.py` → `cache/c5_coeffs.txt`). These are not an independent re-derivation; that would need the
  `j′ = 1` light-ray monomials.
  * (i) `conj R_j(x) = R_j(−x)` exactly for all three `j`. So `Re[cΨ₅/N̂₅]` has a Hermitian extension with σ = 1.
  * (ii) `Res₀R_{±2} = ∓407025i/4` (they cancel in `N̂₅`), and `N̂₅(0) = −3122914500/2873`.
  * (iii) The real and imaginary parts of `N̂₅`'s numerator have gcd `16x² + 100`, so `N̂₅` has **no real zero** (exact).
  * (iv) `κ₁, κ₂, κ₃`, recomputed from the stored `R_j`, are `−192491/1028120`, `1252043/1028120`, `−223923/1028120`, equal to the
    stored values, with `κ₂ + κ₃ = 1` (and `H₅(0) = 1`).
  * (v) The poles and zeros of `R_j` in `|Im ν| < 2`, which the continuum C5 weight uses (§2.2).

### 1.2 Upper end: exactly C5-independent [D]
* `W[c] = Φ` near `t = 0` gives the moment identities `Re Σ_j ∫ c ĥ_j (ν−ij)^{2n} = (−1)ⁿΦ^{(2n)}(0)`, **for every solution `c`**.
* The C5 light side therefore maps `R⁴ → 0`, `g₅ → g₅` (cancelling ∂Ĉ₄'s `g₅`), `|ĝ₃|² → 0`, whatever the gauge (Stage 5 §2.3).
* `A = 17.241` and `B_eff` do not see C5 beyond these, so **72.04 / 92.81 are gauge-independent** [D].
* **For the new continuum weight** the identities hold to `≤ 2.4·10⁻⁷`, `8.3·10⁻⁶`, `1.7·10⁻³` (n = 0, 1, 2) in units of `Λ^{−2n}`. That is
  a relative effect `< 10⁻⁶` on the `g₅` cancellation [DV].

### 1.3 Per-state values: bound over all admissible C5 weights
* **Leading Regge.** Every solution has image `σ_h P_{−1/2}` (Abel injectivity) [D]. All LR / HYB-LR values are independent of the
  choice.
* **Beyond LR** (`c5_binding.py` → `cache/c5_binding.txt`).
  * Write the per-harmonic C5 image as `P(ν−ij) + Γ_j/m²` (Γ in `c₂` units).
  * Then `|δ(FF/N)| ≤ a u³ A₅ G/(m²N)`, with `A₅ = ∫|c|Σ_j|ĥ_j|` and `G = sup|Γ|` [D].
  * At the binding `++` J = 0 gap state (`u = 1`, `m² = Λ²`, `a = 0.1694`, `N = 3.008`), the tolerance `G* = (FF/N)NΛ²/(aA₅)` is:

    | `Δ_g` | `A₅ = 1.80` (XP) | `1.91` (new) | `7.66` (worst Stage-6 gauge) | XP → new gauge change |
    |---|---|---|---|---|
    | 300 | `1.26·10⁵ c₂` | `1.18·10⁵ c₂` | `2.95·10⁴ c₂` | `|Δ(FF/N)| ≤ 3.0·10⁻⁷ G` |
    | 1000 | `1.41·10⁶` | `1.32·10⁶` | `3.3·10⁵` | `≤ 2.7·10⁻⁸ G` |

  * **For comparison:** the C4 derivative piece's tolerance at the same state is `G_U ≈ 9–11 c₂` (Stage 7 §2.2). The C5 piece is
    `Λ⁻²`-suppressed like the C4 value piece.
  * Measured beyond-LR corrections of the k = 4 images are `10–330 c₂` (HLJ_CLOSE, HRS). With `G = 10³ c₂`, any admissible C5 weight
    with `A₅ ≤ 7.7` moves the binding FF/N by `≤ 5·10⁻³` (300) / `5·10⁻⁴` (1000), against a margin of 0.143.
* **Grade:** insensitivity [D], given only that the unmeasured C5 beyond-LR corrections are `≪ 10⁴ c₂` [M, natural scale].
  This agrees with Stage 6's "C5 gauge changes minima by ≤ 0.003" [M].

---

## 2. The k = 4, 5 continuum de-shift for all t

### 2.1 The de-shift identity in shifted form [D]
For a weight `u` (C4: χ, even, real on ℝ; C5: `D := c/N̂₅`, Hermitian):

    W[u](t) = ½ Σ_j ∫_ℝ u(ν) R_j(ν) cos((ν − ij)t) dν.

1. **Shift.** Shift each `j`-term to `μ = ν − ij`.
2. **Assumptions.**
   * `uR_j` is analytic in the strip crossed by the shift (between ℝ and ℝ+ij), apart from singularities cancelled exactly by zeros of
     `R_j` (§2.2).
   * `u` has a pole at `ν = ±i j_max`, i.e. `μ = 0`.
   * `u` decays at `|Re ν| → ∞` uniformly in the strip.
3. **Result.** Then, for **every real t**,

       W[u](t) = Σ_{j=±j_max} (iπ s_j/2)(a₁^{(j)} − a₃^{(j)} t²/2) + ½∫_ℝ E(μ) cos(μt) dμ,

   where:
   * `E = Σ_j [f_j − S_j]`, `f_j(μ) = u(μ+ij)R_j(μ+ij)`, and `S_j` is the Gaussian-regularised principal part of `f_j` at 0;
   * `s_j = −1` (j < 0), `+1` (j > 0), and `a_k` are the Laurent coefficients;
   * a double-pole coefficient `a₂` must vanish.
4. **Consequence.** Imposing "residue part = Φ" exactly gives:
   * `|W(t) − Φ(t)| ≤ ε_W := ∫₀^∞|E_even| dμ` for all `t ≥ 0`;
   * `|image(β) − target(β)| ≤ ε_W P_{−1/2}(cosh β)` for all β (since `|P_{−1/2+iμ}| ≤ P_{−1/2}`);
   * the light-side moments are exact up to `|Δm_n| ≤ ∫₀^∞|E_even|μ^{2n}`.
5. **Equivalently:** exact DS ⇔ `Σ_j u(μ+ij)R_j(μ+ij) ≡ 0` on ℝ∖{0}, with the pole normalisation. This is a third-order linear
   difference equation in steps of 2i with rational coefficients.
6. **Check [DV].** For every weight built, the formula agrees with a direct real-axis quadrature of `W(t)` (no shifting, 384-bit) at
   t = 0 … 10, to all printed digits: e.g. `1.4976578·10⁻²¹` vs `1.4976574·10⁻²¹` at t = 0 for the C4 derivative.

### 2.2 Singularity structure (the "better-conditioned basis") [D + NE]
* **Forced by the identity [D].**
  * χ has zeros at `0` (double), `±i`, `±2i`, which cancel `R_{∓3}`'s poles in the crossed strips.
  * The C4 targets need a pole at `±3i`: simple for the value, triple for the derivative (with `a₂ = a₁ = 0`). This is the
    `e^{3t}F_{−3} → 2Φ` structure HK_CLOSE saw.
* **Allowed and present [NE].** χ may also have simple poles wherever **every** crossed `R_j` vanishes. For C4 these are `±i/2`, `±5i/2`
  and `±i(3+2√2)/2` (the last is a root of `4ν² ∓ 12iν − 1`).
  * An AAA rational fit of the stored XP χ on the real axis finds `±0.5000i`, `±2.4988i`, and a 4-pole cluster at `±2.95i`, i.e. the
    pair {2.914i, 3i}.
  * Without these poles no fit works: `ε_W ≈ 1–10` (first attempt). With them it converges.
* **C5** (two real unknown functions, one equation): I use `D = iν(ν²+1)/(ν²+4)·Π_a i/(ν−ia)·B`, with `a` the zeros of the crossed
  `R_{∓2}`: `−1.19593…` (a quintic root, 125 digits), `3/2`, and `√73 − 15/2`.
* **Basis.** Gaussian RBFs (`h = 0.15`, `σ = 0.45`, centres to ν = 60), solved by least squares on `E` (μ ∈ (0, 80]) plus exact
  constraint rows. A weak gauge row pulls towards the XP weight (C4) or towards min-norm (C5). Arithmetic is 384-bit
  (`ab_solve.py`).

### 2.3 Results (`ab_check.py` → `cache/ab_check_{c4v,c4d}_h015.txt`, `ab_check_c5_p015.txt`) [DV]

| piece | `ε_W` (all t ≥ 0) | residue part | `max|w − w_XP|/max|w_XP|` (XP grid) | amplification `A₀` new / XP | `|w|/max` at ν = 40 / 60 |
|---|---|---|---|---|---|
| C4 value | **`1.2·10⁻¹³`** | `1 + O(10⁻¹⁵)` | **`1.9·10⁻¹⁷`** | 2.152410 / 2.152410 | `1.6·10⁻¹⁶` / `4·10⁻³⁰` |
| C4 derivative | **`4.2·10⁻¹¹`** | t²: `−½`; const `10⁻⁴³`; `a₂ = 8·10⁻⁴²` | **`1.1·10⁻¹⁶`** | 0.173185 / 0.173185 | `4·10⁻¹⁷` / `8·10⁻³¹` |
| C5 | **`2.4·10⁻⁷`** | `1 + O(10⁻¹⁵)` | 0.37 (a different gauge) | 1.915 / 1.797 | `2·10⁻⁵` / `2·10⁻⁶` |

* **Convergence of the C4 value [NE].** As the basis is refined (`h, σ` = `0.5, 1.2` → `0.25, 0.6` → `0.2, 0.5` → `0.15, 0.45`):
  * `ε_W` goes `11 → 1.7·10⁻⁷ → 7.2·10⁻⁹ → 1.2·10⁻¹³` (the first value is the run without the R-zero poles);
  * the distance to XP goes `— → 3.6·10⁻¹⁰ → 1.5·10⁻¹¹ → 1.9·10⁻¹⁷`.
* **Convergence of the C4 derivative:** `1.9·10⁻⁶ → 1.4·10⁻⁶ → 4.2·10⁻¹¹`.
* **Direct check at t = 10.** `W − Φ = 1.5·10⁻¹⁹` (value), `5.4·10⁻¹⁶` (derivative), `5.7·10⁻¹²` (C5).
* **Comparison with the stored weights.** XP fails at `t ≈ 41`, Stage-4 at `t ≈ 12`, and the XP C5 at `t ≈ 20`.
* **Caveat (C5 pole location).** The C5 pole at the quintic root is stored to 125 digits. The left-over residue (`≲ 10⁻¹¹⁵`) gives a term
  `≲ 10⁻¹¹⁵e^{0.81t}`, so the bound as computed holds for `t ≲ 300`. It holds **for all t** if `a₅` is defined as the exact root, which
  is how the functional should be stated (`prl/data/README.md`). The C4 pole `(3+2√2)/2` is algebraic and is evaluated exactly.

### 2.4 What this does for Ω^Q_up and H_K
* **Definition.** Replace the C4/C5 parts of `Ω^Q_up` by the continuum weights (`prl/data/q_deshift_continuum_*.json`).
  * For C4 this is **the same functional** as the stored XP grid weights to `10⁻¹⁶` relative.
  * Every computed per-state value moves by `≤ dA = ∫|w−w_XP|Σ|R_j| ≤ 6·10⁻¹⁶` of the normalisation.
  * For C5 it is a gauge change, bounded in §1.3.
* **The cutoff.** These weights are entire × rational with super-Gaussian decay beyond ν ≈ 60. They are therefore cut-off weights in
  HK_CLOSE's sense. Multiplying by HK's erf plateau (`ν_c = 150`, `a = 6`) changes them by `< e⁻³⁰⁰` relative. `w/N_com` converges
  (`N_com⁻¹ ~ e^{πν}`), so K2ν's kernel is absolutely convergent [D].
* **H_K, k = 4 corner.** Proposition K_c needed `|I_c − I| ≤ ε_c(P_{−1/2} + |D|)` "given the NE structure" of Lemma S4. Here
  `|I_c − I| ≤ ε_W P_{−1/2}` holds **directly**, for all β [DV]. The corner count then gives `J₀ < 7/2` for `Ĉ₄^{ds}`, `∂Ĉ₄^{ds}` without the
  NE input.
  * For `Ĉ₅^{ds}` it gives the same `J₀ < 7/2`, against raw Fubini's `3/2`.
  * K1′, LRC and K3 are untouched, so H_K's net grade is unchanged ([D-sketch] by HK_CLOSE's weakest links).
* **Light side.** `|Δm_n| ≤ ∫|E|μ^{2n}` (§3). The R⁴ normalisation `2g₄` and the graviton pole move by `≤ 10⁻¹²` (C4) [DV].

### 2.5 The HYB / LR region re-evaluated (`lr_scan.py` → `cache/lr_scan.txt`) [DV images + M model]
* **Set-up.** Stage 5's LR scan (Stage 7's `scan`, both δ, three Λ), with β extended from 30 to **45**.
* **Realised images.** Worst-sign envelopes from §2.3, over all 8 sign patterns:
  * value and C5: `|e| ≤ ε_W`;
  * derivative: `|e| ≤ min(2ε_W/t², 2c₀/t² + m₁)`, with `c₀ = 1.5·10⁻²¹`, `m₁ = 2.1·10⁻¹⁰`.

| Λ | δ | min, b̂ ≥ 20 | min, β ≥ 6 | min, β ≥ 12 (to β = 45) | Stage 7 XP (β ≤ 30) | largest worst-sign shift |
|---|---|---|---|---|---|---|
| 100 | 0.2 | +0.2972 | **+0.9866** | **+0.9932** | +0.9645 (β = 30, C5 tail) | `−5·10⁻¹¹` |
| 300 | 0.2 | +0.2260 | +0.9953 | +0.9977 | +0.9960 | `−1.3·10⁻¹⁴` |
| 1000 | 0.2 | +0.3567 | +0.9992 | +0.9999 | +0.9996 | `−5·10⁻¹¹` |
| 100 / 300 / 1000 | 0.05 | +0.2482 / +0.0653 / +0.1910 | +0.9895 / +0.9964 / +0.9993 | +0.9948 / +0.9982 / +0.9999 | +0.9645 / +0.9960 / +0.9996 | `≤ 6·10⁻¹¹` |

* No negative value appears anywhere to β = 45.
* The minima at β ≥ 12 sit at β = 12, so the extended range adds nothing lower.
* Beyond β = 45 the image bound holds for all β. The LR model's own validity is the remaining limit there (HYB [M]).
* **Not affected.** The near-gap / binding states use the exact pure images (HYB, `B₅ = σu³P₀`), so they are unchanged; C5 beyond LR is
  covered by §1.3.

### 2.6 Existence of the exact C4 solution
* **[D] reduction.** Exact DS for C4 ⇔ a solution χ of `Σ_{j=±1,±3} χ(μ+ij)R_j(μ+ij) = 0` (μ ∈ ℝ∖0). The solution must be:
  * even and real;
  * meromorphic in `|Im ν| ≤ 3`, with only the §2.2 singularities;
  * decaying, and normalised by the `±3i` residue (§2.1).
* **[NE] evidence.** The approximants converge geometrically (§2.3) to the XP weight (`ε_W → 10⁻¹³`), with bounded norm (`A₀` constant to
  7 digits).
* **[D-sketch] limit.** A sequence with `ε_W → 0` and uniformly bounded weighted norms has a weak-* limit. `W[·](t)` is continuous for
  each t in that topology, so the limit solves the identity exactly. The uniform bounds are the [NE] input.
* **Gap to [D].** Solvability of the third-order rational difference equation in that class: a Riemann–Hilbert / Nörlund-type statement
  with index count. I did not attempt it (hard).
* **Why it no longer matters.** With `Ω^Q_up` defined by the explicit weights, the theorem needs only ε-DS, which is [DV]. Every place
  "DS" enters (images, moments, H_K corner) is a bound in `ε_W`. **Recommended change:** replace hypothesis "De-shift existence [NE to
  t = 30; C5 [I]]" by "ε-DS for the explicit weights [DV]: `ε_W = 1.2·10⁻¹³ / 4.2·10⁻¹¹ / 2.4·10⁻⁷`".

---

## 3. The moment "disagreement" (`moments.py` → `cache/moments.txt`) [DV]

Absolute moments `|m|_n = ∫|χ|Σ_j|R_j||ν−ij|^{2n}`, value weight, n = 0–3:

| object | `|m|₀ … |m|₃` | matches |
|---|---|---|
| **EE3, Stage 7 XP** (`ν_max = 122.9`; the theorem's) | **2.152, 82.55, 1.362·10⁴, 3.554·10⁶** | V4FINAL §3 (2.15, 82.6, 1.36·10⁴, 3.55·10⁶) |
| EE3, continuum weight | 2.152, 82.55, 1.362·10⁴, 3.554·10⁶ (identical) | — |
| EE3, Stage 4 (`ν_max = 41`) | 2.00, 49.8, 5.09·10³, 9.79·10⁵ | — |
| **sc3** (scalar proxy), Stage 4 | **1.49, 33.2, 3.53·10³, 6.58·10⁵** | FERMION_HEFT §4.4 (1.49, 33, 3.5·10³, 6.6·10⁵) |

* **Which is right.** Both are right. They are moments of different kernels and weights. The ratio grows to 5.4 at n = 3 because the EE3
  XP weight's `|χ|Σ|R_j|` tail (`e^{−~1.3ν}`) is wider than sc3's.
* **The actual error.** RESULT_Q_STATEMENT App. A ("V4F … sc3") and hypothesis row 7 ("the stored sc3 weights … m₄ = −1.4·10⁻²") mix
  them up. FH says it was "a different normalisation or key": it is a different key and a different weight.
* **Signed high moments.** `m_n = Σ_j∫χR_j(ν−ij)^{2n}`, which should equal `δ_{n0}`:

  | weight | `m₄` | `m₅` | `m₆` | `|m_n − δ|`, n ≤ 9 |
  |---|---|---|---|---|
  | sc3, Stage 4 | `−1.4·10⁻²` | `−27` | `−5.2·10⁴` | grows ×1.9·10³ per order (`≈ 1.16 ν_max²`): grid cutoff, as FH says |
  | EE3, Stage 4 | `−2.4·10⁻²` | `−50` | `−9.8·10⁴` | same |
  | **EE3 XP** (theorem) | `−2·10⁻³⁷` | `−4·10⁻³³` | `−6·10⁻²⁹` | **`≤ 2.5·10⁻¹⁶`** |
  | **continuum** (certified bound `∫|E|μ^{2n}`) | `≤ 1.1·10⁻⁵` | `≤ 3.6·10⁻²` | `≤ 1.2·10²` | n = 0…3: `1.2·10⁻¹³, 6.9·10⁻¹³, 2.2·10⁻¹¹, 4.8·10⁻⁹` |

  The XP grid sums cancel their individually large tails exactly: e.g. the `ν > 30` part of `m₉` is `10²⁰`.
* **"Cutoff-controlled".** For the continuum weight, `|m_n − δ|` is set by `E`'s tail at the edge of the basis (ν ≈ 60). It grows like
  `60^{2n}·ε_tail`, which is HK_CLOSE's `ε_c ν_c^{2n}`. It enters a degree-2n EFT term as `Λ^{−2n}`. At n = 8, `1.4·10⁹/Λ¹⁶ ≈ 4·10⁻³¹`
  (300) is negligible.
  * The absolute moments are an intrinsic property of the continuum weight. They are not a grid artefact: the ratios
    `|m|_{n+1}/|m|_n = 38, 165, 261, 307, …, 533` follow the weight's exponential decay.
* **Consequence for RQS row 7.** The sentence "the stored sc3 weights have growing high moments … harmless for Λ ≥ 300" is true of the
  Stage-4 weights. **It does not describe Ω^Q_up**, whose signed moments are exact to `10⁻¹⁶` (XP) or `ε_W`-controlled (continuum).

---

## 4. Independent Abel inversion [DV]

### 4.1 Result D, k = 2 (NF02 ψ₂, EE kernel, min-norm gauge) (`d_abel.py` → `cache/d_abel.txt`)
* **Method.** Own discretisation:
  * `dν = 0.025 / 0.05` (stored: 0.02 / 0.04), same `ν_max`;
  * `2N` offset collocation rows (stored: square system);
  * own FFT type-IV transforms and CGLS from zero. Nothing from `spinfunc/` or `h7i14/` is imported.
* **Inputs.** `R_∓` and `ψ₂` are read as data. The gauge `Y − ψ₂ ⊥ null` in `L²(dν)` was made explicit; it depends on the phase
  convention of `n`, which is now recorded in `prl/data`.

| | Δ_g = 300 | Δ_g = 1000 |
|---|---|---|
| CGLS iterations / relative residual | 19 / `6.3·10⁻¹⁶` | 20 / `5.7·10⁻¹⁶` |
| `max|X_new − X_stored|/max ψ₂`, ν < 2 / 2–10 / 10–25 / 25–Λ / > Λ | `3.8·10⁻⁷ / 9.5·10⁻⁸ / 4.4·10⁻¹⁰ / 3.2·10⁻¹¹ / 1·10⁻¹²` | `1.5·10⁻⁶ / 1.1·10⁻⁷ / 5.3·10⁻¹⁰ / 3.2·10⁻¹¹ / 2·10⁻¹³` |
| `dk₀` new vs stored | `−7.666324·10⁻⁴` (rel. diff `4·10⁻¹¹`) | `−6.830150·10⁻⁵` (`9·10⁻¹¹`) |
| graviton pole V3 | exact to `6·10⁻¹⁶` | `8·10⁻¹⁵` |
| off-grid `|E − Φ|`, new / stored, at t = 10 / 18 / 22 | `1.6·10⁻⁹ / 9·10⁻⁶ / 3.6·10⁻⁴` vs `1.3·10⁻⁹ / 2.8·10⁻⁶ / 7·10⁻⁴` | similar; both fail at `t ≈ 22` (double-precision floor, HK §3) |

* The small-ν difference (`≤ 1.5·10⁻⁶`) is at the level of the two discretisations' `O(dν²)` errors.
* The bound factor `(1+dk₀)^{−½}` agrees to 7 digits.

### 4.2 Result Q
* **C4 value and derivative.** These are reproduced by a different formulation: the shifted difference equation, an analytic RBF basis,
  and 384-bit least squares instead of DCT collocation with SVD and iterative refinement. They agree with the stored XP weights to
  `1.9·10⁻¹⁷` and `1.1·10⁻¹⁶` relative (§2.3).
* **C5.** It is **not** reproduced, and cannot be:
  * the XP weight is rough at grid scale (AAA finds a line of near-real poles), which reflects its stalled null-space content;
  * it is not within 0.3 of any smooth exact solution;
  * the continuum C5 is a different gauge, with the same LR image and moments. Its effect beyond LR is bounded in §1.3.

---

## 5. Exported data (`ads4/prl/data/`, README there)
* **Flat weights.** `nf02_psi.json` (exact rationals, `k₀`); `q3_full_delta0.2.json`.
* **Kernels.** `kernels.json`: EE3 `R_j`, C5 `R_j`, D `R_∓` and `q²` (sympy strings).
* **XP grid weights.** `q_deshift_xp_{C4_value,C4_deriv,C5}.txt`, 64 digits: Ω^Q_up as fixed on 22 Nov.
* **Continuum weights.** `q_deshift_continuum_{c4v,c4d,c5}.json`: formula, 60-digit coefficients and certificate. Also `*_table.txt`.
* **D's Abel weights.** `d_abel_X_{300,1000}.txt`: X, Y, ψ₂ on the stored grids.
* `MANIFEST.txt` (SHA-256 prefixes).

---

## 6. Suggested wording changes (not made here)
1. **UT §3.1 "De-shift" row; RQS hypothesis 7; V4F row 5.** Replace with: "ε-DS for the explicit continuum weights [DV]:
   `|W − Φ| ≤ 1.2·10⁻¹³ / 4.2·10⁻¹¹ / 2.4·10⁻⁷` for all t (C4 value / derivative / C5). C5 coefficients [P + DV] (Stage 5).
   Exact existence [NE + D-sketch], not needed." Delete "C5 spinning [I]".
2. **RQS App. A, last paragraph.** "V4F's moments are those of the EE3 XP weights, FH's of the sc3 Stage-4 weights; no disagreement
   (DESHIFT_CLOSE §3)." **Row 7:** the `m₄, m₅` numbers refer to the Stage-4 sc3 weights, not to Ω^Q_up.
3. **Definition of Ω^Q_up** (UT §1, RQS §1). "… realised by the de-shifted pieces with the continuum weights of
   `prl/data/q_deshift_continuum_*.json` (equal to the Stage-7 C4 weights to `10⁻¹⁶`; C5 in the gauge given there)". This also
   implements HK_CLOSE's required cutoff.
4. **Letter / referee response, M6.** Every item above is addressed.
   * The weights are published.
   * D's Abel inversion and Q's C4 de-shift are independently reproduced.
   * The C5 "[I]" is a stale label, and the bound is shown insensitive to the C5 choice.

## 7. What remains
1. **An analytic existence proof for the exact C4 solution** (§2.6 gap). It is not needed for the theorem as reformulated.
2. **D's k = 2 continuum weight.** I did not build it. The stored double-precision weights fail at `t ≈ 22` (HK §3).
   * HK_CLOSE's Lemma S2 [D] + cutoff covers this analytically.
   * A continuum k = 2 weight needs a basis that resolves the non-analytic ψ₂ target on `[0, Λ]` (about 2000 unknowns). Feasible, not
     done.
3. **The size of the C5 (and C4) beyond-LR corrections Γ.** These are still unmeasured. §1.3 only shows that C5 would need
   `≳ 3·10⁴ c₂`.
4. **Certification.** The `ε_W` are 384-bit midpoint quadratures of smooth, Gaussian-decaying `|E|` on a 0.01 grid, not interval
   bounds. Upgrading to [certified] is mechanical: arb quadrature of `|E|` plus a tail bound.
5. **Git note.** Early versions of some `code/deshift/` files appear in another agent's commit (`efc32f8`). I committed nothing.

## 8. Reproduce
    cd ads4/code/deshift; ulimit -v 3000000; export OMP_NUM_THREADS=1
    AB_TAG=h015 python3 ab_solve.py c4v 0.15 0.45 400 1e-4                 # ~3 min
    AB_TAG=h015 python3 ab_solve.py c4d 0.15 0.45 400 1e-4                 # ~3 min
    AB_C5POLES=1 AB_TAG=p015 python3 ab_solve.py c5 0.15 0.45 400 1e-6     # ~6 min
    python3 ab_check.py c4v h015; python3 ab_check.py c4d h015; python3 ab_check.py c5 p015   # ~4 min each
    python3 lr_scan.py        # ~2 min      python3 moments.py; MOM_CONT=1 python3 moments.py   # ~3 min
    python3 c5_coeffs.py      # ~1 min      python3 c5_binding.py                                # s
    python3 d_abel.py 300 1000   # ~10 s    python3 export_data.py h015 h015 p015               # ~1 min
