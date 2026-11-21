# Result Q, Stage 7: what a measurement of `K` must measure, why it was not measured, and extended-precision de-shift weights

*21 Nov 2026, Stage-7 agent. Unrefereed, not committed by this agent. No existing document or code file was edited. New code:
`ads4/code/resultq/stage7_null.py`, `stage7_kscan.py`, `stage7_xprec.py`, `stage7_xpcheck.py`, `stage7_s5check.py`. New functional file:
`ads4/code/resultq/functionals_q7_xprec.json` (extended-precision weights, §4). Outputs in `resultq/cache/stage7_*` (git-ignored). Read-only
imports: `stage4_hj`, `stage4_deshift`, `stage5_hj`, `stage5_deshift`, `stage5_positivity`, `stage6_model`, `stage6_scan`, `stage6_stress`,
`stage6_gauge`, `stage12_forward_regge`, `common_q`. Nothing in `spinfunc/` was run or imported. Compute: about 1 CPU-hour on ≤ 2 cores, ≤ 1.5 GB.*

Grades: **[P]** proved, **[D]** derived with the stated assumption, **[DV]** derived and checked numerically, **[M]** model / numerical,
**[NE]** numerical evidence, **[I]** open. Flat = D = 4 amplitudes; AdS = the CFT₃ functionals. `FF/N` uses Stage 0's normalisation (as
Stage 6). Functionals: δ = 0.05 = `Q0-full-robust`, δ = 0.2 = `Q3-full-delta0.2`. "Gauge" = a solution of the de-shift Abel equation
`W[χ] = Φ`; H0 … H4, SB as in Stage 6; **XP** = the Stage-7 extended-grid weight (§4).

---

## 0. Answer in brief

1. **`K` was not measured [I].** No available engine reaches the beyond-LR image of `Ψ_{4,ν}[EE]` where the de-shift weights live, at the
   precision the decision needs (§3).
   * **Where.** The derivative weight's response has **72 % of its mass at `ν < 2.5`** and 24 % at `2.5–8.3` (§2).
   * **What is available.** The crystal saddle folds at `ν* ≈ 11.3` for every `(Δ, J)`. Below the fold:
     * the Laplace (`K₂/K₃`) expansions do not exist;
     * the n3 contour is validated only down to `ν = 8.3–10`;
     * n4 thimbles reach `ν = 5`, with errors of `±250` to `±2·10⁴` c₂;
     * nothing reaches `ν < 5`.
   * **What is needed.** About **3 c₂** (δ = 0.05) or **10 c₂** (δ = 0.2) at the `++`, J = 0 gap state.
   * So the requested exact image, and both validations (large-J → Stage 4's LR; scalar → CHMRS skeleton), could not be computed. I did not
     build the `k = 4` (`× R`) integrand, because it would not change this reach (§3.1).
2. **What a measurement of `K` has to deliver, exactly [DV].** The gauge freedom of `Ĉ₄` is spanned by **three** directions:
   * 2 exact null directions + 1 near-null one;
   * supported at `ν ∈ [2.5, 18]` (≥ 94 % of their response weight), around the real zeros 4.55, 8.17, 13.09 of `N̂₄`.
   * Every Stage-6 gauge difference (H2, H3, H4, SB) lies in their span to `≤ 3·10⁻⁷`. The gauge spread of any state's image is therefore
     **three numbers**: the null projections of the exact image on `ν ∈ [2.5, 18]`.
   * The absolute (single-gauge) correction needs the image on all of `ν ≲ 18`, mostly below 2.5. Tolerances per binding state are in §2.2.
3. **Verdict (task 2): not decided; conditional statements only [M].**
   * In the covariant models (`K = 0`) both δ are positive: Stage 6's `+0.0478 / +0.0497 / +0.0491` and `+0.1985 / +0.1977 / +0.1970`.
   * With Stage 6's non-covariant structures at strength `|K|` (§3.3):

     | δ | `|K| = 1` (natural scale `c₄^{++}`), worst gauge | positive in every gauge up to |
     |---|---|---|
     | 0.05 | min `FF/N` = **−0.051** (SB, structure b); `−0.028` (H0/H3, structure a); only XP stays non-negative (+0.0002) | `|K| ≲ 0.5` |
     | 0.2 | min `FF/N` = **+0.068** (SB, b); `+0.098` (H0/H3, a) | `|K| ≲ 1.5` |

   * The binding state is the `++`, J = 0 gap state throughout.
   * **δ = 0.2 remains the value to quote.** It is robust to the natural size `|K| ~ 1` in every gauge; δ = 0.05 is not. Neither is
     established, because `K` itself is [I].
4. **Extended precision (task 3).**
   * **Correction to Stage 6 §2.4 [NE].** The stored weights' large-`t` failure is **truncation at `ν_max = 41`, not double-precision round-off**.
     * A 384-bit solve of the same discrete system cuts the on-grid residual by 10¹⁰.
     * The off-grid `W(t)` is unchanged: `|W/Φ − 1|` at t = 10.1 is `5·10⁻⁷` for the 384-bit solve and `2·10⁻⁶` for the stored weight
       (derivative).
   * **C4 [DV].** On `dν = 1/25`, `ν_max = 123`, the min-norm weights refined to 384 bits satisfy `W = Φ` off-grid to
     **`4.4·10⁻¹⁸` (derivative) and `3.4·10⁻¹⁵` (value) at t = 30**. The light-side moments are exact (`≤ 2·10⁻¹¹`) and `∫χN̂ = 1`.
     * Gauge change from the stored weights: `3.5·10⁻³` / `1.4·10⁻²` (max relative). That is harmless on covariant images (`≤ 8·10⁻¹¹`)
       and shifts the `K`-stress coefficients (`s_a^{4d}`: −0.90 → −0.58).
   * **C5 [NE]: not achieved.** With no spectral gap, refinement stalls at a residual floor of `3·10⁻¹⁹`. The best weight gives
     `|W − 1| = 5·10⁻⁷ / 8·10⁻⁴ / 1.1` at t = 15.3 / 20.2 / 22.7 (stored: `10⁻⁴` at 14).
   * **Stage 5 re-run with realised images [NE + M]:**
     * **stored weights:** `FF/Nrm < 0` from **β = 15.6** (`++`, `J ≈ 4·10⁵–10⁷`, all Λ, both δ), reaching −1 by β ≈ 29;
     * **Stage-7 weights:** positive to β = 30, min (β ≥ 6) **`+0.9645 / +0.9953 / +0.9992`** (Λ = 100 / 300 / 1000);
     * the drop from `+0.9895` at Λ = 100 comes from the unrealised C5 tail at β ≈ 30.

---

## 1. What has to be measured: the gauge space of `Ĉ₄` (`stage7_null.py` → `cache/stage7_null.txt`)

### 1.1 Null space [DV]
Discrete Abel operator of Stage 4 (EE3, rows normalised by `cosh 3t`). The smallest singular values relative to `s_max` are:

| grid | smallest singular values / `s_max` |
|---|---|
| `dν = 0.02`, N = 2048 | `1.8·10⁻¹⁶, 3.6·10⁻¹⁵` (exact null), `2.1·10⁻¹¹` (near-null), then `1.4·10⁻⁷`, `0.15` |
| `dν = 0.04` | `1.0·10⁻¹⁵, 3.6·10⁻¹⁵, 2.1·10⁻¹¹`, then `1.2·10⁻⁴` |

Profiles of the three directions:

| direction | `|Y|²` quantiles (10 / 50 / 90 %) in `ν` | response weight `|χ|Σ|R_j|` in `ν ∈ [2.5, 18]` |
|---|---|---|
| null-0 | 8.5 / 10.6 / 13.0 | 98.4 % |
| null-1 | 4.8 / 6.6 / 8.5 | 98.3 % |
| near-null | 4.3 / 11.6 / 15.4 (peak at 13.55, next to the zero 13.094 of `N̂₄`) | 93.9 % |

* **Analytic remark [D].** A compactly supported `χ` cannot be an exact continuum null vector: `W ≡ 0` forces `e^{±2t}`-decay of Fourier
  transforms, i.e. analyticity in a strip. The null directions are therefore rapidly decaying analytic functions represented on the grid.
  This is consistent with their localisation, and their large-`t` realisation is limited by the same truncation as in §4.

### 1.2 Every gauge difference is a null combination [DV]
`Y_g − Y_H0`, projected on the three directions (`dν = 0.02`), leaves a relative residual outside the span of:

| gauge | value piece | derivative piece | dominant directions |
|---|---|---|---|
| H2 | `2.0·10⁻⁹` | `1.4·10⁻⁸` | uses the near-null direction |
| H3 | `3.3·10⁻⁸` | `2.9·10⁻⁷` | essentially the 2 exact directions only |
| H4 | `4·10⁻¹¹` | `2·10⁻¹¹` | essentially the 2 exact directions only (outside them: `10⁻⁵`) |
| SB | `4·10⁻¹¹` | `6·10⁻¹²` | uses the near-null direction |

H1 differs through the next (`1.4·10⁻⁷`) direction; Stage 6 already found it not grid-stable.

**Consequence.** For any state, the gauge spread of the de-shifted image is `Σ_k c_k^{(g)} P_k[F_X]`, where
`P_k[F] = Σ_j∫χ_null,k R_j F_j` are three linear functionals of the exact raw image on `ν ∈ [2.5, 18]`. "K" in Stage 6 is a two-structure
model of these projections.

## 2. Response of the heavy action to beyond-LR corrections, and tolerances

### 2.1 Definition and amplification [DV]
* Write the per-harmonic forward image as `F_j(ν) = 2[P(ν−ij) + Γ_j(ν)/m²]`, with `Γ` in c₂ units (`m²` × relative correction, as Result D).
* Then `δ(∂B₄) = 2u³ Re Σ_j∫χ_d R_j Γ_j` and `δB₄ = (2u³/Λ²) Re Σ_j∫χ_v R_j Γ_j`. The value piece is `Λ⁻²`-suppressed.
* Stage 6's (a) is `Γ_j = K[ν² − (ν−ij)²]`; (b) is `Γ_j = K(j/3)²ν²`.
* Amplification `A[w] = ∫|χ|Σ_j|R_j| w(ν)`, derivative piece (`dν = 0.02`):

| gauge | `A[1]` | share at `ν <2.5 / 2.5–8.3 / 8.3–18 / >18` | `A[|ν−ij|²]` | spread vs H0 `∫|Δχ|Σ|R_j|` |
|---|---|---|---|---|
| H0 | 0.1714 | 72 / 24 / 4 / 0.1 % | 3.19 | — |
| H2 | 0.1633 | 75 / 24 / 1 / 0.1 % | 2.14 | 0.0113 (97 % in 2.5–18) |
| H3 | 0.1715 | 72 / 24 / 4 / 0.1 % | 3.21 | 0.0002 |
| H4 | 0.1616 | 76 / 24 / 0 / 0 % | 1.76 | 0.0185 (98 % in 2.5–18) |
| SB | 0.2021 | 61 / 20 / 18 / 1.3 % | 9.49 | 0.0664 (94 % in 2.5–18) |
| XP | 0.1732 | — | 4.04 | (other grid) |

### 2.2 Tolerances at the binding states (HYB model of Stage 6) [M]
* **Uniform bound.** `G_U` is the largest uniform `|Γ_j(ν)| ≤ G` (worst sign pattern) keeping `FF ≥ 0`.
* **ν²-growth bound.** `G_Q` is the same for `|Γ_j| ≤ g|ν−ij|²`.
* **Spread bound.** `G_S` is the uniform `|Γ|` that the gauge spread alone needs to flip the state.
* Shown: Δ_g = 300. Δ_g = 100 and 1000 agree to ≤ 2 % except at the `+−` gap states, whose `FF` depends on Δ_g.
* Columns are min … max over H0, H2, H3, H4, SB; the minimum (SB) is listed first.

| state (worst `m`) | `FF/N` | `G_U` δ = 0.05 | `G_U` δ = 0.2 | `G_Q` δ = 0.05 / 0.2 | `G_S` (H4, SB) δ = 0.05 / 0.2 |
|---|---|---|---|---|---|
| `++` J = 0, gap (`m = Λ`) | +0.050 / +0.200 | **2.9 … 3.6** | **8.8 … 11.0** | 0.06–0.33 / 0.19–1.01 | (31, 8.7) / (96, 27) |
| `++` J = 2, `1.05–1.1Λ` | +0.050 / +0.221 | 4.2 … 5.2 | 13.9 … 17.4 | 0.09–0.48 / 0.30–1.60 | (46, 12.7) / (152, 42) |
| `++` J = 4, `1.4–1.45Λ` | +0.31 / +0.50 | 56 … 70 | 74 … 93 | 1.2–6.4 / 1.6–8.6 | (609, 170) / (813, 226) |
| `+−` J = 4, gap `1.013Λ` | +0.079 / +0.228 | 5.2 … 6.5 (Δ_g = 1000: 3.7 … 4.7) | 11.8 … 14.8 | 0.11–0.60 / 0.25–1.36 | (57, 16) / (129, 36) |
| `+−` J = 5, `1.15–1.2Λ` | +0.050 / +0.203 | 8.8 … 11.1 | 28 … 35 | 0.19–1.02 / 0.60–3.26 | (97, 27) / (310, 86) |
| `+−` J = 6, `1.02–1.5Λ` | +0.051 / +0.208 | 10 … 12.6 | 47 … 59 | 0.21–1.16 / 1.0–5.4 | (110, 31) / (515, 143) |
| `+−` J = 7, 8 | ≥ +0.18 / ≥ +0.37 (Δ_g = 1000: +0.088 / +0.29 at J = 8, `1.008Λ`) | ≥ 51 (all Δ_g) | ≥ 112 | ≥ 1.1 / ≥ 2.4 | ≥ (561, 156) |

Reading:
* **Binding state.** The `++`, J = 0 gap state binds every tolerance (`u³ = 1`).
* **Single gauge.** Positivity there needs the EE beyond-LR correction of the `k = 4` forward image, relative to the covariant HYB model, to be
  `|Γ| ≲ 3` c₂ (δ = 0.05) or `≲ 10` c₂ (δ = 0.2), dominated by `ν < 2.5`.
* **Gauge spread.** The spread alone matters only if `|Γ| ≳ 9–31` (δ = 0.05) or `27–96` (δ = 0.2) on `ν ∈ [2.5, 18]`.
* **`c₄` scale.** The covariant finite-J constant HYB already contains is `γ′(0) = (c₄ + ½)/2`, i.e. `Γ ≈ 0.75 μ²` (`++`) and `20.75 μ²` (`+−`).
  * A non-covariant part equal to a fraction ε of it, with that `μ²` growth, is tolerated at the J = 0 gap state up to
    `ε ≈ 0.08–0.44` (δ = 0.05) or `0.25–1.35` (δ = 0.2).
  * At the `+−` J = 4 gap state the same fraction is tolerated only up to `ε ≈ 0.4–3 %` (δ = 0.05) or `1–7 %` (δ = 0.2), because
    `c₄^{+−} = 41`.
  * This framing assumes that the non-covariant part scales with `c₄_h`, which is not known [M].

## 3. Why `K` was not measured, what exists, and natural-strength scans

### 3.1 Reach and precision of the existing machinery (census from EE_DG300_SMALLNU, HYPOTHESES_CLOSE, NEARGAP100) [NE]
* **Fold.** The crystal saddle folds at **`ν* ≈ 11.3–11.5`**, independently of `J` and `Δ` (EE_DG300 §1).
  * The `K₂/K₃` Laplace rows are expansions about that saddle, so they do not exist below the fold.
  * Even above it, `K₂` had the wrong sign at `x = 0.05` for J = 24.
* **n3 contour (`ν_f = 13`).** Scalar-validated to `|Δc₂| ≲ 100` (`ν ≥ 11.7`), `≲ 200` (`ν = 10`) and `±1000` (`ν = 8.3`), at
  `(336, 24/28)`, `2¹¹·2·8` points.
* **n4 thimbles (`ν = 5–10`).** Errors from `±51` (`ν = 9.2`, J = 4, `m ≈ 100`) to `±250` with a 4.7σ phase systematic (`ν = 5.45`), and
  `±3500 … ±2·10⁴` at J = 24/28.
  * Monte-Carlo errors in c₂ scale like `m²` × relative error, so at `m ≈ 300` the J = 4 figures are ×9.
  * **No method below `ν = 5`.**
* **Comparison.** The J = 0 requirement (§2.2) is `≲ 3–11` c₂ with 72 % of the weight at `ν < 2.5`. That is ≳ 10–1000× beyond the existing
  precision, in a window that no method reaches.
* **The `k = 4` integrand.** By Lemma B it is the `k = 2` integrand times the `ν`-independent, O(1) prefactor `R(u′, v′)`. That changes the
  amplitudes, not the saddle/fold structure of the fast exponent `(iν/2)log A + (Δ/2)log B`, so the reach above is unchanged [D-sketch].
  Building it would not have produced a usable `K`, so it was not built.
* **Position-space exact moments** (directgrav). They reach `ν ≲ 5` at `m ≈ 100`, but:
  * they compute Regge moments `Π̂`, not the `Ψ₄` functional image, and the `Π̂ ↔ Ψ_k` relation for spinning blocks is [I] (SPINNING_MELLIN §5);
  * `+−` is noise-limited;
  * J = 0, 2 were never run (2.4 h/state).
* **Validation tasks.** They need an engine and were therefore not executable. What does hold:
  * the de-shift maps every covariant image exactly (Stage 6; re-checked on the XP weights, §4.2);
  * the scalar CLS-D `k = 2` image is exactly covariant by the magic-block identity `N(ν)[g C(ν+i) + g̃ C(ν−i)]` (H0 §24.1), so the scalar
    skeleton is reproduced in every gauge;
  * the `k = 4` analogue of that identity is assumed in EXA [I].

### 3.2 Indirect size indications (proxies, not the object) [NE]
* **Result D's direct crystal values** (`k = 2`, Δ = 400 rows):
  * `++` J = 4: `Γ_c = −18 / +5 / +4` at `ν = 10 / 20 / 40`;
  * `+−` J = 5: C2 `−13 / +2 / +2`, C3H `+119 / +15 / +16` (HYPOTHESES_CLOSE §4.4).
* **directgrav `++` J = 4, 6** (Δ_g = 100): graviton Regge moments equal `c*` × the scalar skeleton, with shape deviation `≤ 3·10⁻³` at
  `ν ≲ 3.5–4.4`, i.e. `|Γ| ≲ 30–50` c₂.
* **Reading.**
  * These proxies are comparable to or larger than the J = 0 / J = 2 tolerances (3–17 c₂), and well below the `++` J = 4 ones (56–93).
  * There are no data at all for J = 0, 2.
  * The proxies neither rule out nor indicate a failure at the binding states.

### 3.3 Minima under Stage 6's structures at natural `|K|` (`stage7_kscan.py` → `cache/stage7_kscan.txt`) [M]
Min `FF/N` over Stage 6's 14 000 states, worst sign, Δ_g = 300. Δ_g = 100 / 1000 agree to ≤ 0.002.

| gauge | δ = 0.05, (a), `|K|` = 0.5 / 1 / 2 | δ = 0.05, (b) | δ = 0.2, (a) | δ = 0.2, (b) |
|---|---|---|---|---|
| H0 | +0.011 / **−0.028** / −0.106 | +0.034 / +0.018 / −0.014 | +0.149 / +0.098 / −0.004 | +0.179 / +0.158 / +0.117 |
| H2 | +0.014 / −0.021 / −0.093 | +0.048 / +0.047 / +0.044 | +0.153 / +0.107 / +0.014 | +0.198 / +0.196 / +0.192 |
| H3 | +0.011 / −0.028 / −0.106 | +0.033 / +0.017 / −0.017 | +0.149 / +0.098 / −0.004 | +0.178 / +0.157 / +0.113 |
| H4 | +0.012 / −0.027 / −0.103 | +0.050 / +0.050 / +0.049 | +0.150 / +0.100 / −0.000 | +0.198 / +0.198 / +0.198 |
| SB | +0.050 (all) | −0.000 / **−0.051** / −0.152 | +0.198 (all) | +0.134 / **+0.069** / −0.063 |
| XP (Stage 7) | +0.025 / +0.000 / −0.050 | +0.043 / +0.036 / +0.021 | +0.168 / +0.135 / +0.070 | +0.191 / +0.181 / +0.162 |

* Derivative stress coefficients `(s_a, s_b)`:

  | H0 | H2 | H3 | H4 | SB | XP |
  |---|---|---|---|---|---|
  | (−0.90, +0.37) | (−0.83, +0.04) | (−0.90, +0.38) | (−0.89, −0.01) | (0, −1.17) | (−0.58, −0.17) |

* `K_tol` (both structures, worst sign):
  * δ = 0.05: 0.50 (SB) … 0.64 (H0) … ≈ 1.0 (XP);
  * δ = 0.2: 1.52 (SB) … 1.96 (H0/H3) … ≈ 3.1 (XP).
* The binding state is always `++`, J = 0 at the gap.
* The new XP gauge is the least sensitive to (a), but only because of how the `ν²`-weighted moment depends on `ν_max`. That is a further sign
  that the `K` coefficients are not intrinsic.

## 4. Extended-precision de-shift weights (task 3) (`stage7_xprec.py`, `stage7_xpcheck.py`, `stage7_s5check.py`)

### 4.1 Method and diagnosis [DV + NE]
* **Method.** Iterative refinement of the same discrete Abel systems:
  * residuals in exact 384-bit fixed-point integer arithmetic (the matrices are DCT/DST-type with a `4N_t`-entry cosine table);
  * corrections from the double-precision pseudo-inverse (tol 1e-11, i.e. the stored gauge);
  * `W(t)` then evaluated **off** the collocation grid in mpmath.
* **On the stored grid** (`dν = 1/50`, `ν_max = 41`), derivative piece:
  * the on-grid residual falls from `2.6·10⁻¹⁶` to `2.2·10⁻²⁶`, floored by the dropped (exact-null) left directions at `5·10⁻²⁶`;
  * the off-grid `|W/Φ − 1|` is essentially unchanged:

    | t | 10.1 | 12.6 | 15.3 |
    |---|---|---|---|
    | 384-bit solve | `5·10⁻⁷` | `3·10⁻⁴` | 0.6 |
    | stored weight | `2·10⁻⁶` (t = 10) | `10⁻³` (t = 12) | |

    The value piece behaves the same way.
* **Diagnosis.** The large-`t` failure is **truncation**: `|Y(ν_max = 41)| ≈ 10⁻¹⁴`. It is not round-off. Stage 6 §2.4's attribution
  (`ε cosh jt` with ε = double precision) is corrected here.
* **Grid tests** (derivative piece):

  | grid | `|W/Φ − 1|` off-grid |
  |---|---|
  | `dν = 1/100`, `ν_max = 41` | worse (`10⁻⁵` at t = 10) |
  | `dν = 1/50`, `ν_max = 82` | `5·10⁻⁷` at t = 25, 0.6 at t = 30 |
  | **`dν = 1/25`, `ν_max = 123`** | **`4.4·10⁻¹⁸` at t = 30**, `1.7·10⁻⁵` at t = 40 |

### 4.2 The kept C4 weights (XP; `functionals_q7_xprec.json`) [DV]
| piece | off-grid `|W/Φ − 1|` at t = 20 / 25 / 30 / 35 | `∫χN̂` | shifted moments `M₀…M₃` | covariant test | gauge change vs stored |
|---|---|---|---|---|---|
| value | `4·10⁻²⁸ / 6·10⁻²² / 3.4·10⁻¹⁵ / 8·10⁻⁹` | 1.000000000000 | `1, 7·10⁻¹⁶, 10⁻¹³, 2·10⁻¹¹` (stored: M₃ = −1.1·10⁻⁵) | `2.4·10⁻¹⁵` | `1.4·10⁻²` |
| derivative | `10⁻³⁰ / 10⁻²⁴ / 4.4·10⁻¹⁸ / 8·10⁻¹²` | 0 | `0, 1, 8·10⁻¹⁵, 9·10⁻¹³` | `8·10⁻¹¹` (finite-difference limited) | `3.5·10⁻³` |

* The covariant test uses the Δ_g = 300 states, since `g` has a pole at `μ = m`, which must exceed `ν_max`.
* XP is the min-norm gauge on the extended grid. It is a valid gauge, with the same light side and the same covariant images.
* The near-gap minima of Stages 3–6 are unchanged in it (covariant models). Its `K`-sensitivity is in §3.3.

### 4.3 C5: not achieved [NE]
* **Runs** (`dν = 1/25`, `ν_max = 82`; general complex gauge, 4096 unknowns): the refinement stalls at the residual in the dropped space.

  | correction tol | residual floor | `|W − 1|` at t = 15.3 | t = 17.9 | t = 20.2 | t = 22.7 |
  |---|---|---|---|---|---|
  | 1e-11 | `10⁻¹⁴` | | | | |
  | 1e-14 | `6·10⁻¹⁸` | | | | |
  | **1e-15 (best)** | `3·10⁻¹⁹` | `4.9·10⁻⁷` | `7.9·10⁻⁵` | `8.4·10⁻⁴` | 1.1 |

* **Why.** C5's spectrum has no gap (Stage 5: about `N − 25` null-ish directions). The large-`t` cancellation needs singular directions far
  below double precision.
* **What would fix it [I]:** an extended-precision SVD, or an analytic construction from the shift equation
  `Σ_j (cĥ_j)(ν + ij) = 0`, whose rational coefficients admit Γ-function solutions.

### 4.4 Stage 5's large-β positivity with the realised images (`cache/stage7_s5check.txt`) [NE + M]
* **Method.** Stage 5's scan (same states, β grid, `Nrm`), with `P₀`, `D` and the C5 image replaced by the Mehler sums weighted by
  `(1 + e(t))`, where `e = W/Φ − 1` of the actual weight.

| weights | min `FF/Nrm`, `b̂ ≥ 20` | β ≥ 6 | first negative |
|---|---|---|---|
| exact (Stage 5) | +0.248 / +0.065 / +0.191 (δ = 0.05); +0.297 / +0.226 / +0.357 (δ = 0.2) | ≥ +0.987 | none |
| **stored** | −1.000 (β = 29) | −1.000 | **β = 15.6**, `++`, J = 3.9·10⁵ … 1.4·10⁷ (all Λ, both δ) |
| **Stage 7** (C4 XP + best C5) | as exact | **+0.9645 / +0.9953 / +0.9992** (Λ = 100 / 300 / 1000) | none to β = 30 |

* **Consequences.**
  * Stage 5's β ≤ 30 positivity holds for the **continuum** functionals and, numerically, for the **Stage-7** weights.
  * The stored double weights, taken literally as functionals, are **not** heavy-positive beyond β ≈ 15.4 (first failure at J ≈ 4·10⁵,
    `m ≈ 3Λ`, Λ = 100).
  * They should be replaced by `functionals_q7_xprec.json` for C4. For C5, it should be stated that the functional is defined through its
    continuum Abel equation.
* **Margin at large β.**
  * At Λ = 100 the remaining margin at β ≈ 30 (+0.9645 vs +0.99) is eaten by the unrealised C5 tail, whose error grows like `e^{2t}`.
  * Beyond β ≈ 31–32 the Stage-7 C5 weight would fail as well [NE, extrapolated].

## 5. Verdict on δ (task 2)
* **In-model (covariant, `K = 0`)** [M, Stage 6, unchanged]: both δ positive.
  * δ = 0.05: min `+0.0478 / +0.0497 / +0.0491`; δ = 0.2: `+0.1985 / +0.1977 / +0.1970`.
  * The binding near-gap minima are as Stage 6. No "measured correction" exists to add.
* **Every reasonable gauge** [M, conditional on the stress model]:
  * δ = 0.2 holds for `|K| ≤ 1` (min **+0.068**) and fails by `|K| ≈ 1.5–2` (SB, H0/H3/H4);
  * δ = 0.05 fails at `|K| = 1` in every gauge but XP (min −0.051), and holds only for `|K| ≲ 0.5`.
  * Equivalently, uniform `|Γ| ≲ 3` c₂ (δ = 0.05) or `≲ 10` c₂ (δ = 0.2) at the J = 0 gap state (§2.2).
* **Unconditional: [I].** The decision needs the EE `k = 4` image at `ν ≲ 18` (mostly `< 2.5`) to ~1–10 c₂.
* **Recommendation unchanged:** quote δ = 0.2 (`A = 17.241`, `B_eff = −19.038 / −19.220 / −19.265`, `κ = 0.110`). Realise the functional
  with the Stage-7 C4 weights.

## 6. Remaining open items [I]
1. **`K`, equivalently the 3 null projections plus the absolute correction:** the exact EE `k = 4` image beyond LR at `ν ≲ 18` for the
   `++` J = 0, 2 and `+−` J = 4 gap states, to ≲ 3 c₂.
   * A new method is needed below the fold. Two candidates:
     * hard-direction Laplace × exact soft-block quadrature, with the `1/m²` coefficient computed directly rather than as a Monte-Carlo
       difference;
     * a spinning magic-block (`Ψ ↔ C`) identity, which would make `K ≡ 0` structurally.
2. C5 weights realised beyond t ≈ 18 (§4.3). Note: β > 30 is not covered by any stage.
3. Whether the `Ĉ₄` null directions survive in the continuum as analytic decaying functions (§1.1 remark). If they do not, the gauge
   question reduces to the absolute correction alone.
4. All items inherited from Stage 6 §5 (flat LP rows at J > 2500, exact skeleton at J ≥ 40, `k = 3` remainder, parity-odd beyond LO, `b₀`,
   `b₅`, `c₁`, (F), (P), A1, A2, K1′/K2/K3, KKSDZ loophole, odd-`k` continuation).

## 7. Graded claims

| claim | grade |
|---|---|
| `Ĉ₄` gauge space = 2 exact + 1 near-null directions, ≥ 94 % of their response weight at `ν ∈ [2.5, 18]`; all Stage-6 gauge differences in their span (`≤ 3·10⁻⁷`) | [DV] |
| gauge spread = 3 null projections of the exact image; absolute correction dominated by `ν < 2.5` (72 %, H0 derivative) | [DV] |
| tolerances `G_U ≈ 3 / 10` c₂ (δ = 0.05 / 0.2) at the `++` J = 0 gap state; table §2.2 | [DV] (response) + [M] (HYB `FF`) |
| no existing engine reaches `ν < 5`, nor `< 10³` c₂ precision at `ν = 5–8`; `K` not measured | [NE] (census) / [I] |
| `k = 4` integrand (`× R`) does not change the reach | [D-sketch] |
| δ = 0.2 positive in every gauge for `|K| ≤ 1` (min +0.068); δ = 0.05 not (min −0.051) | [M] (stress model) |
| stored weights' large-`t` failure is `ν_max`-truncation, not round-off (corrects Stage 6 §2.4) | [NE] |
| XP C4 weights: `W = Φ` to `≤ 4·10⁻¹⁵` at t = 30, valid gauge (moments, `R⁴` normalisation, covariant images) | [DV] |
| C5 extended solve stalls (no spectral gap): realised only to t ≈ 18 | [NE] |
| stored weights, realised: Stage-5 positivity lost from β = 15.6; Stage-7 weights: positive to β = 30 (≥ +0.9645) | [NE + M] |
| δ = 0.2 recommended; constants unchanged | [M] + [D given Q]; robustness [I] |

## 8. Files (`ads4/code/resultq/`)

| file | content | run time |
|---|---|---|
| `stage7_null.py` | null space (two grids), gauge-difference decomposition, response amplifications by `ν` window, per-state tolerances → `cache/stage7_null.{txt,pkl}` | 5 s |
| `stage7_kscan.py` | min `FF/N` under Stage-6 structures (a), (b) at `|K|` = 0.25 … 5, per gauge incl. XP → `cache/stage7_kscan.{txt,pkl}` | 5 s |
| `stage7_xprec.py` | 384-bit iterative refinement of the C4 / C5 Abel systems (env `XP_DNU`, `XP_N`, `XP_TOL`), off-grid `W` → `cache/stage7_xprec_*.pkl`, logs `cache/stage7_xprec_t*.log`; `save` writes **`functionals_q7_xprec.json`** | 25–60 s per piece |
| `stage7_xpcheck.py` | validity of the XP C4 weights (moments, `∫χN̂`, covariant test, amplification, `s_a`, `s_b`) → `cache/stage7_xpcheck.txt` | 3 s |
| `stage7_s5check.py` | realised `e(t)` curves (stored vs Stage 7) and Stage 5's scan with realised images → `cache/stage7_s5check.txt`, `cache/stage7_s5check_curves.pkl` | 3 min |
| `functionals_q7_xprec.json` | C4 value / derivative (`dν = 1/25`, N = 3072) and C5 (`dν = 1/25`, N = 2048) weights as 64-digit strings, with off-grid accuracy and gauge-change records | — |
