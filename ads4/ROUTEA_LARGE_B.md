# H13′ for the route-A′ D-safe functional at large impact parameter

Status: 18 Nov 2026, H13′ agent. **Unrefereed. Not committed.**

**Functional.** Route A′, `Ω^{A′} = ∫ψ₂Ĉ₂ + ∫ψ₃Ĉ₃`, with the NF02 weights and stripped images (`H0_SPINNING_FUNCTIONAL.md` §10). It
replaces the withdrawn B′-uplift `Ω^{NF,SC}`; `NF02_LARGE_B.md` is about that withdrawn functional, apart from its functional-independent
parts (its §7).

**Model [M].** Identical to `spinfunc/f6_routeA2.py`. For each state:
* **flat part:** the stripped `B₂ + B₃` (this is NF02);
* **curvature:** `curv = 2∫ψ₂[P_{−½+iΛp}(cosh β) − J₀(b̂p)] + uσ∫ψ₃[P − J₀]`;
* **remainder:** `S = ∫ψ₂E₂s_hρ₂ + uσ∫ψ₃E₃s_hρ₃`, where `ρ_k` is the exact scalar-skeleton remainder (`k = 3` interpolated). The
  stripped envelope uses `s_h = (1−y)^{±2}`; the "B′ env." variant uses `s_h = 1`;
* **normalisation:** `N_s = |B₂+curv₂| + |B₃+curv₃|`.

`F_hyb := B₂ + B₃ + curv` is computed with exact Mehler–Dirichlet and `J₀` transforms. It reproduces `f6`'s formula to `10⁻¹³`. The
exact action is `F = F_hyb + S`.

**Code** (new files only): `ads4/code/exactskel/h13a_common.py`, `h13a_hyb.py`, `h13a_rem.py`, `h13a_summary.py`, `h13a_cert_curv.py`.
They reuse `h13_common.py`, `h13_cert_m1.py`, and the exact `ρ_k` of `h13_exact.py` (`cache/h13_ex_<Δg>.pkl`). Outputs are in
`exactskel/cache/` (git-ignored): `h13a_summary.txt`, `h13a_rem.txt`, `h13a_rem_table.md`, `h13a_cert_curv.txt`, `h13_cert_m1.txt`.

Grades:
* **[certified]**: Arb or exact, with rigorous enclosures;
* **[D]**: derived;
* **[NE]**: numerics without rigorous truncation bounds;
* **[M]**: model;
* **[I]**: inferred.

## 0. Results

**Main result.** On the whole region asked for, no route-A′ state is non-positive, and the minimum `F/N_s` is about 0.93.
* At large `b̂` the action is dominated by a positive **AdS-only** term `A₁M₁(β)/Λ²`, not by the flat `1/b̂³` tail.
* Its positivity, together with an explicit bound on the remainder of the envelope, is **certified** on a `b̂`-window per `Δ_g`.
* What remains is assumed or numerical: the scalar-skeleton remainder `S` beyond the exact sample, H2/H4 (the transfer to gravitons),
  and H12′.

| item | result | grade |
|---|---|---|
| **Hybrid scan**: every `J`, `m ∈ [1.5, 10]Λ` (+ `20, 50, 100Λ`), `b̂ ∈ [20, 200]` (`Δ_g` = 300, 1000) or `[40, 200]` (100); spot checks `b̂ ∈ [200, 1000]` at 5 masses; both channels | no non-positive state. Min `F_hyb/N_s` = **0.9313** (300) and **0.9311** (1000), both at `+−`, `1.5Λ`, `J = 15`, `b̂ = 20.7` (the corner of the region); **0.9411** (100: `+−`, `1.71Λ`, `J = 95`, `b̂ = 112`). On `b̂ ≤ 1000`: ≥ 0.899. The AdS curvature is **positive on every record** | [NE] (exact ingredients) |
| **Exact skeleton remainder**: 59 states × 2 skeletons (`J ≤ 150`, `m = 1.5–10Λ`) | no non-positive reliable record. Min `F/N_s` = **0.9489** (300), **0.9492** (1000), **0.9669** (100). `|S|/F_hyb ≤ 6.1·10⁻⁴ / 5.4·10⁻⁵ / 3.8·10⁻³`, and `|S|/F_hyb ≤ 122/m²` everywhere. Reliable records reach `b̂` = 150 (1000), 100.5 (300), 70.5 (100) | [NE + M] |
| **Large-`b̂` form** | `F_hyb ≃ A₁M₁(β)/Λ² − A₂/(2b̂³)`, with `A₁ = 377.79 + 50.64σu ≥ 355.3` and `A₂ = −2047.09`. The **AdS-only** term `A₁M₁/Λ²` dominates the flat tail once `b̂ ≳ 4Λ^{2/3}` (85 / 178 / 398). Matches `F_hyb` to `[0.9998, 1.003]` for `b̂ ≥ 100` and `[0.9975, 1.012]` for `b̂ ≥ 50` on every record | form, coefficients, identities [D]; match [NE] |
| **`M₁(β) > 0` for all `β > 0`** | Arb plus an analytic large-`β` bound (`h13_cert_m1`) | **[certified]** |
| **Curvature envelope**: `|curv − A₁M₁(β)/Λ²| ≤ R(b̂)`, explicit (§4.2), every `m ≥ 1.5Λ`, both channels | `curv ≥ A₁M₁/Λ² − R > 0` on **`b̂ ∈ [115.5, 197]` (Δ_g = 100), `[162.6, 597]` (300), `[274.1, 1997]` (1000)** | **[certified]**, given the derived bound |
| ⇒ **`F_hyb > 0`** on those windows | `F_hyb = F_flat + curv`. `F_flat ≥ (3/20)N` is certified (`NOFORWARD_FUNCTIONAL.md` §4) except on `+−`, `J ≥ 400`, `1.5Λ ≤ m < 1.633Λ`. That corner has `b̂ ∈ [490, b̂_max(1.633Λ)]`, with `b̂_max ≈ 1.02Λ`. It is empty at `Δ_g = 100, 300`. At `Δ_g = 1000` it sits inside the window (`b̂ ≈ 490–1020`), where the flat positivity is float-only (referee: min `0.765`, `J ≤ 20000`) | **[certified]** at 100, 300; at 1000 certified except that corner, where it is [NE] (flat) + [certified] (`curv > 0`) |

**Proposed H13′ entry for route A′:**
* [NE] for the hybrid action on the whole region (every `J`);
* [NE] for the exact action on the sample (`b̂ ≤ 150`, at most 70 at `Δ_g = 100`);
* **[certified] positivity of the hybrid action** on `b̂ ∈ [115, 197] / [163, 597] / [274, 1997]`;
* [I] for `|S| < F` off the sample. It is calibrated by `|S|/F_hyb ≤ 122/m²`, a margin of at least 150.

**Not covered rigorously:**
* `b̂ ∈ [20, b̂₀)` (scanned only, [NE]);
* `b̂ > 2Λ` (scan and asymptotic form only).

## 1. What carries over, and what is new

* **Carries over from `NF02_LARGE_B.md`:**
  * conventions (`b̂ = (2J+1)Λ/m`, `cosh β − 1 = b̂²/(2Λ²)` exactly, the twist-gap cap on `b̂`);
  * the exact scalar-skeleton remainders `ρ_k`, which are functional-independent;
  * the method fixes (`mack_mat`, Arb tail weights, and tail order `n_q ≥ 960`).
* **What route A′ changes:**
  * Only `k = 2, 3`, with weights `ψ₂`, `ψ₃` instead of `φ_k`.
  * The stripped remainder envelope.
  * The leading large-`b̂` coefficients `A₁`, `A₂` are **the same** as for the withdrawn functional, because `φ_k′(0)` and `φ_k″(0)` vanish
    for `k ≥ 4` and `φ₃ − ψ₃ = O(p³)`. From `A₃` on they differ (route A′: `A₃ = 4516.3 + 6704.7σu`).

## 2. Hybrid scan (`h13a_hyb.py`, `h13a_summary.py`) [NE]

* About 27 000 (`Δ_g = 300, 1000`) and 18 000 (`Δ_g = 100`) state-helicity records in `b̂ ∈ [20, 200]` (`[40, 200]`), every `J`.
  Plus 10 000–13 000 records per `Δ_g` in `b̂ ∈ [200, 1000]`.
* No non-positive state. The curvature `curv = F_hyb − F_flat` is positive on every record (min `curv/F_hyb = +7·10⁻⁵`, at `b̂ ≈ 20`,
  `Δ_g = 1000`).
* `+−`, `Δ_g = 300`, min `F_hyb/N_s` by band:

  | `m/Λ` \ `b̂` | [20,40) | [40,70) | [70,120) | [120,200] |
  |---|---|---|---|---|
  | [1.5, 2) | 0.931 | 0.982 | 0.974 | 0.935 |
  | [2, 3) | 0.966 | 0.990 | 0.985 | 0.963 |
  | [3, 5) | 0.985 | 0.996 | 0.993 | 0.983 |
  | [5, 10] | 0.995 | 0.999 | 0.998 | 0.994 |

  `++` is ≥ 0.988 everywhere. `Δ_g = 1000` is the same to `10⁻³`. Masses `20–100Λ`: ≥ 0.9996.
* **Refinements.** In 0.005 steps on `m ∈ [1.5, 1.6]Λ` (`Δ_g = 300, 1000`) and 0.01 steps on `[1.6, 1.8]Λ` (100). The minima stay at the
  corner (`m = 1.5Λ`, smallest admissible `b̂`) for 300 and 1000, and at `1.71Λ`, `b̂ ≈ 112` for 100. Below `b̂ = 20` they join the
  already-scanned near-gap region, whose minimum is 0.21–0.23, `+−` `J = 5`.

## 3. Exact scalar-skeleton remainder (`h13a_rem.py`) [NE + M]

* The `ρ_k` come from `h13_exact.py`: Arb Laplace/Barnes, exact Mack polynomials, three `n`-sum configurations, odd `k = 3` by quartic
  interpolation in `k` [M].
* The route-A′ `S` is recomputed from them, for both envelopes.
* **Reliability:** the spread of `S` across the `n`-sum configurations must be ≤ `10⁻³ N_s` (route A′ and old functional).

__REMTABLE__

**Envelope variant.** The B′-envelope variant of `S` is smaller (`|S|/F_hyb ≤ 2.6·10⁻³` at `Δ_g = 100`). The headline uses the stripped
envelope.

**Odd-`k` interpolation sensitivity:** `≤ 2·10⁻⁵ N_s`.

**Pointwise vs integrated remainder.** As for the old functional, the pointwise remainder is up to `3·10² N_s`, but its integral is
`≤ 4·10⁻³ N_s`. The remainder integrates with the same forward-free cancellation as the leading action.

## 4. Large-`b̂` form and certificate

### 4.1 Asymptotic form [D]

**Exact representation.** Write `Φ_k(x) = ∫₀¹ψ_k cos(px)dp`. Mehler–Dirichlet for `P_{−½+iν}`, `J₀(z) = (2/π)∫₀^{π/2}cos(z sin θ)dθ`
and the substitution `s = 2Λ sinh(t/2)` give

`curv = Σ_k w_k (2/π)∫₀^{π/2}D_k(b̂ sin θ)dθ`, `D_k(s) = w(s)Φ_k(τ(s)) − Φ_k(s)`,

with `w₂ = 2`, `w₃ = uσ`, `τ = 2Λ asinh(s/2Λ)`, `w = (1+s²/4Λ²)^{−½}`.

**Leading terms.**
* The large-`Λ` endpoint expansion gives `F_hyb ≃ A₁M₁(β)/Λ² + A₂M₂(β)/(2Λ³) + …`, with `A_n = 2ψ₂^{(n)}(0) + σuψ₃^{(n)}(0)`.
* `M₂/Λ³ = −1/b̂³` exactly, so the second term is the flat tail.
* `M₁(β) = ∫₀^∞2νP_{−½+iν}(cosh β)dν/(e^{2πν}+1) = (√2/π)∫₀^βK(t)[cosh(t/2)/(4sinh²(t/2)) − 1/t²]dt`.
* `M₁` is also the `j = 0` term of the expansion of `D`:
  `(2/π)∫₀^{π/2}[1/s² − w/τ²]dθ = M₁/Λ²`. The proof: change variables to `t = τ/Λ`; the `1/ε` divergences of the AdS and flat kernels cancel
  because `K(0) = √2Λ/b̂`; the flat Hadamard finite part vanishes. It is checked to 15 digits.

**Which term dominates.** `A₁ ≥ 355 > 0` and `M₁ > 0`, so the AdS term is positive. It overtakes the flat tail `−A₂/(2b̂³) > 0` at
`b̂ ≈ 4Λ^{2/3}`. Beyond that, the D-safe action at large `b̂` is an `O(Λ⁻²)` AdS effect, set by the moments at `ν = O(1)`.

**Numerical check [NE].** `F_hyb/(two-term form)` lies in `[0.99983, 1.00294]` (`b̂ ≥ 100`) and `[0.9975, 1.0116]` (`b̂ ≥ 50`) on every
record, and in `[0.99999, 1.0007]` on `b̂ ∈ [200, 1000]`.

### 4.2 Explicit remainder bound and certified window (`h13a_cert_curv.py`) [D + certified]

**Splitting.** For `s ≤ b̂ ≤ 2Λ` (so `ρ := s²/(24Λ²) ≤ 1/6`), split `θ` at `s₁ = b̂ sin θ₁`. Then

`|curv − A₁M₁(β)/Λ²| ≤ R(b̂) = [C_I(1 − s₁²/b̂²)^{−½} + C_sm/Λ² + C_osc]/b̂`,

where:
* **`C_I = (2/π)Z₀`**, with `Z₀ ≥ Σ_k|w_k|max ∫₀^{s₁}|D_k + ψ_k′(0)e₀|ds` and `e₀ = w/τ² − 1/s²`.
  * Arb cells of width 1/25.
  * For `s ≥ 4`: second-order Taylor about the cell centre, using point values of `G` and `G′` and a derived bound on `G″`.
  * For `s < 4`: crude bound.
  * The bounds use the rigorous moments `∫p^n|ψ_k|` (`k = 2`: 10.4, 3.54, 1.52, 0.75; `k = 3`: 12.9, 5.13, 2.35, 1.19).
* **`C_sm = Σ_{k,j≥1}|w_k|max|ψ_k^{(2j+1)}(0)|κ_j(2/π)(π/2)^{2j}s₁^{1−2j}/(2j−1)`.** It comes from the exact endpoint expansion of
  `Φ_k`, with `|e_j| ≤ κ_j s^{−2j}/Λ²`, `κ_j = max(1/8, (1.2^{2j+2} − 1)/4)`, and `sin θ ≥ 2θ/π`.
* **`C_osc = Σ_{k,m≥5}|w_k|max|ψ_k^{(m)}(1)|(1.2^{m+1} + 1)(2/π)(π/2)^{m+1}s₁^{−m}/m`.** These are the `p = 1` endpoint terms; `ψ ∝ (1−p)⁵`.
* `|w₂| = 2`, and `|w₃| ≤ 4/9` for `m ≥ 1.5Λ`.

**Positivity condition.** `A₁,min M₁(β)/Λ² > R(b̂)` must hold on `b̂`-boxes (ratio 1.03).
* `M₁` is bounded below on each `β`-box by Arb lower sums.
* `R` is decreasing in `b̂`.
* `A₁,min = 355.28`.

| `Δ_g` | best `s₁` | `C_IΛ²` | `C_sm` | `C_oscΛ²` | **certified window** |
|---|---|---|---|---|---|
| 100 | 80 | 742 | 43.2 | 364 | **`b̂ ∈ [115.5, 197.0]`** |
| 300 | 120 | 1281 | 27.0 | 252 | **`b̂ ∈ [162.6, 597.0]`** |
| 1000 | 160 | 2693 | 20.0 | 497 | **`b̂ ∈ [274.1, 1997]`** |

**Hence [certified]: `curv > 0`, and therefore `F_hyb = F_flat + curv > 0`, on these windows, for every `m ≥ 1.5Λ` and both channels.**
* This holds wherever the flat certificate applies. At `Δ_g = 1000` the flat `+−`, `J ≥ 400`, `m < 1.633Λ` corner (`b̂ ≈ 490–1020`) is
  float-only.
* The derivation steps are [D]: the exact identities, `asinh x ≥ x − x³/6`, `(1+y)^{−½} ≥ 1 − y/2`, convexity, `sin θ ≥ 2θ/π`, and the
  mean-value and Taylor bounds.
* Every number is an Arb enclosure.
* **Limitation.** The bound is conservative. `C_I` contains no cancellation beyond `D ≈ −A₁e₀`, and `C_osc` drops the oscillation. So
  `b̂₀` sits well above where the numerics show the asymptotic form is already accurate (`b̂ ≈ 50`). Pushing `b̂₀` down to the scan would
  need a sharper bound in region I.

### 4.3 Coverage map for route A′ (`m ≥ 1.5Λ`)

| `b̂` | `Δ_g = 100` | `Δ_g = 300` | `Δ_g = 1000` |
|---|---|---|---|
| `< 20` (`< 40` at 100) | near-gap scans and certificate (`f6`, `CERTIFIED_NEARGAP` model) | same | same |
| `[20, b̂₀)` | hybrid [NE], every `J`; exact `S` [NE] to `b̂ = 70` | hybrid [NE]; exact `S` to 100.5 | hybrid [NE]; exact `S` to 150 |
| `[b̂₀, 2Λ]` | **hybrid [certified]** (`b̂₀ = 115.5`) | **[certified]** (162.6) | **[certified]** (274.1; flat corner [NE]) |
| `> 2Λ` | hybrid [NE] (scan to 1000), asymptotic form [D + NE] | same | same |
| exact `S` beyond the sample | [I], `|S|/F ≤ 122/m²` calibrated | same | same |

## 5. Caveats

* Everything is in the scalar-skeleton model [M]. H4 (transfer to gravitons with the stripped envelope), H1/H2 for `k = 2, 3`, and H12′
  are untouched.
* The `η`-inversion of `Ĉ₂` is not included (`f6`: `≤ 0.086 N_s` near the gap, smaller here). Its symbol is `O(1/Λ²)` relative, so it is
  far below the `≥ 0.93` margins, but it is not computed.
* **Exact remainder coverage.**
  * The Euler–Maclaurin `n`-sum fails near the twist gap: `Δ_g = 100`, `J = 60`, and `Δ_g = 300`, `J = 112`, both at `1.5Λ`.
  * `J = 200` was not run: `mack_mat` takes about 50 min there.
  * The `J = 150` states at `Δ_g = 300, 100` are __J150__.
* **The certificate is for the hybrid action.** Positivity of the exact action on the window still needs `|S| < F` there [I].

## 6. Files

| file | role |
|---|---|
| `h13a_common.py` | route-A′ hybrid action and remainder (as `f6`) |
| `h13a_hyb.py` | hybrid scan → `cache/h13a_hyb_<Dg>_<b0>_<b1>.pkl` |
| `h13a_rem.py` | route-A′ `S` from the cached exact `ρ_k` → `cache/h13a_rem.txt`, `h13a_rem_table.md` |
| `h13a_summary.py` | tables, envelope ratio → `cache/h13a_summary.txt` |
| `h13a_cert_curv.py` | derived remainder bound, Arb constants, certified windows → `cache/h13a_cert_curv.txt` |
| `h13_cert_m1.py` | `M₁ > 0` certificate (functional-independent) → `cache/h13_cert_m1.txt` |
| `h13_exact.py`, `h13_laplace.py`, `mack_mat.py` | exact `ρ_k` (shared with `NF02_LARGE_B.md`) |
