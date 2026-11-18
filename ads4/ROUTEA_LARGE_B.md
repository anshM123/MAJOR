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
| **Exact skeleton remainder**: 62 states × 2 skeletons (`J ≤ 150`, `m = 1.5–10Λ`) | no non-positive reliable record. Min `F/N_s` = **0.9489** (300), **0.9492** (1000), **0.9669** (100). `|S|/F_hyb ≤ 6.1·10⁻⁴ / 5.4·10⁻⁵ / 3.8·10⁻³`, and `|S|/F_hyb ≤ 122/m²` everywhere (`≤ 24/m²` at `b̂ ≥ 150`). Reliable records reach `b̂` = **200.7** (1000), 150.5 (300), 70.5 (100) | [NE + M] |
| **Large-`b̂` form** | `F_hyb ≃ A₁M₁(β)/Λ² − A₂/(2b̂³)`, with `A₁ = 377.79 + 50.64σu ≥ 355.3` and `A₂ = −2047.09`. The **AdS-only** term `A₁M₁/Λ²` dominates the flat tail once `b̂ ≳ 4Λ^{2/3}` (85 / 178 / 398). Matches `F_hyb` to `[0.9998, 1.003]` for `b̂ ≥ 100` and `[0.9975, 1.012]` for `b̂ ≥ 50` on every record | form, coefficients, identities [D]; match [NE] |
| **`M₁(β) > 0` for all `β > 0`** | Arb plus an analytic large-`β` bound (`h13_cert_m1`) | **[certified]** |
| **Curvature envelope**: `|curv − A₁M₁(β)/Λ²| ≤ R(b̂)`, explicit (§4.2), every `m ≥ 1.5Λ`, both channels | `curv ≥ A₁M₁/Λ² − R > 0` on **`b̂ ∈ [115.5, 197]` (Δ_g = 100), `[162.6, 597]` (300), `[274.1, 1997]` (1000)** | **[certified]**, given the derived bound |
| ⇒ **`F_hyb > 0`** on those windows | `F_hyb = F_flat + curv`. `F_flat ≥ (3/20)N` is certified (`NOFORWARD_FUNCTIONAL.md` §4) except on `+−`, `J ≥ 400`, `1.5Λ ≤ m < 1.633Λ`. That corner has `b̂ ∈ [490, b̂_max(1.633Λ)]`, with `b̂_max ≈ 1.02Λ`. It is empty at `Δ_g = 100, 300`. At `Δ_g = 1000` it sits inside the window (`b̂ ≈ 490–1020`), where the flat positivity is float-only (referee: min `0.765`, `J ≤ 20000`) | **[certified]** at 100, 300; at 1000 certified except that corner, where it is [NE] (flat) + [certified] (`curv > 0`) |

**Proposed H13′ entry for route A′:**
* [NE] for the hybrid action on the whole region (every `J`);
* [NE] for the exact action on the sample (`b̂ ≤ 200.7` at `Δ_g = 1000`, `≤ 150.5` at 300, `≤ 70.5` at 100);
* **[certified] positivity of the hybrid action** on `b̂ ∈ [91.2, 2Λ]` at all three `Δ_g` (§4.4; the first pass gave
  `[115, 197] / [163, 597] / [274, 1997]`);
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

**Minima over reliable records** (`cache/h13a_rem.txt`). No reliable record is non-positive.

| `Δ_g` | states (reliable helicity records) | min `F/N_s` (where) | max `|S|/F_hyb` (strip / B′ env.) | max `|S|m²/F_hyb` | `b̂` reached |
|---|---|---|---|---|---|
| 100 | 9 (32/36) | **0.9669** (`+−`, `1.5Λ`, `J = 46`, `b̂ = 62`) | `3.8·10⁻³` / `2.6·10⁻³` | 100 | 70.5 |
| 300 | 26 (100/104) | **0.9489** (`+−`, `1.5Λ`, `J = 16`, `b̂ = 22`) | `6.1·10⁻⁴` / `3.7·10⁻⁴` | 122 | 150.5 |
| 1000 | 27 (108/108) | **0.9492** (same state) | `5.4·10⁻⁵` / `3.3·10⁻⁵` | 122 | **200.7** |

* On the states with `b̂ ≥ 150` (`2Λ`, `J = 150` at 300 and 1000; `1.5Λ`, `J = 150` at 1000), `|S|m²/F_hyb ≤ 24`. The calibrated
  `122/m²` therefore does not grow with `b̂`.
* Per-state table: entries are the worst over the reliable skeletons; "—" means no reliable value.

**Δ_g = 100**

| `m/Λ` | `J` | `b̂` | `++` `F/N_s` | `+−` `F/N_s` | max `|S|/F_hyb` (strip) | max `|S|/F_hyb` (B′ env) | status |
|---|---|---|---|---|---|---|---|
| 1.50 | 30 | 40.7 | +1.0007 | +0.9814 | 3.8e-03 | 2.6e-03 | ok (Δφ=3,2) |
| 2.00 | 40 | 40.5 | +1.0008 | +0.9881 | 2.6e-03 | 1.8e-03 | ok (Δφ=3,2) |
| 1.50 | 46 | 62.0 | +1.0003 | +0.9669 | 1.1e-03 | 9.8e-04 | ok (Δφ=3,2) |
| 1.50 | 60 | 80.7 | — | — | — | — | **unreliable** (n-sum) |
| 3.00 | 60 | 40.3 | +1.0000 | +0.9946 | 3.0e-04 | 1.8e-04 | ok (Δφ=3,2) |
| 2.00 | 70 | 70.5 | +0.9995 | +0.9752 | 1.2e-03 | 9.4e-04 | ok (Δφ=3,2) |
| 5.00 | 100 | 40.2 | +0.9999 | +0.9980 | 1.3e-04 | 1.3e-04 | ok (Δφ=3,2) |
| 3.00 | 104 | 69.7 | +1.0001 | +0.9898 | 1.9e-04 | 1.9e-04 | ok (Δφ=3,2) |
| 4.00 | 120 | 60.2 | +1.0000 | +0.9955 | 6.9e-05 | 5.1e-05 | ok (Δφ=3,2) |

**Δ_g = 300**

| `m/Λ` | `J` | `b̂` | `++` `F/N_s` | `+−` `F/N_s` | max `|S|/F_hyb` (strip) | max `|S|/F_hyb` (B′ env) | status |
|---|---|---|---|---|---|---|---|
| 1.50 | 16 | 22.0 | +1.0002 | +0.9489 | 2.7e-04 | 1.7e-04 | ok (Δφ=3,2) |
| 2.00 | 20 | 20.5 | +0.9938 | +0.9886 | 3.4e-04 | 3.7e-04 | ok (Δφ=3,2) |
| 1.50 | 22 | 30.0 | +0.9999 | +0.9880 | 6.1e-04 | 3.4e-04 | ok (Δφ=3,2) |
| 2.00 | 30 | 30.5 | +0.9999 | +0.9880 | 3.2e-04 | 2.2e-04 | ok (Δφ=3,2) |
| 3.00 | 30 | 20.3 | +0.9936 | +1.0001 | 1.3e-04 | 1.3e-04 | ok (Δφ=3,2) |
| 1.50 | 34 | 46.0 | +1.0000 | +0.9885 | 3.1e-04 | 1.2e-04 | ok (Δφ=3,2) |
| 4.00 | 40 | 20.2 | +0.9951 | +1.0001 | 6.6e-05 | 6.3e-05 | ok (Δφ=3,2) |
| 2.00 | 44 | 44.5 | +0.9998 | +0.9949 | 3.4e-04 | 2.8e-04 | ok (Δφ=3,2) |
| 3.00 | 44 | 29.7 | +1.0001 | +0.9902 | 6.2e-05 | 6.0e-05 | ok (Δφ=3,2) |
| 5.00 | 50 | 20.2 | +0.9963 | +1.0000 | 3.9e-05 | 3.7e-05 | ok (Δφ=3,2) |
| 1.50 | 52 | 70.0 | +1.0000 | +0.9919 | 2.2e-04 | 1.2e-04 | ok (Δφ=3,2) |
| 4.00 | 60 | 30.2 | +1.0000 | +0.9951 | 8.9e-06 | 4.3e-06 | ok (Δφ=3,2) |
| 3.00 | 68 | 45.7 | +0.9999 | +0.9978 | 1.2e-04 | 1.2e-04 | ok (Δφ=3,2) |
| 2.00 | 70 | 70.5 | +0.9999 | +0.9952 | 2.1e-04 | 1.7e-04 | ok (Δφ=3,2) |
| 7.00 | 70 | 20.1 | +0.9978 | +1.0000 | 1.8e-05 | 1.8e-05 | ok (Δφ=3,2) |
| 1.50 | 74 | 99.3 | +1.0000 | +0.9832 | 1.8e-04 | 1.1e-04 | ok (Δφ=3,2) |
| 5.00 | 74 | 29.8 | +1.0000 | +0.9965 | 2.3e-05 | 2.3e-05 | ok (Δφ=3,2) |
| 4.00 | 80 | 40.3 | +1.0000 | +0.9982 | 5.3e-06 | 2.7e-06 | ok (Δφ=3,2) |
| 2.00 | 100 | 100.5 | +1.0000 | +0.9900 | 8.0e-05 | 5.5e-05 | ok (Δφ=3,2) |
| 10.00 | 100 | 20.1 | +0.9988 | +1.0000 | 8.2e-06 | 8.1e-06 | ok (Δφ=3,2) |
| 3.00 | 104 | 69.7 | +1.0000 | +0.9980 | 7.4e-05 | 6.1e-05 | ok (Δφ=3,2) |
| 7.00 | 104 | 29.9 | +1.0000 | +0.9982 | 1.1e-05 | 1.1e-05 | ok (Δφ=3,2) |
| 1.50 | 112 | 150.0 | — | — | — | — | **unreliable** (n-sum) |
| 5.00 | 112 | 45.0 | +1.0000 | +0.9992 | 5.3e-05 | 5.1e-05 | ok (Δφ=3,2) |
| 4.00 | 120 | 60.3 | +1.0000 | +0.9990 | 4.0e-05 | 3.6e-05 | ok (Δφ=3,2) |
| 2.00 | 150 | 150.5 | +1.0000 | +0.9759 | 6.6e-05 | 5.2e-05 | ok (Δφ=3,2) |

**Δ_g = 1000**

| `m/Λ` | `J` | `b̂` | `++` `F/N_s` | `+−` `F/N_s` | max `|S|/F_hyb` (strip) | max `|S|/F_hyb` (B′ env) | status |
|---|---|---|---|---|---|---|---|
| 1.50 | 16 | 22.0 | +1.0000 | +0.9491 | 2.4e-05 | 1.5e-05 | ok (Δφ=3,2) |
| 2.00 | 20 | 20.5 | +0.9936 | +0.9882 | 3.1e-05 | 3.3e-05 | ok (Δφ=3,2) |
| 1.50 | 22 | 30.0 | +1.0000 | +0.9891 | 5.4e-05 | 3.1e-05 | ok (Δφ=3,2) |
| 2.00 | 30 | 30.5 | +1.0000 | +0.9886 | 2.9e-05 | 2.0e-05 | ok (Δφ=3,2) |
| 3.00 | 30 | 20.3 | +0.9936 | +1.0000 | 1.2e-05 | 1.1e-05 | ok (Δφ=3,2) |
| 1.50 | 34 | 46.0 | +1.0000 | +0.9898 | 2.9e-05 | 1.1e-05 | ok (Δφ=3,2) |
| 4.00 | 40 | 20.3 | +0.9950 | +1.0000 | 5.9e-06 | 5.6e-06 | ok (Δφ=3,2) |
| 2.00 | 44 | 44.5 | +1.0000 | +0.9961 | 3.1e-05 | 2.5e-05 | ok (Δφ=3,2) |
| 3.00 | 44 | 29.7 | +1.0000 | +0.9903 | 5.5e-06 | 5.3e-06 | ok (Δφ=3,2) |
| 5.00 | 50 | 20.2 | +0.9963 | +1.0000 | 3.5e-06 | 3.4e-06 | ok (Δφ=3,2) |
| 1.50 | 52 | 70.0 | +1.0000 | +0.9980 | 2.0e-05 | 1.1e-05 | ok (Δφ=3,2) |
| 4.00 | 60 | 30.3 | +1.0000 | +0.9952 | 8.4e-07 | 4.3e-07 | ok (Δφ=3,2) |
| 3.00 | 68 | 45.7 | +1.0000 | +0.9983 | 1.1e-05 | 1.1e-05 | ok (Δφ=3,2) |
| 2.00 | 70 | 70.5 | +1.0000 | +0.9988 | 1.8e-05 | 1.5e-05 | ok (Δφ=3,2) |
| 7.00 | 70 | 20.1 | +0.9978 | +1.0000 | 1.6e-06 | 1.6e-06 | ok (Δφ=3,2) |
| 1.50 | 74 | 99.3 | +1.0000 | +0.9982 | 1.5e-05 | 1.0e-05 | ok (Δφ=3,2) |
| 5.00 | 74 | 29.8 | +1.0000 | +0.9965 | 2.1e-06 | 2.0e-06 | ok (Δφ=3,2) |
| 4.00 | 80 | 40.2 | +1.0000 | +0.9983 | 5.6e-07 | 2.6e-07 | ok (Δφ=3,2) |
| 2.00 | 100 | 100.5 | +1.0000 | +0.9989 | 1.0e-05 | 7.0e-06 | ok (Δφ=3,2) |
| 10.00 | 100 | 20.1 | +0.9988 | +1.0000 | 7.4e-07 | 7.3e-07 | ok (Δφ=3,2) |
| 3.00 | 104 | 69.7 | +1.0000 | +0.9996 | 7.7e-06 | 6.5e-06 | ok (Δφ=3,2) |
| 7.00 | 104 | 29.9 | +1.0000 | +0.9982 | 1.0e-06 | 1.0e-06 | ok (Δφ=3,2) |
| 1.50 | 112 | 150.0 | +1.0000 | +0.9942 | 1.3e-05 | 8.3e-06 | ok (Δφ=3,2) |
| 5.00 | 112 | 45.0 | +1.0000 | +0.9994 | 4.8e-06 | 4.6e-06 | ok (Δφ=3,2) |
| 4.00 | 120 | 60.3 | +1.0000 | +0.9995 | 4.0e-06 | 3.6e-06 | ok (Δφ=3,2) |
| 1.50 | 150 | 200.7 | +1.0000 | +0.9870 | 8.7e-06 | 5.7e-06 | ok (Δφ=3,2) |
| 2.00 | 150 | 150.5 | +1.0000 | +0.9967 | 6.0e-06 | 4.9e-06 | ok (Δφ=3,2) |


**Envelope variant.** The B′-envelope variant of `S` is smaller (`|S|/F_hyb ≤ 2.6·10⁻³` at `Δ_g = 100`). The headline uses the stripped
envelope.

**Odd-`k` interpolation sensitivity:** `≤ 2·10⁻⁵ N_s`.

**Pointwise vs integrated remainder.** As for the old functional, the pointwise remainder is up to `1.7·10³ N_s` (at `b̂ ≈ 200`), but its integral is
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

### 4.4 Tighter constants: the window moves to `b̂₀ = 91.2` (`h13a_cert_curv.py`, `H13A_SIGNED=1`) [certified]

**Two changes, both [D]:**
1. **Region I uses the signed integral.** With `(1 − s²/b̂²)^{−½} = 1 + r`, `0 ≤ r ≤ (s²/2b̂²)(1 − s₁²/b̂²)^{−3/2}`:
   * `|(2/π)∫_I G dθ| ≤ (2/(πb̂))[|∫₀^{s₁}G ds| + (1/2b̂²)(1 − s₁²/b̂²)^{−3/2}∫₀^{s₁}|G|s² ds]`;
   * `∫G` is enclosed by the midpoint rule per cell, with error `h³ sup|G″|/24`.
2. **Oscillatory terms use the near-cancellation of `w Osc(τ)` and `Osc(s)`.**
   * `|w Osc(τ) − Osc(s)| ≤ (s²/8Λ²)Ω(τ) + (s³/24Λ²)Ω₁(τ)`, where `Ω₁` bounds `|Osc′|`.
   * This makes `C_osc` scale like `Λ⁻²`.

**Result: certified `curv > 0`, and hence `F_hyb > 0`, on `b̂ ∈ [91.2, 2Λ]` at `Δ_g = 100`, `300` and `1000`.**
* The previous windows started at 115.5 / 162.6 / 274.1.
* At `s₁ = 80`: `(C_IΛ² + C_sm + C_oscΛ²) ≤ 682` (`Δ_g = 100`) and `≤ 683` (`Δ_g = 300`). The window starts at the same box edge for all
  three `Δ_g`, because after the changes the constants are almost `Λ`-independent.

**What limits further progress.**
* For smaller `s₁` the oscillatory constant `C_oscΛ² ≈ 7·10³ (40/s₁)³` dominates. It is large because `ψ^{(m)}(1)` reaches about `4·10⁷`.
* An integration by parts in `θ` would use the phase `e^{ib̂ sin θ}` and gain about `1/s₁` per step. That is the next tightening; it is not done.
* `b̂₀ ≈ 50` is not reached.

### 4.3 Coverage map for route A′ (`m ≥ 1.5Λ`)

| `b̂` | `Δ_g = 100` | `Δ_g = 300` | `Δ_g = 1000` |
|---|---|---|---|
| `< 20` (`< 40` at 100) | near-gap scans and certificate (`f6`, `CERTIFIED_NEARGAP` model) | same | same |
| `[20, b̂₀)` | hybrid [NE], every `J`; exact `S` [NE] to `b̂ = 70` | hybrid [NE]; exact `S` to 150.5 | hybrid [NE]; exact `S` to 200.7 |
| `[b̂₀, 2Λ]` | **hybrid [certified]** (`b̂₀ = 91.2`, §4.4) | **[certified]** (91.2) | **[certified]** (91.2; flat corner [NE]) |
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
  * The other planned states were not run: `J = 150` at `Δ_g = 100` (`3Λ`; expected to fail, `J/τ ≈ 0.8`) and the `J = 200` states (`2–3Λ`).
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
