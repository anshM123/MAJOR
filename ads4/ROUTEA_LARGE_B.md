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
