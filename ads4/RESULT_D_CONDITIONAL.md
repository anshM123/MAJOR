# Result D as a conditional theorem: the failing band, the band hypothesis, the window/bound trade-off, and the flat-space comparison

*19 Nov 2026, conditional-statement agent. Unrefereed, not committed. No existing document was edited (`RESULT_D.md`, the review,
`LARGEJ_OFFDIAG.md`, `S1_S4_FIXES.md`, `LIT_HELICITY_MIXING.md`, `OFFDIAG_KERNEL.md` were only read). New code:
`ads4/code/spinfunc/t1_band.py`, `t2_lp.py`, `t3_q3check.py`, `t4_flatband.py`, `t1_queue.sh`. Caches: `spinfunc/cache/t1_*`, `t2_*`, `t3_*`
(git-ignored). Project modules (`r5_universal`, `r3_redesign`, `q3_analysis`, `robust/nf_lp`) are imported read-only. Every run used
one thread; at most two processes ran at once.*

Grades: **[P]** proved; **[D]** derived; **[D-sketch]**; **[Q]** quoted; **[DV]** derived and checked numerically; **[NE]** numerical
evidence; **[M]** model; **[I]** assumed or open. Conventions: `Λ = √((Δ_g−2)(Δ_g−1))`; `b̂ = (2J+1)Λ/m`;
`β = arccosh(1 + (2J+1)²/(2m²)) ≈ b̂/Λ`; AdS units (`R_AdS = 1`); `k(0) = 1`.

---

## 0. Summary

1. **The band, determined.** I mapped where the per-state 2×2 test fails, `λ_min < 0`, at Δ_g = 100, 300 and 1000. The model is
   LARGEJ's universal-`g` model with fine quadrature, extended to `β ≤ 2`. For NF02 the failing set is, to grid accuracy, **every
   twist-allowed even-`J` state with `b̂ ≥ b̂_c`, for every sampled `m/Λ ∈ [1.02, 5]`**. No sampled state is positive above the edge, up
   to `β = 2` (100, 300) or `β = 0.5` (1000) [M].

   | `Δ_g` | `b̂_c` (NF02) | `b̂_c/Λ^{2/3}` | `β_c` | inside validated `g` (`β ≤ 0.11`)? |
   |---|---|---|---|---|
   | 100 | 16.2 (±1 between `m/Λ` values) | 0.76 | 0.16 | **no**: the whole band is extrapolated |
   | 300 | 30 (per-`m` edges 30–46) | 0.67 | 0.10 | the onset only (`b̂ ≤ 33`) |
   | 1000 | 62–68 (one grid outlier, 107, at `1.5Λ`) | 0.62–0.68 | 0.062–0.068 | yes, up to `b̂ ≈ 110` |

   * The diagonal stays positive throughout the band: `min F_hh/N_h` = 0.69 (100), 0.88 (300), 0.92 (1000) [M]. **The failure is purely
     the helicity off-diagonal.**
   * The computed graviton corrections (q3 model M, S1_S4) put the first negative state at **`b̂ = 30.0` at Δ_g = 300**, the same as the
     universal model [NE + M].
   * The upper end of the band is not determined. In the model the band runs up to the twist limit (`b̂ ≈ 2Λ` at `m = 2.5–5Λ`), but this
     relies on the small-`β` form `β²g(νβ)` at `β ≈ 0.3–2`, where it is unvalidated [I].
2. **Conditional theorem (§3).** Fix route-A′ weights `ψ`. Take the review §2 hypotheses (A1, A2, V, L0–L3) and per-state positivity
   below the window, `b̂ < b̂_w(ψ)`. Add one band hypothesis: the heavy states with `b̂ ≥ b̂_w` contribute at least `−ε` times the
   functional's own Einstein graviton-pole term. Then
   **`√(t₄²+α_T²) ≤ 240 √((1+ε)·16L_A[ψ](Λ)) / Λ⁴ · (1 ± 0.4 % / ±0.04 %)`** [D given the hypotheses].
   * With NF02 and `ε = 0` this gives **3.13·10⁻⁶ (Δ_g = 300), 2.89·10⁻⁸ (1000), 2.21·10⁻⁴ (100)**.
   * Sufficient for `ε = 0`: the helicity coherence `κ = |ρ₊₊,₊₋|/√(ρ₊₊ρ₊₋)` of the heavy spectral matrix at each `(m, J)` in the band
     stays below `κ_crit(b̂)`. `κ_crit` is ≈ 1 at the band edge, 0.1 at `b̂ ≈ 2b̂_c` and 0.003–0.005 deep in the band (tabulated in §2.2)
     [M].
   * Physical reading: the band (impact parameters `Λ^{2/3} ≲ b̂ ≲ Λ`) must be **helicity-coherence-poor**, as it is if it is dominated
     by Regge-factorised t-channel graviton exchange. In that picture helicity flips cost `(bM)⁻⁴`, and the resulting coherence is
     10–10⁴ times below `κ_crit` (§3.3) [I, heuristic].
   * Plausibility, graded honestly: good for string-like spectra, where the band is populated only at `m ≳ b̂Λ/4`, deep in the Regge
     regime. **Not implied by any axiom**: a spectrum with leading-twist large-spin states near the gap, each with both helicity couplings
     (rank one, so `κ = 1`), violates it.
3. **Trade-off (§4).** Re-weighting inside NF02's LP family pushes the edge out at a cost. The bound below uses the exact `p`-integral;
   the hypothesis region is `[b̂_w, top)`.

   | `Δ_g` | edge `b̂_w` → bound |
   |---|---|
   | 300 | **30 → 3.13**, 36.5 → 3.04, ≈ 43–45 → 3.08, 57.5 → 3.38, 76 → 4.53, ≈ 90 → 9.55 (·10⁻⁶; `d = 12`); `b̂ ≥ 90` infeasible at `d = 8` |
   | 1000 | **62–68 → 2.89**, 68.5 → 2.86, 88.6 → 2.90, 117 → 3.26, ≈ 160 → 5.85 (·10⁻⁸; `d = 12`); `b̂ ≤ 200` infeasible |
   | 100 | **16.2 → 2.21**, 20.5 → 2.13, 29.2 → 2.14, 40.5 → 2.45, 48.2 → 2.87, 57.2 → 4.52 (·10⁻⁴); `b̂ ≤ 69` infeasible |

   A cheap window reaches `≈ 1.3Λ^{2/3}`. Beyond `≈ 2–2.7Λ^{2/3}` no window is feasible. No functional in the family removes the band
   [M].
4. **Flat comparison (§5).**
   * Route A′'s bound is *exactly* a CHLPS-form flat bound with the NF02 functional, with the AdS curvature acting as an IR regulator
     `m_IR = 3/(2R_AdS)`: `|ĝ₃|²Λ⁸ ≤ 3022.34 log(2Λ/3) − 5327.21` (+ O(1/Λ)) [D]. The normalisation of `ĝ₃` is identical to CHLPS's
     [DV: the project reproduces CHLPS (3.7)'s `1598 log − 3211` as `1596.7, −3208.0`].
   * **The log coefficients do not agree, and should not:** 3022.34 vs 24.9, a factor 121.4, i.e. ×11.0 in the bound asymptotically.
     Both equal `16ψ₂′(0)/k` of their own functional, so the coefficient is functional-dependent, not universal.
   * The factor decomposes as 1.05 (CHLPS's extra rules), × 33.6 (no forward terms, forced in AdS by contacts and H9a), × 2.58 (margin
     δ = 0.2), × 1.34 (large-`b̂` certifiability proxy).
   * Transplanting CHLPS's optimum with `m_IR = 3/(2R)` gives 2.23·10⁻⁵ / 3.09·10⁻⁷ / 2.80·10⁻⁹. Result D is weaker by 9.9 / 10.1 / 10.3.
   * Both bounds are conditional on an IR band. For NF02 with a sharp flat cutoff, the band is `b̂ ∈ [1.35Λ^{2/3}, 2.56Λ]` (`t4`) [NE]. The
     AdS graviton band starts about 2× earlier (`0.62–0.76Λ^{2/3}`), and in AdS the diagonal is positive there.
   * Parametrically Result D is `√(log Δ_g)/Δ_g⁴`, stronger than AKT's `log Δ_g/Δ_g⁴`. In AKT's form the effective constant is
     `C ≈ 4.2–4.8·10³`.
5. **Caveat on inputs.** All band numbers use LARGEJ's `g_r4` (r5). `OFFDIAG_KERNEL.md`, work in progress by another agent, finds the
   basis-invariant real off-diagonal to be `β²R(νβ)` with `∫tR = 47/4`, and finds that part of `g_r4` is a basis shear. Since
   `b̂_c ∝ |G₁|^{−1/3}` (§2.1), an O(10–25 %) change in the effective `G₁` moves the edges by 3–8 % [D-sketch]. The same note reports a
   flip-free (2,4) combination `S(κ*)` whose scaling-limit off-diagonal vanishes. If that is adopted, the band and this whole conditional
   layer must be redone, and the band hypothesis may become unnecessary.

---

## 1. Inputs used (latest bullets of `RESULT_D.md` supersede earlier ones)

* **Constants.** A, B_AdS, `16L_A` and the implication "per-state positivity ⇒ bound" are [D]/[DV] (review §0, (a)). `A = 3022.34`,
  `B_AdS = −6552.66`. `16L_A = 7508.5 / 10734.6 / 14339.6` at Δ_g = 100 / 300 / 1000 (reproduced here by `t1_band.LA_exact`).
* **S1/S4** (`S1_S4_FIXES.md`):
  * `Ĉ₃` needs no de-shift [P given kernel (D)].
  * Near-gap graviton corrections break `n = 1` at Δ_g = 300 for even `J ≥ 16–18`, `m/Λ = 1.06–1.25` [NE + M; direct-validated].
  * At Δ_g = 1000 `n = 1` holds for `J ≤ 20` [NE + M] and for `J = 24, 28` [M]; `J = 32` fails [M].
* **LARGEJ_OFFDIAG.** The real off-diagonal is leading-Regge, `Γ_off/m = β²g(νβ)` [NE, `β ≤ 0.11`]. The failure is `m`-independent at
  fixed `b̂`, with `b̂_c ≈ 0.65Λ^{2/3}`. Re-weighting cannot remove it [M + D-sketch].
* **LIT_HELICITY_MIXING.** This is the AdS twin of the flat D = 4 IR band (CHLPS 2201.06602 (3.6); CHMRS 2106.10274 (5.59)). The
  literature states its bounds conditionally on the large-`b` spectrum [Q].

---

## 2. The band, determined from the model (`t1_band.py`, `t3_q3check.py`) [M unless stated]

**Model.** This is r5's universal-`g` 2×2: diagonal = flat NF02 `B₂ + B₃` + leading-Regge curvature (= `F_hyb`); off-diagonal =
`eo2·β²g(νβ)ψ₂ + (k = 3)`. `t1_band` adds composite Gauss–Legendre quadrature (the 64-node rule of r5 is adequate only for `Λβ ≲ 100`;
the two rules agree to ≤ 2 % at `β ≤ 0.3`). States are twist-allowed with r5's convention `α ≥ Δ_g − 2`, which is slightly more
permissive than A1. The grid is `m/Λ ∈ {1.02, 1.05, 1.1, 1.2, 1.3, 1.5, 1.75, 2, 2.5, 3, 4, 5}`, 70 values of `b̂ ∈ [6, 3Λ]`, `β ≤ 2`
(`≤ 0.5` at 1000, where mpmath's `₂F₁` fails beyond).

### 2.1 NF02: where `λ_min < 0` (`cache/t1_band_{100,300,1000}_ext.txt`)

Per `m/Λ`, the cell gives the last positive `b̂` / the first negative `b̂`:

| `m/Λ` | 1.05 | 1.1 | 1.2 | 1.3 | 1.5 | 2 | 3 | 5 |
|---|---|---|---|---|---|---|---|---|
| Δ_g = 100 | none allowed | 15.5 / 19.1 | 14.2 / 17.5 | 13.1 / 16.2 | 14.0 / 16.7 | 14.5 / 16.5 | 15.0 / 16.3 | 13.8 / 18.6 |
| Δ_g = 300 | 27.6 / none (twist limit) | 26.4 / 30.0 | 27.5 / 34.2 | 25.4 / 31.5 | 27.3 / 30.0 | 22.5 / 34.5 | 29.7 / 45.7 | 27.4 / 39.4 |
| Δ_g = 1000 | 46.7 / 61.9 | 51.8 / 62.7 | 57.5 / 67.5 | 62.3 / 68.5 | 62.0 / 107 | 62.5 / 68.5 | 52.3 / 68.3 | 57.0 / 68.2 |

* In every row the negative interval extends from the first negative `b̂` to the largest twist-allowed `b̂` on the grid (e.g.
  `[30, 670]` at 300 and `m = 5Λ`).
* Spacing between neighbouring even-`J` grid states makes the per-`m` edges scatter by ±15 %. The universal edge is
  `b̂_c ≈ 16 / 30 / 65`, i.e. `(0.62–0.76)Λ^{2/3}`.
* **Check with the actually computed graviton corrections** (`t3_q3check`: q3 model M of S1_S4; K₂ rows, aligned off-diagonal, fixed
  `Ĉ₃`; 888 blocks at 300).
  * NF02's first negative state is at **`b̂ = 30.0`**: `J = 16`, `1.06Λ` (`b̂ = 31.1`, −0.16); `J = 18` (34.9, −0.85); `J = 20` (38.3, −1.34).
    Hybrid transfers (`J = 24, 28, 32` at `1.4Λ`) are negative at `b̂ = 35–46`. `J ≤ 14` (`b̂ ≤ 27.6`) is positive.
  * The universal model and model M agree on the edge [NE + M].
  * At 1000 the q3 set reaches only `b̂ ≈ 55` for `J ≤ 28`, and it is all positive. Its `J = 32` row has the broken K₂ diagonal noted in
    S1_S4 §2.6 (−841, meaningless). So the 1000 edge is tested by the universal model only.
* **Scaling of the edge [D-sketch].** In the wedge, `|X|/F ≈ 2ψ₂′(0)|G₁|b̂³/(1023.6Λ²)`, so `b̂_c ∝ (Λ²/|G₁|)^{1/3}`, i.e.
  `b̂_c/Λ^{2/3}` is constant. Measured: 0.76 / 0.67 / 0.62. The residual drift is the `t`-shape of `g` relative to `b̂p`.

### 2.2 Critical helicity coherence `κ_crit = √(F₊₊F₊₋)/|X|` (envelope, min over states in the bin; NF02)

For a PSD heavy density `ρ` at a state, `tr(ρF) ≥ 0` for **every** `ρ` with `|ρ₊₊,₊₋| ≤ κ√(ρ₊₊ρ₊₋)` iff `κ ≤ κ_crit` [D]. This requires
`F_hh ≥ 0`, which holds in the band [M].

| `b̂` bin | 15–20 | 25–30 | 30–35 | 40–50 | 60–70 | 80–100 | 100–120 | 150–200 | 200–300 | 300–500 |
|---|---|---|---|---|---|---|---|---|---|---|
| Δ_g = 100 | 0.49 | 0.13 | 0.069 | 0.026 | 0.011 | 0.0062 | 0.0048 | 0.0029 | 0.0026 | — |
| Δ_g = 300 | 4.4 | 0.95 | 0.60 | 0.20 | 0.081 | 0.028 | 0.016 | 0.0060 | 0.0041 | 0.0030 |
| Δ_g = 1000 | 49 | 10.5 | 7.0 | 2.2 | 0.75 | 0.25 | 0.14 | 0.039 | 0.0125 | 0.0051 |

`κ_crit` tends to ≈ 1/(24|G₁|·(1 + u-terms)) ≈ 1/400–1/300 deep in the band, matching LARGEJ's weight-independent asymptotic
`|X|/F → 24|G₁|` [D-sketch]. Entries at `β > 0.11` (100: all; 300: `b̂ > 33`; 1000: `b̂ > 110`) use the extrapolated `g` [M + I].

### 2.3 Δ_g = 100: a tension, not resolved here

The universal model puts near-gap `J = 10` at `m = 1.1Λ` (`b̂ = 19.1`, `β = 0.19`) in the band, with `λ/N_s = −0.61`. The Δ_g = 100
near-gap chain (H0 §25, Ĉ₂′, contour-validated at `x ≥ 0.1`) has `J = 10` positive (0.49–0.50).
* Either `g` at `β ≈ 0.2` is smaller than extrapolated (which would favour a narrower band), or the Δ_g = 100 chain's off-diagonal is
  affected by the row-sign defect of S1_S4 §2.2(a) (the `(4,100)` and `(10,200)` rows were flipped; the `Δ = 100` rows of `J = 8, 10` are
  unconverged).
* OFFDIAG_KERNEL's exact kernel at `β ≈ 0.2` should decide it. **At Δ_g = 100 the band's location is [M + I], not [NE].**

---

## 3. The conditional theorem

### 3.1 Statement

**Theorem D_cond** [D given the hypotheses; the constants are DV].

*Setting.* Let `Δ_g ∈ {300, 1000}`; Δ_g = 100 is included with the extra caveats of §2.3 and RESULT_D's [I] near-gap status. Let
`Ω^{A′}_ψ = Re∫XΨ_{2,ν}[q²]/(h₋+h₊) + Re∫ψ₃½Λ⁶Ĉ_{3,ν}` be the route-A′ functional of H0 §R. It uses `Ĉ₂′`, `Ĉ₃` without de-shift (S1_S4),
min-norm gauge, and weights `ψ = (ψ₂, ψ₃)` from NF02's LP family with `k(0) = 1`. Let `b̂_w(ψ)` be the window edge of §4.

*Hypotheses.*
* **A1, A2, V, L0, L1, L2, L3** of the review §2: spectrum, tree level, validity/swap, light content, H7, H17, H6-lin. Grades as there:
  [A], [A], [D], [D], [D], [D-sketch + M], [I].
* **P_<**: every twist-allowed heavy state with `b̂ < b̂_w(ψ)` has a PSD per-state action matrix (EE 2×2 ⊕ OO 1×1).
* **H_band(ε)**: the heavy states with `b̂ ≥ b̂_w(ψ)` together contribute to `Ω^{A′}_ψ` at least `−ε` times the graviton-pole term
  `16L_A[ψ]` of the same functional (in the units in which the sum rule reads `Σ_heavy = 16L_A − k|ĝ₃|²Λ⁸`).

*Conclusion.*
`k(0)(1+δ)|ĝ₃|²Λ⁸ ≤ (1+ε)·16L_A[ψ](Λ)`, i.e.

`√(t₄² + α_T²) ≤ 240·√((1+ε)·16L_A[ψ](Λ)) / Λ⁴ · (1 ± 0.4 % (300) / ±0.04 % (1000))`,

with `L_A[ψ] = ∫₀¹ψ₂(t)dt/(t² + 9/(4Λ²))`, `t₄`, `α_T` as in `DICTIONARY_T4.md`, and the band edge `b̂_w` of §2/§4.

*Proof.*
1. Heavy ≥ 0 below the window (P_<).
2. The band term is `≥ −ε·16L_A` (H_band).
3. Hence `−Ω_light ≥ −ε·16L_A`.
4. `−Ω_light = 16L_A − k|ĝ₃|²Λ⁸` (review §0(a), L0–L2).

∎ [D]

**Linear variant.** Suppose the band's negativity is linear in `ĝ₃`, as expected if the off-diagonal heavy density is sourced by the R³
vertex (§3.3): `Σ_band ≥ −2γ√(16L_A)|ĝ₃|Λ⁴`. Then `|ĝ₃|Λ⁴ ≤ √(16L_A)(γ + √(1+γ²))` [D]. For example, `γ = 0.1` gives ×1.105 on the bound.

**Numbers, NF02** (`ψ` = NF02, `b̂_w = b̂_c`):

| `ε` | 0 | 0.1 | 0.5 | 1 | 3 | 10 |
|---|---|---|---|---|---|---|
| Δ_g = 300 (·10⁻⁶) | **3.13** | 3.29 | 3.84 | 4.43 | 6.26 | 10.4 |
| Δ_g = 1000 (·10⁻⁸) | **2.89** | 3.03 | 3.54 | 4.09 | 5.78 | 9.59 |
| Δ_g = 100 (·10⁻⁴) | **2.21** | 2.32 | 2.71 | 3.12 | 4.42 | 7.33 |

### 3.2 Status of each premise

| premise | Δ_g = 300 | Δ_g = 1000 | Δ_g = 100 |
|---|---|---|---|
| P_< near the gap | **[NE + M]**: even `J ≤ 14`, odd `J ≤ 13` (`b̂ ≤ 27.6 < 30`), graviton-corrected (S1_S4); `J = 16` (`b̂ = 31`) lies just inside the band | **[NE + M]** `J ≤ 20`; [M] `J = 24, 28` (all `b̂ ≤ 55 < 62`) | [NE + M] on `x ≥ 0.1`, [I] at `x = 0.05` (RESULT_D); the band starts at `J ≈ 8` |
| P_< away from the gap (`m ≥ 1.3Λ`, `b̂ < b̂_c`) | [M]: hybrid + universal `g` (validated range), positive; graviton diagonal corrections beyond the computed rows [I] | same | [M] (extrapolated `g`) |
| P_< wedge `m ∈ (1.3, 1.5)Λ`, `J > 20` | [M] (review §4 item 2 + this note) | [M] | [M] |
| band location | onset [NE + M] (model M = universal model); above `b̂ ≈ 33` [M + I] | onset [M] (validated `g`); `b̂ > 110` [M + I] | [M + I] (§2.3) |
| **H_band** | **[I]** (§3.3) | **[I]** | **[I]** |
| H_band(0) needs `κ ≤` | ≈ 0.6 at the edge → 0.003 | ≈ 0.75 → 0.005 | ≈ 0.5 → 0.003 |

What the theorem does **not** say: that the bound holds for every CFT satisfying A1–A2. It holds for those CFTs whose band spectrum
satisfies H_band. The band contains the near-gap states `J ≥ 16` (300) and `J ≥ 32` (1000), and everything at `b̂ ≥ b̂_c` for larger `m`.

### 3.3 Interpretations of H_band, from strongest to weakest, and plausibility

Write `ρ(m, J) = Σ_{O at (m,J)} c_O c_Oᵀ` for the heavy spectral 2×2 in the (++, +−) basis. In CHLPS notation this is the matrix
`Im a_J(s)`, whose off-diagonal entry is `Im a^{+++−}_J = Σc⁺⁺c⁺⁻`.

| form | statement | implies | grade |
|---|---|---|---|
| **H-diag** | each band level is helicity-diagonal, `Σ_O c⁺⁺_O c⁺⁻_O = 0` | H_band(0) (the diagonal is positive [M]) | [I] |
| **H-κ** | coherence `κ(m,J) ≤ κ_crit(m,J)` (§2.2) at every band level | H_band(0) | [I] |
| **H-eik** | in the band, `ρ` is that of Regge-factorised t-channel graviton exchange: diagonal ∝ eikonal density, `κ ≈ κ_eik(b̂) = 24|ĝ₃|Λ⁴/(b̂⁴ℓ)` with `ℓ = log(2Λ/(3b̂))` (CHLPS (3.25), `[[log, 24ĝ₃/b⁴],[c.c., log]]`) | H-κ with large margin (below) | [I], heuristic |
| **H-ε** | aggregate: band contribution ≥ `−ε·16L_A` | the theorem | [I] (weakest) |

**H-eik is self-consistent, with margin.** At the bound, `|ĝ₃|Λ⁴ ≤ √(16L_A) = 104 (300) / 120 (1000) / 87 (100)`. With that value,
`κ_eik` compares with `κ_crit` as follows:
* **Δ_g = 300**, `b̂ = 30 / 60 / 100 / 150`: `κ_eik ≈ 1.6·10⁻³ / 1.6·10⁻⁴ / 3.6·10⁻⁵ / 1.8·10⁻⁵`, against `κ_crit ≈ 0.6 / 0.08 / 0.016 / 0.006`.
  The margin is ≈ 300–500.
* **Δ_g = 1000**: the margin is ≳ 10⁴ for `b̂ ≤ 300`.
* **Δ_g = 100**: the margin is ≈ 15–30.
* `ℓ → 0` as `b̂ → 2Λ/3`, where the flat eikonal log is no longer the AdS propagator. So this check covers `b̂ ≲ Λ/2` only. [I/heuristic;
  CHLPS's matrix is quoted from `LIT_HELICITY_MIXING.md` §2.]

**Arguments for plausibility.**
1. *Angular momentum + Regge factorisation.* A helicity flip `|Δh| = 4` on one line makes the t-channel Regge residue vanish like
   `(−t)²`. In impact-parameter space this is `(bM)⁻⁴` relative to the helicity-preserving eikonal, with the coefficient set by the
   non-minimal (R³) vertex. This is CHLPS (3.25)'s structure [Q]. The suppression is generic, not specific to a model [D-sketch, flat
   space].
2. *Where string-like spectra live* [M]. Suppose the leading trajectory has `J_max(m) ≈ 2m²/M² + 2` (type II slope, `M² = 4/α′` ↔ `Λ²`).
   Then `b̂ ≤ (4m²/Λ² + 5)Λ/m`, so the band `b̂ ≥ b̂_w` is populated only at `m ≳ b̂_wΛ/4`: **m ≳ 7.5Λ (300), ≳ 16Λ (1000)**.
   * The equivalent AdS statement: the twist at spin `J` grows like `Δ_g√(J/2) − J`, so near-gap large-`J` states are absent.
   * The band is then reached only deep in the Regge regime `s ≫ M²`, `b ≫ 1/M`, where reggeised graviton exchange is the standard
     expectation (Häring–Zhiboedov 2410.21499 [A]; the project's own flat hypothesis `E(λ, η)`, `FLAT_GRAVITON_BOUND.md` §4).
3. *Sum-rule link.* The flat helicity-violating rule `B⁽³⁾₂` ties smeared `Σc⁺⁺c⁺⁻` to a low-energy side linear in `ĝ₃` (CHLPS (2.35b))
   [Q]. So the integrated off-diagonal heavy density is "of order `ĝ₃`", which is the linear variant. This is not a per-band statement,
   and its CFT uplift does not exist (CLS-D §5) [I].
4. *Precedent.* CHMRS ("only rigorous if one assumes additional information about the spectrum at large impact parameters") and CHLPS
   state their D = 4 bounds under the same kind of condition [Q]. In AdS the condition is milder in kind than in flat space: the
   diagonal is positive in the band (curvature cures it, `M₁ > 0` [certified for the scalar hybrid]), so only helicity coherence must be
   controlled. The flat cutoff band is negative in both channels.

**Arguments against.**
1. *Rank one.* A single non-degenerate heavy operator with both couplings has `κ = 1`. H-κ deep in the band (`κ_crit ~ 10⁻²`) therefore
   needs helicity-selective couplings or cancellations among near-degenerate operators.
2. *Twist gap ≠ Regge slope.* A1 allows near-gap leading-twist states at large `J` (`m ≈ 1.06Λ`, `J ≥ 16` at 300). These are inside the
   band, and no argument above constrains their OPE data. A CFT with a flat leading trajectory at twist ≈ Δ_g and generic helicity
   couplings would violate H_band.
3. *AdS-specific regime.* The upper half of the band, `β ~ 0.3–2` (impact parameter ~ `R_AdS`), is outside both the flat-space
   intuition and the validated model.

**Grade of H_band:** [I]. Plausible (my estimate: likely) for string-like UV completions with a Regge slope `~1/Λ²`. Not supported for
general CFTs satisfying only A1. Untestable with the project's present tools. It has the same epistemic status as the large-`b`
assumptions of CHMRS/CHLPS, and a narrower content (coherence only).

---

## 4. Trade-off: covered window vs bound (`t2_lp.py`, r5 LP caches, `t1_band.py` edges; [M])

Re-weighting is done inside NF02's LP family (`r = 5`, `d = 8` unless noted): flat constraints, margin `δ = 0.2`, proxy 120, `k(0) = 1`,
universal-`g` AdS rows on `β ≤ β_w`, AdS margin `dw`. Only the weights change, so the graviton pole, `Ĉ₃|_E = 0` and (HN) are untouched.
* **"edge"** is the first negative `b̂` of the resulting functional on the fine `t1` map (`t3`'s computed states in brackets where
  checked). The band hypothesis applies to `[edge, top)`.
* **"bound"** is `240√(16L_A)/Λ⁴` with the exact integral. The log form `A log Λ + B_AdS` underestimates it by up to 35 % when `A ≳ 10⁴`.
* **"g valid?"** asks whether the window lies inside `β ≤ 0.11`.

| Δ_g | functional (LP rows, `dw`) | `A` | `B_AdS` | edge `b̂_w` (`/Λ^{2/3}`) | **bound** | g valid? |
|---|---|---|---|---|---|---|
| 300 | **NF02** | 3022.3 | −6552.7 | **30** (0.67) [q3: 30.0] | **3.132·10⁻⁶** | yes |
| 300 | `β ≤ 0.12`, 0.1 | 3131.1 | −7802.7 | 36.5 (0.82) [q3: 34.9] | 3.044·10⁻⁶ | borderline |
| 300 | `β ≤ 0.15`, 0.1 | 3526.8 | −9866.5 | ≈ 43–45 (≈ 1.0) [q3: 43.3] | 3.084·10⁻⁶ | no |
| 300 | `β ≤ 0.2`, 0.1 | 5091.0 | −16981.9 | 57.5 (1.29) [q3: none negative up to 46] | 3.376·10⁻⁶ | no |
| 300 | `β ≤ 0.25`, 0 | 11706 | −45800 | 76.0 (1.70) | 4.530·10⁻⁶ | no |
| 300 | `β ≤ 0.3`, 0, `d = 12` | 63585 | −272727 | ≈ 90 (2.0) (design) | 9.553·10⁻⁶ | no |
| 300 | `β ≤ 0.3`, `d = 8` | — | — | infeasible | — | — |
| 1000 | **NF02** | 3022.3 | −6552.7 | **62–68** (0.62–0.68) | **2.891·10⁻⁸** | yes |
| 1000 | `β ≤ 0.066`, 0.1 | 3051.0 | −7025.6 | 68.5 (0.69) | 2.864·10⁻⁸ | yes |
| 1000 | `β ≤ 0.09`, 0.1 | 3467.4 | −9558.8 | 88.6 (0.89) | 2.901·10⁻⁸ | yes |
| 1000 | `β ≤ 0.12`, 0.1 | 5005.2 | −16509.3 | 117 (1.17) | 3.257·10⁻⁸ | borderline |
| 1000 | `β ≤ 0.16`, 0, `d = 12` | 20725 | −85513 | ≈ 160 (1.6) (design) | 5.846·10⁻⁸ | no |
| 1000 | `β ≤ 0.2`, `d = 8` and `12` | — | — | infeasible | — | — |
| 100 | **NF02** | 3022.3 | −6552.7 | **16.2** (0.76) | **2.209·10⁻⁴** | no |
| 100 | `β ≤ 0.2`, 0.1 | 3106.5 | −7583.6 | 20.5 (0.96) | 2.126·10⁻⁴ | no |
| 100 | `β ≤ 0.3`, 0.1 | 4320.6 | −13516.1 | 29.2 (1.37) | 2.144·10⁻⁴ | no |
| 100 | `β ≤ 0.4`, 0.1 | 7623.2 | −28089.6 | 40.5 (1.90) | 2.446·10⁻⁴ | no |
| 100 | `β ≤ 0.5`, 0.1 | 12707 | −50219 | 48.2 (2.26) | 2.868·10⁻⁴ | no |
| 100 | `β ≤ 0.6`, 0 | 39924 | −168177 | 57.2 (2.68) | 4.518·10⁻⁴ | no |
| 100 | `β ≤ 0.7`, `d = 8` and `12` | — | — | infeasible | — | — |

Notes.
* **Windows on actual computed states.** At 300, the `β ≤ 0.2` and `β ≤ 0.25` functionals are positive on all 888 q3 blocks (graviton
  corrections for `J ≤ 20` + hybrid transfers `J = 24–32`; min 0.1996). NF02 has 30 negative states. The `β ≤ 0.12` / `β ≤ 0.15`
  functionals have 12 / 2 negative states, starting at `b̂ = 34.9 / 43.3` (`t3_q3check_300.txt`) [NE + M]. These match the universal
  edges.
* **Cost.** Moving the edge to `≈ 1.2–1.3Λ^{2/3}` costs 8 % (300) / 13 % (1000) in the bound; at 100 it even improves it, because NF02's LP
  objective was ≈ `A`, not the bound. Beyond ≈ 1.7Λ^{2/3} the cost explodes. **Every functional fails again just past its window**: no
  sampled state is positive above the new edge (`t1` maps) [M].
* The "modest cost" choice at 300 (`b̂ ≤ 57`, bound 3.38·10⁻⁶) puts H_band on `b̂ ≥ 57`, where `κ_crit` starts at ≈ 0.9 (`t1` envelope of
  that functional) instead of 0.6. The gain in plausibility is that the near-gap states `J = 16–28` (`b̂ ≤ 46`) are no longer in the
  band. **The band then begins only at `m ≳ 1.2Λ`, `J ≳ 30` (`b̂ ≈ 57`), and under the string-like spectrum of §3.3 only at `m ≳ 14Λ`.**
  This is the recommended conditional headline at Δ_g = 300 if one prefers H_band to say nothing about near-gap states: **`√(t₄²+α_T²)
  ≤ 3.38·10⁻⁶` (300)**. At 1000 the analogue is `b̂ ≤ 117`, bound **3.26·10⁻⁸**. It keeps every state with `m ≤ 1.05Λ` out of the band (their twist limit is `b̂ ≤ 89`).
* The design for `g/2` (the cached `r5_lp_300_b0.3_dw0.0.pkl` is this functional, A = 10148.8) has, under the full `g`, edge 76 and bound
  4.28·10⁻⁶. It is dominated by the `β ≤ 0.25` row.

---

## 5. Comparison with the sharp flat bound and with AKT

### 5.1 Normalisation [DV]

* **`ĝ₃`.** The project's flat light side is `−B⁽¹⁾₂|_low = 16πG/p² + 2πG|ĝ₃|²p⁶`. This is CHLPS (2.31a) verbatim. The project reproduces
  CHLPS's footnote functional (3.7), `1598 log − 3211`, as `A = 1596.7`, `B = −3208.0` (`NOFORWARD_FUNCTIONAL.md` §0.2). So `ĝ₃` and the
  form `|ĝ₃|²M⁸ ≤ A log(M/m_IR) + B` mean the same thing in both.
* **`M ↔ Λ`.** The flat variable is `p = ν/Λ`, heavy states have `m ≥ Λ`, so `M R_AdS ↔ Λ = √((Δ_g−2)(Δ_g−1)) ≈ Δ_g − 3/2`.
* **`m_IR ↔ 3/(2R_AdS)`, exactly at log accuracy [D].** The AdS graviton pole is `16πG/(ν² + 9/4)` (H7, [D]). So
  `16L_A = A log(Λ/(3/2)) + 16L_fin + O(1/Λ)`, and Result D *is* the CHLPS-form flat bound of NF02 with that cutoff:

  **`|ĝ₃|²Λ⁸ ≤ 3022.34·log(2Λ/3) − 5327.21`** (`B_flat = B_AdS + A log 1.5 = −5327.21`, NOFORWARD's `B`).

  The flat literature's heuristic "log(M R_AdS)" is fixed here to `log(2MR_AdS/3)`.

### 5.2 Log coefficients: 3022.34 vs 24.9

They do not agree: the ratio is **121.4**, i.e. **√121.4 = 11.0 in the bound** as `Λ → ∞`.

* **They are not expected to agree.** For any functional, the log coefficient is `16ψ₂′(0)/k(0)`: the graviton pole integrated against
  the weight's slope at `p = 0`. CHLPS's 24.9 is the SDPB **optimum over a much larger class**. NF02 is one member of a much smaller one.
  "Route A′'s flat limit is NF02" means the coefficient is NF02's, not CHLPS's.
* NF02 is *in the class of* CHLPS's forward-free footnote functional (3.7) (`∫ψ₂B₂ + ∫ψ₃B₃`). That functional's `A ≈ 1597` is smaller,
  but it is not positive on all heavy states: it violates `ψ₃(0) = 0` and the `(1−p)⁵` endpoint condition (NOFORWARD §0.2) [NE].
* Decomposition (from NOFORWARD §0.2; products [NE]):

  | step | `A` | factor |
  |---|---|---|
  | CHLPS optimum: improved `B⁽¹⁻⁴⁾₂`, `k = 3`, forward `B⁽¹⁾₄`, SDPB | 24.9 | — |
  | project LP, forward-term class, margin 0 | 26.1 | 1.05 |
  | **forward-free** `∫ψ₂B₂ + ∫ψ₃B₃` (`r = 5`, `d = 8`), margin 0 | 876 | **33.6** |
  | margin δ = 0.2 (needed against AdS / graviton corrections) | 2258 | 2.58 |
  | + Lemma-M proxy (large-`b̂` certifiability) = NF02 | 3022.34 | 1.34 |

* **Why route A′ must pay the factor 33.6.**
  * Forward terms (`∂_{p²}B₄|₀`, `B₅|₀`) uplift to `Λ¹⁰`-weighted `ν = 0` moments. Their leverage on light-side errors is ~10¹³, and their
    double-trace/contact content is uncontrolled (H9a).
  * No forward-free smeared `B_{k≥4}` combination is contact-free (NOFORWARD obstruction, [D]).
  * CLS-D's family has no helicity-flip `B⁽³⁾` uplift (LIT §1).
  * This is the price of a contact-free AdS uplift, not a defect of the comparison.
* **Constants.** CHLPS `B/A = −1.11`; NF02 `B_flat/A = −1.76`.

### 5.3 Numbers, with CHLPS transplanted by the same identification (`24.9 log(2Λ/3) − 27.6`)

| Δ_g | Result D (NF02, exact `16L_A`) | CHLPS optimum transplanted | ratio | `24.9 log Λ − 27.6` variant |
|---|---|---|---|---|
| 100 | 2.21·10⁻⁴ | 2.23·10⁻⁵ | 9.9 | 2.37·10⁻⁵ |
| 300 | 3.13·10⁻⁶ | 3.09·10⁻⁷ | 10.1 | 3.23·10⁻⁷ |
| 1000 | 2.89·10⁻⁸ | 2.80·10⁻⁹ | 10.3 | 2.90·10⁻⁹ |

The transplanted CHLPS number is **not** an established AdS bound. CHLPS themselves only assert the translation (§5 of 2201.06602). It
would need the forward and helicity-flip rules in AdS, and it carries the same kind of band condition.

### 5.4 The band compared (`t4_flatband.py`) [NE]

* **Flat NF02 with a sharp cutoff** `ε = m_IR/M = 1.5/Λ` (scaling limit, `2∫_ε¹ψ₂J₀(b̂p)dp`). It is negative on
  **`b̂ ∈ [1.35Λ^{2/3}, 2.56Λ]`**, i.e. `[1.77ε^{−2/3}, 3.83ε^{−1}]`, followed by oscillating negative windows. The values are 29.4 / 60.4 /
  134 at Δ_g = 100 / 300 / 1000. Analytic small-`εb̂` estimate: `b̂_c = (1023.6/(2ψ₂′(0)·½))^{1/3}ε^{−2/3} = 1.76ε^{−2/3}` [D-sketch]. This is
  CHLPS's `m_IR^{−2/3} ≲ b ≲ m_IR^{−1}` for this functional.
* **AdS, NF02.** The diagonal is positive at all `b̂` (curvature, `M₁ > 0`). The off-diagonal band starts at `0.62–0.76Λ^{2/3}`, about 2×
  earlier. In the wedge the AdS off-diagonal/base ratio is `≈ 3.7b̂³/Λ²` (LARGEJ), against the flat cutoff deficit `≈ 0.42b̂³/Λ²`, i.e. about
  9× larger.
* **Same scaling, different mechanism.** Flat: a diagonal IR deficit. AdS: a curvature-induced helicity off-diagonal. Hence the AdS
  hypothesis can be phrased as a helicity-coherence condition. The flat one cannot.

### 5.5 AKT 1805.07393 and the other parametric results

* AKT: `|t₄| ≲ ln Δ_gap/Δ_gap⁴` (d = 3, HNEC, unknown constant) [Q].
* Result D is **`√(ln Δ_g)/Δ_g⁴`** times `240√A ≈ 1.3·10⁴` (asymptotically), stronger by `√(ln Δ_g)` [D given the hypotheses]. The √
  comes from the light side being quadratic in `ĝ₃`, and the log from the D = 4 graviton pole, matching CHLPS.
* Written in AKT's form, `√(t₄²+α_T²) ≤ C·ln Δ_g/Δ_g⁴` with **C = 4.8·10³ / 4.4·10³ / 4.2·10³** at Δ_g = 100 / 300 / 1000 (`C` decreases
  like `1/√ln Δ_g`).
* AKT's constant is not known, so the numbers cannot be compared. CEMZ / CHP / Meltzer–Perlmutter give `Δ_g⁻⁴` with unknown constants.
* AKT's quantity `|ñ_s − 4ñ_f|/c_T` is proportional to `t₄`, with a normalisation not converted here [I].
* **Result D (conditional) remains, to my knowledge, the only explicit-constant bound of this type in d = 3** (LIT §3) [Q, literature
  check of LIT_HELICITY_MIXING].

---

## 6. Draft replacement for the obsolete RESULT_D tables (route A′, Ĉ₂′ + Ĉ₃; base: review §3, updated)

### 6a. Headline (replaces §2)

| Δ_g | `Λ` | `16L_A` | **conditional bound** `√(t₄²+α_T²)` (NF02, `ε = 0`) | band hypothesis on | alternative (re-weighted, H_band only beyond) | status |
|---|---|---|---|---|---|---|
| 100 | 98.499 | 7508.5 | **2.21·10⁻⁴** (±2 % H17) | `b̂ ≥ 16` (extrapolated `g`) | `b̂ ≥ 29`: 2.14·10⁻⁴ | conditional on H_band **and** near-gap `n` at `x = 0.05` [I]; band location [M + I] |
| 300 | 298.50 | 10734.6 | **3.13·10⁻⁶** (±0.4 %) | `b̂ ≥ 30` (incl. near-gap `J ≥ 16`) | `b̂ ≥ 57`: **3.38·10⁻⁶** (no near-gap states in the band) | conditional on H_band [I]; P_< [NE + M] near gap, [M] elsewhere |
| 1000 | 998.50 | 14339.6 | **2.89·10⁻⁸** (±0.04 %) | `b̂ ≥ 62` (incl. near-gap `J ≥ 32`) | `b̂ ≥ 117`: **3.26·10⁻⁸** | conditional on H_band [I]; P_< [NE + M] `J ≤ 20`, [M] above |

* `A = 3022.34`, `B_AdS = −6552.66` (NF02, `k(0) = 1`).
* The **unconditional** reading "n = 1 at Δ_g = 300/1000" is **withdrawn** (S1_S4, LARGEJ). D-strong (`Ω^U`/M02, B′ objects) is
  withdrawn (review P1).
* `ε`-dependence: multiply by `√(1+ε)` (§3.1).

### 6b. What is established (replaces §3a–§3c)

| ingredient | status | file |
|---|---|---|
| NF02 flat positivity (stripped), min `𝔉/N = 0.200` | [certified] (one `+−` corner `J ≥ 400`, `m < 1.633M`, float only) | `NOFORWARD_FUNCTIONAL.md` |
| implication per-state positivity ⇒ bound; `A = 3022.34`, `B_AdS = −6552.66`, `16L_A`, bounds | [D]; constants [DV] (review, this note) | review §0(a); `t1_band.LA_exact` |
| AdS IR identification `m_IR = 3/(2R_AdS)`: Result D = CHLPS-form flat bound of NF02, `3022.34 log(2Λ/3) − 5327.21` | [D] | this note §5.1 |
| light side: `−Ĉ₂′|_E = 16πG/(ν²+9/4)`, `Ĉ₃|_E = 0` (H7) | [D] (weakest input: §18.1 kernel) | `H7_I14_LIGHTSIDE.md` |
| validity / swap (I14, H16) | [D] (K1 KKSDZ [Q, loophole], K1′ [D-sketch]) | `H7_I14_LIGHTSIDE.md` |
| `Ĉ₂′` leading-Regge **diagonal** `h₋P_{l−1} + h₊P_{l+1}`, `h₀ ≡ 0`, `++ = +−` | [P given §18.1 (D)] + [DV] | `H7_I14_LIGHTSIDE.md` §18, §23 |
| `Ĉ₃` leading-Regge image `σ_h c₃(ν)P_l`, no de-shift ((HN₃)) | [P given kernel (D)] + [DV ≤ 4·10⁻¹³] | `S1_S4_FIXES.md` §1 |
| EO interference cancels; OO sign `N = tan²(πτ/2) ≥ 0` | [DV]; [NE + P given Q] | `PARITY_ODD_ROUTEA.md` |
| scalar-skeleton near-gap action at Δ_g = 100 (Theorem A) | [certified] (scalar model) | `CERTIFIED_NEARGAP_A.md` |
| scalar-skeleton hybrid at large `b̂`, `b̂ ∈ [91.2, 2Λ]`; `M₁ > 0` (diagonal) | [certified] (scalar model) | `ROUTEA_LARGE_B.md` |
| graviton near-gap corrections: `n = 1` for even `J ≤ 14`, odd `J ≤ 13` (300); `J ≤ 20` (1000) | [NE + M] (model M, K₂ + direct checks) | `S1_S4_FIXES.md` §2 |
| **per-state positivity fails** at 300 for `J ≥ 16–18` near the gap (`b̂ ≥ 30`), direct-validated; at 1000 from `J ≈ 32` | [NE + M] (300); [M] (1000) | `S1_S4_FIXES.md` §2.4–2.6 |
| real helicity off-diagonal is leading-Regge, `β²g(νβ)`, universal (26 rows, `β ≤ 0.11`); fixed-`b̂` direct ladder 1.96 vs 1.94 | [NE] | `LARGEJ_OFFDIAG.md` §2 |
| exact leading-Regge off-diagonal of `Ĉ₂′` (dual-line flip, `β²R(νβ)`, `∫tR = 47/4`); `Ĉ₃` none | [D] + [NE], **unrefereed, in progress** | `OFFDIAG_KERNEL.md` |
| failing band = `{b̂ ≥ b̂_c}` at all sampled `m/Λ`; `b̂_c ≈ 16 / 30 / 62–68` (`0.62–0.76Λ^{2/3}`); diagonal positive in the band | onset [NE + M] (300), [M] (1000, validated `g`), [M + I] (100 and beyond `β = 0.11`) | this note §2 |
| no `(ψ₂, ψ₃)` in NF02's family covers the band; window/bound trade-off (§4); asymptotic `|X|/F → 24|G₁|` | [M]; [D-sketch] | `LARGEJ_OFFDIAG.md` §3; this note §4 |
| flat analogue: D = 4 IR band (CHLPS (3.6)); NF02 flat-cutoff band `[1.35Λ^{2/3}, 2.56Λ]` | [Q]; [NE] | `LIT_HELICITY_MIXING.md`; this note §5.4 |
| **Theorem D_cond** (bound under A1–A2, V, L0–L3, P_<, H_band(ε)) | [D given the hypotheses]; **H_band [I]** | this note §3 |
| helicity-flip `B⁽³⁾`-type CFT sum rule / CEMZ matrix (would replace H_band) | not constructed (CLS-D §5 open) | `LIT_HELICITY_MIXING.md` §4 |

### 6c. Hypothesis table changes (relative to the review's §3b)

* **H13′** splits into **P_<** (graded as in §3.2) and **H_band** ([I], new, load-bearing).
* **H4** for `J ≥ 11`: computed, `n = 1` fails at `b̂ ≥ 30` (300) → absorbed into H_band.
* **New:** H_band (§3.1/§3.3).
* **New (definitional):** the choice of window functional (§4) is part of the statement.

---

## 7. Hand-over

**Files.**
* `t1_band.py`: band map.
  * Env: `T1_DG`, `T1_BMAX`, `T1_NQ` (geometric panels; 0 = r5's 64-node rule), `T1_NL` (linear panels), `T1_NB` (# `b̂` values),
    `T1_MFLS`, `T1_TAG`. Arguments: extra weight pickles from `cache/`.
  * Outputs: `cache/t1_band_<Dg><tag>.txt` and `t1_band_<Dg>_<name>.pkl` (per-state rows: `mfl J bh beta F++/N F+-/N X/sqrtNN lam kappa
    XLam2`).
  * `LA_exact(v, Λ)` gives the exact `16L_A`.
  * Runtimes: 300, `β ≤ 2`: ~4 min build; 1000, `β ≤ 0.5`: ~3 min; 100: ~2 min. Beyond `β ≈ 0.5` at Δ_g = 1000, mpmath's `₂F₁` does
    not converge.
  * Main runs: `t1_band_{100,300,1000}_ext.txt`; quadrature check `t1_band_300_q64.txt` vs `t1_band_300_qfine.txt`.
* `t2_lp.py`: r5's LP with new windows. Env `T2_DG`, `T2_BMAX`, `T2_DW`, `T2_OBJ`; `R3_D=12` for `d = 12`. Outputs
  `cache/t2_lp_<Dg>_b<β>_dw<dw>_<obj>[_d12].pkl/.txt`; ~15 s each. Note that the r5 pickle names do not record `d` or the `g`-scale:
  `r5_lp_300_b0.3_dw0.0.pkl` is the `g/2` design, and `r5_lp_1000_b0.16_dw0.0.pkl` is `d = 12`.
* `t3_q3check.py`: weights on q3's computed-correction blocks (`R3_EXTRA_M` for wedge states); outputs `cache/t3_q3check_<Dg>[_wedge].txt`;
  ~1 min.
* `t4_flatband.py`: flat-cutoff NF02 band (stdout only).
* `t1_queue.sh`: the original queue. Its 1000 line used a `d = 12` pickle and `β ≤ 2` and failed; it was rerun by hand as above.

**Reproduce:**
```
T1_DG=300 T1_BMAX=2.0 T1_NQ=10 T1_NL=40 T1_NB=70 T1_TAG=_ext python3 t1_band.py r5_lp_300_b0.12_dw0.1.pkl r5_lp_300_b0.2_dw0.1.pkl r5_lp_300_b0.25_dw0.0.pkl
T2_DG=100 T2_BMAX=0.3 python3 t2_lp.py
T3_DG=300 python3 t3_q3check.py r5_lp_300_b0.2_dw0.1.pkl
python3 t4_flatband.py
```

**Caveats.**
1. All band edges and `κ_crit` use r4's `g` (median of K₂ rows). OFFDIAG_KERNEL (in progress) reports that the physical real off-diagonal
   is `β²R` (`∫tR = 47/4`) and that r4's `g` contains a basis-shear piece. The swap-in point is `R5.G2` (and `G3 → 0`) in
   `t1_band`/`t2_lp`; edges scale as `|G₁|^{−1/3}`.
2. In the model the diagonal is the scalar hybrid (no graviton diagonal correction beyond q3's rows).
3. The 1000 edge rests on the universal model only, because q3's `J = 32` row is broken.
4. The twist convention is r5's `α ≥ Δ_g − 2`.

**Open, in value order.**
1. Rerun §2 and §4 with the exact off-diagonal kernel (OFFDIAG_KERNEL), especially (a) whether the band closes at some `β ~ 0.3–2` (a
   finite band would make H_band a statement about a bounded region), (b) the Δ_g = 100 tension of §2.3.
2. If OFFDIAG_KERNEL's cancelling (2,4) structure `S(κ*)` is adopted, redo the band from scratch. Theorem D_cond may then hold with
   `ε = 0` and no H_band at the leading-Regge level.
3. A `B⁽³⁾`-type helicity-flip CFT sum rule (CEMZ matrix) would turn H_band into a derived inequality (LIT §4).
4. Graviton diagonal corrections in the window `b̂ < b̂_w` at `m ≥ 1.3Λ` (P_< there is [M]).
5. Convert AKT's `ñ` normalisation to `t₄` for a numerical comparison.

**Rules kept.** No existing document edited, nothing committed, no subagents, ≤ 2 single-thread processes, < 1 GB each.
