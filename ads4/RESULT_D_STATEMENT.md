# Result D — consolidated statement (route A′ with `Ω^EE`; state of 20 Nov 2026)

> **Post-review status (REVIEW_EE_2026-11.md; supersedes §0.3, §2, §4, §5.2, §6 below wherever they differ).**
> No fatal error; constants, bounds and binding states reproduce. Corrections adopted:
> * "Every computed state passes" means **the EE block** only; the parity-odd even-J (OO) block is at leading order only (H_OO).
> * **New hypothesis H_edge:** the light-cone-edge window x ∈ (0.8, Λ/m] is held extrapolation for every computed state (cutting at
>   x = 0.6 moves J = 36/40/31 at 1000 by −0.43/−0.33/−0.23 N_s); not computed [I].
> * **K₂ beyond its convergence criterion (J²/m ≳ 2):** J = 28, 32 at Δ_g = 300 regraded to [M]/[I]; J = 48 at 1000 is at the K₂ limit.
> * **Δ_g = 300, J ≤ 22:** post-fix numbers use the superseded small-ν treatment; J = 22 is +0.123 with t* = 1.1 (10 % margin), not +0.15.
>   The local vs min-norm difference at 300 is 0.02–0.10 N_s (favourable direction), not ≤ 3·10⁻³. J = 28 with the independent thimble
>   repeat: +0.69 (not +1.45 ± 0.58).
> * **Chronology (selection):** the canonical gauge was adopted then dropped after heavy states failed; q² was replaced by EE after the
>   J = 16–20 failures. Any fixed functional gives a valid bound if positive, but J = 18–28 at 300 are in-sample evidence. Out-of-sample:
>   J = 32–48 and odd J = 19–31 at 1000, and the direct small-ν nodes at 300 — all pass. **Ω^EE (min-norm, NF02) is fixed from now on.**
> * **H_LJ is large:** it covers most near-gap spins (J ≳ 50 up to ≈ 344 at 1000; J ≥ 34 and all odd J ≥ 17 at 300), with tolerance
>   1–30 c₂ against raw corrections of 10²–10⁴ c₂ (oscillating parts integrate out), and is out of reach of present numerics for
>   J ≳ 150 (1000) / ≳ 100 (300).
> * Flat limit: [P given the H0 §17 polarisation map] + [NE]; H17 carried over from q² (structure part not recomputed for EE); |c₁| ≤ 1
>   has no argument (c₁ is a free parameter); AKT/novelty comparison softened.
>
> **Result D (conditional), as it can honestly be claimed.** Assume A1 (twist gap, all J), A2 (tree level), V (K1 with KKSDZ's loophole,
> K2, K3) [D], EE admissibility [D], H7 [D], L0 [D], H17 [model carried over], c₁ free. Then per-state positivity of the fixed functional
> Ω^EE implies √(t₄²+α_T²) ≤ 3.13·10⁻⁶ (Δ_g = 300) and ≤ 2.89·10⁻⁸ (Δ_g = 1000) [D; constants DV]. Positivity at Δ_g = 1000 is computed
> [NE + M] on the EE block for near-gap even J ≤ 40 with rows, J = 48 (K₂ limit) and odd J ≤ 31, with the edge window extrapolated;
> [M] elsewhere; OO block leading order only; J ≳ 50 near the gap is H_LJ [I]. At Δ_g = 300, computed for even J ≤ 24 and odd J ≤ 15
> (only J = 24 with the validated small-ν treatment); J ≥ 26 is model-level or open. Δ_g = 100: not established. No condition is a
> spectral assumption beyond A1; all are statements about explicit block integrals.


*20 Nov 2026, consolidation agent. Unrefereed, not committed. Written to replace the obsolete tables (§1–§6) of `RESULT_D.md`,
which still describe the withdrawn B′-uplift functional. No existing document or code file was edited. Sources (later ones win where
they differ): the `RESULT_D.md` update box (all bullets), `REVIEW_ROUTEA_2026-11.md`, `RESULT_D_CONDITIONAL.md`, `EE_FUNCTIONAL.md`,
`EE_FOUNDATIONS.md`, `EE_NEARGAP.md`, `EE_DG300_SMALLNU.md`, `EE_DG1000_LARGEJ.md`, `EE_DG100.md`, `SUBLEADING_REGGE.md`,
`SMOOTH_BOUND.md`, `FOLD_WINDOW_LARGEJ.md`, `H7_I14_LIGHTSIDE.md`, `LIT_HELICITY_MIXING.md`, `H0_SPINNING_FUNCTIONAL.md` §R.
§9 lists the few computations made for this note (odd `J = 19–31` at `Δ_g = 1000`, a phase-rotation precision check, and the
`Δ_g = 300`, `J = 32` small-`ν` chain).*

**Grades.** **[certified]** computer-assisted proof (interval/Arb); **[P]** proved; **[P given X]** proved from the named inputs;
**[D]** derived, not independently checked; **[D-sketch]**; **[Q]** quoted from the literature; **[DV]** derived and checked numerically;
**[NE]** numerical evidence; **[M]** model (an approximation whose accuracy is itself a hypothesis); **[I]** open/assumed;
**[A]** standing assumption of the setting. A chain of inputs carries the grade of its weakest link (review P4).

---

## 0. Headline (read this first)

1. **The implication is sound.** If the functional `Ω^EE` (§2) acts positively (per-state PSD) on every twist-allowed heavy conformal
   block, then **`√(t₄² + α_T²) ≤ 3.13·10⁻⁶` at `Δ_g = 300` and `≤ 2.89·10⁻⁸` at `Δ_g = 1000`** (tree level, AdS units), with
   `±0.4 % / ±0.04 %` from H17 and a linear-`t₄` light term with `|c₁| ≤ 1` (`±2.4 % / ±0.2 %` for `|c₁| ≤ 10`). The step
   "positivity ⇒ bound" is [D]; the constants are [DV] (reproduced by the hostile review).
2. **Result D is not a theorem at any `Δ_g`.** It is a conditional statement. The conditions are of two kinds:
   * *analytic inputs* (validity/swappability K1–K3 with KKSDZ's loophole, the exact graviton pole via EE's (HN), admissibility of `EE`,
     H17, the linear-`t₄` term `c₁`): graded [D], [D-sketch + M] and [I] (§4);
   * *per-state positivity* `P`: [NE + M] on the computed near-gap states, [M] in the Regge/large-`b̂` region, and a named hypothesis
     **`H_LJ(Δ_g)`** on near-gap states of large spin (§5).
3. **Strongest supported statement** (§6):
   * **`Δ_g = 1000`:** the bound `2.89·10⁻⁸` holds given the analytic inputs, the model hypotheses `H_M` (the near-gap chain is accurate on
     the computed states), `H_OO` (the uncorrected even-`J` parity-odd block) and `H_R` (beyond-leading-Regge corrections at `m ≳ 1.3Λ` stay below the curvature margin), and
     **`H_LJ(1000)`**: near the gap (`m ≲ 1.3Λ`), for even `J ≥ 50` (and the even `J = 26, 30, 34, 38, 42–46`, which have no rows) and odd `J ≥ 33`, the non-oscillating beyond-leading-Regge part of
     the graviton correction is above the `SMOOTH_BOUND` tolerance (uniform: `−14 / −5 / −1` `c₂` at `J = 64 / 100 / ≥ 200`; fold window
     `ν < 30`: `−3580 / −1170 / −205` at `J = 64 / 100 / 300`), and the `ν < 50` window is above `−1300 … −90` `c₂` (`J = 64 … 300`).
     Every computed state passes: even `J ≤ 40` and `48`, odd `J ≤ 31`, `m` from the gap to `5Λ`; minimum `+0.2021` (exact-skeleton grid;
     `+0.2009` for the hybrid state at the twist edge).
   * **`Δ_g = 300`:** the bound `3.13·10⁻⁶` holds given the same inputs **plus `H_LJ(300)`**, which is larger: near the gap, even
     `J ≥ 34` (and the uncomputed `J = 26, 30`), odd `J ≥ 17`, and the `ν < 8.3` window at `J = 28, 32` (failure only below
     `Γ_c ≈ −4000 c₂` there; thimble data `±10⁴`). Computed states pass: even `J ≤ 24`, odd `J ≤ 15`; `J = 28, 32` on central values;
     minimum `+0.2052` (`+0.2043` hybrid at the twist edge).
   * **`Δ_g = 100`:** not established. Main-protocol `n = 1` for `J ≤ 10` [NE + M], but the odd extra-row state is marginal and both
     transfer models fail at near-gap `J ≥ 18–24` [M].
4. **Nature of the remaining gap.** `H_M`, `H_R`, `H_LJ` and the even-`J` parity-odd block are statements about explicit conformal-block
   integrals (decidable by computation in principle), not assumptions about the CFT spectrum. The only physics-level inputs are A1 (twist
   gap for all single traces), A2 (tree level), and the KKSDZ Regge bound with its loophole. This is the main gain over the
   `q²` route, whose fallback needed a spectral hypothesis (`H_band`, `RESULT_D_CONDITIONAL.md`).
5. **Literature:** this is exactly a CHLPS-form flat bound of the forward-free NF02 functional with `m_IR = 3/(2R_AdS)`,
   `|ĝ₃|²Λ⁸ ≤ 3022.34 log(2Λ/3) − 5327.21` [D]; about 10× weaker than the transplanted (asserted, not derived) CHLPS optimum, and
   parametrically `√(ln Δ_g)/Δ_g⁴`, stronger than AKT's `ln Δ_g/Δ_g⁴`. No published paper derives explicit-constant `t₄` bounds from
   spinning CFT sum rules (§6.4).

---

## 1. Setting and conventions

* **A1 (spectrum) [A].** `𝒯` is a unitary CFT₃ (parity-even or parity-violating) with a unique stress tensor `T`. **Every exchanged
  single-trace operator other than `T` has twist `Δ − J ≥ Δ_g`, for all `J` (including `J = 0, 1, 2`).** No other light single
  traces. This is the convention every scan uses (`omegau/common.state`). `RESULT_D.md` §1's "`Δ − J ≥ Δ_g − 2` for `J > 2`" is wrong
  (review M1). A1 excludes all known holographic CFT₃s (ABJM included), which have light Kaluza–Klein towers.
* **A2 (tree level) [A].** `O(1/c_T)` truncation, `c_T ≫ Δ_g²`; loops dropped.
* Heavy states `(Δ, J)`: `α = Δ − J − 2 ≥ Δ_g − 2`, `m² = α(α + 2J + 1)`, `Λ := √((Δ_g − 2)(Δ_g − 1))` (`= 98.499 / 298.50 / 998.50`),
  `b̂ = (2J+1)Λ/m`, `cosh β = 1 + (2J+1)²/(2m²)`, `x = ν/m`, `p = ν/Λ`. The gap state of spin `J` has `m_gap = √((Δ_g−2)(Δ_g−1+2J))`.
* Couplings: `⟨TTT⟩` in the basis `(c_T, t₄, α_T)`; in the bulk `t₄² + α_T² = 57600|ĝ₃|²/L⁸`, so `√(t₄²+α_T²) = 240|ĝ₃|` in AdS units
  (`DICTIONARY_T4.md`, reproduced by the foundation referee).
* `n`: positivity threshold (heavy states positive above `Λ/n`); "`n = 1`" means per-state positivity from the gap.

## 2. The functional `Ω^EE`

**Definition** (`EE_FUNCTIONAL.md` §1–§3, `EE_FOUNDATIONS.md` §1–§3, `H0` §R items 1–4, 7; `S1_S4_FIXES` for `Ĉ₃`):

`Ω^EE := ∫₀¹ dp ψ₂(p) Λ⁴ Ĉ′_{EE,Λp} + ∫₀¹ dp ψ₃(p) ½Λ⁶ Ĉ_{3,Λp}`, with

| ingredient | definition | grade |
|---|---|---|
| building block | CLS-D's light-ray-commutator functional `Ψ_{k,ν,j′}[f G]` on `⟨TTTT⟩` (arXiv:2311.04271, eq. `spinningfunctional_definition`) | [Q] |
| `k = 2` part | `j′ = 0`, `f = v′ − u′` (unsubtracted), dual (2,4) structure **`EE = [V_{0,42}²(H₂₄ − 2V_{2,04}V_{4,02})/(X₂₀X₄₀)]²`** (`= V₀⁴(z₂′·z₄′)²` at `x₀ = 0`) | definition; admissible [D] |
| de-shift | `Ĉ′_EE := 𝒜_h⁻¹[Ψ_EE] = c(ν)⁻¹𝒜_g⁻¹[Ψ_EE]`: EE's leading-Regge kernel is `h₋ = g·c(ν+i)`, `h₊ = g̃·c(ν−i)` (CHMRS's scalar `g`), `c(ν) = 32(4ν²+1)(4ν²+9)/((4ν²+25)(4ν²+81))` | [P given H7 §18.1 kernel (D)] |
| smeared form | `∫ψ₂Λ⁴Ĉ′_EE := Re∫dν X(ν)Ψ̂_{EE,ν}` (normalisations as `H0` §R.3–5) (`Ψ̂` forward-normalised by `h₋+h₊`), `X` the solution of the generalised Abel equation, solved in the χ-normalisation `C[aY] − i tanh t·S[bY] = C[ψ₂]/cosh t`, `X = aY` (bounded operator; `X` vanishes automatically at the real zeros `ν = 1.171513, 6.642858` of `h₋+h₊`) | [D] + [NE] (residual `≤ 3·10⁻¹⁵`) |
| **gauge** | the Abel solution has a 2-dim null space (index `w = 5`), supported on `ν ≲ 18`. **Convention: min-norm gauge** (correction orthogonal to the null space). | definition (see §8 item 1) |
| `k = 3` part | `Ĉ₃ := Re[Ψ_{3,ν}[o_H·q]/N̂₃(ν)]`, `j′ = 1`, `f = −1` (unsubtracted), **no de-shift** (its spinning leading-Regge image is exactly `σ_h c₃(ν)P_l`) | [P given kernel reading (D)] + [DV `4·10⁻¹³`] |
| reality | `Re` of the complex sum rule (both parts separately valid for real OPE data); phase convention `logK → logK − iπ(Δ+J)` on every record | definition; legitimacy [D]; strip [P given NE] |
| weights | NF02: `ψ_k(p) = (1−p)⁵ p X_k(p)/k₀`, degree-8 rational polynomials (`robust/functionals_noforward.json`), `k₀ = 27.48509817`, `k(0) := −2∫ψ₂p⁶ + 2∫ψ₃p⁴ = 1` | definition |

**Properties.**
* **Flat limit:** at the `D = 4` scattering crystal `EE ≡ q²` identically [P, sympy, `x` symbolic], so `Ω^EE`'s bulk-point limit is
  the stripped NF02 functional `∫ψ₂B⁽¹⁾₂ + ∫ψ₃B⁽¹⁾₃`; NF02's flat positivity certificate applies (min `𝔉/N = 0.200`) [certified], except
  the `+−`, `J ≥ 400`, `1.5Λ ≤ m < 1.633Λ` corner, which is float-checked only (empty at `Δ_g = 100, 300`; at 1000 inside the large-`b̂`
  window, `b̂ ≈ 490–1020`, where the curvature term is certified positive) [`NOFORWARD_FUNCTIONAL.md`, `ROUTEA_LARGE_B.md`].
* **Leading-Regge heavy action:** helicity-diagonal, `P_{−½+iν}(cosh β)` on both helicities; **the leading-Regge helicity off-diagonal
  vanishes identically for every `(ν, β)`** [P] (the `q²` obstruction of `LARGEJ_OFFDIAG`/`OFFDIAG_KERNEL` is absent).
* **What the chains compute (review M7):** the local first-order de-shift `(1 ± u_EE)/2` at nodes `ν ≥ 15–25`, and the non-local `X` below
  `ν = 25` where direct data exist (z4 chain at 300; EE_DG100). The difference is `≤ 3·10⁻³ N_s` at 300 and negligible at 1000.

## 3. Light side and the bound

**Light side** (tree level; `H7_I14_LIGHTSIDE.md`, `EE_FOUNDATIONS.md` §2, `H0` §R.8):
* Both functionals are unsubtracted, so every `[TT]` double trace sits at a double zero: parity-even (`2 sin²(π(τ−10)/2)`) [D given Q],
  parity-odd (`2 cos²(πτ/2)` at odd `τ`) [P given Q + P-sketch]. Contact terms drop out [D given Q]. The tree-level light side is the
  identity plus the single `T` block.
* **`−Ĉ′_{EE,ν}|_E = 16πG/(ν² + 9/4)` exactly** — (HN) holds in its original CHMRS form for `EE` (`N̂ = c(ν)`), so no rational mismatch
  factor and no double pole [P given Lemma R, Q3, Q4, D2, §18.1 kernel (D, two inputs)] → overall **[D]**.
* **`Ĉ₃|_E = 0`** (bare light-ray commutator, Einstein one-point matrix = identity) [P given Q2, Q8, F]; `Ĉ₃` has no spin-2 Regge moment
  (≤ 7·10⁻⁴) [D + NE].
* `|ĝ₃|²` term: flat images `−B₂ ∋ 2πG|ĝ₃|²p⁶`, `−B₃ ∋ −2πG|ĝ₃|²p⁴`, coefficient `k(0)(1 + δ)` with
  `δk(0)_EE = −7.245·10⁻³ / −7.666·10⁻⁴ / −6.830·10⁻⁵` (h1 convention) and the H17 remainder.
* Possible `E × t₄` cross term `16πG c₁ t₄` (ν-independent): not computed [I] (review M3).

**Sum rule and bound.** With heavy ≥ 0: `Σ_heavy = 16L_A − k(0)(1+δ)|ĝ₃|²Λ⁸ ≥ 0`, `L_A = ∫₀¹ψ₂(t)dt/(t² + 9/(4Λ²))`; `G` cancels.

**`√(t₄² + α_T²) ≤ 240·√(16L_A(Λ)/(k(0)(1+δ)))/Λ⁴`**, `A = 16ψ₂′(0) = 3022.34`, `B_AdS = −6552.66` (`16L_A = A log Λ + B_AdS + O(1/Λ)`).

| `Δ_g` | `Λ` | `16L_A` (exact) | `A log Λ + B_AdS` | **bound** (exact `p`-integral, `k = 1`) | `δk(0)_EE` factor | H17 band (on bound) | `c₁` (per unit, on bound) | status of `P` |
|---|---|---|---|---|---|---|---|---|
| 100 | 98.499 | 7508.5 | 7320.0 | **2.21·10⁻⁴** | ×1.0036 (2.22·10⁻⁴) | ±2 % | ≈ ∓1.7 % | **not established** (§5.3) |
| 300 | 298.50 | 10734.6 | 10670.9 | **3.13·10⁻⁶** | ×1.0004 | ±0.2 % | ∓0.22 % | conditional on `H_M, H_R, H_LJ(300)` (§5.2) |
| 1000 | 998.50 | 14339.6 | 14320.4 | **2.89·10⁻⁸** | ×1.00003 | ±0.02 % | ∓0.017 % | conditional on `H_M, H_R, H_LJ(1000)` (§5.1) |

`A`, `B_AdS` and the bounds are gauge-independent (Lemma N) and identical for `q²` and `EE` (same pole, same `ψ₂`).
The log form underestimates the exact integral by 0.6 % at 300; the table uses the exact integral.

---

## 4. Hypothesis table

| # | hypothesis | grade | supporting doc | comment / what would settle it |
|---|---|---|---|---|
| **A1** | twist gap `Δ − J ≥ Δ_g` for every single trace `≠ T`, all `J`; no other light single traces | [A] | review M1, §2 | setting; excludes ABJM-type spectra |
| **A2** | tree level, `c_T ≫ Δ_g²` | [A] | review §2 | setting |
| **K1** | KKSDZ Regge bound `J₀, J₀^{dDisc} ≤ 1` for spinning correlators; needed `J₀ < 2` (one unit of margin) | [Q], **with KKSDZ's stated loophole** at the first/second-sheet boundary (their §4.2.6) | `H7_I14` §2 (Q6), §7.1 | the loophole is a technical gap in the cited theorem, not closed here |
| K1′ | uniformity of K1 in the angle towards the light-cone edges | [D-sketch] (not written out) | review V, P4 | write out |
| **K2** | kernel power counting at every boundary (Regge corner [Q]; `1 ≈ 3` endpoint and edges [D]; `η → ∞` decay `η^{−½}/log η` of the total de-shifted kernel) | [Q]/[D]; `η → ∞`: [P given K1′] (Prop. K applies verbatim to `EE`: the de-shifted total kernel is `∫ψ₂P_l`); spinning components [D] | `H7_I14` §7, §14; `EE_FOUNDATIONS` §2.4 | — |
| **K3** | per-family Rindler positivity (`P_O` commutes with modular flow) ⇒ Lemma S (swap) | [D] | `H7_I14` §7.2 | elementary; not quoted |
| **V** (validity + swap, I14) | `Ω^EE` annihilates the non-perturbative `⟨TTTT⟩` and commutes with the OPE sum | (a) [P given K1, K2]; (b) [P given K1, K1′, K2, K3] → **overall [D]** | `H7_I14` §7; `EE_FOUNDATIONS` §1.4, §2.3–2.4 | structure-independent in `d = 3` (complete gauge fixing; dual structure is a constant tensor) [D]; for `EE` use the contracted integrand near `x₀₂² = 0` (double pole at the light-cone edge) |
| **Adm** | `EE` is an admissible `j′ = 0` dual structure | **[D]** (weights, transversality, membership in the 13-dim `j′ = 0` space [P]; CLS-D's printed list is a saddle-quotient basis [Q + DV]) | `EE_FOUNDATIONS` §1 | was [D-sketch] in `EE_FUNCTIONAL` |
| **(HN) / H7** | `−Ĉ′_EE|_E = 16πG/(ν²+9/4)` exactly; `Ĉ₃|_E = 0` | (HN) for `EE` in original CHMRS form [P given §18.1 kernel (D)]; pole [P given Lemma R (R4 D-sketch, R5 assumed), Q3, Q4, D2 (D), §18.1 (D)] → **[D]**; `Ĉ₃` part [P given Q2, Q8, F] | `H7_I14` §3–§5, §12, §18–§23; `EE_FOUNDATIONS` §2.1–2.2 | one independent small-`ν` check of `h₋ + h₊` (review M5); a defect would move `B_AdS`, not `A` [M] |
| **L0** | light content: `[TT]` double zeros (even [D given Q]; odd [P given Q + P-sketch]); contacts drop out | [D] | `H0` §15; `PARITY_ODD_ROUTEA` §12 | — |
| **Flat** | flat limit = stripped NF02, certificate applies | [P] (crystal identity `EE ≡ q²`) + [Q + DV] (CLS-D crystal) + [certified] (NF02, one float corner) | `EE_FUNCTIONAL` §1; `NOFORWARD_FUNCTIONAL` | — |
| **Re** | `Re` prescription; β-phase strip `logK − iπ(Δ+J)` | definition; legitimacy [D] (conditions (a)–(c) hold [DV/NE]); strip [P given NE] (verified on 99 shadow entries and > 7000 records); **no `Δ_g ≥ 300` number depends on the strip** (even-integer `Δ`, same-`J` ratios) | review "checked"; `PARITY_ODD_ROUTEA` §12; `EE_NEARGAP` §1 | load-bearing as a definition (treating the imaginary mixing as real breaks `n = 1`) |
| **Gauge / Lemma N** | min-norm gauge (convention); null directions invisible to `A`, `B_AdS`, the flat limit, the graviton pole and the leading-Regge action | gauge: definition; Lemma N [P given NE-observed regularity] + [DV `10⁻¹⁵`]; gauge dependence beyond leading Regge: `≤ 10⁻⁴` at 1000 [NE + M]; at 300 the canonical gauge is **not** controlled [I] | `EE_FOUNDATIONS` §2.5, §3; `EE_NEARGAP` §4; `EE_DG300_SMALLNU` §6 | any gauge gives a valid functional (each solves the same Abel equation), so the choice is a design decision; see §8 item 1 |
| **H17 / L2** | `|ĝ₃|²` light coefficient `k(0)(1+δ)`: `|δ| ≲ 4.4 % / 0.48 % / 0.043 %` (bound ±2 % / ±0.2 % / ±0.02 %), incl. `δk(0)_EE` | [D-sketch + M] (heuristic `|a_j| = D × 36` model; sign unknown) | `H0` §18; review M4 | spinning `R³` exchange to `O(r²)` |
| **c₁ / L3** | any `E × t₄` light term `16πG c₁t₄` has `|c₁| ≲ 10` | **[I]** (not computed); effect on bound `∓0.22 % · c₁` (300), `∓0.017 % · c₁` (1000) | review M3 | Regge residue of the `E × R³` `T`-exchange (Lemma R with `n = 2`) |
| **H_M** | the near-gap chain (model M) is accurate on the computed states: base = exact scalar skeleton or hybrid (flat NF02 + exact leading-Regge curvature); graviton correction `Γ` from `K₂` rows (`Δ_ref = 2·10⁴` reference), direct RQMC/contour/thimble nodes where computed; shape fit / spline in `x`, linear in `1/m` between rows; Hampel despike; curvature double count removed | [M], validated [NE]: direct vs `K₂` at `(336,18)` (+0.466 vs +0.450), `(336,24)` (state values ≤ 1 %), `(1040,40)` (≤ 0.05 `N_s`), `(1600,20)`; ensemble spreads | review S3; `EE_NEARGAP`; `EE_DG1000_LARGEJ` §3–§4; `SMOOTH_BOUND` §4; `FOLD_WINDOW_LARGEJ` §4 | `K₂` usable only for `J²/m ≲ 1.5` |
| **H_OO** | the parity-odd **even-`J`** block (`X^O_same`, couples like `++`) stays positive beyond leading order | leading order `N = tan²(πτ/2) ≥ 0` × `++` density [NE + P given Q]; `O(1/m²)` correction never computed → [M] | review M9; `PARITY_ODD_ROUTEA` | one `p6`-type run at even `J` |
| **H_R** | Regge/large-`b̂` region: beyond-leading-Regge graviton corrections at `m ≳ 1.3Λ` (and everywhere at `b̂ ≳ 20` above the near-gap sets) stay below the leading-Regge + curvature margin | [M]: leading-Regge EE model positive on the whole band (`β ≤ 2.2`), min `0.2263 / 0.2056 / 0.2056`; computed corrections at `m ≥ 1.3Λ`: `≤ 0.005 N_s` (1000, `J ≤ 32`), `≤ 0.02` (1000, `J ≤ 40`), `−0.04 … −0.49` with all states passing (300, `J = 12–28`); scalar exact remainder `|S|/F_hyb ≤ 122/m²` [NE]; `F_hyb > 0` [certified] on `b̂ ∈ [91.2, 2Λ]` (scalar hybrid) | `EE_FUNCTIONAL` §4; `EE_NEARGAP` §6; `ROUTEA_LARGE_B`; review S2, M2 | the spinning `O(r)` Regge kernel (review S2 fix 1) |
| **H_LJ(Δ_g)** | near the gap, at large spin beyond the computed set, the non-oscillating beyond-leading-Regge correction (and the small-`ν` window) is above the stated tolerance | **[I]** (precise form and evidence in §5) | `SUBLEADING_REGGE`; `SMOOTH_BOUND`; `FOLD_WINDOW_LARGEJ`; `EE_DG300_SMALLNU` §8; `EE_DG1000_LARGEJ` §5 | analytic bound on the smooth part (proof of `SMOOTH_BOUND` §1.2 + edge leakage), or a method without the state-side sign problem at `J²/m ≫ 1` |
| H_cont | continuity in `m` between grid masses and at the twist edge | [NE] (hybrid at the edge: 0.2009 at 1000, 0.2043 at 300) | review M8; `EE_DG1000_LARGEJ` §4 | — |
| H3, H9a/b, H10, H14, H15 | odd-`k` moments; forward/smeared contact content; EFT power counting; object identity; `|ĝ₃|²` absent from `Ψ_{5,6,7}` | **not needed** (only unsubtracted `k = 2, 3`; contacts drop out at tree level; object fixed by definition) | review §3b | — |

---

## 5. Per-state positivity: coverage maps

Per-state test: the EE 2×2 (`++`, `+−`) ⊕ OO 1×1 action matrix of each heavy block is PSD; reported as min `λ/N_s`. "Near the gap"
means `m ≲ 1.3Λ`. Labels: **[NE + M]** = graviton correction computed (model M, `H_M`); **[M]** = leading-Regge + curvature model only
(κ-uniform model of `SUBLEADING_REGGE` near the gap; hybrid/band model elsewhere); **[I]** = neither, or the computed data do not
resolve the tolerance. The minimum over all states is always the parity-odd `+−`, `J = 5` state at the gap.

### 5.1 `Δ_g = 1000`

| region (`J`; `m`) | status | value / evidence | source |
|---|---|---|---|
| even `J = 0–24`, `28`, `32`; odd `J = 5–15`; `m/Λ = 1.01 … 5` | **[NE + M]** | min `+0.2021` (`+−`, `J = 5`, `1.01Λ`); diagonal corrections `≤ 0.09`, off-diagonal `≤ 0.001 N_s`; gauge-independent to `10⁻⁴` | `EE_NEARGAP` §2; `SMOOTH_BOUND` §4 (curvature fix) |
| even `J = 32, 36, 40` (own-mass rows), odd `J = 17` | **[NE + M]** | `+0.904 / +0.967 / +0.979 / +0.839`; `J = 40` direct sub-row agrees to `≤ 0.05 N_s` | `EE_DG1000_LARGEJ` §0, §4; `SMOOTH_BOUND` §4 |
| even `J = 48` (fixed-engine rows) | **[NE + M]** | `+0.932` (variants `+0.73 … +0.99`); old "offset" was engine precision loss | `FOLD_WINDOW_LARGEJ` §4 |
| **odd `J = 19–31`** (rows `Δ = 1000 + J`, 35 nodes) | **[NE + M]** (this note, §9.1) | `+0.916 / +0.942 / +0.887 / +0.918 / +0.954 / +0.979 / +0.995` (`J = 19 … 31`); every mass positive; corrections `≤ 0.016 N_s` except `J = 31` (`+0.21`, stabilising) | §9.1 |
| even `J = 26, 30, 34, 38, 42–46` near the gap | [M] (interpolation in `J`; no rows) | neighbours positive with margin `≥ 0.9` | — |
| gap-mass hybrids `J ≤ 15` not on the exact-skeleton grid | [M] | min `+0.2009` (`J = 5`, `1.0048Λ`), `J = 4` `+0.204` | `EE_DG1000_LARGEJ` §4 |
| near gap, even `J = 64`, `100` | [M] + partial [NE] | κ-model `+0.993 / +0.981`; fold window: EE `−202 ± 210` (tol 3580) ✓, C3H `±700` ✓ at 2σ (`J = 64`); EE `−378 ± 243` (tol 1170) ✓, C3H only ~1σ (`J = 100`) | `SUBLEADING_REGGE` §0.3; `FOLD_WINDOW_LARGEJ` §0.2 |
| near gap, even `J ≥ 50` (other), odd `J ≥ 33`, up to `J = 300` | **[M]** | κ-uniform model `≥ +0.90` (`J = 300`: `+0.902`); model dependence `< 10⁻³` | `SUBLEADING_REGGE` §0.3 |
| near gap, `J > 300` (twist-allowed to `J ≈ 344` at `m ≤ 1.3Λ`); wedge `m ≤ 1.5Λ` to the twist cap (`J ≤ 624`) | [M] | band/wedge scans: `J > 20` `≥ 0.21` (hybrid), band model `≥ 0.90` (`J ≤ 300`) | review §4 item 2; `EE_DG1000_LARGEJ` §5 |
| `m ≥ 1.3Λ`, `J ≤ 40` | [NE + M] (rows interpolated/held in `1/m`) | corrections `≤ 0.02 N_s` (`≤ 0.005` for `J ≤ 32`) | `EE_NEARGAP` §6; `EE_DG1000_LARGEJ` §0.3 |
| `m ≥ 1.3Λ`, larger `J`; all `b̂ ≥ 20` above the near-gap sets; `b̂ > 2Λ`; `b̂ < 20` at `m > 1.5Λ` | **[M]** (`H_R`) | leading-Regge EE band model positive on the whole band (min 0.2056; large `b̂` at 0.93–1.00); `F_hyb > 0` [certified] on `b̂ ∈ [274.1, 1997]` and `[91.2, 2Λ]`; `b̂ < 20` at `m ∈ [1.6, 30]Λ` `≥ 0.48` | `EE_FUNCTIONAL` §4; `ROUTEA_LARGE_B`; review §4 items 2–3 |
| parity-odd even-`J` (OO) block, all `J` | [NE] at leading order + [M] (`H_OO`) | `N = tan²(πτ/2) ≥ 0` times the `++` density | review M9 |

**The [I] residue at 1000 = `H_LJ(1000)`:** *for every twist-allowed state with `m ≤ 1.3Λ` and even `J ≥ 50` (the fold items at `J = 64, 100` EE
are measured), even `J ∈ {26, 30, 34, 38, 42, 44, 46}` (no rows; interpolation between positive neighbours) or odd `J ≥ 33`, the non-oscillating part `r` of the graviton beyond-leading-Regge correction
`Γ − (κ-image)` satisfies `r ≥ −T(J)`*, with the `SMOOTH_BOUND` tolerances (`c₂` units, worse helicity, the value that would remove
`0.9λ` from the gap state):

| `J` | 48 | 64 | 100 | 150 | 200 | 300 |
|---|---|---|---|---|---|---|
| uniform smooth `T` | 29 | 14 | 5 | 2 | 1 | 1 |
| fold-localised (`ν < 30`) | 7750 | 3580 | 1170 | 490 | 306 | 205 |
| light-cone-edge-shaped (`A_E`) | 1.2·10⁴ | 2900 | 390 | 80 | 32 | 12 |
| uniform in the `ν < 50` window | −2800 | −1300 | ≈ −450 | ≈ −220 | −130 | −90 |

*Evidence for:* (i) the exact scalar skeleton at `J = 40, 48, 64, 100` has `|rem/N_s| ≤ 4·10⁻⁴`, purely oscillating,
`10²–10⁴` below tolerance [NE, exact]; (ii) graviton bulk smooth part `|r| ≤ 6–12 c₂` at `J ≤ 40`, `r = 0 ± 5` at `J = 48` [NE];
(iii) fold window inside tolerance at `J = 64` (EE, C3H) and `J = 100` (EE) [NE]; (iv) saddle argument: every order of the crystal-saddle
expansion oscillates, so smooth terms need the fold (`ν* ≈ 11.3`, `J`-independent) or the light-cone edge [D-sketch]; (v) every
computed central value is consistent with zero and none indicates a violation. *Against / not known:* the bulk smooth term at
`J > 48` is untestable with present methods (`K₂` unconverged at `(1064, 64)`, `J²/m = 3.9`); C3H at `J = 100` only at ~1σ; `J ≥ 150`
undecidable with the present direct method (state-side sign problem: `(1200, 200)` element errors 6–8.5 %); no graviton data at
`ν < 10` or `x > 0.8`; the `ν < 50` window is decisive for `J ≳ 100` and uncomputed.

### 5.2 `Δ_g = 300`

| region (`J`; `m`) | status | value / evidence | source |
|---|---|---|---|
| even `J ≤ 22`, odd `J ≤ 15`; `m/Λ = 1.02 … 5` | **[NE + M]** in every variant | binding `+0.2052` (`+−`, `J = 5`, `1.02Λ`); worst even `J = 22` `+0.15` (convention A) / `+0.21–0.36` (other variants); odd `J = 15` `+0.707`, odd corrections grow `0.036 → 0.24 N_s` from `J = 9` to `15` | `EE_NEARGAP` §2, §4 |
| even `J = 24` | **[NE + M]** | direct small-`ν` nodes: `+0.82` after the curvature fix (min at 1.5Λ); 1.1Λ state `+1.06 ± 0.14`; fails only if `Γ_c ≲ −8900` below `ν = 8.3` (measured `+350 ± 4900`) | `EE_DG300_SMALLNU` §6; `SMOOTH_BOUND` §4 |
| even `J = 28` | [NE + M] on central values; **[I]** at 1σ in the `ν < 8.3` window | `+0.81` (curvature-fixed); 1.1Λ `+1.45 ± 0.58`; fails at uniform `Γ_c ≈ −4000` below `ν = 8.3`; thimble data `±10⁴` | `EE_DG300_SMALLNU` §6 |
| **even `J = 32`** (`Δ = 336` row, direct n3 nodes `ν ≥ 8.3`) | [NE + M] on central values (this note, §9.3); **[I]** in the `ν < 8.3` window and at `m ∈ [1.102, 1.15)Λ` (not sampled) | `+0.657` at 1.15Λ (Γ held below `ν = 8.3`); `+0.738` (zero below 8.3); fails at uniform `Γ_c ≈ −4200` below `ν = 8.3`; no thimble data | §9.3 |
| even `J = 26, 30` | [M] (interpolation; no rows) | neighbours positive | — |
| near gap, even `J = 64` | [M] + weak [NE] | κ-model `+0.968`; fold window C2 `−51 ± 43` vs tolerance 70: inside at 1σ, not at 2σ | `SUBLEADING_REGGE`; `FOLD_WINDOW_LARGEJ` §0.2 |
| near gap, even `J ≥ 34` (largest near-gap `b̂ ≈ 159` at `J ≈ 103`), odd `J ≥ 17` | **[M]** | κ-uniform model `≥ +0.92` for `J ≤ 150` | `SUBLEADING_REGGE` §0.3 |
| `m ≥ 1.3Λ`, `J = 12–28` | [NE + M] | corrections `−0.04 … −0.49 N_s`, all states pass (worst `+0.47` before / `+0.81` after the curvature fix) | `EE_NEARGAP` §6; `SMOOTH_BOUND` §4 |
| `m ≥ 1.3Λ`, larger `J`; large `b̂`; wedge; `b̂ < 20` at `m > 2Λ` | **[M]** (`H_R`) | as at 1000: band min `0.2056`; `F_hyb > 0` [certified] on `b̂ ∈ [162.6, 597]` and `[91.2, 2Λ]` (scalar hybrid); wedge `J > 20` `≥ 0.27` | `EE_FUNCTIONAL` §4; `ROUTEA_LARGE_B`; review §4 |
| parity-odd even-`J` block | [NE] LO + [M] | as at 1000 | review M9 |

**The [I] residue at 300 = `H_LJ(300)`:** *for every twist-allowed state with `m ≤ 1.3Λ` and even `J ≥ 34` (also `J = 26, 30`) or odd
`J ≥ 17`, the non-oscillating beyond-leading-Regge correction satisfies `r ≥ −T(J)`, and for `J = 28, 32` (and all larger `J`) the
correction in the window `ν < 8.3` is above the failure threshold.* Tolerances (`c₂`):

| `J` | 24 | 28 | 32 | 40 | 64 | 100 | 150 |
|---|---|---|---|---|---|---|---|
| uniform smooth `T` | 25 | — | — | 8 | 3 | 2 | 1 |
| fold-localised (`ν < 30`) | 653 | 440 | — | 186 | 70 | 36 | 27 |
| uniform below `ν = 8.3` (computed / [M]) | −8900 | −4000 | −4200 (§9.3) | ~ −10³ | ~ −3·10² | ~ −10² (`J = 103`) | — |

*Evidence for:* the exact scalar remainder at `J = 24, 40, 64, 100` is `≤ 0.08 c₂`-equivalent, 320–6000× below tolerance [NE, exact];
graviton bulk smooth part `|r| ≤ 20–36` at `J = 24–32` and fold-window `r = +3 … +86` (`≥ −52` at 2σ) at `J = 24, 28` [NE];
direct `Γ` on `ν = 10–33` is `O(10²)`, mostly positive, and does not grow from `J = 24` to `28`. *Against / not known:* the fold
tolerance falls below the measured 2σ range already at `J = 64` (`−51 ± 43` vs 70) and is far below the errors at `J = 100, 150`
(`±150–460` vs 36/27); the `ν < 8.3` window threshold falls like `N_s ∝ b̂⁻³` below thimble resolution for `J ≳ 40`; odd `J ≥ 17` has
no rows, and the odd correction was still growing at `J = 15`; the `Δ_g = 100` transfer models fail at `J ≳ 18–24`
(a different gap, but the same `J²/m` diagonal growth). **`H_LJ(300)` is materially weaker-supported than `H_LJ(1000)`.**

### 5.3 `Δ_g = 100` (status only)

* Main protocol (non-local min-norm weight; direct values at `x = 0.05` thimbles, `0.1–0.15` contour, `0.2–0.3` direct): **`n = 1`
  on even `J ≤ 10`, odd `J = 5, 7, 9`**, binding odd `J = 5` at 1.05Λ, `0.2286` (MC `0.229 ± 0.002`); `x < 0.075` may be scaled by
  `t ∈ [−0.23, 8.8]` [NE + M].
* Not robust: odd extra-row state `J = 9`, 1.1Λ: `+0.03` (MC `0.005 ± 0.044`; `−0.09` with zero extrapolation) [I]; even extra rows
  (`K₂` only, unstable `a_h`) give `J = 10`, 1.28Λ `−0.12` [M, weak]; near-gap `J ≥ 16`: no rows at `m ≈ 110–130`, and **both transfer
  models fail at `J ≥ 18–24`** (diagonal correction `∝ J²/m`) [M]; canonical gauge unusable.
* **Result D at `Δ_g = 100` is not established**; the nominal `2.21·10⁻⁴` (`2.22·10⁻⁴` with `δk(0)`) is conditional on `n = 1` over all
  states, which the model evidence at `J ≳ 18` contradicts. [EE_DG100.md]

---

## 6. Statements

### 6.1 The implication [D]

*Let `𝒯` satisfy A1, A2. Assume V, Adm, H7 (L1), L0, H17 (L2), L3 (`|c₁| ≤ c̄`). If every twist-allowed heavy block has a PSD per-state
action matrix under `Ω^EE` (`P`), then*

`√(t₄² + α_T²) ≤ 240·√(16L_A(Λ)/(k(0)(1+δ)))/Λ⁴ · (1 + O(c̄·Λ⁻²))`,

*i.e. `≤ 3.13·10⁻⁶` (`Δ_g = 300`) and `≤ 2.89·10⁻⁸` (`Δ_g = 1000`), `±0.4 % / ±0.04 %` for `c̄ = 1`.* Proof: heavy `≥ 0` ⇒
`−Ω_light ≥ 0` ⇒ `k(0)(1+δ)|ĝ₃|²Λ⁸ ≤ 16L_A` (review §0(a), reproduced). The constants are [DV].

### 6.2 The honest headline (strongest statement the evidence supports)

> **Result D (conditional).** Under A1, A2 and the analytic inputs V [D] (K1 with KKSDZ's loophole), Adm [D], H7 [D], L0 [D],
> H17 [D-sketch + M] and L3 [I],
> * **at `Δ_g = 1000`**, `√(t₄²+α_T²) ≤ 2.89·10⁻⁸` holds **given `H_M`, `H_OO`, `H_R` and `H_LJ(1000)`**; per-state positivity is
>   verified [NE + M] on every computed state (even `J ≤ 40` and `48`; odd `J ≤ 31`; `m` from the gap to `5Λ`; min `+0.2021`) and in the
>   leading-Regge + curvature model [M] everywhere else (min `+0.2009` at the twist edge);
> * **at `Δ_g = 300`**, `√(t₄²+α_T²) ≤ 3.13·10⁻⁶` holds **given the same plus `H_LJ(300)`**, whose region is larger (near-gap even `J ≥ 34`,
>   odd `J ≥ 17`, the `ν < 8.3` window at `J = 28, 32`) and whose evidence is weaker (fold item at `J = 64` at the tolerance edge);
> * **at `Δ_g = 100`** it is **not established** (model-level failures at near-gap `J ≳ 18`);
> * at any other `Δ_g` nothing has been computed ("`Δ_g ≥ 300`" is an extrapolation, review M10).

What would make it a theorem (in value order): (1) an analytic bound on the non-oscillating beyond-Regge part near the gap (proof of
`SMOOTH_BOUND` §1.2 plus the edge leakage), which removes `H_LJ`; (2) the spinning `O(r)` Regge kernel (removes `H_R`); (3) `c₁`;
(4) rigorous (interval) versions of the near-gap chains (removes `H_M`); (5) the KKSDZ loophole and K1′.

### 6.3 Conditional fallback (generalising `RESULT_D_CONDITIONAL.md`)

Theorem D_cond carries over to `Ω^EE` with the band replaced by the uncovered region. Let `U(Δ_g)` be the `H_LJ` region of §5.1/§5.2,
assume `P` on the complement of `U`, and assume **`H_U(ε)`: the states in `U` contribute to `Ω^EE` at least `−ε·16L_A`**. Then
`√(t₄²+α_T²) ≤ 240√((1+ε)16L_A)/Λ⁴` [D given the hypotheses]: `×√(1+ε)` on the table (`ε = 0.1`: `3.29·10⁻⁶ / 3.03·10⁻⁸`; `ε = 1`:
`4.43·10⁻⁶ / 4.09·10⁻⁸`). Unlike the `q²` route's `H_band` (a helicity-coherence condition on the CFT spectrum, `RESULT_D_CONDITIONAL`
§3.3), `H_U` is implied by `H_LJ`, which is a property of explicit block integrals. The `q²`-route numbers of `RESULT_D_CONDITIONAL`
(band edges `b̂_c ≈ 16 / 30 / 62–68`, `κ_crit`, the re-weighted `3.38·10⁻⁶ / 3.26·10⁻⁸`) remain valid conditional statements *for
`Ω^{A′}[q²]`* but are superseded as the headline: for `Ω^EE` the leading-Regge band is absent [P + M].

### 6.4 Comparison with the literature (`RESULT_D_CONDITIONAL.md` §5, `LIT_HELICITY_MIXING.md`; functional-independent because `Ω^EE`'s flat limit and light side equal `q²`'s)

* **Flat CHLPS form [D].** The AdS graviton pole `16πG/(ν² + 9/4)` makes `m_IR = 3/(2R_AdS)` exact at log accuracy:
  **`|ĝ₃|²Λ⁸ ≤ 3022.34·log(2Λ/3) − 5327.21` (+ `O(1/Λ)`)**, the CHLPS-form bound of NF02 (`ĝ₃` normalisation identical to CHLPS: the
  project reproduces their (3.7) `1598 log − 3211` as `1596.7, −3208.0` [DV]).
* **Log coefficient** 3022.34 vs CHLPS's SDPB optimum 24.9 (2201.06602 (4.4)): ratio 121.4, i.e. ×11.0 in the bound. Decomposition
  [NE]: ×1.05 (CHLPS's extra rules), ×33.6 (no forward terms — forced in AdS by uncontrolled contact content, H9a), ×2.58 (margin
  δ = 0.2), ×1.34 (large-`b̂` proxy). Transplanting CHLPS's optimum with the same `m_IR` gives `3.09·10⁻⁷ / 2.80·10⁻⁹` (300 / 1000), 10×
  stronger — **asserted by CHLPS, not derived**; it would need forward and helicity-flip (`B⁽³⁾`) sum rules in AdS, which CLS-D's family
  does not contain.
* **Scaling:** `√(t₄²+α_T²) ≲ 240√A·√(ln Δ_g)/Δ_g⁴`, stronger by `√(ln Δ_g)` than AKT (1805.07393, `t₄ ≲ ln Δ_g/Δ_g⁴`, d = 3, unknown
  constant; AKT's normalisation not converted [I]). In AKT's form `C ≈ 4.4·10³ / 4.2·10³` at 300 / 1000. CEMZ, CHP, Meltzer–Perlmutter:
  `Δ_g⁻⁴`, parametric only.
* **Novelty [Q, literature check of 19 Nov]:** no published work derives graviton-coupling bounds from spinning CFT sum rules; CLS-D
  derive none and defer the Regge regime; the curvature-induced helicity off-diagonal (and its removal by `EE`) is not in the literature.
  CHMRS/CHLPS state their D = 4 bounds conditionally on the large-impact-parameter spectrum; the present conditions are of a different
  (computational) kind.

---

## 7. Withdrawn or superseded

| item | status | reason | where |
|---|---|---|---|
| **D-safe `Ω^{NF,SC}` (B′ uplift, `k = 2…7`) and its "`n = 1`"** | withdrawn | CFT light-ray sum rules map to the *stripped* flat rules, not B′; the B′ functional is not positive at leading order for gravitons (`−0.72 N_s` at `+−`, `1.1Λ`, `J = 8`) | box bullet 1; `H0` §9–§11 |
| **D-strong (`Ω^U` on M02)** | withdrawn | built on B′ objects, never re-derived for the stripped dictionary | review P1 |
| `RESULT_D.md` §2–§6 tables, hypotheses H0–H17 in their B′ form, "What `n = 1` means" | superseded by this note | describe the withdrawn functional | review P1, P5 |
| literal CLS-D `Ĉ₂` | withdrawn | heavy normalisation vanishes at `ν₀ = 5.979` | box (Δ_g = 100, §24) |
| scalar-`g` `Ĉ₂` (CHMRS de-shift on `q²`) | superseded | violates (HN): pole `16πG(1+ε)/(ν²+9/4)`, `ε = +0.51` at `ν = 1–3`, singular at 5.98 | `H7_I14` §19 |
| **route A′ with `Ĉ₂′` (`q²`) — "`n = 1` at `Δ_g = 300, 1000`"** | withdrawn | leading-Regge helicity off-diagonal: near-gap failures at 300 for `J ≥ 16–18` (direct-validated), modelled failure at 1000 from `J ≈ 32`; band `b̂ ≳ 0.65Λ^{2/3}` not removable within NF02's LP family [M] | `S1_S4_FIXES`; `LARGEJ_OFFDIAG`; `OFFDIAG_KERNEL` |
| `q²` "`n = 1` at `Δ_g = 100`, 0.2206" | withdrawn | used unaligned off-diagonal conventions; the corrected chain gives `n ≠ 1` (`J = 10`, 1.21Λ, `−0.20 / −0.34`) | `EE_DG100` §0 |
| `H0` §22 "`n = 1.5` at 100, `J = 10` `−2.7`" | withdrawn | artefact of the log form | box (§23) |
| `Ĉ₃` with the scalar `k = 3` de-shift | superseded | spurious; `Ĉ₃` needs no de-shift | `S1_S4_FIXES` §1 |
| `l9_neargap` off-diagonal interpolation without basis alignment | fixed | arbitrary basis signs partially cancelled the off-diagonal | `S1_S4_FIXES` |
| "`K₂` `+−` breakdown" at `J ≥ 32` (1000), `J = 32` `−3.63`; 300 `J = 28` at 1.1Λ convention A | superseded | scipy `eval_jacobi` bug for complex argument, fixed by recurrence | `EE_DG1000_LARGEJ` §1 |
| EE `J = 24, 28` failures at 300 (`−0.053 / −0.095`) | superseded | `K₂` node at `x = 0.05` has the wrong sign; spline extrapolation below it | `EE_DG300_SMALLNU` §0 |
| `J = 48` "smooth offset" (`+200 … +650 c₂`), `J = 48` [I] | superseded | spin-`J` engine precision loss `ε²4^J`; fixed engine gives `r = 0 ± 5` and `+0.932` | `FOLD_WINDOW_LARGEJ` §1, §4 |
| any unpatched graviton record at `J ≳ 48` (ee1000 `J = 48` rows, refb `J = 56, 64`, `(1064, 64)`, r1's `J = 72`); `l5_shadow_matrix` shadows at `J ≥ 56` | unreliable | engine precision / unscaled shadow fit | `FOLD_WINDOW_LARGEJ` hand-over |
| hybrid-chain λ values before 20 Nov (curvature counted twice) | superseded | `λ` too high by `≤ 0.004–0.0135 N_s` (1000), up to `0.07` (300 `J = 28`); no verdict changed | `SMOOTH_BOUND` §4 |
| canonical (CHMRS-class) gauge as the definition of `Ω^EE` | reversed | its weight (`|X_can/ψ₂|` up to ≈ 180 at `ν ≈ 5`) makes `n` at 300/100 depend on uncontrolled small-`ν` data | `EE_NEARGAP` §4 vs `EE_FOUNDATIONS` §3.3 |
| ηP / `Ĉ₂′` `δk(0)` values in `H0` §R.9; `B_AdS = −6552.5` | superseded | `δk(0)_EE` as in §3; `B_AdS = −6552.66` | review M4, P2 |
| twist convention "`Δ − J ≥ Δ_g − 2` for `J > 2`" | wrong | code and all scans use `Δ − J ≥ Δ_g` for all `J` | review M1 |
| "for all `Δ_g ≥ 300`" | withdrawn | only 300 and 1000 computed | review M10 |
| "certified window `b̂ ∈ [91.2, 2Λ]`" and Theorem A as graviton statements | re-graded | both concern the scalar-skeleton model | review S2; `CERTIFIED_NEARGAP_A` |
| odd-sector `x = 0.05` failure at `Δ_g = 100` (`J = 9`, `−1.345`) | likely artefact | computed with the spurious `Ĉ₃` de-shift | box (NEARGAP100_SMALLNU) |
| odd `J = 29` at 1000 (`−1.20` in the raw chain) | artefact | an isolated `K₂` outlier in the `Δ_ref = 2·10⁴` C3H reference at `ν = 822` (`|M₂/M₀| = 0.86` vs `≤ 0.03` at neighbours); odd rows are not despiked by the chain | §9.1 |

---

## 8. Contradictions and inconsistencies between documents (flagged)

1. **Gauge.** `EE_FOUNDATIONS` §3.3 *adopts* the canonical gauge as the definition; `EE_NEARGAP` §4, `EE_DG300_SMALLNU` §6,
   `EE_DG100` and the `RESULT_D` box adopt the **min-norm gauge** as the stated convention. Later wins: min-norm. This is a real
   dependence, not just bookkeeping: with the canonical weight the same data give `n = 1.55` (spline) or `−22 … −40` (direct, CV)
   at 300 for `J = 24`. Each gauge is a valid functional, so choosing min-norm is legitimate, but the statement is for the min-norm
   functional and should say so.
2. **`RESULT_D_CONDITIONAL` §2:** "r5's convention `α ≥ Δ_g − 2`, which is slightly more permissive than A1". With `α = Δ − J − 2`,
   `α ≥ Δ_g − 2` *is* A1; the remark is probably about continuous `m` sampling. Minor; no number depends on it.
3. **Odd rows are never despiked.** `EE_NEARGAP` §5 describes a Hampel filter in "the main chain"; in `w2_chain._grid` it acts only on
   even rows (keys `(kind, h)`), not on the odd rows (`P4.grid`, keys `'C2'`, `'C3H'`). This produced the spurious odd `J = 29`
   failure (§9.1). For the previously reported odd `J ≤ 15` (300, 1000) see §9.2.
4. **Superseded intermediate verdicts still standing in older docs:** `EE_NEARGAP` (`J ≥ 40` at 1000 [I] from the `K₂` breakdown; `n` at
   300 [I] from `J = 24, 28`), `EE_DG1000_LARGEJ` (`J = 48` [I]), `SUBLEADING_REGGE` (`J = 48` offset "unexplained") — all superseded
   as listed in §7.
5. **`H0` §R** still lists the pre-`Ĉ₂′`/pre-`EE` items 3, 5, 9, 10 (review P3); §2 here is the current definition.
6. **`Δ_g = 100` tension (q² era):** `RESULT_D_CONDITIONAL` §2.3 found the band model negative at near-gap `J = 10`, 1.1Λ while the chain
   was positive; for `Ω^EE` the leading-Regge off-diagonal vanishes, so this tension does not arise. Not a live contradiction.
7. **`J`-grid gaps:** "even `J ≤ 40`" (1000) and "even `J ≤ 22`" (300) in the summaries mean the computed sets of §5 (`J = 26, 30, 34, 38`
   at 1000 and `26, 30` at 300 have no rows). Minor, but the phrase "all `J ≤ …`" is an interpolation [M].

---

## 9. Computations made for this note (one core; no existing file edited)

### 9.1 Odd `J = 19–31` at `Δ_g = 1000` (Task 1)

* Rows: `cache/ee1000_p3_{EE,C3H}_J{19..31}.pkl` (35 nodes `x = 0.05–0.8`, state `Δ = 1000 + J` and reference `2·10⁴`; `J²/m =
  0.35–0.93`, inside the `K₂` range `≲ 1.5`), finished by `ee1000_queueR1.sh`.
* Chain: `ee1000_chain` via `sm_chainfix.py ee1000` with the curvature double-count fix (`SM_CURV=row`, κ image), shape model
  `K = 4`, `x ∈ [0.05, 0.8]`, `W2_DMIN=1000`, hybrid states `m/Λ ∈ {m_gap, 1.05, …, 5}`. Outputs
  `cache/sm_chainfix_ee1000_fodd_*.txt` (raw) and `cache/sm_chainfix_ee1000_fodd_hampel_*.txt` (with the odd-row despike).
* **Raw chain:** every odd `J` positive except **`J = 29`: `−1.20` at the gap mass**, in every representation variant (`K = 4/5/6`,
  `x ≤ 0.9`, spline (`−3.65`), with LO nodes, 64 `p`-nodes, Bessel image). Cause: the `J = 29` C3H node at `x = 0.8` has
  `Γ = −17933` (neighbours `±2000`; `Im Γ = −99` against `≈ 0` elsewhere), and the shape fit holds `c, s, r` beyond the last node up to
  `x_max = Λ/m ≈ 0.97`. The raw records show the source: the `Δ_ref = 2·10⁴` reference at `ν = 822` has `|M₂/M₀| = 0.86` (C3H) and
  `0.077` (EE), against `≤ 0.03` / `≤ 0.002` at the neighbouring nodes — an isolated `K₂` outlier of the kind `EE_NEARGAP` §5 reports.
* **Fix:** w2's Hampel rule (window ±3, `6·max(MAD, median|Γ|)`, first node kept) applied to the odd rows in-process
  (`scratchpad/odd_despike.py`). It drops `x = 0.8` at `J = 23, 29` and `x = 0.775, 0.8` at `J = 19`, nothing else.

| odd `J` | min `λ/N_s` (main) [`m/Λ`] | uncorrected | max `|dF|/N_s` | range over variants (`K5, K6, spline, x ≤ 0.775, x ≤ 0.7, with LO, 64 p, Bessel, no curv fix`) | `n` |
|---|---|---|---|---|---|
| 17 | **+0.839** [1.05] | +0.844 | 0.008 | +0.838 … +0.844 | 1 |
| 19 | **+0.916** [1.05] | +0.925 | 0.012 | +0.912 … +0.918 | 1 |
| 21 | **+0.942** [1.10] | +0.952 | 0.014 | +0.939 … +0.943 | 1 |
| 23 | **+0.887** [1.023, gap] | +0.902 | 0.015 | +0.866 … +0.889 | 1 |
| 25 | **+0.918** [1.025, gap] | +0.931 | 0.013 | +0.883 … +0.926 | 1 |
| 27 | **+0.954** [1.05] | +0.963 | 0.016 | +0.923 … +0.963 | 1 |
| 29 | **+0.979** [1.75] | +0.978 | 0.004 (0.12 spline) | +0.882 … +0.984 | 1 |
| 31 | **+0.995** [3.0] | +0.983 | 0.21 (positive, gap mass) | +0.947 … +0.996 | 1 |

* **Verdict [NE + M]: `n = 1` for every odd `J ≤ 31` at `Δ_g = 1000`, all masses from the gap to `5Λ`;** the minima are set by the base
  model's oscillation in `m` (uncorrected values `0.84–0.98`), the graviton corrections are `≤ 0.016 N_s` (`≤ 0.12` in the worst variant)
  except the positive `+0.21` at `J = 31` (gap mass), which comes from the large `x = 0.775–0.8` nodes (`±2000–3600`) held beyond the row
  (`x ≤ 0.7`: `+0.004`). Sign-safe either way.
* **Precision (engine fix) [NE]:** a phase rotation of the spin-`J` null-circle sampling (an exact identity) at `J = 31`, `x = 0.5`, with
  the *unpatched* engine: state `(1031, 31)` EE, phases `0 / 0.37 / 1.3`: relative change of `M₀` `≤ 7.5·10⁻¹⁴`, `M₂` `≤ 1.8·10⁻¹¹`,
  `M₄` `≤ 2·10⁻¹²`; recomputation vs the cached record `1.7·10⁻¹¹` (C3H alike; reference `(2·10⁴, 31)` `≤ 4·10⁻¹²`). In `Γ` units
  (`m² ≈ 10⁶`) that is `≲ 10⁻⁴ c₂`, consistent with `ε²4^J ≈ 2·10⁻¹³`. **The pre-fix `J ≤ 31` rows are unaffected.**

### 9.2 Odd rows `J ≤ 15` with the despike (check of item §8.3)

`EE_NEARGAP`'s odd-sector chain (`w2`, rows 200–1600, `J = 5–15`, `SM_CURV=row`) re-run with and without the odd-row despike
(`cache/sm_chainfix_w2_oddcheck_{nodesp,desp}.txt`). The filter drops only `x = 0.8` on the `J = 9` rows 800/1600; the only change is
`J = 9` at 1000, `+0.4541 → +0.4546`. All odd `J ≤ 15` minima are otherwise identical: 1000 `+0.2021 / 0.311 / 0.455 / 0.410 / 0.605 /
0.869`, 300 `+0.2052 / 0.704 / 0.438 / 0.749 / 0.728 / 0.707` (`J = 5 … 15`). **The earlier odd-sector results stand** [NE + M].

### 9.3 `Δ_g = 300`, `J = 32` with direct small-`ν` nodes

`EE_DG300_SMALLNU` §7 left the `J = 32` direct n3 nodes (`ν = 8.3 … 33.3`) computing; they finished (`cache/z2*_pts_*_J32.pkl`) but
were never analysed. Run here: `Z4_EST=best Z4_JS=32 Q3_JEVEN=32 W2_ROWS=336,1040 SM_CURV=row SM_JACFIX=1 sm_chainfix.py z4 300`
(`cache/sm_chainfix_z4_row_jac1_kappa_j32*.txt`):
* `Γ` held below `ν = 8.3`: **`J = 32` min `+0.657`** at 1.15Λ (uncorrected `+0.967`; max `|dF_hh|/N_h` 0.33; off-diagonal 0.013);
  `Γ_c = 0` below 8.3: `+0.738`; `Γ_c = −2000` below 8.3: `+0.388` ⇒ failure at a uniform `Γ_c ≈ −4200 c₂` below `ν = 8.3`
  (cf. `−4000` at `J = 28`). No thimble nodes exist for `J = 32`. Convention A (`K₂` + spline) had `−0.062` — the same small-`ν`
  artefact as `J = 24, 28`.
* Caveat: the lowest chain state is 1.15Λ; the twist-allowed edge is `1.102Λ` (`1.1Λ` is excluded), so `m ∈ [1.102, 1.15)Λ` is not
  sampled [H_cont].
* Grade: [NE + M] on central values; [I] in the `ν < 8.3` window (as `J = 28`).

---

## Hand-over

* **Files.** This note only; scratch scripts in the session scratchpad (`odd_run.sh`, `odd_despike.py` — in-process odd-row Hampel
  patch wrapping `sm_chainfix.py`, takes the mode as arguments —, `odd_rows.py` (row dump), `phase_check.py` (phase-rotation check;
  writes only to the scratchpad), `odd_flags.py`, `odd_w2.sh`). New cache outputs (new names only):
  `cache/sm_chainfix_ee1000_fodd_*`, `cache/sm_chainfix_ee1000_fodd_hampel_*`, `cache/sm_chainfix_z4_row_jac1_kappa_j32*`,
  `cache/sm_chainfix_w2_oddcheck_*`.
* **Reproduce §9.1:** `cd ads4/code/spinfunc; W2_DMIN=1000 EE1000_JODD=17,19,21,23,25,27,29,31 EE1000_JEVEN= EE1000_MODEL=shape
  EE1000_K=4 EE1000_XMAX=0.8 EE1000_XMIN=0.05 SM_CURV=row SM_TAG=fodd_row_kappa python3 sm_chainfix.py ee1000` (raw); the same with
  `python3 <scratch>/odd_despike.py` in place of `python3 sm_chainfix.py` (despiked, mode argument `ee1000` is the default). ≈ 11 s each.
* **Recommended code change (not made):** extend `w2_chain._grid`'s Hampel filter to `q3_analysis.P4.grid` (odd rows), or recompute the
  `J = 29` reference node at `ν = 822` (C3H, EE).
* **For whoever edits `RESULT_D.md`:** replace §1–§6 by §1–§7 of this note; collapse the update box to a pointer here (review P5).
* **Highest-value next steps:** (1) analytic control of the smooth beyond-Regge part near the gap (removes `H_LJ` at 1000 and most of it
  at 300); (2) odd `J ≥ 17` rows at `Δ_g = 300` (the odd correction was still growing at `J = 15`); (3) thimble nodes `ν < 8.3` for
  `J = 28, 32` at 300 with ×10–100 statistics; (4) `c₁`; (5) the spinning `O(r)` Regge kernel (`H_R`); (6) the even-`J` OO block.
