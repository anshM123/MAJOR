# A flat graviton functional for `|ĝ₃|²` with a robust positivity margin, and its AdS4 threshold

Status: written 16 Nov 2026 by a derivation/numerics agent. **Unrefereed.** Code: `ads4/code/robust/` (`bash run_all.sh`; outputs in
`robust/cache/`, git-ignored). It copies the certificate machinery of `ads4/code/flatbound/` (Lemmas P, C, M, S of
`FLAT_GRAVITON_BOUND.md`) and the AdS machinery of `ads4/code/adsfunc/` (`ADS4_FUNCTIONAL.md` t2/t3), generalised to arbitrary
polynomial weights. It imports only `flatbound/common.py` (Clebsch–Gordan bound `cg_abs_bound`, `to_fmpz_poly`), unchanged.

Marks: **[certified]** FLINT/Arb ball arithmetic or exact rationals, given the lemmas named; **[D]** derived (elementary, or from the
named lemmas); **[NE]** floating point (LP solutions, grids, scans); **[M]** statement about the hybrid *model* of
`ADS4_FUNCTIONAL.md` §2.4; **[I]** inferred/assumed. LP outputs are [NE] until a rationalised version is certified.

## 0. Results

Problem. `F_{g3}` (completed, `g₅`-free) is positive on all heavy states [certified, `FLAT_GRAVITON_BOUND.md`], but only marginally
(`m¹⁰F = 19/315` at `+−, J = 4, m = M`; relative margin `≈ 10⁻³` near `m ≈ 21.8M, J ≈ 245`). In the AdS uplift this makes the
leading Regge action negative up to `m ≈ 7.5Λ`, and with an `O(Λ²/m²)` envelope the positivity threshold is `≈ 8.1Λ`, `≈ 22Λ` at
`Δ_g = 1000` (`ADS4_FUNCTIONAL.md` §3(a)), so the conditional AdS bound (`ADS4_LIGHT.md`, `37.84 log Λ − 60.79`) needs heavy
positivity down to `m ≈ Λ`, which is not available (referee M2: the bound weakens like `n⁸` in `|ĝ₃|²` if only `m ≥ nΛ` is usable).

1. **Family and criterion [D].** Weights `ψ₂ = (1−p)³p Σᵢ₌₀⁴ xᵢpⁱ`, `ψ₃ = (1−p)³ Σⱼ₌₀⁴ yⱼpʲ`, forward term
   `−a∂_{p²}B₄^{(1)}|₀ + aB₅^{(1)}|₀`, `a ≥ 0` (`g₄`, `g₅` eliminated). The heavy action per `|c|²` (normal form `𝔉 = m⁴F`, `u = M²/m²`)
   is `𝔉 = B2 + B3 + FW`, with `B2 = ∫ψ₂(2−up²)D_h`, `B3 = σ_hu∫ψ₃D_h`, `FW = a(2𝒥 − c₄ + σ_h)u³`, `c₄ = 1 (++), 41 (+−)`. Robustness
   reference: the component norm `N = |B2| + |B3| + |FW|` (the `N_comp` used in the AdS envelopes of `ADS4_FUNCTIONAL.md`).
   **Uniform margin δ:** `𝔉 ≥ δN` on every heavy state `m ≥ M`, both channels, all allowed `J`.
2. **Pareto curve [NE] (§2).** Maximal margin vs bound constants (LP over the family, `k = 1`, minimise `A`):
   `δ = 0: A = 26.1`; `0.1: 35.5`; `0.2: 48.5`; `0.3: 67.0`; `0.5: 144`; `0.6: 236` (degree 4; degree 6 lowers `A` by 2–3 %).
   `B ≈ −1.1A` throughout. For comparison `F_{g3}`: `A = 37.84`, `B = −45.45`, `δ = 0` (dense grid: min `𝔉/N = 8.1·10⁻³` in the `++` island, `1.2·10⁻³` at the `+−` near-zero).
   Imposing in addition a linear proxy that makes the large-`b̂` (Lemma M) certificate work at `b̂ ≥ 120` costs little:
   `δ = 0.1: A = 37.55, B = −46.93`; `δ = 0.2: A = 51.74, B = −65.32` (§2.2).
3. **Chosen functional `M02` and its certificate [certified] (§3–4).** Exact rational weights (§3, scaled so that `k ≈ 27.49` like
   `F_{g3}`): `ψ₂ = (1−p)³p(711/8 + 129.649p + 561.721p² − 1013.557p³ + 997.013p⁴)`,
   `ψ₃ = (1−p)³(54.971 − 1025.83p − 1276.12p² + 946.687p³ + 15068.662p⁴)`, `a = 1129/250`.
   **Theorem R.** For every heavy state `m ≥ M`, every allowed `J`, both helicity channels: `𝔉 ≥ (9/50)·(|B2| + |B3| + |FW|)`.
   [certified: regions P (`J < 400`, all `m ≥ M`; 2384 exact polynomials, Descartes), C (`J ≥ 400`, `m ≤ 2M`), M (`m ≥ 2M`,
   `b̂ ≥ 120`, Lemma M with exact `1/b̂`, `1/b̂³` terms + remainder), S (`J ≥ 400`, `m ≥ 2M`, `b̂ ≤ 120`, Cauchy/Taylor cells) — the
   regions and lemmas of `FLAT_GRAVITON_BOUND.md` §3, applied to the four sign patterns of `(B2, B3)`.] In particular `𝔉 > 0`
   (Theorem-1 analogue). Floating point: min `𝔉/N = 0.1996` (`++`, `m ≈ 2.7M`, `J = 8`).
   **Low side [certified, exact]:** `g₄, g₅, g₆, g₆′` drop out, `k = 10884217/396000 = 27.48540`, `A = 16x₀/k = 51.7366`,
   `B = −65.3207`, `B_AdS = B − A log(3/2) = −86.2981`; both light-matter coefficients negative on `0 < m_ℓ ≤ M` (Sturm).
   Hence (flat, hypotheses H1–H5 of `FLAT_GRAVITON_BOUND.md` §4.1, same IR caveats as its Theorem 3(i)):
   `|ĝ₃|²M⁸ ≤ 51.74 log(M/m_IR) − 65.32 + Δ_IR`.
   A second certified functional `M01` (`𝔉 ≥ (9/100)N`, all regions PASS) has `A = 37.5457`, `B = −46.9286` — **at least as good as
   `F_{g3}` in both constants** (`37.8387`, `−45.4452`) while having a uniform 9 % margin instead of `≈ 10⁻³`.
4. **Envelope thresholds.** *Flat, from Theorem R [D]:* `𝔉 − C₁uN ≥ (9/50 − C₁u)N > 0` for `m > (50C₁/9)^{1/2}M`:
   `M₀ ≤ 2.36M (C₁ = 1)`, `4.08M (3)`, `7.45M (10)`. *Leading AdS Regge action with the envelope `C₀cosh β/m² + C₁Λ²/m²` [NE] (§5):*
   `M02` is positive on every physical grid point **without** envelope (min `F_R/N ≈ 0.19`, all `Δ_g = 20…1000`), and with envelope
   the last negative point is `m ≈ 2.16Λ (C₁ = 1)`, `3.92Λ (3)`, `6.70Λ (10)`, **independent of `Δ_g`** (`F_{g3}` on the same grids:
   `8.09Λ / 22.1Λ` at `Δ_g = 100 / 1000` for `C₁ = 1`, `10.1Λ / 24.6Λ` for `C₁ = 10`). Hybrid model [M]: `1.92Λ`, `3.58Λ`, `6.14Λ`.
5. **Conditional AdS bound [NE threshold, D otherwise] (§5.4).** With the smearing scale lowered to `Λ_s = Λ/n` so that every heavy
   state lies above the threshold (`n = 2.2` for `C₁ = 1`; the rescaled Regge scan is then positive on every physical point at
   `Δ_g = 100, 1000`), the light side of `ADS4_LIGHT.md` gives `|ĝ₃|²Λ⁸ ≤ n⁸·16L_AdS(Λ/n)/k`, i.e.
   `|ĝ₃|²Λ⁸ ≤ 2.84·10⁴ log Λ − 6.97·10⁴ + …` (`= 6.3·10⁴` at `Δ_g = 100`), against `7.1·10⁸` (`Δ_g ≤ 300`) and `8.9·10¹¹`
   (`Δ_g = 1000`, threshold `22Λ`) for `F_{g3}` with its own `C₁ = 1` threshold, and `113.5` for `F_{g3}` *if* it were positive down to
   `m = Λ`. Functionals optimised directly for the envelope (§2.3, [NE] only, not certified) do better: `≈ 525` at `Δ_g = 100`
   (`C₁ = 1`, `n = 1.75`; `|ĝ₃|²Λ⁸ ≤ 227 log Λ − 531`), `≈ 1.8·10³` (`C₁ = 3`, `n = 2.25`), `≈ 7·10⁴` (`C₁ = 10`, `n = 4`); their rescaled
   Regge actions are positive with the respective envelope on all physical grid points at `Δ_g = 100, 1000`, but their uniform margin
   far above threshold is only `≈ 0.02`.
6. **What is not achieved.** `M₀ ≈ Λ` with `C₁ ≥ 1` is impossible in this criterion: at `m = M` (`u = 1`) the envelope demands
   `𝔉 ≥ C₁N ≥ C₁𝔉`. For `C₁ = 1` the best found is `n ≈ 1.75` (targeted family, §2.3) and `n = 2.2` for the certified M02. Everything AdS-side remains conditional on the inputs listed in
   `ADS4_FUNCTIONAL.md` §4 and `ADS4_LIGHT.md` §4 (heavy `⟨TTTT⟩` formula at `k = 2…5`, spinning bulk-point formula, the envelope
   itself as a model of the `O(Λ²/m²)` corrections, (F), (R), (P)). Nothing here is a theorem about AdS4.

## 1. Setup

### 1.1 Functional, low side [D; certified per functional]
`F = ∫₀¹dp[ψ₂B₂^{(1)} + ψ₃B₃^{(1)}] − a∂_{p²}B₄^{(1)}|₀ + aB₅^{(1)}|₀`. With the low sides of `FLAT_GRAVITON_BOUND.md` §1.3 ((2.31), (2.32),
(3.8) corrected): `F|_high = 16πG L(ε) − πG k|ĝ₃|² + C₀(m_ℓ)a₀ + C₂(m_ℓ)a₂`, `k = −2μ₂₆ + 2μ₃₄ + 4a`, `L(ε) = x₀ log(1/ε) + L_fin + O(ε)`,
`L_fin = −(11/6)x₀ + Σ_{i≥1} xᵢ·6(i−1)!/(i+3)!`, `C₀ = μ₂₂ − 2m_ℓ²μ₂₀ − μ₃₀`,
`m_ℓ⁴C₂ = ∫ψ₂(p² − 2m_ℓ²)(m_ℓ⁴ − 6m_ℓ²p² + 6p⁴) − ∫ψ₃(m_ℓ⁴ − 6m_ℓ²p² + 6p⁴) − 12a` (`μ_{ij} = ∫ψᵢpʲ`). `g₄` never appears
(no `B₄` at `p = 0`), `g₅` cancels between the two forward terms, `g₆, g₆′` do not enter. So `A = 16x₀/k`, `B = 16L_fin/k`
(`ε → 0`), and in AdS (`ADS4_LIGHT.md` §2.6) `L_AdS(Λ) = ∫ψ₂/(p² + 9/(4Λ²)) = x₀ log(2Λ/3) + L_fin + O(1/Λ)`, `B_AdS = B − A log(3/2)`.
`F_{g3}` is `x = (65, 155, 47, 0, 0)`, `y = (162, −3854, 9202, 0, 0)`, `a = 5` (check: `k = 330151/12012`, `A = 37.8387`, `B = −45.4452`,
`B_AdS = −60.7875`, reproduced by `rcommon.py`). For each chosen functional `lowcheck.py` recomputes everything from the amplitude
(2.7a) and (B.1) by residues (independently of `cert_low.py`) and checks: `g₄, g₅, g₆, g₆′` coefficients `0`, `|ĝ₃|²` coefficient
`−πGk`, pole `16πG/p²`, `C₀, C₂ < 0` on `(0, 1]` (exact root counts).

### 1.2 Heavy action [D]
Per state (common factor `16(2J+1)dm²/m²` dropped), generalising `FLAT_GRAVITON_BOUND.md` (0.1) with `∂_{p²}C₄|₀ = (c₄ − 2𝒥)/m¹⁰`,
`C₅|₀ = σ_h/m¹⁰`: `𝔉_h(J,u) = ∫₀¹[ψ₂(2 − up²) + σ_huψ₃]D_h(1 − 2up²)dp + a(2𝒥 − c₄ + σ_h)u³` (`a = 5` gives `Φ₊₊ = 10𝒥`,
`Φ₊₋ = 10𝒥 − 210`). Lemma P (exact polynomial in `u`), the continuous extension `G_h(b̂,u)`, Lemma C and Lemma M hold verbatim for
any weights of the form `W₀(p) + uW₁(p)` with `(1−p)³ | W₀, W₁` and `W₀(0) = 0`.

### 1.3 Robustness criteria (LP constraints)
* (R1) physical states `J` allowed, `u ∈ (0, u₀]`: `𝔉 − (δ + C u)N ≥ 0`, `N = |B2| + |B3| + |FW|`. Uniform margin: `C = 0, u₀ = 1`.
  Envelope-targeted (§2.3): `δ = 0.02`, `C = C₁`, `u₀ = 1/n²`.
* (R2) scaling limit `u → 0` at fixed `b̂`: `2∫ψ₂J₀(b̂p) ≥ δ ρ(b̂)`, `ρ = ρ₀(1 + (b̂/b₀)²)^{−3/2}` with the `F_{g3}` scales
  (`ρ₀ = 0.4489`, `b₀ = 1.865`, tail `2.91/b̂³` per unit `k`): the flat "impact-parameter scale" reference, which (R1) alone does not
  constrain at `u = 0` (there `𝔉/N → 1` wherever `f̂ > 0`).
* (R3) large `b̂ ∈ [150, 10⁴]`, `u ∈ [10⁻⁹, ¼]`: (R1) applied to the two-term Lemma-M asymptotics
  `B2 ≈ −ψ₂″(0)/b̂³`, `B3 ≈ σ_hu[ψ₃(0)/b̂ − ψ₃″(0)/(2b̂³)]`, `FW = a(b̂²u²/2 − u³/2 − (c₄ − σ_h)u³)`. Thus the `1/b̂³` coefficient
  `−ψ₂″(0)` is positive with margin, and the `±ψ₃(0)u/b̂` term is beaten by the forward term (`F_{g3}`: `−ψ₂″(0) = 80`, `ψ₃(0) = 162`;
  `M02`: `274`, `55`).
* (M) matter: `C₀, m_ℓ⁴C₂ ≤ 0` on a 201-point grid in `m_ℓ ∈ [0, 1]` (exact check afterwards).
* Normalisation `k = 1`. (R1) is the maximum of four linear constraints (sign patterns of `B2`, `B3`; `sign FW` is fixed by `J`) and is
  imposed by cutting planes. Grid: 117 values `u ∈ [2.5·10⁻⁵, 1]`, all `J ≤ min(2500, (200/√u − 1)/2)`, ≈ 1.8·10⁵ states.

## 2. LP results [NE]

### 2.1 Uniform margin, no further constraints (`pareto.py`)
Minimise `A` (tie-break `B`) at fixed `δ` (`dsc = δ`):

| `δ` | `A` (deg 4) | `B` | `B_AdS` | `A` (deg 6) | `B` (deg 6) |
|---|---|---|---|---|---|
| 0 | 26.14 | −28.68 | −39.28 | 25.59 | −28.14 |
| 0.05 | 30.47 | −33.77 | −46.12 | 29.77 | −32.94 |
| 0.10 | 35.53 | −39.61 | −54.02 | 34.72 | −38.79 |
| 0.15 | 41.44 | −46.11 | −62.91 | 40.51 | −45.30 |
| 0.20 | 48.46 | −53.90 | −73.55 | 47.34 | −52.73 |
| 0.25 | 56.83 | −63.04 | −86.09 | 55.51 | −61.70 |
| 0.30 | 66.99 | −73.99 | −101.15 | 65.47 | −72.51 |
| 0.40 | 95.58 | −105.08 | −143.84 | 93.43 | −102.65 |
| 0.50 | 143.68 | −157.62 | −215.88 | 140.74 | −152.87 |
| 0.60 | 236.08 | −263.40 | −359.12 | 233.11 | −258.42 |

(`F_{g3}`: `37.84, −45.45`, `δ = 0`.) The `δ = 0` entry, `A = 26.1`, is above but close to the paper's SDPB optimum `24.9 log − 27.6`
(arXiv:2201.06602 (4.4), which uses improved sum rules and more structures); it is a grid statement only. The cost of a margin is
smooth: `A(δ) ≈ 26 e^{3.1δ}` up to `δ ≈ 0.4`. These unconstrained optima have small `−ψ₂″(0)` (≈ 33 in the `k ≈ 27.5` scaling), which
makes the Lemma-M certificate at `b̂ ≥ 120` fail (it needed `b̂ ≳ 600`).

### 2.2 Certifiable family (`mtest.py`, `pareto2.py`)
Add the linear sufficient proxy for the Lemma-M certificate at `B_A = 120` (`lp.py: m_proxy`):
`(1−δ)(−ψ₂″(0))B_A ≥ 1.4·(4/3π)[2(1+δ)(|ψ₂‴(1)| + ‖ψ₂⁗‖₁) + u₁(1+δ)(…ψ₃, p²ψ₂…)] + 1.2(1+δ)²Y²/(2(1−δ)a_min)`, `|ψ₃(0)| ≤ Y = 2`,
`a ≥ a_min = 0.14` (per unit `k`). Results (`B_A = 120, Y = 2`): `δ = 0.1: A = 37.54, B = −46.93`; `δ = 0.2: A = 51.74, B = −65.32`
(`B_A = 200`: `36.43/−42.79`, `49.86/−58.50`). Pareto curve of this family (degree 4, `pareto2.py`):

| `δ` | 0 | 0.05 | 0.10 | 0.15 | 0.20 | 0.25 | 0.30 | 0.35 | 0.40 | 0.50 |
|---|---|---|---|---|---|---|---|---|---|---|
| `A` | 27.47 | 32.11 | 37.54 | 43.99 | 51.74 | 61.10 | 72.55 | 86.66 | 104.4 | 158.1 |
| `B` | −34.11 | −39.93 | −46.93 | −54.93 | −65.32 | −76.41 | −90.83 | −109.4 | −133.0 | −206.6 |
| `B_AdS` | −45.25 | −52.94 | −62.15 | −72.77 | −86.30 | −101.2 | −120.2 | −144.5 | −175.3 | −270.7 |

The proxy costs `≈ 5–10 %` in `A` and gains in `B`; at `δ = 0.1` the result is at least as good as `F_{g3}` in both constants.
Degree 6 (`pareto2_d6`): `A` lower by `≈ 1 %` (`δ = 0.1: 37.29/−45.41`; `δ = 0.2: 50.93/−62.20`), so degree 4 is adequate.

### 2.3 Envelope-targeted family (`target.py`)
For `C₁ ∈ {1, 3, 10}` and `n`, impose (R1) with `δ = 0.02`, `C = C₁`, `u₀ = 1/n²` (smearing scale `Λ_s = Λ/n`, so physical states
have `u ≤ 1/n²`; the grid is augmented by 38 values `u ∈ [0.7u₀, u₀]`, since the constraint binds at the endpoint — a first run
without them was violated just above `m = n` and is discarded), (R2) with `0.02`, (R3), (M), and minimise the exact light-side value
`16L_AdS(Λ/n)` at `Λ = 100`. The would-be bound is `|ĝ₃|²Λ⁸ ≤ n⁸·16L_AdS(Λ/n)` (`k = 1`):

| `C₁` | `n` | `A` | `B_AdS` | `16L_AdS(100/n)` | `n⁸ ×` |
|---|---|---|---|---|---|
| 1 | 1.3 | 35.96 | −57.03 | 100.0 | 816 |
| 1 | 1.5 | 9.05 | −16.29 | 22.29 | 571 |
| 1 | 1.75 | 2.58 | −4.59 | 6.01 | **528** |
| 1 | 2.0 | 1.03 | −1.89 | 2.21 | 566 |
| 3 | 2.0 | 5.52 | −9.73 | 12.29 | 3145 |
| 3 | 2.25 | 1.33 | −2.42 | 2.73 | **1792** |
| 3 | 2.5 | 0.735 | −1.34 | 1.43 | 2175 |
| 10 | 4 | 0.849 | −1.86 | 1.12 | **7.3·10⁴** |
| 10 | 5 | 0.962 | −2.37 | 0.97 | 3.8·10⁵ |

`C₁ = 3` is infeasible for `n ≤ 1.5` and `C₁ = 10` for `n ≤ 3` (first run). For large `n` the objective saturates (`16L_AdS ≈ 1`):
the functional then only has to be positive far above its own smearing scale. These functionals (`T1`, `T3`, `T10` below) are dominated
by `ψ₃` and the forward term, are negative near `m = Λ_s` (e.g. `𝔉 < 0` at `+−, J = 4, u = 1`; irrelevant since `m ≥ nΛ_s`), and are
**not certified**.

## 3. The chosen functionals
`mkfun.py` rationalises the LP solutions (coefficients scaled by `k(F_{g3}) = 27.485…` and rounded to `10⁻³`; `y₀` then raised by
`≤ 2·10⁻³` so that the scalar-matter coefficient `C₀(0) = μ₂₂ − μ₃₀` is strictly negative, since the LP puts it at `0`):

* **M02** (δ_LP = 0.2): `x = (711/8, 129649/1000, 561721/1000, −1013557/1000, 997013/1000)`,
  `y = (54971/1000, −102583/100, −31903/25, 946687/1000, 7534331/500)`, `a = 1129/250`.
  `k = 10884217/396000`, `A = 51.736565`, `B = −65.320678`, `B_AdS = −86.298050`; `C₀ = −(248792918m_ℓ² + 2863)/13860000`,
  `C₂ = −(3234307934m_ℓ⁶ − 3807692525m_ℓ⁴ − 7642880220m_ℓ² + 14856956205)/(180180000m_ℓ⁴)`.
* **M01** (δ_LP = 0.1): `x = (64493/1000, 3779/40, 208699/500, −790783/1000, 884501/1000)`,
  `y = (1001/50, −74971/250, −3284623/1000, −2054369/1000, 2261339/100)`, `a = 4259/1000`. `k = 4951986353/180180000`,
  `A = 37.545657`, `B = −46.928647`, `B_AdS = −62.152101`.

Floating-point dense check (`check.py`; 160 `u` values in `[10⁻⁶, 1]`, `J ≤ min(3000, (200/√u − 1)/2)`) [NE]:

| functional | min `𝔉/N` `++` (at `m/M, J`) | min `𝔉/N` `+−` | min `2f̂/ρ` (per unit `k`) | flat thresholds `n(C)` = sup{`m/M`: `𝔉 < CuN`}, `C = 1 / 3 / 10` (`++`; `+−`) |
|---|---|---|---|---|
| F_g3 | `8.1·10⁻³` (6.9, 30) | `1.2·10⁻³` (21.78, 245) | 0.55 | 7.30 / 8.06 / 9.36 ; 21.78 / 22.89 / 24.06 |
| M01 | `0.09998` (2.33, 6) | `0.09996` (1, 6) | 0.0996 | 2.84 / 4.44 / 6.95 ; 1.22 / 1.91 / 4.91 |
| M02 | `0.1996` (2.70, 8) | `0.19998` (1, 4) | 0.1999 | 1.91 / 3.64 / 5.98 ; 1.22 / 1.91 / 4.02 |

The margins are attained at many points simultaneously (the LP optimum is "flat"), mostly at small `J` and `m ≲ 3M` and in the
scaling region `b̂ ≈ 9`; there is no isolated near-zero. The flat thresholds agree with the AdS Regge-scan thresholds of §5.1 to ≈ 10 %.

## 4. Certification (`rcert.py`, `lowcheck.py`)
**What is certified.** For a functional `(ψ₂, ψ₃, a)` with rational coefficients and a rational `δ′`, the inequality
`𝔉 − δ′N ≥ 0` is equivalent to `𝔉[g_{s₂s₃}] ≥ 0` for the four generalised actions
`W₀ = 2(1−δ′s₂)ψ₂`, `W₁ = −(1−δ′s₂)p²ψ₂ + (1−δ′s₃)σ_hψ₃`, forward `(1−δ′s_F)a(2𝒥 − c₄ + σ_h)u³` (`s_F = sign(2𝒥 − c₄ + σ_h)`,
`= −1` only for `+−, J = 4`), `s₂, s₃ = ±1`, since `min_{s}(s₂B2 + s₃B3) = −|B2| − |B3|`. Each is certified on the four regions of
`FLAT_GRAVITON_BOUND.md` §3 with the generalised weights:

| region | method (unchanged lemmas) | M02, `δ′ = 9/50` | M01, `δ′ = 9/100` |
|---|---|---|---|
| P: `J < 400`, all `u ∈ [0,1]` | Lemma P polynomial (degree ≤ `J+1`), exact Descartes + bisection | PASS (2384 polynomials) | PASS (2384) |
| C: `J ≥ 400`, `u ≥ ¼` | `|D| ≤ 1` resp. `min(binom(J+4,8), (1−p²)⁻⁴)`; interval integrals | PASS | PASS |
| M: `u ≤ ¼`, `b̂ ≥ 120` | Lemma M (exact `1/b̂`, `1/b̂³` terms, remainder `C_R/b̂⁴` by Taylor jets in Arb), Lemma C for `+−`; quadratic minimisation in `u` | PASS (`b̂⁴𝔉 ≥ 1.2·10⁴` vs `C_R ≤ 1.2·10⁴…2.5·10⁴`) | PASS |
| S: `J ≥ 400`, `u ≤ ¼`, `b̂ ≤ 120` | exact Taylor coefficients of `G_h(b̂,u)` + Cauchy bound from Lemma S; 604 cells × 4 patterns per channel; `b̂ ≤ 1` absolute series | PASS (`++`, `+−`) | PASS (`++`, `+−`) |
| low side | exact (`rcert.py low`, `lowcheck.py`) | PASS | PASS |

Region M generalised: with `α = W₁(0)`, `w″(0) = q₀ + q₁u + q₂u²` (`q₂ = 16W₁(0)` for `+−` from `(1−up²)⁻⁸`),
`γᵢ = −qᵢ/2`, `b̂⁴𝔉 ≥ γ₀b̂/(1+δ_s)³ − C_R/(1−δ_s)⁴ − π₁²/(2a_φ(1−η)(1−δ_s)²)` for `u ≤ u₁ = 1/100`
(`π₁ = [α]₋ + [γ₁]₋/b̂₋²`, `η` absorbs the `u³` and `γ₂` terms, `δ_s = 8√u₁/b̂` for `+−`, `0` for `++`), increasing in `b̂`; for
`u ∈ [u₁, ¼]` the forward term dominates (same as `cert_M.py`). Positive terms in `α`, `γ₁` are dropped (conservative).

**Theorem R [certified]** (M02): for all `m ≥ M`, allowed `J`, both channels, `𝔉 ≥ (9/50)N`. (M01: `𝔉 ≥ (9/100)N`, certified the same way.)
Sharpness check (mutation): for M02 region P with `J < 12`, `δ′ = 199/1000` passes and `δ′ = 1/5` fails (6 polynomials), consistent with
the floating-point minimum `0.19998` at `+−, J = 4, m = M`; the certified `9/50` leaves ≈ 10 % slack for regions M, S.
Relies on FLINT/Arb, Lemmas P, C, M, S and (Q1)–(Q4) of `FLAT_GRAVITON_BOUND.md`, exactly as its Theorem 1.

**Not certified:** (R2)/(R3) as separate statements (they are implied at physical points by Theorem R, but the scaling-limit margin
`2f̂ ≥ 0.2ρ` is only [NE]: min `2f̂/ρ = 0.1999` per unit `k` at `b̂ = 9.04`); the envelope-targeted family §2.3; everything AdS.

## 5. AdS uplift (`ads.py`, `adsT.py`, `bound.py`)

### 5.1 Leading Regge action [NE, conditional on the heavy formula as in ADS4_FUNCTIONAL.md]
`F_R = 2∫ψ₂P_{−½+iΛp}(cosh β) + σ_hu∫ψ₃P + 2aΛ²u²W(β) + aσ_hu³P_{−½}(cosh β)`; envelope `(C₀cosh β/m² + C₁u)·Σ|parts|`; physical points
`τ ≥ Δ_g`, `β ≤ 8`; grids of `ADS4_FUNCTIONAL.md` t2 refined (`m/Λ` step 0.004 on `[1,12]`). The code reproduces the published `F_{g3}` numbers
(`7.38Λ`, `8.03Λ` vs `8.05Λ`, `+−: 1.21Λ` at `Δ_g = 100`). Last negative `m/Λ` (`++ / +−`; `0` = positive on all physical points):

| functional | `Δ_g` | min `F_R/N` | none | `C₁ = 1` | `C₁ = 3` | `C₁ = 10` | `C₀ = 1000, C₁ = 1` |
|---|---|---|---|---|---|---|---|
| F_g3 | 100 | `−1.0·10⁻²` | 7.38 / 0 | 8.03 / 1.21 | 8.73 / 1.80 | 10.05 / 3.54 | 8.08 / 1.24 |
| F_g3 | 1000 | `−1.1·10⁻²` | 7.46 / 0 | 8.09 / 22.1 | 8.78 / 23.1 | 10.09 / 24.6 | 8.09 / 22.1 |
| M02 | 20 | `+0.197` | 0 / 0 | 2.15 / 1.18 | 3.86 / 1.84 | 6.30 / 3.96 | (C₀ envelope dominates) |
| M02 | 100 | `+0.194` | 0 / 0 | 2.16 / 1.18 | 3.92 / 1.84 | 6.68 / 4.00 | 2.30 / 1.21 |
| M02 | 1000 | `+0.194` | 0 / 0 | 2.16 / 1.18 | 3.92 / 1.84 | 6.70 / 4.01 | 2.16 / 1.18 |
| M01 | 100 | `+0.098` | 0 / 0 | 3.20 / 1.20 | 5.11 / 1.87 | 7.84 / 4.65 | 3.35 / 1.23 |
| M01 | 1000 | `+0.098` | 0 / 0 | 3.20 / 1.20 | 5.12 / 1.87 | 7.88 / 4.66 | 3.20 / 1.20 |

So the robust functionals have **no `Δ_g`-dependent near-zero**: the leading Regge action is positive everywhere with relative margin
≈ 0.19 (M02), and the envelope thresholds are `≈ (C₁/0.2)^{1/2}Λ`, i.e. `2.2Λ, 3.9Λ, 6.7Λ` for `C₁ = 1, 3, 10` (M02) — compared with
`8–22Λ`, `9–23Λ`, `10–25Λ` for `F_{g3}`. The minimal margin sits at `m ≈ 3.2–3.4Λ`, `β → 0` (`++`), i.e. in the bulk-point region; the
Regge kernel is used there as in `ADS4_FUNCTIONAL.md` t2 (not justified at `ν ~ m`; see 5.2).

### 5.2 Hybrid model [M] (M02, criterion (A) of ADS4_FUNCTIONAL.md t3)
`𝔉_hyb = 𝔉_flat + [Regge − flat limit at fixed b̂]` term by term, 1.9–2.1·10⁴ physical states per `(Δ_g, h)`, `m ≤ 60Λ`:
min `𝔉_hyb/N_comp = 0.200 (++)`, `0.28 / 0.23 (+−, Δ_g = 100 / 1000)`; last negative `m/Λ`: `C₁ = 1: 1.92 / 1.24`, `C₁ = 3: 3.58 / 1.85`,
`C₁ = 10: 6.14 / 4.06` (`++ / +−`, both `Δ_g`). The `m ≈ Λ`, `J = 4` state that was binding for `F_{g3}` is not binding.

### 5.3 Rescaled smearing (`adsT.py`) [NE]
With smearing scale `Λ_s = Λ/n` (physical states `m ≥ Λ = nΛ_s`), the leading Regge action of M02 with envelope is positive on
**every** physical grid point at `Δ_g = 100` and `1000` for: `n = 2.2, C₁ = 1` (also `C₀ = 100`; `C₀ = 1000` needs `n` slightly larger at
`Δ_g = 100`); `n = 4, C₁ = 3`; `n = 6.8, C₁ = 10` (also `C₀ = 1000`). M01: `n = 3.2, C₁ = 1`. The envelope-targeted functionals with
their own `n` (`T1`: `C₁ = 1, n = 1.75`; `T3`: `C₁ = 3, n = 2.25`; `T10`: `C₁ = 10, n = 4`) are positive with their envelope on every physical
grid point at `Δ_g = 100, 1000` (flat check on `m ≥ nΛ_s`: no violation; min `F_R/N = 0.015–0.13`, smallest at `Δ_g = 1000`, `m ≈ 20–40Λ_s`).

### 5.4 Conditional AdS bound constants
Given (i) heavy positivity of the uplifted functional with smearing `Λ/n` on all heavy states (supported by 5.1–5.3 only [NE]/[M], and
conditional on all heavy-side inputs of `ADS4_FUNCTIONAL.md` §4), (ii) the light side of `ADS4_LIGHT.md` with the same inputs
((F), (R), (P), `k_AdS = k(0) + O(Λ⁻²)`, no other light single traces, `c_T ≫ Δ_g²`), one gets [D]
`|ĝ₃|²Λ⁸ ≤ n⁸·16L_AdS(Λ/n)/k = A′ log Λ + B′ + O(n⁹/Λ)`, `A′ = n⁸A`, `B′ = n⁸(B_AdS − A log n)` (`bound.py`, exact `L_AdS`):

| case | `n` | `A′` | `B′` | value `Δ_g = 100` | `Δ_g = 1000` |
|---|---|---|---|---|---|
| F_g3 if positive down to `m = Λ` (ADS4_LIGHT.md; **not supported**) | 1 | 37.84 | −60.79 | 113.5 | 200.6 |
| F_g3, its `C₁ = 1` threshold | 8.03 / 22.1 | 6.5·10⁸ / 2.2·10¹² | | 7.1·10⁸ | 4.8·10¹² |
| F_g3, `C₁ = 10` | 10.05 (`Δ_g = 100`) | 3.9·10⁹ | −1.5·10¹⁰ | 3.6·10⁹ | (threshold 24.6 there) |
| **M02, `C₁ = 1`** | 2.2 | 2.84·10⁴ | −6.97·10⁴ | **6.3·10⁴** | **1.27·10⁵** |
| M02, `C₁ = 3` | 4.0 | 3.39·10⁶ | −1.04·10⁷ | 5.6·10⁶ | 1.3·10⁷ |
| M02, `C₁ = 10` | 6.8 | 2.37·10⁸ | −8.48·10⁸ | 2.9·10⁸ | 7.9·10⁸ |
| M01, `C₁ = 1` | 3.2 | 4.13·10⁵ | −1.16·10⁶ | 7.7·10⁵ | 1.7·10⁶ |
| M01, `C₁ = 3 / 10` (`n` from 5.1, not re-scanned) | 5.2 / 7.9 | 2.0·10⁷ / 5.7·10⁸ | | 2.9·10⁷ / 6.3·10⁸ | 7.3·10⁷ / 1.8·10⁹ |
| targeted T1, `C₁ = 1` (§2.3, [NE] only) | 1.75 | 227.1 | −531.1 | 525 | 1039 |
| targeted T3, `C₁ = 3` | 2.25 | 870.7 | −2292 | 1780 | 3728 |
| targeted T10, `C₁ = 10` | 4 | 5.56·10⁴ | −1.99·10⁵ | 7.3·10⁴ | 1.9·10⁵ |

In CFT terms (`DICTIONARY_T4.md`, `t₄² + α_T² = 57600|ĝ₃|²/L⁸`): with M02 and `C₁ = 1`, `√(t₄² + α_T²) ≤ 240·(6.3·10⁴)^{1/2}/Λ⁴ = 6.4·10⁻⁴`
at `Δ_g = 100` (`Λ = 98.5`; `F_{g3}` with its own threshold: `0.068`; the unsupported `n = 1` value: `2.7·10⁻⁵`).
The gain over `F_{g3}` is a factor `1.1·10⁴` (`Δ_g = 100`) to `3.8·10⁷` (`Δ_g = 1000`) in `|ĝ₃|²`, i.e. `≈ 10²` to `6·10³` in
`√(t₄²+α_T²)`; the price relative to the (unsupported) `n = 1` statement is a factor `≈ 550` (`Δ_g = 100`) in `|ĝ₃|²`.

## 6. Caveats
* The envelope `C₁Λ²/m²·Σ|parts|` is the model of `ADS4_FUNCTIONAL.md` (statement R); `C₁` is not derived for `⟨TTTT⟩`, and the
  `O(ν²/m²)` corrections need not be bounded by the *component* norm (a bound by the integrand norm `∫|Ψ||D|` would be much weaker at
  large `b̂`). The robust margin protects against any correction that is relatively small compared with each of `B2`, `B3`, `FW`.
* At `m ≲ 3Λ` (where the thresholds now sit) the Regge kernel is outside its validity (`ν ~ m`); the hybrid model is a stated
  interpolation. That the thresholds of 5.1 and 5.2 agree to 10–15 % is an internal-consistency check, not a derivation.
* The flat functional's large-impact-parameter (IR) issues are unchanged (`FLAT_GRAVITON_BOUND.md` §4.3): the flat bound with M02 has the
  same `Δ_IR` structure; in AdS the IR is regulated.
* The LP optima depend on the grid; only the rationalised M01/M02 are certified. The `δ = 0` optimum `A = 26.1` is a grid result.

## 7. Code (`ads4/code/robust/`)
| file | content |
|---|---|
| `rcommon.py` | general weights, low-side constants, flat parts by Gauss–Legendre + Jacobi recursion, scaling limit, AdS Regge transforms |
| `lp.py` | grid, LP (HiGHS) with cutting planes for (R1), (R2), (R3), matter, Lemma-M proxy, exact `L_AdS` objective |
| `pareto.py`, `pareto2.py`, `mtest.py`, `target.py` | §2.1, §2.2, §2.3 |
| `mkfun.py` | rationalisation → `cache/functionals.json` |
| `check.py` | dense floating-point check, flat envelope thresholds |
| `rcert.py` | certificates P, C, M, S, low side for `𝔉 − δ′N` (generalised copies of `flatbound/cert_*.py`) |
| `lowcheck.py` | exact low side from the amplitude by residues |
| `ads.py`, `adsT.py`, `bound.py` | §5 |
| `run_all.sh` | everything (`bash run_all.sh cert` for the certificates only) |
