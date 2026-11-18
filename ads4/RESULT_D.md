# Result D — conditional statement (consolidated, 18 Nov 2026; revised after exact near-gap moments and ⟨TTTT⟩ extension)

**Status: a conditional statement, not a theorem.** This file consolidates the AdS4 chain into one precisely stated
claim with every hypothesis listed and graded. Sources: `STATUS.md` (chain), `OMEGA_U.md`, `OMEGA_NF.md`,
`SUPERCONVERGENT_UPLIFT.md`, `FLAT_UNSTRIPPED.md`, `ROBUST_FUNCTIONAL.md`, `ADS4_LIGHT.md`, `SPINNING_4PT.md`,
`SPINNING_TTTT.md`, `SPINNING_TTTT_EXT.md`, `EXACT_SKELETON_NEARGAP.md`, and their referee reports. Nothing here is new
computation. **The two newest inputs (`EXACT_SKELETON_NEARGAP.md`, `SPINNING_TTTT_EXT.md`) are unrefereed**; they will be
reviewed in the batched foundation review.

Labels: **[certified]** computer-assisted proof (Arb/exact); **[D]** derived; **[DVE]** derived and checked end to end on exact
blocks; **[NE]** numerical evidence; **[M]** a modelling convention; **[I]** inferred/assumed.

## 1. Setting

A unitary, parity-even-or-odd CFT₃ with stress tensor `T`, central charge `c_T`, and a gap `Δ_g` to all single-trace operators of
spin `J > 2` (and no other light single traces). Large `c_T` with `c_T ≫ Δ_g²` (tree level). Write
`⟨TTT⟩` in the basis `(c_T, t₄, α_T)` (Einstein + parity-even `t₄` + parity-odd `α_T`); in the bulk,
`t₄² + α_T² = 57600 |ĝ₃|²/L⁸` (`DICTIONARY_T4.md`, unrefereed; three published cross-checks). `Λ = √((Δ_g−2)(Δ_g−1))`.

## 2. Claim D (conditional)

Under hypotheses H1–H12 below, with a heavy-positivity threshold `n` (§4),

`|ĝ₃|² (Λ/n)⁸ ≤ A log(Λ/n) + B_AdS + O(n²/Λ²)`,

i.e. `√(t₄² + α_T²) ≤ 240 √(A log(Λ/n) + B_AdS) · (n/Λ)⁴`, with two versions:

| version | functional | `A`, `B_AdS` | `n` (exact skeleton moments, scanned states) | `√(t₄²+α_T²)` at `Δ_g` = 100 / 300 / 1000 | hypotheses |
|---|---|---|---|---|---|
| **D-safe (headline)** | `Ω^{NF,SC}` built on NF02 (no forward terms; superconvergent replacements) | 3022.3, −6552.7 | **1** at all three `Δ_g` (min margin `F/N_s ≈ +0.20`) | **2.2·10⁻⁴ / 3.1·10⁻⁶ / 2.9·10⁻⁸** | H1–H4, H6–H8, H10–H12, H13′ |
| **D-strong** | `Ω^U` built on M02 (flip-free even kernel) | 51.74, −86.30 | 1.08 / 1.05 / 1.01 (1.09 / 1.06 / 1.01 with `δ = 0.03`) | 4.2·10⁻⁵ / 5.3·10⁻⁷ / 4.1·10⁻⁹ | as D-safe, plus H9 (forward terms) |

Previous values (truncated `1/m` series, tolerance convention): D-strong `n = 1.28` central / 1.40 conservative, 8·10⁻⁵ at
`Δ_g = 100`; D-safe `n = 2.0–2.9`, 3.1·10⁻³ – 1.1·10⁻² at `Δ_g = 100`. D-safe is now chosen as the headline because it avoids the
forward-term hypothesis H9 (leverage up to 10¹³) at a cost of a factor ≈ 5–7 in the bound.

Meaningful only for `Δ_g ≳ 100` (at `Δ_g = 50` contact residues can shift the bound by up to 17 %). Parametrically
`|t₄|, |α_T| ≲ √(log Δ_g)/Δ_g⁴`, the scaling anticipated heuristically in arXiv:1805.07393 §7 and arXiv:2201.06602 (§1, §5),
here with explicit constants and an explicit hypothesis list.

## 3. What is actually established (unconditional, certified, or end-to-end checked)

| ingredient | status | file |
|---|---|---|
| Flat `F_{g3}` positive on all heavy states; flat bound `37.84 log(M/m_IR) − 45.21 + Δ_IR` | **[certified]**, 3 referee rounds | `FLAT_GRAVITON_BOUND.md` |
| JHEP-printed `F_{g3}` not positive (`−6281/315`); v1 (3.3a) misses `5B₅|₀` | **[D]**, re-derived by referee | same |
| Robust flat functionals M01, M02 (margin 9/100, 9/50) | **[certified]** | `ROBUST_FUNCTIONAL.md` |
| Unstripped and stripped fixed-`u` sum rules span the same space; flat bounds identical | **[D]** | `FLAT_UNSTRIPPED.md` |
| No forward-free functional uses smeared `B_{k≥4}`; NF02 from `B₂, B₃` | **[D]**; NF02 **[certified]** except one `+−` corner (float) | `NOFORWARD_FUNCTIONAL.md` |
| Heavy Regge formula `C^± = (2P/m^{2k})(|c_same|² ± |c_opp|²)(1+O(1/m²))` | **[DVE]** `⟨JJJJ⟩` (`k = 4, 6`; odd `k = 5, 7`); **`⟨TTTT⟩`: [DVE] at β = 1, `k = 4–7`, `ν = 0, 1, 3`, both helicities, under the odd kernel `𝔄⊗𝔄` and leading modes (m = 39, 70, 101; `|c₁| ≤ 0.04`); same helicity also `k = 2, 3` and even kernels** | `SPINNING_4PT.md`, `SPINNING_ODDK.md`, `SPINNING_TTTT.md`, `SPINNING_TTTT_EXT.md` |
| Graviton pole in `C⁺₂`: `16πG/(ν²+9/4)` | scalar **[DV]**; spinning [D given F, R] | `ADS4_LIGHT.md` |
| Einstein term absent from superconvergent `Ψ_k`, `k = 3…7` | **[D]** given `J₀(G^E) ≤ 2` (`k = 3`: no assumption) | `SUPERCONVERGENT_UPLIFT.md` |
| Uniform positivity of the scalar-kernel `B₂` part (all `M, β`) | **[certified]** | `UNIFORM_POSITIVITY.md` |
| Exact scalar-skeleton moments near the gap (Mellin/Mack + Arb Barnes lemma), no `1/m` truncation | **[D]** method; Arb vs double `|ΔF| ≤ 2·10⁻⁷`; quadrature stable to `10⁻⁵N_s` | `EXACT_SKELETON_NEARGAP.md` |

## 4. Hypotheses (graded)

| # | hypothesis | grade | evidence / what would settle it |
|---|---|---|---|
| H1 | Two-pair heavy formula for `⟨TTTT⟩` at all needed `k = 2…7`, `ν ≲ m`; MFT constant ¼ | **[DVE] at β = 1, `k = 4–7`, `ν = 0, 1, 3`, both helicities, under `𝔄⊗𝔄`, `𝔖_ff⊗𝔖_ff` (`|c₁| ≤ 0.13`; corrected circle-mode weights) and the flip-free `K* = (D₁+D₂)(D₁−½D₂)` (`|c₁| ≤ 0.05`); `k = 2` (`Ĉ₂`) both helicities under `K*` and `𝔄` [NE]: R = 1.134/1.040/1.021 (same), 1.114/1.030/1.014 (opp) at m = 39/70/101, `m²(R−1)` flat to ~10 %, after cleaning numerical defects of the cached waves (truncation mask, tiny-`t` noise) — cleaning is [I] until a production rerun** | open: `𝔖_ff` at `k = 2` opposite helicity (isotropic-mode noise; hence use `K*` for `Ĉ₂`); β ≠ 1 convergence (β = 0.4 single mass, 𝔄: 1.02–1.06); `ν ~ m` (= H2); ¼ end to end. **Consequence for Claim D: switching the even kernel to `K*` reproduces every threshold and bound exactly, and scaling every even heavy moment by `1 + c₂/m²` with the measured `c₂ ≤ 268` moves no threshold (tightest: D-strong `+−`, J=4, Δ_g=300, m=1.05Λ, margin 0.0408 → 0.0393 N_s > δ)** (`SPINNING_TTTT_EXT.md` §6–7, `EXACT_SKELETON_NEARGAP.md` §9) |
| H2 | B′ dictionary: AdS moments → unstripped sum rules (full `d^J_{44}`, "reading B") for gravitons | **[NE]** | currents (`SKELETON_TEST.md` §8); gravitons: opposite/same double ratio carries 63–98 % of the reading-B term with bounded `m²d_B` (`SPINNING_TTTT_EXT.md`); one point non-discriminating (near a sign change), one unusable (`P_J` near zero); other `b̂` untested |
| H3 | Odd-`k` two-pair moments | [DVE] for `⟨JJJJ⟩` `k = 5, 7`, `⟨TTTT⟩` `k = 5, 7` under `𝔄⊗𝔄`; skeleton odd `k` at `Δ_g ≥ 300` by interpolation in `k` [M] (≤ 2·10⁻³E_k) | `k = 3` at `m ≈ 100` open |
| H4 | Scalar-skeleton remainder transfers to gravitons | [I] — **now the main near-gap gap** | holds for same helicity of currents (`Δφ ≈ 2.5`); opposite-helicity shift is flat physics (reading B). Thresholds use skeletons `Δφ = 2, 3` only |
| H5 | Near-gap tolerance `κ_req`, allowance `δ` | **no longer needed for the skeleton part** (exact moments) | `δ = 0.03` moves D-strong `n` by ≤ 0.01 |
| H6 | (F) flat limit of spinning light side; (P) polynomial contacts | [I] | [Q] for scalars (CHMRS) |
| H7 | (R) graviton-pole residue exact | [D-sketch] | given impact-parameter form |
| H8 | `J₀(G^E) ≤ 2` (tree Einstein Regge intercept) | [Q] | KKSDZ; standard |
| H9 | Double-trace/contact content of the even and forward replacements | [I], load-bearing in D-strong only (leverage up to 10¹³) | absent in D-safe |
| H10 | EFT power counting of contacts `g_f ≲ c_g G_f log Λ`, `c_g = O(1–30)` | [A] | standard |
| H11 | No other light single traces; `c_T ≫ Δ_g²`; double zeros of spinning kernels at double twists | [A]/[I] | — |
| H12 | Parity-odd / indefinite-parity heavy exchanges harmless | **[D] flat side** (f's cuts are diagonal `|c^{+−}|²`, `|c^{++}|²`; all functionals use only f) **and AdS interference-free** (all kernels parity-even, flip responses vanish for `K*`, `𝔖_ff`, `𝔄` ⇒ no `λ_Eλ_O`, no `(c^{++})²` terms); leading formula unchanged [D given H1]; reduced to **H12′** [I]: H1 (and no `O(1/m)`) for the parity-odd `⟨TTO⟩` structure | not tested end to end (needs the ε-structure in the Taylor-mode code, ~2.5 h/mass); tolerance: thresholds unchanged for `|c₁| ≲ 8` (D-strong) / `≲ 19` (D-safe) vs measured parity-even `|c₁| ≤ 0.05` (`PARITY_ODD_HEAVY.md`) |
| H13′ | Heavy positivity outside the scanned set | [NE] inside it, [M] outside | scanned: `J ≤ 10` to `m = 3Λ / 2Λ / 1.5Λ` (Δ_g = 100/300/1000) plus high-`J` samples (`J ≤ 28`, `b̂ ≤ 40–60`), all positive; beyond, `OMEGA_U`'s `κ ≥ 1.3–2.6` [M]. Not certified (sampled states, not intervals) |

## 5. Where it can fail

* **H1 beyond `ν ≲ 3` and β = 1.** The heavy formula is checked end to end for `k = 2–7`, both helicities, at `ν ≤ 3`, β = 1 (`k = 2`
  after cleaning cached-wave defects); the functionals need `ν ~ Λ` (the H2 question) and other β.
* **H4.** Near-gap positivity is computed on scalar skeletons, not on graviton moments.
* **Coverage.** Positivity is on sampled states; the region `b̂ > 60`, `J` in the hundreds, is not covered by exact moments.
* D-strong additionally rests on forward terms whose light-side errors are amplified ~10¹³ (H9).

## 6. Shortest route to an unconditional statement

1. H1: a clean production rerun for `k = 2` (removes the [I] cleaning); `ν ~ m` and β ≠ 1. The even kernel is now `K*` (no effect on any number).
2. H4: near-gap positivity on graviton moments directly (production code exists), at least at the binding states.
3. H13′: Arb-certify near-gap positivity on intervals in `(m, J)`, plus a large-`b̂` envelope.
4. H12′: one end-to-end parity-odd `⟨TTO⟩` point.
5. One batched hostile review of the foundation (this file and its sources).

No unified manuscript is written: Claim D is not yet a theorem.
