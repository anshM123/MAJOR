# Result D — conditional statement (consolidated, 18 Nov 2026)

**Status: a conditional statement, not a theorem.** This file consolidates the AdS4 chain into one precisely stated
claim with every hypothesis listed and graded. Sources: `STATUS.md` (chain), `OMEGA_U.md`, `OMEGA_NF.md`,
`SUPERCONVERGENT_UPLIFT.md`, `FLAT_UNSTRIPPED.md`, `ROBUST_FUNCTIONAL.md`, `ADS4_LIGHT.md`, `SPINNING_4PT.md`,
`SPINNING_TTTT.md`, and their referee reports. Nothing here is new computation.

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

| version | functional | `A`, `B_AdS` | supported `n` | `√(t₄²+α_T²)` at `Δ_g = 100` | at `Δ_g = 1000` | hypotheses |
|---|---|---|---|---|---|---|
| **D-strong** | `Ω^U` built on M02 (with flip-free even kernel) | 51.74, −86.30 | 1.28 central / 1.40 conservative | ≈ 8·10⁻⁵ (range 5.4·10⁻⁵ – 1.5·10⁻⁴) | ≈ 9·10⁻⁹ – 1.3·10⁻⁸ | H1–H12, incl. forward terms (H9) |
| **D-safe** | `Ω^{NF,SC}` built on NF02 (no forward terms; superconvergent replacements) | 3022.3, −6552.7 | 2.0–2.1 central / 2.75–2.9 conservative (for `b̂ ≤ 20`) | 3.1·10⁻³ – 1.1·10⁻² | 5.2·10⁻⁷ – 1.8·10⁻⁶ | H1–H8, H10–H12, plus H13 (`b̂ > 20` positive) |

Meaningful only for `Δ_g ≳ 100` (at `Δ_g = 50` contact residues can shift the bound by up to 17 %). Parametrically
`|t₄|, |α_T| ≲ √(log Δ_g)/Δ_g⁴`, the scaling anticipated heuristically in arXiv:1805.07393 §7 and arXiv:2201.06602 (§1, §5),
here with explicit constants and an explicit hypothesis list.

## 3. What is actually established (unconditional or certified)

| ingredient | status | file |
|---|---|---|
| Flat `F_{g3}` positive on all heavy states; flat bound `37.84 log(M/m_IR) − 45.21 + Δ_IR` | **[certified]**, 3 referee rounds | `FLAT_GRAVITON_BOUND.md` |
| JHEP-printed `F_{g3}` not positive (`−6281/315`); v1 (3.3a) misses `5B₅|₀` | **[D]**, re-derived by referee | same |
| Robust flat functionals M01, M02 (margin 9/100, 9/50) | **[certified]** | `ROBUST_FUNCTIONAL.md` |
| Unstripped and stripped fixed-`u` sum rules span the same space; flat bounds identical | **[D]** | `FLAT_UNSTRIPPED.md` |
| No forward-free functional uses smeared `B_{k≥4}`; NF02 from `B₂, B₃` | **[D]**; NF02 **[certified]** except one `+−` corner (float) | `NOFORWARD_FUNCTIONAL.md` |
| Heavy Regge formula `C^± = (2P/m^{2k})(|c_same|² ± |c_opp|²)`, sign of opposite-helicity term | **[DVE]** for `⟨JJJJ⟩` (`k = 4, 6`); **[DVE]** percent level for `⟨TTTT⟩` at one point (`k=6, ν=0, β=1`) | `SPINNING_4PT.md`, `SPINNING_TTTT.md` |
| Graviton pole in `C⁺₂`: `16πG/(ν²+9/4)` | scalar **[DV]**; spinning [D given F, R] | `ADS4_LIGHT.md` |
| Einstein term absent from superconvergent `Ψ_k`, `k = 3…7` | **[D]** given `J₀(G^E) ≤ 2` (`k = 3`: no assumption) | `SUPERCONVERGENT_UPLIFT.md` |
| Uniform positivity of the scalar-kernel `B₂` part (all `M, β`) | **[certified]** | `UNIFORM_POSITIVITY.md` |

## 4. Hypotheses (graded)

| # | hypothesis | grade | evidence / what would settle it |
|---|---|---|---|
| H1 | Two-pair heavy formula for `⟨TTTT⟩` at all needed `k = 2…7`, `ν ≲ m`; MFT constant ¼ | [I] beyond one point | [DVE] at `k=6, ν=0, β=1`; ¼ single-sourced |
| H2 | B′ dictionary: AdS moments → unstripped sum rules (full `d^J_{44}`, "reading B") for gravitons | [I] | [NE] for currents (`SKELETON_TEST.md` §8) |
| H3 | Odd-`k` two-pair moments behave as continued in `k` | [I] | computation for `⟨JJJJ⟩` running |
| H4 | Scalar-skeleton remainder transfers to gravitons | [I] | holds for same helicity of currents with `Δφ ≈ 2.5`; opposite-helicity shift is flat physics (reading B) |
| H5 | Near-gap tolerance `κ_req` and allowance `δ` | [M] | central 0.3/0.03 calibrated, not derived |
| H6 | (F) flat limit of spinning light side; (P) polynomial contacts | [I] | [Q] for scalars (CHMRS) |
| H7 | (R) graviton-pole residue exact | [D-sketch] | given impact-parameter form |
| H8 | `J₀(G^E) ≤ 2` (tree Einstein Regge intercept) | [Q] | KKSDZ; standard |
| H9 | Double-trace/contact content of the even and forward replacements | [I], load-bearing in D-strong (leverage up to 10¹³) | absent in D-safe (no forward terms; smeared leverage 48–200) |
| H10 | EFT power counting of contacts `g_f ≲ c_g G_f log Λ`, `c_g = O(1–30)` | [A] | standard |
| H11 | No other light single traces; `c_T ≫ Δ_g²`; double zeros of spinning kernels at double twists | [A]/[I] | — |
| H12 | Parity-odd heavy exchanges harmless | not treated | — |
| H13 | (D-safe only) positivity for `b̂ > 20` | [I] | needs a fixed-`b̂` resummation or pointwise relative correction ≲ 4·10⁻³ |

## 5. Where it can fail

* Positivity down to the gap (`n = 1`) is **not** supported for either version; the binding states are low-spin (`J ≤ 6`) near
  `m ≈ 1.3–1.4Λ`, where the CHMRS series is not controlled.
* D-strong rests on forward terms whose light-side errors are amplified ~10¹³ (H9).
* D-safe is 40–140× weaker and has an open large-impact-parameter region (H13).

## 6. Shortest route to an unconditional statement

1. H1, H3: `⟨TTTT⟩` moments at more `(k, ν, β, m)` (production code exists; ~2 h/point).
2. H2: reading-B test for gravitons at fixed `b̂` (same code).
3. Near-gap control (H5, H13): a fixed-`b̂` (impact-parameter) resummation of the heavy moments instead of the `1/m` series.
4. H9: analyse double-trace/contact content of the even replacements, or find forward-free functionals with better `A`.

No unified manuscript is written: Claim D is not yet a theorem.
