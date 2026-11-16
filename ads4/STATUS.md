# Status of the Einstein-gravity programme (living summary; updated as results arrive)

Last update: 16 Nov 2026. **Result D — a theorem on the emergence of local Einstein gravity in AdS4 from
CFT3 data — is not established.** This file lists every link of the chain with its status.

## Target

For CFT3s with large `c_T` and all higher-spin single-trace operators above `Δ_gap`: a rigorous bound on the
non-Einstein cubic graviton coupling of the form `|ĝ_3|² ≤ (function vanishing as Δ_gap → ∞)`, with explicit
assumptions and error terms, valid in the curved AdS4 regime. The flat-space analogue is a physics-level
derivation with numerically checked positivity (arXiv:2201.06602 (4.4)).

## Chain

| # | Link | Status | File |
|---|---|---|---|
| 1 | Scalar-kernel moments: `Q_1 > 0`, `Q_2 < 0` in `d = 3` | **proved, refereed** | `POSITIVITY_LEMMA.md` |
| 2 | Monotonicity of the improvement ratio (all `d > 2`) | **proved, refereed** (classical identity) | `POSITIVITY_LEMMA.md` |
| 3 | Kernels for spinning Regge exchange (harmonic analysis) | **proved, refereed** (mostly classical) | `SPINNING_KERNELS.md` |
| 4 | Which kernel multiplies which graviton coupling | derived from the literature's operators; unrefereed | `REGGE_FROM_BLOCKS.md` |
| 5 | Positivity of the graviton `B_2` weight against the AdS4 scalar kernel, all `M`, `β`; flat-space positivity of that weight | **proved, refereed, independently re-certified** | `UNIFORM_POSITIVITY.md` |
| 6 | Normalised heavy-block Regge formula for stress tensors | **narrow result, refereed (sound; labels corrected):** for one spinning Regge pair in `⟨OφOφ⟩` (`O = J, T`; `Δφ = 1.75`), the leading tensor Regge moment per `λ²`, both helicity components, `k = 4, 6`, `ν ≪ m`, matches exact 3D spinning blocks with no fitted constant. **Both pairs spinning, `⟨JJJJ⟩` (parity-even exchanges), refereed in two rounds — established:** `C^± = (2P/m^{2k})(|c_same|² ± |c_opp|²)(1 + O(1/m²))`, absolute constant fixed from exact Fock-space mean-field data (¼ per helicity quadruple; closed forms identified from exact data), sign of the opposite-helicity term derived and independently re-derived by the referee (same helicity ↔ flat `++`, σ=+1); all 81 modes matched end-to-end with no fitted constant; no `O(1/m)` term (all normalisation factors Weyl-even; under Assumption T). Novelty incremental (cf. arXiv:2311.04271, 2102.08160, 1809.05111). **`⟨TTTT⟩`:** graviton helicity map identified exactly; MFT normalisation of the graviton helicity couplings Weyl-even ⇒ no `O(1/m)` for parity-even exchanges — now from our own exact Fock-space `[TT]` data (closed forms [ID]; Weyl-evenness identity proved symbolically for `l=1`, exact at 4 rational points for `l=2`) and Assumption T; `l=2` absolute constant not yet fixed; unrefereed (`SPINNING_TTTT.md`). **Not established:** `⟨TTTT⟩` end-to-end from exact blocks (symbolic route infeasible; numerical jet route being tested), the coupling→helicity map, (mean-field normalisation for one pair now done: `2P/m^{2k}` form with weight ¼ per helicity pair, even normalisation tested against blocks on the lightcone, odd `l=1` sign measured from the OPE; refereed, labels corrected; leading-twist data agree with Sleight–Taronna arXiv:1804.09334), `k = 2, 3` (the moments the graviton sum rules use; now converged for one spinning pair — all `J` cases and `T` helicity-preserved, `R → 1`, `m²(R−1)` converging, no non-analytic terms seen; `T` flip at `k=2` not determined at β≈0.4; unrefereed, referee running, `SPINNING_LOWK.md`), `ν ~ m`, `TT→φφ`, other `Δφ`. Novelty incremental (cf. arXiv:2311.04271, 2102.08160) | `REGGE_MOMENT_TTTT.md`, `SPINNING_E2E.md`, `REFEREE_SPINNING_E2E.md`, `SPINNING_MFT.md` |
| 7 | First subleading (`1/m²`) corrections, incl. the opposite-sign `B_3` term | no `O(1/m)` term (derived, verified); `1/m²` coefficients measured, spin-dependent part `−2l²P′/P` inferred from flat space, confirmed for currents at all tested β, for `T` only as β→0; analytic derivation not done; opposite-sign term not testable in `⟨OφOφ⟩` | `SPINNING_E2E.md` §6 |
| 8 | Positivity of the complete flat-space graviton functional (all components, all heavy states) and the rigorous flat-space bound | positivity of `F_{g3}` (v1 form and g₅-free completion) on all heavy states, both channels incl. the opposite-sign `B_3` term: **computer-assisted proof (Arb), refereed: correct and sound**; novelty incremental (the paper's plots already indicate it). Bound: (i) given the paper's hypotheses H1–H5, `|ĝ₃|²M⁸ ≤ 37.84 log(M/m_IR) − 45.21 + Δ_IR` with `Δ_IR ≥ 0` explicit; (ii) for `m_IR ≤ M/20000`, under an eikonal hypothesis E(λ,η) restricted to a certified large-impact-parameter region Ω_ε (Lemma N), `≤ 37.8387 log(M/m_IR) − 40.790` (λ=1, η=0); for `M/20000 < m_IR ≤ M/100` only an estimate. Earlier infrared hypothesis H-IR **withdrawn** (referee: contradicts the graviton pole). Findings: `+−` positivity nearly marginal (`m¹⁰F = 19/315`); (3.3a) of v1 omits a `5B₅|₀` term (probable misprint; paper's figures match the completed functional); the JHEP-printed functional is **not** positive (`−6281/315` at `+−, J=4, m=M`; independently re-derived by referee). **Three referee rounds; established** (last write-up fix R3.1 applied as prescribed by the referee). This is a flat-space statement; it is not the AdS4 bound | `FLAT_GRAVITON_BOUND.md` |
| 9 | Uniform remainder near the gap (`m ~ M`) in AdS4 | not started | — |
| 10 | Light side in AdS4 (graviton exchange, double-traces, parity-odd) | not started | — |
| 11 | The AdS4 bound | **does not exist** | — |
| — | Dictionary ĝ₃ ↔ CFT3 data (`t₄² + α_T² = 57600 ĝ₃ĝ₃*/L⁸`, `M L = Δ_gap`) and literature comparison (no rigorous Δ_gap bound on `t₄` in CFT3 found; closest: heuristic `t₄ ≲ ln Δ/Δ⁴` in magnitude, arXiv:1805.07393 §7) | derived and numerically cross-checked against three published cases; **unrefereed** | `DICTIONARY_T4.md` |

## Side results (independent, refereed)

AdS3 covariant entropy inequalities (repository root); ABJM `N^{3/2}` law, finite-`N` bounds, and
`F/c_T → π²/64` (`ABJM_THEOREM.md`); conditional lightcone theorem (`LIGHTCONE_THEOREM.md`); conditional `d = 3`
energy-correlator sum rule and its kernels (`EINSTEIN_PROGRAM.md` §3, `KERNELS_D3.md`).

## Known code issues

* `ads4/code/e2e/cblock.c`: memory corruption when truncation level `N < J` (worked around in `ads4/code/lowk/`; not fixed).

## Retracted claims (recorded so they are not repeated)

* "No positive functional exists in the helicity-flip channel" — wrong channel assignment.
* "The energy-correlator sum rule cannot give a gap-dependent bound" — non sequitur.
