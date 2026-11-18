# Result D — conditional statement (consolidated; revised 18 Nov 2026 after the batched foundation review)

**Status: a conditional statement, not a theorem.** Revised after `REFEREE_FOUNDATION.md` (hostile review of this file and all
inputs added on 18 Nov). The referee independently reproduced every load-bearing number (exact skeleton moments, even-kernel
weights, `N^{K*}`, `A`, `B_AdS`, all six bounds, the 240 dictionary factor, parity-odd interference freedom, bulk-point density ratios),
and found that **the previous headline ("D-safe positive down to the gap, n = 1") was overclaimed**: it holds in a scalar-skeleton
model for a functional whose graviton version is not yet defined. This revision states what is shown at which level.

Sources: `STATUS.md` (chain), `OMEGA_U.md`, `OMEGA_NF.md`, `SUPERCONVERGENT_UPLIFT.md`, `FLAT_UNSTRIPPED.md`, `ROBUST_FUNCTIONAL.md`,
`NOFORWARD_FUNCTIONAL.md`, `ADS4_LIGHT.md`, `SPINNING_4PT.md`, `SPINNING_TTTT.md`, `SPINNING_TTTT_EXT.md`, `EXACT_SKELETON_NEARGAP.md`,
`CERTIFIED_NEARGAP.md`, `DIRECT_GRAVITON_NEARGAP.md`, `PARITY_ODD_HEAVY.md`, `OBJECT_IDENTITY.md`, `FUNCTIONAL_REDEFINITION.md`,
`SPINNING_MELLIN.md`, `DICTIONARY_T4.md`, and their referee reports. Nothing here is new computation.

Labels: **[certified]** computer-assisted proof (Arb/exact); **[D]** derived; **[DVE]** derived and checked end to end on exact
blocks; **[NE]** numerical evidence; **[M]** a modelling convention; **[I]** inferred/assumed; **[A]** standard assumption; **[Q]** quoted.

## 1. Setting

A unitary CFT₃ (parity-even or parity-violating) with stress tensor `T`, central charge `c_T`, and a **twist gap**: every single-trace
operator of spin `J > 2` has twist `Δ − J ≥ Δ_g − 2` (a dimension gap is not enough: it allows states below `Λ` at large `J`). No other
light single traces (H11; note this excludes all known holographic CFT₃s, ABJM included, which have light Kaluza–Klein towers). Large `c_T`
with `c_T ≫ Δ_g²` (tree level). Write `⟨TTT⟩` in the basis `(c_T, t₄, α_T)`; in the bulk `t₄² + α_T² = 57600 |ĝ₃|²/L⁸`
(`DICTIONARY_T4.md`; reproduced by the foundation referee, including the ECG and quasi-topological cross-checks). `Λ = √((Δ_g−2)(Δ_g−1))`.

## 2. Claim D (conditional)

Under hypotheses H0–H17 below, with a heavy-positivity threshold `n`,

`|ĝ₃|² (Λ/n)⁸ ≤ A log(Λ/n) + B_AdS + O(n/Λ)`,

i.e. `√(t₄² + α_T²) ≤ 240 √(A log(Λ/n) + B_AdS) · (n/Λ)⁴` (the tables use the exact `p`-integral, not the log form).

| version | functional | `A`, `B_AdS` | `n` | `√(t₄²+α_T²)` at `Δ_g` = 100 / 300 / 1000 | hypotheses |
|---|---|---|---|---|---|
| **D-safe** | NF02-based, no forward terms; even `k` := `C`-type dispersive functionals (`FUNCTIONAL_REDEFINITION.md`), odd `k` := `Ψ_k/N_k` | 3022.3, −6552.5 | **1 in the scalar-skeleton model** (certified on a grid at Δ_g = 100; scanned at 300, 1000) | 2.21·10⁻⁴ / 3.13·10⁻⁶ / 2.89·10⁻⁸ | H0–H4, H6–H8, H9b, H10–H17 |
| **D-strong** | `Ω^U` on M02 | 51.74, −86.30 | 1.08 / 1.05 / 1.01 in the same model (object-sensitivity range 1.00–1.10 at Δ_g = 100) | 4.2·10⁻⁵ / 5.3·10⁻⁷ / 4.1·10⁻⁹ (up to ×1.4 object/parity sensitivity) | as D-safe, plus H9a |

**What `n = 1` means.** It is the threshold of the near-gap heavy action computed with (i) the flat B′ bulk-point density, which carries
≈ 98 % of the action near the gap, plus (ii) the exact scalar-skeleton remainder (≤ 0.02 `N_s`). Item (ii) is exact; item (i) for gravitons
at `x = ν/m ~ 1` is hypothesis H2/H4. So `n = 1` is **not** shown for gravitons.

Meaningful only for `Δ_g ≳ 100`. Parametrically `|t₄|, |α_T| ≲ √(log Δ_g)/Δ_g⁴`. **Novelty:** the scaling and the AdS₄/CFT₃ translation
with `log(M R_AdS)` are in arXiv:2201.06602 §5 (CLPS), whose footnote functional is of the same forward-free `B₂ + B₃` class as NF02;
D-safe's `A` is ≈ 121× CLPS's optimum (≈ 11× weaker bound). arXiv:2311.04271 lists Regge positivity and the light-side dictionary as
open steps; this chain assumes rather than closes them. The referee found no paper deriving an explicit gap bound on ⟨TTT⟩ in CFT₃ with
explicit constants; the explicit conditional bound and the printed 240 dictionary coefficient appear to be new.

## 3. What is established

### 3a. Flat space and kinematics (unconditional or certified)

| ingredient | status | file |
|---|---|---|
| Flat `F_{g3}` positive on all heavy states; flat bound `37.84 log(M/m_IR) − 45.21 + Δ_IR` | **[certified]**, 3 referee rounds | `FLAT_GRAVITON_BOUND.md` |
| JHEP-printed `F_{g3}` not positive (`−6281/315`); v1 (3.3a) misses `5B₅|₀` | **[D]** | same |
| Robust flat functionals M01, M02 (margin 9/100, 9/50) | **[certified]** | `ROBUST_FUNCTIONAL.md` |
| Unstripped and stripped fixed-`u` sum rules span the same space | **[D]** | `FLAT_UNSTRIPPED.md` |
| No forward-free functional uses smeared `B_{k≥4}`; NF02 from `B₂, B₃` | **[D]**; NF02 **[certified]** as a flat functional except one `+−` corner (float) | `NOFORWARD_FUNCTIONAL.md` |
| Uniform positivity of the scalar-kernel `B₂` part (all `M, β`) | **[certified]** | `UNIFORM_POSITIVITY.md` |
| Even-kernel circle-mode weights `(1/36, 1/16, 1)`; `K*`, `𝔖_ff` flip-free; `N₄^{K*}`, `N₆^{K*}` | **[D]**, reproduced by referee | `SPINNING_TTTT_EXT.md` §6, `EXACT_SKELETON_NEARGAP.md` §9 |
| Parity-odd heavy exchanges: no `λ_Eλ_O` interference (parity invariance of all kernels) | **[D]**, reproduced by referee | `PARITY_ODD_HEAVY.md` |

### 3b. Scalar-skeleton model (exact or certified, but scalar)

| ingredient | status | file |
|---|---|---|
| Exact scalar-skeleton moments `C_{k,ν}` near the gap (Mellin/Mack + Arb Barnes lemma) | **[D]**; independently reproduced by the referee to ≤ 7·10⁻⁶ (Δ_g=100), 5·10⁻⁸ (300), 5·10⁻⁷ (1000) | `EXACT_SKELETON_NEARGAP.md` |
| D-safe action positive on the Δ_g = 100 grid (m/Λ ∈ [1.01, 3], J ≤ 10, both channels; Δφ = 3 and 2), 64-node `p`-rule | **[certified]** for that model (384/384 cases for the 64-node rule; min F/N_s +0.227 for Δφ=3, +0.072 for Δφ=2). **Exact `p`-integral certified** (Bernstein-ellipse Gauss–Legendre error bound) on 400/402 cases incl. J = 11–16 at m/Λ ∈ [1.12, 1.3]; fails only at ++, J=0, m = 1.01Λ (no admissible ellipse; float value +0.997 N_s). Continuity in `m` between grid points [NE] (~30 core-h to certify) | `CERTIFIED_NEARGAP.md` |
| Even-`k` objects `C_{k,ν}` meet: same leading action and B′ density, Einstein-free for `k ≥ 4`, no forward terms, Mellin-exact | **[D]** for external scalars | `FUNCTIONAL_REDEFINITION.md` |

### 3c. Spinning checks (scope-restricted)

| ingredient | status | file |
|---|---|---|
| Heavy formula `C^± = (2P/m^{2k})(|c_same|² ± |c_opp|²)(1+O(1/m²))`, `⟨JJJJ⟩` | **[DVE]** `k = 4–7` | `SPINNING_4PT.md`, `SPINNING_ODDK.md` |
| Heavy formula, `⟨TTTT⟩` | **[DVE] only for β = 1, ν ≤ 3, even J, parity-even exchanges, k ≥ 4, kernels `K*`/`𝔄`**; `k = 2` [NE + I cleaning]; `𝔖_ff` unresolved at low J for k ≤ 4; ¼ not tested end to end | `SPINNING_TTTT.md`, `SPINNING_TTTT_EXT.md`, `DIRECT_GRAVITON_NEARGAP.md` |

Conditional items previously listed here (graviton pole in `C⁺₂`; Einstein term absent from `Ψ_k`) are moved to the hypotheses (H6, H8, H16).

## 4. Hypotheses (graded; referee grades adopted)

| # | hypothesis | grade | evidence / what would settle it |
|---|---|---|---|
| **H0** | **Existence of the spinning (⟨TTTT⟩) functional as defined**: a `C`-type dispersive `Ĉ_k` for even `k` (incl. `Ĉ₂`, which carries the graviton pole and the whole light side) with the properties of `FUNCTIONAL_REDEFINITION.md` (merges the former H14′) | **[I], load-bearing** | sketched (meromorphic-weight light-ray commutators + CHMRS η-inversion), not built. Alternative: keep `Ψ_k` and recompute its flat limit, certify it, recompute its light side |
| H1 | Two-pair heavy formula for `⟨TTTT⟩` at all needed `k`, `ν ≲ m`, β, `J`; MFT constant ¼ | **[NE]** as used; [DVE] in the restricted scope of §3c | needs `ν ~ m` (= H2), β ≠ 1, odd `J` (= H12′), a clean `k = 2` production rerun, ¼ end to end. Finite-`m` sensitivity: scaling even moments by `1 + c₂/m²` (`c₂ ≤ 268`) moves no threshold [NE] |
| H2 | B′ dictionary for gravitons (reading B, full `d^J_{44}`) | [NE] at `x ≤ 0.1`, one `b̂`, even `J`; **[I] at `x ~ 1`, where it is load-bearing** (≈ 98 % of the near-gap action) | spinning moments at `x ~ 1` (`SPINNING_MELLIN.md`: position space cannot reach; Mellin route needs H0-type objects) |
| H3 | Odd-`k` moments | **[M]** (continuation/interpolation in `k`); `⟨TTTT⟩` k = 5, 7 under `𝔄⊗𝔄` [DVE] in the §3c scope | — |
| H4 | Scalar-skeleton remainder transfers to gravitons | [I]; `++` evidence at `x ≲ 0.1` covers ≤ 0.4 % of the norm | as H2 |
| H5 | Near-gap tolerance convention | not needed for the skeleton part | — |
| H6 | (F) flat limit of spinning light side; (P) polynomial contacts; graviton pole `16πG/(ν²+9/4)` in `Ĉ₂` | [I] (scalar: [DV]/[Q]) | — |
| H7 | (R) graviton-pole residue exact | [D-sketch] | — |
| H8 | `J₀(G^E) ≤ 2` | [Q] | KKSDZ |
| H9a | Forward-term double-trace/contact content (D-strong only) | [I], leverage up to 10¹³ | — |
| **H9b** | (P) for the smeared even replacements in D-safe | **[I]**, leverage up to 200 (Δ_g = 100) → 1.1·10⁴ (Δ_g = 1000) | analysable once even `k` are meromorphic-weight `C`-type objects (H0) |
| H10 | EFT power counting of contacts | [A] | — |
| H11 | No other light single traces; `c_T ≫ Δ_g²`; double zeros of spinning kernels at double twists | [A]/[I]; excludes known holographic CFT₃s | — |
| H12′ | Heavy formula (and no `O(1/m)`) for parity-odd exchanges, incl. every odd-`J` `+−` state (the D-safe minimum at all three `Δ_g` is the `+−`, `J = 5` parity-odd state; 0.227 / 0.206 / 0.202) | **shape [DVE]** (`PARITY_ODD_HEAVY.md` §4): parity-odd `⟨TTO⟩` built in the Taylor-mode ⟨TTTT⟩ code (12-dim odd (2,2,J) space; CHL odd structure decomposed to 1e-59); at β = 1, J = 35 (m = 68.1) and J = 19 (m = 37.4): mixed waves W[E,O] vanish identically (pointwise parity-odd to 2e-12), and W[O,O]/W[E,E] is one constant for all k = 2…7, ν ≤ 3 and all kernels (spread ≤ 2e-5 / 2e-4 at k ≥ 5), i.e. the odd exchange has the opposite-helicity heavy-formula shape incl. its 1/m² corrections. **Normalisation N = 3.30 / 3.77 relative to CHL-normalised structures [I], unexplained** (expected 1; varies 14 % between the two masses); positive, so cannot change a sign; if its variation hid an O(1/m) term it would be |c₁| ≈ 14, below the tolerance on odd-J `+−` states (D-safe |c₁| ≤ 56 / 140 / 455 at Δ_g = 100/300/1000; D-strong ≤ 182 at 100) — no threshold moves [NE]. Not run: even-J parity-odd (`++`), other β | parity-odd MFT (Fock) data to fix N; even-J odd structure; other β; small-`m` odd waves at k = 3, 4 |
| H13′ | Heavy positivity off the scanned/certified set | certified grid at Δ_g = 100 incl. the exact `p`-integral (§3b; except m = 1.01Λ, J=0) and J ≤ 16 near the gap. D-safe at large `b̂` (`NF02_LARGE_B.md`): hybrid action (flat + leading AdS Regge) positive for m ∈ [1.5, 10]Λ, b̂ ∈ [20, 200] (Δ_g = 300, 1000; b̂ ∈ [40, 200] at 100) and b̂ ≤ 1000 at sample masses, min F/N_s 0.93 [NE]; exact skeleton remainder on 56 states up to b̂ ≈ 100, |S|/F_hyb ≤ 131/m², min F/N_s 0.95 [NE]; large-b̂ two-term envelope (AdS term A₁M₁/Λ² dominates for b̂ ≳ 4Λ^{2/3}) [D + NE]; **[I] (calibrated) for the remainder at b̂ > 100 and off the sample**; [NE] between grid points, m ∈ [Λ, 1.01Λ], J ≥ 18 near the gap | certify M₁ > 0 (interval arithmetic); derive the envelope remainder; exact remainder at J ≥ 140 (chunked Mack matrices); continuity in `m` |
| H14 | Object identity (`C` vs `Π̂` vs `Ψ_k/N_k`) | for D-safe even `k`: by definition (redefinition) at scalar level, else part of H0; for odd `k` and D-strong: [I], bounded in a bulk-point model [M] | `OBJECT_IDENTITY.md` |
| **H15** | `|ĝ₃|²` absent from the light side of `Ψ₅, Ψ₆, Ψ₇` | [I] | — |
| **H16** | I14 convergence/swappability of the light-ray sum rules; P1(`j′ = 0`); the untracked Jacobian `j_k(ν)`; Einstein term absent from `Ψ_k` (derived given H8, RB₂, microcausality, I14) | [I] | — |
| **H17** | `k_AdS − k(0)` correction to the `|ĝ₃|²` coefficient | [D-sketch]; ≈ 4 % at Δ_g = 100 | — |

## 5. Where it can fail (in order of weight)

1. **H0.** The graviton functional the claim uses is not constructed; its flat limit, light side and contact content are inferred.
2. **H2/H4 at `x ~ 1`.** Near the gap almost all of the action is the flat dictionary at `x → 1`; gravitons are tested only at `x ≲ 0.1`.
3. **H12′.** The binding D-safe state is a parity-odd exchange never tested end to end.
4. **H13′ coverage** for D-safe at large `b̂` and larger `Δ_g`.
5. **H9b** (D-safe) and **H9a** (D-strong) contact content.

## 6. Shortest route (from `REFEREE_FOUNDATION.md` §5, adopted)

1. **Fix the object (H0).** Construct the spinning `C`-type `Ĉ_k` (even `k`, `K*`), or keep `Ψ_k` and recompute its flat limit, certify it,
   and recompute `k(0)`, `A`, `B_AdS` and contacts.
2. **H12′.** One end-to-end parity-odd `+−`, odd-`J` `⟨TTO⟩` point.
3. **H13′ for NF02.** Exact-skeleton scans at Δ_g = 300, 1000, `m ∈ [1.5, 10]Λ`, `b̂ ∈ [20, 200]`, plus a large-`b̂` envelope; extend the
   certificate (continuity in `m`; quadrature error now certified at Δ_g = 100).
4. **H9b.** (P) for the even replacements (becomes analysable with step 1).
5. **H2/H4 at `x ~ 1`.** The real near-gap content for gravitons.

No unified manuscript is written: Claim D is not yet a theorem.
