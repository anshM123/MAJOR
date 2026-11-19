# S1 and S4 of `REVIEW_ROUTEA_2026-11.md`: the definition of `Ĉ₃`, and near-gap graviton corrections for `J = 11–20`

*19 Nov 2026, fix agent. Unrefereed, not committed. Nothing in `RESULT_D.md`, `STATUS.md`, `H0_SPINNING_FUNCTIONAL.md`,
`H7_I14_LIGHTSIDE.md`, the review or `n3_*` was edited. New code only: `ads4/code/spinfunc/q*_*.py` (+ `q2_queue.sh`); caches
`spinfunc/cache/q*` (git-ignored). Project modules are imported read-only. Where a project module's cache file name was hard-wired,
it was redirected to a `q*` file by monkey-patching (`q2_run.py`, `q4_direct.py`); one exception, noted in §4: a timing test at the
start added two shadow entries `(800, 12)`, `(800, 20)` to `l9_shadow_k2.pkl` (new keys, same engine; nothing overwritten).*

Grades: **[P]** proved, **[P given X]**, **[D]** derived, **[DV]** derived and checked numerically, **[NE]** numerical evidence,
**[M]** model, **[I]** assumed/open.

---

## 0. Summary

1. **S4 — `Ĉ₃` fixed: no de-shift.** The leading-Regge kernel of `Ψ_{3,ν}[o_H q]` (`j′ = 1`, `f = −1`) was computed with v3's split
   representation. Its heavy image is exactly `l`-isotypic, `σ_h c₃(ν)P_l` (`++ ≡ −(+−)`), with no `P_{l±1}` [P given the kernel reading
   (D)] + [DV to `4·10⁻¹³`]; `c₃` has no real zero on `(0.02, 300]` [NE]; the kernel's ν-only factor matches the independent crystal
   integrand (`a₁ = −2.50i` both) [DV]. So `Ĉ₃ := Re[Ψ₃[o_H q]/N̂₃]` with **no** de-shift satisfies (HN₃) exactly; `Ĉ₃|_E = 0`, the flat
   limit (NF02, `ψ₃`), `k(0)`, `A`, `B_AdS` and the `n = 1` bounds are untouched. The chains' scalar-`g` shift on `k = 3` was spurious and is
   removed: `J = 10` diagonal correction at 300 drops 0.072 → 0.021; the odd binding state (`+−`, `J = 5`) moves 0.2072 → **0.2052** (300) and
   0.2023 → **0.2021** (1000).
2. **A defect in the existing even-sector chain** [NE]: the off-diagonal of different `Δ` rows is interpolated with inconsistent signs
   (arbitrary basis-vector signs per row). Aligned on the robust imaginary mixing, the `J ≤ 10` off-diagonals at 300 grow
   (`J = 10`: 0.088 → 0.227) and the whole-correction window shrinks to `t* = (−4.6, +3.9)`; `n = 1` for `J ≤ 10` is unaffected.
3. **S1 — `J = 11–20` computed** (`K₂` rows `Δ = 400, 800, 1600`, refined `x`-grid; validated against the direct RQMC integral at
   `Δ = 400` for `J = 16, 18, 20`: the off-diagonal agrees to 10–35 % at `ν ≈ 30–40` and `≤ 10 %` above).
   * **`Δ_g = 300`: `n = 1` does not survive.** Even `J = 18, 20` states at `m/Λ = 1.06–1.25` have negative per-state eigenvalues in every
     variant (`λ/N_s ≈ −0.4 … −1.3`; direct: −0.69 / −0.49; `J = 18` with a direct row at the state's own mass `1.12Λ`: −0.44);
     `J = 16` is marginal (`K₂` −0.16, direct +0.11 or −0.10 depending on the
     `ν ≈ 20` node); `J = 12, 14` survive (0.46, 0.23). The driver is the real `O(1/m²)` helicity off-diagonal at `x ≈ 0.05–0.15`
     against the strongly cancelling near-gap base action at large `J`; it grows roughly as `J^{3.5}/m²`. The crossing (`J ≈ 15–17`) is in
     the near-gap region (`b̂ ≈ 30`), outside the large-`b̂` analysis; single-row extrapolations put failures at `J = 24, 28` up to
     `m ≈ 1.65–2Λ`, i.e. inside the large-`b̂` region too [M]. Within model M, positivity on the computed `J ≤ 20` set holds only from
     `m ≈ 1.2Λ` (direct) / `1.35Λ` (`K₂`).
   * **`Δ_g = 1000`: `n = 1` survives for all `J ≤ 20`** (min 0.2021, binding `+−`, `J = 5`; `t* = (−5.1, +4.3)`, now set by `J = 20`), and for
     `J = 24, 28` on one row with a hybrid base (0.64, 0.38) [M]; **`J = 32` is negative near the gap from the off-diagonal alone**
     (−0.45 at `1.04Λ`; its `K₂` diagonal is broken) [M]. The crossing, `J ≈ 30–32` at `m ≈ 1.02–1.15Λ` (near-gap wedge, `b̂ ≈ 60`), is
     outside every exact-skeleton set and outside the large-`b̂` analysis.
   * Odd `J = 11, 13`: corrections `≤ 0.033 N_s`, harmless.
4. **Status.** The premise "`n = 1` at `Δ_g = 300`" of Result D is **contradicted** in model M ([NE + M]; the model's graviton input is
   validated non-perturbatively where it matters, including a direct row at the failing state's own mass). At `Δ_g = 1000` `n = 1` is
   [NE + M] for `J ≤ 20` and [M] for `J = 24, 28`, but **the same mechanism reaches the margin at `J ≈ 30–32`** ([M], one row, hybrid base,
   off-diagonal only), so `n = 1` at 1000 is **not established either** once `J > 28` near the gap is included.

---

## 1. S4: the definition of `Ĉ₃` (`q1_c3_regge_kernel.py`, `q1b_c3_crystal_check.py`)

### 1.1 The leading-Regge kernel of `Ψ_{3,ν}[o_H q]` (`j′ = 1`, `f = −1`) [D] — same method as `v3_regge_kernel.py`

Frame and conventions as `H7_I14_LIGHTSIDE.md` §18.1 (`x₁ = −y`, `x₃ = y`, `x₂ = x/x²`, `x₄ = −x/x²`, `r = |x||y| → 0`).

* **Dual side (2, 4) at `x₀ = 0`.** By exact row reduction of the 5×5 embedding determinant (`X₀ = (1,0,0)`, `X₄ = (1, x₂², −x₂)`):
  `ε(Z₂,X₂,X₄,X₀,Z₀) = −2x₂² ε₃(z₂, x₂, z₀)` and `ε(Z₄,X₄,X₂,X₀,Z₀) = +2x₂² ε₃(z₄, x₂, z₀)`.
  With the inversion-referred polarisations `z′ = R(x̂)z` of Lemma V (`R = η − 2x̂x̂/x̂²`, `det R = −1`, `Rx̂ = −x̂`),
  `ε₃(z₂, x̂, z₀) = ε₃(z₂′, x̂, z₀)` and `H₂₀ ∝ z₂′·z₀` (as v3). Hence, up to a constant,
  **`o_H(x₀ = 0) = (z₂′·z₀) ε₃(z₄′, x̂, z₀) − (z₄′·z₀) ε₃(z₂′, x̂, z₀)`**, antisymmetric under `2 ↔ 4`, and
  `A₃(x̂, z₀) = (−z₀·x̂)^{−(1−δ)} 16 q·o/(z₀·x̂)⁴` with v3's `q`.
  `ε₃(e, x, z₀) = w_e·x` with a `z₀`-dependent vector `w_e`; it is handled exactly by the third Taylor slot of v3's transform engine
  (`(w·x)|x|^c(−z₀·x)^{−λ}… → (w·∂_p)F_{c,λ}`).
* **Light-ray side (1, 3).** The `j′ = 1` celestial structure is the `j′ = 0` one times the locally constant sign of
  `ε(z₀, z₁, z₃) ∝ sin(θ − φ_m)` (the "`j′ = 0` × locally constant sign" of `l9_points` / `p3`). Only terms odd in `sin` survive the
  θ-average, with moments `⟨c^{2i}|s|^{2j+1−δ}⟩/⟨|s|^{1−δ}⟩ = (½)_i(1−δ/2)_j/(3/2−δ/2)_{i+j}`. The odd terms carry one
  `t_k = ε_k·t̂ = (w_k·y)/(−z₀·y)`; even `t`-pairs reduce as in v3. Without spin the `j′ = 1` side vanishes identically (checked: 0.0).
* **Radial power.** `f = −1` has no `η`: the kernel is `(|x||y|)^{a−2}K₀`, one power of `r` more than `Ĉ₂`'s `η`-kernel (CLS-D's `−r¹P_l`),
  so the transform uses `b = a + 1 = 8` and no gradient. `b` changes only a ν-only constant (`b = 7` gives the same structure; checked).
* **What is assumed (why [D]):** the same three inputs as §18.1 for `Ĉ₂` (the `Ψ⁺` reading, dropped ν-only constants, Lemma V's
  referral), plus the identification of the `j′ = 1` structure with `sign(ε(z₀,z₁,z₃)) × (j′ = 0)`, which is what the crystal code
  implements.

### 1.2 Heavy image: `l`-isotypic and `σ_h`-exact [P given §1.1] + [DV]

**Proof (equivariance; §18.3's argument without the gradient).** For preserving (weight-0) helicity tensors each side's plane-wave
transform is an `SO(2,1)`-equivariant function of `(p, z₀)` of `z₀`-homogeneity `−λ`. `SO(2,1)` acts simply transitively on
`{(p̂, [z₀])}`, so `Â = c_A(−z₀·p)^{−(1−δ)}`, `B̂ = c_B(−z₀·p̄)^{−δ}`, and `∫dφ(−z₀·p)^{−a}(−z₀·p̄)^{a−1} = 2πP_{−a}(cosh β)` gives
**`L = 2πc_A c_B P_{−½+iν}(cosh β)`: no `P_{l∓1}`, for any (2,4) and (1,3) structure, parity-odd ones included.** The reflection that
fixes `(p̂, z₀)` reverses orientation, so it flips the pseudo-scalar `o` and the `j′ = 1` sign while exchanging the helicities:
`c_A(−h₂) = −c_A(h₂)`, `c_B(−h₁) = −c_B(h₁)`. Hence the `++` (`(h₁,h₂) = (+,−)`) and `+−` (`(+,+)`) images are **exactly opposite**. ∎

**Numbers** (`cache/q1_c3_regge_kernel.txt`; 21 values `ν = 0.05 … 300`, 14 `β` each, window `∝ 1/ν` at large `ν`):

| check | result |
|---|---|
| control: v3's `q²`, `j′ = 0` kernel without `η` through the new code is `l`-isotypic (§18.3 item 1) | `|h_∓|/|h₀| ≤ 8·10⁻¹⁶` |
| `Ĉ₃`: fit of `L(β)` to `h₋P_{l−1} + h₀P_l + h₊P_{l+1}` | residual `≤ 4.2·10⁻¹⁴`, **`(|h₋|+|h₊|)/|h₀| ≤ 3.6·10⁻¹³`** at every `ν` |
| helicity: `h₀(+−)/h₀(++)` | **`−1`** to `≤ 10⁻¹⁵` at every `ν`; `max|L₊₋ + L₊₊|/max|L| ≤ 2·10⁻¹⁶` |
| equivariance: `Â(−z₀·p)^{1−δ}`, `B̂(−z₀·p̄)^{δ}` over 5 celestial points | constant to `2·10⁻¹⁵` / `1.5·10⁻¹⁴` |
| total-parity bookkeeping (`(h₁,h₂) → (−h₁,−h₂)`) | `o q` with `j′ = 1`: `+1` (survives on parity-definite blocks); `q²` with `j′ = 1` and `o q` with `j′ = 0`: `−1` (cancel) |
| zeros of the forward normalisation `c₃ = 2πc_Ac_B` on `ν ∈ [0.02, 300]` (250 nodes; relative to the scalar sides) | none: `|c_A/c_A^{sc}| ≥ 7.8`, `|c_B/c_B^{sc}| ≥ 0.018` (minimum at the `ν → 0` end); both → constants at large `ν` (power `−0.002`, `+0.001`) |

### 1.3 End-to-end check against the crystal integrand (analogue of §18.2 check X) [DV]

`q1b_c3_crystal_check.py`: at the heavy reference `Δ_ref = 2·10⁴` (`J = 4`, `l9` caches, `K₂`) the crystal gives
`Ĉ₃/Ĉ₂ = (N₃/N₂)(ν)·(−2ix²/(2−x²))`; the kernel gives the forward ratio `L₃(0)ρ_B(ν)/(L₂(0)ν²)` (identical ν-only rescalings on both
functionals; `ρ_B` restores the θ-moment normalisation). Fits `log Q = c + p log ν + a₁/ν + a₂/ν² + a₃/ν³` on `ν ∈ [40, 1000]`:

| | `p` | `a₁` | `a₂` |
|---|---|---|---|
| crystal (18 nodes) | `0.0000` | `−0.009 − 2.501i` | `50.6` |
| kernel (24 nodes) | `0.0000` | `0.019 − 2.511i` | `47.5` |

Pointwise `Q(ν)/Q(200)` agree to `9·10⁻⁴` at `ν = 49`, `≤ 10⁻⁵` for `ν ≥ 150` (`cache/q1b_c3_crystal_check.txt`). The `a₂` difference
is the crystal's `K₂` truncation (as for `Ĉ₂` in §18.2: 92 vs 82). This independently confirms the `x₀ = 0` reduction of `o_H`, the
`j′ = 1` light-ray side and the radial power `b = 8` (a wrong `b` shifts `p`; a wrong relative sign in `o_H` changes `a₁`).

### 1.4 Decision: `Ĉ₃` needs **no** de-shift [D + DV]

* `Ψ_{3,ν}[o_H q]` has leading-Regge heavy action **`σ_h c₃(ν)P_l`** with `c₃` free of real zeros, and nothing else (§1.2). So
  **`Ĉ₃ := Re[Ψ_{3,ν}[o_H q]/N̂₃(ν)]`, `N̂₃ = c₃` (heavy-limit, "real positive at `++`" convention), with no de-shift**, satisfies the
  spinning (HN₃) exactly: its heavy leading-Regge action is `σ_h P_l`. This is §R.5 ("`k = 3` needs no inversion") and the "`Ĉ₃` is
  unchanged" of the §R definition update, now derived rather than quoted from CLS-D; it closes H0 §24.1's "[I]: its own shift relation
  is not derived" (the shift relation is trivial) and the review's M6 / H20 ("(HN₃) quoted"): **(HN₃) [P given §1.1 (D)] + [DV]**.
* The `k = 3` de-shift in the chains (`l9_analysis.shift`, `p4_odd_analysis.shift`, kept in `v3_heavy_impact`, `v4_*`) used the scalar
  CHMRS `g(ν)`, i.e. it treated `Ψ₃` as if its Regge image were `gP_{l−1} + g̃P_{l+1}`. That is wrong for `Ĉ₃` (no `η`); it was a
  spurious `O(J²/m²)` term (`2αJ(J+1)` type). **It is removed** (`q3_analysis.py`, `Q3_C3SHIFT=none`, the default).
* **Preserved:**
  * `Ĉ₃|_E = 0` [P given Q2, Q8, F, as H7 §4]: it is a statement about the bare commutator `Ψ₃`; division by the ν-only `N̂₃` and `Re`
    keep a zero. (Either de-shift choice kept it; the choice never touched H7.)
  * Flat limit (NF02, `ψ₃`) [D]: the bulk-point image is the stripped `B^{(1)}_3` density with `ε₃ m⁶Ĉ₃ → 2σP_l` normalisation; a local
    de-shift is `O(1/m²)` at fixed `x`, so neither choice changes the flat limit, `k₀ = 27.48509817` or `k(0) = 1`.
  * Light side: `A = 3022.34`, `B_AdS = −6552.66` and the bounds are unchanged — the light side (§15, H7) was always computed with
    no `k = 3` de-shift; the fix makes the heavy chains agree with it. `δk(0)` is the `Ĉ₂′` one only (`ψ̃₂ − ψ₂`), as quoted.
* **The large-`b̂` hybrid** (`ROUTEA_LARGE_B.md`: `σu∫ψ₃[P_l − J₀]`) now rests on the derived (HN₃), not on a quote. No number changes.

### 1.5 Effect on the existing near-gap results (`J ≤ 10`, odd `J = 5, 7, 9`; `q3_analysis.py`) [NE + M]

Same chain as the review's §4 item 5 (C23H, `K₂`, rows as before, `Ĉ₂′` exact `u_T`), all states of the exact-skeleton set, even and odd
sectors in one pass. `scalar` = the old `k = 3` shift (reproduces the review and `v4_odd_impact` exactly: 0.3049 / 0.4332 / 0.6017,
`J = 10` diagonal 0.0724, off-diagonal 0.0878, binding 0.2072, `t* = (−10.2, +8.3)`); `none` = the fixed `Ĉ₃`.

| `Δ_g = 300` | old (`scalar` on `k = 3`) | **fixed (`none`)** |
|---|---|---|
| combined `n`, binding state, min `λ/N_s` | 1; `+−`, `J = 5`, `1.02Λ`: 0.2072 | **1; `+−`, `J = 5`, `1.02Λ`: 0.2052** (uncorrected 0.2054) |
| `t*` (whole computed correction) | (−10.2, +8.3) | **(−10.2, +8.7)** |
| min `λ/N_s` `J = 6 / 8 / 10` (even) | 0.3049 / 0.4332 / 0.6017 | 0.3086 / 0.4275 / 0.6037 |
| min `F/N_s` `J = 5 / 7 / 9` (odd) | 0.2072 / 0.7103 / 0.4592 | 0.2052 / 0.7059 / 0.4459 |
| `J = 10`: max `|ΔF_hh|/N_h`, max `|ΔF₊₋|/√(N₊₊N₊₋)` | 0.0724, 0.0878 | **0.0212**, 0.0878 |

| `Δ_g = 1000` | old | **fixed** |
|---|---|---|
| combined `n`, binding, min | 1; `+−`, `J = 5`, `1.01Λ`: 0.2023 | **1; `+−`, `J = 5`, `1.01Λ`: 0.2021** |
| `t*` | (−36.8, +39.1) | **(−38.7, +39.3)** |
| `J = 10`: max diagonal, off-diagonal | 0.0065, 0.0204 | 0.0018, 0.0204 |

* (These two tables use the chain as it was, i.e. without the off-diagonal row alignment of §2.2(a); with it, `t*` at 300 is
  `(−4.6, +3.9)` and the `J = 8 / 10` off-diagonals are 0.088 / 0.227, while every `n` and minimum above is unchanged.)
* The fix removes most of the diagonal `J`-growth at `J = 10` (0.072 → 0.021) and leaves the off-diagonal (de-shift-independent) as
  it was. In the odd sector it removes the small positive correction of the binding state (0.2072 → 0.2052, i.e. now 0.0002 *below*
  the uncorrected base). `n = 1` at both gaps is unchanged.

---

## 2. S1: near-gap corrections for `J = 11–20` (`q2_run.py`, `q3_analysis.py`, `q4_direct.py`, `q5_tables.py`)

### 2.1 What was computed [NE]

* **Crystal points** (`K₂`, linear chain of H0 §23.2, the `l9_points` / `p3_odd_action` integrands unchanged; `q2_run.py`):
  * even `J = 12, 14, 16, 18, 20` (2×2, `Ĉ₂` and `Ĉ₃` `o_H`), rows `Δ = 400, 800, 1600`;
  * odd `J = 11, 13` (1×1 parity-odd `+−`; the exact-skeleton state sets contain no odd `J = 15, 17, 19` states), rows `Δ = 400, 800, 1600`;
  * each point with its heavy reference `Δ_ref = 2·10⁴` at the same `ν`; for `J ≥ 14` six extra `Ĉ₂` points at the reference itself
    (`x_ref = 0.3 … 0.8`, `q2_refbasis.py`) — `l8_analysis.basis` needs them for the real helicity basis at `Δ_ref` (`J ≤ 12` had them from l8);
  * x-nodes: the l9 set `0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4 … 0.8`, **plus** `0.35, 0.45, …, 0.75` (all even `J ≥ 12`, all rows) and
    `0.075, 0.125, 0.175` (`J = 16–20`, `Δ = 400`). The l9 spacing `0.1` under-resolves the `P_J` oscillation at `J ≳ 16` (period in `x`
    ≈ `π/((J+½)√(1−x²))` ≈ 0.18 at `J = 20`).
* **Phase:** every record uses `logK − iπ(Δ+J)` (PARITY_ODD_ROUTEA §12), asserted on `Im logK mod 2π` (all records: deviation 0.0);
  for the state/reference ratios on these even-integer rows it is a no-op, as expected.
* **Functional:** `Ĉ₂′` (exact `u_T`) and the fixed `Ĉ₃` (no de-shift, §1.4).
* **Near gap:** every exact-skeleton state (`e1`/`e7`, `Δ_f = 2, 3`; `J = 11–20` only `Δ_f = 3`, `m/Λ = 1.04–1.3`) in one pass; even `J`:
  l9_neargap's 2×2 formula; odd `J`: p6's 1×1 formula. Optional hybrid states (`Q3_EXTRA_M`: flat + leading-Regge curvature, scalar
  remainder 0, as `ROUTEA_LARGE_B`'s `F_hyb`) at `m/Λ = 1.4 … 2.5` for `J ≥ 11` [M], to see where positivity returns.
* The `Δ = 200` row was not computed for `J ≥ 12` (`J²/m = 0.7–2`: outside the `K₂` range, H0 §22.3); see §2.2(c) for why it is not
  needed for the conclusion.

### 2.2 Three defects of the existing even-sector chain (they also affect `J ≤ 10`) [NE]

(a) **Off-diagonal sign mismatch between `Δ` rows.** `l8_analysis.basis` fixes each row's real helicity basis vectors only up to sign
  (`realvec` normalises by the largest component). The off-diagonal `Γ_off` of different rows can therefore come with opposite signs,
  and `l9_neargap` interpolates them in `1/m` as they are, i.e. with a partial cancellation. Diagnostic: the imaginary `O(1/m)` mixing
  `m Im Γ_off` (H0 §22.2) is `Δ`-independent asymptotically, so its overlap with the most asymptotic row must be `±1`. Found (`q3`
  header): `(J, Δ) = (4, 100) −0.99, (8, 400) −1.00, (10, 200) −0.94, (20, 400) −0.96` — flipped rows; all others `+0.91 … +1.00`
  (except the unconverged `Δ = 100` rows of `J = 8, 10`: `+0.66, +0.21`). After alignment the real parts are consistent across rows and
  scale as `1/m` (e.g. `J = 10`, `x = 0.2`: `−2.9 / −1.3 / −0.3` at `Δ = 200 / 400 / 1600`). Effect at `Δ_g = 300` (`J ≤ 10`, fixed `Ĉ₃`):
  `J = 8` off-diagonal `0.037 → 0.087`, `J = 10` `0.088 → 0.227`; whole-correction window `t* = (−10.2, +8.7) → (−4.6, +3.9)`;
  `n` unchanged. This is also why the review's "drop the `Δ = 100, 200` rows" variant showed a *larger* `J = 10` off-diagonal (0.170):
  the dropped `Δ = 200` row had the opposite sign.
(b) **`c_off := mΓ_off` is `m`-independent at fixed `x`.** E.g. `J = 16`, `x = 0.1 / 0.15`: `1.39 / 0.91` (`Δ = 400`), `1.32 / 0.91` (800),
  `1.29 / 0.90` (1600); `J = 20`, `x = 0.1`: `2.20 / 2.25` (400 / 800). So the real off-diagonal is `O(1/m²)` in `F` at fixed `x`, and
  holding `Γ_off` (l9) below the lowest row *under*-estimates it by `m_row/m`. `q3` models `c_off` (interpolated in `1/m`, held outside,
  divided by the state's `m`); the l9 convention is kept as a variant (`Q3_OFFSCALE=0`).
(c) **Resolution:** the added nodes (2.1) change the `Δ_g = 300` minima by `≤ 0.06` (`J = 14–20`); the conclusion below is unchanged.

### 2.3 Validation of `K₂` at large `J` against the direct (RQMC) integral (`q4_direct.py`) [NE]

`n1_run.point` (the H0 §24.2 non-perturbative contour integral, read-only) on the `Δ = 400` row (`J²/m = 0.64 / 0.81 / 1.01` for
`J = 16 / 18 / 20`), state and reference, kind `C2`; relative RQMC errors `2–7·10⁻⁴` for `x ≥ 0.075`, `1.3–1.8·10⁻²` at `x = 0.05`
(`ν ≈ 20`, 26 % of the samples masked: at the edge of the method's range).

Real off-diagonal `Γ_off` (`Ĉ₂`, l9 units, same `a_h`; `cache/q4_direct_compare.txt`):

| `J` | source | `x = 0.075` | `0.1` | `0.125` | `0.15` | `0.175` | `0.2` |
|---|---|---|---|---|---|---|---|
| 16 | direct | 1.47 | 2.92 | 3.11 | 2.21 | 0.95 | 0.03 |
| 16 | `K₂` | 2.25 | 3.48 | 3.41 | 2.28 | 0.86 | −0.12 |
| 18 | direct | −2.53 | −4.14 | −3.68 | −1.87 | −0.05 | 0.77 |
| 18 | `K₂` | −3.56 | −4.76 | −3.89 | −1.78 | 0.15 | 0.91 |
| 20 | direct | −3.81 | −4.96 | −3.40 | −0.83 | 0.79 | 0.77 |
| 20 | `K₂` | −5.07 | −5.54 | −3.43 | −0.60 | 0.99 | 0.80 |

* **The large real off-diagonal is a property of the crystal integral, not of the Laplace truncation**: `K₂` overestimates it by
  10–35 % at `x ≤ 0.1` (`ν ≤ 40`) and agrees to `≤ 0.2` (`≤ 10 %`) at `x ≥ 0.125`. (`Ĉ₃`'s off-diagonal enters with a 10× smaller
  weight; no direct `Ĉ₃` points were computed.)
* Diagonal: `K₂ − direct` up to `≈ 240` c₂ units on `x ∈ [0.075, 0.25]` (RQMC errors 5–40) — the same size as at `m ≈ 98`, `J ≤ 10`
  (H0 §24.4); at `x = 0.05` the direct values (errors 200–370) are not usable.
* `K₃` was not run: the off-diagonal's real part first appears at the last computed order (`M₄`), so `K₂ → K₃` would test exactly
  the 10–35 % that the direct integral already measures, at ≈ 10× the cost per point.

### 2.4 `Δ_g = 300` [NE + M]

Per even `J`, worst state of the exact-skeleton set (`m/Λ = 1.04–1.3`); min `λ/N_s` of the per-state 2×2 (uncorrected base in
parentheses), the largest diagonal `|ΔF_hh|/N_h`, the largest off-diagonal `|ΔF₊₋|/√(N₊₊N₊₋)`, and `t*_J` (factor on the computed
correction at which a state of that `J` first fails). Main = `K₂`, fixed `Ĉ₃`, aligned rows, `c_off` model, all nodes
(`cache/q3_analysis_none_all_main300.txt`).

| `J` | main: min `λ/N_s` (uncorr.) | max diag | max off | `t*_J` | direct `Δ = 400` (`x = 0.05` direct / `K₂`) | l9 conventions (hold `Γ_off`, no alignment) | correction only on `x ≥ 0.1` / `≥ 0.15` |
|---|---|---|---|---|---|---|---|
| 4 | 0.2203 (0.2204) | 0.0003 | 0.0015 | (−139, +156) | — | 0.2203 | — |
| 6 | 0.3086 (0.3091) | 0.0018 | 0.013 | (−39, +33) | — | 0.3086 | — |
| 8 | 0.4083 (0.4342) | 0.0051 | 0.088 | (−7.5, +7.5) | — | 0.4275 | — |
| 10 | 0.6037 (0.6046) | 0.021 | 0.227 | (−4.6, +3.9) | — | 0.6037 | — |
| 12 | 0.458 (0.682) | 0.039 | 0.360 | (−2.7, +2.3) | — | 0.519 | 0.433 / 0.592 |
| 14 | 0.234 (0.791) | 0.064 | 0.653 | (−1.6, +1.4) | — | 0.361 | 0.220 / 0.720 |
| **16** | **−0.161** (0.875) | 0.118 | 1.079 | (−1.0, +0.9) | **+0.107 / −0.104** | +0.060 | −0.007 / 0.875 |
| **18** | **−0.726** (0.892) | 0.191 | 1.576 | (−0.7, +0.6) | **−0.694 / −0.694** | −0.403 | −0.195 / 0.729 |
| **20** | **−1.335** (0.924) | 0.206 | 2.152 | (−0.5, +0.4) | **−0.490 / −1.207** | −0.910 | −0.234 / 0.672 |

Odd `J` (1×1, parity-odd `+−`; no off-diagonal): `J = 11`: 0.778 (0.789), max `|ΔF|/N_s` 0.020; `J = 13`: 0.767 (0.788), 0.033.
The odd corrections grow slowly (0.006 / 0.016 / 0.020 / 0.033 at `J = 7 / 9 / 11 / 13`) and are far inside the base margins.

**Reading.**
* **`n = 1` at `Δ_g = 300` does not survive the computed graviton corrections at `J ≥ 18` (and is marginal at `J = 16`), in every
  variant** — `K₂` or direct on the `Δ = 400` row, aligned or l9 conventions, `c_off` or held `Γ_off`. The binding states are
  `m/Λ = 1.06–1.1` (`b̂ = (2J+1)Λ/m ≈ 31–37`): inside the near-gap set, outside the large-`b̂` analysis (`m ≥ 1.5Λ`).
* **Mechanism.** It is the real `O(1/m²)` helicity off-diagonal (`ΔF₊₋`), not the diagonal. At `m ≈ Λ` the near-gap base action is a
  strongly cancelling `ψ`-integral of the oscillating `P_J`, `d̃^J₄₄` (`|base|` ≈ 1 % of `Σ|ψ₂ w env D_h|` at `J = 16–20`;
  `N_s ≈ 0.12 / 0.08 / 0.035 / 0.019` at `J = 10 / 12 / 16 / 20`), while the off-diagonal (`c_off(x)` peaked at `x ≈ 0.08–0.13`, i.e.
  `ν ≈ 25–45`, of fixed sign) does not cancel; its absolute size grows only slowly (`|ΔF₊₋| ≈ 0.003–0.014 / 0.025 / 0.034 / 0.046` at
  `J = 10 / 12 / 16 / 20`). The ratio grows roughly as `J^{3.5}` (×1.8 → ×1.4 per `ΔJ = 2` from `J = 10` to 20).
* **No transfer needed for `J = 18`: a row at the binding states' own mass.** `Δ = 336` (`m = 334 = 1.12Λ`, `q2_queueH.sh`): `K₂` full grid
  plus direct RQMC (`Ĉ₂`) at `x = 0.075 … 0.2` (`cache/q4_direct_compare_336.txt`). Off-diagonal at `x = 0.075 / 0.1 / 0.125 / 0.15`:
  direct `−3.08 / −5.10 / −4.35 / −1.98`, `K₂` `−4.70 / −6.04 / −4.66 / −1.85`; `c_off(0.1) = −1.70` (336) vs `−1.65` (400), confirming the
  `m`-independence down to the near-gap masses. Near gap with rows 336 + 400 (`cache/q5_tables_300_row336_d*.txt`): `J = 18` at
  `m/Λ = 1.06 / 1.12 / 1.2 / 1.3`: **`−0.85 / −0.55 / −0.17 / +0.11` (`K₂`), `−0.71 / −0.44 / −0.11 / +0.14` (direct)**. At `1.12Λ` this is
  the row's own mass: no `m`-transfer and the decisive `x`-window non-perturbative.
* **What the base is.** At these states the base is the exact flat NF02 action to `≥ 99 %` (curvature `≈ 2·10⁻⁴`, scalar remainder
  `≤ 10⁻⁴`, vs `B₂ = 0.081 / 0.035 / 0.018` at `J = 12 / 16 / 20`, `1.07Λ`): nothing in model M other than the crystal corrections
  themselves can move the verdict.
* **Where it sits.** Switching the correction off below `x = 0.1` keeps `J ≥ 16` negative; switching it off below `x = 0.15` makes every
  `J` positive. So the decisive input is the off-diagonal on `x ∈ [0.05, 0.15]` (`ν ≈ 16–50` for the binding states) — exactly the
  window where the direct integral is valid (`ν ≳ 20–25`) and agrees with `K₂` to 10–35 %. Below `ν ≈ 20` (`x ≲ 0.06`) neither method is
  controlled; that node moves `J = 16` across zero (+0.107 vs −0.104) but not `J = 18, 20`.
* **Extent in `m`** (`q5_tables.py`, hybrid base `[M]` for `m/Λ ≥ 1.4`): the failure at `J = 20` extends to `m ≈ 1.25Λ` (direct) /
  `1.35Λ` (`K₂`); at `m = 1.5Λ` (`b̂ = 27`) `J = 16 / 18 / 20` have `+0.55 / +0.33 / +0.23` (`K₂`), i.e. the margin at the edge of the large-`b̂`
  region is shrinking with `J`.
* **Consequence for the bound (illustration only).** With positivity from `m ≥ nΛ`, the bound uses `16n⁸L_A(Λ/n)`: at `Δ_g = 300`,
  `n = 1.2 / 1.3 / 1.45` give `√(t₄²+α_T²) ≤ 6.3·10⁻⁶ / 8.6·10⁻⁶ / 1.3·10⁻⁵` (vs `3.13·10⁻⁶` at `n = 1`). Even that presumes positivity for all
  `J > 20`, which the growth makes doubtful (§2.6).

---

### 2.5 `Δ_g = 1000` [NE + M]

Same chain, rows `Δ = 800, 1600` (interpolated in `1/m`; the states are at `m = 1020–1300`), all nodes
(`cache/q3_analysis_none_all_main1000.txt`):

| `J` | min `λ/N_s` (uncorr.) | max diag | max off | `t*_J` | l9 conventions |
|---|---|---|---|---|---|
| 10 | 0.5876 (0.5882) | 0.0018 | 0.020 | (−39, +40) | 0.5876 |
| 11 (odd) | 0.4105 (0.4112) | 0.0006 | — | (−5163, +634) | 0.4105 |
| 12 | 0.6794 (0.6819) | 0.0021 | 0.039 | (−26, +25) | 0.6793 |
| 13 (odd) | 0.6059 (0.6050) | 0.0020 | — | (−506, +1662) | 0.6059 |
| 14 | 0.7012 (0.7095) | 0.0009 | 0.063 | (−16, +16) | 0.7079 |
| 16 | 0.7250 (0.7503) | 0.0009 | 0.093 | (−10, +10) | 0.7248 |
| 18 | 0.6828 (0.7512) | 0.016 | 0.142 | (−7.1, +6.2) | 0.6825 |
| 20 | 0.6881 (0.8399) | 0.039 | 0.208 | (−5.1, +4.3) | 0.7719 |

* **`n = 1` at `Δ_g = 1000` survives with all `J ≤ 20` corrected**: combined min `λ/N_s = 0.2021` (binding state unchanged: `+−`, `J = 5`,
  parity-odd, `1.01Λ`); whole-correction window `t* = (−5.1, +4.3)`, now set by `J = 20` (it was `(−38.7, +39.3)` with `J ≤ 10`).
* The off-diagonal is `≈ (m₃₀₀/m₁₀₀₀)² ≈ 0.09–0.10` times its `Δ_g = 300` value at the same `J` (e.g. `J = 20`: 0.208 vs 2.15), as the
  `c_off/m` scaling predicts, and grows ×1.6 → ×1.47 per `ΔJ = 2`.
* Alignment matters here too: for `J = 14` (`Δ = 400, 800`) and `J = 20` (`Δ = 800`) the rows come with the opposite sign to `Δ = 1600`
  (overlaps −0.99/−1.00/−0.99); l9's interpolation partially cancels them (`J = 20`: 0.095 instead of 0.208).

### 2.6 Where the corrections cross the margin, and which region that is [NE + M; extrapolation flagged]

`b̂ = (2J+1)Λ/m`; the large-`b̂` analysis (`ROUTEA_LARGE_B.md`) covers `m ≥ 1.5Λ`; the exact-skeleton near-gap sets cover
`m/Λ ≤ 1.3` (`J ≥ 11`), `≤ 2.0 / 1.5` (`J ≤ 10`, 300 / 1000); `m ∈ (1.3, 1.5)Λ` and all `J > 20` near the gap were covered only by
the review's scalar-skeleton hybrid scan. Extra rows here: hybrid base [M] + computed graviton corrections, `m/Λ = 1.4 … 2.5`
(`J ≤ 20`) and, for `J = 24, 28`, the `Δ = 1600` row only (held in `1/m`; justified by the `m`-independence of `c_off` and `Γ_hh` at fixed
`x`, §2.2(b)) with hybrid states at every twist-allowed `m/Λ ≥ 1.02`.

* **`Δ_g = 300`: crossing between `J = 14` and `J = 18`** (computed): `t*_J = 1.4 / 0.9 / 0.6 / 0.4` at `J = 14 / 16 / 18 / 20`. The
  failing states are at `m/Λ = 1.06–1.25` (`b̂ ≈ 27–37`): **inside the near-gap region, not covered by the large-`b̂` analysis**, and not by
  any graviton-specific evidence before this note. At the region boundary `m = 1.5Λ` the corrected margins fall with `J`
  (`+0.55 / +0.33 / +0.23` at `J = 16 / 18 / 20`, `K₂`); extrapolating that trend (≈ −0.1 to −0.15 per `ΔJ = 2`) the failure enters
  `m ≥ 1.5Λ` (`b̂ ≈ 30–35`, the large-`b̂` region, where only the scalar hybrid was checked) around `J ≈ 24–26` [extrapolation].
* **`Δ_g = 1000`: no crossing up to `J = 28`; the crossing is at `J ≈ 30–32`, near the gap** (`cache/q5_tables_1000_hiJ*.txt`,
  `cache/q5b_offonly_1000.txt`). Worst state per `J` (corrected min `λ/N_s`; off-diagonal ratio):

  | `J` | 16 | 18 | 20 | 24 [M] | 28 [M] | 32 [M] |
  |---|---|---|---|---|---|---|
  | min `λ/N_s` (at `m/Λ`) | 0.725 (1.03) | 0.683 (1.02) | 0.688 (1.02) | 0.640 (1.04) | **0.378** (1.02) | **−0.45** (1.04; off-diagonal only) |
  | same, off-diagonal only (diagonal correction off) | — | — | 0.709 | 0.646 | 0.017 | −0.45 |
  | max `|ΔF₊₋|/√(NN)` | 0.09 | 0.14 | 0.21 | 0.34 | 0.96 | 1.45 |

  * `J = 32`: the `K₂` `+−` *diagonal* on the `Δ = 1600` row breaks down (`|Γ_hh|` up to `2·10⁶`; the full eigenvalue −841 is
    meaningless), while the off-diagonal is regular (`|Γ_off| ≤ 3.1`). With the diagonal correction switched off, `J = 32` is negative for
    `m/Λ ≤ 1.13` (`−0.45 / −0.24 / −0.09` at `1.04 / 1.08 / 1.12`, `+0.14` at `1.2`). At `J = 28` the diagonal correction is positive
    (+0.36 at the worst state), so a similar diagonal at `J = 32` could bring it to ≈ 0: **`J = 32` is marginal-to-failing [M]**, and the
    off-diagonal ratio (×1.5 per `ΔJ = 2`) makes `J ≥ 34` fail unless the diagonal grows as fast. `J = 36` (queue G) was stopped
    unfinished.
  * These states (`m ≈ 1.02–1.15Λ`, `b̂ ≈ 57–63`) are in the near-gap wedge `m < 1.5Λ`, `J > 20`, which only the review's scalar hybrid
    scan covered — **not** in the large-`b̂` analysis (`m ≥ 1.5Λ`), where `J = 32` still has `+0.58`.
* **The same rows transferred to `Δ_g = 300`** (`m ≈ 320–650`; far outside the `Δ = 1600` row's `J²/m`, so [M] only):
  `J = 24` fails up to `m ≈ 1.65Λ` (`−2.45` at `1.08Λ`, `−0.05` at `1.6Λ`, `+0.25` at `1.75Λ`), `J = 28` up to `m ≈ 2.0Λ`
  (`−0.49` at `1.75Λ`, `+0.07` at `2.0Λ`). If even roughly right, the graviton corrections break positivity at `Δ_g = 300` inside the
  large-`b̂` region as well (`b̂ ≈ 30–35` at `m = 1.5–1.75Λ`), i.e. S2's concern is realised there.
* **Effect of the `Ĉ₃` fix at `J ≥ 11`:** the spurious `k = 3` shift roughly doubled the diagonal corrections (`Δ_g = 300`, max
  `|ΔF_hh|/N_h` at `J = 12 / 16 / 20`: `0.109 / 0.225 / 0.435` old vs `0.039 / 0.118 / 0.206` fixed), but the minima move by `≤ 0.03` and
  the verdict is the same with either definition (old: `J = 16 / 18 / 20` at `−0.165 / −0.733 / −1.362`).


## 3. Graded summary

| statement | grade |
|---|---|
| Leading-Regge kernel of `Ψ_{3,ν}[o_H q]` (`j′ = 1`, `f = −1`): `o_H(x₀ = 0)` reduction, `j′ = 1` light-ray side, `b = 8` | [D] (same inputs as H7 §18.1) |
| Its heavy image is `l`-isotypic, `σ_h c₃(ν)P_l`, `++ ≡ −(+−)`, for any (2,4)/(1,3) structure | [P given the §1.1 kernel] + [DV] (`≤ 4·10⁻¹³`) |
| `c₃(ν)` has no zero on `ν ∈ (0.02, 300]`, → constant ratio to the scalar sides at large `ν` | [NE] (250-node scan) |
| Kernel vs crystal: same ν-only factor `N₃/N₂` (`a₁ = −2.50i` both; pointwise `≤ 10⁻³`) | [DV] |
| **`Ĉ₃ := Re[Ψ₃[o_H q]/N̂₃]`, no de-shift; (HN₃) exact** | definition; (HN₃) [P given §1.1 (D)] |
| `Ĉ₃|_E = 0`, flat limit, `k₀`, `k(0)`, `A`, `B_AdS`, bounds at `n = 1` unchanged by the fix | [P given H7 §4] / [D] |
| Effect of the fix on `J ≤ 10`: `J = 10` diagonal 0.072 → 0.021 (300); odd binding 0.2072 → 0.2052; `n = 1` at 300, 1000 unchanged | [NE + M] |
| Off-diagonal row-sign mismatch in the l9 even chain (overlaps `−0.94 … −1.00`); fixed by alignment on `m Im Γ_off` | [NE] (defect found); affects every earlier `t*` at 300 |
| `c_off = mΓ_off` `m`-independent at fixed `x` (`Δ = 400–1600`, `J = 12–20`) | [NE] |
| `K₂` off-diagonal at `J = 16–20`, `Δ = 400` (and `J = 18`, `Δ = 336`) vs direct RQMC: 10–50 % high at `x ≤ 0.1`, `≤ 10 %` at `x ≥ 0.125` | [NE] |
| **`Δ_g = 300`: computed graviton corrections break `n = 1` at even `J = 18, 20` (and `J = 16` in `K₂`), `m/Λ = 1.06–1.25`** | **[NE + M]**, robust to every variant tried (direct / `K₂`, alignment, `c_off` / hold, nodes); `J = 16` sign [I] (`ν ≲ 20` node) |
| driver: real `O(1/m²)` helicity off-diagonal at `x ≈ 0.05–0.15` vs the strongly cancelling near-gap base at large `J` | [NE + M] |
| `n` at `Δ_g = 300` | **not 1** in model M; positivity on the `J ≤ 20` set from `m/Λ ≈ 1.21` (direct) / `1.35` (`K₂`); `J > 20` unknown → **[I]** |
| **`Δ_g = 1000`: `n = 1` with all `J ≤ 20` corrected**, min 0.2021 (binding `+−`, `J = 5`), `t* = (−5.1, +4.3)` | **[NE + M]** |
| `Δ_g = 1000`, `J = 24, 28` positive (one row, hybrid base); `J = 32` negative near the gap from the off-diagonal alone (`K₂` diagonal broken) | [M] |
| `n` at `Δ_g = 1000` including `J > 28` near the gap | **[I]**, with a computed [M] failure at `J ≈ 32` |
| odd `J = 11, 13` corrections (`≤ 0.033 N_s` at 300) | [NE + M] |

## 4. Files and hand-over

**Code** (`ads4/code/spinfunc/`, all new; project modules imported read-only):
* `q1_c3_regge_kernel.py` → `cache/q1_c3_regge_kernel.txt`, `q1_c3_table.pkl` (≈ 20 s). `q1b_c3_crystal_check.py` →
  `cache/q1b_c3_crystal_check.txt` (≈ 1 min).
* `q2_run.py` (`K₂` production; cache names redirected to `cache/q2_pts_<kind>_J<J>_k2_<tag>.pkl` (even), `q2_p3_<kind>_J<J>_<tag>.pkl`
  (odd), shadows `q2_shadow_k2_<tag>.pkl`, `q2_p2shadow_<tag>.pkl`), `q2_refbasis.py` (reference-basis points, tag `_c`), queues
  `q2_queue.sh` (A, B), `q2_queueC.sh` (mid nodes), `q2_queueE.sh` (small-`x` nodes + direct), `q2_queueF.sh` / `q2_queueG.sh` (`J = 24–36`,
  `Δ = 1600`); logs `cache/q2_queue*.log`. `cache/q2x_partial200_*` is an incomplete `Δ = 200` row (queue B was stopped there; not used).
* `q3_analysis.py`: loaders, phase strip, `Ĉ₂′`/`Ĉ₃` de-shifts (`Q3_C3SHIFT=none|scalar`), row alignment (`Q3_ALIGN`), `c_off` model
  (`Q3_OFFSCALE`), direct `Δ = 400` rows (`Q3_DIRECT`), small-`x` cut (`Q3_XCUT`), hybrid extra states (`Q3_EXTRA_M`), combined
  even + odd near gap → `cache/q3_analysis_<shift>_<rows>[...]_<tag>.txt`.
* `q4_direct.py` (direct RQMC via `n1_run.point`) → `cache/q4_n1_<kind>_J<J>.pkl`, `q4_direct_compare*.txt`, `q4_direct_grids*.pkl`.
* `q5_tables.py` (per-`(J, m)` tables) → `cache/q5_tables_<Δ_g>_<tag>.txt`; `q5b_offonly.py` (diagonal correction off, used for
  `J = 32`) → `cache/q5b_offonly_1000.txt`. `q6_final.sh` reruns the main variants (≈ 20 min); the `Δ = 336` row runs are
  `Q4_ROW=336 python3 q4_direct.py analyse 18` and `Q3_DIRECT=0|1 Q3_JEVEN=16,18,20 python3 q5_tables.py 300 16`.
* `cache/q2x_partial*`: incomplete rows (stopped queues), excluded from every loader by name.

**Not done / next steps (in order of value):**
1. **Is the failure a model artefact?** Checked here: at `J = 12 / 16 / 20`, `m = 1.07Λ` (300) the base is the exact flat NF02 action
   `B₂ + B₃` to `≥ 99 %` (curvature `cv₂ + cv₃ ≈ 2·10⁻⁴`, transplanted scalar remainder `≤ 10⁻⁴`, against `B₂ = 0.081 / 0.035 / 0.018`). So
   the model's only non-crystal element is negligible; the failure is "flat NF02 action of a near-gap, large-`J` (`b̂ ≈ 30–40`) state vs
   the `O(1/m²)` graviton helicity off-diagonal" — the same competition as the review's S2 (at large `b̂` the margin is an `O(Λ⁻²)`
   effect, and the graviton-specific `O(Λ⁻²)` term has no scalar analogue). What remains model-dependent: the transfer of the `Δ = 400`
   row to `m ≈ 316–390` (`c_off` is `m`-independent on 400–1600 [NE]; below 400 [M]) and the `ν ≲ 20` node. A direct row at the states'
   own masses (`Δ ≈ 320–340`, `J = 16–20`, `x ≥ 0.07`) would remove the first.
2. `ν ≲ 20` (`x ≲ 0.06` at `Δ = 400`): the n3 curved contour for `J = 16` (the marginal `J`).
3. `J = 32` diagonal at `Δ_g = 1000` (the `K₂` `+−` diagonal on `Δ = 1600` breaks down: use the direct integrator or a `Δ ≈ 3200` row);
   `J = 36` (`q2_queueG.sh` was stopped before it; `q2_refbasis.py 36` and the `q2_run.py` line restart it).
4. Graviton corrections in the wedge and at `m ≥ 1.5Λ` for `J > 20` with proper rows (`Δ` matched to `m`), both gaps (review S2).
5. If `Ω^{A′}` is to be kept, the large-`J` failure suggests modifying the functional (e.g. weights that suppress `ν ≈ 20–50` relative to
   the flat action at `J ≳ 15`), not the chain.

**Rules kept:** ≤ 2 compute processes at any time (analysis runs paused a production process), 1 thread each, ≤ 1 GB each; nothing
committed; no protected file edited.
