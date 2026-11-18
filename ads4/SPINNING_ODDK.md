# Odd-`k` tensor Regge moments with both pairs spinning (`⟨JJJJ⟩`): `C⁻₃, C⁻₅, C⁻₇` end to end

Status: written 18 Nov 2026 by a derivation/numerics agent; **not refereed; not committed**. Code: `ads4/code/fourpt/oddk/` (caches in
`oddk/cache/`, git-ignored). Builds on `SPINNING_4PT.md` (the `⟨JJJJ⟩` partial waves, the factorised formula §4.2, the `±` structure §4.3,
the MFT constant §4.6) and serves `FLAT_UNSTRIPPED.md` §0 item 5 (the moments `C⁻₃, C⁻₅, C⁻₇` needed by `Ω^U_{g3}`).

Marks: **[D]** derived; **[DVE]** derived and verified end to end on exact blocks; **[NE]** numerical evidence; **[I]** inferred; **[ND]** not done.

## 0. Summary

SUMMARY_PLACEHOLDER

## 1. The prediction at odd `k` [D]

**Moments.** `Π̂^q_{k,ν} = ∫[dη] 𝒦^ν_{q_1+q_3,q_2+q_4}(ξ) ∫dr r^{k−2} 𝒢_q(r,η)` (`SPINNING_4PT.md` §1) is defined for every real `k` for which it
converges. In the code `k` enters only through the weight `(ρρ̄)^{(k−4)/2}` of the grid sum. No new block was computed: the odd-`k`
moments are new integrals of the exact cached waves of `SPINNING_4PT.md` §5 (the grid points dropped there, §2, are partly recomputed).

**Factorised formula.** The derivation of `SPINNING_4PT.md` §4.2 uses `k` only through the two one-pair Laplace transforms. These give
`γ_a(ν)²m^{−2a}` and the tables `M = R + g_1S` evaluated at `a = k + 2l + 1 = k + 3`, with `g_1 = γ_{a−1}/γ_a`. The tables are rational
in `a` (`e2e/laplace.py`), and no parity of `k` is used anywhere. So the prediction at odd `k` is the same formula evaluated at
`a = 6, 8, 10`:

  `(−1)^J Π̂^q_{k,ν}[s_13²s_24²W] → (π/4)4⁴/b^{(0)} · γ_{k+3}(ν)² m^{−2(k+3)} · Σ_h V V M^{q_1q_3}_{h_1h_3}(ν;k+3) M^{q_2q_4}_{h_2h_4}(−ν;k+3) 𝒦^ν(β)`,

with the `l = 1` vertices `V` of §4.1, which do not depend on `k`. Nothing is fitted.

**The `±` structure.** In the preserved–preserved component the prediction is `P[R_LR_R(Σ_s+Σ_o) − g_1²S_LS_R(Σ_s−Σ_o)]` at any `k`
(§4.3). The even kernel (`D_1⊗D_1`) picks `Σ_+ = Σ_s + Σ_o`, and the odd kernel (`D_5⊗D_5`) picks `Σ_− = Σ_s − Σ_o`, with the physical sign
of §4.5. The two normalisations are regular and smooth across `k` (pair `(1,−1,1,−1)`):

| `k` | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|---|
| `R(ν=0)` | 0.3650 | 0.3291 | 0.3102 | 0.2986 | 0.2907 | 0.2850 |
| `S(ν=0)` | 0.500 | 0.625 | 0.750 | 0.875 | 1.000 | 1.125 |
| `g_1(ν=0)` | 0.6580 | 0.4962 | 0.3981 | 0.3322 | 0.2850 | 0.2495 |
| `g_1²S²/R²` (`ν = 0`; `ν = 1`) | 0.812; 0.890 | 0.888; 0.922 | 0.926; 0.944 | 0.948; 0.958 | 0.961; 0.968 | 0.970; 0.974 |

(`k = 6`, `ν = 0`: `0.961` is the `K_odd/K_even` of `SPINNING_4PT.md` §4.5.)

**`C^±` and the MFT constant.** Normalising the even and odd kernels by `R_LR_R/4` and `g_1²S_LS_R/4` gives
`C^±_{k,ν} = (2P_{−½+iν}(cosh β)/m^{2k})(|c_same|² ± |c_opp|²)`, with the `1/4` per quadruple of §4.6. In that derivation (`tttt/t2_mft_heavy_l1.py`)
the `k`-dependence cancels exactly: `γ_a²` against the functional normalisation, and `m^{−2a}m^{2k} = m^{−6}`. So it holds verbatim at odd `k`
[D, given the factorised formula at that `k`]. **What the AdS moments test is therefore the factorised formula at odd `k`.** The `C⁻`
used by `Ω^U` at odd `k` is its odd-kernel (`Σ_−`) part.

**`ε_odd = ½` is not an AdS-side quantity.** The AdS normalisation is `2P/m^{2k}` for both parities of `k`. The factor `ε_odd = ½` arises only in
the dictionary to the flat sum rules `B′_k` (`FLAT_UNSTRIPPED.md` §0 item 5). There the flat odd density `Z^{(k+1)/2}(…)` and the even one
`(2w−q)Z^{(k+2)/2}(…)` differ by `2w/(…)`, a factor `½` at leading order. The AdS moments cannot test it. It remains [I], tied to the
spinning bulk-point limit.

## 2. Data and the dropped tiny-`t` region (`corner_est.py`, `complete_n2.py`)

* **Waves.** The exact `⟨JJJJ⟩` waves cached by `fourpt/t4_compute_l1.py` were used: `β ≈ 0.4` (`m = 12.4, 22.4, 42.2, 72.0, 101.8`) and `β ≈ 1`
  (`m = 10.6, 20.1, 39.3, 70.0, 100.7`), the `e2e` grid (72 × 62 nodes in `(u,t)`), radial truncation `N = max(4000, 400Δ) ≤ 10⁴`.
* **Combinations.** The waves `W[ĥ_s,ĥ_s]`, `W[ĥ_o,ĥ_o]`, `W[ĥ_o,ĥ_s]` (exact helicity structures) and `W[v_2,v_2]`, which carries both
  `Σ_+` and `Σ_−`.
* **Fits.** Per `(β,m,k,ν)`, the 8-class fit of `fourpt/t5_analysis_l1.py` over 81 modes, unchanged.
* **The dropped region.** The production runs dropped the grid points where the radial series had not converged. These are tiny `t`
  (`r → 1`), 9–25 % of the points. The density there falls steeply towards `t → 0`: it is `O(10⁻⁶)` of the peak at the first dropped column for
  `k = 3`, `m ≈ 72`. Its weight grows as `k` decreases and as `m` grows, because `N` is fixed.
  - Estimate: a power law in `t` from the two smallest accepted nodes, row by row. The bound is the density held at the first accepted node,
    which is monotone. The systematic quoted below is `|R(estimate) − R(bound)|`.
  - Estimated share of the leading mode: `k = 3`: `0.04 %` (m ≈ 12) to `1.2 %` (m ≈ 72); `k = 4`: `≤ 0.14 %`; `k ≥ 5`: `≤ 10⁻⁴`.
* **Why the `lowk` corner method was not used.** `ads4/code/lowk` needs `c_1 + c_2 > 0` for every seed block, and `⟨JJJJ⟩` has seeds with
  `c_1 + c_2 = 0` (equal external dimensions; a log-type crossed-channel singularity). It does not apply.
* **Direct recomputation (`complete_n2.py`).** The dropped points that converge with `N_2 = 30000` radial levels are recomputed with the
  unchanged `fourpt` machinery (seed blocks built one at a time, 1.9 GB). This was done for COMPLETED_CASES; see §3.

## 3. Results

RESULTS_PLACEHOLDER

## 4. Marks and what is not done

MARKS_PLACEHOLDER

## 5. Code (`ads4/code/fourpt/oddk/`)

| file | content | output |
|---|---|---|
| `corner_est.py` | moment with the dropped tiny-`t` region estimated (power law) and bounded (monotone) | — |
| `complete_n2.py` | recomputation of dropped points with `N_2 = 30000` (fourpt machinery, lazy seed blocks) | `cache/modes4c_*.npz`, `cache/complete_*.log` |
| `t1_oddk_l1.py` | odd-`k` (and `k = 4, 6`) end-to-end comparison, `ν = 0, 1`, both `β`; fits in `m` | `out_t1_oddk_l1.txt` |
