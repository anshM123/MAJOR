# Odd-`k` tensor Regge moments with both pairs spinning (`⟨JJJJ⟩`): `C⁻₃, C⁻₅, C⁻₇` end to end

Status: written 18 Nov 2026 by a derivation/numerics agent; **not refereed**; an earlier draft with placeholders was committed, this completed version is not. Code: `ads4/code/fourpt/oddk/` (caches in
`oddk/cache/`, git-ignored). Builds on `SPINNING_4PT.md` (the `⟨JJJJ⟩` partial waves, the factorised formula §4.2, the `±` structure §4.3,
the MFT constant §4.6) and serves `FLAT_UNSTRIPPED.md` §0 item 5 (the moments `C⁻₃, C⁻₅, C⁻₇` needed by `Ω^U_{g3}`).

Marks: **[D]** derived; **[DVE]** derived and verified end to end on exact blocks; **[NE]** numerical evidence; **[I]** inferred; **[ND]** not done.

## 0. Summary

1. **Prediction [D].** The factorised two-pair formula of `SPINNING_4PT.md` §4.2 depends on `k` only through `a = k + 3` in the one-pair
   Laplace tables `M = R + g_1S` and in `γ_a(ν)`. These are rational or meromorphic in `a`, and the derivation uses no parity of `k`. So
   the odd-`k` prediction is the same formula, continued, with nothing fitted.
   * The even and odd kernels pick `Σ_± = Σ_s ± Σ_o`.
   * Normalised by `R_LR_R/4` and `g_1²S_LS_R/4`, they give `C^±_k = (2P/m^{2k})(|c_same|² ± |c_opp|²)` at every `k`. This includes the MFT
     constant `1/4`, whose derivation is `k`-independent.
   * `ε_odd = ½` lives only in the flat dictionary and is not testable by AdS moments [I]. §1.
2. **End to end, `k = 5, 7` [DVE].** Exact `⟨JJJJ⟩` moments over `β ≈ 0.4, 1`, `ν = 0, 1`, `m = 10 … 101`, both `±`:
   * The ratio to the prediction approaches 1 as `1 + c₂/m²`. At the two largest masses `m²(R−1)` is constant: same–same `k = 7`:
     `−102.2, −103.2` (`β ≈ 0.4`) and `−110.1, −111.4` (`β ≈ 1`).
   * `Σ_−` agrees with `Σ_+` to `≤ 2·10⁻³` at `m ≈ 100`. The odd-kernel `C⁻` part is reproduced with sign and normalisation.
   * `O(1/m)` fits give `|c₁| ≤ 0.27` (`k = 5`; `≤ 1.0` where the `m ≈ 100` point is corner-limited) and `≤ 0.04` (`k = 7`). The mutations
     of `SPINNING_4PT.md` §5 gave `|c₁| ≈ 4–6`.
3. **`k = 3` [DVE to `m ≈ 40`; not resolved beyond].**
   * At `m ≤ 42` the ratios are `0.993–0.996`. At `β ≈ 1`, `m = 39.3` the tiny-`t` region was recomputed with `N₂ = 30000`: `Σ_+ = 0.9916`,
     `Σ_− = 0.9953`, same–same `0.9931`, residual systematic `≤ 7·10⁻⁴`.
   * At `m ≥ 70` the tiny-`t` region dropped by the production runs (fixed `N ≤ 10⁴`) carries up to 1–4 % at `k = 3`. Only an estimate
     exists, and it is biased high: at `m = 39` it gave `0.9957`, against `0.9931` once recomputed. So `k = 3` at `m ≥ 70` and its `1/m` fit
     are not resolved.
   * The `N₂` recomputations at `m = 70, 72` were lost in a container restart (about 3 h each) and were dropped.
4. **Odd `k` by continuation in `k` holds for currents [DVE: `k = 5, 7` to `m ≈ 100`; `k = 3` to `m ≈ 40`].** The odd-`k` deviations
   interpolate the even-`k` ones smoothly. Same–same `m²(R−1)` at `β ≈ 1`, `m = 39.3`: `k = 3…7` gives `−10.7, −24.3, −43.9, −70.0, −103.2`,
   monotone and roughly quadratic in `k`. No parity effect is visible at the `10⁻³` level.
5. **Side finding.** The same tiny-`t` region affects the published `k = 4` ratios of `SPINNING_4PT.md` at `m ≳ 70` by `0.3–1 %`; for
   example, same–same at `β ≈ 0.4`, `m = 101.8` is `0.9914` raw and about `1.00 ± 0.01` estimated. `k = 6` is affected at `≤ 10⁻³`.

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
  unchanged `fourpt` machinery (seed blocks built one at a time, 1.9 GB). This was done for `β ≈ 1`, `m = 39.3`: 277 of 606 dropped points, 7700 s on one core. The remaining correction is then `≤ 7·10⁻⁴` at `k = 3`. Against the recomputed value, the power-law estimate on the uncompleted data was `+0.0026` too high, well inside its quoted systematic of `0.0087`. The `m = 70, 72` recomputations (about 3 h each) were killed by a container restart and not repeated.

## 3. Results

Ratios measured/predicted (`R` of the leading class in the 8-class fit; `Σ_±` from `W[v_2,v_2]`). Brackets give the tiny-`t` systematic
in units of `10⁻⁴` (`|R(estimate) − R(bound)|`). Full tables, including `ν = 1`, the raw values, the `W[ĥ_o,ĥ_s]` flip classes and `k = 4, 6`, are in
`out_t1_oddk_l1.txt`.

**`β ≈ 0.4`, `ν = 0`:**

| `m` | `k = 3`: ss / oo / `Σ_+` / `Σ_−` | `k = 5`: ss / oo / `Σ_+` / `Σ_−` | `k = 7`: ss / oo / `Σ_+` / `Σ_−` |
|---|---|---|---|
| 12.4 | 0.9455(12) / 0.359 / 0.9469 / 0.9477 | 0.8021 / 0.042 / 0.8032 / 0.8037 | 0.5959 / −0.264 / 0.5965 / 0.5968 |
| 22.4 | 0.9790(13) / 0.747 / 0.9794 / 0.9798 | 0.9228 / 0.536 / 0.9232 / 0.9235 | 0.8279 / 0.220 / 0.8282 / 0.8286 |
| 42.2 | 0.9961(56) / 0.924 / 0.9962 / 0.9963 | 0.9763 / 0.847 / 0.9764 / 0.9765 | 0.9448 / 0.715 / 0.9449 / 0.9451 |
| 72.0 | not resolved (`±0.02`) | 0.9918(6) / 0.945 / 0.9918 / 0.9919 | 0.9803 / 0.894 / 0.9803 / 0.9804 |
| 101.8 | not resolved (`±0.04`) | 0.9964(21) / 0.973 / 0.9964 / 0.9964 | 0.9900 / 0.946 / 0.9901 / 0.9901 |

**`β ≈ 1`, `ν = 0`** (`m = 39.3` recomputed with `N₂`):

| `m` | `k = 3`: ss / oo / `Σ_+` / `Σ_−` | `k = 5`: ss / oo / `Σ_+` / `Σ_−` | `k = 7`: ss / oo / `Σ_+` / `Σ_−` |
|---|---|---|---|
| 10.6 | 0.9711(17) / 0.437 / 0.9429 / 0.9819 | 0.7901 / −0.056 / 0.7629 / 0.8067 | 0.5399 / −0.403 / 0.5227 / 0.5620 |
| 20.1 | 0.9794(24) / 0.794 / 0.9726 / 0.9857 | 0.9094 / 0.515 / 0.9034 / 0.9226 | 0.7946 / 0.135 / 0.7920 / 0.8171 |
| 39.3 | **0.9931(7) / 0.940 / 0.9916 / 0.9953** | 0.9716 / 0.845 / 0.9705 / 0.9765 | 0.9333 / 0.682 / 0.9336 / 0.9427 |
| 70.0 | not resolved (`±0.03`) | 0.9910(14) / 0.949 / 0.9906 / 0.9926 | 0.9776 / 0.890 / 0.9778 / 0.9809 |
| 100.7 | not resolved (`±0.05`) | 0.9980(51) / 0.977 / 0.9968 / 0.9979 | 0.9890 / 0.945 / 0.9892 / 0.9907 |

**Convergence in `m`** (`m²(R−1)` at the two largest masses; `c₁` from the fit `c₁/m + c₂/m² + c₄/m⁴` over the three largest masses):

| `β`, `ν` | `k = 5` ss: `m²(R−1)`; `c₁` | `k = 5` `Σ_−`: `m²(R−1)`; `c₁` | `k = 7` ss: `m²(R−1)`; `c₁` | `k = 7` `Σ_−`: `m²(R−1)`; `c₁` |
|---|---|---|---|---|
| 0.4, 0 | −42.7, −37.4; +0.25 | −42.2, −37.0; +0.25 | −102.2, −103.2; +0.000 | −101.6, −102.6; +0.001 |
| 0.4, 1 | −38.7, −37.5; +0.07 | −38.2, −36.9; +0.07 | −94.7, −95.6; −0.002 | −94.1, −95.0; −0.002 |
| 1, 0 | −44.1, −20.8; +1.0 (the `m ≈ 100` point is corner-limited, `±0.005`) | −36.5, −21.7; +0.64 | −110.1, −111.4; +0.014 | −93.6, −94.6; +0.007 |
| 1, 1 | −37.2, −31.3; +0.27 | −33.2, −30.5; +0.13 | −99.3, −100.8; −0.003 | −87.9, −89.1; −0.004 |

`k = 3`: the fitted `c₁ ≈ +4` comes from the unresolved `m ≥ 70` points and is not a measurement. The resolved points (`m ≤ 42`) give
`m²(R−1) = −7 ± 10` (`β ≈ 0.4`, estimate) and `−10.7 ± 1` (`β ≈ 1`, recomputed).

**Smoothness in `k`.** Same–same `m²(R−1)` at the largest resolved mass of each `k`:

| | `k = 3` | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|
| `β ≈ 1`, `m = 39.3` (recomputed) | −10.7 | −24.3 | −43.9 | −70.0 | −103.2 |
| `β ≈ 0.4`, `m = 42.2` | −7 ± 10 | −23.7 | −42.2 | −66.8 | −98.3 |

The odd `k` interleave the even `k` monotonically. The opposite-helicity class and the flip classes behave the same way, with larger
`c₂` (for example `−557` at `k = 7`) and no `O(1/m)` term.

## 4. Marks and what is not done

* **[D]** The odd-`k` prediction is the factorised formula continued in `k` (§1). `C^±_k` at odd `k` follows with the same normalisations and
  the `k`-independent MFT constant.
* **[DVE]** For currents, the factorised formula at `k = 5, 7` holds for both `β`, both `ν`, both `±`, to `m ≈ 100`:
  - deviations `≤ 1 %` at `m ≈ 100`, behaving as `c₂/m²`;
  - no `O(1/m)` term at the level `|c₁| ≲ 0.3`;
  - the odd-kernel (`C⁻`) part `Σ_−` reproduced with its sign and normalisation.
* **[DVE]** `k = 3` holds to `m ≈ 40`: `0.993–0.996`, `Σ_−` included.
* **"Odd `k` by continuation in `k`": holds for currents [DVE in the range above].**
* **[ND]**
  - `k = 3` at `m ≳ 70` and its `1/m` fit. This needs the recomputation of the tiny-`t` region (`N₂`; about 3 h per case on one core), or a
    corner method for seeds with `c₁ + c₂ = 0`.
  - `⟨TTTT⟩` at odd `k`. The Taylor-mode caches of `SPINNING_TTTT.md` hold only `m = 39, 70` at `β = 1`, and the moment weight could be
    changed in the same way. Not done.
  - `ν > 1`.
* **[I]**
  - `ε_odd = ½` and the flat-limit (bulk-point) side of the dictionary `C_k → ε_k^{−1}B′_k`.
  - Positivity or validity of `Ω^U_{g3}` in AdS; nothing here addresses it.
  - The light sides.

## 5. Code (`ads4/code/fourpt/oddk/`)

| file | content | output |
|---|---|---|
| `corner_est.py` | moment with the dropped tiny-`t` region estimated (power law) and bounded (monotone) | — |
| `complete_n2.py` | recomputation of dropped points with `N_2 = 30000` (fourpt machinery, lazy seed blocks) | `cache/modes4c_*.npz`, `cache/complete_*.log` |
| `t1_oddk_l1.py` | odd-`k` (and `k = 4, 6`) end-to-end comparison, `ν = 0, 1`, both `β`; fits in `m` | `out_t1_oddk_l1.txt` |
