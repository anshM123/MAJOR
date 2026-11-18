# Odd-`k` tensor Regge moments with both pairs spinning (`⟨JJJJ⟩`): `C⁻₃, C⁻₅, C⁻₇` end to end

Status: written 18 Nov 2026 by a derivation/numerics agent; refereed (`REFEREE_SPINNING_ODDK.md`), revised 18 Nov (§6); not committed. Code: `ads4/code/fourpt/oddk/` (caches in
`oddk/cache/`, git-ignored). Builds on `SPINNING_4PT.md` (the `⟨JJJJ⟩` partial waves, the factorised formula §4.2, the `±` structure §4.3,
the MFT constant §4.6) and serves `FLAT_UNSTRIPPED.md` §0 item 5 (the moments `C⁻₃, C⁻₅, C⁻₇` needed by `Ω^U_{g3}`).

Marks: **[D]** derived; **[DVE]** derived and verified end to end on exact blocks; **[NE]** numerical evidence; **[I]** inferred; **[ND]** not done.

## 0. Summary

(Revised 18 Nov 2026 after `REFEREE_SPINNING_ODDK.md`, §6. The primary numbers now come from a robust extraction: the leading mode, and a
uniform-weight class fit. The `1/|pred|`-weighted 8-class fit inherited from `SPINNING_4PT.md` amplifies the dropped tiny-`t` region by
20–50×; its numbers are kept in `out_t1_oddk_l1.txt` for comparison.)

1. **Prediction [D].** The two-pair factorised formula of `SPINNING_4PT.md` §4.2 depends on `k` only through `a = k + 3` in the one-pair
   Laplace tables `M = R + g_1S` and in `γ_a(ν)`. The derivation uses no parity of `k`, so the odd-`k` prediction is the same formula,
   continued, with nothing fitted.
   * The even and odd kernels pick `Σ_± = Σ_s ± Σ_o`.
   * Normalised by `R_LR_R/4` and `g_1²S_LS_R/4`, they give `C^±_k = (2P/m^{2k})(|c_same|² ± |c_opp|²)` at every `k`. This includes the MFT
     constant `1/4`, whose derivation is `k`-independent.
   * `ε_odd = ½` is a flat-dictionary factor. The heavy moments here (`m ≤ 100`, `ν ≤ 1`) cannot test it; a bulk-point computation at
     `ν ~ m` could. It stays [I]. §1.
2. **End to end, `k = 5, 7` [DVE], currents only.** Exact `⟨JJJJ⟩` moments over `β ≈ 0.4, 1`, `ν = 0, 1`, `m = 10 … 101`, using the leading
   mode:
   * The ratio to the prediction is `1 + c₂/m²`, with `m²(R−1)` flat over the three largest masses. Same helicity: `k = 5`: `−46.6, −47.5,
     −47.7` (`β ≈ 0.4`) and `−66.3, −68.8, −68.9` (`β ≈ 1`); `k = 7`: `−107 … −112` and `−148 … −160`.
   * `|c₁| ≤ 0.02` for every class, `β`, `ν` and parity of `k`, against `c₁ ≈ 4.3` for the referee's `b^{(0)} → b^{(2)}` mutation at `k = 5`.
   * At `m ≈ 100` the same-helicity class is within `0.5 %` (`k = 5`) and `1.6 %` (`k = 7`). The opposite-helicity class is `1.1 %` and
     `2.0 %` off, also converging like `1/m²` (§3).
   * The odd-kernel (`C⁻`) part is tested, with its sign and size. The evidence is the separate convergence of `R_ss` and `R_oo` in the exact
     helicity waves. Under the referee's mutation `S_R → −S_R` (sign of the opposite-sign term), `R_ss` drops to `0.27, 0.10, 0.05`. In
     `W[v_2,v_2]` at `β ≈ 0.4`, `Σ_−` and `Σ_+` coincide trivially (`T_oo/T_ss ≈ 10⁻³`), so they are not independent evidence there.
3. **`k = 3` [DVE to `m ≈ 70`; `m ≈ 100` open at about 0.6 %].**
   * At `β ≈ 1`, `m = 39.3` the dropped tiny-`t` region was recomputed (`N₂ = 30000`): leading mode `0.9872`.
   * Leading mode at `m ≈ 70`: `0.9964` (`β ≈ 1`) and `0.9978` (`β ≈ 0.4`). The tiny-`t` correction there is `+0.19 %` and `+0.12 %`; the
     estimate is biased high, so these are resolved to `≲ 0.3 %`.
   * At `m ≈ 100` the correction is `+0.6 %` (raw `0.9937`, corrected `0.9998`), so that point and the `k = 3` `1/m` fit remain open. At
     `m ≤ 70`, `m²(R−1)` is about `−18` to `−20` (`β ≈ 1`) and `−12` (`β ≈ 0.4`).
4. **Odd `k` by continuation in `k` holds for currents** [DVE: `k = 5, 7`, `m ≤ 100`, `ν ≤ 1`; `k = 3`, `m ≤ 70`].
   * The odd-`k` deviations interpolate the even-`k` ones smoothly. Leading-mode same-helicity `m²(R−1)` at `β ≈ 1`, `m = 39.3`, `k = 3…7`:
     `−19.8, −39.0, −66.4, −102.4, −148.0`. The weighted fit gives `−10.7, −24.3, −43.9, −70.0, −103.2`.
   * Both sets are monotone and smooth, with no parity effect. The `c₂` values themselves depend on the extraction.
5. **Applicability (narrow).** What is tested is the `l = 1` analogue of the project's `C⁻_k`: the `D₅⊗D₅` part of the two-pair Regge moment
   at odd `k`, for currents, at `ν ≤ 1` and `m ≳ 40`. **Not tested:**
   * the graviton `C⁻_k` of `Ω^U` (`𝔄⊗𝔄` kernel, `⟨TTTT⟩`, `a = k + 5`, `ν = Λp` up to `ν ~ m`, `m` near the gap where the odd-`k` `c₂`
     are large);
   * the `Ψ_k` of `Ω^SC` (meromorphic weights, antipodal kernels);
   * whether the odd-`k` moment is the heavy action of a valid dispersive functional (`C⁻₃` is not in the CLS-D family).

   **H3 stays [I] for gravitons.** For currents it is supported at `k = 5, 7` (and `k = 3` to `m ≈ 70`), `ν ≤ 1`, `m ≳ 40`.
6. **Correction to `SPINNING_4PT.md`** (dated note added there). Its `k = 4` weighted class-fit ratios at `m ≳ 70` carry an estimated
   dropped-region shift of `+0.4 %` to `+2.1 %`. Leading-mode `k = 4` moves by `≤ 1.5·10⁻³` and reads `0.996–0.997` at `m ≈ 100`. `k = 6`
   is unaffected.

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
`(2w−q)Z^{(k+2)/2}(…)` differ by `2w/(…)`, a factor `½` at leading order. The heavy moments computed here (`ν ≤ 1`) cannot test it; a
bulk-point computation at `ν ~ m` could. It remains [I], tied to the
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
  - Estimate: a power law in `t` from the two smallest accepted nodes, row by row. A second, larger estimate holds the density at the first
    accepted node; it is called the "bound" in the code, but it is not a bound. It assumes monotone decrease in `t`, and for skipped
    `u`-rows it is a crude extrapolation. The spread between the two is a one-sided systematic estimate.
  - **Extraction matters.** The `1/|pred|`-weighted class fit gives near-zero modes the same weight as the leading mode, and those modes
    carry a larger tiny-`t` share. Per mode, the leading-mode correction at `k = 5`, `m ≈ 100` is `10⁻⁴`, while the weighted fit moves by
    `5·10⁻³`. Leading-mode and uniform-weight extractions are therefore adopted.
  - Estimated share of the leading mode: `k = 3`: `0.04 %` (m ≈ 12) to `1.2 %` (m ≈ 72); `k = 4`: `≤ 0.14 %`; `k ≥ 5`: `≤ 10⁻⁴`.
* **Why the `lowk` corner method was not used.** `ads4/code/lowk` needs `c_1 + c_2 > 0` for every seed block, and `⟨JJJJ⟩` has seeds with
  `c_1 + c_2 = 0` (equal external dimensions; a log-type crossed-channel singularity). It does not apply.
* **Direct recomputation (`complete_n2.py`).** The dropped points that converge with `N_2 = 30000` radial levels are recomputed with the
  unchanged `fourpt` machinery (seed blocks built one at a time, 1.9 GB). This was done for `β ≈ 1`, `m = 39.3`: 277 of 606 dropped points, 7700 s on one core. The remaining correction is then `≤ 7·10⁻⁴` at `k = 3`. Against the recomputed value, the power-law estimate on the uncompleted data was `+0.0026` too high, well inside its quoted systematic of `0.0087`. The `m = 70, 72` recomputations (about 3 h each) were killed by a container restart and not repeated.

## 3. Results (robust extraction; `out_t1_oddk_l1.txt`, sections "ROBUST")

Leading mode `R = meas/pred` for `W[ĥ_s,ĥ_s]` (ss) and `W[ĥ_o,ĥ_o]` (oo), corner-corrected, with raw values in brackets, and the
uniform-weight `Σ_±` of `W[v_2,v_2]`; `ν = 0`. `ν = 1`, `k = 4, 6`, the uniform-weight class ratios and the weighted-fit numbers are in the
output file.

**`β ≈ 0.4`:**

| `m` | `k = 3`: ss / oo / `Σ_+` / `Σ_−` | `k = 5`: ss / oo | `k = 7`: ss / oo |
|---|---|---|---|
| 12.4 | 0.9260 / 0.503 / 0.932 / 0.929 | 0.7657 / 0.400 | 0.5474 / 0.277 |
| 22.4 | 0.9743 / 0.822 / 0.975 / 0.975 | 0.9130 / 0.759 | 0.8106 / 0.665 |
| 42.2 | 0.9926 [0.9923] / 0.9475 / 0.993 / 0.993 | 0.9738 / 0.9256 | 0.9400 / 0.8896 |
| 72.0 | 0.9978 [0.9966] / 0.9821 / 0.999 / 0.998 | 0.9908 / 0.9737 | 0.9787 / 0.9603 |
| 101.8 | 0.9997 [0.9965] / 0.9919 — open (`~0.6 %`) | 0.9954 / 0.9868 | 0.9892 / 0.9799 |

**`β ≈ 1`** (`m = 39.3` recomputed with `N₂`):

| `m` | `k = 3`: ss / oo / `Σ_+` / `Σ_−` | `k = 5`: ss / oo | `k = 7`: ss / oo |
|---|---|---|---|
| 10.6 | 0.8705 / 0.677 / 0.873 / 0.875 | 0.6295 / 0.468 | 0.3665 / 0.263 |
| 20.1 | 0.9552 / 0.889 / 0.953 / 0.957 | 0.8568 / 0.784 | 0.7068 / 0.637 |
| 39.3 | **0.9872 / 0.9685 / 0.986 / 0.988** | 0.9571 / 0.9344 | 0.9044 / 0.8787 |
| 70.0 | 0.9964 [0.9945] / 0.9903 / 0.997 / 0.996 | 0.9860 / 0.9784 | 0.9679 / 0.9590 |
| 100.7 | 0.9998 [0.9937] / 0.9969 — open (`~0.6 %`) | 0.9932 / 0.9895 | 0.9842 / 0.9798 |

At `k = 5, 7` the corner correction of the leading mode is `≤ 10⁻⁴`; raw and corrected values agree to the digits shown.

**Convergence** (`m²(R−1)` at the three largest masses; `c₁` from the fit `c₁/m + c₂/m² + c₄/m⁴`):

| `β`, `ν` | `k = 5` ss | `k = 5` oo | `k = 7` ss | `k = 7` oo |
|---|---|---|---|---|
| 0.4, 0 | −46.6, −47.5, −47.7; `c₁ = +0.004` | −132.7, −136.2, −137.1; `+0.003` | −106.9, −110.6, −111.6; `−0.002` | −196.7, −205.7, −208.2; `−0.006` |
| 0.4, 1 | −42.0, −42.7, −42.9; `−0.002` | −124.0, −127.0, −127.8; `−0.003` | −100.3, −103.6, −104.4; `−0.002` | −186.6, −194.5, −196.7; `−0.005` |
| 1, 0 | −66.3, −68.8, −68.9; `+0.019` | −101.6, −105.8, −106.4; `+0.016` | −147.9, −157.4, −159.9; `−0.009` | −187.7, −201.2, −204.8; `−0.013` |
| 1, 1 | −50.4, −51.8, −52.4; `−0.011` | −80.4, −82.9, −83.7; `−0.011` | −123.6, −130.0, −131.7; `−0.005` | −159.0, −168.3, −170.7; `−0.007` |

The uniform-fit `Σ_±` behave the same way (`|c₁| ≤ 0.06` at `k = 5, 7`). `k = 4, 6` for comparison: `|c₁| ≤ 0.12` (`k = 4`) and `≤ 0.005` (`k = 6`).

**`k = 3`:** `m²(R−1)` (ss, leading mode) is `−19.8, −17.9` at `m = 39, 70` (`β ≈ 1`) and `−13.2, −11.5` at `m = 42, 72` (`β ≈ 0.4`). The
`m ≈ 100` values (`−2.2`, `−2.7`) carry the open `0.6 %` corner correction. The fitted `c₁` (`0.4–0.65` at `ν = 0`, `≤ 0.3` at `ν = 1`) is
dominated by that point and is not a measurement.

**Smoothness in `k`** (same helicity, `m²(R−1)`, `β ≈ 1`, `m = 39.3`, recomputed):

| extraction | `k = 3` | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|
| leading mode | −19.8 | −39.0 | −66.4 | −102.4 | −148.0 |
| weighted 8-class fit | −10.7 | −24.3 | −43.9 | −70.0 | −103.2 |

The odd `k` interleave the even `k` monotonically in both. The `c₂` values themselves are extraction-specific, because the extractions
differ at `O(c₂/m²)`.

## 4. Marks and what is not done

* **[D]** The odd-`k` prediction is the factorised formula continued in `k` (§1). `C^±_k` at odd `k` follows with the same normalisations and
  the `k`-independent MFT constant.
* **[DVE], currents, `ν ≤ 1`.** The factorised formula holds at `k = 5, 7` to `m ≈ 100` and at `k = 3` to `m ≈ 70`, for both `β` and both
  helicity classes (and so both `±`). The approach is `1 + c₂/m²` with `|c₁| ≤ 0.02` (`k = 5, 7`). The referee's mutations break it:
  - `b^{(2)}`: `c₁ ≈ 4.3`;
  - `k = 4/6`-mean tables: `R ≈ 0.92`;
  - `C⁻` sign flip: `R ≈ 0.05`.
* **"Odd `k` by continuation in `k`": holds for currents** in this range.
* **Not resolved:** `k = 3` at `m ≈ 100` (about `0.6 %` corner correction) and its `1/m` fit.
* **[ND]**
  - `⟨TTTT⟩` at odd `k`.
  - The graviton `C⁻_k` of `Ω^U`, i.e. the `𝔄⊗𝔄` kernel and `ν` up to `~m`.
  - `Ψ_k` of `Ω^SC`.
  - `ν > 1`.
  - `m` near the gap.
* **[I]**
  - `ε_odd = ½` and the bulk-point side of `C_k → ε_k^{−1}B′_k`.
  - H3 for gravitons.
  - Validity and positivity of `Ω^U` in AdS.
  - The light sides.

## 5. Code (`ads4/code/fourpt/oddk/`)

| file | content | output |
|---|---|---|
| `corner_est.py` | moment with the dropped tiny-`t` region estimated (power law; a second, larger estimate, not a bound) | — |
| `complete_n2.py` | recomputation of dropped points with `N_2 = 30000` (fourpt machinery, lazy seed blocks) | `cache/modes4c_*.npz`, `cache/complete_*.log` |
| `t1_oddk_l1.py` | odd-`k` (and `k = 4, 6`) end-to-end comparison, `ν = 0, 1`, both `β`: weighted 8-class fit, uniform-weight fit and leading mode; fits in `m` | `out_t1_oddk_l1.txt` |

## 6. Revisions after `REFEREE_SPINNING_ODDK.md` (18 Nov 2026)

* **M1.** A robust extraction (leading mode and uniform-weight class fit) is adopted and reported next to the weighted fit. The verdicts are
  updated:
  - `k = 5, 7`: `|c₁| ≤ 0.02`, `m²(R−1)` flat;
  - `k = 3`: resolved to `≲ 0.3 %` at `m ≈ 70`; only `m ≈ 100` is open.
* **M2.** The `k = 4` shift in `SPINNING_4PT` is `+0.4` to `+2.1 %` (weighted fit), not `0.3–1 %`. A dated correction note was added to
  `SPINNING_4PT.md` §5.
* **M3.** The applicability statement was added (§0 item 5). H3 stays [I] for gravitons.
* **Minor revisions:**
  - m1: `Σ_− ≈ Σ_+` at `β ≈ 0.4` is no longer cited as evidence.
  - m2: the opposite-helicity class is reported as `1–2 %` off at `m ≈ 100`.
  - m3: the corner "bound" is called an estimate.
  - m4: the wording on `ε_odd` is corrected.
  - m5: both extractions are shown in the smoothness table.
