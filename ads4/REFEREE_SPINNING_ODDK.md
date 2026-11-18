# Referee report: `SPINNING_ODDK.md` (odd-`k` two-pair moments of `⟨JJJJ⟩`)

Referee: independent, hostile. 18 Nov 2026. Nothing under review was edited. Scratch code: `ads4/code/referee_oddk/`
(`ref_k5.py`, `ref_fitdep.py`, `ref_corner.py`; logs in `referee_oddk/cache/`, git-ignored). Each script runs in under a minute on one core
from the existing caches.

## Verdicts

* **Correctness: largely correct.** The numbers reproduce, and the main claim (the factorised formula continued to `k = 5, 7`) holds. It
  holds more strongly than reported once a robust extraction is used (M1). Claims (3) and (5) are artefacts of the estimator (M1, M2).
* **Soundness: weakened by the extraction.** The `1/|pred|`-weighted 8-class fit, inherited from `SPINNING_4PT`, amplifies the dropped tiny-`t`
  region by 20–50×. The "unresolved" `k = 3` points and the `k = 4` "side finding" come mostly from that weighting (M1).
* **Labels: several overstated or mis-cited** (m1–m4).
* **Applicability: narrow.** The tested object is the `l = 1` analogue of the project's `C⁻_k`. It is neither the graviton `C⁻_k` of `Ω^U` nor
  the `Ψ_k` of `Ω^SC`. It is tested only at `ν ≤ 1` and `m ≳ 10` (M3).
* **Novelty: modest.** The paper reuses cached waves with a new radial weight and adds a corner estimate and one recomputation.

## What was checked

1. **Reproduction.** I recomputed, with an independent driver, the `k = 5`, `β ≈ 0.4`, `ν = 0` ratios at `m = 42.2, 72.0, 101.8`: `W[ĥ_s,ĥ_s]`
   `R_ss = 0.9763, 0.9918, 0.9964` and `W[ĥ_o,ĥ_o]` `R_oo = 0.8468, 0.9451, 0.9726`. Both are identical to `out_t1_oddk_l1.txt`. The
   one-pair table of §1 (`R, S, g_1, g_1²S²/R²` for `k = 2…7`) is also reproduced.
2. **The `k`-dependence enters only through `a = k+3`.** I confirmed this in `SPINNING_4PT` §4.2, `predict4.py` and `e2e/laplace.py`:
   * `a = 2l+k+1` enters the prefactor `γ_a² m^{−2a}` and the rational `R, S`.
   * `g_1 = γ_{a−1}/γ_a`.
   * The vertices, `b^{(0)}` and `𝒦` do not depend on `k`.
   * The moment integrates over one sheet only, `(ρρ̄)^{(k−4)/2}` with `ρ, ρ̄ > 0`, so no parity of `k` appears.
3. **Mutations at `k = 5`** (`ref_k5.py`, `β ≈ 0.4`, `ν = 0`; the paper ran none at odd `k`):

   | mutation | `R_ss` at `m = 42, 72, 102` | `c₁(A)` |
   |---|---|---|
   | none | 0.976, 0.992, 0.996 | +0.25 |
   | `b^{(0)} → b^{(2)}` | 1.075, 1.050, 1.037 | +4.3 |
   | tables at `k=5` replaced by the mean of the `k=4` and `k=6` tables (a naive "continuation") | 0.904, 0.918, 0.922 | −14 |
   | `S_R → −S_R` (sign of the `g_1²S_LS_R` term, i.e. the `C⁻` sign) | 0.27, 0.096, 0.048 | — |
   | `S → S/2` | 2.45, 2.48, 2.49 | — |
   | `S → −S` on both sides | unchanged | — |

   * The odd-`k` data discriminate the continuation and the relative sign of the opposite-sign term decisively. At `k = 5` the same-helicity
     weight `R² − g_1²S²` is a 95 % cancellation, so `R_ss` tests `S` at the `10⁻⁴` level.
   * The global sign of `S` is invisible, because `S` enters the class sums only quadratically.

## Issues

**M1 (major). The extraction amplifies the dropped region, and the "unresolved" verdicts are artefacts of it.** The weighted 8-class fit
(`w = 1/|pred|`) gives modes with small predicted values the same weight as the leading mode. Those modes carry a larger tiny-`t` share. Below,
"weighted fit" is the paper's estimator; "leading mode" and "uniform-weight fit" are robust alternatives, run in `ref_fitdep.py`:

| case | weighted fit, raw → corrected | uniform-weight fit, raw → corrected | leading mode, raw → corrected |
|---|---|---|---|
| `k=5`, `β≈1`, `m=100.7`, ss | 0.9929 → 0.9980 | 0.9931 → 0.9933 | 0.9931 → 0.9932 |
| `k=4`, `β≈1`, `m=100.7`, ss | 0.9833 → 1.0047 | 0.9949 → 0.9967 | 0.9955 → 0.9963 |
| `k=3`, `β≈1`, `m=70.0`, ss | 0.9688 → 1.0031 | 0.9927 → 0.9968 | 0.9945 → 0.9964 |
| `k=3`, `β≈0.4`, `m=72.0`, ss | 0.9798 → 1.0060 | 0.9956 → 0.9983 | 0.9966 → 0.9978 |

* For the leading mode, the per-mode corner estimate is `1·10⁻⁴` at `k = 5` and `m ≈ 100`. That matches §2's own "`k ≥ 5`: `≤ 10⁻⁴`", which
  contradicts the `(51)` systematic in the §3 table without explanation.
* With the leading mode the convergence is clean:

  | case | `m²(R−1)`, three largest `m` | `c₁(A)` |
  |---|---|---|
  | `k = 5`, `β ≈ 1` | −66.4, −68.7, −69.0 | +0.008 |
  | `k = 5`, `β ≈ 0.4` | −46.7, −47.7, −47.7 | +0.013 |
  | `k = 7`, `β ≈ 1` | — | −0.027 |
  | `k = 7`, `β ≈ 0.4` | — | −0.026 |
  | `k = 4` | — | +0.05 / +0.10 |

  So claim (2) is conservative: `|c₁| ≲ 0.03` against ≈ 4 for the `b^{(2)}` mutation, instead of "`≤ 1.0` where corner-limited".
* `k = 3` at `m ≈ 70` is resolved to `≲ 0.3 %` in the leading mode (raw 0.9945/0.9966, corrected 0.9964/0.9978).
* `k = 3` at `m ≈ 100` remains corner-sensitive at about 0.6 % (raw 0.9937, corrected 0.9998), so its `1/m` fit is still open.
* The ratios also depend on the extraction at `O(c₂/m²)`. At `m ≈ 40`, `R_oo` is 0.845 (weighted fit) against 0.934 (leading mode) at `k = 5`. The
  quoted `c₂` values (e.g. −557) are therefore estimator-specific.
* **Fix:** report leading-mode or uniform-weight ratios alongside the weighted fit; withdraw "k = 3 at m ≥ 70 not resolved (1–4 %)".

**M2 (major, affects `SPINNING_4PT`). Claim (5) is right in substance but misstated.**
* In the published (weighted-fit) `k = 4` ratios at `m ≳ 70`, the shift is estimated at +0.4 % (`β≈0.4`, `m=72`), +1.0 % (`β≈0.4`, `m=101.8`),
  +0.7 % (`β≈1`, `m=70`) and **+2.1 %** (`β≈1`, `m=100.7`). It is not "0.3–1 %".
* The estimator is biased high: at `m = 39`, `k = 3`, it overshot the recomputed shift by about 40 %.
* `SPINNING_4PT`'s "tail estimate" (`2·10⁻⁵ … 3·10⁻³`, `tmoments4.tail_estimate`) only probes an `X`-cut on the accepted points and is blind to
  the convergence-masked points. The paper's numerical-accuracy sentence and weak spot (iii) of §6.2 therefore understate the error.
* For the leading mode, `k = 4` moves by `≤ 1.5·10⁻³`, and the trend is smooth: `m²(R−1) ≈ −26` (`β ≈ 0.4`) and `−38` (`β ≈ 1`).
* `k = 6`, the headline of `SPINNING_4PT`, is unaffected: `≤ 7·10⁻⁴` in the weighted fit and `0` in the leading mode. The one exception is the
  `β ≈ 1`, `ν = 0` `c₁` of the weighted fit, which moves from −0.15 to +0.13, still far below the mutations.
* **`SPINNING_4PT` needs a correction note:**
  * Its `k = 4` class-fit ratios at `m ≳ 70` carry a dropped-region systematic of up to about 1–2 %.
  * The robust (leading-mode) `k = 4` ratios are 0.997 and 0.996 at `m ≈ 100`.
  * The `k = 6` conclusions stand.

**M3 (major, applicability). The test covers neither the graviton `C⁻_k` of `Ω^U` nor `Ψ_k` of `Ω^SC`.**
* It tests the `l = 1` analogue: the `D₅⊗D₅` part of the s-channel Regge moment at odd `k`. That is the right analogue of the project's
  `C⁻_k`, and the `S_R → −S_R` mutation shows the opposite-sign term is genuinely measured.
* `Ω^U` needs the following, none of which is touched here:
  * `𝔄⊗𝔄` with `a = k+5` for `⟨TTTT⟩` (ND, as stated);
  * `ν = Λp` up to `ν ~ m` (here `ν ≤ 1`);
  * `m` near the gap, where the odd-`k` corrections are large: `c₂` grows roughly like `k²`, to −100 (ss) and −560 (oo) at `k = 7`;
  * the [I] statement that the project's odd-`k` moment is the heavy action of a valid dispersive functional (`LIGHT_INPUTS.md` §0 item 3:
    `C⁻₃` is not in the CLS-D family).
* `Ω^SC`'s `Ψ_k` uses the meromorphic weights `A_k`, antipodal `j′ = 1` kernels at `(1,3)` and `𝔄 = D₅D₁` at `(2,4)`. Its heavy action equals
  `N_k ×` the project form only given P1–P3. This paper supports the right-hand factor, for currents only.
* H3 of `RESULT_D` should read: "supported for currents at `k = 5, 7`, `ν ≤ 1`, `m ≳ 40`; gravitons not tested".

**m1.** "Σ₋ matches Σ₊ to `≤ 2·10⁻³`; the `C⁻` part is reproduced with sign and normalisation" cites the wrong evidence.
* At `β ≈ 0.4`, `T_oo/T_ss ≈ 1.6·10⁻³` in `W[v_2,v_2]`, so `R(Σ₋) ≡ R(Σ₊) ≡ R_ss`. This was already `SPINNING_4PT` referee point m4.
* The evidence is the separate convergence of `R_ss` and `R_oo` in the exact-helicity waves, together with the sign mutation above.

**m2.** "Deviations `≤ 1 %` at `m ≈ 100`" (§0 item 2, §4) is false for the opposite-helicity class.
* `k = 7`: `R_oo = 0.945` (weighted fit) or 0.980 (leading mode).
* `k = 5`: 0.973–0.977 (weighted fit) or 0.987–0.989 (leading mode).
* Only the same-helicity class and `Σ_±` are within 1 %. Convergence like `1/m²` does hold for `oo`.

**m3.** The corner "bound" is not a bound.
* It assumes monotone decrease in `t` "as observed".
* For skipped `u`-rows (248 points at `m ≈ 100`) it is a power-law extrapolation weighted by a crude ratio of maximum weights.
* The quoted systematic `|R(est) − R(bound)|` is one-sided.
* Label it as an estimate.

**m4.** "`ε_odd = ½` untestable from AdS moments" is too strong. It is a flat-dictionary factor, and the heavy moments at `m ~ 100`, `ν ≤ 1`
cannot test it. A bulk-point computation at `ν ~ m` could.

**m5.** The smoothness table (claim 4) is estimator-specific.
* Leading mode at `β ≈ 1`, `m = 39.3`, `k = 3…7`: `−19.8, −39.0, −66.4, −102.4, −148.0`.
* Weighted fit (the paper's numbers): `−10.7, −24.3, −43.9, −70.0, −103.2`.
* Both are monotone and smooth, so the conclusion stands.

## Bottom line

The factorised formula continued in `k` is verified end to end for currents at `k = 5, 7` (`ν ≤ 1`, `m ≤ 100`) with no fitted parameter. The
verification includes the sign and size of the opposite-sign (`C⁻`) term, and is robust to the extraction. `k = 3` is verified to about 0.3 % at
`m ≲ 70` by the leading mode.

The "unresolved"/"side-finding" parts should be rewritten around a robust extractor. `SPINNING_4PT` needs a short correction for its `k = 4`
large-`m` class ratios and its tail estimate. H3 remains [I] for gravitons.
