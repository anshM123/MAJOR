# Referee report: `ads4/SPINNING_LOWK.md` and `ads4/code/lowk/`

Referee: independent, hostile. Date: 16 Nov 2026. Object: commit `706738c`. I did not edit anything under review.
Scratch scripts are in `ads4/code/referee_lowk/`. Their logs and pickles are in `ads4/code/referee_lowk/cache/`, which is git-ignored.

## Verdicts

| | Verdict |
|---|---|
| Correctness | **Reproduced; no wrong number found.** I rebuilt one case from scratch and got the published values to all printed digits. Tail exponents and η-derivative tails pass independent checks. The leading formula holds at k = 2, 3 for `l = 1` (both helicities) and for the `l = 2` preserved component, at the stated accuracy, within the factor-2–3 caveat of issue M1. |
| Soundness | **Sound method; error bars are not bounds** (M1). The quoted errors are the scatter between two neighbouring tail settings. An equally good third setting moves results by 1.3–3× the quoted error. |
| Labels | Mostly honest. Some uncertainties are understated and one check is misread (M2, m1–m4). |
| Applicability | Narrow: see the last section. |
| Novelty | Technical and incremental. Lightcone reorganisation plus Darboux tails is a standard asymptotic technique applied to a known problem. |

No fatal issues.

## What I ran

### 1. Reproduction from scratch (`ref_tail.py`, `ONLY=A,B,MUT`)

Case: `l = 1`, `(Δ,J) = (44.565161…, 8)`, m ≈ 42, β ≈ 0.4. I rebuilt the corner-refined grid and the LC seeds A and B, then recomputed the moments.

| quantity | value | published |
|---|---|---|
| k = 2, ν = 0, `R_pres` | 0.990504910 ± 8.4e-7 | identical |
| k = 2, ν = 0, `R_flip` | 0.980559445 ± 6.3e-5 | identical |
| k = 3, ν = 0, `R_pres` | 0.984598889 ± 2.4e-9 | identical |
| k = 3, ν = 0, `R_flip` | 0.970379447 ± 7.9e-8 | identical |

All of `nfill`, `n_rejected` and `valid` are identical to the published build. The script must be given the exact Δ from `cases`. With Δ rounded to 6 decimals the grid shifts, the validation degrades to 1e-5, and R moves by about 4e-7. That is my own pitfall, not a bug in the note.

### 2. Tail-ansatz and fit-window variations

Same grid and the same acceptance mask throughout. The table shows the shift of R relative to setting A, k = 2, ν = 0, (pres / flip).

| setting | ΔR (pres / flip) | corrected Casimir residual (deep corner) |
|---|---|---|
| B | 3e-9 / 1.7e-7 | 2e-10 |
| fit window [1000, 8000] | 1e-8 / 9e-7 | 1e-11 |
| M = 4 | 1e-8 / 1e-6 | 7e-10 |
| log terms at `z̄ = ∞` | 5e-8 / 5.5e-6 | — |
| Q = 3000 | 1.1e-6 / 1.1e-4 (1.3× / 1.8× the quoted error) | 5e-9, as good as A |
| M = 12 | 3.0e-5 / 3.1e-3 (37× / 50×) | 3e-6 |
| window [6000, 8000] | 1.6e-4 / 1.6e-2 (200×) | 5e-6 |

At k = 3 the shifts are: Q3000 flip 2.3e-7 (3× the quoted 7.8e-8); M12 6e-8 / 5e-6; window [6000, 8000] 3.5e-7 / 2.8e-5 (145× / 360×).

The settings that are badly off are exactly those with a large corrected Casimir residual. So the Casimir residual, not |A − B|, is the usable quality gate (see m1). In absolute terms all of these shifts are small (k = 2: `m²ΔR ≤ 0.05` for the acceptable settings), so no conclusion changes.

### 3. Mutation

I changed the main Darboux exponent from `2σ+l+t−m` to `2σ+l+t−m+0.25`, a wrong Casimir exponent.

* **Result:** k = 2 gives `R_pres = 1.50` and `R_flip = 52.9`; k = 3 shifts by 1.2e-3. The mutation is caught.
* **But not by the error budget:** `err_lc = 3|A−B|` is unchanged (2.6e-8 against 4.0e-8), because A and B share the exponent. The `fitres` diagnostic is also identical to the unmutated value.
* **What actually protects the exponents:** the block-level checks. These are the ₂F₁ comparison and the Casimir residual (`casimir_fixed.py`), plus the η-derivative check below.

### 4. η-derivative tails (`fd_check.py`)

This tests the `l ≥ 1` tails and their exponent `−2σ−l` in the deep corner. The note has no test there: its Casimir check uses `l = 0` only.

* Method: a 6th-order finite difference of `S[k,l]` in χ, compared with `S[k,l+1]`.
* Seed: `(−0.375, 0.625)`; r = 0.3 and 0.6; `1−ρ_b = 1e-3 … 1e-7`.
* Result: agreement to ≤ 3e-7, with every `l → l+1` up to `l = 4`.

### 5. Leave-out test of the missing-point estimate (`leaveout2.py`)

I removed accepted nodes next to the missing region and compared the power-law estimate with the true values.

* One node removed: est/true = 1.00–1.19.
* Two nodes removed: est/true = 0.44–1.49.

So the 50 % error assignment is adequate for one-node extrapolation and marginal further out.

### 6. Outer region of (★) (`outer2.py`)

I evaluated the exact shadow integrand along five rays in X or Y, for m = 44 and 182, `l = 1`, from ζ = 1 up to the causal-diamond boundary `ζ_max = m e^{−d(x̂,e)} cosh d(x̂,p̂)`.

* Inside the support, exact/leading decreases monotonically (down to `10^{−200}` at the boundary).
* The density follows `e^{−2ζ}`.
* Outside the support the integrand would be huge (up to `10^{30}` above the typical value), so the bound genuinely depends on the support restriction.

This is numerical support for the [inferred] outer bound, not a proof.

## Issues

**M1 (major) — error bars are fit-scatter, not bounds.**
`3|A−B|` brackets only variations close to A. One alternative setting (Q = 3000) is as good as A by the corrected Casimir residual, yet it moves k = 2 by 1.3–1.8× the quoted error and k = 3 flip by 3×.
*Fix:* add a third, Casimir-certified setting (for example Q = 3000) to the budget, or take max over settings with a corrected Casimir residual ≤ 1e-8. Gate the deep-corner acceptance on that residual, not only on `|A−B| < 1e-2`. Say explicitly that the error is a model-dependence estimate.

**M2 (major, labels) — c₂ uncertainties understated.**
* k = 3: `c₂^pres = −27.88(2)` uses the B2/B3 spread, but B2 has χ² = 3e3, i.e. it is rejected. The A and N fits give −27.82 and −27.92.
* k = 2: −17.08(2), while the A and N fits give −17.12 and −17.02.
* The §7 "±" statistical errors are meaningless where χ² ≫ 1.
* The §7 table is missing the headline row (`l = 1`, β ≈ 0.4, k = 2, ν = 0, pres).

*Fix:* quote −17.08(5) and −27.88(6); restore the row; drop the statistical errors or flag them by χ².

**m1 — the Casimir check is misread.**
`casimir_check.py` computes `1−z` by subtraction. The deep-corner "≤ 1e-4", identical for A and B, is therefore the check's own roundoff (`eps/(1−z)`). With `1−z = ((1−ρ)/(1+ρ))²` the residual is ≤ 1e-8 for A. This strengthens §2.3; correct the text and use the fixed residual (M1).

**m2 — logarithms missed at `z̄ = ∞`.**
§2.2 says the expansion is "logarithmic if c1 = c2". It is also logarithmic when `c1 − c2 ∈ ℤ`. The seed `(−0.375, 0.625)` has `c1 − c2 = 1`, and there the code duplicates an exponent instead of adding a log term. The effect is ≤ 5.5e-6 on the flip component (inside the error). Fix the statement and add the log basis.

**m3 — non-analytic terms.**
* The m^{−5/2} bound is `|c_N| ≲ 6` at 2σ, i.e. about 0.35|c₂|. A term with an O(1) coefficient is allowed.
* `m^{−2k−1} = m^{−5}` is invisible behind `c₄/m⁴ ≈ 300/m⁴`. It is excluded only by the sketch.
* At k = 3 the N-fit finds `c_{7/2} = 39 ± 0.7`, which cannot be told apart from analytic truncation.
* The parity-in-m argument that justifies even-power ansätze (odd `t_n` antisymmetric) should be stated.

Say all of this in §0 item 5.

**m4 — known error sources outside the budgets.**
* The `l = 2` v2 validation deviations are not in the `l = 2` budget. They are comparable to it for the k = 3 flip.
* Quadrature error is excluded and can exceed the quoted error for m ≲ 40.

Add both to the budgets.

**m5 — the fit-residual diagnostic is blind.**
`fitres` is blind to a wrong exponent (§3 above). Do not cite it as validation.

## Applicability

**Established**, for one spinning Regge pair in `⟨OφOφ⟩` with `Δφ = 1.75`, ν ≤ 3, β ∈ {0.4, 1, 2} and m ≤ 100 (≤ 70 at β ≈ 2): the leading formula holds at k = 2, 3, with `m²(R−1)` converging, for:
* currents, both helicity components;
* stress tensors, the preserved component.

Also established, under that numerical evidence and the sketch: there is no detectable O(1/m) or m^{−5/2} term. Assumption T is plausibly unnecessary for ⟨OφOφ⟩.

**Not established:**
* the T helicity-flip component at k = 2 (and k = 3 at β ≈ 0.4) — the component the `B₃`/flip sum rules need;
* `⟨TTTT⟩` and both pairs spinning;
* the MFT/positivity normalisation at k = 2, 3;
* uniformity in ν ~ m or β ≲ 1/m;
* other Δφ;
* a rigorous remainder bound.

STATUS row 6 may move "k = 2, 3" from "not established" to "established for one pair, `l = 1` both helicities and `l = 2` preserved, within the stated ranges", with M1 and M2 fixed.
