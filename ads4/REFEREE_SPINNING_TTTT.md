# Referee report: `SPINNING_TTTT.md` (⟨TTTT⟩, `l = 2`)

Referee: independent, hostile; 16–17 Nov 2026. Read-only use of `ads4/code/fourpt/tttt/` and its caches. The 2-hour production runs were
not repeated. Scratch scripts are in `ads4/code/referee_tttt/`:
* `ref_t13.py`: reproduces the fits, then runs mutations and diagnostics;
* `ref_rcut.py`: weight of the large-`r` grid points;
* `ref_spot.py`, `ref_spot2.py`: recompute the waves at single grid points in float64, long double and with the old `arccosh`.

Outputs are in `referee_tttt/cache/` (git-ignored).

## Verdicts

* **Correctness.** Mostly sound. I reproduced the `t13` numbers exactly from the saved `.npz` caches, for example same–same one-scale
  `1.0828 / 1.0239`, `Σ₋` `0.9215 / 0.9736`, factorisation `0.849 / 0.941`. I recomputed 22 grid points from scratch in float64, and they
  match the caches bit for bit, so the caches were produced with the fixed `arccosh` code. The precision claim, however, is wrong (issue 1).
* **Soundness.** The end-to-end test is a real test: no constant is fitted into the prediction, and the mutations break it. The `1/m²` claim
  and some "prediction-free checks" are overstated (issues 2–4).
* **Labels.** Several are inflated (issues 3–6).
* **Applicability.** This gives the factorised leading heavy `⟨TTTT⟩` moment at `β = 1`, `k = 6`, `ν = 0` to a few percent at `m = 70`,
  for parity-even exchanges only. The `1/4` constant and the absence of an `O(1/m)` term are not established end to end. Neither is any
  statement at the `O(1/m)` precision that the graviton bound needs.
* **Novelty.** The helicity structures and the MFT factors `A_h` are CHL's (2102.08160). The new content is:
  - an independent Fock confirmation of the `A_h`;
  - the identity `κ_h A_h = 2^{3J+3/2}` (symbolic for `l = 1`, four points for `l = 2`);
  - the Taylor-mode spinning partial waves and the numerical end-to-end check.

  There is no conflict with KKSD (1809.05111) or CLS-D (2311.04271), but the document does not compare with them. That comparison is missing.

## What "fit/pred" fits

`t13.fit` makes two fits:
* a least-squares fit of **8 real class amplitudes** (the parity classes) to 625 complex modes, that is 1250 real equations;
* a **1-parameter** one-scale fit.

The prediction itself has no free constant. The design matrix is well conditioned (condition number `1.1` after column normalisation). The
columns have very different weights, though: the `ff_oo` column norm is `1e-6` of the `ss`/`oo` columns, and the `pf`/`fp` columns are
`7e-3`. So the `ff_*` ratios (`−25.3`, `1.41`) carry no information, and the document is right to drop them.

## Mutation tests (saved data, `m = 40 / 70`)

| mutation | effect |
|---|---|
| (M1) swap `λ_same ↔ λ_opp` in `T_pred` | same–same `0.34 / 0.54`, opp–opp `3.01 / 1.90`, `Σ₋` `−0.37 / −1.15`. **Breaks.** The test resolves the `1/A_h` dependence. |
| (M2) flip the sign of the `g₁² S_L S_R` (opposite-sign) term in the design | diagonal one-scale fits `0.13`, residual `0.99`; `Σ₋` `−0.92 / −0.97`, while `Σ₊` is unchanged. **Breaks.** The sign is tested. |
| (M4) swap the predicted helicity classes | one-scale fits `0.05`–`0.14`, residual `0.99`. **Breaks.** |
| (M3) `g₁ → −g₁` | no change, because the `R S` cross terms drop out at `ν = 0`. The sign of `g₁` is **untested**. |

## Issues

1. **[major] Float64 precision is misstated.** The document says float64 is accurate to `≤ 3·10⁻⁵` pointwise (§4.3, §4.5), but `t11` checked
   only `r ≤ 0.64`, and the grid reaches `r = 0.991`.

   Float64 against long double, relative to the largest mode at each point:

   | `r` | error |
   |---|---|
   | 0.86–0.88 | `1e-4` (ss) to `0.65` (os) |
   | 0.92 | `1`–`18` |
   | 0.95 | `2e2`–`1.6e3` |
   | 0.99 | `O(1)` |

   Both masses behave the same way. Long double is itself unverified there. The weights suppress these points:
   * dropping `r > 0.9` changes the ratios by `≤ 2e-4`;
   * `Σ|terms|` from `r > 0.95` is `≤ 4e-6`.

   The band `0.85–0.9` carries about 1.5 % of the moments at `m = 70`, with float64 errors up to about `0.2`. So the contamination is
   `≲ 0.3–0.5 %` on the mixed-helicity and opp–opp waves at `m = 70`. The percent-level conclusions stand. The sub-percent figures do not:
   opp–opp `1.0015` and mixed `0.998` are not resolved.

   The bug fix itself is correct. The old expansion gives `nan` or deviations up to `10²` at `r ≥ 0.87`.
2. **[major] The `1/m²` claim is cherry-picked.** `m²(ratio − 1)` at `m = 40 → 70`:

   | quantity | `m = 40` | `m = 70` | consistent with `c/m²`? |
   |---|---|---|---|
   | same–same | `128` | `117` | yes; also matches the `t14` vertex value `119–124`, which is good independent corroboration |
   | mixed-helicity | `−92` | `−88` | yes |
   | `Σ₋` | `−121` | `−130` | yes |
   | opp–opp | `−66` | `+9` | no: sign change |
   | mixed wave | `−48` | `+40` | no: sign change |
   | `Σ₊` | `−42` | `+54` | no: sign change |

   The vertex `E_opp − 8` is positive (`t14`), yet the opp–opp deviation at `m = 40` is −4.3 %. That is an unexplained effect of about −12 % in
   the opposite direction. §0 ("deviations shrink 3.3–3.5×") should be marked [I] and qualified.
3. **[major] Helicity purity in the moments fails at `m = 40`.** In the mixed-helicity waves, the fitted `ss`/`oo` classes, which are predicted
   to be zero, carry 14 % each of the norm, and `|meas|/|pred| = 1.34`. At `m = 70` this is 1.3 % each, a factor of 11. The mixed wave's
   `pf/fp` classes fit at `0.643` at `m = 40`, and this number is not reported. "Dominant modes 1–8 %" hides this.
4. **[minor] `T^{so} = T^{os}` is not a check.** The document's own table says the two waves are identical by pair exchange. That makes the
   equality a symmetry, not a "prediction-free check" that the Weyl-odd term cancels.
5. **[major, label] The `1/4` constant.** The constant is not circular, but it is single-sourced. It rests on:
   * `E = 8`, a Richardson extrapolation along one ray (`β = 0.4`);
   * the `[ID]` closed form `A_h P_FK/2304`, which was fitted on 39 Fock entries with `N ≤ 6` and `Δ ≤ 10`, then extrapolated to heavy `Δ`.

   The end-to-end run uses unit-normalised `X_h` and **does not test `p[X_h]`**. The §4.5 phrase "through `p[X_h]`, the 1/4 per quadruple" is
   therefore false; the end-to-end run checks `E` only. Rows 7 and 9b should be [NE/ID], not [D].
6. **[minor] Unreported `O(1/m)` term.** `m⁸T^{hh}/λ_h` has an imaginary part of `O(1/m)` (`Im E_opp·m ≈ 46`, constant). It is odd under
   `h → −h` and cancels in the class sums, but it is not reported. The header says "not committed", but commit `70664cf` exists.
7. **[minor] Two masses at one `β`.** The [DVE] label for "no `O(1/m)`" is not justified. The document's own Marks section says [I], and §0
   should say the same.
