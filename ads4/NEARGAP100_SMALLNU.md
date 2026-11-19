# The `x = 0.05` node at `Δ_g = 100` (`ν ≈ 4.5–9`): contour evaluation below the fold (STOPPED at the coordinator's request; hand-over in §8)

*Written 19 Nov 2026 by the numerics agent continuing H0_SPINNING_FUNCTIONAL.md §25 (hand-over §25.9 item 1) and
H7_I14_LIGHTSIDE.md §23.4/§23.6 item 2. Unrefereed, not committed. Updated incrementally; hand-over at the end.*

Grades: [P] proved, [D] derived, [NE] numerical evidence, [M] model, [I] open.

## 1. Saddle structure below the fold, tracked in `ν` (`code/spinfunc/n4_thimble.py`, `track`) [NE]

* The pair `A±` of critical points of the integrand's own scalar exponent `Φ_ν` (`n1_direct.phi_w`) is obtained by a
  symmetry-breaking kick of the symmetric saddle just below the fold, then continued in `ν` (Newton, `Δν ≤ 0.25`). All logs (scalar
  flow function, scalar integrand, graviton integrands) are continued along the `z`-path traced by the saddle, starting from the
  n3 convention (`W_REF` next to the crystal of the calibration `ν`, straight `w`-path). In this tracked convention the
  intersection numbers can only jump at a Stokes point.
* `Δ = 100`, `J = 4`: `Re Φ(A+) = Re Φ(A−)` to all digits at every `ν` (`−28.83 … −35.94` for `ν = 10 … 5`). `Im Φ(A+) − Im Φ(A−)`
  grows from ≈ 0 at the fold to 1.45 / 2.17 / 2.97 / 3.84 at `ν = 8 / 7 / 6 / 5` (two saddles of equal weight with a relative phase
  of < 1 period: the oscillatory side of an Airy-type transition).
* Below the fold the saddles are far from the crystal in `w` (`|w| ≈ 1.2` at `ν = 5`), soft `−Hess` eigenvalues 1.8–13
  (`|λ|`), two hard ones (≈ 1.5·10³, 4.1·10³).

## 2. Intersection numbers at `ν = 10` [NE, weak]

`Δ = 100`, `J = 4`, scalar: thimble sums (2¹¹ points each, 12 % / 25 % RQMC errors) vs the n3 hybrid contour (`ν_f = 13`, 1 %):
`|Y₊ + Y₋ − Y_H|/|Y_H| = 0.096`; all other combinations with `n± ∈ {−1, 0, 1}` are off by ≥ 0.50. So `n₊ = n₋ = 1` (consistent
within the thimble errors; not a precision test).

## 3. The thimble integrals: no sign problem, but a rugged Jacobian; common random numbers [NE]

* **No sign problem.** On the `A±` thimbles (n3_hybrid contour built on `A±` with `ν_f = ν`, `ρ₁ = 0.8`) the weights have nearly constant
  phase: `|Σw| / Σ|w| ≈ 0.99` at `ν = 5` (`Δ = 100`, `J = 0, 4`). On the old `ν_f = 13` contour the same integral is a 50–100 %
  cancellation (`ν = 5.45`: 56 % RQMC error at 2¹³ points); the generalised-thimble step (`N3_TFLOW = 0.3`) does not help (62 %, 25 % of
  the points cut).
* **But the thimble Jacobian is rugged.** `|det ∂ξ/∂ζ|` of the flowed part varies by factors 3–10 over `δζ ≈ 0.02–0.1` (checked: RK4
  tangent propagation = finite differences of the parametrisation, `76.05` vs `76.10`). The flow lines separate along ridges between
  sinks (gradient `|∂Φ|` stays ≈ 1.5–2 on one line and rises to 30 on a neighbour 0.05 away); the spikes are not near a critical
  point and not near `A = 0` (`|A| ≈ 2` on the spike line). Single samples carry up to 10 % of the integral at 2¹⁰ points.
  Consequence: relative RQMC errors 3–25 % per thimble at 2¹³ points (`ν = 5–10`), not improvable by Gaussian proposal adaptation.
* **Cut lines carry negligible weight.** Lines cut by the per-step winding test (`WCUT`, 1–15 % of the points) were re-done with a
  winding-controlled adaptive flow (`n4_thimble._flow_rho_adaptive`, steps halved until every log moves < 0.4 rad): same mean to 7
  digits at 40× the cost. So the default (n3) flow is kept; the graviton log replay is sub-stepped where a segment winds > 0.4 rad
  (`grav_on_paths_sub`).
* **Common random numbers (CRN) rescue the ratio.** The ruggedness is a property of the geometry of the soft block, which is nearly the
  same for the state (`Δ ≈ 100–130`) and the reference (`Δ = 2·10⁴`) at equal `ν` (in the normalised `ξ` coordinates). With the same
  proposal and Sobol seeds, at `ν = 5`, `J = 4`, `A+`, fixed `N(0,1)` proposal (deliberately bad: one replica 30× the others): per-thimble
  relative errors 77 % / 86 %, state/reference ratio 5.8 % (and 2 % without the outlier replica). The production runs therefore compute
  state and reference on common random numbers (`n4_run.run_pair`: proposal adapted on the state, re-used for the reference).

## 4. Scalar validation against the exact Mellin object (`n4_validate.py` → `cache/n4_validate.txt`) [NE]

**Comparison object.** For every block `D` the exact `X_D(ν) = g C(ν+i) + g̃ C(ν−i)` (n3_exact: closed-form part at complex `ν` +
`ρ₂` continued from the exact real-`ν` nodes) is computed from e1 records, or — for blocks not in e1 (the `Δ = 100` row, the extra
rows at the chain's exact `Δ`, and `Δ_ref = 2·10⁴`) — from new exact records `n4_exact.py` (`exactskel.skelexact.ExactState`,
read-only import; reproduces e1's `ρ₂` to 10⁻¹⁵ and `X` to 10⁻¹¹ on an e1 state; ≈ 1–4 s per block).
`d_D = log Y_D − log(q_D X_D)` must be independent of `D` at fixed `(ν, J)` (CLS-D normalisation and the `J`-only constant of `Q_J`).
This holds **as a complex number** in the n3 branch convention: on the validated hybrid points (`ν = 9.197`, `J = 4`, `m/Λ =
1.05–1.3`) `Re d` is constant to 5·10⁻⁴ and `Im d` to 2·10⁻⁴. The test statistic is `d_D − d_ref` with `D_ref = 2·10⁴` — exactly the
state/reference ratio used by the graviton chain — in c2 units `Δc2 = m_D² Re(d_D − d_ref) ex₂/E₂`; errors by jackknife over the 8
common-seed replicas (heavy tails: the jackknife is not reliable at the 2σ level).

| `ν` (`J = 4`) | blocks | `A₊ + A₋` (MAIN): `Re(d − d_ref)`, `Im(d − d_ref)` | `Δc2` | `A₊` alone / `A₋` alone (`Δc2`) |
|---|---|---|---|---|
| 9.200 | 105.0, 109.9, 119.8, 129.6 vs 2·10⁴ | `−0.0018 ± 0.0056`, `+0.027 ± 0.021` (1.05Λ) … `−0.0011`, `+0.017` (1.3Λ) | **−16 ± 51** (all four) | −59 … −72 (±7·10³) / −127 … −146 (±85); `Im` off by ∓0.08, ±0.10 |
| 5.451 | same | `+0.045 ± 0.025`, `−0.106 ± 0.023` (1.05Λ) … `+0.026 ± 0.015`, `−0.064 ± 0.014` (1.3Λ) | **+453 ± 255 … +403 ± 237** | +2230 (±46) / +2303 (±30) |

* At `ν = 9.2` the thimble sum `A₊ + A₋` reproduces the exact state/reference ratio (modulus to 0.2 %, `Δc2 = −16 ± 51`), while
  either thimble alone is wrong (phase off by ±0.1 rad, modulus by 1 %): **this confirms `n₊ = n₋ = 1` below the fold, in the
  tracked branch convention** [NE].
* At `ν = 5.45` the residual is `Δc2 ≈ +430 ± 250` (1.8σ) in the modulus and `−0.064 … −0.106` rad (4.7σ by the jackknife) in the phase;
  all four states show the same pattern `(+470 − 1100i)/m²` (common seeds make the statistical error itself ∝ `1/m² − 1/m_ref²`, so
  this pattern does not distinguish a fluctuation from a systematic). Either thimble alone is off by `Δc2 ≈ +2200`.

* **Odd `J`: no scalar test.** For odd `J` the scalar `n0 = +1` crystal integral is not the identical-scalar functional (whose
  `S⁰ = {n0 = ±1}` sum vanishes for odd `J`), so `d_D − d_ref` is meaningless there (it comes out O(1)); the odd production points
  carry no exact scalar check. The odd graviton integrand (P8.build, `n0 = +1`) is the chain's object as before.

## 5. Production (`n4_run.py`, `n4_prod.py`; caches `n4_raw.pkl`, `n4_pts_<kind>_J<J>.pkl`, `n4_odd_pts_<kind>_J<J>.pkl`) [NE]

* 17 state/reference pairs at the `x = 0.05` node (`ν = 0.05 m`): even `J = 4, 6, 8, 10` on the `Δ = 100` row (`C2`, `C3H`); odd
  `J = 5, 7, 9` on the `Δ = 100` row and on the extra rows (p8_k2 keys: `J = 5` at 1.05/1.1/1.2/1.3Λ, `J = 7, 9` at 1.1/1.2/1.3Λ; at
  1.05Λ only `J = 5` exists above the gap), each with its `Δ_ref = 2·10⁴` partner at the same `ν`.
* Per block: saddles tracked from `ν_c = 10` (the `2·10⁴` track is shared by all `ν` of a `J`); thimbles `A±` (`ρ₁ = 0.8`), 2⁹·2·8 =
  8192 RQMC points each, proposal adapted to the state's scalar + graviton integrands and re-used for the reference (CRN), graviton
  logs replayed along the flow paths (sub-stepped where needed). ≈ 23–35 min per pair (1 core).
* **Stored raw** (for recombination with the non-local `χ` later): `n4_raw.pkl[(mode, Δ, J, ν)]` holds, per thimble, the literal
  crystal integrals `Ψ_ν` (unshifted, at the real node `ν`) of the scalar and of each graviton kind: replica values, means, `logN`
  (graviton: n1_run convention `logN = Φ_g(saddle) + log det L`, log-Jacobian constant included), proposal, diagnostics.
  `n4_(odd_)pts_*` hold the sums `A₊ + A₋` in n1_run / p8 format (key `(Δ, round(ν, 9))`) for the chains. No shift is applied here;
  the chains apply their own (local) de-shift.

* **Status of the run.** Stopped at the coordinator's request after 13 of the 17 pairs: **all odd pairs done** (`Δ = 100` row
  `J = 5, 7, 9`; all 10 extra-row pairs); **the 4 even pairs (`J = 4, 6, 8, 10`) were not computed** (no `n4_pts_*` files).
  `n4_prod.py 0 2` / `1 2` resumes them (skips finished pairs).

## 6. Graviton values at `x = 0.05` (odd sector, `P8.build`, Ĉ₂′ local de-shift; `h7i14/v5_odd100_smallnu.py`) [NE, with caveats]

**CRN errors of the state/reference ratio** (`cache/n4_gamma.txt`; jackknife over the 8 common-seed replicas, in c2 units
`m² × rel. error`): 130–460 c2 for 18 of the 26 (pair, kind) values; 550–770 c2 for J = 7 1.3Λ C2 and J = 9 1.1/1.2Λ C2; 1.7–2.9·10³ c2
for J = 9 on the `Δ = 100` row and J = 9 1.1Λ C3H; **unusable (≥ 100 %)** at J = 9, 1.3Λ. p8's own error column (state and reference
spreads added as if independent) overstates these 5–20×; p8's acceptance rule therefore had to be told the CRN error
(`v5_odd100_smallnu.crn_err`: state record err := |E| × CRN error, reference err := 0; a point with CRN error above `N4_ACCEPT`
(default 0.05) is marked unusable and the chain falls back to `K₂` there).

**Γ at `x = 0.05` (Ĉ₂′ chain, local de-shift; c2 units; direct vs `K₂`, same `a_h`)** — `cache/v5_odd100_smallnu_xrows_acc005.txt`:

| `J` | row | C2 direct / K₂ | C3H direct / K₂ | CRN err C2 / C3H |
|---|---|---|---|---|
| 5 | Δ=100 | −48 / +82 | −271 / −224 | 190 / 320 |
| 5 | 1.05Λ | +124 / +383 | −666 / −660 | 130 / 260 |
| 5 | 1.1Λ | +270 / +461 | −1183 / −1041 | 180 / 360 |
| 5 | 1.2Λ | +500 / +548 | −9612 / −9408 | 240 / 190 |
| 5 | 1.3Λ | +60 / +257 | +1497 / +1598 | 290 / 140 |
| 7 | Δ=100 | +11 / +988 | −1522 / −1459 | 220 / 240 |
| 7 | 1.1Λ | +792 / +1675 | −3445 / −3712 | 130 / 290 |
| 7 | 1.2Λ | +900 / +1687 | −32629 / −32619 | 460 / 320 |
| 7 | 1.3Λ | −890 / +589 | +5813 / +5860 | 770 / 330 |
| 9 | Δ=100 | (K₂ −581; direct CRN err 30 %) | (K₂ −508; 18 %) | 2900 / 1700 |
| 9 | 1.1Λ | (K₂ +2747; 6.5 %) | (K₂ −6538; 17 %) | 760 / 1960 |
| 9 | 1.2Λ | +508 / +2724 | −60593 / −60388 | 550 / 410 (rel. 3.9 % / 3.0 %) |
| 9 | 1.3Λ | (unusable) | (unusable) | — |

(The CRN errors apply to the raw crystal ratio; the Γ entries for rows where Γ is huge (`C3H` at 1.2Λ) carry the same *relative*
error, i.e. ±1.8·10³ on −60.6·10³.)

* **`Ĉ₃` (`C3H`): direct = `K₂` within errors at every accepted node** (ratios 0.94–1.14; at the large values 1.002 / 1.000 / 1.003).
* **`Ĉ₂′` (`C2`): direct is systematically smaller than `K₂`** (same direction as at `x = 0.1, 0.15`, H7 §23.4: `K₂` overstates the
  small-`x` `Ĉ₂` correction), by 200–2200 c2, i.e. 1–5× the CRN errors.
* **The huge `C3H` values at 1.2Λ (−9.6·10³ / −32.6·10³ / −60.6·10³ for J = 5/7/9) are confirmed by the direct crystal values.** They are
  not a `K₂` artefact. Their size comes from the chain's normalisation of the node, not from the crystal: the 1.2Λ node sits at
  `ν = 5.910`, 0.07 from the scalar `ν₀ = √143/2 = 5.979`, where `Ĉ₃`'s local de-shift (scalar `g`) is singular (H7 §23.5). The sign
  flip of `C3H` between 1.2Λ and 1.3Λ (−60.6·10³ → +11.3·10³ at J = 9) fits a pole of the de-shift between `ν = 5.91` and `6.40` [I:
  interpretation, not checked separately].

**`n`, odd sector, `Δ_g = 100` (Ĉ₂′):**

| protocol | `x = 0.05` values | odd min `F/N_s` (state) | `x < 0.075` window | `t*` |
|---|---|---|---|---|
| main (§6 rows), H7 §23.4 MAIN (n3 at 0.1/0.15, `K₂` at 0.05) | `K₂` | 0.2392 (J=5, 1.05Λ) | (−∞, 4.94] | (−6.9, 29.1) |
| **main + n4 direct at `x = 0.05`** (J = 5, 7 accepted; J = 9 falls back to `K₂`) | direct | **0.2426** (J=5, 1.05Λ): **n = 1** | (−∞, 4.94] | (−6.9, 29.1) |
| extra rows, H7 §23.4 | `K₂` | −0.550 (J=9, 1.2Λ) | [−0.93, 0.60] | (−1.55, 0.64) |
| **extra rows + n4 direct** | direct where accepted | **−1.345 (J=9, 1.2Λ): n = 0** | **[−0.93, +0.36]** | (−1.55, 0.43) |

(`N4_ACCEPT = 1.0`, i.e. accepting every n4 point with CRN error ≤ 100 %, gives the same numbers.) With the `x < 0.075` part set to
0, the extra-row minimum is +0.2447 (J = 5, 1.05Λ); with it ×(−1), −0.062 (J = 9, 1.1Λ).

* **The extra-row failure is confirmed by direct contour values, and made worse** (−0.55 → −1.35): it is the J = 9, 1.2Λ node, where
  the direct `C3H` agrees with `K₂` to 0.3 % and the direct `C2` (+508 ± 550) removes part of the `K₂` C2 value (+2724) that had
  partly compensated it. It is a property of the local de-shift next to `ν₀ = 5.98` acting on the correct crystal values, not of the
  `K₂` crystal approximation. [NE for the crystal values (no scalar check is possible for odd `J`, §4); the de-shift diagnosis [I].]
* So the open question at `x = 0.05` is definitional (non-local `χ` instead of the local de-shift near `ν₀`, H7 §23.2/§23.6), not
  numerical. The raw values for that recombination are stored (§5).

## 7. Grades

| statement | grade |
|---|---|
| below the fold the crystal integral = thimble(`A₊`) + thimble(`A₋`) (`n₊ = n₋ = 1`) in the `ν`-tracked branch convention | [NE]: thimble sum vs n3 hybrid at `ν = 10` (9.6 %, within errors); vs exact Mellin at `ν = 9.2` (`Δc2 = −16 ± 51`; single thimbles off in phase by 0.1 rad) |
| the contour method reaches `ν ≈ 4.9–6.4` | [NE], partly: scalar at `ν = 5.45` (`J = 4`, 4 states vs 2·10⁴) `Δc2 = +400–450 ± 240–255` (1.8σ) and phase residual −0.06 … −0.11 rad (4.7σ by the jackknife, which is unreliable with these tails). Not a clean pass; a precision statement below `ν ≈ 9` is not made |
| no exact scalar check at odd `J` | [D-sketch] (odd-`J` scalar `n0 = +1` integral is not the identical-scalar functional) |
| odd graviton values at `x = 0.05`, J = 5, 7 (all rows), J = 9 at 1.2Λ | [NE], CRN errors 130–770 c2; J = 9 at Δ=100/1.1Λ/1.3Λ not usable |
| odd `n = 1` at `Δ_g = 100`, main protocol, with direct `x = 0.05` values | [NE + M + I] (min 0.2426; J = 9 there still `K₂`) |
| odd `n = 0` with the extra rows, now with direct `x = 0.05` crystal values (−1.35, J = 9, 1.2Λ) | [NE + M]; cause: local de-shift next to `ν₀` [I] |
| even sector at `x = 0.05` | **not computed** (stopped) |
| combined `n` at `Δ_g = 100` | main protocol: odd min 0.2426, even unchanged from H0 §25.6b (0.2206, `K₂` at 0.05) → combined `n = 1` [NE + M + I]; with extra rows `n = 0` |

## 8. Hand-over: using the method at `ν ≈ 5–20`

* **Which representation where.** `ν ≥ ν* ≈ 11.3–12` (above the fold): the n3 hybrid contour (`n3_hybrid`, `ν_f = max(ν, 15)` or 13),
  validated at 9.2–19.2 (H0 §25.4). `ν ≈ 9–11`: either (n3 with `ν_f = 13` still works, errors 1–2 %). `ν ≲ 9`: the `A₊ + A₋` thimble
  sum (`n4_thimble`), the only representation without a sign problem there.
* **Recipe** (`n4_run.run_pair(J, ν, [Δ_state, 2e4, …], kinds)`):
  1. `n4_thimble.track(Δ, J, [ν_c = 10, …, ν])`: symmetric saddle from `ν = 20` down to `ν_c`, kick into `A±`, continue in `ν` (`Δν ≤
     0.25`) with all logs (scalar flow function, integrand, graviton kinds incl. `'odd:C2'`, `'odd:C3H'`) continued along the saddle's
     path. `ν_c` must be below the fold of that `(Δ, J)` (the kick raises otherwise). The `2·10⁴` track is slow (5–10 min): reuse it for
     all `ν` of a `J` (`more=`).
  2. `n4_thimble.thimble(rec, Δ, J, kinds, rho1=0.8, m=9, R=8, prop=…)` per saddle; sum with `exp(logN₋ − logN₊)` (`n4_run.combine_kind`).
  3. **Always compute state and reference on common random numbers** (same proposal, same seeds) and quote the error of the RATIO
     (`n4_gamma.jk_ratio`); individual values have 3–30 % errors, ratios 1–5 %.
  4. Validate each new `(ν, J)` region on the scalar (even `J`) with `n4_validate.py` (exact records for any block: `n4_exact.get(Δ, J)`).
* **Known weaknesses / next steps.** (a) The thimble Jacobian is rugged (§3); 2⁹·2·8 points give 3–25 % per thimble. Ideas not tried:
  a flow metric adapted to the soft block, stratifying the 6-sphere of directions, or the scalar as a control variate (exact for even
  `J` via `n4_exact`). (b) The `ν = 5.45` scalar residual (`+430 ± 250` c2, phase −0.1 rad) should be re-run at 4× statistics before
  the method is trusted below `ν ≈ 7` to better than ~500 c2. (c) Even-sector `x = 0.05` points: run `n4_prod.py` (resumes), then
  `N3_SPIN=1 N4_USE=1 N1_RELMAX=0.5 python3 n3_analysis.py` (→ `cache/n3_analysis_spin_n4.txt`; `N4_USE=0` = old behaviour). Note: n1's
  `grid_direct` also adds state and reference errors independently; give it CRN errors as `v5_odd100_smallnu.crn_err` does.
  (d) For the helicity OFF-DIAGONAL at `x ≈ 0.05–0.15` (the large-`J` failure of S1_S4_FIXES.md): the even graviton kinds return the
  full 2×2 matrices (`E`, `reps` per thimble in `n4_raw.pkl`), so the off-diagonal comes for free from the same runs; at large `J`
  check the fold position `ν*(Δ, J)` first (`ν_c` must lie below it).
* **Files.** `code/spinfunc/n4_thimble.py` (saddle tracking, tracked branches, thimble integral, odd evaluator, adaptive flow (off),
  sub-stepped graviton replay), `n4_run.py` (CRN pairs, storage), `n4_prod.py` (job list, 2-worker driver), `n4_exact.py` (exact
  scalar records for any block), `n4_validate.py`, `n4_gamma.py`; `n3_analysis.py` gained `N4_USE` (default off);
  `code/h7i14/v5_odd100_smallnu.py` (copy of v4_odd100 + n4 points + CRN acceptance; `N4_USE=0` reproduces v4). Caches:
  `spinfunc/cache/n4_{raw,scal,exact}.pkl`, `n4_odd_pts_{C2,C3H}_J{5,7,9}.pkl`, `n4_validate.txt`, `n4_gamma.txt`, `n4_prod_w*.log`;
  `h7i14/cache/v5_odd100_smallnu*.txt`. Nothing committed; no protected file edited.
