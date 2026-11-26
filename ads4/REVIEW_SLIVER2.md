# REVIEW_SLIVER2: adversarial referee report on ANALYTIC_SLIVER2 (L ∈ [3.65, 5], Δ_g = 1000)

*26 Nov 2026. Referee agent. No existing file edited. `prl/` not touched (it is only read, through a symlink, for `nf02_psi.json`). Nothing
committed. Review scripts are in `ads4/code/review_s2/` (`rv_*.py`). Logs are in `code/review_s2/logs/`. The reruns use a private copy of
`code/analytic` and `code/analytic2` under `code/review_s2/ws/code/`. The original `analytic2/cache` was moved aside there as
`cache_orig/`, so every rerun output, including the `rmass3` masses, was regenerated from source.*

## 0. Verdict

**No error that breaks any of the four cells was found.** Every lemma is a proof valid on the whole continuum of states in its box, not a
sampled check. I found:
* two gaps (S1, S2), both closed here;
* one wrong justification with a negligible effect (S3);
* one floating-point nit (S4);
* one load-bearing input that is still unreviewed from AS: Lemma X, which carries the right edge.

| item | verdict |
|---|---|
| Lemma F (§3) | **holds**: rigorous over the whole box (affine frames plus Arb Taylor bound with a box Hessian); nit S4 |
| Lemma B (§2, incl. multi-frame min) | **holds**: every family is rigorous for every state in the box, so the per-j minimum is a bound |
| K_cell assembly (§4) | holds: every endpoint and monotonicity direction checked; ι-range and g-box contain all states (exact check, 0/1432 failures) |
| Lemma Loc″ (§5) | **holds**: identity and proof re-derived; ψ₂ = 0 above Λ verified; J(ν) consistent between residual and remainder |
| rmass3 masses (§6) | **hold**: code identical to rt_rmass2/EPSW, regenerated bit-identically, totals consistent with EPSW |
| state relations for L ∈ [3.65, 3.75] | **gap S1, closed here by proof** (the document only scans L ≥ 3.75) |
| G1 bookkeeping above Λ | **wrong argument (S3)**: "\|f_θ\| ≤ \|f\| pointwise" is false. The omitted term is ≤ 7·10⁻¹¹·F |
| contour pieces on [3.65, 3.75] (old code) | hypotheses hold (asserted at d_max, monotone); the right edge 0.32·F rests on Lemma X (AS, prior G4 still open) |
| reruns | see §7 |

## 1. Lemma F (turning point): a proof over the box, not a sample

I re-derived the proof line by line against `tpf3.lemmaF3_log`.

**The identity.** M_iF_i = F_{i+1}diag(R₁, R₂) + [[R₁δ₁, R₂δ₂], [0, 0]] is exact for any affine R. The row-sum bound for E_i is correct.

**How the cell is covered.**
* Each R_{c,i}(p) is an affine function with point coefficients.
* sup|R| uses the expansion |R⁰ + Σg_lΔ_l|² ≤ |R⁰|² + 2Σ|Re(R̄⁰g_l)|D_l + (Σ|g_l|D_l)², which is valid for every p in the box.
* The defect δ = 2w + 1/R − R′ is bounded by Taylor with the integral remainder:
  - the first derivatives are exact at p⁰;
  - the Hessian 2∂²w + 2g_lg_m/R³ is enclosed in Arb over the whole box, with R as a ball, so the bound is infinite if R can vanish;
  - R′ is affine, so it has no Hessian.
* I checked the second derivatives of ε(m, ι) symbolically: ε_m, ε_ι, ε_mm, ε_mι, ε_ιι all agree.
* Q_{i+1} is the Arb inverse of the ball matrix. A singular ball raises an exception, and that i_s is then skipped.
* F₀ has zero slope (g = 0 at i = 0), so its point inverse is exact.
* P_{k_A} is a ball over the box.
* Hence the bound holds for every (τ, a, ι, σ) in the box and every n > K. Nothing in it is sampled.

**Phase I (Theorem E′ stopped at k_A).**
* I re-derived it:
  - the σ-excess sum over m ≥ a + K + 1 is bounded by the integral from a + K;
  - the adiabatic sum runs up to y_{k_A+1} ≤ τ/(a + K − 1), by the convexity step;
  - c₀ ≤ 1 + 2|α|/a².
* Monotone choices in `phase1_log`: κ′ is evaluated at (t_hi, a_lo), q at t_lo, σ at s_hi and c₀ at a_lo. All of these are worst cases.
* The d-ratios are taken at σ_lo, which is valid because each ratio is ≤ 1. I checked the ratio formula d_{k+1}/d_k against d_k = 1 − k(k+1)/(N(N+1)).

**State coverage.**
* The ι-range is [t_lo/a_hi·(1 − 10⁻⁴), t_hi/a_lo]. t is the exact root of t² − (2+e)t + d(1+e) = 0, and it is monotone in d and in e (I checked
  the sign of dt/de).
* It contains the b = 0 state and the helicity-shifted state (n + 4, N + 4) of the b = 4 channel, whose d₀ = d·D/(D + 20).

**S4 (nit, floating point).**
* The Taylor and affine radius is D = (hi − lo)/2 around the float midpoint p⁰.
* In 211 of 1224 boxes, [p⁰ − D, p⁰ + D] misses [lo, hi] by ≤ 1.4·10⁻¹⁵ relative. This was checked in exact rationals.
* No state is lost:
  - the a-boxes have integer (dyadic) endpoints, where p⁰ and D are exact;
  - the true ι lies strictly inside its range, with margins of 10⁻⁴ and 10⁻³;
  - in τ only a single endpoint can be missed, and sup|K| over a cell is a sup of a continuous function.
* Still, D should be inflated by about 10⁻¹².

**Adversarial DV.**

`rv_F.py`: the certificate's own base bound PP_box (Lemma F / E′ / Lemma X, plus helicity) on the real certificate boxes. The boxes are a unit
a-box and a τ-cell [T₀, T₀ + 1.1] within the L-cell. The comparison is with the exact maximum over corner and edge integer states: a = a_lo and
a_hi, b ∈ {0, 4}, n at both L-edges, and τ at both edges.

| L-cell | boxes | min log(bound/exact) | worst box |
|---|---|---|---|
| [3.65, 3.75] | 112 (T₀ ∈ {998 … 1050}, a_lo ∈ {998 … 1045}, σ ∈ {5/2, 7/2}) | **2.53** | τ ∈ [1010, 1011.1], a ∈ [1045, 1046], b = 4 |
| [3.75, 4] | 18 | 2.66 | |
| [4.5, 5] | 18 | 2.49 | |

At the binding box (a ∈ [998, 999], τ ∈ [1031, 1032.1], L ∈ [3.65, 3.75]) the bound is e^{11.23} and the exact value is e^{6.76}.

`rv_defect.py` re-implements the per-step quantities and evaluates, at 50 digits, the exact δ_{c,i}(p), |R_{c,i}(p)| and ‖E_i(p)‖ at the corners
plus 4 random points in each of 5 boxes × 2 L-cells × 2 values of i_s:
* exact/bound for ‖E_i‖ is ≤ 1 − 1.7·10⁻⁶ (at steps with g = 0, ‖E‖ equals the sup|R| part);
* exact/bound for |δ| is ≤ 0.90;
* exact/bound for |R| is ≤ 1, with equality only at i = 0, where g = 0.

**Verdict: holds.**

## 2. Lemma B (β-chain): rigorous, and the minimum over families is legitimate

* For any fixed point frames T_j, the bounds ‖c_{j+1}‖ ≤ ‖T_{j+1}⁻¹𝕄_jT_j‖‖c_j‖ and |x_j| ≤ ‖row₁T_j‖‖c_j‖ hold for every member of the ball
  matrix 𝕄_j.
* The frames are computed from box midpoints only. So each family gives a bound valid uniformly on the box, and the elementwise minimum of
  valid bounds is valid. No family needs to be valid for "the right" state.
* The transform norms are Arb `.upper()`. The inverses are Arb inverses of exact point matrices. The products are rounded to nearest in
  floats, which gives a relative error of about 10⁻¹⁴; this is covered by the final (1 + 10⁻⁹).

**State box.**
* g ∈ [g_lo, g_hi] is checked exactly (`rv_B.py` (1)): 1432 integer states (both channels, both L-edges, every a-range of `a_ranges`
  including [2·10⁵, ∞)), with **0 failures**.
* s ∈ [1, 9] enters as an interval.

**Robustness remark (not an error).**
* `chain_bounds_multi` takes `min(x, y)` with the first family's value as x. A NaN from that family would therefore propagate.
* K_cell then drops a NaN a-range silently, because `v > best` is false.
* No NaN can arise here: den is bounded away from 0 and overflow gives inf, not NaN. A `math.isfinite` assert would make this explicit.

**DV against exact values** (`rv_B.py` (2)).
* Method:
  - the K_cell chain is started from the certificate's own PP_box M₀, M₁;
  - exact values come from the downward β-recurrence (TAIL E3) in Arb at 300 bits, from n-recurrence values at σ = 5/2 and 7/2;
  - stability was checked by perturbing the start by 10⁻⁹; the drift is ≤ 10⁻³ in log.
* Results: max_j log(exact/bound) ≤ **−1.53** over 7 boxes × 2 L-cells (ν ∈ {500, 900, 1003, 1020, 1031, 1045}), b ∈ {0, 4}, j ≤ 109.
* The binding j is always −1 or 0, i.e. the base bound, so the chain itself loses less than that margin.
* (A first attempt used the n-recurrence directly at Re β ≈ −100. It loses all precision, giving a spurious e^{+850}, and was discarded. Arb's
  forward n-recurrence also blows up by wrapping. This is noted only so that later checks avoid the trap.)

**Verdict: holds.**

## 3. Per-cell kernel (§4) and residual bookkeeping (§6)

**Directions checked.**
* Every d-power is taken at d_hi:
  - d³·d^{j−5/2} is increasing;
  - dt/d = g_hi·d is increasing.
* |B| is evaluated at T₀ (decreasing).
* c_j is an Arb ball over the ν-cell.
* `ratio_bound_box` is a sup over j ≥ jmax (convex/decreasing argument re-checked).
* G = max(2, (A+B)/den) is taken at the largest |β| and at d_hi; (A+B)/den is increasing in |β|, d and s.
* The geometric tail needs rb < 0.9, which is asserted.
* The helicity PP formula is evaluated at d₀ = d_hi(1 + 10⁻⁴).
* The residual pairs r with the full K_J (no +1), consistently with Theorem R.

**Sub-cells.**
* The EPSW panels have width ≤ 0.1 for ν ≥ 1 and are aligned on tenths. Each panel starting in [T, T+1) lies in [T, T+1.1], the cell on which
  K_cell is evaluated.
* [0.2, 1) is covered by the [0.4, 1] band and the small-ν bands, each with EPSW's mass of [0, 1].
* The [897.9, 898] panel is counted twice: in rmass3 and in rt_rmass2. This is conservative.

**J(ν) consistency.**
* The residual uses J′ = 2200 for cells with ν_lo ≥ 999. `rt2_loc2` uses J′ for T ≥ 999. The panels beyond Λ are exactly [T, T+1].
* ψ = ψ₂^{NF02}(ν/Λ)/k₀ for ν ≤ Λ and 0 beyond: this is in the `ew_common` docstring, and `inside = False` for panels above Λ.
* So for ν ≥ Λ = 998.4994, f = r, and any J′ is admissible there. **Claim verified.**
* The original [3.65, 3.75] and [3.75, 4] logs agree with a run using `RT2_JP=2200`; the rerun uses that setting explicitly.

**S3 (wrong justification, negligible).**
* §6 and `rt2_resid.py` assert that "above Λ, ψ₂ = 0, so r_θ = f_θ and |f_θ| ≤ |f| = |r| pointwise (0 ≤ θ ≤ 1 on the real axis)", and
  therefore add the θ-correction only for ν < Λ.
* This is false. By TAIL_CERT R.2 and `tc_cutoff.py`, f_θ = 2Re[θ(ν − i)g(ν)], with θ evaluated at the complex point ν − i around the cutoff
  ν_c = Λ + 100 = 1098.5.
* `rv_theta.py` gives θ(1098.5 − i) = 0.49999 + 0.0566i and θ(1080 − i) = 0.9959 + 0.0018i. So 2Re[θg] ≠ θ·2Re g, and wherever Re g ≈ 0,
  |f_θ| > |f|.
* The omitted term is 2∫_{ν≥Λ}|f_θ − f||K_J|. It is at most
  2[sup_{[Λ,1059)}K · ∫|f_θ − f| + sup_{[1059,1150]}K · (∫_{1059}^∞|f_θ| + ∫_{1059}^∞|f|)], with:
  - ∫|f_θ − f| ≤ 1.66·10⁻²³ (tc_cutoff);
  - ∫_{1059}^∞|f_θ| ≤ 2.4·10⁻²⁵⁵ (CONTOUR cb_lines.tail);
  - ∫_{1059}^∞|f| ≤ 5.9·10⁻²⁴⁷ (EPSW).
* The resulting omitted term is ≤ **6.6·10⁻¹¹·F** ([3.65, 3.75]), 7.1·10⁻¹¹·F ([3.75, 4]), 4.5·10⁻¹¹·F ([4, 4.5]) and 5.8·10⁻¹²·F ([4.5, 5]).
  - These use the sup|K_J| of the original caches. On [1059, 1150] that sup is up to 4.7·10³⁹, so the full θ-weighted mass is essential there;
    the "sup|f_θ − f|" route would fail.
* The same pointwise slip affects E_ρ (rt2_loc2 uses EPSW cell masses of f above Λ). There the correction is ≤ 1.66·10⁻²³·e^{(per-cell log
  bound)}, which is irrelevant.
* **Fix:** replace the sentence by the bound above.

**G2 / far term.**
* The weighted tail `wtail` is correct:
  - the completed-square erfc formula was checked;
  - the Y-term ∫e^{wt−t²/4} ≤ e^{w²}(2/(T−2w))e^{−(T−2w)²/4} was checked;
  - ymax = 1 matches CONTOUR's Y = 0 call.
* The full-kernel bound uses |1 − ξ| ≥ 0.49d_lo − 0.255d_max, and its F_{3,ν} carries both ± terms, so the missing factor 2 (compared with
  the band terms) is correct.

**Item 4: rmass3 (EPSW per-unit-cell masses).**
* `rmass3.py` is a line-for-line copy of the panel loop of `rt_rmass2.py`; the diff shows only the range filter and the output.
* Regenerated: `rmass3_1_40.json` is bit-identical. `rmass3_40_898.json` is listed in §7.
* Consistency with EPSW's certified ranges (`ew_cert_1000_re_rep.txt`):

| range | rmass3 sum | EPSW range |
|---|---|---|
| [1, 40) | 6.085038·10⁻⁷ | 6.0850·10⁻⁷ |
| [40, 898) | 1.73921·10⁻⁷ | — |
| [40, 898) + rt_rmass2 [897.9, Λ) | 1.73921·10⁻⁷ + 2.995·10⁻⁹ = 1.7692·10⁻⁷ | **[40, Λ)**: 1.7690·10⁻⁷ |

  The last row agrees up to the duplicated panel and rounding.
* AS's `EPSW_BANDS` label "[40, 898]: 1.7690e-7" is in fact the [40, Λ) range. That is conservative, and the coarse runs only use it as an upper
  bound.

## 4. Lemma Loc″ (§5)

**The proof.**
* `rv_loc.py` verifies the identity ((x−½)² + h)² − x²(x² + 4l) = −2x³ − qx² − 2x(¼+h) + (¼+h)² symbolically (sympy).
* With qx² ≥ (¼+h)², P(x) < 0 for x ≥ J + 1.
* Arb check of the ratio bound on all 2294 (J, ν-cell) pairs where `log_loc4` sets R_t = 1: the maximum is 0.999955 < 1 over j ≤ J + 2·10⁴.
* Since P < 0, the ratio stays < 1 for all larger j.
* The geometric sum with r < 1 is then correct.

**The rest of Loc″.**
* Larger (b₁, θ₀, θ₁) keep r₁ < 0.995.
* The radii |1 − ξ| ≤ (1 − Bρ) + |θ| remain valid for θ₁ ≤ 0.5.
* frac(θ) is decreasing in θ.
* The B-tier weights are unchanged.
* The k = 3 piece uses real cells only.

**J′ = 2200 for ν ≥ 999:** legitimate, see §3.

**Verdict: holds.**

## 5. State relations below L = 4 (item 5)

**S1 (gap, closed here).**
* The document claims the theorem down to L = 3.65.
* But §5 and §8 state the relations as "verified for L ≥ 3.75", and `check_state.py` defaults to L ≥ 4.
* Several hypotheses depend on them:
  - Loc′/Loc″ (0.49d ≤ a/N ≤ 0.51d);
  - REC and TAU_COEF (Dd² ≥ 3.9a²);
  - `state_box` (e ≤ 9·0.51·d/a);
  - rt_pp (N ≥ 40a).

**Proof for every state with L ≥ 3.65** (a ≥ 998, b ∈ {0, 4}, n ≥ 997). This is exact algebra in the style of REVIEW_SLIVER_L5 §3, with
d = a(2N − a + s)/(N(N + s)) and x = a/N.
1. **Lower bound.** a/(Nd) = (N+s)/(2N−a+s) ≥ ½.
2. **First bound on x.** x(2 − x) ≤ d(1 + s/N), with s/N ≤ 9/1995 and d ≤ e^{−3.65} = 0.025991. This gives x ≤ 0.013141.
3. **Bound on s/N.** Hence N ≥ 75948 and s/N ≤ 1.185·10⁻⁴.
4. **Upper bound on a/(Nd).** x ≤ d(1 + s/N)/(2 − x), so **a/(Nd) ≤ 0.50337 ≤ 0.51**.
5. **D relation.** Dd²/a² = (2 − x + s/N)²/(1 + s/N) ≥ **3.9471 ≥ 3.9**.
6. **N relation.** N/a ≥ 76 ≥ 40.

The exact-rational scan (`rv_state2.py`: 70 values of a up to 10¹⁰, both b, 4 values of L in the cell, n and n ± 1) finds extremes
a/(Nd) ∈ [0.50298, 0.50332], Dd²/a² ≥ 3.94785, N/a ≥ 76.44. `check_state.py 3.65` agrees, and gives K₀ = 1.00033 and I_W/bound = 2.10 at
L = 3.65.

**Note.** §5's "Dd² ≥ 3.95a²" is false at L = 3.65 (min 3.9478), but nothing uses 3.95.

**Top line Y = 6.25 (k₀ = 4).**
* The B poles are at ν = ik on the imaginary axis, outside [δ, Λ] × [0, Y] (δ = ¼). The closest pole, 6i, is at distance ≥ ¼ from the top
  line and from the left edge.
* For ψ₃, Y = 7.25 and the pole 7i is likewise at distance ≥ ¼.
* `_digamma_box` shifts by ⌈12.5⌉ + 8 steps; every 1/(z + k) term stays at distance ≥ ½ from 0.
* The chain indices were checked: base = ch[k₀:] and the weight d^{j−¼+k₀}.

## 6. Contour pieces on [3.65, 3.75] (item 6)

**The asserted growth constants.**
* (A+B)/den = 1.9745 ≤ 1.99 and F/den = 0.9300 ≤ 1 at |β| = hypot(1503, 1151) + 2 and d = e^{−3.65}. They are increasing in |β| and d, so
  they hold throughout the cell.
* Every |β| on the contour is ≤ hypot(1503 + 7.25, Λ) ≈ 1812 < 1895.

**The j-sum tails.**
* The ratio is ≈ 9.35 at j = 59 and |ν| ≈ Λ, so rb = 9.35·2·0.02599 ≈ 0.486 < 0.5. This passes, but only narrowly. It is asserted at d_max
  and is monotone in d.
* It is the first assertion to fail below 3.65, as §8 says.

**Endpoints.**
* Every piece takes d at e^{−L_a}. The exponents of d are positive on every piece, and F_cert is the minimum over the β-cell.

**The right edge (0.320·F) rests entirely on AS's Lemma X plus the old absolute-value chain.**
* Lemma X's proof was not re-derived in REVIEW_SLIVER_L5 (G4), and it was not re-derived here either.
* DV (`rv_X.py`): PP(σ ≤ σ_top, τ ∈ [Λ − 1, Λ]) against exact max|Ŝ| over a ∈ [998, 1200], b ∈ {0, 4}, L ∈ {3.65, 3.75}, σ′ ∈ [0, σ_top]. The
  bound exceeds the truth by **e^{8.35} or more**.
* So the bound is safe in practice. But Lemma X is now the single unreviewed analytic input carrying more than half of the [3.65, 3.75] ρ.
* If Lemma X were withdrawn, the right edge could be redone with Lemma F at complex σ, as §8 suggests.

## 7. Reruns from source (item 7)

RERUN_TABLE_PLACEHOLDER

## 8. Per-cell survival

| cell | claimed ρ | rerun ρ | issues affecting it | survives? |
|---|---|---|---|---|
| [4.5, 5] | 0.0267 | RR45 | S3 (≤ 6·10⁻¹²) | **yes** |
| [4, 4.5] | 0.0894 | RR40 | S3 (≤ 5·10⁻¹¹) | **yes** |
| [3.75, 4] | 0.223 | RR375 | S1 (closed), S3 (≤ 7·10⁻¹¹); Lemma X carries 0.068 | **yes** |
| [3.65, 3.75] | 0.545 | RR365 | S1 (closed here: needed), S3 (≤ 7·10⁻¹¹); Lemma X carries 0.32 | **yes, conditional on Lemma X (as for AS)** |

## 9. Files

* `code/review_s2/rv_F.py`: Lemma F, via PP_box on the certificate boxes, against exact values. Logs: `logs/rv_F_*.log`.
* `code/review_s2/rv_defect.py`: per-step δ, |R| and ‖E_i‖ bounds against 50-digit values at interior and corner points. Logs:
  `logs/rv_defect*.log`.
* `code/review_s2/rv_B.py`: exact g-box containment, and the Lemma B chain against exact values from the Arb β-recurrence. Logs:
  `logs/rv_B_*.log`.
* `code/review_s2/rv_loc.py`: the Loc″ identity (sympy) and the ratio bound on all cells. Log: `logs/rv_loc.log`.
* `code/review_s2/rv_state2.py`: proof constants and an exact scan for L ∈ [3.65, 3.75]. Logs: `logs/rv_state2.log` and
  `logs/check_state_3.65.log`.
* `code/review_s2/rv_theta.py`: θ(ν − i) values and the omitted-term bound for S3. Log: `logs/rv_theta.log`.
* `code/review_s2/rv_X.py`: Lemma X at the right edge against exact values. Log: `logs/rv_X.log`.
* `code/review_s2/scripts/runq.sh` and `jobs.txt`: the rerun queue. Logs: `logs/R*.log`, `logs/L*.log`, `logs/C*.log`. Outputs:
  `ws/code/analytic2/cache/`.
