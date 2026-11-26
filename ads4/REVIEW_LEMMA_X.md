# Review of Lemma X (turning-point crossing, b = 0): ANALYTIC_SLIVER §5.5 / `rt_pp.cross`

Adversarial review, 26 Nov 2026. Scripts are in `code/review_lemx/` and logs in `code/review_lemx/logs/`. No existing file was edited.

## Verdict

* **Lemma X holds.** I re-derived the proof step by step, and the certificate's form of it (`rt_pp.cross`, called by `rt_pp.PP0` with c ∈ {0.7, 1, 1.5}) is a valid upper bound.
  * I found five textual imprecisions (N1–N5 below). Each one is absorbed by slack that the code already carries.
  * An independent, fully rigorous re-implementation (`lx_common.bound_indep`) is never larger than `rt_pp.cross`. The minimum of cross − indep is +0.066 nats, over every PP0 box [t−1, t] with t ∈ [949, 1150] and σ ≤ 6.5.
* **On the right edge of [3.65, 3.75] the turning point is never crossed.**
  * There τ = Im β = Λ = 998.49987 and a ≥ 998, so every step has y_k < τ/(a+1) ≤ Λ/999 = 0.99950 < 1.
  * Only phase I (the eigen-coordinates of Theorem E′) and phase II (the window) are active. The phase-III terms enter as non-negative slack (lP = 0.39, lfin = 0.995, G = 0.007).
* **Every hypothesis holds on the whole cell.** This covers all a ≥ 998, all spins, L ∈ [3.65, 3.75], y ∈ [0, Y], and σ ≤ 5.5 including the helicity shift.
* **Numerical margin.** Lemma X's bound sits at least **e^{7.727}** above the exact b = 0 multiplier on the cell. The gap is certified with Arb balls at 12 000 and 16 000 bits. PP with helicity, b = 4, has a margin of at least e^{8.21}.
* **Replacement option.** The reviewed base-multiplier lemma (`code/norming/base_bounds.py`, README §3) could replace Lemma X on this edge. It is 4–6.5 nats sharper there.

## 1. Statement and hypotheses (as used)

S_n(α) = ₃F₂(α, −n, n+1; −N, N+1; 1) with real N = n + a. The spin recurrence is S_{k+1} = S_{k−1} + ε_k α S_k, with S_{−1} = S_0 = 1 and ε_k = 2(2k+1)/((N−k)(N+k+1)).

The normalised multiplier is Ŝ = d^α S_n, where d = 1 − n(n+1)/(N(N+1)). The step parameter is w_k = ε_kα/2 = u_k + i y_k.

**Claim (`rt_pp.cross(σ, t_lo, t_hi, c)`).** The bound is uniform over:
* all real a ≥ 998;
* all n ≥ 1;
* 0 ≤ Re α ≤ σ;
* |Im α| ∈ [t_lo, t_hi] ⊂ [900, 1150].

It reads |Ŝ| ≤ exp(lI + lII + lP + G + lad + lsig + lfin), with δ = c·t_lo^{−2/3}, y_A = 1 − δ and y_B = 1 + δ. The code evaluates it at a = 998. Every term is non-increasing in a; I checked numerically that log cross falls from 14.52 to 13.12 as a goes from 998 to 4500.

Implicit facts the proof needs, each checked:
* (i) y_k is increasing, and y_{n−1} < τ/(a+1).
* (ii) u_k ≤ σ/(a+1).
* (iii) d > 0.
* (iv) |w| = q·y with q = |α|/τ.

## 2. Re-derivation, step by step

| step | claim | status |
|---|---|---|
| Eigen-coordinates | c_{k+1} = P_{k+1}^{−1}P_kΛ_kc_k, with P_{k+1}^{−1}(P_k − P_{k+1}) = (2r_{k+1})^{−1}[[−Δλ₊, −Δλ₋],[Δλ₊, Δλ₋]]. Max-norm growth ≤ \|λ₊\|(1 + (\|Δλ₊\|+\|Δλ₋\|)/(2\|r_{k+1}\|)). \|λ₋\| = 1/\|λ₊\| ≤ 1 for u ≥ 0. | correct; audit P1 slack ≥ 0 (exact to rounding) |
| log\|λ₊\| ≤ u/√(1−y²) (y < 1) | integrate asinh′ from iy to iy+u; \|1+(x+iy)²\| ≥ 1 − y² | correct; audit P2 (equality at σ = 0) |
| d-normalisation | log(d_k/d_{k+1}) = −log(1 − 2(k+1)/P_k) ≥ ε_k/2 | correct |
| Δε_k ≥ ε_k²/2 | Q(2Q + (2k+1)(2k+2)) − (2k+1)²(Q − 2k − 2) = 2Q² + (2k+1)Q + (2k+1)²(2k+2) > 0 | correct (proved; 7991 exact-rational checks) |
| excess Σ(σε_k/2)((1−y_k²)^{−1/2} − 1) ≤ σ log(2/(1+√(1−Y²))) | σε_k/2 = σy_k/τ ≤ σΔy_k/y_k, increasing integrand, left Riemann sum; Y = y_A + Δy_max | correct |
| adiabatic sum ≤ q[asin Y + (q/2)log(1/(1−Y²))] | \|dλ/dw\| ≤ 1 + \|w\|/\|r\| (maximal at the upper y-end in phase I), \|r\| ≥ √(1−y²), ε convex, f increasing | correct |
| step bound Δy ≤ 2τ/a² (code `dyx`) | dε/dk = ((2N+1)² + (2N+1−2m)²)/(m²(2N+1−m)²) ≤ 4/m², m ≥ a+1 | correct (max observed (Δε)(N−k−1)² = 3.997) |
| c₀ ≤ 1 + 2\|α\|/a² | w₀ = α/(N(N+1)) is tiny | correct (c₀ ≈ 1 − 10⁻⁶) |
| entry into the window: N ≤ (2 + e1 + e1/√η)\|c\| | \|λ± − i\| ≤ u + (1−y) + \|r\|, \|r\|² ≤ 1 − y² + u² + 2uy, at y ∈ (y_A − Δy, y_A] | correct; audit E1 slack ≥ 0.034 |
| window step: N_{k+1} ≤ (1+√η+η)N_k | D_{k+1} = D_k + 2(μ−1)T_k, \|μ − 1\| ≤ δ + u = η/2 | correct for y ∈ (y_A, y_B); see N1 for the entry step |
| window length W | τ·Δ(1/y) = (P + ¼)·2/((2k+1)(2k+3)) + ½ ≥ 1 ⇒ #{y_k ∈ (y_A, y_B)} ≤ τ(1/y_A − 1/y_B) + 1 | correct (min observed 1.0069); the code uses +3 |
| conversion to phase III | \|c\| ≤ N(√η(1+e2) + e2)/(2r_B) via S_k − λ∓S_{k−1} = i^k(D + (1 + iλ∓)T_{k−1}) | correct; audit E2 slack ≥ 0.98 |
| phase-III growth | log\|λ₊\| ≤ acosh y + u/√(y²−1−u²); Σ acosh y_k ≤ Σ_{m ≥ a+1} acosh(τ/m) ≤ ∫_a^τ acosh(τ/m)dm = G_a(τ) | correct; audit P4 slack ≥ 4.5e-4; G_a is decreasing in a and increasing in τ |
| σ-term in phase III | Σu_k/r ≤ (σ/r_B)Σ_{m>a}1/m ≤ (σ/r_B)log(τ/a) | correct |

Nits. None is a gap, because each is absorbed by slack the code already has:

* **N1.** The first window step starts from y_K ≤ y_A (the conversion uses y ∈ (y_A − Δy, y_A], as `e1` shows). There \|μ−1\| can reach δ + Δy + u > η/2. The correct factor for that one step is 1 + η′ + η′/√η with η′ = η + 2Δy, which is about 0.03 nats worse. The code's +3 in W is two steps more than the proved +1, so this is covered.
* **N2.** In phase III the integrand f(y) = 1/r + qy/r² decreases in y, and \|dλ/dw\| is maximal at the lower end y_k of each segment. A rigorous sum therefore needs ∫_{y_B}^{y_max} f + f(y_B)·Δy. `lad` omits the boundary term, which is 0.1–0.27 nats for c = 0.7–0.5.
* **N3.** `lmax = y + √(y²−1) + 2u` is not an upper bound for \|λ₊\| near y_B. The rigorous form is y + u + √(y²−1+u²+2uy), about 3% larger at y_B ≈ 1.01. The phase-III adiabatic term also uses √(y²−1) instead of √(y²−1−u²), which is a relative error of at most 0.2%.
* **N4.** The final factor 2·lmax² over-covers the needed \|λ₊\|² + 1 by up to 0.4 nats. The phase-I term q·f_A·τ/a² is redundant, as are σ/(a r_B) and the +3.
* **N5.** `crossing.py`'s DV driver also minimises over c = 0.5. At c = 0.5 the uncorrected `cross` is 0.07–0.17 nats below my rigorous version, at τ ∈ [998.5, 1031]. The certificate path never uses c = 0.5: PP0 uses c ∈ {0.7, 1, 1.5}. Only DV printouts are affected.

`lx_common.bound_indep` repairs N1–N3 and removes the double slack. Against it, `rt_pp.cross` is larger by at least +0.066 nats at every c ∈ {0.7, 1, 1.5}, on a dense scan of the boxes PP0 actually uses (`logs/lx_phase3.log`).

## 3. Where it is used on [3.65, 3.75]

The call chain is `analytic2/rt2_contour.py` → `rt_cert.piece_right(L, p, Q, dress, Y)`.

* The edge is ν = Λ + iy, which gives β₀ = (p − ½) − y + iΛ. The chain bases are `PP(max(stop+k, 0), Λ−1, Λ_up, d0)` and `PP(·+1, …)`.
* PP calls PP0 at σ and σ + 1 for the helicity identity, so σ ≤ 5.5 (ψ₃: 3.5 + 1 + 1).
* For t_hi = 998.5 > 0.95·998, `thmE` returns None. PP0 then returns min over c ∈ {0.7, 1, 1.5} of cross(σ, max(Λ−1, 900), Λ_up, c).

Hypotheses on the cell:
* **τ.** τ ∈ [Λ−1, Λ_up] ⊂ [900, 1150]. ✓
* **σ.** 0 ≤ σ ≤ 5.5. The base covers Re β ∈ [σ_top − ¼, σ_top] ⊂ [0, σ_top], because k = ⌈y₁ − e₀⌉ ensures Re β ≥ 0. ✓
* **Twist and spin.** a ∈ [998, 4500], and indeed all a ≥ 998. Any n. ✓ (The b = 4 shift n+4, N+4 keeps a fixed.)
* **What the lemma does not assume.** It makes no assumption on L, d or N/a. Those (N/a ≥ 40, d²D ≥ 3.9a²) belong to the REC chain, not to Lemma X, and REVIEW_SLIVER2 S1 already proved them for L ≥ 3.65.

The bound actually used is log PP0 = 8.795 (σ = 0) up to 13.329 (σ = 5.5). Without the phase-III terms, the rigorous phase-II-end bound would be ≈ 1.2 nats smaller still.

## 4. Numerical tests (exact ladders)

**Precision.** Arb's interval wrapping makes ball radii grow like ∏(1 + \|ε_kα\|) ≈ e^{2\|α\| log(1/d)}, which is about 10⁴ bits at L ≈ 3.7.
* At 256–2000 bits Arb returns valid but useless balls, and the midpoint drifts to around 10²⁷¹.
* At 12 000 and 16 000 bits the balls are tight and agree with double-precision screening to 6 digits.

**Step audit (`lx_audit.py`, `logs/lx_audit.log`).**
* Scope: 150 trajectory × c combinations.
  * Cell-edge states: a ∈ {998, 998.49, 1000, 1010}, L ∈ {3.65, 3.75}, σ ∈ {0, 2.5, 5.5}, τ ∈ {Λ−1, Λ}.
  * Phase-III states up to τ = 1150, n ≤ 3·10⁵.
* Every intermediate inequality holds: P1/P2 ≥ 0 (to rounding), P3 ≥ 0.0039, E1 ≥ 0.034, W1 ≥ 0.013, W2 ≥ 0.67 steps, E2 ≥ 0.98 and P4 ≥ 4.5e-4.
* The window is entered (up to 10 steps) but never left on the edge.

**Cell scan (`lx_exact.py`, `logs/lx_exact.log`).**
* Grid: 21 twists a ∈ [998, 4500], including 998.1, 998.25, 998.4, 998.49 and 998.5; L ∈ {3.65, 3.70, 3.75}; 9 values of σ ∈ [0, 5.5]; τ ∈ {Λ−1, Λ−½, Λ}.
* Minimum of log PP0 − log\|Ŝ\| is **7.727**, at a = 998, n = 83372 (L = 3.75), σ = 0, τ = Λ.
  * Exact log\|Ŝ\| = 1.067389, from an Arb ball at 12 000 and 16 000 bits.
  * The margin grows with σ (9.17 at σ = 5.5) and with a (8.78 at a = 4500).
* The tightest 8 points are Arb-confirmed, with certified gaps ≥ 7.727–7.741.
* All-spin sweep at a = 998 (n from 1001 to 10⁶, L from 0.3 to 6.2): min gap 7.53.
* b = 4 (PP with helicity, vs exact s = 9): min gap 8.21.

The general lemma at τ up to 1150 (phase III active) has gaps ≥ 12 nats (audit FIN column).

## 5. Can base_bounds replace Lemma X on this edge?

Yes. `base_bounds.power_s1` covers the same boxes: τ ∈ [Λ−1, Λ_up], d ∈ [e^{−3.75}, e^{−3.65}], all a ≥ 998.
* It is spin-uniform and localised by d. Its turning-point parameter is Y ≈ 0.9946 < 1.
* For σ < ¼ it uses exact downward recurrence from σ ≥ ¼.

Comparison on the right-edge boxes (`lx_base_compare.py`):

| σ-box | Lemma X (PP0) | base_bounds | exact max |
|---|---|---|---|
| [0, ¼] | 9.01 | 4.84 | 1.21 |
| [¼, ½] | 9.22 | 3.95 | 1.35 |
| [2¼, 2½] | 10.90 | 5.12 | 2.48 |
| [5¼, 5½] | 13.33 | 6.84 | 4.16 |

Using it would require rewiring `piece_right`:
* call base_bounds on the actual Re β box, so σ_top ≥ ¼, or use its built-in downward steps;
* apply the base_bounds helicity (`base_bound`, b = 4).

The right edge's 0.32·F would then drop by a factor of roughly e^{4}–e^{6}.

## Files

* `code/review_lemx/lx_common.py`: exact Arb and double ladders, plus `bound_indep` (the independent rigorous re-derivation).
* `lx_bound.py`: auxiliary inequalities in exact rationals, monotonicity in a, and a comparison with `rt_pp.cross`.
* `lx_audit.py`: per-step audit.
* `lx_exact.py`: cell scan and Arb confirmation.
* `lx_phase3.py`: N2–N5 and the dense box scan.
* `lx_base_compare.py`: comparison with base_bounds.
* Note: an unrelated commit by another process (b3b67da) swept early versions of these scripts into git. I made no commits.
