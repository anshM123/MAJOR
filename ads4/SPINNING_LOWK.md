# Spinning tensor Regge moments at k = 2, 3: corner completion of the exact d = 3 partial waves and the end-to-end test

Status: written 16 Nov 2026 by a numerics/derivation agent. **Unrefereed.** Code and recorded outputs: `ads4/code/lowk/` (index §9;
`sh run_all.sh` regenerates every number; it needs the e2e block-tensor cache `ads4/code/e2e/data_modes/`). This note continues
`SPINNING_E2E.md` ("E2E"; conventions §1 there, notation unchanged) and addresses its §11 items 1 and 11 and the `**` flags at k = 2.
Nothing in `ads4/code/e2e/` was modified; its modules are imported.

Marks: **[derived]**, **[derived+verified]** (independent numerical check of the derived statement), **[numerical evidence]** (read off from
numbers, no derivation), **[inferred]**. Coefficients obtained by fits are labelled "fitted" and are never called derived.

## 0. Summary

1. **Method [derived+verified].** The slowly convergent corner `ρ_b → 1` of the reflected scalar seed blocks (`z̄ → −∞` at fixed `z`) is
   handled by re-organising the radial series into the lightcone monomial basis `ρ_s^i ρ_b^q` (an exact rearrangement of non-negative
   terms, also for all η-derivatives via Gegenbauer polynomials), and summing the `q`-tail at fixed `i` from its Darboux asymptotics,
   whose exponents follow from the local exponents of the Casimir equation at `z̄' = 1` and `z̄' = ∞`. Checks: the leading lightcone function
   equals the exact collinear block (₂F₁) to `≤ 8·10⁻⁸` down to `1 − ρ_b = 10⁻⁸`; derivative sums agree with an `N = 24000` radial series to
   `≤ 4·10⁻¹³` where the latter converges; Casimir residual `≤ 10⁻⁸` (relative to the largest term) down to `1 − ρ_b = 2·10⁻⁶` (§2).
2. **Currents (`l = 1`), k = 2 and 3: converged [numerical evidence for the ratios; the formula is the derived one of E2E §4.5,
   nothing fitted].** Errors in brackets are **model-scatter estimates** over three tail settings that pass the independent block checks
   (§3), not bounds; quadrature error not included (§4.3). At the largest masses (`m ≈ 100`; β ≈ 0.4 / 1.0, ν = 0):
   * k = 3: `R_pres = 0.997317(2) / 0.99698(1)`, `R_flip = 0.99478(3) / 0.99483(2)`; `m²(R−1)` converges like `c_2 + c_4/m²`
     (β ≈ 0.4, ν = 0: `−23.68, −26.39, −27.45, −27.74, −27.81` (pres), `−40.5, −49.1, −52.8, −53.9, −54.1` (flip), `m = 12 … 102`).
   * k = 2: `R_pres = 0.9984(2) / 0.9986(3)`, `R_flip = 0.9976(40) / 0.9976(16)`; `m²(R−1)` is converged within errors from `m ≈ 40`
     (β ≈ 0.4, ν = 0: `−15.36, −16.49, −16.92, −17.01(10), −16.6(21)` (pres), `−27.7, −32.5, −34.2(22), −33(7), −25(41)` (flip)).
   * all three β (0.4, 1, 2) and ν = 0, 1, 3: the full tables are §5; at β ≈ 2 (only `m ≤ 70`) `m²(R−1)` is still drifting for ν = 3.
3. **Stress tensors (`l = 2`) [partial].** Helicity-**preserved** component converged at k = 3 and consistent at k = 2: k = 3, β ≈ 0.4:
   `R_pres = 0.53930(1), 0.84003(1), 0.95323(4), 0.9838(1), 0.9919(2)` (`m = 12 … 102`; `m²(R−1) = −71.0, −79.9, −83.4, −84.2(4), −84(2)`);
   k = 2: `0.5725(2) … 0.9937(24)`, `m²(R−1) = −65.9, −71.2, −73.0(12), −71(7), −65(24)`; β ≈ 1, m ≈ 101: `R_pres = 0.9961(2)` (k = 3),
   `0.9982(20)` (k = 2). Helicity-**flip** component: k = 3 determined at β ≈ 1, 2 (`R_flip = 0.999(15)` at m ≈ 101, β ≈ 1; `0.978(2)` at
   m ≈ 70, β ≈ 2) but **not** at β ≈ 0.4 (`1.27 ± 0.49` at m ≈ 102); at k = 2 it is **not determined** at β ≈ 0.4 (errors > 1) and only to
   10–30 % at β ≈ 1, 2. Reason: for `l = 2` the deepest corner and the tiny-t region carry up to 10 % of the `(2,2)` mode at k = 2 and cannot be
   evaluated in double precision (cancellation between the seeds, §8 item 1); they are estimated with 50 % error, which dominates.
4. **No `O(1/m)` term at k = 2, 3 [numerical evidence; the argument of E2E §6.1 now holds without Assumption T, see 5].** Fitted `c_1`
   (ansatz `c_1/m + c_2/m² + c_4/m⁴`, `m ≥ 20`), `l = 1`, β ≤ 1: `|c_1| ≤ 0.009` (k = 2, pres), `≤ 0.08 ± 0.25` (k = 2, flip), `≤ 0.010` (k = 3),
   all consistent with zero within their statistical errors at k = 2.
5. **Non-analytic terms / Assumption T [derived (sketch) + numerical evidence].** Through the exact identity (★) of E2E §3.1 the moment
   is a 6-dimensional integral of the elementary shadow integrand; in heavy variables its expansion in `1/m` is uniform in all directions
   (every Cartesian component of `X = mx` is bounded by `e^{β/2} ζ`, `ζ = −P̂·X`, which the exponential `e^{−2ζ}` controls), so the moment has
   an asymptotic expansion in powers of `1/m` with no non-analytic terms — neither `m^{−k−1/2}` (breakdown scale `u ~ 1/m`) nor
   `m^{−2k−1}` (`u ~ 1/m²`) — for every k for which (★) converges (for `Δ_φ = 1.75`, `l ≤ 2`: all k > −2). The non-uniformity at `u → 0` that motivated
   Assumption T is a property of the `(r, η)` representation only. Numerically: the pointwise `1/m²` coefficient of the shadow integrand
   stays bounded (it saturates) when the hyperbolic distances of `x̂, ŷ` from `p̂, p̄̂` grow to 8 (`d(x̂, ŷ)` up to 12.6), `l = 1, 2` (§6);
   fits of the completed k = 2 ratios with an added `m^{−5/2}` term give a coefficient consistent with zero (`−1.3 ± 4.8` against
   `c_2 ≈ −17.0`, `l = 1`, β ≈ 0.4, ν = 0, pres), but only `|c| ≲ 10` (≈ 0.6|c₂|, 2σ) is excluded — an O(1) coefficient is allowed by the data;
   `m^{−5}` at k = 2 and `m^{−7/2}` at k = 3 cannot be separated from the analytic truncation at all. So: no failure at `1/m²`
   [numerical evidence]; beyond that the absence of non-analytic terms rests on the sketch; no rigorous remainder estimate is written.
6. **Fitted `O(1/m²)` coefficients** (`l = 1`, k = 2 / 3, β ≈ 0.4, ν = 0; mid-range of the acceptable fits ± half the range):
   `c_2^pres = −17.03(11) / −27.87(5)`, `c_2^flip = −31(8) / −54.29(26)`. Table §7.
7. **What is not done:** see §8 (most importantly: the `l = 2` deepest corner at k = 2, 3 is limited by double precision; scalar-test caveat).

**Physical reach (applicability).** *Established* [numerical evidence for the ratios; derived leading formula, nothing fitted], for **one**
spinning Regge pair in `⟨OφOφ⟩` (parity-even exchange `(Δ, J)`, unique conserved `⟨OφO⟩` structure), `Δ_φ = 1.75`, ν ≤ 3, β ∈ {0.4, 1, 2},
`m ≤ 100` (`≤ 70` at β ≈ 2): the leading heavy formula of E2E §4.5 holds at k = 2 and 3, with `m²(R−1)` converging, for **currents in both
helicity components** and for the **stress-tensor helicity-preserved component**; no `O(1/m)` term and no `m^{−5/2}` term are detected (the latter
only to `|c| ≲ 10`), and by the sketch of §6 Assumption T is not needed for `⟨OφOφ⟩`. All quoted errors are model-scatter estimates, not bounds.
*Not established:* the **stress-tensor helicity-flip** component at k = 2 (and at k = 3 for β ≈ 0.4) — the component the `B₃`/flip sum rules
need; `⟨TTTT⟩` and configurations with both pairs spinning; the mean-field / positivity normalisation at k = 2, 3 (cf. `SPINNING_MFT.md`);
uniformity in `ν ~ m` or `β ≲ 1/m`; other `Δ_φ`; a rigorous remainder bound for the `1/m` expansion (§6 is a sketch). Accordingly, the STATUS
row for the heavy-block formula may move "k = 2, 3" from "not established" to "established for one pair, `l = 1` both helicities and `l = 2`
preserved, within these ranges".

## 1. The problem and the choice of method

E2E §3.2: in the variables `(u, t)` (`√ρ = s = t/(2(α+ᾱu))`, `√ρ̄ = su`, `u = e^{−ξ}`) the moment integrand behaves like `u^{k−1/2}` for
`u → 0`, where the reflected scalar seed blocks are needed at `ρ_b = ((1−su)/(1+su))² → 1` (`z̄' → 1`), with radial convergence `e^{−4Nsu}`.
With `N ≤ 10⁴` the e2e runs drop the grid points with `u ≲ 30α/(Nt)`; at k = 2 the dropped region carries up to `10⁻¹` of a mode.

Options considered (task list): (i) an expansion of the d = 3 block around `z̄' = 1` needs the connection coefficients of its singular term,
not known in closed form in d = 3; (iii) the d = 3 → d = 2 reduction converges like `ρ_s^{n/2}` and needs its coefficients by an O(N²J²) matching;
(iv) per-point acceleration of the level sums mixes the `ρ_s`- and `ρ_b`-structure. **Chosen:** a hybrid of (i) and (iv) that needs no new
connection data — the lightcone reorganisation makes each `q`-series a function of one variable with *known* singularities and exponents,
whose tail follows from the Darboux expansion with coefficients fitted on `q ∈ [Q/2, Q]` and summed in closed form.

## 2. Scalar seed blocks near `ρ_b → 1` (`lcblock.py`, `lcconv.c`)

### 2.1 Lightcone reorganisation, all η-derivatives [derived]

With `B_{n,j} ≥ 0` the radial coefficients (E2E §2.1; reflected block `g^{c1,c2}`, `c1 + c2 > 0`), `ρ_s = re^{−χ}`, `ρ_b = re^{χ}`, and
`P_j^{(l)} = (2l−1)!! C^{(l+½)}_{j−l}`, `C^λ_m(cos θ) = Σ_k a^λ_k a^λ_{m−k} e^{i(m−2k)θ}`, `a^λ_k = (λ)_k/k! > 0`:

`S[k,l] := Σ_{n,j} B_{n,j} (Δ+n)^k r^n P_j^{(l)}(cosh χ) = (2l−1)!! r^{l−J} Σ_{i,q≥0} (Δ−J+l+i+q)^k C^{(l)}_{i,q} ρ_s^i ρ_b^q`,

`C^{(l)}_{i,q} = Σ_{k'} a^{l+½}_{k'} a^{l+½}_{|q−i|+k'} B_{i+q−J+l, |q−i|+2k'+l}` — every term non-negative, so no cancellation (an earlier version that
took η-derivatives as `(sinh χ)^{−1}(θ_b − θ_s)` lost up to 9 digits near `χ → 0` and was replaced). Check: all fifteen `S[k,l]`, `k + l ≤ 4`,
reproduced to `2·10⁻¹⁴` at a generic point; the leading lightcone function `Σ_q C^{(0)}_{0,q}x^q` equals `4^{−(Δ+J)/2}x^{−(Δ−J)/2}k_{Δ+J}(4x/(1+x)²)`,
`k_β(z) = z^{β/2}₂F₁(β/2+c1, β/2+c2; β; z)`.

### 2.2 Tails in `q` [derived; one input inferred]

At fixed `i`, `H^{(l)}_i(x) = Σ_q C^{(l)}_{i,q} x^q` has on `|x| ≤ 1` only the singular points `x = 1` (`z̄' = 1`) and `x = −1` (`z̄' = ∞`).
The Casimir operator (Dolan–Osborn form) has local exponents `{0, −σ}`, `σ = c1 + c2`, at `z̄ = 1` and `{−c1, −c2}` at `z̄ = ∞` (the `(d−2)` term
is of lower order in both indicial equations). With `1 − z̄ = ((1−x)/(1+x))²` and one extra power per η-derivative:
`H^{(l)}_i(x) = A(x) + (1−x)^{−2σ−l}B(x)` near 1, `(1+x)^{2c_{1,2}−l}`-type at −1 (logarithmic if `c1 − c2` is an integer — `c1 = c2`, and also the seed `(−0.375, 0.625)` with `c1 − c2 = 1`, which the first version missed (referee m2; effect ≤ `5.5·10⁻⁶` on a flip component); log columns are now in the fit, their alternating tails are not summed); [inferred: that the Frobenius
structure survives the expansion in `ρ_s` term by term — tested by §2.3]. For `Δ_φ = 1.75` every `2σ` is a half-integer: no logarithms at 1.
Darboux: `(1−x)^{γ} = Σ (−γ)_q/q! x^q`, hence for `t = 0…4−l`

`q^t C^{(l)}_{i,q} ≃ Σ_{m<M} s_m (2σ+l+t−m)_q/q! + (−1)^q Σ_{c∈{c1,c2}, m<3} a_{c,m} (l−2c+t−m)_q/q!`.

Coefficients: one least-squares problem per `(l, t)` on `q ∈ [Q/2, Q]` (every third `q`) with `I+1` right-hand sides. The series is summed exactly to
`Q`; the tails `T(a,x) = Σ_{q>Q}(a)_q/q! x^q` in closed form: `(1−x)^{−a}I_x(Q+1, a)` for `a > 0` (regularised incomplete beta; the partial sums
are negative-binomial distribution functions), an upward recursion in `a` for `a < 0`, a Pfaff-transformed ₂F₁ at `x < 0` (all checked against
mpmath ₂F₁: `≤ 10⁻¹¹`). The `ρ_s`-sum is truncated at `i = I`; a point is used only if its last term is `< 10⁻¹³` of the sum. Settings:
**A** `(N, I, Q, M) = (10⁴, 1500, 8000, 8)` (used), **B** `(8000, 1200, 6000, 6)` (same `B_{n,j}` — they do not depend on the truncation — but a
different depth, window and number of asymptotic terms; error estimate).

### 2.3 Validation of the block representation (`c1_block_checks.py`, `out_c1_block_checks.txt`) [derived+verified]

| check | blocks | result |
|---|---|---|
| leading lightcone function vs mpmath ₂F₁ (independent) | `(Δ,J,c1,c2) = (14.17,2,.625,.625), (105.36,20,.625,−.375), (115.11,52,1.625,1.625), (109.83,82,.375,.375)` | `1−x = 10⁻²: ≤ 10⁻¹⁴`; `10⁻⁴: ≤ 2.5·10⁻⁸`; `10⁻⁶, 10⁻⁸: ≤ 8.3·10⁻⁸` |
| `S[k,l]`, `k+l ≤ 4`, vs radial series with `N = 24000` where it converges (`1−ρ_b ≥ 4·10⁻³`; tests the extrapolation `q = 8000 → 24000`) | `(14.17,2,.625,−.375), (24.22,10,.625,.625), (24.22,10,1.625,−.375)` | `≤ 4.1·10⁻¹³` |
| Casimir residual / largest term, deep corner, A and B (with `1 − z = ((1−ρ)/(1+ρ))²`; the first version computed `1 − z` by subtraction and its "`≤ 10⁻⁴`" at `1−ρ_b ≈ 3·10⁻⁶` was the check's own roundoff, referee m1) | `(115.11,52,.625,−.375), (105.36,20,.375,.375), (14.17,2,.625,−.375)` | A: `≤ 1.0·10⁻⁸` for `1−ρ_b` down to `2.4·10⁻⁶`, `≤ 2·10⁻¹²` for `1−ρ_b ≥ 2·10⁻³`; B: `≤ 3.3·10⁻⁹` |
| `|S_A/S_B − 1|`, same points | same | `≤ 9·10⁻⁸` (deepest), `≤ 10⁻¹¹` for `1−ρ_b ≥ 10⁻³` |

The deep-corner A/B difference (`4–9·10⁻⁸`) is of the size of the deviation of A from the exact collinear block (`3–8·10⁻⁸`). It is **not** a
bound: A and B share the Darboux exponents, so a wrong exponent would leave `|A − B|` (and the fit residual) unchanged — the referee's mutation
`2σ+l+t−m → +0.25` moves k = 2 by O(1) with `|A − B|` unchanged (m5: the fit-residual diagnostic is not a validation and is not cited as one).
What protects the exponents are the independent checks of this table (₂F₁, `N = 24000`, Casimir) and the referee's η-derivative check of the
`l ≥ 1` tails (finite differences in χ of `S[k,l]` vs `S[k,l+1]`, `1−ρ_b` down to `10⁻⁷`: ≤ `3·10⁻⁷`). [The radial reference with `N = 24000` costs `2.3 GB`; the LC representation of one
seed with setting A takes ~40 s.]

## 3. Corner completion of the moment grids (`corner.py`) and error budget

**Grid.** The e2e grid (72 × 62 Gauss nodes in `(u, t)`) under-resolves the corner at k = 2: below `u = 0.01` it has 8 nodes, too few for the
log-oscillating kernel `P_{−½+iν}(cosh ξ)`, `ξ = −log u`, at ν = 3 (§4.3). The results use the **corner-refined grid** (`refined_base`): the two
u-panels below 0.03 (18 nodes) are replaced by six panels with 46 nodes from `u = 3·10⁻⁶`; above `u = 0.03` the nodes and the e2e radial-series
tensors are kept. All corner nodes, and the e2e points that were dropped, are evaluated with A and B through the *unchanged* spinning machinery
(`spinblock.SpinBlock.eval`: differential operators, `(P,Q)` jets, generated tensor code, Fourier modes); 300 random converged e2e points with
`u < 0.25` are re-evaluated with A as a check.

**Settings.** Three tail settings: A (used), B, and C = `(N, I, Q, M) = (10⁴, 1500, 3000, 8)` (exact sum only to `q = 3000`; added after the
referee found that this setting, as good as A by the corrected Casimir residual, moves k = 2 by 1.3–1.8× and the k = 3 flip by 3× the first
version's errors). The per-point error is `3·max(|A−B|, |A−C|)`. **These errors are model-scatter estimates, not bounds**: they measure the
dependence on the tail setting among settings that pass the independent checks; settings that fail the Casimir check (e.g. `M = 12` or a
window `[6000, 8000]`) move results by up to 200× more and are excluded by the gate below.

**Per point.** With `δ = max_modes max(|A−B|, |A−C|) / max_modes|A|`: the point is accepted if `δ < 10⁻²`, the Casimir residual of every seed
(setting A, corrected formula) is `< 10⁻⁷`, **and** it is connected, within its t-column, to the
radial-series region through accepted points whose size changes by less than `u^{±8}` between neighbours (isolated points inside a rejected
region are where A and B are both corrupted). Its error is `3·max(|A−B|, |A−C|)` (summed in absolute value). Rejected points (deepest corner: the conserved
spinning combination cancels between the seeds by up to ~`(1−ρ_b)^{−(2σ+4)}`, beyond double precision) and points where the lightcone sum has not
converged (tiny `t`, where also the radial series fails) are **estimated** by a power law through the two nearest accepted nodes of their u-row
(tiny t) or t-column (deep corner); the estimate is added to the moment and **50 % of it** is counted as error. [The 50 % is a choice, not derived;
the referee's leave-out test gives estimate/true = 1.00–1.19 for one-node and 0.44–1.49 for two-node extrapolation, so it is adequate for the
former and marginal for the latter. The estimated share is listed for every entry (`est` in `out_t1_lowk_e2e_l*.txt`).]
**Not in the budget** (referee m4): the quadrature error of the grid (§4.3: `≤ 2·10⁻⁵` at k = 2, `≤ 3·10⁻⁸` at k = 3 for the refined grid at m ≈ 70–100;
it can exceed the quoted error for m ≲ 40) and, for `l = 2`, the LC-vs-radial deviations of the validation below (comparable to the quoted
error for the k = 3 flip at m ≈ 42, β ≈ 0.4). Radial-series points: e2e criterion, `≤ 10⁻⁹`.

**Gate statistics.** Per case 47–320 corner points fail the Casimir gate (mostly at tiny t, where also the lightcone sum is marginal) and are
estimated like the other rejected points; for `l = 1` the k = 2, 3 ratios moved relative to the first version by at most 2.2× their new error (k = 3, β ≈ 2, m ≈ 10), in 10 of 84 entries by more than 1×.

**Per ratio.** The relative errors of the modes `(l,−l)`, `(l,l)` are propagated through the 2×2 helicity solve (E2E §5.2) to `R_pres`, `R_flip`.
Validation of the re-evaluated converged points (unchanged by the three-setting rebuild): `max dev ≤ 5·10⁻⁹` for `l = 1` (all cases); for `l = 2` it ranges from `2·10⁻⁶` (β ≈ 1, 2, m ≥ 40) up to `0.2–0.3` (m ≈ 10–20) and `1.5·10²` at the single worst point (m ≈ 12, β ≈ 0.4) — the
`l = 2` tensor at small `s` amplifies `10⁻¹²` errors of the seed jets by up to `10¹⁰`. Its impact on the moments, measured by evaluating **all**
converged points with `u < 0.25` with LC and weighting the deviations (`v2_l2_validation.py`, `out_v2_l2_validation.txt`): `l = 2`, m ≈ 42:
`1.1·10⁻⁵ / 8.3·10⁻⁴` (k = 2, modes `(2,−2)/(2,2)`), `1.1·10⁻⁶ / 2.2·10⁻⁴` (k = 3); m ≈ 102: `2·10⁻⁷ / 1.2·10⁻⁵` (k = 2); `l = 1`: `≤ 3·10⁻¹¹`.
These numbers are not in the `l = 2` budget (the corner points have their own A/B test); they are smaller than the quoted `l = 2` errors
except for the flip component at k = 3, m ≈ 42, β ≈ 0.4, where they are comparable.

## 4. Validation of the completed moments

### 4.1 Scalar partial waves (`v1_scalar.py`, `out_v1_scalar.txt`)

`⟨φ_1φ_2φ_1φ_2⟩`, s-channel, k = 2, 3, compared with the derived scalar formula of E2E §4.1 (nothing fitted); errors as in §3 (three tail settings; the numbers below are unchanged from the two-setting run at the printed precision).

* `(D_1, D_2) = (3, 1.75)` (corner weight `u⁴` at k = 2; the radial series alone nearly converges): on the e2e grid the completed values equal
  the radial-only ones (E2E `t3` method) to `≤ 5·10⁻⁷` (e.g. m ≈ 72, β ≈ 0.4, k = 2, ν = 0: `0.999794878` vs `0.999794844`) [consistency]; the
  corner-refined grid changes them by up to `5·10⁻⁵` (m ≈ 100, β ≈ 1, k = 2, ν = 3) — the quadrature error of the e2e grid at k = 2 (§4.3).
* `(D_1, D_2) = (3.75, 0.5)` (corner weight `u^{3/2}` at k = 2, `u^{5/2}` at k = 3 — the spinning corner weights; a harsher case, since the
  tiny-t region carries far more weight than for the spinning moments). Radial series alone (as in E2E `t3`) vs completed, β ≈ 0.4, ν = 0:

| m | k = 2: radial only | completed (refined grid) | `m²(R−1)` | k = 3: completed | `m²(R−1)` |
|---|---|---|---|---|---|
| 12.4 | 1.098540 | 1.099325(5) | 15.31 | 1.0718088 | 11.07 |
| 22.4 | 1.028018 | 1.028806(5) | 14.39 | 1.0199205 | 9.95 |
| 42.2 | 1.006130 | 1.00795(20) | 14.16(36) | 1.0053633(10) | 9.56 |
| 72.0 | 0.999040 | 1.0034(20) | 17(10) | 1.001832(30) | 9.50(16) |
| 101.8 | 0.995501 | 1.0043(70) | 45(73) | 1.00097(20) | 10.1(21) |

  k = 3: `m²(R−1)` converges (`≈ 9.5`) through m ≈ 102 within errors; k = 2: converged within `±0.4` up to m ≈ 42 and consistent (with large
  errors, from the tiny-t estimate) beyond; the radial series alone is off by up to `10⁻²` (`m²(R−1) = 10.9, −5.0, −46.6` at m = 42, 72, 102).
  β ≈ 1 behaves the same (`out_v1_scalar.txt`). The e2e grid vs the refined grid differ by up to `5·10⁻³` at k = 2, ν = 3 for this scalar (log-oscillating
  kernel), which is why the refined grid is used. [derived formula verified for k = 3 to m ≈ 100 and for k = 2 to m ≈ 40 on a corner-heavy
  scalar; numerical evidence]

### 4.2 k = 4, 6 against the e2e numbers [numerical evidence]

`l = 1`: where E2E's numbers were flagged converged (β ≤ 1, k = 6, all m; k = 4, m ≤ 42) the completed values agree with them to
`≤ 3.5·10⁻⁵` (k = 6, flip, β ≈ 1, m ≈ 100; typically `10⁻⁶`); the larger differences sit exactly at E2E's `*`/`**` entries
(k = 4, `m ≳ 70`: up to `5.9·10⁻³`; β ≈ 2, m ≈ 70: `1.7·10⁻⁴` at k = 6, `8.4·10⁻³` at k = 4) and are the dropped corner that E2E's `tail` estimated
(E2E: k = 4 tail `2·10⁻⁴ … 2·10⁻³` at the largest m). `l = 2`: k = 6 within `≤ 6·10⁻⁴` of E2E (largest at E2E's `*` entry, flip, m ≈ 102, β ≈ 0.4), the preserved component within `≤ 5·10⁻⁶`
(β ≤ 1); k = 4 flip differs by up to `0.14` at m ≈ 102, β ≈ 0.4 (E2E `**`; new error `2·10⁻²`, the `l = 2` flip problem of §8 item 1).

### 4.3 Quadrature (`q1_grid.py`, `out_q1_grid.txt`)

(These grid comparisons were made with the first-version two-setting acceptance; the grid differences do not depend on the error budget.)

Same `(Δ, J)` on the e2e grid (A), the corner-refined grid (B, used) and, for two cases, a grid with ~1.6× the corner nodes (C):

| quantity | `l = 1` (m ≈ 42, 102 at β ≈ 0.4; 101 at β ≈ 1; 70 at β ≈ 2) | `l = 2` (m ≈ 42, 102, β ≈ 0.4) |
|---|---|---|
| `|A − B|`, k = 2 | pres `≤ 1.2·10⁻⁴` (ν = 3), `≤ 1.4·10⁻⁵` (ν ≤ 1); flip `≤ 1.5·10⁻⁴` | pres `≤ 2.9·10⁻⁴`; flip `0.03–0.3` (not determined) |
| `|A − B|`, k = 3 | `≤ 1.4·10⁻⁶` | pres `≤ 9·10⁻⁶`; flip `≤ 10⁻²` |
| `|B − C|` (m ≈ 102, β ≈ 0.4; m ≈ 70, β ≈ 2) | k = 2: `≤ 2.9·10⁻⁶` pres, `≤ 2.2·10⁻⁵` flip; k = 3: `≤ 3·10⁻⁸` | — |

So the e2e grid's quadrature error at k = 2 is `~10⁻⁴` (it under-resolves the ν = 3 kernel near `u = 0`), and the refined grid's is `≤ 2·10⁻⁵`,
below the quoted numerical errors at the masses where they matter (`m ≳ 70`); at `m ≲ 40` the quadrature error of B (`~10⁻⁶`) can exceed the
quoted error, which does not include it.

## 5. Results (`out_t1_lowk_e2e_l{1,2}.txt`, `out_t3_tables.txt`)

Measured = tensor Regge moments of the exact partial wave (corner-completed); predicted = E2E §4.5 with the exact `b^{(1−Δ_φ)}`; nothing fitted.
Errors in units of the last digit (the full numerical budget of §3; quadrature not included, see §4.3).

### 5.1 Currents, `l = 1`

#### l = 1, k = 2: R_pres / R_flip (measured/predicted; error in units of the last digit)
| β ≈ 0.4, (Δ,J) → | (14.17, 2) | (24.30, 4) | (44.57, 8) | (74.96, 14) | (105.36, 20) |
|---|---|---|---|---|---|
| m | 12.4 | 22.4 | 42.2 | 72.0 | 101.8 |
| ν = 0 | 0.900381(1) / 0.8207(2) | 0.966983(3) / 0.9349(4) | 0.99051(1) / 0.9808(7) | 0.99672(2) / 0.9936(13) | 0.9984(2) / 0.9976(40) |
| ν = 1 | 0.917689(1) / 0.83088(2) | 0.973239(1) / 0.93903(5) | 0.992362(2) / 0.9818(2) | 0.99737(1) / 0.9936(4) | 0.9987(1) / 0.9972(24) |
| ν = 3 | 1.08402(1) / 0.9197(1) | 1.02852(2) / 0.9743(2) | 1.00823(4) / 0.9927(3) | 1.0028(2) / 0.9976(7) | 1.0012(14) / 0.9995(52) |
| m²(R−1), ν = 0 | -15.36 / -27.65 | -16.49 / -32.54 | -16.92 / -34.22 | -17.01 / -33.06 | -16.61 / -24.64 |
| m²(R−1), ν = 1 | -12.69 / -26.08 | -13.37 / -30.46 | -13.61 / -32.41 | -13.67 / -33.27 | -13.46 / -29.47 |
| m²(R−1), ν = 3 | 12.95 / -12.37 | 14.25 / -12.86 | 14.67 / -12.93 | 14.57 / -12.34 | 12.10 / -4.70 |

| β ≈ 1.0, (Δ,J) → | (13.40, 5) | (24.22, 10) | (45.86, 20) | (80.48, 36) | (115.11, 52) |
|---|---|---|---|---|---|
| m | 10.6 | 20.1 | 39.3 | 70.0 | 100.7 |
| ν = 0 | 0.888529(1) / 0.82325(1) | 0.963971(2) / 0.93646(3) | 0.989996(4) / 0.9816(1) | 0.99681(3) / 0.9941(3) | 0.9986(3) / 0.9976(16) |
| ν = 1 | 0.934508(1) / 0.844592(1) | 0.980453(1) / 0.945288(2) | 0.994749(2) / 0.98434(1) | 0.99834(3) / 0.9950(2) | 0.9993(3) / 0.9980(12) |
| ν = 3 | 0.77479(1) / 1.05608(1) | 0.94763(3) / 1.02241(2) | 0.9874(1) / 1.0066(1) | 0.9961(9) / 1.0021(5) | 0.9987(80) / 1.0011(45) |
| m²(R−1), ν = 0 | -12.42 / -19.69 | -14.63 / -25.80 | -15.48 / -28.41 | -15.65 / -28.80 | -14.48 / -24.33 |
| m²(R−1), ν = 1 | -7.30 / -17.31 | -7.94 / -22.21 | -8.13 / -24.24 | -8.14 / -24.60 | -7.58 / -20.56 |
| m²(R−1), ν = 3 | -25.09 / 6.25 | -21.26 / 9.10 | -19.44 / 10.29 | -18.90 / 10.45 | -12.83 / 10.81 |

| β ≈ 2.0, (Δ,J) → | (16.60, 11) | (32.36, 23) | (63.87, 47) | (109.83, 82) |
|---|---|---|---|---|
| m | 9.8 | 20.0 | 40.4 | 70.2 |
| ν = 0 | 0.800438(1) / 0.772104(4) | 0.929850(2) / 0.90632(1) | 0.98034(1) / 0.97227(3) | 0.9933(1) / 0.9904(3) |
| ν = 1 | 1.056806(1) / 0.829845(1) | 1.035687(2) / 0.934405(2) | 1.01168(1) / 0.98137(1) | 1.0042(2) / 0.9937(2) |
| ν = 3 | 0.16650(2) / 0.87807(1) | 0.7831(1) / 0.98778(3) | 0.9519(3) / 1.0020(2) | 0.9850(38) / 1.0016(24) |
| m²(R−1), ν = 0 | -19.11 / -21.82 | -28.05 / -37.46 | -32.12 / -45.29 | -33.13 / -47.07 |
| m²(R−1), ν = 1 | 5.44 / -16.29 | 14.27 / -26.23 | 19.09 / -30.43 | 20.58 / -30.94 |
| m²(R−1), ν = 3 | -79.81 / -11.68 | -86.71 / -4.89 | -78.52 / 3.27 | -73.69 / 7.95 |

#### l = 1, k = 3: R_pres / R_flip (measured/predicted; error in units of the last digit)
| β ≈ 0.4, (Δ,J) → | (14.17, 2) | (24.30, 4) | (44.57, 8) | (74.96, 14) | (105.36, 20) |
|---|---|---|---|---|---|
| m | 12.4 | 22.4 | 42.2 | 72.0 | 101.8 |
| ν = 0 | 0.846393(1) / 0.737319(1) | 0.947179(1) / 0.901817(1) | 0.984599(1) / 0.970380(2) | 0.994651(1) / 0.989615(5) | 0.997317(2) / 0.99478(3) |
| ν = 1 | 0.866025(1) / 0.749008(1) | 0.954769(1) / 0.907049(1) | 0.986909(1) / 0.972074(1) | 0.995463(1) / 0.990220(1) | 0.997725(1) / 0.99509(2) |
| ν = 3 | 1.057553(1) / 0.851083(1) | 1.022568(1) / 0.950099(1) | 1.006833(1) / 0.985662(1) | 1.002395(1) / 0.995043(2) | 1.00120(1) / 0.99752(3) |
| m²(R−1), ν = 0 | -23.68 / -40.50 | -26.39 / -49.05 | -27.45 / -52.79 | -27.74 / -53.86 | -27.81 / -54.09 |
| m²(R−1), ν = 1 | -20.66 / -38.70 | -22.60 / -46.43 | -23.33 / -49.77 | -23.53 / -50.73 | -23.58 / -50.94 |
| m²(R−1), ν = 3 | 8.87 / -22.96 | 11.27 / -24.93 | 12.18 / -25.56 | 12.42 / -25.71 | 12.47 / -25.71 |

| β ≈ 1.0, (Δ,J) → | (13.40, 5) | (24.22, 10) | (45.86, 20) | (80.48, 36) | (115.11, 52) |
|---|---|---|---|---|---|
| m | 10.6 | 20.1 | 39.3 | 70.0 | 100.7 |
| ν = 0 | 0.799456(1) / 0.699700(1) | 0.931838(1) / 0.888731(1) | 0.980694(1) / 0.967348(1) | 0.993784(1) / 0.989379(1) | 0.99698(1) / 0.99483(2) |
| ν = 1 | 0.856115(1) / 0.724263(1) | 0.954378(1) / 0.900062(1) | 0.987459(1) / 0.970977(1) | 0.995997(1) / 0.990589(1) | 0.998061(4) / 0.99542(1) |
| ν = 3 | 0.582111(1) / 0.982455(1) | 0.884272(1) / 1.003612(1) | 0.970311(1) / 1.002109(1) | 0.990720(3) / 1.000775(1) | 0.9955(1) / 1.00039(3) |
| m²(R−1), ν = 0 | -22.34 / -33.45 | -27.67 / -45.18 | -29.88 / -50.53 | -30.50 / -52.11 | -30.63 / -52.46 |
| m²(R−1), ν = 1 | -16.03 / -30.72 | -18.52 / -40.58 | -19.41 / -44.92 | -19.64 / -46.17 | -19.69 / -46.44 |
| m²(R−1), ν = 3 | -46.55 / -1.95 | -46.99 / 1.47 | -45.95 / 3.26 | -45.53 / 3.80 | -45.34 / 3.95 |

| β ≈ 2.0, (Δ,J) → | (16.60, 11) | (32.36, 23) | (63.87, 47) | (109.83, 82) |
|---|---|---|---|---|
| m | 9.8 | 20.0 | 40.4 | 70.2 |
| ν = 0 | 0.645885(1) / 0.581859(1) | 0.865664(1) / 0.827279(2) | 0.960846(1) / 0.947692(1) | 0.986403(1) / 0.981627(1) |
| ν = 1 | 0.968816(1) / 0.650313(1) | 1.019637(1) / 0.865391(1) | 1.009717(1) / 0.960872(1) | 1.003765(1) / 0.986435(1) |
| ν = 3 | -0.324554(1) / 0.569521(1) | 0.56827(1) / 0.896693(4) | 0.891677(1) / 0.980039(1) | 0.96433(1) / 0.99410(1) |
| m²(R−1), ν = 0 | -33.91 / -40.04 | -53.72 / -69.06 | -63.96 / -85.45 | -67.01 / -90.54 |
| m²(R−1), ν = 1 | -2.99 / -33.49 | 7.85 / -53.82 | 15.87 / -63.92 | 18.56 / -66.85 |
| m²(R−1), ν = 3 | -126.84 / -41.22 | -172.63 / -41.31 | -176.96 / -32.61 | -175.81 / -29.10 |

### 5.2 Stress tensors, `l = 2`

#### l = 2, k = 2: R_pres / R_flip (measured/predicted; error in units of the last digit)
| β ≈ 0.4, (Δ,J) → | (14.17, 2) | (24.30, 4) | (44.57, 8) | (74.96, 14) | (105.36, 20) |
|---|---|---|---|---|---|
| m | 12.4 | 22.4 | 42.2 | 72.0 | 101.8 |
| ν = 0 | 0.5725(2) / -0.0677(8487) | 0.8575(4) / 1.7637(21862) | 0.9590(7) / 2.7782(40462) | 0.9862(13) / 5.4661(69776) | 0.9937(24) / 7.8756(104207) |
| ν = 1 | 0.5938(1) / 0.3850(2059) | 0.8678(2) / 1.2143(5318) | 0.9621(4) / 1.7504(12224) | 0.9864(8) / 3.3997(28317) | 0.9926(13) / 5.3827(50364) |
| ν = 3 | 0.8230(4) / 0.3919(1172) | 0.9687(11) / 1.2622(3423) | 0.9910(23) / 1.7922(7515) | 0.9937(54) / 1.2698(10932) | 0.9939(102) / 0.3862(16331) |
| m²(R−1), ν = 0 | -65.91 / -164.62 | -71.19 / 381.52 | -73.01 / 3169.39 | -71.44 / 23164.54 | -65.45 / 71280.97 |
| m²(R−1), ν = 1 | -62.63 / -94.83 | -66.06 / 107.05 | -67.59 / 1337.49 | -70.34 / 12446.49 | -76.85 / 45437.11 |
| m²(R−1), ν = 3 | -27.30 / -93.76 | -15.61 / 130.99 | -16.04 / 1411.94 | -32.89 / 1399.26 | -62.91 / -6363.65 |

| β ≈ 1.0, (Δ,J) → | (13.40, 5) | (24.22, 10) | (45.86, 20) | (80.48, 36) | (115.11, 52) |
|---|---|---|---|---|---|
| m | 10.6 | 20.1 | 39.3 | 70.0 | 100.7 |
| ν = 0 | 0.8213(2) / 0.7209(271) | 0.9426(4) / 0.9101(718) | 0.9842(5) / 1.0140(934) | 0.9952(11) / 1.1095(2007) | 0.9982(20) / 1.1900(3081) |
| ν = 1 | 0.88221(4) / 0.7135(21) | 0.9649(2) / 0.8960(103) | 0.9904(3) / 0.9812(139) | 0.9965(7) / 1.0357(509) | 0.9977(13) / 1.0879(1022) |
| ν = 3 | 0.5108(9) / 0.8241(51) | 0.8737(26) / 0.9588(130) | 0.9658(45) / 0.9894(156) | 0.9936(104) / 1.0284(347) | 1.0091(180) / 1.0519(661) |
| m²(R−1), ν = 0 | -19.91 / -31.09 | -23.32 / -36.51 | -24.50 / 21.69 | -23.44 / 537.05 | -18.53 / 1928.68 |
| m²(R−1), ν = 1 | -13.12 / -31.92 | -14.27 / -42.21 | -14.86 / -29.14 | -17.04 / 174.98 | -23.03 / 892.67 |
| m²(R−1), ν = 3 | -54.50 / -19.59 | -51.28 / -16.73 | -52.98 / -16.38 | -31.60 / 139.49 | 92.36 / 526.37 |

| β ≈ 2.0, (Δ,J) → | (16.60, 11) | (32.36, 23) | (63.87, 47) | (109.83, 82) |
|---|---|---|---|---|
| m | 9.8 | 20.0 | 40.4 | 70.2 |
| ν = 0 | 0.8150(2) / 0.7178(65) | 0.9290(3) / 0.8918(72) | 0.9796(8) / 0.9764(206) | 0.9934(17) / 1.0130(393) |
| ν = 1 | 1.0991(1) / 0.7551(11) | 1.0395(2) / 0.9061(13) | 1.0120(7) / 0.9723(21) | 1.0030(22) / 0.9955(44) |
| ν = 3 | 0.2232(28) / 2.1945(111) | 0.7790(33) / 1.3278(143) | 0.9548(109) / 1.0407(396) | 1.0106(274) / 0.9688(865) |
| m²(R−1), ν = 0 | -17.71 / -27.02 | -28.37 / -43.25 | -33.30 / -38.49 | -32.58 / 63.85 |
| m²(R−1), ν = 1 | 9.49 / -23.45 | 15.79 / -37.54 | 19.54 / -45.19 | 14.71 / -22.42 |
| m²(R−1), ν = 3 | -74.38 / 114.39 | -88.37 / 131.09 | -73.82 / 66.47 | 52.01 / -153.91 |

#### l = 2, k = 3: R_pres / R_flip (measured/predicted; error in units of the last digit)
| β ≈ 0.4, (Δ,J) → | (14.17, 2) | (24.30, 4) | (44.57, 8) | (74.96, 14) | (105.36, 20) |
|---|---|---|---|---|---|
| m | 12.4 | 22.4 | 42.2 | 72.0 | 101.8 |
| ν = 0 | 0.53930(1) / 0.2936(245) | 0.84003(1) / 0.7134(563) | 0.95323(4) / 0.9557(1558) | 0.9838(1) / 1.0923(2513) | 0.9919(2) / 1.2689(4901) |
| ν = 1 | 0.560989(2) / 0.3006(54) | 0.85153(1) / 0.7023(115) | 0.95701(2) / 0.9229(442) | 0.98508(4) / 1.0196(875) | 0.9925(1) / 1.1207(2052) |
| ν = 3 | 0.79655(1) / 0.3019(22) | 0.96284(3) / 0.7293(56) | 0.9921(1) / 0.9530(287) | 0.9974(2) / 0.9947(320) | 0.9983(5) / 1.0030(616) |
| m²(R−1), ν = 0 | -71.03 / -108.92 | -79.91 / -143.16 | -83.37 / -79.03 | -84.23 / 478.79 | -84.04 / 2787.29 |
| m²(R−1), ν = 1 | -67.69 / -107.83 | -74.17 / -148.72 | -76.62 / -137.44 | -77.40 / 101.50 | -78.12 / 1251.67 |
| m²(R−1), ν = 3 | -31.37 / -107.63 | -18.56 / -135.23 | -14.03 / -83.86 | -13.53 / -27.74 | -17.43 / 30.68 |

| β ≈ 1.0, (Δ,J) → | (13.40, 5) | (24.22, 10) | (45.86, 20) | (80.48, 36) | (115.11, 52) |
|---|---|---|---|---|---|
| m | 10.6 | 20.1 | 39.3 | 70.0 | 100.7 |
| ν = 0 | 0.742421(3) / 0.5597(5) | 0.91165(1) / 0.8212(18) | 0.97484(2) / 0.9464(25) | 0.9919(1) / 0.9857(76) | 0.9961(2) / 0.9991(148) |
| ν = 1 | 0.810536(1) / 0.57263(3) | 0.93961(1) / 0.8299(3) | 0.98332(1) / 0.9489(3) | 0.99465(3) / 0.9842(15) | 0.9974(1) / 0.9944(39) |
| ν = 3 | 0.34377(1) / 0.7232(1) | 0.8130(1) / 0.9112(3) | 0.9511(2) / 0.9755(4) | 0.9845(5) / 0.9928(10) | 0.9927(8) / 0.9984(23) |
| m²(R−1), ν = 0 | -28.69 / -49.05 | -35.87 / -72.60 | -38.93 / -83.03 | -39.72 / -70.20 | -39.55 / -9.52 |
| m²(R−1), ν = 1 | -21.11 / -47.61 | -24.52 / -69.06 | -25.82 / -79.12 | -26.25 / -77.60 | -26.74 / -56.80 |
| m²(R−1), ν = 3 | -73.11 / -30.84 | -75.93 / -36.04 | -75.72 / -37.98 | -76.26 / -35.49 | -73.83 / -16.25 |

| β ≈ 2.0, (Δ,J) → | (16.60, 11) | (32.36, 23) | (63.87, 47) | (109.83, 82) |
|---|---|---|---|---|
| m | 9.8 | 20.0 | 40.4 | 70.2 |
| ν = 0 | 0.66568(1) / 0.5121(2) | 0.86678(1) / 0.7988(2) | 0.96011(4) / 0.9359(8) | 0.9861(2) / 0.9780(20) |
| ν = 1 | 1.017502(2) / 0.56470(3) | 1.024649(4) / 0.82519(4) | 1.01057(3) / 0.9455(2) | 1.0040(2) / 0.9808(3) |
| ν = 3 | -0.3638(1) / 2.2924(2) | 0.5636(1) / 1.4507(4) | 0.8868(4) / 1.1298(11) | 0.9634(13) / 1.0420(32) |
| m²(R−1), ν = 0 | -32.01 / -46.72 | -53.27 / -80.46 | -65.17 / -104.76 | -68.65 / -108.55 |
| m²(R−1), ν = 1 | 1.68 / -41.68 | 9.86 / -69.90 | 17.27 / -89.09 | 19.58 / -94.81 |
| m²(R−1), ν = 3 | -130.59 / 123.76 | -174.48 / 180.24 | -184.94 / 212.09 | -180.33 / 207.19 |

## 6. Non-analytic terms and Assumption T

**The question.** E2E §6.1 (Assumption T) integrates the pointwise heavy expansion term by term, and notes the non-uniform region `u → 0` of the
`(r, η)` integrand: if the expansion breaks down for `u < u_0`, that region contributes `O(u_0^{k+1/2})`; with `u_0 ~ 1/m` this is `m^{−k−1/2}`
(k = 2: `m^{−5/2}`, between the `1/m²` and `1/m³` orders; k = 3: `m^{−7/2}`), with `u_0 ~ 1/m²` it is `m^{−2k−1}`.

**Argument that there are none [derived (sketch)].** By (★) (E2E §3.1; exact, kinematic) the moment equals
`2∫d³x d³y (|x||y|)^{k−4} K(x̂,ŷ) 𝒯(x,y)` with `𝒯` the elementary shadow integrand in the gauge `x_5 = e`, `n = (1,1,0)`. Put `x = X/m`, `y = Y/m`,
`P̂ = p/m`, `ζ_x := −P̂·X ≥ |X|`. (a) The leading integrand is `e^{−2ζ_x−2ζ_y}|X|^{a_x−3}|Y|^{a_y−3}` times polarisation polynomials; the
`(x̂, ŷ)`-integral converges for `a_x, a_y > ½` (here `a_x = 2Δ_φ+k−1 = 4.5` at k = 2). (b) Every Cartesian component of `X` (in the frame of `e`) is
bounded by `cosh d(X̂,e)|X| ≤ e^{β/2}cosh d(X̂,P̂)|X| = e^{β/2}ζ_x`, and the same for `Y`; the polarisation vectors attached to `ŷ` enter exact and
leading tensors homogeneously. (c) `log(b𝒯/𝒯_0)` is a convergent series `Σ_{j≥1} m^{1−j}P_j(X,Y)` with `P_j` polynomials of degree `≤ j` in the components,
for `|x|, |y|` small (the singular loci `N(±x) = 0`, `𝖵 = 0` are at components `~1`). Hence the `n`-th correction is bounded by `C_n(1+ζ)^{2n}`
times the leading integrand, **uniformly in the directions**, and integrable wherever the leading term is; the remainder after `n` terms is
`O(m^{−n−1})` in the region `ζ ≤ m^{1/3}`. The complement has weight `e^{−2m^{1/3}}`, provided the exact integrand is bounded there by
`C e^{−cζ}` times a polynomial — **this outer bound is not written out** [inferred: the exact integrand is supported in the causal diamond and
carries the same exponentials]. So the moment's expansion is in integer powers of `1/m` to all orders; the `u → 0` non-uniformity is an artefact of
the `(r, η)` representation. With E2E's reflection argument (Lemma S plus the antisymmetry of the `1/m` term) the `O(1/m)` term vanishes at
every such k, now without Assumption T.

**Numerical evidence (`u1_uniformity.py`, `out_u1_uniformity.txt`).** `m² · max_pol|𝒯_sym − 𝒯_lead| / max_pol|𝒯_lead|` at fixed `(X, Y)` with typical
moduli `|X| = 1/cosh d(X̂,P̂)`, 5×5 polarisation pairs on the circle attached to `ŷ`, `m = 4.6·10², 4.6·10³, 4.6·10⁴` (converged columns shown):

| `l` | `d(x̂,p̂)` | `d(ŷ,p̄̂)` | `d(x̂,ŷ)` | `m² × rel. deviation` |
|---|---|---|---|---|
| 1 | 0, 2, 4, 6, 8 | 0 | 0.85 … 8.64 | 7.59, 8.50, 8.45, 8.45, 8.45 |
| 1 | 4, 8 (other azimuth) | 0 | 3.52, 7.52 | 6.21, 6.21 |
| 1 | 0 | 4, 8 | 4.58, 8.58 | 3.73, 3.76 |
| 1 | 6 | 6 | 12.6 | 11.88 |
| 2 | 0, 2, 4, 6, 8 | 0 | 0.85 … 8.64 | 19.27, 18.90, 18.99, 18.99, 18.99 |
| 2 | 4, 8 (other azimuth) | 0 | 3.52, 7.52 | 37.73, 37.75 |
| 2 | 0 | 4, 8 | 4.58, 8.58 | 21.07, 21.35 |
| 2 | 6 | 6 | 12.6 | 76.7 |

The coefficient saturates as the distances grow (no exponential growth), as (b)–(c) predict. **Fits of the completed ratios** (§7): an added
`m^{−k−1/2}` term has a fitted coefficient consistent with zero at k = 2 (`−1.3 ± 4.8`, `l = 1`, β ≈ 0.4, ν = 0, pres), i.e. limited to
`|c| ≲ 10 ≈ 0.6|c_2|` at 2σ (the referee found `≲ 6` with the first-version budget) — **an `m^{−5/2}` term with an O(1) coefficient is not excluded by the data**. At k = 3 the fit finds `c_{7/2} = 39 ± 1.2`
(statistical), but the three-parameter fits with and without it are equally good (χ² 0.09 vs 0.3, one degree of freedom), so it cannot be told
apart from the truncation of the analytic series. `m^{−2k−1} = m^{−5}` at k = 2 is invisible behind `c_4/m⁴ ≈ 300/m⁴` and is excluded only by the
sketch. The even-power ansätze (B2, B3) rest on E2E's parity argument: the pointwise `1/m^n` terms of `b𝒯` are odd under `1 ↔ 3` for odd `n`
(reflection `m → −m`), so they cancel in the moments by Lemma S; with the uniformity above this holds order by order.
**Answer to the task question:** the term-by-term expansion does not fail at k = 2 at order `1/m²` [numerical evidence: `m²(R−1)` converges];
at `1/m³`, `m^{−5/2}` and `m^{−5}` the absence of non-analytic terms rests on the derivation sketch [derived (sketch)], with the data limiting
an `m^{−5/2}` coefficient only to `|c| ≲ 10`.

## 7. Corrections at k = 2, 3 (`t2_corrections.py`, `out_t2_corrections.txt`) [fitted]

Weighted fits on `m ≥ 20` (four masses at β ≈ 0.4, 1; three at β ≈ 2, then exact). `c_2` from `c_2/m² + c_4/m⁴` (B2) and with a `c_6/m⁶` term (B3);
`c_1` from `c_1/m + c_2/m² + c_4/m⁴` (A), and `c_2` with an added `m^{−k−½}` (N). Adopted `c_2`: mid-range of the fits with χ²/dof < 4, ± half the
range (the first version quoted the B2/B3 spread although B2 is rejected at k = 3, χ² ≈ 3·10³ — referee M2). Statistical errors are given only
for `c_1` and `c_N`; they ignore the truncation of the ansatz.

| l | β | k | ν | comp | c_2: B2 [χ²] | B3 [χ²] | A [χ²] | N [χ²] | **adopted c_2** | c_1 (A) | c_{k+½} (N) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 0.4 | 2 | 0.0 | pres | -17.084 [2e-01] | -17.041 [1e-01] | -17.142 [9e-02] | -16.919 [9e-02] | **-17.030±0.11** | +0.0011±4.1e-03 | -1.266±4.8e+00 |
| 1 | 0.4 | 2 | 0.0 | flip | -34.767 [1e-01] | -31.787 [4e-02] | -38.732 [2e-02] | -23.342 [3e-02] | **-31.037±7.7** | +0.0770±2.5e-01 | -87.780±2.9e+02 |
| 1 | 0.4 | 2 | 1.0 | pres | -13.708 [9e-02] | -13.691 [6e-02] | -13.731 [5e-02] | -13.643 [5e-02] | **-13.687±0.044** | +0.0004±2.2e-03 | -0.500±2.5e+00 |
| 1 | 0.4 | 2 | 1.0 | flip | -33.181 [6e-02] | -33.723 [3e-02] | -32.562 [3e-02] | -35.119 [3e-02] | **-33.841±1.3** | -0.0123±7.7e-02 | +14.797±9.0e+01 |
| 1 | 0.4 | 2 | 3.0 | pres | 14.835 [1e-01] | 14.461 [3e-02] | 15.302 [3e-02] | 13.442 [3e-02] | **14.372±0.93** | -0.0093±3.2e-02 | +10.623±3.7e+01 |
| 1 | 0.4 | 2 | 3.0 | flip | -12.948 [5e-02] | -11.852 [2e-02] | -14.339 [2e-02] | -8.840 [2e-02] | **-11.590±2.7** | +0.0275±1.4e-01 | -31.364±1.7e+02 |
| 1 | 0.4 | 3 | 0.0 | pres | -27.865 [1e+03] | -27.899 [8e-02] | -27.823 [1e+00] | -27.918 [2e-01] | **-27.871±0.047** | -0.0008±2.5e-05 | +38.978±1.2e+00 |
| 1 | 0.4 | 3 | 0.0 | flip | -54.256 [2e+01] | -54.440 [3e-02] | -54.031 [2e-01] | -54.545 [6e-02] | **-54.288±0.26** | -0.0044±1.0e-03 | +212.386±4.8e+01 |
| 1 | 0.4 | 3 | 1.0 | pres | -23.619 [1e+03] | -23.638 [7e-02] | -23.596 [9e-01] | -23.649 [2e-01] | **-23.623±0.026** | -0.0005±1.4e-05 | +21.583±6.7e-01 |
| 1 | 0.4 | 3 | 1.0 | flip | -51.075 [2e+02] | -51.241 [4e-02] | -50.874 [3e-01] | -51.336 [9e-02] | **-51.105±0.23** | -0.0040±2.9e-04 | +189.823±1.4e+01 |
| 1 | 0.4 | 3 | 3.0 | pres | 12.531 [3e+01] | 12.551 [4e-02] | 12.506 [1e-01] | 12.563 [5e-02] | **12.534±0.028** | +0.0005±8.9e-05 | -23.202±4.2e+00 |
| 1 | 0.4 | 3 | 3.0 | flip | -25.799 [3e-01] | -25.790 [2e-02] | -25.810 [1e-02] | -25.785 [2e-02] | **-25.798±0.013** | +0.0002±3.8e-04 | -10.266±1.8e+01 |
| 1 | 1.0 | 2 | 0.0 | pres | -15.786 [2e-01] | -15.722 [2e-01] | -15.863 [1e-01] | -15.553 [2e-01] | **-15.708±0.15** | +0.0017±6.1e-03 | -1.704±6.5e+00 |
| 1 | 1.0 | 2 | 0.0 | flip | -29.330 [1e-01] | -28.909 [8e-02] | -29.833 [7e-02] | -27.807 [8e-02] | **-28.820±1** | +0.0108±4.3e-02 | -11.141±4.6e+01 |
| 1 | 1.0 | 2 | 1.0 | pres | -8.195 [1e-01] | -8.144 [5e-02] | -8.254 [5e-02] | -8.014 [5e-02] | **-8.134±0.12** | +0.0013±4.7e-03 | -1.322±5.0e+00 |
| 1 | 1.0 | 2 | 1.0 | flip | -24.956 [2e-01] | -24.726 [1e-01] | -25.229 [1e-01] | -24.126 [1e-01] | **-24.678±0.55** | +0.0059±2.3e-02 | -6.065±2.5e+01 |
| 1 | 1.0 | 2 | 3.0 | pres | -18.794 [6e-03] | -18.602 [5e-03] | -19.029 [5e-03] | -18.090 [5e-03] | **-18.560±0.47** | +0.0051±1.6e-01 | -5.139±1.7e+02 |
| 1 | 1.0 | 2 | 3.0 | flip | 10.710 [3e-03] | 10.509 [5e-05] | 10.932 [9e-05] | 10.007 [7e-05] | **10.470±0.46** | -0.0048±9.3e-02 | +5.136±9.9e+01 |
| 1 | 1.0 | 3 | 0.0 | pres | -30.663 [5e+03] | -30.797 [1e-01] | -30.515 [3e+00] | -30.867 [4e-01] | **-30.691±0.18** | -0.0032±4.5e-05 | +130.530±1.8e+00 |
| 1 | 1.0 | 3 | 0.0 | flip | -52.441 [3e+03] | -52.869 [7e-02] | -51.963 [3e+00] | -53.096 [3e-01] | **-52.529±0.57** | -0.0103±1.8e-04 | +419.861±7.3e+00 |
| 1 | 1.0 | 3 | 1.0 | pres | -19.724 [4e+02] | -19.750 [9e-02] | -19.695 [4e-01] | -19.763 [1e-01] | **-19.729±0.034** | -0.0006±3.3e-05 | +25.385±1.3e+00 |
| 1 | 1.0 | 3 | 1.0 | flip | -46.462 [4e+03] | -46.774 [1e-01] | -46.114 [2e+00] | -46.940 [4e-01] | **-46.527±0.41** | -0.0075±1.1e-04 | +305.958±4.6e+00 |
| 1 | 1.0 | 3 | 3.0 | pres | -45.579 [2e+02] | -45.316 [2e-02] | -45.872 [8e-03] | -45.177 [8e-03] | **-45.525±0.35** | +0.0063±5.0e-04 | -257.384±2.0e+01 |
| 1 | 1.0 | 3 | 3.0 | flip | 3.903 [2e+02] | 4.067 [1e-03] | 3.721 [7e-02] | 4.153 [1e-03] | **3.937±0.22** | +0.0039±2.5e-04 | -160.181±1.0e+01 |
| 1 | 2.0 | 2 | 0.0 | pres | -33.443 [1e-01] | -33.653 [5e-23] | -33.204 [5e-23] | -34.184 [5e-22] | **-33.694±0.49** | -0.0051±1.6e-02 | +5.460±1.7e+01 |
| 1 | 2.0 | 2 | 0.0 | flip | -47.834 [5e-03] | -47.965 [5e-23] | -47.685 [3e-24] | -48.297 [5e-23] | **-47.991±0.31** | -0.0032±4.4e-02 | +3.408±4.7e+01 |
| 1 | 2.0 | 2 | 1.0 | pres | 20.651 [3e-01] | 21.376 [2e-22] | 19.824 [5e-23] | 23.217 [1e-22] | **21.521±1.7** | +0.0176±3.2e-02 | -18.905±3.4e+01 |
| 1 | 2.0 | 2 | 1.0 | flip | -31.793 [2e-01] | -31.141 [5e-23] | -32.536 [3e-22] | -29.486 [8e-25] | **-31.011±1.5** | +0.0158±3.5e-02 | -16.989±3.8e+01 |
| 1 | 2.0 | 2 | 3.0 | pres | -75.865 [3e-02] | -70.892 [0e+00] | -81.537 [3e-24] | -58.259 [3e-24] | **-69.898±12** | +0.1205±7.3e-01 | -129.671±7.9e+02 |
| 1 | 2.0 | 2 | 3.0 | flip | 5.917 [7e-02] | 10.648 [3e-26] | 0.521 [3e-26] | 22.666 [2e-24] | **11.593±11** | +0.1147±4.5e-01 | -123.360±4.8e+02 |
| 1 | 2.0 | 3 | 0.0 | pres | -67.285 [6e+04] | -68.630 [6e-19] | -65.753 [2e-18] | -69.350 [6e-14] | **-67.551±1.8** | -0.0326±1.3e-04 | +1350.683±5.4e+00 |
| 1 | 2.0 | 3 | 0.0 | flip | -90.767 [6e+04] | -93.273 [2e-20] | -87.907 [5e-18] | -94.614 [1e-14] | **-91.260±3.4** | -0.0608±2.4e-04 | +2518.705±1.0e+01 |
| 1 | 2.0 | 3 | 1.0 | pres | 18.475 [3e+04] | 20.012 [1e-20] | 16.723 [9e-20] | 20.834 [1e-15] | **18.779±2.1** | +0.0372±2.1e-04 | -1543.422±8.6e+00 |
| 1 | 2.0 | 3 | 1.0 | flip | -67.194 [2e+04] | -68.398 [1e-18] | -65.821 [5e-19] | -69.042 [1e-18] | **-67.431±1.6** | -0.0292±2.0e-04 | +1209.382±8.5e+00 |
| 1 | 2.0 | 3 | 3.0 | pres | -178.365 [1e+03] | -174.957 [5e-22] | -182.249 [4e-21] | -173.134 [3e-17] | **-177.691±4.6** | +0.0826±2.4e-03 | -3422.548±9.9e+01 |
| 1 | 2.0 | 3 | 3.0 | flip | -29.790 [2e+03] | -27.138 [2e-21] | -32.813 [2e-21] | -25.720 [2e-21] | **-29.266±3.5** | +0.0643±1.4e-03 | -2663.602±5.6e+01 |
| 2 | 0.4 | 2 | 0.0 | pres | -73.584 [2e-01] | -69.441 [5e-02] | -78.807 [2e-02] | -57.990 [3e-02] | **-68.399±10** | +0.1016±2.4e-01 | -119.708±2.9e+02 |
| 2 | 0.4 | 2 | 0.0 | flip | 6019.412 [6e-01] | 45525.841 [1e-01] | -43327.069 [4e-02] | 154141.629 [8e-02] | **55407.280±9.9e+04** | +949.6238±1.2e+03 | -1141098.245±1.5e+06 |
| 2 | 0.4 | 2 | 1.0 | pres | -68.346 [9e-01] | -73.099 [2e-01] | -62.402 [6e-02] | -86.173 [1e-01] | **-74.287±12** | -0.1155±1.3e-01 | +136.916±1.6e+02 |
| 2 | 0.4 | 2 | 1.0 | flip | 2351.146 [1e+00] | 24591.902 [3e-01] | -25278.526 [1e-01] | 85579.469 [2e-01] | **30150.471±5.5e+04** | +541.3953±5.2e+02 | -637446.421±6.3e+05 |
| 2 | 0.4 | 2 | 3.0 | pres | -16.910 [5e-01] | -47.843 [5e-02] | 20.140 [8e-03] | -130.565 [3e-02] | **-55.212±75** | -0.7271±1.0e+00 | +870.075±1.2e+03 |
| 2 | 0.4 | 2 | 3.0 | flip | 1794.843 [2e-01] | -294.448 [2e-01] | 5137.973 [1e-01] | -7093.816 [2e-01] | **-977.921±6.1e+03** | -63.4323±1.9e+02 | +68800.062±2.4e+05 |
| 2 | 0.4 | 3 | 0.0 | pres | -84.709 [7e-02] | -84.630 [5e-02] | -84.830 [4e-02] | -84.576 [5e-02] | **-84.703±0.13** | +0.0024±1.4e-02 | -98.336±7.1e+02 |
| 2 | 0.4 | 3 | 0.0 | flip | 0.296 [5e-01] | 1140.809 [2e-01] | -1535.641 [8e-02] | 1846.507 [1e-01] | **155.433±1.7e+03** | +29.5116±4.8e+01 | -1380631.835±2.4e+06 |
| 2 | 0.4 | 3 | 1.0 | pres | -77.581 [1e+00] | -77.882 [3e-01] | -77.191 [2e-01] | -78.063 [3e-01] | **-77.627±0.44** | -0.0076±7.3e-03 | +357.070±3.6e+02 |
| 2 | 0.4 | 3 | 1.0 | flip | -117.903 [7e-01] | 366.361 [2e-01] | -758.451 [1e-01] | 662.065 [2e-01] | **-48.193±7.1e+02** | +12.4813±1.7e+01 | -577051.616±8.5e+05 |
| 2 | 0.4 | 3 | 3.0 | pres | -12.291 [1e+00] | -13.640 [6e-01] | -10.500 [4e-01] | -14.465 [5e-01] | **-12.482±2** | -0.0351±3.5e-02 | +1600.334±1.7e+03 |
| 2 | 0.4 | 3 | 3.0 | flip | -56.757 [8e-02] | 14.216 [4e-03] | -152.081 [6e-05] | 57.583 [2e-03] | **-47.249±1e+02** | +1.7651±6.3e+00 | -88123.238±3.2e+05 |
| 2 | 1.0 | 2 | 0.0 | pres | -24.877 [2e-01] | -22.196 [5e-02] | -28.004 [2e-02] | -15.210 [4e-02] | **-21.607±6.4** | +0.0666±1.8e-01 | -70.962±2.0e+02 |
| 2 | 1.0 | 2 | 0.0 | flip | 66.861 [6e-01] | 1058.402 [1e-01] | -1063.264 [5e-02] | 3609.282 [9e-02] | **1273.009±2.3e+03** | +23.9020±3.2e+01 | -26066.649±3.6e+04 |
| 2 | 1.0 | 2 | 1.0 | pres | -15.123 [7e-01] | -18.958 [1e-01] | -10.776 [6e-02] | -28.773 [1e-01] | **-19.774±9** | -0.0931±1.1e-01 | +100.056±1.3e+02 |
| 2 | 1.0 | 2 | 1.0 | flip | -21.316 [1e+00] | 367.188 [4e-01] | -466.372 [2e-01] | 1368.299 [3e-01] | **450.964±9.2e+02** | +9.5651±8.6e+00 | -10174.849±9.5e+03 |
| 2 | 1.0 | 2 | 3.0 | pres | -52.580 [8e-01] | -2.237 [4e-01] | -114.195 [2e-01] | 133.335 [3e-01] | **9.570±1.2e+02** | +1.3088±1.7e+00 | -1366.269±1.9e+03 |
| 2 | 1.0 | 2 | 3.0 | flip | -9.966 [1e+00] | 272.775 [2e-01] | -332.994 [8e-02] | 998.368 [2e-01] | **332.687±6.7e+02** | +6.8528±5.8e+00 | -7413.318±6.5e+03 |
| 2 | 1.0 | 3 | 0.0 | pres | -40.025 [6e-02] | -40.051 [5e-02] | -40.016 [6e-02] | -40.058 [6e-02] | **-40.037±0.021** | -0.0002±1.1e-02 | +21.682±4.7e+02 |
| 2 | 1.0 | 3 | 0.0 | flip | -86.395 [4e-01] | -55.702 [1e-01] | -122.088 [7e-02] | -38.736 [1e-01] | **-80.412±42** | +0.7649±1.3e+00 | -30743.130±5.5e+04 |
| 2 | 1.0 | 3 | 1.0 | pres | -26.278 [9e-01] | -26.496 [2e-01] | -26.030 [9e-02] | -26.614 [2e-01] | **-26.322±0.29** | -0.0054±6.0e-03 | +215.843±2.5e+02 |
| 2 | 1.0 | 3 | 1.0 | flip | -82.669 [7e-01] | -75.037 [3e-01] | -91.615 [2e-01] | -70.819 [2e-01] | **-81.217±10** | +0.1932±2.6e-01 | -7603.185±1.1e+04 |
| 2 | 1.0 | 3 | 3.0 | pres | -75.648 [1e-01] | -76.224 [9e-02] | -75.260 [1e-01] | -76.445 [1e-01] | **-75.852±0.59** | -0.0083±7.3e-02 | +512.597±3.2e+03 |
| 2 | 1.0 | 3 | 3.0 | flip | -38.593 [1e+00] | -32.482 [6e-01] | -46.041 [4e-01] | -29.013 [5e-01] | **-37.527±8.5** | +0.1595±1.7e-01 | -6182.199±7.1e+03 |
| 2 | 2.0 | 2 | 0.0 | pres | -34.840 [5e-02] | -31.979 [1e-26] | -38.227 [1e-26] | -24.564 [1e-26] | **-31.396±6.8** | +0.0707±3.3e-01 | -76.111±3.6e+02 |
| 2 | 2.0 | 2 | 0.0 | flip | -32.223 [3e-01] | 127.981 [1e-27] | -224.872 [1e-27] | 546.738 [1e-29] | **160.933±3.9e+02** | +3.9952±7.8e+00 | -4298.285±8.4e+03 |
| 2 | 2.0 | 2 | 1.0 | pres | 20.659 [3e-01] | 11.572 [1e-26] | 31.221 [3e-26] | -11.747 [5e-26] | **9.737±21** | -0.2225±4.3e-01 | +239.355±4.6e+02 |
| 2 | 2.0 | 2 | 1.0 | flip | -46.729 [1e+00] | -7.894 [2e-27] | -92.999 [9e-27] | 93.106 [1e-26] | **0.054±93** | +0.9636±8.5e-01 | -1036.704±9.2e+02 |
| 2 | 2.0 | 2 | 3.0 | pres | -65.780 [8e-01] | 130.605 [6e-26] | -296.651 [6e-27] | 637.662 [3e-26] | **170.506±4.7e+02** | +4.8377±5.4e+00 | -5204.627±5.8e+03 |
| 2 | 2.0 | 2 | 3.0 | flip | 38.113 [2e-01] | -290.458 [6e-28] | 428.356 [3e-28] | -1143.529 [2e-26] | **-357.586±7.9e+02** | -8.1389±1.7e+01 | +8756.248±1.8e+04 |
| 2 | 2.0 | 3 | 0.0 | pres | -69.040 [2e+00] | -70.492 [1e-26] | -67.352 [2e-22] | -71.277 [9e-21] | **-69.314±2** | -0.0356±2.6e-02 | +1473.918±1.1e+03 |
| 2 | 2.0 | 3 | 0.0 | flip | -112.593 [2e-02] | -110.233 [5e-26] | -115.361 [6e-26] | -108.951 [4e-22] | **-112.156±3.2** | +0.0581±3.8e-01 | -2406.909±1.6e+04 |
| 2 | 2.0 | 3 | 1.0 | pres | 19.675 [1e+00] | 20.823 [3e-27] | 18.353 [1e-23] | 21.441 [1e-21] | **19.897±1.5** | +0.0280±2.6e-02 | -1159.349±1.1e+03 |
| 2 | 2.0 | 3 | 1.0 | flip | -95.380 [2e+00] | -97.852 [7e-24] | -92.411 [1e-23] | -99.213 [2e-19] | **-95.812±3.4** | -0.0616±4.2e-02 | +2554.347±1.7e+03 |
| 2 | 2.0 | 3 | 3.0 | pres | -188.229 [1e+00] | -177.142 [3e-24] | -201.082 [8e-24] | -171.156 [3e-23] | **-186.119±15** | +0.2711±2.5e-01 | -11237.058±1.0e+04 |
| 2 | 2.0 | 3 | 3.0 | flip | 222.183 [6e-01] | 203.203 [5e-26] | 244.305 [2e-25] | 192.926 [1e-21] | **218.616±26** | -0.4654±6.2e-01 | +19293.043±2.6e+04 |

## 8. What remains

1. **`l = 2` deepest corner.** The spinning partial wave is built from up to fourth derivatives of seed blocks whose crossed-channel singularity
   is `(1−ρ_b)^{−2σ}` with `σ` up to 3.25; the conserved combination cancels between seeds by up to ~`(1−ρ_b)^{−(2σ+4)}` relative to its value,
   so even jets accurate to `10⁻¹²–10⁻⁸` (§2.3) do not give the `l = 2` tensor where `1−ρ_b ≲ 10⁻³`. Fix: apply the differential operators
   to the lightcone *coefficient arrays* (exact termwise algebra on series, cancellation at coefficient level) and fit the tail of the
   resulting tensor-structure coefficients, whose singularity is the mild physical one; or quad-precision fits. Until then the `l = 2`
   helicity-flip component at k = 2 (and partly k = 3) is not determined (§5.2).
2. **Tiny `t`** (both radial coordinates near 1; the u-channel OPE corner of the frame): neither the radial series (`N = 10⁴`) nor the lightcone
   sum (`I = 1500`) converges; these points are power-law estimates with 50 % error, which dominates the k = 2 error at `m ≈ 100`. Fix: a larger
   lightcone depth (`I ~ 60/(4s)` ≈ 3000–6000; memory ~5×(I×Q)×8 B per seed) or a representation centred at that corner.
3. **Rigour of §6:** the outer-region bound on the exact shadow integrand (`ζ ≳ m^{1/3}`) is not written; the conclusion "no non-analytic terms"
   is a derivation sketch plus numerical evidence.
4. **Analytic `c_2`** at k = 2, 3 (E2E §11 item 2): not attempted; the fitted values of §7 are the input for a check.
5. **Scalar stress test** (§4.1): the corner-heavy scalar `(3.75, 0.5)` at k = 2 is limited at `m ≥ 70` by the tiny-`t` estimate (item 2), with
   errors that cover the deviations; k = 3 is clean.
6. **e2e code issue found:** `ads4/code/e2e/cblock.c` corrupts memory for `N < J` (`CBlock(105.36, 20, ·, ·, 4)` segfaults; `J = 82`, `N = 60`
   aborts later). Worked around here (`N = J + 60` for placeholder seeds); not fixed (e2e code is not to be edited by this task).
7. **Repository hygiene:** an early work-in-progress state of `ads4/code/lowk/` was committed by someone else's commit (including two
   `cache/*.npz` files) before `ads4/code/lowk/cache/` was added to `.gitignore`; those cached files remain tracked until removed with
   `git rm --cached`.
8. **Consequence for E2E:** its 72 × 62 grid has a k = 2 quadrature error of up to `~10⁻⁴` at ν = 3 near `u = 0` (§4.3; scalar `(3, 1.75)`: up to
   `5·10⁻⁵`), which matters for E2E's k = 2 scalar-skeleton numbers (§4.1 there: `m²(R−1)` at ν = 3, m ≈ 120, could shift by `O(1)`); its k = 4, 6
   numbers are not affected at their quoted accuracy (§4.2).
9. Unchanged from E2E §11: `⟨TTTT⟩`, MFT normalisation (see `SPINNING_MFT.md`), `TT → φφ`, other `Δ_φ`, `ν ~ m`.

## 9. Code index (`ads4/code/lowk/`; `sh run_all.sh`)

| file | content | output |
|---|---|---|
| `lcconv.c` | `B_{n,j}` → lightcone coefficients `C^{(l)}_{i,q}` (Gegenbauer, all η-derivatives), C/OpenMP | — |
| `lcblock.py` | `LCBlock`: Darboux tail fits, closed-form tails, bucketed evaluation of `S[k,l]` | — |
| `casimir_check.py` | Casimir residual from the jets | — |
| `c1_block_checks.py` | block-level checks (§2.3) | `out_c1_block_checks.txt` |
| `corner.py` | corner completion of the e2e mode grids, corner-refined grid, two settings, acceptance, estimates, error budget | — |
| `v1_scalar.py` | scalar validation (§4.1) | `out_v1_scalar.txt` |
| `t1_lowk_e2e.py` | completed moments vs the leading formula, `l = 1, 2`, k = 2, 3, 4, 6 | `out_t1_lowk_e2e_l{1,2}.txt` |
| `q1_grid.py` | quadrature-grid check (§4.3) | `out_q1_grid.txt` |
| `u1_uniformity.py` | uniformity of the heavy expansion of the shadow integrand (§6) | `out_u1_uniformity.txt` |
| `t2_corrections.py` | `m²(R−1)`, `c_1`, `c_2`, non-analytic fits (§7) | `out_t2_corrections.txt` |
| `t3_tables.py` | tables of this note | `out_t3_tables.txt` |

Caches (`cache/`, git-ignored): `lcmodes_*.npz` (completed mode grids), `data_t1_l*.pkl`, run logs. Run time on 4 shared cores: the corner
completion is ~10 min per `(l = 1)` case and ~15–20 min per `(l = 2)` case; `v1`, `q1` ~1–2 h each.
