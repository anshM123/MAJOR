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
2. **Currents (`l = 1`), k = 2 and 3: converged and error-controlled [numerical evidence for the ratios; the formula itself is the derived
   one of E2E §4.5, nothing fitted].** At the largest masses (`m ≈ 100`; β ≈ 0.4 / 1.0, ν = 0):
   * k = 3: `R_pres = 0.997317(2) / 0.99698(1)`, `R_flip = 0.99478(2) / 0.99483(2)`; `m²(R−1)` converges like `c_2 + c_4/m²`
     (e.g. β ≈ 0.4, ν = 0: `−23.68, −26.39, −27.45, −27.74, −27.81` (pres), `−40.5, −49.1, −52.8, −53.9, −54.1` (flip), `m = 12 … 102`).
   * k = 2: `R_pres = 0.9984(2) / 0.9986(3)`, `R_flip = 0.9972(30) / 0.9976(16)`; `m²(R−1)` is converged within errors from `m ≈ 40`
     (β ≈ 0.4, ν = 0: `−15.36, −16.49, −16.92, −17.02(6), −16.65(120)` (pres), `−27.7, −32.6, −34.7, −34.8(18), −29(30)` (flip)).
   * all three β (0.4, 1, 2) and ν = 0, 1, 3: the full tables are §5; at β ≈ 2 (only `m ≤ 70`) `m²(R−1)` is still drifting for ν = 3.
3. **Stress tensors (`l = 2`) [partial].** Helicity-**preserved** component converged at k = 3 and consistent at k = 2: k = 3, β ≈ 0.4: `R_pres = 0.53930(1), 0.84002(1),
   0.953220(5), 0.98375(1), 0.99185(2)` (`m = 12 … 102`; `m²(R−1) = −71.0, −79.9, −83.4, −84.3, −84.5`); k = 2: `0.5725(2) … 0.9931(7)`,
   `m²(R−1) = −65.9, −71.3, −73.1, −73.1(21), −72(7)`; β ≈ 1, m ≈ 101: `R_pres = 0.99606(3)` (k = 3), `0.9977(7)` (k = 2). Helicity-**flip**
   component: k = 3 determined at β ≈ 1, 2 (`R_flip = 0.9925(31)` at m ≈ 101, β ≈ 1; `0.9771(4)` at m ≈ 70, β ≈ 2) but only to `±0.09` at β ≈ 0.4;
   at k = 2 it is **not determined** at β ≈ 0.4 (errors > 1) and only to a few % at β ≈ 1, 2. Reason: for `l = 2` the deepest corner and the
   tiny-t region carry up to 10 % of the `(2,2)` mode at k = 2 and cannot be evaluated in double precision (cancellation between the
   seeds, §8 item 1); they are estimated with 50 % error, which dominates.
4. **No `O(1/m)` term at k = 2, 3 [numerical evidence; the argument of E2E §6.1 now holds without Assumption T, see 5].** Fitted `c_1`
   (ansatz `c_1/m + c_2/m² + c_4/m⁴`, `m ≥ 20`): `|c_1| ≤ 0.017` (k = 2) and `≤ 0.011` (k = 3) for β ≤ 1 (`l = 1`), with statistical errors of
   the same size at k = 2.
5. **Non-analytic terms / Assumption T [derived (sketch) + numerical evidence].** Through the exact identity (★) of E2E §3.1 the moment
   is a 6-dimensional integral of the elementary shadow integrand; in heavy variables its expansion in `1/m` is uniform in all directions
   (every Cartesian component of `X = mx` is bounded by `e^{β/2} ζ`, `ζ = −P̂·X`, which the exponential `e^{−2ζ}` controls), so the moment has
   an asymptotic expansion in powers of `1/m` with no non-analytic terms — neither `m^{−k−1/2}` (breakdown scale `u ~ 1/m`) nor
   `m^{−2k−1}` (`u ~ 1/m²`) — for every k for which (★) converges (for `Δ_φ = 1.75`, `l ≤ 2`: all k > −2). The non-uniformity at `u → 0` that motivated
   Assumption T is a property of the `(r, η)` representation only. Numerically: the pointwise `1/m²` coefficient of the shadow integrand
   stays bounded (it saturates) when the hyperbolic distances of `x̂, ŷ` from `p̂, p̄̂` grow to 8 (`d(x̂, ŷ)` up to 12.6), `l = 1, 2` (§6);
   fits of the completed k = 2 ratios with an added `m^{−5/2}` term give a coefficient consistent with zero (`−0.5 ± 2.6` (first-version budget; to be refreshed) against
   `c_2 ≈ −17.1`, `l = 1`, β ≈ 0.4, ν = 0), but only `|c| ≲ 6` (≈ 0.35|c₂|) is excluded — an O(1) coefficient is allowed by the data;
   `m^{−5}` at k = 2 and `m^{−7/2}` at k = 3 cannot be separated from the analytic truncation at all. So: no failure at `1/m²`
   [numerical evidence]; beyond that the absence of non-analytic terms rests on the sketch; no rigorous remainder estimate is written.
6. **Fitted `O(1/m²)` coefficients** (`l = 1`, k = 2 / 3, β ≈ 0.4, ν = 0; fitted, ansatz spread as error): `c_2^pres = −17.08(2) / −27.88(2)`,
   `c_2^flip = −35.1(4) / −54.35(10)`. Table §7.
7. **What is not done:** see §8 (most importantly: the `l = 2` deepest corner at k = 2, 3 is limited by double precision; scalar-test caveat).

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

**Per ratio.** The relative errors of the modes `(l,−l)`, `(l,l)` are propagated through the 2×2 helicity solve (E2E §5.2) to `R_pres`, `R_flip`.
Validation of the re-evaluated converged points: `max dev ≤ 5·10⁻⁹` for `l = 1` (all cases); for `l = 2` it ranges from `2·10⁻⁶` (β ≈ 1, 2, m ≥ 40) up to `0.2–0.3` (m ≈ 10–20) and `1.5·10²` at the single worst point (m ≈ 12, β ≈ 0.4) — the
`l = 2` tensor at small `s` amplifies `10⁻¹²` errors of the seed jets by up to `10¹⁰`. Its impact on the moments, measured by evaluating **all**
converged points with `u < 0.25` with LC and weighting the deviations (`v2_l2_validation.py`, `out_v2_l2_validation.txt`): `l = 2`, m ≈ 42:
`1.1·10⁻⁵ / 8.3·10⁻⁴` (k = 2, modes `(2,−2)/(2,2)`), `1.1·10⁻⁶ / 2.2·10⁻⁴` (k = 3); m ≈ 102: `2·10⁻⁷ / 1.2·10⁻⁵` (k = 2); `l = 1`: `≤ 3·10⁻¹¹`.
These numbers are not in the `l = 2` budget (the corner points have their own A/B test); they are smaller than the quoted `l = 2` errors
except for the flip component at k = 3, m ≈ 42, β ≈ 0.4, where they are comparable.

## 4. Validation of the completed moments

### 4.1 Scalar partial waves (`v1_scalar.py`, `out_v1_scalar.txt`)

`⟨φ_1φ_2φ_1φ_2⟩`, s-channel, k = 2, 3, compared with the derived scalar formula of E2E §4.1 (nothing fitted); errors as in §3.

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
(β ≤ 1); k = 4 flip differs by up to `5·10⁻³` (m ≈ 42, β ≈ 0.4; E2E `*`).

### 4.3 Quadrature (`q1_grid.py`, `out_q1_grid.txt`)

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
| ν = 0 | 0.900380(1) / 0.82062(3) | 0.966982(1) / 0.93477(5) | 0.990505(1) / 0.9806(1) | 0.99672(1) / 0.9933(4) | 0.9984(2) / 0.9972(30) |
| ν = 1 | 0.917688(1) / 0.830898(2) | 0.973238(1) / 0.939064(3) | 0.992361(1) / 0.98190(1) | 0.99736(1) / 0.9937(2) | 0.9987(1) / 0.9973(22) |
| ν = 3 | 1.084005(2) / 0.91982(1) | 1.028496(3) / 0.97440(2) | 1.00823(1) / 0.99281(3) | 1.0028(2) / 0.9975(5) | 1.0013(13) / 0.9991(49) |
| m²(R−1), ν = 0 | -15.36 / -27.66 | -16.49 / -32.59 | -16.92 / -34.65 | -17.02 / -34.81 | -16.65 / -28.86 |
| m²(R−1), ν = 1 | -12.69 / -26.07 | -13.37 / -30.44 | -13.62 / -32.26 | -13.67 / -32.62 | -13.46 / -27.83 |
| m²(R−1), ν = 3 | 12.95 / -12.36 | 14.24 / -12.79 | 14.67 / -12.82 | 14.75 / -12.83 | 13.12 / -9.25 |

| β ≈ 1.0, (Δ,J) → | (13.40, 5) | (24.22, 10) | (45.86, 20) | (80.48, 36) | (115.11, 52) |
|---|---|---|---|---|---|
| m | 10.6 | 20.1 | 39.3 | 70.0 | 100.7 |
| ν = 0 | 0.888528(1) / 0.823242(2) | 0.963971(1) / 0.936452(2) | 0.989995(1) / 0.98163(1) | 0.99681(3) / 0.9941(2) | 0.9986(3) / 0.9976(16) |
| ν = 1 | 0.934508(1) / 0.844593(1) | 0.980452(1) / 0.945291(1) | 0.994747(1) / 0.984345(5) | 0.99834(2) / 0.9950(2) | 0.9993(3) / 0.9980(12) |
| ν = 3 | 0.774783(2) / 1.056069(1) | 0.947628(2) / 1.022406(2) | 0.98749(3) / 1.00668(2) | 0.9962(8) / 1.0022(5) | 0.9988(80) / 1.0011(44) |
| m²(R−1), ν = 0 | -12.42 / -19.69 | -14.63 / -25.80 | -15.48 / -28.43 | -15.66 / -28.90 | -14.49 / -24.54 |
| m²(R−1), ν = 1 | -7.30 / -17.31 | -7.94 / -22.21 | -8.13 / -24.23 | -8.15 / -24.57 | -7.60 / -20.48 |
| m²(R−1), ν = 3 | -25.09 / 6.25 | -21.26 / 9.10 | -19.36 / 10.34 | -18.54 / 10.74 | -11.96 / 11.45 |

| β ≈ 2.0, (Δ,J) → | (16.60, 11) | (32.36, 23) | (63.87, 47) | (109.83, 82) |
|---|---|---|---|---|
| m | 9.8 | 20.0 | 40.4 | 70.2 |
| ν = 0 | 0.800438(1) / 0.772104(1) | 0.929849(1) / 0.906320(2) | 0.98034(1) / 0.97227(2) | 0.9933(1) / 0.9904(3) |
| ν = 1 | 1.056806(1) / 0.829844(1) | 1.035686(1) / 0.934404(1) | 1.01168(1) / 0.98137(1) | 1.0042(2) / 0.9937(2) |
| ν = 3 | 0.166538(5) / 0.878094(3) | 0.78323(2) / 0.98783(1) | 0.9520(3) / 1.0021(2) | 0.9852(38) / 1.0017(23) |
| m²(R−1), ν = 0 | -19.11 / -21.82 | -28.05 / -37.46 | -32.12 / -45.30 | -33.14 / -47.09 |
| m²(R−1), ν = 1 | 5.44 / -16.29 | 14.27 / -26.23 | 19.09 / -30.43 | 20.57 / -30.94 |
| m²(R−1), ν = 3 | -79.81 / -11.67 | -86.68 / -4.87 | -78.36 / 3.38 | -73.17 / 8.31 |

#### l = 1, k = 3: R_pres / R_flip (measured/predicted; error in units of the last digit)
| β ≈ 0.4, (Δ,J) → | (14.17, 2) | (24.30, 4) | (44.57, 8) | (74.96, 14) | (105.36, 20) |
|---|---|---|---|---|---|
| m | 12.4 | 22.4 | 42.2 | 72.0 | 101.8 |
| ν = 0 | 0.846393(1) / 0.737319(1) | 0.947179(1) / 0.901817(1) | 0.984599(1) / 0.970379(1) | 0.994651(1) / 0.989614(1) | 0.997317(2) / 0.99478(2) |
| ν = 1 | 0.866025(1) / 0.749008(1) | 0.954769(1) / 0.907049(1) | 0.986909(1) / 0.972075(1) | 0.995463(1) / 0.990220(1) | 0.997725(1) / 0.99509(2) |
| ν = 3 | 1.057552(1) / 0.851083(1) | 1.022568(1) / 0.950099(1) | 1.006833(1) / 0.985663(1) | 1.002395(1) / 0.995044(1) | 1.00120(1) / 0.99752(3) |
| m²(R−1), ν = 0 | -23.68 / -40.50 | -26.39 / -49.05 | -27.45 / -52.79 | -27.74 / -53.87 | -27.81 / -54.11 |
| m²(R−1), ν = 1 | -20.66 / -38.70 | -22.60 / -46.43 | -23.33 / -49.77 | -23.53 / -50.72 | -23.58 / -50.93 |
| m²(R−1), ν = 3 | 8.87 / -22.96 | 11.27 / -24.93 | 12.18 / -25.55 | 12.42 / -25.71 | 12.47 / -25.71 |

| β ≈ 1.0, (Δ,J) → | (13.40, 5) | (24.22, 10) | (45.86, 20) | (80.48, 36) | (115.11, 52) |
|---|---|---|---|---|---|
| m | 10.6 | 20.1 | 39.3 | 70.0 | 100.7 |
| ν = 0 | 0.799456(1) / 0.699700(1) | 0.931838(1) / 0.888731(1) | 0.980694(1) / 0.967348(1) | 0.993784(1) / 0.989379(1) | 0.99698(1) / 0.99483(2) |
| ν = 1 | 0.856115(1) / 0.724263(1) | 0.954378(1) / 0.900062(1) | 0.987459(1) / 0.970977(1) | 0.995997(1) / 0.990589(1) | 0.998061(4) / 0.99542(1) |
| ν = 3 | 0.582111(1) / 0.982455(1) | 0.884272(1) / 1.003612(1) | 0.970311(1) / 1.002109(1) | 0.990720(3) / 1.000775(1) | 0.9955(1) / 1.00039(3) |
| m²(R−1), ν = 0 | -22.34 / -33.45 | -27.67 / -45.18 | -29.88 / -50.53 | -30.50 / -52.11 | -30.63 / -52.46 |
| m²(R−1), ν = 1 | -16.03 / -30.72 | -18.52 / -40.58 | -19.41 / -44.92 | -19.64 / -46.17 | -19.69 / -46.44 |
| m²(R−1), ν = 3 | -46.55 / -1.95 | -46.99 / 1.47 | -45.95 / 3.26 | -45.53 / 3.80 | -45.33 / 3.95 |

| β ≈ 2.0, (Δ,J) → | (16.60, 11) | (32.36, 23) | (63.87, 47) | (109.83, 82) |
|---|---|---|---|---|
| m | 9.8 | 20.0 | 40.4 | 70.2 |
| ν = 0 | 0.645885(1) / 0.581859(1) | 0.865664(1) / 0.827279(2) | 0.960846(1) / 0.947692(1) | 0.986403(1) / 0.981627(1) |
| ν = 1 | 0.968816(1) / 0.650313(1) | 1.019637(1) / 0.865391(1) | 1.009717(1) / 0.960872(1) | 1.003765(1) / 0.986435(1) |
| ν = 3 | -0.324554(1) / 0.569521(1) | 0.56827(1) / 0.896693(4) | 0.891678(1) / 0.980039(1) | 0.96433(1) / 0.99410(1) |
| m²(R−1), ν = 0 | -33.91 / -40.04 | -53.72 / -69.06 | -63.96 / -85.45 | -67.01 / -90.54 |
| m²(R−1), ν = 1 | -2.99 / -33.49 | 7.85 / -53.82 | 15.87 / -63.92 | 18.56 / -66.85 |
| m²(R−1), ν = 3 | -126.84 / -41.22 | -172.63 / -41.31 | -176.96 / -32.61 | -175.81 / -29.10 |

### 5.2 Stress tensors, `l = 2`

#### l = 2, k = 2: R_pres / R_flip (measured/predicted; error in units of the last digit)
| β ≈ 0.4, (Δ,J) → | (14.17, 2) | (24.30, 4) | (44.57, 8) | (74.96, 14) | (105.36, 20) |
|---|---|---|---|---|---|
| m | 12.4 | 22.4 | 42.2 | 72.0 | 101.8 |
| ν = 0 | 0.5725(2) / 0.0510(9546) | 0.8573(3) / -0.0886(20006) | 0.9590(2) / 1.5571(13059) | 0.9859(4) / 2.1842(22622) | 0.9931(7) / 2.9152(35736) |
| ν = 1 | 0.5938(1) / 0.3909(2096) | 0.8679(2) / 0.7750(4382) | 0.9622(1) / 1.1214(2042) | 0.9869(2) / 1.4066(5032) | 0.9933(4) / 1.8772(10889) |
| ν = 3 | 0.8229(5) / 0.3914(1306) | 0.9685(8) / 0.9809(2265) | 0.9940(6) / 1.0192(1338) | 0.9978(13) / 1.2403(2776) | 0.9971(26) / 1.4595(5291) |
| m²(R−1), ν = 0 | -65.91 / -146.33 | -71.30 / -543.81 | -73.07 / 992.91 | -73.08 / 6141.98 | -71.95 / 19855.07 |
| m²(R−1), ν = 1 | -62.63 / -93.91 | -65.99 / -112.42 | -67.31 / 216.41 | -68.11 / 2108.76 | -69.69 / 9093.67 |
| m²(R−1), ν = 3 | -27.30 / -93.84 | -15.73 / -9.55 | -10.71 / 34.25 | -11.46 / 1246.49 | -30.16 / 4763.42 |

| β ≈ 1.0, (Δ,J) → | (13.40, 5) | (24.22, 10) | (45.86, 20) | (80.48, 36) | (115.11, 52) |
|---|---|---|---|---|---|
| m | 10.6 | 20.1 | 39.3 | 70.0 | 100.7 |
| ν = 0 | 0.8212(1) / 0.7101(186) | 0.9426(1) / 0.8945(182) | 0.9841(2) / 0.9793(301) | 0.9950(3) / 1.0217(670) | 0.9977(7) / 1.0523(1104) |
| ν = 1 | 0.88224(2) / 0.7113(14) | 0.96492(2) / 0.8872(19) | 0.99052(4) / 0.9687(25) | 0.9969(2) / 0.9978(80) | 0.9984(4) / 1.0127(198) |
| ν = 3 | 0.5111(6) / 0.8244(36) | 0.8741(6) / 0.9512(21) | 0.9674(9) / 0.9828(46) | 0.9877(30) / 0.9934(118) | 0.9930(67) / 1.0054(198) |
| m²(R−1), ν = 0 | -19.92 / -32.30 | -23.33 / -42.82 | -24.65 / -32.04 | -24.65 / 106.39 | -23.27 / 530.51 |
| m²(R−1), ν = 1 | -13.12 / -32.16 | -14.24 / -45.80 | -14.67 / -48.45 | -15.13 / -10.76 | -16.26 / 129.26 |
| m²(R−1), ν = 3 | -54.46 / -19.56 | -51.10 / -19.83 | -50.52 / -26.63 | -60.50 / -32.23 | -70.92 / 54.97 |

| β ≈ 2.0, (Δ,J) → | (16.60, 11) | (32.36, 23) | (63.87, 47) | (109.83, 82) |
|---|---|---|---|---|
| m | 9.8 | 20.0 | 40.4 | 70.2 |
| ν = 0 | 0.8149(1) / 0.7127(22) | 0.9290(1) / 0.8898(24) | 0.9794(2) / 0.9666(59) | 0.9930(5) / 0.9945(134) |
| ν = 1 | 1.09909(2) / 0.7558(5) | 1.03963(5) / 0.9062(6) | 1.0124(1) / 0.9705(12) | 1.0042(4) / 0.9902(22) |
| ν = 3 | 0.2243(10) / 2.2008(45) | 0.7797(8) / 1.3337(46) | 0.9451(26) / 1.0914(122) | 0.9830(73) / 1.0142(291) |
| m²(R−1), ν = 0 | -17.72 / -27.51 | -28.38 / -44.05 | -33.65 / -54.52 | -34.67 / -27.15 |
| m²(R−1), ν = 1 | 9.49 / -23.39 | 15.85 / -37.50 | 20.33 / -48.22 | 20.72 / -48.41 |
| m²(R−1), ν = 3 | -74.28 / 114.98 | -88.09 / 133.43 | -89.64 / 149.33 | -83.63 / 70.22 |

#### l = 2, k = 3: R_pres / R_flip (measured/predicted; error in units of the last digit)
| β ≈ 0.4, (Δ,J) → | (14.17, 2) | (24.30, 4) | (44.57, 8) | (74.96, 14) | (105.36, 20) |
|---|---|---|---|---|---|
| m | 12.4 | 22.4 | 42.2 | 72.0 | 101.8 |
| ν = 0 | 0.53930(1) / 0.2903(235) | 0.84002(1) / 0.6406(593) | 0.953220(5) / 0.9172(265) | 0.98375(1) / 0.9960(533) | 0.99185(2) / 1.0190(920) |
| ν = 1 | 0.560990(2) / 0.2995(48) | 0.851533(5) / 0.6867(105) | 0.957021(2) / 0.9054(31) | 0.985101(5) / 0.9722(90) | 0.99253(1) / 0.9955(227) |
| ν = 3 | 0.79654(1) / 0.3018(23) | 0.96281(2) / 0.7242(34) | 0.99222(1) / 0.9156(20) | 0.99758(3) / 0.9722(47) | 0.9988(1) / 0.9944(98) |
| m²(R−1), ν = 0 | -71.03 / -109.43 | -79.92 / -179.53 | -83.38 / -147.63 | -84.31 / -20.84 | -84.52 / 196.80 |
| m²(R−1), ν = 1 | -67.69 / -108.00 | -74.17 / -156.53 | -76.60 / -168.57 | -77.28 / -144.37 | -77.49 / -46.41 |
| m²(R−1), ν = 3 | -31.37 / -107.66 | -18.58 / -137.78 | -13.87 / -150.50 | -12.57 / -144.31 | -12.40 / -58.32 |

| β ≈ 1.0, (Δ,J) → | (13.40, 5) | (24.22, 10) | (45.86, 20) | (80.48, 36) | (115.11, 52) |
|---|---|---|---|---|---|
| m | 10.6 | 20.1 | 39.3 | 70.0 | 100.7 |
| ν = 0 | 0.742420(2) / 0.5596(3) | 0.911649(1) / 0.8209(3) | 0.974839(2) / 0.9456(5) | 0.99189(1) / 0.9826(15) | 0.99606(3) / 0.9925(31) |
| ν = 1 | 0.810536(1) / 0.57261(3) | 0.939613(1) / 0.82977(3) | 0.983322(1) / 0.94866(5) | 0.994666(4) / 0.9833(2) | 0.99741(1) / 0.9921(4) |
| ν = 3 | 0.34378(1) / 0.72318(4) | 0.81305(1) / 0.91125(2) | 0.95119(1) / 0.9756(1) | 0.9846(1) / 0.9921(2) | 0.9924(2) / 0.9961(5) |
| m²(R−1), ν = 0 | -28.69 / -49.06 | -35.87 / -72.72 | -38.94 / -84.24 | -39.81 / -85.31 | -39.96 / -75.89 |
| m²(R−1), ν = 1 | -21.11 / -47.61 | -24.52 / -69.12 | -25.81 / -79.46 | -26.17 / -82.14 | -26.29 / -80.50 |
| m²(R−1), ν = 3 | -73.10 / -30.84 | -75.91 / -36.03 | -75.54 / -37.80 | -75.58 / -38.87 | -76.81 / -39.31 |

| β ≈ 2.0, (Δ,J) → | (16.60, 11) | (32.36, 23) | (63.87, 47) | (109.83, 82) |
|---|---|---|---|---|
| m | 9.8 | 20.0 | 40.4 | 70.2 |
| ν = 0 | 0.665678(2) / 0.51188(4) | 0.866782(2) / 0.79874(4) | 0.96010(1) / 0.9356(2) | 0.98603(2) / 0.9771(4) |
| ν = 1 | 1.017503(1) / 0.56474(1) | 1.024651(2) / 0.82521(1) | 1.010588(3) / 0.94549(3) | 1.00404(1) / 0.9807(1) |
| ν = 3 | -0.36374(2) / 2.2924(1) | 0.56371(1) / 1.4506(1) | 0.88672(5) / 1.1309(3) | 0.9621(2) / 1.0453(7) |
| m²(R−1), ν = 0 | -32.01 / -46.74 | -53.27 / -80.48 | -65.19 / -105.20 | -68.83 / -112.62 |
| m²(R−1), ν = 1 | 1.68 / -41.68 | 9.86 / -69.89 | 17.30 / -89.05 | 19.90 / -95.26 |
| m²(R−1), ν = 3 | -130.59 / 123.75 | -174.46 / 180.19 | -185.06 / 213.90 | -186.69 / 223.37 |

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
`m^{−k−1/2}` term has a fitted coefficient consistent with zero at k = 2 (`−0.5 ± 2.6` (first-version budget; to be refreshed), `l = 1`, β ≈ 0.4, ν = 0, pres), i.e. limited to
`|c| ≲ 6 ≈ 0.35|c_2|` at 2σ — **an `m^{−5/2}` term with an O(1) coefficient is not excluded by the data**. At k = 3 the fit finds `c_{7/2} = 39 ± 0.7`
(statistical), but the three-parameter fits with and without it are equally good (χ² 0.09 vs 0.3, one degree of freedom), so it cannot be told
apart from the truncation of the analytic series. `m^{−2k−1} = m^{−5}` at k = 2 is invisible behind `c_4/m⁴ ≈ 300/m⁴` and is excluded only by the
sketch. The even-power ansätze (B2, B3) rest on E2E's parity argument: the pointwise `1/m^n` terms of `b𝒯` are odd under `1 ↔ 3` for odd `n`
(reflection `m → −m`), so they cancel in the moments by Lemma S; with the uniformity above this holds order by order.
**Answer to the task question:** the term-by-term expansion does not fail at k = 2 at order `1/m²` [numerical evidence: `m²(R−1)` converges];
at `1/m³`, `m^{−5/2}` and `m^{−5}` the absence of non-analytic terms rests on the derivation sketch [derived (sketch)], with the data limiting
an `m^{−5/2}` coefficient only to `|c| ≲ 6`.

## 7. Corrections at k = 2, 3 (`t2_corrections.py`, `out_t2_corrections.txt`) [fitted]

Weighted fits on `m ≥ 20` (four masses at β ≈ 0.4, 1; three at β ≈ 2, then exact). `c_2` from `c_2/m² + c_4/m⁴` (B2) and with a `c_6/m⁶` term (B3);
`c_1` from `c_1/m + c_2/m² + c_4/m⁴` (A), and `c_2` with an added `m^{−k−½}` (N). Adopted `c_2`: mid-range of the fits with χ²/dof < 4, ± half the
range (the first version quoted the B2/B3 spread although B2 is rejected at k = 3, χ² ≈ 3·10³ — referee M2). Statistical errors are given only
for `c_1` and `c_N`; they ignore the truncation of the ansatz.

| l | β | k | ν | comp | c_2 (B2) | c_2 (B3) | c_1 (A) | c_{k+½} (N) |
|---|---|---|---|---|---|---|---|---|
| 1 | 0.4 | 2 | 0.0 | flip | -35.451±1.6e-01 | -34.783±1.6e-01 | +0.0167±7.4e-02 | -18.979±8.6e+01 |
| 1 | 0.4 | 2 | 1.0 | pres | -13.711±1.1e-03 | -13.697±1.1e-03 | +0.0004±1.7e-03 | -0.420±2.0e+00 |
| 1 | 0.4 | 2 | 1.0 | flip | -32.971±2.4e-02 | -32.772±2.4e-02 | +0.0051±4.0e-02 | -5.684±4.7e+01 |
| 1 | 0.4 | 2 | 3.0 | pres | 14.836±1.7e-02 | 14.787±1.7e-02 | -0.0013±2.3e-02 | +1.437±2.6e+01 |
| 1 | 0.4 | 2 | 3.0 | flip | -12.834±8.4e-02 | -12.822±8.4e-02 | +0.0006±9.3e-02 | -0.486±1.1e+02 |
| 1 | 0.4 | 3 | 0.0 | pres | -27.864±5.9e-06 | -27.899±5.9e-06 | -0.0008±1.5e-05 | +39.120±6.9e-01 |
| 1 | 0.4 | 3 | 0.0 | flip | -54.254±2.0e-04 | -54.450±2.0e-04 | -0.0047±2.2e-04 | +223.634±1.1e+01 |
| 1 | 0.4 | 3 | 1.0 | pres | -23.619±4.3e-06 | -23.638±4.3e-06 | -0.0005±1.1e-05 | +21.699±5.2e-01 |
| 1 | 0.4 | 3 | 1.0 | flip | -51.074±4.4e-05 | -51.238±4.4e-05 | -0.0040±1.6e-04 | +187.272±7.6e+00 |
| 1 | 0.4 | 3 | 3.0 | pres | 12.531±1.7e-05 | 12.551±1.7e-05 | +0.0005±6.4e-05 | -23.228±3.0e+00 |
| 1 | 0.4 | 3 | 3.0 | flip | -25.798±6.7e-05 | -25.786±6.7e-05 | +0.0003±2.2e-04 | -13.762±1.0e+01 |
| 1 | 1.0 | 2 | 0.0 | pres | -15.788±2.8e-03 | -15.730±2.8e-03 | +0.0015±5.5e-03 | -1.523±5.9e+00 |
| 1 | 1.0 | 2 | 0.0 | flip | -29.364±2.6e-02 | -29.060±2.6e-02 | +0.0077±3.2e-02 | -7.952±3.4e+01 |
| 1 | 1.0 | 2 | 1.0 | pres | -8.198±1.9e-03 | -8.156±1.9e-03 | +0.0010±7.3e-05 | -1.085±4.7e+00 |
| 1 | 1.0 | 2 | 1.0 | flip | -24.945±1.0e-02 | -24.681±1.0e-02 | +0.0067±2.2e-02 | -6.923±2.4e+01 |
| 1 | 1.0 | 2 | 3.0 | pres | -18.682±6.6e-02 | -18.081±6.6e-02 | +0.0149±2.5e-03 | -15.573±1.6e+02 |
| 1 | 1.0 | 2 | 3.0 | flip | 10.787±4.1e-02 | 10.935±4.1e-02 | +0.0036±8.5e-02 | -3.788±9.1e+01 |
| 1 | 1.0 | 3 | 0.0 | pres | -30.663±7.8e-06 | -30.797±7.8e-06 | -0.0032±4.2e-05 | +130.557±1.7e+00 |
| 1 | 1.0 | 3 | 0.0 | flip | -52.440±3.6e-05 | -52.870±3.6e-05 | -0.0104±1.4e-04 | +420.256±5.7e+00 |
| 1 | 1.0 | 3 | 1.0 | pres | -19.724±6.0e-06 | -19.750±6.0e-06 | -0.0006±3.1e-05 | +25.419±1.3e+00 |
| 1 | 1.0 | 3 | 1.0 | flip | -46.462±2.0e-05 | -46.774±2.0e-05 | -0.0075±1.1e-04 | +305.884±4.5e+00 |
| 1 | 1.0 | 3 | 3.0 | pres | -45.579±6.3e-05 | -45.316±6.3e-05 | +0.0064±2.4e-06 | -258.073±6.2e-03 |
| 1 | 1.0 | 3 | 3.0 | flip | 3.903±3.6e-05 | 4.067±3.6e-05 | +0.0040±1.4e-06 | -160.612±3.5e-03 |
| 1 | 2.0 | 2 | 0.0 | pres | -33.445±1.1e-02 | -33.656±1.1e-02 | -0.0051±1.6e-02 | +5.520±1.7e+01 |
| 1 | 2.0 | 2 | 0.0 | flip | -47.843±3.9e-02 | -47.988±3.9e-02 | -0.0035±4.1e-02 | +3.794±4.4e+01 |
| 1 | 2.0 | 2 | 1.0 | pres | 20.648±2.0e-02 | 21.365±2.0e-02 | +0.0174±3.2e-02 | -18.695±3.4e+01 |
| 1 | 2.0 | 2 | 1.0 | flip | -31.794±2.6e-02 | -31.137±2.6e-02 | +0.0159±3.5e-02 | -17.133±3.8e+01 |
| 1 | 2.0 | 2 | 3.0 | pres | -75.655±4.9e-01 | -70.155±4.9e-01 | +0.1333±7.2e-01 | -143.417±7.8e+02 |
| 1 | 2.0 | 2 | 3.0 | flip | 6.055±3.1e-01 | 11.170±3.1e-01 | +0.1240±4.4e-01 | -133.386±4.7e+02 |
| 1 | 2.0 | 3 | 0.0 | pres | -67.285±7.9e-05 | -68.630±7.9e-05 | -0.0326±1.3e-04 | +1350.713±5.3e+00 |
| 1 | 2.0 | 3 | 0.0 | flip | -90.767±2.6e-04 | -93.273±2.6e-04 | -0.0608±2.3e-04 | +2518.857±9.5e+00 |
| 1 | 2.0 | 3 | 1.0 | pres | 18.475±1.4e-04 | 20.012±1.4e-04 | +0.0372±2.1e-04 | -1543.390±8.5e+00 |
| 1 | 2.0 | 3 | 1.0 | flip | -67.193±9.1e-05 | -68.398±9.1e-05 | -0.0292±2.0e-04 | +1209.439±8.4e+00 |
| 1 | 2.0 | 3 | 3.0 | pres | -178.364±1.3e-03 | -174.955±1.3e-03 | +0.0826±2.4e-03 | -3424.075±9.8e+01 |
| 1 | 2.0 | 3 | 3.0 | flip | -29.790±6.2e-04 | -27.137±6.2e-04 | +0.0643±1.3e-03 | -2664.439±5.6e+01 |
| 2 | 0.4 | 2 | 0.0 | pres | -73.730±4.4e-01 | -72.847±2.8e+00 | +0.0217±6.2e-02 | -25.621±7.6e+01 |
| 2 | 0.4 | 2 | 0.0 | flip | 2025.982±3.1e+03 | 11918.968±1.8e+04 | +242.1542±4.0e+02 | -287945.261±5.0e+05 |
| 2 | 0.4 | 2 | 1.0 | pres | -67.843±1.9e-01 | -68.774±1.4e+00 | -0.0230±3.1e-02 | +26.854±3.8e+01 |
| 2 | 0.4 | 2 | 1.0 | flip | 421.695±5.0e+02 | 4091.264±4.2e+03 | +91.5983±9.6e+01 | -105965.415±1.2e+05 |
| 2 | 0.4 | 2 | 3.0 | pres | -8.904±1.4e+00 | -14.252±1.0e+01 | -0.1471±2.3e-01 | +162.698±2.8e+02 |
| 2 | 0.4 | 2 | 3.0 | flip | 122.736±3.3e+02 | 2506.709±2.3e+03 | +58.4569±5.2e+01 | -68704.447±6.3e+04 |
| 2 | 0.4 | 3 | 0.0 | pres | -84.728±1.2e-02 | -84.797±9.6e-02 | -0.0015±2.2e-03 | +78.272±1.1e+02 |
| 2 | 0.4 | 3 | 0.0 | flip | -127.877±6.5e+01 | 93.785±4.4e+02 | +5.2220±9.7e+00 | -259069.149±5.0e+05 |
| 2 | 0.4 | 3 | 1.0 | pres | -77.555±4.3e-03 | -77.637±4.0e-02 | -0.0019±9.1e-04 | +94.512±4.5e+01 |
| 2 | 0.4 | 3 | 1.0 | flip | -172.454±7.7e+00 | -119.860±7.6e+01 | +1.3003±1.8e+00 | -61183.689±8.6e+04 |
| 2 | 0.4 | 3 | 3.0 | pres | -12.039±2.6e-02 | -11.886±2.1e-01 | +0.0033±4.9e-03 | -172.857±2.4e+02 |
| 2 | 0.4 | 3 | 3.0 | flip | -154.932±4.7e+00 | -129.680±3.9e+01 | +0.6876±8.9e-01 | -30373.501±4.4e+04 |
| 2 | 1.0 | 2 | 0.0 | pres | -25.108±2.3e-01 | -24.481±2.3e-01 | +0.0159±5.0e-02 | -16.731±5.5e+01 |
| 2 | 1.0 | 2 | 0.0 | flip | -22.117±6.2e+01 | 247.619±6.2e+01 | +6.5992±1.1e+01 | -7123.980±1.2e+04 |
| 2 | 1.0 | 2 | 1.0 | pres | -14.827±7.8e-02 | -15.456±7.8e-02 | -0.0154±2.4e-02 | +16.353±2.7e+01 |
| 2 | 1.0 | 2 | 1.0 | flip | -48.678±5.1e+00 | 21.105±5.1e+00 | +1.7297±1.4e+00 | -1830.057±1.5e+03 |
| 2 | 1.0 | 2 | 3.0 | pres | -50.468±1.8e+00 | -67.066±1.8e+00 | -0.3855±5.0e-01 | +422.757±5.5e+02 |
| 2 | 1.0 | 2 | 3.0 | flip | -28.973±9.5e+00 | -21.932±9.5e+00 | +0.3092±1.9e+00 | -256.004±2.2e+03 |
| 2 | 1.0 | 3 | 0.0 | pres | -40.034±5.2e-03 | -40.220±5.2e-03 | -0.0042±1.7e-03 | +181.157±7.1e+01 |
| 2 | 1.0 | 3 | 0.0 | flip | -88.298±9.4e-01 | -84.602±9.4e-01 | +0.0955±2.5e-01 | -3729.906±1.1e+04 |
| 2 | 1.0 | 3 | 1.0 | pres | -26.271±1.3e-03 | -26.340±1.3e-03 | -0.0016±6.6e-04 | +66.892±2.7e+01 |
| 2 | 1.0 | 3 | 1.0 | flip | -83.136±9.9e-02 | -83.277±9.9e-02 | -0.0015±2.4e-02 | +122.188±1.0e+03 |
| 2 | 1.0 | 3 | 3.0 | pres | -75.406±2.3e-02 | -75.675±2.3e-02 | -0.0071±1.0e-02 | +270.307±4.2e+02 |
| 2 | 1.0 | 3 | 3.0 | flip | -38.441±1.1e-01 | -39.463±1.1e-01 | -0.0233±3.3e-02 | +997.126±1.4e+03 |
| 2 | 2.0 | 2 | 0.0 | pres | -35.360±3.6e-01 | -35.156±3.6e-01 | +0.0050±9.1e-02 | -5.398±9.8e+01 |
| 2 | 2.0 | 2 | 0.0 | flip | -56.957±1.2e+01 | -9.661±1.2e+01 | +1.1689±2.6e+00 | -1257.524±2.8e+03 |
| 2 | 2.0 | 2 | 1.0 | pres | 21.774±2.0e-01 | 20.837±2.0e-01 | -0.0229±7.2e-02 | +24.596±7.7e+01 |
| 2 | 2.0 | 2 | 1.0 | flip | -51.590±2.5e+00 | -48.227±2.5e+00 | +0.0840±4.4e-01 | -90.422±4.7e+02 |
| 2 | 2.0 | 2 | 3.0 | pres | -90.006±5.6e+00 | -79.807±5.6e+00 | +0.2505±1.4e+00 | -269.492±1.5e+03 |
| 2 | 2.0 | 2 | 3.0 | flip | 152.009±2.6e+01 | 20.082±2.6e+01 | -3.2556±5.7e+00 | +3502.602±6.2e+03 |
| 2 | 2.0 | 3 | 0.0 | pres | -69.062±1.2e-02 | -70.782±1.2e-02 | -0.0420±4.0e-03 | +1741.487±1.7e+02 |
| 2 | 2.0 | 3 | 0.0 | flip | -113.249±2.7e-01 | -116.579±2.7e-01 | -0.0816±7.6e-02 | +3383.159±3.1e+03 |
| 2 | 2.0 | 3 | 1.0 | pres | 19.726±7.5e-03 | 21.328±7.5e-03 | +0.0392±2.2e-03 | -1626.281±8.9e+01 |
| 2 | 2.0 | 3 | 1.0 | flip | -95.316±6.6e-02 | -98.608±6.6e-02 | -0.0811±1.5e-02 | +3363.983±6.2e+02 |
| 2 | 2.0 | 3 | 3.0 | pres | -188.486±9.9e-02 | -187.412±9.9e-02 | +0.0262±3.7e-02 | -1085.618±1.5e+03 |
| 2 | 2.0 | 3 | 3.0 | flip | 224.872±4.7e-01 | 228.346±4.7e-01 | +0.0852±1.3e-01 | -3531.325±5.3e+03 |

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
