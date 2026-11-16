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
   `≤ 4·10⁻¹³` where the latter converges; Casimir residual `≤ 10⁻⁷` (relative to the largest term) for `1 − ρ_b ≥ 10⁻⁴` (§2).
2. **Currents (`l = 1`), k = 2 and 3: converged and error-controlled [numerical evidence for the ratios; the formula itself is the derived
   one of E2E §4.5, nothing fitted].** At the largest masses (`m ≈ 100`; β ≈ 0.4 / 1.0, ν = 0):
   * k = 3: `R_pres = 0.997317(2) / 0.99698(1)`, `R_flip = 0.99478(2) / 0.99483(2)`; `m²(R−1)` converges like `c_2 + c_4/m²`
     (e.g. β ≈ 0.4, ν = 0: `−23.68, −26.39, −27.45, −27.74, −27.81` (pres), `−40.5, −49.1, −52.8, −53.9, −54.1` (flip), `m = 12 … 102`).
   * k = 2: `R_pres = 0.9984(2) / 0.9986(3)`, `R_flip = 0.9972(30) / 0.9976(16)`; `m²(R−1)` is converged within errors from `m ≈ 40`
     (β ≈ 0.4, ν = 0: `−15.36, −16.49, −16.92, −17.02(6), −16.65(120)` (pres), `−27.7, −32.6, −34.7, −34.8(18), −29(30)` (flip)).
   * all three β (0.4, 1, 2) and ν = 0, 1, 3: the full tables are §5; at β ≈ 2 (only `m ≤ 70`) `m²(R−1)` is still drifting for ν = 3.
3. **Stress tensors (`l = 2`) [partial].** @@L2SUMMARY@@
4. **No `O(1/m)` term at k = 2, 3 [numerical evidence; the argument of E2E §6.1 now holds without Assumption T, see 5].** Fitted `c_1`
   (ansatz `c_1/m + c_2/m² + c_4/m⁴`, `m ≥ 20`): `|c_1| ≤ 0.017` (k = 2) and `≤ 0.011` (k = 3) for β ≤ 1 (`l = 1`), with statistical errors of
   the same size at k = 2.
5. **Non-analytic terms / Assumption T [derived (sketch) + numerical evidence].** Through the exact identity (★) of E2E §3.1 the moment
   is a 6-dimensional integral of the elementary shadow integrand; in heavy variables its expansion in `1/m` is uniform in all directions
   (every Cartesian component of `X = mx` is bounded by `e^{β/2} ζ`, `ζ = −P̂·X`, which the exponential `e^{−2ζ}` controls), so the moment has
   an asymptotic expansion in powers of `1/m` with no non-analytic terms — neither `m^{−k−1/2}` (breakdown scale `u ~ 1/m`) nor
   `m^{−2k−1}` (`u ~ 1/m²`) — for every k for which (★) converges (here all k > −1). The non-uniformity at `u → 0` that motivated
   Assumption T is a property of the `(r, η)` representation only. Numerically: the pointwise `1/m²` coefficient of the shadow integrand
   stays bounded (it saturates) when the hyperbolic distances of `x̂, ŷ` from `p̂, p̄̂` grow to 8 (`d(x̂, ŷ)` up to 12.6), `l = 1, 2` (§6);
   fits of the completed k = 2 ratios with an added `m^{−5/2}` term give a coefficient consistent with zero (`−0.5 ± 2.6` against
   `c_2 = −17.1`, `l = 1`, β ≈ 0.4, ν = 0). So the term-by-term expansion does **not** fail at k = 2 at order `1/m²` or `1/m³`; the only
   evidence for that is the (sketch) derivation plus the absence of any detectable drift — a rigorous remainder estimate is not written.
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
`H^{(l)}_i(x) = A(x) + (1−x)^{−2σ−l}B(x)` near 1, `(1+x)^{2c_{1,2}−l}`-type at −1 (logarithmic if `c1 = c2`); [inferred: that the Frobenius
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
| Casimir residual / largest term, deep corner, A and B | `(115.11,52,.625,−.375), (105.36,20,.375,.375), (14.17,2,.625,−.375)` | `1−ρ_b ≈ 3·10⁻⁶ (t = 1): ≤ 10⁻⁴`; `1−ρ_b ≥ 10⁻⁴: ≤ 7·10⁻⁸`; `≥ 10⁻³: ≤ 7·10⁻¹⁰` |
| `|S_A/S_B − 1|`, same points | same | `≤ 9·10⁻⁸` (deepest), `≤ 10⁻¹¹` for `1−ρ_b ≥ 10⁻³` |

The deep-corner A/B difference (`4–9·10⁻⁸`) is of the size of the deviation of A from the exact collinear block (`3–8·10⁻⁸`), which supports
`|A − B|` as error proxy; a safety factor 3 is used (§3). [The radial reference with `N = 24000` costs `2.3 GB`; the LC representation of one
seed with setting A takes ~40 s.]

## 3. Corner completion of the moment grids (`corner.py`) and error budget

**Grid.** The e2e grid (72 × 62 Gauss nodes in `(u, t)`) under-resolves the corner at k = 2: below `u = 0.01` it has 8 nodes, too few for the
log-oscillating kernel `P_{−½+iν}(cosh ξ)`, `ξ = −log u`, at ν = 3 (§4.3). The results use the **corner-refined grid** (`refined_base`): the two
u-panels below 0.03 (18 nodes) are replaced by six panels with 46 nodes from `u = 3·10⁻⁶`; above `u = 0.03` the nodes and the e2e radial-series
tensors are kept. All corner nodes, and the e2e points that were dropped, are evaluated with A and B through the *unchanged* spinning machinery
(`spinblock.SpinBlock.eval`: differential operators, `(P,Q)` jets, generated tensor code, Fourier modes); 300 random converged e2e points with
`u < 0.25` are re-evaluated with A as a check.

**Per point.** With `δ = max_modes|A−B| / max_modes|A|`: the point is accepted if `δ < 10⁻²` **and** it is connected, within its t-column, to the
radial-series region through accepted points whose size changes by less than `u^{±8}` between neighbours (isolated points inside a rejected
region are where A and B are both corrupted). Its error is `3|A − B|` (summed in absolute value). Rejected points (deepest corner: the conserved
spinning combination cancels between the seeds by up to ~`(1−ρ_b)^{−(2σ+4)}`, beyond double precision) and points where the lightcone sum has not
converged (tiny `t`, where also the radial series fails) are **estimated** by a power law through the two nearest accepted nodes of their u-row
(tiny t) or t-column (deep corner); the estimate is added to the moment and **50 % of it** is counted as error. [The 50 % is a choice, not derived;
the estimated share is listed for every entry (`est` in `out_t1_lowk_e2e_l*.txt`).] Radial-series points: e2e criterion, `≤ 10⁻⁹`.

**Per ratio.** The relative errors of the modes `(l,−l)`, `(l,l)` are propagated through the 2×2 helicity solve (E2E §5.2) to `R_pres`, `R_flip`.
Validation of the re-evaluated converged points: `max dev ≤ 2·10⁻⁹` (`l = 1`, all cases), @@L2VALID@@ (`l = 2`).

## 4. Validation of the completed moments

### 4.1 Scalar partial waves (`v1_scalar.py`, `out_v1_scalar.txt`)

@@SCALAR@@

### 4.2 k = 4, 6 against the e2e numbers [numerical evidence]

`l = 1`: where E2E's numbers were flagged converged (β ≤ 1, k = 6, all m; k = 4, m ≤ 42) the completed values agree with them to
`≤ 3.5·10⁻⁵` (k = 6, flip, β ≈ 1, m ≈ 100; typically `10⁻⁶`); the larger differences sit exactly at E2E's `*`/`**` entries
(k = 4, `m ≳ 70`: up to `5.9·10⁻³`; β ≈ 2, m ≈ 70: `1.7·10⁻⁴` at k = 6, `8.4·10⁻³` at k = 4) and are the dropped corner that E2E's `tail` estimated
(E2E: k = 4 tail `2·10⁻⁴ … 2·10⁻³` at the largest m). @@L2CONTROL@@

### 4.3 Quadrature (`q1_grid.py`, `out_q1_grid.txt`)

@@GRID@@

## 5. Results (`out_t1_lowk_e2e_l{1,2}.txt`, `out_t3_tables.txt`)

Measured = tensor Regge moments of the exact partial wave (corner-completed); predicted = E2E §4.5 with the exact `b^{(1−Δ_φ)}`; nothing fitted.
Errors in units of the last digit (the full numerical budget of §3; quadrature not included, see §4.3).

### 5.1 Currents, `l = 1`

@@L1TABLES@@

### 5.2 Stress tensors, `l = 2`

@@L2TABLES@@

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
| 1 | 0, 2, 4, 6, 8 | 0 | 0.85 … 8.64 | 7.47, 8.50, 8.45, 8.45, 8.45 |
| 1 | 0 | 4, 8 | 4.58, 8.58 | 3.73, 3.76 |
| 1 | 6 | 6 | 12.6 | 11.88 |
| 2 | 0, 2, 4, 6, 8 | 0 | 0.85 … 8.64 | 19.27, 18.90, 18.99, 18.99, 18.99 |
| 2 | 0 | 4, 8 | 4.58, 8.58 | 21.07, 21.35 |
| 2 | 6 | 6 | 12.6 | 76.7 |

The coefficient saturates as the distances grow (no exponential growth), as (b)–(c) predict. **Fits of the completed ratios** (§7): an added
`m^{−k−1/2}` term has a fitted coefficient consistent with zero at k = 2 (`−0.5 ± 2.6`, `l = 1`, β ≈ 0.4, ν = 0, pres; `|c_N| ≲ 0.3|c_2|`); at k = 3 the
three-parameter fits with and without it are equally good (χ² 0.09 vs 0.3, one degree of freedom) — the data cannot separate `m^{−7/2}` from
the truncation of the analytic series, i.e. there is no evidence for it either way beyond the derivation. **Answer to the task question:** the
term-by-term expansion does not fail at k = 2 at order `1/m²` [numerical evidence: `m²(R−1)` converges], nor at `1/m³` [derived (sketch); the
`m^{−5/2}` coefficient is bounded only to `≲ 5` by the data]; no non-analytic term was found or is expected.

## 7. Corrections at k = 2, 3 (`t2_corrections.py`, `out_t2_corrections.txt`) [fitted]

Weighted fits on `m ≥ 20` (four masses at β ≈ 0.4, 1; three at β ≈ 2, then exact). `c_2` from `c_2/m² + c_4/m⁴` (B2) and with a `c_6/m⁶` term (B3);
`c_1` from `c_1/m + c_2/m² + c_4/m⁴` (A). Statistical errors only; the B2/B3 spread is the better uncertainty.

@@CORRTABLE@@

## 8. What remains

@@REMAINS@@

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
