# Spinning tensor Regge moments at k = 2, 3: corner completion of the exact d = 3 partial waves and the end-to-end test

Status: written 16 Nov 2026 by a numerics/derivation agent. Unrefereed. Code and recorded outputs: `ads4/code/lowk/` (index §9; `sh run_all.sh`
regenerates every number; it needs the e2e block-tensor cache `ads4/code/e2e/data_modes/`). This note continues `SPINNING_E2E.md`
("E2E"; conventions §1 there, notation unchanged) and closes its §11 item 1 and the `**` flags at k = 2.

Marks: **[derived]**, **[derived+verified]**, **[numerical evidence]**, **[inferred]**. "Verified" means an independent numerical check
of the derived statement; "numerical evidence" means read off from numbers without derivation. No constant is fitted anywhere in a
statement marked derived; fitted coefficients (§7) are labelled as such.

@@SUMMARY@@

## 1. The problem and the choice of method

E2E §3.2: the moment integrand in the variables `(u, t)` (`√ρ = s = t/(2(α+ᾱu))`, `√ρ̄ = su`, `u = e^{−ξ}`) behaves like `u^{k−1/2}` for
`u → 0`, and there the reflected scalar seed blocks are needed at `ρ_b := ((1−su)/(1+su))² → 1` (`z̄' → 1` of the reflected block, i.e.
`z̄ → −∞` at fixed `z`), where the radial series converges like `ρ_b^N = e^{−4Nsu}`. With `N ≤ 10⁴` the e2e runs drop the grid points with
`u ≲ 30α/(Nt)`; at k = 2 the dropped region carries up to 10⁻¹ of a mode.

Options considered (task list): (i) a crossed-lightcone expansion of the d = 3 block at `z̄ = ∞` needs the connection coefficients of the
`z̄' → 1` singular term, which are not known in closed form for d = 3; (iii) the d = 3 → d = 2 dimensional reduction converges only like
`ρ_s^{n/2}` and needs its coefficients computed by an O(N²J²) matching; (iv) per-point acceleration of the level sums `Σ_n L_n` mixes the
`ρ_s`- and `ρ_b`-dependence. **Chosen:** a hybrid of (i) and (iv) that is exact in structure and needs no new connection data:
reorganise the radial series into the lightcone (monomial) basis `ρ_s^i ρ_b^q`; at fixed `i` the `q`-series is a one-variable power series
whose singularities on the unit circle and their exponents are known from the Casimir equation, so its tail beyond `q = Q` follows from the
Darboux (transfer) expansion; its coefficients are fitted on `q ∈ [Q/2, Q]` and the tail is summed in closed form (incomplete beta function).
Nothing about the physics formula is used; the representation is tested against independent evaluations (§2.3).

## 2. Scalar seed blocks near `ρ_b → 1` (`lcblock.py`, `lcconv.c`)

### 2.1 Lightcone reorganisation [derived]

With `G̃(r, η) = Σ_{n,j} B_{n,j} r^n P_j(η)` (E2E §2.1, `g = (4r)^Δ G̃`, reflected block `g^{c1,c2}`, `c1 + c2 > 0` so that `B_{n,j} ≥ 0`),
`ρ_s = re^{−χ}`, `ρ_b = re^{χ}` (`η = cosh χ`) and `P_j(cos θ) = Σ_{k=0}^{j} a_k a_{j−k} e^{i(j−2k)θ}`, `a_k = (½)_k/k!`:

`G̃ = (ρ_sρ_b)^{−J/2} Σ_{i,q≥0} C_{i,q} ρ_s^i ρ_b^q`,   `C_{i,q} = Σ_{k'=0}^{min(i,q)} a_{k'} a_{|q−i|+k'} B_{i+q−J, |q−i|+2k'}`

(an exact rearrangement of a series of non-negative terms; `C_{i,q} = C_{q,i}`). Checks: reconstruction of `G̃` at a generic point
`1·10⁻¹⁵`; symmetry exact; `H_0(x) := Σ_q C_{0,q} x^q` equals the collinear block `4^{−(Δ+J)/2} x^{−(Δ−J)/2} k_{Δ+J}(4x/(1+x)²)`,
`k_β(z) = z^{β/2} ₂F₁(β/2+c1, β/2+c2; β; z)`, to `10⁻¹⁵` where the series converges.

### 2.2 Tails in `q` [derived; the joint-analyticity input is inferred]

At fixed `i`, `H_i(x) = Σ_q C_{i,q} x^q` has on `|x| ≤ 1` only the singularities `x = 1` (`z̄' = 1`) and `x = −1` (`z̄' = ∞`). The local
exponents of the Casimir equation are `{0, −σ}`, `σ = c1 + c2`, at `z̄' = 1` and `{−c1, −c2}` at `z̄' = ∞` (the `(d−2)` term is subleading
in both indicial equations; checked by hand in §2 of the code docstring). With `1 − z̄' = ((1−x)/(1+x))²` this gives
`H_i(x) = A_i(x) + (1−x)^{−2σ} B_i(x)` near `x = 1` with `A_i, B_i` analytic, and `(1+x)^{2c_{1,2}}`-type behaviour at `x = −1` (logarithmic if
`c1 = c2`). [inferred: that the expansion of the two Frobenius pieces in powers of `ρ_s` preserves this form term by term; tested by (a)–(c)
of §2.3.] For `Δ_φ = 1.75` every `2σ` is a half-integer, so there are no logarithms at `x = 1`. Darboux: the Taylor coefficients of
`(1−x)^{m−2σ}` are `(2σ−m)_q/q!`; hence, for `t = 0…4` (the `ρ_b`-derivatives needed by the jets),

`q^t C_{i,q} ≃ Σ_{m<M} s^{(t)}_{i,m} (2σ+t−m)_q/q! + (−1)^q Σ_{c∈{c1,c2}} Σ_{m<3} a^{(t)}_{i,c,m} (−2c+t−m)_q/q!`   (`q → ∞`).

The coefficients are fitted on `q ∈ [Q/2, Q]` (every third `q`, one least-squares problem per `t` with `I+1` right-hand sides); the series is
summed exactly up to `Q` and the tails `Σ_{q>Q}(a)_q/q! x^q` in closed form: `(1−x)^{−a} I_x(Q+1, a)` for `a > 0` (regularised incomplete beta:
the partial sums are negative-binomial distribution functions), upward recursion in `a` for `a < 0`, and a Pfaff-transformed ₂F₁ for the
alternating family (checked against mpmath ₂F₁: `≤ 10⁻¹¹`). The `ρ_s`-sum is truncated at `i = I`; points where its last term exceeds
`10⁻¹³` of the sum (tiny `t`) are not used. Derivative sums: `S[k,l] = r^{−Δ}(θ_s+θ_b)^k ∂_η^l[r^ΔG̃]` with `∂_η = (sinh χ)^{−1}(θ_b − θ_s)`
(free of cancellations where `1 − ρ_b ≪ χ`, i.e. in the corner; not used near `u ~ 1`, where `χ → 0`); they are fed unchanged into
`spinblock.pq_jets` and the generated tensor code of E2E. Two settings: **A** `(N, I, Q, M) = (10⁴, 1500, 8000, 8)`, **B** `(8000, 1200, 6000, 6)`
(B uses the same `B_{n,j}`, which do not depend on the truncation, but a different window, depth and number of asymptotic terms).

### 2.3 Validation of the block representation (`c1_block_checks.py`, `out_c1_block_checks.txt`) [derived+verified]

@@BLOCKCHECKS@@

## 3. Corner completion of the moment grids (`corner.py`) and error budget

The e2e grids (72 × 62 Gauss nodes in `(u, t)`, cached tensors at the converged points) are kept. Every dropped point is evaluated with A and
B through the unchanged spinning machinery (differential operators, `(P,Q)` jets, tensor code, Fourier modes). Per point, the relative
difference `|A−B|/max_modes|A|` is computed:

* accepted if `< 0.1`; then `3|A − B|` is that point's error (summed in absolute value over points, no cancellation assumed). The factor 3
  is a safety factor: at block level `|A − B|` tracks the deviation of A from a `Q = 14000` reference within a factor ~2;
* otherwise (deepest corner, where the conserved spinning combination cancels between seeds by many orders of magnitude and double
  precision is insufficient) and at tiny `t` (lightcone sum not converged): the point is replaced by a power-law extrapolation in `u` from the
  two smallest accepted nodes of its `t`-column (a `t^{a−1}` law for columns without accepted nodes), and **the whole replaced contribution
  is counted as error** ("drop");
* radial-series points: the e2e criterion `X − 10 log X > 25` bounds their truncation error by `~10⁻⁹`;
* re-evaluation of 300 random converged points with `u < 0.25` with A reproduces the cached radial values to `@@VALID@@` (all cases).

The errors of the two transverse modes are propagated through the 2×2 helicity solve (E2E §5.2) to `R_pres`, `R_flip` (`t1_lowk_e2e.py`).
The quadrature error of the `(u, t)` grid is estimated separately (§4.3).

## 4. Validation of the completed moments

@@VALIDATION@@

## 5. Results: k = 2, 3 (`out_t1_lowk_e2e_l{1,2}.txt`, `out_t3_tables.txt`)

@@RESULTS@@

## 6. Non-analytic terms and Assumption T

@@NONANALYTIC@@

## 7. Corrections: `O(1/m)` and `O(1/m²)` at k = 2, 3 (`t2_corrections.py`, `out_t2_corrections.txt`)

@@CORRECTIONS@@

## 8. What remains

@@REMAINS@@

## 9. Code index (`ads4/code/lowk/`; `sh run_all.sh`)

| file | content | output |
|---|---|---|
| `lcconv.c` | `B_{n,j}` → lightcone coefficients `C_{i,q}` (C, OpenMP) | — |
| `lcblock.py` | `LCBlock`: lightcone reorganisation, Darboux tail fits, closed-form tails, derivative sums `S[k,l]` | — |
| `casimir_check.py` | Casimir residual from the jets | — |
| `c1_block_checks.py` | block-level checks (§2.3) | `out_c1_block_checks.txt` |
| `corner.py` | corner completion of the e2e mode grids, two settings, error budget, moments with errors | — |
| `v1_scalar.py` | scalar validation (§4.1) | `out_v1_scalar.txt` |
| `t1_lowk_e2e.py` | completed moments vs the leading formula, `l = 1, 2`, `k = 2, 3, 4, 6` | `out_t1_lowk_e2e_l{1,2}.txt` |
| `q1_grid.py` | quadrature-grid check (§4.3) | `out_q1_grid.txt` |
| `u1_uniformity.py` | uniformity of the heavy expansion of the shadow integrand (§6) | `out_u1_uniformity.txt` |
| `t2_corrections.py` | `m²(R−1)`, `c_1`, `c_2`, non-analytic fits (§7) | `out_t2_corrections.txt` |
| `t3_tables.py` | tables of this note | `out_t3_tables.txt` |

Caches (`cache/`, git-ignored): `lcmodes_*.npz` (completed mode grids), `data_t1_l*.pkl`, run logs.
