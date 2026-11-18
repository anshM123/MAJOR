# Flat functionals for `|ĝ₃|²` without forward-limit terms: an obstruction for smeared `B_{k≥4}`, and a certified `B₂ + B₃` functional

Status: written 17 Nov 2026 by a derivation/numerics agent; **revised 18 Nov 2026 after `ads4/REFEREE_NOFORWARD.md`** (changes in §7). Code: `ads4/code/robust/nf_*.py` (`bash nf_run_all.sh`;
outputs in `robust/cache/nf/`, git-ignored). Certified coefficients are in the tracked file
`ads4/code/robust/functionals_noforward.json`. `functionals_certified.json` (M01, M02) is untouched. The scripts reuse
`rcommon.py`, `lp.py`, `rcert.py` and `lowcheck.py` from `ROBUST_FUNCTIONAL.md`, read-only. `nf_cert.py` patches two `rcert` evaluation
helpers at run time to use centred interval evaluation (§4).

Motivation (`OMEGA_U.md` §3 I10, `SUPERCONVERGENT_UPLIFT.md`): in the AdS uplift, the forward terms `−a∂_{p²}B₄|₀ + aB₅|₀` become
`Λ¹⁰`-weighted `ν = 0` moments. Their leverage on light-side errors is `~10¹³`, and their double-trace and contact content is not
analysed. The question here is whether a flat `|ĝ₃|²` functional can do without them.

Marks:
* **[D]** derived;
* **[certified]** exact rationals or FLINT/Arb;
* **[NE]** floating point (LP, grids).

Notation: `ROBUST_FUNCTIONAL.md`, with `u = M²/m²`, `N = |B2| + |B3|`, and `q = p²`.

## 0. Results

1. **Obstruction [D; certified exact check].** Let `K ≥ 4` and `ψ_k ∈ L²[0,1]` for `4 ≤ k ≤ K`. Suppose
   `Σ_k ∫₀¹ψ_k(p)(−B^{(1)}_k|_low)dp` receives no contribution from any EFT contact term, for every symmetric polynomial
   `f(s,u) = f(u,s)`. Then `ψ_k = 0` for all `k`.
   * So, with no `p = 0` terms, contacts cannot be eliminated from smeared `B₄ … B_K`.
   * **Basis caveat.** In the stripped basis, the only contact-free forward-free functionals are `∫ψ₂B₂ + ∫ψ₃B₃`: `B₂^{(1)}`, `B₃^{(1)}`
     receive no contact term at all. In the unstripped `B′` basis (same span, `FLAT_UNSTRIPPED.md` item 2), which is the one relevant to
     the AdS uplift, `B′₂` and `B′₃` *do* see contacts. The same functionals are then combinations of `B′₂ … B′₇` with smeared weights,
     via `B₂ = B′₂ − 4qB′₃ + 2q²B′₄ − 3q³B′₅ + ½q⁴B′₆ − ½q⁵B′₇` (and similarly `B₃`). Contacts cancel only in that combination, so the AdS uplift of NF02 still involves `C₂ … C₇`,
     but with no forward (`ν = 0`) terms.
   * Eliminating only `g₄, g₅, g₆, g₆′` (the truncation printed in (2.7a)) is easy: with `ψ₄…ψ₇` of degree ≤ 16 there are 64 free
     directions. But it is not legitimate. A generic contact of degree `n` in `(s,u)` enters `−B_k|_low` as `c·p^{2(n−k+4)}` for every
     `k ≤ n+4` (special ones enter fewer, e.g. `su` only `B₄`, `B₅`). Since `p` runs up to `M`, higher contacts contribute at `O(1)`.
   * Forward terms evade the obstruction: at `p = 0` only finitely many contacts survive (`∂_{p²}B₄|₀`, `B₅|₀` see only `g₅`).
2. **The forward-free class `∫ψ₂B₂ + ∫ψ₃B₃` exists but is expensive [NE].**
   * The paper's footnote functional (arXiv:2201.06602 (3.7), `1598 log − 3211`) is of this form; its low side reproduces: `k = 8.006`,
     `A = 1596.7`, `B = −3208.0`. It **violates both necessary conditions below**: `ψ₃(0) = 0.0017 ≠ 0`, and only `(1−p)²` decay. It is
     therefore not positive on all heavy states (consistent with `EINSTEIN_OBSTRUCTIONS.md` A.4: negative at `m ≈ M`).
   * Two conditions are necessary [D]:
     * `ψ₃(0) = 0`. Otherwise, at fixed `u`, `𝔉 → σ_huψ₃(0)/b̂` as `b̂ → ∞`, with opposite signs in the two channels.
     * The weights must vanish like `(1−p)^5` at `p = 1`. The `+−` backward peak of `d̃^J_{44}` at `m = M` gives `~J^{6−2r}(−1)^J`
       against the forward `γ/(8J³)` (asymptotic argument, confirmed numerically).
   * LP optimum (`r = 5`, degree 8): `A = 876` (margin `δ = 0`), `1561` (`δ = 0.1`), `2258` (`δ = 0.2`), with `B ≈ −1.6A … −1.75A`.
   * Adding the linear proxy that makes the large-`b̂` certificate work at `b̂ ≥ 120` gives `1175 / 2055 / 3022` at `δ = 0 / 0.1 / 0.2`.
   * Compare `F_{g3}` (`37.84`), M01 (`37.55`, margin `9/100`) and M02 (`51.74`, margin `9/50`). **Dropping the forward terms costs a
     factor `≈ 33` in `A` at `δ = 0` (`≈ 860` against the forward-term LP optimum `26.1`; `≈ 23` against `F_{g3}`), and `≈ 45–60` at the
     same robust margin.**
   * At `Λ = 100` the AdS-style value `A log(2Λ/3) + B_AdS` is `≈ 6.1·10³` (NF02) against `131` (M02): a factor `≈ 47` in `|ĝ₃|²`,
     `≈ 6.8` in `√(t₄²+α_T²)`.
3. **Certified functional NF02 [certified in the regions stated] (§4).** The weights are
   `ψ₂ = (1−p)⁵p X(p)` and `ψ₃ = (1−p)⁵p Y(p)`, with `X`, `Y` of degree 8 and exact rationals (in `functionals_noforward.json`,
   sha256 `b8bd1237…20db4a`). There are no forward terms.
   * **Low side (exact):** `k = 88865643959/3233230000 = 27.48510`, `A = 3022.3377`, `B = −5327.2096`, `B_AdS = −6552.662`. `g₄…g₆′`
     and all higher contacts are absent. Both light-matter coefficients are negative on `0 < m_ℓ ≤ M`. These were recomputed from the
     amplitude by residues (`lowcheck.py`) and by an exact Descartes check.
   * **Heavy side:** `𝔉 ≥ (3/20)N` is certified on:
     * all states with `J < 400` (all `m ≥ M`);
     * all `m ≥ 2M` (`b̂ ≥ 120` by Lemma M; `b̂ ≤ 120` by the S cells);
     * `++` with `J ≥ 400`, `M ≤ m ≤ 2M`;
     * `+−` with `J ≥ 400` and `m ≥ (8/3)^{1/2}M = 1.633M`.
   * **Not certified:** `+−`, `J ≥ 400`, `M ≤ m < 1.633M` (the backward region near the gap). This region is **float-checked
     only**: here up to `J ≤ 3000`, and by the referee up to `J ≤ 20000`, with min `𝔉/N = 0.765` in this region [NE]. Why the
     certificate cannot be closed with the present lemmas, and the cost of a successor that could be, is in §4.1.
   * Dense check over all states: min `𝔉/N = 0.2000` in both channels (attained at `m = M`, `J = 4, 5`) [NE]. The flat envelope
     thresholds `n(C)` for `C = 1, 3, 10` are `1.47, 2.00, 3.43` (`++`) and `1.40, 2.00, 3.14` (`+−`); M02 has `1.91, 3.64, 5.98`.
   * All certificate logs print the sha256 above. All verdicts PASS: low side, P, M, S`++`, S`+−`, C`++`, C`+−` (`u ≤ 3/8`).

**Answer to the question.**
* Smeared `B₄…B₇` cannot help: no forward-free combination of them cancels all contacts (item 1, proved).
* In the stripped basis a forward-free, contact-free functional must be `∫ψ₂B₂ + ∫ψ₃B₃`.
* NF02 is such a functional, with margin `3/20`. It is **certified on all heavy states except `+−`, `J ≥ 400`, `M ≤ m < 1.633M`**, where
  it is float-checked only (min `𝔉/N = 0.765` up to `J = 20000`, referee). So full certification is **not** achieved.
* A successor that the present lemmas could certify needs `(1−p)^{12}` weights, at about `2–2.5×` the `A` (§4.1).
* Either way, the forward-free class is weaker than the forward-term functionals by a factor `≈ 45–60` in `A` at the same margin.

## 1. The obstruction

### 1.1 Contact data of `B^{(1)}_k` [D; checked exactly]
* For a polynomial (contact) amplitude, the integrand of (2.26) has poles only at `s = 0, q`. Put `h(s) := f(s, −q)` and
  `z := s(s−q)`, so that `(2s−q)ds = dz`. Then:
  * even `k`: `−B_k|_contact = [z^{k/2−2}] (h(s) + h(q−s))`;
  * odd `k`: `−B_k|_contact = [z^{(k−5)/2}] ((h(q−s) − h(s))/(2s−q))`.
* Both expansions exist because the brackets are symmetric under `s ↔ q−s`. The constant is 1, checked against the residues of
  `cert_low.py`'s method for all monomials of degree ≤ 6 and `k = 4…7` (`nf_contacts.py`).
* For `k = 2, 3` the integrand has no poles, so `B₂`, `B₃` are contact-free.
* By homogeneity, a degree-`n` contact contributes `c_k·q^{n−k+4}` (zero for `k > n+4`).

### 1.2 Proof (revised: the induction pairs `(B_{2ℓ+4}, B_{2ℓ+5})`)
*Data.* Write `f = P(t, w)`, `w = su`, with `P` an arbitrary polynomial; on `u = −q`, `h(s) = P(q−s, −qs)`. Since `z = −qs + s²`,
`[z^ℓ]` of a function regular at `s = 0` is a combination of its `s`-derivatives at 0 of order `≤ ℓ`, with top coefficient
`(−1/q)^ℓ/ℓ!`. The level-`ℓ` sum rules are `B_{2ℓ+4}` (from `F = h(s) + h(q−s)`) and `B_{2ℓ+5}` (from `G = F⁻/(2s−q)`,
`F⁻ = h(q−s) − h(s)`). Both use only data of order `≤ ℓ`. In terms of the top data `a_ℓ := h^{(ℓ)}(0)` and `c_ℓ := (−1)^ℓh^{(ℓ)}(q)`:
* `[z^ℓ]F = (−1/q)^ℓ/ℓ! · (a_ℓ + c_ℓ) + (orders < ℓ)`;
* `[z^ℓ]G = (−1/q)^ℓ/ℓ! · (c_ℓ − a_ℓ)/(−q) + (orders < ℓ)`, since `G^{(ℓ)}(0) = F⁻^{(ℓ)}(0)/(−q) + (lower)`.

*Isolating the two families.* Fix `L` larger than every order in play.
* `P = w^ℓ t^L g(t)` (`g` arbitrary): `h(s) = (−qs)^ℓ(q−s)^L g(q−s)`. All `s = 0` data of order `< ℓ` vanish and
  `a_ℓ = ℓ!(−q)^ℓq^L g(q)`. All `s = q` data of order `< L` vanish. So this `P` reaches only levels `≥ ℓ`, and at level `ℓ` it enters
  through `a_ℓ` alone. The functions `q^{L+ℓ}g(q)` give exponents `p^{2(L+ℓ)+2j}`, `j ≥ 0`.
* `P = t^ℓ w^L g(w)`: `h(s) = (q−s)^ℓ(−qs)^L g(−qs)`. The `s = 0` data of order `< L` vanish, and `c_ℓ = ℓ!(−q²)^L g(−q²)`. The
  functions are `q^{2L}g(−q²)`, with exponents `p^{4L+4j}`.

*Induction.* Let `ℓ_max` be the top level present. Suppose all weights of levels `> ℓ` already vanish.
* Pair the test polynomials of the previous step with the level-`ℓ` weights; levels `< ℓ` see nothing. The functional equals
  `(−1/q)^ℓ/ℓ!·∫[ψ_{2ℓ+4} + ψ_{2ℓ+5}/q]·a_ℓ dp` resp. `(−1/q)^ℓ/ℓ!·∫[ψ_{2ℓ+4} − ψ_{2ℓ+5}/q]·c_ℓ dp`.
* By Müntz–Szász in `L²[0,1]`, both sets `{2(L+ℓ)+2j}` and `{4L+4j}` have `Σ1/λ = ∞`, so the families are complete. Hence
  `ψ_{2ℓ+4} + ψ_{2ℓ+5}/q = 0` and `ψ_{2ℓ+4} − ψ_{2ℓ+5}/q = 0` a.e. The two data families enter with opposite relative signs, so the 2×2
  system `[[1, 1/q], [1, −1/q]]` has determinant `−2/q ≠ 0`, and both weights vanish.
* If the top level contains only `B_{2ℓ+4}` (`K` even), the `a`-family alone gives `ψ_K = 0`.

Descending to `ℓ = 0` gives `ψ_k = 0` for all `k ≥ 4`. ∎

*Checks.* `nf_contacts.py` verifies the structure for `ℓ ≤ 2` exactly: with `P = w^ℓt⁸` and `P = t^ℓw⁸`, the contributions vanish in
`B₄ … B_{2ℓ+3}`, and `B_{2ℓ+5}/B_{2ℓ+4} = +1/q` resp. `−1/q`. This matches the referee's independent residue check.
* **Exact check** (`nf_contacts.py`, FLINT rational rank; `ψ₄…ψ₇` arbitrary polynomials of degree ≤ 16, i.e. 68 unknowns, against
  all contact monomials `t^a(su)^b` of degree ≤ `N`):

  | contacts up to degree `N` | 2 (`g₄…g₆′`) | 4 | 8 | 12 | 16 | 20 | 24 | 28 | 32 |
  |---|---|---|---|---|---|---|---|---|---|
  | null-space dimension | 64 | 59 | 47 | 35 | 23 | 11 | 3 | 0 | 0 |

  `VERDICT nf_contacts: PASS`. The rank grows by 3 per degree for 4 weights: the `s = q` data are functions of `q²` only. This is why
  the completeness (Müntz) step is needed and why low-degree truncations look deceptively free.

### 1.3 Consequences for the LP class [D]
* No forward term means no `2a𝒥u³` growth, so:
  * `ψ₃(0) = 0` is required (item 2 of §0);
  * the large-`b̂` action is `−[2ψ₂″(0) + σ_huψ₃″(0)]/(2b̂³) + O(b̂⁻⁴)`, so `−2ψ₂″(0) > |ψ₃″(0)|` is needed (with margin);
  * at `m → M` the `+−` backward peak (`d̃^J_{44}(−1) = binom(J+4,8)`) must be suppressed by `(1−p)^r` with `r ≥ 5`.
* `k = −2μ₂₆ + 2μ₃₄` loses the `4a` that supplies 20 of `F_{g3}`'s 27.5. The cubic coupling must then be measured by `μ₃₄`, i.e. through
  `B₃`, whose heavy density has the channel-odd sign. This is the origin of the cost.
* Numerical note: the flat parts at `u ≈ 1`, `J ≳ 100` must be evaluated with the weights in **factored** form `(1−p)^r p X(p)`. The
  expanded monomial form loses all relative accuracy near `p = 1`, where `P^{(0,8)}_{J−4} ~ J⁸`: our first LP runs were wrong for
  this reason; caught by an mpmath comparison.
  This does not affect any certificate of `ROBUST_FUNCTIONAL.md`, which are exact. Its float scans there were dominated by the
  forward term at `u ≈ 1`.

## 2. LP results [NE] (`nf_lp.py`)
Weights `(1−p)^r p × (shifted Legendre basis, degree d)`; constraints:
* (R1) `𝔉 ≥ δN` on ≈ 2.1·10⁵ states (`u ∈ [10⁻⁴, 1]`, `J ≤ 1200` for `u ≥ ¼`), by cutting planes;
* (R2) scaling limit `2f̂ ≥ δρ`;
* (R3) `1/b̂³` coefficient with margin, for `u ≤ 1`;
* matter;
* normalisation `k = 1`.

Objective: minimise `A`.

| `r, d` | `δ = 0` | 0.05 | 0.10 | 0.15 | 0.20 | 0.30 |
|---|---|---|---|---|---|---|
| 3, 6 | 1012.6 / −1630 | — | 2199 / −3434 | 2657 / −4201 | 3124 / −4981 | — |
| 5, 6 | 947.6 / −1472 | 1378 / −2378 | 1737 / −2948 | 2088 / −3662 | 2478 / −4462 | 3356 / −6578 |
| 5, 8 | **875.6 / −1357** | 1225 / −1993 | 1561 / −2607 | 1901 / −3239 | 2258 / −3924 | 3084 / −5649 |
| 6, 8 | 910.0 / −1419 | 1200 / −1938 | 1532 / −2542 | 1844 / −3107 | 2195 / −3778 | 3035 / −5576 |
| 5, 12 | 864.3 / −1341 | | 1466 / −2412 | | 2154 / −3692 | |
| 11, 10 | 2589 / −4553 | | 3306 / −5674 | | 4310 / −7243 | |
| 5, 8 + Lemma-M proxy (`B_A = 120`) | 1175 / −2099 | | 2055 / −3587 | 2517 / −4297 | **3022 / −5327** | |

(entries `A / B`; "—": cutting planes did not converge). The `δ = 0` optimum converges slowly with degree (`≈ 860`), i.e. `≈ 33×` the
forward-term class LP optimum (`δ = 0`: 26.1, `ROBUST_FUNCTIONAL.md` §2.1); `≈ 23×` against `F_{g3}`'s 37.84.

## 3. NF02
* `nf_mkfun.py`: the LP solution (`r = 5, d = 8, δ = 0.2`, proxy) is scaled to `k ≈ 27.485` and rounded to `10⁻³`. `Y₀` is raised until
  the scalar-matter coefficient is strictly negative. The result is copied to `functionals_noforward.json`.
* The cofactors `X`, `Y` have monomial coefficients up to `~10⁹`, i.e. strong cancellation. All interval evaluations therefore use
  centred forms (§4).
* Dense floating-point check (`nf_check.py`; 144 `u` values, `J ≤ 3000` for `u ≥ ¼`, factored weights):
  * min `𝔉/N = 0.20000`: `++` at `m = M, J = 4`; `+−` at `m = M, J = 5`;
  * flat envelope thresholds `n(C) = sup{m : 𝔉 < CuN}` for `C = 1, 3, 10`: `1.470, 2.000, 3.433` (`++`) and `1.397, 2.000, 3.138` (`+−`).
* The certified margin `3/20` leaves 25 % slack.

## 4. Certification (`nf_cert.py`, `nf_lowcheck.py`)
* Same reduction as `ROBUST_FUNCTIONAL.md` §4: four sign patterns, with `W₀ = 2(1−δ′s₂)ψ₂`,
  `W₁ = −(1−δ′s₂)p²ψ₂ + (1−δ′s₃)σψ₃`, and no forward part.
* Interval evaluation of the weights and of their Taylor jets on a `p`-cell uses the exact rational Taylor shift to the cell midpoint
  and then Horner in the small offset. The monomial interval Horner over-estimates by `~10⁷` for these coefficients.

| region | method | NF02, `δ′ = 3/20` |
|---|---|---|
| P: `J < 400`, all `m ≥ M`, both channels | exact Lemma-P polynomials, Descartes (2384 polynomials) | PASS (2384 polynomials) |
| M: `m ≥ 2M`, `b̂ ≥ 120` | Lemma M, no forward term: `b̂⁴𝔉 ≥ γ(u)b̂/(1+δ_s)³ − C_R(u)/(1−δ_s)⁴`, `u` split at `1/100, 1/25, 1/10, 1/4` | PASS (tightest: `+−`, `u ∈ [1/10,1/4]`, `b̂⁴𝔉 ≥ 4.3·10³` vs `C_R ≤ 2.27·10⁶`) |
| S: `J ≥ 400`, `m ≥ 2M`, `b̂ ≤ 120` | `rcert.region_S` (Taylor/Cauchy cells, 604 per pattern), forward part 0 | PASS (`++`, `+−`; 2416 cells each) |
| C, `++`: `J ≥ 400`, `M ≤ m ≤ 2M` | Lemma M with the Legendre kernel (`m ≥ M`), `b̂ ≥ (2·400+1)√u ≥ 400.5`, `u` in six slices | PASS (`b̂⁴𝔉 ≥ 9·10⁶`) |
| C, `+−`: `J ≥ 400`, `1.633M ≤ m ≤ 2M` (`u ≤ 3/8`) | Lemma M + Lemma C (`(1−up²)^{−8}`), `u` in two slices | PASS |
| C, `+−`: `J ≥ 400`, `M ≤ m < 1.633M` | — | **not certified** (§4.1); float only: min `𝔉/N = 0.765` (`J ≤ 20000`, referee) |
| low side | exact (residues + Descartes) | PASS |

### 4.1 The uncovered `+−` region: can it be closed?
* **Present route.** Lemma C writes `d̃^J_{44} = (1−up²)^{−8}Σc_jP_{J+j}`, and Lemma M is applied to `w₋ = Ψ₋(1−up²)^{−8}`.
  * Put `δ = 1−p` and `ε = 1−u`. Near `p = 1`, `1−up² ≈ 2δ + ε`, and `Ψ₋ ~ δ⁵`, so `w₋ ~ δ⁵/(2δ+ε)⁸`. This reaches `~ε^{−3}` at
    `δ ~ ε`, and `w₋` is not even bounded uniformly as `u → 1`.
  * Hence `C_R(u) → ∞` for any `J`: it is `6·10⁶` at `u ≤ 3/8`, `2·10⁷` on `[3/8, 1/2]`, `1.2·10⁸` on `[1/2, 5/8]`, against `γb̂ ≈ 10⁷`.
  * Raising `JP` does not help, since the failure is at fixed `u`, not at large `J`. Finer Cauchy cells do not help either: region S
    is not the obstacle.
* **What a certifiable successor needs.**
  * Uniform `C⁴` bounds on `w₋` up to `u = 1` need `Ψ₋ ∝ (1−p)^r` with `r ≥ 12`. Then `w₋ = δ^{r−8}(δ/(2δ+ε))⁸·(…)` has derivatives of
    order ≤ 4 bounded uniformly in `ε`; with `r = 11` the fourth derivative is `~1/δ` and `‖w⁗‖₁ ~ log(1/ε)`.
  * The interval jets would also have to be evaluated in this factored form near the corner `(p, u) → (1, 1)`. That code does not
    exist yet.
  * Cost (LP with the Lemma-M proxy, `δ = 0.2`, [NE]):

    | weights | `A` | relative to NF02 (3022) |
    |---|---|---|
    | `(1−p)^{11}`, `d = 10` | 5468 | 1.8× |
    | `(1−p)^{12}`, `d = 10` | 7441 | 2.5× |
    | `(1−p)^{12}`, `d = 12` | 6135 | 2.0× |

    So a fully certifiable forward-free functional costs `≈ 2–2.5×` NF02's `A`, i.e. `≈ 120×` M02, before any certificate is attempted.
* **A dedicated backward-endpoint lemma** would avoid the CG factor and could certify NF02 itself. It would be a Hilb/Mehler-type
  expansion of `P^{(0,8)}_n` near `x = −1` with explicit remainder, giving the backward contribution `O(J^{6−2r})` against the forward
  `γ/(8J³)` with explicit constants. It is **not done**.
* **Status: NF02 is not fully certified, and no successor has been certified.**

## 5. Remarks
* The forward-free functional removes the `Λ¹⁰` forward moments (`C⁺₄`, `C⁻₅` at `ν = 0`) from the AdS uplift. Its AdS version uses only
  `C⁺₂` and `C⁻₃`, the two moments whose light side `OMEGA_U.md` treats as best controlled. The price is the factor `≈ 47` in the
  would-be `|ĝ₃|²` bound at `Λ = 100` (`≈ 6.8` in `√(t₄²+α_T²)`). The AdS heavy-side analysis (Regge/hybrid scans) of NF02 is not done here.
* Other families (fixed-`t` `B^{(2)}`, the improved sum rules of the paper §4) are not covered by item 1. The paper's improved
  `B₂^{imp}`, `B₄^{imp}` are themselves built from forward limits.
* Novelty: the obstruction is elementary (it is the reason the literature uses forward limits or improved sum rules). The exact
  statement and the certified forward-free functional with a margin are new here.

## 6. Code and reproduction (`ads4/code/robust/`)
| file | content |
|---|---|
| `nf_contacts.py` | contact data of `B_k`, closed form vs residues, proof-structure check (§1.2), exact null-space ranks |
| `nf_lp.py` | LP (factored weights, Legendre basis, cutting planes, optional Lemma-M proxy) |
| `nf_mkfun.py` | rationalisation |
| `nf_check.py` | dense floating-point check |
| `nf_cert.py`, `nf_lowcheck.py` | certificates (fingerprint printed); reads `functionals_noforward.json` first |
| `functionals_noforward.json` | **tracked**: exact NF02 |
| `nf_run_all.sh` | everything |

## 7. Revisions after `ads4/REFEREE_NOFORWARD.md`
* **M1.** §1.2 rewritten. Both parities use data of order `≤ ℓ`; the pair is `(B_{2ℓ+4}, B_{2ℓ+5})`; the two data families enter with
  relative signs `±1/q`, the system is nonsingular, and Müntz applies. The structure is now checked exactly in `nf_contacts.py`.
* **M2.** The headline and answer now state that `+−`, `J ≥ 400`, `M ≤ m < 1.633M` is float-checked only (referee: min 0.765 up to
  `J = 20000`). The attempt to close it is reported in §4.1: not closed. A `(1−p)^{12}` successor would cost `≈ 2–2.5×` in `A`.
* **Minor fixes.**
  * "≈ 23×" corrected to `≈ 33×` against the LP optimum.
  * "Every `k ≤ n+4`" restricted to generic contacts.
  * The footnote functional is now said to violate both necessary conditions.
  * The `nf_cert.py` docstring is updated.
  * Basis caveat added: in the `B′` basis NF02 is a contact-cancelling combination of `B′₂ … B′₇`, with no forward terms.
