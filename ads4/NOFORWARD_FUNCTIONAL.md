# Flat functionals for `|ĝ₃|²` without forward-limit terms: an obstruction for smeared `B_{k≥4}`, and a certified `B₂ + B₃` functional

Status: written 17 Nov 2026 by a derivation/numerics agent. **Unrefereed.** Code: `ads4/code/robust/nf_*.py` (`bash nf_run_all.sh`;
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
   * So, with no `p = 0` terms, contacts cannot be eliminated from smeared `B₄ … B_K`. This holds in the stripped basis and
     equivalently in the unstripped `B′` basis, which spans the same space (`FLAT_UNSTRIPPED.md` item 2).
   * The only contact-free forward-free functionals in these families are `∫ψ₂B₂ + ∫ψ₃B₃`. `B₂^{(1)}` and `B₃^{(1)}` receive no
     contact term at all.
   * Eliminating only `g₄, g₅, g₆, g₆′` (the truncation printed in (2.7a)) is easy: with `ψ₄…ψ₇` of degree ≤ 16 there are 64 free
     directions. But it is not legitimate. A contact of degree `n` in `(s,u)` enters `−B_k|_low` as `c·p^{2(n−k+4)}` for every
     `k ≤ n+4`, and `p` runs up to `M`, so every higher contact contributes at `O(1)`.
   * Forward terms evade the obstruction: at `p = 0` only finitely many contacts survive (`∂_{p²}B₄|₀`, `B₅|₀` see only `g₅`).
2. **The forward-free class `∫ψ₂B₂ + ∫ψ₃B₃` exists but is expensive [NE].**
   * It contains the paper's footnote functional (arXiv:2201.06602 §3, `1598 log − 3211`), which we reproduce: `k = 8.006`,
     `A = 1596.7`, `B = −3208.0`.
   * Two conditions are necessary [D]:
     * `ψ₃(0) = 0`. Otherwise, at fixed `u`, `𝔉 → σ_huψ₃(0)/b̂` as `b̂ → ∞`, with opposite signs in the two channels.
     * The weights must vanish like `(1−p)^5` at `p = 1`. The `+−` backward peak of `d̃^J_{44}` at `m = M` gives `~J^{6−2r}(−1)^J`
       against the forward `γ/(8J³)` (asymptotic argument, confirmed numerically).
   * LP optimum (`r = 5`, degree 8): `A = 876` (margin `δ = 0`), `1561` (`δ = 0.1`), `2258` (`δ = 0.2`), with `B ≈ −1.6A … −1.75A`.
   * Adding the linear proxy that makes the large-`b̂` certificate work at `b̂ ≥ 120` gives `1175 / 2055 / 3022` at `δ = 0 / 0.1 / 0.2`.
   * Compare `F_{g3}` (`37.84`), M01 (`37.55`, margin `9/100`) and M02 (`51.74`, margin `9/50`). **Dropping the forward terms costs a
     factor `≈ 23` in `A` at `δ = 0`, and `≈ 45–60` at the same robust margin.**
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
   * **Not certified:** `+−`, `J ≥ 400`, `M ≤ m < 1.633M` (the backward region near the gap). There the Clebsch–Gordan route of
     Lemma C introduces `(1−up²)^{−8}`, and with `(1−p)⁵` weights its remainder blows up as `u → 1`. Dense floating point, including
     that region with `J ≤ 3000`: no violation [NE].
   * Dense check over all states: min `𝔉/N = 0.2000` in both channels (attained at `m = M`, `J = 4, 5`) [NE]. The flat envelope
     thresholds `n(C)` for `C = 1, 3, 10` are `1.47, 2.00, 3.43` (`++`) and `1.40, 2.00, 3.14` (`+−`); M02 has `1.91, 3.64, 5.98`.
   * All certificate logs print the sha256 above. All verdicts PASS: low side, P, M, S`++`, S`+−`, C`++`, C`+−` (`u ≤ 3/8`).

**Answer to the question.** A forward-free, contact-free flat functional for `|ĝ₃|²` with a robust margin exists. It necessarily uses
only smeared `B₂`, `B₃`, and it is weaker than the forward-term functionals by a factor `≈ 45–60` in `A` at the same margin. Smeared
`B₄…B₇` cannot help (item 1).

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

### 1.2 Proof
* Write `f = P(t, w)`, `w = su`, with `P` an arbitrary polynomial. On `u = −q`: `h(s) = P(q−s, −qs)`. The `z`-coefficient
  `[z^ℓ]` is a fixed combination of `h^{(i)}(0)` and `h^{(i)}(q)` for `i ≤ ℓ` (even `k`) or `i ≤ ℓ+1` (odd `k`). The top derivative
  enters with a nonzero coefficient (`(−1/q)^ℓ/ℓ!`).
* `h^{(i)}(0)` involves `∂_w^iP(q, 0)·(−q)^i` plus lower data, and `∂_w^iP(t, 0)` is an arbitrary new polynomial in `t`.
* `h^{(i)}(q)` involves `∂_t^iP(0, −q²)·(−1)^i` plus lower data. The new piece is an arbitrary polynomial in `q²`.
* Taking `P = w^i t^L g(t)` or `P = t^i w^L g(w)` with `L` large isolates each kind of datum. Then the data are an arbitrary
  polynomial in `p²` times `p^{2L}`, respectively an arbitrary polynomial in `p⁴` times `p^{2i+4L}`.
* Induct downward from the top `k`:
  * The highest odd `k = 2L+5` alone carries the order-`(L+1)` `s = 0` datum. So `ψ_{2L+5}·(coefficient)·p^{2(L+1)}` is orthogonal
    to all `p^{2j}`, and `ψ_{2L+5} = 0`.
  * The pair `(2L+4, 2L+3)` then carries the order-`L` data at `s = 0` and `s = q` with independent sign combinations, which forces
    both weights to vanish.
* The density statements are the Müntz–Szász theorem in `L²[0,1]`: the exponent sets `{2j + c}` and `{4j + c}` have `Σ 1/λ = ∞`. ∎
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

(entries `A / B`; "—": cutting planes did not converge). The `δ = 0` optimum converges slowly with degree (`≈ 860`), i.e. `≈ 23×` the
forward-term class (`δ = 0`: 26.1, `ROBUST_FUNCTIONAL.md` §2.1).

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
| C, `+−`: `J ≥ 400`, `M ≤ m < 1.633M` | — | **not certified**: `C_R ≥ 2·10⁷` at `u ∈ [3/8, 1/2]` against `γb̂ ≈ 1.2·10⁷`, growing without bound as `u → 1` (would need weights `∝ (1−p)^{11}`, which costs a further factor ≈ 1.5–2 in `A`, table §2) |
| low side | exact (residues + Descartes) | PASS |

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
| `nf_contacts.py` | contact data of `B_k`, closed form vs residues, exact null-space ranks |
| `nf_lp.py` | LP (factored weights, Legendre basis, cutting planes, optional Lemma-M proxy) |
| `nf_mkfun.py` | rationalisation |
| `nf_check.py` | dense floating-point check |
| `nf_cert.py`, `nf_lowcheck.py` | certificates (fingerprint printed); reads `functionals_noforward.json` first |
| `functionals_noforward.json` | **tracked**: exact NF02 |
| `nf_run_all.sh` | everything |
