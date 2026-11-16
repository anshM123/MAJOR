# Flat dispersive sum rules for the unstripped graviton helicity amplitude, and the `|ĝ₃|²` bound they give

Status: written 16 Nov 2026 by an automated analyst; not refereed. Code: `ads4/code/flatunstripped/` (`bash run_all.sh`, ≈ 1 min;
`--scan` adds the floating-point scan; outputs in `cache/`, git-ignored). It imports nothing from other code directories; it reads the
certified coefficients of M01/M02 from `ads4/code/robust/functionals_certified.json`. Context: `ADS4_FUNCTIONAL.md` §7,
`SKELETON_TEST.md` §8, `FLAT_GRAVITON_BOUND.md` (Theorem 1, 3), `ROBUST_FUNCTIONAL.md` (M01, M02, Theorem R).

Marks: **[D]** derived, with exact symbolic verification in `sumrules.py`; **[certified]** Arb ball arithmetic or exact rationals;
**[NE]** floating point; **[I]** inferred/assumed. Units `M = 1`, `q = p²`, heavy `w = m²`, `Z := 1/(w(w−q))`, `P := d^J_{00}(x) = P_J(x)`,
`d := d^J_{44}(x)` (full Wigner), `x = 1 − 2q/w`.

## 0. Results

1. **Unstripped sum rules [D].** With `𝓜(s,u) := t⁴f(s,u)` — the helicity amplitude `M(1⁺2⁻3⁻4⁺)` in the CM phase convention of
   arXiv:2201.06602 App. A, where `|⟨23⟩[14]|⁴ = t⁴` — define at fixed `u = −q`, for `k ≥ 2`:
   `B′_k(q) = ∮ds/(4πi) (2s−q)/[s(s−q)]^{(k+2)/2} [𝓜(s,u) + 𝓜(q−s,u)]` (`k` even),
   `B′_k(q) = ∮ds/(4πi) [s(s−q)]^{−(k+1)/2} [−𝓜(s,u) + 𝓜(q−s,u)]` (`k` odd).
   * **Regge:** the kernel is `O(s^{−k−1})` and `𝓜 = O(s^{J₀})` with `J₀ < 2` (the paper's (2.23)/(2.25), smeared), so `B′_k`
     converges for every `k ≥ 2`. Tree-level Einstein exchange (`J₀ = 2`) is exactly marginal at `k = 2` and appears as the
     `16πG/q` pole on the low side.
   * **Heavy side:** on both cuts the prefactors cancel (`Im 𝓜 = 16π d` on the `s`-cut, `16π P` on the `t`-cut). Per state,
     dropping `16(2J+1)/π`, the raw densities are `(2w−q)Z^{(k+2)/2}(|c⁺⁺|²P + |c⁺⁻|²d)` for even `k` and
     `Z^{(k+1)/2}(|c⁺⁺|²P − |c⁺⁻|²d)` for odd `k`. These are the full `d₀₀`, `d₄₄` with a common prefactor; the `B₃`-type sign flip
     between channels is the odd-kernel sign. The leading term `2/m^{2k}` matches the AdS `C^±_k` with `k`-parity locked to the
     kernel parity.
   * **Low side:** from (2.7a) and (B.1), e.g. `−B′₂|_low = 16πG/q + 2πG|ĝ₃|²q³ + 4g₄q² + 5g₅q³ + (6g₆ − g₆′)q⁴ + …` and
     `−B′₃|_low = 2πG|ĝ₃|²q² + 2g₄q + 3g₅q² + …` (§3). **Every unstripped sum rule sees contact couplings**; none is superconvergent.
2. **The unstripped and stripped families span the same space [D, exact].** Exact identities of kernels, hence of heavy densities
   and low sides (all three checked independently):
   * `B^{(1)}_k = B′_k + 2q²B′_{k+2} + ½q⁴B′_{k+4} − 4qB′_{k+1} − 3q³B′_{k+3} − ½q⁵B′_{k+5}` (`k` even);
   * `B^{(1)}_k = B′_k + 2q²B′_{k+2} + ½q⁴B′_{k+4} − qB′_{k+1} − ½q³B′_{k+3}` (`k` odd);
   * conversely `B′₂ = B₂ + 4qB₃ + 2q²B₄ + 3q³B₅ + ½q⁴B₆ + ½q⁵B₇`, `B′₃ = B₃ + qB₄ + 2q²B₅ + ½q³B₆ + ½q⁴B₇` (`B_i = B^{(1)}_i`).

   So stripping is a change of basis with polynomial-in-`q` coefficients. It costs five extra units of Regge spin (`k → k+5`), and in
   the unstripped basis the contacts cancel between `B′_k … B′_{k+5}`.
3. **Bound from unstripped sum rules [D + certified].** `F_{g3}` (completed), M01 and M02, rewritten in the unstripped basis
   (weights on `B′₂ … B′₇` and forward terms `−a∂_qB′₄|₀ + 5aB′₅|₀`), are functionals of unstripped sum rules.
   * **Low sides** recomputed *from the unstripped low sides*: `g₄, g₅, g₆, g₆′` cancel exactly; `k(0)`, `A` and `B` coincide with
     the stripped values, and the matter terms equal the stripped ones pointwise.
   * **Heavy actions** are identical as functions of `(m, J, h)` (exact density identity). So `FLAT_GRAVITON_BOUND.md` Theorem 1 and
     `ROBUST_FUNCTIONAL.md` Theorem R apply verbatim: positivity on all heavy states, with certified uniform margins 9/100 (M01) and
     9/50 (M02).

   | functional (unstripped basis, `B′₂…B′₇` + forward `B′₄, B′₅`) | `A′` | `B′` (`ε → 0`) | heavy positivity |
   |---|---|---|---|
   | `F_{g3}` | 37.8387 | −45.4452 | certified, margin ≈ 0 (near-marginal) |
   | M01 | 37.5457 | −46.9286 | certified, `𝔉 ≥ (9/100)N` |
   | M02 | 51.7366 | −65.3207 | certified, `𝔉 ≥ (9/50)N` |

   An LP over unstripped sum rules with polynomial weights is the stripped LP in another basis (item 2), so with matching weight
   degrees its optimum equals the stripped optimum; no new LP is needed. **The flat bound from unstripped sum rules is the same as
   the stripped one:** `|ĝ₃|²M⁸ ≤ 37.84 log(M/m_IR) − 45.45` (`F_{g3}`) and `≤ 37.55 log − 46.93` (M01), with the hypotheses and IR caveats
   of `FLAT_GRAVITON_BOUND.md` Theorem 3.
4. **`F′` (`F_{g3}` with `d̃₄₄ → d₄₄`) is not a combination of fixed-`u` sum rules [D].** Its light side is therefore not defined by any
   sum rule. The heavy density of any (smeared + forward) combination has the form `E(w,q)(P+d) + O(w,q)(P−d)`, with `E`, `O` built
   from `Z^n` (`Z` is symmetric under `w ↔ q−w`). `F′` would need `E, O ∝ 1/w⁴`, but `1/w⁴` is not a power series in `Z`: its
   expansion in `√Z` has an odd term `−2q Z^{5/2}`. Equality of heavy actions for all `J` forces pointwise equality of the
   `p`-densities, by completeness of `{P_J}` on `[1−2/m², 1]`. The floating-point positivity of `F′` (`ADS4_FUNCTIONAL.md` §7.3) is
   thus not a statement about any sum rule.
5. **The naive transcription of `Ω_{g3}` is not a valid functional [D + certified].** If the AdS moments tend exactly to the unstripped
   sum rules, `ε_kΛ^{2k}C^{(−)^k}_{k,Λp} → B′_k` (the natural spinning analogue of the CHMRS bulk-point formula, whose scalar density
   `(2m²−ν²)P/(m^k(m²−ν²)^{(k+2)/2})` has exactly the `B′_k` form [I]), then `Ω_{g3}` as written in `ADS4_FUNCTIONAL.md` has flat limit
   `F_U = ∫[ψ₂B′₂ + ψ₃B′₃] − 5∂_qB′₄|₀ + 5B′₅|₀`. This `F_U` fails on both counts:
   * **Low side:** it keeps the contacts, `+ (60526/3465)g₄ − (2746189/360360)g₅ + (729983/90090)g₆ − (358681/180180)g₆′ + …`
     (`ε → 0`). Since `g₄ ≥ 0` enters with a positive coefficient, `F_U` bounds nothing.
   * **Heavy side:** it is **not positive**: `m⁴F_U(++, J = 2, m = 213/200) = −3.33384834423(±3·10⁻¹²)` and `(m = 6/5) = −1.12774359442`
     [certified, Arb integration]. In the `+−` channel the scan is positive (min `3.6·10⁻⁵` at `m ≈ 21.7, J = 244`) [NE].
   * **Correct uplift [D given the dictionary, I as a statement about AdS].** The uplift of `F_{g3}` must carry the identity weights
     `φ₂ = ψ₂`, `φ₃ = ψ₃ − 4x²ψ₂`, `φ₄ = 2x⁴ψ₂ − x²ψ₃`, `φ₅ = −3x⁶ψ₂ + 2x⁴ψ₃`, `φ₆ = ½x⁸ψ₂ − ½x⁶ψ₃`, `φ₇ = −½x¹⁰ψ₂ + ½x⁸ψ₃`:

     `Ω^U_{g3} = ∫₀¹dx Σ_{j=2}^{7} φ_j(x) ε_j Λ^{2j} C^{(−)^j}_{j,Λx} − 5Λ¹⁰∂_{ν²}C⁺_{4,ν}|₀ + (25/2)Λ¹⁰C⁻_{5,0}`

     (`ε_j = 1` even, `½` odd). Its flat limit is `F_{g3}` exactly, including the light side, so the contacts cancel. `ADS4_FUNCTIONAL.md`
     §7.2's "no combination of the `C^±` has flat limit `F_{g3}`" holds only for `k`-diagonal combinations at leading order in `1/m²`.
     With the exact densities and `ν`-dependent coefficients across `k = 2…7` it is false.
6. **What is assumed.** The flat statements 1–4 use only the hypotheses of `FLAT_GRAVITON_BOUND.md` §4.1. These include smeared
   `B′_k` for `k ≤ 7` (and `≤ 9` for the forward pieces of `B^{(1)}_4`, `B^{(1)}_5` if those are smeared), and unsmeared forward limits
   `∂_qB′₄|₀`, `B′₅|₀`. Item 5's AdS part is [I]: it needs both "reading B" for gravitons and the exact spinning bulk-point limit
   `C_k → ε_k^{−1}B′_k`.

## 1. Definitions and heavy densities [D]

On real kinematics `|⟨ij⟩| = |[ij]| = √|s_ij|`, and in the CM spinors of App. A of arXiv:2201.06602 `⟨23⟩⁴[14]⁴ = t⁴` (phase 1).
`𝓜 := t⁴f` is analytic with the cuts of `f`. At fixed `u = −q` its `s`-cut carries `+−` states and its `t`-cut (`s = q − m²`) carries
`++` states:
* `Im 𝓜(m², −q) = (q−m²)⁴·(16π/m⁸)Σ(2J+1)|c⁺⁻|²d̃₄₄ = 16πΣ(2J+1)|c⁺⁻|²d₄₄`, because `(1−q/m²)⁴d̃₄₄ = d₄₄`;
* `Im` of `𝓜(q−s,u) = s⁴f(t,u)` at `s = m²` is `m⁸(16π/m⁸)Σ|c⁺⁺|²P_J = 16πΣ|c⁺⁺|²P_J`.

The bracket `𝓜(s) ± 𝓜(q−s) = t⁴f(s,u) ± s⁴f(t,u)` is the unstripped analogue of the paper's `f(s,u) ± f(t,u)` (2.26). Reflecting
`s → q − s` turns each sum rule into `∮Keff(s)f(s,u)ds/(4πi)` with
* stripped: `Keff = ±2K^s_k`;
* unstripped: `Keff = ±2K′_k t⁴` (`+` even, `−` odd).

These effective kernels are what the identities of §2 compare. The Regge degree is `deg K′_k = −k−1` (`sumrules.py` U2).

## 2. Identities [D, `sumrules.py` U2–U3]

Write `a = 1/w`, `b = 1/(w−q)`. Then `ab = Z`, `a − b = −qZ`, `a + b = (2w−q)Z`, and `a⁴ ± b⁴` are polynomials in `Z` with factors
`(2w−q)` (using `(2w−q)² = 4/Z + q²`):
* `a⁴ + b⁴ = q⁴Z⁴ + 4q²Z³ + 2Z²`;
* `a⁴ − b⁴ = −q(2w−q)(q²Z⁴ + 2Z³)`.

Since the stripped densities are `(2w−q)Z^{(k−2)/2}(a⁴P + b⁴d)` (even) and `Z^{(k−3)/2}(a⁴P − b⁴d)` (odd), splitting `a⁴P ± b⁴d` into
`(P+d)`- and `(P−d)`-parts gives the identities of §0 item 2. `sumrules.py` checks them three ways:
* as identities of effective kernels (rational functions of `s`);
* as identities of heavy densities;
* as identities of low sides, computed independently by residues, for the EFT (2.7a) and for matter (B.1).

The inverse identities were obtained by solving for the coefficients.

## 3. Low sides of `B′_k` [D, from (2.7a), (B.1)]

Through the contacts kept in (2.7a); higher contacts add higher powers of `q`:

| `k` | `−B′_k|_low` (EFT) |
|---|---|
| 2 | `16πG/q + 2πG|ĝ₃|²q³ + 4g₄q² + 5g₅q³ + (6g₆ − g₆′)q⁴` |
| 3 | `2πG|ĝ₃|²q² + 2g₄q + 3g₅q² + (4g₆ − g₆′)q³` |
| 4 | `2g₄ + (4πG|ĝ₃|² + 5g₅)q + (9g₆ − 3g₆′)q²` |
| 5 | `g₅ + (3g₆ − g₆′)q` |
| 6 | `2g₆` |
| 7 | `0` |

Matter: `−B′₂ = (−2a₀m⁶ − 3a₀m⁴q − 2a₂m⁶ + 9a₂m⁴q + 6a₂m²q² + 6a₂q³)/m⁴` and `−B′₃ = (−(a₀+a₂)m⁴ + 6a₂m²q + 6a₂q²)/m⁴` (`m = m_ℓ`).
Compare the stripped values: `−B^{(1)}_2 = 16πG/q + 2πG|ĝ₃|²q³` and `−B^{(1)}_3 = −2πG|ĝ₃|²q²` have no contacts.
**Superconvergence is a property of the stripped basis**: the factor `t⁴` gives every contact term Regge spin `≥ 4`, and the
unstripped kernels `O(s^{−k−1})` with `k ≤ 3` do not suppress it.

## 4. The three certified functionals in the unstripped basis [D + certified]

For `F = ∫_ε¹[ψ₂B^{(1)}_2 + ψ₃B^{(1)}_3] − a∂_qB^{(1)}_4|₀ + aB^{(1)}_5|₀`:
* the smeared part becomes `∫Σ_jφ_jB′_j`, with `φ_j` as in §0 item 5 (degrees: `B′₂` 6–8 … `B′₇` 16–18);
* the forward part becomes `−a∂_qB′₄|₀ + 5aB′₅|₀`, using `∂_qB^{(1)}_4|₀ = ∂_qB′₄|₀ − 4B′₅|₀` and `B^{(1)}_5|₀ = B′₅|₀`.

`sumrules.py` U4 recomputes the low side from the `B′_j` low sides and finds:

| | contacts | `k(0)` | `A` | `B` |
|---|---|---|---|---|
| `F_{g3}` | 0 | `330151/12012` | 37.838686 | −45.445232 |
| M01 | 0 | `4951986353/180180000` | 37.545657 | −46.928647 |
| M02 | 0 | `10884217/396000` | 51.736565 | −65.320678 |

These agree with `FLAT_GRAVITON_BOUND.md` and `ROBUST_FUNCTIONAL.md`. The heavy actions are identical (§2), so the existing certificates
cover them.

## 5. `F′` and `F_U`

* **`F′`** [D, U6]: `1/w⁴ = Z² − 2qZ^{5/2} + 2q²Z³ + …` in `r = √Z`, while every sum-rule combination has `E, O ∈ span_q{Z^n}`
  (times `(2w−q)` for `E`). So `F′` is not a sum-rule functional, and its "light side" is undefined.
* **`F_U`** [D, U5; certified, `cert_FU_negative.py`]:
  * normal-form heavy action, `++`: `∫[ψ₂(2−uq) + uψ₃](1−uq)^{−2}P_J + (10𝒥 − 20)u³`;
  * `+−`: `∫[ψ₂(2−uq) − uψ₃](1−uq)²P^{(0,8)}_{J−4} + (10𝒥 − 190)u³`;
  * negative at `(++, J = 2, m ∈ {213/200, 6/5})`, certified;
  * low side with contacts (§0 item 5).

  Scan [NE] (`scan_FU.py`, `m ≤ 100`):
  * `F_U ++`: min `−3.33` at `m ≈ 1.065, J = 2`;
  * `F_U +−`: min `3.60·10⁻⁵` at `(21.7, 244)`;
  * reference `F_{g3}`: `++` positive (`b̂³`-normalised min 1.25 at `m ≈ 6.15, J = 26`), `+−` min `2.68·10⁻⁵` at `(21.7, 244)`.

## 6. Status

| statement | status |
|---|---|
| `𝓜 = t⁴f`, its cut discontinuities (full `d₀₀, d₄₄`, common prefactor) | [D] |
| definition, Regge convergence (`k ≥ 2`, `J₀ < 2`) of `B′_k` | [D] |
| low sides of `B′_k` (EFT + matter) | [D, exact] |
| stripped ↔ unstripped identities (kernels, densities, low sides) | [D, exact] |
| `F_{g3}`, M01, M02 as unstripped functionals; contacts cancel; `A, B` | [D, exact] |
| their heavy positivity | [certified] (existing Theorem 1 / Theorem R + exact identity) |
| equality of unstripped and stripped LP optima (same span) | [D] |
| `F′` not a sum-rule combination | [D] |
| `F_U` keeps contacts; `F_U < 0` at `(++, 2, 213/200)`, `(++, 2, 6/5)` | [D]; [certified] |
| `F_U` scans | [NE] |
| `C_k → ε_k^{−1}B′_k` exactly; `Ω^U_{g3}` is the correct AdS uplift | [I] |

## 7. Code (`ads4/code/flatunstripped/`)

* `sumrules.py`: U1–U6, exact sympy; ≈ 20 s.
* `cert_FU_negative.py`: Arb integration.
* `scan_FU.py`: [NE].
* `run_all.sh`: driver.
* `cache/`: logs.
