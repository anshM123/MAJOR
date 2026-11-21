# Result Q, Stage 0: flat D = 4 two-sided bounds on the quartic graviton contacts, and the minimal sum-rule set

*21 Nov 2026, Stage-0 agent. Unrefereed, not committed. No existing document or code was edited. New code: `ads4/code/resultq/`
(caches in `resultq/cache/`, git-ignored). Downloads (data only): scratchpad `…/scratchpad/resultq/{bcp,chz}` (arXiv sources of
1906.07195 and 2412.08689); the CHLPS source `…/scratchpad/a2201/gravity_4d.tex` was already there.*

Grades as in `RESULT_D_STATEMENT.md`: **[certified] [P] [P given X] [D] [D-sketch] [Q] [DV] [NE] [M] [I]**. A chain has the grade of
its weakest link. Every LP number here is floating point on a finite grid: **[NE]**. Nothing in this note is certified.

Conventions: CHLPS arXiv:2201.06602 (`FLAT_GRAVITON_BOUND.md` §1). `f_low = 8πG/(stu) + 2πG su|ĝ₃|²/t + g₄ + g₅t + g₆t² − g₆′su + …`,
`M` = mass of the lightest heavy state (spin-4 gap), `M = 1` unless written, `u = 1/m²`, `𝒥 = J(J+1)`. IR: `log(M/m_IR)`, translated
to AdS₄ as in Result D by `m_IR ↔ 3/(2R_AdS)`, i.e. `log(M/m_IR) → log(2Λ/3)` with `Λ = M R_AdS` [M, as in D].

---

## 0. Answer in brief

1. **Counting.** In D = 4 there are exactly **three real 8-derivative (`R⁴`-type) graviton contacts**: `g₄` (MHV, real;
   `= 8πG(α₄+α₄′)`) and `ĝ₄` (all-plus, complex; `= 8πG(α₄−α₄′+iα̃₄)`) [D; agrees with CHLPS (2.10)–(2.11) Q]. Among the forward
   terms, `g₄` appears only in `B₄⁽¹⁾|₀` and `ĝ₄` only in `B₄⁽⁴⁾|₀`; smeared `B₂⁽¹⁾, B₃⁽¹⁾` see no contact at all (§1).
2. **Minimal sum-rule set for a two-sided `g₄` bound** (all four needed, §2):
   **`S_min = { ∫ψ₂B₂⁽¹⁾ (smeared),  B₄⁽¹⁾|₀,  ∂_{p²}B₄⁽¹⁾|₀,  B₅⁽¹⁾|₀ }`**. Forward content: **`k = 4` (value and first derivative),
   `k = 5` (value)**. Adding smeared `∫ψ₃B₃⁽¹⁾` is optional and lowers `A` by 24 % (no light matter) or by a factor 2 (light matter
   allowed). `ψ₃(0) = 0` costs nothing. The roles of the four members:
   * `B₂` is the only carrier of `G` [P];
   * a `p = 0` term is unavoidable: forward-free smeared functionals cannot isolate `g₄` [P given `NOFORWARD_FUNCTIONAL.md` §1.2];
   * `∂B₄|₀` is the only admissible term that grows with `J` and has the same sign in both helicity channels: it is needed for
     positivity at `J → ∞`, fixed `m` [D], and the LP without it is infeasible [NE];
   * `B₅|₀` cancels `g₅`, which has no sign [D].
3. **The two-sided bound (flat, tree level, CHLPS hypotheses incl. the IR/eikonal one):**

   **`¼|ĝ₃|²M⁸ ≤ g₄M⁶/(8πG) ≤ A log(M/m_IR) + B − κ|ĝ₃|²M⁸`,  and  `|ĝ₄| ≤ g₄`**

   * Lower bound `¼` and `|ĝ₄| ≤ g₄`: [P given the forward sum rules]; the `|ĝ₄|` statement needs `B₄⁽⁴⁾|₀` in addition.
   * Upper bound [NE]:

   | family (`d = 12`, `(1−p)³` weights, no light matter unless stated) | margin `δ` | `A` | `B` | `κ` | dense off-grid check |
   |---|---|---|---|---|---|
   | `S_min` + `∫ψ₃B₃` | 0.05 | **11.45** | **−12.98** | 0.069 | pass (min 𝔉/N = 0.019) |
   | `S_min` alone (no `B₃`) | 0.05 | **15.00** | **−15.97** | 0.086 | pass (0.050) |
   | `S_min` + `B₃`, light matter allowed | 0.05 | 14.27 | −15.58 | 0.065 | pass (0.043) |
   | `S_min` + `B₃`, optimum | 0 | 10.11 | −11.48 | 0.059 | fails off-grid (marginal; limit value only) |
   | `S_min` + `B₃`, light matter, optimum | 0 | **12.38** | **−13.57** | 0.057 | limit value; **CHLPS (4.4): 12.3, −13.5** |

   In AdS variables, `g₄Λ⁶/(8πG) ≤ A log(2Λ/3) + B` (AdS units), i.e. `≈ 35 / 48 / 61` at `Λ = 100 / 300 / 1000` (robust full set).
4. **Validation.**
   * CHLPS's example `F_{g₄}`, with the `B₅|₀` completion, reproduces exactly: `A = 99/7 = 14.143`, `B = −15.507`, `κ = 0.04192`
     (paper: 14.1, −15.5, 0.04). It is positive on a dense grid (min 𝔉/N = 0.010 at `+−`, `J = 4`, `m = M`) [DV].
   * The CHLPS optimal bound `12.3 log − 13.5` (improved sum rules, `∂^{≤2}B₄`, light matter) is reproduced to 0.7 % in `A`
     (12.38, −13.57) by `S_min + B₃` alone [NE].
   * The `|ĝ₃|²` closure of the same LP with light matter lies at `|ĝ₃|²M⁸ ∈ (100, 110)` for `log = 5.30`. CHLPS (4.4) gives
     24.9·5.30 − 27.6 = 104 [NE].
   * The lower-bound slope `¼` is reproduced; CHLPS fit 0.26 [NE].
5. **Go/no-go: GO for Stage 1, conditional.**
   * The gate is passed: forward content `k ≤ 5`, `∂_{p²}^{≤1}`. It cannot be reduced further.
   * The condition is not reduced by Stage 0 either. The forward pieces carry the signal, and their AdS light side must be known
     **exactly**: an `O(G)` residue `(ν²+9/4)^{−4}` at `ν = 0` in `C⁺₄` has leverage `~Λ⁸` on the bound (`~10¹³` at `Δ_g = 100`,
     `~10²¹` at 1000; same scaling as `ADS4_LIGHT.md` §2.1).
   * So Stage 2's gate (residual leverage ≤ 10²) can only be met by an exact tree-level evaluation (T-exchange + contact Witten
     diagrams with vanishing (E)-type residues), not by estimates. This is H9a in a new role. D avoided it; Q cannot (item 2) [D].
   * A G-free sharper test (`g₆/g₄`, item 6) would need `k = 6` forward terms: **no-go** by the gate's own criterion.
6. **Virasoro–Shapiro.**
   * Type II VS has `g₄/(8πG) = 2ζ(3) = 2.404` [DV], `ĝ₃ = g₅ = ĝ₄ = 0`. Its first level is `J = 4` (`+−`) and `J = 0` (`++`) at
     `m = M`, with unit positive residue [DV].
   * It is **inside**, far from the upper boundary: `2ζ(3)/(A log(2Λ/3)+B) = 0.068 / 0.050 / 0.039` at `Λ = 100 / 300 / 1000`
     (robust full set; CHLPS optimum: 0.063 / 0.047 / 0.036).
   * VS would touch the boundary only at `Λ ≈ 5–6`. The bound grows like `log Λ` and VS does not; the upper bound is log-dominated
     and cannot be near-saturated by any theory with `Λ`-independent `g₄M⁶/G`.
   * By contrast, VS sits at **86 %** of the G-free forward bound `g₆M⁴ ≤ g₄` (`g₆M⁴/g₄ = ζ(5)/ζ(3) = 0.863`), with `g₆′/g₆ = 1`
     [DV]. That needs `B₆|₀` (`k = 6`).
7. **ABJM (IIA regime).**
   * Dictionary [Q]: BCP (1906.07195) and CHZ (2412.08689) give `L⁸/ℓ_s⁸ = 4π⁴λ²` (`R⁴/ℓ_s⁴ = 2π²(λ−1/24)`). The flat limit of the
     localization-fixed `R⁴` term is the IIA `ζ(3)ℓ_s⁶stu/32`.
   * With maximal SUSY the `R⁴`/Einstein ratio is the same for every component, hence `g₄/(8πG) = 2ζ(3)/M⁶`, `M² = 4/α′` [D given Q].
     The first massive level has twist `≈ 2(2π²λ)^{1/4}` (CHZ §2.2) [Q], so `Λ ≈ 4.22λ^{1/4}`.
   * Then `g₄Λ⁶/(8πG) = 2ζ(3)[1 + O(λ^{−1/4})]`: **ABJM's `λ^{−3/2}` is exactly the `Δ_g^{−6}` scaling, with coefficient 2ζ(3)**.
     It lies at 4–7 % of a Result-Q-type bound for `λ = 3·10⁵ … 3·10⁹`. **Pass, but a weak test** [D given Q; the AdS bound itself is [I]].
   * The genus-one term (`+(2π²/3)g_s²`, finite `μ = N/k⁵`) is outside A2. The M-theory regime is outside A1′/A2.
8. **Novelty.**
   * The flat `g₄` bounds themselves are **known** (CHLPS 2022). New here:
     * the minimal-set analysis and necessity arguments;
     * the observation that the CHLPS optimum is saturated by `S_min + B₃` to < 1 %;
     * the AdS-liftability gate assessment;
     * the explicit VS/ABJM placement in AdS variables.
   * No paper bounding AdS₄ graviton quartic contacts / ⟨TTTT⟩ contacts from dispersive CFT sum rules was found (abstract-level search).
     The closest analogue is AdS₅: Caron-Huot–Coronado–Zahraee et al. 2207.01615, 2412.00249 (scalar N = 4 correlators + localization,
     two-sided bounds, VS comparison) [Q, abstracts].

---

## 1. Quartic contacts in D = 4 and where they appear

**Counting [D].**
* In spinor form the Weyl tensor is `Ψ_{ABCD} ⊕ Ψ̄_{ȦḂĊḊ}`. Quartic scalars without derivatives:
  * `Ψ⁴`: the invariants of a binary quartic are generated by `I` (degree 2) and `J` (degree 3), so at degree 4 only `I²` exists:
    one complex coefficient;
  * `Ψ²Ψ̄²`: only `IĪ` (undotted and dotted indices cannot be contracted with each other): one real coefficient;
  * `Ψ³Ψ̄`: none.
* Total: **3 real**, i.e. `(R²)², (RR̃)², R²RR̃` ↔ `α₄, α₄′, α̃₄` ↔ `g₄ ∈ ℝ`, `ĝ₄ ∈ ℂ` (CHLPS (2.10)–(2.11) [Q]; this agrees with the
  general-D classification of Chowdhury et al. 1910.14392, where D = 4 identities remove the remaining structures [Q, abstract only]).
* `g₄` enters the MHV amplitude `f` as a constant; `ĝ₄` enters the all-plus `h` as `½ĝ₄(s²+t²+u²)²`. The `+++−` amplitude `g` has no
  8-derivative contact (its first contact is `ĝ₆″`).
* **SUSY (type II, ABJM): `ĝ₄ = 0`, `α₄ = α₄′`** (only MHV survives) [Q].

**Which sum rules see which contacts** (low sides re-derived in `FLAT_GRAVITON_BOUND.md` §1.3; general rule from `NOFORWARD_FUNCTIONAL.md`
§1.1: a degree-`n` contact enters `−B_k⁽¹⁾` as `c·p^{2(n−k+4)}`, so **`∂^j_{p²}B_k|₀` sees only contacts of degree `n = k − 4 + j`**):

| term | low side (tree, no matter) | contacts seen | forward order |
|---|---|---|---|
| `∫ψ₂B₂⁽¹⁾` | `16πG/p² + 2πG|ĝ₃|²p⁶` | none | smeared |
| `∫ψ₃B₃⁽¹⁾` | `−2πG|ĝ₃|²p⁴` | none | smeared |
| `B₄⁽¹⁾|₀` | `2g₄` | **`g₄` only** | `k = 4`, `∂⁰` |
| `∂_{p²}B₄⁽¹⁾|₀` | `4πG|ĝ₃|² + g₅` | `g₅` | `k = 4`, `∂¹` |
| `B₅⁽¹⁾|₀` | `g₅` | `g₅` | `k = 5`, `∂⁰` |
| `∂²B₄|₀`, `∂B₅|₀`, `B₆|₀` | `g₆ ± g₆′`, `2g₆` | `g₆, g₆′` | beyond the gate |
| `B₄⁽⁴⁾|₀` (fixed-`t`, `h`) | `ĝ₄` (CHLPS (2.40)) | **`ĝ₄`** | `k = 4`, `∂⁰`, non-MHV |
| smeared `B_{k≥4}` | all contacts of degree `≥ k − 4`, at `O(1)` | infinitely many | — |

Heavy actions per `|c^h_{J,m²}|²`, normal form `𝔉 = m⁴F` (`common_q.py`; `FLAT_GRAVITON_BOUND.md` (0.1), §1.2):
* `B₄|₀ → 2u²`;
* `∂B₄|₀ → (c₄ − 2𝒥)u³`, with `c₄ = 1 (++), 41 (+−)`;
* `B₅|₀ → σ_h u³`;
* the smeared parts as in `rcommon.py`.

## 2. The minimal family and why each member is needed

Family: `F = ∫ψ₂B₂ + ∫ψ₃B₃ − a(∂B₄|₀ − B₅|₀) − c B₄|₀`, `a ≥ 0`, normalised `c = ½`. Heavy positivity gives
`−F|_low = 16πG L(ε) − πGk|ĝ₃|² − g₄ ≥ 0`, with `k = −2μ₂₆ + 2μ₃₄ + 4a`, `L = ∫_ε¹ψ₂/p²`. Hence `A = 2ψ₂′(0)`, `B = 2L_fin`, `κ = k/8`.
The heavy action is `𝔉 = B₂-part + B₃-part + a(2𝒥 − c₄ + σ_h)u³ − u²`.

* **(N1) `∫ψ₂B₂` is necessary [P].** It is the only term whose low side contains `G` without `ĝ₃`. Without it no inequality relates
  `g₄` to `G`; only ratios such as `|g₅| ≤ g₄/M²` survive.
* **(N2) A `p = 0` term containing `g₄` is necessary [P given `NOFORWARD_FUNCTIONAL.md` §1.2].** The Müntz induction there uses only
  test contacts `w^ℓt^Lg(t)`, `t^ℓw^Lg(w)` of arbitrarily large degree. So a forward-free combination of smeared `B₄…B_K` that sees
  only finitely many contacts has `ψ_k = 0` for all `k ≥ 4`, and cannot see `g₄` alone. Within `k ≤ 5`, `∂ ≤ 1` the only forward term
  carrying `g₄` is `B₄|₀`. (`g₄ ≥ 0` itself follows from `B₄|₀` alone.)
* **(N3) `∂B₄|₀` is necessary [D; LP infeasible without it, NE].**
  * At fixed `m` and `J → ∞` every smeared part tends to 0. `P_J(1−2up²)` concentrates in `p ≲ m/J`, where `ψ₂ ~ p` gives `O(J⁻²)`;
    elsewhere it oscillates with amplitude `O(J^{−1/2})`. The `B₄|₀` term contributes the `J`-independent `−u² < 0`.
  * A positive term growing with `J` in both channels is therefore required. Among the forward terms with `k ≤ 5`, `∂ ≤ 1`, only `∂B₄|₀`
    qualifies (`2𝒥u³` in both channels); `∂B₅|₀` grows as `σ_h𝒥`, with opposite signs. Positivity then needs `𝒥 ≳ m²/(2a)`, which
    is automatic in the scaling region (`b̂ = 2J/m` fixed, `u → 0`), where `−u²` is subleading.
* **(N4) `B₅|₀` is necessary [D].** `∂B₄|₀` brings `+a g₅` into the low side. `g₅` has no sign (VS: 0; heterotic: < 0), and the only
  admissible term that removes it is `B₅|₀` with coefficient `−a`. Without it one bounds `g₄ + a g₅`, not `g₄`.
* **(O) `∫ψ₃B₃` is optional [NE].**
  * Without light matter: `A = 12.93` without `B₃` vs `10.11` with it (δ = 0); `15.00` vs `11.45` (δ = 0.05).
  * With light matter: `25.7` vs `12.4` (δ = 0).
  * `ψ₃(0) = 0` (needed in the forward-free class) changes nothing here.
  * AdS-style weights `(1−p)⁵`: `A = 11.16` (δ = 0), `12.59` (δ = 0.05, with `ψ₃(0) = 0`); `15.53` without `B₃`. That is a +10 % cost.
* **Lower side [P given the forward sum rules].** `B₄|₀`: `2g₄ = ⟨2(|c⁺⁺|²+|c⁺⁻|²)/m⁸⟩ ≥ 0`. Combining `∂B₄|₀` and `B₅|₀`:
  `4πG|ĝ₃|² = ⟨(−2𝒥|c⁺⁺|² + (42 − 2𝒥)|c⁺⁻|²)/m¹⁰⟩ ≤ 2⟨|c⁺⁻|²/m¹⁰⟩ ≤ 2g₄/M²` (`+−` has `J ≥ 4`, so `42 − 2𝒥 ≤ 2`). Hence
  **`g₄/(8πG) ≥ ¼|ĝ₃|²M²`** (Maldacena et al. 1407.5597 (6.13) [Q via CHLPS]; CHLPS fit 0.26).
  * The LP `slope.py` (`ψ₃` smeared, `ψ₂ = 0`, which is forced for a tangent through the origin, since `ψ₂ ≠ 0` implies
    `ψ₂′(0) > 0` and a log term) gives `0.2500` for `d₃ = 0…8` [NE]: `B₃` does not improve it.
* **`ĝ₄` [P given `B₄⁽⁴⁾|₀`].** By Cauchy–Schwarz, `|ĝ₄| = |⟨(c⁺⁺)²/m⁸⟩| ≤ ⟨|c⁺⁺|²/m⁸⟩ ≤ g₄` (CHLPS (2.41)). So the upper bound on
  `g₄` bounds all three `R⁴` couplings two-sidedly: `α₄, α₄′ ∈ [0, g₄/(8πG)]`, `α̃₄² ≤ 4α₄α₄′`.
  * This needs the fixed-`t` all-plus rule `B₄⁽⁴⁾`, a "non-shock" rule that CLS-D do not construct (`LIT_HELICITY_MIXING.md` §1) [I for AdS].

## 3. LP results (`lp_g4.py`, `dense_check.py`; all [NE])

* **Method.** As in `robust/lp.py` and `nf_lp.py`.
  * Weights: `ψ₂ = (1−p)^r Σ x_i p^{i+1}`, `ψ₃ = (1−p)^r Σ y_j p^j`.
  * Positivity `𝔉 ≥ δN` (`N = |B₂|+|B₃|+|FW|+u²`) on ≈ 2·10⁵ physical states (`u ∈ [2.5·10⁻⁵, 1]`, `J ≤ 2500`, `b̂ ≤ 200`), by cutting
    planes over the sign patterns.
  * Scaling limit: `2∫ψ₂J₀(b̂p) ≥ δρ(b̂)`, `b̂ ≤ 400`.
  * Lemma-M asymptotic rows for `b̂ ∈ [150, 10⁴]`, `u ≤ ¼`.
  * Optional light-matter rows (as CHLPS App. C).
  * Objective: minimise `A` (tie-break `B`).
* **Dense off-grid check.**
  * Set-up: 160 `u` values disjoint from the LP grid, `J ≤ 3000` / `b̂ ≤ 400`, and the scaling limit to `b̂ = 2000`.
  * All δ = 0.05 functionals pass; the δ = 0 optima fail by up to `−0.3N` at `J ≈ 2900`, `b̂ ≈ 50`. So δ = 0 numbers are limits,
    not functionals.
  * The coverage gap is `u < 10⁻⁴` with `b̂ ∈ (25, 150)`, which only the scaling-limit margin controls [NE].
* **Convergence in degree** (δ = 0 / 0.05, no matter):
  * `A = 10.68 / 12.04` (d = 4), `10.34 / 11.69` (6), `10.20 / 11.54` (8), `10.13 / 11.47` (10), `10.11 / 11.45` (12);
  * with `(1−p)²`: `10.03` (d = 8, δ = 0); with `(1−p)⁵`: `11.16`;
  * with matter: `12.89 → 12.58 → 12.46 → 12.40 → 12.38`.
* **`|ĝ₃|²` dependence** (`lp_g4.py --lref L --G3 g`, δ = 0.05, d = 8, no matter). `U(g)` is the best upper bound on `g₄M⁶/(8πG)` at
  `|ĝ₃|²M⁸ = g`:

  | `log(M/m_IR)` | `U(0)` | `U(20)` | `U(50)` | closure (`U < ¼g`, actually `U < 0`) |
  |---|---|---|---|---|
  | `log 200 = 5.30` (`Λ = 300`) | 47.7 | 46.2 | 43.8 | between `g = 60` and 70 |
  | `log 666.7 = 6.50` (`Λ = 1000`) | 61.7 | 60.2 | 58.1 | between `g = 70` and 80 |

  * The boundary is nearly flat (`κ ≈ 0.07`) until it collapses at the forward-term `|ĝ₃|²` bound (with light matter:
    `(100, 110)` vs CHLPS 104).
  * In AdS language the same forward terms that bound `g₄` would bound `|ĝ₃|²Λ⁸ ≲ 60–80`. That is about 150× below Result D's
    forward-free `≈ 10⁴`, but it is exactly the H9a-dependent route D chose not to take [D].
  * D's bound on `ĝ₃` does not move the `g₄` window appreciably: `κ|ĝ₃|²` ≤ a few units inside the closure.

## 4. Go/no-go gate (plan §4 of `UNIFIED_SCOPING.md`)

| criterion | result | grade |
|---|---|---|
| forward terms beyond `k = 5`? | **no**: `k = 4` (`∂⁰, ∂¹`), `k = 5` (`∂⁰`) suffice and are necessary | [D] + [NE] |
| `∂^{≥2}_{p²}`? | **no**; `∂²B₄` improves `A` by < 1 % (CHLPS's optimum vs ours) | [NE] |
| `ĝ₄` sector | needs fixed-`t` `B₄⁽⁴⁾|₀` (all-plus amplitude): not in CLS-D's shock family; ABJM has `ĝ₄ = 0`, so not needed for the test | [I] |
| leverage of the forward pieces (`NOFORWARD`/`ADS4_LIGHT` measure) | an `O(G)` residue `16πG(ν²+9/4)^{−4}` in `C⁺₄(0)` shifts `g₄M⁶/(8πG)` by `≈ 2c(4/9)⁴Λ⁸ = 0.039Λ⁸` (up to the `O(1)` flat-limit normalisation `ε₄`): `4·10¹⁴` (Λ = 100) against a signal `≈ 35`, i.e. leverage `~10¹³`; for `∂B₄|₀`, `~Λ¹⁰` | [D-sketch] |
| can the forward pieces be avoided? | **no** (N2) | [P given NOFORWARD §1.2] |

**Verdict: GO (conditional).**
* The formal gate is passed, and the forward content is minimal and fixed.
* The real risk moves entirely to Stages 1–2. The subtracted `k = 4, 5` CLS-D functionals at `ν = 0` (including `∂_{ν²}`) must have
  light sides that are **exactly** known at tree level:
  * `T`-exchange: zero, or computed;
  * `[TT]` content: entirely T-exchange + contact Witten diagrams;
  * nothing else under A1′.
* No numerical light-side error budget can survive a `10¹³` leverage. If Stage 2 cannot show exactness, record the negative result
  ("forward-free spinning sum rules cannot bound quartic contacts in AdS₄; forward ones need exact `ν = 0` light sides").
* Caveat [I]: (N3) is a flat-space statement. In AdS the large-`J` region at fixed `m` is the Regge/large-`b̂` region, where D needed
  curvature (`H_R`). Whether `∂C⁺₄(0)` keeps its `2𝒥` growth with the right sign there is Stage 1/3 work.

## 5. Virasoro–Shapiro placement (`vs_check.py`)

* `f_VS = 8πG/(stu)·Γ(1−s)Γ(1−t)Γ(1−u)/(Γ(1+s)Γ(1+t)Γ(1+u))` (`M² = 4/α′`), giving `g₄/(8πG) = 2ζ(3) = 2.40411` (to 10⁻⁷),
  `g₅ = 0`, `g₆ = g₆′ = 2ζ(5)·8πG` [DV].
* The residue at `s = M²` is `t`-independent and equal to `+8πG`, i.e. pure `+−`, `J = 4`; by symmetry the `t = M²` residue is pure
  `++`, `J = 0`. VS is a legitimate point with gap `M` [DV].

| bound | Λ = 100 | 300 | 1000 | VS on boundary at |
|---|---|---|---|---|
| robust `S_min + B₃` (11.45, −12.98) | 35.1 (VS 6.8 %) | 47.7 (5.0 %) | 61.5 (3.9 %) | Λ = 5.75 |
| robust `S_min` (15.00, −15.97) | 47.0 (5.1 %) | 63.5 (3.8 %) | 81.6 (2.9 %) | 5.10 |
| optimum `S_min + B₃` (10.11, −11.48) | 31.0 (7.8 %) | 42.1 (5.7 %) | 54.3 (4.4 %) | 5.92 |
| CHLPS (4.4), with matter | 38.2 (6.3 %) | 51.7 (4.7 %) | 66.5 (3.6 %) | 5.47 |

* **VS is allowed, far from the upper boundary.** It is allowed on the lower side too: `g₄ > 0 = ¼|ĝ₃|²`, `ĝ₄ = 0 ≤ g₄`.
* The distance is controlled by the IR log. A G-normalised bound cannot be near-saturated by a theory whose `g₄M⁶/G` is
  `Λ`-independent. In CHLPS's plots the superstring also sits well inside the `(|ĝ₃|², g₄)` region [Q].
* Where VS **is** near a boundary: the G-free ratio `g₆M⁴/g₄ = ζ(5)/ζ(3) = 0.863 ≤ 1` (forward `B₆|₀` vs `B₄|₀`), and
  `g₆′/g₆ = 1 ∈ [0, 3]` (CHLPS (4.10)). These need `k = 6` forward terms (outside the gate).

## 6. AdS₄/CFT₃ dictionary and ABJM

* **Observable.** Each of the 3 bulk `R⁴` contacts gives a contact Witten diagram in `⟨TTTT⟩`. Under A1′/A2 it contributes only to
  `[TT]` double-trace data at `O(1/c_T)`; in Mellin space it is a polynomial of top degree 4.
  * Its flat-space (bulk-point) limit reproduces `g₄⟨23⟩⁴[14]⁴` (resp. the `h` structure) [D-sketch: Penedones' formula / CLS-D crystal
    for spinning correlators; the `⟨TTTT⟩` normalisation constant, a Γ-function factor depending on the Mellin degree, is not computed] [I].
  * The number of independent leading-derivative AdS₄ quartic graviton contacts equals the flat count 3 (2 parity-even, 1 parity-odd)
    [D-sketch: leading term of each contact ↔ polynomial flat S-matrix].
  * Dimensionless target: `g₄R⁶/(8πG)·Λ⁶ ≤ A log(2Λ/3) + B` [I until Stages 1–3].
* **ABJM numbers** (BCP 1906.07195 §4, eqs. (class), (IIAAMP), (IIARes1); CHZ 2412.08689 eq. (1.1)) [Q]:
  * `L⁸/ℓ_s⁸ = 4π⁴λ²`, i.e. `L²/α′ = π√(2λ)`;
  * `g_s² = 512λ²/(3c_T)`;
  * the 't Hooft expansion `M = c_T⁻¹[B¹_SG M_SG + λ^{−3/2}(… + B⁴₄M₄) + …]`, `B⁴₄/B¹_SG = 105ζ(3)/(64√2π³)`, fixed from localization (the
    finite-`μ` result (IIAResAgain), whose `μ → 0` limit matches the flat IIA `ζ(3)ℓ_s⁶stu/32`);
  * massive string states at level `δ`: `Δ(Δ−3) = R²m² = 4δR²/α′` and twist `τ = 2ν^{1/4}√δ + O(1)`, `ν = R⁴/ℓ_s⁴ = 2π²(λ − 1/24)`
    (CHZ eq. (class) and §2.2), so the lowest string twist is `Δ_g ≈ 2(2π²λ)^{1/4} = 4.22λ^{1/4}`.
* **Mapping [D given Q].**
  * Maximal SUSY: the 4D graviton MHV amplitude equals the 10D superamplitude restricted to 4D (zero modes on CP³; KK modes do not couple
    to two zero modes at tree level, which is A1′). So `f = 8πG/(stu)·(1 + ζ(3)ℓ_s⁶stu/32 + …)`, i.e. `g₄/(8πG) = 2ζ(3)(α′/4)³ = 2ζ(3)/M⁶`
    with `M² = 4/α′`, the first massive level (which contains the `J = 4`, `+−` state).
  * Hence, in Result-Q variables, `g₄R⁶Λ⁶/(8πG) = 2ζ(3)·(1 + O(1/Λ))`, `Λ = MR ≈ 2^{5/4}π^{1/2}λ^{1/4} = 4.22λ^{1/4}`. **The `λ^{−3/2}`
    suppression is precisely `Λ^{−6}`; the ratio to the gap-normalised bound is the VS ratio of §5.**
  * Λ = 100, 300, 1000 ↔ `λ ≈ 3.2·10⁵, 2.6·10⁷, 3.2·10⁹`; ABJM sits at 6.8 %, 5.0 %, 3.9 % of the robust bound.
* **Is ABJM inside a Result-Q-type bound for its gap?** Yes, with a factor 15–25 to spare at the flat level. The `O(1)` shifts are
  negligible against that margin [D given Q]:
  * `Λ` vs `Δ_g` (twist of the `J = 4` member ≈ `ML − 2.5`);
  * `λ → λ − 1/24`;
  * curvature corrections (CHZ's first curvature correction to `R⁴`, `O(λ^{−1/2})` relative).
  * It would fail only for `Λ ≲ 6` (`λ ≲ 4`), outside the IIA regime.
* **Finite `μ = N/k⁵` [Q].** The `R⁴` coefficient becomes `2ζ(3) + (2π²/3)g_s²` (genus one), with `g_s⁴ = 32π²μ`. This is a string
  loop: outside A2's tree-level positivity [I]. (Numerically it would reach the bound only at `g_s² ≈ (3/2π²)(bound − 2ζ(3)) ≈ 5–9`.)
  The M-theory regime (fixed `k`, no string gap) is outside A1′/A2.
* **Honest assessment.** The ABJM check is non-trivial in kind (`R⁴ ≠ 0`, unlike `t₄`) but weak in strength. It tests the sign
  (`g₄ > 0`) and the `Λ⁻⁶` scaling, with a log-dominated margin. A sharp test needs a log-free observable. Options:
  * a G-free ratio (`g₆/g₄`, where VS sits at 86 %; CHZ fix the ABJM `D⁴R⁴` at finite curvature [Q, abstract]), which needs `k = 6`
    forward terms;
  * an AdS bound whose `B` is far more negative.

## 7. Literature and novelty (abstract level except CHLPS, BCP, CHZ §1)

| paper | relation |
|---|---|
| CHLPS 2201.06602 (read §2–§4) | flat D = 4 two-sided `g₄` bounds: example `14.1 log − 15.5 − 0.04|ĝ₃|²`, optimum `12.3 log − 13.5` (improved sum rules, `∂^{≤2}B₄`, matter), lower `g₄/8πG ≥ 0.26|ĝ₃|²M²`, `|ĝ₄| ≤ g₄`; VS/heterotic/bosonic placed in their plots. **Reproduced (§0 item 4).** Their AdS₄ remark (`log(MR)`) is asserted, not derived. |
| Bern–Kosmopoulos–Zhiboedov 2103.12728 | D = 4, G-free bounds on ratios (`R³` vs `R⁴`, islands, low-spin dominance); string/loop models in small islands. Complementary (no G). |
| Guerrieri–Penedones–Vieira 2102.02847; Guerrieri–Murali–Penedones–Vieira 2212.00151 | D = 9–11, SUSY, lower bound on the `R⁴` coefficient α vs gap; string/M-theory near the minimal α. Different D, SUSY, non-perturbative. |
| Chowdhury et al. 1910.14392 | classification/counting of four-graviton contacts in all D (used for the count) |
| Caron-Huot, Coronado, Trinh, Zahraee 2207.01615; + 2412.00249 | **closest analogue**: AdS₅ N = 4 scalar correlators, dispersive sum rules + integrability + localization, two-sided bounds, VS flat-space connection. Not AdS₄, not gravitons/⟨TTTT⟩. |
| BCP 1906.07195; CHZ 2412.08689 | the ABJM `R⁴` data used in §6 (read the relevant sections) |
| Xu 2606.19283 | 10D SUSY dispersive bootstrap of VS: bounds on leading coefficients normalised by G, island around VS with a nonlinear ansatz |
| Agmon–Chester–Pufu (N = 8 numerical bootstrap, M-theory `R⁴`) | quoted in `UNIFIED_SCOPING.md`, not re-checked; no gap, SUSY input |

**Novelty verdict.**
* The flat bound is not new.
* New, at the level of a technical note:
  * (i) the minimal-set theorem-sketch (N1–N4) and the fact that `S_min + B₃` already saturates the CHLPS optimum;
  * (ii) the AdS-lift gate analysis with the leverage estimate;
  * (iii) the clean statement that ABJM's IIA `R⁴` is VS in gap units and sits at 4–7 % of the bound.
* A genuinely new physics result requires Stages 1–3, i.e. the first dispersive CFT₃ bound on `⟨TTTT⟩` quartic contacts. None was found
  in the literature [I, abstract-level search only; arXiv API queries listed in the hand-over].

## 8. Grades of the claims

| claim | grade |
|---|---|
| 3 real `R⁴` couplings in D = 4; `ĝ₄ = 0` for SUSY | [D] + [Q] |
| forward-order table (§1) | [D] (from CHLPS low sides + `NOFORWARD` §1.1) |
| (N1) | [P] |
| (N2) | [P given NOFORWARD §1.2] |
| (N3), (N4) | [D]; (N3) also [NE] (LP infeasible) |
| `g₄/(8πG) ≥ ¼|ĝ₃|²M²`; `|ĝ₄| ≤ g₄` | [P given the forward sum rules] (which carry CHLPS's forward-limit/loop caveat) |
| upper bounds `A, B, κ` | [NE] (grid LP + dense check; flat IR hypothesis as in `FLAT_GRAVITON_BOUND.md` Thm 3) |
| reproduction of CHLPS `F_{g₄}` constants | [DV] (exact); positivity [NE] |
| reproduction of CHLPS optimum 12.3/−13.5 and `|ĝ₃|²` closure ≈ 104 | [NE] |
| AdS form `log(2Λ/3)` | [M] (as in D) |
| leverage `~Λ⁸` of forward pieces | [D-sketch] |
| GO verdict | [D] (conditional on exact `ν = 0` light sides: [I]) |
| VS values, first-level content | [DV] |
| VS placement | [NE] (inherits the bound's grade) |
| ABJM `g₄/(8πG) = 2ζ(3)/M⁶`, `Λ ≈ 4.22λ^{1/4}` | [D given Q] |
| ABJM inside a Result-Q bound | [D given Q] for the flat image; the AdS bound itself [I] |
| novelty | [I] (abstract-level search) |

## 9. Hand-over

**Files** (`ads4/code/resultq/`; run with `OMP_NUM_THREADS=2`):

| file | content |
|---|---|
| `common_q.py` | grid of heavy actions (imports `robust/rcommon.py` read-only), AdS log integral, matter rows, `heavy_action` |
| `lp_g4.py` | LP. `python3 lp_g4.py d2 d3 r [--delta δ] [--no-b3] [--no-fwd] [--psi3zero] [--matter] [--lref L --G3 g]`; appends JSON to `cache/lp_g4.jsonl` |
| `dense_check.py` | off-grid check of stored functionals (`python3 dense_check.py <line numbers>`) |
| `val_fg4.py` | CHLPS `F_{g₄}` constants + positivity |
| `slope.py` | lower-bound slope LP |
| `vs_check.py` | VS coefficients, residues, placement table |
| `summarize.py`, `save_functionals.py` | list LP runs; write `cache/functionals_q.json` (reference functionals: weights, `a`, `c`, `A`, `B`, `κ`) |
| `batch1.sh` | the degree/variant scan (≈ 10 min) |

Grids take seconds to minutes and are cached as `cache/grid_*.npz` (≈ 20 MB each; delete freely).

**Next steps (Stage 1).**
1. Identify the CLS-D subtracted `k = 4` (`j′` even) and `k = 5` (`j′` odd) spinning functionals whose crystal images are
   `B₄⁽¹⁾|₀`, `∂_{p²}B₄⁽¹⁾|₀` and `B₅⁽¹⁾|₀`. Check flip-freeness (as `EE_FOUNDATIONS.md` did for `k = 2`).
2. Before any heavy work, decide (E) at `ν = 0` for these: either exact vanishing of the `T`-exchange (Regge-spin counting, `J₀ ≤ 2`
   with one unit of margin is not enough at `k = 4` unless the subtraction is exact), or an exact Witten-diagram evaluation. The leverage
   (§4) leaves no room for estimates.
3. AdS large-`J` behaviour of `∂C⁺₄(0)` (does the flat `2𝒥u³` growth survive curvature with the right sign?).
4. Sharper ABJM test: consider whether a `k = 6` forward term is liftable after all (G-free ratio `g₆/g₄`, VS at 86 %).

**Not done.**
* No certificate of any upper-bound functional; the δ = 0.05 ones pass the dense check only.
* No full-text reading of BKZ, GPV, Xu, or of 2412.00249.
* The `⟨TTTT⟩` Mellin normalisation of the `R⁴` contact.
* `ĝ₄` LP (needs `B₄⁽⁴⁾`; the bound `|ĝ₄| ≤ g₄` is analytic).
* AdS-curvature corrections to the ABJM `R⁴` coefficient (CHZ) — irrelevant at the present margin.

**arXiv queries used** (export API, 21 Nov 2026): ids 1906.07195, 2207.01615, 2412.00249, 2212.00151, 2102.02847, 1910.14392,
2103.12728, 2412.08689, 2606.19283, 2011.05728; searches "dispersive CFT sum rules", "ABJM AND gap AND bounds", "AdS bulk locality",
"graviton AND quartic AND causality AND bounds", title "Virasoro-Shapiro", "graviton AND AdS4 AND dispersive" (no AdS₄ graviton-contact
bound found).
