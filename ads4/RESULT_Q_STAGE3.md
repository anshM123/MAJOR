# Result Q, Stage 3: the validity input V₄, and heavy-side positivity of the forward functionals near the gap

*21 Nov 2026, Stage-3 agent. Unrefereed, not committed. No existing document or code was edited. New code:
`ads4/code/resultq/stage3_*.py` (+ `stage3_run.sh`, `stage3_after.sh`); outputs in `resultq/cache/stage3_*` (git-ignored). Read-only
imports: `exactskel/exactmom.py` (exact CHMRS Mellin moments), `resultq/common_q.py`, `stage12_forward_regge.py`, `robust/rcommon.py`.
`cache/lp_g4.jsonl` received three appended lines (runs of the existing `lp_g4.py`, `δ = 0.1, 0.2, 0.3`). The CLS-D source (2311.04271, scratchpad `dict/src`) was re-read for the definition of `Ψ⁺`, the subtraction factors and the
commutator/double-commutator rewriting. Compute: ≈ 2.5 CPU-hours on 2 cores.*

Grades: **[P]** proved, **[D]** derived with the stated assumption, **[DV]** derived and numerically verified, **[M]** numerical/model,
**[I]** open. (Stage 0/12 also used [Q], [NE], [D-sketch]; kept where quoted.) "Flat" = D = 4 amplitudes; "AdS" = the CFT₃ functionals.

---

## 0. Answer in brief

1. **V₄ (Task 1): holds for the tree correlators of Lemma A, grade [D]. No counterexample.** (§1)
   * **Lemma B [P].** `f⁺₄;₂ = f⁺₂;₂·R` and `f⁻₅;₂ = f⁻₃;₂·R`, with one common `R = ((v′−u′)² − 1)/(4u′v′)`.
     * `R` is bounded at every boundary of the Result-D K2 table other than the light-cone edge, and `→ 16η²r²` at the Regge corner.
     * It is singular only at the light-cone edge `u′ → 0` (`C(r)/s²`) and at the OPE corner (`−(1/z+1/z̄)/2`).
     * It is independent of `x₀`, so `EE`'s `x₀₂²`-edge is inherited unchanged from `k = 2`.
     * The forward pieces have no `ν`-window, so `H_edge` does not arise for them.
   * **Lemma C [P].** The poles of `f` carry the same Wightman `iε` in both orderings of the commutator, and `f` is analytic off the real
     axis along every light ray. Consequences:
     * Lemma A's pointwise vanishing holds as a distributional identity;
     * the vacuum-annihilation contour steps survive with `f`;
     * at the edge, the poles produce only the finitely many light-family terms (identity, `T`, leading `[TT]`), i.e. the analytic
       continuation in `Δ` [D].
   * **MFT (Task 1a).**
     * Direct (commutator form): **exactly 0**, for every `ν` and any `f` [P]. The light-ray-polarised `⟨T₃T₁⟩` numerator is
       `O(ε²)`, so the commutator is `O(ε⁴)`; the numerical pairing goes `3·10⁻⁵ → 2·10⁻¹¹` as `ε = 0.1 → 0.003`.
     * Swapped (s-channel) form: reduced to "identity + `[TT]₀` terms on the edge `u′ = 0`, all heavy `[TT]_{n,J}` exactly 0", but
       **not evaluated** (needs the spinning kernel's edge Taylor coefficients; no engine exists) [I].
   * **New caveat, beyond V₄ proper (non-perturbative validity/swap).**
     * In the large-`η` Regge corner (`r → 0` at fixed `ρ = 2rη`), `R → 4ρ²/(1−ρ²)²` gives no extra suppression.
     * The **raw** forward `Ĉ₄(0)`, `∂Ĉ₄(0)` have kernel growth `η^{5/2}` there. Under KKSDZ's uniform Regge bound they then need
       `J₀ < ½` ✗ (same failure as raw `Ĉ₂` in Result D) [D].
     * Raw `Ĉ₅(0)` is fine (`J₀ < 3/2`). A de-shifted pure-moment `Ĉ₄` would need `J₀ < 7/2`.
     * **The 4-shift de-shift is therefore load-bearing for validity** [I]. For the tree correlators themselves the explicit
       impact-parameter decay suffices [D-sketch].
2. **Stage 3 (Task 2): heavy positivity survives near the gap in the hybrid-exact scalar-skeleton model [M]; the graviton-specific part is [I].** (§2)
   * Model: exact CHMRS Mellin moments for the forward `ν = 0` pieces (all orders in `1/m`, `Δ_φ = 3`); hybrid flat + curvature for the
     smeared parts. Coverage: 298 / 274 / 283 (state, channel) cases at `Δ_g = 100 / 300 / 1000`, `J ≤ 24`, `m` from the gap to `4Λ`.
   * **Minimum heavy image of the Stage-0 robust functional** (`S_min + B₃`, `δ = 0.05`):
     **`min FF/N = +0.0478 / +0.0499 / +0.0500`**.
     * Location: `+−`, `J = 5`, `m = 1.2Λ`, `b̂ = 9.2` at 100 and 300; the `++`, `J = 0` gap state at 1000.
     * With the exact `k = 2` skeleton remainder (checked at 100): `+0.048`, i.e. ±0.005.
     * This is the flat LP margin itself: AdS curvature costs ≤ 0.003 `N`.
     * `b̂ ≥ 20` rows: `≥ +0.83 / +0.51 / +0.33`, consistent with Stage 12's leading-Regge coverage.
   * **Robustness is the issue, not the sign.**
     * The margin `0.05 N` is below the size of the graviton corrections Result D measured near the gap for `Ω^EE` (`0.1–0.5 N_s`,
       different normalisation). These cannot be computed at `ν = 0` with any existing engine [I].
     * A margin-0.2 re-optimisation (`δ = 0.2`, dense check passed, min `0.185`) has AdS min `+0.1985 / +0.2003 / +0.2001`, at a cost
       of `A`: 11.45 → 17.24.
   * **Lower bound:** per-state slope `1/(4c_L) = 0.2539 / 0.2548 / 0.2519 ≥ ¼` (binding: the `+−`, `J = 4` gap state). So
     **`¼|ĝ₃|²Λ⁸ ≤ g₄Λ⁶/(8πG)` survives [M]**, up to the uncomputed `b₀` [I].
   * **AdS constants** (functional unchanged, exact pole):
     * `A = 11.448` (exact);
     * `B_eff = 2L_AdS − A log(2Λ/3) = −12.804 / −12.938 / −12.970` (flat −12.979);
     * `κ = 0.0686` (+ `O(Λ⁻²)` and `b₀` [I]);
     * `2L_AdS = 35.10 / 47.66 / 61.45`.
   * `b₀`, `b₅`: not computed (not cheap) [I].
3. **ABJM:** `2ζ(3) = 2.404` sits at **6.8 % / 5.0 % / 3.9 %** of the AdS upper bound (robust `δ = 0.05`), and at 4.5 / 3.3 / 2.6 % of the
   `δ = 0.2` one. It is inside, far from the boundary, and above the lower bound (`ĝ₃ = 0`). This is unchanged from Stage 0/12: a weak,
   log-dominated test [D given Q; bound M].

---

## 1. Task 1: V₄ (validity of the subtracted `k = 4, 5` sum rules on `J₀ < 4` tree correlators)

### 1.1 What V₄ has to cover
`Ψ_{k,f}[𝒢] = Ψ⁺[f𝒢] − (1↔3)` (CLS-D eq. `spinningfunctional_definition` + §3.4.2). The vanishing argument (Lemma A) uses the
**commutator (Wightman) form** `∫ K f ⟨T₄[T₃,T₁]T₂⟩` (pointwise zero: 1, 3 lie on two null rays from `x₀`, hence `x₁₃² > 0`). The sum
rule used on the heavy side is the **double-commutator form** `Ψ⁺[f𝒢] ∝ ∫ K f ⟨[T₄,T₃][T₁,T₂]⟩`, expanded in the s-channel `(12)(34)`.
Passing from one to the other needs (a) the light-ray contour arguments ("light-transformed operators annihilate the vacuum") with `f`
inside, and (b) integrability at every boundary. New relative to Result D's `k = 2, 3` (K2): the factor `1/(u′v′)` with poles on
`x₁₂² = 0, x₃₄² = 0` (the edge of the double-commutator support), `x₁₄² = 0, x₂₃² = 0`.

### 1.2 Lemma B: one common factor (`stage3_v4.py` (1)–(3)) [P]
**`f⁺₄;₂ = f⁺₂;₂·R` and `f⁻₅;₂ = f⁻₃;₂·R` with the same `R(u′,v′) = ((v′−u′)² − 1)/(4u′v′)`** [P, sympy]. In the CHMRS frame
(`r = |x||y|`, `η`, `s = 1 − 2rη + r²`, `s̄ = 1 + 2rη + r²`): `R = 16r²(η²(1+r²)² − 4r²)/(s s̄)²` [P]. Its limits:

| region | `R` | consequence |
|---|---|---|
| Regge corner, fixed angle (`r → 0`, `η` fixed) | `16η²r²` | +2 units of Regge margin over `k = 2, 3` (CLS-D's `J₀ < 4`) [P] |
| single light-ray endpoints (`x₃ → x₀`; `x₃ →` null infinity) | finite limits (300 random configurations, linear approach) | `k = 4, 5` endpoint rows = `k = 2, 3` rows × bounded factor [DV] |
| light-cone edge `u′ → 0` (`s → 0`, the support edge of `[T₁,T₂]`, `[T₄,T₃]`) | `C(r)/s²`, `C = (1−r²)²(r⁴+6r²+1)/(1+r²)² > 0` | the only new singularity; `sup|s²R| = 1.079` on the region [P + DV] |
| OPE corner (`u′ = zz̄ → 0`, `v′ → 1`; `r → 1`) | `−(1/z + 1/z̄)/2` | one extra power, light s-channel terms only (§1.4) [P] |
| `EE`'s `x₀₂²`, `x₀₄²` edges | `R` depends on `(r, η)` only, not on `x₀` | the `x₀`-edge integrability is that of `k = 2` (K2) times a bounded factor [P given K2 (D)] |
| **large-`η` Regge corner** (`r → 0` at fixed `ρ = 2rη < 1`) | **`4ρ²/(1−ρ²)²`: bounded, no gain** | see §1.5: the only place where `k = 4` has no more margin than `k = 2` [P] |

At `ν = 0` (and for `∂_{ν²}` at 0) the forward functionals have **no `ν`-window**, so the heavy-side `H_edge` window `x = ν/m ∈ (0.8, Λ/m]`
of Result D does not arise for them at all [P]. It enters Result Q only through the smeared `k = 2, 3` parts, unchanged from Result D
(and suppressed by the `(1−p)³` weights at `p → 1`).

### 1.3 Lemma C: the `iε` of the subtraction poles (`stage3_v4.py` (4)) [P]
* **Same `iε` in both orderings.** In the two Wightman orderings of `⟨4[3,1]2⟩` (4312, 4132) the pairs `(12), (34), (14), (23)` — the
  poles of `f` — have identical relative order; only `(13)` differs, and `x₁₃² > 0` on the whole domain. Giving `f` the Wightman `iε`
  of the ordering it multiplies is therefore unambiguous, and `f𝒢` in both orderings are boundary values from one common domain. The
  pointwise commutator identity of Lemma A holds as an identity of distributions including the poles of `f`.
* **Tube analyticity.** Along a light ray `x₃ = x₀ + βz₃` every `x_{3j}²` is linear in `β` with real coefficients, so `f` is a rational
  function of `β`, analytic off the real axis and bounded at `|β| → ∞` (it tends to a constant, §1.2 row 2). Equivalently `x² ≠ 0` for
  `Im x ∈ ±V⁺`. Hence the contour arguments that delete the vacuum-annihilated terms work with `f` exactly as without it (convergence at
  the ends is the endpoint row of §1.2).
* **The double commutator acts on `f` too.** In `[4,3][1,2] = 4312 − 4321 − 3412 + 3421` the `(12)` and `(34)` prescriptions flip
  independently. A term without `x₁₂, x₃₄` (e.g. the Wick term `⟨41⟩⟨32⟩`) acquires `W·Δ₁₂Δ₃₄f`, a distribution supported on
  `x₁₂² = x₃₄² = 0` (`u′ = 0`). This is the position-space form of CLS-D's "finitely many non-vanishing double traces".

### 1.4 The light-cone edge and the OPE corner [D]
* Near `u′ → 0` the s-channel OPE converges (light-cone limit, `v′` away from 0 and ∞). The block expansion therefore commutes with the
  distributional pairing there.
* For one block, the double commutator is `(2i sin πa)²|x₁₂²|^a_+|x₃₄²|^a_+` with `a` shifted by −1 by `f`. The analytic continuation
  in `Δ` of the absolutely convergent integral (what CLS-D's per-block formula computes) **equals** the `iε` distribution
  (Gel'fand–Shilov: `sin(πλ)x₊^λ` is entire in `λ`) [P, standard].
* Only light s-channel terms are non-integrable at the edge: identity (`τ = 0`), `T` (`τ = 1`), the leading `[TT]` family (`τ = 2`, now a
  single instead of a double zero). Heavy states (`τ ≥ Δ_g`) are integrable with margin `≈ Δ_g`.
* So the edge produces **no boundary term beyond the finitely many light-family terms**, which are part of the light side by
  definition [D]. The OPE corner is the same mechanism with one power (`R ~ −(1/z+1/z̄)/2`) [D].

### 1.5 The large-`η` Regge corner: power counting (new) [D]
Leading-Regge kernels (Stage 12 §1.2, §1.4): raw `Ĉ_k ~ r^{a}κ(η)` with `a = 2` for `k = 4, 5` (`a = 0` for `k = 2, 3`). `κ(η)` is the
growth of the raw heavy image: `P_{5/2}` (`k = 4`, `η³ ×` stripped), `P_{3/2}` (`k = 5`), `ηP_{l±1} ~ η^{3/2}` (raw `k = 2`), `P_{−1/2}`
(`k = 3`, and any de-shifted pure moment). Corner region `r ≤ c/(2η)`, KKSDZ uniform bound K1′ `|dDisc| ≤ Br^{1−J₀}`:
`∫dη κ(η)∫₀^{c/2η}r^{a+1−J₀}dr ~ ∫dη κ(η)η^{−(a+2−J₀)}`.

| functional | `κ(η)` | corner converges under K1′ iff | `J₀ ≤ 1` (non-pert.) | check against Result D |
|---|---|---|---|---|
| raw `Ĉ₂` | `η^{3/2}` | `J₀ < −½` | **fails** | D: "shifted building blocks not absolutely integrable; swap with the total de-shifted kernel" ✓ |
| `Ĉ₃` (no de-shift needed) | `η^{−½}log η` | `J₀ < 3/2` | ok | D uses raw `Ĉ₃` ✓ |
| de-shifted `Ĉ′_EE` | `η^{−½}log η` | `J₀ < 3/2` | ok | D's Prop. K / K1′ corollary ✓ |
| **raw `Ĉ₄(0)`, `∂_{ν²}Ĉ₄(0)`** | **`η^{5/2}`** (×log) | **`J₀ < ½`** | **fails** | new |
| raw `Ĉ₅(0)` | `η^{3/2}` | `J₀ < 3/2` | ok (margin ½) | new |
| de-shifted (pure-moment) `Ĉ₄` | `η^{−½}log η` | `J₀ < 7/2` | ok (margin 2½) | needs the 4-shift de-shift [I] |

For the **tree correlators** of Lemma A the uniform bound is not needed: the exchange diagrams have an explicit impact-parameter profile
(`T` exchange: `∝ Π_{Δ−1}(β) ~ η^{−2}`), giving `∫dη η^{5/2}η^{−2}η^{−(4−J₀)} = ∫η^{−3/2}` for Einstein (`J₀ = 2`) — convergent
[D-sketch]. The `E×W³` and contact terms are assumed to have profiles at least as fast-decaying [D-sketch].

### 1.6 The MFT test (Task 1a)
* **Direct (commutator form) [P + DV].** For MFT, `⟨4[3,1]2⟩ = ⟨T₄T₂⟩⟨[T₃,T₁]⟩`. With the light-ray polarisations forced by `L[T₁], L[T₃]`
  (`z₁ ∝ x₁₀`, `z₃ ∝ x₃₀`) the numerator `(z₃·z₁)x₁₃² − 2(z₃·x₁₃)(z₁·x₁₃)` of `⟨T₃T₁⟩` **vanishes identically** for `x₁₃ = βz₃ − αz₁`,
  and equals `ε²(2(e·z₁)(e·z₃) + z₁·z₃)` with the `iε` shift [P, sympy]. So `⟨[T₃,T₁]⟩ = O(ε⁴)/(x₁₃²)⁵`. Paired with a smooth test function
  and any finite-order pole (`f₄` has at most `1/x₁₃²`), it tends to **0**: all moments `∫sⁿ(s+i)^{−p}ds`, `n ≤ p−2`, vanish.
  * Numerically, `|ε⁴∫φ[(t+iε)^{−p} − (t−iε)^{−p}]| = 3·10⁻⁵ → 2·10⁻¹¹` (`p = 5`, `ε = 0.1 → 0.003`), against an `O(1)` control without
    the `ε⁴`. With `f₄`'s pole (`p = 6`): `2·10⁻⁴ → 2·10⁻¹⁰`.
  * **Residual: 0 exactly, for any `f`, every `ν`** (`Ψ₄`, `∂_{ν²}Ψ₄`, `Ψ₅`).
  * This test is genuine but **weak**: it does not probe the `1/(u′v′)` poles (the two orderings' `W_s`, `W_t` terms cancel pointwise
    and identically).
* **Swapped (s-channel) form [I, not computed].** By Lemma C, the double commutator of `f·𝒢_MFT` collapses onto the edge `u′ = 0`:
  * the `W_u = ⟨42⟩⟨31⟩` term is the zero distribution above;
  * the identity gives `δ‴(x₁₂²)δ‴(x₃₄²)`-type terms;
  * `W_t = ⟨41⟩⟨32⟩` gives `δ(x₁₂²)δ(x₃₄²)`, which is exactly the sum over the leading `[TT]_{0,J}` family. The heavy `[TT]_{n,J}` are at
    exact double zeros (zero contribution).

  So the swapped MFT check is "identity edge term + `[TT]₀` edge term = 0": a one-dimensional integral along `u′ = 0` of the Taylor
  coefficients (to third order in `x₁₂²`, `x₃₄²` separately) of the spinning `EE`/celestial kernel. No existing engine evaluates the
  kernel off the heavy-block saddle (all engines are crystal/thimble evaluators of heavy blocks), so this was **not computed**. Residual
  and convergence: **not available**.

### 1.7 Verdict on V₄
* **V₄ holds for the tree correlators of Lemma A, grade [D]:**
  * edges, endpoints, Regge corner at fixed angle and `EE`'s `x₀`-edges: [D] via Lemmas B and C (B, C themselves [P]);
  * large-`η` Regge corner: [D-sketch] (explicit impact-parameter decay of the tree exchanges).
* **No counterexample:**
  * the `1/(u′v′)` poles carry the Wightman `iε` consistently;
  * they generate only the light-family edge terms that the light side already contains;
  * the MFT direct evaluation is exactly 0.
* **Gate: GO stays conditional, not a no-go.**
* **New caveat for the non-perturbative sum rule (I14(a)/(b) analogue, not V₄ proper):**
  * for the **raw** forward `Ĉ₄(0)` and `∂_{ν²}Ĉ₄(0)`, the KKSDZ uniform Regge bound does **not** give absolute convergence in the
    large-`η` corner (it needs `J₀ < ½`; KKSDZ give `J₀ ≤ 1`), exactly as for raw `Ĉ₂` in Result D;
  * raw `Ĉ₅(0)` is fine (`J₀ < 3/2`).
  * The pure-moment (de-shifted) `Ĉ₄` of Stage 12 §1.4 is therefore **load-bearing for validity**, not only a convenience: with it the
    corner needs `J₀ < 7/2`. Its construction (4-shift Abel problem) is [I].
  * Equivalently: the raw ν = 0 functional needs a Regge bound with impact-parameter decay, which is not a theorem.
  * The same power counting applies to CLS-D's raw `k = 4` functionals at any `ν` (`|P_{5/2+iν}| ~ η^{5/2}`).

---

## 2. Task 2: heavy-side positivity near the gap (`stage3_fwexact.py`, `stage3_neargap.py`, `stage3_smeared_check.py`)

### 2.1 Model (hybrid-exact skeleton) [M]
Per heavy state `(Δ, J)` with `Δ − J ≥ Δ_g`, `m² = (Δ−J−2)(Δ+J−1)`, `Λ = √((Δ_g−2)(Δ_g−1))`, `u = Λ²/m²`, `b̂ = (2J+1)√u`, per helicity
channel (`++`: even `J`; `+−`: `J ≥ 4`; the parity-odd even-`J` block couples like `++` at this order):
* **Smeared `k = 2, 3` parts** (Result D's hybrid):
  * Stage-0 flat parts with exact Wigner functions;
  * plus the leading-Regge curvature `∫ψ₂(2−up²)[P_{−½+iΛp}(cosh β) − J₀(b̂p)]` and `σu∫ψ₃[…]`.
* **Forward parts at `ν = 0`:**
  * flat graviton normal forms (`2u²`, `(c₄−2𝒥)u³`, `σu³`);
  * plus the **exact scalar-skeleton curvature** from the CHMRS Mellin representation (`exactmom.py`, `Δ_φ = 3`, no `1/m` expansion):
    * `B₄|₀ → u²ex₄(0)`;
    * `∂B₄|₀ → (c₄−2𝒥)u³ + u²[Λ²∂_{ν²}ex₄|₀ − u(5−2𝒥)]`;
    * `B₅|₀ → σu³ex₅(0)/2` (odd `k` by `Γ`-ratio continuation: a definition, [M]).
  * `ν = 0` is a removable `0/0` of the residue formula. It is reached by a cubic fit in `ν²` through `ν = 0.1 … 0.4`; an independent
    fit through `0.15 … 0.45` agrees to `≤ 2·10⁻¹¹` (value) and `≤ 10⁻⁹` (slope). `N₀ = 200` vs 400: `≤ 10⁻¹¹`.
* **Heavy action** `FF = B₂ + B₃ − a(∂B₄ − B₅) − cB₄` (`c = ½`), normalised by Stage 0's
  `N = |B₂| + |B₃| + |a(∂B₄−B₅)| + c|B₄|`.
* **Variants**, reported side by side:
  * `flat`: the Stage-0 action at the AdS state's `u`;
  * `LR`: Stage 12's leading-Regge forward pieces;
  * `mult`: the `∂B₄` curvature applied multiplicatively.
* **States**: `J ∈ {0, 2, 4–14, 16, 18, 20, 22, 24}`, `m/Λ ∈ {gap, 1.03, 1.06, 1.1, 1.15, 1.2, 1.3, 1.5, 1.75, 2, 2.5, 3, 4}`. At
  `Δ_g = 300, 1000`, `J ≥ 16` is kept only for `m ≥ 1.5Λ` (near the gap those states have `b̂ ≥ 25`, inside Stage 12's leading-Regge
  coverage), and `J = 18, 22` are dropped there.
  * Counts: 189 / 181 / 187 states, i.e. 298 / 274 / 283 (state, channel) cases at `Δ_g = 100 / 300 / 1000`. Exact-moment precision: the two `ν²`-fits agree to `≤ 1.4·10⁻⁸ / 1.7·10⁻⁷ / 2.6·10⁻⁶` (value) and `≤ 10⁻⁸` (slope); the largest cancellation ratio is `1.2·10⁴`.

**What the model does not contain [I]:** the graviton-specific (beyond-scalar-skeleton) corrections. For `Ω^EE` these were measured at
`0.1–0.5 N_s` near the gap (Result D, `K₂`/direct rows). Here they cannot be computed: at `ν = 0` the crystal saddle degenerates into the
Regge limit, and every existing graviton engine is a crystal/thimble evaluator at `ν ≳ 4.5`. Also missing: the raw-vs-de-shifted
difference of `Ĉ₄` (§1.5) and the off-diagonal helicity mixing beyond leading Regge.

### 2.2 Results
`cache/stage3_neargap_q0.txt` (Stage-0 functionals), `cache/stage3_neargap_q3.txt` (re-optimised margins). Minimum of `FF/N`:

| functional (flat `A`, `B`) | `Δ_g` | flat at AdS states | **exact skeleton** (location) | mult | LR | `b̂ < 10` / `10–20` / `≥ 20` (exact) | `c_max` |
|---|---|---|---|---|---|---|---|
| `S_min+B₃` robust, `δ = .05` (11.45, −12.98) | 100 | +0.0500 | **+0.0478** (`+−` `J=5`, `1.2Λ`, `b̂ = 9.2`) | +0.0495 | +0.0561 | +0.048 / +0.073 / +0.835 | 0.562 |
| | 300 | +0.0500 | **+0.0499** (`+−` `J=5`, `1.2Λ`) | +0.0501 | +0.0567 | +0.050 / +0.075 / +0.506 | 0.559 |
| | 1000 | +0.0500 | **+0.0500** (`++` `J=0`, gap) | +0.0500 | +0.0568 | +0.050 / +0.074 / +0.326 | 0.558 |
| `S_min` robust (15.00, −15.97) | 100 / 300 / 1000 | +0.051 | **+0.0460 / +0.0504 / +0.0510** (`+−` `J=5`, `1.2Λ`) | | | | 0.60 / 0.61 / 0.59 |
| `r = 5`, `ψ₃(0) = 0` (12.59, −14.66) | 100 / 300 / 1000 | +0.050 | **+0.0474 / +0.0499 / +0.0500** (`+−` `J=5`, `1.15Λ` / `++` `J=2`, `1.06Λ`) | | | | 0.56 |
| re-optimised `δ = 0.1` (13.02, −14.76) | 100 / 300 / 1000 | +0.0999 | **+0.0976 / +0.0997 / +0.0999** (`+−` `J=5`, `1.2Λ`) | | | | 0.63 |
| re-optimised `δ = 0.2` (17.24, −19.28) | 100 / 300 / 1000 | +0.2000 | **+0.1985 / +0.2003 / +0.2001** | | | | 0.81 / 0.80 / 0.80 |

* **No computed image is negative** in any variant, at any `Δ_g` [M]. `c_max` is the largest `B₄|₀` coefficient keeping every computed
  image `≥ 0`. It exceeds the normalisation `c = ½` everywhere, so the bound is not weakened.
* **Size of the AdS effects.**
  * The forward exact skeleton equals leading Regge to `|ex₄(0)/(2P₀) − 1| ≤ 1.0·10⁻³ / 1.2·10⁻⁴ / 1.0·10⁻⁵`. The `∂_{ν²}` piece differs
    from `−2Λ²D + 5uP₀` by up to 10 %, the finite-`J` constants that the LR model treats crudely.
  * So near the gap the curvature of the forward pieces is small and well described. The exact-minus-flat shift of the minimum is
    `−0.002 / −0.0001 / 0`.
* **Smeared parts** (`stage3_smeared_check.txt`, `Δ_g = 100`, `Q0-full-robust`). The exact `k = 2` skeleton remainder beyond the hybrid
  changes `FF/N` at the binding states by `+0.003` (`+−` `J=5`, `1.2Λ`), `−0.0035` (`J=5`, `1.1Λ`), `−0.005` (`++` `J=0`, gap),
  `−0.003` (`+−` `J=8`, `1.1Λ`), and `+0.19` (`+−` `J=4` at its gap, from the `x → 1` region). The `k = 3` remainder was not computed
  (odd-`k` continuation [M]).
* **Off-grid in the flat limit:** `dense_check.py` on the `δ = 0.1, 0.2` functionals passes (min `FF/N = 0.079`, `0.185`, no negatives;
  scaling limit to `b̂ = 2000` ≥ 0) [NE].
* **Verdict on positivity.**
  * In the model, every heavy state is positive: near-gap and small-`b̂` states by this computation, `b̂ ≥ 20` by Stage 12. The binding
    states are the flat ones (`+−`, `J = 5`, `m ≈ 1.2Λ`; `++`, `J = 0` at the gap).
  * The `δ = 0.05` functionals pass with essentially their flat margin, but that margin (`0.05 N`) is **not robust** against the
    uncomputed graviton-specific corrections at `J = 4–8`, `m ≈ 1.05–1.5Λ`.
  * The `δ = 0.2` functional tolerates corrections up to `≈ 0.2 N` at a 51 % cost in `A`.
  * Which one is legitimate depends on the graviton correction at `ν = 0`, which is [I].

### 2.3 Lower bound and the AdS constants
**Lower bound.** The functional `c_L B₄|₀ − (∂B₄|₀ − B₅|₀)` is heavy-positive iff `c_L ≥ max_states (∂B₄ − B₅)/B₄`. With the exact light
sides of Stage 12 (`−Ĉ₄(0) = 2g₄ + b₀G|ĝ₃|²`, `−∂Ĉ₄ = 4πG|ĝ₃|² + g₅`, `−Ĉ₅ = g₅`) this gives
`g₄Λ⁶/(8πG) ≥ [1/(4c_L) − b₀/(16π)]|ĝ₃|²Λ⁸`.

| `Δ_g` | `c_L`: flat at AdS states | **exact skeleton** | LR | slope `1/(4c_L)` (exact) | binding state |
|---|---|---|---|---|---|
| 100 | 0.9252 | **0.9848** | 0.708 | **0.2539** | `+−`, `J = 4`, gap (`m = 1.0396Λ`) |
| 300 | 0.9739 | **0.9813** | 0.732 | **0.2548** | `+−`, `J = 4`, gap |
| 1000 | 0.9921 | **0.9925** | 0.744 | **0.2519** | `+−`, `J = 4`, gap |

* Flat: `c_L = 1` at `u = 1`, which is unreachable in AdS because `u_gap(J=4) = (Δ_g−1)/(Δ_g+7)`.
* **The `¼` lower bound survives [M]** (slightly stronger), given `b₀` [I]. Graviton-specific corrections at the `J = 4` gap state are [I].

**AdS constants of the upper bound** (the functional is unchanged because it is positive; `A` is exact; `B_eff` uses the exact pole
`2L_AdS = 2∫ψ₂/(p² + 9/(4Λ²))`):

| functional | `Δ_g` (`Λ`) | `A` | `B_eff = 2L_AdS − A log(2Λ/3)` | `κ` | `2L_AdS` (bound at `ĝ₃ = 0`) | ABJM `2ζ(3)/2L_AdS` |
|---|---|---|---|---|---|---|
| `S_min+B₃`, `δ = 0.05` | 100 (98.50) | 11.448 | −12.804 | 0.0686 | **35.10** | 6.8 % |
| | 300 (298.50) | 11.448 | −12.938 | 0.0686 | **47.66** | 5.0 % |
| | 1000 (998.50) | 11.448 | −12.970 | 0.0686 | **61.45** | 3.9 % |
| `δ = 0.2` (robust variant) | 100 / 300 / 1000 | 17.241 | −19.038 / −19.220 / −19.265 | 0.1099 | **53.11 / 72.04 / 92.81** | 4.5 / 3.3 / 2.6 % |
| `δ = 0.1` | 100 / 300 / 1000 | 13.021 | −14.557 / −14.708 / −14.745 | 0.0802 | 39.93 / 54.22 / 69.90 | 6.0 / 4.4 / 3.4 % |

* `κ`: Stage-0 value. The `∂Ĉ₄` coefficient `4πG|ĝ₃|²` is exact [P given (F), (P)].
* The smeared parts' `|ĝ₃|²` coefficients carry `O(Λ⁻²)` AdS corrections (by analogy with Result D's `δk(0) ~ 10⁻²…10⁻⁴`) [M], plus the
  unknown `b₀` with leverage `≤ 2·10⁻⁴…1.5·10⁻²` per unit (Stage 12 §2.4) [I].
* Two-sided statement at the computed `Δ_g` [M + I]:
  `0.254|ĝ₃|²Λ⁸ ≲ g₄Λ⁶/(8πG) ≤ 11.448 log(2Λ/3) + B_eff − 0.0686|ĝ₃|²Λ⁸`.

### 2.4 ABJM
* ABJM in the IIA regime (Stage 12 §3, [DV given Q]): `g₄Λ⁶/(8πG) = 2ζ(3) = 2.404` with `Λ = 4.2156λ^{1/4}`, and `ĝ₃ = 0` (maximal SUSY).
* **Placement:** at 6.8 % / 5.0 % / 3.9 % of the AdS upper bound (robust `δ = 0.05` functional) at `Δ_g = 100 / 300 / 1000`
  (`λ ≈ 3.0·10⁵ / 2.5·10⁷ / 3.1·10⁹`). At 4.5 / 3.3 / 2.6 % of the `δ = 0.2` bound.
* It sits above the lower bound (`0 ≤ 2ζ(3)`, trivially, since `ĝ₃ = 0`).
* **Inside, far from both boundaries.** The test checks the sign and the `Λ⁻⁶` scaling only. The margin is log-dominated, so no
  near-saturation is possible for any theory with `Λ`-independent `g₄Λ⁶/G` (Stage 0 §5).
* A1 excludes ABJM's light KK spectrum, so the placement is a consistency check of the flat image, not an application of a theorem [D].

---

## 3. What remains open [I]
1. **Graviton-specific near-gap corrections of the forward pieces at `ν = 0`** (beyond the scalar skeleton), at the binding states
   `J = 4–8`, `m ≈ 1.05–1.5Λ`, `b̂ ≈ 5–10`. At `ν = 0` the crystal saddle degenerates into the Regge limit, and no existing graviton engine
   applies. This decides whether `δ = 0.05` (A = 11.45) or a larger margin (`δ = 0.2`: A = 17.24) is legitimate.
2. **The pure-moment (de-shifted) `Ĉ₄`** (4-shift Abel problem). It is now needed for the non-perturbative validity/swap of the forward
   `k = 4` sum rule (§1.5), not only for convenience. Alternative: a Regge bound with impact-parameter decay.
3. **Swapped-form MFT check** (identity + `[TT]₀` edge terms; §1.6), and an explicit Einstein `⟨TTTT⟩` evaluation of `Ψ₄`.
4. The `k = 3` skeleton remainder; the odd-`k` (`Ĉ₅`) skeleton is a definition by continuation [M]; parity-odd even-`J` block beyond
   leading order (as Result D's `H_OO`).
5. `b₀` (`W³×W³` constant in `Ĉ₄(0)`), `b₅`, Result D's `c₁`; the spinning (F) and (P) of Stage 12 (identification `g̃₄ = g₄`, ABJM).
6. Heavy states not in the grid:
   * near gap, `J ≥ 26`: `b̂ ≥ 50`, covered by Stage 12's leading-Regge model [M];
   * `m > 4Λ`: flat regime, covered by the Stage-0 dense check [NE];
   * the `H_LJ`-type beyond-leading-Regge corrections at large `J` near the gap, as in Result D.
7. Everything inherited: A1, A2, KKSDZ loophole, K1′, Stage-12 (i) (microcausality of tree correlators), `J₀(E×W³) ≤ 3`.

## 4. Graded claims
| claim | grade |
|---|---|
| `f⁺₄ = f⁺₂R`, `f⁻₅ = f⁻₃R`, `R = ((v′−u′)²−1)/(4u′v′)`; frame form; all boundary limits | [P] (sympy) |
| `R` bounded at single-ray endpoints | [DV] |
| same `iε` for `f`'s poles in both orderings; `f` analytic off the real axis on light rays | [P] |
| edge: no boundary terms beyond the light-family terms (analytic continuation = `iε` distribution) | [D] |
| `EE` `x₀`-edges and forward `ν = 0` (no `H_edge` window) inherited / absent | [D] / [P] |
| **V₄ for Lemma A's tree correlators** | **[D]** (large-`η` corner [D-sketch]) |
| MFT direct (commutator form) value = 0 for any `f`, all `ν` | [P] + [DV] |
| MFT swapped (s-channel) form | [I] (not computed) |
| raw forward `Ĉ₄` not absolutely convergent in the large-`η` corner under the KKSDZ uniform bound (needs `J₀ < ½`); `Ĉ₅` ok (`J₀ < 3/2`); de-shifted `Ĉ₄` ok (`J₀ < 7/2`) | [D] (power counting on Stage 12's leading-Regge kernels; reproduces Result D's raw-`Ĉ₂` / `Ĉ₃` verdicts) |
| exact skeleton forward moments `ex₄(0)`, `∂_{ν²}ex₄(0)`, `ex₅(0)` on 557 states | [DV] (moments; odd `k` a definition [M]) |
| near-gap positivity of the Stage-0 functionals at `Δ_g = 100, 300, 1000` (min `+0.048 / +0.050 / +0.050`) | [M] (hybrid-exact skeleton; graviton-specific part [I]) |
| `δ = 0.2` re-optimisation: `A = 17.24`, flat dense check, AdS min `+0.1985` | [NE] (flat) + [M] (AdS) |
| lower-bound slope `≥ 0.2519` (¼ survives) | [M], `b₀` [I] |
| AdS `A` exact, `B_eff`, `κ` | `A` [P given the pole]; `B_eff` [DV]; `κ` [M] + `b₀` [I] |
| ABJM at 3.9–6.8 % (robust), 2.6–4.5 % (`δ = 0.2`) | [D given Q]; bound [M] |

## 5. Files (`ads4/code/resultq/`)

| file | content | run time |
|---|---|---|
| `stage3_v4.py` | Lemma B (`R`, frame form, all boundary limits), endpoint numerics, `iε` ordering table, light-ray-polarised `TT` identity, `ε⁴` pairing → `cache/stage3_v4.txt` | seconds |
| `stage3_fwexact.py` | exact scalar-skeleton `ex₄(0)`, `∂_{ν²}ex₄(0)`, `ex₅(0)` on the near-gap state grid (`Δ_φ = 3`, `N₀ = 200`, two independent `ν²`-fits) → `cache/stage3_fw_<Δg>_Df3_w*.pkl`, logs `cache/stage3_fw_w*.log` | ≈ 2 h on 2 cores (`stage3_run.sh`) |
| `stage3_neargap.py` | heavy images (flat / exact / mult / LR), minima by region, lower-bound `c_L`, `c_max`, `2L_AdS`, `B_eff` → `cache/stage3_neargap.txt`, `.pkl`; env `S3_FUNC=functionals_q3.json` for the re-optimised margins | ≈ 1 min per functional |
| `stage3_smeared_check.py` | exact `k = 2` skeleton remainder vs the hybrid on the binding states (`Δ_g = 100`) → `cache/stage3_smeared_check.txt` | minutes |
| `stage3_run.sh`, `stage3_after.sh` | queues | — |
| `cache/functionals_q3.json` | LP functionals with margins `δ = 0.1, 0.2, 0.3` (`lp_g4.py 12 12 3 --delta δ`) + `Q0-full-robust` | — |
