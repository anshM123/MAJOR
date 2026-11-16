# `Ω′`: replacing `C⁻₃` by a `k = 3` bare light-ray commutator functional

Status: written 16 Nov 2026 by the `LIGHT_INPUTS.md` agent. Updated the same day:
* §6: the normalisation `N(ν)` and the flip check;
* §7: the `(2,4)` side;
* **revised after `REFEREE_OMEGA_PRIME.md`** (M1–M3, m1–m5; list in §9); new §8 on reality and `ν`-evenness.

Code: `ads4/code/lightinputs/` `t6_normalisation.py`, `t6b_flip_scan.py`, `t7_shadow_reality.py` (logs in `cache/`, git-ignored). **Not refereed in this
revised form.** It builds on:
* `LIGHT_INPUTS.md` §3 (the bare commutator; (I1) is false for the project's `C⁻₃`);
* `ADS4_FUNCTIONAL.md` §1;
* `SPINNING_E2E.md` §3.3, §4.4 (frame modes, light-ray mode map, `M`-tables);
* `SPINNING_4PT.md` §4.2–4.6 (two-pair heavy formula);
* `FLAT_GRAVITON_BOUND.md` §1.

Marks:
* **[D]** derived; **[D-sketch]** outlined; **[DV]** derived and checked by exact sympy evaluation;
* **[Q]** quoted: KKSDZ = 1904.05905, CLS-D = 2311.04271, CLPSD = 2201.06602;
* **[I]** inferred; **[ND]** not done.

**Premises used throughout** (they are part of every [D]/[DV] label below that concerns the functional's heavy action):
* **(P1) light-ray mode map.** At leading Regge order the `(1,3)` part of a bare-commutator functional is the antipodal mode combination
  `w_q = (−1)^{q₃}d_{(q₁+q₃−1)/2}(δ)/d₀(δ)` of `SPINNING_E2E.md` §3.3 [Q+D].
* **(P2) `k_r = 3`.** In the CHMRS normalisation (kernel `r^{k−2}`) the bare `T` commutator has leading moment index `3`. This is [D-sketch]: it is fixed by the
  convergence threshold `J₀ < J₁+J₃−1 = 3` (`LIGHT_INPUTS.md` §2.3).
* **(P3) two-pair heavy formula.** `SPINNING_4PT.md` §4.2: established for currents, conditional for stress tensors (`SPINNING_TTTT.md`).

## 0. Results

1. **Definition.** `C̃⁻_{3,ν} := Ψ̃_{3,ν}/N(ν)`, which is **real and even in `ν`** (item 4).
   * `Ψ̃_{3,ν}` is the bare light-ray commutator functional `∫K·⟨T[L[T](x₀,z₃), L[T](x₀,z₁)]T⟩`, normalised so that its leading `(1,3)` kernel is the `d₀`-normalised
     antipodal combination (P1).
   * Its kernel is the CLS-D parity-odd celestial structure (`j′ = 1`) at `(1,3)`.
   * At `(2,4)` it is the unique covariant three-point lift of the project's odd kernel `𝔄 = D₅D₁` (§7).
   * It is **not** shown to be a member of CLS-D's family. Their duals are light transforms of `⟨TTO_{3/2+iν,3}⟩` structures, at most 23 of them. The `𝔄` lift lives in
     the 25-dimensional space of covariant kernels for `L[O]` of continuous spin `−½−iν` (referee M1).
   * It is a valid superconvergence sum rule by itself, because the commutator vanishes as an operator (§7).
2. **Heavy action, leading order (`ν ≪ m`) [DV given P1–P3].**
   * Form: `C̃⁻₃[Δ,J] = (2P_{−½+iν}(cosh β)/m⁶)(|c_same|² − |c_opp|²)(1 + O(1/m²))`. It is *identical* to the leading heavy action of the project's `C⁻₃`, with the same
     helicity combination and the same pres–pres harmonic.
   * The normalisation `N(ν) = −(32/7)(2ν+i)(4ν²+9)/(4ν²+201)` is exact (§6). It is nonzero on real `ν`, in the `d₀` convention. The convention gap is in item 5.
   * **No helicity-flip interference [D, all `ν`, `a`, `l`; DV for `l = 1, 2`].** The flip responses vanish at both pairs **by symmetry** (referee's argument, §6).
   * This does not make the action sign-definite. `|c_same|² − |c_opp|²` is indefinite, exactly as in `Ω`. Positivity comes only from the full functional (Theorem 1).
3. **Flat limit [D given P1–P3 and the `ADS4_FUNCTIONAL.md` §1.2 dictionary].**
   * The leading heavy action equals that of `C⁻₃`, so the flat limit is the same: `ε₃Λ⁶C̃⁻_{3,Λp} → B₃^{(1)}(p)`, with `σ_h = ±1`.
   * Consequently `Ω′` has the flat limit `F_{g3}`. **Theorem 1 of `FLAT_GRAVITON_BOUND.md` applies**, and the leading-order AdS heavy analysis of
     `ADS4_FUNCTIONAL.md` carries over verbatim.
   * Consistency check: in `D = 4` the odd fixed-`u` sum rules are only the `B^{(1)}_k` (CLPSD [Q]).
   * CLS-D's (200) is **no longer cited** as a premise (M1).
4. **Light side and reality [D given planar `J₀ < 3`, HM/CDP, (F); evenness DV at leading-Regge kernel level (`t7`), D-sketch for the full functional].**
   * `−C̃⁻₃|_light = λ²_TTT Ψ̃₃[G_T]/N(ν) ∝ t₄²+α_T² = 57600|ĝ₃|²`, with **exactly** no Einstein, linear-`t₄`/`α_T`, contact or `[TT]` term at `O(1/c_T)`.
   * **`C̃⁻_{3,ν}` is real and even in `ν`.** This follows from two exact identities (§8), which together give `C̃_{−ν} = C̃_ν = C̃̄_ν`:
     * a shadow identity, `Ψ̃_{−ν} = [δ/(1−δ)]Ψ̃_ν`, which `N` shares exactly;
     * a reality identity, `Ψ̃_{−ν}[G] = −Ψ̃̄_ν[G]` and `N(−ν) = −N̄(ν)`.
   * Hence its light side is even in `ν`, so no `ν³` or `ν` term can appear (M3). Then `−C̃⁻₃|_light = −4πG|ĝ₃|²(ν⁴ + f̃₁ν² + f̃₀) + G|ĝ₃|²h(ν)`, with `h` even and decaying.
   * So `k(0) = 330151/12012` and `A = 37.8387` are unchanged. The `k = 3` term shifts `k_AdS` only at `O(Λ⁻²)` (through `f̃₁`, [ND]) and at `O(Λ⁻⁵)` (through `h`).
   * Hence the claim "`B_AdS`, `C₁` are unchanged by the `k = 3` replacement" holds, given the full-functional step of §8.
   * `Ω′` uses `C̃⁻₃ = Re C̃⁻₃` (m1).
5. **Remaining gap (M2): the light-ray Jacobian `j(ν)`.**
   * `C̃` is defined through its `d₀`-normalised leading kernel. Here `d₀(δ) = −i2^δΓ(1−δ)/(Γ((1−δ)/2)Γ((3−δ)/2))` is finite and nonzero on ℝ [D, DV `t7`].
   * The constant `j(ν)` relating this kernel to the leading moment of the light-ray integral itself is **not tracked** (`SPINNING_E2E.md` §3.3). It comes from
     the change of variables (light-ray parameters, `x₀`, `z₀`) → (`r`, `η`, frame angles) and the `z₀`-integral.
   * `C̃` is a rescaled sum rule wherever `j(ν)` is finite and nonzero. That is plausible (Γ-ratios), but **not shown**.
   * A real zero of `j` at isolated `ν` would make the subleading Regge moments of `C̃` singular there. It would not affect the leading heavy action.
   * The zeros and phase of `N` quoted below refer to the `d₀` convention. The former "`T`-pole" remark is convention-dependent and has been withdrawn.
6. **Effect on the inputs.** For the `k = 3` term of `Ω′`, input (E) needs only the planar chaos bound and HM/CDP data. Still load-bearing:
   * the forward terms (RB₃/RB₄, [I]);
   * (R) for `C⁺₂` ([D-sketch given IP]);
   * `j(ν) ≠ 0, ∞` (item 5, [I]);
   * P2 ([D-sketch]) and P3 (conditional).

## 1. Sum rule and double traces

* **Sum rule.** `Ψ̃₃[G] = 0` for every planar correlator with `J₀ < 3`. This covers planar holographic CFTs, by the chaos bound [Q KKSDZ]. The commutator vanishes as an operator,
  so any kernel at `(2,4)` gives a valid sum rule (§7).
* **Double traces.** `[TT]` families enter at `O(γ²)`, because `⟨0|[T,T]|[TT]⟩ = O(γ)` (`LIGHT_INPUTS.md` §3.2).
* **Parity.** The `(1,3)` celestial structure is parity-odd. A parity-even `(2,4)` kernel would give zero leading action on parity-even exchanges, so the `(2,4)` kernel is
  chosen odd (`𝔄`).

## 2. Heavy action at leading order

* **Formula.** By P3, `Π̂^q → 𝒩γ_a²m^{−2a}Σ_h c^{h₁h₂}c̄^{h₃h₄}M^{q₁q₃}_{h₁h₃}(ν)M^{q₂q₄}_{h₂h₄}(−ν)𝒦^ν_{μ_L(h),μ_R(h)}(β)`.
  * The `β`-harmonic is fixed by the helicity flow (`μ = 0` preserved, `±2l` flipped).
  * The kernel enters only through the responses `A(h₁,h₃) = Σ_q w_q M^q_{h₁h₃}`.
* **Structure.** Both the antipodal `(1,3)` kernel and `𝔄` are purely odd, `A(+,−) = −A(−,+)` (§6). So the two-pair product weighs `h₁h₂`, i.e. `+` for same and `−` for opposite.
  The physical sign follows the covariant convention of `SPINNING_4PT.md` §4.5.
* **Result.** The pres–pres action is `N(ν)` × that of `C⁻₃`. The MFT normalisation is a property of the states and is untouched. Flip pieces vanish (§6).

## 3. Flat limit

The flat limit is inherited from `C⁻₃` through the identity of the leading heavy actions (§2), and then from the `ADS4_FUNCTIONAL.md` §1.2 dictionary (label, power `m^{−2k}`,
sign `σ_h`, `ε₃ = ½`; [D|H*] there). Two consistency checks, neither needed for the conclusion:
* In `D = 4` the only odd fixed-`u` sum rules are CLPSD's `B^{(1)}_k` (2.26) [Q]. This is because `g` is fully crossing symmetric and `h(s,−p²) ~ s`.
* The light side matches flat `−B₃|_low = −2πG|ĝ₃|²p⁴`: no Einstein, no `g_i` (§4).

## 4. `Ω′` and its light side

**Definition.** `Ω′_{g3} := ∫₀¹dp [ψ₂Λ⁴C⁺_{2,Λp} + ½ψ₃Λ⁶ Re C̃⁻_{3,Λp}] − 5Λ¹⁰∂_{ν²}C⁺_{4,ν}|₀ + (5/2)Λ¹⁰C⁻_{5,0}`. By §8, `Re C̃⁻₃ = C̃⁻₃`.

**Positivity.**
* At leading order `Ω′|_heavy = Ω|_heavy` for every heavy state [DV given P1–P3]. So Theorem 1 (flat) and the leading-Regge AdS analysis apply unchanged.
* The `O(1/m²)` terms of `C̃⁻₃` differ from those of `C⁻₃`. Neither set is computed.

**Light side of the `k = 3` term.**
* `λ²_TTT Ψ̃₃[G_T] ∝ ε†[H₁,H₃]ε = −2i|c|²sin4φ₁₃σ_z`, with `|c|² = (t₄²+α_T²)/16 = 3600|ĝ₃|²`.
* By (F) and evenness, `−C̃⁻₃|_light = −4πG|ĝ₃|²(ν⁴ + f̃₁ν² + f̃₀) + G|ĝ₃|²h(ν)`, with `h` even and decaying.

**Effect on the bound.**
* `k(0)` is unchanged.
* The `k = 3` contribution to `k_AdS − k(0)` is `225f̃₁/(14Λ²) + 7f̃₀/(3Λ⁴)` [ND coefficients].
* A decaying `h` enters at relative `O(Λ⁻⁵)`, since `½Λ⁶∫ψ₃h(Λp)dp ≈ 81Λ⁵∫h` against `k(0)πGΛ²x²` with `x² = |ĝ₃|²Λ⁸`.
* No `O(Λ⁻¹)` shift arises from this term (§8). `A`, `B_AdS`, `C₁` are as in `ADS4_LIGHT.md` §4, which is still conditional on (R) and on the forward terms.

## 5. Status table

| item | status |
|---|---|
| `Ψ̃₃` is a valid sum rule (bare commutator, any `(2,4)` kernel), planar `J₀ < 3` | [D] given [Q] KKSDZ; operator statement (§7) |
| leading heavy action `= N(ν)·(2P/m⁶)(|c_s|²−|c_o|²)`, `N = −(32/7)(2ν+i)(4ν²+9)/(4ν²+201)` | [DV given P1–P3] |
| no flip interference at either pair | [D] by symmetry, all `ν, a, l`; [DV] `l = 1, 2` |
| flat limit `= ψ₃B₃^{(1)}`, Theorem 1 applies | [D given P1–P3 + `ADS4_FUNCTIONAL.md` dictionary] |
| `C̃⁻₃` real and even in `ν`; no odd-`ν` light terms | [DV] at leading-kernel level (`t7`); [D-sketch] full functional (§8) |
| light side: no Einstein/linear/contact/`[TT]` term; top `−4πG|ĝ₃|²ν⁴` | [D given planar `J₀ < 3`, HM/CDP, (F)] |
| `d₀` finite and nonzero on ℝ | [D, DV] |
| light-ray Jacobian `j(ν)` finite and nonzero on ℝ | **[I], not tracked** |
| `f̃₁, f̃₀`; `O(1/m²)` heavy terms | [ND] |

## 6. `N(ν)` and the flip check (`t6_normalisation.py`, `t6b_flip_scan.py`)

**Method.** This uses the exact `M`-tables of `ads4/code/e2e/laplace.py` (imported from the existing cache), giving the response `M^{q₁q₃}_{h₁h₃} = R + g₁S`
(`g₁ = γ_{a−1}/γ_a`) of each frame mode to a heavy transverse-helicity plane wave.
* `(1,3)` kernel: the weights `w_q` of P1. `t6` now builds them with exact integers; the imported `e2e/clsd_kernel.py` gives `±1.0` floats for `q₃ < 0` (m4), and the two
  agree.
* Project kernels: `D₅` (`l = 1`), `𝔄 = D₅D₁` (`l = 2`).
* All at `a = 2l + k + 1`, `k = 3`.

**Results (exact) [DV given P1–P3].** The referee reproduced them independently by direct `H²` quadrature, to `5·10⁻⁸`.

| | `l = 1` (`a = 6`) | `l = 2` (`a = 8`) |
|---|---|---|
| project `A(±,∓)` | `±(5i/8)g₁` (∝ the tabulated `D₅` entry) | `±(7i/64)g₁(4ν²+201)/(4ν²+81)` (∝ the tabulated `𝔄` entry) |
| bare commutator `A(±,∓)` | `±g₁(iν/2 − ¼)` | `±g₁(4ν²+9)(1−2iν)/(2(4ν²+81))` |
| pres even part / flip responses | 0 / 0 | 0 / 0 |
| ratio | `(2/5)(2ν + i)` | `N(ν) = −(32/7)(2ν+i)(4ν²+9)/(4ν²+201)` |

* The `l = 1` column is a **parallel computation for currents**, not a check of `N` for stress tensors (m3). It only checks the `D₅` table entry.
* **Zeros (in the `d₀` convention).** `|N|² ∝ (4ν²+1)(4ν²+9)²/(4ν²+201)²` is strictly positive on ℝ. `|N(0)| = 288/1407 ≈ 0.205`, and `|N|` grows like `(64/7)ν`.
* **Flip-freedom is structural [D] (referee).**
  * At `(1,3)` with `z₁,z₃ = ŷ ± X(φ)`, the flip tensor `(A² − (X·ε)²)^l` is even under `φ → φ+π`, while the weight `sgn sin(−φ)|sin φ|^{−δ}` is odd.
  * At `(2,4)`, `𝔄` is antisymmetric under `z₂ ↔ z₄` while the flip tensor is symmetric.
  * So `A(±,±) = 0` for every `ν`, `a` and `l`.
  * Individual frame modes do have flip responses: 20 of the 25 modes at `l = 2` (`t6b`).

## 7. The `(2,4)` side: any covariant kernel is allowed, in particular `𝔄` [D]

* **Why the kernel is free.** `⟨T(x₄,z₄)[L[T](x₀,z₃), L[T](x₀,z₁)]T(x₂,z₂)⟩ = 0` holds pointwise, being an operator statement (KKSDZ; planar `J₀ < 3`). So *any* kernel at `(2,4)`
  gives a valid sum rule. The referee checked that the two possible sources of contact terms are controlled: the shared endpoints are controlled by `J₀ < 3`, and
  collinear `z₁ ∥ z₃` is suppressed by `|z₁₃|^{2+δ}`.
* **Covariant kernels with spectral parameter `ν`.** These are three-point structures of `(T₂, T₄, X)` with `X` at `x₀`, polarisation `z₀ = ℓ`. At leading Regge order they are
  evaluated at `(x₂, x₀, x₄) = (−x̂, 0, x̂)`.
  * The stabiliser is the `SO(2)` of rotations about `x̂`, and it acts transitively on `z₀` mod scale.
  * So the structures are all `SO(2)`-covariant bidegree-`(2,2)` polynomials in the frame `(x̂, n_ℓ, n_ℓ^⊥)`: the 25 frame modes, 13 parity-even and 12 parity-odd.
  * The relevant count is the generic one, `(2l₁+1)(2l₂+1) = 25`, because `X ∼ L[O]` has *continuous* spin `−½−iν` (m2). The earlier analogy with `(2,2,J ≥ 4)` was inapposite:
    `O` itself has `J = 3`, with 23 structures.
  * The covariant extension of a given leading kernel is unique, because three points have no cross-ratio.
* **Asymmetry with `(1,3)`.** At `(1,3)` the light transforms force antipodal polarisations, so `𝔄` cannot be used there. At `(2,4)` there is no such constraint.
* **Choice.** The `(2,4)` kernel is the lift of `𝔄`. So the `(2,4)` factor of `N` is 1 and there is no flip at `(2,4)` (§6).

## 8. Reality and `ν`-evenness of `C̃⁻₃` (`t7_shadow_reality.py`) [DV at the leading-kernel level; D-sketch for the full functional]

At leading Regge order the frame kernel of `Ψ̃_ν` is `Σ_q w_q(ν)F^ν_{q₁+q₃}(ξ) ⊗ 𝔄`. Here `𝔄` is a real, `ν`-independent `(2,4)` kernel with `μ_R = 0`, so its harmonic is `F^ν_{μ_L}`.
Write `F^ν_μ = P^{|μ|}_{−½+iν}/(½+iν)_{|μ|}`. The associated Legendre function `P^μ_{−½+iν}` is even in `ν`, so `F^{−ν}_μ = F^ν_μ(δ)_{|μ|}/(1−δ)_{|μ|}`. Exact identities (`t7`, all odd modes, `l = 1, 2`):

1. **Shadow.** `w_q(−ν)(δ)_{|μ|}/(1−δ)_{|μ|} = [δ/(1−δ)]w_q(ν)`, the Pochhammer identity of the `sgn|sin|^{−δ}` Fourier coefficients. Hence `Ψ̃_{−ν} = [δ/(1−δ)]Ψ̃_ν`.
2. **Reality.** `w_q(−ν) = −w̄_{−q}(ν)` for real `ν`. On real correlators (`𝒢_{−q} = 𝒢̄_q`) this gives `Ψ̃_{−ν}[G] = −Ψ̃̄_ν[G]`.
3. **`N`.** `N(−ν)/N(ν) = δ/(1−δ)` and `N(−ν) = −N̄(ν)`, exactly.

Hence `C̃_{−ν} = C̃_ν` (by 1 and 3) and `C̃_{−ν} = C̃̄_ν` (by 2 and 3): **`C̃⁻₃` is real and even in `ν`**.
* Its heavy action is real and even, as §2 requires. That is a consistency check.
* Its light side `λ²Ψ̃[G_T]/N` is a real even function of `ν`. So the `ν³` term feared in M3 is absent, and `Re` in `Ω′` is redundant.
* **Full functional [D-sketch].**
  * Identity 2 lifts trivially, by complex conjugation of the full kernel.
  * Identity 1 lifts if the covariant extensions at `ν` and `−ν` are related by the shadow transform in `z₀`. They are: the celestial structure is unique for `j′ = 1`,
    the `(2,4)` lift is unique (§7), and shadow transforms preserve covariance.
  * Hence the full functionals satisfy 1 up to the same constant. This step was not written out in detail.

## 9. Revisions after `REFEREE_OMEGA_PRIME.md`

* **M1.** The functional is no longer called a CLS-D member. The flat limit is re-derived from the identity of the leading heavy actions plus the `ADS4_FUNCTIONAL.md` dictionary.
  CLS-D (200) is dropped as a premise. P1 and P2 are now part of the labels.
* **M2.** The untracked Jacobian `j(ν)` is stated as a gap. `d₀` is computed in closed form (finite, nonzero). The zeros and phase of `N` are now attributed to the convention, and the
  `T`-pole remark is withdrawn.
* **M3.** `C̃⁻₃` is shown real and even (§8). So no odd-`ν` light term arises, and `B_AdS`, `C₁` are unchanged, given the full-functional sketch. `Ω′` is written with `Re C̃`.
* **m2.** Counting corrected (continuous spin, 13 + 12).
* **m3.** The `l = 1` column is relabelled as a parallel computation.
* **m4.** The weights are now exact in `t6`.
* **m5.** "No sign-indefinite terms" is replaced by "no flip interference"; `|c_s|² − |c_o|²` is indefinite as in `Ω`.
* **Referee's stronger statement adopted.** Flip-freedom is structural, for all `ν`, `a`, `l`.
