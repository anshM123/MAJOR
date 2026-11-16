# `Ω′`: replacing `C⁻₃` by the `k = 3` bare light-ray commutator

Status: written 16 Nov 2026 by the `LIGHT_INPUTS.md` agent. Mostly analytic. **Updated the same day with the computed normalisation `N(ν)` and the
flip check (§6; code `ads4/code/lightinputs/t6_normalisation.py`, `t6b_flip_scan.py`).**
**Not refereed.** It builds on `LIGHT_INPUTS.md` §3 (bare commutator; (I1) false for the project's `C⁻₃`), `ADS4_FUNCTIONAL.md` §1,
`SPINNING_4PT.md` §4.3–4.6 and `FLAT_GRAVITON_BOUND.md` §1.

Marks: **[D]** derived; **[D-sketch]** outlined; **[Q]** quoted (CLS-D = 2311.04271; CLPSD = 2201.06602, LaTeX read; KKSDZ = 1904.05905);
**[DV]** derived and checked by exact (sympy) evaluation of the `laplace.py` tables; **[I]** inferred; **[ND]** not done.

## 0. Results

1. **Definition.** `C̃⁻_{3,ν} := Ψ̃_{3,ν}/N(ν)`.
   * `Ψ̃_{3,ν}` is a member of the CLS-D `k = 3`, `J_e = 2` family with constant `f⁻_{3;2}` (the *bare* commutator `⟨T[L[T],L[T]]T⟩`).
   * It has the parity-odd celestial structure at the `(1,3)` pair (`j′ = 1`) and a parity-odd dual structure at the `(2,4)` pair.
   * `N(ν)` is fixed so that the leading heavy action equals that of the project's `C⁻₃`.
2. **Heavy action, leading order (`ν ≪ m`) [DV at the `(1,3)` pair, §6; D for the `(2,4)` choice, §7].**
   * Form: `C̃⁻₃[Δ,J] = (2P_{−½+iν}(cosh β)/m⁶)(|c_same|² − |c_opp|²)(1 + O(1/m²))`. This is the **same** helicity combination and the **same** `ν`-kernel
     (the pres–pres harmonic `P = 𝒦_{0,0}`) as `C⁻₃`.
   * Normalisation [DV]: `N(ν) ∝ (2ν + i)(4ν² + 9)/(4ν² + 201)` at `a = 8` (`l = 2`, `k = 3`). It is **nonzero for every real `ν`** (`|N|` is smallest at `ν = 0`).
     Its zeros sit at `ν = −i/2, ±3i/2`, off the real axis. So `Ω′` is well defined on the whole real `ν` axis.
   * No helicity flip [DV]: at the `(1,3)` pair the bare-commutator flip responses `A(+,+)`, `A(−,−)` vanish identically. Individual frame modes do have odd flip
     responses (`t6b`); they cancel in the antipodal combination. The project's `𝔄` kernel also has zero flip response. So there is no
     `c_same c̄_opp` piece: the `(2,4)` kernel can be *chosen* equal to `𝔄` (§7, [D]), and `𝔄` has zero flip response [DV].
3. **Flat limit [D given Q].** In `D = 4` the odd-`k` fixed-`u` sum rules are exactly the `B^{(1)}_k` of CLPSD (2.26) (label `defB`) [Q]. CLS-D map every bare-commutator functional
   to a lowest-subtracted odd fixed-`u` flat sum rule [Q]. Hence `ε₃Λ⁶C̃⁻_{3,Λp} → B₃^{(1)}(p)`, with `σ_h = +1` for `++` and `−1` for `+−`, as before.
   * Consequence: `Ω′` has the same flat limit `F_{g3}` as `Ω`. **Theorem 1 of `FLAT_GRAVITON_BOUND.md` applies unchanged** to the flat leading heavy action.
   * The AdS leading-order heavy analysis of `ADS4_FUNCTIONAL.md` (thresholds ≈7.5Λ etc.) carries over **verbatim**, since the leading heavy action is identical.
   * What changes: the `O(1/m²)` heavy corrections of the `k = 3` term (not computed, for either functional).
4. **Light side [D given planar `J₀ < 3`].**
   * `−C̃⁻₃|_light = λ²_TTT Ψ̃₃[G_T]/N(ν) = (t₄²+α_T²)·X(ν) = 3600|ĝ₃|²X(ν)`.
   * It contains **exactly** no Einstein term, no linear `t₄`/`α_T` term, no contact term (`g₄, g₅, g₆, …`) and no `[TT]` family at `O(1/c_T)`.
   * By (F), `X(ν) = −(4πG/3600)(ν⁴ + f̃₁ν² + f̃₀) + h(ν)`, where `h` is any decaying remainder.
   * The top coefficient is the old one, so **`k(0) = 330151/12012` and `A = 37.8387` are unchanged** [D given (F)]. `f̃₁` replaces `f₁` in the `Λ⁻²` term of
     `k_AdS` [ND].
   * A non-polynomial `h` is now harmless: it enters `Ω` at relative order `Λ⁻⁵` (§4).
5. **Effect on the inputs.** In `Ω′`, input (E) for the `k = 3` term no longer depends on any Regge bound for the Einstein truncation or on (I1). It needs only:
   * the planar chaos bound [Q];
   * HM/CDP one-point data [Q].

   This removes the `~10⁴` leverage item of `ADS4_LIGHT.md` §2.1. Still load-bearing:
   * the forward terms `∂_{ν²}C⁺₄|₀` and `C⁻₅|₀` (RB₃/RB₄, [I]);
   * (R) for `C⁺₂` ([D-sketch given IP]);
   * the new normalisation `N(ν)`: computed, nonzero on real `ν` [DV, §6], and complete, because the `(2,4)` factor is 1 by choice (§7) [D].

## 1. Definition and existence

**Two-pair structure.**
* At the `(1,3)` pair, CLS-D's functional applies `L[T](x₀,z₁)` and `L[T](x₀,z₃)`, which in the CHMRS frame means antipodal polarisations (`LIGHT_INPUTS.md` §3.1).
* The celestial structure `⟨P̃₁P̃₃P_{δ,j′}⟩` with `δ = ½+iν` is `∝ |sin(θ₀−φ)|^{−δ}`, times `sgn sin(θ₀−φ)` for `j′ = 1` [Q, `SPINNING_E2E.md` §3.3].
* At the `(2,4)` pair the dual structure `(⟨T₄L[O]T₂⟩^{(a)})^{−1}` is free. It ranges over the `⟨TTO_{3/2+iν,3}⟩` structures; CLS-D count 2 odd functionals in
  general `d` for `j′ = 1`. In `d = 3` the celestial group is trivial, so the count collapses [I: the number of independent `d = 3` members was not determined]. Superseded by §7: any covariant `(2,4)` kernel is admissible, and we choose `𝔄`.

**Choice.** Choose the `(2,4)` structure parity-odd. A parity-even `(2,4)` structure combined with the odd `(1,3)` structure has zero leading action on
parity-even exchanges (the `h₁`-odd weight sums to zero over parity images; `SPINNING_E2E.md` §4.4: "`S` cancels for parity-even blocks").

**Sum rule and double traces.**
* The sum rule `Ψ̃₃[G] = 0` holds for every planar correlator with `J₀ < 3`, which includes all planar holographic CFTs by the chaos bound [Q KKSDZ].
* `[TT]` families enter at `O(γ²)`, because `⟨0|[T,T]|[TT]⟩ = O(γ)` (`LIGHT_INPUTS.md` §3.2).

## 2. Heavy action at leading order [D-sketch]

The machinery is the two-pair heavy formula (`SPINNING_4PT.md` §4.2), **established for currents and conditional for stress tensors** (`SPINNING_TTTT.md`):
`Π̂^q → 𝒩γ_a²m^{−2a}Σ_h c^{h₁h₂}c̄^{h₃h₄}M^{q₁q₃}_{h₁h₃}(ν)M^{q₂q₄}_{h₂h₄}(−ν)𝒦^ν_{μ_L(h),μ_R(h)}(β)`.
The `β`-harmonic is fixed by the state's helicity flow (`μ = 0` for preserved, `±2l` for flipped), not by the kernel mode. The kernel enters only
through the `M`-tables, i.e. through the rational functions `R` and `g₁S`.

* **Pres–pres part.**
  * The `(1,3)` kernel is parity-odd, so only the odd parts survive: `M^{(13)}_{h,−h} → h·S̃_L(ν)` for both pairs.
  * The product therefore weighs `h₁h₂`. As in `SPINNING_4PT.md` §4.3/4.5, this gives `+` for same and `−` for opposite, with the physical sign fixed by the covariant
    convention of §4.5 [D-sketch].
  * Hence `Ψ̃₃|_{pres–pres} = 𝒩′ g₁²S̃_LS̃_R P_{−½+iν}(cosh β)(|c_same|² − |c_opp|²)`. This has the same form as `C⁻₃`, whose normalisation is `g₁²S_LS_R/4`.
  * So `N(ν) = 𝒩′S̃_LS̃_R/(S_LS_R/4)`, up to the common constant. The odd-mode `S̃` tables (`l = 2`, odd `q₁+q₃`, antipodal weights `(−1)^{q₃}c_M`) are
    **computed in §6 [DV]**; `N` has no real zeros.
  * The MFT normalisation (`|c|²` per helicity class, `SPINNING_4PT.md` §4.6) is untouched, because it is a property of the states, not the kernel.
* **Flip pieces.**
  * Could the bare commutator also pick up states with a helicity flip at one or both pairs? These give `c_same c̄_opp`-type interference times `𝒦_{0,±4}`, `𝒦_{±4,±4}`,
    which are not sign-definite.
  * They are absent at leading order [D given Q]. The flat dictionary (CLS-D) maps `Ψ̃₃` to a combination of `D = 4` odd fixed-`u` sum rules, and the only such
    sum rule is `B₃^{(1)}`, built from the MHV `f`, which has no helicity flip at the `u`-channel pairs (§3).
  * In the overlap regime `1/Λ ≪ β ≪ 1`, `𝒦_{·,±4} → J₄(νβ)` is not suppressed. So a `β`-independent flip coefficient at leading order would contradict the flat limit
    and must vanish.
  * Caveat [I]: this uses the CLS-D dictionary (derived for `Δ, ν ≫ 1` with `ν < m`) in the overlap with the Regge regime `ν ≪ m`. That is the same use
    `ADS4_FUNCTIONAL.md` makes of it for `C^±`.

## 3. Flat limit: `ψ₃B₃^{(1)}` is reproduced [D given Q]

* **Odd sum rules in `D = 4`.** CLPSD (2.26) (label `defB`) plus the fixed-`t` rules below it list the complete basis of fixed-`u` sum rules. The odd-`k` ones are only `B^{(1)}_k = ∮ds [s(s−p²)]^{−(k−3)/2}[−f(s,−p²) + f(p²−s,−p²)]`.
  * The all-plus `g` and one-minus `h` are crossing symmetric. So `g` has no `s↔t`-odd part at fixed `u`, and `h(s,−p²) ~ s` gives no convergent odd rule.
  * So at `k = 3` the lowest-subtracted odd fixed-`u` space is one-dimensional at each `p`.
* **CLS-D map.** CLS-D (200) map the `j′`-odd bare-commutator functionals to lowest-subtracted odd flat sum rules `𝒞^{−(i)}_{4g;3,−ν²}`. Their `M⁻_{4g}` matrix is written
  for generic `D`; at `D = 4` all of its entries reduce to multiples of `B₃^{(1)}`.
* **Result.** `Ψ̃₃ → c(ν)·B₃^{(1)}`, with heavy density `(|c⁺⁺|²d̃₀₀ − |c⁺⁻|²d̃₄₄)/m⁶`, i.e. `σ_h = ±1`. After normalising by `N(ν)`, `ε₃Λ⁶C̃⁻_{3,Λp} → B₃^{(1)}(p)` with
  `ε₃ = ½`, exactly the dictionary of `ADS4_FUNCTIONAL.md` §1.2.
* **Assumed [I]:** `c(ν) ≢ 0` for the chosen dual structure; equivalently `N ≢ 0`.

## 4. `Ω′` and its light side

**Definition.** `Ω′_{g3} := ∫₀¹dp [ψ₂Λ⁴C⁺_{2,Λp} + ½ψ₃Λ⁶C̃⁻_{3,Λp}] − 5Λ¹⁰∂_{ν²}C⁺_{4,ν}|₀ + (5/2)Λ¹⁰C⁻_{5,0}`.

**Positivity.**
* At leading order in `1/m`, `Ω′|_heavy = Ω|_heavy` for every heavy state (§2). So the flat statement (Theorem 1: `F_{g3}[m,J,h] > 0`, both channels, computer-assisted)
  and the AdS leading-Regge analysis of `ADS4_FUNCTIONAL.md` apply unchanged [D given §2–3].
* The uniform remainder near the gap involves the `O(1/m²)` terms of `C̃⁻₃`, which differ from those of `C⁻₃`. Neither set is computed; this remains an open gap of
  `ADS4_FUNCTIONAL.md`, now with a different functional.

**Light side of the `k = 3` term [D given planar `J₀ < 3`, HM/CDP].**
* `λ²_TTT Ψ̃₃[G_T]` is `∝ ε†[H₁,H₃]ε = −2i|c|²sin4φ₁₃σ_z`, with `|c|² = (t₄²+α_T²)/16 = 3600|ĝ₃|²` (`t₄²+α_T² = 57600|ĝ₃|²`, `DICTIONARY_T4.md`).
* So it has no Einstein, linear or contact term. This matches the flat `−B₃|_low = −2πG|ĝ₃|²p⁴` exactly, which also has no `16πG/p²`-type term and no `g_i` (`FLAT_GRAVITON_BOUND.md` §1.3).
* By (F), `−C̃⁻₃|_light = −4πG|ĝ₃|²(ν⁴ + f̃₁ν² + f̃₀) + |ĝ₃|²G h(ν)`.

**Effect on the bound.** In `−Ω′_light = πGΛ²[16L_AdS − k_AdS x²] + …`:
* `k(0) = 330151/12012` is unchanged.
* The `C⁻₃` contribution `225f₁/(14Λ²) + 7f₀/(3Λ⁴)` to `k_AdS − k(0)` becomes the same expression with `f̃₁, f̃₀`.
* A decaying `h` contributes `(½)Λ⁶∫ψ₃h(Λp)dp ≈ 81Λ⁵∫₀^∞h`, i.e. `81G x²Λ^{−3}∫h` against the leading `k(0)πGΛ²x²`, which is relative `O(Λ^{−5})`. So property (P) is not needed for this term.
* In `Ω` (old `C⁻₃`), contacts and Einstein residues could enter the `k = 3` term. In `Ω′` they cannot.

**Constants.** `A = 1040/k(0)`, `B_AdS` and `C₁` come from `C⁺₂` and `k(0)`. They are unchanged, and still conditional on (R) and on the forward terms
(`ADS4_LIGHT.md` §4).

## 5. What remains

| item | status |
|---|---|
| `C̃⁻₃` exists, parity-odd `(2,4)` structure, sum rule valid for planar `J₀ < 3` | [Q] CLS-D, KKSDZ; number of `d = 3` members [I] |
| leading heavy action `= (2P/m⁶)(|c_s|²−|c_o|²)` after normalisation by `N(ν)` | [DV] at `(1,3)` (§6): `N ∝ (2ν+i)(4ν²+9)/(4ν²+201)`, no real zeros, no flip; `(2,4)` kernel `= 𝔄` by choice [D, §7], so `N₂₄ = 1` and no flip at `(2,4)` |
| flat limit `= ψ₃B₃^{(1)}`, Theorem 1 applies | [D given Q] (CLPSD basis, CLS-D map); `c(ν) ≢ 0` [I] |
| light side: no Einstein/linear/contact/`[TT]` term; top coefficient `−4πG|ĝ₃|²ν⁴` | [D given planar `J₀ < 3`, HM/CDP, (F)] |
| `f̃₁, f̃₀`; `O(1/m²)` heavy terms | [ND] |

Next computation, single-core and light: the odd-mode `M`-tables `S̃` for `l = 2` at antipodal weights. These use the `laplace.py` machinery of `ads4/code/e2e/` and give `N(ν)`
in closed form, as products of the rational `S` tables and `g₁`.

## 6. Computation of `N(ν)` and the flip check (`lightinputs/t6_normalisation.py`, `t6b_flip_scan.py`; logs `cache/t6*.log`)

**Method [DV].** This uses the exact `M`-tables of `ads4/code/e2e/laplace.py` (imported; they were cached before). They give the Fourier–Laplace response
`M^{q₁q₃}_{h₁h₃} = R + g₁S`, `g₁ = γ_{a−1}/γ_a`, of each frame mode to a heavy transverse-helicity plane wave.
* The bare commutator's `(1,3)` kernel is the antipodal combination `w_q = (−1)^{q₃}d_{(q₁+q₃−1)/2}(δ)/d₀`, with `q₁+q₃` odd and `δ = ½+iν`.
  These are the Fourier coefficients of `sgn(sin ψ)|sin ψ|^{−δ}`, re-checked numerically in `e2e/clsd_kernel.py`.
* The project's odd kernels are `D₅` (`l = 1`) and `𝔄 = D₅D₁` (`l = 2`), with mode content `(±1,∓1)` and `(±2,∓2)`.
* Response per helicity flow: `A(h₁,h₃) = Σ_q w_q M^q_{h₁h₃}`, at `a = 2l+k+1` with `k = 3`.

**Results (exact):**

| | `l = 1` (`a = 6`) | `l = 2` (`a = 8`) |
|---|---|---|
| project `A(±,∓)` | `±(5i/8)g₁` (∝ the tabulated `∓i(a−1)g₁/2`) | `±(7i/64)g₁(4ν²+201)/(4ν²+81)` (∝ tabulated `𝔄` entry) |
| bare commutator `A(±,∓)` | `±g₁(iν/2 − ¼)` | `±g₁(1/648)(−32iν³+16ν²−72iν+36)/(ν²/20.25+1)` |
| pres even part / flip responses | 0 / 0 for both kernels | 0 / 0 for both kernels |
| `N(ν) = A_CLSD/A_proj` | `(2/5)(2ν + i)` | `−(32/7)(2ν+i)(4ν²+9)/(4ν²+201)` |

* **Structure.** Both kernels are purely odd (`A(+,−) = −A(−,+)`), so the two-pair product weighs `h₁h₂`. That gives the `σ_h` structure `|c_s|² − |c_o|²`
  of `C⁻₃`, with the same harmonic `P_{−½+iν}(cosh β)`. The heavy actions are therefore proportional, with proportionality factor `N(ν)` [DV, given the
  two-pair heavy formula of `SPINNING_4PT.md` §4.2 and the `(2,4)` choice].
* **The `l = 1` row is a consistency check.** `N` is linear and the project entry is the known `D₅` row.
* **Zeros of `N`.** `|N(ν)|² ∝ (4ν²+1)(4ν²+9)²/(4ν²+201)²` is strictly positive on real `ν`. At `ν = 0` it is `(32/7)·9/201·1 ≈ 0.205` in these units; it grows
  like `(64/7)ν` at large `ν`. **`Ω′` does not fail anywhere on the real axis** [DV].
* **Phase.** `N` is complex, `arg N = arctan(1/2ν) + const`. The normalised `C̃⁻₃ = Ψ̃₃/N` has a real leading heavy action. Its light side is `λ²Ψ̃₃[G_T]/N`, which is
  analytic on real `ν`.
* **Remark.** `N` vanishes at `ν² = −9/4`, the `T`-exchange pole location. So `1/N` has a pole there. On the real axis this has no effect; whether it
  cancels against a zero of `Ψ̃₃[G_T]` was not examined.

**Flip check [DV].**
* `t6b` lists, for every frame mode, the parity-odd parts of its preserved and flip responses. For `l = 2`, 20 of 25 modes have a nonzero odd flip response.
  So flip-freedom is a property of the specific kernels, not of the tables.
* For both the antipodal bare-commutator kernel and the project's `D₅`/`𝔄`, `A(+,+) = A(−,−) = 0` identically in `ν` and `a`.
* This replaces the earlier "[D given Q + overlap matching]" argument by a direct computation at the `(1,3)` pair.

**What remains [ND].**
* The `O(1/m²)` heavy terms.
* `f̃₁`, `f̃₀`.
* The two-pair heavy formula for stress tensors is itself conditional (`SPINNING_TTTT.md`).

## 7. The `(2,4)` side: any covariant kernel is allowed, in particular `𝔄` [D]

* **Why the kernel is free.** The bare-commutator sum rule rests on `⟨T(x₄,z₄)[L[T](x₀,z₃), L[T](x₀,z₁)]T(x₂,z₂)⟩ = 0`. This holds **pointwise** in `x₂, x₄, z₂, z₄`
  whenever the integral converges (KKSDZ; planar `J₀ < 3`). So integrating against *any* kernel at the `(2,4)` pair gives a valid sum rule. CLS-D's dual
  structure `(⟨T₄L[O]T₂⟩^{(a)})^{−1}` is one convenient covariant choice; it is not forced.
* **Which kernels are covariant.** Such a kernel (spectral parameter `ν`, i.e. a three-point structure of the `(T₂, T₄, X_ν)` representations with `X` at `x₀`, polarisation `z₀ = ℓ`)
  is, at leading Regge order, that structure evaluated at `(x₂, x₀, x₄) = (−x̂, 0, x̂)`. (In the CHMRS frame `x₀ = |y|n̂ → 0`.)
  * The stabiliser of three points on a timelike line in `d = 3` is the `SO(2)` of rotations about `x̂`.
  * The null polarisation `z₀` (mod scale) is a point `χ` of the celestial circle at `x₀`, on which this `SO(2)` acts transitively.
  * Hence the three-point structures, both parities, are exactly the `SO(2)`-covariant polynomials of bidegree `(2,2)` in `(z₂, z₄)` written in the frame
    `(x̂, n_ℓ, n_ℓ^⊥)`, i.e. arbitrary combinations of the 25 frame modes `k_{q₂q₄}`.
  * This is the familiar count of `(2l₁+1)(2l₂+1) = 25` structures for spins `(2, 2, J ≥ 4)` in `d = 3`. The homogeneity in `z₀` only fixes the factor `(−x̂·ℓ)^{−½+iν}`.
* **The asymmetry with `(1,3)`.** At `(1,3)` the light transforms force antipodal polarisations, which is why `𝔄` cannot be used there (`LIGHT_INPUTS.md` §3.1).
  At `(2,4)` the operators are not light-transformed, so no such constraint arises.
* **Conclusion.** Choose the `(2,4)` kernel to be the structure whose Regge restriction is `𝔄 = D₅D₁`, exactly the project's `C⁻` kernel at that pair. Then:
  * the `(2,4)` factor of `N` is identically 1, so **`N(ν)` of §6 is complete**;
  * `𝔄` has zero flip response (`A(±,±) = 0`, §6 [DV]), so **no `c_same c̄_opp` (flip) pieces arise at `(2,4)` either**;
  * none arise at `(1,3)` (§6).

  So there is nothing sign-indefinite to quantify at leading order.
* **Effect of a different `(2,4)` kernel.** Any other choice multiplies `N` by `S̃_R/S_R` and may introduce flip terms. This one avoids both.
* **Not checked.** That the chosen structure is "nice" beyond leading Regge order (its subleading Regge moments enter only the `O(1/m²)` heavy terms, which are [ND] anyway).
