# H0: the spinning (`⟨TTTT⟩`) functional — route decision and first milestone

Status: 18 Nov 2026, derivation/numerics agent. M1 (§0–§8) committed (ce46a54); M2/M3 (§9–§11) committed (7211186); §R–§14 committed (78eca6a), §15–§16 (29e98d1), §17–§18 (1c7a156), §19 (a2618ed); §20.1–§20.7 (ce8480a); **§20.8 (step 2, round 2) and §21 (round 3, partial) unrefereed, not committed.** Targets hypothesis H0 of `RESULT_D.md` (with H2, H4, H14,
H9b) and finding F1 of `REFEREE_FOUNDATION.md`. Code: `ads4/code/spinfunc/` (`run_all.sh`, single core, ≈ 1–2 h with M3b; outputs in
`spinfunc/cache/`, git-ignored). It imports `omegau/`, `exactskel/` (e1/e7 caches) and `objid/` **read-only**. Nothing else was modified.

Source read for this round: Chang–Landau–Simmons-Duffin, arXiv:2311.04271 (**CLS-D**), LaTeX source (`draft.tex`, §2–§4.6, App. B),
downloaded to the session scratchpad. Equation labels below are CLS-D's LaTeX labels.

Marks: **[Q]** quoted; **[D]** derived; **[DV]** derived and checked (exact algebra, sympy); **[D-sketch]**; **[NE]** numerical evidence;
**[M]** model; **[I]** assumed/inferred.

> **Update (M2 and M3, §9–§11; same day, unrefereed, not committed). Several M1 statements below are superseded.**
>
> 1. **In `d = 3` the bulk-point images of CLS-D-class functionals are the STRIPPED flat sum rules `B^{(1)}_k` (reading B̃, `d̃^J_{44}`), not
>    the unstripped `B′_k` (reading B) that `RESULT_D.md`'s near-gap model and M1 used.** Three independent arguments (§9):
>    * an explicit `D = 4` crystal computation with CLS-D's polarisation-sum prescription: every flip-free kernel has channel ratio
>      `(1−x²)^{−4}` [DV];
>    * CLS-D's unsubtracted `k = 2, 3` graviton sum rules have truncating EFT series, which fits the contact-free stripped `B^{(1)}_{2,3}`
>      and excludes the contact-carrying `B′_{2,3}` [D from Q + `FLAT_UNSTRIPPED.md`];
>    * `f^±_{k;2}|crystal` is exactly the stripped kernel [DV].
> 2. **Route A′ is therefore redefined (§10): `Ω^{A′} := ∫ψ₂Ĉ₂ + ∫ψ₃Ĉ₃`**, with only CLS-D's two unsubtracted functionals and the NF02
>    weights `ψ₂, ψ₃`.
>    * Its flat limit is NF02 itself, in the stripped basis, so the certificate applies [Q + D].
>    * `k = 4…7`, odd-`k` modelling (H3), H9b and H15 are not needed.
>    * Only `k = 2` needs an `η`-inversion (`j = 1`).
> 3. **The B′-uplift (`Ω^{NF,SC}`, `k = 2…7`) is not positive at leading order once its objects are CLS-D-class (stripped) [NE].**
>    * The flat action has min `F/N_s = −0.72` at `+−`, `m/M = 1.1`, `J = 8` (AdS-reachable at `Δ_g = 100`) and `−3478` at `m = M`, both
>      with the project's `Ψ_k` (`R_k`) and with CLS-D weights (`f7`).
>    * This supersedes §4.1 and the `PsiG` columns of §4.2. The `R_k` of §3 stand, but relative to the stripped `B^{(1)}_k`.
>    * The `A′` (`R = 1`, B′) rows of §4.2–4.4 model a functional that does not exist. The tolerance numbers of §4.3 remain valid as
>      sensitivity statements of the near-gap model.
> 4. **Near gap, route A′** (`f6`, same exact-skeleton model, stripped densities) [NE + M]:
>    * `n = 1` at `Δ_g = 100 / 300 / 1000`;
>    * min `F/N_s = 0.2256 / 0.2054 / 0.2021`, and `0.2231 / 0.2052 / 0.2021` including the `η`-inversion estimate.
> 5. **M3** [D-sketch + NE]:
>    * The graviton-pole light term is unchanged by an exact (Mehler–Fock) `η`-inversion.
>    * `k(0)` shifts by `−1.5·10⁻³ / −1.6·10⁻⁴ / −1.5·10⁻⁵`.
>    * So `A`, `B_AdS` scale by `1.0015 / 1.00016 / 1.00002`, and the bounds by `≤ 1.0008`.
>    * **The headline numbers (`A = 3022.3`, `B_AdS = −6552.5`, `2.21·10⁻⁴ / 3.13·10⁻⁶ / 2.89·10⁻⁸`) are unchanged at the quoted precision for
>      route A′**, in the same model and under the same hypotheses minus those listed in §11.

## R. The route-A′ functional `Ω^{A′}`: self-contained statement (for review; M2b/M3b/M4′ details in §12–§14)

**Setting.**
* Unitary CFT₃ with stress tensor `T` (`Δ = 3`, `J = 2`), large `c_T`, twist gap `Δ_g`; `Λ := √((Δ_g−2)(Δ_g−1))`.
* Correlator `⟨TTTT⟩`, cross-ratios `u′ = x₁₂²x₃₄²/(x₁₃²x₂₄²)`, `v′ = x₁₄²x₂₃²/(x₁₃²x₂₄²)`.
* Heavy states `(Δ, J)` with `m² = α(α+2J+1)`, `α = Δ − J − 2` (`omegau/common.state`), `u := Λ²/m²`, `η = cosh β` as there.

**1. Building block** [Q, CLS-D arXiv:2311.04271 eq. `spinningfunctional_definition` and §3.5.2].
`Ψ^{(a),(c)}_{k,ν,j′}[f G]` integrates `f(u′,v′)·⟨T(x₄)[L[T](x₀,z₃), L[T](x₀,z₁)]T(x₂)⟩` against:
* the inverse ("dual") structure `(a)` of `⟨T L[𝒪] T⟩`, where `𝒪` has dimension `3/2 + iν`, spin `J₁+J₃−1 = 3` and transverse spin `j′`;
* the celestial structure `(c)` (unique for each `j′`; in `d = 3`, `j′ ∈ {0, 1}`).

It vanishes on every physical correlator for which the integral converges: the two light-ray operators are on a common null plane, so
their commutator vanishes. Regge convergence: `J₀ < 2` for both functionals below.

**2. The two functionals.**

| | `Ĉ^{CLS}_{2,ν}` | `Ĉ^{CLS}_{3,ν}` |
|---|---|---|
| transverse spin | `j′ = 0` | `j′ = 1` |
| subtraction factor | `f = v′ − u′` (`f⁺_{2;2}`, unsubtracted) | `f = −1` (`f⁻_{3;2}`, unsubtracted) |
| `(2,4)` structure | `(a₊)`: the `j′ = 0` dual-structure combination whose bulk-point image is `Q⊗Q` | `(a₋)`: image `O⊗Q` |
| leading Regge kernel (`r^{k−2}` × weight × harmonic) | `½ η r⁰ P_{−½+iν}(η)` | `−r¹ P_{−½+iν}(η)` |

Covariantly, `(a₊) = q²` with **`q := (2H₂₀H₄₀ + V_{0,42}² H₂₄)/(X₂₀X₄₀)`** (photon-level, CLS-D's `H`, `V`, `X` of eq. `CFT_building_blocks_VH`;
§17, [DV]). Its crystal image is `−Q/(4ν²(m²−ν²)²)` for all `x`, where
`Q := (𝒲₂·e₂)(𝒲₄·e₄) − (𝒱₂·e₂)(𝒱₄·e₄)` (∝ CLS-D's photon `𝒦^{(1)}_{j′=0}` at `D = 4`), `O := (𝒲₂·e₂)(𝒱₄·e₄) − (𝒱₂·e₂)(𝒲₄·e₄)`, and
`𝒱`, `𝒲` are CLS-D eq. `kernel_VWvec_definition`. Both images are flip-free at the crystal for all `x` (§12) [DV].

**3. Normalisation.** `Ĉ_k := Ĉ^{CLS}_k/N̂_k(ν)`, with `N̂_k` fixed by the leading Regge heavy action. On a heavy block, in the per-state
normalisation of `omegau/common.py` (positive state factors divided out), `ε_k m^{2k}Ĉ_{k,ν} → 2P_{−½+iν}(cosh β)` at leading Regge order.
Here `ε₂ = 1`, `ε₃ = ½`, and the `k = 3` term carries `σ = +1 (++)`, `−1 (+−)`.

**4. Weights.**
* `ψ₂(p) = (1−p)⁵ p X(p)/k₀` and `ψ₃(p) = (1−p)⁵ p Y(p)/k₀`, `p ∈ [0, 1]`.
* `X`, `Y` are degree-8 polynomials with 9 exact rational coefficients each (`ads4/code/robust/functionals_noforward.json`, key `NF02`).
* `k₀ = 27.48509817…`, chosen so that `k(0) := −2∫₀¹ψ₂p⁶dp + 2∫₀¹ψ₃p⁴dp = 1` (`omegau/nf_common.NF`).

**5. The functional.**

`Ω^{A′} := ∫₀^∞ dp ψ̃₂(p) Λ⁴ Ĉ_{2,Λp} + ∫₀¹ dp ψ₃(p) ½Λ⁶ Ĉ_{3,Λp}`,

where `ψ̃₂ = 𝔍[ψ₂]` is the **exact `η`-inversion**:

`Φ(t) = ∫₀^Λ ψ₂(ν/Λ) cos(νt) dν`,
`Φ̃(t) = Φ(t)/cosh t − ½ tanh t · sech^{½}t ∫₀^t Φ(s) sech^{½}(s) ds`,
`ψ̃₂(ν/Λ) = (2/π)∫₀^∞ Φ̃(t) cos(νt) dt`.

* **What the inversion guarantees [D].** `∫ψ̃₂ Ĉ^{CLS}_2` has exactly the leading-Regge `η`-kernel `∫ψ₂ P_{−½+iΛp}(η)dp` of the pure-harmonic
  object, for every `η ≥ 1`.
* **Size and support.** `ψ̃₂ − ψ₂` is `O(Λ⁻²)` in the interior. `ψ̃₂` has a tail at `p > 1` of relative mass `2·10⁻⁹` (`Δ_g = 100`).
* **`k = 3` needs no inversion** (`j = 0`).

**6. Per-state heavy action** (used in all scans): `F(m, J, h) = Σ_{k=2,3} ∫dp ψ_k(p) u^{k−2} ε_k m^{2k}Ĉ_{k,Λp}[state]`, with `ψ₂ → ψ̃₂`.

**7. Flat limit [Q + DV, §9, §12].**
* As `Λ → ∞` at fixed `m/Λ`, `J`, `F → ∫ψ₂B^{(1)}_2 + ∫ψ₃B^{(1)}_3`. These are CLPS's (arXiv:2201.06602, (2.26a,b)) stripped sum rules for the
  MHV `f(s,u)`.
* Densities: `(2−y)(|c⁺⁺|²P_J + |c⁺⁻|²d̃^J_{44})` and `u σ(|c⁺⁺|²P_J − |c⁺⁻|²d̃^J_{44})`, with `y = p²/(m/Λ)²` and `d̃ = d/(1−y)⁴`.
* This is NF02, with flat min `F/N_s = 0.200`, certified (`NOFORWARD_FUNCTIONAL.md`; one `+−` corner float).

**8. Light side (tree level, `1/c_T`) [D given Q / H7, §14–§15].**
* Both functionals are unsubtracted, so they have double zeros at every `[TT]` twist (even integers `≥ 2`). The light side is therefore
  the identity (0) plus the single `T` block.
* `Ĉ₂` gives the graviton pole `16πG/(ν² + 9/4)` (marginal `J₀ = 2`) and the `|ĝ₃|²` term. `Ĉ₃` is Einstein-free (`A₃ = f⁻_{3;2}`,
  `OMEGA_PRIME.md`) and gives `|ĝ₃|²`.
* Flat images: `−B₂ = 16πG/p² + 2πG|ĝ₃|²p⁶`, `−B₃ = −2πG|ĝ₃|²p⁴`. There are no contact terms.

**9. Bound.**
* `k(0)_{A′} |ĝ₃|² Λ⁸ ≤ 16 n⁸ L_A`, with `L_A = ∫₀¹ ψ₂(t) dt/(t² + 9/(4L²))` and `L = Λ/n`.
* `L_A` uses `ψ₂`, not `ψ̃₂`: the inversion leaves the graviton-pole term exactly invariant (§13, checked to `5·10⁻¹⁴`).
* `√(t₄²+α_T²) ≤ 240 √(|ĝ₃|²Λ⁸)/Λ⁴`.
* `k(0)_{A′} = 1 + δk(0)`, with exact `δk(0) = −1.507·10⁻³ / −1.638·10⁻⁴ / −1.463·10⁻⁵` at `Δ_g = 100 / 300 / 1000`.
* `A = 3022.3/k(0)_{A′}`, `B_AdS = −6552.5/k(0)_{A′}`.

**10. What it rests on.**
* CLS-D's bulk-point factorisation and its `d = 3` specialisation (§9, §12).
* The graviton-level covariant realisation of `Q⊗Q`, `O⊗Q` ([D-sketch], §12).
* Convergence and swappability (§14).
* H4 (finite-`m` remainder transfer), H12′ (binding state `+−`, `J = 5`, parity-odd), H13′ (off the scanned set), H11, and tree level.

## 0. Summary (M1; see the update box above for what is superseded)

1. **Route decision: A, in a concrete form ("A′").** The spinning C-type functional the project has been calling "not constructed" is,
   up to one step, already in the literature:
   * CLS-D define light-ray-commutator sum rules `Ψ^{(a),(c)}_{k,ν,j′}` for external spinning operators (stress tensors included), with
     **meromorphic** subtraction factors `f^±_{k;J_e}(u′,v′)` (`subtraction_factor_evenjprime`, `…oddjprime`). Their scalar version
     *is* CHMRS's `C̃_{k,ν}` (CLS-D footnote to `scalarfunctional_nokernelexpression`) [Q].
   * They prove (saddle point, general `d`, `ν, m ≫ 1` at fixed `ν/m < 1`) that the heavy action of these functionals equals the flat
     graviton sum rule (shock amplitudes; the lowest-subtracted fixed-`u` sum rules of CLPS) [Q].
   * **The project's `Ψ_k` are exactly these CLS-D functionals with the subtraction factor replaced by the weights `A_k`**
     (`SUPERCONVERGENT_UPLIFT.md` item 1); `A₃`, `A₇` coincide with CLS-D's, `A₄, A₅, A₆` do not. The replacement was made to get pure
     harmonics in the Regge limit (CLS-D's factors carry `η^j`, `j = 1,0,3,2,1,0` for `k = 2…7`). It is what breaks meromorphy (H9b) and
     the flat limit (F1).
   * **Route A′:** `Ĉ_k :=` CLS-D's `Ψ_{k,ν,j′; f^±_{k;2}}` (same light-ray operators, same dual/celestial structures as the project's
     `Ψ_k`), followed by the CHMRS `η`-inversion in `ν` (Legendre recursion), for **all** `k = 2…7`, including `Ĉ₂` (CLS-D's unsubtracted
     `f⁺_{2;2} = v′−u′`). Route B is dominated (§1).
2. **Milestone M1 executed: the graviton bulk-point dictionary of both objects, exactly in `x = ν/m`** (§3). Using CLS-D's factorisation
   (the subtraction factor enters the bulk-point action only through its value at the "scattering crystal"):
   * the crystal is the CHMRS Regge-frame point `η = 1`, `r = r* = (1−s)/(1+s)`, `s = √(1−x²)`: all 9 printed CLS-D saddle values
     `f^±_{k;J_e}|crystal` are reproduced **[DV]**;
   * cross-check: the same method reproduces `OBJECT_IDENTITY.md`'s independently derived (and refereed) scalar density ratio
     `(1−r*²)^{k+2}/(1+r*²)` for `k = 2, 4, 6` **[DV]**;
   * **graviton result [D given CLS-D]:** relative to B′ (= CLS-D's flat sum rule = route A′ at leading order), the bulk-point density of
     the project's `Ψ_k` is `R_k(x)`:

     | `k` | `R_k` (in `r = r*`) | in `s = √(1−x²)` | `x = 0.5 / 0.7 / 0.9 / 0.99` |
     |---|---|---|---|
     | 3, 7 | 1 | 1 | 1 |
     | 4 | `(1−r²)²/(1+r²)` | `8s²/((1+s)²(1+s²))` | 0.985 / 0.920 / 0.620 / 0.120 |
     | 5 | `(1+6r²+r⁴)/(1−r²)²` | `(1+s⁴)/(2s²)` | 1.042 / 1.235 / 2.73 / 25.1 |
     | 6 | `(1−r²)⁴/(1+r²)` | `128s⁴/((1+s)⁶(1+s²))` | 0.975 / 0.869 / 0.443 / 0.023 |

     `R_k − 1 = O(x⁴)` (`R₄ = 1 − 3x⁴/16 − …`, `R₅ = 1 + x⁴/2 + …`, `R₆ = 1 − 5x⁴/16 − …`).
   * **The scalar-analogue model of `OBJECT_IDENTITY.md` was wrong for gravitons**, qualitatively at `k = 5` (it has `Ψ₅`'s density → 0 as
     `x → 1`; the graviton `Ψ₅` density diverges like `1/(2(1−x²))`), and by four powers of `(1−r²)` at `k = 4, 6`.
3. **Consequences [NE + M]** (§4):
   * **F1 sharpened: the flat limit of the functional as currently defined is not positive.** With the graviton densities, the D-safe flat
     action is negative for `+−`, `m/M < 1.0003`, `J = 16–24` (min `F/N_s = −0.0162` at `m = M`, `J = 20`; quadrature-independent, 600–6000
     nodes). The cause is the `k = 5` term. Under the twist gap (`m_min(J)/Λ = √(1+2J/(Δ_g−1))`) these states exist only for
     `Δ_g ≳ 6·10⁴`, so the `Δ_g = 100–1000` numbers are not affected at leading order, but the current functional is **not** the uplift of
     a positive flat functional, and the NF02 certificate does not apply to it. **Route A′'s flat limit is NF02 exactly**, so the
     certificate applies [Q + D].
   * **Near gap, `Δ_g = 100 / 300 / 1000`** (bulk-point-factor model on the exact scalar skeleton, all 1745 / 1223 / 1042 records):
     `n = 1` for both objects. Min `F/N_s`: route A′ `0.2267 / 0.2057 / 0.2022` (= `RESULT_D.md`, by construction); current `Ψ_k` with
     graviton densities `0.2281 / 0.2077 / 0.2042`. Largest object shift `0.035 / 0.052 / 0.121 N_s`.
   * **Where the near-gap action lives** (route A′; worst-sign and coherent tolerances, §4.3): the region `x > 0.7` carries a median
     `≤ 0.5 %` of `Σ|density|`; a coherent relative error of the graviton leading density on `x > 0.7` of `±22 %` (`Δ_g = 1000`, worst
     of three sign patterns) to `±180 %` (`Δ_g = 100`) does not move `n = 1`. The load-bearing window is **`x ∈ [0.3, 0.7]`**: there a
     coherent density error of `≈ 4.5–8 %` (`x > 0.3`) or `≈ 8–30 %` (`x > 0.5`) would break `n = 1`. The graviton-vs-scalar remainder (H4) may differ coherently by a
     factor `≈ (−1.2, 2.6)` on `x > 0.3` at `Δ_g = 100`, and by factors `> 10` at `Δ_g ≥ 300`.
   * **`η`-inversion cost** (large-`ν` symbol estimate [D-sketch + NE]): `≤ 0.087 / 0.0088 / 0.0012 N_s` (`∝ 1/Λ²`); `n = 1` unchanged.
4. **Impact on the hypotheses** (§6): H0 goes from "object not constructed" [I] to "object identified in the literature (general `d`);
   `d = 3` specialisation, `η`-inversion analyticity, convergence and light side to be shown" [D-sketch]. H2 at leading order in `1/m` and
   fixed `x < 1` becomes [Q] for route A′ (conditional on the `d = 3` specialisation, M2). For the current `Ψ_k` the B′ assumption is
   **false** at `x ~ 1` for `k = 4, 5, 6` (densities `R_k`·B′). H4 is unchanged ([I]) but localised: `x ∈ [0.3, 0.7]`, and at
   `O(1/m²)` it is a computable fluctuation correction around the CLS-D crystal (M6).
5. **The check suggested for this round ("verify `Ĉ₄` against the existing `Ψ₄`/`Π̂` data at `x ≲ 0.1`") cannot discriminate.** All
   candidate objects agree to `O(x⁴)`: `|R_k − 1| ≤ 5·10⁻⁵` at `x = 0.1`, against `10⁻³` data precision. The leading-Regge identity of `Ĉ_k`
   with `Ψ_k` holds by construction; it was not tested numerically here.

## 1. Route decision

| | Route A (spinning `C`-type `Ĉ_k`) | Route B (keep `Ψ_k/N_k`) |
|---|---|---|
| object | CLS-D `Ψ_{k,ν,j′}` with CLS-D's meromorphic `f^±_{k;2}`, then `η`-inversion | CLS-D `Ψ_{k,ν,j′}` with the project's `A_k` (non-meromorphic `r^{k−3}` for even `k`) |
| literature | construction, validity (commutator vanishing; Regge convergence `k > J₀` for even `j′`, `k − 1 > J₀` for odd `j′`) and bulk-point dictionary to flat graviton sum rules: CLS-D [Q]; `η`-inversion: CHMRS App. superconvergence (scalar) [Q] | same functional class; nothing on non-meromorphic weights |
| flat limit | **NF02 exactly** at leading order (CLS-D) → certified (`NOFORWARD_FUNCTIONAL.md`) | B′ × `R_k` (§3): **not** NF02, **not positive** (§4.1); needs new weights and a new certificate. And for even `k` the flat image is not a meromorphic dispersive sum rule (`R₄`, `R₆` contain `√(1−x²)`), so it is not a standard flat functional at all |
| light side (`k(0)`, `A`, `B_AdS`, graviton pole, contacts, `[TT]`) | CLS-D/CHMRS-type: meromorphic weight ⇒ `dDisc` form, finitely many subtracted double-trace families, polynomial contacts ⇒ **(P)/H9b analysable** [D-sketch]; flat-limit light side by the same dictionary (CLS-D: "it follows that the contributions of light states must match") | even `k`: (P) **unanalysable** (no `dDisc` rewriting; H9b stays [I] forever); the flat image of the light side is not a standard EFT sum rule |
| exact near-gap value (H2/H4 at `x ~ 1`, finite `m`) | not available; but: (i) leading order in `1/m` is [Q]; (ii) `O(1/m²)` is a fluctuation integral around the crystal with algebraic slow factors (M6); (iii) Mellin route possible in principle (meromorphic weights) | not available even for scalars (non-meromorphic weight: no Mellin form, `OBJECT_IDENTITY.md` §2.4) |
| cost to a statement at the level of `RESULT_D.md` | M2–M5 below: ≈ 5–8 weeks (§7) | re-optimise flat weights, certify, recompute light side with non-standard flat sum rules; (P) still open; ≈ 4–6 weeks for a weaker endpoint |

**Decision: route A′.**
* **Why A is feasible now.** The two hard components were thought to be missing: the spinning functional and its flat limit. Both are in
  CLS-D. What is left is the `d = 3` specialisation (CLS-D work in general `d`, with `D ≥ 8` for their amplitude basis), the
  `η`-inversion (an explicit difference equation in `ν`), convergence at the single-ray endpoints, and the light side.
* **Why B is dominated.** Its only advantage, a pure-harmonic Regge kernel without inversion, was bought by changing the subtraction
  factor. That change is what destroys:
  * meromorphy (H9b);
  * the identity of the flat limit with NF02 (F1);
  * and, as found here, flat positivity itself.

  Route B would have to repair all three. Route A′ keeps the CLS-D factor and moves the Regge-side harmonic purity into the
  `ν`-inversion, at an `O(1/Λ²)` cost (§4.4).
* **What neither route provides cheaply:** exact finite-`m` graviton moments at `x ∈ [0.3, 0.7]` (H4). Position space needs 50–90 digits
  there (`SPINNING_MELLIN.md` §1.1). The realistic next step is the `1/m` expansion around the CLS-D crystal (M6), not exact Mellin.

## 2. Route A′: definition and properties

**Definition.** For `k = 2…7`, `j′ ≡ k (mod 2)`, with `(a)` the dual structure at `(2,4)` realising the project's flip-free kernel (`K*` or
`𝔖_ff` for even `k`, `𝔄` for odd `k`) and `(c)` the celestial structure at `(1,3)` of the project's `Ψ_k`:

`Ψ^{CLS}_{k,ν} := Ψ^{(a),(c)}_{k,ν,j′}[f^±_{k;2}(u′,v′) G]`  (CLS-D `spinningfunctional_definition` with `subtraction_factor_…jprime`, `J_e = 2`),

with `f⁺_{2;2} = v′−u′`, `f⁻_{3;2} = −1`, `f⁺_{4;2} = ((v′−u′)³ − (v′−u′))/(4u′v′)`, `f⁻_{5;2} = −((u′−v′)²−1)/(4u′v′)`,
`f⁺_{6;2} = (v′−u′)/(u′v′)`, `f⁻_{7;2} = −1/(u′v′)`.

* **Leading Regge kernel [DV, `f1` item (3)].** `f_{k;2} = c_k η^j r^{k−3}(1 + O(r²))`, with `(c_k, j) = (½,1), (−1,0), (8,3), (−16,2),
  (128,1), (−256,0)`. The project's `A_k` are `r^{k−3}(1 + O(r²))` with no `η`.
* **`η`-inversion [D-sketch].** In `d = 3`, `ηP_ν = ½(1 + 1/(2iν))P_{ν−i} + ½(1 − 1/(2iν))P_{ν+i}` (`P_ν ≡ P_{−½+iν}(η)`). For a
  `ν`-smeared functional `∫dν φ(ν) X_ν` (`X`: pure-harmonic leading kernel, the object the D-safe weights `φ_k` refer to), the same leading
  kernel is obtained from `∫dν φ̃(ν) Ψ^{CLS}_{k,ν}` with `φ̃` solving
  `φ(ν) = ½φ̃(ν+i)(1 + 1/(2i(ν+i))) + ½φ̃(ν−i)(1 − 1/(2i(ν−i)))` (iterated `j` times). Equivalently, `φ̃` is the Mehler–Fock transform of
  `K_φ(η)/η^j`, which exists because `η ≥ 1`. At large `ν`, `φ̃ = φ + j[φ″/2 − φ′/(2ν) + φ/(2ν²)] + O(ν⁻⁴)`. This needs analyticity of
  `Ψ^{CLS}_{k,ν}` in `|Im ν| ≤ j` (as for CHMRS `C₂`) [I].
* **Bulk-point limit.** At the crystal `η = 1`, so the inversion is the identity up to `O(1/ν²)`. The bulk-point action is CLS-D's flat
  sum rule (B′) [Q], with an inversion correction estimated in §4.4.
* **Heavy Regge action.** Same leading kernel as the project's `Ψ_k` (normalisation `N̂_k(ν) = c_k N_k(ν)`), hence the same
  `(2P/m^{2k})(|c_s|² ± |c_o|²)` form and the same [DVE] scope of `SPINNING_TTTT_EXT.md` [D-sketch]. Subleading-in-`r` structure differs
  (affects `Π̂_{k+2}`-level terms only).
* **Validity / convergence [Q, conditional].** Commutator vanishing for any meromorphic `f` (CLS-D §3.3). Regge convergence for
  `k > J₀` (even `j′`) and `k − 1 > J₀` (odd `j′`), CLS-D §3.5.2. Endpoint convergence for `3 ≤ J₀ < k` (the `|ĝ₃|²` question, H15) is not
  addressed by CLS-D [I] (M4).
* **Light side.** `Ĉ₂` uses the unsubtracted `f⁺_{2;2}`: zero on all double traces; the graviton pole arises from the marginal `J₀ = 2`, as
  for CHMRS `C₂` [D-sketch]. For `k ≥ 3` the weights are meromorphic, so the CHMRS/CLS-D double-trace pattern applies and (P) is
  analysable [D-sketch] (M5).

## 3. Milestone M1: the graviton bulk-point dictionary (`f1_bulkpoint.py`)

**The factorisation used [Q].** CLS-D §3.3, "Saddle point analysis", and §3.5.1: in the bulk-point limit the quickly varying part of the
integrand is identical for scalar and spinning operators. The integral localises on the scattering crystal (`scattering_crystal`), and
"since the subtraction factors do not have any `ν` and `m` dependence", `f(u′,v′)` enters only through its crystal value. **Hence two
functionals with the same `(a)`, `(c)` and weights `f₁`, `f₂` have bulk-point actions in the ratio `f₁/f₂|crystal`** [D from Q]. This
holds for any weight that is smooth and nonzero near the crystal. `sgn(v′−u′)r^{k−3}` qualifies: `sgn = +1` on the moment domain,
`SUPERCONVERGENT_UPLIFT.md` `s3`.

**Normalisation.** Both `Ψ_k/N_k` and `Ĉ_k/N̂_k` are normalised by their leading Regge kernel. The Regge and bulk-point regimes overlap at
`x → 0` (CHMRS, scalar case; CLS-D expect the same for spin). So the density of `Ψ_k` relative to `Ĉ_k` is
`R_k(x) = [A_k/f_{k;2}](η=1, r*(x)) / [A_k/f_{k;2}](η=1, r → 0)`.

**Checks [DV] (`cache/f1_bulkpoint.txt`):**
1. **Crystal ↔ frame.** With the CHMRS frame `u′ = (1−2ηr+r²)²/(16r²)`, `v′ = (1+2ηr+r²)²/(16r²)` at `η = 1`, `r = r*(x)` and CLS-D's
   orientation (`u′ ↔ v′`), every printed crystal value `f^±_{k;J_e}|crystal` is reproduced to `≤ 4·10⁻¹⁴` relative: `J_e = 0`, `k = 2, 4, 6`;
   `J_e = 2`, `k = 2…7`. E.g. `u′v′ = s⁴/x⁸`, `|v′−u′| = (2−x²)/x²`, `(u′−v′)² − 1 = 4s²/x⁴`.
2. **Scalar cross-check.** `A = r^{k+1}` against the scalar `f^{scalar}_k` (= CHMRS `C̃_k`) gives exactly `(1−r*²)^{k+2}/(1+r*²)`
   (`k = 2, 4, 6`). This is `OBJECT_IDENTITY.md` §2.5's `Ψ/C` ratio, derived there from CHMRS's large-`ν` inversion kernels and confirmed by
   the referee (`rf4`). Two independent derivations agree.
3. **Graviton `R_k`:** the table of §0 item 2. The difference from the scalar analogue is a factor `(u′v′)^{±1}` between the `J_e = 0` and
   `J_e = 2` reference weights. At the crystal that factor is `∝ (1−r²)⁴`, which removes four powers of `(1−r²)` at even `k`. At `k = 5` the
   project's `A₅ = (u′+v′)/(32u′v′)` (pure `r²` at leading Regge order) replaces CLS-D's `η²`-carrying `f⁻_{5;2}`. At the crystal
   `(u′−v′)² − 1 → 0` as `x → 1`, hence the divergence.

**Scope [I].** CLS-D's derivation is for general `d` at `ν, m ≫ 1`, `ν < m`. The `d = 3` case (`S^{d−3} = S⁰` locus, `SO(d−2)` trivial,
`j′ ∈ {0,1}`) is not separately checked here (M2). Continuity to finite `ν` is CLS-D's stated expectation, proven for scalars by CHMRS.

## 4. Impact

### 4.1 Flat limit (`f2_flat.py`, `f2b_check.py`) [NE]

Min flat `F/N_s` over `J ≤ 80` (NF02 weights; object `A` = NF02 itself, B′ identity `≤ 10⁻¹³`):

| object | `++` | `+−` |
|---|---|---|
| A′ (`R = 1`) | +0.2000 (`m = M`, `J = 4`) | +0.2000 (`m = M`, `J = 5`) |
| current `Ψ_k`, graviton densities (`PsiG`; `C⁺₂` taken as B′) | +0.2540 | **−0.0162** (`m = M`, `J = 20`) |
| same, `C⁺₂` as `Π̂₂` (`PsiG2`) | +0.2303 | −0.0158 |
| `OBJECT_IDENTITY` scalar analogue (`PsiS`) | +0.2715 | +0.1782 |

* **Negative states.** `+−`, `J = 16–24`, `m/M < 1.000033 … 1.000284` (bisection per `J`). Identical at 600, 2000 and 6000 quadrature
  nodes and for `J ≤ 120`.
* **Mechanism.** At `m = M`, `J = 20` the `k = 5` part goes from `−0.020 N_s` (B′) to `−0.288 N_s`.
* **AdS reach.** Under the twist gap a spin-`J` state has `m ≥ m_min(J) = √((Δ_g−2)(Δ_g+2J−1))`. The window is therefore reached only for
  `Δ_g > 6.3·10⁴` (`J = 18`), `7.9·10⁴` (`J = 20`), `1.4·10⁵` (`J = 22`), `4.7·10⁵` (`J = 16, 24`).

### 4.2 Near gap (`f3_ads.py`) [NE + M]

* **Model** (as `objid/o2`): object moment = `R_X(k,x)` × the exact-skeleton model moment, applied to density, curvature and exact
  remainder. `A′` reproduces e7 (node sum vs `F_ex`: `≤ 2·10⁻⁹ N_s`).
* **Scope:** all e1 states (both skeletons, `J ≤ 10`, high-`J` samples).

| `Δ_g` | object | `n` (`δ = 0` / `0.03`) | min `F/N_s` [state] | max `|ΔF|/N_s` vs A′ |
|---|---|---|---|---|
| 100 | A′ / PsiG / PsiG2 | 1 / 1 / 1 | 0.2267 / 0.2281 / 0.2381 [`+−`, 1.05Λ, `J = 5`] | — / 0.035 / 0.058 |
| 300 | A′ / PsiG / PsiG2 | 1 / 1 / 1 | 0.2057 / 0.2077 / 0.2171 [`+−`, 1.02Λ, `J = 5`] | — / 0.052 / 0.067 |
| 1000 | A′ / PsiG / PsiG2 | 1 / 1 / 1 | 0.2022 / 0.2042 / 0.2077 [`+−`, 1.01Λ, `J = 5, 4`] | — / 0.121 / 0.123 |

`√(t₄²+α_T²)` is unchanged (`n = 1`) in this model for both objects. But only for A′ is the flat limit the certified NF02 and the light side
(`A`, `B_AdS`) the one computed (§1).

### 4.3 Where H2/H4 are load-bearing (`f3_ads.py`, `f5_coherent.py`; route A′) [NE + M]

* **Share in `x`.** The share of `Σ_{k,i}|density + curvature|` from `x > 0.5 / 0.7 / 0.9` is (median over states)
  `0.004 / 0.000 / 0.000` (`Δ_g = 100`) and `0.124 / 0.005 / 0.000` (`Δ_g = 1000`). The exact scalar remainder is `≤ 0.020 / 0.006 / 0.003 N_s`.
* **Worst-sign sufficient tolerances** (any sign per node; minimum over states). A relative error `δ` of the graviton leading density,
  or of the graviton remainder relative to the scalar one, on `x > x₀`, does not move `n = 1` if:

  | `x₀` | density `δ` (100 / 300 / 1000) | remainder `δ` (100 / 300 / 1000) |
  |---|---|---|
  | 0 | 0.003 / 0.003 / 0.002 | 0.044 / 0.27 / 1.35 |
  | 0.5 | 0.049 / 0.013 / 0.008 | 0.15 / 0.40 / 1.7 |
  | 0.7 | 0.70 / 0.086 / 0.032 | 5.3 / 1.4 / 3.1 |
  | 0.9 | 182 / 27 / 2.2 | 265 / 82 / 22 |

  The small-`x₀` density numbers are dominated by `k`-cancellations at high `J` (`J = 20–28`) and are very conservative.
* **Coherent deviations** (`F → F + ε Σ_{k, x>x₀} c_k T_k`, with `c_k = 1`, `(−1)^k` or `[k≥4]`). Interval of `ε` keeping `n = 1`:

  | `T`, `x₀` | `Δ_g = 100` | 300 | 1000 |
  |---|---|---|---|
  | density, `x > 0.3` (`c = 1`) | (−0.079, +0.065) | (−0.049, +0.079) | (−0.048, +0.076) |
  | density, `x > 0.5` | (−0.32, +0.22) | (−0.15, +0.10) | (−0.12, +0.077) |
  | density, `x > 0.7` | (−1.8, +2.3) | (−0.60, +0.63) | (−0.54, +0.62) |
  | remainder, `x > 0.3` | (−2.3, +1.7) | (−16, +49) | (−144, +340) |
  | remainder, `x > 0.5` | (−1.2, +1.9) | (−10, +6.7) | (−83, +52) |

* **Reading.**
  * The statement "98 % of the near-gap action is the flat B′ density **at `x → 1`**" should read "**at `x ≲ 0.7`**".
  * Beyond `x ≈ 0.7` the graviton density may be off coherently by `≈ 20–50 %` (`Δ_g = 300–1000`, depending on the sign pattern; the
    `k ≥ 4`-only pattern is the tightest, `(−0.26, +0.22)` at `Δ_g = 1000`), or by a factor `≈ 2` (`Δ_g = 100`).
  * The load-bearing window for both H2 (`≈ 5–30 %` accuracy needed) and H4 (`O(1)` factor at `Δ_g = 100`, loose at `Δ_g ≥ 300`) is
    `x ∈ [0.3, 0.7]`.
  * The current position-space reach is `x ≲ 0.1–0.17`.

### 4.4 Cost of the `η`-inversion (`f4_etainv.py`) [D-sketch + NE]

* **Change of weights.** With the large-`ν` symbol of §2, `|φ̃/φ − 1|` is `≤ 10⁻²` on most of `p ∈ [0.1, 0.9]` but reaches `0.4–0.7` at
  `p = 0.9`, `Δ_g = 100`. These are the `(1−p)⁵` zeros, `φ″/φ ∼ 20/(Λ(1−p))²`.
* **Action shift.** `max|ΔF|/N_s = 0.087 / 0.0088 / 0.0012` (`Δ_g = 100 / 300 / 1000`); nodes with `ν ≥ 3` only: `0.021 / 0.0058 / 0.0011`.
  `n = 1` unchanged; min `F/N_s` `0.2256 / 0.2056 / 0.2022`.
* **Reliability.** The symbol is unreliable at `ν = O(1)` and in the `p → 1` boundary layer. The exact difference-equation solution is M3.
* **Alternative.** Keep `φ_k` on `Ψ^{CLS}` with no inversion. Then the bulk-point action is NF02 exactly, but the Regge-side term
  `2ε_kσ^k P_{−½+iν}(cosh β)` becomes a combination of `P_{−½+i(ν±im)}`. Which is cheaper is to be decided in M3.

## 5. What the `x ≲ 0.1` data can and cannot test

* **All objects agree there.** `Ψ_k`, `Ĉ_k`, `Π̂` and `C` agree at `x ≲ 0.1` up to `O(x⁴)` in the bulk-point density
  (`|R_k − 1| ≤ 5·10⁻⁵` at `x = 0.1`; `|Π̂/C − 1| ≤ 5·10⁻⁵`) and up to common `O(1/m²)` terms. The `⟨TTTT⟩` position-space moments have
  `10⁻³` precision. So **no existing data can distinguish route A′ from the current objects**; the comparison is by construction
  (same leading Regge kernel).
* **What the data do test.** They test the common leading Regge formula (H1) and reading B (H2) at small `x`. The latter is now implied
  at leading order by CLS-D's block ↔ partial-wave map at the crystal (`Blockvspartialwave`).
* **What would discriminate.** An `x ≈ 0.5` test of the `k = 5` density (`R₅ − 1 = 4 %`) or of `k = 4, 6` at `x ≈ 0.7` (`8–13 %`).

## 6. Grades and effect on `RESULT_D.md`

| statement | grade |
|---|---|
| CLS-D functionals with `f^±_{k;J_e}` are valid spinning sum rules; scalar version = CHMRS `C̃`; bulk-point action = flat graviton (shock) sum rules, general `d` | [Q] |
| the project's `Ψ_k` = CLS-D functionals with `f → A_k`; `A₃`, `A₇` = CLS-D's; `A₄, A₅, A₆` ≠ | [D] |
| crystal = CHMRS frame `(η=1, r*)`; 9 CLS-D crystal values reproduced | [DV] |
| weight-ratio rule `f₁/f₂|crystal` for the bulk-point action | [D from Q] |
| scalar cross-check reproduces `OBJECT_IDENTITY` (`k = 2, 4, 6`) | [DV] |
| graviton densities `R_k` of the current `Ψ_k` | [D given CLS-D, `d = 3` specialisation [I]] |
| flat limit of the current `Ω^{NF,SC}` negative at `+−`, `m/M < 1.0003`, `J = 16–24`; AdS-reachable only for `Δ_g ≳ 6·10⁴` | [NE] (given `R_k`) |
| near-gap `n = 1` for both objects at `Δ_g = 100/300/1000`; margins of §4.2 | [NE + M] (bulk-point-factor model) |
| tolerances / load-bearing window `x ∈ [0.3, 0.7]` | [NE + M] |
| route A′ leading Regge kernel = project's (after `η`-inversion); `η`-inversion cost `≤ 0.09 N_s` (`Δ_g = 100`) | [D-sketch] / [D-sketch + NE] |

**Hypotheses (proposed regrading; `RESULT_D.md` not edited):**
* **H0** [I] → **[D-sketch]**: the object exists in the literature (CLS-D) up to the `η`-inversion. Open: `d = 3` specialisation (M2);
  inversion analyticity and normalisation (M3); endpoint convergence (M4); light side (M5).
* **H14 / F1:**
  * For the current `Ψ_k` the flat limit is now known exactly (§3). It is neither NF02 nor positive (§4.1), so F1 is **confirmed and
    sharpened**.
  * Under route A′ the flat limit is NF02 [Q + D], so the H14′ part of H0 is reduced to M2.
* **H2:**
  * For route A′, at leading order in `1/m` and fixed `x < 1`: [Q] (CLS-D), pending M2.
  * For the current `Ψ_k`: the B′ assumption is **false** for `k = 4, 5, 6` at `x ≳ 0.5`. In the model this does not change `n`.
* **H4:** still [I]. It is load-bearing only on `x ∈ [0.3, 0.7]`, with the tolerances of §4.3. The `x → 1` region is irrelevant to the
  threshold at the stated tolerances.
* **H9b:** analysable under route A′ (meromorphic weights) [D-sketch]; for the current even `Ψ_k` it stays [I].
* **H15** (`|ĝ₃|²`-free `Ψ_{5,6,7}`): under A′, CLS-D's Regge counting gives convergence only for `J₀ < k` (even `j′`), `J₀ < k − 1` (odd
  `j′`); single-ray endpoints still [I] (M4).

## 7. Plan for the remaining milestones

| # | milestone | content | validation | estimate |
|---|---|---|---|---|
| M2 | `d = 3` specialisation of CLS-D | `SO(d−2)` trivial, `j′ ∈ {0,1}`, `S⁰` crystal locus; 4D helicity form of the 7 + 3 lowest-subtracted sum rules; identify the `(a)` combinations that give `K*`/`𝔖_ff`/`𝔄` and the MHV `B_k^{(1)}[f]` | reproduce the project's leading Regge pair responses (`laplace.py` M-tables, `N_k`) from CLS-D's large-`ν` dual structures; reproduce reading B (`d^J_{44}`) at the crystal | 3–5 days |
| M3 | exact `η`-inversion | solve the difference equation for `j = 1, 2, 3` (or Mehler–Fock of `K_φ/η^j`); analyticity strip of `Ψ^{CLS}_ν`; support of `φ̃` beyond `p = 1`; `N̂_k(ν)` zeros; recompute `k(0)`, `A`, `B_AdS`, H17 | the large-`ν` symbol of §4.4; scalar `C₂` = CHMRS identity | ≈ 1 week |
| M4 | validity | endpoint convergence of `Ψ^{CLS}` with `f^±_{k;2}` for `J₀ ≤ 2` (Einstein-free for `k ≥ 4`) and `J₀ ≤ 4` (H15); I14 swappability | KKSDZ / CLS-D counting; scalar CHMRS case | 1–2 weeks (to [D-sketch]) |
| M5 | light side | graviton pole in `Ĉ₂` (H6/H7) via the unsubtracted `f⁺_{2;2}`; contacts; `[TT]` double-trace content (H9b) | flat dictionary of contacts (`B̂_{k,T}`-type, polynomial in `T`) | 2–3 weeks |
| M6 | finite-`m` correction at fixed `x` (H4 at `O(1/m²)`) | expand the CLS-D crystal integral to next order; the fast part is spin-independent, so the spin dependence is in derivatives of algebraic slow factors | the scalar version must reproduce the CHMRS `1/m²` series and exactskel at `x ≤ 0.6`; `⟨TTTT⟩` position-space data at `x ≤ 0.1` | 2–4 weeks |
| M7 | exact Mellin evaluation of `Ĉ_k` on heavy `⟨TTTT⟩` blocks | spinning Mellin in a structure basis; only needed if M6's accuracy at `x ∈ [0.3, 0.7]` misses the §4.3 targets | — | months |

**Shortest useful sequence:**
* **M2 → M3.** This redefines D-safe on route A′ with the same `n = 1` numbers (now for a defined functional whose flat limit is the
  certified NF02).
* **M6.** This gives H4 on the load-bearing window.
* **M4, M5.** These close H0's validity and light-side clauses.

## 8. Files (`ads4/code/spinfunc/`)

| file | content | output (`cache/`) |
|---|---|---|
| `f1_bulkpoint.py` | sympy: CLS-D crystal ↔ CHMRS frame; leading Regge forms; scalar cross-check; graviton `R_k` | `f1_bulkpoint.txt` |
| `spf_common.py` | `R_X(k,x)` for objects A′ / PsiG / PsiG2 / PsiS; near-gap integrand pieces for any `Δ_g`; state loader (e1 + e7) | — |
| `f2_flat.py`, `f2b_check.py` | flat limit with graviton densities; robustness; twist-gap reach | `f2_flat.txt`, `f2b_check.txt` |
| `f3_ads.py` | near-gap model at `Δ_g = 100/300/1000`; worst-sign tolerances; shares in `x` | `f3_ads_<Δ_g>.txt`, `f3_<Δ_g>_recs.pkl` |
| `f4_etainv.py` | `η`-inversion cost (large-`ν` symbol) | `f4_etainv.txt` |
| `f5_coherent.py` | coherent deviation tolerances | `f5_coherent.txt` |
| `g1_crystal4d.py` | M2: explicit `D = 4` crystal; helicity content of CLS-D kernels; flip-free subspace; channel ratio | `g1_crystal4d.txt` |
| `f6_routeA2.py` | M2/M3: near-gap model of `Ω^{A′} = ∫ψ₂Ĉ₂ + ∫ψ₃Ĉ₃` (stripped), `η`-inversion shift, `δk(0)` | `f6_routeA2_<Δ_g>.txt` (logs `f6_*.log`) |
| `f7_flat_stripped.py` | M2: flat limits under the stripped dictionary (route A′, B′-uplift with CLS-D or project weights) | `f7_flat_stripped.txt` |
| `g2_structures.py` | M2b: `(2,4)` structures at the `D = 4` crystal; flip-free family; `Q²`, `O·Q` | `g2_structures.txt` |
| `h1_etainv_exact.py` | M3b: exact `η`-inversion (Abel/ODE), checks V1–V3, `δk(0)`, tails | `h1_etainv_exact_*.txt`, `h1_<Δ_g>.pkl`, logs `h1_*.log` |
| `h2_neargap_exact.py` | M3b: exact near-gap shift of route A′ | `h2_neargap_exact_*.txt`, logs `h2_*.log` |
| `g3_homogeneity.py` | M2c: crystal value of `(2V_{0,42})²X₂₄/(X₂₀X₄₀)`; consistency of CLS-D's printed photon images vs kernel | `g3_homogeneity.txt` |
| `g4_images.py` | M2c: independent crystal images of `H₂₀H₄₀/(X₂₀X₄₀)`, `H₂₄/X₂₄` vs CLS-D's printed ones | `g4_images.txt` |
| `g5_covariant_Q.py` | M2c: `q = (2H₂₀H₄₀ + V²H₂₄)/(X₂₀X₄₀)` ↦ `−Q` exactly; CLS-D dual coefficients + corrected image reproduce printed `𝒦^{(1)}` | `g5_covariant_Q.txt` |
| `h3_H17_model.py` | M5′b: H17 (`k_AdS − k(0)`) for route A′ in the `nf5` heuristic model; bound factors | `h3_H17_model.txt` |
| `m1_scalar_scaling.py` | M6: `1/m` scaling of the exact scalar remainder at fixed `x` across `Δ_g` (calibration) | `m1_scalar_scaling.txt` |
| `m2_tolerance_c2.py` | M6: H4 tolerance for route A′ in units of `c₂/m²` (sign patterns, windows) | `m2_tolerance_c2_*.txt` |
| `l1_tps.py` | M6: truncated multivariate power series (7 vars, degree 6) and exact Gaussian moments (Laplace engine) | — |
| `l1_crystal.py` | M6 step 1: exact gauge-fixed CLS-D scalar integrand around the crystal; Laplace expansion to `O(1/Δ²)` | — |
| `l2_scalar_validate.py` | M6 step 1: Laplace vs exact scalar skeleton (e1), three gaps, Richardson | `l2_*.pkl`, `l2_scalar_validate_*.txt`, `l2_richardson.txt` |
| `l3_shadow_laplace.py` | M6 step 2 (piece): shadow coefficient from the `x₅` saddle to `O(1/Δ²)`, scalar check | `l3_shadow_laplace.txt` |
| `l4_tto_even.py` | M6 step 2 (piece): parity-even conserved `⟨TT𝒪⟩` structures in `d = 3`, TPS-evaluable, numerical null space | `l4_tto_even.txt` |
| `run_all.sh` | runs all of the above | — |

## 9. Milestone M2: CLS-D in `d = 3` / 4D graviton helicities (`g1_crystal4d.py`, `f7_flat_stripped.py`)

### 9.1 The functionals in `d = 3`

* **Transverse structures.** `SO(d−2) = SO(1)`, so the transverse spin takes only `j′ ∈ {0, 1}`. These are the `1↔3`-even/odd celestial
  structures on `S¹`, which the project already uses (`SUPERCONVERGENT_UPLIFT.md`). `j′ ≥ 2` structures do not exist.
* **Saddle and weight rule.** The crystal locus is `S^{d−3} = S⁰` (two points). The CLS-D Hessian is non-degenerate at `d = 3` for
  `0 < x < 1`. The weight-ratio rule of §3 contains no `d`. **So the `R_k` of §3 hold in `d = 3`** [D given Q].
* **The two functionals route A′ needs:**
  * `Ĉ₂`: `j′ = 0`, `f⁺_{2;2} = v′−u′`; unsubtracted; Regge convergence `J₀ < 2`; marginal for Einstein, which gives the graviton pole.
  * `Ĉ₃`: `j′ = 1`, `f⁻_{3;2} = −1`; unsubtracted; convergence `J₀ < 2`.

  Both are CLS-D's "unsubtracted sum rules for gravitons", whose "action on all double-traces will remain zero" and whose EFT series
  truncates [Q, CLS-D §3.5.2].
* **The `(2,4)` dual structures.** The `j′ = 0` structures are `H₂₀^{2−n}H₄₀^{2−n}H₂₄ⁿ`, `n = 0, 1, 2` (CLS-D §4.3). At the crystal they
  map to `{P², PH, H²}`, with `P ≅ 𝒱₂𝒱₄` (`H₂₀H₄₀`, eq. `H20H40_kernel_toflat`) and `H ≅ 𝒱₂𝒱₄ + 𝒲₂𝒲₄` (`H₂₄` at `D = 4`, eq.
  `H24_kernel_toflat`).

### 9.2 Explicit `D = 4` crystal [DV, `g1`]

* **Kinematics.** The CLS-D shock frame at `D = 4` has real CM kinematics, with `cos θ₁₃ = 1 − 2x²` to `10⁻¹²`.
* **Shock polarisations.** The polarisation `s̄` of particles 1 and 3 is a pure helicity state (self-dual field strength), with
  coefficient `√2 x`. It is the same for every helicity of 2 and 4, so it is a common factor.
* **`H₂₄` is pure flip.** At the crystal, `H ≅ 𝒱𝒱 + 𝒲𝒲` vanishes on every physically helicity-preserving `(2,4)` configuration. This
  agrees with `tr(F₂F₄) = 0` for opposite all-incoming labels and independently checks the helicity bookkeeping.
* **Flip-free subspaces.**
  * Photons: one-dimensional, `𝒱𝒱 − 𝒲𝒲` (∝ CLS-D's `𝒦^{(1)}_{j′=0}` at `D = 4`).
  * Gravitons: two-dimensional, `aP² + bPH + cH²` with `a/4 + b/2 + c = 0`. This mirrors the Regge-frame flip-free family
    `(D₁+D₂)(αD₁+βD₂)` ∋ `K*`, `𝔖_ff`, under `H ↔ D₁+D₂ = z′_in·z′_out` [D-sketch] and `P ↔ D₂` [I].
* **Channel ratio.** CLS-D polarisation-sum the probe polarisations (`sumrule_bulk_polarizationmapall`), so the block sees the
  conjugate of the vector the kernel sees. With that pairing, every flip-free kernel has

  **(opposite-helicity channel)/(same-helicity channel) `= (1−x²)^{−2}` (photons), `(1−x²)^{−4}` (gravitons)**,

  to machine precision at `x = 0.1 … 0.9`.
* **Reading.** The bulk-point density is `|c_s|²P_J + |c_o|²d̃^J_{44}`, with `d̃ = d/(1−x²)⁴`, times a common factor. That is the stripped
  `B^{(1)}_k` (`FLAT_GRAVITON_BOUND.md` (2.33)/(2.35a); reading B̃). It is not `|c_s|²P + |c_o|²d` (B′, reading B), and not
  `P + (1−x²)⁴d`, which is what the unconjugated pairing would give.

### 9.3 Two further, convention-independent arguments

1. **Truncation [D from Q].**
   * CLS-D: the unsubtracted `k = 2, 3` graviton sum rules vanish on all double traces, and their flat EFT series truncate.
   * In the flat basis, the stripped `B^{(1)}_2`, `B^{(1)}_3` are contact-free: `−B₂ = 16πG/q + 2πG|ĝ₃|²q³`, `−B₃ = −2πG|ĝ₃|²q²`.
   * `B′₂`, `B′₃` see `g₄, g₅, g₆, g₆′` (`FLAT_UNSTRIPPED.md` §3, [D]).
   * So the images of `Ĉ₂`, `Ĉ₃` are the stripped rules.
2. **Kernel form [DV, `f1`].** `f⁺_{k;2}|crystal ∝ (2−y)/(1−y)^{(k−2)/2}` and `f⁻_{k;2}|crystal ∝ 1/(1−y)^{(k−3)/2}`. These are exactly the
   `B^{(1)}_k` kernels `(2s−q)/[s(s−q)]^{(k−2)/2}` and `[s(s−q)]^{−(k−3)/2}`. CLS-D's flat sum rules are written for structure-coefficient
   functions, i.e. they are stripped-type by construction (eq. `photon_flatspace_sumrules`, `graviton_flatspace_sumrules_evenk`).

### 9.4 What this does to the existing chain

* **The B′ (reading-B) dictionary is not the dictionary of any sum-rule object here.**
  * It was assumed for the project's moments (`FLAT_UNSTRIPPED.md` item 5, "[I]").
  * Its numerical support (`SKELETON_TEST.md` §8; `SPINNING_TTTT_EXT.md` §3, B over B̃ at `x = 0.1`, 4 % discrimination) concerns
    position-space **Regge moments `Π̂^K`**, which are not sum rules.
  * At `O(x²)` Regge moments and CLS-D functionals are different objects: the same leading Regge kernel, but `(1,3)` polarisations
    summed with kernel weights in `Π̂`, versus fixed to the shock direction in CLS-D. The data do not contradict §9.2–9.3. Neither do
    they test it.
* **The current `Ω^{NF,SC}` is not positive at leading order** (`f7`). Its objects (`C⁺₂` and the `Ψ_k`) are CLS-D-class
  light-ray-commutator sum rules, so their images are stripped (`R_k`-modified for `k = 4, 5, 6`). The B′-uplift weights `φ_k` on
  stripped images give:

  | flat `F/N_s` (`+−`, min over `J ≤ 80`) | `m/M = 1.0` | 1.01 | 1.05 | 1.1 | 1.2 | 1.3 |
  |---|---|---|---|---|---|---|
  | route A′ (`ψ₂B₂ + ψ₃B₃` = NF02) | +0.200 | +0.202 | +0.226 | +0.296 | +0.377 | +0.447 |
  | B′-uplift, CLS-D weights (`R = 1`) | −3479 (`J = 80`) | −17.0 | −0.91 | −0.72 (`J = 8`) | +0.011 | +0.317 |
  | B′-uplift, project `Ψ_k` (`R_k`) | −72286 | −53.9 | −1.04 | −0.68 (`J = 8`) | +0.014 | +0.316 |

  `++` stays positive (min `+0.31`). The `J = 8`, `m = 1.1Λ` state is above `m_min(8) = 1.078Λ` at `Δ_g = 100`, so the failure is
  AdS-reachable. **So "n = 1 for D-safe" in `RESULT_D.md` holds only for the B′ model, not for the functional as defined.** The repair is
  route A′, which is not a B′-uplift.

## 10. Route A′ redefined: `Ω^{A′} = ∫ψ₂Ĉ₂ + ∫ψ₃Ĉ₃` (`f6_routeA2.py`)

**Definition.** `Ω^{A′} := ∫₀¹dp [ψ₂(p) Λ⁴ Ĉ_{2,Λp} + ψ₃(p) ½Λ⁶ Ĉ_{3,Λp}]`, with the NF02 weights and the two unsubtracted CLS-D
functionals of §9.1. `Ĉ₂` carries the `η`-inversion with `j = 1`; `Ĉ₃` has `j = 0` and needs none.

**Properties:**
* **Flat limit = NF02** (stripped basis), so the certificate of `NOFORWARD_FUNCTIONAL.md` applies [Q + D, given §9].
* **No forward terms.**
* **Light side.**
  * Unsubtracted, so zero on all double traces.
  * Flat image contact-free, so (P)/H9b is not needed.
  * No `k ≥ 4` objects, so H15 (`|ĝ₃|²` in `Ψ_{5,6,7}`) and H3 (odd-`k` modelling beyond `k = 3`) drop out.
  * Graviton pole in `Ĉ₂` only, with weight `ψ₂`, the same as the project's `C⁺₂`. Hence `A`, `B_AdS` are the same numbers
    [D-sketch; flat light side by CLS-D's matching expectation].
* **Still open:** endpoint convergence of the two unsubtracted functionals at `J₀ ≤ 2` (I14-type), the exact `η`-inversion, and H4.

**Near-gap model [M].** Per state and node, on the e1/e7 states:
* **Flat part:** stripped `B₂ + B₃` (= NF02, identical to e7's flat part).
* **Curvature:** `ψ₂·2(P_{−½+iν}(cosh β) − J₀) + ψ₃·uσ(P − J₀)`.
* **Remainder:** the exact scalar `ρ₂` and the continued `ρ₃`, with the stripped envelope (`e_k(1−y)^{±2}`; variant: the B′ envelope).
* **Normalisation:** `N_s = |B₂+curv₂| + |B₃+curv₃|`.

| `Δ_g` | `n` (`δ = 0`, `0.03`) | min `F/N_s` (stripped env.) [state] | B′ env. | with `η`-inversion (symbol) | max `|rem|/N_s` | old `Ω^{NF,SC}` B′ model |
|---|---|---|---|---|---|---|
| 100 | 1, 1 | 0.2256 [`+−`, 1.05Λ, `J = 5`] | 0.2266 | 0.2231 | 0.0077 | 0.2267 |
| 300 | 1, 1 | 0.2054 [`+−`, 1.02Λ, `J = 5`] | 0.2056 | 0.2052 | 0.0051 | 0.2057 |
| 1000 | 1, 1 | 0.2021 [`+−`, 1.01Λ, `J = 5`] | 0.2021 | 0.2021 | 0.0030 | 0.2022 |

* **Why the numbers barely move.** The near-gap action is ≈ 98 % the NF02 flat action, which is identical in both. The curvature and
  remainder terms differ, but they are small. So **`n = 1` and the margins carry over to route A′, as a model statement**. They now refer
  to a functional that is defined and whose flat limit is certified.
* **What the model still assumes.** H4 (graviton remainder = scalar remainder × stripped envelope), H1 for `k = 2, 3` under
  CLS-D's `(2,4)` structure, H12′ (the binding state is still `+−`, `J = 5`, parity-odd), and H13′ off the scanned set.

## 11. Milestone M3: `η`-inversion and the headline numbers (`f6_routeA2.py`) [D-sketch + NE]

* **The graviton-pole light term is unchanged exactly** [D-sketch].
  * Do the inversion in `η`-space: `φ̃ = MF[K_φ(η)/η]` (Mehler–Fock; it exists because `η ≥ 1`).
  * Then `∫φ̃(ν) ηP_ν(η) dν = K_φ(η)` exactly, with no contour shift and no analyticity of `φ` needed. The difference equation of §2 needs
    analyticity, which `ψ₂ ∝ p(1−p)^r·poly` lacks at `ν = 0`, so it is not used.
  * The graviton term is a leading-Regge `η`-density (I9). So `∫φ̃ Ĉ^{CLS}_2|_grav = ∫ψ₂ C_2|_grav`, and the integral `L_A` that gives `A`
    and `B_AdS` is unchanged.
* **Shift of `k(0)`.** With the large-`ν` symbol in the interior (reliable there: the `p⁶` weight suppresses small `ν`),
  `δk(0) = −2∫δψ₂p⁶ = −1.50·10⁻³ / −1.64·10⁻⁴ / −1.46·10⁻⁵` (`Δ_g = 100 / 300 / 1000`). NF02 is normalised to `k(0) = 1`, and the check
  reproduces it to `2·10⁻¹¹`.
* **Effect on `A`, `B_AdS` and the bounds.** `A`, `B_AdS` scale as `1/k(0)`: × `1.0015 / 1.00016 / 1.00002`. The bounds scale as
  `k(0)^{−1/2}`: × `1.00075 / 1.00008 / 1.00001`. So **`A = 3022.3` (→ 3026.8 at `Δ_g = 100` if the inversion is included), `B_AdS = −6552.5`
  (→ −6562.4), and `√(t₄²+α_T²) ≤ 2.21·10⁻⁴ / 3.13·10⁻⁶ / 2.89·10⁻⁸` are unchanged at the quoted precision.**
* **Near-gap shift.** Symbol estimate `≤ 0.086 / 0.013 / 0.003 N_s`, dominated by small `ν`, where the symbol is least reliable; `n = 1`
  unchanged. An exact Mehler–Fock evaluation is still to do (M3b; it needs conical functions at `ν ≤ Λ`; est. 2–3 days).
* **Alternative (no inversion).** Use `ψ₂` directly on `Ĉ^{CLS}_2`. Then the flat limit is NF02 exactly with no `δk(0)`. But the Regge
  side becomes `η`-weighted: the curvature term becomes a combination of `P_{−½+i(ν±i)}(cosh β)`, an `O(β)` change, and the graviton
  pole moves to `ν = ±i/2` (so `A` unchanged, `B_AdS` changed at `O(1)`). Not computed; the inverted version is preferred.

### 11.1 Regraded table (supersedes §6 where they differ)

| statement | grade |
|---|---|
| `d = 3`: `j′ ∈ {0,1}`; crystal non-degenerate; `R_k` hold (now relative to stripped `B^{(1)}_k`) | [D given Q] |
| `D = 4` crystal: `s̄` pure helicity; `H₂₄` pure flip; flip-free subspace (photon 1-dim, graviton 2-dim) | [DV] |
| channel ratio `(1−x²)^{−4}` ⇒ CLS-D-class images are the **stripped** `B^{(1)}_k` (reading B̃) | [DV] (given CLS-D's polarisation-sum pairing) + [D from Q] (truncation; kernel form) |
| project kernels `K*`, `𝔖_ff` ↔ crystal flip-free family (`H ↔ D₁+D₂`, `P ↔ D₂`) | [D-sketch] / [I] |
| the B′ dictionary of `RESULT_D` describes no sum-rule object; its support (B over B̃ at `x = 0.1`) is for `Π̂^K` | [D] (scope) |
| current `Ω^{NF,SC}` (B′-uplift of CLS-D-class objects) flat limit negative near the gap (`−0.72 N_s` at `1.1Λ`, `J = 8`), AdS-reachable | [NE] given the above |
| route A′ `= ∫ψ₂Ĉ₂ + ∫ψ₃Ĉ₃`: flat limit = NF02 (certified); contact-free; no `k ≥ 4` | [Q + D] (given §9) |
| route A′ near-gap `n = 1` at `Δ_g = 100/300/1000`, margins `0.226/0.205/0.202` | [NE + M] |
| graviton-pole light term unchanged by Mehler–Fock inversion; `A`, `B_AdS` same numbers | [D-sketch] |
| `δk(0) = −1.5·10⁻³ … −1.5·10⁻⁵`; bounds × `≤ 1.0008` | [D-sketch + NE] (large-`ν` symbol) |

**Hypotheses (proposed; `RESULT_D.md` not edited):**
* **H0:**
  * The object is now specified: two CLS-D functionals, `Ĉ₂` (`j′ = 0`, `f = v′−u′`, `η`-inverted) and `Ĉ₃` (`j′ = 1`, `f = −1`).
  * Its flat limit is the certified NF02 [Q + D].
  * Remaining [I]:
    * the exact `(2,4)` structure in `d = 3`, i.e. which flip-free combination (CLS-D large-`ν` dual-structure algorithm, their App. D);
    * endpoint convergence at `J₀ ≤ 2`;
    * the exact Mehler–Fock inversion;
    * the light-side matching (CLS-D expectation, [Q-expectation]).
  * Proposed grade: [D-sketch].
* **H14 / F1:** resolved for route A′; for the B′-uplift `Ω^{NF,SC}`, F1 is **fatal** (flat limit not positive near the gap).
* **H2:**
  * For route A′: the leading density is the stripped `B^{(1)}_{2,3}` [DV + Q], given CLS-D.
  * The project's B′/reading-B assumption is retired for sum-rule objects.
* **H3, H9b, H15:** not needed for route A′.
* **H4:** unchanged [I]; the transfer model is now with the stripped envelope. Variant spread `≤ 0.001 N_s` in the margin.
* **H17:** gains the `δk(0)` term above (`≤ 1.5·10⁻³`).

### 11.2 Plan update (replaces §7 where they differ)

| # | milestone | estimate |
|---|---|---|
| M2b | CLS-D App. D large-`ν` dual structures in `d = 3` for `⟨TT L[O]⟩`: fix the `(2,4)` combination of `Ĉ₂`, `Ĉ₃` that is flip-free at the crystal for all `x` and reduces to `K*`/`𝔄` at leading Regge order; prove the `H ↔ D₁+D₂`, `P ↔ D₂` map; reproduce the M-table responses | 3–5 days |
| M3b | exact Mehler–Fock inversion for `ψ₂` (conical functions, `ν ≤ Λ`); exact `δk(0)`, near-gap shift; tail of `φ̃` beyond `p = 1` | 2–3 days |
| M4′ | endpoint convergence of the unsubtracted `Ĉ₂`, `Ĉ₃` (`J₀ ≤ 2` only; H15 gone) | ≈ 1 week |
| M5′ | light side of `Ĉ₂`, `Ĉ₃` (graviton pole, `|ĝ₃|²`); double traces vanish by construction | 1–2 weeks |
| M6 | `1/m` correction at the crystal for `Ĉ₂`, `Ĉ₃` (H4), on `x ∈ [0.3, 0.7]` | 2–4 weeks |
| — | update `RESULT_D.md` (coordinator): D-safe := `Ω^{A′}`; B′-uplift statements to be withdrawn | — |

## 12. Milestone M2b: the `(2,4)` structures of `Ĉ₂`, `Ĉ₃` (`g2_structures.py`)

**Exact facts at the `D = 4` crystal [DV, `x = 0.05 … 0.95`, machine precision; CLS-D polarisation-sum pairing].** The photon-level images
are `P := 𝒱₂𝒱₄` (`H₂₀H₄₀`) and `H := 𝒱₂𝒱₄ + 𝒲₂𝒲₄` (`H₂₄`).
1. `H` vanishes on preserved helicities, `H_flip = 2P_flip`, and `P_flip² = P_same P_opp`.
2. Hence, in the graviton `j′ = 0` span `aP² + bPH + cH²`, the flip-free condition is **`a + 2b + 4c = 0`, independent of `x`**. This gives
   a 2-dimensional family, like the Regge-frame flip-free family `(D₁+D₂)(αD₁+βD₂) ∋ K*, 𝔖_ff`.
3. Every member of the family has preserved weight `a·P²_pres`. **So the bulk-point density is `a ×` (stripped `B^{(1)}_2` density) for
   any flip-free choice**; the choice only fixes the normalisation `a`. Members with `a = 0` are null at the bulk point.
4. Even sector: `Q := H − 2P = 𝒲𝒲 − 𝒱𝒱` is exactly flip-free (photon), and so is `Q²` (graviton; `a = 4`, `b = −4`, `c = 1`). Channel ratio
   (opposite/same) `= (1−x²)^{−4}` exactly.
5. Odd sector: `O := 𝒲₂𝒱₄ − 𝒱₂𝒲₄` is flip-free, with channel ratio `−(1−x²)^{−2}` (photon). The graviton `O·Q` is flip-free with
   `−(1−x²)^{−4}`, the `B^{(1)}_3` sign structure. The symmetric `𝒲𝒱 + 𝒱𝒲` is pure flip. CLS-D's single `𝒲₂𝒱₄` has helicity-odd
   flips that cancel in the `2↔4`-antisymmetric combination.

**Covariant realisation.**
* **Photons [Q + DV].** CLS-D's photon dual structure `(1)` at `j′ = 0` maps at the crystal to `(K_ν/2)[(𝒱𝒱+𝒲𝒲) − (D−2)𝒱𝒱]`
  (eq. `JJO_kernels_flatspace`), which is `(K_ν/2)Q` at `D = 4`, **for every `x`**, with `x`-independent coefficients. So the photon
  flip-free kernel is a fixed covariant CFT structure.
* **Gravitons [D-sketch].** The graviton `j′ = 0` structures `H₂₀^{2−n}H₄₀^{2−n}H₂₄ⁿ` carry the same `V`-power bookkeeping per factor.
  Hence `Q²` and `O·Q` are images of fixed covariant combinations. CLS-D do not print the graviton kernels; an explicit large-`ν` dual
  computation (their App. D algorithm) would make this [D] (est. 2–3 days).
* **This is now the definition of `(a₊)`, `(a₋)` in §R.**

**Relation to the project's `K*`, `𝔖_ff`, `𝔄` [D-sketch / I].**
* **What matches.** Both families are 2-dimensional. The natural identification `D₁+D₂ = z′_in·z′_out ↔ Q` pairs the two exactly
  flip-free photon structures in their regimes: Lemma L flip of `D₁+D₂` is `0` for all `ν`; `Q`'s flip is `0` for all `x`. They agree in the
  overlap (`x → 0` ↔ `ν → ∞`), where both are pure-preserved and helicity-symmetric. Under it, the Regge family maps into the crystal family
  for any image of `D₂`.
* **What is not shown.** A proof needs the Lemma-L transform of the covariant `(a₊)`; it was not done.
* **Why this does not block route A′.** It is not needed to define `Ω^{A′}`. It is needed only to carry over the existing `K*`/`𝔄`
  heavy-formula evidence (H1) at finite `β`.

## 13. Milestone M3b: exact `η`-inversion (`h1_etainv_exact.py`, `h2_neargap_exact.py`)

**Method [D].**
* Mehler–Dirichlet: `P_{−½+iν}(cosh ξ) = (√2/π)∫₀^ξ cos(νt) dt/√(cosh ξ − cosh t)`. So the Regge `η`-kernel is an Abel transform of
  `Φ(t) = ∫φ(ν)cos(νt)dν`.
* Integrating by parts, `cosh ξ · A[Φ̃] = A[cosh t Φ̃ + ½ sinh t H̃]` with `H̃′ = Φ̃`. So `ηK̃ = K` is a first-order linear ODE with the
  explicit solution of §R item 5.
* No analyticity of `φ` in `ν` is needed. The `ν`-difference equation of §2 would need it, and `ψ₂ ∝ p` lacks it at `ν = 0`.

**Numerics.** `Φ` in closed form (terminating integration by parts, `ψ₂` polynomial; Gauss–Legendre for `Λt < 400`); the ODE integral by
integration by parts plus cumulative Simpson; `Δt = π/(16Λ)`, `t ≤ 60`.

**Checks:**

| check | `Δ_g = 100` | 300 | 1000 |
|---|---|---|---|
| V1: interior `δφ/φ` vs large-`ν` symbol (`p = 0.2…0.8`) | agree to 0.3–1.2 % of `δφ` | `≤ 0.15 %` | `≤ 10⁻⁴` |
| V2: Abel identity `cosh ξ K̃ = K` | `4.7·10⁻¹⁰` (`ξ = 0.05`) … `5·10⁻⁴` (`ξ = 5`, check-quadrature-limited) | `7.7·10⁻¹⁰` … | `7·10⁻⁹` … |
| **V3: graviton-pole identity** `∫ψ̃₂ g̃ = ∫ψ₂ g` (all `ν`, incl. `ν → 0`) | **`5.3·10⁻¹⁴`** (without inversion: `1.5·10⁻²`) | **`2.6·10⁻¹⁴`** (`8.5·10⁻³`) | **`1.9·10⁻¹⁴`** (`5.7·10⁻³`) |
| max `|ψ̃₂ − ψ₂|/max ψ₂` on `[0,1]` (attained at small `ν`) | `7.6·10⁻³` | `2.5·10⁻³` | `7.3·10⁻⁴` |
| tail mass `∫_{p>1}|ψ̃₂|/∫|ψ₂|` | `2.0·10⁻⁹` | `4.7·10⁻¹²` | `2.4·10⁻¹⁴` |

V3 is the end-to-end statement that **the graviton-pole light term, and hence `L_A`, `A`, `B_AdS`, are exactly unchanged** by the inversion.
It assumes the graviton term is a leading-Regge `η`-density (I9).

**Results:**

| `Δ_g` | `δk(0)` exact (symbol) | `A`, `B_AdS` factor `1/k(0)` | bound factor `k(0)^{−½}` | near-gap `n`; min `F/N_s` without / with exact inversion | max `|δF|/N_s` (median) |
|---|---|---|---|---|---|
| 100 | `−1.5074·10⁻³` (`−1.5036·10⁻³`) | 1.00151 | 1.00075 | 1; 0.2256 / **0.2231** [`+−`, 1.05Λ, `J = 5`] | `0.089` (`2.5·10⁻³`) |
| 300 | `−1.6377·10⁻⁴` (`−1.6372·10⁻⁴`) | 1.00016 | 1.00008 | 1; 0.2054 / **0.2052** | `0.013` (`4·10⁻⁴`) |
| 1000 | `−1.4632·10⁻⁵` (`−1.4632·10⁻⁵`) | 1.00001 | 1.00001 | 1; 0.2021 / **0.2021** | `0.003` (`6·10⁻⁵`) |

* **Headline numbers, route A′ (in the exact-skeleton model, same hypotheses):**
  * `A = 3022.3 → 3026.9` and `B_AdS = −6552.5 → −6562.4` at `Δ_g = 100`; both unchanged to `2·10⁻⁴` at `Δ_g ≥ 300`.
  * Bounds `√(t₄²+α_T²) ≤ 2.21·10⁻⁴ / 3.13·10⁻⁶ / 2.89·10⁻⁸`, unchanged at the quoted precision (factor `≤ 1.00075`).
* **Grades.** [D] for the inversion. [DV] for its numerics. [NE + M] for the near-gap numbers. The tails at `p > 1` (`≤ 2·10⁻⁹`) are not
  in the near-gap sum [I].

## 14. Milestone M4′: convergence and the light side of `Ĉ₂`, `Ĉ₃` (as far as clean) [Q + D-sketch]

1. **Regge corner [Q].** CLS-D §3.5.2: `f⁺` changes the effective Regge spin by `k − 2J_e + 1`, `f⁻` by `k − 2J_e`. So `Ĉ₂`
   (`f = v′−u′`) and `Ĉ₃` (`f = −1`) converge for `J₀ < 2`. The full non-perturbative correlator has `J₀ ≤ 1`, so both are valid sum rules
   there. Einstein exchange (`J₀ = 2`) is marginal: the graviton pole in `Ĉ₂` (as for CHMRS `C₂`).
2. **Elsewhere (light-cone and single-ray endpoints) [D-sketch].** Both weights are polynomials in `(u′, v′)`, with no poles. So away from
   the Regge corner they cannot worsen the convergence of the bare commutator `⟨T[L[T],L[T]]T⟩`; that convergence is quoted from KKSDZ /
   CLS-D and not re-derived [I-minor]. The objection of `SUPERCONVERGENT_UPLIFT.md` M2 (endpoints for `3 ≤ J₀ < k`) concerned
   `k = 5, 6, 7` and **does not arise**: route A′ has no `k ≥ 4`.
3. **Light side at tree level [D-sketch].**
   * The `[TT]` double twists `T∂^ℓ□ⁿT` have `τ = Δ − J ∈ {2, 4, 6, …}`. The unsubtracted functionals carry the double-commutator factor
     `∝ sin²(π(τ − 2)/2)` (CLS-D: "action on all double-traces will remain zero").
   * MFT (identity + double traces) is a physical correlator, so the identity contributes `0` as well.
   * So **the tree-level light side of `Ω^{A′}` is exactly `Ĉ₂[G_T] + Ĉ₃[G_T]`**: the single stress-tensor block, with its full
     `⟨TTT⟩ = (c_T, t₄, α_T)` dependence.
   * Contact interactions only shift double-trace data and drop out. So H9b, H10 (contact power counting) and H15 are not needed at tree
     level.
4. **What remains:**
   * the marginal (`J = 2`) evaluation of `Ĉ₂[G_T]` (graviton pole, `|ĝ₃|²` coefficient `k_AdS`) and `Ĉ₃[G_T]`: one explicit
     spinning-block evaluation, which is M5′;
   * the swappability of the heavy sum with the functional (I14, as in CHMRS) [I].

### 14.1 Regraded table (M2b/M3b/M4′; supersedes §11.1 where they differ)

| statement | grade |
|---|---|
| crystal: `a + 2b + 4c = 0` flip-free family (`x`-independent); universal preserved weight; `Q²`, `O·Q` flip-free with exact stripped ratios `±(1−x²)^{−4}` | [DV] |
| photon `Q` = CLS-D dual structure `(1)` at `D = 4` for all `x` | [Q + DV] |
| graviton `Q²`, `O·Q` realised by fixed covariant `j′ = 0, 1` structures | [D-sketch] |
| identification with the project's `K*`/`𝔖_ff`/`𝔄` at finite `β` (needed only to import H1 evidence) | [D-sketch] / [I] |
| exact `η`-inversion formula | [D] |
| graviton-pole term (`A`, `B_AdS`) invariant under the inversion | [D-sketch] (I9) + [DV] (V3, `≤ 5·10⁻¹⁴`) |
| `δk(0)` exact `−1.507·10⁻³ / −1.638·10⁻⁴ / −1.463·10⁻⁵`; bounds × `≤ 1.00075` | [DV] (flat-limit light side) |
| route A′ near-gap `n = 1`, min `F/N_s 0.2231 / 0.2052 / 0.2021` with exact inversion | [NE + M] |
| `Ĉ₂`, `Ĉ₃` converge for `J₀ < 2`; polynomial weights add no endpoint divergence | [Q] / [D-sketch] |
| tree-level light side = the `T` block only (double traces, contacts, identity drop out) | [D-sketch] |

**Hypotheses (proposed):**
* **H0:**
  * The object is fully specified (§R). Its flat limit is certified NF02, and its inversion is exact.
  * Remaining: the graviton covariant structure ([D-sketch]), I14, and the `T`-block light side (M5′).
* **H9b, H10, H15:** not needed at tree level under route A′.
* **H6, H7, H16:** reduce to `Ĉ₂[G_T]`, `Ĉ₃[G_T]`.
* **H17:** gains `δk(0)` (exact).

**Next:**
| step | content | estimate |
|---|---|---|
| M5′ | `Ĉ₂[G_T]`, `Ĉ₃[G_T]`, i.e. `k_AdS`, the pole residue, and the `t₄`/`α_T` linear terms, via the CLS-D crystal plus the Regge-side `T`-exchange density | 1–2 weeks |
| M2c | explicit graviton dual structures (CLS-D App. D) | 2–3 days |
| M2d | Lemma-L transform of `(a₊)`, `(a₋)` (imports H1) | 3–5 days |
| M6 | H4 | 2–4 weeks |

## 15. Milestone M5′: the light side of `Ω^{A′}` and the headline constants [Q + D-sketch]

**Sources.** CHMRS (arXiv:2106.10274) LaTeX: `sec_light.tex` (eqs. `eftformula`, `MGravityRegge`, `BknuGravity`,
`BknuFlatSpaceLimit`, `gravitonpropagator`), `app_omega_light.tex`. CLS-D eqs. `Blockvspartialwave`, `heavyaction_definition_generalrep`
(the spinning double-commutator factor).

### 15.1 What the light side is (tree level, `O(1/c_T)`)

1. **Double traces drop out.**
   * CLS-D: for a symmetric-traceless exchange the double commutator gives `2 sin(π(τ − Δ₁−J₁−Δ₂−J₂)/2)·sin(π(τ − Δ₃−J₃−Δ₄−J₄)/2)` [Q]. For
     external `T` (`Δ+J = 5`) this is `2 sin²(π(τ−10)/2)`: a double zero at every even `τ`, and `≠ 0` at `τ = 1` (`T` itself).
   * The `[TT]` families `T∂^ℓ□ⁿT` have `τ ∈ {2, 4, 6, …} + O(1/c_T)` (spin components lower `τ = 6 + 2n` by `0…4`).
   * The weights `f₂ = v′−u′` and `f₃ = −1` are **polynomials** in `(u′, v′)`. Multiplying by `u′` shifts the effective `s`-channel twist by `+2`,
     so the even-twist zeros are kept. (A subtraction factor with `1/(u′v′)` poles would cancel zeros on the first trajectories; that is
     the CHMRS `C_{k,ν}` situation.)
   * So each `[TT]` contributes `∝ γ² = O(1/c_T²)` [D given Q].
   * **Caveat: parity-odd double traces.** In `d = 3`, `[TT]` also contains parity-odd double traces (`ε`-contractions, e.g.
     `ε T∂T`), with **odd** `τ`. Their double zeros need the parity-odd `⟨TT𝒪⟩` structures, which carry an odd power of `|x₁₂|`, to
     shift the commutator phase by `π/2`, giving `sin(π(τ±1−10)/2)`. This is expected but not shown: CLS-D derive the factor only for the
     parity-even structures (their footnote: "interesting to determine `τ̃_ρ` for general `ρ`") [D-sketch]. It is tied to H12′.
2. **Identity drops out.** `L[T]|0⟩ = 0`, so `⟨T₄T₂⟩⟨[L[T₃], L[T₁]]⟩ = 0`, and the cross terms vanish the same way [D]. (Equivalently, MFT is a
   physical correlator and Ψ annihilates it.)
3. **Hence, at tree level, `−Ω^{A′}|_light = Ĉ₂[p_T G^s_T] + Ĉ₃[p_T G^s_T]`.** This is the single stress-tensor block, quadratic in `⟨TTT⟩ =
   (c_T, t₄, α_T)`.
   * Contact interactions change only double-trace data, so they drop out. **H9b, H10 (contact power counting) and H15 are not needed.**
   * Loop-level (`1/c_T²`) double-trace effects are outside the tree-level statement, as everywhere in `RESULT_D`.

### 15.2 Evaluating it: CHMRS's residue formula

* **The formula [Q].** If a functional is a pure spin-`k` Regge moment up to moments of Regge spin `> n`, with `n` the Regge spin of the
  tree-level EFT correlator `𝒢^EFT`, then `−ω|_light = (2πi/4)·[𝒢^EFT_+(r,η)]_{r^{1−k}}` exactly: the `r^{1−k}` Laurent coefficient at the
  `u`-channel Regge limit. The proof (App. `omega_light`) uses only contour deformation and crossing antisymmetry of `ω`.
* **The scalar example [Q].** For `C_{2,ν}` on graviton exchange the result is exactly `−C_{2,ν}|_grav = 8πG/(ν² + (d/2)²)`, checked
  numerically by CHMRS against the Mellin form.
* **Applicability to `Ĉ₂`, `Ĉ₃` [D-sketch].**
  * After the exact inversion, `Ĉ₂`'s leading kernel is the pure `Π_{2,η}` with weight `K_φ(η)`.
  * The next terms are `O(r²)`: `f₂ = (η/2)r^{−1}(1 + r²)`, and the light-ray kernel expands in even steps. So `Ĉ₂ = Π̂₂ + O(Π̂₄)`, and likewise
    `Ĉ₃ = Π̂₃ + O(Π̂₅)` (`f₃ = −1`).
  * The CHMRS proof carries over to the kernel-contracted spinning correlator, assuming its Laurent expansion near `r = 0` is single-valued
    with integer powers (tree level) [I-minor].
* **Einstein part (`n = 2`).** The formula applies (`k + 1 = 3 > 2`):
  * `−Ĉ₂|_E = (2πi/4)[𝒢^E]_{r^{−1}}`, i.e. the leading Regge density of Einstein exchange, which is the transverse graviton propagator
    `h_{d−1}`.
  * The Einstein shock is helicity-preserving and helicity-independent [Q, KKSDZ/CEMZ]. So under the normalised flip-free kernel `Q⊗Q` it
    gives, with CHMRS eq. `gravitonpropagator` (MF transform of `h_{d−1}` ↔ `1/(ν² + (d/2)²)`),
    **`−Ĉ₂|_E = 16πG/(ν² + 9/4)` exactly.** The normalisation `16πG` is fixed by the flat limit, `−B₂ ∋ 16πG/p²`.
  * `−Ĉ₃|_E = [𝒢^E]_{r^{−2}} = 0`: the spin-2 Regge term is `s↔t` symmetric at leading order. This agrees with `OMEGA_PRIME.md`'s
    pointwise statement for `A₃ = f⁻_{3;2}`.
  * **So I9 is upgraded.** It is no longer "light side ≈ its leading-Regge term" but "light side = a Regge residue, exactly", conditional
    only on the helicity-independent Einstein Regge density (H7) and on applicability above.
* **Non-Einstein part (`t₄`, `α_T`; `|ĝ₃|²`).**
  * The flat MHV term `2πG|ĝ₃|²su/t` makes `𝓜` grow like `s⁴`, so `n = 4 > k + 1`. The residue formula then also receives the `Π̂₄`
    (`Π̂₅`) admixtures of `Ĉ₂` (`Ĉ₃`).
  * At large `ν` the result is the flat one (CHMRS eq. `BknuFlatSpaceLimit`: `Res_{s=0}[kernel × 𝓜_flat](1 + O(1/ν²))`):
    `−Ĉ₂|_{ĝ₃} → 2πG|ĝ₃|²(ν/Λ)⁶·Λ^{−2}…`, `−Ĉ₃|_{ĝ₃} → −2πG|ĝ₃|²(…)⁴`, which is `k(0)`.
  * The `O(1/ν²)` AdS corrections are H17 (`k_AdS − k(0)`). They now need to be recomputed for `Ĉ₂`, `Ĉ₃`: the `Π̂₄` coefficient of `Ĉ₂` and
    the `r^{−3}` Laurent coefficient of the `R³`-type `T`-exchange (M5′b).
  * A term linear in `t₄` is absent from the flat MHV `f`. At `O(1/Λ²)` it is not excluded [I], as in `RESULT_D`.

### 15.3 The headline constants

* **`A`.** `A = 16x₀/k(0)`: the coefficient of `log Λ`, i.e. the large-`ν` tail `16πG/ν²` (flat) times `ψ₂′(0)`, divided by `k(0)`.
  **Exact** given §15.2's Einstein statement, up to `k(0)`.
* **`B_AdS`.** `B_AdS = 16L_fin/k(0)`: the finite part of `∫ψ₂(t)dt/(t² + 9/(4L²))`, which uses the exact `ν = O(1)` shape `1/(ν² + 9/4)`.
  **Exact** given §15.2, up to `k(0)`.
  * Sensitivity: a different shape `1/(ν² + c²)` would shift `B_AdS` by `−A log(c/(3/2))`, e.g. `≈ 300` for `δ log c = 0.1`.
  * This is why the CHMRS-type exactness matters. The shape is now [D given H7], not an approximation.
* **`k(0)`.** The flat value plus `O(1/Λ²)` (H17, previously `≈ 4 %` at `Δ_g = 100` for the B′ objects, not yet redone for `Ĉ₂`, `Ĉ₃`), plus the
  `η`-inversion `δk(0)` of §13.
* **Verdict.** **`A = 3022.3`, `B_AdS = −6552.5` and the bounds survive under route A′.**
  * Their graviton-pole content is exact given (a) the CHMRS residue formula for `Ĉ₂`, `Ĉ₃` [Q + D-sketch] and (b) the helicity-independent
    Einstein Regge density (H7) [Q + D-sketch].
  * Their `|ĝ₃|²` normalisation carries the same `O(1/Λ²)` (H17) caveat as before.

### 15.4 I14 (swappability)

Both the sum rule and the residue formula need the heavy sum `Σ p_{Δ,J}Ĉ[G_{Δ,J}]` to converge and to commute with the functional
integral. This is the I14 premise, as in CHMRS. For route A′ it is needed only for two functionals that converge for `J₀ < 2` (§14).
**It is not proved here [I].** CHMRS make the same assumption for the scalar `C₂`.

## 16. Milestone M2c: graviton `(2,4)` structures — status [D-sketch], with an open consistency check (`g3_homogeneity.py`)

> **Resolved in §17:** the discrepancy is a typo in CLS-D's printed `H₂₀H₄₀` crystal image. Their kernel is right, and there is no flip leakage.

**The argument.**
* Products of structures evaluate to products at the crystal, so the graviton `j′ = 0` basis
  `(−2V_{0,42})^{J_F−4+2n}H₂₀^{2−n}H₄₀^{2−n}H₂₄ⁿ` maps to `(common)·(c_P P)^{2−n}(c_H H)ⁿ`.
* `Q² = (H − 2P)²` is therefore the image of a **fixed** covariant combination iff `c_H/c_P` is `x`-independent.
* CLS-D's printed photon kernel `𝒦^{(1)}_{j′=0} ∝ 𝒱𝒱 + 𝒲𝒲 − (D−2)𝒱𝒱` (constant coefficients, all `x`) implies exactly that, because it is a fixed
  dual combination of the photon structures. Taking that formula as given, **`Q²` (and by the same argument `O·Q`) is the image of a fixed
  covariant graviton structure [D given Q]**.

**Independent check (`g3`).**
* At the crystal, `(2V_{0,42})²X₂₄/(X₂₀X₄₀) = 1` exactly for all `x` [DV]. That is consistent with CLS-D's fast factor
  (`bulkpoint_integrand_factor`).
* But combining it with CLS-D's separately printed images `H₂₀H₄₀/(X₂₀X₄₀) = 𝒱𝒱/(4ν²(m²−ν²))` and
  `H₂₄/X₂₄ = −(𝒱𝒱+𝒲𝒲)/(ν²(m²−ν²)²)` (`D = 4`), and their printed dual coefficients `16(d−1)/(d−2)`, `4/(d−2)`, gives
  `𝒦^{(1)} ∝ (d−1)(m²−ν²)𝒱𝒱 − (𝒱𝒱 + 𝒲𝒲)`.
* This differs from their printed `𝒦^{(1)}` by a factor `(1−x²)` on `𝒱𝒱`. Both readings agree at `x → 0` (Regge).
* If the printed images (rather than the printed kernel) were right, no fixed photon structure would be flip-free for all `x`. The
  flip-free direction would rotate with `x`, and route A′'s `(2,4)` structure would carry `O(x²)` helicity-flip interference at the bulk
  point.

**Status.** M2c is **not closed**. It needs an independent large-`ν` computation of the photon dual structures (App. D:
continuous-spin basis, light-transform replacement rule `V → −(V + H/V_{0,42})`, the pairing at the standard configuration), which
settles which printed formula holds. Estimate 2–3 days. Until then: `Q²`, `O·Q` realisable at all `x` [D given CLS-D eq.
`JJO_kernels_flatspace`]; possible `O(x²)` flip leakage [I]; the stripped-dictionary conclusion of §9 is unaffected (it holds for any
flip-free kernel, and the truncation argument of §9.3 is independent).

### 16.1 Regraded table (M5′, M2c; supersedes §14.1 where they differ)

| statement | grade |
|---|---|
| `Ĉ₂`, `Ĉ₃` vanish on all parity-even `[TT]` (double zeros at even `τ`; polynomial `f` keeps them); identity drops out | [D given Q] (CLS-D commutator factor) |
| same for parity-odd `[TT]` (odd `τ`; needs the `π/2` phase of odd structures) | [D-sketch] |
| tree-level light side = `Ĉ₂[p_T G_T] + Ĉ₃[p_T G_T]`; H9b, H10, H15 not needed | [D given Q] |
| CHMRS residue formula applies to `Ĉ₂` (`Π̂₂ + O(Π̂₄)`), `Ĉ₃` (`Π̂₃ + O(Π̂₅)`) | [Q] (formula) + [D-sketch] (spinning applicability) |
| `−Ĉ₂|_E = 16πG/(ν² + 9/4)` exactly; `Ĉ₃|_E = 0` | [D given H7 + the above] (was I9, a leading-Regge assumption) |
| `A`, `B_AdS` graviton-pole content exact; `A = 3022.3`, `B_AdS = −6552.5` survive | [D given H7, I14, residue-formula applicability] |
| `|ĝ₃|²` coefficient = flat `k(0)` + `O(1/Λ²)` (H17, to be redone for `Ĉ₂`, `Ĉ₃`) + `δk(0)` | [Q-analogue + D-sketch] |
| graviton `Q²`, `O·Q` realised by fixed covariant structures | [D given CLS-D's printed `𝒦`]; open check (§16) |
| I14 swappability | [I] |

## 17. Milestone M2c resolved: CLS-D's printed `H₂₀H₄₀` image has a typo; the kernel is right; no flip leakage (`g4_images.py`, `g5_covariant_Q.py`) [DV]

**Independent evaluation at the `D = 4` crystal.** Built from CLS-D's own definitions: embedding `X`, `Z`, `H_{ij}`, `V_{i,jk}`, `X_{ij}`; the
scattering-crystal points; the polarisation maps `polarvec_125/346_CFTtoflat`. Results for all `x = 0.1…0.97` and both `S⁰` points:
* `H₂₄/X₂₄` = CLS-D's printed image `−(𝒱𝒱+𝒲𝒲)/(ν²(m²−ν²)²)` **exactly**.
* `H₂₀H₄₀/(X₂₀X₄₀) = 𝒱𝒱/(4ν²(m²−ν²)²)`. CLS-D print `𝒱𝒱/(4ν²(m²−ν²))`, **a typo (missing `1/(m²−ν²)`)**: the ratio of the two is
  `1/(1−x²)` at every `x`.
* `(2V_{0,42})²X₂₄/(X₂₀X₄₀) = 1` exactly.
* With the corrected image, CLS-D's printed dual coefficients `16(d−1)/(d−2)`, `4/(d−2)` (eq. `JJO_kernel_dual`) reproduce their printed
  kernel `𝒦^{(1)}_{j′=0} ∝ 𝒱𝒱 + 𝒲𝒲 − (D−2)𝒱𝒱` with `x`-independent normalisation. **So the printed kernel is right.** The
  `(1−x²)` discrepancy of §16 came entirely from the image typo.

**Explicit covariant flip-free structure.**
* `q := (2H₂₀H₄₀ + V_{0,42}²H₂₄)/(X₂₀X₄₀)` maps to `−Q/(4ν²(m²−ν²)²)` (`Q = 𝒲𝒲 − 𝒱𝒱`) **exactly** (`≤ 10⁻¹⁵`, all four entries) at the
  `S⁰` point of the `g1` shock frame. The other point is its `x₃`-reflection, which also flips the shock helicity of 1 and 3, so kernel and
  channel labels flip together.
* So the photon flip-free structure is a fixed covariant structure, flip-free at the bulk point for all `x`.
* The graviton `(2,4)` structure of `Ĉ₂` is **`q²`**, i.e. `4S₀ + S₁ + S₂/16` in CLS-D's basis
  `S_n = (−2V_{0,42})^{J_F−4+2n}H₂₀^{2−n}H₄₀^{2−n}H₂₄ⁿ`. Products of structures evaluate to products, so `q² ↦ Q²/(…)`, which is flip-free
  with the exact stripped channel ratio `(1−x²)^{−4}` (§12) [DV].
* Odd sector: `O·Q` with `O` the `2↔4`-antisymmetrised CLS-D `𝒦_{j′=1}`. This is a single-structure kernel, so there is no
  relative-normalisation issue [D given Q].

**Consequences.**
* The leakage scenario of §16 (`O(x²)` helicity-flip interference) **does not occur**, so no near-gap quantification is needed.
* The (2,4) part of H0 is now explicit and checked: [DV] (photon), [D] (graviton even), [D given Q] (odd).
* The relation to the project's Regge-frame `K*`/`𝔖_ff` (§12) is still only needed to import the H1 evidence: [D-sketch]/[I].

## 18. Milestone M5′b: H17 (`k_AdS − k(0)`) for `Ĉ₂`, `Ĉ₃` — structure and model estimate [D-sketch + M]

**What the correction is.**
* The `|ĝ₃|²` light side has Regge spin `n = 4`, so by CHMRS's residue formula (§15.2):
  * `−Ĉ₂|_{ĝ₃} = (2πi/4)([𝒢^{R³}]_{r^{−1}} + c₄^{(2)} ⊛ [𝒢^{R³}]_{r^{−3}})`;
  * `−Ĉ₃|_{ĝ₃} = (…)([…]_{r^{−2}} + c₅^{(3)} ⊛ […]_{r^{−4}})`.
* Here `c₄^{(2)}`, `c₅^{(3)}` are the `Π̂₄` / `Π̂₅` admixtures of `Ĉ₂`, `Ĉ₃`. They come from the subleading terms of the subtraction factor
  (`f₂ = (η/2)r^{−1}(1 + r²)`: relative `+r²`) and of the spin-3 light-ray kernel (not derived for spinning externals).
* At large `ν` the sum is the flat `Res_{s=0}` (`−B₂ ∋ 2πG|ĝ₃|²p⁶`, `−B₃ ∋ −2πG|ĝ₃|²p⁴`). The AdS terms are `O(1/ν²)` relative:
  `−Ĉ₂|_{ĝ₃} ∝ ν⁶ + a₄ν⁴ + a₂ν² + a₀`, by analogy with CHMRS's exact scalar results such as `−C₂|_{g₃} = ν² − (2Δ_φ − (d−4)/2)(2Δ_φ + (d−4)/2)`.
* The exact `a_j` need the subleading Regge Laurent coefficients of the spinning `R³` (`t₄`, `α_T`) `T`-exchange in AdS₄, plus `c₄`, `c₅`.
  **Not computed** (estimate 1–2 weeks: spinning exchange Witten diagram in the Regge limit to `O(r²)`, or its Mellin/impact-parameter form).

**Model estimate (`h3`; the same heuristic as `omegau/nf5_light.py`, `|a| = D × 36` per power lowered) [M]:**

| `Δ_g` | route A′ `|δk|` | (old B′-uplift, same model) | `A`, `B_AdS` factor | bound factor |
|---|---|---|---|---|
| 100 | 0.044 | (0.042) | [0.958, 1.046] | [0.979, 1.023] |
| 300 | 0.0048 | (0.0046) | [0.995, 1.005] | [0.998, 1.002] |
| 1000 | 0.00043 | (0.00041) | [0.9996, 1.0004] | [0.9998, 1.0002] |

* Route A′ is the same size as the old functional in this model. `|ĝ₃|²` enters only through `ψ₂p⁶` and `ψ₃p⁴`, rather than through `B′₂…B′₄`.
* With the exact inversion `δk(0)` (§13) added, the headline bounds `2.21·10⁻⁴ / 3.13·10⁻⁶ / 2.89·10⁻⁸` carry an H17 uncertainty of
  `≈ ±2 % / ±0.2 % / ±0.02 %` [M]. That is unchanged from `RESULT_D`'s budget; the sign is unknown.
* **H17 stays [D-sketch + M].** It becomes [D] once the `a_j` are computed.

### 18.1 Regraded table (supersedes §16.1 where they differ)

| statement | grade |
|---|---|
| CLS-D printed `H₂₀H₄₀` crystal image off by `1/(m²−ν²)`; printed `H₂₄` image, dual coefficients and `𝒦^{(1)}` consistent once corrected | [DV] |
| photon flip-free `(2,4)` structure = fixed covariant `q = (2H₂₀H₄₀ + V²H₂₄)/(X₂₀X₄₀)`, image `−Q` exactly, all `x` | [DV] |
| graviton `(2,4)` structure of `Ĉ₂` = `q²` (flip-free, stripped ratio `(1−x²)^{−4}`); no flip leakage | [DV] (crystal) / [D] (covariance) |
| odd `(2,4)` structure of `Ĉ₃` = `O·q`-type | [D given Q] |
| `Ĉ₂`, `Ĉ₃` vanish on parity-even `[TT]`; parity-odd `[TT]` | [D given Q]; [D-sketch] |
| `−Ĉ₂|_E = 16πG/(ν² + 9/4)`, `Ĉ₃|_E = 0` exactly; `A = 3022.3`, `B_AdS = −6552.5` survive (graviton-pole content) | [D given H7, I14, residue-formula applicability] |
| `|ĝ₃|²` coefficient: `k(0)`; AdS correction H17 `≈ 4.4 % / 0.48 % / 0.043 %` (bounds `±2 / 0.2 / 0.02 %`) | [D-sketch + M] |
| I14 | [I] |

## 19. Milestone M6 (H4 for route A′): scalar calibration, graviton structure, explicit targets, plan

**Scope of this round.** The full graviton fluctuation expansion around the CLS-D crystal was **not** done. It needs, in about 10 integration
variables around a complex saddle:
* the spin-2 dual structure `q²` and the celestial structures;
* the heavy `⟨TT𝒪_{Δ,J}⟩` conserved structures beyond leading large-`Δ` order;
* the polarisation sums and the shadow-coefficient saddle.

Delivered below: the scalar calibration (from the exact skeleton), the structure of the graviton correction, the H4 targets in calibrated
units, and a plan.

### 19.1 Scalar calibration from the exact skeleton (`m1_scalar_scaling.py`) [NE, exact data]

* **Method.** The 64 `p`-nodes are shared by `Δ_g = 100, 300, 1000`, so the exact remainder `ρ_k(x; m, J)` (Mellin/Barnes, `Δ_φ = 3`) can be
  compared at identical `x` and `(m/Λ, J)` for `m ≈ 100, 300, 1000`.
* **Result: on the load-bearing window `x ∈ [0.3, 0.7]`, `ρ_k = c₂(x, J, k)/m² + c₄/m⁴ + …` with no odd terms.**
  * The local exponent is `1.97–1.98` (`100 → 300`) and `2.00` (`300 → 1000`) for every state checked (`m/Λ = 1…1.3`, `J = 0…6`, `k = 2, 3, 4`).
  * `m²ρ` converges, e.g. `−41.2, −42.2, −42.3` (`k = 2`, `J = 0`, `x = 0.51`); the drift `∝ 1/m²` points to a `c₄/m⁴` term.
  * Sizes: `|c₂| ≈ 20–160`, depending on `k`, `J`, `x`; e.g. `k = 2`: `−42` (`J = 0`), `+22` (`J = 2`), `−19` (`J = 4`).
* **What this calibrates.** Any fluctuation expansion of the crystal integral must give **zero at `O(1/m)`** and these `c₂(x, J, k)` at
  `O(1/m²)`. This fixes the target table for the scalar Laplace computation (plan step 1).
* **Interpretation.** The bulk-point corrections are AdS-curvature effects, `O(1/(mR)²)`. This agrees with CLS-D's "`1 + O(J²/m²)`", and with
  the project's `⟨TTTT⟩` Regge-moment fits, where `|c₁| ≤ 0.04` (`SPINNING_TTTT_EXT.md`).

### 19.2 Structure of the graviton correction [D-sketch]

* **Factorised integrand.** In CLS-D's representation, `Ψ⁺[G^s] = ∫ e^{f_{ν,m}} g_scal · S`, where:
  * the fast factor `f_{ν,m}` is spin-independent (CLS-D);
  * `S` is the ratio of slow factors: (the `q²`/`O·q` dual structure) × (celestial `j′` structure) × (heavy `⟨TT𝒪⟩` structures with
    polarisations `x₁₀`, `x₃₀` for 1, 3) × (the polarisation sum for 2, 4), divided by the scalar analogues.
* **Laplace expansion.** It gives, per helicity channel `h`, the stripped density × `S*` × `[1 + (c₂^{scal} + δc₂^{spin}(x, J, h))/m² + O(m⁻⁴)]`.
* **The spin-specific coefficient.** `δc₂^{spin}` is built from:
  * `½A^{ab}∂_a∂_b log S`, `A^{ab}∂_a log g_scal ∂_b log S`, and the cubic/quartic terms of `f` contracted with `∂ log S`, where `A` is the
    inverse Hessian, `O(1/m)` per entry;
  * the `O(1/Δ²)` corrections of the heavy conserved-structure basis (CLS-D §4.1 gives the leading order only);
  * the helicity-dependent normalisation of the shadow coefficient at `O(1/Δ²)`.
* **Odd terms.** These are expected to cancel as for the scalar (curvature expansion). For Regge-moment objects this is supported:
  `|c₁| ≤ 0.04` [NE].
* **Common factors.** A helicity- and `k`-independent part of `δc₂` is a positive overall rescale, like `directgrav`'s `c*`. It cannot change
  a sign. Only the `k`- and helicity-dependent part matters.

### 19.3 H4 targets for route A′ in calibrated units (`m2_tolerance_c2.py`) [NE + M]

The model replaces the graviton density by density × `(1 + c w(x)/m²)` on `x > x₀`, with sign patterns:
* `all`: same `c` for `k = 2, 3` and both helicities;
* `hel`: `c₊₊ = −c₊₋`;
* `k`: `c₂ = −c₃`.

The table gives the interval of `c` keeping `n = 1` (route-A′ functional, f6 model; all states, both skeletons):

| `Δ_g` | `x > 0.3`: `all` / `hel` / `k` | `x > 0.5`: `all` / `k` | all `x`: `k` |
|---|---|---|---|
| 100 | (−1110, +962) / (−1110, +962) / **(−686, +629)** | (−4360, +2820) / (−970, +921) | (−2530, +2410) |
| 300 | (−4980, +7730) / same / (−3350, +5390) | (−15500, +10100) / (−3920, +3380) | (−2.0·10⁴, +1.9·10⁴) |
| 1000 | (−5.5·10⁴, +7.9·10⁴) / same / (−3.5·10⁴, +5.0·10⁴) | — / (−3.2·10⁴, +2.7·10⁴) | (−2.1·10⁵, +2.1·10⁵) |

**H4 for route A′ therefore holds at `Δ_g = 100` if `|δc₂^{spin}(x)| ≲ 630` on `x > 0.3`** (worst pattern, any `x`-dependence of fixed sign);
`≲ 3·10³` at `Δ_g = 300`; `≲ 3·10⁴` at `Δ_g = 1000`.

**Available evidence on `δc₂^{spin}`** (Regge-moment objects, `x ≲ 0.1` only) [NE]:
* `⟨TTTT⟩` total `c₂ ≤ 310` (`k = 2…7`, both helicities, `ν ≤ 3`; `SPINNING_TTTT_EXT.md`), which includes the common factor;
* the non-common shape residual is `≤ 3·10⁻³` at `m ≈ 102`, i.e. `|δc₂| ≲ 30` (`DIRECT_GRAVITON_NEARGAP.md`, `++`);
* scalar `|c₂| ≤ 160` on the window.

**Margins against the target.**
* At `Δ_g = 100`: a factor `≳ 2` against the full `⟨TTTT⟩` `c₂`, and `≳ 20` against the measured spin-specific residual, **if** `δc₂` does not
  grow by more than these factors between `x ≈ 0.1` and `x ≈ 0.7`.
* At `Δ_g ≥ 300`: a factor `≳ 10–100`.
* Growth `∝ (1−x²)^{−n}` (saddle degeneration) gives at most `×4` for `n = 2` at `x = 0.7`. So H4 is **plausible with margin, not shown** [I].

### 19.4 Plan to close M6

| step | content | validation | estimate |
|---|---|---|---|
| 1 | scalar Laplace expansion of the CLS-D crystal integral to `O(1/m²)`: `f_{ν,m}` (explicit in CLS-D) to 4th order, inverse Hessian, scalar slow factors (subtraction factor, dual/celestial structures, shadow coefficient and heavy-block normalisation to `O(1/Δ²)`) | must give `c₁ = 0` and reproduce §19.1's `c₂(x, J, k)` to `≲ 1 %` on `x ∈ [0.3, 0.7]` | 2–3 days |
| 2 | graviton slow factors: `q²`, `O·q` (§17), celestial `j′ = 0, 1`, heavy `⟨TT𝒪⟩` conserved structures to `O(1/Δ²)`, polarisation sums; `δc₂^{spin}(x, J, h)` | `c₁ = 0`; `x → 0` limit vs the Regge-moment `c₂` data (modulo the `Π̂` vs `Ĉ` difference) | 4–6 days |
| 3 | feed `δc₂^{spin}` into the `h2_neargap_exact` pipeline; `n` at `Δ_g = 100/300/1000` | `|δc₂|` vs §19.3 targets | 1 day |

### 19.5 Regraded table (supersedes §18.1 where they differ)

| statement | grade |
|---|---|
| exact scalar remainder on `x ∈ [0.3, 0.7]` is `c₂(x,J,k)/m² + O(m⁻⁴)`, no `1/m`; `|c₂| ≈ 20–160` | [NE] (exact data, three gaps) |
| graviton correction structure `1 + (c₂^{scal} + δc₂^{spin})/m²`, odd terms absent | [D-sketch]; odd-term absence [NE] (Regge-moment `c₁`) + [I] |
| H4 target for route A′: `|δc₂^{spin}| ≲ 630 / 3·10³ / 3·10⁴` on `x > 0.3` (`Δ_g = 100/300/1000`, worst pattern) | [NE + M] |
| H4 (`δc₂^{spin}` within target on `x ∈ [0.3, 0.7]`) | **[I]**, plausible: evidence `|δc₂| ≲ 30` (non-common) and `≤ 310` (total) at `x ≲ 0.1` |


## 20. M6 step 1: Laplace expansion of the CLS-D crystal integral to `O(1/m²)` — scalar, **validated** (`l1_tps.py`, `l1_crystal.py`, `l2_scalar_validate.py`) [DV + NE]

**Result.** The exact CLS-D integral expanded around the scattering crystal to relative `O(1/Δ²)` reproduces the exact scalar skeleton
(e1, Mellin) on `x ∈ [0.3, 0.7]`, `J = 0, 2, 4, 6`, `m/Λ = 1…1.5`, at `Δ_g = 100, 300, 1000`.
* The mismatch of the `1/m²` coefficient is `|Δc₂| ≤ 0.21` (rms `0.05`) in the units of §19.1, where `|c₂|` is up to `200`. The median
  relative mismatch is `0.05 %`.
* No `1/m` mismatch (`c₁`), modulo a `ν`-only factor (see "What is fitted").
* So the step-1 criterion of §19.4 (`c₁ = 0`, `c₂` to `≲ 1 %`) is met, apart from the qualifications below.

### 20.1 What is computed

**The integral.**
* CLS-D `Ψ⁺_{2,ν}[f G^s_{Δ,J}]` with `f = (v′−u′)/(u′v′)`, in `d = 3`, `Δ_φ = 3` (eq. `scalarfunctional_plus_gaugefixed_nozintegral`).
* The delta functions are solved (`t̃₀ = 0`; `|y₀|` from CLS-D's formula, with the Jacobian re-derived:
  `(u_yv_y)^{1/2}(u_yv_y−y⊥²)^{−1/2}(u_yv_y(1−x̃₀²)−y⊥²)^{−1/2}`). In `d = 3`, `S^{d−3} = {n₀ = ±1}`; both give the same value.
* So it is an exact 7-dimensional integral over `(x̃₀, u_x, v_x, u_y, v_y, x⊥, y⊥)`.

**The split.**
* Exponent `λF = (iν/2) log A + (Δ/2) log B`: CLS-D's fast factor with `m → Δ`, exact.
* `G`: the log of all other `z`-dependent factors, exact:
  * the measure, `(−x²)^{−2d}` and the Faddeev–Popov determinant;
  * the remainders of the dual structure (`Δ = 3/2 + iν`, `J = −1`) and of the celestial structure (`δ = 1/2 + iν`, `δ_i = 2`);
  * the subtraction factor;
  * the scalar parts of both 3-pt structures, and their spin-`J` contraction `Q_J = (Z²W²)^{J/2} P_J(Ẑ·Ŵ)`.
* `1/β_{Δ,J}` is exact (CLS-D eq. `shadow_coeff_scalar`; the gauge `x₅ = e` leaves no `x₅` integral).

**Backbone checks against CLS-D [DV, to round-off].**
* `∇F = 0` at the crystal.
* `det(−∂²λF)` equals CLS-D's `Det M_{ν,m,y₀}` exactly, with `m = Δ`.
* `λF` equals CLS-D's closed form `f_{ν,m}` at off-saddle points.
* `f|crystal` equals CLS-D eq. `subtraction_scalar_atsaddle`, and `|y₀| = a`.
* `e^{λF}|crystal` equals CLS-D's `e^{−πν}2^{iν}` up to the `ν`-only branch factor `e^{−πν}`.

**The Laplace expansion (`l1_tps.py`).**
* Exact truncated multivariate power series in 7 variables to total degree 6 (1716 monomials). There are no finite differences.
* Change of variables `z = z* + Lξ` with `Lᵀ(−∂²λF)L = 1` (complex-symmetric Cholesky of `−H`, refined twice).
* The relative correction is `1 + K₁ + K₂`, with
  * `K₁ = E[A₂ + A₁²/2]`,
  * `K₂ = E[A₄ + A₂²/2 + A₁A₃ + A₁²A₂/2 + A₁⁴/24]`,
  * `A₁ = (λF)₃ + G₁`, `A₂ = (λF)₄ + G₂`, `A₃ = (λF)₅ + G₃`, `A₄ = (λF)₆ + G₄` (subscript = homogeneous degree in `ξ`), and Gaussian
    moments taken exactly (up to degree 12).
* `K₁Δ` and `K₂Δ²` are real and `Δ`-independent at fixed `ν/Δ` (checked).
* The neglected term is `K₃ ∝ 1/Δ³`.

### 20.2 The exact comparison object

**From CLS-D to CHMRS.**
* CLS-D `Ψ_{2,ν}` (`f = (v′−u′)/(u′v′)`) is the `J = −1` light-ray/inversion functional `Φ_{−1,ν}[(v′−u′)G/(u′v′)]` (Kologlu et al.; CHMRS
  app. on superconvergence), up to a `ν`-only factor `N(ν)`.
* CHMRS's magic block identity then gives `X(ν) = g(ν) C₂(ν+i) + g(−ν) C₂(ν−i)`, with
  `g(ν) = ((2−d)/2+iν)/(iν) · (γ₇(ν+i)/γ₇(ν))²` and `C₂` = CHMRS `C_{2,ν}`.
* `C₂(ν ± i)` at fixed `m` is `x → x ± i/m`. It is built from the exact e1 value at the node plus the shift terms of the leading
  `E₂(x)P_J(1−2x²)/m⁴`, up to third order; the error is `O(m⁻⁴)`.

**Normalisation.** e1 is the heavy density `C₂[G]/(q_{Δ,J}·2 sin²)` (CHMRS eq. `definitionofheavydensity`). So the comparison is
`Y/(q_{Δ,J} X)`, with `Y = e^{logY₀}(1 + K₁ + K₂)`; this removed a `4^Δ Δ^{−17/2}` mismatch. After that normalisation, `Y₀/q` agrees
with CLS-D's leading heavy action to `O(1/Δ)` [NE].

**What is fitted.** Each gap has about 470 (state, node) points.
* One `N(ν)` per node `ν`.
* One `n(J)`: a `J`-only structure normalisation. The fitted `n(J)` are identical at all three gaps (`0.0937, 0.0171, 0.0035` for
  `J = 2, 4, 6`).
* The phase of `Y/(qX)` is state-independent to `3·10⁻¹¹` rad.
* A `ν`-only factor is a reweighting `ψ(p) → ψ(p)N(Λp)` of the functional, which is irrelevant for H4. The light side is computed
  exactly anyway. What H4 needs, and what is validated, is the dependence on the state `(m, J)` at fixed `ν`.

### 20.3 Numbers (`cache/l2_scalar_validate_*.txt`, `cache/l2_richardson.txt`)

Residuals, weighted fit, in `c₂` units (`Δc₂ = m² ε · ex/E₂`), rms over all points:

| `Δ_g` | `Y₀` only | `Y₀(1+K₁)` | `Y₀(1+K₁+K₂)` |
|---|---|---|---|
| 100 | 296 | 42.4 | 11.6 |
| 300 | 1.05·10³ | 40.8 | 4.2 |
| 1000 | 3.56·10³ | 40.3 | 1.27 |

* **Reading of the table.** Each column scales exactly as its truncation predicts:
  * `Y₀` only: `∝ m`, a missing `1/m` term.
  * `+K₁`: `m`-independent, a missing `1/m²` term of size `≈ 40` (c₂ units).
  * `+K₂`: `∝ 1/m` (`11.6 : 4.2 : 1.27`). This is the `K₃/Δ³` truncation, not a `c₂` error. Fitting `Δc₂ = d₀ + d₁/m + d₂/m²` per matched point over the three gaps gives:
  * **`d₀` (the mismatch of the `1/m²` coefficient itself):** rms `0.054`, max `0.21`, median `0.05 %` of `|c₂|`. This is on 443 points
    with `|P_J(1−2x²)| > 0.05`. The worst ratio is `2.9 %`, at a point with `|c₂| ≈ 5`.
  * **The truncation term** is `|d₁| ~ 10³` (so `K₃/m³`); at `Δ_g = 100` it is about `10/m²`.
* **Without the near-zero cut.** The 23 points closest to zeros of `P_J(1−2x²)` (where `log Q_J` is singular at the saddle, so the
  expansion is non-uniform) reach `|d₀| ≤ 7`.
  * Fix if needed: treat `Q_J` as a polynomial prefactor, not in the exponent.
* **The `1/m` term.** `Y₀` alone misses a `1/m` term (residual `∝ m`). `K₁` supplies it (residual becomes `m`-independent). So the CLS-D crystal reproduces
  `c₁ = 0` of the exact data, modulo `N(ν)`.

### 20.4 Consequences and grades

* **[DV + NE] The `1/m` expansion around the CLS-D crystal is quantitatively correct through `O(1/m²)` for the scalar `k = 2` functional,
  state by state.**
  * The engine (`l1_tps`, `l1_crystal.laplace`) is generic: any spinning slow factor enters only through `G`.
  * So step 2 is a change of `build()`'s slow factors plus the structure normalisations; the expansion machinery is validated.
* **[NE] Usable accuracy for step 3.**
  * At `Δ_g = 100` the `O(1/m²)` truncation is exact to `|Δc₂| ≲ 0.2`; the `O(1/m³)` remainder is `≈ 10` in `c₂` units.
  * Both are far below the §19.3 tolerance `630`.
* **Caveat [I]:** the proportionality `Ψ^{CLS-D}_{2,ν} ∝_ν Φ_{−1,ν}[(v′−u′)·/(u′v′)]` is used, not derived. It is confirmed a posteriori: a
  state-dependent discrepancy would have shown up.

### 20.5 Step 2 status: one ingredient done (§20.7); graviton structures not started

The graviton slow factors still to be built are listed below.
* **On the `(2,4)` side:** the dual structure with the covariant flip-free `q²` (`Ĉ₂`) and `O·q` (`Ĉ₃`) of §17, with `T₂, T₄`
  polarisations, and the polarisation sums.
* **On the `(1,3)` side:** the celestial `j′ = 0, 1` structures with light-ray polarisations `z₁ = x_{01}`, `z₃ = x_{30}`.
* **The heavy `⟨TT𝒪_{Δ,J}⟩` conserved structures:** two parity-even and one parity-odd in `d = 3`, with their exact
  `Δ`-dependent coefficients.
* **The spinning shadow coefficient (matrix).** In the gauge `x₅ = e` it is needed exactly. Alternatively, CLS-D's `x₆` / `r`-tensor
  route, expanded to `O(1/Δ²)` by the same engine (3 variables). **This route is now validated on the scalar `β_{Δ,J}` (§20.7).**
* **The 2×2 (`J` even) heavy OPE-coefficient matrix.** This makes "`δc₂` per helicity" basis-dependent at `O(1/Δ)`. The
  basis-independent statement is matrix positivity, so step 3 must use the 2×2 action.

H4 for route A′ therefore remains **[I]** (§19.5), now with a validated method and an `O(1/m³)` error budget `≈ 10` (c₂ units) at
`Δ_g = 100`.

### 20.6 Regraded table (supersedes §19.5 where they differ; rows for §20.7–§20.8 included; superseded by §21.8 where they differ)

| statement | grade |
|---|---|
| Laplace expansion of the CLS-D crystal integral to `O(1/m²)` reproduces the exact scalar `c₂(x, J)` (`k = 2`) to `|Δc₂| ≤ 0.21` (median `0.05 %`), and `c₁` (mod `ν`-only) | [DV backbone + NE] (three gaps, 443 points) |
| `O(1/m³)` truncation of that expansion | `≈ 10³/m` in `c₂` units [NE] |
| shadow coefficient from the `x₅` saddle (CLS-D route) to `O(1/Δ²)` | [NE] (scalar: `Δ⁻³` residual, `J`-only constant) |
| parity-even conserved `⟨TT𝒪_{Δ,J}⟩` (`d = 3`), exact `Δ`-dependent coefficients at general complex position | [NE] (counts `1, 1, 2, 2` for `J = 0, 2, 4, 6`; matches `SPINNING_TTTT.md`) |
| graviton `δc₂^{spin}(x, J, h)` | not computed; remaining list in §20.8 |
| H4 for route A′ | **[I]** (unchanged; §19.3 targets) |

### 20.7 Step 2, first ingredient: the shadow coefficient from the `x₅` saddle to `O(1/Δ²)` (`l3_shadow_laplace.py`) [NE, validated]

* **Purpose.** For `⟨TTTT⟩` the shadow coefficient is a matrix over `⟨TT𝒪⟩` structures with no closed form.
  * CLS-D get it at leading order from the saddle of the shadow integral at `x₅ = e`, with `x₃, x₄` at the crystal and `x₆ = (−e)^+`.
  * Here that integral is expanded to relative `O(1/Δ²)` with the l1 engine: 3 live variables, plus 4 inert Gaussian dummies.
  * It is tested on the scalar, where `β_{Δ,J}` is known.
* **Result.**
  * `x₅ = e` is a saddle (`|∇| ~ 10⁻¹³`), with `K₁Δ = 13/8` and `K₂Δ² = 2.195 + J(J+1)/2`.
  * `|β_Laplace|/|β_exact| = 2^{1−2J}`, a `J`-only constant from the index-contraction conventions.
  * After `K₂` the residual is `2.8·10⁻⁶ / 1.0·10⁻⁷ / 2.7·10⁻⁹` at `Δ = 100/300/1000` (`J = 0`). It scales as `Δ⁻³`, so it is the
    `K₃` truncation.
  * The result is identical for crystal points `ν/Δ = 0.4` and `0.6`.
  * One error was found and fixed on the way: a missing `x₃₄^{−Δ}` power of `⟨Õφ₃φ₄⟩`.
* **Consequence.** For step 2 the spinning shadow matrix can be obtained to `O(1/Δ²)` by projecting this saddle expansion onto the
  `⟨𝒪†(x₆) T₃T₄⟩` structure basis. No exact weight-shifting computation is needed.

### 20.8 Round 2 of step 2 (unrefereed, not committed): prefactor fix, conserved `⟨TT𝒪⟩` structures

**(a) `P_J`-zero prefactor fix (`l1_tps.expand_pref`, `l1_crystal.laplace_pref`, `l2_scalar_validate.py` with `L2_PREF/L2_USEPREF/L2_LINEAR`) [NE].**
* What changed: `Q_J` is now kept as a polynomial prefactor, with `E[q e^A]` taken to `ε⁴`; the exponent no longer contains `log Q_J`.
* Agreement with the log form: identical to `10⁻¹³` away from zeros. Both are the same `O(1/Δ²)` truncation.
* Effect on the outliers: the Richardson numbers are unchanged, including with sign-carrying linear residuals. So the 23 points near
  zeros of `P_J(1−2x²)` (`|d₀| ≤ 7.8`) were **not** caused by `log Q_J` (my §20.3 diagnosis was wrong).
* Their `Δc₂` is non-monotone across the gaps, e.g. `3.0 / 15.4 / 1.5` at `J = 4`, `x = 0.575`. I suspect, but have not shown, that
  the exact e1 values near the zeros of the leading term are less accurate there [I].
* Size: `≤ 8` in `c₂` units, against the tolerance `630`. Away from the zeros, §20.3 stands.
* The prefactor engine is what step 2 needs: the summed polarisation/structure contractions enter as a prefactor `q`, which can vanish.

**(b) Parity-even conserved `⟨T T 𝒪_{Δ,J}⟩` in `d = 3`, general complex Lorentzian position, TPS-evaluable (`l4_tto_even.py`) [NE].**
* Basis: CPPR monomials `V₁^{a₁}V₂^{a₂}V₅^{a₅}H₁₂^{b₁₂}H₁₅^{b₁₅}H₂₅^{b₂₅}/X…` with the CLS-D `V`, `H` conventions.
* Conditions imposed:
  * conservation of `T₁` and `T₂`, using the `d = 3` Todorov operator `C[f] = (3/2)∂_x·∂_z f − ½(z·∂_x)□_z f`, applied by exact TPS
    differentiation;
  * `1↔2` symmetry;
  * the `d = 3` linear relations are quotiented out numerically.
* Counts (null-space gap `10⁻¹⁵` vs `10⁻¹`, at `Δ = 7.3, 23.7, 101.1`):

| `J` | basis | independent in `d = 3` | conserved and `1↔2` symmetric |
|---|---|---|---|
| 0 | 3 | 3 | 1 |
| 2 | 11 | 10 | 1 |
| 4, 6 | 14 | 13 | 2 |

* The `14 → 13`, 2-conserved count agrees with `SPINNING_TTTT.md` §1 and arXiv:1708.05718.
* Self-test that caught a bug: the operator first failed on `⟨JJ⟩` and `⟨TT⟩` two-point functions. Cause: scalar + TPS was added by
  numpy broadcast. Fixed by converting all inputs to TPS constants.
* `conserved_even(Δ, J)` returns the exact (numerical) `Δ`-dependent coefficient vectors, and `structures()` evaluates them on TPS
  points. They are therefore ready for the crystal integrand.

**Step 2 remaining list (concrete).**
1. **Parity-odd conserved structure.**
   * Either import `fourpt/tttt/odd3.py`'s words or construct `ε`-structures here. In `d = 3`, `ε(Z₁, Z₂, X₁, X₂, X₅)`-type pseudo-scalars
     times even monomials, then the same null-space step.
   * odd3 works in a frame at infinity with sympy, so a covariant re-expression is needed. A same-point comparison with odd3 would be
     the validation.
2. **Light-transformed externals.**
   * After the δ-function rewrite, `T₁, T₃` have polarisation `z₁ ∝ x_{01}`, `z₃ ∝ x_{30}`. Their celestial `j′ = 0, 1` structure is
     contracted with the transverse polarisation; the scalar celestial factor is in `l1_crystal`.
   * `T₂, T₄` are contracted with the `(2,4)` dual structure. That means `q²` (`Ĉ₂`) or `O·q` (`Ĉ₃`) of §17 (`g5_covariant_Q.py`), with
     the polarisation sum over the dual basis (CLS-D §4).
3. **Shadow coefficient matrix `β_{ab}(Δ,J)`.** Use the `l3` route projected onto the 2 (+1) conserved `⟨𝒪†(x₆) T₃T₄⟩` structures. In the
   `x₅ = e` gauge the integrand needs `β⁻¹`.
4. **The integrand.** `build()` gets `G` (scalar parts) plus the prefactor
   `q = Σ_{a,b} λ_a λ_b × [⟨T₁T₂𝒪⟩_a · β⁻¹ · ⟨𝒪̃T₃T₄⟩_b contracted]`. Each `(a,b)` gives one entry of the 2×2 (3×3 with odd) action
   matrix.
5. **Validation at `x → 0`** against the `⟨TTTT⟩` Regge-moment `c₂` (`SPINNING_TTTT_EXT.md`). This is modulo the `Π̂` vs `Ĉ` difference;
   the scalar analogue (`Ĉ₂` vs `Π̂₂`) has to be quantified first with the scalar engine.
6. **Step 3.** Feed the 2×2 matrix (minimum eigenvalue per state) into `h2_neargap_exact`, and get `n` at `Δ_g = 100/300/1000` vs §19.3.

### 20.9 Round 3 of step 2: hand-over recipe for the parity-even 2×2 path (executed in §21)

Priority order (coordinator): even 2×2 action end to end through step 3 first; the odd sector afterwards.
1. **Shadow matrix `β_{ab}` (even).** Extend `l3_shadow_laplace.build`:
   * Replace `(a·W)^J` by `f_b(x₅, z₅ → a; x₃, z₃, x₄, z₄)`, with `a = I(x₅₆)·z₆`. This is `l4.structures` for `⟨𝒪̃(x₅) T₃ T₄⟩` at
     dimension `3 − Δ`, using the conserved combination `conserved_even(3−Δ, J)`.
   * Why the substitution is exact: contracting a traceless spin-`J` tensor with `(a·z₅)^J` for null `a` is evaluation at `z₅ = a`.
   * Expand in `x₅` around `e`, at fixed numerical `(z₃, z₄, z₆)`.
   * Fit `LHS_b(z₃, z₄, z₆) = Σ_a β_{ba} R_a(z₃, z₄, z₆)` over about 10 random polarisation triples. Here `R_a` is `l4.structures` for
     `⟨𝒪†(x₆) T₃T₄⟩` at dimension `Δ`, conserved basis.
   * Validation: residual `∝ Δ⁻³`; `β` independent of the crystal point and of the polarisations.
2. **Externals.**
   * `T₁, T₃`: light-ray polarisations `z₁ = x₀₁`, `z₃ = x₃₀` (null by the δ-functions), times the celestial `j′ = 0` (`Ĉ₂`) or
     `j′ = 1` (`Ĉ₃`) structure. Take its transverse polarisation from CLS-D §4.2.
   * `T₂, T₄`: contracted with the dual structure `q²` (`Ĉ₂`) or `O·q` (`Ĉ₃`), using `g5_covariant_Q.py`'s covariant forms. The
     polarisation sum is the CLS-D §4.3 dual-basis sum: sum over a null basis of `z₂, z₄` with the dual metric.
3. **Assemble `build()`.**
   * `G`: the scalar parts (`l1_crystal`) with the `T` dimensions and spins in the denominators.
   * Prefactor `q_{ab}`: `⟨T₁T₂𝒪⟩_a(z₅ → ·) β⁻¹ ⟨𝒪̃T₃T₄⟩_b`, evaluated on the TPS.
   * Run `expand_pref` per `(a, b)` to get the 2×2 action matrix `M_{ab}(ν; Δ, J)` to `O(1/Δ²)`.
   * Leading-order check: `M` maps onto the stripped flat helicity densities of §9/§17 (channel ratio `(1−x²)^{−4}`).
4. **`x → 0` comparison** with the Regge-moment `c₂` (`SPINNING_TTTT_EXT.md`), after quantifying `Ĉ₂` vs `Π̂₂` on the scalar engine.
5. **Step 3.** Per state, the minimum eigenvalue of the `ψ`-weighted 2×2 action, through `h2_neargap_exact`; `n` at `Δ_g = 100/300/1000`
   vs §19.3.

## 21. M6 step 2, round 3 (§20.9 recipe executed): the graviton `Ĉ₂` crystal action, parity-even `2×2`, to `O(1/m²)` (`l5_*.py`, `l6_spin_action.py`, `l7_*.py`, `g6_odd_structure.py`) — partial; hand-over in §21.9

*Written 18 Nov 2026 by the numerics agent that took over M6; unrefereed, not committed. Done: items 1 and 2–3 (LO). Not done: the `δc₂` extraction (unstable, §21.4), items 4–5, and the odd sector.*

### 21.1 Machinery (`l5_spin.py`; self-tests in `cache/l5_spin_tests.txt`) [DV, round-off]

* **Index contractions in `d = 3`.** Spin-`l` traceless tensors are contracted from their values on the null circle
  `w(φ) = (1, cos φ, sin φ)`, using `⟨T,S⟩ = Σ_{n=−l}^{l} T_n S_{−n}/κ_n`.
  * `T_n` are the Fourier coefficients of `T(w(φ))`, from `2l+1` samples.
  * `κ_n = (−1)^{l+n} 2^{−l} C(2l, l+n)`. Derivation: `T = Σc_k w_k^{⊗l}` and `(w·w′)^l = (cos(φ−φ′) − 1)^l` [D].
  * Checks: `⟨(a·z)^l, (b·z)^l⟩ = (a·b)^l` at random complex null `a, b` to `10⁻¹⁶–10⁻⁷` (`l ≤ 20`), `10⁻⁴` at `l = 28`.
  * On the crystal TPS, the contraction reproduces the scalar `Q_J` of `l1_crystal` (the scalar spin-`J` exchange), coefficient by
    coefficient, up to a `J`-only constant. The error is `7·10⁻¹⁶ / 10⁻¹³ / 2·10⁻¹¹ / 10⁻¹⁰` at `J = 4 / 12 / 20 / 28`, so the crystal
    configuration is benign at high `J`.
* **Building blocks.** `V`, `H`, `X` are batched over TPS and agree with `l4_tto_even` to `10⁻¹⁶`.
  * They are linear in each polarisation, so sampled polarisations cost three basis evaluations.
  * The conserved `⟨TT𝒪⟩` numerators agree with `l4` to `4·10⁻¹⁵`.
* **Conserved coefficients `conserved(Δ, J)`.**
  * `l4.conserved_even` runs on a reduced TPS (6 variables, degree 3, which is all the Todorov operator needs), with point-wise normalised
    denominators so that it works at `|Δ| ~ 10³`.
  * The null-space basis is rotated to a real one. The null space is real, which is checked: rank 2 of `[Re, Im]`.

### 21.2 Item 1: shadow matrix `β_{cb}(Δ, J)`, parity-even (`l5_shadow_matrix.py`, `cache/l5_shadow_matrix.txt`) [NE, validated]

* **Definition.** `∫d³x₅ ⟨𝒪(x₆,z₆)𝒪(x₅,·)⟩·⟨T₃T₄𝒪̃(x₅,·)⟩_c = Σ_b β_{cb}⟨T₃T₄𝒪(x₆,z₆)⟩_b`, with `x₆ = (−e)⁺` and `x₃, x₄` at the
  crystal.
  * The `x₅` integral is Laplace-expanded around `e` to relative `O(1/Δ²)` (the l3 route).
  * The `𝒪`-contraction is evaluation at `z₅ = I(x₅₆)z₆` (exact).
  * `𝒪̃` side: conserved structures at `3 − Δ`. Right side: conserved structures at `Δ`, in the same basis as the functional uses for
    `⟨T₁T₂𝒪⟩`.
  * `β` is fitted over 6 random complex null polarisation triples `(z₃, z₄, z₆)`.
* **Validation (recipe criteria), all three ∝ `Δ⁻³`:**

| `J` | fit residual, order 2, at `Δ = 100 / 300 / 1000` | `|β(ν/Δ=0.4) − β(0.6)|/|β|` | `|β(pol. set 1) − β(set 2)|/|β|` |
|---|---|---|---|
| 0 | `2.9·10⁻⁴ / 1.1·10⁻⁵ / 2.9·10⁻⁷` | `≤ 2·10⁻¹³` | `4.3·10⁻⁵ / 1.3·10⁻⁶ / 3.2·10⁻⁸` |
| 2 | `3.3·10⁻⁴ / 1.3·10⁻⁵ / 3.5·10⁻⁷` | `9.8·10⁻⁵ / 4.0·10⁻⁶ / 1.1·10⁻⁷` | `3.9·10⁻⁴ / 1.4·10⁻⁵ / 3.8·10⁻⁷` |
| 4 | `6.1·10⁻⁵ / 2.4·10⁻⁶ / 6.7·10⁻⁸` | `4.1·10⁻⁴ / 1.7·10⁻⁵ / 5.1·10⁻⁷` | `9.4·10⁻⁴ / 2.3·10⁻⁵ / 7.1·10⁻⁷` |
| 6 | `4.5·10⁻⁶ / 1.7·10⁻⁷ / 4.6·10⁻⁹` | `5.1·10⁻⁴ / 2.0·10⁻⁵ / 7.6·10⁻⁷` | `8.8·10⁻⁴ / 3.0·10⁻⁵ / 9.2·10⁻⁷` |

  At order 0 / 1 the residuals scale as `Δ⁻¹ / Δ⁻²` (e.g. `J = 4`: `1.6·10⁻² / 8.1·10⁻⁴` at `Δ = 100`). So the polarisation dependence of
  the shadow integral is exactly that of a combination of the conserved structures, up to the expected truncation.
* **Use.** A `β` error of `4·10⁻⁴` at `Δ = 100` would be `≈ 4` in `c₂` units. In §21.4 it is removed with the other `1/m³` terms by
  the multi-`Δ` fit.

### 21.3 Items 2–3: the integrand, and the leading-order check (`l6_spin_action.py`) [NE, validated at LO]

**Integrand** (`Ψ⁺` of CLS-D eq. `spinningfunctional_definition`, same gauge and δ-function solution as `l1_crystal`):
* **Unchanged from the scalar** (validated, §20): the measure, `(−x²)^{−2d}`, the FP determinant, the scalar part of the dual
  (`Δ̃_{2,4} = 0`, `J = −1`; `Δ_T = Δ_φ = 3`), and the celestial `j′ = 0` structure with `δ_i = 2`.
  * The light-transform weight `(−α)^{−Δ−J}` and the polarisation scaling `z₁ = −αx₁₀` combine to `(−α)^{−Δ}`. So the celestial factor
    is the scalar one [D].
* **Subtraction factor:** `f⁺_{2;2} = v′ − u′` (instead of the scalar `(v′−u′)/(u′v′)`).
* **`(2,4)` structure:** `q(z₂,z₄)²/V_{0,42}⁴` times the scalar dual.
  * `q = (2H₂₀H₄₀ + V_{0,42}²H₂₄)/(X₂₀X₄₀)` (§17).
  * `V⁻⁴` turns the `J = −1` homogeneity into that of the `J = 3` light-ray operator: it changes the `X₀` degree by `−4`, the `Z₀`
    degree by 0 and the `X₂`, `X₄` degrees by 0. This equals CLS-D's `S_n = (−2V)^{J_F−4+2n}…` with `q² = 4S₀ + S₁ + S₂/16` [D].
* **Polarisations:** `T₁, T₃` take `z₁ = x₁₀`, `z₃ = x₃₀`; `T₂, T₄` and `𝒪` are contracted via §21.1.
* **Block `G_{(a′b′)}`:** `⟨T₁T₂𝒪(e)⟩_{a′} · β⁻¹ · ⟨𝒪̃(e)T₃T₄⟩_c`.
  * The Δ-parts of the 3-pt denominators equal CLS-D's fast exponent `(Δ/2)log B` exactly. The remainders are slow.
  * All polarisation dependence is one prefactor `P_{a′c}`, expanded with `expand_pref`.
  * `M_{a′b′} = Σ_c P_{a′c}(β⁻¹)_{b′c}`, a `2×2` matrix for even `J ≥ 4` and `1×1` for `J = 0, 2`.

**Leading-order check** (`cache/t_lo.pkl`; `J = 4, 6`; `Δ = 1000, 4000`; `x = 0.2 … 0.8`):
1. The order-0 `M(x)` lie in a 2D span: the third singular value of the `(x × 3)` stack is `2.2·10⁻⁷ → 1.5·10⁻⁸` (`J = 4`) and
   `1.2·10⁻⁵ → 7·10⁻⁷` (`J = 6`), i.e. `∝ Δ⁻²`. The span has two rank-1 directions `aaᵀ`, `bbᵀ`, which are x-independent.
2. In that basis:
   * one density ÷ the scalar `k = 2` action is `(1−x²)²` exactly (`J`-independent; e.g. `1.6384, 1.3688, 1, 0.5929, 0.2304` at
     `x = 0.2…0.8`, vs `(1−x²)²` normalised: `1.6384, 1.3689, 1, 0.5929, 0.2304`);
   * the other ÷ the first is `d^J₄₄(θ)/((1−x²)⁴P_J(cos θ))` times an x-independent constant (`J = 4`: `1.1110 ± 0.0005`, all `x`).
3. So **`M` at leading order is `ρ₊₊ aaᵀ + ρ₊₋ bbᵀ`, with the stripped flat helicity densities of §9/§17 (channel ratio `(1−x²)^{−4}`)**.
   This is the recipe's item-3 criterion.
   * It also reproduces, as an output, the graviton/scalar envelope ratios `(1−y)^{±2}` that `h2_neargap_exact`/`f6` assume.
   * Helicity labels are assigned by these shapes, with no external input.
   * `M` is symmetric to `10⁻⁵`.

### 21.4 Item 5 input: the `O(1/m²)` graviton-specific correction — NOT extracted reliably (`l7_*.py`, `cache/l7_*`) [NE: negative / diagnostic]

**Data.** `cache/l7_J{4,6}.pkl` (`l7_dc2_scan.py`; `J = 8, 10, 12, 0, 2` were still running at hand-over, resumable).
* Grid: physical `x = ν/m` at `0.05, 0.1, 0.2, …, 0.8`; `Δ = 100, 200, 400, 800, 1600, 3200`.
* Per point: the `2×2` graviton `Ĉ₂` action at orders 0/1/2, and the scalar `k = 2` action (`l1_crystal.laplace_pref`) at the same `(ν, Δ, J)`.

**What is robust [NE]:**
1. The leading-order structure of §21.3 holds on the whole grid.
   * The order-0 third singular value is `∝ Δ⁻²`, e.g. `J = 4`: `3.3·10⁻⁷ … 3.1·10⁻¹⁰` for `Δ = 100 … 3200` (unnormalised rows).
   * The helicity shapes approach the stripped densities as `1/Δ`.
2. **A ν-only factor.** The graviton/scalar ratio contains `N(ν) ≈ exp[(−0.6 + 21.5i)/ν + (115…123 + 32…46i)/ν² + …]`.
   * It is the same for `J = 4, 6` and both helicities, to the fit accuracy. This is what a ν-only normalisation difference does: the
     dual/celestial Γ-factors of the `J = 3` versus `J = −1` duals were dropped.
   * It is the analogue of the scalar `N(ν)` that §20.2 fitted away. It is irrelevant for H4 only after it is absorbed into `N̂₂(ν)`.
   * It is non-perturbative at small `ν` (`|21.5/ν| ≈ 4` at `x = 0.05`, `m = 100`). So a normalisation at small `x_ref` contaminates every `x`.
3. **Helicity mixing at `O(1/m)` (new).**
   * In the `x`-independent helicity basis (order-0 rank-1 vectors, at each `Δ`), the order-2 matrices have an off-diagonal element
     `M̃₊,₋ = μ(x, J)/m × O(diag)`. The order-1/2 stacks leave the 2D span at `O(1/m)`.
   * The profile `m·|M̃₊,₋|/max|diag|` is identical at `Δ = 800` and `3200`, so this is a clean `1/m` law: `J = 4`: `0.2 … 12.8` on
     `x = 0.2 … 0.8`; `J = 6`: up to `24` at `x = 0.4`.
   * Mechanism (my reading [I]): `q` is flip-free only at the crystal (§17). Its gradient times the cubic fluctuation term gives a helicity
     flip at `O(1/m)`.
   * It enters positivity at second order. Its `c₂`-equivalent `m²|M̃₊,₋|²/|M̃₊₊M̃₊₋|` (relative to the smaller density) is `5 … 484`
     (`J = 4`, `x = 0.3 … 0.7`, the maximum at `x = 0.6`) and `≤ 691` (`J = 6`, at `x = 0.4`). It is up to `2.7·10³` at `x = 0.8`.
     These numbers are **comparable to the §19.3 tolerance `630` (`Δ_g = 100`)**. They are pointwise in `x`. The `ψ`-integrated `2×2`
     `F` per state has not been formed.

**What is not robust [negative].** The per-helicity `δc₂^{spin}(x)` and the helicity-relative `δc₂^{+−} − δc₂^{++}` are unstable.
* The fit `log R = Σ_{k≤K} e_k(x)/m^k` over `m ≈ 97 … 3200` moves `δc₂` by `1–2·10³` between `K = 3, 4, 5` and between `Δ`-subsets
  (`cache/l7_hel_relative.txt`).
* The helicity ratio, which is free of ν-only factors and of the scalar, still has an imaginary, `J`-independent, `x`-growing `e₁(x)`
  (`≈ +0.9i … +5.5i` on `x = 0.1 … 0.8` at `K = 4`).
* So the data are not yet described by a short `1/m` series with real coefficients. **No `δc₂^{spin}` number is quoted.**
* Hypothesis tested and refuted: summing the two `S⁰` points. The `n₀ = −1` integrand equals the `n₀ = +1` one to round-off (checked at
  `J = 4`, `Δ = 400, 1600`, `x = 0.3, 0.6`, all orders).
* Untested candidates:
  * a relative phase convention between the `T₁T₂` and `T₃T₄` sides (`β` fitted with complex polarisations);
  * the `Δ`-dependence of the real null-space basis rotation, which is arbitrary per `Δ` but only enters through `β`;
  * a non-analytic ν-only piece in the per-helicity ratio not cancelled by the normalisation at `x_ref`;
  * genuinely large `1/m³–1/m⁴` coefficients. `m²(log order 2 − log order 1)` is `≈ −7·10³` at `Δ = 100`, so the `K₂` term alone is
    `O(70 %)` at `m = 100`. The Laplace series of the spinning integrand converges much more slowly than the scalar one. That alone
    could explain the fit instability.

### 21.5 Item 4 (`x → 0` vs the `⟨TTTT⟩` Regge-moment `c₂`): not done [ND]

* The Regge-moment `c₂` (`SPINNING_TTTT_EXT.md`, `|c₂| ≤ 310` at `ν ≤ 3`) is an absolute correction to the leading heavy formula for `Π̂_k`.
* Here only ratios are computed: graviton/scalar at fixed `(ν, Δ, J)`, normalised per state. The small-`x` region is dominated by the
  non-perturbative ν-only factor of §21.4.
* The prerequisite `Ĉ₂` vs `Π̂₂` comparison on the scalar engine would need a Regge-moment integrand, which the crystal engine does not
  provide. Not attempted.

### 21.6 `Ĉ₃` ingredient: the covariant parity-odd `(2,4)` structure (`g6_odd_structure.py`, `cache/g6_odd_structure.txt`) [DV at the crystal image]

* **The `d = 3` celestial `j′ = 1` structure.** It is the `j′ = 0` one times the 1D sign function `ε(z₁,z₃,z₀)/|…|^{1/2}`, which is
  locally constant (`±1`). So it carries no fluctuation corrections; at fixed `n₀` it is a constant [D].
* **The covariant odd photon structure.**
  * The candidates are the ε-dualised `j′ = 1` monomials of CLS-D eq. `3pt_convention_spinning` (`U_{0,24}`-type).
  * **`o := ε(Z₂,X₂,X₄,X₀,Z₀)H₄₀/√(X₂₄X₂₀X₄₀³) + ε(Z₄,X₄,X₂,X₀,Z₀)H₂₀/√(X₂₄X₄₀X₂₀³)`** has crystal image exactly
    `2i·O/(4ν²(m²−ν²)²)`, with `O = 𝒲₂𝒱₄ − 𝒱₂𝒲₄`, on all four entries (spread `≤ 5·10⁻¹⁵`) at `x = 0.1 … 0.9`.
  * The `x`-dependence is the same as that of `q ↦ −Q/(4ν²(m²−ν²)²)`. So `o·q` has the B₃ flat limit with no extra `x`-factor.
  * It is odd under the `S⁰` reflection, as required for the `j′ = 1` sum.
  * The `V`-type completion `V₀₄₂[ε(Z₂,…)V₄ − V₂ε(Z₄,…)]` has the same image. The two differ beyond the crystal, so `Ĉ₃`'s `O(1/m²)`
    correction depends on this choice (a functional-design freedom).
  * Note: the naive integer-power normalisation `/(X₂₄X₄₀)` gives a wrong `x`-dependence. The `√` powers are required by the
    homogeneity.
* **Not done:** the `Ĉ₃` integrand. It is `l6` with `q² → o·q`, `f = −1`, and the same celestial. Also not done: its scan, and its
  leading-order check (same helicity vectors as `Ĉ₂`; ratio to `Ĉ₂` `= σ_h · y/(2−y)` × const).

### 21.7 Item 5 (step 3, `n` at `Δ_g = 100/300/1000`): not done

There is no reliable `δc₂^{spin}` (§21.4), so nothing was fed into `h2_neargap_exact`. **No `n` and no `δc₂`-vs-tolerance statement is
made.** The mixing of §21.4 (item 3 there) is a new, positivity-relevant `O(1/m²)` effect of tolerance size. It has to be included as the
minimum eigenvalue of the `ψ`-integrated `2×2` `F`, not per helicity.

### 21.8 Regraded table (supersedes §20.6 where they differ)

| statement | grade |
|---|---|
| `d = 3` spin-`l` contraction by null-circle Fourier coefficients; reproduces the scalar `Q_J` on the crystal TPS (`J ≤ 28`) | [D] + [DV] |
| conserved parity-even `⟨TT𝒪⟩` coefficients at `|Δ| ≤ 3200`, real basis | [NE] |
| spinning shadow matrix `β_{cb}(Δ,J)` (even) to `O(1/Δ²)`: residual, crystal-point and polarisation dependence all `∝ Δ⁻³` (`J = 0…6`) | [NE], validated |
| covariant `Ĉ₂` integrand (`q²/V⁴` dual, `f⁺_{2;2}`, conserved `⟨TT𝒪⟩·β⁻¹·⟨𝒪̃TT⟩`) | [D] (homogeneity, structure) |
| LO: `Ĉ₂` on even blocks `= ρ₊₊aaᵀ + ρ₊₋bbᵀ`, `ρ₊₋/ρ₊₊ ∝ d^J₄₄/((1−x²)⁴P_J)`, `ρ₊₊/scalar = (1−x²)²` (`J = 4, 6`; `Δ ≤ 3200`) | [NE], validated (recipe item 3) |
| ν-only factor `N(ν)` in graviton/scalar, `J`- and helicity-independent | [NE] |
| `O(1/m)` helicity mixing of the `Ĉ₂` heavy action; second-order `c₂`-equivalent up to `~5·10²` on `x ∈ [0.3,0.7]` (`J = 4, 6`) | [NE] (clean `1/m` scaling); mechanism [I] |
| graviton-specific `δc₂^{spin}(x, J, h)` | **not extracted** (fits unstable at `±2·10³`; §21.4) |
| `Ĉ₃` odd `(2,4)` structure `o` (image `∝ O`, B₃-compatible normalisation) | [DV] at the crystal image; integrand not built |
| `x → 0` comparison with the Regge-moment `c₂` | [ND] |
| parity-odd `⟨TT𝒪⟩` sector (incl. the binding state `+−`, `J = 5`) | [ND] |
| H4 for route A′ | **[I]** (unchanged). New risk: the `O(1/m)` mixing is of tolerance size pointwise. |

### 21.9 Hand-over (agent context exhausted; next steps in order)

1. **Stabilise the `1/m` extraction.**
   * Fit on `Δ ≥ 400` only, with `Δ` up to `10⁴`. Cost: `conserved()` and the shadow fit are `≈ 1 min` per `(Δ, J)`; one point is
     `5–15 s`.
   * First check whether the helicity ratio's imaginary `e₁(x)` survives with real (not complex) random polarisations in the `β` fit, and
     with the `S0/S2/S4` orders separately. `ord2 − ord1` is `O(70 %)` at `Δ = 100`, so the spinning Laplace series may need `K₃`.
     `l1_tps.DEG = 6` supports only `K₂`. `K₃` needs `DEG = 8` (memory: `N = 3432`, `_I ≈ 1.2·10⁵`).
2. **Mixing.** Form `F = ∫ψ₂M(Λp)` per state, as a `2×2`, on the h2 nodes: interpolate `μ(x, J)/m` in `x`, using `cache/l7_J*.pkl`.
   Then report `λ_min(F̃)/N_s` against the diagonal. This is decisive for H4 independently of item 1.
3. **`Ĉ₃`.** Copy `l6.build_spin` with `D = q² → o·q`. `o` comes from `g6` (ε by cofactors with TPS; `x₀`-dependent rows), with
   `f = −1`. LO check as in §21.6. Then the `k`-pattern difference `δc₂^{(3)} − δc₂^{(2)}`, where the ν-only factor does **not** cancel
   automatically (dual normalisations differ).
4. Odd `J` / parity-odd `⟨TT𝒪⟩` (binding state). An ε-structure null space in `l4`, compared with `fourpt/tttt/odd3.py` at matched points.

## 22. M6 round 4: the `Ĉ₂` normalisation and phase fixed, the `O(1/m)` mixing explained, and the `2×2` near-gap test (`l8_*.py`) — `Δ_g ≥ 300` done for `Ĉ₂`; `Δ_g = 100` open

*Written 18 Nov 2026 by the numerics agent continuing M6 after §21.9. Unrefereed, not committed. Updated incrementally; §22.6 is the hand-over.*

### 22.1 Diagnosis: why the §21.4 fits were unstable [NE, `l8_analysis.py`, `l8_points.py`]

* **The ν-only factor sits inside the Laplace terms, not in the prefactor.** In the real saddle helicity basis (below), the
  graviton diagonal elements have `m·K₁ ≈ 20 + 22i/x` for both helicities and every `Δ` (`J = 4`: `m·Im K₁ = 217i` at `x = 0.1` for
  `Δ = 100…3200`), and `m²K₂ ≈ −(22/x)²/2`, i.e. the Laplace series is expanding an exponential `exp(22i/ν + 82/ν² + …)`.
  * `l7` summed `S0 + S2 + S4` and then took the log. At `ν = 30` the dropped `(22/ν)³/6` is `≈ 6 %`, i.e. `≈ 600` in `c₂` units at
    `m = 100` and different at every `(x, m)`. This, and the normalisation at `x_ref = 0.05` (`ν = 5`), produced the `10³`-level fit
    instability and the "`O(70 %)`" `K₂`.
  * **Fix 1 (cumulant form):** `log E[p e^A] = log p₀ + (p₂+p₄)/p₀ − ½(p₂/p₀)² + O(ε⁶)`, applied to each helicity-projected element.
    `l8_points` stores the orders separately.
* **Fix 2: exact removal of every ν-only factor.** `D_h(ν; Δ) := L_h(ν; Δ) − L_h(ν; Δ_ref)`, with `L_h = log(M̃_hh/(Y_scal W_h))`
  (graviton helicity element ÷ validated scalar crystal action ÷ the §21.3 LO shape) and a heavy reference state `Δ_ref = 2·10⁴` at
  the **same `ν`**.
  * This is literally §R.3's definition of `N̂₂(ν)` (the `m → ∞` limit at fixed `ν`), including the Laplace-truncation of the ν-only
    factor, which cancels. No analytic dual normalisation is needed: any `J = 3` vs `J = −1` dual/celestial Γ-factor is ν-only and
    drops out identically.
  * For the record: the scalar-type ratio `[Γ((Δ+J)/2)²/Γ(Δ+J−2)]_{J=3}/[…]_{J=−1}` at `Δ = 3/2+iν` gives `exp(−6i/ν + …)`, not the
    fitted `exp(22i/ν + 82/ν²)`. So the spin-3 dual normalisation is not the scalar formula at `J = 3`; CLS-D give the spinning dual
    coefficients `α^{(a)}_{(b)}(ν)` only at large `ν` (their app. on dual structures). Not needed here.
  * Checks: (i) reference independence: `Δ_ref = 2·10⁴` vs `4·10⁴` differ by a per-state constant to `≤ 3·10⁻⁶` (`J = 4`, all `x`);
    (ii) the ν-only factor is **helicity-independent**: `L₊₋ − L₊₊` at the reference is constant in `ν = 20…2200` to `3·10⁻⁸`, with
    relative phase `0` to `2·10⁻⁶π`. So `C₊₊C₊₋ > 0`: no relative sign or phase between the helicity channels.
* **Reality of the ingredients [NE].** The shadow matrix `β` is real (`1.5·10⁻⁶ / 2.9·10⁻⁸` at `Δ = 400 / 1600`). The saddle helicity
  vectors (order-0 rank-1 decomposition, §21.3) are real up to `7·10⁻⁶ / 1.3·10⁻⁸` at `Δ = 400 / 3200` (`∝ Δ⁻²`). So the basis
  `T₀(Δ)` used below is real, and it is common to `Ĉ₂` and `Ĉ₃` (same saddle, same `⟨TT𝒪⟩` structures; only the dual and `f` differ).

### 22.2 The imaginary `1/m` term and the `O(1/m)` mixing: a genuine imaginary part of `Ψ_ν[q²]`, removed by `Re` [D for the prescription; NE for the numbers; I for the mechanism]

* **What survives Fix 1–2.** After the exact ν-only removal, the action is real up to a clean `i e₁ⁱ(x)/m`:
  * diagonal: `m·Im D₊₊ ∈ [−2.3, 0]`, `m·Im D₊₋ ∈ [0, +2.3]` on `x = 0.05…0.8` (`J = 4, 6`), identical for `Δ = 400…6400`;
  * **off-diagonal: the whole `O(1/m)` helicity mixing of §21.4 is imaginary** relative to the (real) diagonal. In pipeline units
    (`(2−y)·offp/m`, see §22.3), `m·Im` is `−0.3, −0.9, −2.1, −4.2, −7.2, −11.5, −17.1` at `x = 0.2…0.8` (`J = 4`), identical for
    `Δ = 400…6400`; the real part falls like `1/m²` (`m·Re = −0.7, −0.3, −0.2, −0.1` at `x = 0.8`, `Δ = 400…3200`; but `+4.7` at `Δ = 100`). Same for
    `J = 6…10`; the real part at `Δ = 100` is up to `m·|Re| ≈ 50` (`J = 6`), see §22.3. The §21.4 "`c₂`-equivalent up to 691" is `m²|κ|²` of this imaginary element (e.g. `J = 4`, `x = 0.7`: `18.8² = 353`, as there).
* **Prescription [D].** OPE coefficients are real in the real structure basis, so `Σλ_aλ_b M_ab(ν) + (light) = 0` holds for the real and
  imaginary parts separately: both are valid sum rules. Define `Ĉ₂ := Re[Ĉ^{CLS}_2/N̂₂(ν)]` (amendment to §R.3). Unchanged by this:
  the crystal image and flat limit (the leading action is real), the light side (`−Ĉ₂|_E = 16πG/(ν²+9/4)` is real), the double zeros
  on `[TT]`, convergence. Under `Re` the imaginary diagonal term enters as `−(e₁ⁱ)²/2 ≥ −3` (`c₂` units) and the imaginary mixing drops
  out entirely.
* **Mechanism [I].** `q²` is a fixed (ν-independent) combination of continuous-spin structures. Complex conjugation maps `Ψ_ν` to the
  functional at `−ν`, and the shadow relation brings it back only if the structure is a shadow eigenvector. For a single scalar dual
  (`J = −1`) it is; for `q²` it need not be, and the mismatch is an `O(1/ν)` admixture of another (also valid) functional. This
  explains an `i·O(1/m)` term that vanishes at the crystal. Not derived.
* **Consequence for §21.4/§21.8.** The "`O(1/m)` helicity mixing of tolerance size" is withdrawn as a positivity risk: it is the
  imaginary part. The real off-diagonal is `O(1/m²)`, its second-order effect `O(1/m⁴)`.

### 22.3 `δc₂^{spin}(x, J, h)` extracted (`l8_refscan.py`, `l8_dc2.py`; `cache/l8_dc2_J*.txt`, `cache/l8_summary.txt`) [NE]

**Grid.** `J = 0, 2, …, 12`; `Δ = 100, 200, 400, 800, 1600, 3200, 6400`; `x = 0.05, 0.1, 0.2, …, 0.8`; reference `Δ_ref = 2·10⁴` at
every `ν = x·m(Δ)`. About 10 s per point, one core.

**Definitions** (all relative to the Regge point `x → 0` of the same state; that per-state constant is a convention for `Ĉ₂` alone, and
fixed only jointly with `Ĉ₃`, i.e. it is the "`k`" pattern of §19.3):
* `++`: `δc₂^{spin,++}(x) = m²[Re D₊₊(x) − Re D₊₊(0)]` (graviton ÷ scalar, so the exact scalar remainder that the pipeline already
  contains is divided out).
* `+−`: the pipeline models the `+−` remainder by transplanting the absolute scalar remainder `ρ₂` onto `d^J₄₄` (§19, `h2`: `fac`).
  The graviton's own `+−` relative correction is `E₋(x) = m²[Re D₊₋(x) − Re D₊₋(0)] + s(x)`, with `s` the scalar crystal's relative
  correction (`= m²ρ₂/P_J`, poles at zeros of `P_J`, which cancel in `E₋`). So `δc₂^{spin,+−} = E₋ − m²[ρ₂/d^J₄₄ − ρ₂(0)]`.
* `x → 0` constants: from the absolute form `w·v = c·w + b x² + d x⁴` on `x = 0.05, 0.1, 0.2` (`w` = LO shape), which is regular
  through zeros of `P_J`, `d^J₄₄`.

**Stability (the §21.9 item-1 criterion).** At fixed `x`, `m²·Re D` is now `Δ`-independent up to a `1/m` drift on `Δ = 400…6400`,
e.g. `J = 4`, `+−`, `x = 0.3`: `1946, 1937, 1932, 1930` at `Δ = 400, 800, 1600, 3200` (§21.4 had `±2·10³` jumps). `Δ = 6400` is noisier
(`≈ 10⁻⁵` relative in the basis/shadow fit, amplified by `m²`) and is not used in fits. The scalar crystal `s(x)` converges the same
way (`J = 4`, `x = 0.5`: `−213, −204, −198, −195, −194, −194`).

**Numbers: asymptotic absolute corrections** (`a₀ + a₁/m` fitted on `Δ = 400…3200`; `A₊ = P_J δc₂^{++}`, `A₋ = d^J₄₄ E₋`, i.e. `c₂` units
times the LO shape normalised to 1 at `x = 0`; max over `x ∈ [0.3, 0.7]`). The absolute forms are the exact (not linearised) differences
`A₊ = P_J(m²+s)(e^{δc₂^{++}/m²}−1)` and `A₋ = d^J₄₄[(m²+s)e^{m²ΔD₊₋/m²}−m²]`, which stay regular through zeros. Grid points with
`|P_J| < 0.05` (or `|d^J₄₄| < 0.02` for `+−`) are excluded, because the log form loses the sign of the scalar remainder exactly at a zero
(`J = 10`, `x = 0.4`: `P₁₀ = 0.002`):

| `J` | `A₊` asympt. | `A₊` at `m ≈ 100` (direct) | `A₋` asympt. | `A₋` at `m ≈ 100` | real mixing `m·|Re|` (pipeline units): `Δ = 100 / 400 / 3200` | imaginary mixing `m·|Im|` |
|---|---|---|---|---|---|---|
| 0 | 275 | 256 | — | — | — | — |
| 2 | 131 | 184 | — | — | — | — |
| 4 | 86 | 321 | 3 | 70 | 4.7 / 0.7 / 0.1 | 17.1 |
| 6 | 85 | 783 | 375 | 191 | 49.7 / 4.8 / 0.6 | 7.7 |
| 8 | 112 | 1752 | 467 | 2.2·10⁴ | 7.8 / 2.8 / 0.4 | 6.5 |
| 10 | 47 | 5.0·10³ | 660 | 2.6·10³ | 15.5 / 4.2 / 0.9 | 9.4 |
| 12 | unusable | | | | | |

* Relative values `δc₂ = A/shape` are larger near zeros of `P_J` and `d^J₄₄` (e.g. `J = 4`, `++`: `+777` at `x = 0.6`, `−3130` at `x = 0.8`);
  there the absolute correction is what enters `F`.
* The `J = 10` asymptotic fit (`Δ = 400…3200`) has residuals `≈ 1.2·10³`: `Δ = 400` is not yet asymptotic at `J = 10`.
* **`m ≈ 100` is outside the controlled regime for `J ≥ 6`.** The direct crystal values at `Δ = 100` differ from the asymptotic form by
  up to `1.5·10³`, and their `x → 0` constant uses `ν ≤ 20`. `J = 12` is not usable below `Δ ≈ 400` (`J²/m ≳ 1`, the per-state
  constant jumps by `10⁴–10⁵`). This is a limit of the Laplace expansion at `J²/m = O(1)`, not a normalisation problem.
* For `+−` the graviton's own correction `E₋` is small for `J = 4` (`≤ 20` asymptotically): the transplanted scalar remainder of the
  pipeline is the larger effect, and `δc₂^{spin,+−}` mostly removes it.
* `J = 0, 2` (`1×1`): `δc₂^{++}` grows towards large `x` (`J = 0`: `6, 16, 33, 64, 125, 274, 792` at `x = 0.2…0.8`), with the same sign
  for `Ĉ₂` at all `x`.

### 22.4 Step 3: the `ψ`-integrated `2×2` per state and `n` (`l8_neargap.py`; `cache/l8_neargap_{100,300,1000}*.txt`) [NE + M]

**What is computed.** Per state (`m/Λ`, `J`, `Δ_f`) of the exact-skeleton set (`load_states`):
* `F = diag(base₊₊, base₊₋) + ΔF`, where `base` is the `m2`/`h2` per-helicity action without inversion (the inversion shifts it by
  `≤ 3·10⁻³ N_s`, §13).
* `ΔF_hh = Σ_nodes ψ₂w (2−y)·[A₊ or A₋/(1−y)⁴ − (transplanted model remainder)]/m²`, and
  `ΔF₊₋ = Σ ψ₂w (2−y)·Re(off)/m`.
* `λ_min` of `N^{−1/2}FN^{−1/2}`, `N = diag(N_s)`. `n` as in `h2`.
* `t*` is the scale factor on the whole computed correction `ΔF` at which `n = 1` first fails, for each sign.
* Two models for the `m`-dependence: `asym` (`a₀ + a₁/m`, fitted on `Δ = 400…3200`) and `direct` (the crystal values interpolated in
  `1/m` through `Δ = 100…3200`). `absmix` is `asym` with `|Re + i Im|` of the off-diagonal, i.e. what §21.4 implied.
* Coverage: even `J ≤ 10` (both helicities; `1×1` for `J = 0, 2`). Odd `J` (the parity-odd sector, incl. the binding `+−`, `J = 5`) and
  `J > 10` are left unchanged. `Ĉ₃` is not corrected. `x > 0.8` carries `2.4·10⁻³` of the `ψ₂` weight and is held at `x = 0.8`.

| `Δ_g` | coverage | `n` (`asym` / `direct`) | min `λ/N_s` (`asym` / `direct`); uncorrected | `t*` (`asym`) | `t*` (`direct`) | `absmix`: `n`, `t*` |
|---|---|---|---|---|---|---|
| 1000 | `J ≤ 10` even | **1 / 1** | 0.2021 / 0.2021 (binding: `+−`, `J = 5`, uncovered); 0.2021 | (−14.6, +335) | (−13.7, +290) | 1, (−1.3, +1.5) |
| 300 | `J ≤ 10` even | **1 / 1** | 0.2054 / 0.2054 (`J = 5`); 0.2054 | (−2.7, +4.6) | (−4.4, +6.7) | 1.27 (fails) |
| 100 | `J ≤ 4` | **1 / 1** | 0.2256 / 0.2200 (`J = 4`, `m/Λ = 1.04`); 0.2256 | (−12.7, +16.7) | (−12.6, +8.1) | 1, (−5.6, +6.2) |
| 100 | `J ≤ 8` | 1 / — | 0.103 (`J = 8`, `m/Λ = 1.08`) / overflow | (−0.5, +1.2) | — | 1.21 (fails) |
| 100 | `J ≤ 10` | **1.5** / — | −2.70 (`J = 10`, `m/Λ = 1.1`) / overflow | (−0.5, +0.2) | — | 1.6 (fails) |

**Reading.**
* **`Δ_g = 1000, 300`: `n = 1` survives the computed `Ĉ₂` correction (even `J ≤ 10`), with margin factors `t* ≳ 14` and `≳ 3`.**
  The binding state stays `+−`, `J = 5`. Both `m`-models agree.
  * The margin at 300 is set by `J = 10`, whose fits are the least asymptotic (§22.3). Treat `t* ≈ 3` as a lower bound on the
    robustness, not a precise number.
* **`Δ_g = 100`: not decided.**
  * For `J ≤ 4`, where the expansion is controlled at `m ≈ 100`, `n = 1` with `t* ≥ 8`.
  * For `J = 6…10` at `m/Λ ≈ 1.05–1.1` (`m ≈ 105`, `J²/m ≈ 0.3–0.95`), the crystal expansion is not controlled. The `direct` values
    overflow (`(m²+s)e^{δ/m²}` with `δ ~ 10⁶` at `J = 6`). The `asym` extrapolation from `Δ ≥ 400` breaks `n = 1` at `J = 10`
    (`λ_min/N_s = −2.7`) and leaves no margin at `J = 8`.
  * This is a real risk signal, but it is not a computed failure: the `asym` form is used far outside its fit range. A controlled
    statement at `m ≈ 100`, `J ≥ 6` needs either `K₃` (`l1_tps.DEG = 8`) or exact spinning data.
* **The imaginary mixing would have mattered.** Treating it as real (`absmix`) breaks `n = 1` at `Δ_g = 300` and `100`. The `Re`
  prescription of §22.2 is therefore load-bearing. It is [D] as a statement about valid sum rules, but it changes the definition of
  `Ω^{A′}` (§R.3), which the referee should check.
* **In `c₂` language.** The asymptotic graviton-specific absolute corrections on `x ∈ [0.3, 0.7]` are `≤ 470` (`J ≤ 8`) and `≤ 660`
  (`J = 10`, `+−`), against the §19.3 tolerance `630 / 3.4·10³ / 2.7·10⁴`. They are comparable to the tolerance at `Δ_g = 100` and far
  inside it at `300` and `1000`, consistent with the `t*` above. The uniform-`c` tolerance is not the right comparison pointwise;
  the `2×2` `n`-test is.

### 22.5 Regraded table (supersedes §21.8 where they differ)

| statement | grade |
|---|---|
| ν-only normalisation `N̂₂(ν)` removed exactly by a heavy reference state at the same `ν` (`= §R.3` definition); reference-independent to `3·10⁻⁶`; helicity-independent to `3·10⁻⁸` | [D] (definition) + [NE] |
| §21.4 instability caused by truncating an exponentiated ν-only factor `e^{22i/ν+82/ν²}` (cumulant form fixes it) | [NE] |
| `β`, saddle helicity vectors real; `C₊₊C₊₋ > 0` | [NE] |
| residual imaginary part `i e₁ⁱ(x)/m` (diagonal `|e₁ⁱ| ≤ 2.3`; off-diagonal = the whole `O(1/m)` mixing of §21.4); real off-diagonal `O(1/m²)` | [NE] (`Δ = 400…6400`, `J = 4…10`) |
| `Ĉ₂ := Re[Ĉ^{CLS}₂/N̂₂]` is a valid sum rule with the same image, flat limit and light side | [D] (amends §R.3) |
| mechanism of the imaginary part (`q²` not a shadow eigenstructure) | [I] |
| `δc₂^{spin}(x, J, h)`, `J = 0…10`, stable `1/m²` extraction on `Δ = 400…3200`; asymptotic absolute size `≤ 470` (`J ≤ 8`), `≤ 660` (`J = 10`) on `x ∈ [0.3, 0.7]` | [NE] |
| same at `m ≈ 100` for `J ≥ 6`; `J = 12` below `Δ ≈ 400` | not controlled (Laplace series at `J²/m = O(1)`) |
| `n = 1` with the `Ĉ₂` correction (`2×2`, `Re` mixing) at `Δ_g = 1000` (`t* ≳ 14`) and `300` (`t* ≳ 3`), even `J ≤ 10` | [NE + M] |
| `n = 1` at `Δ_g = 100` | **open**: holds for `J ≤ 4` (`t* ≥ 8`); `J = 6–10` outside the controlled regime, and the extrapolation fails at `J = 10` |
| `Ĉ₃` correction, parity-odd sector (binding `+−`, `J = 5`), `J > 10`, the `k`-relative per-state constant | [ND] |
| **H4 for route A′** | **[I]**. Upgraded evidence at `Δ_g ≥ 300` (`Ĉ₂`, even `J ≤ 10`). New concrete risk at `Δ_g = 100`, `J ≥ 6`, `m ≈ Λ`. |

### 22.6 Hand-over (next steps, in order)

1. **`Δ_g = 100`, `J = 6…10`, `m ≈ 100–130`.** Either:
   * add `K₃`: `l1_tps.DEG = 8` (memory `N = 3432`), and `expand_pref` to `ε⁶`; then redo `l8_points` at `Δ = 100…400` only (the
     cumulant gains `p₆/p₀ − p₂p₄/p₀² + p₂³/(3p₀³)`); or
   * test the series at `Δ = 100` against an independent evaluation (exact spinning data are not available).
   Until then, `Δ_g = 100` is open.
2. **`Ĉ₃`** (§21.9 item 3). Copy `l6.build_spin` with `q² → o·q` (`g6`) and `f = −1`, and rerun `l8_refscan`/`l8_dc2` with the same
   machinery. Its reference subtraction fixes `N̂₃` the same way. Check whether `Ĉ₃` also needs `Re`. Then add its diagonal and
   off-diagonal to `l8_neargap` (the per-state `x → 0` constant becomes meaningful only for the `Ĉ₂`/`Ĉ₃` pair).
3. **Parity-odd `⟨TT𝒪⟩`** (odd `J`, binding `+−`, `J = 5`) and `J ≥ 12` (needs `Δ ≥ 800`).
4. Referee point: the `Re` amendment of §R.3, and whether `Ω^{A′}`'s light-side residue formula (§15) was derived for the complex or
   the real functional. It is real on the `T` block, so the bounds are unaffected, but state it.
5. Files: `l8_points.py` (cached points, orders separate), `l8_refscan.py` (+ `cache/l8_refscan_run.sh`), `l8_analysis.py` (basis,
   cumulant, `D_h`, mixing), `l8_dc2.py` (+ `summary()`), `l8_neargap.py`; caches `cache/l8_pts_J*.pkl`, `l8_dc2.pkl`, `l8_dc2_J*.txt`,
   `l8_summary.txt`, `l8_neargap_*.txt`. Added to `run_all.sh`.

## 23. M6 round 5: `K₃`, a linear-form chain, a CLS-D vs CHMRS shift term, and `Ĉ₃` (`l9_*.py`) — `Δ_g ≥ 300` confirmed; `Δ_g = 100` still open, for a sharper reason; hand-over in §23.7

*Written 18 Nov 2026 by the numerics agent continuing M6 after §22.6. Unrefereed, not committed. Updated incrementally.*

### 23.1 The next Laplace order `K₃` and its accuracy at `J²/m = O(1)` (`l9_tps.py`, `l9_scalar_k3.py`) [DV + NE]

* **Engine.** `l1_tps.DEG` is now read from `L1_DEG` (default 6, unchanged). With `L1_DEG=8`, `l9_tps.expand_pref6` gives the
  `ε⁶` term of `E[q e^A]` from the complete Bell polynomials `E₁…E₆`. Expectations of products up to degree 18 are exact: they are paired in
  a degree-9 product space by the Hermite inner product `E[ab] = Σ_α α! â_α b̂_α`. Self-tests (`python3 l9_tps.py`): `S₀…S₄` equal
  `l1_tps.expand_pref` to `3·10⁻¹⁶`; `E[e^{tξ}]` and `E[e^{cξ³}]` orders exact to `10⁻¹⁶`. Cost `≈ 5 s` per scalar point, `≈ 80–120 s` per
  graviton point (`J = 6–10`). The shadow matrix is also computed at `K₃` (fit residual `5·10⁻⁷` at `Δ = 100`, `J = 6`).
* **Accuracy check on the exact scalar skeleton** at the states that matter (`Δ_g = 100`, `m/Λ = 1.0–1.2`, i.e. `m ≈ 98–125`, `J ≤ 10`,
  `x ∈ [0.3, 0.7]`, 161 points; same fit as l2: `a_ν + b_J`; `cache/l9_scalar_k3_100.txt`). Residual `|Δc₂|` (c₂ units), rms / max:

| form | order | `J = 0` | `J = 4` | `J = 6` | `J = 8` | `J = 10` |
|---|---|---|---|---|---|---|
| linear (prefactor) | `K₁` | 22 / 79 | 52 / 123 | 40 / 80 | 63 / 130 | 160 / 609 |
| linear | `K₂` | 9 / 36 | 14 / 30 | 9 / 19 | 20 / 47 | 68 / 225 |
| linear | `K₃` | 1.6 / 6.7 | 3.4 / 8.9 | 2.7 / 5.2 | 5.4 / 10 | **19 / 61** |
| cumulant (log) | `K₃` | 27 / 109 | 1224 / 6587 | 59 / 261 | 158 / 596 | 136 / 468 |

* **Reading.** On `x ∈ [0.3, 0.7]` the Laplace series **in linear (prefactor) form converges at `m ≈ 100` up to `J = 10`**, by a factor `≈ 3`
  per order; after `K₃` the error is `≤ 61` in c₂ units (tolerance scale `630`). The **cumulant (log) form**, which doc 22 used for the
  graviton helicity elements, **does not converge** there (it is singular near zeros of `P_J` and degrades with order). So the doc-22
  `Δ = 100` blow-ups (`2·10⁴`, overflow at `J = 6`) were artefacts of the log form and of the `x → 0` extrapolation through points next to
  zeros of `P_J` (`J = 6`: `x = 0.2`; `J = 10`: `x = 0.1`), not of the expansion itself.

### 23.2 Linear-form chain; and a shift term between the CLS-D and CHMRS functionals (`l9_analysis.py`) [NE + D-sketch]

* **Chain.** Per element, the state series is multiplied order by order by `exp(−κ^{ref})` (the cumulant series of the heavy reference at the
  same `ν`), and truncated at `K`. This removes the ν-only factor exactly to the computed order, in linear form, and is regular through
  zeros of `D_++ = P_J`, `D_+− = d^J₄₄`. Per-state helicity constants `a_h` come from `Ĉ₂`; they are common to `Ĉ₃` (§23.3). Output:
  absolute corrections `Γ_{k,h}(x; m)` in c₂ units × LO shape, and the off-diagonal.
* **Check at `Δ_g = 300, 1000`.** With the same base as doc 22 the chain reproduces doc 22's verdict: `n = 1`, binding `+−`, `J = 5`.
* **New finding: the crystal computes the CLS-D functional, the base model and the light side use the CHMRS one.** In the linear chain the
  scalar crystal itself has `|Γ_s| ≈ 900` at `x ≈ 0.1–0.2` for `J = 10` (and `≈ 450` for `J = 6`), `m`-independent in c₂ units. The exact
  CHMRS remainder at the same `(m, J, x)` is `|m²ρ₂| ≤ 85`. The difference is CHMRS's magic-block relation (§20.2),
  `Ψ^{CLS}_ν ∝ g(ν)C₂(ν+i) + g(−ν)C₂(ν−i)`. Its `(g(ν) − g(−ν))/(g + g̃) ∝ i/ν` term gives `≈ −(α/x)∂_x(LO) ≈ 2αJ(J+1)` near `x → 0`.
  * Subtracting this shift from the LO density (exact complex shift `x → x ± i/m`, scalar `g` of `l2`) brings the scalar crystal to
    `|Γ_s| ≲ 100` at `m ≥ 800` [NE]. At `m ≤ 400` it does not work below `x ≈ 0.3` (residual `∝ 1/m²` in c₂ units, `±1.5·10³` at
    `m = 98`). Reason: the truncated crystal does not reproduce the non-perturbative `1/ν` structure of `g(ν)` at `ν ≲ 30`.
  * **ψ-integrated size.** The shift changes the per-state action by up to `1.8 / 0.27 / 0.033 N_s` (`J = 10`; `Δ_g = 100 / 300 / 1000`).
    This is the difference between `∫ψ₂ Ψ^{CLS}_ν` and `∫ψ₂ C^{CHMRS}_ν`: two different valid functionals (a contour shift of `ν` turns one
    into the other with reweighted `ψ` plus a pole term from `g`'s `1/ν`).
  * **Consequence.** §R defines `Ω^{A′}` with `Ĉ^{CLS}`, but the near-gap base (exact skeleton `ρ_k`), the light side (§15, CHMRS residue
    formula) and doc 22 (graviton ÷ scalar crystal × exact CHMRS scalar) all use, implicitly, the CHMRS-type functional `∫ψ C_ν`. That is a
    consistent choice for the scalar parts. For the graviton it requires the spinning `Ĉ^{CLS} ↔ C`-type relation, which is not
    derived: the scalar `g` is assumed [I]; the graviton's ν-only factor differs from the scalar's, §22.1. **§R must state which functional
    `Ω^{A′}` is.** The CHMRS-type choice (`∫ψ C`, de-shifted crystal) is the one used below.

### 23.3 `Ĉ₃` integrand built; leading-order check passed for both odd structures (`l9_points.py`, kinds `C3H`, `C3V`) [NE]

* **Integrand.** As `Ĉ₂` (§21.3), with:
  * dual `q²/V⁴ → o·q/V⁴`, `f = −1`;
  * celestial `j′ = 1` = the `j′ = 0` one times a locally constant sign. Both `S⁰` points give the same value, since `o` and the sign are
    both odd.
  * `ε(·)` is the 5×5 embedding determinant by cofactors on TPS (checked against `numpy.det` to `4·10⁻¹⁶`).
* **The two crystal-equivalent completions:**
  * `o_H = ε(Z₂,X₂,X₄,X₀,Z₀)H₄₀/√(X₂₄X₂₀X₄₀³) + ε(Z₄,X₄,X₂,X₀,Z₀)H₂₀/√(X₂₄X₄₀X₂₀³)` (§21.6);
  * `o_V = V₀,₄₂[ε(Z₂,…)V₄,₀₂/√(X₂₄X₂₀X₄₀³) − V₂,₀₄ ε(Z₄,X₄,X₂,X₀,Z₀)/√(X₂₄X₄₀X₂₀³)]`. Its homogeneity was re-derived; it has the same
    degrees as `q`.
* **LO check** (`J = 4`, `Δ = 800` and `1600` at equal `ν`, order 0, helicity basis of `Ĉ₂`):
  * `Ĉ₃` is diagonal in the `Ĉ₂` helicity basis: off/diag `≤ 9·10⁻⁵`.
  * `r_{+−}/r_{++} = −1.000` (to `≤ 2·10⁻³`, both structures), i.e. the `σ_h` pattern.
  * `(2−y)·r_{++} = −2i x²` to `0.3 %`, with `r_h := Ĉ₃/Ĉ₂`. The `x²` at fixed `ν` scales as `m⁻²` between `Δ = 800` and `1600`
    (ratio `4.000`). So the per-state constants are common to `Ĉ₂` and `Ĉ₃`, and `Ĉ₃` carries the relative `u m²` of §R.6.
  * `o_H` and `o_V` agree at LO; their `O(1/m)` imaginary parts differ by a factor `2` (`r_{+−}/r_{++}` phase `0.8 %` vs `1.6 %` at `x = 0.3`).
* **Which structure defines the functional.** Both are valid `j′ = 1`, `f = −1` light-ray functionals with the same crystal image and flat
  limit, so the choice is a design freedom. **`o_H` is adopted** (the first, §21.6 form). Its `O(1/m²)` difference from `o_V` was
  not computed (see §23.7).

### 23.5 Independent assessment of the `Re` prescription (§22.2; amends §R.3) [D, with three conditions]

**Verdict.** `Ĉ₂ := Re[Ĉ^{CLS}₂/N̂₂(ν)]` is a legitimate sum rule with the same flat limit and the same Regge convergence, provided
conditions (a)–(c) hold. (a) and (b) hold for the parity-even sector computed so far. (c) is [D-sketch] at the same level as §15.

**Why it is a valid sum rule.** For each `ν`, `Σ_O Σ_{ab} λ_a λ_b Ĉ_ν[G_{O,ab}] + (light) = 0` is a complex identity. Dividing by any
state-independent `N̂(ν)` and taking the real part gives a real identity, provided `λ_aλ_b` are real. `Re[Ĉ_ν/N̂]` equals
`½[Ĉ_ν/N̂ + (Ĉ_ν/N̂)*]`. On real `G`, `(Ĉ_ν)*` is the light-ray functional with conjugated kernel, i.e. the `−ν` celestial weight with the
(real-coefficient) structure `q²`. That is again a CLS-D functional, with the same `J₀ < 2` convergence. `ψ` is real, so `Re` commutes with
`∫ψ`. The flat limit is unchanged: the leading crystal action is real after `N̂` (§22.2), and the imaginary part is `O(1/m)`.

**Conditions.**
* **(a) Real OPE-coefficient basis.** `λ_aλ_b` must be real in the structure basis in which `Re` is taken. This holds for the parity-even
  2×2: real conserved basis, real `β`, real saddle helicity vectors (§22.1).
  * For the parity-odd structure the natural `ε`-form may carry a factor `i` in Lorentzian signature. Then the real-coefficient basis is
    `i·ε`, and in the even–odd mixed entries "Re" becomes "Im" of the naive entries. **This must be fixed before the odd sector
    (binding `+−`, `J = 5`) is computed.** (`Ĉ₃`'s own `−i` relative to `Ĉ₂`, §23.3, is ν-only and absorbed in `N̂₃`.)
* **(b) `N̂` state-independent.** `N̂₂(ν)` is defined as the `m → ∞` limit at fixed `ν` (§R.3), so it is ν-only by construction.
  * Numerically it is a single reference `Δ_ref = 2·10⁴` (reference-independent to `3·10⁻⁶`; helicity-independent, relative phase
    `≤ 2·10⁻⁶π`).
  * The per-state helicity constants `c_h(Δ)` are real (ratios of real structure normalisations), hence a congruence. They do not mix `Re`
    and `Im`.
  * Subtlety: the **phase** of `N̂` is a convention ("heavy limit real positive"). `Re[e^{iθ(ν)}Ĉ/N̂]` is a *different* valid functional.
    Heavy and light sides must use the same `θ`.
* **(c) Light side with the same phase.** The bound uses `−Ĉ₂|_E = 16πG/(ν²+9/4)`, real, from the CHMRS residue formula in the
  heavy-Regge normalisation, so `θ = 0` on both sides [D-sketch, as §15].
  * If the `T`-block action carried a residual ν-dependent phase `φ(ν)` relative to the heavy limit, `Re` would multiply the graviton pole
    by `cos φ(ν)`. That changes `L_A` and the constants `A`, `B_AdS`, **not positivity**.
  * §23.2 adds a sharper version of this point: heavy side, light side and `N̂` must refer to the **same** functional, CLS-D (`Ψ_ν`) or CHMRS
    (`C_ν`). They differ by the shift relation, whose `g(ν)` is complex. `Re` of one is not `Re` of the other.

**Is it load-bearing?** Yes. Treating the imaginary mixing as real breaks `n = 1` at `Δ_g = 300` (§22.4 `absmix`). The prescription is
legitimate given (a)–(c), but it is part of the definition of `Ω^{A′}` and must be stated in §R together with the CLS/CHMRS choice.

### 23.4 `Δ_g = 100`: the graviton series at `J = 10`, `m ≈ 98` is not converged at `K₃`; `n` stays open (`l9_run.py`, `l9_neargap.py`) [NE]

* **`K₃` graviton points** (`Ĉ₂`, `Δ = 100`, `x = 0.05…0.8` plus `0.15`, `0.25`, each with its `Δ_ref = 2·10⁴` reference at `K₃`).
  * `J = 10`: done.
  * `J = 8`, `6`: running at hand-over (`cache/l9_run_C2_a.log`; resumable: `L1_DEG=8 L9_DS=100 python3 l9_run.py C2 10 8 6`).
* **Convergence, normalisation-free** (change `K₂ → K₃` of `m² Re log M̃_hh`, relative to `x = 0.5`; `J = 10`, `m = 98`, points away from zeros
  of `P₁₀`):
  * `++`: `+689 / +250 / +478 / +607` at `x = 0.3 / 0.6 / 0.7 / 0.8`;
  * `+−`: `+716 / +406 / −332 / +359`;
  * `x ≤ 0.25`: `10³–10⁴`.
  For comparison, on the scalar at the same `(m, J)` the `K₂ → K₃` change is `≈ 70` rms (§23.1). **The spinning series moves by roughly the
  §19.3 tolerance (`≈ 630`) at the last computed order, so `K₃` does not control `J = 10` at `m ≈ 100`.**
* **The per-state normalisation at `m ≈ 100` cannot be taken from small `x`.** The fitted Regge constant `a_h` (from `ν ≥ 15`, `x ≤ 0.25`)
  changes sign between `K₂` and `K₃` (`J = 10`). At `ν ≲ 30` the truncated crystal misses the non-perturbative `1/ν` structure of `g(ν)`
  (§23.2).
  * With `Ĉ₂` alone this constant is load-bearing: it fixes the k-relative normalisation against the uncorrected `Ĉ₃` base.
  * With `Ĉ₂ + Ĉ₃` both corrected it is a congruence and drops out of `n`, but the `x < 0.3` values still enter through their `ψ` weight.
* **`n` with what is controlled** (`Ĉ₂`, linear chain, de-shifted, `K₂` rows `Δ = 100…3200`, direct interpolation in `1/m`;
  `cache/l9_neargap_C2_K2_deshift.txt`):

| `Δ_g` | `n` | min `λ/N_s` (binding) | `t*` | `J = 10` max `|ΔF_hh|/N_s` | status |
|---|---|---|---|---|---|
| 1000 | 1 | 0.2021 (`+−`, `J = 5`, uncovered) | (−30.6, +38.8) | 0.008 | consistent with §22.4 (`t* ≥ 14`) |
| 300 | 1 | 0.2054 (`+−`, `J = 5`, uncovered) | (−10.7, +8.4) | 0.036 | consistent with §22.4 (`t* ≥ 3`) |
| 100 | 1 | 0.1552 (`J = 4`, `m/Λ = 1.04`) | (−0.8, +1.6) | 0.88 | **not controlled** (`J ≥ 6` at `m ≈ Λ`: series and normalisation, above) |

  For the record, without the de-shift (literal CLS-D functional, §23.2): `J = 10` `|ΔF|/N_s` up to `1.8 / 0.27 / 0.033`, `t*` at
  `Δ_g = 100` `(−0.3, +11.5)`, `n = 1` at all three gaps.
* **Verdict for `Δ_g = 100`: open.** It is not decided either way.
  * The crystal (Laplace) expansion of the spinning functional does not converge fast enough at `J²/m ≈ 1` for `J ≥ 8–10`.
  * Its small-`x` region (`ν ≲ 30`) is outside the expansion at `m ≈ 100` for every `J`.
  * No computed number breaks `n = 1`. The doc-22 `−2.70` at `J = 10` came from the log form and the `x → 0` extrapolation, and is
    superseded.
  * A decision needs a non-perturbative evaluation at those few states: direct numerical integration of the 7-dimensional CLS-D integrand
    on the steepest-descent contour, or exact spinning Mellin data. Alternatively, the honest scope of Result D at `Δ_g = 100` is the
    `J ≤ 4` sector.

**Combined `Ĉ₂ + Ĉ₃` per-state `2×2`** (`o_H`, `K₂`, de-shifted, rows `Δ = 100 / 400 / 1600` for both functionals;
`cache/l9_neargap_C23H_K2_100_300_1000.txt`). With both corrected, the per-state normalisation is a congruence (§23.2).

| `Δ_g` | `n` | min `λ/N_s` (binding) | `t*` | `J = 10` max `|ΔF_hh|/N_s`, `|ΔF_+−|/√(NN)` |
|---|---|---|---|---|
| 1000 | **1** | 0.2021 (`+−`, `J = 5`, uncovered) | (−37.0, +38.9) | 0.006, 0.020 |
| 300 | **1** | 0.2054 (`+−`, `J = 5`, uncovered) | (−5.1, +5.0) | 0.041, 0.184 |
| 100 | 1 (not controlled) | 0.1855 (`J = 4`, `m/Λ = 1.04`) | (−0.6, +1.6) | 1.18, 0.64 |

### 23.6 Regraded table (supersedes §22.5 where they differ)

| statement | grade |
|---|---|
| `K₃` engine (Hermite-exact Gaussian moments, `ε⁶`) | [DV] (self-tests to `10⁻¹⁶`) |
| scalar crystal Laplace series at `m ≈ 100`, `J ≤ 10`, `x ∈ [0.3, 0.7]`: converges in linear form (`K₃` error ≤ 61 c₂ units); log form does not | [NE] (exact skeleton, 161 points) |
| doc-22 `Δ = 100` blow-ups and the `−2.70` at `J = 10` | artefacts (log form; `x → 0` extrapolation through zeros of `P_J`); **withdrawn** |
| crystal = CLS-D `Ψ_ν`; base model and light side = CHMRS `C_ν`; they differ by the magic-block shift, up to `1.8 / 0.27 / 0.03 N_s` at `J = 10` | [NE] (scalar: de-shift validated at `m ≥ 800`); **§R must fix which functional `Ω^{A′}` is** |
| graviton uses the scalar `g(ν)` in the de-shift | [I] |
| `Ĉ₃` integrand (`o_H`, `o_V`); LO: diagonal in the `Ĉ₂` helicity basis, `σ_h` pattern, relative `u m²`, common per-state constants | [NE] |
| `o_H` vs `o_V` at `O(1/m²)` | [ND] (design freedom; `o_H` adopted) |
| `Re` prescription | [D] given (a) real OPE basis, (b) ν-only `N̂` with a fixed phase convention, (c) the light side in the same convention and the same (CLS-D or CHMRS) functional (§23.5) |
| `n = 1` at `Δ_g = 1000` / `300`, `Ĉ₂ + Ĉ₃` corrected (even `J ≤ 10`), `t* ≈ ±37` / `±5` | [NE + M] (`K₂`, CHMRS-type functional) |
| `n` at `Δ_g = 100` | **open**: graviton series at `J = 10`, `m ≈ 98` moves by `≈ 250–720` c₂ units from `K₂` to `K₃`; normalisation and `x < 0.3` (`ν ≲ 30`) are outside the expansion. No computed number breaks `n = 1`. `J ≤ 4` fine (§22.4) |
| odd `J` / parity-odd sector (binding `+−`, `J = 5`) | [ND]; handed to a separate agent |
| **H4 for route A′** | **[I]**. Evidence at `Δ_g ≥ 300` strengthened (`Ĉ₂ + Ĉ₃`). `Δ_g = 100` undecided by the crystal method |

### 23.7 Hand-over

1. **Running at hand-over** (2 cores, nohup, resumable):
   * `l9_run.py C2 10 8 6` at `K₃`, `Δ = 100` (`J = 8` refs, then `J = 6`; log `cache/l9_run_C2_a.log`);
   * `l9_run.py C3H` at `K₂`, rows `Δ = 200, 800, 3200` (log `cache/l9_run_c3k2.log`).
   * When done, run K₂-vs-K₃ comparisons via `l9_analysis.grid(J, 3, kinds, Ds=[100.0])`, and `l9_neargap.py C23H_K2` (and `C2_K3`).
2. **Deciding `Δ_g = 100` needs a non-perturbative evaluation at a few states** (`J = 6, 8, 10`, `m/Λ = 1.0–1.3`; the series is not
   controlled there).
   * Recommended: direct numerical integration of the 7-variable CLS-D integrand on the real `ξ`-contour `z = z* + Lξ` (`l1_crystal`'s
     `L`). Use `l6`/`l9_points.build_spin` evaluated pointwise rather than as TPS, with the steepest-descent Gaussian as importance
     weight (QMC, about `10⁵–10⁶` points per `(state, ν)`).
   * Validate it first on the scalar against e1 at the same `(m, J)`, including `x < 0.3`.
   * Both `Ĉ₂` and `Ĉ₃` are needed (so that the normalisation is a congruence), and the CHMRS de-shift at the same `ν`.
3. **Fix the functional (§R).** Decide CLS-D vs CHMRS-type (§23.2) and state the `Re` convention (§23.5). If CLS-D is kept, the base model
   (exact skeleton `ρ_k`) and the light side (§15) must be re-derived through the shift relation.
4. `o_V` vs `o_H` at `O(1/m²)` (kind `C3V` in `l9_points`).
5. Files: `l1_tps.py` (env `L1_DEG`, default unchanged), `l9_tps.py`, `l9_scalar_k3.py`, `l9_points.py` (kinds `C2`/`C3H`/`C3V`, `K₂`/`K₃`),
   `l9_run.py`, `l9_analysis.py` (`L9_DESHIFT`, `L9_NUFIT`), `l9_neargap.py`. Caches: `cache/l9_*`. Added to `run_all.sh`.

## 24. M6 round 6: the definition of `Ω^{A′}` fixed (CHMRS-type), and a non-perturbative evaluation of the crystal integral (`n1_*.py`) — `x ≥ 0.25` decided, `Δ_g = 100` still open through `x < 0.25`; hand-over in §24.7

*Written 18 Nov 2026 by the numerics agent continuing after §23.7. Unrefereed, not committed. Updated incrementally.*

### 24.1 Definition: why the literal CLS-D functional cannot be `Ω^{A′}`, and what replaces it (`n1_definition.py`, `cache/n1_definition.txt`) [D given Q + NE]

* **The CHMRS magic-block coefficient has a real zero.** With `g(ν) = ((2−d)/2+iν)/(iν)·(γ₇(ν+i)/γ₇(ν))²` (l2, CHMRS), `g(−ν) = g(ν)*` for
  real `ν`, and `g + g̃ = 2 Re g` changes sign at **`ν₀ = 5.9791`** (`+2.18` at `ν = 2`, `+0.43` at 5, `−0.36` at 7, `→ −2`).
  `ν·Im g/(g+g̃) → 23/4` at large `ν` (validated at `ν ≥ 30` by §20: it is the `O(1/m²)` state dependence l2 reproduces to `0.2` c₂).
* **Consequence for the literal CLS-D functional [D given the identity].** At fixed `ν`, `m → ∞`, `C(ν ± i) → C_heavy`, so the heavy-limit
  normalisation of `Ψ^{CLS}_{2,ν} = N(ν)[g C(ν+i) + g̃ C(ν−i)]` is `N(ν)(g+g̃)C_heavy`. `N` cannot have a pole at `ν₀` (else `Ψ^{CLS}` itself,
  a convergent integral, would be infinite on every block with `C(ν₀+i) ≠ C(ν₀−i)`). So **the heavy limit of `Ψ^{CLS}` vanishes at `ν₀`**, and
  `Ψ^{CLS}/N̂(ν)` of §R.3 has a pole at `p₀ = ν₀/Λ = 0.0607 / 0.0200 / 0.0060` (`Δ_g = 100 / 300 / 1000`). With any regular normalisation
  instead, the heavy action of the k = 2 part changes sign at `p₀` (positive Regge-limit action lost for `p < p₀`). Either way the literal
  CLS-D `∫ψ₂Ĉ^{CLS}₂` is not an admissible `Ω^{A′}` at finite `Δ_g`. (Caveat: rests on the exact form of `g`, quoted from CHMRS and validated only
  through its large-`ν` asymptotics.)
* **Why §R.5's `η`-inversion does not repair it.** The identity implies the CLS-D Regge kernel `[gP_{l−1} + g̃P_{l+1}]/(g+g̃)`, `l = −½+iν`, i.e.
  `ηP_l + i(σ − 1/(4ν))(P_{l−1} − P_{l+1})` with `σ = Im g/(g+g̃) ≈ 23/(4ν)`; §R.5 inverted only the `ηP_l` part. The remainder,
  `≈ −(23/(2ν))∂_ν`, is exactly the CLS→CHMRS "shift term" of §23.2 (up to `1.8 N_s` at `J = 10`, `Δ_g = 100`).
* **Definition adopted (amends §R.1–§R.5).** `Ω^{A′} := ∫₀¹dp ψ₂(p) Λ⁴ Re Ĉ_{2,Λp} + ∫₀¹dp ψ₃(p) ½Λ⁶ Re Ĉ_{3,Λp}`, with **`Ĉ_{2,ν}` the CHMRS-type
  functional**: the linear combination of CLS-D functionals `Ψ_{2,ν′}` (`ν′` on the shifted lines `Im ν′ = ±1` plus the pole term of `g` at
  `ν = 0`) that inverts the magic-block identity, i.e. **the same combination for every correlator**; on scalar correlators it equals CHMRS's
  `C_{2,ν}` (heavy limit exactly `2P_{−½+iν}(cosh β)/m⁴`), on `⟨TTTT⟩` it is that combination of the spinning `Ψ_{2,ν′}[q²]`. Equivalently:
  §R.5's `η`-inversion with the exact kernel `gP_{l−1} + g̃P_{l+1}` in place of `½ηP_l`. `Ĉ₃` likewise (CHMRS-type `C₃`, as in the base
  model; its own shift relation is not derived — the `k = 3` analogue of the statement below is [I]).
  * **What this buys.** The base model (exact CHMRS skeleton `ρ_k`), Theorem A (`CERTIFIED_NEARGAP_A.md`) and the light side (§15, CHMRS residue
    formula for `C_{2,ν}`) are all statements about this functional — no re-derivation needed.
  * **The graviton heavy action is then the scalar-`g` de-shift of the crystal, by definition** (it is the same combination of `Ψ_ν′`). The
    previous [I] ("graviton uses the scalar `g`", §23.6) is **removed from the heavy side** and moves to the light side: that the graviton-
    pole residue of this combination on `⟨TTTT⟩` is `16πG/(ν²+9/4)` (it is if the spinning `Ψ_ν` has the same Regge kernel `gP_{l−1}+g̃P_{l+1}`
    up to a ν-only factor; [I], same level as §15's spinning carry-over). That affects `A`, `B_AdS` only, **not `n`**.
  * **Residual heavy-side issue.** The de-shift used in l9/l8 (`shift()`, applied to the LO density) is the local form of the inverse; it is
    ill-conditioned where `g + g̃ → 0`, i.e. for `x ≲ 0.1` at `m ≈ 100` (`ν ≈ ν₀`). There the crystal does not determine `Ĉ₂` without the
    exact (non-local in `ν`) inverse.
* **Effect on the near-gap numbers.** Unchanged for the base model, Theorem A and §23's `C23H` "de-shifted" rows (they already use this
  object). The literal-CLS-D numbers of §23.4 ("without the de-shift: `|ΔF|/N_s` up to `1.8 / 0.27 / 0.033`") describe a functional that is
  not admissible; they are withdrawn as an option.

### 24.2 Non-perturbative evaluation of the crystal integral: method (`n1_pw.py`, `n1_direct.py`) [NE]

* **Point-wise backend.** `n1_pw` re-executes `l1_crystal`, `l5_spin`, `l6_spin_action`, `l9_points` with the TPS engine replaced by point
  arrays, so the **exact 7-variable integrand** (scalar `build(pref=True)`; graviton `build_spin`, kinds `C2`, `C3H`) is evaluated at
  arbitrary complex points. Logs and fractional powers are continued from the saddle along the ray (two-step path continuation,
  levels `ξ/2, ξ`); per-point winding is recorded.
* **Contour.** `u_x, v_x, u_y, v_y → c_i e^{w_i}` (log coordinates: the linear contour reaches `u = 0` at `4σ`), `x̃₀ = tanh w₀`
  (`x̃₀ ∈ (−1,1)`), others linear; the centre `w_c` is the saddle of `Φ = λF + G + (J/2)log(Z²W²) + log Jacobian` (the `J`-growing part of
  the integrand; `P_J` and all graviton structures stay as bounded prefactors), `L` from `−∂²Φ`. Real `ξ`, `w = w_c + Lξ`.
* **Masks.** Points with `Re(Φ − Φ_c) > 2` (past a branch wall: there the integrand jumps to `e^{+45…+105}`) or winding `> 3` are dropped;
  inside `|ξ| ≲ 6` the integrand is `≤ e^{−8}` of its peak in the worst direction, so the walls cost `≲ 10⁻⁴` (dropped fraction of the
  sampling measure `≤ 10⁻³`, mostly at `x ≤ 0.25`).
* **Integration.** Randomised QMC (scrambled Sobol, `2^{13–14}` points × 8 scramblings, antithetic), Gaussian proposal of scale `1.2`;
  error = replica spread. (Tensor and Smolyak Gauss–Hermite rules fail: corner nodes cross the walls; Smolyak converges only at
  `10⁻³`.) Cost ≈ 8 s (scalar) / 55 s (graviton 2×2 + scalar on the same samples) per point.
* **Structural finding: five of the seven directions are `ν`-soft.** In log coordinates the Hessian of `λF` has two eigenvalues `∝ Δ`
  (`x⊥`, `y⊥`) and five `∝ ν` (`u_x, v_x, u_y, v_y, x̃₀`; e.g. `ν = 10`, `Δ = 400`: `±4.98i`, `0.50 ± 4.96i`, `−51 + 13i`). The crystal
  expansion at fixed `x` is therefore an expansion in `1/ν = 1/(xm)` in those directions, and the slow factors' gradients (`O(1)`, `J`-
  enhanced) shift the saddle by `O(J/ν)`. This is why the Laplace series fails for `x ≲ 0.3` at `m ≈ 100` for every `J` (§23.4) and why the
  "ν-only factor" `exp(22i/ν + 82/ν²)` is non-perturbative.
* **Range of the method.** Converges for `ν ≳ 20–25` (relative error `10⁻³–10⁻²` per element at `2^{14}·16` points). **Fails for
  `ν ≲ 15`** at every `Δ` tested (100, 400, 2·10⁴): the real-`ξ` contour through the `ν`-soft directions crosses walls within `2–3σ` (masked
  fraction → 100 %), and adaptive Gaussian proposals do not help. A non-perturbative value at `ν ≲ 15` needs a curved (steepest-descent)
  contour in the five soft directions; not done.

### 24.3 Scalar validation against the exact skeleton (`n1_scalar_validate.py`, `cache/n1_scalar_validate.txt`) [NE]

Same points and fit as l9_scalar_k3 (`Δ_g = 100`, `m/Λ = 1.0–1.2`, `J = 0…10`, exact e1 Mellin data, `X = gC(ν+i) + g̃C(ν−i)`, fit
`a_ν + b_J`), plus `x ∈ [0.2, 0.3) ∪ (0.7, 0.85]` with `ν ≥ 20`. **Status: see the numbers block at the end of this subsection (filled when
the run finished; if absent, the run did not finish — resumable, `cache/n1_chain3.sh`).** Note the comparison object `X` itself is built
from the LO density with third-order shifts (`O(m⁻⁴)`, §20.2), so it is a test of the direct integral to `≈ 10` c₂ units, not better.

### 24.4 Graviton `Ĉ₂`, `Ĉ₃` (o_H): direct vs Laplace on the load-bearing window (`n1_run.py`, `n1_ratio_check.py`, `n1_analysis.py`) [NE]

* **Points.** l9's `Δ = 100` row (`m = 97.9–98.4`, the row that carries ≥ 80 % of the `1/m` interpolation weight for the `m/Λ ≈ 1.04–1.3`
  states), `x = 0.25, 0.3, 0.4, …, 0.8`, `J = 4, 6, 8, 10`, both kinds, plus the heavy references `Δ_ref = 2·10⁴` at the same `ν`
  (`2^{14}·16` resp. `2^{12}·16` RQMC points; `cache/n1_pts_*.pkl`). The graviton 2×2 and the scalar are integrated on the same samples.
* **Normalisation-free check** (state/reference ratio of the helicity-diagonal elements, `m²(v_K/v_direct − 1)`, c₂ units relative to the
  local value; `cache/n1_ratio_check.txt`). `K₂` misses the direct value by `≈ 50–900` (up to `1.5·10³` at `x = 0.8`, `J = 6`), `K₃` (C2,
  `J = 8, 10`) by `≈ 20–700`; RQMC errors `5–60` (`≤ 275` next to zeros of `P_J`). **So neither Laplace order is accurate pointwise at
  `m ≈ 98` for `J ≥ 6`, as §23.4 suspected, but the error is relative to small local values:**
* **In the units that enter `F`** (absolute `Γ`, i.e. c₂ × LO shape, l9 normalisation): `max |Γ_direct − Γ_K₂|` on `x ≥ 0.25` is
  `40 / 88 / 117 / 77` (`J = 4 / 6 / 8 / 10`, `C2 ++`), `≤ 103` for `C3H`, `≤ 47` for `+−`, `≤ 1` off-diagonal — against `|Γ|` up to
  `1.5–2.1·10³` there and `2–3.9·10³` at `x < 0.25`, and the §19.3 tolerance scale `630`. The `P_J`-oscillation makes the LO shapes small on
  `x ≥ 0.25` for `J ≥ 8`; the weight sits at `x < 0.25`.

### 24.5 `n` at `Δ_g = 100` with the direct values (`cache/n1_analysis.txt`) [NE + M]

`C23H_K2` (both functionals corrected, de-shifted = the §24.1 definition, `Re`), rows `Δ = 100 / 400 / 1600` (+200 for `J = 10`), even
`J ≤ 10` covered, odd `J` and `J > 10` unchanged (parity-odd sector: separate agent):

| rows | `n` | min `λ/N_s` (binding) | `J = 4 / 6 / 8 / 10` min `λ/N_s` | `t*` |
|---|---|---|---|---|
| l9 (Laplace `K₂`) | 1 | 0.1855 (`J = 4`, `m/Λ = 1.04`) | 0.1855 / 0.2146 / 0.2304 / 0.4886 | (−0.7, +1.6) |
| **direct on `x ≥ 0.25` (`Δ = 100` row), `K₂` below** | **1** | **0.1770** (`J = 4`, `m/Λ = 1.04`) | 0.1770 / 0.2046 / 0.2249 / 0.4886 | (−0.7, +1.6) |

* Replacing the Laplace values by non-perturbative ones on `x ≥ 0.25` moves `λ_min/N_s` by `≤ 0.010` (`J = 6`). RQMC error on these
  numbers: `≲ 0.003 N_s` (from the replica spreads, `|ΔΓ| ≲ 30` c₂ on the relevant shapes).
* **What is still not controlled: `x < 0.25` (`ν ≲ 25` at `m ≈ 100`).** There the corrections enter with 35–50 % of `Σ|ψ density|` for
  `J ≥ 6`, the Laplace chain gives `|ΔF_hh|/N_s` up to `1.0` (`J = 10`), and §24.2's direct method does not reach `ν ≲ 15`. `t* = (−0.7, +1.6)`
  says that a 70 % error of opposite sign in the total correction would break `n = 1` (at `J = 4`, `m/Λ = 1.04`).
* **Verdict: `n = 1` at `Δ_g = 100` is supported, not decided [NE + M + I]:** it holds with every computed number, and the controlled
  part (`x ≥ 0.25`) is now non-perturbative; the uncontrolled part is the `ν`-soft region `x < 0.25`, where the margin is `t* ≈ −0.7`.
  `Δ_g = 300, 1000`: unchanged (§23.4, `n = 1`, `t* ≈ ±5 / ±37`).

### 24.6 Regraded table (supersedes §23.6 where they differ)

| statement | grade |
|---|---|
| `g + g̃` (CHMRS magic block, `a = 7`, `d = 3`) vanishes at `ν₀ = 5.979`; `ν Im g/(g+g̃) → 23/4` | [DV] (closed form, numerics) |
| heavy-limit normalisation of the literal CLS-D `Ψ_{2,ν}` vanishes at `ν₀`; `Ψ/N̂` (§R.3) has a pole at `p₀ = 0.061 / 0.020 / 0.006` | [D] given the CHMRS identity (only its large-`ν` form is validated, §20) |
| `Ω^{A′}` := CHMRS-type `Ĉ₂` (inverse magic-block combination of CLS-D functionals, the same for every correlator), `Ĉ₃` CHMRS-type | **definition (amends §R.3–§R.5)**; base model, Theorem A, §15 light side refer to it |
| graviton heavy action of that functional = scalar-`g` de-shift of the crystal | [D] (by definition; was [I] in §23.6). Local de-shift ill-conditioned at `x ≲ 0.1` (`ν ≈ ν₀`) |
| graviton light side of that functional `= 16πG/(ν²+9/4)` | [I] (needs the spinning `Ψ` to share the scalar Regge kernel `gP_{l−1}+g̃P_{l+1}`); affects `A`, `B_AdS`, not `n` |
| `Ĉ₃` shift relation | [I] (not derived) |
| crystal Hessian: 2 directions `∝ Δ`, 5 `∝ ν` (log coordinates) — the expansion parameter at small `x` is `1/ν` | [NE] |
| direct RQMC crystal integral (re-centred log contour, path-continued logs, wall masks) | [NE]: agrees with Laplace at `Δ = 1000` to `6·10⁻⁵`; converges for `ν ≳ 20–25`; **fails for `ν ≲ 15`** |
| scalar validation of the direct integral vs exact skeleton | see §24.3 block |
| `K₂`/`K₃` vs direct at `m ≈ 98`, `x ∈ [0.25, 0.8]`, `J = 4–10`: absolute `|ΔΓ| ≤ 117` c₂ (`K₂`) | [NE] |
| `n = 1` at `Δ_g = 100` with direct values on `x ≥ 0.25`: min `λ/N_s = 0.177`, `t* = (−0.7, +1.6)` | [NE + M]; **`x < 0.25` uncontrolled → `n` at `Δ_g = 100` remains [I]** |
| `n = 1` at `Δ_g = 300 / 1000` | [NE + M] unchanged (§23) |
| **H4 for route A′** | **[I]**; at `Δ_g = 100` reduced to the `ν`-soft region `x < 0.25` |

### 24.7 Hand-over

1. **The remaining obstacle is `ν ≲ 20` (`x < 0.25` at `m ≈ 100`), and it is a contour problem, not a truncation problem.** Five of seven
   directions are `ν`-soft. Next step: a curved steepest-descent contour in the soft block (`u_x, v_x, u_y, v_y, x̃₀`) — e.g. integrate the
   two `Δ`-stiff directions by Gauss–Hermite at each soft point and treat the 5-dim soft integral with a Picard–Lefschetz flow (gradient
   flow of `Re(νf)` from the real `w`-plane for a fixed time), which removes the walls. Validate on the scalar at `ν = 5–20` against exact
   Mellin data (the e1 nodes exist down to `ν ≈ 0.03`).
2. **Alternatively decide `n` without small-`x` crystal values**: since both `k` are corrected, the per-state constant is a congruence; what
   matters at `x < 0.25` is the graviton-specific `x`-dependence. A graviton/scalar ratio at fixed `(ν, Δ, J)` on common samples has much
   smaller errors (the soft structure is shared) — worth testing at `ν = 15–25` first (`n1_run` stores the scalar on the same samples, `Es`).
3. **Definition**: amend §R.3–§R.5 as in §24.1 (CHMRS-type `Ĉ₂`; `η`-inversion with the exact kernel `gP_{l−1}+g̃P_{l+1}`). The spinning
   light-side statement (16πG/(ν²+9/4) for the inverse combination) is the new [I]; `k = 3` shift relation open.
4. Files: `n1_pw.py` (point-wise backend), `n1_direct.py` (contour, masks, RQMC/Smolyak/GH rules), `n1_run.py` (graviton points),
   `n1_scalar_validate.py`, `n1_ratio_check.py`, `n1_analysis.py` (direct rows into l9's chain and l9_neargap), `n1_definition.py`;
   caches `cache/n1_*` (pickles `n1_pts_<kind>_J<J>.pkl`, `n1_scal_100.pkl`). Added to `run_all.sh`. Points stored before the
   log-Jacobian constant fix lack the key `ljc` and are corrected in `n1_analysis.Hrep`.
