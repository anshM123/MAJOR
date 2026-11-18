# H0: the spinning (`⟨TTTT⟩`) functional — route decision and first milestone

Status: 18 Nov 2026, derivation/numerics agent. M1 (§0–§8) committed (ce46a54); M2/M3 (§9–§11) committed (7211186); §R, M2b/M3b/M4′ (§12–§14) committed (78eca6a); **M5′/M2c (§15–§16) unrefereed, not committed.** Targets hypothesis H0 of `RESULT_D.md` (with H2, H4, H14,
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

Here `Q := (𝒲₂·e₂)(𝒲₄·e₄) − (𝒱₂·e₂)(𝒱₄·e₄)` (∝ CLS-D's photon `𝒦^{(1)}_{j′=0}` at `D = 4`), `O := (𝒲₂·e₂)(𝒱₄·e₄) − (𝒱₂·e₂)(𝒲₄·e₄)`, and
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
