# Short analytic lemmas behind Lemma A, H_K and the light side: C′, K1′, K2ν, LRC, and the collinear set

*23 Nov 2026, lemmas agent. Unreviewed, not committed. Answers REFEREE_REPORT §M7 (problems 1–3) and the "what remains" items 1–2 of
HK_CLOSE §9. No existing file was edited.*
* *New code:* `ads4/code/lemmas/lem1_cprime.py`, `lem2_k1prime.py`, `lem3_k2nu_lrc.py`, `lem4_collinear.py`.
* *Outputs:* `ads4/code/lemmas/cache/lem{1,2,3,4}_*.txt` (git-ignored).
* *Read-only imports:* `resultq/ein_dbar.py` (EIN's D-function engine and physical-route continuation), `resultq/v4c_wu.py` (`Kfun`),
  `resultq/stage4_hj.py` (exact `R_j`).
* *Literature read from source (scratchpad copies):* KKSDZ 1904.05905 §4.2–4.4; Caron-Huot 1703.00278 §2.2–2.3.
* *Compute:* about 5 CPU-minutes, at most 2 processes, one thread each, < 1 GB.

Grades: **[P]** proved; **[D]** derived; **[D-sketch]**; **[DV]** derived and checked numerically; **[NE]**; **[M]**; **[Q]** quoted;
**[A]** assumption; **[I]** open. A chain carries the grade of its weakest link.

---

> **Update (R4ARC_HEFT.md, unreviewed): L6′ [I] → [D] with a free-standing proof; H_EFT(a), spinning [D-sketch] → [D]; H_EFT(b) stays [A], with no assumption beyond CHMRS. Result Q's upper end goes from [I] in one link to [D-sketch].**
> * **L6′.** In the og/CLS-D gauge the collinear variable y = 2sin²(Δφ/2) enters only through X₁₃ and an explicit y⁻¹ in f₂, f₄. Counting Lemma Y [P] (exhaustive over 997,836 parity-even and one-ε structures, worst margin 0): every term of a spin-2 contact diagram is O(y⁰ log y). So |L_{R⁴}| ≤ C Δφ⁻²(1+|log Δφ|): p = 2 against the criterion p < 7/2, margin 3/2. It uses neither H_EFT(a) nor LRC\*, and the same bound holds for all contacts and heavy-exchange double-trace parts.
> * **Premise: dispersiveness** (finitely many non-zero s-channel residues δ₁₂ = −n). [D] for integer Mellin shifts, and CLS-D's statement [Q]; every holographic ⟨TTTT⟩ tree term qualifies. In non-dispersive toy classes (half-integer free-field links, 1243) the collinear growth reaches the edge (slope 4.7 → 3.5), so LEMMAS's worry was real for non-holographic terms. LEMMAS §4.1 item 4 and LRC §6.5 item 3 (η^{5/2}, p = 5/2) are superseded.
> * **H_EFT(a) [D] for dispersive terms.** It holds pointwise in angle, hence also ν-paired: only finitely many double-trace residues at bounded Mellin arguments enter. "False pointwise in angle" (F_NORM B3, §3.4, §4) applies only to non-dispersive toy classes. Checks on T′ = σ⊗T: σ-exchanges with J_h = 0, 2, 4 and τ_h = 41–81 at k = 2, 4; the term-wise expansion matches the double-trace value to 10⁻¹⁸–10⁻³⁷. An unexplained ≈ −3 % leading-order ratio for J_h = 2 stays at [NE]; route II needs only the scalar case.
> * **Residual:** neither the actual AdS₄ R⁴ ⟨TTTT⟩ nor a spinning χ-exchange was evaluated. The physical Δφ scan is not yet asymptotic at Δφ = 0.04 (slopes 2.5–2.8), so p = 2 rests on the analytic checks. Refs to add to the bib: Paulos 1107.1504, Costa–Gonçalves–Penedones 1404.5625, Penedones 1011.1485.

## 0. Summary

| item | before | after | key statement |
|---|---|---|---|
| **Lemma C** ("one analytic h on (ℂ∖cut)²") | [D], mis-stated (M7) | **withdrawn**, replaced by C′ | — |
| **Lemma C′** | [DV on examples] (EIN §2.5) | **[D]** given s-channel OPE convergence [Q]; closed form **[DV]** | The swapped-form integrand is one single-valued analytic function on Ω = {Re t₁, Re t₃ > 0; t₁, t₃ ∉ (0, ½]} ⊃ (region between hairpin and C). Reason: on Re t ≥ 0 neither z nor z̄ ever reaches [1, ∞) (an exact identity), so the s-channel OPE converges there. |
| **Lemma A** | [D] on the old C | **[D given C′]**, i.e. [D]. The Einstein instance is **[D], not numerically evaluated** | Lemma A's proof and LIGHTSIDE_DIST's Cauchy argument use only C′. |
| **K1′** | [D-sketch] ("write it out") | **[D]**; scalar examples **[DV]** | `|dDisc 𝒢| ≤ 2B(c)` on the whole Regge box, with an explicit `B(c)` and **no η-dependence**: edge exponent 0. Uniformity holds at the Rindler value J₀ = 1; a bound with the true intercept J₀ < 1 is not uniform. |
| **K2ν** | [D-sketch] | **[D]**, celestial factor **[DV]** | `|K_ν(x)| ≤ K̄(x)` for all real ν. γ = 0 in position normalisation; γ = π in heavy normalisation (from 1/N_com). The heavy-normalised `R_j` are polynomial. |
| **LRC, leading form** | [D] (H7 §14 / §18) | **[D + DV]** for scalar k = 2, 4 and EE k = 2; **[D-sketch]** for EE k = 4 | Heavy image = diagonal (per-harmonic) image of the position kernel `r^a ηⁿP_l` (exact to 10⁻³¹ / 10⁻¹⁶). HK_CLOSE's LRC formula is the *heavy* image; the position kernel's coefficients are `R_j N_com/E_b(ν−ij)`. |
| **LRC, uniform remainder** (Prop K_c step 2) | asserted: "extra factor r" | **corrected; open [I]** | In the large-η part of the Regge box the corrections are O(ρ), ρ = re^β, not O(r) [P + DV]. With K1′ alone the remainder counts like the raw kernel (J₀ < ½). Closing it needs LRC* or LC* (§3.4). **H_K stays [D-sketch]**; its freedom from J₀ < ½ is conditional on LRC* or LC*. |
| **a = b integrability** | [I] (LIGHTSIDE_DIST; EIN §3.6) | partly closed | Kernel bound `|C_ν| ≤ A₁₃^{5/4}K̄`, `K̄ ≈ 2 log(1/Δφ) + 8.69` [D + DV]. Criterion: growth `Δφ^{−p}` with p < 7/2 [D]. The J₀ < 4 part is identically 0 off a = b [D]. The joint-corner collinear limit **is** the large-η edge of the Regge box, `η ≥ 1/sin(Δφ/2)` [P]. Borderline (arc) exponent **[I]**: it is the same corner as LRC. |

**The one thing that turned out hard.** The uniform corner remainder of the de-shifted kernels is a real gap (§3.4), not a write-up.
* The leading-Regge expansion of the per-ν kernel is accurate at fixed η, not uniformly over the Regge box. The ν-dependent bases deviate
  by a relative O(ρ) there, with ρ = re^β ≤ c.
* The de-shift cancellation is proven only for the leading term. HK_CLOSE's statement that the remainder "carries an extra factor r" is
  true at fixed η, false towards the light-cone edges.
* With K1′ alone, Proposition K_c therefore does not yet remove the J₀ < ½ requirement.
* The same corner controls the borderline a = b question (§4).

---

## 1. Lemma C′ (answers M7 problem 1)

### 1.1 Statement

**Setting.** The og_free / EIN gauge:
* points: `x₀ = 0`, `x₂ = −e`, `x₄ = e`, `x₃ = t₃z(b)`, `x₁⁺ = −t₁z(a)`, with celestial angles a ≠ b and `y = 2sin²(Δφ/2) ∈ (0, 2]`;
* effective invariants: `X₁₂ = 2t₁−1`, `X₃₄ = 2t₃−1`, `X₁₃ = 2yt₁t₃`, `X₁₄ = 2t₁+1`, `X₂₃ = 2t₃+1`, `X₂₄ = 4`.

Write `s_i = 1/t_i`. Then

    U = zz̄ = (2−s₁)(2−s₃)/(8y),   V = (1−z)(1−z̄) = (2+s₁)(2+s₃)/(8y).

Let **Ω = {(t₁, t₃) : Re t_i > 0, t_i ∉ (0, ½]}**. In the s-variables this is `{Re s_i > 0, s_i ∉ [2, ∞)}`. Ω contains:
* the region between the hairpin and C, `D′ × D′` with `D′ = {|t − ½| < ½} ∖ (0, ½]`;
* every hairpin-family ellipse `E_κ`;
* `C ∖ {0}` in each variable.

> **Lemma C′.** Let `G = ⟨T₄T₃T₁T₂⟩` be a four-point function of a unitary CFT₃, or any tree-level (Witten-diagram) term of one,
> with the light-ray polarisations fixed. Then `h = t₁ʷt₃ʷ f_k G` (w the light-transform weight) is **one single-valued analytic
> function on Ω**. Its boundary values from the two sides of the slits, `t_i ∈ (0, ½) ± i0`, are the four Wightman orderings of
> `⟨[T₄,T₃][T₁,T₂]⟩`. On Ω the s-channel sheet is fixed by
>
>     log U = L(t₁, t₃) := Log(2 − s₁) + Log(2 − s₃) − log 8y          (principal logs)
>
> Equivalently, every `X_ij` takes its principal branch, with `X₁₃^{−p} := (2y)^{−p}t₁^{−p}t₃^{−p}`. Everything else (`log(1−z)`,
> `Li₂`, blocks in ρ) is principal.

**What C′ does not claim.**
* It claims nothing about `Re t < 0`. There the "wall" `z ∈ [1, ∞)` is reached (§1.3 (A)), so the OPE argument gives nothing.
* Whether h is single-valued on the larger `(ℂ∖(−∞, ½])²` is not decided. No test found a failure there (§1.3 (D)). It is not needed.

### 1.2 Derivation

**Lemma W (no wall) [P].** For `Re s₁ ≥ 0`, `Re s₃ ≥ 0` (t ≠ 0) and every `y ∈ (0, 2]`, neither root of `λ² − (1+U−V)λ + U` lies
in `[1, ∞)`.

*Proof.*
* Multiplying the quadratic by 8y gives the exact identity
  `8y[λ² − (1+U−V)λ + U] = (s₁ + 4λ − 2)(s₃ + 4λ − 2) − 8(2−y)λ(λ−1)` (sympy, `lem1` (A)).
* For λ ≥ 1 and `Re s_i ≥ 0`: `|s_i + 4λ − 2| ≥ 4λ − 2 > 0`. So the modulus of the product is at least `(4λ−2)² = 16λ(λ−1) + 4`.
* That is strictly larger than `8(2−y)λ(λ−1)`, since `0 ≤ 2 − y ≤ 2`.
* So no λ ≥ 1 is a root. The bound is uniform in y and holds down to the collinear limit y → 0. ∎

**Lemma C′ from Lemma W.**
1. **Simple connectivity and no zeros.** On Ω, `U ≠ 0` (U vanishes only at s = 2, i.e. t = ½), `V ≠ 0` (s = −2), and both are finite.
   Ω is a product of two simply connected sets (a half-plane minus a slit from an interior point to ∞).
2. **The ρ-map.** By Lemma W, z and z̄ lie in `ℂ∖[1, ∞)` and are nonzero. The map `ρ(w) = w/(1+√(1−w))²` (principal root) sends
   `ℂ∖[1,∞)` onto the open unit disc. So `|ρ(z)|, |ρ(z̄)| < 1`, and on any compact `K ⊂ Ω` they are bounded by some `ρ_K < 1`.
3. **OPE convergence [Q].** The s-channel (12)(34) expansion converges absolutely and uniformly on compacts of `|ρ|, |ρ̄| < 1`. This holds
   with ρ and ρ̄ independent complex numbers, for any external spins and complex polarisations:
   * KKSDZ §4.2 (first-sheet convergence, citing Mack and Pappadopulo et al.);
   * Caron-Huot §2.2 ("this holds whether or not ρ and ρ̄ are complex conjugate");
   * the radial-quantisation Cauchy–Schwarz bound on `⟨Ψ₃₄|r^D e^{iθM}|Ψ₁₂⟩`.

   Its terms are `c_{h,h̄}ρ^hρ̄^{h̄}`, with h − h̄ ∈ ℤ in a fixed frame. The only multivaluedness is through
   `(ρρ̄)^{min(h,h̄)}`, i.e. through `log ρρ̄ = log U − 2Log[(1+√(1−z))(1+√(1−z̄))]`. The second log is principal, because each factor
   has real part > 1.
4. **The continued log.** `L` is analytic on Ω, since `Log(2−s)` has its cut exactly at `s ∈ [2,∞)`. It satisfies `e^L = U`, and it is
   real and equal to the Euclidean `log U` at `t₁, t₃ → +∞`.
   * So the OPE series with `log U := L` is a series of functions analytic on Ω (symmetric in z ↔ z̄, hence analytic in (U, V) across
     z = z̄), converging uniformly on compacts. By Weierstrass, h is analytic and single-valued on Ω.
   * The prefactors are principal powers of `2t−1`, t and `2t+1`, all analytic on Ω. `f_k` is rational, with poles only at `X_ij = 0`,
     i.e. outside Ω.
   * `L = Log X₁₂ + Log X₃₄ − Log X₁₃ − Log X₂₄` with principal per-factor logs on Ω. The difference can only jump where `(2t−1)/t` crosses
     ℝ₋, i.e. on the slit.
5. **Identification.**
   * On `t₁, t₃ ∈ (½, ∞)` all pairs are spacelike and h is the common value of the four orderings (EIN's physical route, step 1, lies in Ω).
   * Approaching a slit point `t₁ ∈ (0, ½)` from `Im t₁ ≷ 0` gives `X₁₂ → −|X₁₂| ± i0`. This is the Wightman iε of "1 left/right of 2"
     (og_free §2.2 convention [D]; Lemma C of Stage 3 for the same iε on f's poles [P]).
   * The four orderings of the double commutator are therefore the four one-sided boundary values of the single function h. The real
     double-commutator integral equals the hairpin integral, and Cauchy in each variable inside Ω moves it to C. Endpoints and edges are
     unchanged: `t = 0` is endpoint convergence (k > J₀), and `t = ½` is Gel'fand–Shilov (Stage 3 §1.4) [D].
   ∎

**Tree-level terms (the Einstein instance).** The O(1/c_T) Einstein correlator is not itself a unitary correlator. Two routes give C′ for it:
* **(i) Witten diagrams [Q].** A finite sum of tree Witten diagrams is analytic on the s-channel sheet away from `z, z̄ ∈ {0, 1, ∞}`. The
  bulk-point singularity z = z̄ sits only on the sheet reached by crossing [1, ∞) (Maldacena–Simmons-Duffin–Zhiboedov 1509.03612;
  Gary–Giddings–Penedones). Lemma W says Ω never reaches that sheet.
* **(ii) Large-N expansion [D].** Take the O(1/c_T) coefficient of the unitary correlator's convergent OPE (anomalous-dimension logs are
  powers of L).

Grade: **[D given Q]**.

### 1.3 Checks (`lem1_cprime.py` → `cache/lem1_cprime.txt`)

**(A) Wall algebra and sampling.**
* The identity holds exactly (sympy). The strict separation `(4λ−2)² − 16λ(λ−1) = 4` holds.
* Sampling: 4·10⁵ random points with `Re s ∈ [0, 3]`, plus 4·10⁵ hugging `Re s ∈ [0, 0.02]`, with all y. The minimum distance of z, z̄
  from [1, ∞) is **0.500**, and no point is on the wall.
* Control `Re s < 0`: 4.7 % of samples lie within 10⁻² of the wall. An explicit wall point: `s₁ = s₃ = −0.3253`, `y = 0.05`, gives
  `z̄ = 3` exactly.

**(B) The C′ closed form against EIN's tracked physical route** (`ein_dbar.continue_ellipse`). For D̄₁₁₁₁ the closed form is
`Φ_C′ = [2Li₂(z) − 2Li₂(z̄) + L·(log(1−z) − log(1−z̄))]/(z−z̄)`.

| Δφ | κ = 1 (C) | κ = 0.5 | κ = 0.2 (near the hairpin) |
|---|---|---|---|
| 0.1 | 2.6·10⁻¹⁰ | 7.2·10⁻¹⁰ | 1.8·10⁻⁹ |
| 1.0 | 8.9·10⁻¹² | 2.2·10⁻¹¹ | 2.9·10⁻¹¹ |
| 2.5 | 6.0·10⁻¹² | 2.0·10⁻¹² | 2.9·10⁻¹¹ |

* Entries are the maximum relative difference over 576 node pairs, both continuation orders.
* Every EIN piece is `ΠX^e·D_U^mD_V^nΦ` with integer powers. So this checks the C′ sheet for all of EIN's Witten family **[DV]**, and it
  explains EIN's order and κ independence.

**(C) Path tests** (`ein_dbar.continue_path`).
* Routes inside Ω reproduce Φ_C′ to **3·10⁻¹⁶**:
  * via `Re s = 0.03` (t near the imaginary axis);
  * through `s ∈ (1, 2)` twice (t through (½, 1));
  * a large loop around t = ½ inside Re s > 0.
* Crossing the slit at `s = 3` (t = ⅓) changes the value by exactly the s-channel monodromy `L → L − 2πi`, to **4.7·10⁻¹⁶**.

**(D) Outside Ω.**
* Two paths around discriminant points `a*` (z = z̄) with `Re a* < 0`, which pass within 3·10⁻³ of the wall, give the same value to
  10⁻¹⁶.
* So no monodromy was produced in `(ℂ∖(−∞, ½])²`. Single-valuedness there is not decided and not needed.

**(E) EIN's "joint interior" route.** The route is `(t₁, t₃) → (0.8∓0.4i)`, then straight lines in t to near-endpoint nodes.
* When exactly one leg crosses the real t-axis inside (0, ½), the value is `Φ_C′ + 1 × (2πi monodromy term)`, to 1.5·10⁻¹⁷.
* When both legs cross with opposite orientation, the two monodromies cancel.
* So the disagreement EIN saw "at endpoint nodes" is what a straight-in-t leg to the opposite arc produces: it crosses the
  double-commutator segment. It is not a failure of analyticity on the region between the hairpin and C.

### 1.4 What Lemma A and LIGHTSIDE_DIST actually use

| step | input | needs C′? |
|---|---|---|
| V4FINAL §2.2 Step 1: commutator = 0 at fixed a ≠ b | microcausality at real points | no |
| Step 2: t-integrals converge | k > J₀ (endpoint), K2 edges, Gel'fand–Shilov at t = ½ | no |
| Step 3a: vacuum annihilation, `∫dt₁⟨…T₁(t₁)|0⟩ = 0`, f inside | single-variable half-plane analyticity (spectral condition) + Stage 3 Lemma C ("f carries the Wightman iε") [P] | no (the old Lemma C's correct part) |
| Step 3b: real double-commutator integral = hairpin = C × C | one analytic function on the swept region | **yes: C′** |
| Step 4 + EIN/LSD pair identities node by node on C × C | value at a C × C node independent of route | **yes: C′** (`C∖{0} ⊂ Ω`) |
| LSD §1.2: MB representation with principal per-`X_ij` branches = physical continuation | the sheet of h on C × C | **yes: C′** (`L` = sum of principal per-factor logs, §1.2 item 4) |
| LSD §1.2: K entire in δ₁₂ | `|U/V| ≤ 1` | `|2−s| ≤ |2+s|` iff Re s ≥ 0, so it holds on all of Ω̄ |
| LSD: term-wise residues are the same on the hairpin as on C | each OPE term analytic on Ω | per-term C′ (trivial) |
| LSD: Cauchy in δ₁₂ at a fixed node | Mellin-space only | no (needs Mellin convergence on C × C, [DV] for the family) |

So **Lemma A is [D given C′]** and LIGHTSIDE_DIST's decomposition `Ω_k = D_LR + 𝒜_k` is **[D given C′ + Mellin convergence]**. Nothing
used the false "analytic on (ℂ∖cut)²".

**Einstein instance (M7 problem 2).** Every input of Lemma A holds for the Einstein tree correlator:
* microcausality [D];
* J₀ = 2 < 4 [Q, Costa–Hansen–Penedones];
* C′, via route (i) [D];
* relabelling [P].

So the instance is **[D]**. What is missing is a numerical evaluation, not a derivation.

### 1.5 Grade and wording

* Lemma C (old): **withdrawn**.
* **Lemma C′: [D]** (given OPE convergence [Q] and the og_free iε convention [D]). Its closed-form consequence is **[DV]** on EIN's Witten
  family.
* **Lemma A: [D given C′] = [D].** The Einstein instance is **[D], numerically unevaluated**.

> *Letter:* "The swapped-form integrand is a single analytic function on the region swept between the hairpin and the contour
> (Lemma C′): for Re t₁, Re t₃ ≥ 0 neither s-channel cross-ratio can reach [1, ∞), so the (12)(34) OPE converges there, and its four
> one-sided boundary values on the double-commutator segment are the four Wightman orderings. Lemma A uses nothing more; it therefore
> holds, as a derivation, for the Einstein term itself, which we have not evaluated numerically."

---

## 2. K1′: the Rindler Regge bound, uniform over the Regge box

### 2.1 Where K1′ is used, and the precise statement

**Where it is used.** H7 §14 (Corollary) and HK_CLOSE §1 / Prop K_c use it in the form
`|dDisc 𝒢| ≤ B r^{1−J₀}` uniformly in η on the box `{ρ = re^β ≤ c}`, with J₀ ≤ 1.

**Frame.** Caron-Huot 1703.00278 eq. (2.5): `x₄ = (ρ, ρ̄) = −x₃`, `x₁ = (1, 1) = −x₂` in light-cone coordinates. The Regge sheet is
`0 < ρ < 1 < ρ̄`. Put `ρ_b := 1/ρ̄`. The project's box is `B_c = {0 < ρ ≤ c, 0 < ρ_b ≤ c}`, so `r = √(ρρ_b)`, `η = (ρ+ρ_b)/2r`. The
large-η corner is the light-cone edge `ρ_b/ρ → 0` (or `ρ/ρ_b → 0`).

> **K1′ (statement).** Let `𝒢` be ⟨TTTT⟩ stripped and normalised by two-point functions with generic polarisations (KKSDZ eq.
> `J0definition`), on the Regge sheet. For `0 < c < 1` and all normalised polarisations:
>
>     |𝒢(ρ, ρ̄ ∓ i0)| ≤ M(ρ, ρ_b) ≤ B(c) := sup_{B̄_c, pol} M,      |dDisc 𝒢| ≤ 2B(c),
>
> where `M = ½[𝒢^E_A + 𝒢^E_B]`, built from KKSDZ's first-sheet, Rindler-positive comparison correlators at the same cross-ratios.
> * `B(c) < ∞`, and it is **independent of η**: in the project's variables `|dDisc 𝒢(r, η)| ≤ 2B(c)` for all `r ≤ ce^{−β}` and all β.
> * Identical scalars: `B(c) = 𝒢_E(c, c)`.
> * Edges: `M → M_identity` as `ρ_b → 0` or `ρ → 0`, with `|M − M_id| ≤ M(c,c)(min(ρ, ρ_b)/c)^{τ_min/2}`. τ_min is the lowest non-zero
>   twist in the comparison OPEs.

This is H7 §14's K1′ with J₀ = 1 and `B = 2B(c)`. Every H_K count is done at J₀ ≤ 1, so this is all that is needed.

### 2.2 Derivation
1. **Cauchy–Schwarz [Q].** By Rindler positivity (Casini) and Cauchy–Schwarz,
   `|⟨O₄O₁O₂O₃⟩|² ≤ ⟨Ō₂Ō₃O₂O₃⟩⟨O₄O₁Ō₄Ō₁⟩` (KKSDZ eq. `finalCSbound`). The two factors are first-sheet correlators at the **same**
   cross-ratios. AM–GM turns the geometric mean into M.
   * For identical scalars this is Caron-Huot's `|𝒢(ρ,ρ̄)| ≤ 𝒢(ρ, 1/ρ̄)` (eq. `inequality`), from t-channel positivity.
   * dDisc is a combination of 𝒢_E and the two continuations with phase weights of modulus ≤ 1 (Caron-Huot eq. `dDisc_full`). Hence
     `|dDisc| ≤ 2M`.
2. **Uniform bound on the closed box [D].**
   * The comparison correlators are first-sheet. Their s-channel OPE converges uniformly for `|ρ|, |ρ_b| ≤ c < 1` [Q, as §1.2 step 3].
     So they are continuous on the compact closed box `B̄_c`, including its edges `ρ_b = 0` or `ρ = 0`, where only the identity survives.
   * They are polynomial in the polarisations. KKSDZ show that the configuration-dependent polarisations and scale factors stay finite
     at r → 0, and become aligned (`z̃₁ → cz₁`) at the light-cone edges.
   * Taking the supremum over a compact set of normalised polarisations gives `B(c) < ∞`.
   * For identical scalars the s-channel coefficients are positive, so `𝒢_E` is monotone in each variable and the supremum sits at the
     corner (c, c) [P given Q].
3. **Explicit edge dependence (spinning) [D].** Work in KKSDZ variables `rσ = 1/α₁`, `r(1−σ) = −1/α₂`, where σ → 0, 1 are the edges.
   * The identity piece of `⟨Ō₂O₂Ō₃O₃⟩` with aligned polarisations scales as `σ^{J−Δ}(1−σ)^{−Δ−J}r^{−2Δ}`.
   * The partner scales as `(1−σ)^{J−Δ}σ^{−Δ−J}r^{−2Δ}`.
   * Their geometric mean is `(σ(1−σ))^{−Δ}r^{−2Δ}`. That is exactly the generic-polarisation normalisation `(x₁₂²)^{−Δ}`.
   * So the normalised bound has **exponent 0 at both edges**: no η-dependence. The non-identity terms are suppressed by
     `(edge distance)^{τ/2}` (KKSDZ eq. `f0general`).
4. **The loophole is outside the box [P].** KKSDZ's loophole sits at the first/second-sheet boundary, `x₁₄² = 0`, i.e. ρ̄ = 1. On `B_c`
   one has `ρ̄ = 1/ρ_b ≥ 1/c > 1`. The loophole concerns the light-cone rows (K2 edges), not K1′.
5. **Why J₀ = 1 [DV, §2.3].** A bound `B r^{1−J₀}` with the true intercept J₀ < 1 is in general **not** uniform towards the edges. Only
   the Rindler form (r⁰) is.

### 2.3 Checks (`lem2_k1prime.py` → `cache/lem2_k1prime.txt`; exact unitary correlators, scalar externals)

**2d Ising ⟨σσσσ⟩.** Exact form: `𝒢 = (1+√(ρρ̄))/((1−ρ²)(1−ρ̄²))^{1/4}`. On the Regge sheet,
`|𝒢_R| = (√ρ + √ρ_b)/((1−ρ²)(1−ρ_b²))^{1/4}`. The inequality `|𝒢_R| ≤ 𝒢_E` is equivalent to `(1−√ρ)(1−√ρ_b) ≥ 0`.

| quantity | c = 0.5 | c = 0.9 |
|---|---|---|
| sup over the box of `|𝒢_R|/𝒢_E(c,c)` | 0.943 | 0.9986 |
| sup over the box of `|dDisc|/2𝒢_E(c,c)` | 0.289 | 0.147 |
| dDisc ≥ 0 everywhere | yes | yes |
| edge `ρ_b = 10⁻⁹`: max over ρ of `|𝒢_R|` | 0.760 | 1.437 |
| sup of `|𝒢_R|/r^{1−J₀}` at the true J₀ = ½ | ≥ 1.6·10² | ≥ 2.6·10² |

* Pointwise `|𝒢_R| ≤ 𝒢_E(ρ, ρ_b)` everywhere on the box.
* At the edge the value stays finite, `→ √ρ/(1−ρ²)^{1/4}`, so it is uniform.
* The true-J₀ ratio is attained at `ρ = 10⁻⁹`, the grid edge, and grows without bound as the grid extends. **The true-intercept form is
  not uniform.**

**GFF, Δ = 0.5, 1.3, 3.** `sup |𝒢_R|/𝒢_E(c,c) ≤ 1.0000`. dDisc is ≥ 0 up to float cancellation. The edge and Regge-corner limits are
→ 1 (identity).

**Not checked numerically.** The spinning version. Its derivation is the same Cauchy–Schwarz with a supremum over polarisations.

### 2.4 Grade and wording
**K1′: [D-sketch] → [D]** (inputs [Q]: Rindler positivity / Cauchy–Schwarz, first-sheet OPE convergence, KKSDZ's polarisation limits).
The scalar examples are [DV].

> *Letter:* "Rindler positivity bounds the Regge-sheet correlator by first-sheet correlators at the same cross-ratios, and these are
> continuous on the closed Regge box; hence |dDisc 𝒢| ≤ 2B(c) uniformly in the transverse angle, including towards the light-cone
> edges (KKSDZ's loophole lies outside the box). The uniform statement holds at the Rindler value J₀ = 1, which is all the interchange
> argument uses."

---

## 3. K2ν and LRC for Ĉ′_EE, Ĉ₄^{ds}, ∂Ĉ₄^{ds}

### 3.1 K2ν: statement and derivation

> **K2ν.** For the j′ = 0 per-ν functionals of `Ĉ′_EE` (k = 2) and `Ĉ₄^{ds}`, `∂Ĉ₄^{ds}` (k = 4), the per-ν position kernel satisfies,
> for every real ν and every point x of the CLS-D domain,
>
>     |𝒦_ν(x)| ≤ K̄(x)    (γ = 0 in position normalisation),
>
> where `K̄` is the same kernel with every complex power replaced by its modulus (`δ → ½`). Hence:
> * heavy normalisation: `|𝒦_ν/N_com| ≤ C ν⁻⁹e^{πν} K̄` (γ = π);
> * leading-Regge heavy coefficients: `|R_j(ν)| ~ ν²` (k = 2) and `ν⁶` (k = 4), polynomial.
>
> `K̄` is integrable against `|dDisc 𝒢|` away from the Regge corner, by K2's endpoint and edge rows. Those depend only on real parts of
> exponents, which are ν-independent.

**Derivation [D].**
1. **The fibre integral.** The per-ν kernel is an integral over the light-ray parameters, the circle `x₀ = |y|n̂_θ` and the celestial
   point z₀. Its integrand is a ν-independent factor F times complex powers `B_i^{−e_i(δ)}`, with `δ = ½ + iν` and `Re e_i ∈ {½, 1}`.
   * F contains: `f_k`, the measure `t³t³`, the polarisation polynomials `τ₁₃`, `q²/V⁴` and the EE structure (H7 §18.1).
   * The B_i are real bases on the real domain:
     * celestial: `A₁₀, A₃₀ = 2−2cos(·) ≥ 0`, giving `C_ν = A₁₃^{(2+δ)/2}∫dφ̄(A₁₀A₃₀)^{−δ/2}`;
     * dual: `−z₀·x₂₀`, `X₂₀`, `X₄₀`.
2. **Fixed signs.** The bases keep a fixed sign. The celestial ones are trivially ≥ 0. For the dual ones this was checked on the whole
   Regge box in the exact frame (§3.2 (b), all θ, φ, β up to 16), and outside the box CLS-D's domain keeps 2, 4 spacelike from x₀.
3. **Modulus.** Hence `|B^{−e(δ)}| = |B|^{−Re e}` for every real ν, and `|𝒦_ν| ≤ K̄` follows by the triangle inequality under the fibre
   integral.
4. **Worst case.** If a base did change sign somewhere, `|B^{−iν}| ≤ e^{π|ν|}`. With two ν-dependent powers that gives `γ ≤ 2π`, still
   finite. Prop K_c needs only a finite γ, because `c(ν)` decays super-Gaussianly.
5. **Heavy side.** `N_com = 2¹¹π²|Γ(11/4+iν/2)|⁴ ~ e^{−πν}ν⁹`. The `R_j` are the exact rational functions of stage4_hj.

**Checks** (`lem3` (a), (d)).

| Δφ | 10⁻³ | 10⁻² | 0.1 | 0.5 | 2 | 3 |
|---|---|---|---|---|---|---|
| `K̄(Δφ)` | 22.50 | 17.90 | 13.29 | 10.11 | 8.14 | 7.42 |
| `max_ν |K_ν|/K̄` (ν ∈ {0.5,…,80}) | 0.20 | 0.24 | 0.25 | 0.44 | 0.66 | 0.74 |

* `K̄ = 2 log(1/Δφ) + 8.69 + o(1)`, checked to 3·10⁻³ at Δφ = 10⁻⁴.
* Near the corner `K_ν = a + bΔφ^{−2iν}` (V4DPRIME), with **`|a| = |b| = √(π/ν)`** to 4 digits at ν = 2, 10 and 40. So `|C_ν|` decays like
  `ν^{−1/2}`, as V4FINAL quoted.
* The quadrature agrees with `v4c_wu.Kfun` to 2·10⁻⁴, and with the corner expansion to 10⁻⁴.
* Heavy degrees, from the exact `R_j`: `p = 2` for sc1 and EE1; `p = 6` for sc3 and EE3.

**Grade: K2ν [D-sketch] → [D]**; the celestial factor is [DV].

### 3.2 LRC, leading form

> **LRC-a.**
> * **Leading term.** In the CHMRS Regge frame, as r → 0 at fixed η, the per-ν position kernel is `r^a c_k(ν) ηⁿ 𝒦₀,ν(η)(1 + O(r))`:
>   * `(a, n) = (0, 1)` for k = 2, because `f₂ = η(1+r²)/2r` exactly;
>   * `(a, n) = (2, 3)` for k = 4, because `f₄ = f₂R` and `R → 16r²η²`.
>
>   `𝒦₀,ν` is the f-stripped kernel. In the scalar structure it is `P_{−½+iν}(η)`, by the split representation [P].
> * **Heavy map.** The plane-wave (heavy) map sends the harmonic `P_{−½+iν′}(η)` to `E_b(ν′)P_{−½+iν′}(cosh β)`, with `b = 7` (k = 2) and
>   `b = 9` (k = 4). Hence `ηⁿP_l = Σ_j a_j(l)P_{l+j}` (Legendre recurrence) has heavy coefficients
>   `R_j ∝ a_j(l)E_b(ν−ij)s(ν−ij)/N_com`, with s = 1 for scalars and `s(w) = 32(4w²+1)(4w²+9)/((4w²+25)(4w²+81))` for EE at k = 2.

**Check against the stored exact R_j** (`lem3` (d)). Spread of `R_j N_com/(a_j E_b(ν−ij) s(ν−ij))` over j and ν ∈ {0.37, 2.9, 11.3, 41}:

| key | result |
|---|---|
| sc1 (scalar, k = 2) | 1.8·10⁻³¹ |
| EE1 (EE, k = 2) | 2.4·10⁻¹⁶ |
| sc3 (scalar, k = 4) | 3.7·10⁻³¹ |
| EE3 (EE, k = 4) | **not** per-harmonic: O(1) deviation. The ratios `R_j^{EE3}/R_j^{sc3}` are j-dependent rational functions (printed in the cache), with poles at `ν = ±3i/2` inside the de-shift strips |

So for EE at k = 4 the rank-3 contraction (stage4_hj) mixes harmonics. The position kernel is not `η³ ×` (a per-harmonic image) there.

**Correction to HK_CLOSE §4.1 (LRC).**
* The formula "per-ν kernel = `r^a N_com Σ_j R_j P_{l+j}(η)`" is the **heavy image** read as a function of η. The position kernel has
  coefficients `R_j N_com/(E_b(ν−ij)s(ν−ij))`.
* For the scalar structures and EE at k = 2 the corner argument nevertheless goes through. `E_b` has no zeros, and its poles lie at
  `|Im ν′| ≥ b − ½ ≥ 6.5`, outside the strips `|Im| ≤ j_max`. So the contour shifts of HK Lemma S2/S4 apply verbatim to
  `w_c a_j = χ_cR_j/(E_b s)(ν−ij)`.
* After the shift, the position spectral function is the heavy one divided by `E_b(ν′)s(ν′)` on the real line. The leading-Regge
  position profile of the cutoff functional is therefore `≲ P_{−½}(η)` whenever the heavy image is [D].
* For EE at k = 4 the same needs the position-space rank-3 kernel, which is not computed [D-sketch].

**Grade: LRC-a [D] → [D + DV]** for scalar k = 2, 4 and EE k = 2; **[D-sketch]** for EE k = 4.

### 3.3 LRC, the uniform remainder: what is true

**Geometry of the corrections [P + DV]** (`lem3` (b)).

Exact frame: `ŷ = (1,0,0)`, `x̂ = (cosh β, sinh β, 0)`, `x₂ = x/x² = −x₄`, `x₀ = |y|n̂_θ`, `n̂_θ = (0, cos θ, sin θ)`.

* **The dual invariant:** `X₂₀|x|² = −(1 − 2r sinh β cos θ − r²)`. Its maximal relative deviation from the leading-Regge value is
  `2r sinh β ≈ ρ`.
* **The z₀-bases:** `|y z₀·n̂|/(−z₀·x̂/|x|) ≤ re^β = ρ`.
* **Measured** (721 × 721 grid over θ and z₀). The relative correction equals ρ, independently of r:
  * `0.300` at ρ = 0.3, for β = 5 and 12 (r = 2·10⁻³ and 1.8·10⁻⁶);
  * `0.900` at ρ = 0.9;
  * `0.496` at ρ = 0.5, β = 16.

  All signs stay fixed.
* **The f-factor (Stage 3 Lemma B) in radial variables [P]:**
  `R = [4(ρ+ρ̄)²(1+ρρ̄)² − 64ρ²ρ̄²]/((1−ρ²)(1−ρ̄²))²`. Its relative correction is O(ρ²) (2.03·10⁻² at ρ = 0.1; 0.78 at ρ = 0.5,
  independent of r), and it is **ν-independent**.

**Consequences.**
1. **The f-part is harmless [P].** A ν-independent bounded factor `g(ρ, ρ̄)` commutes with the ν-integral:
   `|g·K_c^{LR}| ≤ ‖g‖_∞|K_c^{LR}|`. The same holds for the measure.
2. **The ν-dependent part is O(ρ), not O(r) [P].** At fixed η, ρ = O(r), so H7 §18.1's "corrections O(r)" is right there. On the box,
   ρ ≤ c is bounded but not small towards the light-cone edges.
3. **Without more input the remainder counts as raw [D].** Use the only bound available pointwise in ν:
   `|remainder_ν| ≤ C(1+|ν|)ρ·K̄^{LR}_{raw}`, with `K̄^{LR}_{raw} ~ r^aη^{n−½}`.
   * In (ρ, ρ̄) the edge exponent of this term is `ρ̄^{(a−n+½)/2} = ρ̄^{−1/4}` (k = 2 and k = 4 alike). A factor ρ is neutral at the edge.
   * With Stage 3's measure (`dr dη ∝ dρdρ̄/ρ̄` at the edge) and K1′, integrability needs `−¼ + (1−J₀)/2 > 0`, i.e. **J₀ < ½**.
   * HK_CLOSE's count "`∫η^{J₀−5/2}`" assumed an extra r. With the correct extra ρ it becomes `∫η^{J₀−3/2}`.
4. **No modulus-only rescue [D].** The de-shift cancellation of the `ηⁿ` prefactor is an interference effect. A modulus bound on the
   ν-Fourier transform of the weight reproduces only the raw growth `η^{n−½}` (§3.4 note). Bounded O(ρ) shifts of the "time" variable
   therefore cannot be absorbed by absolute bounds.

**Verdict.** Prop K_c step 2 ("subleading terms carry an extra factor r") holds at fixed η, not uniformly. **The uniform remainder is open
[I].** HK_CLOSE's margins (J₀ < 3/2 for k = 2, J₀ < 5/2 for k = 4) are established for the leading-Regge part of the corner only. The same
gap is present in H7 §14's Corollary for D's Ĉ′_EE, which also used only the leading-Regge kernel on the whole box.

### 3.4 What would close it

Either of two inputs suffices. Each is graded [I].

**LRC\* (kernel side).** The de-shift cancellation survives the O(ρ) corrections. Precisely: at the edge ρ̄ → 0 the exact per-ν kernel
expands as `Σ_j ρ̄^{(a+½−j−iν)/2}F_j(ρ; ν)(1 + O(ρ̄))`, with `F_j(ρ; ·)` analytic in the strips `0 ≤ ±Im ν ≤ |j|`,
exponentially bounded on their lines uniformly for ρ ≤ c, and with no poles there beyond those of the leading-Regge coefficient.
* Then the contour shifts of HK Lemma S2/S4 apply at every fixed ρ, and the edge exponent of the cutoff kernel is the de-shifted one,
  `ρ̄^{(a+½)/2}`.
* Plausibility: the ν-dependence enters only through positive bases `B_i^{LR}(1+O(ρ))`. A shift `ν → ν + ij` multiplies by
  `(1+O(ρ))^{∓j} ≤ (1−c)^{−3}`.
* What is missing: the decomposition into ν-analytic edge pieces for the exact (non-leading) kernel. Natural route: the conformal Casimir
  equation of the per-ν kernel in (ρ, ρ̄), whose indicial exponents at ρ̄ = 0 are exactly the leading-Regge ones.

**LC\* (correlator side).** A uniform light-cone suppression on the Regge box of the dDisc that enters Lemma S (heavy-projected, Rindler
double commutators): `|dDisc| ≤ B·(min(ρ, ρ̄)/max(ρ, ρ̄))^κ` with **κ > ¼**.
* The remainder then counts as `−¼ + κ + (1−J₀)/2 > 0`, which holds at J₀ = 1 for both k = 2 and k = 4.
* At fixed ρ, KKSDZ's asymptotic light-cone expansion gives `κ = τ_*/2 ≥ ½`. That is the same assumption that closes K1's loophole and
  that K2's edge rows already use for validity.
* LC\* is its uniform version: the Regge and light-cone limits commute.
* Caveat [DV]: for identical scalars the full dDisc does **not** decay at the edge, because of the identity (Ising: dDisc → `(1−√(ρ/2))/(1−ρ²)^{1/4}`).
  So LC\* must be formulated for the heavy-projected, light-ray-polarised object, where the identity is absent.

*Note to point 4 of §3.3.* In split form,
`∫w(ν)P_{−½+iν}(η)dν = ∫Dz₀(B_xB_y)^{−½}ŵ(log B_x/B_y)`.
With `|ŵ(T)| ≤ Ce^{−s|T|}` the z₀-integral is `≍ η^{−½}` for every s, so `ηⁿ ×` it gives `η^{n−½}`: raw.

### 3.5 Grades and wording for §3

| item | grade |
|---|---|
| K2ν | **[D-sketch] → [D]**; celestial factor [DV] (γ = 0 position, π heavy) |
| LRC-a (leading form) | **[D + DV]** scalar k = 2, 4 and EE k = 2; **[D-sketch]** EE k = 4 |
| LRC uniform remainder | **[I]** (was implicitly [D] in HK_CLOSE) |
| Prop K_c | [D given K1′, K2ν, LRC-a, K3, DS] **+ (LRC\* or LC\*)** |
| H_K | stays **[D-sketch]**. [D given H_ip] remains valid, since the raw count needs J₀ < ½ and covers the remainder. HK_CLOSE's "no longer needs J₀ < ½" holds **only given LRC\* or LC\*** |

> *Letter (replaces HK_CLOSE §8's Tonelli bullet):* "For the leading-Regge part of the corner, the de-shift identity forces strip
> analyticity of the weights and the cutoff spoils the cancellation only at O(ε_c), so the heavy sum converges for J₀ < 3/2; the
> per-ν kernels themselves are bounded uniformly in ν (K2ν). The sub-leading corner terms are suppressed by ρ = re^β rather than r,
> and controlling them uniformly towards the light-cone edges needs one further input — either that the de-shift cancellation persists
> beyond leading Regge order, or a uniform light-cone suppression of the double discontinuity. We therefore grade the interchange
> [D-sketch]; under the stronger hypothesis J₀ < ½ it is derived."

---

## 4. The collinear set a = b (LIGHTSIDE_DIST [I]; EIN §3.6; M7 problem 3)

### 4.1 Statements
1. **Kernel bound and criterion [D].**
   * `|C_ν(φ₁, φ₃)| ≤ A₁₃^{5/4}K̄(Δφ)`, with `K̄ = 2 log(1/Δφ) + 8.69 + o(1)`, uniformly in ν (§3.1).
   * Hence a light function with `|L(a,b)| ≤ CΔφ^{−p}` near a = b pairs absolutely with C_ν, uniformly in ν and then against any
     `w_c ∈ L¹`, **iff p < 7/2**. The measure transverse to the diagonal is dΔφ; the integrand is `Δφ^{5/2−p}log(1/Δφ)`.
2. **J₀ < 4 tree terms [D].**
   * Their symmetrised light function is identically zero on a ≠ b (Lemma A, now [D given C′]).
   * D_LR integrates over a ≠ b only, as a Lebesgue integral (H_LT → D_LR, V4FINAL §1).
   * So there is nothing to integrate and no distribution at a = b. "The collinear set carries nothing" is [D] for these terms, not
     only [DV] at Δφ ≥ 0.025.
   * Individual families are never paired separately. If they were:
     * the free-scalar family law Δφ⁻³ (OG §3.6) gives p = 3 < 7/2, integrable with margin ½;
     * EIN's Witten-family local slopes (2.6–3.4) refer to T′ (Δ = 5), whose own celestial kernel differs. They are not covered by this
       criterion.
3. **The joint-corner collinear limit is the Regge-box edge [P].** At the joint endpoint `t₁, t₃ → 0` with `t₃/t₁ = λ` fixed,
   `η² = (z+z̄)²/(4zz̄) → (1+λ)²/(2yλ) ≥ 2/y = 1/sin²(Δφ/2)`.
   * Checked: `|η| = 2.000000·10³` at Δφ = 10⁻³, t = 10⁻⁷ (`lem4`).
   * So Δφ → 0 pushes the whole joint corner to η → ∞: the light-cone edge of the u-channel Regge box of §3.3–3.4.
4. **The borderline (arc) light function [I].**
   * Its a = b behaviour is the large-η behaviour of the borderline term's Regge-corner contribution, `η ≈ 2/Δφ`.
   * Integrability needs that contribution to grow slower than `η^{7/2}`.
   * The raw k = 4 per-ν growth `η^{5/2}` would give p = 5/2, margin 1. But that identification is heuristic [D-sketch]: the arc is a
     Laurent coefficient, not a per-ν kernel.
   * LIGHTSIDE_DIST's three toy values (P1–P3) have different polarisations and cannot fix an exponent.

### 4.2 Grade and wording
* The criterion and the J₀ < 4 part: **[I] → [D]**.
* The geometric identification: **[P]**.
* The borderline arc's integrability: **stays [I]**. It is now located: the same light-cone edge of the Regge corner as LRC\* / LC\*.

> *Letter:* "The light function of every J₀ < 4 term vanishes identically off the collinear set, and D_LR integrates only there, so
> the collinear set carries nothing. For the R⁴ term the collinear limit at the joint endpoint coincides with the light-cone edge of
> the Regge limit (η ≥ 1/sin(Δφ/2)); its integrability against the celestial kernel (growth slower than Δφ^{−7/2}) is not yet
> established."

---

## 5. Proposed hypothesis-table rows (UT §3.1; not edited here)

| # | grade now | proposed |
|---|---|---|
| K1′ | [D-sketch] | **[D]** (LEMMAS_CLOSE §2); uniform at J₀ = 1, explicit B(c); loophole outside the box |
| Lemma A + D_LR | [D] (on old Lemma C) | **[D given C′]**, C′ [D] + closed form [DV] (§1). Einstein instance **[D], numerically unevaluated** |
| H_K | [D-sketch] | **[D-sketch]**. K1′ and K2ν are now [D]. The open link is the **uniform corner remainder (LRC\*/LC\*) [I]**; [D given H_ip] stands |
| K2ν (new row from HK_CLOSE) | [D-sketch] | **[D]**, γ = 0 position / π heavy |
| LRC | [D] | leading form [D + DV] (EE k = 4: [D-sketch]); HK_CLOSE's "extra factor r" corrected to ρ |
| a = b (LIGHTSIDE_DIST) | [I] | J₀ < 4 part [D]; borderline arc [I], same corner as LRC\* |

## 6. Reproduce

    cd ads4/code/lemmas
    ulimit -v 3000000
    OMP_NUM_THREADS=1 python3 lem1_cprime.py      # ~3 min -> cache/lem1_cprime.txt
    OMP_NUM_THREADS=1 python3 lem2_k1prime.py     # s      -> cache/lem2_k1prime.txt
    OMP_NUM_THREADS=1 python3 lem3_k2nu_lrc.py    # ~1 min -> cache/lem3_k2nu_lrc.txt
    OMP_NUM_THREADS=1 python3 lem4_collinear.py   # s      -> cache/lem4_collinear.txt

| file | content |
|---|---|
| `lem1_cprime.py` | wall identity (sympy), wall sampling, Φ_C′ vs `ein_dbar.continue_ellipse` (κ = 1, 0.5, 0.2; both orders), path/slit/monodromy tests, outside-Ω tests, EIN-type route |
| `lem2_k1prime.py` | Regge-box bound on 2d Ising σ and GFF; edge limits; non-uniformity of the true-intercept form |
| `lem3_k2nu_lrc.py` | K_ν vs K̄ (any ν), corner expansion; exact-frame geometry of the x₀-corrections; f-factor in radial variables; diagonal-map test of the stored R_j (sc1, EE1, sc3, EE3) |
| `lem4_collinear.py` | η at the joint corner vs 1/sin(Δφ/2); kernel-bound margins |
