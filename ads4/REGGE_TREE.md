# REGGE_TREE: is H_K valid for a tree-level result? (referee report 3, M2)

*23 Nov 2026. Unreviewed; nothing committed; no existing file edited.*
* *New code:* `ads4/code/reggetree/regge_tree.py` (about 1 CPU-min, 1 process, < 200 MB).
* *Output:* `ads4/code/reggetree/cache/regge_tree.txt` (git-ignored).
* *Literature:* arXiv sources of 2106.10274 (CHMRS), 2311.04271 (CLS-D, v3), 1904.05905 (KKSDZ) and 1209.4355 (CGP) were
  downloaded and read at the passages quoted. Equation and section numbers were taken from the arXiv PDFs.

Grades: [P] proved; [D] derived; [D-sketch]; [DV] derived and checked numerically; [NE] numerical evidence; [M] model;
[A] assumption; [Q] quoted; [I] open. A chain carries the grade of its weakest link.

---

> **Update (INDEP_CHECK2.md, independent reproduction with new code): all five targets reproduced, nothing load-bearing refuted. Two substantive qualifications: Theorem W's J\* = 50 → 64 once the k = 3 piece is included, and Lemma LR∞'s r(β) table uses a cruder constant than the stated formula.**
> * **Magic-block identity:** exact with l and b symbolic to O(r⁶), so it is an identity in b. It is exact to O(r⁸) for every b = −4…4, and holds to 4·10⁻⁴⁸ at 50 digits to r⁴⁰ (including b = 0.37, 7, 10 and complex ν); controls fail at O(1). **The monopole constant D is reproduced to 8 digits, and its t-independence is now derived:** 1/c = 1/32 + (15/2)/(4μ²+1) − (9/2)/(4μ²+9), so D = (15/2)·defect[1/(4ν²+1)].
> * **Lemma LR∞:** reproduced; the true r(β) is ≤ 8·10⁻⁷ on β ∈ [1, 32]. **Discrepancy:** the stated formula (with ∫|ψ|) gives r(1) = 0.016, while the table (0.577) used Σ|c_n|/(n+1). Both are rigorous; the text must name the constant used. With the tighter constant the bound holds down to β = 0.5 (r = 0.28), which would give J\* = 324 at the gap instead of 857 (needs a certified ∫|ψ₃| and a denser grid).
> * **Order of limits (route N):** sound, and the quotes are verbatim or faithful. Minor fixes:
>   * the CHMRS ε quote is in §1, not §1.1.1;
>   * the "technical assumption" is eq. (5.3), not §5.3;
>   * CLS-D restates the logic in §2.2.2 but never discusses the order of limits, so do not say "the logic of CLS-D";
>   * f⁻ = −1.
> * **Theorem W:** the 2×2 criterion d ≤ H is correct and sharp (200 000 random tests). All J\* table entries are reproduced. The Gershgorin factor is ≤ 1.5 only if f₊₊ ≈ f₊₋ (otherwise up to 3; 1.50–1.53 on the samples). **Main qualification: d omits the k = 3 piece of the corrections. Including it moves J\* 50 → 64 at (3, 10, 10), 30 → 33, and 64 → 80**; the J = 50, 51 states fail at d/H ≈ 1.24. This must be stated as a named premise or used as the J\*. The cache stores the 4a/4b shapes with 1/m² instead of 1/m⁴ (documentation slip only).
> * **⟨TTO⟩ counts:** verified against DKKPS's TTTT.tex table and by an independent helicity count (ℓ ≤ 11).

## 0. Verdict

**The verdict is (a), in this sense: the Result never needs a tree-level H_K.**
* The sum rule is applied to the exact correlator at finite c_T. There R1 (J₀ ≤ 1) holds and H_K is [D] as proved, with margin ½ at
  k = 2, 3, 5 and 3/2 at k = 4.
* (Pos) makes the heavy side ≥ 0 at every c_T. (Pos) is a statement about (Δ, J, parity) and does not refer to c_T.
* Only the light side is expanded in 1/c_T: Ω[light](c_T) = L₁/c_T + o(1/c_T). Then L₁ ≤ 0, and that is the Result.

This is exactly the logic of CHMRS §1.1.3, eqs. (1.12)–(1.14), and of CLS-D (§2.1, §3.4). Tree level enters only through the light side,
which is already an assumption (A2 + L7(b), [A]). The heavy side is never expanded.
* It need not be analytic in 1/c_T.
* CHMRS §4.4 show that an eikonalised heavy sector satisfies the k = 2 sum rule while its term-wise O(G) heavy sum diverges.

**Grades are unchanged.** H_K stays [D] and the implication stays [D-sketch]. Two things are needed:
* A2's [A] must be restated so that it includes the uniformity of the light-side remainder (§5).
* The order of limits must be written down. It is written nowhere in the Letter or the Supplement.

**Where the referee is right.**
1. The order of limits is not stated anywhere.
2. "Regge spin k must exceed the CFT's Regge intercept J₀ ≤ 1", placed next to "A2: tree level", invites the term-wise reading.
3. Under that term-wise reading, the uniform corner count fails at J₀^planar = 2 (Table 0, column U).

**Where the referee is wrong or incomplete.**
1. *The term-wise reading is not the logic.* The tree-level heavy sum (single traces above the gap) never has to converge.
2. *The referee missed k = 5.* Under the uniform planar count, k = 5 (Ĉ₅^{ds}, c = 3/2, a = 2) also fails, not only k = 2, 3.
   k = 4 survives only because Q's weights have f_pos(±3i/2) = 0; generic weights (c = 1) fail as well.
3. *The literature attribution is inaccurate.* "CHMRS and CLS-D state a planar Regge-boundedness assumption (J₀ < 2 at tree level)" is
   not what they do.
   * Both apply the non-perturbative sum rule with J₀ ≤ 1 (§2).
   * CHMRS (2.8) is a flat-space S-matrix statement. CHMRS §5.7 derive twice-subtracted flat dispersion relations *from* the CFT sum rules.
   * The planar intercept appears only in KKSDZ §4.5, eq. (4.80) (J₀^planar ≤ 2 from chaos). There KKSDZ also note that the 1/N
     expansion does not commute with the Regge limit.
4. *The intercept quoted for the heavy sum is mislabelled.* "J₀ = 2 (Einstein)" is the planar intercept of the **Wightman function**,
   i.e. of elastic graviton exchange. The **heavy** tree dDisc is empty for pure Einstein gravity.
   * In a stringy Regge-pole model it has fixed-η intercept 2 − 9/(2Λ²), with amplitude ∝ (2 − j). Its k = 2 Regge moment is the
     transverse graviton propagator, ∝ η⁻² at large impact parameter (§3, CHMRS (4.45)).

**Per-k table.** Corner count c < a + 1 − J₀ (HK §1; INDEP_CHECK §2.4). Columns: **N** = finite-c_T route, the one used; **U** = term-wise,
uniform planar bound J₀ = 2; **T** = term-wise, with the tree-level impact-parameter tail of the heavy moment.

| k | piece | a | c | N: needs J₀ < | N: margin under J₀ ≤ 1 | U: needs c < a − 1 | T, fixed η: needs j_tree < | T: available (strings) | T: fixed-η margin | T, corner: needs c < a + 1 | T: corner margin |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 2 | Ĉ′_EE (cut, de-shifted) | 0 | −½ (×log) | 3/2 | **½** | fails (and fixed η diverges) | 2 | 2 − 9/(2Λ²) | 4.5·10⁻⁶ (1000), 5.0·10⁻⁵ (300) [M] | ✓ | 3/2 |
| 3 | Ĉ₃ | 0 | −½ (×log) | 3/2 | **½** | fails (and fixed η diverges) | 2 | same | same [M] | ✓ | 3/2 |
| 4 | Ĉ₄^{ds}, ∂Ĉ₄^{ds} (Q, f_pos = 0) | 2 | ½ | 5/2 | **3/2** | ✓ margin ½ | 4 | 2 − O(Λ⁻²) | 2 | ✓ | 5/2 |
| 4 | generic weights | 2 | 1 | 2 | 1 | fails (borderline) | 4 | — | 2 | ✓ | 2 |
| 5 | Ĉ₅^{ds} (raw) | 2 | 3/2 | 3/2 | **½** | **fails** | 4 | 2 − O(Λ⁻²) | 2 | ✓ | 3/2 |

Grades:
* **N:** [D]; its inputs are as in Table S2 (R1 [Q], R1′ [D + DV]).
* **U:** [P] (algebra) with [DV] (`regge_tree.txt` §C).
* **T, fixed-η:** [A] (needs J₀^planar < 2, beyond KKSDZ (4.80)), with [M] for the stringy value.
* **T, corner:** [D-sketch], given the CHMRS (4.45) tail, which is itself derived from the sum rule. As a justification of a term-wise
  sum rule this is circular, so the T column is [A] overall. In the Regge-pole model it is [P + DV].

**What is lost?** Under the finite-c_T route, nothing.
* If one insisted on the term-wise route with only the uniform bound, item (1) and the k = 2, 3, 5 parts of Ω^Q_up would lose their
  derivation; the k = 4 part (Q weights) would not.
* No modified functional rescues the k = 2 parts. The graviton pole lives at k = 2, so dropping k = 2 removes item (1).
* Making c < −1 at k = 2 would mean changing the de-shift target, i.e. re-optimising the functional [I]. It is unnecessary.

---

## 1. Trace: what J₀ in the H_K criterion is

1. **Statement (supplement.tex:483-487; HK_CLOSE §1).** R1′/K1′ is a uniform bound on the Regge box
   B_c = {ρ, 1/ρ̄ ≤ c}: |dDisc 𝒢| ≤ B r^{1−J₀}.
   * It bounds the double commutator of the **exact, full** correlator, in the Rindler-reflection-positive polarisation configuration.
   * It is KKSDZ's J₀^dDisc, eq. (4.70), J₀^dDisc ≤ 1, obtained from Rindler positivity and Cauchy–Schwarz (§4.2.3–4.2.4, (4.55)).
2. **How it reaches the heavy sum.** R3 (per-block Rindler positivity) gives Σ_O ∫|K||dDisc G_O| ≤ ∫|K| dDisc 𝒢.
   * So J₀ is the growth of a *majorant* of the heavy sum: the full correlator, light parts included.
   * The count is ∫dη η^c ∫₀^{c₀/2η} r^{a+1−J₀} dr, which converges iff c < a + 1 − J₀ (INDEP_CHECK §2.4, reproduced).
   * Equivalently (HK §1, "heavy-measure form"), ∫e^{γβ}dμ_k < ∞ for γ < a + 1 − J₀, where β is the AdS impact parameter:
     e^β ≈ (Δ+J)/(Δ−J), as in CHMRS §5 (β ≈ log((Δ+J)/(Δ−J))).
   * In CHMRS's language the r-integral up to r_max(η) = η − √(η²−1) is exactly their k-th Regge moment Π_{k,η} (CHMRS (2.39)). The
     "corner" is therefore the large-impact-parameter tail of the heavy Regge moments.
3. **Which J₀.** It is a property of the correlator *at whatever c_T the sum rule is applied*.
   * Applied to the exact correlator (route N), J₀ ≤ 1 [Q] and H_K is [D].
   * Applied to the O(1/c_T) correlator, the relevant majorant is dDisc𝒢^{(1)}.
4. **The tree-level heavy sum alone.** At O(1/c_T) and k = 2, 3 the light s-channel part is the T block plus light single traces. Double
   traces sit at double zeros (L1). This is a finite sum and does not grow, so the heavy tree sum and dDisc𝒢^{(1)} have the same Regge
   behaviour.
   * **Fixed η.** Conformal Regge theory (CGP §4, eq. (56): 𝒜 ≈ 2πi∫dν α(ν)σ^{1−j(ν)}Ω_{iν}(ρ)) gives intercept j(0).
     * j(0) = 2 for strict Einstein, where the scattering is elastic: "In this strict λ → ∞ limit the spin j(ν) = 2 and the scattering is
       elastic" (CGP §5.2.1). The heavy dDisc then vanishes and the heavy sum is empty.
     * j(0) = 2 − 9/(2Λ²) for the stringy trajectory j(ν) = 2 − 2(ν²+9/4)/(ML)², with ML ≃ Λ (cf. N = 4: j(ν) = 2 − (4+ν²)/(2√λ), CGP §5.2, eq. (95)).
     * The inelastic (dDisc) part carries a factor ∝ 2 − j(ν). This is what makes its k = 2 moment O(1) and equal to the graviton
       propagator.
     * Our weights reach ν → 0 (ψ₂(p) ~ p, the IR log), so the fixed-η margin is the minimum of 2 − j(ν): 9/(2Λ²) = 4.5·10⁻⁶ at
       Δ_g = 1000.
   * **Corner / large impact parameter.** The k = 2 heavy moment equals the transverse graviton propagator h₂(η−1) ∝ η⁻²
     (CHMRS (4.45), for "any holographic CFT", from the sum rule).
     * In ν-language, the tail is set by the first singularity of α(ν) off the real axis. That is the stress-tensor pole at ν = ±3i/2,
       the same point where f_pos(±3i/2) = 0 matters for k = 4 (EESTEP §0).
     * Along the shifted contour, Re j(x − iμ) ≤ 2 for μ < 3/2, so the effective corner exponent is j − ½ − μ → 0.
     * Hence the corner criterion becomes c < a + 1, as in column T.
5. **"J₀ = 2 < 4" (LEMMAS §1.4; RESULT_Q_STAGE3 §1.5).** This is the growth of the *light* tree correlator (Einstein exchange). It is
   used for the k = 4 light-side lemma (Lemma A, D_LR), where its impact-parameter profile η⁻² is used explicitly. It is a light-side
   statement and consistent with the above.

## 2. Literature: what CHMRS, CLS-D and KKSDZ actually assume

**CHMRS (2106.10274).**
* §1.1.1: "In weakly coupled AdS EFTs, where tree-level interactions are suppressed by ε (which we assume to be much smaller than any power
  of 1/Δ_gap)" (§1, p. 3). Also: "we make no assumptions about the primaries with τ > Δ_gap, other than that they obey the usual unitarity
  constraints. While the physics below Δ_gap is assumed to be weakly coupled, the physics above Δ_gap can be strongly coupled."
* §1.1.3, (1.12)–(1.14): ω|light = Σ_{τ≤Δgap}, ω|heavy = Σ_{τ>Δgap}, "−ω|light = ω|heavy". "It is easy to see that ω|light = O(ε) …
  If we manage to construct a functional ω that is non-negative for τ > Δgap, we conclude from (1.13) that −ω|light ≥ 0." The heavy side
  is used only through its sign.
* §1.1.3: "it is fully rigorous … boundedness in the Regge limit (with Regge intercept ≤ 1) is a consequence of unitarity and the OPE at
  the non-perturbative level [Caron-Huot 2017]."
* §4.1: "We assume that the light OPE data agrees with a tree-level AdS effective field theory (EFT) with a derivative expansion." This
  is our L7(b) [A].
* §4.4, (4.45): "a universal prediction for the k = 2 moment of the double-discontinuity of heavy states in the Regge limit in any
  holographic CFT": 8πG (vol S^{d−2}/2^{d−2}) h_{d−1}(γ−1) = ⟨(2/m⁴)δ(γ−η)/…⟩(1 + O(1/Δgap)).
* §4.4, after (4.47), on the eikonal model: "Naively expanding the integrand in small Γ_tree would result in a divergent integral
  ∼∫dΔ Δ. Instead … we do not expand the exponentials … It is a nice consistency check … that the eikonal model satisfies the sum
  rules." The heavy side need not have a convergent O(G) term; CHMRS never expand it.
* §2.1, (2.8): "lim M(s,u)/s² = 0 (u < 0)", spin-2 convergence of the **flat-space S-matrix** (review).
* §5.7: "The present paper establishes a slightly stronger, nonperturbative result: … twice-subtracted dispersion relations hold. … The
  argument is that dispersive CFT sum rules with Regge spin-2 converge when acting on any physical correlator."
* Their own loophole, §2.4 footnote 9: "(G − 1) vanishes in the u-channel collinear limit even on the second sheet, which has never been
  rigorously proven." This is the analogue of KKSDZ §4.2.6.
* §5.3 (conclusion of the eq. (5.3) discussion) is the "technical assumption" on towers of couplings, our L7(b)'s "no conspiracy".

**CLS-D (2311.04271 v3).**
* §3.1, below (3.1): "we also need the condition J₁+J₃−1 > J₀, where J₀ is the Regge intercept [KKSDZ]."
* §3.4: "for non-perturbative CFTs we only know that J₀ ≤ 1 [17, 21, 24]"; "the condition from convergence in the Regge limit becomes
  J₀ < k, which is indeed true for any k ≥ 2."
* §3.5.2: "we must have k > J₀ for even j′ and k − 1 > J₀ for odd j′. For non-perturbative CFTs, we have J₀ ≤ 1. So, we should choose
  k ≥ 2 for even j′ and k ≥ 3 for odd j′." For J_e = 2 the unsubtracted factors are k = 2 (f⁺) and k = 3 (f⁻ = 1). This is the origin
  of a = 0 for both k = 2 and k = 3.
* §2.2.2: "We leave an exploration of the Regge regime for four-point functions of spinning operators to future work."
* No planar intercept appears anywhere. The large-N application follows CHMRS (CLS-D §1, footnote 1: "One must additionally check other regimes
  to derive a rigorous CFT bound").

**KKSDZ (1904.05905).**
* (4.55), (4.70), (4.73): J₀ ≤ 1 and J₀^dDisc ≤ 1 non-perturbatively. §4.2.6 assumes the asymptotic light-cone expansion (4.66) ("we
  don't have a general proof").
* **§4.5, (4.80):** "In large-N perturbation theory, the bound on chaos implies that J₀^planar ≤ 2 … saturated in holographic theories
  with a large gap, where J₀^planar = 2 comes from tree-level graviton exchange in the bulk."
* §4.5, (4.81)–(4.83): "The four-point function has an expansion in powers of 1/N at fixed values of cross-ratios, but this expansion
  does not commute with taking the Regge limit … term-wise integration is reliable for the orders at which it converges."
* §5.3: the planar superconvergence sum rules are applied term-wise (single traces). This is justified there because their
  [ℰ, ℰ] needs J₀ < 3 (§4.4: "3 > J₀ = J₀^dDisc"), and J₀^planar = 2 < 3.
* For our k = 2, 3 functionals the term-wise condition would be J₀^planar < 2, which (4.80) does not give. So KKSDZ's large-N
  discussion does **not** cover a term-wise use at k = 2, 3. It does support route N: exact sum rule at finite N, expand only what
  converges.
* §4.2.6's caveat (first/second-sheet boundary) is unrelated to M2. It concerns R1 at finite c_T and lies outside B_c (Table S2, R1).

**CGP (1209.4355).**
* §4 (below (55)): "At leading order in 1/N the CFT is described by the dual AdS theory at tree-level, which implies the absence of poles in the
  complex J-plane from multi-particle states."
* §5.2, eqs. (95), (98): j(ν) = 2 − (4+ν²)/(2√λ)(1 + …); j(0) = 2 − 2/√λ − 1/λ + …
* §5.2.1, graviton exchange: α(ν) ∝ 1/(4+ν²) (the d = 4 analogue of our 1/(ν²+9/4)).

**Bottom line.** Every one of these papers uses the non-perturbative J₀ ≤ 1 for the sum rule and treats 1/c_T only on the light side. No
paper states a tree-level heavy-sum Regge assumption for k = 2.

## 3. The term-wise route, quantified (not used; for the record)

**Model [M]** (`regge_tree.py` §B). D_ε(r, β) = (ε/r)e^{−2εℓ}K_{εℓ}(β), with ℓ = log(1/r), ε = 2/Λ², and K the heat kernel on H₂. This is the
conformal-Regge form with j(ν) = 2 − ε(ν²+9/4), an inelastic residue ∝ 2 − j(ν), and the normalisation that reproduces the graviton
pole.

* **Positivity [P].** D_ε ≥ 0, because e^{−εℓ(ν²+9/4)}/(ν²+9/4) = ∫_{εℓ}^∞dτ e^{−τ(ν²+9/4)}.
* **Fixed-η intercept → 2 − 9ε/4 [DV].** The local value is (2 − J)/ε = 9/4 + 3/(2τ) + …, measured 2.2546 at τ = 320. It is reached only
  at ℓ ≳ 1/ε ≈ 5·10⁵ (Δ_g = 1000).
* **k = 2 Regge moment [P + DV].** Π₂(η) = ∫_{εℓ₀}^∞ e^{−2τ}K_τ(β)dτ ≤ G₂(β) = Q₁(cosh β)/2π ~ 1/(6πη²).
  * The heat-kernel identity is checked to 10⁻¹² (β ≤ 4) and 4·10⁻⁷ (β = 8).
  * Π₂/G₂ = 0.980 (ε = 0.05, β = 0.5) and 1.00000000 (β ≥ 6, or ε = 2/Λ²).

  So the model obeys CHMRS (4.45) with tail η⁻².
* **Corner integrals [DV]** (§C; η up to 10⁸).
  * Majorant N (J₀ = 1): converges for all k; margins ½, ½, 3/2, ½.
  * Majorant U (J₀ = 2): diverges at k = 2, 3 (already at fixed η, r → 0) and at k = 5 (I = 2.1, 25, 250, 2500 at η_max = 10², 10⁴,
    10⁶, 10⁸). It also diverges at k = 4 with generic weights (logarithmically). It converges at k = 4 with Q's weights.
  * Tree moment T: converges for all k, with margins 3/2, 3/2, 5/2, 2 (generic k = 4), 3/2.
* **Einstein limit.** As ε → 0, D_ε → 0 pointwise while Π₂ → G₂ ≠ 0. Pure Einstein gravity, with no heavy states, cannot satisfy the tree
  k = 2 sum rule: the graviton term is unbalanced. This is the CEMZ / CHMRS statement, not a defect of H_K.

**Term-wise H_K in assumption form.** It would hold under the following assumption:

> **(TR)** "The O(1/c_T) double commutator in the Rindler-positive configuration has finite k = 2 Regge moments Π₂(η) for every η (that
> is, j_tree(ν) < 2 on the real ν axis), with Π₂(η) ≤ Cη^{−p}, p > ½." [A]

* Under (TR) the margins are p − ½ at k = 2, 3, 5 and p + ½ at k = 4.
* p = 2 is CHMRS (4.45) and the Regge-pole model.
* (TR)'s fixed-η part is the AdS image of CHMRS (2.8) restricted to the heavy sector. Its corner part is a planar analogue of LC\*,
  which LRC_CLOSE §5 graded [A], beyond the literature.
* We do not recommend stating (TR). Route N needs neither part.

## 4. The order-of-limits argument (route N) in full

1. *(Exact sum rule.)* For every finite c_T, H_K holds for the exact ⟨TTTT⟩ [D, R1–R3 etc.]:
   Σ_{Δ−J ≥ Δ_g} Tr[P_O M_Ω(O)] + Σ_{light s.t.} Tr[P_O M_Ω(O)] = −Ω[light](c_T).
   * The heavy sum includes multi-trace states above the gap; spins 1 and 3 cannot occur at any c_T.
2. *(Pos.)* M_Ω(Δ, J) ⪰ 0 on every heavy (Δ, J, parity), and on light J = 0, 2 single traces. Hence Ω[light](c_T) ≤ 0 for every c_T.
3. *(Light side.)* Ω[light](c_T) is the action on the identity, T and states with Δ − J < Δ_g other than light single traces.
   * **k = 2, 3.** It is c_T⁻¹ × (the T-block action with the exact ⟨TTT⟩ coefficients), plus light double- and multi-trace terms.
     * Those terms are O(1/c_T²): the double zeros give sin²(πγ/2), and λ²_{TT[φφ]}, λ²_{TT[Tφ]} = O(1/c_T²).
   * **k = 4, 5.** The first double-trace families enter at O(1/c_T) through their tree data (L4–L7), with O(1/c_T²) corrections.
   * The large-spin sums of these remainders converge: γ_{n,J} ~ e^{−β}·Δ_g^p/c_T on light families, against images e^{cβ} with c ≤ 3/2
     and the a-suppression [D-sketch].
   * Hence Ω[light] = L₁/c_T + O(Δ_g^p/c_T²), given A2 as restated in §5 [A].
4. *(Limit.)* L₁ ≤ c_T·O(Δ_g^p/c_T²) → 0. L₁ ≤ 0 is the Result, with L₁ the tree-level light side computed in L1–L8.

No step uses the O(1/c_T) heavy sum, its convergence, or J₀^planar.

**Side remark on M4(b).** Step 3 shows that for item (1) the inequality constrains the T-block coefficients themselves (the exact t₄, α_T).
The error is the O(1/c_T²) light remainder relative to the O(1/c_T) T term [D-sketch].
* That remainder contains one-loop light-sector pieces (for example [φφ] double traces).
* We have not determined the power p, so we cannot say how large c_T must be for the 10⁻⁸ statement to hold for the exact t₄.
* The referee's "at leading order in 1/c_T" is the safe wording; keep it.

## 5. Proposed wording

**letter.tex:109 (A2).**
> \item[A2] (weak coupling below the gap) $c_T\to\infty$ at fixed $\Delta_g$: the OPE data of states with $\Delta-J<\Delta_g$ (other than
> light single traces) are those of a tree-level bulk effective theory up to $o(1/c_T)$, uniformly in spin. The sum rules are used at
> finite $c_T$; heavy states are never expanded in $1/c_T$.

**letter.tex:133-135 (Mechanism).** Replace "(… its Regge spin $k$ must exceed the CFT's Regge intercept $J_0\le1$~\cite{KKSDZ,CaronHuotInversion}), the total vanishes." with:
> ($k$ must exceed the Regge intercept $J_0$ of the exact correlator; $J_0\le1$ in every unitary CFT~\cite{KKSDZ,CaronHuotInversion}),
> the total vanishes. We use it at finite $c_T$, as in Refs.~\cite{CHMRS,CLSD}: only the light part is expanded in $1/c_T$, and the heavy
> part enters only through its sign. The planar intercept ($J_0=2$ for graviton exchange~\cite{KKSDZ}) never enters.

**letter.tex:183-185 (Functionals).**
> The interchange of the OPE sum and the integrals is absolutely convergent for $J_0<\tfrac32$ ($k=2,3,5$) and $J_0<\tfrac52$ ($k=4$)
> [D], hence at every finite $c_T$, with margins $\tfrac12$ and $\tfrac32$. The bounds are the $c_T\to\infty$ limit of the resulting
> inequality: its heavy side is non-negative at every $c_T$ by (Pos), and its light side is $O(1/c_T)$ (A2). Applied term by term to the
> $O(1/c_T)$ correlator, the same count would need a planar intercept below $\tfrac32$, which fails; we never do this.

**supplement.tex, Table S2.**
* **Row A2:** "$c_T\to\infty$ at fixed $\Delta_g$; light data ($\Delta-J<\Delta_g$, except light single traces) = tree-level EFT $+\,o(1/c_T)$
  uniformly in spin; sum rules used at finite $c_T$, heavy side never expanded | D, Q | [A]".
* **Row R1:** append "for the exact correlator at every $c_T$; the planar intercept ($=2$, KKSDZ (4.80)) is not used".
* **Row $H_K$:** append "at every finite $c_T$ (order of limits: Sec.~\ref{sec:HK})".

**supplement.tex, end of Sec. \ref{sec:HK}** (after "Conclusion"), new paragraph:
> \emph{Order of limits [D given A2].} $H_K$ is used only for the exact correlator at finite $c_T$, where R1 holds. (Pos) is a statement
> about $(\Delta,J)$ and parity, so it applies at every $c_T$ to every state with $\Delta-J\ge\Delta_g$, single- or multi-trace, and to
> light single traces; the sum rule then gives $\Omega[{\rm light}](c_T)\le0$. At $k=2,3$ the light side is $c_T^{-1}$ times the action
> on the $T$ block with the exact $\langle TTT\rangle$ coefficients, plus light double and multi traces, which enter at $O(c_T^{-2})$
> (double zeros). At $k=4,5$ the first double-trace families enter at $O(c_T^{-1})$ with their tree-level data (L4–L7). Given A2,
> $\Omega[{\rm light}]=c_T^{-1}L_1+o(c_T^{-1})$, so $L_1\le0$, which is the Result. The heavy side is never expanded and need not be
> analytic in $1/c_T$: in an eikonalised heavy sector the $k=2$ sum rule holds while its term-wise $O(G)$ heavy sum diverges
> \cite[Sec.~4.4]{CHMRS}. This is the logic of \cite[Sec.~1.1.3]{CHMRS} and \cite{CLSD}.
>
> \emph{Term-wise use (not made).} Applying the sum rules to the $O(1/c_T)$ correlator itself needs two things:
> (i) $J_0^{\rm planar}<2$ at $k=2,3$ ($<4$ at $k=4,5$), which Einstein gravity saturates \cite[eq.~(4.80)]{KKSDZ} and which a stringy
> leading trajectory satisfies by $2-j(0)\simeq9/(2\Lambda^2)$ [M];
> (ii) in the corner, $c<a-1$ under a uniform planar bound, which fails at $k=2,3,5$.
> With the tree-level large-impact-parameter tail of the heavy $k=2$ Regge moment, $\propto h_2(\eta-1)\sim\eta^{-2}$
> \cite[eq.~(4.45)]{CHMRS}, the corner margins would be $\tfrac32$ ($k=2,3,5$) and $\tfrac52$ ($k=4$).

**Response to the referee (M2), one paragraph.**
> We thank the referee. The sum rule is applied to the exact correlator at finite $c_T$, where $J_0\le1$ and $H_K$ holds as proved; (Pos)
> is $c_T$-independent and is used only for the sign of the heavy side; only the light side is expanded in $1/c_T$. This is precisely
> CHMRS's procedure (their eqs. (1.12)–(1.14), and Sec. 4.4 on the eikonal heavy sector). No tree-level Regge bound on the heavy sum is
> needed or used. The term-wise reading would indeed fail the uniform corner count at $k=2,3$, and also at $k=5$. We now state the order
> of limits explicitly (A2; Supplement, Sec. on $H_K$) and have reworded the two sentences.

## 6. Open points

* **The light-side remainder.** Its uniformity in spin (§4 step 3) is [D-sketch] inside [A]. Its size (the power p) is not computed. Only
  the c_T → ∞ statement is affected; the quantitative M4(b) claim is.
* **(TR) is not established** [A], and is not needed.
* **Whether de-shift targets with c < −1 at k = 2 exist** with comparable bounds is open [I], and not needed.
* **The finite-c_T heavy spectrum includes multi-trace states.** (Pos) evidence is (Δ, J)-generic, so nothing changes. The Letter should
  say "every state with Δ − J ≥ Δ_g" rather than "heavy single traces" wherever the latter appears.

## 7. Reproduce

    cd ads4/code/reggetree
    ulimit -v 3000000; OMP_NUM_THREADS=1 python3 regge_tree.py     # ~1 min -> cache/regge_tree.txt
