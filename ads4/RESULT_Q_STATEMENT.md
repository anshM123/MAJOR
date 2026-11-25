# Result Q: two-sided bound on the R⁴ coupling `g₄` (consolidated statement, 22 Nov 2026)

*Consolidation of 22 Nov 2026. This is a writing pass with no new computation; every number is quoted from the source named next to
it. This file replaces the version of commit `17dc133`, whose eight stacked update boxes superseded its own body. That version is in
git history; Appendix A summarises the sequence.*

**Status.** This is a conditional result.
* The implication "per-state positivity of the two fixed functionals ⇒ the bounds" is [D], given the analytic inputs of §4.
* Per-state positivity is computed [NE + M] on the binding and near-binding channels, and assumed elsewhere.
* The decisive below-fold input, `H_RS` for the binding state at `ν < 6`, is **[I]**.

**Sources and precedence** (later wins):
1. the previous `RESULT_Q_STATEMENT.md` (body and boxes);
2. `RESULT_Q_STAGE0–7`, `RESULT_Q_REVISION`, `RESULT_Q_V4PRIME/V4DPRIME`;
3. `RESULT_Q_V4FINAL` (V4F);
4. `REVIEW_RESULTQ_2026-11`;
5. `RESULT_Q_FWD_CLOSE` (FWD);
6. `REVIEW_UNIFIED_2026-11` (RU; authoritative on grades);
7. `HRS_CHECK` (HRS);
8. `RESULT_Q_OG_CHECK` (OG);
9. `RESULT_Q_FERMION_HEFT` (FH).

Note: HRS, OG and FH are later than RU and refine three of its items. Each such refinement is marked.

**Grades.**
* **[P]**: proved. **[D]**: derived. **[D-sketch]**: derivation sketched. **[DV]**: derived and checked numerically.
* **[Q]**: quoted. **[NE]**: numerical evidence. **[NE + M]**: computed in the chain/HYB model. **[M]**: model.
* **[A]**: assumption. **[I]**: open.

---

> **Update (HRS_CLOSE.md, unreviewed): H_RS for Q's binding state [I] → [D-sketch], i.e. [D] given premise RD.** The argument is analytic and covers every ν window, including ν < 1.
> - At every order of the 1/Δ expansion, the engine's integrand is one fixed multivalued function times a polynomial. Integration by parts reduces every term to a rational multiple of one base integral, the same on any cycle (s-/φ-odd terms vanish; there are no α-odd terms).
> - So the crystal cycle = c_Γ(ν) × engine order by order, and c_Γ cancels in the state/reference ratio.
> - The hypotheses of the argument are checked on the engine's pieces to a⁴–a⁵ [DV]: lattice-shifted exponents shared by scalar/EE2/EE4, no interior poles, no integer exponents, no boundary layers.
> - **Remaining premise RD** (the crystal integral is Regge-dominated, with the far region e^{−O(Δ)}): [DV] for the scalar, [NE] for EE at ν ≥ 10, [D-sketch] for EE below 10.
> - **New K = 5 images:**
>   - truncation is ≤ 0.01 c₂ per harmonic (300), with |Δs_d| ≤ 0.015 (300) / 0.0015 (1000);
>   - the binding FF/N is unchanged (+0.1427 / +0.1428 at δ = 0.2), with the margin ≥ 85× the bound.
> - **Covered:** ++ J = 0, 2. OO J = 0 is covered via Γ_OO = Γ_++ (direct application [D-sketch]). **Not covered:** Result D's (336, 28/32) window, where the obstacle is the J²/Δ resummation, not the cycle.

> **Update (RESULT_Q_EINSTEIN_CHECK.md, unreviewed): Lemma A + D_LR stays [D]; the Einstein `⟨TTTT⟩` itself is still not evaluated [I].**
> * **New evidence [DV].** The mechanism was checked on interacting spinning correlators: `T′ = σ ⊗ T_free`, where `G_σ` is an exact sum of AdS₄ Witten diagrams.
>   * The diagrams: contacts with `J₀ = 0, 2, 4`, and a massive spin-2 exchange (`τ = 2`, `J₀ = 2`) in the s, t and u channels.
>   * The pointwise zero holds at `Ψ₄` and `Ψ₅` to 10⁻¹⁴–2·10⁻¹², and in the collinear approach down to Δφ = 0.025.
>   * Controls fail at O(1).
> * **Not exercised:** the graviton's `τ = 1`, its own tensor structures, and the exact Einstein coefficients at external Δ = 3.
> * **Definition of the light side (clarified, [D-sketch]).** The light side is **D_LR on the `J₀ < 4` terms plus term-wise (F)+(P) on the `J₀ ≥ 4` terms**. It is not D_LR applied to the full correlator.
>   * In the test family the swapped-form zero holds whenever `J₀ ≤ k`, because the joint Regge corner then has an integrable exponent `k − 1 − J₀`.
>   * So the R⁴ value `−Ĉ₄ = 2g₄` cannot sit in the fixed-angle integrand. Applied to the full correlator, D_LR would give it 0.
>   * This is how `g̃₄` is already defined (§1 table); the Einstein zero is unaffected.
> * **Weight control qualified.** The light-transform-weight control is one-sided: lower weights also give 0. It therefore tests convergence of the light-ray integral, not light-transform covariance (RESULT_Q_OG_CHECK §3.4).

> **Update (REPRO_ENGINE.md, independent): the engine ratio ρ(ν) and Q's binding `s_d` are reproduced by a different route, with no discrepancy beyond stated precision.**
> * **Route.** A direct 7-D quadrature of the exact pointwise integrand on D at finite Δ. It uses no 1/Δ expansion, and the conserved ⟨TTO⟩ coefficients were re-derived symbolically.
> * **ρ agreement:** `≤ 3·10⁻⁷` (Δ_g = 300) and `≤ 3·10⁻⁸` (Δ_g = 1000) relative to the engine's K = 5 ratio on ν < 4.
> * **`s_d`:** `+0.50967` / `+0.50776`, against `+0.5096` / `+0.5077`.
> * **FF/N:** `+0.14265` / `+0.14282`.
> * **Sensitivity.** Reaching the fail point `s_d = 1.776` needs a worst-sign shape error of `2.1·10⁻⁴` / `1.9·10⁻⁵` in ρ, i.e. ≈ 2.5 % of its beyond-leading-Regge deviation. That is ~600× the measured direct-vs-engine difference.
> * **New grade.** The engine's algebra and truncation for the binding state are [DV] on D.
> * **Not tested:** premise RD (crystal cycle ≡ real section D), because both routes integrate over D. H_RS stays [D-sketch].
> * **Record note.** Production K = 4 is missing some a⁴ pieces (HRS_CLOSE §3.5). Its effect on `s_d` is +0.0004, inside ±0.015.

> **Update (RD_CLOSE.md, unreviewed): premise RD at ν < 10 is reduced to one cycle-class premise P_cyc, which EE shares with the scalar.**
> * **EE = a·Q·scalar exactly [P + DV].** Q is single-valued rational, every pole of Q is already singular for the scalar, and Q is regular at m₁₄ = 0 and at s → ±1.
> * **Far-region suppression on D′** (the real section cut at the zeros of B): `O(Δ^N e^{−cΔ})`, with `c = 2 artanh(r/2)` (≥ 0.336) [D]. The integrated check gives fitted exponents within 2 % of analytic for EE and scalar [DV].
> * **On the actual crystal thimbles** the far weight at r ≥ 0.3 is pointwise `≤ e^{−66}` (Δ_g = 300) and `≤ e^{−277}` (1000). After the thimble volume factor (≲ e^{40} at 300) this is `≲ e^{−26}`, against failure thresholds `e^{−8.5}` / `e^{−10.8}` [NE + D-sketch]. The r ≥ 0.2 figure (e^{−37}) does not survive the volume factor and must not be quoted (REFEREE_REPORT M10).
> * **Net:** RD for EE at ν < 10 is [D given P_cyc], and H_RS for Q's binding state is [D given P_cyc].
> * **What P_cyc says:** the physical cycle lies in the class of D′. For the scalar it is supported through its periods (engine = exact skeleton).
> * **Side finding.** At ν ≲ 2 the scalar thimble has an angular end on which EE does not decay. This plausibly explains HRS_CHECK's signal loss below ν ≈ 2 [D-sketch].

> **Update (C1_COMPUTE.md, unreviewed): c₁ is one uncomputed universal constant, not a free parameter. Exact (non-linearised) bounds hold for every c₁.**
> * **What vanishes, in H7's prescription:**
>   * the k = 3 part [P given Q2, Q8, F];
>   * the u-channel E×R³ T-exchange Regge residue [D]. At O(t₄) it is helicity-flip, and EE kills every flip configuration; checked across all 16 polarisation assignments.
> * **What c₁ is.** The remaining contact-like piece is ν- and Δ_g-independent [D, given J₀(E×W³) = 2 (D-sketch) and polynomiality].
>   * Its value is [I].
>   * Its natural size is ~10⁻³–10⁻¹ [M].
>   * No bound on it exists [I].
> * **Cubic bound, exact:** `√(t₄²+α_T²) ≤ R₀(x + √(1+x²))`, with `x = 2.160·10⁻³|c₁|` (Δ_g = 300) or `1.670·10⁻⁴|c₁|` (1000) [D].
>   * |c₁| = 100 gives 3.88·10⁻⁶ / 2.94·10⁻⁸.
>   * |c₁| = 1000 gives 1.42·10⁻⁵ / 3.41·10⁻⁸.
> * **Quartic upper end, exact (given κ > 0):** `2L_AdS + 6.81·10⁻⁵c₁²` (300) / `+ 5.44·10⁻⁷c₁²` (1000), taking the supremum over unknown t₄, α_T [D].
>   * It no longer depends on Result D, which removes RU issue 14.
>   * **The `0.57|c₁|` / `0.059|c₁|` figures used above are superseded.**

> **Update (LIGHTSIDE_DIST.md, unreviewed): M3 resolved structurally. Result Q's light side is ONE functional, and it is the one the heavy blocks use [D].**
> * **The functional.** `Ω_k` is CLS-D's term-wise swapped double-commutator form on the s-channel OPE, with no separate collinear or ν-boundary terms.
> * **Decomposition, by Cauchy in Mellin space:** `Ω_k = D_LR + 𝒜_k`, where `𝒜_k` is the Regge arc [D].
>   * The arc is 0 for Regge spin `n < k`; that is Lemma A, now read as "no arc below the borderline" [DV].
>   * At `n = k` the arc is fixed by the top Regge coefficient alone (P) [D + DV].
> * **Checks on EIN's Witten family [DV]:**
>   * term-wise values = whole-function values to 1e-11–3e-5;
>   * the borderline contact gives a finite light value, e.g. −8.7642e-3 at P1, equal to its top-Regge term to ≤ 7e-7;
>   * a borderline exchange (the H7 analogue) matches its top-Regge term to 8e-8.
> * **Consequences:**
>   * The piecewise "D_LR on J₀ < 4 + (F)+(P) on J₀ ≥ 4" light side is now the *derived* value of Ω₄.
>   * `g̃₄ = −½·𝒜₄[light] = −(heavy sum)` given H_K [D].
>   * Result D's H7 pole is the arc of the same functional at `n = k = 2` [D].
>   * No number changes.
> * **Corrections to the Einstein-check box above:**
>   * Its claim that J₀ = k values sit at a = b or in the ν integral is **withdrawn**.
>   * Its "J₀ = 4 contact under Ψ₄" is subcritical in that product family (effective n = J₀ − 2), so it is not an R⁴ control.
> * **Remaining:**
>   * (F), the normalisation `−𝒜₄[R⁴] = 2g₄` (now [D], see the F_NORM box);
>   * the primaries-vs-Mellin-residue rearrangement [D-sketch];
>   * a = b integrability against C_ν [I].

> **Update (HK_CLOSE.md, unreviewed): H_K no longer needs J₀ < ½. With an analytic cutoff on the de-shift weights it follows under the known bound J₀ ≤ 1, with margin ½.**
> * **No-go [P].** Any de-shift weight supported on a bounded ν-range, including every stored grid weight, regrows at the raw Regge rate by Paley–Wiener. The interchange then needs J₀ < ½, the old H_ip.
> * **The fix: Proposition K_c [D given K1′, K2ν, LRC, K3, DS].**
>   * Each de-shift weight is multiplied by an entire erf cutoff c(ν): ν_c = 150, a = 6 for Q; ν_c = Λ+100, a = 10 for D.
>   * The de-shift identity forces strip-analyticity [D], so the cutoff does not reintroduce growth.
>   * Tonelli then gives an absolutely convergent heavy sum for J₀ < 3/2 (k = 2) and J₀ < 5/2 (k = 4).
>   * Checks: exactly solvable toys [DV]; real Q weights [NE].
> * **Every computed number is unchanged.** ε_c ≲ e^{−180} (Q) and ≲ 10⁻⁴¹ (D).
> * **Ĉ₅^{ds} and Ĉ₃ need nothing new** (plain Fubini).
> * **Net grade.** H_K stays [D-sketch] by weakest link, but it is no longer an independent hypothesis. It reduces to K1′ (existing), K2ν (new: per-ν kernel growth ≤ e^{γν} away from the corner, [D-sketch]), LRC and DS.
> * **Required change.** The cutoff must be added to the functional definitions (§1 above; RESULT_D_STATEMENT §2; RESULT_Q_STATEMENT).
> * **Side note.** The C5 de-shift is poor beyond t ≈ 15, which affects per-state positivity at β ≳ 15 (the HYB region).

> **Update (F_NORM.md, unreviewed): (F) [D-sketch] → [D]. The normalisation `−𝒜₄[R⁴] = 2g₄` (= `−Ĉ₄^{ds}|_{R⁴}`) is confirmed with factor 2, and the all-plus `ĝ₄` gives 0.**
> * **Sign correction.** "𝒜₄[R⁴] = 2g₄" in the LIGHTSIDE box above and in LIGHTSIDE_DIST §4 drops a minus; the correct form is `𝒜₄[R⁴] = −2g₄`.
> * **Method: calibration.** An auxiliary heavy scalar coupled through `χ(W⁺)²` has Regge spin 3, below the borderline. Lemma A therefore makes the arc of the R⁴ contact it generates equal minus its block value. That value is `2g₄^χ` by the flat sum rule; the identity is exact and independent of m.
> * **Grade** [D], given Lemma A (term-wise), ν-paired H_EFT(a) at leading order, and the CLS-D J = 0 heavy dictionary, all already on Q's list.
> * **Checks [DV]:**
>   * CHMRS contact functional reproduced exactly;
>   * Penedones top coefficient vs the heavy-density coupling, extrapolated ratio 0.9999990;
>   * calibration identity to 5e-10;
>   * Lemma A on heavy σ-exchanges, ≤ 1.7e-4;
>   * (P) polynomial ratios to 1e-7.
> * **Negative finding.** At fixed celestial angle the calibration fails (an O(1) residue). H_EFT(a) is false pointwise in angle and must be used in ν-paired form.
> * **Net.** `g̃₄ = g₄ + O(Λ⁻²)`, given H_EFT(b) [A]. This is CHMRS's standing assumption (2106.10274 §5.1) and is also assumed by CLS-D. No numbers change.
> * **(P), the g₅ cancellation: [D-sketch], with tension.** Exact cancellation requires the ++ J ≤ 3 heavy-image correction to vanish as m → ∞, but the measured s_d ≈ 0.51 is m-independent. The consistent reading is RU's double counting, which gives lower slope 0.255 (CHLPS 0.26).
>   * The δ = 0.2 / 0.05 upper bounds survive either reading.
>   * Deciding test, not run: the engine ++ J = 0 derivative-piece correction at m = 3–10Λ (F_NORM §2.6, §6).

> **Update (LEMMAS_CLOSE.md, unreviewed).**
> * **Lemma C′: → [D]**, with a closed-form check [DV] at ≤ 2·10⁻⁹.
>   * On Ω = {Re t₁, Re t₃ > 0} ∖ (0, ½] the swapped integrand is one single-valued analytic function, by an exact identity keeping z, z̄ off [1, ∞) [P].
>   * Lemma A is [D given C′]. The Einstein instance is [D] (not numerically evaluated), which answers M7.
> * **K1′: → [D]**, uniform on the Regge box at J₀ = 1, with no edge dependence.
> * **K2ν: → [D].** `|𝒦_ν| ≤ K̄` for all real ν.
> * **LRC leading form:** [D + DV] for scalar k = 2, 4 and EE k = 2; [D-sketch] for EE at k = 4.
> * **Correction to the H_K box above: the LRC uniform remainder is [I].**
>   * Over the large-η part of the Regge box the corrections are O(ρ = re^β), not O(r). Prop K_c's "extra factor r" step therefore fails uniformly, and with K1′ alone the remainder again needs J₀ < ½.
>   * The same gap sits in H7 §14's corollary for Ĉ′_EE.
>   * It closes given either of two inputs, both [I]:
>     * LRC\*: the de-shift cancellation survives the O(ρ) corrections;
>     * LC\*: uniform light-cone suppression of dDisc with exponent κ > ¼.
>   * **H_K stays [D-sketch].**
> * **a = b integrability:**
>   * the kernel bound is [D] (integrable iff growth p < 7/2);
>   * J₀ < 4 terms carry nothing at a = b [D];
>   * the borderline R⁴ arc in that corner is [I], tied to LRC\*/LC\*.

> **Update (LRC_CLOSE.md, unreviewed): LRC\* holds. H_K no longer needs J₀ < ½ or LC\*; it has no [I] link left.**
> * **The exact kernel.** The exact per-ν kernel is a d = 3 continuous-spin partial wave (dimension 1, spin −½+iν) in the (13)(42) channel, checked against direct quadrature to 3·10⁻⁸.
> * **The light-cone-edge corrections:**
>   * they resum to `₂F₁(½, h; h+½; ρ²)`, with no O(ρ¹) term;
>   * the result is analytic for Im ν < 3/2, so the de-shift cancellation survives to all orders [D + DV, scalar-type kernels].
> * **Criteria under the known J₀ ≤ 1:**
>   * Ĉ′_EE (k = 2): J₀ < 3/2;
>   * Ĉ₄^{ds}, ∂Ĉ₄^{ds}: **J₀ < 2**, which corrects HK_CLOSE's 5/2 (that needs f_pos(±3i/2) = 0);
>   * Ĉ₅^{ds}, Ĉ₃: J₀ < 3/2.
> * **H_K: [D-sketch]**, explicitly [D given K1′, K2ν, K3, DS, LRC-a, LRC\*]. Its weakest links are:
>   * the EE step (spinning partial wave → scalar ones with shifted labels, standard CPPR/weight-shifting, not checked end to end) [D-sketch];
>   * LRC-a for EE at k = 4 [D-sketch];
>   * DS for Q [NE].
> * **H7 §14 corollary (Result D):** it now holds for the exact kernel, with margin ½ [D-sketch through the EE step].
> * **LC\* is not needed.** The fixed-ρ light-cone expansion (Hartman–Jain–Kundu) is itself [A] for spinning operators.

> **Update (DESHIFT_CLOSE.md, unreviewed): M6 addressed. No bound value or binding margin changes.**
> * **C5 j′ = 1 coefficients.** The "[I]" label was stale: Stage 5 computed them [P + DV].
>   * The upper end 72.04 / 92.81 and the leading-Regge per-state values are exactly independent of the C5 weight [D].
>   * Beyond leading Regge, moving the binding state would need corrections of ≳ 3·10⁴ c₂ [D + M].
> * **Explicit continuum de-shift weights for C4 and C5, valid for all t ≥ 0 [DV].** The identity error is bounded uniformly by ε_W:
>   * C4 value: 1.2e-13;
>   * C4 derivative: 4.2e-11;
>   * C5: 2.4e-7.
>   * By comparison, the stored weights fail at t ≈ 41 (C4) and ≈ 20 (C5).
>   * The new weights decay super-Gaussianly, so they are already cut-off weights in HK_CLOSE's sense.
>   * Ω^Q_up can now be *defined* by these published weights, so DS is replaced by ε-DS [DV]. Exact C4 existence is [NE + D-sketch] and no longer needed.
> * **"C5 poor beyond t ≈ 15": fixed.** Worst-sign leading-Regge scan to β = 45: the Λ = 100 dip goes from +0.9645 to ≥ +0.9866; nothing is negative.
> * **Moments.** There is no real disagreement: RESULT_Q_STATEMENT mislabels V4FINAL's EE3 moments as "sc3". For the theorem's XP weights, |m_n − δ| ≤ 2.5e-16 (n ≤ 9).
> * **Independent Abel inversion:**
>   * D (NF02, k = 2) reproduced [DV];
>   * Q's C4 reproduced to 1e-16 by a different method;
>   * Q's stored C5 is not reproducible, but it is superseded by the new weight.
> * **Data published** in prl/data/ (weights, kernels, README, checksums).

> **Update (PCYC_CLOSE.md, unreviewed): P_cyc, as H_RS needs it, is an identity. The physical integration chain IS D′. H_RS [D given P_cyc] → [D].**
> * **Definition [P] + [DV].** CLS-D (arXiv:2311.04271) define the functional as an absolutely convergent integral over real Lorentzian points in the causal order 4>5>3>0>1>2. In the gauge x₃=−x₁=y, x₄=−x₂=x/(−x²), x₅=e this region is exactly D′, the real section cut at the light cones of x₅ (the zeros of B). Checked with an independent Minkowski code: 0 mismatches in 2·10⁵ configurations.
> * **Integrand [D] + [DV].** On D′ all 25 factor bases are positive, A, B > 0 and Q is real. So the physical integrand = point-independent constant × the project's integrand, and the constant cancels in ρ.
> * **End to end, scalar [DV].** ∫_{D′}, exact in Δ with the cut, against the exact CHMRS Mellin skeleton: |ε| ≤ 2·10⁻⁸ (Δ = 20), ≤ 6·10⁻⁹ (Δ = 30) for ν ≤ 5, and ≤ 1·10⁻⁹ to ν = 10 at Δ = 50. At these Δ, 82 % / 47 % / 3 % of ∫|Φ| lies beyond raw radius 0.3, outside the convergence region of the 1/Δ series. So this checks the chain itself, non-perturbatively.
> * **H_RS [D]** = RD_CLOSE Prop. D′ (engine = ∫_{D′} + O(Δ^N e^{−0.336Δ}); one [D-sketch] sub-step: the convergence radius of the x₀-dependent slow factors) + the identification above. REPRO_ENGINE's direct ∫_{D′} = engine K = 5 to 3·10⁻⁷ [DV].
> * **The crystal (saddle-point) cycle is only a numerical method; its class is no longer a premise.** Crystal = D: [P] in the Regge-scaled system; [I] at finite Δ (gap: no Stokes jump for a ∈ (0, 1/Δ], n± = 1 through the fold). A direct EE test T = ⟨Q⟩_crys/⟨Q⟩_{D′} − 1 finds no discrepancy (|T|/σ ≤ 1.7), but it is weak: (|T|+2σ)/threshold ≥ 4.8 above the fold and 16–2000 below it [NE, weak].
> * **Retraction.** RD_CLOSE's census statement "no critical point at raw radius 0.1–39" is withdrawn. There are 0–276 critical points per run at 0.3 ≤ r < 1.5, many with |B| ≈ 1, at heights 30–94 above the crystal (12 × 4000 starts, not saturated) [NE]. This does not touch D′.

> **Update (NEARGAP2.md, unreviewed): under A1′, ε* = 0.038034 leaves zero flat margin on a massless scalar; use ε = 0.040.**
> * With ε*, the repaired functional Ω^Q_up + ε*Ω^EE_θ on a light J = 0 has flat margin ≈ 3.3 (m/Λ)², which vanishes as m → 0. For Δ_L = O(1) the sign is set by O(Λ⁻²) AdS terms that nothing controls, so it is **[I]** (Lemma C can sign it only for Δ_L ≳ 30 at 1000 and ≳ 15 at 300).
> * **Adopted fix [D-sketch]: ε = 0.040.** The flat margin at m → 0 becomes 0.022 N, ×100 above every computed curvature effect. The upper end under a pure higher-spin gap becomes **164.51 (Δ_g = 1000) / 125.71 (300)**, up from 160.99 / 123.08; the c₁² coefficient has to be recomputed for ε = 0.040.
> * The published 72.04 / 92.81 are unaffected when no scalar lighter than 0.46Λ couples to TT.

> **Update (C1_VALUE.md, unreviewed): c₁ = 0 given one named premise F₂ [D-sketch given F₂; F₂ is [A]]. The shockwave route gives exactly zero and cannot see c₁ in principle. There is still no unconditional value or rigorous bound [I].**
> * **Shock route [P].** In d = 3 the preserved-helicity elements of the Costa–Hansen–Penedones impact-parameter operators (1707.07689) are ¼□(□+ν²+¼) (t₄-type) and ½(□+ν²+¼) (t₂-type). Both vanish identically on every H₂ eigenfunction, for every ν and every L, so the "polynomial + pole" remainder is zero. A ν-independent term is a δ-function at zero impact parameter, which shock computations at L > 0 cannot fix [D]. ADS4_LIGHT's t₄[V(ν)−V(−3i/2)]/(ν²+9/4) is identically zero in d = 3, so c₁ does not originate there.
> * **Structure [D] (one [D-sketch] input: a regular second-order anti-self-dual solution).** In the preserved⊗preserved configuration (legs 1, 2 self-dual, 3, 4 anti-self-dual), the whole O(ĝ₃) tree (s- and t-channel exchanges plus the quartic seagull) collapses at integrand level to 3ĝ₃∫TrΨ₁Ψ₂Φ₃₄. This is a factorised energy-channel T exchange, ⟨T⁺T⁺T⁺⟩_{W₊} × [⟨T⁺T⁻T⁻⟩_E − κ·CS contact]/⟨TT⟩. So c₁ is fixed by ⟨TT⟩ and ⟨TTT⟩ data alone (hence universal), and there is no independent seagull. Checked against Jain–John–Mehta–Nizami–Suresh (2106.00016) [Q + DV]. This corrects C1_COMPUTE §3.5 and §7.1 (third bullet).
> * **Premise F₂ [A].** The k = 2 functional sees the s-channel T family only through light-ray one-point matrices: diagonal for Einstein, purely helicity-flipping at O(t₄). This is the k = 2 analogue of H7's factorisation F, which is [D] for k = 3. Given F₂, the O(t₄) term flips helicity on the (2,4) line and EE annihilates it [P], so c₁ = 0. Not written out: F₂ for the k = 2 weight f = v′ − u′, and the map from Euclidean momentum-space helicity to Lorentzian light-ray helicity for the intermediate T.
> * **Without F₂.** D moves by < 1 % iff |c₁| < 4.61 (300) / 59.6 (1000). Q moves by < 1 % iff |c₁| < 103 (300) / 1306 (1000). Letter wording: "uncomputed universal constant, fixed by ⟨TTT⟩ data, zero given F₂", not "free parameter". No numbers change.

> **Update (R4ARC_HEFT.md, unreviewed): L6′ [I] → [D] with a free-standing proof; H_EFT(a), spinning [D-sketch] → [D]; H_EFT(b) stays [A], with no assumption beyond CHMRS. Result Q's upper end goes from [I] in one link to [D-sketch].**
> * **L6′.** In the og/CLS-D gauge the collinear variable y = 2sin²(Δφ/2) enters only through X₁₃ and an explicit y⁻¹ in f₂, f₄. Counting Lemma Y [P] (exhaustive over 997,836 parity-even and one-ε structures, worst margin 0): every term of a spin-2 contact diagram is O(y⁰ log y). So |L_{R⁴}| ≤ C Δφ⁻²(1+|log Δφ|): p = 2 against the criterion p < 7/2, margin 3/2. It uses neither H_EFT(a) nor LRC\*, and the same bound holds for all contacts and heavy-exchange double-trace parts.
> * **Premise: dispersiveness** (finitely many non-zero s-channel residues δ₁₂ = −n). [D] for integer Mellin shifts, and CLS-D's statement [Q]; every holographic ⟨TTTT⟩ tree term qualifies. In non-dispersive toy classes (half-integer free-field links, 1243) the collinear growth reaches the edge (slope 4.7 → 3.5), so LEMMAS's worry was real for non-holographic terms. LEMMAS §4.1 item 4 and LRC §6.5 item 3 (η^{5/2}, p = 5/2) are superseded.
> * **H_EFT(a) [D] for dispersive terms.** It holds pointwise in angle, hence also ν-paired: only finitely many double-trace residues at bounded Mellin arguments enter. "False pointwise in angle" (F_NORM B3, §3.4, §4) applies only to non-dispersive toy classes. Checks on T′ = σ⊗T: σ-exchanges with J_h = 0, 2, 4 and τ_h = 41–81 at k = 2, 4; the term-wise expansion matches the double-trace value to 10⁻¹⁸–10⁻³⁷. An unexplained ≈ −3 % leading-order ratio for J_h = 2 stays at [NE]; route II needs only the scalar case.
> * **Residual:** neither the actual AdS₄ R⁴ ⟨TTTT⟩ nor a spinning χ-exchange was evaluated. The physical Δφ scan is not yet asymptotic at Δφ = 0.04 (slopes 2.5–2.8), so p = 2 rests on the analytic checks. Refs to add to the bib: Paulos 1107.1504, Costa–Gonçalves–Penedones 1404.5625, Penedones 1011.1485.

> **Update (C1_FINAL.md, unreviewed): c₁ = 0 EXACTLY [D] (same grade as the implication it feeds; was [D-sketch given F₂ [A]]). The headline numbers are unchanged, and the c₁ caveat and every c₁ error term are removed.**
> * **Route A, F₂ proved [D + DV]:**
>   1. With light rays based at x₀ = ∞, the k = 2 weight f = v′−u′ is affine [P], so the functional contains only E₃E₁, E₃E₁⁽¹⁾ and E₃⁽¹⁾E₁.
>   2. The first moment E⁽¹⁾ = ½{E, −i n̄·D} has the same helicity pattern as E on the T family, for every parity-even ⟨TTT⟩ (free boson and free fermion span the structures). Ward identities exclude a diagonal O(t₄) or O(α_T) element.
>   3. The EE (2,4) contraction is helicity-diagonal through first order [P].
>
>   So the O(t₄) T-block term vanishes pointwise for every ν, after the de-shift. The Euclidean-vs-Lorentzian helicity point is resolved [D-sketch, not load-bearing].
> * **Route B (direct):** the k = 2 EE action on the explicit energy-channel exchange is zero to ≤ 1.75·10⁻¹¹ at 8 points × 3 structures; the controls (flip-sensitive contraction, Einstein term, W×W term) are O(1).
> * **Consequences:** Result D is √(t₄²+α_T²) ≤ 2.89126·10⁻⁸ (1000) / 3.13206·10⁻⁶ (300), with a symmetric t₄ interval. Result Q's upper end is 92.8135 / 72.0406 with no c₁² term. Referee M5 / RU issue 14 are closed. The whole O(t₄) light term of Ĉ′_EE vanishes, not just the constant, and an O(α_T) constant vanishes too.
> * **Inherited links:** Ω_OPE = Lemma R on tree terms [D]; the [TT] double zeros at the borderline n = k; and the choice of the EE dual structure (Ω^EE and Ω^Q_up both use EE at k = 2). An independent Lemma-R evaluation of the Witten diagrams and the 2106.00016 cross-check were not done.

> **Update (QBIND_COROLLARY.md, unreviewed):**
> * **The quartic functional Ω^Q_up's binding states are error-controlled at Δ_g = 1000.**
>   * The binding state ++ J=0 (m = Λ): +0.14282, lower bound +0.1405, **L = 62 [NE±]** (was [NE+M] given H_RS). Its five direct nodes agree within 1.1σ.
>   * OO J=0: L = 58.
>   * ++/OO J=2 channel minimum: L = 114.
>   * +− J=5 channel minimum: L = 90.
>   * +− gap states J = 4, 5, 6: [NE±] given K₂ edge rows (L ≥ 71). With strict charging of x > 0.955 they fall to L = 0.4 / 1.2 / 1.1, so **+− J=4 at the gap is [NE+M] under the strict protocol** (a channel minimum at 0.259, not the binding state).
>   * The largest budget term is the C5 convention (10⁴ c₂ vs measured 10–330 c₂).
>   * Nothing negative.
> * **Corollary (energy-flux isotropy) [D from the cubic bound; CDP 1707.03007 kinematics [Q], angular algebra [DV]]:** ⟨𝓔(θ)⟩ = (E/2π)[1 + A cos(4θ−φ)] with **A ≤ τ/4 ≤ 7.23·10⁻⁹** at Δ_g = 1000 (ANEC allows A ≤ 1, saturated by free fields). Equivalently |n_s − n_f| ≤ 7.23·10⁻⁹ c_T and |p_T| ≤ 1.11·10⁻¹⁰ c_T. No κ or light-scalar condition is needed. Literature note: Buchel et al. 0911.4257 eq. (3.8) prints the d=3 constant as 1/8; it should be 1/4.

## 0. Headline

**Theorem-form statement (conditional).** Assume A1 and A2 at `Δ_g ∈ {300, 1000}`, the analytic inputs of §4, and per-state positivity of the upper
functional `Ω^Q_up` and the lower functional `Ω^Q_low` on every twist-allowed heavy block. Then the CFT-defined R⁴ coefficient `g̃₄`
(§1; equal to the bulk `g₄` given (F)+(P)) satisfies

    −O(Λ⁻²) ≤ g̃₄ Λ⁶/(8πG) ≤ 2L_AdS − κ|ĝ₃|²Λ⁸ + ε_{c₁},
    2L_AdS = 17.241·log(2Λ/3) + B_eff,   B_eff = −19.220 (Δ_g = 300) / −19.265 (Δ_g = 1000),
    i.e.  2L_AdS = 72.04 (300) / 92.81 (1000).

* **Cubic terms.**
  * **Case `t₄ = α_T = 0`** (CFT-checkable; then `ĝ₃ = 0`): the `κ` term drops, and so does `ε_{c₁}` (∝ `t₄`).
  * **Case general `ĝ₃`:** the upper bound needs **`κ ≥ 0`** [D-sketch + M], with `κ = 0.110`. Then `ĝ₃ ≠ 0` only tightens it.
  * If `κ` could be negative, Result D's `|ĝ₃|²Λ⁸ ≤ 1.07·10⁴` would allow a shift of order `10³`, about 15× `2L_AdS` (RU issue 7).
  * The lower bound for `ĝ₃ ≠ 0` reads `g̃₄Λ⁶/8πG ≥ |ĝ₃|²Λ⁸/(4c_L) − O(Λ⁻²)` with `c_L > 0`. It uses the sign of the W³W³ light term
    of the lower functional (§2.4).
* **`c₁`** (free) is the coefficient of the `E × W³` light term (`k = 2` and `k = 3` parts).
  * `|ε_{c₁}| = 2|c₁t₄|Λ²∫ψ₂`, with `∫ψ₂ = 1.016` for the δ = 0.2 weights.
  * With Result D's `t₄` bound: `|ε_{c₁}| ≤ 0.57|c₁|` (0.79 % of `2L_AdS`, 300) and `≤ 0.059|c₁|` (0.06 %, 1000).
  * Importing D's bound makes Q inherit D's `H_LJ`, `H_M`, `H_R`, `H_edge` and `H_RS`. Carrying `c₁t₄` explicitly avoids that.
* **Heavy side.**
  * Full-set minimum of FF/N for the δ = 0.2 functional: **`+0.1426 ± 0.0017`** (SB gauge `+0.1409`).
  * Binding pair: `++` J = 0 and OO J = 0 at the gap, at both `Δ_g` [NE + M given `H_RS`].
  * The binding margin in the derivative piece is `s_fail/s_d = 1.78/0.51 ≈ 3.5`.
* **Not claimed.**
  * `Δ_g = 100` or any other `Δ_g`.
  * Loop level.
  * The δ = 0.05 functional as a headline (its margin is `+0.0060`, 10 % of δ, and depends on `σ_hU`).
  * **Any cubic closure `|ĝ₃|²Λ⁸ ≤ …` from Q** (§5.3).
* **ABJM:** `g₄Λ⁶/8πG = 2ζ(3)` lies at 3.3 % (300) / 2.6 % (1000) of the bound. It violates A1, so this is a consistency check only (§6).

---

## 1. Setting, normalisations, functionals

* **A1, A2, `c_T`, `G`, `Λ`, `m`** as in `RESULT_D_STATEMENT.md` §1:
  * free real scalar `c_T = 1`;
  * `8πG = 256/c_T`;
  * `Λ = √((Δ_g−2)(Δ_g−1)) = 298.50 / 998.50`;
  * `m² = (Δ−J−2)(Δ+J−1)`;
  * gap state `m_gap² = (Δ_g−2)(Δ_g−1+2J)`.

  These conventions are identical in `resultq/stage6_model.py` and Result D (RU §1).
* **Couplings.** `g₄` is the flat R⁴ Wilson coefficient in Stage 0's normalisation: the type II Virasoro–Shapiro amplitude has
  `g₄/(8πG) = 2ζ(3)/M⁶`, `M² = 4/α′`. `ĝ₃` is the CHLPS `R³` normalisation, with `√(t₄²+α_T²) = 240|ĝ₃|` (AdS units).
* **`g̃₄` (CFT-defined).** After Lemma A removes every tree term with Regge intercept `J₀ < 4`, `g̃₄` is defined by the light side of the
  de-shifted `k = 4` value functional, `−Ĉ₄^{ds}|_light = 2g̃₄`. It equals the bulk `g₄` given (F)+(P) [D-sketch].
* **Upper functional `Ω^Q_up`.**
  * The flat functional is δ = 0.2, `Q3-full-delta0.2` of `code/resultq/functionals_q3.json`: `A = 2ψ₂′(0) = 17.241`, `κ = 0.1099`.
  * It is realised in AdS by the de-shifted `Ĉ₄^{ds}`, `∂Ĉ₄^{ds}`, `Ĉ₅^{ds}`, with the C4/C5 weights of Stage 7
    (`functionals_q7_xprec.json`), plus the smeared `k = 2, 3` parts `∫ψ₂Ĉ′_EE`, `∫ψ₃Ĉ₃` (Result D's objects with Q's weights).
  * `s_d` and FF/N are quoted in gauge H0, with H2/H3/H4/SB/XP as variants (spread `≤ 0.016`).
* **Lower functional `Ω^Q_low`.** The value piece: heavy image `P_{−½}(cosh β) > 0` at leading Regge, flat
  `2(|c⁺⁺|²+|c⁺⁻|²)/m⁸ ≥ 0`, light side `2g̃₄` (RU issue 17; `RESULT_Q_REVISION` §3.3).
* **Fixed as of 22 Nov 2026.**

## 2. Light side

### 2.1 Lemma A + D_LR (tree terms with `J₀ < 4` contribute zero)
* **Pointwise Lemma A [D]** (V4F §2.2). For every tree term with `J₀ < 4` (Einstein, `E × W³`, MFT), the integrand integrated over
  `t₁, t₃` and symmetrised under 1↔3 vanishes identically at `φ₁ ≠ φ₃`.
  * Inputs: microcausality, CLS-D endpoint convergence (`J₀ < k`), Lemmas B and C, `J₀(E×W³) ≤ 3`.
  * So this light side is the zero function in every order of the ν-, angle- and `t`-integrations.
  * The earlier ν-inside/outside circularity is resolved.
* **D_LR [D] (consequence [P]).** CLS-D's domain `x₁₃² > 0`, with the operators on their null contours, fixes the cone-preserving
  prescription. On the light-transform manifold, `W_u·f₄ ≡ 0`. The timelike-tube regularisation is not a regularisation of CLS-D's
  integral.
* **Checks at `O(1/C_T)` [DV]** (OG, FH).
  * **Correlators.** `Ψ₄` and `Ψ₅` were evaluated in swapped form on the exact connected `⟨TTTT⟩` of:
    * the free conformally coupled scalar: all 25 (2,4) components, `k = 2–5`, residual `≤ 8·10⁻¹⁵`; per ν (ν = 1, 3) `≤ 5·10⁻¹³`;
    * the free Dirac/Majorana fermion: `≤ 3.8·10⁻¹⁵`; per ν `≤ 1.4·10⁻¹⁵`.

    Hence they vanish on every boson + fermion mixture, including the `t₄ = 0` mixture that has Einstein's `⟨TTT⟩`.
  * **Generic exponents and logs.** GFF `[φφ]₀,₂` and `[ψ̄ψ]₀,₂` families and their Δ-derivatives: `≤ 3·10⁻¹¹` / `1·10⁻¹⁰`.
  * **Controls fail at O(1).**
    * Scalar: non-primary `T` 0.5–1.7; wrong light-transform weight 0.4–1.8; no T-shift 0.03–0.04; wrong symmetrisation 2.00.
    * Fermion: weight, primarity and frame controls 0.04–1.8. The no-T-shift control is degenerate there; the branch rule is pinned by
      the 2- and 3-point validation instead (`κ = i` gives `−1`).
  * **Corner.** The collinear-corner law `Δφ^{−1/2}` (scalar) confirms V4DPRIME's twist-1 homogeneity lemma [DV].
  * **Limitation.** The fermion is blind to the (13) link and the collinear corner in the CLS-D gauge (two Wick classes vanish
    identically [P]). The (13)-channel/corner part of D_LR is therefore checked on the scalar only.
* **Grade.** Lemma A + D_LR [D]; mechanism [DV] at `O(1/C_T)` on two free theories. **The Einstein instance is unchecked [I].** No
  light-ray-usable Einstein `⟨TTTT⟩` exists in the repository. Its leverage on `g₄` is `10¹⁰–10²²` (RU issue 9).

### 2.2 Contacts and heavy-exchange double traces (`J₀ ≥ 4`): (F)+(P) and `H_EFT`
* **(F)+(P) [D-sketch].** The spinning flat limit and polynomial light side: `−Ĉ₄^{ds}|_{R⁴} = 2g₄`, and the `g₅` cancellation. The R⁴
  normalisation control in position space is not reachable: no `J₀ = 4` spinning correlator exists (OG §5.1, FH §6).
* **`H_EFT` splits in two (FH §4.5, refining RU issue 9 and OG §5.2).**
  * **(a)** The term-wise (P) value of a heavy exchange's double-trace part equals its actual (OPE-side) double-trace value.
    * **[DV] in the scalar-Mellin toy** (`J_h = 0, 2, 4, 6`; `Δ_h = 20, 50, 100`). The error at the asymptotic optimum or numerical
      floor is `≤ 1·10⁻⁹` (`Δ_h = 20`) to `≤ 1.5·10⁻⁶` (`Δ_h = 100`) of the heavy image, `≥ 10⁵×` below the bound's sensitivity.
    * Spinning and position-space versions: [D-sketch].
  * **(b)** The heavy-sum / contact-tower interchange ("no conspiracy"): **[A]**, untested. The absolute moments grow super-geometrically,
    so the term-wise series is at best asymptotic.
  * Caution: "double traces = −heavy block" holds only for `J_h ≤ 3`, where it is Lemma A, checked `D = −H` to the numerical floor. For
    `J_h ≥ 4`, `Arc/H ≈ −30 … −35` (`J_h = 4`) and `≈ −2.5·10³ … −3.4·10³` (`J_h = 6`), tending to the flat values `−35`, `−3423`.
    Result Q does not use the `J_h ≤ 3` shortcut.

### 2.3 The `|ĝ₃|²` sector and `c₁`
* **W³W³.** The u-channel diagram has Regge spin 2 and is annihilated by `k = 4`. The s- and t-channel diagrams act like contacts:
  `−Λ²∂_{ν²}Ĉ₄(0)|_{W³W³} = 4πG|ĝ₃|²` [DV algebra in the scalar-Mellin form; spinning D-sketch] (`RESULT_Q_REVISION` §0.3).
* **`κ = 0.110`.** Flat moments of the δ = 0.2 functional, plus H17-type `O(Λ⁻²)` terms. **`κ ≥ 0`: [D-sketch + M]** (RU §4). The earlier
  table grade [I] is superseded.
* **`c₁`:** free (§0). It includes the `Ĉ₃` (`k = 3`) `E × W³` term as well as the `k = 2` one (RU issue 14).

### 2.4 Lower bound
* At `ĝ₃ = 0`: `g̃₄Λ⁶/8πG ≥ −O(Λ⁻²)` **[D-sketch]**, limited by (F)+(P) and `H_K` (RU issue 17). It uses only the value piece. It is
  unaffected by `σ_hU` and by `s_d`. "≥ 0" really means "≥ −O(Λ⁻²) × (EFT data)", from `b₀`, `s_v` and the contact tower.
* For `ĝ₃ ≠ 0`: `≥ |ĝ₃|²Λ⁸/(4c_L)`, with `c_L = 0.509` (measured `s_d`; slope `0.49`) or `0.98` (`s_d` double counted; slope `0.255`).
  * The slope is positive in both readings [D-sketch + M]: W³W³ sign as in §2.3, plus HYB.
  * This item is not in RU §5's wording. It is added here because the lower functional's `|ĝ₃|²` term enters once `ĝ₃ ≠ 0`.

---

## 3. Heavy side: per-state positivity

### 3.1 Which channels can bind (FWD §1, tolerances `s_fail` = the `s_d` at which a channel's FF reaches 0)

| channel | `s_fail` (δ = 0.2) | `s_fail` (δ = 0.05) | measured `s_d` (gap; 300 / 1000) | sign | source |
|---|---|---|---|---|---|
| **`++` J = 0** | **+1.78** | +0.58 | `+0.509 ± 0.016` / `+0.508 ± 0.002` | unfavourable | BF2 |
| **OO J = 0** | +1.78 | +0.58 | `+0.509 ± 0.016` / `+0.510 ± 0.002` | unfavourable (= `++`) | FWD §2 |
| `++` / OO J = 2 | +2.81 | +0.83 | `++` `+0.511 ± 0.010` / `+0.508 ± 0.001`; OO `+0.511 ± 0.010` / `+0.510 ± 0.001` (1.05Λ) | unfavourable | BF2; FWD §2 |
| `+−` J = 4 | +2.38 / +2.10 | +1.04 / +0.75 | `−0.516 ± 0.039` / `−0.510 ± 0.004` | favourable | BF3 |
| `+−` J = 5 (parity-odd) | +5.72 | +1.78 | `−0.539 ± 0.033` / `−0.512 ± 0.004` | favourable | PO §2 |
| `+−` J = 6 | +9.5 / +9.4 | +1.99 | `−0.564 ± 0.012` / `−0.515 ± 0.003` | favourable | FWD §3 |
| `+−` J = 7, 8, 9, 11 | ≥ 20 | ≥ 8.9 | 1000 (used): `−0.48 … −0.52`; 300: `−0.43 … −0.59` | favourable | FWD §3 |
| `++` J = 6, 8 | +36 / +66 | +28 / +52 | `+0.511 ± 0.008` / `+0.54 ± 0.41` | unfavourable | FWD §3 |

* The correction is helicity-odd and universal: `Γ_K^h ≈ σ_h U(ν)`, with `σ_{++} = σ_OO = +1`, `σ_{+−} = −1`, and `|s_d| ≈ 0.51 ± 0.08`.
* Unmeasured channels (`++` J ≥ 10, `+−` J = 10 and ≥ 12, OO J ≥ 4, masses `> 1.475Λ`) have `s_fail ≥ 11`, or are `(Λ/m)⁶`-suppressed.
* **Full-set minimum.**
  * δ = 0.2: `+0.1426 ± 0.0017` (SB `+0.1409`). The 300 state alone is `+0.1427 ± 0.0017`; the 1000 state `+0.1426 ± 0.0002`.
  * δ = 0.05: `+0.0060 ± 0.0013` (SB `+0.0047`).
  * Both `Δ_g`, all gauges [NE + M] (FWD §0; unchanged by HRS).

### 3.2 `H_RS` for the binding state (HRS; refines RU issue 3)

* **Weight.** `s_d = 0.04 Σ_ν χ_d(ν) Re Γ_K(ν)`.
  * Per window, engine `s_d` (H0) is `+0.10` [0,1), `+0.57` [1,2), `−0.05` [2,3), `−0.10` [3,4), `−0.01` [4,8), `≈ 0` above.
  * `ν < 2` carries 68 % of the weight, and `ν < 4` 89 %.
* **Crystal-cycle test** (direct n3/n4 nodes, `J = 0` at the gap at 300, 1000, `Δ = 100`; `J = 2` at 1.05Λ). Crystal = engine beyond
  leading order:
  * **[NE] at `ν ≥ 10`**: `Γ_exact` up to `−209` per harmonic reproduced to 0.3–1 %; e.g. `ν = 20.02`: `−0.3 ± 1.5` against
    `Γ_exact = −73.4`;
  * **[NE, weak] at `6 ≤ ν < 10`**;
  * **[I] at `2 ≤ ν < 6`**: consistent, but the errors (`±39 … ±1800` per harmonic, with ×3 seed inflation) are 3–30× the δ = 0.2
    thresholds;
  * **untested at `ν < 2`**.
* **Propagated.** On `ν ∈ [5, 25]` at 300: `Δs_d = −0.10 ± 0.35` and `ΔFF/N = +0.011 ± 0.039` (×3 errors), against the margin
  `s_fail − s_d = 1.27`. On `[2, 25]`: no information.
* **To break δ = 0.2.** A uniform per-harmonic error of 30 (ν ∈ [0,1)) or 17 ([1,2)) would break it.
* **`H_RS` for Q's binding state stays [I].** `s_d = +0.51` is [NE given `H_RS`].

### 3.3 Other heavy-side items

| item | content | grade | source |
|---|---|---|---|
| `H_OO^Q` | OO J = 0, 2 measured, equal to `++` (helicity lemma at `k = 4`: order 0 to `10⁻⁹`; beyond, `≤ 1.3·10⁻³ c₂` per harmonic; `|Δs_d| ≤ 5·10⁻⁶`) | [NE + M], i.e. **[D given `H_RS`]**: OO and `++` use the same engine (RU issue 11) | FWD §2; PO §1 |
| `H_fwd^Q` | unmeasured forward channels; the `Δ_g = 300` J-drift of `s_d` (`−0.07` at J = 7, 2.4× the truncation error, cause [I]), J ≥ 7 quoted from 1000 (a substitution; tolerance ×17 at δ = 0.05) | [NE + M] given `H_RS`; drift [I], harmless at δ = 0.2 | FWD §3; RU issue 12 |
| `H_smear^Q` | graviton corrections to the smeared `k = 2, 3` parts at the binding states: `≤ 3·10⁻⁵` of B₂ against a tolerance `≈ 24 %` | [NE + M]; edge `x > 0.8` and `k = 3` below the fold bounded by tolerances | REVISION §1 |
| `H_LJ^Q` | large-J near-gap corrections stay below Q's tolerance: `≥ 1.1·10⁴ c₂` (`+−` J = 8, 300), `≥ 1.6·10⁵` (J ≥ 16), against measured `10²–10⁴`. The forward term `a(2J(J+1) − c₄)u³` dominates. | [M], margin ×10–10³ | RU issue 13 |
| HYB | scalar skeleton + leading Regge accuracy away from measured states; `β > 30` | [M] (was [NE]) | RU §4 |

---

## 4. Hypotheses (single list)

| # | hypothesis | role | grade | evidence | what would upgrade it |
|---|---|---|---|---|---|
| 1 | **A1**, **A2** | setting | [A] | — | — |
| 2 | **K1** (KKSDZ 1904.05905, with its loophole), **K1′**, K2/K3 for the de-shifted functionals (validity `J₀ < 7/2`) | Regge | [Q] + [D-sketch] | V4F §5 | closing the loophole; writing out K1′ |
| 3 | **`H_K`**: heavy OPE sum ↔ (corner, ν) interchange for the heavy-normalised `Ĉ₄^{ds}`, `∂Ĉ₄^{ds}`, `Ĉ₅^{ds}` (same statement as Result D's) | the only load-bearing Regge input | **[D-sketch]** ([D given `H_ip`], raw corner `J₀ < ½`) | V4F §2.4 | Proposition K for the heavy-normalised kernel |
| 4 | **Lemma A + D_LR** | `J₀ < 4` light side `= 0` | **[D]**; mechanism [DV] at `O(1/C_T)` (free scalar + free fermion); **Einstein instance [I], unchecked** | V4F §2.2; OG; FH | one swapped-form `Ψ₄` on Einstein `⟨TTTT⟩` (needs O(G) `[TT]` data or a light-ray-usable representation) |
| 5 | **(F)+(P)** spinning | EFT light values (`2g₄`, `g₅` cancellation) | [D-sketch] | V4F §5 | a `J₀ = 4` spinning position-space correlator (AdS contact diagram) |
| 6a | **`H_EFT`(a)**: term-wise (P) of heavy-exchange double traces = actual value | `J₀ ≥ 4` light side | [DV] scalar-Mellin toy; [D-sketch] spinning | FH §4 | spinning toy |
| 6b | **`H_EFT`(b)**: no conspiracy in the heavy sum / contact tower | same | **[A]** | — | — |
| 7 | **De-shift** existence: continuum C4/C5 solutions for all β | heavy images beyond `t = 30` | [NE] to `t = 30`; C5 spinning coefficients [I]. The stored sc3 weights have growing high moments (`m₄ = −1.4·10⁻²`, `m₅ = −27`, ratio `≈ 1.9·10³`, a ν-grid cutoff effect), harmless for `Λ ≥ 300` (ratio `≤ 0.02` per order). | Stage 7; FH §4.4 | — |
| 8 | Adm, H7 (inherited from D) | `k = 2` smeared part | [D] | Result D §4 | — |
| 9 | **`κ ≥ 0`** | only if `ĝ₃ ≠ 0` (upper) | [D-sketch + M] | REVISION §3.3 | spinning W³W³ |
| 9′ | `c_L > 0` (lower slope) | only if `ĝ₃ ≠ 0` (lower) | [D-sketch + M] | REVISION §3.3 | same |
| 10 | **`c₁`** | `E × W³` light term | free parameter, carried | RU issue 14 | Regge residue of `E × R³` |
| 11 | **`H_RS`** | every `s_d`, incl. the binding state | **[I]** for `ν < 6`; [NE] `ν ≥ 10`; [NE, weak] 6–10 | HRS | an EE `k = 4` crystal value at `ν < 4` to ~10–20 per harmonic (lower Δ, variance control for the `R` factor, or an analytic continuation argument) |
| 12 | `H_OO^Q` | OO J = 0, 2 | [D given `H_RS`] | FWD §2 | — |
| 13 | `H_fwd^Q` | unmeasured channels, 300 J-drift | [NE + M]; drift [I] | FWD §3 | — |
| 14 | `H_smear^Q` | smeared parts | [NE + M] | REVISION §1 | — |
| 15 | `H_LJ^Q` | large-J near-gap | [M], margin ×10–10³ | RU issue 13 | — |
| 16 | HYB accuracy; `β > 30` | rest of the heavy action | [M] | RU §4 | — |
| — | `σ_hU` double counting | **not needed** for the δ = 0.2 headline. It decides δ = 0.05 and the (withdrawn) closure. | [I]; the flat-limit argument favours double counting | RU issue 8 | — |
| — | `H_νout`, H_LT, V₄/V₄′/V₄″/V₄‴ | **deleted/replaced** by item 4 | — | V4F | — |

---

## 5. Caveats

1. **The weakest links.**
   * `H_RS` at `ν < 6` for the binding state [I]. An engine error of ×3.5 in the derivative piece at J = 0 would break δ = 0.2, and
     nothing independent excludes it.
   * Then `H_EFT`(b) [A], the unchecked Einstein instance of Lemma A, `H_K` and (F)+(P) [D-sketch].
2. **`g̃₄` vs `g₄`.** Only `g̃₄` is CFT-defined. Its identification with the bulk coupling is (F)+(P).
3. **Cubic closure withdrawn** (RU issue 8).
   * The values "`|ĝ₃|²Λ⁸ ≤ 120–198 / 154–257`" are not claimed.
   * The `Δ_g`-independent lower slope 0.49 (with the measured `σ_hU` correction) is 1.9× CHLPS's flat optimum `0.26`, which an AdS
     bound whose sum rules reduce to the flat ones cannot exceed in the flat limit. This favours the double-counting reading.
   * At most a remark: `|ĝ₃|²Λ⁸ ≤ 2L_AdS/(κ + ¼) ≈ 200 / 257` [D-sketch + M]. No ratio `|ĝ₃|²/g₄` is claimed.
   * Result D's `1.07·10⁴ / 1.43·10⁴` stands under fewer hypotheses.
4. **δ = 0.05.** Its minimum `+0.0060 ± 0.0013` is positive but depends on `σ_hU`. It is not the headline.
5. **Rounding.** `B_eff(1000) = −19.2647` (REVIEW_RESULTQ §1) has been quoted as "−19.27" in earlier documents; to two decimals it is
   −19.26. The effect on `2L_AdS = 92.81` is `< 0.01`.

## 6. ABJM consistency check

* **Dictionary** (Stage 0; arithmetic re-checked by RU §1):
  * `M² = 4/α′`, `L²/α′ = π√(2λ)`, so `Λ ↔ ML = 2^{5/4}π^{1/2}λ^{1/4} ≈ 4.22 λ^{1/4}`. The identification `Λ ↔ ML` holds only at leading
    order in `1/(ML)`.
  * The type IIA `R⁴` term `ζ(3)α′³/32` gives `g₄M⁶/8πG = 2ζ(3) = 2.404`, independent of `λ`.
* **Placement.** `2ζ(3)/2L_AdS = 3.3 %` (300) and `2.6 %` (1000), recomputed: `0.0334`, `0.0259`.
* **`ĝ₃ = 0`** because 10d type IIA supersymmetry forbids an `R³` term. It does *not* follow from "maximal SUSY": ABJM at `k > 2` has
  `N = 6`.
* **`c_T` does not enter the placement.** Do not confuse flavour-current coefficients with `c_T`:
  * for `N = 2` theories, `∂²_m F` gives a flavour-current two-point coefficient, not `c_T`;
  * the identity `c_T = (64/π²)∂²_m F` (free real scalar `c_T = 1`) relies on the mass coupling to a current in the stress-tensor
    multiplet. It is used only in `ABJM_THEOREM`'s separate `c_T` result, consistently with `8πG = 256/c_T`.
* **Outside A1** (light KK modes); only the flat image (F) is used; log-dominated. **Consistency check, not evidence.**

## 7. Literature

* **CHMRS 2106.10274** (Caron-Huot–Mazáč–Rastelli–Simmons-Duffin) is the **methodological precedent**: two-sided `log(Δ_gap)` bounds
  on contact couplings in AdS₄ with dynamical gravity, from uplifted flat sum rules, for scalar correlators. Q's light-side formula and
  `H_EFT` are of CHMRS's (5.3) type.
* **CLS-D 2311.04271** (Chang–Landau–Simmons-Duffin; JHEP 04 (2025) 016) build the spinning light-ray sum rules and the crystal
  dictionary, and leave the uplift of flat bounds to AdS as future work. **Result Q is a first explicit-constant, conditional execution
  of that programme for `g₄`**. Method priority: CLS-D and CHMRS.
* **CHLPS 2201.06602:** the flat D = 4 two-sided graviton bounds (log IR). Stage 0 reproduces the flat ¼ from the same forward rules.
* **To be discussed (not only cited):** Kundu 2104.11238 and the `N = 6` bootstrap (Binder–Chester–Jerdee–Pufu 2011.05728).
* **Novelty [I] at search level.** No `⟨TTTT⟩` sum-rule bound on `R⁴` was found (abstract-level search plus one extended search, RU §1).

## 8. Open problems (value order)

1. `H_RS` at `ν < 4` for the binding state.
2. One swapped-form `Ψ₄` on Einstein `⟨TTTT⟩` (routes: O(G) `[TT]` data in every helicity structure, days; momentum space, impractical).
   Also the (13)-channel/corner part of D_LR on a correlator with a spin-1/spin-2 single-trace exchange.
3. `H_EFT`(b); spinning `H_EFT`(a).
4. `H_K`.
5. (F)+(P) in position space; the R⁴ normalisation control.
6. Spinning W³W³ (would firm up `κ ≥ 0` and `c_L > 0`); `σ_hU`.
7. The 300 J-drift of `s_d`.
8. `c₁`.

---

## Appendix A. History (superseded content, in order)

1. **Stages 0–7** built the forward functionals.
   * δ = 0.05 (`A = 11.45`) was replaced by δ = 0.2 (`A = 17.241`) as the quoted bound.
   * The `Δ_g = 100` column was dropped.
2. **REVIEW_RESULTQ: major revision.**
   * `κ`, the ¼ slope and `4πG|ĝ₃|²` were graded [I].
   * V₄ (swapped form, `[TT]₀` tail) was graded [D-sketch].
3. **REVISION.**
   * `H_smear^Q` → [NE + M].
   * The W³W³ residue in `∂Ĉ₄(0)` was found to be zero (Mellin).
   * Cubic terms were restored at [D-sketch + M], and a closure `|ĝ₃|²Λ⁸ ≤ 120–198 / 154–257` was proposed. **Later withdrawn (RU).**
4. **The V₄ sequence.**
   * V₄′ was false at O(1) for the raw `Ψ₄,ν`; the tube contact term grows like `ε^{−1/4}`.
   * V₄″ was false [DV], because of `e^{πν/4}` growth.
   * V₄‴ worked per ν with H_LT and H_νout, but left a circularity.
   * **V4F replaced all of them** by the pointwise Lemma A [D] and D_LR [D], deleted H_νout, and isolated `H_K`.
5. **PARITY_ODD_CLOSE / FWD.**
   * `+−` J = 5, OO J = 0, 2 and `+−` J = 6–11 were measured.
   * `H_OO^Q` [I] → [NE + M]; `H_fwd^Q` coverage [I] → [NE + M].
   * The minimum was unchanged.
6. **RU.**
   * Sound with corrections.
   * New or renamed: `H_EFT`, `H_RS` [I], `H_LJ^Q`; `H_K` [D-sketch].
   * Cubic closure withdrawn; `κ ≥ 0` required for `ĝ₃ ≠ 0`; `c₁t₄` carried; HYB [M].
7. **HRS.** `H_RS` regraded by window; the binding state is still [I]. The "512 c₂" item is the resolved `k = 2` comparison artefact
   (BF3 §2), not an open `k = 4` discrepancy.
8. **OG, FH.**
   * Lemma A mechanism [DV] at `O(1/C_T)` on free scalar and fermion; the Einstein instance is still [I].
   * `H_EFT` split into (a) [DV toy] and (b) [A].
   * The OG reading "`H_EFT` is an assumption only for `J_h ≥ 4`" is refined by FH: (b) is assumed for all `J_h`.

**Numbers unchanged throughout:** `2L_AdS = 72.04 / 92.81`; full-set min `+0.1426 / +0.1427` (δ = 0.2).

**Superseded statements of the previous body:**
* hypotheses 3 ("V₄ incl. the `[TT]₀` tail") and 9 ("`H_OO^Q` assumed");
* the regrading table's "`κ` [I]" and "Lemma A [D-sketch]";
* "Next steps" 1–3 (all done: smeared `k = 2` correction, swapped-form `Ψ₄` on MFT, W³W³ residue).

**Unresolved between documents.** The absolute de-shift moments of the sc3 weights are `2.2, 83, 1.4·10⁴, 3.6·10⁶` in V4F but
`1.49, 33, 3.5·10³, 6.6·10⁵` in FH ("different normalisation or key"). No bound depends on them.
