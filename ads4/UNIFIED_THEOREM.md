# Einstein gravity through quartic order in AdS₄/CFT₃ with a large gap: a conditional theorem

**Ansh Mishra and Aryan Senthilkumar**

*22 Nov 2026. Unrefereed. This document unifies Result D (`RESULT_D_STATEMENT.md`) and Result Q (`RESULT_Q_STATEMENT.md`). Its starting
point is `REVIEW_UNIFIED_2026-11.md` (RU) §5's recommended wording, updated with `HRS_CHECK.md`, `RESULT_Q_OG_CHECK.md`,
`RESULT_Q_FERMION_HEFT.md` and `RESULT_D_CLEANUP.md`. No new computation was made for it, and every number is quoted from the two
statements, which name its source.*

**This is a conditional theorem.**
* The implication "per-state positivity of three fixed functionals ⇒ the bounds" is [D], given a list of analytic inputs, themselves
  graded [Q], [D], [D-sketch] or [A].
* **Per-state positivity is not established in general.** It is computed [NE + M] on finite sets of heavy states and assumed elsewhere,
  through named hypotheses, two of which are open [I].

**Grade legend.**
* **[P]**: proved. **[certified]**: computer-assisted proof.
* **[D]**: derived. **[D-sketch]**: derivation sketched. **[DV]**: derived and checked numerically.
* **[Q]**: quoted from the literature.
* **[NE]**: numerical evidence. **[NE + M]**: computed inside the near-gap chain/HYB model.
* **[M]**: model. **[A]**: assumption. **[I]**: open.

A chain carries the grade of its weakest link.

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

> **Update (DG_SCAN.md, unreviewed): Δ_g dependence.**
> * **Bound values as explicit functions of Λ [DV].** The functionals are held fixed, and the forms reproduce every quoted digit at 100, 300 and 1000:
>   * cubic: `√(t₄²+α_T²) ≤ 1.3194·10⁴·√(ln(Λ/8.7415) + 6.384/Λ + O(lnΛ/Λ²))/Λ⁴`;
>   * quartic upper end: `17.2407·ln(2Λ/3) − 19.2782 + 10.314/Λ + (860.3·ln(2Λ/3) − 2508.6)/Λ² + …`.
> * **Where the constant drift comes from.** The shift −19.220 → −19.265 is the finite-Λ tail of the AdS pole, and the constant tends to the flat value −19.2782.
> * **Error terms in Λ:**
>   * c₁ (Result D): `∓19941/(Λ²√(16L_A))` per unit c₁;
>   * H17: `±212.5/Λ²` [M].
> * **Positivity threshold not computed.** The model bracket is `140 ≲ Δ_g* ≲ 300` [M]. The sign-free large-J window closes only from Δ_g ≈ 275, with margin ×1.1 at 300. So 300 is the lowest gap where the current method closes.
> * **Proposed corollary.** Under H_mono (corrections fall monotonically at fixed (J, m/Λ)) and H_cont (continuity in m/Λ), the bounds hold for all `Δ_g ≥ 300` [D-sketch + M]. Nothing is claimed below 300.
> * **Q's binding state is Δ_g-stable.** The law `s_d ≈ 0.50 + 1440/Δ²` was mis-fitted (REPRO D1). The cached data give `s_d ≈ 0.5076 + 147/Δ²` [NE/M], which is more stable than stated.
> * **Record correction.** RESULT_D_STATEMENT's c₁ entry at Δ_g = 100 should read `∓2.37 %`, not `≈ ∓1.7 %`, by the same formula as the 300/1000 entries.

> **Update (REPRO.md, independent reproduction with fresh code): every headline number reproduces within its stated precision.**
> * **Bounds and constants (fully independent):**
>   * cubic bound `3.13206·10⁻⁶` / `2.89126·10⁻⁸`;
>   * `2L_AdS` `72.0406` / `92.8135`;
>   * `B_eff` and the DG_SCAN Λ-forms.
> * **Error terms (fully independent):** c₁, H17 and δk(0).
> * **Binding states:**
>   * Q binding FF/N `+0.14270` / `+0.14282`. Partly independent: the engine ratio ρ(ν), the de-shift weights χ and `R_j` were taken as given.
>   * D binding states: hybrid `+0.20427` (300) and `+0.20090` (1000).
> * **ABJM:** placement `3.34 %` / `2.59 %`, with `8πG = 256/c_T` confirmed.
> * **Not independently reproduced:**
>   * the analytic-radial engine ρ(ν), which carries all of `s_d`;
>   * the Abel de-shift weights;
>   * D's crystal corrections;
>   * the NF02 flat certificate (grid check only);
>   * the full 14000-state Q scan.
> * **Record fixes applied:**
>   * H17 band `±0.24 %` (was `±0.2 %`);
>   * twist edge `1.00499Λ` (was `1.0048Λ`);
>   * the `s_d` vs Δ_g law corrected.
> * **Remaining minor items (none binding):** D4, D5, D7 and D9 in REPRO §2.

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

> **Update (HLJ_CLOSE.md, unreviewed): H_LJ(300) is not closed. Its Regge window for J ≥ 34 stays [I], and H_LJ(1000) stays [D-sketch + M].**
> * **Diagnosis.** The large a¹ term is a J-independent EE normalisation, of size 16–21/Δ relative.
>   * The genuinely J-dependent part is the helicity structure. It is small in the window (25–60 c₂ per order), but its 1/Δ series at fixed J/Δ is asymptotic, not convergent.
>   * No resummation was achieved [I].
> * **Analytic normalisation `ρ_an`.** It has no fitted constant, agrees with the old fits within their errors, and Prop. RS bookkeeping holds at J = 28, 32 [DV].
> * **Window states at Δ_g = 300, `ν < 8.3`: [M] → [NE + M].**
>   * (336, 28): `λ/N_s = +0.759`; stresses give 0.751–0.772.
>   * (336, 32): `+0.837`.
>   * The binding margin is the structure error, at 3–6× the engine's structure terms. Every measured discrepancy has the safe sign.
>   * (336, 28) at +0.759 is now the lowest computed large-J value at 300, below the previously quoted J = 32 +0.783.
> * **High J (J = 64, 100).** The engine's last term is 20–55 × T_B, so it cannot decide H_win either way.
> * **Unexplained:** a negative structure offset of 20–330 c₂ at ν = 5–10 [I].
> * **Net for Result D:** the only remaining [I] is the H_LJ(300) Regge window for J ≥ 34. The headline at Δ_g = 1000 stands (PRL draft switch = open).

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

> **Update (NEARGAP2.md, unreviewed): Result D's binding state is now error-bounded [NE + M] → [NE±]; light J = 0, 2 states are positive up to just below the gap.**
> * **Binding state** (J = 5, `+−`, twist edge 1.0049925Λ, Δ_g = 1000), on a K₂ row at its own mass (Δ = 1005, 36 nodes): λ/N_s = +0.20088. Error budget, sign-free with ×3 on every model term: 2.6·10⁻³, so **λ/N_s ≥ +0.198, leverage L = 78** (130 with the λ-level truncation). Certified flat floor minus the Lemma-C bound gives ≥ +0.197 with no quadrature input. 9 direct RQMC nodes agree with K₂ (λ moves < 10⁻⁶); the parity-odd analytic-radial engine agrees on ν = 12.5–38. At Δ_g = 300: ≥ +0.146, L = 3.5 [NE±], with no direct check.
> * **Near-gap J = 4, 6** (own rows): [NE±], L ≥ 22 / 5.0 sign-free. **J = 8**: [NE±] on the λ-level budget only (L = 4.3–8.4). **J = 12, 16**: [NE + M]; the budget does not close (L ≈ 1–3), and the next Laplace order is needed. The published chain differs from the own-row values by at most 0.015 N_s.
> * **H_cont, J = 15 twist edge**: computed, +0.6900 [NE + M] (was [M]). **H_edge**: J = 18, 19, 24, 25 shift by ≤ 0.032 [NE + M]. **J = 29 is not resolved:** one K₂ outlier at x = 0.8 drives the plain w2 odd chain to −3.62, while the direct nodes give +1.39 ± 0.03. It stays [M], and is carried into the high-spin task.
> * **Light J = 0, 2 (A1′):** positive at every Δ_L up to Δ_g − 1 in the leading-Regge model (F/N ≥ 0.788 / 0.638), and Lemma-C certified there for Δ_L ≥ 4.5. Own rows at Δ_L = 990 give L = 4362 / 913 [NE±]; Δ_L ≳ 0.7Δ_g gives [NE + M]. Whether the leading-Regge form is the true AdS image of a light block at Δ_L = O(1) remains [D-sketch].
> * **Correction to POS_ERRORS §3.1:** odd J has only the parity-odd `+−` channel, a 1×1 block; the sign-free tolerance at the binding state is 3.6·10⁴ c₂, not 6.6·10³.

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

> **Update (INDEP_CHECK.md, independent reproduction from the statements with new code; no project code imported or read): the H_K EE reduction is reproduced, and Lemma Y is reproduced for the count and at the correlator level.**
> * **EE kernel = Σ_{b=−4..4} shifted scalar waves:** holds to 2.6·10⁻¹² at 81 cases; end to end to 10⁻⁹–5·10⁻⁸. **Corner exponents** (w2g −0.50, w2s ≈ −0.4, w4g −0.50, w4p +1.00) and the **thresholds** J₀ < 3/2 (k = 2, 3, 5) and J₀ < 5/2 (k = 4) are re-derived. The w4g exponent stays at −0.500 out to β = 13 in 45-digit arithmetic, so the c = ½ input is conservative. **f_pos(±3i/2) = 0** for Q's published k = 4 continuum weights (contour ≤ 10⁻⁶⁸; decay-rate test against a control with a pole).
> * **Lemma Y:** the count of 997,836 = 35,637 × 28 is reproduced, worst margin 0. The closed form σ = (g+e+f) − (z+e′+f′) + (k₁₃−k₂₄)/2 makes it [P] for any P·P content. Five explicit bulk spin-2 contacts (833 polynomial × D-function terms each) behave as A + B log y down to y = 1.2·10⁻⁶.
> * **Not independently reproduced:** the step from the correlator to the light-function bound for L6′ (it stays [D] on the original work); the non-dispersive counterexample slopes; HK_RESIDUAL §1's k = 4 kernel normalisation; K1′, K2ν and LRC-a.
> * **Minor corrections:** the failure of EESTEP §1.3's negative-b series is truncation, not round-off (not load-bearing); HK_CLOSE and RESULT_Q_STAGE3 differ by one factor of η in the raw k = 2 corner entry (both say "fails"); one displayed line in R4H §2.3 is garbled, and the closed form replaces it.

> **Update (LARGESPIN_ANALYTIC.md, unreviewed): no large-spin hypothesis is replaced by a proof. H_R′ is upgraded [M] → [D-sketch + M]; H_LJ(1000) and H_bulk keep their grades. A possible convention problem (N2) affects H_LJ's leverage.**
> * **Common missing input:** the beyond-leading-Regge correction inside the Regge window (ν ≲ 3/β). It needs the ⟨TTTT⟩ heavy-block Regge formula at subleading order (the spinning analogue of CHMRS (3.17)–(3.20)), which is not derived.
> * **Lemma I [P]:** for an image ∫P_{−½+iν}(cosh t)σ(dt), the action is ∫Ω_LR(t)σ(dt). This gives explicit tolerance profiles for H_LJ (≥ 0.062 of the ring weight at t = β) and tolerated sizes of candidate correction shapes (ring amplitude a ≥ 80/m, ν² a ≥ 37/m², index-shifted a ≥ 10.8/m, …). The remaining model input R_win (non-ring window correction): exact scalar ≈ 10 c₂ (×38); graviton J = 100 ≈ 550 c₂ (leverage 1.1–1.5); no data at J ≥ 150.
> * **H_R′ [D-sketch + M] for b̂ ≥ b̂_A(m).** Lemma C′ [P + DV], a sign-free curvature bound using the decay of Φ_k, makes the hybrid floor positive for b̂ ≤ 119 at every m/Λ from 1.3 to 20 (the old Lemma C reached b̂ < 10). With ROUTEA, the leading-Regge hybrid is rigorous for all b̂ ≤ 2Λ when m ≥ 1.5Λ. The AG saddle test passes on 5 states. 5 ≤ b̂ < b̂_A stays [NE + M].
> * **H_bulk [M], restated:** the bulk correction may reach ≥ 3.9 % of the flat envelope (9.8 % at 1.3Λ) before any state fails. There is no analytic route.
> * **N2 (open, may weaken H_LJ):** Ω^EE's k = 2 part acts on raw images with the de-shift weight X, not ψ₂, and |X−ψ₂|/ψ₂ = 0.3–0.9 at ν ≤ 10. HLJ's T_B used ψ₂ inside the window. With |X|, the leverage at (1100, 100) drops to ≈ 1.1, and the sign-free ring-amplitude tolerance (3.2) is below the measured excess (2.7–5.9). Which convention the fw/pw Γ data use must be checked. **Until then quote the H_LJ(1000) leverage as 1.1–4.8 (not 2.7–4.8).**
> * **Corrections:** ROUTEA's certified window [91.2, 2Λ] holds only for m ≥ 1.5Λ (UT3 §5.2 and PE §2.4 quote it for m ≥ 1.3Λ). HLJ §5.2's [D-sketch] (subleading terms are P_{l+n}(cosh β) of the same β) does not by itself make window corrections harmless (N3). PE's "no sign-free argument can work" is true for δ but false for the curvature.

> **Update (N2_SUBREGGE.md, unreviewed): N2's premise is right (X, not ψ₂, multiplies the window data), but the claimed leverage collapse does not happen. A dir–cv estimator systematic at (1100, 100) is now the dominant uncertainty.**
> * **What the rows are [DV].** Window rows Γ are raw images / (R₋+R₊), and X multiplies them; Re ΣX𝒜[P_β] = Σψ₂P_β holds to 2·10⁻¹². Γ has poles at the real zeros of R₋+R₊ (ν = 1.171, 6.643), which explains HLJ_CLOSE's "local Γ singular at every J". H_win should be stated on the regular aΓ_c with weight |Y|.
> * **T_B with the exact weight:** (1100, 100): 1409 → 1439 (|X|) / 1152 (|Y|). Over all sample states the change is −28 % to +3 %. LARGESPIN's "1.1" and "3.2" were artefacts (ψ₂ + |X−ψ₂| overcounts; a different floor and structure were switched in; signed ring tolerance 31 vs measured 4.5–6.6). The D_h excess is an a_h normalisation, already in de-shift form (cross-checked against k = 3, which shares a_h). The published Γ − κ contains a stripped-κ artefact of ±120–460 c₂.
> * **Corrected H_LJ(1000) leverage at (1100, 100)** (cv estimator, raw-form refit): **3.0–4.2 central; ≥ 4.1 at 2σ with the full ring class.** (364, 64): 20–36. **With the dir estimator: 0.8–1.4.** This dir–cv systematic is unresolved and dominant.
> * **New explicit hypothesis R1 [I] for EE beyond leading Regge:** ring-type window corrections enter in de-shift form. It holds exactly for the a_h mechanism and for the scalar. In stripped form they would cost 27–228× the ring action at J = 150–344.
> * **Subleading Regge:** CHMRS (3.20) is checked against the exact scalar skeleton (9 states). All its terms sit on the state's own ring, with cost 6.5·10⁻⁴ (J = 100) to ≤ 6·10⁻⁵ (J ≥ 150) of the floor, so margins are ×17 to ×10⁴. The graviton window data fit the same ring classes, with non-ring remainder consistent with zero. The **spinning formula is not derived [I]** (obstacle: the EE = Σ_b identity decomposes the kernel, not the spinning block; it needs ⟨TTO⟩ structure ratios). At each finite 1/m² order the corrections are ring-located or index-shifted at the state's β [D-sketch].
> * **Corrections:** the "leverage 1.1" in the LARGESPIN box above is superseded. Quote "3.0–4.2 (cv), ≥ 4.1 at 2σ; 0.8–1.4 with dir". HWIN_THIMBLE: refit the D_h class in de-shift form, after which the (364, 64) remainder is 7–13 c₂.

> **Update (HIGHSPIN.md, unreviewed): J = 29 is resolved and positive; near-gap J = 12, 16 move to [NE±]; no physical state came out negative.**
> * **J = 29 twist edge (1.0286Λ):** λ/N_s = **+0.972 ± 0.018** [NE + M] (was [M]). Method: a full direct n1 row (34 nodes) plus GK15 quadrature of 16 direct edge nodes, with no interpolation. Budget 0.70, so L = 1.39. The one model input is the ν < 50 window: flipping the state needs |Γ| ≥ 3.2·10³ c₂ there with the wrong sign, while the measured values at ν = 51–77 are ≤ 116 c₂. Masses 1.04–1.10Λ: +0.956 … +0.995, L = 2.8 → 9.5.
> * **Diagnosis of the −3.62 [DV]:** the K₂ method fails, not the physics. The Δ_ref = 2·10⁴ reference's leading partial wave has its first zero at x₀ = 0.8008, right at the bad node (|r₂/r₀| = 0.86 vs ≈ 0.02 elsewhere). w2's Hampel filter skips odd J, so the bad node is kept. The +1.39 is spline aliasing. The published +0.948 removed the node, and the shape chain is within 0.025 of the direct result.
> * **The same reference-zero problem in published chains:** J = 31 (+1.208 → +1.02; this one drops, though it stays above the published minimum +0.988), J = 32 (+0.904 → +0.997), and the J = 40 and J = 48 outliers explained. (Correction: not all shifts are favourable; J = 31 drops.)
> * **J = 12:** +0.8528, L = 3.6 / 3.9 / 2.5 (twist edge, 1.02Λ, 1.05Λ). **J = 16:** +0.9402, L = 1.9 / 1.7 / 1.4 / 1.2. Both now **[NE±]**, from full direct rows (38 nodes) that remove the truncation term. Odd J = 21, 23, 27 and even J = 20, 22, 28, 32 edges: [M] → [NE + M]. Not run: J = 26, 30, 48, 10, 14. **H_LJ(1000): unchanged [D-sketch + M].**
> * **Cautions:** the K₂ λ-level Laplace step cannot close a ×3 budget for J ≳ 17, although the actual K₂ error is only 4–5 % of the step and always favourable. The ν < 50 window is now the binding model input at large J. Lowest computed values: J = 21 at 1.03Λ, +0.868; J = 20 twist edge, +0.83 ± 0.14. Code bugs worked around: `p4_odd_analysis.grid` and `l9_analysis.grid` drop direct nodes whose x is not round to 6 decimals.

> **Update (DIRCV.md, unreviewed): the dir–cv discrepancy is resolved. Both estimators are biased, in opposite directions, by one per-wave systematic. H_LJ(1000) leverage at (1100, 100): 2.3–3.4 central, 1.3–2.0 at 2σ (ring classes A + C); with the full ring class (R1), 7–9 central and ≥ 2.7 at 2σ. A leverage of ≈ 1 is excluded at ≥ 2σ in every fit [NE].**
> * **Errors:** the quoted r errors were 1-dof least-squares residual errors. Replica-jackknife errors are 2–20× larger: cv −160 ± 213 / +327 ± 216, dir −1310 ± 526 / −810 ± 529, paired gap 2.4σ. The gap comes entirely from the node at ν = 13.3, where seed 4242 is a 1.4σ fluctuation (4-seed scalar bias −0.0075 ± 0.0021).
> * **Bias:** the CV mean is exact (≤ 10⁻¹⁰), but the graviton carries λ = 0.5 ± 0.25 of the scalar's per-wave bias (three tests: 0.55, 0.58, 0.43). dir assumes λ = 0 and cv assumes λ = 1, so each is disfavoured at ≈ 2σ. The unbiased P_J-cv twin agrees with the pooled result. The k = 3 "miss" disappears with jackknife errors. No heavy tails were found (Hill index 2.6–3.2).
> * **Pooled fits** (7× statistics at ν = 13.3 and 16.3): A + C, 4 nodes: r = −441 ± 180, L 2.3–2.8 (2σ 1.3). A + C, 6 nodes: +147 ± 77, L 2.7–3.4 (2σ 2.0). A + B + C: −102 ± 122, L 7.2–9.0 (2σ 2.9). J ≥ 150 is not feasible (no basis rows; windows at or below the fold), so it stays [M].
> * **Corrections:** this supersedes the leverage quoted in the N2 and LARGESPIN boxes ("3.0–4.2 (cv), ≥ 4.1 at 2σ; dir 0.8–1.4", "1.1–4.8"). Quote "2.3–3.4 (1.3–2.0 at 2σ); ≥ 2.7 at 2σ given R1". The cv window values in HWIN carry a systematic of ≈ 0.5·|Γ_dir − Γ_CV|. Next cheap step: 4× P_J-twin statistics at ν ≈ 13 (≈ 1.5 CPU-h), giving σ_r ≈ 90 with no λ assumption.

> **Update (DIRCV.md §7, unreviewed): with 4× more P_J-twin statistics the error on r did not shrink. The twin has a heavy tail at ν = 13.3 (tail index ≈ 2.3; one replica sits ≈ 8× from the typical one). Without a λ assumption, the 4-node A+C fit reaches leverage 0.8 at 2σ, so leverage ≤ 1 is no longer excluded in every fit.**
> * **P_J-cv, no λ, worse helicity:**
>   * A+C, 4 nodes: r = +534 ± 387, L = 1.8–2.3 central, 0.8 at 2σ.
>   * A+C, 6 nodes: +289 ± 124, L = 2.7–3.4, 1.7 at 2σ.
>   * A+B+C: +181 ± 231, L = 4.7–6.0, 1.6 at 2σ.
> * **λ-corrected numbers unchanged** (combined λ ≈ 0.6 ± 0.25): A+C 2.3–3.4 central (1.3–2.0 at 2σ); A+B+C 7–9 central (2.9 at 2σ). Every central leverage is ≥ 1.8.
> * **The earlier 56-replica P_J-cv cross-check was not robust** (e.g. `++` 4-node −393 ± 174 → +23 ± 380); the bracketed numbers in the DIRCV box above are superseded.
> * **Next:** add window nodes at ν = 11–15 so the fit leans less on ν = 13.3, or use a variance-bounded twin at that node. More plain samples are inefficient.

> **Update (DIRCV.md §8, unreviewed): 5 new window nodes at (1100, 100) (ν = 11.0–15.0, avoiding the D_h zero at 13.09). The A+C fit no longer hinges on ν = 13.3: dropping that node moves r by +6 to −149, inside the errors. Weighted fits (χ²/dof 0.8–1.3) are now required. λ = 0.7 ± 0.25.**
> * **Leverage, worse helicity, central / 2σ:**
>
>   | fit | λ-corrected | λ-free P_J-cv |
>   |---|---|---|
>   | A+C, 9 nodes | 2.2–2.7 / 1.4 | 2.0–2.6 / 1.0–1.3 |
>   | A+C, 11 nodes | 3.1–4.5 / 1.8 | 2.2–2.9 / 1.7 |
>   | A+B+C, 9 nodes | 2.7–3.8 / 1.4–1.7 | 1.6–2.1 / 0.5–0.6 |
>   | A+B+C, 11 nodes | 4.2–5.6 / 2.3–2.5 | 5.1–7.4 / 2.1–3.2 |
>
> * **Every central leverage is ≥ 1.6, and every A+C fit is ≥ 1.0 at 2σ.** The only 2σ entry below 1 is λ-free A+B+C on 9 nodes, where the ν²D coefficient rests on noisy P_J nodes. The 2σ floor is set by P_J noise at ν < 13 (fold side). Grade [NE].
> * The one-contour twin does not reduce variance. P_J-dir is unusable (point drops at ν < 13; only the control variate cancels them).

> **Update (SMALLNU_WINDOW.md, unreviewed): the ν < 50 window is now measured. The 10⁴ c₂ convention is kept only below the lowest measured node (6–18 % of the window weight). Nothing came out negative. J = 29 and J = 21 → [NE±] (thin); J = 12, 16 leverage roughly doubles; J = 23, 27, 31 stay [NE+M] (now limited by the K₂ bulk truncation, not the window).**
> * **Method:** direct n1 points at ν ≥ 20 and n3 points at ν = 10–17.5 on each state's own row, with jackknife errors. n1 fails below ν ≈ 20 (crystal-saddle fold); n1 and n3 agree at ν = 25 and 30. Reference zeros cannot occur in the window.
> * **J = 29 twist edge:** measured Γ ≤ 310 c₂ on ν = 25–47; ν = 10–20 is consistent with 0 (σ up to 730). Signed window: L 0.61 → **1.23**, λ ≥ +0.18. The remaining model input is |Γ| ≤ 10⁴ on ν < 10; the state tolerates 4.5·10⁴ there, 16× the largest known value. Masses 1.04–1.10Λ: L 2.2–5.7. **A sign-free ×3 envelope does not close (L = 0.94)**, because of Monte Carlo noise at ν = 10–20.
> * **J = 21:** a new 34-node direct bulk row puts the twist edge at +0.998 ± 0.020, L 0.58 → 1.20 (1.7–6.2 at 1.03–1.10Λ). **J = 16:** L → 4.4 / 3.7 / 2.6 / 1.7. **J = 12:** L → 7.4 / 7.5 / 3.6.
> * **J = 23 / 27 / 31:** window cost now 0.15 / 0.29 / 0.46, but the K₂ λ-level truncation blocks them (with calibrated truncation L ≈ 1.5 / 1.4). **J = 32:** the window costs 0.26 measured, vs 2.31 under the convention. J = 20 and J = 32 have no budget.
> * **Cautions:** the window weight grows steeply with J (1.4 → 23·10⁻⁵ N_s/c₂ from J = 12 to 32), so for J ≳ 32 it is a live input again. Lemma C/C′ cannot bound Γ (it is beyond leading Regge). HIGHSPIN's odd-J window weight omits a gauge factor (8.52 vs 8.73·10⁻⁵). EE_DG1000_LARGEJ's "direct RQMC cannot reach the precision" is false on ν = 25–50.

> **Update (LARGESPIN_W.md, unreviewed): Route W gives a conditional theorem with an explicit J\*, not an unconditional one. H_LJ(1000) stays [D-sketch + M]; do NOT read it as "proved for J ≥ 50".**
> * **Theorem W** (implication [P] given P0–P4; costs [DV] on sampled states). Premises:
>   * P0: the leading-Regge hybrid floor.
>   * P1: parity, M = 𝔫U†gU with g even in 1/m at fixed β; [D-sketch].
>   * P2: beyond-hybrid O(1/m²) corrections are ring class with |κ_s| ≤ Λ_amp·|κ_s(scalar CHMRS)|, and index-shifted harmonics have |κ| ≤ K_sh. Magnitudes are [NE] at J ≤ 100 and [M] beyond.
>   * P3: R1, i.e. de-shift form; [I] at O(1/m²) for the ν² classes.
>   * P4: O(1/m⁴) ≤ K₄, plus the AG beyond-all-orders hypothesis; [M].
>
>   If the summed Lemma-I cost d ≤ H = f₊₊f₊₋/(f₊₊+f₊₋), then M_Ω ⪰ 0 for every coupling ratio.
> * **J\*** over the samples (even J 26–344, odd J 33–301, m ≤ 1.3Λ):
>   * with (Λ_amp, K_sh, K₄) = (3, 10, 10) and R1: J\* = 50 at margin ×1, 100 at ×2;
>   * with (1, 10, 1): 30 / 38;
>   * if R1 fails (stripped form), there is no J\*: it fails at J = 300 and 344.
> * **New pieces:**
>   * W1a [D+DV]: the 1/m parity comes from a Weyl reflection T(−x,−y; Δ,J) = T(x,y; 3−Δ, −1−J), holding to 10⁻⁴⁰.
>   * W1c [P], congruence lemma: ν-independent normalisation and coupling-basis effects cannot change positivity. This disposes of the measured +21/m term and N2's a_h excess.
>   * W2 [D] / [D-sketch]: every finite order is ring class, so there is no window-smooth class.
>   * R1 at O(1/m) is no longer needed given P1.
>   * The leading term c(ν) = 32γ₃²γ₄/(γ₁γ₂γ₇) is exact.
>   * The (1100, 100) B fits correspond to Λ_amp ≈ 2–3 and sit ×10–100 inside tolerance.
> * **Decisive residuals:**
>   * the P2 magnitudes, which need the next-order EE-contracted expansion of the integrand at fixed β;
>   * R1 at O(1/m²) for the ν² classes at J ≳ 200. A sufficient, cheap check: the s = ±1 harmonics of the leading-r EE kernel share one polarisation tensor (`eestep/ee_quad.py`).

> **Update (LARGESPIN_W.md §9, unreviewed): neither premise of Theorem W closed.**
> * **R1's sufficient condition is false [DV].** The leading-r EE kernel has rank exactly 3 in polarisation (singular values 1, 0.24–0.37, 3·10⁻⁴–4·10⁻³; the next is at the 10⁻¹¹–10⁻¹⁶ floor). At leading order the heavy symbols see only the de-shift combination. At O(1/m²) the two other channels can enter in non-de-shift form; P1 forbids that at O(1/m). R1 stays [I], now pinned to these two kernel channels. If they are index-shifted harmonics with rational coefficients, they are tolerated even in stripped form (K_sh ≥ 52 at J = 344).
> * **P2 magnitudes not computed.** Three steps remain:
>   1. conserved ⟨TTO_{Δ,J}⟩ structures in d = 3 for general J;
>   2. second-order heavy symbols at fixed β;
>   3. contraction with the three kernel channels.
>
>   Λ_amp ≈ 2–3 stays [NE] (from J = 100); K_sh and K₄ stay [M].
> * **J\*:** with P3, J\* = 50 (100 at ×2). Without P3, positivity holds on 50 ≤ J ≤ 250 (and J = 301) but fails at J = 300 and 344. So R1 at O(1/m²) is load-bearing only for J ≳ 300.

> **Update (LARGESPIN_L.md, unreviewed): Route L gives no unconditional theorem either. It proves the leading-Regge floor at large impact parameter, gives a conditional Theorem L with J\* ≈ 0.52m (857 at the twist edge), and reduces the residual to two CFT-data-free statements about the EE functional's own plane-wave images.**
> * **Lemma LR∞ [P + DV]:** F_LR ≥ (1 − r)·A₁M₁/Λ² + 2047/(2b̂³) > 0 for all β ≥ 1, all m ≥ Λ and both helicities, with r ≤ 0.58 (β ≥ 1) down to 1.6·10⁻³ (β ≥ 3). The constants are float on a β grid to 64; the tail is [D]. With Lemma C′ and ROUTEA, the leading-Regge floor is rigorous on every b̂ for m ≥ 1.5Λ. This closes the old "b̂ > 2Λ: hybrid [NE]" entry and supplies Route W's P0 for J ≥ 857.
> * **Spin-Gamma representation:** [P] for the scalar heavy factor; [D-sketch] for ⟨TTO⟩⟨OTT⟩. Every heavy block is a positive average of plane waves, with J entering only through a Gamma(J+1) average.
>   * Corollary [D]: R1 holds automatically for the heavy-factor part of any correction. **R1 can fail only through the EE kernel's own plane-wave images, a property of the functional, not of the CFT.**
>   * Lemma W [P given Paley–Wiener on H₂]: a correction of exponential type T has its density within T of the ring.
> * **Theorem L [D given P_L]:** for J ≥ J\*(m) = ⌈(1.0421m − 1)/2⌉, the per-state matrix is ⪰ 0 for every coupling vector.
>   * Premise P_L, part (R1): de-shift form.
>   * Premise P_L, part (Loc): the de-shifted correction's density is ≤ tol(1)/4 at t ≥ 1 and ≤ 10⁻⁷–4·10⁻⁹ at t < 1.
>
>   P_L reduces to two open statements: **(PW-R1)** the EE functional's plane-wave images are in de-shift form, and **(PW-Loc)** they carry no density at t < 1 [I].
> * **Negative findings:**
>   * Corrections do not decay like O(1/J) at fixed twist; they tend to a J-independent lightcone amplitude (scalar ≤ 1.2·10⁻⁴ on ν ≤ 10).
>   * The de-shift tolerance has its minimum of 6.7·10⁻⁴ at J ≈ 350–500, so large J is not easier.
>   * In raw form the tolerance falls like 0.7/J.
>   * ANEC positivity holds only for the sum over states.
>   * Graviton non-ring remainder 2–5·10⁻⁴ (relative) at J = 29–100, 40× the scalar; the sign-free L ≈ 1.3–3 in the hardest band if it stays J-independent.
> * **Inventory:**
>   * J ≥ J\*(m): [D given P_L], with a rigorous floor.
>   * 150 ≤ J < 857 near the gap: the hardest band, still numeric/model.
>   * Route W covers J ≥ 50 given scalar-size coefficients and R1. Both routes stall on the same object: the EE-kernel polarisation contraction beyond leading order.

> **Update (P2_SPINNING.md, unreviewed): steps 1–2 done; step 3 done only inside a labelled model [M]. P2 is not closed, and H_LJ(1000) stays [D-sketch + M]. Do not read this as "Λ_amp = 1 derived".**
> * **Step 1 [DV + Q]: conserved ⟨TTO_{Δ,J}⟩ in d = 3, solved exactly over ℚ(Δ,J).** Generic even J: 2 parity-even + 1 parity-odd structure; odd J: 1 parity-odd. This matches Dymarsky–Kos–Kravchuk–Poland–Simmons-Duffin 1708.05718; conservation checked to 4·10⁻²⁷. The odd-J "+−" structure is parity-odd.
> * **Step 2 [D + DV]: heavy symbols to O(1/m²).** The leading symbols are pure helicity. The shadow side is exactly m → −m. O(1/m) is a basis drift plus a Wigner phase. **P1 is now [D given (S)]**, where (S) is that the shadow coefficients are diagonal in the helicity basis with positive leading entries [D-sketch]; check against 1809.05111.
> * **Step 3:** the CHMRS ∂_p mechanism reproduces CHMRS's O(1/m²) scalar series [DV] (ν⁰P off by exactly 8η, from normalisation). By-product [NE]: N2's unexplained "−8…−9·P" scalar residual is −8ηP at all 8 states. **Model [M]** (helicity-diagonal kernel): the spinning symbol adds nothing to the ν² classes, giving Λ_amp(B) = Λ_amp(C2) = 1, Λ_A = 3.17, Λ_C3 = 2.03. In the conjugate pairing, Λ_C1 ≈ 205. The (1100,100) and (364,64) fits bracket the model (inconclusive).
> * **Theorem W rerun (model constants):** J\* = 50 with R1 or in stripped form at K_sh = 10, K₄ = 10; 30 at K₄ = 1; fails at J = 300 and 344 if K_sh ≥ 30 without R1. **Now binding:** K₄ (below J = 50) and K_sh (J ≥ 300), both [M]. In the model, R1 is not needed at margin ×1. Helicity leakage into the non-de-shift channels is ≲ 10² at O(1/m²), possibly 1.8·10³ in the even basis; estimate κ ≲ 30 [M] against a tolerance of 52–59 at J = 300–344.
> * **Remaining:**
>   * the true three-channel kernel contraction at O(1/m²) (overlaps PW_IMAGES);
>   * K_sh from the kernel's O(r) terms;
>   * K₄ to O(1/m⁴) (cheap);
>   * the (S) check;
>   * the J = 0 engine O(a²) comparison.

> **Update (PW_IMAGES.md, unreviewed): R1 is replaced by an exact identity plus one bounded monopole term. PW-R1 and PW-Loc are false as literally stated, but neither is needed any more. Theorem W keeps J\* = 50 (100 at ×2) WITHOUT P3. Theorem L keeps J\* = 857 with P_L reduced to magnitudes.**
> * **Exact de-shift identity for the EE kernel.** Every b-shifted scalar block in EE = Σ_b satisfies CHMRS's magic-block identity with b-independent coefficients:
>   2η(1+r²)Ĝ¹_b(ν) = (1+i/2ν)Ĝ⁰_b(ν+i) + (1−i/2ν)Ĝ⁰_b(ν−i).
>   It is [P] to r⁸, [DV] to r⁴⁰ at 50 digits, [D] at all orders. So every plane-wave image, in any gauge and with any insertion, is exactly 𝒜_g[F] (CHMRS's scalar de-shift), with F even and analytic for |Im μ| ≤ 3/2. End to end against ee_quad: 2·10⁻⁷–1.1·10⁻⁶. **Route W's stripped (R2) worst case cannot occur.**
> * **PW-R1 is refuted** with respect to the EE operator 𝒜 = 𝒜_g∘c. It needs F(±i/2) = 0, which O(1/m²) radial insertions violate [P]: |y|² gives a ratio of 0.914, r gives 0.500. Whether the full ⟨TTTT⟩ sum cancels this is [I].
> * **The failure is exactly a monopole term, and it is bounded [DV].** The exact cost is ∫ψ₂F/c + D·F(i/2), with D = −9.98·10⁻⁴ (1000) / −1.10·10⁻² (300), ≈ −990/Λ². D is constant in impact parameter, and the equivalent cost kernel stays below 10⁻³ (no e^{t/2} growth).
> * **PW-Loc is false at O(1/m²):** index-shifted harmonics carry density 0.271 per unit coefficient at t < 1 (closed form [P]). It is no longer needed, since every finite-order structure now has an exact cost.
> * **Theorem W without P3 [P given P0, P1, P2, P2′, P4] + [DV]:** the scalar-shaped worst case is tolerated at ×24–29 for J ≥ 150 (stripped form gave ×1.8–7). J\* = 50 at (3, 10, 10); J\* = 100 at ×2. P2′ (the radial-shift family exhausts the O(1/m²) multipliers) is [D-sketch].
> * **Theorem L [D given P_L′]:** J\* = 857 at the gap. The monopole tolerance gives a tolerated Λ_amp of 1240 (J = 857), 555 (10⁴), 166 (10⁵) and 43 (10⁶).
> * **Still needed by both:** the spinning O(1/m²) magnitudes (Λ_amp ≈ 2–3 [NE] for J ≤ 100, [M] above) and the beyond-all-orders premise. Untested: polarisation-tilt O(1/m²) terms.

> **Update (P2_SPINNING.md §7, unreviewed): K₄ ≤ 3.8 in the model; the (S) check FAILS, so the P1 upgrade is withdrawn and P1 stays [D-sketch]; a new gap in P4 (the flat ν⁶ tower) is found. Model J\* = 38.**
> * **K₄:** the scalar seed matches CHMRS's O(1/m⁴) series exactly (ν⁴, ν⁶ coefficients at two β). In the helicity-diagonal model the spinning symbol adds only −(5/24)(η²−1)ν⁴P, so K₄ ≤ 3.6 + 0.21 ≤ 3.8 [M]. There is no model-independent bound.
> * **New gap in P4:** CHMRS's O(1/m⁴) term contains −(η²−1)ν⁶P/72, which is the same for every Δφ (part of the flat-space tower). Costed as a Regge term it is equivalent to K₄ ≈ 35 (J = 26) to ≈ 3000 (J = 120–200), and every sampled J fails. **P4 needs this tower's hybrid subtraction made explicit:** the identification is [DV], the subtraction is [I].
> * **(S) fails:** Karateev–Kravchuk–Simmons-Duffin 1809.05111 shadow matrices mix structures, so the two parity-even structures can mix. P1 stays [D-sketch]. What survives: the leading shadow matrix is diagonal and positive [D], and in the model the O(1/m) off-diagonal is constant, so it cannot affect positivity.
> * **−8η factor:** the first-order term of the moment's normalisation, exact at O(1/m²) and O(1/m⁴). It is harmless (a rescaling); its sign relative to N2's residual is [I].
> * **J\* (model constants, K₄ = 4):** K_sh = 10 gives 38 with R1 or stripped (θ = 1), 50 at ×2 with R1. K_sh = 30 gives 38 with R1 and none stripped. **This holds only if the flat ν⁶ tower is assigned to the hybrid floor, now a named premise.** K_sh is uncomputed [M] and binding at J ≥ 300.
> * The J = 0 engine O(a²) comparison was not done (it needs the full three-channel contraction).

> **Update (REGGE_TREE.md, unreviewed): referee-3 M2 is a presentation defect (the order of limits was never stated), not a derivation gap. Done as CHMRS do it, the Result never needs H_K at tree level; all grades stay.**
> * **H_K's J₀** is the exponent of the uniform bound on dDisc of the **exact** correlator (KKSDZ (4.70)). **Route N [D]:**
>   1. Apply the sum rule at finite c_T. H_K is [D] under J₀ ≤ 1, with margin ½ at k = 2, 3, 5 and 3/2 at k = 4.
>   2. (Pos) makes the heavy side ≥ 0 at every c_T, multi-trace states included.
>   3. Expand only the light side in 1/c_T.
>
>   The tree-level heavy sum is never expanded; CHMRS §4.4 show it need not be analytic in 1/c_T. The order-of-limits step is [D-sketch] given A2. A2 must be restated as tree-level EFT light data up to o(1/c_T), uniformly in spin.
> * **Term-by-term reading (not made):** with a uniform planar J₀ = 2, the corner count fails at k = 2, 3 **and 5**, and for generic k = 4 weights (Q's k = 4 weights pass with margin ½). "J₀ = 2 (Einstein)" is the Wightman/elastic intercept; the heavy tree dDisc is empty for pure Einstein, and for strings its fixed-η intercept is 2 − 9/(2Λ²) with an η⁻² k = 2 moment (CHMRS (4.45)).
> * **Literature [Q]:**
>   * CHMRS (§1.1.1, (1.12)–(1.14), §4.1, §4.4) and CLS-D (§3.4, §3.5.2) use J₀ ≤ 1 non-perturbatively and assume no tree-level J₀ < 2.
>   * KKSDZ §4.5 gives J₀^planar ≤ 2 (4.80); their planar sum rules need J₀ < 3.
> * **Model check:** with the uniform J₀ = 1 bound, the corner integrals converge for every k; with the uniform J₀ = 2 bound they diverge at k = 2, 3, 5; with a stringy Regge-pole model they converge for every k.
> * **Wording:** proposed replacements for letter.tex:109 (A2), :133–135 and :183–185, plus supplement paragraphs "Order of limits" and "Term-wise use (not made)", are in REGGE_TREE §5. Say "every state with Δ−J ≥ Δ_g", not "heavy single traces".
> * **Open:** the size of the light-side remainder O(Δ_g^p/c_T²) is not computed, so keep "at leading order in 1/c_T".

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

> **Update (LRINF_TIGHT.md, unreviewed): Lemma LR∞ is certified down to β₀ = 0.37, so Theorem L's J\* at the twist edge goes 857 → 223. The 345–857 band is now inside Theorem L′ at every J (not only sampled J), and Theorems W (J\* = 64 with k = 3, sampled to 344) and L overlap at the gap on 223 ≤ J ≤ 344.**
> * **Constants [certified]:** ∫|ψ₂| = 10.386029907129 and ∫|ψ₃| = 12.835622283960 (Arb, ±7·10⁻²⁰; Sturm-counted sign changes). The old table used Σ|c_n|/(n+1) = 7.9·10⁵ / 1.9·10⁶, which explains the 0.577. A₁,min and A₂′ are exact rationals.
> * **r(β) [certified]:** monotone pieces evaluated in Arb on cells covering [0.1, 64], plus a new explicit tail for β ≥ 64 with r ≤ 2.03·10⁻⁷ (was [D]). sup r = 0.0162 for β ≥ 1 and 0.280 for β ≥ 0.5. **β₀ = 0.37 for r < 1.** F_LR > 0 is certified for β ≥ 0.245. Step 3 (the Hadamard identity) is now [P]; Lemma LR∞ is [P + certified].
> * **New J\*(m) = ⌈(0.37211·m − 1)/2⌉ [certified]:**
>
>   | m | old J\* | new J\* |
>   |---|---|---|
>   | twist edge | 857 | 223 |
>   | 1.1 m_gap | 987 | 250 |
>   | 1.3 m_gap | 1275 | 307 |
>   | 1.5Λ | 780 | 279 |
>   | 2Λ | 1041 | 372 |
>   | 5Λ | 2602 | 929 |
>
>   At J\* = 223 the certified floor is 0.57 of the numerical floor.
> * **Tolerances at J\* = 223 (certified floor) [DV]:**
>   * monopole Q_tol = 4.8·10⁻³;
>   * tolerated Λ_amp: 181 for the full scalar-shaped vector, 1025 for the monopole alone; ≥ 181 up to J = 10⁴, 133 to 10⁵, 32 to 10⁶;
>   * K_tol for sh±1 ≥ 449.
> * **Floor coverage:** the floor is rigorous at every b̂ for m ≥ 1.5Λ. For m ∈ [1.3, 1.5)Λ the [NE] hole shrinks to b̂ ∈ (119, 245). For m < 1.3Λ the floor is [NE] below b̂ = 245. Theorem W's P0 is now certified at its sampled gap states with J ≥ 139.
> * **Still uncovered:** non-sampled 64 < J < 223, J < 64, and each theorem's premises (P_L′ magnitudes [M], beyond-all-orders). This is not an unconditional theorem.

> **Update (P2_CONTRACT.md, unreviewed): the helicity model is replaced by an exact three-channel computation at all 47 Theorem W sample states. All of P2 (Λ_amp, K_sh, P2′, P3) is now computed [DV]. P4 (beyond-hybrid O(1/m⁴) and beyond-all-orders) is the only magnitude premise left, and it is still open [I]: the state-wise scalar O(1/m⁴) remainder costed to ν ≈ m exceeds tolerance at every sampled J ≥ 80 except 301.**
> * **Engine:** an exact leading-Regge EE plane-wave image engine with any polarisation and insertion (all three channels). It reproduces the heavy images (ratio 1.000000), the R∓ normalisation, and CHMRS ∂_p to 10⁻¹⁵.
> * **Computed constants:**
>   * Λ_B (block) = Λ_C2 = 1 exactly.
>   * The kernel's own O(r) term is ring class with κ_B = ±2, giving **Λ_B(+−, odd J) ≈ 2.0 (worse than the model)**; Λ_B(++) ≤ 0.09.
>   * C1(+−) ≈ 200× scalar (harmless); Λ_A ≈ 6–8 (a rescaling).
>   * Off-diagonal ++/+− at O(1/m²) = 0. The O(1/m) image after drift removal = 0, so P1′ is [DV] at leading-Regge kernel level and the (S) failure cannot act.
>   * **The monopole cancels between block and kernel (∓12.8 vs ±12.8), so PW-R1 holds for ⟨TTTT⟩ at O(1/m²).**
>   * Exact O(1/m²) cost/H: ≤ 0.41 (J = 26), ≤ 0.05 (J ≥ 64), ≤ 0.005 (J ≥ 200).
>   * κ_B(+−) ≈ 3.8 matches the (1100, 100) B fits [NE].
> * **Flat ν⁶ tower:** subtracted exactly at O(y): the block reproduces the flat partial waves E₂d⁰⁰ and E₂d⁴⁴, and with the kernel's O(r) term gives exactly the hybrid's stripped flat part. At O(y²) the ++ case matches up to an unexplained +(3/16)J₀y²; +− at this order is [D-sketch].
> * **P4 remainder:** scalar O(1/m⁴) minus the flat tower gives K₄ ≲ 3 for J ≤ 50 (exceptions J = 46: 26.6; J = 26, 34 marginal), 10–16 at J = 64, and 20–7·10³ at J ≥ 80. The cause is the AdS/flat kinematic mismatch when costed to ν ≈ Λ ≈ m, where the fixed-β series is invalid. **P4 needs a ν-window statement instead of the truncated full-range series.**
> * **Theorem W [P given P0, P1′, P4] + [DV]:** J\* = 26 (K₄ = 0), 30 (K₄ = 1), 38 (K₄ = 4), 50 (K₄ = 10), 100 (K₄ = 35), at θ = 1. With the state-wise computed scalar remainder there is no J\*: only J = 30, 33, 38, 42, 50, 64 and 301 pass.
> * **Open:**
>   * P4 (a ν-window bound, plus the spinning O(1/m⁴) symbols and the kernel O(r²) term);
>   * the (3/16)J₀y² residual;
>   * P1′ beyond the leading-Regge kernel;
>   * the J = 0 engine comparison;
>   * behaviour between samples.

> **Update (P4_WINDOW.md, unreviewed): P4 is reformulated as a ν-window statement and reduced to three premises (W, AG′, EDGE). It is NOT proved: the beyond-all-orders part remains an assumption [A], now with explicit constants. Do not write "proved for J ≥ 26 / 42".**
> * **The full-range P4 failure is a series artefact [DV].** On exact scalar heavy moments, the true remainder costs ≤ 1.5·10⁻³ of the budget (K₄-units 0.99 at J = 64, 8.1 at J = 100). The fixed-β series is asymptotic and worsens past x ≈ 0.3–0.4. In the Regge window ν < 3/β it converges (successive ratio ≤ 6.5·10⁻³ at all 47 states) and matches the exact cost to 0.3 %.
> * **Lemma P4W [P] (Cauchy).** cost(G) ≤ window sup-bound + boundary term + an edge term at ν = Λ (suppressed by κ⁻⁶) + the exact monopole. Premises:
>   * (W) a real-axis window bound;
>   * (AG′) exact − flat is a pair of analytic travelling waves on [3/β, Λ], with local decay 0.97β/√(1−x²);
>   * (EDGE) an amplitude bound at ν = Λ.
>
>   The 1/m series is never used outside the window. Checked on the exact scalar: bound 4.0–7.9 % of budget.
> * **Theorem W, 47 states [DV]:** computed part ≤ 0.49 H everywhere (k = 3 included); every state where the full-range series failed now passes. With the scalar-calibrated edge, J\* = 26 (edge multiple tolerated up to 114; spinning O(1/m⁴) coefficient tolerated up to 20.8).
>   * **If the +− edge carries the Ω^U helicity factor (1−x₁²)⁻²**, then (26, 1.026Λ), (30, 1.030Λ), (34, 1.033Λ) and (38, 1.037Λ) fail. J\* becomes 42 (edge multiple 1), 51 (×3), 100 (×10).
>   * The remaining inequality there is |ρ₊₋(x₁)| ≤ 0.15–0.63; the exact scalar has 0.27–0.69.
> * **Theorem L′ range (16 states, J = 223 to 10⁶):** P4 is discharged with ε_R ≤ 0.086 of the certified half-floor (0.28 at Λ_amp = 10). EDGE is not binding (tolerance ≥ 600× the scalar envelope).
> * **(3/16)J₀y² residual resolved:** it is a genuine, harmless beyond-hybrid O(1/m⁴) term.
> * **Remaining:**
>   * (AG′) beyond all orders [A]: local rate [NE] at x ≤ 0.86, [M] at x₁ > 0.94;
>   * the graviton edge amplitude and its +− factor [M];
>   * the spinning O(1/m⁴) window coefficients [M] (margin ×5–10⁴);
>   * the boundary factor 2 [D-sketch];
>   * sampled states only; float only.
>
>   Sensitivity: with threshold 10/κ the edge multiple tolerated drops from 114 to 18.7.

> **Update (EXT_FOLD.md, unreviewed): the external-machine points are folded in. J = 7 and J = 27 → [NE±]; J = 23 and J = 31 are [NE±] at their higher masses but stay [NE+M] at the twist edge; J = 29 still does not close sign-free (L 0.94). No state came out negative.**
> * **J = 7:** twist edge +0.2683, L = 6.7; 1.01–1.05Λ: L = 7.2–25. **[NE±]** (its own direct row; shadow coefficient computed and checked bit-for-bit).
> * **J = 27:** twist edge +0.976 ± 0.019, L 0.28 → **1.64** (signed 2.02); 1.03–1.10Λ: L = 1.93–3.21. **[NE±].**
> * **J = 23:** twist edge +0.923, L 0.45 → 0.85, stays [NE+M]. It is now limited by the degree-7 shape-fit representation term, which Monte Carlo noise amplifies; a noise-aware term would give L ≈ 1.52, but that is a protocol change. 1.03–1.10Λ: L = 1.16–2.33, **[NE±]**.
> * **J = 31:** it has no GK15 edge run, so the tail beyond x = 0.9 must be covered sign-free; twist edge to 1.05Λ: L = 0.31–0.83, stays [NE+M]. 1.075/1.10Λ: L = 1.53/1.91, **[NE±]**. A GK15 edge window on [0.8, 0.970] is needed.
> * **J = 29 sign-free:** not closed. The 16 B seeds give only 23 independent streams: in n3, replica r of seed s uses stream s + r, so consecutive seeds share 7 of 8 replicas; independent seeds must be spaced by ≥ 8. The streams are heavy-tailed (|z| ≈ 11), and the binding inputs are the n1 node at ν = 20 and the 10⁴ convention below ν = 10. Even with σ → 0 on ν = 10–17.5, L = 0.97. **SMALLNU's "≈16× n3 points would close it" is false.**
> * **Checks:**
>   * External points reproduce local recomputation to 1.5·10⁻¹² / 4.5·10⁻⁷.
>   * Direct-vs-K₂ differences are always favourable (median 58–89 c₂).
>   * The K₂ error is 0.6–8 % of the last Laplace step over seven states (HIGHSPIN's "4–5 %" is corrected).
>   * Job A's seed is irrelevant: n1 uses fixed Sobol seeds.
>   * Job C (J = 400) failed: a float SVD failure in the conserved-structure solver at large J.
> * **Code issues:** hs_odd.py's non-GK edge lacks a sign-free tail term, and its degree-7 representation probe is noise-amplified on direct rows.

> **Update (EXT_FOLD addendum, batch 2):** the J = 31 GK15 edge points from the external machine replace the unmeasured edge tail. J = 31's twist edge goes +0.805 → +0.862 ± 0.023, L 0.31 → 0.70 (×3 window envelope) / 0.85 (signed); it stays [NE+M], limited by the noise-amplified representation probe (0.377) and the edge quadrature (0.213 + 0.167). **J = 31 at 1.05Λ → [NE±]** (L 1.22); 1.04Λ is borderline (1.00). Nothing negative.

> **Update (AG_PROOF.md, unreviewed): (AG′) is NOT proved, and in its literal local-rate form it is false near the gap edge. It is replaced by a strictly weaker pair of premises, (AG_U) + (EDGE-R), that suffices. With them Theorem W's sampled J\* = 26 (no longer 42 under the old premises), and Theorem L′ holds on the uniform rate alone.**
> * **Negatives:**
>   * **N1 [P + DV]:** the flat partial wave itself violates (AG′)'s pointwise bound near the edge (branch point at ν = m; ×4.1 excess at J = 26). The edge integral the lemma actually needs still holds (0.72–0.82 of (AG′)'s). P4_WINDOW's 15/κ-vs-10/κ sensitivity is an artefact.
>   * **G [NE]:** the exact scalar continued to complex ν (Arb, entire series) has an effective edge rate of 0.64–0.86 κ_loc at x₁ = 0.94–0.96, against the assumed 0.97. In the bulk it is 0.98–1.06. The uniform rate holds everywhere with margin (≥ 1.11κ₀; ≥ 2.3κ₀ at the edge).
>   * **N2 [P + DV]:** the plane-wave route cannot reach the near-gap edge. A positive average keeps only the worst rate, and at J ≤ 50, 14–35 % of the measure has |p| < Λ. It is harmless for J ≥ 100.
> * **Proved [P]:**
>   * **Lemma F:** an exact flat-wave split W₊ = (2/πi)C_J e^{i(J+1)θ}₂F₁(½, J+1; J+3/2; e^{2iθ}) with an explicit bound, so the uniform rate 0.97β holds for every explicit piece of the hybrid.
>   * **Lemma P4W-χ:** an entire erfc cutoff splits the complement cost into (EDGE-R), one real-axis integral E_χ[G] = ∫_{p_w}^1 ψ₂χG with no continuation and no rate, and (AG_U), a Cauchy part needing only the uniform rate 0.97Λβ. (AG′)+(EDGE) ⇒ (AG_U)+(EDGE-R).
> * **Constants [DV]:**
>   * The exact scalar edge integral is ≤ 0.30 S_E, 7–100× below (AG′)'s charge.
>   * Theorem W (47 states): everything except the edge integral uses ≤ 0.489 H. With the plain envelope the tolerated edge integral is ≥ 158 S_E, so **sampled J\* = 26**. With the Ω^U +− transfer the tolerance is ≥ 0.406 S_E,t at (26, 1.026Λ); J\* = 26 if the graviton +− edge is scalar-like, 38 if tolerance ≥ 1.2 is demanded.
>   * Theorem L′: the uniform rate alone suffices (ε_R ≤ 0.133; edge tolerated to 18.9× the transferred scalar envelope).
> * **Remaining premises:**
>   * **(AG_U) beyond all orders [A]:** [D] to all orders; near the edge it only asks for the split at rate κ₀.
>   * **(EDGE-R) for the graviton [M]:** decided by real-axis graviton values near the edge at J ≤ 51.
>   * The exact scalar data stop at x₁ = 0.962.
>   * Unchanged: P0, P1′, W, the boundary factor 2 [D-sketch], the +− transfer convention [M]; sampled states only, float only.

> **Update (DSKETCH_UP.md, unreviewed):**
> * **c₁: universality and ν-independence → [D].**
>   * The regular second-order ASD solution is explicit [DV] (exact sympy at two kinematic points), and exists for all momenta [D].
>   * The remnant is Φ₃₄ = −b₊W₊ with b₊ ≠ 0, which agrees with the JJ21 check.
>   * J₀(E×W³) = 2 and polynomiality are [D]; the Mellin → ν map T ↦ (A²+ν²)/2A is exact by Barnes' first lemma [P+DV].
>   * So "(Pos) ⇒ (1)" is [D] when no light single trace couples to TT.
> * **Light single traces: NEGATIVE finding.** In CHMRS's exact Mellin representation, the per-ν action of a single light block grows like e^{cν} (c → π/2) above the block's mass. **The leading-Regge/hybrid "image" of a light block is not its AdS image: the UT3 §5.5 row "LR form = true AdS image" is refuted as stated.**
>   * Its light value is the arc at infinity: a polynomial in ν² of degree n − k, [D given (E′)]. (E′) is CHMRS's light-side formula carried to the spinning k = 2, 3 functionals, and is [D-sketch].
>   * For J = 0 the AdS effect is m² → μ² = (Δ_L−2)(Δ_L−1) + δ.
>   * Ω^EE on light states is robust: J = 0 positive iff μ² > −0.545Λ².
>   * LS6 (ε = 0.040) is positive iff μ² > −692.5 / −7748.7, against a natural |δ| ≈ 21 [M]. ε\* is expected to fail for Δ_L ≲ 6, so 0.040 is necessary.
> * **Caution (unchecked):** weights reaching ν above a near-gap state's mass (e.g. de-shift tails to Λ+100) may meet the same exponential growth. This needs a re-check.

> **Update (EDGE_R.md, unreviewed): (EDGE-R) for the graviton is NOT closed; it stays [M] at all 30 Theorem W states with J ≤ 51. No violation was found, but no graviton E_χ was computed at any W state. External compute is the route.**
> * **What decides it:** at the gap-mass states the value comes from x ∈ [0.94, x₁]. For the exact scalar the parts above and below ν = 0.94Δ cancel to a net 0.11–0.29 S_E2, so sign-free covers cost 4–10× the tolerance and the oscillating values are needed.
> * **Why n1 fails at ν/Δ ≳ 0.94:** the spin prefactor grows at complex angle and puts a second hump on the real contour (e^{14} above the saddle). Re-centring cannot fix it.
> * **n3:**
>   * scalar: works to the edge (within 1.6σ of the exact Arb scalar up to ν/Δ = 0.970);
>   * graviton k = 2 (EE): agrees with n1 where both work, and gives values beyond 0.94;
>   * **graviton k = 3 (C3H): broken in n3** (state values ×200–1000 too small, cause not found).
> * **Exact scalar at J = 26–38 (new):** 0.285 / 0.201 / 0.108 / 0.151 S_E2 at the states' own cutoffs. **Correction to AG_PROOF: "exact scalar ≤ 0.30 S_E over all cutoffs" is wrong.** At J = 26 the grid maximum is 0.405, essentially at the tolerance of 0.406. At the tabulated cutoff a scalar-like graviton passes with margin 1.42.
> * **Units correction:** the plain tolerance (≥ 158 S_E) and the transferred tolerance (≥ 0.406 S_E,t) are the same number. AG_PROOF's "plain margin ≥ 520" compares with the scalar edge, not the graviton.
> * **External route:** ext/jobs_batch3.txt (1586 n1 + 150 n3 jobs, about 90 CPU-h); ext/jobs_batch3_n1.txt holds the n1 part only. The 20 states with m ≥ 1.13Λ need only n1 and could reach [NE±]. The gap-mass states also need the n3 C3H fix and 2–4× less bulk MC noise.

> **Update (NU_TAIL_CHECK.md, unreviewed): the exponential growth above a state's mass affects NO computed number, but the published grid-weight DEFINITION at Δ_g = 1000 is defective, and the fix is free.**
> * **No effect on computed numbers [P, by inspection]:** every computed per-state number (crystal/D′, K₂ rows, direct n1/n3, GK edge nodes, windows, EDGE-R, hybrid floors) uses ν ≤ Λ ≤ m. The light-side constants are unaffected, and so is everything built from the continuum weights or the exact min-norm X: X and X^Q decay like e^{−(π/2)(ν−Λ)} beyond Λ, and tails beyond m are ≤ 4.4·10⁻⁷ (D) and ≤ 1.6·10⁻⁴ (Q) of a state's action.
> * **Problem (c):** the stored grids `prl/data/d_abel_X_1000.txt` and `prl/data/k2/xq_abel_X_1000.txt` run to Λ+312 and carry a non-decaying double-precision floor (|X| ≈ 10⁻¹⁵–10⁻¹⁶) beyond Λ. Read literally with θ(Λ+100, 10), that floor times the exact image growth (Γ ≈ 10¹⁶ at Λ+100) gives tails up to 650× (D) and 1300× (Q) a near-gap state's action. The sign is set by round-off, so (Pos) for that grid-defined object is untested at m − Λ ≲ 60–100. HK_CLOSE §5 and §6 claims are false for that object, though true for the computed values. Δ_g = 300 grids are unaffected (they end at Λ+29).
> * **Fix (changes no number):** define the functional by the continuum weights; or read the grids as zero beyond Λ+10; or use the cutoff (Λ+60, 6) for the k = 2 pieces (same leakage, tail 5·10⁻⁹ / 4.6·10⁻⁷).
> * **Growth numbers [DV]:** exact scalar Arb series. Γ(m) = 1.8·10³ and Γ(Λ+100) = 1.6·10¹⁶ at Δ_g = 1000, J = 0; the rate per unit ν is 0.13–0.19 at x = 1, rising to π/2 only as x → ∞.

> **Update (GAP_SCALING.md, unreviewed): extending positivity to all Δ_g ≥ 1000 by scaling is PARTIALLY viable. The certified layers carry to every larger gap; the beyond-leading-Regge corrections carry over only under tested scaling hypotheses.**
> * **Weights:** the flat weights ψ_k are exactly Λ-independent in gap units. The k=2 de-shift weight is ψ₂ + Λ⁻²·(explicit symbol) for p ≥ 0.1 [DV]. The leading-Regge action is exactly ∫ψP at every Λ.
> * **Certified layers:**
>   * The curvature is g(J, m/Λ)/Λ² [DV].
>   * The **Lemma C bound is monotone decreasing in Λ [P]**, so "certified flat floor − Lemma C" only improves with the gap; it is positive near the twist edge up to J = 13 / 17 / 23 / 35 / 61 at Δ_g = 1000 / 2000 / 5000 / 2·10⁴ / 10⁵.
>   * **Lemma LR∞ is certified uniformly for 998.5 ≤ Λ ≤ 10⁶** (454 Arb Λ-balls, remainder ratio ≤ 0.9988), so the Theorem L′ floor and J\*(m) hold for every gap in that range [certified].
> * **Binding state at Δ_g = 2000 [NE±]:** own-mass row Δ = 2005 gives λ/N_s = +0.20041, budget 6.3·10⁻⁴, L = 320 (it was 78 at 1000). The truncation scales like Λ⁻² and the window like Λ⁻⁴; the graviton correction profile converges (|Γ| ≤ 41 c₂).
> * **What does not carry over automatically:**
>   * Sign-free budgets at fixed β do not improve (pointwise corrections do not shrink).
>   * EDGE-R's +− transfer factor grows like (Δ_g/J)², so at fixed spin it gets harder.
>   * The graviton edge integral is uncomputed at every gap.
>   * There is a seam strip between the 1000 and new twist edges with no 1000 counterpart (J up to ≈ 344 / 424 / 552 at Δ_g = 2000 / 3000 / 5000).
>   * The structural premises (AG_U, EDGE-R, P1′, W) carry over as premises.
> * **Plan (3–5 days):**
>   1. Write up Lemma C monotonicity and uniform LR∞.
>   2. Own-mass K₂ rows at Δ_g = 2000, 5000 for every error-budgeted spin.
>   3. Direct n1/n3 nodes at 2000 for J = 12, 16, 21, 27, 29 (external machine).
>   4. Parameterise Δ_g in the large-spin code and re-run the P2/P4W-χ costs at the chart-B images.
>
>   Target: error-controlled positivity at Δ_g = 1000, 2000 and 5000, plus a transfer to all Δ_g ≥ 1000 under named, tested scaling hypotheses.

> **Update (GAP_LARGESPIN.md, unreviewed): the large-spin floor is now certified for EVERY Δ_g ≥ 1000 (no upper limit), and the large-spin theorems' margins grow with the gap.**
> * **Lemma LR∞ for all Δ_g ≥ 1000 [P + certified]:** a scaling lemma makes Λ⁴·(each bound piece) non-increasing in Λ (envelope R̂ = min(R_∞, V/s⁴), V₂ ≤ 684625, V₃ ≤ 789636), plus a 13-step monotone Λ-grid below 2604.28. Thus r < 1 on β ≥ 0.37 and F_LR > 0 on β ≥ 0.245; J\*(m) = ⌈(0.37211m−1)/2⌉, with J\*/Λ → 0.22387 as Δ_g → ∞. This closes the tail above Λ = 10⁶.
>   * Thresholds improve with the gap: β₀ = 0.2375 / 0.1675 / 0.1175 for Δ_g ≥ 2606 / 5000 / 10⁴.
>   * Gap-unit floor: Λ²F_LR ≥ f_U(β) for all Δ_g ≥ 1000 [certified].
> * **ROUTEA curvature window uniform in the gap [D + certified]:** curv > 0 on b̂ ∈ [114.06, 2Λ] for m ≥ Λ (new exact identity ∫₀^{s₁}g₁ = −s₁³r(s₁)/24). **Lemma C′ is monotone in Λ [P]**, so b̂ ≤ 119 carries to every gap. **For m ≥ 1.3Λ the hybrid floor is rigorous on every b̂ ≤ 2Λ at every gap**, which closes the old (119, 245) hole for m in [1.3, 1.5)Λ.
> * **Theorems L′ / W at Δ_g = 1000 / 2000 / 5000 [DV, premises carried with Δ_g-independent AdS-unit constants]:**
>   * Λ²D = −995.3 / −999.1 / −1001.3.
>   * L′ remainder: max d/(F/2) = 0.19 / 0.035 / 0.005; edge tolerance 12.6 / 228 / 9050.
>   * **Theorem W: d_χ/H = 0.489 / 0.246 / 0.100; the tolerated graviton edge integral at the binding state is 0.406 / 3.45 / 57 S_E,t (∝ Λ^3.1).** So the EDGE-R premise gets much easier with the gap.
>   * Caveat: the 5000 C2/vector tolerance entries are not numerically converged (finite-difference β-derivative).
> * **Remaining:** the premises (P_L′ magnitudes, AG_U, EDGE-R, P1′, W) carry over only under the hypothesis of Δ_g-independent AdS-unit constants. The trend between computed gaps is [D-sketch]. The small-J twist-edge seam and β₀ < 0.10 are not addressed.

> **Update (GAP_UNIFORM.md, unreviewed): there is no uniform-in-gap proof for the near-gap corrections, and the literature has none either. The needed inequality is pinned down, and J ≤ 7 is positive at every gap given one explicit "no-growth" premise with large safety factors.**
> * **Literature [Q]:** no explicit uniform-in-Δ_gap bound exists.
>   * 2205.01495: the uplift is one sentence, with an unquantified O(1/Δ_gap⁴) in eq. (30).
>   * CHMRS §5: non-effective ("sufficiently large M", (5.7), (5.11), (5.26)–(5.29)).
>   * Knop–Mazáč AdS₂: their uniformity comes from pointwise kernel positivity against dDisc ((3.25)–(4.9)), which does not transfer to Ω^EE.
>   * CLS-D: leading order only.
>   * van Rees–Zhao 2210.15683: uniformity "not immediate".
> * **Lemma T [P]:** in chart A the node coefficients, the flat parts and N_f are exactly Λ-independent and the Lemma-C bound decreases, so positivity at a larger gap needs only that Γ/m² does not grow. Checked on 366 states at 1000/2000/5000.
> * **Prop. S [D]:** at fixed ξ = ν/Δ the crystal integral is ∫e^{ΔΦ_ξ}g_J with Φ_ξ exactly Δ-independent, so H_mono-A is ordinary fixed-phase Laplace asymptotics.
> * **Lemma R [D-sketch]:** the D′ thimble coefficients are Δ-independent in chart A (removes the finite-Δ Stokes worry). **Theorem E [D-sketch]:** each fixed near-gap state is positive for all Λ ≥ Λ\*(J, μ). Prop. D′ controls only ξ < 0.107.
> * **Coverage:**
>   * **J ∈ {0, 2, 4, 5, 6, 7} positive at every Δ_g ≥ 1000** [P given B_F + windows]. Premise B_F: sup|Γ_Δg| ≤ F_J·(sup|Γ| at 1000), with F_J = 11.8 (J=7), 14.9 (6), 43 (4), 743 (5), 799 (2), 3.4·10⁴ (0). B_F is [NE]: graviton 1000 → 2000 ratio 0.97–1.06; exact scalar 1000 → 10⁴ drifts −2.6 % to +1.2 %.
>   * J = 8 fails with the certified floor (F = 0.85).
>   * J ≥ 12: the sign-free form fails already at 1000.
> * **Main missing step:** explicit Laplace remainder constants C₃(ξ, J). They need interval sup-norm bounds of the spinning amplitudes on a fixed complex polydisc (checked once at Δ₀, valid for all Δ ≥ Δ₀), which needs a pointwise complex evaluator the project lacks.

> **Update (GAP_C3.md, unreviewed): the Laplace remainder constant C₃ is explicit and certified (Arb Taylor-model evaluator on a complex polydisc around the crystal; Lemma L3 [P]), but it is ~10⁸–10¹² too large to use at Δ_g = 1000. It gives a certified threshold Δ₁ ≈ 5·10⁵–10⁷ above which K₂ truncation is accurate enough. The premise B_F is still needed on 1000 ≤ Δ_g ≲ 10⁶.**
> * **Certified:** the scalar k=2 integrand with Q_J at J ∈ {0, 2, 4–8} × ξ ∈ {0.1 … 0.85} (35 nodes), giving Δ₁ = 4.7·10⁵ … 10⁷; graviton EE and C3H at J = 0, 2, giving Δ₁ = 1.0–1.8·10⁶.
> * **Eventual positivity:** GAP_UNIFORM's Λ\* is now explicit, Δ_g\* ≲ Δ₁, at the certified nodes. This holds given:
>   * Lemma R (contour) [D-sketch];
>   * the chain propagation factor Π ≤ 3 [M];
>   * the K₂ part of Γ for all m [NE];
>   * continuity in ξ [D-sketch].
> * **The obstruction is the method:** the genuine cubic size of Φ_ξ (S₃ ≈ 10–12) caps the radius at ~0.04 against a Gaussian width of 0.032 at Δ = 1000. Ideal sup norms give the same Δ₁. The true next Laplace term is ~2·10⁻⁵ relative at Δ = 1000, so K₂ is accurate and the certificate is pessimistic by ~10¹⁰.
> * **Next lever:** a certified holomorphic Morse chart (rough estimate Δ₁ ≈ 3·10³–1.3·10⁴), ξ as a Taylor-model variable, and the J ≥ 4 two-structure / odd-J graviton amplitudes.

> **Update (GAP_NEARGAP.md, unreviewed): near-gap states at Δ_g = 2000 (and partly 5000), on own-mass K₂ rows. Nothing came out negative. Chart A (same J, same m/Λ) scaling holds exactly for J ≤ 16, and leverages grow 4–9× at 2000.**
> * **[NE±] at 2000:**
>   * light J = 0, 2: L = 40612 / 5712.
>   * J = 4: L = 106 (twist edge 103).
>   * J = 5 (binding): +0.2004, L = 320–334; at 5000, L = 1922–2066.
>   * J = 6: L = 22.
>   * J = 7: L = 27.5.
>   * **J = 8 sign-free: L = 5.2** (was λ-level only at 1000).
>   * J = 12, 16: λ-level L = 8.0 / 4.4 (J = 12 at 5000: sign-free 6.4).
> * **[NE+M] at the 2000 twist edges for J = 21, 27, 29:** λ-level L = 0.4–1.1. This is an artefact: the Δ_ref = 2·10⁴ reference's partial-wave zero falls inside the 2000 rows for J ≥ 16. Use Δ_ref ≥ 5·10⁴ at 2000 and ≳ 1.3·10⁵ at 5000.
> * **H_mono-A confirmed numerically:** Γ(x) is Λ-independent at fixed (J, x) to 30 c₂ (J ≤ 6) and ≤ 330 c₂ (J ≤ 29) from 1000 to 2000, and 9 / 248 c₂ (J = 5 / 12) from 1000 to 5000. Budget terms scale by 0.25 ± 0.01 at 2000 and 0.040–0.043 at 5000, as Λ⁻² predicts.
> * **Batch 4** (`ext/jobs_batch4.txt`): 345 n1 jobs (247 main for J = 12, 16, 21, 27, 29 at 2000 + J; 98 optional seam jobs). About 1–2 h on 16 workers; pilot validated.
> * **Still open:** seam spins J = 9–11, 13–15, 17–20, 22–26, 28, 30–32; the parity-odd OO states J = 4–10 at 2000.

> **Update (GAP_MORSE.md, unreviewed): a certified holomorphic Morse chart plus Laplace remainder (Lemma M [P], data [certified]) lowers the certified threshold Δ₁ by ×16–47 relative to GAP_C3, to Δ₁\*(J) ≈ 0.9·10⁵ (J = 0), 1.1·10⁵ (2), 1.8·10⁵ (4), 1.3·10⁵ (5), 2.2·10⁵ (6), 2.5·10⁵ (7), 4.2·10⁵ (8), and graviton EE/C3H at J = 0, 2: 2.2/3.0·10⁵. It still does not reach 1000–5000.**
> * **Gap coverage:** 1000 ≤ Δ_g < Δ₁\* (about 2–2.6 decades) rests on premise B_F, checked at the [NE±] rows at 1000, 2000 (and 5000) and by the exact scalar to 10⁴. Above Δ₁\* positivity is [P given GAP_C3's premises]: the TAIL, Π ≤ 3, the K₂ part bounded, and ξ-continuity.
> * **Interval coverage:** ξ-intervals work (one certificate per interval, ×4–8 cost); full [0.1, 0.85] coverage is ~2 CPU-h per spin and was not run. ξ = 0.1 sets the threshold.
> * **Why 1000 is out of reach:** Δ₁ ∝ r_c⁻², and Δ = 1000 needs r_c ≈ 0.24, against a certified radius of 0.023–0.052.
> * **Defect:** GAP_C3's `gc_ball.log` used float series coefficients outside the enclosure. It is numerically irrelevant, and fixed in GAP_MORSE's copy.

> **Update (GAP_SEAM.md, unreviewed): every state run at Δ_g = 2000 is positive.**
> * **Sign-free [NE±] at 2000** (Δ_ref = 5·10⁴): seam J = 9 (L = 8–9), 10 (1.8–2.0), 11 (2.6), 13 (1.03, thin). Parity-odd OO J = 4, 6, 8, 10 (L = 5.7–167); the OO value comes from Γ_OO = Γ_++, checked to ≤ 0.044 c₂.
> * **λ-level only at 2000:** J = 14, 15, 16, 17, 19 and the J = 21 rerun (λ-level L ≥ 1.7–5.2); sign-free fails for J ≥ 14.
> * **Δ_ref artefact:** confirmed for the λ-level truncation and fixed by Δ_ref = 5·10⁴ (J = 21's λ-step falls ×2.5–2.7, toward the direct values). It is not the cause of the J ≥ 17 sign-free inflation, which comes from the extra x = 0.025 node in the a-fit; use the same x grid at both gaps.
> * **Chart A holds:** budget ratios 0.23–0.29 (predicted 0.25), small-ν term 0.059–0.065; leverages grow ×4–14.
> * **Cautionary:**
>   * The seam strip holds the lowest values of J = 11, 13, 15, 17 at 2000 (+0.30 to +0.53; hybrid dips, still ≥ ×1.5 the binding state). Chart A cannot supply these states.
>   * At 1000 the OO J = 4 twist edge (+0.2059) is below the previously quoted minimum (+0.2150).
>   * GAP_NEARGAP's printed J = 21 value (+0.967) should read +1.0009.
> * **Remaining:**
>   * J = 18, 20, 22–26, 28, 30–32 not run;
>   * J = 27 and 29 should be rerun with Δ_ref = 5·10⁴;
>   * sign-free closure for J ≥ 14 needs direct rows (batch 4);
>   * there is no direct check at 2000.

> **Update (GAP_SEAM2.md, unreviewed): the Δ_g = 2000 reruns of J = 27, 29 with Δ_ref = 5·10⁴ did not lift their twist edges. They stay [NE+M]; no state is negative (all ≥ +0.91).**
> * **Why the expected doubling failed:** the λ-level truncation step fell as predicted (chart-A ratio 0.22–0.25), but two other terms took over. The edge term (the row stops at 0.8 while |Γ| ≈ 10⁴ c₂ at x = 0.85–0.9) and the small-ν window (larger when the first node is 0.05) dominate.
> * **Chart A holds cleanly** once the node set and the reference are matched: E_trunc ratio 0.26–0.27 and λ-step ratio 0.22–0.25, against 0.25 predicted. Γ(x) at 2027 matches 1027 to ≤ 40 c₂ for x ≤ 0.6.
> * **New lower state:** J = 27 at 1.0175Λ, +0.910, inside the seam strip. J = 29 at 1.075/1.10Λ and J = 27 at 1.075/1.10Λ are λ-level [NE±].
> * **J = 27, 29 twist edges at 2000 need the batch-4 direct rows (or a K₃ order).** J = 18, 20, 22, 24 not done; resume commands are in GAP_SEAM2 §3. Hybrid minima for J = 22/24 at 2000: +0.780 / +0.826.

> **Update (NEARGAP_FILL.md, unreviewed): at Δ_g = 1000 every near-gap spin J ≤ 16 is now [NE±] (twist edge and 1.02–1.05Λ). J = 23's twist edge → [NE±], thin (L = 1.18). J = 31's twist edge stays [NE+M]. Nothing negative.**
> * **Method:** full direct n1 rows (38 nodes) on own-mass rows plus the NEARGAP2 budget and a 3σ MC term. J = 9: L = 38; 10: 5.6; 11: 16; 13: 6.8; 14: 2.07; 15: 1.75 at the twist edge (thinnest 1.52, J = 14 at 1.05Λ). The K₂ error is 0.4–7 % of the last Laplace step, always favourable. The Δ_ref zero lies outside every row for J ≤ 16.
> * **J = 23/31:** a consistency correction (the representation term evaluated on the same windowed λ as the leverage; the window had been counted twice) lifts J = 23's twist edge to L = 1.18 and J = 27 from 1.64 to 1.81. J = 31's twist edge stays at L = 0.76: genuine half-grid under-resolution plus window and edge terms. J = 31 at 1.04Λ → [NE±].
> * **Batch 5** (`ext/jobs_batch5.txt`): 666 n1 jobs for J = 17–20, 22, 24–26, 28, 30, 32 plus J = 23/31 mid-points; ≈ 1.8 h on 16 workers. J = 26 and 30 need even shadows computed locally first.
> * **Cautions:** for J ≥ 13 the base is the hybrid quadrature (the certified floor does not close). The 10⁴ c₂ small-ν convention (~100× the measured values) is still the largest or second-largest budget term at J ≤ 14.

> **Update (EDGE_R2.md, unreviewed): the batch-3 n1 rows do not close (EDGE-R) at any of the 8 gap-mass Theorem W states. The ++ entry closes at (26, 1.026) and (30, 1.030) [NE+M]; the binding +− entry stays open everywhere. Theorem W's J\* = 26 is unchanged and still conditional on (EDGE-R).**
> * **Data-handling finding:** **even-J records from the external machine are in a different ⟨TTO⟩ conserved-structure basis** (an SVD null-space basis, machine-dependent, differing by a 2×2 transformation). Repaired by joint diagonalisation (`code/edger2/e2_basis.py`); the EE fit residual is 4e-7–3e-4 and C3H comes out diagonal unfitted. **Every even-J external job needs this repair**, or the external runs should ship their basis. Odd J is unaffected.
> * **Warning (not an established violation):** the ++ central values at (46, 1.045) and (50, 1.049) are 2.0× / 2.6× the EDGE-R tolerance. Their 3σ budgets reach below it, and the bulk is aliased on the 0.025 grid, so they are [NE] at best. Recheck on a denser grid (~20 n1 nodes per state). This concerns the premise's tolerance, not the sign of the state's functional value.
> * **To close +−:** the n3 k=2 edge jobs, the n3 C3H fix (the k=3 edge remainder alone exceeds tolerance at all 8 states), lower MC, and a 0.0125 bulk grid. A combined-weight version of Lemma P4W-χ that keeps the C2/C3H cancellation (J = 26: D₂+D₃ = 0.05 vs 0.27) is a possible route [M].
> * **Code issues:** `er_jobs.state_nodes` puts the edge-window start one shift late; the `ng_val._clean` tolerance fails after hs_core re-keying (J = 12/16 checked unaffected).

> **Update (BATCH5_FOLD.md, unreviewed): batch 5 is folded in at Δ_g = 1000. 25 of 31 near-gap spins J ≤ 32 are now error-controlled at the twist edge and every computed mass up to 1.05Λ (was 19). Nothing negative.**
> * **New [NE±]:** J = 17 (L 1.53), 18 (1.83), 19 (1.84), 20 (1.41), 25 (1.65); **J = 23 no longer thin** (L 1.62; the new mid-points drop the representation term from 0.258 to 0.020); J = 31 on the signed window only (1.14, like J = 29).
> * **Twist edges still [NE+M]:** even J = 22 (0.73), 24 (0.76), 26 (0.95), 28 (0.66), 30 (0.47), 32 (0.41), limited by the flat sign-free edge term beyond the last valid node plus the ×3 window. With the odd-J edge envelope J = 22–28 would reach L ≈ 1.0–1.4 (a protocol change; diagnostic only).
> * **Basis repair** of even-J external records: per-row rotation/reflection solved exactly, fit residual 2e-12–2e-5, held-out nodes ≤ 7e-5, Γ vs local recomputation ≤ 0.34 c₂. (EDGE_R2's e2_basis fails for J = 18, 20, 24, 28 because of a per-helicity sign.)
> * **Cautions:**
>   * J = 20 drops +0.965 → +0.835 (it matches HIGHSPIN's +0.83 ± 0.14 and still closes); J = 28 drops 0.04.
>   * The hs_odd shape fit is ill-conditioned at J = 17 and 19 (folded with the own-row method instead).
>   * `ng_bind` was despiking real injected direct edge nodes (exempted; values move ≤ 0.001).
>   * Held states (J = 17, 19, 22, 26, 30 above the twist) rest on the own-row hold.
> * **Count variants:** 25 (as stated); 23 without the signed-window-only J = 29, 31; 18 excluding held masses; 16 if every mass must close sign-free with the ×3 envelope.

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

> **Update (INDEP_CHECK3.md, independent): c₁ = 0 partly reproduced and not refuted; batch 5 reproduced.**
> * **c₁ steps:**
>   * (i) the affine k=2 weight is reproduced [P], with the finite-r error ∝ 1/r.
>   * (ii) the first-moment helicity pattern is reproduced with an independent 30-digit free boson/fermion engine [D+DV]: the momentum Ward identity holds to 3e-30, the flip part of M is exactly 0, and M = −½{Γ,H} to 3e-25. The parity-odd α_T Ward argument re-derives cleanly but is [D] only.
>   * (iii) the EE contraction is helicity-diagonal [P]; the d=3 orientation reversal is checked harmless.
>   * The Route-B zero is reproduced to 6.9e-24, with controls O(0.1–4).
> * **c₁ grade refinement:** "[D] given the [TT] double-zero statement (parity-odd odd-twist [TT] double zeros [D-sketch], H0 §15.1)", the same caveat as the Einstein term. Route B is not independent of Route A.
> * **Batch 5:** the basis repair is re-derived independently (J = 18, 26). The unfitted node agrees to 2.6e-7 against MC 3.6e-4, and the mapped records match to 1e-11. Values, budgets and L for J = 18, 23, 25 are bit-identical on re-run. "25 of 31" is confirmed, and so are the stricter counts 23 / 18 / 16. Minor: BATCH5 §1.1's failed-node list omits J = 18 x = 0.9447 and J = 20 x = 0.9432 (the chain drops them correctly).

> **Update (K_CERT.md, unreviewed): δK₀, the Einstein light value and τ₀(K) are certified (Arb); the H17 remainder is NOT proved (stays [D-sketch + M]), so there is no certified K_min yet. The 2.89×10⁻⁸ headline needs only a weak lower bound on H17.**
> * **Certified:**
>   * K₀ = 1 + 3.16e-11.
>   * δK₀ = −7.666324482972e-4 (300) and −6.830149983909e-5 (1000).
>   * K_c = 0.999233 / 0.999932.
>   * τ₀(K_c) = 3.133256e-6 / 2.891361e-8.
>   * The Einstein light value of the actual weight is 16L_A(1 ± 1e-9).
>   * At 2000 / 5000: τ₀(K_c) = 1.928332e-9 / 5.326183e-11.
> * **Conditions on H17 for the headline:**
>   * at 1000, τ₀ ≤ 2.90e-8 needs H17 ≥ −5.95e-3 (14× the model band), and |ĝ₃|Λ⁴ ≤ 120 needs H17 ≥ −4.1e-3 (9.7×);
>   * at 300, ≤ 3.15e-6 needs 2.2×, and ≤ 3.14e-6 is not robust.
>
>   The quartic upper end does not involve K.
> * **H17 reduced:** Λ²·H17 = −0.8106·a₄ + 4.6867·b₂ exactly (certified weights), where a₄, b₂ are the |ĝ₃|² light-function coefficients of the two functional pieces. b₂ (the Ĉ₃ piece) carries 79 % of the band. The headline ≤ 2.90e-8 holds while |b₂| ≤ 1247. **Computing b₂ (and a₄) gives a certified K_min directly.**
> * δK₀ is convention-dependent: it is 0 in the Ĉ′-primary convention, so the tables cover both ends.

> **Update (K_B2A4.md, unreviewed): H17 computed in closed form, so the model band is no longer needed. b₂ = −35/2, a₄ = −581/4 (leading), Λ²H17 = +35.7 at every gap [D]; the arithmetic to K and τ₀ is [certified]. K = 1 − 32.38/Λ².**
> * **Method:** at x₀ = ∞ the light-ray factorisation (C1_FINAL F₂, H7 F) plus a Fourier transform to the s-channel rest frame makes each piece's T-block light value an exact celestial-circle integral, evaluated in closed form. No Regge expansion and no Witten diagram are used.
>   * k = 3: −2πG|ĝ₃|²p⁴[1 − (35/2)ν⁻² + (569/16)ν⁻⁴] exactly.
>   * k = 2: rational, with poles at the zeros of R₋+R₊ that X_θ cancels; H17 is integrated from the exact function.
> * **Checks:** the same reduction reproduces exactly the Lemma-R Einstein function, including 16ν⁴ − 440ν² + 273. The flat top coefficient is exactly 1/(8·240²) (= K₀). Kernels match direct quadrature to 1e-10–1e-12. The value agrees with K_CERT's linear form to 3e-3. The physical shifted Einstein value settles δK₀ as physical.
> * **Certified K and headline:**
>   * K = 0.99963650 (300), 0.99996752 (1000), 0.99999189 (2000), 0.99999870 (5000).
>   * **τ₀ = 2.89131·10⁻⁸ (≤ 2.8914·10⁻⁸) at 1000** and 3.13262·10⁻⁶ at 300.
>   * τ₀ ≤ 2.90·10⁻⁸ holds with margin +5966/Λ² in K, and |ĝ₃|Λ⁴ ≤ 120 holds.
> * **Caveats:**
>   * Inherited [D] links: light-ray factorisation, the [TT] double zeros (parity-odd [D-sketch]), and Ω_OPE = Lemma R.
>   * k = 3 has two steps without end-to-end checks: the j′ = 1 structure sign, and N̂₃ = the J = 0 heavy action.
>   * The k = 2 α_T² part is [D-sketch].
>   * A factor-10 error in b₂ would still keep the headline.

> **Update (UNCOND_INVENTORY.md, analysis only): no finite Δ\* makes "Ω^EE ⪰ 0 on every heavy state" unconditional with known methods. The fully rigorous reachable layer is T0: the leading-Regge (hybrid) action ≥ c·N_s on every heavy state at Δ\* = 1000 (≈ 14–24 agent-hours, ≈ 85 %).**
> * **Obstructions to the full action:**
>   * (O1) a contour theorem (TAIL / Lemma R) needed by every Laplace certificate;
>   * (O2) AG_U beyond all orders at large b̂, where the value comes from cancellation;
>   * (O3) uniformity across J ~ √Δ_g;
>   * (O4) the light-cone fold, removable by redefining the weights to support p ≤ p_max < 1.
> * **Targets:**
>   * T1: the full action for J ≤ 8 at Δ\* ≈ 10⁶, given CT; 35–60 h; 30–40 % given CT, with CT itself ≲ 20 %.
>   * T∞: every heavy state at Δ\* ≈ 10⁶–10⁷, given CT; 120–200 h; ≲ 20–25 %.
> * **Infeasible:** any certified statement about the beyond-leading-Regge correction δ below Δ_g ≈ 2·10⁴. The optimal-truncation accuracy at Δ = 1000 (7e-3–1e-5) is comparable to the tolerances (3.5e-2–5.6e-5), and no exact non-asymptotic spinning representation exists. So 1000 ≤ Δ_g < 10⁵–10⁶ stays at the B_F premise. Below Δ_g ≈ 140–300 positivity is expected to be false (DG_SCAN model).
> * **Plan:**
>   * Phase 1: T0 (M1 + a continuum functional at every Δ_g, M20) in parallel with a capped CT probe (8–10 h; go/no-go).
>   * Phase 2 (T1) only on go.
>   * Cheap extras: the M19 implication-side [D-sketch] links and the M18 A1′ wording.

> **Update (LOWNU_SCOPE.md, ANGULAR_ROUTE.md, 25 Nov; unreviewed): low-ν unconditional functional is partially viable.**
> * **LF1/LF2 (flat):** the dilation identity is exact [P]. Dilated NF02 is certified flat-positive (margin 3/20) for a ≤ 0.612. A new LP that only needs positivity for u′ ≤ a² gains about 19× and passes a dense check [NE], but is not certified. τ at Δ_g=1000 with K(L): LP 1.56e-4 (a=0.05), 1.11e-5 (a=0.1); dilated NF02 2.9e-3 / 2.1e-4.
> * **LF3 (AdS remainder): blocked on the contour route.** Prop. D′ is a far-region absolute bound only. The near-region constant is not explicit. The usable ξ_max is ≈ 0.088 at Δ=1000 (≤ 0.118 from polydisc analyticity). The prefactor grows with spin, so J-uniformity fails.
> * **LF4:** K(L) is exact via the Q₆ identity. The e^{πν}-weighted tail forces L ≲ 0.05–0.06Λ, giving a best estimate of τ ≈ 0.8–1.6e-4. 10⁻⁶ needs Δ_g ≈ 3000.
> * **Angular route (ANGULAR_ROUTE.md):** the celestial kernels are continuous dual Hahn polynomials (verified exactly), so completely monotone smearing is positive on every harmonic. The decisive check is whether cubic sensitivity survives (the LP's sensitivity comes from a ψ₃ with opposite helicity signs). Heavy positivity remains conditional unless that check or the plane-wave image check succeeds (estimate < 25% by Thursday).

> **Update (T0_CERT.md, 25 Nov; unreviewed): T0 proved and certified. The leading-Regge (hybrid) action of Ω^EE is ≥ 0.1999·N_s on every heavy state for Δ_g ≥ 1000 (≥ 0.199·N_s for Δ_g ≥ 300).**
> * **Coverage:** ++ for all even J, +− for all J ≥ 4 (odd J included), every m ≥ Λ, every b̂ > 0, every real Δ_g ≥ 300. The 2×2 matrix has λ_min ≥ c·N_s for even J. No physics input is used; OO is [certified given helicity lemma (A)].
> * **Method:** an algebraic reduction to four sign-pattern conditions, each "flat margin + signed curvature ≥ 0". A new curvature lemma [P] in gap units replaces ROUTEA's [D] derivation. The certificate has three parts, using Arb/exact Descartes with proved inter-cell bounds:
>   * curvature > 0 on 16 ≤ b̂ ≤ 2Λ, with min 4.92/Λ²;
>   * b̂ > 2Λ via LR∞ plus a new flat-Regge remainder lemma [P] (≤ 4.6% of the floor);
>   * b̂ ≤ 16, where the flat action exceeds the possible negative curvature T = 2.9e-5.
>
>   Reproduces in about 20 CPU-min (`code/t0cert/run_all.sh`).
> * **Open:** M20, i.e. the realised published weights. For those there is only a corollary at Δ_g = 300/1000, resting on ε_W [DV] and limited to β ≲ 200. The flat certificate's corner code is unreviewed (inherited). The beyond-leading-Regge correction δ is untouched.

> **Update (ANGULAR_GO.md, 25 Nov; unreviewed): the angular route is a GO, conditional and with a ceiling τ ∝ 1/Δ_g.**
> * **Sensitivity.** Cubic sensitivity survives through the fourth-harmonic k=3 sum rule: \|ĝ₃\|² sits in heavy slot \|m\|=4. Einstein sits at harmonic 0, which no heavy state populates. The k=3 functional alone is forced to zero.
> * **Sign and light side.** The global sign is fixed physically: heavy must sit opposite Einstein. 303/5120 is confirmed against K_B2A4.
> * **Scalar first moment.** M_mid/H = i(Δ²−7)/6 is reproduced in both parities. The literal monomial rule is the retarded-time u = α/2 moment. The normalisation step to iT_Δ/6 is [D-sketch].
> * **Ceiling.** Δ-free harmonic positivity cannot depend on the gap. The gap enters only through T_Δ, giving τ ≈ 22/Δ_g (0.0224 at 1000). Bandwidth is capped at ν = O(1), because the CDH kernels alternate beyond ν ≈ 1.9m: the authors' weight works, stretched by 1.5 it fails. The route needs no LF3 remainder and no de-shift tail bound.
> * **Open.**
>   * The all-J tridiagonal identity (~1 day).
>   * The b=0, 4 ladder ratios (pipeline data to J≈30, then WZ proof; ~40% by Thursday).
>   * Physical matching: the rest-frame reduction for spinning heavy states, k=3 flip entries at the |m|=4 projection, the shared |c|² normalisation in 47266, swap and convergence, double-trace and contact terms, and the parity-odd α_T piece.

> **Update (INDEP_CHECK_T0.md, 25 Nov; independent): T0 confirmed, with minor documentary fixes. No number, constant or [certified] grade changes.**
> * **Reproduction:** byte-identical outputs from all 9 logs.
> * **Algebra:** the reduction is verified. Part F and the hybrid share the same ψ, and the flat part agrees to 5e-11.
> * **Curvature lemma:** re-derived. An independent ₂F₁ test gives h_eff ≤ 191, against the certified H_k ≈ 1.2e4.
> * **Rigour:** every inter-cell bound is proved, the domains are complete, the 584 polynomials are exact, and the LR∞ anchoring is valid for all larger Λ.
> * **Hostile scan:** 1534 random states per gap. min H/N_s = 0.19999969 (300), 0.19999998 (1000), 0.20000001 (5000). The true margin over c = 0.1999 is about 1e-4.
> * **Issues:**
>   * **I1:** the OO block also needs M17 (EO cancellation), κ > 0 at every Δ ([DV] per row), and the N_s normalisation, so OO is conditional; ++ and +− are unaffected.
>   * **I2:** Part L uses LR∞ / R̂ for 300 ≤ Δ_g < 1000, outside GAP_LARGESPIN's stated range. Valid, but say so.
>   * **I3:** the cached t0_checks.txt says c = 199/1000, a stale run.
> * **Unchanged scope:** M20 (the realised weight, only at 300 and 1000) and beyond-leading-Regge δ are not covered.

> **Update (LADDER_ALLJ.md, 25 Nov; independent audit): the §7 first-moment operator and the §8 all-J ladder pass every pipeline test, and both were re-derived independently.**
> * **First moment.** M^orig = (i/6)(W−18)H holds exactly for J = 0, 2, 4–8, both parities, symbolic Δ, all m. The +18 = 3(Δ+3) translation − 3(Δ−3) normalisation; the authors' −i(Δ−3)/2 sign is confirmed from first principles. In Wigner-normalised midpoint form the operator is (i/6)W, tridiagonal with off-diagonal −b_m, for all J. −K̂_z = W − 18 is re-derived for symbolic J, m, Δ [P given the standard K action].
> * **Ladder.** All 158 cases pass exactly (J ≤ 40, both parities, Δ = J+10/3 and J+1001+1/7): F₁ = F₂ = F₃ = 0, the closed forms, and the recurrence at every m. h = −2 is the mirror image F^(−2)_m = σF^(+2)_{−m}.
> * **Conservation step.** The Todorov contraction (spin −2) with K̂ gives exactly −½ × the recurrence for symbolic J, m, Δ. The endpoint terms decay as O(\|α\|⁻²).
> * **Grade.** Both are [P] given the two operator identities [z·K, E] = 6iE⁽¹⁾ and D_z^μ[K_μ, E] = 0 (authors' [D]), with consequences checked exactly against the pipeline.
> * **Note.** The commutator sum rule needs the translation term (structure of Ỹ_m), which §9 includes.

## 1. Setup and normalisations

| item | definition |
|---|---|
| theory | `𝒯`: a unitary CFT₃ (parity-even or parity-violating) with a unique stress tensor `T` |
| `c_T` | normalised so that a free real scalar has `c_T = 1` |
| Newton constant | `8πG := 256/c_T` (AdS units `L = 1`; from `c_T = 32L²/(πG)`) |
| **A1** (gap) | every single-trace operator other than `T` has twist `Δ − J ≥ Δ_g`, for all `J` (including `J = 0, 1, 2`) |
| **A2** (tree level) | `O(1/c_T)` truncation, `c_T ≫ Δ_g²`; loops dropped |
| gap scale | `Λ := √((Δ_g − 2)(Δ_g − 1))`: `298.50` (`Δ_g = 300`), `998.50` (`Δ_g = 1000`) |
| heavy kinematics | `m² = (Δ − J − 2)(Δ + J − 1)`; gap state `m_gap² = (Δ_g − 2)(Δ_g − 1 + 2J)` |
| cubic couplings | `⟨TTT⟩` in the basis `(c_T, t₄, α_T)`; bulk dictionary `√(t₄² + α_T²) = 240|ĝ₃|` (AdS units; `ĝ₃` in CHLPS normalisation) |
| quartic coupling | `g₄`: flat R⁴ coefficient in Stage 0's normalisation (type II Virasoro–Shapiro: `g₄/(8πG) = 2ζ(3)/M⁶`, `M² = 4/α′`). `g̃₄`: the CFT-defined light-side coefficient `−Ĉ₄^{ds}|_light = 2g̃₄` after Lemma A; `g̃₄ = g₄` given (F)+(P). |
| `c₁` | coefficient of a possible `ν`-independent `E × W³` light term `16πG c₁t₄`, collecting the `k = 2` (`Ĉ′_EE`) and `k = 3` (`Ĉ₃`) parts. **Not computed; a free real parameter.** |

RU §1 checked that the two results use the same conventions: gap, `Λ`, `G`, twist convention, and the `k = 2` object `Ĉ′_EE` in the
min-norm gauge.

**The fixed functionals** (fixed as of 22 Nov 2026; not to be changed).
* **`Ω^EE`** (Result D). It is built from CLS-D's light-ray functional (2311.04271):
  * `Ω^EE = ∫₀¹dp ψ₂ Λ⁴Ĉ′_{EE,Λp} + ∫₀¹dp ψ₃ ½Λ⁶Ĉ_{3,Λp}`;
  * the `k = 2` part uses the dual structure `EE`, de-shifted by `c(ν)⁻¹𝒜_g⁻¹`, in the **min-norm gauge**;
  * the `k = 3` part `Ĉ₃` (`o_H·q`, `j′ = 1`) has no de-shift;
  * NF02 weights `ψ₂, ψ₃` from `code/robust/functionals_noforward.json`, with `k₀ = 27.48509817` and `k(0) = 1`.

  Full definition: `RESULT_D_STATEMENT.md` §2.
* **`Ω^Q_up`** (Result Q, upper).
  * The δ = 0.2 functional `Q3-full-delta0.2` of `code/resultq/functionals_q3.json` (`A = 17.241`, `κ = 0.1099`).
  * It is realised by the de-shifted `Ĉ₄^{ds}`, `∂Ĉ₄^{ds}`, `Ĉ₅^{ds}` with the Stage 7 weights (`functionals_q7_xprec.json`), plus the
    smeared `k = 2, 3` parts with Q's weights.
* **`Ω^Q_low`** (Result Q, lower): the de-shifted `k = 4` value piece (`RESULT_Q_STATEMENT.md` §1, §2.4).

---

## 2. The theorem

> **Theorem (conditional; Einstein gravity through quartic order at `Δ_g ∈ {300, 1000}`).**
> Let `𝒯` satisfy A1 and A2 (§1). Assume the analytic inputs of §3:
> * the KKSDZ Regge bound with its loophole, and K1′;
> * `H_K`;
> * Adm, H7, L0;
> * Lemma A with D_LR;
> * (F)+(P);
> * `H_EFT`;
> * existence of the de-shift;
> * H17.
>
> Assume also **per-state positivity** of `Ω^EE`, `Ω^Q_up` and `Ω^Q_low` on every twist-allowed heavy block. Positivity is computed
> [NE + M] on the finite sets of §3.2. Elsewhere it is assumed through `H_M`, `H_R`, `H_edge`, `H_OO`, `H_RS`, `H_LJ(Δ_g)`, `H_LJ^Q`,
> `H_fwd^Q`, `H_smear^Q` and HYB. Let `c₁` be the (uncomputed, free) coefficient of the `E × W³` light term. Then:
>
> 1. **(cubic)** `√(t₄² + α_T²) ≤ 240·√(16L_A/(k(0)(1+δ)))/Λ⁴ · (1 + O(c₁Λ⁻²))`, with `L_A = ∫₀¹ψ₂(t)dt/(t² + 9/(4Λ²))`. In AdS
>    units:
>    * **`≤ 3.13·10⁻⁶` at `Δ_g = 300`** and **`≤ 2.89·10⁻⁸` at `Δ_g = 1000`**;
>    * `∓0.22 %` / `∓0.017 %` per unit `c₁`, and `±0.24 %` / `±0.021 %` from H17.
>
>    Equivalently, `|ĝ₃|²Λ⁸ ≤ 1.07·10⁴ / 1.43·10⁴`.
> 2. **(quartic)** Add one of two hypotheses:
>    * **`κ ≥ 0`** [D-sketch + M] for the upper end and **`c_L > 0`** [D-sketch + M] for the lower end;
>    * or simply **`t₄ = α_T = 0`**.
>
>    Then
>
>    `−O(Λ⁻²) ≤ g̃₄Λ⁶/(8πG) ≤ 17.241·log(2Λ/3) − 19.220 (Δ_g = 300) / − 19.265 (Δ_g = 1000) + ε_{c₁}`,
>
>    i.e. the upper end is `72.04 / 92.81 + ε_{c₁}`, with `|ε_{c₁}| = 2.03·|c₁t₄|Λ²`.
>    * Using item 1 to bound `t₄`: `|ε_{c₁}| ≤ 0.57|c₁|` (300) and `≤ 0.059|c₁|` (1000).
>    * For `ĝ₃ ≠ 0` and `κ ≥ 0`, the term `−κ|ĝ₃|²Λ⁸` (`κ = 0.110`) only tightens the upper end.
>    * The CFT-defined `g̃₄` equals the bulk `g₄` given (F)+(P).
>
> **Grade.**
> * The implication "positivity ⇒ (1), (2)" is **[D]** given the analytic inputs; the constants are [DV].
> * **Positivity itself is not established.** It holds at the grades of §3: two items are [I] (`H_RS` for Q's binding state; the Regge
>   window of `H_LJ(300)` for `J ≥ 34`), and the analytic inputs `H_K` and (F)+(P) are [D-sketch].
> * The Einstein instance of Lemma A has not been checked.
>
> Nothing is claimed at any other `Δ_g`, at `Δ_g = 100`, at loop level, or for the ratio `|ĝ₃|²/g₄` (no cubic closure from (2)).

**Remarks.**
1. **Dependence structure.**
   * Item 1 uses only `Ω^EE`, and so does not need Lemma A, (F)+(P), `H_EFT` or anything forward.
   * Item 2 needs `Ω^Q_up` and `Ω^Q_low`. It needs item 1 only through `ε_{c₁}`.
   * If item 1 is used to bound `ε_{c₁}`, item 2 at `Δ_g = 300` inherits `H_LJ(300)`, `H_M`, `H_R`, `H_edge` and `H_RS`. Carrying `c₁t₄`
     explicitly avoids that.
2. **Two deviations from RU §5's wording.**
   * The lower end for `ĝ₃ ≠ 0` additionally needs the sign `c_L > 0` of the lower functional's `|ĝ₃|²` slope (measured 0.255–0.49,
     `RESULT_Q_REVISION` §3.3). RU's wording covered only `κ`.
   * `B_eff(1000)` is written `−19.265`. The value `−19.2647` (REVIEW_RESULTQ §1) has been rounded to "−19.27" in earlier documents; to two
     decimals it is −19.26. The effect is `< 0.01` on `92.81`.
3. **Where the margins are.**
   * D's binding states: `+0.2038` (300) and `+0.2009` (1000) in units of `N_s`, against a functional built with δ = 0.2.
   * Q's binding pair (`++`/OO J = 0 at the gap): FF/N `+0.1426 ± 0.0017`, with derivative-piece margin `s_fail/s_d = 1.78/0.51 ≈ 3.5`.

---

## 3. Consolidated hypothesis table

"D"/"Q" mark which result needs the item. Grades follow RU §4, updated by HRS_CHECK (HRS), RESULT_Q_OG_CHECK (OG),
RESULT_Q_FERMION_HEFT (FH) and RESULT_D_CLEANUP (CL).

### 3.1 Setting and analytic inputs

| # | hypothesis | D | Q | grade | evidence | what would upgrade it |
|---|---|---|---|---|---|---|
| A1 | twist gap `Δ − J ≥ Δ_g`, all `J`; unique `T` | ✓ | ✓ | [A] | setting | — |
| A2 | tree level, `c_T ≫ Δ_g²` | ✓ | ✓ | [A] | setting | — |
| K1 | KKSDZ Regge bound (1904.05905), `J₀ ≤ 1` for spinning correlators | ✓ | ✓ | [Q], with KKSDZ's loophole (first/second-sheet boundary, their §4.2.6) | `H7_I14_LIGHTSIDE` §2 | closing the loophole |
| K1′ | uniformity of K1 towards the light-cone edges | ✓ | ✓ | [D-sketch] | review V | write it out |
| K2, K3 | kernel power counting; per-family Rindler positivity ⇒ swap lemma; validity `J₀ < 7/2` of the de-shifted `k = 4, 5` functionals | ✓ | ✓ | [Q]/[D]; `J₀ < 7/2` [D-sketch] | `H7_I14` §7; V4FINAL §5 | — |
| **H_K** | heavy OPE sum ↔ (corner, `ν`) interchange for the heavy-normalised de-shifted kernels (`Ĉ′_EE`; `Ĉ₄^{ds}`, `∂Ĉ₄^{ds}`, `Ĉ₅^{ds}`) | ✓ | ✓ | **[D-sketch]** ([D given `H_ip`]: raw per-`ν` heavy-sum convergence, `J₀ < ½`, not implied by KKSDZ) | RU issue 6; V4FINAL §2.4 | Proposition K for the heavy-normalised kernel, or `H_ip` |
| Adm | `EE` is an admissible `j′ = 0` dual structure | ✓ | ✓ | [D] | `EE_FOUNDATIONS` §1 | — |
| H7 | `−Ĉ′_EE|_E = 16πG/(ν² + 9/4)`, `Ĉ₃|_E = 0` | ✓ | ✓ | [D] | `H7_I14` §3–§5, §18–§23 | an independent small-`ν` check of `h₋ + h₊` |
| L0 | `[TT]` double zeros at `k = 2, 3`; contacts drop out | ✓ | — | [D] | `H0` §15; `PARITY_ODD_ROUTEA` §12 | — |
| **Lemma A + D_LR** | tree terms with `J₀ < 4` (Einstein, `E × W³`, MFT) are annihilated pointwise by `Ψ₄`, `Ψ₅`; the collinear set carries nothing | — | ✓ | **[D]**. MFT instance [P]. Mechanism **[DV] at `O(1/C_T)`** on the exact free-scalar and free-fermion `⟨TTTT⟩` (residual `≤ 8·10⁻¹⁵`; controls fail at O(1)). Also **[DV] on spinning Witten-diagram correlators** with a `J₀ = 2` spin-2 exchange (EIN). **Einstein `⟨TTTT⟩` itself not evaluated [I].** | V4FINAL §2.2; OG; FH; EIN | one swapped-form `Ψ₄` on Einstein `⟨TTTT⟩` at one `ν` |
| **(F)+(P)** | spinning flat limit and polynomial light side (`−Ĉ₄^{ds}|_{R⁴} = 2g₄`, `g₅` cancellation) | — | ✓ | [D-sketch] | V4FINAL §5 | a `J₀ = 4` spinning position-space correlator |
| **H_EFT** | (a) term-wise (P) of heavy-exchange double traces = their actual value; (b) no conspiracy in the heavy sum / contact tower | — | ✓ | (a) [DV] in the scalar-Mellin toy (`J_h = 0–6`), [D-sketch] spinning; **(b) [A]** | FH §4.5 (refines RU issue 9 and OG §5.2) | (a) a spinning toy; (b) — |
| De-shift | continuum C4/C5 solutions exist for all `β` | (C2: [D]) | ✓ | [NE] to `t = 30`; C5 spinning `j′ = 1` coefficients [I] | Stage 7; V4FINAL §5 | — |
| H17 | `|ĝ₃|²` light coefficient `k(0)(1+δ)` | ✓ | (κ) | [D-sketch + M], carried over from `q²` | `H0` §18; REVIEW_EE m4 | spinning `R³` exchange to `O(r²)` |
| **κ ≥ 0** | sign of Q's upper `|ĝ₃|²` coefficient (`κ = 0.110`) | — | ✓ if `ĝ₃ ≠ 0` | [D-sketch + M] | REVISION §3.3; RU issue 7 | spinning W³W³ |
| c_L > 0 | sign of the lower functional's `|ĝ₃|²` slope (0.255–0.49) | — | ✓ if `ĝ₃ ≠ 0` | [D-sketch + M] | REVISION §3.3 | same |
| **c₁** | `E × W³` light term (`k = 2` and `k = 3`) | ✓ | ✓ | **free parameter, carried explicitly** | RU issue 14 | Regge residue of the `E × R³` `T` exchange |

### 3.2 Per-state positivity: what is computed, and the hypotheses that cover the rest

| # | hypothesis | D | Q | grade | evidence | what would upgrade it |
|---|---|---|---|---|---|---|
| — | **computed set, D** | ✓ | — | [NE + M]. 1000: even `J ≤ 24`, 28, 32, 36, 40, 48; odd `J ≤ 31`; min `+0.2009`. 300: even `J ≤ 24`, 28, 32; odd `J ≤ 15`; min `+0.2038`. One value per state: `RESULT_D_STATEMENT` §5. | CL; HC; HC2; `EE_DG1000_LARGEJ` | — |
| — | **computed set, Q** | — | ✓ | [NE + M] given `H_RS`. Measured: `++` J = 0, 2, 4, 6, 8; OO J = 0, 2; `+−` J = 4–9, 11. Every other channel needs `|s_d| ≥ 17–30×` the measured `0.51` to fail. | BF2; BF3; PO; FWD | — |
| **H_RS** | analytic-radial engine image = crystal-cycle integral for EE beyond leading order | (window 28/32; OO beyond order 0) | ✓ (every `s_d`) | **[I]** where it matters: `ν < 6` for Q's binding state (`ν < 4` carries 89 % of the `s_d` weight; `ν < 2` is untested). [NE] at `ν ≥ 10`; [NE, weak] at 6–10; scalar and EE LO [DV]. At `J = 28` the entry-resolved test fails at 3–16σ. | HRS (refines RU issue 3); RU issue 2 | an EE `k = 4` crystal value at `J = 0`, `ν < 4`, to ~10–20 per harmonic |
| H_M | the near-gap chain is accurate on the computed states (incl. `K₂` validity) | ✓ | — | [M], validated [NE] on rows with direct data. `K₂` fails at `J²/m ≥ 2.2`; all such states are now on direct sub-rows (incl. `(1048, 48)`: `+0.915 ± 0.08`). The 1.6–2.2 band is untested. | HC §5; CL §0.3 | direct sub-rows in the band; interval chains |
| H_edge | light-cone-edge window `x ∈ (0.8, Λ/m]` | ✓ | (tolerance `≥ 760 c₂`) | [NE + M] where computed (300: `J = 4–24`; 1000: `J = 4–14`, 16, 17, 31, 36, 40); [M] elsewhere | HC §2; CL §0.2 | the remaining rows |
| H_OO | parity-odd even-`J` block | ✓ | ✓ (J = 0, 2) | [D given `Γ_OO = Γ_++`]. The identity is [P given the helicity lemma] at order 0 and [NE] beyond (`≤ 0.05 / 0.12 c₂`; `|Δs_d| ≤ 5·10⁻⁶`). It reduces to EE `++` and `H_RS`. | HC2; PO; FWD; RU issue 11 | a derivation beyond order 0 |
| H_R | beyond-leading-Regge at `m ≳ 1.3Λ`, large `b̂` | ✓ | (HYB) | [M] | `EE_NEARGAP` §6; `ROUTEA_LARGE_B` | spinning `O(r)` Regge kernel |
| **H_LJ(1000)** | near-gap large spin: `H_AG ∧ H_win (refined) ∧ H_edge-amp` | ✓ | — | **[D-sketch + M]**. The sign-free `H_win` is violated at `(1100, 100)`; the refined form rests on HLJ §5.2 [D-sketch]. | HLJ_ANALYTIC; HWIN_THIMBLE; RU issue 4 | spinning `O(r)` kernel at `β ≤ 0.53` |
| **H_LJ(300)** | same at 300, plus the `(336, 28/32)` window `ν < 8.3` | ✓ | (only if `t₄` is imported) | **[D-sketch/M] bulk; [I] Regge window for `J ≥ 34`** (no graviton data); `(336, 28/32)` window [M] (no validated direct datum below `ν = 8.3` at `J = 32`) | RU issue 4; HC2 §3; RU issue 2 | graviton window data at 300; a validated below-fold method |
| H_LJ^Q | large-`J` near-gap corrections below Q's tolerance (`≥ 1.1·10⁴ c₂`) | — | ✓ | [M], margin ×10–10³ | RU issue 13 | `s_d` at larger `J` |
| H_fwd^Q | unmeasured forward channels; the 300 J-drift (J ≥ 7 quoted from 1000) | — | ✓ | [NE + M] given `H_RS`; drift [I], harmless at δ = 0.2 | FWD; RU issue 12 | — |
| H_smear^Q | smeared graviton corrections at the binding states (`≤ 3·10⁻⁵` of B₂ vs tolerance `≈ 24 %`) | — | ✓ | [NE + M] | REVISION §1 | — |
| HYB | scalar skeleton + leading Regge elsewhere; `β > 30` | (hybrid) | ✓ | [M] | RU §4 | — |
| σ_hU | double counting of the helicity-odd correction | — | δ = 0.05 and cubic closure only | [I]; **not used by the theorem** | RU issue 8 | — |

---

## 4. What is proved, derived, numerical, assumed

* **Proved / certified / derived and checked.**
  * The constants `A = 3022.34`, `B_AdS = −6552.66`, `16L_A = 10734.61 / 14339.59` and the cubic bounds [DV]. They were reproduced
    independently by REVIEW_EE, and recomputed for this note from `16L_A` and `Λ`.
  * `2L_AdS = 72.041 / 92.814` [DV] (RU §1).
  * NF02's flat positivity certificate [certified] (one float-checked corner).
  * The crystal identity `EE ≡ q²` on CLS-D's physical polarisation plane [P]. Vanishing of `Ω^EE`'s leading-Regge helicity
    off-diagonal [P].
  * On the light-transform manifold `W_u·f₄ ≡ 0` [P]. The MFT instance of Lemma A [P].
  * The Lemma A mechanism on the free scalar and the free fermion at `O(1/C_T)` [DV].
  * The crystal saddle sits on a helicity eigen-polarisation [P], so `Γ_OO = Γ_++` at order 0 [P given the helicity lemma (D, DV)].
  * The ABJM arithmetic and the `c_T ↔ G` normalisation [DV].
* **Derived [D].**
  * The two implications (theorem items 1 and 2) given their inputs.
  * Adm, H7, L0; Lemma A + D_LR (Einstein instance unchecked); V(a).
* **Derived in sketch [D-sketch].**
  * `H_K`; K1′; validity `J₀ < 7/2`; (F)+(P); spinning `H_EFT`(a); H17, `κ ≥ 0` and `c_L > 0` ([D-sketch + M]).
  * The lower end `g̃₄ ≥ −O(Λ⁻²)`.
  * The refined `H_win` and the Regge-window structure behind `H_LJ`.
* **Numerical [NE], [NE + M].**
  * Per-state positivity on the computed sets (§3.2).
  * `H_edge` where computed; the `Γ_OO = Γ_++` identity beyond order 0.
  * `H_RS` at `ν ≥ 10`; `H_smear^Q`; `H_fwd^Q`; de-shift existence to `t = 30`; `H_EFT`(a) in the toy.
* **Model [M].** `H_M`, `H_R`, HYB, `H_LJ^Q`, `H_edge` on uncomputed rows, the `(336, 28/32)` window, and the large-spin and
  large-`b̂` regions.
* **Quoted [Q].** The KKSDZ Regge bound, with its loophole. CLS-D's definition of the light-ray functional.
* **Assumed [A].** A1, A2, `H_EFT`(b).
* **Open [I].**
  * `H_RS` below `ν ≈ 6` for Q's binding state.
  * The `H_LJ(300)` Regge window for `J ≥ 34`.
  * The Einstein instance of Lemma A (derived, never checked).
  * The C5 spinning de-shift coefficients.
  * `c₁` (carried as a free parameter rather than assumed).

## 5. Weakest links, ranked

1. **`H_RS` at `ν < 6` for Q's binding state (`++`/OO J = 0 at the gap)** [I]. Result Q.
   * `ν < 4` carries 89 % of the `s_d` weight; `ν < 2` is untested; no available crystal-cycle method reaches the needed precision.
   * A ×3.5 engine error in the derivative piece would break the δ = 0.2 functional at the gap.
   * HRS_CHECK confirmed the engine at `ν ≥ 10` and found nothing against it, but cannot test where it matters.
2. **The `H_LJ(300)` Regge window for `J ≥ 34`, and the `(336, 32)` window `ν < 8.3`** [I] / [M]. Result D at 300, and Q at 300 if
   `t₄` is imported.
   * There are no graviton window data at `J ≥ 34`.
   * The sign-free tolerance (`180–390 c₂`) is of the size of the only measured window values (`≈ 450` at `J = 24, 28`). It holds with
     margin only if the window correction stays leading-Regge shaped [D-sketch].
   * `J = 32` has no validated datum below `ν = 8.3`.
3. **The Einstein instance of Lemma A** [D, unchecked]. Result Q.
   * Its leverage on `g₄` is `10¹⁰–10²²`.
   * The mechanism is [DV] on two free theories, and the controls show what the zero depends on: primarity of `T`, the light-transform
     weight, the T-shift, microcausality. Einstein satisfies all of these [D].
   * But Einstein's (13) channel (`T` + `[TT]` with anomalous dimensions) is probed only by the free scalar, because the fermion is
     blind to it in this gauge.
4. **`H_EFT`(b), no conspiracy in the heavy sum / contact tower** [A]. Result Q. Untested; the term-wise series is at best asymptotic.
5. **`H_K`** [D-sketch]. Both results. The main analytic hypothesis; it has no criterion short of `H_ip`.
6. **(F)+(P)** [D-sketch]. Result Q. It also carries the identification `g̃₄ = g₄`. There is no position-space R⁴ normalisation control.
7. **`H_LJ(1000)`** [D-sketch + M]. Result D at 1000. The sign-free `H_win` is violated at `(1100, 100)`; the refined form's class
   tolerance is computed at 1000 only; `ν < 10` is unmeasurable.
8. **KKSDZ's loophole and K1′** [Q] + [D-sketch]. Both.
9. **`κ ≥ 0` and `c_L > 0`** [D-sketch + M]. Q, only for `ĝ₃ ≠ 0`. Not needed if `t₄ = α_T = 0`.
10. **Model-level items** [M]: `H_M` (the `K₂` band `J²/m = 1.6–2.2`), `H_R` (the `N₀` pathology at 300), `H_edge` on uncomputed rows,
    HYB, `H_LJ^Q`, and the C5 spinning de-shift coefficients [I]. RU does not list the last among the weakest, but its grade is [I].

## 6. What is NOT claimed

* **No other `Δ_g`.** Only 300 and 1000 were computed; "`Δ_g ≥ 300`" is not claimed.
* **Not `Δ_g = 100`.** Positivity is contradicted at model level at near-gap `J ≳ 18`.
* **No loop level** (A2).
* **No cubic closure** `|ĝ₃|²Λ⁸ ≤ …` from item 2, and no ratio `|ĝ₃|²/g₄`.
  * The earlier "`120–198 / 154–257`" is withdrawn. Its lower slope 0.49 is 1.9× CHLPS's flat optimum 0.26, inconsistent with the flat
    limit (RU issue 8).
  * At most a remark: `|ĝ₃|²Λ⁸ ≤ 2L_AdS/(κ + ¼) ≈ 200 / 257` [D-sketch + M].
* **The Einstein instance of Lemma A is unchecked.** Only free-field proxies at `O(1/C_T)` were evaluated.
* **Per-state positivity is not established in general**, and neither result is a theorem without it.
* **Not the δ = 0.05 functional.** Its margin `+0.0060` depends on `σ_hU`.
* **Nothing about theories that violate A1**, e.g. ABJM (§7).
* **Novelty only at search level** [I].

## 7. ABJM consistency check

* **Dictionary** (Stage 0; re-checked by RU §1):
  * `M² = 4/α′` and `L²/α′ = π√(2λ)`, so `Λ ↔ ML = 2^{5/4}π^{1/2}λ^{1/4} ≈ 4.22 λ^{1/4}`. The identification holds only at leading order
    in `1/(ML)`.
  * The type IIA `R⁴` term `ζ(3)α′³/32` gives `g₄M⁶/(8πG) = 2ζ(3) = 2.404`, independent of `λ`.
* **Placement.** `2ζ(3)` sits at **3.3 %** (`Δ_g = 300`) and **2.6 %** (`Δ_g = 1000`) of the upper end of item 2 (`2.404/72.04`,
  `2.404/92.81`).
* **`ĝ₃ = 0`.** This follows from **10d type IIA supersymmetry** (no `R³`), not from "maximal SUSY": ABJM at `k > 2` has `N = 6`.
* **No `c_T` confusion.** `c_T` does not enter the placement, because `2ζ(3)` is `λ`-independent and `Λ` depends only on `λ`.
  * The mass-derivative identity `c_T = (64/π²)∂²_m F` relies on the mass coupling to a current in the stress-tensor multiplet. It is used
    only in `ABJM_THEOREM`'s separate `c_T` result, with the same `c_T` (free real scalar `= 1`, `8πG = 256/c_T`).
  * In an `N = 2` theory `∂²_m F` gives a **flavour-current** coefficient, not `c_T`. The two must not be confused.
* **Status.** ABJM violates A1 (light Kaluza–Klein modes). The check uses only the flat image (F), and the bound is log-dominated.
  **A consistency check, not evidence for the theorem.**

## 8. Literature positioning

* **CLS-D (Chang–Landau–Simmons-Duffin, arXiv:2311.04271, v3 Feb 2025; JHEP 04 (2025) 016).**
  * They construct the spinning light-ray dispersive sum rules on `⟨TTTT⟩` and the dictionary to flat-space graviton sum rules in the
    bulk-point limit.
  * They list "the remaining steps for uplifting flat space bounds to AdS" as **future work**.
  * Results D and Q are precisely that programme. They are presented as **a first explicit-constant, conditional execution of CLS-D's
    uplift** for the cubic and quartic graviton couplings.
  * Priority on the method is CLS-D's. The new ingredients are:
    * the `EE` structure, which removes the curvature-induced helicity off-diagonal;
    * the de-shifted `k = 4, 5` functionals;
    * the pointwise Lemma A;
    * the near-gap computations.
* **CHMRS (Caron-Huot–Mazáč–Rastelli–Simmons-Duffin, arXiv:2106.10274)** is the **methodological precedent**: two-sided
  `log(Δ_gap)` bounds on contact couplings in AdS₄ with dynamical gravity, from uplifted flat sum rules, for scalar correlators. Q's
  light-side formula and `H_EFT` are of CHMRS's (5.3) type, and D's de-shift uses CHMRS's scalar kernel `g`. Credit is due to both
  CLS-D and CHMRS.
* **CHLPS (arXiv:2201.06602)** gives the flat D = 4 graviton bounds.
  * D's CHLPS-form rewriting, `|ĝ₃|²Λ⁸ ≤ 3022.34·log(2Λ/3) − 5327.21`, is about 10× weaker than the transplanted CHLPS optimum.
  * That optimum would need forward and helicity-flip sum rules in AdS, which are absent; it is asserted by CHLPS, not derived in AdS.
* **Parametric `Δ_g⁻⁴` bounds** on `⟨TTT⟩` deviations: CEMZ (1407.5597), Meltzer–Perlmutter (1712.04861), Costa–Hansen–Penedones
  (1707.07689), AKT (1805.07393; `t₄ ≲ ln Δ_g/Δ_g⁴`, unknown constant). Item 1 has a different parametric form, `√(ln Δ_g)/Δ_g⁴`, and
  explicit constants; the constants are not comparable.
* **To be discussed:** Kundu (2104.11238); the `N = 6` bootstrap (Binder–Chester–Jerdee–Pufu, 2011.05728).
* **Novelty [I] at search level.** No published bound on `R³`/`R⁴` couplings from `⟨TTTT⟩` sum rules was found (abstract-level search
  of 19 Nov plus one extended search of 22 Nov).

## 9. Open problems

1. **`H_RS` below the fold.** An EE `k = 4` crystal value at `J = 0`, `ν < 4`, to ~10–20 per harmonic. Routes: lower `Δ` with a
   higher-order engine, variance control for the Lemma-B factor `R`, or an analytic proof that the Γ-continuation of the `R`-weighted
   monomials equals the crystal-cycle integral.
2. **A validated below-fold method at large `J`.** This would give graviton window data at `Δ_g = 300`, `J ≥ 34`, and `(336, 32)` at
   `ν < 8.3`. The thimble sums currently fail the scalar `Im` check.
3. **The spinning `O(r)` Regge kernel** at near-gap `β ≤ 0.53`. It would make `H_win` [D] at both `Δ_g` and remove most of `H_R`.
4. **One swapped-form `Ψ₄` on Einstein `⟨TTTT⟩`** (O(G) `[TT]` data in every helicity structure). Also the (13)-channel/corner part of
   D_LR on a correlator with a spin-1 or spin-2 single-trace exchange (the free fermion is blind to it).
5. **`H_K`** (Proposition K for the heavy-normalised kernel) and the KKSDZ loophole.
6. **(F)+(P) in position space**, including an R⁴ normalisation control from a `J₀ = 4` spinning AdS contact diagram. Spinning
   `H_EFT`(a); any handle on `H_EFT`(b).
7. **`c₁`:** the Regge residue of the `E × R³` `T` exchange (`k = 2` and `k = 3`). It may vanish at leading Regge for `EE`.
8. **Spinning W³W³** (firms up `κ ≥ 0`, `c_L > 0`); `σ_hU`; the 300 J-drift of `s_d`.
9. **H17 for `EE`** (the `Π̂₄` admixture of `Ĉ′_EE`).
10. **`H_edge` on the remaining rows; `K₂` in the `J²/m = 1.6–2.2` band; interval versions of the near-gap chains.**
11. **Other `Δ_g`**, and `Δ_g = 100` in particular.
12. **Release checklist of RU §5:**
    * statements rewritten: done (22 Nov);
    * `H_RS` check: done, inconclusive at `ν < 6`;
    * spinning O(G) Einstein `Ψ₄`: not done (free-field proxies only);
    * direct sub-row for `(1048, 48)`: done (`+0.915 ± 0.08`);
    * BF_SPINJ entry-by-entry residuals: computed by RU itself, regrade adopted. The BF_SPINJ authors have not re-reported.

---

### Pointers

* Result D: `RESULT_D_STATEMENT.md` (one value per state, history of superseded values).
* Result Q: `RESULT_Q_STATEMENT.md`.
* Review: `REVIEW_UNIFIED_2026-11.md`.
* Checks: `HRS_CHECK.md`, `RESULT_Q_OG_CHECK.md`, `RESULT_Q_FERMION_HEFT.md`, `RESULT_D_CLEANUP.md`.
