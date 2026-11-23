# LARGESPIN_W: Route W (weight-shifting) to the large-spin region of Ω^EE at Δ_g = 1000. A conditional asymptotic theorem with explicit J\*, and where the route stalls

*23 Nov 2026. Unreviewed, not committed. No existing file was edited. `SMALLNU_WINDOW.md`, `code/smallnu/`, `LARGESPIN_L.md` and
`code/lsl/` were not touched.*
* *New code:* `ads4/code/lsw/`. Caches are in `cache/`; its `.gitignore` ignores everything except itself.
  * `lsw_common.py`: helpers (read-only imports of `n2sub/n2_common.py` and `largespin/ls_core.py`).
  * `lsw_weyl.py`: the Weyl-reflection parity check (§2.2).
  * `lsw_costs.py`: exact Lemma-I costs of every correction structure (§4.1).
  * `lsw_theorem.py`, `lsw_classwise.py`, `lsw_jstar.py`: Theorem W's criterion, the tolerances and J\* (§4).
  * `run_all.sh` reproduces everything.
* *Read-only inputs:* `prl/data/d_abel_X_1000.txt` and `kernels.json` (via n2_common); NF02 weights and hybrid pieces (via ls_core).
* *Literature, fetched from arXiv source:*
  * CHMRS 2106.10274: `sec_heavy.tex` §3.1–3.3, `app_heavy.tex`, `sec_sumrules.tex` (eqs. `moment appendix`, `regge subleading`,
    `gammafactors`, `magicmoment`);
  * Costa–Penedones–Poland–Rychkov 1109.6321: eqs. `diffops`, `diffops2`, `main`, and the d = 3 discussion;
  * Karateev–Kravchuk–Simmons-Duffin 1706.07813: source fetched; only its general statement is used.
* *Compute:* ≈ 0.6 CPU-h, one process, `ulimit -v 3000000`, `OMP_NUM_THREADS=1`.

Grades are as in UNIFIED_THEOREM_v3 §0. **"c₂ units"**: a coefficient κ of a structure S means that the image correction is
`δI = (κ/m²)·S(ν; β)`, relative to the leading raw image with ring weight 2 (LARGESPIN §3). **"Floor"**:
`f_h = N_flat,h/5 + min(0, curv_h)` per helicity. That is the certified flat floor (CERT_CLOSE) plus the exact leading-Regge curvature.

---

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

## 0. Summary

| # | claim | grade | explicit constants / J\* | residual |
|---|---|---|---|---|
| W0 | **Exact spinning representation.** The Ω^EE action on a ⟨TTTT⟩ heavy block with couplings λ is the CHMRS symmetric-frame integral (`moment appendix`). In it, the scalar product of three-point functions T(x, y) is replaced by `Σ_{μμ′} C_μ C′_{μ′} D_μ D_{μ′}[seed]`. Here the D are CPPR's (Δ, J)-independent first-order operators on the external points, the seeds are scalar products with shifted external dimensions, and C, C′ are coupling vectors in the CPPR differential basis. The project's crystal integrand (`repro_engine/integrand.py`) is this object at J = 0. | [Q] (CPPR eq. `main`) + [D] | — | — |
| W1a | **Why CHMRS's Regge series is even in 1/m.** `T(−x, −y; Δ, J) = T(x, y; 3−Δ, −1−J)` (the Weyl reflection (Δ, J) → (d−Δ, 2−d−J)). With `Δ+J−1 = m e^{β/2}` and `Δ−J−2 = m e^{−β/2}` this is exactly m → −m at fixed β. So `F(ε) = T(εX, εY; m = 1/ε)` is even. | **[D] + [DV]** (identity to 10⁻⁴⁰; odd coefficient 10⁻³⁴) | — | — |
| W1b | **Unequal-dimension seeds** (which CPPR produces) are even **iff Δ₁₂ = Δ₃₄**. The odd O(1/m) coefficient of log T is `≈ 2.2·(Δ₁₂ − Δ₃₄)` at the test point. It is antisymmetric under exchanging the left and right structures, and zero on diagonal pairs. | [DV] | — | whether the odd seed terms cancel after the D-operators for conserved identical T: **[D-sketch]** |
| W1c | **Congruence lemma.** If `M(m) = 𝔫(m) U(m)† g(m) U(m)`, with 𝔫 > 0 and U invertible and ν-independent (normalisation and coupling-basis maps), then `M ⪰ 0 ⇔ g ⪰ 0`. So every ν-independent O(1/m) normalisation or basis effect is harmless, at any size and uniformly in the ratios. Examples are the measured +21/m ν-flat term and the a_h (D_h-class) excess of N2. | **[P]** | none needed | — |
| W2 | **Ring structure at every order.** Each order of the block-side expansion is a finite sum of `η^i ν^{2l} ∂^n P(η)`, with ν-polynomial coefficients from γ-ratios. That is ring class A–C, with no window-smooth (F/G) class. The kernel's own r-expansion adds index-shifted harmonics `P_{l±1}, P_{l±2}` at O(1/m²). | [D] for the Regge-moment kernel + CPPR insertions; [D-sketch] for the EE kernel | — | coefficient magnitudes: not derived (§3.5) |
| W3 | **R1 (de-shift form) at O(1/m).** It is not needed: by W1+W1c there is no genuine O(1/m) term. **R1 at O(1/m²) is load-bearing only at J ≳ 200:** in stripped form, the ν² ring classes B and C2 are tolerated only at 1.4–7.5 c₂ (J = 250–344), against 290–4500 in de-shift form. | [D-sketch] (O(1/m)); [DV] (costs) | §4.2 table | **R1 at O(1/m²) for the ν² classes: [I]**, with a precise sufficient condition (§3.4) |
| **TW** | **Theorem W (conditional).** Assume (P0) the leading-Regge hybrid floor; (P1) parity (W1); (P2) the beyond-hybrid O(1/m²) correction is ring class with \|κ_s\| ≤ Λ_amp·\|κ_s^scalar\|, index-shifted pieces ≤ K_sh; (P3) de-shift form; (P4) O(1/m⁴) remainder ≤ K₄ (ν⁴ ring units). Then `M_Ω(O) ⪰ 0` for **every coupling vector** (2×2 ++/+−, and 1×1 for odd J) on every sampled H_LJ(1000) state with J ≥ J\*. | **[P given P0–P4]** (Lemma I + 2×2 criterion) + [DV] (sampled costs) | **J\* = 50 for Λ_amp = 3, K_sh = 10, K₄ = 10**, margin ×1. **J\* = 80–100** at margin ×2. J\* = 30 for Λ_amp = 1, K₄ = 1. J\* = 26 if the ν² classes carry the beyond-hybrid β² suppression. | P1 [D-sketch]; P2 structure [D-sketch], magnitude [NE] (J ≤ 100) / [M] (J ≥ 150); P3 [I] (bites only at J ≥ 200–300); P4 [M]; beyond-all-orders (H_AG) [NE + M] unchanged |
| C1 | Checks of the leading term. (i) `c(ν) = 32 γ₃²γ₄/(γ₁γ₂γ₇)` exactly (the Pochhammer reading of the EE normalisation as γ-factors of shifted seeds). (ii) The cost engine reproduces N2's scalar-subleading costs (J = 26: 0.118 vs 2 × 0.059; J = 100: 6.3·10⁻³ vs 2 × 3.16·10⁻³). | [P] (i); [DV] (ii) | — | — |
| C2 | Checks against graviton window data. The fitted non-ring constant r is consistent with 0, as TW predicts at every finite order: (1100, 100) r = −75…−92 ± 105–120 and +46…57 ± 86–164 (DIRCV §8); (364, 64)@300 r = +1 ± 1 / +1 ± 5. The fitted B coefficients (3.7–5.0 c₂/ν², J = 100) correspond to **Λ_amp ≈ 2–3** relative to the scalar 11/6. They sit ×10 (stripped) to ×100 (de-shift) inside the class-B tolerance. | [NE] | Λ_amp ≈ 2–3 | J ≥ 150: no data |

**Bottom line.**
* The route gives a **theorem-shaped reduction**, not an unconditional theorem.
* Weight-shifting (CPPR) plus the Weyl reflection plus a congruence lemma remove the O(1/m) problem entirely. These are the ring-amplitude
  pieces at O(1/m) whose stripped form would be fatal at J ≥ 200 (a_tol ≤ 0.2–0.7 in (a/m) units, §4.2). The O(1/m) removal is:
  * [P] for every ν-independent normalisation or coupling-basis effect;
  * [D-sketch] for the full spinning integrand.
* What remains is a finite list of O(1/m²) coefficients and the O(1/m⁴) remainder. Their exact costs are tabulated (Lemma I) on all
  H_LJ sample states. With scalar-sized coefficients (×3, as the J = 100 graviton fits suggest) the criterion closes for **J ≥ 50**.
* **Where it stalls.**
  1. The magnitudes of the spinning O(1/m²) coefficients were not computed (§6, item 1). That needs the EE-kernel polarisation
     contraction at next order.
  2. R1 at O(1/m²) for the ν² classes was not proved (§3.4). It matters only at J ≳ 200, where stripped-form ν² terms of scalar size
     consume up to 50–100 % of the floor.

---

## 1. What weight-shifting provides in d = 3 [Q + DV]

1. **Structures.** For identical conserved T in d = 3, a spin-J primary has the following ⟨TTO⟩ structures:
   * even J ≥ 4: two parity-even (++, +−) and one parity-odd (OO);
   * odd J ≥ 5: one, +−;
   * J = 0, 2: ++ and OO.

   This is the helicity count `|h₁ − h₂| ≤ J` with Bose symmetry, as used by the project (UT3 §2.5; A1P §2.1). Nothing new here.
2. **CPPR (1109.6321, eq. `main`).** The conformal partial wave for couplings in the differential basis is
   `𝒟_left 𝒟_right W_𝒪`.
   * `𝒟_left = H₁₂^{n₁₂} D₁₂^{n₁₀} D₂₁^{n₂₀} D₁₁^{m₁} D₂₂^{m₂}`.
   * W_𝒪 is the **scalar** partial wave for the **same** exchanged 𝒪(Δ, J), with external dimensions shifted:
     `Δ₁ → Δ₁ + m₁ + n₂₀ + n₁₂` and `Δ₂ → Δ₂ + m₂ + n₁₀ + n₁₂`.
   * The D_ij are fixed, (Δ, J)-independent, O(d+1, 1)-invariant first-order operators on (P_i, Z_i). They are quoted in
     `cache`-independent form in the source (eqs. `diffops`, `diffops2`).
   * In d = 3 the differential basis is over-complete but spanning (CPPR §5). Parity-odd structures are generated with the ε-operators
     D̃.
3. **Consequence in CHMRS's symmetric frame [D].**
   * CHMRS's exact heavy-block formula (`moment appendix`) is built from the Lorentzian shadow representation
     `∫⟨φ₁φ₂𝒪(x₅,n)⟩⟨𝒪^S(x₅,n)φ₃φ₄⟩`.
   * The D's act on external points only, so they pass through the x₅, n integral. Hence the spinning integrand is
     `Σ C_μC′_{μ′} D_μD_{μ′}[T_seed]`.
   * Equivalently, it is `T_scalar(Δφ = 3) × Σ C_μC′_{μ′} E_{μμ′}`, where E is a ratio of CPPR monomials
     `V₁^{2−a−b}V₂^{2−a−c}V₅^{−b−c}H₁₂^aH₁₅^bH₂₅^c/X₁₂²` (left), times the same for the right.
   * The project's exact crystal integrand (`repro_engine/integrand.py`, J = 0) has exactly this form, with the EE kernel and the light-ray
     x₀ integral. So **W0 is a re-reading of existing project objects, not a new computation.**

## 2. Parity: what is odd in 1/m, and why it does not matter

### 2.1 Why parity is the crux
* At O(1/m), a ring-amplitude term in **stripped** form is tolerated only at a ≤ 0.18–0.69 (J = 200–344), in units of (a/m)P
  (§4.2, column aA(R2)).
* The measured O(1/m) term is a = +21 (ν-flat; PE/HLJ_CLOSE).
* So any genuine O(1/m) ring term that is not of de-shift (𝒜) form would end H_LJ at J ≥ 200. R1 was introduced (N2 §2.5) for exactly
  this reason.

### 2.2 W1a: the Weyl reflection [D + DV] (`lsw_weyl.py` → `cache/weyl.txt`)

**Proof sketch** (complete for d = 3, scalar externals of equal dimension). In CHMRS's T, under (x, y) → (−x, −y):
* the ratio `R = (e−x)²(e−y)²/((e+x)²(e+y)²)` inverts;
* `𝒱_e(−x, −y) = −𝒱_{−e}(x, y)`;
* the prefactors (`(x²y²)^{−Δ̃}`, `(1 + 2x·y + x²y²)^{…}`) are even;
* using `(e+x)²(e+y)² = (e−x)²(e−y)²/R`, the flipped T equals T with `(Δ+J)/2 → 1 − (Δ+J)/2` and `J → −1−J`.

This is (Δ, J) → (3−Δ, −1−J). It sends `τ̄ = Δ+J−1 → −τ̄` and `τ = Δ−J−2 → −τ`, i.e. **m → −m at fixed β**. Hence
`F(ε) := T(εX, εY; m = 1/ε, β)` is even, which is CHMRS's "only even powers" remark, now explained.

**Checks (40-digit mpmath):**
* the identity holds at 4 random points to ≤ 2·10⁻⁴⁰;
* the odd O(1/m) coefficient of log F is −1.7·10⁻³⁴.

### 2.3 W1b: seeds with unequal external dimensions [DV], and the spinning statement [D-sketch]
* Adding the seed factors `|x₁₅²/x₂₅²|^{Δ₁₂/2}|x₃₅²/x₄₅²|^{Δ₃₄/2}` gives the following odd O(1/m) coefficient of log F:

  | (Δ₁₂, Δ₃₄) | (1, 0) | (1, −1) | (2, 2) |
  |---|---|---|---|
  | odd O(1/m) coefficient of log F | +2.200 | +4.400 | 10⁻³⁴ |

* **The odd part is ∝ (Δ₁₂ − Δ₃₄).** It is X, Y-dependent (a genuine ring-displacement insertion), antisymmetric under left ↔ right, and
  zero whenever the left and right seeds carry the same shift.
* In the CPPR basis the diagonal pairs (μ, μ) have Δ₁₂ = Δ₃₄. So the diagonal of g_μμ′ is even, and the odd part of g is antisymmetric at
  the seed level.
* The D's commute with the flip. They are O(3, 2)-invariant, built from dot products, and of even degree in each Z, so x → −x, Z → −Z
  leaves them invariant.
* **Not shown:** that the full g(m) (seeds plus D's, contracted with the EE kernel tensor) has an antisymmetric odd part. If it does,
  then for couplings with a common phase (the parity-even sector, real λ_++ and λ_+−) the odd part drops out of the Hermitian form
  exactly. In the parity-mixing sector it can contribute only through `Im(λ_μ λ̄_{μ′})`. **Grade: [D-sketch].**
* The EE kernel is invariant under the flip and under the Klein swap (13)(24) [D, by inspection of the published integrand]:
  * EE is built from V, X and H, with no ε-tensor;
  * `V_{0,42}²` and `V_{2,04}V_{4,02}` are 2↔4 symmetric;
  * `A = Nn²c₁₃/(m₂₄m₄₀m₂₀c₁₀c₃₀)` is (13)(24)-invariant, since `c₁₀ ↔ −c₃₀` and Nn → −Nn.

### 2.4 W1c: the congruence lemma [P]
* **Statement.** For any invertible U and any 𝔫 > 0, `λ†Mλ = 𝔫·(Uλ)†g(Uλ)`. So M ⪰ 0 on ℂⁿ ⇔ g ⪰ 0.
* **Consequence for the O(1/m) terms.** If g = g₀ + g₂/m² + … (W1), every O(1/m) term of M, written in the helicity basis normalised
  at leading order, has the form `𝔫(DV + V†D)`, where D = diag(leading images) and V = U₀⁻¹U₁ is constant. Such a term is an
  infinitesimal congruence.
  * Its diagonal is a normalisation of the raw leading image, which is exactly 𝒜-form.
  * Its off-diagonal is a ν-independent mixture of the two helicities' leading images.
  * Positivity is untouched for every coupling vector.
* **This is what the data show** [NE]:
  * the O(1/m) term is ν-flat and J-independent (a = +21, J = 0, 28–100);
  * the D_h-class excess at (1100, 100) is an a_h normalisation (N2-e, C3H σ_h test: +4230 ± 239 vs +4168 predicted).

## 3. Ring structure and R1

### 3.1 Block side [D]
* After rescaling x, y → x/m, y/m, E_μμ′ and the seeds are rational with convergent Taylor expansions at fixed (X, Y). Each order is
  polynomial in (X, Y) up to powers of X², Y².
* x^μ → ∂_p^μ and |x|^{2j} → shifts of a in the Fourier–Laplace identity
  `∫|x|^{a−d}P(η_{xy})e^{x·p} = 2^{a−1}π^{(d−2)/2}γ_a(ν)|p|^{−a}P(η_{py})`. The ratios γ_{a+2j}/γ_a are polynomials in ν².
* So each order is `Σ η^i ν^{2l} P^{(n)}(η)` at the state's own η = cosh β. That is ring class (A: P; B: ν²P; C: ∂_βP), with H₂ density a
  finite-order distribution on t = β (Lemma I; N2 §3.3).
* Negative powers of X², if any survive, give γ_{a−2}/γ_a, i.e. P/(ν² + c²) (class E). Class E is also tabulated in LARGESPIN; its
  costs are of the index-shifted size.

### 3.2 Kernel side [D-sketch]
* EESTEP §1.3 gives `Ĝ_b = Σ rⁿ Σ_m B^{(b)}[n, m] P_{l+m}`. Since r = |x||y| ~ 1/m², the kernel's O(r) terms are O(1/m²). They carry
  index-shifted harmonics P_{l+m} with |m| ≤ 2 and coefficients rational in ν.
* By N3, these are not ring class (their density lies inside the disc). They cost per §4.2 (column "sh", stripped).

### 3.3 Where 𝒜-form comes from
* Write the EE leading-r kernel, after polarisation contraction, as `Σ_s w_s(ν)P_{l+s}`. The block-side transform of harmonic l+s is
  `F_s(ν − is)` (the Fourier–Laplace identity is analytic in ν).
* The raw image is then `Σ_s w_s(ν)F_s(ν − is)`. It is exactly of de-shift form `𝒜[F]` iff F_{+1} = F_{−1} ≡ F as functions.
* At leading order this holds: the project's exact `h∓ = g·c(ν ± i)` [P given H7 §18.1].

### 3.4 R1 at O(1/m²): open, with a precise sufficient condition [I]
* **Sufficient condition:** the s = ±1 components of the EE leading-r kernel carry the same polarisation tensor, up to ν-dependent
  scalar factors. Then the spinning block-side transform factorises, F₊ = F₋ to all orders, and R1 holds for every block-side
  correction.
* This is checkable with `eestep/ee_quad.py`, by decomposing L_ν per polarisation set (EESTEP §1.4). It was not done here.
* The kernel-side terms (§3.2) are not 𝒜-form, but they are tolerated at K_sh ≥ 52 (stripped) everywhere (§4.2).

### 3.5 Coefficient magnitudes: not derived
* For the scalar, CHMRS (3.20) gives `κ_A = 86.75/6 − 8.0 = 6.46` (the −8.0 is N2's measured normalisation residual), `κ_B = 11/6`,
  `κ_C1 = 1/6`, `κ_C2 = 1/6`, `κ_C3 = −11.75/6`.
* The ν² coefficients `(2d + 3k − 1)/6` are **independent of Δφ**. This follows from (3.20):
  `(7d/2 + 3a − 1 − 6Δφ)` with `a = 2Δφ + k − d/2`.
  * So every weight-shifted seed has the same top-ν coefficients. These are CHMRS's "flat-space" coefficients (eq. `QPolynomials`).
  * The beyond-hybrid ν² content is therefore the difference between the CHMRS and hybrid ν² terms. It vanishes as β → 0, so it is
    O(β²)·O(1) [D-sketch]. This is the "ν² × β²" scenario of §4.3.
* The ν⁰ coefficients do depend on the seed. For the scalar formula, the ν⁰ coefficient of ηP (before the −8.0 residual) is −7.5 / −9.5 / +14.5 / +78.5 / +196.5 for Δφ′ = 1…5. After the D-operator
  corrections the spinning values are **unknown**.

## 4. Theorem W and its explicit constants

### 4.0 Statement
Let O be a heavy primary in H_LJ(1000) (near gap, m ≤ 1.3Λ), with couplings λ ∈ ℂ² for even J (++, +−) or ℂ¹ for odd J (+−). Write
the per-state action as `Q(λ) = Re Σ λ̄_h λ_{h′} M_{hh′}`. Assume:
* **(P0) Leading-Regge floor.** `M^{LR} = diag(f_h + slack)`, f_h > 0.
  * Certified flat part.
  * Curvature: [P + DV] for b̂ ≤ 119 (Lemma C′); computed positive ([NE]) on the rest of H_LJ.
* **(P1) Parity.** `M = 𝔫U†gU` with 𝔫 > 0, U ν-independent and invertible, and g even in 1/m at fixed β (W1).
* **(P2) O(1/m²) structure.** In the leading-normalised helicity basis, `U₀†g₂U₀ = Σ_s κ_s^{hh′} S_s + Σ_sh κ_sh^{hh′} P_{l+sh}`.
  * s runs over {A: ηP, B: ην²P, C1: P′(η), C2: (η²−1)ν²P′(η), C3: (η²−1)P′(η)}.
  * sh ∈ {±1, ±2}.
  * Coefficient bounds: `|κ_s^{hh′}| ≤ Λ_amp|κ_s^{sc}|` and `|κ_sh| ≤ K_sh`.
* **(P3) De-shift form** of the P2 ring terms (R1).
* **(P4) Remainder.** The O(1/m⁴) part is bounded by K₄ in units of `{η²ν⁴P, ηsinhβ ν⁴∂_βP}/m⁴`, and nothing beyond all orders
  (H_AG).

**Conclusion.** Let `d := Λ_amp Σ_s|κ_s^{sc}|c_s + K_sh Σ c_sh + K₄ Σ c₄`, where the c are the exact Lemma-I action costs
(`lsw_costs.py`) per unit coefficient. Let `H := f₊₊f₊₋/(f₊₊ + f₊₋)` for even J, or f₊₋ for odd J. If d ≤ θH with θ ≤ 1, then
`Q(λ) ≥ (1 − θ)·min_h f_h·|λ|²/2 ≥ 0` **for all λ**.

**Proof.**
1. Each structure has an H₂ density supported at t = β (ring class) or given explicitly (P_{l+sh}).
2. Lemma I [P] makes each entry's correction equal to its cost, computed exactly (sign-free).
3. Each entry is therefore bounded by d. For a 2×2 Hermitian matrix with diagonal ≥ f_h − d and off-diagonal ≤ d, positivity holds iff
   `(f₊₊ − d)(f₊₋ − d) ≥ d²`, i.e. iff d ≤ H.
4. The 1×1 case is immediate.
5. The O(1/m) terms drop out by P1 and W1c. ∎

**Grades.**
* The implication is [P] given P0–P4.
* The costs are [DV]: floating-point quadrature on the published grid, for sampled states. Behaviour between the sampled J and masses is
  [NE] (monotone trends).
* For a parity-mixing 3×3 block (++, +−, OO), Gershgorin replaces H by min f/3. That lowers tolerances by ≤ 1.5.

### 4.1 Per-unit costs (`cache/costs_1000.txt`, `theorem_1000.txt`)
Quadrature over the full range ν ∈ [0, Λ] (Lemma I is exact on [0, Λ]), using the published X and ψ₂. Uniform bound |κ_s| ≤ K on each of
A–C3, sign-free. Minimum over sampled masses:

| J | 26 | 34 | 50 | 64 | 100 | 150 | 200 | 250 | 300 | 344 |
|---|---|---|---|---|---|---|---|---|---|---|
| K_tol(A–C3, de-shift) | 0.68 | 1.6 | 5.3 | 9.2 | 56 | 63 | 133 | 251 | 272 | 1490 |
| K_tol(A–C3, stripped) | 0.68 | 1.6 | 5.3 | 9.0 | 21 | 5.3 | 3.2 | 3.6 | 2.1 | **1.0** |
| K_tol(index-shifted, stripped) | 5.8e4 | 9.9e4 | 1.1e4 | 5.2e3 | 3.1e3 | 1.1e3 | 197 | 89 | 59 | **52** |
| K₄ tol (stripped) | 1.1 | 2.3 | 13 | 16 | 45 | 142 | 354 | 911 | 853 | 1160 |
| O(1/m) ring amplitude a_tol, de-shift / stripped | 104/181 | 104/372 | 104/48 | 103/23 | 100/16 | 89/3.3 | 74/**0.69** | 60/**0.32** | 48/**0.21** | 40/**0.18** |

The uniform-K column is dominated by C2 and B (the ν² classes) at J ≤ 50, and by their stripped-form costs at J ≥ 150.

### 4.2 Class-wise tolerances (`cache/classwise_1000.txt`; H-normalised, minimum over masses)

| J | A | B (R1/R2) | C1 | C2 (R1/R2) | C3 | index-shifted | ×scalar vector (R1/R2) | ×scalar, ν² × β² (R1/R2) |
|---|---|---|---|---|---|---|---|---|
| 26 | 1.1e5 | 19/19 | 92 | 0.7/0.7 | 3.5e4 | 2.2e5 | **3.3/3.3** | 415/417 |
| 50 | 1.1e5 | 96/71 | 340 | 5.5/5.6 | 3.7e4 | 3.9e4 | 23/22 | 1140/1270 |
| 100 | 1.1e5 | 410/48 | 1350 | 69/37 | 4.0e4 | 1.1e4 | 140/23 | 2240/256 |
| 150 | 9.8e4 | 990/61 | 3020 | 69/6.5 | 4.3e4 | 3.9e3 | 225/17 | 2150/87 |
| 200 | 8.3e4 | 1600/22 | 5340 | 150/4.1 | 4.5e4 | 774 | 419/7.1 | 2450/32 |
| 250 | 6.8e4 | 2500/7.6 | 8270 | 290/7.5 | 4.7e4 | 337 | 682/3.5 | 2680/15 |
| 300 | 5.5e4 | 3600/5.1 | 1.2e4 | 310/3.9 | 4.9e4 | 221 | 812/2.3 | 2450/8.1 |
| 344 | 4.6e4 | 4500/5.0 | 1.5e4 | 2900/1.4 | 5.0e4 | 195 | **1520/1.8** | 3130/5.3 |

**Reading.**
* The ν⁰ ring classes (A, C1, C3) and the index-shifted harmonics are harmless at any plausible size: ≥ 10⁴ in de-shift form, ≥ 30–200 in
  stripped form.
* **Everything is decided by the ν² classes (B, C2):**
  * at small J, by their full-range cost. This is where the Regge series is used up to ν ~ m, and the flat partial wave (already exact
    in the floor) resums it;
  * at J ≥ 200, by **R1**. In de-shift form the scalar vector is tolerated ×419–1520. In stripped form it is tolerated only ×1.8–7.1.

### 4.3 J\* (`cache/jstar_1000.txt`)
J\* is the smallest sampled J₀ such that every sampled state with J ≥ J₀ (even J = 26–344, odd J = 33–301; m ≤ 1.3Λ) satisfies
d ≤ θH.

| premise constants (Λ_amp, K_sh, K₄) | form | θ = 1 | θ = ½ |
|---|---|---|---|
| (1, 10, 1) | R1, full ν² | **30** | 38 |
| **(3, 10, 10)** | **R1, full ν²** | **50** | **100** |
| (10, 30, 10) | R1, full ν² | 64 | none (fails at 300, 344 through K_sh) |
| (3, 10, 10) | R1, ν² × β² | 50 | 80 |
| (1, 10, 1) | R1, ν² × β² | **26** | 33 |
| (3, 10, 10) | R2 (no R1), full ν² | **none**: fails at 300, 344 | none |
| (3, 10, 10) | R2, ν² × β² | 50 | none (fails at 300, 344) |

* Below J ≈ 50 the binding term is K₄ (the full-range ν⁴ costs). That is an artefact of using the Regge series at ν ~ m. It is not a
  physical failure: the exact rows at J ≤ 48 exist (UT3 §5.2).
* At J ≥ 300 the binding terms are K_sh, and the ν² classes in stripped form.
* **Hoped-for J\* ≤ 100: obtained at margin ×2, for scalar-sized coefficients ×3 in de-shift form.**

## 5. Checks

### 5.1 Leading term [P/DV]
1. **c(ν) as a γ-factor ratio.**
   * `(4ν² + (2j−1)²) = 16|(2j−1)/4 + iν/2|²`. Hence `c(ν) = 32|Γ(5/4+z)²Γ(7/4+z)/(Γ(1/4+z)Γ(3/4+z)Γ(13/4+z))|²` with z = iν/2,
     i.e. `c = 32 γ₃²γ₄/(γ₁γ₂γ₇)` in CHMRS's `γ_a = |Γ((a−½)/2 + iν/2)|²`.
   * Checked numerically at ν = 0.3, 2.0, 7.1 to 10⁻¹⁵.
   * γ₇ is the scalar Δφ = 3, k = 2 magic normalisation. So the EE leading normalisation is a product of Fourier–Laplace factors of
     shifted seeds (a = 1, 2, 3, 3, 4), as the CPPR picture predicts [P identity; D-sketch interpretation].
2. **The theorem's leading term is the project's own.**
   * The floors are `ls_core.hyb_fast`: certified NF02 + curvature, with 2Ω_LR(β) = base (e.g. (1100, 100): 1.81·10⁻⁴).
   * Lemma-I costs via the grid ψ₂ reproduce `OmLR` to ≤ 1.4·10⁻⁴ relative (β = 0.05, 0.18, 0.5).
   * Stripped costs reproduce N2's R2/R1 = −6.4 at β = 0.18.
3. **Scalar subleading reproduced.** The CHMRS-shaped correction costs 0.118 / 6.3·10⁻³ / 2.9·10⁻⁴ of H at J = 26 / 100 / 344
   (de-shift). With H ≈ f_min/2 these equal N2 §3.3's 5.9·10⁻² / 3.2·10⁻³ / 1.4·10⁻⁴ of the LS floor.

### 5.2 Subleading against graviton window data [NE]
* **TW's prediction for the rows.** After the D_h (normalisation) class, the window correction is ring class (D_h, ν²D_h, xD_h′ shapes),
  plus O(1/m⁴). The non-ring constant r is 0 at every finite order.
* **(1100, 100), DIRCV §8, 11 nodes, full ring class.** r = −75…−92 ± 105–120 (λ-corrected) and +46…57 ± 86–164 (λ-free): consistent
  with 0.
  * The fitted B coefficients 3.7–5.0 c₂/ν² give **Λ_amp ≈ 2–3** against the scalar 11/6.
  * The class-B tolerance at J = 100 is 410 (de-shift) / 48 (stripped). The data are inside by ×10–100.
* **(364, 64) at Δ_g = 300.** r = +1 ± 1 / +1 ± 5 and B = 0.11 (N2 §2.4): consistent. No new 300 computation was made.
* **J ≥ 150:** no data (DIRCV §4). So P2's magnitude there is [M].

### 5.3 Parity checks
* §2.2 to 10⁻³⁴.
* Data: the O(1/m) term is ν-flat (a = +21 at J = 0 and 28–100; PE §2.3 / HLJ_CLOSE §1), as W1c predicts for a congruence.

## 6. Where the route stalls, and what completes it

1. **The spinning O(1/m²) coefficients (P2 magnitude).** Needed: expand the EE-contracted symmetric-frame integrand to second order.
   This is the J > 0 generalisation of `repro_engine/integrand.py` at fixed β, with `[−n·𝒱]^J` exponentiated as in CHMRS. Integrate
   the insertions with the p-derivative trick; the x₀ light-ray integral reduces through EESTEP's b-decomposition.
   * Output: κ_s^{hh′}(β), for 5 structures × 3 entries.
   * Expected effort: several days of symbolic work. Validate at J = 0 against the engine's order-a² images and at (1100, 100) against
     the DIRCV B fits.
2. **R1 at O(1/m²) (P3).** Check whether the s = ±1 harmonics of the leading-r EE kernel carry a common polarisation tensor (§3.4).
   This is cheap with `eestep/ee_quad.py`. If it holds, P3 is [D] for all block-side terms, and the stripped-form columns of §4 become
   irrelevant.
3. **Parity of the full spinning integrand (P1).** Show that the odd seed terms ∝ (Δ₁₂ − Δ₃₄) stay antisymmetric after the D-operators.
   Alternatively, check numerically that the J = 0 engine's order-a¹ image is a congruence (it is ν-flat [NE]).
4. **P4 / beyond all orders.** Unchanged (H_AG [NE + M]). TW does not split the window and the bulk: Lemma I acts on the full range. So
   the residual is exactly "the expansion is asymptotic, with these ring classes, uniformly in ν ≤ Λ". At J < 50 this is genuinely
   false as a full-range statement (the flat resummation), which is why J\* ≥ 50 there.

## 7. Effect on the Δ_g = 1000 positivity inventory (suggested; not applied)

| region | before (LARGESPIN, N2, DIRCV) | after Route W |
|---|---|---|
| near gap even J ≥ 50, odd J ≥ 33 (H_LJ) | [D-sketch + M]; residual R_win + R1 (stripped-form ring terms would end J ≥ 200) | **[P given P0–P4] for every sampled state with J ≥ 50, uniformly in coupling ratios**, with (Λ_amp, K_sh, K₄) = (3, 10, 10). P1 [D-sketch] and W1c [P] remove the O(1/m) part of R1. P2 is [NE] at J ≤ 100 (Λ_amp ≈ 2–3) and [M] beyond. **R1 survives only at O(1/m²), for the ν² classes, at J ≳ 200 [I].** Overall grade unchanged: **[D-sketch + M]**, now with explicit constants and J\* |
| J ≥ 150 | no data, [M] | still no data. The criterion is explicit: Λ_amp ≤ 7 (stripped) or ≤ 400 (de-shift) at J = 200, and ≤ 1.8 / 1500 at J = 344 |
| H_R′ (m ≥ 1.3Λ) | [D-sketch + M] | not tabulated here. `lsw_costs.py` extends by adding states |

**Suggested letter wording.**

> "At large spin the per-state action is an average of the leading-Regge action over the state's impact-parameter density. Every
> correction generated by the weight-shifted Regge expansion is supported on the state's own impact parameter, and its cost is computed
> exactly. Positivity then holds for every three-point coupling ratio once the second-order coefficients are within a few times their
> scalar values, which is what the J = 100 data indicate. This is proved conditionally; the coefficients themselves are not derived."

Do not write "proved for J ≥ 50".

## 8. Caveats
1. Only Δ_g = 1000, only the k = 2 piece (the k = 3 piece has weight u ≤ 0.59, as in LARGESPIN). States are sampled (40 states, J = 26–344).
2. Costs use floating-point quadrature: the Mehler–Dirichlet P with 600 nodes, and 4th-order finite differences for ∂_βP. The costs of
   the ν⁴ structures at small β carry ≈ 10⁻³ relative error; nothing is certified.
3. "Stripped" (R2) uses the published grid X, which is correct to t ≈ 22–25 (UT3 §2.2). The ν⁴ structures probe it hardest.
4. The O(1/m⁴) representatives are two shapes, not CHMRS's full O(1/m⁴) basis. The ancillary `regge_moments_order8.m` was not processed.
5. The 2×2 criterion assumes the same bound d on all entries. With entry-wise bounds, H can be replaced by the exact 2×2 condition.

---

## 9. Follow-up (coordinator request): R1 at O(1/m²), and the P2 magnitudes

*Same day. Additional code: `code/lsw/lsw_r1rank.py` → `cache/r1rank.txt`. It imports `eestep/ee6_lrprofile.py`, `ee3_deshift.py`,
`ee_quad.py` and `ee_blocks.py` read-only. Compute ≈ 0.1 CPU-h.*

### 9.1 The sufficient condition of §3.4 is false: the leading-r EE kernel has exactly three polarisation channels [DV]

**Test.** Let `L_ν[w₁, w₃, z₂, z₄](B)` be the r⁰ Cauchy coefficient of `K^EE/(N r²)`, i.e. EESTEP §1.4's leading-Regge position kernel,
computed with eestep's own code. One polarisation tensor (`L = 𝒯(pols; B)·ℓ(ν, B)`) holds iff the matrix `M[pol set, ν]` has rank 1 at
fixed B. The test used 7 random null-polarisation sets × ν ∈ {0.5, 1.5, 3, 6, 10}.

| B | singular values / σ₀ |
|---|---|
| 0.3 | 1, 0.24, 3.0·10⁻⁴, **2·10⁻¹¹**, 8·10⁻¹³ |
| 1.0 | 1, 0.37, 3.3·10⁻³, **4·10⁻¹⁵**, 3·10⁻¹⁶ |
| 3.0 | 1, 0.26, 4.2·10⁻³, **6·10⁻¹⁶**, 3·10⁻¹⁷ |

* **The rank is exactly 3 at every B.** The third singular value is small but 10⁸–10¹³ above the floor.
  So `L_ν = Σ_{i=1}^{3} 𝒯_i(pols; B) ℓ_i(ν, B)`: three tensor channels with three independent ν-functions.
* **Consequence.**
  * The factorisation argument of §3.3 cannot give R1.
  * At leading order the heavy symbols see only the de-shift combination `𝒜[cP]`: the project's exact
    `h∓ = g·c(ν ± i)` holds for both helicities with the same R∓ (EE_FUNCTIONAL §2). So at least one channel combination, call it
    ℓ_𝒜, carries the de-shift form. The other two directions are invisible to the leading symbols.
  * At O(1/m²), the spinning block's subleading polarisation tensors can project onto the other two channels. Their contributions
    `κ_i ℓ_i(ν; β)/m²` are **not of de-shift form** in general.
  * By P1 they cannot appear at O(1/m). That is consistent with the ν-flat O(1/m) data.
* **Status of R1 at O(1/m²): [I]**, now localised. Two things are missing:
  1. the identity of ℓ₂, ℓ₃: their harmonic content P_{l+s} and the rational ν-coefficients;
  2. the size of their coefficients, i.e. the contraction of 𝒯₂, 𝒯₃ with the O(1/m²) part of the heavy symbols.

  If ℓ₂ and ℓ₃ are index-shifted harmonics with rational coefficients, they fall in the "index-shifted (stripped)" class of §4.2,
  tolerated at K_sh ≥ 52 (J = 344) to ≥ 10⁴ (J ≤ 100). If they contain ν²-weighted pieces, the stripped ν² columns of §4.2 apply.

### 9.2 P2 magnitudes: not computed. What exactly remains
The next-order EE-contracted expansion was **not attempted beyond its set-up**. The work in §9.1 shows that it requires three ingredients
that do not exist in the project for J > 0.

1. **Conserved ⟨TTO_{Δ,J}⟩ structures for general J in d = 3**, in the helicity basis, with their (Δ, J)-dependent coefficients.
   * `repro_engine/tto_conserved.py` solves J = 0 only.
   * The general-J solve is the same linear algebra over ℚ(Δ, J) on the 14 CPPR monomials `V₁^{2−a−b}V₂^{2−a−c}V₅^{J−b−c}H₁₂^aH₁₅^bH₂₅^c`,
     with conservation imposed by the Todorov operator.
   * Estimate: about 1 day, sympy.
2. **The heavy symbols to O(1/m²) at fixed β.**
   * Evaluate those structures in CHMRS's symmetric frame at x₁ = −y/m, x₂ = −(x/m)/(−x²/m²), x₅ = e, with n fixed and
     `V₅^{J}` exponentiated through `J + ½ = m sinh(β/2)` (as in §2.2).
   * Taylor-expand the ratio monomials to second order. They are rational, so this is exact and mechanical.
   * The shadow-side structures follow by the Weyl reflection m → −m (§2.2), which also controls the parity (P1).
3. **Contraction with the three kernel channels and the p-derivative integration.**
   * Contract the symbols with 𝒯₁, 𝒯₂, 𝒯₃ of §9.1. These need to be extracted as explicit tensors, e.g. by fitting `L_ν[pols]` on a
     basis of degree-(2, 2, 2, 2) trigonometric polynomials in the polarisation angles. That is a cheap linear solve with ee6's code.
   * Then apply the polynomial insertions as ∂_p derivatives in the Fourier–Laplace identity: CHMRS's mechanism, validated for the
     scalar in N2 §3.2.
   * Output: κ_s^{hh′}(β) for the five ring classes, κ₂, κ₃ for the non-𝒜 channels, and the leading-order check `h∓ = g c(ν ± i)`.

**Validation targets** (unchanged):
* the J = 0 engine's order-a² images;
* the (1100, 100) ring-class fits: Λ_amp ≈ 2–3, r ≈ 0;
* (364, 64) at Δ_g = 300: B = 0.11, r = +1 ± 5.

Estimated total: several days. This is beyond the present budget, so **Λ_amp, K_sh and K₄ for ⟨TTTT⟩ remain [NE] (Λ_amp ≈ 2–3 from
the J = 100 fits) and [M] (K_sh, K₄).**

### 9.3 Updated J\* for Theorem W (`cache/jstar_1000.txt`, plus the run in §9.1's session)
With (Λ_amp, K_sh, K₄) = (3, 10, 10) or (2, 10, 10) and margin θ = 1, over the sampled H_LJ states:

| premise on form | sampled J covered | J\* |
|---|---|---|
| P3 (de-shift form) assumed, as in §4.3 | all J ≥ 50 | **50** (×2 margin: 100) |
| **P3 dropped (worst case: all O(1/m²) ring terms stripped)** | **50 ≤ J ≤ 250** (including J = 51, 101, 201, 301); fails at J = 300 and 344 | no J\*; a **window** [50, 250] (301 also passes) |
| P3 dropped, ν² classes with the beyond-hybrid β² suppression (§3.5, [D-sketch]) | J = 50–250 at θ = 1 (Λ_amp = 3); J ≥ 26 at Λ_amp = 1 | 50 (window) |

**Net.**
* R1 at O(1/m²) is now the one structural premise that separates "J\* = 50" from "the window 50 ≤ J ≤ 250". It is load-bearing only for
  J ≳ 300, i.e. near the top of H_LJ (m → 1.3Λ).
* It is localised to two extra kernel channels that are invisible at leading order (§9.1).
* Theorem W's grade is unchanged: [P given P0–P4] + [DV]. P3 is **[I]** (the cheap route is closed by §9.1). P2's magnitudes are **[NE]
  (J ≤ 100) / [M] (J ≥ 150)**.
