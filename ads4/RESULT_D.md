# Result D — conditional statement (consolidated; revised 18 Nov 2026 after the batched foundation review)

**Status: a conditional statement, not a theorem.** Revised after `REFEREE_FOUNDATION.md` (hostile review of this file and all
inputs added on 18 Nov). The referee independently reproduced every load-bearing number (exact skeleton moments, even-kernel
weights, `N^{K*}`, `A`, `B_AdS`, all six bounds, the 240 dictionary factor, parity-odd interference freedom, bulk-point density ratios),
and found that **the previous headline ("D-safe positive down to the gap, n = 1") was overclaimed**: it holds in a scalar-skeleton
model for a functional whose graviton version is not yet defined. This revision states what is shown at which level.

> **Update (H0 milestones M2/M3, `H0_SPINNING_FUNCTIONAL.md` §9–§11; not yet refereed).** Specialising Chang–Landau–Simmons-Duffin
> (arXiv:2311.04271) to d = 3 shows that light-ray-commutator sum rules for `⟨TTTT⟩` map at the scattering crystal to the **stripped** flat
> sum rules `B⁽¹⁾_k` (`d̃^J_44`), not to the B′ rules (full `d^J_44`) used below [D from Q + DV: unsubtracted k = 2, 3 rules vanish on
> double traces, as stripped `B₂, B₃` do and B′ do not; `f^±_{k;2}` at the crystal equals the stripped kernel for k = 2…7; explicit 4D
> crystal check]. Consequences:
> * **The D-safe functional as defined below (`Ω^{NF,SC}`, k = 2…7, B′ uplift) is not positive at leading order** for gravitons: flat action
>   −0.72 `N_s` at `+−`, m = 1.1Λ, J = 8 (reachable at Δ_g = 100) [NE given the identification]. Its "n = 1" holds only in the B′ model and is
>   **withdrawn** for the functional as defined. The reading-B support (`Π̂` data at x = 0.1) concerns position-space moments, not sum rules.
> * **Repair, route A′:** D-safe := `∫ψ₂Ĉ₂ + ∫ψ₃Ĉ₃`, the two unsubtracted CLS-D functionals (j′ = 0, f = v′−u′; j′ = 1, f = −1) with the NF02
>   weights. Flat limit = NF02 in the stripped basis, so its certificate applies [Q + D]. No k ≥ 4 objects and no forward terms: H3, H9b,
>   H15 are not needed. Near-gap model (exact scalar skeleton, stripped envelopes): **n = 1** at Δ_g = 100/300/1000, min F/N_s
>   0.2256/0.2054/0.2021 (0.2231 at 100 with an η-inversion estimate) [NE + M].
> * η-inversion (only Ĉ₂): leaves the graviton-pole light term unchanged under I9; δk(0) = −1.5·10⁻³/−1.6·10⁻⁴/−1.5·10⁻⁵, so A, B_AdS and
>   the bounds are unchanged at quoted precision (A ×1.0015 at most) [D-sketch + NE].
> * Still open for route A′: exact (2,4) dual structure in d = 3 (M2b), exact Mehler–Fock inversion (M3b), endpoint convergence (M4′),
>   light side (M5′), H4 on x ∈ [0.3, 0.7] (M6). H12′/H13′ results below were computed for the B′-uplift functional; the H12′ shape result
>   is functional-independent, the H13′ scans must be redone for route A′.
>
> * **Update (M2b/M3b/M4′, `H0_SPINNING_FUNCTIONAL.md` §R, §12–§14):** self-contained definition of Ω^{A′} in §R. Exact η-inversion
>   (Abel form, no analyticity in ν) [D, DV]: graviton-pole identity holds to 5·10⁻¹⁴, so A, B_AdS unchanged given I9; exact δk(0) =
>   −1.507·10⁻³/−1.638·10⁻⁴/−1.463·10⁻⁵; n = 1, min F/N_s 0.2231/0.2052/0.2021; bounds unchanged (factor ≤ 1.00075) [NE + M].
>   Flip-free (2,4) structures: stripped density holds for every flip-free choice [DV at the D = 4 crystal]; graviton covariant form
>   [D-sketch] (M2c). Ĉ₂, Ĉ₃ converge for J₀ < 2 [Q]; weights polynomial, no endpoint problem [D-sketch]; they vanish on all [TT]
>   double twists, so the tree-level light side is exactly Ĉ₂[G_T] + Ĉ₃[G_T] (contacts drop out) [D-sketch] — H9b, H10, H15 not needed
>   at tree level; H6, H7, H16 reduce to evaluating that single block (M5′).
>
> * **Update (M5′/M2c, `H0_SPINNING_FUNCTIONAL.md` §15–§16):** parity-even [TT] double traces drop out (double zero of CLS-D's
>   2 sin²(π(τ−10)/2), kept by polynomial weights) [D given Q]; parity-odd odd-twist [TT] zeros [D-sketch] (ties to H12′). By the CHMRS
>   Regge-residue formula (spinning carry-over [D-sketch]): −Ĉ₂|_E = 16πG/(ν² + 9/4) **exactly**, Ĉ₃|_E = 0, conditional on H7 — I9 upgraded
>   from a leading-Regge approximation to an exact residue. **A = 3022.3, B_AdS = −6552.7 and the bounds survive**; only k(0) is approximate
>   (H17 for Ĉ₂, Ĉ₃ to redo, M5′b). I14 swappability still [I]. M2c open: CLS-D's printed photon kernel and printed crystal images
>   disagree by (1−x²) on the 𝒱𝒱 term; if the images are right, route A′'s (2,4) structure has O(x²) helicity-flip leakage at the bulk
>   point (the stripped-dictionary conclusion is unaffected).
>
> * **Update (M2c/M5′b, §17–§18):** M2c resolved [DV]: the (1−x²) discrepancy is a missing factor 1/(m²−ν²) in CLS-D's printed
>   H₂₀H₄₀ crystal image; their printed photon kernel is correct, no helicity-flip leakage. Explicit covariant structure: Ĉ₂ uses q² with
>   q = (2H₂₀H₄₀ + V₀,₄₂²H₂₄)/(X₂₀X₄₀) (= 4S₀ + S₁ + S₂/16 in CLS-D's basis), flip-free with stripped ratio (1−x²)⁻⁴ [DV at crystal, D
>   covariantly]. H17 for Ĉ₂, Ĉ₃: exact coefficients need the spinning R³ exchange in AdS₄ to subleading Regge order (not computed); model
>   estimate |δk| = 4.4 % / 0.48 % / 0.043 %, i.e. bounds ±2 % / ±0.2 % / ±0.02 % (same budget as before) [D-sketch + M].
>
> * **Update (M6 partial, §19):** H4 stays **[I]** for route A′, now with explicit targets. Scalar calibration [NE]: on x ∈ [0.3, 0.7]
>   the exact remainder is c₂/m² + c₄/m⁴ with no 1/m term (local exponent 1.97–2.00), |c₂| ≈ 20–160. Tolerance [NE + M]: n = 1 survives
>   a graviton-specific density correction (1 + c/m²) on x > 0.3 for |c| ≲ 630 (worst, k-antisymmetric pattern) at Δ_g = 100, ≈ 3.4·10³
>   at 300, ≈ 2.7·10⁴ at 1000; a helicity- and k-common part is a positive rescale and cannot flip a sign. Evidence at x ≲ 0.1 only:
>   ⟨TTTT⟩ total c₂ ≤ 310, spin-specific residual ≲ 30. Closing M6: scalar Laplace expansion to O(1/m²) calibrated on the exact skeleton,
>   then graviton slow factors (~1–1.5 weeks).
>
> * **Update (H13′ for route A′, `ROUTEA_LARGE_B.md`):** hybrid action (flat NF02 + leading AdS curvature) positive on every state,
>   m ∈ [1.5, 10]Λ (+ 20/50/100Λ), b̂ ∈ [20, 200] at Δ_g = 300/1000, [40, 200] at 100, b̂ ≤ 1000 at 5 masses; min F_hyb/N_s 0.931
>   [NE]. Exact remainder on 62 states × 2 skeletons (J ≤ 150; b̂ reliably to 70/150/200 at Δ_g = 100/300/1000): min F/N_s 0.949,
>   |S|/F_hyb ≤ 122/m² [NE]. Large b̂: F_hyb ≃ A₁M₁(β)/Λ² − A₂/(2b̂³) with A₁ ≥ 355 (AdS term dominates for b̂ ≳ 4Λ^{2/3}); M₁ > 0 for
>   all β **[certified]**; explicit Arb-certified remainder bound gives **F_hyb > 0 [certified] on b̂ ∈ [115.5, 197] / [162.6, 597] /
>   [274.1, 1997]** (Δ_g = 100/300/1000; flat part float-only in the NF02 corner +−, J ≥ 400, m < 1.633Λ at Δ_g = 1000). Still [NE]/[I]:
>   b̂ ∈ [20, b̂₀), b̂ > 2Λ, |S| < F off the sample, η-inversion not included.
>
> * **Update (M6 step 1, `H0_SPINNING_FUNCTIONAL.md` §20):** scalar Laplace expansion of CLS-D's Ψ₂ around the scattering crystal to
>   O(1/m²) validated against the exact CHMRS skeleton [DV backbone + NE]: c₂ reproduced to ≤ 0.21 absolute (median 0.05 %) on 443 points
>   of x ∈ [0.3, 0.7], c₁ = 0 reproduced, residuals scale as the truncation predicts (rms 11.6 / 4.2 / 1.27 at Δ_g = 100/300/1000); 1/m³
>   remainder ≈ 10 in c₂ units at Δ_g = 100 vs tolerance 630. Uses CLS-D ≡ J = −1 inversion functional up to a ν-only factor [I].
>   Spinning shadow-coefficient route to O(1/Δ²) checked on the scalar. Graviton slow factors (step 2) not yet done: H4 still [I].
>
> * **Update (M6 step 2, `H0_SPINNING_FUNCTIONAL.md` §21):** graviton crystal machinery for Ĉ₂ built and validated at leading order
>   [NE]: shadow matrix β (residual ∝ Δ⁻³, crystal- and polarisation-independent), and the 2×2 action reproduces the stripped flat helicity
>   densities exactly (ρ₊₊/ρ_scalar = (1−x²)², ρ₊₋/ρ₊₊ ∝ d^J₄₄/((1−x²)⁴P_J)). Ĉ₃ (2,4) structure o found (same crystal normalisation as q).
>   **The O(1/m²) graviton correction δc₂ is not yet extracted** (the 1/m fits are unstable at the 10³ level). **New risk:** an O(1/m)
>   helicity mixing whose second-order effect is pointwise 5–691 in c₂ units on x ∈ [0.3, 0.7] — comparable to the Δ_g = 100 tolerance 630;
>   its integrated effect per state is not computed. H4 stays **[I]**.
>
> * **Update (M6 round 4, §22):** fit instability fixed (cumulant form; ν-only factors removed exactly by a heavy reference state, which is
>   §R.3's N̂₂). Remaining O(1/m) term, incl. the whole §21 helicity mixing, is imaginary; **prescription Ĉ₂ := Re[Ĉ^{CLS}₂/N̂₂]** (real and
>   imaginary parts are separately valid sum rules for real OPE data; mechanism [I]) — **load-bearing, flagged for review**: treating the
>   imaginary mixing as real breaks n = 1 at Δ_g = 300 and 100. Graviton-specific asymptotic corrections on x ∈ [0.3, 0.7] ≤ 470 (J ≤ 8),
>   ≤ 660 (J = 10). Per-state min eigenvalue of the ψ-integrated 2×2 (Ĉ₂ only, even J ≤ 10) [NE + M]: **n = 1 at Δ_g = 1000 and 300**
>   (min λ/N_s unchanged 0.2021 / 0.2054; correction may be scaled by (−14.6, +335) / ≈(−2.7, +4.6) before failing). **Δ_g = 100 undecided:**
>   n = 1 for J ≤ 4 (min 0.2256/0.2200), 0.103 for J ≤ 8, but J = 10 at m = 1.1Λ gives −2.7 (n = 1.5) — from extrapolating a Δ ≥ 400 fit to
>   m ≈ 105, a risk signal, not a computed failure. Not done: Ĉ₃, parity-odd sector (binding state +−, J = 5), J > 10. H4 stays [I].
>
> * **Update (certification for route A′, `CERTIFIED_NEARGAP_A.md`):** **Theorem A [certified]** — at Δ_g = 100 the route-A′ near-gap
>   action of the scalar-skeleton model (without η-inversion) has exact `p`-integral > 0 on the full grid (m/Λ ∈ [1.01, 3], J ≤ 10, both
>   channels; J = 11–16 sample) except ++, J = 0, m = 1.01Λ in each skeleton (229/230 for Δφ = 3, 178/179 for Δφ = 2); min certified
>   F/N_s +0.2251 / +0.2240 (+−, J = 5, 1.05Λ). With the exact η-inversion shift [DV, not certified]: min +0.2226 / +0.2214. Δ_g = 300: one
>   certified spot check; 1000 not attempted. Certified large-b̂ window now b̂ ∈ [91.2, 2Λ] at all three Δ_g. Continuity in m, m < 1.01Λ
>   not certified.
>
> * **Update (M6 round 5, §23):** with Ĉ₂ and Ĉ₃ (o_H structure) both corrected, even J ≤ 10, both helicities, per-state 2×2: **n = 1 at
>   Δ_g = 300 and 1000** (min λ/N_s 0.2054 / 0.2021; correction may be scaled by ±5 / ±37) [NE + M]. **Δ_g = 100 undecided**: the §22
>   negative J = 10 value and blow-ups were artefacts of the log form (withdrawn); in linear form K₃ converges for the scalar, but the
>   graviton series is not converged at m ≈ 100 (K₂→K₃ shifts ≈ 250–720 c₂ units vs tolerance 630), and x < 0.3 is outside the expansion.
>   No computed number breaks n = 1; deciding it needs a non-perturbative evaluation (J = 6, 8, 10; m/Λ = 1.0–1.3). **Definitional issue:**
>   the crystal computes the CLS-D functional Ψ_ν, while the base model and light side use CHMRS's C_ν (related by
>   Ψ ∝ g(ν)C(ν+i) + g(−ν)C(ν−i)); their ψ-integrated gap at J = 10 is up to 1.8 / 0.27 / 0.033 N_s at Δ_g = 100/300/1000 — §R must fix
>   which one Ω^{A′} is. Re-prescription judged legitimate [D] under three conditions (real OPE data in that basis — for the parity-odd
>   sector "Re" may become "Im" in mixed entries; ν-only N̂ with fixed phase; same phase and functional on the light side).
>
> * **Update (Δ_g = 100, §24):** **definition fixed** [D given the CHMRS identity]: Ω^{A′} is the CHMRS-type functional (Ĉ₂ = inverse
>   magic-block combination of CLS-D functionals, i.e. η-inversion with the exact kernel g P_{l−1} + g̃ P_{l+1}; Ĉ₃ likewise) — literal
>   CLS-D has a vanishing heavy normalisation at ν₀ = 5.979 and is withdrawn. Base model, Theorem A and the light side already refer to
>   this object; the remaining [I] (graviton shift relation) moves to the light side (A, B_AdS), not to n. **Direct (QMC) evaluation of the
>   exact 7-variable crystal integrand**, validated on the scalar vs the exact skeleton (rms 1.9 c₂ units on x ∈ [0.3, 0.7]; 8.4 on
>   [0.2, 0.85]; tolerance 630), works for ν ≳ 20. With direct graviton values on x ≥ 0.2 (Ĉ₂ + Ĉ₃, even J ≤ 10): **n = 1 at Δ_g = 100,
>   min λ/N_s 0.181** (J = 4, m = 1.04Λ; QMC error ≈ 0.005–0.02) [NE + M on x ≥ 0.2]; x < 0.2 (ν ≲ 20, 25–35 % of the weight for J ≥ 6)
>   covered only by the K₂ chain → **[I] overall**. No computed number breaks n = 1 at Δ_g = 100.
>
> * **Update (parity-odd sector, `PARITY_ODD_ROUTEA.md`):** covariant conserved odd ⟨TT𝒪⟩ built and matched to CHL's X^O at 9 states [DV];
>   no symmetric odd structure at J = 1, 3 (odd-J +− states start at J = 5). Real OPE coefficient; even–odd interference cancels
>   identically in the S⁰ sum at every Laplace order [DV], so Re acts on EE ⊕ OO. Leading order: odd exchanges give the stripped +− density
>   with a positive per-state normalisation N(Δ,J) = tan²(π(Δ−J)/2) (the ratio of the even/odd double-commutator factors sin²/cos²; with
>   those restored N_full = 1; the earlier "constant 2.31756" and Regge 3.30/3.77 are this function, PARITY_ODD_ROUTEA.md §12;
>   irrelevant for positivity) [NE]. **Odd-sector n = 1 at Δ_g = 1000 and 300** with J = 5, 7, 9 corrected (binding +−, J = 5: min
>   F/N_s 0.2023 / 0.2072; correction may be scaled by ≳ −210 / −57) [NE + M]. Δ_g = 100: n = 1 with K₂ (0.2217) but **not controlled**
>   (same small-ν non-convergence as the even sector). H12′ for route A′: shape and sign [NE].
>
> * **Combined even + odd (`PARITY_ODD_ROUTEA.md` §9):** sectors are block-diagonal; **n = 1 at Δ_g = 300 and 1000**, binding state +−, J = 5
>   (min 0.2072 / 0.2023), whole computed correction may be scaled by (−5.1, +5.0) / (−37, +39) [NE + M]; even J ≥ 12, odd J ≥ 11 uncorrected.
>
> * **Update (Δ_g = 100, x < 0.2; §24.8):** graviton/scalar ratio route fails at ν ≈ 15 (outliers don't cancel). Tolerance [NE + M]: with
>   direct values on x ≥ 0.2, **n = 1 at Δ_g = 100 iff the graviton correction on x < 0.2 is t × the K₂ value with t ∈ [−0.20, +1.25]**
>   (common scale factor; t = 0 → min 0.226, t = 1 → 0.181). Evidence consistent with t ≈ 1 ± 0.25, not shown. The x < 0.2 correction is
>   load-bearing; remaining route: a Picard–Lefschetz curved contour at ν = 5–20, validated on the scalar vs exact Mellin nodes.
>
> * **Update (odd sector at Δ_g = 100, direct QMC; PARITY_ODD_ROUTEA.md §11):** On x ≥ 0.2 the direct pointwise evaluation (`p8_odd_pw.py`) matches the K₂ Laplace chain (to ≤10⁻⁴ at Δ = 1600; within K₂ truncation at Δ = 400) [NE]. In the §24.5 protocol the odd sector gives n = 1 at Δ_g = 100, binding state +−, J = 5, 1.05Λ, min F/N_s = 0.2262 ± 0.0008; combined even+odd n = 1, min 0.1814 (even J = 4) [NE + M]. **This n = 1 is not robust:** replacing K₂ by K₃ on x < 0.2, or adding rows m = 1.05–1.3Λ, gives n = 0 through J = 9 (−1.74 and −0.48). The driver is x = 0.05–0.15 (ν ≈ 5–16, near the de-shift singularity ν₀ ≈ 6), where Laplace orders disagree by 10³–10⁴ c₂. The binding J = 5 state is positive in every variant (0.197–0.251). The J = 7 crystal sign at m ≈ Λ is also unstable under K₂ → K₃. **n at Δ_g = 100 stays [I]** in both sectors; Δ_g = 300, 1000 unchanged. New [NE] finding: non-even-integer Δ rows carry an untraced overall phase e^{−iπ(Δ−100)} in the state/reference normalisation, removed in `p8_analysis` (`P8_PHASE=1`); all earlier rows were even integers, so earlier results are unaffected.
>
> * **Update (phase, double-trace zeros, N; PARITY_ODD_ROUTEA.md §12):** (i) The e^{−iπ(Δ−100)} phase of §11 is a branch artefact: at the crystal the shadow-coefficient right-hand side takes log of the negative separations m36, m46 with +iπ each, against l3's own "+|·|" convention, so the computed β = e^{iπ(Δ+J)}·β_real [P]; verified to 0.0 on all 99 cached shadow entries, and the crystal normalisation phase is ν-only to 4·10⁻¹² [NE]. The even sector has the identical phase [NE]. Stripping it restores real β, so the three Re-prescription legitimacy conditions hold [P given NE]. **Prescription for every future row, both sectors, state and reference: logK − iπ(Δ+J)**, with an assertion on Im logK; n1_analysis/l9_analysis need it before any non-even-integer Δ row, and graviton/scalar ratios need it always. All earlier results are unchanged (Δ ∈ 2ℤ rows; N is a same-state ratio). (ii) Parity-odd [TT] double traces: the CLS-D double-commutator factor on odd structures is 2cos²(πτ/2), with double zeros at odd τ, where odd double traces sit at leading order, so their tree-level contribution is O(1/c_T²) [P given Q + P-sketch for odd τ]. (iii) The odd/even normalisation is not a constant: **N(Δ,J) = tan²(π(Δ−J)/2)** at 14 states including 2 out-of-sample and the 2 Regge points [NE]. It equals sin²/cos² of the double-commutator factors, so N_full = 1 [NE + P given Q]. Positive per state, so no route-A′ number changes.
>
> * **Update (light side: H7 and I14; H7_I14_LIGHTSIDE.md):** H7 upgraded [D-sketch] → **[D]**. Ĉ₃|_E = 0 because Ĉ₃ (f = −1) is the bare light-ray commutator, whose T-block light side is ∝ [H(n₁), H(n₃)] and the Einstein part of H is the identity [P given Q2, Q8, F]. −Ĉ₂|_E = 16πG/(ν² + 9/4): the leading Regge term of Einstein exchange is the exact AdS shock (H₂ propagator of dimension 2) with transverse polarisations, and the ν-dependent factors cancel between heavy and light sides because k = j = 2 [D, given (HN), the heavy-normalisation premise already used throughout, plus D2, D3]. The same route reproduces the CHMRS scalar formula 8πG/(ν²+d²/4) analytically. A double pole and the old "rational mismatch factor" are excluded. I14 upgraded [I] → **[P given K1–K3]**: both functionals need J₀ < 2; KKSDZ prove J₀ ≤ 1 for spinning correlators (with their stated loophole), so j′ = 0 is within the proven bound; swapping via per-family Rindler positivity + Cauchy–Schwarz + Fubini. Remaining gap: decay of the total de-shifted Ĉ₂ kernel at large η is [D-sketch]. Checks in code/h7i14/v1_lightside_checks.py (all pass). B_AdS corrected to −6552.66 (was quoted −6552.5; 2.5·10⁻⁵ relative, nothing changes). Open: D3 explicit computation; check (HN) for the spinning Ψ; whether Ĉ₃ has a spin-2 Regge moment (H17 only).
>
> * **Update (H7/I14 gaps; H7_I14_LIGHTSIDE.md §12–§17):** D3 proved: for kernels homogeneous of degree a−3 on the forward cone, Fourier = Laplace × e^{iπa/2} per pair, with identical angular/polarisation integrals, so N̂₂(ν) cancels identically [P]; numerical checks to ≤10⁻¹². So **−Ĉ₂|_E = 16πG/(ν²+9/4) is [P given (HN) + cited literature + D2]**. **(HN) is the remaining physics gap**: it holds by structure up to admixtures in span{P_{l−1}, P_l, P_{l+1}} [D]; Re removes the conjugation-odd ones [P], but a helicity-even P_l piece would survive. It is [I] at ν ≲ 15 and untestable with the crystal code there. A model failure (fixed transverse tensor structure) would shift **B_AdS by ~10–20%, not A** [M]. Fix: define Ĉ₂ with the spinning de-shift coefficients h_{−,0,+}(ν) from the leading-Regge kernel of Ψ_{2,ν}[q²], which makes (HN) true by construction. I14: large-η decay of the total kernel is η^{−1/2}/log η (correcting the earlier /log²η), proved given K1′, with numerics to β = 400; Fubini holds with a spare η^{−1/2}. Ĉ₃ has no spin-2 Regge moment (≤7·10⁻⁴ share at Δ = 1600) [D + NE].
>
> * **Update ((HN) resolved by redefining Ĉ₂ → Ĉ₂′; H7_I14_LIGHTSIDE.md §18–§22):** The leading-Regge kernel of Ψ_{2,ν}[q²] was computed in closed form [D]; it reproduces CHMRS's scalar g to 10⁻¹⁴ and the independent crystal ν-only factor exp(22i/ν + ~90/ν²) [DV]. It has **no helicity-even P_l piece (h₀ ≡ 0)** and ++ = +− [P given kernel]. But h₋/h₊ differs from the scalar g/g̃ (O(1) at ν ≲ 8, ≈ −180/ν³ at ν ≥ 20), so **the old scalar-g Ĉ₂ violates (HN)**: its graviton pole would be 16πG(1+ε)/(ν²+9/4) with ε = +0.51 at ν = 1–3, singular at ν₀ = 5.98, −0.0017 at ν = 20. **New definition Ĉ₂′ := two-term de-shift with the computed spinning h₋, h₊**; normalisation-independent, (HN) by construction, pole exact. **A = 3022.34, B_AdS = −6552.66, NF02/k₀ and the bounds 2.21e-4 / 3.13e-6 / 2.89e-8 are unchanged.** Heavy side: Δ_g = 300 and 1000 change by ≤ 0.003 / 10⁻⁴ (n = 1 stands); Theorem A (scalar skeleton) is unaffected. **At Δ_g = 100 the change is concentrated at the ν ≈ 4.9 and 9.8 nodes**: using Ĉ₂′ at all nodes, the direct-row minimum drops from 0.181 to 0.055 (J = 8) and the Laplace rows give n = 1.09; using it only at ν ≥ 15–20 changes the minimum by ≤ 0.002. Since the local de-shift is invalid at those nodes for both functionals, n at Δ_g = 100 stays [I], and the small-ν computation must now use Ĉ₂′. Open: non-local weight χ for Ĉ₂′; zeros/poles of h∓ near the real axis.
>
> * **Update (Δ_g = 100 small-x control via Picard–Lefschetz contour; H0_SPINNING_FUNCTIONAL.md §25):** A hybrid contour (linear disk of radius 3.5 through the saddle of Φ at ν_f, then downward PL flow with Re Φ = −ρ²/2, so Im Φ is constant and no branch locus is crossed) evaluates the exact integrand down to ν ≈ 9 (`n3_hybrid.py`). Scalar validation against exact Mellin nodes: 152 points, ν = 9.2–19.15, residual rms 19, max 82 c₂, i.e. ≥ 8× inside the 630 tolerance [NE]. Below ν ≈ 9 accuracy is unusable. New structural finding: the crystal saddle disappears in a fold at ν* ≈ 11.3 and a σ-conjugate saddle pair takes over (Airy-type), which is the origin of the "ν-soft" failures. **With Ĉ₂′ (direct values at x = 0.1, 0.15; K₂ only at x = 0.05): n = 1 at Δ_g = 100, λ_min/N_s = 0.2206 (J = 4, 1.04Λ)**; per J = 4/6/8/10: 0.221 / 0.279 / 0.392 / 0.502 [NE + M]. The J = 8 value 0.055 from the previous update came from K₂ at all small-x nodes; it is superseded. The graviton correction factor t is 0.70–0.98, and the remaining uncontrolled node x = 0.05 (ν ≈ 4.9) can be scaled by t ∈ [−32.6, +3.62] times its K₂ value without losing n = 1. **Status of n = 1 at Δ_g = 100: [NE + M] on x ≥ 0.1, [I] at the single node x = 0.05** (wide margin), plus a definitional caveat: at ν < ν₀ the chain uses the local de-shift, and the exact non-local inverse is not implemented. Old scalar-g Ĉ₂ results (n = 1, 0.2067) kept for comparison only. The odd sector J = 5, 7, 9 has also been evaluated on the contour at x = 0.1, 0.15 (`cache/n3_odd_pts_*`), not yet run through p8_analysis. Still running: a Ĉ₂′ tolerance bisection over all of x < 0.2.
>
> * **Update (Ĉ₂′ consolidation; H7_I14_LIGHTSIDE.md §23):** (i) Closed form h∓ = N_com(ν)·R∓(ν), with R∓ rational and N_com a ratio of Γ-functions without zeros/poles for |Im ν| < 11/2 [P given §18.1]; matches the v3 table to 2·10⁻¹¹. **No zeros or poles of h₋, h₊, h₋+h₊ on the real axis (0, 300]**, apart from the pole at ν = 0 that the scalar g also has; the kernel is finite there. The old de-shift singularity is ν₀ = √143/2 exactly; for Ĉ₂′ that zero moves off the axis to 6.3458 − 0.7923i. (ii) Non-local weight χ via the generalised Abel form [D]: the η-inversion / graviton-pole identity holds to ≤ 8·10⁻¹² at Δ_g = 100/300/1000, so A and B_AdS are unchanged; flat limit (NF02, k₀) unchanged; δk(0) is ≈ 4.7× the old ηP value but Λ²δk(0) → 0; bound at Δ_g = 100 moves 2.21 → 2.22·10⁻⁴ (inside the H17 budget), 300/1000 unchanged [NE + M]. **New definitional issue [I]:** the χ equation has a 3-dim null space (from the winding of ρ = h₊/h₋, 7π vs 3π for the scalar g), living at ν ≲ 18; it is invisible to the leading heavy action and the graviton term, but whether it drops out of the full functional depends on analyticity of Ĉ′_ν in |Im ν| ≤ 1. The old scalar-g definition already had a 1-dim null space (not noted in H0 §24.1). (iii) Parity-odd with Ĉ₂′: n = 1 at Δ_g = 300 (0.2072) and 1000 (0.2023), and combined n = 1 at both [NE + M]. (iv) **Δ_g = 100 with Ĉ₂′, odd + even, main protocol** (direct x ≥ 0.2, contour at x = 0.1, 0.15, K₂ only at x = 0.05): odd n = 1, min 0.2392; combined n = 1, binding even J = 4 at 0.2206, x < 0.075 window [−32.6, +3.62] [NE + M]. **Not robust:** with extra rows m = 1.05–1.3Λ (K₂ on all x < 0.2 there) odd J = 9 at 1.2Λ goes to −0.55 and the window shrinks to [−0.93, +0.60]. The whole failure is the x = 0.05 node at ν ≈ 5.4–6.4. **n at Δ_g = 100 stays [I]; the one remaining computation is contour values at ν ≈ 5–6.5 with the non-local χ.**
>
> * **Update (Ĉ₂′ tolerance bisection, Δ_g = 100, even sector; `cache/n3_analysis_spin.txt`):** with direct values wherever computed, scaling the whole graviton correction on x < 0.2 by t keeps n = 1 for **t ∈ [−0.82, +1.80]**; the measured direct/K₂ factors at x = 0.1, 0.15 are 0.70–1.0 (Ĉ₂′) and 0.69–1.25 element by element, inside the window. On x < 0.075 alone the window is [−32.6, +3.62]. [NE + M]
>
> * **Update (uniqueness of the Ĉ₂′ weight; H7_I14_LIGHTSIDE.md §24):** the 3-dim null space of χ is not fixed by any analyticity, growth or scalar-limit condition [D + NE]; it is fixed by **stating the gauge as part of the definition** (min-norm: X − ψ₂ ⊥ null space; now in H0 §R). Lemma N [P given observed regularity]: null vectors annihilate the Â_h-image of every strip-analytic density, so A, B_AdS, the flat limit and the leading Regge heavy action are gauge-independent. Beyond leading Regge the gauge dependence is [I] but numerically small: n = 1 at Δ_g = 300/1000 gauge-independent to 4 digits; at Δ_g = 100 (main protocol) the combined minimum 0.2206 stays in [0.1915, 0.2256] over natural gauges (in [0.2196, 0.2212] without K₂ imaginary parts). The freedom cannot (admissibly) fix the odd extra-row failure at x = 0.05. Side finding: K₂ imaginary parts at x ≤ 0.1 are 4–7× the exact leading-Regge value, so non-local chains must not use them. The old scalar-g definition had one null direction, unique in the CHMRS class (now stated in H0 §R).
>
> * **Update (batched hostile review of the route-A′ chain; REVIEW_ROUTEA_2026-11.md):** No fatal error; A, B_AdS and the bounds reproduce exactly (2.209·10⁻⁴ / 3.132·10⁻⁶ / 2.891·10⁻⁸) and the step positivity ⇒ bound is valid [D]. **Grade correction: the Δ_g = 300 and 1000 bounds are conditional, not [NE + M].** Per-state positivity (n = 1) is [NE + M] only near the gap for even J ≤ 10 and odd J = 5, 7, 9; for gravitons elsewhere it is [M] + [I]. Serious findings: **S1** near-gap graviton corrections grow ≈ ×2.5 per two units of J at Δ_g = 300 (off-diagonal 0.001 → 0.088 from J = 4 to 10), so J ≥ 11 is load-bearing, not "safe by margin" (crude extrapolation reaches the margin at J ≈ 15–16; ~10× smaller at 1000); **S2** the large-b̂ / large-J region (incl. the "certified" window b̂ ∈ [91.2, 2Λ]) is scalar-skeleton model only, graviton corrections beyond leading Regge never computed (the reviewer's hybrid scan of m ∈ (Λ, 1.5Λ), J > 20 is positive, min 0.2056 / 0.2009); **S3** the model M (scalar skeleton + stripped flat densities + K₂ graviton correction interpolated in 1/m) is now stated; ~10–29% of the near-gap correction at 300 comes from unconverged Δ = 100/200 rows (dropping them keeps n = 1); **S4** Ĉ₃ is defined inconsistently (no de-shift in §R.5, scalar k = 2 de-shift in the chains); removing the shift keeps n = 1 but changes the J = 10 correction 0.072 → 0.021. Minor: the twist condition is Δ − J ≥ Δ_g for all J (§1's "Δ_g − 2 for J > 2" is wrong); a linear-in-t₄ light-side term c₁ is unaddressed (|c₁| ≲ 10 assumed, [I]); the Ĉ₂′ δk(0) = −7.66·10⁻⁴ / −6.83·10⁻⁵ replaces the ηP values; H7 is [D] overall because the §18.1 kernel is [D]; "Δ_g ≥ 300" means the two computed values 300 and 1000. **Supported statement:** under the hypotheses listed in the review §2, √(t₄²+α_T²) ≤ 3.13·10⁻⁶ (Δ_g = 300) and 2.89·10⁻⁸ (Δ_g = 1000), ±0.4% / ±0.04%, *conditional on per-state positivity for gravitons at near-gap J ≥ 11 and at large b̂*. Replacement tables for route A′ are in the review §3. Fixes in progress: S1 and S4.
>
> * **Update (S1/S4 fixes — n = 1 at Δ_g = 300 FAILS at large J; S1_S4_FIXES.md):** (S4) Ĉ₃'s spinning leading-Regge kernel is exactly l-isotypic (σ_h c₃(ν)P_l, ++ = −(+−), no P_{l±1}) [P given kernel reading (D) + DV ≤ 4·10⁻¹³], matched by the crystal integrand; **Ĉ₃ := Re[Ψ₃/N̂₃] with no de-shift** ((HN₃) exact). The scalar k = 3 shift in the chains was spurious and is removed; light side, flat limit and A, B_AdS unchanged; binding state 0.2072 → 0.2052 (300), 0.2023 → 0.2021 (1000). Chain defect found and fixed: l9_neargap interpolated off-diagonal Γ_off between Δ rows without aligning arbitrary basis signs (overlaps −0.94 to −1.00), partially cancelling them; aligned, the J = 10 off-diagonal at 300 is 0.227 (was 0.088). (S1) Near-gap graviton corrections for even J = 12–20 (2×2) and odd J = 11, 13, Ĉ₂′ + fixed Ĉ₃, K₂ rows Δ = 400/800/1600, validated by direct RQMC at Δ = 400 and 336: **at Δ_g = 300 the per-state minimum is 0.46 / 0.23 / −0.16 / −0.73 / −1.33 at J = 12/14/16/18/20 (K₂); direct rows give J = 18: −0.44 at its own mass 1.12Λ, J = 20: −0.49. So n = 1 at Δ_g = 300 is contradicted [NE + M].** Mechanism: the real O(1/m²) helicity off-diagonal at x ≈ 0.05–0.15 against a near-gap base action that is small at large J (the base is the exact flat NF02 action to ≥ 99%; B₂ = 0.018 at J = 20); ratio grows ≈ J^3.5/m². Crossing at J ≈ 15–17, m ≈ 1.06–1.25Λ (b̂ ≈ 30), inside the near-gap region. With positivity only from m ≈ 1.2–1.35Λ the bound would weaken to 6.3–8.6·10⁻⁶, and only if J > 20 is positive (single-row transfers [M] say J = 24/28 are negative up to m ≈ 1.65–2Λ). **At Δ_g = 1000: n = 1 holds for J ≤ 20 [NE + M] (min 0.2021) and J = 24, 28 [M], but J = 32 is −0.45 at 1.04Λ from the off-diagonal alone [M]** (crossing J ≈ 30–32, b̂ ≈ 60) → n = 1 at 1000 is [I] with a modelled failure. **Consequence: with the current weights (ψ₂, ψ₃ = NF02), Result D is not established at any Δ_g; the Δ_g = 300/1000 bounds above are withdrawn pending either a verification overturning this or a redesign of the weights.** The x = 0.05 question at Δ_g = 100 is moot until then.
>
> * **Update (x = 0.05 at Δ_g = 100, two-thimble contour; NEARGAP100_SMALLNU.md; stopped early, now secondary):** below the fold the crystal integral is the sum of the two σ-conjugate thimbles (intersection numbers 1, 1) [NE]; shared-seed state/reference runs give ratios to 1–5%. Scalar validation: clean at ν = 9.2 (Δc₂ = −16 ± 51); at ν = 5.45 not clean (+400–450 ± 250 c₂, phase residual ≈ −0.1 rad), so the method is not certified below ν ≈ 7 to better than ~500 c₂. Odd sector at x = 0.05 computed: main protocol n = 1 (min 0.2426); with extra rows n = 0, J = 9 at 1.2Λ now −1.345, driven by Ĉ₃ (direct = K₂ to 0.3%) at ν = 5.91, next to ν₀ = √143/2 — **but these runs still applied the scalar-g de-shift to Ĉ₃, which S1_S4_FIXES.md has since shown is spurious (Ĉ₃ needs no de-shift), so this failure is likely an artefact of the old Ĉ₃ definition [I].** Even x = 0.05 points not computed. Tools (n4_*.py) kept for reuse at ν ≈ 5–11.
>
> * **Update (large-J off-diagonal: genuine, structural; LARGEJ_OFFDIAG.md):** No bug [NE; r7_checks]. The real ++/+− off-diagonal is a **leading-Regge effect**: Γ_off/m = β²·g(νβ) with one universal function g (peak ≈ −1.3 at νβ ≈ 3.7), on which 26 K₂ rows (J = 6–32, Δ = 336–1600, β ≤ 0.11) collapse to ±0.1 [NE]; direct RQMC ladder (336,18) → (672,36) at fixed b̂ ≈ 33 grows 1.96× (leading-Regge prediction 1.94×) [NE]. Consistent with the flat limit: at fixed b̂ and m/Λ it falls as 1/Λ². The NF02 flat base in the wedge is the universal tail 1023.6/b̂³, so the ratio ≈ 3.7·b̂³/Λ² crosses 1 at **b̂_c ≈ 0.65·Λ^{2/3}** (≈ 30 at Δ_g = 300, ≈ 65 at 1000), and the 2×2 is negative on the whole band b̂ ≳ b̂_c at every m/Λ ∈ [1.1, 5] [M]. **Re-weighting within the NF02 LP family cannot fix it [M]**: covering b̂ ≤ 60 at 300 costs little (bound 3.32·10⁻⁶) but b̂ ≤ 90 at 300 and b̂ ≤ 200 at 1000 are infeasible, while the twist-allowed band reaches b̂ ≈ Λ–2Λ. Reason [D-sketch given universal g]: where AdS curvature dominates, off-diagonal/diagonal → 24|G₁| ≈ 300 for any weights (G₁ = ∫t g(t)dt ≈ −12.5). **Status: route A′ with Ĉ₂′ + Ĉ₃ is obstructed at Δ_g = 300 and 1000 [M]**, with the caveat that the infeasibility uses g at β ≈ 0.2–0.3, beyond the validated β ≤ 0.11 [I]. Open ways out: (i) the exact off-diagonal at β ≈ 0.2–0.5 (leading-Regge kernel for helicity-changing structures) may differ from the extrapolation; (ii) an additional k = 2 structure whose off-diagonal cancels Ĉ₂′'s.
>
> * **Update (literature on helicity mixing; LIT_HELICITY_MIXING.md):** No published work derives graviton-coupling bounds from spinning CFT sum rules, and none treats a curvature-induced helicity off-diagonal. CLS-D derive no bounds and defer the Regge/large-impact-parameter regime (§2.2.2); positivity is stated as PSD over 3-pt structures but never tested. In flat D = 4 (CHLPS 2201.06602) the rule Ĉ₂ maps to (B⁽¹⁾_k) is helicity-diagonal; c⁺⁺c⁺⁻ enters only a separate fixed-t helicity-flip rule B⁽³⁾. **The obstruction is the AdS twin of the known flat D = 4 IR band**: their IR cutoff makes the functional negative for m_IR^{−2/3} ≲ b ≲ m_IR^{−1} (their (3.6), Fig. 5), which with m_IR → 1/R_AdS is exactly our b̂ ≳ Λ^{2/3}; they accept it (CHMRS: "only rigorous if one assumes additional information about the spectrum at large impact parameters"). Known d = 3 bounds on t₄ are parametric only (CEMZ, CHP, Meltzer–Perlmutter, AKT 1805.07393: t₄ ≲ ln Δ/Δ⁴); the sharp flat expectation |ĝ₃|²M⁸ ≤ 24.9 log(M/m_IR) − 27.6 matches the form of Result D's bound. **Honest fallback statement: Result D conditional on the heavy spectrum in the band b̂ ∈ [~0.65Λ^{2/3}, O(Λ)]**, the same kind of condition as the flat D = 4 literature. Only literature-backed fix: a new (helicity-flip, B⁽³⁾-type) CFT sum rule outside CLS-D's family.
>
> * **Update (conditional Result D; RESULT_D_CONDITIONAL.md):** **Theorem D_cond** [D given the review's §2 hypotheses, per-state positivity below a window edge b̂_w, and the new hypothesis H_band(ε): states with b̂ ≥ b̂_w contribute ≥ −ε·16L_A]: √(t₄²+α_T²) ≤ 240·√((1+ε)·16L_A)/Λ⁴. With NF02 (ε = 0): **2.21·10⁻⁴ / 3.13·10⁻⁶ / 2.89·10⁻⁸ at Δ_g = 100 / 300 / 1000**, band edge b̂_c ≈ 16 / 30 / 62–68 (≈ 0.62–0.76 Λ^{2/3}) [M]; at 300 the computed graviton corrections (888 states) give the first negative state at b̂ = 30.0, the same edge. In the band the diagonal stays positive; the failure is purely helicity off-diagonal. Trade-off (re-weighted NF02 family): edge 57.5 → 3.38·10⁻⁶ at 300 (positive on all 888 computed states; band then starts at m ≳ 1.2Λ, J ≳ 30); edge 117 → 3.26·10⁻⁸ at 1000 (keeps all m ≤ 1.05Λ out of the band); no window beyond ≈ 2–2.7Λ^{2/3} is feasible. H_band(0) ⇔ small helicity coherence κ of the heavy spectral matrix in the band (κ_crit ≈ 0.6–0.75 at the edge, ≈ 0.1 at 2b̂_c); plausible if the band is dominated by Regge-factorised graviton exchange (flip cost (bM)⁻⁴) and for string-like spectra; not implied by the axioms [I]. **Flat comparison:** route A′ is exactly the CHLPS-form flat bound of NF02 with m_IR = 3/(2R_AdS): |ĝ₃|²Λ⁸ ≤ 3022.34·log(2Λ/3) − 5327.21 [D]; the log coefficient is 121× CHLPS's SDPB optimum 24.9 (×33.6 from having no forward terms or helicity-flip rule, ×2.6 from the margin, ×1.3 from the large-b̂ proxy), so the bound is ≈ 10× weaker than the (asserted, not derived) transplanted CHLPS optimum. Parametrically √(ln Δ_g)/Δ_g⁴, stronger than AKT's ln Δ_g/Δ_g⁴. Open tension at Δ_g = 100: model puts near-gap J = 10 at 1.1Λ in the band (−0.61) while the near-gap chain has +0.50. All band numbers use the extrapolated g; the exact kernel is being computed.
>
> The tables below still describe the B′-uplift functional and will be rewritten at the next batched review.

Sources: `STATUS.md` (chain), `OMEGA_U.md`, `OMEGA_NF.md`, `SUPERCONVERGENT_UPLIFT.md`, `FLAT_UNSTRIPPED.md`, `ROBUST_FUNCTIONAL.md`,
`NOFORWARD_FUNCTIONAL.md`, `ADS4_LIGHT.md`, `SPINNING_4PT.md`, `SPINNING_TTTT.md`, `SPINNING_TTTT_EXT.md`, `EXACT_SKELETON_NEARGAP.md`,
`CERTIFIED_NEARGAP.md`, `DIRECT_GRAVITON_NEARGAP.md`, `PARITY_ODD_HEAVY.md`, `OBJECT_IDENTITY.md`, `FUNCTIONAL_REDEFINITION.md`,
`SPINNING_MELLIN.md`, `DICTIONARY_T4.md`, and their referee reports. Nothing here is new computation.

Labels: **[certified]** computer-assisted proof (Arb/exact); **[D]** derived; **[DVE]** derived and checked end to end on exact
blocks; **[NE]** numerical evidence; **[M]** a modelling convention; **[I]** inferred/assumed; **[A]** standard assumption; **[Q]** quoted.

## 1. Setting

A unitary CFT₃ (parity-even or parity-violating) with stress tensor `T`, central charge `c_T`, and a **twist gap**: every single-trace
operator of spin `J > 2` has twist `Δ − J ≥ Δ_g − 2` (a dimension gap is not enough: it allows states below `Λ` at large `J`). No other
light single traces (H11; note this excludes all known holographic CFT₃s, ABJM included, which have light Kaluza–Klein towers). Large `c_T`
with `c_T ≫ Δ_g²` (tree level). Write `⟨TTT⟩` in the basis `(c_T, t₄, α_T)`; in the bulk `t₄² + α_T² = 57600 |ĝ₃|²/L⁸`
(`DICTIONARY_T4.md`; reproduced by the foundation referee, including the ECG and quasi-topological cross-checks). `Λ = √((Δ_g−2)(Δ_g−1))`.

## 2. Claim D (conditional)

Under hypotheses H0–H17 below, with a heavy-positivity threshold `n`,

`|ĝ₃|² (Λ/n)⁸ ≤ A log(Λ/n) + B_AdS + O(n/Λ)`,

i.e. `√(t₄² + α_T²) ≤ 240 √(A log(Λ/n) + B_AdS) · (n/Λ)⁴` (the tables use the exact `p`-integral, not the log form).

| version | functional | `A`, `B_AdS` | `n` | `√(t₄²+α_T²)` at `Δ_g` = 100 / 300 / 1000 | hypotheses |
|---|---|---|---|---|---|
| **D-safe** | NF02-based, no forward terms; even `k` := `C`-type dispersive functionals (`FUNCTIONAL_REDEFINITION.md`), odd `k` := `Ψ_k/N_k` | 3022.3, −6552.5 | **1 in the scalar-skeleton model** (certified on a grid at Δ_g = 100; scanned at 300, 1000) | 2.21·10⁻⁴ / 3.13·10⁻⁶ / 2.89·10⁻⁸ | H0–H4, H6–H8, H9b, H10–H17 |
| **D-strong** | `Ω^U` on M02 | 51.74, −86.30 | 1.08 / 1.05 / 1.01 in the same model (object-sensitivity range 1.00–1.10 at Δ_g = 100) | 4.2·10⁻⁵ / 5.3·10⁻⁷ / 4.1·10⁻⁹ (up to ×1.4 object/parity sensitivity) | as D-safe, plus H9a |

**What `n = 1` means.** It is the threshold of the near-gap heavy action computed with (i) the flat bulk-point density, which carries
≈ 98 % of the action near the gap, plus (ii) the exact scalar-skeleton remainder (≤ 0.02 `N_s`). Item (ii) is exact. Item (i) for gravitons
is H2/H4: at leading order it is the CLS-D scattering-crystal theorem (arXiv:2311.04271) for route-A′ functionals (H0); the
load-bearing window is `x = ν/m ∈ [0.3, 0.7]` (not `x → 1`), where the `O(1/m)` graviton correction is not computed. So `n = 1` is
**not** shown for gravitons.

Meaningful only for `Δ_g ≳ 100`. Parametrically `|t₄|, |α_T| ≲ √(log Δ_g)/Δ_g⁴`. **Novelty:** the scaling and the AdS₄/CFT₃ translation
with `log(M R_AdS)` are in arXiv:2201.06602 §5 (CLPS), whose footnote functional is of the same forward-free `B₂ + B₃` class as NF02;
D-safe's `A` is ≈ 121× CLPS's optimum (≈ 11× weaker bound). arXiv:2311.04271 lists Regge positivity and the light-side dictionary as
open steps; this chain assumes rather than closes them. The referee found no paper deriving an explicit gap bound on ⟨TTT⟩ in CFT₃ with
explicit constants; the explicit conditional bound and the printed 240 dictionary coefficient appear to be new.

## 3. What is established

### 3a. Flat space and kinematics (unconditional or certified)

| ingredient | status | file |
|---|---|---|
| Flat `F_{g3}` positive on all heavy states; flat bound `37.84 log(M/m_IR) − 45.21 + Δ_IR` | **[certified]**, 3 referee rounds | `FLAT_GRAVITON_BOUND.md` |
| JHEP-printed `F_{g3}` not positive (`−6281/315`); v1 (3.3a) misses `5B₅|₀` | **[D]** | same |
| Robust flat functionals M01, M02 (margin 9/100, 9/50) | **[certified]** | `ROBUST_FUNCTIONAL.md` |
| Unstripped and stripped fixed-`u` sum rules span the same space | **[D]** | `FLAT_UNSTRIPPED.md` |
| No forward-free functional uses smeared `B_{k≥4}`; NF02 from `B₂, B₃` | **[D]**; NF02 **[certified]** as a flat functional except one `+−` corner (float) | `NOFORWARD_FUNCTIONAL.md` |
| Uniform positivity of the scalar-kernel `B₂` part (all `M, β`) | **[certified]** | `UNIFORM_POSITIVITY.md` |
| Even-kernel circle-mode weights `(1/36, 1/16, 1)`; `K*`, `𝔖_ff` flip-free; `N₄^{K*}`, `N₆^{K*}` | **[D]**, reproduced by referee | `SPINNING_TTTT_EXT.md` §6, `EXACT_SKELETON_NEARGAP.md` §9 |
| Parity-odd heavy exchanges: no `λ_Eλ_O` interference (parity invariance of all kernels) | **[D]**, reproduced by referee | `PARITY_ODD_HEAVY.md` |

### 3b. Scalar-skeleton model (exact or certified, but scalar)

| ingredient | status | file |
|---|---|---|
| Exact scalar-skeleton moments `C_{k,ν}` near the gap (Mellin/Mack + Arb Barnes lemma) | **[D]**; independently reproduced by the referee to ≤ 7·10⁻⁶ (Δ_g=100), 5·10⁻⁸ (300), 5·10⁻⁷ (1000) | `EXACT_SKELETON_NEARGAP.md` |
| D-safe action positive on the Δ_g = 100 grid (m/Λ ∈ [1.01, 3], J ≤ 10, both channels; Δφ = 3 and 2), 64-node `p`-rule | **[certified]** for that model (384/384 cases for the 64-node rule; min F/N_s +0.227 for Δφ=3, +0.072 for Δφ=2). **Exact `p`-integral certified** (Bernstein-ellipse Gauss–Legendre error bound) on 400/402 cases incl. J = 11–16 at m/Λ ∈ [1.12, 1.3]; fails only at ++, J=0, m = 1.01Λ (no admissible ellipse; float value +0.997 N_s). Continuity in `m` between grid points [NE] (~30 core-h to certify) | `CERTIFIED_NEARGAP.md` |
| Even-`k` objects `C_{k,ν}` meet: same leading action and B′ density, Einstein-free for `k ≥ 4`, no forward terms, Mellin-exact | **[D]** for external scalars | `FUNCTIONAL_REDEFINITION.md` |

### 3c. Spinning checks (scope-restricted)

| ingredient | status | file |
|---|---|---|
| Heavy formula `C^± = (2P/m^{2k})(|c_same|² ± |c_opp|²)(1+O(1/m²))`, `⟨JJJJ⟩` | **[DVE]** `k = 4–7` | `SPINNING_4PT.md`, `SPINNING_ODDK.md` |
| Heavy formula, `⟨TTTT⟩` | **[DVE] only for β = 1, ν ≤ 3, even J, parity-even exchanges, k ≥ 4, kernels `K*`/`𝔄`**; `k = 2` [NE + I cleaning]; `𝔖_ff` unresolved at low J for k ≤ 4; ¼ not tested end to end | `SPINNING_TTTT.md`, `SPINNING_TTTT_EXT.md`, `DIRECT_GRAVITON_NEARGAP.md` |

Conditional items previously listed here (graviton pole in `C⁺₂`; Einstein term absent from `Ψ_k`) are moved to the hypotheses (H6, H8, H16).

## 4. Hypotheses (graded; referee grades adopted)

| # | hypothesis | grade | evidence / what would settle it |
|---|---|---|---|
| **H0** | **Existence of the spinning (⟨TTTT⟩) functional as defined**: a `C`-type dispersive `Ĉ_k` (all `k = 2…7`, incl. `Ĉ₂`, which carries the graviton pole and the whole light side) with meromorphic weights and flat limit NF02 (merges the former H14′) | **[D-sketch], load-bearing** (`H0_SPINNING_FUNCTIONAL.md`): the object exists in the literature up to one step — Chang–Landau–Simmons-Duffin, arXiv:2311.04271, build spinning light-ray sum rules for `⟨TTTT⟩` with meromorphic subtraction factors and prove their bulk-point (“scattering crystal”) action equals the flat graviton sum rules [Q]. Route A′: `Ĉ_k` := CLS-D's functional with its own `f^±_{k;2}`, then the CHMRS η-inversion in ν. Its flat limit is NF02 exactly [Q + D]. The project's `Ψ_k` are CLS-D's functionals with the subtraction factor replaced by `A_k` (same for k = 3, 7; different for 4, 5, 6) [D]; with that replacement the flat limit is **not** positive (`+−`, m/M < 1.0003, J = 16–24, min −0.016 N_s; only reachable for Δ_g ≳ 6·10⁴) [NE] | M2 d = 3 helicity specialisation (3–5 d); M3 exact η-inversion, recompute k(0), A, B_AdS (~1 wk); M4 endpoint convergence incl. H15 (1–2 wk); M5 light side (2–3 wk) |
| H1 | Two-pair heavy formula for `⟨TTTT⟩` at all needed `k`, `ν ≲ m`, β, `J`; MFT constant ¼ | **[NE]** as used; [DVE] in the restricted scope of §3c | needs `ν ~ m` (= H2), β ≠ 1, odd `J` (= H12′), a clean `k = 2` production rerun, ¼ end to end. Finite-`m` sensitivity: scaling even moments by `1 + c₂/m²` (`c₂ ≤ 268`) moves no threshold [NE] |
| H2 | Flat bulk-point dictionary for gravitons | for route A′ at leading order in 1/m, fixed x < 1: **[Q]** (CLS-D), pending the d = 3 specialisation (M2). Crystal = CHMRS Regge point η = 1, r* = (1−s)/(1+s), s = √(1−x²) [DV: all 9 printed CLS-D values reproduced; scalar ratio of `OBJECT_IDENTITY.md` reproduced]. For the current `Ψ_k` the B′ density is **false** at x ≳ 0.5 for k = 4, 5, 6 (R₄ = (1−r²)²/(1+r²), R₆ = (1−r²)⁴/(1+r²), R₅ = (1+6r²+r⁴)/(1−r²)² diverges as x → 1) [D given CLS-D]; in the bulk-point-factor model n = 1 is unchanged either way (min F/N_s A′ 0.2267/0.2057/0.2022) [NE + M] | M2 |
| H3 | Odd-`k` moments | **[M]** (continuation/interpolation in `k`); `⟨TTTT⟩` k = 5, 7 under `𝔄⊗𝔄` [DVE] in the §3c scope | — |
| H4 | `O(1/m)` corrections to the graviton bulk-point density (formerly: scalar-skeleton remainder transfers to gravitons) | [I], **localised** [NE + M]: the load-bearing window is x ∈ [0.3, 0.7]; x > 0.7 carries ≤ 0.5 % of Σ\|density\| (tolerates ±22–180 % coherent error). On x > 0.3 a coherent density error of ≈ 4.5–8 % breaks n = 1; at Δ_g = 100 the graviton remainder may be −1.25× to 2.67× the scalar one. η-inversion shift ≤ 0.087 / 0.0088 / 0.0012 N_s (∝ 1/Λ²) [D-sketch + NE] | M6: 1/m fluctuation correction around the CLS-D crystal on x ∈ [0.3, 0.7] (2–4 wk); exact spinning Mellin only if insufficient |
| H5 | Near-gap tolerance convention | not needed for the skeleton part | — |
| H6 | (F) flat limit of spinning light side; (P) polynomial contacts; graviton pole `16πG/(ν²+9/4)` in `Ĉ₂` | [I] (scalar: [DV]/[Q]) | — |
| H7 | (R) graviton-pole residue exact: −Ĉ₂′\|_E = 16πG/(ν²+9/4), Ĉ₃\|_E = 0 | For the spinning de-shift **Ĉ₂′** (H7_I14_LIGHTSIDE.md §18–§22): (HN) holds by construction, so the pole is **[P given Lemma R, Q3, Q4, D2, and the Regge-kernel derivation [D]]**; h₀ ≡ 0 [P given kernel]. For the old scalar-g Ĉ₂, (HN) fails (ε = O(1) at ν ≲ 8) — superseded. Ĉ₃\|_E = 0 [P given Q2, Q8, F]; no Π₂ Regge moment [D + NE] | — |
| H8 | `J₀(G^E) ≤ 2` | [Q] | KKSDZ |
| H9a | Forward-term double-trace/contact content (D-strong only) | [I], leverage up to 10¹³ | — |
| **H9b** | (P) for the smeared even replacements in D-safe | **[I]**, leverage up to 200 (Δ_g = 100) → 1.1·10⁴ (Δ_g = 1000); **analysable under route A′** (meromorphic weights) | M3/M5 of `H0_SPINNING_FUNCTIONAL.md` |
| H10 | EFT power counting of contacts | [A] | — |
| H11 | No other light single traces; `c_T ≫ Δ_g²`; double zeros of spinning kernels at double twists | [A]/[I]; excludes known holographic CFT₃s | — |
| H12′ | Heavy formula (and no `O(1/m)`) for parity-odd exchanges, incl. every odd-`J` `+−` state (the D-safe minimum at all three `Δ_g` is the `+−`, `J = 5` parity-odd state; 0.227 / 0.206 / 0.202) | **shape [DVE]** (`PARITY_ODD_HEAVY.md` §4): parity-odd `⟨TTO⟩` built in the Taylor-mode ⟨TTTT⟩ code (12-dim odd (2,2,J) space; CHL odd structure decomposed to 1e-59); at β = 1, J = 35 (m = 68.1) and J = 19 (m = 37.4): mixed waves W[E,O] vanish identically (pointwise parity-odd to 2e-12), and W[O,O]/W[E,E] is one constant for all k = 2…7, ν ≤ 3 and all kernels (spread ≤ 2e-5 / 2e-4 at k ≥ 5), i.e. the odd exchange has the opposite-helicity heavy-formula shape incl. its 1/m² corrections. **Normalisation N = 3.30 / 3.77 relative to CHL-normalised structures [I], unexplained** (expected 1; varies 14 % between the two masses); positive, so cannot change a sign; if its variation hid an O(1/m) term it would be |c₁| ≈ 14, below the tolerance on odd-J `+−` states (D-safe |c₁| ≤ 56 / 140 / 455 at Δ_g = 100/300/1000; D-strong ≤ 182 at 100) — no threshold moves [NE]. Not run: even-J parity-odd (`++`), other β | parity-odd MFT (Fock) data to fix N; even-J odd structure; other β; small-`m` odd waves at k = 3, 4 |
| H13′ | Heavy positivity off the scanned/certified set | certified grid at Δ_g = 100 incl. the exact `p`-integral (§3b; except m = 1.01Λ, J=0) and J ≤ 16 near the gap. D-safe at large `b̂` (`NF02_LARGE_B.md`): hybrid action (flat + leading AdS Regge) positive for m ∈ [1.5, 10]Λ, b̂ ∈ [20, 200] (Δ_g = 300, 1000; b̂ ∈ [40, 200] at 100) and b̂ ≤ 1000 at sample masses, min F/N_s 0.93 [NE]; exact skeleton remainder on 56 states up to b̂ ≈ 100, |S|/F_hyb ≤ 131/m², min F/N_s 0.95 [NE]; large-b̂ two-term envelope (AdS term A₁M₁/Λ² dominates for b̂ ≳ 4Λ^{2/3}) [D + NE]; **[I] (calibrated) for the remainder at b̂ > 100 and off the sample**; [NE] between grid points, m ∈ [Λ, 1.01Λ], J ≥ 18 near the gap | certify M₁ > 0 (interval arithmetic); derive the envelope remainder; exact remainder at J ≥ 140 (chunked Mack matrices); continuity in `m` |
| H14 | Object identity (`C` vs `Π̂` vs `Ψ_k/N_k`) | for D-safe even `k`: by definition (redefinition) at scalar level, else part of H0; for odd `k` and D-strong: [I], bounded in a bulk-point model [M] | `OBJECT_IDENTITY.md` |
| **H15** | `|ĝ₃|²` absent from the light side of `Ψ₅, Ψ₆, Ψ₇` | [I] | — |
| **H16** | I14 convergence/swappability of the light-ray sum rules; P1(`j′ = 0`); the untracked Jacobian `j_k(ν)`; Einstein term absent from `Ψ_k` | **[P given K1 (+ uniformity K1′), K2, K3]** for route A′: K1 = KKSDZ Regge bound J₀ ≤ 1 [Q, with their loophole]; K2 large-η decay η^{−1/2}/log η of the *total* de-shifted Ĉ₂ kernel [P given K1′] (scalar-type), [D] (spinning components); endpoint/edge rows and K3 per-family Rindler positivity [D]. j′ = 0 within the proven bound. P1 not used; j_k absorbed in N̂_k; Einstein-absence = Ĉ₃\|_E = 0 (H7_I14_LIGHTSIDE.md §7, §14) | — |
| **H17** | `k_AdS − k(0)` correction to the `|ĝ₃|²` coefficient | [D-sketch]; ≈ 4 % at Δ_g = 100 | — |

## 5. Where it can fail (in order of weight)

1. **H0.** The graviton functional the claim uses is not constructed; its flat limit, light side and contact content are inferred.
2. **H2/H4 at `x ~ 1`.** Near the gap almost all of the action is the flat dictionary at `x → 1`; gravitons are tested only at `x ≲ 0.1`.
3. **H12′.** The binding D-safe state is a parity-odd exchange never tested end to end.
4. **H13′ coverage** for D-safe at large `b̂` and larger `Δ_g`.
5. **H9b** (D-safe) and **H9a** (D-strong) contact content.

## 6. Shortest route (from `REFEREE_FOUNDATION.md` §5, adopted)

1. **Fix the object (H0).** Construct the spinning `C`-type `Ĉ_k` (even `k`, `K*`), or keep `Ψ_k` and recompute its flat limit, certify it,
   and recompute `k(0)`, `A`, `B_AdS` and contacts.
2. **H12′.** One end-to-end parity-odd `+−`, odd-`J` `⟨TTO⟩` point.
3. **H13′ for NF02.** Exact-skeleton scans at Δ_g = 300, 1000, `m ∈ [1.5, 10]Λ`, `b̂ ∈ [20, 200]`, plus a large-`b̂` envelope; extend the
   certificate (continuity in `m`; quadrature error now certified at Δ_g = 100).
4. **H9b.** (P) for the even replacements (becomes analysable with step 1).
5. **H2/H4 at `x ~ 1`.** The real near-gap content for gravitons.

No unified manuscript is written: Claim D is not yet a theorem.
