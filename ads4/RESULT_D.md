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
>   from a leading-Regge approximation to an exact residue. **A = 3022.3, B_AdS = −6552.5 and the bounds survive**; only k(0) is approximate
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
| H7 | (R) graviton-pole residue exact | [D-sketch] | — |
| H8 | `J₀(G^E) ≤ 2` | [Q] | KKSDZ |
| H9a | Forward-term double-trace/contact content (D-strong only) | [I], leverage up to 10¹³ | — |
| **H9b** | (P) for the smeared even replacements in D-safe | **[I]**, leverage up to 200 (Δ_g = 100) → 1.1·10⁴ (Δ_g = 1000); **analysable under route A′** (meromorphic weights) | M3/M5 of `H0_SPINNING_FUNCTIONAL.md` |
| H10 | EFT power counting of contacts | [A] | — |
| H11 | No other light single traces; `c_T ≫ Δ_g²`; double zeros of spinning kernels at double twists | [A]/[I]; excludes known holographic CFT₃s | — |
| H12′ | Heavy formula (and no `O(1/m)`) for parity-odd exchanges, incl. every odd-`J` `+−` state (the D-safe minimum at all three `Δ_g` is the `+−`, `J = 5` parity-odd state; 0.227 / 0.206 / 0.202) | **shape [DVE]** (`PARITY_ODD_HEAVY.md` §4): parity-odd `⟨TTO⟩` built in the Taylor-mode ⟨TTTT⟩ code (12-dim odd (2,2,J) space; CHL odd structure decomposed to 1e-59); at β = 1, J = 35 (m = 68.1) and J = 19 (m = 37.4): mixed waves W[E,O] vanish identically (pointwise parity-odd to 2e-12), and W[O,O]/W[E,E] is one constant for all k = 2…7, ν ≤ 3 and all kernels (spread ≤ 2e-5 / 2e-4 at k ≥ 5), i.e. the odd exchange has the opposite-helicity heavy-formula shape incl. its 1/m² corrections. **Normalisation N = 3.30 / 3.77 relative to CHL-normalised structures [I], unexplained** (expected 1; varies 14 % between the two masses); positive, so cannot change a sign; if its variation hid an O(1/m) term it would be |c₁| ≈ 14, below the tolerance on odd-J `+−` states (D-safe |c₁| ≤ 56 / 140 / 455 at Δ_g = 100/300/1000; D-strong ≤ 182 at 100) — no threshold moves [NE]. Not run: even-J parity-odd (`++`), other β | parity-odd MFT (Fock) data to fix N; even-J odd structure; other β; small-`m` odd waves at k = 3, 4 |
| H13′ | Heavy positivity off the scanned/certified set | certified grid at Δ_g = 100 incl. the exact `p`-integral (§3b; except m = 1.01Λ, J=0) and J ≤ 16 near the gap. D-safe at large `b̂` (`NF02_LARGE_B.md`): hybrid action (flat + leading AdS Regge) positive for m ∈ [1.5, 10]Λ, b̂ ∈ [20, 200] (Δ_g = 300, 1000; b̂ ∈ [40, 200] at 100) and b̂ ≤ 1000 at sample masses, min F/N_s 0.93 [NE]; exact skeleton remainder on 56 states up to b̂ ≈ 100, |S|/F_hyb ≤ 131/m², min F/N_s 0.95 [NE]; large-b̂ two-term envelope (AdS term A₁M₁/Λ² dominates for b̂ ≳ 4Λ^{2/3}) [D + NE]; **[I] (calibrated) for the remainder at b̂ > 100 and off the sample**; [NE] between grid points, m ∈ [Λ, 1.01Λ], J ≥ 18 near the gap | certify M₁ > 0 (interval arithmetic); derive the envelope remainder; exact remainder at J ≥ 140 (chunked Mack matrices); continuity in `m` |
| H14 | Object identity (`C` vs `Π̂` vs `Ψ_k/N_k`) | for D-safe even `k`: by definition (redefinition) at scalar level, else part of H0; for odd `k` and D-strong: [I], bounded in a bulk-point model [M] | `OBJECT_IDENTITY.md` |
| **H15** | `|ĝ₃|²` absent from the light side of `Ψ₅, Ψ₆, Ψ₇` | [I] | — |
| **H16** | I14 convergence/swappability of the light-ray sum rules; P1(`j′ = 0`); the untracked Jacobian `j_k(ν)`; Einstein term absent from `Ψ_k` (derived given H8, RB₂, microcausality, I14) | [I] | — |
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
