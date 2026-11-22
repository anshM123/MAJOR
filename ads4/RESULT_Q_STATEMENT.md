# Result Q: current statement (post-review, 2026-11)

Authoritative. Supersedes the headline claims of RESULT_Q_STAGE0–7, BELOW_FOLD2/3 and PARITY_ODD_CLOSE wherever they differ. Review: REVIEW_RESULTQ_2026-11.md, verdict **major revision**; nothing fatal found for the upper bound at ĝ₃ = 0.

> **Revision update (RESULT_Q_REVISION.md):**
> - **H_smear^Q → [NE + M].** The smeared corrections at the binding states are ≤ 3·10⁻⁵ of B₂, against a tolerance of about 24 %. The minimum FF/N is unchanged at +0.1427. The edge window x > 0.8 and k = 3 below the fold are still [I], but bounded by tolerances.
> - **V₄ → V₄′ (the single load-bearing input of the upper bound) [D-sketch].** V₄′ states absolute convergence at the corner where the light-cone edge meets the Regge corner, and Mellin equivalence of the spinning functional. The infinite [TT]₀ J-sum is never performed. In the scalar Mellin proxy the O(G) residues cancel identically [DV].
> - **The W³W³ residue in ∂Ĉ₄(0) is zero [D-sketch].** The u-channel is annihilated, and the s- and t-channel light sides are polynomial in ν². This assumes a spinning analogue of CHMRS's light-side formula, which is a new hypothesis.
> - **Cubic terms restored at grade [D-sketch + M].** κ = 0.110, and the lower slope is 0.49 (0.255 if double-counted), so ≥ ¼ either way.
> - **Implied cubic bound.** Result Q implies |ĝ₃|²Λ⁸ ≤ 120–198 (Δ_g = 300) and 154–257 (Δ_g = 1000). This is a *conditional strengthening* of Result D (1.07·10⁴ / 1.43·10⁴), not a replacement: Result D needs none of Q's extra hypotheses.

> **V₄′ test (RESULT_Q_V4PRIME.md):**
> - **V₄′ is false for the raw per-ν Ψ₄,ν, already at O(1).** The u-channel identity ⟨T₄T₂⟩⟨T₃T₁⟩ in MFT gives a contact term on f₄'s 1/x₁₃² pole that grows like ε^{−1/4}, and its value depends on the regularisation prescription. Result Q does not use the raw functional.
> - **For the de-shifted Ĉ₄^{ds}, the 3-unit contour shift makes it go like ε^{5/4} → 0 [D-sketch].** The identity and the J-resummed [TT]₀ term cancel exactly [P].
> - **Ψ₅ (k = 5) has no such problem [P / D + NE].**
> - **O(G) corrections are corner-convergent by power counting [D-sketch].** The exact Einstein cancellation in the swapped form is [I].
> - **V₄′ is replaced by V₄″:** on tree-level ⟨TTTT⟩ the total de-shifted kernels, with the ν-integral inside, converge absolutely at the corner, and the swapped form equals the commutator form there. Grade [D-sketch], with the O(1) regular part [P]. Hypothesis 5 (the de-shift exists and is analytic in |Im ν| ≤ 3) is now load-bearing for the light side too.

> **V₄″ test (RESULT_Q_V4DPRIME.md): V₄″ is false [DV]; V₄‴ replaces it.**
> - **The cure fails.** With the actual de-shift weights, the u-channel contact term still diverges (exponent −0.15 to −0.23), and its size depends on the ν-truncation. There are two separate obstructions:
>   - the ν-integrand grows like e^{πν/4};
>   - the two branches do not join at ν = 0.
> - **Cone-preserving regularisation removes it.** With the cone-preserving prescription t → t + iε the contact term is exactly 0. This prescription is the one consistent with L[T]|Ω⟩ = 0.
> - **V₄‴ [D-sketch]:** per ν (with the ν-integral outside), swapped form = commutator form at the corner. The O(1) part is exactly 0 under the light-ray prescription, and the O(G) part converges per ν with margin r^{1/2}. The O(G) homogeneity lemma is [D], checked on the exact structures.
> - **Two new hypotheses:**
>   - **H_LT:** the CLS-D functional uses the light-ray regularisation at the collinear point.
>   - **H_νout:** the per-ν light side, summed over families, is polynomially bounded in the heavy normalisation.
> - **Hypothesis 5 reworded:** "χR_j is analytic in |Im ν| < 3 up to edge poles" [NE]. χ itself has poles at ±i/2.
> - **Open:** heavy-side validity uses the ν-integral inside, while the light side now uses it outside. This circularity is concrete and unresolved.

> **Final light-side status (RESULT_Q_V4FINAL.md), superseding the V₄′/V₄″/V₄‴ boxes above:**
> - **Pointwise Lemma A [D].** For tree terms with J₀ < 4, the t-integrated, 1↔3-symmetrised integrand vanishes at every φ₁ ≠ φ₃. So the light side is the zero function in every order of integration, and the ν-circularity is resolved: Stages 3–7 computed ν-outside images on both sides, consistently.
> - **H_LT → definition D_LR [D].** CLS-D's domain and KS-D's contour fix the cone-preserving prescription. On the light-transform manifold, W_u·f₄ ≡ 0 [P]. The tube regularisation is not a regularisation of CLS-D's integral.
> - **H_νout deleted.** For J₀ ≥ 4 terms the light side is (P) and converges absolutely.
> - **The remaining interchange (heavy sum vs corner/ν) is H_K** [D-sketch], the same class as Result D's I14(b). It becomes [D] given H_ip (J₀ < ½ at the raw corner).
>
> **Final upper-bound hypotheses:**
> 1. A1, A2, a local EFT, and no conspiracy in the contact tower.
> 2. K1′/K2/K3 plus H_K.
> 3. Pointwise Lemma A with D_LR [D].
> 4. Spinning (F) + (P) [D-sketch].
> 5. The continuum de-shift exists [NE to t = 30; the C₅ spinning coefficients are I].
> 6. Inheritances from Result D.
> 7. H_smear^Q [NE + M].
> 8. H_fwd^Q.
> 9. H_OO^Q.
> 10. Accuracy.
>
> The numbers are unchanged.

> **Hostile review corrections (REVIEW_UNIFIED_2026-11.md), superseding conflicting text below:**
> - The bound at ĝ₃ = 0 is **sound with corrections**. The implication positivity ⇒ bound is [D].
> - **Light side:** Einstein gravity's contribution is zero at [D] (pointwise Lemma A + D_LR; no circularity remains). Two caveats:
>   - no spinning O(G) evaluation of Ψ₄ has been done (to do);
>   - the J₀ ≥ 4 part of the remainder (contacts **and heavy-exchange double traces**) needs a new named hypothesis **H_EFT [A]/[I]** (CHMRS (5.3)-type).
> - **ĝ₃ ≠ 0 requires κ ≥ 0 [D-sketch + M].** The bound is stated for κ ≥ 0 or t₄ = α_T = 0. If κ could be negative, Result D's |ĝ₃|² bound would allow a shift ~15× the bound.
> - **The implied cubic closure is withdrawn.**
>   - The "|ĝ₃|²Λ⁸ ≤ 120–257" values are dropped from claims.
>   - The Δ_g-independent lower slope 0.49 is 1.9× CHLPS's flat optimum of 0.26, which cannot survive the flat limit. This favours the σ_h·U double-counting reading.
>   - At most a remark remains (≈ 200/257 [D-sketch + M]); no ratio |ĝ₃|²/g₄ is claimed.
> - **c₁t₄** is carried explicitly: ε_{c₁} ≤ 0.57|c₁| (300) / 0.059|c₁| (1000) using Result D's t₄ bound. Importing D's bound makes Q inherit H_LJ.
> - **New and regraded items:**
>   - H_K is [D-sketch].
>   - **H_RS is [I]** (every s_d depends on it).
>   - **H_LJ^Q** is new, [M]: large-J tolerances ≥ 1.1·10⁴ c₂, a margin of ×10–10³.
>   - The lower bound reads −O(Λ⁻²) ≤ g̃₄Λ⁶/(8πG) [D-sketch].
> - **ABJM:** the arithmetic is confirmed and c_T does not enter the placement. ĝ₃ = 0 there follows from 10d IIA supersymmetry, not "maximal SUSY". Consistency check only (outside A1).
> - **Literature:** CLS-D 2311.04271 leaves the AdS uplift as future work, so Results D and Q are a first conditional execution of that programme; credit CLS-D and CHMRS 2106.10274. No prior ⟨TTTT⟩ sum-rule bound on R³/R⁴ couplings was found.

> **Forward channels and parity-odd J = 0, 2 measured (RESULT_Q_FWD_CLOSE.md):**
> - **H_OO^Q [I] → [NE+M].** OO J = 0, 2 (gap, 1.05Λ, 1.3Λ; Δ_g = 300, 1000): s_d = +0.509…+0.511, equal to ++ within 5·10⁻⁶ in every gauge. The helicity lemma holds for the k = 4 functional (order 0 exact to 10⁻⁹; beyond: ≤ 1.3·10⁻³ c₂, falling like Δ^−2.3…−2.5).
> - **H_fwd^Q, channel coverage [I] → [NE+M].**
>   - +− J = 6 (the last channel with failure threshold within 4× the measured size) has s_d = −0.515…−0.564, which is favourable.
>   - +− J = 7, 8, 9, 11 are favourable. ++ J = 6, 8 are unfavourable but have thresholds ≥ 28.
>   - Caveat: at Δ_g = 300 there is a J-drift of s_d (−0.07 at J = 7, 2.4× the truncation error), cause open. J ≥ 7 is therefore quoted from Δ_g = 1000.
> - **Still assumed [I]:** real-section ≡ crystal-cycle equivalence for EE beyond LO (shared with Results D and Q). Also unmeasured: ++ J ≥ 10, +− J = 10 and ≥ 12, and states above 1.475Λ (all with thresholds ≥ 11 or (Λ/m)⁶-suppressed).
> - **Untouched:** the σ_h·U double-counting question (it decides the sign at δ = 0.05, not at δ = 0.2).
> - **Full-set minimum unchanged:** +0.1426 ± 0.0017 (δ = 0.2), +0.0060 ± 0.0013 (δ = 0.05). Binding pair: ++ J = 0 and OO J = 0 at the gap. **Numbers in the headline unchanged.**

## Headline (conditional)

At tree level, with twist gap Δ_g = 300 or 1000, and at ĝ₃ = 0:

    0 ≤ g̃₄ Λ⁶/(8πG) ≤ 17.24 · log(2Λ/3) − 19.22 (Δ_g = 300) / −19.27 (Δ_g = 1000)   (+ c₁ t₄ term)

- The functional is the δ = 0.2 functional, with A = 17.241.
- **Heavy positivity minimum:** +0.143 on the measured forward channels, with the binding state ++ J = 0 at the gap. The smeared pieces use the scalar skeleton plus leading Regge only, so graviton corrections to them are [I].
- **Status:** the Δ_g = 100 column is dropped, and the |ĝ₃|² terms (κ, the ¼ lower slope, 4πG|ĝ₃|²) are [I].
- **ABJM:** g₄Λ⁶/(8πG) = 2ζ(3) sits at 3.3 % / 2.6 % of the bound. This is a log-dominated consistency check on a theory outside A1, not evidence for the bound.

## Regradings applied

| claim | grade |
|---|---|
| Lemma A / V₄ (the swapped form and the infinite [TT]₀ tail are untested) | D-sketch |
| −Ĉ₄^{ds}\|_{R⁴} = 2g₄; the g₅ cancellation | D-sketch |
| κ, the ¼ slope, 4πG\|ĝ₃\|² in AdS (W³W³ is an exchange, not a contact) | I |
| validity J₀ < 7/2 of the de-shifted functionals | D-sketch |
| full-set minimum +0.143 | M (smeared) + NE (forward, 5 channels, given real-section ≡ crystal for EE) |
| ABJM placement | D-sketch given Q, outside A1 |

## Hypotheses

1. **Setup:** A1 (twist gap for all J) and A2 (tree level). The light sector is a local EFT, with no conspiracy in the contact tower.
2. **Regge:** the KKSDZ Regge bound plus K1′/K2/K3 for the de-shifted functionals.
3. **V₄:** including the [TT]₀ large-J tail; microcausality; J₀(E×W³) ≤ 3.
4. **Spinning functionals:** (F) + (P).
5. **De-shift:** the continuum C4/C5 de-shift solutions exist for all β.
6. **Inherited from Result D:** H7, EE admissibility, c₁, and the t₄ bound — or carry c₁t₄ explicitly.
7. **H_smear^Q:** graviton corrections to the smeared parts stay below about 24 % of B₂.
8. **H_fwd^Q:** covers the forward channels not measured, and real-section ≡ crystal for EE beyond leading order.
9. **H_OO^Q:** the parity-odd J = 0, 2 channels equal ++.
10. **Accuracy:** HYB accuracy elsewhere; β > 30.

## Literature

The novelty claim stands at [I] (based on an abstract-level search). The closest real precedent is CHMRS 2106.10274: two-sided log(Δ_gap) bounds on scalar contact couplings in AdS4 with gravity. Kundu 2104.11238 and the N = 6 bootstrap papers (Binder–Chester–Jerdee–Pufu 2011.05728) need to be discussed, not just cited.

## Next steps (cheapest first, largest effect on the grades)

1. The smeared k = 2 EE correction at ++ J = 0, 2, from the existing `bg_eeq` images (H_smear^Q).
2. The swapped-form Ψ₄ on MFT ⟨TTTT⟩, including the [TT]₀ large-J tail (V₄).
3. The W³W³ residue in ∂Ĉ₄(0) (the |ĝ₃|² sector).
