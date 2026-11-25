# CHANGES v16 (25 Nov 2026)

* T0 (T0_CERT.md, independently confirmed in INDEP_CHECK_T0.md) folded in.
  * **Abstract:** "check numerically in AdS" becomes "prove ... at leading order in the Regge expansion, in anti-de Sitter space".
  * **Assumption (iii):** now states the certified leading-Regge floor (≥ 1/5 of N_s, every heavy state, every Δ_g ≥ 300). For the cubic functional, (iii) only asks that beyond-Regge corrections stay within that margin.
  * **"Is the functional positive?":** the second certified statement is now the uniform floor 0.1999 (Δ_g ≥ 1000) / 0.199 (≥ 300), with a one-line proof description.
  * **Supplement:** new subsection "T0: a certified leading-Regge floor on every heavy state", with its scope (target weights; the OO conditions; nothing beyond leading Regge).

## v16b (26 Nov 2026)
* **Abstract rewritten for impact.** It adds "with a certified margin" (T0) and one closing sentence on the angular route ("two-percent bound at every spin ... given one light-ray identity"). **The letter body does not yet contain the angular-route paragraph**; it must be added before submission (pending ANGULAR_MATCH).
* **Cover letter.** The "why PRL" paragraph is rewritten (COVER_LETTER.md/.txt).
* **Length.** The abstract is 800 characters; check the PRL length limit (a 560-character variant is in the chat log).

## v16c (26 Nov 2026): scope edits
* **Introduction.** Adds the physical observable up front: stress-tensor energy flux isotropic to 1e-8, where causality allows O(1). The positivity sentence is updated for T0: the certified leading-Regge floor of one fifth of the flat margin.
* **"Why it matters".** Adds the bulk-language scale: |ĝ₃| = M_*⁻⁴ gives M_* ≥ Λ/120^{1/4} ≈ 0.30Λ.
* **Discussion.** Adds one outlook sentence: the ingredients are not special to d = 3; AdS₅/CFT₄ with N = 4 SYM checks is next.
* **Cover letter.** The method paragraph is updated for T0 and the bulk-scale sentence.

## v16d (26 Nov 2026): angular route in the body; length
* **Abstract.** Cut to 569 characters (the limit tracked in NOTES is 600). Its last sentence now matches the audited result: "Without that property, energy-detector algebra and crossing still give a two-percent bound at every spin."
* **Letter.** New paragraph "Without the positivity assumption.": τ ≤ (160/9)/(Δ_g−2)·[…]^{−1/2} = 0.0178 at 1000, with the inputs named.
* **Trims for length.** "Large spin" is compressed; the "Larger gaps" transfer sentences are compressed; one sentence about examined states is shortened. Estimated APS length ≈ 3670 (limit 3750).
* **Supplement.** New section "An argument without positivity: the angular route" (sum rule, exact ingredients, result, crossing proof of 𝒟 = 0, inputs), updated per ANGULAR_AUDIT.

## v16e (26 Nov 2026): batch 4, free-field check, all-J proofs
* **Batch 4 (BATCH4_FOLD.md).** All `% BATCH4` markers are resolved: supplement summary item (2), Table `tab:gap2000` (J = 12, 16, 20, 21, 23, 27, 29 from direct rows), the *Other near-gap states* paragraph and the coverage-table rows; the letter's "Larger gaps" sentence (22 spins, 7 direct, every computed state closes except the two lightest of J = 29).
* **Supplement, angular section.** Adds the free-field check (endpoint rule 3/2 + 3/2 = 3 exactly) and a new subsection "The finite-c_T spectral estimate" (ANGULAR_PROOFS: all-J collinear identities, kernel envelopes, two-scale crossing lemma). The duplicate formula statement is replaced by a cross-reference.
* **Length.** Estimated APS length ≈ 3670.

## v17 (26 Nov 2026): fixes from REFEREE_REPORT_6 (fresh-eyes hostile read; verdict "major revision")
* **Abstract (M1).** Rewritten to make conditionality explicit ("Given a positivity property that we prove in flat space and certify at leading Regge order in AdS …"), to add "at leading order in 1/c_T", and to replace the ambiguous "two-percent bound at every spin" with "bound the deviation by 0.018". 545 characters. ABJM is dropped from the abstract (P12).
* **T0 wording (M2).** One consistent statement: ≥ 0.1999 N_s (0.199 for Δ_g ≥ 300), curvature costs < 0.1% of the flat margin, for the exact de-shifted weights; the OO block and the published-weight scope are stated.
* **Supplement vs T0 (M3).** The S7 intro, S7.C, P0, the Table S9 reading, Table S11 and Table S17 now cite T0.
* **Angular route (M4).** New Result item (3) with Eq. (angular). The Letter paragraph "Without assumption (iii)" now lists the same inputs as the Supplement: the Ward and conservation identities as operator identities (ladder checked to J = 40), existence and light-cone convergence of the mixed commutator, no scalar with Δ ≤ 1, plus (i) and (ii) at leading order in 1/c_T; it notes the non-optimal 1/Δ_g scaling. Added to the Intro, the Discussion and Fig. 1(a) (dotted curve).
* **Implication (1) (M5).** "Derived for parity-invariant theories; for parity-violating ones one step … sketched."
* **Large spin (M6).** "Two conditional statements", with premises named: decay [A], model magnitudes, and the edge integral [M]; "less than half of the margin" is attributed correctly.
* **Assumptions (M7).** "Three assumptions, together with standard inputs graded in the SM, the strongest being … no conspiracy (quartic only)."
* **Larger gaps (M8).** The uncovered band 33 ≤ J ≲ 0.026Λ is now stated in the text and in the Fig. 1 caption. The Fig. 2 caption claim about Δ_g⁻² shrinking is removed.
* **S1 legend (M9).** Regenerated from the current Letter. It maps T0, the angular route, the new positivity, large-spin and larger-gaps wording, and the parity caveat; stale phrases are removed.
* **Other fixes.**
  * Script names are removed from the captions (P9).
  * The Fig. 2 caption now says "medium blue", not a dashed line, and mentions the T0 floor (P7).
  * "Conformal-collider bound, which holds without a gap" replaces "causality alone" (P3).
  * Result (2a) wording: "weakens it by" (P4).
  * M_* "0.3 of the gap" (P14).
  * The ⟨TTT⟩ free-field statement is restricted to the parity-even part (T13).
  * The certified τ₀ ≤ 2.8914e-8 is quoted (T4).
  * The quartic scan is "consistent but pre-asymptotic" (T11).
  * Supplement: the garbled "give hence" (T1), the I3 pointer (T2), the free-field Δ_φ remark (T12), "hostile audit" → "separate check", the date stamp is removed, "byte for byte" is reworded (P10), and H17 value [D] vs asymmetry band [D-sketch+M] (T9).
* **Not done (author decisions):**
  * M10: publishing the code and notes (`code/t0cert/`, the angular proofs) in a citable archive.
  * The broader jargon rewrite of S6–S8 (P9).
  * The selection-bias sentence in the Letter (P11): length.
* **Length.** Estimated ≈ 3610–3690 APS words (limit 3750).
