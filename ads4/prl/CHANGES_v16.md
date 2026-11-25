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
