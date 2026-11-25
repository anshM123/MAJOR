# CHANGES v13 → v14 (25 Nov 2026): c₁ = 0, batch 5, Ω^Q_up binding states, energy-flux corollary

Sources: `ads4/C1_FINAL.md`, `ads4/BATCH5_FOLD.md`, `ads4/QBIND_COROLLARY.md`, the newest update boxes at the top of
`ads4/UNIFIED_THEOREM.md`, and `ads4/INDEP_CHECK3.md` (independent check: reproduces c₁ = 0 and batch 5, refutes nothing; its one
grade refinement is adopted). Nothing was recomputed except the figures (same scripts) and the arithmetic checks listed in §5.
**Abstract and acknowledgments / AI-use disclosure byte-identical** (checked by diff against v13). Tone and structure of v13 kept.

## 1. Letter

| place | change |
|---|---|
| Setting and assumptions | "The quoted numbers are at c₁ = 0, where c₁ is one universal constant …" removed; the paragraph now ends with κ > 0 (or τ = 0) and the light-scalar condition. |
| Why it matters | three sentences added (QBIND_COROLLARY §2.5, adapted): energy flux (E/2π)[1 + A cos(4θ − φ)], A ≤ τ/4 [CDPrakash]; collider bounds require A ≤ 1, saturated by a free scalar or fermion; here A ≤ 7.2·10⁻⁹, flux isotropic to 1 part in 10⁸, ⟨TTT⟩ that of equal numbers n_s, n_f of free real scalars and Majorana fermions to \|n_s − n_f\| ≤ 7.2·10⁻⁹ c_T. "between the two" → "between the two bounds". |
| Fig. 1 caption | ", at c₁ = 0" removed; `fig/make_fig.py` annotation "(c₁ = 0, K = 1)" → "(K = 1)"; figure regenerated (numbers unchanged). |
| Functionals | "and a term 16πGc₁t₄. We derive that c₁ is a single universal number … sketch an argument that it vanishes under one further assumption; the bounds are exact functions of c₁" → "We derive that no term linear in t₄ or α_T survives: at k = 2 the weight is affine in the light-ray times, so the stress-tensor exchange enters only through the energy flux and its first moment; both are helicity-diagonal for Einstein gravity and helicity-flipping at linear order in t₄ and α_T, and the EE structure annihilates helicity flips. A direct numerical evaluation gives zero to 10⁻¹¹ of the Einstein term." |
| Result | "and for every value of c₁" removed; Eq. (2) is now τ ≤ τ₀ := 240Λ⁻⁴√(16L_A/K) (one display row, as before); "and x = 1.67·10⁻⁴\|c₁\|" removed; Eq. (3) loses the 5.946·10⁴c₁²/(κΛ⁴) term; "At Δ_g = 1000 and c₁ = 0" → "At Δ_g = 1000". The sentence "The dependence on c₁ is exact …; less than 1 % unless \|c₁\| exceeds about 60 (cubic) or 1300 (quartic)" removed. |
| Is the functional positive? | "near the gap on every other spin up to 16 and on several larger spins" → "near the gap on 25 of the 31 spins up to 32 that have near-gap states (all of them up to 21 and the odd spins from 23 to 31)". New sentence: "With error control, Ω^Q_up is positive on its binding states near the gap and at masses 1.3–3Λ." "Elsewhere, including Ω^Q_up outside masses 1.3–3Λ, …" → "Elsewhere the values are computed in models or interpolated." `% BATCH5` marker removed (resolved). The strict-protocol caveat on the +− J = 4 gap state is in the Supplement only. |
| Fig. 2 caption | "at J = 29 the lightest part (pale)" → "at J = 29 and 31". |
| Fig. 2 (`fig/make_fig2.py`) | NEPM adds J = 17, 18, 19, 20, 25 (twist edge to 1.05Λ); J = 31 twist edge → 1.04Λ drawn pale (signed window, like J = 29). The partial [NE±] mass ranges of the even J = 22–30 (held or λ-level, above [NE + M] twist edges) are deliberately not drawn: at this log scale the bars merge and would suggest covered twist edges. Docstring updated. |

`% BATCH4` (Letter, "Larger gaps") kept.

## 2. Supplement

| place | change |
|---|---|
| §S1 legend | [NE±] list → J = 0, 2, 4–21, 23, 25, 27, 29, 31 from the twist edge, thinnest/held/signed listed; "25 of the 31 spins" mapped; Ω^Q_up sentence mapped (binding pair, J = 2 and +− J = 5 minima [NE±] L ≥ 58; +− gap J = 4, 5, 6 [NE±] given K₂ edge rows, J = 4 [NE + M] strict); "no term linear in t₄ or α_T survives" → row L3 and §S5.E, **[D] given the [TT] double-zero statement of L1 (parity-odd, odd-twist [TT] double zeros [D-sketch])**; "direct numerical evaluation" = the numerical form of the same pointwise statement [DV] (per INDEP_CHECK3: not an independent route); energy-flux statement → §S5.F. "The quoted numbers at c₁ = 0" removed. "What is here" and "Independent checks" updated (fourth check). |
| Table S2 | L1 grade: "+ the parity-odd, odd-twist [TT] double zeros [D-sketch]" (INDEP_CHECK3; the same input carries the Einstein term). L3 rewritten: "c₁ = 0 [D] given L1's double zeros (parity-odd, odd-twist [D-sketch]) and L4; F₂ [D + DV]; O(α_T) step [D]; shock part [P]" (was "c₁ = 0 [D-sketch given F₂], F₂ [A]; value [I], carried exactly"). |
| §S4 Full statement | "for every value of c₁", x, c₁² terms removed; (1) τ ≤ τ₀ = 3.13206·10⁻⁶ / 2.89126·10⁻⁸, symmetric interval \|t₄\| ≤ τ₀; (2a) 72.0406 / 92.8135 with no c₁² term; (2b) 125.71 / 164.51. The "c₁ terms are exact … thresholds" paragraph replaced by the c₁ = 0 form of the two optimisations. Net grades: implication (1) "[D] … given the [TT] double-zero statement of L1 (its parity-odd, odd-twist part [D-sketch]), which carries both the Einstein term and c₁ = 0". I3: "c₁ form" → "no O(t₄) or O(α_T) term (c₁ = 0)". |
| §S5.E (was "The constant c₁ and the exact bounds") | rewritten as "The O(t₄) light term: c₁ = 0, and the bounds": the chain (i) L4, (ii) T block only / double zeros, (iii) F₂ in three steps — affine weight with the closed form [P]; E⁽¹⁾\|_T = ½{E, −i n̄·D}, Ward/rotation argument for α_T [D + DV]; EE contraction = identity [P] — conclusion [D], numerical form ≤ 1.75·10⁻¹¹ with O(1) controls, independent re-implementation to 6.9·10⁻²⁴; the old (a) shock, (b) k = 3, (c) ASD closure kept as consistency statements; b₊ remnant reconciled [D-sketch, not load-bearing]; "not done": an independent Lemma-R evaluation. Items (1), (2), Eq. (S-repair) as exact numbers; **one-line remark** on a hypothetical c₁ (59.6 / 1306 thresholds). **Table S3** (was "Exact bounds versus \|c₁\|") replaced by the grade table of the chain (numbering of all later tables unchanged). |
| new §S5.F "Energy-flux corollary" | QBIND_COROLLARY §2.5 paragraph, adapted (c₁ = 0): CDP Eq. (4.39) matrix [Q], 4θ rotation, A_ε ≤ τ/4 [D], Buchel et al. Eq. (3.8) 1/8 → 1/4 note, A ≤ 7.23·10⁻⁹, max/min ≤ 1 + 1.45·10⁻⁸, \|n_s − n_f\| ≤ 7.23·10⁻⁹ c_T, \|p_T\| ≤ 1.11·10⁻¹⁰ c_T; status of the cubic bound. |
| §S5.H light single traces (subsections after S5.E shift by one letter: de-shift S5.G, light S5.H, Λ-forms S5.I; all cross-references are `\ref`s) | "+16I_DΛ²c₁t₄" and "(plus the c₁ term)" removed. |
| Table S4 (Λ-forms) | columns x/\|c₁\| and c₁² coeff. removed. |
| §S6 "How it is computed" | direct rows now J = 7, 9–32; Gauss–Kronrod edges J = 21, 23, 25, 27, 29, 31. |
| Table S5 | J = 23/27/31 row rewritten (J = 23 L 1.18 → 1.62, no longer thin; J = 31 twist edge +0.920, L 0.90 sign-free / 1.14 signed; 1.04–1.10 L 1.35–2.13); **new row J = 17–20, 22, 24–26, 28, 30, 32** (twist-edge values/L, masses, held / λ-level suffixes, grades); "other" row trimmed (the superseded J = 20 Monte-Carlo edge and the J = 32 chain note removed); **new row Ω^Q_up at 1000, binding states** (QBIND §0: L = 62 (25), OO 58, J = 2 114, +− J = 5 90, +− gap J = 4/5/6 L = 71/148/108 given K₂ edge rows, strict 0.4/1.2/1.1 → J = 4 [NE + M] strict). Table split into three parts (float height). |
| §S6 new paragraph "Direct rows for J = 17–32" | structure-basis repair (anchored O(2) angle + reflection; residuals; unfitted node; discarded basis artefact; independently reproduced), J = 17, 19 own-row route (ill-conditioned shape basis), outlier-filter exemption, reference-zero rule, failed-node rule, J = 20 / 28 shifts, the not-used odd-J edge envelope for the even twist edges. |
| Table S6 | 1000 row: direct-row minima for J = 17–32 added; the old J = 23 and J = 31 direct values moved there. |
| Table S7 | H_M (no longer lists J = 17–32), H_edge (from direct nodes for J ≤ 32; the flat edge term sets the even twist-edge grade), H_low (J = 12, 16–32), **new H_hold [M]** (own-row hold for J = 11, 13, 15, 17, 19, 22, 26, 30). Split into two parts (float height). |
| Table S11 | [NE±] row and [NE + M] row rewritten for batch 5 (25 of 31; λ-level, signed, held lists); Theorem-W row: "J = 26–32 direct rows and error budgets". |
| Table S17 | J = 17, 19 → [NE±] at 1000; J = 23–31 row: "29, 31 signed at the twist edge"; "other near-gap states" row with the batch-5 grades. All `% BATCH5` markers resolved; `% BATCH4` kept (4). |
| §S8.A | "(K = 1, c₁ = 0)" → "(K = 1)". |
| §S10 / Table S18 | four rows from INDEP_CHECK3 (c₁ steps; pointwise zero 6.9·10⁻²⁴; basis repair; J = 18, 23, 25 and the count — our fold chain re-run, not re-implemented); row pointers ("last four/five rows") corrected; the not-reproduced list adds the fold protocol of J = 17–32; the Barnes map / ASD solution now "consistency statements". |
| §S12 Scope | "a value of c₁ without F₂" removed. |

## 3. Grades (no upgrade beyond the deliverables)

* c₁ = 0: [D-sketch given F₂ [A]] → **[D] given L1's double zeros (parity-odd, odd-twist [TT] double zeros [D-sketch]) and L4**
  (C1_FINAL's [D], with INDEP_CHECK3's refinement as instructed). F₂ [A] → [D + DV].
* **New caveat on L1** (INDEP_CHECK3 §1.5; H0 §15.1): the parity-odd, odd-twist [TT] double zeros are [D-sketch]. This input also
  carries the Einstein term, so the net grade of implication (1) is now stated "given" it. *Flag for the authors:* strictly, by the
  weakest-link rule, this makes the parity-odd-double-trace part of item (1) [D-sketch]; the Letter's "we derive" is kept, as
  instructed ("the Letter can still say c₁ = 0 in plain words").
* Batch 5: exactly BATCH5_FOLD's grades (twist edges of even J = 22–32 stay [NE + M]; held / λ-level / signed suffixes kept).
* Ω^Q_up binding: [NE + M] given H_RS → [NE±] with H_RS as an input; +− J = 4 gap [NE + M] under the strict protocol.
* Corollary: [D] from the cubic bound, kinematics [Q] checked [DV]; status of the cubic bound.

## 4. Length and build

`wc_prl8.py letter.tex 3 0.9375 1.03846` (acknowledgments excluded): prose 2938, inline formulas 232 (conservative 295), 3 display rows
(48), figures 344: **TOTAL 3562 (3625 conservative)**, headroom **188 (125)**. (The same tool gives 3489 for the v13 file as it stood
before this pass; CHANGES_v13 quoted 3461.) Letter: 6 pages, text ends on page 4, references on 5–6. Supplement: 34 pages. Both compiled
into `prl/build/` (pdflatex, bibtex, pdflatex ×2): 0 errors, 0 undefined references/citations, 0 overfull boxes, no float-too-large.

## 5. Numbers verified

* τ/4 = 2.89126·10⁻⁸/4 = 7.228·10⁻⁹; \|n_s − n_f\|/c_T = \|t₄\|/4 ≤ τ/4; \|p_T\|/c_T = 3\|α_T\|/(8π⁴) ≤ 1.113·10⁻¹⁰; (1+A)/(1−A) ≤ 1 + 1.45·10⁻⁸.
* Bounds 3.13206·10⁻⁶ / 2.89126·10⁻⁸, 72.0406 / 92.8135, 125.71 / 164.51 (C1_FINAL §0, REPRO; the c₁ = 0 values of the v13 formulas).
* Batch-5 count recounted from BATCH5_FOLD's table: J = 0, 2, 4–21 (20), 23, 25, 27, 29, 31 (5) = 25; J ≤ 32 minus J = 1, 3 = 31.
* QBIND table values transcribed and cross-checked against the UNIFIED_THEOREM box.
* CDPrakash: arXiv:1707.03007 (title/authors checked on arXiv), DOI 10.1007/JHEP11(2017)171 resolves (doi.org → Springer); entry
  already in `refs_verified.bib`, no bib change.

## 6. Left open / for the authors

* `% BATCH4` markers: Letter 1, Supplement 4 (Δ_g = 2000 direct rows).
* BATCH5_FOLD's strictest count drops "λ-level (8, 20)", but Table S5 has J = 8 also sign-free (L = 1.1–1.8); the Supplement lists only
  J = 20 as λ-level-only. Source inconsistency, not resolved here.
* Fig. 2 does not draw the partial [NE±] ranges of even J = 22–30 (see §1).
