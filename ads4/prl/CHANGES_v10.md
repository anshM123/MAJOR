# Changes v9 → v10 (24 Nov 2026): positivity at every gap Δ_g ≥ 1000

This version answers the referees' main objection, "positivity at one gap only". It draws on six new analyses, all unreviewed:
`ads4/GAP_SCALING.md`, `GAP_LARGESPIN.md`, `GAP_UNIFORM.md`, `GAP_C3.md`, `GAP_NEARGAP.md` and `GAP_MORSE.md`; their summaries are the
update boxes at the top of `ads4/UNIFIED_THEOREM.md`, newest last. Every number below was checked against the deliverables
themselves, not only the boxes.

No physics was recomputed. The only new arithmetic: √(16L_A) and τ₀ at 2000 / 5000 / 10⁴ / 10⁵, printed by `fig/make_fig.py`
from the published weights, which agree with DG_SCAN's Λ-forms. Nothing was committed, and no file outside `ads4/prl/` was edited.

## 1. The all-gap statement, as now written

The Letter states the following picture in plain words; the Supplement's new §S8 gives the grades.

* **Bounds.**
  * A theory with gap Δ′_g also satisfies the assumptions at every smaller gap. So, given positivity at 1000, the headline numbers
    bound every theory with Δ_g ≥ 1000. This is a new one-line remark, graded [D] (§S8.A).
  * Positivity at a theory's own gap gives that gap's sharper bound: τ₀ ∝ √(ln Λ)/Λ⁴, i.e. 2.89·10⁻⁸ / 1.93·10⁻⁹ / 5.33·10⁻¹¹ at
    1000 / 2000 / 5000.
  * In gap units the cubic constant *grows* slowly: |ĝ₃|Λ⁴ ≤ √(16L_A) ≈ 55√(ln(Λ/8.74)), which is 120 / 128 / 139 / 146 / 168 at
    1000 / 2000 / 5000 / 10⁴ / 10⁵.
  * **The brief's guess "bounds in gap units only get stronger" is wrong.** The Δ_g = 1000 numbers are the weakest of the family in
    absolute terms, not in gap units.
  * The quartic bound stays a Δ_g = 1000 statement, because Ω^Q_up was not examined at other gaps. None of the GAP analyses treats Ω^Q.
* **Certified at every Δ_g ≥ 1000:**
  * the flat floor;
  * Lemma C, whose bound is monotone in Λ [P];
  * Lemma LR∞ with no upper limit, its improved thresholds and J\*(m) [P + certified];
  * the ROUTEA window, uniform in the gap [D + certified], and Lemma C′, monotone in the gap [P];
  * hence the hybrid floor is rigorous on every b ≤ 2Λ for m ≥ 1.3Λ, and on b ≥ 114.06 for every m ≥ Λ.
* **Near-gap states.**
  * Error-controlled at 1000 and 2000; at 5000 for J = 5 and 12.
  * Nothing is negative.
  * Budget terms scale like Λ⁻² (×0.25 ± 0.01 at 2000; ×0.040–0.043 at 5000) for J ≤ 16.
  * The ratio of margin to budget grows 4–9× for J ≤ 8. J = 8 is now sign-free.
  * J = 12 and 16 close on the λ-level budget from K₂ alone.
  * J = 21, 27, 29 stay thin or [NE + M] at 2000, pending direct rows. Part of this is the reference-zero artefact.
* **Transfer between and beyond the computed gaps.**
  * Lemma T [P] and premise B_F [NE]. For J ≤ 7 the implication is proved: [P given B_F + windows], F_J ≥ 11.8.
  * For the other near-gap states: [D-sketch given H_mono-A].
  * For near-gap large spin: chart B with H_mono-B [M].
* **Certified eventual positivity** above Δ₁\*(J) = 0.9·10⁵ (J = 0) … 4.2·10⁵ (J = 8).
  * Grade: [P given TAIL, Π ≤ 3 [M], a bounded K₂ part [NE], and ξ-continuity].
  * The graviton amplitude is certified only at J = 0, 2 (2.2 / 2.8–3.0·10⁵). For J = 4–8 the certified data are those of the
    scalar integrand with the same phase.
  * The Letter therefore says "for the lowest spins, given further technical premises", not "for spins up to 8".
* **Large spin.** The theorems' margins grow with the gap: Theorem W's d_χ/H is 0.489 / 0.246 / 0.100, and its EDGE-R tolerance is
  0.406 / 3.45 / 57 S_E,t. This holds given premise constants carried to each gap in AdS units [M].

## 2. Letter (`letter.tex`)

| place | change |
|---|---|
| header | v10 comment, and a list of the `% SEAM` / `% BATCH4` marks |
| Introduction | New sentence: "We do this at Δ_g = 1000 and carry it, assuming that near-gap corrections do not grow with the gap, to every larger gap." |
| Setting | The reason spins 1 and 3 are absent is shortened to a citation (length). |
| Setting (end) | "Of the eight-derivative couplings we bound only g₄" is removed; its citation (CGGHJM) moves to the Discussion's "all-plus R⁴ coefficient" (length). |
| Why it matters | The last sentence (ABJM below the bound) is removed, because the ABJM paragraph and the abstract already say it (length). |
| Functionals | "Applying the sum rules term by term … we never do this" is removed. Mechanism already says the planar intercept never enters (length). |
| Result, last paragraph | "we claim nothing at other gaps" is replaced by four statements: monotonicity in the gap; the sharper per-gap bound, τ₀ ∝ √(ln Λ)/Λ⁴; the gap-unit constant growing to 128 (2000) and 139 (5000); and positivity of Ω^EE carrying over to every Δ_g ≥ 1000, while Ω^Q_up was examined only at 300 and 1000. |
| Is the functional positive? | The connector sentence "For large spin, positivity follows from the conditional theorems described next" is removed, and "Everywhere else" becomes "Elsewhere" (length). |
| **new paragraph "Larger gaps."** | Placed after "Large spin", ≈ 200 words. It covers: the certified statements at every gap; the error-controlled rechecks at 2000 (two spins at 5000); budgets ∝ Δ_g⁻² up to spin 16; margin/budget ×4–9 up to spin 8; J = 21, 27, 29 thin pending direct rows; the no-growth premise and where it was checked (1000, 2000, 5000; exact scalar to 10⁴); "for spins up to 7 we prove that it suffices, with a safety factor above ten"; the certified remainder bound above 10⁵–4×10⁵ for the lowest spins, "given further technical premises"; and the large-spin margins at 5000 (a tenth instead of a half, edge tolerance ×140). |
| ABJM | "Our bound is dominated by the infrared logarithm" is removed (the Discussion says it); the sentence is merged (length). |
| Discussion | "We do not claim … results at other gaps" becomes "… positivity at larger gaps without the premise on near-gap corrections, results below Δ_g = 1000". In the next steps, "positivity over a range of gaps" becomes "a proof of the uniform-decay assumption and of that premise". The de Sitter sentence is shortened. |
| Fig. 1 caption | "The solid cubic curve above Δ_g = 1000 holds given the premise on near-gap corrections stated in the text; faded lines: the implication alone." |
| Fig. 2 caption | "In these units every larger gap looks the same or better. Grades and sources: Supplemental Material, Secs. S7, S8." |
| acknowledgments | The AI-use disclosure is **unchanged, verbatim**. |
| abstract | **Unchanged** (580 characters). Accuracy does not require a change: the abstract states the Δ_g = 1000 result, which remains true. Writing "Δ_g ≥ 1000 … ≤ 120 in units of the gap" would be wrong, because the gap-unit constant grows. |

## 3. Fig. 1 (`fig/make_fig.py`)

* Panel (a): the "this work" cubic curve is drawn solid (full opacity, lw 1.3) for Δ_g ≥ 1000. Small open circles mark 2000 and 5000.
* The legend gains "solid: Δ_g ≥ 1000, see text", and "lines: implication only" becomes "faded: implication only".
* Panel (b), the quartic, is unchanged and stays faded beyond 1000.
* The script prints the constants at 2000, 5000, 10⁴ and 10⁵ as well (√(16L_A) = 128.17 / 138.54 / 145.90 / 168.06).
* Fig. 2 was not redrawn; only its caption changed. An inset was judged not cheap enough: the plane is drawn at Δ_g = 1000, and the
  all-gap statement is about the same picture in gap units.

## 4. Supplement (`supplement.tex`)

* **New §S8 "Positivity at every gap Δ_g ≥ 1000"** (`sec:allgap`), with six new tables, S12–S17. The old §S8–§S12 become §S9–§S13.
  * **Intro and summary.** States the net grade: "at best the grade it has at 1000, plus the no-growth premise below Δ₁\*; not (Pos)
    proved for all Δ_g ≥ 1000; nothing about Ω^Q_up".
  * **S8.A Monotonicity in the gap, and the constants.** The new remark [D]; Eq. (g3gap), √(16L_A) at five gaps; τ₀ at five gaps; the
    quartic as an implication only.
  * **S8.B Table S12 (uniform statements).**
    * Weights: the S(p)/Λ² symbol, the inner profile, δk₀.
    * Flat floor.
    * Lemma C monotone, with its J-reach 13 / 17 / 23 / 35 / 61.
    * LR∞ for every gap: both certificates; thresholds 0.2375 / 0.1675 / 0.1175; uniform floor f_U.
    * J\* = 223 / 447 / 1119 / 2238.
    * The uniform ROUTEA window [114.06 / 110.74 / 107.51, 2Λ]; Lemma C′ monotone; the consequence for the hybrid floor.
    * Lemma T; the fixed-phase structure [D]; cycle rigidity and non-effective eventual positivity [D-sketch].
  * **S8.C Table S13 (near-gap states at 2000 / 5000, one row per mass).** Values and leverages at tw₂₀₀₀ / tw₁₀₀₀ (/ tw₅₀₀₀); grade
    at 2000. The reading covers (i) no negatives; (ii)–(iv) per spin group; the reference-zero artefact (x₀ values, Δ_ref ≥ 5·10⁴);
    the batch-4 jobs; (v) chart-A scaling of every budget term and of Γ(x); (vi) the seam, with the open seam spins and the OO
    states, and the seam as the x → 1 window of Lemma T.
  * **S8.D Transfer.**
    * Premise B_F and the Proposition [P given B_F + windows], with its proof line.
    * Table S14: F_J; the J = 8 and 12 failures.
    * Evidence for B_F.
    * Other near-gap states [D-sketch given H_mono-A], with its tests and the incremental form (≥ +0.07 / +0.38).
    * Near-gap large spin in chart B under H_mono-B, with the scalar numbers. Why chart A fails at the twist edge: the +− transfer
      factor goes 408 → 1550 → 9400.
  * **S8.E Certified eventual positivity.**
    * Lemma M and the certified Morse chart, with premises (i)–(iv).
    * Table S15: Δ₁(J, ξ), scalar J = 0–8 and graviton EE/C3H J = 0, 2.
    * Why 1000 is out of reach (r_c ≈ 0.24 needed, 0.014–0.052 certified).
    * GAP_C3's first lemma is superseded, together with its float-coefficient defect.
  * **S8.F Large-spin theorems vs Δ_g.** Table S16: Λ²D; L′ Λ_amp 714 / 2851 / 17804; max d/(F/2) 0.190 / 0.035 / 0.0050; edge
    multiples 12.6 / 228 / 9050; W d_χ/H 0.489 / 0.246 / 0.100; EDGE-R tolerance 0.406 / 3.45 / 57. The text adds the ×8.5 / ×140
    margins under H_W, and notes that the 5000 C2/vector entries are not converged.
  * **S8.G Coverage.** Table S17 gives the status at every gap by region, followed by a "Not claimed" list.
* **Changed at Δ_g = 1000 as a consequence** (GAP_LARGESPIN U3/U5 hold at 1000 too):
  * §S7.B "Coverage of the floor": the hole 119 < b < 245.2 at 1.3 ≤ m/Λ < 1.5 is closed, and the floor for m < 1.3Λ is rigorous at
    b ≥ 114.06 (was b ≥ 245.2).
  * Theorem W's P0 is [D + certified] at b ≥ 114.06: every sampled state with J ≥ 80, and (64, 1.062Λ).
  * Table S11's intermediate-b row is updated to match.
  * These are upgrades. Flag them to the authors, because they change a 1000 statement.
* **Scope sentences.**
  * §S1 maps the Letter's new plain phrases to grades, adds §S8 to "What is here", and adds §S8's analyses to the provenance
    paragraph's unreviewed list.
  * §S5.G (Λ-forms): "values of the implication, not bounds" now reads "bounds where positivity holds: Ω^EE ≥ 1000 in the sense of
    §S8; Ω^Q_up at 1000", plus the monotonicity remark.
  * The Table S4 caption is changed to match.
  * The full statement points to §S8.
  * The Theorem L′ statement points to §S8.F for larger gaps.
  * §S9, the conditional remark on ≥ 300, is now marked as superseded for Δ_g ≥ 1000.
  * §S12 "What is not claimed" now reads "any Δ_g < 1000 other than 300; at Δ_g > 1000 anything beyond §S8".
  * "Corrections to earlier drafts" gains an "In v10" item.

## 5. Marks for incoming results

| mark | where | what will change |
|---|---|---|
| `% SEAM` | letter.tex header and the "Larger gaps" paragraph; supplement Table S13 (comment above the tabular), §S8.C (vi), after Table S17 | Results for seam spins J = 9–11, 13–15, 17–20 (also 22–26, 28, 30–32) and the parity-odd OO states J = 4–10 at Δ_g = 2000, from `ads4/GAP_SEAM.md` when it exists: new rows in Table S13, reading (vi), the Table S17 row "other near-gap states", and possibly "We repeated the error-controlled near-gap checks" in the Letter. |
| `% BATCH4` | letter.tex header and the "Larger gaps" paragraph (the "Spins 21, 27 and 29 remain thin" sentence); supplement Table S13 (comment), §S8.C (iv), §S8.D (other near-gap states), after Table S17 | External direct rows for J = 12, 16, 21, 27, 29 at Δ_g = 2000 (`ads4/ext/jobs_batch4.txt`, 247 main and 98 optional seam jobs for J = 20, 23). They replace the K₂-only entries. If they close J = 21, 27, 29 as at 1000, delete the Letter's "remain thin" sentence and upgrade the grades in Tables S13 and S17. |

## 6. Inconsistencies and cautions found

1. **Gap units (brief vs DG_SCAN).** The brief suggested that bounds in gap units get stronger with the gap. They do not:
   |ĝ₃|Λ⁴ ≤ √(16L_A) grows like √(ln Λ), and the quartic 2L_AdS like ln Λ. Only the absolute τ₀ falls, like √(ln Λ)/Λ⁴. Both
   files now say this.
2. **"Certified above ~10⁵" (brief, UT box) vs GAP_MORSE §0.** That statement holds only under premises (Π ≤ 3 [M], K₂ bounded [NE],
   TAIL and ξ-continuity [D-sketch]). It covers spins J ∈ {0, 2, 4, …, 8} only, and covers the *graviton* amplitude only at J = 0, 2.
   The Letter says "for the lowest spins, given further technical premises".
3. **Lemma T vs B_F wording.** The UT box for GAP_UNIFORM says Lemma T needs only that "Γ/m² does not grow", i.e. the relative
   correction. U7's premise B_F is stated in c₂ units with a Δ_g-independent F_J, which is stronger than Lemma T requires. The
   Supplement states B_F as written in GAP_UNIFORM and notes that Lemma T would allow more. The Letter's "relative to its flat value,
   does not grow" is the weaker, Lemma-T form; "a safety factor above ten" refers to F_J ≥ 11.8 (≥ 10.6 after the windows).
4. **J = 8 at 2000.** The UT box quotes "sign-free L = 5.2", which is the value at the chart-A image of the 1000 twist edge. At the
   2000 twist edge, L = 3.9. Table S13 gives both.
5. **"Budgets ∝ Λ⁻²"** holds for J ≤ 16. For J = 21 and 27 the factor is only ×0.43–0.70, because of the reference-zero artefact
   (GAP_NEARGAP N1). The Letter says "up to spin 16".
6. **GAP_SCALING §3.1 says Theorem L′'s tolerances were "not recomputed".** GAP_LARGESPIN §2 recomputed them, so the Supplement uses
   GAP_LARGESPIN. Its 1000 value of max d/(F/2) is 0.190 with the uniform floor, against 0.133 with LRINF's floor (Table S10); both
   are given.
7. **The 1000 statements upgraded by GAP_LARGESPIN U3** (the hole 119 < b < 245.2, and W's P0) were [NE] in v9. They are now
   [D + certified], on an unreviewed analysis.
8. **`ads4/UNIFIED_THEOREM.md` was not edited.** The brief's "(update boxes at the top of …)" was read as describing where the new
   results are summarised; the rules allow edits only under `ads4/prl/`.
9. **The monotonicity remark (§S8.A) is new wording, not taken from a deliverable.** It is a one-line logical consequence of A1′/A2.
   The authors should confirm it.

## 7. Build and length

* Build: `(cd fig && python3 make_fig.py)`, then pdflatex / bibtex / pdflatex ×2 per file into the git-ignored `prl/build/`.
  * **Letter:** 6 pages, 0 errors, 0 undefined references or citations, 0 overfull boxes, 7 underfull notices. The text now ends near
    the top of page 5; references fill pages 5–6.
  * **Supplement:** 30 pages (was 25), 0 errors, 0 undefined, 0 overfull.
  * BibTeX: the usual `jnrlst` notice and "empty year in SM".
* Length (`scratchpad/tools/wc_prl7.py letter.tex 3 0.9375 1.0385`):
  * prose 3056 words and 236 inline formulas (298 conservative);
  * displays 48; figures 344;
  * **TOTAL ≈ 3684 (3746 conservative)**, against a limit of 3750.
* Abstract: 580 rendered characters (unchanged).
