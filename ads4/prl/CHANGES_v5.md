# CHANGES v4 → v5 (HIGHSPIN places filled; large-spin theorems; front matter)

*23 Nov 2026. Writing pass only: no physics computation was redone. Sources (all unreviewed except INDEP_CHECK, an independent
check): `ads4/UNIFIED_THEOREM.md` update boxes and the deliverables they cite, each read directly — `HIGHSPIN.md`, `SMALLNU_WINDOW.md`,
`N2_SUBREGGE.md`, `DIRCV.md` (§8 latest), `LARGESPIN_ANALYTIC.md`, `LARGESPIN_W.md`, `LARGESPIN_L.md`, `P2_SPINNING.md` (its §7 was
still being written), `PW_IMAGES.md`, `INDEP_CHECK.md`. Files edited: `letter.tex`, `supplement.tex`, `macros.tex`,
`refs_verified.bib`, `NOTES.md`; new: this file. `fig/` untouched (no change needed). Nothing outside `ads4/prl/` was touched; nothing
committed. Line numbers refer to the files as of this writing.*

## 1. HIGHSPIN places (all filled; `\HSdraftfalse`, macros.tex:76)

| v4 place | v5 content | where | source |
|---|---|---|---|
| Letter (i) near-gap J ≥ 12 "[NE + M]" | [NE±] near gap for J = 0, 2, 4, 6, 8, 12, 16, 21, 29 (+ binding J = 5, ≥ +0.198, L = 78); leverage 1.7–22 for J ≤ 16 (J = 8 λ-level), 1.2 for J = 21, 29 at the twist edge, ≥ 1.7 above; J = 7, 9–11, 13–15, 17–20, 22–28, 31, 32, 36, 40, 48 [NE + M]; even J = 26, 30, 34–46 [M] | letter.tex:252–296 | NEARGAP2 §0; HIGHSPIN §4, §6; SMALLNU_WINDOW §0, §3 |
| Letter (ii) J = 29 "−3.6 vs +1.39 [M]" | +0.972 ± 0.018; −3.62 = K₂ node at the heavy reference's partial-wave zero (x₀ = 0.8008); +1.39 = spline aliasing; J = 31, 32 → +1.02, +1.00; J = 29 [NE±] only with the signed measured window (sign-free envelope L = 0.94) | letter.tex:262–272 | HIGHSPIN §2; SMALLNU_WINDOW §3.1 |
| Letter (iii) H_LJ(1000) | [D-sketch + M]; sign-free violation at (1100, 100) is a normalisation of the leading term [DV], harmless by congruence [P]; leverage ≥ 1.6 central in every fit, ≥ 1.0 at 2σ for every A + C fit, **0.5–0.6 at 2σ in the λ-free A + B + C 9-node fit** (stated) | letter.tex:283–288 | DIRCV §8; N2_SUBREGGE; LARGESPIN_W W1c |
| (new) ν < 50 window | measured; convention kept below the lowest node | letter.tex:268–271; supplement.tex:885–895 | SMALLNU_WINDOW |
| (new) H_R′, leading-Regge floor | H_R′ [D-sketch + M] for b ≥ b_A(m) ≈ 8.5–15, [NE + M] below (weak point (iii)); floor proved at every b for m ≥ 1.5Λ, at b ≥ 1.04Λ for every m [P + DV; certified] | letter.tex:273–276, 289–291 | LARGESPIN_ANALYTIC §2, §5.2; LARGESPIN_L §3 |
| Supplement Table S5 (J ≥ 12 row) | three rows: J = 12, 16 [NE±]; J = 21, 29 [NE±] (thin); J = 23, 27, 31 [NE + M]; "other" row adds J = 20, 32 and the direct-edge list | supplement.tex:776–779 | HIGHSPIN §2–4; SMALLNU §3 |
| Table S6 (J = 12, 16, 29) | published minima kept, direct values in parentheses; J = 29, 31, 32 updated | supplement.tex:800–802 | HIGHSPIN §2, §4 |
| Table S7 (H_M, H_edge, H_LJ(1000)) | H_M with the K₂ calibration (4–5 % of the last step); H_edge computed for J = 18–25, 27–29, 32; new row H_low; H_R′ upgraded (ROUTEA window m ≥ 1.5Λ only); H_LJ(1000) with Table S8 | supplement.tex:815–820 | HIGHSPIN §3–4; SMALLNU §0; LARGESPIN_ANALYTIC §5.2, N4; DIRCV §8 |
| Supplement §S1 provenance "being recomputed" | replaced (recomputed; pointers) | supplement.tex:69–75 | — |

New supplement paragraphs in §S6: the J = 29 edge and the reference-zero artefact (supplement.tex:869–883), the window ν < 50
(885–895), H_LJ(1000) at (1100, 100) + **Table S8** (897–927). "How it is computed" mentions direct rows (745–747).

## 2. Large spin (task 2)

* **Letter** "Large spin." paragraph (letter.tex:300–321): conditional reduction, not a proof; the exact magic-block identity
  ([P] to r⁸, [DV] to r⁴⁰, [D] all orders) and the monopole term; Theorem W (J\* = 50, 100 at margin 2; m ≤ 1.3Λ; implication [P],
  costs [DV] on sampled spins ≤ 344); Theorem L (J ≳ 0.52m, 857 at the gap; leading-Regge action proved positive beyond the AdS
  radius; Λ_amp tolerated to 43 at J ≤ 10⁶); open premises named: spinning O(1/m²) magnitudes (≈ 2–3 [NE] at J ≤ 100, [M] above) and
  beyond-all-orders [M]. `% P2FOLLOW` comment at letter.tex:317–319.
* **Supplement §S7 "Large spin: two conditional theorems"** (supplement.tex:929–1125; old §S7–§S10 are now §S8–§S11): A. the identity
  (Eq. S-magic), consequence, monopole cost D = −9.98332·10⁻⁴, Lemma I, congruence lemma, Weyl reflection, CPPR weight-shifting;
  B. Lemma LR∞ with r(β) table and the floor coverage (the LARGESPIN_ANALYTIC "Lemma C′" is renamed "extended Mehler–Dirichlet curvature
  bound" to avoid a clash with §S5.B's Lemma C′); C. Theorem W (premises P0, P1, P2, P2′, P4 with grades, proof, **Table S9** J\* and
  class tolerances, evidence on the premises, the open (P4) subtraction point); D. Theorem L (statement, **Table S10** J\*(m) and
  tolerances, reading); E. coverage **Table S11** and what is not claimed. `% P2FOLLOW` at the (S) check, Table S9, the model constants
  and the (P4) paragraph.
* New reference **CPPR** (Costa, Penedones, Poland, Rychkov, *Spinning Conformal Blocks*, JHEP 11 (2011) 154, arXiv:1109.6321),
  verified on the arXiv API and INSPIRE (recid 930033): refs_verified.bib:763–780; cited at supplement.tex:971.

## 3. "Why it matters" (task 3) and outlook (task 4)

* New paragraph (letter.tex:87–99): if positivity and the inputs hold, up to c₁, τ ≤ 2.89·10⁻⁸ at Δ_g = 1000, i.e. < 10⁻⁸ of the
  conformal-collider range τ ≤ 4 (macro `\CollFrac`, macros.tex:78–81, follows the H_LJ(300) switch); a number, not a scaling law;
  quartic order reached; ABJM's 2ζ(3) below the R⁴ upper end as a consistency check against string theory, "not a test of positivity".
* Intro (letter.tex:72–79): "positivity is proved in the flat limit; in AdS it is supported numerically, with explicit error budgets on
  the low-lying heavy states, and at large spin it is reduced to a conditional theorem".
* Discussion (letter.tex:350–358): not-claimed list adds "unconditional positivity at large spin"; next steps rewritten (spinning
  O(1/m²) coefficients first); **de Sitter sentence only here** (final paragraph), not in the abstract.
* Abstract (letter.tex:33–41): "(proved in the flat limit, tested numerically in AdS)" → "(proved in the flat limit; in AdS no
  computed state is negative)". Nothing else; no large-spin or de Sitter claim. 563 characters rendered.
* Result paragraph (letter.tex:243–245): independent reproduction of the EE-kernel reduction and the counting lemma (INDEP_CHECK).
  Supplement reproduction section (supplement.tex:1159–1203): four INDEP_CHECK rows and the "not covered" list.

## 4. Front matter (coordinator message, authors' data)

Both files: `\affiliation{Independent researcher}` for each author, `\email{ansh.mishra2025@gmail.com}`,
`\email{aryan.senthilk@gmail.com}` (letter.tex:22–28, supplement.tex:27–33); the TODO(authors) comments are gone. **No
acknowledgments**: the `acknowledgments` environment and its placeholder were removed from letter.tex. The word counter now counts to
`\bibliography`.

## 5. Other

* Figure environment moved from before "Result" to before "Functionals" (letter.tex:163) to avoid a stuck float; figure unchanged.
* Supplement record corrections (supplement.tex:1258–1266) and not-claimed list (1245–1248) updated.
* NOTES.md: header, file table (stale v3 duplicate removed), build status, length (v5), switches, provenance (v5 table), referee
  status (front matter done), open items, **§7 v5 source inconsistencies (nine items)**, §8 fallback status.

## 6. Length and build

* APS words (tool `scratchpad/tools/wc_prl5.py letter.tex 3 0.9375`; v4's tool dropped the `\HLJstatus` macro text, ≈ 45 words):
  **≈ 3110** (≈ 3210 conservative); limit 3750. Abstract 563 characters. Text ends ≈ 0.3 down the right column of page 4
  (≈ 3.65 pages). Letter 5 pages, Supplement 20 pages.
* Both compile (pdflatex; `(cd build && BIBINPUTS=..:. bibtex $f)`; pdflatex × 2) with no errors, no undefined references or
  citations, no overfull boxes, no stuck floats. BibTeX notices: `jnrlst` style notice; "empty year in SM" (intended).

## 7. What remains

1. P2 follow-up (running): K₄ to O(1/m⁴), the (S) check against 1809.05111, the O(1/m⁴) hybrid subtraction of premise (P4); update
   the `% P2FOLLOW` places.
2. Not addressable by writing: the spinning O(1/m²) coefficients (P2 magnitudes); error budgets for J = 7, 9–11, 13–15, 17–20,
   22–28, 31, 32, 36–48 and near-gap 345 ≤ J < 857; J ≥ 150 window data; c₁ / F₂; a re-optimised Q.
3. Reviews: every analysis used here except REPRO, REPRO_ENGINE, INDEP_CHECK is unreviewed; independent checks of PW_IMAGES'
   identity and Lemma LR∞ (interval constants) would be cheap and valuable.
4. INSPIRE sweep before submission; `refs.bib` (unused legacy) still stale.
