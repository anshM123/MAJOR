*Revised 24 Nov 2026, **v10** (positivity at every gap Δ_g ≥ 1000, after GAP_SCALING, GAP_LARGESPIN, GAP_UNIFORM, GAP_C3,
GAP_NEARGAP, GAP_MORSE; all unreviewed; no physics recomputed): the Letter gains a paragraph "Larger gaps" (certified layers at every
gap; near-gap rechecks at 2000, two spins at 5000, budgets ∝ Δ_g⁻² to J = 16; the no-growth premise, proved sufficient for J ≤ 7 with
F_J ≥ 11.8, checked at 1000/2000/5000 and on the exact scalar to 10⁴; certified remainder above 10⁵–4·10⁵ for the lowest spins given
premises; large-spin margins grow), the Result's last paragraph (monotonicity in the gap [D]; per-gap constant |ĝ₃|Λ⁴ ≤ 128 / 139 at
2000 / 5000, i.e. the gap-unit constant GROWS like √ln Λ while τ₀ falls like √ln Λ/Λ⁴; Ω^Q_up only at 300/1000), and scope edits
(Intro, Discussion, Fig. 1/2 captions; Fig. 1 cubic curve solid above 1000). Abstract and the AI-use disclosure unchanged. Supplement:
new **§S8 "Positivity at every gap Δ_g ≥ 1000"** (Tables S12–S17); the 1000 hybrid-floor hole 119 < b < 245.2 and W's P0 upgraded via
the uniform ROUTEA window. Marks **`% SEAM`** (GAP_SEAM.md: seam spins and OO states at 2000) and **`% BATCH4`** (external direct rows
J = 12–29 at 2000). Change log: **`CHANGES_v10.md`** (incl. inconsistencies §6).*

*Revised 24 Nov 2026, **v9** (cheap fixes after `REFEREE_REPORT_4.md`; no new physics): title "How Close to Einstein Must
Holographic Gravity Be?"; abstract leads with the assumption and |ĝ₃|Λ⁴ ≤ 120 (τ ≤ 2.89·10⁻⁸ as "i.e."), "proved in the flat-space
limit for the cubic bound and supported numerically in AdS", "numbers are for c₁ = 0"; Introduction lists the large-spin premises
(decay assumption + uncomputed graviton edge integral) instead of "one technical assumption"; Result stated without light single
traces in T×T (implication (1) [D]), light single traces as a one-sentence extension resting on (E′), (2b) moved there; M6 wording
(J = 29 twist edge only with the measured sign; "to about one part in 10⁶"; "we derive"; "separate code … EE-kernel reduction and
the counting lemma"); Theorem L′ "remainder checked on sampled states up to 5Λ"; identity "prove through eighth order, check to
fortieth"; Fig. 2: Theorem W drawn as its 47 sampled states (dots) over a hatched "between samples" band, J = 29 twist edge pale.
Supplement: new title; version narrative removed from the main flow (S1 "What is here", provenance, S2–S7 inline "v7/v8"
remarks) and collected in a final section "Corrections to earlier drafts" (`sec:corr`); S5.F states which representation each
light-state statement refers to; S2.B mentions the ×100 bump of the continuum k = 2 weight at Λ+30–33; S2.C says the grids are
tabulations, not functionals. Change log: **`CHANGES_v9.md`**.*

*Revised 24 Nov 2026, **v8** (writing pass after AG_PROOF, DSKETCH_UP, EDGE_R, NU_TAIL_CHECK and the J = 31 addendum of EXT_FOLD;
title/abstract of commit 65c1a44 kept): (i) **large spin**: (AG′) is replaced by the strictly weaker pair (AG_U) (uniform decay rate,
[A] beyond all orders) + (EDGE-R) (one real-axis edge integral of the graviton, [M], not closed) via Lemma F and Lemma P4W-χ [P];
Theorem W: sampled J\* = 26 given (EDGE-R) (v7's "42" superseded; Table S9 rewritten with the EDGE_R corrections: exact scalar at
(26, 1.026Λ) 0.285 at the tabulated cutoff, margin 1.42, grid maximum 0.405 vs tolerance 0.406; "plain margin ≥ 520" not a graviton
statement); Theorem L′ on the uniform rate alone (ε_R ≤ 0.133, EDGE tolerance 18.9×); Fig. 2's dashed J = 42 line removed; every
`% AGPROOF` mark resolved. (ii) **c₁** universality and ν-independence [D] (L3, §S5.E). (iii) **light single traces**: the
leading-Regge "AdS image" is refuted; light value = Regge polynomial given (E′) [D-sketch]; LS6 ⇔ μ² > −692.5 / −7748.7; ε\* expected
to fail. (iv) **functional definition**: k = 2 pieces defined by the continuum weights; k = 2 grids read as zero beyond Λ + 10
(`data/truncate_grids.py` writes `*_trunc.txt` copies; published files unchanged; data/README §4, §7; MANIFEST). Change log:
**`CHANGES_v8.md`** (incl. inconsistencies §3).*

# PRL draft: notes, build, length, provenance, referee responses, fallback plan, open items

*Revised 23 Nov 2026, **v7** (writing pass after INDEP_CHECK2, LRINF_TIGHT, P2_CONTRACT, P4_WINDOW, EXT_FOLD incl. its J = 31 addendum,
an outside critique and the authors' instruction that the Letter read like a PRL): the Letter's main text has **no grade tags** (one
sentence points to the legend, Supplement §S1, which now maps every plain phrase to its grade); internal labels (A1′/A2, (Pos), H_LJ,
W/L′, AG′, [NE±], K₂) are gone from the main text and Fig. 2; the positivity paragraph ("Is the functional positive?") is six sentences
around Fig. 2; the large-spin paragraph describes the two conditional theorems and their one technical assumption in words. Supplement
§S7 rewritten: Lemma LR∞ certified to β ≥ 0.37, Theorem L′ with J\*(m) = ⌈(0.37211 m − 1)/2⌉ (223 at the gap), computed O(1/m²)
corrections (monopole cancels), window Lemma P4W and the assumption (AG′), Theorem W in window form (47 sampled states, J\* = 26, or 42
with the Ω^U +− edge factor); the v6 "345 ≤ J < 857 uncovered" statement is gone (W and L′ overlap at the gap). Fig. 2 redrawn
accordingly. Change log: **`CHANGES_v7.md`** (incl. the PRLs read for voice, and inconsistencies §4). Marks: **`% AGPROOF`** at every
place that changes if (AG′) is proved (list in CHANGES_v7 §3); v6's `% P2CONTRACT` marks are resolved.*

*Revised 23 Nov 2026, **v6** (writing pass after `REFEREE_REPORT_3.md` and `ads4/REGGE_TREE.md`): order of limits stated (A2,
Mechanism, Functionals; Supplement Table S2 rows A2/R1/H_K and two paragraphs at the end of Sec. S5.D); the Letter's "Status of
positivity" + "Large spin" (≈ 830 words) condensed to one paragraph (≈ 310 prose words) with a new **Fig. 2** (coverage of (Pos) in
the (J, m/Λ) plane, `fig/make_fig2.py`); Theorem W/L scope and premises stated (P0 [NE], P1 [D-sketch] after the (S) failure, flat ν⁶
tower subtraction [I], Λ_amp from one state); "why it matters" reframed (Λ⁻⁴ scaling, leading order in 1/c_T); B_eff(Λ) fixed; Fig. 1
markers/key (m14). Change log: **`CHANGES_v6.md`**. Marks: `% P2CONTRACT` (constants the running three-channel contraction may
update); the v5 `% P2FOLLOW` marks are resolved. No physics recomputed except B_eff/2L_AdS (from `fig/make_fig.py`'s integrals).*

*Revised 23 Nov 2026, v5 (writing pass after HIGHSPIN, SMALLNU_WINDOW, N2_SUBREGGE, DIRCV §8, LARGESPIN_ANALYTIC/_W/_L,
P2_SPINNING, PW_IMAGES, INDEP_CHECK): the HIGHSPIN places are filled and `\HSdraftfalse`; the Letter gains a "Why it matters"
paragraph, a "Large spin" paragraph (conditional Theorems W and L) and one de Sitter outlook sentence in the Discussion; the
Supplement gains §S7 "Large spin: two conditional theorems" (old §S7–§S10 become §S8–§S11), Table S8 (H_LJ(1000) leverage) and the independent-check rows; front
matter filled (authors' message), acknowledgments removed. Change log: **`CHANGES_v5.md`**. No new physics computation in v5.*

*Previous header (v4): Revised 23 Nov 2026 (v4 writing pass after `REFEREE_REPORT_2.md`: the Letter is narrowed to item (1) at Δ_g = 1000 plus item (2a)'s
upper end, rewritten for a general reader, and given one figure; everything removed from the Letter is in the Supplement). The only new
arithmetic in v4: the figure's curves (same integrals as `code/repro/bounds.py`, recomputed from `prl/data/`), λ(Λ = ML) = 2.51·10⁷ /
3.15·10⁹, the |ĝ₃| ratio to the transplanted CHLPS optimum (10.1 / 10.3), and the Q3\* data file (exact division, fingerprint
reproduced). Unrefereed. Nothing committed. Source of truth: `ads4/UNIFIED_THEOREM_v3.md` (UT3). Change logs: `CHANGES_v3.md`
(v2 → v3), **`CHANGES_v4.md` (v3 → v4, referee 2 point by point)**.*

Files in `ads4/prl/`:

| file | content |
|---|---|
| `letter.tex` | the Letter (REVTeX 4.2, `aps,prl,twocolumn,floatfix`), v4 |
| `supplement.tex` | Supplemental Material (one column): grades and composites, definitions/data, inputs table, **full statement + I1–I6 + lower-end outlook (new §S4)**, derivation sketches, Λ-forms (now with the displayed forms), positivity evidence, conditional remark on Δ_g ≥ 300, reproduction, ABJM |
| `macros.tex` | shared macros: H_LJ(300) switch, `\QuartRep`, `\GthreeHead` (new), HIGHSPIN marker `\HS` |
| `fig/make_fig.py`, `fig/fig1_bounds.pdf` | Fig. 1 (new): bounds vs Δ_g from the published data; `fig1_bounds_preview.png` is a raster preview |
| `fig/make_fig2.py`, `fig/fig2_coverage.pdf` | Fig. 2 (v6): coverage of (Pos) for Ω^EE at Δ_g = 1000 by grade in the (J, m/Λ) plane; boundaries transcribed from the deliverables (docstring); `fig2_coverage_preview.png` |
| `refs_verified.bib` | bibliography; v4: JHEP year convention, key BCPW → BCP, CMT, HaringZhiboedov, HKT, `SM` |
| `refs.bib` | legacy, unused (not updated in v4) |
| `data/` | published functional data; v4 adds `q3star_factored.json` + `make_q3star.py`, README errata §5, MANIFEST lines |
| `build/` | all build outputs (git-ignored) |

| `CHANGES_v3.md` … `CHANGES_v10.md` | change logs v2 → v3, v3 → v4 (referee 2), v4 → v5 (HIGHSPIN filled, large spin), v5 → v6 (referee 3), v6 → v7 (plain prose, large-spin theorems), v7 → v8 (AG_U + EDGE-R, c₁ [D], light single traces, k = 2 definition), v8 → v9 (referee 4 cheap fixes), v9 → v10 (positivity at every gap Δ_g ≥ 1000, Supplement §S8) |
| `data/truncate_grids.py`, `data/*_trunc.txt`, `data/k2/*_trunc.txt`, `data/TRUNC_SHA256SUMS` | v8: the k = 2 grids cut at ν ≤ Λ + 10 (the reading the Supplement now prescribes); published grids unchanged |

(v5: the stale v3 copy of this table that followed here was removed; `letter.tex` is v5, `supplement.tex` has the new §S7,
`macros.tex` has `\HSdraftfalse` and `\CollFrac`, `refs_verified.bib` adds CPPR 1109.6321.)

---

## 1. Build

```
cd ads4/prl
ulimit -v 3000000; export OMP_NUM_THREADS=1
(cd fig && python3 make_fig.py)          # Fig. 1 (prints the headline constants as a check)
(cd fig && python3 make_fig2.py)         # Fig. 2 (v6; prints the twist-edge / Theorem-L crossings as a check)
(cd data && python3 make_q3star.py)      # optional: regenerates q3star_factored.json, asserts the certified fingerprint
for f in letter supplement; do
  pdflatex -interaction=nonstopmode -output-directory=build $f.tex
  (cd build && BIBINPUTS=..:. bibtex $f)
  pdflatex -interaction=nonstopmode -output-directory=build $f.tex
  pdflatex -interaction=nonstopmode -output-directory=build $f.tex
done
```

Status (24 Nov, **v10**): both compile with no errors, no undefined references or citations, no overfull boxes (Letter: 7
underfull-box notices). `build/letter.pdf` 6 pages: the text ends near the top of page 5 (the new "Larger gaps" paragraph, ≈ 0.25
page), references fill pages 5–6. `build/supplement.pdf` 30 pages (new §S8 with Tables S12–S17). `fig/make_fig.py` (v10) draws the
cubic curve solid above 1000 and prints √(16L_A) at 2000 / 5000 / 10⁴ / 10⁵ (128.17 / 138.54 / 145.90 / 168.06). Fig. 2 not redrawn.

Status (24 Nov, **v8**): both compile with no errors, no undefined references or citations, no overfull boxes, no stuck or oversized
floats (Letter: 7 underfull-box notices; Supplement: hyperref "token not allowed" for math in section titles). `build/letter.pdf` 5
pages, text ends ≈ 0.45 down the right column of page 4; `build/supplement.pdf` 25 pages. Optional data step (v8):
`(cd data && python3 truncate_grids.py)` rewrites the four `*_trunc.txt` copies and `TRUNC_SHA256SUMS` (deterministic; < 1 s).

Status (23 Nov, **v7**): both compile with no errors, no undefined references or citations, no overfull boxes, no stuck floats
(Letter: 7 underfull-box notices; Supplement: hyperref "token not allowed" for math in section titles). `build/letter.pdf` 5 pages,
text ends ≈ 0.55 down the right column of page 4; `build/supplement.pdf` 23 pages. `make_fig2.py` prints Theorem L′'s twist-edge J\*
(223) as a check.

Status (23 Nov, **v6**): both compile with no errors, no undefined references or citations, no overfull boxes, no
stuck floats (Letter: 5 underfull-box notices; Supplement: Table S2 set in `\scriptsize` after the A2/R1/H_K rows grew).
`build/letter.pdf` 5 pages, text ends ≈ 0.55 down the right column of page 4 (≈ 3.8 pages); `build/supplement.pdf` 21 pages.

Status (v5, for the record): both compile with **no errors, no undefined references or citations, no overfull boxes, no stuck floats**
(the Letter has 5 underfull-box notices, cosmetic). `build/letter.pdf` **5 pages** (text ends ≈ 0.3 down the right column of
page 4, i.e. ≈ 3.65 pages incl. the title block; references fill the rest of page 4 and page 5); `build/supplement.pdf` **20 pages**.
In v5 the figure environment was moved to just before "Functionals" (it had become a stuck float after the new paragraphs).
(v4 status: Letter 4 pages, text ending ≈ 0.95 down page 3; Supplement 15 pages.) BibTeX: one
`jnrlst` style notice and "empty year in SM" (the Supplemental-Material pointer, intended). The JHEP "missing number" warnings are gone.
Figure: `cd fig && python3 make_fig.py` (needs mpmath, matplotlib) before compiling the Letter.

## 2. Length

**v10** (same tool and arguments; it counts the acknowledgments too): abstract **580** rendered characters (unchanged); prose 3056
words, 236 inline formulas (conservative 298), 3 display rows (48), Fig. 1 180, Fig. 2 164; **TOTAL ≈ 3684 (≈ 3746 conservative)**;
limit 3750. To make room for "Larger gaps" (≈ 200 words) v10 cut ≈ 90 words elsewhere (CHANGES_v10 §2: spins 1, 3 reason, g₄-only
sentence, Why-it-matters ABJM sentence, the term-by-term sentence in Functionals, a connector in the positivity paragraph, ABJM
IR-log clause, de Sitter sentence shortened). Little room is left (≈ 65 words APS, ≈ 5 conservative).

**v9** (`scratchpad/tools/wc_prl7.py letter.tex 3 0.9375 1.0385`): abstract **580** rendered characters (pdftotext, whitespace
normalised); prose 2850 words, 215 inline formulas (conservative 268), 3 display rows (48), Fig. 1 180, Fig. 2 164;
**TOTAL ≈ 3457 (≈ 3510 conservative)**; limit 3750. Letter 5 pages, Supplement 25 pages; 0 errors, 0 undefined references,
0 overfull boxes.

**v8** (same tool, `scratchpad/tools/wc_prl7.py letter.tex 3 0.9375 1.0385`): abstract **588** rendered characters (the new abstract
of 65c1a44, unchanged); prose 2733 words, 209 inline formulas (conservative 261), 3 display rows (48), Fig. 1 180, Fig. 2 164;
**TOTAL ≈ 3334 (≈ 3386 conservative)**; limit 3750, aim ≤ 3350 met.

**v7** (tool `scratchpad/tools/wc_prl7.py letter.tex 3 0.9375 1.0385`: wc_prl6 plus expansion of `\HLJstatusLetter`): abstract
**563** characters (unchanged); prose 2715 words, 212 inline formulas (conservative 264), 3 display rows (48), Fig. 1 180, Fig. 2 164;
**TOTAL ≈ 3319 (≈ 3371 conservative)**; limit 3750. Plain prose costs words (v6: 3207), offset by moving numbers to the Supplement.

**v6** (tool `scratchpad/tools/wc_prl6.py letter.tex 3 0.9375 1.0385`: as wc_prl5, for any number of figures):
abstract **563** characters rendered; prose 2501 words, 314 inline formulas (conservative 423), 3 display rows (48), Fig. 1 180,
Fig. 2 164 (aspect 3.375/3.25 = 1.0385); **TOTAL ≈ 3207 (≈ 3316 conservative)**; limit 3750. Text ends ≈ 3.8 pages.


**v5** (tool `scratchpad/tools/wc_prl5.py letter.tex 3 0.9375`: as wc_prl4, plus it expands the prose macro `\HLJstatus`, which
wc_prl4 silently dropped (≈ 45 words), and counts to `\bibliography`, since there are no acknowledgments now):

| item | value |
|---|---|
| abstract | **563 characters** rendered incl. spaces (555 without the spaces pdftotext inserts in ⟨TTTT⟩); limit 600 |
| main text, APS counting | **≈ 3110** (inline formula = 1 word) / **≈ 3210** (conservative); prose 2567 words, 315 inline formulas, 3 display rows (48), Fig. 1 180, no table; limit 3750 |
| rendered | text ends ≈ 0.3 down the right column of page 4 (≈ 3.65 pages; target ≤ 3.95) |

### v4 (for the record)

| item | value | method |
|---|---|---|
| abstract | **551 characters** rendered incl. spaces (545 after removing the spaces pdftotext inserts inside ⟨TTTT⟩) | pdftotext → whitespace-normalised `len()`; limit 600 |
| main text, APS counting | **≈ 2395** (inline formula = 1 word) / **≈ 2480** (conservative: long inline formulas counted as ⌈chars/6⌉ words) | prose 1928 words from `\maketitle` to the acknowledgments, captions included; 239 inline formulas; 3 display rows × 16 = 48; Fig. 1 single-column, aspect 3.375/3.6 = 0.94 → 150/0.94 + 20 = 180; no table |
| tool | `scratchpad/tools/wc_prl4.py letter.tex 3 0.9375` | display rows checked in the PDF (Eqs. (1)–(3) one row each) |
| rendered | text ends ≈ 0.95 down page 3 of the REVTeX `prl,twocolumn` layout, i.e. ≈ 2.95 pages incl. the title block | page images |

Target was ≤ 3300 APS words and ≤ 3.7 pages: met with ≈ 900 words of headroom for the HIGHSPIN text (v3: ≈ 3080 by the old tool, which
undercounted, ≈ 4.2 pages).

## 3. Switches

* **H_LJ(300)**: `\HLJclosedfalse` (open; headline Δ_g = 1000). Unchanged. The abstract, the introduction's headline and the
  discussion use `\CubicHead`, `\DgHead`, `\QuartHead`, `\GthreeHead`, so they follow the switch; the Result paragraph prints both gaps.
  Fig. 1 marks 1000 filled, 300 open, and does not follow the switch (re-run `make_fig.py` with the markers swapped if it is flipped).
* **HIGHSPIN**: **`\HSdraftfalse` (v5: every place filled; list in CHANGES_v5.md §1)**. The v4 text follows. Places in v4 (each with a `% HIGHSPIN:` comment): `letter.tex` Status-of-positivity items (i) near-gap
  J ≥ 12, (ii) J = 29 edge, (iii) H_LJ(1000); `supplement.tex` Tables S5 (J ≥ 12 row), S6 (J = 12, 16, 29), S7 (H_M, H_edge, H_LJ(1000)).
  The v3 Letter's Table I placeholders are gone with the table (its content is in Tables S5/S7). When HIGHSPIN lands: update those
  places, UT3 §5.2/§5.6/§10, and set `\HSdraftfalse`. If any state comes back negative, use §8 below.

## 4. Provenance of the numbers (letter items new or changed in v3)

| letter item | source |
|---|---|
| A1′ (higher-spin gap), spins 1, 3 absent from T × T | A1PRIME §2.1, §6 |
| light Regge spins 3 / 4; drop-out at k = 4, 5; block action at k = 2, 3 | A1PRIME §2–§3, §5 |
| Ω^EE positive on light states at every mass; Ω^Q negative on scalars below 0.46Λ | A1PRIME §4.2–§4.3 (exact polynomials) |
| repair 0.04 Ω^EE; 125.71 / 164.51 | NEARGAP2 §0 ("Repaired Q"); UT3 §3 |
| c₁² coefficients 1.29·10⁻⁴ / 1.03·10⁻⁶ (repaired) | this pass: `57600(I_Q+εI_D)²/((κ+ε/8)Λ⁴)` (A1PRIME §4.3 formula; reproduces A1PRIME's ε* values 1.2598·10⁻⁴ / 1.0062·10⁻⁶) |
| 1 % tolerances 4.6 / 60 (item 1), 99 / 1262 (item 2, either form) | C1_VALUE §6 (published form: 103 / 1306); this pass (repaired form) |
| c₁ universal; shock part identically 0; energy-channel T exchange; zero given F₂ | C1_VALUE §0, §2–§4 |
| ε_W ≤ 2.7·10⁻¹³ / 9.1·10⁻¹¹ [certified]; 2.4·10⁻⁷; k = 2 ≤ 1.1·10⁻⁹ relative | HK_RESIDUAL §3.3; DESHIFT §2.3; CERT_CLOSE §3.2 |
| R⁴ light function Δφ⁻² log Δφ, criterion Δφ^{−7/2} | R4ARC_HEFT §0, §2 |
| H_EFT(a) pointwise for dispersive terms [D; DV] | R4ARC_HEFT §3 |
| g₅ not cancelled; s_d mass-independent Λ → 20Λ | EESTEP_CLOSE §5; HK_RESIDUAL §4 |
| H_K [D]; EE kernels = finite sums of shifted scalar partial waves; J₀ < 3/2, 5/2; margins ½, 3/2 | EESTEP_CLOSE §0–§3, §6; HK_RESIDUAL §0–§1 |
| I1: DV on free-field ⟨TTTT⟩ | HK_RESIDUAL §2 |
| flat ≥ 1/5 (Ω^EE) [certified], attained to 10⁻⁸; ≥ 0.18 (Ω^Q) on most states, min 0.1835 | CERT_CLOSE §0–§2 |
| binding ≥ +0.198, ×78; ≥ +0.146, ×3.5 at 300; J = 0, 2, 4, 6, 8 [NE±] | NEARGAP2 §0; POS_ERRORS §3.1 (J = 0, 2) |
| large mass m ≥ 1.3Λ, b < 5: ≥ 8.7×; no m_* | POS_ERRORS §2.4 |
| Mehler–Dirichlet curvature bound [P] | POS_ERRORS §2.1 (Lemma C) |
| light J = 0, 2 positive to Δ_g − 1 (leading Regge) | NEARGAP2 §0 |
| H_RS: physical region = D′; O(Δ^N e^{−0.336Δ}); scalar 10⁻⁸ at Δ = 20 | PCYC_CLOSE §0–§1; RD_CLOSE Prop. D′ |
| ABJM inside A1′, A2; 2.6 % / 3.3 %; one digit | A1PRIME §7 |

**v4 additions.**

| letter item | source |
|---|---|
| Fig. 1 curves; τ₀ = 2.89126·10⁻⁸ / 3.13206·10⁻⁶; 2L_AdS = 92.8135 / 72.0406; (2b) 164.51 / 125.71 | `fig/make_fig.py` from `data/nf02_psi.json`, `data/q3star_factored.json` (prints them); = UT3 §3, REPRO |
| \|ĝ₃\|Λ⁴ ≤ 119.75 / 103.61; transplanted CHLPS optimum √(24.9 ln(2Λ/3) − 27.6) = 11.59 / 10.21; ratio 10.3 / 10.1 | `make_fig.py`; RESULT_D.md (route A′ comparison, "≈ 10× weaker", log coefficient 121×) |
| energy-flux disk ⇔ \|ĝ₃\|Λ⁴ ≤ Λ⁴/60 | τ ≤ 4, τ = 240\|ĝ₃\| (DICTIONARY_T4 §0) |
| τ = 240\|ĝ₃\|: t₄ half [DV], α_T half [D] with sign convention [I] | DICTIONARY_T4.md:27–29 |
| (1100, 100): −2154 ± 43 vs 1409, leading-Regge (D_h) shaped excess; refined form | HWIN_THIMBLE (via RESULT_D.md update; fw −2123 ± 172); HLJ_CLOSE §5–6; UT3 §5.2 |
| J = 29: −3.62 (plain chain, one K₂ outlier at x = 0.8) / +1.39 ± 0.03 (direct nodes); other odd J ≳ 27 possibly alike | NEARGAP2 §0, §1.6 (lines 36, 42, 230–235) |
| J = 12, 16 minima +0.6759, +0.7191 | supplement Table S6 (NG2) |
| Q: +0.143, s_d/s_fail = 0.51/1.776 (factor 3.5); ≥ 18× (11/0.59 = 18.6); [NE±] on 1.3–3Λ at 1000 | UT3 §5.4; PE §2.5 |
| c₁ 1 % thresholds at 1000: 59.6 (1), 1306 (2a), 1262 (2b) | UT3 §3 |
| ABJM: (Λ/ML)⁶ shift ≈ −3 % at 1000; λ = 2.51·10⁷ / 3.15·10⁹ for Λ = ML | UT3 §8; `(Λ/2^{5/4}π^{1/2})⁴` (referee 2 m12) |
| CHLPSpw = D ≥ 5 graviton partial waves; "record new bounds on the central charges of holographic conformal theories" | arXiv API abstract of 2205.01495 (checked 23 Nov 2026) |

**v5 additions and replacements** (rows above that they supersede: the J = 29 row and the J = 12, 16 minima).

| letter item | source |
|---|---|
| binding J = 5: +0.2009, ≥ +0.198, L = 78 [NE±] (unchanged) | NEARGAP2 §0 |
| [NE±] near gap: J = 0, 2 (L ≥ 6), 4 (≥ 22), 6 (≥ 5.0), 8 (λ-level 4.3–8.4; sign-free 1.1–1.8) | NEARGAP2 §0 (Letter: "from 1.7 to 22 for J ≤ 16, J = 8 on a less conservative budget") |
| J = 12: +0.8528 / +1.0006 / +0.6872, L 7.4 / 7.5 / 3.6; J = 16: +0.9402 / +0.7681 / +0.7509 / +1.0009, L 4.4 / 3.7 / 2.6 / 1.7 | HIGHSPIN §4 (values); SMALLNU_WINDOW §0, §3.2 (leverages with the measured window) |
| J = 21: +0.998 ± 0.020, L = 1.20 (×3) / 1.26 (signed); 1.7–6.2 at 1.03–1.10Λ | SMALLNU_WINDOW §0, §3.3 |
| J = 29: +0.972 ± 0.018; L = 1.23 (signed window), 0.94 (sign-free ×3 envelope); 2.2–5.7 at 1.04–1.10Λ; −3.62 = reference-zero node, +1.39 = spline aliasing | HIGHSPIN §2; SMALLNU_WINDOW §3.1 |
| J = 31: twist edge +1.20 → +1.02; J = 32 (1.0315Λ) +0.904 → +0.997 | HIGHSPIN §2.1 (K = 4 values; K6: +1.208 → +1.019, +0.785 → +0.945) |
| J = 23, 27, 31 [NE + M], limited by K₂ truncation (1.48 / 2.95 / 9.47 N_s); calibrated L ≈ 1.5 / 1.4 | SMALLNU_WINDOW §3.3; HIGHSPIN §3, §4 |
| ν < 50 window measured; convention kept on 6–18 % of the weight, ×3.6 above the largest known value | SMALLNU_WINDOW §0, §1, §4.4 |
| H_LJ(1000) at (1100, 100): central ≥ 1.6; A+C ≥ 1.0 at 2σ; λ-free A+B+C 9 nodes 0.5–0.6 at 2σ | DIRCV §8 (table and "Net") |
| D_h excess = a_h normalisation [DV]; congruence lemma [P] | N2_SUBREGGE §1 / UT box; LARGESPIN_W W1c |
| H_R′ [M] → [D-sketch + M] at b̂ ≥ b̂_A(m) = 8.5–15; ROUTEA window only m ≥ 1.5Λ | LARGESPIN_ANALYTIC §0.1, §5.2, N4 |
| leading-Regge floor rigorous at all b̂ for m ≥ 1.5Λ; b̂ ≥ 1.0421Λ for every m ≥ Λ (Lemma LR∞; r ≤ 0.577 …) | LARGESPIN_L §0.1 L1, L1′, §3 |
| magic-block identity for every b-shifted block, [P] to r⁸, [DV] to r⁴⁰, [D] all orders; monopole D = −9.98332·10⁻⁴ | PW_IMAGES §0.1 PW1–PW4, §2, §4 |
| Theorem W: J\* = 50 (100 at ×2) at (Λ_amp, K_sh, K₄) = (3, 10, 10) without P3; tolerance table | PW_IMAGES §5; LARGESPIN_W §4, §9 |
| Theorem L: J\*(m) = ⌈(1.0421m − 1)/2⌉, 857 at the gap; Λ_tol 1240 … 43 (J = 857 … 10⁶) | PW_IMAGES §6; LARGESPIN_L §3.3 |
| Λ_amp ≈ 2–3 [NE] from the J = 100 fits; model Λ_B = Λ_C2 = 1, Λ_A = 3.17, Λ_C3 = 2.03 [M] | LARGESPIN_W C2; P2_SPINNING §0, §3 |
| P4 gap: flat ν⁶ term −(η²−1)ν⁶P/72 at O(1/m⁴), subtraction [I] (P2FOLLOW) | P2_SPINNING §7.1 (in progress when read) |
| EE-kernel reduction and Lemma Y reproduced from statements | INDEP_CHECK §0 |
| τ₀/4 < 10⁻⁸ (conformal-collider comparison) | arithmetic: 2.89·10⁻⁸/4 = 7.2·10⁻⁹ |

All other letter numbers: unchanged from v2/v3 (UT2 §12, UT3 §12).

## 5. Referee reports

**Second report (`REFEREE_REPORT_2.md`, reject for PRL as written).** Point-by-point response with file:line: `CHANGES_v4.md`. Summary:
§6 items 1 (narrow/rewrite + figure), 4 (positioning), 5 (lower end dropped from abstract/discussion), 8 (housekeeping) done; M4 grade
slips fixed; M1 weak points now named in the main text; fallback written (§8 below). Not addressable by writing: §6 items 2 (close
positivity at one Δ_g), 3 (auditable derivations, independent reproduction of the EE reduction and Lemma Y), 6 (c₁), 7 (re-optimised Q).

**First report: status after v3** (kept for the record).

| point | status |
|---|---|
| M1 theorem/[D] too strong | addressed (v2); v3 net grades: (1) [D-sketch]; (2) [D-sketch] + [A] via H_EFT(b); **no [I] input**; positivity not established. |
| M2 hypotheses as labels | addressed (v2); v3 adds the [NE±] grade with explicit leverage, and Table S5 separates error-bounded from model-only states. |
| M3 k = 4 light side | addressed; **remaining items closed in v3**: a = b integrability (L6′ [D]) and spinning H_EFT(a) ([D] for dispersive terms). |
| M4 lower end | addressed (no constant). |
| M5 c₁ | addressed; v3: universal, shock route identically zero, zero given F₂ [A]; sensitivity thresholds stated. **Remaining:** F₂ or a direct value. |
| M6 de-shift, weights | **fully addressed in v3**: X^Q published; continuum k = 2 weights for both functionals; C4 ε_W interval-certified. Remaining minor: C5 ε_W not certified (not load-bearing); Q3\* factored coefficients not a separate data file. |
| M7 Lemma C / Einstein instance | addressed (v2). |
| M8 H_K | **closed in v3**: H_K [D] (EE step derived; f_pos(±3i/2) = 0; criteria 3/2, 5/2). |
| M9 certificate; (1100, 100) | **certificate closed in v3** (≥ N/5 everywhere); (1100, 100) still stated; H_LJ(1000) awaiting HIGHSPIN. |
| M10 H_RS | **closed in v3**: the physical chain is D′ by CLS-D's definition; P_cyc no longer a premise; census claim withdrawn. |
| M11 Δ_g = 300 | addressed (v2); binding state at 300 now [NE±], L = 3.5. |
| M12 title/scope | addressed; v3 replaces A1 by A1′, which ABJM-IIA satisfies at leading order. |
| M13 ABJM | v3: an in-scope consistency check, explicitly not evidence for positivity. |
| M14 novelty/credit | addressed (v2); v3 adds Paulos, Costa–Gonçalves–Penedones, Penedones, Jain et al. |
| M15 jargon, derivations | addressed (v2); v3 adds derivation sketches for the EE kernel, L6′, H_EFT(a), c₁, light states, D′. |
| minor 11 front matter | **done in v5** (authors' data via the coordinator): both authors "Independent researcher", e-mails ansh.mishra2025@gmail.com, aryan.senthilk@gmail.com; no acknowledgments (section removed). |
| overall verdict | positivity is still numerical; v3 makes part of it error-bounded and states exactly which states are model-only. |

## 6. Open items before submission

1. ~~HIGHSPIN~~ **done in v5**: places filled, `\HSdraftfalse`; no state negative, so §8's fallback was not triggered.
2. ~~Front matter~~ **done in v5**: affiliations "Independent researcher", e-mails in both files; acknowledgments removed (none).
2c. **v8: (AG_U) + (EDGE-R)** replace (AG′) (AG_PROOF: (AG′) not proved and pointwise false near the gap edge; `% AGPROOF` marks
   resolved, CHANGES_v8 §1). Open: (AG_U) beyond all orders [A]; (EDGE-R) for the graviton [M] — needs real-axis graviton values on
   x ∈ [0.94, x₁] at the J ≤ 51 gap-mass states (external job list `ads4/ext/jobs_batch3.txt`; the n3 C3H evaluator is broken, EDGE_R
   §2.4); the 20 states with m ≥ 1.13Λ need only n1. If a graviton E_χ exceeds Table S9's tolerance at (26, 1.026Λ), the sampled J\*
   rises (38 if ≤ 1.2 S_E,t). Also open (v8): (E′) for light single traces (DSKETCH_UP §3.3) and the spinning δ(Δ_L) for LS6.
2b. (superseded by 2c) **v7: (AG′)** (beyond-all-orders analytic structure of exact − flat on ν ≥ 3/β; Supplement §S7.D) was the one structural
   assumption of both large-spin theorems; a proof was attempted (`ads4/AG_PROOF.md`). The magnitude premises (W: spinning O(1/m⁴) window coefficients [M]; EDGE: graviton edge amplitude [M], +− factor;
   P_L′(i) above J = 344 [M]) and W's P0 below b = 245 [NE] remain. Also open: the J = 23 and 31 twist edges (noise-aware representation
   term), J = 29 sign-free (resolve the ν = 20 node), certification (Arb) of Lemma P4W's constants, the J = 0 engine comparison.
2a. ~~**P2 follow-up**~~ (done: P2_CONTRACT, P4_WINDOW; v6's `% P2CONTRACT` marks resolved in v7) (running when v5 was written): K₄ to O(1/m⁴), the (S) check against 1809.05111, and the O(1/m⁴) hybrid
   subtraction behind premise (P4) of Theorem W (P2_SPINNING §7.1 found the flat ν⁶ term). Update every `% P2FOLLOW` place
   (letter.tex Large-spin paragraph; supplement §S7.C, Table S9 and the "Evidence on the premises" paragraph).
3. ~~Publish the Q3\* factored coefficients~~ **done in v4** (`data/q3star_factored.json`).
4. Optional: certify C5's ε_W (root isolation of the quintic pole).
5. Before submission, the INSPIRE sweep recommended in LITERATURE.md §6 (v4 added CMT 1710.03199, Häring–Zhiboedov 2202.08280,
   Hartman–Kundu–Tajdini 1610.05308, verified on arXiv and INSPIRE).
6. Unreviewed deliverables: every analysis behind v3–v5 except REPRO / REPRO_ENGINE / INDEP_CHECK is unreviewed (referee 2 M3).
   v5: INDEP_CHECK reproduced the H_K EE-kernel reduction and Lemma Y (count and correlator level) from the statements; not covered:
   the correlator → light-function step of L6′, the non-dispersive counterexample. Still wanted: full write-ups of H_K, L6′, the R⁴
   calibration and c₁; an independent check of PW_IMAGES' identity (cheap: sympy to r⁸) and of Lemma LR∞; interval evaluation of
   LR∞'s constants (now floating point). **v7:** INDEP_CHECK2 reproduced the identity, D, LR∞ and W's algebra; LR∞'s constants are
   now certified (LRINF_TIGHT); LRINF_TIGHT, P2_CONTRACT, P4_WINDOW, EXT_FOLD are unreviewed.
7. `refs.bib` (unused legacy) still has the old JHEP fields; delete it or mirror `refs_verified.bib` before submission.
8. UT3 still carries the stale "natural size 10⁻³–10⁻¹" for c₁ (§3); the paper uses no size estimate. (v8: universality is [D] again,
   now with the explicit ASD solution of DSKETCH_UP, so UT3's "universality [D]" agrees with the paper.) UT3 was not edited.

9. **v10 (all gaps):** open for the all-gap statement: the seam spins J = 9–11, 13–15, 17–20, 22–26, 28, 30–32 and the parity-odd
   OO states J = 4–10 at Δ_g = 2000 (`% SEAM`, GAP_SEAM.md); the batch-4 direct rows for J = 12, 16, 21, 27, 29 at 2000 (`% BATCH4`;
   J = 21, 27, 29 are [NE + M] at their 2000 twist edges until then); a proof of the no-growth premise B_F / H_mono-A on
   1000 ≤ Δ_g < Δ₁\* (≈ 2 decades); certified graviton amplitudes at J ≥ 4 and ξ-interval coverage for the eventual-positivity
   certificate; Ω^Q_up at any gap other than 300/1000; the monotonicity remark of §S8.A (new wording, [D]) should be confirmed by the
   authors; GAP_* analyses unreviewed.

## 7. Source inconsistencies

**v10**: listed in `CHANGES_v10.md` §6 (gap-unit constants grow, not shrink; "certified above 10⁵" only under premises and for the
graviton only at J = 0, 2; Lemma T vs B_F wording; J = 8 L = 5.2 vs 3.9; budgets ∝ Λ⁻² only to J = 16; L′ tolerances recomputed in
GAP_LARGESPIN although GAP_SCALING says not; the v9 [NE] floor hole upgraded; UNIFIED_THEOREM.md not edited).


**v8**: listed in `CHANGES_v8.md` §3 (AG_PROOF's "exact scalar ≤ 0.30 S_E over all cutoffs" and "plain margin ≥ 520" corrected by
EDGE_R; the 1.35 vs 1.42 margin; UT box for DSKETCH_UP vs its §0 table; HK_CLOSE §5/§6 vs NU_TAIL_CHECK; k2/README "θ changes nothing
computed"; NOTES' own v7 open item 8).

**v7**: eleven items, listed in `CHANGES_v7.md` §4 (J\*(1.3Λ) = 242 not 241 in LRINF_TIGHT §6; 40 vs 47 W states; the J\* = 64
vs 26/42 lineage; EXT_FOLD verdict vs addendum for J = 31; magic-identity residual normalisations; Theorem L tolerance 43 (numerical
floor) vs 32 (certified); P_L′ status; K₂ calibration 4–5 % vs 0.6–8 %; CHMRS/CLS-D citation labels; the 3/16 residual; a stale UT
box).

**v6** (found while answering referee 3; the paper's choice is stated):
1. Letter v5 wrote B_eff = −19.265 for both gaps; with the published weights B_eff = −19.2646 (1000), −19.2197 (300); the single value
   gives 72.00 at 300, not 72.04. Fixed (Letter Eq. (3) text; Supplement already had both).
2. Letter v5 spin lists overlapped (J = 26, 36, 40 both "computed in models" and "interpolated"); Table S7 said "even J = 26, 30,
   34–46 without rows". Resolved from POS_ERRORS §3.2/HIGHSPIN: computed (chain) 36, 40, 48; interpolated 26, 30, 34, 38, 42, 44, 46.
3. Supplement graded Theorem W's P0 "curvature [P + DV] for b ≤ 119"; that bound needs m ≥ 1.3Λ, so P0 is [NE] in W's region.
4. Supplement "residual ≤ 8·10⁻⁴⁸ (also b = 7, 10)": 7.9·10⁻⁴⁸ for b ≤ 4; b = 7: 2.0·10⁻⁴⁷; b = 10: 1.1·10⁻⁴⁷ (`pw1_magic.txt`).
5. Letter v5 "|c₁| ≳ 1260 (item 2)" is the (2b) value; (2a) is 1306.
6. Letter v5 "measured … only at J ≤ 100" for Λ_amp: it is one state, (1100, 100).
7. P1 "[D given (S)]" (v5 supplement) is superseded: (S) fails (P2_SPINNING §7.2); P1 [D-sketch].
8. REGGE_TREE: the referee's term-wise corner failure also hits k = 5, not only k = 2, 3 (stated in the Supplement's "Term-wise use").

**v5: found between the deliverables read for this pass (not silently resolved; the paper's choice is stated).**
1. The coordinator's brief said "near-gap J = 4–8 [NE±]"; J = 7 is [NE + M] (HIGHSPIN §6, NEARGAP2). The paper lists J = 0, 2, 4, 6, 8
   (+ binding 5), 12, 16, 21, 29.
2. HIGHSPIN §0 / negative finding 2 and its UT box say removing the reference-zero node *raises* J = 31, 32 ("all shifts favourable"),
   but J = 31's twist edge goes +1.208 → +1.023 (down; still positive and above the published minimum +0.988). HIGHSPIN §0 also pairs
   "+1.208 → +1.023" although §2.1 gives K = 4: +1.200 → +1.023 and K6: +1.208 → +1.019. Paper: "+1.20 → +1.02", no "favourable".
3. J = 31 twist edge has three values: published +1.20, reference-zero-free shape chain +1.02, HIGHSPIN's hs_odd spline-edge run
   +0.786 (SMALLNU: +0.841 after window removal). Paper quotes +1.02 and notes "+0.84 with a spline-edge representation".
4. Name clashes: LARGESPIN_ANALYTIC's "Lemma C′" (curvature bound to b̂ ≈ 119) ≠ the Supplement's Lemma C′ (LEMMAS_CLOSE); LARGESPIN_L's
   "Lemma W" (Paley–Wiener locality) ≠ the Supplement's Lemma W (no roots in [1, ∞)). The Supplement calls the former "an extended
   Mehler–Dirichlet curvature bound" and does not use the latter.
5. DIRCV's first UT box still says "leverage ≤ 1 excluded at ≥ 2σ in every fit"; DIRCV §7 and §8 withdraw it (λ-free fits 0.8, then
   0.5–0.6 at 2σ). Paper follows §8.
6. LARGESPIN_L §8 inventory calls near-gap 150 ≤ J < 857 "still numeric/model"; PW_IMAGES/Theorem W′ covers sampled J = 50–344 at
   m ≤ 1.3Λ (conditionally). Near-gap 345 ≤ J < 857 (m ≈ 1.3–1.65Λ) is in neither theorem's sample; the paper assigns it to H_R′
   [D-sketch + M] (Table S11).
7. Theorem W's P4 as stated in LARGESPIN_W/PW_IMAGES bounds the whole O(1/m⁴) term; P2_SPINNING §7.1 (in progress) finds a flat ν⁶
   term that must be excluded. The paper states P4 as "the part not already resummed in the floor" and flags the subtraction as [I]
   (P2FOLLOW).
8. HIGHSPIN's caution "the ν < 50 window is now the binding model input at large J" is superseded by SMALLNU (binding input for
   J = 23–31 is the K₂ truncation); HIGHSPIN's odd-J window weight omitted a gauge factor (8.52 vs 8.73·10⁻⁵).
9. The update boxes in UNIFIED_THEOREM.md run oldest → newest (PW_IMAGES last), not newest first as the brief said.

**Earlier (v3/v4):**

UT3 §11: UT2's twenty items (17 resolved; #9 RD fit 2 % vs 3 % and #20 CHMRS equation numbers still open in the record; #1 procedural)
and fifteen new ones (N1–N15). Items that changed letter text in v3: H_K k = 4 criterion back to 5/2 for our weights; certificate now
1/5 [certified]; Q flat minimum 0.1835 (not 0.1966); certified ε_W; ABJM inside A1′; H_RS by definition; L6′ [D]; H_EFT(a) restricted
negative finding. One grading tension is left as found: PCYC_CLOSE calls H_RS "[D]" with one [D-sketch] sub-step (UT3 N13); the letter
writes "[D, one step D-sketch]".

## 8. Fallback plan if HIGHSPIN returns a negative state (referee 2, M1(d))

**Status v8: not triggered** (no new state values; AG_PROOF found (AG′) false in its pointwise form, but its replacement (AG_U) +
(EDGE-R) keeps J\* = 26 given (EDGE-R); the live trigger is now a graviton edge integral above Table S9's tolerance).
**Status v7: not triggered** (EXT_FOLD: no state negative; the lowest new value is J = 7's twist edge, +0.268). The P2 trigger in
the next paragraph is obsolete: the O(1/m²) corrections are computed and use ≤ 0.41 H. A new trigger: a proof that (AG′) fails, or a
graviton edge amplitude above Theorem W's tolerance (Table S9), would remove large-spin support without making a state negative.
**Status v5: not triggered.** HIGHSPIN and SMALLNU_WINDOW found no negative physical state (both J = 29 negatives were chain
artefacts). The plan below **remains valid** as the response to any future negative state of Ω^EE at Δ_g = 1000 (e.g. from the
untested premises of Theorems W/L, the near-gap band 345 ≤ J < 857, or J ≥ 150 window data), with one update: the trigger list should
include "a P2 computation giving Λ_amp beyond the Theorem-W tolerance (≈ 23 at J = 50, ≈ 27–31 at J ≥ 150, Table S9)", which would
not by itself make a state negative but would remove the large-spin support. The "If HIGHSPIN closes positively" paragraph at the end
applies: weak points (i)–(iii) are kept in the Letter, rewritten with the new numbers, because the budgets for J = 21, 29 are thin and
H_LJ(1000)/H_R′ remain [D-sketch + M].

Worked out from UT3 §3 ("Dependence"), §5.2–§5.4 and §5.6. Trigger: any heavy state of `Ω^EE` at Δ_g = 1000 comes back with
`λ_min < 0` (the J = 29 edge resolves negative, a near-gap J ≥ 12 state is negative, or the refined H_LJ(1000) fails on a measured
window so that some large-J near-gap state is negative). Then (Pos) for `Ω^EE` is **false** at 1000, and:

| statement | fate | why |
|---|---|---|
| item (1) at 1000 (τ ≤ 2.89·10⁻⁸; \|ĝ₃\|Λ⁴ ≤ 120) | **withdrawn** | its only positivity input is (Pos) for Ω^EE |
| item (1) at 300 | **withdrawn** | already rests on an open window [I]; beyond-Regge corrections at fixed (J, m/Λ) are larger at 300 (∝ Λ⁻² or faster), so a failure at 1000 is expected to recur at 300 for the same class; no positive evidence could rescue it |
| item (2b) (164.5 / 125.7) | **withdrawn** | uses Ω^Q_up + 0.040 Ω^EE and needs (Pos) for Ω^EE |
| Λ-forms of item (1), the conditional remark on Δ_g ≥ 300 | **withdrawn** | values of a now-vacuous implication |
| **item (2a) at 1000** | **survives unchanged** | uses only Ω^Q_up; the exact c₁ form needs only κ > 0, not item (1) |
| item (2a) at 300 | survives formally, not headlined | Q's evidence at 300 is weaker (large mass not converted, L < 1 at m ≥ 3Λ, b ≥ 5) |

**What the Letter would then claim (precisely).** *Assume A1′, A2, the inputs of Table S2 used by item (2) (H_K, L4–L9, LS rows), and
(Pos) for Ω^Q_up on heavy states and on light spin-0, 2 states. If no scalar single trace with m < 0.4602Λ appears in T × T, and either
τ = 0 or κ > 0 [D-sketch + M], then at Δ_g = 1000*

    g̃₄Λ⁶/8πG ≤ 92.81 + 5.44·10⁻⁷ c₁²   (and ≤ 72.04 + 6.81·10⁻⁵ c₁² at Δ_g = 300),

*with g̃₄ = g₄ + O(Λ⁻²) given H_EFT(b) [A]. Implication [D-sketch].* Its positivity evidence (UT3 §5.4): binding pair (++/OO, J = 0 at
the gap) +0.1427 / +0.1428 with s_d = 0.508 against s_fail = 1.776 [NE + M given H_RS]; measured channels min +0.1426 ± 0.0017; unmeasured
channels s_fail ≥ 11 (≥ 18× the measured |s_d|) [NE + M]; 1.3Λ ≤ m ≤ 3Λ, J²/Δ ≤ 1: [NE±], L ≥ 3.95; m > 3Λ: L ≥ 1.65 [NE + M];
J²/Δ > 1 and the large-J near-gap tail: [M] (H_LJ^Q, tolerance ≥ 1.1·10⁴ c₂, margin ×1–10³); HYB [M]. Flat: ≥ 9/50 certified on most
regions, two regions [NE]. Note: with item (1) gone, τ is no longer bounded by the paper, so the sign κ > 0 [D-sketch + M] becomes
load-bearing whenever ĝ₃ ≠ 0 (κτ²Λ⁸/57600 can exceed the upper end by orders of magnitude), and must be stated as such.

**Checks to run before using the fallback.** (a) Evaluate Ω^Q_up on the failing Ω^EE state (same K₂/direct rows, Q's weights): Q's
smeared k = 2, 3 parts use the same EE structure, so a failure mechanism of Ω^EE at near-gap large spin could also act on Ω^Q_up, though
Q's sign-free tolerance there is much larger (≥ 1.1·10⁴ c₂ vs T_B = 1409 c₂ for Ω^EE’s window at (1100, 100)). (b) Do **not** re-tune Ω^EE to absorb the failure before
submission: that adds a further round of post-hoc selection (referee 2, M1(d)); if it is done, it must be disclosed and re-checked on all
computed states, at both gaps.

**Framing.** Title: "A conditional bound on an R⁴ graviton coupling in AdS₄/CFT₃ from ⟨TTTT⟩ light-ray sum rules". Abstract: only
(2a) at 1000, with its condition. The cubic programme becomes a reported negative finding in the Supplement (the state, its value and
grade, and that the fixed Ω^EE is not positive at 1000). PRL suitability would drop further (the surviving result is the
IR-log-dominated g₄ window's upper end, in which ABJM sits at 2.6 %); the referee's PRD/JHEP route becomes the default.

**If HIGHSPIN closes positively** (J ≥ 12 at [NE±], J = 29 positive with a resolved edge, H_LJ(1000) upgraded): fill the `\HS` places,
set `\HSdraftfalse`, and remove weak points (i)–(iii) from the Letter's list accordingly; nothing else changes.
**If it comes back weak but not negative** (budget still not closed): the Letter stands as written; the weak points stay listed.

## AI-use disclosure (authors, 24 Nov)
The Letter ends with an acknowledgments block containing only the authors' AI-use disclosure, verbatim:
"AI (Claude) assisted with the formulation and drafting of this manuscript. The authors take full responsibility for the work presented, and have thoroughly reviewed all output."
No other acknowledgments. Keep this sentence in every future version (and in the arXiv version).
