# CHANGES_jargon: Supplement internal language and independence claims (Referee 6, P9 and P10)

Only `supplement.tex` was edited. `letter.tex` and every other file are unchanged. Prose only: no claim, number or grade changed.

## 1. Glossary (new, Sec. S1, after the grades paragraph)

The glossary is a compact in-text table (not a float) with 16 entries, each one line of plain language. It defines:
NF02 and Q3\* (and "flat floors"); twist edge, tw_Δ and tw_1000; binding state; hybrid action, hybrid floor, hybrid model and Γ; node and row;
K₂ row, direct row and own (own-mass) row; chain; held (H_hold); error budget and leverage L; λ-level budget; window and signed window;
edge term; reference row and "reference zero"; chart A and chart B; seam; H17 (and so Λ²H17). The new label `sec:grades` on S1 lets
later sections point to it.

## 2. Renamed or removed in S6–S8

| was | now |
|---|---|
| "crystal integral" (S8 summary, Table S12 "Fixed-phase structure") | "light-ray integral" |
| "physical chain" (S8, Table S12 "Cycle rigidity"); "$H_{\rm RS}$: the physical integration chain"; "tests the chain" | "physical integration cycle"; "$H_{\rm RS}$: from the definition to the expansion"; "tests these steps" |
| "fold chain", "fold protocol" (S1, S10, Table S18 last row) | "re-run of our code"; "assembly of the direct rows into values and budgets" |
| "seam strip" | "seam" (glossary term) |
| "own-mass rows" in the S8 text | "own rows" (glossary term) |
| "hybrid contour", "hybrid thimble contour" (these clashed with "hybrid action") | "mixed contour (real, then steepest descent)" |
| "external records" (structure-basis repair) | "the even-J direct rows were computed in a different basis …" |
| "the plain odd-J chain … held its value over the edge"; "published shape chain" | "the odd-J chain … continued its value over the edge"; "the published chain from the shape fit" |
| "the chain's despiking filter" | "the outlier filter used on chain rows" |
| "Protocol notes", "adopted protocol", "production s_d", "the protocol used for" | "Conventions", "the conventions used here", "the s_d used here" |
| "the protocol is that of Sec. S6 (the same budget code …)" | "the error budget is that of Sec. S6 (same code …)" |
| "engine envelope" (Table S5), "the engine reproduces" (S7.C) | "beyond-Regge envelope H_bulk"; "the code that computes these corrections reproduces" |
| "as the earlier Monte-Carlo-limited direct edge indicated" | "consistent with a less precise direct edge evaluation" |
| "that protocol change is not used" | "this alternative is not used" |
| "now certified", "A new lemma", "(new exact identity …)", "no longer links" | "certified", "A lemma", "(exact identity …)", "not links of the argument" |
| T0's symbol H (it clashed with Theorem W's 2×2 budget H) | H_hyb (T0 subsection only) |
| "mass entries … hold the own row"; "It holds the lowest values" | "are held at the own row"; "It contains the lowest values" |
| "reached by no transfer route" (S8 summary) | "reached by neither chart" |
| "fold of the saddle" | kept, with "the ν below which it leaves the real integration region" added |
| "A chain carries the grade of its weakest link" (S1) | "An argument carries the grade of its weakest step" |

The S6 opening ("How it is computed") and the paragraph on the J = 17–32 direct rows were rewritten in plainer sentences.

## 3. Lighter pass elsewhere

* **File names and project history.** Removed: `code/t0cert/`, `fig/make_fig.py`, `fig/make_fig2.py`, `ANGULAR_ROUTE.md`, `LADDER_ALLJ.md`,
  `ANGULAR_MATCH.md`, `ANGULAR_AUDIT.md`, `ANGULAR_PROOFS.md`, `code/angular/`, `code/ladder/`, `code/angmatch/`, `code/angproofs/`,
  the `README.md` errata pointer, "Stage-7 grid weights" (now "k = 4,5 grid weights"), and "audits".
  Code pointers now read "the reproduction code (`reproduce.py` and the scripts it calls)". `reproduce.py` calls the T0 certificate and the angular,
  angmatch and angproofs scripts.
* `Q3-full-delta0.2` is now "the stored linear-programming solution with margin δ = 0.2 (`q3_full_delta0.2.json`)". The fingerprint
  `80e0834a…33f2` is kept, and the text now says what it is (a hash of the coefficients, by which the flat certificate identifies its input).
* The errata pointer is now in words: the `note` field of `nf02_psi.json` and the ε_W values stored in the continuum C4 files are superseded by the
  Supplement's values. The files are unchanged so that their hashes stay valid.
* The published data-file names in Table S1 (`data/`) are kept.

## 4. Independence claims (P10)

* **S10** now opens with "Who performed the checks, and what 'independent' means here". Every check was a separate run of the AI system
  acknowledged in the Letter, made after the computation it checks, and none was made by a person outside the project. A check is called
  "independent" only if it wrote its own code from the written definitions and the published data, without importing or running project scripts.
  A run of our own code is called a "re-run": byte-for-byte agreement shows reproducibility, not correctness. The runs are listed as (a)–(f), each
  with its disclosed exceptions, taken from the reproduction notes:
  * (a) κ and the H17-band model were read from project code;
  * (c) the Theorem-W floors were read as data, and the LR∞ routine was read only after the run's own result was fixed;
  * (d) the J = 18, 23, 25 values are a re-run of our code, and the local node and structure routines were used read-only;
  * (e) T0 = a re-run of the certificate code, plus a hand re-derivation and a separately written evaluator;
  * (f) the angular-route check imported our kernel code read-only, so it is not independent of that code.
  "third/fourth checker" is now run (b)/(d).
* **S1** "Independent checks" is now "Checks by separate code" and points to S10. It now says that the J = 17–32 values were a re-run.
* **S7.B (T0)** says the certificate re-run had identical outputs and that this shows reproducibility only. The other "independent check" mentions
  in S4–S7 now read "separate run/code (Sec. S10)".

## 5. Verification

* **Numbers.** A script compared every number token outside `\texttt`, `\ref`, `\cite` and `\label`, in order (difflib), between the old and the
  new tex. There are **no deletions or changes**, only insertions: glossary text, "δ = 0.2", restated J values in S1 and S10, and "k = 4,5". The one
  removed token is the "7" of "Stage-7".
* **Tables.** All 21 original tabulars have the same number sequence, apart from that "Stage-7" → "k = 4,5" cell.
* **Grade tags.** The multiset of all 659 grade tags ([P], [D], [NE±], …) is identical.
* **Build.** `pdflatex` ×2: 0 errors, 0 undefined references or citations, and 2 overfull boxes (the same 2 as before). Citations are unchanged, so
  bibtex was not re-run. **Pages: 36 before, 37 after** (the glossary table adds about half a page).

## 6. For the authors to confirm

* The S10 sentence "no check was made by a person outside the project" and the "AI system acknowledged in the Letter" attribution are inferred
  from the reproduction notes (REPRO, INDEP_CHECK*, INDEP_CHECK_T0, ANGULAR_AUDIT). Please confirm this wording.
* The Letter (l. 248–249) still says "Separate code reproduces every constant, and a separate check rebuilt …". This matches S10, but the Letter was
  not edited.
