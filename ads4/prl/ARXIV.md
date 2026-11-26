# arXiv submission metadata (draft for the authors to approve)

Updated 27 Nov 2026 for Letter v17: positivity of the cubic functional is now established on every heavy state at Δ_g = 1000
(Theorem FP, Supplemental Sec. S11), so the cubic bound no longer assumes it; Fig. 2 replaced by the certificate cover.
Copy each field into the arXiv form as given. Items marked **[AUTHORS]** need a decision.

## Title

How Close to Einstein Must Holographic Gravity Be?

## Authors

Ansh Mishra, Aryan Senthilkumar

(arXiv author field: `Ansh Mishra, Aryan Senthilkumar`. Both are independent researchers; arXiv has no affiliation field, so
affiliations appear only in the PDF.)

## Abstract (plain text, about 690 characters; arXiv limit 1920)

A plain-text version of the Letter's abstract, written for listings (`10^8` instead of `$10^8$`; arXiv also renders `$...$`).

```
Holography predicts that quantum field theories with a large gap behave like Einstein gravity, but only as a scaling law, with no numbers. We supply the numbers for three-dimensional theories: with a gap of 1000 or more, the graviton's three-point couplings match Einstein gravity to one part in 10^8 at leading order in 1/c_T. The key input is causality, that energy detectors commute, together with a positivity property that we prove in flat space and establish in anti-de Sitter space, on every heavy state at a gap of 1000, by a computer-assisted proof. Under a further positivity assumption the first stringy correction is bounded from above, with ABJM theory inside the bound.
```

## Comments

```
4 pages + references, 2 figures; 43 pages of Supplemental Material. Code and data: https://doi.org/10.5281/zenodo.24571131
```

* The Letter PDF is 6 pages in total: the text and acknowledgments end near the top of page 5, and the references fill the rest of pages 5 and 6 (REVTeX
  `prl` style, build of 27 Nov 2026). The Supplement is 43 pages (build of 27 Nov 2026). If you would rather count the whole PDF, write
  "6 pages, 2 figures, 43 pages of Supplemental Material".
* Code and data: Zenodo, DOI 10.5281/zenodo.24571131 (filled in).
* If the Supplement goes into the same arXiv submission (recommended: append it to the Letter PDF or upload it as an ancillary
  file), say so: "Supplemental Material included".

## Categories

* **Primary: quant-ph** (the authors' choice).
* **Cross-lists: hep-th, gr-qc.**
* Note: hep-th is the conventional primary category for AdS/CFT and conformal-bootstrap work, and arXiv moderators may reclassify
  a quant-ph primary to hep-th.

## Other fields

* MSC class: none.
* ACM class: none.
* Journal-ref: blank (fill in after publication).
* DOI: blank (fill in after publication).
* Report number: none.
* License: **[AUTHORS]** choose one (arXiv's non-exclusive license is the default; CC BY 4.0 is compatible with APS publication).

## Checks before upload

* The acknowledgments must contain the AI-use disclosure verbatim (NOTES.md, "AI-use disclosure"): "AI (Claude) assisted with
  the formulation and drafting of this manuscript. The authors take full responsibility for the work presented, and have
  thoroughly reviewed all output." It is in `letter.tex`, so it is in the arXiv version.
* Upload the TeX source (`letter.tex`, `macros.tex`, `refs_verified.bib` or the generated `letter.bbl`, `fig/fig1_bounds.pdf`,
  `fig/fig2_certificate.pdf`; for the Supplement also `supplement.tex` and `sec_fullpos.tex`) rather than a PDF only. arXiv runs BibTeX, but including `build/letter.bbl` is safer.
* The `% BATCH4` comment in `letter.tex` (and the four in `supplement.tex`) would be visible in the arXiv source. Strip TeX
  comments before upload (e.g. with `arxiv_latex_cleaner`) or decide to leave them.
* The `\cite{SM}` entry reads "See Supplemental Material at [URL will be inserted by publisher]". For arXiv you may want
  "See Supplemental Material, appended" or the ancillary-file path. **[AUTHORS]**
* **Code and data link.** Theorem FP rests on published certificates (inputs, outputs, replay scripts) under `ads4/code/`
  (Supplemental Sec. S11.I lists the directories). The large coefficient tables are not in git (SHA-256 sums are); archive them with
  the code on Zenodo so the replay instructions work from the DOI.
