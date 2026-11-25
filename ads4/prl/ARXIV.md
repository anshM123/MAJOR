# arXiv submission metadata (draft for the authors to approve)

Prepared 25 Nov 2026 for Letter v14 (Fig. 1 / Fig. 2 redrawn in the v15 presentation pass; no claim changed).
Copy each field into the arXiv form as given. Items marked **[AUTHORS]** need a decision.

## Title

How Close to Einstein Must Holographic Gravity Be?

## Authors

Ansh Mishra, Aryan Senthilkumar

(arXiv author field: `Ansh Mishra, Aryan Senthilkumar`. Both are independent researchers; arXiv has no affiliation field, so
affiliations appear only in the PDF.)

## Abstract (plain text, 576 characters; arXiv limit 1920)

The Letter's abstract, unchanged except that `$10^8$` is written as `10^8` (arXiv renders `$...$` with MathJax, so `$10^8$` also
works; plain text is safer in listings and e-mail digests).

```
Holography predicts that quantum field theories with a large gap behave like Einstein gravity, but only as a scaling law, with no numbers. We supply the numbers for three-dimensional theories: with a gap of 1000 or more, the graviton's three-point couplings match Einstein gravity to one part in 10^8 at leading order, and the first stringy correction is bounded from above, with ABJM theory inside the bound. The key input is causality, that energy detectors commute, together with a positivity property we prove in flat space and check numerically in anti-de Sitter space.
```

## Comments

```
4 pages + references, 2 figures; 34 pages of Supplemental Material. Data (functional weights) and code: <URL placeholder>
```

* The Letter PDF is 6 pages in total: the text and acknowledgments end on page 4, and the references fill pages 5 and 6 (REVTeX
  `prl` style, build of 25 Nov 2026). The Supplement is 34 pages. If you would rather count the whole PDF, write
  "6 pages, 2 figures, 34 pages of Supplemental Material".
* **[AUTHORS]** Fill in the data and code URL (e.g. a Zenodo DOI or a public GitHub tag). Use a permanent archive (Zenodo) if you
  can: the Letter says that every number can be reproduced from published data and code, so the link should not change.
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
  `fig/fig2_coverage.pdf`) rather than a PDF only. arXiv runs BibTeX, but including `build/letter.bbl` is safer.
* The `% BATCH4` comment in `letter.tex` (and the four in `supplement.tex`) would be visible in the arXiv source. Strip TeX
  comments before upload (e.g. with `arxiv_latex_cleaner`) or decide to leave them.
* The `\cite{SM}` entry reads "See Supplemental Material at [URL will be inserted by publisher]". For arXiv you may want
  "See Supplemental Material, appended" or the ancillary-file path. **[AUTHORS]**
