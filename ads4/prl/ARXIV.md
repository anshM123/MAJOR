# arXiv submission metadata (draft for the authors to approve)

Updated 26 Nov 2026 for Letter v19 (prose rewrite; concept figure added as Fig. 1). Positivity of the cubic functional is
established on every heavy state at Δ_g = 1000 (Theorem FP, Supplemental Sec. S11), so the cubic bound does not assume it.
Copy each field into the arXiv form as given. Items marked **[AUTHORS]** need a decision.

## Title

How Close to Einstein Must Holographic Gravity Be?

## Authors

Ansh Mishra, Aryan Senthilkumar

(arXiv author field: `Ansh Mishra, Aryan Senthilkumar`. Both are independent researchers; arXiv has no affiliation field, so
affiliations appear only in the PDF.)

## Abstract (plain text, about 800 characters; arXiv limit 1920)

A plain-text version of the Letter's abstract, written for listings (`10^8` instead of `$10^8$`; arXiv also renders `$...$`).

```
Holography predicts that conformal field theories with a large central charge c_T and a large gap behave like Einstein gravity, with corrections that fall as powers of the gap. For three-dimensional theories we determine the constants. With a gap of 1000 or more, the graviton's three-point couplings agree with Einstein gravity to one part in 10^8 at leading order in 1/c_T. The bounds follow from the commutativity of energy detectors, a consequence of causality, together with a positivity property. We prove it in flat space and, at a gap of 1000, establish it in anti-de Sitter space on every heavy state by a computer-assisted proof. Under a further positivity assumption the first stringy correction is bounded from above, and the ABJM value lies below the bound. Without positivity, energy-detector commutativity and crossing still bound the deviation tau by about 17.8/Delta_g at any gap; read in reverse, a measured deviation tau places a higher-spin state, and so the bulk string scale, below twist about 2 + 17.8/tau.
```

## Comments

```
4 pages + references, 3 figures; 43 pages of Supplemental Material. Code and data: https://doi.org/10.5281/zenodo.24571131
```

* The Letter PDF is 5 pages in total: the text, acknowledgments and data-availability statement end on page 4, and the references fill pages 4 and 5 (REVTeX
  `prl` style, build of 26 Nov 2026). The Supplement is 43 pages (build of 26 Nov 2026). If you would rather count the whole PDF, write
  "5 pages, 3 figures, 43 pages of Supplemental Material".
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
* Upload `arxiv_submission.tar.gz` (built from `letter.tex`, `macros.tex`, `letter.bbl`, the three figures and the compiled
  Supplement, which is appended after the references). It has no TeX comments.
* The `\cite{SM}` entry reads "See Supplemental Material at [URL will be inserted by publisher]". For arXiv you may want
  "See Supplemental Material, appended" or the ancillary-file path. **[AUTHORS]**
* **Code and data link.** Theorem FP rests on published certificates (inputs, outputs, replay scripts) under `ads4/code/`
  (Supplemental Sec. S11.I lists the directories). The large coefficient tables are not in git (SHA-256 sums are); archive them with
  the code on Zenodo so the replay instructions work from the DOI.
