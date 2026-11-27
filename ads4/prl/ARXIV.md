# arXiv submission metadata (draft for the authors to approve)

Updated 27 Nov 2026 for Letter v2 (`letter_v2.tex`: thermal sum rules, spectral-to-gap theorem, Hall viscosity). Copy each field
into the arXiv form as given. Items marked **[AUTHORS]** need a decision.

## Title

A Higher-Spin Gap Constrains the Shear and Hall Spectral Weight of Three-Dimensional Conformal Field Theories

## Authors

Ansh Mishra, Aryan Senthilkumar

(arXiv author field: `Ansh Mishra, Aryan Senthilkumar`. arXiv has no affiliation field, so affiliations appear only in the PDF;
the PDF currently says "Independent researcher" **[AUTHORS]**.)

## Abstract (plain text, about 1304 characters; arXiv limit 1920)

A plain-text version of the Letter's abstract, written for listings.

```
The stress-tensor three-point function of a three-dimensional conformal field theory has two coefficients that vanish for Einstein gravity, the parity-even t_4 and the parity-odd alpha_T. We derive a Hall-channel counterpart of the known shear sum rule. Together, the two sum rules equate the thermal spectral weight of the stress tensor in the shear and Hall channels, measured from its Einstein value, to (t_4 + i alpha_T)/32, in theories without light scalars that have thermal expectation values; the conformal-collider bound then places every such theory in a disc of radius 1/8. At leading order in a large central charge, a gap Delta_g in the higher-spin spectrum shrinks the radius to about 0.56/(Delta_g - 2), using detector commutativity and no positivity assumption, so that a spectral deviation delta requires a higher-spin state below twist 2.25 + 5/(9 delta). This scaling is three powers of Delta_g weaker than the parametric expectation. At Delta_g = 1000, positivity of a dispersive functional, which we certify on every heavy state, reduces the radius to 9 x 10^-10, with an estimate of the light-state remainder that we sketch. In the effective theory of the dual black brane, the parity-odd cubic coupling produces a Hall viscosity without an axion, 4 pi |eta_H|/s = (3/40) |alpha_T|.
```

## Comments

```
5 pages, 2 figures; 57 pages of Supplemental Material (appended). Code and data: https://doi.org/10.5281/zenodo.24571131
```

* The Letter PDF is 5 pages (REVTeX `prl` style, build of 27 Nov 2026); the Supplement is 57 pages and is appended after the
  references in the arXiv tarball, so the arXiv PDF has 62 pages.
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
  thoroughly reviewed all output." It is in `letter_v2.tex`, so it is in the arXiv version.
* Upload `arxiv_submission_v2.tar.gz` (built from `letter_v2.tex` as `main.tex`, `macros.tex`, `main.bbl`, the two figures
  `fig/fig_disc.pdf` and `fig/fig_tau.pdf`, and the compiled Supplement, appended after the references).
* The `\cite{SM}` entry reads "See Supplemental Material at [URL will be inserted by publisher]". For arXiv you may want
  "See Supplemental Material, appended" or the ancillary-file path. **[AUTHORS]**
* **Code and data link.** Theorem FP rests on published certificates (inputs, outputs, replay scripts) under `ads4/code/`
  (Supplemental Sec. S11.I lists the directories). The large coefficient tables are not in git (SHA-256 sums are); archive them with
  the code on Zenodo so the replay instructions work from the DOI.
