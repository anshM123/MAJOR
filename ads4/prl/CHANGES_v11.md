# CHANGES v10 → v11 (24 Nov 2026): cheap fixes after `REFEREE_REPORT_5.md`

No physics was recomputed except the checks of the k = 2 weights at Δ_g = 2000 and 5000 (M5, read-only from the project cache). The
acknowledgments / AI-use disclosure are **unchanged, byte for byte** (referee M6 is for the authors). `% SEAM` and `% BATCH4` marks
are kept (Letter "Larger gaps", Supplement §S8.C, Tables S13 and S17, and a new mark at §S8's summary item (2)).

## 1. Report 5, item by item

| item | change |
|---|---|
| **M2.1** Intro "carry it … to every larger gap" | "We do this at Δ_g = 1000. The evidence extends to larger gaps, at no better strength, if near-gap corrections do not grow with the gap and the large-spin premises do not change." |
| **M2.2** "We argue below that positivity … carries over" | "What supports positivity of Ω^EE at larger gaps is described below"; "gives a sharper bound" → "would give". |
| **M2.3** "The evidence … carries over to every Δ_g ≥ 1000" | "Above Δ_g = 1000 the evidence for Ω^EE is at best as strong as at 1000 and needs further premises" (§S8's net statement). |
| **M2.4/5** "We repeated the checks at 2000"; "21, 27, 29 remain thin" | "No state is negative, and up to spin 16 the checks close … Spins 21, 27 and 29 are not yet error-controlled at their lowest mass (at 1000 this needed direct evaluations), and the other near-gap spins have not been recomputed." §S8 summary item (2) says the same (23, 31 not recomputed; seam and OO without rows). |
| **M2.6** "checked at 1000, 2000 and 5000" | "At 2000 and 5000 we tested this only on the saddle-point evaluations used for the checks, and up to 10⁴ on an exactly computable scalar analogue." Legend (§S1) maps this to the own-mass K₂ rows and their structural gap independence. |
| **M2.7** "a certified bound … replaces it" | "supplements it for the lowest spins, given further premises"; the legend says the third premise is B_F for an explicit function. |
| **M2.8** "margins grow" | "If the large-spin premises do not change with the gap in AdS units, the theorems' margins grow." The "tenth of its margin / 140 times" numbers moved out (Table S16 has them). |
| **M2.9** Fig. 2 caption "every larger gap looks the same or better" | → "At larger gaps the near-gap corrections shrink like Δ_g⁻², and the large-spin regions move to proportionally larger J." |
| **M2.10** Fig. 1 caption "holds given the premise" | → "Above Δ_g = 1000 the solid cubic curve is supported, given the premises stated in the text, at the strength it has at 1000". |
| Discussion (M2) | "positivity at larger gaps without the premise …" → "unconditional positivity at large spin or at larger gaps". |
| **M3** unassigned band | New Table S17 row "near gap above the seam, 33 ≤ J < J_W" (J_W = chart-B image of J = 26: 52 at 2000, 132 at 5000): no transfer route, not covered, hypothesis only. New §S8.D paragraph "A band with no transfer route"; §S8 summary item (6) and net statement; both "Not claimed" lists. Letter: "Near-gap states with spins from 33 to about 0.026Λ are covered by no argument yet." |
| **M5** functional at 2000/5000 | New `data/gap/`: `d_abel_X_{2000,5000}_trunc.txt.gz` (the min-norm grid weights the near-gap runs used, `code/gapscale/cache/weights_*.npz`, ν ≤ Λ+10), `gap_weights_checks.txt`, `SHA256SUMS`; script `data/export_gap_weights.py` (read-only import of `code/deshift/d_abel.py`). Checks: discrete residual ≤ 5·10⁻¹³; off-grid de-shift residual ≤ 6·10⁻⁶ for t ≤ 14, failing beyond t ≈ 18 as for the 1000 grid; V3 to 10⁻¹⁴; Λ²dk₀ = −68.081 / −68.077; floor ≤ 1.5·10⁻¹⁴ beyond Λ+10. **Not done:** the continuum (strip-analytic) fit with an all-t ε_W, which H_K uses, at 2000/5000 (≲ 2.5 CPU-h per gap with `cert/k2_cont.py`). Table S12 row *Weights* states exactly this ("[DV] at 300–5000 (grids; continuum form with ε_W only at 300, 1000); asymptotic beyond"); Table S1 gains a row; `data/README.md` §8, `MANIFEST.txt`. Letter: "All weights are published (at Δ_g = 2000 and 5000 as tables only)"; "the leading-Regge ones wherever the de-shifted weight exists (we constructed it up to 5000)". |
| **m1** length | "Large spin" 316 → ≈ 225 words (L′ "leading margin certified on every state", W's "leading margin certified only at larger spins" and "given the edge bound below" left to §S7); "Larger gaps" rewritten with fewer numbers; Discussion de Sitter sentence and two connectors cut; Intro "sharp swampland" clause shortened; "we have not computed it" (said twice elsewhere) cut in Functionals. |
| **m2** prose | "Three ingredients make it tractable"; "Spins 23 and 31 are covered only above their lowest mass, and at J = 29 …"; the dangling κ clause split ("…are only sketched; the sign of κ also rests on a model"). |
| **m4 / R4-M1.5** κ | Abstract: "and, given a sign condition, one quartic coupling is bounded above" (paid for by "flat-space limit" → "flat limit"); Discussion: "(the quartic one if κ > 0 or τ = 0)". ABJM stays in the abstract (authors' instruction, v9). |
| **m5** x ∝ \|c₁\|/Λ² | → "x = 1.67×10⁻⁴\|c₁\|" at Δ_g = 1000 (Supplement §S4). |
| **m7** Table S17 J ≤ 7 row | adds "at the sampled masses (twist edge, 1.02, 1.05, 1.1, 1.3Λ), between them H_cont [NE]". |
| **m8** r_c | "0.023–0.052 (scalar) and 0.014–0.032 (graviton prefactors)". |
| **m9** Fig. 1 | Transplanted flat optimum drawn without the paper's status markers (`make_fig.py`, `status=False`); "(2b): light scalars allowed; needs both functionals"; key "solid: Δ_g ≥ 1000, given premises (text)". Regenerated; constants unchanged. |
| **m11** references | Added, each checked on its arXiv abstract page (title, authors, journal DOI; INSPIRE not reachable): Kravchuk–Simmons-Duffin 1805.00098 (light transform, Mechanism), Penedones–Silva–Zhiboedov 1912.11100 (with J₀ ≤ 1), Hartman–Jain–Kundu 1509.00014 (Regge causality, Intro). |
| **m12** literature on uniformity | §S8 opening: one "Literature" sentence citing CHMRS, CHLPSpw (2205.01495), Knop–Mazáč 2203.11170 and van Rees–Zhao 2210.15683 (arXiv only). Not in the Letter (length). |
| §6.2 caveat 1 | Table S12 Lemma T row: "Λ-independent" holds for the flat weights; the realised k = 2 weight differs by S/Λ². |
| Legend (§S1) | Re-mapped every new Letter phrase of "Larger gaps" and the Intro sentence to its grade; "narrowly at the lowest spins" → "narrowly". |

## 2. Not done (needs computation or an authors' decision)

* **M1, M4**: positivity at 1000 (near-gap J ≤ 25, EDGE-R, Ω^Q binding states), c₁, Ω^Q at 2000 — computations.
* **M3 fix by computation** (K₂ rows for 33 ≤ J ≤ 51 at 2000, Δ_ref ≥ 5·10⁴) — not run; the band is now stated as uncovered.
* **M5 continuum weights + ε_W at 2000/5000** — not fitted (see above); the construction of Ω^EE at general Λ is stated as "grids [DV] to 5000, asymptotic expansion beyond".
* **M6 / acknowledgments**: untouched (authors). The "unreviewed" labels in the Supplement are unchanged.
* **m3** title, **m4** ABJM in the abstract, **m6** "certified" leading-Regge for the exact vs published weight (≈ 10 words; not added for length), **m10** Fig. 2 inset in J/Λ, **m13** stripping source comments (marks still needed for SEAM/BATCH4).

## 3. Length and build

* APS count (`scratchpad/tools/wc_prl7.py letter.tex 3 0.9375 1.0385`, run on a copy without the acknowledgments, as APS counts):
  prose 2940, 240 inline formulas (conservative 302), 3 display rows 48, figures 344: **3572** (conservative 3634); limit 3750,
  **headroom 178** (116 conservative). With the acknowledgments counted: 3600 / 3662. v10 was 3656 / 3718 on the same basis.
* Abstract: 599 rendered characters by the NOTES §2 method (pdftotext, whitespace-normalised; ≈ 593 without the spaces pdftotext
  inserts in ⟨TTTT⟩ and √); limit 600. Tight: any further abstract words need a cut.
* `build/letter.pdf` 6 pages: the text, acknowledgments included, ends at the bottom of page 4; references fill pages 5–6.
  `build/supplement.pdf` 31 pages. Both: 0 errors, 0 undefined references/citations, 0 overfull boxes (Letter: 8 underfull notices).
  BibTeX: the known `jnrlst` notice and "empty year in SM".
