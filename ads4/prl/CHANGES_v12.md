# CHANGES v11 → v12 (24 Nov 2026): NEARGAP_FILL, GAP_SEAM, GAP_SEAM2 folded in

No physics was recomputed. Every number below was checked against the deliverable's cache (read-only):
`code/ngfill/cache/nf_summary.txt` and `nf_reprsum.txt` (J = 9–15 states, held-row estimate, K₂ calibration, J = 23/27/31 on λ_c),
`code/gapseam/cache/table.txt` and `gsm_bind_OO468.txt` (seam rows, OO block), `code/gapseam2/cache/gs2_bind_J27e/J27k/J29k.txt`, and
`code/gapnear/cache/gn_bind_J27.txt` (to reconcile the 2·10⁴ values, §4). The acknowledgments / AI-use disclosure are **unchanged,
byte for byte**. `ads4/UNIFIED_THEOREM.md` was not edited: its GAP_SEAM, GAP_SEAM2 and NEARGAP_FILL boxes are already at the end of
the box list (newest last) and agree with the deliverables (checked; see §4 item 7).

## 1. Letter

| place | change |
|---|---|
| Positivity paragraph (referee-5 M1) | "on thirteen further near-gap spins up to J = 31" → "near the gap on every other spin up to 16 and on spins 21, 23, 27, 29 and 31"; "Spins 23 and 31 are covered only above their lowest mass" → "Spin 31 is covered only above its lowest mass, spin 23 at its lowest mass only narrowly". `% BATCH5` mark added. |
| "Larger gaps" | The 2000 sentence now says what was computed (every spin up to 21 except 18, 20, and 27, 29), that up to 21 the checks close, "for spin 12 and from 14 on only with a signed estimate of the truncation error", that budgets shrink like Δ_g⁻² up to 16, that 27 and 29 are open at their lowest mass, and (new) "Several spins have their lowest values just above the new gap, with no counterpart at 1000" (the seam strip, GAP_SEAM N2). `% SEAM` marks removed; `% BATCH4` kept and narrowed (27/29 twist edges; sign-free closure for J ≥ 12). |
| Fig. 2 (`fig/make_fig2.py`) | Bars added for J = 9, 10, 11, 13, 14, 15 (twist edge to 1.05Λ); J = 23 now from its twist edge (to 1.10); J = 31 from 1.04Λ (was 1.05). Docstring sources updated; `# BATCH5` note at `NEPM`. Caption unchanged (still true). Regenerated. |
| Header comment | v12 summary and the mark list. |

Not changed: abstract (599 rendered characters), Intro, Discussion ("error control for the remaining near-gap states" still true),
acknowledgments.

## 2. Supplement

| place | change |
|---|---|
| §S6 intro | direct rows now "J = 7, 9–16, 21, 23, 27, 29, 31"; GK edges said for odd J ≥ 21 only (J = 9–15 rows end at x = 0.95 with the edge term). |
| **Table S5** | New row "J = 9–11, 13–15" (twist-edge values and L, min L per spin, NU_LOW = 25 variant, lowest value +0.3485 ≥ +0.327, hybrid base for J ≥ 13, held-row estimate). Row J = 23/27/31 rewritten: the **consistency correction** stated (representation term evaluated on λ_c, from which ν < 50 is removed and bounded separately; on λ it contained each probe's extrapolation into that window, i.e. counted it twice — J = 23: +0.120 of the K7 probe's +0.184; every probe and ×3 kept); J = 27 L = 1.81 (1.64 on λ), J = 23 L = 1.18 (thin; binding terms named), J = 31 twist edge L = 0.76 stays [NE + M] (real half-grid under-resolution), 1.04Λ L = 1.09; WLS 1.01 reported, not used. Row "other": the still-unbudgeted spins listed; J = 15 own-row clause removed. New row: OO J = 4–10 at 1000 [NE± given H_OO]. The table is now split over two floats ("TABLE S5 (continued)", one number) because it no longer fits a page. |
| Table S6 | J = 9, 10, 11, 13, 14, 15 direct values in parentheses (J = 14's chain minimum lies between budgeted masses; J = 15 published +0.6751, own K₂ +0.6900, direct +0.6943); OO J = 4 +0.215 → +0.206 at its twist edge. |
| Table S7 | H_M: where no own/direct row exists (mainly J = 17–20, 22, 24–26, 28, 30, 32; `% BATCH5`), calibration 0.4–8 % on 13 states; H_low: J = 9–15 measured ≤ 61 c₂ at ν ≈ 25, 50; H_OO: budgeted for J ≤ 10. Kept short so the table fits a page. |
| Table S11 | [NE±] row: J = 0, 2, 4–16, 21, 23, 27, 29 from the twist edge (mass ranges stated), J = 31 from 1.04Λ, OO J = 4–10; thin entries named (14, 15, 21, 23 with the λ_c correction, 29 signed). [NE + M] row: J = 17–20, 22, 24, 25 and the J = 31 twist edge (`% BATCH5`). |
| §S8 summary (2) | Rewritten for the 2000 seam/OO results. **Also fixes a v11 bug**: the words "(3) Between and beyond the computed gaps the" sat inside the `% SEAM / BATCH4` comment line and were missing from the PDF. |
| §S8.C and **Table S13** | Δ_ref per row stated. New rows J = 9, 10, 11, 13, 14, 15, 17, 19 (2000 twist edge, 1000 twist edge, 2000 hybrid minimum where it is elsewhere), J = 16 and 21 with Δ_ref = 5·10⁴ (J = 21 twist edge now λ-level [NE±], L 1.7), J = 27, 29 from GAP_SEAM2 (new J = 27 low +0.9098 at 1.0175Λ; twist edges [NE + M]; λ-level from 1.04/1.075Λ), OO J = 4–10 (sign-free). Caption says which value is tabulated (windowed central for J ≥ 9 and OO, own-row K₂ for J ≤ 8, 12). Readings (i)–(vi) rewritten: lowest values; sign-free vs λ-level; the reference artefact (λ-step fixed by 5·10⁴; at 27/29 edge and window terms take over; values move −0.02…−0.03); grid dependence of J = 27/29; OO identity and the J = 4 minimum; chart A incl. the x = 0.025 node effect; seam strip lows (+0.30…+0.53, J = 27 +0.910), so every larger gap needs its own seam rows. `% SEAM` resolved; `% BATCH4` kept. |
| §S8.D | "For J = 16, 21 the K₂ budget alone exceeds the margin at any gap" → "For J ≥ 14 the sign-free K₂ budget alone exceeds the margin at 1000 and 2000 (at 2000 only the λ-level budget closes)". |
| Table S17 | Near-gap rows split: J = 9–11, 13 + OO; J = 12, 14–16; J = 17, 19 (`% BATCH5`); J = 21; J = 23, 27, 29, 31 (`% BATCH4`); other J ≤ 32 (named; `% BATCH5`). Band row: chart-B images "mostly [NE + M]". |
| §S1 | Legend maps the new Letter phrases ("every other spin up to 16 and spins 21, 23, 27, 29 and 31"; "spin 23 … only narrowly" = L 1.18 after the correction; "signed estimate of the truncation error" = λ-level; "up to spin 21 the checks close …"; "lowest values just above the new gap" = reading (vi)). Provenance: NEARGAP_FILL rows and the λ_c evaluation listed as unreviewed. |
| §S13 (corrections) | v12 entries: J = 9–15 and the J = 23 twist edge upgrades, J = 27 L, J = 31 at 1.04; the earlier "noise amplification" description of the J = 23/31 representation term corrected; OO J = 4 minimum; calibration 0.4–8 % over 13 states; J = 21 at 2000; the reference-vs-node explanation of the budget scaling; the "+0.967" misprint in a source (not used here). |

## 3. Length and build

* APS count (`scratchpad/tools/wc_prl7.py letter.tex 3 0.9375 1.0385`, on a copy without the acknowledgments): prose 3006, 240 inline
  formulas (conservative 302), 3 display rows 48, figures 344: **3638** (conservative 3700); limit 3750, **headroom 112** (50
  conservative). v11: 3572.
* Abstract: not edited, 599 rendered characters.
* `build/letter.pdf` 6 pages (text now ends two lines into page 5, then acknowledgments and references); `build/supplement.pdf` 32
  pages. Both: 0 errors, 0 undefined references/citations, 0 overfull boxes, no "float too large" (Letter: 8 underfull notices).
  BibTeX: the known `jnrlst` notice and "empty year in SM".

## 4. Inconsistencies found (reported, not hidden)

1. **Two value conventions at 2000.** GAP_NEARGAP tabulates the own-row K₂ value; GAP_SEAM/SEAM2 the windowed central value λ_c. They
   agree to ≤ 0.001 for J ≤ 16 but differ by up to 0.03 at J = 21–29 (J = 27 tw₁₀₀₀ with 2·10⁴: K₂ +1.0305, central +1.0239). GAP_SEAM2's
   "2·10⁴ dense = GAP_NEARGAP's run" column is the central value of that run, not GAP_NEARGAP's printed number. Table S13's caption now
   says which is shown.
2. **GAP_NEARGAP's J = 21 "1000 same protocol" value** +0.967 should read +1.0009 (GAP_SEAM N5; its own json). The PRL never used it
   (Table S13's 1000 column is the graded direct value +0.998).
3. **OO J = 4 at 1000**: the twist edge +0.2059 is below the previously quoted minimum +0.2150 (that grid started at 1.01Λ). Still above
   the binding +0.2009. Table S6 corrected.
4. **J = 27/29 grades at 2000 depend on the x grid** (GAP_SEAM2 N2): J = 27 is λ-level from 1.04Λ only with edge nodes to x = 0.9 (from
   1.075Λ on the 1000 grid); J = 29 at 1.05/1.06Λ was λ-level with 2·10⁴ on the dense grid but is [NE + M] with 5·10⁴ on the 1000 grid.
   Table S13 grades the best variant and says so in reading (iv).
5. **J = 14's chain minimum** (+0.6982, near 1.035Λ) lies between the direct-row masses 1.03 (+0.7264) and 1.05 (+0.7946); covered
   by H_cont as for J = 16 before. GAP_SEAM's 2000 rows do include 1.035 (+0.6992, λ-level).
6. **v11 Supplement bug**: summary item (3)'s opening words were inside a `%` comment (fixed, §2).
7. UNIFIED_THEOREM.md: boxes for GAP_SEAM, GAP_SEAM2, NEARGAP_FILL already exist (added by the analysis commits), newest last, and
   match the deliverables; not edited here (this pass edits only `ads4/prl/`). The NEARGAP_FILL box's "0.4–7 %" is the six new states;
   the 13-state range is 0.4–8 %.

## 5. Marks left

* `% BATCH4` (external direct rows at 2000): Letter "Larger gaps"; Supplement §S8 summary (2), Table S13 comment, §S8.D, Table S17.
* `% BATCH5` (external direct rows at 1000, `ads4/ext/jobs_batch5.txt`): Letter positivity paragraph; Supplement Table S5 (rows J = 23/
  27/31 and "other"), Table S7 (H_M), Table S11, Table S17; `fig/make_fig2.py` (`NEPM`).
* `% SEAM`: none left (resolved).

## 6. Not done

* Rows at 2000 for J = 18, 20, 22–26, 28, 30–32 (GAP_SEAM2 stopped early); J = 23, 31 never run at 2000.
* Batch 4 and batch 5 (external).
* The Letter's Intro sentence "a sample of other low-lying states" left as is (still true; length).
