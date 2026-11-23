# CHANGES v5 → v6 (after REFEREE_REPORT_3.md; order of limits from REGGE_TREE.md)

*23 Nov 2026. Writing pass. The only arithmetic redone: 2L_AdS and B_eff at Δ_g = 300 and 1000, from the published weights with the
integrals of `fig/make_fig.py` (= `code/repro/bounds.py`). Sources read: `REFEREE_REPORT_3.md`, `ads4/REGGE_TREE.md`, the update boxes
of `ads4/UNIFIED_THEOREM.md` (through REGGE_TREE), and for Fig. 2 and the premises `POS_ERRORS.md` (§2.4, §3), `NEARGAP2.md`,
`HIGHSPIN.md`, `SMALLNU_WINDOW.md`, `DIRCV.md` §8, `LARGESPIN_W.md`, `LARGESPIN_L.md`, `PW_IMAGES.md`, `P2_SPINNING.md` §7,
`code/pwimg/cache/pw1_magic.txt`. Files edited: `letter.tex`, `supplement.tex`, `macros.tex`, `NOTES.md`, `fig/make_fig.py` (and its
outputs). New: `fig/make_fig2.py`, `fig/fig2_coverage.pdf`, `fig/fig2_coverage_preview.png`, this file. Nothing outside `ads4/prl/`
touched; `P2_CONTRACT.md` / `code/p2con/` not touched; nothing committed. Line numbers refer to the files as of this writing
(L = letter.tex, S = supplement.tex).*

## 1. Major issues

| referee 3 | action | where |
|---|---|---|
| **M2** order of limits / J₀ for a tree-level result | REGGE_TREE §5 wording implemented. A2 restated (tree-level EFT light data up to o(1/c_T), uniformly in spin; sum rules at finite c_T; heavy side never expanded). Mechanism: "Regge intercept J₀ of the exact correlator; J₀ ≤ 1 in every unitary CFT … used at finite c_T as in CHMRS, CLS-D … planar intercept never enters". Functionals: "hence at every finite c_T … bounds are the c_T → ∞ limit … term by term would need a planar intercept below 3/2, which fails; we never do this". "Heavy" defined as every state with Δ − J ≥ Δ_g, single- or multi-trace (the Letter never said "heavy single traces"; the Supplement's (Pos) now says so too). "At leading order in 1/c_T" in abstract, Why it matters, Result, Discussion. Supplement: Table S2 rows A2, R1 (exact correlator at every c_T; planar intercept (KKSDZ 4.80) not used), H_K (at every finite c_T; order of limits [D given A2]); new paragraphs "Order of limits [D given A2]" and "Term-wise use (not made)" (incl. k = 5 and the circularity of the CHMRS (4.45) tail); I2 and the full statement point to them | L:113–117, 141–150, 196–201; S:175–176, 227–228, 236, 276–277, 319, 558–580 |
| **M1** large-spin coverage and premises | Letter now names the three near-gap regions: Theorem W on 40 sampled states with 50 ≤ J ≤ 344 (m ≤ 1.3Λ); 345 ≤ J < 857 (m ≈ 1.3–1.65Λ) and m ≥ 1.3Λ at b_A ≤ b < 1.04Λ in neither theorem (H_R′ [D-sketch + M]); Theorem L for J ≥ 0.52m up to 10⁶. Premises of W listed: leading-Regge floor [NE] there, the three constants (3, 10, 10), the flat ν⁶ subtraction [I], two structural premises [D-sketch] (P1 after the (S) failure, P2′), beyond all orders [M]. "≈ 2–3× scalar" attributed to one state (1100, 100). Theorem L: "scalar-shaped O(1/m²) correction up to 43×", floor "[P + DV]" with grid-evaluated constants. Intro: "on parts of the large-spin region it is reduced to conditional statements whose premises are measured on one state or modelled". Supplement: P0 regraded [NE] in W's region (m22); P1 [D-sketch] with the KKSD 1809.05111 finding; P4 presupposes the ν⁶ assignment [I] and the Regge-costed alternative (K₄ ≈ 35–3·10³, all J fail) stated; W's grade sentence (40 states, odd 303–343 extrapolated); Λ_amp from one state, (364, 64) inconclusive, 2σ leverage 0.5–0.6 in one fit; model K₄ ≤ 3.8 and model J\* = 38; Theorem L's tolerance is for the scalar-shaped vector, shift-by-two without light-cone normalisation only K ≈ 0.4 at J = 10⁶; Table S11 split (345–856 row, J > 10⁶ row) | L:76–81, 289–303; S:1041–1062, 1097–1108, 1119–1123, 1155–1172 |
| **M3** thin [NE±] states | J = 29: "+0.972 (Monte Carlo error 0.018; error budget 0.79; margin/budget 1.2 using the sign of the measured window, 0.94 without)"; "window measured only for ν ≳ 10–20, and for J ≳ 32 it is again a live input". Table S5: "± 0.018 (Monte Carlo error only)", grade "J = 29 only with the signed window, borderline [NE±]/[NE + M] on the sign-free ×3 rule" | L:281–285; S:811 |
| **M4(a)** "proved in the flat limit" | Removed from the abstract. Intro: "certified in the flat limit for the functional of item (1)". Status: Ω^Q_up ≥ 0.18 certified except two asymptotic regions ([NE]); Functionals paragraph keeps "negative on scalars lighter than 0.4602Λ [P] … origin of the condition (2a)" | L:32–41, 79–80, 275–278, 218–220 |
| **M4(b)** tree level vs finite c_T | "at leading order in 1/c_T" (abstract, Why it matters, Result, Discussion); "We do not bound the 1/c_T corrections, so at finite c_T this does not yet constrain t₄ at the 10⁻⁸ level". Supplement order-of-limits paragraph: remainder size Δ_g^p/c_T² not computed; Not-claimed list adds the size of 1/c_T corrections | L:91–99, 225–226, 333–337; S:565–568, 1293–1294 |
| **M4(c)** "10⁻⁸" mostly Λ⁻⁴ scaling | Why it matters rewritten: collider range τ ≤ 4 stated, "most of the factor ∼ 10⁸ between the two is the Λ⁻⁴ scaling already implied by causality [CEMZ, AKT]. What is new is the coefficient: in units of the gap \|ĝ₃\|Λ⁴ ≤ 120". `\CollFrac` no longer printed (macro kept, comment updated) | L:90–103; macros.tex:79–80 |
| **M4(d)** de Sitter sentence | **Kept, rewritten as an open question in one clause**: "Whether the flat-space core of the argument has a de Sitter analogue is open; it would need a positivity principle there, which is not established." Justification: the v5 sentence ("suggests analogous constraints … during inflation") implied a carry-over the paper gives no argument for (the AdS mechanism is a de-shift against AdS curvature); the new sentence claims nothing beyond the absence of a dS positivity principle, so it needs no dS literature and cannot be read as a result. Dropping it entirely would also be fine (it costs ≈ 30 words) | L:353–355 |
| **M5** importance (one gap, loose constants, c₁) | Not addressable by writing. The Discussion's next steps now include "positivity over a range of gaps" and the near-gap band 345 ≤ J < 857 | L:350–353 |
| **M6** unreviewed v5 layer | Supplement provenance: order-of-limits argument and premise updates are added and unreviewed | S:83–85 |
| **M7** readability | v5 "Status of positivity" (566 words) + "Large spin" (264 words) replaced by one paragraph of ≈ 310 prose words (≈ 360 counting each inline formula as a word, incl. the `\HLJstatus` sentence, which was shortened in macros.tex): what is certified, error-budgeted (key numbers only), conditional, model, and the uncovered band. All leverage values, spin lists, ring classes, bias corrections, the monopole, the J = 29 artefact story, the J = 31/32 shifts and the Ω^Q channel numbers are only in the Supplement (Tables S5–S8, S11, Sec. S6–S7), where they already were; the gauge disclosure is kept as one clause. **New Fig. 2** (below) | L:271–316; macros.tex:47–50 |

## 2. Minor issues

| item | action | where |
|---|---|---|
| m1 length | recounted (§4) | — |
| m2 B_eff(Λ) | Eq. (3) now has B_eff(Λ): "−19.265 at Δ_g = 1000 and −19.220 at 300 (it tends to the flat −19.278)". Recomputed from the published weights (make_fig.py integrals): 2L_AdS = 92.81351 / 72.04059; B_eff = 2L_AdS − 17.2407 ln(2Λ/3) = −19.26464 / −19.21966 (with A_Q = 17.240713: −19.26473 / −19.21973, the Table S4 values); the v5 single value gives 71.995 at 300. Supplement full statement says "B_eff(Λ) … depends on Λ" | L:243–248; S:287 |
| m3 overlapping spin lists | Disjoint lists (Supplement Table S11 and Table S7): [NE + M] J = 7, 9–11, 13–15, 17–20, 22–25, 27, 28, 31, 32, 36, 40, 48; interpolated [M] even 26, 30, 34, 38, 42, 44, 46; odd 33–49 H_LJ(1000). The Letter no longer lists spins beyond the [NE±] set | S:854, 1163 |
| m4 item (ii) mass range | "near-gap states with 50 ≤ J ≤ 344 (m ≤ 1.3Λ)" | L:292–293 |
| m5 "proved … floating-point constants" | "proved with grid-evaluated constants [P + DV]" | L:299–300 |
| m6 K grade | "K = 1 + O(Λ⁻²) (the shift is [DV] up to a ±0.02 % band [D-sketch + M])" | L:204–205 |
| m7 "So" + [D] | "An exact identity ([P] to eighth order …, [D] at all orders) makes every …" | L:289–291 |
| m8 Theorem L "any amplitude" | "a scalar-shaped O(1/m²) correction"; Supplement adds the light-cone-normalisation caveat | L:300–301; S:1119–1123 |
| m9 "only at J ≤ 100" | "fitted at one state, (Δ, J) = (1100, 100)" | L:297–298 |
| m10 J\* constants | Letter states the three constants (3, 10, 10) and the open subtraction; J\* itself is kept out of the Letter (only the range 50–344) | L:292–296 |
| m11 c₁ thresholds | "\|c₁\| ≳ 60 (item 1), 1306 (2a) or 1262 (2b)" | L:257 |
| m12 x formula | "x = 19941\|c₁\|/(Λ²√(16KL_A))" | L:234 |
| m13 composite grades | Intro legend: "[X + Y]", "[D given X]", weakest-link rule, pointer to Sec. S1; [NE] added to the list | L:81–89 |
| m14 Fig. 1 | annotation carries (c₁ = 0, K = 1); 300 markers are open diamonds; key: "Δ_g = 1000: claimed / Δ_g = 300: one hypothesis open / lines: implication only"; band label "nothing claimed"; caption updated | fig/make_fig.py; L:173–176 |
| m15 two scope points | H_bulk: "(one model input, tolerated ×76)"; the 1.1–1.3Λ band at J ≤ 48: "masses 1.1–1.3Λ there … [NE + M]/[M]"; Table S7 H_cont row says it is the only cover of that band | L:282, 287–288; S:861 |
| m16 DKKPS | cited for J = 1, 3; "odd J ≥ 5 couple only through a parity-odd structure" | L:118–120 |
| m17 weighted (2a) | "A scalar below 0.4602Λ with φW² coupling g₀ does not invalidate (2a); it adds (0.4304 − 2.0321 m²/Λ²)\|g₀\|²Λ⁴/8πG to the right side" (Supplement Sec. S5.F formula) | L:248–250 |
| m18 ABJM orders | "through the orders computed (supergravity and R⁴) only the stress-tensor multiplet is exchanged" | L:323–325 |
| m19 HKT | "Averaged null energy is positive [HM, HKT], and the light-ray operators that measure it commute at spacelike separation [KKSDZ, KKSDZlightray]" | L:56–59 |
| m20 residual | "≤ 7.9·10⁻⁴⁸ for b = 0…4; 2.0·10⁻⁴⁷ at b = 7 and 1.1·10⁻⁴⁷ at b = 10" (`pw1_magic.txt`) | S:985–986 |
| m21 end-to-end accuracy | "with refined quadrature (… default up to 7.6·10⁻³ at ν = 6)" | S:990–991 |
| m22 P0 grade | P0 [NE] in Theorem W's region (curvature bounds need m ≥ 1.3Λ, b ≤ 119, or b ≥ 1.04Λ) | S:1041–1043 |
| κ (§6) | "flat value 0.1099, corrections O(Λ⁻²)" | L:239 |

## 3. Fig. 2 (new, single column)

`fig/make_fig2.py` → `fig/fig2_coverage.pdf` (+ `fig2_coverage_preview.png`). Ω^EE at Δ_g = 1000 in the plane of spin J (log axis
at J + ½, so J = 0 is shown) and m/Λ: a top strip for m → ∞ (flat limit, certified), an upper log panel 1.3 ≤ m/Λ ≤ 20 and a lower
linear panel 1.0 ≤ m/Λ ≤ 1.3 (near gap). Ordinal one-hue blue ramp (dataviz reference sequential steps 700/550/400/250, validated with
`validate_palette.js --ordinal`: all checks pass), plus neutral grey with hatching for "hypothesis only"; direct labels on every region.
Regions: [NE±] bars at J = 0, 2 (to 1.2Λ), 4, 5, 6, 8, 12, 16 (to 1.05Λ), 21, 29 (to 1.10Λ) and the strip m ≥ 1.3Λ, b < 5; conditional
theorems: W (50 ≤ J ≤ 344, m ≤ 1.3Λ) and L (J + ½ ≥ 0.52105 m, J ≤ 10⁶); [NE + M]/[M]: J ≤ 32 near gap (and to 5Λ), computed even
J = 36, 40, 48, and 5 ≤ b < b_A(m); hypothesis only: odd 33–49 and interpolated even 34–46 near the gap (inside H_LJ(1000)),
m ≥ 1.3Λ with b_A ≤ b < 1.04Λ (bracket "345 ≤ J < 857: no theorem"); twist edge drawn; star at (1100, 100) "premise measured here
only". Sources per region are in the script's docstring. Choices to note: (i) the [NE±] mass ranges for J = 4, 5, 6, 8 are drawn to
1.05Λ (NEARGAP2 §5 item 4: "only masses up to 1.05Λ were budgeted"); (ii) the b < 5 strip is drawn to 20Λ, where it was computed
(POS_ERRORS §2.4; above 20Λ only J ≤ √Δ is engine-covered); (iii) J > 2·10⁵ near the gap is above the plotted masses, so the J > 10⁶
"hypothesis only" region is stated in the caption ("Theorem L up to J = 10⁶, off scale") rather than drawn.

## 4. Length and build

* **Length** (`scratchpad/tools/wc_prl6.py letter.tex 3 0.9375 1.0385`, wc_prl5 extended to several figures): prose 2501 words,
  314 inline formulas (conservative 423), 3 display rows (48, checked in the PDF), Fig. 1 180, Fig. 2 164 (aspect 1.0385):
  **≈ 3207 APS words (≈ 3316 conservative)**, limit 3750 (v5: 3110 / 3210 with one figure). Abstract **563** rendered characters.
  Text ends ≈ 0.55 down the right column of page 4 (≈ 3.8 pages). Letter 5 pages, Supplement 21 pages.
* **Build** (NOTES.md §1; into git-ignored `build/`): both compile with 0 errors, no undefined references or citations, no overfull
  boxes, no "float too large"/stuck floats (Table S2 is now `\scriptsize`; it overflowed by 111 pt after the A2/R1/H_K rows grew).
  Notices: 5 underfull boxes (Letter), BibTeX `jnrlst` notice and "empty year in SM" (intended), hyperref "token not allowed in PDF
  string" for math in Supplement section titles (cosmetic).

## 5. What remains (not addressable by writing)

1. **Running: `P2_CONTRACT`** (three-channel EE contraction of the spinning O(1/m²) symbols). It fixes Λ_amp, K_sh and the monopole
   F(i/2) for ⟨TTTT⟩, and with the ν⁶ subtraction would settle K₄. Update every `% P2CONTRACT` place: letter.tex:296, 300 (Theorem W
   constants (3, 10, 10), range 50–344; Theorem L tolerance 43), supplement.tex:1090 (Table S9 comment), 1102 (model constants),
   1142 (Table S10 comment), 1164 and 1169 (Table S11 Theorem W/L rows). If Λ_amp exceeds the Theorem-W tolerance (≈ 23 at J = 50, 27–31 at
   J ≥ 150) the large-spin support is removed (NOTES.md §8 trigger).
2. The explicit hybrid subtraction of the flat ν⁶ tower at O(1/m⁴) (P4) [I].
3. Error budgets for near-gap J = 7, 9–11, 13–15, 17–20, 22–28, 30–48 and the band 345 ≤ J < 857; J ≥ 150 window data.
4. c₁ (or F₂); a re-optimised Q with light-scalar positivity; positivity over a range of gaps (M5).
5. The size of the light-side remainder O(Δ_g^p/c_T²) (REGGE_TREE §6): needed before any statement at finite c_T.
6. Reviews: REGGE_TREE, PW_IMAGES, LARGESPIN_W/L, SMALLNU_WINDOW, DIRCV, P2_SPINNING §7 are unreviewed; interval constants for Lemma
   LR∞; an independent check of the Theorem W/L costs at 2–3 spins.

## 6. Inconsistencies found (also NOTES.md §7, v6)

1. Letter v5 B_eff single value (−19.265 at both gaps) gives 71.995 at 300 vs the quoted 72.04 (exact 72.0406).
2. Letter v5 lists put J = 26, 36, 40 in both "computed in models" (22–28, 36, 40) and "interpolated" (26, 34–46); Table S7 likewise.
   POS_ERRORS §3.2/HIGHSPIN: 36, 40, 48 have chain rows; 26, 30, 34, 38, 42, 44, 46 are interpolated.
3. Supplement P0 grade "[P + DV] for b ≤ 119" applied to Theorem W's region, where m ≤ 1.3Λ and the bound does not apply; P0 is [NE].
4. "≤ 8·10⁻⁴⁸ (also b = 7, 10)": the b = 7, 10 residuals are 2.0·10⁻⁴⁷ and 1.1·10⁻⁴⁷.
5. Letter v5 "|c₁| ≳ 1260 (item 2)" is the (2b) threshold; (2a) is 1306.
6. Supplement v5 "P1 [D given (S)]" was overtaken by P2_SPINNING §7.2 ((S) not supported by 1809.05111).
7. REGGE_TREE vs referee: the term-wise corner failure also occurs at k = 5 (and for generic k = 4 weights); the literature attribution
   "CHMRS/CLS-D state J₀ < 2 at tree level" is not borne out (REGGE_TREE §2) — the Letter cites CHMRS/CLS-D only for using the sum
   rule at finite c_T.
8. Minor: v5 Letter said near-gap budgets reach "1.05–1.1Λ"; J = 0, 2 were budgeted to 1.2Λ (POS_ERRORS §3.1) — the Letter now says
   "1.05–1.2Λ".
