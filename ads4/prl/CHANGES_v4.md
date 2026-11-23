# CHANGES v3 → v4 (response to `REFEREE_REPORT_2.md`)

*23 Nov 2026. Writing pass only: no physics computation was redone except the figure's curves (the same exact integrals as
`code/repro/bounds.py`, recomputed from `prl/data/`), the Q3\* export (exact division; certified fingerprint reproduced) and two
conversions (λ for Λ = ML; the |ĝ₃| ratio to the transplanted CHLPS optimum). Source of truth unchanged: `ads4/UNIFIED_THEOREM_v3.md`
(UT3). Files edited: `letter.tex` (rewritten), `supplement.tex`, `macros.tex`, `refs_verified.bib`, `NOTES.md`, `data/README.md`,
`data/MANIFEST.txt`. New: `data/q3star_factored.json`, `data/make_q3star.py`, `fig/make_fig.py`, `fig/fig1_bounds.pdf`,
`fig/fig1_bounds_preview.png`, this file. Nothing outside `ads4/prl/` was touched; nothing committed. Line numbers refer to the files
as of this writing.*

**Principle.** The Letter is narrowed to item (1) at Δ_g = 1000 (headline) and item (2a)'s upper end (secondary), with (2b) in one
sentence. Everything removed from the Letter is in the Supplement (§S1 "What moved here from the Letter", supplement.tex:52).

## 1. Major issues

| issue | what was done | where |
|---|---|---|
| **M1** positivity; main text understates weak points | New "Status of positivity" paragraph names, in plain words: (i) near-gap J ≥ 12 budget does not close [NE + M]; (ii) the J = 29 edge state's sign depends on the edge model (−3.6 vs +1.39 ± 0.03), possibly other odd J ≳ 27; (iii) H_LJ(1000) [D-sketch + M], sign-free form violated at (1100, 100) (−2154 ± 43 vs 1409; verified in HWIN_THIMBLE via RESULT_D/UT2/UT3 and HLJ_CLOSE §5–6: a leading-Regge-shaped excess; HLJ_CLOSE's engine "neither supports nor violates"), holds only in refined form; (iv) b ≥ 5 model-only, no m_*; (v) post-hoc choice of EE and gauge, canonical-gauge failure at J = 24 (300). The `\HS` markers and `% HIGHSPIN:` comments are on (i)–(iii). [NE±] is described as "error budget with ×3 safety factors on model terms; not a bound". Fallback for a negative HIGHSPIN state: NOTES.md §8. | letter.tex:231–260; NOTES.md §8 |
| (M1 fix 2) abstract | Abstract now says "positivity … (proved in the flat limit, tested numerically in AdS)" and includes light spin-0, 2 states; listing modelled regions does not fit in 600 characters, so the list is in the main text. | letter.tex:29–37 |
| **M2** lower end | Dropped from the abstract, Result and Discussion; one sentence in the Letter says it is not claimed and why. Kept in the Supplement as a conditional outlook with its grade, the decaying leading-Regge margin βe^{−β/2}, no constant, c_L slope not reproduced. | letter.tex:224–226, 294; supplement.tex:312–327 (§S4.A); not-claimed list (§S10) |
| **M3** unreviewed derivations | Not addressable by writing. Stated in NOTES.md §6 item 6 as the top open item (write-ups of H_K, L6′, calibration, c₁; independent reproduction of the EE reduction and Lemma Y). | NOTES.md §6 |
| **M4(a)** c₁ universality [D] | → [D-sketch] in the Letter, Table S2 (L3) and §S5.E. | letter.tex:157; supplement.tex:220, 554 |
| **M4(b)** arc statement drops qualifier | The Letter no longer states the arc grade separately; Table S2 L4 keeps "[D given C′ + Mellin convergence]". | supplement.tex L4 row |
| **M4(c)** L6′ | Letter: "[D] … rests on an analytic counting lemma; the direct numerical scan is still pre-asymptotic (local slopes 2.5–2.8 against the predicted 2)". Table S2 L6′ and I4 say the same, plus "AdS₄ R⁴ correlator not evaluated". | letter.tex:159–164; supplement.tex:224, 303 |
| **M4(d)** H_EFT(a) | Table S2 L7 and I5 now carry the unexplained ≈ −3 % spin-2 ratio [NE] (calibration uses only a scalar). The Letter no longer claims "[D; DV] for all dispersive terms". | supplement.tex:225, 307 |
| **M4(e)** "bound"/"error-bounded" | Replaced by "error budget"/"error-budgeted" everywhere; §S1 and §S6 say what the budget covers and that only the flat floor and the Lemma-C curvature term are bounds ("flat floor minus Lemma C alone (a bound): ≥ +0.19993 N_s"). | letter.tex:65–67, 236–238; supplement.tex:42–49, §S6 method paragraph, Tables S5, S7, S8 |
| **M4(f)** "Theorem" | Now "Result (conditional; implication derived at sketch level, [D-sketch])" in the Letter and §S4. No "Theorem" remains in either document. | letter.tex:188; supplement.tex:251 |
| M-physics: parity-odd half of τ = 240\|ĝ₃\| | Stated with its grade from DICTIONARY_T4: t₄ half [DV]; α_T half [D] with sign convention [I]; only the modulus enters. | letter.tex:103–106; supplement.tex:78–81 |
| A1′ and ABJM (m12) | A1′ now holds for ABJM-IIA by the spin count alone (supergravity spins ≤ 2; lightest J ≥ 4 is a string state). The "only the stress-tensor multiplet is exchanged" argument is presented separately, as what satisfies the (2a) scalar condition (and makes rows LS vacuous). | letter.tex:263–275; supplement.tex:933–946 |
| **M5** positioning | CHLPSpw (2205.01495) verified on the arXiv API: "Graviton partial waves and causality in higher dimensions", D ≥ 5, "record new bounds on the central charges of holographic conformal theories". The Letter now describes it so, and names CHLPS's CFT₃ remark ("without a derivation") as precedents for item (1). The comparison is "≈ 10× weaker in \|ĝ₃\| (≈ 10² in \|ĝ₃\|²)" (computed 10.1 / 10.3 at 300 / 1000; shown in Fig. 1a). The g₄ upper end is stated to be IR-log dominated, and ABJM's 2.6 % to be a one-digit consistency check (with the ≈ 3 % Λ ↔ ML shift). | letter.tex:56–62, 271–275, 282–291; supplement.tex:§S5.H (Λ-forms) |
| **M6** broad audience | Letter rewritten: physics-first introduction; a plain-language "Mechanism" paragraph (commuting energy detectors → sum rule; heavy side positive in the flat limit; light side fixed by bulk couplings; what AdS curvature spoils; the three ingredients EE, de-shift, k = 4 arc); one figure. Moved out: Inputs I1–I6 (→ §S4), c₁ threshold details (→ §S4, §S5.E), Λ-forms (→ §S5.H, displayed), Δ_g ≥ 300 corollary (→ §S7, reworded), Reproduction (→ §S8, already tabulated), Table I (its rows are in Tables S5/S7). Jargon removed from the main text: de-shift explained once, no D_LR, s_fail, "binding pair", N_s, Mehler–Dirichlet, "dual-line flip", "light-ray one-point matrices". | letter.tex:42–296 |
| (M6) figure | `fig/fig1_bounds.pdf` (single column, 3.375″ × 3.6″): (a) \|ĝ₃\|Λ⁴ ≤ √(16L_A) vs Δ_g against the transplanted CHLPS flat optimum, (b) (2a), (2b) upper ends vs Δ_g with 2ζ(3); 1000 filled, 300 open, Δ_g < 300 shaded "not claimed"; lines labelled as values of the implication. Colours: first three slots of the dataviz reference palette (all-pairs CVD-validated), each series also with its own line style and a direct label. A (t₄, α_T) disk plot was rejected: R₀ is 10⁻⁸ of the energy-flux disk radius, so both cannot share axes; the caption gives the disk in these units (Λ⁴/60). | letter.tex:174–186; fig/make_fig.py |
| **M7** (2a) condition | Abstract: "if no scalar lighter than 0.46Λ couples to TT". Result (2a): "if no scalar single trace with m < 0.4602Λ appears in T × T". (2b) states that it needs (Pos) for Ω^EE; the Letter says (2a) is independent of Ω^EE positivity and (2b) is not. Re-optimising Q is listed as a next step. | letter.tex:35–36, 203–216, 219–220, 294–296 |
| **M8** corollary, "explicit dependence on Λ" | Removed from the abstract. In the Letter: one sentence, "values of the implication … other gaps are not claimed". In the Supplement the corollary is now "A conditional remark on Δ_g ≥ 300" with a "What it assumes" paragraph: an [I] input at 300; H_LJ(Δ_g) at every gap, close to the conclusion; proxy support only. | letter.tex:226–229; supplement.tex:857–880 |

## 2. Minor issues

| item | what was done | where |
|---|---|---|
| m1 length | Recounted (NOTES.md §2): ≈ 2395 APS words (≈ 2480 counting long inline formulas as several words), 3 display rows, figure 180, no table; text ends ≈ 0.95 down page 3. Abstract 551 characters. | NOTES.md §2 |
| m2 front matter | "TODO: email" no longer prints; affiliations/acknowledgments are visible placeholders "[… to be supplied …]" with `% TODO(authors)` comments. **Authors must supply them.** | letter.tex:20–24, 298–301; supplement.tex:26–30 |
| m3 JHEP years | All 41 JHEP entries now use `volume = {year}`, `number = {issue}` (apsrev4-2 special-cases SISSA journals and prints "J. High Energy Phys. 2009 (10), 079"). "missing number" warnings gone. | refs_verified.bib (header comment, all JHEP entries) |
| m4 abstract wording | "an R⁴ coefficient defined by the CFT obeys g̃₄Λ⁶/8πG ≤ …"; light spin-0, 2 states included; scalar condition precise. | letter.tex:29–37 |
| m5 notation | R → τ (both files); k(0) → K₀ (flat value of K), k₀ → 𝒩_D; Abel a, b → r_±; Q's c → c_Q; c(ν) → ϱ(ν); footnote-style clause "Λ … not a cosmological constant". Left as is: k for the Regge spin; the cutoff width a and Q's a (flagged in place); the Regge-box parameter c in R1′. | letter.tex:96–97, 103; supplement.tex:52–57, 78–81, 106–135 |
| m6 grade legend | Letter: compact legend, composites avoided except [NE + M], [D-sketch + M]. Supplement §S1 defines the composites ([X + Y], [X; DV], [D given X], "[D, one step D-sketch]" = net [D-sketch]). | letter.tex:64–68; supplement.tex:42–50 |
| m7 Table I | Table I removed from the Letter; every row is in Tables S5/S7 (light J = 0, 2 row graded [P given LR]; [NE±]; [NE + M]; J = 8 "λ-level budget"; large-mass row carries the ×76 bulk-ν tolerance). | supplement.tex Table S5 |
| m8 "≥ 19×" | → "≥ 18×" (11/0.59 = 18.6). | letter.tex:256; supplement.tex:826 |
| m9 A1′ wording | "Spins 1 and 3 never appear: identical conserved T's in d = 3 have no ⟨TTO_J⟩ structure for J = 1, 3"; same in Table S2 A1′. | letter.tex:93–94; supplement.tex A1′ row |
| m10 parity violation | New paragraph in §S2.C: flat side [D] (MHV-only functionals, no cross terms); AdS side, parity-invariant dual structures ⇒ even in the odd OPE coefficient, mixed waves vanish ([D; DV] for the earlier functionals, carry-over to EE and the de-shifted pieces [D-sketch]); odd-coupling normalisation [I], unused. | supplement.tex:160–169 |
| m11 natural size of c₁ | Withdrawn (rested on the corrected seagull picture); thresholds kept. | supplement.tex:559–563 |
| m12 ABJM | A1′ decoupled (above); λ = 2.51·10⁷ / 3.15·10⁹ for Λ = ML (old values noted as ML = Δ_g); CPY cited for IIA only alongside BCP; CHMRSswamp Sec. 3.8, Eqs. (3.42)–(3.43) for [0, 3.000] and 2ζ(3). Bib key BCPW → BCP. | supplement.tex:933–955; letter.tex:268–270 |
| m13 references | Added, verified on arXiv API and INSPIRE: Córdova–Maldacena–Turiaci 1710.03199 (cited at the light-scalar condition), Häring–Zhiboedov 2202.08280 (Regge input R1), Hartman–Kundu–Tajdini 1610.05308 (ANEC from causality). | refs_verified.bib (end); letter.tex:52, 172; supplement.tex:209, 666 |
| m14 truncated quartic form | Letter no longer uses it; §S5.H calls it asymptotic and on the unsafe side by 0.02–0.06. | supplement.tex §S5.H |
| m15 flat certificate | "the code for one corner lemma is unreviewed" next to [certified] in the Letter and §S6. | letter.tex:234–235; supplement.tex:751 |
| m16 corollary J = 0 | Already "[D given RD] + [NE]" at J = 0 in §S7 (iii); unchanged. | supplement.tex §S7 |
| m17 "plausibly ≤ 0.80" | → "≤ 0.80 at large mass [M]". | supplement.tex Table S6 |
| m18 Regge caveat | R1 row and I1 say what KKSDZ Sec. 4.2.6's caveat is (the light-cone-limit step's loophole at the first/second-sheet boundary ρ̄ = 1) and why it lies outside the Regge box. | supplement.tex:209, 293–296 |
| m19 data hygiene | `data/README.md` §3 now quotes the certified ε_W 2.673·10⁻¹³ / 9.058·10⁻¹¹; new §5 errata table for the stale `note` in `nf02_psi.json` ("3/20 except …" → ≥ 1/5 everywhere [certified]) and the midpoint ε_W inside the C4 JSONs (files left byte-identical so the published hashes stay valid); **Q3\* published** as `data/q3star_factored.json` (exact X, Y, expanded ψ₂, ψ₃, remainders, fingerprint `80e0834a…33f2` reproduced by `data/make_q3star.py`, standard library only); MANIFEST and Table S1 updated. | data/README.md §1, §3, §5, §6; data/MANIFEST.txt; supplement.tex:172–185 |

## 3. Supplement: structural changes

* §S1: grade composites; "What moved here from the Letter"; notation changes.
* §S2: τ and its grades; Λ footnote; K₀, 𝒩_D, r_±, ϱ, c_Q; Q3\* file; parity paragraph; Table S1 (Q3\* row, errata pointer).
* §S3 Table S2: A1′ (m9), R1 caveat (m18) + Häring–Zhiboedov, L2 (K₀), L3 [D-sketch], L6′ pre-asymptotic, L7(a) −3 %; net-grade
  paragraph points to §S4 and the lower end.
* **§S4 (new) Full statement**: the v3 Letter's theorem as "Result (conditional; [D-sketch])", items (1), (2a), (2b), (2′), exact c₁
  forms and thresholds, H17 remark, net grades, dependence, inputs I1–I6 verbatim in substance; **§S4.A Outlook: the lower end**.
* §S5.E: c₁ universality [D-sketch]; natural size withdrawn; τ. §S5.G: CMT. §S5.H: the Letter's Λ-forms (Eqs. S11–S12), AKT remark,
  CHLPS comparison (121× in the log coefficient, ≈ 10× in |ĝ₃|), "values of the implication, not bounds".
* §S6: error budget ≠ bound; corner code unreviewed; "Error-budgeted" captions; m17; m8.
* §S7: corollary → conditional remark, with "What it assumes".
* §S9 ABJM: A1′ by spin count; exchange argument separately for (2a); λ; equation numbers; CPY/BCP.
* §S10: not-claimed list adds the lower end and the remark.

## 4. What remains (authors / pending computations)

1. **HIGHSPIN** (another agent): fill the `\HS` places — letter.tex:241–249 (items (i)–(iii)), supplement Tables S5 (J ≥ 12 row),
   S6 (J = 12, 16, 29), S7 (H_M, H_edge, H_LJ(1000)) — and set `\HSdraftfalse`. If any state is negative, apply the fallback in
   NOTES.md §8 (only item (2a) survives; precise statement there).
2. **Front matter**: affiliations, e-mails, acknowledgments/funding (authors).
3. Referee §6 items not addressable by writing: (2) close (Pos) for Ω^EE at one gap across the heavy spectrum (near-gap J ≥ 12 at
   [NE±], J = 29, an analytic large-spin/large-b argument); (3) full derivations and independent reproduction of the H_K EE reduction,
   Lemma Y, the R⁴ calibration and the c₁ structure; (6) the value of c₁ or premise F₂; (7) a Q functional optimised with light-scalar
   positivity (removes the 0.46Λ condition and the ε repair); Ω^Q_low evidence if a lower end is ever to be claimed.
4. The referee's PRL verdict hinges on item 2 above; NOTES.md §8 records that under the fallback the PRD/JHEP route becomes the default.
5. Housekeeping left: `refs.bib` (unused) not updated; UT3 still has "universality [D]" and the c₁ natural-size estimate (not edited:
   outside the allowed files); INSPIRE sweep before submission; C5 ε_W certification (optional).
6. If the H_LJ(300) switch is ever flipped, regenerate Fig. 1 with the marker roles swapped (the text follows the macros).
