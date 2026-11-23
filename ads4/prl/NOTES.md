# PRL draft: notes, build, length, provenance, referee response, open items

*Revised 23 Nov 2026 (writing pass, no new computation). Unrefereed. Nothing committed. Source of truth: `ads4/UNIFIED_THEOREM_v2.md`
(UT2), which consolidates `UNIFIED_THEOREM.md` and all fourteen of its Update boxes; UT2 §12 traces every number to a source file and
section, and UT2 §11 lists the inconsistencies found among the sources.*

Files in `ads4/prl/`:

| file | content |
|---|---|
| `letter.tex` | the Letter (REVTeX 4.2, `aps,prl,twocolumn`) |
| `supplement.tex` | Supplemental Material (REVTeX 4.2, one column): definitions, data, input table, derivation sketches, positivity evidence, corollary, reproduction, ABJM |
| `macros.tex` | shared macros, including the H_LJ(300) switch |
| `refs_verified.bib` | bibliography used by both documents (`\bibliography{refs_verified}`); `refs.bib` is no longer used |
| `data/` | published functional data (README, MANIFEST); referenced by Supplement Table S1 |
| `build/` | all build outputs |

---

## 1. Build

```
cd ads4/prl
ulimit -v 3000000; export OMP_NUM_THREADS=1
for f in letter supplement; do
  pdflatex -interaction=nonstopmode -output-directory=build $f.tex
  (cd build && BIBINPUTS=..:. bibtex $f)
  pdflatex -interaction=nonstopmode -output-directory=build $f.tex
  pdflatex -interaction=nonstopmode -output-directory=build $f.tex
done
```

Status (23 Nov): both compile with **no errors, no undefined references or citations, no overfull boxes, no hyperref
PDF-string warnings**. `build/letter.pdf` 5 pages; `build/supplement.pdf` 10 pages. BibTeX: only "missing number" warnings (issue
numbers absent in `refs_verified.bib`) and one `jnrlst` style notice.

## 2. Length

| item | value | method |
|---|---|---|
| abstract | **≈ 518 characters** rendered (pdftotext, incl. spaces; pdftotext adds spaces inside ⟨T T T T⟩, so the true count is ≈ 510) | limit 600 |
| main text, PRL counting | **≈ 2790**: prose 2427 words (inline math = 1 word; the `\HLJstatus` macro adds ≈ 60), 9 display lines × 16 = 144, Table I (single-column, 8 rows × 13 + caption + header) = 158 | `scratchpad/tools/wc_prl.py` |
| main text, rendered | ≈ 3200 words (pdftotext, Introduction → end of Discussion, incl. table and equation tokens) | upper estimate |

Both are below the 3750-word limit.

## 3. H_LJ(300) switch

`\HLJclosedfalse` in `macros.tex` (open; headline Δ_g = 1000). HLJ_CLOSE.md is now final and did **not** close the Regge window of
H_LJ(300) for J ≥ 34 (it stays [I]); it moved the (336, 28/32) windows from [M] to [NE + M] (+0.759 / +0.837), which the documents now
use. If the window is ever closed: flip the switch, set `\HLJgrade` and `\HLJsource`, and re-check the corollary paragraph (letter,
"Dependence on the gap"; Supplement §S6). The closed position prints a neutral placeholder instead of a TODO.

## 4. Provenance of the numbers

All numbers are in UT2 §12 with file and section. Letter-specific items:

| letter item | source |
|---|---|
| abstract 2.89·10⁻⁸, 92.8 | UT2 §3 (RD-S §3; DG §1.2; REPRO §1) |
| R₀ 3.13·10⁻⁶ / 2.89·10⁻⁸; x/|c₁| 2.16·10⁻³ / 1.67·10⁻⁴; 19941 | C1 §4.1 |
| B_eff −19.220 / −19.265; 72.04 / 92.81 | DG §1.2 |
| 5.946·10⁴ = 57600 I_Q² (I_Q = 1.0160366); 6.8·10⁻⁵ / 5.4·10⁻⁷ at κ = 0.1099 | C1 §1.1, §4.2 (541017 = 57600 I_Q²/κ) |
| 2.032 = 2 I_Q | C1 §1.1; DG §1.3 |
| |c₁| = 100: 3.88·10⁻⁶ / 2.94·10⁻⁸, 72.72 / 92.82 | C1 §4 tables |
| H17 ±0.24 % / ±0.021 %; δk(0) = −68.08/Λ² | DG §0; REPRO D2 |
| ε_W 1.2·10⁻¹³ / 4.2·10⁻¹¹ / 2.4·10⁻⁷ | DESHIFT §2.3 |
| cutoff effect < 10⁻⁴⁰ | HK §0 (e^{−180} Q, 10⁻⁴¹ D) |
| H_K margins ½ / 1 (J₀ < 3/2, J₀ < 2) | LRC §0, §6.1 |
| R1′ `|dDisc 𝒢| ≤ 2B(c)` | LEMMAS §2 |
| light-side decomposition; checks | LIGHTSIDE §0–§1 |
| −𝒜₄[R⁴] = 2g₄ | F_NORM §0, §2.3 |
| Table I rows | RD-S §5; HLJ §4; RQ-S §3.1; NOFORWARD §0, §4 (certified 3/20) |
| +0.2038 / +0.2009; J = 24 canonical-gauge failure; in-sample J = 18, 20, 24, 28 | RD-S §2, §5 |
| (1100, 100) violation | RD-S §4; HLJ §5; REFEREE M9 |
| s_d +0.5096 / +0.5077; s_fail 1.776; ≥ 19× | HRS §0; REPRO §2; RQ-S §3.1 (s_fail ≥ 11; 11/0.59 = 18.6) |
| far region ≲ e^{−26} vs e^{−8.5} | RD §3.2 (e^{−66} at r ≥ 0.3, volume ≲ e^{40}); UT RD box |
| s_d direct +0.50967 / +0.50776 | REPRO_ENGINE §0 |
| Λ-forms 1.3194·10⁴, 8.7415, 6.384; 17.2407, −19.2782, 10.314, 860.3, −2508.6 | DG §0, §1 |
| corollary hypotheses | DG §2.5, §3; UT DG box |
| ABJM 2.6 %, 80 % of 3.000 | UT §7; LIT §3 (2.404/3.000 = 0.801) |
| |ĝ₃|Λ⁴ ≤ 104 / 120 | √(16L_A) = 103.6 / 119.7 (RD-S §3; my arithmetic) |
| 3022.34 ln(2Λ/3) − 5327.21; ×10 weaker than CHLPS | RD-S §0 item 4, §6.4 |

## 5. Referee report: what was done

| point | status |
|---|---|
| M1 "theorem"/[D] too strong | **addressed**: "Theorem (conditional)" with every input listed precisely (I1–I6; Supplement Table S2) and net grades per item: (1) [D-sketch]; (2) [I] in one technical link, [A] via H_EFT(b), [D-sketch] elsewhere; positivity not established. Three layers separated. |
| M2 hypotheses are labels; table mixes inputs/evidence | **addressed**: precise statements in the main text; Supplement Table S2 (inputs) and Tables S5–S7 (evidence) are separate; letter Table I is evidence only; positivity defined (M_Ω(𝒪) ⪰ 0, normalisation); corollary lists all its assumptions. |
| M3 k = 4 light side, meaning of g̃₄ | **addressed** by LIGHTSIDE_DIST + F_NORM: one functional, Ω₄ = D_LR + 𝒜₄, same object on heavy blocks; −𝒜₄[R⁴] = 2g₄ by calibration [D]; g̃₄ = g₄ + O(Λ⁻²) given H_EFT(b) (CHMRS/CLS-D standard); H7's value is the arc at n = k = 2 (item 1 sentence). **Remaining:** a = b integrability of the R⁴ arc [I]; spinning H_EFT(a). |
| M4 lower end | **addressed**: "−O(Λ⁻²), no computed constant" stated in abstract, theorem, discussion; graded with item (2). |
| M5 c₁ | **addressed**: one universal constant, u-channel and k = 3 parts vanish; exact bounds for all c₁ (both items); "uncomputed", not "free". **Remaining:** value of c₁. |
| M6 de-shift, weights | **addressed**: continuum weights for all t (ε-DS [DV]); C5 "[I]" was stale; upper end C5-independent; data published; Abel weights re-solved independently. **Remaining:** X^Q (Q's k = 2 Abel weight) is not in `data/`; no continuum k = 2 weight for Ω^EE; ε_W not interval-certified. |
| M7 Lemma C / Einstein instance | **addressed**: Lemma C withdrawn, C′ [D]; Lemma A [D given C′]; Einstein instance "derived, not evaluated numerically"; collinear set for J₀ < 4 terms [D]. |
| M8 H_K buried | **addressed**: H_K stated precisely (I2), derived by Tonelli from J₀ ≤ 1 with margins ½ / 1, with the analytic cutoff in the definition; the Paley–Wiener no-go for grid weights disclosed. **Remaining:** the EE step (EESTEP, in progress elsewhere). |
| M9 1000 headline relies on refinement; certificate | **addressed**: (1100, 100) violation and refined form in main text; certificate stated as 𝔉 ≥ (3/20)N except one corner (and the "0.200 [certified]" error corrected). |
| M10 H_RS | **addressed**: "[D given P_cyc], preliminary (previously [I])"; P_cyc untested by the reproductions (main text); far-region number corrected to r ≥ 0.3 with volume factor; coverage stated; contiguity argument summarised in Supplement. **Remaining:** RD_CLOSE's divisor table is not reproduced in the Supplement; HRS_CLOSE / RD_CLOSE still unreviewed. |
| M11 Δ_g = 300 inconsistencies | **addressed**: "lowest gap" sentence qualified (only with the [D-sketch] refinement); corollary stated only with all hypotheses, including positivity at 300 (open); HLJ_CLOSE reconciled (+0.759 / +0.837). |
| M12 title/scope | **addressed**: title "Conditional bounds on R³ and one R⁴ graviton coupling …"; g₄ defined by the MHV amplitude; ĝ₄ (all-plus) and g₅ stated unbounded; A1 restricted to single traces in T × T; "no known theory satisfies A1" stated. |
| M13 ABJM | **addressed**: out of abstract; short normalisation check citing BCP (`BCPW` key, 1906.07195) and CPY; 80 % of CHMRS's 10d window; 2.6 % reflects the 4d IR log; outside A1; c_T/flavour parenthesis deleted. |
| M14 novelty/credit | **addressed**: positioning paragraph (KKSDZ, BHM, CEMZ, AHKT, CHP, MP, AKT, AHKTshock, CHMRSdisp, CHMRSswamp, CaronHuotVanDuong, CHMRS, Kundu, CHLPS, CHLPSpw, BKZ, CLS-D, CPY, BCP, BCJP); TODO and search-level language removed; LITERATURE.md §4 replacements applied (rows 1–14; row 13's sentence dropped with the paragraph). |
| M15 jargon, derivations, disclosure | **addressed**: "engine", "crystal", "fold", "HYB", "chain" removed from the main text (defined once in Supplement §S5); s_d, P_cyc, H_RS, D_LR, binding, de-shift defined; Supplement has derivation sketches (arc decomposition, C′, F calibration, H_K/LRC*, c₁ exact bounds, de-shift identity, Λ-forms); chronology and gauge sensitivity disclosed in the main text. **Remaining:** derivations of Adm, H7 and L0 and the explicit positivity matrix are only cited, not sketched. |
| minors 1–10, 12–16 | addressed (δ ≈ 0.2 nominal; units/c_T → ∞; α_T citations CDPrakash, GPY; H17 weighted norm; 19×; J = 0, 2 at 300 superseded (Supplement); truncated Λ-form unsafe; +0.2054 protocol recorded; KKSDZ §4.2.4/4.2.6 + Caron-Huot; refs_verified; "uncomputed"; J = 15 ≈ +0.69; hyperref warnings gone; grade notation; 1.00499Λ). |
| minor 11 front matter | **remaining**: affiliations, emails, acknowledgments are TODO (not in the repository). |
| referee's overall verdict (long-form venue) | **not addressed** by writing: positivity is still numerical on finite sets inside a model. The letter states this plainly. |

## 6. Open items before submission

1. EESTEP: when `EESTEP_CLOSE.md` lands, update the H_K grade in letter I2 and Net grades, Supplement Table S2 (LRC*, H_K rows) and
   §S4.D ("EE step (pending)"), and UT2 §3–§4 (search for "EESTEP").
2. Publish X^Q (Q's k = 2 min-norm Abel weight) in `data/` (outside this task's authorised files).
3. Front matter (affiliations, emails, acknowledgments).
4. Optional: RD_CLOSE's divisor table and HRS_CLOSE's Lemmas 1–2 in checkable form (referee M10).
5. Before submission, the INSPIRE sweep recommended in LITERATURE.md §6.

## 7. Source inconsistencies

See UT2 §11 (20 items). The ones that changed letter text: certificate 0.15 vs 0.200; ε_{c₁} linear → exact; sign of (F); H_K k = 4
criterion 5/2 → 2; far-region number; DS for Q [NE] → [DV]; (336, 28) +0.759 as lowest large-J value at 300; "≥ 19×"; X^Q not published.
