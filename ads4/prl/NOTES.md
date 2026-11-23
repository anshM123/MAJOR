# PRL draft: notes, build, length, provenance, referee responses, fallback plan, open items

*Revised 23 Nov 2026 (v4 writing pass after `REFEREE_REPORT_2.md`: the Letter is narrowed to item (1) at Δ_g = 1000 plus item (2a)'s
upper end, rewritten for a general reader, and given one figure; everything removed from the Letter is in the Supplement). The only new
arithmetic in v4: the figure's curves (same integrals as `code/repro/bounds.py`, recomputed from `prl/data/`), λ(Λ = ML) = 2.51·10⁷ /
3.15·10⁹, the |ĝ₃| ratio to the transplanted CHLPS optimum (10.1 / 10.3), and the Q3\* data file (exact division, fingerprint
reproduced). Unrefereed. Nothing committed. Source of truth: `ads4/UNIFIED_THEOREM_v3.md` (UT3). Change logs: `CHANGES_v3.md`
(v2 → v3), **`CHANGES_v4.md` (v3 → v4, referee 2 point by point)**.*

Files in `ads4/prl/`:

| file | content |
|---|---|
| `letter.tex` | the Letter (REVTeX 4.2, `aps,prl,twocolumn,floatfix`), v4 |
| `supplement.tex` | Supplemental Material (one column): grades and composites, definitions/data, inputs table, **full statement + I1–I6 + lower-end outlook (new §S4)**, derivation sketches, Λ-forms (now with the displayed forms), positivity evidence, conditional remark on Δ_g ≥ 300, reproduction, ABJM |
| `macros.tex` | shared macros: H_LJ(300) switch, `\QuartRep`, `\GthreeHead` (new), HIGHSPIN marker `\HS` |
| `fig/make_fig.py`, `fig/fig1_bounds.pdf` | Fig. 1 (new): bounds vs Δ_g from the published data; `fig1_bounds_preview.png` is a raster preview |
| `refs_verified.bib` | bibliography; v4: JHEP year convention, key BCPW → BCP, CMT, HaringZhiboedov, HKT, `SM` |
| `refs.bib` | legacy, unused (not updated in v4) |
| `data/` | published functional data; v4 adds `q3star_factored.json` + `make_q3star.py`, README errata §5, MANIFEST lines |
| `build/` | all build outputs (git-ignored) |

---|---|
| `letter.tex` | the Letter (REVTeX 4.2, `aps,prl,twocolumn`) |
| `supplement.tex` | Supplemental Material (REVTeX 4.2, one column): definitions, data, input table, derivation sketches, light single traces, positivity evidence, corollary, reproduction, ABJM |
| `macros.tex` | shared macros: the H_LJ(300) switch, the repair headline `\QuartRep`, the HIGHSPIN marker `\HS` (switch `\HSdrafttrue`) |
| `refs_verified.bib` | bibliography used by both documents; v3 adds Paulos, CGPspin, PenedonesMellin, JJMNS (verified against arXiv and INSPIRE) |
| `refs.bib` | legacy, unused; the four v3 entries are mirrored for completeness |
| `data/` | published functional data (README, MANIFEST; `k2/` with `X^Q` and the continuum k = 2 weights, SHA256SUMS) |
| `build/` | all build outputs (git-ignored) |
| `CHANGES_v3.md` | change log v2 → v3 and the HIGHSPIN placeholder list |

---

## 1. Build

```
cd ads4/prl
ulimit -v 3000000; export OMP_NUM_THREADS=1
(cd fig && python3 make_fig.py)          # Fig. 1 (prints the headline constants as a check)
(cd data && python3 make_q3star.py)      # optional: regenerates q3star_factored.json, asserts the certified fingerprint
for f in letter supplement; do
  pdflatex -interaction=nonstopmode -output-directory=build $f.tex
  (cd build && BIBINPUTS=..:. bibtex $f)
  pdflatex -interaction=nonstopmode -output-directory=build $f.tex
  pdflatex -interaction=nonstopmode -output-directory=build $f.tex
done
```

Status (23 Nov, v4): both compile with **no errors, no undefined references or citations, no overfull boxes, no stuck floats**.
`build/letter.pdf` **4 pages** (text ends ≈ 0.95 down page 3; references fill page 4); `build/supplement.pdf` 15 pages. BibTeX: one
`jnrlst` style notice and "empty year in SM" (the Supplemental-Material pointer, intended). The JHEP "missing number" warnings are gone.
Figure: `cd fig && python3 make_fig.py` (needs mpmath, matplotlib) before compiling the Letter.

## 2. Length (v4)

| item | value | method |
|---|---|---|
| abstract | **551 characters** rendered incl. spaces (545 after removing the spaces pdftotext inserts inside ⟨TTTT⟩) | pdftotext → whitespace-normalised `len()`; limit 600 |
| main text, APS counting | **≈ 2395** (inline formula = 1 word) / **≈ 2480** (conservative: long inline formulas counted as ⌈chars/6⌉ words) | prose 1928 words from `\maketitle` to the acknowledgments, captions included; 239 inline formulas; 3 display rows × 16 = 48; Fig. 1 single-column, aspect 3.375/3.6 = 0.94 → 150/0.94 + 20 = 180; no table |
| tool | `scratchpad/tools/wc_prl4.py letter.tex 3 0.9375` | display rows checked in the PDF (Eqs. (1)–(3) one row each) |
| rendered | text ends ≈ 0.95 down page 3 of the REVTeX `prl,twocolumn` layout, i.e. ≈ 2.95 pages incl. the title block | page images |

Target was ≤ 3300 APS words and ≤ 3.7 pages: met with ≈ 900 words of headroom for the HIGHSPIN text (v3: ≈ 3080 by the old tool, which
undercounted, ≈ 4.2 pages).

## 3. Switches

* **H_LJ(300)**: `\HLJclosedfalse` (open; headline Δ_g = 1000). Unchanged. The abstract, the introduction's headline and the
  discussion use `\CubicHead`, `\DgHead`, `\QuartHead`, `\GthreeHead`, so they follow the switch; the Result paragraph prints both gaps.
  Fig. 1 marks 1000 filled, 300 open, and does not follow the switch (re-run `make_fig.py` with the markers swapped if it is flipped).
* **HIGHSPIN**: `\HSdrafttrue`. Places in v4 (each with a `% HIGHSPIN:` comment): `letter.tex` Status-of-positivity items (i) near-gap
  J ≥ 12, (ii) J = 29 edge, (iii) H_LJ(1000); `supplement.tex` Tables S5 (J ≥ 12 row), S6 (J = 12, 16, 29), S7 (H_M, H_edge, H_LJ(1000)).
  The v3 Letter's Table I placeholders are gone with the table (its content is in Tables S5/S7). When HIGHSPIN lands: update those
  places, UT3 §5.2/§5.6/§10, and set `\HSdraftfalse`. If any state comes back negative, use §8 below.

## 4. Provenance of the numbers (letter items new or changed in v3)

| letter item | source |
|---|---|
| A1′ (higher-spin gap), spins 1, 3 absent from T × T | A1PRIME §2.1, §6 |
| light Regge spins 3 / 4; drop-out at k = 4, 5; block action at k = 2, 3 | A1PRIME §2–§3, §5 |
| Ω^EE positive on light states at every mass; Ω^Q negative on scalars below 0.46Λ | A1PRIME §4.2–§4.3 (exact polynomials) |
| repair 0.04 Ω^EE; 125.71 / 164.51 | NEARGAP2 §0 ("Repaired Q"); UT3 §3 |
| c₁² coefficients 1.29·10⁻⁴ / 1.03·10⁻⁶ (repaired) | this pass: `57600(I_Q+εI_D)²/((κ+ε/8)Λ⁴)` (A1PRIME §4.3 formula; reproduces A1PRIME's ε* values 1.2598·10⁻⁴ / 1.0062·10⁻⁶) |
| 1 % tolerances 4.6 / 60 (item 1), 99 / 1262 (item 2, either form) | C1_VALUE §6 (published form: 103 / 1306); this pass (repaired form) |
| c₁ universal; shock part identically 0; energy-channel T exchange; zero given F₂ | C1_VALUE §0, §2–§4 |
| ε_W ≤ 2.7·10⁻¹³ / 9.1·10⁻¹¹ [certified]; 2.4·10⁻⁷; k = 2 ≤ 1.1·10⁻⁹ relative | HK_RESIDUAL §3.3; DESHIFT §2.3; CERT_CLOSE §3.2 |
| R⁴ light function Δφ⁻² log Δφ, criterion Δφ^{−7/2} | R4ARC_HEFT §0, §2 |
| H_EFT(a) pointwise for dispersive terms [D; DV] | R4ARC_HEFT §3 |
| g₅ not cancelled; s_d mass-independent Λ → 20Λ | EESTEP_CLOSE §5; HK_RESIDUAL §4 |
| H_K [D]; EE kernels = finite sums of shifted scalar partial waves; J₀ < 3/2, 5/2; margins ½, 3/2 | EESTEP_CLOSE §0–§3, §6; HK_RESIDUAL §0–§1 |
| I1: DV on free-field ⟨TTTT⟩ | HK_RESIDUAL §2 |
| flat ≥ 1/5 (Ω^EE) [certified], attained to 10⁻⁸; ≥ 0.18 (Ω^Q) on most states, min 0.1835 | CERT_CLOSE §0–§2 |
| binding ≥ +0.198, ×78; ≥ +0.146, ×3.5 at 300; J = 0, 2, 4, 6, 8 [NE±] | NEARGAP2 §0; POS_ERRORS §3.1 (J = 0, 2) |
| large mass m ≥ 1.3Λ, b < 5: ≥ 8.7×; no m_* | POS_ERRORS §2.4 |
| Mehler–Dirichlet curvature bound [P] | POS_ERRORS §2.1 (Lemma C) |
| light J = 0, 2 positive to Δ_g − 1 (leading Regge) | NEARGAP2 §0 |
| H_RS: physical region = D′; O(Δ^N e^{−0.336Δ}); scalar 10⁻⁸ at Δ = 20 | PCYC_CLOSE §0–§1; RD_CLOSE Prop. D′ |
| ABJM inside A1′, A2; 2.6 % / 3.3 %; one digit | A1PRIME §7 |

**v4 additions.**

| letter item | source |
|---|---|
| Fig. 1 curves; τ₀ = 2.89126·10⁻⁸ / 3.13206·10⁻⁶; 2L_AdS = 92.8135 / 72.0406; (2b) 164.51 / 125.71 | `fig/make_fig.py` from `data/nf02_psi.json`, `data/q3star_factored.json` (prints them); = UT3 §3, REPRO |
| \|ĝ₃\|Λ⁴ ≤ 119.75 / 103.61; transplanted CHLPS optimum √(24.9 ln(2Λ/3) − 27.6) = 11.59 / 10.21; ratio 10.3 / 10.1 | `make_fig.py`; RESULT_D.md (route A′ comparison, "≈ 10× weaker", log coefficient 121×) |
| energy-flux disk ⇔ \|ĝ₃\|Λ⁴ ≤ Λ⁴/60 | τ ≤ 4, τ = 240\|ĝ₃\| (DICTIONARY_T4 §0) |
| τ = 240\|ĝ₃\|: t₄ half [DV], α_T half [D] with sign convention [I] | DICTIONARY_T4.md:27–29 |
| (1100, 100): −2154 ± 43 vs 1409, leading-Regge (D_h) shaped excess; refined form | HWIN_THIMBLE (via RESULT_D.md update; fw −2123 ± 172); HLJ_CLOSE §5–6; UT3 §5.2 |
| J = 29: −3.62 (plain chain, one K₂ outlier at x = 0.8) / +1.39 ± 0.03 (direct nodes); other odd J ≳ 27 possibly alike | NEARGAP2 §0, §1.6 (lines 36, 42, 230–235) |
| J = 12, 16 minima +0.6759, +0.7191 | supplement Table S6 (NG2) |
| Q: +0.143, s_d/s_fail = 0.51/1.776 (factor 3.5); ≥ 18× (11/0.59 = 18.6); [NE±] on 1.3–3Λ at 1000 | UT3 §5.4; PE §2.5 |
| c₁ 1 % thresholds at 1000: 59.6 (1), 1306 (2a), 1262 (2b) | UT3 §3 |
| ABJM: (Λ/ML)⁶ shift ≈ −3 % at 1000; λ = 2.51·10⁷ / 3.15·10⁹ for Λ = ML | UT3 §8; `(Λ/2^{5/4}π^{1/2})⁴` (referee 2 m12) |
| CHLPSpw = D ≥ 5 graviton partial waves; "record new bounds on the central charges of holographic conformal theories" | arXiv API abstract of 2205.01495 (checked 23 Nov 2026) |

All other letter numbers: unchanged from v2/v3 (UT2 §12, UT3 §12).

## 5. Referee reports

**Second report (`REFEREE_REPORT_2.md`, reject for PRL as written).** Point-by-point response with file:line: `CHANGES_v4.md`. Summary:
§6 items 1 (narrow/rewrite + figure), 4 (positioning), 5 (lower end dropped from abstract/discussion), 8 (housekeeping) done; M4 grade
slips fixed; M1 weak points now named in the main text; fallback written (§8 below). Not addressable by writing: §6 items 2 (close
positivity at one Δ_g), 3 (auditable derivations, independent reproduction of the EE reduction and Lemma Y), 6 (c₁), 7 (re-optimised Q).

**First report: status after v3** (kept for the record).

| point | status |
|---|---|
| M1 theorem/[D] too strong | addressed (v2); v3 net grades: (1) [D-sketch]; (2) [D-sketch] + [A] via H_EFT(b); **no [I] input**; positivity not established. |
| M2 hypotheses as labels | addressed (v2); v3 adds the [NE±] grade with explicit leverage, and Table S5 separates error-bounded from model-only states. |
| M3 k = 4 light side | addressed; **remaining items closed in v3**: a = b integrability (L6′ [D]) and spinning H_EFT(a) ([D] for dispersive terms). |
| M4 lower end | addressed (no constant). |
| M5 c₁ | addressed; v3: universal, shock route identically zero, zero given F₂ [A]; sensitivity thresholds stated. **Remaining:** F₂ or a direct value. |
| M6 de-shift, weights | **fully addressed in v3**: X^Q published; continuum k = 2 weights for both functionals; C4 ε_W interval-certified. Remaining minor: C5 ε_W not certified (not load-bearing); Q3\* factored coefficients not a separate data file. |
| M7 Lemma C / Einstein instance | addressed (v2). |
| M8 H_K | **closed in v3**: H_K [D] (EE step derived; f_pos(±3i/2) = 0; criteria 3/2, 5/2). |
| M9 certificate; (1100, 100) | **certificate closed in v3** (≥ N/5 everywhere); (1100, 100) still stated; H_LJ(1000) awaiting HIGHSPIN. |
| M10 H_RS | **closed in v3**: the physical chain is D′ by CLS-D's definition; P_cyc no longer a premise; census claim withdrawn. |
| M11 Δ_g = 300 | addressed (v2); binding state at 300 now [NE±], L = 3.5. |
| M12 title/scope | addressed; v3 replaces A1 by A1′, which ABJM-IIA satisfies at leading order. |
| M13 ABJM | v3: an in-scope consistency check, explicitly not evidence for positivity. |
| M14 novelty/credit | addressed (v2); v3 adds Paulos, Costa–Gonçalves–Penedones, Penedones, Jain et al. |
| M15 jargon, derivations | addressed (v2); v3 adds derivation sketches for the EE kernel, L6′, H_EFT(a), c₁, light states, D′. |
| minor 11 front matter | **remaining**: affiliations, emails, acknowledgments (TODO). |
| overall verdict | positivity is still numerical; v3 makes part of it error-bounded and states exactly which states are model-only. |

## 6. Open items before submission

1. **HIGHSPIN** (in progress elsewhere): fill the `\HS` places (§3) and set `\HSdraftfalse`; if a state is negative, §8.
2. Front matter (affiliations, e-mails, acknowledgments): placeholders "[… to be supplied …]" with `% TODO(authors)` comments in both files.
3. ~~Publish the Q3\* factored coefficients~~ **done in v4** (`data/q3star_factored.json`).
4. Optional: certify C5's ε_W (root isolation of the quintic pole).
5. Before submission, the INSPIRE sweep recommended in LITERATURE.md §6 (v4 added CMT 1710.03199, Häring–Zhiboedov 2202.08280,
   Hartman–Kundu–Tajdini 1610.05308, verified on arXiv and INSPIRE).
6. Unreviewed deliverables: every analysis behind v3/v4 except REPRO / REPRO_ENGINE is unreviewed (referee 2 M3). Highest value:
   independent reproduction of the H_K EE-kernel reduction and of Lemma Y; full write-ups of H_K, L6′, the R⁴ calibration and c₁.
7. `refs.bib` (unused legacy) still has the old JHEP fields; delete it or mirror `refs_verified.bib` before submission.
8. UT3 still carries the stale "natural size 10⁻³–10⁻¹" for c₁ (§3) and "universality [D]" (L3d); the paper now uses [D-sketch] and no
   size estimate. UT3 was not edited (outside this task's files).

## 7. Source inconsistencies

UT3 §11: UT2's twenty items (17 resolved; #9 RD fit 2 % vs 3 % and #20 CHMRS equation numbers still open in the record; #1 procedural)
and fifteen new ones (N1–N15). Items that changed letter text in v3: H_K k = 4 criterion back to 5/2 for our weights; certificate now
1/5 [certified]; Q flat minimum 0.1835 (not 0.1966); certified ε_W; ABJM inside A1′; H_RS by definition; L6′ [D]; H_EFT(a) restricted
negative finding. One grading tension is left as found: PCYC_CLOSE calls H_RS "[D]" with one [D-sketch] sub-step (UT3 N13); the letter
writes "[D, one step D-sketch]".

## 8. Fallback plan if HIGHSPIN returns a negative state (referee 2, M1(d))

Worked out from UT3 §3 ("Dependence"), §5.2–§5.4 and §5.6. Trigger: any heavy state of `Ω^EE` at Δ_g = 1000 comes back with
`λ_min < 0` (the J = 29 edge resolves negative, a near-gap J ≥ 12 state is negative, or the refined H_LJ(1000) fails on a measured
window so that some large-J near-gap state is negative). Then (Pos) for `Ω^EE` is **false** at 1000, and:

| statement | fate | why |
|---|---|---|
| item (1) at 1000 (τ ≤ 2.89·10⁻⁸; \|ĝ₃\|Λ⁴ ≤ 120) | **withdrawn** | its only positivity input is (Pos) for Ω^EE |
| item (1) at 300 | **withdrawn** | already rests on an open window [I]; beyond-Regge corrections at fixed (J, m/Λ) are larger at 300 (∝ Λ⁻² or faster), so a failure at 1000 is expected to recur at 300 for the same class; no positive evidence could rescue it |
| item (2b) (164.5 / 125.7) | **withdrawn** | uses Ω^Q_up + 0.040 Ω^EE and needs (Pos) for Ω^EE |
| Λ-forms of item (1), the conditional remark on Δ_g ≥ 300 | **withdrawn** | values of a now-vacuous implication |
| **item (2a) at 1000** | **survives unchanged** | uses only Ω^Q_up; the exact c₁ form needs only κ > 0, not item (1) |
| item (2a) at 300 | survives formally, not headlined | Q's evidence at 300 is weaker (large mass not converted, L < 1 at m ≥ 3Λ, b ≥ 5) |

**What the Letter would then claim (precisely).** *Assume A1′, A2, the inputs of Table S2 used by item (2) (H_K, L4–L9, LS rows), and
(Pos) for Ω^Q_up on heavy states and on light spin-0, 2 states. If no scalar single trace with m < 0.4602Λ appears in T × T, and either
τ = 0 or κ > 0 [D-sketch + M], then at Δ_g = 1000*

    g̃₄Λ⁶/8πG ≤ 92.81 + 5.44·10⁻⁷ c₁²   (and ≤ 72.04 + 6.81·10⁻⁵ c₁² at Δ_g = 300),

*with g̃₄ = g₄ + O(Λ⁻²) given H_EFT(b) [A]. Implication [D-sketch].* Its positivity evidence (UT3 §5.4): binding pair (++/OO, J = 0 at
the gap) +0.1427 / +0.1428 with s_d = 0.508 against s_fail = 1.776 [NE + M given H_RS]; measured channels min +0.1426 ± 0.0017; unmeasured
channels s_fail ≥ 11 (≥ 18× the measured |s_d|) [NE + M]; 1.3Λ ≤ m ≤ 3Λ, J²/Δ ≤ 1: [NE±], L ≥ 3.95; m > 3Λ: L ≥ 1.65 [NE + M];
J²/Δ > 1 and the large-J near-gap tail: [M] (H_LJ^Q, tolerance ≥ 1.1·10⁴ c₂, margin ×1–10³); HYB [M]. Flat: ≥ 9/50 certified on most
regions, two regions [NE]. Note: with item (1) gone, τ is no longer bounded by the paper, so the sign κ > 0 [D-sketch + M] becomes
load-bearing whenever ĝ₃ ≠ 0 (κτ²Λ⁸/57600 can exceed the upper end by orders of magnitude), and must be stated as such.

**Checks to run before using the fallback.** (a) Evaluate Ω^Q_up on the failing Ω^EE state (same K₂/direct rows, Q's weights): Q's
smeared k = 2, 3 parts use the same EE structure, so a failure mechanism of Ω^EE at near-gap large spin could also act on Ω^Q_up, though
Q's sign-free tolerance there is much larger (≥ 1.1·10⁴ c₂ vs T_B = 1409 c₂ for Ω^EE’s window at (1100, 100)). (b) Do **not** re-tune Ω^EE to absorb the failure before
submission: that adds a further round of post-hoc selection (referee 2, M1(d)); if it is done, it must be disclosed and re-checked on all
computed states, at both gaps.

**Framing.** Title: "A conditional bound on an R⁴ graviton coupling in AdS₄/CFT₃ from ⟨TTTT⟩ light-ray sum rules". Abstract: only
(2a) at 1000, with its condition. The cubic programme becomes a reported negative finding in the Supplement (the state, its value and
grade, and that the fixed Ω^EE is not positive at 1000). PRL suitability would drop further (the surviving result is the
IR-log-dominated g₄ window's upper end, in which ABJM sits at 2.6 %); the referee's PRD/JHEP route becomes the default.

**If HIGHSPIN closes positively** (J ≥ 12 at [NE±], J = 29 positive with a resolved edge, H_LJ(1000) upgraded): fill the `\HS` places,
set `\HSdraftfalse`, and remove weak points (i)–(iii) from the Letter's list accordingly; nothing else changes.
**If it comes back weak but not negative** (budget still not closed): the Letter stands as written; the weak points stay listed.
