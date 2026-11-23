# PRL draft: notes, build, length, provenance, referee response, open items

*Revised 23 Nov 2026 (v3 writing pass; the only new arithmetic is the `c₁²` coefficient of the ε = 0.040 repair). Unrefereed. Nothing
committed. Source of truth: `ads4/UNIFIED_THEOREM_v3.md` (UT3), which integrates EESTEP_CLOSE, HK_RESIDUAL, CERT_CLOSE, POS_ERRORS,
A1PRIME, PCYC_CLOSE, NEARGAP2, C1_VALUE and R4ARC_HEFT on top of `UNIFIED_THEOREM_v2.md` (UT2). UT3 §12 (new numbers) and UT2 §12 (older
numbers) trace every number to a file and section; UT3 §11 resolves UT2's inconsistency list and adds new items. What changed in the
letter and supplement, and the list of HIGHSPIN placeholders: `CHANGES_v3.md`.*

Files in `ads4/prl/`:

| file | content |
|---|---|
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
for f in letter supplement; do
  pdflatex -interaction=nonstopmode -output-directory=build $f.tex
  (cd build && BIBINPUTS=..:. bibtex $f)
  pdflatex -interaction=nonstopmode -output-directory=build $f.tex
  pdflatex -interaction=nonstopmode -output-directory=build $f.tex
done
```

Status (23 Nov, v3): both compile with **no errors, no undefined references or citations, no overfull boxes**. `build/letter.pdf`
**6 pages** in the REVTeX draft layout (main text ends near the top of page 5; references fill the rest); `build/supplement.pdf` 13 pages.
Underfull-box notices only (letter: item I2 and Table I). BibTeX: only "missing number" warnings (issue numbers absent from the bib) and
one `jnrlst` style notice.

## 2. Length

| item | value | method |
|---|---|---|
| abstract | **≈ 560 characters** rendered (pdftotext, incl. spaces; pdftotext inserts extra spaces inside math, so the true count is ≈ 545) | limit 600 |
| main text, PRL counting | **≈ 3080**: prose 2778 words (inline math = 1 word), 9 display lines × 16 = 144, Table I (single-column, 9 rows × 13 + caption + header) = 155 | `scratchpad/tools/wc_prl.py` (same tool as v2) |
| rendered | text occupies ≈ 4.3 REVTeX draft pages before the references | pdftotext per page |

Below the 3750-word limit (v2 was ≈ 2730). The growth (≈ +350) is the light-single-trace paragraph, the error-bounded positivity
sentences and the ABJM rewrite.

## 3. Switches

* **H_LJ(300)**: `\HLJclosedfalse` (open; headline Δ_g = 1000). Unchanged: no new document closes the Regge window for J ≥ 34 at 300.
* **HIGHSPIN**: `\HSdrafttrue` prints a superscript "HS" at every place whose number or grade awaits `ads4/HIGHSPIN.md` (near-gap
  J ≥ 12 at 1000, the J = 29 edge, H_LJ(1000)). Each such place also carries a `% HIGHSPIN:` comment. When the high-spin results land:
  update those places (list in `CHANGES_v3.md` §3), UT3 §5.2/§5.6/§10, and set `\HSdraftfalse`.

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

All other letter numbers: unchanged from v2 (see the v2 table in git history of this file, or UT2 §12).

## 5. Referee report: status after v3

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

1. **HIGHSPIN** (in progress elsewhere): fill the `\HS` places (CHANGES_v3.md §3) and set `\HSdraftfalse`.
2. Front matter (affiliations, emails, acknowledgments).
3. Publish the Q3\* factored coefficients as a data file (outside this task's authorised files); currently derived from
   `q3_full_delta0.2.json` by exact division.
4. Optional: certify C5's ε_W (root isolation of the quintic pole).
5. Before submission, the INSPIRE sweep recommended in LITERATURE.md §6.
6. Unreviewed deliverables: every analysis behind v3 except REPRO / REPRO_ENGINE is unreviewed.

## 7. Source inconsistencies

UT3 §11: UT2's twenty items (17 resolved; #9 RD fit 2 % vs 3 % and #20 CHMRS equation numbers still open in the record; #1 procedural)
and fifteen new ones (N1–N15). Items that changed letter text in v3: H_K k = 4 criterion back to 5/2 for our weights; certificate now
1/5 [certified]; Q flat minimum 0.1835 (not 0.1966); certified ε_W; ABJM inside A1′; H_RS by definition; L6′ [D]; H_EFT(a) restricted
negative finding. One grading tension is left as found: PCYC_CLOSE calls H_RS "[D]" with one [D-sketch] sub-step (UT3 N13); the letter
writes "[D, one step D-sketch]".
