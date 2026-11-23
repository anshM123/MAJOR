# CHANGES v2 → v3 (letter, supplement, macros, bibliography)

*23 Nov 2026. Source of truth for every change: `ads4/UNIFIED_THEOREM_v3.md` (UT3), which cites the deliverable and section for each
item. Pre-v3 copies of the six prl files were not kept in the repository (the v2 state is in git HEAD; the letter and supplement were last committed in 4cdc4bb).*

## 1. Letter (`letter.tex`)

| place | v2 | v3 | why (source) |
|---|---|---|---|
| header comments | `% EESTEP` placeholders; source UT2 | removed; source UT3; HIGHSPIN note | EESTEP_CLOSE, HK_RESIDUAL landed |
| abstract | "large gap"; 92.8 unconditional on light states; "tested numerically on finite sets" | "large higher-spin gap"; 92.8 "if no light scalar couples to TT (164.5 otherwise)"; "tested numerically" | A1PRIME §6; NEARGAP2 (ε = 0.040); positivity now partly error-bounded (NEARGAP2, POS_ERRORS) |
| "What is established" | positivity "computed on finite sets inside a model" | proved in the flat limit, error-bounded on some states, model elsewhere; grade legend adds [certified], [NE±] | CERT_CLOSE; NEARGAP2; POS_ERRORS |
| Setting | A1 (all J); "no known holographic CFT₃ satisfies A1: ABJM has light KK modes" | **A1′** (spin ≥ 4 gap; spin 0, 2 allowed; spins 1, 3 absent from T × T); ABJM sentence removed; g₄ defined with light fields kept in the amplitude | A1PRIME §2, §6, §7 |
| Sum rules | "over …" | adds "over real Lorentzian points"; light spin-0, 2 single traces deferred to their paragraph | PCYC_CLOSE §1 |
| Quartic functionals | "margin δ ≈ 0.2" | "nominal margin δ = 0.2" | CERT_CLOSE §2 (true flat min 0.1835) |
| Cutoff and data | Q continuum weights 1.2·10⁻¹³, 4.2·10⁻¹¹, 2.4·10⁻⁷ [DV]; Ω^EE Abel solutions | all de-shift weights continuum: k = 4 ≤ 2.7·10⁻¹³, 9.1·10⁻¹¹ **[certified]**; k = 5 2.4·10⁻⁷; k = 2 ≤ 1.1·10⁻⁹ relative [DV] | HK_RESIDUAL §3.3; CERT_CLOSE §3 |
| Light side, R⁴ | — | new sentence: R⁴ light function ~ Δφ⁻² log Δφ at coincident angles, integrable [D]; cites Paulos, CGPspin, PenedonesMellin | R4ARC_HEFT §2 |
| Light side, new paragraph | — | **light single traces**: Regge spins 3 / 4; drop-out at k = 4, 5; block action at k = 2, 3 [D-sketch]; Ω^EE positive on all [P]; Ω^Q negative on scalars < 0.46Λ, repair 0.04 Ω^EE [D-sketch] | A1PRIME §2–§4; NEARGAP2 |
| c₁ paragraph | universal number, uncomputed | fixed by ⟨TT⟩, ⟨TTT⟩ [D]; shock part identically zero [P]; energy-channel T exchange (cites JJMNS); zero given one premise [A] (F₂); otherwise uncomputed | C1_VALUE §0–§4 |
| Theorem | A1; (Pos) on heavy states | **A1′**; (Pos) also on light spin-0, 2 states where they enter; item (2) needs "no scalar lighter than 0.46Λ couples to TT", otherwise **125.71 + 1.29·10⁻⁴c₁² / 164.51 + 1.03·10⁻⁶c₁²** with (Pos) for Ω^EE | A1PRIME §4.3; NEARGAP2; c₁² coefficient computed in this pass |
| c₁ sentence | "natural size moves (1) by ≲ 0.02 %; even |c₁| = 100 gives …" | 1 % thresholds: 4.6 / 60 (item 1), 99 / 1262 (item 2, either form) | C1_VALUE §6; this pass |
| Net grades | (1) [D-sketch] through H_K; (2) [D-sketch] H_K, g₅; one link open [I] | (1) [D-sketch] through c₁ ν-independence and (if present) light states' k = 2, 3 action; (2) [D-sketch] in g₅, light states, repair; [D-sketch + M] signs; [A] H_EFT(b); **no input open** | EESTEP / HK_RESIDUAL (H_K [D]); R4ARC_HEFT (L6′ [D]); A1PRIME |
| I1 | [D] | [D; DV on free-field ⟨TTTT⟩] | HK_RESIDUAL §2 |
| I2 | J₀ < 3/2, **J₀ < 2** (k = 4), margins ½, 1; Net [D-sketch], EE step unchecked | EE kernels exact finite sums of shifted scalar partial waves [P; DV]; **J₀ < 5/2** (k = 4; no pole at ±3i/2 [P]); margins ½, **3/2**; **Net [D]** | EESTEP_CLOSE §1–§3; HK_RESIDUAL §1 |
| I4 | g₅ cancellation [D-sketch] "in tension"; R⁴ arc integrability [I] | integrability [D]; g₅ not cancelled (s_d mass-independent Λ → 20Λ [NE]); upper end holds either way | R4ARC_HEFT; EESTEP_CLOSE §5; HK_RESIDUAL §4 |
| I5 | (a) [DV scalar; D-sketch spinning], "fails pointwise in angle" | (a) pointwise for dispersive terms [D; DV] | R4ARC_HEFT §3 |
| Table I | flat "certified ≥ 0.15 except one corner"; Ω^EE min only; H_RS "[D given P_cyc], preliminary" | flat ≥ 1/5 [certified], Q ≥ 0.18 [certified] on most states, min 0.1835; binding ≥ +0.198 / +0.146 [NE±]; J = 0, 2, 4, 6, 8 [NE±]; large-mass strip [NE±]; light J = 0, 2 row; H_RS row dropped (now [D], in text); HIGHSPIN markers | CERT_CLOSE; NEARGAP2; POS_ERRORS; PCYC_CLOSE |
| Positivity paragraph | min values; chronology; H_LJ(1000); H_RS preliminary given P_cyc, far region e^{−26} | flat certificates; binding error-bounded (×78; ×3.5 at 300); near-gap J ≤ 8 (even) and large-mass strip [NE±]; Mehler–Dirichlet bound [P]; "no mass above which all spins are covered"; chronology unchanged; **H_RS: real section cut at the shadow point's light cones = CLS-D's causal domain [P]; O(Δ^N e^{−0.336Δ}) [D, one step D-sketch]; scalar check at Δ = 20 [DV]**; far-region/P_cyc text removed | NEARGAP2; POS_ERRORS; PCYC_CLOSE |
| "Normalisation check" | ABJM outside A1, sign/normalisation only | **"ABJM consistency check"**: IIA satisfies A1′ and A2 at leading order (no light single traces in T × T at tree level); 2.6 % / 3.3 %, one digit; not evidence for positivity | A1PRIME §7 |
| Discussion | next steps incl. proof of P_cyc and the R⁴ arc at coincident angles | next steps: certified positivity, c₁, a quartic functional with light-state positivity imposed | PCYC_CLOSE; R4ARC_HEFT; A1PRIME §9 |

Length: ≈ 2730 → ≈ 3080 PRL-counted words (limit 3750); abstract ≈ 560 rendered characters (limit 600); 6 rendered REVTeX draft pages
(text ≈ 4.3 pages, then references). See NOTES.md §2.

## 2. Supplement (`supplement.tex`)

* **§S1** grades: [certified], [NE±] added, ordering updated; provenance: "EE step being checked" removed, HIGHSPIN marker explained.
* **§S2** normalisations: Λ is the spin-4 gap; g₄ in the minimal light-exchange completion. Cutoff: both k = 2 pieces use
  (Λ + 100, 10) (CERT_CLOSE §3.2; v2 applied Q's (150, 6) to Q's k = 2 piece). Ω^EE: continuum k = 2 weight added as a second
  realisation. Ω^Q_up: **Q3\*** (factored weights, definitional fix of the stored-weight rounding defect), nominal δ = 0.2, true flat
  minimum 0.1835. Positivity blocks: odd J ≥ 5 is 1 × 1 (`+−`), J = 0, 2 have no `+−`, light states included under A1′.
  **Table S1**: `X^Q` and the continuum k = 2 weights added (`data/k2/`, SHA-256 prefixes); caption no longer says X^Q is missing.
* **§S3 Table S2**: A1 → A1′; R1′ [D + DV]; DS row rewritten (k = 4, 5 analytic structure [P], values C4 [D + certified]); LRC-a [D + DV]
  incl. EE k = 4; LRC\* [D + DV] incl. EE; new row f_pos [P]; H_K **[D]**, criteria 3/2 and **5/2**; L3 (universal, shock zero,
  F₂); L4 incl. k = 4 checks; L6′ **[D]**; L7(a) [D]; L8 g₅ [D-sketch + NE]; **new rows LS1–LS6** (light single traces). Net grades and
  "what changed" rewritten.
* **§S4.A** collinear set: the [I] is replaced by the L6′ derivation (gauge, Lemma Y, dispersiveness, non-dispersive counter-example,
  residuals). k = 4 borderline/drop-out checks added (A1PRIME §5).
* **§S4.C**: H_EFT(a) [D] for dispersive terms; the "fails at fixed angle" finding restricted to non-dispersive toy classes; (P) for g₅:
  deciding test results at J = 0 and J = 2, exact cancellation excluded; C5 piece not run [I].
* **§S4.D** H_K: spinning R1′ check; strip analyticity by construction for the continuum weights; **EE kernel reduction** (replaces "EE
  step (pending)"); f_pos(±3i/2) = 0; conclusion **[D]**, J₀ < 5/2 for k = 4, margins ½ and 3/2.
* **§S4.E** c₁: shock part identically zero (CHP operators in d = 3) and blind to c₁; ASD closure (no independent seagull); F₂ and
  c₁ = 0 given F₂; 1 % thresholds; **new Eq. (S-repair)** for Ω^Q_up + εΩ^EE with ε = 0.040.
* **§S4.F** de-shift: certified ε_W for C4; k = 2 shifted form and continuum ε_W; X^Q re-solve; "not certified" disclaimer replaced.
* **§S4.G (new)** light single traces under A1′: occurrence, Regge spins, drop-out, flat sign polynomials, consequences for both items,
  repair, KK towers.
* **§S5** positivity: method (Lemma C [P]; K₂ rows; [NE±] definition); flat certificates (NF02 ≥ 1/5 everywhere; Q3\* ≥ 9/50 on listed
  regions); **Table S5 rewritten** as "error-bounded and computed sets" (binding, near gap, J ≥ 12 with HS, large mass, 300, Q, light);
  Table S6 J = 15 = +0.6900, HS markers on J = 12, 16, 29; **Table S7** covering hypotheses: H_RS [D], P_cyc row removed, H_M/H_edge/H_R′
  restricted, H_bulk and Pos_light added, σ_hU "consistent reading". H_RS paragraph rewritten (physical chain = D′; Prop. D′; end-to-end
  scalar check; crystal cycle a numerical method; census claim withdrawn). Record notes: 300 own row; s_d reference artefact; odd-J 1 × 1.
* **§S6** corollary: A1 → A1′; H_mono at J = 0 [D given RD] + [NE].
* **§S7** reproduction: X^Q re-solve; certified items; Q flat minimum record fix.
* **§S8** ABJM rewritten as a consistency check inside A1′, A2 at leading order (CPY/BCP exchange content; KK orthogonality; 𝒩 = 6 caveat;
  one-digit placement; not evidence).
* **§S9** not-claimed list extended; record corrections list rewritten for v3.

## 3. HIGHSPIN placeholders remaining (fill from `ads4/HIGHSPIN.md`)

Every item below prints the marker `\HS` (superscript "HS") while `\HSdrafttrue`, and has a `% HIGHSPIN:` comment in the source. Current
content is the grade as of NEARGAP2 / POS_ERRORS / HLJ_CLOSE; nothing was invented.

| # | file:line (at writing) | what is there now | what goes there |
|---|---|---|---|
| 1 | `letter.tex` Table I caption (≈ l. 306) | "HS: being recomputed." | delete when all items are filled |
| 2 | `letter.tex` Table I, row "Ω^EE, 1000" (≈ l. 311) | "others: computed sets, even J ≤ 48, odd ≤ 31" ([NE + M]) | grade (and leverage) of near-gap J ≥ 12 at Δ_g = 1000 |
| 3 | `letter.tex` Table I, row "H_LJ(1000)" (≈ l. 316) | [D-sketch + M] | H_LJ(1000) grade / leverage |
| 4 | `letter.tex` Positivity paragraph (≈ l. 335) | "Other near-gap spins … computed on finite sets inside a model" | status of near-gap J ≥ 12 at 1000 |
| 5 | `letter.tex` Positivity paragraph (≈ l. 343) | "holds only in a refined form [D-sketch + M]" | H_LJ(1000) grade |
| 6 | `supplement.tex` Table S5 caption and row "Ω^EE, 1000, J ≥ 12" (≈ l. 633, 641) | J = 12, 16 own-row values; [NE + M], L ≈ 1–3 | grade / leverage for J ≥ 12 (and any new J) |
| 7 | `supplement.tex` Table S6 caption and Δ_g = 1000 row (≈ l. 654, 663–665): J = 12, 16, 29 | +0.6759, +0.7191; "29: edge window unresolved" | updated J = 12, 16 values; J = 29 edge value |
| 8 | `supplement.tex` Table S7 caption; rows H_M, H_edge, H_LJ(1000) (≈ l. 672, 679, 680, 683) | [M] / J = 29 unresolved [M] / [D-sketch + M], leverage 2.7–4.8 | H_M for J ≥ 12; J = 29 edge grade; H_LJ(1000) grade / leverage |
| 9 | `UNIFIED_THEOREM_v3.md` §5.2 (4 rows), §5.6 (H_LJ(1000), H_M, H_edge, net paragraph), §10 | ⟦HIGHSPIN⟧ with current grades | as above |

When done: set `\HSdraftfalse` in `macros.tex`, and remove the "HS: being recomputed" sentences from the three captions.

## 4. Bibliography

* `refs_verified.bib`: added **Paulos** (1107.1504, JHEP 10 (2011) 074), **CGPspin** (Costa–Gonçalves–Penedones 1404.5625, JHEP 09 (2014)
  064), **PenedonesMellin** (1011.1485, JHEP 03 (2011) 025), **JJMNS** (Jain–John–Mehta–Nizami–Suresh 2106.00016, JHEP 09 (2021) 041).
  Each verified on 23 Nov 2026 against the arXiv API (title, authors, journal_ref, DOI) and INSPIRE (publication info, DOI, authors;
  recids in the file). Costa–Hansen–Penedones 1707.07689 was already present as `CHP` and was re-verified.
* `refs.bib` (unused): the same four entries mirrored, with a note.

## 5. Macros (`macros.tex`)

* Source pointer → UT3. `\EESTEP` removed (EE step closed). New `\QuartRep` (repaired upper end 164.5 / 125.7, follows the H_LJ(300)
  switch). New HIGHSPIN marker `\HS` with switch `\HSdrafttrue`. H_LJ(300) switch unchanged (open).
