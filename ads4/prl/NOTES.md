# PRL draft: notes, provenance of every number, TODOs, H_LJ switch

*Drafted 22 Nov 2026 (writing pass only; no new computation). Unrefereed. Nothing committed; no existing file edited.*

Files in `ads4/prl/`:

| file | content |
|---|---|
| `letter.tex` | the Letter (REVTeX 4.2, `aps,prl,twocolumn,superscriptaddress`) |
| `supplement.tex` | Supplemental Material (REVTeX 4.2, one column) |
| `macros.tex` | shared macros, including the **H_LJ(300) switch** (read by both `.tex` files; added so that the switch is one line for both documents) |
| `refs.bib` | bibliography; every entry is a work named in the repo docs; unconfirmed details are marked `TODO(bib)` |
| `NOTES.md` | this file |
| `build/` | all build outputs; `build/.gitignore` ignores everything there except itself |

**Source of truth:** `UNIFIED_THEOREM.md` (UT) including its six Update boxes (HRS_CLOSE, EINSTEIN_CHECK, DG_SCAN, REPRO,
REPRO_ENGINE, RD_CLOSE), which supersede the UT body where they conflict. Abbreviations below: RD-S = `RESULT_D_STATEMENT.md`,
RQ-S = `RESULT_Q_STATEMENT.md`, DG = `DG_SCAN.md`, RP = `REPRO.md`, RPE = `REPRO_ENGINE.md`, HRS = `HRS_CLOSE.md`,
RDC = `RD_CLOSE.md`, EIN = `RESULT_Q_EINSTEIN_CHECK.md`, ABJM = `ABJM_THEOREM.md`. `HLJ_CLOSE.md` and `code/spinfunc/lj_*`
were not read.

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

* TeX was not installed initially. A minimal `apt-get install --no-install-recommends texlive-latex-base texlive-publishers
  texlive-latex-recommended texlive-fonts-recommended` (≈ 74 MB download) succeeded and provides `pdflatex`, `bibtex` and
  `revtex4-2.cls`.
* Status at hand-off: both documents compile with **no errors, no undefined references, no overfull boxes**
  (`build/letter.pdf`: 5 pages; `build/supplement.pdf`: 7 pages). BibTeX warnings only for missing fields of entries marked
  `TODO(bib)` (empty year in OsbornPetkou; missing issue numbers).
* The closed position of the switch (`\HLJclosedtrue`) was also test-compiled in a scratch copy: compiles; only overfull
  warnings come from the placeholder `[TODO(HLJ)]` in narrow table cells, which disappear once the grade is filled in.

## 2. Length

| item | value | method |
|---|---|---|
| abstract, fallback switch (Δ_g = 1000) | **561 characters** (rendered, from `pdftotext`, incl. spaces) | limit 600 |
| abstract, closed switch (Δ_g = 300) | **560 characters** | same |
| main text, PRL convention | **≈ 3450 word-equivalents**: prose 2548 words (inline math = 1 word; the three switch macros add ≈ 115 words, not in this count), 9 display-equation lines × 16 = 144, Table I 23 rows × 26 = 598, caption 41 | `scratchpad/tools/wc_prl.py` + manual add of macro text |
| main text, rendered | 3647 words (`pdftotext`, Introduction → acknowledgments, incl. table text and equation tokens) | upper estimate |

Both estimates are below the ≈ 3750-word PRL limit.

## 3. H_LJ(300) switch (instructions)

`H_LJ(300)` = Result D's near-gap large-spin hypothesis at Δ_g = 300. Its Regge-window part for `J ≥ 34` is **[I]** in UT (§2
Grade, §3.2, §4 "Open", §5 item 2) and RD-S §5.3. Its outcome (HLJ_CLOSE.md, other agent) was unknown when this was written.

**The switch is one line in `macros.tex`:**
```
\HLJclosedfalse   % <<< SWITCH      (current: fallback, headline Δ_g = 1000)
\HLJclosedtrue    % <<< SWITCH      (window closed, headline Δ_g = 300)
```

What flips automatically (both documents):
* the headline Δ_g in the abstract (`\DgHead`), the two headline bounds (`\CubicHead`: 2.89·10⁻⁸ ↔ 3.13·10⁻⁶; `\QuartHead`:
  92.81 ↔ 72.04) and the ABJM placement (`\ABJMHead`: 2.6 % ↔ 3.3 %);
* the status sentence `\HLJstatus` at the end of the letter's introduction and in Supplement §S1:
  * (a) closed: *"The Regge window of H_LJ(300) (near-gap spins J ≥ 34), open in the consolidated statement, has since been closed
    at grade X (HLJ_CLOSE.md, §…); this is why we headline Δ_g = 300."*
  * (b) open: *"The Regge window of H_LJ(300) (near-gap spins J ≥ 34) is still an open hypothesis [I]: there are no graviton window
    data at these spins, and the sign-free tolerance exceeds the measured window corrections only after their leading-Regge-shaped
    part is removed (a [D-sketch] step), and then by a factor of just 1.1. Every use of Ω^EE positivity at Δ_g = 300 is therefore conditional on an open input, and we headline Δ_g = 1000."*
* the grade printed in Table I (letter) and Table S2 (supplement) for the window (`\HLJwindowGrade`: `[I]` ↔ `\HLJgrade`);
* the corollary's inheritance clause (`\HLJcorollary`) and the "weakest links" sentence (`\HLJweakest`).

What must also be done **by hand** if the window closes (all marked `TODO(HLJ)`):
1. set `\HLJgrade` (the grade X) and `\HLJsource` (section of HLJ_CLOSE.md) in `macros.tex`;
2. check whether HLJ_CLOSE also covers the separate `(336, 28/32)` window below ν = 8.3, graded [M] (RD-S §5.3, UT §3.2). If it
   does, edit the H_LJ(300) rows of Table I and Table S2 and `\HLJweakest` by hand; if it does not, leave them;
3. check the "Weakest links" paragraph and the "What would make it stronger" list in the Discussion (first item: graviton window
   data at 300 for J ≥ 34), and drop that item if it is done;
4. the bodies of the letter keep both Δ_g = 300 and 1000 numbers in the theorem regardless of the switch; nothing else changes.

Things that do **not** depend on H_LJ(300): every Δ_g = 1000 statement; Result Q (item 2) at Δ_g = 300 as long as `c₁t₄` is carried
explicitly (UT §2 Remark 1). The letter says so.

Option *not* used: in fallback mode one could add "for Δ_g ≥ 1000 without H_LJ(300)" by applying DG §2.5's general statement
("positivity at Δ_g₀ implies positivity at Δ_g ≥ Δ_g₀") with Δ_g₀ = 1000. DG does not state this instance explicitly, so it is not in
the draft.

## 4. Every number in the letter, with its source

### Abstract and headline
| number | source |
|---|---|
| √(t₄²+α_T²) ≤ 2.89·10⁻⁸ (1000) / 3.13·10⁻⁶ (300) | UT §2 item 1; RD-S §0, §3 table; RP §1 (2.89126·10⁻⁸ / 3.13206·10⁻⁶) |
| upper end 92.81 (1000) / 72.04 (300) | UT §2 item 2; RQ-S §0; RP §1 (92.8135 / 72.0406) |
| lower end −O(Λ⁻²) | UT §2 item 2; RQ-S §2.4 |
| Λ ≃ Δ_g − 3/2 | DG header (`Λ = Δ_g − 3/2 − 1/(8Δ_g) + O(Δ_g⁻²)`) |
| ABJM 2.6 % / 3.3 % | UT §7; RQ-S §6; RP §1 (2.59 % / 3.34 %); DG §1.2 table |
| "one uncomputed light-side term" = c₁ | UT §1 table, §2 |

### Setting
| number | source |
|---|---|
| c_T = 1 for a free real scalar; c_T = 32L²/(πG); 8πG = 256/c_T | UT §1 table; RD-S §1; RP §1 last row |
| A1 incl. J = 0, 1, 2; A2 c_T ≫ Δ_g² | UT §1 table |
| m² = (Δ−J−2)(Δ+J−1) | UT §1 table |
| Λ = √((Δ_g−2)(Δ_g−1)) = Δ_g − 3/2 − 1/(8Δ_g) + O(Δ_g⁻²); 298.50 / 998.50 | UT §1 table; DG header |
| √(t₄²+α_T²) = 240|ĝ₃| (CHLPS normalisation) | UT §1 table; RD-S §1 |
| g₄/(8πG) = 2ζ(3)/M⁶, M² = 4/α′ | UT §1 table; RQ-S §1 |

### Method
| number / formula | source |
|---|---|
| EE = [V_{0,42}²(H₂₄ − 2V_{2,04}V_{4,02})/(X₂₀X₄₀)]² | RD-S §2 table |
| helicity off-diagonal vanishes identically [P] | RD-S §2 Properties; UT §4 |
| EE ≡ q² on CLS-D's plane [P]; flat-limit identification [P given map] + [NE]; flat min 0.200 [certified], one corner float-checked | RD-S §2 Properties; UT §4 |
| Ω^EE formula; k(0) = 1 | UT §1 "fixed functionals"; RD-S §2 |
| L0 [D]; −Ĉ′_EE|_E = 16πG/(ν²+9/4) (H7 [D]) | UT §3.1; RD-S §3 |
| Ĉ′_EE := c(ν)⁻¹𝒜_g⁻¹[Ψ_EE], Abel equation, min-norm gauge | RD-S §2 table |
| de-shift existence: C2 [D]; C4/C5 [NE] (to t = 30) | UT §3.1 "De-shift" row |
| δ = 0.2; log coefficient 17.241 | UT §1 (A = 17.241); RQ-S §1 |
| lower functional: P_{−½}(cosh β) > 0 | RQ-S §1 |
| Lemma A (J₀ < 4: Einstein, E×W³, MFT) [D]; D_LR [D]; domain x₁₃² > 0 | RQ-S §2.1; UT §3.1 |
| light side = D_LR on J₀ < 4 + term-wise (F)+(P) on J₀ ≥ 4 [D-sketch]; full-correlator D_LR would give R⁴ zero | UT EINSTEIN_CHECK box; EIN §0 item 5, §5 |
| −Ĉ₄^{ds}|_light = 2g̃₄ | UT §1 table; RQ-S §1 |

### Theorem (letter Eqs. (3)–(4))
| number | source |
|---|---|
| E×W³ light term 16πG c₁t₄; c₁ free | UT §1 table |
| L_A = ∫ψ₂/(p² + 9/(4Λ²)) | UT §2 item 1 |
| bound 240√(16L_A/(k(0)(1+δ)))/Λ⁴(1 + O(c₁Λ⁻²)) | UT §2 item 1. (The letter writes δ₃ for UT's δ, to avoid a clash with the margin δ = 0.2; per RD-S §3, δ = δk(0) + H17 remainder.) |
| 3.13·10⁻⁶ / 2.89·10⁻⁸ | UT §2; RP §1 |
| ∓0.22 % / ∓0.017 % per unit c₁ | UT §2; DG §1.2 (0.216 % / 0.0167 %); RP |
| ±0.24 % / ±0.021 % (H17) | UT §2 (corrected values, REPRO box); DG §1.2; RP D2 |
| |ĝ₃|²Λ⁸ ≤ 1.07·10⁴ / 1.43·10⁴ | UT §2; RP §1 |
| 17.241 ln(2Λ/3) + B_eff; B_eff = −19.220 / −19.265 | UT §2 item 2; DG §1.2 (−19.2197 / −19.2647) |
| 72.04 / 92.81 | UT §2 |
| |ε_c₁| = 2.03|c₁t₄|Λ² | UT §2 (2.03); DG §1.3 (2.032) |
| |ε_c₁| ≤ 0.57|c₁| / 0.059|c₁| | UT §2; DG §1.2 (0.567 / 0.0586) |
| κ = 0.110 (only tightens) | UT §2 |
| κ ≥ 0, c_L > 0 [D-sketch + M]; slope 0.255–0.49 | UT §2 Remark 2; UT §3.1 |
| item (1) needs no Lemma A/(F)+(P)/H_EFT/forward; item (2) needs (1) only via ε_c₁ | UT §2 Remark 1 |
| no cubic closure | UT §2 last line, §6 |

### Λ-forms (letter Eqs. (5)–(6))
| number | source |
|---|---|
| 1.3194·10⁴, 8.7415, 6.384, O(lnΛ/Λ²) | UT DG_SCAN box; DG §0 item 1, §3 |
| 17.2407, −19.2782, 10.314, 860.3, 2508.6, O(Λ⁻³lnΛ) | UT DG_SCAN box; DG §0, §1.1 |
| "exact closed forms reproduce every quoted digit" | DG §0 item 1, §1.2 "Reproduction"; RP §1 |
| truncated cubic form error ≲ 5·10⁻⁵; quartic ≤ 10⁻³ | DG §1.3 (+4.5·10⁻⁵ / +3.6·10⁻⁶; −0.0009 / −0.00002) and RP `code/repro/cache/bounds.txt` (+4.45·10⁻⁵ / +3.1·10⁻⁶; 1.3194e4-form +3.0·10⁻⁵ / −1.1·10⁻⁵; quartic −0.00103 / −0.00014). The bounds quoted in the letter cover both. |
| B_eff drift −19.220 → −19.265 → flat −19.2782 | UT DG_SCAN box; DG §0 |
| ∓19941/(Λ²√(16L_A)); ±212.5/Λ² [M]; 487.7√(16L_A)|c₁|/Λ² | UT DG_SCAN box; DG §0 item 1, §1.3 |
| "√lnΛ/Λ⁴"; "replacing Λ by Δ_g is not innocuous at 300" | DG §1.3, §3 (the −2.2 % of DG is −1.99 % per RP D6; the letter quotes no number, the supplement −2.0 %) |

### Table I (letter)
All grades from UT §3.1–3.2 as updated by the boxes: H_RS → [D given P_cyc] (RDC box); Lemma A mechanism on Witten diagrams (EIN
box); light-side definition (EIN box); H_mono/H_cont (DG box; DG §2.5; RD-S §4 H_cont row). K2/K3 from UT §3.1. "J₀ ≤ 1", "J₀ < 7/2",
"J ≥ 34", "(336, 28/32)", "ν < 8.3": UT §3.1–3.2, RD-S §5.3.

### Status of positivity
| number | source |
|---|---|
| computed sets (1000: even J ≤ 24, 28, 32, 36, 40, 48; odd J ≤ 31. 300: even J ≤ 24, 28, 32; odd J ≤ 15); twist edge to 5Λ | UT §3.2; RD-S §0 item 2, §5 |
| +0.2038 at 1.0166Λ; +0.2009 at 1.00499Λ | UT §2 Remark 3; RD-S §5.2–5.3; twist edge 1.00499Λ per REPRO box (record fix) |
| flat 0.200 | RD-S §2; DG §2.3 |
| J = 32 (300): +0.783 | RD-S §5.3 |
| s_d = +0.5096 / +0.5077 | UT HRS_CLOSE box / HRS §4 (K = 5); RPE §0 (+0.50967 / +0.50776) |
| s_fail = 1.78; margin ≈ 3.5 | UT §2 Remark 3; RP §2 (1.7758/0.5096 = 3.48) |
| FF/N +0.1427 / +0.1428 | UT HRS_CLOSE box; HRS §0; RP; RPE |
| full-set minimum +0.1426 ± 0.0017 | UT §2 Remark 3; RQ-S §0, §3.1 |
| measured channels (++ J = 0–8 even; OO J = 0, 2; +− J = 4–9, 11); unmeasured need ≥ 17–30× | UT §3.2 "computed set, Q" |
| contiguity argument [D], hypotheses [DV] | UT HRS_CLOSE box; HRS §0, §2 |
| c = 2 artanh(r/2) ≥ 0.336 [D] | UT RD_CLOSE box; RDC §0, §2.3 |
| ≤ e^{−37} / e^{−171} vs e^{−8.5} / e^{−10.8} [NE + D-sketch] | UT RD_CLOSE box; RDC §0, §3.2 |
| EE = a·Q·scalar, no new divisor [P + DV] | UT RD_CLOSE box; RDC §1 |
| P_cyc, [DV] for scalar via periods | UT RD_CLOSE box; RDC §4.1 |
| |Δs_d| ≤ 0.015 / 0.0015; ≥ 85× | UT HRS_CLOSE box; HRS §6 |

### Independent reproduction
| number | source |
|---|---|
| 3.13206·10⁻⁶ / 2.89126·10⁻⁸; 72.0406 / 92.8135 | UT REPRO box; RP §1 |
| ρ agreement ≤ 3.4·10⁻⁷ / 3·10⁻⁸ | RPE §0, §1 (the UT box rounds the first to "3·10⁻⁷") |
| s_d +0.50967 / +0.50776 | UT REPRO_ENGINE box; RPE §0 |
| ~600× | UT REPRO_ENGINE box; RPE §0, §5 |
| not reproduced: de-shift weights, D's crystal corrections, NF02 certificate (grid only), full Q scan | UT REPRO box; RP §3 |

### Weakest links
| number | source |
|---|---|
| free scalar ≤ 8·10⁻¹⁵ (free fermion ≤ 3.8·10⁻¹⁵ in supplement) | UT §3.1 Lemma A row; RQ-S §2.1 |
| Witten diagrams 10⁻¹⁴–2·10⁻¹²; J₀ = 2 spin-2 exchange, τ = 2; controls fail at O(1); not exercised τ = 1 etc. | UT EINSTEIN_CHECK box; EIN §0 |
| ⟨TTTT⟩ only in momentum space (Raju; Albayrak–Kharel) | EIN §2.1 |

### H_LJ(300) status sentence (b)
| number | source |
|---|---|
| J ≥ 34; no graviton window data | UT §3.2, §5 item 2; RD-S §5.3 |
| ×1.1 after removing the D_h (leading-Regge-shaped) class [D-sketch]; full measured size ≈ 450 c₂ would give ×0.62 | DG §2.2 (P2), §2.3 table, §2.4; RD-S §5.3 ("refined remainder ≈ 150–250 at the level of T_B ≈ 250") |

### Gaps above 300 / corollary
| number | source |
|---|---|
| only 300 and 1000 computed; monotone margins; baseline → 0.200; T_B ∝ Λ^{1.0–1.2} | DG §0 item 2, §2.5 |
| H_mono, H_cont statements; grade [D-sketch given H_mono + H_cont + H_LJ(Δ_g)], H_mono [M] | DG §2.5 Conclusion, §3 grade table; UT DG_SCAN box |
| ×1.1 at 300 → ×2 at ≈ 490 | DG §3 Corollary; DG §2.3 (P2 column) |
| 300 is the lowest gap where the method closes, marginally; Δ_g = 100 fails at model level | DG §2.4; UT DG_SCAN box; UT §6 |

### ABJM
| number | source |
|---|---|
| M² = 4/α′, L²/α′ = π√(2λ), ML = 2^{5/4}π^{1/2}λ^{1/4} ≈ 4.22λ^{1/4} | UT §7; RP §1 (4.21563) |
| ζ(3)α′³/32 ⇒ 2ζ(3) = 2.404 | UT §7; RP §1 (Virasoro–Shapiro check) |
| 3.3 % / 2.6 % | UT §7 |
| ĝ₃ = 0 from 10d IIA SUSY; N = 6 at k > 2 | UT §7 |
| c_T not in the placement; ∂²_m F flavour-current caveat; c_T = (64/π²)∂²_m F needs a stress-tensor-multiplet current | UT §7; ABJM "Stress tensor versus flavour current"; RQ-S §6 |

### Discussion
| number | source |
|---|---|
| 1.32·10⁴√ln(Λ/8.74)/Λ⁴ | DG §1.3 "Leading behaviour" |
| 17.24 ln(Λ/4.59) + O(1/Λ) | DG §1.3 table (17.2407 ln(Λ/4.5889)) |
| |ĝ₃|²Λ⁸ ≤ 3022.34 ln(2Λ/3) − 5327.21 + O(1/Λ); ~10× weaker than CHLPS optimum | UT §8; RD-S §0 item 4, §6.4 (×11.0) |
| AKT ln Δ_g/Δ_g⁴, constants not comparable | UT §8; RD-S §6.4 |
| novelty at search level (abstract-level search plus one extended search) | UT §8 |

### Supplement-only numbers (beyond the above)
* Table S1/S2 grades and evidence: UT §3.1–3.2, RD-S §4, RQ-S §4, updated by boxes (HRS, RDC, EIN, DG). σ_hU: UT §3.2.
* Table S3: DG §1.2 (all bound constants, δk(0) factor, H17 band, c₁, ε_c₁, ABJM); RD-S §5.2–5.3 (D states, J = 48 +0.915 ± 0.08);
  RP §1 (hybrid +0.2043); RQ-S §0/§3.1 (+0.1426 ± 0.0017); HRS §4, §6 (s_d ± truncation, ×85/×850); RDC §3.2; DG §2.3 (×1.1, ×4.8).
  `s_fail/s_d ≈ 3.5` at 1000 is my arithmetic 1.7758/0.5077 = 3.50 (DG's 3.55 used the superseded s_d law).
* Table S4 (implication at Δ_g = 300, 500, 1000, 2000, 5000): DG §1.2 table, rows ≥ 300 only.
* §S4 Λ-forms: DG §1.1, §1.3 (constants A = 3022.3377, B_∞ = −5327.2096, B_AdS = −6552.662, b₁ = 19293.4, b₂L = −13549, b₂ = −14705,
  Λ₀^D = 8.7415; A_Q = 17.240713, B_∞ = −19.278176, b₁ = 10.3137, b₂L = 860.35, b₂ = −2508.6, Λ₀^Q = 4.5889; drift split 0.0450 =
  0.0242 + 0.0199 + 0.0009; H17 424.97/Λ²; δk(0) = −68.08/Λ², exact −68.3/−68.1, factor 1 + 34.0/Λ²); "−2.0 %" for Λ⁴ → Δ_g⁴ from RP D6.
* Table S5 (D per-state values): RD-S §5.2, §5.3 tables. Caveats on J = 15 (1000, ≈ +0.69) and J = 16 (300, +0.805 at 3.17Λ):
  RP D4, D5.
* Q channels (s_fail values): RQ-S §3.1 table.
* H_RS details: HRS §0, §2; RDC §0–§4 (scalar periods ≤ 6·10⁻⁹ / 1·10⁻⁹: RDC §4.1(d)); angular corner: RDC §3.4.
* Lemma A checks: RQ-S §2.1 (25 components, ≤ 8·10⁻¹⁵, fermion ≤ 3.8·10⁻¹⁵); EIN §0, §3 (Witten, Δφ = 0.025, one-sided weight
  control, J₀ ≤ k, exponent k − 1 − J₀).
* Corollary support: DG §2.5 (Λ⁻², Λ⁻⁴, p = 2.83 (1.4–3.2), excess 0.0043/0.0009/0.00016), DG §2.6 (easier/harder list), s_d law
  0.5076 + 147/Δ² (UT DG_SCAN box; RP D1). Δ_g = 100 failure at J ≳ 18: UT §6 (see inconsistency 15). Window closes from ≈ 275: DG §2.4.
* Table S6 (reproduction): RP §1, §3; RPE §0, §1, §5. Record fixes: UT REPRO box; RPE E1.
* What is not claimed: UT §6; RQ-S §5.3 (withdrawn closure, 0.49 vs 0.26, 200/257); RQ-S §0 (δ = 0.05, +0.0060).

## 5. TODO list

**Front matter (letter and supplement)**
* `TODO: affiliation` (both authors; letter and supplement), `TODO: email` (both authors; letter).
* `TODO: acknowledgments and funding` (letter).
* Supplemental Material URL placeholder "[URL will be inserted by publisher]" (standard PRL wording; nothing to do).

**H_LJ(300)** (`macros.tex`): the switch, `\HLJgrade`, `\HLJsource`; see §3 above.

**Content**
* Letter, Discussion: `% TODO: expand if space allows`: UT §8 and RQ-S §7 say Kundu (2104.11238) and the N = 6 bootstrap
  (Binder–Chester–Jerdee–Pufu, 2011.05728) should be *discussed, not only cited*. The draft only cites them. No repo document contains
  the discussion, so it was not written.

**Bibliography** (`refs.bib`, `TODO(bib)`): given names/initials for every entry (the repo records surnames only), and:

| key | arXiv | confirmed in repo | missing |
|---|---|---|---|
| HPPS | 0907.0151 | authors (LITERATURE.md) | title, journal |
| HofmanMaldacena | 0803.1467 | authors (DICTIONARY_T4.md) | title, journal |
| Buchel | 0911.4257 | authors (DICTIONARY_T4.md) | title, journal |
| OsbornPetkou | — | authors only (DICTIONARY_T4.md, "Osborn–Petkou normalisation") | everything else |
| CEMZ | 1407.5597 | authors | title, journal |
| AHKT | 1610.09378 | authors | title, journal |
| CHP | 1707.07689 | authors | title, journal |
| MeltzerPerlmutter | 1712.04861 | authors | title, journal |
| AKT | 1805.07393 | authors | title, journal |
| KKSDZ | 1904.05905 | authors, title (EINSTEIN_OBSTRUCTIONS.md) | journal |
| CHMRS | 2106.10274 | authors, title (LITERATURE.md) | journal |
| CHLPS | 2201.06602 | authors, JHEP 05 (2023) 122 | title |
| CLSD | 2311.04271 | authors, JHEP 04 (2025) 016 | title |
| CaronHuotLi | 2102.08160 | authors, title as paraphrased in EIN §2.1 | exact title, journal |
| Raju | 1201.6452 | author, title (EIN §2.1) | journal |
| AlbayrakKharel | 1908.01835 | authors, title, JHEP 12 (2019) 135 | — |
| ABJM | 0806.1218 | authors | title, journal |
| Kundu | 2104.11238 | author | title, journal |
| BCJP | 2011.05728 | authors | title, journal |
| ACP | 1711.07343 | authors | title, journal |
| BCPW | 1906.07195 | **author list inconsistent** (LITERATURE.md: Binder–Chester–Pufu–Wang; ABJM_THEOREM.md: Binder–Chester–Pufu) | authors, title, journal |
| Fmax | 1205.4142 | arXiv number only (ABJM_THEOREM.md) | authors, title, journal |

## 6. Where the source documents were inconsistent or unclear (and what the draft does)

1. **Grade of H_RS.** UT body (§2 Grade, §3.2, §4, §5 item 1) and RQ-S body (§0 Status, §3.2, §4 #11, §5) say **[I]**; the HRS_CLOSE
   box says [D-sketch] ([D given RD]); the REPRO_ENGINE box says "H_RS stays [D-sketch]"; the later RD_CLOSE box says **[D given P_cyc]**.
   Draft: [D given P_cyc], with P_cyc [DV] for the scalar via periods and "not proved", and the source marked unreviewed (†). Hence the
   UT body's "two items are [I]" is no longer true: of positivity inputs only the H_LJ(300) window (and the 300 J-drift of H_fwd^Q,
   "harmless") remain [I].
2. **"Δ_g ≥ 300".** UT §6 says "'Δ_g ≥ 300' is not claimed" and RD-S App. A.3 lists "Δ_g ≥ 300" as withdrawn; the DG_SCAN box (later)
   proposes it as a conditional corollary [D-sketch + M]. Draft: stated only as a separate, conditional corollary.
   The DG_SCAN box lists only H_mono and H_cont as extra hypotheses, while DG §3 (grade table) also requires **H_LJ(Δ_g)**; the draft
   includes H_LJ(Δ_g).
3. **s_d law.** DG §2.7/§3 give `s_d ≈ 0.500 + 1440/Δ²` and Q margins ×2.76 … ×3.55 derived from it; RP D1 shows the law misses its own
   anchors; the corrected `0.5076 + 147/Δ²` is in the UT DG_SCAN and REPRO boxes. Draft uses the corrected law and measured margins
   (3.48 at 300; 3.50 at 1000 by my arithmetic), and none of DG's sub-300 margins.
4. **Q binding value at 1000.** RQ-S §3.1 says "the 1000 state +0.1426 ± 0.0002"; HRS §6, RP and RPE give +0.1428 (+0.14282). Draft
   uses +0.1428 for the binding state and RQ-S's +0.1426 ± 0.0017 for the full-set minimum.
5. **ρ agreement at 300.** UT REPRO_ENGINE box "≤ 3·10⁻⁷"; RPE §0/§1 "≤ 3.4·10⁻⁷". Draft: 3.4·10⁻⁷.
6. **RD far-region fit.** UT RD_CLOSE box "within 2 %"; RDC §0/§2.4 "within 3 %". Draft: 3 %.
7. **Accuracy of the truncated Λ-forms at 1000.** DG §1.3 quotes −0.00002 (quartic, all terms) and +3.6·10⁻⁶ (cubic); RP's recomputation
   (`code/repro/cache/bounds.txt`) gives −0.00014 and +3.1·10⁻⁶ (and +3.0·10⁻⁵ / −1.1·10⁻⁵ for the rounded 1.3194·10⁴ form). The RP
   summary table compares only 300. Draft quotes conservative bounds (≲ 5·10⁻⁵ relative, ≤ 10⁻³ absolute).
8. **Λ⁴ → Δ_g⁴.** DG §3: −2.2 % at 300; RP D6: −1.99 %. Letter gives no number; supplement says −2.0 %.
9. **H17 at 1000.** RD-S §3 table "±0.02 %" vs ±0.021 % elsewhere (the UT corrected value). Draft: ±0.021 %.
10. **s_d error at 300.** RQ-S §3.1 "±0.016"; HRS "±0.015". Draft: ±0.015.
11. **Flat-limit grade.** UT §4 lists the crystal identity EE ≡ q² as [P]; RD-S §2 lowers "Flat" to [P given H0 §17 map] + [NE] because
   in the full null-circle contraction EE/q² → 1 only as m → ∞. Draft: identity [P], identification [P given the polarisation map] + [NE].
12. **"Every other channel needs 17–30×"** (UT §3.2) refers to *unmeasured* channels; measured ++/OO J = 2 have s_fail/s_d ≈ 5.5
   (RQ-S §3.1). Draft words it as "every unmeasured channel".
13. **Δ_g = 100 failing spins.** UT §6 "near-gap J ≳ 18"; DG §2.1 "J ≥ 16 (rows held); 18, 24, 28 (Γ ∝ m)"; RD-S §5.4 "J ≥ 18–24".
   Draft: J ≳ 18 (supplement only; letter says "near-gap large spins").
14. **D per-state records.** RP D4: RD-S's J = 15 (1000) value +0.869 is not a minimum over the physical mass range (hybrid ≈ +0.69 at
   the twist edge); RP D5: J = 16 (300) +0.884 is the near-gap minimum only (hybrid ++ +0.805 at 3.17Λ). Neither is binding; both are
   flagged in the supplement.
15. **BCPW author list** differs between LITERATURE.md and ABJM_THEOREM.md (see §5).
16. **DG §1.2 vs RD-S δk(0) at 100** (RP D7) and **c₁ at 100** (∓2.37 %): irrelevant here, nothing at 100 is quoted.
17. **Corollary below 300.** DG §2.4 gives a model bracket `140 ≲ Δ_g* ≲ 300` [M]. Deliberately omitted from both documents, because
   quoting it would read as a claim below 300; the draft only says the threshold was not computed.
18. **ε_c₁ coefficient.** UT §2: 2.03; DG: 2.032; RQ-S §0: "2|c₁t₄|Λ²∫ψ₂ with ∫ψ₂ = 1.016". Consistent; letter uses 2.03, supplement 2.032.
19. **Unclear:** none of the six update documents is marked reviewed; REPRO and REPRO_ENGINE are independent reproductions but not
   referee reports. The draft marks items resting on HRS_CLOSE, RD_CLOSE, EINSTEIN_CHECK and DG_SCAN with † ("not yet reviewed").
