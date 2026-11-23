# Third referee report

**Manuscript:** "Conditional bounds on R³ and one R⁴ graviton coupling in AdS₄/CFT₃ from ⟨TTTT⟩ light-ray sum rules", v5
(`ads4/prl/letter.tex`, `supplement.tex`, `macros.tex`, `refs_verified.bib`, `fig/`, build in `prl/build/` of 23 Nov 17:01;
switches `\HLJclosedfalse` (headline Δ_g = 1000) and `\HSdraftfalse`).

*Written 23 Nov 2026 as an internal pre-submission review. No existing file was edited and nothing was committed. References are
`file:line` in the current text, relative to `ads4/prl/` unless they say otherwise. I formed the view in §§1–3 and §6 before reading
`REFEREE_REPORT.md` and `REFEREE_REPORT_2.md`, then checked report 2 item by item (§4). Lines marked `% P2FOLLOW` were judged on
whether the claims stand whatever the running P2 follow-up returns.*

**What I checked myself** (scratch scripts in `/tmp/claude-0/-home-user-MAJOR/60eea3a7-07d6-5453-af55-c4346d86f9b5/scratchpad/`):
* **Headline constants.** I recomputed them from `data/nf02_psi.json` and `data/q3star_factored.json` in 40-digit mpmath (`chk1.py`).
* **Magic-block identity.** I reran it at fresh random complex ν (`chk_magic.py`):
  * the exact sympy check, n ≤ 4, b = 0…4;
  * the 50-digit check, n ≤ 30;
  * a control with wrong coefficients;
  * a cross-check that the script's recursion equals the production block code `eestep/ee_blocks.block_coeffs_b`.

  I imported the project modules read-only, with no bytecode, and wrote only to the scratchpad.
* **Positivity, large spin and window numbers.** I read the J = 29, binding-state, Theorem W/L and monopole numbers in the stored
  outputs (`code/highspin/cache/hs_odd_J29_final.json`, `code/neargap2/cache/ng_bind_J5.txt`,
  `code/pwimg/cache/pw{1,1b,3,4,5}_*.txt`). I traced the window numbers to `SMALLNU_WINDOW.md` §0.
* **Length.** I counted it from the source and from the rendered PDF.

---

## 1. Summary and recommendation

**What the Letter claims.**
* **The implication.** The implication is graded [D-sketch]. Assume:
  * a higher-spin gap Δ_g in T × T (A1′) and tree level (A2);
  * a list of analytic inputs;
  * per-state positivity (Pos) of fixed light-ray functionals on all heavy states and on light spin-0, 2 states.

  Then three bounds follow.
  * **(1)** τ = √(t₄²+α_T²) ≤ 2.89·10⁻⁸, i.e. |ĝ₃|Λ⁴ ≤ 120, at Δ_g = 1000.
  * **(2a)** g̃₄Λ⁶/8πG ≤ 92.8, if no scalar lighter than 0.46Λ couples to TT.
  * **(2b)** g̃₄Λ⁶/8πG ≤ 164.5 without that condition.

  Each holds exactly for every value of an uncomputed universal constant c₁.
* **Positivity.** Positivity is certified in the flat limit for Ω^EE. In AdS it is supported by error-budgeted numerics on a finite set of
  states.
* **Large spin (new in v5).** It is "reduced to a conditional theorem": the magic-block identity, Theorem W (J\* = 50) and Theorem L
  (J ≥ 857).

**Assessment.** The Letter is now honest in most places, and much better organised than v3.
* Every constant I traced reproduces.
* The J = 29 problem of report 2 is resolved without a negative state, and the main text says how.
* The magic-block identity for the b-shifted EE blocks is a clean, checkable, exact result. I reproduced it.

But the physical content is still positivity of a functional on an infinite spectrum, and it is still not established.
* **Near the gap at Δ_g = 1000, only ten spins are error-budgeted:** J = 0, 2, 4, 5, 6, 8, 12, 16, 21, 29. J = 8 is error-budgeted only
  on a less conservative budget, and J = 21 and 29 are thin (leverage 1.2).
* **Near-gap 7 ≤ J ≤ 48 otherwise** is model or interpolation.
* **Near-gap 50 ≤ J ≤ 344** is a conditional statement. Its key premise (Λ_amp ≤ 3) is measured on one state, (1100, 100).
* **Near-gap 345 ≤ J < 857** (twist-edge masses 1.3–1.65Λ) is covered by neither theorem, and the record calls it "the hardest band".
* **J ≥ 857** is conditional on the same unknown spinning coefficients.

**A new physics concern (M2).** The interchange H_K, the one analytic input on which every bound rests, is proved under the
non-perturbative Regge bound J₀ ≤ 1. The bounds, however, are tree-level statements, where the planar intercept is
2 − O(Δ_g⁻²) and the authors' own k = 2, 3 criterion (J₀ < 3/2) fails.

**Beyond positivity:**
* the implication is [D-sketch];
* c₁ is uncomputed;
* the bounds hold at two gap values only;
* they are about 10× (|ĝ₃|) or 10² (|ĝ₃|²) weaker than the transplanted flat optimum;
* the R⁴ window is IR-log dominated;
* the record behind v5 (PW_IMAGES, LARGESPIN_W/L, SMALLNU_WINDOW, DIRCV) was written the same day and is unreviewed.

**PRL criteria.**
* *Importance.* Real, for specialists. Explicit-constant graviton bounds from ⟨TTTT⟩ in AdS are an open target named by CLS-D, and the
  EE structure, de-shift and k = 4 arc are genuine technical ideas.
* *Broad interest.* Limited. The result is conditional, at a single gap value, with loose constants, and the scaling was known (CEMZ,
  AKT). About 28 % of the main text is a positivity ledger (§3, m1).
* *Claims established at the level stated.* The conditional implication is established at [D-sketch], modulo M2. The numbers that make
  the paper interesting are not established, because they rest on (Pos).

**Would an editor send it out?** Perhaps; it is close to a coin flip. The topic is current and the abstract is careful. Against it, the
abstract is hedged three times, the result is at one gap value, and both authors list "Independent researcher" with a 20-page
supplement of internal grades. Some editors will judge it "specialised and conditional" and decline without review.

**Would a strong referee accept it?** No. A referee from the CHMRS/CLS-D community would:
* ask for positivity across the full spectrum, at least at the [NE±] level with an analytic large-spin argument, as CHMRS supply for
  scalars;
* ask how the tree-level sum rule is justified (M2);
* ask for c₁;
* recommend a long paper.

**Recommendation: reject for PRL in the present form; transfer to PRD (or a JHEP long paper) strongly encouraged.** This is unchanged in
substance from report 2. v5 improves the evidence (J = 29, J = 21, the window, the large-spin structure) but does not change its
category. A PRL resubmission becomes realistic if four things happen:
1. M2 is settled.
2. The spinning O(1/m²) coefficients are computed, so that Theorems W and L have derived rather than measured premises.
3. The near-gap band 7 ≤ J < 857 is error-budgeted.
4. The main text is cut back to the physics.

---

## 2. Major issues

### M1. Large spin: the "conditional theorems" cover less than the Letter suggests, and the Letter omits premises

**Where:** letter.tex:78-80 (Introduction), :272-276, :281-289 (Status, (ii)–(iii)), :300-317 (Large spin); supplement.tex:1002-1124.

**Claim:** "at large spin it is reduced to a conditional theorem" (:79-80). "Theorem W: near-gap states with m ≤ 1.3Λ and J ≥ J\* = 50
(J\* = 100 with a margin of 2) give a positive semidefinite matrix for every coupling ratio, provided …" (:308-312). "Theorem L: for
J ≳ 0.52m … absorbs any O(1/m²) amplitude up to 43 times the scalar one for J ≤ 10⁶" (:312-315).

**Problems.**

1. **The theorems leave a gap.** For near-gap states, m_gap(J)/Λ = √(1+2J/(Δ_g−1)) exceeds 1.3 at J ≥ 345 (my arithmetic, and Table S11
   of the authors).
   * So "m ≤ 1.3Λ, J ≥ 50" is exactly 50 ≤ J ≤ 344, and Theorem L starts at J = 857.
   * Near-gap 345 ≤ J < 857 is covered by neither. It rests on H_R′ [D-sketch + M], with a leading-Regge floor that is only [NE] for
     m < 1.5Λ, b ∈ (119, 1.04Λ) (supplement.tex:995-996, :1116).
   * LARGESPIN_L §0 calls 150 ≤ J < 857 near the gap "the hardest band, still numeric/model". The de-shift tolerance has its minimum,
     6.7·10⁻⁴, at J ≈ 350–500.
   * The authors know this (NOTES.md §7 item 6), but the Letter never says it. Item (ii) (:281) says "larger near-gap spins (even J ≥ 50,
     odd J ≥ 33) rest on … H_LJ(1000)", which a reader takes to mean all larger spins.
2. **Theorem W is a statement about 40 sampled states, and P0 is omitted.**
   * The record states it "on every sampled H_LJ(1000) state with J ≥ J\*" (LARGESPIN_W §0 TW; PW_IMAGES TW′). The samples are even
     J = 26–344 and odd J = 33–301, at three masses. Between samples it is [NE], and odd 303–343 are extrapolated.
   * Premise P0, the leading-Regge floor, is not proved in Theorem W's region. The extended Mehler–Dirichlet bound covers m ≥ 1.3Λ only
     (LARGESPIN_ANALYTIC §2; LARGESPIN_L L1′: "m ∈ [Λ, 1.3Λ), b̂ < 119: [NE]"). The Letter does not list P0 at all.
   * The supplement grades P0 "curvature [P + DV] for b ≤ 119" (supplement.tex:1006-1007), which is wrong for m < 1.3Λ.
3. **The Letter's premise list for Theorem W is incomplete.**
   * J\* = 50 holds for (Λ_amp, K_sh, K₄) = (3, 10, 10). The Letter states only Λ_amp ≤ 3.
   * The index-shifted bound K_sh ≤ 10 (absolute, in c₂ units) and the value K₄ ≤ 10 are missing.
   * (P4) has an open sub-point [I]: the flat ν⁶ term of the O(1/m⁴) series must be subtracted. "If it were costed as a Regge term over
     0 ≤ ν ≤ Λ it would exceed every tolerance" (supplement.tex:1056-1060; NOTES §7 item 7).
   * "The O(1/m⁴) part is bounded" (:311) hides an [I].
4. **Λ_amp ≈ 2–3 comes from one state.** It comes from the class-B fits at (1100, 100), with the (364, 64) fit at Δ_g = 300 "bracketing
   the model (inconclusive)". "measured as ≈ 2–3 times the scalar ones only at J ≤ 100 [NE]" (:316) reads as a range of spins. The DIRCV
   fits that produce it have 2σ leverage down to 0.5–0.6 in one fit.
5. **Theorem L's tolerance is for a specific shape.** "Any O(1/m²) amplitude up to 43 times the scalar one" is the tolerance for the
   **scalar-shaped** vector (CHMRS (3.20) coefficients times one common factor, worst radial multiplier). It also assumes index shifts
   "normalised as the light-cone limit requires" (supplement.tex:1072-1074).
   * Without that normalisation, a shift by two is tolerated only at K ≈ 0.4 c₂ at J = 10⁶ (`pwimg/cache/pw5_theoremL.txt:45`).
   * The tolerance falls with J (1240 → 43 from J = 857 to 10⁶), because the monopole cost is J-independent while the floor decays like
     βe^{−β/2}. **Nothing is said for J > 10⁶.**
6. **Neither theorem reduces the uncertainty in kind.** Both convert "H_LJ / H_R′ hold" into "the spinning O(1/m²) coefficients are at
   most a few times the scalar ones". That is useful structure, but the premise is the same unknown physics.

**Why it matters.** The new "Large spin" paragraph is the main v5 addition, and its framing will be read as "large spin is handled up to
a mild premise". The record (LARGESPIN_W: "Do not write 'proved for J ≥ 50'"; PW_IMAGES: "Neither becomes unconditional") is more
careful than the Letter.

**Fix.** Rewrite :300-317 in about 150 words. The paragraph should:
* name the three regions near the gap: 50 ≤ J ≤ 344 (Theorem W, on sampled states); 345 ≤ J < 857 (no theorem, H_R′ [D-sketch + M]);
  J ≥ 857 (Theorem L, up to J = 10⁶);
* list all premises, P0 (computed, [NE] near the gap) and the K_sh, K₄ values among them, plus the [I] subtraction in P4;
* say that Λ_amp is measured at one state;
* replace "any O(1/m²) amplitude" with "a scalar-shaped O(1/m²) correction".

In the Introduction (:79-80), replace "at large spin it is reduced to a conditional theorem" by "on parts of the large-spin region it is
reduced to conditional statements whose premises are measured on one state or modelled".

**Fatal for PRL?** In combination with M5–M7, yes. As a wording issue alone, no.

### M2. The Regge input: a non-perturbative J₀ ≤ 1 is used for a tree-level interchange

**Where:** letter.tex:133-135, :183-185, :109; supplement.tex:217, :225, :298-306, :475-543; UT3 §4.1 (R1, H_K).

**Claim.** "its Regge spin k must exceed the CFT's Regge intercept J₀ ≤ 1 [KKSDZ, Caron-Huot]" (:134-135). "The interchange of the OPE
sum and the integrals is absolutely convergent for J₀ < 3/2 (k = 2, 3, 5) and J₀ < 5/2 (k = 4) [D], with margin under J₀ ≤ 1" (:183-185).

**Problem.**
* **The bound J₀ ≤ 1 is non-perturbative**, a property of the full correlator at finite c_T.
* **The Result is a tree-level statement** (A2, :109; "bounds refer to c_T → ∞ at fixed Δ_g", UT3 §1). Its light side is the tree-level
  graviton exchange, and its heavy sum is the O(1/c_T) single-trace sum.
* **At O(1/c_T) the relevant intercept is the planar one.**
  * It is J₀ = 2 for Einstein exchange; the record itself says "J₀ = 2 < 4" for the Einstein tree correlator (LEMMAS_CLOSE §1.4).
  * It is J₀ = 2 − O(Δ_g⁻²) after stringy Reggeization.
  * That is why the graviton is a borderline "arc" at k = 2 (UT3 L1: "H7 = the Regge arc at n = k = 2").
* **The authors' own corner count fails at that intercept.** The count requires c < a + 1 − J₀ with c = −½, a = 0 for k = 2, 3
  (supplement.tex:487-490, :515). That is J₀ < 3/2, which J₀_planar ≈ 2 violates.
* **The order of limits is not discussed anywhere in the Letter or Supplement.** One could apply the sum rule at finite c_T and expand in
  1/c_T, but that needs uniformity of the heavy sum in c_T. The Regge tail that the large-c_T expansion fails to control is exactly what
  J₀_planar measures.
* **The k = 4 criterion is not at risk.** J₀ < 5/2 is above 2.

**Why it matters.**
* H_K is the input that turns the per-state positivity of a fixed functional into a bound. It is graded [D] with "margins ½ and 3/2".
* If its tree-level version needs J₀_planar < 3/2, it fails for every stringy holographic CFT, which are the theories the Letter is
  about (ABJM included). Item (1) and the k = 2, 3 parts of Ω^Q_up would then lose their derivation.
* The flat-space analogue (CHLPS) works with spin-2 subtracted sum rules that converge for Regge growth < s², and the graviton sits on
  the boundary. CHMRS and CLS-D state a planar Regge-boundedness assumption of this kind (J₀ < 2 at tree level); please check and cite the
  exact statements. The margin at k = 2 is then O(Δ_g⁻²), not ½.

**Fix.** Do one of the following:
* **(a)** Show that the cut-off EE/Q weights give a corner exponent c < a + 1 − 2 at k = 2, 3, i.e. c < −1, so that absolute convergence
  holds at J₀_planar < 2.
* **(b)** Replace absolute convergence by an argument that needs only k > J₀_planar.
* **(c)** State the tree-level sum rule as an explicit assumption [A], as CHMRS and CLS-D effectively do. Then grade H_K "[D]
  non-perturbatively; tree-level use [A]", and change the sentences at :134-135 and :183-185.

Either way, say which J₀ is meant, and whether KKSDZ §4.2.6's caveat or their large-N discussion covers it.

**Fatal?** Not necessarily: option (c) is a one-paragraph change. But until it is addressed the grade of the implication is overstated,
and a CHMRS/CLS-D referee will ask this first.

### M3. The new [NE±] states are thinner than the main text conveys

**Where:** letter.tex:257-271; supplement.tex:869-895.

**Problems.**
* **The ±0.018 on J = 29 is Monte-Carlo only.** "Direct evaluation gives +0.972 ± 0.018" (:267) quotes the jackknife MC error
  (`hs_odd_J29_final.json`: `main` 0.97159, `jk` 0.0183).
  * The error budget is 0.70 (HIGHSPIN λ-level) + 0.09 (window) = 0.79, i.e. L = 1.23.
  * The interpolation variants in the same file range from 0.67 (`K2/eshape`) to 1.44 (`dir/espl`).
  * A reader will take ±0.018 as the uncertainty of the sign. It is not.
* **"Now measured" overstates.** "The low-frequency window (ν < 50) … is now measured on these states" (:269-270). Below the lowest
  measured node (ν < 10 at J = 29, ν < 20 at J = 21), the old 10⁴ c₂ convention is kept (SMALLNU_WINDOW §0). At J = 29 it carries 0.052
  of the 0.093 window cost.
* **J = 29's grade is borderline.** J = 29 is [NE±] only with the signed window, where the sign-free ×3 envelope gives L = 0.94 (stated
  at :270-271, which is good). By the definition at :82-83 ("×3 safety factors on all model terms"), J = 29 at the twist edge is
  borderline [NE±]/[NE + M]. The grade should carry "(signed window)".
* **The window becomes live again at higher spin.** SMALLNU_WINDOW §0 item 5: "For J ≳ 32 the window is again a live input even with
  measurements". That applies to every interpolated and H_LJ spin, and the Letter does not say it.

**Fix.**
* Write "+0.972 (MC ± 0.018; error budget 0.79, leverage 1.2 with the signed window)".
* Write "measured on ν ≥ 10–20; conventional below".
* Add one clause on J ≳ 32.

**Fatal?** No.

### M4. Framing in the abstract, the Introduction and "Why it matters" is stronger than the record

**Where:** letter.tex:37, :78, :88-98.

**Problems.**
* **(a) "proved in the flat limit" (abstract :37; Intro :78).**
  * For Ω^Q_up, the functional of item (2), flat positivity is certified on four regions. It is [NE] on two others (float minima 0.997
    and 0.304; UT3 §5.1; letter :254 itself says "all but two asymptotic regions").
  * On light scalars below 0.46Λ, Ω^Q_up is *proved negative*.
  * The abstract's parenthesis is true for Ω^EE only.
* **(b) "every holographic CFT₃ with Δ_g = 1000 has stress-tensor three-point couplings within τ ≤ 2.89·10⁻⁸ … less than 10⁻⁸ of the
  range τ ≤ 4" (:88-92).**
  * The bound is the O(c_T⁰) (tree-level) value. Loop corrections to t₄ are O(1/c_T): in the free-field basis c_T t₄/4 = n_B − n_F
    (DICTIONARY_T4 §0), so light-field loops generically give δt₄ ~ #/c_T.
  * For the tree-level statement to bound the actual t₄ at the 10⁻⁸ level one needs c_T ≳ 10⁸ at Δ_g = 1000. A2's "c_T ≫ Δ_g²" = 10⁶ is
    not enough.
  * The sentence should say "at leading order in 1/c_T".
* **(c) The "10⁻⁸" is mostly the known Λ⁻⁴ scaling.** The fraction τ₀/4 = 7·10⁻⁹ is dominated by Λ⁻⁴ ≈ 10⁻¹², i.e. by the
  CEMZ/AKT parametric scaling, which the Letter itself credits (:53-55).
  * The new content is the O(10²) coefficient |ĝ₃|Λ⁴ ≤ 120.
  * Comparing with the collider disk advertises the scaling as if it were the result.
  * State the bound in gap units in "Why it matters", as the Discussion does (:338-342).
* **(d) The de Sitter / inflation sentence (:356-358)** is speculative and unreferenced. Drop it, or cite the dS positivity literature
  and say why the AdS mechanism (de-shift against AdS curvature) would carry over.

**Fatal?** No; wording.

### M5. Importance for PRL: one gap value, loose constants, and an uncomputed constant

**Where:** abstract; letter.tex:248-249, :338-356; Fig. 1.

**Problems.**
* **One gap value.**
  * The headline holds at Δ_g = 1000 only. Δ_g = 300 has an open [I] hypothesis, and "other gaps are not claimed" (:249).
  * A bound that exists at one value of the gap is a data point, not a constraint on a family of theories.
  * The Λ-forms are values of the implication; Fig. 1 draws them faded.
* **Loose constants.**
  * Item (1) is ≈ 10× weaker in |ĝ₃| than the transplanted flat optimum, and ≈ 10² weaker in |ĝ₃|² (stated, :342-345).
  * The R⁴ upper end is dominated by the 4d IR log: 17.24 ln(2Λ/3) = 112 of 92.8.
  * The one theory in scope, ABJM-IIA, sits at 2.6 %, so the check can only fail grossly (stated).
* **An uncomputed constant.**
  * c₁ is a universal number fixed by a tree-level ⟨TTTT⟩ computation (L3d). It is not computed, and c₁ = 0 rests on an assumption F₂
    [A].
  * Since it is universal, a referee will ask why the one Witten-diagram calculation that fixes it (graviton exchange with one W³
    vertex) was not done before the bounds were advertised.
* **A sketch-level implication.** The implication is [D-sketch], with weakest links in the c₁ form, the light single traces and the
  ε-repair.

**Fix.**
* Extend positivity to a range of gaps, or prove that it is monotone in Δ_g. This is the single change that would most raise breadth.
* Compute c₁.
* Otherwise accept that this is a long-paper result.

**Fatal for PRL?** Yes, in combination with M1 and M6–M7.

### M6. The v5 material is unreviewed, same-day work, and the Supplement still sketches the load-bearing steps

**Where:** supplement.tex:39-75 (provenance), :929-1125; CHANGES_v5 §7 item 3.

**Problem.**
* PW_IMAGES, LARGESPIN_W/L, SMALLNU_WINDOW, DIRCV §8 and HIGHSPIN are all dated 23 Nov and unreviewed. The P2 follow-up that fixes
  Theorem W's K₄ and (S) premises was still running when the Letter was written.
* INDEP_CHECK independently reproduced the EE reduction and the Lemma Y count, which is real progress on report 2's M3. Not reproduced:
  * the L6′ correlator-to-light-function step;
  * K1′, K2ν and LRC-a;
  * the k = 4 kernel normalisation;
  * everything in v5's large-spin section.
* I reproduced the magic-block identity (§5). It is the easiest piece; the costs, floors and tolerances behind J\* are not
  independently checked.

**Fix.**
* An independent reproduction of Lemma LR∞ (with interval constants) and of the Theorem W costs at a few J.
* A written-out long-paper derivation of H_K and the c₁ structure.

**Fatal?** Not by itself.

### M7. Readability has regressed: the positivity ledger is back in the main text

**Where:** letter.tex:252-320.

**Problem.**
* v4 moved the ledger to the Supplement. v5 put back about 830 words, by my count, about 28 % of the prose.
* Those words are mostly status bookkeeping:
  * "Status of positivity" (566 words) and "Large spin" (264 words);
  * leverage values to two digits;
  * Monte-Carlo noise, "minimal ring class", "bias correction", "sign-free envelope";
  * the monopole, the "magic-block relation", P0–P4 premises in prose;
  * three lists of spins.
* A general reader cannot use this, and a specialist wants it in a table.
* The figure plots two curves that are essentially ln Λ. It does not show the most informative object, the coverage of (Pos) in the
  (J, m/Λ) plane by grade.

**Fix.**
* Replace :252-320 by about 250 words: what is certified, what is error-budgeted, what is conditional, and what is open.
* Replace or augment Fig. 1 with a (J, m/Λ) coverage map at Δ_g = 1000, coloured by grade, with the two theorems' regions and the
  uncovered band 345 ≤ J < 857 marked.
* Move every leverage value to Tables S5, S8 and S11.

**Fatal for PRL?** In combination with M1 and M5, yes.

---

## 3. Minor issues

**Length.**
* **m1. Length is compliant. Density is the problem, not the count.**
  * *My count.* About 2970 prose tokens if every math token counts (detex of :43-359), and about 1825 if each inline formula counts as
    one word. To that add:
    * the caption, ≈ 80–120;
    * three display rows, 48;
    * the figure, single column at aspect 0.9375, 150/0.9375 + 20 = 180.
  * *APS-style total:* ≈ 3100–3300, consistent with the authors' 3110–3210 and under 3750.
  * *Rendered length.* The rendered body is 3576 PDF tokens. The text ends ≈ 3.65 pages. The abstract has 563 rendered characters
    (limit 600).
  * *Verdict.* Within limits with about 450–600 words of headroom. M7's cut would add headroom.

**Numbers and wording.**
* **m2. B_eff depends on Λ (:226-231).** The text writes B_eff = −19.265 and quotes 72.04 at Δ_g = 300.
  * With the published weights, B_eff = −19.2646 at 1000 but −19.2197 at 300 (my `chk1.py`, matching report 2).
  * 17.2407 ln(2Λ/3) − 19.265 gives 71.995 at 300, not 72.04.
  * Write "B_eff(Λ) = −19.265 (1000), −19.220 (300)", or keep only the exact 2L_AdS.
* **m3. Overlapping spin lists (:277-280).** J = 26, 36 and 40 appear both as "computed inside explicit models" and inside "even J = 26,
  30, 34–46 … interpolated". HIGHSPIN lists J = 26 as not run; the underlying UT3 rows (§5.2) carry the same clash. Make the lists
  disjoint.
* **m4. Item (ii) needs the mass range (:281).** Write "near-gap 50 ≤ J ≤ 344 (m ≤ 1.3Λ)". The H_LJ region is m ≤ 1.3Λ, and twist-edge
  states leave it at J = 345 (M1.1).
* **m5. "proved inequalities with floating-point constants" (:274) is self-contradictory.** In Lemma LR∞ the constants are evaluated in
  floating point on a β grid (step 0.01–0.5 to β = 64), and the tail is [D] (supplement.tex:984-988). By the paper's weakest-link rule,
  "for all β ≥ 1" is [D] (or [P + DV] on the grid), not [P].
  * The same applies to "which we prove positive" (:313).
  * Write "[P + DV] (inequalities proved, constants on a grid; tail [D])".
* **m6. "K = 1 + O(Λ⁻²) [DV]" (:188).** The O(Λ⁻²) H17 band is [D-sketch + M] (UT3 L2; supplement.tex:310). Write "δK₀ [DV]; H17 band
  [D-sketch + M]". The effect is a factor 1.00003 at Δ_g = 1000, so no number changes.
* **m7. The "So" in "So every beyond-leading-Regge correction is exactly the de-shift of a scalar image" (:305-306)** rests on the
  identity at all orders [D] and on the EESTEP representation [D + DV]. Attach "[D]".
* **m8. Theorem L's "any O(1/m²) amplitude" (:314).** Write "a scalar-shaped O(1/m²) correction (index shifts with light-cone
  normalisation)" (M1.5).
* **m9. "only at J ≤ 100" (:316).** Write "at one state, (Δ, J) = (1100, 100)" (M1.4).
* **m10. Theorem W's J\* depends on all three constants (:308-312).** State "for (Λ_amp, K_sh, K₄) = (3, 10, 10)". Mention the open
  P4 subtraction, or keep J\* out of the Letter until P2FOLLOW finishes (M1.3). As an implication Theorem W stands whatever P2 returns;
  the number J\* = 50 may not.
* **m11. The c₁ threshold for item 2 (:238).** "|c₁| ≳ 1260 (item 2)" is the (2b) value; for (2a) it is 1306 (UT3 §3). Both are
  conservative, but name the item.
* **m12. Give the c₁ coefficient (:218).** "x ∝ |c₁|/Λ²": give x = 19941|c₁|/(Λ²√(16KL_A)). It is short and exact.
* **m13. Composite grades undefined in the Letter.** The Letter uses "[P + DV; certified]", "[implication P; costs DV …]",
  "[D given these magnitudes]" and "[D-sketch + NE]"-type composites. The legend (:80-85) does not define "+", "given" or
  "implication". Add one clause, or point to §S1.
* **m14. Fig. 1.**
  * The annotation "√(t₄²+α_T²) ≤ 2.89·10⁻⁸" should carry "c₁ = 0, K = 1" (the caption does).
  * The "not claimed" shading is only below 300, although everything except the two markers is unclaimed (:248-249). Shade or hatch
    the whole line, or label it "implication only".
  * The 300 markers deserve a symbol that says "one hypothesis open".
* **m15. Two scope points belong in the Letter.**
  * The [NE±] strip m ≥ 1.3Λ, b < 5 (:259) depends on one model input, the bulk-ν envelope H_bulk [M], tolerated ×76. Say so in a
    clause.
  * The mass band 1.1Λ < m < 1.3Λ at J ≤ 48 is covered only by H_cont [NE] (Table S7). It appears nowhere in the Letter's weak points.
* **m16. Cite DKKPS (arXiv:1708.05718) at :111-112** for the absence of ⟨TTO_J⟩ structures at J = 1, 3. Note that odd J ≥ 5 occur only
  parity-odd; that is why the binding state is parity-odd.
* **m17. Light-scalar condition (2a).** The binary condition "no scalar … appears" is stronger than needed. The Supplement's weighted
  form (upper end + Σ max(0, 0.4304 − 2.0321 m²/Λ²)|g₀|²Λ⁴/8πG; UT3 §3) makes item (2a) usable when light scalars couple weakly through
  φW², as they do generically. One sentence in the Letter would turn a condition into a bound.
* **m18. ABJM (:328-329).** "only the stress-tensor multiplet is exchanged in T × T at tree level" holds through the orders computed
  (supergravity + R⁴; UT3 §8). Say so; at higher orders in α′, KK modes can couple to TT through φW²-type terms.
* **m19. HKT is mis-cited for commutativity (:56).** HKT proves ANEC; commutativity of ANEC operators is KKSDZ (and Kologlu et al.'s
  light-ray OPE). Cite HKT only for ANEC.
* **m20. A trivial overstatement in the record (not in the Letter).** The identity residual "≤ 8·10⁻⁴⁸ (also b = 7, 10)" is ≤ 7.9·10⁻⁴⁸
  for b = 0…4; for b = 7 it is 2.0·10⁻⁴⁷ (`pw1_magic.txt`).
* **m21. End-to-end accuracy of the magic identity (supplement.tex:960).** "2·10⁻⁷–1.1·10⁻⁶" is the refined run (`pw1b_e2e.txt`). The
  default run gives up to 7.6·10⁻³ at ν = 6 (`pw1_magic.txt` (c)). Say "with refined quadrature".
* **m22. Recheck the supplement's P0 grade (supplement.tex:1006-1007)** for m < 1.3Λ (M1.2).

---

## 4. Report 2: status of each item in the current text

| report-2 item | status | evidence (current v5 text) |
|---|---|---|
| **M1** positivity not established; main text understates | **partly** | J = 29 resolved and explained in the main text (:264-268); J = 12, 16, 21, 29 [NE±] (:260-262); canonical-gauge failure disclosed (:290-291); weak points (i)–(iv) listed. Not closed: J = 7–48 mostly [NE + M]/[M]; 345 ≤ J < 857 is in no theorem (M1); the new large-spin theorems rest on one measured state; J = 21, 29 thin (M3). Abstract: "no computed state is negative", but the modelled regions are not named (report 2's fix 2 only partly done). |
| M1(b) J = 29 in the main text | **addressed** | :264-271 |
| M1(c) gauge / in-sample selection | **addressed (disclosed)** | :290-291 |
| M1(d) fallback for a negative HIGHSPIN state | **addressed** | NOTES.md §8 (not triggered; trigger list updated) |
| **M2** lower end advertised without evidence | **addressed** | Dropped from the abstract and Discussion; "A lower bound on g̃₄ is not claimed" (:245-247); outlook in §S4.A |
| **M3** unreviewed derivations, only sketched | **partly** | INDEP_CHECK reproduces the EE-kernel reduction and the Lemma Y count from the statements (:243-245; §S10). Still sketched: H_K corner count, L6′ correlator step, R⁴ calibration, c₁ ASD closure. v5 adds a large unreviewed same-day layer (M6) |
| M4(a) c₁ universality [D] | **addressed** | [D-sketch] (:190) |
| M4(b) arc qualifier | **addressed** | Not stated separately in the Letter; Table S2 L4 keeps the qualifier |
| M4(c) L6′ numerics pre-asymptotic | **addressed** | :191-195 |
| M4(d) H_EFT(a) −3 % | **addressed** | Supplement L7, I5; the Letter no longer claims it |
| M4(e) "bound" for [NE±] | **addressed** | "error budget … not a bound" (:82-83, :258) |
| M4(f) "Theorem" label | **addressed for the Result** | "Result (conditional; … [D-sketch])" (:207). New: "Theorem W/L" (:308, :312) are honest as implications but see M1 |
| **M5** positioning, strength, CHLPSpw | **addressed** | CHLPSpw correctly described, CHLPS CFT₃ remark named (:60-65); "≈ 10× in \|ĝ₃\|, ≈ 10² in \|ĝ₃\|²" (:342-344); IR-log dominance (:345-346); ABJM one-digit (:333-334). New overstatement in "Why it matters" (M4(b, c)) |
| **M6** not written for a broad audience | **partly; regressed in v5** | v4's Mechanism paragraph and figure kept (:130-161, Fig. 1). v5 re-adds ≈ 830 words of ledger (M7); jargon is back ("leverage", "twist edge", "sign-free envelope", "ring class", "monopole", "bias correction"). No coverage figure |
| **M7** (2a) condition; repair couples the items | **addressed (physics open)** | Abstract precise (:39-40); (2a)/(2b) dependence stated (:238-239); re-optimised Q still a next step (:356) |
| **M8** corollary and "explicit dependence on Λ" | **addressed** | Removed from the abstract; "other gaps are not claimed" (:248-249); remark in §S8 |
| m1 length | **addressed** | ≈ 3100–3300 APS words; text ends ≈ 3.65 pp (§3 m1) |
| m2 front matter | **addressed** | Affiliations and e-mails (:22-28); no acknowledgments (acceptable) |
| m3 JHEP years | **addressed** | Rendered "J. High Energy Phys. 2009 (10), 079" (build/letter.bbl) |
| m4 abstract wording | **mostly addressed** | g̃₄ "defined by the CFT", light states, 0.46Λ. New: "proved in the flat limit" overstates for Ω^Q_up (M4(a)) |
| m5 notation | **mostly addressed** | R → τ; K₀, r_±, ϱ, c_Q in the Supplement. k (Regge spin) and K (coefficient) remain; harmless |
| m6 grade legend | **partly** | Supplement §S1 defines the composites; the Letter uses new undefined ones (m13) |
| m7 Table I | **moot** | Table removed; rows in Tables S5/S7 |
| m8 "≥ 19×" | **addressed** | "at least 18×" (:295) |
| m9 A1′ wording | **addressed** | :111-112 |
| m10 parity violation | **addressed (Supplement)** | §S2.C paragraph (per CHANGES_v4) |
| m11 natural size of c₁ | **addressed** | Withdrawn |
| m12 ABJM λ convention / citations | **addressed** | supplement.tex:1231-1232; CPY alongside BCP; CHMRSswamp equations cited |
| m13 missing references | **addressed** | CMT (:203), HKT (:56), Häring–Zhiboedov (Supplement). New misuse of HKT (m19) |
| m14 truncated quartic form | **addressed** | Not in the Letter |
| m15 corner code unreviewed | **addressed** | :255-256 |
| m16 corollary J = 0 grade | **addressed** | Supplement |
| m17 "plausibly" | **addressed** | Supplement Table S6 |
| m18 Regge caveat | **addressed as asked** | supplement.tex:217, :298-301. But see M2: the more important tree-level vs non-perturbative intercept question was not raised before and is not addressed |
| m19 data hygiene | **addressed** | `data/README.md` §3, §5 errata; `q3star_factored.json` published |
| Report 2 §6 item 2 (close (Pos) at 1000) | **partly** | See M1, M3 |
| Report 2 §6 item 6 (compute c₁ or derive F₂) | **not addressed** | c₁ still [I] without F₂ |
| Report 2 §6 item 7 (re-optimise Q) | **not addressed** | Listed as a next step |

**Report 1:** report 2's §4 table remains accurate for report 1's items, with these updates.
* R1-M4 (lower end) is now addressed.
* R1-M11 (corollary) is addressed: it is a Supplement remark.
* R1-M15 (presentation) has regressed again in v5 (M7).
* R1-M3's explicit spinning R⁴ ⟨TTTT⟩ contact is still not evaluated (L6′ "not evaluated", supplement Table S2).

---

## 5. Spot checks of load-bearing numbers

| # | number (Letter) | traced to | my check | verdict |
|---|---|---|---|---|
| 1 | τ₀ = 2.89·10⁻⁸ / 3.13·10⁻⁶; \|ĝ₃\|Λ⁴ ≤ 120 / 104 (:218-220) | UT3 §3; `nf02_psi.json` | Recomputed: L_A = 896.224 / 670.913; √(16L_A) = 119.748 / 103.608; τ₀ = 2.891263·10⁻⁸ / 3.132055·10⁻⁶; k₀ = 27.48509817 | ✓ |
| 2 | 2L_AdS = 92.81 / 72.04; B_eff = −19.265 (:226-231) | `q3star_factored.json` | 92.81351 / 72.04059; B_eff = −19.2646 (1000) but −19.2197 (300) | ✓ for 2L_AdS; **B_eff is Λ-dependent** (m2) |
| 3 | c₁² coefficients 5.4·10⁻⁷ / 6.8·10⁻⁵; x/\|c₁\| = 1.67·10⁻⁴ (:219, :230-231) | UT3 §3 | 5.946·10⁴/(κΛ⁴) = 5.4426·10⁻⁷ / 6.8143·10⁻⁵; 19941/(Λ²√(16L_A)) = 1.67025·10⁻⁴ / 2.16006·10⁻³ | ✓ |
| 4 | (2b) 164.51 (:234) | UT3 §3 | 2L_AdS + 2·0.04·L_A = 164.5115 / 125.7136 | ✓ |
| 5 | 2ζ(3) at 2.6 % (:330-331) | A1PRIME §7 | 2.5903 % (1000), 3.337 % (300) | ✓ |
| 6 | **Binding state** +0.2009, ≥ +0.198, L = 78 (:257-259) | NEARGAP2 §0; `neargap2/cache/ng_bind_J5.txt` | own row +0.20088; total budget 2.59·10⁻³ (truncation 1.91·10⁻³, ν < 50 convention 3.1·10⁻⁴, edge 3.7·10⁻⁴); lower +0.19829; L = 77.7; certified-floor route +0.19993 | ✓ (note: AdS correction at this state is −2·10⁻⁵; the state is essentially flat) |
| 7 | **J = 29** +0.972 ± 0.018; old −3.6 / +1.4 (:264-268) | HIGHSPIN §2; `highspin/cache/hs_odd_J29_final.json` | `main` = 0.97159 (dir/gk), `jk` = 0.01830; budget E_tot = 0.697, L = 1.39 (HIGHSPIN) → 1.22 with the measured window (`L_lowmeas`); interpolation variants 0.67–1.44 | ✓ value; ± is MC only (M3) |
| 8 | J = 21, 29 leverage 1.2 at the edge, ≥ 1.7 above; sign-free J = 29 0.94 (:261-262, :270-271) | SMALLNU_WINDOW §0 | J = 29: 1.23 (signed), 2.2–5.7 above, 0.94 sign-free ×3; J = 21: 1.20, 1.7–6.2 | ✓ |
| 9 | **Theorem W** J\* = 50 (100 at margin 2) (:309) | PW_IMAGES TW′; `pwimg/cache/pw4_theoremW.txt:64-65, :80` | "(3, 10, 10), θ = 1: J\* = 50; failing [26, 30, 33, 34, 38, 42, 46]"; θ = ½: J\* = 100 | ✓, for (Λ_amp, K_sh, K₄) = (3, 10, 10) on sampled states only (M1) |
| 10 | **Theorem L** 43× at J ≤ 10⁶; J\* = 857 (:312-315) | `pwimg/cache/pw5_theoremL.txt:14-16, :44-46` | J = 857: m/Λ = 1.648, β = 1.000, tolerated scalar-shaped vector 1.24·10³; J = 10⁶: 43.3; shift-by-two without light-cone normalisation K_tol = 0.399 | ✓ for the scalar-shaped vector (M1.5) |
| 11 | **Magic-block identity** [P] to r⁸, [DV] to r⁴⁰ (:303-305) | PW_IMAGES PW1; `pw1_magic.txt` | **Rerun at fresh ν:** recursion = `ee_blocks.block_coeffs_b` to 3·10⁻¹⁴; 50-digit residual ≤ 6.8·10⁻⁴⁸ (n ≤ 30, b = 0…4); sympy exact, n ≤ 4, b = 0…4: 0 nonzero of 12–20 coefficients; control with coefficients (1, 1): residual 6·10³ | ✓ |
| 12 | Monopole D = −9.98332·10⁻⁴ (supplement.tex:968) | `pwimg/cache/pw3_cost.txt:34` | −9.9833211·10⁻⁴, spread 6.4·10⁻¹⁵ over t ∈ [0, 12]; grid halving 1.4·10⁻⁵ | ✓ (sup of the cost kernel is 0.32 at t = 0; the "≤ 10⁻³ for t ≥ 0.03" qualifier is needed) |
| 13 | m_gap(J) crosses 1.3Λ at J = 345 | (Δ_g−2)(Δ_g+2J−1) | √(1+2·344/999) = 1.2996; √(1+2·345/999) = 1.3003 | supports M1.1 |
| 14 | 0.4602Λ; ε\* = 0.038034 | A1PRIME | √(0.430430/2.032073) = 0.46024; 0.430430/11.317027 = 0.038034 | ✓ |

---

## 6. Physics checks

* **Dictionary normalisations.**
  * *8πG.* C_T^Einstein = 3L²/(π³G) and C_T^free scalar = 3/(32π²) give c_T = 32L²/(πG), so 8πG = 256/c_T. ✓
  * *τ = 240|ĝ₃|.* 2(d−1)d(d+1)(d+2) = 240 at d = 3 (720 at d = 4) for the t₄ half, [DV]. The α_T half is [D] with an unfixed sign,
    stated at :122-123. ✓
  * *Energy-flux disk.* τ ≤ 4 ⇔ |ĝ₃| ≤ 1/60 (caption). ✓
  * *Free fields.* The free-field values t₄ = ±4 (scalar / Majorana) are standard. ✓
* **No flavour-current / c_T confusion.** c_T enters only through 8πG. Spin-1 and spin-3 single traces, including R- and flavour
  currents, are correctly excluded from T × T by the ⟨TTJ⟩ structure count (supplement.tex:642-643). No current two-point coefficient is
  used anywhere in the Letter. ✓
* **A1′ versus ABJM.**
  * A1′ holds for IIA-ABJM by the spin count alone: supergravity KK spins ≤ 2, and the first massive string level has J ≤ 4 at twist
    ≈ ML − 4.
  * ML = 2^{5/4}π^{1/2}λ^{1/4} follows from L²/α′ = π√(2λ). ✓
  * The scalar condition of (2a) needs the separate exchange argument (CPY/BCP), at the orders computed (m18). Massive string scalars sit
    at m ≈ Λ > 0.46Λ and are heavy. ✓
* **The logic from positivity to the bounds.**
  * P_O ⪰ 0 together with (Pos) makes the heavy side ≥ 0, so Ω[light] ≤ 0.
  * The light side of Ω^EE is K|ĝ₃|²Λ⁸ − 16L_A − (c₁ term), and its quadratic in τ gives Eq. (2) exactly.
  * For item (2), maximising over (t₄, α_T) is correct given κ > 0.
  * **Caveat M2:** H_K, the step that lets one write Eq. (1) with an absolutely convergent heavy sum, is proved under the
    non-perturbative J₀ ≤ 1 and applied at tree level, where J₀ ≈ 2.
* **The large-spin "conditional theorems".**
  * *Theorem W.* As implications both are correct and simple: Lemma I plus the 2×2 criterion d ≤ f₊₊f₊₋/(f₊₊+f₊₋). I checked the
    algebra, including (f₊₊−d)(f₊₋−d) ≥ d² ⇔ d ≤ H. The premises are stated completely in the Supplement, apart from the P0 grade for
    m < 1.3Λ (m22). They are stated incompletely in the Letter (M1.3).
  * *Theorem L.* Same as Theorem W; it is a statement about J ≤ 10⁶ and a scalar-shaped correction.
  * *What is presented as proved but is not:*
    * Lemma LR∞ "for all β ≥ 1" uses grid-evaluated constants and a [D] tail (m5);
    * "proved in the flat limit" for Ω^Q_up (M4(a));
    * Theorem W's statement for unsampled spins (M1.2).
  * *What is genuinely proved:*
    * the b-shifted magic-block identity to O(r⁸), which I reproduced;
    * the congruence lemma (trivial and correct);
    * the closed-form densities of index-shifted harmonics.
* **c₁ and F₂.**
  * The bounds are exact in c₁ and correct as functions of it.
  * c₁ = 0 needs F₂ [A], and is otherwise uncomputed.
  * The < 1 % thresholds are right: |c₁| < 59.6 for (1) and 1306 / 1262 for (2a) / (2b).
  * Being universal, c₁ is computable from one tree-level W³-vertex ⟨TTTT⟩ diagram. Its absence is the most visible hole after (Pos)
    (M5).
* **H_EFT(b).**
  * g̃₄ = g₄ + O(Λ⁻²) needs the standard CHMRS/CLS-D assumption [A], stated at :195-197. The abstract says "an R⁴ coefficient defined by
    the CFT". ✓
  * Note that H_EFT(b) is also implicitly what makes the tree-level sum rule meaningful (M2).
* **κ > 0.**
  * κ is needed only if τ ≠ 0. The flat value is 0.1099, so the sign is plausible at Δ_g = 1000 unless the O(Λ⁻²) coefficient is
    ≳ 10⁵.
  * It is graded [D-sketch + M]. The Letter could say "flat value 0.11, corrections O(Λ⁻²)" in one clause.
  * The sign matters: with |ĝ₃|Λ⁴ up to 120, κ|ĝ₃|²Λ⁸ can reach ≈ 1.6·10³ ≫ 92.8.
* **The (2a) light-scalar condition.**
  * Ω^Q_up is proved negative on flat light scalars with m < 0.46024Λ (exact polynomial −0.4304 + 2.0321m²), so a condition is needed.
  * As worded it is binary and strong: any φW² coupling, however small, violates it. The weighted formula of UT3 §3 is the right
    statement (m17).
  * (2b) removes the condition at the cost of (Pos) for Ω^EE, and the Letter says so. ✓
* **Tree level versus finite c_T (M4(b)).** The 10⁻⁸ statement is about the O(c_T⁰) part of t₄. A2's c_T ≫ Δ_g² does not make
  1/c_T corrections to t₄ (generically of order (n_B − n_F)/c_T) smaller than τ₀ at Δ_g = 1000.

---

## 7. What would most raise PRL acceptance odds from here (ranked by impact ÷ effort)

1. **Settle the tree-level Regge input (M2).** *Effort: days to state it honestly as an assumption [A]; weeks if a tree-level H_K proof
   with J₀_planar < 2 is attempted.* A CHMRS/CLS-D referee will raise this first, and an unanswered version undermines the grade of every
   bound.
2. **Rewrite the main text around the physics (M4, M7, M1 wording).** *Effort: 1–2 days, no computation.*
   * Cut Status and Large spin to about 250 words.
   * Add a (J, m/Λ) coverage figure by grade, with the 345–856 band and both theorems' regions marked.
   * Fix the "proved in the flat limit" and "every holographic CFT₃ … 10⁻⁸" framing.
   * List Theorem W's full premises, or drop J\* from the Letter until P2FOLLOW returns.

   This is the highest return per hour, and it decides whether an editor sends the paper out.
3. **Compute the spinning O(1/m²) coefficients (P2 step 3: the three-channel EE contraction) and K₄ with the ν⁶ subtraction.** *Effort:
   high, 1–3 weeks plus review.* This turns Theorems W and L into statements with derived premises, the only route by which "large spin"
   stops being a hypothesis.
4. **Close the near-gap band.** *Effort: moderate compute (≈ 1–2 CPU-h per spin for direct bulk rows, by SMALLNU's and HIGHSPIN's
   figures), plus a structural argument for J ≳ 32 where the window is live again.*
   * Error-budget J = 7, 9–11, 13–15, 17–20, 22–28, 30–48.
   * Give the 345 ≤ J < 857 band (m ≈ 1.3–1.65Λ) at least sampled [NE±] budgets.
5. **Compute c₁.** *Effort: moderate–high, one spinning AdS₄ tree diagram with a W³ vertex.* It removes "up to one uncomputed universal
   constant" from the abstract, and c₁ is universal, so one theory suffices.
6. **Independent review of the v5 layer.** *Effort: low–moderate.* Interval constants for Lemma LR∞ and recomputed Theorem W/L costs at
   2–3 spins (CHANGES_v5 §7 item 3 already flags these as cheap). The authors' grade discipline is only as good as its audit.
7. **Positivity over a range of gaps.** *Effort: high.* For example, show that the de-shifted corrections scale monotonically in Λ, so
   that (Pos) at 1000 implies (Pos) for all Δ_g ≥ 1000. This converts a data point into a constraint on a family, and is the single change
   that most raises broad interest.
8. **Re-optimise Ω^Q with light-scalar positivity imposed.** *Effort: high.* It removes the 0.46Λ condition and the ε-repair, and makes
   item (2) clean.

The realistic path is a long paper now (PRD or JHEP), containing everything above that is done, with a short PRL only if items 1, 3 and 4
(or 7) succeed.
