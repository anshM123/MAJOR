# Fifth referee report

**Manuscript:** "How Close to Einstein Must Holographic Gravity Be?", v10 (`ads4/prl/letter.tex`, `supplement.tex`, `macros.tex`,
`refs_verified.bib`, `fig/`).
* Builds: `prl/build/letter.pdf` (6 pp.) and `supplement.pdf` (30 pp.), dated 24 Nov 11:23. Both are newer than every source file.
* Switches: `\HLJclosedfalse` (headline Δ_g = 1000), `\HSdraftfalse`.
* Pending marks: `% SEAM` and `% BATCH4` (seam spins and external direct rows at Δ_g = 2000). I treat both as pending. I do not count
  their absence against the paper, except where the Letter's wording already assumes their outcome.

*Written 24 Nov 2026 as an internal pre-submission review. I edited no existing file and committed nothing.*
* References are `file:line` in the current text, relative to `ads4/prl/`. "L" means letter.tex and "S" means supplement.tex.
* Scratch scripts are in `/tmp/claude-0/-home-user-MAJOR/60eea3a7-07d6-5453-af55-c4346d86f9b5/scratchpad/r5/`. One of them,
  `rerun_bind.py`, imports project code in-process. It runs behind a write guard that refuses every write outside the scratch
  directory, and `git status` was clean afterwards.

*Method.*
* I read the Letter, the Supplement and the figures cold, and wrote §§1–3 and §6 before opening the earlier reports.
* I then checked report 4 item by item (§4).
* I checked the new all-gap material against `ads4/UNIFIED_THEOREM_v3.md`, against the dated boxes of `ads4/UNIFIED_THEOREM.md`
  (through GAP_MORSE), and against the deliverables `GAP_SCALING`, `GAP_LARGESPIN`, `GAP_UNIFORM`, `GAP_C3`, `GAP_NEARGAP` and
  `GAP_MORSE`, including their caches.

---

## 1. Summary and recommendation

**What the Letter claims (v10).**
* **Setting.** A unitary CFT₃ with a unique stress tensor and no light single trace in T×T, at leading order in 1/c_T. Assume that two
  fixed light-ray functionals on ⟨TTTT⟩ are positive on every heavy state.
* **Cubic bound.** At Δ_g = 1000, |ĝ₃|Λ⁴ ≤ 120, i.e. τ = √(t₄²+α_T²) ≤ 2.89·10⁻⁸.
* **Quartic bound.** g̃₄Λ⁶/8πG ≤ 92.8. This needs κ > 0 or τ = 0, and no scalar lighter than 0.46Λ.
* **c₁.** Both bounds are exact in one universal constant c₁; the numbers quoted are for c₁ = 0.
* **New in v10:**
  * (i) A remark: a theory with a larger gap also has a gap of 1000, so the 1000 numbers bound it.
  * (ii) The per-gap bound τ₀ ∝ √(ln Λ)/Λ⁴. In gap units it grows slowly: 120 / 128 / 139 at Δ_g = 1000 / 2000 / 5000.
  * (iii) An argument that positivity of Ω^EE "carries over to every Δ_g ≥ 1000". It rests on certified uniform-in-gap layers,
    near-gap rechecks at 2000 (two spins at 5000), a "no-growth" premise, a certified remainder bound above Δ_g ≈ 10⁵, and large-spin
    margins that grow with the gap.

**Assessment.**
* **The implication and its constants are correct.**
  * I recomputed every gap-dependent constant from the published exact weights by a different quadrature (§5).
  * I reran the binding-state error budget at Δ_g = 2000 from the project's rows, and it reproduces to the printed digit.
  * The monotonicity remark is right (§6.1).
* **v10 is a real step on the breadth objection.**
  * The certified layers are now gap-uniform: the flat floor, the monotone Lemma C, LR∞ for every Λ ≥ 998.5, and the ROUTEA window.
  * The structural Lemma T is genuinely useful.
  * The 2000 rechecks show the expected Λ⁻² shrinkage on every budget term for J ≤ 16.
  * All of this is honest work, and S8 grades it honestly (S:1403-1413, S:1702-1705).
* **The category of the result has not changed.** The physical content is still positivity of a fixed functional on an infinite AdS
  spectrum, and that positivity is not established at *any* gap:
  * Below J ≈ 25 near the gap, ten spins are still model-level at 1000.
  * The graviton's edge integral (EDGE-R) is still uncomputed.
  * The binding states of the quartic functional are still [NE + M].
  * (AG_U) is still an assumption.
  * By the Supplement's own net statement (S:1411-1413), the all-gap extension inherits this status. It adds a no-growth premise
    below Δ₁\* ≈ 10⁵ and model-level premises at large spin.
* **The Letter states the extension more strongly than S8 does** (M2), and leaves one band of states unassigned (M3).

**PRL criteria.**
* **Importance / broad interest.** The topic is current: explicit-constant graviton bounds from ⟨TTTT⟩ are what CLS-D list as missing.
  The ideas are genuine (EE structure, de-shifting, the k = 4 arc). But the result is still conditional, and its weaknesses are
  unchanged:
  * the constant is ≈ 10× weaker than the transplanted flat optimum that CHLPS already remark on;
  * c₁ is uncomputed;
  * the quartic bound is set by the infrared logarithm and holds at one gap;
  * the one theory checked (ABJM) sits at 2.6 % of the bound.

  The new large-gap direction adds the least where it is most rigorous. As Δ_g → ∞ the statement approaches the flat bound uplifted
  by the CHMRS/CHLPS argument, which those authors already assert. The value added is finite-gap control, and that is exactly where
  positivity is weakest.
* **Correctness.** I found no error in the implication, in the constants, or in the monotonicity remark. The physical premise is not
  established.
* **Claims at the level stated.** The title, the abstract and the "Large spin" paragraph are now accurate. The new "Larger gaps"
  paragraph, Fig. 2's caption and the Introduction's new sentence are not (M2).
* **Clarity.** The writing is improved, but "Larger gaps" is dense specialist prose with nine numbers (m2).

**Would an editor send it out?** Probably yes, though no more readily than v9. I estimate 45–60 %.
* For sending it out:
  * the honest question title;
  * an abstract that leads with the condition;
  * a current topic.
* Against:
  * the abstract's "If … are positive" announces a conditional result in its first word;
  * a 6-page Letter with a 30-page Supplement that reads as a grade ledger;
  * independent authors;
  * an AI-use acknowledgment.

**What would a strong referee decide?** From the CHMRS/CLS-D community: reject for PRL, recommending PRD or JHEP. They would grant
that the "one gap" objection, in its literal form, is now answered. They would then say:
* the binding objection is positivity in AdS at the headline gap, not the number of gaps;
* the "every gap" statement rests on the same unproven core plus premises;
* c₁ and the graviton edge integral should be computed before an explicit-constant claim is advertised;
* the supplementary material should be rewritten as derivations.

A more sympathetic referee might ask for major revision. However, the paper cannot meet "claims established" without new results on
positivity at Δ_g = 1000.

**Has the "one gap" objection been answered adequately?** Partly, and adequately for what can be done without solving positivity.
* *The literal objection is answered.* Report 4 M3 said "a constraint at one value of the gap is a data point, not a constraint on a
  family". The remark (L:250; S:1416-1419) is correct: given (Pos) at 1000, τ ≤ 2.89·10⁻⁸ holds for every theory with Δ_g ≥ 1000.
  This should have been in the paper from the start. It answers the "family" point in absolute units.
* *The scaling-law version is answered only conditionally.* The coefficient of the Λ⁻⁴ law at each gap (|ĝ₃|Λ⁴ ≤ 128, 139, …) needs
  (Pos) at each gap. That rests on:
  * for J ≤ 7, [P given B_F + window conventions], with B_F [NE] (sampled masses only);
  * for J = 8–32 near the gap, [D-sketch given H_mono-A [M]];
  * for the near-gap seam, chart B under H_mono-B [M];
  * at large spin, premise constants carried in AdS units [M];
  * an unassigned band (M3).

  This is "at best the 1000 grade plus premises" (S:1411), as the Supplement says. It is not a proof that (Pos) at Δ_g implies (Pos)
  at Δ′_g > Δ_g. Report 4 asked for either that proof or "positivity over a range of gaps"; v10 gives the second at two gaps (partly
  three), with a transfer premise.
* *Net.* The one-gap objection is no longer fatal on its own. What remains fatal is report 4's M2, which v10 does not touch at 1000.

**Recommendation: reject for PRL in the present form; transfer to PRD (or a JHEP long paper) strongly encouraged.**
* This is unchanged from reports 2–4, but the reasons have narrowed. Breadth is now mostly a wording and coverage problem (M2, M3).
  Establishing positivity at 1000 is the substantive one (M1).
* M2, M3, M5 and the minor items must be fixed whatever the venue.

---

## 2. Major issues

### M1. Positivity in AdS is still not established at any gap (carried over from report 4 M2)

**Where:**
* L:79-85, :260-272, :278-299;
* Fig. 2 (L:322-334);
* S Tables S5, S7, S11 (S:890-955, :1361-1381);
* S Sec. S8 net statement (S:1411-1413).

**State at Δ_g = 1000** (unchanged in v10 except the two upgrades noted below):
* *Error-controlled ([NE±], a ×3 budget on model terms; an estimate, not a bound):*
  * J = 5;
  * J = 0, 2, 4, 6, 7, 8, 12, 16, 21, 27 from their twist edges;
  * J = 29 only with a signed window;
  * J = 23 at m ≥ 1.03Λ and J = 31 at m ≥ 1.05Λ;
  * the strip m ≥ 1.3Λ, b < 5.
* *Model level ([NE + M]):*
  * near-gap J = 9–11, 13–15, 17–20, 22, 24, 25;
  * the twist edges of 23 and 31;
  * the band 1.1–1.3Λ (H_cont);
  * 5 ≤ b < b_A;
  * every binding state of Ω^Q_up.
* *Conditional theorems:*
  * Theorem W on 47 sampled states, needing AG_U [A], EDGE-R [M] (graviton uncomputed; scalar analogue 0.405 against 0.406) and model
    magnitudes;
  * Theorem L′ (AG_U, EDGE [M], P_L′(i) [M] above J = 344).
* *Hypothesis only:*
  * between samples;
  * intermediate impact parameter at m ≥ 1.3Λ;
  * J > 10⁶.
* **v10's upgrades at 1000** come from unreviewed GAP_LARGESPIN (S:1137-1140, :1460-1462):
  * the hybrid-floor hole 119 < b < 245.2 at 1.3 ≤ m/Λ < 1.5 is closed;
  * Theorem W's P0 is rigorous at b ≥ 114.06.

  Both are welcome. Neither touches a beyond-leading-Regge premise.

**Problem.** The physical content of the paper is (Pos). It is now argued to transfer to larger gaps "at best at the grade it has at
1000" (S:1411). So the extension multiplies the scope of an unproven statement; it does not prove it anywhere. This is the structural
objection of reports 2 (M1), 3 (M1, M5) and 4 (M2).

**Fix (computation, not writing).** In order of value:
* error-control near-gap J ≤ 25 over 1.0–1.3Λ at 1000;
* compute the graviton EDGE-R at the Theorem W states;
* put the Ω^Q_up binding states under error control.

**Fatal for PRL?** Yes, with M4.

### M2. The Letter states the all-gap extension more strongly than Sec. S8 and the record allow

Every item below is a wording problem. Taken together, though, they are the paper's new headline, and each claims more than S8 does.

1. **L:86-87, Introduction.** "We do this at Δ_g = 1000 and carry it, assuming that near-gap corrections do not grow with the gap, to
   every larger gap."
   * "It" has no clear antecedent: the bounds? positivity? the checks?
   * The transfer needs more than one premise. The large-spin theorems need their premise constants carried in AdS units [M] (S:1643-1645,
     Table S16 caption). The Discussion's disclaimer (L:362-363) omits this too.
   * This repeats report 4 M4.1's "one technical assumption" pattern in new text.
2. **L:253-254.** "We argue below that positivity of Ω^EE carries over to every Δ_g ≥ 1000."
   * By the Letter's own legend (S:77), "argue" means [D-sketch].
   * The actual grades are weaker: [P given B_F] with B_F [NE] for J ≤ 7 at sampled masses; [D-sketch given H_mono-A [M]]; [M] in chart B.
3. **L:304, "The evidence for Ω^EE carries over to every Δ_g ≥ 1000."** This is stated as a fact, but it is conditional on the no-growth
   premise and on the carried premise constants. Write "carries over, given …".
4. **L:306-307, "We repeated the error-controlled near-gap checks at Δ_g = 2000."** This is only partly true:
   * At 2000 the twist edges of J = 21, 27 and 29 are [NE + M], not error-controlled (S Table S13, S:1512-1520; λ-level leverage
     0.4–1.1, sign-free 0.0–0.2).
   * J = 23 and 31 were not recomputed.
   * The seam spins and the parity-odd OO states have no rows (S:1541-1545).
5. **L:309-310, "Spins 21, 27 and 29 remain thin at 2000".** "Thin" is the 1000 grade ([NE±], L ≈ 1.2–1.6).
   * At 2000 they do not close at all at the twist edge.
   * Write "are not yet error-controlled at their lowest mass at 2000 (at 1000 this needed direct evaluations)".
   * Part of the shortfall is the reference-zero artefact (S:1530-1533). That explains it; it does not change the grade.
6. **L:315-316, "We checked this at 1000, 2000 and 5000".** At 2000 and 5000 the check uses the K₂ Laplace rows.
   * Their chart-A gap independence is structural (Prop. S: the phase Φ_ξ is exactly Δ-independent; S:1467-1468).
   * The "ratio 0.97–1.06" (Table S14) therefore tests the asymptotic function more than the true correction.
   * The non-asymptotic tests are the exact scalar (to 10⁴) and the 1000 direct rows.
   * The first graviton test at 2000 is batch 4, which is pending.
   * Say so, or wait for BATCH4.
7. **L:317-319, "Above Δ_g ≈ 10⁵–4×10⁵ a certified bound on the saddle-point remainder replaces it for the lowest spins".** It does not
   replace the premise. Premise (iii) of the statement (S:1611) is sup_m|Γ_{K₂,m}| bounded, i.e. B_F for the explicit K₂ function,
   graded [NE]. Further limits:
   * the certificate covers five ξ-nodes, and ξ-continuity is [D-sketch];
   * the graviton amplitude is certified only at J = 0, 2 (S:1613-1615).

   Write: "a certified remainder bound reduces it, for spins 0 and 2 (and for a scalar analogue up to spin 8), to the same premise for
   an explicit function".
8. **L:319-321, "At large spin the theorems' margins grow".** This holds only with Δ_g-independent premise constants in AdS units [M]
   (S:1643-1645; GAP_LARGESPIN U6). It also holds only at the chart-B images of the 47 sampled states. Add "if the premises' constants
   do not change with the gap".
9. **Fig. 2 caption (L:332), "In these units every larger gap looks the same or better."** Not true in the (J, m/Λ) plane drawn.
   * J\* of Theorem L′ moves from 223 to 447 to 1119.
   * The Theorem W samples move to J ≥ 52 and J ≥ 132.
   * The certified sign-free near-gap reach grows only as J_LC ∝ Λ^{1/3} (13 / 17 / 23 / 35 / 61; `code/gapscale/cache/lcscan.txt`),
     while the near-gap band to be covered grows ∝ Λ (to J\* ≈ 0.224Λ).
   * Chart A at fixed J makes the EDGE-R transfer *harder*, ∝ (Δ_g/J)² (S:1599-1601).
   * Delete the sentence, or replace it with a correct one, e.g. "At larger gaps the near-gap corrections shrink like Δ_g⁻²; the
     large-spin regions move to proportionally larger J".
10. **Fig. 1 caption (L:175-176), "The solid cubic curve above Δ_g = 1000 holds given the premise …".** "Holds" reads as "is a bound".
    The curve is a bound given (Pos) at that gap, for which the premise is one input among several. Write "is supported, given the
    premise …, at the grade it has at 1000".

**Fix.** Rewrite L:86-87, :253-255, :303-321, the Fig. 1 and Fig. 2 captions and L:362-363 along these lines. Net cost ≈ +20 words,
to be paid by cuts (m1).

**Fatal?** No. But "every larger gap" is now the most prominent new claim, and a referee will test it first.

### M3. A band of near-gap states at larger gaps is assigned to no transfer route

**Where:** S Table S17 (S:1671-1698); S:1586-1602; GAP_LARGESPIN §0, table row "the twist-edge seam at small J, non-sampled J, J < 26:
not addressed here".

**Problem.** At Δ_g = 2000, consider the near-gap states (above the seam) with 33 ≤ J ≲ 51. At 5000 the band is 33 ≤ J ≲ 131.
* Table S17 covers "other near-gap states, J ≤ 32" by H_mono-A inheritance.
* It covers Theorem W at "their chart-B images", which start at J = 52 / 132 (the images of J = 26).
* Neither route covers the band in between:
  * *Chart A* maps these states to 1000 states with J = 33–51, which are Theorem W territory. There the Supplement itself says chart A
    "is not usable": the +− EDGE-R transfer factor grows like (Δ_g/J)² (S:1599-1601).
  * *Chart B* maps them to 1000 near-gap states with J ≈ 16–25. Those are mostly [NE + M] at 1000, and "in chart B pointwise
    corrections do not shrink, so sign-free per-node budgets do not improve" (S:1600-1601).
* As far as I can tell, this band carries no grade at any Δ_g > 1000. It grows with the gap.

**Fix.** Either add a Table S17 row with its honest grade (I expect [NE + M] or "hypothesis only"), and mark it in any figure that claims
all gaps; or compute K₂ rows there at 2000 (cheap, by GAP_NEARGAP's protocol; with Δ_ref ≥ 5·10⁴).

**Fatal?** No, but it contradicts the Letter's "carries over to every Δ_g ≥ 1000" (L:304) as written.

### M4. Importance and breadth for PRL (report 4 M3, now narrowed)

**Where:** abstract; L:95-103, :241-255, :350-368; Fig. 1; S Sec. S5.G.

**What remains after v10:**
* **c₁ is uncomputed.** It is universal ([D]); c₁ = 0 rests on F₂ [A]. One tree-level spinning Witten diagram with a W³ vertex fixes
  it. The abstract's last sentence advertises the gap.
* **Loose constants.**
  * |ĝ₃|Λ⁴ ≤ 120 is 10.3× the transplanted flat optimum (S:816-818). The Letter attributes this to the missing forward and
    helicity-flip sum rules (L:353-356).
  * In gap units the constant now grows like 55√(ln(Λ/8.74)) (S:1423). The "number rather than a scaling law" of L:102 is a slowly
    growing function, as the Letter now says at L:251-253.
* **The quartic bound is a Δ_g = 1000 statement only** (L:254-255; S:1428-1429). It is set by the infrared logarithm (17.24 ln(2Λ/3) −
  19.3).
* **The ABJM check sits at 2.6 %** of the bound and can only fail grossly.
* **Where the extension adds least.** As Δ_g → ∞ the large-gap statement approaches what the flat-space bound plus the CHMRS/CHLPS uplift
  ("up to corrections suppressed by 1/Δ_gap", quoted at GAP_UNIFORM §1.1) already asserts, with a 10× better constant. The new
  explicit content is at finite gap, where positivity rests on B_F (1000 ≤ Δ_g < Δ₁\* ≈ 10⁵, "about two decades", S:1618).

**Fix.** Compute c₁. Extend Ω^Q_up to 2000 (K₂ rows at its binding states are cheap) if the quartic bound is to stay in the
headline. Otherwise this remains a long-paper result.

**Fatal for PRL?** Yes, together with M1.

### M5. The functional at gaps other than 300 and 1000 is neither published nor shown to exist with the required properties

**Where:**
* L:185-186: "All weights are published as explicit functions";
* S Table S1 (S:264-274): only 300 and 1000;
* S Table S12, first row (S:1439-1442): "weights [DV] at 300–5000";
* S:1420-1424 and Fig. 1, which quote bounds at 2000, 5000, 10⁴, 10⁵.

**Problem.** The cubic number at Δ_g needs only ψ₂ (L_A, K₀). But the theorem behind it needs more at Δ_g itself:
* the interchange H_K needs the de-shifted k = 2 weight to exist and be strip-analytic, or the continuum weight with an ε_W bound (row
  DS, S:296);
* the "certified at every gap" leading-Regge statements use "the leading-Regge action is exactly ∫ψ_kP at every Λ" (S:1439-1440),
  which holds only if the de-shift identity holds at that Λ.

What the paper has:
* the 2000 and 5000 weights exist only in project caches (`code/gapscale/cache/weights_{2000,5000}.npz`);
* no continuum weight or ε_W is stated at any Λ ∉ {300, 1000};
* for Λ > 5000 the weight has not been constructed; there is only the expansion Re X = ψ₂ + S/Λ² + O(Λ⁻⁴) [DV].

**Fix.**
* Publish the 2000 and 5000 continuum weights with their ε_W.
* State the construction of Ω^EE at general Λ and the grade of its existence.
* Restrict the "every gap" wording to "[DV] at 300–5000; asymptotic expansion beyond".

**Fatal?** No; a half-day of writing plus an export.

### M6. The unreviewed record has grown, and the AI-use disclosure does not match the Supplement

**Where:**
* S:114-123 (provenance: "none of these is externally refereed"; every S8 analysis "unreviewed");
* S:1401 ("None of the analyses quoted here is reviewed or independently reproduced");
* acknowledgments L:370-373.

**Problems.**
* **Six new unreviewed analyses.** v10 adds GAP_SCALING, GAP_LARGESPIN, GAP_UNIFORM, GAP_C3, GAP_NEARGAP and GAP_MORSE. They are
  unreviewed, and two of them already change 1000 statements: the floor hole and P0.
* **The external review report 4 asked for (M7) has not happened.** The Supplement grew from 25 to 30 pages and still ends in a version
  history ("Corrections to earlier drafts", S:1841-1885).
* **The acknowledgment contradicts the Supplement.** It says the authors "have thoroughly reviewed all output". The Supplement says in
  many places that analyses are unreviewed. An editor or referee who reads both will notice.
* **The disclosure is narrower than what the record shows.** It covers "formulation and drafting of this manuscript". The record
  (deliverables headed "gap-near agent", "gap-Morse agent", …) shows AI agents doing computations, code, certificates and derivations.

**Fix.**
* Make the disclosure match reality. For example: "AI tools (Claude) were used extensively in this research, including code,
  numerical computations, computer-assisted certificates, derivations and drafting; the authors take responsibility for all results."
* Repeat it in the Supplement.
* Replace "unreviewed" by a neutral statement of what has been independently reproduced, or drop "thoroughly reviewed all output".
* Obtain one external read of S7–S8 and of the c₁ closure.

**Fatal?** Not by itself. It is an editorial risk, and an integrity point that must be right before submission.

---

## 3. Minor issues

**Length (APS counting).**
* **m1. Compliant on the standard count, with no headroom left.**
  * *The authors' count.* I reproduce it: `wc_prl7.py letter.tex 3 0.9375 1.0385` gives 3056 prose + 236 formulas + 48 + 344 =
    **3684** (conservative 3746).
  * *My count* (`r5/wc5.py`):
    * acknowledgments excluded, as APS does; each inline formula one word; citations not counted;
    * that gives **3656**;
    * with displays counted at 5–6 rows (Eqs. (2) and (3) are tall), 3688–3704.
  * *A hostile reading* counts each citation bracket as a word and tokenises long formulas (e.g. counting relations). It gives
    **3830–3880, i.e. 2–3 % over.**
  * *Abstract:* 556–580 rendered characters, against a limit of 600.
  * *Pages.* The text ends at the top of page 5 and the references fill pages 5–6. PRL counts words, not pages, so 6 pages is not a
    violation, but it looks long to an editor.
  * *Verdict.* Within the limit on the standard convention, with ≈ 50–90 words of headroom. Any referee-driven addition, including the
    SEAM/BATCH4 updates and M2, must be paid for by cuts.
  * *Candidate cuts:*
    * "Large spin" (316 words; mean sentence 32 words) to about 200;
    * "Larger gaps" (209 words) to about 130;
    * Discussion (240 words; mean sentence 34) to about 180.

**Prose** (meant to read as plain, human prose with few numbers).
* **m2.** Most of the Letter now reads well. Problems:
  * *"Larger gaps" (L:303-321)* is the densest paragraph in the paper. It carries nine numbers ("Δ_g⁻²", "four- to ninefold",
    "10⁵–4×10⁵", "140 times", …) and undefined jargon: "flat-space floor", "leading-Regge floor", "near-gap", "saddle-point
    remainder", "graviton edge integral". Keep one number (the Λ⁻² shrinkage) and push the rest to S8.
  * *Long sentences:* L:88-92 (≈ 60 words), L:250-255 (≈ 65), L:279-290, L:361-364 ("We do not claim …" lists nine items).
  * *L:268, "Three of these spins are covered only above their lowest mass; at J = 29 …".* The three are 23, 29 and 31
    (`fig/make_fig2.py:69-70`), so a reader counts four. Write "Spins 23 and 31 are covered only above their lowest mass, and at
    J = 29 …".
  * *L:246-247.* "the sign of κ, which also rests on a model, are only sketched": the relative clause dangles.
  * *L:155, "Three ingredients overcome it".* Since positivity is not established, write "make it tractable".
  * *L:309-310, "until direct evaluations, needed at 1000 too, are included".* Clumsy; see M2.5.

**Title and abstract.**
* **m3. Title.** "How Close to Einstein Must Holographic Gravity Be?" is acceptable and no longer overclaims. It does not say AdS₄/CFT₃,
  and the answer is conditional. Consider "… in AdS₄?" or a subtitle.
* **m4. Abstract.**
  * *ABJM.* "a bound ABJM theory obeys" is still in the abstract. Reports 1 (M13) and 4 (M1.6) asked to move it to the body. It is a
    2.6 % consistency check.
  * *κ.* "one quartic coupling is bounded above" omits κ > 0 (or τ = 0).
  * *What is right.* Leading with the condition, giving |ĝ₃|Λ⁴ ≤ 120, restricting the flat-space proof to the cubic bound, and stating
    c₁ = 0 are all good.
  * *Δ_g ≥ 1000.* Given the monotonicity remark, the abstract could say "a gap Δ_g ≥ 1000" for the absolute τ ≤ 2.89·10⁻⁸. The authors
    are right not to write "≥ 1000 … 120 in units of the gap" (CHANGES_v10 §2). Phrase it so that only τ carries "≥".

**Wording stronger than the Supplement's grades** (beyond M2).
* **m5.** L:223, "x ∝ |c₁|/Λ²". This is still a proportionality (report 4 m5). Give x = 1.67·10⁻⁴|c₁| at Δ_g = 1000 (S:340).
* **m6.** L:263-264, "positivity of the leading Regge part of Ω^EE on every state whose impact parameter exceeds about a third of the
  AdS radius" is certified for the exact minimal-norm weight. For the defined continuum weight it holds up to ε_W ≈ 10⁻⁶ [DV] (report 4
  m9, still unstated in the Letter). Trivial, but it is the only place where "certified" applies to an object other than the one
  defined.
* **m7.** S Table S17, row "near gap J = 0, 2, 4–7: every larger gap [P given B_F + windows]" (S:1681). The Proposition covers only the
  sampled masses (twist edge, 1.02, 1.05, 1.1, 1.3Λ; S:1557). Between them, H_cont [NE]. Say so in the row.
* **m8.** S:1615-1616: "certified radius r_c ≈ 0.014–0.052". GAP_MORSE §0 gives 0.023–0.052 (scalar) and 0.014–0.032 (graviton
  prefactors). Give both ranges.

**Figures.**
* **m9. Fig. 1** (report 4 m7, not addressed).
  * The transplanted flat optimum still carries the paper's own marker code: a filled "claimed" dot at 1000 and an open "one
    hypothesis open" diamond at 300 (`fig/make_fig.py`, `series()` draws both for every series; visible in the preview). Use plain
    markers for the comparison curve.
  * "(2b): light scalars allowed" should say "needs both functionals".
  * The open circles at 2000 and 5000 are fine.
* **m10. Fig. 2** is unchanged since v9 apart from its caption.
  * The caption's new sentence is wrong (M2.9).
  * If the all-gap claim stays, an inset or a second x-axis in J/Λ would show what transfers in chart B and what does not (seam;
    M3 band).

**References.**
* **m11.** Report 4 m13 is not addressed. I would add:
  * Kravchuk–Simmons-Duffin, arXiv:1805.00098, for the light transform (L:131);
  * Penedones–Silva–Zhiboedov, arXiv:1912.11100;
  * Hartman–Jain–Kundu, arXiv:1509.00014, with the Regge-causality references (L:56).
* **m12.** Literature on uniformity in the gap. S8 should cite what GAP_UNIFORM §1 found:
  * CHMRS §5's "sufficiently large M" uplift;
  * CHLPS's one-sentence uplift with an unquantified O(1/Δ_gap⁴);
  * Knop–Mazáč's pointwise-kernel mechanism (arXiv:2203.11170);
  * van Rees–Zhao (arXiv:2210.15683).

  The Letter's "Larger gaps" paragraph could then say in one clause that no explicit uniform-in-gap bound exists in the literature.
  That is a fair point in the authors' favour.

**Source hygiene.**
* **m13.** Internal comments remain in the source: version notes, internal file names, `% SEAM`, `% BATCH4`, `% Front matter supplied
  by the authors (v5)`. Strip them before arXiv (report 4 m12).

---

## 4. Report 4: status of each item in v10

| report-4 item | status | evidence (v10) |
|---|---|---|
| **M1.1** title promises a derivation of Einstein gravity | **addressed** | "How Close to Einstein Must Holographic Gravity Be?" (L:22); no derivation claimed (m3: scope not in the title) |
| **M1.2** "certified in flat space" said of both functionals | **addressed** | "Positivity is proved in the flat-space limit for the cubic bound" (L:41) |
| **M1.3** "hold up to one universal constant" ambiguous | **addressed** | "The numbers are for c₁ = 0, where c₁ is one universal constant we have not computed" (L:42) |
| **M1.4** headline number is scaling-dominated | **addressed** | Abstract leads with \|ĝ₃\|Λ⁴ ≤ 120 (L:39) |
| **M1.5** "we show" before the condition; quartic conditions | **partly** | The condition comes first (L:36). The quartic sentence still omits κ > 0 (or τ = 0) (m4). The scalar condition is implied by "a gap … in the single traces of T×T" |
| **M1.6** ABJM in the abstract | **not addressed** | "a bound ABJM theory obeys" (L:40-41) |
| **M2** positivity not established in AdS | **not addressed** (no new 1000 near-gap or EDGE-R data) | v10 adds two 1000 upgrades at leading-Regge level only (floor hole; P0 at b ≥ 114.06; S:1137-1140). See M1 |
| **M3** one gap value | **partly** | Monotonicity remark (L:250; S:1416-1419), correct. Per-gap constants (L:251-253). Transfer argument (S8) resting on B_F [NE] and [M] premises. See §1 and M2, M3 |
| **M3** loose constants | **not addressed** (unchanged by nature) | 10.3× transplanted optimum (S:816-818); quartic IR-dominated |
| **M3** c₁ uncomputed | **not addressed** | L:194-196; S row L3 |
| **M4.1** Intro "one technical assumption" | **addressed** for large spin at 1000 (L:82-85) | **regressed** in the new all-gap sentence (L:86-87; M2.1) |
| **M4.2** J ≥ 26 condition in the same sentence | **addressed** | "positivity on the 47 sampled near-gap states of spin 26 ≤ J ≤ 344, given the edge bound below" (L:287-288) |
| **M4.3** Theorem L′ "every state" | **addressed** | "the leading margin is certified on every state and the remainder was checked on sampled states up to mass 5Λ" (L:289-290). The reason ε_R does not grow with m is still not given (report 4 m10) |
| **M4.4** Fig. 2 draws Theorem W as solid | **addressed** | Dots at the 47 states plus a hatched band labelled "between theorem-1 samples" (Fig. 2; caption L:329-331) |
| **M4.5** "prove order by order" | **addressed** | "we prove through eighth order …, check numerically to fortieth order and derive at all orders" (L:281-282) |
| **M5** light single traces in the headline | **addressed** (option (a)) | Result assumes no light single trace (L:213); A1′ is an extension "that rests on the sketched step" (L:235-237); "given the light-side formula we sketch … also bounds their couplings" (L:206-207); representation dependence explicit (S:761-763, :791-798) |
| **M6** "establish … the low-lying heavy states" | **addressed** | "we check it numerically, with error control, on the most constraining state and a sample of other low-lying states" (L:80-82) |
| **M6** J = 29 counted under error control | **addressed** | "at J = 29 the lowest mass is within error control only with the measured sign of one low-frequency contribution" (L:268-269); m2 notes an ambiguity |
| **M6** "exactly its flat value" | **addressed** | "to about one part in 10⁶" (L:160) |
| **M6** "We show" (H_K) | **addressed** | "We derive that the OPE sum may be exchanged" (L:186) |
| **M6** "Independent code … most delicate steps" | **addressed** | "Separate code reproduces every constant, and a separate check rebuilt the EE-kernel reduction and the counting lemma" (L:247-248) |
| **M7** unreviewed record; Supplement as ledger | **not addressed; grown** | Six more unreviewed analyses; 25 → 30 pp.; version history kept at S:1841-1885; no external read (M6) |
| m1 length | **addressed but marginal** | 3684 / 3746 against 3750 (was 3327). Headroom ≈ 50–90 words (§3 m1) |
| m2 prose | **partly** | Grade tags gone; "Larger gaps" is new dense prose; long sentences remain |
| m3 "its" antecedent | **addressed** | Abstract rewritten |
| m4 "small AdS mass shift" | **addressed** | "an AdS mass shift we have not computed" (L:204-205) |
| m5 x ∝ \|c₁\|/Λ² | **not addressed** | L:223 |
| m6 binding state is essentially flat | **not addressed** (optional) | LR curvature 7.8·10⁻⁷ N_s at 1000, 1.9·10⁻⁷ at 2000 (my rerun, §5 #3) |
| m7 Fig. 1 markers; (2b) label | **not addressed** | `fig/make_fig.py` `series()`; preview |
| m8 Fig. 2 J = 27, 29, 31 bars; J = 29 distinct | **partly** | J = 29's lightest part drawn pale (caption L:328-329); overlap unchanged |
| m9 functional beyond Λ (truncated grid not admissible; bump) | **addressed** | S:155-169, :198-203 ("not separate functionals … only the continuum weight is the definition"; bump at Λ+30–33 stated) |
| m10 Theorem L′ mass column | **partly** | The Letter now names the sampled masses; no argument for m > 5Λ |
| m11 order of limits / uniformity in spin | **addressed** | "uniformly in spin" in (ii) (L:113-114); CHMRS cited at S row A2 |
| m12 internal comments in source | **not addressed** | More of them now (`% SEAM`, `% BATCH4`, v10 headers) |
| m13 reference additions | **not addressed** | None of the three added |
| m14 AI disclosure | **partly** | Added in the acknowledgments (L:370-373), but too narrow and in tension with "unreviewed" (M6) |
| §7.1 title / abstract / Intro | **done** (except ABJM and κ) | — |
| §7.2 restrict headline to no light single traces | **done** | — |
| §7.3 graviton EDGE-R | **not done** | Still uncomputed at any gap (S:1602, :1649) |
| §7.4 error-control the rest of near-gap at 1000; Ω^Q binding states | **not done** | — |
| §7.5 compute c₁ | **not done** | — |
| §7.6 external review; Supplement as derivations | **not done** | — |
| §7.7 positivity over a range of gaps | **partly done** | 2000 (5000 for J = 5, 12) [NE±] for J ≤ 16; transfer premise B_F; certified layers uniform; no proof of the implication (Pos at Δ ⇒ Pos at Δ′) |
| §7.8 AG_U | **not done** | — |
| §7.9 re-optimise Ω^Q | **not done** | — |

**Reports 1–3.** Report 4's §4 stands for them, with these changes in v10:
* **R1-M12 / R4-M1.1 (title):** fixed.
* **R1-M13 / R4-M1.6 (ABJM in the abstract):** still open.
* **R2-M4, R3-M4(a), R4-M1.2 (grade slips in the abstract):** fixed.
* **R2-M5, R3-M5 (importance):** narrowed (M4), not resolved.
* **R3-M6, R4-M7 (unreviewed layer):** larger.

---

## 5. Spot checks of load-bearing numbers

I ran nine checks; checks 3–7 and 9 concern the new all-gap material.

| # | number (Letter / Supplement) | traced to | my check | verdict |
|---|---|---|---|---|
| 1 | \|ĝ₃\|Λ⁴ ≤ √(16L_A) = 119.75 / 128.17 / 138.54 / 145.90 / 168.06 at Δ_g = 1000 / 2000 / 5000 / 10⁴ / 10⁵; "120 / 128 / 139" (L:224, :253; S:1424; Table S4) | `data/nf02_psi.json` | `r5/consts.py`: mpmath adaptive quadrature of ψ₂/(p²+9/4Λ²) from the exact rationals (not the authors' moment recursion). K₀ = k₀ = 27.4850981709. Results 119.748 / 128.169 / 138.536 / 145.896 / 168.058, and 103.608 at 300. The fit 54.98√(ln(Λ/8.7415)+6.384/Λ) is within 0.01 at every gap. √(16ψ₂′(0)) = 54.976 | ✓ |
| 2 | τ₀ = 2.89·10⁻⁸ / 1.93·10⁻⁹ / 5.33·10⁻¹¹ / 3.50·10⁻¹² / 4.03·10⁻¹⁶ (S:1420-1421); 2L_AdS = 92.81 / 104.77 / 120.57; B_eff; (2b) 164.51 | same, plus `q3star_factored.json` | 2.8913·10⁻⁸ / 1.9283·10⁻⁹ / 5.3262·10⁻¹¹ / 3.5036·10⁻¹² / 4.0336·10⁻¹⁶. 2L_AdS = 92.8135 / 104.7694 / 120.5709 (72.0406 at 300). B_eff = −19.2647 (1000), drifting to −19.2781. 2ψ₂^Q′(0) = 17.240713. (2b) 164.511. ABJM 2ζ(3)/92.8135 = 2.590 % | ✓ |
| 3 | **Binding state at Δ_g = 2000:** +0.2004, budget 6.3·10⁻⁴, L = 320 (twist edge), 334 at the 1000 image; 78 at 1000 (S Table S13; GAP_SCALING) | `code/gapnear/gn_bind.py` on the cached K₂ rows | `r5/rerun_bind.py`: read-only in-process rerun of `gn_bind.run()` (NEARGAP2's `ng_bind.run_state` unchanged), behind a write guard. At 2000: λ = +0.20041; trunc 4.87·10⁻⁴ + interp 1.7·10⁻⁸ + low 1.98·10⁻⁵ + edge 1.20·10⁻⁴ = 6.27·10⁻⁴; **L = 319.8** (λ-level 564.6); 334.4 at m/Λ = 1.00499. At 1000: L = 77.7 (λ-level 130.2). Trunc ratio 0.255, low ratio 0.064: the Λ⁻² and Λ⁻⁴ laws. LR curvature 1.9·10⁻⁷ N_s | ✓. **The binding state is essentially flat**, so its large leverage is expected |
| 4 | **LR∞ certified uniformly:** "454 Arb Λ-balls on [998.5, 10⁶], sup r ≤ 0.9988" (S:1448-1449) | `code/gapscale/cache/lrball_part1.log`, `lrball_resume.txt` | Parsed both: 488 rows, of which **34 are NaN** (Λ = 1652–128005, in `part1`), superseded by the resume run. The **454 non-NaN balls cover [998.50, 10⁶] without gaps**. Max sup r(β ≥ 0.37) = 0.9988 on the first ball [998.5, 999.5], falling to 0.42 by Λ ≈ 1550. Beyond 10⁶: `gaplarge/cache/gl2_allgaps.txt` gives the scaling-lemma bounds; f_U > 0 on every β-cell ≥ 0.245, and f_U ≥ 0.99 A₁M₁ from β = 0.34 | ✓. Margin 0.12 % at Δ_g ≈ 1000. The NaN rows in the log should be explained in GAP_SCALING (the README line at :245 is terse) |
| 5 | J\*(m) = ⌈(0.37211m−1)/2⌉: 223 / 447 / 1119 / 2238 at the twist edges of 1000 / 2000 / 5000 / 10⁴; J\*/Λ → 0.22387 (S:1453-1454) | formula | `r5/jstar.py`: 2 sinh(0.185) = 0.372114; self-consistent twist-edge J\* = 223 / 447 / 1119 / 2238; asymptote 0.223864. Also 2 sinh(0.1225) = 0.245613 (b ≥ 0.2456Λ). Twist-edge masses in Table S13 (1.0020, 1.0025, …, 1.0144 at 2000; 1.0010, 1.0024 at 5000) all reproduce | ✓ |
| 6 | **Δ₁\* thresholds** 0.9·10⁵ (J = 0) … 4.2·10⁵ (J = 8); graviton 2.2 / 2.8–3.0·10⁵; target 0.9 G_cert/(Πμ²Λ₀²) (S:1607-1640) | `code/gapmorse/cache/gm_report.txt`, `gm_table.txt`; GAP_MORSE §0 | Table S15 equals the cache to the printed digit. rel\*: J = 5 gives 0.9·3.01·10⁴/(3·1.69·998.5²) = 5.36·10⁻³ (cache 0.00536); J = 0 gives 0.0354. μ = 1.3 is the conservative end of the sampled masses. **r_c**: the Supplement's "0.014–0.052" merges scalar (0.023–0.052) and graviton (0.014–0.032) radii (m8) | ✓ (not rerun: the certificates take CPU-hours) |
| 7 | **F_J of premise B_F:** 11.8 / 14.9 / 43 / 743 / 799 / 3.4·10⁴; J = 8: 0.85 (S Table S14) | GAP_UNIFORM §3.4 | G_cert/sup\|Γ\|: 2340/197.5 = 11.85; 3600/241.8 = 14.9; 1.10·10⁴/254.5 = 43.2; 3.01·10⁴/40.5 = 743; 6.34·10⁴/79.3 = 799; 1.99·10⁵/5.9 = 3.4·10⁴; 316/370.2 = 0.854 | ✓ |
| 8 | Lemma C reach 13 / 17 / 23 / 35 / 61 at 1000 / 2000 / 5000 / 2·10⁴ / 10⁵ (S:1445-1446) | `code/gapscale/cache/lcscan.txt` | m/Λ = 1.00, certified floor, worse helicity: min(13, 15), min(17, 21), min(23, 27), min(35, 39), min(61, 61). Matches. The scan's own fit says J_LC ∝ Λ^{1/3} | ✓; this growth law is relevant to M2.9 |
| 9 | Large-spin margins: W d_χ/H 0.489 / 0.246 / 0.100 (∝ Λ⁻¹); EDGE-R tolerance 0.406 / 3.45 / 57 (∝ Λ^{3.1}); "×8.5 / ×140"; L′ max d/(F/2) 0.190 / 0.035 / 0.0050 (≈ Λ^{−2.3}) (L:319-321; Table S16) | GAP_LARGESPIN §0, §4 | Exponents: ln(3.45/0.406)/ln(1998.5/998.5) = 3.09; ln(57/0.406)/ln(5.006) = 3.07; d_χ ratios 0.503 and 0.2045 against Λ⁻¹ ratios 0.4996 and 0.1998; L′ exponents 2.44 (2000) and 2.26 (5000) | ✓ (arithmetic of the quoted values). All [DV] given Δ_g-independent premise constants [M] |

Not rechecked (verified in report 4 and unchanged): 0.46Λ scalar threshold; c₁ 1 % thresholds 59.6 / 1306; 8πG = 256/c_T;
τ = 240|ĝ₃|; ∫|ψ₂| = 10.3860299071.

---

## 6. Physics checks

**6.1 The monotonicity remark** ("a theory with a larger gap also has a gap of 1000, so the 1000 numbers bound it"; L:250; S:1416-1419).
**Correct, and correctly used.**
* *(i) A1′ and A2 are inherited downward.* Single traces of spin ≥ 4 with twist ≥ Δ′_g also have twist ≥ 1000. The states with twist
  < 1000 are among those with twist < Δ′_g, whose data are tree-level EFT up to o(1/c_T) uniformly in spin.
* *(ii) (Pos) at 1000 is a property of the functional.* It concerns (Δ, J, parity), not the theory (S:236-238). It therefore applies to
  whatever states the larger-gap theory has at twist ≥ 1000, including its EFT double traces in [1000, Δ′_g).
* *(iii) Those double traces do not spoil the light side.* This is the point a careful referee will raise. They now sit on the heavy side
  of the 1000 sum rule, but that is already true in a theory whose gap is exactly 1000: [TT]_{n,ℓ} with twist ≥ 1000 are always
  "heavy" in this framework (L:115-116).
  * At k = 2, 3 they sit at double zeros, so they contribute at O(1/c_T²).
  * At k = 4 the light side is the arc of the EFT amplitude (L4 / Lemma A: D_LR = 0 on tree terms). The arc does not depend on where
    the light/heavy split is placed.
  * The quartic remark is therefore also consistent at leading order.
* *(iv) Scope of the remark.* It is an absolute-units statement. In the larger theory's own gap units it gives |ĝ₃|Λ′⁴ ≤
  120(Λ′/Λ₀)⁴, which is very weak. The Letter uses it in exactly this way and keeps it distinct from per-gap positivity (L:250-255).
  The Supplement grades it [D] and asks the authors to confirm it (CHANGES_v10 §6.9). I confirm it.
* *(v) One inheritance is worth stating.* The quartic conditions become conditions at Λ(1000), i.e. "no scalar below 0.46Λ(1000)".
  This is a weaker condition, and favourable.

**6.2 "Positivity at every gap ≥ 1000": logic and strength.**
* *Lemma T [P] (S:1463-1466; GAP_UNIFORM §3.1) is correct as a statement about the error-budget structure.*
  * In chart A the node coefficients, the flat parts and N_f are Λ-independent, and the Lemma-C bound decreases in Λ.
  * Hence the tolerance on sup|Γ| (in c₂ units) grows at least like Λ².
  * Positivity at a larger gap then needs only that the *relative* correction Γ/m² does not grow. The Letter's wording at L:314-315
    matches this weaker form, which is right.
  * Caveat 1: the realised k = 2 weight differs from ψ₂ by S(p)/Λ², so "exactly Λ-independent" holds for the flat weights. GAP_UNIFORM
    §3.1 states this caveat; S Table S12 does not. Add one clause.
  * Caveat 2: Lemma T is a statement about the *chain's* representation of the correction (nodes, windows). It is not a statement about
    the exact functional action. The window conventions are [M].
* *Premise B_F [NE] (S:1551-1554) is stronger than Lemma T needs.* It caps the c₂-unit correction by F_J times its 1000 value; Lemma T
  would allow growth by a further (Λ/Λ₀)². The Supplement says this, and "safety factor above ten" (L:317) refers to F_J ≥ 11.8.
  Correct.
* *Evidence for B_F.* At 2000 and 5000 the evidence is K₂-based, and in chart A the K₂ coefficients are structurally gap-independent
  (Prop. S). The real tests of B_F are:
  * the exact scalar to 10⁴ (drift ≤ 3.1 %), which is a convincing analogue;
  * the 1000 direct-minus-K₂ difference (≤ 13 c₂, 50–2000× below the allowance).

  The graviton test at 2000 is BATCH4 (M2.6). Physically the premise is very plausible: corrections approach flat space as 1/m² at fixed
  (J, x), with Λ-independent coefficients. The risk lies in non-uniformity: Stokes phenomena, which Lemma R addresses only at sketch
  level, and the light-cone edge.
* *Certified eventual positivity (S:1605-1618).*
  * The Morse-chart remainder lemma is a nice piece of work. It is [P] in structure, and its data are certified at five ξ-nodes.
  * It replaces the *truncation* part of the premise, not the premise: (iii) sup_m|Γ_{K₂,m}| bounded remains (M2.7).
  * It covers J ∈ {0, 2, 4, …, 8} only, with graviton data at J = 0, 2. The binding state's graviton amplitude (J = 5) is not certified.
  * It starts at ~10⁵, so 1000 ≤ Δ_g ≲ 10⁵ (two decades) rests on B_F.
* *Large spin at other gaps.*
  * The chart-B tolerances grow (Table S16), which is good news for EDGE-R (×8.5, ×140).
  * Everything is conditional on carrying the premise constants in AdS units [M]. The graviton edge integral is uncomputed at every gap.
  * The band of M3 is unassigned.
* *Net, my reading.* The right statement is the Supplement's: "positivity of Ω^EE at Δ_g > 1000 has at best the grade it has at 1000,
  plus the no-growth premise below Δ₁\*" (S:1411-1412). Its weakest links are [NE] (B_F) and [M] (H_mono-A, H_mono-B, the carried
  constants). The Letter should use that sentence's strength (M2).

**6.3 Normalisations.** Unchanged, and correct:
* 8πG = 256/c_T;
* τ = 240|ĝ₃| (t₄ half [DV], α_T half [D], with only the modulus used);
* Λ = √((Δ_g−2)(Δ_g−1)), not a cosmological constant (L:117);
* the collider disk τ ≤ 4;
* Λ⁻² scalings: the budget ratio 0.25 ± 0.01 against (998.5/1998.5)² = 0.2496, and 0.040–0.043 against 0.0399 at 5000. Checked.

**6.4 Order of limits.** Unchanged and correct. The sum rule is used at finite c_T with J₀ ≤ 1 and margins ½ and 3/2, and (Pos) concerns
(Δ, J). The new material respects this, since c_T → ∞ is taken at each fixed Δ_g. Two remarks:
* Uniformity in the gap of the o(1/c_T) light remainder is not needed, because each gap is a separate statement.
* "Bounds at every gap" are therefore pointwise in Δ_g, not uniform in a double limit. That is fine, but the Letter should not suggest
  more.

**6.5 A1′ and ABJM.** Unchanged and sound at the level claimed.
* In the IIA regime KK modes have J ≤ 2, and the string scale sets the spin-4 gap.
* The scalar condition holds at the orders computed, at sketch level.
* The monotonicity remark now makes "a bound ABJM theory obeys" apply to every IIA-regime ABJM theory with string-scale gap ≥ 1000. This
  needs λ ≳ 3·10⁹ in the convention ML ≈ Δ_g (S:1843-1844). Placement remains the issue (m4).

**6.6 c₁.** Unchanged: universal [D], c₁ = 0 given F₂ [A], value [I].
* At larger gaps its effect shrinks: x/|c₁| = 3.9·10⁻⁵ at 2000 (Table S4).
* The abstract's last sentence is accurate.
* It is the cheapest substantive improvement left (one spinning Witten diagram).

**6.7 H_EFT(b).** Unchanged. It is the standard CHMRS/CLS-D assumption [A], stated at L:199-201. The new remark does not affect it.

**6.8 κ > 0.** κ is needed only if τ ≠ 0. Its sign is [D-sketch + M]. It matters: with |ĝ₃|Λ⁴ up to 120, κ|ĝ₃|²Λ⁸ ≈ 0.11·1.4·10⁴ ≫ 92.8.
The Letter states the condition in (2a) (L:226); the abstract omits it (m4).

**6.9 Large-spin premises at 1000 (AG_U, EDGE-R).** Unchanged.
* (AG_U) is an assumption beyond all orders. The plane-wave representation is shown *not* to supply it at the edge for J ≤ 50.
* EDGE-R for the graviton is uncomputed. The scalar analogue gives 0.405 against 0.406 at (26, 1.026Λ).
* The Letter now says both plainly (L:290-296).
* v10's upgrade of P0 to rigorous at b ≥ 114.06 removes one [NE] premise at the sampled states with J ≥ 80 and at (64, 1.062Λ). That is
  real but small.

**6.10 Functional definition beyond Λ.** The v8 fix stands: continuum weights times θ(Λ+100, 10); grids are tabulations (S:198-203). At
other gaps the functional is not published (M5).

**6.11 AI-use disclosure.** Placement in the acknowledgments is appropriate for APS journals; it has no place in the main text. Its
content is narrower than the record and conflicts with the Supplement's "unreviewed" labels (M6). Add the same disclosure to the
Supplement.

---

## 7. What still stands between this and PRL acceptance (ranked by impact ÷ effort)

1. **Fix the all-gap wording (M2), the Fig. 2 caption, the Fig. 1 markers (m9), and the abstract's ABJM and κ (m4).** *Effort: half a
   day of writing.*
   * It costs ≈ 20–40 words, paid by cutting "Large spin" to about 200 words.
   * Necessary for any venue. It does not change the PRL verdict.
2. **Fix the AI disclosure and the "reviewed / unreviewed" tension (M6). Strip source comments (m13).** *Effort: 1 hour.*
3. **Publish the 2000 and 5000 weights; state the construction and grade of Ω^EE at general Λ (M5).** *Effort: half a day.*
4. **Fold in BATCH4 and SEAM, and fill the M3 band.**
   * *Effort, BATCH4:* ≈ 13–35 CPU-h; 1–2 h on 16 workers.
   * *Effort, SEAM:* in progress.
   * *Effort, M3 band:* K₂ rows for 33 ≤ J ≤ 51 at 2000 with Δ_ref ≥ 5·10⁴, a few CPU-h per spin.

   This makes "we repeated the checks at 2000" true as written, and gives the first graviton test of B_F.
5. **Ω^Q_up at 2000 (binding states, K₂ rows).** *Effort: days.* It lets the quartic bound join the all-gap statement, or justifies its
   absence.
6. **Compute c₁.** *Effort: 1–3 weeks; one spinning AdS₄ tree diagram with a W³ vertex.* It removes the abstract's last caveat.
7. **Compute the graviton EDGE-R at the Theorem W states.** *Effort: ≈ 90 CPU-h external, plus fixing the k = 3 (C3H) bug, days to
   weeks.* Turns Theorem W's only [M] premise into data. At larger gaps the tolerance grows like Λ^{3.1}, so one computation at 1000
   settles the most.
8. **Error-control the rest of the near-gap region at 1000.** *Effort: ≈ 50–100 CPU-h.* Covers J = 9–11, 13–15, 17–20, 22, 24, 25 to
   1.3Λ, the twist edges of 23 and 31, and J = 29 sign-free. Together with item 7, this is what M1 needs. It also carries to every larger
   gap through Lemma T.
9. **Close B_F on 1000 ≤ Δ_g < Δ₁\*.** *Effort: weeks.* Options: higher Laplace orders in ball arithmetic, a degree-12–16 C/Arb Taylor
   model (GAP_MORSE "What remains" 1), or K₂ plus direct rows on a gap ladder. It turns the J ≤ 7 transfer from "[P given B_F]" into a
   certificate.
10. **External review of S5 (c₁, H_K, the arc, E′) and S7–S8; rewrite the Supplement as derivations without version history.** *Effort:
    1–2 weeks plus one external reader.*
11. **A proof of (AG_U), or a large-spin argument without it.** *Effort: high, possibly out of reach.*

**Realistic path.** Do items 1–3 now (a day), then submit a long paper (PRD or JHEP) with items 4–6 as they land. With v10, the
one-gap objection is no longer what blocks a PRL. PRL becomes plausible if items 7 and 8 succeed: positivity at error-control level on
the whole near-gap region at 1000, and the graviton edge integral computed. Since Lemma T carries both to every larger gap, the
"every gap ≥ 1000" claim would then rest on a single explicit premise (B_F). Item 6 (c₁) should be done in addition.
