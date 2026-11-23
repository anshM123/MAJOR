# Second referee report

**Manuscript:** "Conditional bounds on R³ and one R⁴ graviton coupling in AdS₄/CFT₃ from ⟨TTTT⟩ light-ray sum rules"
(`ads4/prl/letter.tex`, `supplement.tex`, `macros.tex`, `refs_verified.bib`; build in `prl/build/`; switches `\HLJclosedfalse`
(headline Δ_g = 1000) and `\HSdrafttrue` (HIGHSPIN placeholders visible))

*Written 23 Nov 2026 as an internal pre-submission review. No existing file was edited and nothing was committed. File references are
`path:line`, relative to `ads4/`. I formed the view in §§1–3 before reading the first report (`prl/REFEREE_REPORT.md`), then checked its
points against the current text (§4).*

**What I checked myself.**
* **Headline constants.** I recomputed them from the published data (`prl/data/nf02_psi.json`, `q3_full_delta0.2.json`), in mpmath at
  40 digits. Script: `/tmp/claude-0/-home-user-MAJOR/60eea3a7-07d6-5453-af55-c4346d86f9b5/scratchpad/ref2_spot.py`.
* **Independent reproduction script.** I ran a copy of `code/repro/bounds.py`, redirected to the scratchpad so that the cache in the repo
  was not touched.
* **Lemma W.** I verified its algebraic identity (supplement:292-298) by hand.
* **Literature.** I checked titles, authors and abstracts of nine references against the arXiv API: 2311.04271, 2106.10274, 2102.08951,
  1805.07393, 2104.11238, 1904.05905, 1904.05892, 2201.06602, 2205.01495. I also checked three candidate additions (§3, m13).
* **Length.** I counted it with the authors' tool and by rendered layout.

---

## 1. Summary and recommendation

**Summary.** The Letter builds three fixed light-ray functionals on ⟨TTTT⟩ in a large-c_T CFT₃ with a higher-spin gap Δ_g (A1′), at tree
level (A2):
* Ω^EE, at k = 2, 3, with a new dual structure "EE";
* Ω^Q_up and Ω^Q_low, de-shifted k = 4, 5 pieces.

It proves an **implication**, graded **[D-sketch]** overall: per-state positivity (Pos) of these functionals on the heavy spectrum
(and on light spin-0, 2 states) gives the following.
* √(t₄²+α_T²) ≤ 3.13·10⁻⁶ at Δ_g = 300 and 2.89·10⁻⁸ at Δ_g = 1000. Equivalently |ĝ₃|Λ⁴ ≤ 104 / 120.
* −O(Λ⁻²) ≤ g̃₄Λ⁶/8πG ≤ 72.04 / 92.81 if no scalar lighter than 0.46Λ couples to TT, and 125.7 / 164.5 otherwise. This holds exactly
  for every value of an uncomputed universal constant c₁. The identification g̃₄ = g₄ is an assumption, [A].

**Arithmetic.** Every constant I traced reproduces from the published weights (§5). The bookkeeping of what is proved, derived,
numerical and assumed is unusually careful. Since the first report, the following were added:
* the cutoff/H_K argument (now [D]);
* certified flat floors;
* published continuum weights;
* an exact treatment of c₁;
* an honest A1′ under which ABJM-IIA is in scope;
* a much better positioning paragraph.

**What is not there.** The physical content is still positivity, and positivity is still not established at either value of Δ_g.
* **Δ_g = 1000.** Only a handful of states are error-budgeted [NE±]: the binding J = 5 state, near-gap J ≤ 8, and a strip m ≥ 1.3Λ,
  b < 5. Near-gap J ≥ 12, large spin and large impact parameter rest on model hypotheses. The sign-free form of one of them is
  violated at (1100, 100). One computed state (J = 29 at the edge) has a sign that depends on how the edge window is modelled.
* **Δ_g = 300.** A covering hypothesis is open, [I].
* **The lower functional Ω^Q_low** has no positivity evidence at all beyond leading Regge.
* **Choice of functional.** The functional and its gauge were chosen after failures had been seen. In the canonical gauge a computed
  state fails badly.

**Status of the derivation.** The implication itself is [D-sketch]. Most of its links come from internal analyses written in the last two
days and not yet reviewed. The Supplement gives sketches, not derivations a referee can check.

**PRL criteria.**
* *Importance.* Explicit-constant AdS graviton bounds from ⟨TTTT⟩ are a recognised open target (CLS-D).
* *Strength of the bounds.* They are weak: about 10× the transplanted flat optimum in |ĝ₃| (about 100× in |ĝ₃|²), and an IR-log
  dominated g₄ window in which the only known theory sits at 2.6 %.
* *Breadth of interest.* Limited by the ledger-style presentation.
* *Are the central claims established at the level stated?* The conditional implication is, at [D-sketch]. The numbers that make it
  interesting are not, because they depend on positivity.

**Recommendation: reject for PRL in the present form, with transfer to PRD (or JHEP) as a long paper strongly encouraged.** This is not a
minor-revision case. A PRL version is realistic only if (Pos) for Ω^EE is closed at one Δ_g across the whole heavy spectrum, at least at
the [NE±] level, with an analytic large-spin and large-impact-parameter argument. The Letter would then also need to be cut down to that
clean statement and rewritten for a general audience (§6, items 1–2). If HIGHSPIN comes back weak (M1(d)), the item-(1) headline at
1000 loses its support and the Letter's structure does not survive.

---

## 2. Major issues

### M1. Positivity, the physical content, is not established, and the main text understates its weakest points

**Where.**
* letter.tex:84-89 ("What is established");
* letter.tex:304-321 (Table I);
* letter.tex:323-355 (Positivity);
* supplement.tex:602-730;
* UT3 §5.2–§5.6.

**Claim.** "Positivity itself is not established" (letter.tex:87). The Letter states this honestly, but the evidence is presented as
stronger than it is.
* **(a) Only six regions are error-budgeted at 1000.** These are the J = 5 binding state, near-gap J = 0, 2, 4, 6, 8, and the strip
  m ≥ 1.3Λ, b < 5. J = 8 is error-budgeted only on the λ-level budget; its sign-free leverage is 1.1–1.8 (UT3 §5.2; NG2 §0).
  * Near-gap J ≥ 12 is [NE + M], and the budget does not close (L ≈ 1–3).
  * Near-gap large spin rests on H_LJ(1000), [D-sketch + M]. Its sign-free form is *violated* at (1100, 100) (−2154 ± 43 against a
    tolerance of 1409). It "holds" only after a refinement that removes a leading-Regge-shaped class.
  * Large impact parameter (b ≥ 5) rests on H_R′, H_bulk and the hybrid model, all [M].
  * There is no mass above which all spins are covered (letter.tex:334-335).
* **(b) One computed state has a model-dependent sign, and the main text does not say so.** The J = 29 twist-edge state at 1000 comes out
  −3.62 in the plain chain and +1.39 ± 0.03 with direct nodes (NEARGAP2.md:42-43, :230-235). NEARGAP2 calls it "representation-dominated"
  and says the same may hold for other odd J ≳ 27. A value above 1 in units of N_s means the AdS correction is larger than the flat value.
  The main text mentions this state only in a TeX comment (letter.tex:4). Table S6 (supplement.tex:664) says "edge window unresolved".
* **(c) Gauge sensitivity and in-sample selection.** In the canonical Abel gauge the J = 24 state at 300 fails, at n = 1.55 by spline and
  −22 to −40 from direct data (supplement.tex:617-621). So the AdS corrections the functional must survive are not uniformly small. They
  are O(1) and gauge-dependent for some states, and the gauge was chosen after the failure was seen. This is disclosed (letter.tex:338-341),
  but it undercuts the implicit message that "AdS corrections are perturbations of a certified flat margin".
* **(d) HIGHSPIN.** The \HS items are near-gap J ≥ 12, the J = 29 edge and H_LJ(1000). If they come back weak, the outcome depends on how
  weak.
  * *Budget still does not close (status quo, [NE + M]).* The claims stand formally, because the theorem assumes (Pos). The Δ_g = 1000
    headline then rests on model-only evidence for every near-gap state with J ≥ 12. That is the same kind of situation the authors use to
    demote 300.
  * *Any state negative (e.g. J = 29 resolves negative, or the refined H_LJ fails).* (Pos) for Ω^EE is then false at 1000, and the
    following collapse:
    * item (1) at 1000;
    * the 164.5 form of item (2), which needs (Pos) for Ω^EE (letter.tex:249-250);
    * the Δ_g ≥ 300 corollary;
    * the abstract's 2.89·10⁻⁸ headline.

    Only item (2a) at 1000 would survive. Its own large-mass evidence is [NE±] only on 1.3Λ ≤ m ≤ 3Λ and [M] elsewhere (UT3 §5.4). Fixing
    the failure by re-tuning Ω^EE would add a further round of post-hoc selection.
  * *So the Letter's architecture is fragile to the pending computation.* The authors should write down a fallback before submission.

**Why it matters.** The comparable literature presents bounds whose functionals are checked across the full heavy spectrum, including
large-spin asymptotics. CHMRS 2106.10274 is the example, with scalar blocks, and its abstract says "we prove". Here the comparable check is
incomplete at both Δ_g, with known stress points. A reader of the abstract ("tested numerically") will not infer this.

**Fix.**
1. Close (Pos) for Ω^EE at Δ_g = 1000. That means near-gap J ≥ 12 at [NE±], and an analytic large-spin / large-b argument (e.g.
   light-cone asymptotics of the de-shifted kernel plus a certified compact interval).
2. Until then, put the J = 29 status and the canonical-gauge failure in the main text, and say in the abstract which regions are only
   modelled.
3. Prepare the fallback statement for the case that HIGHSPIN fails.

**Fatal for PRL?** Yes, in the present form.

### M2. The lower end has neither a constant nor any positivity evidence, yet it is advertised as a window

**Where.**
* letter.tex:36-37 (abstract);
* letter.tex:155-157 (definition of Ω^Q_low);
* letter.tex:238-251 (theorem);
* letter.tex:405-407 (Discussion: "the window between −O(Λ⁻²) and …");
* Table I (letter.tex:304-321);
* supplement.tex:132-134 and Table S5.

**Problem.**
* **No positivity evidence.** Item (2)'s lower end needs (Pos) for Ω^Q_low. Neither the Letter nor the Supplement contains any evidence for
  it beyond leading Regge. Table I and Table S5 report only Ω^Q_up. I found none in UT3 §5.4 either. RESULT_Q_STATEMENT.md:269-270 gives
  only the leading-Regge image P_{−½}(cosh β) > 0 and the flat image.
* **The margin vanishes where corrections are largest.** The leading-Regge margin P_{−½}(cosh β) decays like βe^{−β/2} at large β (large
  J at fixed m). So the lower functional has its smallest margin exactly where beyond-Regge corrections are least controlled.
* **No constant.** The lower end itself is "−O(Λ⁻²) × (EFT data)" with no computed constant (UT2 §3). In the R ≠ 0 case it also needs
  c_L > 0, graded [D-sketch + M]. The c_L slope (0.255 or 0.49) was not reproduced (supplement.tex:786).

**Why it matters.** The abstract and Discussion present "−O(Λ⁻²) ≤ g₄Λ⁶/8πG ≤ 92.8" as a two-sided result. What is established is a
conditional "g̃₄ ≳ 0" whose positivity input is untested.

**Fix.**
* Either give Ω^Q_low evidence at the binding states and at large β (at least the same [NE + M] level as Ω^Q_up),
* or drop the lower end from the abstract and Discussion and keep it in the theorem with an explicit "(Pos) for Ω^Q_low: no evidence
  beyond leading Regge".

**Fatal?** No, but it must be fixed.

### M3. The implication rests on recent, unreviewed internal derivations that the Supplement only sketches

**Where.**
* supplement.tex:50-52 (Provenance: "none of these is externally refereed");
* NOTES.md §6 item 6 ("every analysis behind v3 except REPRO / REPRO_ENGINE is unreviewed");
* letter.tex:270-301 (Inputs).

**Problem.** In this paper "[D]" means "argument written out, *not independently checked*" (UT3 §0). The load-bearing links upgraded in
v3 come from deliverables dated 22–23 Nov:
* H_K [D], which depends on the exact EE-kernel reduction to nine b-shifted scalar partial waves, f_pos(±3i/2) = 0, and the
  J₀ < 3/2, 5/2 criteria (EESTEP_CLOSE, HK_RESIDUAL);
* L6′ [D], which depends on Lemma Y, "exhaustive over 997 836 structures" (R4ARC_HEFT);
* the ASD closure and universality of c₁ (C1_VALUE);
* the R⁴ calibration −𝒜₄[R⁴] = 2g₄ (F_NORM, R4ARC_HEFT);
* H_RS via "D′ is CLS-D's causal domain" (PCYC_CLOSE).

The Supplement gives each as one paragraph. Lemma W (supplement.tex:292-298) is complete, and I checked its identity. Lemma C′ and the
calibration are readable. The H_K corner count, the EE reduction, Lemma Y and the ASD closure cannot be checked from the text. Only the
light-side constants and Ω^Q's s_d have been reproduced independently. The beyond-Regge correction at Ω^EE's binding state and the Ω^Q
scan have not (letter.tex:357-361).

**Why it matters.** The paper's selling point is grade discipline. A grade a referee cannot audit is a claim, not a result. The overall
[D-sketch] is then a lower bound on the uncertainty, not an estimate of it.

**Fix.**
* Write the full derivations of H_K (corner count plus EE reduction), L6′, the calibration and the c₁ structure, in a long companion paper
  or a much longer Supplement.
* Give, for each [DV]/[certified] claim, the script that checks it.
* Have the H_K EE reduction and Lemma Y independently reproduced, as was done for the constants.

**Fatal for PRL?** Not by itself, but together with M1 it means the result is not "established at the level stated" in the sense a PRL
referee needs.

### M4. Grade slips: places where the Letter claims more than UT3, or than its own weakest-link rule, allows

The paper asks to be judged on its grades, so these matter.
* **(a) Universality of c₁, letter.tex:213-214 ("[D]").** UT3 L3d grades the ASD closure "[D] (one input, a regular second-order ASD
  solution, [D-sketch])". By the paper's own rule (letter.tex:94) that is [D-sketch]. The same slip is in supplement.tex:192 (L3,
  "universality [D]"). UT3 copies it from C1_VALUE.md:26.
* **(b) Arc statement, letter.tex:179-180.** "depends only on X's top Regge coefficient [D; DV]" drops the qualifier "[D given C′ + Mellin
  convergence]" (UT3 L4; supplement.tex:193 has it).
* **(c) Integrability of the R⁴ light function, letter.tex:191-195.** The Letter says "[D]" and omits that the physical Δφ scan gives local
  slopes 2.5–2.8 at Δφ = 0.04. The analytic count predicts 2, and the criterion is < 3.5. It also omits that the AdS₄ R⁴ ⟨TTTT⟩ itself was
  never evaluated (UT3 L6′; supplement.tex:282-284). The numerics do not yet confirm the count. Say so in the Letter.
* **(d) H_EFT(a), letter.tex:296-298.** "pointwise in angle for all dispersive (holographic) tree terms [D; DV]". An unexplained ≈ −3 %
  leading-order ratio for a spin-2 heavy exchange remains [NE] (UT3 L7(a); supplement.tex:345-346). The calibration uses only a scalar,
  so the R⁴ normalisation is not directly affected. But "[D; DV] for all dispersive terms" is contradicted by the one spinning case
  checked. Restrict the claim, or explain the −3 %.
* **(e) The word "bound", letter.tex:88 and :329-331.**
  * The text says "error-bounded numerically" and "we bound every AdS correction to it with factor-3 safety margins … λ_min ≥ +0.198".
  * [NE±] is a ×3 safety factor on *model* terms (K₂ Laplace rows, the chain, the bulk-ν envelope). It is not a bound.
  * What is genuinely bounded is the flat floor minus Lemma C (≥ +0.19993 N_s). The beyond-Regge part is an estimate.
  * Replace "bound"/"error-bounded" by "error-budgeted", and say what the budget covers.
* **(f) "Theorem" label, letter.tex:223.** It is attached to an implication whose net grade is [D-sketch] (letter.tex:262). "What is
  established" (letter.tex:84-86) says the implication "is derived". Use "Proposition (conditional, [D-sketch])", or equivalent, and say
  "derived at sketch level".

**Fix.** Correct each item. None of them changes a number.

**Fatal?** No. The pattern (a), (c), (d), (e) is worth fixing before any submission, because a referee who finds one will distrust the rest.

### M5. Importance and positioning: the bounds are weak, and the closest explicit-constant precedent is mis-cited

**Where.**
* letter.tex:66-67 ("Flat four-dimensional graviton bounds with an infrared logarithm are in Refs. [CHLPS, CHLPSpw]");
* letter.tex:408-413;
* letter.tex:415-422.

**Problem 1: mis-citation and a missing precedent.** CHLPSpw (2205.01495) is the D ≥ 5 paper. Its abstract says "in spacetime dimensions
D ≥ 5 where the S-matrix does not suffer from infrared ambiguities". It also "record[s] new bounds on the central charges of holographic
conformal theories": explicit-constant, gap-dependent bounds on the higher-d analogue of t₄ (DICTIONARY_T4.md:54, :237 quotes
|a−c|/c ≤ 23/Δ_gap² in d = 4, uplift asserted). CHLPS (2201.06602) remarks that its bound translates to holographic CFT₃. These are the
closest precedents for item (1). The novelty sentence (letter.tex:415-416) is technically restricted to d = 3 and to "derived from
⟨TTTT⟩", but it must name them.

**Problem 2: strength.** The Discussion writes |ĝ₃|²Λ⁸ ≤ 3022.34 ln(2Λ/3) − 5327.21 and then says item (1) is "about ten times weaker" than
the transplanted CHLPS optimum.
* The factor 10 is in |ĝ₃|, i.e. in R.
* In |ĝ₃|², the quantity just displayed, it is about 100×.
* The log coefficient is 121× CHLPS's 24.9 (RESULT_D.md:157).
* For g₄, the window [0, 92.8] is set by the 4d IR log (17.24 ln(2Λ/3)). The one theory in scope sits at 2.6 %, which the authors call a
  check that "could only fail grossly" (UT3 §8).

**Why it matters.** PRL asks for importance. What is new is explicit constants, and those constants are loose. The paper should own that
quantitatively. It should not imply that the cubic result is near-optimal.

**Fix.**
* Correct the CHLPSpw description.
* Add one sentence on CHLPSpw's holographic central-charge bounds and on CHLPS's CFT₃ remark.
* State the comparison as "≈ 10× in |ĝ₃| (≈ 10² in |ĝ₃|²)".
* Say plainly that the g₄ window is IR-log dominated.

**Fatal?** Contributes to the PRL verdict; easy to fix.

### M6. Not written for a broad audience

**Where.** Throughout the main text, especially letter.tex:84-94, :144-219, :262-301, :323-361.

**Problem.** Roughly 40 % of the Letter is grades, inputs and evidence bookkeeping. A general physics reader cannot decode the
vocabulary: "de-shift", "Regge arc", D_LR, "dual-line flip", "light-ray one-point matrices", "binding pair", s_d, s_fail, "twist edge",
"refined form", H_LJ, [NE±], "N_s", "Mehler–Dirichlet curvature bound", "celestial kernel". There is no figure. The physical idea is never
explained in plain words: why commuting light-ray operators give a positive heavy side, what the de-shift buys, and why the R⁴ term sits at
a Regge "borderline".

**Why it matters.** Broad interest is a PRL criterion, and the current text reads as an internal status report.

**Fix.**
* Move the Inputs list (letter.tex:270-301), the Reproduction paragraph (:357-361) and the c₁-threshold details (:253-260) to the
  Supplement.
* Spend the space on one paragraph of mechanism and one figure. Two good candidates:
  * the bound versus Δ_g, with the ABJM point and the shaded "values of the implication only" region;
  * a (J, m/Λ) map of the positivity evidence, coloured by grade.

**Fatal for PRL?** In combination with M1, yes.

### M7. Item (2)'s headline needs a physical condition that A1′ does not supply, and its repair couples the two items

**Where.** letter.tex:37-38 (abstract), :201-210, :238-250, :253-256; supplement.tex:555-567.

**Problem.**
* **The headline functional fails on a generic class of states.** Ω^Q_up is *negative* on light scalars with m < 0.4602Λ [P]. Light moduli,
  the dilaton or axions coupling through φW² or φWW̃ are generic in string compactifications. So the 92.8 headline needs "no scalar below
  0.46Λ couples to TT". That condition is physical, not implied by A1′.
* **The abstract misstates the condition.** It says "if no light scalar couples to TT" (letter.tex:37). The actual condition is on
  scalars below 0.46Λ.
* **The repair has a cost.**
  * The repair (Ω^Q_up + 0.04 Ω^EE) is [D-sketch].
  * Its flat margin on a massless scalar is only 0.022N.
  * At ε* = 0.038 the sign is [I], and ε = 0.040 is chosen to clear it.
  * Positivity on light states with Δ_L = O(1) beyond leading Regge is itself [D-sketch] (Pos_light, supplement.tex:691).
  * The repaired form needs (Pos) for Ω^EE, so the 164.5 headline inherits all of M1.
* **The quartic functional is not optimised for this.** It was not built with light-state positivity imposed; the authors list this as a
  next step (letter.tex:427-428).

**Fix.**
* State the condition precisely in the abstract.
* Say in the theorem that the 164.5 form depends on (Pos) for Ω^EE.
* Ideally, re-optimise Q with light-scalar positivity imposed.

**Fatal?** No.

### M8. The Δ_g ≥ 300 corollary and "explicit dependence on Λ" oversell values of the implication

**Where.** letter.tex:38-39 (abstract: "with explicit dependence on Λ"), :364-383; supplement.tex:732-751.

**Problem.**
* **The Λ-forms are not bounds.** They are values of the implication with the functionals fixed. Positivity evidence exists only at 300
  and 1000.
* **The corollary leans on an open input and on its own conclusion.** It assumes (Pos) at 300, *including* the window that is [I]
  (letter.tex:379). It also assumes H_LJ(Δ_g) for *every* Δ_g ≥ 300, which is nearly the conclusion. Its support is a proxy analysis.
* **The abstract invites the wrong reading.** Advertising "explicit dependence on Λ" leads readers to treat Eqs. (4)–(5) as bounds for all
  Δ_g.

**Fix.**
* Drop "explicit dependence on Λ" from the abstract, or qualify it.
* Move the corollary to the Supplement as a remark.

**Fatal?** No.

---

## 3. Minor issues

**Length and format**

* **m1. Length.**
  * *Authors' count.* `NOTES.md` §2 reports ≈ 3080 words by its tool. I reproduce that: 2779 prose words (inline math = 1 word), 9 display
    rows (144) and Table I (155), total 3078.
  * *The tool undercounts.* Table I renders as about 18 lines, not 9 rows, because every row wraps in the single column; that is about
    +100. Eq. (2) and Eq. (3) occupy more display rows than counted. The text has 374 inline formulas (about 3100 rendered math characters).
    Counting long formulas as several words adds roughly 300–400.
  * *Rendered length.* Text ends at about 4.1–4.3 REVTeX `prl,twocolumn` pages, including the full-width title and abstract block. That is
    at, not below, the four-page equivalent of 3750 words.
  * **Verdict.** Nominally compliant under the literal APS rule (≈ 3200–3600 by reasonable conventions), with **no headroom**. The pending
    HIGHSPIN text, any requested definitions (M6) and a figure must be paid for by cuts. I recommend cutting about 500 words (the Inputs
    list, Reproduction and c₁-threshold material). The real problem is density, not raw count.
  * *Abstract.* About 530–560 rendered characters, under the 600 limit.
* **m2. Front matter.** Affiliations, emails and acknowledgments are still TODO (letter.tex:22-26, :431; supplement.tex:26, :28). "TODO:
  email" prints in the PDF footnote.
* **m3. References print without years for every JHEP entry.** For example [1] reads "J. High Energy Phys. 10, 079". apsrev4-2 drops the year
  for these entries (build/letter.bbl:63-64). Use `volume = {2009}`, `number = {10}` (JHEP's convention), or otherwise force the year.
  BibTeX also warns "missing number" for all JHEP entries.

**Abstract and statement**

* **m4. Abstract wording.**
  * It bounds "g₄", but the theorem bounds g̃₄ and g̃₄ = g₄ needs H_EFT(b) [A]. Say "the CFT-defined R⁴ coefficient", or add "(given a
    standard EFT assumption)".
  * (Pos) is also assumed on light spin-0, 2 states.
  * "if no light scalar couples to TT" should read "if no scalar lighter than 0.46Λ couples to TT" (M7).
* **m5. Notation clashes.**
  * R := √(t₄²+α_T²) (letter.tex:116) clashes with R³/R⁴ in the title and text. Rename it, e.g. τ or |t|.
  * k is used for Regge spin, for k(0) and for k₀; K for the light coefficient.
  * c is used for c = ½ in Q, c(ν), c₁, c_L, c_T and c₂.
  * b is used for the impact parameter and for the Abel coefficient.
  * Λ (gap scale) clashes with the cosmological constant, which is harmless but worth a footnote.
* **m6. Grade legend.** "[D; DV]", "[P; DV]", "[Q + DV]", "[D-sketch + M]", "[D-sketch + NE]" and "[D, one step D-sketch]" are used but not
  defined (legend at letter.tex:89-94). State the composition rule once.
* **m7. Table I (letter.tex:304-321).**
  * The caption says "[NE + M] unless stated", but the "light J = 0, 2" row is [P given LR form] + [NE] (UT3 §5.5).
  * The J = 8 entry needs "λ-level budget only".
  * The "large mass" [NE±] row depends on the bulk-ν envelope H_bulk [M], tolerated up to ×76. Say so.
* **m8. "≥ 19×" (letter.tex:350; supplement.tex:701).** 11/0.59 = 18.6, so write "≈ 19×" or "≥ 18×".
* **m9. A1′ wording (letter.tex:102-105).** "(at tree level only T × T enters ⟨TTTT⟩)" is confusing, since only T × T ever enters the
  s-channel. Say instead that spins 1, 3 have no ⟨TTO_J⟩ structure.
* **m10. Parity violation.** The theorem covers parity-violating theories (letter.tex:98). The Letter and Supplement should state, with a
  grade, why an indefinite-parity heavy operator gives no ++/OO interference in Ω^EE and Ω^Q: the functionals are parity-invariant, so the
  action is even in λ_O. PARITY_ODD_HEAVY.md §0 argues this for older functionals and leaves the odd-coupling normalisation [I]. Say whether
  this carries over to EE / de-shifted Q.
* **m11. Natural size of c₁ (letter.tex:257-258, "10⁻³–10⁻¹ [M]").**
  * The estimate comes from C1_COMPUTE §3.5. That section's picture (an independent seagull, by Maldacena–Pimentel analogy) was corrected by
    C1_VALUE (UT3 N11).
  * The estimate itself takes "O(1) read generously as 30–60" (C1_COMPUTE.md, §3.5).
  * Either re-derive it from the energy-channel picture or drop it. The thresholds (|c₁| < 4.6 / 60, 99 / 1262) are the useful content.

**ABJM and literature**

* **m12. ABJM paragraph (letter.tex:386-398).**
  * *A1′ holds without the extra claim.* A1′ for IIA-ABJM needs only that supergravity fields have spin ≤ 2. The [D-sketch] "no light single
    traces in T × T at tree level" is a separate statement. It is what removes the LS rows and satisfies item (2a)'s scalar condition. The
    sentence chains them as if A1′ needed it.
  * *λ convention.* λ ≈ 2.6·10⁷ / 3.2·10⁹ (supplement.tex:811) corresponds to ML = Δ_g. With Λ = ML, as the Letter states, it is
    2.51·10⁷ / 3.15·10⁹ (my run of `repro/bounds.py`). Pick one convention.
  * *Citations.* CPY is the 𝒩 = 8 M-theory paper, so cite it for IIA only alongside BCP. Give the equation or figure in CHMRSswamp for the
    ten-dimensional window "[0, 3.000]". UT3 §11 #20 notes that equation numbers are still owed.
* **m13. Possibly missing references** (verified to exist, with titles and authors confirmed on the arXiv API).
  * Córdova–Maldacena–Turiaci, arXiv:1710.03199 (interference bounds on ⟨TTO⟩ in the conformal collider), relevant to the light-scalar
    discussion.
  * Häring–Zhiboedov, arXiv:2202.08280 (gravitational Regge bounds), relevant to I1.
  * Hartman–Kundu–Tajdini, arXiv:1610.05308 (ANEC from causality).

**Smaller grade and wording points**

* **m14. Quartic form in the Discussion (letter.tex:406).** "17.24 ln(Λ/4.59) + O(1/Λ)" is the truncated form, which is 0.06 below the true
  upper end at 300, i.e. on the unsafe side. That is acceptable with the "+O(1/Λ)", but add "(asymptotic)".
* **m15. Flat certificate.** Its corner lemma (Lemma B) code is unreviewed (UT3 §5.1, "new corner code unreviewed"). Mention this next to
  "[certified]" (letter.tex:310, :325).
* **m16. Corollary wording (letter.tex:380-381).** "if at fixed (J, m/Λ) the beyond-Regge corrections grow at most like m … [M]". At J = 0
  this is [D given RD] + [NE] (UT3 §7). Fine; say "(at J = 0: derived)".
* **m17. Table S6 (supplement.tex:659-660).** "plausibly ≤ 0.80 at large mass" is not a grade. Write "[M]".
* **m18. Regge-intercept caveat (letter.tex:273-277).** "with the caveat of their Sec. 4.2.6": say in one clause what the caveat is and why
  it lies outside the Regge box.

**Published data**

* **m19. Data hygiene (`prl/data/`).**
  * *Stale ε_W in the README.* `README.md` §3 still quotes ε_W = 1.2·10⁻¹³ / 4.2·10⁻¹¹ (midpoint values), not the certified
    2.673·10⁻¹³ / 9.058·10⁻¹¹ that the Letter quotes (UT3 N8).
  * *Stale note in the NF02 file.* `nf02_psi.json`'s `note` field still says "certified margin 3/20 except region (+-, J>=400, 3/8<u<=1)",
    which contradicts "≥ 1/5 everywhere".
  * *Q3\* is not published as a file.* `q3_full_delta0.2.json` stores the weights whose flat action the authors say is **negative** at +−,
    m = M, J ≥ 1990 (supplement.tex:119-122). The defining Q3\* coefficients exist only "by exact polynomial division". Publish them; the
    published data should define the functional the theorem uses.

---

## 4. First report: status of each point in the current text

| first-report item | status | evidence (current text) |
|---|---|---|
| M1 "theorem"/[D] too strong; give net grades; separate layers | **partly** | Net grades given (letter.tex:262-267). Layers separated (letter.tex:84-94; supplement.tex:46-48). Still labelled "Theorem" (letter.tex:223), and "What is established" says "derived" (letter.tex:84-86). New slips: M4(a, c, d, e). |
| M2 hypotheses are labels; Table I mixes inputs and evidence; corollary omits an assumption | **mostly addressed** | Inputs I1–I6 are statements, with a displayed H_K (letter.tex:278-280). Table I is "Evidence … (not an input)" (letter.tex:305). (Pos) is defined (letter.tex:138-142). The corollary lists H_LJ(Δ_g) (letter.tex:381-382). I5(b) "no conspiracy" is still a label (letter.tex:298). |
| M3 k = 4 light side / meaning of g̃₄ | **mostly addressed** | One term-wise functional; Ω₄ = D_LR + 𝒜₄; borderline arc (letter.tex:175-199; supplement.tex:225-264). k = 4 checks on a contact and on crossed and direct exchanges [DV]. Calibration gives 2g₄ (supplement.tex:316-337). Not done: the requested explicit spinning R⁴ ⟨TTTT⟩ contact (not evaluated; supplement.tex:282-283), and Lemma A on Einstein (not evaluated). |
| M4 lower end overstated | **partly** | "no computed constant" is stated (letter.tex:251). Abstract and Discussion still present a two-sided window (letter.tex:36-37, :405-407). New problem: no (Pos) evidence for Ω^Q_low (M2 above). |
| M5 c₁ uncomputed; use the exact quadratic | **addressed (value still open)** | Exact root and supremum forms (letter.tex:228-260; supplement.tex:459-499); "for every value of c₁". Value [I]; c₁ = 0 only given F₂ [A]. The natural-size estimate is stale (m11). |
| M6 de-shift existence, weights, reproducibility | **mostly addressed** | Continuum weights published, C4 ε_W certified, X^Q published, Abel weights re-solved independently (letter.tex:159-167; supplement.tex:145-166, :501-523). Remaining: Q3\* not a file, stale README and NF02 note (m19), C5 ε_W not certified (stated). |
| M7 Lemma C / Einstein instance / collinear corner | **addressed** | Lemma C′ [D] via Lemma W [P], whose identity I checked (supplement.tex:286-314). Einstein instance "derived but not evaluated" (letter.tex:184). Collinear corner L6′ [D] (but see M4(c)). |
| M8 H_K buried / needs J₀ < ½ | **addressed (unreviewed)** | I2 states H_K, with criteria 3/2 and 5/2 and margins ½ and 3/2 (letter.tex:278-287; supplement.tex:364-432). Rests on 22–23 Nov deliverables (M3). |
| M9 1000 relies on the refinement; (1100, 100); uncertified corner | **addressed for the text; the physics remains** | The violation and the "refined form" are in the main text (letter.tex:341-343). Flat certificate now covers everything (letter.tex:310, :325). J = 29 is not in the main text (M1(b)). |
| M10 H_RS preliminary; far-region number | **addressed** | Physical chain = CLS-D's causal domain [P]; O(Δ^N e^{−0.336Δ}) [D, one step D-sketch]; scalar check at Δ = 20 (letter.tex:350-355; supplement.tex:705-723). Far-region number and P_cyc text removed; census claim withdrawn. |
| M11 "lowest gap"; corollary; HLJ_CLOSE reconciliation | **partly** | "Lowest gap" is now qualified (supplement.tex:749-750). Most-fragile state updated to (336, 28) +0.759 (letter.tex:312). Corollary still in the Letter, anchored at 300 with an [I] input (letter.tex:378-383). See M8. |
| M12 title / scope / define g₄ / A1 too strong | **addressed** | New title (letter.tex:17-19). g₄ defined via M(1⁺2⁻3⁻4⁺) (letter.tex:117-123). A1 → A1′ (letter.tex:102-105). |
| M13 ABJM in abstract; 10d flat limit; digits; c_T digression | **addressed** | Not in the abstract. Ten-dimensional flat limit and the "one digit" caveat (letter.tex:392-398). c_T/flavour parenthesis removed from the Letter; only a vestigial clause remains (supplement.tex:529-530). BCJP described (letter.tex:421-422). |
| M14 novelty and credit | **mostly addressed** | TODO removed; KKSDZ, BHM, CHMRS (disp, swamp), BKZ, Kundu, Caron-Huot and BCJP positioned (letter.tex:52-82, :415-422). New: CHLPSpw mis-described, and its holographic bounds not acknowledged (M5). |
| M15 jargon; derivations; chronology | **partly** | Derivation sketches and data added; chronology in the main text (letter.tex:338-341). Jargon and ledger density remain (M6); several derivations are not checkable (M3). |
| minor 1 δ = 0.2 nominal | addressed | "nominal margin δ = 0.2"; flat min 0.1835 (letter.tex:155, :310) |
| minor 2 units; c_T → ∞ | addressed | letter.tex:108, :115 |
| minor 3 parity-odd citation | addressed | CDPrakash, GPY (letter.tex:113-115) |
| minor 4 H17 in the theorem | addressed | letter.tex:258-260 |
| minor 5 "17–30×" | addressed, rounding | "≥ 19×", but 18.6× (m8) |
| minor 6 superseded small-ν at 300 | addressed | not claimed in the Letter; flagged in Tables S5/S6 |
| minor 7 truncated quartic form | partly | warned at letter.tex:375, but the Discussion uses the truncated form (m14) |
| minor 8 300 protocol | addressed | supplement.tex:725-726 |
| minor 9 Regge citation | addressed | KKSDZ Sec. 4.2.4 plus Caron-Huot (letter.tex:130, :273) |
| minor 10 bibliography TODOs | addressed | refs_verified.bib:17; new problem: JHEP years missing (m3) |
| minor 11 front matter | **not addressed** | TODOs (m2) |
| minor 12 "free parameter" | addressed | "universal constant … uncomputed" |
| minor 13 J = 15 value | addressed | +0.6900 own row (supplement.tex:663) |
| minor 14 hyperref warnings | addressed | log: no hyperref/overfull warnings |
| minor 15 grade notation | partly | [NE + M] and [D given X] defined; other composites not (m6) |
| minor 16 twist edge 1.0048 | addressed | 1.00499 / 1.0049925 throughout |
| overall recommendation (reject for PRL) | **unchanged in substance** | The decisive item (positivity on an infinite spectrum) is improved at the binding state but not closed (M1). |

One correction to the first report. Its M13 said ABJM has light R-singlet single traces in T × T (the breathing mode). Under A1′ that is
irrelevant: A1′ allows light spin-0, 2 states. At tree level the current text's [D-sketch] argument against any light exchange is
plausible. The breathing-mode claim itself was never re-derived (UT3 N14).

---

## 5. Spot checks of numbers

| # | number (letter) | traced to | my recomputation | verdict |
|---|---|---|---|---|
| 1 | R₀ = 3.13·10⁻⁶ / 2.89·10⁻⁸; \|ĝ₃\|Λ⁴ ≤ 104 / 120 (letter.tex:235, :408) | UT3 §3; `nf02_psi.json` | 16L_A = 10734.606 / 14339.589; R₀ = 3.132055·10⁻⁶ / 2.891263·10⁻⁸; √(16L_A) = 103.61 / 119.75; k(0) = 1 + 3·10⁻¹¹ | ✓ |
| 2 | 2L_AdS = 72.04 / 92.81; B_eff = −19.220 / −19.265 (letter.tex:246-247) | UT3 §3; `q3_full_delta0.2.json` | 72.04059 / 92.81351; −19.2197 / −19.2647 | ✓ |
| 3 | x/\|c₁\| = 2.16·10⁻³ / 1.67·10⁻⁴; 19941; c₁² coefficient 6.8·10⁻⁵ / 5.4·10⁻⁷; 5.946·10⁴ (letter.tex:232-248) | supplement.tex:459-476 | 1920·I_D = 19941.017; 2.160065·10⁻³ / 1.670256·10⁻⁴; 57600 I_Q²/κ = 541017.0, giving 6.81453·10⁻⁵ / 5.44276·10⁻⁷; 57600 I_Q² = 59462 | ✓ |
| 4 | repaired upper end 125.71 + 1.29·10⁻⁴c₁² / 164.51 + 1.03·10⁻⁶c₁² (letter.tex:250) | supplement.tex:480-484; UT3 §3 | 125.7136 + 1.29379·10⁻⁴ / 164.5115 + 1.03335·10⁻⁶; ε* = 0.430430/11.317027 = 0.038034; margin at ε = 0.04: 0.0223 | ✓ |
| 5 | 1 % thresholds \|c₁\| < 4.6 / 60 and 99 / 1262 (letter.tex:256-257) | C1_VALUE §6; UT3 §3 | 4.607 / 59.57 (item 1); 102.8 / 1305.9 (2a); 98.57 / 1261.8 (2b) | ✓ |
| 6 | ABJM 2.6 % / 3.3 %; ML = 4.2156 λ^{1/4} (letter.tex:393-394; supplement.tex:792-793) | A1PRIME §7 | 2ζ(3)/2L_AdS = 3.337 % / 2.590 %; 2^{5/4}π^{1/2} = 4.21563; λ(Λ = ML) = 2.51·10⁷ / 3.15·10⁹ (text says 2.6·10⁷ / 3.2·10⁹; m12) | ✓ (λ convention) |
| 7 | Λ-forms, Eqs. (4)–(5): 1.3194·10⁴, 8.7415, 6.384, 17.2407, −19.2782, 10.314, 860.3, −2508.6 (letter.tex:366-374) | supplement.tex:569-583 | form vs exact at 300: quartic 72.0396 vs 72.0406 (10⁻³), cubic 3.13215·10⁻⁶ vs 3.13206·10⁻⁶ (3·10⁻⁵); Λ₀ = 8.7415, Λ₀^Q = 4.5889; `repro/bounds.py` (copy) reproduces all its coefficients | ✓ |
| 8 | twist edges 1.0049925Λ (1000), 1.0165851Λ (300) for J = 5 (supplement.tex:639, :644) | m²_gap = (Δ_g−2)(Δ_g−1+2J) | 1.004990 / 1.016585 | ✓ |
| 9 | light-scalar root 0.46Λ (letter.tex:209, :239) | supplement.tex:549-551 | root of −0.430430 + 2.032073 m² is m = 0.46024 | ✓ |
| 10 | 8πG = 256/c_T; t₄ = −240α₃/L⁴ (letter.tex:99, :116) | DICTIONARY_T4 §0 | C_T^{Einstein} = 3L²/(π³G), C_T^{scalar} = 3/(32π²), so c_T = 32L²/(πG) and 8πG = 256/c_T; 2(d−1)d(d+1)(d+2) = 240 at d = 3 | ✓ (parity-odd half of 240\|ĝ₃\| is [D]/[I]: m5/M4) |

**Physics checks requested.**
* **Dictionary.** 8πG = 256/c_T is correct with c_T(free real scalar) = 1. The d = 3 Hofman–Maldacena t₄ normalisation (constant ¼,
  t₄ = ±4 for free scalar / Majorana, collider disc t₄²+α_T² ≤ 16) is the standard one. t₄ = −240 Re ĝ₃/L⁴ is [DV]. The α_T half is
  [D]/[I] (DICTIONARY_T4.md:28), so "R = 240|ĝ₃|" should carry that grade.
* **Flavour currents versus c_T.** The Letter never uses a current two-point coefficient. c_T enters only through 8πG = 256/c_T, and the
  ABJM placement is c_T-free. The project's ABJM c_T code uses (64/π²)∂²_m F, which is valid only because the m_± currents sit in the 𝒩 ≥ 6
  stress-tensor multiplet. ABJM_THEOREM.md:186-189 states this correctly. Spin-1, 3 states are correctly excluded from T × T
  (supplement.tex:527-530).
* **A1′ and ABJM.** A1′ is satisfied by IIA-ABJM:
  * the lightest J ≥ 4 single traces in T × T are first-massive-level string states, at twist ≈ ML − 4;
  * supergravity KK modes have spin ≤ 2;
  * A2 holds as g_s → 0.

  The identification Λ ↔ ML is leading-order (−11 % / −3 % shifts, stated). See m12 on the logical chaining.
* **Logic from positivity to the bounds.** Eq. (1) plus P_O ⪰ 0 and (Pos) gives −Ω[light] ≥ 0. With the stated light side, the quadratic
  inequalities in (t₄, α_T) give Eqs. (2)–(3) exactly. The maximisation over (t₄, α_T) in item (2) is correct, given κ > 0.
* **c₁, F₂, H_EFT(b), κ > 0.**
  * c₁ is handled exactly and correctly.
  * F₂ is needed only for c₁ = 0, which the Letter does not use.
  * H_EFT(b) [A] is needed only for g̃₄ = g₄. It is correctly flagged, but the abstract writes g₄ (m4).
  * κ > 0 and c_L > 0 are [D-sketch + M] and needed only if ĝ₃ ≠ 0. That matters: κR²Λ⁸/57600 can be up to about 0.11 × 1.4·10⁴, far above
    92.8, so a negative κ would destroy item (2) for R ≠ 0.
  * There is no evidence for (Pos) of Ω^Q_low (M2).

---

## 6. What would most increase the chance of PRL acceptance (ranked by impact ÷ effort)

1. **Narrow and rewrite.** *High impact, low–moderate effort.*
   * Make the Letter one clean statement: item (1) at Δ_g = 1000, with item (2a)'s upper end as a secondary result.
   * Move the inputs ledger, c₁ thresholds, corollary and reproduction to the Supplement.
   * Add one figure and one paragraph of mechanism (M6).
   * Fix the grade slips (M4) and the abstract (m4, M7, M8).
   * This alone also creates the length headroom (m1).
2. **Close (Pos) for Ω^EE at Δ_g = 1000.** *Decisive impact, high effort.*
   * Near-gap J ≥ 12 at [NE±].
   * Resolve J = 29.
   * An analytic large-spin / large-impact-parameter argument replacing H_LJ(1000), H_R′ and H_bulk.
   * A certified compact region (M1).
   * Without this, PRL is unlikely. With it, the paper is at or above the CHMRS standard.
3. **Make the derivations auditable.** *High impact, moderate effort.* Full write-ups and scripts for H_K (EE reduction, corner count),
   L6′ (Lemma Y), the R⁴ calibration and the c₁ structure. Get an independent reproduction of the EE reduction and Lemma Y (M3).
4. **Fix positioning and strength statements.** *Moderate impact, low effort.* Correct CHLPSpw; acknowledge the holographic central-charge
   bounds and the CHLPS CFT₃ remark; "≈ 10× in |ĝ₃|, ≈ 10² in |ĝ₃|²"; the IR-log-dominated g₄ window (M5).
5. **The lower end.** *Moderate impact, low effort to drop, moderate to support.* Give Ω^Q_low positivity evidence, or drop the lower end
   from the abstract and Discussion (M2).
6. **Compute c₁ or derive F₂.** *Moderate impact, moderate–high effort.* This removes "up to one uncomputed universal constant" from the
   abstract.
7. **Re-optimise the quartic functional with light-scalar positivity imposed.** *Moderate impact, high effort.* This removes the 0.46Λ
   condition and the ε-repair (M7).
8. **Housekeeping.** *Low impact, low effort.* Front matter, JHEP years, data README, NF02 note, a Q3\* data file, the grade legend,
   notation (m2, m3, m5, m6, m19).
