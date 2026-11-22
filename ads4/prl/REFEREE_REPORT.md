# Referee report

**Manuscript:** "Einstein gravity through quartic order in AdS₄/CFT₃: conditional bounds on the R³ and R⁴ couplings"
(`ads4/prl/letter.tex`, `supplement.tex`, `macros.tex`; switch position `\HLJclosedfalse`, headline Δ_g = 1000)

*Report written 22 Nov 2026 as an internal pre-submission referee. I only read files; nothing was edited or committed. I ran one
cheap Python spot check of the headline arithmetic (≤ 1 s, `ulimit -v 3000000`, one thread). I also ran three web searches, using
only paper titles, to confirm literature I cite below. File references are `path:line` relative to `ads4/`.*

**What I checked.** Every number in the abstract, the theorem and the Λ-forms reproduces from the stored constants. I recomputed:
* the cubic bound, 3.13206·10⁻⁶ / 2.89126·10⁻⁸;
* the ABJM placement, 0.03337 / 0.02590;
* ε_{c₁}/|c₁|, 0.567 / 0.0586;
* the c₁ shift, 0.216 % / 0.0167 %, and the H17 band, 0.238 % / 0.021 %;
* the Λ-form constants 1.3194·10⁴, 8.7415 and 4.5889, and ML = 4.2156 λ^{1/4}.

The arithmetic is not the problem. The problems are how the result is stated, the logic of the k = 4 light side, the scope, and the
positioning.

---

## 1. Summary of claims

The authors construct three fixed light-ray functionals on ⟨TTTT⟩ in a parity-even or parity-violating CFT₃:
* Ω^EE (k = 2, 3, with a new "EE" dual structure);
* Ω^Q_up and Ω^Q_low (de-shifted k = 4, 5).

The setting is two standing assumptions:
* A1: every single trace other than T has twist ≥ Δ_g;
* A2: tree level, c_T ≫ Δ_g².

They claim:

1. **An implication** ("conditional theorem", `prl/letter.tex:164-203`). Assume per-state positivity of the three functionals on every
   twist-allowed heavy block, plus about 20 analytic inputs graded [Q] through [A] (Table I). Then:
   * √(t₄²+α_T²) ≤ 240√(16L_A/(k(0)(1+δ₃)))/Λ⁴ · (1+O(c₁Λ⁻²)), numerically 3.13·10⁻⁶ (Δ_g = 300) and 2.89·10⁻⁸ (Δ_g = 1000);
   * −O(Λ⁻²) ≤ g̃₄Λ⁶/8πG ≤ 17.241 ln(2Λ/3) + B_eff + ε_{c₁}, i.e. 72.04 / 92.81.

   Here c₁ is an uncomputed light-side coefficient, and g̃₄ is a CFT-defined R⁴ coefficient that equals the bulk g₄ given (F)+(P).
2. **Closed-form Λ-dependence** of both bounds, with explicit error terms (`prl/letter.tex:205-228`).
3. **Numerical evidence for positivity** [NE + M] on finite sets of heavy states at Δ_g = 300 and 1000. The margins quoted are:
   * Ω^EE binding value +0.2038 / +0.2009 against the flat 0.200;
   * Ω^Q binding value +0.1427 / +0.1428, controlled by s_d ≈ 0.51 against s_fail = 1.78.

   The decisive below-fold input H_RS is reduced to a cycle-class premise P_cyc.
4. **A conditional corollary** extending positivity, and hence the bounds, to all Δ_g ≥ 300 (`prl/letter.tex:320-335`).
5. **An ABJM consistency check**: the type IIA value 2ζ(3) sits at 2.6 % (1000) or 3.3 % (300) of the upper end.

## 2. PRL criteria

**Importance.** Explicit-constant, two-sided AdS bounds on graviton self-couplings from ⟨TTTT⟩ would be a natural and important
completion of the CLS-D programme, which CLS-D explicitly left open. If the statement were unconditional, or conditional only on a
small number of crisp, physically motivated assumptions, it would be of PRL interest. I do not doubt the importance of the target.

**What is actually delivered.**

1. *The implication is the routine part.* "Positivity ⇒ bound" is the standard dispersive step once the functionals and the light
   side are fixed. The new work is in three places:
   * constructing functionals with the right flat limit (EE, de-shift);
   * the light side (H7, L0, Lemma A, (F)+(P));
   * per-state positivity.
2. *Positivity is the entire physical content, and it is not established.* The heavy spectrum is infinite in J and in mass. The
   computations cover:
   * finite J sets up to 5Λ;
   * at Δ_g = 300, only even J ≤ 32 and odd J ≤ 15.

   Everything else rests on model hypotheses graded [M], [D-sketch + M] or [I]. Even the headline Δ_g = 1000 depends on H_LJ(1000)
   [D-sketch + M], whose unrefined form is *violated* at (1100, 100) (see M9).
3. *The implication itself is not [D].* It contains [D-sketch] steps (H_K, K1′, J₀ < 7/2, (F)+(P), the lower end) and one [A]
   (H_EFT(b)). The de-shift of Ω^Q_up has an open [I] component. By the authors' own rule ("a chain carries the grade of its weakest
   link", `prl/letter.tex:71-72`), item (2) is [A]/[I] and item (1) is [D-sketch].
4. *A1 excludes every known holographic CFT₃*, ABJM included (`RESULT_D_STATEMENT.md:81`). The theorem applies to a class with no
   known member.
5. *The bound is weak compared with what one would want.* It is ≈ 10× weaker than the transplanted flat optimum, and only one of the
   quartic structures is bounded (M12).

**Verdict on the PRL bar.** A "conditional theorem" can meet the PRL bar when the condition is a single clean, plausible statement and
the rest is proved. Two examples are a result conditional on a spectral assumption, and the CEMZ results conditional on A1 and a Regge
bound. Here the condition is an infinite family of block inequalities, supported numerically inside models, plus a ledger of more
than 20 inputs at mixed grades.

That is a research programme reported honestly, not a letter-sized result. It would meet the bar if any of these held:
* (a) positivity were reduced to a finite, certified computation plus a proven asymptotic regime. Concretely: certify the flat limit
  with no float corner, prove large-spin/large-mass positivity analytically, and certify the remaining compact region with interval
  arithmetic;
* (b) the result were narrowed to a statement whose inputs are all at least [D]. Item (1) at Δ_g = 1000 is the obvious candidate,
  once H_K and H_LJ(1000) are upgraded;
* (c) the scope were widened so that it says something about actual theories, e.g. A1 weakened to operators in T×T with a controlled
  light sector, or a supersymmetric version.

The methodological precedent, CHMRS, appeared as a long-form paper. I think that is the right venue for the present material too.

---

## 3. Major issues

Labels: **W** = wrong, **O** = overstated relative to the backing document, **U** = unclear or underspecified. Each item says where
the problem is, why it matters and what would fix it.

### M1. "Conditional theorem" and "the implication is [D]" are stronger than the inputs allow (O)
* **Where:** `prl/letter.tex:29-30` (abstract), `:63-72`, `:164`, `:196`; `UNIFIED_THEOREM.md:172-176`.
* **Problem.** The letter says "the implication is [D] given the inputs". The inputs include:
  * H_K [D-sketch];
  * K1′ [D-sketch];
  * J₀ < 7/2 [D-sketch];
  * (F)+(P) [D-sketch];
  * H_EFT(b) [A];
  * de-shift existence [NE], with C5 spinning coefficients [I] (`prl/letter.tex:248`).

  The lower end of item (2) is itself [D-sketch] (see M4). With the authors' weakest-link rule, the honest net grades are item (1)
  [D-sketch] and item (2) [A], with an [I] in the definition of Ω^Q_up. "[D] given inputs" is formally true of the last step only,
  but readers will read "Theorem … [D]" as "proved modulo A1, A2 and positivity".
* **Why it matters.** The word "theorem" sets expectations that this result does not meet.
* **Fix.**
  * Call it a "conditional result" or "claim".
  * In the statement itself, give the net grade of each item.
  * Separate the three layers: physics assumptions (A1, A2, K1); unproven analytic steps (H_K, (F)+(P), H_EFT, K1′, de-shift);
    and the positivity hypothesis.

### M2. Hypotheses are labels, not statements; Table I mixes inputs with evidence (U)
* **Where:** `prl/letter.tex:231-263`, `:164-166`, `:325-331`.
* **Problem 1: a reader cannot tell what is assumed.** Rows such as these state no mathematical claim:
  * "K2, K3: kernel power counting; Rindler positivity ⇒ swap lemma";
  * "(F)+(P): spinning flat limit, polynomial light side";
  * "H_EFT(b): no conspiracy in heavy sum";
  * "H_K: interchange".

  Positivity is never defined in the letter beyond "PSD per state". The block normalisation N_s, the matrix and the set of
  "twist-allowed heavy blocks" are not given.
* **Problem 2: the table mixes two kinds of row.** The lower half (H_M, H_R, HYB, H_edge, H_LJ, H_cont, H_mono) is not a hypothesis of
  the theorem, which assumes positivity outright. It is the evidence structure for positivity, yet it sits in the same table with the
  same "used by" column.
* **Problem 3: the abstract conflates the two.** It says positivity is "assumed elsewhere through named hypotheses"
  (`prl/letter.tex:34-35`). The theorem actually assumes positivity directly.
* **Problem 4: the corollary omits an assumption.** It says "assume in addition H_mono … and H_cont" (`prl/letter.tex:326-330`), but
  the grade bracket then adds H_LJ(Δ_g), which is not among the assumptions. "In addition to the hypotheses of the theorem at 300"
  (`prl/supplement.tex:303`) already includes positivity at 300.
* **Fix.**
  * Give a numbered list of precise statements (one displayed formula each), at least for H_K, (F)+(P) and H_EFT(a, b), in the letter
    or in the first page of the Supplement.
  * Split Table I into "inputs of the implication" and "evidence for positivity".
  * List every assumption of the corollary in its statement.

### M3. The k = 4 light side and the meaning of g̃₄ (U, possibly W)
* **Where:** `prl/letter.tex:97-99`, `:146-153`; `prl/supplement.tex:70`, `:294-298`; `RESULT_Q_EINSTEIN_CHECK.md:59-66`,
  `:340-344`, `:395-407`.
* **Problem.** EIN finds that the swapped-form pointwise zero holds for every J₀ ≤ k. In particular, a J₀ = 4 contact (the R⁴-type
  control) is annihilated pointwise by Ψ₄. EIN states that "where its (F)+(P) value would sit (an a = b distribution, or the ν
  integration) is not determined" [I].
  * So the functional *as an integral* (D_LR) assigns zero to the R⁴ term.
  * The letter repairs this by *defining* the light side piecewise: D_LR on J₀ < 4 terms, plus the term-wise flat value on J₀ ≥ 4
    terms. It calls this a "clarification" [D-sketch].
* **Why it is not a clarification.** The sum rule Σ_heavy Ω = −Ω[light] holds for *one* functional applied to the full correlator.
  The following are not shown:
  * that the same functional, i.e. the same prescription at the collinear set and in the ν/corner integration, acts on the heavy
    blocks;
  * that this prescription yields 2g₄ on the R⁴ (or heavy-exchange double-trace) terms while still annihilating the J₀ < 4 tree terms
    at the collinear set.

  Otherwise g̃₄ is just "minus the heavy sum" by definition. Its identification with a bulk coupling then rests entirely on (F)+(P)
  [D-sketch] and H_EFT (toy [DV] in the scalar-Mellin case, [A] for (b)).
* **The same issue touches item (1).** EIN itself (`RESULT_Q_EINSTEIN_CHECK.md:396-397`) notes that the Einstein k = 2 pole of
  Result D is also a J₀ = k term whose value "lives at the collinear set or in the ν integration". The letter calls H7 "exact" [D] and
  says item (1) needs no light-side subtleties (`prl/letter.tex:129-130`, `:197`). One sentence is owed on why the H7 value is the
  value of the same prescription used for the heavy states.
* **Fix.**
  * Define Ω^Q as a distribution, including its collinear and ν-boundary terms.
  * Show on one explicit spinning J₀ = 4 position-space contact that it returns 2g₄. EIN's Witten-diagram engine is half-way there.
  * Show that the same definition is used for heavy blocks.
  * Until then, state plainly in the abstract that g̃₄ is *defined* as the k = 4 heavy moment and that g̃₄ = g₄ is a hypothesis.

### M4. The lower end "−O(Λ⁻²)" is graded and presented too strongly (O)
* **Where:** `prl/letter.tex:32`, `:183-186`, `:196`, `:360`.
* **Backing.** `RESULT_Q_STATEMENT.md:178-180` grades it [D-sketch] and says it "really means ≥ −O(Λ⁻²) × (EFT data), from b₀, s_v
  and the contact tower". `UNIFIED_THEOREM.md:263` lists "the lower end g̃₄ ≥ −O(Λ⁻²)" among the [D-sketch] items.
  `prl/supplement.tex:216` confirms there is no computed constant.
* **Problem.** The abstract and Discussion present a two-sided window "with explicit constants", and the theorem grades the
  implication [D]. The lower end has neither an explicit constant nor a [D] grade, and its size depends on unknown EFT data.
* **Fix.**
  * Grade the lower end [D-sketch] in the theorem.
  * Either compute the constant in terms of named EFT data, or state the lower end as "g̃₄ ≥ 0 up to O(Λ⁻²) corrections controlled
    by …".
  * Drop "explicit constants" for the lower end.

### M5. The uncomputed c₁ makes the headline numbers incomplete (U)
* **Where:** `prl/letter.tex:166-168`, `:172-180`, `:191-192`; `RESULT_D_STATEMENT.md:44-47`.
* **Problem 1: c₁ is not a parameter of the theory.** It is a definite number fixed by the functional: the Regge residue of the E×R³ T
  exchange. Calling it "free" hides that it is simply uncomputed.
* **Problem 2: c₁ is unbounded.** "There is no argument for |c₁| ≤ 1", so the quoted ∓0.22 % per unit is not a size estimate. The
  letter's linearised form "(1 + O(c₁Λ⁻²))" is valid only for |c₁| ≪ Λ²√(16L_A)/19941 ≈ 460 (Δ_g = 300).
* **What does survive.** The light term is linear in t₄ while the |ĝ₃|² term is quadratic, so the exact inequality is a quadratic with
  a finite solution for every c₁. The bound is therefore not void, but it grows like |c₁| for large c₁.
* **Fix.** Either:
  * (a) compute c₁. REVIEW_EE m5 argues c₁ = 0 is plausible at leading Regge, and proving that would remove an entire error column;
    or
  * (b) state the exact quadratic-root bound, valid for all c₁, instead of the linearisation.

### M6. De-shift existence, the definition of Ω^Q_up, and reproducibility of the weights (O/U)
* **Where:** `prl/letter.tex:132-137`, `:248`; `RESULT_Q_STATEMENT.md:249`, `:362-363`; `REPRO.md:150-152`, `:164`.
* **What is and is not established.**
  * The k = 2 de-shift is [D], with a 2-dimensional null space fixed by a gauge choice.
  * The k = 4, 5 continuum solutions are [NE] only "to t = 30".
  * The C5 spinning j′ = 1 coefficients are [I].
  * The stored weights have high moments that grow from a ν-grid cutoff.
  * Two documents disagree on those moments by factors up to 5.
  * The Abel weights were never reproduced independently.
* **Why it matters.** A "fixed functional" whose definition contains an [I] component, and whose numerical realisation exists in a
  single implementation, is not yet a fixed functional. The letter mentions [NE] but does not flag that Ω^Q_up is incompletely
  defined.
* **Fix.**
  * Close the C5 spinning coefficients, or show the bound is insensitive to them.
  * Publish the weights (ψ₂, ψ₃ for NF02; Q3-full-δ0.2; the Abel solutions) as data in the Supplement.
  * Have the Abel inversion reproduced independently.

### M7. Lemma A is graded [D], but its input Lemma C is now known to be mis-stated (O); the Einstein-instance grading is incoherent (U)
* **Where:** `prl/letter.tex:144-147`, `:244-245`, `:309-317`; `RESULT_Q_EINSTEIN_CHECK.md:193`, `:390`, `:422-423`.
* **Problem 1: Lemma C.** EIN finds that Lemma C's "one analytic function h(t₁, t₃)" is false for interacting correlators: they are
  not single-valued on (ℂ∖cut)². Only the weaker statement "analytic between the hairpin and C" holds, and it is verified [DV] on
  examples, not derived. Lemma A's [D] rests on the old Lemma C.
* **Problem 2: the Einstein instance.** It is listed as [I] while Lemma A is [D].
  * If the general derivation covers Einstein, the instance is [D] and "[I]" just means "not numerically checked". Say so.
  * If it does not, Lemma A is not [D] in the generality item (2) uses, and item (2) is [I].
* **Problem 3: the collinear corner.** EIN did not compute whether the Δφ⁻³ collinear growth is integrable against the celestial
  kernel (`RESULT_Q_EINSTEIN_CHECK.md:339-340`). "The collinear set carries nothing" (`prl/letter.tex:148-149`) is therefore [DV]
  only for the cancelling total at Δφ ≥ 0.025.
* **Fix.** Grade Lemma A as [D given Lemma C′], derive C′ from the s-channel OPE (EIN §6 item 3), and reword the Einstein-instance
  entry.

### M8. H_K, the main analytic input, is buried (U)
* **Where:** `prl/letter.tex:114-115`, `:242`; `prl/supplement.tex:63`; `UNIFIED_THEOREM.md:299`.
* **Problem.** The Supplement and UT say H_K is [D given H_ip], where H_ip means per-ν heavy-sum convergence with J₀ < ½, which is
  "not implied by KKSDZ". UT calls H_K "the main analytic hypothesis; it has no criterion short of H_ip".
  * A Regge intercept below ½ is not expected in general: non-perturbatively only J₀ ≤ 1 is known.
  * So H_K must be proved by the other route ("Proposition K"), which does not exist.

  The letter mentions only "sketched [D-sketch]".
* **Fix.** State H_K precisely in the letter, say that the only sufficient criterion currently known (J₀ < ½) is stronger than the
  known Regge bound, and say what Proposition K would need to show.

### M9. Positivity at the headline Δ_g = 1000 rests on the same kind of step that makes 300 "[I]" (U/O)
* **Where:** `prl/letter.tex:257`, `:266-274`, `macros.tex:48-54`; `RESULT_D_STATEMENT.md:249-255`; `HLJ_CLOSE.md:279-281`;
  `REPRO.md:119-121`.
* **The 1000 headline relies on a refinement.** The status sentence explains that H_LJ(300) is open because the sign-free tolerance
  is met only after the D_h class is removed, a [D-sketch] step. But the sign-free H_win is *violated at (1100, 100)* at Δ_g = 1000
  (C2 +− = −2154 ± 43 against T_B = 1409). So H_LJ(1000) also rests on exactly that [D-sketch] refinement, with margin ×4.8 instead of
  ×1.1.
* **Parts of the 1000 window have no data at all.** ν < 10 is unmeasurable at every J, and for J ≥ 200 the whole window lies below
  ν ≈ 10. HLJ_CLOSE §5 adds that the engine cannot decide H_win at (1100, 100).
* **None of this is in the letter's main text.** The supplement notes the violation (`prl/supplement.tex:100`).
* **Two related points.**
  * The flat certificate is called "[certified]" with "one corner float-checked only" (`prl/letter.tex:123`). A certificate with an
    uncertified corner is not a certificate.
  * REPRO finds that D's binding values equal the flat NF02 values to 10⁻⁴. The genuinely AdS content of the +0.2038 is a −0.0005
    correction that was not reproduced. That is fine, but it shows the binding state is not where the risk lies; the risk is the
    large-spin region discussed above.
* **Fix.**
  * Say explicitly in the letter that the 1000 headline also needs the D_h-shape refinement, and quote the (1100, 100) violation.
  * Grade the flat certificate "[certified] except the corner (+−, J ≥ 400, 1.5Λ ≤ m < 1.633Λ) [NE]", or certify that corner.

### M10. H_RS: the upgrade rests on two unreviewed same-day analyses, and the far-region number is the wrong one (O)
* **Where:** `prl/letter.tex:282-296`, `:253`; `RD_CLOSE.md:231-238`; `REPRO_ENGINE.md:35-42`; `HRS_CLOSE.md:311-313`.
* **The upgrade.** In the consolidated statement, H_RS for Q's binding state is [I], ranked weakest link #1. It moves to
  [D given P_cyc] on the strength of HRS_CLOSE and RD_CLOSE. Both are unreviewed and dated the same day as the draft. The Supplement
  summarises them but contains neither the proofs (Lemmas 1–2, Proposition RS) nor the divisor analysis, so an outside referee cannot
  check them.
* **Far-region weight.** The letter quotes "≤ e^{−37} (300) against e^{−8.5}". That is a *pointwise* maximum at r ≥ 0.2. RD_CLOSE
  itself says the thimble volume factor can be up to e^{40} at Δ = 300 and that the conclusion holds "for r ≥ 0.3"; at r ≥ 0.2,
  e^{−37+40} > e^{−8.5}. The supporting number is e^{−66} at r ≥ 0.3, i.e. ≲ e^{−26} after the volume factor.
* **P_cyc is not tested by either reproduction.** Both integrate over the real section. The letter's main text omits this; it appears
  only in the Supplement (`prl/supplement.tex:360-362`).
* **Coverage.** The closure covers only ++ J = 0, 2 and OO J = 0. The other measured Q channels (+− J = 4–11, ++ J = 6, 8) use the
  same engine without the bookkeeping. The margins there are large, but this should be stated.
* **Fix.**
  * Put the contiguity argument and the divisor table into the Supplement in checkable form.
  * Quote the r ≥ 0.3 far-region number, with the volume factor.
  * State in the main text that P_cyc is untested by the reproductions.
  * Until HRS_CLOSE and RD_CLOSE are reviewed, present H_RS as "[D given P_cyc] (preliminary)", with the previous [I] grade mentioned.

### M11. Internal inconsistencies about Δ_g = 300 and the range corollary (W/O)
* **"Lowest gap" vs "open".** The letter says "Three hundred is the lowest gap at which the present method closes"
  (`prl/letter.tex:333-334`; `prl/supplement.tex:326-327`). In the same letter, the H_LJ(300) window is [I] and that is the reason
  for headlining 1000 (`macros.tex:48-54`). DG_SCAN says the sign-free window closes only from Δ_g ≈ 440 with the measured 450 c₂;
  ≈ 275 needs the [D-sketch] refinement (`DG_SCAN.md:39`, `:195`).
* **The corollary revives a withdrawn claim.** "Δ_g ≥ 300" is "not claimed" (`UNIFIED_THEOREM.md:310`) and "withdrawn"
  (`RESULT_D_STATEMENT.md:480`). The letter revives it as a corollary, from an unreviewed proxy analysis (DG_SCAN). In fallback mode
  the corollary is anchored at a Δ_g whose positivity input is open. Its extra assumption H_LJ(Δ_g) for *every* Δ_g ≥ 300 is close to
  assuming the conclusion.
* **HLJ_CLOSE (unfinished; PLACEHOLDER-SUMMARY etc., `HLJ_CLOSE.md:19`, `:168`, `:203`, `:278`) changes the 300 picture.**
  * The (336, 28) gap state comes out at **+0.759 / +0.765**, below the letter's "most fragile computed state at 300 is J = 32
    (+0.783)" (`prl/letter.tex:272-273`).
  * The binding margin of that window is a *structure* error of only ≈ +215 c₂, i.e. 3–6×. The letter's quoted failure thresholds are
    of order −5·10³ c₂ (`HLJ_CLOSE.md:213`, `:221`).
  * The window grade moves from [M] to [NE + M].

  The switch rightly stays "open" (HLJ_CLOSE §6: J ≥ 34 stays [I]). But the "(336, 28/32) window [M]" entries and the "most fragile"
  sentence will need editing once HLJ_CLOSE is final. `macros.tex:21-25` anticipates this only for the closed position.
* **Fix.**
  * Drop the "lowest gap" sentence, or qualify it as "with the [D-sketch] refinement".
  * Move the corollary to the Supplement as a remark, or anchor it at Δ_g = 1000 (NOTES §3 "option not used").
  * Reconcile with HLJ_CLOSE before submission.

### M12. The title and scope claim more than is bounded (O)
* **Where:** title (`prl/letter.tex:14-15`), Discussion `:356-361`.
* **Only one quartic coefficient is bounded.** "Einstein gravity through quartic order" suggests every coupling through eight
  derivatives is bounded. In D = 4 there are several independent eight-derivative four-graviton contact structures. They enter
  helicity-preserving and helicity-violating amplitudes differently (see e.g. Chowdhury et al., arXiv:1910.14392, for the
  classification). The letter bounds one coefficient, g₄, normalised through the type II Virasoro–Shapiro amplitude and presumably the
  MHV function f(s,u). The authors themselves note that forward and helicity-flip sum rules are absent in AdS
  (`prl/letter.tex:365`).
* **No known theory is covered.** A1 excludes every known holographic CFT₃. A1 is also stronger than needed: only single traces
  appearing in T×T (bosonic R-singlets) matter for ⟨TTTT⟩.
* **Fix.**
  * Retitle, e.g. "Conditional bounds on R³ and one R⁴ coupling from ⟨TTTT⟩ light-ray sum rules in AdS₄/CFT₃".
  * Define g₄ by an explicit amplitude, e.g. M(1⁺2⁻3⁻4⁺) = ⟨23⟩⁴[14]⁴ f(s,u), f = 8πG/(stu) + … + g₄ + …, and list which quartic
    structures are left unbounded.
  * Weaken A1 to operators in T×T, and discuss whether any theory satisfies it.

### M13. ABJM: outside the hypotheses, and in the abstract (O/U)
* **Where:** abstract `prl/letter.tex:35`; `:338-352`; `prl/supplement.tex:374-387`.
* **Problem 1: misleading prominence.** ABJM violates A1, and even restricted to T×T it has light R-singlet single traces (e.g. the
  breathing mode). A placement number for it in the abstract reads as support for the bound. The letter itself says it is "not
  evidence".
* **Problem 2: the flat limit is 10d.** The flat limit of AdS₄×CP³ is ten-dimensional, because the internal radius equals the AdS
  radius. The 4d coefficient "g₄" is the 10d R⁴ term restricted to 4d momenta, and the 4d sum rule's light side would also receive
  the light KK exchanges.
* **Problem 3: the Λ ↔ ML dictionary.** It identifies the J = 0 gap mass with the string mass. The lowest string states have spin up
  to 4 and twist ≈ ML − 4, so the identification holds only at leading order. That is stated, but the 3.3 % / 2.6 % are quoted to two
  digits.
* **Problem 4: an irrelevant digression.** The c_T / flavour-current parenthesis (`prl/letter.tex:347-350`) does not enter the
  placement. It reads as a defence against an objection nobody raised. ABJM_THEOREM.md is likewise irrelevant to the letter except
  for the normalisation 8πG = 256/c_T.
* **A better comparison exists.** The R⁴ coefficient of ABJM(ABJ) has been obtained within the CFT by localisation plus the
  stress-tensor-multiplet bootstrap: Chester–Pufu–Yin (arXiv:1804.00949) in the M-theory regime, and Binder–Chester–Jerdee–Pufu
  (2011.05728) in the IIA regime.
* **Fix.**
  * Remove ABJM from the abstract.
  * Keep a two-line footnote-level sanity check of sign and normalisation of g₄, citing those CFT determinations.
  * Delete the c_T / flavour parenthesis.

### M14. Novelty and credit (U, and partly W)
* **Where:** `prl/letter.tex:41-55`, `:366-370`.
* **Problem 1: a TODO in the text.** "remains to be discussed in detail. % TODO" cannot go to a referee.
* **Problem 2: an unconvincing novelty statement.** "Literature searches (abstract level, plus one extended search) … novelty is
  asserted only at search level" is honest, but it is not a statement a PRL referee can accept. The authors must position the work
  themselves.
* **Problem 3: "no published bound on R³ … from ⟨TTTT⟩ sum rules" must be reconciled with:**
  * KKSDZ (1904.05905). It is cited only for the Regge bound, but its superconvergence (shock-commutativity) sum rules are
    light-ray sum rules on ⟨TTTT⟩-type correlators, used there to constrain stress-tensor couplings in large-gap theories.
  * Belin–Hofman–Mathys, "Einstein gravity from ANEC correlators" (arXiv:1904.05892, JHEP 08 (2019) 032). They show from ANEC
    commutativity, in large-N, large-gap CFTs, that the dual is Einstein gravity (a = c in 4d). This is a ⟨TTTT⟩-based argument for
    Einstein ⟨TTT⟩.

  The difference (explicit constants, AdS₄, R⁴) is real, but it has to be stated against these.
* **Problem 4: Kundu is not described.** Kundu (2104.11238) is about scalar couplings φ²□ᵏφ² in AdS, from the CFT Regge limit:
  positivity, monotonicity and log-convexity, with gravity modifying only the k = 2 case. Say that in one sentence rather than "to be
  discussed".
* **Problem 5: missing citations.**
  * Caron-Huot–Mazáč–Rastelli–Simmons-Duffin, "Dispersive CFT sum rules" (2008.04931): the functional framework.
  * CHMRS, "Sharp boundaries for the swampland" (2102.08951): gravity bounds with an IR log.
  * Bern–Kosmopoulos–Zhiboedov (2103.12728): 4d four-graviton R⁴ coefficients from dispersive bounds.
  * Caron-Huot's Lorentzian inversion paper (1703.00278), for the Regge bound alongside KKSDZ.
  * The N = 6 bootstrap (BCJP) should be described, not only cited.
* **Fix.** Write a proper positioning paragraph covering these. Remove the search-level language.

### M15. Presentation and reproducibility below PRL standard (U)
* **Where:** throughout; `prl/supplement.tex` as a whole.
* **Problem 1: undefined internal vocabulary.** The letter relies on terms a PRL reader cannot decode: "engine", "crystal-cycle",
  "fold", "near-gap chains", "hybrid model/HYB", "K₂", "s_d", "s_fail", "P_cyc", "twisted local systems", "swapped form", "D_LR",
  "Regge window", "D_h class", "c₂", "N_s". Roughly half the main text is a hypothesis ledger.
* **Problem 2: the Supplement is a status report, not a derivation.** It gives grades and numbers but none of the following:
  * the EE structure's admissibility;
  * the H7 pole;
  * L0;
  * Lemma A;
  * the de-shift equation and its solution;
  * the functional coefficients;
  * the certificate;
  * the positivity matrix.

  A referee cannot verify the core claims.
* **Problem 3: in-sample states and gauge fragility are not disclosed.** The functional's structure (q² → EE) and gauge
  (canonical → min-norm) were chosen *after* heavy-state failures were seen. J = 18, 20, 24, 28 at 300 are in-sample, and in the
  canonical gauge the J = 24 state fails badly, with n = 1.55 or −22 … −40 (`RESULT_D_STATEMENT.md:129-135`). Any fixed functional
  gives a valid implication, but this weakens the *evidence* for positivity, and the letter does not say so. "In-sample" appears only
  in a Supplement caption.
* **Fix.**
  * Define every term used, or remove it from the main text.
  * Add derivations and the functional data (machine-readable) to the Supplement.
  * Add one sentence on the chronology and the gauge sensitivity.

---

## 4. Minor issues

1. **"Margin δ = 0.2" is nominal.** REPRO D9 finds the flat minimum of Q3-full-δ0.2 is 0.1966 on a dense grid (`REPRO.md:104-109`).
   Say "δ ≈ 0.2 (grid LP)" (`prl/letter.tex:140`).
2. **Units of t₄.** "√(t₄²+α_T²) ≤ … in AdS units" (`prl/letter.tex:178`). t₄ and α_T are dimensionless CFT data; the units refer to
   ĝ₃. Also state that the bound concerns the c_T → ∞ limit at fixed Δ_g.
3. **Parity-odd ⟨TTT⟩ structure.** Cite a source for the α_T notation; HofmanMaldacena/Buchel do not cover it in d = 3.
4. **H17.** It allows separate t₄², α_T² coefficients at O(Λ⁻²), in which case the bound is on a weighted norm
   (`RESULT_D_STATEMENT.md:190`). Mention this in the theorem.
5. **"Every unmeasured channel needs |s_d| ≥ 17–30×"** (`prl/letter.tex:282`). RQ-S §3.1 gives s_fail ≥ 11 for unmeasured channels,
   i.e. ≥ 19–22×. Reconcile the lower end.
6. **Superseded small-ν treatment at 300.** J = 0, 2 at Δ_g = 300 still use it (margin ≈ 1). Note this in the letter's coverage
   sentence (`prl/letter.tex:268-269`).
7. **Truncated quartic Λ-form.** With only the 1/Λ term it gives 72.017 at 300 against the exact 72.041. That is on the unsafe side,
   as DG_SCAN itself notes (`DG_SCAN.md:289`). The letter's "≤ 10⁻³ with all terms shown" is correct; keep it that way and do not
   quote the shorter form as a bound.
8. **HLJ_CLOSE's 300 binding value.** It reports +0.2054 for the Δ_g = 300 binding state under its chain protocol
   (`HLJ_CLOSE.md:178`), against the adopted +0.2038. Harmless, but record which protocol defines the quoted minimum.
9. **Regge-bound citation.** "KKSDZ … with its stated loophole": give the section, and add Caron-Huot 1703.00278.
10. **Bibliography.** It has 25 `TODO(bib)` markers. The BCPW author list is inconsistent between the repo's documents (NOTES §5).
    The Fmax entry has only an arXiv number.
11. **Front matter.** Affiliations, emails and acknowledgments are TODO (`prl/letter.tex:18-22`, `:381`;
    `prl/supplement.tex:20`, `:22`). `macros.tex:33-34` would print "[TODO(HLJ)]" if the switch were flipped.
12. **"c₁ … free parameter".** Use "uncomputed coefficient" (see M5) in the theorem and in Table I (`prl/letter.tex:251`).
13. **Record caveats.** REPRO D4/D5 (J = 15 at 1000; J = 16 at 300) appear in the Supplement only. They are fine, but quote the J = 15
    value as "≈ +0.69 (model)" rather than +0.869 in Table S5.
14. **Hyperref warnings.** Unicode tokens in PDF strings (from the section titles with math). Cosmetic.
15. **Grade notation.** "[NE + D-sketch]" and "[D given P_cyc]" are not in the grade legend (`prl/letter.tex:69-71`); add them.
16. **Twist edge.** "+0.2009 at 1.00499Λ" is now right (record fix). Make sure RESULT_D_STATEMENT's "1.0048Λ" does not reappear
    anywhere in the Supplement.

---

## 5. Recommendation

**Reject for PRL in its present form, with encouragement to resubmit a long-form version to JHEP (or PRD).** I do not consider this a
"major revision" case for PRL.

* **The calculations are sound.** The light-side constants are correct and independently reproduced, and the bookkeeping of what is
  and is not established is unusually honest. That honesty is exactly what shows the problem: the deliverable is an implication whose
  decisive hypothesis (positivity on an infinite heavy spectrum) is supported only numerically inside models.
* **The implication itself is not [D].** It contains [D-sketch] and [A] steps (M1, M8).
* **Its central new quantity g̃₄ is defined by fiat** at the point where the functional's integral gives zero (M3).
* **Its scope covers no known theory** (M12).
* **It cannot be checked from the Supplement** (M15).
* **The novelty claim is not yet made against the closest literature** (M14).

A JHEP paper with full derivations, the functional data, and a statement built around item (1) at Δ_g = 1000 would be a valuable
contribution. A future letter would become appropriate if M3 is settled by explicit computation, c₁ is computed, and positivity is
certified at one Δ_g.

**The three fixes with the highest leverage:**
1. **Restate the result honestly and precisely** (M1, M2, M4, M7, M12).
   * Call it a conditional result, not a theorem.
   * Give precise displayed statements of H_K, (F)+(P) and H_EFT.
   * Give each item its net grade: item (1) [D-sketch]; item (2) [A]/[I]; lower end [D-sketch] with no explicit constant.
   * Separate inputs from positivity evidence.
   * Scope the title to "R³ and one R⁴ coupling".
2. **Make the k = 4 light side a theorem about one functional, and compute c₁** (M3, M5).
   * Define Ω^Q as a distribution, including collinear and ν-boundary terms.
   * Show on an explicit spinning J₀ = 4 contact that it gives 2g₄, and that the same definition acts on the heavy blocks.
   * Compute the E×R³ Regge residue c₁, or replace the linearisation by the exact quadratic bound.
3. **Fix the positioning and the ABJM framing** (M13, M14).
   * Add a real literature paragraph: KKSDZ superconvergence, Belin–Hofman–Mathys, CHMRS 2008.04931 and 2102.08951,
     Bern–Kosmopoulos–Zhiboedov, an accurate description of Kundu, and BCJP / Chester–Pufu–Yin for ABJM.
   * Remove the TODO and the search-level language.
   * Take ABJM out of the abstract.

---

### Sources consulted for the literature checks
* [Belin, Hofman, Mathys, "Einstein gravity from ANEC correlators", arXiv:1904.05892](https://arxiv.org/pdf/1904.05892)
* [Kundu, "Swampland Conditions for Higher Derivative Couplings from CFT", arXiv:2104.11238](https://arxiv.org/abs/2104.11238)
* [Bern, Kosmopoulos, Zhiboedov, "Gravitational effective field theory islands, low-spin dominance, and the four-graviton amplitude", arXiv:2103.12728](https://arxiv.org/abs/2103.12728v2)
* Cited from memory, not re-checked online: Chowdhury et al. 1910.14392, Chester–Pufu–Yin 1804.00949, CHMRS 2008.04931 and
  2102.08951, Caron-Huot 1703.00278. The authors should verify the precise content attributed to these.
