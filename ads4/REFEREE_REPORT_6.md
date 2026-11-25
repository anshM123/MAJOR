# Referee report (round 6): "How Close to Einstein Must Holographic Gravity Be?"

Letter: `prl/letter.tex` (built `prl/build/letter.pdf`, 6 pp.). Supplemental Material: `prl/supplement.tex` (built `prl/build/supplement.pdf`, 36 pp.).
Line numbers refer to the `.tex` sources as they stand on 2026-11-25.

## Verdict: major revision

The submission makes the program of Chang, Landau and Simmons-Duffin concrete for AdS4/CFT3. It derives (i) an explicit-constant implication "positivity of a fixed functional ⇒ τ ≤ 2.89×10⁻⁸ at Δ_g = 1000", (ii) an upper bound on the R⁴ coefficient, and (iii) an assumption-light "angular" bound τ ≤ 0.0178. The machinery (de-shifting, the EE structure, a certified leading-Regge floor, the energy-detector ladder) is new and impressive, and the grading of claims is far more careful than is usual. Still, I cannot recommend publication in the present form, for four reasons:

- **The abstract reads as more than the paper shows.** The headline 10⁻⁸ statement is conditional on assumption (iii), which in AdS is established only at leading Regge order and otherwise rests on numerics, model premises and a decay assumption.
- **The three headline statements are worded inconsistently** across the abstract, the body, the Result and the Discussion:
  - The certified floor is "one fifth of the flat margin", "at least one fifth of its normalisation" and "0.1999 / 0.199" in different places.
  - Its scope is "every heavy state" in the Letter but "target weights only, β ≲ 200 for the published weights, OO block conditional" in the Supplement.
  - The angular bound has one list of inputs in the Letter and a different one in the Supplement, and it is absent from the Result and the Discussion.
- **The Supplement contradicts itself.** Several of its grades predate the certified floor and still call the same object [NE]. Its Sec. S1 "legend", which the Letter says grades every statement, maps about ten Letter phrases that no longer exist and omits the newest claims.
- **Key proofs point to unpublished internal files** (`ANGULAR_*.md`, `code/t0cert/`, `code/angular/`), and the Supplement is written in project-internal language that an outside reader cannot follow.

All of this can be fixed, and the physics may well warrant a Letter. But the paper must say plainly what is proved, what is assumed and what is modelled, in the same words everywhere.

---

## 1. Fatal or major scientific problems

### M1. The abstract presents a conditional result as if positivity were established
**Where:** `letter.tex` lines 21–26.

> "with a gap of 1000, the graviton's three-point couplings match Einstein gravity to one part in 10⁸ … The input is causality plus a positivity property we prove in flat space and, with a certified margin, at leading Regge order in AdS."

**Problem.** A reader will conclude that the positivity input is proved. It is not. In AdS it is *assumed* (iii). The Supplement's Sec. S13 says plainly that "an unconditional bound" lies outside the results. The evidence beyond leading Regge order is error-controlled numerics on a finite set of states ([NE±], an estimate, not a bound), model values ([NE+M], [M]), and two large-spin "theorems" whose premises include an assumption ((A_GU), [A]) and model inputs ([M]).

The abstract also omits other conditions:
- the result holds at leading order in 1/c_T;
- the quartic bound needs κ > 0 or τ = 0, plus the scalar condition;
- "two-percent bound at every spin" is ambiguous. Two percent of what? τ ≤ 0.0178 is about 0.45% of the collider bound τ ≤ 4. And the bound is not "at every spin": the *contributions* are non-negative at every spin.

**Fix.** Make the conditionality explicit. Suggested abstract (about 550 characters rendered):

> Holographic CFTs with a large higher-spin gap are expected to be dual to Einstein gravity, but only parametrically. In three dimensions we supply constants. If a positivity property holds that we prove in flat space, certify at leading Regge order in AdS and support numerically beyond it, then a gap of 1000 forces the stress-tensor three-point couplings to their Einstein values to 3×10⁻⁸, and bounds the first stringy correction, with ABJM theory inside. Without that property, commutativity of energy detectors and crossing still give τ ≤ 0.018.

Also replace "large-gap quantum field theories" (line 21). A large gap alone is not enough; large c_T is also required.

### M2. The certified leading-Regge floor ("T0") is stated in four mutually inconsistent ways
**Where:** `letter.tex` lines 64–65, 91–94, 253–258, 276–277; `supplement.tex` lines 1207–1222 (Sec. S7.B).

1. Intro, l. 64–65: "its leading Regge part keeps **one fifth of the flat margin** on every heavy state". This is wrong. The flat floor is itself 1/5 of the normalisation N_s, and the AdS floor is 0.1999 N_s. Curvature therefore keeps about 99.95% of the flat margin, not 20% of it. As written, the sentence understates the result by a factor of five and contradicts l. 257–258 ("curvature costs essentially nothing").
2. Setting (iii), l. 92–93: "at least **one fifth** of its natural normalisation on every heavy state, at every gap Δ_g ≥ 300". This is false as stated: the certified constant is 0.199 < 1/5 for Δ_g ≥ 300, and 0.1999 < 1/5 for Δ_g ≥ 1000. "Natural normalisation" is also undefined in the Letter.
3. Positivity paragraph, l. 254–256: "at least 0.1999 … on every heavy state of every spin, mass and impact parameter, for every gap Δ_g ≥ 1000". The Supplement's own *Scope* (l. 1219–1221) is narrower:
   - "the target (de-shifted) weights; the realised published weights are covered at Δ_g = 300 and 1000 for β ≲ 200 only";
   - "the OO block additionally needs the EO cancellation and κ > 0".
   So "every heavy state" excludes the parity-odd even-J (OO) block unless further inputs are granted. "EO cancellation" is undefined anywhere, and κ here is apparently not the quartic κ (see m9). Nor is "every impact parameter" certified for the functional actually published.
4. Larger gaps, l. 276–277: "the leading-Regge ones [hold] wherever the de-shifted weight exists (we construct it up to 5000)". This conflicts with "for every gap Δ_g ≥ 1000" at l. 256. The existence of the minimal-norm weight is itself graded [D + NE] (Table S2, row DS).

**Fix.** Use one sentence everywhere, for example:

> "We certify that the leading-Regge part of Ω^EE is at least 0.1999 N_s (0.199 N_s for Δ_g ≥ 300) on every heavy state in the parity-even and odd-J sectors, for the exact de-shifted weights at every gap where they exist; the flat value is N_s/5, so curvature costs less than 0.1% of the margin."

Intro l. 64: replace "keeps one fifth of the flat margin" with "keeps all but 0.1% of the flat margin". Setting l. 92: replace "at least one fifth of its natural normalisation" with "at least 0.199 of its flat-space normalisation N_s (the flat floor is 0.2)". Also state the OO-block and published-weight caveats once, in the Letter or in a footnote.

### M3. The Supplement still grades the same leading-Regge floor as [NE] in several places
**Where:** `supplement.tex`:
- l. 1157–1158 (S7 intro): "Theorem W also [rests] … on a leading-Regge floor that is [NE] at its lower spins".
- l. 1251–1252 (S7.C): "for m < 1.3Λ at b ≥ 114.06, [NE] below".
- l. 1351 (S7.F, premise P0): "[NE] below; hybrid vs F_LR [NE]".
- l. 1404 (Table S9 reading): "the floor at m < 1.3Λ below b = 245 is [NE]".
- Table S11, l. 1480: "hybrid floor rigorous on every b ≤ 2Λ (Lemma-C′ part at sampled masses)".
- Table S17: "hybrid floor … [D + certified] (Lemma-C′ part [P + DV] at sampled masses)".

**Problem.** S7.B claims the hybrid (leading-Regge) action is certified ≥ 0.1999 N_s for *every* m ≥ Λ and every b̂ > 0. Either T0 supersedes these rows, or it does not cover them (for example, if the difference is target versus published weights, or hybrid versus F_LR). A referee cannot tell which. If T0 supersedes them, Theorem W's premise P0 becomes certified and the theorem's grade should improve. If it does not, the Letter's "every heavy state" (M2) is wrong.

**Fix.** Rewrite S7 intro, S7.C, S7.F (P0), the Table S9 reading, Table S11 and Table S17 so that they cite T0, and state exactly which object (target or published weight; hybrid or F_LR) each grade refers to.

### M4. The angular bound: its inputs are stated inconsistently, and "no positivity assumption" needs qualifying
**Where:** `letter.tex` lines 301–310, 25–26 (abstract), 57–62 (Intro), 209–248 (Result), 324–339 (Discussion); `supplement.tex` Sec. S11, l. ~1975–2066.

1. **Different input lists.** The Letter (l. 308–310) names three inputs: "planar regularity of multi-trace states, no scalar with Δ ≤ 1 in T×T, and light-cone convergence of detector commutators". The Supplement's *Inputs* (l. 2019–2021) names four:
   - "(1) the first-moment and ladder identities **as physical operator identities**";
   - "(2) the contact coefficient";
   - "(3) fixed-window planar regularity";
   - "(4) light-cone convergence **and existence of the mixed commutator**".

   The Letter presents (1) as established fact ("fix this ladder exactly", l. 304–305), and it never mentions the existence of the mixed commutator or the J₀ < 2 requirement (S11).
2. **Hidden use of (i) and (ii).** The bound needs the twist gap for J ≥ 6 ("twist ≥ g, J ≥ 6"), which is assumption (i). It also needs L_g → 0 from planar factorisation, which is (ii) and holds only at leading order in 1/c_T. Neither the Letter paragraph nor the abstract says this. Since the Letter states that "our results rest on three assumptions" (l. 81), a reader cannot tell which of (i)–(iii) the angular route uses.
3. **"With no positivity assumption" is accurate only narrowly.** The route uses ANEC (the "fourth-harmonic ANEC rule"), reflection positivity and unitarity. These are theorems or axioms, not assumptions, but the phrase invites the objection. Write "without assumption (iii)".
4. **"At every spin".** The all-J statement rests on "the all-m helicity norms (checked for J ≤ 30)" (S11.A, end) and a pipeline check "for J ≤ 40" (S11). Is the all-J helicity-norm formula proved, or quoted? If it is checked only to J ≤ 30–40, then "at every spin" in the abstract and l. 305 is not shown.
5. **Grade and provenance.** The Letter gives the angular result no grade word, and Sec. S1, which "grades every statement" (Letter l. 72–73), does not map it at all. In S11 the chain is "confirmed" by "an independent hostile audit" (l. 2010), and the "derivations, audits and check scripts" are `ANGULAR_ROUTE.md`, `LADDER_ALLJ.md`, `ANGULAR_MATCH.md`, `ANGULAR_AUDIT.md`, `ANGULAR_PROOFS.md`, `code/angular/`, `code/ladder/`, `code/angmatch/`, `code/angproofs/` (l. 2021–2022, 2066). None of these is in the published `data/` directory. For a result advertised in the abstract, the proof must be in the Supplement, or in published files.
6. **Prominence.** The angular result appears in the abstract but nowhere in the Intro's list of results (l. 57–60), the Result (l. 209–248), Fig. 1 or the Discussion (which mentions only "the algebra of energy detectors" as an ingredient). Either it is a headline result, and should be stated in the Result and drawn in Fig. 1(a) (a 1/Δ_g curve fits easily in the plotted range), or it should leave the abstract.
7. **Scaling.** τ ≲ (160/9)/Δ_g is three powers of Δ_g weaker than the Λ⁻⁴ scaling that the Letter itself credits to causality (l. 121–122). The Letter should say that this bound has an explicit constant but non-optimal scaling, so that readers do not compare 0.018 with 2.9×10⁻⁸ as two "bounds" of the same kind.

**Fix.** Replace l. 307–310 with:

> "… this gives, under (i) and (ii) and at leading order in 1/c_T but without assumption (iii), τ ≤ (160/9)(Δ_g−2)⁻¹[1−49/(2Δ_g+9)²]^{−1/2}, i.e. τ ≤ 0.0178 at Δ_g = 1000. Besides (i) and (ii) it uses: that the special conformal Ward identity and light-transform conservation hold as operator identities on the states involved; existence and light-cone convergence of the mixed detector commutator (J₀ < 2); and the absence of scalars with Δ ≤ 1 in T×T. The bound scales as 1/Δ_g, weaker than the Λ⁻⁴ expected from causality, but it has an explicit constant."

Then add the result to the Result and to Fig. 1(a); give it a grade in S1; put the all-J proofs in the Supplement; and remove "hostile audit" and the `.md` filenames.

### M5. Item (1) is called "derived", but the Supplement's own weakest-link rule makes it a sketch for parity-violating theories
**Where:** `letter.tex` l. 240 ("The cubic implication (1) is derived"), l. 105–108; `supplement.tex` l. 76, 301, 360, 645 (the L1 row: "the parity-odd, odd-twist [TT] double zeros [D-sketch]").

**Problem.** S1 says that "a chain carries the grade of its weakest link". Item (1) is graded [D] "given the [TT] double-zero statement of L1 (its parity-odd, odd-twist part [D-sketch])". So for CFTs where the parity-odd sector matters, which are exactly those with α_T ≠ 0, implication (1) is net [D-sketch]. The Letter advertises the α_T bound as a new element (l. 334) while calling (1) "derived".

The Supplement also grades the AdS-side parity block structure for parity-violating theories as [D-sketch] (S2.C, *Parity-violating theories*), and the α_T sign convention as [I].

**Fix.** Either prove the parity-odd double-zero statement, or write in the Letter: "The cubic implication (1) is derived for parity-invariant theories; for parity-violating theories one step (double zeros of parity-odd, odd-twist double traces) is sketched." Also restate the net grade in S3/S4.

### M6. The large-spin "theorems" are presented more strongly than their premises allow
**Where:** `letter.tex` l. 265–272, 95–97, Fig. 2 caption l. 292–295; `supplement.tex` S7.F, Table S9.

1. The Letter says that one explicit term "uses less than half of the margin on near-gap states up to spin 344", and that "this gives two theorems, conditional on the uniform decay stated in (iii)". But the "less than half" is d_χ ≤ 0.489H (S7.F). That is the *entire* budget except the edge integral: window, boundary, cutoff complement and top terms. It is not the cost of "one explicit term", which is the monopole term, and that term cancels at O(1/m²) (S7.D(iv)).
2. The excluded edge integral (EDGE-R) is a model premise [M] for the graviton. At the thinnest sampled state the exactly computable scalar analogue reaches 0.405 against a tolerance of 0.406 over the lemma's cutoff grid (Table S9 reading). That is essentially zero margin.
3. Theorem W further depends on P0, partly [NE] (but see M3); on P1′ ([D-sketch] beyond the leading-Regge kernel); on a boundary factor of 2 ([D-sketch]; "a factor-2 error would raise d/H to ≈ 0.9 at J = 300, 344"); on a "+− transfer [that] is a model convention [M]"; and on "between samples [NE]". Theorem L′ is graded [D given …], not [P], with P′_L(i) [M] above J = 344.

Calling these "theorems", with only the decay assumption named in the main-text paragraph, overstates them.

**Fix.** Replace l. 268–272 with:

> "… plus one explicit term. With a contour argument for higher orders this gives two conditional statements: positivity on 47 sampled near-gap states of spin 26 ≤ J ≤ 344, and on states of spin below 10⁶ whose impact parameter exceeds about a third of the AdS radius. Both assume a uniform decay of beyond-Regge corrections (iii) and take a few magnitudes from models; for the first, the computed costs use < 0.49 of the margin, and the remainder must absorb an edge integral of the graviton that we estimate only in a model."

Use "Propositions" rather than "Theorems" in the Supplement, or give "Theorem" an explicit premise list each time it is cited.

### M7. The claim "three assumptions" is incomplete
**Where:** `letter.tex` l. 81 ("Our results rest on three assumptions") vs. l. 210–211 ("and the analytic inputs listed in the SM").

**Problem.** The "analytic inputs" include:
- L7(b), graded [A] (no conspiracy in the heavy sum and contact tower), which the quartic bound uses;
- R1 (J₀ ≤ 1), [Q], quoted with a stated KKSDZ loophole;
- κ > 0 [D-sketch + M] and the scalar condition;
- an H17 band graded [D-sketch + M] for item (1)'s constant (Table S2, row L2) but [D] in I3. That is itself inconsistent; see m4.

**Fix.** At l. 81 write: "Our results rest on three assumptions, (i)–(iii), together with standard inputs listed and graded in the Supplemental Material, the strongest of which is an assumption that light double-trace data contain no conspiracy (used only for g̃₄ = g₄)." Alternatively, promote L7(b) to an explicit (iv) for the quartic bound.

### M8. Larger gaps: the Letter and Fig. 1 overstate the support
**Where:** `letter.tex` l. 172–175 (Fig. 1 caption), l. 274–284, l. 296 (Fig. 2 caption); `supplement.tex` S8.D, S8.G, Table S17.

1. The Fig. 1 caption says that above 1000 "the solid cubic curve is supported, under the no-growth assumption …, at the strength of the evidence at 1000". The Supplement says instead:
   - a near-gap band 33 ≤ J < J_W (J_W = 52 at 2000, 132 at 5000) "carries no grade at any Δ_g > 1000" and "grows with the gap";
   - the large-spin premises are "carried in AdS units [M]".

   The Letter never mentions the uncovered band. The Supplement's S1 still maps the phrase "near-gap states with spins from 33 to about 0.026Λ", so it was removed from the Letter but not from the legend.
2. The no-growth test "at 2000 and 5000" partly uses K₂ (Laplace) rows whose "chart-A gap independence is largely structural" (S1, l. ~93–95). A test that is gap-independent by construction cannot test gap independence. The independent evidence is the exact scalar analogue and the direct rows at 2000. Say so.
3. Fig. 2 caption, l. 296: "At larger gaps the near-gap corrections shrink like Δ_g⁻²". This is measured only for J ≤ 16 (S8.C); at near-gap large spin, in chart B, "pointwise corrections do not shrink".

**Fix.** Add to l. 284: "One near-gap band (33 ≤ J ≲ 0.026Λ) is reached by neither transfer, so at Δ_g > 1000 positivity there is a hypothesis." In the Fig. 1 caption, replace "supported … at the strength of the evidence at 1000" with "supported at most at the strength of the evidence at 1000, under the no-growth assumption and with one near-gap band uncovered". In the Fig. 2 caption, write "for J ≤ 16 the near-gap corrections shrink like Δ_g⁻²".

### M9. Sec. S1 is out of date
**Where:** `supplement.tex` l. 25–110 (Sec. S1); `letter.tex` l. 72–73 ("the Supplemental Material grades every statement against this legend (Sec. S1)").

**Problem.** S1 maps about ten Letter phrases that are no longer in the Letter (checked by grep):
- "One edge integral of the graviton, estimated in a model and met by an exactly computable scalar analogue" (l. 70);
- "the other magnitudes the theorems use are model estimates with large margins";
- "a hypothesis supported by a sketched argument" (l. 75);
- "The remainder was checked on sampled states up to mass 5Λ" (l. 83);
- "The states examined while the functionals were being chosen" (l. 84);
- "Twenty spins up to 29 … the checks close for every computed spin up to 21" (l. 88). The Letter now says "twenty-two spins … every computed state closes except the two lightest of spin 29", a different count;
- "with a safety factor above ten" (l. 96);
- "a certified remainder bound supplements it for the lowest spins" (l. 97);
- "the strip just above each new gap" (l. 99);
- "near-gap states with spins from 33 to about 0.026Λ" (l. 99–100);
- "with the large-spin premises fixed in AdS units, the theorems' margins grow".

Conversely, S1 does not map the T0 sentences, "computer-certified bounds … one fifth", the angular route (l. 301–310), or "Separate code reproduces every constant". The Letter's promise is therefore false at present.

**Fix.** Regenerate S1 from the current Letter text, one row per claim sentence, and add rows for T0 and the angular route. Remove the rows for phrases that were deleted.

### M10. Reproducibility: certified and proved claims rely on code and notes that are not published
**Where:** `supplement.tex` l. 1222 (`code/t0cert/`), l. 2021–2022, 2066 (`ANGULAR_*.md`, `code/angular/`, `code/ladder/`, `code/angmatch/`, `code/angproofs/`), l. 257 (`README.md` errata); `letter.tex` l. 175, 298 (`fig/make_fig.py`, `fig/make_fig2.py`). The published `prl/data/` directory holds only weights, kernels and two scripts.

**Problem.** A "certified" (computer-assisted) proof is checkable only if its code is available. The T0 certificate is a headline claim of the abstract.

**Fix.** Publish `code/t0cert/` and the angular-route proofs and code with the Supplement (or in a cited Zenodo or GitHub archive), and cite them as such. Remove internal `.md` file names, or turn them into published appendices.

---

## 2. Presentation problems that would hurt at PRL

**P1. The evidence summary is repeated about five times.** It appears in the Intro (l. 63–67), Setting (iii) (l. 91–98), "Is the functional positive?" (l. 251–263), "Large spin" (l. 265–272) and the Fig. 2 caption (l. 288–298). The sentence "at Δ_g = 300 positivity … is an open hypothesis" appears in both the Fig. 1 caption and l. 263.
*Fix:* reduce the Intro version to one sentence ("We prove positivity in flat space, certify it at leading Regge order in AdS, and give graded evidence beyond (Fig. 2)"), and reduce (iii) to its definition plus "evidence: below". This saves about 150 words.

**P2. Paragraph (iii) is overlong** (l. 90–101): a run-on of about 120 words that mixes four grades. See P1.

**P3. "Causality alone allows order-one anisotropy"** (l. 61–62) is incorrect. Causality *with* a large gap already forces anisotropy ~Λ⁻⁴, as l. 121–122 itself says. What allows A ≤ 1 is the conformal-collider (ANEC) bound with no gap.
*Fix:* "… where the conformal-collider bound, which holds without a gap, allows order-one anisotropy."

**P4. Result (2a), l. 237–240, contradicts itself.** "(2a) … if no scalar lighter than 0.46Λ couples" and then "a scalar lighter than 0.46Λ does not invalidate it, but adds a term proportional to its squared coupling."
*Fix:* "a scalar lighter than 0.46Λ weakens (2a) by a term proportional to its squared coupling (Supplemental Material); (2b) removes this dependence."

**P5. Undefined or ambiguous symbols in the Letter:**
- κ (l. 99, 223): define it as "the coefficient of −τ²Λ⁸/57600 on the light side of Ω^Q_up (κ = 0.11 in flat space)";
- the relation ν = Λp between celestial frequency and the weight variable (l. 189 vs. 214);
- 𝒜₄ (l. 162, 197);
- "natural normalisation" (l. 93);
- "EE" (l. 153, used before it is defined);
- "load-bearing" (l. 183);
- ψ₂ is called "the published k = 2 weight" (l. 215), but the Supplement has two k = 2 weights: the flat polynomial ψ₂(p) of NF02, and the de-shifted χ(ν). Write "ψ₂ is the flat k = 2 weight (an exact-rational polynomial, published)".

**P6. "All weights are published (at Δ_g = 2000 and 5000 as tables)"** (l. 182). The Supplement says only the continuum weights define the functional, because grid weights fail the interchange by its own no-go (S2.C, S5.D). Continuum k = 2 weights are published only at 300 and 1000; at 2000 and 5000 only truncated grids exist (Table S1, `gap/`).
*Fix:* publish the continuum weights at 2000 and 5000, or write "at 2000 and 5000 as grid tabulations of the continuum weight".

**P7. Fig. 2 does not match its caption.**
- The caption cites "(dashed line)" for the second theorem, but the rendered figure shows a filled medium-blue band labelled "theorem 2"; I see no dashed line.
- The caption names "Dark blue" and "Lightest blue" but not the medium blue of theorem 2.
- "Every spin up to 21" includes J = 1, 3, which have no states; write "every spin up to 21 that has near-gap states".
- The figure does not show T0, the main *certified* AdS statement. Add a sentence: "Throughout the plane the leading-Regge part is certified ≥ 0.1999 N_s; shading grades the beyond-Regge corrections."

**P8. Fig. 1.** The angular bound is not drawn (see M4.6). "Faded lines: the implication alone" is opaque; write "faded: value of the implication where positivity is not examined". Remove the script name from the caption.

**P9. Internal language and file names.**
- *In the Letter:* `fig/make_fig.py` (l. 175) and `fig/make_fig2.py` (l. 298). Remove both.
- *In the Supplement:* "chart A/B", "K₂ rows", "the chain", "hybrid", "own row", "seam", "crystal integral", "fold chain", "Stage-7 grid weights" (l. 268), "Q3-full-delta0.2", "fingerprint 80e0834a…33f2", "NF02", "H_hold", "H_R′", "Λ²H17", "Ω^EE has been fixed since 20 Nov 2026, Ω^Q since 22 Nov" (l. 928), "independent hostile audit" (l. 2010), "fourth checker", the `README.md` errata, and the `.md` files (M10).

   The Supplement needs a glossary at minimum, and preferably a rewrite of S6–S8 for outside readers. Date-stamped project history does not belong in a journal supplement.

**P10. Independence claims are unexplained.** "reproduced independently" (l. 258), "Separate code …" (l. 242) and S10's "two reproducers … a third checker … a fourth checker" do not say who or what performed the checks. The acknowledgments say an AI system assisted. Say whether the "independent" checks were separate AI runs with no shared code, and avoid "independent" where the evidence is a byte-for-byte rerun (S7.B, l. 1218: "matched every output byte for byte"). Matching byte for byte shows the same code was rerun, not that it was independently reimplemented.

**P11. Selection bias is disclosed only in the Supplement.** The functionals were chosen after seeing heavy-state failures (S6 Remark (i)). One sentence in the Letter would pre-empt the obvious referee question: "The EE structure and minimal-norm gauge were chosen after failures of the standard choices were seen; the evidence quoted excludes those in-sample states."

**P12. The ABJM check is weak as a test.** The abstract claims "ABJM theory lies inside our bound". Any theory with modest positive g₄ lies at about 3% of an upper bound of 92.8. The body says this honestly ("not evidence for positivity"); the abstract should not present it as a result.

**P13. Section headings.** "Why it matters" actually states the headline result; rename it "Main result in physical terms". Rename "Is the functional positive?" to "Evidence for positivity".

**P14.** "Defining a scale by |ĝ₃| = M_*⁻⁴ … cannot set in below about a third of the gap" (l. 125–126): 0.30 is "about 0.3", not "about a third". Minor, but PRL readers will check.

**P15.** "J₀ < 3/2 (k = 2, 3, 5) or J₀ < 5/2 (k = 4)" (l. 184–185) looks paradoxical: k = 5 is more restrictive than k = 4. Add "(set by the analytic structure of each weight)".

---

## 3. Typos, LaTeX, references, notation

| # | Location | Problem | Fix |
|---|---|---|---|
| T1 | `supplement.tex` l. 2006–2007 | "explicit kernel envelopes give hence μ((L,∞)) = …" is garbled | "… and explicit kernel envelopes give μ((L,∞)) = O(L⁻³) …" |
| T2 | `supplement.tex` l. 339 | "(Sec. S3, I3)": I1–I6 are listed in S4, not S3 | "(I3 in Sec. S4)" |
| T3 | `supplement.tex` l. 1936 (Table S18) | τ₀ = 3.13206×10⁻⁶ / 2.89126×10⁻⁸ disagrees with S4 (3.13262×10⁻⁶ / 2.89131×10⁻⁸). These are the K = 1 values (l. 709) but are unlabelled | label the column "τ₀ at K = 1" |
| T4 | `letter.tex` l. 221 | "τ₀ ≤ 2.90×10⁻⁸" while the Supplement certifies 2.8914×10⁻⁸ | quote "τ₀ ≤ 2.8914×10⁻⁸ (certified)", or explain the rounding |
| T5 | `supplement.tex` S7.B vs S7.F | H means the hybrid action in S7.B and the 2×2 budget f₊₊f₊₋/(f₊₊+f₊₋) in S7.F | rename one |
| T6 | `supplement.tex` l. 1221 vs S2.C | κ is the quartic light coefficient (0.109908), and also "κ > 0" for the OO block in T0; "EO cancellation" is undefined | define, or rename (e.g. κ_OO) |
| T7 | `supplement.tex` S7.B | b̂ and Λ_min are undefined | define b̂ = b/Λ (if that is meant) and give Λ_min |
| T8 | `supplement.tex` l. 1538 | §S8.B "What holds uniformly in the gap" contains only a float, so in the PDF it reads as an empty section | add a one-line lead-in referencing Table S12 |
| T9 | `supplement.tex` H17 | graded [D] in I3 (l. 379–380) but [D-sketch + M] in Table S2 row L2 (l. 302) and in S3/S4 (l. 322, 361) | the closed-form value (D) and the t₄/α_T asymmetry band (M) appear to be different objects; give them different names |
| T10 | `letter.log` | "A float is stuck (cannot be placed) on input line 200"; harmless in the current build but fragile | move Fig. 1 earlier or use `[t]` |
| T11 | `letter.tex` l. 197 | "a direct numerical scan follows the predicted trend": the Supplement says the scan is pre-asymptotic (slopes 2.5–2.8 vs. 2) | "a direct numerical scan is consistent with it but pre-asymptotic" |
| T12 | `supplement.tex` S11, free-field check | "the bound itself does not apply there: no gap, Δ_φ < 1": the free fermion's ψ̄ψ has Δ = 2; the failure there is only the missing gap | "(no gap; and for the scalar, Δ_{φ²} = 1)" |
| T13 | `letter.tex` l. 119–120 | "⟨TTT⟩ is that of equal numbers n_s, n_f …": α_T ≠ 0 cannot be produced by free fields | "the parity-even part of ⟨TTT⟩ is that of …, and α_T ≤ τ" |
| T14 | both | British spellings (neighbouring, normalised, localisation); APS will convert | none needed |
| T15 | `letter.blg` | "empty year in SM" | harmless for the SM placeholder |

---

## 4. Length

- **Abstract:** about 566 characters rendered (571 in source), under the 600-character limit. My suggested replacement in M1 is about 550.
- **Main text:**

  | Component | Word equivalents |
  |---|---|
  | Body text, excluding captions, display equations and references (inline math counted as one word each) | about 2,860 |
  | Fig. 1 caption | 135 |
  | Fig. 2 caption | 146 |
  | Three display equations, about 3 rows at 16 words each | about 50 |
  | Fig. 1 (aspect 0.94, so 150/0.94 + 20) | about 180 |
  | Fig. 2 (aspect 1.125, so 150/1.125 + 20) | about 155 |
  | Acknowledgment | about 30 |
  | **Total** | **about 3,550** |

  That is roughly 95% of the 3,750 limit: compliant but at risk. Any fix that adds text, for example M4 (adding the angular result to the Result) or M6 and M7, must be paid for. Cutting the repeated evidence summaries (P1) frees about 150–200 words, and shortening the Fig. 2 caption about 40 more.

---

## Summary of the top 10 items

1. **M1.** The abstract implies positivity is proved; it is assumption (iii). It also omits leading order in 1/c_T, and "two-percent bound at every spin" is ambiguous. A replacement abstract is supplied.
2. **M2.** The T0 floor is stated inconsistently. "Keeps one fifth of the flat margin" (Intro) is wrong: it keeps about 99.95%. "At least one fifth … Δ_g ≥ 300" is false (0.199). "Every heavy state" ignores the Supplement's scope (target weights only; published weights for β ≲ 200 only; OO block needs extra inputs). "Every gap ≥ 1000" conflicts with "wherever the weight exists (≤ 5000)".
3. **M3.** The Supplement contradicts T0: S7 intro, S7.C, S7.F (P0), the Table S9 reading, Table S11 and Table S17 still grade the leading-Regge floor [NE] or "at sampled masses".
4. **M4.** The angular route's inputs are listed differently in the Letter (3) and the Supplement (4; operator identities and existence of the mixed commutator are missing from the Letter). The route uses (i), (ii) and leading order in 1/c_T without saying so. "At every spin" rests on helicity norms checked to J ≤ 30. It is ungraded in S1, its proofs sit in unpublished `.md` files, it is missing from the Result, Fig. 1 and the Discussion, and its 1/Δ_g scaling is not flagged as non-optimal. Say "without assumption (iii)" instead of "no positivity assumption".
5. **M5.** "The cubic implication (1) is derived" conflicts with the Supplement's weakest-link rule: the parity-odd double-zero step is [D-sketch], which matters exactly for α_T.
6. **M6.** The large-spin "theorems": "less than half of the margin" is the whole non-edge budget, not "one explicit term". The excluded edge integral is a model premise at about zero margin for the scalar proxy (0.405 vs. 0.406). Other premises are [D-sketch] or [M]. Call them conditional propositions and list the premises.
7. **M9.** The S1 legend is stale: about ten mapped Letter phrases no longer exist, and T0 and the angular route are unmapped. The Letter's promise that the Supplement "grades every statement" is currently false.
8. **M8.** Larger gaps: the uncovered band 33 ≤ J < J_W is not mentioned in the Letter. The no-growth test partly uses rows whose gap independence is structural. "Corrections shrink like Δ_g⁻²" holds only for J ≤ 16.
9. **M10 and P9.** Certified and proved claims depend on unpublished code and notes (`code/t0cert/`, `ANGULAR_*.md`). Script names appear in the Letter's captions. The Supplement is full of internal language ("chart A", "K₂ rows", "hostile audit", "fixed since 20 Nov 2026").
10. **M7, P3 and P4.** "Three assumptions" omits the [A] input L7(b) and other graded inputs. "Causality alone allows order-one anisotropy" is incorrect. Result (2a) says a light scalar "does not invalidate" a bound that assumes there is none.
