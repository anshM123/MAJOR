# Hostile PRL referee report on the covariant entropy-cone Letter (28 Nov 2026)

The referee read only `prl_hrt/letter.pdf` and `prl_hrt/supplement.pdf`, plus the arXiv literature. Its check scripts are in this directory. Downloaded papers in `lit/` are not committed.

## Recommendation

**Reject for PRL and transfer to JHEP** (CMP or ATMP would also suit the vacuum theorem). If the editor keeps it, the minimum is a major revision addressing M1–M6.

- **Confidence:** about 75% that it fails PRL's importance and broad-interest bar as written; about 80% that the vacuum results are correct.
- **Advice to the authors:** the referee would not tell them to submit to PRL now.
- **Single most useful change:** rebuild the Letter around the vacuum theorem alone, with
  1. an abstract under 600 characters that states the regulator and the scope;
  2. an accessible account of the mechanism: positive current, cut-and-paste on the quotient, Bell pairs;
  3. a real argument for relevance beyond 2+1 vacuum. The best candidates are the regulator-independent finite part G and the corollary that a rotating wormhole has the entropy vector of a time-symmetric wormhole.

  Move Theorems 4 and 5, the smooth-slice proposition and the NEC material to a companion long paper.

## Validity

The referee spot-checked the whole unconditional chain: Lemma S3.1, Prop. S4.3, Cor. S5.1, Thm. S5.2, Prop. S6.1, Thms. S7.6–S7.7, Lemmas S8.1–S8.3, Prop. S8.4, Thm. S9.1 and Prop. S10.1. Nothing it found invalidates the chain.

Independent checks:
- **Prop. S4.3 endpoint integrals:** confirmed to O(Y⁻²).
- **Separation Lemma:** minimum ratio 1.1995 ≥ 1 on random acausal circles.
- **Explicit strip, Eqs. (S27)–(S28):** the identities are exact.
- **Rotating-BTZ HRT solver** (π₁ = Z only):
  - 400 configurations at N = 3 and 80 at N = 5, with no violation; the minimum facet value is −3.6e−16;
  - G is ε-independent to 2e−14;
  - the common-class constraint changes the minimiser in about 10% of regions;
  - the control functionals fail, as they should.

## Major issues

### M1. A special case with an expected answer

The answer was strongly expected. GGHH's numerics already pointed this way, as do the majorization papers 2508.21823 and 2601.09989, which are not cited. The unconditional result is confined to 2+1 vacuum gravity, which has no local degrees of freedom. The Letter needs a relevance argument beyond that setting.

### M2. Abstract and headline over-read the results

- **(a)** "We settle it unconditionally": only a special case is settled.
- **(b)** "For every topology" contradicts the stated exclusions: cusps, conical defects, non-orientable quotients, and pure AdS₃ with |P| ≤ 2.
- **(c)** "At every sufficiently small cutoff" holds only for the horocyclic regulator.
  - A Fefferman–Graham cutoff adds O(ε²) terms that are not of Bell-pair form.
  - Exact membership can then fail for saturated inequalities.
  - The threshold ε₁ also depends on the regulator and on the state.
- **(d)** Calling (JM) "the maximin hypothesis that underlies Wall's covariant proof" is inaccurate; see M3.
- **(e)** "The NEC is necessary" is Callan–He–Headrick (1204.2309). The paper's own thin-shell solver does not exhibit it, and the single example shown is on a near-null slice in simply connected planar Vaidya.
- **(f)** The strip theorem is conditional and should not be in the headline.
- **(g)** The abstract is about 1340 characters; PRL's limit is 600.

### M3. Theorem 4's hypotheses are stronger than, and different from, Wall's

- (JM) restricts to homotopy classes, not homology, and is required in every finite cover.
- (R2) demands a modulus of continuity uniform over all covers and all multi-classes.
- "Exactly as conditional as Wall" (Remark S12.1) must be withdrawn.
- Precompactness (K) fails generically for NEC wormholes.
- "Holds in vacuum" needs the qualifier "given (R1) and (R4)".
- "Positivity locally without maximin" (Letter ll. 312 and 361) contradicts SM ll. 2009–2013.
- Letter l. 303 writes L(h) where the SM's Θ(h) is meant.
- "Any spatial topology" needs the list of restrictions.

### M4. Theorem 5 is close to vacuous given (H4)

The reduction lands in Czech–Dong's simply connected setting. Given the assumed maximin content, the graph theorem is the classical circular-split (Kalmanson) fact. Calling it an extension of Bao–Mezei is an overclaim. Remove it from the Letter.

### M5. Theorem 3 is a corollary

It follows from Theorem 2, closedness of the cone, and Bao et al.

### M6. Spot-checks of the vacuum proof

The proof is judged essentially complete. Each item below needs a sentence or a fix.

- **Lemma S3.1:** add a sentence each on funnel cells and on the reduced-class minimum.
- **Prop. S4.3:** the measure-zero limit set is why cusps and extremal horizons are excluded; say so.
- **Theorem H scope:** state that Wall's achronal slices and homology regions off a Cauchy slice are excluded.
- **Prop. S10.1 (smooth slices):** not needed for the cone result; move it out of the Letter.
- **Numerics:** uncertified and partial.
  - Word lengths are only searched up to 3–4.
  - Closed geodesics never entered a minimiser on the tori.
  - The explicit graph was built only on the disc.
  - The code location is a placeholder.

## Minor issues (abridged)

1. Cite Maxfield 1412.0687 in the text for the length formula. Cite Bonahon and Glorieux 1503.09067 for currents.
2. Refs. [53]–[56] are cited nowhere.
3. **Error in SM Remark "Cutoff scheme" (ll. 565–568).** G does *not* depend on the horoball normalisation; only t_p does. The referee confirmed this numerically: rescaling a horoball by 0.1–10 changes G by less than 2e−14. G is a canonical, conformal-frame-independent finite part, and the paper should say so.
4. Cite 2508.21823 and 2601.09989.
5. Make the remark about "the recent review" neutral.
6. Clarify AdS₃ = PSL(2,R) and the positive lift.
7. Write "non-elementary", not "non-abelian". Time symmetry also needs a reflection-symmetric boundary slice.
8. Theorem 2(ii): state that the divergence uses Q(b_p) ≥ 0 and Q(b_p) ≠ 0 for some p.
9. Define "balanced at every entangling point" at first use and relate it to superbalance.
10. "Without any maximin input" should read "given (R1) and (R4)".
11. Put the fixed boundary slice into the Theorem 4 statement.
12. Mention the possible ν-atoms at peripheral axes in the Letter.
13. Write out the isotropic m(m+1)(kφ)² term of Eq. (5) explicitly.
14. Affiliations, acknowledgments, SM URL and code location are placeholders.
15. Give the real fractions for control 2: the pants N = 4* run fails in only 6 of 30.
16. The status table should cite Callan–He–Headrick for NEC necessity.
17. The symbol c is used for both cut classes and the central charge.

## Presentation

- **Readability:** the Letter is jargon-dense. The core idea fits in three sentences:
  1. AdS₃ lengths are averages of two hyperbolic lengths.
  2. Each hyperbolic length counts crossings with a positive weight.
  3. Static min-cut arguments need only positive crossing weights, plus one Bell pair per endpoint.
- **Best paragraph:** "Why the common class was the obstacle".
- **Length:** about 3680 word-equivalents, at the limit, with five theorems each only sketched.
- **Figures:** Fig. 2 is overloaded; panel (b) is the useful one.
- **Supplement:** honest and well organised, but at about 27,000 words it is a monograph.
