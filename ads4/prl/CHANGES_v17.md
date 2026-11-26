# v17 — positivity of Omega^EE on every heavy state at Delta_g = 1000 (Theorem FP)

Status: drafted while the final checks run; NOT pushed until the portable-GMP replay of the coupled-cap certificate
(172 shards) and the referee of sliver_fullcov/PROOF.md both report clean (ANGULAR_ROUTE §127).

## Letter
- Abstract: the cubic bound no longer rests on an assumed positivity property; positivity is established in AdS on every heavy
  state at Delta_g = 1000 by a computer-assisted proof. The quartic bound is stated as resting on a further positivity assumption.
- Intro: positivity "certified at leading Regge order ... graded evidence beyond" -> established on every heavy state for the
  cubic functional at 1000; quartic functional remains an assumption. Legend: "certified" includes directed rounding arithmetic.
- (iii): for the cubic functional at Delta_g = 1000 it is not an assumption (flat proof, T0 floor, full-action certificate).
- "Why it matters": under (i) and (ii); gap Delta_g >= 1000.
- Result: (iii) proved for item (1) at Delta_g = 1000.
- "Is the functional positive?" rewritten around the certificate cover; "Large spin" paragraph removed (its two conditional
  statements are now a cross-check, mentioned in one sentence).
- "Larger gaps": the certified bound at 1000 holds for every larger gap; sharper bounds at larger gaps need positivity there.
- Discussion: (i),(ii) for the cubic bound; new element "a proof that a dispersive functional is positive on the entire heavy
  spectrum of AdS at finite gap".
- Fig. 1 caption: above 1000 the solid curve assumes positivity at its own gap.
- Fig. 2 replaced: fig/fig2_certificate.pdf (make_fig2_fp.py), the certificate cover in the (L, twist) plane.
- Word count essentially unchanged (2970 -> ~2980, text only).

## Supplement
- New Sec. S11 (sec_fullpos.tex, \input before the angular route): Theorem FP, the exact action (kept from the earlier draft),
  coordinates and the cover (Table tab:fpcover, completeness argument with verified a* thresholds), fixed spin, high twist,
  L >= 3.65 (tail, T', T''), 2.3 <= L <= 3.65 (signed transfer, norming, base multipliers, three circles), 0.05 <= c <= 0.9 (direct
  prefix, coupled cap), inputs/reproduction/code.
- S1 mapping: (iii) for Omega^EE at 1000 is Theorem FP; positivity-paragraph mapping rewritten; Fig. 2 -> Table tab:fpcover.
- Full statement: (Pos) for Omega^EE at 1000 is Theorem FP, not an assumption.
- Sec. S6 header: superseded at 1000 for Omega^EE; remains evidence elsewhere and for Omega^Q_up.
- Scope: "an unconditional bound" replaced by the precise exclusions (A1', A2, inputs; other gaps; Omega^Q_up).

## Not changed
- Acknowledgment (as instructed). Quartic bound and its assumptions. Delta_g = 300 statements.
