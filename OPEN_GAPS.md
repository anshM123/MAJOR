# OPEN GAPS

## G1. The link to the HRT formula — Theorem H, with smooth-slice regularity closed in revision 8

Since revision 5 of PROOFS.md the HRT entropy is defined in §6.3 by the standard prescription (minimal
extremal curve anchored on `∂A`, subject to the spacelike homology condition on a bulk slice) and Theorem H
proves that it equals the minimum over topological cut classes. The former assumptions H1 and H2 are now the
two halves of that proof. What the statement still depends on:

* **The definition of admissible slice.** Slices satisfy (S1)–(S3). Proposition S+ (PROOFS §6.6b,
  [SMOOTH_SLICES.md](SMOOTH_SLICES.md)) now proves that the identical minimizing geodesic cut lies on a
  compact smooth spacelike Cauchy slice of the conformal completion, with the correct homology. Its induced
  conformal metric has a positive uniform margin relative to any fixed smooth auxiliary metric on the
  quotient. Thus Theorem H also holds with smooth uniformly spacelike slices for **all regions, including
  endpoints**. The original pleated construction remains nonsmooth and asymptotically lightlike at rotating
  axes; the new proof replaces the slice while fixing the cut.
* **One cited topological fact** in Lemma T, Step 3 (a compact surface with prescribed vertices, not one of
  four small exceptions, has a triangulation with those vertices). The spiralling ("spinning") step, which in
  revision 5 was a second citation, is proved in the text since revision 6 (Lemma T, Step 4).
* **Equivalence with "acausal Cauchy slice of the spacetime with boundary"**: proved in revision 7 (PROOFS
  §6.3b, Theorems C and E), with an explicit conformal model and the signed causal criterion. Awaiting
  independent review.
* **Smooth slices: proved in revision 8**, using a smooth conformal double. Anchored geodesic arcs reflect
  into compact smooth spacelike circles; folding proves that the double is globally hyperbolic and the
  doubled cut acausal. Bernal–Sánchez Theorem 1.1 / Remark 4.14 extends a compact acausal ribbon, and its
  restriction is the desired boundary Cauchy slice. A peripheral-homology argument proves compactness, and
  a boundary-tangent timelike flow transports the homology region while fixing the cut. Boundary traces
  may change away from the endpoints, exactly as §6.3 already permits. The uniform constant depends on the
  geometry and cut; a universal-cover fixed-frame or moduli-independent constant is not asserted.
* **Achronal slices with null pieces** (Wall's class): not treated. The lower bound `≥` would need the
  Jordan-curve argument on an achronal disc, which looks unproblematic, but it is not written.
* **The cited description of the spacetimes** (§6.1): proper discontinuity on the invisible domain, structure
  of the asymptotic boundaries.
* **Independent review** of §§6.3–6.9: see G6.

## G2. Matter (non-vacuum 2+1 spacetimes), multiply connected

Not covered. The proof uses that each chiral length is a Liouville integral of intersection numbers. A
possible replacement:

> **Reduction (sketch, not a theorem).** Let `M` be a globally hyperbolic 2+1 spacetime with universal cover
> `M̃`, and let `D(a,b)` be the regulated length of the extremal geodesic of `M̃` between lifted boundary points
> `a, b`. Suppose (i) `D(a,c) + D(b,d) ≥ D(a,d) + D(b,c)` whenever `a,b,c,d` are cyclically ordered on the
> boundary circle of a Cauchy surface of `M̃` ("strong subadditivity on the universal cover"), and (ii) the
> length of the closed extremal geodesic in the class of `γ` equals
> `lim_n [D(a, γ^{n+1} b) − D(a, γ^n b)]`. Then the mixed second difference of `D` is a positive
> `π_1`-invariant measure `ν` on pairs of boundary points, candidate lengths are `ν`-integrals of intersection
> numbers, and Lemma A gives all contraction-map inequalities.

Missing: a proof of (i) and (ii) from the null curvature condition; a careful version of the telescoping
argument that expresses differences such as `D(a, γ b) − D(a, b)` (balanced only modulo the deck group) as
`ν`-integrals; treatment of non-unique extremal geodesics.

## G3. Higher dimensions

Not addressed. There is no analogue of the left/right decomposition, and extremal surfaces are not geodesics.

## G4. Inequalities without contraction maps

Theorem 2 covers inequalities that admit a contraction map. Whether every valid RT inequality admits one is
claimed in arXiv:2506.18086 "up to certain disclaimers" and is not universally accepted (arXiv:2601.09987:
"we will not assume this"). Not used and not resolved here.

## G5. Edge cases excluded by assumption

* Extremal horizons (parabolic holonomy, cusps): excluded by "all non-trivial elements hyperbolic".
* Conical defects / point particles (elliptic holonomy): excluded.
* Regions whose endpoints do not lie on a common boundary Cauchy slice: excluded.

## G6. Independent verification

All reviews are by automated referees of the same model family as the author. None is a human expert and
nothing is formalised.

* §§1–5 (Lemma A, Lemma B, Theorem 1, Lemma G, Lemma F, Theorem 3 for whole-boundary parties): three passes,
  no error found after repairs.
* §§6.1–6.9, revision 5 (first audit of the HRT identification). Verdict: no fatal error; Lemma 6.4,
  Proposition P and Theorem H correct for the entropy defined in §6.3. Gaps found and repaired in revision 6:
  the unsigned separation criterion was wrong inside `D̃` (explicit rotating-BTZ example) → Lemma 6.3a;
  candidates must be simple (explicit self-crossing candidate) → definition changed; uniqueness of the
  invariant geodesic for equal chiral lengths and all causal types → added; `D̃` is not convex → claim
  removed; citations for "spinning" inadequate → replaced by a proof (Lemma T, Step 4). The repairs
  themselves have **not** been re-read by that referee.
* The same audit's findings that are **not** repairable by rewording:
  - The geometric core (acausal pleated discs over an acausal curve) is prior work of Mazzoli–Viaggi,
    arXiv:2206.06946. Now cited; novelty claim reduced accordingly.
  - "No admissibility hypothesis" overstated the result at that revision: smooth slices were then known
    only for whole-boundary regions. The regularity restriction is now removed by Proposition S+
    (§6.6b), with the original geometric assumptions and boundary-trace freedom retained.
  - A slice-independent lower bound (homology region not assumed to lie in a complete slice) is not written.
* Referee's own tests: exact rational test of Lemma 6.4, 20,000 cases, 0 violations; rotating three-boundary
  wormholes with fans spiralling both ways around an axis, 8 group translates, ≈36,000 pairs, 0 violations;
  fans to depth 25 at 300 digits, 0 violations (ratio → 1 from above at the axis).
* Closed-geodesic case and Cauchy-slice assertion (second, targeted audit, prompted by the objection that an
  axis cannot be an edge of an invariant locally finite triangulation; read revision 6 and the first version
  of §6.3b). Verdict: **the revised Lemma T resolves the objection; could not be broken.** Lemma T Step 4
  (a)–(e) checked line by line, including the pieces annulus-with-one-marked-point, pants bounded by three
  closed geodesics, loop edges. Proposition P: no error (continuity and injectivity at axis points and at the
  ideal endpoints; different spiralling choices on the two sides). Theorem H `≤`: no error. Lemma K, Lemma L:
  correct. Theorems C, E and the identification with the Headrick–Hubeny–Lawrence–Rangamani condition
  (their paper read by the referee, including that their own Lemma 14 uses non-smooth slices): correct
  modulo gaps, all filled in revision 7 — closure of the lifted slice; order of two steps; proper
  discontinuity on the boundary (now an explicit cited assumption); scope `Γ` non-trivial. The fills were
  **not** re-read by a referee. Not reviewed by this referee: the `≥` half of Theorem H (reviewed by the
  first), Mazzoli–Viaggi.
* That referee's tests: rotating BTZ with fans on both sides and all four spiralling choices, 55,244,160
  pairs, 0 violations (control with a crossing triangle: 114 / 20,000); rotating one-holed-torus wormhole with
  a non-separating closed geodesic plus an arc, group translates, 60 digits, 4 × 456,489 pairs, 0 violations;
  Lemma K formula against direct evaluation, agreement to `7·10^{-9}`.
* Historical open items after those audits: smooth/uniformly spacelike slices with arcs; Wall's achronal
  class; the trivial-group formulation in §6.3b. Revision 8 resolves the first by Proposition S+, with a
  direct global-AdS version of that regularity construction. Wall's achronal class remains separate.
* **Revision 8 targeted regularity audit:** the new proof checks smooth conformal descent, folding with
  infinitely many seam contacts, both directions of causal inextendibility, opposite-copy acausality,
  compact ribbon extension, transverse restriction, peripheral compactness and homology transfer. Exact
  rotating-strip and endpoint identities are checked symbolically in `smooth_slices/`. These are automated
  checks of the written proof and local algebra, not formal or human-expert verification. The earlier
  standing geometric inputs and the rest of §§1–6 are not re-audited by this addition.
* Numerical tests: EXPERIMENTS.md. Lemma A is tested directly only on the cylinder (`π_1 = Z`).

## G6a. Referee's list: implications asserted but not proved

Seriousness: F = fatal for the affected statement if false; S = changes scope; C = cosmetic.

| # | Location | Item | Seriousness |
|---|---|---|---|
| 1 | §1.1 | convex cocompactness, `dim Λ < 1`, `Λ` null: cited only, citations unchecked | C |
| 2 | §1.2 | uniqueness of the reduced representative of a class | C |
| 3 | §1.3(c) | minimality for multiply covered closed geodesics | C |
| 4 | §1.3(d) | isotopy to geodesic position for ideal arcs | C |
| 5–8 | Lemma A | each piece's neighbourhood lies in one cell; small horoballs meet the arrangement only in rays; closure of `U'_r` at marked points; parallel-copy step for proper powers | C (all checked true by the referee) |
| 9–10 | §3 | unfolding Crofton to the quotient; existence of `ε_0` | C |
| 12 | §4 | intersection numbers preserved by non-smooth homeomorphisms | C |
| 15 | Lemma G | several terminals per party in the graph-model definition | C |
| 16 | Lemma F | normal coordinates determine the class (standard, not re-proved) | C |
| 19 | §5.4 | six missing steps | S (Theorem 3 for intervals) |
| 20 | §6.1 | proper discontinuity, convexity, structure of the asymptotic boundaries: cited | S if misquoted |
| 21 | §6.1 | lifted slice is the graph of an increasing homeomorphism; equivariant extension | C |
| 23 | §6.2(a) | the long invariant geodesic lies in the spacetime and the short one does not | F for H if false; argument now included |
| 25 | (6.1) | identification of the spacetime cut-off with horocycles | C |
| 26 | H1 | spacelike-homologous extremal surfaces are embedded cut classes | S |
| 27 | H2 | every embedded class is admissible | proved by Proposition P / Theorem H; smooth uniformly spacelike version added by Proposition S+ |
| 28 | H3 | the homology condition of arXiv:2407.07165 implies the colouring hypothesis of Lemma 6.4 | C |
| 31 | §5 | closedness and graph characterisation of the cone: cited | C |

(Items 11, 13, 14, 17, 18, 22, 24, 29, 30 of the referee's list were errors or omissions in revision 2 and are
repaired in revision 3.)

## G6b. Where a counterexample could conceivably live (referee's list)

* **H2 (historical target):** Proposition P / Theorem H prove admissibility, and Proposition S+ supplies
  smooth uniformly spacelike homology slices. A counterexample within the stated scope would therefore
  have to invalidate an identified step or standing input of those proofs.
* **H1:** an admissible extremal configuration that is not an embedded cut class.
* **Marking:** a twist mismatch between spacetime classes and the identification of the two surfaces. Theorem 1
  would still hold for the mis-marked functional, but that functional would not be HRT.
* **Unbalanced inequalities at finite cut-off:** the constant in Theorem 1(2) can be negative; the statement is
  only that the difference diverges to `+∞`.
* **Out of scope by assumption:** cusps, conical defects, non-orientable quotients, endpoints not on one slice,
  matter.

## G7. Secondary target (SYK Gaussianity)

Not started. See APPROACHES.md.
