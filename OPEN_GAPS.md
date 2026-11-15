# OPEN GAPS

## G1. Hypothesis H (identification of the proved functional with the HRT entropy)

Theorem 2 is proved for the functional `min_c ½(ℓ_L(c) + ℓ_R(c))` over all topological cut classes. That this
equals the HRT entropy is taken from the literature. A derivation from Wall's maximin theorems is sketched in
DERIVATIONS.md §5; it is at the level of rigour of those theorems (existence and extremality of maximin
surfaces), not a self-contained proof. Specifically unproved here:

* existence of a maximin surface for every boundary region in these spacetimes;
* that the null congruence from each geodesic component reaches every Cauchy slice in a curve of the same
  class (used to show maximin ≤ length of every class).

What *is* proved towards H: geodesics of an embedded class are pairwise spacelike separated (PROOFS.md
Lemma 6.3'); non-embedded collections never beat embedded ones (Lemma 6.4). What remains assumed (H2): a
pairwise spacelike-separated collection of geodesics bounding `A` topologically also bounds a *spacelike*
region with `A` (extension of an acausal set to a Cauchy surface).

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

* Two rounds by one automated referee (same model family as the author; not a human expert; nothing is
  formalised). Round 1: no fatal error, repairs listed. Round 2: no fatal error; the referee regards as
  established Lemma A, Facts 3.1–3.2, Lemma B, Theorem 1, Lemma G, Lemma F, Theorem 3 for whole-boundary
  parties, the length formula and Lemma 6.4; as conditional the two statements about HRT; as not established
  Theorem 3 for regions with endpoints, H1 and H2.
* A third pass checked each round-2 repair (12 items): all made correctly, no new error, status table judged
  accurate. Its remaining wording points (definition of the domain `D̃`, "limit set", hypotheses of Lemma 6.3',
  two additions to Lemma 6.4, one justification in Theorem 3) are applied in revision 4 and were not re-read.
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
| 27 | H2 | every embedded class is admissible | S — **the main open item** |
| 28 | H3 | the homology condition of arXiv:2407.07165 implies the colouring hypothesis of Lemma 6.4 | C |
| 31 | §5 | closedness and graph characterisation of the cone: cited | C |

(Items 11, 13, 14, 17, 18, 22, 24, 29, 30 of the referee's list were errors or omissions in revision 2 and are
repaired in revision 3.)

## G6b. Where a counterexample could conceivably live (referee's list)

* **H2:** an embedded class whose geodesics do not bound a spacelike region. Then the proved functional would
  be below the HRT entropy for some region and the HRT conclusion could fail. Judged unlikely given Lemma 6.3'.
* **H1:** an admissible extremal configuration that is not an embedded cut class.
* **Marking:** a twist mismatch between spacetime classes and the identification of the two surfaces. Theorem 1
  would still hold for the mis-marked functional, but that functional would not be HRT.
* **Unbalanced inequalities at finite cut-off:** the constant in Theorem 1(2) can be negative; the statement is
  only that the difference diverges to `+∞`.
* **Out of scope by assumption:** cusps, conical defects, non-orientable quotients, endpoints not on one slice,
  matter.

## G7. Secondary target (SYK Gaussianity)

Not started. See APPROACHES.md.
