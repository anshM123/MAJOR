# OPEN GAPS

## G1. Hypothesis H (identification of the proved functional with the HRT entropy)

Theorem 2 is proved for the functional `min_c ½(ℓ_L(c) + ℓ_R(c))` over all topological cut classes. That this
equals the HRT entropy is taken from the literature. A derivation from Wall's maximin theorems is sketched in
DERIVATIONS.md §5; it is at the level of rigour of those theorems (existence and extremality of maximin
surfaces), not a self-contained proof. Specifically unproved here:

* existence of a maximin surface for every boundary region in these spacetimes;
* that the null congruence from each geodesic component reaches every Cauchy slice in a curve of the same
  class (used to show maximin ≤ length of every class).

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

* The proof has been reviewed by one independent automated referee (report summarised in the manuscript once
  received). It has not been checked by a human expert and is not formalised.
* Numerical tests cover the conclusion of Theorem 1 in two topologies; Lemma A's construction itself is not
  exercised in code.

## G7. Secondary target (SYK Gaussianity)

Not started. See APPROACHES.md.
