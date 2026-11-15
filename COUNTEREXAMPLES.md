# COUNTEREXAMPLES, FAILED CONJECTURES, ADVERSARIAL TESTS

## 1. A natural lemma that is false as stated

**False lemma.** "Let `γ_1,…,γ_L` be cut classes and `f` a contraction map. The homotopy class of the
cut-and-paste output `∂U_r` is the same for every minimal-position arrangement of representatives of the
`γ_l`."

Why it fails: a triangle move replaces a triangular cell with bit string `x` by one with `x ⊕ 111` (on the
three curves involved). Take the coordinate function `m = f(·)_r` equal to 1 on `000, 100, 010, 001` and 0 on
`110, 011, 101, 111`. Before the move `U_r` is connected inside the disc (the triangle joins the three
weight-1 sectors) and `∂U_r` consists of three arcs, each going around a weight-2 sector. After the move the
three weight-1 sectors are disjoint and `∂U_r` consists of three arcs, each going around a weight-1 sector.
The six endpoints on the boundary of the disc are matched differently. A single coordinate of a contraction
map can have this pattern (the contraction property constrains only the weighted sum over `r`).

Consequence: the "obvious" proof (cut and paste separately in `X_L` and `X_R` and compare) does not work.
The proof in PROOFS.md never needs this lemma.

## 2. Functionals that are not of Crofton type do violate the inequalities

`F_max(A) = min_c max(ℓ_L(c), ℓ_R(c))` violates monogamy of mutual information in 923/3000 random cylinder
configurations and in 299/600 four-holed-sphere configurations, and the 5-party cyclic inequality in 559/3000
(EXPERIMENTS.md). So "any reasonable combination of the two chiral lengths with a common class" is not
enough; linearity in the measure is used.

## 3. Regularisation artefact

Comparing finite parts (cut-off term dropped) for an inequality that is not balanced at some endpoint gives
meaningless signs: strong subadditivity "failed" in 1830/3000 single-metric cylinder trials. With the cut-off
kept there are none. Theorem 1 is accordingly stated with a cut-off and an `o(1)`.

## 4. Invalid test set-ups caught (not counterexamples)

* Random Schottky groups for the four-holed sphere: essential curves entered a minimal cut in ≤ 7% of trials,
  so "zero violations" was nearly vacuous.
* Wrong-side Fenchel–Nielsen gluing: non-discrete representation (trace 1.195 element); the single-metric
  control then showed 325/500 spurious violations.

## 5. Failure-mode analysis for the two halves of the HRT identification

The former Hypothesis H had two halves. For each: what could have gone wrong, what was tried, what was found.

**H1 — "every admissible HRT candidate is an embedded cut class."**

| Conceivable failure | Outcome |
|---|---|
| Candidate geodesics cross when drawn on the surface | Excluded: they are disjoint curves on an embedded disc, so their endpoints on the boundary circle are unlinked (Jordan curve theorem); the left coordinate preserves cyclic order. |
| Candidate uses the *inner* closed geodesic (length `½|ℓ_L − ℓ_R|`, shorter than the outer one) | Excluded: it is not in the spacetime (PROOFS §6.2(a)). If it were admitted, the minimum would change and the inequalities would not follow from this proof. |
| Topological homology satisfied but no spacelike region (the gap between "homologous" and "spacelike-homologous") | This is H2, below. For H1 only the existence of *some* slice is used. |
| Marking mismatch: spacetime classes identified with surface classes up to a Dehn twist | Excluded by construction of the marking through the equivariant boundary map (§6.1); a mismatch would be visible as violations in the twisted `k = 2` rows of E6 only if it changed the minimum, so this rests on the proof, not on numerics. |
| Self-intersecting or multiply covered closed geodesics as candidates | Do not lie on an embedded slice; and Lemma 6.8 shows non-embedded collections never beat the minimum. |

**H2 — "every embedded cut class bounds a region on an admissible slice."**

| Conceivable failure | Outcome |
|---|---|
| Two geodesics of the class (or two group translates) are timelike or null related | Excluded by the Separation Lemma (points on hulls of complementary arcs of an acausal circle are spacelike related; equality only at shared ideal points). 240,000 random tests + 360,000 on a rotating pleated slice, none. |
| Pairwise spacelike curves that lie on no common acausal surface | Was the main worry (the user's "pairwise separation is not sufficient"). Resolved constructively: the surface is built (Proposition P) as a union of totally geodesic ideal triangles of an invariant triangulation containing the cut; the lemma applies to *all* pairs of points of the surface, not only to points on the curves. |
| The surface is not embedded / not a disc / has the wrong boundary | Proved by an explicit homeomorphism from the Klein disc (piecewise projective on triangles). **Not tested numerically.** |
| Null degeneration | **Occurs, in a limited sense.** At a closed geodesic with `ℓ_L ≠ ℓ_R` the faces of the pleated surface tilt to the null planes through the geodesic (`null_tilt_check.py`). Every pair of points is still strictly spacelike related, so the surface is acausal; it is not uniformly spacelike, and not smooth. |
| Slice is acausal but not a Cauchy slice of the spacetime-with-boundary | Not proved here. An argument is sketched in PROOFS §6.3 Remark (iv); it relies on the cited product structure of globally hyperbolic spacetimes with timelike boundary. |
| A closed geodesic of the cut cannot be an edge of a locally finite invariant ideal triangulation (its stabiliser would have to preserve the adjacent triangle) | **True, and raised as an objection.** The construction never uses such a triangulation: on the collapsed surface the closed geodesic is a vertex, and after straightening the triangles spiral onto it; no triangle is adjacent to it. Lemma T, Step 4, proves disjointness, local finiteness off the axes and covering directly. One topological fact (triangulations with prescribed vertices) is cited. |

**Net:** no counterexample to H1 or H2 was found. The one real degeneration (null tilt at rotating closed
geodesics) limits the theorem to definitions of HRT that use acausal slices, which is the convention of
Headrick–Hubeny–Lawrence–Rangamani (Cauchy slices "required to be acausal") and is contained in Wall's
("complete achronal slice", null permitted).

## 6. First-phase results (quantum information; archived in archive_qi/)

* Ordinary-Petz recovery bound for conditional mutual information: violated. Independently rediscovered here
  (ratio `I/(−log F)` down to ≈ 0.69 at second order for qubits, verified with an independent 80-digit
  implementation), but **already known**: certified counterexamples were recorded in September 2026, and the
  general-channel form was disproved in arXiv:2309.02074. No novelty claimed.
* Zhang–Yeung inequality for von Neumann entropies: no violation found by three different searches. Open.
* Transpose-degradable versus degradable channels: no separation found by SDP in dimensions 2×2, 2×4, 4×2,
  3×3. Suggests, but does not prove, that transpose-degradable implies degradable.
