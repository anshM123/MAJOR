# Covariant holographic entropy inequalities for vacuum rotating wormholes in 2+1 dimensions

**Ansh Mishra and Aryan Senthilkumar**

*Draft of 15 November 2026.*

*Note on methods and status.* The proofs and numerical tests were developed with AI assistance (Claude) and
checked by automated adversarial review and numerical falsification tests. They have not yet been checked by
independent human experts and are not formalised. The status of every statement is given in §10. Complete
proofs are in PROOFS.md; this text gives the statements and the ideas.

## Abstract

For time-symmetric states the Ryu–Takayanagi (RT) entropies obey an infinite family of linear inequalities,
the holographic entropy cone. Whether the covariant Hubeny–Rangamani–Takayanagi (HRT) entropies obey the same
inequalities is open; in 2+1 bulk dimensions it was known for simply connected bulks and for `π_1 = Z`, with
only numerical evidence for multiply connected bulks (arXiv:2407.07165). We settle the question for vacuum
solutions of arbitrary topology: rotating BTZ, rotating multi-boundary wormholes, wormholes with handles.
These are labelled by two hyperbolic structures `(X_L, X_R)` on one surface. We prove two things. (i) The HRT
entropy, defined by the usual prescription — the minimal extremal curve anchored on the endpoints of the
region and bounding, together with it, an acausal region of a bulk slice — equals the minimum over
topological cut classes `c` of `½(ℓ_{X_L}(c) + ℓ_{X_R}(c))`, the same class on both sides. The spacelike
homology region is constructed explicitly as a pleated surface, using a separation lemma for convex hulls of
arcs of an acausal circle in the boundary of `AdS_3`. (ii) That functional satisfies every holographic entropy
inequality admitting a contraction map, for arbitrary boundary regions, and lies in the full RT cone when the
parties are unions of whole boundaries. Consequently the HRT entropies of these spacetimes obey all
contraction-map inequalities, with no time-symmetry, no energy condition and no admissibility hypothesis. The
mechanism for (ii) is that of Czech and Dong; the new step is to run it on the quotient surface, where a
crossing inequality holds for every homotopy class of test curve and can be integrated against the Liouville
measures of both chiral surfaces.

## 1. Introduction

**The question.** When the bulk has a time-reflection-symmetric slice, the entropy `S(A)` of a boundary region
is the area of a minimal surface on that slice and satisfies strong subadditivity, monogamy of mutual
information, and all further inequalities provable by contraction maps (Bao et al., arXiv:1505.07839). In a
general spacetime `S(A)` is the area of the minimal *extremal* surface (HRT). Strong subadditivity and
monogamy survive (Wall, arXiv:1211.3494). For the higher inequalities the October 2026 review of Grimaldi and
Headrick states: "The validity of higher HEIs for the HRT formula remains an open problem."

**What was known in 2+1 dimensions.** Czech and Dong (arXiv:1905.03787) proved all contraction-map
inequalities for simply connected bulks, using a positive measure on kinematic space. Their extension to
multiply connected bulks via the universal cover is incorrect: "an infinite number of surfaces enter on both
sides of the inequality" (Grado-White, Grimaldi, Headrick, Hubeny, arXiv:2407.07165). Those authors gave an
argument for `π_1 = Z`, explained why it fails for free groups, and tested 1,877 inequalities numerically on
rotating three-boundary and one-boundary-torus vacuum wormholes without finding a violation; they conclude
"The validity of the RT inequalities for multiply-connected 2+1 dimensional bulk spacetimes remains an open
question." Zhao (arXiv:2602.04888) proves a covariant graph-model theorem in any dimension conditional on the
non-emptiness of "exposed regions" for every pair of HRT surfaces, an assumption not established for any
class of spacetimes.

**Results.** Let `M` be a locally `AdS_3` spacetime determined by two convex-cocompact Fuchsian
representations `ρ_L, ρ_R` of the fundamental group of an oriented surface `S` of finite type with boundary,
with the same marking. Let the parties be disjoint unions of intervals of an acausal boundary Cauchy slice
and of whole asymptotic boundaries.

> **Theorem H.** The HRT entropy of a region `A` (Definition in §7) equals
> `(1/4G_N) · min_c ½(ℓ^ε_{X_L}(c) + ℓ^ε_{X_R}(c))`, the minimum over topological cut classes for `A`; the
> infimum in the definition is attained.

> **Theorem 2.** The HRT entropies satisfy every holographic entropy inequality admitting a contraction map:
> exactly, at every sufficiently small cut-off, if the inequality is balanced at each endpoint of the
> regions; and with a difference of the two sides diverging to `+∞` otherwise.

> **Theorem 3'.** If the parties are unions of whole asymptotic boundaries, the vector of HRT entropies lies in
> the RT holographic entropy cone.

Theorems 2 and 3' follow from Theorem H and from a statement with no physics in it:

> **Theorem 1.** Let `X_1,…,X_k` be hyperbolic structures with funnel ends on one marked surface and
> `w_j > 0`. The functional `A ↦ min_c Σ_j w_j ℓ^ε_{X_j}(c)`, the same class in every structure, satisfies every
> inequality admitting a contraction map.

**Why the common class is the difficulty.** If the two chiral lengths were minimised independently the entropy
vector would be a sum of two RT vectors and there would be nothing to prove. The spacetime geodesic is
labelled by one group element acting on both chiral factors, so the minimisation is over a common class.
Cutting and pasting the `X_L`- and `X_R`-geodesics separately does not help: the two outputs need not be
homotopic (§4.3).

## 2. Setting

`X = H²/Γ` is a complete hyperbolic surface all of whose ends are funnels. `∂_∞X` is its ideal boundary, a
union of circles; `P ⊂ ∂_∞X` is a finite set of marked points cutting it into elementary pieces, each assigned
to a party or to the purifier. A *cut* for a region `A` is an embedded multi-curve — closed curves and arcs
ending at the marked points in `∂A` — bounding a region whose closure meets `∂_∞X` exactly in `A`; `𝒞(A)` is the
set of cuts up to isotopy and deletion of trivial or repeated components. A *test class* is the homotopy
class of a closed curve, or of a proper arc with ends in `∂_∞X \ P` free to move within their pieces;
`i(τ, c)` is the geometric intersection number.

## 3. The topological contraction lemma

**Lemma A.** Let `Σ_l α_l S(A_l) ≥ Σ_r β_r S(B_r)` admit a contraction map `f`, and let `γ_l ∈ 𝒞(A_l)` be
arbitrary. There exist `δ_r ∈ 𝒞(B_r)` such that, for every test class `τ`,
`Σ_r β_r i(τ, δ_r) ≤ Σ_l α_l i(τ, γ_l)`.

*Idea.* Realise the `γ_l` as geodesics; label each cell of the arrangement by the bit string of its
memberships in the homology regions; let `U_r` be the union of cells on which the `r`-th coordinate of `f` is
1; resolve `∂U_r` at crossings and cap pairs of rays at marked points to obtain a cut `δ_r`. A geodesic
representative of `τ` meets each geodesic of the arrangement minimally; at each crossing point the adjacent
bit strings differ exactly in the coordinates of the `γ_l` containing that geodesic, and the contraction
property bounds the weighted number of `r` for which the point lies on `∂U_r`. ∎

The lemma is topological; this is what allows one set of `δ_r` to serve several metrics. The mechanism
(auxiliary geometry, crossing counts, positive measure) is that of Czech and Dong.

## 4. Length as an integral of intersection numbers

**Lemma B (Crofton with cut-offs).** Let `m̄` be the Liouville measure on complete geodesics of `X`. Almost
every geodesic runs from funnel to funnel and defines a test class. Regulate arcs by removing horoballs
`B_p(ε)` at marked points; below a threshold `ε_0`, `ℓ^ε(c) = ℓ^{ε_0}(c) + e(c) log(ε_0/ε)` exactly, where `e(c)`
is the number of arc-ends. With `𝒢_ε` the geodesics disjoint from all horoballs,
`ℓ^ε(c) = ∫_{𝒢_ε} i(τ(ĝ), c) dm̄ + κ e(c) + o(1)`, `κ = 1 − log 2`.

**Proof of Theorem 1.** Take minimisers `γ_l`, apply Lemma A, integrate over `𝒢_ε(X_j)` for each `j`, and use
Lemma B: `Σ_l α_l F^ε(A_l) − Σ_r β_r F^ε(B_r) ≥ κ(Σ_j w_j)Δ − o(1)` with
`Δ = Σ_l α_l|∂A_l| − Σ_r β_r|∂B_r| ≥ 0`. The left side equals `C + (Σ_j w_j)Δ log(ε_0/ε)`; for `Δ = 0` this forces
`C ≥ 0`. ∎

**4.3 The naive argument fails.** Two minimal-position arrangements of the same classes differ by triangle
moves, and a triangle move can change how `∂U_r` pairs the six boundary points of a small disc. So outputs of
cut-and-paste performed in `X_L` and in `X_R` need not be homotopic. Theorem 1 uses only one arrangement.

## 5. The full cone for whole-boundary parties

**Lemma G (graph realisation).** For a finite positive combination `μ = Σ w_k τ_k` of test classes,
`A ↦ min_c i(μ, c)` is the min-cut function of a finite weighted graph (vertices: intersection points of the
geodesic representatives and their ideal ends, the latter labelled by parties; edges: the segments).

**Theorem 3.** For parties that are unions of whole boundary circles, `(min_c Σ_j w_j ℓ_{X_j}(c))_A` lies in the
holographic entropy cone. (Push the Liouville measure forward to a countable combination of test classes,
truncate, apply Lemma G, and pass to the limit using a finite filling family of arcs to confine the
minimisers.) For regions with endpoints the corresponding statement is not proved.

## 6. Vacuum 2+1 gravity

`AdS_3 = PSL(2,R)` with `ds² = −det dx`. The pair `(ρ_L, ρ_R)` defines `Γ ⊂ PSL(2,R)²`, acting on the
invisible domain `D̃` of its limit set; `M = D̃/Γ` has one asymptotic region per end of `S` (Mess; Barbot;
Bonsante–Krasnov–Schlenker). A boundary Cauchy slice, with the limit set, is an acausal circle `C` in
`∂AdS_3 = RP¹ × RP¹`, the graph of an equivariant orientation-preserving homeomorphism; it identifies the
marked ideal boundaries of `X_L` and `X_R`, hence their cut classes, including boundary twists.

**Lengths.** `(ρ_L(γ), ρ_R(γ))` translates a spacelike geodesic of `D̃` by `½(ℓ_L(γ) + ℓ_R(γ))`, and this is
the only invariant spacelike geodesic of `D̃` translated non-trivially (the other invariant geodesic, the
inner horizon, is lightlike-related to the limit set). For boundary points `p = ab^T`, `q = cd^T`, the
regulated geodesic length is `log N + 2 log(1/ε)`, `N = ω(a,c)ω(b,d)`: half the sum of the two horocycle
distances. Hence the geodesics of class `c` have regulated length `½(ℓ^ε_{X_L}(c) + ℓ^ε_{X_R}(c))`.

## 7. The HRT entropy is the minimum over cut classes

**Definition.** A *bulk slice* anchored on an invariant acausal circle `C' ⊃` limit set is a `Γ`-invariant set
`Σ̃ ⊂ D̃` such that `Σ̃ ∪ C'` is an embedded closed disc with boundary `C'` and any two distinct points of `Σ̃` are
spacelike separated. An *HRT candidate* for `A` consists of a boundary slice `C'` through the endpoints of
`A`; a bulk slice anchored on it; a finite disjoint union `γ` of complete geodesics of `M` lying in the
quotient slice, each closed or anchored at endpoints of `A`, one end at each endpoint; and a region `R` of the
quotient slice with frontier `γ` meeting the boundary in the region corresponding to `A`. The HRT entropy is
`(1/4G_N)` times the infimum of the regulated lengths of candidates. This is the minimal-extremal-curve
prescription with the spacelike homology condition, for acausal slices.

**Separation Lemma.** Let `C'` be an acausal circle cut at two points `r, r'` into closed arcs `α, β`. Every
point of the convex hull of `α` is spacelike separated from every point of the convex hull of `β`, except
that the two hulls share the geodesic from `r` to `r'`.

*Proof.* Move `r, r'` to `(∞,∞)` and `(0,0)` by Möbius maps in the two factors. Then `C'` is the graph of an
increasing homeomorphism fixing `0` and `∞`; positive combinations of points `[[uφ(u), u],[φ(u), 1]]` of `α`
(`u, φ(u) ≥ 0`) are matrices `X` with non-negative entries, and those of `β` are matrices
`Y = [[y_{11}, −y_{12}], [−y_{21}, y_{22}]]` with `y_{ab} ≥ 0`. So
`tr(adj(X)Y) = x_{22}y_{11} + x_{11}y_{22} + x_{12}y_{21} + x_{21}y_{12} ≥ 2√(x_{11}x_{22}y_{11}y_{22}) ≥ 2√(det X det Y)`,
which says that `X/√det X` and `Y/√det Y` are not timelike or lightlike related; equality forces both onto
the geodesic `rr'` and equal. ∎

**Pleated slices.** Let `c` be a cut class, with lifted geodesics `𝓕` in `H²` (chords for the arcs, axes of
hyperbolic elements for the closed components).

*An obstruction.* An axis `A` with stabiliser `⟨h⟩` cannot be an edge of a triangle in any invariant family of
ideal triangles with disjoint interiors: `h` would have to preserve the triangle adjacent to `A` on that
side, and `h` fixes only the two endpoints of `A`. So a cut containing a closed geodesic does **not** extend
to a locally finite invariant ideal triangulation, and none is used.

*Lemma T (spun family).* There is an invariant family `𝒯` of ideal triangles with disjoint interiors and
vertices at lifted boundary points and endpoints of the axes, such that: every arc-chord of `𝓕` is an edge;
no axis of `𝓕` meets any triangle; the triangles together with the axes fill `H²`; the family is locally
finite away from the axes, and triangles accumulate only on the axes. *Construction.* Over each elementary
boundary arc use a dyadic family of triangles. Cut the rest along `𝓕` and the boundary chords into pieces
`Y`; collapse each closed geodesic in `∂Y` to a point and take a topological triangulation of the collapsed
surface with vertices at the marked points and the collapsed points (the one cited fact). In the universal
cover `Ỹ` every lifted edge joins two "corners" (ideal vertices or boundary axes). Choose for every boundary
axis one of its two endpoints, equivariantly — possible because the stabiliser of an axis fixes both — and
replace each lifted edge by the chord between the chosen points, each lifted triangle by the ideal triangle
on them. The choice map preserves cyclic order, so the chords do not cross and the triangles have disjoint
interiors. An end of a chord at an axis endpoint is asymptotic to the axis, which gives local finiteness away
from the axes by proper discontinuity; and a point of `Ỹ` off the axes that lay in no triangle would be
separated from a fixed triangle by infinitely many distinct chords crossing one compact segment disjoint
from the axes, contradicting that local finiteness. In the quotient, the triangles spiral around each closed
geodesic of `c` from both sides.

*Proposition P (the slice).* Replace each triangle of `𝒯` by the totally geodesic ideal triangle of `AdS_3` on
the corresponding points of `C`, and each axis by its spacetime geodesic; let `Z` be the union. Any two
pieces have vertex sets in complementary closed arcs of `C`, so by the Separation Lemma any two points of `Z`
are spacelike separated, and two pieces meet only along a common edge. The piecewise-affine map from the
Klein disc (in a projective chart, affine on each triangle and on each axis) is continuous — at a point of an
axis because the triangles near it have all their vertices near the two endpoints of the axis — and
injective, so `Z ∪ C` is an embedded closed disc. `Z` is invariant because `𝒯` is, and it contains the
spacetime geodesics of `c`. The triangles and axes on one side of the cut form the homology region.

**Proof of Theorem H.** `≥`: on any bulk slice, the lifted components of a candidate are disjoint arcs of a
closed disc with endpoints on its boundary circle, hence have unlinked endpoints (Jordan curve theorem); the
left coordinate preserves cyclic order, so the corresponding geodesics of `X_L` form an embedded multi-curve,
and the parity of the number of separating chords reconstructs an invariant homology region. `≤`: the
pleated slice. ∎

**A feature of the construction.** At a closed component with `ℓ_L ≠ ℓ_R` (a geodesic carrying angular
momentum) the faces of the pleated slice accumulate and tilt towards the two lightlike planes through the
geodesic: the slice is acausal, each face is spacelike, but it is not uniformly spacelike there. This is
forced for any invariant slice that is totally geodesic on one side of the geodesic, because the group acts
on its normal plane by a boost. If the homology region were required to be smooth and uniformly spacelike, a
further argument (replacing the slice near such geodesics by a boost-twisted smooth strip) would be needed;
it is not given here.

## 8. Numerical tests

These test conclusions, not proofs. Details and commands: EXPERIMENTS.md.

| Test | Setting | Result |
|---|---|---|
| Separation Lemma | random acausal circles; positive combinations over complementary arcs; points of non-crossing ideal triangles | 0 violations in 160,000 cases; control with crossing triangles finds timelike pairs in 16% of configurations |
| Lightlike tilt at rotating closed geodesics | Gram determinant on the planes of spinning triangles | constant for `ℓ_L = ℓ_R`; decays to 0 for `ℓ_L ≠ ℓ_R` |
| Theorem 1, cylinder, interval regions, 1–3 metrics | strong subadditivity, monogamy, 5- and 7-party cyclic | 0 violations in 3,000 trials per case (750 for 7-party); control functional violates in 19–31% |
| Theorem 1, four-boundary wormhole, whole boundaries, Schottky groups, exhaustive enumeration | monogamy, strong subadditivity | 1 metric: 0 / 40; 2 metrics with relative Dehn twists: 0 / 40 |
| Same topology, Fenchel–Nielsen construction, 1–3 metrics | monogamy, strong subadditivity | 0 / 600 per case |
| Lemma A itself, cylinder, exact integers | four inequalities | dominating output exists in 7,800 / 7,800; for a false inequality it fails in 1,077 / 2,400 |
| Lemma G on the disc (referee) | min over classes versus graph max-flow | equal in 5,464 / 5,464 |

Limitations: double precision; some sample filtering; Lemma A is tested directly only for `π_1 = Z`; the
embedding part of the pleated-slice construction is not tested numerically; three earlier test versions were
vacuous or used a non-discrete representation and are documented, not counted.

## 9. What is new and what is not

Not new: contraction maps; cut-and-paste in an auxiliary geometry with crossing counts and a positive measure
(Czech–Dong); the left/right length formula (Maxfield; Grado-White et al.); pleated surfaces in `AdS_3`
(Mess); the numerical evidence for the conclusion.

New, to the best of our knowledge after the search described in LITERATURE.md: the crossing inequality on the
quotient surface (Lemma A) and the Crofton formula with horoball cut-offs (Lemma B) in this context; the
multi-metric theorem; the graph realisation and the full-cone statement for whole-boundary parties; the
Separation Lemma; the construction of the spacelike homology region for every topological cut class; and the
resulting theorems that the HRT entropy of vacuum multi-black-hole spacetimes is the minimum over cut classes
and obeys all contraction-map inequalities in every topology.

Scope: vacuum only. Matter, higher dimensions, extremal horizons, conical defects and the black-hole
information problem are not addressed.

## 10. Status of statements

| Statement | Status |
|---|---|
| Lemma A, Lemma B, Theorem 1 | proved; three automated referee passes, no error found |
| Lemma G, Theorem 3 (whole-boundary parties) | proved; referee: correct after repairs, which are applied |
| Theorem 3 for regions with endpoints | **not proved** |
| Length formula, marking | proved |
| Separation Lemma | proved; tested numerically |
| Invariant spun family of ideal triangles containing a cut (Lemma T) | proved given one cited topological fact (triangulations with prescribed vertices); spiralling step proved |
| Pleated slice; Theorem H | proved given the above; independent review: see OPEN_GAPS.md |
| Theorems 2 and 3' (HRT) | proved; no admissibility hypothesis; for the definition of §7 with acausal slices |
| Same statements for smooth, uniformly spacelike homology regions | **not proved** |
| Non-vacuum spacetimes, higher dimensions | not covered |

## References

Bao, Nezami, Ooguri, Stoica, Sully, Walter, arXiv:1505.07839. Wall, arXiv:1211.3494. Czech, Dong,
arXiv:1905.03787. Grado-White, Grimaldi, Headrick, Hubeny, arXiv:2407.07165 and arXiv:2502.09894. Zhao,
arXiv:2602.04888. Grimaldi, Headrick, arXiv:2610.33905. Maxfield, arXiv:1412.0687. Headrick, Takayanagi,
arXiv:0704.3719. Hubeny, Rangamani, Takayanagi, arXiv:0705.0016. Mess, *Lorentz spacetimes of constant
curvature* (1990). Barbot, *Causal properties of AdS-isometry groups I, II*. Bonsante, Krasnov, Schlenker,
arXiv:math/0610429. Bonsante, Seppi, arXiv:2004.14414. Aminneborg, Bengtsson, Holst, *Class. Quantum Grav.* 15
(1998). Santaló, *Integral Geometry and Geometric Probability*. Penner, *Comm. Math. Phys.* 113 (1987).
Epstein, *Acta Math.* 115 (1966). Farb, Margalit, *A Primer on Mapping Class Groups*. Hatcher, *Topology
Appl.* 40 (1991). Fomin, Shapiro, Thurston, *Acta Math.* 201 (2008). Bonahon, *Ann. Fac. Sci. Toulouse* 5
(1996). Bao, Furuya, Naskar, arXiv:2506.18086.
