# Exploratory note: the generalised-measure theorem and four-dimensional black branes

Status: exploratory. Nothing here is refereed. The comparison with Bousso–Kaya (2025),
Akers–Hernández-Cuenca–Rath (2021) and Bao–Furuya–Naskar (2025) is delegated to the literature audit
(`ads4/PROBLEM_CANDIDATES.md`) and is **not** asserted here from memory.

## 1. The abstract statement (what the AdS3 proof actually used)

Let `X` be a surface as in PROOFS.md §1 with universal cover `H²`, deck group `Γ`, ideal boundary circle `S¹`,
and let `ν` be a `Γ`-invariant positive Radon measure on the space of unordered pairs of distinct points of
`S¹` (a geodesic current, allowed to charge chords ending in the domain of discontinuity) such that the
regulated intersection `i_ε(ν, c)` with every cut class is finite. Put `F_ν(A) := min_{c ∈ 𝒞(A)} i_ε(ν, c)`.

**Theorem 1′ (sketch of proof only).** `F_ν` satisfies every contraction-map inequality.

Reason: Lemma A gives, for arbitrary input classes, output classes with
`Σ_r β_r i(τ, δ_r) ≤ Σ_l α_l i(τ, γ_l)` for **every** test class `τ`; integrate over `τ` against `ν`. Theorem 1 of
PROOFS.md is the case `ν = Σ_j w_j L_{X_j}` (Liouville currents). What Theorem 1 used beyond this: the Crofton
formula (to identify `i(L_X, c)` with hyperbolic length) and null sets of non-generic chords; for a general
`ν` with atoms the generic-position step must be redone with intersection numbers of classes. Not written.

So the mathematical generalisation is essentially free. The question is whether it says anything new about
gravity.

## 2. When does an HRT functional have such a representation?

Let `D(a,b)` be the regulated HRT length (or area per unit transverse volume) for the region between boundary
points `a, b` of the universal cover. A positive measure with `D(a,b) = ν{chords separating a from b}` exists
iff the mixed second difference is non-negative,

`D(a,c) + D(b,d) − D(a,d) − D(b,c) ≥ 0` for cyclically ordered `a,b,c,d`,

i.e. iff strong subadditivity holds for adjacent intervals on the cover; then `ν = −∂_a∂_b D` (the density of
conditional mutual information: the "kinematic space" measure of Czech–Lamprou–McCandlish–Sully). For HRT,
strong subadditivity is Wall's maximin theorem, which assumes the null curvature condition. So:

> the positive-measure representation of single-interval HRT lengths is *equivalent* to strong
> subadditivity on the universal cover, and is therefore available exactly when the usual energy-condition
> hypotheses hold. It is not an additional structure special to vacuum.

Two further inputs are needed to turn `F_ν` into the HRT entropy of a general region, and both are
gravitational statements, not consequences of the measure:
(a) *composite regions:* `S(A) = min` over non-crossing pairings of sums of single-interval values, which
needs the minimal extremal curves of nested or disjoint intervals to lie on one Cauchy slice (the H2-type
statement; in vacuum AdS3 this is PROOFS.md Proposition P; in general it is maximin nesting);
(b) *closed components:* the length of the closed extremal curve in a class `g` must equal the `ν`-translation
length `lim_n [D(a, g^{n+1}b) − D(a, g^n b)]`.

## 3. Planar-symmetric time-dependent AdS4 black branes

Metric `ds² = g_{ab}(x^c) dx^a dx^b + e^{2σ(x^c)} dy²`, `x^a = (t, z, x)` (e.g. Vaidya-AdS4). For strips
`{x ∈ I} × R_y`, `y`-invariant extremal surfaces are geodesics of the reduced 2+1 metric `ĝ = e^{2σ} g`. The
reduced manifold has fundamental group trivial (`x ∈ R`) or `Z` (`x` periodic), one- or two-sided.

* `x ∈ R`: single strips give `D` with the crossing inequality (maximin in the 4D spacetime); composite
  regions are minima over pairings. This is the hypothesis set of the simply connected theorem of
  Czech–Dong (arXiv:1905.03787), whose kinematic-space argument uses only these two properties. **Expected
  verdict: a corollary, not a new theorem** (to be confirmed by the audit: whether their hypotheses are
  stated so as to include a reduced, non-Einstein 2+1 metric).
* `x` periodic: the universal cover is the planar spacetime; winding arcs are single-strip values on the cover;
  closed curves are `x`-invariant extremal cross-sections. Statement (b) is then the physical claim
  "entropy density of the whole circle = slope of the single-strip entropy at infinite width". It holds at
  `t = 0` of the eternal brane. In a two-sided brane perturbed by a shock it is consistent with the
  Liu–Suh picture (at fixed time the infinitely wide strip runs along the old horizon, giving the old
  entropy density, which is the fine-grained entropy density of one side), but **it is not proved here**.
  The `π_1 = Z` case of 2+1 HRT is treated in arXiv:2407.07165; whether their hypotheses cover the reduced
  metric is again for the audit.

**Conclusion for the question asked.** For AdS4 black branes with planar symmetry the generalised-measure
theorem would at most re-derive, in a symmetric sector, inequalities already expected from existing 2+1
theorems. Its genuinely new content is in 2+1 dimensions (matter, non-abelian `π_1`), conditional on (a) and
(b). It is **not** a result about four-dimensional gravitational dynamics and should not be the primary
target on that basis.

## 4. The obstruction in four dimensions, and what it points to

For regions without a transverse symmetry the HRT surface is a genuine 2-surface in a 3+1 spacetime and
there is no reduction to geodesics. The AdS3 proof worked because of two facts with no 4D analogue:

1. **A positive Crofton representation.** On a hyperbolic 3-space slice, area is the invariant integral of
   the number of intersections with lines. For a spacelike 2-surface in `AdS_4` not contained in a
   time-symmetric slice, the natural candidates (totally geodesic 2-planes of signature (1,1) meeting the
   surface in points) have non-compact isotropy — the boosts fixing a point of the surface and its normal
   plane — so the invariant measure of the planes meeting any piece of surface is infinite. In `AdS_3` the same
   divergence is avoided only because the isometry group is a product, `SO(2,2) ≃ SL(2,R) × SL(2,R)` up to
   finite index, and lengths split into two Riemannian (hyperbolic) Crofton integrals. `SO(2,3)` is simple.
2. **A purely topological description of candidates** (cut classes on a surface) on which the contraction
   lemma acts independently of the metric.

This isolates a sharply posed question that is genuinely four-dimensional and, as far as this note knows,
open (status to be confirmed by the audit):

> **Question V4.** In the vacuum state of a holographic CFT3 (bulk = pure `AdS_4`), let the regions lie on an
> arbitrary, non-time-symmetric boundary Cauchy slice. Do the HRT entropies satisfy all RT inequalities
> (beyond strong subadditivity and monogamy)?

Here there is no matter, no topology and no horizon: everything is fixed by the conformal geometry of the
boundary slice. A positive answer would need a Lorentzian replacement for the Crofton formula; a negative
one would be a new constraint distinguishing time-dependent from static holography in `d ≥ 3`. The 2+1
analogue is true (simply connected case). No proof strategy is available yet; recorded as a candidate, to be
ranked with the others.

## 5. Failure of one construction versus a no-go theorem

Section 4 shows that **one** construction fails: the isometry-invariant measure on totally geodesic
(1,1)-planes of `AdS_4` is infinite on the planes meeting any piece of spacelike surface, so there is no
invariant Lorentzian Crofton formula of that type. This is **not** a no-go theorem for positive-measure
representations in general, for three reasons.

1. *Non-invariant measures are not excluded.* In the reduction of §2 the measure `ν = −∂_a∂_b D` is built from
   the entropies themselves and is not invariant under anything. An analogue for 2-surfaces would be a
   positive measure on some auxiliary space of "test objects" whose intersection numbers with HRT surfaces
   reproduce areas; nothing here rules that out.
2. *For finitely many regions, existence of a representation is equivalent to the conjecture itself.* By the
   graph-model characterisation of the RT cone (Bao et al., arXiv:1505.07839), an entropy vector lies in the
   cone iff it is the vector of minimal cuts of a weighted graph, i.e. iff it has a positive combinatorial
   "measure" representation. So asking whether *some* positive representation exists for each finite family of
   regions is exactly Question V4, not a stronger or weaker statement.
3. *What a genuine no-go would be.* Either an explicit configuration violating one RT inequality (which
   would settle V4 negatively), or a proof that no assignment of test objects can be made consistently for
   all regions at once (a statement about a continuum representation, strictly stronger than failure of the
   inequalities and not implied by §4).

What §4 does establish is narrower and should be quoted as such: the specific mechanism of the AdS3 proof —
splitting the area functional into Riemannian Crofton integrals via a product structure of the isometry
group — is unavailable in `AdS_4`.
