# Finite orbifold quotients preserve the vacuum full-cone theorem

Research checkpoint, 28 November 2026. This result uses the standing vacuum
Theorem H and the equivariant form of the arrangement construction in
`PROOFS.md` and `fullcone/FULLCONE.md`, at repository revision
`acda6f59a240537f521583630ff9e9fac71e145f`. It adds no maximin assumption.

## 1. Precise statement

Let \(\widetilde M\) be a smooth vacuum spacetime in the class covered by
Theorem H and the full-cone arrangement theorem. Let a finite group \(F\),
of order \(q\), act by orientation- and time-orientation-preserving
isometries. Assume the action respects the conformal completion and its
double: each group element extends smoothly and conformally to the double,
preserves its time orientation, and commutes with boundary reflection.
Assume also that its induced actions on the two marked chiral hyperbolic
surfaces are isometries in the common topological marking. Take the quotient
\(M=\widetilde M/F\), including its fixed timelike axes as orbifold particle
worldlines. Require no fixed points on the conformal boundary.

Let \(A_1,\ldots,A_N\) be disjoint boundary regions on a boundary Cauchy slice
of \(M\), and let \(\widetilde A_i\) be their full preimages. Use compatible
quotient/lifted endpoint cutoffs. A quotient HRT candidate is an embedded
regular-locus spacelike extremal multicurve with a homology region on an
acausal orbifold Cauchy slice; cone points are filled for that homology.

**Theorem.** For every party union \(A\),

\[
 S_M(A)=q^{-1}S_{\widetilde M}(\widetilde A).
 \tag{1}
\]

There is one finite positive-capacity graph representing all the quotient
entropies simultaneously. Thus the quotient vector lies in the full static
RT cone, with precisely the endpoint/UV convention of the cover theorem.
In particular, regulated vectors are positive-graph vectors for all
sufficiently small compatible cutoffs, for which the cover's endpoint
Bell capacities are nonnegative.
The conclusion holds for any finite number of parties and arbitrary boundary
Cauchy time profiles. The quotient HRT minimum is attained on an acausal
Cauchy slice that is smooth and uniformly spacelike in orbifold charts of
the conformal completion; the minimizing curves avoid every particle
worldline. The uniform margin depends on the geometry and the cut, as in
the cover theorem.

This is a theorem for finite orbifold quotients, not for arbitrary cone
angles, spinning particles, or collisions.

## 2. An equivariant arrangement graph exists

All endpoint lifts form a finite \(F\)-invariant set. For each required
boundary union choose a minimizing cut class in the cover. In the
arrangement construction include that class and all of its \(F\)-translates.
The resulting finite arrangement in the reference chiral hyperbolic surface
is \(F\)-invariant. Adding these translates does not invalidate the graph
proof: every selected optimum remains an arrangement cut, and the lower
bound applies to every cut of the refined arrangement.

The total renormalized intersection functional is \(F\)-invariant, since
the group acts by isometries of each chiral metric. Choose the auxiliary
test-geodesic measure equivariantly. Equivalently, replace any realizing
measure by its finite average \(q^{-1}\sum_{g\in F}g_*\nu\). Its intersection
pairing with any class is unchanged. Endpoint terms can likewise be chosen
constant on endpoint orbits, because the cutoff is lifted from the quotient.

Consequently each arrangement-edge capacity equals the capacity of each of
its translates. Add the endpoint Bell edges equivariantly. The resulting
finite graph \(\mathcal G\) carries an \(F\)-action by graph automorphisms;
all boundary labels are \(F\)-invariant. The complete cover entropy vector
for the lifted parties is its mincut vector.

## 3. A minimum cut can be chosen invariant

For a set \(U\) of graph vertices, let \(c(U)\) be its undirected cut
capacity. Nonnegative edge capacities imply

\[
 c(U\cup V)+c(U\cap V)\le c(U)+c(V).
 \tag{2}
\]

If \(U,V\) are both minimum cuts with the same fixed terminal labels, their
union and intersection are admissible and each has capacity at least the
minimum \(m\). Equation (2) forces both to have capacity exactly \(m\).

Starting from any minimum cut \(U\), form

\[
 U_F=\bigcup_{g\in F}gU.
 \tag{3}
\]

Every \(gU\) is a minimum cut with the same labels. Finite iteration of (2)
shows that \(U_F\) is a minimum cut and is \(F\)-invariant. This argument
does not assume a unique minimum or strict inequalities.

## 4. From the invariant graph cut to an invariant topological cut

First ignore the added endpoint Bell edges, whose contribution is fixed by
the boundary union. The selected arrangement cells of \(U_F\) form an
invariant, possibly nonsmooth homology region. Smooth its frontier
equivariantly in small disks around arrangement vertices, using an
\(F\)-invariant auxiliary Riemannian metric and orbitwise equal radii.
A regular-neighborhood construction chooses the same resolution at all
vertices in an orbit and is automatically invariant at vertices with
nontrivial stabilizers. At a fixed vertex, filling a small invariant disk
and rounding the unselected sectors, or the complementary construction,
provides such a resolution. Endpoint half-disks are treated with their
marked endpoint fixed; the group acts freely on these endpoints.

The resolutions give an embedded cut class \(c\) with the required trace.
Its isotopy class is invariant. Inessential loops and redundant parallel
boundary pairs are removed in their full group orbits; this is the usual
reduction of a cut class and preserves the homology trace.

The standing cover construction uses a measure absolutely continuous with
respect to Liouville measure on proper arc geodesics: almost every test
geodesic has both endpoints outside the limit set and the marked endpoint
orbits, meets the finite arrangement transversely in finitely many points,
and hits no arrangement vertex. Finite averaging preserves these facts.
For each such proper geodesic, sufficiently small smoothing disks and
endpoint horoballs are disjoint from it and do not change its crossings
with the frontier. Merely avoiding a vertex would not suffice for an
arbitrary recurrent complete geodesic; the properness supplied by the
cover's test measure is used here. The intersection number of its test class with \(c\)
is therefore no larger than its crossings with the original selected
arrangement edges. The class is constant for all sufficiently small
resolutions. Integrating the pointwise inequality gives

\[
 \mathcal L(c)\le \operatorname{cap}(U_F)=S_{\widetilde M}(\widetilde A).
 \tag{4}
\]

Here \(\mathcal L\) includes the same endpoint terms and the factor
\(1/(4G_N)\) as the entropy functional. The cover's minimum-over-classes
theorem gives the reverse inequality.
Thus \(c\) is a minimizing cover class.

Every reduced simple class has its unique chiral geodesic representatives,
and Theorem H gives its unique spacetime geodesic multicurve \(K\) in the
vacuum domain. Isometries carry these representatives to the representatives
of the transformed classes. Since \(c\) is invariant, \(K\) is invariant.
Theorem H supplies acausality, embeddedness, and a homology Cauchy slice.

## 5. An invariant smooth Cauchy slice through the curve

Work first in the smooth conformal double of the cover. The finite group
action extends to that double and commutes with boundary reflection. It
preserves the causal cones; isometry of an arbitrarily selected conformal
representative is unnecessary. If an invariant representative is desired,
average a reflection-invariant representative as
\(\bar g_F=q^{-1}\sum_g g^*\bar g\). The summands are pointwise positive
conformal multiples, so the result is a smooth Lorentz metric with the
same cones, invariant under both the group and reflection. The
doubled compact curve is invariant and acausal. The cover's smooth-slice
theorem supplies a smooth spacelike Cauchy surface \(S\) containing it.

Bernal--Sanchez, arXiv:gr-qc/0512095, Theorem 1.2, supplies a smooth Cauchy
temporal function \(T\) with \(S=T^{-1}(0)\). Define

\[
 T_F(x)=q^{-1}\sum_{g\in F}T(gx).
 \tag{5}
\]

Every summand has a past-directed timelike differential (equivalently the
appropriate directed timelike gradient). The timelike covector cone is
convex, so \(T_F\) is temporal. Along any future-inextendible causal curve,
each of the finitely many summands tends to \(+\infty\); along any
past-inextendible one, each tends to \(-\infty\). Thus \(T_F\) is Cauchy.
Since \(K\) is invariant and \(T\) vanishes on it, \(T_F\) vanishes on it.

The level surface \(S_F=T_F^{-1}(0)\) is an invariant smooth spacelike
Cauchy surface through \(K\). Restrict it to the original conformal half.
Transversality to the timelike conformal boundary and the Cauchy property
follow exactly as in Proposition S+ of the cover theorem.
The restricted surface is compact by that proposition's compactness
argument. Strict positivity on its compact auxiliary unit tangent bundle
therefore supplies a uniform spacelike margin; this descends in orbifold
charts because the covering group is finite.

Transport a homology region for \(K\) from the original slice to this one
along a complete future timelike flow tangent to the conformal boundary.
The transport fixes \(K\) and preserves its marked boundary trace. No
equivariance of that transport is needed for the next argument.

## 6. The homology region is invariant and the curve misses every axis

**Coloring lemma.** On a connected surface with nonempty boundary, a
specified embedded reduced cut multicurve and its boundary trace determine
the homology-region coloring uniquely. Indeed, the symmetric difference of
two such colorings has no jump across any cut component, so its indicator is
locally constant on the entire connected surface. It vanishes in a boundary
collar where the traces agree, hence vanishes everywhere.

Apply this to the region \(R\subset S_F\) just constructed. Both \(gR\) and
\(R\) have frontier \(K\) and the same invariant boundary trace, so the
coloring lemma implies \(gR=R\).

Suppose \(p\in K\) were fixed by a nonidentity deck element. On the invariant
spacelike plane \(T_pS_F\), its linear action is a nontrivial finite rotation.
Indeed it preserves the future normal and the spatial orientation; an
identity spatial derivative would be an identity ambient derivative, and
an isometry fixing a point with identity derivative is the identity on the
connected spacetime.
Embeddedness makes the tangent cone of \(R\) at \(p\) one of the two closed
half-planes bounded by \(T_pK\). No nontrivial rotation preserves a
half-plane. This contradicts invariance of \(R\).

This explicitly includes isotropy order two: a rotation by \(\pi\) preserves
the tangent line but exchanges its two sides. Thus the angle-\(\pi\) case
does not create a projected cut ending at a particle.

The quotient \(S_F/F\) is an acausal Cauchy orbifold surface, smooth in
orbifold charts and orthogonal to each timelike particle axis. Its homology
region \(R/F\) includes or excludes the filled cone points normally. The
projected multicurve \(K/F\) is regular, embedded, spacelike, extremal, and
admissible. Its length is \(q^{-1}\operatorname{Length}(K)\).
One can also see the quotient's causality directly: \(T_F\) descends to a
temporal function, and quotient causal curves lift locally across finite
orbifold charts. A lift of an inextendible quotient curve is inextendible,
since any endpoint upstairs would project to an endpoint downstairs. Thus
the descended temporal function is Cauchy.

The candidate convention is substantive: a cut is a codimension-one
boundary of a homology region, with no particle endpoint or additional
homology boundary. If one instead defines an orbifold-extremal cut by
requiring a smooth extremal lift, the same invariant-half-plane argument
forces it to avoid every axis, so the two conventions agree. No claim is
made here for a different prescription that treats the particle itself as
an allowed endpoint, excises a tube, or assigns an entropy term to that
internal boundary.

## 7. Equality with quotient HRT and the quotient graph

Any admissible regular-locus quotient candidate lifts to an admissible cover
candidate of \(q\) times its length: finite properness preserves the Cauchy
property on lifting, and homology regions lift. Therefore

\[
 S_{\widetilde M}(\widetilde A)\le qS_M(A).
\]

The projected minimizing curve of Sections 4--6 proves the opposite
inequality. This proves (1), with an attained minimum.

Finally form \(\mathcal Q=\mathcal G/F\). Its vertices are vertex orbits;
between distinct orbits assign the sum of capacities of all original edges
between them, divided by \(q\). Loops are irrelevant. Every cut of
\(\mathcal Q\) corresponds to an invariant cut of \(\mathcal G\), with
capacity divided by \(q\). Section 3 then implies

\[
 S_{\mathcal Q}(A)=q^{-1}S_{\mathcal G}(\widetilde A)=S_M(A).
\]

One graph works for every union. This proves the full RT-cone assertion.

## 8. Geometric scope and why it does not include collisions

A geometrically finite orientable hyperbolic orbifold with cone angles
\(2\pi/n_i\), integers \(n_i\ge2\), has a finite torsion-free normal cover:
apply Selberg's lemma to its finitely generated Fuchsian group and take the
normal core. For two compatible chiral orbifold holonomies, restriction to
that subgroup gives the smooth vacuum cover used above. The disk with a
single order-\(n\) point has the elementary finite cover by the disk and is
included separately. Retain the cover theorem's other global assumptions;
in particular, cusp/extremal degenerations are a separate extension.

These cone points have positive deficit and zero intrinsic spin. This is
unconditional matter in the finite-orbifold class, including arbitrary
allowed chiral moduli and boundary times. It is not restricted to a static
boundary time slice.

A finite time-orientation-preserving isotropy group in a Lorentz three-plane
is conjugate to a subgroup of spatial rotations; it fixes one timelike line.
Consequently the singular set locally consists of one timelike axis, not
distinct axes meeting at a collision vertex. Collisions are therefore not
covered by this finite-orbifold quotient theorem. Generic cone angles are
also not limits automatically covered by a dense orbifold-angle set:
\(\{2\pi/n:n\ge2\}\) is not dense in \((0,2\pi)\).

## Primary reference

A. N. Bernal and M. Sanchez, *Further results on the smoothability of Cauchy
hypersurfaces and Cauchy time functions*, Lett. Math. Phys. 77 (2006),
183--197, arXiv:gr-qc/0512095, Theorems 1.1--1.2.

The new deduction is finite-group descent of the repository's existing
vacuum arrangement and HRT theorems. Its validity is conditional on those
standing theorems being correct; it does not independently re-audit them.

## Exact combinatorial check

`check_equivariant_mincut.py` independently enumerates every cut in 321
small symmetric graphs, using integer cover capacities and exact rational
quotient capacities. Cyclic group orders are 2, 3, 4 and 6, with fixed
vertices and stabilized edges. It checks unrestricted minimum = invariant
minimum = group-order times quotient minimum. An explicit degenerate
order-two graph has four minima, two asymmetric; taking the union of their
translates recovers an invariant minimum. All checks pass; the compact
result is in `check_equivariant_mincut.json`. This validates the algebraic
normalization and degeneracy handling, not the geometric hypotheses.
