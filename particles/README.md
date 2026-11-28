# Point particles, cone cuts, and extremal limits

Research extension, 28 November 2026, based on repository revision
`acda6f59a240537f521583630ff9e9fac71e145f` and its standing vacuum Theorem H,
smooth-slice theorem, and exact full-cone arrangement construction.

## What this work establishes

The general collision theorem is not established here, and the proposed
arbitrary-angle chiral-cone HRT identification is false in an explicit
static example. There is an actual matter extension with no new maximin assumption, an exact
all-time positivity theorem for the explicit two-particle interval kernel,
and a quantitative treatment of the extremal limit.

| Result | Status and scope |
|---|---|
| Finite spinless orbifold particles | HRT/full-cone descent from a compatible smooth vacuum cover; arbitrary finite party number and allowed boundary times |
| Smooth homology slice through quotient HRT | Proved in orbifold charts, including avoidance of every particle axis |
| Filled-class cone Crofton | False at every positive-mass cone angle; an explicit positive-measure counterexample is supplied |
| Arbitrary-angle marked cone functional | Full RT-cone theorem for the relaxed common-class length infimum |
| Cone angles below pi | Exact attained geodesic-length functional, via uniform collar exclusion and funnel replacement |
| Naive arbitrary-angle chiral-cone HRT identity | Explicitly false on a tilted boundary slice of a static spinless cone |
| Ordinary common cusp ends | Exact positive-graph theorem for the geometric cut functional, with unmarked cusp colors |
| Compact extremal BTZ | Full-cone theorem for the nonextremal limiting entropy functional, with a uniform error bound |
| Unwrapped extremal BTZ | Actual bounded-interval HRT full-cone theorem using the explicit Poincare patch |
| Two-massless-particle formation kernel | Exact positive kinematic measure at all BTZ boundary times |
| Two-particle multiparty matching functional | Exact all-party RT-cone theorem and a 372-facet numerical audit |
| Physical HRT through collisions for every union | Still requires simultaneous admissibility/acausality; not inferred from endpoint noncrossing |

All HRT claims are at the classical area/length order. The endpoint cutoff
or positive endpoint Bell terms are retained for full-cone membership.
Endpoint-balanced inequalities also apply directly to the renormalized
finite parts. The standing vacuum theorems are inputs to the new deductions;
this addition is not an independent verification of the entire repository.

## The unconditional matter result

Let a finite orientation- and time-orientation-preserving group \(F\) act
isometrically on a smooth vacuum spacetime \(\widetilde M\) covered by the
standing theorem. Assume compatibility with the marked chiral surfaces and
the smooth conformal double, and no boundary fixed points. In the quotient
\(M=\widetilde M/F\), include the timelike fixed axes as ordinary filled
spinless orbifold particles. For every boundary region,

\[
S_M(A)=\frac{1}{|F|}\,S_{\widetilde M}(\widetilde A),
\]

where \(\widetilde A\) is the full preimage and cutoffs are compatible.
One finite positive-capacity graph represents every party union.

The proof does not assume the cover minimizer is unique or symmetric.
Submodularity of graph cuts provides an invariant minimum. An invariant
frontier class then gives an invariant geodesic cut. Averaging a Cauchy
temporal function constructs an invariant smooth slice through it. Its
homology-side coloring is invariant, so a nontrivial stabilizer cannot fix
a point of the curve: a rotation cannot preserve the tangent half-plane.
This includes the order-two case.

The resulting particle angles are \(2\pi/n_i\), \(n_i\ge2\). These are
genuine positive-deficit particles, but the angles are not dense in
\((0,2\pi)\). Finite Lorentz isotropy fixes a single timelike axis; this
theorem therefore does not include collision vertices.

Full proof: `FINITE_ORBIFOLD_THEOREM.md`.

## The cone issue is a topology issue, not a measure-zero exception

On a cone with any \(0<\theta<2\pi\), a globally minimizing anchored cut
can be crossed twice by an open family of regular embedded test geodesics,
while their geometric intersection number in the filled disk is zero.
Neither geodesic hits the tip. Thus removing cone-hitting geodesics does not
restore the original intersection formula, even below angle \(\pi\).

Particles should be filled for **homology**, while their marks remain in
the **isotopy classes**. The finite-cover theorem supplies one rigorous
way around the resulting analytic problem.

Proof and figure: `CONE_INTERSECTION_AUDIT.md`, `figures/cone_bigon.pdf`.

### A replacement proof for every cone angle

The Crofton failure can be bypassed. Replace a radius-\(a\) cone neighborhood
by a smooth negatively curved funnel, with its homology color free. The
replacement admits a marking-preserving 1-Lipschitz map back to the old
punctured cone. Consequently every marked class length is bounded below by
its old infimum. Any fixed regular comparison curve misses all sufficiently
small collars, giving the matching upper bound as \(a\to0\).

Choose one nearly optimal common class for each required region, and then
nearly optimal representatives in each metric. This proves convergence of
the common-class infimum without interchanging an uncontrolled infinite
minimum and a limit. The smooth replacements have positive graph models;
the endpoint constants are independent of \(a\). Closedness of the RT cone
therefore proves full-cone membership of the relaxed marked-cone
functional at every angle in \((0,2\pi)\), and mathematically beyond that
range as well.

Below angle \(\pi\), every embedded complete regular geodesic avoids a
uniform collar,
\[
r<\operatorname{arctanh}\!\left(\tanh R_{\rm collar}\cos(\theta/2)\right).
\]
Replacing the cones inside those collars preserves all reduced simple
geodesic lengths exactly. No relaxation or limit is then needed.

These are mathematical functional theorems. They become physical entropy
theorems wherever the HRT/common-class identification is independently
established. The proofs are in `ALL_ANGLE_RELAXED_FUNCTIONAL.md` and
`SMALL_ANGLE_FUNCTIONAL_THEOREM.md`; `COMMON_SIDE_ROUNDING.md` records a
further local improvement lemma and its precise limitation.

## Exact positivity in the formation example

For the explicit equal-time interval kernel of arXiv:1706.07390, write
\(X=e^{Ra}<1<Y=e^{Rb}\), \(Z=e^{Rt}\), and
\(E=\coth(\pi R/2)\). The crossing expression factors as

\[
D_{\mathrm c}=\log\frac{F_1F_2}{XYZ^2},\qquad
\begin{cases}
F_1=E(X-Z)(Y-Z)+Z(Y-X),\\
F_2=E(XZ-1)(YZ-1)+Z(Y-X).
\end{cases}
\]

Both factors are strictly positive on every shortest-arc chart. Their
bilinear determinants equal \(Z^2\). Moreover, if
\(h=D_{\mathrm c}-D_{\mathrm d}\), then \(h_a\ge0\) and \(h_b\le0\).
Concave smoothing of the branch minimum therefore proves

\[
\partial_a\partial_b\min(D_{\mathrm c},D_{\mathrm d})\ge0
\]

as a distribution, including singular switching loci. Injection and
antipodal seams also contribute no negative measure. The finite disk graph
theorem then places the noncrossing-matching functional in the full RT cone
at every time, for every party number.

This settles the algebraic/current part for the explicit collision kernel.
It does not silently identify every minimizing matching with a physical
HRT multicurve. `COLLISION_ADMISSIBILITY.md` records the further geometric
analysis and the precise remaining condition.

The numerical audit evaluates 2,169,504 facet values from 5,832 matching
vectors. All 372 five-party facets pass within \(10^{-10}\); the minimum
is \(-1.42\times10^{-13}\). The analytic theorem, rather than this sample,
establishes the all-time assertion for the matching functional.

## A physical counterexample to the assumed chiral identity

In static conical AdS with \(\theta=3\pi/2\), choose
\(\Delta\phi=11\pi/16\) and \(\Delta t=9\pi/16\). These endpoints fit on
a spacelike boundary Cauchy curve. The direct spacetime geodesic stays in
an embedded angular sector, misses the particle, is globally acausal, and
admits a Cauchy homology slice. Its length is
\[
L_{\rm HRT}=\log\!\left[4\sin(5\pi/8)\sin(\pi/16)\right]
=-0.3271718686\ldots.
\]
The proposed chiral cone classes have sweeps \(5\pi/4\) and \(\pi/8\).
The first has no regular hyperbolic-cone geodesic. Relaxing its length to
the cone-hitting infimum gives instead
\[
F_{\rm rel}=\log\!\left[4\sin(\pi/16)\right]
=-0.2479982767\ldots.
\]
The exact mismatch is \(-\log\sin(5\pi/8)>0\), in identical endpoint
normalization. The second embedded branch is more expensive in both
prescriptions. Its mismatch is different, so a common endpoint counterterm
cannot fix both class lengths.

This disproves the proposed large-angle identification, not an entropy
inequality. The all-angle mathematical functional theorem remains true,
but a physical particle theorem needs a different argument. See
`TILTED_CONE_COUNTEREXAMPLE.md` and its 80-digit checker/figure.

## What changes at extremality

For \(\alpha=(r_+-r_-)/2\to0\), fixed
\(\beta=(r_++r_-)/2>0\), boundary circumference \(P\), and a union with
\(2k\) endpoints,

\[
0\le F_{\alpha,\beta}(A)-F_{0,\beta}(A)
\le k\frac{\alpha^2P^2}{6}+P\alpha.
\]

Only two same-boundary embedded arc windings are possible; the candidate
cut set is finite. Closedness of the RT cone gives the limiting full-cone
theorem. However, the genuine hyperbolic-parabolic quotient has no invariant
AdS geodesic in its peripheral class. A positive limiting core cost is not
an attained extremal horizon circle. The unwrapped exterior has no compact
core requirement and admits an actual bounded-interval HRT theorem.

Full proofs and conventions: `CUSPS_EXTREMAL_BTZ.md`.

## Reproduction and file map

From this directory, with the dependencies in `requirements.txt`:

```bash
python code/check_equivariant_mincut.py
python code/check_cone_bigon.py
python code/check_extremal_limit.py
python code/check_tilted_cone.py
python code/bilocal_collision.py \
  --facets data/facets_n5_all372.json \
  --output collision_reproduction
```

The first three checkers write their diagnostic outputs beside their
scripts; archived run outputs are in `data/`. The final command writes its
JSON and figures to the named directory. No computation contacts a remote
service. The facet catalogue includes its original author/source URL and
subset ordering.

The exact finite-group checker covers 321 examples, including asymmetric
degenerate minima, fixed vertices, and stabilized edges. The extremal
checker includes exact matrix and horizon-acceleration identities, complete
finite annular-cut enumeration, and 155 convergence-bound instances. The
cone checker verifies the analytic counterexample at angles on both sides
of \(\pi\). None of these finite tests replaces its associated proof.

`OPEN_GAPS_UPDATE.md` gives the precise status update supported by the work.
The source audit in `SPACETIME_AUDIT.md` explains the distinct issues posed
by spin, elliptic particles, massless particles, and collisions. The
presented results do not justify the blanket headline that arbitrary
particle collisions never leave the static cone.
