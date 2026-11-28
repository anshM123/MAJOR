# Cone points, homotopy, and the Crofton step

Research audit, 28 November 2026. This note addresses the proposed point-particle extension of `PROOFS.md` §§1–4 and `fullcone/FULLCONE.md`.

## 1. Decisive result: filling the cones in the intersection classes fails at every positive-mass angle

**Proposition.** For every cone angle \(0<\theta<2\pi\), the complete hyperbolic cone

\[
C_\theta:\qquad ds^2=dr^2+\sinh^2r\,d\phi^2,
\qquad \phi\in\mathbb R/\theta\mathbb Z
\]

has a globally minimizing regular anchored geodesic \(c\), and an open positive-Liouville-measure family of regular embedded complete geodesics \(g\), for which

\[
\#(g\cap c)=2,
\qquad i_{\mathrm{filled}}([g],[c])=0.
\]

Neither curve hits the cone. Thus deleting the measure-zero set of cone-hitting geodesics does not repair the identity `(B0)` in the repository. Restricting the angle below \(\pi\) does not repair that identity when intersection classes are taken on the filled surface.

### Explicit construction

Choose

\[
\max(0,\theta-\pi)<\beta<\theta/2,
\qquad
0<\eta<\min\{\beta/2,(\pi-\theta+\beta)/2,2\pi-\theta\}.
\]

Let \(c\) have ideal endpoints at angular coordinates \(0,\beta\), using the short angular sector. Let \(g\) have physical ideal endpoints \(\eta,\beta-\eta\), but use the lift whose angular interval is

\[
(\beta-\eta,\theta+\eta).
\]

Its angular sweep is

\[
\delta_g=\theta-\beta+2\eta<\pi,
\qquad \delta_g<\theta.
\]

It is therefore a regular embedded complete geodesic. A complete regular hyperbolic geodesic spanning an angular interval \((a,b)\) of length \(\delta=b-a<\pi\) has polar equation

\[
\tanh r\cos(\phi-\mu)=\cos(\delta/2),
\qquad \mu=(a+b)/2.
\tag{1}
\]

The first intersection of \(g\) with the lift of \(c\) in \((0,\beta)\) lies in \((\beta-\eta,\beta)\). The second intersection, with the next lift of \(c\), lies in \((\theta,\theta+\eta)\). At each overlap, the corresponding developed geodesic endpoint pairs alternate, so there is exactly one transverse intersection. No further lifts of \(c\) have an angular interval meeting that of \(g\).

After projection to the cone, the two intersections remain distinct: their angular coordinates lie respectively in \((\beta-\eta,\beta)\) and \((0,\eta)\). Both physical endpoints of \(g\) lie in the same boundary cap cut off by \(c\). On the filled topological disk, \(g\) can be moved entirely into that cap, disjoint from \(c\). Hence its filled geometric intersection number is zero.

The cut \(c\) is globally minimizing. Any of its subsegments has angular difference smaller than \(\beta<\theta/2\). Among lifts of its endpoints, this is the smallest absolute angular difference. For fixed radial endpoints the hyperbolic distance is increasing with angular difference on \([0,\pi]\), by

\[
\cosh d=\cosh r_1\cosh r_2-\sinh r_1\sinh r_2\cos\Delta\phi.
\]

A filled-space minimizing path cannot pass through a cone of angle below \(2\pi\): one of the two sectors at the tip has angle below \(\pi\), and smoothing into that sector shortens it. Thus no omitted cone-hitting path is shorter. The complete anchored version follows by the cutoff limit.

All endpoint inequalities above are strict. The two intersections persist under independent small perturbations of the endpoints of \(g\). In a developed endpoint chart the Liouville density is positive,

\[
d\mu=\frac{da\,db}{4\sin^2((b-a)/2)}.
\]

Consequently the failure occurs on an open set of strictly positive measure, in a compact part of the cut. It is not a cutoff-end artifact.

### Concrete angles on both sides of the proposed threshold

| Cone angle | Cut sweep | Endpoint displacement | Test-geodesic sweep |
|---|---|---|---|
| \(\theta=\pi/2\) | \(\beta=\pi/6\) | \(\eta=\pi/48\) | \(3\pi/8\) |
| \(\theta=3\pi/2\) | \(\beta=5\pi/8\) | \(\eta=\pi/32\) | \(15\pi/16\) |

The companion code verifies the two transverse intersections at high precision and draws the curves in a topological disk rendering. The argument itself is analytic and valid for every \(\theta\in(0,2\pi)\).

## 2. Filled homology and punctured isotopy are different operations

The physically natural condition that particles are included in the homology region does **not** require forgetting the particles when labeling curve classes.

The consistent alternative is:

- Curves and test classes lie in the regular surface \(S^\circ=S\setminus Q\), with isotopies avoiding the particle marks \(Q\).
- A homology region is a region of the filled surface \(S\), so it may contain any subset of particles.
- A closed cut component bounding a disk in the filled surface may be deleted when minimizing. Deleting it can only decrease intersection counts with every test class. This deletion is an operation on admissible cuts; it must not be used to identify arbitrary arc classes while retaining an ambiguous intersection number.

For the construction above,

\[
i_{\mathrm{punctured}}([g],[c])=2.
\]

The apparent bigon contains the cone point and is not an empty bigon in \(S^\circ\).

### Minimal-position lemma in punctured classes

**Lemma.** Suppose the curves in question have regular geodesic representatives, and one is an embedded cut component. Their transverse intersections realize the geometric intersection number in the punctured surface (with the usual proper/ideal endpoint control).

**Proof.** A nonminimal transverse configuration has an empty bigon, or the corresponding monogon for a removable self-intersection. An empty bigon in the punctured surface contains no cone point. If its corner angles are \(\alpha,\beta>0\), Gauss–Bonnet gives

\[
-\operatorname{Area}(B)+(\pi-\alpha)+(\pi-\beta)=2\pi,
\]

which is impossible. Ideal-end versions follow by truncation in the hyperbolic ends. The same calculation rules out an empty monogon. The ordinary bigon criterion therefore applies. The conclusion does not require an angle threshold once regular representatives exist. \(\square\)

This repairs the *local minimal-position obstruction*, but does not by itself establish the full geodesic-current/Crofton representation: one must also justify the appropriate global space and measure of complete test geodesics, escape/properness where used, and the near-end renormalization.

## 3. What the angle \(\pi\) actually controls

The threshold controls realization of arbitrary punctured simple classes by regular geodesics, rather than the filled-topology Crofton identity.

For cone angles strictly below \(\pi\), both angular sectors at a cone-hitting limiting curve have angle below \(\pi\). One can shorten while staying on whichever side is required by the punctured class. This is the local reason minimizing sequences cannot collapse into the cone for essential simple curves. Bonsante–Schlenker prove geodesic realization of laminations in the marked-point complement. Tan–Wong–Zhang prove unique regular realizations for essential nonperipheral simple closed curves and essential normal arcs on compact cone surfaces; their \(\pi\)-angle addenda explicitly contain doubled cone-to-cone or boundary-to-cone segments.

For angles above \(\pi\), a simple punctured class can have no regular geodesic representative. On \(C_\theta\), choose ideal endpoints with

\[
0<\beta<\theta-\pi.
\]

The long-side simple arc class has angular sweep \(\theta-\beta>\pi\). No regular complete geodesic has this sweep: equation (1) requires a sweep strictly below \(\pi\). Its infimum is approached by paths through an arbitrarily small neighborhood of the cone; the limiting two radial rays pass through the singular point. Smoothing to the other side would shorten it, but changes the punctured class.

Thus:

| Claim | \(0<\theta<\pi\) | \(\pi<\theta<2\pi\) |
|---|---|---|
| Filled-space minimizing path avoids cone | Yes | Yes |
| Filled-class geodesic/Crofton minimal position for all test geodesics | **False** | **False** |
| Regular geodesics are in punctured minimal position | Yes | Yes |
| Every essential simple punctured class has a regular geodesic | Available under the cited hypotheses | False in general |

## 4. The topological contraction lemma with particle marks

Lemma A can be adapted with punctured isotopy and filled homology. Choose an auxiliary smooth negatively curved reference surface with additional funnel ends representing the particles. These auxiliary ends have **free colors**: no party boundary condition is imposed there. Physical asymptotic ends and their marked endpoints retain their original party assignment.

Place the input cuts in geodesic position for the reference metric. The cell-coloring and contraction-map construction is unchanged. Near an auxiliary particle end, the output color is whatever the contraction map assigns to the input cell color; this is allowed because the filled homology region may include or exclude that particle. No output cut must terminate at a particle. Any resulting closed component which is null-homotopic after filling the particles may be removed. This decreases every punctured intersection count and preserves the required physical boundary trace.

The proof of the contraction inequality then uses punctured intersection numbers in the auxiliary negatively curved surface, exactly as in the repository proof. This is a topological lemma. Its successful extension does not establish the physical HRT formula or the required positive measure in each chiral cone metric.

## 5. A direct disk alternative

For a single Riemannian cone disk, globally minimizing boundary-to-boundary distances satisfy both cyclic crossing inequalities by ordinary cut-and-paste. If \(a,b,c,d\) occur in cyclic order, the two minimizing paths \(a\!c\) and \(b\!d\) must intersect. Reconnecting their segments gives

\[
D(a,c)+D(b,d)\ge D(a,d)+D(b,c),
\qquad
D(a,c)+D(b,d)\ge D(a,b)+D(c,d).
\]

No curvature sign enters this argument. Positive sums of such boundary-distance matrices retain the inequalities and admit the finite positive graph construction already proved in `strips/FINITE_GRAPH.md`, after the same explicit endpoint completion.

This is a genuine mathematical bypass of the cone Crofton obstruction on a disk. The physical premise must still be established: HRT for the spacetime must equal the noncrossing matching minimum built from those interval data. One must not replace a particle-marked common class by independent left/right minima without proof.

## 6. Branched CAT(-1) cover: promising, not completed here

For finitely many cone metrics on the same marked bordered surface, choose an integer \(q\) with \(q\theta_{j,p}\ge2\pi\) at every particle and in every metric. A finite cyclic branched cover can have ramification \(q\) at all particle points; boundary monodromy absorbs the relation because the surface has boundary. The pullback metrics have nonpositive cone curvature and are locally CAT(-1).

Liouville-current length identities for nonpositively curved cone metrics are available in the literature (see Constantine below). A possible proof would build a deck-equivariant graph upstairs and descend an invariant minimum cut. This requires a separate argument; it is not proved merely by lifting the metrics.

One useful local lemma is available. Suppose a descended finite frontier meets a particle. Choose a pairing of adjacent incident rays bounding a sector of angle below \(\pi\) in one metric; such a sector exists because the total angle is below \(2\pi\). Resolve that pairing at a fixed small radius in this metric to obtain a strict length reduction. In every other metric, implement the same topological pairing in arbitrarily small collars; the excess over the cone-hitting broken path tends to zero. Thus the weighted sum strictly decreases. No side needs to have angle below \(\pi\) in *all* metrics. Applying this lemma to an invariant graph minimum requires proving that its descended finite frontier and metricwise representatives carry the same local pairing data.

## Primary sources

1. F. Bonsante and J.-M. Schlenker, *AdS manifolds with particles and earthquakes on singular surfaces*, [arXiv:math/0609116](https://arxiv.org/abs/math/0609116), Lemma 1.1 / Lemma 2.2 and Remark 2.1. Their topology keeps the particle marks in the isotopy problem. The spacetime particles in that paper are spinless.
2. S. P. Tan, Y. L. Wong and Y. Zhang, *Generalizations of McShane's identity to hyperbolic cone-surfaces*, [arXiv:math/0404226](https://arxiv.org/abs/math/0404226), Theorem 4.1 and Addendum 4.2. The boundary-arc theorem there uses normal arcs to geometric boundary, so its application to fixed ideal endpoints needs the usual end-extension argument.
3. D. Constantine, *Marked length spectrum rigidity in nonpositive curvature with singularities*, [author PDF](https://dconstantine.faculty.wesleyan.edu/files/2019/08/MLSnonpos.pdf). Relevant to the proposed branched-cover repair; its hypotheses and the extension to asymptotic ends must be checked in a complete proof.

The all-angle counterexample and the direct elementary arguments above are derived in this note. They do not claim a counterexample to the holographic entropy cone: they disprove the proposed filled-topology intermediate identity.
