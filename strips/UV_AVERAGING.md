# Relative-cutoff averaging and the asymptotic HRT comparison

This note supplies the UV step in the symmetry proof without assuming that
finite-cutoff HRT surfaces converge, and without applying convexity to
renormalized perimeter. The reference trace is the trace of the invariant
reduced-maximin surface itself.

There are two logically separate results:

1. An exact finite-cutoff spatial inequality with a boundary-mismatch penalty.
2. Sufficient asymptotic conditions under which both that penalty and the
   cutoff mismatch in the null-focusing comparison tend to zero.

The first result is a BV statement. The second uses same-anchor extremal
asymptotics plus a uniformly spacelike conformal slice and a regular local
null representative. Bare NEC does not imply all of those boundary
regularity conditions.

## 1. Notation and the exact relative inequality

The bulk dimension is \(D=d+1\) and the transverse dimension is \(m=d-2\). Work first
with a finite transverse torus of coordinate volume V. Let Σ be an invariant
spacelike slice containing \(X=\pi^{-1}\gamma\), where γ is a reduced spatial minimum on
\(\Sigma/T^m\). At this stage, no full-dimensional spatial minimality is assumed.
The needed property is only that invariant admissible cuts obtained by
changing γ in the interior cannot have smaller reduced length.

Choose an invariant defining function z. Let

\[
 \Sigma_\epsilon=\Sigma\cap\{z\ge\epsilon\},\qquad
 B_\epsilon=\Sigma\cap\{z=\epsilon\}.
\]

Other infrared boundaries, if present, carry the same fixed trace for all
competitors. The quotient domain under consideration is compact; an infrared
limit can be treated separately. Let \(R_X\) be X's homology region on Σ, and let
\(a_\epsilon\) be its characteristic trace on \(B_\epsilon\). This trace is invariant.

For any finite-perimeter homology region R in \(\Sigma_\epsilon\), let \(b_\epsilon\) be its boundary
trace on \(B_\epsilon\), and write \(P_\epsilon(R)\) for its interior perimeter, excluding the
prescribed outer boundary. Define the mismatch set

\[
 E_\epsilon=\{a_\epsilon\ne b_\epsilon\}\subset B_\epsilon.
\]

**Relative-cutoff averaging lemma.** For almost every sufficiently small
regular cutoff ε,

\[
 \boxed{
 \operatorname{Area}(X\cap\Sigma_\epsilon)
 \le P_\epsilon(R)+
 \mathcal H^{d-1}_{h_\Sigma}(E_\epsilon).
 }
 \tag{1}
\]

The lemma requires no smallness of \(E_\epsilon\). It compares against nonsymmetric R.

### Proof

First replace R in a layer \(\epsilon<z<\epsilon+\eta\) by the reference region \(R_X\), retaining R
for \(z>\epsilon+\eta\). Call the resulting region \(R_\eta\). Its trace at \(z=\epsilon\) is exactly \(a_\epsilon\).
The BV gluing formula, or direct counting of the interface, gives

\[
 P_\epsilon(R_\eta)
 \le P_\epsilon(R)
 +\mathcal H^{d-1}(E_{\epsilon+\eta})
 +P(R_X;\epsilon<z<\epsilon+\eta).
\]

Here \(E_{\epsilon+\eta}\) is the mismatch between the traces from the two sides of the
gluing hypersurface. For fixed ε>0, the metric is regular in this thin layer.
The last term tends to zero as \(\eta\to0\), and the BV trace theorem gives a sequence
\(\eta\to0\) for which the mismatch integrals converge to
\(\mathcal H^{d-1}(E_\epsilon)\). Thus it suffices to prove the desired spatial
lower bound for a competitor with exact trace \(a_\epsilon\).

Average \(\chi_{R_\eta}\) over the torus, obtaining an invariant f with the same trace.
Convexity of total variation and weighted coarea imply

\[
 P_\epsilon(R_\eta)
 \ge V\int_0^1P_\lambda(\{f>s\};\sigma_\epsilon)\,ds,
 \qquad \sigma_\epsilon=\Sigma_\epsilon/T^m.
\]

For almost every s, the reduced level region has exactly the trace of \(R_X\) at
\(z=\epsilon\). Glue it to \(R_X\)'s own exterior collar \(0<z<\epsilon\). The result is an admissible
invariant homology region with the same asymptotic anchor as γ. There is no
interface perimeter because the traces agree. By the assumed spatial
minimality of γ, its interior weighted perimeter is at least that of
\(\gamma\cap\{z\ge\epsilon\}\). Hence

\[
 P_\epsilon(R_\eta)\ge
 \operatorname{Area}(X\cap\Sigma_\epsilon).
\]

Let \(\eta\to0\) along the trace-convergent sequence. This proves (1). ∎

This argument does not invoke any extremal-surface theorem at the artificial
timelike cutoff boundary. It uses an ordinary Riemannian variational problem
on the already chosen slice and gluing to the reference surface's own collar.

## 2. A quantitative boundary-band criterion

Assume the conformal spatial metric \(\bar h_\Sigma=z^2h_\Sigma\) extends to a
positive-definite metric at the entangling edges, with bounded coefficients
and inverse. Suppose \(E_\epsilon\) lies in finitely many edge bands whose total
coordinate width normal to the invariant edge planes is at most \(\delta_\epsilon\). On
\(B_\epsilon\) the physical (d−1)-dimensional measure scales as \(\epsilon^{-(d-1)}\), so

\[
 \frac{\mathcal H^{d-1}(E_\epsilon)}{V}
 \le C\epsilon^{-(d-1)}\delta_\epsilon.
 \tag{2}
\]

Consequently, the sufficient trace-control condition is

\[
 \boxed{\delta_\epsilon=o(\epsilon^{d-1}).}\tag{3}
\]

Under (3), the relative averaging error tends to zero. In particular,
\(\delta_\epsilon=O(\epsilon^d|\log\epsilon|^q)\), for fixed q, gives an error
\(O(\epsilon|\log\epsilon|^q)\). The more intrinsic sufficient condition is simply
\(\mathcal H^{d-1}(E_\epsilon)/V\to0\); bands are a convenient way to verify it.

## 3. Why two same-anchor extremal surfaces provide the needed scale

Let Y be a smooth competing extremal surface with the same asymptotic anchor
as X. In common boundary-adapted normal graph coordinates, their embeddings
obey

\[
 U_Y(z,y)-U_X(z,y)=o(z^{d-1}).
 \tag{4}
\]

The usual smooth or polyhomogeneous extremal expansion gives the stronger
\(O(z^d|\log z|^q)\) estimate. Coefficients below the first undetermined order are
fixed by the common anchor and the same ambient metric. Derivative estimates
of the corresponding order are understood, uniformly on the finite torus
and the finite set of edge components.

A direct Lorentzian primary source is Jonathan Sorce, *Holographic
entanglement entropy is cutoff-covariant*,
[arXiv:1908.02297](https://arxiv.org/abs/1908.02297), §3.1, Lemma 2 and the
discussion following eq. (3.32). Sorce uses d for **bulk dimension**. In the
notation of this note, his result fixes the embedding through order
\(z^{D-2}=z^{d-1}\), hence the difference (4). With smooth Taylor remainders this
is \(O(z^d)\). The result compares extremal surfaces in a single ambient
spacetime; it does not require claiming that all divergent coefficients are
state independent between different spacetimes.

For standard Einstein FG expansions, the same leading exponent can be seen
from the normal graph equation. Its homogeneous radial principal part is

\[
 \partial_z^2 U-\frac{d-1}{z}\partial_zU+\cdots=0,
\]

whose indicial roots are 0 and d. The anchor fixes the zero-root data and the
lower local recursion, leaving the first freely adjustable coefficient at
order \(z^d\). This is only an explanation of the scale; equation (4) is the
precise asymptotic input used below.

Equation (4) does **not** directly apply to a null representative, which is
usually not extremal. The next step propagates it through the null projection.

## 4. Propagating the estimate to the maximin slice

Assume X is already known to be fully extremal, as established by invariant
stationarity in the symmetry proof. Let Σ contain X. In a conformal collar,
write Σ as \(F=0\), where F is a regular time-graph defining function. Require
that Σ be uniformly spacelike for \(\bar g=z^2g_D\). Equivalently, after a
fixed normalization, dF remains uniformly inside the timelike covector cone.

Because F vanishes on X and has bounded derivatives, (4) gives

\[
 F(Y(z,y))=o(z^{d-1}).\tag{5}
\]

Normalize the relevant null normals k to Y by a regular conformal temporal
coordinate, for example \(dT(k)=\pm1\). Uniform spacelikeness supplies a uniform
transversality bound

\[
 |dF(k)|\ge c>0.\tag{6}
\]

Null geodesic paths are unchanged by conformal rescaling. The conformal
metric and the asymptotically regular normal fields have bounded local
coefficients. The implicit-function theorem applied to F along each local
null generator therefore yields an intersection with Σ at conformal flight
parameter

\[
 s(z,y)=o(z^{d-1}),\tag{7}
\]

and every coordinate displacement is of this order. In the usual stronger
expansion class, these bounds are \(O(z^d|\log z|^q)\). The derivative bounds in
the extremal expansion keep the local projection a graph; its Jacobian is
uniformly bounded and tends to the identity in the radial and tangential
coordinates.

Provided that the actual null representative agrees with this local sheet
in the asymptotic collar, its intersection with \(z=\epsilon\) differs from X's
intersection only in normal-coordinate bands of width

\[
 \delta_\epsilon=o(\epsilon^{d-1}).\tag{8}
\]

Substitution into (1)–(3) gives

\[
 \operatorname{Area}(X_\epsilon)
 \le\operatorname{Area}(\widetilde Y_\epsilon)+o(V).
 \tag{9}
\]

The constants may depend on the chosen competitor Y. Uniformity over the
entire infinite collection of possible extremal competitors is unnecessary:
the comparison is made with each one separately.

### Why an arbitrary smooth spacelike slice is insufficient

In pure Poincaré AdS, take the vertical extremal plane Y at \(t=x=0\) and a slice
Σ given by \(t=az\), \(0<|a|<1\). The normal null sheets have \(t=\pm x\) at fixed z, so their
intersections with Σ shift by \(x=\pm az\). The mismatch is \(O(z)\), not \(o(z^{d-1})\).
This slice is conformally smooth and uniformly spacelike. What fails is the
additional property used in (5): it does not contain a same-anchor extremal
reference surface with the required common asymptotic expansion.

In the symmetry proof Σ does contain X, and X is already extremal before the
representative comparison is invoked. That fact removes this obstruction.

## 5. Sufficient conditions for the actual representative to be local

It is not correct to identify the local normal-sheet intersection with an
arbitrary globally defined representative without a locality argument. A
concrete sufficient geometric criterion is as follows:

- The relevant conformal completion has a smooth causal collar and closed
  causal relation J in the compact causal region under consideration. Each
  anchor has a basis of causally convex, strongly causal neighborhoods, so
  causal curves with nearby endpoints cannot make a remote excursion and
  bypass the local normal neighborhood.
- The conformal closure of Y is a compact embedded spacelike C² submanifold
  with boundary equal to the specified entangling edges. Polyhomogeneous
  expansions with the required C² regularity are sufficient.
- Each anchor point p is causally unrelated to every other point of that
  closure: \(J^\pm(p)\cap\overline Y=\{p\}\). A conformally regular acausal
  homology/Cauchy slice containing \(\overline Y\) is a sufficient stronger
  condition.
- The representative is the usual boundary of the corresponding future or
  past domain, assembled from the appropriate null normal sheets; its
  homology trace agrees with A at conformal infinity. It is proper, with
  conformal boundary exactly the specified anchors and no extra boundary
  accumulation. The containing slice is Cauchy in the relevant causal
  region, so a generator cannot meet it first locally and then a second
  time deeper in the bulk.

Here is the localization argument. If remote generators could replace the
local sheet arbitrarily close to an anchor p, there would be target points
\(q_n\to p\) and source points \(y_n\) in a fixed closed portion of \(\overline Y\) away
from p, with \(q_n\) causally related to \(y_n\). Compactness gives \(y_n\to y\ne p\). Closedness
of J then makes p causally related to y, contradicting the third condition.
Thus all sufficiently nearby relevant sources lie in a common small local
neighborhood of p.

Inside that neighborhood, regular spacelikeness and the normal exponential
map give the usual short null collar without focal or cut points. This is a
local tubular-neighborhood statement for a smooth spacelike submanifold; one
may extend the conformal metric and the surface smoothly across the boundary
to formulate it. The intersections in (7) occur at a displacement much
smaller than z itself. Their local footpoints therefore stay on the physical
\(z>0\) portion of Y, and the surface's boundary cannot intrude into this short
normal construction. The appropriate local null sheet is consequently part
of the global representative. This establishes the locality used in §4.

For geometries where these sufficient conditions have not been checked,
retain the following explicit, weaker boundary hypothesis instead:

> In a collar of each anchor, the representative is given by the short
> local null-normal projection of Y to the conformally uniformly spacelike
> slice Σ, with bounded projection derivatives.

Neither this condition nor its sufficient compactness/acausality criterion
is being inferred solely from NEC. They concern the regular boundary behavior
of the representatives already used in the maximin framework.

## 6. Focusing with the same radial cutoff

Even after locality is established, a representative generator can cross
\(z=\epsilon\), so the unregulated area comparison must be aligned with the cutoff.
The flight estimate in §4 supplies that alignment.

Take a monotone upper envelope of the \(o(z^{d-1})\) displacement estimate. There
is a quantity \(\rho_\epsilon=o(\epsilon^{d-1})\) such that any generator contributing to
\(\widetilde Y\cap\{z\ge\epsilon\}\) originates in
\(Y\cap\{z\ge\epsilon-\rho_\epsilon\}\). For sufficiently small ε,
\(\rho_\epsilon\ll\epsilon\). Focusing on the selected generators, with caustic generators discarded,
then gives

\[
 \operatorname{Area}(\widetilde Y_\epsilon)
 \le\operatorname{Area}(Y\cap\{z\ge\epsilon-\rho_\epsilon\}).
\]

The asymptotic area element on Y is bounded by a constant times
\(z^{-(d-1)}\,dz\,d^my\). Hence the additional radial layer contributes at most

\[
 C V\rho_\epsilon\epsilon^{-(d-1)}=o(V).
\]

Therefore the correctly aligned focusing estimate is

\[
 \boxed{
 \operatorname{Area}(\widetilde Y_\epsilon)
 \le\operatorname{Area}(Y_\epsilon)+o(V).
 }\tag{10}
\]

In the \(O(z^d|\log z|^q)\) expansion class, both the radial-layer error here and
the seam error in §2 are \(O(V\epsilon|\log\epsilon|^q)\).

## 7. The asymptotic symmetry conclusion

Combining (9) and (10) yields

\[
 \operatorname{Area}(X_\epsilon)
 \le\operatorname{Area}(Y_\epsilon)+o(V).
\]

X and Y are both extremal with the same asymptotic anchor, so their local
divergences agree and their area difference has a finite regulator-covariant
limit. Dividing by V and passing to that limit proves

\[
 \boxed{
 \operatorname{Area}^{\rm ren}(X)/V
 \le\operatorname{Area}^{\rm ren}(Y)/V.
 }
\]

There is no need to separately define a renormalized area for the generally
nonextremal representative: the two cutoff inequalities suffice. Since X is
already fully extremal, it is an HRT representative in the specified regular
asymptotic class.

This closes the UV averaging step under a concrete boundary-control condition.
It does not prove exact HRT or SSA statements at a finite radial wall, and it
does not identify truncated continuum HRT surfaces with solutions of a
different finite-boundary variational problem.

## Primary sources and relationship to earlier notes

- J. Sorce, *Holographic entanglement entropy is cutoff-covariant*,
  [arXiv:1908.02297](https://arxiv.org/abs/1908.02297), §3.1, Lemma 2, eqs.
  (3.16) and (3.32), for same-anchor extremal asymptotics and finite area
  differences. Its symbol d denotes bulk dimension.
- C. R. Graham and E. Witten, *Conformal Anomaly Of Submanifold Observables In
  AdS/CFT Correspondence*, [hep-th/9901021](https://arxiv.org/abs/hep-th/9901021),
  and C. R. Graham, *Volume and Area Renormalizations for Conformally Compact
  Einstein Metrics*, [math/9909042](https://arxiv.org/abs/math/9909042), for the
  standard local minimal-submanifold expansion and the first undetermined
  coefficient. The Lorentzian same-spacetime statement needed here is
  supplied directly by Sorce.
- `SYMMETRY.md` supplies the exact KK area identity, compact-torus
  BV argument, and full extremality of the invariant lift. The present note
  replaces a vague appeal to a compatible UV variational regulator with the
  explicit relative-cutoff inequality (1) and boundary estimates.
