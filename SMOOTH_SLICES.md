# Smooth uniformly spacelike homology slices for wormhole HRT cuts

**G1 regularity upgrade.** This note proves Proposition S+ below for the
geometries and HRT candidates of `PROOFS.md` §§6.1–6.3. It replaces the
smooth-slice restriction to whole-boundary regions by a statement covering
all finite interval regions. The original geodesic cut is kept pointwise.

**Base revision:** `329d7430eee1440b3fd0ad974b273fcca95db886`.
**Proof date:** 27 November 2026.
**Inputs retained:** the spacetime description and proper boundary action
in §§6.1, 6.3b; the acausal pleated realization in Proposition P; and
Theorems C and E. The external extension input is Bernal–Sánchez,
Theorem 1.1 and Remark 4.14 [BS]. The argument below verifies its hypotheses
after a conformal doubling construction. This is a proof of the G1
regularity upgrade, not a new audit of all the other results in the paper.

## 1. Statement and meaning of uniformity

Let $M=\widetilde D/\Gamma$ be an orientable nonextremal locally
$\mathrm{AdS}_3$ wormhole in §6.1, with convex-cocompact marked left and
right holonomies. Write $\overline M$ for its spacetime completion by the
timelike conformal boundary. A region $A$ is a finite union of intervals
and whole components of a boundary Cauchy curve. For a cut class
$c\in\mathcal C(A)$, let $\gamma_c$ denote the corresponding finite union
of simple spacelike geodesic arcs and closed geodesics in $M$.

**Proposition S+ (smooth homology realization, including arc ends).**
There exists a compact, smooth, spacelike Cauchy hypersurface with boundary
$\overline\Sigma\subset\overline M$ such that:

1. $\overline\gamma_c\subset\overline\Sigma$, including every anchored
   endpoint. The physical geodesics themselves are unchanged.
2. Its boundary is a smooth spacelike Cauchy curve of each boundary
   cylinder, through the points of $\partial A$.
3. The cut bounds a homology region on $\overline\Sigma$ with the
   corresponding boundary region $A'$. The intervals of $A'$ have the
   same endpoints, marking, and boundary domain of dependence as $A$.
4. For any fixed smooth representative $\bar g$ of the quotient conformal
   metric and any smooth auxiliary Riemannian metric $h$, there is
   $\epsilon>0$ such that

$$
\bar g(v,v)\geq\epsilon h(v,v)
\qquad\text{for every }v\in T\overline\Sigma.                 \tag{1}
$$

The interior lift is $\Gamma$-invariant and satisfies (S1)–(S3). The
homology region is a smooth surface with the usual corners at the
anchoring points where its two boundary pieces meet.

The statement also applies to an already given admissible candidate:
every finite simple geodesic cut on an acausal Cauchy slice can be kept
fixed while replacing that slice by one with these regularity properties.
For the empty cut, use an unconstrained smooth Cauchy slice in the same
construction.

Equation (1) is the precise uniformity relevant here: a positive margin
on the completed quotient slice, for the fixed spacetime and cut. It
also gives the corresponding physical inequality after multiplying both
metrics by the inverse square of a defining function. It does not assert
one constant for all moduli or for all translates measured in a fixed
non-equivariant frame on the universal cover.

The boundary freedom is exactly that in `PROOFS.md` §6.3. The curve $C'$
may change away from $\partial A$, while retaining its interval labels.
A merely continuous original boundary trace cannot in general be the
boundary of a surface smooth up to conformal infinity, so retaining that
trace pointwise would be a different assertion.

## 2. A smooth conformal double exists

Work first in the conformal cover used in §6.3b. Its metric is

$$
g_E=-d\theta^2+d\chi^2+\sin^2\chi\,d\phi^2,
\qquad 0\leq\chi\leq\pi/2.
$$

Reflect the spatial hemisphere through its equator. The doubled metric
is smooth on $S^2\times\mathbb R$, and the reflection
$\chi\mapsto\pi-\chi$ preserves time orientation. The boundary diamonds
are open in the timelike equator. Consequently two copies of
$\mathcal E_\Gamma$, glued along their boundary diamonds, form an open
subset of this conformal double.

It remains to check compatibility with the holonomy, rather than assume
that a coordinate reflection descends. Introduce the ambient quadratic
form

$$
Q(X,w)=-\det X+w^2,
\qquad X\in M_2(\mathbb R),\quad w\in\mathbb R.               \tag{2}
$$

The positive null rays, in the appropriate Einstein cover, carry the
Lorentzian conformal structure. On the continuous positive matrix lift
of the invisible domain, a point $\det X=1$ is represented by
$[(X,1)]$. The reflected copy uses $[(X,-1)]$, and the common boundary is
$w=0$. This is a statement in the positive lift/causal cover, not a
claim that all of $PSL(2,\mathbb R)$ embeds by the formula $(X,1)$.

The transformations

$$
(X,w)\longmapsto(g_LXg_R^{-1},w),\qquad
\iota[(X,w)]=[(X,-w)]                                      \tag{3}
$$

preserve $Q$ and commute. Their lifts agree with the existing holonomy
action on $\mathcal E_\Gamma$ and with spatial reflection on the other
half. Thus the holonomy extends smoothly and conformally across the
boundary. Any overall matrix-sign choice in the positive-lift
description commutes with the same reflection.
More explicitly, the lifted reflection fixes the entire lifted equator
pointwise. Its commutator with a chosen lifted holonomy is a deck
transformation fixing an equatorial point, and is therefore the identity.

The action on the doubled domain is free and proper. Indeed, folding a
compact subset to the original half gives a compact subset there, so
properness follows from the proper action on $\mathcal E_\Gamma$,
boundary included. Freeness follows in the same way. Therefore

$$
\widehat M=D\overline M
$$

is a smooth Lorentzian conformal manifold without boundary, with the two
original copies as closed halves. Choose a smooth metric $g_d$ in this
conformal class. Since $\iota$ is conformal,
$g_d+\iota^*g_d$ is a positive scalar multiple of $g_d$ and is invariant
under $\iota$. We henceforth use this reflection-invariant representative.
No causal cone has been altered.

## 3. Folding proves global hyperbolicity

Let $F:\widehat M\to\overline M$ be the fold. In collar coordinates
adapted to reflection,

$$
F(y,z)=(y,|z|).
$$

**Folding lemma.** A future causal locally Lipschitz curve in the double
folds to a future causal locally Lipschitz curve in the original half.
Nonconstant segments stay nonconstant, and an inextendible curve folds
to an inextendible curve.

*Proof.* Away from the seam, the derivative of $F$ is a time-orientation
preserving isometry. On the set where the normal coordinate $z$ of an
absolutely continuous curve is zero, $z'=0$ almost everywhere. There
the folded tangent equals the original boundary-tangent causal vector.
This also covers curves with infinitely many contacts with the seam.
A folding fiber is one point or two isolated points, so a connected
nonconstant curve cannot fold to a constant.

If the folded future tail had an endpoint $p$, the original tail would
also have an endpoint. For $p$ on the boundary, convergence to $p$ in
either copy is convergence to the same glued point. For $p$ in the
interior, the tail eventually lies in one of the two disjoint copies of
a neighborhood of $p$, and converges to that copy of $p$. It could then
be continued timelike in the double. This contradicts inextendibility.
The past argument is identical. $\square$

Let $\overline\Sigma_0$ be the original pleated Cauchy slice supplied by
Proposition P and Theorem C. For an already given admissible candidate,
take its original acausal Cauchy slice instead. Its double

$$
\widehat\Sigma_0=F^{-1}(\overline\Sigma_0)
$$

is Cauchy: every inextendible causal curve folds to one meeting
$\overline\Sigma_0$. It is also acausal. A causal curve between distinct
points of $\widehat\Sigma_0$ would fold either to a causal curve between
distinct points of $\overline\Sigma_0$, or to a nonconstant closed causal
curve. Both are impossible in the original Cauchy spacetime. Thus the
intersection is unique. The ordinary Cauchy-surface characterization of
global hyperbolicity applies to the boundaryless double.

This is the same global doubling mechanism used in [AFS], §5.2. Their
general construction first modifies the metric near a timelike boundary.
The locally AdS conformal double above is already smooth, so that
modification is unnecessary here.

## 4. Geodesic ends become smooth compact curves

Take positive rank-one endpoint representatives $p,q$ of one lifted arc,
and set

$$
N=N(p,q)=\operatorname{tr}(\operatorname{adj}(p)q)>0.
$$

The physical geodesic is

$$
X(t)=\frac{e^tp+e^{-t}q}{\sqrt N}.
$$

Near the endpoint $p$, use $u=e^{-t}$ and rescale the positive null ray:

$$
\boxed{[(X(t),1)]=[(p+u^2q,\,u\sqrt N)].}                   \tag{4}
$$

The identity $\det(p+u^2q)=u^2N$ verifies the null-cone equation exactly.
Formula (4) is smooth for both signs of $u$, and reflection is exactly
$u\mapsto-u$. The derivative of its displayed null representative is
$(2uq,\sqrt N)$ and has $Q$-norm $N$ for every $u$. Its endpoint tangent is represented by $(0,\sqrt N)$,
whose $Q$-norm is $N>0$. It is transverse and conformally orthogonal to
the timelike boundary. Near $q$, the same argument with $v=e^t$ gives
$[(q+v^2p,v\sqrt N)]$.

Thus the compactified arc and its reflected copy join to a smooth
embedded spacelike circle. Closed components of the cut give two
disjoint circles. Since a cut has finitely many components and exactly
one arc end at each anchoring point, its double

$$
K=F^{-1}(\overline\gamma_c)
$$

is a compact smooth embedded spacelike one-manifold, with no branch
points. It is acausal because $K\subset\widehat\Sigma_0$. In particular,
this checks causality between different components and different mirror
copies, not just the tangent norm of each separate curve.

## 5. A compact acausal ribbon and the extension theorem

Orient each circle of $K$. The normal bundle has signature $(1,1)$ and
is time-oriented. The orientation of the ambient three-manifold and of
the circle fixes a spacelike orientation as well. Project a smooth
future timelike field orthogonally to the circle and choose the
orthogonal unit spacelike normal $n$ with this orientation. This gives
a smooth global section, even when parallel transport has boost
holonomy.

For sufficiently small $\delta_0>0$ the normal exponential map gives a
compact embedded ribbon with boundary,

$$
H_{\delta_0}=
\{\exp_x(r n_x):x\in K,\ |r|\leq\delta_0\}.
$$

Its tangent plane along $K$ is positive definite. Compactness and
openness of spacelikeness make the whole ribbon spacelike after
shrinking $\delta_0$.

**Compact ribbon lemma.** Some nested ribbon $H_\delta$ is globally
acausal in $\widehat M$.

*Proof.* Otherwise, for $\delta_j\downarrow0$ choose distinct causal
pairs $x_j,y_j\in H_{\delta_j}$. Subsequence limits $x,y$ lie in $K$.
Closedness of $J^+$ in a globally hyperbolic spacetime gives
$y\in J^+(x)$. Acausality of $K$ forces $x=y$. Near this point the
**fixed** smooth spacelike ribbon $H_{\delta_0}$ is a zero level of a
smooth temporal coordinate. Strong causality supplies a smaller
causally convex neighborhood. For large $j$ a causal curve between
$x_j$ and $y_j$ must stay in that neighborhood, where a temporal
coordinate cannot take equal values at its two endpoints. Contradiction.
$\square$

The fixed-ribbon qualification matters: this is not a claim that
arbitrary small perturbations of a noncompact acausal pleated surface
remain acausal.

Now apply [BS], Theorem 1.1, to $H_\delta$: a compact, spacelike, acausal
submanifold with boundary in a globally hyperbolic spacetime lies in a
smooth spacelike Cauchy hypersurface. Remark 4.14 explicitly permits the
disconnected prescribed submanifold needed here. Obtain

$$
\widehat\Sigma\supset H_\delta\supset K,
$$

a smooth spacelike Cauchy hypersurface of $\widehat M$. If $K$ is empty,
one may apply the same theorem to any sufficiently small compact
spacelike disc in a causally convex neighborhood. Reflection symmetry
of $\widehat\Sigma$ is not needed.

## 6. Restriction gives the required boundary Cauchy slice

Set

$$
\overline\Sigma=\widehat\Sigma\cap\overline M.               \tag{5}
$$

At the seam, a spacelike two-plane cannot equal the timelike tangent
two-plane of $\partial\overline M$. The two hypersurfaces therefore
intersect transversely. Equation (5) is consequently a smooth spacelike
surface with boundary, containing the complete compactified cut.

Every inextendible causal curve of $\overline M$ is inextendible in
$\widehat M$. If it had an endpoint in the double, that endpoint would
belong to the closed original half. An interior endpoint permits a
timelike continuation inside the half; a boundary endpoint permits a
timelike continuation tangent to its timelike boundary. Both contradict
the original inextendibility. The curve therefore meets
$\widehat\Sigma$ exactly once, necessarily in $\overline\Sigma$.
Thus $\overline\Sigma$ is Cauchy and acausal. Its boundary is a Cauchy
curve of each boundary cylinder by the same argument.

Everything was constructed on the quotient. Lifting automatically
gives $\Gamma$-invariance. Theorem E supplies (S1)–(S3), including the
correct signed spacelike-separation condition. In particular no
averaging of embedded surfaces or matching of accumulating pleats is
required.

## 7. Compactness follows from the peripheral relation

Uniformity requires compactness of the completed quotient slice. It
does not follow by assuming that the affine map $\Psi$ in Proposition P
is equivariant. The following argument avoids that assumption.

**Compactness lemma.** Every connected oriented acausal Cauchy surface
$\overline S$ in the completion of the marked wormholes of §6.1 is
compact.

*Proof.* There are finitely many boundary cylinders, so

$$
\partial\overline S=B_1\sqcup\cdots\sqcup B_b,
\qquad b\geq1.
$$

Theorem E identifies the lifted slice with a graph over a closed disc
with the limit-set points of its boundary removed. This domain is
contractible by radial contraction to an interior point. Hence it is
the simply connected cover and
$\pi_1(\overline S)=\Gamma$.

Each $B_i$ represents a conjugate of the marked peripheral element
$h_i$ or its inverse. For nontrivial $\Gamma$, it is a generator rather than a proper power:
a Cauchy circle in a boundary cylinder has degree one in its spatial
circle factor. The compact marked oriented surface underlying the
holonomies has the relation

$$
\prod_{j=1}^{g}[a_j,b_j]\,h_1\cdots h_b=1.
$$

Consequently $\sum_i[h_i]=0$ in $\Gamma_{\rm ab}$. For suitable
$\varepsilon_i=\pm1$, the nonzero element
$\sum_i\varepsilon_i[B_i]\in H_1(\partial\overline S;\mathbb Z)$
maps to zero in $H_1(\overline S;\mathbb Z)$.

If $\overline S$ were noncompact, ordinary relative homology would give

$$
H_2(\overline S,\partial\overline S;\mathbb Z)=0.             \tag{6}
$$

For example, oriented Poincaré–Lefschetz duality identifies this group
with $H_c^0(\overline S;\mathbb Z)=0$. Equivalently, in a locally
finite triangulation the orientation multiplicity of a relative
two-cycle is constant across all interior edges. A nonzero constant on
a connected noncompact surface has noncompact support and cannot be
represented by a finite ordinary chain. This proves (6) directly.

The exact sequence of the pair would then make

$$
H_1(\partial\overline S;\mathbb Z)
\longrightarrow H_1(\overline S;\mathbb Z)
$$

injective, contrary to the nonzero peripheral relation. Thus
$\overline S$ is compact. The disc and annulus cases use respectively
one trivial peripheral element and two inverse peripheral elements.
$\square$

Apply the lemma to $\overline\Sigma$. The continuous positive function
$v\mapsto\bar g(v,v)$ on its compact $h$-unit tangent bundle has a
positive minimum. This is precisely (1), including the conformal
boundary and the rotating closed components.

## 8. Homology is preserved pointwise along the cut

Choose a smooth future timelike vector field $V$ on $\overline M$
tangent to its timelike boundary. Local choices exist and their
positive partition-of-unity sum stays in the convex future timelike
cone. Rescale by a smooth positive function to make its norm bounded
in a complete auxiliary Riemannian metric. The resulting field is
complete, preserves the boundary, and has inextendible timelike orbits.

Every orbit meets $\overline\Sigma_0$ and $\overline\Sigma$ exactly once.
Let

$$
\Phi:\overline\Sigma_0\longrightarrow\overline\Sigma
$$

be the orbit-intersection map. It is bijective. Its continuity follows
from the implicit-function theorem at the smooth target surface,
which is transverse to $V$, followed by restriction to the possibly
only Lipschitz original surface. The source is compact by the preceding
lemma and the target is Hausdorff, so $\Phi$ is a homeomorphism.

Every point of $\overline\gamma_c$ already lies on both slices.
Uniqueness of an orbit's intersection therefore gives
$\Phi|_{\overline\gamma_c}=\mathrm{id}$. The map preserves each
boundary component and fixes the endpoints. It sends each oriented
boundary interval to the interval on the new boundary Cauchy curve
with the same marked endpoints; orientation is preserved by the
timelike flow. Since both boundary curves are increasing graphs in
the boundary null coordinates, these are exactly the intervals
specified in §6.3.

Transport the old homology region by $\Phi$. A homeomorphism preserves
closure, interior, and frontier. The transported region consequently
has interior frontier $\gamma_c$ and boundary at infinity $A'$. The
curve $\gamma_c$ is smooth and properly embedded in the new smooth
slice, and its tangent is transverse to the boundary at each endpoint
by (4). Thus the region is smooth with the usual endpoint corners.
This proves all parts of Proposition S+. $\square$

## 9. Exact smooth strips around rotating closed components

The global proof already establishes existence. An explicit model also
identifies the rotating pleat degeneration as a feature of that
construction, rather than an obstruction from holonomy.

Let $a=\ell_L>0$, $b=\ell_R>0$, and set

$$
L=\frac{a+b}{2},\qquad \eta=\frac{a-b}{2},\qquad
q=\frac\eta L=\frac{a-b}{a+b},\qquad |q|<1.
$$

Use the unit-speed long axis and a rotating spacelike normal

$$
c(s)=\begin{pmatrix}0&e^s\\-e^{-s}&0\end{pmatrix},\qquad
n(s)=\begin{pmatrix}e^{qs}&0\\0&-e^{-qs}\end{pmatrix}.
$$

Then

$$
\boxed{\mathcal X(s,r)=\cosh r\,c(s)+\sinh r\,n(s)}          \tag{7}
$$

has determinant one and obeys the exact holonomy relation
$h\mathcal X(s,r)=\mathcal X(s+L,r)$. Direct polarization of
$g=-\det(dX)$ gives

$$
\boxed{\mathcal X^*g
=dr^2+\bigl[1+(1-q^2)\sinh^2r\bigr]ds^2.}                  \tag{8}
$$

It is spacelike for all $s,r$. In particular, a thin strip descends
smoothly through the rotating closed geodesic.

For $\Delta=s-s'$, its exact pairwise separation is

$$
\begin{aligned}
\tfrac12N(\mathcal X(s,r),\mathcal X(s',r'))
&=\cosh r\cosh r'\cosh\Delta
 -\sinh r\sinh r'\cosh(q\Delta)\\
&=\cosh(r-r')\cosh(q\Delta)\\
&\quad+\cosh r\cosh r'
       [\cosh\Delta-\cosh(q\Delta)]>1
\end{aligned}                                                   \tag{9}
$$

for distinct parameters. There is no time-winding ambiguity: its
matrix entries satisfy
$(X_{12}-X_{21})/2=\cosh r\cosh s>0$, so its conformal time lies
in one interval of length $\pi$. Thus (9) is genuine acausality on
the relevant causal cover, and it descends for the cyclic holonomy.

There is also a quantitative invariant margin. The equivariant field

$$
T(s)=\operatorname{diag}(e^{qs},e^{-qs})
$$

is a unit timelike observer tangent to AdS along the strip. With the
associated Riemannian metric $g_R=g+2T^\flat\otimes T^\flat$,

$$
g(v,v)\geq\frac{1-q^2}{1+q^2}\,g_R(v,v)
=\frac{2ab}{a^2+b^2}\,g_R(v,v),
\qquad v\in T\mathcal X.                                  \tag{10}
$$

For a general wormhole the strip is taken sufficiently thin to stay
inside the invisible domain and avoid other translates. Global
acausality of the finite prescribed collection follows from the
compact ribbon lemma, not from disjointness of separately constructed
strips alone.

## 10. Explicit collars around arc endpoints

In a Poincaré patch set

$$
g=z^{-2}(dz^2+dx^2-dt^2).
$$

For endpoints $(0,0)$ and $(X,T)$ with $d^2=X^2-T^2>0$, the branch
near the first endpoint is

$$
x_\gamma(z)=\frac X2\left(1-\sqrt{1-\frac{4z^2}{d^2}}\right),
\quad
t_\gamma(z)=\frac T2\left(1-\sqrt{1-\frac{4z^2}{d^2}}\right). \tag{11}
$$

Both functions are smooth and even in $z$, and are $O(z^2)$. Let a
chosen smooth boundary trace be $t=b(x)$, with $b(0)=0$ and
$|b'|\leq\kappa<1$ in a compact coordinate interval. Choose a smooth
bump $\chi$ supported in the interval and equal to one near zero.
Define

$$
\boxed{F(x,z)=b(x)+\chi(x)
 \bigl[t_\gamma(z)-b(x_\gamma(z))\bigr].}                    \tag{12}
$$

For small $z$, the graph $t=F(x,z)$ has the prescribed boundary trace
and contains the exact geodesic. Its derivatives satisfy

$$
F_x=b'+O(z^2),\qquad F_z=O(z),\qquad
F_x^2+F_z^2\leq\kappa^2+Cz^2.
$$

Choosing the collar width so that $Cz^2\leq(1-\kappa^2)/2$ gives

$$
\bar g|_{t=F}=dx^2+dz^2-dF^2
\geq\frac{1-\kappa^2}{2}(dx^2+dz^2).                        \tag{13}
$$

Finitely many distinct endpoints admit disjoint such neighborhoods
and a common positive minimum width. These formulas answer the local
boundary construction directly; the double supplies the global
extension without gluing derivatives across infinitely many pleats.

An independent rotating-BTZ check is also useful. With $z=1/r$, its
conformal metric is exactly

$$
\bar g=
\frac{dz^2}{(1-r_+^2z^2)(1-r_-^2z^2)}
-[1-(r_+^2+r_-^2)z^2]dt^2
-2r_+r_-z^2\,dt\,d\phi+d\phi^2.                            \tag{14}
$$

All coefficients are even in $z$. Writing
$\bar g=A(z^2)dz^2+h_{ab}(z^2)dy^a dy^b$ and using the conserved
momenta $p_a$ of a unit-speed spacelike geodesic gives

$$
\frac{dy^a}{dz}
=\pm z\frac{\sqrt A\,h^{ab}p_b}
 {\sqrt{1-z^2h^{cd}p_cp_d}}.                                \tag{15}
$$

The right side is odd and smooth near zero. Thus $y^a(z)$ is even and
smooth, consistent with (4). Its conformal tangent norm is
$A/(1-z^2h^{cd}p_cp_d)>0$. The global proof uses the ambient
construction, so it needs no additional assumption about a global BTZ
coordinate system for the whole wormhole.

## 11. Consequence for the HRT minimum and the paper

Let $S^\varepsilon_{\rm sm}(A)$ denote the infimum with the additional
requirement that the completed homology slice is smooth and uniformly
spacelike in the sense of (1). Every such candidate is admissible in
the original definition, so
$S^\varepsilon_{\rm sm}(A)\geq S^\varepsilon(A)$.

Choose a minimizing cut class from Theorem H. Proposition S+ realizes
its identical geodesic union on a smooth uniformly spacelike Cauchy
slice with the same homology. The regulator and every regulated length
are unchanged. Therefore the reverse inequality holds, and

$$
\boxed{
S^\varepsilon_{\rm sm}(A)=S^\varepsilon(A)
=\frac{1}{4G_N}
 \min_{c\in\mathcal C(A)}
 \frac{\ell^\varepsilon_{X_L}(c)+\ell^\varepsilon_{X_R}(c)}{2}.
}                                                               \tag{16}
$$

The minimum is attained in the smooth class. This closes the
smooth/uniformly-spacelike part of G1 for interval regions as well as
whole-boundary regions. Theorem 2 inherits the same upgrade because
its input is precisely this entropy functional.

For the trivial group (global AdS), the same argument starts with the
whole conformal hemisphere and its globally hyperbolic double
$S^2\times\mathbb R$. The original acausal graph is Cauchy because
conformal time is unbounded along inextendible causal curves. After
extension and restriction, vertical timelike curves identify the new
slice with a graph over the closed hemisphere, giving (S1)–(S3)
directly. Its boundary circle has trivial image in the fundamental
group. Thus this case does not require the nontrivial-group formulation
of Theorem E.

The existing geometric scope is unchanged: orientable vacuum locally
AdS3 quotients, hyperbolic nonextremal holonomies, and finite regions
on a common acausal boundary Cauchy curve. Other questions recorded
in `OPEN_GAPS.md`, including the achronal class with null pieces and
homology regions not required to lie in a complete Cauchy slice, are
separate from this regularity result.

## 12. Verification and sources

The proof is analytic. The scripts in `smooth_slices/` independently
check the exact local algebra for (7)–(15) using SymPy. They do not
certify the global topology or replace the causal arguments above.

- **[BS]** A. N. Bernal and M. Sánchez, *Further results on the
  smoothability of Cauchy hypersurfaces and Cauchy time functions*,
  Letters in Mathematical Physics **77** (2006), 183–197.
  [arXiv:gr-qc/0512095](https://arxiv.org/abs/gr-qc/0512095).
  Theorem 1.1, Lemma 3.9, and disconnected-case Remark 4.14.
- **[AFS]** L. Aké, J. L. Flores, and M. Sánchez, *Structure of globally
  hyperbolic spacetimes-with-timelike-boundary*, Revista Matemática
  Iberoamericana **37** (2021), 45–94.
  [arXiv:1808.04412v4](https://arxiv.org/abs/1808.04412),
  [DOI:10.4171/rmi/1201](https://doi.org/10.4171/rmi/1201).
  Sections 3.2 and 5.2; the doubling and Cauchy-flow arguments.
- **[Repo]** `anshM123/MAJOR`, `PROOFS.md` at
  [the base revision](https://github.com/anshM123/MAJOR/blob/329d7430eee1440b3fd0ad974b273fcca95db886/PROOFS.md),
  §§6.1–6.3b, Proposition P, and Theorem H. These supply the standing
  geometric description, marked peripheral subgroups, and the initial
  acausal Cauchy realization.

The reduction to a smooth conformal double, reflected-arc formula, and
its application to the prescribed HRT cut are the additions in this
note. General smooth Cauchy extension and the doubling method are
credited to the references above; no separate novelty claim is made
for those established theorems.
