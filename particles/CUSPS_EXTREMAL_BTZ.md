# G5: cusps, extremal BTZ, and the distinction between HRT and its limit

Research note, 28 November 2026. This note does not modify the repository's standing proofs.

## Status

1. **Proved:** quantitative continuity of the complete finite-party entropy vector in the nonextremal-to-extremal BTZ limit, with no escape of minimizing cut classes. Consequently the *limiting BTZ entropy prescription* lies in the full RT cone.
2. **Proved:** actual HRT for bounded intervals on the **unwrapped extremal BTZ exterior** has the disk matching/positive-graph representation, using a direct Poincare patch map and the existing disk homology argument. No horizon component or maximin theorem is needed.
3. **Proved obstruction:** a genuine hyperbolic-parabolic quotient has no invariant bulk geodesic in its peripheral class. Thus the extremal horizon cost is not the length of an attained closed extremal curve. The compact black-hole result in item 1 must not be relabelled as an attained-HRT theorem for every region.
4. **Mathematical extension:** for hyperbolic surfaces with cusps, the Crofton/current argument survives after cusp-peripheral zero-infimum components are removed, provided the same cusp labels are unmarked in every metric. This does not alone establish the spacetime Theorem H for mixed hyperbolic-parabolic ends.

The slogan "cusps are just horoball cutoffs" misses the difference between a parabolic translation-length infimum and a realized geodesic axis.

## 1. The extremal geometry and its holonomy

Set the AdS radius to one. Write

\[
\alpha=\frac{r_+-r_-}{2},\qquad
\beta=\frac{r_++r_-}{2},\qquad \phi\sim\phi+P.
\]

The ordinary nonextremal family has \(\alpha,\beta>0\). Fix \(\beta>0\) and take \(\alpha\downarrow0\). At the limit \(r_+=r_-=r_0=\beta\).

Caputa--Jejjala--Soltanpanahi, arXiv:1309.7852v2, equation (2.5), give the exact exterior-to-Poincare map

\[
w_- =\frac{e^{2r_0(\phi-t)}}{2r_0},\qquad
w_+=\phi+t-\frac{r_0}{r^2-r_0^2},\qquad
U=\frac{e^{r_0(\phi-t)}}{\sqrt{r^2-r_0^2}}.
\tag{1}
\]

It takes the exterior \(r>r_0\) to the half-Poincare patch \(U>0,w_->0\), with

\[
ds^2=\frac{dU^2+dw_+dw_-}{U^2}.
\]

The spatial identification acts as

\[
(w_+,w_-,U)\longmapsto
(w_++P,e^{2r_0P}w_-,e^{r_0P}U).
\tag{2}
\]

One chiral boundary action is a nontrivial translation, hence parabolic; the other is a dilation, hence hyperbolic. Setting one diagonal eigenvalue parameter to zero without a compensating conjugation instead produces an identity factor and describes a different quotient/near-horizon degeneration. The self-dual near-horizon quotient must not be identified with (2).

The same source, equations (2.8)--(2.14), explicitly tracks how the spatial identification changes in the near-horizon limit.

### The massless and ordinary-cusp cases are different

Taking **both** \(\alpha,\beta\to0\), with the correct parabolic normalization in both factors, gives massless BTZ. Its exterior metric is

\[
ds^2=-r^2dt^2+\frac{dr^2}{r^2}+r^2d\phi^2,
\qquad \phi\sim\phi+P.
\]

At \(t=0\), setting \(y=1/r\) identifies the spatial metric with
\((d\phi^2+dy^2)/y^2\), a parabolic cylinder. The deep end is a cusp, its peripheral length has infimum zero, and there is no positive-length horizon/core geodesic. The boundary-arc formula becomes
\(D_{0,0}=\log(|u_nv_n|/(\epsilon_p\epsilon_q))\).
The two-parameter version of (9) is

\[
0\le F_{\alpha,\beta}(A)-F_{0,0}(A)
\le kP^2(\alpha^2+\beta^2)/6+P(\alpha+\beta).
\]

This again proves full-cone membership for the limiting functional. Calling the result an attained HRT prescription requires specifying whether a homology region may extend down the infinite cusp. If only compact homology chains are permitted, a missing core component cannot silently be replaced by a zero-length curve. Ordinary cusp ends in non-elementary hyperbolic surfaces are treated separately in Section 6.

## 2. No closed peripheral geodesic in a hyperbolic-parabolic quotient

Use the matrix model \(\mathrm{AdS}_3=\{X\in M_2(\mathbb R):\det X=1\}\). The action is \(T(X)=LXR^{-1}\), with

\[
L=\begin{pmatrix}e^a&0\\0&e^{-a}\end{pmatrix},\qquad
R=\begin{pmatrix}1&s\\0&1\end{pmatrix},\qquad a>0,\ s\ne0.
\]

Its two generalized eigenspaces are

\[
V_+=\operatorname{span}(E_{11},E_{12}),\qquad
V_-=\operatorname{span}(E_{21},E_{22}),
\]

with distinct eigenvalues \(e^a,e^{-a}\). Each restriction is one nontrivial Jordan block. Every invariant subspace decomposes under the two spectral projections, since those projections are polynomials in T.

An invariant two-plane W therefore has dimensions \((2,0)\), \((0,2)\), or \((1,1)\) in these spaces. In the first two cases it is one of the row planes \(V_\pm\), on which determinant vanishes identically. In the last case, the only invariant line in each Jordan block is its eigenline, giving

\[
W=\operatorname{span}(E_{12},E_{22}),
\]

again a rank-one/determinant-zero plane. Thus **no invariant two-plane intersects \(\det X=1\)**.

Every AdS geodesic is the intersection with a two-plane through the origin (including the degenerate-plane description of null geodesics). It follows that there is no T-invariant bulk geodesic. The same proof applies to every nonzero power, since it remains hyperbolic-parabolic.

A closed geodesic in a quotient lifts to a geodesic invariant under its deck transformation. Consequently exact extremal BTZ has no closed geodesic in the peripheral winding class. In particular,

\[
\frac12(\ell_L+\ell_R)=P\beta
\]

at extremality is **not** the length of a realized closed geodesic, even though it is the positive limit of the nonextremal horizon lengths.

A direct horizon check gives the same warning. In regular ingoing coordinates the extremal metric is

\[
ds^2=-(r^2-2r_0^2)dv^2+2dv\,dr-2r_0^2dv\,d\psi+r^2d\psi^2.
\]

For a horizon circle tangent to \(\partial_\psi\), at \(r=r_0\),

\[
\nabla_{\partial_\psi}\partial_\psi
=-r_0(\partial_v+\partial_\psi)\ne0.
\]

It is a null normal acceleration, not zero. The horizon circle is therefore not an extremal curve. The absence of a bifurcation circle is geometrically substantive.

## 3. Exact boundary-arc lengths and their extremal limit

Let two boundary endpoints have differences \(\Delta\phi,\Delta t\), with a spacelike common boundary Cauchy curve. For winding n set

\[
u_n=\Delta\phi+\Delta t+nP,\qquad
v_n=\Delta\phi-\Delta t+nP.
\]

The arc length in units of the AdS radius, with a fixed compatible boundary endpoint prescription, is

\[
D_{\alpha,\beta}^{\epsilon}(n)
=\log\!\left[
\frac{\sinh(\alpha|u_n|)}{\alpha}
\frac{\sinh(\beta|v_n|)}{\beta}
\frac1{\epsilon_p\epsilon_q}
\right].
\tag{3}
\]

The two factors have the same orientation sign for a spacelike anchor pair. Formula (3) is the standard rotating BTZ geodesic formula; its \(\alpha=0\) limit agrees with arXiv:1309.7852v2, equations (3.9)--(3.11). Different endpoint renormalization conventions add only endpoint-local constants and do not affect the following estimates when held fixed across the family.

For an exact endpoint prescription, use the chiral developing coordinate
\(f_\alpha(u)=(e^{2\alpha u}-1)/(2\alpha)\), continuously extended by \(f_0(u)=u\). Give the horoball at \(f_\alpha(u_p)\) Euclidean diameter \(\epsilon_p f'_\alpha(u_p)\), and likewise in the beta factor. The exact horocycle-to-horocycle length is

\[
\ell_L(p,q)
=\log\frac{(f_\alpha(u_q)-f_\alpha(u_p))^2}
 {\epsilon_p\epsilon_q f'_\alpha(u_p)f'_\alpha(u_q)}
=2\log\frac{\sinh(\alpha|u_q-u_p|)}\alpha
-\log(\epsilon_p\epsilon_q).
\]

Half the sum with \(\ell_R\) gives (3), with no finite-radius error term. This is an exact compatible horocycle regulator, or equivalently the renormalized HRT length with its stated endpoint terms. A raw constant-r cutoff has additional vanishing local asymptotic errors and should not be substituted into exact regulated formulas without saying so.

At extremality,

\[
D_{0,\beta}^{\epsilon}(n)
=\log\!\left[
|u_n|\frac{\sinh(\beta|v_n|)}{\beta}
\frac1{\epsilon_p\epsilon_q}
\right].
\tag{4}
\]

For a constant-time interval of width \(\ell\) on the unwrapped exterior,

\[
S(\ell)=\frac{c}{6}\log\!\left[
\frac{\ell\sinh(r_0\ell)}{r_0\epsilon^2}\right],
\qquad c=\frac3{2G_N}.
\tag{5}
\]

This is the sum of one vacuum chiral contribution and one thermal chiral contribution. It is not the zero-temperature vacuum interval formula on a compact spatial circle.

## 4. No cut-class escape for a compact BTZ exterior

Fix finitely many marked endpoints on one boundary circle of a nonextremal BTZ annulus. A smooth periodic spacelike time graph satisfies \(|t'|<1\), so both lifted null boundary coordinates \(\phi\pm t(\phi)\) strictly increase by P per period.

Take endpoints with \(0<\Delta\phi<P\). An embedded same-boundary geodesic arc has only the two winding possibilities \(n=0,-1\). Indeed, if \(n\ge1\), a lift with boundary endpoints \(p,q+nP\) and its P-translate have alternating endpoints; their geodesics intersect, and the projected arc self-intersects. The case \(n\le-2\) is identical after reversal. This also follows from the classification of simple proper arcs in an annulus: each cuts off a boundary-parallel disk, and there are only the two choices of boundary segment between its endpoints.

A reduced cut has a perfect matching of the endpoint set, one of two arc classes for each pair, and at most one core component. Homological compatibility and disjointness select a subset of this finite list. This subset is topological and does not depend on \(\alpha>0\).

Hence for a region with 2k endpoints, there are at most

\[
2^{k+1}(2k-1)!!
\]

candidate reduced cuts. The estimate is deliberately loose; it establishes finiteness, which is what the limit needs. The two-sided exterior region assignment can be included in the core-side coloring, but the theorem below is stated for parties on one exterior only.

For either simple arc class, \(|u_n|,|v_n|<P\). The horizon/core cost is

\[
H_{\alpha,\beta}=P(\alpha+\beta).
\tag{6}
\]

The elementary bound

\[
0\le\log\frac{\sinh x}{x}\le\frac{x^2}{6}
\tag{7}
\]

follows termwise from \((2j+1)!\ge6^j j!\), hence \(\sinh x/x\le e^{x^2/6}\). Each arc therefore satisfies

\[
0\le D_{\alpha,\beta}^{\epsilon}-D_{0,\beta}^{\epsilon}
\le\frac{\alpha^2P^2}{6}.
\tag{8}
\]

For the formal limiting cut functional F, minimizing the finite list gives

\[
\boxed{
0\le F_{\alpha,\beta}(A)-F_{0,\beta}(A)
\le k\frac{\alpha^2P^2}{6}+P\alpha.
}
\tag{9}
\]

The last term can be removed when an \(\alpha=0\) minimizer has no core component. Equation (9) controls every required party union simultaneously. There is no unproved interchange of a limit with an infimum over infinitely many windings.

### Full-cone corollary

Choose one sufficiently fine compatible endpoint prescription throughout a compact parameter interval \(0\le\alpha\le\alpha_0\). The marked endpoints remain distinct in the parabolic normalization \(f_\alpha(x)=(e^{2\alpha x}-1)/(2\alpha)\to x\), so the endpoint horoballs/counterterms can be chosen uniformly. Equivalently, retain fixed sufficiently large positive endpoint Bell terms.

More explicitly, the pairwise horoball-disjointness condition is
\(\epsilon_p\epsilon_q< [\sinh(\alpha|u_p-u_q+nP|)/\alpha]^2\).
Its right-hand side is at least \(|u_p-u_q+nP|^2\), whose minimum over distinct marked endpoint lifts is positive and independent of alpha. Thus one fixed sufficiently small epsilon works even as the hyperbolic chiral factor degenerates to the parabolic one. The local endpoint Bell constants likewise vary continuously on this compact parameter interval.

The nonextremal vacuum full-cone theorem gives

\[
\mathbf F_{\alpha,\beta}\in\mathcal C_N^{\rm RT}\qquad(\alpha>0).
\]

The RT cone is closed. By (9),

\[
\boxed{\mathbf F_{0,\beta}\in\mathcal C_N^{\rm RT}.}
\tag{10}
\]

Every endpoint-balanced RT inequality passes directly to the renormalized limit. This is an unconditional statement about the **near-extremal limiting entropy functional** for arbitrary finite collections of intervals on a single BTZ exterior.

It does not establish that every formal core-containing minimizer is an attained extremal surface of the exact extremal quotient. Section 2 proves that such a conclusion would be false under the repository's existing complete-geodesic definition of HRT. A choice of limiting/generalized homology prescription is necessary if one wants to call (10) the entropy of the exact compact extremal black hole for all regions. In particular the whole-boundary entropy \(P\beta/(4G_N)\) cannot be supplied by an ordinary peripheral closed extremal curve there.

For a single constant-time interval the limiting thermal prescription is explicitly

\[
F_{0,\beta}([0,\ell])
=\min\{D_{0,\beta}(\ell),D_{0,\beta}(P-\ell)+P\beta\}.
\tag{11}
\]

For intervals close enough to the whole circle, the core-containing branch wins. This makes the distinction observable already for one party; it is not a technicality that appears only at large party number.

## 5. An actual HRT theorem on the unwrapped extremal exterior

Remove the periodic spatial identification and consider bounded unions of intervals on a common spacelike boundary curve. Equation (1) is now a global isometry from the entire exterior to \(U>0,w_->0\).

Each pair of finite boundary endpoints has a unique spacelike geodesic in the Poincare patch. Along it, each light-cone boundary coordinate lies between its endpoint values (the usual affine-plus-tanh geodesic parametrization). Since both endpoint values of \(w_-\) are positive, the geodesic remains in the extremal exterior.

The endpoint order in both chiral coordinates is the same. The disk argument applies to the finite endpoint set: noncrossing matchings give acausal geodesic cuts, realized on a spacelike homology slice; a finite ideal-triangle filling stays in \(w_->0\). The latter is a projective half-space condition: the defining linear functional is positive on the ideal vertices and hence on their convex ideal triangles. Completion of this slice inside the half-Poincare patch uses the same conformal-double/spacelike-extension construction as the repository's smooth-slice note. No compact horizon/core component is available or needed for bounded regions on the line. Conversely, every admissible HRT cut projects to a noncrossing matching, with no contractible geodesic loop.

Thus

\[
S_{\rm HRT}(A)=\frac1{4G_N}
\min_{\Pi\in\mathrm{NC}(\partial A)}
\sum_{\{p,q\}\in\Pi}D_{0,\beta}(p,q).
\tag{12}
\]

A useful explicit global-hyperbolicity check for the slice extension is available in the conformally flat double. Write \(w_-=X-T>0\). The graph

\[
T=-\log(1+e^{-X})
\]

has slope strictly between zero and one, lies in \(T<X\), and is a Cauchy graph of the null half-space; adding the unrestricted transverse U-coordinate gives a spacelike Cauchy surface. To check the last claim, put \(f=T+\log(1+e^{-X})\). This increases strictly along future causal curves. A future end either meets \(T=X\), where \(f=\log(1+e^X)>0\), or has \(T\to+\infty\), where \(X>T\) forces \(f>0\). Toward the past, \(T\to-\infty\); the null coordinate \(X-T\) is nondecreasing toward the past and bounded below by its positive value at any chosen point. Thus \(f<0\) sufficiently far in the past, including the potentially dangerous \(X\to-\infty\) direction. Every inextendible causal curve crosses the graph once. The function \(f\) itself is not asserted to have all levels Cauchy.

The reflected finite acausal ribbons can therefore be extended in a globally hyperbolic conformal double. The conformal boundary U=0 is an ordinary reflection seam, exactly as in the existing smooth-slice construction. The prescribed boundary Cauchy curve can be retained: its union with the finite geodesic cut is acausal by the disk argument, and the completion is performed away from its boundary collar.

Finally, either apply the finite disk graph theorem directly, or take \(\alpha\to0\) in the finite Catalan matching problem. On a fixed finite endpoint set let \(U_*\) be the largest relevant null-coordinate difference; then

\[
0\le F_{\alpha,\beta}(A)-F_{0,\beta}(A)
\le k\alpha^2U_*^2/6.
\]

With the endpoint terms retained, the full entropy vector lies in \(\mathcal C_N^{\rm RT}\). This establishes a genuine arbitrary-party HRT result for the extremal **unwrapped exterior**, not merely a horizon-entropy convention.

## 6. Ordinary hyperbolic cusp surfaces

Let X be a finite-type complete hyperbolic surface with cusp ends and surviving conformal funnel ends, and mark endpoints only on the latter. Suppose cusp ends carry no separately prescribed party trace.

A cusp is \(\{y\ge Y\}/(x\sim x+1)\) with metric \((dx^2+dy^2)/y^2\). Its peripheral class has length infimum zero, approached by horocycles \(y=Y'\to\infty\), and has no closed geodesic. A cusp-peripheral component of a homology cut can be removed while preserving the trace on every surviving conformal boundary: reverse the color of the cusp tail. Accordingly, remove these components from the class system, just as null-homotopic and paired parallel components are removed in the vacuum proof.

All remaining closed classes have hyperbolic holonomy and a unique closed geodesic; marked-endpoint proper arcs also have unique geodesic representatives. Because curvature is still -1 everywhere in X, the usual linked-lift minimal-intersection proof is unchanged. This differs fundamentally from cone points, which are filled positive-curvature points at finite distance.

The Liouville Crofton formula is local and remains valid. Geodesics with an end at a parabolic fixed point form a measure-zero set. For a geometrically finite second-kind surface, the limit set has Lebesgue measure zero, so almost every complete geodesic has both ends in the surviving ordinary set and defines an allowed proper test class. A nontrivial test class crossing between distinct ordinary-set lifts has finite Liouville mass by the surviving funnel-geodesic transversal, exactly as in PROOFS Section 5.3. The finite-mass endpoint classes in a single funnel have the unchanged local calculation of FULLCONE Section 1. The endpoint-divergence classes at marked ordinary points likewise have the same local horoball calculation as before.

There is no escape of a nonperipheral bounded-length minimizing cut into the cusp. A curve entering \(y\ge Y'\) from \(y=Y\) and returning has length at least \(2\log(Y'/Y)\). A cut with bounded truncated length must therefore stay within a fixed cusp height. After deleting peripheral components, all components return; bounded candidates consequently lie in a compact core with fixed endpoint collars, where discreteness gives finitely many relevant homotopy classes.

Use the **exact finite-arrangement proof** of FULLCONE Section 3, rather than its finite-current approximation proof: no new filling-family lemma at punctures is needed. More explicitly:

1. Bounded-length properness just proved gives a minimizing reduced class for each of the finitely many required regions. Choose their geodesic representatives in one reference cusp metric.
2. Their finite union has finitely many crossings and cells. Every component has bounded cusp depth, so sufficiently deep in each cusp there is just one cell. That cell is an ordinary internal graph vertex: the cusp carries no terminal condition.
3. Transfer the positive finite mass of each ordinary-ended test class to its normalized Liouville measure in the reference metric, as in FULLCONE equation (3.1). Minimal intersection is unchanged. Consequently each arrangement-edge weight is finite and nonnegative, bounded by the length of an entropy-minimizing cut containing it.
4. The frontier of any union of selected cells can be resolved away from the cusp ends. The reference geodesic of its reduced class has no larger intersection with any test class. This proves the graph lower bound. Each selected entropy minimizer proves the reverse bound.
5. The endpoint Bell edges are precisely those of the standing theorem, because every marked endpoint lies on an unchanged ordinary conformal end.

Thus, for a positive weighted sum of several marked cusp metrics with the same cusp set, the regulated cut vector is exactly a finite graph mincut vector for sufficiently fine compatible endpoint cutoffs. Cusp-peripheral loops have zero intersection with every ordinary-ended test class and are removable in every metric, so their deletion loses no cost. This proves a full-cone statement for the resulting **geometric cut functional**. The statement here assumes genuine surviving funnel ends (and is therefore non-elementary); the elementary parabolic cylinder has already been handled as a BTZ limit in Section 1.

Establishing a spacetime HRT theorem still requires the corresponding Theorem H. In particular, mixed hyperbolic-parabolic ends fail the "same removable cusp class" condition: their other chiral length is positive, and Section 2 shows that this cost has no invariant geodesic realization.

## 7. Reproducible checks

Run `python check_extremal_limit.py` in this directory. It writes `cusp_validation.json` and checks:

- the exact hyperbolic-parabolic linear action and determinant-null invariant planes;
- the exact nonzero null acceleration of an extremal horizon circle;
- complete finite candidate enumeration for five disjoint intervals on a nonconstant periodic spacelike boundary graph: the number of compatible reduced cuts for \(k=1,\ldots,5\) interval components is \(2,6,20,70,252\);
- the uniform bound (9) in 155 finite-party/parameter instances, including all 31 nonempty unions;
- MMI and the five-party cyclic inequality on the sample vectors (diagnostics, not substitutes for the theorem);
- a large interval for which the formal extremal core-containing branch is strictly cheaper, illustrating why attained-HRT and limiting entropy must be distinguished.

The length/entropy normalization is explicit: the script works in length units and all entropy values are obtained by division by \(4G_N\).

## 8. Precise G5 update supported by this work

- Bounded-interval HRT on the unwrapped extremal BTZ exterior: **closed** by the explicit patch map and disk argument.
- Compact BTZ, arbitrary finite interval parties, entropy defined as the nonextremal limit: **closed**, with the uniform estimate (9).
- Ordinary cusp-surface positive-length/current/graph functional with unmarked common cusp ends: **closed at the mathematical functional level**.
- Compact extremal BTZ, attained closed HRT curve for every homology class: **false for the peripheral class**.
- General mixed-parabolic multiboundary spacetime Theorem H: **not closed by the cusp argument**; its prescription and end structure require a separate treatment.

## Primary references

- P. Caputa, V. Jejjala, H. Soltanpanahi, *Entanglement Entropy of Extremal BTZ*, arXiv:1309.7852v2, equations (2.3)--(2.14), (3.9)--(3.11), Appendix A. https://arxiv.org/html/1309.7852v2
- M. Banados, M. Henneaux, C. Teitelboim, J. Zanelli, *Geometry of the 2+1 Black Hole*, arXiv:gr-qc/9302012. https://arxiv.org/abs/gr-qc/9302012
- F. Bonsante, K. Krasnov, J.-M. Schlenker, *Multi Black Holes and Earthquakes on Riemann Surfaces with Boundaries*, arXiv:math/0610429, for cusp/geodesic-boundary end data. https://arxiv.org/abs/math/0610429
- N. Bao et al., *The Holographic Entropy Cone*, arXiv:1505.07839, for closedness and graph characterization. https://arxiv.org/abs/1505.07839

The cusp-functional extension uses standard hyperbolic-surface facts explicitly listed in Section 6; it is not a claim that the cited black-hole paper already proves this entropy-cone extension.
