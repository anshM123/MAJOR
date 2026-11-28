# A continuous-angle cone theorem for the common-class length functional

Research checkpoint, 28 November 2026.

## Status and scope

This note proves a curve-functional extension of the repository's positive-graph construction to hyperbolic cone surfaces whose cone angles lie in `(0, pi)`. The angles need not be rational multiples of `2 pi`. The proof removes the positive-curvature singularities from the calculation: small cone neighborhoods are replaced by smooth negatively curved funnels that no admissible simple geodesic visits. The replacement preserves the relevant marked lengths exactly.

The theorem does **not** identify this functional with HRT in a particle spacetime. That is a separate spacetime theorem, including common-class marking, acausality, and homology. In particular this note does not cover collision vertices, generic spinning particles, or assert a fixed chiral-pair description across a collision.

The new ingredients below are elementary collar and funnel-replacement lemmas, together with the extension of renormalized Crofton to smooth negative curvature using a local boundary subtraction. The finite arrangement graph argument is the one in `fullcone/FULLCONE.md`, restated to make the changed hypotheses explicit.

## 1. Statement

Let `S` be a compact connected oriented surface with nonempty physical boundary, let `Q` be a finite set of interior marked points, and let `P` be finitely many marked endpoints on the physical boundary. Give `S \ Q` finitely many marked complete hyperbolic cone metrics `h_j`, with conformally compact hyperbolic ends at the physical boundary and cone angles

\[
0<\theta_{j,q}<\pi \qquad(j=1,\ldots,m;\ q\in Q).
\]

The filled metric completions include the points `Q`. Markings identify the particle points, physical boundaries and marked endpoints. We require the standard geodesic realization property for reduced simple particle-marked classes; for angles below `pi` this is the usual cone-surface realization theorem, discussed in Section 4 below.

Cut classes use **isotopy in `S \ Q`**, with physical endpoints fixed. Their homology regions lie in the **filled surface `S`**. Thus particle points can lie on either side, but the winding class of an arc around a particle is not discarded. No cut ends on a particle. Remove contractible and particle-peripheral closed components when reducing a cut; these cannot improve a minimum because particles carry no boundary label.

Assign disjoint physical boundary pieces to finitely many parties and a purifier. Let `C(A)` be the cut classes with boundary trace `A`. For positive weights `w_j`, let

\[
F^\epsilon(A)=\min_{c\in\mathcal C(A)}
  \sum_{j=1}^m w_j\ell_j^\epsilon(c),
\]

using fixed compatible local hyperbolic endpoint cutoffs. Then, for every fixed finite party system:

1. The minimum is attained, and the relevant cut classes are finite below every fixed renormalized length bound.
2. There are a finite nonnegative graph-cut function `G(A)` and endpoint coefficients `b_p(A)=[p in boundary A]` such that

\[
F^\epsilon(A)=G(A)+\sum_{p\in P}t_p(\epsilon)b_p(A),
\qquad
t_p(\epsilon)=\Big(\sum_jw_j\Big)\log\frac1\epsilon+O(1).
\tag{1}
\]

3. For sufficiently small `epsilon`, every `t_p(epsilon)` is nonnegative. Adding an edge between the two party labels adjacent to `p` realizes the endpoint term. Therefore `F^epsilon` is exactly the min-cut vector of one finite nonnegative graph and belongs to the full static RT cone.

The same conclusion applies to whole-boundary parties, where no endpoint counterterms occur. Auxiliary particle funnels in the proof are unlabelled; they are graph internal vertices, not additional parties or a fixed purifier.

## 2. Uniform exclusion of a small cone collar

**Lemma 1.** Suppose a cone has angle `theta < pi` and an embedded radial chart of radius `R`,

\[
h=dr^2+\sinh^2r\,d\phi^2,
\qquad \phi\in\mathbb R/\theta\mathbb Z.
\]

Every embedded complete regular geodesic avoids the open collar

\[
r<\rho(R,\theta):=
\operatorname{arctanh}\big(\tanh R\cos(\theta/2)\big)>0.
\tag{2}
\]

Here complete includes closed geodesics and complete arcs with endpoints at physical conformal infinity.

**Proof.** If a regular geodesic enters this collar, extend its local portion in the radial chart to its radial minimum `r_0`. A radial segment hitting the cone is excluded by regularity. In the developed hyperbolic chart its equation is

\[
\tanh r\cos(\phi-\phi_0)=\tanh r_0.
\]

The two points at angular positions `phi_0 +/- theta/2` have the same radius

\[
r_* = \operatorname{arctanh}
       \frac{\tanh r_0}{\cos(\theta/2)}<R.
\]

They are identified in the cone and the two tangent directions are distinct. Thus the geodesic self-intersects before leaving the radial chart, a contradiction. The strict positivity of `cos(theta/2)` is precisely why this uniform argument stops at `theta=pi`. `QED`

Choose, independently for every metric and particle, a replacement radius `a` strictly below (2). Every component of every regular simple cut lies entirely outside the modification.

## 3. Explicit replacement by a smooth negatively curved funnel

**Lemma 2.** For any `a>0`, the punctured cone collar `0<r<a` can be replaced by a smooth complete funnel collar with the following properties:

- the metric is unchanged for `r>=a` and in a neighborhood just inside `a`;
- its curvature is smooth and pinched between two negative constants;
- sufficiently far inside the new end it is exactly a hyperbolic funnel;
- there is an orientation-preserving, marking-preserving `1`-Lipschitz diffeomorphism from the new open collar to the old punctured collar, equal to the identity outside the modification.

**Construction.** Keep the angular coordinate of period `theta` and use a new radial coordinate `t in R`. Choose a nonzero smooth bump `eta>=0` supported in `(a/3,2a/3)`. Put

\[
E_-:=\int e^{-s}\eta(s)\,ds,
\qquad E_+:=\int e^s\eta(s)\,ds.
\]

Because the support is positive, `E_+>E_->0`. Pick

\[
E_+^{-1}<c<E_-^{-1},\qquad H(s)=c\eta(s).
\]

Define, with the integral zero above the support,

\[
f(t)=\sinh t+\int_t^a\sinh(s-t)H(s)\,ds.
\tag{3}
\]

For `t>=2a/3`, `f(t)=sinh t`. For `t<=a/3`,

\[
f(t)=A e^t+B e^{-t},\quad
A=\tfrac12(1-cE_-)>0,\quad
B=\tfrac12(cE_+-1)>0.
\tag{4}
\]

For `t in [a/3,a]`, both `sinh t` and the integral in (3) are nonnegative, with the former positive. Thus `f>0` everywhere. The new metric is

\[
\widehat h=dt^2+f(t)^2d\phi^2.
\]

It is complete towards `t=-infinity`, matches the original metric smoothly, and

\[
K_{\widehat h}=-\frac{f''}{f}=-1-\frac{H}{f}\le-1.
\tag{5}
\]

Since `H` is compactly supported and `f` is positive, the curvature is also bounded below by a finite negative constant. Equation (4) is

\[
f(t)=2\sqrt{AB}\cosh(t-t_*),
\qquad t_* =\tfrac12\log(B/A),
\]

so the new asymptotic end is exactly hyperbolic.

To construct the domination map, note that `f(t)>=sinh t` for `t>=0`. Let `m=min_R f>0` and choose `0<c_0<min(a/3,arsinh m)`. Choose a smooth increasing diffeomorphism

\[
r:\mathbb R\longrightarrow(0,\infty)
\]

with `r(t)=t` for `t>=c_0`, `0<r(t)<=c_0` for `t<=c_0`, and `0<r'(t)<=1`. Such a map is obtained by integrating a smooth function `q in (0,1]` that equals `1` near and above `c_0` and has integral `c_0` on `(-infinity,c_0]`. The transition can be placed sufficiently far to the left to achieve this integral. Then

\[
(r')^2dt^2+\sinh^2(r(t))d\phi^2
\le dt^2+f(t)^2d\phi^2.
\tag{6}
\]

For `t<=c_0` the angular inequality uses `sinh r<=sinh c_0<m`; for `t>=c_0` it uses `r=t` and (3). This proves the claimed `1`-Lipschitz marking map. `QED`

No orbifold angle or finite cover enters this construction. It also does not require the cone angles of different metrics to agree: the same topological end is used, with its own angular normalization in each metric.

## 4. Exact preservation of the simple marked lengths

Apply Lemma 2 inside a collar smaller than (2). Denote the resulting smooth metric on the surface with extra funnel ends by `hhat_j`.

Let `c` be a reduced simple particle-marked class, and let `c_j` be its regular embedded geodesic representative in `h_j`. By Lemma 1, `c_j` misses every modified collar. Hence it is still a complete geodesic for `hhat_j`. Completeness and negative curvature of `hhat_j` imply that it is the unique geodesic representative of the same class, including for prescribed physical ideal endpoints. Consequently

\[
\widehat\ell_j(c)=\ell_j(c),
\qquad
\widehat\ell_j^\epsilon(c)=\ell_j^\epsilon(c).
\tag{7}
\]

The domination map (6) also gives `lhat_j(c)>=l_j(c)` directly at the level of infima, since the maps preserve punctured classes and are the identity near the physical endpoint cutoffs.

The geodesic realization input for `theta<pi` is standard: Bonsante--Schlenker, Lemma 2.2 and Remark 2.1, and Tan--Wong--Zhang, Theorem 4.1. The latter states the boundary-arc result for normal arcs to geometric boundary; it is not literally a fixed-ideal-endpoint statement. The extension here follows by truncating the prescribed ideal rays in the unchanged hyperbolic end, applying length minimization in the fixed simple class, and taking the limit. Compactness on the core is standard; a cone-hitting minimizing limit can be shortened in the side prescribed by the punctured class because **both** complementary cone sectors have angle below `pi`. The local collar argument above then prevents regular embedded minimizers from approaching a cone as the endpoint truncation is removed. Negative curvature away from the points rules out empty geodesic bigons and proves uniqueness.

For a cut multicurve, disjointness is preserved under geodesic realization by the same no-empty-bigon argument in the punctured surface. Thus (7) holds for the whole cut, with the same marked multicurve class in every metric.

The auxiliary ideal boundary circles carry no prescribed label. A particle-peripheral closed cut in the modified surface separates only a free end from the physical boundary and can be deleted, changing the free end's color. Its positive length decreases the cost in every metric. More general components enclosing only unlabelled ends can likewise be removed when redundant. Consequently the infimum over modified-surface cuts with free auxiliary-end colors is exactly the original filled-homology infimum over reduced particle-marked cuts.

The modified surfaces are smooth negatively curved asymptotically hyperbolic surfaces of finite type. On such a surface, each essential marked class has a unique geodesic representative. For a fixed finite endpoint set, there are only finitely many reduced simple classes below a fixed renormalized length bound: outside a compact core the end segments are fixed by their endpoint and first core entrance; a length bound bounds the associated deck translates in the proper universal cover. Closed components obey the corresponding ordinary length-spectrum finiteness, and a disjoint reduced multicurve has topologically bounded component number. This proves attainment and finiteness for the common-class minimum by (7).

## 5. Crofton and escape for the smooth replacement metrics

The replacement is useful because the necessary global flow statement is available in precisely this smooth setting. Lefeuvre, *Boundary rigidity of negatively-curved asymptotically hyperbolic surfaces*, arXiv:1805.05155, Section 2.1.2, Lemma 3.1 and Remark 3.1, provides:

1. the past- and future-trapped tails have zero Liouville measure;
2. outside those tails, geodesics have two endpoints at conformal infinity and form the scattering-geodesic space;
3. in transverse coordinates along a geodesic segment the Liouville current is `sin(theta) dtheta dt`, up to the orientation/normalization convention.

With the normalization in the repository this gives

\[
\ell_h(\sigma)=\int \#(g\cap\sigma)\,d\overline\mu_h(g)
\tag{8}
\]

for compact rectifiable test segments, with multiplicity. It also gives the standard Crofton formula for a finite union of such segments. In particular, this proof does **not** assume that a Gromov-boundary homeomorphism preserves a measure-zero limit set. The a.e. escape statement is an independent smooth negative-curvature theorem.

For a geodesic cut representative and an a.e. proper scattering geodesic, the actual intersection number equals the topological minimum in the punctured/funnel topology. An empty bigon would contradict negative curvature by Gauss--Bonnet. Properness and generic transversality ensure a finite number of intersections with every fixed finite cut arrangement.

The test geodesics may end on auxiliary funnels. This is necessary: their Crofton contributions encode the metric's response to particle winding. Their endpoints do not create new parties or fix the auxiliary-end homology color.

## 6. Local endpoint subtraction in variable negative curvature

The global metric need not be hyperbolic for the endpoint part of `FULLCONE.md` to work. A small neighborhood of each physical ideal endpoint is still exactly hyperbolic.

Choose disjoint local hyperbolic charts around the marked endpoints. For one endpoint `p`, place it at infinity in the upper half-plane. Choose `A_p>0` so that the entire outside-half-disk neighborhood

\[
U_p=\{z\in\mathbb H^2:|z|>A_p\}
\]

is contained in that chart and its ideal boundary contains no other marked endpoint. Define `T_p` to be the local geodesics with ideal endpoints

\[
a<-A_p,\qquad b>A_p.
\]

Every such semicircle lies in `U_p`; it is also a complete global geodesic by uniqueness in negative curvature. Topologically it cuts off just the endpoint `p`, so for every cut class

\[
\#(g\cap c)=e_p(c)\in\{0,1\}\quad(g\in T_p).
\tag{9}
\]

When `e_p(c)=1`, the terminal part of the cut is a vertical line `x=x_0`, with `|x_0|<=A_p`. Its crossing with `g=(a,b)` has height

\[
y=\sqrt{(x_0-a)(b-x_0)}.
\]

The local Liouville density is the standard hyperbolic one. Put `u_0=A_p+x_0`, `v_0=A_p-x_0`. For `Y^2>u_0v_0`, the exact crossing integral is

\[
\int_{T_p}\#(g\cap c^Y)\,d\overline\mu_h(g)
=e_p(c)\left\{
\frac12\log\big[(Y^2+u_0^2)(Y^2+v_0^2)\big]
-\log(2A_pY)\right\}.
\]

Therefore the elementary integral in `FULLCONE.md` gives

\[
\int_{T_p}\#(g\cap c^Y)\,d\overline\mu_h(g)
=e_p(c)\big(\log Y-\log(2A_p)+o(1)\big).
\tag{10}
\]

The constant is independent of `x_0`, hence independent of the cut class. Removing the disjoint families `T_p` leaves a positive measure

\[
\mu'_h=\overline\mu_h\big|_{\mathcal G\setminus\bigcup_p T_p}.
\]

By (8), (9), monotone convergence and (10),

\[
\ell_h^{\rm ren}(c)
=\int i(\tau(g),c)\,d\mu'_h(g)
-\sum_p e_p(c)\log(2A_p).
\tag{11}
\]

The positive integral is finite for every reduced cut class. Restoring the exact local horoball cutoff gives

\[
\ell_h^\epsilon(c)
=\int i(\tau(g),c)\,d\mu'_h(g)
+\sum_p e_p(c)\log\frac1{2A_p\epsilon}
\tag{12}
\]

once the cutoff lies in the terminal local charts. The choice is uniform over all cut classes and all replacement radii: `|x_0|<A_p`, since otherwise the vertical geodesic would remain entirely in the local unmarked boundary neighborhood and have its second endpoint there. A geodesic with neither endpoint `p` cannot reach the horoball `y>A_p` without one of its ideal endpoints being in this unmarked local neighborhood. Thus `Y>A_p` is a uniform terminal cutoff, and the error in (10) is uniformly `O(A_p^2/Y^2)`. As in the repository, changing endpoint normalization changes only the endpoint constants.

Unlike the original globally hyperbolic argument, it is not necessary to remove the *entire* topological endpoint class. The local family `T_p` contains all its divergent mass; the finite residual mass in that class is harmless and remains in the positive measure. This is the reason the same construction survives a compact change of curvature.

Push `mu'_h` to the countable set of relative test classes `tau`, dropping classes with zero intersection against every cut under consideration. Class masses are finite: classes whose geodesics enter a compact core are bounded by Crofton for finitely many compact cross-sections; classes confined to an asymptotic collar have the usual local endpoint-coordinate estimate, with the only divergent endpoint families removed in (10). Equivalently, finiteness follows whenever a class has positive intersection with one reduced cut, from (11); classes invisible to every cut can be omitted.

## 7. One finite graph for all the common-class minima

For the modified metrics put

\[
\nu(\tau)=\sum_j w_j\mu'_j(\tau),
\qquad
G(A)=\min_{c\in\mathcal C(A)}i(\nu,c).
\]

Equation (12) is (1), with

\[
t_p(\epsilon)=\sum_jw_j\log\frac1{2A_{j,p}\epsilon}.
\]

Choose a smooth constant-curvature reference metric on the same surface with the particle ends replaced by funnels. For every test class with positive finite `nu(tau)`, choose any probability measure `pi_tau` absolutely continuous with respect to reference Liouville measure and supported on reference geodesics in that class. Even if a full reference class has infinite mass, a compact positive-measure subfamily supplies such a probability measure. Set

\[
\widehat\nu=\sum_\tau\nu(\tau)\pi_\tau.
\]

For each physical party union `A`, choose a minimizing class `c_A` and put its reference geodesic representative into one finite arrangement. Define an edge capacity to be the `nuhat`-measure of test geodesics crossing that arrangement edge, counted with multiplicity. It is finite because the edge lies on one selected minimizing cut with finite `i(nu,c_A)`.

The arrangement cells form the vertices. Physical boundary cells receive their party labels. Cells meeting an auxiliary funnel have **no additional prescribed label**. Because no cut endpoint is placed there, the whole auxiliary boundary circle lies in one cell.

Each selected `c_A` gives a graph cut of capacity `G(A)`. Conversely, resolve the frontier of any admissible union of cells into an embedded cut. Reference geodesics are Liouville-a.e. proper, transverse to the arrangement and avoid its vertices. For each such test geodesic, the smoothing disks and endpoint collars can be chosen small enough to miss it. The topological output class is independent of these radii, and

\[
i(\tau,\delta)\le
\sum_{e\subset\partial U}\#(g\cap e).
\]

Integration gives `i(nu,delta)<=capacity(U)`. Deleting peripheral or redundant components preserves the physical trace and only decreases the intersection count. Therefore the graph minimum equals `G(A)` for every union. Adding the nonnegative endpoint edges for small `epsilon` proves the theorem.

## 8. What this establishes and what remains separate

The established result is a full finite-graph theorem for the positive weighted common-class length functional of any finite collection of marked hyperbolic cone metrics with cone angles below `pi`. It covers a continuous open angle range and every number of parties. It keeps particle winding in the class labels and fills particles only for homology.

For a spacetime whose HRT entropy is independently proved to equal this functional, the full RT-cone conclusion follows immediately. Showing that equality, with an admissible acausal homology slice, is a separate Lorentzian task. Particle collisions cannot be included by merely treating their pre- and post-collision chiral data as a single fixed pair.

For angles in `[pi,2pi)`, the uniform embedded-geodesic collar fails and some punctured simple classes have cone-hitting infima. A finite chosen set of regular attained minimizing classes can still be preserved by the same funnel replacement and domination argument. The companion `ALL_ANGLE_RELAXED_FUNCTIONAL.md` proves the full graph-cone conclusion for the **relaxed common-class infimum functional at arbitrary angles**, by shrinking the replacement collars and using closedness of the fixed-party graph cone. Identifying that relaxed functional with a regular HRT minimum still requires a proof of joint regular attainment and the Lorentzian common-class formula.

## Primary references

1. F. Bonsante and J.-M. Schlenker, *AdS manifolds with particles and earthquakes on singular surfaces*, [arXiv:math/0609116](https://arxiv.org/abs/math/0609116), Lemma 2.2 and Remark 2.1.
2. S. P. Tan, Y. L. Wong and Y. Zhang, *Generalizations of McShane's identity to hyperbolic cone-surfaces*, [arXiv:math/0404226](https://arxiv.org/abs/math/0404226), Theorem 4.1 and Addendum 4.2.
3. T. Lefeuvre, *Boundary rigidity of negatively-curved asymptotically hyperbolic surfaces*, Comment. Math. Helv. 95 (2020), 129--166, [arXiv:1805.05155](https://arxiv.org/abs/1805.05155), Section 2.1.2, Lemma 3.1 and Remark 3.1. The latter explicitly states that both trapped tails have zero Liouville measure.
