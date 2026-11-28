# Explicit compactness and regular-maximizer lemmas

This note isolates the existence input in `SYMMETRY.md`. The first
lemma proves attainment of the reduced maximin **value** in a compact
achronal class. The second converts a smooth spacelike maximizing slice into
an extremal minimizing cut without a Jacobi nondegeneracy assumption. Neither
lemma proves that every achronal maximizer can be replaced by a smooth
spacelike one in every singular spacetime.

## 1. A compact weighted quotient problem

Let S be a compact smooth two-dimensional spatial domain with the regulated
boundary anchors fixed. Work in a globally hyperbolic slab with an orthogonal
temporal splitting

\[
 g=-N(t,x)^2dt^2+h_t(x),\qquad t_-(x)\le t\le t_+(x),
\]

where the slab is compact, N>0, h_t is positive definite, and the fiber-volume
weight λ>0 is smooth. All these data and their needed derivatives are bounded,
and N, λ, and the smallest eigenvalue of h_t are bounded away from zero.
The two barriers and the fixed boundary trace are chosen so that the
admissible class is nonempty and its members are complete achronal graphs
over S. Any fixed boundary collars used to control the anchor asymptotics
are included in this class.

Let F be the class of continuous achronal graph functions f with
t_-≤f≤t_+ and the prescribed boundary/collar trace. Include null graphs in
the closure. For any fixed projected admissible cut C⊂S, define

\[
 L_C(f)=\int_C\lambda(f(x),x)
 \sqrt{h_{f(x)}(\dot C,\dot C)-N(f(x),x)^2
              ((f\circ C)')^2}\,ds.
 \tag{1}
\]

C may be a finite union of arcs and closed curves with the required relative
homology. Its parametrization is immaterial. Set

\[
 L(f)=\inf_C L_C(f).
 \tag{2}
\]

**Compactness lemma.** F is compact in the uniform topology, L is upper
semicontinuous, and therefore \(\max_{f\in F}L(f)\) is attained.

**Proof.** The barriers give uniform boundedness. On the compact slab the
lightcones have a uniform nonzero opening, so all achronal graphs have a
common spatial Lipschitz bound. For example, after fixing a reference
Riemannian metric h₀ on S, the spacelike/achronal differential inequality gives
\(|df|_{h_0}\le C\) almost everywhere, with C controlled by the maximum of h_t
and minimum of N. A uniform limit of achronal graphs is achronal because the
timelike relation is open. The barrier and trace conditions are closed.
Arzelà–Ascoli therefore gives compactness.

Fix a Lipschitz projected cut C, parameterized by h₀ arclength on each
component. If f_n→f uniformly, then f_n∘C→f∘C uniformly and their derivatives
converge weak-* in L∞ after taking a subsequence. The integrand in (1) is
concave in the time derivative:

\[
 v\longmapsto\lambda\sqrt{a-N^2v^2},\qquad
 |v|\le\sqrt a/N.
\]

Its coefficients converge uniformly with f_n. Upper semicontinuity of the
integral of this concave integrand gives

\[
 \limsup_n L_C(f_n)\le L_C(f).
\]

One can derive this directly from the tangent-line representation of a
concave function, or by weak upper semicontinuity after an arbitrarily small
interior regularization of the square root. Its continuity at null slopes
then permits removal of that regularization. Infima of upper-semicontinuous
functions are upper semicontinuous, so

\[
 \limsup_n L(f_n)\le L(f).
\]

The maximum is attained by the extreme-value theorem. ∎

### What this does and does not prove

The conclusion does not require an interchange of max and min. It is simply
compactness followed by upper semicontinuity. The positivity and smoothness
of λ make the weighted argument as well behaved as the unweighted one.

On a smooth strictly spacelike graph, the induced weighted metric is
uniformly Riemannian. The spatial minimum exists by the ordinary
finite-perimeter or shortest-curve compactness argument. On a completely
general limiting null graph, the induced form is degenerate. The compactness
lemma alone does not prove that its minimizing cut is smooth or even that a
classical curve realizes the infimum. Nor does it prove that the maximizing
graph is strictly spacelike. These are precisely the regularity issues to be
separated from the symmetry argument.

The requirement of a compact slab must be obtained from actual causal or
area barriers. Putting in arbitrary time or radial walls and then ignoring
their variational constraints would not establish an HRT theorem for the
unregulated spacetime.

## 2. Extremality at a regular maximizing slice

The following gives a useful stronger replacement for an abstract stable
maximin hypothesis in the regular case.

**Regular-maximizer lemma.** Suppose the maximum in (2) is attained on a
smooth spacelike graph f₀. Suppose the following standard spatial-regularity
conditions hold:

1. Minimizing cuts have the prescribed finite set of distinct anchors and
   do not meet any other constraint boundary; suitable fixed anchor collars
   have been supplied. Any frozen collars are already extremal, or are
   removed in a justified limiting argument before claiming full
   extremality.
2. The family M of minimizing cuts is compact in C¹ on interior compact sets,
   including convergence of their lengths and quadratic tangent integrals.
   Minimizers for sufficiently small smooth metric perturbations have the
   same subsequential compactness to members of M.
3. Near the union of minimizing cuts, both signs of compactly supported
   normal-time variations of the slice remain admissible. In particular,
   the auxiliary slab barriers are not active constraints there.

Condition 3 applies at every interior point of a minimizing cut. If an open
boundary collar is frozen as an auxiliary regulator, one must either already
know extremality in that collar or remove it in the limiting argument; a
frozen nonextremal collar would only give a constrained variational problem.

Then at least one minimizing cut on f₀ is fully extremal for the weighted
Lorentzian metric \(\widehat g=\lambda^2g\).

In a compact regular two-dimensional Riemannian problem with barrier-controlled
endpoints, condition 2 follows from the geodesic equation, bounded total
length, and a positive injectivity radius. Only actual minimizing curves in
the nearby smooth metrics need this C¹ compactness. Arbitrary wiggly
near-minimizing competitors need not have it and are not being assumed to.

**Proof.** Use Gaussian normal coordinates for \(\widehat g\) near the
maximizing slice, so that

\[
 \widehat g=-d\tau^2+q_\tau,\qquad \tau=0\text{ on }f_0.
\]

Every C∈M is a geodesic cut for q₀. Two distinct minimizing cuts cannot cross
transversely: union/intersection of their homology regions gives two
admissible cuts with the same total length, and rounding the resulting
corners strictly lowers that total. Tangential contact of geodesics forces
local coincidence, by uniqueness of the geodesic initial-value problem.
Thus, away from the fixed anchors, the union Λ of all minimizing components
is a geodesic lamination. It has a well-defined continuous unoriented tangent
line field. The continuity follows also by taking a convergent subsequence
of geodesic tangent vectors; a different limiting tangent would force a
crossing.

For a q₀-unit tangent T to this lamination, define the temporal first-variation
density

\[
 b(x)=\frac12\partial_\tau q_\tau(T,T)\big|_{\tau=0}.
\]

b is continuous on Λ away from the anchors and does not depend on the
orientation of T. A normal slice deformation τ=εv(x) changes the length of a
fixed cut C to first order by

\[
 \frac d{d\epsilon}L_C(\epsilon)\bigg|_0
 =\int_Cv\,b\,ds_{q_0}.
 \tag{3}
\]

Assume for contradiction that no C∈M is fully extremal. Since its spatial
first variation already vanishes, this means that b is not identically zero
on any C. Hence

\[
 Q(C)=\int_Cb^2\,ds_{q_0}>0\qquad(C\in M).
\]

Compactness and continuity give a strictly positive uniform lower bound on
Q. The same remains true with a sufficiently small fixed anchor collar
removed: choose a cutoff χ supported away from the boundary, using
compactness of M and monotone exhaustion of each curve's interior, such that

\[
 \inf_{C\in M}\int_C\chi b^2\,ds_{q_0}>0.
\]

The continuous function χb on the closed interior portion of Λ extends to a
continuous function on the slice and can be approximated uniformly there by
a smooth compactly supported v. Choosing the approximation sufficiently
accurate yields a common positive first-order improvement:

\[
 \int_Cvb\,ds_{q_0}\ge\delta>0\qquad(C\in M).
 \tag{4}
\]

The deformed slice remains spacelike and admissible for small ε, by condition
3. Let C_ε be an actual minimizing cut in its smooth induced metric q_ε.
Those metrics converge smoothly to q₀ and remain uniformly elliptic. Along
every sequence ε→0, their minimizing cuts have a subsequence converging to
some C₀∈M. The metric Taylor expansion gives, uniformly for these curves,

\[
 L_{q_\epsilon}(C_\epsilon)
 =L_{q_0}(C_\epsilon)+\epsilon B(C_\epsilon)+O(\epsilon^2),
 \quad B(C)=\int_Cvb\,ds_{q_0},
\]

where b for a general tangent is the same smooth quadratic form
\(\tfrac12\partial_\tau q_\tau\). By compactness,
B(C_ε)→B(C₀)≥δ along such subsequences. Since
\(L_{q_0}(C_\epsilon)\ge L(f_0)\), for all sufficiently small ε>0,

\[
 L(f_\epsilon)\ge L(f_0)+\tfrac12\epsilon\delta,
\]

contradicting maximality. Therefore some C∈M has b≡0. Its spatial and temporal
first variations both vanish, so it is an extremal curve of \(\widehat g\).
∎

This lemma is dimensionally special only in its convenient geodesic
lamination regularity: after reduction, the spatial cuts are curves. It
allows degeneracies and continuous families. It does not prove regularity
of a null or nonsmooth maximizing slice.

Combining this lemma with the BV/Haar spatial-minimality lemma and the full
null-representative comparison gives the HRT symmetry theorem directly.
There is no need to show that every spatial minimizing cut on the maximizing
slice is extremal.

## 3. A concrete singularity barrier for planar black branes

Consider a planar Schwarzschild–AdS_{d+1} region with

\[
 ds^2=-F(r)dt_s^2+\frac{dr^2}{F(r)}
 +r^2(dx^2+dy_1^2+\cdots+dy_{d-2}^2),
 \qquad F(r)\sim-\mu r^{-(d-2)}\quad(r\to0).
\]

Inside the horizon, proper time to the spacelike singularity satisfies
\(\tau\propto r^{d/2}\). The Kasner exponents are

\[
 p_{t_s}=-\frac{d-2}{d},\qquad
 p_x=p_{y_i}=\frac2d.
\]

Since \(\lambda=r^{d-2}\), the reduced effective spatial metric at constant τ has
leading behavior

\[
 \widehat q_{t_st_s}\propto\tau^{2(d-2)/d},\qquad
 \widehat q_{xx}\propto\tau^{4(d-1)/d}.
 \tag{5}
\]

Both powers are strictly positive for d≥3. In the exact stationary metric,
both coefficients increase away from the singularity in a sufficiently
small collar. The same conclusion holds for a general Kasner collar with
uniform C¹ asymptotic control that preserves these derivative signs; leading
power laws without derivative control are insufficient. A
nonconstant slice τ=f(t_s,x) also has a nonpositive tilt contribution to its
induced metric. Replacing its singular portion by

\[
 f\longmapsto\max(f,\tau_0)
\]

raises the effective spatial coefficients and removes the negative tilt
term on the clamped portion. It therefore cannot decrease any reduced cut
length. With a suitable collar extension this gives an achronal barrier
separated from the singularity, exactly the invariant analogue of Wall's
Kasner argument.

For a smooth spacetime whose relevant singularities have such uniform
Schwarzschild-like Kasner collars, this supplies a local interior barrier.
The remaining outer and time barriers and regularity hypotheses must still
be established before applying the compact-slab lemma. It is not necessary
to assume an invariant HRT surface
to construct this barrier: it is a statement about the entire reduced
maximin slice class.

For Vaidya geometries, a final stationary black-brane region has this collar.
To claim a global theorem through a quench, one must also control any
time-dependent singularity collar and, for an ideal thin shell, the
shell/singularity junction. A smooth positive-energy shell can be used where
the ordinary smooth maximin assumptions hold; passing to a distributional
thin shell requires the usual junction or limiting argument. The stationary
collar alone should not be described as checking all these global issues.

## 4. Relation to the existing maximin existence proof

The compactness mechanism is the quotient/weighted version of §3.4 in
Aron Wall, *Maximin Surfaces, and the Strong Subadditivity of the Covariant
Holographic Entanglement Entropy*,
[arXiv:1211.3494](https://arxiv.org/abs/1211.3494). The singularity clamp is
the specialization of Theorem 11. Wall explicitly separates these existence
arguments from the regularity/stability issue for general achronal maximizing
slices in §3.5.

Accordingly, the correct clean theorem is one of the following:

- Under the ordinary stable reduced-maximin existence and regularity
  hypotheses, HRT symmetry follows from BV averaging and representatives.
- Under the explicit compact-slab and regular-spacelike-maximizer hypotheses
  above, compactness proves attainment and the lamination lemma proves the
  extremal minimizing cut; the same HRT symmetry conclusion follows.

Neither formulation makes full HRT existence alone imply reduced-maximin
existence, and neither hides a new assumption that HRT itself is symmetric.
