# Exact branch admissibility in the symmetric two-particle collision

28 November 2026. Model: Aref'eva--Khramtsov--Tikhanovskaya,
[arXiv:1706.07390v3](https://arxiv.org/abs/1706.07390), Sections 2.2--3.3.
The statements in Sections 1--5 concern individual equal-time interval
branches. They do not establish simultaneous physical admissibility of
an arbitrary multicomponent matching.

## 1. Notation and developed geodesics

Let \(R>0\), \(\mu=\pi R\), \(E=\coth(\mu/2)\). Consider a shortest
boundary arc crossing particle 1,

\[
 -\pi<a<0<b<\pi,\quad b-a\le\pi,
 \quad c=\frac{R(a+b)}2,\quad d=\frac{R(b-a)}2.
\]

Thus \(|c|<d\le\mu/2\). Both endpoints have the same BTZ time \(t\ge0\).
Write \(C=\cosh Rt\), \(S=\sinh Rt\), and use ambient signature
\((-++-)\). A null boundary representative is

\[
 p(x)=(-\cosh Rx,C,\sinh Rx,S).
\]

The particle holonomy is \(u_1=\bigl(\begin{smallmatrix}1&0\\-2E&1
\end{smallmatrix}\bigr)\). Its \(\#\)-action on an ambient vector is

\[
 \begin{aligned}
 x_0'&=x_0+2Ex_2+2E^2(x_0+x_1),\\
 x_1'&=x_1-2Ex_2-2E^2(x_0+x_1),\\
 x_2'&=x_2+2E(x_0+x_1),\qquad x_3'=x_3.
 \end{aligned}
\]

The crossing curve is recovered from the developed geodesic
\(G(p(a),p(b)^\#)\): keep its left segment up to the wedge face, and map
the other segment by the inverse holonomy. This is precisely the image
construction of the source's equations (3.12)--(3.21).

## 2. Direct appearance and crossing disappearance without a triangle inequality

Define

\[
 P(C)=(C-\cosh(c-d))(C-\cosh(c+d)),\qquad
 Q(C)=C\cosh c-\cosh d,
\]

and

\[
 F(C)=E P(C)+\sinh d\,Q(C).
\]

The direct branch crosses the symmetry plane \(x_2=0\). At that point,
\(x_0+x_1\) has the sign of \(Q\). The regular exterior of the wedge
at that plane has \(x_0+x_1\ge0\). Consequently the direct branch first
becomes available at

\[
 C_{\rm cr}=\frac{\cosh d}{\cosh c},\qquad Q(C_{\rm cr})=0.
\]

This is the source's equation (3.30), obtained directly from the cut faces.

For the crossing branch, set
\(u_a=C-\cosh(c-d)\), \(u_b=C-\cosh(c+d)\), and

\[
 w_a=\sinh(c-d)-Eu_a<0,\qquad
 w_q=\sinh(c+d)+Eu_b>0.
\]

These inequalities hold for \(t\ge0\) and endpoints in the stated half
circles. The left wedge plane is \(w=x_2-E(x_0+x_1)=0\). An unnormalized
intersection point of the image geodesic with it is

\[
 Y=w_q p(a)-w_a p(b)^\#.
\]

A direct calculation gives

\[
 Y_0+Y_1
 =u_a\sinh(c+d)-u_b\sinh(c-d)+2Eu_au_b
 =2F(C).
\]

The proper wedge face lies on the side \(x_0+x_1\le0\); its endpoint
on particle 1 has \(x_0+x_1=x_2=0\). Hence the crossing branch reaches
the particle precisely when

\[
 F(C_o)=0.
\]

At this point its normalized intersection is on the fixed null line
\(X=(-s,s,0,1)\). This criterion comes from the actual geodesic embedding,
not a comparison of plotted lengths.

There is exactly one root \(C_o>1\). To see this, put
\(H=\cosh d-\cosh c>0\). Then

\[
 F(1)=H(EH-\sinh d)<0,
\]

because \(H\le\cosh d-1=\sinh d\tanh(d/2)\) and
\(E\tanh(d/2)<1\). The quadratic \(F\) has positive leading
coefficient, so precisely one of its two roots is above one.

By reflecting the setup assume \(c\ge0\). Let
\(m=\cosh(d-c)\), \(M=\cosh(d+c)\). Then

\[
 m\le C_{\rm cr}\le M,
 \qquad F(C_{\rm cr})=E P(C_{\rm cr})\le0,
 \qquad F(M)=\sinh d\,Q(M)\ge0.
\]

Thus

\[
 \boxed{C_{\rm cr}\le C_o\le M.}
\]

Equality \(C_o=C_{\rm cr}\) occurs for the symmetric interval \(c=0\).
For a strictly asymmetric interval both branch coexistence inequalities
are strict.

## 3. Exact length comparison at both physical thresholds

Use compatible radial cutoffs and omit their common additive divergence.
The exponential of the direct renormalized length is

\[
 K_{\rm dir}=4\sinh^2d.
\]

Expanding \(-2\langle p(a),p(b)^\#\rangle\), or equivalently the
source's equation (3.40), gives

\[
 K_{\rm cross}
 =4\left[\sinh^2d+E^2P(C)+2E\sinh d\,Q(C)\right].
\]

At direct appearance,

\[
 K_{\rm cross}(C_{\rm cr})-K_{\rm dir}
 =4E^2P(C_{\rm cr})\le0.
\]

At crossing disappearance, use \(EP=-\sinh d\,Q\) to obtain

\[
 \boxed{K_{\rm cross}(C_o)-K_{\rm dir}
 =4E\sinh d\,Q(C_o)\ge0.}
\]

This proves the comparison claimed in equation (3.39) of the source without
using its purported Lorentzian triangle inequality. A general triangle
inequality for spacelike sides would not justify that step.

Moreover \(K_{\rm cross}\) is strictly increasing for \(C\ge1\):

\[
 \frac{dK_{\rm cross}}{dC}
 =8E\left[E(C-\cosh c\cosh d)+\sinh d\cosh c\right]>0.
\]

For a convenient lower bound at \(C=1\), use \(\cosh c\le\cosh d\)
and \(E\tanh d\le1\); the derivative increases with \(C\). Strictness
follows from \(|c|<d\). Therefore there is exactly one length switch, and
it lies in the interval where both branches are available. The formal
minimum of these two length expressions never selects the branch before
its physical appearance or after its disappearance.

This uses the source's classification of the relevant short-arc branches
and does not independently classify every winding extremal curve.
At a contact threshold, a limiting curve can meet the null particle
worldline. The expressions above remain continuous limiting lengths; they
are not a claim that the curve at the contact instant belongs to a
regular-locus-only candidate space. A prescription restricted to regular
curves must distinguish its infimum from attainment there.

## 4. Independent audit of the positive interval kernel

Put \(A=e^{Ra}<1<B=e^{Rb}\), \(T=e^{Rt}\). The exact ratio factors are

\[
 R_1=1+\frac{E(A-T)(B-T)}{T(B-A)},\qquad
 R_2=1+\frac{E(AT-1)(BT-1)}{T(B-A)},
\]

with \(L_{\rm cross}-L_{\rm dir}=h=\log R_1+\log R_2\). Both factors
are positive: their common lower bound is
\(1-E\tanh(R(b-a)/4)>0\).

The exact endpoint derivatives are

\[
 \begin{array}{ll}
 (R_1)_A=\dfrac{E(B-T)^2}{T(B-A)^2},&
 (R_1)_B=-\dfrac{E(A-T)^2}{T(B-A)^2},\\[6pt]
 (R_2)_A=\dfrac{E(BT-1)^2}{T(B-A)^2},&
 (R_2)_B=-\dfrac{E(AT-1)^2}{T(B-A)^2}.
 \end{array}
\]

Thus \(h_a\ge0\), \(h_b\le0\). Writing the two crossing factors as
bilinear numerators, each logarithm has positive mixed derivative because
the numerator's \(AB\)-determinant is exactly \(T^2\).

For switch regularity, one need not assume a nonsingular \(h=0\) wall.
Let \(f\) be a smooth concave approximation to \(\min(0,h)\), with
\(0\le f'\le1\), \(f''\le0\). Then

\[
 \partial_a\partial_b[L_{\rm dir}+f(h)]
 =(1-f')\partial_a\partial_bL_{\rm dir}
 +f'\partial_a\partial_bL_{\rm cross}
 +f''h_ah_b\ge0.
\]

Distributional convergence proves positivity of the minimum kernel even
at degenerate switch points. Where the switch wall is regular, the extra
term is the positive measure \(-h_ah_b\,\delta(h)\).

The resulting phase ordering is also exact: a nested outer interval
cannot be in the direct phase while its inner interval is in the crossing
phase.

## 5. Why the Fuchsian two-cusp vacuum domain is insufficient

The special symmetric model has diagonal holonomies, so it is natural to
try to embed it into the vacuum domain of a Fuchsian surface with two
cusps and one funnel. This cannot cover the full post-quench boundary.

The parabolic limit point of particle 1 is
\(P=(-1,1,0,0)\), at global time \(-\pi/2\) and angle zero. The physical
boundary point at BTZ time \(t>0\), angle zero, is

\[
 Q=(-1,\cosh Rt,0,\sinh Rt).
\]

It is timelike to \(P\): the global spatial angle is unchanged and global
time has advanced. Equivalently
\(\langle P,Q\rangle=\cosh Rt-1>0\). Since \(P\) belongs to the
diagonal group's limit set, \(Q\) is outside its invisible domain. The
smooth vacuum cusp theorem therefore does not identify the generic
post-quench HRT family.

## 6. Status of simultaneous realization and the Czech--Dong result

[Czech--Dong, arXiv:1905.03787](https://arxiv.org/abs/1905.03787),
Section 2.1, organizes HRT phases as nonalternating colorings; Appendix A
uses SSA to show that spatial ordering minimizes the corresponding sums
of interval entropies. It does not prove an acausal-slice realization lemma
for independently chosen interval minimizers in a singular point-particle
geometry. Its known contraction-inequality result, including all the
five-party inequalities, must be credited; checking those inequalities
numerically here is not a first theorem of their covariant validity.

The actual missing physical implication remains

\[
 \begin{gathered}
 \text{optimal noncrossing matching of the exact interval kernel}\\
 \Downarrow\\
 \text{one embedded extremal multicurve}\\
 \text{on an acausal homology slice}.
 \end{gathered}
\]

SSA and phase ordering alone are algebraic statements and do not prove
that implication. The exact image-geodesic description above reduces it
to a concrete finite wedge/causality problem in this special model, rather
than the unrestricted BBS collision-surgery problem.

## 7. A reflection-restricted route, preserved as an open proposition

There is extra structure specific to the symmetric collision which is
absent in general BBS collision surgery. The gluing map on each paired
\(W\) or \(V\) face preserves \(x_0,x_1,x_3\) and sends
\(x_2\mapsto-x_2\). The physical spacetime therefore has a spatial
reflection exchanging two copies of a right half-domain. That half-domain
is bounded by three timelike totally geodesic planes with inward unit
spacelike normals

\[
 n_0=(0,0,1,0),\qquad
 n_1=(-E,E,1,0),\qquad
 n_2=(\sinh\mu,0,-\cosh\mu,0).
\]

Its inequalities are \(\langle X,n_j\rangle\ge0\). The folding map
reflects the negative-\(x_2\) half into the positive half. Piecewise it is
a time-orientation-preserving AdS isometry; along identification faces it
is continuous. Hence it sends causal curves to causal curves. If a
proposed multicurve's folded image is acausal in the ambient AdS patch,
its physical image is acausal as well. The converse need not hold.

For boundary points of the half-circle write
\(p(u)=(-\cosh u,C,\sinh u,S)\), \(0<u<\mu\). The three positive wall
weights are

\[
 w_0=\sinh u,\quad
 w_1=\sinh u+E(C-\cosh u),\quad
 w_2=\sinh(\mu-u).
\]

A reflection-symmetric full geodesic folds to an orthogeodesic from
\(p(u)\) to one of these walls. Its unfolded endpoints are
\(p(u)\) and \(J_jp(u)\), where
\(J_jX=X-2\langle X,n_j\rangle n_j\). Its full length factor is
\(4w_j^2\). The minimum wall is therefore

- \(H_0\) for \(u<Rt\),
- \(H_1\) for \(Rt<u<\mu-Rt\),
- \(H_2\) for \(u>\mu-Rt\),

when \(Rt<\mu/2\); after that time only the direct \(H_0,H_2\) branches
remain. Equalities are particle-contact limits.

For reflection-invariant boundary labels, the formal positive graph has
an invariant minimum cut by submodularity. An invariant noncrossing
matching cannot contain an asymmetric opposite-half component together
with its reflected distinct component: their endpoints alternate. Thus
opposite-half components in such a matching are individually symmetric,
while same-half components occur in mirror pairs.

This makes a physical theorem for reflection-invariant parties plausible.
Several local acausality comparisons become transparent. For example,
with \(J=J_1\),

\[
 \langle p(u),Jp(v)\rangle
 =\langle p(u),p(v)\rangle-2w_1(u)w_1(v)<0.
\]

Thus the half-boundary arc and its reflected image form compatible
acausal endpoint sets. Ordered orthogeodesics to the same wall are
noncrossing chords of that set. A same-half direct chord whose endpoints
are on the same side of the orthogeodesic's boundary endpoint is likewise
noncrossing and acausal with it in the developed patch.

**Unfinished proposition.** Complete all mixed-wall comparisons, including
the \(H_1\)--\(H_0,H_2\) pairs and contact limits, and show that the
folded family extends to a slice whose double is a physical Cauchy
homology surface. These steps have not been packaged into a checked proof
here. This section is a research route, not an additional established
collision entropy-cone theorem. For arbitrary, non-reflection-invariant
parties, folding can create extra interleavings and cannot be assumed to
preserve acausality of the selected family.
