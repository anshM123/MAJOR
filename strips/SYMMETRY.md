# HRT symmetry reduction for translational strip regions

## Status and precise conclusion

The symmetry step can be proved without assuming uniqueness of HRT surfaces,
without excluding symmetry-breaking competitors, and without a Jacobi
nondegeneracy assumption. The mechanism is:

1. Invariant cuts on a fixed invariant Cauchy slice attain the **unrestricted**
   spatial minimum. This follows by averaging homology-region indicators and
   applying BV coarea.
2. A stable reduced maximin curve lifts to a fully extremal surface.
3. Every competing full-dimensional extremal surface has a null representative
   on that invariant slice with no larger area. Spatial minimality then proves
   that the invariant surface is a full HRT surface.

There is also a useful correction to the proposed argument: the effective
three-dimensional metric **does obey null convergence**, provided the full
metric does and the translational fibers are spacelike and nondegenerate. This
holds with anisotropic fiber metrics and off-diagonal Kaluza–Klein connections.

The theorem below retains the usual existence and regularity inputs of
maximin. It does not derive stable maximin existence in every singular NEC
spacetime, and it does not identify full HRT existence alone with the stronger
input needed for the reduced maximin construction.

## 1. Geometric setup

Let the bulk dimension be D=m+3, with \(m=d-2\)≥1. Initially take a free spacelike
\(T^m\) action, with coordinate volume V, preserving the bulk metric and the
boundary entangling region. Locally write

\[
 ds_D^2=g_{ab}(x)dx^a dx^b+
 G_{ij}(x)(dy^i+A_a^i(x)dx^a)(dy^j+A_b^j(x)dx^b),
 \tag{1}
\]

where g is Lorentzian on the three-dimensional orbit space B and G is positive
definite. The orbit directions need not have equal warp factors, and A need not
vanish. The local calculations also apply to a principal torus bundle; the
global strip examples ordinarily use the product bundle.

Set

\[
 \lambda=(\det G)^{1/2},\qquad u=\log\lambda,\qquad
 \widehat g=\lambda^2g=(\det G)g.
 \tag{2}
\]

An allowed boundary region A is a union of complete fibers over intervals on a
boundary Cauchy curve. Its entangling planes and its domain of dependence must
be invariant. An arbitrary boundary Cauchy slice is not sufficient: if it is
given by t=f(x,y) and f(a,y) varies along an edge x=a, the anchor itself breaks
translations and no invariant bulk surface can have that anchor.

The local symmetry theorem does not require a disc orbit topology. Disc
topology, or another topology covered by the graph representation theorem, is
a separate input to the proposed entropy-cone application. Translation
invariance alone does not impose it.

### Cauchy slices

For a base spacelike Cauchy slice σ, its lift Σ=π⁻¹σ is spacelike by (1). Under
the usual global assumptions it is also Cauchy. A full causal vector projects
to a nonzero base causal vector because G>0. A full inextendible causal curve
cannot project to a curve with an ordinary endpoint: in a compact base chart a
local temporal coordinate bounds both base spatial velocities and the fiber
velocity dy+A dx, so the full curve also has an endpoint if its projection
does.

For a compact torus action, invariant Cauchy slices exist whenever a smooth
Cauchy temporal function exists. Average such a function τ over normalized
Haar measure. Future timelike covectors form a convex cone, so the averaged
function is temporal. Along any inextendible causal curve, τ(g·c(s)) tends to
the appropriate infinity for every g. Monotonicity in s and compactness of the
torus make this divergence uniform in g, so the average remains Cauchy. This
also provides the standard Cauchy structure on the quotient. **It does not
assert that averaging a maximizing slice preserves its maximin value.**

For a noncompact \(\mathbb R^m\) action, an appropriate causal orbit space must be assumed
or established separately; spacelike orbits alone do not prove this global
property.

## 2. Exact area reduction

Let γ be a spacelike curve in B, with arbitrary parameter s. The induced metric
of its invariant lift has block form

\[
 h=\begin{pmatrix}
 g(\dot\gamma,\dot\gamma)+G(A\dot\gamma,A\dot\gamma)
 & (GA\dot\gamma)^T\\
 GA\dot\gamma&G
 \end{pmatrix}.
\]

The Schur complement gives

\[
 \det h=\det G\,g(\dot\gamma,\dot\gamma).
\]

Consequently,

\[
 \boxed{\frac{\operatorname{Area}(\pi^{-1}\gamma)}{V}
 =\int_\gamma\lambda\,ds_g
 =\operatorname{Length}_{\widehat g}(\gamma).}
 \tag{3}
\]

This identity is valid even when the full metric has time-fiber and
space-fiber cross terms. It applies componentwise to disconnected cuts and
closed components.

For \(G=e^{2\phi}I_m\) it becomes the proposed formula
\(\widehat g=e^{2mφ}g\).

## 3. The effective metric obeys null convergence

Define

\[
 P_a=G^{-1}\partial_aG,\qquad F^i=dA^i.
\]

For a g-null vector k, its horizontal lift

\[
 K=k^a(\partial_a-A_a^i\partial_{y^i})
\]

is a full null vector. The horizontal Ricci identity for (1) is

\[
 \operatorname{Ric}_D(K,K)
 =\operatorname{Ric}_g(k,k)
 -k^ak^b\nabla_a\nabla_bu
 -\frac14\operatorname{Tr}(P_k^2)
 -\frac12G_{ij}(\iota_kF^i)_a(\iota_kF^j)^a,
 \tag{4}
\]

with base indices raised using g. To see the scalar-matrix terms directly,
temporarily put A=0. The additional mixed Christoffels are

\[
 \Gamma^a_{ij}=-\tfrac12g^{ab}\partial_bG_{ij},\qquad
 \Gamma^i_{aj}=\tfrac12(G^{-1}\partial_aG)^i{}_j.
\]

Their contribution to horizontal Ricci is

\[
 -\tfrac12\operatorname{Tr}(G^{-1}\nabla_a\partial_bG)
 +\tfrac14\operatorname{Tr}(P_aP_b)
 =-\nabla_a\nabla_bu-\tfrac14\operatorname{Tr}(P_aP_b).
\]

For nonzero A, use the horizontal frame
\(E_a=\partial_a-A_a^i\partial_i\). Its bracket
\([E_a,E_b]=-F^i_{ab}\partial_i\) supplies the additional horizontal Ricci
term \(-\tfrac12G_{ij}F^i_{ac}F^j_b{}^c\), giving (4).

In base dimension three, the conformal Ricci formula for \(\widehat g=e^{2u}g\)
gives on null vectors

\[
 \operatorname{Ric}_{\widehat g}(k,k)
 =\operatorname{Ric}_g(k,k)-k^ak^b\nabla_a\nabla_bu+(ku)^2.
\]

Combining these equations yields the key identity

\[
 \boxed{\operatorname{Ric}_{\widehat g}(k,k)
 =\operatorname{Ric}_D(K,K)+(ku)^2
 +\frac14\operatorname{Tr}(P_k^2)
 +\frac12G_{ij}(\iota_kF^i)_a(\iota_kF^j)^a.}
 \tag{5}
\]

Each added term is nonnegative:

- \(P_k\) is similar to the real symmetric matrix
  \(G^{-1/2}(\partial_kG)G^{-1/2}\), so its trace square is nonnegative.
- Each one-form \(\iota_kF^i\) annihilates k. The inverse Lorentzian form is
  positive semidefinite on that null-annihilator subspace. Contracting the
  resulting positive-semidefinite matrix with G>0 gives a nonnegative number.

Thus full NCC implies effective NCC. In the isotropic, untwisted case,

\[
 \boxed{\operatorname{Ric}_{\widehat g}(k,k)
 =\operatorname{Ric}_D(k,k)+m(m+1)(kφ)^2\ge0.}
 \tag{6}
\]

The geometric assumption is NCC, \(\operatorname{Ric}_D(K,K)\ge0\). Matter NEC
implies it through the classical Einstein equations, because both the
cosmological and trace terms vanish upon null contraction. Matter NEC by
itself is not a replacement for NCC in arbitrary higher-derivative gravity.

The single-circle Einstein-frame version of this inheritance is explicitly
discussed in Parikh–van der Schaar, *Derivation of the Null Energy Condition*,
§IV, eq. (23), [arXiv:1406.5163](https://arxiv.org/abs/1406.5163). Equation (5)
also includes the varying fiber shape and KK gauge fields.

An independent exact check is recorded in `check_kk_ricci.py`: it constructs
full coordinate Ricci from metric two-jets for m=1 and m=2 examples with
nonzero gauge fields and varying fiber metrics. Both the complete horizontal
Ricci identity and its null-conformal consequence have exactly zero residual.
This is a check of signs and coefficients, not a substitute for the derivation.

### Asymptotic qualification

For an AdS_D Poincaré asymptotic region,
\(\widehat g\sim z^{-2(m+1)}(-dt^2+dx^2+dz^2)\), up to powers of the AdS radius.
It is generally **not asymptotically AdS₃**. Its causal structure agrees with
the quotient's, and it obeys NCC, but an AdS₃ theorem whose hypotheses concern
the precise AdS boundary expansion cannot be imported verbatim. The relevant
length functional is the original area density and inherits its strip
divergences.

## 4. The spatial symmetry lemma

**Lemma.** On an invariant spacelike slice Σ=π⁻¹σ with invariant boundary trace
A, the unrestricted infimum of the areas of cuts homologous to A equals V
times the weighted-length infimum of cuts on σ. In particular, an invariant
minimizer exists whenever the reduced minimum is attained.

**Proof.** Use a compact spatial regulator, or the usual limiting
finite-perimeter formulation, and let R⊂Σ be any admissible homology region.
Its boundary has trace A on the asymptotic/regulator boundary. Write its
relative perimeter as \(P(R)=|D\chi_R|(\Sigma)\), excluding the fixed boundary
piece. Average its indicator over the torus:

\[
 f=\int_{T^m}\chi_{gR}\,d\mu(g).
\]

Because the torus acts isometrically and total variation is convex,

\[
 |Df|(\Sigma)\le\int_{T^m}|D\chi_{gR}|(\Sigma)\,d\mu(g)=P(R).
 \tag{7}
\]

The averaged function is invariant, lies in [0,1], and has the same boundary
trace χ_A. Writing it as f(x), the inverse-metric and volume formulas from
(1) give

\[
 |Df|(\Sigma)=V\int_\sigma\lambda\,|df|_{h_\sigma}
 \,d\operatorname{vol}_{h_\sigma}.
 \tag{8}
\]

Weighted coarea gives

\[
 \int_\sigma\lambda|df|\,d\operatorname{vol}
 =\int_0^1P_\lambda(\{f>s\})\,ds.
 \tag{9}
\]

For almost every s∈(0,1), the set {f>s} has the same prescribed boundary
trace A, hence its boundary is an admissible reduced cut. Let L_σ(A) be the
minimum weighted length of such a cut. Equations (7)–(9) yield

\[
 \boxed{P(R)\ge V L_\sigma(A).}\tag{10}
\]

Conversely, lift a reduced minimizing region. Equation (3) gives equality.
This proves the lemma. ∎

The use of regions is essential: HRT homology means that a cut and A bound
such a region. The averaging does not impose an extra restriction on allowed
relative-homology cuts. It includes disconnected competitors, winding
competitors, and nonsymmetric competitors. It also avoids any regularity
assumption on those competitors beyond finite perimeter.

This proof does not rely on a fixed-y slicing inequality, which would need
additional care in the presence of the connection A. It proves equality of
the full spatial minima directly.

## 5. Extremality of invariant lifts

**Lemma.** If γ is an extremal curve for \(\widehat g\), then \(X=\pi^{-1}\gamma\) is extremal
for full-dimensional area.

**Proof.** The mean-curvature vector H_X is invariant under the torus because
X and the ambient metric are invariant. It is horizontal: vertical vectors
are tangent to X, and the normal bundle is orthogonal to the entire fiber.
Invariant normal variations are therefore arbitrary normal variations of γ,
lifted horizontally. Equation (3) says that their first area variation
vanishes. If H_X were nonzero at any point, an invariant normal variation
supported over a sufficiently small segment of γ would detect it, contrary
to this vanishing. Hence H_X=0 pointwise. The full first variation then
vanishes for all normal variations, including y-dependent ones. ∎

Equivalently, compact-group averaging of any test variation leaves its first
variation unchanged, because H_X and the induced measure are invariant.

## 6. Symmetry theorem for HRT

**Theorem (translational HRT reduction).** Let a classical asymptotically AdS
Einstein spacetime satisfying NCC admit the regular free spacelike torus
action (1), preserving the boundary region A. Assume:

1. The ordinary causal/homology hypotheses needed for the null-representative
   construction hold, with a common boundary prescription. For asymptotic
   anchors assume the regular same-anchor expansion, conformally spacelike
   containing slice, and local null-representative behavior specified in
   `UV_AVERAGING.md`; that note gives sufficient boundary conditions and
   a direct cutoff proof.
2. The reduced weighted maximin problem
   \[
   L_*(A)=\sup_{\sigma}\inf_{\gamma\sim A,\,\gamma\subset\sigma}
   \operatorname{Length}_{\widehat g}(\gamma)
   \tag{11}
   \]
   has a stable maximizer represented by a smooth spacelike curve γ_* on an
   invariantly liftable Cauchy slice σ_*; the standard local maximin
   extremality argument is applicable. A sufficient regularity regime is a
   smooth spacelike maximizing slice together with a stable minimizing cut.

Then X_*=π⁻¹γ_* is a full HRT surface, and

\[
 \boxed{\frac{S_{\rm HRT}(A)}{V}
 =\frac{1}{4G_D}\operatorname{Length}_{\widehat g}(\gamma_*).}
 \tag{12}
\]

**Proof.** Stable reduced maximin is extremal for \(\widehat g\), by the usual
two-normal-directions maximin argument. The previous lemma makes its lift
fully extremal. By the spatial symmetry lemma, X_* minimizes area against
**all** cuts on Σ_*=π⁻¹σ_*, not merely invariant cuts.

Let Y be any other full-dimensional extremal surface with the same anchors
and homology. Its representative \(\widetilde Y\) on Σ_* is obtained from the
appropriate future/past null sheets, with generators discarded at caustics.
The initial expansion is zero. Raychaudhuri and full NCC make expansion
nonpositive along these sheets, so

\[
 \operatorname{Area}(\widetilde Y)\le\operatorname{Area}(Y).
\]

The representative has the required relative homology. Therefore

\[
 \operatorname{Area}(X_*)
 \le\operatorname{Area}(\widetilde Y)
 \le\operatorname{Area}(Y).
 \tag{13}
\]

Since X_* itself is extremal, it realizes the least extremal area and is HRT.
No strict focusing, uniqueness, or exclusion of continuous families is
needed for these weak inequalities. ∎

For asymptotic anchors, the displayed area comparison is shorthand for a
finite area-difference comparison, not an inequality between infinities.
`UV_AVERAGING.md` proves
$\operatorname{Area}(X_{*,\epsilon})\le
\operatorname{Area}(Y_\epsilon)+o(V)$ directly: first add the exact BV
boundary-mismatch penalty, then control both it and the focusing cutoff
mismatch by the common extremal boundary expansion. The finite same-anchor
limit gives the required ordering. No renormalized area of the generally
nonextremal representative is needed.

The null-representative comparison is Wall's maximin/HRT mechanism; see
*Maximin Surfaces, and the Strong Subadditivity of the Covariant Holographic
Entanglement Entropy*, [arXiv:1211.3494](https://arxiv.org/abs/1211.3494),
Corollary 6(e), Theorem 15(c). The added ingredient here is the unrestricted
spatial minimality supplied by BV averaging.

### Exactly what remains an input

The proof establishes symmetry once a stable reduced maximin exists. The
class of invariant slices is nonempty, but nonemptiness is not attainment of
the maximin supremum. Averaging an arbitrary full maximizing slice is not a
proof of attainment, and no minimax interchange is used here.

Wall's original proof establishes existence in specified classes of
spacetimes and explicitly discusses a stability/regularity lacuna for
general achronal maximizing slices (§3.5). Its Theorem 12 proves a stable
choice when the maximizing slice is smooth and spacelike; Theorems 10–11
give existence under additional global hypotheses. The same compactness
strategy can be applied to the invariant slice class, which is closed under
uniform convergence, but singularity avoidance and stability still have to
be checked in the geometry class being claimed. The later restricted-maximin
work, [arXiv:1901.03879](https://arxiv.org/abs/1901.03879), also explicitly
formulates a smooth-region condition.

Thus the scientifically correct scope is “under standard stable-maximin
existence and regularity hypotheses,” not “NEC alone proves existence in
every asymptotically AdS spacetime.” This qualification already exists in
the usual covariant HRT framework.

Nor should an arbitrary hard radial cutoff be silently treated as an
asymptotic AdS boundary: Wall notes finite-distance boundary subtleties. The
spatial BV lemma is exact on compact regulators, while the covariant
representative theorem should be applied with its own valid boundary
prescription and only then compared in a common asymptotic cutoff limit.

## 7. Noncompact translations and entropy density

### Periodic definition

For a product \(\mathbb R^m\) translational solution, impose transverse periods, assuming
the quotient remains in the admissible causal class. Every period gives the
same reduced problem. The full HRT area, including all periodic nonsymmetric
competitors, is V times the same reduced length. Consequently the
torus-regulated entropy density is independent of the periods, and its
decompactification limit is exactly (12).

The torus does not artificially remove competing winding cuts: the spatial
BV lemma compared against every cut that bounds an admissible homology
region in the quotient.

### Direct Følner limit for arbitrary competitors

There is also a direct version on a noncompact invariant spatial slice.
Fix a compact regulated base, let \(Q_L=[-L,L]^m\), and average a homology-region
indicator only over Q_L:

\[
 f_L(x)=\frac1{|Q_L|}\int_{Q_L}\chi_R(x,y)\,dy.
\]

The horizontal derivatives are
\(D_a=\partial_a-A_a^i\partial_i\). Hence

\[
 \partial_af_L=\operatorname{avg}_{Q_L}(D_a\chi_R)
 +A_a^i b_{L,i},\qquad
 b_{L,i}=\frac1{|Q_L|}\int_{Q_L}\partial_i\chi_R.
\]

The face-flux term obeys \(|b_{L,i}|\le L^{-1}\) for almost every x and
almost every L, in the appropriate BV trace sense. Convexity of the
horizontal covector norm, followed by weighted coarea, gives

\[
 \frac{P(R;\sigma\times Q_L)}{|Q_L|}
 \ge L_\sigma(A)-\frac{C}{L}.
 \tag{14}
\]

C is finite on the regulated base and depends on the bounded connection and
the weighted base volume. Thus

\[
 \liminf_{L\to\infty}\frac{P(R;\sigma\times Q_L)}{|Q_L|}
 \ge L_\sigma(A).
\]

For a block-diagonal metric A=0 the face-flux error vanishes. The same
argument holds along suitable Følner sets.

To extend the **covariant** comparison to arbitrary aperiodic extremal
surfaces, one must also control transverse flux of the null representatives.
A sufficient condition is that the relevant generators remain in a fixed
regulated causal slab with bounded transverse displacement C₀. Then

\[
 \operatorname{Area}(\widetilde Y\cap Q_L)
 \le\operatorname{Area}(Y\cap Q_{L+C_0}),
\]

and existence of the entropy density plus negligible side corrections
recovers (13) after division by |Q_L|. Such a displacement bound follows on
compact causal slabs with smooth coefficients and uniformly positive G.
Without this infrared control, the theorem should be stated for the standard
periodic density rather than asserting equality for every conceivable
aperiodic exhaustion of an infinite-area surface.

## 8. Relation to existing symmetry arguments

Appendix B of Grado-White–Grimaldi–Headrick–Hubeny,
*Testing holographic entropy inequalities in 2+1 dimensions*,
[arXiv:2407.07165v2](https://arxiv.org/html/2407.07165v2), develops a
dimension-independent HRT symmetry argument. Its Corollary B.3 covers a
**finite** isometry group. The appendix assumes, for every HRT surface, that
the principal eigenvalues of the two null diagonal Jacobi operators are not
both zero. It therefore excludes certain continuous degenerate families.

The proof above is a different route for the translational/torus setting. It
does not require that Jacobi genericity assumption, but keeps stable reduced
maximin existence and the ordinary null-representative framework explicit.

## 9. Consequence for the intended strip-cone application

For every union of invariant strips, (12) supplies exactly the same curve
minimization data that would be obtained by assuming invariant HRT surfaces,
now with a comparison against all nonsymmetric competitors. Full
strong-subadditivity follows from the original Einstein/NCC spacetime.
Effective NCC is additionally available from (5).

This note does not prove the graph/Crofton representation for those curve
data. That is a separate step, including its topology, homology,
renormalization, and any requirement of a common minimizing slice for
different endpoint pairings. The symmetry theorem cannot repair a failure
of those later hypotheses by itself.
