# Static graph models for covariant strip entropies

Research note, 28 November 2026. Based on `anshM123/MAJOR` revision
`acda6f59a240537f521583630ff9e9fac71e145f`.

## Result and its precise scope

The reduction and full-cone argument close in the **regular disk-orbit
maximin class**. The result applies at arbitrary strip widths and boundary
times; it uses neither a membrane scaling limit nor completeness of
contraction maps. The proof supplies an invariant HRT representative by
comparison with every nonsymmetric competitor, rather than assuming symmetry
of the minimizer.

The literal statement “translation invariance and NEC alone, on any boundary
Cauchy slice” is broader than the theorem proved here. The entangling edges
must be invariant, the relevant spatial orbit space must have disk topology,
and the usual maximin existence/regularity and null-representative hypotheses
must hold. In particular, this note does not deduce stable reduced-maximin
existence from the bare existence of a full HRT surface. The proof of symmetry
keeps that global input explicit. The invariant compactness argument and a
black-brane singularity barrier are detailed in the symmetry note.

### Theorem

Let $(M,g_D)$ be a classical asymptotically locally AdS Einstein spacetime of
dimension $D=d+1=m+3$, $m=d-2\ge1$, satisfying null convergence. Assume:

1. **Regular translations.** There is a free spacelike translation action in
   $m$ directions, with positive-definite orbit metric. Transverse entropy
   densities are defined by admissible periodic quotients and their
   decompactification limit. The metric may have anisotropic orbit sizes and
   nonzero time-orbit or space-orbit cross terms.
2. **Invariant boundary data.** The finitely many parties are disjoint unions
   of bounded parallel strips over a boundary Cauchy curve in the orbit space.
   Their entangling edges and domains are translation invariant. The quotient
   curve may have arbitrary spacelike time profile $t=t(x)$.
3. **Disk homology.** The relevant orbit-space Cauchy surface is a disk or
   half-plane. The cuts have only the specified boundary endpoints; there are
   no extra internal boundaries or relative-homology endpoints. The cuts are
   finite and proper. This is an independent topological hypothesis.
4. **Maximin admissibility.** The weighted maximin problems in the invariant
   slice class have regular stable solutions for every interval between
   marked endpoints and every required strip union. The full spacetime's
   null representatives and positive weighted joint-maximin construction
   apply to the finite laminar families of those intervals. At the anchors,
   same-anchor extremal asymptotics and conformally regular spacelike
   slices hold, and the representatives agree with their short local null
   projections. Sufficient boundary acausality and compactness conditions
   for this last property are proved in `UV_AVERAGING.md`. All areas use
   one consistent asymptotic boundary prescription.

Let $s^{\rm ren}(A)=S^{\rm ren}_{\rm HRT}(A)/V_y$ be the entropy density after
local endpoint subtraction, and let $p$ run over all marked strip endpoints.
Define $b_p(A)=1$ when $p$ belongs to the entangling boundary of $A$, and zero
otherwise. Then there are finite constants $a_p\ge0$ and **one finite graph
with nonnegative capacities**, $G$, such that, simultaneously for every union
of parties,

$$
\boxed{\ s^{\rm ren}(A)+\sum_p a_p b_p(A)=S_G(A).\ }       \tag{T1}
$$

Consequently, for every local endpoint-counterterm prescription
$C_p\ge a_p$,

$$
\boxed{\ s^{\rm ct}(A):=s^{\rm ren}(A)+\sum_pC_p b_p(A)
=S_G(A)+\sum_p(C_p-a_p)b_p(A)\ \in\mathcal C_N^{\rm RT}.\ } \tag{T2}
$$

Here the last membership refers to the complete entropy vector as $A$ ranges
over the $2^N-1$ nonempty party unions. The same exact membership holds for
any finite regulated interval data for which the matching formula and the
cyclic SSA/weak-monotonicity conditions hold exactly. Every endpoint-balanced
RT inequality holds for $s^{\rm ren}$ itself.

Here endpoint-balanced means that the coefficients $q_A$ of the inequality
obey $\sum_A q_A b_p(A)=0$ for every marked endpoint $p$.

The conclusion is **containment in the full RT cone for every finite number
of parties**. It does not assert that every RT-cone ray is realizable by
static strips in the same disk geometry. All entropies here are the classical
area/HRT contribution.

## Proof architecture

### 1. The correct effective metric, including boosts

Write the invariant metric in Kaluza--Klein form,

$$
g_D=g_{ab}dx^a dx^b+
G_{ij}(dy^i+A^i_a dx^a)(dy^j+A^j_b dx^b),
\qquad G>0.
$$

Put $\lambda=\sqrt{\det G}$, $u=\log\lambda$, and
$\widehat g=\lambda^2g$. A block-determinant calculation gives, exactly,

$$
\frac{\operatorname{Area}(\pi^{-1}\gamma)}{V_y}
=\int_\gamma\lambda\,ds_g
=\operatorname{Length}_{\widehat g}(\gamma).               \tag{T3}
$$

For isotropic fibers, $G=e^{2\phi}I_m$, this is
$\widehat g=e^{2m\phi}g$ as proposed. Cross terms cancel in the Schur
complement; they need not be excluded.

There is also a useful strengthening of the proposed reduction. **The
effective metric preserves null convergence.** Set
$P_a=G^{-1}\partial_aG$, $F^i=dA^i$. For a base-null vector $k$ and its
horizontal lift $K$,

$$
\operatorname{Ric}_{\widehat g}(k,k)
=\operatorname{Ric}_D(K,K)+(ku)^2+
\frac14\operatorname{Tr}(P_k^2)
+\frac12G_{ij}(\iota_kF^i)_a(\iota_kF^j)^a\ge0.            \tag{T4}
$$

The matrix term is a trace square of a matrix similar to a real symmetric
matrix. The gauge-field covectors annihilate the null vector, so their
Lorentzian norm matrix is positive semidefinite. In the isotropic untwisted
case the extra contribution is $m(m+1)(k\phi)^2$.

This is an Einstein-frame reduction. It is generally not asymptotically
AdS$_3$: near the Poincare boundary $\widehat g$ has conformal factor
$z^{-2(d-1)}$. Maximin is therefore applied to the original area-density
functional with the stated causal and existence hypotheses, not by silently
assuming every AdS$_3$ asymptotic theorem applies.

### 2. Symmetry inheritance is a variational theorem

On any invariant spacelike slice $\Sigma$, let $R$ be an arbitrary admissible
homology region, with no assumption of symmetry. In a transverse torus
regulator average its characteristic function over the isometry group:

$$
f=\int_{T^m}\chi_{gR}\,d\mu(g).
$$

Convexity of total variation and isometric invariance imply

$$
|Df|(\Sigma)\le P(R).
$$

The averaged function has the same boundary trace. Its level sets are
invariant admissible homology regions, so weighted BV coarea gives

$$
P(R)\ge V_y\int_0^1P_\lambda(\{f>s\})\,ds
\ge V_y\min_{\gamma\sim A,\,\gamma\subset\Sigma/T^m}
\operatorname{Length}_{\widehat g}(\gamma).                \tag{T5}
$$

Thus the invariant spatial minimum is the unrestricted spatial minimum.
This includes disconnected, winding, and nonsymmetric competitors.
The BV statement is proved for finite unrenormalized perimeters with fixed
invariant boundary trace. To pass to infinity, use the reduced minimizer's
own trace as reference. If a competitor's trace differs on $E_\epsilon$,
BV gluing and the same averaging give the exact relative inequality

$$
\operatorname{Area}(X_{*,\epsilon})
\le P_\epsilon(R)+\mathcal H^{d-1}_\Sigma(E_\epsilon).
$$

The shared extremal expansion and local null projection imply an edge-band
width $o(\epsilon^{d-1})$, making the mismatch area $o(V_y)$. The same
estimate aligns the focusing comparison at the two radial cutoffs.
`UV_AVERAGING.md` proves these statements, giving a direct asymptotic
comparison without assuming convergence of cutoff minima. Convexity is
never applied to a subtracted perimeter.

Now take a stable weighted maximin curve $\gamma_*$ and lift it to $X_*$. Its
mean-curvature vector is invariant and horizontal. Stationarity under all
invariant normal variations forces that vector to vanish pointwise; hence
$X_*$ is fully extremal, including against noninvariant variations. For any
other full extremal surface $Y$, its null representative $\widetilde Y$ on
the lifted maximin slice satisfies Wall's focusing comparison,

$$
\operatorname{Area}(X_*)\le\operatorname{Area}(\widetilde Y)
\le\operatorname{Area}(Y).                                \tag{T6}
$$

The first inequality is (T5), the second is Raychaudhuri plus NCC. Therefore
$X_*$ is an invariant HRT representative. No uniqueness or Jacobi
nondegeneracy assumption enters this comparison. Existence and regularity
of the stable weighted maximin remain the explicit global input.
For asymptotic anchors, (T6) is understood with the vanishing cutoff errors
just established, followed by the finite same-anchor area-difference limit.

### 3. Strip unions are noncrossing matching minima

Let $d_{ab}$ be the renormalized single-strip entropy density between marked
endpoints $p_a,p_b$. For a union $A$,

$$
s^{\rm ren}(A)=
\min_{\Pi\in\operatorname{NC}(\partial A)}
\sum_{\{a,b\}\in\Pi}d_{ab}.                               \tag{T7}
$$

For the lower bound, the invariant HRT surface projects to disjoint arcs on
the spatial disk. A closed component can be removed, reversing the homology
color inside it, and decreases area. Each remaining extremal arc is an
admissible single-strip candidate, so its entropy density is at least
$d_{ab}$. Its endpoints form a noncrossing matching.

For the upper bound, the intervals associated with any noncrossing matching
are nested or disjoint. The positive weighted joint-maximin construction for
this finite laminar family supplies single-strip HRT surfaces on one common
slice in the full spacetime. These common-slice representatives need not
themselves be symmetric: their areas equal the already identified
single-strip values. Their disjoint union, with any coincident components
canceled, is extremal; the symmetric difference of their homology regions
has boundary trace $A$. Its entropy is no greater than the matching cost.
Together with the lower bound this proves (T7), including phase transitions
between matchings.

### 4. Positive interval differences produce a static graph

Put the endpoints in cyclic order by including the exterior rays and any
formal purifier in one boundary atom. For finite interval data $D$, define

$$
w_{ij}=\frac12\left(D_{ij}+D_{i+1,j+1}
-D_{i,j+1}-D_{i+1,j}\right).                              \tag{T8}
$$

These are half conditional mutual informations. The nonwrapping cases use
SSA; wrapping cases use weak monotonicity. Adjacent cases are ordinary
mutual informations. Telescoping gives the exact inverse,

$$
D_{ab}=\sum_{i\in I(a,b),\ j\notin I(a,b)}w_{ij}.          \tag{T9}
$$

When every $w_{ij}\ge0$, put all endpoint chords in an auxiliary disk. For
each atom pair choose a generic test chord of weight $w_{ij}$. Subdivide the
endpoint chords at their intersections and give each segment the total
weight of test chords crossing it. The dual cell graph has nonnegative
capacities.

Every noncrossing matching is a graph cut with exactly its cost in (T9).
Conversely, resolve any graph-cut frontier into disjoint arcs, remove closed
loops, and straighten the remaining arcs. This does not increase their
intersection number with any test chord. Therefore every graph cut costs at
least a noncrossing matching. The graph reproduces (T7) for **all** party
unions simultaneously.

An ordinary Bell graph on boundary atoms reproduces only single intervals;
its disconnected-union entropies are generally wrong. The internal cells of
the arrangement graph are essential. The graph-to-RT characterization then
gives the full cone, with no list of inequalities and no contraction-map
completeness assumption.

### 5. The UV completion is explicit and positive

Form (T8) from the finite data $d$, setting $d_{aa}=0$, and call the result
$w^0$. Local endpoint divergences cancel in every nonadjacent weight, so
SSA and weak monotonicity give $w^0_{ij}\ge0$ there. For at least three
endpoints set

$$
a_p=\max\{0,-w^0_{p-1,p}\},\qquad
D^*_{ab}=d_{ab}+a_a+a_b\quad(a\ne b).                     \tag{T10}
$$

This changes only $w^0_{p-1,p}$, adding $a_p$. All weights are now
nonnegative. Every matching uses each selected endpoint once, proving (T1).
The single-edge Bell graph between the party labels adjacent to $p$ has
entropy vector $b_p$. Adding its nonnegative capacity $C_p-a_p$ proves (T2).
For two endpoints it suffices that $d_{01}+C_0+C_1\ge0$.

This is the exact counterpart of the endpoint decomposition in
`fullcone/FULLCONE.md`; it works with power-law strip divergences as well as
logarithmic ones. The bare finite part is not asserted to be a graph vector:
the vacuum AdS$_4$ strip already has a negative finite term proportional to
$-1/w$.

For a flat boundary with the usual AdS radius $L$, the leading term per
endpoint is $C_p(\epsilon)=L^{d-1}/[4G_D(d-2)\epsilon^{d-2}]$.
It tends to positive infinity. For each fixed state and finite endpoint set,
the condition $C_p\ge a_p$ therefore holds at sufficiently small cutoff in
the exact endpoint prescription. A state-independent bound on that cutoff
is neither needed nor asserted.

A literal radial truncation generally has additional vanishing terms,
$D^\epsilon_{ab}=d_{ab}+C_a(\epsilon)+C_b(\epsilon)+o(1)$. Without further
regulator control this proves asymptotic cone membership, not an exact sign
for every saturated inequality at every finite cutoff. The exact statement
(T2) uses a specified endpoint-counterterm prescription. Numerical radial
truncations are tested directly in the accompanying example.

## Numerical illustration and reproducibility

The accompanying thin-shell Vaidya-AdS$_4$ calculation uses

$$
ds^2=z^{-2}\big[-(1-\Theta(v)z^3)dv^2-2dv\,dz+dx^2+dy^2\big],
$$

with AdS radius and final horizon radius set to one. It solves the strip
area functional, enumerates the competing shell-crossing solutions found
by continuation, includes the admissible static branch, and minimizes over
noncrossing pairings for every party union.

The single-strip solver uses the standard reflection-symmetric family in
$x$. This numerical restriction is separate from the transverse-translation
theorem. The matching formula is an analytic input; the computation does not
independently establish symmetry or the existence of common maximin slices.

Two configurations are supplied: contiguous widths
$(0.10,0.14,0.18,0.16,0.12)$ and separated width-$0.12$ strips with gap
$0.0741$. The first 3, 4, and all 5 parties give six cases. All 372
five-party facet rows are evaluated at 551 boundary times, with absent
parties set to the empty system. The complete five-party catalog and its
purifier-permutation expansion are included.

The minimum residual is approximately $-7.11\times10^{-14}$ in area-density
units, with no value below the stated $-2\times10^{-8}$ numerical threshold.
The largest tested interval width is $0.8964$; an analytic short-strip bound
excludes the deep radially returning branches in the symmetric single-shell
family. Small competing time branches near saturation are still included.
Quadrature refinement and a separate multistart root audit are documented.
The separated pair loses mutual information at
$t/z_h\simeq0.2336006115$.

These are sampled floating-point checks, not an interval-arithmetic proof of
all-time positivity. The analytic graph theorem is what covers every time
and every number of parties under its hypotheses. A separate exact-rational
arrangement program verifies 756 graph-mincut/noncrossing-matching equalities
with 3 through 8 endpoints, including zero capacities and a case where the
boundary-only Bell graph fails.

## Relation to the literature and proposed headline

The scoped primary-source search found no exact prior theorem with this
higher-dimensional strip-sector scope. Czech--Dong (arXiv:1905.03787)
establishes the three-dimensional contraction-certifiable result.
Bao--Mezei (arXiv:1811.00019) transfers static inequalities to higher-dimensional
thermalization in the leading large-region/late-time membrane regime.
Current general covariant graph constructions retain extra geometric
conditions: cooperating time sheets in arXiv:2502.09894, or exposed regions
in the revised arXiv:2602.04888v2. The full comparison and source caveats are
recorded in the novelty note.

The result supports the precise headline:

> **Covariant strip entropies admit static graph models in every dimension.**

The text should immediately specify the regular disk-orbit maximin class
and the endpoint prescription. It establishes exact finite-width and
arbitrary-time containment; it does not by itself settle all orbit
topologies or supply a reverse realization of every RT-cone ray by strips.

## Files and proof dependencies

- `SYMMETRY.md`: metric reduction, Ricci identity, BV spatial symmetry,
  extremal comparison, compactness and infrared qualifications.
- `UV_AVERAGING.md`: the exact trace-mismatch inequality and the vanishing
  seam/focusing errors needed at conformal infinity.
- `FINITE_GRAPH.md`: full finite graph proof, physical matching lemma, exact
  endpoint padding, and regulator distinctions.
- `NUMERICS.md`: equations, branch treatment, facet conventions, numerical
  results, figure captions, and reproduction commands.
- `NOVELTY.md`: primary-source audit dated 28 November 2026.
- `code/`: exact graph and Ricci checks; facet expansion.
- `numerics/`, `data/`, `figures/`: executable Vaidya calculation, reusable
  data, validation reports, and publication-size figures.

Primary proof inputs are Wall, arXiv:1211.3494, Theorems 15, 17 and 18, and
the RT graph characterization of Bao et al., arXiv:1505.07839. The weighted
interval differences agree with Czech--Dong eq. (2.11). The present work's
new role is their combination with an unrestricted symmetry comparison and
the explicit finite graph/endpoint construction for higher-dimensional
strip dynamics. The novelty search is evidence about published scope, not
a guarantee of priority or journal acceptance.
