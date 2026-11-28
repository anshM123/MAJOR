# A finite graph theorem for ordered interval data

This note supplies the full-cone step for the strip project. It separates three
logically different statements: an exact finite-data graph theorem; the physical
matching formula that connects those data to HRT; and the treatment of UV
counterterms. None of them requires differentiability of the interval entropy,
a continuum Crofton density, or completeness of contraction maps.

The source `fullcone/FULLCONE.md` at repository revision
`acda6f59a240537f521583630ff9e9fac71e145f` uses a hyperbolic Liouville measure and an
exact logarithmic horoball identity. That exact identity does not transfer
verbatim to the higher-dimensional strip metric. The finite construction below
implements its graph-realization idea directly from endpoint data instead.

## 1. Data and statement

Place $n\geq2$ marked points $p_0,\ldots,p_{n-1}$ in cyclic order on a disc
boundary. Let $Z_i$ be the open boundary atom from $p_i$ to $p_{i+1}$; all indices
are modulo $n$. A real symmetric matrix $D$ has $D_{ii}=0$. Define, for $i\ne j$,

$$
w_{ij}=\frac12\left(
D_{ij}+D_{i+1,j+1}-D_{i,j+1}-D_{i+1,j}
\right).                                                     \tag{1}
$$

Assume every $w_{ij}\geq0$. This is the finite circular-split condition, with
the adjacent cases included. For a selection $A$ of boundary atoms let
$\partial A$ be the marked points where membership changes, and define

$$
F_D(A)=\min_{\Pi\in\mathrm{NC}(\partial A)}
       \sum_{\{a,b\}\in\Pi}D_{ab},                          \tag{2}
$$

where $\mathrm{NC}(\partial A)$ denotes noncrossing perfect matchings in the
disc. The empty matching has value zero. Atoms can be grouped into arbitrarily
many disconnected parties, with one label designated the purifier.

**Finite graph theorem.** There is one finite undirected graph with nonnegative
edge capacities whose min-cut entropy for every union of parties is exactly
$F_D(A)$. Consequently the complete vector of these values lies in the full RT
holographic entropy cone.

The last implication uses the graph characterization of the RT cone in Bao et
al., arXiv:1505.07839. The construction and equality are proved below. The
conditional-mutual-information weights (1) are already present in
Czech–Dong, arXiv:1905.03787, eq. (2.11); they are not a new identity claimed
here. The direct graph is a finite specialization of the repository's current
arrangement argument.

## 2. Exact inversion by telescoping

For $a\ne b$, let $I(a,b)$ be the atoms on the oriented boundary arc from
$p_a$ to $p_b$. Then

$$
\boxed{D_{ab}=
\sum_{i\in I(a,b),\ j\notin I(a,b)} w_{ij}.}                 \tag{3}
$$

To prove it, cyclically relabel so $a<b$. Twice the right side is

$$
\sum_{i=a}^{b-1}\sum_{j=b}^{a-1\ ({\rm cyclic})}
 (D_{ij}+D_{i+1,j+1}-D_{i,j+1}-D_{i+1,j}).
$$

The inner difference first telescopes in $i$ to
$D_{a,j}-D_{b,j}-D_{a,j+1}+D_{b,j+1}$, then in the complementary cyclic
$j$-interval to
$D_{ab}-D_{bb}-D_{aa}+D_{ba}=2D_{ab}$. No regularity assumption or limiting
argument enters. In particular $D_{ab}\geq0$ follows from the weights.

One should distinguish (3) from a graph model for all unions. The complete
Bell-pair graph on the atoms, with capacities $w_{ij}$, reproduces **single
intervals**, but generally gives the wrong entropy for a disconnected union.
For four atoms and $A=Z_0\cup Z_2$, that graph gives $B(A)=\sum_{i\in A,j\notin A}w_{ij}$,
whereas (2) gives

$$
F_D(A)=B(A)+2\min(w_{02},w_{13}).                             \tag{4}
$$

The arrangement graph below supplies the additional internal degrees of freedom.

## 3. Constructing the graph

Use straight chords in an auxiliary round disc. Let $\mathcal A$ be the union
of every chord $[p_a,p_b]$, not only the winning matching for one region. This
is a finite embedded arrangement after subdividing at its intersections.
For each unordered pair $\{i,j\}$ choose a straight **test chord** $\tau_{ij}$
with one endpoint in $Z_i$ and one in $Z_j$. Choose these chords generically so
none passes through an arrangement vertex or overlaps an arrangement edge.
There are finitely many forbidden conditions, so such choices exist.

Give test chord $\tau_{ij}$ the weight $w_{ij}$. For each arrangement edge $e$,
set

$$
c(e)=\sum_{i<j}w_{ij}\,\#(\tau_{ij}\cap e).                 \tag{5}
$$

All capacities are finite and nonnegative. Take the complementary cells of
$\mathcal A$ as graph vertices; an arrangement edge gives a graph edge between
its two incident cells with capacity (5). The boundary atom $Z_i$ lies in the
closure of a unique cell. For $n\geq3$ these boundary cells are distinct,
because the adjacent-endpoint chord cuts off a cap containing only $Z_i$.
For $n=2$ the one chord gives the two cells directly. Label each boundary cell
by its atom's party. Merge graph vertices carrying the same party label.

Every endpoint chord $[p_a,p_b]$ crosses a test chord exactly when its two test
endpoints lie on opposite boundary arcs. Thus (3) says that its total
intersection cost is precisely $D_{ab}$.

### Upper bound

A noncrossing matching $\Pi$ consists of arrangement chords and bounds a
checkerboard homology region whose boundary trace is $A$. Select the cells in
that region. Its graph-cut edges are exactly the matching's subdivided edges.
Using (5) and (3), its capacity is $\sum_{\{a,b\}\in\Pi}D_{ab}$. Hence

$$
\operatorname{mincut}_A(G)\leq F_D(A).
$$

### Lower bound

Conversely, let a graph cut select any collection of cells with the required
boundary labels. Its frontier is a finite collection of arrangement edges.
At interior vertices, resolve the frontier into disjoint curves according to
the adjacent cell colors. At each marked boundary point, pair the incident
rays, leaving one arc end exactly when the two neighboring atoms have different
membership. All such surgery can be performed in tiny discs or half-discs
avoiding **every** test chord: there are only finitely many test chords, and
they miss all arrangement vertices and all marked points.

Delete closed frontier components, which carry no boundary endpoint and can be
removed by reversing the color of their innermost disc. The remaining curves
give a noncrossing matching $\Pi$. Straightening each arc to its chord cannot
increase intersections with any test chord: the minimum is one for alternating
endpoint pairs and zero otherwise. Therefore

$$
\sum_{\{a,b\}\in\Pi}D_{ab}
\leq\sum_{e\text{ in the graph cut}}c(e).
$$

Taking minima proves equality. This proves the theorem for every atom selection
simultaneously. Identifying all same-party terminals preserves these cuts, so it
also proves the theorem for arbitrary unions of atoms assigned to each party.

## 4. Which entropy inequalities ensure positive weights?

For a pure state on the cyclic atoms $Z_0,\ldots,Z_{n-1}$, set

$$
D_{ab}=S(Z_a\cdots Z_{b-1}),\qquad D_{aa}=0.
$$

Purity makes $D$ symmetric. For $i<j$, equation (1) becomes

$$
\boxed{w_{ij}=	frac12 I(Z_i:Z_j\mid Z_{i+1}\cdots Z_{j-1})\geq0.} \tag{6}
$$

For adjacent atoms the conditioning system is empty; their weight is half an
ordinary mutual information. This is why adjacent positivity uses subadditivity,
not just the four-distinct-endpoint crossing inequality.

For finitely many endpoints on a line, let $Z_0,\ldots,Z_{n-2}$ be the bounded
atoms between consecutive endpoints. Put the two exterior rays, together with
any reference system purifying a mixed global state, into the final atom
$Z_{n-1}$. This produces an abstract pure cyclic system. Only entropies of
bounded connected intervals are needed to specify $D$; any interval containing
the exterior atom is assigned the entropy of its bounded complement. Thus a
thermal state causes no problem at the level of finite interval data.

When expressed solely in the original unpurified state, some of the cyclic
positivity conditions are weak monotonicity. For $a<b<c<d$, the two required
crossing inequalities are

$$
D_{ac}+D_{bd}\geq D_{ad}+D_{bc},\qquad
D_{ac}+D_{bd}\geq D_{ab}+D_{cd}.                              \tag{7}
$$

The first is strong subadditivity for adjacent intervals $AB,B C$; the second
is weak monotonicity for three consecutive regions, equivalently SSA in a
purification. A continuum mixed derivative away from the diagonal captures
the first family locally but, by itself, does not supply all global/cyclic
conditions or the adjacent mutual informations.

There is also a direct bulk proof of weak monotonicity that does not assume a
geometric purification. Choose the common maximin slice on which the HRT
surfaces for the disjoint intervals $[a,b]$ and $[c,d]$ are both minimal,
using Wall's Corollary 17(h). The HRT surfaces for $[a,c]$ and $[b,d]$ have
null representatives on this slice with no greater area. In the disc their
alternating endpoint arcs cross. Cut and rejoin them into candidates for
$[a,b]$ and $[c,d]$, discarding any closed components. The operation does not
increase total area. Minimality of the first two HRT surfaces on that slice
then gives the second inequality in (7). Thus the abstract reference system
above is only an entropic shorthand; no extra asymptotic boundary or bulk
purification is needed for this step.

For infinite transverse directions, these statements apply first in a finite
periodic regulator and then to entropy per transverse volume. Positivity
survives the density limit if the required entropy densities exist.

## 5. The physical matching formula: what has to be proved separately

Suppose the symmetry step has supplied a translation-invariant least-area HRT
representative for every bounded strip union. In addition assume:

1. the translational action is free and spacelike;
2. the relevant orbit-space Cauchy topology is a disc with one boundary line;
3. there are no additional internal boundary components or relative-homology
   endpoints at which a reduced surface may terminate;
4. the usual stable maximin/HRT existence and nesting statements apply.

Then the reduced HRT entropy obeys (2), with $D_{ab}$ the single-strip HRT
entropy density between its two endpoints.

### Lower bound does not require componentwise minimality

An invariant union surface reduces to disjoint extremal curves. Its global
acausality and an invariant timelike-flow projection put their spatial traces
on the disc without intersections. A closed curve bounds a disc; deleting that
component and reversing the homology color inside removes positive area while
leaving the boundary region unchanged. Thus a least-area representative has
only arcs and determines a noncrossing endpoint matching.

Each arc separately bounds its bounded boundary interval on a Cauchy slice, so
its invariant lift is an admissible single-strip extremal candidate. Its area
density is therefore at least $D_{ab}$. Summing gives the lower bound in (2).
It is not necessary to assert that every component of the union is already the
single-interval minimizer.

The timelike-flow projection can be made equivariant after the finite torus
regulator: average a future timelike vector field over the compact translation
group (the future timelike cone is convex), then rescale invariantly to make
it complete. Projection of an invariant acausal surface to an invariant
reference Cauchy surface is injective and equivariant. Its quotient traces are
therefore embedded arcs on the spatial disc.

### Upper bound uses a finite laminar common-slice lemma

The bounded intervals associated with a noncrossing matching form a laminar
family: any two are disjoint or nested. Wall's Theorem 17 proves common-slice
minimality for a nested pair, and Corollary 17(h) for a disjoint family. The
same positive weighted joint-maximin proof applies to any finite laminar
family: minimize all areas on one slice, then maximize their positively
weighted sum. Pairwise uncrossing gives disjoint components or coincident
components. Local slice variations at disjoint components establish individual
extremality; coincident components are grouped with the sum of their positive
weights. The null representatives then establish individual global HRT
minimality exactly as in Theorem 17(f). The generic/stability and existence
conditions are the same as in that theorem, not consequences of NEC alone.

Choose the single-strip HRT surfaces for all those intervals on the resulting
common slice. Their symmetric difference homology region has boundary exactly
the target strip union, because the matching endpoints alternate its boundary
membership. The disjoint surface union is extremal and has the sum of the
single-strip entropies. This gives an admissible candidate attaining the
matching cost, and proves the upper bound.

Translation invariance by itself does not establish the disc-topology
hypothesis. Horizons, extra asymptotic ends, internal boundaries, or nontrivial
orbit-space topology can introduce further components and homology classes.
Likewise, a finite strip party must have invariant entangling edges: a fully
arbitrary boundary Cauchy slice with $y$-dependent edge times is not an invariant
boundary-value problem.

## 6. UV renormalization and exact positive endpoint padding

Let $d_{ab}=d_{ba}$ be finite renormalized interval data, with the convention
$d_{aa}=0$. Define $w^0_{ij}$ by (1) with $d$ in place of $D$. For four
distinct endpoints all local endpoint divergences cancel, so SSA/weak
monotonicity imply $w^0_{ij}\geq0$ for nonadjacent atoms. The adjacent weights
can be negative; renormalized entropy is not itself an ordinary quantum entropy.

For $n\geq3$, put

$$
a_p=\max\{0,-w^0_{p-1,p}\},\qquad
D^*_{ab}=d_{ab}+a_a+a_b\quad(a\ne b),\qquad D^*_{aa}=0.        \tag{8}
$$

Adding endpoint terms changes only the adjacent split weights:

$$
w^*_{p-1,p}=w^0_{p-1,p}+a_p\geq0,
\qquad w^*_{ij}=w^0_{ij}\quad\text{for nonadjacent }i,j.       \tag{9}
$$

Therefore $D^*$ has an exact positive graph realization. Every matching uses
each selected endpoint once, so

$$
G(A):=F_{D^*}(A)=F_d(A)+\sum_{p\in\partial A}a_p.             \tag{10}
$$

Writing $b_p(A)=[p\in\partial A]$, this gives the exact endpoint-counterterm
statement

$$
\boxed{
F^{\rm ct}=F_d+\sum_p C_p b_p
=G+\sum_p(C_p-a_p)b_p\in\mathcal C_N
\quad\text{whenever every }C_p\geq a_p.
}                                                            \tag{11}
$$

Each $b_p$ is a Bell-pair graph ray between the labels on the two sides of the
endpoint (zero if the two labels coincide). This proves a positive graph
decomposition for any local endpoint-counterterm convention with sufficiently
large positive UV terms. The $n=2$ case is the single condition
$d_{01}+C_0+C_1\geq0$.

**The finite renormalized remainder need not lie in the RT cone.** A single
vacuum AdS4 strip already has a negative finite area term proportional to
$-1/\ell$, whereas a one-party nonnegative graph cut cannot be negative. The
finite positive graph core in (10) includes the explicit finite endpoint
padding; it should not be identified with a bare subtraction-scheme finite part.

### A literal radial truncation needs one more regulator statement

For a general geometry one typically only has

$$
D^{\varepsilon}_{ab}
=d_{ab}+C_a(\varepsilon)+C_b(\varepsilon)+r^{\varepsilon}_{ab},
\qquad r^{\varepsilon}_{ab}\longrightarrow0.                  \tag{12}
$$

The hyperbolic horoball identity in `FULLCONE.md` sets the analogous remainder
exactly to zero. Equation (12) alone cannot establish exact finite-cutoff cone
membership when a limiting nonadjacent weight vanishes. A small negative
remainder can make that weight negative.

An explicit algebraic control is $n=4$, $D^0_{ab}=2U$ for $a\ne b$ and
$D^0_{aa}=0$. It has adjacent split weights $U$ and the two nonadjacent weights
zero. Subtract $\varepsilon$ times the circular split metric separating
$\{p_1,p_2\}$ from $\{p_3,p_0\}$. Every pair distance remains positive for
$U>\varepsilon$, but

$$
D^{\varepsilon}_{02}+D^{\varepsilon}_{13}
-D^{\varepsilon}_{03}-D^{\varepsilon}_{12}=-2\varepsilon<0.   \tag{13}
$$

Thus arbitrarily small remainders can spoil an exactly saturated SSA condition.
This is an obstruction to that inference, not a physical HRT counterexample.

The valid alternatives are:

- use exact common-regulator interval data known to obey all conditions (1);
- use the exact local endpoint-counterterm convention (11);
- or state the standard radial-truncation conclusion asymptotically.

For the last alternative, if $r_\varepsilon=\max_{a\ne b}|r^\varepsilon_{ab}|$,
the minimum of finitely many matching costs obeys

$$
|F_{D^\varepsilon}(A)-F^{\rm ct}(A)|
\leq\frac{|\partial A|}{2}\,r_\varepsilon.                  \tag{14}
$$

Consequently the distance of the finite entropy vector from the RT cone is
$O(r_\varepsilon)$ for the matching prescription $F_{D^\varepsilon}$; every
endpoint-balanced RT inequality holds for the finite
renormalized data, and every inequality with a nonzero positive Bell-pair
divergence becomes positive for sufficiently small cutoff. An exactly
saturated endpoint-balanced inequality needs more information to determine the
sign of its finite-cutoff remainder.

A second finite-cutoff issue is logically separate from the sign of the split
weights. Truncating a selected continuum HRT surface for a union is not by
definition the same as minimizing the sums of truncated single-interval
surfaces. At a continuum phase degeneracy, the two matching sums agree after
renormalization but their vanishing cutoff remainders can differ. A chosen
continuum winner need not be the smaller truncated branch. The same issue can
occur between degenerate single-interval HRT representatives. The exact
finite-data theorem applies to $F_{D^\varepsilon}$; identifying this with an
independently chosen literal truncated HRT vector requires a compatible branch
prescription or a nondegeneracy argument. For a fixed finite family, ordinary
asymptotic counterterm subtraction still gives an $o(1)$ difference between
the literal vector and its exact endpoint completion, provided that every
selected HRT surface has the assumed common local counterterm expansion.

One must not silently invoke finite-cutoff HRT SSA as a general theorem.
`Radial cutoffs and holographic entanglement` distinguishes ordinary extremal
surfaces from the restricted-maximin prescription precisely because finite
cutoffs change the causal problem. Whatever regulator is chosen must be stated.

### Optional strictness and a phase-stable radial-cutoff refinement

Under smooth common-maximin-slice regularity, all four-distinct-endpoint
inequalities in (7) are strict in the disc class. This can be proved without
strict null focusing. For SSA choose a common slice on which $[a,d]$ and
$[b,c]$ are minimal; for weak monotonicity choose one on which $[a,b]$ and
$[c,d]$ are minimal. On that slice, minimize separately for $[a,c]$ and
$[b,d]$. Their intrinsic minima are no larger than their global maximin
values. The two smooth minimizing curves have alternating endpoints, hence
intersect. They cannot meet tangentially: uniqueness of the Riemannian geodesic
initial-value problem would make them the same curve, incompatible with their
distinct boundary endpoints. Rejoining and rounding the transverse corner
strictly shortens the total length. Comparing the resulting candidates with
the two prescribed minima yields the strict version of (7).

This argument assumes that the required common slice admits a smooth,
positive-definite reduced metric and regular minimizing curves. When those
facts are not established, use nonnegative weights and the exact endpoint
completion instead. Wall's Theorem 18(d) supplies the same strict-SSA conclusion
under its saturation hypotheses: saturation requires the two overlapping HRT
surfaces to be minimal on a common slice without crossing, which alternating
connected disc arcs cannot do.

For a fixed finite endpoint set, strictness gives a positive minimum of all
nonadjacent split weights. Then every sufficiently small remainder in (12)
preserves their signs, while the adjacent weights are eventually positive from
their UV divergences. If, additionally, each single-interval HRT curve and each
winning union matching is unique, the literal truncated vector equals
$F_{D^\varepsilon}$ for sufficiently small $\varepsilon$. Indeed, the continuum
matching proof forces every winning union component to be the corresponding
unique single-interval HRT curve; the finitely many positive gaps between a
winning matching and its competitors survive a sufficiently small remainder.
The exact finite graph theorem then proves exact radial-cutoff cone membership
in this nondegenerate, smooth-slice regime. This optional refinement is not
needed for the exact counterterm-regulated theorem.

## 7. Scope of the resulting cone statement

Given symmetry, the disc matching formula, and exact regulated positivity,
the theorem proves

$$
\{\text{covariant strip entropy vectors in this class}\}
\subseteq\mathcal C_N^{\rm RT}.
$$

This is sufficient to prove every static holographic entropy inequality for
these dynamics. It does not establish that every vector or extreme ray of the
full RT cone is realizable by static strips in the same geometrical class.
Calling the result equality of those two restricted realization sets would
require that reverse-realization theorem separately.

For continuum HRT data with the usual local UV expansion, the unconditional
regulator formulation furnished by this note is: **the positively
endpoint-completed vector (11) lies in the full RT cone, and every
endpoint-balanced RT inequality is satisfied independently of this endpoint
completion.** An exact common regulator whose interval data obey (1) and whose
union entropies obey (2) also gives exact full-cone membership by the same finite
graph theorem. Literal radial truncation has the additional qualifications
spelled out above.

## 8. Reproducible finite-graph validation

Run `python verify_finite_graph.py --max-n 8 --output finite_graph_validation.json`
from this directory. The script uses only the Python standard library. It
constructs the full chord-arrangement dual with exact rational intersections,
assigns integer crossing capacities, and computes unrestricted graph mincuts
with an integer maximum-flow algorithm. An independent noncrossing-matching
dynamic program gives the comparison values.

For three through eight endpoints, three weight families (dense unit weights,
random weights including zeros, and sparse crossing weights), and every boundary
selection up to complementation, all **756 exact equalities** pass. The three-
and four-endpoint cases additionally enumerate every interior cell coloring,
independently checking the maximum-flow output. A four-atom control has complete
Bell-pair-graph cut zero but arrangement-graph cut and matching minimum ten,
verifying that the construction tests the disconnected-union issue as well as
the split reconstruction. These finite checks support the implementation; the
proof of the theorem is the construction and surgery argument in Sections 2–3.

## Primary references

- N. Bao et al., *The Holographic Entropy Cone*, arXiv:1505.07839. Graph
  characterization; static geometry-to-graph construction.
- B. Czech and X. Dong, *Holographic Entropy Cone with Time Dependence in Two
  Dimensions*, arXiv:1905.03787, especially eq. (2.11), §§2.3–2.5, and §7.
- A. C. Wall, *Maximin Surfaces, and the Strong Subadditivity of the Covariant
  Holographic Entanglement Entropy*, arXiv:1211.3494, Theorem 17 and Corollary
  17(h), Theorem 18. These supply the cited nesting and SSA input under their
  stated maximin hypotheses.
- B. Grado-White, D. Marolf, and S. J. Weinberg, *Radial Cutoffs and Holographic
  Entanglement*, arXiv:2008.07022, JHEP 01 (2021) 009. Its distinction between
  ordinary finite-boundary HRT and restricted maximin explains why a regulator
  hypothesis cannot be omitted.
- Repository `fullcone/FULLCONE.md`, revision
  `acda6f59a240537f521583630ff9e9fac71e145f`, §§1–3: the current-based graph
  construction and Bell-pair endpoint decomposition being generalized here.
