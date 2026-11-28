# SEARCH4D: covariant HRT vs. static RT cone in AdS4 (search log)

Started 2026-11-28. All work in this directory; nothing outside it is edited.
Scripts run with `python3 -I`. Facets: `../strips/data/facets_n5_all372.json`
(orbits: [0] MMI 3 singletons (20), [1] MMI with composite (45), [2] SA (15),
[3..7] the five genuinely 5-party orbits of sizes 90, 72, 60, 60, 10).

## 0. Plan

(a) SO(2)-symmetric concentric discs/annuli in spherically symmetric
Vaidya-AdS4 (Poincare). HRT surfaces reduce to curves in the orbit space
(v,z,rho) with effective metric G = (rho^2/z^4)(-f dv^2 - 2 dv dz + drho^2).
Pieces: tubes between two edges, caps from an edge to the axis rho=0.
Entropy of any party union = min over admissible configurations of pieces.

Before numerics: determine abstractly which entropy vectors this "matching
formula with a free axis" can produce (Section 1).

## 1. Structural result: the SO(2) class cannot violate (abstract part)

**Reduced problem.** For SO(2)-symmetric regions (concentric discs/annuli,
or polar caps/bands in global AdS4) in any SO(2)-invariant spacetime, and
assuming the HRT pieces are SO(2)-invariant (averaging argument as in the
strips note), every HRT surface is a non-crossing union of curves in the
orbit space: tubes between two edges (area L_ij), caps from an edge to the
axis (c_i).  The orbit space is a disc whose boundary circle is the asymptotic
boundary ray plus the axis, which is a *free arc*: curves may end on it at any
point, any number of times, at no extra cost.  Homology is automatic (each
non-crossing configuration 2-colours consistently).  So for any union X,

  S(X) = min over non-crossing configurations {tubes, caps} on the active edges,
         caps never nested under a tube   ("free-arc matching formula").

**Lemma A (abstract, numerically checked).** Put the axis as one extra point 0
on the circle in front of e_1 and define d(0,e_i)=c_i, d(e_i,e_j)=
L'_ij := min(L_ij, c_i+c_j).  If d is Kalmanson for the circular order
(0,e_1,...,e_m) then the free-arc matching vector lies in the RT cone for any
assignment of annuli to parties; raw tube areas above c_i+c_j are irrelevant.
Proof sketch: Kalmanson + triangle = circular-decomposable (Christopher-Farach-
Trick / Chepoi-Fichet); expand the axis into k copies at mutual distance 0,
apply the disc (2+1) statement (kinematic-space / split decomposition gives a
planar graph whose cuts reproduce the matching formula, the `ensure no endpoint
counterterm needed' because trivial-split weights are non-negative), then
minimise over the free memberships X of the copy sub-arcs, which is the same
as making the copy vertices bulk vertices.  Tests (`abstract_axis.py`,
`abstract_axis_run1.txt`, 1500 data sets x 30 labellings each, m=5..10 edges,
all 357 balanced facets): disc (no axis) Kalmanson, circular-decomposable with
axis, LP vertices of the K1,K2,K4 cone (316 with K0 violated), and 'inflate'
(16258 raw tube areas pushed above c_i+c_j or to infinity): min margin
>= -6e-14 (round-off) in every orbit.  Control with unconstrained random data:
margins down to -5 (all orbits violated), so the test has teeth.

**Lemma B (physics: these conditions are SSA instances).** With i<j<k<l edges,
D_r the disc of radius r, [a,b] the annulus:
 * K1  L'_ik + c_j >= L'_ij + c_k : SSA for X=[i,k], Y=D_j (X u Y = D_k, X n Y=[i,j]).
 * K2  L'_ik + c_j >= c_i + L'_jk : weak monotonicity for the same pair.
 * K4  L'_ik + L'_jl >= max(L'_ij + L'_kl, L'_il + L'_jk) : SSA / WM for X=[i,k], Y=[j,l].
 * K0  c_i + c_j >= L'_ij : by definition of L'.
Physical entropies: S(D_r)=c_r, S([a,b])=L'_ab.  SSA and WM hold for HRT in
any dimension under NEC (Wall's maximin).  Hence Lemma A applies.

**Consequence.** No covariant violation of ANY RT-cone inequality (any n) can
come from concentric annuli in any SO(2)-symmetric NEC spacetime: thin/thick
Vaidya shells, several shells, Poincare or global (polar caps/bands are the
same disc-with-one-free-arc problem), arbitrary boundary times of the edges.
The axis is NOT an obstruction: it acts like zero-cost bulk vertices.

**Two free arcs** (e.g. two-sided global black hole: orbit-space circle =
bdry1, S-axis, bdry2, N-axis; axis-to-axis chords = closed spheres such as
the horizon): same mechanism with the metric closure d' (paths through free
points).  Each Kalmanson condition on d' is again an SSA/WM instance; e.g. for
(N, a, S, b): S(Q1(a) u Q2(b)) + S(bdry1) >= S(P2(b)) + S(Q1(a)) gives
L'_ab + H >= c'_N(b) + c'_S(a) (purity used for the union).  Checked with the
copies engine (`circle_engine.py`, `abstract_twoarc.py`): circular-decomposable
two-arc data and perturbed/inflated data that pass all physical SSA/WM
instances (all unions of asymptotic arcs): no facet violation.

**Adversarial check of the logic** (`abstract_hillclimb.py`, hill_one.txt,
hill_two.txt): 40 one-arc and 30 two-arc chains, 300-400 moves each,
minimising the worst balanced facet margin subject to ALL physical SSA/WM
slacks >= 0: never below -3e-15.

**All symmetry reductions of Vaidya are covered.** Planar Vaidya-AdS4 has
spatial isometry ISO(2); its one-parameter subgroups are translations (strips:
strips theorem) and rotations about a point (concentric annuli: above).
Global Vaidya has SO(3); one-parameter subgroups are rotations about an axis
(polar caps/bands: above).  So every configuration of Vaidya-AdS4 whose HRT
problem reduces to curves gives entropy vectors in the RT cone.  A violation in
Vaidya-AdS4, if it exists, needs regions without a continuous symmetry, where
HRT surfaces with three or more boundary components (pants) appear and the
pairwise matching structure fails.  (Pure pairwise-tube matching among
non-concentric discs would violate even MMI; Wall's theorem is what forces
the multi-boundary surfaces in.)

The numerics below are therefore a *control* of this prediction and of the
code, not a search with a realistic chance of success.
