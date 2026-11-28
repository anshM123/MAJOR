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
planar graph whose cuts reproduce the matching formula; no endpoint
counterterm is needed for balanced facets, and the Kalmanson conditions are
invariant under endpoint shifts d(p,q) -> d(p,q)+l_p+l_q, which restore the
triangle inequality for renormalised areas), then
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

## 2. Physical numerics: solver and validation

Code: `vaidya4.py` (G-geodesic equations in h-arclength, near-boundary
expansion rho = rho_p - z^2/(2 rho_p) + alpha z^3, v = t - z + beta z^3 with
free (alpha, beta); renormalised area per 2 pi with the O(z_s^2) end
corrections derived analytically: integrand = rho_p/z^2 + z(m rho_p/2 - 2 alpha)),
`hrt_pieces.py` (all branches: coarse vectorised RK4 over a shooting grid,
triangle preimage search with local batch refinement of stretched cells,
Newton refinement with the DOP853 integrator, continuation from the vacuum
for caps), `vaidya_scan.py` (pools of 8 edges, every tube computed from BOTH
ends, K-slacks, facets).

Validation:
 * pure AdS4 disc: A_ren/2pi = -1 to 1e-8 (exact: -1).
 * tube area computed from the inner end and from the outer end agree to
   ~1e-8 after the analytic end correction (before it: 1e-5 discrepancy,
   exactly the predicted zs^2 (m rho/4 - alpha) term).
 * Liu-Suh early-time universal growth, thin shell (dlt=0.01), R=4:
   Delta(A/2pi) = kappa t^2 with kappa = 0.998 (t=0.25), 0.988 (t=0.5);
   prediction Delta S = (pi/(d-1)) E |dA| t^2  ->  Delta(A/2pi) = M R t^2/4 = 1.
   Saturation for R=4 reached at t ~ 4.2 with A/2pi = 4.1719 (constant after).
   Pseudo-arclength continuation of the R=4 branch shows t(v*) monotone up to
   t = 3.70 (tip z* = 1.67); the last stretch to saturation is exponentially
   stretched (surface hugs the critical surface) and not resolved: pools only
   use edges whose cap was found.
 * static control (all edges at t=-3, pre-shell, pure AdS): K-slacks >= -4e-16,
   all facet orbits >= -7e-13 (round-off).

Bugs found and fixed on the way (each produced spurious "violations" that
were always accompanied by SSA (K-slack) violations, i.e. impossible data):
 1. annuli thinner than the shooting grid (ratio 1.003) -> missed tube; pools now
    require ratio >= 1.07 and gap >= 0.03.
 2. Newton tolerance below the integrator noise floor -> converged roots
    rejected; fixed (accept at 1e-8).
 3. coarse 181^2 grid missed near-fold tube branches (two-end mismatch); now
    321^2 grid + branch propagation by continuation between neighbouring
    targets.

**Admissibility of the matching formula (why HRT = min over non-crossing
configurations).** A crossing configuration is never strictly better than an
uncrossed one (uncrossing uses exactly K1, K2, K4), so the minimum can be
taken over configurations whose pieces are HRT surfaces of nested or disjoint
discs/annuli (pieces with L_ij > c_i + c_j are always beaten).  Such pieces are
mutually spacelike by entanglement-wedge nesting (Wall, NEC), so they lie on a
common achronal slice and the union is an admissible homologous extremal
surface.  Remaining assumption: the minimal HRT surface of an SO(2)-invariant
region is SO(2)-invariant (maximin + averaging, as in the strips note).

**Generalisation.** Nothing used roundness or concentricity: any *chain of
nested* closed boundary curves (non-round, non-concentric, on any spacelike
boundary slice) whose HRT pieces each have at most two boundary components
gives the same free-arc matching formula and hence RT-cone vectors.  What can
evade the argument is (i) pieces with >= 3 boundary components (pants),
(ii) laminar families that are trees rather than chains (disjoint discs inside
a disc), where the pairwise planar order is lost.

## 3. Proof of Lemma A (graph construction; replaces the sketch above)

Let d be the augmented matrix on the circle (F, e_1, ..., e_m), F the axis,
d(e_i,F) = c_i, d(e_i,e_j) = L'_ij = min(L_ij, c_i + c_j).
1. *Endpoint shifts.* d(p,q) -> d(p,q) + l_p + l_q (l_F = 0) changes every
   configuration cost by sum_p l_p b_p(X), so balanced facets are unchanged;
   Kalmanson inequalities are unchanged; for l_p large d becomes a metric
   (triangle inequalities through F hold because L' <= c_i + c_j).
2. *Split decomposition.* A Kalmanson metric for a circular order is circular
   decomposable: d = sum_S a_S delta_S, a_S >= 0, S ranging over arcs of
   consecutive points (Christopher-Farach-Trick 1996; Chepoi-Fichet 1998).
3. *Graph.* Put the points on a convex circle, draw for every split S a
   straight chord with endpoints on the boundary gaps separating S from its
   complement, weight a_S.  Let G be the planar dual: one vertex per face,
   edge capacity a_S across each chord segment; all faces touching one
   asymptotic interval are glued to that interval's boundary vertex
   (capacity infinity); faces touching the axis arc F stay *bulk* vertices.
4. *Cut = matching.*  A min cut of G is a non-crossing system of curves with
   endpoints at the marked edges where membership changes, or anywhere on F;
   a curve from p to q (or to F) crosses every chord separating its endpoints,
   so costs >= d(p,q) (resp. c_p).  Straight segments attain this, and
   non-crossing configurations give non-crossing segments.  Hence
   S_G(X) = min over non-crossing free-arc configurations of sum d.
5. *Raw tubes.* The physical formula uses L >= L'.  Any configuration using
   a d-chord (p,q) with L' = c_p + c_q < L_pq is matched by a physical one of
   no larger cost: if no tube encloses (p,q), replace it by caps p, q; if
   (r,s) encloses it, K2 gives d(r,s) + c_q >= c_r + d(q,s), i.e. the
   physical configuration {cap r, cap p, tube (q,s)} is no worse, and the
   nesting depth drops; iterate.  So S_phys = S_G (numerically confirmed by the
   'inflate' test, 16258 inflated tubes).
Therefore S_phys(X) = S_G(X) for one graph G with non-negative capacities,
simultaneously for all unions X: the vector is in the RT cone for every n.
Two free arcs: same construction with two bulk-face families; axis-to-axis
cut curves cost d(N,S) = sum of a_S over splits separating N and S.
Physical input: only SSA/WM for discs and annuli (Lemma B), i.e. Wall's
maximin theorem, plus SO(2)-invariance of the minimal HRT pieces.

## 4. Scan results (runs/SUMMARY.txt, runs/v_*.json)

Searched: Vaidya-AdS4 Poincare, M=1 (z_h=1), shell profile tanh with width
dlt = 0.02 (near-thin), 0.1, 0.3; boundary slices t = t0 + kappa (r - r_min),
t0 in {0, 0.3, 0.6, 1.0, 1.5, 2.5}, kappa in {0, +0.6, -0.6} (31 of the 54
combinations, every t0 with kappa=0 for dlt=0.02 and 0.3, tilted slices at
dlt=0.02, 0.1, 0.3); edge times therefore span -0.1 .. 3.8, i.e. pre-shell,
during collapse, and after saturation.  Each run: 8 random edges with radii
0.25-2.4 (ratios >= 1.07), all 8 caps, all 56 directed tube shootings,
169 branches added by propagation, 400 random sub-pools of 5-9 edges x 10
random assignments of the annuli to A..E,O: 4000 entropy vectors per run,
124000 in total, all 372 facets evaluated.

  orbit        min margin over 124000 vectors
  MMI(3)       -2.8e-14
  MMI(comp)    -3.4e-14
  n5 #90       -5.2e-14
  n5 #72       -3.6e-14
  n5 #60a      -9.9e-14
  n5 #60b      -6.8e-14
  n5 #10       -4.3e-14
  (SA not evaluated on renormalised areas; it holds trivially with the 2 pi rho/eps terms)
  min K-slack (SSA/WM instances of Lemma B): -4.4e-16

Error budget: renormalised areas reproducible to <= 2.2e-7 (two-end
mismatch), so margins are certain to about 1e-5; the observed margins are
exact zeros up to round-off (the facet is tight on many matching vectors), as
Lemma A predicts.  Caps: one branch in every case sampled (no swallow-tail for
R <= 2.4 at these times); no edge lacked a cap.  13 dominant tubes were found
from one end only (all on tilted slices, near the existence fold); since a
miss of a lower branch from both ends would show up as a K-slack (SSA)
violation and none occurred after the fixes, the data are consistent.
Controls: static pre-shell pool (all margins >= -7e-13), negative-energy
(NEC-violating) shell at t0 = 0.5, 1.0 (no SSA violation for these annuli;
facets hold, consistent with Lemma A, which needs only SSA).

## 5. Verdict

NO VIOLATION in the searched space, and none is possible there: for
rotationally symmetric (concentric) configurations the covariant vector is
exactly a min-cut vector of one planar graph (Sections 1, 3), using only
SSA/WM of discs and annuli (Wall) and SO(2)-invariance of the minimal pieces.
The same holds for polar caps/bands in global AdS4, any number of shells,
two-sided geometries (two free arcs), and nested non-round chains of curves
whose HRT pieces have <= 2 boundary components.  Together with the strips
theorem this exhausts every symmetry-reduced (curve) problem in Vaidya-AdS4.
A 3+1 violation, if it exists, requires HRT surfaces with >= 3 boundary
components or tree-like (non-chain) laminar families: genuinely 2D extremal
surfaces with no continuous symmetry.  That needs a Lorentzian 2D extremal-
surface solver with topology enumeration; not attempted here.
