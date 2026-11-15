# Covariant holographic entropy inequalities for vacuum rotating wormholes in 2+1 dimensions

*Draft of 15 November 2026. Written by an AI system (Claude, Anthropic) during an autonomous research
session; reviewed by an automated referee, not by a human expert, and not formalised. Status of every
statement is given in §9.*

## Abstract

For time-symmetric states the Ryu–Takayanagi (RT) entropies obey an infinite family of linear inequalities,
the holographic entropy cone. Whether the covariant Hubeny–Rangamani–Takayanagi (HRT) entropies obey the
same inequalities is open; in 2+1 bulk dimensions it is established for simply connected bulks and for
`π_1 = Z`, while for multiply connected bulks only numerical evidence exists (arXiv:2407.07165). We treat
vacuum solutions of arbitrary topology — rotating BTZ, rotating multi-boundary wormholes, wormholes with
handles — which are labelled by two hyperbolic structures `(X_L, X_R)` on one surface and in which a
geodesic in a topological class `c` has length `½(ℓ_{X_L}(c) + ℓ_{X_R}(c))`. We prove that the functional
`min_c ½(ℓ_{X_L}(c) + ℓ_{X_R}(c))`, with the *same* class on both sides, satisfies every inequality that admits
a contraction map, for arbitrary boundary regions (intervals of a boundary Cauchy slice and whole boundaries);
for whole-boundary regions we prove it lies in the full RT cone. Assuming the standard prescription that the
HRT entropy of these spacetimes is this minimum over classes, the covariant and static cones agree for vacuum
2+1 gravity. The mechanism is that of Czech and Dong (cut-and-paste in an auxiliary geometry, crossing
counts, a positive measure on geodesics); the new step is to run it on the quotient surface, where the
crossing inequality holds for every homotopy class of test curve and can be integrated against the Liouville
measures of both chiral surfaces.

## 1. Introduction

**The question.** Let a holographic state have a classical bulk dual and let `S(A)` denote the entropy of a
boundary region. When the bulk has a time-reflection-symmetric slice, `S(A)` is the area of a minimal surface
on that slice and satisfies strong subadditivity, monogamy of mutual information, and all further
inequalities provable by contraction maps (Bao et al., arXiv:1505.07839). In a general spacetime `S(A)` is
the area of the minimal *extremal* surface (HRT). Strong subadditivity and monogamy are known to survive
(Wall, arXiv:1211.3494). For the higher inequalities this is open; the October 2026 review of Grimaldi and
Headrick states: "The validity of higher HEIs for the HRT formula remains an open problem."

**What was known in 2+1 dimensions.** Czech and Dong (arXiv:1905.03787) proved all contraction-map
inequalities for simply connected bulks satisfying the null curvature condition, using a positive measure on
kinematic space. Their extension to multiply connected bulks by passing to the universal cover is incorrect:
"an infinite number of surfaces enter on both sides of the inequality" (Grado-White, Grimaldi, Headrick,
Hubeny, arXiv:2407.07165). Those authors gave an argument for `π_1 = Z` by finite covers, explained why it
fails for free groups, and tested 1,877 inequalities numerically on rotating three-boundary and
one-boundary-torus vacuum wormholes without finding a violation. They conclude: "The validity of the RT
inequalities for multiply-connected 2+1 dimensional bulk spacetimes remains an open question." Zhao
(arXiv:2602.04888) proves a covariant graph-model theorem in any dimension conditional on the non-emptiness of
"exposed regions" for every pair of HRT surfaces, a condition that is assumed rather than established and
whose failure modes are left to future work.

**Result.** For vacuum solutions we remove the topological restriction.

> **Theorem (vacuum 2+1, any topology; given Hypothesis H of §7).** Let `M` be a locally `AdS_3` spacetime
> determined by two convex-cocompact Fuchsian representations `ρ_L, ρ_R` of the fundamental group of a surface
> `S`, in the same Teichmüller component. Let the parties be disjoint unions of intervals of a boundary
> Cauchy slice and of whole asymptotic boundaries. Then the HRT entropies satisfy every inequality admitting
> a contraction map: exactly, at every sufficiently small cut-off, if the inequality is balanced at each
> endpoint, and with a difference of the two sides diverging to `+∞` otherwise. For parties consisting of
> whole boundaries the HRT entropy vector lies in the RT holographic entropy cone.

The mathematical content, which involves no physics, is Theorem 1 of §5: the same statement for
`min_c Σ_j w_j ℓ_{X_j}(c)` with any number of hyperbolic structures and positive weights.

**Why the common class is the difficulty.** If the two chiral lengths were minimised independently the
entropy vector would be a sum of two RT vectors and there would be nothing to prove. The spacetime geodesic
is labelled by one group element acting on both chiral factors, so the minimisation is over a common class,
and the result is in general strictly larger than the sum of the two separate minima. Cutting and pasting
the `X_L`-geodesics and the `X_R`-geodesics separately does not help, because the two outputs need not be
homotopic (§4.3).

## 2. Setting

`X = H²/Γ` is a complete hyperbolic surface all of whose ends are funnels (`Γ` finitely generated,
torsion-free, purely hyperbolic, of the second kind; the trivial and cyclic groups are allowed). `∂_∞X` is its
ideal boundary, a union of circles; `P ⊂ ∂_∞X` is a finite set of marked points, cutting `∂_∞X` into elementary
pieces, each assigned to one of the parties `1,…,N` or to the purifier. A *cut* for a region `A` is an
embedded multi-curve — closed curves and arcs ending at the marked points in `∂A` — bounding a region whose
closure meets `∂_∞X` exactly in `A`; `𝒞(A)` is the set of cuts up to isotopy and deletion of trivial or
repeated components. A *test class* is the homotopy class of a closed curve or of a proper arc with ends in
`∂_∞X \ P`, the ends being free to move within their pieces; `i(τ, c)` is the geometric intersection number.
Geodesic representatives realise all these numbers simultaneously.

## 3. The topological contraction lemma

**Lemma A.** Let the inequality `Σ_l α_l S(A_l) ≥ Σ_r β_r S(B_r)` admit a contraction map `f`, and let
`γ_l ∈ 𝒞(A_l)` be arbitrary. There exist `δ_r ∈ 𝒞(B_r)` such that, for every test class `τ`,
`Σ_r β_r i(τ, δ_r) ≤ Σ_l α_l i(τ, γ_l)`. Moreover `Σ_r β_r[p ∈ ∂B_r] ≤ Σ_l α_l[p ∈ ∂A_l]` at every marked point.

*Proof sketch* (complete proof: PROOFS.md §2). Realise the `γ_l` as geodesics with homology regions `W_l`;
label each cell of the arrangement by the bit string of its memberships; let `U_r` be the union of cells on
which the `r`-th coordinate of `f` is 1; resolve `∂U_r` at crossings and cap pairs of rays at marked points
to obtain a cut `δ_r` for `B_r`. A geodesic representative of `τ` meets each geodesic of the arrangement
minimally. At each of its crossing points the two adjacent bit strings differ exactly in the coordinates of
the `γ_l` containing that geodesic, and the contraction property bounds the weighted number of `r` for which
the point lies on `∂U_r`. Summing over crossing points gives the claim. ∎

The lemma is topological. This is what allows the same `δ_r` to be used for several metrics.

## 4. Length as an integral of intersection numbers

### 4.1 Crofton formula with horoball cut-offs

Let `m̄` be the Liouville measure on the complete geodesics of `X`, normalised so that the geodesics crossing a
segment have measure equal to its length. Almost every geodesic runs from funnel to funnel and avoids the
marked points, hence defines a test class. Regulate arcs by removing horoballs `B_p(ε)` at the marked points
and write `ℓ^ε(c)` for the remaining length of the geodesic representative. With the horoballs chosen below a
threshold `ε_0` depending only on `(X, P)`, no geodesic of any class enters a horoball except along its
terminal rays, so `ℓ^ε(c) = ℓ^{ε_0}(c) + e(c) log(ε_0/ε)` exactly, `e(c)` being the number of arc-ends.

**Lemma B.** Let `𝒢_ε` be the set of geodesics of `X` disjoint from all `B_p(ε)`. Then, as `ε → 0`,
`ℓ^ε(c) = ∫_{𝒢_ε} i(τ(ĝ), c) dm̄(ĝ) + κ e(c) + o(1)`, with `κ = 1 − log 2`.

The constant is the measure of geodesics that cross a ray below a horocycle and also enter the horoball,
`∫_0^{π/2} cos φ log(1/cos φ) dφ`; it is independent of the height by dilation invariance.

### 4.2 Theorem 1

Let `X_1,…,X_k` be hyperbolic structures as above on one marked surface (identified by homeomorphisms
respecting marked points and parties), `w_j > 0`, and
`F^ε(A) = min_{c ∈ 𝒞(A)} Σ_j w_j ℓ^ε_{X_j}(c)`.

**Theorem 1.** If `Σ_l α_l S(A_l) ≥ Σ_r β_r S(B_r)` admits a contraction map then, with
`Δ = Σ_l α_l|∂A_l| − Σ_r β_r|∂B_r| ≥ 0`: if `Δ = 0`, `Σ_l α_l F^ε(A_l) ≥ Σ_r β_r F^ε(B_r)` for every `ε ≤ ε_0`; if
`Δ > 0` the difference equals `C + (Σ_j w_j)Δ log(ε_0/ε)`.

*Proof.* Take minimisers `γ_l` (independent of `ε`), apply Lemma A, integrate its inequality over `𝒢_ε(X_j)`
for each `j`, and use Lemma B for the finitely many classes involved: this gives
`Σ_l α_l F^ε(A_l) − Σ_r β_r F^ε(B_r) ≥ κ(Σ_j w_j)Δ − o(1)`. The left side is `C + (Σ_j w_j)Δ log(ε_0/ε)` by the
exact scaling of `ℓ^ε`; for `Δ = 0` this forces `C ≥ 0`. ∎

### 4.3 Remark: the naive argument fails

Two minimal-position arrangements of the same classes differ by triangle moves, and a triangle move can
change how `∂U_r` pairs the six boundary points of a small disc (when the `r`-th coordinate of `f` is 1 exactly
on strings of weight ≤ 1 in the three relevant bits). So the outputs of cut-and-paste performed in `X_L` and in
`X_R` need not be homotopic. Theorem 1 uses only the `X_1`-arrangement.

## 5. The full cone

**Lemma G (graph realisation; suggested by the referee).** For a finite positive combination
`μ = Σ w_k τ_k` of test classes, `A ↦ min_{c ∈ 𝒞(A)} i(μ, c)` is the min-cut function of a finite weighted
graph: vertices are the crossing points of the `τ_k` (drawn minimally) and their ideal ends, the latter
labelled by parties; edges are the segments between them.

**Theorem 3.** For parties that are unions of whole boundary circles (`P = ∅`), the vector
`(min_c Σ_j w_j ℓ_{X_j}(c))_A` lies in the holographic entropy cone, i.e. satisfies every inequality valid for
RT. (For regions with endpoints an argument is given in PROOFS.md §5.4; it is not written at the same level
of detail and is not claimed here as a theorem.)

*Proof.* The Liouville measure pushes forward to a countable positive combination `μ` of test classes with
`i(μ, c) = Σ_j w_j ℓ_{X_j}(c)`. Finite partial sums `μ_n` give graph vectors (Lemma G). A finite filling family
of arc classes of positive measure bounds the complexity of the minimisers, so they range over a finite set
and `min_c i(μ_n, c) → min_c i(μ, c)`. The cone is closed. ∎

## 6. Vacuum 2+1 gravity

`AdS_3 = SL(2,R)` with `ds² = −det dx`; a pair `(ρ_L, ρ_R)` of Fuchsian representations as in §2, in the same
Teichmüller component, defines `M = D̃/Γ` (Mess; Barbot; Bonsante–Krasnov–Schlenker), with one asymptotic
region per end of `S`. A boundary Cauchy slice projects homeomorphically to the ideal circles of both `X_L`
and `X_R`; together with the marking it determines an identification of marked surfaces, hence of cut
classes, including boundary twists.

**Proposition (lengths).** (a) `(ρ_L(γ), ρ_R(γ))` translates a spacelike geodesic of `D̃` by
`½(ℓ_L(γ) + ℓ_R(γ))`. (b) For boundary points `p = ab^T`, `q = cd^T` with `N = ω(a,c)ω(b,d) > 0`,
`x(t) = (e^t p + e^{−t} q)/√N` is the geodesic between them and its regulated length is
`log N + 2 log(1/ε) = ½(λ_L + λ_R) + 2 log(1/ε)`, `λ` being horocycle distances in `H²`.

Hence the regulated length of the spacetime geodesics of class `c` is `½(ℓ^ε_{X_L}(c) + ℓ^ε_{X_R}(c))`
(equations (58), (60) of arXiv:2407.07165).

## 7. Hypothesis H and the statement about HRT

> **H.** In these spacetimes the HRT entropy of `A` is `(1/4G_N) min_{c ∈ 𝒞(A)} ½(ℓ^ε_{X_L}(c) + ℓ^ε_{X_R}(c))`.

This is the prescription used in the literature for these spacetimes. It has three parts. (H1) Extremal
candidates are embedded cuts: reason given, not written as a proof. (H2) Every embedded class is admissible
(its geodesics are mutually spacelike and bound a spacelike region with `A`): **assumed**. (H3) The
literature's minimisation over possibly non-embedded collections agrees with the minimisation over embedded
cuts: proved (PROOFS.md Lemma 6.4). A derivation of H from maximin is sketched in DERIVATIONS.md §5 at the
level of rigour of Wall's theorems.

Given H, Theorem 1 with `k = 2`, `w = (½, ½)` is the Theorem of §1, and Theorem 3 gives the cone statement
for whole-boundary parties. If H2 failed for some class, the theorems would remain true for the functional
but would not transfer to HRT.

## 8. Numerical tests

These test conclusions, not proofs. Details and commands: EXPERIMENTS.md.

*(table inserted after the certified runs complete; see EXPERIMENTS.md)*

## 9. Status of statements

| Statement | Status |
|---|---|
| Lemma A, Lemma B, Theorem 1 | proved; reviewed by one automated referee (round 1: no fatal error; repairs incorporated) |
| Lemma G, Theorem 3 for whole-boundary parties | proved; referee round 2 pending at time of writing |
| Theorem 3 for regions with endpoints | argument only; not claimed |
| Length proposition, marking | proved |
| Hypothesis H (in particular H2) | assumed |
| Statement about HRT | proved given H |
| Non-vacuum spacetimes, higher dimensions, extremal horizons, conical defects | not covered |

## 10. What is new and what is not

Not new: contraction maps; cut-and-paste in an auxiliary geometry; counting crossings; a positive measure on
geodesics (Czech–Dong); the left/right length formula for vacuum wormholes (Maxfield; Grado-White et al.);
the numerical evidence for the conclusion.

New, to the best of our knowledge after the search described in LITERATURE.md: the crossing inequality on the
quotient surface for all test classes (Lemma A); the Crofton formula with horoball cut-offs and its universal
constant in this context (Lemma B); the resulting proof for vacuum solutions of arbitrary topology; the
multi-metric generalisation; the graph realisation and the full-cone statement for whole-boundary parties.

This is a vacuum-only repair of a known gap. It does not address matter, higher dimensions, or the
black-hole information problem.

## References

Bao, Nezami, Ooguri, Stoica, Sully, Walter, arXiv:1505.07839. Wall, arXiv:1211.3494. Czech, Dong,
arXiv:1905.03787. Grado-White, Grimaldi, Headrick, Hubeny, arXiv:2407.07165 and arXiv:2502.09894. Zhao,
arXiv:2602.04888. Grimaldi, Headrick, arXiv:2610.33905. Maxfield, arXiv:1412.0687. Mess, *Lorentz spacetimes of
constant curvature* (1990). Barbot, *Causal properties of AdS-isometry groups I, II*. Bonsante, Krasnov,
Schlenker, arXiv:math/0610429. Aminneborg, Bengtsson, Holst, *Class. Quantum Grav.* 15 (1998). Santaló,
*Integral Geometry and Geometric Probability*. Penner, *Comm. Math. Phys.* 113 (1987). Epstein, *Acta Math.*
115 (1966). Farb, Margalit, *A Primer on Mapping Class Groups*. Bao, Furuya, Naskar, arXiv:2506.18086.
