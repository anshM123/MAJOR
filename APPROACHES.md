# APPROACHES

Target: contraction-map inequalities for HRT in vacuum 2+1 gravity of arbitrary topology (PROBLEM.md).

The structural fact used throughout: in a vacuum solution labelled by two hyperbolic structures `(X_L, X_R)` on
the same surface, the spacetime geodesic in a topological class `c` has length `½(ℓ_L(c) + ℓ_R(c))`, with the
**same** class on both sides. If the two minimisations were independent the result would be trivial (a sum of
two RT entropy vectors lies in the cone). The whole difficulty is the common class.

## Approach 1 — Crofton / intersection numbers. **Succeeded.** (PROOFS.md §§2–6)

Idea: do not compare lengths of cut-and-paste pieces at all. Perform the cut-and-paste once, on one auxiliary
hyperbolic structure, and record only the *topological classes* `δ_r` of the outputs. Show the purely
topological inequality

`Σ_r β_r i(τ, δ_r) ≤ Σ_l α_l i(τ, γ_l)` for every test curve or arc `τ`

(Lemma A), by counting crossings of a geodesic representative of `τ` with the arrangement. Then integrate
over `τ` against the Liouville measure of `X_L` and of `X_R` separately (Crofton formula, Lemma B). Because
hyperbolic length is the Liouville integral of intersection numbers, both chiral lengths obey the inequality
with the same classes `δ_r`.

What made it work: the output classes are metric-independent, and the inequality is linear in the measure.

## Approach 2 — keep one class by showing cut-and-paste is metric-independent. **Failed; obstruction identified.**

Idea: take geodesic representatives separately in `X_L` and in `X_R`, cut and paste in each, and prove the two
outputs are homotopic. Two minimal-position arrangements of the same classes differ by isotopy and triangle
(Reidemeister III) moves, so it would suffice that a triangle move does not change the homotopy class of
`∂U_r`.

Obstruction (DERIVATIONS.md §4): restrict the coordinate function `m = f(·)_r` to the three bits flipped
around a triangle. If `m = 1` exactly on strings of weight ≤ 1, then before the move `∂U_r` pairs the six
boundary points of a small disc around the triangle "around the weight-2 sectors", and after the move
"around the weight-1 sectors". These are different planar matchings, and nothing in the contraction property
forbids this pattern for a single coordinate `r`. So the homotopy class of the output can depend on the
arrangement. Approach 1 avoids the issue by never needing the `X_R`-arrangement.

## Approach 3 — covering spaces (Grado-White–Grimaldi–Headrick–Hubeny). **Fails for non-abelian π₁ (their result).**

Pass to finite covers in which the relevant geodesics become minimal; works for π₁ = ℤ because the Cayley
graph is a line; fails for free groups because the fraction of the cover near the cut does not vanish
("because of the hyperbolicity of G"). Not pursued.

## Approach 4 — geodesic currents (Bonahon). **Equivalent to Approach 1 for closed surfaces; not used.**

`ℓ_X(c) = i(L_X, c)` with `L_X` the Liouville current, and closed geodesics are dense in the space of currents,
so an inequality valid for every closed test curve extends to `L_{X_L} + L_{X_R}`. For surfaces with funnels
the Liouville measure is carried by geodesics that exit through the funnels and is *not* a limit of closed
geodesics, and boundary-anchored arcs need a cut-off. Approach 1 is the direct version that handles both
(test arcs with ends in the ideal boundary; horoball cut-offs with the universal deficit `κ = 1 − log 2`).

## Approach 5 — beyond vacuum: a Crofton measure for a general 2+1 spacetime. **Open; reduction only.** (OPEN_GAPS.md)

Lemma A needs only that candidate lengths are integrals of intersection numbers against a positive measure.
In a general 2+1 spacetime one can try to define that measure on pairs of boundary points of the universal
cover as the mixed second difference of anchored geodesic lengths. Positivity is exactly a crossing
inequality ("strong subadditivity on the universal cover"). This would reduce the general 2+1 problem to that
statement plus a statement about closed geodesics. Not proved.

## Secondary target — Gaussianity of the one-time-point SYK partition function. **Not started.**

`z` is a homogeneous polynomial of degree `p = N/q` in independent Gaussians, i.e. an element of the `p`-th
Wiener chaos. Planned steps: (i) product formula ⇒ `z² − E z² − :z²:` is a sum over partial contractions;
(ii) its `L²` norm is bounded by the fourth cumulant; (iii) Nourdin–Peccati ⇒ total-variation distance to a
Gaussian `≤ (2/√3)·sqrt(κ₄/σ⁴)`; (iv) rigorous asymptotics of the exact fourth moment (eq. 2.23 of
arXiv:2103.16754) to show `κ₄/σ⁴ → 0` for `q > 2`.
