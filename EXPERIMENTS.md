# EXPERIMENTS

All numerical work is evidence only. None of it is used in a proof. Commands are run from `code/`.

## E1. Contraction-map certificates (`contraction.py`)

Exhaustive back-tracking search for a contraction map `f:{0,1}^L→{0,1}^R`, followed by an independent
integer-arithmetic verification of the contraction property on all pairs and of the boundary conditions.

```bash
python3 contraction.py
```

| inequality | L, R | result |
|---|---|---|
| strong subadditivity | 2, 2 | map found, verified |
| monogamy of mutual information (cyclic, k=1) | 3, 4 | map found, verified |
| 5-party cyclic inequality (k=2) | 5, 6 | map found, verified |
| 7-party cyclic inequality (k=3) | 7, 8 | map found, verified |
| reversed monogamy (deliberately false control) | 4, 3 | no map exists |

## E2. Cylinder topology, π₁ = ℤ, boundary intervals (`annulus_test.py`)

Hyperbolic cylinder with `n` marked points on one ideal circle, parties = random unions of the `n` arcs, the
other end and at least one arc assigned to the purifier. Cut classes are enumerated exactly (non-crossing
matchings × position of the hole). `k` independent metrics (core length log-uniform in [0.3, 12], random
point positions), random positive weights, the same class used in all metrics.

```bash
python3 annulus_test.py 1 3000 25.0     # seed, trials, log(1/eps) per arc end
```

| inequality | k=1 | k=2 | k=3 |
|---|---|---|---|
| strong subadditivity | 0 / 3000 | 0 / 3000 | 0 / 3000 |
| monogamy | 0 / 3000 | 0 / 3000 | 0 / 3000 |
| 5-party cyclic | 0 / 3000 | 0 / 3000 | 0 / 3000 |
| 7-party cyclic (seed 0, finite parts) | 0 / 750 | 0 / 750 | 0 / 750 |

(entries: violations / trials; the smallest value of LHS − RHS is −2·10⁻¹³, i.e. rounding.)

Controls, same sampler:

| functional / inequality | violations |
|---|---|
| `min_c max_j ℓ_j(c)` (not of Crofton type), monogamy, k=2 | 923 / 3000 |
| `min_c max_j ℓ_j(c)`, 5-party cyclic, k=2 | 559 / 3000 |
| reversed monogamy, sum functional, k=2 | 2677 / 3000 |

**A mistake caught here.** The first run dropped the cut-off term and compared finite parts only; it reported
1830 / 3000 "violations" of strong subadditivity even for a single metric. That was a regularisation artefact:
strong subadditivity is not balanced at an endpoint where pieces of `A` and `C` touch, so the true difference
diverges to `+∞` there and the finite parts carry no sign. With the cut-off included (`log(1/ε) = 25` per arc
end) there are no violations. Inequalities that are balanced at every endpoint (monogamy, cyclic) are
insensitive to this.

## E3. Four-holed sphere, π₁ = F₃, whole-boundary regions (`fourholed_test.py`)

This is the first multiply connected case with free non-abelian fundamental group, and with whole-boundary
regions the inequality of Theorem 1 is exact. `k = 2`, weights (½, ½) is the HRT entropy of a rotating
four-boundary vacuum wormhole.

* Hyperbolic structures: Fenchel–Nielsen construction (two pairs of pants glued along a neck of length
  `l_e ∈ [0.3, 4]` log-uniform, twist uniform in `[−2 l_e, 2 l_e]`, boundary lengths uniform in `[1, 6]`).
* Essential simple closed curves: all `x₁x₂`, `x₂x₃` for braid-group images of the standard generators up to
  braid word length 7; every matrix carries a running bound on its rounding error, and a curve whose trace is
  not known to `10⁻⁷` is discarded.
* Cut classes for a set `T` of boundaries: a set of boundary geodesics plus at most one essential curve.

```bash
python3 fourholed_test.py 3 600 7 fn    # seed, trials, braid depth, sampler
```

| functional | violations of monogamy / strong subadditivity | trials in which an essential curve is part of a minimal cut |
|---|---|---|
| sum, k=1 (Ryu–Takayanagi) | 0 / 600 | 99% |
| sum, k=2 (HRT, vacuum wormhole) | 0 / 600 | 100% |
| sum, k=3 | 0 / 600 | 100% |
| max, k=2 (control) | 299 / 600 | 100% |

Stability: for 150 two-metric samples the minimal cuts found with braid depth 4 and depth 7 agree to
`4·10⁻¹³`, so the enumeration depth is not truncating the minimum.

**Two mistakes caught here.**
1. A first version, using randomly generated Schottky groups, reported zero violations but an essential curve
   entered a minimal cut in at most 7% of trials: the test was nearly vacuous (boundary geodesics alone always
   satisfy the inequalities). It also mishandled overflowed matrices. Those numbers are not reported.
2. The first Fenchel–Nielsen gluing put the second pair of pants on the wrong side of the neck: the
   representation was not discrete (an elliptic element with trace 1.195 appeared), and the single-metric
   control then "violated" monogamy in 325 / 500 trials. Four gluing conventions were compared; exactly one
   gives hyperbolic traces for every enumerated simple curve and the correct fourth boundary length in 100/100
   samples. All reported numbers use that one.

## E4. Direct falsification test of Lemma A (`lemmaA_test.py`)

The earlier tests check a consequence of the theorem for Liouville measures. This one tests statement (A1)
itself, on the cylinder with marked points, in exact integer arithmetic. Test classes: arcs between any two
elementary pieces winding up to two full turns, and arcs from each piece to the other end. Intersection
numbers are computed by linking of lifts.

* **(L2), the lemma itself.** For random party assignments and **arbitrary** (uniformly random, not minimal)
  cut classes `γ_l ∈ 𝒞(A_l)`, search exhaustively over all tuples `(δ_r) ∈ Π 𝒞(B_r)` for one satisfying
  `Σ_r β_r i(τ,δ_r) ≤ Σ_l α_l i(τ,γ_l)` for **all** test classes simultaneously.
* **(L1).** For every single test class `τ`, check the inequality for `F_τ(A) = min_c i(τ,c)`.

```bash
python3 lemmaA_test.py 0 400
```

| inequality | (L2) tuples with no dominating `δ` | (L1) single-test-class violations |
|---|---|---|
| strong subadditivity | 0 / 2400 | 0 / 23010 |
| monogamy | 0 / 2400 | 0 / 22436 |
| 5-party cyclic | 0 / 2400 | 0 / 42022 |
| 7-party cyclic | 0 / 600 | 0 / 17193 |
| reversed monogamy (control) | **1077 / 2400** | 0 / 22326 |

The control shows that (L2) is a discriminating test (it fails 45% of the time for a false inequality) and
that (L1) is not (a single test class cannot detect the false inequality). Limitation: the cylinder has
`π_1 = Z`; Lemma A is not tested directly on a surface with non-abelian fundamental group.

## E5. Lemma 6.3' (spacelike separation), `spacelike_check.py`

20,000 random configurations of four boundary points with the same cyclic order in both chiral factors: the
closed formula for `min ½ tr(x^{-1}y)` agrees with direct evaluation in all cases; it is `> 1` for all 40,000
unlinked pairings and `≤ 1` for all 20,000 linked pairings.

## E6. Certified four-boundary test (`fourholed_certified.py`)

Addresses three questions about E3: is the representation discrete and faithful, are all homotopy classes
included, is the homology condition enforced.

* **Discrete and faithful, no relations to check.** Each hyperbolic structure is a classical Schottky group:
  six pairwise disjoint intervals and three Möbius maps pairing adjacent ones. By the ping-pong lemma the
  group is free of rank 3 and discrete. The script verifies the disjointness margin and that each generator
  maps the endpoints of its interval to those of the partner. The fourth boundary class is defined as
  `(g_3 g_2 g_1)^{-1}`; the side pairing has four boundary cycles and Euler characteristic −2, so the quotient
  is a four-holed sphere.
* **All classes.** Every closed geodesic is a cyclically reduced word; a word `x_1…x_n` crosses the fundamental
  domain `n` times, the `k`-th time between two specific sides, so its length is at least the sum of the
  corresponding side-to-side distances. A branch-and-bound over all cyclically reduced words with this lower
  bound below a threshold `Θ` therefore visits every closed geodesic, simple or not, of weighted length `≤ Θ`.
  `Θ` is set to the largest class minimum found among short words, so no unvisited word can lower any
  minimum. The lower bound was calibrated: it holds for all 3,918 words of length ≤ 5 with the correct product
  convention and fails for 936–1020 of them with the reversed one.
* **Homology.** The mod-2 homology class of a word is its vector of exponent parities; a cut for a set `T` of
  boundaries is a collection of closed curves with total class `Σ_{i∈T}[∂_i]`. Non-embedded collections never
  beat embedded ones (Lemma 6.4), so `F(T)` is the cheapest decomposition of that class.
* **Time dependence.** For `k ≥ 2` the metrics differ (independent jitter of the interval data, independent
  pairing parameters, and in 60% of cases a relative Dehn twist of the marking), so `X_L ≠ X_R`: no
  time-reflection symmetry, and for this non-abelian group no stationarity.
* **Non-vacuity.** Base geometries are kept only if a non-boundary class enters a minimal cut.

```bash
python3 -u fourholed_certified.py 3 40
```

Results: see the table at the end of this section (filled in when the run completes).

## E7. What has not been tested numerically

* Boundary intervals on surfaces with non-abelian fundamental group (done at scale for `k = 2` in
  arXiv:2407.07165: 1,877 inequalities, no violation; not repeated here).
* Inequalities with more than 7 parties.
* Lemma A on a surface with non-abelian fundamental group (E4 tests it on the cylinder only). The classes `δ_r`
  are found by exhaustive search, not by executing the construction in the proof.
* Theorem 3 for regions with endpoints beyond the inequalities listed.
