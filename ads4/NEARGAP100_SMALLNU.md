# The `x = 0.05` node at `Δ_g = 100` (`ν ≈ 4.5–9`): contour evaluation below the fold (IN PROGRESS)

*Written 19 Nov 2026 by the numerics agent continuing H0_SPINNING_FUNCTIONAL.md §25 (hand-over §25.9 item 1) and
H7_I14_LIGHTSIDE.md §23.4/§23.6 item 2. Unrefereed, not committed. Updated incrementally; hand-over at the end.*

Grades: [P] proved, [D] derived, [NE] numerical evidence, [M] model, [I] open.

## 1. Saddle structure below the fold, tracked in `ν` (`code/spinfunc/n4_thimble.py`, `track`) [NE]

* The pair `A±` of critical points of the integrand's own scalar exponent `Φ_ν` (`n1_direct.phi_w`) is obtained by a
  symmetry-breaking kick of the symmetric saddle just below the fold, then continued in `ν` (Newton, `Δν ≤ 0.25`). All logs (scalar
  flow function, scalar integrand, graviton integrands) are continued along the `z`-path traced by the saddle, starting from the
  n3 convention (`W_REF` next to the crystal of the calibration `ν`, straight `w`-path). In this tracked convention the
  intersection numbers can only jump at a Stokes point.
* `Δ = 100`, `J = 4`: `Re Φ(A+) = Re Φ(A−)` to all digits at every `ν` (`−28.83 … −35.94` for `ν = 10 … 5`). `Im Φ(A+) − Im Φ(A−)`
  grows from ≈ 0 at the fold to 1.45 / 2.17 / 2.97 / 3.84 at `ν = 8 / 7 / 6 / 5` (two saddles of equal weight with a relative phase
  of < 1 period: the oscillatory side of an Airy-type transition).
* Below the fold the saddles are far from the crystal in `w` (`|w| ≈ 1.2` at `ν = 5`), soft `−Hess` eigenvalues 1.8–13
  (`|λ|`), two hard ones (≈ 1.5·10³, 4.1·10³).

## 2. Intersection numbers at `ν = 10` [NE, weak]

`Δ = 100`, `J = 4`, scalar: thimble sums (2¹¹ points each, 12 % / 25 % RQMC errors) vs the n3 hybrid contour (`ν_f = 13`, 1 %):
`|Y₊ + Y₋ − Y_H|/|Y_H| = 0.096`; all other combinations with `n± ∈ {−1, 0, 1}` are off by ≥ 0.50. So `n₊ = n₋ = 1` (consistent
within the thimble errors; not a precision test).

## 3. Why the thimbles of `A±` do not give a usable estimator (so far) [NE]

(see §3 below; being written)
