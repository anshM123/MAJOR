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

## 3. The thimble integrals: no sign problem, but a rugged Jacobian; common random numbers [NE]

* **No sign problem.** On the `A±` thimbles (n3_hybrid contour built on `A±` with `ν_f = ν`, `ρ₁ = 0.8`) the weights have nearly constant
  phase: `|Σw| / Σ|w| ≈ 0.99` at `ν = 5` (`Δ = 100`, `J = 0, 4`). On the old `ν_f = 13` contour the same integral is a 50–100 %
  cancellation (`ν = 5.45`: 56 % RQMC error at 2¹³ points); the generalised-thimble step (`N3_TFLOW = 0.3`) does not help (62 %, 25 % of
  the points cut).
* **But the thimble Jacobian is rugged.** `|det ∂ξ/∂ζ|` of the flowed part varies by factors 3–10 over `δζ ≈ 0.02–0.1` (checked: RK4
  tangent propagation = finite differences of the parametrisation, `76.05` vs `76.10`). The flow lines separate along ridges between
  sinks (gradient `|∂Φ|` stays ≈ 1.5–2 on one line and rises to 30 on a neighbour 0.05 away); the spikes are not near a critical
  point and not near `A = 0` (`|A| ≈ 2` on the spike line). Single samples carry up to 10 % of the integral at 2¹⁰ points.
  Consequence: relative RQMC errors 3–25 % per thimble at 2¹³ points (`ν = 5–10`), not improvable by Gaussian proposal adaptation.
* **Cut lines carry negligible weight.** Lines cut by the per-step winding test (`WCUT`, 1–15 % of the points) were re-done with a
  winding-controlled adaptive flow (`n4_thimble._flow_rho_adaptive`, steps halved until every log moves < 0.4 rad): same mean to 7
  digits at 40× the cost. So the default (n3) flow is kept; the graviton log replay is sub-stepped where a segment winds > 0.4 rad
  (`grav_on_paths_sub`).
* **Common random numbers (CRN) rescue the ratio.** The ruggedness is a property of the geometry of the soft block, which is nearly the
  same for the state (`Δ ≈ 100–130`) and the reference (`Δ = 2·10⁴`) at equal `ν` (in the normalised `ξ` coordinates). With the same
  proposal and Sobol seeds, at `ν = 5`, `J = 4`, `A+`, fixed `N(0,1)` proposal (deliberately bad: one replica 30× the others): per-thimble
  relative errors 77 % / 86 %, state/reference ratio 5.8 % (and 2 % without the outlier replica). The production runs therefore compute
  state and reference on common random numbers (`n4_run.run_pair`: proposal adapted on the state, re-used for the reference).

## 4. Scalar validation against the exact Mellin object (`n4_validate.py` → `cache/n4_validate.txt`) [NE]

**Comparison object.** For every block `D` the exact `X_D(ν) = g C(ν+i) + g̃ C(ν−i)` (n3_exact: closed-form part at complex `ν` +
`ρ₂` continued from the exact real-`ν` nodes) is computed from e1 records, or — for blocks not in e1 (the `Δ = 100` row, the extra
rows at the chain's exact `Δ`, and `Δ_ref = 2·10⁴`) — from new exact records `n4_exact.py` (`exactskel.skelexact.ExactState`,
read-only import; reproduces e1's `ρ₂` to 10⁻¹⁵ and `X` to 10⁻¹¹ on an e1 state; ≈ 1–4 s per block).
`d_D = log Y_D − log(q_D X_D)` must be independent of `D` at fixed `(ν, J)` (CLS-D normalisation and the `J`-only constant of `Q_J`).
This holds **as a complex number** in the n3 branch convention: on the validated hybrid points (`ν = 9.197`, `J = 4`, `m/Λ =
1.05–1.3`) `Re d` is constant to 5·10⁻⁴ and `Im d` to 2·10⁻⁴. The test statistic is `d_D − d_ref` with `D_ref = 2·10⁴` — exactly the
state/reference ratio used by the graviton chain — in c2 units `Δc2 = m_D² Re(d_D − d_ref) ex₂/E₂`; errors by jackknife over the 8
common-seed replicas (heavy tails: the jackknife is not reliable at the 2σ level).

| `ν` (`J = 4`) | blocks | `A₊ + A₋` (MAIN): `Re(d − d_ref)`, `Im(d − d_ref)` | `Δc2` | `A₊` alone / `A₋` alone (`Δc2`) |
|---|---|---|---|---|
| 9.200 | 105.0, 109.9, 119.8, 129.6 vs 2·10⁴ | `−0.0018 ± 0.0056`, `+0.027 ± 0.021` (1.05Λ) … `−0.0011`, `+0.017` (1.3Λ) | **−16 ± 51** (all four) | −59 … −72 (±7·10³) / −127 … −146 (±85); `Im` off by ∓0.08, ±0.10 |
| 5.451 | same | `+0.045 ± 0.025`, `−0.106 ± 0.023` (1.05Λ) … `+0.026 ± 0.015`, `−0.064 ± 0.014` (1.3Λ) | **+453 ± 255 … +403 ± 237** | +2230 (±46) / +2303 (±30) |

* At `ν = 9.2` the thimble sum `A₊ + A₋` reproduces the exact state/reference ratio (modulus to 0.2 %, `Δc2 = −16 ± 51`), while
  either thimble alone is wrong (phase off by ±0.1 rad, modulus by 1 %): **this confirms `n₊ = n₋ = 1` below the fold, in the
  tracked branch convention** [NE].
* At `ν = 5.45` the residual is `Δc2 ≈ +430 ± 250` (1.8σ) in the modulus and `−0.064 … −0.106` rad (4.7σ by the jackknife) in the phase;
  all four states show the same pattern `(+470 − 1100i)/m²` (common seeds make the statistical error itself ∝ `1/m² − 1/m_ref²`, so
  this pattern does not distinguish a fluctuation from a systematic). Either thimble alone is off by `Δc2 ≈ +2200`.

* **Odd `J`: no scalar test.** For odd `J` the scalar `n0 = +1` crystal integral is not the identical-scalar functional (whose
  `S⁰ = {n0 = ±1}` sum vanishes for odd `J`), so `d_D − d_ref` is meaningless there (it comes out O(1)); the odd production points
  carry no exact scalar check. The odd graviton integrand (P8.build, `n0 = +1`) is the chain's object as before.

## 5. Production (`n4_run.py`, `n4_prod.py`; caches `n4_raw.pkl`, `n4_pts_<kind>_J<J>.pkl`, `n4_odd_pts_<kind>_J<J>.pkl`) [NE]

* 17 state/reference pairs at the `x = 0.05` node (`ν = 0.05 m`): even `J = 4, 6, 8, 10` on the `Δ = 100` row (`C2`, `C3H`); odd
  `J = 5, 7, 9` on the `Δ = 100` row and on the extra rows (p8_k2 keys: `J = 5` at 1.05/1.1/1.2/1.3Λ, `J = 7, 9` at 1.1/1.2/1.3Λ; at
  1.05Λ only `J = 5` exists above the gap), each with its `Δ_ref = 2·10⁴` partner at the same `ν`.
* Per block: saddles tracked from `ν_c = 10` (the `2·10⁴` track is shared by all `ν` of a `J`); thimbles `A±` (`ρ₁ = 0.8`), 2⁹·2·8 =
  8192 RQMC points each, proposal adapted to the state's scalar + graviton integrands and re-used for the reference (CRN), graviton
  logs replayed along the flow paths (sub-stepped where needed). ≈ 23–35 min per pair (1 core).
* **Stored raw** (for recombination with the non-local `χ` later): `n4_raw.pkl[(mode, Δ, J, ν)]` holds, per thimble, the literal
  crystal integrals `Ψ_ν` (unshifted, at the real node `ν`) of the scalar and of each graviton kind: replica values, means, `logN`
  (graviton: n1_run convention `logN = Φ_g(saddle) + log det L`, log-Jacobian constant included), proposal, diagnostics.
  `n4_(odd_)pts_*` hold the sums `A₊ + A₋` in n1_run / p8 format (key `(Δ, round(ν, 9))`) for the chains. No shift is applied here;
  the chains apply their own (local) de-shift.
