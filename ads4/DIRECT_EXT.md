# DIRECT_EXT: larger realised coefficient table and direct sign cells above c = 0.258

Agent report, 26 Nov. Nothing committed. Outputs are in `code/sliver_direct_ext/`. The working copy is in the session scratchpad at `direct_ext/`.

## 1. Larger realised moment table and coefficient vector (done)

* **Pipeline.** The author's `run_full_gauss.py` → `make_coefficients.py` (code/sliver_direct), run with **K = 3300, D = 2600** on 3 shards.
  * The algorithms are unchanged.
  * Two hard-coded *inputs* had to be changed (diffs are saved):
    * `make_coefficients.py`: `K=1254;D=1600` → `K=3300;D=2600`.
    * `log_constants.cpp`: output scale 10^1600 → 10^2600 and MPFR precision 6528 → 10528 bits, which follows the 4D+128 rule of gaussian_fast_adjacent.
  * Why the second change is needed: `log_constants.txt` (G0, E0) is written at a fixed scale of 10^1600. Feeding it unchanged into a D = 2600 run silently gives a wrong A0 (26.73 instead of 20.77), and the A0 assertion catches it.
  * The new constants agree with the old ones after rescaling, to the last unit.
* **How D was chosen.** Radius growth was modelled from the raw CDH error recurrence plus the binom² factor. The model matches the old table: about 10^1208 of growth at k = 1252, against the observed 10^1224. The prediction was radius ≈ 10^(−D+2083+~55) at k = 3300, so D = 2600 was chosen.
* **Results.**
  * `S_3300_2600.txt`: PASS_FULL_INVENTORY, 3564 centres, max width 1.8e70 units (that is, 1.8e-2530).
  * `direct_coefficients.txt`: **3299 pairs (N = 3298), D = 2600**.
  * Max radius: **A: 10^-466.2, G: 10^-543.6**, reached at k = 3298. At k = 1200 the radius is 10^-1349.
  * A0 = 20.77189301533 (assertion passes).
  * This is enough for need_k(c) up to c = 0.56 (K = 3230).
* **Timing.** About 0.7 s per centre at K = 3300 (0.18 s at the old K = 1254). The inventory took 854 s wall on 3 workers; coefficients took 23 s.
* **Consistency checks.**
  * All 2 × 1253 old coefficient intervals intersect the new ones.
  * Max midpoint difference from the old table is 1e-454.
  * Rerunning the author's `direct_sign_tight.py .223 .258` (b = 0, 4) on the new table reproduces the author's exact rational `physical_lower` (9.970877e-6 and 9.415199e-6). The files are `cells/xcheck_*`.
* **Saved files.** `tables/*.xz` and SHA256SUMS_tables.

## 2. New sign cells with the unchanged `direct_sign_triangle.py`/`.cpp`

| c-cell | b | K | (mc, mh), terms | wall time | lower bound | status |
|---|---|---|---|---|---|---|
| [0.258, 0.262] | 0 | 1397 | (280, 1682), 236883 | 42 s | 1.4037e-5 | PASS_CONDITIONAL |
| [0.258, 0.262] | 4 | 1397 | (280, 1683), 237164 | 42 s | 1.1126e-5 | PASS_CONDITIONAL |
| [0.262, 0.264] | 0 | 1405 | (241, 1688), 204611 | 37 s | 1.0939e-5 | PASS_CONDITIONAL |
| [0.262, 0.263] | 4 | 1401 | (224, 1569), 176850 | 35 s | 1.7369e-6 | PASS_CONDITIONAL |
| [0.262, 0.264] | 4 | 1405 | | 37 s | −9.3e-5 | FAIL_INCONCLUSIVE |
| [0.264, 0.266] | 0 | 1413 | | 35 s | −1.3e-4 | FAIL_INCONCLUSIVE |
| [0.264, 0.265] | 0 | 1409 | | 38 s | −4.9e-6 | FAIL_INCONCLUSIVE |
| [0.258, 0.29] | 0, 4 | 1512 | | 147–150 s | center ~1e12, L1 ~1e15 | FAIL_INCONCLUSIVE |
| [0.29, 0.3125] (tight, K = 1460) | 0 | 1460 | | 137 s | center 6.6e32 | FAIL_INCONCLUSIVE |

* **Cells that need the new table.** All of the PASS cells use K ≥ 1397, beyond the old 1253-coefficient table.
* **Verification.** `verify_ext_cells.py` is the author's `verify_extension.py` logic applied to these cells. It reports PASS for the 4 PASS cells. Its checks:
  * exact transfer of each coefficient ball from the new table;
  * spin caps, real and complex pole gaps, and the need_k prefix condition;
  * the Cauchy majorant and triangle tail, recomputed;
  * the final rational inequality;
  * contiguity from 0.258.

  The ledger is `ext_ledger.json`.
* **New coverage.** c ∈ [0.258, 0.263] for both b, plus [0.263, 0.264] for b = 0 only. The combined top is L = −ln(1 − 0.263) ≈ 0.3052.
* **The tight-variant run.** `direct_sign_tight_ext.py` is the author's tight script with K and the cell endpoints turned into inputs; its `tail(u, K) < 3e-6` assertion is kept. It was tried once and failed for the reason in §3.

## 3. Why it stops at c ≈ 0.263: the h-cap of the rectangle, not the table

* **How the drivers define the rectangle.** Both drivers (cheb and triangle) take the (c, h) rectangle to be h ∈ [0, c_hi/((1001−b)(1002+b))].
  * This contains all physical states only while that cap is the binding constraint.
  * For c > 0.2509 the retained constraint a ≥ 998, that is h ≤ (1−√c)²/998² (README §3, used in full_prefix.py), is smaller.
  * The rectangle therefore includes unphysical points with a < 998.
* **What happens at those points.** There the prefix numerator itself is not positive. These are non-rigorous mpmath probes at 1000 digits (`probe3.py`), corner (c_hi, h_cap):
  * b = 0: 1.5e-5 at c = 0.262; 7.3e-6 at 0.265; −1.5e-5 at 0.266; 2.8e-2 at 0.270; 1.6e1 at 0.274; −4.4e4 at 0.278; 3.6e14 at 0.29.
  * b = 4: 1.30e-5 at 0.262; −1.4e-4 at 0.265.
  * The values do not depend on K: the tight K gives the same numbers (`probe4.py`).
* **Consequence.** Using a larger table, narrower cells or a smaller K cannot get past about c ≈ 0.264 with the drivers as written.
* **What would work (not run; needs your decision).** Replace the cap with h_cap = min(c_hi/((1001−b)(1002+b)), (1−√c_lo)²/998²).
  * This is a one-line change to the h-range in the driver; the C++ is unchanged.
  * It is rigorous because √((n+a)(n+a+s)) − √(n(n+s)) ≥ a.
  * I prepared it as a separate driver, `direct_sign_triangle_cap2.py`, but the permission system blocked it as a modification of the author's code. It was not created or run.
  * Non-rigorous probes with that cap (`probe2.py`) show a smooth, positive numerator over the whole rectangle for b = 0. Corner values (c_lo / c_hi, h at the min-cap):

    | cell | corner values |
    |---|---|
    | [0.29, 0.30] | 1.41e-5 / 1.45e-5 |
    | [0.35, 0.36] | 1.34e-5 / 1.38e-5 |
    | [0.45, 0.46] | 1.23e-5 / 1.25e-5 |
    | [0.55, 0.56] | 1.10e-5 / 1.12e-5 |

  * **Complex pole check.** The driver's smallest h-radius is 2, so it needs 1.5·h_cap·K² < 1. With need_k and the min-cap this holds up to c_hi ≈ 0.51 (h_cap·K² = 0.659 at c = 0.50, 0.673 at 0.52).
  * Above about 0.51 one would need a smaller K (a tight tail allowance) or a driver radius option below 2.

## 4. Cost extrapolation (with the cap change)

* **Scaling.** Per-cell time is roughly proportional to K × terms. Measured: about 40 s at K ≈ 1400 with 2e5 terms; 150 s at K = 1512 with 7.6e5 terms.
* **Estimate to c = 0.5.** K rises to 2700, terms are perhaps 5e5–1.5e6 (the majorant log10 M grows), so roughly 3–10 min per cell and channel. With about 15–25 cells of width 0.01–0.02 for each b, that is about 2–8 CPU-hours, or 1–3 h wall on 3 cores.
* **Above 0.51.** This additionally needs the pole or radius issue above resolved.
* The coefficient table needs no further work for c ≤ 0.56.
