# FULLPOS_STRIPSCAN: uncertified positivity-margin scan of the low-twist bridge strip (Δ_g = 1000)

*26 Nov 2026. Research agent (FULLPOS). **Uncertified [DV]**, kept separate from FULLPOS_CERT.md. Code: `ads4/code/fullpos/fp_stripscan.py`.
Outputs: `cache/fp_stripscan.txt`, `cache/fp_stripscan_edge.txt`.*

## Setup
* **States.** b ∈ {0, 4}; twist a = Δ − J − 2 ∈ {998, 1500, 3000, 10⁴, 3·10⁴, 10⁵}.
* **Choice of J.** J is chosen (even for b = 0) so that L = −log(1−r), r = C/T, is closest to each target in
  {0.04, 0.1, 0.3, 1, 3, 10, 20, 28}. Here 1 − r = (2aJ + a(a+1))/((J+a)(J+a+1)), so J ≈ 2a·e^L for large L.
* **H_ang.** The exact finite ladder sum, evaluated in Arb from the certified K = 1000 moments; H^{(2)} is its k = 2 part.
* **H_hyb and N_s.** T0's `hyb_fast` (float).
* **Normalised margin.** H_ang/N_s (and H^{(2)}/N_s), to be compared with 0.1999.

## Where it stops
The moment code handles n = J − b ≤ 1000, so the largest reachable L is:

| a | 998 | 1500 | 3000 | 10⁴ | 3·10⁴ | 10⁵ |
|---|---|---|---|---|---|---|
| L at n = 1000 | 0.289 | 0.175 | 0.065 | 0.0084 | 0.0011 | 0.0001 |

Every target L ≥ 0.3 needs n ≥ 1030 (a = 998) up to n ~ 10¹⁶ (L = 28), so all of them are skipped. The scan covers only the
low-L end of the strip. It says **nothing** about L ∈ [0.3, 28], the actual bridge.

## Results (every state reachable)
| b | a | L target | J | L | H_ang | H^(2) | H_hyb | N_s | H_ang/N_s | H^(2)/N_s |
|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 998 | 0.04 | 246 | 0.0400 | 3.249327e-05 | 3.109946e-05 | 3.249124e-05 | 3.249124e-05 | 1.00006 | 0.95716 |
| 0 | 998 | 0.1 | 446 | 0.1004 | 1.994082e-05 | 1.886976e-05 | 1.993877e-05 | 1.993877e-05 | 1.00010 | 0.94639 |
| 0 | 1500 | 0.04 | 370 | 0.0400 | 3.171518e-05 | 3.109819e-05 | 3.171316e-05 | 3.171316e-05 | 1.00006 | 0.98061 |
| 0 | 1500 | 0.1 | 670 | 0.1003 | 1.935249e-05 | 1.887814e-05 | 1.935045e-05 | 1.935045e-05 | 1.00011 | 0.97559 |
| 0 | 3000 | 0.04 | 740 | 0.0400 | 3.127892e-05 | 3.112464e-05 | 3.127692e-05 | 3.127692e-05 | 1.00006 | 0.99513 |
| 4 | 998 | 0.04 | 246 | 0.0400 | 2.970853e-05 | 3.110218e-05 | 2.970749e-05 | 3.249513e-05 | 0.91425 | 0.95713 |
| 4 | 998 | 0.1 | 445 | 0.1001 | 1.781860e-05 | 1.889090e-05 | 1.781707e-05 | 1.996175e-05 | 0.89264 | 0.94636 |
| 4 | 1500 | 0.04 | 370 | 0.0400 | 3.048248e-05 | 3.109945e-05 | 3.048089e-05 | 3.171487e-05 | 0.96114 | 0.98060 |
| 4 | 1500 | 0.1 | 669 | 0.1001 | 1.841758e-05 | 1.889229e-05 | 1.841578e-05 | 1.936523e-05 | 0.95106 | 0.97558 |
| 4 | 3000 | 0.04 | 740 | 0.0400 | 3.097067e-05 | 3.112495e-05 | 3.096878e-05 | 3.127735e-05 | 0.99019 | 0.99513 |

**Edge of reach (n = 1000):**

| b | a | L target | J | L | H_ang | H^(2) | H_hyb | N_s | H_ang/N_s | H^(2)/N_s |
|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 998 | edge of reach | 1000 | 0.2885 | 1.522111e-05 | 1.459897e-05 | 1.521906e-05 | 1.521906e-05 | 1.00013 | 0.95926 |
| 0 | 1500 | edge of reach | 1000 | 0.1745 | 1.656628e-05 | 1.619352e-05 | 1.656424e-05 | 1.656424e-05 | 1.00012 | 0.97762 |
| 0 | 3000 | edge of reach | 1000 | 0.0646 | 2.295555e-05 | 2.281874e-05 | 2.295353e-05 | 2.295353e-05 | 1.00009 | 0.99413 |
| 0 | 10000 | edge of reach | 1000 | 0.0083 | 1.844692e-04 | 1.844517e-04 | 1.844669e-04 | 1.844669e-04 | 1.00001 | 0.99992 |
| 4 | 998 | edge of reach | 1004 | 0.2899 | 1.396510e-05 | 1.458520e-05 | 1.396333e-05 | 1.520355e-05 | 0.91854 | 0.95933 |
| 4 | 1500 | edge of reach | 1004 | 0.1754 | 1.580263e-05 | 1.617439e-05 | 1.580077e-05 | 1.654430e-05 | 0.95517 | 0.97764 |
| 4 | 3000 | edge of reach | 1004 | 0.0650 | 2.260851e-05 | 2.274508e-05 | 2.260657e-05 | 2.287972e-05 | 0.98815 | 0.99412 |
| 4 | 10000 | edge of reach | 1004 | 0.0084 | 1.825949e-04 | 1.826123e-04 | 1.825928e-04 | 1.826276e-04 | 0.99982 | 0.99992 |

## Summary
* **All 18 scanned states are comfortably positive.**
  * min H_ang/N_s = **0.8926** (b = 4, a = 998, J = 445, L = 0.10);
  * min H^{(2)}/N_s = **0.9464** (same state);
  * everything is far above 0.1999.
  * For b = 0, H_ang/N_s ≈ 1.0001 everywhere (the k = 3 part has the same sign there).
* **Trend.** At a = 998, b = 4, H_ang/N_s goes 0.914 → 0.893 → 0.919 at L = 0.04 → 0.10 → 0.29. There is no sign of an approaching
  zero, but the range is short.
* **Agreement with the hybrid.** H_ang agrees with T0's hybrid to ≤ 1.3·10⁻⁴ relative on every scanned state, so the margin is
  essentially T0's.
* **Not done.** L ≥ 0.3 at a = 998, and L ≥ 0.1 at larger a, needs n > 1000. It can be reached in two ways:
  * more certified moments (cost ~Õ(n²) at ~n/2 digits);
  * the per-frequency kernel ∫g·𝒦_fin dν (FULLPOS_CERT §4.4). That costs O(n·#nodes) per state at ~0.75ν digits, feasible to n ~ 5·10³
    (L ≈ 1 at a = 998), not to L = 28.
* **CONTOUR note (ANGULAR_ROUTE §43).** Nothing in FULLPOS_CERT assumes an L ≥ 28 (or any L-threshold) tail; the tail stage was shown to
  be unusable at any accessible K. No flag is needed.
