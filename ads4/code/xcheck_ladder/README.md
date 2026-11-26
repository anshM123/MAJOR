# xcheck_ladder: exact ladder formula (Theorem FP) vs the direct rows

Referee-requested consistency check, 26 Nov 2026. **Question:** does the exact finite-mass action of Ω^EE
(`prl/sec_fullpos.tex`, Eqs. fp-target, fp-H) reproduce the values of the same functional computed earlier by the
direct light-ray/K₂/RQMC/Gauss–Kronrod rows (`prl/supplement.tex`, Tables tab:computed, tab:Dstates, tab:gap2000), within
their stated error budgets? Nothing outside this directory was modified. The scripts import `code/largespin/ls_core.py` read-only,
for a cross-check of N_s only, and read the cached certified moments `code/fullpos/cache/moments_1000_K1000.json`
read-only.

**Note on git.** I made no commits. Commits `b767952` and `b1ee978` ("Cross-check (in progress) …"), which picked up
in-progress files from this directory, were made by some other process.

## 1. Normalisation map

* **Exact side.** Eq. (fp-target) gives Ω_O = 𝓑(Δ,J) Σ_b |a_{O,b,0}|² H_ang,b(Δ,J) with 𝓑 > 0. The matrix M_Ω(O) is
  diagonal in the structure basis: H_ang,0 sits on ++ and on parity-odd even J (OO), and H_ang,4 on +− and on
  parity-odd odd J. Each entry carries the positive factor 𝓑·(map a_{b,0} ↔ OPE coefficient).
* **Direct side.** Every direct value is λ_min/N_s (UT_v2 §2.5, supplement l. 273). It is the lowest eigenvalue of
  N^{-1/2} M_Ω N^{-1/2} with N = diag(N_{s,h}) (`spinfunc/ee1000_common.lam2`), where
  `N_{s,h} = |B2 + c2| + |B3 + c3|`, the flat NF02 parts plus the leading-Regge curvature (`ee1000_common.hyb`,
  `ls_core.hyb_fast`). This N_s is the "flat normalisation"; the hybrid H_hyb/N_s ≤ 1 by construction.
* **Link.** Eq. (fp-lift) defines N^phys_{s,O} by the *same* Eq. (fp-target), with H_ang,b → N_{s,b}. The positive
  factor 𝓑|a_{O,b,0}|² is therefore common to M_Ω and to the normalisation, and it cancels. Eq. (fp-LR), checked by
  `fullpos/fp_norm.py` and `fp_compare.py`, fixes the relative normalisation of H_ang,b and T0's per-state action
  (H_hyb) to exactly 1. Hence

  **λ_min/N_s (exact) = min over the block's structures of H_ang,b(Δ,J) / N_{s,b}(m)**, with:
  * even J: b = 0 (++) and b = 4 (+−);
  * odd J (parity-odd 1×1): b = 4;
  * OO (parity-odd even J): b = 0, normalised by N_s(++). This matches the supplement's "normalised ++ diagonal".

  Kinematics: T = (Δ−1)(Δ−2), C = J(J+1), m² = T − C = (Δ−J−2)(Δ+J−1), m/Λ as printed. The twist edge is Δ = J + Δ_g
  exactly. At Δ_g = 2000 and 5000 the "tw₁₀₀₀ image" is taken at the same m/Λ, exactly √((998+2J+1)/999); the table
  prints it to 4 decimals.
* **Checks of the map.** At J ≤ 13 the exact and direct values agree to ≤ 1.4·10⁻³ N_s (2·10⁻⁴ at 2000). The
  *beyond-hybrid corrections* also agree in sign and size, e.g. Δ_g = 1000: J=9 −0.00084 vs −0.00086,
  J=10 +0.00168 vs +0.00154, J=11 −0.00070 vs −0.00068; Δ_g = 2000: J=9 −0.00036 vs −0.00040. A wrong factor (2, 𝓑,
  |a|²) or a wrong N_s would be visible at once. My own mpmath N_s (`xl_norm.py`) equals `ls_core.hyb_fast` to
  3·10⁻¹¹ relative on every Δ_g = 1000 state.
* **Residual ambiguities.** The parity-odd coupling normalisation, graded [I] in the supplement, cancels in the
  per-block ratio. For the 2×2 even-J blocks the direct rows include small measured off-diagonals; the exact formula
  has none (Eq. fp-target). This affects λ_min at most at second order. I found no unresolved normalisation ambiguity.

## 2. How the exact value is computed (independent code)

* `xl_common.py` re-implements Eqs. (fp-U, fp-zy, fp-moments, fp-H) and the single-sequence form (fp-single) from the
  TeX, in Arb (python-flint):
  * the CDH recurrences are exact (fmpq_poly);
  * Q_k is exact rational from `nf02_psi.json`;
  * fp-single agrees with fp-H to 10⁻³⁵⁰.
* **Moments, three independent sources at Δ_g = 1000** (`xl_moments.py`):
  1. **cert:** the certified closed-form Arb moments, read only.
  2. **quad:** moments of a discrete measure Σ wᵢ g(νᵢ) δ(ν−νᵢ), using Gauss–Legendre (24 or 32 nodes/panel) on the
     published repaired continuum weight × θ(Λ+100, 10), up to Λ+70.
  3. **grid:** the published truncated grid `d_abel_X_1000_trunc.txt`, with χ_re = Re X/(R₋+R₊).

  Power moments are accumulated in Arb at 900 or 1400 bits. Because H is linear in the measure, the quad route equals
  ∫g K_T dν exactly, where K_T is the exact per-ν finite-T ladder kernel. The ladder's 10^k cancellation therefore does
  not amplify quadrature error. The conditioning ∫|g K_T|/|H^(2)| is 3 (J=5), 3·10² (J=18), 1·10³ (J=29) and 1.5·10⁴
  (J=100) (`out/kernel_routeB.txt`, which evaluates K_T node by node).
* **Δ_g = 2000 and 5000.** No continuum weight is published at these gaps (supplement Sec. T0: "At Δ_g ≥ 2000 the
  realised statement needs continuum weights there"). Only the grids `prl/data/gap/d_abel_X_{2000,5000}_trunc.txt.gz`
  exist. They are the same weights the 2000/5000 direct rows used. At 1000 the grid route differs from the continuum by
  ≤ 1.6·10⁻⁵ N_s.
* **N_s.** `xl_norm.py` implements the hybrid parts independently, in mpmath with ₂F₁ conical functions, and is
  cross-checked against `ls_core.hyb_fast`.

**Numerical stability** (`out/precision.txt`, `out/xcheck_1000.csv` columns `exact_*`):

| check | result |
|---|---|
| cert vs quad (GL24 / GL32), all 61 states at 1000 | ≤ 7·10⁻⁹ N_s |
| quad at 900 vs 1400 bits | identical to 10⁻¹³ |
| grid vs cert | ≤ 1.6·10⁻⁵ N_s |
| contribution of ν > Λ | ≤ 1.5·10⁻⁵ N_s |
| contribution of ν > m (route B) | ≤ 6·10⁻⁸ in H |

Ladder cancellation, the largest |term|/|H|:

| state | ratio |
|---|---|
| J=18 | 2·10⁸ |
| J=29 | 2·10¹⁵ |
| J=32 | 10²⁰ |
| J=100 | 2·10⁶⁵ |
| J=300 | 10¹⁷⁵ |

64-bit arithmetic fails from J = 29 on, and 1024 bits are needed at J = 300. Every number reported here is stable from
128 bits (J ≤ 32) or 256 bits (J = 100) up to 2048 bits.

## 3. Results

Full table: `out/results.csv` (116 comparisons), also as Markdown in `out/results_table.md`; the summary is in
`out/summary.txt`. In the columns, budget = value/L (or the stated budget), and stat = the jackknife/Monte Carlo error
where quoted. The beyond-hybrid corrections are d_exact = exact − hybrid and d_direct = direct − hybrid.

**Headline numbers**

| state | exact | direct |
|---|---|---|
| Δ_g=1000, J=5 twist edge (binding) | +0.200868 | +0.20088 ± 0.0026 |
| J=11 | +0.348482 | +0.3485 |
| J=18 | +0.74889 | +0.746 ± 0.41 (stat 0.018) |
| J=29 twist edge | +0.99839 | +0.972 ± 0.784 (MC 0.018) |
| J=29 at Δ_g=2000, its twist edge 1.0144Λ | +0.92111 | +0.9407 |
| J=29 at Δ_g=2000, the 1000-image 1.0286Λ | +0.99966 | +1.0357 |

**Agreement within budgets.** All 116 differences lie well inside the stated budgets:

| Δ_g | max \|diff\|/budget | worst state |
|---|---|---|
| 1000 | 0.084 | J=31, 1.04Λ |
| 2000 | 0.215 | J=27, image mass |
| 5000 | 0.44 | J=5 twist edge; 0.20015 vs printed 0.2002, so only the printed rounding |

No exact value is ≤ 0. The smallest is +0.20015, J = 5 at the 5000 twist edge.

**Agreement in detail at J ≤ 13.**

| Δ_g | max \|diff\| (N_s) |
|---|---|
| 1000 | 1.4·10⁻³ |
| 2000 | 1.8·10⁻⁴ |

At this level the finite-mass (beyond-leading-Regge) corrections of the two methods match in sign and magnitude. This is
the non-trivial part of the test.

**Tension with the quoted statistical errors alone (not with the budgets).** Where a jackknife or Monte Carlo error is
quoted, these states differ by more than 2 stat:

| J (Δ_g = 1000) | diff | in units of stat |
|---|---|---|
| 15, twist edge | −0.0026 | −4.4 |
| 16, twist edge | +0.0027 | +2.2 |
| 30, twist edge | +0.076 | +4.2 |
| 31, twist edge | +0.070 | +2.2 |
| 32, twist edge | +0.045 | +2.5 |

* The exact beyond-hybrid correction is ≤ 0.005 N_s for every J ≤ 32 near the gap.
* The direct rows at J ≥ 17 imply "corrections" of up to ±0.075 N_s. These are therefore almost entirely numerical
  systematics of the direct method: edge term, window, missing nodes. That is consistent with their budgets being large
  (L ≈ 0.4–1.8).
* At Δ_g = 2000 the direct values of J = 20, 21, 23, 27 and 29 all lie above the exact ones, by 0.01–0.065. The largest
  is J = 20 at every mass: +0.054 to +0.065, about 0.1–0.2 of the budget.

**Artefacts.** The exact values settle the artefacts flagged in Secs. S6–S8:

| state | earlier value | exact | verdict |
|---|---|---|---|
| J=18 | unmapped-basis −0.12 | +0.749 | artefact |
| J=29 | chain −3.62 | +0.998 | artefact |
| J=29 | spline +1.39 | +0.998 | artefact |
| J=29 | published chain +0.948 | +0.998 | — |
| J=20 | chain +0.965 | +0.839 | the direct row's +0.835 is right |
| J=31 | +1.20 → +1.02 | +0.990 | — |
| J=32 | +0.904 → +0.997 | +1.000 | — |

**(1100, 100)** (`out/window_1100_100.txt`). Exact H/N_s = +1.000001 (++) and +0.981007 (+−). This is within
2·10⁻⁶ N_s of the hybrid, a total beyond-hybrid correction of about ±2 c₂. The supplement's window numbers
(−2154 ± 43 c₂ vs T_B = 1409; leverages of Table tab:dircv) are per-ν fit residuals of Γ on ν = 11–25, not values of the
action, so they cannot be compared number by number. The exact result does confirm their conclusion: the excess has no
effect on the sign, and the state is positive with the full hybrid margin.

**Minimum over masses** (`out/minscan_1000_*.txt`). The minimum of λ/N_s over masses from the twist edge to 5Λ, for
J ≤ 48 at Δ_g = 1000, against Table tab:Dstates: see §4.

## 4. Minimum over masses

(filled from `out/minscan_*.txt`; see the end of this file)

## Files

| file | content |
|---|---|
| `xl_common.py` | the exact formula: CDH polynomials, Q_k, H_ang (fp-H), H₂ single-sequence (fp-single), weights, moment accumulation |
| `xl_moments.py` | Arb power moments of the GL / continuum and grid measures → `out/pmom_*.json` |
| `xl_norm.py` | independent N_s and hybrid (mpmath) |
| `xl_states.py` | transcription of every direct value used, with budgets (value/L), stat errors and notes |
| `xl_eval.py` | evaluates all states and writes `out/xcheck_<Dg>.csv` |
| `xl_report.py` | merges into `out/results.csv` and `out/summary.txt` |
| `xl_extra.py` | precision/cancellation (`prec`), (1100,100) (`win`), per-ν kernel route B (`kernel`), mass scan (`scan`) |
| `out/*.log` | run logs |

Reproduce, in about 10 min, excluding the mass scan:

```
python3 xl_moments.py cont 1000 24 120 1400; python3 xl_moments.py cont 1000 32 120 1400
python3 xl_moments.py cont 1000 24 120 900;  python3 xl_moments.py grid 1000 0 120 1400
python3 xl_moments.py grid 2000 0 120 1400;  python3 xl_moments.py grid 5000 0 120 1400
python3 xl_eval.py 1000 2000 5000; python3 xl_report.py
python3 xl_extra.py prec; python3 xl_extra.py win; python3 xl_extra.py kernel
for j in 0 1 2; do python3 xl_extra.py scan 48 $j 3; done
```
