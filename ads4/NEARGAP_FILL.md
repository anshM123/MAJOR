# NEARGAP_FILL: error control at Δ_g = 1000 for the near-gap spins still model-level (referee report 5, M1)

*24 Nov 2026, near-gap fill agent. Unreviewed, not committed by this agent. Other agents' commits (9d1797b, 9550bf3) swept in intermediate
snapshots of `code/ngfill/` and `ext/jobs_batch5.txt`; the working tree holds the final versions. No existing file was edited. The one new
file outside `code/ngfill/` is `ext/jobs_batch5.txt`. `GAP_SEAM.md` and `code/gapseam/` were not touched.*
*One write escaped the guard: Python's bytecode cache `code/ngfill/__pycache__/`, from a scratch import of `nf_core` before the guard was
active. It was deleted.*
* *New code:* `ads4/code/ngfill/`.
  * `nf_guard.py`: process-local write guard. Every write outside `ngfill/cache/` raises. It registers itself as `ng_guard`, `hs_guard`,
    `sn_guard`, `ef_guard` and `gn_guard`.
  * `nf_core.py`: read-only set-up. It runs `highspin/hs_core` → `neargap2/ng_val` → `ng_chain` → the w2/q3 chain, and adds `ngfill/cache`
    to the row loader, the direct-point directory and the shadow candidates.
  * `nf_rows.py` = `neargap2/ng_rows.py` with every cache path moved. `nf_points.py` = `highspin/hs_points.py` with the cache moved.
    `nf_shadow.py`: the odd ⟨TTO⟩ shadow without a row.
  * `nf_state.py`: the HIGHSPIN §4 own-row protocol (`neargap2/ng_bind.run_state`, unchanged) for both parities. `nf_nodes.py`: direct
    vs K₂ node tables. `nf_summary.py`.
  * `nf_repr.py`, `nf_reprcut.py`, `nf_reprsum.py`: the J = 23 / 31 representation analysis (§3).
  * `nf_batch5.py`: writes `ext/jobs_batch5.txt`. `nf_b5pilot.py`: the 1-node pilot.
  * `nf_lane.sh`, `nf_pause.sh`: compute lanes, and short analyses with one worker SIGSTOPped.
* *Caches:* `ads4/code/ngfill/cache/`; its `.gitignore` ignores everything but itself.
* *Read-only imports:*
  * `highspin/` (`hs_core`, `hs_points`, `hs_odd` via runpy), `neargap2/` (`ng_val`, `ng_chain`, `ng_bind`), `smallnu/` (`sn_core`,
    `sn_odd`, `sn_window`), `extfold/ef_core` (its merged point sets), `gapnear/gn_batch4` (pure functions);
  * `spinfunc/` (`w1_rows`, `p3_odd_action`, `p2_shadow_odd`, `l9_points`, `l5_shadow_matrix`, `w0_ee_common`, `hc_points`, `hc_nodes`,
    `ee1000_common` → `w2_chain`, `q3_analysis`, `p4_odd_analysis`); `ext/point.py` (pilot).
* *Compute:* ≈ 6.3 CPU-h, with at most two computing processes (two lanes 11:55–14:58 UTC; analyses ≈ 0.3 h, run with one worker SIGSTOPped or after the lanes).
  `ulimit -v 3000000`, `OMP_NUM_THREADS=1`.

Grades as in POS_ERRORS / NEARGAP2 / HIGHSPIN: **[NE±]** means an explicit error budget, every model input with a stated safety factor
(×3 on model terms, 3σ on Monte Carlo), and the conclusion surviving it. **[NE + M]**, **[M]**, **[I]**. `L = value / budget`.
"Sign-free" and "λ-level" are as in NEARGAP2 §1.2. Γ is in `c₂ = m² ×` relative correction. Λ = 998.5.

---

## 0. Result

**Result table** (Δ_g = 1000; budgets are sign-free unless marked; a mass marked "twist" is the twist edge Δ = 1000 + J):

| J | state (m/Λ) | value λ/N_s (old → new) | budget | L | grade old → new |
|---|---|---|---|---|---|
| **9** | twist 1.0090 | +0.6755 (published chain) → **+0.6757 ± 0.0001** (own direct row) | 0.018 | **38** | [NE + M] → **[NE±]** |
| 9 | 1.02 / 1.03 / 1.05 | +0.4859 / +0.4543 / +0.5294 | 0.017 / 0.020 / 0.032 | **29 / 23 / 17** | [NE + M] → **[NE±]** |
| **10** | twist 1.0100 | +0.5852 → **+0.5892 ± 0.0004** | 0.106 | **5.6** | [NE + M] → **[NE±]** |
| 10 | 1.02 / 1.03 / 1.05 | +0.8058 / +0.9186 / +0.8972 | 0.125 / 0.149 / 0.161 | **6.5 / 6.2 / 5.6** | [NE + M] → **[NE±]** |
| **11** | twist 1.0110 | +0.3489 → **+0.3485 ± 0.0002** | 0.022 | **16** | [NE + M] → **[NE±]** |
| 11 | 1.02 / 1.03 / 1.05 | +0.4099 / +0.5453 / +0.9997 | 0.022 / 0.027 / 0.049 | **19 / 20 / 20** | [NE + M] → **[NE±]** |
| **13** | twist 1.0129 | +0.4094 → **+0.4076 ± 0.0003** | 0.060 | **6.8** | [NE + M] → **[NE±]** |
| 13 | 1.02 / 1.03 / 1.05 | +0.6034 / +0.9994 / +1.0002 | 0.066 / 0.093 / 0.081 | **9.1 / 10.8 / 12.3** | [NE + M] → **[NE±]** |
| **14** | twist 1.0139 | +0.9965 → **+1.0009 ± 0.0006** | 0.483 | **2.07** | [NE + M] → **[NE±]** |
| 14 | 1.02 / 1.03 / 1.05 | +1.0011 / +0.7264 / +0.7946 | 0.460 / 0.396 / 0.523 | **2.2 / 1.8 / 1.5** (λ-level 2.6 / 3.0 / 3.9) | [NE + M] → **[NE±]** |
| **15** | twist 1.0149 | +0.6900 (own K₂ row; published +0.6751) → **+0.6943 ± 0.0005** | 0.397 | **1.75** | [NE + M] → **[NE±]** |
| 15 | 1.02 / 1.03 / 1.05 | +0.9998 / +1.0005 / +0.9493 | 0.403 / 0.232 / 0.116 | **2.5 / 4.3 / 8.2** | [NE + M] → **[NE±]** |
| **23** | twist 1.0228 | +0.923 ± 0.016 (λ_c +0.901), unchanged | 0.762 (×3 window; was 1.056) | 0.85 → **1.18** (signed 1.27) | [NE + M] → **[NE±] (thin)**, §3 |
| 23 | 1.03 / 1.05 / 1.075 / 1.10 | unchanged | | 1.16 / 1.69 / 2.09 / 2.33 → **1.87 / 3.56 / 5.22 / 6.22** | [NE±] (strengthened) |
| **31** | twist 1.0306 | +0.862 ± 0.023 (λ_c +0.881), unchanged | 1.158 | 0.70 → **0.76** (signed 0.94) | **[NE + M] unchanged** |
| 31 | 1.04 | +0.869 | 0.811 | 1.00 → **1.09** (signed 1.45) | [NE + M] (borderline) → **[NE±]** |
| 31 | 1.05 / 1.075 / 1.10 | unchanged | | 1.22 / 1.55 / 1.90 → 1.35 / 1.75 / 2.18 | [NE±] (strengthened) |
| 27 (control) | twist 1.0267 | +0.976 | 0.551 | 1.64 → 1.81 | [NE±] (unchanged; nothing got worse) |
| 17–20, 22, 24–26, 28, 30, 32 | twist and low masses | not computed here | — | — | unchanged ([NE + M]; J = 26, 30 [M]); **batch 5 prepared**: 666 n1 jobs, ≈ 29 core-h, ≈ 1.8 h on 16 workers (§4) |


**Inventory before this note** (Δ_g = 1000, near-gap J ≤ 32, twist edge and masses to ≈ 1.05–1.10Λ; from HIGHSPIN, SMALLNU, EXT_FOLD
+ addendum, NEARGAP2, POS_ERRORS and the UNIFIED_THEOREM boxes):
* **[NE±]:** J = 0, 2, 4, 5, 6, 7, 8 (λ-level), 12, 16, 21, 27, 29 (signed window); J = 23 at ≥ 1.03Λ; J = 31 at ≥ 1.05Λ.
* **Not [NE±]:**
  * J = 9, 10, 11, 13, 14, 15 [NE + M] (J = 15 on its own K₂ row, the rest on the published 1/m chain);
  * J = 17, 18, 19, 20, 22, 24, 25, 28, 32 [NE + M];
  * J = 26, 30 [M] (no rows);
  * J = 23's twist edge [NE + M]; J = 31 twist edge – 1.04Λ [NE + M].
  * (J = 1, 3 have no state: odd J has only the parity-odd `+−` structure, which starts at J = 5.)

**Negative findings and cautions (stated up front).**
1. **No state came out negative.** The lowest new value is J = 11 at its twist edge, `+0.3485` (lower bound `+0.327`). Next are J = 13
   (`+0.408`) and J = 9 at 1.03Λ (`+0.454`), all well above the binding J = 5 state (`+0.20`).
2. **J = 31's twist edge does not close** (`L = 0.76` ×3, `0.94` signed), with any defensible representation term. Its half-sample probe
   is a real (5–6σ) under-resolution, not noise (§3).
3. **The J = 23 closure is thin** (`L = 1.18`). It uses a *corrected evaluation* of the representation term: on λ_c rather than on λ,
   because the probes' extrapolation below the first bulk node lies in a window that the budget bounds separately. That is a protocol
   correction, argued in §3; it is not a new tolerance.
4. **The "jackknife-weighted" (WLS) idea does not work for the reason expected.** The direct-row node errors are nearly homoscedastic
   (36 c₂ throughout), so WLS gives only `L = 1.01` at J = 23. It is reported, not used.
5. **J = 14 and 15 are the thinnest of the J ≤ 16 fill** (`L ≥ 1.52`, `≥ 1.75`).
   * J = 14 is limited by the ν < 50 convention (0.23) and the edge beyond x = 0.925 (0.25; the x = 0.95 node failed). At 1.05Λ it is
     limited by the 1600 K₂ partner row's truncation (0.34).
   * J = 15 is limited by the edge (0.30, x = 0.95 failed).
   * With the measured small-ν nodes standing for ν ∈ [25, 50) (HIGHSPIN's `NU_LOW = 25` variant) they become `L ≥ 1.97` and `≥ 2.12`.
6. **The certified-floor base route does not close for J ≥ 13** (`(N_f/5 − Lemma C)/N_s` = +0.117 at J = 13, < 0 at J = 14, 15): the
   Lemma-C curvature bound is loose there. The base is then the hybrid (exact flat + leading-Regge curvature by quadrature), exactly as for
   J = 12, 16 in HIGHSPIN. For J = 9, 10, 11 the certified route alone gives `λ/N_s ≥ +0.13 … +0.19` before the budget.
7. **Held rows (odd J = 11, 13, 15).** These J have no K₂ row *above* the own row (the published rows are 400/800 and 336/400), so a
   state with m > m_row is held at the own row. An estimate of that hold (§2.3) is `≤ 1.1·10⁻³ N_s`, against budgets ≥ 0.02. It is
   reported, not added; adding it changes no `L` by more than 1 %.

**Net for referee M1.** At Δ_g = 1000 every near-gap spin J ≤ 16 is now [NE±] at its twist edge and low masses: J = 0, 2, 4–16, with
J = 8 at λ level as before. Above that, J = 21, 27, 29 and now J = 23 are [NE±] at every computed mass (J = 23 thin, J = 29 signed), and
J = 31 is [NE±] from 1.04Λ. Of the referee's model-level near-gap list (J = 9–11, 13–15, 17–20, 22, 24, 25, plus the twist edges of 23
and 31), **J = 9, 10, 11, 13, 14, 15, the J = 23 twist edge and J = 31 at 1.04Λ move to [NE±]**. Still open:
* J = 17–20, 22, 24, 25 [NE + M];
* J = 26, 28, 30, 32 [NE + M] / [M];
* the J = 31 twist edge.

All of them are covered by batch 5.

---

## 1. Protocol (J ≤ 16)

The protocol is HIGHSPIN §4's, the one that moved J = 12 and 16 to [NE±], applied to both parities as EXT_FOLD did for J = 7.
* **Rows.** For each J ∈ {9, 10, 11, 13, 14, 15}, a full **direct n1 row** sits on the own-mass row Δ = 1000 + J (the twist edge):
  * x = 0.025 … 0.95 (38 nodes), `hc_points` n1 adaptive, `HC_M = HC_MREF = 11`, RQMC 2¹¹·2 points × 8 replicas;
  * state and Δ_ref = 2·10⁴ on common random numbers.

  Only x = 0.95 failed (replica error > 20 %), at J = 14 and 15. That node is dropped and covered by the edge term.
* **K₂ own rows** (`nf_rows.py` = NEARGAP2's `ng_rows`, L1_DEG = 6). They are coarse (x = 0.05, 0.1, 0.2, …, 0.9, 10 nodes) for
  J = 9, 10, 11, 13, 14; J = 15 uses NEARGAP2's existing 36-node row. They give the K₂-only own-row value, the K₂ λ-level budget, the
  direct-vs-K₂ calibration and the ⟨TTO⟩ shadows that the direct records need.
  * Coarse rather than dense: the direct row overrides every K₂ node anyway, and at ≈ 10–27 min per 10-node row a dense K₂ row would have
    cost as much as the direct row it is checked against.
  * For odd J = 9, 11 the shadow was computed first on its own (`nf_shadow.py`, identical function), so that the direct rows could run before the K₂ rows.
* **Reference zero.** x₀ = 2.76 / 2.27 / 2.20 / 1.84 / 1.64 / 1.58 for J = 9 / 10 / 11 / 13 / 14 / 15 on Δ = 1000 + J, so the
  Δ_ref = 2·10⁴ zero is outside every row. No node is affected, and no other Δ_ref is needed for J ≤ 16 at 1000. GAP_NEARGAP's problem
  at 2000 comes from the rows at 2000 being twice as heavy (x₀ halves).
* **State and budget** (`nf_state.py`): `neargap2/ng_bind.run_state`, unchanged (the EE_NEARGAP w2 chain, hybrid states, NEARGAP2 §1.2
  budget), with the direct row injected through `ng_val`'s loaders. The terms, all sign-free, ×3 on every model term:
  * `E_trunc`: 0 on the own row (direct). At m > m_row it comes from the next K₂ row up (1600 for J = 9, 10, 14).
  * `E_interp`.
  * `E_low`: |Γ| ≤ 10⁴ c₂ on ν < 50, with the min-norm gauge factor.
  * `E_edge`: beyond the last valid node, at max(10⁴, 10 max|Γ| of the last 4 nodes).
  * `E_mc` = 3σ: the jackknife over the 8 replicas, every injected node left out together.

  The masses are the twist edge, 1.02, 1.03 and 1.05Λ: the grid used for J = 12, 15, 16 in NEARGAP2 / HIGHSPIN.
* **Variant `NU_LOW = 25`:** the direct nodes at ν ≈ 25 and 50 (x = 0.025, 0.05) stand for the window ν ∈ [25, 50), with 10⁴ c₂ only
  below 25. The measured |Γ| there is ≤ 60 c₂ (σ ≤ 70; `nf_nodes_J*.txt`), ×100 below the convention.

## 2. Results (J ≤ 16)

### 2.1 State table (`cache/nf_summary.txt`, from `nf_state_J*.json`)

| J | m/Λ | published model | K₂ own row (L sign-free / λ-level) | **direct own row** λ ± jk | budget: low / edge / trunc / interp / mc | lower bound | **L** (NU_LOW = 25) | certified-floor route |
|---|---|---|---|---|---|---|---|---|
| 9 | 1.0090 | +0.6755 | +0.6755 (2.13 / 10.69) | **+0.6757** ± 0.0001 | 0.013 / 0.004 / 0.000 / 0.0001 / 0.0004 | +0.658 | **37.96** (83.49) | +0.188 |
| 9 | 1.0200 | +0.4858 | +0.4858 (2.27 / 12.11) | **+0.4859** ± 0.0001 | 0.010 / 0.000 / 0.006 / 0.0001 / 0.0002 | +0.469 | **29.41** (50.59) | +0.191 |
| 9 | 1.0300 | +0.4542 | +0.4541 (2.40 / 12.88) | **+0.4543** ± 0.0001 | 0.009 / 0.000 / 0.010 / 0.0001 / 0.0002 | +0.435 | **23.27** (34.52) | +0.192 |
| 9 | 1.0500 | +0.5292 | +0.5292 (2.59 / 13.13) | **+0.5294** ± 0.0001 | 0.010 / 0.000 / 0.022 / 0.0001 / 0.0002 | +0.497 | **16.60** (21.50) | +0.191 |
| 10 | 1.0100 | +0.5852 | +0.5885 (0.51 / 2.92) | **+0.5892** ± 0.0004 | 0.077 / 0.028 / 0.000 / 0.0005 / 0.0012 | +0.483 | **5.55** (11.67) | +0.146 |
| 10 | 1.0200 | +0.8026 | +0.8050 (0.68 / 3.95) | **+0.8058** ± 0.0004 | 0.079 / 0.015 / 0.029 / 0.0005 / 0.0013 | +0.681 | **6.46** (11.91) | +0.149 |
| 10 | 1.0300 | +0.9174 | +0.9176 (0.76 / 4.31) | **+0.9186** ± 0.0003 | 0.079 / 0.009 / 0.058 / 0.0015 / 0.0010 | +0.770 | **6.16** (10.01) | +0.152 |
| 10 | 1.0500 | +0.8961 | +0.8962 (0.94 / 4.74) | **+0.8972** ± 0.0003 | 0.067 / 0.002 / 0.090 / 0.0013 / 0.0009 | +0.736 | **5.58** (7.97) | +0.158 |
| 11 | 1.0110 | +0.3489 | +0.3484 (0.60 / 4.59) | **+0.3485** ± 0.0002 | 0.015 / 0.006 / 0.000 / 0.0001 / 0.0007 | +0.327 | **15.99** (32.94) | +0.176 |
| 11 | 1.0200 | +0.4099 | +0.4097 (0.61 / 4.81) | **+0.4099** ± 0.0003 | 0.019 / 0.002 / 0.000 / 0.0002 / 0.0008 | +0.388 | **18.98** (52.22) | +0.171 |
| 11 | 1.0300 | +0.5450 | +0.5449 (0.62 / 4.85) | **+0.5453** ± 0.0003 | 0.026 / 0.000 / 0.000 / 0.0002 / 0.0010 | +0.518 | **19.93** (63.29) | +0.161 |
| 11 | 1.0500 | +0.9984 | +0.9988 (0.67 / 4.88) | **+0.9997** ± 0.0006 | 0.047 / 0.000 / 0.000 / 0.0004 / 0.0018 | +0.950 | **20.21** (65.40) | +0.130 |
| 13 | 1.0129 | +0.4094 | +0.4066 (0.21 / 1.92) | **+0.4076** ± 0.0003 | 0.035 / 0.024 / 0.000 / 0.0009 / 0.0008 | +0.348 | **6.81** (11.68) | +0.117 |
| 13 | 1.0200 | +0.6045 | +0.6018 (0.21 / 2.20) | **+0.6034** ± 0.0004 | 0.053 / 0.011 / 0.000 / 0.0014 / 0.0011 | +0.537 | **9.14** (21.57) | +0.074 |
| 13 | 1.0300 | +0.9990 | +0.9969 (0.22 / 2.55) | **+0.9994** ± 0.0006 | 0.086 / 0.002 / 0.000 / 0.0024 / 0.0017 | +0.907 | **10.79** (33.29) | -0.004 |
| 13 | 1.0500 | +0.9978 | +0.9977 (0.27 / 2.86) | **+1.0002** ± 0.0005 | 0.077 / 0.000 / 0.000 / 0.0023 / 0.0016 | +0.919 | **12.29** (39.51) | +0.022 |
| 14 | 1.0139 | +0.9965 | +0.9955 (0.21 / 0.63) | **+1.0009** ± 0.0006 | 0.227 / 0.250 / 0.000 / 0.0037 / 0.0019 | +0.518 | **2.07** (3.14) | -0.136 |
| 14 | 1.0200 | +0.9966 | +0.9959 (0.23 / 0.76) | **+1.0011** ± 0.0006 | 0.215 / 0.154 / 0.085 / 0.0035 / 0.0017 | +0.541 | **2.18** (3.29) | -0.119 |
| 14 | 1.0300 | +0.7237 | +0.7234 (0.25 / 0.99) | **+0.7264** ± 0.0004 | 0.175 / 0.058 / 0.161 / 0.0003 / 0.0012 | +0.331 | **1.84** (2.70) | -0.099 |
| 14 | 1.0500 | +0.7919 | +0.7921 (0.29 / 1.36) | **+0.7946** ± 0.0004 | 0.167 / 0.015 / 0.340 / 0.0002 / 0.0012 | +0.271 | **1.52** (1.97) | -0.064 |
| 15 | 1.0149 | +0.6751 | +0.6899 (0.14 / 1.30) | **+0.6943** ± 0.0005 | 0.096 / 0.298 / 0.000 / 0.0014 / 0.0016 | +0.297 | **1.75** (2.12) | -0.140 |
| 15 | 1.0200 | +0.9742 | +0.9937 (0.15 / 1.41) | **+0.9998** ± 0.0008 | 0.138 / 0.260 / 0.000 / 0.0023 / 0.0023 | +0.597 | **2.48** (3.30) | -0.285 |
| 15 | 1.0300 | +0.9781 | +0.9944 (0.16 / 1.60) | **+1.0005** ± 0.0007 | 0.131 / 0.095 / 0.000 / 0.0024 / 0.0022 | +0.769 | **4.32** (7.34) | -0.255 |
| 15 | 1.0500 | +0.9322 | +0.9439 (0.20 / 1.99) | **+0.9493** ± 0.0006 | 0.107 / 0.005 / 0.000 / 0.0021 / 0.0018 | +0.833 | **8.16** (24.61) | -0.163 |

### 2.2 Checks

* **Direct vs K₂, node by node** (`cache/nf_nodes_J*.txt`; the common nodes of the coarse K₂ rows, J = 15 its 36-node row):

  | J | 9 | 10 | 11 | 13 | 14 | 15 |
  |---|---|---|---|---|---|---|
  | median / max \|direct − K₂\| (c₂) | 28 / 84 | 26 / 103 | 32 / 155 | 50 / 154 | 29 / 123 | 78 / 359 |
  | median \|Δ\|/σ | 1.5 | 1.2 | 1.4 | 2.4 | 0.9 | 2.0 |

  The largest single deviations (4–12σ) are at x = 0.05–0.4. There K₂ overstates |Γ| by 50–150 c₂: for example J = 15, x = 0.25, C3H
  `−182` (K₂) vs `−50 ± 11` (direct); J = 11, x = 0.2, C2 `−67` vs `−9`. The pattern is the one HIGHSPIN and EXT_FOLD found: K₂ has
  real errors of `10²` c₂ that oscillate with x and largely cancel in λ (net shifts below, all favourable).
* **K₂ calibration extended to six more states** (twist edge; actual K₂ error vs its last Laplace step):

  | J | 9 | 10 | 11 | 13 | 14 | 15 |
  |---|---|---|---|---|---|---|
  | λ_direct − λ_K₂ | +0.0002 | +0.0009 | +0.0001 | +0.0010 | +0.0067 | +0.0027 |
  | λ_K₂ − λ_K₁ | +0.012 | +0.030 | +0.014 | +0.029 | +0.102 | +0.096 |
  | ratio | 2 % | 3 % | 0.4 % | 3 % | 7 % | 3 % |

  Together with EXT_FOLD's seven states, the actual K₂ error is **0.4–8 % of the last step over 13 states, always with direct > K₂
  (favourable)**. (The coarse K₂ rows carry some interpolation in this comparison; J = 15's row is dense.)
* **K₂-only budgets on the new own rows.** NEARGAP2's protocol alone (no direct points) gives sign-free `L` = 2.1–2.6 (J = 9),
  0.5–0.9 (10), 0.6–0.7 (11), 0.2–0.3 (13, 14), 0.14–0.2 (15); λ-level 10.7–13 (9), 2.9–4.7 (10), 4.6–4.9 (11), 1.9–2.9 (13),
  0.6–1.4 (14), 1.3–2.0 (15). Only J = 9 would have closed on K₂ alone (sign-free). The direct rows are what close J = 10–15.
* **Published chain vs own direct row:** `|Δλ| ≤ 0.004` at J = 9–14. At J = 15 it is `+0.019 … +0.026`, because the published chain holds
  the 400 row; the shift is favourable. NEARGAP2 already saw `+0.015` with its own K₂ row. **H_M is confirmed at J = 9–14.**
* **Small ν, measured:** direct |Γ| at ν ≈ 25 and 50 is ≤ 61 c₂ (σ ≤ 70) at every J. `G_low = 10⁴` is ×100 above the measurement, yet
  it is the largest or second-largest budget term at every J ≤ 14 (the others being the edge term and, above the own row, the
  partner row's truncation).

### 2.3 The held-row estimate (odd J = 11, 13, 15)

The 1/m model holds a state above the highest row at that row. For J = 11 and 13 (published rows 400 and 800) and J = 15 (336 and 400),
the own row is the highest, so states at 1.02–1.05Λ are held at m_row = 1.011–1.015Λ. Estimate: the 1/m step from the own row to the
state, as a fraction f of the step from the next-lower row to the own row, times the λ shift that step produced (direct own row minus the
published chain, which held the lower row), ×3.

| J | m/Λ | f | Δλ (own − lower) | E_hold | `L` with E_hold |
|---|---|---|---|---|---|
| 11 | 1.02 / 1.03 / 1.05 | 0.034 / 0.070 / 0.141 | −0.0001 / +0.0003 / +0.0013 | ≤ 0.0005 | 19.0 / 19.9 / 20.0 |
| 13 | 1.02 / 1.03 / 1.05 | 0.026 / 0.062 / 0.132 | −0.0012 / +0.0005 / +0.0025 | ≤ 0.0010 | 9.1 / 10.8 / 12.2 |
| 15 | 1.02 / 1.03 / 1.05 | 0.003 / 0.009 / 0.022 | +0.026 / +0.022 / +0.017 | ≤ 0.0011 | 2.48 / 4.31 / 8.09 |

This is negligible. It is also the H_mono-A behaviour (Γ at fixed (J, x) nearly Δ-independent) that GAP_NEARGAP measured over
1000 → 2000. The cleaner fix would be a direct row at Δ = 1050 + J, ≈ 45 min each (not run).

### 2.4 Verdict (J ≤ 16)

**J = 9, 10, 11, 13, 14, 15 near the gap (twist edge to 1.05Λ) are [NE±]** (were [NE + M]). The sign-free leverages are:

| J | 9 | 10 | 11 | 13 | 14 | 15 |
|---|---|---|---|---|---|---|
| min `L` | 16.6 | 5.6 | 16.0 | 6.8 | 1.52 | 1.75 |
| min `L`, `NU_LOW = 25` | 21.5 | 8.0 | 32.9 | 11.7 | 1.97 | 2.12 |

* The values are on full direct own-mass rows, with no Laplace truncation on the own row.
* The model inputs left are those of every [NE±] near-gap state (J = 5, 12, 16): the direct n1 integral as the crystal integral; the
  hybrid base (for J ≥ 13 the certified-floor route does not close); `G_low = 10⁴` on ν < 50; `G_edge` beyond the last node.
* With J = 0, 2, 4–8, 12 and 16 this makes **every near-gap spin J ≤ 16 [NE±] at Δ_g = 1000**.
* Since Lemma T (GAP_UNIFORM) carries [NE±] at 1000 to every larger gap under the no-growth premise, these spins now enter the all-gap
  statement at the same grade as J = 12, 16.
  * *Caveat:* at 2000, GAP_SEAM (another agent, in progress) reports the seam spins J = 9–13 sign-free and 14–21 at λ level.


---

## 3. J = 23 and J = 31: a noise-aware representation term (`nf_repr.py`, `nf_reprcut.py`, `nf_reprsum.py`)

**Reproduction first.** `nf_repr.py J ols` reruns `highspin/hs_odd.py` and `smallnu/sn_odd.py` unchanged (runpy / import) on EXT_FOLD's
merged point set `all`. It reproduces EXT_FOLD to every printed digit: J = 23 twist edge `+0.923 ± 0.016`, λ-level budget `0.925`
(mc 0.047, quad 0.023, miss 0.303, repr 0.553), `L` ×3 / signed `0.85 / 0.90`; J = 31 `+0.862`, `0.828`, `L` `0.70 / 0.85`; J = 27
`1.64 / 2.02`.

**Three candidate fixes were tested. Only one survives as a principled change.**

1. **Jackknife-weighted (WLS) shape fit** (`nf_repr.py J wls`): weights `1/σᵢ²` from each node's replica jackknife, held fixed in the
   jackknife. *Premise false:* the node errors on these direct rows are nearly homoscedastic (J = 23: median σ 36 c₂; x < 0.3: 35;
   x > 0.6: 38). WLS moves the degree-7 probe only from `+0.184` to `+0.134`, and the twist edge from `L = 0.85` to `1.01` (signed
   `1.07`). That is not a robust closure, and it is not the mechanism. **Rejected as the fix** (kept as a sensitivity check).
2. **Model selection by χ².** With the jackknife node errors, every degree from K = 4 up fits at the noise level. For J = 23, C2 / C3H
   give χ²/dof = 13.7/19, 14.1/19 (K4), 9.6/13, 9.3/13 (K6) and 7.0/10, 6.6/10 (K7). The K6 → K7 step improves χ² by 5.3 for 6 extra
   parameters (p ≈ 0.5). So the data do not support K7. But they do not *exclude* it either. Dropping the probe on that ground alone
   would be a judgement call, so this was not used.
3. **Evaluate the representation term on the quantity the budget bounds (adopted).** Each probe is split by x-window
   (`nf_reprcut.py`, with a jackknife on every difference):
   * On the first bulk node, x = 0.05 (ν ≈ 50), every shape fit continues as a polynomial *extrapolation* down to x = 0.
   * In the final budget that window (ν < 50) is **removed** from the central value λ_c and covered separately by SMALLNU's measured
     envelope (×3) or signed window (`sn_odd`, `lamc = main − X_low`).
   * The protocol's `E_repr = 3 max|λ(probe) − λ(K6)|` is evaluated on λ, so it includes each probe's extrapolation below its first node.
     For J = 23, that piece is **+0.120 of the K7 probe's +0.184** (K6's own continuation is +0.022). It is pure polynomial extrapolation
     of a 24-parameter fit into a window that the budget covers anyway. Counting it in E_repr counts the window twice.
   * Evaluated on λ_c, the same probes, with the same ×3 rule and nothing dropped, give a K7 deviation of `+0.086 ± 0.106` (0.8σ: noise).
     E_repr goes `0.553 → 0.258`.

   This is a consistency correction, not a tolerance change: every probe is kept, ×3 is kept, and no noise is subtracted. It is applied
   to every odd direct-bulk row that has the data (J = 23, 27, 31), and it lowers E_repr in every case (none gets worse).

| J | m/Λ | λ_c | mc / quad / miss | E_repr: protocol (on λ) → on λ_c | window ×3 / signed | `L` ×3: protocol → **on λ_c** | `L` signed | binding probe on λ_c (Δ ± jk) |
|---|---|---|---|---|---|---|---|---|
| **23** | **1.0228 (twist)** | +0.901 | 0.047 / 0.023 / 0.303 | 0.553 → 0.258 | 0.131 / 0.081 | 0.85 → **1.18** | 0.90 → 1.27 | K7 +0.086 ± 0.106 |
| 23 | 1.03 | +0.923 | 0.046 / 0.011 / 0.103 | 0.501 → 0.202 | 0.133 / 0.083 | 1.16 → **1.87** | 1.24 → 2.08 | K7 +0.067 ± 0.091 |
| 23 | 1.05 / 1.075 / 1.10 | +0.998 / +0.950 / +0.998 | ≤ 0.045 / ≤ 0.003 / ≤ 0.001 | 0.406 / 0.302 / 0.282 → 0.097 / 0.029 / 0.013 | 0.134 / 0.113 / 0.109 | 1.69 / 2.09 / 2.33 → **3.56 / 5.22 / 6.22** | → 4.32 / 6.74 / 8.21 | K7 ≤ +0.032 |
| 27 (control) | 1.0267 (twist) | +0.998 | 0.058 / 0.064 / 0.136 | 0.089 → 0.031 | 0.263 / 0.149 | 1.64 → **1.81** | 2.02 → 2.29 | half1 +0.010 ± 0.013 |
| 27 | 1.03 – 1.10 | +0.96 … +1.00 | | 0.088 … 0.070 → 0.031 … 0.018 | | 1.93 … 3.21 → **2.17 … 3.88** | → 2.87 … 5.63 | |
| **31** | **1.0306 (twist)** | +0.881 | 0.070 / 0.213 / 0.167 | 0.377 → 0.280 | 0.427 / 0.221 | 0.70 → **0.76** | 0.85 → 0.94 | **half1 +0.093 ± 0.014** (half0 +0.092 ± 0.019) |
| 31 | 1.04 | +0.886 | 0.062 / 0.094 / 0.014 | 0.309 → 0.232 | 0.410 / 0.215 | 1.00 → **1.09** | 1.29 → 1.45 | half1 +0.077 ± 0.012 |
| 31 | 1.05 / 1.075 / 1.10 | +0.925 / +0.942 / +0.979 | | 0.273 / 0.185 / 0.134 → 0.196 / 0.115 / 0.067 | 0.403 / 0.363 / 0.338 | 1.22 / 1.55 / 1.90 → **1.35 / 1.75 / 2.18** | → 1.90 / 2.61 / 3.43 | |

**Verdict.**
* **J = 23 twist edge: [NE + M] → [NE±] (thin, `L = 1.18` ×3 envelope; 1.27 signed).** All of J = 23 is now [NE±] (1.03Λ:
  1.87). The binding terms are now the missing-GK-node envelope (x = 0.941 failed and nodes above 0.945, 0.303) and the residual K7
  probe (0.258, a 0.8σ number ×3). Batch 5's second GK window on [0.9, Λ/m] and the bulk mid-points address exactly these two.
* **J = 31 twist edge: stays [NE + M]** (`L = 0.76` ×3, `0.94` signed). Its representation term is not noise:
  * both every-other-node half-samples sit `+0.09` above the full row, agree with each other, and are 5–6σ from it;
  * degrees K4–K8 agree on λ_c within 0.02;
  * so the half-sample probe measures the under-resolution of the *half* grid (node spacing 0.05 against an oscillation period
    ≈ 0.05–0.1 at J = 31). No reweighting removes it.

  The other blockers are the window (×3 envelope 0.427, of which ν < 20 at the 10⁴ convention 0.189; n1 cannot go below ν ≈ 20 and n3 is
  broken for C3H) and the edge quadrature/missing nodes (0.213 + 0.167). Even with E_repr → 0 the twist edge would give `L ≈ 1.0`.
* **J = 31 at 1.04Λ: [NE + M] (borderline 1.00) → [NE±]** (`L = 1.09` ×3, 1.45 signed); 1.05–1.10Λ were already [NE±] and strengthen
  (1.35–2.18).
* **J = 27** (control): `L` 1.64 → 1.81 at the twist edge; nothing gets worse.
* **Would a referee accept it?** The λ_c version should be acceptable. It removes a double count of a window that is bounded elsewhere,
  keeps every probe and the ×3 factor, and is applied uniformly. The WLS version is a defensible estimator, but here it neither
  helps much nor addresses the mechanism. A χ²-based dropping of K7 is defensible statistics, but it looks like choosing the probe after
  seeing it. It should not be the basis of a grade.
* **Caveat.** J = 23's twist edge is thin (`L = 1.18`) and rests on the one remaining K7 deviation being counted at ×3 although it is
  consistent with zero. A sceptical reader can still call it [NE + M]. If so, batch 5's bulk mid-points (59 nodes, so K7 has 2.5× more
  data per parameter) are the decisive test.

---

## 4. Batch 5: the external job list for J ≥ 17 (`ads4/ext/jobs_batch5.txt`, `nf_batch5.py`)

**What is still model-level at J ≥ 17** (Δ_g = 1000; after HIGHSPIN, SMALLNU, EXT_FOLD and §3 here):

| J | current state | what blocks it | what batch 5 adds |
|---|---|---|---|
| 17, 19, 25 (odd) | [NE + M]: own-mass K₂ row (ee1000, Δ = 1000 + J) and 5 equispaced direct edge nodes x = 0.8 … 0.9 (HC / NEARGAP2); no direct bulk | the K₂ λ-level step (0.3–0.6 N_s at these J) does not close ×3; the edge is a spline through 5 nodes per ≈ 2 periods (aliasing, HIGHSPIN §2); J = 25's reference zero x₀ = 0.933 lies *inside* the edge | full direct bulk row + GK15 edge + n1 window (J = 25: + flanking nodes x₀ ∓ 0.0125) |
| 18, 20, 22, 24, 28, 32 (even) | [NE + M]: K₂ rows at 336–1600 (none at 1000 + J); edge nodes on non-own rows | K₂ truncation; 1/m interpolation; J = 20's 1600-row edge noise; J = 20, 32 have no budget at all | full 38-node own direct row + GK15 edge + window (+ reference-zero flanks for J ≥ 24) |
| 26, 30 (even) | [M]: no K₂ row exists (J-interpolated) | everything | as above; the analysis also needs the even ⟨TTO⟩ shadow at (1000 + J, J), computed locally (`nf_rows.py` does that) |
| 23 (twist edge) | [NE±] thin after §3 (`L = 1.18`) | missing GK nodes x ≥ 0.941 (0.303); residual K7 probe (0.258) | bulk mid-points x = 0.1875 … 0.7875 (25) + a second GK15 window on [0.9, Λ/m] |
| 31 (twist edge) | [NE + M] (`L = 0.76`) | half-sample probe +0.09 (real under-resolution of the half grid); window ×3 0.427; edge quadrature 0.213 + missing 0.167 | as J = 23 |

J = 21, 27 and 29 are [NE±] already. J = 29's remaining sign-free gap is the ν < 20 window, which n1 cannot reach and n3 cannot
evaluate (C3H broken), so it is not in the batch.

**Reference zero at Δ_g = 1000 (the check asked for).** With Δ_ref = 2·10⁴ the first zero of the reference's leading partial wave sits at
x₀ = 2.76 (J = 9), 2.27 (10), 2.20 (11), 1.84 (13), 1.64 (14), 1.58 (15), 1.44 (16). That is outside every row, so Δ_ref = 2·10⁴ is
safe for all of §1–2. For J ≥ 17 on the own rows Δ = 1000 + J it moves in:
* x₀ = 1.39 / 1.28 / 1.24 / 1.15 / 1.05 / 1.02 at J = 17 / 18 / 19 / 20 / 22 / 23, outside the row;
* 0.960 / 0.933 / 0.886 / 0.822 at J = 24 / 25 / 26 / 28, in the edge window;
* 0.767 / 0.748 / 0.718 at J = 30 / 31 / 32, in the bulk.

Δ_ref is hard-wired (`spinfunc/hc_points.py`, `ext/point.py`). Changing it needs a new point script, which was not created here (no
existing file may be edited). Batch 5 instead follows batch 4: two flanking nodes x₀ ∓ 0.0125, so that the node nearest x₀ can be dropped
(`hs_core.refbad` / `ZWIN`). If the reference is ever made configurable, Δ_ref = 5·10⁴ puts x₀ ≥ 1.8 for every J ≤ 32 at 1000.

**Content** (one job per ν, n1 only, `HC_ADAPT=1, HC_M=HC_MREF=11, L1_DEG=6`, seed irrelevant):
* *odd J = 17, 19, 25:* bulk x = 0.05 … 0.15 (step 0.0125) + 0.175 … 0.775 (step 0.025) (34 nodes, = job A); x = 0.025, 0.0375; the
  GK15 nodes of [0.8, Λ/m_row] with x ≤ 0.945 (10–11); n1 window ν = 20, 25, 32.5, 40, 47; J = 25 also its reference-zero flanks.
* *even J = 18 … 32:* the 38-node row x = 0.025 … 0.95 (HIGHSPIN §4); the same GK15 edge nodes (11), so that the even edge can be
  integrated without the spline aliasing HIGHSPIN §2 found; the n1 window; reference-zero flanks for J ≥ 24.
* *J = 23, 31:* bulk mid-points (25) and the second GK15 window on [0.9, Λ/m] with x ≤ 0.945 (8 / 9 nodes).

**Size and time.** **666 jobs**: J = 17: 51, 18: 54, 19: 52, 20: 54, 22: 54, 23: 33, 24: 56, 25: 54, 26: 56, 28: 56, 30: 56,
31: 34, 32: 56. The one-node pilot (below) took 2.8 min, and bulk nodes at J = 13 take 1.3–1.4 min here. At ≈ 2.6 min per job this is
**≈ 29 core-h, ≈ 1.8 h on 16 workers** (0.7 GB each). The file is ordered by J, so the lowest spins come back first.

**Pilot** (`nf_b5pilot.py`, `ext/point.py` run unchanged behind the write guard): job `N5_J24_edge`, ν = 908.1325 (x = 0.8884,
ν/Δ = 0.887). It took 169 s (41 + 40 s for the state, 44 + 43 s for the reference). The largest replica relative error was
6.8·10⁻³ (state) / 1.4·10⁻⁴ (reference), so it is valid (threshold 0.2). The output is in `cache/b5pilot/`, in the same layout as
`run_ext.py`.

**How the returned points will be analysed** (not run here):
* odd J: `hs_odd.py` + `sn_odd.py` via an `ef_merge`/`ef_oddrun`-style merge, with E_repr on λ_c (§3);
* even J: `nf_state.py` (the `ng_bind` own-row protocol of §1) with the GK edge nodes injected as extra nodes. A GK-quadrature even-J edge
  (an even `hs_odd.value`) would be new code;
* J = 26, 30: first `nf_rows.py even EE,C3H J 1000+J 0.5` (one K₂ node: it creates the even shadow and the row scaffold).

Expected outcome, by analogy with J = 21 and 27 (direct bulk + GK edge + n1 window, `L` 1.2–1.8):
* J = 17–22 should close;
* J = 24–32 are thinner, because the window weight grows ≈ J³ (SMALLNU) and is then limited by the ν < 20 convention, which no n1 job can
  remove.

## 5. Suggested record changes (not applied)

1. **POS_ERRORS / HIGHSPIN §6 / NEARGAP2 §0 inventory:** "J = 7, 9–11, 13–15 near gap [NE + M]" → **J = 9, 10, 11, 13, 14, 15
   [NE±]** (own direct rows; `L ≥ 16.6 / 5.6 / 16.0 / 6.8 / 1.52 / 1.75`; twist-edge values `+0.676 / +0.589 / +0.349 / +0.408 / +1.001
   / +0.694`).
2. **EXT_FOLD §3.2 / UNIFIED_THEOREM box:**
   * J = 23 twist edge [NE + M] → [NE±] (thin, `L = 1.18`), with the representation term evaluated on λ_c (§3).
   * J = 31 at 1.04Λ → [NE±] (`L = 1.09`).
   * The phrase "a noise-aware term would give L ≈ 1.52" should be read with §3: the mechanism is double counting of the extrapolated
     window, not noise per se. A WLS fit gives only 1.01.
3. **hs_odd.py's E_repr** (code, not changed) is evaluated on λ while the reported `L` uses λ_c. Future odd-J budgets should use the same
   window cut in the probes (`nf_reprcut.py` shows how).
4. **Referee report 5, M1:** "near-gap J = 9–11, 13–15 … model level" is superseded at Δ_g = 1000. The open near-gap list at 1000 is now
   J = 17–20, 22, 24–26, 28, 30, 32 and the J = 31 twist edge (batch 5).
5. **EXT_FOLD §3.5 calibration:** "0.6–8 % of the step over 7 states" → "0.4–8 % over 13 states (J = 9–16, 21–31), always favourable".

## 6. Caveats

1. The direct n1 integrals are the project's reference standard [NE], not certified. Every [NE±] here, as in NEARGAP2 / HIGHSPIN, takes
   "direct = crystal integral" as input.
2. `G_low = 10⁴ c₂` on ν < 50 and `G_edge` beyond the last valid node are conventions, with measured headroom ≥ ×100 (small ν) and ×10
   (edge).
3. For J ≥ 13 the base is the hybrid quadrature. The certified flat floor minus Lemma C does not close there; the same holds at J = 12, 16.
4. The masses are a grid (twist, 1.02, 1.03, 1.05Λ); between them H_cont applies as before. Masses above 1.05Λ were not budgeted here.
   For J = 11, 13, 15 they are held at the own row (§2.3).
5. The K₂ rows of J = 9, 10, 11, 13, 14 are coarse (10 nodes). They enter only the K₂-only comparison and the calibration, never the
   [NE±] values.
6. J = 23's closure (`L = 1.18`) still counts a 0.8σ K7 deviation at ×3; the [NE±] is thin (§3).
7. Batch 5 is prepared and piloted (1 node), not run. Its outcome for J ≥ 24 is uncertain: the ν < 20 window is beyond n1, and n3's C3H
   evaluator is broken.

## 7. Reproduce

```bash
cd ads4/code/ngfill            # every compute command under: ulimit -v 3000000; OMP_NUM_THREADS=1 (nf_lane.sh / nf_pause.sh do this)
python3 nf_points.py 13 1013 dense            # full direct own row x = 0.025 ... 0.95 (38 nodes; ~45-60 min); likewise 9/1009,
                                              # 10/1010, 11/1011, 14/1014, 15/1015
python3 nf_rows.py odd EE,C3H 13 1013 0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9     # coarse own K2 row (also 9, 11 odd; 10, 14 even)
python3 nf_shadow.py 9 1009 11 1011           # odd shadows without a row (only needed if the direct row is analysed first)
python3 nf_state.py 13                        # -> cache/nf_state_J13.{txt,json} (twist, 1.02, 1.03, 1.05; K2-only, direct, NU_LOW=25)
python3 nf_nodes.py 13                        # -> cache/nf_nodes_J13.txt
python3 nf_summary.py                         # -> cache/nf_summary.txt (§0/§2 tables, held-row estimate, K2 calibration)
python3 nf_repr.py 23 ols 1.03,1.05,1.075,1.1 # hs_odd + sn_odd unchanged on extfold's merged points -> cache/repr_ols/ (EXT_FOLD repro)
python3 nf_repr.py 23 wls 1.03,1.05           # the WLS variant -> cache/repr_wls/
python3 nf_reprcut.py 23 1.03,1.05,1.075,1.1  # probes split by x-window -> cache/nf_reprcut_J23.{txt,json}  (also 27, 31)
python3 nf_reprsum.py                         # -> cache/nf_reprsum.txt (§3 table)
python3 nf_batch5.py                          # -> ../../ext/jobs_batch5.txt
python3 nf_b5pilot.py N5_J24_edge 24 1024 908.132512   # 1-node pilot through ext/point.py -> cache/b5pilot/
```
Exact job lists and times: `cache/lane_{A,B,C,D,E}.log`; analysis logs `cache/run_*.log`.
