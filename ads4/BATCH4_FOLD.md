# BATCH4_FOLD: the batch-4 external direct rows folded into the Δ_g = 2000 near-gap budgets (J = 12, 16, 21, 27, 29; seam J = 20, 23)

*25 Nov 2026, batch-4 fold agent. Unreviewed, not committed. No existing file was edited; `prl/` was only read (§6 gives the replacement
texts for the `% BATCH4` markers, to paste by hand). `LADDER_ALLJ.md` / `code/ladder/` were not touched.*
* *New code:* `ads4/code/b4fold/` (write guard `b4_guard.py`: every write outside `b4fold/cache/` raises). *Caches:* `b4fold/cache/`, with a
  `.gitignore` that ignores everything but itself.
* *Compute:* ≈ 1.4 CPU-h, at most two computing processes, `ulimit -v 3000000`, one thread, no bytecode written.
* *Input:* `ext/results/ext_results_batch4.json` (3380 records; the 690 batch-4 records `N4_*`, `N4S_*` were decoded with `ext/unpack.dec`
  into `b4fold/cache/unpacked/`, the rest of the file — earlier batches — was not used). **All 345 jobs (247 main + 98 seam) are present.**

Grades and conventions are BATCH5_FOLD's: **[NE±]**, [NE + M], [M], [I]; `L = value / budget`, ×3 on model terms, 3σ on Monte Carlo;
"value" is λ_c/N_s (the windows ν < 50 and x > x_last removed from the central value and bounded separately). Λ = Λ(2000) = 1998.5.

---

## 0. Result

**Per-state table at Δ_g = 2000.**
* "old" = the current K₂-only entry at 2000 (supplement Table `tab:gap2000` = GAP_NEARGAP for J = 12, GAP_SEAM for J = 16, 21, GAP_SEAM2
  for J = 27, 29): value / sign-free L (λ-level L). J = 20, 23 had no entry at 2000.
* "new": λ_c/N_s (full λ in brackets where it differs by ≥ 0.001) ± jackknife σ; budget sign-free, ×3 measured window envelope,
  representation/interpolation term on λ_c; **L** (alternative in brackets: λ-level for the own-row route, signed window for hs_odd).
* "tw" = 2000 twist edge (Δ = 2000 + J, the own row's mass); "tw₁₀₀₀" = the chart-A image of the 1000 twist edge; "min" = the 2000 hybrid
  minimum (`b4_hscan`, GAP_SEAM2). "(h)": held (§2.3).

| J | m/Λ | route | old (K₂ only) | **new** λ_c/N_s ± jk | budget | **L** (alt) | grade old → new |
|---|---|---|---|---|---|---|---|
| 12 | 1.0060 tw | own row | +0.4094 / 0.7 (5.4) | **+0.4095** ± 0.0002 | 0.026 | **15.9** | [NE±] λ-level → **[NE±]** |
| 12 | 1.0119 tw₁₀₀₀ / 1.02 / 1.05 (h) | own row | +0.8496 / 1.1 (8.0); — / 1.3 (11.6); — / 1.6 (18.2) | **+0.8499 / +0.9996 / +0.6868** | 0.027 / 0.018 / 0.002 | **31.8 / 54 / 301** | [NE±] λ-level → **[NE±], held** |
| 16 | 1.0080 tw | own row | +0.9997 / 0.4 (4.1) | **+0.9996** ± 0.0007 | 0.214 | **4.68** | [NE±] λ-level → **[NE±]** |
| 16 | 1.0159 tw₁₀₀₀ (h) | own row | +0.9399 / 0.4 (5.1) | **+0.9448** ± 0.0015 | 0.073 | **13.0** | [NE±] λ-level → **[NE±], held** |
| 16 | 1.02 / 1.025 min / 1.03 / 1.05 (h) | own row | — / 0.4 (4.6); not run; — / 0.5 (5.1); — / 0.6 (6.4) | **+0.7713 / +0.7271 / +0.7531 / +0.9998** | 0.043 / 0.026 / 0.022 / 0.010 | **18.1 / 27.7 / 34.0 / 101** | [NE±] λ-level → **[NE±], held** |
| **20** | **1.0100 tw** | own row (seam) | not computed | **+0.9651** (+0.9672) ± 0.0084 | 0.714 | **1.35** | — → **[NE±]** |
| 20 | 1.015 min / 1.0198 tw₁₀₀₀ / 1.03 / 1.05 (h) | own row | not computed | **+0.8179 / +0.8948 / +0.9988 / +0.9991** | 0.369 / 0.275 / 0.061 / 0.101 | **2.22 / 3.25 / 16.3 / 9.9** | — → **[NE±], held** |
| **21** | **1.0105 tw** | own row (odd) | +0.9973 / 0.1 (1.7) | **+1.0118** ± 0.0021 | 0.772 | **1.31** (λ-lvl 1.89) | [NE±] λ-level, thin → **[NE±]** |
| 21 | 1.0208 tw₁₀₀₀ (h) | own row (odd) | +0.9979 / 0.2 (2.3) | **+1.0101** ± 0.0020 | 0.405 | **2.49** (4.24) | [NE±] λ-level → **[NE±], held** |
| 21 | 1.03 min / 1.05 / 1.075 / 1.10 (h) | own row (odd) | — / λ 2.7 … 3.9 | **+0.8909 / +1.0081 / +0.9681 / +0.9575** | 0.217 / 0.107 / 0.057 / 0.034 | **4.10 / 9.4 / 17.1 / 27.8** | [NE±] λ-level → **[NE±], held** |
| **23** | **1.0114 tw** | hs_odd/GK | not computed | **+1.0446** (+1.0449) ± 0.0074 | 0.738 | **1.41** (signed 1.42) | — → **[NE±]** |
| 23 | 1.0228 tw₁₀₀₀ / 1.025 min / 1.03 / 1.05 / 1.075 / 1.10 | hs_odd/GK | not computed | **+0.9322 / +0.9237 / +0.9515 / +1.0206 / +0.9651 / +1.0106** | 0.220 / 0.190 / 0.153 / 0.080 / 0.042 / 0.037 | **4.24 / 4.86 / 6.24 / 12.8 / 22.9 / 27.5** | — → **[NE±]** |
| **27** | **1.0134 tw** | hs_odd/GK | +0.9670 / 0.06 (0.46) | **+1.0028** (+1.0041) ± 0.0072 | 0.834 | **1.20** (signed 1.22) | [NE + M] → **[NE±]** |
| **27** | **1.0175 min** | hs_odd/GK | +0.9098 / 0.07 (0.56) | **+0.9443** (+0.9455) ± 0.0067 | 0.441 | **2.14** (2.19) | [NE + M] → **[NE±]** |
| 27 | 1.0267 tw₁₀₀₀ / 1.03 | hs_odd/GK | +0.9983 / 0.08 (0.77); +0.9984 / 0.08 (0.85) | **+1.0365 / +1.0367** | 0.171 / 0.130 | **6.05 / 7.97** | [NE + M] → **[NE±]** |
| 27 | 1.04 / 1.05 / 1.075 / 1.10 | hs_odd/GK | +0.992 … +0.999 / λ 1.07 … 1.79 | **+1.0298 / +0.9970 / +1.0311 / +1.0013** | 0.082 / 0.070 / 0.059 / 0.052 | **12.6 / 14.3 / 17.4 / 19.4** | [NE±] λ-level, thin → **[NE±]** |
| **29** | **1.0144 tw** | hs_odd/GK | +0.9143 / 0.05 (0.16) | **+0.9407** (+0.9413) ± 0.0072 | 1.319 | **0.71** (signed 0.72) | [NE + M] → **[NE + M]** |
| **29** | **1.015 min** | hs_odd/GK | +0.9133 / 0.05 (0.16) | **+0.9404** (+0.9409) ± 0.0072 | 1.220 | **0.77** (0.78) | [NE + M] → **[NE + M]** |
| 29 | 1.0286 tw₁₀₀₀ / 1.04 / 1.05 / 1.06 | hs_odd/GK | +0.9960 / 0.06 (0.25) … +0.9983 / 0.09 (0.70) | **+1.0357 / +1.0087 / +1.0261 / +1.0403** | 0.337 / 0.175 / 0.134 / 0.113 | **3.08 / 5.78 / 7.65 / 9.18** | [NE + M] → **[NE±]** |
| 29 | 1.075 / 1.10 | hs_odd/GK | +0.9867 / 0.11 (1.01); +0.9926 / 0.14 (1.42) | **+1.0262 / +1.0290** | 0.094 / 0.071 | **10.9 / 14.6** | [NE±] λ-level, thin → **[NE±]** |

**Count at Δ_g = 2000** (near-gap spins computed at 2000; a spin counts as error-controlled when its 2000 twist edge and every computed mass are
[NE±], any suffix, as in BATCH5_FOLD):
* *Before:* 20 spins computed (J = 0, 2, 4–17, 19, 21, 27, 29); sign-free [NE±] at every mass: J = 0, 2, 4–11, 13 (**11**); λ-level only:
  J = 12, 14–17, 19, 21 (8); [NE + M] at the twist edge: J = 27, 29.
* **Now:** 22 spins computed (+ J = 20, 23). **Sign-free [NE±] at every computed mass: J = 0, 2, 4–13, 16, 20, 21, 23, 27 (17)**; λ-level only:
  J = 14, 15, 17, 19 (4, unchanged K₂ entries); **[NE + M]: J = 29 at its 2000 twist edge and hybrid minimum only** (L = 0.71 / 0.77;
  [NE±] from 1.0286Λ up). Error-controlled spins (any suffix): 18 → **21 of 22**.
* Stricter readings: no signed-window-only state is used (none needed: every L above is on the ×3 envelope); dropping held masses changes
  nothing relative to before (every 2000 mass above the twist edge sits above every row of its J, in the old K₂ entries as well, §2.3).

**Negative findings and cautions (stated first).**
1. **No state is negative.** Lowest new values: J = 12 at its twist edge +0.4095 (the hybrid's own dip, +0.4090), J = 16 at its hybrid minimum 1.025Λ +0.727 (hybrid +0.724), J = 20 at its hybrid
   minimum +0.818 (hybrid +0.763), J = 21 at 1.03Λ +0.891 (hybrid +0.883), J = 23 at 1.025Λ +0.924, J = 29 at its twist edge +0.941.
2. **J = 29's 2000 twist edge (and its hybrid minimum 1.015Λ) stays [NE + M]**, L = 0.71 / 0.77 (×14 / ×5 better than K₂-only). Binding
   term: the five GK nodes x > 0.945 that n1 cannot reach, covered by the light-cone envelope 3A((1 − x_l)/(1 − x))² (E_miss = 0.97), then
   representation 0.15 and quadrature 0.15. The signed window does not help (the window is 0.025 at 2000). Under EDGE_R2's calibrated
   pointwise factor c = 1.21–2.06 in place of 3 (protocol change, not used for the grade), L = 0.93–1.27.
3. **Even-J external records are in a different ⟨TTO⟩ structure basis** (raw ‖E_loc − E_ext‖/‖E‖ = 0.31 / 1.2 / 1.1 at J = 12 / 16 / 20,
   state row). Every even number above uses the repaired records (§1.2); the repair reproduces local recomputations to ≤ 2·10⁻⁶ (fitted) and
   1.5·10⁻⁵ (unfitted node), and Γ to ≤ 0.25 c₂ against σ = 18–85 c₂.
4. **hs_odd's shape fit is ill-conditioned at J = 21 on the 2000 row** (as BATCH5 found at J = 17, 19): the degree-6 basis has condition
   number 7.9·10⁵ (above the 4.6·10⁵ maximum where hs_odd was used at 1000), degree 7 1.3·10⁸, degree 8 1.7·10¹¹; the K7 probe gives
   +0.27 and K8 −525 against K6 +1.04, so the protocol's E_repr = 2.3 and L = 0.40. As BATCH5 did for J = 17, 19, **J = 21 is folded by
   the own-row route** (§2.1). hs_odd's K4–K6 values (+1.011 … +1.041) agree with the own row (+1.012) to ≤ 0.03. At J = 23, 27, 29 the
   degree-6 condition numbers are 4.5·10⁵ / 1.9·10⁵ / 1.8·10⁵ and the probes behave (largest |Δλ_c| 0.023 / 0.014 / 0.051).
5. **J = 23: one GK edge node failed** (x = 0.9139, replica error 0.27). With 9 valid GK nodes hs_odd falls back to its spline edge
   ('espl', L = 4.95 at the twist edge), which carries **no** sign-free term beyond the last node. We use instead the GK edge with the failed
   node and the five unreachable ones in E_miss (B5_GKMIN = 9, the switch BATCH5 introduced): L = 1.41. The more conservative value is graded.
6. **Every 2000 mass above the twist edge is held on the own row** (no row above Δ = 2000 + J exists at these J); the K₂-only entries held the
   same row. For the hs_odd masses this is the protocol, as at 1000.
7. **K₂ against direct at 2000 (GAP_NEARGAP caveat 1):** the direct values are higher than the K₂ values on the same rows by +0.0002 … +0.026
   (every shift favourable), ≤ 3 % of the K₂ λ-level truncation charge 3|K₂ − K₁|. Pointwise, Γ_dir − Γ_K₂ has rms 70–137 c₂ on x ≤ 0.775
   with χ²/n = 2.2–3.0 against the jackknife: the K₂ truncation is visible node by node but averages out in λ (§4).

---

## 1. Inputs and the structure-basis repair

### 1.1 Points (`b4_merge.py` → `cache/pts/all`, log `cache/b4_merge.txt`)
* 345 jobs, 690 records (state + reference, EE + C3H per job); no key collisions; every node of `jobs_batch4.txt` present.
* Failed by the replica-error rule (> 20 %, `hs_core.failed_x`): **only J = 23 x = 0.9139** (GK edge). Unlike the 1000 rows, x = 0.95 does not
  fail on the even 2000 rows (ν/Δ = 0.949).
* GAP_NEARGAP's two local pilot jobs (J = 29, ν = 20 and x = 0.9306) reproduce the external records to 1·10⁻¹⁴ – 4·10⁻¹⁰ (odd J: no basis
  ambiguity; n1 deterministic).
* Shadows: gapnear's (2012, 12), (2016, 16), (2021/2027/2029, odd) used as the K₂ rows used them; new (2020, 20) and (2023, 23, 6) computed
  with the project's functions (`b4_shadow.py`), which reproduce gapnear's (2012, 12) and (2021, 21, 6) bit for bit. The even seam row 2020
  has no K₂ row: 5 K₂ EE state records (x = 0.2–0.8, `b4_k2.py`) feed the l8 helicity basis, as BATCH5's `b5_k2`.

### 1.2 Even J: external records mapped into the local structure basis (`b4_points.py`, `b4_basis.py`, `b4_basischk.py`)
BATCH5_FOLD §1.2's anchored O(2)+reflection fit, unchanged: one local recomputation per even row (x = 0.3, EE, state and reference,
ext/point.py's code path, 16–28 s per point here), G and H fitted on the 8 replicas + mean, every record mapped replica by replica.

| J | raw ‖E_loc − E_ext‖/‖E‖ (state / ref) | O(2) fit residual (state / ref); general 2×2 | 2nd node, unfitted (x = 0.6, EE + C3H; state / ref) | Γ(mapped) − Γ(local) at the local nodes |
|---|---|---|---|---|
| 12 | 0.31 / 0.21 | 1.0·10⁻⁶ / 1.6·10⁻⁹; 8.4·10⁻⁷ / 9.1·10⁻¹⁰ | — | ≤ 0.25 c₂ (σ 18–24) |
| 16 | 1.22 / 0.76 | 1.1·10⁻⁸ / 1.8·10⁻⁸; 9.9·10⁻⁹ / 1.7·10⁻⁸ | 1.5·10⁻⁵ / 2.5·10⁻⁷ | ≤ 0.09 c₂ (σ 20–85), C3H ≤ 0.04 |
| 20 | 1.08 / 0.003 | 3.2·10⁻⁸ / 6.0·10⁻¹²; 3.2·10⁻⁸ / 3.8·10⁻¹² | — | ≤ 0.004 c₂ (σ 53–68) |

* After the mapping the rows are helicity-diagonal: |c_off| ≤ 5.4 c₂ below x = 0.92 (J = 12, 16, 20).
* The reference (20000, 20) needs the same transformation as in batch 5 (up to an overall sign), i.e. the external machine's basis is stable
  across batches.

## 2. Protocols

### 2.1 Routes (all BATCH5_FOLD's, with Δ_g a parameter)
* **Even J = 12, 16, 20** (`b4_state.py` = b5_state with Dg = 2000; `neargap2/ng_bind.run_state` unchanged; `b4_window.py` → `b4_sneven`
  = `smallnu/sn_even.main` with Λ(2000), FULL normalisation): direct own row x = 0.025 … 0.95 (+ rz flanks, + ν = 20–47 small-ν nodes)
  injected into the own row (J = 12, 16: over gapnear's K₂ row 2000 + J, all of whose nodes are overridden; J = 20: direct only). Budget:
  trunc (0 on the own row), interpolation (on λ_c), the flat edge 10 × max|Γ| of the last four nodes beyond x_last, 3σ MC, the ν < 50 window
  replaced by the measured ×3 envelope.
* **Odd J = 23, 27, 29** (`b4_odd.py` = b5_odd with `b4_hsodd` = hs_odd and `b4_snodd` = sn_odd, Dg = 2000): direct bulk shape fit (K = 6)
  + GK15 edge on [0.8, Λ/m_row] with the missing nodes under 3A((1 − x_l)/(1 − x))², representation on λ_c, measured window.
* **Odd J = 21:** own-row route (§0 item 4), = BATCH5's route for J = 17, 19 (`b4_state.py` + the odd branch of `b4_window.py`).
* Δ_g = 2000 redirections (re-typed from `gapnear/gn_bind.py`, in-process): gapnear's 2000 K₂ rows in the row loader (Δ_ref = 2·10⁴, the
  reference the direct points use; GAP_SEAM's 5·10⁴ rows are not mixed in), the min-norm gauge factor of `gapscale/cache/weights_2000.npz`,
  `sn_window.DG = 2000`.

### 2.2 Choices made here (all fixed before grading; none changes a sign)
1. **Bulk lower end x = 0.025 for hs_odd** (hs_odd: 0.05). At 1000, x = 0.05 is ν = 50, where the separately bounded window starts; at 2000
   it is ν = 100, so batch 4 added the 'blo' nodes x = 0.025, 0.0375 (GAP_NEARGAP §4). With XB = 0.05 (`odd_xb05`) the twist-edge values
   move by ≤ 0.011 and L by ≤ 0.06 (J = 27: 1.14 vs 1.20; J = 29: 0.71 vs 0.71; J = 21 stays ill-conditioned, 0.39).
2. **Reference-zero rule** (BATCH5 §2.2 item 3) on the even rows and on the odd own row: direct nodes with |x − x₀(Δ_ref = 2·10⁴)| < 0.012
   dropped (J = 12: 0.95; 16: 0.725; 20: 0.575; 21: 0.575), **and the K₂ node at the same x is not allowed to fall back into the injected row**
   (`b4_core._drop_own`; without this the J = 16 row carried a non-finite K₂ node at x₀). hs_odd keeps its ZWIN = 0.002 rule on the GK nodes
   only, as at 1000.
3. **Hampel exemption for the injected direct rows** (BATCH5 §2.2 item 2), unchanged.
4. **J = 23: GK edge with GKMIN = 9** (§0 item 5); **J = 21: own-row route** (§0 item 4).

### 2.3 What binds
| state | binding terms (N_s) |
|---|---|
| J = 12 tw | flat edge x > 0.925 0.024 (total 0.026) |
| J = 16 tw | flat edge x > 0.95 (G = 1.2·10⁵) 0.203 |
| J = 20 tw | flat edge x > 0.95 (G = 2.8·10⁵) 0.536; interpolation 0.143 |
| J = 21 tw (own row) | flat edge x > 0.933 (G = 8.6·10⁴) 0.405; K₂ truncation on the row's K₂ nodes 0.825–0.9 between the GK nodes 0.307; interpolation 0.046 |
| J = 23 tw | E_miss 0.635 (5 unreachable + 1 failed GK node) |
| J = 27 tw | E_miss 0.733 |
| J = 29 tw | E_miss 0.966; repr 0.154; quad 0.152 |

The window ν < 50 costs ≤ 0.025 everywhere (×1/16 of 1000, chart A: E_low ∝ Λ⁻⁴).

## 3. Odd-J diagnostics (`b4_cond.py`, variants)
* Degree scan at the twist edges (K = 4 … 8; condition number; λ): J = 21 4.2·10³ … 1.7·10¹¹, +1.015 / +1.011 / +1.041 / +0.268 / −525;
  J = 23 3.6·10³ … 3.0·10⁹, +1.032 / +1.034 / +1.045 / +1.051 / −0.68; J = 27 3.3·10³ … 1.8·10⁸, +1.016 / +1.017 / +1.004 / +0.996 / +1.27;
  J = 29 3.3·10³ … 2.8·10⁷, +0.964 / +0.957 / +0.941 / +0.958 / +1.03. (The protocol's probes are K4, K5, K7 and two half-samples; K8 is
  diagnostic only.)
* hs_odd at J = 21 for the record (not graded): twist edge +1.041, L = 0.40 (E_repr 2.32); 1.05 / 1.075 / 1.10Λ L = 1.35 / 3.26 / 9.34.
* **Reference-zero rule applied to the odd bulk too** (`B4_RZODD=1`, diagnostic; drops x = 0.575 / 0.525 / 0.425 / 0.40 at J = 21 / 23 /
  27 / 29): twist-edge L = 0.44 (J = 21, hs_odd) / 1.10 (27) / 0.68 (29), against 0.40 / 1.20 / 0.71; above the twist edges L moves by
  ≤ 15 %. No grade changes. (J = 23's variant ran with the default GKMIN, i.e. on the espl edge: 6.22 vs 4.95.)
* **Bulk from x = 0.05** (`B4_XB=0.05`): §2.2 item 1.
* **Own-row route at J = 29** (cross-check, `state/b4_state_J29`): twist edge λ_c +0.965 (hs_odd +0.941), but L = 0.26: the flat edge term
  (G = 2.3·10⁵ beyond x = 0.931) is 1.85 and the K₂ truncation on the own row's K₂ nodes 0.825–0.9 is 1.82. At 1.0286Λ L = 0.69 (λ-level 1.77)
  against hs_odd's 3.08. The hs_odd route is the designated one for J = 29 and the better one; neither closes the twist edge.
* The published-rows-only chain (no 2000 + J row; ng_bind's `lam_model`, printed in the state logs) gives −0.91 … −0.16 for J = 29 at 2000:
  the rows below 2000 do not represent the 2000 states at this spin. It is a diagnostic of the old H_M hold, not a state value; every graded
  value uses the own row.

## 4. K₂ against direct at Δ_g = 2000 (`b4_k2cmp.py` → `cache/b4_k2cmp.txt`)
| J | masses | λ_dir − λ_K₂ | K₂ charge 3\|K₂ − K₁\| | ratio |
|---|---|---|---|---|
| 12 (own row) | tw … 1.05 | −0.0003 … +0.0009 | 0.029–0.053 | ≤ 0.022 |
| 16 (own row) | tw … 1.05 | +0.0002 … +0.0030 | 0.11–0.16 | ≤ 0.020 |
| 21 (own row) | tw … 1.10 | +0.0019 … +0.0068 | 0.44–0.86 | ≤ 0.008 |
| 21 / 27 / 29 (hs_odd bulk, same edge) | tw … 1.10 | +0.003 … +0.017 / +0.012 … +0.019 / +0.018 … +0.026 | 0.43–1.03 | ≤ 0.025 |

So the K₂-only entries at 2000 were right in value (they are all within 0.04 of the direct values, the K₂ chain with Δ_ref = 5·10⁴ included),
and their truncation charge was ×30–200 too pessimistic: "calibrated" L (SMALLNU's 5–15 % of the step) was itself conservative. Chart A holds
for the direct values: Γ on 2000 + J agrees with 1000 + J to O(100 c₂) in the bulk (node tables in `cache/state/`, `cache/odd_main/`).

## 5. Dense-grid recheck of the J = 46 / 50 ++ EDGE-R warning (optional)
DENSE_PLACEHOLDER

## 6. Replacement texts for the `% BATCH4` markers
Line numbers refer to the current `prl/supplement.tex` / `prl/letter.tex`. Every `% BATCH4` comment line is deleted together with its replacement.

### 6.1 `supplement.tex` ≈ 1485–1488, Summary item (2)
Replace from `(2) At $\Delta_g=2000$ every near-gap spin` through `(Sec.~\ref{sec:allgap}.C).` and delete the `% BATCH4: update summary item (2) …` line, with:
```latex
(2) At $\Delta_g=2000$ every near-gap spin $J\le21$ except $18$, $J=23$, $27$, $29$ and the parity-odd OO states $J=4$--$10$
were computed on own-mass rows, for $J=12$, $16$, $20$, $21$, $23$, $27$, $29$ with direct (Monte Carlo) evaluations of the crystal
integral, including the seam strip (two spins also at $5000$); nothing is negative; $J\le13$, $J=16$, $20$, $21$, $23$, $27$ and the OO
states close [NE$\pm$] sign-free at every computed mass, $J=14$, $15$, $17$, $19$ on the $\lambda$-level budget only, and $J=29$ from
$1.03\Lambda$ on ([NE\,+\,M] at its $2000$ twist edge and hybrid minimum, $L=0.71$, $0.77$); the direct values exceed the $K_2$ values on the
same rows by at most $0.026$, $\le3\%$ of the $K_2$ truncation charge; for $J\le16$ every budget term scales like $\Lambda^{-2}$ or
faster; $J=18$, $22$, $24$--$26$, $28$, $30$--$32$ are not computed at $2000$ (Sec.~\ref{sec:allgap}.C).
```

### 6.2 `supplement.tex` ≈ 1583–1587, Table `tab:gap2000` (caption, comment, rows)
(a) In the caption replace `Values at $2000$ for $J\ge9$ and OO are the windowed central values of Sec.~\ref{sec:pos} ($\nu<50$ cut and budgeted separately), for $J\le8$ and $12$ the own-row $K_2$ values; the two differ by $\le0.001$ for $J\le16$ and by $\le0.03$ at $J=21$--$29$.}` with:
```latex
Values at $2000$ for $J\ge9$ and OO are the windowed central values of Sec.~\ref{sec:pos} ($\nu<50$ cut and budgeted separately), for
$J\le8$ the own-row $K_2$ values; the two differ by $\le0.001$ for $J\le16$. For $J=12$, $16$, $20$, $21$, $23$, $27$, $29$ the $2000$
entries come from direct own rows ($\Delta_{\rm ref}=2\times10^4$; even $J$ repaired into the local structure basis), with the protocol
of Sec.~\ref{sec:pos} for direct rows (sign-free budget, measured $\times3$ window; odd $J\ge23$ with the direct bulk and the
Gauss--Kronrod edge); there no bracket is given (no truncation term); mass entries above each $2000$ twist edge hold the own row.}
```
(b) Delete the two comment lines `% BATCH4: the J = 12-29 entries at 2000 are K2-only; …` and `% plus 98 optional seam jobs) replace them. …`.

(c) Rows. Replace the J = 12 rows `tw$_{2000}$` and `tw$_{1000}$` (keep the `tw$_{5000}$` row) with:
```latex
12 & tw$_{2000}=1.0060$ & --- & $+0.4095$ / $15.9$ & $+0.4086$ / $4.4$ ($31$) & [NE$\pm$]\\
12 & tw$_{1000}=1.0119$ & $+0.8528$ / $7.4$ & $+0.8499$ / $31.8$ & $+0.8482$ / $6.4$ ($47$) & [NE$\pm$]\\
```
Replace the two J = 16 rows with:
```latex
16 & tw$_{2000}=1.0080$ & --- & $+0.9996$ / $4.7$ & --- & [NE$\pm$]\\
16 & tw$_{1000}=1.0159$ & $+0.9402$ / $4.4$ & $+0.9448$ / $13.0$ & --- & [NE$\pm$]\\
16 & $1.025$ (min) & (model) & $+0.7271$ / $27.7$ & --- & [NE$\pm$]\\
```
Insert after the second J = 19 row:
```latex
20 & tw$_{2000}=1.0100$ & --- & $+0.9651$ / $1.35$ & --- & [NE$\pm$]\\
20 & $1.015$ (min) & --- & $+0.8179$ / $2.2$ & --- & [NE$\pm$]\\
20 & tw$_{1000}=1.0198$ & $+0.835$ / $1.41$ & $+0.8948$ / $3.3$ & --- & [NE$\pm$]\\
```
Replace the J = 21, 27 and 29 rows (from `21 & tw$_{2000}=1.0105$` through `29 & $1.075$--$1.10$ …`) with:
```latex
21 & tw$_{2000}=1.0105$ & --- & $+1.0118$ / $1.31$ & --- & [NE$\pm$]\\
21 & tw$_{1000}=1.0208$ & $+0.998$ / $1.20$ & $+1.0101$ / $2.5$ & --- & [NE$\pm$]\\
21 & $1.03$ (min)--$1.10$ & $L=1.7$--$6.2$ & $+0.8909$ / $4.1$; $L\le28$ & --- & [NE$\pm$]\\
23 & tw$_{2000}=1.0114$ & --- & $+1.0446$ / $1.41$ & --- & [NE$\pm$]\\
23 & tw$_{1000}=1.0228$ & $+0.926$ / $1.62$ & $+0.9322$ / $4.2$ & --- & [NE$\pm$]\\
23 & $1.025$ (min) & --- & $+0.9237$ / $4.9$ & --- & [NE$\pm$]\\
23 & $1.03$--$1.10$ & $L=2.9$--$6.4$ & $L=6.2$--$27$ & --- & [NE$\pm$]\\
27 & tw$_{2000}=1.0134$ & --- & $+1.0028$ / $1.20$ & --- & [NE$\pm$]\\
27 & $1.0175$ (min) & --- & $+0.9443$ / $2.1$ & --- & [NE$\pm$]\\
27 & tw$_{1000}=1.0267$ & $+0.976$ / $1.81$ & $+1.0365$ / $6.1$ & --- & [NE$\pm$]\\
27 & $1.03$--$1.10$ & $L=2.9$--$3.9$ (from $1.04$) & $L=8.0$--$19$ & --- & [NE$\pm$]\\
29 & tw$_{2000}=1.0144$ & --- & $+0.9407$ / $0.71$ & --- & [NE\,+\,M]\\
29 & $1.015$ (min) & --- & $+0.9404$ / $0.77$ & --- & [NE\,+\,M]\\
29 & tw$_{1000}=1.0286$ & $+0.972$ / $1.24$ (signed) & $+1.0357$ / $3.1$ & --- & [NE$\pm$]\\
29 & $1.04$--$1.10$ & $L=4.3$--$5.7$ (from $1.075$) & $L=5.8$--$15$ & --- & [NE$\pm$]\\
```
(d) Optional, in the protocol paragraph before the table, after `… ($\Delta_{\rm ref}=5\times10^4$ … (reading (iv)).`, add:
```latex
For $J=12$, $16$, $20$, $21$, $23$, $27$, $29$ the own row at $2000$ is a direct (Monte Carlo) row, as at $1000$; its values exceed the
$K_2$ values by $+0.0002$ to $+0.026$.
```

### 6.3 `supplement.tex` ≈ 1731–1734, end of *Other near-gap states*
Replace `For $J\ge14$ the sign-free $K_2$ budget alone exceeds the margin at $1000$ and $2000$ (at $2000$ only the $\lambda$-level budget closes, Table~\ref{tab:gap2000}), so a sign-free transfer needs direct rows there, as at $1000$. % BATCH4` with:
```latex
For $J\ge14$ the sign-free $K_2$ budget alone exceeds the margin at $1000$ and $2000$. With direct rows at $2000$ ($J=12$, $16$, $20$,
$21$, $23$, $27$, $29$) the sign-free budget closes at every computed mass except the $J=29$ twist edge and hybrid minimum
(Table~\ref{tab:gap2000}), and the direct values lie within $3\%$ of the $K_2$ truncation charge of the $K_2$ values, so at $2000$ the
$K_2$ rows are checked directly; beyond $2000$ a sign-free transfer for $J\ge14$ still needs direct rows, as at $1000$.
```

### 6.4 `supplement.tex` ≈ 1846–1853, coverage table rows
Replace the row `near gap $J=12$, $14$--$16$ & … $5000$; transfer as above\\` with:
```latex
near gap $J=12$, $16$ & [NE$\pm$] at $1000$ and sign-free at $2000$ (direct rows, seam included; $L\ge4.7$); $J=12$ sign-free at
$5000$; transfer as above\\
near gap $J=14$, $15$ & [NE$\pm$] at $1000$ (direct rows), $\lambda$-level at $2000$ ($K_2$ only, seam included); transfer as above\\
```
Replace the row `near gap $J=21$ & … transfer as above\\` with:
```latex
near gap $J=20$, $21$ & [NE$\pm$] at $1000$ ($J=21$ thin) and sign-free at $2000$ at every computed mass (direct own rows, seam
included; $L=1.35$, $1.31$ at the $2000$ twist edges); transfer as above\\
```
Replace the row `near gap $J=23$, $27$, $29$, $31$ & … transfer as above\\` with:
```latex
near gap $J=23$, $27$, $29$, $31$ & at $1000$ as in Table~\ref{tab:coverLS} (direct rows; $29$, $31$ signed at the twist edge);
at $2000$ (direct rows) $J=23$, $27$ [NE$\pm$] at every computed mass ($L=1.41$, $1.20$ at the twist edges), $J=29$ [NE\,+\,M] at its
twist edge and at $1.015\Lambda$ ($L=0.71$, $0.77$) and [NE$\pm$] from $1.03\Lambda$; $31$ not computed at $2000$; transfer as above\\
```
In the row `other near-gap states, $J\le32$ ($J=18$, $20$, $22$, $24$--$26$, $28$, $30$, $32$) & … ([NE$\pm$] for $J=18$, $20$, $25$; …` delete `$20$, ` in both lists (giving `($J=18$, $22$, $24$--$26$, $28$, $30$, $32$)` and `([NE$\pm$] for $J=18$, $25$; …`). Delete the comment line `% BATCH4: the rows "J = 12, 14-16", "J = 21", "J = 23, 27, 29, 31" change …`.

### 6.5 `letter.tex` ≈ 276–281
Replace from `At $\Delta_g=2000$ we` through `as the approach to flat space predicts.` and delete the two `% BATCH4` / `% (J = 27, 29 …` comment lines, with:
```latex
At $\Delta_g=2000$ we
repeated the error-controlled near-gap checks for twenty-two spins up to $29$, seven of them on direct Monte Carlo
evaluations of the crystal integral: no state is negative, the checks close for every computed state except the two lightest
states of spin $29$, for four spins only with a signed estimate of the truncation error, the direct evaluations agree with the
saddle-point ones within a few per cent of the latter's own error estimate, and up to spin 16 the error budgets shrink like
$\Delta_g^{-2}$, as the approach to flat space predicts.
```

## 7. Reproduce
```bash
cd ads4/code/b4fold        # every command via lane.sh: ulimit -v 3000000; OMP_NUM_THREADS=1; PYTHONDONTWRITEBYTECODE=1
python3 b4_merge.py; python3 b4_shadow.py; python3 b4_k2.py 20 0.2,0.35,0.5,0.65,0.8
python3 b4_points.py 12 2012 603.138342; python3 b4_points.py 16 2016 604.329728; python3 b4_points.py 20 2020 605.518769
B4_KINDS=EE,C3H python3 b4_points.py 16 2016 1208.659456
python3 b4_basis.py; python3 b4_basischk.py 12,16,20; python3 b4_hscan.py 12,16,20,21,23
python3 b4_state.py 12 1.0119407,1.02,1.05; python3 b4_window.py 12        # likewise 16 (1.0158898,1.02,1.025,1.03,1.05), 20 (1.015,1.0198235,1.03,1.05)
python3 b4_state.py 21 1.0208046,1.03,1.05,1.075,1.1; python3 b4_window.py 21
python3 b4_odd.py 27 main 1.0175,1.0266713,1.03,1.04,1.05,1.075,1.1; python3 b4_odd.py 29 main 1.015,1.0286195,1.04,1.05,1.06,1.075,1.1
B5_GKMIN=9 python3 b4_odd.py 23 gk9 1.0227639,1.025,1.03,1.05,1.075,1.1;  python3 b4_odd.py 21 main 1.0208046,1.03,1.05,1.075,1.1  (diagnostic)
python3 b4_cond.py 21,23,27,29; python3 b4_sum.py; python3 b4_k2cmp.py     # -> cache/b4_sum.txt (every number in §0), cache/b4_k2cmp.txt
```

## 8. Caveats
1. The direct n1 integrals are [NE] (the reference standard), not certified; the even-J results rest on the anchored basis repair (one local
   node per row; J = 16 also an unfitted second node).
2. Model inputs kept: 10⁴ c₂ below ν = 20; the flat edge G beyond x_last (even, own-row) or the light-cone envelope (hs_odd); the own-row hold
   above every 2000 twist edge (H_mono-A), as in the K₂-only entries.
3. Only the listed masses were evaluated; between them H_cont as before. J = 14, 15, 17, 19 at 2000 are unchanged (K₂ only, λ-level).
4. Choices introduced here: the hs_odd bulk lower end 0.025, the K₂-fallback drop at a removed own-row node, GKMIN = 9 at J = 23, the own-row
   route at J = 21. Each is argued in §0/§2, none is tuned to Γ, and the alternatives are reported.
