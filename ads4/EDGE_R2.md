# EDGE_R2: (EDGE-R) at Theorem W's gap-mass states from the batch-3 n1 rows — the ++ entry closes at J = 26, 30; the binding +− entry stays open at all 8 states

*24 Nov 2026. Unreviewed, not committed. No existing file was edited.*
* *New code:* `ads4/code/edger2/` (`run_all.sh`, ≈ 1.5 CPU-h used, at most 2 processes, `ulimit -v 3000000`, `OMP_NUM_THREADS=1`).
* *Caches:* `edger2/cache/`, whose `.gitignore` ignores everything except itself. The write guard `e2_guard.py` confines writes to that directory.
* *Input:* `ext/results/ext_results_batch3_partial.json`, decoded with `ext/unpack.dec` (no external pickle was loaded). It holds the n1 jobs of the 8 states named below; none of their n3 jobs.
* *Read-only modules:* `edger/` (er_core, er_jobs), `highspin/hs_core` → `neargap2/ng_val` → the EE chain (`spinfunc/` w2/q3/l9, ee1000), `agproof/` (ag_exact, chi.json), `p4win/`, `lsw/`.

Grades as in AG_PROOF / EDGE_R. Units: S_E,t (EDGE_R §1). cost_E = |cost₂| + |cost₃| per entry, where cost_k = ∫χ w_k Γ_k. The tolerance is AG_PROOF's ε_E tol = (1 − d_χ)H.

---

## 0. Result

**No violation was found, and none can be excluded. (EDGE-R) is closed at none of the 8 states.**
* **The ++ entry closes at (26, 1.026Λ) and (30, 1.030Λ)** [NE + M].
* **The binding +− entry is open at all 8 states.** Two things keep it open: the unmeasured edge x > x_l, and the error budget of the measured part.
* Theorem W's sampled J\* is unchanged: J\* = 26, conditional on (EDGE-R).

| state (J, m/Λ); x₁; last n1 node x_l | entry | determined part: \|D₂\| + \|D₃\| (D₂, D₃) ± budget [mc / repr / quad] | bounded remainder: edge x > x_l (C2 + C3H) + gap | tolerance | upper bound | verdict |
|---|---|---|---|---|---|---|
| **(26, 1.026)**; 0.975; 0.928 | +− | 0.27 (+0.16, −0.11) ± [0.23 / 0.41 / 0.40] | 2.16 (1.17 + 0.99) + 0.42 | **0.406** | 3.90 | **OPEN** |
| | ++ | 0.010 ± [0.06 / 0.004 / 0.000] | 0.000 + 0.001 | 0.406 | **0.075** | **closed [NE + M]** |
| **(30, 1.030)**; 0.971; 0.928 | +− | 1.80 (−1.48, +0.32) ± [2.10 / 0.86 / 0.81] | 3.52 (1.90 + 1.62) | **0.629** | 9.08 | **OPEN** |
| | ++ | 0.034 ± [0.12 / 0.06 / 0.00] | 0.000 | 0.629 | **0.22** | **closed [NE + M]** |
| **(34, 1.033)**; 0.968; 0.926 | +− | 1.83 (+0.82, −1.01) ± [3.59 / 5.10 / 2.81] | 7.66 (4.21 + 3.45) + 5.37 | **0.886** | 26.4 | OPEN |
| | ++ | 0.19 ± [1.35 / 0.47 / 0.00] | 0.001 + 0.005 | 0.886 | 2.03 | OPEN (MC, repr) |
| (33, 1.033) odd; 0.968; 0.926 | +− | 2.60 (−1.93, −0.66) ± [5.40 / 1.04 / 3.11] | 4.32 (2.44 + 1.88) | 1.77 | 16.5 | OPEN |
| (38, 1.037); 0.964; **0.899**† | +− | 7.30 (−2.07, +5.23) ± [3.88 / 11.6 / 5.44] | 37.7 (21.7 + 16.0) + 7.16 | 1.22 | 73.0 | OPEN |
| | ++ | 0.61 ± [1.77 / 0.68 / 0.02] | 0.015 + 0.014 | 1.22 | 3.10 | OPEN (MC, repr) |
| (42, 1.041); 0.960; 0.929 | +− | 7.85 (−6.06, +1.79) ± [3.73 / 7.90 / 6.47] | 9.45 (5.32 + 4.13) + 42.0 | 1.59 | 77.4 | OPEN |
| | ++ | 0.38 ± [0.62 / 0.87 / 0.13] | 0.001 + 0.09 | 1.59 | 2.10 | OPEN (MC, repr) |
| (46, 1.045); 0.957; 0.918 | +− | 32.2 (+0.48, −31.7) ± [37.2 / 120 / 4.5] | 30.5 (17.0 + 13.5) | 2.05 | 225 | OPEN |
| | ++ | **4.01** (−3.86, −0.15) ± [5.61 / 7.22 / 0.01] | 0.007 | 2.05 | 16.9 | OPEN; **central value 2.0 × tol** ‡ |
| (50, 1.049); 0.953; 0.915 | +− | 16.8 (+15.5, −1.30) ± [23.8 / 51.6 / 4.5] | 44.3 (25.2 + 19.1) + 182 | 2.57 | 323 | OPEN |
| | ++ | **6.60** (+6.57, +0.03) ± [6.91 / 2.39 / 0.38] | 0.013 + 0.81 | 2.57 | 17.1 (lower bound 1.30) | OPEN; **central value 2.6 × tol** ‡ |

All numbers are in S_E,t.
* *mc*: 3σ jackknife over the 8 RQMC replicas.
* *repr*: 3 × the bulk shape-degree spread (K = 6 against 4, 5, 7).
* *quad*: 3\|K15 − G7\| of the measured edge sum.
* *edge*: sign-free, with \|Γ\| ≤ 3A((1 − x_l)/(1 − x))² on the n3-only GK nodes. These nodes were never run.
* *gap*: measured-range GK nodes that have no usable value, covered by 3 × the largest \|Γ\| of their four nearest valid nodes. A node has no usable value if it lies in the reference-zero window, has a helicity off-diagonal > 2 % after the repair (§1.1), or has a relative error > 0.2.
* *Upper bound*: Σ_k(\|D_k\| + budget + remainder).
* *Lower bound*: Σ_k max(0, \|D_k\| − budget − remainder). It is 0 everywhere except (50) ++, where it is 1.30 < 2.57.
* † J = 38: the nodes at x = 0.906–0.932 fail the 2 % off-diagonal test (n1 MC failure), so 10 GK nodes are unmeasured.

**‡ Prominent warning, not an established violation.** The ++ central values at (46, 1.045Λ) and (50, 1.049Λ) exceed the tolerance, by 1.95× and 2.56×. Their 3σ budgets reach below the tolerance, so this is not a violation. But these values come from a bulk that is aliased: at J ≥ 46 the 0.025 grid is about one oscillation wavelength of D_h at x ≈ 0.6. The grade is [NE] at best. Before any claim about Theorem W at J ≈ 46–50 these must be rechecked on a denser grid (§4).

**Per-state verdict of (EDGE-R):**
* **Open at all 8 states.** The binding entry is +− everywhere.
* **Narrowed at (26, 1.026) and (30, 1.030):** the ++ entry is closed there [NE + M]; only +− remains.
* Elsewhere, both entries are open.

## 1. Method

### 1.1 A finding the analysis depends on: the external even-J records are in a different structure basis [DV]
* **What goes wrong without a fix.** Injected as is (`hc_nodes.to_record` with the local shadow b), the batch-3 even-J records are **not diagonal in the chain's helicity basis**: off/diag is ≈ 0.1–2. For local direct records it is ≲ 5·10⁻⁵ (checked on HIGHSPIN's J = 16 row). The resulting Γ_+− is wrong by ×10–10³. At J = 26 it gives rms 8·10³ c₂ on x < 0.6, where the true value is ≈ 10².
* **Cause.** `l5_spin.conserved` returns an SVD basis of the ⟨TTO⟩ null space. `hc_points`' own docstring requires the points to be in the same basis as the local rows. The external machine's basis differs by a general 2×2 transformation, not only by column signs (all 16 sign patterns were tried). Odd J is 1×1, so a basis change is a constant there and cancels. That is why EXT_FOLD's odd rows, and J = 33 here, are unaffected.
* **Repair (`e2_basis.py`).** In any basis E(ν) = U C(ν) Vᵀ, with C helicity-diagonal and U, V constant.
  * U and V are fitted from the external EE records themselves by joint diagonalisation over all nodes, separately for the state and for Δ_ref = 2·10⁴. The reference nodes straddle D_++'s zero, which separates the two helicities.
  * The records are mapped into the local basis, E′ = T U⁻¹E V⁻ᵀ(bT)ᵀ, and the chain then runs unchanged.
  * The per-helicity constants are absorbed in a_h, which is common to k = 2 and k = 3 because the same U, V are applied to C3H.
  * The helicity labels are the permutation that minimises rms Γ. The alternatives are ×10²–10³ worse.
* **Checks.**
  * The fitted EE residual is off/‖C‖ ≤ 4·10⁻⁷ (median, J = 26) to 3·10⁻⁴ (J = 50). For the reference it is ≤ 2·10⁻⁹.
  * **C3H was not fitted and is diagonal too**: median 5·10⁻⁶ (J = 26) to 3·10⁻⁴. This confirms the structure-basis explanation.
  * Γ_++ and Γ_+− come out at the size of HIGHSPIN's J = 16 and J = 29 rows, and with their sign pattern: C3H ≈ C2 for ++, C3H ≈ −C2 for +−.
* **Consequence for other jobs.** Every even-J external job (P3 included) needs this repair. Alternatively, the external runs should ship `conserved(Δ, J)` / `conserved(3 − Δ, J)`, or the shadows computed there.

### 1.2 Γ and the costs (`e2_k2.py`, `e2_gamma.py`)
* **The chain.** Γ_{k,h}(x) (c₂, κ subtracted) comes from the project's chain, unchanged: hs_core → ng_val injection → w2/ee1000 EE mode → l9 grid, with a_h, the de-shift and the Jacobi fix.
* **Inputs it needed.** No K2 rows exist at J = 26–50, so the chain's l8 helicity basis needed K2 EE records. These are local, cost ≈ 30 s each, and are not row nodes:
  * the state row at x = 0.2–0.8;
  * the reference at x_ref = 0.3–0.8 for J ≥ 34 (J = 26, 30 have w1 refb);
  * the shadows b(Δ, J).
* **Every node is a direct n1 node.**
* **Costs.** As in EDGE_R §1 (er_core weights, the W state's own cutoff):
  * bulk [x_w, a): degree-6 shape fit of the 0.025 grid;
  * edge: the job list's two GK15 windows with K15 weights on the measured nodes.
* **Two code issues met and worked around** (read-only modules, not edited):
  * (i) `er_jobs.state_nodes` reports the edge start one shift late when its loop ends without a break. At J = 26 the nodes lie on [0.818, …], not [0.824, …]. `e2_gamma.layout` reconstructs the real windows.
  * (ii) `ng_val._clean` recognises injected nodes to \|Δν\| < 10⁻⁵, but hs_core's re-keying moves ν by up to ≈ 5·10⁻⁴. Non-round GK nodes could therefore be Hampel-despiked; J = 26, x = 0.8961 was dropped this way. `e2_gamma` uses a copy with 10⁻³. HIGHSPIN's J = 12, 16 rows (round x) are not affected; this was checked.
* **The reach of the n1 rows.** The rows reach ν/Δ < 0.93, i.e. x_l = 0.915–0.931, not x ≤ 0.945. The rest of [x_l, x₁] is the job list's n3-only nodes: 5–8 per state, none of which have been run.

### 1.3 The remainder: envelope calibrated on the exact scalar (`e2_calib.py`) [DV]
The exact scalar (Arb) was run at the same 7 even J, with the graviton's own windows and its own x_l.
* **The pointwise factor the envelope A((1 − x_l)/(1 − x))² needs is c_pt = 1.21–2.06.** It is largest at J = 50. The c = 3 used here covers the scalar with margin ≥ 1.46.
* The scalar's \|Γ\| grows like (1 − x)^−q, with q = 2.4–3.0 on [0.8, x_l] and 2.3–2.8 beyond. The graviton's measured q is 1.9–3.4. A variant with the measured q is in `cache/report.txt`; it raises the +− remainder by ≤ 25 %.
* **The part beyond x_l is not a correction.** For the scalar, E_rem ≈ −E_det at J ≥ 38. For example, at J = 46: E_det = +4.25, E_rem = −4.13, E_full = +0.125 S_E2. The signed remainder is 0.01–1.5 × the sign-free c = 1 envelope. So x > x_l must be computed with its phase. It can be neither neglected nor extrapolated (as in EDGE_R §4).
* GK15 on the two windows reproduces the dense integral to < 10⁻⁴ S_E2. The quad term in the table is therefore an artefact of zeroing the missing nodes, and it disappears once the windows are complete.

## 2. Grades
* **Determined parts.** [NE±] at J = 26, 30. [NE] at J ≥ 33, where the 0.025 bulk grid aliases near x ≈ 0.6–0.8 and the repr spread may understate the representation error.
* **Remainder bounds.** [M]: a growth envelope, calibrated on the scalar only.
* **++ at J = 26, 30 is closed [NE + M]; the M is very weak.** The ++ weights vanish at the edge (env₊₊ = 2 − y; ψ_k ~ (1 − p)⁵). The unmeasured nodes would break the tolerance only if \|Γ_++\| exceeded ≈ 6.7·10⁷ c₂ there. That is 4300× (J = 26) and 3600× (J = 30) the last measured amplitude, i.e. a growth exponent q ≳ 8 against the measured 2.1–2.8.
* The off-diagonal entry is not evaluated (as in EDGE_R §1; its sign is basis-convention dependent).

## 3. Remarks
* **The C2/C3H cancellation is lost by the premise.** For +−, Γ_C3H ≈ −Γ_C2, and the chain's actual correction is cost₂ + cost₃. The premise as stated charges \|cost₂\| + \|cost₃\|. The combined D₂ + D₃ is smaller at several states: +0.05 against 0.27 at J = 26, and −0.19 against 1.83 at J = 34. A combined-weight version of Lemma P4W-χ (one premise for ψ₂G₂ + uψ₃G₃) would target that number. This is a suggestion, not checked [M].
* **The ++ / +− asymmetry is structural.** The edge x > 0.93 matters only for +−: env₊₋ = (2 − y)/(1 − y)⁴. For ++, the remaining work is bulk-MC work.

## 4. What would close each state

| state | +− entry needs | ++ entry needs |
|---|---|---|
| (26, 1.026) | (a) the n3 k = 2 edge nodes (the `P1_J26_*_n3a/b` jobs, 10 nodes × 2 seeds); (b) **the n3 C3H fix**: the k = 3 edge alone is 0.99 > 0.406; (c) MC ×2–3 lower (M = 13 on x ≥ 0.6) and a denser bulk on [0.5, a) (repr 0.41) | done |
| (30, 1.030) | (a) + (b) (k = 3 edge 1.62 > 0.63); (c) MC ×4 (2.10 → ≲ 0.3): n3 or M + 4 on x ≥ 0.8, where n1 has heavy tails | done |
| (34), (33), (38), (42) | (a) + (b) (the k = 3 edge alone is 3.5 / 1.9 / 16 / 4.1); (c) MC ×5–10 and the off-diagonal-failed gap nodes rerun; n1 fails the 2 % test at 2–6 nodes in [0.87, 0.93], so n3 is needed there too; (d) a 0.0125 bulk grid on [0.5, a) | MC ×2–3 and a denser bulk; no edge data needed |
| (46), (50) | (a)–(d); MC ×20–50: not reachable by raising n1's M, so n3 on x ≥ 0.8 is needed for both kinds, again requiring the C3H fix | **the decisive check (‡):** a 0.0125 bulk grid on [0.3, a) and M + 2. With central values at 2–2.6 × tol, this could turn into a violation of the ++ premise at J = 46–50. |

* **Every state's +− entry needs a working n3 C3H evaluator**, because the k = 3 edge remainder alone exceeds the tolerance everywhere. The n3 k = 2 jobs alone cannot close any state.
* **Any further even-J external output must go through the §1.1 repair**, or ship its structure basis.
* **Measured costs.** n1: 2.4 min per node at J = 26, M = 11; ×4 per +2 in M. n3: 8–10 min per node.
* **Size of the ++ recheck at J = 46, 50.** About 20 extra n1 nodes per state, ≈ 1–1.5 CPU-h each at M = 11. These can run on the external machine, with the repair applied afterwards.

## 5. Files
| file | content |
|---|---|
| `edger2/e2_merge.py` | JSON → `cache/pts/` point files; node inventory `cache/merge.txt` |
| `edger2/e2_k2.py`, `run_k2.sh` | K2 basis records + shadows → `cache/k2/` |
| `edger2/e2_basis.py` | structure-basis repair (§1.1) → `cache/ptsfix/` |
| `edger2/e2_gamma.py` | Γ, the costs and their budgets per state → `cache/gamma_J*.{txt,json}` |
| `edger2/e2_calib.py` | exact-scalar calibration → `cache/calib.{txt,json}` |
| `edger2/e2_report.py` | the table → `cache/report.{txt,json}` |
