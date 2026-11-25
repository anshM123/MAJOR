# BATCH5_FOLD: the batch-5 external direct rows folded into the Δ_g = 1000 near-gap budgets (J = 17–32)

*25 Nov 2026, batch-5 fold agent. Unreviewed, not committed by this agent. `git status` shows `code/b5fold/b5_core.py` and `b5_state.py`
as modified: another agent's snapshot commit (978b1cd, "WIP: batch-5 fold and c1 agents (snapshots)") swept in intermediate versions, and the working
tree holds the final ones. No existing file was edited; `ads4/prl/`, `code/edger2/` and `code/ngfill/` were only read.*
* *New code:* `ads4/code/b5fold/` (guard `b5_guard.py`; every write outside `b5fold/cache/` raises). *Caches:* `b5fold/cache/`, with a
  `.gitignore` that ignores everything but itself.
* *Compute:* ≈ 1.8 CPU-h. At most two computing processes at a time, `ulimit -v 3000000`, one thread, no bytecode written.

Grades follow NEARGAP2 / HIGHSPIN / SMALLNU / EXT_FOLD / NEARGAP_FILL: **[NE±]**, [NE + M], [M], [I]. `L = value / budget`, with ×3 on
model terms and 3σ on Monte Carlo. Suffixes:
* *(λ-level)*: closes only with NEARGAP2's λ-level truncation, as J = 8 does.
* *(signed window, thin)*: closes only with SMALLNU's signed window, as J = 29 does.
* *held*: the mass lies above every row of its J, so the 1/m model holds the own row. The hs_odd protocol does this at every odd-J
  mass, and NEARGAP_FILL §2.3 at J = 11, 13, 15.

---

## 0. Result

**Per-state table** (Δ_g = 1000):
* "old" is the published or model value, or for J = 23 and 31 the EXT_FOLD / NEARGAP_FILL value. "new" is the full λ/N_s; λ_c in
  brackets is the value the leverage uses, with the windows ν < 50 and x > x_last removed and bounded separately.
* Budget: sign-free, ×3 measured window envelope, representation/interpolation term evaluated on λ_c (NEARGAP_FILL §3).
* "alt": λ-level for the own-row route; signed window for the hs_odd route.
* "tw" marks the twist edge (Δ = 1000 + J); "(h)" marks a held state.

| J | m/Λ | route | λ/N_s: old → **new** (λ_c) ± jk | budget | **L** (alt) | grade old → new |
|---|---|---|---|---|---|---|
| 17 | 1.0169 tw | own row | +1.003 → **+1.001** (+1.005) ± 0.002 | 0.656 | **1.53** (λ-lvl 1.53) | [NE + M] → **[NE±]** |
| 17 | 1.02 / 1.03 / 1.05 (h) | own row | +1.002 / +0.999 / +0.841 → **+1.002 / +1.004 / +0.848** | 0.49 / 0.18 / 0.05 | **2.04 / 5.56 / 18.4** | [NE + M] → **[NE±], held** |
| 18 | 1.0179 tw | own row | +0.740 → **+0.746** (+0.744) ± 0.002 | 0.407 | **1.83** (λ-lvl 1.83) | [NE + M] → **[NE±]** |
| 18 | 1.02 / 1.03 / 1.05 | own row | +0.743 / +0.909 / +0.994 → **+0.749 / +0.914 / +0.999** | 0.41 / 0.50 / 0.82 | **1.81 / 1.83 / 1.21** (λ-lvl 2.0 / 4.0 / 5.6) | [NE + M] → **[NE±]** |
| 19 | 1.0188 tw | own row | +0.995 → **+1.006** (+1.007) ± 0.001 | 0.547 | **1.84** (λ-lvl 1.84) | [NE + M] → **[NE±]** |
| 19 | 1.02 / 1.03 / 1.05 (h) | own row | +0.995 / +0.926 / +0.918 → **+1.006 / +0.937 / +0.928** | 0.51 / 0.20 / 0.09 | **1.97 / 4.64 / 10.1** | [NE + M] → **[NE±], held** |
| 20 | 1.0198 tw | own row | **+0.965 → +0.835** (+0.835) ± 0.002 | 0.594 | **1.41** (λ-lvl 1.41) | [NE + M] → **[NE±]** |
| 20 | 1.02 / 1.03 | own row | +0.969 / +0.990 → **+0.839 / +1.004** | 0.60 / 0.86 | **1.40 / 1.17** (λ-lvl 1.42 / 2.06) | [NE + M] → **[NE±]** |
| 20 | 1.05 | own row | +0.991 → **+0.975** | 1.157 | 0.84 (**λ-lvl 2.89**) | [NE + M] → **[NE±] (λ-level)** |
| **22** | **1.0218 tw** | own row | +0.901 → **+0.996** (+1.003) ± 0.008 | 1.365 | **0.73** (λ-lvl 0.73) | [NE + M] → **[NE + M]** |
| 22 | 1.03 / 1.04 / 1.05 (h) | own row | +0.904 / +0.880 / +0.845 → **+1.002 / +0.974 / +0.932** | 0.68 / 0.36 / 0.23 | **1.48 / 2.74 / 4.07** | [NE + M] → **[NE±], held** |
| **23** | **1.0228 tw** | hs_odd/GK | +0.923 → **+0.926** (+0.904) ± 0.012 | 0.557 | **1.62** (signed 1.79; unchanged protocol **1.45**) | [NE±] thin (1.18) → **[NE±]** |
| 23 | 1.03 / 1.05 / 1.075 / 1.10 | hs_odd/GK | +0.945 / +1.021 / +0.970 / +1.018 → **+0.948 / +1.023 / +0.973 / +1.021** | 0.32 / 0.19 / 0.16 / 0.16 | **2.92 / 5.18 / 5.85 / 6.37** | [NE±] → **[NE±]** (strengthened) |
| **24** | **1.0237 tw** | own row | +0.945 → **+1.024** (+1.025) ± 0.010 | 1.341 | **0.76** (λ-lvl 0.76) | [NE + M] → **[NE + M]** |
| 24 | 1.03 / 1.04 / 1.05 | own row | +0.954 / +0.902 / +0.937 → **+1.021 / +0.953 / +0.982** | 1.28 / 1.33 / 1.65 | 0.80 / 0.71 / 0.59 (**λ-lvl 1.13 / 1.78 / 2.06**) | [NE + M] → **[NE±] (λ-level)** |
| 25 | 1.0247 tw | hs_odd/GK | +0.903 → **+0.960** (+0.930) ± 0.018 | 0.564 | **1.65** (signed 1.91) | [NE + M] → **[NE±]** |
| 25 | 1.03 / 1.04 / 1.05 | hs_odd/GK | +0.966 / +0.972 / +0.951 → **+1.028 / +1.034 / +1.011** | 0.41 / 0.31 / 0.28 | **2.44 / 3.26 / 3.53** | [NE + M] → **[NE±]** |
| **26** | **1.0257 tw** | own row | — (no row) → **+0.985** (+0.978) ± 0.012 | 1.031 | **0.95** (λ-lvl 0.95; 'measured' window 1.07) | [M] → **[NE + M]** |
| 26 | 1.03 / 1.04 / 1.05 (h) | own row | — → **+0.956 / +0.956 / +0.986** | 0.77 / 0.45 / 0.40 | **1.24 / 2.12 / 2.44** | [M] → **[NE±], held** |
| **28** | **1.0276 tw** / 1.03 | own row | +0.990 / +0.990 → **+0.950 / +0.948** (+0.935) ± 0.014 | 1.42 / 1.44 | **0.66 / 0.65** (λ-lvl 0.66 / 0.79) | [NE + M] → **[NE + M]** |
| 28 | 1.04 / 1.05 | own row | +1.019 / +1.021 → **+0.995 / +1.000** | 2.02 / 2.67 | 0.49 / 0.37 (**λ-lvl 1.18 / 1.21**) | [NE + M] → **[NE±] (λ-level)** |
| **30** | **1.0296 tw** / 1.03 (h) / 1.04 (h) | own row | — → **+0.900 / +0.903 / +0.938** (λ_c +0.875 / +0.878 / +0.921) | 1.87 / 1.85 / 1.02 | **0.47 / 0.48 / 0.90** | [M] → **[NE + M]** |
| 30 | 1.05 (h) | own row | — → **+0.941** (+0.927) | 0.779 | **1.19** | [M] → **[NE±], held** |
| **31** | **1.0306 tw** | hs_odd/GK | +0.862 → **+0.920** (+0.896) ± 0.032 | 0.999 | **0.90** (**signed 1.14**) | [NE + M] (0.76) → **[NE±] (signed window, thin)** |
| 31 | 1.04 / 1.05 / 1.075 / 1.10 | hs_odd/GK | +0.869 / +0.907 / +0.926 / +0.964 → **+0.920 / +0.954 / +0.963 / +0.995** | 0.66 / 0.57 / 0.50 / 0.46 | **1.35 / 1.64 / 1.88 / 2.13** | [NE±] → **[NE±]** (strengthened) |
| **32** | **1.0315 tw** | own row | +0.913 → **+0.955** (+0.941) ± 0.017 | 2.321 | **0.41** (λ-lvl 0.41) | [NE + M] → **[NE + M]** |
| 32 | 1.04 / 1.05 / 1.075 | K₂ rows 1040 / 1600 only (m > m₁₀₄₀ = 1.0396Λ; the own row is not used) | +0.947 / +0.954 / +0.999 (unchanged) | trunc 60–90 | ≈ 0 | [NE + M] unchanged |

**Negative findings and cautions (stated first).**
1. **No state came out negative.** The lowest values are J = 18 at its twist edge (+0.746), J = 20 (+0.835), and J = 30 (+0.900,
   λ_c +0.875). Every [NE + M] row above has a positive value; what fails there is the *budget* (lower bound < 0), not the sign.
2. **Even-J external records were in a different ⟨TTO⟩ structure basis** (as EDGE_R2 found on batch 3). Injected as they were, J = 18
   gave λ = −0.12 and λ_c = −14.3, with Γ₊₋ and c_off at 10⁴–10⁶ c₂. That run was a basis artefact, not physics: it was discarded, and
   every even number above uses the repaired records (§1.2). Odd J (1×1) is unaffected.
3. **J = 20's value drops by 0.13**: published +0.965, direct +0.835 at the twist edge. This is unfavourable, but it confirms HIGHSPIN's
   MC-limited edge estimate (+0.83 ± 0.14), and the state still closes (L = 1.41). J = 28 drops by 0.04 (+0.990 → +0.950). All the
   other shifts are favourable (+0.004 … +0.10).
4. **hs_odd's shape fit is ill-conditioned at J = 17 and 19.** Unchanged, the protocol gives λ = +0.50 ± 0.26 (J = 17) and
   L ≈ 0.0 / 0.13. The degree-6 shape basis has condition number 6.0·10⁷ (J = 17) and 1.7·10⁷ (J = 19), against ≤ 4.6·10⁵ wherever
   it was used before (J = 21–31). Its continuation over [0.775, 0.8) blows up to 3·10⁴ c₂. Degrees 4, 5 and both half-samples give
   +1.00. These two states are therefore folded with the own-row protocol that NEARGAP_FILL used for odd J ≤ 15 (§2.1). Both routes give
   the same value (λ_c, shape degree 5 vs own row: J = 17 +1.002 vs +1.005; J = 19 +0.998 vs +1.007).
5. **ng_bind's Hampel despike hit the direct own rows.** It removed real edge nodes, up to 5 per row (J = 26: 0.8226, 0.8362, 0.923,
   0.925, 0.9387); NEARGAP_FILL's J = 14 lost 0.925 this way. The HC protocol exempts injected direct rows (ng_val's `_clean`), and that
   exemption is applied here (§2.2). It changes J = 26's twist-edge L from 0.49 to 0.95 and leaves the values unchanged (Δλ_c ≤ 0.001).
6. **Even twist edges J ≥ 22 do not close.** The binding terms are the sign-free edge beyond the last valid node (G = 10 × max|Γ| of the
   last 4 nodes, flat, x > 0.926–0.939: 0.6–1.4 N_s) and, from J = 26, the ×3 window envelope (0.38–0.79). Under the odd-J
   light-cone envelope the edge term would halve, giving L ≈ 1.2 / 1.2 / 1.4 / 1.0 for J = 22 / 24 / 26 / 28. That is a protocol
   change, used for no grade (§3.3).
7. **Held states.** J = 17, 19 (all masses above the twist edge), J = 22 (1.03–1.05: its only other row is 336) and J = 26, 30 (no other
   row) are held at the own row. No neighbouring row exists for a NEARGAP_FILL §2.3-style hold estimate, so the hold is a stated model
   input (H_mono-A). They are graded "[NE±], held"; the strict count below excludes them.

**Twist edges of J = 23 and 31 (the decision asked for).**
* **J = 23: [NE±], no longer thin.** With the 25 bulk mid-points the K7 representation probe is resolved: E_repr on λ goes 0.553 → 0.087
  and on λ_c 0.258 → 0.020, and the binding probe is now half0 +0.007 ± 0.004.
  * Twist edge **L = 1.62** (×3, λ_c); signed 1.79; even the *unchanged* protocol gives 1.45 (was 0.85).
  * The second GK window [0.9, Λ/m] as the edge quadrature (HS_XCUT = 0.9) gives λ_c +0.908 vs +0.904 and L = 2.31: an independent
    check of the edge.
  * Binding term: the missing GK nodes x > 0.945 and the failed x = 0.941 (E_miss 0.349).
* **J = 31: twist edge [NE±] only on the signed window (L = 1.14), like J = 29. The ×3 envelope does not close (L = 0.90).**
  * The half-sample probe that NEARGAP_FILL found to be real under-resolution shrinks with the mid-points: λ_c probe +0.093 → +0.019.
  * Remaining: window ×3 0.427 (ν < 20 at 10⁴ c₂: 0.189, unreachable by n1), GK − G7 quadrature 0.213, missing GK nodes 0.207.
  * The second GK window confirms it: L ×3 0.93, signed 1.19.
  * 1.04–1.10Λ strengthened (1.35–2.13).

**Count: error-controlled near-gap spins J ≤ 32 at Δ_g = 1000** (twist edge and every computed mass to 1.05Λ at [NE±], any suffix):
* Before batch 5: J = 0, 2, 4–16, 21, 23, 27, 29 (**19**).
* **Now: J = 0, 2, 4–21, 23, 25, 27, 29, 31 = 25 of the 31 spins that have a near-gap state** (J = 1, 3 have none). Added: J = 17, 18,
  19, 20, 25, and J = 31 (signed window, thin).
* **Still [NE + M] at the twist edge:** J = 22, 24, 26, 28, 30, 32. Above it they are [NE±] held / λ-level, except J = 28 at 1.03Λ,
  J = 30 at 1.03–1.04Λ and J = 32 above the twist edge (K₂ rows only).
* Stricter readings:
  * no signed-window-only spins (drop 29, 31): 23;
  * also no held masses (drop 17, 19, and 11, 13, 15 from NEARGAP_FILL; the hs_odd masses above the twist edge, which hold the own
    row by construction, are kept as in EXT_FOLD): 18;
  * sign-free ×3 at every mass, also dropping λ-level (8, 20): 16.

---

## 1. Inputs and the structure-basis repair

### 1.1 Points (`b5_merge.py` → `cache/pts/{all,old,base}`, log `cache/b5_merge.txt`)
* All 666 batch-5 jobs are present (`ext/unpacked/out/N5_*`). They were merged with the existing direct points: spinfunc / neargap2 /
  highspin / smallnu caches, and for J = 23 and 31 extfold's `all` set.
  * 208 key collisions, all among old sources and bit-identical.
  * No batch-5 key collides.
* Failed by the replica-error rule (> 20 %, `hs_core.failed_x`) in the external basis: x = 0.95 on every even row, and
  * GK nodes 0.9439 (J = 19), 0.9417 (J = 22), 0.9402 (J = 24), 0.9425 (J = 31 edge2);
  * the rz flank 0.9453 (J = 25); 0.9472 / 0.9722 (J = 24);
  * a few interior nodes (J = 20, 24 x = 0.725; J = 28 0.75 / 0.925; J = 30 0.825; J = 32 0.80–0.875, 0.8847, 0.9191).

  The element-wise rule is basis-dependent, so on the mapped records (§1.2) the failed set is recomputed. Interior failures mostly sit
  where one matrix element crosses zero: Γ and σ are normal (e.g. J = 22 x = 0.875: σ ≈ 120 c₂). They are dropped as the protocol says.

### 1.2 Even J: external records mapped into the local ⟨TTO⟩ structure basis (`b5_points.py`, `b5_basis.py`, `b5_basischk.py`)
**Mechanism.**
* hc_points' even record is E = Σ F[a] G̃[c], with F linear in `sol = l5_spin.conserved(Δ, J)` and G̃ linear in `conserved(3 − Δ, J)`.
* `conserved` returns U[:, :2]/colmax, where U is an orthonormal SVD basis of a 2-dim null space. U is fixed only up to Q ∈ O(2)
  (degenerate null singular values; LAPACK-dependent).
* Hence, exactly, **E_ext = Gᵀ E_loc H** with G = diag(c_loc) Q diag(1/c_ext), and the same form for H.
* G and H each depend on one angle and one reflection bit. They are separate for the state Δ and for Δ_ref = 2·10⁴, and shared by the
  EE and C3H kinds.

**Determination: anchored, not fitted to the physics.**
* One batch-5 node per even row (x = 0.3, EE, state and reference) was recomputed **locally** with ext/point.py's code path (n1 is
  deterministic up to the adaptive proposal), ≈ 1 min per node.
* The two angles were fitted on the 8 replicas plus the mean. Every record was then mapped, E_loc = G⁻ᵀ E_ext H⁻¹ (replicas one by one;
  `err` recomputed from them as hc_points does).
* **Checks:**
  * a fully general 2×2 G, H fit gives the same residual (so the O(2) form is right);
  * a second, *unfitted* local node (x = 0.6, EE and C3H) at J = 18 and 26 is predicted;
  * at the chain level, Γ of the injected row with the local record substituted is unchanged.

| J | raw ‖E_loc − E_ext‖/‖E‖ (state / ref) | fit residual (state / ref) | 2nd node, unfitted (state / ref) | Γ(mapped) − Γ(local) at the local nodes (c₂; σ) | EDGE_R2 joint diag. (`e2_basis`, read-only): off/‖C‖ median, EE state / ref; Γ vs mapping |
|---|---|---|---|---|---|
| 18 | 1.5 / 1.1 | 9·10⁻⁸ / 8·10⁻⁹ | 2·10⁻⁶ / 2·10⁻⁸ (EE + C3H) | ≤ 0.02 (σ 19–32) | 2·10⁻⁷ / 2·10⁻¹¹; row did not build (sign of a₊₊) |
| 20 | 0.49 / 0.003 | 5·10⁻⁶ / 2·10⁻¹² | — | ≤ 0.15 (σ 35–38) | 2·10⁻⁷ / 3·10⁻¹¹; row did not build |
| 22 | 0.12 / 0.04 | 2·10⁻⁷ / 3·10⁻¹⁰ | — | ≤ 0.01 | 2·10⁻⁷ / 6·10⁻¹¹; diag. Γ agree ≤ 6·10⁻⁴ c₂ (52 nodes) |
| 24 | 1.4 / 1.1 | 1·10⁻⁵ / 2·10⁻¹⁰ | — | ≤ 0.004 | 6·10⁻⁷ / 1·10⁻¹⁰; row did not build |
| 26 | 0.15 / 1.7 | 3·10⁻⁷ / 1·10⁻⁹ | 7·10⁻⁵ / 5·10⁻⁹ (EE + C3H) | ≤ 0.11 (σ 32–53) | 4·10⁻⁷ / 1·10⁻¹⁰; diag. Γ agree ≤ 3·10⁻³ c₂ (53 nodes) |
| 28 | 0.34 / 0.05 | 1·10⁻⁶ / 5·10⁻¹¹ | — | ≤ 0.01 | 3·10⁻⁶ / 1·10⁻⁹; row did not build |
| 30 | 0.02 / 1.05 | 2·10⁻⁵ / 4·10⁻⁹ | — | ≤ 0.34 (σ 46–50) | 1·10⁻⁶ / 3·10⁻⁹; diag. Γ agree ≤ 0.03 c₂ (54 nodes) |
| 32 | 1.1 / 0.21 | 6·10⁻⁷ / 1·10⁻¹⁰ | — | ≤ 0.02 | 6·10⁻⁶ / 5·10⁻¹⁰; diag. Γ agree ≤ 6.3 c₂ (54 nodes) |

* After the mapping the rows are helicity-diagonal: |c_off| ≤ 10 c₂ below x = 0.92, against 10⁵–10⁶ unmapped.
* **EDGE_R2's joint diagonalisation agrees** wherever its output builds a row (J = 22, 26, 30, 32). There the diagonal Γ is identical
  to ≤ 6 c₂ (median 10⁻⁴). The off-diagonal differs by ≤ 260 c₂ near failed edge nodes, because e2_basis forces the EE off-diagonal to
  ~0.
* For J = 18, 20, 24, 28 its arbitrary per-helicity sign leaves a negative row normalisation (a₊₊ < 0), and only one node survives in
  this chain. So no comparison was possible there. The anchored local mapping is used throughout.
* Other even-row inputs (no K₂ row exists at Δ = 1000 + J for even J):
  * the even shadows b(1000 + J, J), via `b5_shadow.py`, the nf_rows function. The same call reproduces (1040, 32) bit for bit, and
    agrees bit for bit with EDGE_R2's (1026, 26) and (1030, 30);
  * 5 K₂ EE state records per row (x = 0.2–0.8, `b5_k2.py`) for the chain's l8 helicity basis. They are not row nodes.

## 2. Protocols

### 2.1 Routes
* **Even J** (`b5_state.py` = ngfill/nf_state, `neargap2/ng_bind.run_state` unchanged; `b5_window.py` = `smallnu/sn_even.main`
  unchanged, FULL normalisation):
  * the direct own row, x = 0.02 … 0.95, plus the GK nodes as extra row nodes;
  * budget: trunc (0 on the own row; the 1600 K₂ row above it), interpolation, the edge beyond the last valid node, 3σ MC;
  * the ν < 50 window replaced by the measured ×3 envelope, with 10⁴ c₂ only below ν = 20.
* **Odd J = 23, 25, 31** (`b5_odd.py`: hs_odd via an unchanged copy `b5_hsodd.py` + `sn_odd.main` + the λ_c representation split of
  nf_reprcut): exactly NEARGAP_FILL §3 / EXT_FOLD.
* **Odd J = 17, 19:** the own-row route (as NEARGAP_FILL for J = 9–15 and EXT_FOLD for J = 7). The window uses sn_odd's envelope and
  signed formulas (copied) on ng_bind's λ_c. hs_odd numbers are kept as diagnostics (`cache/odd_main`, `odd_kb5`).
* Cross-check, J = 25 by the own-row route: λ_c +0.931 vs +0.930. Its flat even-style edge term (1.4) would not close it (L = 0.58);
  the designated hs_odd route gives 1.65.

### 2.2 Choices made here (all stated before grading; none changes a sign)
1. **Representation/interpolation on λ_c** (NEARGAP_FILL §3), applied to both routes. For even J the effect is ≤ 0.02 in budget.
2. **Hampel exemption for injected direct rows** (§0 item 5). The as-run (despiked) budgets are kept in `cache/state_asrun/`. Twist-edge
   λ_c moves ≤ 0.001.
3. **Reference-zero rule (J ≥ 24 even):** direct nodes strictly closer to x₀(++) than the batch's flanking rz nodes (|x − x₀| < 0.012)
   are dropped. This is fixed by the job design, not by Γ: J = 24 0.95; J = 26 0.875, 0.8875; J = 28 0.8117, 0.8224, 0.825;
   J = 30 0.775; J = 32 0.725. For odd J = 25, hs_odd's ZWIN = 0.002 applies (no GK node within).

### 2.3 What blocks the even twist edges (`b5_edgeenv.py`, diagnostic)

| J | 22 | 24 | 26 | 28 | 30 | 32 |
|---|---|---|---|---|---|---|
| edge term (flat 10× G, x > x_last) | 1.12 | 0.93 | 0.60 | 0.81 | 0.99 | 1.38 |
| … under the odd-J envelope 3A((1 − x_l)/(1 − x))² | 0.57 | 0.45 | 0.25 | 0.32 | 0.40 | 0.39 |
| window ×3 (10⁴ convention) | 0.21 (0.83) | 0.29 (1.07) | 0.38 (1.33) | 0.50 (1.60) | 0.64 (1.97) | 0.79 (2.40) |
| L: protocol → odd-J edge envelope (not used for grades) | 0.73 → 1.23 | 0.76 → 1.20 | 0.95 → 1.42 | 0.66 → 1.01 | 0.47 → 0.68 | 0.41 → 0.71 |

n1 cannot reach beyond x ≈ 0.94 (ν/Δ ≳ 0.94 fails), and below ν = 20 (crystal fold). J = 30 and 32 are also limited by the window and by
interpolation (0.21 at J = 30).

## 3. Reproduce
```bash
cd ads4/code/b5fold      # every command via lane.sh: ulimit -v 3000000; OMP_NUM_THREADS=1; PYTHONDONTWRITEBYTECODE=1
python3 b5_merge.py; python3 b5_shadow.py check; python3 b5_shadow.py 18,20,22,24,26,28,30,32
python3 b5_k2.py J 0.2,0.35,0.5,0.65,0.8                                   # even J
python3 b5_points.py J 1000+J <nu of x=0.3>   (cache/loc_nodes.txt); B5_KINDS=EE,C3H python3 b5_points.py 18 1018 609.798983 (and 26)
python3 b5_basis.py; python3 b5_basischk.py 18,...,32; python3 b5_e2chk.py 18,...,32
python3 b5_state.py J [masses]; python3 b5_window.py J                     # even J, and odd 17, 19 (25 cross-check)
python3 b5_odd.py J main [masses]                                          # 23, 25, 31 (17, 19 diagnostics); B5_SET=old ... old (repro)
HS_XCUT=0.9 B5_GKMIN=8 python3 b5_odd.py 23 x09 ...; HS_KB=5 python3 b5_odd.py 17 kb5 ...
python3 b5_edgeenv.py 18,...,32; python3 b5_sum.py                        # -> cache/b5_sum.txt (all numbers above)
```
**Reproduction check.** With set `old`, b5_odd reproduces NEARGAP_FILL exactly. J = 23 twist edge: L 0.85 → 1.18, signed 1.27, value
+0.923 / λ_c +0.901. J = 31: 0.70 → 0.76, signed 0.94.

## 4. Caveats
1. The direct n1 integrals are [NE] (the reference standard), not certified.
2. The even-J results rest on the structure-basis mapping. It is anchored on one local node per row and checked to ≤ 0.34 c₂ in Γ
   (σ ≥ 19). Only J = 18 and 26 have a second, unfitted node.
3. Model inputs kept:
   * 10⁴ c₂ below ν = 20;
   * the flat edge G beyond x_last (even, own-row route) or the light-cone envelope (hs_odd);
   * for held states, the own-row hold (no hold estimate was possible for J = 17, 19, 22, 26, 30).
4. The masses are a grid; between them H_cont as before. J = 32 above its twist edge is still on K₂ rows (1040 / 1600).
5. Choices introduced here: the Hampel exemption, the hs_odd → own-row switch for J = 17 and 19, and the rz-drop threshold. Each is
   argued in §2, none is tuned to Γ, and the as-run alternatives are in the cache.
