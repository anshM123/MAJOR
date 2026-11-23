# NEARGAP2: the near-gap rows of Result D on their own-mass rows, with error budgets; light J = 0, 2 single traces at finite Δ_L

*23 Nov 2026. Unreviewed, not committed by this agent (another agent's commits 508ea02 / 3a3d6c8 swept in intermediate snapshots of this
file and of `code/neargap2/`; the working tree holds the final versions). No existing file was edited; `PCYC_CLOSE.md` and `code/pcyc/`
were not touched.*
* *New code:* `ads4/code/neargap2/` — `ng_guard.py` (process-local write guard: every write outside `neargap2/cache/` raises),
  `ng_rows.py` (own-mass K₂ rows), `ng_points.py` (direct n1 RQMC points), `ng_chain.py` (the EE_NEARGAP chain, unchanged, with the new
  rows visible), `ng_bind.py` (per-state error budget), `ng_val.py` (direct vs K₂, node by node and at the state), `ng_engine5.py` +
  `ng_engcmp.py` (parity-odd analytic-radial engine at the binding state), `ng_rowdump.py`, `ng_light.py` (light states),
  `ng_edgew.py`, `ng_summary.py`, `ng_queue.sh`, `ng_pause.sh`.
* *Caches:* `ads4/code/neargap2/cache/` (its `.gitignore` ignores everything but itself).
* *Read-only imports:* `spinfunc/` (`w1_rows`, `w2_chain` via `ee1000_common`, `q3_analysis`, `p4_odd_analysis`, `l9_analysis`,
  `hc_points`, `hc_nodes`, `po_ee_odd`, `bg_eeq`), `prl/data/q3_full_delta0.2.json`; read-only data: the existing K₂ rows and shadows in
  `spinfunc/cache/`, `x2_chi_EE_can_*.pkl`, `u3_chi_EE_*.pkl`.
* *Compute:* ≈ 5.5 CPU-h, at most 2 processes of this note at any time (short analyses ran while one worker was SIGSTOPped, `ng_pause.sh`),
  `ulimit -v 3000000`, `OMP_NUM_THREADS=1` everywhere.

Grades as in UNIFIED_THEOREM_v2 §0 and POS_ERRORS: [P], [certified], [D], [D-sketch], [DV], [NE], [NE±] (numerical evidence with an
explicit error budget, every model input with a stated safety factor, conclusion surviving it), [NE + M], [M], [I].
`c₂ = m² ×` relative correction; `λ/N_s` as in UT_v2 §2.5.

---

> **Update (HIGHSPIN.md, unreviewed): J = 29 is resolved and positive; near-gap J = 12, 16 move to [NE±]; no physical state came out negative.**
> * **J = 29 twist edge (1.0286Λ):** λ/N_s = **+0.972 ± 0.018** [NE + M] (was [M]). Method: a full direct n1 row (34 nodes) plus GK15 quadrature of 16 direct edge nodes, with no interpolation. Budget 0.70, so L = 1.39. The one model input is the ν < 50 window: flipping the state needs |Γ| ≥ 3.2·10³ c₂ there with the wrong sign, while the measured values at ν = 51–77 are ≤ 116 c₂. Masses 1.04–1.10Λ: +0.956 … +0.995, L = 2.8 → 9.5.
> * **Diagnosis of the −3.62 [DV]:** the K₂ method fails, not the physics. The Δ_ref = 2·10⁴ reference's leading partial wave has its first zero at x₀ = 0.8008, right at the bad node (|r₂/r₀| = 0.86 vs ≈ 0.02 elsewhere). w2's Hampel filter skips odd J, so the bad node is kept. The +1.39 is spline aliasing. The published +0.948 removed the node, and the shape chain is within 0.025 of the direct result.
> * **The same reference-zero problem in published chains:** J = 31 (+1.208 → +1.02), J = 32 (+0.904 → +0.997), and the J = 40 and J = 48 outliers explained. All shifts are favourable.
> * **J = 12:** +0.8528, L = 3.6 / 3.9 / 2.5 (twist edge, 1.02Λ, 1.05Λ). **J = 16:** +0.9402, L = 1.9 / 1.7 / 1.4 / 1.2. Both now **[NE±]**, from full direct rows (38 nodes) that remove the truncation term. Odd J = 21, 23, 27 and even J = 20, 22, 28, 32 edges: [M] → [NE + M]. Not run: J = 26, 30, 48, 10, 14. **H_LJ(1000): unchanged [D-sketch + M].**
> * **Cautions:** the K₂ λ-level Laplace step cannot close a ×3 budget for J ≳ 17, although the actual K₂ error is only 4–5 % of the step and always favourable. The ν < 50 window is now the binding model input at large J. Lowest computed values: J = 21 at 1.03Λ, +0.868; J = 20 twist edge, +0.83 ± 0.14. Code bugs worked around: `p4_odd_analysis.grid` and `l9_analysis.grid` drop direct nodes whose x is not round to 6 decimals.

> **Update (SMALLNU_WINDOW.md, unreviewed): the ν < 50 window is now measured. The 10⁴ c₂ convention is kept only below the lowest measured node (6–18 % of the window weight). Nothing came out negative. J = 29 and J = 21 → [NE±] (thin); J = 12, 16 leverage roughly doubles; J = 23, 27, 31 stay [NE+M] (now limited by the K₂ bulk truncation, not the window).**
> * **Method:** direct n1 points at ν ≥ 20 and n3 points at ν = 10–17.5 on each state's own row, with jackknife errors. n1 fails below ν ≈ 20 (crystal-saddle fold); n1 and n3 agree at ν = 25 and 30. Reference zeros cannot occur in the window.
> * **J = 29 twist edge:** measured Γ ≤ 310 c₂ on ν = 25–47; ν = 10–20 is consistent with 0 (σ up to 730). Signed window: L 0.61 → **1.23**, λ ≥ +0.18. The remaining model input is |Γ| ≤ 10⁴ on ν < 10; the state tolerates 4.5·10⁴ there, 16× the largest known value. Masses 1.04–1.10Λ: L 2.2–5.7. **A sign-free ×3 envelope does not close (L = 0.94)**, because of Monte Carlo noise at ν = 10–20.
> * **J = 21:** a new 34-node direct bulk row puts the twist edge at +0.998 ± 0.020, L 0.58 → 1.20 (1.7–6.2 at 1.03–1.10Λ). **J = 16:** L → 4.4 / 3.7 / 2.6 / 1.7. **J = 12:** L → 7.4 / 7.5 / 3.6.
> * **J = 23 / 27 / 31:** window cost now 0.15 / 0.29 / 0.46, but the K₂ λ-level truncation blocks them (with calibrated truncation L ≈ 1.5 / 1.4). **J = 32:** the window costs 0.26 measured, vs 2.31 under the convention. J = 20 and J = 32 have no budget.
> * **Cautions:** the window weight grows steeply with J (1.4 → 23·10⁻⁵ N_s/c₂ from J = 12 to 32), so for J ≳ 32 it is a live input again. Lemma C/C′ cannot bound Γ (it is beyond leading Regge). HIGHSPIN's odd-J window weight omits a gauge factor (8.52 vs 8.73·10⁻⁵). EE_DG1000_LARGEJ's "direct RQMC cannot reach the precision" is false on ν = 25–50.

## 0. Verdict

| item | result | grade |
|---|---|---|
| **Binding state** (`J = 5`, `+−`, parity-odd, twist edge `1.0049925Λ`, Δ_g = 1000) | On its **own-mass K₂ row** (`Δ = 1005`, 36 nodes `x = 0.025 … 0.9`): `λ/N_s = +0.20088` (hybrid `+0.20090`, published chain `+0.20088`). **Error budget (sign-free, ×3 on the model terms): `2.6·10⁻³`; `λ/N_s ≥ +0.1983`, leverage `L = 78`.** With the λ-level truncation estimate: budget `1.5·10⁻³`, `L = 130`. Base: certified flat floor minus the Lemma-C curvature bound gives `≥ +0.19993 N_s` [certified + P], so `λ/N_s ≥ +0.197` with no quadrature input. Checks: direct n1 RQMC nodes on the same row agree with K₂ node by node (9 nodes, `|Δ| ≤ 13 c₂`, all but one within 1.2σ) and move `λ` by `< 10⁻⁶`; the parity-odd analytic-radial engine agrees with K₂ on `ν = 12.5–38` to `19–99 c₂`, inside its own last-order term. | **[NE±]** (was [NE + M]) |
| Binding state at Δ_g = 300 (J = 5, `+−`, `1.0165851Λ`) | Own row `Δ = 305`: `+0.20377` (published `+0.20384`, RD-S `+0.2038`). Budget `5.8·10⁻²` (dominated by the sign-free `ν < 50` window with `G = 10⁴ c₂`: `3.6·10⁻²`), `λ/N_s ≥ +0.146`, `L = 3.5` (λ-level `4.3`). No direct or engine check at 300. | **[NE±]**, `L = 3.5` |
| Why the small-ν engine route does not bind | At Δ_g = 1000 the near-gap states (J = 4–15, twist edge) carry `≈ 0.1 %` of their sign-free weight at `ν < 20` and `≈ 0.5 %` at `ν < 50`; `ν < 50` costs `3·10⁻⁴ N_s` at the binding state even for a uniform `|Γ| = 10⁴ c₂`. The decisive window is `x ∈ [0.05, 0.9]`, where K₂ converges (`J²/m ≤ 0.26`). POS_ERRORS §3.1 also treated odd J as a 2×2 (`++`, `+−`) block; at odd J only the parity-odd `+−` exists, and the 1×1 sign-free tolerance at the binding state is `3.6·10⁴ c₂` (not `6.6·10³`). | [D] (weights) + [NE] |
| **H_M, J = 4, 6 near gap** | Own rows `Δ = 1004, 1006`: chain unchanged to `≤ 10⁻⁴` (`+0.20413`, `+0.22953` at the twist edges). Budget sign-free: `L ≥ 22` (J = 4), `≥ 5.0` (J = 6); λ-level `L ≥ 80 / 28`. J = 6: 5 direct nodes on the own row agree with K₂ (≤ 1.9σ), state shift `−2·10⁻⁶`. | **[NE±]** |
| H_M, J = 8 | Own row `Δ = 1008`: `+0.42973` (published chain `+0.42930`). Sign-free budget `L = 1.1–1.8`; λ-level `L = 4.3–8.4`. 4 direct nodes agree with K₂ (≤ 2.6σ); state shift `+7·10⁻⁵` (twist), `−1.6·10⁻⁴` (1.02Λ). | [NE±] on the λ-level budget only; sign-free marginal |
| H_M, J = 12 | Own row `Δ = 1012`: twist edge `+0.84968` (published `+0.84326`); `1.05Λ` `+0.68574`. Sign-free budget fails (`L = 0.3`); λ-level `L = 1.6–3.1`. | [NE + M] (own row removes the 1/m interpolation; truncation not bounded with ×3 margin) |
| H_M, J = 16 | Own row `Δ = 1016`: twist edge `+0.93421` (published `+0.93981`), `1.02Λ` `+0.76316`, `1.03Λ` `+0.74581`. λ-level `L = 1.0–1.3`; sign-free `L = 0.1`. | [NE + M] (own row; not error-bounded) |
| **H_cont, J = 15 twist edge** (`1.0149042Λ`) | **Computed: `+0.6900` on its own row `Δ = 1015`** (36 nodes); `+0.6751` with the published rows (Δ = 336/400, held at the 400 row); hybrid `+0.6946`. λ-level budget `L = 1.3` (`≥ +0.16`); sign-free does not close. Not binding (×3.4 above the binding state). | **[NE + M]** (was [M]) |
| **H_edge** (direct n1 edge nodes `x = 0.8 … 0.9`, injected as HC does) | Computed on even J = 18 (rows 800 + 1600), 24 (500 + 1600) and odd J = 19, 25, 29 (rows 1000 + J). Shift of the twist-edge state: J = 18 `−0.0066`, J = 19 `+0.0041`, J = 24 `+0.024`, J = 25 `+0.032` (jackknife `≤ 4·10⁻⁴`). **J = 29: not resolved** — the edge window carries `≈ 1.2 N_s` sign-free at the measured `|Γ| ≈ 2.8·10³ c₂`; in the w2 chain the K₂ node at x = 0.8 is an outlier (C3H `−1.8·10⁴` vs direct `+2.4·10³ ± 0.4·10³`) that drives the state to `−3.6`, and with the direct nodes it is `+1.39 ± 0.03`, i.e. the value is set by the edge representation. Not run: J = 20, 22, 26–28, 30–32, odd 21, 23, 27, J = 48. | J = 18, 19, 24, 25: [NE + M] (was [M]); J = 29: [M], representation-dominated |
| bj engine at even J = 4–16 | **Not run.** It would decide only `ν < 12–20`, which carries `≈ 0.1 %` of the weight of every near-gap state at Δ_g = 1000 (row 2); the compute went to own-mass rows. The odd-J analogue at the binding state was run (row 1). | — |
| **Light J = 0, 2 under Ω^EE** (hybrid: exact flat image + leading-Regge AdS curvature) | Positive at every `Δ_L` from 2.6 (J = 0) / 4.5 (J = 2) to `Δ_g − 1`, both Δ_g: `F/N ≥ 0.788` (J = 0), `≥ 0.638` (J = 2); the curvature moves `F/N` by `≤ 1.8·10⁻⁴`. **Lemma C certifies the leading-Regge image** (flat [P] minus the sign-free curvature bound `> 0`) for `Δ_L ≥ 4.5` (J = 0, 1000), `≥ 3.5` (J = 0, 300), and every computed J = 2 state (`Δ_L ≥ 4.5`, both Δ_g). | flat [P]; leading-Regge image [P given the LR form] + [NE]; beyond-LR: below |
| Light J = 0, 2 just below the gap, beyond leading Regge (K₂ rows) | **Own rows at `Δ_L = 990`** (1000): `λ/N = +1.0000`, `L = 4362` (J = 0), `913` (J = 2). Published rows (1/m-interpolated; `x > 0.8`, which runs past `x = 1` for light states, sign-free): `Δ_L ≥ 0.9Δ_g`: `L ≥ 378`; `0.7Δ_g`: `L ≥ 28` (1000). At 300 (rows held at Δ = 400): `L ≥ 13` (J = 0), `≥ 8.8` (J = 2, `Δ_L ≥ 270`), `1.4` (J = 2, `Δ_L = 200`). | `Δ_L = 0.99Δ_g`: **[NE±]**; `≳ 0.7Δ_g`: [NE + M] (was [I]) |
| **Repaired Q on light J = 0** (`Ω^Q_up + 0.038034 Ω^EE_θ`) | Flat margin is `F/N ≈ 3.3 (m/Λ)²`, **zero at `m → 0` by construction of ε\***. Leading-Regge model: positive at every computed `Δ_L ≥ 2.6`, but the margin there (`2·10⁻⁵`, 1000; `2·10⁻⁴`, 300) is the AdS curvature itself, and Lemma C cannot sign it below `Δ_L ≈ 30` (1000) / `15` (300). For `1 < Δ_L < 2` the heavy-state convention gives `m² < 0` and a flat value `< 0` (of size `10⁻⁶`). **The sign for `Δ_L = O(1)` is decided by `O(Λ⁻²)` AdS terms that no computation controls.** Remedy: `ε = 0.040` gives a flat margin `0.022 N` at `m → 0` (×100 above every computed curvature effect) for `+3.5` (1000: 160.99 → 164.51) / `+2.6` (300: 123.08 → 125.71) on the upper end. | `Δ_L ≥ 30`: [P given LR] + [NE]; `Δ_L = O(1)` with ε\*: **[I]**; with `ε = 0.040`: [D-sketch] |

**Failures.** No physical state is found negative. One chain variant is: the J = 29 twist edge at Δ_g = 1000 in the plain w2 odd chain
(`−3.62`), caused by a single K₂ outlier at x = 0.8 that direct RQMC replaces (`+1.39 ± 0.03`); the edge window decides that state (§1.6).
The sign-free budgets do not close from J = 12 on, and the repaired Q functional has no robust margin on very light scalars.

---

## 1. Job 1: near-gap rows at Δ_g = 1000

### 1.1 Where the weight is

`W(ν-window)` = cost in `λ/N_s` of a uniform sign-free `|Γ| = 1 c₂` (k = 2 / k = 3 parts), hybrid states at the twist edge, 1024 nodes
(`ee1000_common.hyb`, the chain's weights):

| state | `ν < 20` | `20 ≤ ν < 50` | `50 ≤ ν < 0.5m` | `0.5m ≤ ν < 0.8m` | `0.8m ≤ ν < 0.9m` | `ν ≥ 0.9m` | `G_tol = λ/ΣW` |
|---|---|---|---|---|---|---|---|
| J = 4 `+−` | 2.6·10⁻⁹ | 1.3·10⁻⁸ | 1.3·10⁻⁶ | 1.5·10⁻⁶ | 1.8·10⁻⁷ | 2.2·10⁻⁸ | 6.8·10⁴ |
| J = 4 `++` | 3.0·10⁻⁸ | 1.6·10⁻⁷ | 9.8·10⁻⁶ | 2.7·10⁻⁶ | 2.2·10⁻⁸ | 7·10⁻¹¹ | 1.6·10⁴ |
| **J = 5 `+−`** | **4.7·10⁻⁹** | 2.4·10⁻⁸ | 2.4·10⁻⁶ | 2.8·10⁻⁶ | 3.3·10⁻⁷ | 3.7·10⁻⁸ | **3.6·10⁴** |
| J = 6 `+−` | 9.6·10⁻⁹ | 4.9·10⁻⁸ | 4.8·10⁻⁶ | 5.8·10⁻⁶ | 6.5·10⁻⁷ | 6.9·10⁻⁸ | 2.0·10⁴ |
| J = 8 `+−` | 1.7·10⁻⁷ | 8.6·10⁻⁷ | 8.4·10⁻⁵ | 9.8·10⁻⁵ | 1.1·10⁻⁵ | 1.0·10⁻⁶ | 2.2·10³ |
| J = 15 `+−` | 1.5·10⁻⁶ | 7.6·10⁻⁶ | 7.4·10⁻⁴ | 8.3·10⁻⁴ | 7.4·10⁻⁵ | 5.1·10⁻⁶ | 4.2·10² |

**Consequences.**
* The window that the 1/Δ engine reaches (`ν ≲ 12–20`) carries `≈ 0.1 %` of the weight at the binding state, and at every tabulated J; even `ν < 50` is `≈ 0.5 %`.
  An engine value there, with any error bar, cannot change a near-gap verdict at Δ_g = 1000. (At Δ_g = 300 the same window carries
  `≈ 1 %` (J = 5, 15); it was decisive there only at J ≥ 24, where `N_s` is small.)
* The decisive region is the bulk `x = ν/m ∈ [0.05, 0.9]`. There K₂ (the Laplace series of the crystal integral at the saddle) is in its
  convergent regime (`J²/m ≤ 0.26` for J ≤ 16), and its own next order is a truncation estimate.
* So the cheapest error-bounded route is not "engine at small ν" but **K₂ rows at the state's own mass** (which removes the chain's 1/m
  interpolation between rows) with the Laplace truncation as error bar, plus direct RQMC nodes as a check.
* **A slip in POS_ERRORS §3.1.** `large_m_map.state` gives every `J ≥ 4` both helicities and an off-diagonal. At odd J the only
  structure is the parity-odd `+−` (UT_v2 §2.5; A1PRIME §2.1), so the J = 5 row of POS_ERRORS's table (`G_tol = 6.6·10³`, `L = 1.8`)
  contains a non-existent `++` channel. For the 1×1 block, `G_tol = 3.6·10⁴ c₂` at the twist edge, and the E1′ leverage
  (`g(5) = 2850 c₂`) is `≈ 12`, not 1.8 [NE]. This alone does not make the state error-bounded (E1′'s bulk envelope is a model);
  §1.3 does.

### 1.2 Method (`ng_rows.py`, `ng_chain.py`, `ng_bind.py`)

* **Rows.** K₂ crystal points (w1_rows, unchanged: `L1_DEG = 6`, odd J → parity-odd `p3` integrand with `E` in place of `q`, even J →
  w0 EE integrand; C3H for the k = 3 part; state and `Δ_ref = 2·10⁴` at the same ν), at the twist-edge dimension `Δ = 1000 + J`,
  `x = 0.025, 0.05, …, 0.9` (36 nodes, both kinds). Cost 7–8 s per point at J = 5, 18 s at J = 12, 21 s at J = 15.
* **Chain.** The EE_NEARGAP main chain (w2 in EE mode via `ee1000_common`: Jacobi fix, EE helicity basis, local de-shift, Hampel despike,
  hybrid states) runs unchanged; the only change is that its row loaders also see `neargap2/cache`. A state at the twist edge then sits
  on its own row (`t = 0` in the 1/m interpolation). `NG_NOOWN=1` reproduces the published chain: `+0.20088` at the J = 5 twist edge
  (REPRO: hybrid `+0.20090`, correction `≤ 10⁻³`).
* **Error budget of a state** (all sign-free unless stated; `λ_c` = central value with the two outer windows removed):
  * `E_trunc` = 3 × Σ_nodes `|w|·|Γ_{K=2} − Γ_{K=1}|` on the own row's x-range (the K = 2 Laplace term, node by node, sign-free).
    For even J the order is set inside `pred` (the K = 2 records are needed to load the rows).
  * `E_trunc^λ` = 3 × `|λ_{K=2} − λ_{K=1}|` (the signed effect of the last order kept; an alternative, less conservative estimate).
  * `E_interp` = 3 × max(`|λ_spline − λ_linear|`, `|λ − λ_every-other-node|`).
  * `E_low` = `G_low × W(ν < 50)` with `G_low = 10⁴ c₂`, each node weight multiplied by `max(1, |Re X_min/ψ₂|)` (≤ 6.8 at ν < 5, ≤ 1.24 on
    5–12, 1 above 25) so that it holds in the functional's own min-norm gauge, not only in the chain's local form.
  * `E_edge` = `G_edge × W(x > x_max)` with `G_edge = max(10⁴, 10 × max |Γ| of the row's last four nodes)`.
  * Base: `λ_hyb` (flat with exact Wigner functions + LR curvature, quadrature) and, independently, the certified flat floor `N_f/5`
    (CERT_CLOSE) minus the Lemma-C bound on the curvature (POS_ERRORS §2.1).
  * Off-diagonal (even J): the same with `c_off = mΓ_off`; the bound on `λ` is the Frobenius norm of the normalised error matrix.
  * Leverage `L = λ_c / (sum of the four errors)`.

### 1.3 The binding state (J = 5, `+−`, Δ_g = 1000)

**Row values** (`cache/ng_rowdump_J5_D1005.txt`; `Γ` in `c₂`, K = 2 [K = 1, K = 0]): C2 rises from `−8` (x = 0.025) to `+22` (x = 0.4)
and falls to `+0.1` (x = 0.9); C3H from `−4` to `−39` (x = 0.45) to `−0.4`. The Laplace series converges by about ×10 per order (e.g.
x = 0.4, C2: `+2150 → −190 → +22`). `|Γ| ≤ 40 c₂` everywhere on the row.

**State values and budget** (`cache/ng_bind_J5.txt`):

| `m/Λ` | hybrid | published chain | own row (K = 2) | K = 1 | `E_trunc` | `E_trunc^λ` | `E_interp` | `E_low` | `E_edge` | lower bound | `L` (sign-free / λ-level) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| **1.0049927** (twist edge) | +0.20090 | +0.20088 | **+0.20088** | +0.20059 | 1.9·10⁻³ | 8.7·10⁻⁴ | 2·10⁻⁷ | 3.1·10⁻⁴ | 3.7·10⁻⁴ | **+0.1983** / +0.1993 | **78** / 130 |
| 1.01 | +0.20212 | +0.20210 | +0.20210 | +0.20181 | 1.9·10⁻³ | 8.9·10⁻⁴ | 2·10⁻⁷ | 3.2·10⁻⁴ | 2.4·10⁻⁴ | +0.1996 / +0.2007 | 81 / 139 |
| 1.02 | +0.20562 | +0.20559 | +0.20560 | +0.20528 | 2.0·10⁻³ | 9.5·10⁻⁴ | 2·10⁻⁷ | 3.3·10⁻⁴ | 1.3·10⁻⁴ | +0.2032 / +0.2042 | 85 / 146 |
| 1.05 | +0.22621 | +0.22620 | +0.22620 | +0.22581 | 2.1·10⁻³ | 1.2·10⁻³ | 2·10⁻⁷ | 3.8·10⁻⁴ | 2·10⁻⁵ | +0.2237 / +0.2246 | 89 / 146 |

* Base at the twist edge: flat `+0.20090` (`N_f/N_s = 1.0000`), LR curvature `+7.8·10⁻⁷ N_s`, Lemma-C bound `6.8·10⁻⁵ N_s`
  (`b̂ = 10.95`). Certified floor route: `(N_f/5 − Lemma C)/N_s = +0.19993` [certified + P]. With the row correction (`dF/N_s = −2.2·10⁻⁵`)
  and the budget: **`λ/N_s ≥ 0.19993 − 0.00002 − 0.0026 = +0.197`**.
* **Direct check on the same row** (`ng_points.py` = hc_points n1, adaptive, `2¹¹·2·8` points, shared seeds; `ng_val.py`,
  `cache/ng_val_J5_1005.txt`; `Γ` in `c₂`, jackknife over 8 replicas):

  | x | C2: K₂ / direct | C3H: K₂ / direct |
  |---|---|---|
  | 0.1 | +9.4 / +22.1 ± 3.9 | −8.1 / −7.9 ± 6.9 |
  | 0.3 | +20.3 / +17.2 ± 16.5 | −31.1 / −21.3 ± 12.5 |
  | 0.5 | +19.8 / +11.7 ± 10.1 | −36.9 / −27.4 ± 11.8 |
  | 0.7 | +6.4 / +1.8 ± 11.0 | −14.7 / −12.3 ± 12.4 |
  | 0.8 / 0.825 / 0.85 / 0.875 / 0.9 | +1.7 / +1.0 / +0.6 / +0.3 / +0.1 vs −1.3 / −1.6 / −1.7 / −1.5 / −0.8 (± 7.1 … 2.1) | −4.6 … −0.4 vs −3.3 … +0.1 (± 8.0 … 2.3) |

  One node is off by 3.3σ (x = 0.1, C2, `+13 c₂`); the rest agree within 1.2σ. Replacing the K₂ nodes by the direct ones changes `λ/N_s`
  by `< 10⁻⁶ ± 10⁻⁶` (twist edge and 1.02Λ). The edge window `x ∈ [0.8, 0.9]` is `|Γ| ≤ 2 c₂` directly; `G_edge = 10⁴` for `x > 0.9`
  is ×5000 above it.
* **Engine at small ν** (`ng_engine5.py`: po_ee_odd EE2, depth 4, K = 3, state 1005 and reference 2·10⁴; `ng_engcmp.py`,
  `cache/ng_engcmp.txt`). In the normalisation-free form `m²[D(ν) − D(ν₀)]`, `D = log(state/ref)`, `ν₀ = 0.05m`:

  | ν | 12.5 | 15.1 | 17.6 | 20.1 | 25.1 | 37.6 |
  |---|---|---|---|---|---|---|
  | engine (Kt = 3) | 25128 | 24371 | 23497 | 22496 | 20102 | 11771 |
  | K₂ | 25029 | 24304 | 23446 | 22454 | 20070 | 11752 |
  | engine − K₂ | +99 | +67 | +51 | +42 | +32 | +19 |
  | engine last order (Kt3 − Kt2) − its value at ν₀ | +593 | +575 | +555 | +531 | +475 | +279 |

  The two independent expansions agree to `≤ 99 c₂` where both exist, well inside the engine's own last-order term (whose ν-variation on
  this range is `≈ 300 c₂`). Below `ν = 12.5` only the engine exists; its truncation is `≈ 600 c₂`, almost all ν-flat (a normalisation).
  This supports `G_low = 10⁴ c₂` by a factor ≥ 15; the window costs `3·10⁻⁴ N_s` at that size anyway.
* **What remains a model input**, with its safety factor:
  1. *K₂ equals the crystal integral* (the Laplace series of the real-section integral at the saddle): checked against direct RQMC at 9
     nodes (above) and against the analytic-radial engine at 6 [NE]; RD for the engine (UT_v2 §5.2) is not used.
  2. *Truncation*: next order ≤ 3 × the last one kept, node by node. The observed ratio per order is ≈ 0.1, so the factor is ≈ 30 in
     practice.
  3. *`|Γ| ≤ 10⁴ c₂` at `ν < 50` and at `x > 0.9`*: ≥ ×15 (engine) and ≥ ×5000 (direct edge) above the measured values.
  4. *x-interpolation*: ×3 on the spline–linear and subsample differences (`2·10⁻⁷`).
  5. *Local vs min-norm de-shift at odd J* (POS_ERRORS D9): covered by the gauge factor in `E_low` (the weights differ only at ν ≤ 25).

**Verdict (binding state).** `λ/N_s(J = 5, +−, twist edge) = +0.2009`, with every AdS correction bounded: `≥ +0.198` (hybrid base) or
`≥ +0.197` (certified flat floor and Lemma C), `L = 78` against a sign-free budget with ×3 on every model term. **[NE±].** The same holds
on the J = 5 states to 1.05Λ (`L ≥ 78`). This replaces POS_ERRORS §3.1's "not error-bounded" for this state.

### 1.4 Other near-gap rows (own-mass K₂ rows, `cache/ng_bind_J*.txt`)

| J (own Δ) | `m/Λ` | hybrid | published chain | own row | `E_trunc` / `E_trunc^λ` | `E_low` | `E_edge` | lower bound (sign-free / λ-level) | `L` | grade |
|---|---|---|---|---|---|---|---|---|---|---|
| 4 (1004) | 1.0040 (twist) | +0.20413 | +0.20413 | +0.20413 | 6.9·10⁻³ / 8·10⁻⁵ | 2.3·10⁻³ | 2·10⁻⁴ | +0.1948 / +0.2016 | 22 / 80 | **[NE±]** |
| 4 | 1.01, 1.02 | +0.21028, +0.22039 | same | same | ≤ 6.7·10⁻³ | ≤ 2.2·10⁻³ | ≤ 10⁻⁴ | ≥ +0.2012 | ≥ 23 | [NE±] |
| 6 (1006) | 1.0060 (twist) | +0.22946 | +0.22954 | +0.22953 | 3.8·10⁻² / 3.3·10⁻⁴ | 6.7·10⁻³ | 1.2·10⁻³ | +0.1833 / +0.2213 | 5.0 / 28 | **[NE±]** |
| 6 | 1.01, 1.02 | +0.25070, +0.31000 | +0.25078, +0.31007 | same | ≤ 3.8·10⁻² | ≤ 6.6·10⁻³ | ≤ 10⁻³ | ≥ +0.2057 | ≥ 5.6 | [NE±] |
| 8 (1008) | 1.0080 (twist) | +0.42903 | +0.42930 | +0.42973 | 0.28 / 2.5·10⁻² | 3.1·10⁻² | 3.8·10⁻² | +0.080 / +0.336 | 1.2 / 4.6 | [NE±] (λ-level) |
| 8 | 1.02, 1.05 | +0.57870, +0.43458 | +0.57835, +0.43441 | +0.57829, +0.43448 | ≤ 0.45 / ≤ 4.9·10⁻² | ≤ 4.2·10⁻² | ≤ 4.3·10⁻² | ≥ +0.039 / ≥ +0.38 | ≥ 1.1 / ≥ 4.3 | [NE±] (λ-level) |
| 12 (1012) | 1.0119 (twist) | +0.84949 | +0.84326 | +0.84968 | 3.0 / 0.21 | 0.14 | 0.18 | < 0 / +0.319 | 0.3 / 1.6 | [NE + M] |
| 12 | 1.02, 1.05 | +1.00000, +0.68671 | +0.99511, +0.68545 | +0.99786, +0.68574 | ≤ 3.0 / ≤ 0.21 | ≤ 0.14 | ≤ 0.10 | — / ≥ +0.46 | — / ≥ 2.2 | [NE + M] |
| 16 (1016) | 1.0159 (twist) | +0.93879 | +0.93981 | +0.93421 | 8.3 / 0.48 | 0.32 | 0.13 | < 0 / −0.010 | 0.1 / 1.0 | [NE + M] |
| 16 | 1.02, 1.03, 1.05 | +0.76702, +0.74987, +1.00000 | +0.76687, +0.74823, +0.99573 | +0.76316, +0.74581, +0.99424 | ≤ 6.9 / ≤ 0.47 | ≤ 0.28 | ≤ 0.08 | — / ≥ +0.008 | — / ≥ 1.0 | [NE + M] |
| 15 (1015, odd) | 1.0149 (twist) | +0.69459 | +0.67511 | +0.68995 | 4.6 / 0.29 | 0.10 | 0.15 | < 0 / +0.162 | 0.1 / 1.3 | [NE + M] |
| 15 | 1.02, 1.03, 1.05 | +1.0, +1.0, +0.94824 | +0.97416, +0.97812, +0.93223 | +0.99371, +0.99440, +0.94394 | ≤ 6.5 / ≤ 0.42 | ≤ 0.14 | ≤ 0.15 | — / ≥ +0.29 | — / ≥ 1.4 | [NE + M] |

**Direct checks on the own rows** (`cache/ng_val_J6_1006.txt`, `ng_val_J8_1008.txt`; n1 RQMC, `Γ` in `c₂`, 2×2 helicity entries):
* J = 6, x = 0.2, 0.5, 0.8, 0.85 (+0.9, where the K₂ node was despiked): all 16 entries within 1.9σ (largest: `+10 ± 5`, C3H `++`, x = 0.5;
  `+61 ± 38`, C3H `++`, x = 0.85). Injecting the direct nodes moves the state by `−2·10⁻⁶ ± 6·10⁻⁶` (twist edge) and `−3·10⁻⁶` (1.02Λ).
* J = 8, x = 0.2, 0.5, 0.8, 0.85: 13 of 16 entries within 1.3σ; three `++` entries at 1.8–2.6σ
  (`−20 ± 8` C3H and `−13 ± 6` C2 at x = 0.2, `−56 ± 30` C3H at x = 0.8), all with direct below K₂. State shift
  `+6.7·10⁻⁵ ± 4.8·10⁻⁵` (twist), `−1.6·10⁻⁴ ± 1.0·10⁻⁴` (1.02Λ).
* So "K₂ = crystal integral" is checked at J = 5, 6, 8 at the level of `≤ 2·10⁻⁴ N_s`; J = 4, 12, 15, 16 have no direct check here.

**Reading.**
* **H_M is confirmed where it was used.** Own row vs published chain (1/m interpolation between rows 400/800/1600): `≤ 10⁻⁴` at
  J = 4–6, `4·10⁻⁴` at J = 8, `+6·10⁻³` at J = 12 (published chain pessimistic), `−5.6·10⁻³` at J = 16 (optimistic), `+1.5·10⁻²` at
  J = 15 (published: held 400 row). All shifts are `≤ 0.015 N_s`, against margins `≥ 0.2`. [NE]
* **Error bounding stops between J = 8 and J = 12.** The sign-free tolerance falls steeply with J (`G_tol = 6.8·10⁴ → 4.2·10² c₂`
  from J = 4 to 15, §1.1; `N_s ∝ b̂⁻³`), while the K₂ terms themselves stay O(10²–10³) `c₂`; from J = 8 on, only the signed λ-level truncation closes (the corrections cancel across x, as
  POS_ERRORS §3.3 found), and at J = 12, 15, 16 even that has only ×1.0–1.6. The K₂ λ-level step K = 1 → 2 is
  `3·10⁻⁵ / 1·10⁻⁴ / 8·10⁻³ / 7·10⁻² / 0.10 / 0.16` at J = 4 / 6 / 8 / 12 / 15 / 16: the Laplace series slows down as `J²/m` grows
  (0.016 → 0.25). Error-bounding J ≥ 12 would need the next Laplace order (`L1_DEG = 8` records, M6) on the own rows.
* `G_low = 10⁴ c₂` is generous; at J ≥ 8 it costs 0.03–0.32 `N_s`. A measured small-ν bound (e.g. the J = 5 engine, ≲ 10² `c₂`
  away from the normalisation) would cut `E_low` ×10–100, but it was measured only at J = 5.

### 1.5 H_cont: the J = 15 twist edge

* **Own row** `Δ = 1015` (`ng_rows.py`, 36 nodes `x = 0.025 … 0.9`, EE + C3H, 21 s per point): twist edge **`+0.68995`**, `1.02Λ`
  `+0.99371`, `1.03Λ` `+0.99440`, `1.05Λ` `+0.94394` (§1.4 table). The state is computed on a row at its own mass; the remaining model
  inputs are those of every K₂ row; the K = 1 → 2 step moves it by `0.10`, so the λ-level budget gives only `≥ +0.16` (`L = 1.3`).
  **H_cont's exception at J = 15 is removed: the twist edge is a computed state, `+0.690` [NE + M], ×3.4 above the binding state.**
* Published rows (odd J = 15 rows exist only at Δ = 336, 400; a 1000-gap state is held at the 400 row), same chain, mass grid from the
  twist edge: `+0.6751` at `1.0149042Λ` (hybrid `+0.6946`), `+0.9742` (1.02), `+0.9322` (1.05), **`+0.8035` (1.07)**, `+0.9375` (1.10).
  This is the "computed" version of POS_ERRORS D6 in the published protocol [NE + M]. The held-row budget does not close (`L ≈ 0.1`).

### 1.6 H_edge

Direct n1 points (`ng_points.py` = hc_points n1, adaptive, `HC_M = HC_MREF = 11`, 8 replicas, shared seeds; RESULT_D_CLEANUP §2) at
`x = 0.8, 0.825, 0.85, 0.875, 0.9` on the rows that carry the 1000 near-gap states, injected into the chain exactly as `hc_chain`'s
`HC_INJ` does (`hc_nodes.to_record`; the injected nodes exempt from the Hampel rule, `HC_DESPIKE_INJ = 0`; nodes with a non-finite
entry dropped); `ng_val.py`, `cache/ng_val_J*.txt`. Odd J = 19–29 sit on the single ee1000 row `Δ = 1000 + J` (EE_DG1000_LARGEJ), read
here through the w2 odd chain (not the ee1000 shape chain that produced the published odd-J values).

| J (rows) | twist edge `m/Λ` | without → with edge nodes | shift | other masses (shift) | largest direct `|Γ|` in the window (`c₂`) | K₂ vs direct at x = 0.8 |
|---|---|---|---|---|---|---|
| 18 (800, 1600) | 1.01786 | +0.7399 → +0.7333 ± 0.0003 | −0.0066 | 1.03: −0.0048; 1.05: −0.0018 | 1.3·10³ | `|Δ| ≤ 172` (800); the 1600 K₂ node at 0.8 is Hampel-dropped |
| 19 (1019) | 1.01884 | +0.9954 → +0.9995 ± 0.0003 | +0.0041 | 1.03: +0.0017; 1.05: −0.0002 | 5.6·10³ | −131 / +122 (≈ 3σ) |
| 24 (500, 1600) | 1.02374 | +0.9454 → +0.9699 ± 0.0002 | +0.024 | 1.03: +0.017; 1.05: +0.006; 1.1: +0.001 | 1.8·10³ | row 500 (`J²/m = 1.15`): `++` off by `−1006 / −1054`; row 1600: `|Δ| ≤ 210` |
| 25 (1025) | 1.02472 | +0.8918 → +0.9234 ± 0.0004 | +0.032 | 1.03: +0.026; 1.05: +0.008 | 3.7·10³ | −428 / +472 |
| 29 (1029) | 1.02862 | **−3.62 → +1.39 ± 0.03** | (+5.0) | 1.04: −1.83 → +1.24; 1.06: −0.21 → +1.15 | 2.8·10³ | C3H: K₂ −17 933 vs direct +2370 ± 401 |

Sign-free weight of the window (`ng_edgew.py`, hybrid weights, twist edge): `W[0.8, 0.9] = 1.2·10⁻⁴` (J = 18, `+−`), `1.9·10⁻⁴`
(19), `3.0·10⁻⁴` (24, 25), `4.3·10⁻⁴` (29) per `c₂`; `W(0.9, Λ/m] = 0.7–2.1·10⁻⁵`. At the measured `|Γ|` the sign-free cost of
`[0.8, 0.9]` is `0.35` (J = 18) to `1.2 N_s` (J = 19–29): the window is decided by the oscillation of `Γ` (signs alternate node to node
with a period ≈ 0.05–0.1 in x), not by its size. At J = 5 the same window costs `< 10⁻⁴ N_s`.

**Reading.**
* J = 18, 19, 24, 25: the direct edge values move the near-gap states by `−0.007 … +0.032 N_s`, jackknife `≤ 4·10⁻⁴`; the two
  sizeable moves are favourable (as RESULT_D_CLEANUP found at 300). With the remaining spline representation between 5 nodes per
  window, this is [NE + M], the same grade as the computed rows at lower J (was [M]).
* **J = 29 is not resolved, and is the one place where a chain variant gives a negative value.** On the (1029, 29) row the K₂ series is
  slow (`J²/m = 0.83`; K = 0 → 1 → 2 steps of `1.8·10⁴ → 4·10³ → 9·10²` `c₂` at x = 0.75) and its x = 0.8 C3H node is off by
  `2·10⁴ c₂`. The w2 odd chain has no Hampel filter on odd rows (hc_chain and the ee1000 shape chain do), so that node turns the
  twist-edge state to `−3.62`; with the direct nodes it is `+1.39 ± 0.03`. A value `> 1` means the edge window dominates the
  correction. The published odd-J values (`+0.89 … +0.98`, EE_DG1000_LARGEJ / HC) come from the ee1000 shape chain, whose window stops at
  x = 0.8. **So `H_edge` at J = 29 (and plausibly at the other odd J ≳ 27 on ee1000 rows) is representation-dominated: the sign of the state depends on how `x ∈ [0.8, Λ/m]` is modelled
  between nodes spaced at a quarter of the oscillation period.** It stays [M]. The K₂ node at x = 0.8 on rows with `J²/m ≳ 0.8` should
  not be trusted (also J = 24, row 500: `++` off by `10³ c₂`).
* Not run (time): even J = 20, 22, 26, 28, 30, 32; odd 21, 23, 27; J = 48 (fw chain). J = 22 at 1000 still draws on the Δ = 336 row only,
  where `x = 0.9` is within 10 % of Δ (the RQMC fails there; RESULT_D_CLEANUP §0.4).

---

## 2. Job 2: light single traces (A1′) at finite Δ_L

### 2.1 Set-up (`ng_light.py`)
* A light single trace of spin J ∈ {0, 2} enters the k = 2, 3 sum rules with its block action (A1PRIME §3), i.e. as a state of that
  dimension and spin. Only the `++` structure exists for J < 4, so its image is a 1×1 entry and only its sign matters.
* **Mass convention.** The heavy states use `m² = α(α + 2J + 1)`, `α = Δ − J − 2`. Light states are parametrised the same way, so
  `m > 0` needs `Δ_L > J + 2` (J = 0: `Δ_L > 2`; J = 2: `Δ_L > 4`). Below that, only the `m → 0` limit is reported.
* **Images.** Exact flat image (NF02 weights; Q3 weights for Q; `d^J_{00} = P_J(1 − 2y)`, `y = p²Λ²/m²`) plus the leading-Regge AdS
  curvature `2ψ₂(P_{−½+iν}(cosh β) − J₀(νB)) + uψ₃(…)`, `cosh β = 1 + (2J+1)²/2m²`, exactly as the heavy hybrid (`ee1000_common.hyb`).
  `P` by the Mehler–Dirichlet integral (agrees with mpmath `hyp2f1` to `7·10⁻⁸` where that converges; it does not for `β ≳ 2`, ν ~ 10³).
  The node count resolves `νB ≤ Λ(2J+1)/m`; doubling it changes the curvature by `≤ 10⁻¹³ N`.
* **Flat check** (exact polynomials of A1PRIME §4.2–4.3, `m/Λ` = 0.1, 0.4602, 0.675, 0.99): `F_D⁰`, `F_Q⁰`, `F_D²` reproduced to `10⁻⁶`
  (e.g. `F_D²(0.675) = 0.448298`, `F_Q⁰(0.460237) = 0.000000`).

### 2.2 Ω^EE on light states (`cache/ng_light.txt`)

| | J = 0, Δ_g = 1000 | J = 0, 300 | J = 2, 1000 | J = 2, 300 |
|---|---|---|---|---|
| range of Δ_L | 2.6 … 999 | 2.6 … 299.7 | 4.5 … 999 | 4.5 … 299.7 |
| min flat `F/N` (at `m → 0`) | 0.7884 | 0.7884 | 0.6379 | 0.6379 |
| min AdS (LR) `F/N` | 0.7884 | 0.7884 | 0.6379 | 0.6379 |
| max `|LR curvature|/N` | 1.4·10⁻⁵ | 1.4·10⁻⁵ | 1.6·10⁻⁵ | 1.8·10⁻⁴ |
| Lemma C certifies `flat − |curv| > 0` for | Δ_L ≥ 4.5 | Δ_L ≥ 3.5 | Δ_L ≥ 4.5 (all) | Δ_L ≥ 4.5 (all) |

* `F/N → 1` from `m ≳ 0.3Λ` (J = 0) and `m ≳ 0.5Λ` (J = 2): both pieces positive.
* The curvature is tiny: `β ≈ (2J+1)/m` is small unless `m = O(1)`, and where it is not (`Δ_L ≲ 5`) the curvature integrand oscillates
  at frequency `νB ≈ Λ(2J+1)/m` and averages out against the smooth weights. Lemma C is sign-free and fails only at `Δ_L ≲ 4` (J = 0), where `ν(B − β)` saturates.
* **Grade.** The flat image is [P] (A1PRIME). The leading-Regge AdS image is positive at every Δ_L [P given the LR form + quadrature
  NE] for `Δ_L ≥ 4.5`, and [NE] (computed, not bounded) for `2 < Δ_L < 4.5` at J = 0. Whether the leading-Regge form is the AdS image
  of a light block at `Δ_L = O(1)` is not established: for heavy states it is the leading term of a `1/Δ` expansion, which a light state
  does not have. A1PRIME's statement (image = `m → 0` flat value × `(1 + O(Δ_L²/Λ²))`) remains [D-sketch]. For Ω^EE this does not
  matter: the margin is `F/N ≥ 0.64`, so any relative AdS correction below 60 % is harmless.

### 2.3 Beyond leading Regge, just below the gap (`cache/ng_bind_light{1000,300}.txt`)

The chain's K₂ rows for J = 0, 2 (Δ = 400, 1600; the 1000 near-gap rows) give `Γ(x)` for `x ≤ 0.8`; a light state at `m < Λ` interpolates
them in 1/m (Δ_g = 1000) or holds the 400 row (Δ_g = 300, `m < 398`). Its x-range extends to `Λ/m > 1`; `x > 0.8` is covered sign-free
with `G = 10⁴ c₂`.

| Δ_g | J | `m/Λ` (Δ_L) | `λ/N` | budget | `L` |
|---|---|---|---|---|---|
| 1000 | 0 | 0.70 (700) / 0.90 (900) / 0.95 / 0.99 / 0.999 | +1.0000 | 3.1·10⁻³ / 4.5·10⁻⁴ / 3·10⁻⁴ / 2.5·10⁻⁴ / 2.4·10⁻⁴ | 326 … 4130 |
| 1000 | 2 | same | +0.9999 … +1.0000 | 3.6·10⁻² / 2.7·10⁻³ / 1.6·10⁻³ / 1.2·10⁻³ / 1.1·10⁻³ | 28 … 900 |
| 300 | 0 | 0.66 (200) / 0.83 (250) / 0.90 / 0.97 / 0.997 | +0.9999 | 7.9·10⁻² / 3.0·10⁻² / 2.4·10⁻² / 2.0·10⁻² / 1.8·10⁻² | 13 … 55 |
| 300 | 2 | same | +0.9978 … +1.0003 | 0.71 / 0.21 / 0.11 / 7.3·10⁻² / 6.2·10⁻² | 1.4 … 16 |

**Own rows at `Δ_L = 990`** (Δ_g = 1000; `x = 0.025 … 0.9`; `Λ/m = 1.0101`): J = 0 `λ/N = +1.00000`, budget `2.3·10⁻⁴`, `L = 4362`;
J = 2 `+1.00001`, budget `1.1·10⁻³`, `L = 913` (λ-level: `4634`, `1372`). With the row at the state's own mass, the only model inputs left
are those of the binding state (K₂ = crystal, `|Γ| ≤ 10⁴ c₂` at `ν < 50` and at `x ∈ (0.9, 1.01]`). → **[NE±]** for light J = 0, 2 at
`Δ_L = 0.99Δ_g`.

**Grade.** Light J = 0, 2 states with `Δ_L ≳ 0.7Δ_g` at 1000 (`≳ 0.8Δ_g` at 300): positive with explicit budgets, inside the chain model
(1/m interpolation of the rows; `x > 0.8` sign-free) — **[NE + M]**, replacing A1PRIME's [I] "AdS action of J ≤ 2 states with twist just
below Δ_g". Lower Δ_L at 300 (`Δ_L = 200`, J = 2: `L = 1.4`) is held from the 400 row and not converted.

### 2.4 The repaired Q functional on light scalars

`Ω^Q_up + ε*Ω^EE_θ`, `ε* = 0.4304302/11.3170270 = 0.0380339`; for a light J = 0 the k = 4, 5 parts drop out (A1PRIME §3), so its image is
`F_Q⁰ + ε*F_D⁰` (smeared k = 2, 3 parts only).

* **Flat:** `F_R(m) = (2.032073 + ε*·20.771893) m² = 2.8221 m²` (m in Λ units): zero at `m = 0` by the choice of ε\*, `F/N ≈ 3.3 m²`.
* **Leading-Regge AdS curvature** (`cache/ng_light.txt`, column R): `+1.6·10⁻⁵ N` (1000) and `+1.7·10⁻⁴ N` (300) at `Δ_L = 2.6–8`,
  i.e. **larger than the flat margin `3.3(m/Λ)²` for `Δ_L ≲ 4`** (both Δ_g). It happens to be positive at every computed Δ_L
  (min `F/N = +1.8·10⁻⁵`, 1000; `+2.0·10⁻⁴`, 300), but Lemma C cannot sign it for `Δ_L < 30` (1000) / `< 15` (300), and the leading-Regge
  form is itself not established there (§2.2).
* **`1 < Δ_L < 2`** (allowed for a scalar in d = 3): `m² < 0` in the heavy convention and `F_R < 0` at order `m²/Λ² ~ 10⁻⁶`. Whether this
  is physical is again an `O(Λ⁻²)` question.
* **Conclusion.** With ε\* the repaired functional is positive on light scalars with `Δ_L ≳ 30` (1000) / `≳ 15` (300) [P given LR]
  and computed positive down to `Δ_L = 2.6` in the LR model, but its sign for `Δ_L = O(1)` rests on `O(Λ⁻²)` AdS terms: **[I]**.
* **Remedy [D, arithmetic].** Take `ε` slightly above ε\*. With `ε = 0.040` the flat margin at `m → 0` is
  `(0.040 − 0.0380339)·11.317 = 0.0222` per unit `|g₀|²` (`F/N ≈ 0.022`), ×100 above the largest computed curvature effect, so only a
  relative AdS correction `> 2 %` could flip it (A1PRIME's [D-sketch] expects `O(Δ_L²/Λ²)`). Cost (A1PRIME §4.3, `+2(ε − ε*)L_A`):
  upper end **164.51** (1000; was 160.99) and **125.71** (300; was 123.08). The `c₁²` coefficient and the κ term change in proportion
  (not recomputed here).

---

## 3. Changes to POS_ERRORS's inventory (Δ_g = 1000, Result D)

| # | item | before | after |
|---|---|---|---|
| D3 | `H_M`, binding state | [M], validated [NE]; "not error-bounded" (§3.1) | **[NE±]**, `L = 78`, lower bound `+0.198` (§1.3) |
| D3 | `H_M`, J = 4, 6 | [NE + M] | **[NE±]** (`L ≥ 22`, `≥ 5`) |
| D3 | `H_M`, J = 8 | [NE + M] | [NE±] on the λ-level truncation (`L ≥ 4.3`); sign-free `L ≥ 1.1` |
| D3 | `H_M`, J = 12 | [NE + M] | [NE + M], now on an own-mass row (no 1/m interpolation); λ-level `L ≥ 1.6` |
| D3 | `H_M`, J = 15, 16 | [NE + M] | [NE + M] on own-mass rows (λ-level `L = 1.0–2.0`) |
| D8 | value below ν ≈ 8.3 at J = 5 | [M] | inside the sign-free `E_low` (engine-supported) |
| D9 | odd J, local de-shift only | [M], `L ≥ 100` | at J = 5: inside `E_low` via the min-norm gauge factor |
| D6 | J = 15 twist edge | [M], hybrid ≈ +0.69 | **computed `+0.690`** (own row) / `+0.675` (published rows) [NE + M] |
| D4 | `H_edge` | [M] on even J = 18–32, odd 19–29, J = 48 | **J = 18, 19, 24, 25: [NE + M]** (shifts −0.007 … +0.032); **J = 29: [M], representation-dominated** (§1.6); J = 20–22, 26–28, 30–32, 21, 23, 27, 48: [M] (not run) |
| — | binding state at Δ_g = 300 | [NE + M] (`+0.2038`) | **[NE±]**, own row `Δ = 305`, `L = 3.5` (`≥ +0.146`) |
| A1PRIME §6 | AdS light J ≤ 2 just below the gap | [I] | **[NE±]** at `Δ_L = 0.99Δ_g` (own rows, `L ≥ 913`); [NE + M] for `Δ_L ≳ 0.7Δ_g` (1000) / `≳ 0.8Δ_g` (300) (§2.3) |
| A1PRIME §6 | AdS light J ≤ 2 at `Δ_L ≪ Λ`, Ω^EE | [D-sketch] | unchanged as a statement about the true image; the LR image is positive with `F/N ≥ 0.64` ([P given LR] for `Δ_L ≥ 4.5`) |
| A1PRIME §4.3 | repaired Q on light J = 0 | "positive on all light states" [P, flat] | flat margin 0 at `m → 0`; `Δ_L = O(1)`: **[I]**; use `ε = 0.040` |

**Letter wording (suggested).** "At Δ_g = 1000 the lowest state (J = 5, parity odd, at the twist edge) has λ/N_s = 0.201; on a crystal row
at its own mass, with the Laplace truncation, the small-ν and light-cone windows bounded sign-free with factor-3 or larger safety
margins, and checked against direct Monte-Carlo integration, λ/N_s ≥ 0.198." Do not extend "error-bounded" beyond J ≤ 8.

## 4. Caveats
1. "K₂ = crystal integral" is checked, not derived: direct RQMC at 9 (J = 5), 5 (J = 6) and 4 (J = 8) nodes, and the analytic-radial
   engine at 6 (J = 5). The own rows at J = 4, 12, 15, 16, 0, 2 and Δ_g = 300 have no direct check here.
2. The ×3 factors are choices. At J = 5 the observed convergence (×0.1 per order) makes the effective factor ≈ 30.
3. `G_low = G_edge = 10⁴ c₂` are hypotheses with large measured headroom at J = 5 (engine, direct), not derived bounds.
4. Only masses up to 1.05Λ were budgeted; beyond, the states sit between the own row and the 1600 row (1/m interpolation, as before).
5. The light-state analysis uses the heavy-state parametrisation `m² = α(α + 2J + 1)` and the leading-Regge form at small Δ_L (§2.2).
6. The `c₁²` coefficient of the ε = 0.040 upper end is not recomputed.

## 5. Reproduce
```bash
cd ads4/code/neargap2        # every command under: ulimit -v 3000000; OMP_NUM_THREADS=1
python3 ng_rows.py odd EE,C3H 5 1005 dense                 # own row, binding state (~17 min)       -> cache/w1_p3_*_J5_ng.pkl
python3 ng_rows.py even EE,C3H 4 1004 dense                # likewise J = 6, 8, 12, 16 (1000 + J); 0, 2 at 990
python3 ng_rows.py odd EE,C3H 15 1015 <x list>             # J = 15 twist edge
python3 ng_points.py 5 1005 0.8,0.825,0.85,0.875,0.9,0.1,0.3,0.5,0.7     # direct n1 nodes        -> cache/hc_pts_*_J5_D1005.pkl
python3 ng_engine5.py 20000 3 4; python3 ng_engine5.py 1005 3 4          # parity-odd engine EE2  -> cache/ng_po_eeq2_*.pkl
python3 ng_rows.py odd EE 5 1005 0.0125,0.015,0.0175,0.02,0.0375; python3 ng_engcmp.py
python3 ng_bind.py 1000 J5 5:1005:1.0049927,1.01,1.02,1.05             # budgets -> cache/ng_bind_<tag>.{txt,json}
python3 ng_val.py 1000 5 1005 1.0049927,1.02                             # direct vs K2 -> cache/ng_val_J5_1005.txt
python3 ng_rowdump.py 5 1005                                             # Gamma per order
python3 ng_light.py                                                      # light states -> cache/ng_light.{txt,json}
python3 ng_points.py 19 1019 0.8,0.825,0.85,0.875,0.9                    # edge nodes (also 25/1025, 29/1029, 18/800+1600, 24/500+1600)
python3 ng_val.py 1000 18 800,1600 1.0178588,1.03,1.05                    # edge injection -> cache/ng_val_J18_800_1600.txt
python3 ng_edgew.py 1000 19:1.0188417:5650 ...                           # sign-free weight of the edge windows
python3 ng_summary.py                                                    # every budget and state shift -> cache/ng_summary.txt
python3 ng_bind.py 1000 light1000 0:400:<mfl list> 2:400:<mfl list>      # and 300 light300 ...
python3 ng_bind.py 1000 light990 0:990:0.989985 2:990:0.989982          # light states on their own rows
python3 ng_rows.py odd EE,C3H 5 305 dense; python3 ng_bind.py 300 J5_300 5:305:1.0165851,1.02,1.05   # binding state at 300
```
The exact job lists, in run order, are in `cache/lane_A.log` and `cache/lane_B.log`.
