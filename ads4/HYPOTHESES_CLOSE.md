# Closing the computational hypotheses of Result D (`H_edge`, `H_OO`, `Δ_g = 300` small-`ν`, `K₂` beyond its criterion)

*20–21 Nov 2026, hypothesis-closing agent. Unrefereed, not committed. No existing document or code file was edited. New code
`ads4/code/spinfunc/hc_*.py` (+ `hc_queue*.sh`), caches `spinfunc/cache/hc_*` (one file per state; no shared reference keys).
The functional `Ω^EE` (min-norm gauge, NF02 weights) and the analysis protocol (sm_chainfix chains: curvature double count removed,
`SM_CURV=row`, κ image; Jacobi fix; Hampel despike, odd rows included as in `st_odd_despike.py`) are FIXED; nothing was tuned on
the signs found here.*

Grades: **[P]** proved, **[D]** derived, **[NE]** numerical evidence, **[M]** model, **[I]** open.

**STATUS: in progress** (see hand-over at the end).

## 1. Methods and validation

### 1.1 Direct points (`hc_points.py`, `hc_nodes.py`)
* **Edge and bulk nodes (`HC_METHOD=n1`).** The n1 real re-centred contour (`n1_run.point` integrand for even `J`; `p8_odd_pw`'s
  point-wise `p3_odd_action` for odd `J`, kind `EE` = `w0.make_build_odd_EE`), with
  * the `Q_J` overflow rescaling of z2 (`set_qsc` at the n1 saddle; without it every heavy-reference point at `ν ≲ 650` is masked,
    which is why `ee1000_direct` had to use `K₂` references),
  * the scalar `Q_J` by the Legendre recurrence and the spin-`J` Fourier floor `fw_frame.patch_S` (FOLD_WINDOW §1),
  * a moment-matched Gaussian proposal (`n1_direct.adapt`, pilot seed 5) and RQMC seeds `2000 + r`, `r = 0…7`, **shared by the state and
    its reference** (`Δ_ref = 2·10⁴`, same `ν`). At `(1040, 40)`, `x = 0.9` the adaptive proposal lowers the element errors ×8 at equal
    cost (fixed `s = 1.2`: 28–40 %; adaptive: 4–5 % at `2⁹·2·8` points; `0.9–1.4 %` on the helicity diagonals at `2¹³·2·8`).
  * `fw_frame.patch_conserved_clean` is **not** applied (needed only at `J ≥ 150`): it rotates the SVD basis of the conserved `⟨TTO⟩`
    null space, and the points must be in the structure basis of the `K₂` rows they are injected into.
* **Small-`ν` nodes (`HC_METHOD=n3`).** z2's n3 hybrid contour (`ν_f = 13` below 13), proposal adapted on the state and re-used for the
  reference (common random numbers), `Z2_NUSPLIT=5` (n3 also at `ν = 8.3, 10`; no thimbles), `2¹⁰·2·8` points.
* **Records.** `M0 = E b⁻ᵀ` (even) / `E/b_OO` (odd), `M2 = M4 = 0`, `logC = logN`, `logK` of the shadow; the shadow `b` is the one whose
  stored fit residual equals the `shadow_resid` of the row's `K₂` records (fw_fix's rule), so direct and `K₂` records share the
  structure/shadow basis. One cache file per state (`cache/hc_pts_<kind>_J<J>_D<Δ>.pkl`), holding the state and its own references.
* **Validation (state values).** `(1040, 40)`, `x = 0.8`, EE: the adaptive point agrees element by element with the independent
  `ee1000_n1` point (fixed proposal, `2¹⁴` points): ratios `1.0079 / 0.9977 / 0.9962 / 0.9993` against errors `0.64 % / 0.18–0.5 %`
  (≤ 1.2σ), identical `logN` [NE]. So the hc points are in the same structure basis and normalisation as the project's direct points.

### 1.2 Chains (`hc_chain.py`)
* `sm_chainfix`'s modes re-run in-process (`ee1000`, `w2`, `z4`), curvature double count removed (`SM_CURV=row`, κ image), Jacobi fix,
  odd-row Hampel despike as `st_odd_despike.py`. The hc records are injected into the record loaders of the named rows (both kinds,
  state and reference); a `K₂` node at the same `ν` is replaced; the helicity basis is computed **without** the injected records.
* **Reproduction:** with nothing injected, `ee1000` (`J = 36, 40, 17, 31`, `XMAX = 0.8`) gives `+0.9669 / +0.9794 / +0.8389 / +0.9945`,
  the statement's numbers exactly; `z4 300` (`J = 24, 28`, rows 336/500/1600, `Z4_EST=best`) gives `+0.8204 / +0.8133` (review §2).
* **Edge representation (ee1000 mode, `HC_EDGE_MODEL=piece`).** On `x ≤ 0.8` the protocol's model is kept (shape fit `K = 4` on the
  row's nodes `x ≤ 0.8`); on `x > 0.8` a separate shape fit `Γ = c D_h + s xD_h′/(2J+1) + r` with constant (`HC_EK=0`) or linear
  coefficients on the nodes `x ≥ 0.8`, held beyond the last node. Variants (`hc_edge.py`): `A` = direct only for `x > 0.8` (the `K₂`
  node at 0.8 kept), `B` = direct from 0.8 on (main), `B_lin` / `B_spl` = linear / natural-spline interpolation of the edge nodes,
  `B_nd` = injected nodes exempt from the Hampel filter, `glob` = one `K = 4` shape fit over `[0.05, 0.9]`; jackknife over the 8
  replicas of the injected nodes. The `+−` window weight sits at `x ≤ 0.875` (`J = 40` gap: 61 / 32 / 11 % of the window in
  `(0.8, 0.825] / (0.825, 0.85] / (0.85, 0.875]`, the rest `≤ 3 %`, with sign changes from `ψ₂`), so nodes `0.8 … 0.9` cover it.

## 2. `H_edge`: the light-cone-edge window `x ∈ (0.8, Λ/m]` with direct values

### 2.1 Node values (`Γ`, `c₂` units, curvature image subtracted; `hc_chain` node report)
Direct edge nodes oscillate with the `D_h` period and **grow** towards the edge, much faster than the protocol's held extrapolation
(constant `c, s, r` at `x = 0.8`):

| state | key | `x = 0.8` (direct / `K₂`) | 0.825 | 0.85 | 0.875 | 0.9 | held model at 0.825 … 0.9 |
|---|---|---|---|---|---|---|---|
| `(1040, 40)` | C2 `+−` | `+4379 / +1367` | `−6117` | `+11509` | `−17578` | (despiked) | `−2834, +3984, −4221` |
| `(1040, 40)` | C2 `++` | `+5864 / +5475` | `−7021` | `+4622` | `+6494` | — | `−5773, +3882, +1837` |
| `(1036, 36)` | C2 `+−` | `−2596 / −3853` | `+3860` | `−3981` | `+4736` | `+13136` | `+3769, −3596, +1635, +2202` |
| `(1036, 36)` | C2 `++` | `+2486 / +594` | `−3442` | `+6363` | `−12634` | `+14615` | `−907, +2034, −4484, +5875` |
| `(1031, 31)` | C2 (`+−`) | `+1792 / +2744` | `−3853` | `+5975` | `−8856` | `+14453` | `−3277, +3334, −3369, +3700` |

(C3H follows C2 with the sign of `σ_h`; `K₂` at `x = 0.825–0.875` exists only on the `J = 40` row and is 15–25 % below direct at
0.85–0.875 and ×3 off at 0.8.) Statistical errors (jackknife over replicas) are `≲ 1 %` of these values; the dominant uncertainty is
the representation between nodes (half-period of the oscillation ≈ node spacing 0.025 at `J ≥ 36`).

### 2.2 Per-state effect, `Δ_g = 1000` (min `λ/N_s`; gap state and minimum over `m/Λ ∈ [m_gap, 5]`)

| `J` | protocol (held) gap / min | B (direct, shape edge) gap / min | other representations, gap | jackknife (B, gap) |
|---|---|---|---|---|
| 40 | `+1.018 / +0.979` | `+1.007 / +0.979` | A `+1.013`, lin `+1.004`, spline `+0.993`, glob `+0.974` | `±0.0014` |
| 36 | `+0.994 / +0.967` | `+1.055 / +0.967` | A `+1.043`, lin `+0.853`, spline `+0.843`, glob `+0.986` | `±0.0032` |
| 31 (odd) | `+1.200 / +0.995` | `+1.037 / +0.988` | A `+1.077`, lin `+1.016`, spline `+1.029`, glob `+1.111` | `±0.0039` |
| 17 (odd) | `+0.992 / +0.839` | `+0.997 / +0.840` | A, lin, spline, glob: `+0.997 … +0.998` | `±0.0008` |

### 2.3 Per-state effect, `Δ_g = 300` (z4 chain, rows 336/500/1600, `Z4_EST=best`, natural spline in `x` as the protocol; edge nodes
exempt from the Hampel filter, `x = 0.9` excluded at `J = 24` (RQMC error > 100 % there: `ν = 300` is too close to `Δ = 336`))

| `J` | variant | gap state (`1.078Λ`) | 1.1Λ | min over `m` |
|---|---|---|---|---|
| 24 | protocol (held) | `+1.021` | `+1.013` | `+0.8204` [1.5] |
| 24 | B (direct 0.8–0.875) | `+1.011 ± 0.0007` (jk) | `+1.009` | `+0.8204` [1.5] |
| 24 | A (K₂ kept at 0.8) / B with Hampel on | `+1.028` / `+1.008` | `+1.016` / `+1.008` | `+0.8204` |
| 24 | linear interpolation in `x` (protocol / B) | `+1.040 / +1.014` | `+1.053 / +1.042` | `+0.8101` [1.5] |

| 16 | protocol (held) | `+0.767` [gap 1.053] | `+0.707` | `+0.7074` [1.1] |
| 16 | B (direct 0.8–0.9) | `+0.750` | `+0.703 ± 0.0001` (jk) | `+0.7032` [1.1] |
| 16 | A / B with Hampel on | `+0.761 / +0.755` | `+0.705 / +0.703` | `+0.7048 / +0.7034` |
| 16 | linear interpolation (protocol / B) | `+0.758 / +0.743` | `+0.702 / +0.698` | `+0.7020 / +0.6982` |

### 2.4 The binding state, `J = 5` (odd, `+−`), both `Δ_g`
* Direct edge nodes `x = 0.8, 0.85, 0.9` on all four `J = 5` rows (`Δ = 200, 400, 800, 1600`; `x = 0.95` has `~100 %` RQMC error on
  the 200/400 rows and is excluded): `|Γ| ≤ 50 c₂` (C2 `+43, +17, +3` at 1600; `≤ 6` at 200/400), against a failure threshold for a
  uniform edge value of `≈ −6·10⁴ c₂` (300) / `−5·10⁵ c₂` (1000).
* w2 chain, all `J = 5` states: **`+0.2052` (300, 1.02Λ) and `+0.2021` (1000, 1.01Λ), unchanged to `< 5·10⁻⁵`** in every variant
  (protocol, direct, direct + Hampel, linear interpolation).

Edge node values at `(336, 24)` (C2 `++ / +−`): `x = 0.8`: `−279 / +1935` (K₂ `+1316 / +1282`), 0.825: `−1844 / −953`, 0.85:
`+3418 / −2489`, 0.875: `−3566 / +5769`; the held protocol value is `+1316 / +1282` throughout. At `(400, 16)`: `x = 0.8 … 0.9`:
`++` `+556, +878, −818, −1555, +2580`, `+−` `−1009, −559, +1286, +2216, −899` (held `+244 / −891`).

* **`J = 31`**: the favourable `+0.21` gap correction the review flagged (S3) came from the held `x = 0.775–0.8` nodes; with direct
  edge values it is `+0.03 … +0.12`. The gap state stays at `≥ +1.016` in every representation.
* **`J = 36`**: the representation spread is large (`+0.84 … +1.06`) because the `+−` edge values alternate in sign from node to node
  (period ≈ 0.05 ≈ 2 × spacing); linear/spline interpolation is at the Nyquist limit there. Finer nodes (`0.8125, 0.8375, 0.8625`) are
  queued (`hc_queueA2.sh`). Every representation is positive (worst `+0.84`).

## 3. `H_OO`: the even-`J` parity-odd block beyond leading order

### 3.1 What it is and which states it affects
* At even `J` the conserved `⟨TT𝒪_{Δ,J}⟩` space is [parity-even (2 at `J ≥ 4`: `++`, `+−` helicity combinations), parity-odd (1: the
  symmetric odd structure = CHL's `X^O_same`)] (PARITY_ODD_ROUTEA §1: one symmetric odd structure at `J = 0, 2, 4, 5, 6, …`, none at
  `J = 1, 3`). The functional's `S⁰` sum cancels the even–odd interference identically (§3 [DV]), so every even-`J` heavy block is
  `EE (2×2) ⊕ OO (1×1)`. **The OO block is a separate per-state condition for every even-`J` state, at every mass** (near gap and
  beyond, `J = 0, 2` included); at odd `J` the single structure is parity-odd (`X^O_opp`, `+−`) — that is the already computed odd
  sector, not the OO block.
* At leading order `OO = N_same × (++ density)`, `N = tan²(πτ/2) ≥ 0` (§4/§12), so `λ_OO/N_s = (base_++ + dF_OO)/N_s,++`, with
  `dF_OO` the beyond-LO correction of the odd structure (never computed before this note).

### 3.2 Rows and LO check (`hc_oo_rows.py`, `hc_oo.py`)
* `K₂` rows of the 3×3 `p3_odd_action` action with `w0`'s EE integrand (kind `EE`) and `o_H q` (C3H), rows `Δ = 400, 1600` (the EE
  rows of `J ≤ 10`), 19 nodes `x = 0.05 … 0.8` at `J = 4, 6`, 11 nodes at `J = 8, 10`, reference `2·10⁴`; OO = entry `[−1, −1]`.
* `Γ_OO` in the l9/p4 conventions with the `++` shapes (`a_O` per row from `x ≤ 0.25`, `ν ≥ 15`; spinning de-shift `(1 ± u_EE)/2`;
  `C3H` no de-shift; curvature image subtracted as for `++`).
* **LO check [NE]:** the order-0 ratio `C3H_OO/C2_OO × (2 − y)` is `x`-independent to `≤ 0.3 %` on every row (`1.986` at `Δ = 400`,
  `1.996` at 1600), i.e. the odd structure couples to `Ĉ₃` exactly like `++` (`env₃ = env₂/(2 − y)`, `σ = +1`).
* `Γ_OO` is of the same size as the EE `++` correction: `|Γ_OO| ≤ 70 / 150 / 270 / 430 c₂` (`J = 4 / 6 / 8 / 10`, row 400, both kinds).

### 3.3 Folded into the chain (EE_NEARGAP main chain, `SM_CURV=row`, all exact-skeleton + hybrid states `m/Λ ≤ 5`)

| `Δ_g` | `J` | min `λ_OO/N_s` [`m/Λ`] | LO (uncorrected) | max `|dF_OO|/N_++` | max `|dF_++|/N_++` (EE) | EE-block min `λ` |
|---|---|---|---|---|---|---|
| 300 | 4 | **+0.2301** [1.02] | +0.2308 | 0.0008 | 0.0006 | +0.2204 |
| 300 | 6 | **+0.3071** [1.25] | +0.3091 | 0.0044 | 0.0039 | +0.3073 |
| 300 | 8 | **+0.4774** [1.65] | +0.4797 | 0.0164 | 0.0142 | +0.4293 |
| 300 | 10 | **+0.6022** [2.0] | +0.6046 | 0.0454 | 0.0371 | +0.6026 |
| 1000 | 4 | **+0.2150** [1.01] | +0.2150 | `< 10⁻⁴` | `< 10⁻⁴` | +0.2103 |
| 1000 | 6 | **+0.3093** [1.25] | +0.3093 | 0.0002 | 0.0001 | +0.2508 |
| 1000 | 8 | **+0.5032** [1.75] | +0.5033 | 0.0005 | 0.0004 | +0.4344 |
| 1000 | 10 | **+0.6046** [2.0] | +0.6047 | 0.0015 | 0.0015 | +0.5853 |

* Near the gap `dF_OO` has the sign and size of the EE `++` correction (ratio `1.0–1.25`; e.g. 300, `J = 10`, 1.05Λ: `−0.0438` vs
  `−0.034`); `λ_OO` sits within `≤ 0.008` of the EE `++` diagonal at every state.
* **Grade:** OO block positive with margin `≥ 0.21` on every computed state, `J = 4–10`, `Δ_g = 300, 1000`, `m/Λ` from the gap to 5
  **[NE + M]** (M = the same K₂ chain model as the EE block). For `J ≥ 12` the OO block is **[M]**: it tracks `base_++ + dF_++` with a
  correction ratio `≈ 1–1.25` at `J ≤ 10`; where `dF_++` is large (300, `J ≥ 18`: `|dF_++|/N ≈ 0.5–1` in convention A) a 25 %
  difference is `O(0.1–0.2) N_s` — a `J = 18` OO row (`Δ = 336, 400`) is queued as a comparative test (`hc_queueOO2.sh`).

## Hand-over (kept current)

* **Queues** (`ads4/code/spinfunc/`, all resumable: re-launch with `setsid nohup ./hc_queueX.sh > cache/hc_queueX.log 2>&1 &`;
  finished points are skipped; at most 2 workers run at a time):
  * `hc_queueA.sh` (edge nodes, done), `hc_queueOO.sh` (OO rows J = 4–10, done);
  * `hc_queueA2.sh`: finer edge nodes on `(1040, 40)`, `(1036, 36)`, `(1031, 31)`; bulk nodes `(336, 28)`, `(336, 32)`; small-`ν`
    nodes `(400, 18)`; spot nodes `(400, 4)`, `(400, 5)`;
  * `hc_queueB.sh`: small-`ν` nodes `(336, 22)`, `(400, 20)`, `(336, 18)`;
  * `hc_queueOO2.sh` (after B): OO rows `J = 18`, `Δ = 336, 400`.
* **Code:** `hc_points.py` (direct points; `HC_METHOD=n1|n3`), `hc_nodes.py` (records, shadow matching, injection), `hc_chain.py`
  (chains with injection; modes `ee1000`, `w2`, `z4`; node report), `hc_edge.py` (edge variant table, ee1000), `hc_oo_rows.py`,
  `hc_oo.py` (OO block), `hc_sum.py`, `hc_jk.py` (summaries / jackknife of chain outputs).
* **Caches:** `cache/hc_pts_<kind>_J<J>_D<Δ>.pkl` (+ `.lock`), `cache/hc_raw_J<J>_D<Δ>.pkl`, `cache/hc_p3_<kind>_J<J>.pkl`,
  `cache/hc_p2shadow_<kind>_J<J>.pkl`, `cache/hc_exact.pkl` (z3's exact scalar records, redirected), `cache/hc_chain_*`, `cache/hc_edge_J*.txt`,
  `cache/hc_oo_*.txt/.pkl`, logs `cache/hc_queue*.log`.
* **Reproduce:** edge (1000): `python3 hc_edge.py 40 1040 1000` (likewise 36/1036, 31/1031, 17/1017); edge (300, `J = 24`):
  `HC_TAG=t HC_INJ=24:336 HC_IXMAX=0.88 HC_DESPIKE_INJ=0 Z4_EST=best Z4_JS=24 Q3_JEVEN=24 Q3_JODD=5 Q3_EXTRA_M=1.078,1.1,… W2_ROWS=336,500,1600 python3 hc_chain.py z4 300`
  then `python3 hc_sum.py z4 300 t J=24`; OO: `HC_OOJ=4,6,8,10 Q3_EXTRA_M=1.02,…,5 python3 hc_oo.py 300 1000`.
* Analysis runs were made with one production worker SIGSTOPped (scratch `paused.sh`), to stay at 2 busy processes.
