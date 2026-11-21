# Closing the computational hypotheses of Result D (`H_edge`, `H_OO`, `Δ_g = 300` small-`ν`, `K₂` beyond its criterion)

*20–21 Nov 2026, hypothesis-closing agent. Unrefereed, not committed. No existing document or code file was edited. New code
`ads4/code/spinfunc/hc_*.py` (+ `hc_queue*.sh`), caches `spinfunc/cache/hc_*` (one file per state; no shared reference keys).
The functional `Ω^EE` (min-norm gauge, NF02 weights) and the analysis protocol (sm_chainfix chains: curvature double count removed,
`SM_CURV=row`, κ image; Jacobi fix; Hampel despike, odd rows included as in `st_odd_despike.py`) are FIXED; nothing was tuned on
the signs found here.*

Grades: **[P]** proved, **[D]** derived, **[NE]** numerical evidence, **[M]** model, **[I]** open.

**STATUS: all queued computations finished (21 Nov).** Results per hypothesis (details and variants in §2–§5):

| hypothesis | what was computed | result | grade now |
|---|---|---|---|
| **`H_edge`** (edge window `x ∈ (0.8, x_edge]` held) | direct n1 nodes `x = 0.8–0.9` (8 nodes, step 0.0125, at `J = 31, 36, 40`; 5 at `J = 17, 16, 24`; 3–4 on the four `J = 5` rows) on the rows carrying the gap states: 1000 `J = 5, 17, 31, 36, 40`; 300 `J = 5, 16, 24`; chains re-run | every state passes; gap-state change `+0.005 … −0.14 N_s` (K1 edge model), `J = 31`'s favourable `+0.21` held artefact reduced to `+0.07`; minima `+0.9795 / +0.957 / +0.988 / +0.840 / +0.2021` (1000: `J = 40/36/31/17/5`), `+0.820 / +0.703 / +0.2052` (300: `J = 24/16/5`) | **closed [NE + M]** on the computed rows (M = edge representation between nodes; all representations positive, worst `+0.86`); **[M]** (assumption, supported) on the other rows |
| **`H_OO`** (even-`J` parity-odd block) | identified (`X^O_same`, the 1×1 block of every even-`J` state, `J = 0, 2, 4, …`); K₂ OO rows `J = 4, 6, 8, 10` (Δ = 400, 1600) and `J = 18` (Δ = 336, 400); folded into the EE_NEARGAP chain | LO coupling `= ++` confirmed (`C3H/C2·(2−y)` constant to 0.3 %); min `λ_OO/N_s` `+0.230 / +0.307 / +0.477 / +0.602` (300) and `+0.215 / +0.309 / +0.503 / +0.605` (1000) at `J = 4/6/8/10`; `J = 18` (300): `+0.222` (superseded small-`ν` treatment), OO correction `1.6×` the EE `++` one | **closed [NE + M]** for `J = 4–10` (both `Δ_g`); `J = 18` positive in the superseded treatment only; **`J ≥ 12` [M]** — "OO tracks `++`" fails quantitatively at `J = 18` |
| **`Δ_g = 300`, `J ≤ 22`, validated small-`ν`** | z4 protocol with direct n3 nodes `ν = 8.3–40` on the rows of `J = 18` (336, 400), 20 (400), 22 (336); spot nodes `J = 4, 5`; full sweep incl. every `J`'s twist-edge mass | min `λ/N_s`: `J = 18 / 20 / 22` **`+0.872 / +0.863 / +0.841`** (convention A had `+0.472 / +0.465 / +0.046`); `J = 4 / 5` `+0.2137 / +0.2038` (unchanged); jackknife `±0.03–0.08`; failure only for a uniform `Γ_c ≲ −1·10⁴ … −2·10⁴ c₂` below `ν = 8.3` | **[NE + M]** for `J = 4, 5, 18, 20, 22`; other `J ≤ 16` (even) and odd `J ≤ 15` still convention A (positive, min `+0.2038`), which was pessimistic wherever tested |
| **`K₂` beyond criterion, `(336, 28)`, `(336, 32)`** | complete direct n1 sub-rows `x = 0.125–0.8` (16 nodes each); z4 chains | `K₂` wrong by up to ×5 in amplitude and in phase from `x = 0.2`; with direct bulk `J = 28` `+0.835` (`+0.794` with the repeat thimbles), `J = 32` `+0.783` (its twist-edge state `+0.284 → +0.826`) | bulk **[NE + M]**; `ν < 8.3` window **[I]** (failure at uniform `Γ_c ≈ −5400 / −4250 c₂`) |

Nothing here is a proof: every "closed" means "computed in the project's chain model with direct crystal values replacing the
extrapolation/expansion in question" ([NE + M]).

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
  node at 0.8 kept), `B` = direct from 0.8 on, `B_K1` = B with linear-in-`x` edge coefficients (main once 8 edge nodes exist),
  `B_lin` / `B_spl` = linear / natural-spline interpolation of the edge nodes, `B_desp` = injected nodes subject to the Hampel filter
  (all other variants exempt the direct nodes), `glob` = one `K = 4` shape fit over `[0.05, 0.9]`; jackknife over the 8
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
| 40 (5 edge nodes) | `+1.018 / +0.979` | `+1.007 / +0.979` | A `+1.013`, lin `+1.004`, spline `+0.993`, glob `+0.974` | `±0.0014` |
| **40 (8 edge nodes, `Δx = 0.0125`)** | `+1.018 / +0.979` | **K1 `+1.006 / +0.9795`** (edge-fit residual 2–10 %) | B (K0, residual 26–41 %) `+0.995`, A `+1.001`, lin `+1.005`, spline `+1.003`, Hampel on `+1.008`, glob `+0.915` (min over `m` `+0.908`) | `±0.0013` |
| 36 (5 edge nodes) | `+0.994 / +0.967` | `+1.055 / +0.967` | A `+1.043`, lin `+0.853`, spline `+0.843`, glob `+0.986` | `±0.0032` |
| **36 (8 edge nodes)** | `+0.994 / +0.967` | **K1 `+0.957 / +0.957`** [gap] (edge-fit residual 1–9 %) | B (K0, residual 25–35 %) `+1.036`, A `+1.022`, lin `+0.858`, spline `+0.876`, Hampel on `+1.007`, glob `+0.986` | `±0.0032` (B) |
| 31 (odd, 5 edge nodes) | `+1.200 / +0.995` | `+1.037 / +0.988` | A `+1.077`, lin `+1.016`, spline `+1.029`, glob `+1.111` | `±0.0039` |
| **31 (odd, 8 edge nodes)** | `+1.200 / +0.995` | **K1 `+1.056 / +0.988`** (residual 8 %) | B (K0, residual 39 %) `+1.051`, A `+1.090`, lin `+1.000`, spline `+1.014`, glob `+1.193` | `±0.0041` (B) |
| 17 (odd) | `+0.992 / +0.839` | `+0.997 / +0.840` | A, lin, spline, glob: `+0.997 … +0.998` | `±0.0008` |

* **Main representation = K1** (edge shape fit with linear-in-`x` coefficients on 8 nodes `x = 0.8 … 0.9`, spacing 0.0125): it is
  the only edge model whose residual (1–10 %) is at the level of the in-row fits; K0 (constant coefficients) misses the growth of the
  amplitude towards the edge (residual 25–40 %), and linear/spline interpolation alias the `D_h` oscillation (half-period ≈ 0.025).
  Jackknife of K1 at the gap: `±0.0012 / ±0.0035 / ±0.0042` (`J = 40 / 36 / 31`).
* **`J = 31`**: the favourable `+0.21` gap correction the review flagged (S3) came from the held `x = 0.775–0.8` nodes; with direct
  edge values it is `+0.065` (K1; `+0.01 … +0.10` over the non-global representations). The gap state stays `≥ +1.000`.
* **`J = 36`**: the gap state drops from `+0.994` to `+0.957` (K1) and becomes the minimum over `m`; linear/spline interpolation give
  `+0.86 / +0.88` (aliased); every representation is positive.
* **`J = 40`**: `−0.012` at the gap (K1); the minimum (`+0.9795`, 5Λ) is not in the edge window.
* **`J = 17`**: `+0.005`; edge values `|Γ| ≤ 2800 c₂`, failure threshold for a uniform edge value `≈ −5600 c₂`.
* Changing `x = 0.8` from `K₂` to direct also moves states without an edge window (`J = 31`, 1.4Λ: `+0.9945 → +0.9877`): `K₂` at
  `x = 0.8` is off by `≈ 1000 c₂` (`J = 31`: `+2744` vs direct `+1792`; `J = 40` `+−`: `+1367` vs `+4379`) — an `H_M` (K₂ accuracy) item,
  small here.

### 2.3 Per-state effect, `Δ_g = 300` (protocol chains: z4 for `J = 24` (rows 336/500/1600, `Z4_EST=best`), w2 for `J = 16`; natural
spline in `x` as the protocol; edge nodes exempt from the Hampel filter; `x = 0.9` excluded at `J = 24` (RQMC error > 100 % there:
`ν = 300` is too close to `Δ = 336`))

| `J` | variant | gap state | 1.1Λ | min over `m` |
|---|---|---|---|---|
| 24 | protocol (held) | `+1.021` [1.078] | `+1.013` | `+0.8204` [1.5] |
| 24 | B (direct 0.8–0.875) | `+1.011 ± 0.0007` (jk) | `+1.009` | `+0.8204` [1.5] |
| 24 | A (K₂ kept at 0.8) / B with Hampel on | `+1.028` / `+1.008` | `+1.016` / `+1.008` | `+0.8204` |
| 24 | linear interpolation in `x` (protocol / B) | `+1.040 / +1.014` | `+1.053 / +1.042` | `+0.8101` [1.5] |
| 16 | protocol (held) | `+0.767` [1.053] | `+0.707` | `+0.7074` [1.1] |
| 16 | B (direct 0.8–0.9) | `+0.750` | `+0.703 ± 0.0001` (jk) | `+0.7032` [1.1] |
| 16 | A / B with Hampel on | `+0.761 / +0.755` | `+0.705 / +0.703` | `+0.7048 / +0.7034` |
| 16 | linear interpolation (protocol / B) | `+0.758 / +0.743` | `+0.702 / +0.698` | `+0.7020 / +0.6982` |

Edge node values at `(336, 24)` (C2 `++ / +−`): `x = 0.8`: `−279 / +1935` (K₂ `+1316 / +1282`), 0.825: `−1844 / −953`, 0.85:
`+3418 / −2489`, 0.875: `−3566 / +5769`; the held protocol value is `+1316 / +1282` throughout. At `(400, 16)`: `x = 0.8 … 0.9`:
`++` `+556, +878, −818, −1555, +2580`, `+−` `−1009, −559, +1286, +2216, −899` (held `+244 / −891`).

### 2.4 The binding state, `J = 5` (odd, `+−`), both `Δ_g`
* Direct edge nodes `x = 0.8, 0.85, 0.9` on all four `J = 5` rows (`Δ = 200, 400, 800, 1600`; `x = 0.95` has `~100 %` RQMC error on
  the 200/400 rows and is excluded): `|Γ| ≤ 50 c₂` (C2 `+43, +17, +3` at 1600; `≤ 6` at 200/400), against a failure threshold for a
  uniform edge value of `≈ −6·10⁴ c₂` (300) / `−5·10⁵ c₂` (1000).
* w2 chain, all `J = 5` states: **`+0.2052` (300, 1.02Λ) and `+0.2021` (1000, 1.01Λ), unchanged to `< 5·10⁻⁵`** in every variant
  (protocol, direct, direct + Hampel, linear interpolation).

### 2.5 Grade of `H_edge`
* **Computed states (the binding and representative ones):** `Δ_g = 1000`: `J = 5` (odd, binding), 17 (odd), 31 (odd), 36, 40 at the
  gap mass and every mass up to 5Λ; `Δ_g = 300`: `J = 5`, 16, 24. The held extrapolation on `x ∈ (0.8, x_edge]` is replaced by direct
  values (RQMC errors `0.1–1 %` per element; jackknife on `λ/N_s` `≤ 0.004`) and the states are re-run in the protocol's chains.
  **Every one of them passes**; change of the gap state `+0.005 … −0.14 N_s` (K1 / protocol spline), minima after the change
  `+0.9795 / +0.957 / +0.988 / +0.840 / +0.2021` (1000: `J = 40 / 36 / 31 / 17 / 5`) and `+0.820 / +0.703 / +0.2052` (300: `J = 24 / 16 / 5`).
  Grade **[NE + M]** — NE for the node values, M for the representation of the edge oscillation between nodes (the K1 edge fit; the
  alternative representations move the gap state by up to `0.1 N_s` at `J = 36` and stay positive, worst `+0.86`; the one-polynomial
  `glob` fit, which the K1 residuals show to be inadequate there, gives `+0.915` at `J = 40`).
* The review's truncation estimate (`x ≤ 0.6` cut: `−0.43 / −0.33 / −0.23` at `J = 36 / 40 / 31`) over-states the effect: with direct
  edge data the change is `−0.04 / −0.01 / −0.14`. The one sizeable shift (`J = 31`) removes the favourable held artefact flagged in
  review S3 rather than revealing a deficit.
* **Not computed:** the other rows (1000: even `J ≤ 32`, 48, odd `J = 7–15, 19–29`; 300: even `J ≤ 22` except 16, `J = 28, 32`, odd
  `J = 7–15`) still use the held extrapolation. There `H_edge` remains an assumption, now supported by the computed rows ([M]): the
  edge values grow towards the edge like the `+−` envelope, but the integrated effect on the computed states is `≤ 0.04 N_s` except
  where a favourable held artefact was removed (`J = 31`, `−0.14`), against margins `≥ 0.63` (`≥ 0.20` at `J ≤ 9`, where the edge
  contribution is `< 10⁻³ N_s`).

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
* **`J = 18` at 300** (OO rows `Δ = 336, 400`, 11 nodes; the comparative test). LO check `C3H/C2 × (2−y) = 1.99–2.00` (one node
  2.05–2.11 at `x = 0.15`). In the same `K₂` protocol as the EE block (convention A, local; `K₂` + spline below `x = 0.05`):
  `Γ_OO` at `x = 0.05` is `−902 / −658 c₂` (rows 336 / 400) — the same small-`x` `K₂` pattern that was an artefact for the EE `++`
  (`K₂ −1331` vs direct `+212` at `(336, 24)`) — and **min `λ_OO/N_s = +0.222`** at the gap `1.0586Λ` (`dF_OO/N = −0.78` vs EE `++`
  `−0.49`; EE-block min `+0.472`); with min-norm weights on the `K₂` nodes (`W2_GAUGE=min`): `+0.239`. So at `J = 18` the OO
  correction is `≈ 1.6×` the EE `++` one, and the OO block, not the EE block, is the binding block of `J = 18` in this treatment.
* **Grade:** OO block positive with margin `≥ 0.21` on every computed state, `J = 4–10`, `Δ_g = 300, 1000`, `m/Λ` from the gap to 5
  **[NE + M]** (M = the same K₂ chain model as the EE block); `J = 18` at 300: positive (`+0.22 / +0.24`) **in the superseded
  small-`ν` treatment only** [NE + M(superseded)]; the validated treatment would need direct small-`ν` nodes of the odd structure at
  even `J` (the hc n3 path supports the 3×3 `p8` integrand, not run). For `J ≥ 12` the OO block remains **[M]**: the OO/`++` correction
  ratio grows from `1.0–1.25` (`J ≤ 10`) to `1.6` (`J = 18`), so "OO tracks `++`" is **not** a safe extrapolation where `dF_++` is
  large (300, `J ≳ 18`; 1000 large `J`). `H_OO` is narrowed, not closed.

## 4. `Δ_g = 300`, `J ≤ 22` with the validated small-`ν` treatment (z4 protocol)

Method: z4_chain (EE_DG300_SMALLNU §6) with the hc n3 nodes (`ν ≈ 8.3, 10, 11.7, 13.3, 16.7, 20.9, 25, 33.4`, `2¹⁰·2·8` points, CRN,
state and reference in one per-state file) injected on the rows that carry the near-gap states, non-local min-norm `X` below
`ν = 25`, local de-shift above; curvature fix; Jacobi fix; estimator `best` (per node and kind, scalar control variate where it lowers
the jackknife error); below the lowest node (`ν = 8.34`) `Γ_NL` held (variants: zero, linear, constant stress); jackknife over the 8
common-seed replicas (`Z4_JK`). States: all exact-skeleton + hybrid states incl. the twist-edge mass (review m2).

### 4.1 `J = 22` (row `Δ = 336`, the only row; states `m/Λ ≥ 1.072` = gap)

| variant | gap `1.072Λ` | 1.1Λ | min over `m` |
|---|---|---|---|
| convention A (`K₂` + spline below `x = 0.05`, local `ψ₂` weights) — superseded | `+0.046` | `+0.123` | `+0.046` [gap] |
| **direct small-`ν`, `best`, held below 8.3** | **`+0.842 ± 0.084`** (jk) | `+0.852` | **`+0.842`** [gap] |
| `dir` / `cv` estimator | `+0.802 / +0.842` | `+0.819 / +0.852` | gap |
| zero / linear below 8.3 | `+0.897 / +0.745` | `+0.901 / +0.767` | — |
| `Γ_c := −1000 / −2000 / −4000 / −8000` below 8.3 | `+0.81 / +0.72 / +0.54 / +0.19` | | failure at `≈ −1.0·10⁴ c₂` |

* Node values (`Γ_c`, `++`, `c₂`): `ν = 8.34`: `−1230 − 492i` (5 % RQMC error); `10.0`: `−58`; `11.7`: `+87`; `13.3`: `+70`; `16.7`: `+28`;
  `20.9`: `−44`; `25.0`: `−118`; `33.4`: `−109` — `O(10²)` above `ν = 10`, as at `J = 24, 28`.
* The convention-A failure mode is the same `K₂`-at-`x = 0.05` artefact as at `J = 24, 28`; with direct values the `J = 22` margin
  goes from `+0.05` (gap state; `+0.12` at 1.1Λ) to `+0.84`. **New:** the gap state `1.072Λ` (not sampled before, review m2) was the
  worst state in convention A (`+0.046`, i.e. `t* ≈ 1.05`).

### 4.2 `J = 20` (row `Δ = 400`, which carries every near-gap state `m ≤ 1.33Λ`; gap `1.0648Λ`, first sampled mass 1.0656Λ)

| variant | gap state | 1.1Λ | min over `m` |
|---|---|---|---|
| convention A — superseded | `+0.465` | `+0.534` | `+0.465` [gap] |
| **direct small-`ν`, `best`, held below 8.36** | **`+0.866 ± 0.050`** (jk) | `+0.943` | **`+0.866`** [gap] |
| `dir` / `cv` | `+0.856 / +0.865` | `+0.929 / +0.944` | gap |
| zero / linear below 8.36 | `+0.903 / +0.808` | `+0.980 / +0.887` | `+0.866` [4Λ] / gap |
| `Γ_c := −1000 / −2000 / −4000 / −8000` below 8.36 | `+0.84 / +0.78 / +0.65 / +0.40` | | failure at `≈ −1.45·10⁴ c₂` |

* Node values (`Γ_c`, `++`): `ν = 8.36`: `−1229` (5 % RQMC error); `9.95`: `−182`; `11.7`: `−26`; `13.3`: `+15`; `16.7`: `+39`;
  `19.9`: `+14`; `29.8`: `−33`; `39.8`: `−46`. The `ν ≈ 8.3` node is `≈ −1230` at both `J = 20` and 22 (−395 / −895 at `J = 24 / 28`):
  the n3 contour's least reliable node (scalar check at `J = 24`: `−301 ± 1020 c₂`); the zero/const variants bracket it.

### 4.3 `J = 18` (rows `Δ = 336` and `400`, both with direct nodes; gap `1.0585Λ`; states below `1.119Λ` held at the 336 row)

| variant | gap `1.0586Λ` | worst near-gap state | min over `m` |
|---|---|---|---|
| convention A — superseded | `+0.482` | `+0.472` [1.072] | `+0.472` |
| **direct small-`ν`, `best`, held below 8.35** | **`+0.993 ± 0.037`** (jk) | **`+0.910 ± 0.027`** [1.08] | `+0.872` [4Λ, no small-`ν` dependence] |
| `dir` / `cv` | `+0.969 / +0.993` | `+0.898 / +0.910` | `+0.872` |
| zero / linear below 8.35 | `+1.015 / +0.952` | `+0.927 / +0.879` | `+0.872` |
| `Γ_c := −1000 / −4000 / −8000` below 8.35 | `+0.96 / +0.81 / +0.60` | `+0.885 / +0.761 / +0.589` | failure at `≈ −2·10⁴ c₂` |

* Node values (`Γ_c`, `++`): row 336: `ν = 8.35`: `−823`; `10.0`: `+2`; `11.7`: `+99`; `13.4`: `+108`; `16.7`: `+86`; `20.9`: `+15`;
  `25.0`: `−35`; `33.4`: `−108`. Row 400: `8.36`: `−958`; `9.95`: `−100`; `11.7`: `+28`; `13.3`: `+61`; `16.7`: `+76`; `19.9`: `+45`;
  `29.9`: `−21`; `39.8`: `−64`.

### 4.4 `J = 4` (row 400, spot nodes) and the full `Δ_g = 300` sweep, even `J = 0–22`, odd `J = 5–15`
* `J = 4`: direct `Γ_c` (`++`) `−18 / +5 / +4 c₂` at `ν = 10.0 / 19.9 / 39.8` (failure would need `≲ −4·10⁴` below `ν = 25`, §4 method
  note); z4 with these nodes: gap state `1.0133Λ` **`+0.2137`**, identical (`±2·10⁻⁴`) to convention A. Closed [NE + M].
* `J = 5` (odd, binding; row 400, spot nodes, n3 with the odd `p8` integrand): direct `Γ` (`+−`) C2 `−13 / +2 / +2`, C3H `+119 / +15 /
  +16 c₂` at `x = 0.025 / 0.05 / 0.1` (`ν = 10 / 20 / 40`), against `K₂`/spline `−7 / +4 / +23` and `−40 / −40 / −39`; failure would need
  `≲ −2·10⁵ c₂` uniformly below `ν = 25`. With the direct nodes in the row (local de-shift; the odd sector has no non-local
  implementation): binding state **`+0.2038`** at `1.0166Λ`, unchanged. Closed [NE + M].
* **Sweep** (one z4 run: direct small-`ν` + non-local min-norm `X` for `J = 18, 20, 22`; convention A — `K₂` rows, local `ψ₂` weights
  — for every other `J`, which has no direct small-`ν` nodes; all exact-skeleton + hybrid states **including every `J`'s twist-edge
  mass**; `cache/hc_chain_z4_t3_all.txt`):

| `J` | 0 | 2 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 | 18* | 20* | 22* |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| min `λ/N_s` | +0.999 | +1.000 | +0.214 | **+0.2038** | +0.307 | +0.704 | +0.426 | +0.438 | +0.603 | +0.645 | +0.628 | +0.728 | +0.653 | +0.630 | +0.707 | +0.872 | +0.863 | +0.841 |
| at `m/Λ` | 1.0 | 1.7 | 1.013 (gap) | 1.017 (gap) | 1.25 | 1.1 | 1.052 (gap) | 1.03 | 2.0 | 1.036 (gap) | 1.055 (gap) | 1.1 | 1.046 (gap) | 1.068 | 1.1 | 4.0 | 1.065 (gap) | 1.071 (gap) |

  (`*` = direct small-`ν` nodes.) `n = 1` at 300 on all these `J`; binding state `J = 5` at its twist edge `1.0166Λ`, `+0.2038`
  (hybrid; uncorrected `+0.2043`, as the statement's `+0.2043`). The twist-edge masses lower several minima that the old grid did not
  sample (`J = 11: 0.749 → 0.645`, `12: 0.634 → 0.628`, `14: 0.690 → 0.653`, `15: 0.707 → 0.630`, `22: 0.123 → 0.046` in convention A).
* **Grade.** `J = 18, 20, 22` (and `J = 4`): **[NE + M]** with the validated small-`ν` treatment (margins `+0.84 … +0.87`, jackknife
  `±0.03–0.08`; failure needs a uniform `Γ_c ≲ −1.0·10⁴ / −1.45·10⁴ / −2·10⁴ c₂` below `ν = 8.3` at `J = 22 / 20 / 18`, against
  `O(10²–10³)` measured at `ν ≥ 8.3`). Even `J = 6–16` and odd `J = 5–15`: still convention A ([NE + M] in the superseded treatment, which
  at `J = 18–28` was always **pessimistic**, by `0.4–0.8 N_s`); the corrections there are `≤ 0.25 N_s` against margins `≥ 0.20`.

## 5. `K₂` beyond its convergence criterion: direct bulk nodes on `(336, 28)`, `(336, 32)`

Method: n1 adaptive direct nodes (shared seeds, `2¹²·2·8` points, element errors `0.06–0.16 %`) at `x = 0.2, 0.3, 0.4, 0.5, 0.6`,
compared node by node with the `K₂` row (`Γ`, curvature image subtracted; same `a_h`, which the injection leaves unchanged to `10⁻⁵`).

**`(336, 28)`, `J²/m = 2.35` [NE]:**

| `x` | C2 `++` direct / `K₂` | C2 `+−` | C3H `++` | C3H `+−` |
|---|---|---|---|---|
| 0.2 | `−46 / −61` | `−99 / −81` | `−58 / −224` | `+115 / +250` |
| 0.3 | `−88 / −157` | `−148 / −159` | `−112 / −289` | `+176 / +292` |
| 0.4 | `−106 / −257` | `−187 / −233` | `−144 / −383` | `+229 / +358` |
| 0.5 | `−61 / −402` | `−236 / −403` | `−111 / −535` | `+299 / +543` |
| 0.6 | `+258 / −412` | `−174 / −729` | `+219 / −519` | `+264 / +890` |

**`(336, 32)`, `J²/m = 3.08` [NE]:**

| `x` | C2 `++` direct / `K₂` | C2 `+−` | C3H `++` | C3H `+−` |
|---|---|---|---|---|
| 0.2 | `+153 / −86` | `+143 / −104` | `+196 / +40` | `−180 / +3` |
| 0.3 | `+157 / +105` | `+177 / +54` | `+201 / +271` | `−219 / −203` |
| 0.4 | `+67 / +397` | `+177 / +391` | `+103 / +552` | `−227 / −554` |
| 0.5 | `−265 / +306` | `−67 / +592` | `−282 / +316` | `+42 / −660` |
| 0.6 | `−279 / −917` | `−580 / −552` | `−400 / −1147` | `+684 / +750` |

* `K₂` is off by `15–170 c₂` at `x = 0.2–0.3` and by `150–670 c₂` at `x = 0.4–0.6`, on `(336, 28)` always on the same side on these
  nodes (the negative `P_J` lobes; `K₂` overshoots the lobe amplitude); on `(336, 32)` by `50–660 c₂` with both signs (phase error
  of the oscillation). This is the non-convergence the order ratios `|M₄/M₂| = 0.79 / 1.20` suggested: **`K₂` is not usable on these
  rows at `x ≥ 0.2` at the `10²`-`c₂` level [NE]**, confirming the review's regrading (S1).
* **A partial replacement is meaningless:** the `K₂` row oscillates with the `P_J` period (`Δx ≈ 0.1` at `x ≈ 0.5`), and replacing
  only the nodes `0.2, 0.3, …, 0.6` (one lobe sign) turns the `K₂` lobe error into a spurious non-oscillating shift (`dF_++/N` at the
  gap `+0.86 → +4.2`). Hence complete direct sub-rows `x = 0.125 … 0.8` (step 0.05, plus 0.125/0.15/0.175; 16 nodes) were computed.
* **`(336, 28)`, complete direct sub-row** (`cache/hc_t4f_28_nodes.txt`; `Γ` with the row's own `a_h`, which now includes the direct
  nodes at `x ≤ 0.25`): `K₂` overshoots the `P_J` lobes by a factor 2–5 on `x = 0.25–0.55` (C2 `++`: direct `+32, −44, +55, −60, +55, −22`
  vs `K₂` `+113, −157, +203, −257, +324, −402` at `x = 0.25 … 0.5`), and is out of phase at `x = 0.55–0.8` (`0.65`: `−607` vs `+61`;
  `0.75`: `−69` vs `−1660`; `0.8`: `−2503` vs `−326`).
  z4 chain (rows 336/500/1600; the 336 row now direct on `x ≤ 0.8`: z2 small-`ν` nodes on `ν ≤ 33`, hc nodes on `x = 0.125–0.8`; the
  500 and 1600 rows are `K₂` with `J²/m = 1.58 / 0.49`, inside the criterion):

  | thimble nodes `ν = 5.0, 6.67` | `K₂` bulk (as before) min / gap `1.091Λ` | **direct bulk** min / gap | jackknife (direct, hc nodes) |
  |---|---|---|---|
  | main run | `+0.813` [1.5] / `+1.538` | **`+0.835` [1.5] / `+1.728`** | `±0.002` at the minimum, `≤ ±0.026` over states |
  | independent repeat (`s2`) | `+0.660` [gap] / `+0.660` | **`+0.794` [1.5] / `+0.847`** | |

  **`J = 28` at 300 with direct bulk: `n = 1` [NE + M]**; the `K₂` bulk error was unfavourable (it lowered the gap state by `0.19` (main) and `0.19` (repeat)
  and the repeat-thimble minimum from `+0.85` to `+0.66`). The `ν < 8.3` window caveat (thimble values `±10⁴`, failure at a uniform
  `Γ_c ≈ −4000`) is unchanged: **`J = 28` stays [I] at the 1σ level below `ν ≈ 8`**, but is no longer [M] in the bulk.
* **`(336, 32)`, complete direct sub-row** (`cache/hc_t4f_32_nodes.txt`): `K₂` is wrong in amplitude and phase from `x = 0.2` on
  (C2 `++` at `x = 0.2 … 0.5`: direct `+107, −115, +117, −99, +41, +78, −258` vs `K₂` `−86, +11, +105, −252, +397, −461, +306`; `0.75`:
  `+996` vs `−2159`). z4 chain (rows 336/1040, z2 n3 small-`ν` nodes `ν ≥ 8.3`, no thimbles, `Γ` held below 8.3; twist-edge mass
  `1.102Λ` added):

  | | `K₂` bulk | **direct bulk** |
  |---|---|---|
  | min `λ/N_s` [`m/Λ`] | `+0.284` [1.102, the twist edge — never sampled before] | **`+0.783` [1.15]** (jk `±0.016`) |
  | gap state `1.102Λ` / 1.15Λ | `+0.284 / +0.657` | `+0.826 / +0.783` |

* **The `ν < 8.3` window with direct bulk** (`Z4_NUMIN=8`, uniform `Γ_c` below `ν = 8.3`): `J = 28`: `0 / −2000 / −4000 / −8000` →
  `+0.797 / +0.556 / +0.228 / −0.428` (failure at `≈ −5400 c₂`); `J = 32`: `+0.885 / +0.515 / +0.056 / −0.861` (failure at
  `≈ −4250 c₂`). Measured `Γ_c` at `ν = 8.3`: `−1175` (`J = 28`), `−895 … −1230` at `J = 20–28`; thimble values at `ν = 5–6.7`
  (`J = 28`) `±10⁴` (two seeds disagree by `≈ 10⁴`).
* Side result: with the independent thimble repeat at `ν = 5.0, 6.67` (`Z3_PREFER=x,s2`) and the twist-edge mass `1.091Λ`, `J = 28`'s
  `K₂`-bulk minimum is `+0.660` at the gap (`+0.691` at 1.1Λ, as the review found); with direct bulk `+0.794`.

### 5.1 Grade (`K₂` beyond its criterion, `J²/m ≳ 2`)
* **`K₂` on `(336, 28)`, `(336, 32)` is wrong at `x ≥ 0.2` by up to `×5` in amplitude and in phase [NE]** — the review's S1 regrading
  was right. Replacing it by a complete direct sub-row, `J = 28` and `J = 32` at `Δ_g = 300` **pass** (`+0.835 / +0.783`;
  repeat-thimble `+0.794` at `J = 28`) **[NE + M]** in the bulk; the `K₂` error had been *unfavourable* (`J = 32` at its twist edge:
  `+0.284 → +0.826`).
* Remaining for both: the window `ν < 8.3` — **[I]** at the level of the available data (failure needs a uniform
  `Γ_c ≲ −5400 / −4250 c₂` there; nodes at 8.3 are `≈ −10³`, the thimble nodes below are not resolved to `10⁴`).
* `J = 48` at 1000 (`J²/m = 2.2`, the other "`K₂`-limit" state) was not tested here.

## Hand-over (kept current)

* **Done:** every queue finished (`hc_queueA`, `A2`, `A3`, `B`, `OO`, `OO2`; logs `cache/hc_queue*.log`). No process is running.
  To recompute anything, the queue scripts are resumable (`setsid nohup ./hc_queueX.sh > cache/hc_queueX.log 2>&1 &` from
  `ads4/code/spinfunc/`; finished points are skipped).
* **Not done / open (value order):**
  1. OO block at `J ≥ 12` with the validated small-`ν` treatment (direct n3 nodes of the 3×3 `p8` integrand at even `J`;
     `hc_points.setup(odd=True)` would have to be forced for even `J` and the z4 NL term extended to the OO entry). `J = 18` at 300 is
     the binding block (`+0.22`) in the superseded treatment.
  2. The `ν < 8.3` window at 300 for `J = 28, 32` (and the `ν = 8.3` node itself, `≈ −10³` at 5 % RQMC error at every `J = 18–28`).
  3. `H_edge` on the uncomputed rows (1000: even `J ≤ 32`, 48, odd `J = 7–15, 19–29`; 300: even `J ≤ 22` except 16, 28, 32);
     same recipe: `HC_METHOD=n1 HC_ADAPT=1 HC_M=12–13 python3 hc_points.py J Δ 0.8,0.8125,…,0.9`, then `hc_edge.py` (1000) or
     `hc_chain.py w2/z4` with `HC_INJ`, `HC_DESPIKE_INJ=0`.
  4. Direct small-`ν` nodes for even `J = 6–16` and odd `J = 7–15` at 300 (still convention A); `J = 48` at 1000 against direct.
* **Code** (`ads4/code/spinfunc/`): `hc_points.py` (direct points; `HC_METHOD=n1|n3`, `HC_ADAPT`, `HC_M`, `HC_MREF`; per-state
  caches), `hc_nodes.py` (records, resid-matched shadows, injection), `hc_chain.py` (chains with injection: modes `ee1000`, `w2`,
  `z4`; env `HC_INJ`, `HC_IXMIN/IXMAX`, `HC_Z3`, `HC_DESPIKE_INJ`, `HC_EDGE_MODEL`, `HC_EK`, `HC_EREP`, `HC_JK`, `HC_ROWDUMP`; node
  report), `hc_edge.py` (edge variant table, ee1000 mode), `hc_t3.sh` (small-`ν` variants + jackknife), `hc_oo_rows.py`, `hc_oo.py`
  (OO block), `hc_sum.py`, `hc_jk.py` (summaries).
* **Caches:** `cache/hc_pts_<kind>_J<J>_D<Δ>.pkl` (+ `.lock`), `cache/hc_raw_J<J>_D<Δ>.pkl`, `cache/hc_p3_<kind>_J<J>.pkl`,
  `cache/hc_p2shadow_<kind>_J<J>.pkl`, `cache/hc_exact.pkl` (z3's exact scalar records, redirected), `cache/hc_chain_<mode>_<tag>.*`,
  `cache/hc_edge_J*.txt`, `cache/hc_oo_<tag>.*`, `cache/hc_t4f_{28,32}_nodes.txt`.
* **Reproduce (each ≤ a few minutes):**
  * edge, 1000: `python3 hc_edge.py 40 1040 1000` (also `36 1036`, `31 1031`, `17 1017`);
  * edge, 300 `J = 24`: `HC_TAG=t HC_INJ=24:336 HC_IXMAX=0.88 HC_DESPIKE_INJ=0 Z4_EST=best Z4_JS=24 Q3_JEVEN=24 Q3_JODD=5
    Q3_EXTRA_M=1.078,1.1,1.15,1.2,1.3,1.4,1.5,1.75,2 W2_ROWS=336,500,1600 SM_CURV=row python3 hc_chain.py z4 300; python3 hc_sum.py z4 300 t J=24`;
    `J = 16`: same in `w2` mode with `Q3_JEVEN=16 HC_INJ=16:400`; `J = 5`: `w2 300 1000` with `Q3_JODD=5 HC_INJ=5:200,400,800,1600 HC_IXMAX=0.92`;
  * OO: `HC_OOJ=4,6,8,10 Q3_EXTRA_M=1.02,1.05,1.1,1.15,1.2,1.3,1.4,1.5,1.75,2,2.5,3,4,5 python3 hc_oo.py 300 1000` (`HC_OOJ=18 … hc_oo.py 300`);
  * small `ν`: `./hc_t3.sh 22 22:336`, `./hc_t3.sh 20 20:400`, `./hc_t3.sh 18 18:336,400`; sweep: `HC_TAG=t3_all Z4_EST=best
    Z4_JS=18,20,22 HC_Z3='18:336,400;20:400;22:336' Q3_JEVEN=4,6,…,22 Q3_JODD=5,…,15 Q3_EXTRA_M=<gap masses>,1.1,…,5 python3 hc_chain.py z4 300`;
  * `K₂` check: `scratchpad t4full.sh` equivalent: `Z4_EST=best Z4_JS=28 Q3_JEVEN=28 Q3_JODD=5 W2_ROWS=336,500,1600 HC_INJ=28:336
    HC_IXMIN=0.12 HC_IXMAX=0.81 HC_DESPIKE_INJ=0 python3 hc_chain.py z4 300` (`J = 32`: `W2_ROWS=336,1040`).
* Analysis runs were made with one production worker SIGSTOPped, to stay at 2 busy processes.
