# SMALLNU_WINDOW: the `ν < 50` window of the large-J near-gap states at Δ_g = 1000, measured

*23 Nov 2026. Unreviewed, not committed by this agent. No existing file was edited; `DIRCV.md` and `code/dircv/` were not touched.*
* *New code:* `ads4/code/smallnu/`: `sn_guard.py` (process-local write guard: every write outside `smallnu/cache/` raises; it also
  stands in for `hs_guard` / `ng_guard`), `sn_points.py` (direct n1 / n3 RQMC points at small ν), `sn_core.py` (Γ at small ν in the chain's own
  convention and in the baseline normalisation, jackknife), `sn_window.py` (window weights in NEARGAP2's `E_low` convention, measured
  envelopes, signed window matrix), `sn_odd.py` / `sn_even.py` (budgets), `sn_oddrun.py` (HIGHSPIN's `hs_odd.py`, unchanged, run on a
  row with a direct bulk row computed here), `sn_summary.py`, `sn_run.sh` (compute lane), `sn_pause.sh` (short analyses with the lane
  SIGSTOPped).
* *Caches:* `ads4/code/smallnu/cache/` (its `.gitignore` ignores everything but itself); `n1chk/`, `n3chk/` hold the cross-check points.
* *Read-only imports:* `highspin/` (`hs_core`, `hs_odd`), `neargap2/` (`ng_bind`, `ng_val`, `ng_chain`), `spinfunc/` (`hc_points`,
  `hc_nodes`, `z2_points`, `p4_odd_analysis`, `l9_analysis`, `q3_analysis`, `ee1000_common`), their caches, and HIGHSPIN's result files
  (`highspin/cache/hs_odd_J*.json`, `hs_even_J*.json`).
* *Compute:* ≈ CPUH CPU-h, at most two computing processes of this note at any time (`sn_run.sh` lanes; analyses ran with one lane
  SIGSTOPped), `ulimit -v 3000000`, one thread everywhere.

Grades as in UNIFIED_THEOREM_v2 §0, POS_ERRORS, NEARGAP2 and HIGHSPIN: [P], [certified], [D], [D-sketch], [DV], [NE], [NE±] (explicit
error budget, every model input with a stated safety factor, conclusion surviving it), [NE + M], [M], [A], [I]. `c₂ = m² ×` relative
correction. `L = value / budget`. "Window" = `ν < 50` of the state (`ν = x m_state`; the chain uses its rows at fixed `x`).

---

## 0. Verdict

**The `10⁴ c₂` convention is replaced by measured values on `ν ≥ 10` (J = 16, 29, 32), `ν ≥ 12.5` (J = 12), `ν ≥ 16` (J = 20) and
`ν ≥ 20` (J = 21, 23, 27, 31).** `10⁴ c₂` is kept only *below* the lowest measured node, where it now carries 6 % (`ν < 10`) to 18 %
(`ν < 20`) of the window's weight, and is ×3.6 above the largest value known there (§4.4).

| state (Δ_g = 1000) | measured Γ on the window (`c₂`, baseline normalisation) | window cost: `10⁴` → envelope ×3 / measured | budget (with window) | `L` old → new | grade old → new |
|---|---|---|---|---|---|
| **J = 29, twist edge `1.02862Λ`** | `|Γ| ≤ 310` on ν = 25–47 (n1, σ ≤ 200); ν = 10–20: `≤ 1113` (n3/n1, σ 230–730, all within 2.1σ of 0); `−56 … +116` at 51–77 (HIGHSPIN) | `0.873` → `0.327` (×3) / **`0.093`** (signed: direct values on `[10, 50)`, `±0.041`; `ν < 10` at `10⁴`: `0.052`) | **`0.79`** (HIGHSPIN λ-level `0.697` + window `0.093`) | `0.61` (window at `10⁴`) → **`1.23`** (signed); `0.94` with the ×3 envelope | [NE + M] → **[NE±]** (thin; sign-free ×3 does not close at the twist edge, §3.1) |
| J = 29, 1.04 / 1.05 / 1.06 / 1.075 / 1.10Λ | same row | `0.80 … 0.61` → `0.30 … 0.25` (×3) / `0.085 … 0.070` (signed) | `0.43 / 0.35 / 0.30 / 0.23 / 0.17` | `0.83 … 1.38` → **`2.2 / 2.8 / 3.4 / 4.3 / 5.7`** (×3 envelope: `1.5 / 1.7 / 2.0 / 2.3 / 2.8`) | [NE + M] → **[NE±]** |
| **J = 16, twist edge `1.01589Λ`** | `|Γ| ≤ 104` on ν = 25–51; ν = 10–20: `≤ 648` (σ 110–350) | `0.305` → **`0.040`** (×3) / `0.021` (measured) | `0.21` sign-free | `1.92` → **`4.4`** | [NE±] → [NE±] (strengthened) |
| J = 16, 1.02 / 1.03 / 1.05Λ | same row | `0.27` → `0.035` | `0.21 / 0.29 / 0.59` sign-free | `1.69 / 1.43 / 1.20` → **`3.7 / 2.6 / 1.7`** (λ-level `5.2 / 7.1 / 10.0`) | [NE±] (strengthened; 1.05Λ limited by the 1600 K₂ row's truncation) |
| **J = 12, twist edge `1.01194Λ`** | `|Γ| ≤ 55` on ν = 25–51; ν = 12.5–20: `≤ 288` | `0.137` → **`0.018`** (×3) | `0.12` sign-free | `3.56` → **`7.4`** | [NE±] (strengthened) |
| J = 12, 1.02 / 1.05Λ | same row | `0.135 / 0.094` → `0.017 / 0.013` | `0.13 / 0.19` | `3.90 / 2.50` → **`7.5 / 3.6`** | [NE±] (strengthened) |
| **J = 21, twist edge `1.0208Λ`** (new: direct bulk row, 34 nodes) | `|Γ| ≤ 121` on ν = 25–47 (σ ≤ 106); ν = 20: `250 ± 393` (in the normalisation of the direct-bulk row; `≤ 275` in that of the K₂ row, §4.2) | `0.347` → **`0.102`** (×3) / `0.068` (signed; `ν < 20` at `10⁴`: `0.063`) | **`0.82`** (direct-bulk λ-level `0.720`, of which representation `0.55`, + window ×3) | `0.58` (K₂ bulk, `10⁴`) → **`1.20`** (×3; signed `1.26`) | [NE + M] → **[NE±]** (thin) |
| J = 21, 1.03 / 1.05 / 1.075 / 1.10Λ | same row; values `+0.885 / +1.005 / +0.968 / +0.959` | `0.29 … 0.22` → `0.09 … 0.07` | `0.51 / 0.35 / 0.15 / 0.16` | `0.65 … 0.83` → **`1.7 / 2.9 / 6.2 / 6.2`** (signed `1.8 / 3.1 / 7.5 / 7.2`) | [NE + M] → **[NE±]** |
| J = 23 (`1.0228Λ` … 1.10Λ) | `|Γ| ≤ 349` on ν = 25–47; ν = 20: `622 ± 512` | `0.40` → `0.15` / `0.08` | `1.84 + 0.15` (K₂ truncation `1.48`) | `0.40` → `0.45` (twist); calibrated truncation `1.5` | [NE + M] unchanged (window no longer binding; the K₂ bulk is) |
| J = 27 (`1.0267Λ` … 1.10Λ) | `|Γ| ≤ 361` on ν = 25–47; ν = 20: `793 ± 659` | `0.71` → `0.29` / `0.15` | `3.25 + 0.29` (truncation `2.95`) | `0.25` → `0.28`; calibrated `1.4` | [NE + M] unchanged |
| J = 31 (`1.0306Λ` … 1.075Λ) | `|Γ| ≤ 465` on ν = 25–47; ν = 20: `1068 ± 807` | `1.04` → `0.46` / `0.22` | `9.9 + 0.46` (truncation `9.5`) | `0.08` → `0.08` | [NE + M] unchanged (K₂ diverging, `J²/m = 0.93`) |
| J = 20 (rows 800 / 1600; `1.0198Λ`, 1.03, 1.05Λ) | `|Γ| ≤ 315` on state-ν 16–45 (row 1600) and 38–48 (row 800); row 800 at state-ν 26: `702 ± 404` | `0.56` → `0.23` (×3) / `0.17` (measured; `0.16` of it below ν_min) J20UPD | no budget exists (HIGHSPIN: value `+0.83 ± 0.14`) | — | [NE + M] unchanged |
| **J = 32 (row 1040; `1.0315Λ`, 1.03, 1.05Λ)** | `|Γ| ≈ 450–550` on ν = 25–40 (C2 and C3H, both helicities, σ 90–170); ν = 10–20: `≤ 1262` (σ 370–830); `|c_off| ≤ 5` | **`2.31`** → `0.73` (×3) / **`0.26`** (measured; the measured window matrix itself is `+0.027` to first order, favourable) | no budget exists (value `+0.997 … +1.006`) | — | [NE + M] unchanged; **the `10⁴` convention alone would have exceeded the state value ×2.3** |

**What this note establishes.**
1. **No state is negative, and no measured small-ν value is large.** On every row the direct Γ on `ν = 25–50` is `O(10²) c₂`:
   `≤ 55` (J = 12), `≤ 104` (J = 16), `≤ 121` (J = 21, direct-bulk normalisation), `≤ 310` (J = 29), `≤ 349–465` (J = 23–31, K₂-row
   normalisation), `≈ 500` (J = 32), i.e. ×20–180 below `10⁴` (the absolute level shifts by `≈ 10²` with the row normalisation, §4.2). Below ν = 20 the values
   are compatible with the same size; the direct method's errors there are `σ = 100–800 c₂` (§2.2).
2. **J = 29 moves to [NE±]** (twist edge `L = 1.23`, thin; `L ≥ 2.2` at 1.04–1.10Λ). The one remaining model input is `|Γ| ≤ 10⁴ c₂` on
   `ν < 10`, which the state tolerates up to `4.5·10⁴ c₂` (×4.5; ×16 above the largest value known there). With a *sign-free* ×3
   envelope the twist edge does **not** close (`L = 0.94`): the envelope there is set by Monte-Carlo noise of the n3 / near-fold n1
   nodes, not by Γ (bare `|Γ| + 3σ`: `L = 1.15`). This is reported as a negative finding (§3.1).
3. **J = 12, 16 stay [NE±] with ×2 larger leverage**; the window, which dominated their budgets, now costs `≤ 0.04 N_s`.
4. **For J = 21, 23, 27, 31 the window is no longer the binding input** (×3 envelope `0.12–0.46` vs `10⁴`: `0.35–1.04`), but their
   budgets with K₂ bulk do not close (λ-level truncation `1.2–9.5 N_s`). **J = 21 was given a direct bulk row here** (34 n1 nodes,
   ≈ 1.2 CPU-h): `λ/N_s = +0.998 ± 0.020` at the twist edge (K₂ `+0.980`), and with the measured window it moves to **[NE±]**
   (`L = 1.20`, ×3 envelope; `1.7–6.2` at 1.03–1.10Λ). J = 23, 27, 31 stay [NE + M] (they need the same direct bulk rows).
5. **The window grows steeply with J at fixed `ν`**: `W(ν < 50) = 1.4 / 3.1 / 3.5 / 5.6 / 8.7 / 23·10⁻⁵ N_s/c₂` at J = 12 / 16 / 21 / 20 /
   29 / 32 (twist edges, gauge factor included; even J with the off-diagonal). At J = 32 the old convention costs `2.3 N_s`; the
   measured window `0.26`, of which `0.14` is still the `ν < 10` convention. **For `J ≳ 32` the window is again a live input** even
   with measurements, because `σ(Γ)` of the direct method at `ν < 20` is `≥ 400 c₂` at Δ ≈ 10³.

---

## 1. The window and the convention it replaces

* NEARGAP2 §1.2: `E_low = G_low × W(ν < 50)`, `G_low = 10⁴ c₂`, each node weight × the min-norm gauge factor `max(1, |Re X_min/ψ₂|)`;
  HIGHSPIN used the same `10⁴` (and quoted `10 ×` the three lowest direct nodes as an alternative). Its support was the J = 5 engine
  (`≲ 10² c₂` away from the normalisation) and HC's Δ_g = 300 n3 nodes (`O(10²)` above ν = 10) — never a measurement at large J and
  Δ ≈ 10³.
* The window carries a small fraction of the weight but `N_s` falls steeply with J, so at `J ≳ 20` the convention dominated every
  budget: HIGHSPIN's J = 29 (`0.85` of a total `1.55`), J = 16 (`0.32` of `0.49`).
* Weight distribution inside the window (all states alike, gauge factor included): `ν < 8.3`: 4.7 %, `ν < 10`: 6.0 %, `ν < 12.5`:
  8.4 %, `ν < 20`: ≈ 18 %, `ν < 25`: 27 % (density ∝ ν; the gauge factor `≤ 6.84` at ν < 5).

---

## 2. Methods (`code/smallnu/`)

* **Write guard.** `sn_guard.py` (copy of `hs_guard` with root `smallnu/cache/`); it also registers itself as `ng_guard` and `hs_guard`,
  so importing `neargap2/` and `highspin/` modules installs no second guard. Every compute and analysis process imports it first.
* **Direct points** (`sn_points.py` → `spinfunc/hc_points.py`, unchanged; cache moved by its own `HC_CACHE` switch): state and
  `Δ_ref = 2·10⁴` at the same ν on common random numbers, 8 RQMC replicas.
  * **n1** (real re-centred contour, adaptive Gaussian proposal, `2¹¹·2·8` points; the HIGHSPIN/NEARGAP2 settings) at ν ≥ 20.
    **n1 fails below ν ≈ 20 at Δ ≈ 10³** (J = 29: replica errors 77–505 % on the records and up to 38 % of the points masked at
    ν = 12.5 and 15; `cache/n1chk/`, `cache/lane_b.log`): the fold of the crystal saddle (ν* ≈ 11.4, EE_DG300_SMALLNU §1) is felt up to ν ≈ 15–20.
  * **n3** (z2's hybrid contour, `ν_f = 13` below 13; HYPOTHESES_CLOSE's small-ν protocol, `2¹⁰·2·8` points) at ν = 10 … 17.5. It
    works, but its replica errors at Δ ≈ 10³ are 0.4–1.6 % on the records, i.e. `σ(Γ) ≈ 80–790 c₂` (J = 12 … 32) after the common-random-number
    cancellation (`m² ≈ 10⁶`).
* **Γ in the chain's own convention** (`sn_core.py`). The direct records are injected through `highspin/hs_core` (= the chain's
  rows, `neargap2/ng_val` loader wrappers, `hc_nodes` records, κ image subtracted for odd J) and read off the row exactly where the
  chain reads its K₂ nodes. **Normalisation.** Each row's constant `a` is fitted on its C2 nodes `x ≤ 0.25, ν ≥ 15`, and each record
  carries `exp(−c₀)` with `c₀` the log-ratio at the row's first node. A new small-ν node enters that fit and becomes the first node:
  the effective normalisation moves by `1–2.5·10⁻⁴`, which shifts Γ on the *whole* row by `≈ 10²–2.5·10² c₂ × D_h` (a near-pure
  renormalisation of the image, positivity-neutral to first order). All values here are transformed back to the **baseline
  normalisation** (HIGHSPIN's row): one ratio per helicity measured on the old nodes, `Γ′ = m²[(Γ/m² + S)/r − S]`,
  `S = sh_k/env_k + σ_k D_h` from the chain module's own functions; checked to reproduce every old node to `≤ 0.3 c₂`.
  Jackknife over the 8 replicas (all injected nodes together), each replica transformed with its own baseline.
* **Window weights** (`sn_window.py`): NEARGAP2's `E_low` structure exactly (`ng_bind.absdens`, `errmat`, `gauge_factor` — the
  min-norm/local de-shift factor `≤ 6.84` at ν < 5, `≤ 1.67` on 5–8, `≤ 1.04` above), 1024 hybrid nodes. (HIGHSPIN's odd-J `W_low`
  had no gauge factor: `8.52` vs `8.73·10⁻⁵ N_s/c₂` at J = 29.) The window is `ν_state < 50`; rows are used at fixed `x = ν/m` (the
  chain's convention), so a node at row-ν `ν_r` sits at state-ν `ν_r m_state/m_row`.
* **Envelopes.**
  * *Node bound* `g_i = |Γ_i| + 3σ_i`, max over C2/C3H (odd), per entry `++`, `+−` and `|c_off| = m|Γ_off|` (even).
  * *Measured envelope ×k* (`sn_window.envelope`): piecewise constant, on `[x_i, x_{i+1}]` `k·max(g_i, g_{i+1})`, beyond the last node its
    bound; **below the lowest reliable node `G_lo = 10⁴ c₂`** (the old convention, now carrying only 5–18 % of the window weight);
    `k = 3` (main) and `k = 10`.
  * *Signed (measured) window*: the direct values integrated over `[ν_min, 50)` (natural spline in x, held beyond the last node, rows
    1/m-interpolated as in the chain), error `3 × Σ|w_i|σ_i` (fully correlated) `+ 3 × |spline − linear|`; for even J the window's
    2×2 matrix enters through its spectral norm (the largest first-order effect on λ of either sign) plus the Frobenius norm of the
    error matrix; `G_lo = 10⁴` below `ν_min`.
* **Budgets.** All other terms are HIGHSPIN's (λ-level for odd J: MC, quadrature, missing GK nodes, representation, K₂ truncation;
  sign-free and λ-level for J = 12, 16 own rows). The odd-J central value of HIGHSPIN carried the bulk shape fit's continuation on
  ν < 50 (`X_low = +0.0075` at J = 29); it is removed here (`λ_c = main − X_low`) and the window is covered by the envelope (or
  replaced by the direct values in the signed variant), as `ng_bind`'s `lam_central` does for even J.

---

## 3. Results

### 3.1 J = 29 (row 1029, odd, `+−`) — `cache/sn_odd_J29.txt`

| ν | 10 | 12.5 | 15 | 17.5 | 20 | 25 | 30 | 35 | 40 | 45 | 51.4 (HS) | 64.2 (HS) | 77.0 (HS) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| method | n3 | n3 | n3 | n3 | n1 | n1 | n1 | n1 | n1 | n1 | n1 | n1 | n1 |
| C2 | −681 ± 634 | −173 ± 338 | +50 ± 361 | +66 ± 412 | +19 ± 233 | +230 ± 98 | +215 ± 91 | +169 ± 132 | +251 ± 174 | +152 ± 115 | −56 ± 59 | −91 ± 14 | −81 ± 32 |
| C3H | +1113 ± 528 | −87 ± 398 | −235 ± 381 | −393 ± 403 | −761 ± 727 | −283 ± 200 | −276 ± 137 | −249 ± 123 | −307 ± 165 | −107 ± 128 | +86 ± 80 | +115 ± 23 | +116 ± 34 |

(Γ − κ in `c₂`, jackknife σ over 8 replicas.) Γ is smooth and `O(250)` on ν = 25–45 with C2 ≈ −C3H (the pattern of the whole row),
changes sign near ν ≈ 48 (the state's leading partial wave has its first zero near ν ≈ 42); below ν = 20 every value is within
2.1σ of zero, with σ = 230–730.

| `m/Λ` | `λ_c` (window removed) | signed `λ` | HIGHSPIN budget (no window) | window: `10⁴` / ×1 / ×3 / ×10 / signed | `L`: `10⁴` / ×3 / ×10 / **signed** |
|---|---|---|---|---|---|
| **1.02862** | +0.964 | **+0.972** | 0.697 | 0.873 / 0.144 / 0.327 / 0.968 / **0.093** | 0.61 / 0.94 / 0.58 / **1.23** |
| 1.04 | +0.949 | +0.956 | 0.342 | 0.795 / 0.133 / 0.303 / 0.898 / 0.085 | 0.83 / 1.47 / 0.77 / **2.24** |
| 1.05 | +0.970 | +0.977 | 0.269 | 0.770 / 0.131 / 0.296 / 0.876 / 0.085 | 0.93 / 1.71 / 0.85 / **2.76** |
| 1.06 | +0.988 | +0.995 | 0.215 | 0.746 / 0.128 / 0.290 / 0.857 / 0.082 | 1.03 / 1.96 / 0.92 / **3.35** |
| 1.075 | +0.980 | +0.986 | 0.150 | 0.682 / 0.120 / 0.272 / 0.803 / 0.078 | 1.18 / 2.32 / 1.03 / **4.33** |
| 1.10 | +0.989 | +0.995 | 0.104 | 0.610 / 0.109 / 0.248 / 0.733 / 0.070 | 1.38 / 2.81 / 1.18 / **5.70** |

* `λ_c` = HIGHSPIN's main value minus its bulk-fit continuation on the window (`X_low = +0.0075 … +0.0061`). The signed window adds
  the direct values on `[10, 50)`: `+0.0076 ± 0.041` (twist edge) — the continuation HIGHSPIN used was right to `0.0001`.
* The signed budget at the twist edge: HIGHSPIN's `0.697` + MC/interpolation of the window `0.041` + `ν < 10` at `10⁴` `0.052` =
  `0.789`; **`λ ≥ +0.18`, `L = 1.23`**. The `ν < 10` input is tolerated up to a uniform `|Γ| = 4.5·10⁴ c₂` (`1.2–2.2·10⁵` at 1.04–1.10Λ).
* **Sign-free envelopes do not close at the twist edge**: ×3 gives `L = 0.94`, ×10 `0.58`. Their node bounds (`|Γ| + 3σ`) on ν = 10–20
  are `1.3–2.9·10³ c₂` — set by σ, not by Γ. To close a sign-free ×3 envelope there would need σ ≲ 150 c₂ at ν = 10–20, i.e. ≈ 16×
  the n3 points (≈ 2.5 CPU-h per node at J = 29); not attempted.

### 3.2 J = 12, 16 (own rows 1012, 1016, even) — `cache/sn_even_J12.txt`, `sn_even_J16.txt`

| row | ν | 10 | 12.5 | 15 | 16 / 17.5 | 20 | 25.4 (HS) | 30 | 35 / 37.5 | 42 / 45 | 50.5 (HS) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1016 | C2 `++` | −362 ± 257 | −175 ± 113 | −96 ± 118 | −87 ± 139 | −27 ± 69 | +20 ± 20 | +37 ± 12 | +48 ± 20 | +53 ± 31 | +47 ± 39 |
| 1016 | C3H `+−` | +648 ± 202 | +102 ± 121 | +55 ± 122 | −8 ± 124 | −400 ± 347 | −39 ± 70 | −54 ± 55 | −77 ± 44 | −88 ± 39 | −82 ± 37 |
| 1012 | C2 `++` | — | −130 ± 82 | — | −70 ± 88 | +1 ± 55 | −4 ± 20 | — | +7 ± 7 | +10 ± 4 | +7 ± 10 |
| 1012 | C3H `+−` | — | +107 ± 108 | — | +41 ± 77 | −288 ± 259 | +3 ± 59 | — | −11 ± 38 | −20 ± 29 | −21 ± 23 |

(`+−` C2 and `++` C3H follow the same pattern; `|c_off| ≤ 1`.) Budgets: HIGHSPIN §4's terms unchanged, `E_low` replaced.

| state | HIGHSPIN `E_low` (`10⁴`) | ×1 / ×3 / ×10 / measured | other terms (sign-free / λ-level) | `L` sign-free old → ×3 (×10, measured) | λ-level old → ×3 |
|---|---|---|---|---|---|
| J = 12, 1.01194 | 0.142 | 0.013 / 0.018 / 0.034 / 0.012 | 0.098 / 0.098 | 3.56 → **7.40** (6.48, 7.76) | 3.56 → 7.40 |
| J = 12, 1.02 | 0.140 | 0.013 / 0.017 / 0.034 / 0.012 | 0.116 / 0.055 | 3.90 → **7.50** (6.68, 7.82) | 5.12 → 13.8 |
| J = 12, 1.05 | 0.098 | 0.010 / 0.013 / 0.024 / 0.009 | 0.178 / 0.014 | 2.50 → **3.61** (3.40, 3.68) | 6.15 → 25.7 |
| J = 16, 1.01589 | 0.317 | 0.025 / 0.040 / 0.097 / 0.021 | 0.172 / 0.172 | 1.92 → **4.43** (3.49, 4.86) | 1.92 → 4.43 |
| J = 16, 1.02 | 0.283 | 0.022 / 0.035 / 0.085 / 0.019 | 0.171 / 0.111 | 1.69 → **3.72** (3.00, 4.04) | 1.95 → 5.23 |
| J = 16, 1.03 | 0.267 | 0.022 / 0.034 / 0.083 / 0.019 | 0.259 / 0.072 | 1.43 → **2.56** (2.20, 2.70) | 2.21 → 7.07 |
| J = 16, 1.05 | 0.281 | 0.023 / 0.037 / 0.089 / 0.020 | 0.554 / 0.063 | 1.20 → **1.69** (1.56, 1.74) | 2.91 → 10.0 |

### 3.3 J = 21, 23, 27, 31 (rows 1000 + J, odd) — `cache/sn_odd_J{21,23,27,31}.txt`

n1 nodes at ν = 20, 25, 32.5, 40, 47 (n1 fails below 20; no n3 run on these rows). On ν = 25–47 all four rows show the J = 29 pattern
(C2 `+107 … +331`, C3H `−152 … −361`, σ 49–178); at ν = 20 C3H is `−474 / −622 / −793 / −1068` with σ `442 / 512 / 659 / 807`.

| J | twist-edge `λ_c` | HIGHSPIN budget (no window; K₂ truncation) | window `10⁴` → ×3 / signed (`ν < 20` at `10⁴`) | `L` strict: `10⁴` → ×3 / signed | `L` with the calibrated truncation (5 % of the step, HIGHSPIN §4): ×3 / signed |
|---|---|---|---|---|---|
| 21 | +0.992 | 1.357 (1.204) | 0.347 → 0.118 / 0.069 (0.063) | 0.58 → 0.67 / 0.70 | 3.0 / 3.5 (1.03–1.10Λ: 4.7–6.9 / 6.1–9.1) |
| 23 | +0.898 | 1.843 (1.480) | 0.404 → 0.149 / 0.083 (0.074) | 0.40 → 0.45 / 0.47 | 1.5 / 1.7 (2.5–4.7 / 3.0–6.2) |
| 27 | +0.996 | 3.251 (2.951) | 0.714 → 0.286 / 0.151 (0.130) | 0.25 → 0.28 / 0.30 | 1.4 / 1.7 (1.9–2.5 / 2.6–3.3) |
| 31 | +0.841 | 9.889 (9.472) | 1.040 → 0.461 / 0.222 (0.189) | 0.08 → 0.08 / 0.08 | — (truncation step itself 3.2 `N_s`) |

(J = 31's value is HIGHSPIN's `hs_odd` run with spline edge, `+0.786` + window removal; the published refzero-corrected shape chain
gives `+1.02`. The K₂ series on (1031, 31) is not converging, `J²/m = 0.93`.)

**J = 21 with a direct bulk row** (`cache/hs_odd_J21.txt` — HIGHSPIN's `hs_odd.py`, unchanged, run by `sn_oddrun.py`; 34 n1 nodes
`x = 0.05 … 0.775`, ≈ 1.2 CPU-h; GK15 edge nodes of HIGHSPIN). Direct vs K₂ in the bulk moves the twist edge `+0.980 → +0.998 ± 0.020`
(1.03Λ `+0.868 → +0.885`, all masses `+0.017 … +0.021`, favourable, as at J = 12, 16, 29). λ-level budget without the window (HIGHSPIN
convention): MC `0.059`, quadrature `0.022`, missing GK nodes `0.088`, representation `0.550` (the degree-7 shape fit gives `+0.815`;
degrees 4, 5 and both half-samples `+1.000 … +1.009`), truncation 0 → `0.720`, `L = 1.39` (1.03 / 1.05 / 1.075 / 1.10Λ: `2.1 / 4.0 /
13 / 12`). The window (in this row's own normalisation, `cache/sn_odd_J21.txt`): C2 `+37 … +76`, C3H `−74 … −121` on ν = 25–47,
`−34 ± 80 / −250 ± 393` at ν = 20; cost `10⁴` `0.347` → ×3 `0.102`, ×10 `0.193`, signed `0.068`. **Twist edge `L = 1.20` (×3) / `1.26`
(signed) / `1.08` (×10)**; 1.03Λ `1.71 / 1.81`; 1.05Λ `2.87 / 3.14`; 1.075Λ `6.24 / 7.47`; 1.10Λ `6.15 / 7.19`. **[NE±]**, thin at the twist
edge, where the representation term (×3 on the degree-7 fit) dominates.

### 3.4 J = 20 (rows 800, 1600) and J = 32 (row 1040), even — `cache/sn_even_J20.txt`, `sn_even_J32.txt`

* **J = 32** (state on the 1040 row): n1 at ν = 20–48, n3 at 10–17.5. `Γ ≈ +450 … +550` for C2 (`++`, `+−`) and C3H `++` (`−` for C3H
  `+−`) on ν = 25–40 (σ 90–170), `+68 … +148` at ν = 48 (the K₂ node at 51.9: `|Γ| = 56 … 207`); ν = 10–17.5 all within 1.6σ of 0
  except C3H `+−` at ν = 10 (`+1198 ± 628`). `|c_off| ≤ 5`. Window at the twist edge: `10⁴` **`2.305`** (this alone exceeds the state's
  `λ/N_s ≈ 1.0` ×2.3) → ×1 `0.32`, ×3 `0.73`, ×10 `2.18`; **measured `0.257`** = `|dF_win|₂ 0.038` (first-order effect on the hybrid
  `+0.027`, favourable) + error `0.081` + `ν < 10` at `10⁴` `0.138`. 1.03 / 1.05Λ: `0.259 / 0.235`.
* **J = 20** (state at `t = 0.43` between rows 800 and 1600): J20TEXT

---

## 4. Checks

### 4.1 Two contours at the same ν
* J = 29, ν = 25: n1 `C2 +230 ± 98, C3H −283 ± 200` vs n3 `−4 ± 441, −347 ± 409` (`cache/n3chk/`): differences `0.5σ / 0.1σ`.
* J = 16, ν = 30: n3 `C2 ++ −93 ± 183, +− −79 ± 172, C3H ++ +83 ± 173, +− −17 ± 149` vs n1 `+37 ± 12, +40 ± 11, +90 ± 31, −54 ± 55`
  (0.7 / 0.7 / 0.04 / 0.2σ). The two contours agree; n3 is ×3–15 noisier at equal cost.
* n1 vs n3 at ν < 20 is not possible (n1 fails); the n3 values at 10–17.5 join the n1 values at 20–25 smoothly within their errors
  on every row (§3).

### 4.2 Normalisation and the transform
* Injecting new nodes into a row moves its fitted normalisation by `r − 1 = −3·10⁻⁶ … −2.6·10⁻⁴` (J = 12 … 32), i.e. Γ on the whole row
  by up to `−260 c₂ × D_h` in the injected grid. In the baseline normalisation (the transform of §2) every old node is reproduced to
  `≤ 0.5 c₂` (odd rows, x < 0.94) and `≤ 10⁻⁹ c₂` (even rows). The quoted envelopes are therefore those of HIGHSPIN's rows, i.e. of the
  state values whose budgets are recomputed. A uniform shift of this kind is a multiplicative change of the image and does not change
  the sign of λ.
* **The absolute small-ν level of Γ depends on the row's normalisation at the `10²`-`c₂` level.** J = 21: with the K₂ bulk the
  small-ν C2 / C3H values are `+107 … +268 / −152 … −474`; after the direct bulk row is injected (which refits `a` on direct instead of
  K₂ nodes, `δa/a ≈ 1.9·10⁻⁴`) the *same* points read `+37 … +76 / −74 … −250` — a uniform `∓190 c₂ × D_h` shift. Envelopes must therefore
  be quoted in the normalisation of the state value they are combined with, which is what is done throughout (J = 21: direct-bulk
  normalisation for the direct-bulk value; J = 23–31: K₂-row normalisation for HIGHSPIN's K₂ values). The J = 21 K₂-bulk numbers of
  §3.3 come from the first `sn_odd.py 21` run, before the bulk row existed (its cache file now holds the direct-bulk version).

### 4.3 Reference zeros and near-degenerate references
* The ratio's reference is `Δ_ref = 2·10⁴` at the same ν: at ν < 50 its `x_ref < 0.0025`, `(2J+1)x_ref ≤ 0.17`, far from its partial
  wave's first zero (`2.405`); HIGHSPIN's "reference-zero" failure cannot occur in the window. The *state's* partial wave has its first
  zero inside the window for J ≳ 27 (ν₀ ≈ 2.405 m/(2J+1): 45 / 42 / 39 / 38 at J = 27 / 29 / 31 / 32), but Γ is not divided by it
  (Γ = m²(G/(a env) − D_h)); the direct values pass smoothly through ν₀.

### 4.4 What is known below the lowest node (the remaining `10⁴` input)
* Direct n3 at ν = 10 (Δ ≈ 10³): `|Γ| ≤ 648` (J = 16), `≤ 1113` (J = 29), `≤ 1198` (J = 32), each with σ 200–790.
* HC's Δ_g = 300 rows (Δ = 336 / 400, n3, `m² ≈ 10⁵`, so ×10 smaller errors in `c₂`): `ν = 8.3`: `−823 … −1230`; `ν = 10`: `−182 … +2`;
  `ν = 11.7–40`: `|Γ| ≤ 120` (J = 18, 20, 22; HYPOTHESES_CLOSE §4; z4 convention).
* The bj engine's order-1 remainder at J = 28 (POS_ERRORS §2.3): `≤ 2.8·10³` (ν < 2), `620` (2–6), `150` (6–12) `c₂`.
* So `10⁴` on `ν < 10` is ×3.6 above the largest value known there; the J = 29 twist edge tolerates `4.5·10⁴` (×16).

---

## 5. Routes 2 and 3

* **Route 2 (analytic-radial / bj engine) was not run.** At J = 5 (NEARGAP2 §1.3) its last-order term was already `≈ 600 c₂` (mostly a
  normalisation); POS_ERRORS §2.3 shows the engine's spin-J remainder at `J²/Δ ≈ 0.04–0.5` is `10²–10³ c₂` below ν = 12 and that it says
  nothing at ν ≥ 20 for J > 0 (depth truncation `10³–10⁵ c₂`). At J = 29 (`J²/Δ = 0.82`) its truncation would exceed the direct
  envelope on the whole measured range; its only use would be ν < 10, where §4.4's engine numbers are already quoted. The compute went
  to direct nodes.
* **Route 3 (Lemma C / C′) does not bound Γ.** Lemma C bounds `|P_{−½+iν}(cosh β) − J₀(νB)| ≤ (1 − sech(β/2)) + ν(2 sinh(β/2) − β)`,
  i.e. the *leading-Regge curvature*, which is in the state's base (odd J: κ subtracted from Γ). In `c₂` units it is
  `≈ 435 + 8.3ν` at J = 29 (`136 + 1.5ν` at J = 16, `528 + 11ν` at J = 32) — comparable to the measured Γ, but it is not Γ: Γ is the
  beyond-leading-Regge remainder δ, for which LARGESPIN_ANALYTIC (§3.2, §5.3, table row (a)) shows that no sign-free per-ν argument can work. No explicit
  analytic constant for Γ on ν < 50 follows.

---

## 6. Caveats

1. The direct n1/n3 integrals are the project's reference standard [NE], not certified. n1 and n3 agree where both run (§4.1).
2. Below ν ≈ 20 the direct errors at Δ ≈ 10³ are 100–800 `c₂`: the window there is measured to "small compared with `10⁴`", not to
   the size of Γ itself. Sign-free envelopes are noise-dominated on ν < 20 (J = 29: bare `|Γ| + 3σ` already gives `L = 1.15`).
3. The signed window (J = 29 main statement) uses a natural spline in x through the direct nodes (×3 on spline − linear) and a fully
   correlated MC error; between nodes Γ is assumed smooth on the node spacing (Δν = 2.5–7.5 against an oscillation period ≈ 110 at J = 29).
4. `|Γ| ≤ 10⁴ c₂` on `ν < 10` (J = 29, 16, 32), `< 12.5` (J = 12), `< 16–26` (J = 20), `< 20` (J = 21–31) remains a model input
   (stated safety factor ×3.6 vs the largest known value, §4.4). Below ν ≈ 8 there is no direct value at Δ ≈ 10³.
5. Only the masses of HIGHSPIN's tables were evaluated; between them, `H_cont` as before. All non-window budget terms are HIGHSPIN's,
   unchanged (for J = 21 recomputed with the same code, §3.3).
6. The even-J "measured" window uses the spectral norm of the window matrix, not its signed effect in the full chain matrix.

---

## 7. Suggested corrections to existing documents (not applied)

1. **HIGHSPIN §0 / §2.2 / §6 and the update boxes (NEARGAP2, POS_ERRORS, RESULT_D_STATEMENT, EE_DG1000_LARGEJ):** J = 29 near gap
   "[NE + M], `ν < 50` window the one model input" → "**[NE±]** (SMALLNU_WINDOW): direct Γ on ν = 10–50 (`|Γ| ≤ 310 c₂` on 25–47,
   within 2.1σ of 0 below 20); twist edge `L = 1.23` (signed window), `L ≥ 2.2` at 1.04–1.10Λ; remaining input `|Γ| ≤ 10⁴` on ν < 10
   (tolerance `4.5·10⁴`)."
2. **HIGHSPIN §4 / §6, J = 12, 16:** leverages `3.6 / 3.9 / 2.5` → `7.4 / 7.5 / 3.6` (J = 12) and `1.9 / 1.7 / 1.4 / 1.2` →
   `4.4 / 3.7 / 2.6 / 1.7` (J = 16).
3. **NEARGAP2 §1.2 `E_low` / HIGHSPIN §1 "Budgets":** `G_low = 10⁴ c₂` on `ν < 50` should be replaced by the measured envelope on
   `ν ≥ 10–20` and kept only below; HIGHSPIN's odd-J `W_low` omits the gauge factor (`8.52` vs `8.73·10⁻⁵` at J = 29, +2.5 %).
4. **HIGHSPIN §0 "Negative findings" 4 / POS_ERRORS caution "the ν < 50 window is now the binding model input at large J":** true for the
   convention, false for the measured window at J = 12–31: the binding input for J = 23–31 is the K₂ bulk truncation, and for J = 21 (direct bulk row) and J = 29 the representation /
   missing-node terms of the bulk and edge, not the window.
5. **EE_DG1000_LARGEJ §3 item 4 ("direct RQMC cannot reach `m² × 10⁻⁴` precision at `m ≈ 10³`"):** with common random numbers it reaches
   `σ = 4–200 c₂` (`≲ 2·10⁻⁴` of `m²`) on ν = 25–50 (n1) and `100–800 c₂` on ν = 10–20 (n3) at J = 12–32.
6. **For J = 32 (and any J ≳ 30) the window must be quoted with its measured value**: under the convention it costs `2.3 N_s`.

## 8. Reproduce
```bash
cd ads4/code/smallnu      # compute lanes: bash sn_run.sh <jobfile> <log>  (ulimit -v 3000000, one thread); jobs: cache/lane_{A,B}.jobs
python3 sn_points.py n1 29 1029 20,25,30,35,40,45          # n1 small-nu nodes (explicit nu); likewise 21/1021, 23/1023, 27/1027, 31/1031
python3 sn_points.py n3 29 1029 10,12.5,15,17.5            # n3 nodes below nu = 20; likewise 16/1016, 12/1012 (12.5,16), 32/1040, 20/800, 20/1600
python3 sn_points.py n3 29 1029 25 n3chk                   # cross-checks in separate caches (n1chk/, n3chk/)
python3 sn_odd.py 29;  python3 sn_odd.py 21; ... 31        # odd-J windows and budgets  -> cache/sn_odd_J<J>.txt/.json
python3 sn_even.py 16 1016; python3 sn_even.py 12 1012; python3 sn_even.py 32 1040; python3 sn_even.py 20 800,1600
python3 sn_points.py n1 21 1021 <34 nu of x = 0.05 ... 0.775>; python3 sn_oddrun.py 21 1.03,1.05,1.075,1.1   # J = 21 direct bulk row
python3 sn_summary.py
```
