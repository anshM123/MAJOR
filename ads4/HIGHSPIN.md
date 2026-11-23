# HIGHSPIN: the high-spin near-gap states of Result D at Δ_g = 1000 — J = 29 resolved, J = 12, 16 error-bounded, the edge sweep, H_LJ(1000)

*23 Nov 2026. Unreviewed, not committed by this agent (another agent's commits db03df3 / 0e58356 swept in intermediate snapshots of this
file and of `code/highspin/`; the working tree holds the final versions). No existing file was edited; `R4ARC_HEFT.md` and
`code/r4heft/` were not touched. One write escaped the guard: Python's bytecode cache `code/highspin/__pycache__/` (from importing
`hs_core` in scratch scripts; git-ignored); it was deleted.*
* *New code:* `ads4/code/highspin/` — `hs_guard.py` (process-local write guard: every write outside `highspin/cache/` raises),
  `hs_core.py` (rows with direct nodes injected, reference-zero diagnostic, κ image, shape model, GK15 quadrature), `hs_points.py` (direct
  n1 RQMC points), `hs_odd.py` (odd-J state with direct bulk / GK edge and budget), `hs_even.py` (even-J own-row state with direct nodes,
  through `neargap2/ng_bind`), `hs_edgeval.py` (edge injection on non-own-mass even rows, `ng_val` protocol), `hs_nodes.py` (node tables),
  `hs_refzero.py` (reference-zero survey), `hs_ee1000.py` (published shape chain with/without the reference-zero node), `hs_lane.sh`,
  `hs_run.sh`.
* *Caches:* `ads4/code/highspin/cache/` (its `.gitignore` ignores everything but itself). Existing direct points of NEARGAP2 / HC were
  *copied* there (`hs_points.py J D seed`), not moved.
* *Read-only imports:* `neargap2/` (`ng_val`, `ng_chain`, `ng_bind`), `spinfunc/` (`hc_points`, `hc_nodes`, `ee1000_common`, `ee1000_chain`,
  `sm_curvfix`, `w2_chain`, `q3_analysis`, `p4_odd_analysis`), and their caches.
* *Compute:* ≈ 7.2 CPU-h, at most two computing processes of this note at any time (short analyses ran while one worker was SIGSTOPped,
  `hs_run.sh`; one accidental 40-s overlap of a third worker was killed), `ulimit -v 3000000`, one thread everywhere.

Grades as in UNIFIED_THEOREM_v2 §0 and POS_ERRORS: [P], [certified], [D], [D-sketch], [DV], [NE], [NE±] (explicit error budget, every
model input with a stated safety factor, conclusion surviving it), [NE + M], [M], [I]. `c₂ = m² ×` relative correction. Leverage
`L = value / budget`. "Sign-free" and "λ-level" budgets as NEARGAP2 §1.2.

---

> **Update (SMALLNU_WINDOW.md, unreviewed): the ν < 50 window is now measured. The 10⁴ c₂ convention is kept only below the lowest measured node (6–18 % of the window weight). Nothing came out negative. J = 29 and J = 21 → [NE±] (thin); J = 12, 16 leverage roughly doubles; J = 23, 27, 31 stay [NE+M] (now limited by the K₂ bulk truncation, not the window).**
> * **Method:** direct n1 points at ν ≥ 20 and n3 points at ν = 10–17.5 on each state's own row, with jackknife errors. n1 fails below ν ≈ 20 (crystal-saddle fold); n1 and n3 agree at ν = 25 and 30. Reference zeros cannot occur in the window.
> * **J = 29 twist edge:** measured Γ ≤ 310 c₂ on ν = 25–47; ν = 10–20 is consistent with 0 (σ up to 730). Signed window: L 0.61 → **1.23**, λ ≥ +0.18. The remaining model input is |Γ| ≤ 10⁴ on ν < 10; the state tolerates 4.5·10⁴ there, 16× the largest known value. Masses 1.04–1.10Λ: L 2.2–5.7. **A sign-free ×3 envelope does not close (L = 0.94)**, because of Monte Carlo noise at ν = 10–20.
> * **J = 21:** a new 34-node direct bulk row puts the twist edge at +0.998 ± 0.020, L 0.58 → 1.20 (1.7–6.2 at 1.03–1.10Λ). **J = 16:** L → 4.4 / 3.7 / 2.6 / 1.7. **J = 12:** L → 7.4 / 7.5 / 3.6.
> * **J = 23 / 27 / 31:** window cost now 0.15 / 0.29 / 0.46, but the K₂ λ-level truncation blocks them (with calibrated truncation L ≈ 1.5 / 1.4). **J = 32:** the window costs 0.26 measured, vs 2.31 under the convention. J = 20 and J = 32 have no budget.
> * **Cautions:** the window weight grows steeply with J (1.4 → 23·10⁻⁵ N_s/c₂ from J = 12 to 32), so for J ≳ 32 it is a live input again. Lemma C/C′ cannot bound Γ (it is beyond leading Regge). HIGHSPIN's odd-J window weight omits a gauge factor (8.52 vs 8.73·10⁻⁵). EE_DG1000_LARGEJ's "direct RQMC cannot reach the precision" is false on ν = 25–50.

## 0. Verdict

| state (Δ_g = 1000) | value `λ/N_s` (old → new) | budget | `L` | grade (old → new) |
|---|---|---|---|---|
| **J = 29, twist edge `1.02862Λ`** | w2 chain `−3.62` / with 5 direct edge nodes `+1.39` / published `+0.948` → **`+0.972 ± 0.018`** (direct bulk row + GK15 direct edge) | `0.70` λ-level (MC 0.055, quadrature 0.11, missing/reference-zero GK nodes 0.48, representation 0.05); `ν < 50` window separately | **1.39**; with `ν < 50` sign-free at `G = 10⁴ c₂`: 0.63; at `G = 10 ×` measured (`1.2·10³`): 1.22 | [M] → **[NE + M]**: positive on every representation tried; error-bounded except the `ν < 50` window, which needs a uniform `|Γ| ≥ 3.2·10³ c₂` there to fail (direct `|Γ| ≤ 116 c₂` at ν = 51–77, ×28) |
| J = 29, 1.04 / 1.05 / 1.06 / 1.075 / 1.10Λ | `+0.956 / +0.977 / +0.995 / +0.986 / +0.995` (± 0.017) | `0.34 / 0.27 / 0.21 / 0.15 / 0.10` | `2.8 / 3.6 / 4.6 / 6.6 / 9.5` (with `ν < 50` at `10⁴`: `0.86 / 0.96 / 1.05 / 1.21 / 1.42`) | [M] → [NE + M] (`ν < 50` only); `ν < 50` tolerance `7.9·10³ … 1.5·10⁴ c₂` |
| **J = 12, twist edge `1.01194Λ`** | own-row K₂ `+0.8497` → **`+0.8528 ± 0.0004`** (full direct row, 38 nodes) | `0.24` sign-free (`ν < 50` at `10⁴`: 0.14; edge `x > 0.925`: 0.094) | **3.6** (6.2 with direct nodes standing for `ν ∈ [25, 50]`) | [NE + M] → **[NE±]** |
| J = 12, 1.02 / 1.05Λ | `+1.0006 / +0.6872` | `0.26 / 0.28` sign-free | `3.9 / 2.5` (λ-level `5.1 / 6.2`) | [NE + M] → **[NE±]** |
| **J = 16, twist edge `1.01589Λ`** | own-row K₂ `+0.9342` → **`+0.9402 ± 0.0012`** (full direct row) | `0.49` sign-free (`ν < 50`: 0.32; edge: 0.16) | **1.9** (3.6 with `ν ≥ 25` direct) | [NE + M] → **[NE±]** |
| J = 16, 1.02 / 1.03 / 1.05Λ | `+0.7681 / +0.7509 / +1.0009` | `0.45 / 0.53 / 0.84` sign-free (`0.39 / 0.34 / 0.34` λ-level) | `1.7 / 1.4 / 1.2` (λ-level `2.0 / 2.2 / 2.9`; with `ν ≥ 25` direct `3.1 / 2.3 / 1.6`) | [NE + M] → **[NE±]** (1.03, 1.05Λ draw on the 1600 K₂ row; their margin is set by its truncation) |
| J = 31, 32 (published shape chain) | J = 31 twist edge `+1.208 → +1.023`; J = 32 `1.0315Λ` `+0.904 → +0.997` (reference-zero node removed, K = 4) | — | — | [NE + M] unchanged; published values carried a K₂ node at the reference zero (§2). Both shifts favourable |
| edge sweep (§3): odd J = 21, 23, 27 | twist edge `+0.984 / +0.880 / +0.956` → `+0.980 / +0.883 / +0.959` (GK15 direct edge); minima `+0.868` (J = 21, 1.03Λ), `+0.883`, `+0.928` | edge terms `≤ 0.33`; bulk K₂ step 0.4–1.0 not closed at ×3 | — (≈ 2–6 with the §4 calibration) | edge [M] → [NE + M]; states [NE + M] |
| edge sweep (§3): even J = 20, 22, 28, 32 | edge shifts `−0.13 ± 0.14` (J = 20, MC-limited), `−0.017`, `−0.016`, `+0.093` | — | — | edge [M] → [NE + M]; states [NE + M] |
| J = 10, 14 | not run (no own-mass K₂ rows; building them and the direct rows did not fit the budget) | | | [NE + M] unchanged |
| J = 26, 30 (even), J = 48 edge | not run: no K₂ rows exist at J = 26, 30 (interpolated in J, inside `H_LJ(1000)`); J = 48's edge needs the fw structure basis (`fw_frame.patch_conserved_clean`) | | | [M] unchanged |
| **`H_LJ(1000)`** | — | — | — | **[D-sketch + M] unchanged** (§5): no compute route to [NE±] within budget; what is model is listed precisely |

**Negative findings (reported prominently).**
1. **No physical state was found negative.** Every representation of every state computed here is positive except the two chain
   artefacts at J = 29 explained in §2 (`−3.62`, and `−1.17` in the ee1000 shape chain *without* its Hampel step).
2. **The J = 29 discrepancy is a method failure that affects other rows.** The K₂ outlier at `(1029, 29)`, x = 0.8 is the **heavy
   reference** (`Δ_ref = 2·10⁴`) passing through the first zero of its leading partial wave at `ν₀ ≈ 2.405·Δ_ref/(2J+1)` (§2). The same
   happens on 9 other chain rows for `J ≥ 18` (J = 31 and J = 32 in the *published* chains, uncorrected; J = 40, 48 were caught as
   unexplained outliers). Removing the node raises the published J = 31, 32 near-gap values (favourable).
3. **The K₂ λ-level Laplace step is not a truncation estimate at large J.** At the J = 19–31 near-gap states the K = 1 → 2 step of
   `λ/N_s` is `0.3–3.2` (×3 budgets of 1–10 `N_s`); direct rows show the actual K₂ error is only **4–5 % of that step** (J = 12: 0.003
   vs 0.069; J = 16: 0.008 vs 0.160; J = 29: 0.046 vs 1.28). So the published K₂ near-gap values for `J ≳ 17` are accurate (to
   ≲ 0.05 `N_s`) but **cannot be error-bounded from K₂ alone**; only direct rows bound them.
4. **The `ν < 50` window does not close sign-free at large J** under NEARGAP2's `G_low = 10⁴ c₂` (`W(ν < 50) ≈ 8.5·10⁻⁵ N_s/c₂` at
   J = 29). It is the one remaining model input at J = 29 (tolerance `3.2·10³ c₂` vs measured ≤ 116).

---

## 1. Methods (`code/highspin/`)

* **Write guard.** `hs_guard.py` (copy of `ng_guard` with root `highspin/cache/`): every `open(…, 'w'/'a'/'x'/'+')` and every
  `os.replace/rename/remove/unlink` outside `highspin/cache` raises. It registers itself as `ng_guard`, so importing `neargap2` modules
  does not install the neargap2 guard on top.
* **Direct points** (`hs_points.py` = `ng_points.py` with the cache moved): `hc_points` n1, adaptive re-centred contour, RQMC `2¹¹·2`
  points × 8 replicas, state and `Δ_ref = 2·10⁴` on common random numbers. Dense rows (`x = 0.025 … 0.95`) or the 15 Gauss–Kronrod nodes
  of an edge window (`gk15:a:b[:cap]`). A node whose state or reference record has replica error > 20 % is "failed"
  (`hs_core.failed_x`); this happens for `ν/Δ ≳ 0.94` (J = 29: every node with x ≥ 0.95; J = 16: x = 0.95).
* **Injection** (`hs_core.grid`): through `neargap2/ng_val`'s loader wrappers (= `hc_chain`'s `HC_INJ`): direct records (`M₂ = M₄ = 0`)
  override K₂ records at the same ν. **Bug found and worked around:** `p4_odd_analysis.grid` / `l9` rebuild ν from `x` rounded to
  6 decimals and match the reference record within `10⁻⁷ν`; a direct node at a non-round `x` (here the GK nodes) is then *silently
  dropped*. `hs_core` re-keys direct records to `round(x, 6)·m` (shift `≤ 5·10⁻⁷ m`, phase `≤ 10⁻⁵`). NEARGAP2's nodes were round.
* **Odd J** (`hs_odd.py`). `Γ` in the published odd-J convention (κ-curvature image subtracted, SMOOTH_BOUND §4; base = hybrid with its
  curvature). Bulk `x ∈ [0.05, 0.8)`: shape fit (`ee1000_chain.ShapeFit`, copied: `Γ = c(x)D_h + s(x)xD_h′/(2J+1) + r(x)`, degree 6) of
  the direct nodes (or of the K₂ nodes, for comparison); below x = 0.05 its continuation, and the window `ν < 50` also sign-free. Edge
  `x ∈ [0.8, Λ/m_row]`: **GK15 quadrature of the direct values at the GK nodes — no interpolation** — with the Gauss-7 subset as the
  quadrature error; for `m > m_row` the same nodes with the state's weights cut at `p = 1` (`ψ_k ∝ (1 − p)⁵`, so the integrand is C⁴
  there). Reproduces the published configuration exactly (J = 29 twist edge `+0.9477` vs `+0.948`).
* **Even J** (`hs_even.py`): `neargap2/ng_bind.run_state` unchanged (w2 conventions, own-mass row), every direct node injected,
  `E_mc = 3σ` over the replicas; `NG_NULOW = 25` variant uses the direct nodes at ν = 25, 50 for the window `ν ∈ [25, 50]`.
* **Budgets.** NEARGAP2 §1.2 conventions: ×3 on every model term (truncation, representation, quadrature), 3σ on Monte Carlo,
  sign-free `|Γ| ≤ 10⁴ c₂` on `ν < 50` (× the min-norm gauge factor for even J) and `G_edge = max(10⁴, 10 max|Γ|)` beyond the last node.
  Odd-J edge: a failed GK node is covered sign-free with `|Γ| ≤ 3A((1 − x_l)/(1 − x))²` (A = max|Γ| of the last three valid nodes: the
  measured `(1 − x)⁻²` growth of the edge amplitude, ×3); a GK node within `0.002` of the reference zero with `3 × max|Γ|` of its four
  nearest valid nodes.

---

## 2. J = 29: diagnosis and value

### 2.1 Why K₂ fails at (1029, 29), x = 0.8: the reference passes through a partial-wave zero [D + DV]

Every chain row is a **ratio**: the state record over the heavy reference `Δ_ref = 2·10⁴` at the same ν, multiplied back by the
reference's leading shape `D_h(x_ref)`, `x_ref = ν/m_ref` (`p4_odd_analysis.grid`, `l9_analysis.grid`). The reference sits in the Regge
region, `D_h(x_ref) ≈ J₀((2J+1)x_ref)`, whose first zero is at

    ν₀ ≈ 2.405 m_ref/(2J+1),     x₀(row Δ) = ν₀/m(Δ, J)      (exact Jacobi zero: hs_core.ref_zero_x).

Near ν₀ the reference's Laplace series is a series in `r₂/r₀` with `r₀ ∝ D_h(x_ref) → 0`, so the K₂ truncation of the cumulant form of
`1/(r₀ + r₂ + r₄)` (`pred`) fails there, however well the *state's* series converges. The direct method does not expand, but divides
by a small Monte-Carlo number: its error is amplified by `1/|D_h(x_ref)|`.

* **(1029, 29):** `x₀ = 0.8008`. The reference record at the K₂ node x = 0.8 has `|r₂/r₀| = 0.86` (C3H) and `0.077` (EE), against the
  smooth `≈ 19/ν` trend (`0.004–0.03`, C3H) and `≤ 0.001` (EE) at every other node. K₂ C3H there: `K = 0, 1, 2` =
  `−1.24·10⁵, −7.2·10⁴, −1.79·10⁴ c₂` (not converging) vs direct `+2.7·10³ ± 400` (κ-subtracted). The direct errors at x = 0.8 are ×4
  those of the neighbours, and at the GK node x = 0.8007 (`7·10⁻⁵` from `x₀`) the direct value is unusable (`C2 = +3.8·10⁴ ± 4.5·10³`).
* **Why the plain w2 odd chain gives −3.62:** w2's Hampel filter acts on even-J keys only, so the outlier stays, and the spline holds the
  last node's value (`≈ −1.8·10⁴ c₂`, non-oscillating) over the whole edge (`W[0.8, 0.972] = 4.5·10⁻⁴ N_s/c₂`).
* **Why +1.39 with the five direct edge nodes:** the spline through 5 nodes per ≈ 2 periods of an oscillation whose amplitude grows
  like `(1 − x)⁻²` (direct: `|Γ| ≈ 3·10³` at 0.82, `10⁴` at 0.92, `3·10⁴` at 0.94) aliases it into a smooth component and then holds the
  x = 0.9 value to the end. With 16 direct edge nodes the same spline gives `+1.44`, a shape fit `+0.72`, GK15 `+0.97`: **the spline
  value is an artefact; the w2 odd chain is not right at J = 29.**
* **The published value** (`+0.948`) came from the ee1000 shape chain **with the x = 0.8 node removed** (`sm_chainfix_ee1000_fodd_hampel_*`;
  the same chain without that step gives `−1.17`, `…fodd_K6_row_kappa`), shape-fitted to x ≤ 0.775 and continued with the exact `D_h`.
  That representation is close to right: its continuation of the shape fit beyond x = 0.8 is within `0.025 N_s` of the GK15 direct edge
  (twist edge `+0.022`, 1.04Λ `+0.009`), and its K₂ bulk is `0.046` below the direct bulk.

**Survey** (`cache/hs_refzero.txt`, every K₂ row of the 1000 chain, J = 17–40, 48; reference excess `|r₂/r₀|` minus its smooth trend):

| row | `x₀` | nearest node | excess | in the chain |
|---|---|---|---|---|
| (1600, 18) | 0.813 | 0.800 | 0.047 | dropped by Hampel (NEARGAP2: "the 1600 K₂ node at 0.8 is Hampel-dropped") |
| (1600, 20) | 0.733 | 0.750 | 0.056 | kept |
| (1600, 24) | 0.614 | 0.625 | 0.082 | dropped by Hampel |
| (1027, 27) | 0.862 | — (edge) | — | edge only |
| (1600, 28) | 0.528 | 0.550 | 0.041 | kept |
| **(1029, 29)** | **0.8008** | 0.800 | **0.83** | kept in w2; dropped in the published shape chain |
| **(1031, 31)** | **0.7477** | 0.750 | **0.54** | **kept in w2 and in the published shape chain** |
| **(1040, 32)** | **0.713** | 0.725 | **0.12** | **kept** |
| (1040, 40) | 0.572 | 0.575 | 0.44 | dropped (EE_DG1000_LARGEJ's "genuine outlier in both kinds") |
| (1048, 48) | 0.474 | 0.475 | — | (fw rows) RESULT_D_CLEANUP's "reference matrix nearly singular at x = 0.475" |

Published ee1000 shape chain with the reference-zero node removed (`cache/hs_ee1000_refzero*.txt`, K = 4 main / K = 6):
J = 31 twist edge `+1.200 → +1.023` (K6: `+1.208 → +1.019`), 1.05Λ `+1.135 → +1.034`; **J = 32 `1.0315Λ` `+0.904 → +0.997`**
(K6: `+0.785 → +0.945`), 1.05Λ `+0.937 → +0.994`; J = 40 unchanged. The large K-spread of J = 32 in EE_DG1000_LARGEJ's ensemble
(`+0.79 … +0.90`) was mostly this node.

**Remedy (for future chains).** Drop K₂ nodes with reference excess `|r₂/r₀| > 0.1` (`hs_core.refbad` — a property of the *reference*
record, decided before looking at the state's Γ, so not a post-hoc outlier cut), or move `Δ_ref` so that `ν₀` lies outside the row.

### 2.2 The value (`hs_odd.py 29`, `cache/hs_odd_J29_final.txt`)

Data: the K₂ row (34 nodes `x = 0.05 … 0.8`), **a full direct row** (34 nodes `x = 0.05 … 0.775`, 8 replicas), and **16 direct edge
nodes** (the 15 GK15 nodes of `[0.8, 0.97218]` plus 0.825 … 0.9); the four nodes with `x ≥ 0.95` failed (RQMC, `ν/Δ ≥ 0.95`).

*Direct vs K₂ in the bulk* (`cache/hs_nodes_J29_D1029.txt`): median `|Δ| = 50` (C2) / `66` (C3H) `c₂`, max `614 / 629` (x = 0.75),
1.4–1.7σ median, while the K₂ Laplace steps `|K2 − K1|` are `10²–3.5·10³ c₂`. At ν = 51–77 the direct `|Γ| ≤ 116 c₂`.

| `m/Λ` | published | K₂ bulk + GK edge | **direct bulk + GK edge** (± jk) | G7 | direct bulk + spline edge | K₂ bulk, K = 1 |
|---|---|---|---|---|---|---|
| **1.02862 (twist)** | +0.948 | +0.926 | **+0.972 ± 0.018** | +0.935 | +1.441 | −0.350 |
| 1.04 | +0.923 | +0.914 | **+0.956 ± 0.017** | +0.925 | +1.279 | −0.242 |
| 1.05 | +0.941 | +0.935 | **+0.977 ± 0.017** | +0.955 | +1.215 | −0.185 |
| 1.06 | +0.958 | +0.953 | **+0.995 ± 0.017** | +0.981 | +1.165 | −0.136 |
| 1.075 | +0.949 | +0.947 | **+0.986 ± 0.017** | +0.981 | +1.080 | −0.061 |
| 1.10 | +0.958 | +0.957 | **+0.995 ± 0.016** | +0.995 | +1.024 | +0.028 |

Budget (λ-level, ×3 on model terms, 3σ on MC; twist edge): MC `0.055`; GK15 − G7 `0.110`; missing GK nodes `0.48` (the reference-zero node
x = 0.8007: `0.20`; the failed `x ≥ 0.95` nodes: `0.28`); bulk representation (shape degree 4, 5, 7, two half-samples) `0.051`; total
**`0.70`, lower bound `+0.275`, `L = 1.39`**. At 1.04 / 1.05 / 1.06 / 1.075 / 1.10Λ: `L = 2.8 / 3.6 / 4.6 / 6.6 / 9.5`.
Window `ν < 50` (`W = 8.5·10⁻⁵ N_s/c₂` at the twist edge): at NEARGAP2's `G_low = 10⁴ c₂` it costs `0.85` (`L → 0.63`); at
`G = 10 ×` the largest direct `|Γ|` of the three lowest nodes (`1.16·10³`) `0.10` (`L = 1.22`). A failure needs a *uniform*
`|Γ| ≥ 3.2·10³ c₂` of the wrong sign on `ν < 50` (`7.9·10³ … 1.5·10⁴` at the higher masses).

**Verdict J = 29.** The state is **positive: `λ/N_s = +0.972 ± 0.018` at the twist edge, `≥ +0.956` on the computed masses**. The
correct chain is neither the w2 odd chain (`−3.62`: reference-zero K₂ node held over the edge) nor its spline with five direct nodes
(`+1.39`: aliasing); it is a shape/GK representation, which the published ee1000 shape chain approximates (`+0.948`, within 0.024 of
the direct value). Everything but the `ν < 50` window is error-bounded (`L = 1.4–9.5`); the window is the one model input (×28 above
the measured small-ν values, not ×1 above `10⁴`). **Grade [NE + M]** (was [M], "representation-dominated"); [NE±] under a measured
small-ν envelope.

---

## 3. The other edge states (sweep)

Direct n1 edge nodes on every row that carries one of the requested near-gap states (`cache/hs_odd_J{21,23,27}.txt`,
`cache/hs_edgeval_J{20,22,28,32}.txt`). Odd J: the 15 GK nodes of `[0.8, Λ/m_row]` (nodes above x = 0.945 skipped: weight
`≤ 2·10⁻⁶ N_s/c₂`), own rows `Δ = 1000 + J`, bulk K₂ (shape fit, as published). Even J: `x = 0.8, 0.825, …, 0.9` on the chain's rows
(not at the state's mass; 1/m interpolation, w2 spline chain, NEARGAP2's `ng_val` protocol), edge nodes only (`HS_IXMIN = 0.7`).

| J (rows) | `J²/m` | reference zero in row? | twist edge: published → with direct edge | other masses (shift) | budget / grade |
|---|---|---|---|---|---|
| 21 (1021) | 0.43 | no (x₀ = 1.12) | `+0.984 → +0.980 ± 0.007` (GK; G7 `+0.988`) | 1.03Λ `+0.869 → +0.868` (min); 1.05, 1.075, 1.1: `≤ 0.001` | edge: quadrature 0.022, missing 0.09; bulk K₂ step `λ_K2 − λ_K1 = 0.40` (×3 = 1.2) |
| 23 (1023) | 0.52 | no (x₀ = 1.02) | `+0.880 → +0.883 ± 0.0004` (G7 `+0.891`) | 1.03Λ `+0.903 → +0.905`; 1.05–1.1: `≤ 0.001` | edge: 0.023 + missing 0.31 (x = 0.941 failed); bulk step 0.49 |
| 27 (1027) | 0.71 | **in the edge**, x₀ = 0.862 (no GK node within 0.002) | `+0.956 → +0.959 ± 0.002` (G7 `+0.937`) | 1.05Λ `+0.928 → +0.928` (min); others `≤ 0.003` | edge: 0.064 + 0.14; bulk step 0.98 |
| 20 (800, 1600) | 0.50 / 0.25 | 1600: x₀ = 0.733 (K₂ 0.75 kept) | `+0.965 → +0.830 ± 0.138` (jk) | 1.03Λ `0.000`; 1.05Λ `−0.043 ± 0.055` | the 1600-row direct values carry `±1.7–9.5·10³ c₂` replica errors (vs `±50–180` on the 800 row, which agrees with them node by node to `≲ 100 c₂`): the twist-edge shift is not resolved beyond 1σ; ≥ `+0.41` at 3σ |
| 22 (336 only) | 1.45 | no | `+0.901 → +0.884 ± 0.001` | 1.03Λ `−0.013`; 1.05Λ `−0.007` | 0.875 failed (`ν/Δ = 0.87`); K₂ on this row is off by up to `900 c₂` at x = 0.8 (`J²/m = 1.45`); a 336 row for a 1000-gap state is itself a model (held in 1/m) |
| 28 (500, 1600) | 1.58 / 0.49 | 1600: x₀ = 0.528 (K₂ 0.55 kept) | `+0.990 → +0.974 ± 0.003` | 1.05, 1.075Λ `≤ 0.0004` | 1600: x = 0.825 failed |
| 32 (1040; + 336, 1600) | 0.99 | 1040: x₀ = 0.713 (K₂ 0.725 **kept**) | w2 chain `+0.913 → +1.006 ± 0.001` | 1.05Λ `+0.016`; 1.075Λ `−0.002` | favourable; the published ee1000 value (+0.904) also changes by `+0.09` when the reference-zero node is removed (§2.1) |
| 26, 30 | — | — | **not run: no K₂ rows exist** (the states are J-interpolated, inside `H_LJ(1000)`) | | [M] |
| 48 | 2.2 | x₀ = 0.474 | not run (the fw chain's structure basis needs `fw_frame.patch_conserved_clean` on the points and its own injection path) | | `H_edge` [M] |

**Reading.**
* **The edge is small wherever it was computed with a resolved representation.** Odd J = 21, 23, 27: the GK15 direct edge moves the
  published shape-chain values by `−0.004 … +0.003` (J = 29: `+0.024`, but its bulk moved too). Even J = 22, 28: `−0.017, −0.016`;
  J = 32: `+0.09` (favourable). J = 20's twist-edge shift `−0.13 ± 0.14` is Monte-Carlo-limited (1600-row direct noise).
* **Every computed state stays positive**; the lowest is J = 21 at 1.03Λ (`+0.868`) and J = 20 at its twist edge (`+0.83 ± 0.14`).
  J = 21 at 1.03Λ lies between the published mass grid points (`m_gap = 1.0208`, then 1.05): it is below the published J = 21
  minimum (`+0.939` at 1.1Λ) because the hybrid itself dips there (`+0.883`); an `H_cont` observation, not an edge effect.
* **Grades.** Edge: [M] → **[NE + M]** for J = 20, 21, 22, 23, 27, 28, 32 (the representation inside the window is now direct
  quadrature for odd J; spline between 5 nodes for even J). The states themselves stay **[NE + M]**: the bulk is K₂ and its λ-level
  Laplace step (0.40–0.98 for J = 21–27) does not close a ×3 budget. With the calibration of §4 (actual K₂ error 4–5 % of the step,
  three direct rows), the ×3-calibrated truncation is `0.06 / 0.07 / 0.15` and the twist-edge leverages would be ≈ `6 / 2 / 2`
  (J = 21 / 23 / 27) — an [NE + M] statement, the M being the calibration.

---

## 4. J = 12, 16: closing the budgets with direct rows

The next Laplace order (`L1_DEG = 8`, M6 records) would have cost ≈ 3–7 CPU-h per row; direct n1 rows on the whole decisive window
cost 38 min (J = 12) and 55 min (J = 16) and remove the Laplace truncation altogether. So: **full direct rows** on the own-mass rows
`Δ = 1012, 1016`, `x = 0.025 … 0.95` (38 nodes; J = 16: x = 0.95 failed), injected at every node (`hs_even.py`).

| state | own-row K₂ (NEARGAP2) | direct row | sign-free budget: trunc / interp / `ν < 50` / edge / MC | `L` sign-free (λ-level) | with `ν ∈ [25, 50]` direct |
|---|---|---|---|---|---|
| J = 12, 1.01194 (twist) | +0.8497 (`L` 0.3 / 1.6) | **+0.8528 ± 0.0004** | 0 / 0.003 / 0.142 / 0.094 / 0.001 → 0.24 | **3.6** (3.6) | 6.2 |
| J = 12, 1.02 | +0.9979 | +1.0006 ± 0.0008 | 0.065 / 0.001 / 0.140 / 0.047 / 0.002 → 0.26 | 3.9 (5.1) | 6.5 |
| J = 12, 1.05 | +0.6857 | +0.6872 ± 0.0003 | 0.175 / 0.001 / 0.098 / 0.001 / 0.001 → 0.28 | 2.5 (6.2) | 3.4 |
| J = 16, 1.01589 (twist) | +0.9342 (`L` 0.1 / 1.0) | **+0.9402 ± 0.0012** | 0 / 0.013 / 0.317 / 0.155 / 0.004 → 0.49 | **1.9** (1.9) | 3.6 |
| J = 16, 1.02 | +0.7632 | +0.7681 ± 0.0010 | 0.064 / 0.008 / 0.283 / 0.096 / 0.003 → 0.45 | 1.7 (2.0) | 3.1 |
| J = 16, 1.03 | +0.7458 | +0.7509 ± 0.0009 | 0.201 / 0.004 / 0.267 / 0.051 / 0.003 → 0.53 | 1.4 (2.2) | 2.3 |
| J = 16, 1.05 | +0.9942 | +1.0009 ± 0.0009 | 0.535 / 0.002 / 0.281 / 0.015 / 0.003 → 0.84 | 1.2 (2.9) | 1.6 |

(The remaining truncation at `m > m_row` comes from the published 1600 K₂ row, which the 1/m interpolation mixes in.)

* **Direct vs K₂ node by node** (`cache/hs_nodes_J12_D1012.txt`, `…J16_D1016.txt`): median `|Δ| ≈ 50 c₂` (J = 12) / `≈ 100 c₂`
  (J = 16), up to `200 / 480`, typically 2–4σ: K₂ has real errors of that size, but they oscillate with x and cancel in λ (state shifts
  `+0.003 / +0.008`).
* **Small ν, measured:** direct `|Γ| ≤ 50 c₂` (J = 12) and `≤ 150 c₂` (J = 16) at ν = 25–76, against K₂ `90–400`. So `G_low = 10⁴` is
  ×70–200 above the measurement; it is still the dominant budget term.
* **Edge:** direct values to x = 0.925 (`|Γ| ≤ 2.7·10³` J = 12); `x > 0.925` sign-free at `G_edge ≥ 2.7·10⁴`.

**Verdict.** **J = 12 and J = 16 near the gap are [NE±]** (were [NE + M]): sign-free `L ≥ 2.5` (J = 12) and `≥ 1.2` (J = 16, `≥ 1.9`
at the twist edge, `≥ 2.0` λ-level), with no Laplace truncation on the own row. J = 10, 14 were not run.

**Calibration of the K₂ truncation estimate.** `λ_K2 − λ_K1` vs `λ_direct − λ_K2` at the twist edge: J = 12 `0.069` vs `+0.003`; J = 16
`0.160` vs `+0.008`; J = 29 `1.276` vs `+0.046` (GK edge on both). The actual K₂ error is **4–5 % of the last Laplace step**, with the
same (favourable) sign in all three. NEARGAP2's `3 × |K2 − K1|` is therefore 60–80× conservative near the gap. [NE, three points]

---

## 5. H_LJ(1000): what is model, precisely

No compute route to [NE±] fits the budget (and the engine cannot decide `H_win` at `J²/Δ ≳ 9`, HLJ_CLOSE §5). The grade stays
**[D-sketch + M]**. Its content, item by item (RESULT_D_STATEMENT §5.2, HWIN_THIMBLE):

| piece | statement | grade | what would make it [NE±] |
|---|---|---|---|
| region | twist-allowed `m ≤ 1.3Λ`, even `J ≥ 50`, even `J ∈ {26, 30, 34, 38, 42, 44, 46}` (no rows: J-interpolated), odd `J ≥ 33` | — | — |
| `H_AG` | on `ν ∈ [ν_B, Λ]` the beyond-model correction is two travelling waves in class AG | [D] to all orders of the saddle expansion (J = 50–344); [NE + M] beyond all orders | direct rows at 2–3 large J (cost ∝ J: ≥ 3 CPU-h per row at J = 100) |
| `H_win` (refined) | `|Γ − κ − c₀D_h| ≤ T_B` on `ν < ν_B` (`T_B ≈ 10³ c₂`) | the `D_h`-class step [D-sketch]; data only at `ν ≥ 13` (J = 64, 100) [NE]; **`ν < 10` unmeasurable at any J**; sign-free form violated at (1100, 100) | a method below the fold (thimbles fail the scalar check, HYPOTHESES_CLOSE2) |
| `H_edge-amp` | edge amplitude `≲ 3·10⁶ / 9·10⁶` (J = 48 / 64) vs the scalar envelope | [M], margin ×3 (J = 48), ≥ 19 (J ≥ 64); graviton `x > 0.8` unmeasured at large J | GK edge nodes as in §2 (≈ 45 min per row at J = 29; more at large J, and n1 fails for `ν/Δ ≳ 0.94`) |
| state values | κ-uniform model `≥ +0.90` | [M] | — |

**What this note adds to it.** (i) The J = 29 value with a direct row is within `+0.02 … +0.05` of the published shape-chain value,
i.e. the shape representation of the oscillating Γ (the basis of the `H_AG` / κ-model arguments) is accurate at `J²/m = 0.8` [NE];
(ii) the edge amplitude grows like `(1 − x)⁻²` at J = 29 (`3·10³ → 3·10⁴ c₂` from x = 0.82 to 0.94) — this is the behaviour
`H_edge-amp` must bound at large J, and it is weighted by `ψ ∝ (1 − p)⁵`, so it stays harmless (`0.28 N_s` sign-free at J = 29);
(iii) the reference-zero mechanism (§2.1) affects large-J rows in the bulk (`x₀ = 0.47` at J = 48, `0.35` at J = 64, `0.22` at
J = 100 on Δ ≈ 1000 + J rows): any future large-J K₂ row must drop that node. (iv) The K₂ calibration (§4) suggests the published
K₂ values at `J²/m ≲ 1` are accurate to `≲ 0.05 N_s`; it does not reach `J²/m ≳ 1.5`, where K₂ diverges (EE_DG1000_LARGEJ).

---

## 6. Revised positivity inventory for Result D at Δ_g = 1000

| region / state | grade before (NEARGAP2 / POS_ERRORS) | grade now | basis |
|---|---|---|---|
| flat limit | [certified] | [certified] | CERT_CLOSE |
| `m ≥ 1.3Λ`, `b̂ < 5` | [NE±] (`L ≥ 8.7`) | unchanged | POS_ERRORS §2.4 |
| binding J = 5 `+−` near gap | [NE±] (`L = 78`) | unchanged | NEARGAP2 §1.3 |
| J = 0, 2, 4, 6 near gap | [NE±] | unchanged | NEARGAP2 |
| J = 8 near gap | [NE±] (λ-level, `L ≥ 4.3`) | unchanged | NEARGAP2 |
| **J = 12 near gap (to 1.05Λ)** | [NE + M] | **[NE±]**, `L ≥ 2.5` sign-free | §4 |
| **J = 16 near gap (to 1.05Λ)** | [NE + M] | **[NE±]**, `L ≥ 1.2` sign-free (`≥ 1.9` twist, `≥ 2.0` λ-level) | §4 |
| J = 7, 9–11, 13–15 near gap | [NE + M] (J = 15 own row) | unchanged | NEARGAP2; K₂ truncation at λ level ≈ 4–5 % of the step (§4 calibration) |
| J = 17–28 near gap (K₂ rows) | [NE + M]; `H_edge` [M] except J = 18, 19, 24, 25 | [NE + M]; **`H_edge` now computed for J = 20, 21, 22, 23, 27, 28** (and 32) [NE + M] | K₂ not error-bounded at the λ level (steps 0.3–1.1 `N_s`), accurate to ≲ 5 % of the step where checked; J = 26 has no row ([M]) |
| **J = 29 near gap** | [M] (representation-dominated) | **[NE + M]**: `+0.972 ± 0.018`, `L = 1.4–9.5` with the `ν < 50` window as the one model input (tolerance `≥ 3.2·10³ c₂`, measured ≤ 116) | §2 |
| J = 31, 32 near gap | [NE + M] (`+0.988`, `+0.904`) | [NE + M]; published values carried a reference-zero K₂ node; without it `+1.02`, `+0.997` (favourable) | §2.1 |
| J = 36, 40, 48 | [NE + M]; J = 48 edge [M] | unchanged | — |
| even J = 26, 30, 34–46 (no rows), even ≥ 50, odd ≥ 33 | [M] / [D-sketch + M] (`H_LJ(1000)`) | unchanged | §5 |
| `m ≥ 1.3Λ`, `b̂ ≥ 5` | [NE + M] / [M] (`H_R`) | unchanged | POS_ERRORS |

**Bottom line at Δ_g = 1000.** Error-bounded ([NE±] or better): flat limit, `b̂ < 5` at `m ≥ 1.3Λ`, and near the gap J = 0, 2, 4, 5, 6, 8,
**12, 16**. Computed with explicit model inputs ([NE + M]): the other near-gap J ≤ 32, now including **J = 29** (one model input: the
`ν < 50` window). Model only: `H_LJ(1000)` and `H_R` at large `b̂`. Nothing computed is negative.

---

## 7. Suggested corrections to existing documents (not applied)

1. **NEARGAP2 §0 / §1.6 and POS_ERRORS' update box:** "J = 29: not resolved … with the direct nodes +1.39" → "J = 29 resolved
   (HIGHSPIN §2): the w2 `−3.62` is a K₂ node at the reference's partial-wave zero (x = 0.8008) held over the edge; the `+1.39` is spline
   aliasing; with a direct bulk row and GK15 direct edge, `+0.972 ± 0.018`, [NE + M] (`ν < 50` window only)."
2. **NEARGAP2 §0, H_M rows J = 12, 16:** [NE + M] → [NE±] (direct rows; `+0.8528`, `+0.9402` at the twist edges).
3. **EE_DG1000_LARGEJ §4 / RESULT_D_STATEMENT §5.2, J = 31, 32:** the published near-gap values include a K₂ node at the reference
   zero (J = 31: x = 0.75; J = 32, row 1040: x = 0.725); without it J = 31 `+1.208 → +1.02`, J = 32 `+0.904 → +0.997`; the J = 32
   K-ensemble spread was mostly this node. EE_DG1000_LARGEJ's J = 40 "genuine outlier" and RESULT_D_CLEANUP's J = 48 "nearly singular
   reference at x = 0.475" are the same mechanism.
4. **The published odd-J chain description** should state that the x = 0.8 node of (1029, 29) was removed (the "hampel" variant);
   without it the same chain gives `−1.17`.
5. **NEARGAP2 §1.2 truncation rule:** `3|λ_K2 − λ_K1|` is 60–80× conservative near the gap (actual K₂ error 4–5 % of the step at
   J = 12, 16, 29); it is a valid bound but should not be read as an estimate.
6. **POS_ERRORS §3.2 "model-only states" list:** remove "the edge window of … odd J = 19–29" for J = 29 (computed; and see §3 here for
   21, 23, 27).
7. **Code (not changed):** `p4_odd_analysis.grid` / `l9_analysis.grid` drop any node whose `x` is not a 6-decimal number (reference
   matched within `10⁻⁷ν` after rounding x); w2's Hampel filter does not act on odd-J rows.

## 8. Reproduce
```bash
cd ads4/code/highspin          # every compute command under: ulimit -v 3000000; OMP_NUM_THREADS=1 (hs_lane.sh / hs_run.sh do this)
python3 hs_points.py 29 1029 seed                                    # copy NEARGAP2's direct points (likewise 18, 19, 25, 31, 17)
python3 hs_points.py 29 1029 gk15:0.8:0.972178                       # GK15 edge nodes (~50 min)
python3 hs_points.py 29 1029 0.05,0.0625,...,0.775                   # direct bulk row (cache/lane_A.log has the exact list; ~90 min)
python3 hs_odd.py 29 1.04,1.05,1.06,1.075,1.1                        # -> cache/hs_odd_J29.txt/.json (final copy: hs_odd_J29_final.*)
python3 hs_refzero.py; EE1000_JEVEN=32,40 EE1000_JODD=29,31 EE1000_K=6 python3 hs_ee1000.py refzero   # §2.1
python3 hs_points.py 12 1012 0.025,0.05,...,0.95; python3 hs_even.py 12 1012 1.0119407,1.02,1.05 J12    # §4 (also 16 1016)
NG_NULOW=25 python3 hs_even.py 12 1012 1.0119407,1.02,1.05 J12_nu25; python3 hs_nodes.py 12 1012
python3 hs_points.py 27 1027 gk15:0.8:0.974022:0.945; python3 hs_odd.py 27         # §3 (21, 23 likewise)
python3 hs_points.py 20 1600 0.8,0.825,0.85,0.875,0.9; python3 hs_edgeval.py 20 800,1600 1.0198198,1.03,1.05   # §3
python3 hs_edgeval.py 28 500,1600 1.027654,1.05,1.075; python3 hs_edgeval.py 32 1040 1.0315315,1.05,1.075
HS_IXMIN=0.7 python3 hs_edgeval.py 22 336 1.0217776,1.03,1.05      # edge nodes only (HC's small-nu 336 nodes are in the copied file)
```
Exact job lists and times: `cache/lane_A.log`, `cache/lane_B.log`.

## 9. Caveats
1. The direct n1 integrals are the project's reference standard, validated elsewhere [NE]; they are not certified.
2. `G_low = 10⁴ c₂` (and the measured-envelope alternative) on `ν < 50` remain hypotheses; at J = 29 the window is the binding input.
3. The J = 29 edge budget uses an extrapolated `(1 − x)⁻²` envelope for the four failed GK nodes (weight `≤ 1.7·10⁻⁶ N_s/c₂`).
4. Only the masses listed were evaluated; between them, `H_cont` as before.
5. The K₂ calibration (4–5 % of the step) rests on three states.
