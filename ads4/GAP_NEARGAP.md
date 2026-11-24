# GAP_NEARGAP: the error-budgeted near-gap states of Ω^EE at Δ_g = 2000 (and 5000) on own-mass K₂ rows; chart-A test; batch 4

*24 Nov 2026, gap-near agent. Unreviewed, not committed by this agent. GAP_SCALING.md §6.1, plan item 2 (and the preparation of item 3). No
existing file was edited; the one new file outside `code/gapnear/` is `ext/jobs_batch4.txt`. `GAP_LARGESPIN.md` / `code/gaplarge/` were
not touched.*
* *New code:* `ads4/code/gapnear/`.
  * `gn_guard.py`: process-local write guard (every write outside `gapnear/cache/` raises; it registers itself as `ng_guard`, `hs_guard`,
    `sn_guard`, `ef_guard`).
  * `gn_rows.py`: own-mass K₂ rows, odd J = `gapscale/gs_rows.py`, even J = `neargap2/ng_rows.py`'s even branch, both re-typed with every
    cache path moved to `gapnear/cache/rows/`.
  * `gn_bind.py`: NEARGAP2's budget (`neargap2/ng_bind.run_state`, unchanged, in-process) at Δ_g = 1000 / 2000 / 5000, plus the Γ(x)
    tables (`--hmono`).
  * `gn_chartA.py`, `gn_hsum.py`, `gn_table.py`: comparisons. `gn_refdiag.py`, `gn_stepnorm.py`: the K₂-step diagnostics of §2.4.
  * `gn_batch4.py`: writes `ext/jobs_batch4.txt`. `gn_pilot.py`: runs one batch-4 job locally (`ext/point.py`, unchanged).
  * `gn_lane.sh`, `gn_pause.sh`: compute lanes, and short analyses with one lane SIGSTOPped.
* *Caches:* `ads4/code/gapnear/cache/`; its `.gitignore` ignores everything but itself.
* *Read-only imports:* `neargap2/` (`ng_chain`, `ng_bind`), `spinfunc/` (`w1_rows`, `p3_odd_action`, `p2_shadow_odd`, `l9_points`,
  `l5_shadow_matrix`, `w0_ee_common`, `ee1000_common` → `w2_chain`, `q3_analysis`, `p4_odd_analysis`), `ext/point.py` → `hc_points`.
  Read-only data: the K₂ rows of `spinfunc/cache`, `neargap2/cache`, `gapscale/cache/rows`, and `gapscale/cache/weights_{2000,5000}.npz`.
* *Compute:* ≈ 7.8 CPU-h (two lanes, 04:28–≈08:00), `ulimit -v 3000000`, one thread, no bytecode written.

Grades: **[P]**, **[certified]**, **[D]**, **[DV]**, **[NE±]** (explicit error budget, every model input with a stated safety factor,
conclusion surviving it), **[NE]**, **[NE + M]**, **[M]**, **[I]**. `L = value / budget`. "Sign-free" and "λ-level" budgets as NEARGAP2
§1.2. `c₂ = m² ×` relative correction. Λ = √((Δ_g − 2)(Δ_g − 1)): 998.5 / 1998.5 / 4998.5.

---

## 0. Table and verdict

**Same protocol at every gap:** own-mass K₂ row at the state's twist edge `Δ = Δ_g + J` (light J = 0, 2: `Δ_L = 0.99Δ_g`), dense
`x = 0.025 … 0.9`, NEARGAP2's `ng_bind` budget unchanged (×3 on the K₂ truncation and interpolation; `|Γ| ≤ 10⁴ c₂` sign-free on
`ν < 50` with the min-norm gauge factor of the gap's own weight; `G_edge = max(10⁴, 10 × max|Γ|)` beyond `x = 0.9`). "L" is the
sign-free leverage, "(λ)" the λ-level one. The 1000 *graded* column is the value that carries the 1000 grade (for J ≥ 7 it rests on
direct rows, not on K₂); the "1000 same protocol" column is this note's rerun of `ng_bind` at 1000 (it reproduces NEARGAP2 to every
digit, §1.3), which is the like-for-like reference for chart A.

| state (J, m/Λ) | Δ_g = 1000 graded: value / L | 1000, same protocol (K₂): value / L (λ) | **2000** own row K₂: value / L (λ) | 5000: value / L (λ) | grade at 2000 |
|---|---|---|---|---|---|
| 0 light, 0.98999 (Δ_L = 0.99Δ_g) | +1.0000 / 4362 [NE±] | +1.0000 / 4362 (4634) | **+1.0000 / 40612** (46824) | — | **[NE±]** |
| 2 light, 0.98998 | +1.0000 / 913 [NE±] | +1.0000 / 913 (1372) | **+1.0000 / 5712** (12640) | — | **[NE±]** |
| 4, 1.00200 (2000 twist edge) | — (seam) | — | **+0.2021 / 103** (853) | — | **[NE±]** |
| 4, 1.00400 (1000 twist edge) / 1.01 / 1.02 | +0.2041 / 22 (λ 80) [NE±] | +0.2041 / 21.9 (79.9); 23.2; 25.4 | +0.2041 / **105.5** (925); 112; 123 | — | [NE±] |
| 5, 1.00250 (2000 twist edge) | — (seam) | — | **+0.2004 / 320** (565) (GAP_SCALING) | — | **[NE±]** |
| 5, 1.00499 / 1.01 / 1.02 / 1.05 | +0.2009 / 78 (λ 130) [NE±] | +0.2009 / 77.7 (130); 81; 85; 89 | +0.2009 / **334** (608); 352; 370; 394 | PENDING5000 | [NE±] |
| 6, 1.00300 (2000 twist edge) | — (seam) | — | **+0.2144 / 20.1** (236) | — | **[NE±]** |
| 6, 1.00599 / 1.01 / 1.02 | +0.2295 / 5.0 (λ 28) [NE±] | +0.2295 / 5.0 (27.8); 5.6; 7.3 | +0.2295 / **22.0** (279); 24.7; 32.5 | — | [NE±] |
| 7, 1.00350 (2000 twist edge) | — (seam) | — | **+0.2298 / 26.0** (93) | — | **[NE±]** |
| 7, 1.00698 / 1.01 / 1.02 / 1.05 | +0.2683 / 6.7 [NE±] (direct row 1007) | +0.2682 / 4.3 (7.6) (rows 800/1600, no 1007 K₂ row) | +0.2684 / **27.5** (103); 28.8; 33.5; 53 | — | [NE±] |
| 8, 1.00399 (2000 twist edge) | — (seam) | — | **+0.2849 / 3.9** (19.3) | — | **[NE±]** (sign-free) |
| 8, 1.00798 / 1.02 / 1.05 | +0.4297 / 1.2 (λ 4.6) [NE±, λ-level only] | +0.4297 / 1.2 (4.6); 1.1 (4.3); 1.8 (8.4) | +0.4292 / **5.2** (24.1); 4.5 (22.4); 7.9 (51) | — | **[NE±] sign-free** (was λ-level only) |
| 12, 1.00599 (2000 twist edge) | — (seam) | — | **+0.4094 / 0.7 (λ 5.4)** | — | **[NE±] λ-level** |
| 12, 1.01194 / 1.02 / 1.05 | +0.8528 / 7.4 [NE±] (direct row) | +0.8497 / 0.3 (1.6); 0.3 (2.2); 0.4 (3.1) | +0.8496 / 1.1 (**8.0**); 1.3 (11.6); 1.6 (18.2) | PENDING5012 | [NE±] λ-level (K₂ alone; 1000 needed a direct row) |
| 16, 1.00797 (2000 twist edge) | — (seam) | — | **+0.9993 / 0.4 (λ 3.7)** | — | **[NE±] λ-level** |
| 16, 1.01589 / 1.02 / 1.03 / 1.05 | +0.9402 / 4.4 [NE±] (direct row) | +0.9342 / 0.1 (1.0); 0.1 (1.0); 0.1 (1.1); 0.1 (1.3) | +0.9422 / 0.4 (**4.4**); 0.4 (4.6); 0.5 (5.1); 0.6 (6.4) | — | [NE±] λ-level |
| 21, 1.01045 (2000 twist edge) | — (seam) | — | **+1.0043 / 0.1 (λ 0.9)** | — | **[NE + M]** (calibrated L 3.0) |
| 21, 1.02080 / 1.03 / 1.05 / 1.075 / 1.10 | +0.998 ± 0.020 / 1.20 [NE±] (direct bulk) | +0.967 / 0.07 (0.17); … (0.24–0.75) | +1.0047 / 0.2 (**1.1**); 0.2 (1.2); 0.2 (1.4); 0.3 (1.7); 0.3 (2.0) | — | [NE±] λ-level, thin (L 1.1–2.0) |
| 27, 1.01342 (2000 twist edge) | — (seam) | — | **+1.0121 / 0.0 (λ 0.4)** | — | **[NE + M]** (calibrated L 0.7) |
| 27, 1.02667 / 1.03 / 1.04 / 1.05 / 1.075 / 1.10 | +0.976 ± 0.019 / 1.64 [NE±] (direct bulk) | +0.996 / 0.02 (0.09); … (0.10–0.41) | +1.0305 / 0.1 (**0.6**); 0.7; 0.8; 1.0; 1.2; 1.3 | — | [NE + M] to 1.05Λ; [NE±] λ-level at 1.075–1.10Λ (thin) |
| 29, 2000 twist edge / 1.02862 / 1.04–1.10 | +0.972 ± 0.018 / 1.24 signed (×3: 0.94) [NE±, thin] | +0.912 / 0.01 (0.05) (x = 0.8 reference-zero node dropped) | PENDING29 | — | PENDING29G |

(Masses other than the twist edges at 1000 also have graded values: e.g. J = 12 at 1.02 / 1.05Λ L = 7.5 / 3.6, J = 16 at
1.02–1.05Λ L = 3.7 / 2.6 / 1.7, J = 21 at 1.03–1.10Λ L = 1.7–6.2, J = 27 at 1.03–1.10Λ L = 1.93–3.21, J = 29 at 1.04–1.10Λ L = 2.27–5.86;
SMALLNU_WINDOW, EXT_FOLD.)

**Verdict.**
1. **No state is negative at Δ_g = 2000 (or 5000).** The lowest value is the binding state, `+0.2004` at its 2000 twist edge. The lowest
   at J ≥ 8 is J = 12 at its 2000 twist edge, `+0.409` (a dip of the hybrid itself, `+0.409`, not a correction; λ-level `L = 5.4`).
2. **J ≤ 8 and the light states: [NE±] at 2000 with ×4–9 larger leverage than at 1000**, on K₂ alone (sign-free `L = 3.9 … 4·10⁴`).
   J = 8 is now sign-free [NE±] (`L = 3.9–7.9`; it was λ-level only at 1000). J = 7 needs no direct row any more (`L = 26–53` on its K₂ row,
   against 4.3 / 6.7 at 1000).
3. **J = 12, 16: [NE±] on the λ-level budget from K₂ alone** (`L = 4.4–18`), as J = 8 was at 1000; sign-free they do not close
   (`L = 0.4–1.6`). At 1000 they needed a direct row for any budget. The direct rows of batch 4 would make them sign-free [NE±].
4. **J = 21, 27 (and 29): K₂ alone does not close at the twist edges** (λ-level `L = 0.4–1.1`), as at 1000; these stay **[NE + M]** until
   the batch-4 direct rows are folded in (at 1000 they were [NE±] only through direct rows). J = 21 at ≥ 1.02Λ and J = 27 at ≥ 1.075Λ
   close thinly on λ-level (`L = 1.1–2.0`).
5. **Chart A holds on every budget term for J ≤ 16, exactly**: `E_trunc`, the λ-level K₂ step and `E_edge` scale by `0.25 ± 0.01`
   (`(Λ₁₀₀₀/Λ₂₀₀₀)² = 0.2496`), `E_low` by `0.065` (≈ Λ⁻⁴). **It fails for J = 21, 27 (×0.43–0.70)** — a K₁-order artefact of the
   heavy reference, not a physics trend: the reference's leading partial wave has its first zero *inside* every 2000 row with
   J ≥ 16 (x₀ = 0.72 / 0.56 / 0.44 / 0.41 at J = 16 / 21 / 27 / 29), and the K = 1 ratio degrades around it (§2.4). The K₂ values Γ(x)
   themselves satisfy H_mono-A to `≤ 30 c₂` (J ≤ 6), `≤ 110` (J = 8, 12), `≤ 170 / 300 / 330` (J = 21 / 27 / 16, x ≤ 0.8).
6. **The seam** (2000 states between the 2000 and 1000 twist edges) is closed by direct computation for every J of the set (the 2000
   twist edges in the table). The other seam spins (J = 9–11, 13–15, 17–20, 22–26, 28, 30–32) need rows; §5 lists them, and batch 4 carries
   an optional seam group for J = 20, 23.
7. **Batch 4** (`ext/jobs_batch4.txt`): **247 main n1 jobs (J = 12, 16, 21, 27, 29) + 98 optional seam jobs (J = 20, 23) = 345**.
   Pilot here: 3.2–3.6 min per node at Δ = 2029 (J = 29), so ≈ 13–15 CPU-h main (≈ 0.9 h on 16 workers), ≈ 20 CPU-h with the seam
   group (≈ 1.3 h); worst case at 6 min per node 35 CPU-h (≈ 2.2 h).

**Negative findings (reported prominently).**
* **N1. The Δ_ref = 2·10⁴ reference is too close at Δ_g = 2000 for J ≥ 16.** Its leading partial wave vanishes at
  `x₀ = 2.405…·m_ref/((2J+1) m_row)`, inside the row for J ≥ 16 at 2000 (J ≥ 5 at 5000: x₀ = 0.874 on the 5005 row). The K₂ nodes nearest
  x₀ are flagged by HIGHSPIN's reference-excess rule (`refbad`, applied here) or removed by the even chain's Hampel filter (J = 16: x = 0.725),
  and the K = 1 truncation of the ratio is visibly worse over a band around x₀ (§2.4). Consequence: at J = 21, 27 the K₂ step — the
  protocol's truncation estimate — *does not shrink* like Λ⁻², and the K₂-only budgets at 2000 are pessimistic by a factor ≈ 2.
  Remedy for any further 2000 / 5000 rows: `Δ_ref ≥ 5·10⁴` (x₀ > 1 for J ≤ 29 at 2000); it needs the module constant `DREF` of
  `p4_odd_analysis` / `l9_analysis` / `hc_points` made a parameter. **Batch 4 keeps Δ_ref = 2·10⁴** (the direct rows divide by the
  reference without expanding it; only the node nearest x₀ is unusable, and two flanking nodes are added).
* **N2. The δ = λ_K2 − λ_hyb test of chart A is not decisive beyond J = 8**: at J = 12–27 the beyond-hybrid shift of the K₂ chain is
  smaller than its own K₁ → K₂ step at both gaps (e.g. J = 16: `δ = −4.6·10⁻³ → +3.3·10⁻³` with steps `0.16 → 0.05`), so its "ratio" is
  noise. Chart A is tested on the budget terms and on Γ(x), not on δ, for J ≥ 12.
* **N3.** J = 21, 27 (and J = 29) twist edges are **not** error-bounded at 2000 by this note: K₂ alone does not close, as at 1000.

---

## 1. Protocol

### 1.1 Rows (`gn_rows.py`)
* K₂ crystal points (`w1_rows` engines unchanged, `L1_DEG = 6`): odd J → parity-odd `p3` integrand with `E` (`w0_ee_common`), 1×1 `+−`;
  even J → `l9_points` with the w0 EE integrand, 2×2. Kinds EE (k = 2) and C3H (k = 3); state and `Δ_ref = 2·10⁴` at the same ν.
* Rows: `Δ = 2000 + J` for J = 4, 6, 7, 8, 12, 16, 21, 27, 29 (J = 5: GAP_SCALING's `2005` row, read-only), `Δ = 1980` for the light
  J = 0, 2, `Δ = 5005`, `5012` at 5000; 36 nodes `x = 0.025 … 0.9`, 144 points per row. Shadows (`⟨TTO⟩` normalisation) were
  computed at the new (Δ, J) with the engines' defaults and stored in `cache/rows/gn_*shadow*.pkl`.
* Cost per point (one core): 5 s (J = 4) … 9 s (J = 7) … 11 s (J = 16), 23 s (J = 21), 30 s (J = 27), 32 s (J = 29), 7 s (J = 5 at 5005).

### 1.2 Budget (`gn_bind.py` → `neargap2/ng_bind.run_state`, unchanged)
NEARGAP2 §1.2 exactly: central value with the two outer windows removed; `E_trunc = 3 Σ|w||Γ_K2 − Γ_K1|` (sign-free, own-row range);
`E_trunc^λ = 3|λ_K2 − λ_K1|`; `E_interp` (spline vs linear, every-other-node); `E_low = 10⁴ c₂ × W(ν < 50) ×` gauge factor;
`E_edge = G_edge × W(x > 0.9)`; base: hybrid and the certified floor `N_f/5 − Lemma C`. Redirections, all in-process:
* the row loader also reads `gapnear/cache/rows` and `gapscale/cache/rows` (hidden for the Δ_g = 1000 reruns);
* the gauge factor at 2000 / 5000 uses the min-norm weight `gapscale/cache/weights_<Dg>.npz`, as `gapscale/gs_bind.py`;
* `q3_analysis.odd_grid / even_grid`: identical copies with the row cut `D ≤ 3200` raised to `6000`, and HIGHSPIN's reference-zero rule
  (`refbad`: drop a K₂ node whose reference record has excess `|r₂/r₀| > 0.1`, decided on the reference before Γ is looked at) applied to
  the new rows and to the 1000 ee1000 rows of odd J ≥ 17. It dropped exactly one node: x = 0.8 of (1029, 29), the node HIGHSPIN identified.
  (On the 2000 rows no node exceeded 0.1; the even chain's Hampel filter had already removed x = 0.725 of (2016, 16), excess 0.31.)

### 1.3 Reproduction at 1000 (`cache/gn_bind_base1000.json`)
With the new rows hidden, `gn_bind` reproduces NEARGAP2's table to every printed digit (J = 0, 2, 4, 5, 6, 8, 12, 16; e.g. J = 5
`+0.20088`, `L = 77.7 / 130.2`; J = 16 `L = 0.1 / 1.0`), EXT_FOLD's "published rows alone, `L = 4.3`" for J = 7, and GAP_SCALING's 2000
binding state (`+0.20041`, `L = 319.8 / 564.6`). For J = 21, 27, 29 the same chain on the ee1000 rows gives K₂ values `+0.967 / +0.996 /
+0.912` at the twist edges with `L ≤ 0.2` — the [NE + M] situation that the direct rows resolved at 1000.

### 1.4 Grading
As NEARGAP2 / HIGHSPIN: [NE±] if the sign-free budget closes, "[NE±] λ-level" if only the λ-level one does (J = 8's 1000 grade);
[NE + M] otherwise. "Calibrated" L (SMALLNU / EXT_FOLD: actual K₂ error 0.6–8 % of the step, budgeted as 15 % of the step) is quoted as
[NE + M] information only. The K₂ = crystal-integral input is checked against direct RQMC at 1000 (J = 5, 6, 7, 8, 12, 16, 21–31), not at
2000 (batch 4 does that for J = 12–29).

---

## 2. Results

### 2.1 J ≤ 8 and the light states
All close sign-free with large margins (table). Per state the change from 1000 is a pure scaling of the budget (§3); the values move by
`≤ 5·10⁻⁴` at fixed m/Λ (J = 8: `+0.42973 → +0.42921`). At the 2000 twist edges (the seam, §5) every state is [NE±]: J = 4 / 5 / 6 / 7 / 8
`L = 103 / 320 / 20 / 26 / 3.9`. Certified-floor route at the 2000 twist edges (`(N_f/5 − Lemma C)/N_s`): `+0.19993 … +0.19999` (J = 4–7).

### 2.2 J = 12, 16
* K₂ budgets at 2000 (twist edge of 1000, chart A): J = 12 `E_trunc 0.75`, `E_trunc^λ 0.053`, `E_low 0.009`, `E_edge 0.045`; J = 16
  `2.04 / 0.15 / 0.021 / 0.037`. λ-level `L = 8.0 / 4.4` (1000: 1.6 / 1.0). Sign-free `L = 1.1 / 0.4` (1000: 0.3 / 0.1).
* Own-row K₂ vs the 1000 direct value: J = 12 `+0.8496` (1000 direct `+0.8528`); J = 16 `+0.9422` (1000 direct `+0.9402`, 1000 K₂
  `+0.9342`).
* 2000 twist edges: J = 12 `+0.4094` (hybrid `+0.4090`), λ-level `L = 5.4`; J = 16 `+0.9993`, `L = 3.7`.
* Calibrated (15 % of the step) `L = 8.6–69` [NE + M information].

### 2.3 J = 21, 27 (and 29)
* J = 21 (row 2021, `J²/m = 0.22`): twist edge of 1000 `+1.0047`, λ-level budget `0.94` (`E_trunc^λ 0.78`, `E_edge 0.11`), `L = 1.1`;
  1.03 / 1.05 / 1.075 / 1.10Λ `L = 1.2 / 1.4 / 1.7 / 2.0`; 2000 twist edge `+1.0043`, `L = 0.9`. Calibrated `L = 3.0–16`.
* J = 27 (row 2027, `J²/m = 0.36`): twist edge of 1000 `+1.0305`, `L = 0.6` (λ-level `0.97 + E_edge 0.52`); 1.03–1.10Λ `L = 0.7–1.3`;
  2000 twist edge `+1.0121`, `L = 0.4` (`E_edge = 1.22`: the last four nodes reach `|Γ| ≈ 1.1·10⁴ c₂`, so `G_edge = 1.1·10⁵`). Calibrated
  `L = 0.7–11`.
* J = 29: PENDING29TEXT
* Compared with the published-row model (the 1000 + J row held): the own row moves the states by `+0.013` (J = 21) / `+0.031` (J = 27) at
  the 1000 twist-edge mass — favourable, as every direct-vs-K₂ shift at 1000 was.

### 2.4 Why the K₂ step does not shrink at J ≥ 21 (`gn_refdiag.py`, `gn_stepnorm.py`)
* The raw K = 2 term of the ratio (`a₂ + a₁e₂ + a₀e₄`, c₂ units) is the same on the 1021 and 2021 rows node by node (e.g. EE x = 0.4:
  `3318 / 3042`; C3H `3237 / 3702`): the *state's* Laplace series is Λ-independent in c₂, as chart A requires.
* The reference's own second-order term `e₂` is not: at x = 0.6 (C3H) `3.0·10⁴ → 1.45·10⁵ c₂` from row 1021 to 2021 — the 2000 row
  probes the reference at twice the x_ref, next to its partial-wave zero (x₀ = 0.563).
* At the chain level the K₁ rows pick this up: the C3H K₂ − K₁ step of row 2021 is mostly a row renormalisation `c D_h` with `c = +8657 c₂`
  (row 1021: `−1158`); after removing it the rms step is `957` vs `739 c₂`. Same at J = 27 (`c = +9903` vs `−4695`; residual rms C2
  `3223` vs `1898`). At J = 16 (x₀ = 0.72, near the row's end) the split is identical on both rows.
* So the protocol's truncation estimate is inflated by the reference, not by the state; the K₂ values themselves agree across gaps
  (§3.2). With `Δ_ref = 5·10⁴` the zero would leave the rows (N1).

---

## 3. Chart-A check (`gn_chartA.py` → `cache/chartA.txt`; `gn_hsum.py` → `cache/hsum.txt`)

### 3.1 Budget terms at fixed (J, m/Λ), ratio 2000 / 1000 (expected `(Λ₁₀₀₀/Λ₂₀₀₀)² = 0.2496` for c₂-type corrections, `≈ 0.0623` for the
small-ν window, whose weight per c₂ falls as Λ⁻⁴)

| J | E_trunc | λ-level K₂ step | E_low | E_edge | δ = λ_K2 − λ_hyb | L: 1000 → 2000 (λ-level) |
|---|---|---|---|---|---|---|
| 0, 2 (light) | 0.223, 0.256 | 0.218, 0.243 | 0.065 | 0.250 | (|δ| ≤ 5·10⁻⁷) | 4362 → 40612; 913 → 5712 |
| 4 | 0.253 | 0.250–0.277 | 0.065 | 0.250 | 0.262 | 22 → 106 (80 → 925) |
| 5 | 0.256 | 0.252 | 0.065 | 0.250 | 0.19–0.21 | 78 → 334 (130 → 608) |
| 6 | 0.253 | 0.254 | 0.065 | 0.250 | 0.239 | 5.0 → 22 (28 → 279) |
| 7 | 0.26 (vs 800/1600 rows) | 0.252 | 0.065 | — (1000: no own row) | 0.19 | 4.3 → 27.5 (7.6 → 103) |
| 8 | 0.253 | 0.253 | 0.065 | 0.250 | 0.288 at the twist edge (small δ elsewhere) | 1.2 → 5.2 (4.6 → 24) |
| 12 | 0.252 | 0.253 | 0.065 | 0.252 | not resolved (N2) | 0.3 → 1.1 (1.6 → 8.0) |
| 16 | 0.247–0.251 | 0.28–0.31 | 0.065 | 0.280 | not resolved (N2) | 0.1 → 0.4 (1.0 → 4.4) |
| 21 | 0.61–0.70 | 0.48–0.64 | 0.065 | (1000 row ends at 0.8) | not resolved | 0.07 → 0.2 (0.17 → 1.1) |
| 27 | 0.43–0.49 | 0.43–0.50 | 0.065 | — | not resolved | 0.02 → 0.1 (0.09 → 0.6) |
| 29 | PENDING29CA | | | | | |

**Verdict on chart A.** For J ≤ 16 every budget term scales exactly as predicted: the c₂-type terms by 0.25 (to ±0.01; J = 16's step
0.28–0.31 shows the next O(J²/m) term), the small-ν window by 0.065 ≈ Λ⁻⁴, and where δ is resolved (J ≤ 8) δ by 0.19–0.29. Leverages
therefore grow ×4–9. **[NE] for chart A at J ≤ 16, on the graviton, Δ_g = 1000 → 2000.** For J = 21, 27 the K₂ budget terms scale only
by 0.43–0.70 (N1, §2.4): the protocol's own error estimate violates chart A through the reference, while the quantity it bounds (Γ)
satisfies it (§3.2). Chart A for the *values* at J ≥ 21 awaits the direct rows (batch 4).

### 3.2 H_mono-A on the rows: Γ(x) in c₂ at fixed (J, x), K₂ (`cache/gn_hmono_*.txt`)

| J | rows | max |Γ| on the 2000 row | max |Γ₂₀₀₀ − Γ₁₀₀₀| (x ≤ 0.8) | max |Γ₂₀₀₀ − Γ₁₆₀₀| | K₂ − K₁ step (max, 1000 / 2000 row) |
|---|---|---|---|---|---|
| 0 / 2 | 990, 1600, 1980 | 6 / 77 | 0.5 / 4.8 | 0.1 / 0.4 | 44 / 39; 236 / 231 |
| 4 | 1004, 1600, 2004 | 255 | 15 | 15 | 593 / 602 |
| 5 | 1005, 1600, 2005 | 41 | 8.6 | 1.2 | 350 / 358 |
| 6 | 1006, 1600, 2006 | 234 | 31 | 61 | 1232 / 1245 |
| 7 | 1600, 2007 | 198 | — | 21 | — / 1135 |
| 8 | 1008, 1600, 2008 | 370 | 50 | 3.4 | 2029 / 2048 |
| 12 | 1012, 1600, 2012 | 985 | 107 | 22 | 3761 / 3790 |
| 16 | 1016, 1600, 2016 | 1998 | 328 (1931 incl. x = 0.85–0.9) | 389 | 4251 / 4642 |
| 21 | 1021, 2021 | 2786 | 167 | — | 7218 / 5821 |
| 27 | 1027, 2027 | 10816 | 297 | — | 6904 / 10210 |
| 29 | PENDING29HM | | | | |

* **Γ converges to an m-independent profile at fixed (J, x)** — H_mono-A [NE] at every computed spin. The 1000 → 2000 differences are
  O(J²/m) (they are ≤ 5 % of max|Γ| at J ≤ 12 and the 1600 → 2000 differences are ×2–15 smaller), and the K₂ step per node in c₂ is the
  same on both rows for J ≤ 12 — which is why every K₂ budget term falls as Λ⁻².
* At J = 16 the two edge nodes x = 0.85, 0.9 differ by up to 1.9·10³ c₂ between 1016 and 2016 (a sign change at x = 0.85, both rows with
  K₂ steps ≥ 10³): the K₂ edge is not resolved on either row, so this is not a test of H_mono-A there; the direct rows of batch 4 reach
  x = 0.95.
* At J = 21, 27 the K₂ *values* agree to ≤ 170 / 300 c₂ over x ≤ 0.8 (e.g. J = 21 C3H, x = 0.8: `−1758 / −1741`), while the K₂ − K₁ step
  grows (§2.4).

### 3.3 5000
PENDING5000TEXT

---

## 4. Batch 4 (`ext/jobs_batch4.txt`, generated by `gn_batch4.py`)

**Content.** Direct n1 nodes on the 2000 own rows `Δ = 2000 + J`, on the x grids of the 1000 protocol (chart A), so that the 1000
analysis scripts (`highspin/hs_odd.py`, `hs_even.py`, `smallnu/sn_odd.py`, `sn_even.py`, through `extfold`-style wrappers with Dg = 2000)
can be rerun unchanged once Δ_g is parametrised (they hard-code `Dg = 1000` and `D = 1000 + J`; §6).

| group | J = 12, 16 (even) | J = 21, 27, 29 (odd) | source of the 1000 protocol |
|---|---|---|---|
| bulk | full direct row x = 0.025 … 0.95, 38 nodes | x = 0.05 … 0.15 (step 0.0125) + 0.175 … 0.775 (step 0.025), 34 nodes | HIGHSPIN §4; ext jobs A / SMALLNU §3.3 |
| blo | — | x = 0.025, 0.0375 (at 2000 the bulk's first node is at ν ≈ 100, not 50) | new (chart-A gap between window and bulk) |
| edge | (inside the row) | GK15 on `[0.8, Λ/m_row]` (`Λ/m = 0.9897 / 0.9868 / 0.9858`), the 10 nodes with x ≤ 0.945 | HIGHSPIN §2–3 (n1 fails for ν/Δ ≳ 0.94; the 5 skipped nodes are covered by HIGHSPIN's E_miss envelope) |
| rz | J = 16: x₀ ± 0.0125 | x₀ ± 0.0125 (x₀ = 0.563 / 0.436 / 0.406) | new (N1: the node nearest x₀ is unusable) |
| lo | n1 at ν = 20, 25, 32.5, 40, 47 | same | SMALLNU (n1 part only) |
| jobs | 43 (J = 12), 45 (J = 16) | 53 each | **247 main** |
| **optional seam group S** | J = 20: 38 + 2 (rz, x₀ = 0.581) + 5 = 45 | J = 23: 34 + 2 + 10 + 2 (x₀ = 0.513) + 5 = 53 | **98** (§5) |

* **n1 only.** No n3 jobs: the n3 C3H evaluator is broken (EDGE_R finding 4), and the small-ν window, where SMALLNU used n3 at
  ν = 10–17.5, matters ×16 less at 2000 (`E_low ∝ Λ⁻⁴`, §3.1): below ν = 20 the 10⁴ c₂ convention costs `≤ 0.006 N_s` (J = 21) to
  `0.01` (J = 27) — `≈ 18 %` of `E_low`, table — so n3 nodes are not needed. n1 seeds are irrelevant (fixed Sobol seeds; `4242` kept for the
  format). Settings: `ext/point.py` defaults (`HC_ADAPT = 1`, `HC_M = HC_MREF = 11`, `L1_DEG = 6`), as batches 1–3.
* **Pilot here** (`gn_pilot.py`, `cache/pilot/`, two J = 29 jobs on the 2029 row, exactly as `run_ext.py` runs them):
  * ν = 20 (`N4_J29_lo`): 214 s; replica errors 0.31 % (EE) / 0.60 % (C3H) on the state, 0.31 / 0.60 % on the reference, the same as the
    1029 node of SMALLNU (0.33 / 0.63 %). n1 works at ν = 20 on Δ ≈ 2000 (the crystal fold is felt at the same ν). Crude CRN error of the
    state/reference ratio in c₂: `680` (EE) / `1260` (C3H), against `780 / 1670` for the 1029 node by the same estimator (the chain's
    jackknife at 1029 was `≈ 730`): the per-node σ in c₂ does not grow with m², because the relative ratio error falls ×4.5.
  * x = 0.9306 (highest kept GK node, `N4_J29_edge`): 192 s; replica errors 0.24 % (state), 0.5 % (reference), `nmask 0`; crude σ ≈ 1.5·10⁴ c₂
    (the edge values at 1029 were 1–3·10⁴ c₂ with ≈ 10⁴ errors, EDGE_R): usable, as at 1000.
* **Cost.** 3.2–3.6 min per node at J = 29 (one core, this machine); J = 12, 16 nodes are ≈ 1.5–2× cheaper at 1000 (HIGHSPIN: 38 / 55
  min per 38-node row). Estimate 2.5 min (J = 12, 16) and 3.5 min (J = 21–29) per node: **main ≈ 88 × 2.5 + 159 × 3.5 ≈ 780 min ≈ 13 CPU-h
  → ≈ 50 min on 16 workers**; with group S (+98 × 3 min ≈ 5 CPU-h) **≈ 18 CPU-h → ≈ 70 min**. At the brief's pessimistic 6 min per node:
  25 / 35 CPU-h → 1.6 / 2.2 h on 16 workers. Memory as batches 1–3 (≤ 0.7 GB per job).
* **What batch 4 buys.** For J = 12, 16: removing the Laplace truncation turned 1000's `L = 0.3 / 0.1` into `3.6 / 1.9` (HIGHSPIN), then
  `7.4 / 4.4` with the window; at 2000, where every other term is ×4–16 smaller, the same is expected to give sign-free `L ≳ 10`. For
  J = 21, 27, 29: the direct-bulk + GK15-edge protocol gave `L = 1.20 / 1.64 / 1.24` at 1000; its MC, quadrature and representation
  terms are expected to shrink ×2–4 at 2000 (the pilot's per-node σ in c₂ is unchanged while the weights per c₂ fall ×4), i.e. `L ≈ 3–5`.
* **Folding needs** (not in this note): (i) the odd shadows `(2000 + J, J)` exist in `gapnear/cache/rows/gn_p2shadow.pkl` (computed with the
  rows); the even ones in `gn_shadow_k2.pkl`; J = 20, 23 have none yet (the first `hc_nodes.to_record` will compute them, as EXT_FOLD did
  for (1007, 7)); (ii) `hs_odd.py` / `sn_odd.py` hard-code `Dg = 1000`, `D = 1000 + J` and `ref_zero_x(J, 1000 + J)`; (iii) the
  normalisation fit of a row moves with injected nodes (SMALLNU §4.2) — quote windows in the direct-row normalisation.

---

## 5. The seam strip between the 1000 and 2000 twist edges

At Δ_g = 2000 the twist edge of spin J is `m/Λ = 1 + J/1999 + …`; chart A has a 1000 counterpart only from `1 + J/999`. The strip
`[tw₂₀₀₀(J), tw₁₀₀₀(J))` has width `≈ J/2000` (0.005 at J = 10, 0.016 at J = 32).

**Closed here** (a state on its own 2000 row at the 2000 twist edge, plus the chart-A image at the 1000 twist edge; between them the own
row is held, H_cont as everywhere): J = 4, 5, 6, 7, 8 [NE±] sign-free; J = 12, 16 [NE±] λ-level; J = 21, 27 [NE + M] (K₂ alone; batch 4
makes them direct); J = 29 PENDING29SEAM. The values at the 2000 twist edges are in the table (`+0.2021 … +1.0121`).

**Open** (no 2000 row; at 2000 these states sit on the 1600 row held, an H_M extrapolation):

| J | strip in m/Λ at 2000 | width | J²/m | 1000 grade (twist edge) | needed | cost |
|---|---|---|---|---|---|---|
| 9, 10, 11 | [1.0045, 1.0090), [1.0050, 1.0100), [1.0055, 1.0110) | 0.0045–0.0055 | 0.04–0.06 | [NE + M] | own K₂ row at 2000 + J (λ-level closes, as J = 8, 12) | 15–30 min each, local |
| 13, 14, 15 | [1.0065, 1.0129), [1.0070, 1.0139), [1.0075, 1.0149) | 0.006–0.007 | 0.08–0.11 | [NE + M] (J = 15 own row) | own K₂ row (λ-level, as J = 12, 16) | 30–50 min each, local |
| 17, 18, 19, 20 | [1.0085, 1.0169) … [1.0100, 1.0198) | 0.008–0.010 | 0.14–0.20 | [NE + M] (J = 20: no budget) | own K₂ row; direct row for sign-free (J = 20: group S) | ≈ 45–60 min K₂; ≈ 2–3 CPU-h direct |
| 22, 23, 24, 25, 26 | [1.0109, 1.0218) … [1.0129, 1.0257) | 0.011–0.013 | 0.24–0.33 | [NE + M] (J = 23 twist edge; J = 26: no row) | own K₂ row + direct bulk + GK edge (J = 23: group S) | ≈ 1 h K₂ + ≈ 3 CPU-h direct each |
| 28, 30, 31, 32 | [1.0139, 1.0276) … [1.0159, 1.0315) | 0.014–0.016 | 0.39–0.50 | [NE + M] / [M] | as J = 27, 29 | ≈ 1.2 h K₂ + ≈ 3 CPU-h direct each |
| 33 … J_B ≈ 52 | — | — | ≥ 0.53 | H_LJ(1000) [D-sketch + M] | chart B (GAP_SCALING §5) | — |

* Which seam states need *direct* evidence: on the evidence of §2, every J ≥ 17 (K₂ alone does not close a ×3 budget there, at 1000 or
  2000); J ≤ 16 needs only K₂ rows. The two weakest of the 1000 near-gap set (J = 20: value `+0.83 ± 0.14`, no budget; J = 23: twist edge
  [NE + M]) are batch 4's group S.
* **Plan entries** (GAP_SCALING §6.1 item 2): local K₂ rows at 2000 + J for J = 9–11, 13–15, 17–20, 22–26, 28, 30–32 (19 rows,
  ≈ 15–20 CPU-h, two lanes ≈ 9 h), each with Δ_ref = 5·10⁴ once `DREF` is a parameter (N1); then an ext batch with direct rows for
  J = 17–19, 22, 24–26, 28, 30–32 (≈ 11 × 50 nodes ≈ 30 CPU-h ≈ 2 h on 16 workers).
* At 5000 the strip is `[1 + J/4999, 1 + J/999)`, ×1.25 wider; the 5000 twist edges have no rows except J = 5 and 12 (§3.3).

---

## 6. What remains (for GAP_SCALING §6.1)
1. Run batch 4; fold it with the HIGHSPIN/SMALLNU/EXT_FOLD scripts parametrised in Δ_g (Dg, D = Dg + J, the reference-zero x₀ per row).
2. Rows at 2000 for the open seam spins (§5); `DREF` → parameter, Δ_ref = 5·10⁴ for every new 2000 / 5000 row.
3. 5000: own rows for J = 4, 6–8 (cheap; chart A is expected to hold as at 2000), J ≥ 8 with Δ_ref ≥ 10⁵ (x₀ = 0.87 at J = 5 already).
4. Parity-odd OO J = 4–10 (GAP_SCALING item 2's last clause) were not run.

## 7. Reproduce
```bash
cd ads4/code/gapnear                      # every command under: ulimit -v 3000000; OMP_NUM_THREADS=1; PYTHONDONTWRITEBYTECODE=1
bash gn_lane.sh cache/lane_A.jobs cache/lane_A.log & bash gn_lane.sh cache/lane_B.jobs cache/lane_B.log &   # rows (≈ 7.5 CPU-h) + pilot
#   e.g. python3 gn_rows.py even EE,C3H 12 2012 dense; python3 gn_rows.py odd EE,C3H 21 2021 dense; ... 0/2 at 1980; 5 5005; 12 5012
python3 gn_bind.py base1000 1000:0:990:0.989985 1000:2:990:0.989982 1000:4:1004:tw,1.01,1.02 1000:5:1005:tw,1.01,1.02,1.05 \
   1000:6:1006:tw,1.01,1.02 1000:7:1007:tw,1.01,1.02,1.05 1000:8:1008:tw,1.02,1.05 1000:12:1012:tw,1.02,1.05 \
   1000:16:1016:tw,1.02,1.03,1.05 1000:21:1021:tw,1.03,1.05,1.075,1.1 1000:27:1027:tw,1.03,1.04,1.05,1.075,1.1 \
   1000:29:1029:tw,1.04,1.05,1.06,1.075,1.1
python3 gn_bind.py light 2000:0:1980:0.989985 2000:2:1980:0.989982;  python3 gn_bind.py J4 2000:4:2004:tw,tw1000,1.01,1.02
python3 gn_bind.py J5 2000:5:2005:tw,tw1000,1.01,1.02,1.05;  python3 gn_bind.py J7 2000:7:2007:tw,tw1000,1.01,1.02,1.05
python3 gn_bind.py J6J8 2000:6:2006:tw,tw1000,1.01,1.02 2000:8:2008:tw,tw1000,1.02,1.05
python3 gn_bind.py J12J16 2000:12:2012:tw,tw1000,1.02,1.05 2000:16:2016:tw,tw1000,1.02,1.03,1.05
python3 gn_bind.py J21 2000:21:2021:tw,tw1000,1.03,1.05,1.075,1.1;  python3 gn_bind.py J27 2000:27:2027:tw,tw1000,1.03,1.04,1.05,1.075,1.1
python3 gn_bind.py J29 2000:29:2029:tw,tw1000,1.04,1.05,1.06,1.075,1.1
python3 gn_bind.py J5k 5000:5:5005:tw,2000tw,tw1000,1.01,1.02,1.05;  python3 gn_bind.py J12k 5000:12:5012:tw,tw1000,1.02,1.05
python3 gn_bind.py --hmono 0,2 (likewise 4; 5,6,7,8,12,16; 21; 27; 29);  python3 gn_hsum.py > cache/hsum.txt
python3 gn_chartA.py light,J4,J5,J6J8,J7,J12J16,J21,J27,J29,J5k,J12k > cache/chartA.txt;  python3 gn_table.py <same tags> > cache/table.txt
python3 gn_refdiag.py 21 1021,2021; python3 gn_stepnorm.py 21 1021,2021 (and 16, 27)
python3 gn_pilot.py n1 29 2029 20 4242 N4_J29_lo; python3 gn_pilot.py n1 29 2029 1886.599298 4242 N4_J29_edge
python3 gn_batch4.py                      # -> ../../ext/jobs_batch4.txt
```

## 8. Caveats
1. Everything at 2000 / 5000 here is K₂; "K₂ = crystal integral" is checked at 1000 (direct rows), not at 2000. Batch 4 is that check.
2. The ×3 factors and the 10⁴ c₂ conventions are NEARGAP2's choices, unchanged; the gauge factor at 2000 / 5000 comes from the
   min-norm *grid* weights (GAP_SCALING §1.2), not the published continuum weights.
3. Only the listed masses were evaluated; between them, H_cont as before.
4. The reference-zero rule is HIGHSPIN's; at 2000 it flagged nothing, and the contamination of §2.4 is below its threshold but visible in
   the K₁ rows. It inflates, never deflates, the truncation estimate as used here (the K₂ values agree with 1000).
5. Graded values at 1000 for J ≥ 7 are quoted from NEARGAP2 / HIGHSPIN / SMALLNU_WINDOW / EXT_FOLD, not recomputed.
