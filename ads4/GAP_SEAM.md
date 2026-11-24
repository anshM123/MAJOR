# GAP_SEAM: the seam / near-gap spins J = 9–11, 13–15, 17, 19 and the parity-odd OO block J = 4–10 at Δ_g = 2000, with a parameterised Δ_ref

*24 Nov 2026, gap-seam agent. Unreviewed, not committed by this agent. Another agent's commit (f2c01f6, "Gap Morse") swept in an intermediate
snapshot of `code/gapseam/`; the working tree holds the final versions. No existing file was edited. `GAP_MORSE.md` / `code/gapmorse/`
were not touched.*
* *New code:* `ads4/code/gapseam/`. `cache/.gitignore` ignores everything but itself.
  * `gsm_guard.py`: write guard; every write outside `gapseam/cache/` raises. It registers itself as `ng_guard`, `gn_guard`, `hs_guard`, ….
  * `gsm_rows.py`: own-mass K₂ rows at Δ = 2000 + J (copy of `gapnear/gn_rows.py`).
    * **`Δ_ref` is a parameter** (`GSM_DREF`, default 5·10⁴).
    * It reuses gapnear's state records (read-only) when they exist, so the J = 16, 21 reruns only compute the new reference points.
    * `refb` mode: reference-basis points for even J. `oo` mode: the 3×3 parity-odd action at even J.
  * `gsm_bind.py`: NEARGAP2's budget (copy of `gapnear/gn_bind.py`).
    * `p4_odd_analysis.grid` and `l9_analysis.grid` are re-typed with **a per-row reference `Δ_ref(J, D)`**. Every other step is the
      module's own function, called through the module.
    * Also new: the OO (++-diagonal) 1×1 value and budget, and `--hmono` with a per-row reference.
  * Other scripts:
    * `gsm_ooid.py`: the Γ_OO = Γ_++ check at 2000.
    * `gsm_hscan.py`: hybrid minima over m/Λ.
    * `gsm_table.py`: the tables.
    * `gsm_lane.sh`, `gsm_pause.sh`: compute lanes, and short analyses run with one lane SIGSTOPped.
* *Read-only imports:* `neargap2/` (`ng_chain`, `ng_bind`) → `spinfunc/` (`ee1000_common`, `w2_chain`, `q3_analysis`, `p4_odd_analysis`,
  `l9_analysis`, `l8_analysis`), and the row engines (`w1_rows`, `p3_odd_action`, `p2_shadow_odd`, `l9_points`, `l5_shadow_matrix`,
  `w0_ee_common`).
  * Read-only data: the rows of `spinfunc/cache`, `neargap2/cache`, `gapscale/cache/rows` and `gapnear/cache/rows`, and
    `gapscale/cache/weights_2000.npz`.
* *Compute:* ≈ 7.5 CPU-h, on two lanes (08:05–12:00 and 08:05–11:27), `ulimit -v 3000000`, one thread, at most 2 processes.
  * ≈ 1 360 new K₂ points; 13–26 s per point at Δ = 2000 + J and at Δ_ref = 5·10⁴ alike.

Grades are as in GAP_NEARGAP:
* **[NE±]**: the sign-free budget closes (L > 1).
* **[NE±] λ-level**: only the λ-level budget closes.
* **[NE + M]**: neither budget closes.

`L = value / budget`; "(λ)" is the λ-level leverage. The budget is NEARGAP2's, unchanged:
* ×3 on the K₂ truncation and on the interpolation;
* `|Γ| ≤ 10⁴ c₂` sign-free on ν < 50, with the gap's min-norm gauge factor;
* `G_edge = max(10⁴, 10 × max|Γ|)` beyond the row's last node.

Λ = 998.5 / 1998.5 at Δ_g = 1000 / 2000. "tw" is the 2000 twist edge; "tw1000" is the 1000 twist edge. Every new 2000 row uses
**Δ_ref = 5·10⁴**.

---

> **Update (GAP_SEAM2.md, unreviewed): the Δ_g = 2000 reruns of J = 27, 29 with Δ_ref = 5·10⁴ did not lift their twist edges. They stay [NE+M]; no state is negative (all ≥ +0.91).**
> * **Why the expected doubling failed:** the λ-level truncation step fell as predicted (chart-A ratio 0.22–0.25), but two other terms took over. The edge term (the row stops at 0.8 while |Γ| ≈ 10⁴ c₂ at x = 0.85–0.9) and the small-ν window (larger when the first node is 0.05) dominate.
> * **Chart A holds cleanly** once the node set and the reference are matched: E_trunc ratio 0.26–0.27 and λ-step ratio 0.22–0.25, against 0.25 predicted. Γ(x) at 2027 matches 1027 to ≤ 40 c₂ for x ≤ 0.6.
> * **New lower state:** J = 27 at 1.0175Λ, +0.910, inside the seam strip. J = 29 at 1.075/1.10Λ and J = 27 at 1.075/1.10Λ are λ-level [NE±].
> * **J = 27, 29 twist edges at 2000 need the batch-4 direct rows (or a K₃ order).** J = 18, 20, 22, 24 not done; resume commands are in GAP_SEAM2 §3. Hybrid minima for J = 22/24 at 2000: +0.780 / +0.826.

## 0. Table and verdict

**1000 column.** This is the same protocol (`ng_bind`) at Δ_g = 1000, on the rows that exist there: own rows `1000 + J` only for J = 15,
16, 17, 19, 21; for J = 9–14 the chain interpolates the 400 / 800 / 1600 rows. It reproduces `gapnear/cache/gn_bind_base1000.json` to
every digit where both exist. The graded 1000 values (RESULT_D_STATEMENT §5.2, minima over masses, all [NE + M]) are J = 9 / 10 / 11 /
13 / 14 / 15 / 17 / 19 = `+0.4542 / +0.5886 / +0.3482 / +0.4076 / +0.6982 / +0.869 (own row tw +0.6900) / +0.840 / +0.916`.

### 0.1 EE block (odd J: parity-odd `+−` 1×1; even J: 2×2 minimum)

| state J (m/Λ: tw, tw1000, others) | 1000: value / L (λ) at tw1000; other masses | **2000 own row (Δ_ref 5·10⁴): value / L (λ)** at tw; tw1000; others | grade 2000 |
|---|---|---|---|
| 9 (1.00449; 1.00897; 1.02, **1.03** (hybrid min), 1.05) | +0.6763 / 1.5 (3.6); +0.4862 / 1.7 (4.9), **+0.4544 / 2.0 (6.3)**, +0.5293 / 2.3 (8.7) | +0.9275 / **8.0** (47); +0.6763 / 8.2 (52); +0.4862 / 8.8 (60), **+0.4544 / 9.3 (64)**, +0.5294 / 10.2 (67) | **[NE±]** |
| 10 (1.00499; 1.00996; **1.0075**, 1.02, 1.05) | +0.5857 / 0.4 (1.2); —, +0.8030 / 0.6 (1.9), +0.8960 / 0.9 (3.7) | +0.5820 / **1.8** (14); +0.5879 / 2.0 (17); **+0.5759 / 1.9 (16)**, +0.8049 / 2.7 (23), +0.8970 / 3.8 (29) | **[NE±]** |
| 11 (1.00549; 1.01095; 1.02, 1.05) | +0.3488 / 0.5 (1.7); +0.4099 / 0.5 (2.0), +0.9984 / 0.6 (3.3) | **+0.3349 / 2.6** (19); +0.3490 / 2.6 (21); +0.4102 / 2.6 (23), +0.9998 / 2.8 (25) | **[NE±]** |
| 13 (1.00648; 1.01293; 1.02, 1.05) | +0.4077 / 0.2 (1.0); +0.6028 / 0.2 (1.1), +0.9975 / 0.2 (2.1) | **+0.3027 / 1.03** (7.4); +0.4084 / 1.03 (9.0); +0.6040 / 1.04 (10.2), +0.9997 / 1.3 (14) | **[NE±] (sign-free thin, L = 1.03)**; λ-level 7–14 |
| 14 (1.00698; 1.01392; 1.02, 1.035, 1.05) | +0.9962 / 0.1 (0.7); +0.9963 / 0.2 (0.8), +0.6976 / 0.2 (1.1), +0.7920 / 0.2 (1.4) | +0.6842 / 0.5 (**5.2**); +0.9991 / 0.6 (7.0); +0.9994 / 0.7 (8.5), +0.6992 / 0.7 (10.1), +0.7935 / 0.8 (10.9) | **[NE±] λ-level** |
| 15 (1.00748; 1.01490; 1.02, 1.05) | +0.6932 / 0.1 (1.3); +0.9972 / 0.1 (1.4), +0.9451 / 0.2 (2.0) | **+0.3659** / 0.6 (**5.0**); +0.6946 / 0.6 (5.9); +0.9999 / 0.6 (6.5), +0.9481 / 0.8 (9.5) | **[NE±] λ-level** |
| 16, Δ_ref 5·10⁴ (1.00797; 1.01589; 1.02, 1.03, 1.05) | +0.9344 / 0.1 (1.0); +0.7635 / 0.1 (1.0), +0.7466 / 0.1 (1.1), +0.9949 / 0.1 (1.3) | +0.9997 / 0.4 (**4.1**); +0.9399 / 0.4 (5.1); +0.7676 / 0.4 (5.4), +0.7501 / 0.5 (6.1), +0.9996 / 0.6 (7.1) | [NE±] λ-level (Δ_ref 2·10⁴: 3.7; 4.4; 4.6, 5.1, 6.4) |
| 17 (1.00847; 1.01687; 1.02, 1.05) | +0.9991 / 0.1 (0.4); +0.9989 / 0.1 (0.5), +0.8417 / 0.2 (0.9) | **+0.5329** / 0.4 (**3.3**); +1.0000 / 0.4 (4.1); +1.0000 / 0.4 (4.3), +0.8440 / 0.5 (6.4) | **[NE±] λ-level** (1000: not bounded) |
| 19 (1.00946; 1.01884; 1.02, 1.0375, 1.05) | +0.9918 / 0.1 (0.3); +0.9920 / 0.1 (0.3), +0.8626 / 0.2 (0.5), +0.9212 / 0.2 (0.6) | +0.8552 / 0.2 (**1.7**); +0.9989 / 0.2 (2.5); +0.9989 / 0.2 (2.5), +0.8663 / 0.3 (3.7), +0.9240 / 0.3 (4.3) | **[NE±] λ-level** (thin at tw) |
| 21, Δ_ref 5·10⁴ (1.01045; 1.02080; 1.03, 1.05, 1.075, 1.10) | +1.0009 / 0.07 (0.17) [graded direct +0.998 ± 0.020, L 1.20]; +0.8821 / (0.24), +0.9983 / (0.37), +0.9588 / (0.60), +0.9492 / (0.75) | +0.9973 / 0.1 (**1.7**); +0.9979 / 0.2 (2.3); +0.8813 / (2.7), +0.9987 / (3.2), +0.9605 / (3.6), +0.9512 / (3.9) | **[NE±] λ-level at every mass, the twist edge included** (Δ_ref 2·10⁴: tw 0.9 [NE + M]; 1.1–2.0 elsewhere) |

### 0.2 Parity-odd OO block, even J (1×1 = normalised `++` diagonal, Γ_OO = Γ_++; §3)

| J | 1000: value / L (λ) | **2000: value / L (λ)** | grade 2000 |
|---|---|---|---|
| 4 | tw1000 1.00400: +0.2059 / 22.7 (33); 1.01 (HC min): +0.2150 / 24.2 (35); 1.02: +0.2311 / 27 (39) | tw 1.00200: **+0.2029 / 106** (172); tw1000: +0.2059 / 108 (175); 1.01: +0.2151 / 116 (186); 1.02: +0.2312 / 129 (206) | **[NE±]** |
| 6 | tw1000 1.00599: +0.3551 / 8.0 (17); 1.02: +0.3567 / 8.8 (18); 1.25 (HC min): +0.3093 / 20.7 (32) | tw 1.00300: +0.3548 / **34.8** (86); tw1000: +0.3552 / 35.5 (87); 1.02: +0.3567 / 38.6 (93); 1.25: **+0.3093** / 87.6 (145) (held row) | **[NE±]** |
| 8 | tw1000 1.00798: +0.5884 / 2.7 (8.0); 1.02: +0.5783 / 2.9 (8.4); 1.05: +0.5606 / 3.4 (9.3); 1.75 (HC min): +0.5033 / 43 (57) | tw 1.00399: +0.5924 / **11.4** (40); tw1000: +0.5887 / 11.6 (40); 1.02: +0.5786 / 12.3 (41); 1.05: +0.5608 / 14.2 (46); 1.75: **+0.5033** / 167 (232) (held row) | **[NE±]** |
| 10 | tw1000 1.00996: +0.9378 / 1.4 (5.6); 1.02: +0.9276 / 1.5 (5.8); 1.05: +0.8960 / 1.7 (6.2); 2.0 (HC min): +0.6046 / 31 (45) | tw 1.00499: +0.9441 / **5.7** (28); 1.0075: +0.9415 / 5.8 (28); tw1000: +0.9391 / 5.9 (28); 1.02: +0.9288 / 6.1 (29); 1.05: +0.8970 / 6.9 (31); 2.0: **+0.6046** / 143 (217) (held row) | **[NE±]** |
| (14, 16: by-products) | 14: +0.996 / 0.4–0.5 (2.2–2.5); 16: +0.994 / 0.3 (1.5–1.7) | 14: +0.9997 / 1.5–1.9 (11–12); 16: +0.9994–0.9997 / 1.00–1.24 (7.0–8.1) | 14: [NE±]; 16: [NE±] (razor-thin at tw: L = 1.002 / 1.005); identity not checked at 2000 for J ≥ 12 |

### 0.3 Verdict
1. **No state is negative.** The new lowest values are hybrid dips inside the seam strip at the 2000 twist edges:
   * J = 13: `+0.3027` (L = 1.03, λ 7.4);
   * J = 11: `+0.3349` (L 2.6);
   * J = 15: `+0.3659` (λ 5.0).

   All are ×1.5–1.8 above the binding state (`+0.2004`). They are dips of the hybrid itself (`+0.3028 / +0.3351 / +0.3658`; `gsm_hscan.py`),
   not of the corrections. **The seam strip is lower than its 1000 counterpart for J = 11, 13, 15, 17** (e.g. J = 15: `+0.366` at the
   2000 edge vs `+0.693` at the 1000 edge), so the strip is not a chart-A image and needs its own rows, as GAP_SCALING §5(iv) anticipated.
2. **J = 9, 10, 11, 13 are sign-free [NE±] at 2000 from K₂ alone** (L = 8–10 / 1.8–3.8 / 2.6–2.8 / 1.03–1.3). At 1000 they were [NE + M]
   (same protocol, no own row: L = 0.2–2.3).
3. **J = 14, 15, 16, 17, 19, 21 are [NE±] λ-level at 2000** (λ-level L ≥ 5.2 / 5.0 / 4.1 / 3.3 / 1.7 / 1.7). Sign-free they do not close
   (L = 0.1–0.8). At 1000, K₂ gave λ-level L = 0.2–2.0 for these spins (J = 15, 16: 1.0–2.0, graded [NE + M] there). **J = 21's twist edge moves from [NE + M] to [NE±] λ-level**
   through the new reference alone (§2).
4. **OO block J = 4, 6, 8, 10: sign-free [NE±] at 2000** (L ≥ 106 / 35 / 11 / 5.7), with L ×3.9–4.8 larger than at 1000.
   * The 1000 hc_oo minima (`+0.2150 / +0.3093 / +0.5032 / +0.6046`) are reproduced to 10⁻⁴, and are now error-budgeted at 1000 too
     (L = 24 / 21 / 43 / 31).
   * The identity Γ_OO = Γ_++ holds on the 2000 rows to ≤ 0.044 c₂ (§3). The budget carries it as a 10 c₂ sign-free term.
5. **Δ_ref artefact (GAP_NEARGAP N1): confirmed for the λ-level truncation step, which is now chart-A-clean. It was not the cause of
   the sign-free E_trunc inflation** (§2).
   * At J = 21 the λ-level K₂ step falls ×2.5–2.7 with Δ_ref = 5·10⁴ (tw1000: `0.781 → 0.297`). Its 2000/1000 ratio becomes
     **0.24** (2·10⁴: 0.48–0.64). The predicted ratio is 0.25.
   * At J = 16 the zero sat at the row's end (x₀ = 0.72): the λ-level step falls ×1.1–1.2. The despiked node x = 0.725 comes back, and Γ at
     x = 0.70–0.75 now agrees with the 1016 row to ≤ 60 c₂ (was up to 330 c₂).
   * **The sign-free E_trunc at J ≥ 17 still scales by only 0.46–0.70.** This does not depend on Δ_ref. Its cause is the node set: the
     1000 rows 1017 / 1019 / 1021 start at x = 0.05, the 2000 rows at x = 0.025. That node is inside the per-row a-fit (x ≤ 0.25, ν ≥ 15)
     and renormalises the K = 1 row. Dropping it on 2021 (`GSM_XMIN=0.045`, diagnostic only) halves E_trunc, and the ratio becomes
     ≈ 0.29. Where both gaps have the same node set (J = 15, 16), every term scales by 0.25–0.26.
6. **Chart A holds** on every budget term whose definition is the same at both gaps (§4):
   * the λ-level step: 0.23–0.29 (with Δ_ref = 5·10⁴, for every J computed);
   * E_low: 0.059–0.065 ≈ Λ⁻⁴;
   * E_trunc: 0.24–0.28 (J = 9–16; J = 13: 0.18);
   * Γ(x) on the 2000 rows agrees with the 1600 / 1015 / 1017 rows to ≤ 4–19 c₂ (J = 9, 10) and ≤ 80–140 c₂ (J = 15, 17) over x ≤ 0.8,
     except isolated nodes where a coarse 1000 row is interpolated.

**Negative / cautionary findings (reported prominently).**
* **N1. Sign-free budgets do not close from J = 14 on, even at 2000 with Δ_ref = 5·10⁴** (L = 0.1–0.8). The K₂ step per node is a property
  of the state's own Laplace series (≈ 1–4·10³ c₂ at every gap). It falls only as Λ⁻², which is not enough at J ≥ 14 with a ×3 safety
  factor. J = 13 closes sign-free by only 3 % (lower bound `+0.009`).
* **N2. The seam strip contains the lowest near-gap values of J = 11, 13, 15, 17 at Δ_g = 2000** (`+0.30 … +0.53`, against
  `+0.35 … +1.00` at the 1000 twist edges, §0.3 item 1). They are positive, and the K₂ corrections there are small
  (|λ_K2 − λ_hyb| ≤ 6·10⁻⁴), but chart A cannot supply these states.
* **N3. The OO values at the HC minimum masses for J = 6, 8, 10** (1.25 / 1.75 / 2.0Λ) lie beyond the own row. The 2000 chain *holds*
  the 2000 + J row there (at 1000 it interpolated 1000 + J ↔ 1600). The corrections are ≤ 3·10⁻⁴ of the value there.
* **N4. At 1000 the OO twist edge of J = 4 (`1.0040Λ`: `+0.2059`) is below HYPOTHESES_CLOSE's quoted minimum `+0.2150` [1.01].** That
  grid did not include the twist edge. It is positive and bounded (L = 22.7).
* **N5. Small inconsistency in GAP_NEARGAP's table.** Its "1000, same protocol" entry for J = 21 at the 1000 twist edge reads `+0.967 / 0.07 (0.17)`.
  Its own `gn_bind_base1000.json` (reproduced here) has `+1.0009 / 0.07 (0.17)`. The leverages agree; only the printed value differs.
* **J = 18, 20 were not run** (budget). J = 22–26, 28, 30–32 remain open (§5).

---

## 1. Protocol

### 1.1 Rows
* K₂ crystal points (`L1_DEG = 6`, the w1 engines unchanged):
  * odd J: the parity-odd `p3` integrand with `E`;
  * even J: `l9_points` with the w0 EE integrand.
* Kinds are EE and C3H. Each row has 36 nodes, `x = 0.025 … 0.9`, at Δ = 2000 + J, with the reference **Δ_ref = 5·10⁴ at the same ν**.
* Reference zero `x₀ = 2.405 m_ref / ((2J + 1) m_row)`: with 5·10⁴ it lies at 1.4–3.5 for J ≤ 21, outside every row. With 2·10⁴ it was
  0.56 (J = 21) / 0.72 (J = 16).
* **Even J needs a reference helicity basis at the new Δ_ref.** `l8_analysis.basis` requires ≥ 4 points with x_ref ≥ 0.15. The `refb`
  mode adds six EE points at `(5·10⁴, x_ref = 0.3 … 0.8)`, i.e. ν = 1.5–4·10⁴ (≈ 15 s each).
* New shadows: `(2000 + J, J)` and `(5·10⁴, J)`. Shadow residuals at 5·10⁴ are ≤ 6·10⁻⁸ (EE / C3H). There were no non-finite records,
  no phase-strip failure, and no `refbad` node on any new row.
* Rows computed:
  * J = 9, 10, 11, 13, 14, 15, 17, 19: full rows, 144 points each;
  * J = 16, 21: only the 72 reference points; the state records are gapnear's;
  * refb for J = 10, 14, 16;
  * OO points (§3).

### 1.2 Per-row reference in the chain (`gsm_bind.py`)
`p4_grid` / `l9_grid` are line-by-line copies of `p4_odd_analysis.grid` / `l9_analysis.grid` with three changes:
* `DREF → dref(J, D)`;
* the row filter `k[0] < DREF` becomes `< 10⁴`;
* `Bref` is built per row.

They are installed where the chain calls them: `ng_chain._P4grid` and `w2_chain._A_grid`. Through those, the w2 Hampel filter, the q3 row
filter / alignment, the ee1000 Jacobi fix and de-shift, and ng_chain's Laplace-order control all apply unchanged. The `refbad` rule and
`ref_zero_x` use the row's own reference.

**Reproduction check** (`cache/gsm_bind_chk.txt`): with all references at 2·10⁴ the copy reproduces gn_bind to every printed digit:
* J = 5 at 2000: `+0.20041`, L = 319.8 / 564.6;
* J = 12: `+0.40891`, 0.7 / 5.4;
* J = 16: `+0.99937`, 0.4 / 3.7;
* and all 1000 entries of `gn_bind_base1000.json`.

### 1.3 States
* The 2000 twist edge (tw), the 1000 twist edge (tw1000), and 1.02 / 1.05Λ.
* The hybrid minima over m/Λ ∈ [tw, 1.15] from `gsm_hscan.py`:
  * J = 9 at 1.03 (`+0.4545`), where RESULT_D's graded +0.4542 sits;
  * J = 10 at 1.0075;
  * J = 14 at 1.035;
  * J = 19 at 1.0375.
* For J = 11, 13, 15, 17 the minimum *is* the twist edge at both gaps.
* At 1000 the same set minus the masses below its twist edge.

---

## 2. The Δ_ref artefact (J = 16, 21)

| | J = 16, Δ_ref 2·10⁴ → 5·10⁴ | J = 21, Δ_ref 2·10⁴ → 5·10⁴ |
|---|---|---|
| x₀ on the 2000 row | 0.723 → 1.807 | 0.563 → 1.408 |
| λ-level step 3\|λ_K2 − λ_K1\| (tw; tw1000; 1.05) | 0.161 → 0.141; 0.148 → 0.125; 0.132 → 0.119 | **0.858 → 0.316; 0.781 → 0.297; 0.647 → 0.274** |
| ratio of that step to 1000 (tw1000 … 1.05) | 0.28–0.31 → **0.25–0.26** | 0.48–0.64 → **0.24** |
| λ-level L (tw; tw1000; others) | 3.7 → 4.1; 4.4 → 5.1; 4.6–6.4 → 5.4–7.1 | **0.9 → 1.7**; 1.1 → 2.3; 1.2–2.0 → 2.7–3.9 |
| sign-free E_trunc (tw1000) | 2.04 → 2.12 | 6.26 → 5.98 |
| E_interp (tw) | 2.5·10⁻³ → 2.7·10⁻⁴ | 2.8·10⁻² → 1.7·10⁻² (the two every-other-node variants now differ by 0.010 instead of 0.016) |
| value (tw1000) | +0.9418 → +0.9399 (1000 direct row: +0.9402) | +1.0051 → +0.9979 (1000 direct: +0.998 ± 0.020) |
| Hampel-despiked nodes on 2016 | 0.725, 0.85, 0.875 → 0.85, 0.875 | none → none |

Γ near the zero (`cache/gsm_hmono_J16.txt`), C2++ at x = 0.70 / 0.75:
* 1016 row: `+280 / −220`;
* 2016, 2·10⁴: `+101 / −271`;
* 2016, 5·10⁴: `+221 / −216`.

The K₂−K₁ steps there fall by 10–60 % (C3H++ x = 0.70: 1684 → 652). The x = 0.85 sign difference between 1016 and 2016 is unchanged by
the reference: it is the unresolved K₂ edge (GAP_NEARGAP §3.2).

**Why the sign-free E_trunc does not follow** (`cache/gsm_hmono_J21x.txt`, `gsm_bind_J21x.txt`, with `GSM_XMIN=0.045`):
* The per-node |Γ_K2 − Γ_K1| on 2021 is 1–2·10³ c₂ at x = 0.1–0.7 *for both references*. On 1021 it is 80–900.
* Removing the x = 0.025 node of 2021, which 1021 does not have, brings the 2021 steps back to the 1021 size (x = 0.3: 1978 → 347;
  x = 0.7: 1455 → 152).
* E_trunc then falls 5.98 → 3.02 (tw1000), and the λ-level step is unchanged (0.297 → 0.312).

So:
* **the reference zero inflated the λ-level step**, and Δ_ref ≥ 5·10⁴ removes that inflation;
* **the sign-free inflation at J ≥ 17 is the x = 0.025 node entering the K = 1 a-fit**: a node-set difference between the 1000 and 2000
  rows, not a reference or physics effect.

The grades here use the full protocol (all nodes), as GAP_NEARGAP does. The diagnostic changes no grade, since sign-free L stays < 1 at
J = 21 either way (0.3).

---

## 3. The parity-odd OO block, J = 4–10 at 2000

**Protocol.** The protocol at 1000 is HYPOTHESES_CLOSE §3: K₂ OO rows of the 3×3 action, folded into the chain;
`λ_OO/N_s = (base_++ + dF_OO)/N_s,++`; graded [NE + M] with no budget. PARITY_ODD_CLOSE §1 and HYPOTHESES_CLOSE2 §1.4 showed
Γ_OO = Γ_++ node by node (≤ 0.05 / 0.12 c₂, falling as Δ⁻⁴). So here:
* `λ_OO/N_s` is the normalised `++` diagonal of the EE chain on the own rows;
* its 1×1 budget is the `++` entry of every NEARGAP2 term, plus `E_id`, a sign-free term for |Γ_OO − Γ_++| ≤ 10 c₂ on every node;
* J = 4, 6, 8 use gapnear's 2004 / 2006 / 2008 rows (Δ_ref 2·10⁴, x₀ = 1.4–2.7, outside); J = 10 uses the new 2010 row.

At 1000 the same construction reproduces hc_oo's minima, `+0.2150 / +0.3093 / +0.5032 (here +0.5033) / +0.6046`, to 10⁻⁴.

**Identity check at 2000** (`gsm_ooid.py` → `cache/gsm_ooid.txt`):
* K₂ points of the 3×3 action on 2004–2010 at x = 0.1, 0.2, 0.3, 0.5, 0.7, 0.85, with kinds EE and C3H. That is 47 points, state side
  only. On the reference side κ is constant to 10⁻¹³ (PARITY_ODD_CLOSE).
* T0 is the row's own order-0 helicity basis, and `κ_k = M_k[OO] / (T0⁻¹ M_k T0⁻ᵀ)_++`.
* The order-0 κ is constant in x to ≤ 7·10⁻¹⁴ (EE) / ≤ 3·10⁻⁸ (C3H).
* The higher orders κ_k/κ₀ − 1 are ≤ 3·10⁻⁸ (EE) / ≤ 7·10⁻⁶ (C3H).
* The c₂-unit deviation `m²|ρ/ρ_ref − 1|` (ρ the K₂-summed ratio) is ≤ 10⁻⁶ (EE) and ≤ 0.003 / 0.004 / 0.004 / **0.044 c₂** (C3H; J = 4 /
  6 / 8 / 10).
* `G_ID = 10 c₂` is therefore ×230 conservative. It costs 1·10⁻⁵ … 7·10⁻⁴ N_s.

**Result** (§0.2): sign-free [NE±] at every state, with L ×3.9–4.8 the 1000 values.
* The OO twist edges at 2000 are `+0.2029` (J = 4, the lowest OO value), `+0.3548`, `+0.5924`, `+0.9441`.
* The chart-A ratios of every term are 0.252–0.261 (E_low 0.065).
* At 1.25 / 1.75 / 2.0Λ (held own row) the corrections are below 3·10⁻⁴ and L = 88–167.

---

## 4. Chart A (`cache/table.txt`; ratios 2000 / 1000 at fixed (J, m/Λ); expected 0.2496 for c₂-type terms, ≈ 0.062 for E_low)

| J | E_trunc | λ-level step | E_low | L growth sign-free (λ) | note |
|---|---|---|---|---|---|
| 9, 10, 11 | 0.24–0.28 | 0.25–0.28 | 0.059–0.065 | ×4.2–5.5 (×7.5–14) | 1000 rows 200–1600 interpolated (not own rows) |
| 13 | 0.18 | 0.28–0.29 | 0.065 | ×5–6 (×7–9) | 1000: rows 400 / 800 only |
| 14 | 0.26 | 0.26–0.27 | 0.065 | ×4–6 (×8–11) | |
| 15 (like-for-like 1015) | 0.254–0.255 | 0.252 | 0.065 | ×4–6 (×4.5–4.8) | same node set at both gaps |
| 16 (like-for-like 1016), Δ_ref 5·10⁴ | 0.257–0.258 | 0.252–0.261 | 0.065 | ×4 (×5.1–5.5) | 2·10⁴: λ-step 0.28–0.31 |
| 17 / 19 (1017 / 1019) | 0.46 / 0.58–0.61 | **0.25–0.27 / 0.23–0.24** | 0.065 | ×2–3 (×7–10) | E_trunc: x = 0.025 node (§2) |
| 21 (1021), Δ_ref 5·10⁴ | 0.58–0.68 | **0.24** | 0.065 | ×2–3 (×5–13) | 2·10⁴: λ-step 0.48–0.64 |
| OO 4–10 | 0.252–0.261 | 0.253–0.258 | 0.065 | ×3.9–4.8 | |

E_edge is not comparable for J = 9–14, 17, 19, 21: the 1000 rows end at x = 0.8, the 2000 rows at 0.9. Where it is comparable (J = 15, 16,
OO 4–8) its ratio is 0.25–0.28.

H_mono-A (Γ in c₂ at fixed (J, x), K₂; `cache/gsm_hmono_*.txt`), max |Γ₂₀₀₀ − Γ₁₀₀₀-side| over x ≤ 0.8:

| J | rows compared | max difference | max |Γ| | where |
|---|---|---|---|---|
| 9 | 1600 vs 2009 | 12 c₂ | 620 | |
| 10 | 1600 vs 2010 | 4–19 c₂ | 1053 | |
| 14 | 1600 vs 2014 | ≤ 6 c₂ for x ≤ 0.7 | | 710–816 c₂ only at x = 0.75, a node despiked on the 16-node 1600 row and interpolated |
| 11 / 13 | 800 vs 2011 / 2013 | 350–620 c₂ | 1.3–3.5·10³ | different J²/m (0.15 vs 0.06), and the 11-node 800 row at x = 0.75 |
| 15 / 17 | 1015 / 1017 vs 2015 / 2017 | 81–142 c₂ | 2–3·10³ | |
| 21 | 1021 vs 2021 (5·10⁴) | ≤ 45 c₂ on 0.1 ≤ x ≤ 0.8 | 1.5–1.8·10³ | 163 c₂ at x = 0.05, where the K₂ step is 6–7·10³ on both rows |

---

## 5. What remains
1. **Seam spins not run:** J = 18, 20 (even; they need refb + rows, ≈ 1 CPU-h each), and J = 22–26, 28, 30–32 (≈ 1–1.3 CPU-h each with
   Δ_ref 5·10⁴).
   * On the evidence here, K₂ alone gives λ-level [NE±] up to at least J = 21 at 2000 once Δ_ref = 5·10⁴.
   * The J = 27, 29 twist edges (GAP_NEARGAP: [NE + M], λ-level L = 0.4–0.6 with 2·10⁴) should be rerun with 5·10⁴ references:
     72 points each, ≈ 45 min. The J = 21 result suggests they would roughly double their λ-level L.
2. **Sign-free closure for J ≥ 14 at 2000** needs direct rows (batch 4 covers J = 12, 16, 21, 27, 29 and seam J = 20, 23), or a K₃ order.
   * The node-set effect of §2 suggests using the same x grid at both gaps, or excluding x = 0.025 from the a-fit, in any future chart-A
     comparison.
3. **Batch 4** should use Δ_ref ≥ 5·10⁴ for any new reference points on the 2000 rows. Its direct rows divide by the reference without
   expanding it, so the 2·10⁴ choice there matters only at the node nearest x₀.
4. The seam strip's own grade rests on K₂ = crystal integral at 2000 (as in GAP_NEARGAP); there is no direct check at 2000 yet.
5. OO beyond J = 10: by-products show J = 14 sign-free (L 1.5–1.9) and J = 16 sign-free but razor-thin (L 1.00–1.24). The identity was not checked at 2000 for J ≥ 12. J ≥ 18 is not run.

## 6. Reproduce
```bash
cd ads4/code/gapseam          # every command: ulimit -v 3000000; OMP_NUM_THREADS=1; PYTHONDONTWRITEBYTECODE=1 (gsm_lane.sh / gsm_pause.sh set them)
bash gsm_lane.sh cache/lane_A.jobs cache/lane_A.log & bash gsm_lane.sh cache/lane_B.jobs cache/lane_B.log &   # rows, refb, OO points (≈ 7.3 CPU-h)
#   e.g. python3 gsm_rows.py odd EE,C3H 9 2009 dense;  python3 gsm_rows.py refb EE 10;  python3 gsm_rows.py even EE,C3H 10 2010 dense
#        GSM_NOREF=1 python3 gsm_rows.py oo EE,C3H 4 2004 0.1,0.2,0.3,0.5,0.7,0.85
python3 gsm_bind.py chk 2000:5:2005:tw 2000:12:2012:tw 2000:16:2016:tw                  # reproduction of gn_bind (all Δ_ref = 2e4)
python3 gsm_bind.py base1000 1000:9:1009:tw,1.02,1.03,1.05 1000:10:1010:tw,1.02,1.05,2.0 1000:11:1011:tw,1.02,1.05 \
  1000:13:1013:tw,1.02,1.05 1000:14:1014:tw,1.02,1.035,1.05 1000:15:1015:tw,1.02,1.05 1000:16:1016:tw,1.02,1.03,1.05 \
  1000:17:1017:tw,1.02,1.05 1000:18:1018:tw,1.02,1.05 1000:19:1019:tw,1.02,1.0375,1.05 1000:20:1020:tw,1.02,1.05 \
  1000:21:1021:tw,1.03,1.05,1.075,1.1 1000:4:1004:tw,1.01,1.02 1000:6:1006:tw,1.02,1.25 1000:8:1008:tw,1.02,1.05,1.75
python3 gsm_bind.py J9 2000:9:2009:tw,tw1000,1.02,1.03,1.05:50000      # likewise J10 (…,1.0075,1.02,1.05,2.0), J11, J13, J14 (…,1.035,…),
                                                                       # J15, J17, J19 (…,1.0375,…), each with :50000
python3 gsm_bind.py J16ref 2000:16:2016:tw,tw1000,1.02,1.03,1.05 2000:16:2016:tw,tw1000,1.02,1.03,1.05:50000
python3 gsm_bind.py J21ref 2000:21:2021:tw,tw1000,1.03,1.05,1.075,1.1 2000:21:2021:tw,tw1000,1.03,1.05,1.075,1.1:50000
GSM_XMIN=0.045 python3 gsm_bind.py J21x 2000:21:2021:tw,tw1000,1.05:50000 2000:21:2021:tw1000:20000   # §2 diagnostic
python3 gsm_bind.py OO468 2000:4:2004:tw,tw1000,1.01,1.02 2000:6:2006:tw,tw1000,1.02,1.25 2000:8:2008:tw,tw1000,1.02,1.05,1.75
python3 gsm_ooid.py 4,6,8,10;  python3 gsm_hscan.py 9,10,11,13,14,15,17,18,19,20,4,6,8
python3 gsm_bind.py --hmono J16 16:1016/20000,2016/20000,2016/50000   # also J9, J21, seam, like (see cache/gsm_hmono_*.txt)
python3 gsm_table.py base1000 J9,J10,J11,J13,J14,J15,J16ref,J17,J19,J21ref,OO468 > cache/table.txt
```

## 7. Caveats
1. Everything at 2000 is K₂ (crystal integral). "K₂ = crystal integral" is checked at 1000 by direct rows, not at 2000.
2. Δ_ref = 5·10⁴ is new to the project. Its records pass the phase-strip assertion, give finite Laplace series and shadow residuals
   ≤ 6·10⁻⁸, and reproduce the 2·10⁴ results at nodes away from x₀ to ≤ 10–40 c₂. They are otherwise unvalidated, and there is no direct
   check at 5·10⁴.
3. The 1000 comparisons for J = 9–14 are not like-for-like (no own rows at 1000). The 1000 L values there are lower bounds on what own
   rows would give.
4. "λ-level" grades use NEARGAP2's `3|λ_K2 − λ_K1|`, a signed truncation estimate, as J = 8 / 12 / 16 were graded before.
5. The OO grades rest on Γ_OO = Γ_++ [NE], checked here at 2000 on 6 nodes per J and budgeted at ×230 its measured size.
