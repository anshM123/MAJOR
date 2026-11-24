# GAP_SEAM2: J = 27, 29 at Δ_g = 2000 rerun with Δ_ref = 5·10⁴ on the 1000 x grid (J = 18, 20, 22, 24 not finished)

*24 Nov 2026, gap-seam follow-up agent. Unreviewed, not committed. No existing file edited. `NEARGAP_FILL.md` / `code/ngfill/` not touched.
Code: `ads4/code/gapseam2/` (copies of `gapseam/` scripts; `cache/.gitignore` ignores everything but itself; `gs2_guard` refuses
writes outside `gapseam2/cache`). Compute ≈ 2.1 CPU-h (two lanes, 12:08–13:09, `ulimit -v 3000000`, one thread). The run was stopped
early, so items 2 and 3 of the task were not done (§3).*

## 0. Table (EE block, odd J, parity-odd `+−` 1×1; value / sign-free L (λ-level L); ratios are 2000 / 1000 at fixed m/Λ)

Columns:
* "2·10⁴ dense" is GAP_NEARGAP's run (`gapnear/cache/gn_bind_J27/J29.json`).
* "**5·10⁴, 1000 grid ∪ edge**" (`J27e`) uses the 1027 grid, `x = 0.05, 0.0625, …, 0.15, 0.175, …, 0.8`, plus `0.825 … 0.9`.
* "5·10⁴, 1000 grid exactly" (`J27k`, `J29k`) uses the 1027 / 1029 grid and stops at x = 0.8. This is the like-for-like variant.

Grades follow GAP_SEAM.

| J, m/Λ | 1000 same protocol | 2000, 2·10⁴ dense | **2000, 5·10⁴, 1000 grid ∪ edge** | 2000, 5·10⁴, 1000 grid exactly | grade 2000 (best variant) | chart-A ratios E_trunc, E_trunc^λ, E_low, E_edge (exact-grid variant; edge-variant E_trunc / E_trunc^λ in brackets) |
|---|---|---|---|---|---|---|
| 27, tw 1.01342 | — | +0.9932 / 0.05 (0.41) | +0.9670 / 0.06 (**0.46**) | +0.9726 / 0.06 (0.25) | [NE + M] | — |
| 27, 1.0175 (new hybrid min, `+0.9125`) | — | not run | +0.9098 / 0.07 (0.56) | +0.9140 / 0.07 (0.27) | [NE + M] | — |
| 27, tw1000 1.02667 | +1.0041 / 0.0 (0.1) | +1.0239 / 0.06 (0.63) | +0.9983 / 0.08 (0.77) | +1.0006 / 0.08 (0.38) | [NE + M] | 0.260, 0.243, 0.246, 0.236 [0.287, 0.283] |
| 27, 1.03 | +1.0030 / 0.0 (0.1) | +1.0237 / 0.06 (0.69) | +0.9984 / 0.08 (0.85) | +1.0004 / 0.08 (0.43) | [NE + M] | 0.260, 0.244, 0.246, 0.236 |
| 27, 1.04 | +0.9926 / 0.0 (0.1) | +1.0158 / 0.07 (0.84) | +0.9920 / 0.10 (**1.07**) | +0.9933 / 0.10 (0.55) | [NE±] λ-level (edge variant only) | 0.260, 0.248, 0.246, 0.236 |
| 27, 1.05 | +0.9611 / 0.0 (0.2) | +0.9842 / 0.07 (0.99) | +0.9620 / 0.11 (1.30) | +0.9625 / 0.11 (0.73) | [NE±] λ-level (edge variant only) | 0.260, 0.250, 0.246, 0.236 |
| 27, 1.075 | +0.9967 / 0.0 (0.3) | +1.0196 / 0.09 (1.20) | +0.9989 / 0.15 (1.62) | +0.9989 / 0.15 (1.25) | [NE±] λ-level | 0.260, 0.252, 0.246, 0.236 |
| 27, 1.10 | +0.9715 / 0.1 (0.4) | +0.9918 / 0.11 (1.34) | +0.9740 / 0.20 (1.79) | +0.9738 / 0.19 (1.64) | [NE±] λ-level | 0.260, 0.253, 0.246, 0.236 |
| 29, tw 1.01440 | — | +0.9480 / 0.04 (0.63) | (edge nodes incomplete) | +0.9143 / 0.05 (**0.16**) | [NE + M] | — |
| 29, 1.015 (hybrid min, `+0.9190`) | — | not run | — | +0.9133 / 0.05 (0.16) | [NE + M] | — |
| 29, tw1000 1.02862 | +0.9991 / 0.0 (0.0) | +1.0278 / 0.04 (0.83) | — | +0.9960 / 0.06 (0.25) | [NE + M] | 0.272, 0.220, 0.246, 0.171 |
| 29, 1.04 | +0.9667 / 0.0 (0.1) | +0.9952 / 0.05 (0.95) | — | +0.9673 / 0.07 (0.37) | [NE + M] | 0.270, 0.222, 0.246, 0.161 |
| 29, 1.05 | +0.9823 / 0.0 (0.1) | +1.0103 / 0.05 (1.04) | — | +0.9836 / 0.08 (0.52) | [NE + M] on this grid (2·10⁴ dense: λ-level 1.04) | 0.269, 0.226, 0.246, 0.135 |
| 29, 1.06 | +0.9961 / 0.0 (0.1) | +1.0239 / 0.06 (1.11) | — | +0.9983 / 0.09 (0.70) | as above | 0.268, 0.231, 0.246, 0.125 |
| 29, 1.075 | +0.9837 / 0.0 (0.2) | +1.0103 / 0.07 (1.18) | — | +0.9867 / 0.11 (1.01) | [NE±] λ-level (thin) | 0.266, 0.239, 0.246, 0.109 |
| 29, 1.10 | +0.9888 / 0.0 (0.3) | +1.0138 / 0.08 (1.28) | — | +0.9926 / 0.14 (1.42) | [NE±] λ-level | 0.263, 0.248, 0.246, 0.088 |

J = 18, 20, 22, 24: **not run** (§3). Their Δ_g = 1000 same-protocol baselines are in `cache/gs2_bind_base1000.json`:
* J = 18 at tw: `+0.7465 / 0.1 (0.3)`;
* J = 20 at tw: `+0.9543 / 0.1 (0.2)`;
* J = 22 at tw: `+0.9786 / 0.1 (0.2)`;
* J = 24 at tw: `0.0 (0.1)`.

For J = 22 / 24 the baseline also has a 1.02 entry, which lies below the 1000 twist edge and is invalid.

## 1. Verdict

**Negative findings (prominent).**
* **N1. GAP_SEAM's prediction fails: the λ-level L does not double.**
  * J = 27 at the twist edges: tw 0.41 → 0.46, tw1000 0.63 → 0.77 (edge variant). On the exact 1000 grid the values are *lower*: 0.25 / 0.38.
  * J = 29 on the exact grid: tw 0.63 → 0.16, tw1000 0.83 → 0.25.
  * **The J = 27 and J = 29 twist edges stay [NE + M].**
  * The λ-level truncation step does fall as predicted: J = 27 tw1000 0.970 → 0.549; J = 29 0.939 → 0.548. But it is no longer the only large term:
    * **E_edge** (x > xmax, with `G_edge = 10 × max|Γ|` of the last four nodes, |Γ| ≈ 10⁴ c₂ at x = 0.85–0.9) dominates when the row stops at 0.8: 2.0–5.1 at the twist edges;
    * **E_low** rises ×3.7 (0.05 → 0.18–0.22), because the low window is `ν < x_first·m` and x_first moves from 0.025 to 0.05 (ν < 101 instead of 50).
* **N2. The "same x grid" prescription has a cost at 2000.**
  * It makes chart A exact for the truncation terms: E_trunc 0.26–0.27, E_trunc^λ 0.22–0.25, against 0.25 predicted. With 2·10⁴ on the dense grid they were 0.43–0.70.
  * It also raises E_low and, at x ≤ 0.8, E_edge, so the leverage *drops*. Grades depend on the variant:
    * J = 27 is λ-level [NE±] from 1.04Λ with the edge nodes, but only from 1.075Λ on the exact grid;
    * J = 29 at 1.05 / 1.06Λ moves from λ-level [NE±] (2·10⁴ dense, L 1.04 / 1.11) to [NE + M] (exact grid, L 0.52 / 0.70).
* **N3. The values shift with the reference.** Every 2000 value moves by −0.02 … −0.03 (J = 27) and −0.03 (J = 29) between 2·10⁴ and 5·10⁴, toward the 1000 values. Examples:
  * J = 27 tw1000: +1.0239 → +0.9983 (1000: +1.0041);
  * J = 29 tw: +0.948 → +0.914.

  All stay positive, ≥ +0.91.
* **N4. New lowest J = 27 state at 2000.** The hybrid minimum is at 1.0175Λ inside the seam strip: `+0.9098` (hybrid `+0.9125`), not evaluated before. J = 29's minimum is at 1.015Λ: `+0.913`. Both are positive, [NE + M].

**Positive findings.**
1. No state is negative. The lowest value is `+0.910` (J = 27, 1.0175Λ).
2. **Chart A holds on every like-for-like budget term once the node set and the reference are matched.** On the 1000 grid with 5·10⁴ the ratios are:
   * J = 27: E_trunc 0.260, E_trunc^λ 0.243–0.253, E_edge 0.236;
   * J = 29: 0.263–0.272, 0.220–0.248;
   * E_low: 0.246, which is the x-window effect of N1.

   With 2·10⁴ on the dense grid E_trunc^λ was ×0.50. The diagnostic `J27x20` (2·10⁴, x = 0.025 dropped) gives E_trunc^λ 0.878, so the λ-step reduction is mostly the **reference** (0.878 → 0.549) and only partly the node (1.026 → 0.878). The sign-free E_trunc drop is mostly the **node** (19.5 → 13.8 → 11.7).
3. **H_mono-A** (`cache/gs2_hmono_J27.txt`). Γ(x) on 2027/5·10⁴ agrees with 1027 to ≤ 40 c₂ for x ≤ 0.6 and ≤ 240 c₂ at 0.7–0.8. The per-node |K₂ − K₁| steps are equal node by node (e.g. x = 0.2: 555 vs 546; 0.3: 358 vs 350; 0.7: 239 vs 93). With 2·10⁴ they were inflated to 2–2.8·10³ at x = 0.2–0.4.
   The per-node step in c₂ is Λ-independent, so the Λ⁻² scaling of E_trunc comes only from the weights. With a ×3 factor this is insufficient at J ≥ 27, which confirms GAP_SEAM N1: **J = 27, 29 need the batch-4 direct rows (or K₃) at their twist edges.**

## 2. Protocol
* The rows are GAP_SEAM's `gsm_rows` / `gsm_bind`, unchanged except as listed below. The state records of 2027 / 2029 are gapnear's (read-only).
  * New: state points at x = 0.0625, 0.0875, 0.1125, 0.1375, and every reference point at Δ_ref = 5·10⁴ on `k1000e`. That is 86 points for J = 27; for J = 29, 83 of 86, missing C3H x = 0.85–0.9, so J = 29 has no edge variant.
  * 37–50 s per point.
  * New shadows (5·10⁴, 27 / 29): residuals ≤ 9·10⁻¹⁴.
* `gs2_bind.py` = `gsm_bind.py` plus `GS2_XMAX`, which drops 2000-row nodes above XMAX. The exact variant uses `GS2_XMAX=0.8001`.
  * The x = 0.025 node is absent because no 5·10⁴ reference exists there.
  * `refbad` dropped nothing on 2027 / 2029 (x₀ = 1.09 / 1.01). At 1000 it drops x = 0.8 of 1029, as before.
* Masses: GAP_NEARGAP's, plus the 2000 hybrid minima from `gs2_hscan.py` (`cache/gs2_hscan.txt`):
  * J = 22 / 24 minima at their 2000 twist edges: +0.7795 / +0.8260;
  * J = 27 at 1.0175;
  * J = 29 at 1.015.
* Reproduction check: `gs2_bind` reproduces GAP_SEAM's J = 19 twist edge exactly (`+0.8552`, L 0.2 / 1.7; `cache/gs2_bind_chk19.txt`).

## 3. Not done / how to resume
* The lanes were stopped at 13:09 at the caller's request for an early hand-back.
  * Done before the stop: J = 18 refb, 1 of 6 points; J = 29 edge, 83 of 86 points.
  * Not started: J = 20 refb and rows, J = 22, J = 24.
* To resume (cached points are skipped):
  * lane A: `bash gs2_lane.sh cache/lane_A.jobs cache/lane_A.log`, then drop its J27 line;
  * lane B: likewise.
  * Planned grid for J = 18 / 20 / 22 / 24: `d05` (x = 0.05 … 0.9; no own 1000 rows exist, the 1000-side rows 336–1600 start at 0.05). About 140 points, ≈ 45–55 min per J.
* Then run the bind, e.g. `python3 gs2_bind.py J18 2000:18:2018:tw,1.0175,tw1000,1.02,1.05:50000` (J = 20: `…,1.015,…`).

## 4. Reproduce
```bash
cd ads4/code/gapseam2     # ulimit -v 3000000; OMP_NUM_THREADS=1; PYTHONDONTWRITEBYTECODE=1
python3 gs2_rows.py odd EE,C3H 27 2027 k1000e;  python3 gs2_rows.py odd EE,C3H 29 2029 k1000e
python3 gs2_bind.py base1000 1000:27:1027:tw,1.03,1.04,1.05,1.075,1.1 1000:29:1029:tw,1.04,1.05,1.06,1.075,1.1 (+ J18,20,22,24)
python3 gs2_bind.py J27e 2000:27:2027:tw,1.0175,tw1000,1.03,1.04,1.05,1.075,1.1:50000
GS2_XMAX=0.8001 python3 gs2_bind.py J27k <same>;  GS2_XMAX=0.8001 python3 gs2_bind.py J29k 2000:29:2029:tw,1.015,tw1000,1.04,1.05,1.06,1.075,1.1:50000
GSM_XMIN=0.045 python3 gs2_bind.py J27x20 2000:27:2027:<same masses>:20000          # diagnostic
python3 gs2_bind.py --hmono J27 27:1027/20000,2027/20000,2027/50000;  python3 gs2_table.py base1000 J27e,J27k,J29k > cache/table.txt
```
Caveats: as GAP_SEAM §7. Everything at 2000 is K₂, and Δ_ref = 5·10⁴ has no direct check.
