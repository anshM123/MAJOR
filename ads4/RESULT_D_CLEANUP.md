# Result D clean-up: direct small-`ν` nodes at `Δ_g = 300` for `J = 6–15`, `H_edge` on the remaining low-`J` rows, `J = 48` check

*22 Nov 2026, clean-up agent. Unrefereed; not committed by this agent (intermediate drafts of this file were picked up by other commits). No existing document, code file or cache was edited: new code
`ads4/code/spinfunc/cl_*.py` / `cl_*.sh`, new caches `spinfunc/cache/cl_*` only (a write guard, `cl_guard.py`, refuses any write to a
non-`cl_` cache file in every `cl_` process). Background: RESULT_D.md (update bullets), HYPOTHESES_CLOSE.md (HC), HYPOTHESES_CLOSE2.md
(HC2), PARITY_ODD_CLOSE.md (PO), REVIEW_EE_2026-11.md (the review), RESULT_D_STATEMENT.md.*

**Fixed and not tuned on anything found here:** the functional `Ω^EE` (min-norm gauge, NF02); the z4 protocol (direct n3 nodes, non-local
min-norm `X` below `ν = 25` for even `J`, `SM_CURV=row`, Jacobi fix, estimator `best`, value held below the lowest node, Hampel despike);
for odd `J` the row-injection protocol of PO §3 (z4's non-local term does not act on odd `J`); for the edge, the HC protocol (direct n1
nodes `x = 0.8 … 0.9`, exempt from the Hampel rule, natural spline in `x` as the protocol's row model). Two protocol choices of this note
are stated up front: **(i)** a direct node whose state RQMC error is `> 10 %` (`max|err|/max|E|`) is dropped (this removes the four
`x = 0.9` nodes at `Δ = 336`/`200`, where `ν` is within 10–15 % of `Δ`, and the `J = 5` `x = 0.95` nodes, exactly the class HC excluded by
hand); **(ii)** every chain is run on a dense near-gap mass grid (twist edge, then `1.015 … 1.12` in steps of `0.005`, `1.13 … 1.2` in
steps of `0.01`, then `1.25 … 5`), because odd `J ≥ 11` and even `J ≥ 12` have no exact-skeleton states and the older sparse grids missed
the minimum (`J = 15`: PO's `+0.829` at `1.0586Λ` is not the minimum; it is `+0.789` at `1.07Λ`).

Grades: **[P]** proved, **[D]** derived, **[DV]** derived + checked, **[NE]** numerical evidence, **[NE + M]** direct crystal values in the
chain model, **[M]** model, **[I]** open.

## 0. Summary

### 0.1 Task 1 — direct small-`ν` nodes at `Δ_g = 300`, `J = 6–15` (near-gap states `m ≤ 1.3Λ`, twist edge included)

Old = the superseded convention A (`K₂` + natural spline below `x = 0.05`, local `ψ₂` weights; HC's sweep `hc_chain_z4_t3_all`).
New = direct n3 nodes `ν ≈ 8.37 / 9.96 / 11.67 / 13.35 / 16.73 / 19.92` (`x = 0.021 … 0.05`) on the row that carries the near-gap states,
state and reference `2·10⁴` on common random numbers. `n` = number of negative eigen-directions over all states of that `J` (gap to `5Λ`).
Spread = jackknife over the 8 common-seed replicas, then the range over the estimator / Hampel / lowest-node variants.

| `J` | block | old (convention A) min [`m/Λ`] | **new min** [`m/Λ`] | jackknife | variant range (near-gap min) | `n` | edge shift (§0.2) | grade |
|---|---|---|---|---|---|---|---|---|
| 5 (odd, binding) | `+−` | `+0.2038` [1.0166] (HC: spot nodes, row 400 only) | **`+0.2038`** [1.0166] (rows 200 + 400 direct; `0.20384 → 0.20377`) | `±0.0000` | `+0.2038` | 0 | `0.0000` | [NE + M] |
| 6 | EE / OO | `+0.307` [1.25] | **`+0.3075` / `+0.3078`** [1.25] (PO; reproduced bit-for-bit) | `±0.0001` | `+0.307 … +0.308` | 0 | `0.0000` | [NE + M] |
| 7 (odd) | `+−` | `+0.704` [1.10] (same grid `+0.7025` [1.095]) | **`+0.7033`** [1.095] (rows 200 + 400 direct; row 400 only: `+0.7024`) | `±0.0002` | `+0.7033` | 0 | `0.0000` | [NE + M] |
| 8 | EE | `+0.426` [1.052] | **`+0.4271`** [1.055] | `±0.0009` | `+0.4271` (dir/cv/no-Hampel); conv. A same grid `+0.4258` | 0 | `−0.0003` | [NE + M] |
| 8 | OO | `+0.5408` [1.08] (same grid); all `m` `+0.4774` [1.65] (HC) | **`+0.5493`** [1.085] near gap; all `m` `+0.4787` [1.65] (LO base) | `±0.0006` | `+0.549 … +0.550`; uniform `Γ_c = −8000` below `ν = 8.37`: `+0.520` | 0 | — | [NE + M] |
| 9 (odd) | `+−` | `+0.4375` [1.03] | **`+0.4561`** [1.03] (rows 200 + 400 direct; row 400 only: `+0.4448 ± 0.0003`) | `±0.0006` | `+0.4561 … +0.4562` | 0 | `−0.0055` | [NE + M] |
| 10 | EE | `+0.6631` [1.28] (all `m` `+0.6026` [2.0]) | **`+0.6697`** [1.28] (all `m` `+0.6033` [2.0]) | `±0.0009` | `+0.6697` | 0 | `+0.0001` | [NE + M] |
| 10 | OO | `+0.6594` [1.28] (all `m` `+0.602`) | **`+0.6751`** [1.28] (all `m` `+0.6038` [2.0]) | `±0.0009` | `+0.674 … +0.676`; `Γ_c = −8000`: `+0.648` | 0 | — | [NE + M] |
| 11 (odd) | `+−` | `+0.6447` [1.036, twist edge] | **`+0.6728`** [1.036] | `±0.0024` | `+0.6728 … +0.6737` | 0 | `−0.0041` | [NE + M] |
| 13 (odd) | `+−` | `+0.7280` [1.10] | **`+0.7764`** [1.10] | `±0.0047` | `+0.7764 … +0.7788` | 0 | `−0.0018` | [NE + M] |
| 15 (odd) | `+−` | `+0.6304` [1.07] (same grid; HC sweep `+0.630` [1.068]) | **`+0.7892`** [1.07] (PO: `+0.829` at 1.0586, sparse grid) | `±0.015` | `+0.7892 … +0.7892` | 0 | `−0.0038` | [NE + M] |

Odd `J = 5, 7, 9`: the near-gap states lie between the `Δ = 200` and `400` rows (weight of the 200 row `≈ 0.3` at the gap), so both rows got direct nodes (§1.4).

**Verdict (task 1).** Every state of even `J = 6, 8, 10` and odd `J = 7, 9, 11, 13, 15` at `Δ_g = 300` is positive with direct small-`ν`
nodes, with jackknife errors `≤ 0.015` and variant spreads `≤ 0.003`; no negative direction anywhere (`n = 0`). Convention A was again
**pessimistic** wherever the small-`ν` window matters (`J = 9 / 11 / 13 / 15`: `+0.019 / +0.028 / +0.048 / +0.159`), and irrelevant where
it does not (`J = 6, 7, 8, 10`: `|Δ| ≤ 0.007`, minima away from the gap or `N_s` large). Every direct node passes the scalar check
against the exact Mellin object (`≤ 1.2σ`; §1.3). The binding state at `Δ_g = 300` is unchanged: **`J = 5`, `+0.2038`** (hybrid, twist edge
`1.0166Λ`), now with direct small-`ν` nodes on both rows it draws on and direct edge nodes. The next-lowest states are `J = 4`
(`+0.2137`, HC), `J = 6` (`+0.3075`), `J = 8` EE (`+0.4271`) and `J = 9` (`+0.4561`).

### 0.2 Task 2 — `H_edge` on the remaining rows (shift of the near-gap minimum, `N_s` units)

Direct n1 nodes `x = 0.8, 0.825, 0.85, 0.875, 0.9` (adaptive proposal, shared seeds, `2¹¹·2·8` points) on every row carrying a near-gap
state; "shift" = (with edge nodes) − (same chain, same small-`ν` nodes, held edge); "max" = largest |shift| over the near-gap states.

| `Δ_g` | `J` (rows) | min without → with edge [`m/Λ`] | shift of min | max shift near gap | edge jackknife |
|---|---|---|---|---|---|
| 300 | 4 (400) | `+0.2137 → +0.2137` [1.0133] | `0.0000` | `0.0000` | — |
| 300 | **5** (200/400/800/1600, `x ≤ 0.92`) | **`+0.2038 → +0.2038`** [1.0166] | `0.0000` | `< 5·10⁻⁵` | `±0.0000` |
| 300 | 6 (400) | `+0.3075 → +0.3075` | `0.0000` | `0.0000` | — |
| 300 | 7 (200/400) | `+0.7033 → +0.7033` | `0.0000` | `0.0004` | — |
| 300 | 8 (400) | `+0.4271 → +0.4268` | `−0.0003` | `−0.0006` | — |
| 300 | 9 (200/400) | `+0.4561 → +0.4506` | `−0.0055` | `−0.0056` | `±0.0006` |
| 300 | 10 (400) | `+0.6697 → +0.6698` | `+0.0001` | `+0.0017` | — |
| 300 | 11 (400) | `+0.6728 → +0.6687` | `−0.0041` | `−0.0041` | `±0.0022` |
| 300 | 12 (400) | `+0.6712 → +0.6706` [1.055] | `−0.0006` | `−0.0006` | — |
| 300 | 13 (400) | `+0.7764 → +0.7746` | `−0.0018` | `−0.0053` | `±0.0046` |
| 300 | 14 (400) | `+0.7495 → +0.7487` [1.0458] | `−0.0008` | `−0.0008` | — |
| 300 | 15 (336/400) | `+0.7892 → +0.7853` [1.07] | `−0.0038` | `−0.0156` [1.049] | `±0.015` |
| 300 | 18 (336/400) | `+0.9104 → +0.9224` [1.08] | `+0.0120` | `+0.0201` | `±0.0012` |
| 300 | 20 (400) | `+0.8631 → +0.9051` [1.0648] | `+0.0420` | `+0.0420` | `±0.0036` |
| 300 | 22 (336) | `+0.8411 → +0.8440` [1.0711] | `+0.0029` | `+0.0029` | — |
| 1000 | 4 (400/1600) | `+0.2041 → +0.2041` [1.004] | `0.0000` | `0.0000` | — |
| 1000 | 6 (400/1600) | `+0.2296 → +0.2296` [1.006] | `0.0000` | `0.0000` | — |
| 1000 | 7 (200–1600) | `+0.2684 → +0.2684` [1.007] | `0.0000` | `0.0000` | — |
| 1000 | 8 (400/1600) | `+0.4304 → +0.4309` [1.008] | `+0.0005` | `+0.0006` | — |
| 1000 | 9 (200–1600) | `+0.4546 → +0.4542` [1.03] | `−0.0005` | `−0.0016` | — |
| 1000 | 10 (400/1600) | `+0.5853 → +0.5886` [1.01] | `+0.0032` | `+0.0032` | `±0.0006` |
| 1000 | 11 (400/800) | `+0.3491 → +0.3482` [1.011] | `−0.0008` | `−0.0008` | — |
| 1000 | 12 (400/800/1600) | `+0.6755 → +0.6759` [1.055] | `+0.0004` | `+0.0069` | `±0.0000` |
| 1000 | 13 (400/800) | `+0.4109 → +0.4076` [1.013] | `−0.0033` | `−0.0039` | `±0.0001` |
| 1000 | 14 (400/800/1600) | `+0.6973 → +0.6982` [1.035] | `+0.0009` | `+0.0016` | — |
| 1000 | 16 (400/800/1600) | `+0.7221 → +0.7191` [1.025] | `−0.0029` | `−0.0064` | `±0.0001` |

**Verdict (task 2).** With direct edge values on every row that carries a near-gap state at `Δ_g = 300` (`J = 4–22`) and on the low-`J`
rows at `Δ_g = 1000` (`J = 4–16`), every state stays positive; the largest unfavourable shift of any near-gap state is **`−0.016 N_s`**
(`J = 15`, twist edge, which is not its minimum), the largest unfavourable shift of a per-`J` minimum `−0.0055` (`J = 9`, 300). The two
sizeable shifts are favourable (`J = 20`: `+0.042`, `J = 18`: `+0.012`), the same pattern as HC's `J = 31`. Together with HC (1000:
`J = 5, 17, 31, 36, 40`; 300: `J = 5, 16, 24`) `H_edge` is now computed on **every row that carries a near-gap state at `Δ_g = 300` for
`J ≤ 24`** and at `Δ_g = 1000` for `J ≤ 17` and `J = 31, 36, 40`.

### 0.3 Task 3 — `J = 48` at `Δ_g = 1000` against direct values

**The consistency check failed, the state still passes.** Direct n1 nodes on `(1048, 48)` (fw structure basis, shared seeds) disagree with the
fixed-engine `K₂` row that FOLD_WINDOW §4 used: they agree at `x ≤ 0.2` (C2 `|Δ| ≲ 80 c₂`, C3H `≤ 180 c₂`), but from
`x ≈ 0.22` on `K₂` overshoots the `P_J` lobes by ×2–4 (`x = 0.4–0.55`: `|Δ| ≈ 400–750 c₂`) and is out of phase at `x ≥ 0.7` (`|Δ|` up to
`5·10³`, `1.9·10⁴` at `x = 0.8` `+−`) — the pattern HC found at `(336, 28)` (`J²/m = 2.35`; here `2.2`). So **`K₂` is not accurate on this
row**, against EE_DG1000_LARGEJ's criterion `J²/m ≲ 2.2`. With a complete direct sub-row (`x = 0.125 … 0.8`, step `0.025`, 27 nodes) in the
unchanged J = 48 chain: **min `λ/N_s = +0.915 ± 0.08`** (jackknife) at `1.15Λ` (was `+0.932`), `n = 1`; variants `+0.67` (natural spline,
the aliasing variant, was `+0.73`) `… +0.99` (`K = 6`); max `|dF_hh|/N_h` `0.37` (was `0.065`). [NE + M].

### 0.4 Side findings
* **Failed edge nodes.** At `Δ = 336` (`J = 15, 18, 22`) and `Δ = 200` (`J = 9`) the `x = 0.9` node (`ν = 0.9 m ≈ 0.89 Δ`) has RQMC errors
  `23–121 %`. Kept, a single such node moves the `J = 18` gap state to **`−6.50`** (one node at `Γ ≈ −5.8·10⁶ c₂`) — an artefact of the
  failed point, not a measurement; the HC rule (drop it) is now automatic (`CL_MAXREL = 0.1`). The same class at `J = 5`, `x = 0.95`
  inflates the `J = 5` jackknife to `±0.05` without moving the central value beyond it.
* **Twist-edge states at `Δ_g = 1000`.** The low-`J` chains at 1000 had never sampled the twist-edge masses. They lower several per-`J`
  minima (old → new, both without edge nodes): `J = 4` `+0.2103 → +0.2041` (1.004Λ), `J = 6` `+0.2508 → +0.2296`, `J = 7` `+0.3114 →
  +0.2684`, `J = 11` `+0.4099 → +0.3491`, `J = 13` `+0.6046 → +0.4109`, `J = 16` `+0.7477 → +0.7221`. All remain above the binding state
  (`J = 5`, `+0.2021` at 1.01Λ; `+0.2009` hybrid at 1.0048Λ, EE_DG1000_LARGEJ), so `n = 1` at 1000 is unaffected, but the margins quoted
  per `J` in RESULT_D_STATEMENT were taken on a grid that started at `1.01Λ`.
* **`K₂` at `x = 0.05` (convention A's anchor) vs direct**, odd `J`: `C3H` `K₂` `−263 / −400 / −690` vs direct `−47 / −63 / −145 c₂`
  (`J = 9 / 11 / 15`); `C2` `−21 / −37 / −333` vs `+34 / +32 / +84`. Convention A's pessimism comes from this node and the spline below it.

## 1. Task 1: methods, checks, details
### 1.1 Points
* Even `J = 8, 10`: `CL_MODE=oo cl_points.py J 400 0.021,0.025,0.0293,0.0335,0.042,0.05` = `oo_points.py` (unchanged code path:
  `hc_points.setup(odd=True)` forced, 3×3 `[EE(2), OO(1)]` p8 integrand, n3 hybrid contour, `ν_f = 13` below 13, `2¹⁰·2·8` points,
  reference `2·10⁴` on CRN) → `cache/cl_oopts_{EE,C3H}_J<J>_D400.pkl`. One file per state.
* Odd `J = 7, 9, 11, 13`: `CL_MODE=hc HC_METHOD=n3 cl_points.py J 400 …` = `hc_points.py` n3 path (odd 1×1 integrand) →
  `cache/cl_pts_*_J<J>_D400.pkl`; plus the `Δ = 200` row for `J = 9, 5, 7` (`x = 0.042, 0.05, 0.059, 0.067, 0.084, 0.1`, `ν = 8.33 … 19.84`).
* `J = 6` (oo_pts, PO) and `J = 15` (hc_pts, PO) were not recomputed; their chains were re-run here on the dense mass grid.

### 1.2 Chains
* Even `J`: `cl_teven.sh J 400 oo` runs `oo_chain.py` (via `cl_chain.py oo`, unchanged source, redirections only) with the oo_t.sh variant
  set: `best` (main), `dir`, `cv`, no Hampel, `zero`/`lin` below the lowest node, convention A, uniform `Γ_c = −1000 … −8000` below
  `ν = 8.37`, jackknife `jk0–7`. run0 is the EE block (z4 with the even 2×2 block of the 3×3 points), run1 the OO block.
* Odd `J`: `cl_todd.sh J 400 …` runs `hc_chain.py z4 300` (via `cl_chain.py hc`) with `HC_INJ=J:400`, `HC_IXMAX = 0.2` (the small-`ν`
  nodes only): `inj` (main, Hampel on, as PO), `inj_nodesp`, `inj_drop` (the `ν = 8.37` node dropped), `convA` (no injection), `injjk0–7`.
* Reproduction: `cl_chain.py` reproduces PO's `J = 15` (`hc_chain_z4_po15_inj`) and `J = 6` (`oo_chain_po6_best`) outputs line by line
  (identical per-`J` blocks), and the EE-block path through `hc_chain z4` with the 3×3 points (`CL_Z3OO`) gives PO's `J = 6` `+0.3075`.

### 1.3 Node checks [NE]
Scalar check against the exact Mellin object (z3 `scalar_check`, `cache/cl_check.txt`), every small-`ν` node of `J = 7–13`: `Δc₂` within
`0.2–1.2σ` (`±200–600 c₂` at `ν = 8.37`, `±11–90` above); graviton RQMC errors `7.6–8.0 %` at `ν = 8.37`, `1.6 %` at 9.96, `≤ 0.6 %` above
(the same as PO/HC2). Direct `Γ` values (`c₂`, curvature subtracted), e.g. OO `Γ_c` at `J = 8`: `−166 −96i, −24, +6, +21, +33, +32`;
`J = 10`: `−311 −151i, −45, +1, +23, +37, +33` (`ν = 8.37 … 19.92`); odd `C3H` at `J = 11`: `+639, +188, +73, +5, −64, −63`.

### 1.4 The `Δ = 200` row for odd `J = 5, 7, 9`
The odd rows of `J ≤ 9` are `Δ = 200, 400, 800, 1600`; a near-gap state at `m ≈ 1.02–1.1Λ` (`m ≈ 305–330`) is interpolated in `1/m`
between the 200 row (`m ≈ 198`) and the 400 row (`m ≈ 398`), with weight `≈ 0.3` on the 200 row at the gap. The 200 row's lowest `K₂` node
`x = 0.05` sits at `ν ≈ 9.9`, deep in the region where `K₂` fails. Direct values there (`J = 9`, C2 / C3H): `x = 0.05`: `−75 / +131` vs
`K₂` `−875 / −434`; `x = 0.1` (`ν ≈ 20`): `−1 / −24` vs `−74 / −235`. Every 200-row node passes the scalar check (`≤ 1.1σ`). Effect on the
near-gap minimum (both rows direct vs 400 only): `J = 9` `+0.4448 → +0.4561`, `J = 7` `+0.7024 → +0.7033`, `J = 5` `+0.20377` (both). The
`Δ = 200` row with direct nodes is the main value in §0.1.

## 2. Task 2: method and notes
* Points: `CL_MODE=hc HC_METHOD=n1 HC_ADAPT=1 HC_M=11 HC_MREF=11 cl_points.py J D 0.8,0.825,0.85,0.875,0.9` (= `hc_points.py` n1 path;
  even `J`: w0 2×2 EE integrand + C3H; odd `J`: 1×1). Errors `≤ 1 %` except the dropped `x = 0.9` nodes listed in §0.4.
* Chains: 300 — z4 (`cl_teven.sh … edge`: `HC_INJ=J:rows HC_IXMIN=0.7 HC_DESPIKE_INJ=0` against the same chain without it; odd `J`:
  `cl_todd.sh … edge`, small-`ν` + edge vs small-`ν` only); 1000 — the EE_NEARGAP main chain `hc_chain w2 1000` (`cl_t1000.sh`), all
  edge rows of the `J` injected. Variants: Hampel on the injected nodes (`edge_desp`; changes the near-gap states by `≤ 0.005`, except `J = 13`: `0.009` and
  `J = 20`: `0.020`, min `+0.885`); edge jackknife where a shift exceeds `0.003`.
* **Code fix needed for even `J < 11` (in the wrapper, not in hc_chain):** hc_chain injects `C3H` edge nodes only through `q3._q2_even`,
  which the w2 loader uses for `J ≥ 11`; for even `J < 11` the `C3H` rows come from `W2._q3_load` (l9 caches), so `C3H` edge nodes would
  silently not be injected. `cl_chain.py` injects them there too (`CL_C3HLOW=1`); the injection log lists both kinds for every run.
* Edge node values grow towards the edge and oscillate with the `D_h` period, as in HC: e.g. `(400, 20)` C2 `++` `+988, −1021, −1052,
  +3151, −2721` and `+−` `−11, +1685, −297, −3419, +3972` (`x = 0.8 … 0.9`; held value `+1069 / −623`); `(400, 8)` `++` `−200 … +640`.
* Representation: the protocol's natural spline through the nodes (spacing `0.025`, at least ~4 nodes per half-period at `J ≤ 22`);
  HC's shape-fit edge models exist only in ee1000 mode and were not needed here. Grade of the shifts: [NE + M].

## 3. Task 3: `J = 48` at `Δ_g = 1000`
* **Points.** `cl_points.py 48 1048 …` with `CL_FWPATCH=1` (adds `fw_frame.patch_conserved_clean`, the structure basis of the fw rows),
  n1 adaptive, `2¹²` (4 check nodes) / `2¹¹` / `2¹⁰` points per replica, `R = 8`; relative errors `0.03–1.7 %`. The `x = 0.475` node is left
  out: its `2·10⁴` reference has `E₀₀/E₀₁ = 82` against `1.3–1.4` at every other node (the reference matrix is nearly singular at that
  `ν`; the chain's off-diagonal becomes non-finite there). The `K₂` node at `0.475` stays (removing it as well changes nothing, `+0.9148`).
* **Chain.** `cl_fw48.py` runs `fw_chain48.py` unchanged (ee1000 shape chain `K = 4`, `x ∈ [0.05, 0.8]`, hybrid states, κ image); the
  direct nodes enter through fw_gamma's own direct-node path (replacing the `K₂` node at the same `ν`); the generalised helicity basis is
  built from the pure `K₂` order-0 records (as hc_chain builds it without injected nodes; with the direct nodes in the basis window the
  basis fit fails). Reproduction: `FW_NODIRECT=1` gives `+0.9319`, FOLD_WINDOW's number.
* **Node comparison** (`cache/cl_j48.txt`, `Γ` in `c₂`, curvature image subtracted, C2 `++` direct / `K₂`): `x = 0.2`: `+237 / +212`;
  `0.3`: `−214 / −402`; `0.4`: `−250 / +179`; `0.45`: `−14 / +587`; `0.5`: `+408 / +1097`; `0.6`: `+1424 / +1142`; `0.7`: `+2186 / +3248`;
  `0.75`: `+2091 / +5203`; `0.8`: `−6600 / −2328`; C2 `+−` at `0.8`: `−5071 / +14099`. Jackknife errors `30–60 c₂` (`≤ 0.6`), `120–300` above.
* **Per state** (main): gap `1.047Λ` `+0.979`, `1.1Λ` `+0.927`, `1.15Λ` `+0.915`, `1.2Λ` `+0.918`; jackknife replicas `+0.870 … +0.960`.
* Partial replacement is again meaningless: the four check nodes alone give `+0.340` (gap), as HC warned for oscillating rows.

## 4. Regradings
| item | before | now |
|---|---|---|
| **`Δ_g = 300`, even `J = 6, 8, 10`, odd `J = 7, 9, 11, 13, 15` (review S4)** | [NE + M] in convention A (superseded small-`ν` treatment), `J = 12–16`, 18–22, 4, 5 already direct | **[NE + M] in the validated small-`ν` treatment for every even `J = 4–22` and odd `J = 5–15`**, near gap and to `5Λ`; min `+0.2038` (`J = 5`). Caveats unchanged: below the lowest direct node (`ν ≈ 8.3`) the value is held [M] (stress: for `J = 8, 10` a uniform `Γ_c = −8000 c₂` there still leaves `+0.52 / +0.65`); odd `J` use the local de-shift (z4's non-local term is wired for even `J` only) [M]. |
| `H_OO`, `J = 8, 10` | [NE + M] with `K₂` (convention A) rows | **[NE + M] with direct small-`ν` OO nodes** (min `+0.4787` [1.65Λ, LO base] / near gap `+0.549`; `+0.6038` / `+0.675`) |
| **`H_edge`** | computed on 1000 `J = 5, 17, 31, 36, 40`, 300 `J = 5, 16, 24`; elsewhere an assumption [M] | **[NE + M] on every row carrying a near-gap state at 300 for `J = 4–22` (and 24, 16, 5 from HC)** and at 1000 for `J = 4–16` (+ 17, 31, 36, 40); largest unfavourable change of a near-gap state `−0.016 N_s`, of a per-`J` minimum `−0.006`. Still [M] on the rows not computed (§5). |
| `J = 5` binding state at 300 | `+0.2038` (spot nodes on row 400 only) | `+0.2038` with direct small-`ν` nodes on rows 200 + 400 and direct edge nodes on all four rows (`x ≤ 0.92`) [NE + M] |
| per-`J` minima at `Δ_g = 1000`, low `J` | taken on a mass grid from `1.01Λ` | the twist-edge states lower `J = 4, 6, 7, 11, 13, 16` by `0.006–0.19` (§0.4); all `≥ +0.204`, binding state unchanged |
| `J = 48` at `Δ_g = 1000` | [NE + M], `+0.932` on the fixed-engine `K₂` row (FOLD_WINDOW §4) | **[NE + M] with a direct sub-row**: `+0.915 ± 0.08` (`+0.67 … +0.99`); the `K₂` row itself is regraded **not accurate** at `x ≳ 0.22` (errors `10²–10⁴ c₂`), so `K₂` rows with `J²/m ≈ 2.2` are not usable at the `10² c₂` level [NE] |

Nothing here changes the implication [D] or the bounds `3.13·10⁻⁶` (300) / `2.89·10⁻⁸` (1000); every computed state is positive and the
binding states are unchanged (`J = 5`: `+0.2038` at 300, `+0.2021` / `+0.2009` (hybrid) at 1000).

## 5. What remains open
1. **`Δ_g = 300`, `ν < 8.3` below the lowest direct node** [M]: held in every chain. Irrelevant at `J ≤ 22` (failure needs a uniform
   `Γ_c ≲ −10⁴ c₂` there; `J = 8, 10` tolerate `−8000` with margin `≥ 0.52`), decisive at `J ≥ 28` (HC2 §3: `J = 32` the most fragile
   state, no validated method below the fold) — unchanged by this note.
2. **Odd `J` small-`ν` treatment uses the local de-shift** [M]; the non-local min-norm term is implemented for even `J` only (PO §3).
3. **`H_edge` not computed** [M]: `Δ_g = 300` rows of `J = 28, 32` (`Δ = 336`, where `x = 0.9` is within 10 % of `Δ` and fails — the
   edge window there needs a larger-`Δ` row or a different contour) and every `J ≥ 26` row; `Δ_g = 1000` even `J = 18–32`, odd
   `J = 19–29` and `J = 48` (the ee1000 rows `Δ = 1000 + J`).
4. **`K₂` rows near `J²/m ≈ 2` at `Δ_g = 1000`.** `J = 48` (`2.2`) is now on a direct sub-row; `J = 40` (`1.5`) was validated by a
   25-node direct sub-row in EE_DG1000_LARGEJ (≤ `0.05 N_s`, old engine); the even rows in between (`J = 42–46`) have no rows at all
   (`H_LJ`). Whether `K₂` is adequate at `J²/m ≈ 1.6–2` at the `10² c₂` level is open [I].
5. Unchanged from HC2/RESULT_D_STATEMENT: the `ν < 5` window at 300 for `J ≥ 28`, a derivation of `Γ_OO ≡ Γ_++`, `H_LJ`, `H_R`, the
   analytic inputs; `Δ_g = 100` not established.

## 6. Reproduce
All from `ads4/code/spinfunc/`, `OMP_NUM_THREADS=1`, `ulimit -v 3000000`; at most two busy processes (analyses go through `cl_pause.sh`,
which SIGSTOPs one point worker and serialises analyses with a lock). Every `cl_` process loads `cl_guard.py` (no write outside `cache/cl_*`).
* **Points** (resumable; finished nodes skipped): `./cl_queueA.sh` (small `ν`, `J = 8, 9, 10, 11, 7, 13`, row 400), `./cl_queueA2.sh`
  (row 200, `J = 9, 5, 7`), `./cl_queueB.sh` (edge rows, 300 then 1000), `./cl_queueC.sh` (`J = 48` check nodes), `./cl_queueC2.sh` +
  `./cl_queueC3.sh` (`J = 48` sub-row). Single node set: `CL_MODE=oo python3 cl_points.py 8 400 0.021,0.025,0.0293,0.0335,0.042,0.05`;
  `CL_MODE=hc HC_METHOD=n3 python3 cl_points.py 9 400 …`; edge `CL_MODE=hc HC_METHOD=n1 HC_ADAPT=1 HC_M=11 HC_MREF=11 python3
  cl_points.py 12 400 0.8,0.825,0.85,0.875,0.9`. Caches `cache/cl_oopts_*`, `cache/cl_pts_*` (+ `cl_*raw*`), logs `cache/cl_q*.log`.
* **Chains** (each run 0.5–1 min): even `J`: `./cl_teven.sh 8 400 oo edge` → `cache/cl_oochain_cl8_<variant>.*` (OO + EE block) and
  `cache/cl_chain_z4_cl8_{ebase,edge,…}.*`; `J = 18–22`: `CL_Z3SRC=hc ./cl_teven.sh 18 336,400 edge`; odd `J`: `./cl_todd.sh 11 400 400
  small edge` (`J = 9, 7` with both rows: `cl_todd_r2.sh 9 200,400 200,400 small edge`, generated by `cl_driver5.sh`); `J = 5`:
  `cl_driver5.sh` (tags `cl5r2_*`) and `./cl_t5b.sh`; 1000: `./cl_t1000.sh 13 400,800` → `cache/cl_chain_w2_k13_*`. The drivers
  `cl_driver{,2,3,4,5,1000}.sh` ran all of these in the order of this note.
* **Summaries:** `python3 cl_sum.py oo cl8 best dir cv nodesp zero lin convA c-1000 c-2000 c-4000 c-8000`;
  `python3 cl_sum.py hc z4 300 cl9r2 inj convA inj_nodesp inj_drop --jk=injjk --ref=convA --J=9`;
  `python3 cl_sum.py hc z4 300 cl20 edge ebase edge_desp --jk=edgejk --ref=ebase --J=20`;
  `python3 cl_sum.py hc w2 1000 k13 edge base edge_desp --jk=edgejk --ref=base --J=13`; per-state lists with `--states=1`.
* **Node validation:** `python3 cl_check.py 8:400:oo 10:400:oo 9:400 11:400 7:400 13:400 9:200 5:200 7:200` → `cache/cl_check.txt`.
* **Reproduction checks:** `cl_chain.py` on PO's settings reproduces `hc_chain_z4_po15_inj` and `oo_chain_po6_best` (identical `J` blocks);
  `FW_NODIRECT=1 CL48_TAG=k2 python3 cl_fw48.py` reproduces `fw_chain48_k2` (`+0.9319`).
* **`J = 48`:** `python3 cl_j48.py` (node comparison → `cache/cl_j48.txt`); chain with the direct sub-row: `CL48_TAG=full python3 cl_fw48.py` (`CL48_JK=i` jackknife, `EE1000_K=5|6`,
  `EE1000_MODEL=spline`, `SM_CURV=none` variants; `FW_NODIRECT=1` = published chain) → `cache/cl_fw48_<tag>.txt`.
