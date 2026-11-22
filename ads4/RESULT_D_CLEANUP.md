# Result D clean-up: direct small-`ν` nodes at `Δ_g = 300` for `J = 6–15`, `H_edge` on the remaining low-`J` rows, `J = 48` check

*22 Nov 2026, clean-up agent. Unrefereed, not committed. No existing document, code file or cache was edited: new code
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
| 6 | EE / OO | `+0.307` [1.25] | **`+0.3075` / `+0.3078`** [1.25] (PO; reproduced bit-for-bit) | `±0.0001` | `+0.307 … +0.308` | 0 | `0.0000` | [NE + M] |
| 7 (odd) | `+−` | `+0.704` [1.10] | **`+0.7024`** [1.095] | `±0.0001` | `+0.7024 … +0.7025` | 0 | `+0.0000` | [NE + M] ⟨r2⟩ |
| 8 | EE | `+0.426` [1.052] | **`+0.4271`** [1.055] | `±0.0009` | `+0.4271` (dir/cv/no-Hampel); conv. A same grid `+0.4258` | 0 | `−0.0003` | [NE + M] |
| 8 | OO | `+0.477` [1.65] (all `m`, HC) | **`+0.5493`** [1.085] near gap; all `m` `+0.4787` [1.65] (LO base) | `±0.0006` | `+0.549 … +0.550`; uniform `Γ_c = −8000` below `ν = 8.37`: `+0.520` | 0 | — | [NE + M] |
| 9 (odd) | `+−` | `+0.4375` [1.03] | **`+0.4448`** [1.03] | `±0.0003` | `+0.4448 … +0.4449` | 0 | `−0.0032` | [NE + M] ⟨r2⟩ |
| 10 | EE | `+0.6631` [1.28] (all `m` `+0.6026` [2.0]) | **`+0.6697`** [1.28] (all `m` `+0.6033` [2.0]) | `±0.0009` | `+0.6697` | 0 | `+0.0001` | [NE + M] |
| 10 | OO | `+0.6594` [1.28] (all `m` `+0.602`) | **`+0.6751`** [1.28] (all `m` `+0.6038` [2.0]) | `±0.0009` | `+0.674 … +0.676`; `Γ_c = −8000`: `+0.648` | 0 | — | [NE + M] |
| 11 (odd) | `+−` | `+0.6447` [1.036, twist edge] | **`+0.6728`** [1.036] | `±0.0024` | `+0.6728 … +0.6737` | 0 | `−0.0041` | [NE + M] |
| 13 (odd) | `+−` | `+0.7280` [1.10] | **`+0.7764`** [1.10] | `±0.0047` | `+0.7764 … +0.7788` | 0 | `−0.0018` | [NE + M] |
| 15 (odd) | `+−` | `+0.6304` [1.07] (same grid; HC sweep `+0.630` [1.068]) | **`+0.7892`** [1.07] (PO: `+0.829` at 1.0586, sparse grid) | `±0.015` | `+0.7892 … +0.7892` | 0 | `−0.0038` | [NE + M] |

⟨r2⟩: see §1.4 (the odd `J = 7, 9` near-gap states also draw on the `Δ = 200` row).

**Verdict (task 1).** Every state of even `J = 6, 8, 10` and odd `J = 7, 9, 11, 13, 15` at `Δ_g = 300` is positive with direct small-`ν`
nodes, with jackknife errors `≤ 0.015` and variant spreads `≤ 0.003`; no negative direction anywhere (`n = 0`). Convention A was again
**pessimistic** wherever the small-`ν` window matters (`J = 9 / 11 / 13 / 15`: `+0.007 / +0.028 / +0.048 / +0.159`), and irrelevant where
it does not (`J = 6, 7, 8, 10`: `|Δ| ≤ 0.007`, minima away from the gap or `N_s` large). Every direct node passes the scalar check
against the exact Mellin object (`≤ 1.2σ`; §1.3). The binding state at `Δ_g = 300` is unchanged: **`J = 5`, `+0.2038`** (hybrid, twist edge
`1.0166Λ`). The lowest new values are `J = 8` EE (`+0.427`) and `J = 9` (`+0.445`).

### 0.2 Task 2 — `H_edge` on the remaining rows (shift of the near-gap minimum, `N_s` units)

Direct n1 nodes `x = 0.8, 0.825, 0.85, 0.875, 0.9` (adaptive proposal, shared seeds, `2¹¹·2·8` points) on every row carrying a near-gap
state; "shift" = (with edge nodes) − (same chain, same small-`ν` nodes, held edge); "max" = largest |shift| over the near-gap states.

| `Δ_g` | `J` (rows) | min without → with edge [`m/Λ`] | shift of min | max shift near gap | edge jackknife |
|---|---|---|---|---|---|
| 300 | 4 (400) | `+0.2137 → +0.2137` [1.0133] | `0.0000` | `0.0000` | — |
| 300 | **5** (200/400/800/1600, `x ≤ 0.92`) | **`+0.2038 → +0.2038`** [1.0166] | `0.0000` | `< 5·10⁻⁵` | `±0.0000` |
| 300 | 6 (400) | `+0.3075 → +0.3075` | `0.0000` | `0.0000` | — |
| 300 | 7 (200/400) | `+0.7024 → +0.7024` | `0.0000` | `0.0003` | — |
| 300 | 8 (400) | `+0.4271 → +0.4268` | `−0.0003` | `−0.0006` | — |
| 300 | 9 (200/400) | `+0.4448 → +0.4417` | `−0.0032` | `−0.0032` | ⟨jk⟩ |
| 300 | 10 (400) | `+0.6697 → +0.6698` | `+0.0001` | `+0.0017` | — |
| 300 | 11 (400) | `+0.6728 → +0.6687` | `−0.0041` | `−0.0041` | ⟨jk⟩ |
| 300 | 12 (400) | `+0.6712 → +0.6706` [1.055] | `−0.0006` | `−0.0006` | — |
| 300 | 13 (400) | `+0.7764 → +0.7746` | `−0.0018` | `−0.0053` | ⟨jk⟩ |
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
(`J = 15`, twist edge, which is not its minimum), the largest unfavourable shift of a per-`J` minimum `−0.0041` (`J = 11`, 300). The two
sizeable shifts are favourable (`J = 20`: `+0.042`, `J = 18`: `+0.012`), the same pattern as HC's `J = 31`. Together with HC (1000:
`J = 5, 17, 31, 36, 40`; 300: `J = 5, 16, 24`) `H_edge` is now computed on **every row that carries a near-gap state at `Δ_g = 300` for
`J ≤ 24`** and at `Δ_g = 1000` for `J ≤ 17` and `J = 31, 36, 40`.

### 0.3 Task 3 — `J = 48` at `Δ_g = 1000` against direct values

⟨J48⟩

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
  `cache/cl_pts_*_J<J>_D400.pkl`. ⟨r2: plus the `Δ = 200` row for `J = 9, 5, 7` (`x = 0.042 … 0.1`, same `ν`).⟩
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
⟨r2⟩

## 2. Task 2: method and notes
* Points: `CL_MODE=hc HC_METHOD=n1 HC_ADAPT=1 HC_M=11 HC_MREF=11 cl_points.py J D 0.8,0.825,0.85,0.875,0.9` (= `hc_points.py` n1 path;
  even `J`: w0 2×2 EE integrand + C3H; odd `J`: 1×1). Errors `≤ 1 %` except the dropped `x = 0.9` nodes listed in §0.4.
* Chains: 300 — z4 (`cl_teven.sh … edge`: `HC_INJ=J:rows HC_IXMIN=0.7 HC_DESPIKE_INJ=0` against the same chain without it; odd `J`:
  `cl_todd.sh … edge`, small-`ν` + edge vs small-`ν` only); 1000 — the EE_NEARGAP main chain `hc_chain w2 1000` (`cl_t1000.sh`), all
  edge rows of the `J` injected. Variants: Hampel on the injected nodes (`edge_desp`; changes `≤ 0.004` except `J = 20`: `+0.885`); edge
  jackknife where a shift exceeds `0.003`.
* **Code fix needed for even `J < 11` (in the wrapper, not in hc_chain):** hc_chain injects `C3H` edge nodes only through `q3._q2_even`,
  which the w2 loader uses for `J ≥ 11`; for even `J < 11` the `C3H` rows come from `W2._q3_load` (l9 caches), so `C3H` edge nodes would
  silently not be injected. `cl_chain.py` injects them there too (`CL_C3HLOW=1`); the injection log lists both kinds for every run.
* Edge node values grow towards the edge and oscillate with the `D_h` period, as in HC: e.g. `(400, 20)` C2 `++` `+988, −1021, −1052,
  +3151, −2721` and `+−` `−11, +1685, −297, −3419, +3972` (`x = 0.8 … 0.9`; held value `+1069 / −623`); `(400, 8)` `++` `−200 … +640`.
* Representation: the protocol's natural spline through the nodes (spacing `0.025`, at least ~4 nodes per half-period at `J ≤ 22`);
  HC's shape-fit edge models exist only in ee1000 mode and were not needed here. Grade of the shifts: [NE + M].

## 3. Task 3: `J = 48` at `Δ_g = 1000`
⟨J48 details⟩

## 4. Regradings
⟨regrade⟩

## 5. What remains open
⟨open⟩

## 6. Reproduce
⟨repro⟩
