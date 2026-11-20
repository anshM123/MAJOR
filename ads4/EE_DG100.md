# `Ω^EE` at `Δ_g = 100`: per-state positivity near the gap (n), with direct small-ν crystal values and the non-local weight

**Status: complete for the briefed scope (final numbers below); n at Δ_g = 100 remains [I].**

*19–20 Nov 2026, EE-Δ_g=100 agent. Unrefereed, not committed by me (the coordinator's WIP commit e5b9ce1 holds the pre-restart state).
No existing document or project code was edited. New code (`ads4/code/spinfunc/`): `y1_rows.py` (+ `y1_queue.sh`) K₂ EE rows,
`y4_direct.py` (n1 / n3 / p8 EE points), `y5_n4.py` (n4 two-thimble EE points), `y2_chain100.py` (the chain), `y3_runs.sh`
(variants), `y_queue2.sh`, `y_queue3.sh`; caches `spinfunc/cache/y*`. Project modules, the 300/1000 sweep's `w0_ee_common` and
its `w1_*`/`u5_*` caches were imported/read only.*

Grades: [P] proved, [D] derived, [DV] derived + checked, [NE] numerical evidence, [M] model, [I] open.

## 0. Verdict

| protocol (`Ω^EE`, NONLOCAL weight, min-norm gauge) | even `n` / min `λ/N_s` (state) | odd `n` / min (state) | combined | grade |
|---|---|---|---|---|
| **main** (rows 100/200/400/…, direct on the Δ = 100 row: n4 `x = 0.05`, n3 `0.1, 0.15`, n1/p8 `0.2–0.3`) | **1 / +0.2398** (`J = 4`, 1.04Λ) | **1 / +0.2286** (`+−`, `J = 5`, 1.05Λ) | **n = 1, min +0.2286** | [NE + M] (J ≤ 10) |
| same, MC over direct-node errors (40 draws) | +0.237 ± 0.010 (lowest draw +0.193) | +0.229 ± 0.002 | n = 1 in every draw | [NE] |
| **+ odd extra rows** `m = 1.05–1.3Λ` (non-even-integer Δ, phase stripped; n4 direct at `x = 0.05`) | 1 / +0.2398 | **1 / +0.0315** (`J = 9`, 1.1Λ); MC +0.005 ± 0.044 (range −0.07 … +0.10); +0.195 with p8's `a`-fit nodes | **n = 1 at central value, marginal** | [NE + M + I] |
| + even extra rows (K₂ only, unstable `a_h`) | **n = 1.35, −0.12** (`J = 10`, 1.28Λ) | as above | n ≠ 1 in this [M] variant | [M] (weak) |
| canonical gauge (x2 `X_can`) | n = 1.35, −4.45 (`J = 10`, 1.1Λ); MC −9 ± 6 | +0.228; MC +0.03 ± 0.67 | **unusable at Δ_g = 100** | [NE] |
| local de-shift `(1 ± u_EE)/2` (comparison only) | +0.2397 | +0.2283 | n = 1 | (singular at 6.64; MC even lowest +0.07) |
| `q²` (`Ĉ′₂`) through the identical chain + its direct values (incl. new even n4 C2 at `x = 0.05`) | **n = 1.30–1.35, −0.195 (local) / −0.336 (nonlocal)** (`J = 10`, 1.21Λ) | +0.221 | n ≠ 1 | [NE + M] |

* **Main answer [NE + M + I]: `n = 1` at `Δ_g = 100` for `Ω^EE` on the computed set (even `J = 0–10`, odd `J = 5, 7, 9`), main protocol,
  min-norm gauge; binding state odd `+−`, `J = 5`, `m = 1.05Λ`, `λ/N_s = 0.2286` (even: `J = 4`, 1.04Λ, 0.2398).** Joined with the
  leading-Regge band model (253 states, `β ≤ 2.2`, min 0.2263 at 1.05Λ, `b̂ ≈ 10`, EE_FUNCTIONAL §4 [M]): combined min **0.2263** (band)
  / 0.2286 (near gap, computed corrections). The nominal bound `2.21·10⁻⁴` (`2.22·10⁻⁴` with `δk(0)`) would follow from `n = 1` on ALL states; that is not established (next bullets).
* **Extra rows (task 2) [NE + M + I]:** the robustness test that broke `q²` (odd `J = 9`, 1.2Λ, −0.55 → −1.35 with direct values) does
  **not** produce a clear failure for EE: odd `J = 9` at 1.1Λ stays at +0.03 (consistent `a_h` rule) / +0.19 (p8's rule), local +0.23,
  K₂-only +0.22; but it is **marginal**: MC over the n4 errors 0.005 ± 0.044, `t*_J` upper edge 1.03, and the extrapolation below
  `x = 0.05` moves it to +0.006 (hold) / −0.09 (zero). The K₂-only even extra rows give `J = 10` at 1.28Λ −0.12 [M, unreliable rows].
  **So `n = 1` at `Δ_g = 100` is supported by the main protocol but not robust: [I].**
* **What remains [I]:** (o) the odd `J = 9`, 1.1Λ extra-row state (needs n4 `x = 0.05` at ~4× statistics, a direct node below 0.05 or a
  finer `x` grid at `x ≤ 0.1`, and a stable `a_h` determination: a 5 % change of `a₊₋` on that row moves `λ/N_s` by 0.16); (i) near-gap states with `J ≥ 11` (exact-skeleton `J = 11–28` at `m ≥ 1.11Λ`): only the 300/1000 sweep's K₂ EE rows
  (`Δ ≥ 336`) exist, transferred to `m ≈ 110–130` by a model; **both transfer models fail at `J ≥ 18–24`** (§2.5) — a model-level
  counter-indication, not a computed failure; (ii) even extra rows have K₂ only (no direct `x = 0.05`), and their per-row `a_h` is
  unstable (§2.4); (iii) the canonical gauge cannot be evaluated with present data; (iv) everything inherited ([D]/[I] items of
  EE_FUNCTIONAL/EE_FOUNDATIONS, chain model M, beyond-leading-Regge gauge question).
* **Side result [NE + M]: the `q²` functional fails `n = 1` at `Δ_g = 100` in the S1_S4-corrected chain** (aligned off-diagonals,
  `c_off`, `Ĉ₃` without de-shift): already with K₂ only (n = 1.35, `J = 10` −0.33, `J = 8` −0.05) and with direct values (−0.20/−0.34).
  The driver is the helicity off-diagonal (`max|dF₊₋|/√NN` 0.87–1.2 at `J = 10`); for EE it is ≤ 0.05. The earlier "q² n = 1 at 100"
  (H0 §25.6b, 0.2206) used the unaligned l9 conventions.

## 1. Method

### 1.1 Chain (`y2_chain100.py`)
* S1_S4 / q3_analysis conventions: K₂ linear chain, rows linear in `1/m` (held outside), cubic spline in `x`, `c_off = mΓ_off`, row
  alignment by the sign of `q²`'s `m Im Γ_off` (computed in the same basis; EE's own imaginary mixing vanishes at leading Regge),
  `Ĉ₃` (o_H q) without de-shift (S1_S4 §1.4), phase strip `logK → logK − iπ(Δ+J)` on every record (asserted), combined even + odd pass
  over all exact-skeleton states of `Δ_g = 100` (`load_states(100)`), same base model (flat `B₂+B₃` + curvature + scalar remainder).
* **Validation [NE]:** on the `q²` rows with K₂ only it reproduces `q3_analysis.run([100])` to all printed digits (n, minima, per-J
  `t*`, alignment overlaps (4,100) −0.99, (8,400) −1.00, (10,200) −0.94 as S1_S4 §2.2(a)). The local `k = 2` term on the 1200-node fine
  quadrature equals the chain's own quadrature to `10⁻⁴`.
* EE rows: K₂ (`y1`: Δ = 100 all J, Δ = 200 for J = 10, Δ = 400 even; `w1` (other agent, read-only): odd Δ = 200/400/800/1600, even 1600);
  helicity basis per row from the `q²` records (u6 convention); `a_h` fitted on the three smallest `x` with `ν ≥ 9` on the near-gap rows
  (l9's rule would use the `ν = 4.9` node, where the local `u_EE` de-shift inside the fit is near its pole; variant `afit0` below).

### 1.2 The non-local weight (task 1)
`∫ψ₂Ĉ′_EE = Re∫X Ψ̂` exactly (Abel identity, u3), with `Ψ̂` the forward-normalised crystal value. The `k = 2` term is evaluated as
`Re∫X[G/a_h − Â_h(f_LO)] + ∫ψ₂ f_LO` (the LO part exactly by the identity, H7 §24.8 item 3), in the χ-normalisation:
`X·Ψ̂ = Y·a(ν)Ψ̂`, `X·Â_h[f] = Y(ĥ₋f(ν+i) + ĥ₊f(ν−i))`, `ĥ = R∓/n` (u3's `n`). The interpolated object is
`Γ_NL = m² a(ν)(G/a_h − Â_h[f_LO])/env`, finite through the real zeros 1.17, 6.64 (nothing divided by `h₋+h₊`); `Y` from u3 (min-norm) or
x2 (canonical), integrated on a 1200-node graded quadrature in `p`. The raw crystal ratios at the real nodes are used directly; no local
de-shift enters the nonlocal number. (The reference's own shift `h₋C_ref(ν+i)+h₊C_ref(ν−i)` vs `(h₋+h₊)C_ref` differs by `O(uν J²/m_ref²)
~ 10⁻⁶u` and is neglected.) **`Im Ψ̂`**: direct nodes keep their imaginary part; at K₂ nodes with `ν < 20` it is replaced by the exact
leading-Regge value (H7 §24.3: K₂ Im 4–7× too large). Variant `imK2` keeps K₂ Im: no change at 4 digits.
* `ψ₂` of the chain and u3's `φ` agree to `6·10⁻⁸`; `max|X − aY|/max|X| = 0` (q²), consistent for EE.

### 1.3 Direct values computed here [NE]
| node | method | EE: kinds / points | errors (Γ units, `a(ν)`-weighted for EE) |
|---|---|---|---|
| `x = 0.05` (`ν ≈ 4.9`, Δ = 100 row) | n4 two thimbles, CRN (`y5_n4.py`) | even J = 4, 6, 8, 10 (EE, C2, C3H on common samples); odd J = 5, 7, 9 | EE 56–580 c2 (J = 9: 579); C3H J = 6, 9, 10 > 1000 → K₂ kept (C3H direct = K₂ at accepted nodes, NEARGAP100 §6) |
| odd extra rows `x = 0.05` (`ν = 5.17–6.40`) | n4 (`odd:EE`) | J = 5 (4 rows), 7 (3), 9 (3) | J = 9: 282 (1.1Λ); others see §2.3 |
| `x = 0.1, 0.15` | n3 hybrid contour (`y4 n3even/n3odd`) | even 4–10, odd 5–9 | 6–110 c2 (even J = 8, `x = 0.15`, `P₈` zero: 822 → K₂) |
| `x = 0.2, 0.25, 0.3` | n1 (even) / p8 (odd) | even 4–10, odd 5–9 | 3–40 c2; `x ≥ 0.4` K₂ (accurate there, H0 §24.4) |
Acceptance: Γ-unit error ≤ 600 c2 (`Y2_GERR`; variant 300). 264 direct node values used in the main run.

## 2. Results

### 2.1 Main protocol, per J (nonlocal min-norm | local), `cache/y2_chain100_main.txt`
| J | 0 | 2 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|---|
| min `λ/N_s` NL | 0.9999 | 1.0002 | **0.2398** (1.04) | **0.2286** (1.05) | 0.3050 (1.25) | 0.7045 (1.1) | 0.4606 (1.65) | 0.6851 (1.13) | 0.3915 (1.24) |
| local | 0.9999 | 1.0004 | 0.2397 | 0.2283 | 0.3052 | 0.7027 | 0.4075 | 0.6438 | 0.5070 |
* Off-diagonal `max|dF₊₋|/√NN ≤ 0.045` (q² same chain: up to 1.2) — the leading-Regge off-diagonal of EE vanishes [P], confirmed near the gap.
* Whole-correction window `t* = (−21.1, +2.16)`; **window of the `x < 0.075` region (the n4 node), its whole correction × t:
  `n = 1` for `t ∈ [−0.23, +8.80]`** (lower edge: the correction there is load-bearing and positive). Extrapolation below `x = 0.05`:
  spline / hold / zero give min 0.2286 / 0.2287 / 0.2288. K₂-only (no direct): 0.2212.
* MC over direct-node errors (40 draws): even 0.237 ± 0.010 (range 0.193–0.242), odd 0.229 ± 0.002 — `n = 1` in all draws (nonlocal);
  the local form's even minimum spreads to 0.07 (amplification by `u_EE` near 6.64).

### 2.2 Canonical gauge [NE]
`Y_can` (x2): main protocol `n = 1.35`, min −4.45 (`J = 10`, 1.1Λ), MC even −9 ± 6, odd +0.03 ± 0.67; with odd extra rows worse.
`X_can/ψ₂ ≈ −160…−190` at `ν = 3–5` and `+66` at 10 multiplies the beyond-Regge data (n4 errors 50–600 c2, K₂ extrapolation below
`x = 0.05`) by O(100): **the canonical gauge is unusable at `Δ_g = 100` with present data** (not a statement about the functional's
positivity: errors dominate). It requires either the exact spinning magic-block identity (then the gauge is irrelevant) or ~10× more
precise small-ν data.

### 2.3 Odd extra rows (task 2) [NE + M] (`cache/y2_chain100_xrows_odd.txt` and variants `*_xrows`)
Rows `m = 1.05` (J = 5), 1.1, 1.2, 1.3Λ: K₂ EE (`y1`) on all nodes, n4 EE at `x = 0.05` (`ν = 5.17 / 5.42 / 5.91 / 6.40`; accepted where the
`a(ν)`-weighted error ≤ 600 c2: J = 9 at 1.1Λ 282 c2), p8/n3 `C3H`/K₂ as in the q² protocol; phase `logK − iπ(Δ+J)` stripped.

| variant (nonlocal, min-norm) | odd min `λ/N_s` (state) | comment |
|---|---|---|
| **default** (`a_h` from `x = 0.1, 0.15, 0.2` on every near-gap row) | **+0.0315** (`J = 9`, 1.1Λ) | `t*_J = (−∞, 1.03)`; `x < 0.075` window `[0.70, 16]` |
| MC over direct errors (40 draws) | +0.005 ± 0.044 (−0.07 … +0.10) | about half the draws within 0.03 of 0; some negative |
| `a_h` from `x = 0.15, 0.2, 0.25` (p8 rule) | +0.1946 (`J = 9`, 1.1Λ); MC 0.18 ± 0.04 | `a₊₋` differs by 5 % (1.12 vs 1.18·10⁻⁹) |
| extrapolation below `x = 0.05`: hold / zero | +0.006 / **−0.090** | the region `ν < 5.4` carries `X/ψ₂ ≈ 1–2.6` |
| K₂ Im kept / acceptance 300 c2 | +0.030 / +0.0315 | |
| K₂ only (no direct) | +0.221 (J = 5) | |
| local de-shift | +0.227 (J = 5); MC lowest −0.08 | |
| canonical gauge | −20 (`J = 9`) | unusable (§2.2) |
| `q²`, same chain, its own direct values | +0.224 (odd) — the q² odd sector does not fail here; q² fails in the even sector (§0) | |

Reading: the extra-row `x = 0.05` node at `J = 9` is again the sensitive spot, but with direct EE values and the non-local weight the
result is a margin of `≈ 0 – 0.2 N_s` (variant-dependent), not the robust negative value found for `q²` with the local de-shift.
### 2.4 Even extra rows (K₂ only) [M]
K₂ rows EE, C3H, C2 at `m = 1.05/1.1/1.2/1.3Λ` (`J = 4`), `1.1/1.2/1.3Λ` (`J = 6, 8, 10`) were computed (`y1 evenx`). Their per-row `a_h`
is unstable (fits on `x = 0.15–0.25`: `J = 10` at 1.1Λ gives `a₊₊ < 0`; with the `ν ≥ 9` rule `a₊₋(J = 8)` jumps 1.4 → 8.8 → 4.0 → 1.6·10⁻⁹
across 1.0–1.3Λ, while odd rows are smooth), i.e. K₂ is not reliable on these rows at `x ≤ 0.2` (as H0 §24.4 found for `q²` at `m ≈ 98`).
Reported as a variant only: with them (EE only, `cache/y2_chain100_xrows_all.txt`) even `J = 10` at 1.28Λ is −0.12 (nonlocal; local
+0.24), `J = 8` 0.33, `J = 6` 0.23, `J = 4` 0.24 — a weak [M] counter-indication pending direct `x ≤ 0.2` values on these rows.

### 2.5 Even `J = 12–28`, odd `J = 11, 13` near the gap [M] → [I]
Only the 300/1000 sweep's K₂ EE rows exist (Δ = 336/400/500/800/1600; `J²/m ≈ 0.4–2.3`), while the `Δ_g = 100` states sit at `m ≈ 110–130`.
* rows held below the lowest row (chain default, Γ m-independent): `J = 11–14` positive (0.61 / 0.45 / 0.43 / 0.09), `J = 16` −0.19,
  `J = 18` −1.50, `J = 20` −1.34, `J = 24` −3.30, `J = 28` −2.98 (1.2–1.3Λ); driver the DIAGONAL correction (up to 4 N_s), not the off-diagonal.
* `Γ ∝ m` below the lowest row (relative correction `∝ J²/m`): `J ≤ 16` positive (≥ 0.53), `J = 18` −0.07, `J = 20` +0.005, `J = 24` −0.62,
  `J = 28` −0.57.
Neither transfer is controlled (K₂ at `J²/m ≈ 3–7` is outside its range; no row at the states' masses). **`n` for `J ≥ 16` near the gap at
`Δ_g = 100` is [I], with a model-level counter-indication at `J ≥ 18–24`** — the `J²/m` diagonal growth EE_FUNCTIONAL §5 flagged.

## 3. Graded summary
| statement | grade |
|---|---|
| chain reproduces q3_analysis at Δ_g = 100 (q² rows, K₂) exactly | [NE] (code check) |
| `q²` (`Ĉ′₂`) in the S1_S4-corrected chain at Δ_g = 100: n ≠ 1 (J = 10, 1.21Λ; K₂ and direct) | [NE + M] |
| nonlocal evaluation in the χ-normalisation, finite at the real zeros; LO part exact | [D] + [DV] (quadrature 10⁻⁴) |
| **EE, main protocol, min-norm: n = 1, min 0.2286 (odd J = 5, 1.05Λ); even 0.2398 (J = 4, 1.04Λ); J = 0–10 + odd 5, 7, 9** | **[NE + M]** |
| robustness of the main protocol: MC over direct errors, extrapolation, Im rule, `a_h` rule, K₂-only, local vs nonlocal: n = 1 in all | [NE + M] |
| odd extra rows (n4 direct at x = 0.05): n = 1 at central values, J = 9 at 1.1Λ marginal (+0.03 / +0.19 by `a_h` rule; MC ±0.04–0.05; −0.09 with zero extrapolation) | [NE + M]; sign [I] |
| even extra rows (K₂ only): J = 10 at 1.28Λ −0.12 | [M] (rows unreliable) |
| x = 0.05 node: n = 1 for t ∈ [−0.23, 8.8] × its correction | [NE + M] |
| canonical gauge at Δ_g = 100 | unusable with present data [NE]; n there [I] |
| joined with the leading-Regge band (min 0.2263): combined n = 1 for the computed set | [M] (band) + [NE + M] |
| near gap J ≥ 16 (even), J ≥ 11 beyond K₂ transfer | [I]; [M] failure at J ≥ 18–24 in both transfer models |
| bound 2.21·10⁻⁴ (2.22·10⁻⁴ with δk(0)) at Δ_g = 100 | conditional on n = 1 over ALL states → **not established** while §2.5 is open |

## 4. Hand-over
* **Data (all complete, survive in cache/):** `y1_pts_{EE,C2,C3H}_J*_k2.pkl` (K₂ rows), `y1_p3_EE_J{5,7,9}.pkl`, `y4_{n1,n3,n3_odd,p8}_pts_EE_J*.pkl`,
  `y5_pts_{EE,C2,C3H}_J{4,6,8,10}.pkl`, `y5_odd_pts_EE_J{5,7,9}.pkl`, `y5_raw.pkl` (raw thimble records, all 17 n4 jobs).
* **Outputs:** `cache/y2_chain100_<tag>.txt` (+ `y2_recs_<tag>.pkl`): `main`, `can`, `ex_hold`, `ex_zero`, `imK2`, `k2only`, `jhi`,
  `jhi_scale`, `val_q2_k2` (validation), `q2_dir`; after the container restart (all production data had survived): `xrows_odd` (odd extra rows),
  `xrows_all` (EE, odd + even extra rows), `can_xrows`, `ex_{hold,zero}_xrows`, `imK2_xrows`, `afitp8_xrows`, `afit0`, `gerr300_xrows`,
  `k2only_xrows` (`y3_runs.sh`, log `cache/y3_runs2.log`). `y2_chain100_xrows.txt` is a pre-restart file with partial even extra rows and the
  old `a_h` rule (superseded; its EE odd numbers equal `afitp8_xrows`).
* **Rerun:** `OMP_NUM_THREADS=1 Y2_MC=40 Y2_TOL=1 python3 y2_chain100.py main` (≈ 10 min); env switches documented in the module docstring
  (`Y2_GAUGE`, `Y2_XROWS`, `Y2_XEVEN`, `Y2_EXTRAP`, `Y2_JHI`, `Y2_BELOW`, `Y2_DIRECT`, `Y2_NUIM`, `Y2_AFITNU`, `Y2_AFITXR`, `Y2_GERR`, `Y2_MC`, `Y2_TOL`).
* **Next steps (value order):** (1) rows at the near-gap masses for `J = 12–24` at Δ_g = 100 (`Δ ≈ 120–140`, direct n1/n3 at `x ≥ 0.1`,
  `J²/m ≈ 1–5` excludes K₂) — decides §2.5; (2) even extra rows with direct `x ≤ 0.2` (n3/n4) to fix their `a_h`; (3) J = 9 Δ = 100 n4 at
  4× statistics (its error 579 c2 is the largest accepted); (4) canonical gauge only after the exact magic-block identity or ~10× data.
