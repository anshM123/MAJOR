# Closing `H_OO` at `J ≥ 12` and the `ν < 8.3` window at `Δ_g = 300` (continuation of HYPOTHESES_CLOSE.md)

*21 Nov 2026, hypothesis-closing agent (2). Unrefereed, not committed. No existing document or code file was edited. New code
`ads4/code/spinfunc/oo_*.py` (+ `oo_queue*.sh`, `oo_t.sh`, `oo_window.sh`). Caches are under `spinfunc/cache/`:
* `oo_*` files, one file per state. The reference points of a state live only in that state's file, and no reference key is shared
  (asserted when files are merged).
* `hc_chain_z4_oowin_*`, written by the unchanged `hc_chain.py`.
* the tagged thimble caches `z2oo*_pts_*`.

Fixed and not tuned on any sign found here:
* the functional `Ω^EE` (min-norm gauge, NF02);
* the z4 protocol of EE_DG300_SMALLNU §6 and HYPOTHESES_CLOSE §4. It consists of direct n3 nodes, the non-local min-norm `X` below
  `ν = 25`, `SM_CURV=row`, the Jacobi fix, the estimator `best`, the value held below the lowest node, and the Hampel despike.*

Grades: **[P]** proved, **[D]** derived, **[NE + M]** direct crystal computation in the chain model, **[NE]** numerical evidence,
**[M]** model/extrapolation, **[I]** open.

## 0. Summary

| item | result | grade |
|---|---|---|
| **OO block, `J = 12–22`, `Δ_g = 300`, direct small-`ν` nodes** (3×3 integrand) | min `λ_OO/N_s` over near-gap states and twist-edge masses: **`+0.886 ± 0.002 / +0.998 ± 0.004 / +1.007 ± 0.009 / +0.994 ± 0.032 / +0.936 ± 0.059 / +0.842 ± 0.085`** at `J = 12 / 14 / 16 / 18 / 20 / 22` (jackknife). `J = 18` was `+0.222` in the superseded treatment; that number is reproduced exactly, and §1.4 shows it was an artefact. | **[NE + M]** for `m ≤ 1.3Λ` |
| **New structural fact** | the forward-normalised OO correction **equals the `++` correction node by node**. `K₂` rows (`J = 12, 14, 18, 22`, 16–19 nodes): `|Γ_OO − Γ_++| ≤ 0.05 c₂` (C2), `≤ 0.12 c₂` (C3H). Direct small-`ν` records at the same estimator (`J = 18–22`): equal within `≤ 1 c₂` at `ν ≥ 10` and within `1–12 c₂` at the noisy `ν = 8.3` nodes (RQMC error there `≈ 10³ c₂`). Hence `λ_OO/N_s` is the EE block's normalised `++` diagonal, which is `≥` the EE 2×2 minimum. | identity **[NE]**; the implication "EE-block positivity ⇒ OO positivity" is then **[D]** |
| Validation `J = 4–10` | `+0.2301 / +0.3071 / +0.4774 / +0.6022` reproduced exactly; with direct OO spot nodes at `J = 4`: `+0.2306` | [NE + M] |
| **`ν < 8.3` window, `(336, 28)` / `(336, 32)`** | Not decided. `J = 28` improved: a pooled thimble node at `ν = 5.0` (two new/independent seeds) passes the scalar check (`≤ 1.7σ`) with `Re Γ_c = −1300 ± 980`, and failure now needs a uniform `Γ_c ≲ −7.3·10³` on `ν < 5`; gap state `+0.93 ± 0.24`. `J = 32`: no validated node. | `J = 28`: [NE + M] to `ν = 5`, [M]/[I] below; `J = 32`: **[I]** |
| … why | 98 % of the window's weight lies at `ν < 6` and 58 % at `ν = 3–5`. Only Lefschetz thimbles reach there, and they fail the project's own scalar check: the new `J = 32` thimble node misses the exact scalar by `Δc₂ = +32 600 ± 4900` (`6.7σ` in `Re`, `19σ` in `Im`). On that node the two estimators give opposite signs (`best`: twist-edge state `−1.63 ± 0.24`; `dir`: `+1.53`). An independent second seed gives `+0.884`, and the pooled seeds `+0.939`, but both still fail the `Im` check (`3.4σ`, `11σ`): a systematic thimble phase bias. The `J = 28` seeds disagree by `1.1·10⁴ c₂`, but only because of the old `2⁹` main seed; the two `2¹⁰` seeds agree. The n3 contour below the fold has 70–120 % errors at `ν = 5`. | |
| … best achievable bound | `J = 28`: failure needs a uniform `Γ_c ≲ −1.0·10⁴` (main thimble seed) or `≲ −7.7·10³ c₂` (repeat seed) on `ν < 5`, and `≲ −5400` on `ν < 8.3` without thimbles. `J = 32`: `≲ −4250` on `ν < 8.3` without thimbles. With the two new (unvalidated) seeds pooled, the central value passes (min `+0.865 ± 0.12`), but failure needs only `≲ −1.4·10³` on `ν < 5`, so **`J = 32` is the most fragile state at 300**. | |
| **Item 3 (EE block), even `J = 12, 14, 16`, direct small-`ν` nodes** (the even 2×2 block of the same 3×3 points) | min `λ/N_s` near the gap: **`+0.673 ± 0.006 / +0.750 ± 0.007 / +0.884 ± 0.022`**, against convention A `+0.630 / +0.653 / +0.707`. Convention A was again pessimistic. | **[NE + M]** for `J = 12, 14, 16`. `J = 6–10` and odd `J = 7–15` were not computed. |

The binding state at `Δ_g = 300` is unchanged: `J = 5`, `+0.2038` (hybrid, twist edge). Nothing here goes negative except on a node
that fails validation (§3.3).

## 1. Item 1: the OO block at `J ≥ 12` with the validated small-`ν` treatment

### 1.1 Method
* **Direct points (`oo_points.py`).** hc_points' n3 path is unchanged:
  * n3 hybrid contour, `ν_f = 13` below `ν = 13`;
  * proposal adapted on the state and re-used for the reference `Δ_ref = 2·10⁴`, with the same Sobol seeds;
  * `Q_J` overflow rescaling, `Q_J` recurrence, Fourier floor `patch_S`;
  * `2¹⁰·2·8` points.

  `hc_points.setup(odd=True)` is **forced at even `J`**. The graviton evaluators are `n3_run_odd.OddEval` on p8's point-wise copy of
  `p3_odd_action.build`:
  * kind `EE` = `w0.make_build_odd_EE`;
  * kind `C3H` = `o_H q`;
  * the integrand is 3×3 [even(2), odd(1)], with OO = entry `[−1, −1]`.

  The scalar is the even-`J` scalar on the same samples, so the scalar control variate applies.

  Nodes:
  * `ν ≈ 8.35, 10.0, 11.7, 13.3, 16.7, 20.9/19.9, 25.0/29.9, 33.4/39.8` on `(336, 18)`, `(400, 18)`, `(400, 20)`, `(336, 22)`;
  * `ν ≈ 8.36 … 19.9` (6 nodes) on `(400, 12)` and `(400, 14)`, and 8 nodes on `(400, 16)`;
  * spot nodes `ν = 10, 20, 40` on `(400, 4)`.

  Cost: ≈ 10 min per node (state + reference).
* **`K₂` OO rows (`oo_rows.py`).** These are hc_oo_rows with one file per state row and the EE-row `x` grid:
  * 19 nodes `x = 0.05 … 0.8` at `J = 16, 18, 20, 22`. At `J = 18` this adds the 8 missing nodes to the old 11-node rows.
  * 16 nodes at `J = 12, 14`.
* **Records (`oo_chain.oo_direct`).** `M0 = E b_odd⁻ᵀ`, where `b_odd` is the `p2_shadow_odd` matrix whose stored residual equals the
  `shadow_resid` of the row's `K₂` OO records. `M2 = M4 = 0` and `logC = logN`. A direct node replaces a `K₂` node at the same `ν`.
  The estimators are `dir`, `cv` and `best`; `best` picks, per node, whichever gives the smaller jackknife error of the OO element.
  `OO_JK = i` leaves replica `i` out.
* **Chain (`oo_chain.py`).**
  * hc_oo's construction (OO rows in place of the even rows; `++` keys = `Γ_OO`; `+−` and off-diagonal keys = 0) is moved into the
    **z4 chain**, set up exactly as `hc_chain.py z4`.
  * For every OO row with direct nodes, the complex `Γ_c,OO` (curvature-subtracted) and `Γ_NL = a_χ Γ_c` are put into
    `z4_chain.NLROWS[(J, D)]['++']`. The chain's `++` density below `ν = 25` is therefore the non-local min-norm
    `Re[(Y/φ)Γ_NL,OO]` (z4's term, extended to the OO entry), and the local values at direct nodes `ν < 14` are flattened as z4 does.
  * The Hampel filter is applied to the OO rows (variant: none).
  * A single-row `J`'s duplicate row carries the same NL row (`OO_DUPNL=1`). This matters only at `J = 22`, for `m > 1.118Λ`:
    `+0.871` vs `+0.912` at 1.15Λ.
  * In the same process, run0 is the unmodified z4 chain (EE block), with the EE direct nodes taken from `hc_pts` (`J = 18–22`) or
    from the even 2×2 block of `oo_pts` (`OO_EESRC=oo`, used for `J = 12–16`).
  * Per state, `λ_OO/N_s = (base_++ + dF_OO)/N_s,++`.

### 1.2 Validation [NE]
* **`J = 4–10` (task requirement).** `oo_chain` with `K₂` OO rows only (`OO_DESPIKE=0`) reproduces hc_oo **exactly**:
  * `λ_OO/N_s = +0.2301 / +0.3071 / +0.4774 / +0.6022` at `J = 4 / 6 / 8 / 10`;
  * EE-block minima `+0.2204 / +0.3073 / +0.4293 / +0.6026`;
  * with the Hampel filter on: identical (`cache/oo_chain_val410*.txt`).

  With the direct OO spot nodes at `J = 4` (best): **`+0.2306`**. The nodes are `Γ_c = +2 / +25 / +21 c₂` at `ν = 10 / 20 / 40`, and
  direct/`K₂` = `0.996 / 0.999` at `ν = 20 / 40`.
* **Superseded `J = 18` number reproduced.** Convention A (`OO_DIRECT=0 Z4_DIRECT=0 Z4_NL=0`) on the old 11-node rows gives
  **`+0.2224`** at 1.0586Λ, with `dF_OO/N = −0.778` vs EE `++` `−0.495` and EE block `+0.4718`. These are hc_oo's numbers.
* **Regression of the EE block.** run0 reproduces the hc z4 values: `J = 18 / 20 / 22` near-gap `+0.910 / +0.863 / +0.841`, all-`m`
  `+0.872`. Taking the EE nodes from the even 2×2 block of the 3×3 points (`OO_EESRC=oo`) instead of hc's w0 even integrand changes
  every state by `≤ 5·10⁻⁴` (`J = 18–22`). The even-block matrix elements agree with hc's to 4–5 digits on the same seeds.
* **Per-node checks (`oo_check.py`, `cache/oo_check.txt`), every direct node at `J = 4, 12–22`:**
  * scalar `d_D − d_ref` (exact Mellin object) within `≤ 1.2σ` (`Δc₂`: `±10³` at `ν = 8.3`, `±50–190` above);
  * OO-element ratio errors: `0.6–2.4 %` (`ν = 8.3`), `0.2–0.4 %` (`ν = 10`), `≤ 0.1–0.5 %` above;
  * direct/`K₂` of the full state/reference OO ratio: `1.0004–1.0016` at `x = 0.1` (normalisation and shadow basis consistent),
    `1.008–1.015` at `x = 0.05`. A 1 % ratio error is `≈ 10³ c₂`, which is the known small-`x` `K₂` error, as for EE `++`;
  * order-0 coupling `C3H/C2·(2 − y)` → `1.975–1.999` at `ν ≈ 33–40`. At `ν ≤ 20` it is `1.6–1.9`, which is beyond-LO, the same
    for `++`.

### 1.3 Result: the OO block at `J = 12–22`
Main run `cache/oo_chain_f_*` (`oo_t.sh f 12,14,16,18,20,22`, `OO_EESRC=oo`). States: every exact-skeleton and hybrid state, plus the
twist-edge masses `1.0395 / 1.0459 / 1.0523 / 1.0586 / 1.0649 / 1.0712Λ` (`J = 12 … 22`). Jackknife over the 8 common-seed replicas of
the direct nodes.

| `J` (rows) | **min `λ_OO/N_s`, `m ≤ 1.3Λ`** [`m/Λ`] | `best / dir / cv / no Hampel` | held / zero / linear below 8.3 | convention A | EE `++` diag at that state |
|---|---|---|---|---|---|
| 12 (400) | **`+0.886 ± 0.002`** [1.3] | `+0.886 / +0.885 / +0.886 / +0.886` | `+0.886 / +0.888 / +0.884` | `+0.851` | `+0.886` |
| 14 (400) | **`+0.998 ± 0.004`** [1.3] | `+0.998 / +0.995 / +0.998 / +0.998` | `+0.998 / +1.001 / +0.992` | `+0.860` | `+0.998` |
| 16 (400) | **`+1.007 ± 0.009`** [1.3] | `+1.007 / +0.994 / +1.010 / +1.007` | `+1.007 / +1.013 / +0.986` | `+0.766` | `+1.007` |
| 18 (336, 400) | **`+0.994 ± 0.032`** [1.09] | `+0.994 / +0.972 / +0.994 / +0.993` | `+0.994 / +1.007 / +0.956` | `+0.505` (11-node rows: `+0.222`) | `+0.994` |
| 20 (400) | **`+0.936 ± 0.059`** [1.0649, twist edge] | `+0.936 / +0.920 / +0.937 / +0.936` | `+0.936 / +0.979 / +0.870` | `+0.483` | `+0.928` |
| 22 (336) | **`+0.842 ± 0.085`** [1.0712, twist edge] | `+0.842 / +0.802 / +0.842 / +0.842` | `+0.842 / +0.897 / +0.745` | `+0.070` | `+0.841` |

* **Stress test: uniform `Γ_c` below the lowest node (`ν ≈ 8.35`).** For `Γ_c = −1000 / −2000 / −4000 / −8000`:
  * `J = 22`: `+0.81 / +0.72 / +0.54 / +0.19`; failure at `≈ −1.0·10⁴ c₂`;
  * `J = 20`: `+0.91 / +0.84 / +0.70 / +0.43`; failure at `≈ −1.4·10⁴`;
  * `J = 18`: `+0.97 / +0.92 / +0.81 / +0.60`; failure at `≈ −2·10⁴`;
  * `J = 16`: `+0.99 / +0.95 / +0.88 / +0.73`;
  * `J = 14`: `+0.99 / +0.97 / +0.93 / +0.82`;
  * `J = 12`: `+0.88 / +0.88 / +0.86 / +0.84`.

  These are the EE block's numbers, by §1.4.
* **Beyond `1.3Λ`.** The OO rows exist only on the lowest row (400 or 336), so states above `m_row` are held at that row.
  * All-`m` minima: `+0.700 / +0.777 / +0.813 / +0.874 / +0.867 / +0.842` at `2.5 / 3 / 3 / 4 / 4 / 1.07 Λ`.
  * These are the LO bases (`dF_OO ≈ 0` there) and match the EE block's all-`m` minima: `+0.673 / +0.750 / +0.810 / +0.872 / +0.863 /
    +0.841`.
  * Graded **[M]**: held row. By §1.4 the EE rows at 800/1600 would give the same OO values.
* **Direct OO node values** (`Γ_c`, `c₂`, best; real part):
  * `J = 12`: `−504, −85, −21, +8, +27, +23` (`ν = 8.4 … 19.9`);
  * `J = 14`: `−703, −120, −37, 0, +23, +17`;
  * `J = 16`: `−837, −78, +19, +58, +73, +51, −2, −45`;
  * `J = 18`: `−824, +2, +99, +108, +86, +15, −35, −108` (336);
  * `J = 20`: `−1229, −182, −26, +15, +39, +14, −33, −46`;
  * `J = 22`: `−1227, −58, +92, +70, +28, −44, −118, −109`.

  These are identical to hc's EE `++` values (§1.4). They are `O(10²)` above `ν = 10`. At `J = 18–22` they replace a `K₂` value at
  `x = 0.05` of `−650 … −1420 c₂`.

### 1.4 Why "OO tracks `++`" failed at `J = 18`, and the identity `Γ_OO ≡ Γ_++`
* **`K₂` rows.** OO (entry `[−1, −1]` of the 3×3 p3 action) and EE `++` (w1/u5 rows) compared node by node:
  * rows `(336, 18)`, `(400, 18)`, `(336, 22)` (19 nodes) and `(400, 12)`, `(400, 14)` (16 nodes);
  * `|Γ_OO − Γ_++| ≤ 0.05 c₂` in C2 (printed `0.0` at every node) and `≤ 0.12 c₂` in C3H, against `|Γ| ≤ 1460 c₂`;
  * only the per-row normalisation differs: `a_O/a_++ = 2.20 / 2.09 / 1.74 / 1.77 / 1.14`. This is CHL's `N_same`, constant along a
    row: PARITY_ODD_ROUTEA §2 ("exactly proportional … at every Laplace order") and §12.
* **Direct points.** With the same estimator (`dir`), the OO `Γ_c` equals the EE `++` `Γ_c` from hc's independent w0 integrand on the
  same seeds, at every direct node of `J = 18, 20, 22`: within `≤ 1 c₂` at `ν ≥ 10`, within `1–12 c₂` in `Re` and `≤ 12 c₂` in `Im` at `ν ≈ 8.35`. Example, `(336, 18)`, `ν = 8.35 … 33.4`: OO `−478−1243i, −113, +39, +57, +53,
  −29, −66, −113` vs `++` `−477−1255i, −113, +39, +57, +53, −29, −66, −113` (`cache/oo_chain_m3_dircmp.txt`).
* **What failed before.** The OO rows of HYPOTHESES_CLOSE had 11 nodes at `J = 8, 10, 18`; the EE rows had 19 (`J ≥ 16`) or 9
  (`J ≤ 10`, l9/u5 caches). The natural spline of an oscillating `Γ` aliases differently on different node sets. That alone produced
  the "1.6×" at `J = 18`: on identical 24-node rows, convention A gives `dF_OO/N = dF_++/N = −0.495` and `λ_OO = +0.505`, against
  `−0.778` / `+0.222` on the 11-node OO rows. It also produced the `1.0–1.25` ratios at `J ≤ 10`: on the common nodes the `J = 4 / 10`
  rows differ by `≤ 27 / 48 c₂`, an `a`-fit offset `∝ D_h` that is equal in C2 and C3H.
* **Consequence.** Wherever the OO and EE rows carry the same nodes and estimator, `dF_OO = dF_++`. So `λ_OO/N_s` equals the EE
  block's normalised `++` diagonal `(nFn)₀₀`, which is `≥ λ_min(nFn)`, the EE-block value (Rayleigh quotient).
  * In the runs above the residual differences `|λ_OO − (nFn)₀₀| ≤ 0.03` come from estimator choices per element and from the Hampel
    filter (w2's acts on all keys, the OO filter only on `++`). Examples: `J = 14` near the gap, `+1.025` vs `+0.999`; `J = 20`,
    `+0.936` vs `+0.928`.
  * Every OO state is at or above the EE block's minimum at the same state, except where these representation differences act. The
    worst such case is `J = 22` at 1.15Λ with `DUPNL=1`: `+0.871` vs EE `+0.937`, still positive.
* **Grade.** The identity is **[NE]**: `K₂`, 5 rows, `J = 12–22`, to `≤ 10⁻⁴` relative, and all-orders direct nodes at `ν ≥ 8.3`, to
  the RQMC noise shared on common seeds. It is not derived. A proof would be a statement that CHL's `X^O_same` and the even `++`
  combination have proportional helicity amplitudes under both `Ĉ'_EE` and `Ĉ₃`; PARITY_ODD_ROUTEA §2 already conjectured this
  ("`N` is a normalisation of the whole block") [I]. **Given the identity, H_OO is not an independent hypothesis: OO positivity follows
  from EE-block positivity state by state [D].**

### 1.5 Grade of H_OO now
* `J = 4–10` (both `Δ_g`, from HYPOTHESES_CLOSE) and **`J = 12–22` at `Δ_g = 300`, near gap (`m ≤ 1.3Λ`)**: **[NE + M]** with the
  validated small-`ν` treatment. Margins: `+0.84 … +1.01`; jackknife `±0.002 … ±0.085`. Failure would need a uniform `Γ_c ≲ −10⁴ c₂`
  below `ν = 8.3` at `J = 22`, more at lower `J`.
* Every other even-`J` state (`Δ_g = 1000`, `J ≥ 24` at 300, `m > m_row`): **[D given the identity [NE]]** — the OO block is the EE
  block's `++` diagonal and inherits the EE block's grade.
* The `J = 18` binding OO value `+0.222` of HYPOTHESES_CLOSE is withdrawn. It was an 11-node aliasing artefact plus convention A.

## 2. Item 3 (partial): the EE block at even `J = 12, 14, 16` with direct small-`ν` nodes
The even 2×2 block of the 3×3 oo_points is the EE integrand on the same samples, verified to `≤ 5·10⁻⁴ N_s` against hc's EE nodes at
`J = 18–22`. So the OO points give direct small-`ν` EE nodes on `(400, 12)`, `(400, 14)` and `(400, 16)`. They come from run0 of the
same main run (z4, `best`, held below 8.36, jackknife).

| `J` | convention A (HYPOTHESES_CLOSE sweep) | **direct small-`ν`** near-gap min [`m/Λ`] | all `m` (rows 400/800/1600) |
|---|---|---|---|
| 12 | `+0.628` (sweep) / `+0.630` (this run) | **`+0.673 ± 0.006`** [1.0586] | `+0.673` |
| 14 | `+0.653` / `+0.653` | **`+0.750 ± 0.007`** [1.0459 = twist edge] | `+0.750` |
| 16 | `+0.707` / `+0.707` | **`+0.884 ± 0.022`** [1.10] | `+0.810` [3Λ, base] |

Convention A was pessimistic by `0.04 / 0.10 / 0.18 N_s`. Grade **[NE + M]** for `J = 12, 14, 16` (EE block). **Not done:**
* even `J = 6, 8, 10`, which stay on convention A and are positive (min `+0.307` at `J = 6`);
* odd `J = 7–15`, also convention A (min `+0.630`, `J = 15`), with no direct nodes. Odd `J` has no OO block, and its `+−` structure is
  the 1×1 p8 integrand (hc_points odd path).

## 3. Item 2: the `ν < 8.3` window at `(336, 28)` and `(336, 32)`

### 3.1 Where the window's weight sits
The weight of a uniform `Γ_c` in `dF` is `∝ ∫ Re X_min(ν) env dν`. From `cache/u3_chi_EE_300.pkl`:

| `ν` interval | `[0, 2)` | `[2, 3)` | `[3, 4)` | `[4, 5)` | `[5, 6)` | `[6, 7)` | `[7, 8.3)` |
|---|---|---|---|---|---|---|---|
| `∫ Re X dν` (total 18.0) | 1.50 | 3.47 | 5.53 | 4.82 | 2.32 | 0.24 | 0.12 |

98 % of the window lies below `ν = 6`, 85 % below 5, and 58 % in `[3, 5]`. `Re X_min ≈ 0` on 6–8.3. So the `ν = 6.67` and `8.3` nodes
are irrelevant to the window. The decision needs `Γ_c` on `ν ≈ 2–5`, i.e. deep below the fold (`ν* ≈ 11.3`).

### 3.2 Methods tried
* **n3 contour below the fold** (`ν_f = 13`, `J = 28`, `ν = 5.0`, `2⁹·2·8` points). Relative errors are 66–120 % (scalar 68 % /
  120 %), so the method is useless there, as expected from its growth from 2–5 % at `ν = 8.3`.
* **n4 thimbles (A₊ + A₋) at `ν = 5.0`.** About 1–2 core-hours per node at `2¹⁰·2·8` points. 2–6 % of the points are masked (walls).
  * The scalar check against the exact Mellin object (z3) is the protocol's validation:

    | node | seed | `d − d_ref` (Re, Im) | `Δc₂` |
    |---|---|---|---|
    | `J = 28`, `ν = 5.0` | main (`2⁹`) | `+0.23 ± 0.23`, `−0.19 ± 0.14` | `+20 900 ± 20 800` |
    | `J = 28`, `ν = 5.0` | s2 (`2¹⁰`) | `−0.049 ± 0.019`, `−0.009 ± 0.077` | `−4460 ± 1720` (2.6σ) |
    | **`J = 32`, `ν = 4.99` (new, tag `oo`)** | 4242 (`2¹⁰`) | **`+0.381 ± 0.057`, `−0.476 ± 0.025`** | **`+32 600 ± 4900` (6.7σ / 19σ)** |
    | `J = 32`, `ν = 8.32` (n3) | — | `−0.038 ± 0.021` | `−1830 ± 1000` |

  * The new `J = 32` thimble node **fails** the scalar validation decisively. On it:
    * the `best` estimator (`cv` for EE, chosen per the protocol rule) gives `Γ_c,++ = −18 700 ± 2000 c₂`;
    * the plain `dir` estimator gives `+4950 − 12 300i`;
    * the `cv` factor is the large complex scalar mismatch `e^{−d}`.
  * **Report of a negative state.** With `best`, the thimble node held below `ν = 5.0` (`cache/hc_chain_z4_oowin_32_x_oo_n0_*`):
    * `J = 32` at its twist edge 1.102Λ is **`−1.63 ± 0.24`** (jackknife), at 1.15Λ `−1.12`, and at 1.2Λ `−0.65`;
    * with `Γ_c := 0` below 5.0 it is `+0.22`; with `dir` it is `+1.53`.

    Because the node fails the scalar test at 6.7σ / 19σ, and the estimators disagree in sign, this is **not a measurement of the
    window** but a failure of the method. It is reported, not used. The `J = 28` seeds tell the same story: Γ_c(5.0) is `+10 070`
    (main) vs `−1110` (s2), and `NL` density `+12 400` vs `−2000`.
  * **An independent second `J = 32` seed** (6161, tag `oo2`, `2¹⁰·2·8`) was run to test reproducibility.
    * Scalar check: `Re` `+0.061 ± 0.119`, which is fine, but `Im` `−0.345 ± 0.101`, which is `3.4σ`. The phase error has the same
      sign and size as seed 1 (`−0.48`), so it is a systematic thimble phase bias, like the `0.1–0.3` rad of EE_DG300 §4.
    * `Γ_c,++(4.99) = −1066 − 9254i` (best = `dir`), against `−18 713` for seed 1: a seed-to-seed spread of `1.8·10⁴ c₂` in `Re`.
    * With this seed held below 5.0, `J = 32` is positive: twist edge `+0.884`, min `+0.825 ± 0.30` (jk) at 1.15Λ.
* **Pooled seeds (`oo_merge.py`).** Independent seeds are pooled replica-wise with point weights (`Z3_PREFER=x,oom`).
  * `J = 32`, seeds `oo + oo2`: scalar check still fails (`Re +0.225 ± 0.072`, `3.1σ`; `Im −0.446 ± 0.041`, `11σ`).
  * `Γ_c,++(4.99) = +1972 − 11 054i`.
  * Held below 5.0: twist edge `+0.939`, min `+0.865 ± 0.123` (jk) at 1.2Λ. `Γ_c := 0 / −4000` on `ν < 5` gives `+0.258 / −0.474` at
    the twist edge, i.e. failure at a uniform `Γ_c ≈ −1.4·10³ c₂` on `ν < 5` given this (unvalidated) node.
  * The spline between the `ν = 5.0` node and 8.3 is itself decisive: with `Γ_c := 0` below 5.0 every seed gives only
    `+0.22 … +0.47` at the twist edge, against `+0.885` with no thimble node and `Γ_c := 0` below 8.3.
  * **`J = 28`, a third seed** (`oo3`, 6161, `2¹⁰·2·8`):
    * `Γ_c,++(5.0) = −866 − 4380i`, against s2 `−1108 + 3610i` (`Re` consistent) and main `+10 073` (`2⁹`, the outlier);
    * scalar check `Re −0.049 ± 0.071`, `Im −0.196 ± 0.090` (`2.2σ`).
  * **Pooled s2 + oo3** (`Z3_PREFER=x,oom2`; the two `2¹⁰` seeds):
    * scalar check **passes at `≤ 1.7σ`**: `Re −0.060 ± 0.035`, `Im −0.106 ± 0.069`, `Δc₂ = −5490 ± 3160`. The scalar test therefore
      validates the method only to `≈ ±5·10³ c₂`;
    * `Γ_c,++(5.0) = −1300 ± 980 − (590 ± 1590)i` (jackknife);
    * with this node held below 5.0: `J = 28` min `+0.799 ± 0.011` [1.5Λ], gap state 1.091Λ `+0.927 ± 0.235`;
    * uniform `Γ_c` on `ν < 5` of `0 / −2000 / −5000 / −10 000`: gap state `+0.96 / +0.70 / +0.31 / −0.35`, i.e. failure at
      `≈ −7.3·10³ c₂`.
  * Pooling all three `J = 28` seeds (main weight 0.2) gives `Γ_c = +2650 − 480i`, gap `+1.29`; the scalar `Im` is off by `2.5σ`.
    The main seed's heavy tail (OO-element error 62 %) degrades the pool, so s2 + oo3 is the reference value.

### 3.3 Best achievable statement
All runs use direct bulk (`HC_INJ`, `x = 0.125–0.8`), `best`, and rows 336/500/1600 (`J = 28`) or 336/1040 (`J = 32`).

| | thimble data used | failure threshold (uniform `Γ_c` below the lowest node) | min `λ/N_s` [state]; gap/twist-edge state (held below the lowest node) |
|---|---|---|---|
| `J = 28` | none (`ν ≥ 8.3`, HYPOTHESES_CLOSE §5) | `≈ −5400 c₂` on `ν < 8.3` | `+0.797` [1.5Λ] (`Γ_c := 0` below 8.3) |
| `J = 28` | main seed at `ν = 5.0, 6.67` | `≈ −1.04·10⁴` on `ν < 5` | `+0.835` [1.5Λ]; gap 1.091Λ `+1.728` |
| `J = 28` | s2 seed | `≈ −7.7·10³` on `ν < 5` | `+0.794 ± 0.007` [1.5Λ]; gap `+0.847 ± 0.163` (jk) |
| `J = 28` | oo3 seed | `≈ −8.1·10³` on `ν < 5` | `+0.810` [1.5Λ]; gap `+1.176` |
| **`J = 28`** | **s2 + oo3 pooled (scalar `≤ 1.7σ`)** | **`≈ −7.3·10³` on `ν < 5`** | **`+0.799 ± 0.011` [1.5Λ]; gap `+0.927 ± 0.235`** |
| `J = 32` | none | `≈ −4250` on `ν < 8.3` | `+0.885` [1.2Λ] (`Γ_c := 0` below 8.3) |
| `J = 32` | seed 1 (`oo`; scalar off `6.7σ / 19σ`) | — | `best`: **`−1.63 ± 0.24`** [1.102Λ]; `dir`: `+0.962` [4Λ], twist edge `+1.53` |
| `J = 32` | seed 2 (`oo2`; scalar `Im` off `3.4σ`) | `≈ −2.6·10³` on `ν < 5` | `+0.825 ± 0.30` [1.15Λ]; twist edge `+0.884` |
| `J = 32` | seeds 1 + 2 pooled (`Im` off `11σ`) | `≈ −1.4·10³` on `ν < 5` | `+0.865 ± 0.123` [1.2Λ]; twist edge `+0.939` |

* **`J = 32` is the most fragile state at 300.** The `ν = 8.3` node (validated, `Δc₂` within `1.8σ`) is `Γ_c ≈ −1700 c₂`. The pooled
  thimble central value passes (`+0.865 ± 0.12`), but on `ν < 5` it tolerates only a uniform `Γ_c ≳ −1.4·10³ c₂` (`−2.6·10³` with
  seed 2 alone).
* **`J = 28`: improved but not closed.**
  * The `ν = 5.0` node is now resolved to `±10³ c₂` (two independent `2¹⁰` seeds that agree in `Re`) and passes the scalar check at
    `≤ 1.7σ`. Its value, `Re Γ_c = −1300 ± 980`, is of the size of the validated `ν = 8.3` node (`−1175`).
  * Failure needs a uniform `Γ_c ≲ −7.3·10³ c₂` on `ν < 5` (85 % of the window's weight), ≈ 6× the measured `ν = 5` value.
  * Grade: **[NE + M]** down to `ν = 5`. Below 5 it is **[M]**, a held extrapolation that is not measured, and that sub-window stays
    [I] in the strict sense.
* **`J = 32`: [I].** There is no validated node below 8.3: both thimble seeds fail the `Im` scalar check, and their `Re` values differ
  by `1.8·10⁴`.
* **Verdict for the window as posed: [I].** No validated method reaches the window's weight (`ν ≈ 2–5`) at the needed precision. That precision
  is `|Γ_c| ≲ 2000 c₂`, i.e. relative state/reference errors `≲ 1–2 %` on a validated thimble sum.
* The available thimble sums are not validated at `J = 28` or `J = 32` (scalar off by `2.6σ` and `6.7σ / 19σ`). Their seed-to-seed
  spread is `≈ 10⁴ c₂`.
* Brute force would need `≳ 25–100×` the statistics per node (1–4 core-days per node) **and** a fix of the systematic thimble failure
  (masked points, A₋ convergence). The latter dominates at `J = 32`.
* What is supported: every validated datum (`ν ≥ 8.3`) is `O(10³)` at `ν = 8.3` and `O(10²)` above. A failure needs the
  beyond-leading-Regge `Γ_c` to average `≲ −5400 / −4250 c₂` over `ν ≈ 2–5`, i.e. ≈ 4–5× the most negative validated value and
  of the same sign.

## 4. Remaining open items (value order)
1. **The `ν < 5` window at 300 for `J ≥ 28`** [I], with `J = 32` the most fragile state (no validated node below 8.3). It needs a validated method below the fold: a thimble sum that passes the scalar
   check at `ν ≈ 3–5` (fix masking/A₋ tracking; check the saddle set at `J = 32`), or an analytic small-`ν` bound on the
   beyond-Regge `Γ_c`. The same window governs all `J ≥ 24` at 300 (EE_DG300 §8).
2. **A derivation of `Γ_OO ≡ Γ_++`** (CHL `X^O_same` vs even `++` helicity amplitudes under `Ĉ'_EE`, `Ĉ₃`). This would turn the
   OO block from [D given NE] into [D] everywhere.
3. Direct small-`ν` EE nodes for even `J = 6–10` and odd `J = 7–15` at 300 (still convention A, positive, min `+0.307` / `+0.630`).
   Recipe: `oo_points.py J 400 x…` (even; it gives both blocks) and `HC_METHOD=n3 hc_points.py J 400 x…` (odd).
4. `H_edge` on the uncomputed rows (HYPOTHESES_CLOSE item 3), and `J = 48` at 1000 against direct values.

## 5. Hand-over
* **Code** (`ads4/code/spinfunc/`):
  * `oo_points.py` — direct 3×3 points at even `J`;
  * `oo_rows.py` — `K₂` OO rows, per state;
  * `oo_chain.py` — OO block in the z4 chain, plus run0 = EE block. Env: `OO_J`, `OO_EST`, `OO_JK`, `OO_DIRECT`, `OO_DESPIKE`,
    `OO_DUPNL`, `OO_EESRC=hc|oo`, `OO_Z3EE`, `OO_RUN0`, `OO_TAG`, and all `Z4_*` / `Q3_*`;
  * `oo_check.py` — per-node validation;
  * `oo_t.sh` — variants and jackknife; `oo_sum.py` — summary;
  * `oo_window.sh` — window stress runs through `hc_chain.py z4`;
  * `oo_merge.py` — pooling of thimble seeds;
  * queues `oo_queueD/D2/D3.sh` (direct), `oo_queueK.sh` (`K₂`), `oo_queueT/T2.sh` (thimbles).
* **Caches:**
  * `cache/oo_pts_{EE,C3H}_J<J>_D<D>.pkl` (+ `oo_raw_*`) and `cache/oo_p3_{EE,C3H}_J<J>_D<D>.pkl`;
  * `cache/oo_p2shadow_*`, `cache/oo_exact.pkl`;
  * `cache/oo_chain_<tag>.{txt,pkl,log}`. Tags: `f_*` = main, `m3_*` = `J = 18–22` with hc EE nodes, `val410*` = validation,
    `val18A` = superseded `J = 18`;
  * `cache/oo_check.txt`;
  * thimble caches: `cache/z2oo_pts_*_J32` and `z2oo2_pts_*_J32` (seeds 4242 / 6161), `z2oo3_pts_*_J28` (6161), and the pools
    `z2oom_pts_*_J{28,32}` and `z2oom2_pts_*_J28`;
  * `cache/hc_chain_z4_oowin_<J>_<src>_n<numin>_<variant>.*`.
* **State at hand-over:** no process of mine is running. The fourth `J = 28` seed (`oo4`, 7171) was not run (stopped for time).
  `oo_queueT2.sh` resumes it: finished nodes are skipped. Re-pool afterwards with `python3 oo_merge.py 28 4.999254944 oom2 s2,oo3,oo4`.
* **Reproduce (each a few minutes):**
  * validation: `OO_J=4,6,8,10 OO_DESPIKE=0 Q3_EXTRA_M=1.02,1.05,1.1,1.15,1.2,1.3,1.4,1.5,1.75,2,2.5,3,4,5 python3 oo_chain.py`;
  * main: `./oo_t.sh f 12,14,16,18,20,22 OO_EESRC=oo "OO_Z3EE=12:400;14:400;16:400;18:336,400;20:400;22:336"`, then
    `python3 oo_sum.py f`;
  * window: `JK=1 ./oo_window.sh 32 x,oo "0 -2000 -4000 -8000 -15000"` and `./oo_window.sh 28 x,s2 "…"`, read with
    `HC_MMAX=1.2 python3 hc_sum.py z4 300 oowin_<J>_<src>_n0_<variant> J=<J>`.
* Analysis runs were made with the thimble worker SIGSTOPped, to stay at 2 busy processes.
