# `Ω^EE` at `Δ_g = 1000`: the large-`J` near-gap region (even `J ≥ 32`, odd `J ≥ 17`)

*20 Nov 2026, EE-1000-large-J agent. Unrefereed. No existing document or code file was edited. New code
(`ads4/code/spinfunc/`):
* `ee1000_common.py`: EE-mode setup and the Jacobi fix.
* `ee1000_rows.py`: K₂ rows.
* `ee1000_direct.py`: direct RQMC points.
* `ee1000_chain.py`: hybrid-state chain and shape model.
* `ee1000_w2rerun.py`: re-run of EE_NEARGAP's main chain.
* Checks and tables: `ee1000_checks.py`, `ee1000_directcmp.py`, `ee1000_rawcmp.py`, `ee1000_shapecheck.py`, `ee1000_band.py`,
  `ee1000_smallnu.py`, `ee1000_tables.py`.
* Scripts: `ee1000_ensemble.sh`, `ee1000_queue*.sh`, `ee1000_run.sh`.

Caches are in `spinfunc/cache/ee1000_*`. Project modules are imported read-only, and every patch is applied inside the process only.*

Grades: **[P]** proved, **[D]** derived, **[NE]** numerical evidence, **[M]** model, **[I]** open — as in `EE_NEARGAP.md`.

**STATUS (20 Nov, final for the even sector):**
* Even `J = 32, 36, 40` and odd `J = 17`: final, positive.
* `J = 48`: direct sub-row done (EE + C3H, 24 nodes). The verdict is **[I]** (§3a).
* Even `J ≥ 56`: [I] (§3).
* Odd `J = 19–31`: K₂ rows running in `ee1000_queueR1.sh` (EE queue: 31, 27, 21, 29; C3H queue: 23, 19, 25). Not yet analysed (§7).

## 0. Summary

1. **The "K₂ `+−` breakdown" of `EE_NEARGAP.md` §5 (and S1_S4's `J = 32` on 1600) was a numerical bug, not a failure of the Laplace
   series** [NE; checked against mpmath at 50 digits, §1].
   * The local de-shift evaluates `((1+c)/2)⁴ P^{(0,8)}_{J−4}(c)` at complex `c = 1 − 2(x ± i/m)²` with `scipy.special.eval_jacobi`.
     For `n = J − 4 ≳ 24` at `x ≳ 0.5` the complex branch of that function returns garbage.
   * Fixed in-process by a three-term recurrence.
   * Re-run of EE_NEARGAP's main chain with the fix (`cache/ee1000_w2rerun_main_fix.txt`):
     * 1000: unchanged for even `J ≤ 24` and odd `J ≤ 15`.
     * `J = 28`: max `|dF|/N` goes from 0.46 to 0.047.
     * `J = 32` (rows 1040 + 1600): `+0.954`, was `−3.63`. EE_NEARGAP's 1040-only `+0.967` had hidden a spurious `+33 N_s` `+−`
       diagonal at `1.05Λ`.
     * **For the 300 agent: the 300 `J = 28` state at `1.1Λ` (convention A) moves from `−0.095` to `−0.310`.**
2. **Per-`J` verdict at `Δ_g = 1000`, even `J ≥ 32`, own-mass rows** (§4; `cache/ee1000_tables_f1.txt`, `_fd40.txt`).

   | `J` | row (`J²/m`) | K₂ status | min `λ/N_s`, `m/Λ ∈ [m_gap, 5]` | max `|dF_hh|/N_h` | grade |
   |---|---|---|---|---|---|
   | 32 | 1040 (0.99) + 1600 | converged | `+0.904` (main); `+0.79 … +0.91` over the ensemble | 0.10 (0.21) | [NE + M] |
   | **36** | 1036 (1.25) | converged | **`+0.967`** in every variant | 0.044 (0.07; spline 0.16) | **[NE + M]** |
   | **40** | 1040 (1.54) | marginal, checked against direct | **`+0.979`** (K₂ and direct sub-row alike); `≥ +0.905` over all non-aliased variants | 0.074 (K₂) / 0.048 (direct) | **[NE + M]** |
   | 48 | 1048 (2.20) | **not converged**: `|M₄| ≈ |M₂|`; shape-fit residual 20–33 %; K₂ gives `+4 … +6 N_s`, which is meaningless | **undetermined**: with direct data on `x ∈ [0.1, 0.7]` the integrated correction depends on the representation (`−5 … +40 N_s`), see §3a | — | **[I]** |
   | 56, 64 | — (3.0, 3.9) | **divergent** (`J = 64`: `‖M_K − M_dir‖` grows K₀ → K₂, `|M₄/M₀| > 1`); direct RQMC errors `2·10⁻³–2·10⁻²` ↔ `Γ` errors `4·10³–3·10⁴` | not computable here | — | **[I]** |

   * Odd `J = 17` (row 1017): `+0.839`; corrections `≤ 0.008` [NE + M].
   * Odd `J = 19–31`: queued. K₂ is converged there (`J²/m ≤ 0.93`).
3. **There is no `J`-crossover above which the leading-Regge band model "takes over" near the gap. The opposite holds** [NE + M]:
   * Near the gap the weight of any beyond-Regge correction relative to `N_s` grows like `J³`, since `N_s ≈ 1024/b̂³` while the
     correction weight is `∝ 1/m²`.
   * K₂'s expansion parameter is `J²/m`.
   * So the band model is *least* reliable at large `J` near the gap.
   * The crossover is in `m` at fixed `J`. For `J ≤ 40` the measured corrections are `≤ 0.02 N_s` at `m ≥ 1.3Λ`, `≤ 0.01` at `1.5Λ`
     and `≤ 0.004` at `2Λ`.
   * For larger `J`, K₂ is usable only where `J²/m ≲ 1.5` (`m ≳ J²/1.5`, e.g. `≳ 2.7Λ` at `J = 64`). A K₂ row for `(2900, 64)` is
     queued to test that.
4. **The small-`ν` window `ν < 50`** (below the first reliable node; K₂ has the wrong sign there per EE_DG300_SMALLNU, and direct RQMC
   cannot reach `m² × 10⁻⁴` precision at `m ≈ 10³`) is uncomputed.
   * A uniform `Γ = C` there flips the gap state at `C* = −4600 / −2800 / −1300 / −740 / −480 / −270 / −130 / −90` for
     `J = 40 / 48 / 64 / 80 / 96 / 128 / 200 / 300` (`cache/ee1000_smallnu_50.txt`).
   * For comparison, the direct values at 300 (`ν = 12–21`, `J = 24, 28`) are `+440 … −30`.
   * So for `J ≲ 48` this window is harmless unless `Γ ≲ −3000` [M]. For `J ≳ 100` it is a genuine [I].
5. **`n` at `Δ_g = 1000`.**
   * **`n = 1` on every computed state**: even `J ≤ 40`, odd `J ≤ 17`, `m/Λ` from each `J`'s gap mass to 5.
   * The minimum is `+0.2009` at the `J = 5` gap-mass hybrid (`m = 1.0048Λ`, not in the exact-skeleton grid) [M]; it is `+0.2021` on
     the exact-skeleton states [NE + M].
   * **What remains [I] for the `Δ_g = 1000` statement:**
     * (a) even `J ≥ 48` at `m ≲ J²/1.5`: K₂ not converged; at `J = 48` a direct sub-row does not fix it either (§3a);
     * (b) odd `J = 19–31` (queued; expected unproblematic);
     * (c) the window `ν < 50` for all large `J`, decisive for `J ≳ 100`;
     * (d) `J > 64` at any `m` with `m < J²/1.5`.
   * The leading-Regge band model alone is positive on all of these (`λ₀/N_s ≥ 0.90` for `J ≤ 300`, §5) [M], but the beyond-Regge
     correction there is not controlled.

## 1. The Jacobi bug

`cache/ee1000_checks_jacobi.txt`. The table gives the relative error of `P^{(0,8)}_n(c)` at `c = 1 − 2(x + i/1040)²`, against mpmath at
50 digits.

| `n` | scipy, `x = 0.5` | scipy, `x = 0.7` | scipy, `x = 0.8` | recurrence (all `x ≤ 0.9`) |
|---|---|---|---|---|
| 20 | `9·10⁻⁹` | `1.5·10⁻⁵` | `5.6·10⁻⁵` | `≤ 2·10⁻¹⁵` |
| 24 | `9·10⁻⁷` | `1.3·10⁻³` | `5·10⁻²` | `≤ 2·10⁻¹⁵` |
| 28 | `6·10⁻⁴` | `0.35` | `1.1` | `≤ 7·10⁻¹⁵` |
| 36 | `7·10⁻²` | `5·10⁴` | `2·10⁶` | `≤ 3·10⁻¹⁵` |
| 60 | `2·10⁸` | `5·10¹⁷` | `7·10²⁰` | `≤ 5·10⁻¹⁵` |

* Only the `+−` de-shift term is affected; the real-argument shapes and numpy's Legendre are accurate.
* `C3H` (no de-shift) never was affected.
* `(1040, 40)`: `Γ_{+−}` goes from `7·10³ … 4·10¹⁰` to `|Γ| ≲ 2·10³` (`cache/ee1000_checks_deshift.txt`).
* The fix is applied to `Q.A.shift` and `Q.P4.shift` (`ee1000_common`). It also enters `n1_analysis`, which calls `A.shift`.

## 2. Rows and direct points computed

| item | content | grade |
|---|---|---|
| K₂ rows, EE + C3H, state + `Δ_ref = 2·10⁴` | `(1036, 36)` 35 nodes `x = 0.05–0.8`; `(1040, 40)` 41 nodes `0.025–0.9`; `(1048, 48)` 35 nodes; odd `(1017, 17)` 35 nodes; refb `J = 48, 56, 64` | [NE] |
| direct RQMC, EE state (K₂ reference) | `(1040, 40)`: 25 nodes `x = 0.2–0.8` (step 0.025); `(1064, 64)`: `x = 0.3, 0.6`; `(1048, 48)`: running, EE + C3H, `x = 0.1–0.7` | [NE] |
| K₂ state points for the comparison | `(1064, 64)` at `x = 0.3, 0.45, 0.6` | [NE] |

States: hybrid states (flat NF02 + leading-Regge curvature, scalar remainder 0, as q3 `Q3_EXTRA_M`). The masses are
`m/Λ ∈ {m_gap(J), 1.05, 1.075, 1.1, 1.15, 1.2, 1.3, 1.4, 1.5, 1.75, 2, 2.5, 3, 4, 5}`, twist-allowed only, with
`m_gap = √(998(999+2J))/Λ` (`Δ = 1000 + J`). EE's own helicity basis, S1_S4 row alignment, `logK − iπ(Δ+J)` on every record, Hampel
despike, 256 `p`-nodes (64 checked; they agree to `≤ 0.002`).

## 3. K₂ convergence and direct validation

* **`(1040, 40)`, `J²/m = 1.54`** (`cache/ee1000_directcmp.txt`):
  * Direct and K₂ `Γ_{EE,h}` have the same oscillation.
  * K₂ is larger by 5–15 % at `x ≤ 0.45` and by 20–40 % at `0.5–0.65`. At `x ≥ 0.7` the `+−` values differ outright (direct
    `69 / −56 / 1520` against K₂ `−684 / −1727 / −1235` at `0.7 / 0.75 / 0.8`).
  * RQMC error: 20–60 at `x ≤ 0.6`, 70–280 above.
  * **Integrated:** with direct on `x ∈ [0.2, 0.8]` the gap-state corrections are `+0.007 / +0.048` (`++ / +−`), against K₂'s
    `+0.016 / +0.074`, in the main shape model. Over the shape variants `K = 4–6` they are `+0.007 … +0.096`.
  * The K₂ systematic at `J = 40` is therefore `≲ 0.05 N_s`, against a margin of 0.98 [NE].
* **Order ratios** (median `|M₂|`, `|M₄|` relative to the local envelope of `|M₀|`; `ee1000_rawcmp.py`):
  * `J = 36`: `0.020, 0.007`.
  * `J = 40`: `0.024, 0.014`.
  * `J = 48`: `0.016–0.029, 0.017–0.018`. `M₄ ≈ M₂`, so the series no longer decreases.
  * `J = 64`, `x = 0.3 / 0.6`: `|M₂/M₀| = 0.24 / 2.0`, `|M₄/M₀| = 1.4 / 4.2`. The deviation from direct grows `0.46 → 0.55 → 1.14`
    and `1.1 → 3.7 → 6.4` from K₀ to K₂ (divergent).
  * Direct RQMC relative errors at `J = 64`: `3.6·10⁻³ / 2.3·10⁻²`, which is `Γ ≈ 4·10³ / 2.6·10⁴` c₂. That is useless at the needed
    `≲ 10²`.
  * **Criterion: K₂ is usable for `J²/m ≲ 1.5` [NE].** Above that, near the gap, neither K₂ nor this direct method controls `Γ`.
* `x > 0.8` (K₂ `Γ` grows to `10⁴` at `J = 40`) and `ν < 50` (`x < 0.05`; K₂ `Γ(ν = 26) = −1500` at `J = 40`, the wrong sign per
  the 300 direct values) are excluded from the main fit. They are variants `x09` and `withlo` in the tables.

### 3a. `J = 48` with a direct sub-row (`cache/ee1000_directcmp.txt`, `ee1000_chain_d48*.txt`, `ee1000_tables_d48.txt`) [NE]

* Direct EE and C3H state values at 24 nodes `x = 0.1–0.7` (step 0.025, K₂ reference). The RQMC error is 10–150 c₂.
* Direct minus K₂ is up to `~10³` (e.g. `x = 0.25`, `++`: `−432` vs `−1209`).
* The direct `Γ` has a clear **non-oscillating positive component** of `≈ +500` c₂ (`++` and `+−`) on top of the `P_J` oscillation.
  With the 0.025 grid the oscillation is sampled about 2.6 times per period. The shape model does not represent it: residuals are
  30–50 % even on the direct data, so the class `c D + s xD′ + r` with smooth `c, s, r` fails at `J = 48`, `J²/m = 2.2`.
* The gap-state correction therefore depends on the representation:
  * shape, `K = 4`, direct window only: `+5.2 / +0.02 N_s` (min `λ/N_s` `+0.68`, at `1.75Λ`);
  * `K = 3`: min `+0.96`;
  * spline: `−4.7 / +39.7` (its natural-spline extrapolation below `x = 0.1` alone gives `−17.8`);
  * half-node fits: `±10²`.
  * With the K₂ end nodes the result is `+9.3 / +2.1`.
* **Verdict for `J = 48`: undetermined [I].** The positive offset is suggestive (it adds positively), but not established. Settling it
  needs a denser direct row (`Δx ≲ 0.01`, and x down to 0.05 with a direct reference), or a method that does not need `Γ` point-wise.

## 4. Representation of the oscillating `Γ` (aliasing) and the ensemble

* At `J ≥ 40`, `Γ` oscillates with the `P_J` period, `≈ 2π/(2J+1)` in `x` at small `x` and shorter at large `x`. A cubic spline on the
  0.025 grid aliases:
  * spline variant at `J = 40`: `+1.09 N_s`;
  * `J = 48`: `+12.7 N_s`;
  * at `J = 40` the spline also jumps the Hampel-dropped node 0.575, a genuine outlier in both kinds.
* The main model is a least-squares **shape fit** per row and key: `Γ(x) = c(x)D_h(x) + s(x)·xD_h′(x)/(2J+1) + r(x)`, with `c, s, r`
  polynomials of degree `K = 4`. The oscillation is carried by the exact LO shape, and `c, s, r` are held beyond the last node.
  * Residuals: `≤ 5 %` (K₂, `J ≤ 40`), `≤ 1 %` (direct, `J = 40`); C3H at `J = 15`: `5·10⁻⁵`.
* Per-bin contributions to `dF` are `±0.1–0.5 N_s` and cancel. The interior bins agree between all models at `J ≤ 40`. The spread comes
  from the end regions, i.e. `x < 0.05` and `x > x_last`.
* Ensemble (`ee1000_ensemble.sh`): main (`K = 4`, `0.05 ≤ x ≤ 0.8`), `K = 5`, `K = 6`, `x ≤ 0.9`, spline, with the LO nodes, and 64
  `p`-nodes.

| `J` | uncorr. | main | K5 | K6 | x≤0.9 | spline | with LO | t*_J (main) |
|---|---|---|---|---|---|---|---|---|
| 17 (odd) | `+0.844` | `+0.839` | `+0.839` | `+0.839` | `+0.839` | `+0.841` | `+0.839` | `(+142, −∞)` |
| 32 | `+0.967` | `+0.904` | `+0.855` | `+0.788` | `+0.904` | `+0.913` | `+0.904` | `(+10, −36)` |
| 36 | `+0.967` | `+0.967` | `+0.967` | `+0.967` | `+0.967` | `+0.968` | `+0.967` | `(+∞, −23)` |
| 40 (K₂) | `+0.979` | `+0.979` | `+0.979` | `+0.979` | `+0.980` | (aliased) | `+0.914` | `(+∞, −13.5)` |
| 40 (direct `x ≥ 0.2`) | `+0.979` | `+0.979` | `+0.980` | `+0.980` | `+0.980` | (aliased) | `+0.905` | `(+∞, −21)` |

* Lower-`J` gap-mass states, which the exact-skeleton grid misses (hybrids, existing rows, `cache/ee1000_tables_low.txt`):
  * all positive;
  * the minimum is `J = 5` at `1.0048Λ`: `+0.201`;
  * `J = 4`: `+0.204`;
  * `J = 15` at `1.015Λ`: `+0.675` (its corrections are held from the 336/400 rows) [M].

## 5. Band model and small-`ν` thresholds

* **Band model** (`ee1000_band.py`; `cache/ee1000_band.log`): the uncorrected `λ₀/N_s` on the discrete states, `J ≤ 300`,
  `b̂ ≤ 2Λ`. These are the leading-Regge model M′ values [M].
  * `J ≤ 15`: min `+0.2009`.
  * `16–31`: `+0.746`.
  * `32–64`: `+0.962`.
  * `65–128`: `+0.967`.
  * `129–300`: `+0.902`.
  * 128 `p`-nodes are not enough at `J ≥ 200`: they gave a spurious `+0.06` at `J = 300`. 512 nodes are used for `J ≥ 100`.
* **Small-`ν` window** (§0.4) [M].

## 6. Graded summary

| statement | grade |
|---|---|
| EE_NEARGAP's "K₂ `+−` breakdown" (`J ≥ 32`) is scipy's complex `eval_jacobi`; fixed by recurrence | [NE] (vs mpmath) |
| with the fix, EE_NEARGAP's 1000 chain: `J = 28` corrections 0.047 (not 0.46), `J = 32` with both rows `+0.954`; `n = 1`, min `+0.2021` | [NE + M] |
| with the fix, 300 `J = 28` at `1.1Λ` (convention A): `−0.310` (was `−0.095`) | [NE + M] (for the 300 agent) |
| `J = 36`, `40` at 1000: positive at every mass, min `+0.967 / +0.979`, beyond-Regge corrections `≤ 0.07` | [NE + M] |
| `J = 40`: K₂ vs direct on `x ∈ [0.2, 0.8]` agree in the integrated correction to `≤ 0.05 N_s` | [NE] |
| K₂ usable for `J²/m ≲ 1.5`; not converged at `J = 48` near the gap, divergent at `J = 64` | [NE] |
| odd `J = 17`: `+0.839`, corrections `≤ 0.008` | [NE + M] |
| no large-`J` crossover to the band model near the gap: correction weight `/N_s ∝ J³` | [NE + M] |
| `n = 1` at 1000 on all computed states (even `J ≤ 40`, odd `J ≤ 17`; `m_gap … 5Λ`) | [NE + M] |
| `J = 48` near the gap: K₂ not converged; direct sub-row (0.025 grid) does not determine the integrated correction (`−5 … +40 N_s`) | **[I]** |
| even `J ≥ 48` near the gap (`m ≲ J²/1.5`); odd `J = 19–31` (rows running); `ν < 50` window at large `J` | **[I]** |

## 7. Hand-over

* **Running** (resumable; after a restart re-run exactly these two commands, and existing points are skipped):
  `nohup ./ee1000_queueR1.sh EE 31,27,21,29 > cache/ee1000_queueR1_EE.log 2>&1 &` and
  `nohup ./ee1000_queueR1.sh C3H 23,19,25 > cache/ee1000_queueR1_C3H.log 2>&1 &`.
  * Their `J = 48` direct part is done.
  * They are now on the odd rows: 35 nodes, about 70 min each, 3–4 rows per queue.
* Optional afterwards: `./ee1000_queueA6.sh`, the K₂ row `(2900, 64)`, which checks that the band takes over at `m ≳ J²/1.5`.
* **Analyse `J = 48` with direct:** `W2_DMIN=1000 EE1000_DSTATE=1 EE1000_JEVEN=48 ./ee1000_ensemble.sh fd48`. This uses the direct
  state values where they exist (EE and C3H) and K₂ elsewhere; restrict with `EE1000_XMIN=0.1 EE1000_XMAX=0.7`.
* **Odd rows:** `W2_DMIN=1000 EE1000_JODD=17,19,...,31 ./ee1000_ensemble.sh fodd`.
* Direct-vs-K₂ point tables: `python3 ee1000_directcmp.py "48:1048"`, `python3 ee1000_rawcmp.py "48:1048"`.
* Rules used: at most 2 cores. `ee1000_run.sh cmd` pauses one queue process while an analysis runs.
* **Next steps:**
  1. A method for `J²/m ≳ 2` near the gap: higher-precision direct with state/reference common random numbers and the scalar control
     variate of EE_DG300_SMALLNU §3, or a resummed (`J²/m`-exact) Laplace expansion.
  2. The `ν < 50` window at large `J`.
