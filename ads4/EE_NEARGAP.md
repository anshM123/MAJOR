# `Ω^EE` near the gap and beyond: per-state positivity sweep with `EE` crystal rows

*19 Nov 2026, EE-near-gap agent. Unrefereed, not committed. No existing document was edited. New code (`ads4/code/spinfunc/`):
`w0_ee_common.py` (EE integrand factories, TPS and point-wise), `w1_rows.py` (+ `w1_queue*.sh`; K2 rows: EE/C2/C3H, even and odd,
reference-basis points), `w2_chain.py` (the S1_S4 q3 chain with EE rows; variants by env), `w3_direct.py` (direct RQMC for EE/C3H),
`w4_tables.py` (per-state tables). Caches `spinfunc/cache/w*` (`w1_pts_*`, `w1_p3_*`, `w1_shadow_*`, `w3_n1_*`, `w2_chain_<tag>.txt`,
`w2_recs_<tag>.pkl`; `w1x_outliertest_*` = K2 outlier test points, kept out of the loaders). Project modules imported read-only
(cache names redirected / loaders wrapped inside the process).*

Grades: **[P]**, **[D]**, **[NE]** numerical evidence, **[M]** model, **[I]** open — as in `EE_FUNCTIONAL.md`.

## 0. Summary

1. **`Δ_g = 1000`: `n = 1` for `Ω^EE` on every computed state** — even `J ≤ 32`, odd `J ≤ 15`, `m/Λ = 1.01 … 5` [NE + M].
   Min `λ/N_s = +0.2021`, binding state unchanged (`+−`, parity-odd, `J = 5`, `1.01Λ`). Graviton diagonal corrections are
   `≤ 0.09 N_s` (`J ≤ 28`); the off-diagonal is `≤ 0.001 N_s`. `J = 32` uses a 35-node row at its own near-gap mass (`Δ = 1040`,
   `J²/m = 0.99`). `J ≥ 40` is **not covered**: the K2 `+−` diagonal on `(1040, 40)` breaks down for `x ≥ 0.5` (`Γ` up to `10⁸`), the same
   failure S1_S4 saw for `J = 32` on 1600 [I].
2. **`Δ_g = 300`: `n = 1` holds for all even `J ≤ 22` and odd `J ≤ 15` in every variant** (worst `J = 22`: `+0.15`, `J = 18`: `+0.48`;
   binding `+−`, `J = 5`, `1.02Λ`: `+0.2052`) [NE + M].
   **For `J = 24, 28` the verdict depends on the uncontrolled small-`ν` window.**
   * In the q3 convention of S1_S4/u6 (natural-spline extrapolation of `Γ` below the lowest node, local `ψ₂` weights), both fail at
     the lowest twist-allowed mass `1.1Λ`: `−0.053 / −0.095`, with `dF_hh ≈ −1.0 N_s` and `t*_J = 0.9`. They pass from `≈ 1.15Λ`.
   * They pass if any one of the following is changed: the min-norm `X_min` weight (`+0.043 / +0.087`); `Γ` held below the lowest node
     (`+0.27 / +0.48`); `Γ` set to zero below it (`+0.70 / +0.93`); the correction cut at `ν ≲ 25` (`+0.91 / +0.93`).
   * The whole effect comes from `ν ≲ 25`, mostly from `ν < 17`, where there is no crystal data at all.
   * **So `n = 1` at 300 is [I] for `J ≥ 24`.** It is [NE + M] for `J ≤ 22`.
3. **The limiting term is the graviton diagonal correction concentrated at small `ν`, not a `J²/m` bulk effect.** At the worst
   near-gap state the bin `x < 0.075` (`ν ≲ 25`) carries `−0.005 / −0.011 / −0.009 / −0.013` of `dF₂ = −0.008 / −0.012 / −0.009 / −0.013`
   (`J = 16 / 18 / 20 / 24`, `++`). This absolute size is roughly `J`-independent, while `N_s ≈ 1024/b̂³` falls. That produces the
   apparent growth `|dF_hh|/N_h ≈ 0.014, 0.036, 0.079, 0.135, 0.22, 0.50, 0.48, 0.82, 1.01, 1.08` (`J = 8 … 28`, even).
   * The margin is crossed between `J = 22` (`t* = 1.2`) and `J = 24` (`0.9`), at `m ≈ 1.1Λ`, `b̂ ≈ 41–44`.
   * At 1000 the same window has `≈ (300/1000)⁴` of the weight and is irrelevant.
4. **Validation [NE].**
   * Direct RQMC at `(336, 18)`: K2 is adequate at the binding mass (`+0.466` direct vs `+0.450` K2, `q²` basis).
   * Direct state values at `(336, 24)`, `x = 0.075 – 0.4`, with K2 references (the n1 contour fails at the reference for `ν ≲ 35`)
     confirm the K2 state actions to `≤ 1 %`, and leave `J = 24` at `−0.063 / −0.066`. They cannot reach `ν < 25`.
   * K2 overestimates the negative correction at `ν = 25` (`J = 18`: K2 `−364` vs direct `−65` c₂ units).
   * `(1600, 20)` direct: 1000-gap `J = 20` goes `+0.937 → +0.952`.
5. **Join with the band (§6) [M].** The leading-Regge band model is positive for `EE` on the whole twist-allowed band, with the
   off-diagonal `≡ 0` (EE_FUNCTIONAL §4). Beyond leading Regge, the computed diagonal corrections at `m ≥ 1.3Λ` are:
   * at 300: `−0.04 … −0.49 N_s` (`J = 12 … 28`); every state passes, worst `+0.47` (`J = 24, 28`, `1.3Λ`);
   * at 1000: `≤ 0.005 N_s`.

   Remaining failing or uncovered regions:
   * at 300, `J ≥ 24`, `m ≲ 1.15Λ` (in the q3 convention only, [I]);
   * at 300, `J ≥ 32` near the gap (no row at the state mass);
   * at 1000, `J ≥ 40` near the gap (K2 `+−` breakdown at large `x`).
6. **Off-diagonal.** In `EE`'s own order-0 helicity basis the `EE` off-diagonal is `≤ 0.03 N_s` everywhere computed. The u6 convention
   (`q²`'s basis) shows `0.16 / 2.4 N_s` at `J = 24 / 28` (`Δ = 1600`); this is an order-0 basis mismatch amplified by `m` (§3).
   **Odd sector.** The `p3` integrand was run with `E` in place of `q`, `Ĉ₃` unchanged. Corrections at 300: `J = 9 / 11 / 13 / 15` →
   `0.036 / 0.053 / 0.090 / 0.24 N_s`. At 1000 they are `≤ 0.03`.
7. **Gauge (coordinator request).** Results are in §4.
   * 1000: all gauges agree to `10⁻⁴`; `n = 1` in both gauges.
   * 300, `J = 18 / 20`: the canonical gauge reproduces the min-norm/local verdict once `Γ` is held below the first node
     (`+0.54 / +0.59`).
   * 300 with the spline extrapolation, which `X_can/ψ₂ ≈ −170 … +55` amplifies: `J = 18 / 20` (`q²`-basis run) go to `−0.63 / −0.25`,
     and `J = 28` (EE basis) to `−1.54`.
   * `n` at 300: min-norm gauge `n = 1` (spline: marginal `+0.04` at `J = 24`; hold: `+0.31`). Canonical gauge: `n = 1` with hold,
     `n = 1.55` with the spline, so **[I]**.

## 1. What was computed

| item | content | grade |
|---|---|---|
| EE integrand | `w0`: `E = V_{0,42}²(H₂₄ − 2V_{2,04}V_{4,02})/(X₂₀X₄₀)`, `D = E⊗E`, l9 (even) and p3 (odd, full conserved set) integrands; reproduces u5's records exactly (`ΔM0 = ΔM2 = 0.0`); point-wise copy = TPS at the saddle to `2·10⁻¹⁴` | [DV] |
| even `EE` K2 rows | `J = 0 … 16` at the `q²` keys (400/1600 for `J ≤ 10`; 400/800/1600 for `J = 12–16`); u5's `J = 18, 20`; `J = 24, 28, 32` on 1600; near-gap-mass rows `(336, 22)`, `(336, 24)`, `(336, 28)`, intermediate `(500, 24)`, `(500, 28)`, 1000-gap `(1040, 32)` (35 nodes), `(1040, 40)` (19 nodes); `J = 24` on 1600 with 35 nodes. New `J` with their `C3H` rows. | [NE] |
| odd `EE` K2 rows | `J = 5, 7, 9` (200/400/800/1600), `11, 13` (400/800), `15` (336/400, + C3H). `J = 19` started, not finished. | [NE] |
| reference-basis points | `EE` (and `q²` where missing) at `Δ_ref`, `x_ref = 0.3–0.8`, `J = 4 … 44` | — |
| direct RQMC (`EE`) | `(336, 18)` `x = 0.075–0.2`; `(336, 24)` `x = 0.075–0.4`; `(1600, 20)` `x = 0.05–0.2` (state + reference each) | [NE] |
| chain | q3: S1_S4 conventions (aligned rows, `c_off = mΓ_off`, `Ĉ₃` without de-shift, `logK − iπ(Δ+J)`, local de-shift with `u_EE`), EE-own helicity basis, Hampel despiking. States: all exact-skeleton states, plus hybrid states (flat + leading-Regge curvature, scalar remainder 0) at `m/Λ = 1.02 … 5` for every `J` with rows | [M] |

Every record is stripped with `logK → logK − iπ(Δ+J)`. The check `Im logK ≡ π(Δ+J) mod 2π` holds on all records (`> 7000`).

## 2. Per-`J` results (`cache/w2_chain_final_xc0.txt` / `final_xc0.075`; J = 32 from `j32_1040`)

Worst state per `J` (min `λ/N_s` over all near-gap and hybrid states) and the most negative diagonal correction `dF_hh/N_h`. The
off-diagonal correction is `≤ 0.03 N_s` throughout.
Columns: **A** = q3 convention (main): spline extrapolation, local `ψ₂` weights. **B** = `Γ` held below the lowest node (local weights).
**C** = min-norm weights, spline. **D** = correction off for `x < 0.075`.

| `J` | 300 A [`m/Λ`] | 300 diag (A) | 300 B | 300 C | 300 D | 1000 A [`m/Λ`] | 1000 diag (A) |
|---|---|---|---|---|---|---|---|
| 0 / 2 | `+0.9999 / +1.000` | `≤ 0.0006` | | | | `+1.000` | `≤ 10⁻⁴` |
| 4 | `+0.2204` [1.02] | `0.0005` | | | | `+0.2103` [1.01] | `0.0000` |
| **5 (odd)** | **`+0.2052`** [1.02] | `0.0009` | `+0.2052` | `+0.2052` | `+0.2052` | **`+0.2021`** [1.01] | `0.0000` |
| 6 | `+0.3074` [1.25] | `0.004` | | | | `+0.2508` [1.01] | `0.0001` |
| 7 (odd) | `+0.7043` [1.10] | `0.011` | | | | `+0.3114` [1.01] | `0.0002` |
| 8 | `+0.4297` [1.05] | `0.014` | | | | `+0.4345` [1.05] | `0.002` |
| 9 (odd) | `+0.4376` [1.03] | `0.036` | | | | `+0.4541` [1.03] | `0.0007` |
| 10 | `+0.6027` [2.0] | `0.036` | | | | `+0.5881` [1.01] | `0.001` |
| 11 (odd) | `+0.7504` [1.04] | `0.053` | | | | `+0.4099` [1.02] | `0.002` |
| 12 | `+0.6348` [1.06] | `0.079` | | | | `+0.6808` [1.06] | `0.003` |
| 13 (odd) | `+0.7287` [1.10] | `0.090` | | | | `+0.6045` [1.02] | `0.003` |
| 14 | `+0.6916` [1.05] | `0.135` | | | | `+0.7080` [1.04] | `0.003` |
| 15 (odd) | `+0.7084` [1.05] | `0.239` | `+0.749` | `+0.754` | `+0.708` | `+0.8696` [1.2] | `0.026` |
| 16 | `+0.7113` [1.10] | `0.222` | | | `+0.811` | `+0.7483` [1.03] | `0.005` |
| 18 | `+0.4759` [1.07] | `0.497` | `+0.542` | `+0.499` | `+0.829` | `+0.7431` [1.02] | `0.008` |
| 20 | `+0.4947` [1.07] | `0.481` | `+0.594` | `+0.527` | `+0.867` | `+0.8672` [4.0] | `0.010` |
| 22 | `+0.1470` [1.10] | `0.820` | `+0.363` | `+0.213` | `+0.896` | `+0.8449` [1.05] | `0.086` |
| **24** | **`−0.0534`** [1.10] | **`1.012`** | `+0.267` | `+0.043` | `+0.910` | `+0.9098` [5.0] | `0.031` |
| **28** | **`−0.0950`** [1.10] | **`1.083`** | `+0.480` | `+0.087` | `+0.934` | `+0.9349` [5.0] | `0.02` |
| 32 | no row at the 300 state mass | — | | | | `+0.967` (1040 row, 35 nodes; diag `+0.05`) | — |
| 40 | — | — | | | | **[I]** (K2 `+−` breakdown, §5) | — |

**`n`:**
* 300: `n = 1.11` in A (the `J = 24, 28` states at `1.1Λ` fail, `1.15Λ` passes); `n = 1` in B, C and D. Grade **[I]**.
* 1000: `n = 1` in every variant, min `+0.2021` (`+−`, `J = 5`, `1.01Λ`). Grade [NE + M] for `J ≤ 32`.

**Where the margin is crossed (A).** `t*_J` (factor on the correction at first failure) at 300 is `4.5 / 2.0 / 2.1 / 1.2 / 0.9 / 0.9`
for `J = 16 / 18 / 20 / 22 / 24 / 28`. Crossing is at `J ≈ 23`, `m ≈ 1.1Λ`, `b̂ ≈ 42`.
* Profile in `m/Λ` for `J = 24 / 28` (rows 336/500/1600):
  * `1.1: −0.05 / −0.10`; `1.15: +0.12 / +0.07`; `1.2: +0.29 / +0.27`; `1.3: +0.51 / +0.49`;
  * `1.5: +0.75 / +0.76`; `2: +0.94 / +0.93`; `≥ 2.5: ≥ 0.91`.
* That failing region `m ∈ [m_twist ≈ 1.08Λ, 1.15Λ]` is **inside** the leading-Regge band. The band model covers it but drops exactly
  this beyond-Regge diagonal, so it does not rescue these states.

## 3. Helicity basis [NE + D]

* On `Δ = 1600` at `J = 24 / 28`, the `q²` basis (u6) leaves an `EE` off-diagonal of `0.16 / 2.44 N_s`. `EE`'s own rank-1 order-0
  basis gives `0.002 / 0.0008`. The diagonals are unchanged (`0.200`).
* Reason: `Γ_off = m·Re(…)`, so any order-0 basis mismatch is amplified by `m`.
* The chain replaces the crystal's order-0 diagonal by the model base. That is consistent only in the basis where the functional's own
  order-0 matrix is diagonal: `EE`'s basis for `Ω^EE`.
* For `J ≤ 20` the two bases agree on the verdict, to within `0.03`.

## 4. Gauge of the `Ĉ′_EE` weight at small `ν` (`W2_GAUGE`; coordinator request)

The node weights `p2 = ψ₂·w` of the beyond-leading-Regge terms are rescaled by `Re X_g/ψ₂` at `ν = pΛ ≤ 25`. These terms are the
graviton correction `dF` and the scalar remainder. The flat `B₂` and the leading-Regge curvature are gauge-invariant by Lemma N and are
not rescaled. `X_min` is taken from u3, `X_can` from x2.

| `Δ_g` | `J` | local `ψ₂` (spline / hold) | min-norm (spline / hold) | canonical (spline / hold) |
|---|---|---|---|---|
| 300 | 5 (odd, binding) | `+0.2052` | `+0.2052` | `+0.2053` |
| 300 | 15 (odd) | `+0.708 / +0.749` | `+0.754 / +0.754` | `+0.619 / +0.619` |
| 300 | 18 | `+0.476 / +0.542` | `+0.499 / +0.556` | `+0.826 / +0.543` (q²-basis run: `−0.634`) |
| 300 | 20 | `+0.495 / +0.594` | `+0.527 / +0.612` | `+0.839 / +0.595` (q²-basis run: `−0.247`) |
| 300 | 22 | `+0.147 / +0.363` | `+0.213 / +0.396` | `+0.685 / +0.378` |
| 300 | 24 | `−0.053 / +0.267` | `+0.043 / +0.310` | `+0.560 / +0.288` |
| 300 | 28 | `−0.095 / +0.480` | `+0.087 / +0.539` | **`−1.544` (1.2Λ)** / `+0.514` |
| 1000 | all computed | as local, to `10⁻⁴` | same | same |

* The spline columns depend on the natural cubic spline's extrapolation of `Γ` below the lowest computed node (`x = 0.05`, `ν ≈ 17`).
  There is no crystal data below that node.
* `X_can/ψ₂` is `≈ −170` at `ν = 5` and `+55` at `ν = 10`. The canonical gauge therefore makes the extrapolated region decide.
  The `q²`-basis runs at `J = 18, 20` failed only through `Γ_off` there; with a cut at `ν ≤ 15` they give the same `−0.61`.
* With `Γ` held below the lowest node, all three gauges agree to within `≈ 0.05` and give `n = 1` at 300.
* **Verdict.** `n` in the min-norm gauge: 300 `n = 1`, marginal in the spline variant; 1000 `n = 1`. `n` in the canonical gauge:
  300 **[I]** (`n = 1.55` with the spline, `n = 1` with hold); 1000 `n = 1`.
* Recommendation: keep the min-norm gauge as the stated convention for the heavy side, unless the `ν ≲ 20` beyond-Regge correction is
  computed.

## 5. Validation of K2 for `EE`; x-resolution and outliers

* **Direct RQMC at `(336, 18)`** (`J²/m = 0.97`).
  * The off-diagonal agrees with K2 to `≤ 0.03`.
  * The diagonal `Γ_hh` differs by up to ~300 c₂ units at `x = 0.075–0.2`. At `x = 0.075` direct gives `−65`, K2 gives `−364`
    (RQMC error 23).
  * Chain at 300, `J = 18`: K2 `+0.450`; with the direct row `+0.466`; with direct state values and K2 reference (`ν ≥ 20`) `+0.528`.
* **Direct at `(336, 24)`** (`J²/m = 1.73`).
  * The n1 contour fails at the heavy **reference** for `ν ≲ 35`. Direct/K2 is `0.15 / 0.83` at `ν = 25 / 33`, with 52738 / 13729
    masked samples. The pure direct row is therefore unusable (`a_h` ×3.3).
  * State points agree with K2 to `≤ 1 %`. One exception is a K2 outlier at `x = 0.3`, `++` (direct/K2 `0.872`; K2 `M₂/M₀ = −0.62`).
  * Mixed chain (direct state, K2 reference; `W2_DSTATE`): `J = 24` at 300 gives `−0.063` (`ν ≥ 30`) and `−0.066` (`ν ≥ 20`).
* **Direct at `(1600, 20)`.**
  * At `x = 0.05 / 0.1` direct gives `+486 / −656`, K2 gives `−33 / −6` (errors 120–185).
  * Effect at 1000 is small: `+0.937 → +0.952`.
* **x-resolution** (`J = 24`, `Δ = 1600`, 19 vs 35 nodes).
  * The 19-node grid (spacing 0.05 against a `P_J` period of 0.11–0.13) aliases the oscillating `Γ_hh`.
  * `dF/N` at 300, 1.1Λ, with the 1600 row held: 19 nodes `0.198`; dense with outliers removed `0.036`.
  * On 1600 rows at `J ≥ 24`, 19-node diagonal corrections are uncertain by O(100 %) of themselves. They are small there.
* **K2 outliers.**
  * Isolated nodes reach `Γ` of `10⁴–7·10⁴` against neighbours of `~300`. Example: `J = 24`, 1600, `x = 0.6125–0.6312`, four extra
    nodes, erratic node to node.
  * They are not visible in `|M₂/M₀|` in the structure basis.
  * The main chain drops them with a Hampel filter (window ±3, 6× max(MAD, median|Γ|)). Dropped nodes are listed in the chain output:
    mostly `x = 0.8`, plus `J = 24`/1600/0.625.
  * Kept, they turn `J = 24` into a spurious `−7 N_s`.
* **K2 `+−` breakdown at large `J`, large `x`.**
  * `(1040, 40)`: `Γ_{+−}` is `7·10³ … 10⁸` at `x = 0.5–0.7`.
  * The same holds on `(1600, 32)` (S1_S4). `J = 32` at 1000 is therefore taken from the 1040 row alone; adding the broken 1600 row
    gives a spurious `−3.6`.
  * `J ≥ 40` at 1000 needs another method [I].

## 6. Joined map (near-gap chain + leading-Regge band) [M]

* **Leading-Regge band** (EE_FUNCTIONAL §4, model M′).
  * Off-diagonal `≡ 0`.
  * NF02 is positive on the whole twist-allowed band (`β ≤ 2.2`, `b̂ ≤ 2Λ`).
  * Min `0.2056` (300 / 1000), at the near-gap flat corner. Large `b̂` sits at `0.93–1.00`.
* **Beyond-Regge diagonal at `m ≥ 1.3Λ`** (hybrid states with the computed rows; interpolated in `1/m` between the near-gap-mass row,
  500/800 and 1600; held above the largest row).
  * 300, most negative `dF_hh/N_h` at `1.3Λ`, for `J = 12 / 14 / 16 / 18 / 20 / 24 / 28`: `−0.04 / −0.07 / −0.12 / −0.19 / −0.26 /
    −0.46 / −0.49`.
  * 300, `J = 24 / 28`: `1.5Λ` `−0.22`; `2Λ` `−0.06 / −0.04`; `≥ 3Λ` `≤ 0.01`.
  * 300, every one of these states passes; the worst is `+0.49` (`1.3Λ`).
  * 1000: `≤ 0.005 N_s` at `m ≥ 1.3Λ` for `J ≤ 32`.
* **Combined per-state map.**
  * `Δ_g = 1000`: positive on every computed state, `m/Λ ∈ [1.01, 5]`, even `J ≤ 32`, odd `J ≤ 15`. Beyond `J = 32` the band model
    covers `m ≥ 1.3Λ`, where the computed corrections are `≤ 0.5 %`; that coverage is [M].
  * `Δ_g = 1000`, uncovered: `J ≥ 40` near the gap (`m ≲ 1.3Λ`) [I].
  * `Δ_g = 300`: positive on every computed state with `J ≤ 22`.
  * `Δ_g = 300`, `J = 24, 28`: fail for `m ≲ 1.15Λ` in convention A only.
  * `Δ_g = 300`, `J ≥ 32`: not computed near the gap.
* **How far the 300 failure would extend [M].** The small-`ν` correction is `≈ −0.01` in absolute size and roughly independent of `J`
  and `m`, while `N_s` falls with `b̂`. In convention A every `J ≥ 24` near the gap fails. Further from the gap it fails while
  `N_s(b̂, m) ≲ 0.01`; the computed `1.3Λ` margins (`+0.49`) show it does not reach `1.3Λ` for `J ≤ 28`.

## 7. Graded summary

| statement | grade |
|---|---|
| w0 EE integrand = u5 (TPS) and = point-wise copy (direct) | [DV] |
| `Ω^EE` per-state positive at 1000 on every computed state (`J ≤ 32` even, `≤ 15` odd); min `+0.2021`, binding `+−`, `J = 5` | [NE + M] |
| `Ω^EE` per-state positive at 300 for even `J ≤ 22`, odd `J ≤ 15`; binding `+−`, `J = 5`, `+0.2052` | [NE + M] |
| `J = 24, 28` at 300 near the gap: fail in the q3 convention (`−0.05 / −0.10`), pass under every small-`ν` variant; decided at `ν ≲ 25` (mostly `ν < 17`) | [NE + M] for each variant; verdict [I] |
| `n = 1` at 300 | **[I]** |
| `n = 1` at 1000 (computed `J ≤ 32`) | [NE + M]; `J ≥ 40` near the gap [I] |
| K2 for `EE` validated at `(336, 18)` and `(1600, 20)` (verdict level); direct state at `(336, 24)` confirms K2 state values `≤ 1 %` | [NE] |
| EE off-diagonal beyond leading Regge `≤ 0.03 N_s` (own basis); u6's `q²`-basis off-diagonal is a basis artefact | [NE + D] |
| canonical gauge: no effect at 1000; at 300 decided by the extrapolation below `ν ≈ 17`; with `Γ` held, all gauges agree | [NE + M]; canonical `n` at 300 [I] |
| 19-node grids alias `Γ_hh` at `J ≳ 24`; isolated K2 outliers; K2 `+−` breakdown at `J ≥ 32–40`, `x ≥ 0.5` | [NE] |

## 8. Hand-over

* **Reproduce the main runs:**
  `Q3_JEVEN=0,2,…,28 Q3_JODD=5,…,15 Q3_EXTRA_M=1.02,…,5 W2_ROWS=200,336,400,500,800,1600 W2_BASIS=ee W2_MODES=ee W2_TAG=x python3 w2_chain.py 300 1000`.
  For `J = 32` at 1000 use `W2_ROWS=1040`. A run takes ~5 min on one core.
* **Variants:**
  * `W2_MODES=q2,ee`, `W2_BASIS=q2`;
  * `W2_GAUGE=none|min|can` (+ `W2_GNU`);
  * `W2_XLOW=spline|hold|zero`, `Q3_XCUT=0.075`;
  * `W2_ROWS`, `W2_DMIN`, `W2_INTERP=linear`, `W2_XSUB=19`, `W2_DROPX=J:Δ:x/x`, `W2_DESPIKE=0|1`;
  * `W2_DSTATE=1` (+ `W2_DS_NUMIN`), `W2_DIRECT=1`.
* **Tables:** `python3 w4_tables.py <tag>`.
* **Rows:** `W1_TAG=_x python3 w1_rows.py even|odd EE,C3H J-list Δ-list [x-list|keys]`. `w1_rows.py refb J-list` is needed for
  `W2_BASIS=ee`.
* **Direct:** `python3 w3_direct.py run J Δ x-list [EE|C3H]`; analyse with
  `W3_NUMIN=30 W2_BASIS=ee python3 w3_direct.py analyse "J:Δ,…"`.
* **Next steps, in order of value:**
  1. The `EE` beyond-Regge correction at `ν ≲ 25` near the gap at 300. Options are the n3 curved contour or thimble machinery, or a
     direct method that holds at the heavy reference for `ν ≲ 35`. It decides `n` at 300 and the gauge question.
  2. Dense (≥ 35-node), outlier-screened rows at the state masses for `J ≥ 32` at 300.
  3. A method for the `+−` diagonal at large `x` for `J ≥ 40` (1000), e.g. K3 or direct.
  4. Odd `J = 17, 19` rows. `J = 19` at 336 was started, its cache is partial, and it is not used.
