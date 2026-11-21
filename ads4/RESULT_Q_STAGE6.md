# Result Q, Stage 6: full-region heavy coverage, cross-validation of the two heavy models, and gauge robustness of the de-shifted forward pieces

*21 Nov 2026, Stage-6 agent. Unrefereed, not committed by this agent. No existing document or code file was edited. New code:
`ads4/code/resultq/stage6_model.py`, `stage6_scan.py`, `stage6_scan_ext.py`, `stage6_fwexact.py`, `stage6_gauge.py`, `stage6_canon.py`,
`stage6_stress.py`, `stage6_wcheck.py`. New functional file: `ads4/code/resultq/functionals_q6_canonical.json` (a candidate gauge, **not adopted**, §2.2).
Outputs in `resultq/cache/stage6_*` (git-ignored). Read-only imports: `common_q.py`, `stage3_fwexact.py` (→ `exactskel/exactmom.py`),
`stage4_deshift.py`, `stage4_hj.py`, `stage5_deshift.py`, `stage5_hj.py`, `stage5_positivity.py`, `stage12_forward_regge.py`. Compute: about 2 CPU-hours
on ≤ 2 cores, ≤ 3 GB.*

Grades: **[P]** proved (incl. exact computer algebra), **[D]** derived with the stated assumption, **[DV]** derived and checked numerically,
**[M]** model / numerical, **[NE]** numerical evidence, **[I]** open. Flat = D = 4 amplitudes; AdS = the CFT₃ functionals. All minima below are
`FF/N` with Stage 0's normalisation `N = |B₂| + |B₃| + |a(∂B₄ − B₅)| + c|B₄|`, the one the LP margin δ refers to (Stage 5 used a larger `Nrm`;
both are given where they differ). Functionals: δ = 0.05 = `Q0-full-robust`; δ = 0.2 = `Q3-full-delta0.2`. Both forward pieces de-shifted throughout.

---

## 0. Answer in brief

1. **Coverage is now complete within the model hierarchy [M].**
   * One model, **HYB** (§1.1), covers every physical state: all integer J, the twist gap imposed, `m` from the gap to `10⁴Λ`, both channels.
   * HYB = flat-exact Stage-0 action (exact Wigner functions, flat graviton normal forms) + the leading-Regge AdS curvature of every piece. It is
     Stage 3's hybrid with the forward parts treated the same way, and it reduces to Stage 5's LR model at large J.
   * Previously uncovered regions are now covered:
     * J ≥ 25 with `m ≤ 4Λ`;
     * `m > 4Λ` with `b̂ < 20`;
     * off-grid `m` between Stage 3's 12 grid points;
     * J = 18, 22 and near-gap J ≥ 16 at Δ_g ≥ 300;
     * `J > 3000` with `b̂ ≤ 200`;
     * the scaling corner `m ≥ 1000Λ`.
2. **Cross-validation [M].**
   * HYB vs Stage 3's exact-skeleton model (EXA), same states: `≤ 0.005` (298 states, Δ_g = 100), `≤ 6·10⁻⁴` (300), `≤ 5·10⁻⁵` (1000).
   * HYB vs EXA on **39 new exact CHMRS states** outside Stage 3's box (J = 26–32 at 2.6–8Λ; J ≤ 24 at 6–20Λ; Δ_g = 100): `≤ 5·10⁻⁵`.
   * HYB vs LR in the overlap `b̂ ∈ [10, 30]`:
     * at `m ≥ 4Λ`: `≤ 0.018` (`b̂ ≥ 20`) and `≤ 0.077` (`b̂ ∈ [10, 20)`);
     * at Stage 5's own edge minima: `±0.012`.
   * LR fails near the gap (`m < 4Λ`) even at `b̂ ∈ [20, 30]`: differences up to 0.6. So Stage 5's "`b̂ ≥ 20` = LR-valid" holds only for `m ≳ 4Λ`.
3. **Full-region minima [M]** (§1.5): **δ = 0.05: +0.0478 / +0.0497 / +0.0491; δ = 0.2: +0.1985 / +0.1977 / +0.1970** (Δ_g = 100 / 300 / 1000).
   No negative heavy image anywhere. The new region `m > 4Λ, b̂ < 20` binds for δ = 0.2 (`+−`, J = 26, `m = 6.27Λ`), slightly below the nominal 0.2.
4. **Flat-space finding [NE].** The Stage-0 flat actions dip at small `u`, `b̂ ≈ 55–85` (`J ≈ 2700–4300`, `m ≈ 90–130Λ`), a band the Stage-0 LP grid
   (`J ≤ 2500`) never constrained:
   * δ = 0.05: **+0.025** (flat);
   * δ = 0.2: **+0.185** (flat).

   AdS curvature lifts them to `≥ +0.087` and `≥ +0.2185` at Δ_g = 1000, so they bind nowhere at the three Δ_g, but they set the Δ_g → ∞ floor.
5. **Gauge [DV + M]** (§2).
   * Every valid gauge gives the **same** heavy image on any covariant (shift-structured) raw image. Checked to `≤ 10⁻¹⁰` with 7 gauges of `Ĉ₅`
     and 5 of `Ĉ₄`, at the binding states, including beyond-LR finite-`m` skeleton functions. So in LR, HYB and EXA the minima are
     **gauge-independent (spread 0)**.
   * Gauge dependence needs a non-covariant beyond-LR part. With such terms of strength `K` (finite-J-constant scale, §2.3):
     * the `Ĉ₅` gauge is **immaterial**: it changes the minima by `≤ 0.003` at `|K| = 5`, being `u⁴/Λ²`-suppressed;
     * the `Ĉ₄` **derivative** piece is the sensitive one. Positivity survives `|K| ≲ 0.5–0.64` (δ = 0.05) or `≲ 1.5–2.0` (δ = 0.2), worst
       structure × gauge. The binding state is the `++`, J = 0 gap state.
   * A candidate canonical gauge ("shift-blind": extra condition `μ₁ = M₁`) exists and is stable, but it trades sensitivities rather than
     removing them, and for `Ĉ₄` it uses a near-null direction. **Not adopted; worst case over gauges reported.**
6. **Verdict** (§3).
   * Both δ are heavy-positive over the full region in the covariant models.
   * **Neither is gauge-robust** at the natural size `|K| ~ 1–5` of the non-covariant finite-J terms. δ = 0.2 tolerates about 3× more.
   * **Recommended: δ = 0.2 as the quoted (conservative) bound**: `A = 17.241`, `B_eff = −19.038 / −19.220 / −19.265`, `κ = 0.110`. δ = 0.05 is the
     in-model bound: `A = 11.448`, `B_eff = −12.804 / −12.938 / −12.970`, `κ = 0.0686`.
   * **No LP re-optimisation was needed** (no negative state). ABJM sits at **4.5 / 3.3 / 2.6 %** (δ = 0.2) and **6.8 / 5.0 / 3.9 %** (δ = 0.05) of
     the bound. Constants and placements are unchanged from Stages 3–5.
7. **New caveat [NE].** The stored double-precision de-shift weights realise the pure image only for `t ≲ 10` (`Ĉ₄` value: `|W − 1| = 7.5·10⁻⁴` at
   t = 10, 0.41 at t = 12) and `t ≲ 14` (`Ĉ₅`). Stage 4/5's β ≤ 30 positivity therefore refers to the exact (continuum) functionals. Realising
   them at `β > 10` needs an extended-precision solve [I].

---

## 1. Coverage map and closure

### 1.1 The model (`stage6_model.py`) [M]
* **State variables.** For a state `(Δ, J)`:
  * `m² = (Δ−J−2)(Δ+J−1)`, `u = Λ²/m² ≤ u_gap(J) = (Δ_g−1)/(Δ_g+2J−1)`, `b̂ = (2J+1)√u`;
  * `cosh β = 1 + b̂²/(2Λ²)`. **β depends on `b̂` alone**, so every AdS-curvature correction is a 1D function of `b̂`. This is exact
    bookkeeping, not an approximation.
* **LR (Stage 12/5).**
  * Smeared parts: `B₂ = 2I₂₀(β) − uI₂₂(β)`, `B₃ = σuI₃(β)`.
  * Forward parts: `B₄ = 2u²P₀`, `∂B₄ = −2u²Λ²D + c₄u³P₀`, `B₅ = σu³P₀`.
* **HYB = FLAT_EXACT + (LR − LR_flat).**
  * `FLAT_EXACT` is the Stage-0 action at `(u, J)`: `common_q.parts` with exact Wigner functions; `B₄ = 2u²`, `∂B₄ = (c₄ − 2JJ)u³`, `B₅ = σu³`.
  * `LR_flat` is the β → 0 limit of LR at fixed `b̂`: `P_{−½+iΛp} → J₀(b̂p)`, `P₀ → 1`, `Λ²D → b̂²/4`.
  * For the smeared parts this is exactly Stage 3's hybrid. The forward parts now get the same treatment, which removes the need for exact moments
    state by state.
  * Tables: 7 000–10 000 `b̂` nodes per (functional, Λ), with cubic splines. The same exact polynomial cos transforms and Mehler quadrature as
    Stage 5 are used for both the AdS and the flat columns.
* **EXA (Stage 3).** HYB smeared parts plus exact CHMRS scalar-skeleton forward moments `ex₄(0)`, `∂ex₄`, `ex₅(0)`.
  * Check: `stage6_model.hyb(..., fw=exact)` reproduces Stage 3's 298 `exact` values to `2.9·10⁻⁹` [DV].
* **Channels.** `++` (even J), `+−` (J ≥ 4). The parity-odd even-J block has the `++` image at this order (Stage 3 §2.1), which is not
  re-examined here [M].

### 1.2 Map before Stage 6 → closure

| region | before | now (model, file) |
|---|---|---|
| R1: J ≤ 24, `m ≤ 4Λ` | Stage 3 EXA on 12 `m`-values + gap (J ∈ {0, 2, 4–14, 16, …, 24}; at Δ_g ≥ 300 J ≥ 16 only for `m ≥ 1.5Λ`, J = 18, 22 dropped) | every integer J, 100 `m`-values in `[1, 2]`, dense above; every gap state J ≤ 300 (HYB, `stage6_scan`) + EXA on Stage 3's states |
| R2: J ≥ 25, `m ≤ 4Λ`, `b̂ < 20` (J = 25–39 at `m = 2.55–4Λ`) | **none** | HYB + EXA at J = 26, 28, 32 (2.6Λ, 4Λ) |
| R3: `m > 4Λ`, `b̂ < 20` (all J < `10m/Λ`, unbounded) | **none** in AdS (flat dense check only, [NE]) | HYB up to `m = 1000Λ` + EXA at J ≤ 24 (6, 10, 20Λ) and J = 26–32 (8Λ) |
| R4: overlap `b̂ ∈ [10, 30]` | not matched | both models on the same ≈ 285 000 states; differences §1.4 |
| R5: `b̂ ≥ 20`, J ≤ 3000 | LR (Stage 5; J continuous, gap not imposed) | HYB (finite-J corrections added) + LR alongside |
| A′: J > 3000, `b̂ ≤ 200`, `m ∈ [20, 600]Λ` | LR only; **LR drops the O(u) finite-J terms that dominate the flat dips there (§1.6)** | HYB with J ≤ 40 000 (`stage6_scan_ext`) |
| B′: J > 3000 with `b̂ > 200` or `m > 600Λ` | LR (coarse β grid at small β) | LR, dense in `b̂` (1200 extra β), gap imposed |
| scaling corner `m ≥ 1000Λ` (`b̂ < 6` at J ≤ 3000) | flat LP scaling row | `FF/N = 1.0000` at `m = 10³, 3·10³, 10⁴Λ` (direct); `S_AdS(b̂)/|S_flat(b̂)| ≥ 0.94` on the full `b̂` range |

**Not covered:** β > 30 (Stage 5's cut); `m > 10⁴Λ` beyond the scaling argument; the beyond-LR graviton-specific (non-skeleton) corrections. The
last are [I] throughout and are the subject of §2.3.

### 1.3 New exact CHMRS states (`stage6_fwexact.py`, Δ_g = 100) [DV]
* 39 states: J ∈ {0, 2, 4, 5, 6, 8, 12, 16, 20, 24} at 6, 10, 20Λ; J ∈ {26, 28, 32} at 2.6, 4, 8Λ.
* Numerics: two-fit agreement `≤ 2·10⁻¹¹` (value), `≤ 10⁻⁹` (slope); cancellation `≤ 1.1·10⁴`. Run time 1.4 s (J = 8) to 270 s (J = 32).
* J = 40 was stopped: Chebyshev overflow warnings in `mack.py`, so the result was not trusted and not used.
* EXA − HYB on these states: `≤ 9·10⁻⁸` (δ = 0.05, 33 cases), `≤ 5.3·10⁻⁵` (δ = 0.2, 59 cases).
* On Stage 3's states: `≤ 0.0033 / 0.0004 / 0.00005` (δ = 0.05) and `≤ 0.0050 / 0.0006 / 0.00005` (δ = 0.2). The maximum is always at `+−`, J = 4,
  just above the gap.
* So HYB's forward curvature (LR on the flat-exact base) is an excellent proxy for the exact skeleton everywhere except the J = 4 gap corner,
  where EXA is used [DV for the moments, M for the identification].

### 1.4 Cross-validation in the overlap (`stage6_scan`, all region-A states with `10 ≤ b̂ < 30`)
`max |FF/N_HYB − FF/N_LR|`, Δ_g = 100 / 300 / 1000:

| band | `m < 4Λ` | `4Λ ≤ m < 30Λ` | `m ≥ 30Λ` |
|---|---|---|---|
| δ = 0.05, `b̂ ∈ [10, 20)` | 0.51 / 0.60 / 0.70 (LR invalid) | 0.063 / 0.076 / 0.077 | 0.023 / 0.033 / 0.035 |
| δ = 0.05, `b̂ ∈ [20, 30)` | 0.076 / 0.39 / 0.59 (LR invalid) | 0.011 / 0.015 / 0.017 | 0.008 / 0.016 / 0.018 |
| δ = 0.2, `b̂ ∈ [10, 20)` | 0.51 / 0.54 / 0.60 | 0.026 / 0.026 / 0.026 | 0.003 / 0.003 / 0.004 |
| δ = 0.2, `b̂ ∈ [20, 30)` | 0.047 / 0.41 / 0.63 | 0.011 / 0.013 / 0.013 | 0.006 / 0.009 / 0.009 |

* Mean differences are `≤ 6·10⁻⁴`.
* **At Stage 5's own minima** (`b̂ ≈ 20.8`, `m ≈ 19–27Λ`):
  * δ = 0.05: HYB **+0.214 / +0.072 / +0.049** vs LR +0.205 / +0.060 / +0.037;
  * δ = 0.2: HYB +0.280 / +0.210 / +0.200 vs LR +0.289 / +0.221 / +0.212.
* The LR column reproduces Stage 5:
  * `LR-S5` (Stage-5 convention: no gap, J continuous) gives +0.2065 at Λ = 100, vs Stage 5's +0.206;
  * region A's LR gives +0.0367 at Λ = 1000, vs Stage 5's +0.042 in `Nrm` on a coarser grid.
* **Finding [M].** LR is accurate to `≲ 0.02` for `m ≳ 4Λ`, `b̂ ≥ 20`. It is **not** valid near the gap even at `b̂ ∈ [20, 30]` (J = 10–13,
  `m ≈ 1.01–1.12Λ`), where it overestimates by up to 0.6. Those states are inside Stage 3's EXA box and are positive there (≥ +0.34).

### 1.5 Full-region minima (`stage6_scan.py`, `stage6_scan_ext.py`) [M]
Minimum `FF/N` over all physical states (gap imposed, both channels). The binding state's location is given; "EXA" is the exact-skeleton value
where it is lower.

| Δ_g | δ = 0.05: min (location) | δ = 0.2: min (location) |
|---|---|---|
| 100 | **+0.0478** EXA (`+−`, J = 5, `1.20Λ`, `b̂ = 9.2`); HYB dense +0.0489 (`+−`, J = 6, `1.475Λ`) | **+0.1985** EXA (`+−`, J = 5, `1.20Λ`); HYB dense +0.1994 (`+−`, J = 5, `1.187Λ`) |
| 300 | **+0.0497** HYB (`+−`, J = 6, `1.475Λ`, `b̂ = 8.8`); EXA +0.0499 | **+0.1977** HYB (`+−`, J = 26, `6.27Λ`, `b̂ = 8.5`; region R3) |
| 1000 | **+0.0491** HYB (`++`, J = 278, `26.8Λ`, `b̂ = 20.8`; Stage 5's edge); near-gap +0.0498 | **+0.1970** HYB (`+−`, J = 26, `6.27Λ`, `b̂ = 8.5`; region R3) |

By region (HYB; δ = 0.05 | δ = 0.2; Δ_g = 100 / 300 / 1000):

| region | δ = 0.05 | δ = 0.2 |
|---|---|---|
| R1 Stage-3 box | +0.0489 / +0.0497 / +0.0498 | +0.1994 / +0.1999 / +0.1999 |
| R2 J > 24, `m ≤ 4Λ`, `b̂ < 20` | +0.543 (`+−`, J = 31, 4Λ) | +0.705 (`++`, J = 28, 4Λ) |
| R3 `m > 4Λ`, `b̂ < 20` | +0.115 / +0.058 / +0.051 (`+−`, J ≈ 102–111, `m ≈ 14–15Λ`, `b̂ = 14.5`) | **+0.2036 / +0.1977 / +0.1970** (`+−`, J = 26, `6.27Λ`) |
| R5 `b̂ ≥ 20`, J ≤ 3000 | +0.214 / +0.072 / +0.049 | +0.280 / +0.210 / +0.200 |
| A′ J > 3000, `b̂ ≤ 200` | +0.829 / +0.428 / +0.087 (`+−`, J ≈ 3000–3100, `m ≈ 69–104Λ`) | +0.792 / +0.425 / +0.2185 |
| B′ LR remainder | +0.941 / +0.877 / +0.703 | +0.926 / +0.844 / +0.642 |
| gap states J ≤ 300 | +0.0500 (`++`, J = 0) | +0.2000 (`++`, J = 0) |
| β ≥ 2 (LR) | +0.951 / +0.983 / +0.995 | — (Stage 5: ≥ +0.94) |

**Error budget at the minima [M].**
* Forward model: HYB vs EXA `≤ 0.005`.
* Smeared parts: the exact `k = 2` skeleton remainder is `±0.005` at the binding states (Stage 3, Δ_g = 100). The `k = 3` remainder is not computed.
* LR vs HYB at the large-`b̂` edge: `≤ 0.012`.
* All minima exceed their budgets: δ = 0.05 by ≈ 0.04, δ = 0.2 by ≈ 0.19.
* **Excluded from the budget:** graviton-specific beyond-skeleton corrections (§2.3, §4).

### 1.6 Flat-space dips outside the Stage-0 LP grid [NE]
* **What.** Both flat functionals have shallow minima at small `u`, `b̂ ≈ 55–85` (`stage6_scan_ext` docstring):
  * δ = 0.05: +0.025 at `m = 90.6Λ`, `+−`, J = 2689; +0.039 at `130Λ`, J = 4274;
  * δ = 0.2: +0.185 at `90.6Λ`, `++`, J = 2980 (this is Stage 3's dense-check minimum).
* **Why the LP missed them.** The scaling function `2∫ψ₂J₀(b̂p)` is tiny there, so the O(u) finite-J terms decide the sign. The Stage-0 LP grid
  stopped at J = 2500 (`b̂ ≤ 55` at these `u`), and its asymptotic rows start at `b̂ = 150`.
* **AdS lift.** In AdS the curvature at `β = b̂/Λ ≈ 0.06–0.9` lifts them: at Δ_g = 1000 to +0.087 (δ = 0.05) and +0.2185 (δ = 0.2).
* **Consequences.**
  * They do not bind at the three Δ_g.
  * At larger Δ_g the δ = 0.05 minimum would approach the flat +0.025.
  * A future LP should add rows with J up to `~100 m/Λ` at `u ≲ 10⁻⁴` [I, cheap].
  * HYB and LR agree there to `≤ 0.002` at J > 3000. LR is adequate there once the AdS curvature is included; the flat dip is a finite-`u`
    effect that LR also carries at these J.

## 2. Gauge robustness

### 2.1 Gauge families (`stage6_gauge.py` → `cache/stage6_gauge.txt`)
* **`Ĉ₅`** (3-shift, complex `c`). Gauges:
  * Stage 5's G0 (min-norm), G1 (`1+(ν/5)²`), G2 (Im×10), G3 (Re×10);
  * new: G4 smooth (min `‖Y′‖² + 10⁻²‖Y‖²`), G5 amplification-L2 (min `∫(|c|Σ_j|ĥ_j|)²`), G6 amplification-L1 (min `∫|c|Σ_j|ĥ_j|`, LP over the
    1002-dim null space).
* **`Ĉ₄`** (4-shift, value and derivative). Gauges:
  * H0 min-norm, tol `10⁻¹¹` (Stage 4; = amplification-L2, since the amplification density of χ is `|Y|`);
  * H1 tol `10⁻⁶`;
  * H2 `1+(ν/5)²`;
  * H3 smooth;
  * H4 amplification-L1.
* **Checks.** All gauges solve the Abel equation: residual `≤ 4·10⁻¹⁴` for `Ĉ₅` and `≤ 8·10⁻¹³` for `Ĉ₄`, with `W(t ≤ 6) = Φ ± 10⁻⁸`. The exception is
  H1, which is **not grid-stable**: valid at `dν = 0.04`, residual `2·10⁻⁷` and `W(6)` off by `10⁻²` at `dν = 0.02`.
* **Reproduction.** Stage 5's G0–G3 numbers are reproduced: `∫|c| = 1.150 / 1.112 / 1.139 / 3.401`, distances 0 / 0.17 / 0.29 / 2.57.
* **Amplification `A₀ = ∫|c|Σ_j|ĥ_j|`:**
  * `Ĉ₅`: 1.62 (L1) … 1.80 (G0) … 7.66 (G3);
  * `Ĉ₄`: 1.78–2.15 (value), 0.16–0.17 (derivative).
* **Null-space dimension.** `Ĉ₄` has 2 exact null directions plus one near-null (s.v. `2·10⁻¹¹`; Stage 4 counted it as a third). `Ĉ₅` has
  493 / 1002 null directions at `N = 512 / 1024`.
* **Covariant test [DV].** For the flat scalar-skeleton function `g(μ²) = (2−y)(1−y)⁻³P_J(1−2y)`, `y = μ²/m²` (all orders in `1/m`, entire up
  to `|μ| = m`), the quantity `Re Σ_j∫cĥ_j g((ν−ij)²)` equals `g(0)` (value, `Ĉ₅`) or `g′(0)` (derivative) **in every valid gauge**. The worst
  spread over gauges:
  * `3·10⁻¹⁵` (`Ĉ₅`);
  * `1.4·10⁻¹³` (`Ĉ₄` value);
  * `7·10⁻¹²` (`Ĉ₄` derivative).

  This holds at five states: `+−` J = 5 at `1.2Λ`; J = 4 at the gap; J = 8 at `1.1Λ`; J = 24 at `1.3Λ`; J = 32 at `3Λ` (Δ_g = 100).
* **Consequence [D + DV].** Every model with a covariant raw image (LR exactly, Stages 4/5 [P]; HYB; EXA, whose per-harmonic skeleton is
  analytic for `|μ| < m`) gives **identical minima in all gauges**. The Stage 5 / Stage 4 gauge question lives entirely in the non-covariant
  beyond-LR part.

### 2.2 A candidate canonical gauge: shift-blind (SB) (`stage6_canon.py`) [DV], not adopted
* **Construction.** The first structure that distinguishes gauges is a term read with the unshifted argument, `ν² − (ν−ij)²`. Its image is
  `∝ s_a := μ₁ − M₁`, with:
  * `μ₁ = Re Σ_j∫cĥ_jν²` (unshifted);
  * `M₁ = Re Σ_j∫cĥ_j(ν−ij)²` (shifted; 0 for value and `Ĉ₅`, 1 for the derivative).

  The SB gauge is min-norm subject to the Abel equation **and** `s_a = 0`. This is the "vanishing subleading (shift-violating) image" condition.
* **Results** (refinement-stable: `Ĉ₄` `5·10⁻⁵` and `1.2·10⁻⁴` for `dν` .04 → .02; `Ĉ₅` `9·10⁻³` for .08 → .04):

  | piece | `A₀` before → after | `s_b` before → after | covariant test |
  |---|---|---|---|
  | `Ĉ₄` value | 2.00 → 2.31 | −2.6 → **+20.8** | passes (`≤ 10⁻¹²`) |
  | `Ĉ₄` derivative | 0.17 → 0.20 | −0.17 → **−1.17** | passes (`≤ 10⁻¹⁰`) |
  | `Ĉ₅` | 1.80 → 2.05 | +2.8 → **+18.1** | passes (`≤ 10⁻¹⁴`) |

  Here `s_b` is the j-weighted probe, which SB does not constrain.
* **Why not adopted.**
  1. It does not remove sensitivity, it moves it: `s_b` grows ×7 to ×8, and the second structure then binds (§2.3).
  2. For `Ĉ₄` the condition is met only through the near-null direction, whose large-`t` accuracy is poor (`|W − 1| = 1.7·10⁻³` at t = 8, vs
     `6·10⁻⁶` for H0).
  3. With 2–3 null directions, `Ĉ₄` can cancel at most 2–3 non-covariant structures, while the general O(1/m²) class is larger.
* **Conclusion: no canonical gauge emerges cleanly; I report the worst case over gauges.** The SB weights are saved in
  `functionals_q6_canonical.json` for a later decision.

### 2.3 Minima under non-covariant stress, per gauge (`stage6_stress.py` → `cache/stage6_stress.txt`) [M]
* **Allowed size of the stress.**
  * LR is exactly covariant at every order in β (Stages 4/5, [P]). So a non-covariant part can only be beyond LR: O(1/m²) at fixed β, of the size
    of the **finite-J constants**. For scale, the scalar skeleton's is 5 and the graviton normal forms have `c₄ = 1` (`++`) and 41 (`+−`).
  * The unit-strength "no shift" model (full `g(ν²)` for every harmonic, first version in the cache log) violates LR at O(β²). It is therefore
    excluded and not used.
* **Stress structures**, with per-harmonic image `F₀ = 2`:
  * (a) `dF_j = F₀K[ν² − (ν−ij)²]/m²`;
  * (b) `dF_j = F₀K(j/j_max)²ν²/m²`.
* **Shifts in the heavy action:**
  * `δB₄ = 2Ku³s^{4v}/Λ²`;
  * `δ(∂B₄) = 2Ku³s^{4d}`: the derivative functional reads `Λ² ×` the image, so this one is **not** `1/Λ²`-suppressed;
  * `δB₅ = σKu⁴s⁵/Λ²`.
* **States:** 14 000 cases per (δ, Δ_g):
  * exact states;
  * gap states J ≤ 60;
  * region-A states with `m ≤ 30Λ`, J ≤ 60, `b̂ < 30`.

| coefficient (a, b) | H0 | H2 | H3 | H4 | SB |
|---|---|---|---|---|---|
| `s^{4d}` (derivative) | (−0.904, +0.369) | (−0.828, +0.038) | (−0.904, +0.385) | (−0.889, −0.006) | (0, −1.168) |
| `s^{4v}` (value) | (+5.57, +10.5) | (+6.71, +3.78) | (+5.62, +9.42) | (+5.85, +3.25) | (0, +20.8) |

`Ĉ₅` coefficients `s⁵ = (+5.2 … +8.0, +2.8 … +4.0)`, and (+53, −81) for G3.

**Tolerated strength** `K_tol` (largest `|K|` with all states ≥ 0; worst sign), Δ_g = 100 / 300 / 1000:

| | (a) shift-blind | (b) j-weighted | best gauge for (a) / (b) |
|---|---|---|---|
| δ = 0.05, H0 (Stage-4 gauge) | **0.64** (`K < 0`) | 1.55 (`K > 0`) | SB: ∞ / 0.50; H4: 0.65 / 102–139 (H1, not grid-stable: 1.0 / –) |
| δ = 0.05, worst over gauges | 0.64 | 0.50 (SB) | |
| δ = 0.2, H0 | **1.97** | 4.8 | SB: ∞ / 1.53; H4: 2.0 / 312–366 |
| δ = 0.2, worst over gauges | 1.96 | 1.52 (SB) | |

* The values are the same at all three Δ_g to 2 % (the derivative stress is Λ-independent).
* The binding state is always the **`++`, J = 0 gap state** (`u = 1`, where `u³` is maximal). The second most sensitive is `+−`, J = 6 at `1.475Λ`.
* Min `FF/N` at `|K| = 5`, over all gauges: δ = 0.05 from **−0.454** to +0.048; δ = 0.2 from **−0.458** to +0.198.
* **Gauge spread of the minimum at fixed `|K| = 5`: about 0.5 (both δ).** At `K = 0` the spread is 0.
* **The `Ĉ₅` gauge changes the minima by `≤ 0.003`** even at `|K| = 5` (G3 the extreme). Stage 5's 17–29 % gauge dependence of `Ĉ₅` is immaterial
  for positivity [M].

**Reading [M].**
* The minima are gauge-independent unless the beyond-LR forward images contain non-covariant O(1/m²) finite-J terms.
* If such terms exist with `|K| ≳ 0.5–0.6` (δ = 0.05) or `≳ 1.5–2` (δ = 0.2), some gauge of the `Ĉ₄` derivative piece fails at the J = 0 gap state.
* The natural scale is `|K| ~ 1–5`, since `c₄^{++} = 1` itself is a finite-J constant. So **gauge robustness is not established for either δ**.
  δ = 0.2 has about 3× the tolerance.
* Equivalent statement: Stages 3/12/6 assume that the de-shifted `∂Ĉ₄` reproduces the flat normal-form constants `c₄_h` covariantly (`K = 0`).
  That assumption is the open input, and it is the same [I] as the graviton-specific near-gap correction (Stage 4 §5).

### 2.4 Large-`t` precision of the stored weights (`stage6_wcheck.py` → `cache/stage6_wcheck.txt`) [NE]

| `|W/Φ − 1|` at t = | 6 | 8 | 10 | 12 | 14 |
|---|---|---|---|---|---|
| `Ĉ₄` value, H0 (`dν = .02`) | `3·10⁻⁹` | `6·10⁻⁶` | `7.5·10⁻⁴` | **0.41** | 38 |
| `Ĉ₄` derivative, H0 | `8·10⁻¹¹` | `2·10⁻⁸` | `5·10⁻⁷` | `8·10⁻⁴` | 0.11 |
| `Ĉ₅` G0 (`dν = .04`) | `10⁻¹¹` | `5·10⁻¹⁰` | `10⁻⁷` | `3.5·10⁻⁶` | `3·10⁻⁴` |

* **Cause.** A double-precision residual `ε` of the `cosh(jt)`-normalised rows becomes `|W − Φ| ~ ε cosh jt`.
* **Consequence.** The heavy image at β uses `W` on `[0, β]`, so the stored `Ĉ₄` value weight realises `P_{−½}` only for β ≲ 10. That is
  `J ≲ m·e⁵`, i.e. about 150m.
* Stage 4/5's β ≤ 30 positivity, and the validity counting, hold for the exact (continuum) weight. The Abel equation is exact there, and Stage 4/5's
  grades stand for that object.
* Realising the functional numerically beyond β ≈ 10 needs an extended-precision (or analytic large-`t`) solution [I].
* This was not visible in Stage 4/5's checks, which stopped at t = 8.

## 3. Verdict and constants

* **Heavy positivity over the full physical region [M]:** δ = 0.05 min +0.0478 / +0.0497 / +0.0491; δ = 0.2 min +0.1985 / +0.1977 / +0.1970.
  * This holds in the covariant (LR / HYB / EXA) models with both forward pieces de-shifted, at Δ_g = 100 / 300 / 1000.
  * **Both δ are supported in-model at all three Δ_g.** Minima are within 1.5 % of δ.
* **Not supported robustly [I]:** neither δ against non-covariant beyond-LR finite-J terms of natural size (§2.3).
  * δ = 0.2 tolerates `|K| ≈ 1.5–2`; δ = 0.05 only `≈ 0.5–0.6`.
  * Together with Stage 4's crude near-gap estimate (`0.007–0.10N` at Δ_g = 100), this makes **δ = 0.2 the value to quote** and δ = 0.05 the
    in-model value.
* **LP re-optimisation: not needed.** No newly covered state is negative. The only sub-δ minima are the flat off-grid dips (§1.6) and the R3 state
  (`+−`, J = 26, `6.27Λ`) at `0.985δ`. Both are positive and outside any binding role for the bound.
  * Adding them as LP rows would change `A` by an amount of order the 1.5 % shortfall [not run].
  * The cheap improvement for later is LP rows at J up to `100m/Λ` for `u ≲ 10⁻⁴` (§1.6).
* **Constants (unchanged; the functionals are unchanged).**

  | | `A` | `B_eff` (Δ_g = 100 / 300 / 1000) | `κ` | ABJM `2ζ(3)` as fraction of bound |
  |---|---|---|---|---|
  | **δ = 0.2 (quoted)** | **17.241** | **−19.038 / −19.220 / −19.265** | **0.110** | **4.5 / 3.3 / 2.6 %** |
  | δ = 0.05 (in-model) | 11.448 | −12.804 / −12.938 / −12.970 | 0.0686 | 6.8 / 5.0 / 3.9 % |

  * `A` is exact given the pole; `B_eff = 2L_AdS − A log(2Λ/3)` [DV]; `κ` is [M] with `b₀` [I].
  * The lower-bound slope `≥ 0.2519` is unchanged [M].
  * ABJM is inside, far from both boundaries. This is a weak, log-dominated test [D given Q; bound M]. A1 excludes ABJM's light KK spectrum, so
    this is a consistency check, not an application.

## 4. Task 4 (graviton-specific ν = 0 correction at J = 4–8): not attempted [I]
No engine evaluates the graviton (non-skeleton) forward image at `ν = 0` beyond leading Regge. §2.3 sharpens what is needed: the shift structure of
the O(1/m²) finite-J terms of the de-shifted `∂Ĉ₄` at the `++`, J = 0 gap state and the `+−`, J = 5–6 states near `1.2–1.5Λ`. Equivalently, `K` in
§2.3: `|K| < 0.5` would establish δ = 0.05, and `|K| < 1.5` would establish δ = 0.2.

## 5. Remaining open items [I]
1. The non-covariant beyond-LR finite-J terms of the de-shifted `∂Ĉ₄` (strength `K`, §2.3). They decide gauge robustness and δ. This is the same
   item as the graviton-specific near-gap correction.
2. An extended-precision de-shift solve so that the stored weights realise the pure image beyond β ≈ 10 (`Ĉ₄`) / 14 (`Ĉ₅`) (§2.4).
3. A canonical gauge: SB is a candidate (saved), not adopted (§2.2).
4. Flat LP rows for J > 2500 at `u ≲ 10⁻⁴` (dips +0.025 / +0.185 flat, §1.6). Relevant for Δ_g ≫ 1000.
5. Exact skeleton at J ≥ 40 (`mack.py` Chebyshev overflow); the `k = 3` skeleton remainder; the parity-odd block beyond leading order.
6. Inherited from Stages 0–5: `b₀`, `b₅`, `c₁`, spinning (F) and (P), A1, A2, K1′/K2/K3 and KKSDZ's loophole (validity), microcausality of tree
   correlators, the odd-`k` CHMRS continuation.

## 6. Graded claims

| claim | grade |
|---|---|
| β depends on `b̂` only (`cosh β = 1 + b̂²/2Λ²`); HYB = FLAT_EXACT + (LR − LR_flat) reproduces Stage 3 EXA to `3·10⁻⁹` with exact forward input | [P] / [DV] |
| HYB vs EXA `≤ 0.005` (Stage-3 states), `≤ 5·10⁻⁵` (39 new exact states, J ≤ 32, `m ≤ 20Λ`) | [DV] (moments) + [M] |
| HYB vs LR in the overlap: `≤ 0.018` (`b̂ ≥ 20`, `m ≥ 4Λ`); LR invalid near the gap even at `b̂ ∈ [20, 30]` (up to 0.6) | [M] |
| full-region minima δ = 0.05: +0.0478 / +0.0497 / +0.0491; δ = 0.2: +0.1985 / +0.1977 / +0.1970; no negative state | [M] |
| flat Stage-0 dips +0.025 / +0.185 at `b̂ ≈ 55–85`, J ≈ 2700–4300 (outside the LP grid); AdS-lifted at Δ_g ≤ 1000 | [NE] / [M] |
| all valid gauges map covariant images (incl. all-orders-in-1/m scalar skeleton) identically (`≤ 10⁻¹¹`) | [D] (Taylor / Abel) + [DV] |
| SB gauge exists, refinement-stable; trades `s_a` for `s_b`; not adopted | [DV] / judgement |
| `Ĉ₅` gauge immaterial for positivity (`≤ 0.003`); `Ĉ₄`-derivative gauge spread ≈ 0.5 at `|K| = 5`; `K_tol` 0.5–0.64 (δ = 0.05), 1.5–2.0 (δ = 0.2) | [M] (stress model) |
| stored weights realise `W = Φ` only to t ≈ 10 (`Ĉ₄` value), ≈ 14 (`Ĉ₅`) | [NE] |
| constants `A`, `B_eff`, `κ` and ABJM placement unchanged; δ = 0.2 recommended | [M] + [D given Q]; robustness [I] |

## 7. Files (`ads4/code/resultq/`)

| file | content | run time |
|---|---|---|
| `stage6_model.py` | HYB / LR / EXA heavy actions; `b̂`-tables (cached `cache/stage6_tab_*`); flat tables; gap states; exact-moment loader | (library) |
| `stage6_scan.py` | region A (J ≤ 3000, 459 `m`-values, gap states J ≤ 300), Stage-5-grid region B, EXA, overlap statistics, scaling limit → `cache/stage6_scan_{0.05,0.2}.{txt,pkl}`, `stage6_scan.log` | 25 min (tables) + 2 min |
| `stage6_scan_ext.py` | region A′ (J ≤ 40 000, `b̂ ≤ 200`, `m ∈ [20, 600]Λ`, HYB) and dense LR region B′ → `cache/stage6_scan_ext_{0.05,0.2}.txt` | 20 min |
| `stage6_fwexact.py` | 39 new exact CHMRS states (Δ_g = 100) → `cache/stage6_fw_100_Df3_w0.pkl` | 40 min |
| `stage6_gauge.py` | 7 `Ĉ₅` + 5 `Ĉ₄` (value, derivative) gauges, probes, covariant test → `cache/stage6_gauge.{txt,pkl}` | 6.5 min |
| `stage6_canon.py` | SB gauge, checks; writes **`functionals_q6_canonical.json`** (candidate, not adopted) → `cache/stage6_canon.{txt,pkl}` | 15 s |
| `stage6_stress.py` | K-stress per gauge, `K_tol`, minima at `|K| = 5` → `cache/stage6_stress.{txt,pkl}` (`cache/stage6_stress.log` holds the discarded unit-strength first version) | 15 s |
| `stage6_wcheck.py` | large-`t` accuracy of the stored weights → `cache/stage6_wcheck.txt` | 10 s |
