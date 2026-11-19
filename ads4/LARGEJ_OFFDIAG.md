# The large-`J` graviton off-diagonal of `Ω^{A′}`: a leading-Regge effect, and what it means for the route-A′ family

*19 Nov 2026, large-J agent. Unrefereed, not committed. Nothing in `RESULT_D.md`, `STATUS.md`, `H0_SPINNING_FUNCTIONAL.md`,
`H7_I14_LIGHTSIDE.md`, `S1_S4_FIXES.md`, the review or `n3_*` was edited. New code only: `ads4/code/spinfunc/r*_*.py` (+ `r1_queue{A,B}.sh`);
caches `spinfunc/cache/r*` (git-ignored). Project modules are imported read-only (cache names redirected by monkey-patching, as q2/q4 did).*

Grades: **[P]** proved, **[D]** derived, **[D-sketch]**, **[NE]** numerical evidence, **[M]** model, **[I]** open/assumed.
Conventions as `S1_S4_FIXES.md`: `Λ = √((Δ_g−2)(Δ_g−1))`, `b̂ = (2J+1)Λ/m`, `β = arccosh(1 + (2J+1)²/(2m²))` (the state's AdS impact
parameter, `b̂ ≈ Λβ` at small `β`), `x = ν/m`, `p = ν/Λ`; `Γ_off` in l9 units, so the `F₊₋` density per unit `ψ₂`-weight is `eo2·Γ_off/m`
(`eo2 = √(env₊₊env₊₋) ≈ 2`).

---

## 0. Summary

1. **No conflict with the flat limit, and no bug — but the failure is genuine** (verdict on the three options of the task):
   * **(a) is half right.** The fit "`ratio ∝ J^{3.5}/m²`" was made at `m ≈ Λ`; rewriting it as `(b̂/2)^{3.5} m^{1.5}` at fixed `b̂` is not
     legitimate. What the data actually say (§2) is that the real off-diagonal is a **leading-Regge object**:
     **`Γ_off/m = β² g(νβ)`** with one universal function `g` (peak `g ≈ −1.3` at `t = νβ ≈ 3.7`, zeros at `t ≈ 1.7, 6.5`), the same on 26
     rows (`J = 6–32`, `Δ = 336–1600`, `β = 0.008–0.11`) to `±0.1` [NE]. Hence at fixed `b̂` and fixed `Δ_g` the off-diagonal/base ratio is
     **independent of `m`** (it neither grows like `m^{1.5}` nor decays like `1/m²`), and at fixed `b̂`, `m/Λ` it **vanishes like `Λ⁻²` in the
     flat limit** (`≈ 3.7 b̂³/Λ²` on the wedge) — consistent with NF02's flat certificate [D given the collapse]. Equivalently `c_off(x; J) = (2J+1)² g((2J+1)x)`, which is
     both S1_S4's "`c_off` is `m`-independent at fixed `x`" and its "`J^{3.5}`" growth.
   * **(b) in the weak sense.** NF02's flat base has no zero: in the near-gap wedge it is already the universal tail `1023.6/b̂³`
     (`= −A₂/(2b̂³)`) to `≤ 2 %` for `b̂ ≥ 25` [NE] — "small" only in absolute terms. The off-diagonal integrates to a nearly `b̂`-independent
     `|ΔF₊₋| ≈ 3.8·10³/Λ²` (asymptotically `2ψ₂′(0)|G₁|/Λ²`, `G₁ = ∫t g(t)dt ≈ −12.5`) [NE + M]. So the ratio is `≈ 3.7 b̂³/Λ²`, and
     it crosses 1 at **`b̂_c ≈ 30` (`Δ_g = 300`) and `≈ 65` (`Δ_g = 1000`), for every `m/Λ`** (`b̂_c ≈ 0.65Λ^{2/3}`) [M].
   * **(c) no bug found.** The 2×2 test is `F ⪰ 0` (≡ det ≥ 0 and trace ≥ 0, basis independent for a real congruence); the off-diagonal
     normalisation is the same congruence as the diagonal's (§1.3); reference ratios are real positive, the antisymmetric part is
     `≤ 3 %`, and the collapse of rows with independent basis signs and normalisations onto one `g` is itself a strong consistency check [NE].
2. **Fixed-`b̂` ladder at `Δ_g = 300`** (`b̂ ≈ 33`): `(336, 18) → (672, 36)` with direct RQMC values: `Γ_off` grows by **1.96** (leading-Regge
   prediction 1.94; "`O(1/m²)` at fixed `(ν,β)`" would give 0.50), so the `F₊₋` density and the off-diagonal/base ratio are `m`-independent:
   `|F₊₋|/√(NN) = 1.53 / 1.25 / 1.21` at `m/Λ = 1.12 / 2.24 / 4.49` (universal model, §2.4) [NE + M]. The K₂ ladders at fixed
   `J/m` (`(8,400) → (16,800) → (32,1600)` etc.) give `Γ_off ∝ m^{s}` with **`s = +0.98 … +1.10`** (C2), i.e. an `m`-independent
   off-diagonal density at fixed `(ν, β)` [NE]. `(1344, 72)` is out of reach: the K₂ engine overflows (`M₂ ≈ M₀`) and the direct method
   was not run at `J = 72` [—].
3. **Consequences** [M, with the universal `g`]: with NF02 the per-state 2×2 is negative on the whole band `b̂ ≳ 30` (`Δ_g = 300`) / `≳ 65`
   (`1000`), **at every `m/Λ` where the band is twist-allowed (`1.1Λ … 5Λ` sampled)** (e.g. `λ/N_s = −2.1 / −7.7 / −26` at `b̂ ≈ 42 / 60 / 86`, `Δ_g = 300`). So positivity
   from `m ≥ nΛ` does **not** rescue a bound: the failure is in `b̂`, not in `m`. The `16n⁸L_A(Λ/n)` illustrations of S1_S4 §2.4 are moot.
4. **Redesign within route A′** (`ψ₂, ψ₃` re-optimised in NF02's own LP family, flat constraints and margin of NF02 kept, `k(0) = 1`,
   `Ĉ₃|_E = 0` and the graviton pole untouched because only the weights change) [M]:
   * On the **computed states only** (`J ≤ 20` near gap + hybrids to `2.5Λ` at 300; the universal model up to `b̂ ≈ 65` at 1000) `n = 1` is
     restored cheaply: `A = 3301 / B_AdS = −8728` (300, margin 0.1) → bound `3.04·10⁻⁶` (NF02: `3.13·10⁻⁶`); `A = 3051 / −7026` (1000)
     → `2.86·10⁻⁸` (NF02: `2.89·10⁻⁸`). **But these functionals just move the crossing** (300: `b̂_c ≈ 30 → 43`).
   * Covering the **band** costs fast-growing `A`, then becomes **infeasible**: `Δ_g = 300`: `b̂ ≤ 36 / 60 / 75 / 90` → `A = 3131 / 5091 /
     1.1·10⁴ / 6.4·10⁴` (`d = 12`), bound `3.03 / 3.32 / 4.2 / 9.1 ·10⁻⁶`; `b̂ ≤ 90` infeasible at `d = 8`. `Δ_g = 1000`: `b̂ ≤ 65 / 120 / 160`
     → bound `2.86 / 3.25 / 5.8 ·10⁻⁸`; **`b̂ ≤ 200` infeasible even at `d = 12`**. The twist-allowed band reaches `b̂ ~ Λ` (300) – `2Λ`.
   * **Why (structural, [D-sketch] given the universal `g`):** in the curvature-dominated regime `Λ^{2/3} ≪ b̂ ≪ Λ` both the diagonal
     (`≈ A₁M₁(0)/Λ² = 2ψ₂′(0)/(24Λ²)`) and the off-diagonal (`≈ 2ψ₂′(0)G₁/Λ²`) are proportional to `ψ₂′(0)`, so `|X|/F → 24|G₁| ≈ 300`
     **for every weight pair** with `ψ₂′(0) ≠ 0` (the `ψ₃` terms are `u`-suppressed, ≤ 15 %). `ψ₂′(0) = 0` would remove `A = 16ψ₂′(0)` (the
     4D IR log) and is excluded by flat positivity in this family (the LP's minimal `A` at margin 0 is 875 > 0) [NE].
   * **A third sum rule** could only help if its own crystal off-diagonal has a `G₁` that can cancel `Ĉ₂`'s uniformly in `u`; route A′
     has no such partner (`Ĉ₃`'s is `u`-weighted, sign-unstable `G₁⁽³⁾ ∈ [−6.5, +4]`), and any `k ≥ 4` rule reintroduces contacts/forward
     terms (H3, H9b, H15) [I, not computed]. **A weight that suppresses `x ≈ 0.05–0.15` is not available**: weights are functions of `ν`,
     and the dangerous window is `p ≈ 3.7/b̂`, which sweeps all small `p` as `b̂` varies.
5. **Status.** Per-state positivity of `Ω^{A′}` fails on a Regge-invariant band `b̂ ∈ [≈0.65Λ^{2/3}, O(Λ)]` at every `Δ_g` [NE + M]; within
   route A′ no choice of `(ψ₂, ψ₃)` repairs it [M + D-sketch], up to the validity of the small-`β` form of `g` at `β ≈ 0.2–0.3` (the infeasibility
   at `Δ_g = 300/1000` sits there; for `Δ_g ≳ 5·10³` it sits inside the validated range `β ≤ 0.11`). **In model M (universal `g`, diagonal graviton corrections and scalar remainder dropped beyond the computed
   states), Result D cannot be obtained from `Ω^{A′}` — NF02 or re-optimised — with per-state positivity at any `n`** [M]; the decisive open
   input is the off-diagonal at `β ≈ 0.2–0.3` (`ν ≈ 12–20`).

---

## 1. Checking the S1_S4 result

### 1.1 What "fixed `b̂`" means here [D]
`b̂ = (2J+1)Λ/m` is the flat-space impact parameter in units of `1/M` (`M ↔ Λ`). The flat limit is `Δ_g → ∞` at **fixed `J` and `m/Λ`** (so
fixed `b̂`); the leading-Regge limit is `m → ∞` at fixed `(ν, β)`, i.e. **fixed `J/m` at fixed `Δ_g`** (fixed `b̂` too, `m/Λ → ∞`). The task's
ladder `(336, 18), (672, 36), (1344, 72)` is the second: `b̂ = 33.07 / 32.55 / 32.29`, `m/Λ = 1.12 / 2.24 / 4.49`, `β = 0.111 / 0.109 / 0.108`.
S1_S4's two scalings are consistent with both limits once written in these variables:
* at fixed `(x, J)` (rows 400–1600): `c_off = mΓ_off` constant ⇒ `F₊₋` density `∝ 1/m²` — the flat direction, `∝ 1/Λ²` at fixed `m/Λ`;
* at `m ≈ Λ`, varying `J`: ratio `∝ J^{3.5}` — this is a `b̂`-dependence (`b̂ ≈ 2J`), not a `J`-at-fixed-`b̂` one.

### 1.2 The flat base in the wedge (option b) (`r4`/inline, NF02 weights) [NE]
`2∫ψ₂J₀(b̂p)dp` (the scaling-limit flat action; `B₂` near gap agrees to `≥ 99 %`, S1_S4 §2.4): `b̂³·F = 962 / 985 / 1046 / 1033 / 1018 / 1026 /
1023 / 1024` at `b̂ = 20 / 25 / 30 / 35 / 40 / 50 / 60 ≥ 80`. **No zero, no anomalous dip**: the wedge is already on the universal tail
`−A₂/(2b̂³) = 1023.5/b̂³` of ROUTEA_LARGE_B. `N_s ≈ 0.12 / 0.08 / 0.035 / 0.019` (J = 10/12/16/20) is this tail.

### 1.3 Normalisation and criterion (option c) (`r7_checks.py`, `cache/r7_checks.txt`) [NE + D]
* **Criterion.** `λ_min(N^{−½}FN^{−½}) ≥ 0 ⇔ F ⪰ 0 ⇔ det F ≥ 0 and tr F ≥ 0`; for real OPE vectors only the symmetric part of `Re F`
  enters. This is what q3 tests. [D]
* **Congruence.** Diagonal: `M̃_jj s_j/(M̃^{ref}_jj a_j)`; off-diagonal: `M̃_01 √(s₊s₋)/√(M̃^{ref}_00 M̃^{ref}_11 a₊a₋)` (via
  `sq = √|M̃^{ref}_00/M̃^{ref}_11|`), i.e. `D F D` with one real `D` per state, common to `Ĉ₂` and `Ĉ₃` (the `Ĉ₃` off-diagonal uses the `Ĉ₂`
  reference in `sq`, which is the correct common `D`). [D, code read]
* Checks: `M̃^{ref}_{++}/M̃^{ref}_{+−}` (order 0, `Δ_ref = 2·10⁴`) with positive real part on every node of every row used (`Re ∈ [36, 660]`;
  `|Im/Re| ≤ 5·10⁻⁴` on `Δ ≤ 800`, `≤ 6·10⁻²` on `Δ = 1600` except one large-`ν` node of `J = 20` with 0.47), so `|·|` in `sq` hides no
  sign; `|M₀₁ − M₁₀|/|M₀₁ + M₁₀| ≤ 3·10⁻²`; `a₊₊, a₊₋ > 0` on every row (`J = 8–32`). The order-0 matrix is diagonal in the chosen real
  basis to `≤ 2·10⁻²` relative at `Δ = 1600` (larger, up to 0.24, only at nodes next to zeros of `D_h` on the `Δ = 336/400` rows); any such
  residual is part of the physical off-diagonal in that basis and is basis-consistent. [NE]
* **Row-sign alignment (S1_S4 §2.2a)** only matters for interpolation between rows of one `J`; single-row tests (direct `Δ = 336` row at
  the binding state's own mass) are unaffected, and the per-state eigenvalue depends on `|F₊₋|` only.
* The strongest check is §2.2: rows with independent basis signs, `a_h` fits and `c0` constants collapse onto one function of `νβ`.

## 2. The off-diagonal is a leading-Regge object

### 2.1 Fixed-`J/m` ladders from the existing K₂ rows (`r2_ladder.py`, `cache/r2_ladder.txt`) [NE]
`Γ_off(ν)` at common `ν`; `Γ_off ∝ m^s` at fixed `(ν, β)` (the density `∝ m^{s−1}`):

| ladder (`J/m`) | `s` (rms over the common ν-window), C2 | C3H |
|---|---|---|
| (6,400)→(12,800)→(24,1600) (0.015) | +0.34 → **+0.71** | +0.49 → +0.99 |
| (8,400)→(16,800)→(32,1600) (0.020) | +0.84 → **+0.98** | +2.7 → +1.5 |
| (10,400)→(20,800) (0.025) | **+1.10** | +2.0 |
| (14,800)→(28,1600) (0.0175) | **+1.09** | +0.44 |

The C2 profiles coincide in shape, with zeros at the same `ν` (e.g. `ν ≈ 160` on the 0.020 ladder). `s = 1` means an `m`-independent density:
the real off-diagonal does **not** vanish in the Regge limit. (C3H is noisier and enters with weight `uψ₃`, `≲ 15 %`.)

### 2.2 Universal function (`r4_collapse.py`, `cache/r4_collapse.txt`, `cache/r4_g.pkl`) [NE]
`f/β² = Γ_off/(mβ²)` against `t = νβ` for all 26 rows with `Δ ≥ 336`:

| `t` | 1 | 1.5 | 2 | 2.5 | 3 | 3.5 | 4 | 4.5 | 5 | 6 | 7 | 8 | 10 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| median `g(t)` | 0.80 | 0.24 | −0.30 | −0.78 | −1.10 | −1.28 | −1.23 | −1.09 | −0.85 | −0.25 | 0.11 | 0.08 | −0.11 |
| spread | 0.03 | 0.09 | 0.07 | 0.09 | 0.09 | 0.08 | 0.07 | 0.09 | 0.11 | 0.11 | 0.06 | 0.21 | 0.42 |

Rows span `β = 0.008–0.111`, `m = 334–1598`, `J = 6–32`, `J²/m ≤ 1.0`; e.g. `(8,400)`, `(16,800)`, `(32,1600)` and `(18,336)` agree to
`≈ 0.1` at the peak. So **`F₊₋` density `= eo2 β² g(νβ) ψ₂ + (k = 3)`**, an `O(β²) = O(b̂²/Λ²)` AdS-curvature effect at leading Regge order
(the same order as the positive scalar curvature `A₁M₁(β)/Λ²`). The `ν`-window of S1_S4 (`x ≈ 0.05–0.15` at `m ≈ Λ`, `J = 16–20`) is the peak
`t ≈ 2.5–5`. Note `Γ_off = m β² g` gives `c_off(x; J) = (2J+1)² g((2J+1)x)` at small `x`.

* **Regge structure [D-sketch].** The leading-Regge heavy images of the diagonal are `l`-isotypic (`P_l`, S1_S4 §1.2, H7 §18). The
  off-diagonal connects structures whose transverse helicities differ, so equivariance allows an associated-Legendre-type image vanishing
  like `β²` at small `β`; the measured `β²g(νβ)` is of that form. Not derived here [I].

### 2.3 Non-perturbative ladder point `(672, 36)` (`r1_run.py`, `r6_ladder_direct.py`, `cache/r6_ladder_direct.txt`) [NE]
Same chain as q4 (`n1_analysis.grid_direct`: K₂ helicity basis of the row, `a_h` fitted on the six direct nodes, `Re`, l9 units), on
`ν = 25.0, 33.4, 41.7, 50.1, 58.4, 66.8` at both states (`t = νβ = 2.7–7.4`). The `(672, 36)` points are new: K₂ (13 nodes `ν = 16.7 … 167`, C2 + C3H,
state + reference, ~155 s each; queue A was stopped before `ν = 184 … 267`, not needed here), reference basis, and direct RQMC C2 at the six nodes (state `≈ 285 s`, reference `≈ 75 s`; RQMC relative
errors `1.6·10⁻⁴ – 3.7·10⁻³`; the reference needed `point_g`, see §5). `J²/m = 0.97` and `1.94`.

| `ν` | 25.0 | 33.4 | 41.7 | 50.1 | 58.4 | 66.8 |
|---|---|---|---|---|---|---|
| `Γ_off` C2, `(336,18)` direct | −3.06 | −5.08 | −4.33 | −1.97 | 0.13 | 0.82 |
| `Γ_off` C2, `(672,36)` direct | −5.17 | −9.42 | −8.94 | −4.97 | −0.41 | 2.15 |
| `Γ_off/(mβ²)`, `(336,18)` direct | −0.75 | −1.24 | −1.06 | −0.48 | 0.03 | 0.20 |
| `Γ_off/(mβ²)`, `(672,36)` direct | −0.65 | −1.18 | −1.12 | −0.62 | −0.05 | 0.27 |
| median `g(νβ)` (§2.2, K₂ rows; read off at `t ≈ 2.75 … 7.4`) | ≈ −0.9 | ≈ −1.27 | ≈ −1.05 | ≈ −0.6 | ≈ −0.05 | ≈ 0.1 |
| `Γ_off` C2, `(672,36)` K₂ | −8.02 | −10.97 | −9.42 | −4.82 | −0.08 | 2.35 |

* **rms ratio `Γ_off(672,36)/Γ_off(336,18)` = 1.96 (direct) / 1.88 (K₂)**, against **1.94** predicted by the leading-Regge form `mβ²g(νβ)` and
  **0.50** for an off-diagonal that is `O(1/m²)` at fixed `(ν, β)`: `s = +0.97` (direct). The non-perturbative values collapse onto `g` to
  `±0.1`. **(a) in the sense "the off-diagonal/base ratio dies at fixed `b̂` as `m` grows" is excluded** [NE].
* `Ĉ₃` (K₂ only; its own C3H direct was not run): `Γ_off` equal at the two states to `≤ 0.5` (i.e. `s ≈ 0` for this kind on this ladder);
  it enters with `uψ₃`, `u = 0.80 → 0.20`, and is subdominant.
* The **diagonal** corrections are also not suppressed along this ladder (`Γ_hh/m²` at `(672,36)` ≈ 0.9 × that at `(336,18)`, direct), unlike
  the `β ≈ 0.02–0.03` K₂ ladders of §2.1 (`∝ m⁻¹`). They are `≤ 0.2 N_s` at the binding states (q3) and are dropped in the universal model
  of §2.4/§3 [M].
* `(1344, 72)`: not computed. The K₂ engine overflows at `J = 72` (`M₂ ≈ M₀`, `M₄ ≈ M₀` at the state; `l1_tps` overflow warnings,
  `cache/r1_pts_*_J72_k2.pkl`, one node) and a direct row (`≈ 20 min` per state point, untested contour at `J = 72`) did not fit the budget;
  the universal model gives `|F₊₋|/√(NN) = 1.21` there (§2.4).

### 2.4 Consistency with the near-gap numbers (`r5_universal.py check`) [NE + M]
Universal model: `F_hh` = flat + leading-Regge curvature (= ROUTEA_LARGE_B's `F_hyb`; scalar remainder and graviton diagonal corrections
dropped — `≤ 0.2 N_s` at the near-gap binding states (q3); they are `∝ m⁻¹` along the small-`β` K₂ ladders (r2) but not along the
`β ≈ 0.11` direct ladder (§2.3), so this is a model choice [M]), `F₊₋` from `g₂, g₃`. It reproduces
q3's near-gap off-diagonal ratios (`Δ_g = 300`: `|X|/√(NN) = 2.13` at `b̂ = 37`, `m = 1.1Λ`; q3: 2.15 at `J = 20`, `1.07Λ`).
First negative `b̂` per `m/Λ`:

| `m/Λ` | 1.1 | 1.2 | 1.3 | 1.5 | 2 | 3 | 5 |
|---|---|---|---|---|---|---|---|
| `Δ_g = 300` | 30 | 31 | 35 | 30 | 30.5 | 32 | 46 (grid) |
| `Δ_g = 1000` | 66 | 94 (grid) | 65 | 73 | 72.5 | 66 | 66 |

**Off-diagonal/base ratio at fixed `b̂` vs `m` (the task's ladder) and vs `Δ_g` (the flat limit)** (`r5_universal.py states`, NF02):

| `Δ_g` | state `(m/Λ, J)` | `b̂` | `β` | `|F₊₋|Λ²` | `|F₊₋|/√(N₊₊N₊₋)` | `λ/N_s` |
|---|---|---|---|---|---|---|
| 300 | (1.119, 18) = (336, 18) | 33.07 | 0.111 | 3959 | **1.53** | −0.53 |
| 300 | (2.243, 36) = (672, 36) | 32.55 | 0.109 | 3333 | **1.25** | −0.25 |
| 300 | (3.0, 50) | 33.67 | 0.113 | 3391 | **1.40** | −0.40 |
| 300 | (4.491, 72) = (1344, 72) | 32.29 | 0.108 | 3251 | **1.21** | −0.21 |
| 1000 | (1.119, 18) | 33.07 | 0.033 | 3959 | 0.137 | +0.86 |
| 3000 | (1.119, 18) | 33.07 | 0.011 | 3959 | 0.015 | +0.99 |
| 1000 | (1.04, 32) | 62.50 | 0.063 | 4714 | 1.11 | −0.11 |
| 1000 | (2.08, 64) | 62.02 | 0.062 | 4208 | 0.98 | +0.02 |
| 1000 | (4.16, 130) | 62.74 | 0.063 | 4181 | 1.01 | −0.01 |

So at fixed `b̂` the ratio is `m`-independent to `±15 %` (the residual is `eo2` and the `m/Λ`-dependence of the flat base), and at fixed
`(b̂, m/Λ)` it falls exactly as `Λ⁻²` (`|F₊₋|Λ²` constant). (q3 with the computed `J = 18` rows at `1.12Λ`: `λ/N_s = −0.44` (direct) /
`−0.55` (K₂), consistent with −0.53.)

`Δ_g = 300`, `β ≤ 0.3`: `λ/N_s = −2.1 / −7.7 / −26` at `b̂ = 42 / 60 / 86`. The absolute off-diagonal is nearly constant,
`|ΔF₊₋|Λ² ≈ 3.7–3.9·10³` for `b̂ = 37–86` (asymptote `2ψ₂′(0)|G₁| + uψ₃′(0)|G₁⁽³⁾| ≈ 4.6·10³ ± 0.3·10³`), against the curvature
`A₁M₁(0) ≈ 17.9` and the flat tail `1023.6Λ²/b̂³`.

## 3. Redesign within route A′ (`r3_redesign.py`, `r5_universal.py lp`) [M]

**Family and constraints.** NF02's own LP (`robust/nf_lp.py`, read-only): `ψ_k = (1−p)⁵p Σ_{i≤d} v_i L_i(2p−1)`, `k(0) = 1`, flat heavy states
with margin `δ = 0.2` (cutting planes on NF02's `(u, h, J)` grid), scaling limit, `b̂ → ∞` tail, scalar matter, Lemma-M proxy 120. With no AdS
rows it reproduces NF02 exactly (`A = 3022.34`, `B_AdS = −6552.66`). AdS rows: per state, `F(v)` is linear in the weights; 2×2 PSD
imposed by angle cuts at the minimal eigenvector (exact in the limit), margin `d_w` on `N_s` of NF02. Only the weights change, so the
graviton pole identity, `Ĉ₃|_E = 0`, (HN), the flat limit (= the new flat functional, positive by construction) and the definitions of
`A = 16ψ₂′(0)`, `B_AdS = 16L_fin − A log 1.5` are kept. Bound `= 240√(A log Λ + B_AdS)/Λ⁴`.

| AdS rows | `A` | `B_AdS` | bound 300 | bound 1000 | `n = 1` where |
|---|---|---|---|---|---|
| none (NF02) | 3022.3 | −6552.7 | 3.12·10⁻⁶ | 2.89·10⁻⁸ | — (fails, §2.4) |
| q3 computed states, 300 (J ≤ 20, 873 blocks, actual corrections), `d_w = 0 / 0.1` | 3269.7 / 3300.6 | −8587 / −8728 | 3.03 / 3.04·10⁻⁶ | — | those states; universal check: fails from `b̂ ≈ 43` |
| universal, 300, `β ≤ 0.12` (`b̂ ≤ 36`), `d_w = 0.1` | 3131.1 | −7802.7 | 3.03·10⁻⁶ | — | `b̂ ≤ 36` |
| universal, 300, `β ≤ 0.2` (`b̂ ≤ 60`), 0.1 | 5091.0 | −16982 | 3.32·10⁻⁶ | — | `b̂ ≤ 60` |
| universal, 300, `β ≤ 0.25` (`b̂ ≤ 75`), 0 | 11706 (d=8) / 10622 (d=12) | −45800 / −41321 | 4.37 / 4.19·10⁻⁶ | — | `b̂ ≤ 75` |
| universal, 300, `β ≤ 0.3` (`b̂ ≤ 90`), 0 | infeasible (d=8) / 63585 (d=12) | — / −2.7·10⁵ | — / 9.05·10⁻⁶ | — | `b̂ ≤ 90` |
| universal, 1000, `β ≤ 0.066` (`b̂ ≤ 66`, covers `J ≤ 34` near gap), 0.1 | 3051.0 | −7025.6 | — | 2.86·10⁻⁸ | `b̂ ≤ 66` |
| universal, 1000, `β ≤ 0.12` (`b̂ ≤ 120`), 0.1 | 5005.2 | −16509 | — | 3.25·10⁻⁸ | `b̂ ≤ 120` |
| universal, 1000, `β ≤ 0.16`, 0 (d=12) | 20725 | −85513 | — | 5.80·10⁻⁸ | `b̂ ≤ 160` |
| universal, 1000, `β ≤ 0.2` (`b̂ ≤ 200`), 0 | infeasible (d = 8 and d = 12) | | | | |

Robustness of the infeasibility (`d = 8`, margin 0) [M]: **starting positivity at `m ≥ 1.5Λ` or `m ≥ 2Λ` (i.e. `n = 1.5, 2`) does not help** —
infeasible at 300 (`β ≤ 0.3`) with only `m/Λ ∈ {1.5, 2, 3, 5}` and with only `{2, 3, 5}` in the AdS rows, and at 1000 (`β ≤ 0.2`) with only
`{1.5, 2, 3, 5}`. With `g` halved
(a crude allowance for the unvalidated `β ≥ 0.2` part): 300 / `β ≤ 0.3` becomes feasible at `A = 1.0·10⁴` (bound `4.1·10⁻⁶`), 1000 / `β ≤ 0.2`
stays infeasible.

(The q3-state redesign improves the 300 bound slightly because NF02's objective was `≈ A`, not the bound; with the bound as objective the
no-AdS optimum is lower still — the comparison that matters is feasibility, not the third digit.)

**Reading.** Positivity can be bought on any finite window of `b̂`, at a cost that explodes as the window approaches `b̂ ~ 3Λ^{2/3}`, and is
infeasible beyond. The twist-allowed band extends to `b̂ ≈ Λ` (and `≈ 2Λ` at larger `m/Λ`). Since the failure is `m`-independent at fixed `b̂`,
the window cannot be avoided by starting positivity at `m ≥ nΛ`. **No member of the route-A′ family gives a Result-D bound with per-state
positivity** [M; asymptotic argument D-sketch]. Caveat: the infeasibility at `Δ_g = 300/1000` uses `g` at `β = 0.2–0.3`, beyond the
validated `β ≤ 0.11` (the peak there is at `ν ≈ 15`, below the direct method's range; only the n3 contour reaches it) [I].

**Asymptotic argument [D-sketch, inputs NE].** For `b̂ ≫ t`-support of `g`, `ΔF₊₋ = (b̂/Λ)²∫ψ₂ eo2 g(b̂p)dp → 2ψ₂′(0)G₁/Λ²` (`G₁ = ∫t g dt
≈ −12.3 … −13.4` over `t ≤ 6–12`); diagonal `→ 1023.6/b̂³·(…) + A₁M₁(β)/Λ²` with `A₁ = 2ψ₂′(0) + σuψ₃′(0)` and `M₁(0) = 1/24`. Once
the curvature dominates (`b̂ ≳ 4Λ^{2/3}`), `|F₊₋|/F_hh → 24|G₁|·(1 + O(u ψ₃′/ψ₂′)) ≈ 300`, independent of the weights. Before that,
`|F₊₋|/F ≈ 2ψ₂′(0)|G₁|b̂³/(1023.6Λ²)` with the tail coefficient itself fixed by `ψ₂″(0)`-type data, so a redesign can only shift `b̂_c ∝ Λ^{2/3}`
by an `O(1)` factor, as the LP shows. Escape routes not excluded: `g` changing sign/decaying at `β ≳ 0.2` (then the obstruction moves to
larger `Δ_g` only); an extra sum rule whose off-diagonal cancels `Ĉ₂`'s (needs its crystal computation; `k ≥ 4` reintroduces H3/H9b/H15).

## 4. Graded summary

| statement | grade |
|---|---|
| Real helicity off-diagonal = `eo2 β² g(νβ)ψ₂ + (k = 3)` with universal `g` (26 K₂ rows, `β ≤ 0.11`, `J²/m ≤ 1`, ±0.1) | [NE] |
| Fixed-`b̂` ladder, direct RQMC: `Γ_off(672,36)/Γ_off(336,18) = 1.96` (Regge form 1.94; `O(1/m²)` form 0.50) | [NE] |
| ⇒ off-diagonal/base at fixed `b̂`: `m`-independent at fixed `Δ_g`; `→ 0` as `Λ⁻²` at fixed `b̂`, `m/Λ` (flat limit consistent) | [D given the collapse] |
| `J^{3.5}/m²` → "`(b̂/2)^{3.5}m^{1.5}` growth" | **withdrawn** (misreading of a fixed-`m ≈ Λ` fit) [D] |
| NF02 flat base in the wedge `= 1023.6/b̂³` (±2 %, `b̂ ≥ 25`), no zero | [NE] (+ form [D], ROUTEA_LARGE_B) |
| 2×2 criterion and off-diagonal normalisation consistent; no sign hidden by `sq`; antisymmetric part ≤ 3 % | [D] + [NE] |
| With NF02, `n = 1` fails on `b̂ ≳ 30` (300) / `≳ 65` (1000) at **every** `m/Λ ∈ [1.1, 5]` | [NE + M] (universal model reproduces q3) |
| Positivity from `m ≥ nΛ` cannot rescue the bound | [M] |
| Redesigned `(ψ₂, ψ₃)` restore `n = 1` on the computed states (300: `A = 3301`, bound `3.04·10⁻⁶`; 1000: `A = 3051`, `2.86·10⁻⁸`) | [M] |
| ... but only by moving the crossing; band coverage infeasible at `b̂ ≈ 90` (300, d = 8; d = 12 gives `9·10⁻⁶`) and `b̂ ≈ 200` (1000) | [M]; `β ≥ 0.2` input [I] |
| Weight-independent asymptotic obstruction `|F₊₋|/F → 24|G₁| ≈ 300` | [D-sketch] given the universal `g` |
| Third sum rule; leading-Regge derivation of `g` | [I], not computed |

## 5. Files and hand-over

**Code** (`ads4/code/spinfunc/`, new):
* `r1_run.py` (+ `r1_queueA.sh`, `r1_queueB.sh`): K₂ (`k2`, `k2c2`, `refbasis`) and direct (`direct`; `point_g` = `n1_run.point` without the
  stacked scalar row, whose relative exponent overflows at `Δ_ref = 2·10⁴` for `J = 36` and masked every sample) points for arbitrary
  `(Δ, J)` on a `ν`-grid; caches `cache/r1_pts_<kind>_J<J>_k2.pkl`, `cache/r1_shadow_k2.pkl`, `cache/r1_n1_<kind>_J<J>.pkl`.
* `r2_ladder.py`: fixed-`J/m` ladders from existing rows (also merges r1 caches into q3's loader; imported by r3–r7).
* `r4_collapse.py`: the universal `g(t)` (→ `cache/r4_g.pkl`).
* `r5_universal.py`: universal-`g` 2×2 model; modes `check`, `eval`, `lp` (env `R5_DG`, `R5_BMAX`, `R5_DW`, `R3_D`).
* `r3_redesign.py`: per-state linear maps from q3's model (actual computed corrections) + NF02's LP with AdS angle cuts.
* `r6_ladder_direct.py`: `(336,18)` vs `(672,36)` with direct values. `r7_checks.py`: normalisation checks.

**Restart / reproduce:** `bash r1_queueA.sh` (resumes the `(672,36)` K₂ grid), `bash r1_queueB.sh`; `python3 r2_ladder.py`;
`python3 r4_collapse.py`; `python3 r6_ladder_direct.py`; `R5_DG=300 R5_BMAX=0.3 python3 r5_universal.py check|lp` (+ `R5_STATES="1.1189:18,..." ... states`,
`R5_MFLS`, `R5_GSCALE`, `R3_D=12`); `R3_DG=300 R3_JEVEN=0,2,...,20 R3_DW=0.0,0.1 python3 r3_redesign.py` (~70 s). Each analysis run is < 2 min.

**Next steps (value order):** (1) the off-diagonal at `β ≈ 0.2–0.5` (n3 contour at `ν ≈ 10–20`, or a leading-Regge kernel computation of
the off-diagonal à la `q1`/`v3` with helicity-changing structures) — this decides whether the obstruction at `Δ_g = 300/1000` is
exactly as modelled; (2) derive `g` from the Regge kernel (should be a weight-2 associated-Legendre image); (3) if Result D is to be
kept, a functional outside route A′ whose `Ĉ₂`-type off-diagonal moment `G₁` is cancelled (or a non-per-state positivity argument that
uses the OPE structure of the `(++, +−)` couplings).
