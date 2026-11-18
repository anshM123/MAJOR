# H13′ for D-safe (NF02 / `Ω^{NF,SC}`) at large impact parameter: exact-skeleton scans and a large-`b̂` envelope

> **Superseded functional (18 Nov 2026, later).** This note was written for the **B′-uplift D-safe functional `Ω^{NF,SC}` (k = 2…7)**.
> That functional is **withdrawn** (`H0_SPINNING_FUNCTIONAL.md` §9–§11; box at the top of `RESULT_D.md`). D-safe is now route A′,
> `∫ψ₂Ĉ₂ + ∫ψ₃Ĉ₃`. H13′ for route A′ is in **`ROUTEA_LARGE_B.md`**.
>
> Still valid here, because it does not depend on the functional:
> * the exact scalar-skeleton moments `ρ_k` (route A′ reuses them);
> * the method fixes of §3.1, plus those of §7 below;
> * the certified positivity of `M₁(β)` (§7);
> * the derived identities of §4.1 for `M₁` and `M₂`.
>
> The positivity statements, the tables of `F/N_s` and the H13′ grade below refer to the withdrawn functional only.

Status: 18 Nov 2026, H13′ agent. **Unrefereed. Not committed.** This answers the H13′ row of `RESULT_D.md` §4 and step 3 of §6 for the
D-safe functional (NF02 weights, no forward terms, even kernel `K*`). Only new files were written: those under `ads4/code/exactskel/`
and this note.

Code: `ads4/code/exactskel/h13_*.py`, `mack_flint.py`, `mack_mat.py`, `t_h13_*.py`. Caches and logs are in `exactskel/cache/`
(git-ignored). Every number below is printed by one of:
* `h13_summary.py` → `cache/h13_summary.txt`;
* `h13_asym.py` → `cache/h13_asym.txt`;
* the `t_h13_*` checks → `cache/t_h13_*.txt` and `cache/t_h13_checks.log`.

Grades:
* **[certified]**: Arb/exact with a rigorous enclosure;
* **[D]**: derived;
* **[NE]**: floating-point or Arb numerics without a rigorous truncation bound;
* **[M]**: a model convention;
* **[I]**: inferred or assumed.

As everywhere in `EXACT_SKELETON_NEARGAP.md`, "exact" means the exact **scalar-skeleton** model. Nothing here touches H2/H4: the
graviton dictionary, and the transfer of the scalar remainder to gravitons and to `+−`.

## 0. Results

**Bottom line.** No non-positive D-safe state was found anywhere in `m ∈ [1.5, 10]Λ`, `b̂ ∈ [20, 200]` (`Δ_g = 300, 1000`) or
`b̂ ∈ [40, 200]` (`Δ_g = 100`), in either channel or skeleton. The minimum is `F/N_s ≈ 0.93`, at the corner `m = 1.5Λ`, `b̂ ≈ 20`, `+−`.
Above `b̂ ≈ 50` the action has a derived, positive asymptotic form. It is an AdS term `∝ Λ⁻²`, not the flat `1/b̂³` tail.

**What is proven, and what is not.**
* The positivity of the full action rests on two pieces: a dense scan of the hybrid action, and an exact skeleton remainder sampled up
  to `b̂ ≈ 100`.
* The remainder is `|S| ≤ 131 F_hyb/m²` on the sample, i.e. `≤ 0.6 %` of the action.
* For `b̂ ∈ (100, 200]`, and beyond, the remainder is **not** computed exactly. There it is a calibrated assumption [I], with a
  `≥ 150×` margin.

| item | result | grade |
|---|---|---|
| **1a. Hybrid action** `F_hyb` (flat NF02 + exact leading AdS Regge curvature), every integer `J`, 53–80 masses in `[1.5, 10]Λ` (+ `20, 50, 100Λ`), both channels | positive everywhere. Min `F_hyb/N_s` = **0.931** (`Δ_g` = 300, 1000; `+−`, `1.5Λ`, `J = 15`, `b̂ = 20.7`) and **0.941** (`Δ_g = 100`; `+−`, `1.68Λ`, `b̂ = 108`). On `b̂ ∈ [200, 1000]` at 5 masses: ≥ 0.899 | [NE] (ingredients exact) |
| **1b. Exact scalar-skeleton remainder** `S`, Arb, both skeletons: 56 reliable states (of 59), `J ≤ 120`, `b̂ ≤ 100.5`, `m ≤ 10Λ` | no non-positive record. Min `F/N_s` = **0.949** (300, 1000), **0.967** (100). `|S|/F_hyb ≤ 5.9·10⁻⁴` (300), `5.3·10⁻⁵` (1000), `3.2·10⁻³` (100). `|S|/F_hyb ≤ 131/m²` on all of them | [NE + M(odd k)] on the sample |
| **2. Large-`b̂` envelope** | `F_hyb ≃ A₁M₁(β)/Λ² + A₂M₂(β)/(2Λ³)`, with `A₁ ≥ 355 > 0`, `A₂ = −2047`, `M₂ = −(cosh β−1)^{−3/2}/(2√2) < 0`, `M₁(β) = ∫₀^∞2νP_{−½+iν}(cosh β)dν/(e^{2πν}+1) > 0`. Matches `F_hyb` to `[0.9998, 1.003]` for `b̂ ≥ 100` (`[0.9975, 1.012]` for `b̂ ≥ 50`) on every record | expansion and signs [D]; `M₁ > 0` and the remainder factor [NE] |
| **2′. Remainder envelope** | `F ≥ F_hyb(1 − C/m²)`, `C = 131` (calibrated on the sample) ⇒ positive with factor ≥ 150 to spare at `m ≥ 1.5Λ`, `Δ_g ≥ 100` | [I] outside the sample |
| **3. `Δ_g = 100`, `b̂ ∈ [40, 200]`** | hybrid min 0.941. Exact (`b̂ ≤ 70.5`) min 0.967, `|S|/F ≤ 3.2·10⁻³`. Consistent with `Δ_g = 300, 1000` and with the earlier `J ≤ 28` sample | [NE] |
| **Method** | (a) exact Mack polynomials ~200–1000× faster (`mack_flint`, `mack_mat`, bit-identical). (b) **Bug-class found:** the `n`-sum tail of `skelexact3.Laplace3` (double-precision Gauss weights; Euler–Maclaurin at `J` not `≪ τ`) gives `O(1)` errors at `J ≳ 40`. It does not affect the existing `J ≤ 28` results (`≤ 2·10⁻⁹` change at `J = 20`). Fixed with Arb weights plus a 3-configuration spread estimate; 3 states remain unreliable and are excluded | [D]/[NE] |

**Proposed H13′ entry (D-safe):**
* [NE] for `F_hyb` on `m ∈ [1.5, 10]Λ`, `b̂ ∈ [20, 200]` at `Δ_g = 300, 1000` and `[40, 200]` at `Δ_g = 100`, every `J`;
* [NE] for the exact `F` on 56 sampled states up to `b̂ ≈ 100`;
* [D + NE] for the large-`b̂` leading envelope (all `m ≥ 1.5Λ`, `b̂ ≥ 100`);
* [I] for `|S| < F` at `b̂ > 100` and off the sample (calibrated `|S|/F ≤ 131/m²`, margin ≥ 150).

The earlier statement "the `OMEGA_U` κ-envelope applies to M02 only; [I] for D-safe at `b̂ > 20`" is replaced by the above.
**Still untouched:** H2/H4 (gravitons at `x ~ 1`, `+−` transfer) and H12′.


## 1. Conventions (checked in the code)

* `m = (m/Λ)·Λ` and `Λ = √((Δ_g−2)(Δ_g−1))`.
* **Impact parameter** `b̂ = (2J+1)/(m/Λ) = (2J+1)Λ/m`, the convention of `e7_nf.py` and `nf_common.action_nf`. The `b̂ ≈ 2J/m` of the
  task is this up to the `+1`. So `J ≈ ½ b̂ m/Λ`, and the region `m ≤ 10Λ`, `b̂ ≤ 200` needs `J` up to 1000.
* **AdS impact parameter** `β`, with `cosh β = η = 1 + (2J+1)²/(2m₂)`, so `2 sinh(β/2) ≃ b̂/Λ`. Within `b̂ ≤ 200`:
  * `β ≤ 0.2` at `Δ_g = 1000`;
  * `β ≤ 0.66` at `Δ_g = 300`;
  * `β ≤ 1.8` at `Δ_g = 100`.
* **Twist gap.** It bounds the impact parameter: `b̂ ≲ Λ(m/Λ − Λ/m)`, i.e. `β ≤ 2 log(m/(Δ_g−2))`. At fixed `m` there is no limit
  `b̂ → ∞`: large `b̂` requires large `m`. At `Δ_g = 100` the physical `b̂ ≤ 200` region starts only at `m ≳ 2.1Λ` (`b̂_max = 82` at
  `1.5Λ`, `148` at `2Λ`).
* **Action** (`h13_common.py`). The exact D-safe action of a state is `F = F_hyb + S`:
  * `F_hyb = B₂ + B₃ + Σ_k curv_k` is the flat stripped NF02 action (exact Legendre/Jacobi, Gauss rule exact for every `J`) plus the
    leading AdS Regge curvature `curv_k = u^{k−2}2ε_kσ^k[∫φ_kP_{−½+iΛp}(cosh β) − ∫φ_kJ₀(b̂p)]` (exact Mehler–Dirichlet transforms,
    `nf_common.TransformsNF`);
  * `S = Σ_k ∫u^{k−2}σ^kφ_kE_kρ_k` is the exact scalar-skeleton remainder, with the same `ρ_k` as `e8_states_arb.py`;
  * `N_s = |B₂+curv₂| + |B₃+curv₃|`, as in `e7_nf.py`.

  Reproduced to `10⁻¹⁵` on cached `e7` records (`F_h`, `N_s`, `S`).
* **Kernel `K*`.** The normalised heavy moments of a flip-free even kernel are the kernel-independent preserved densities, and `e10`'s
  normalised response ratio is 1 on every node (`EXACT_SKELETON_NEARGAP.md` §9). So `F` is the `K*` action, and equally the `𝔖_ff` one.
* **Skeletons.** `Δφ = 3` and `2`, both helicity channels. `++` uses even `J` only. `+−` uses the scalar `ρ` transferred (reading B [I], as before).

## 2. Task 1a: the hybrid action on the whole region, every `J` (`h13_hyb.py`) [NE; ingredients exact]

**Grid.**
* Every integer `J` with `b̂` in range, both helicities (`++` even `J`; `+−` with `J ≥ 4`).
* Masses: `m/Λ ∈ [1.5, 10]` at 53 values (step 0.1 to 3, 0.25 to 5, 0.5 to 10), plus refinements near the minima: 0.005 steps on
  `[1.5, 1.6]` (`Δ_g = 300, 1000`), 0.02 on `[3, 3.5]` (300), 0.01 on `[1.6, 1.8]` (100).
* `b̂` ranges: `[20, 200]` (`Δ_g = 300, 1000`) and `[40, 200]` (100). Extensions: `b̂ ∈ [200, 1000]` at `m/Λ = 1.5, 2, 3, 5, 10`
  (`2, 3, 5, 7, 10` at `Δ_g = 100`), and `m/Λ = 20, 50, 100` at `b̂ ∈ [20, 200]` with every 7th `J`.
* About 13 000–34 000 (state, helicity) records per `Δ_g` and range.

**Result: no non-positive state anywhere. `F_hyb/N_s ≥ 0.93`.**

| `Δ_g` | range | `++` min `F_hyb/N_s` (where) | `+−` min `F_hyb/N_s` (where) |
|---|---|---|---|
| 300 | `b̂ ∈ [20,200]`, `m ≤ 10Λ` | 0.9886 (`3.25Λ`, `J = 32`, `b̂ = 20.0`) | **0.9312** (`1.5Λ`, `J = 15`, `b̂ = 20.7`) |
| 1000 | same | 0.9885 (`3.25Λ`, `J = 32`, `b̂ = 20.0`) | **0.9311** (`1.5Λ`, `J = 15`, `b̂ = 20.7`) |
| 100 | `b̂ ∈ [40,200]`, `m ≤ 10Λ` | 0.9998 (`1.62Λ`, `J = 32`, `b̂ = 40.1`) | **0.9410** (`1.68Λ`, `J = 90`, `b̂ = 107.7`) |
| 300 | `b̂ ∈ [200,1000]` | 1.0000 | 0.9186 (`1.5Λ`, `J = 186`, `b̂ = 248.7`) |
| 1000 | `b̂ ∈ [200,1000]` | 1.0000 | 0.8990 (`1.5Λ`, `J = 624`, `b̂ = 833`) |
| 100 | `b̂ ∈ [200,1000]` (physical part) | 1.0000 | 0.9721 (`3Λ`, `J = 396`, `b̂ = 264`) |
| 300 / 1000 | `m/Λ = 20, 50, 100`, `b̂ ∈ [20,200]` (every 7th `J`) | ≥ 0.9997 | ≥ 0.9996 |

* **Where the minima sit.**
  * At `Δ_g = 300` and `1000` the minimum is the corner of the region: `m = 1.5Λ` at the smallest admissible `b̂`. The 0.005-step
    refinement in `m` does not lower it (0.9341 at `1.505Λ`).
  * It continues smoothly into the already-scanned `b̂ < 20` region, where the D-safe minimum is `0.20` (`+−`, `J = 5`, near the gap).
  * At `Δ_g = 100` the `+−` minimum is interior (`1.68Λ`, `b̂ ≈ 108`). The 0.01-step refinement moves it from 0.9413 to 0.9410.
* **Per band** (`cache/h13_summary.txt` §1, `+−`, `Δ_g = 300`): `F_hyb/N_s` is 0.93–0.995 on `b̂ ∈ [20,40)`, ≥ 0.98 on `[40,70)`,
  ≥ 0.97 on `[70,120)` and ≥ 0.935 on `[120,200]`. It rises with `m`, reaching ≥ 0.994 for `m ≥ 5Λ`. `++` is ≥ 0.989 everywhere.
* **What `F/N_s ≈ 1` means.** At large `b̂` both components `B₂+curv₂` and `B₃+curv₃` are positive, apart from the `σ = −1` sign of the
  `B₃` part in `+−`. So `F ≈ N_s`, and the margin against any correction of relative size `ε` is `1 − ε`.

## 3. Task 1b and Task 3: the exact scalar-skeleton remainder at large `J` (`h13_exact.py`) [NE + M]

### 3.1 What made `J ~ 100` reachable, and a numerical problem found on the way

1. **Mack polynomials.** `mack.mack_P_fast` takes 138 s at `J = 30` (measured) and scales like `J⁶–J⁷` (an extrapolated ~15 min at `J = 40`).
   * `mack_flint.py` is the same algorithm on FLINT `fmpq_mpoly`. Output is identical (bit-for-bit, `J ≤ 13`, 3 values of `Δ`). It is
     ~200× faster.
   * `mack_mat.py` reorganises the Costa–Gonçalves–Penedones sum as one integer matrix product per `r`, after noting that the four
     Pochhammer indices always sum to `2J`. Identical to `mack_flint` for `J ≤ 21`, 3 values of `Δ`.
   * Timings: `J = 50` in 2.5 s, `J = 70` in 18 s, `J = 104` in 210 s.
   * States are placed on `Δ ∈ (1/64)ℤ` (`h13_common.exact_state`), so the rationals stay small.
2. **The `n`-sum tail of `skelexact3.Laplace3` is inaccurate at large `J`.** This is new, and relevant to anyone extending the Laplace
   method. It does **not** affect the existing `J ≤ 28` results.
   * The Euler–Maclaurin tail used double-precision Gauss–Legendre nodes and weights. At `J ≳ 40` the `n`-sum cancels by many orders
     of magnitude: the per-`n` contributions grow like `n^{~25}` and oscillate.
   * `10⁻¹⁶` weight errors then gave `O(1)` errors at large `ν`, invisible in the Arb radius. Example: `Δ_g = 300`, `Δφ = 2`,
     `m = 1.5Λ`, `J = 52` gave `ρ₂ ≈ −23` at `p ≥ 0.85` and `S/N_s = +0.12`, against `−5·10⁻⁵` once fixed.
   * **Fix (`h13_laplace.Laplace3H`):** tail nodes and weights in Arb (`arb.legendre_p_root`). Everything else is a verbatim copy of
     `Laplace3`.
   * **Residual problem.** The Euler–Maclaurin tail itself becomes unreliable when `J` is not small compared with `τ = Δ − J`. The
     error grows with `N₀` beyond ~400 and is large for small `N₀`; higher-order Euler–Maclaurin corrections make it worse
     (`t_h13_nsum.py`, `cache/t_h13_nsum.txt`).
   * **Error estimate.** Each state is therefore computed with 3 (earliest records: 2) `(N₀, n_q, L)` configurations. The spread of
     `S` is the "`n`-sum err" column, and records with spread `> 10⁻³ N_s` are flagged **unreliable** and excluded.
   * **Small `N₀` also fails, at `Δ_g = 1000`.** With `J ≳ 74`, `N₀ = 120–200` gave spreads up to `10¹⁴N_s` while `N₀ = 320`
     was stable. The final configuration set is `(480,480,½)`, `(320,480,¼)`, `(320,320,½)`, and the flagged records were redone with
     it (`H13_REDO=1`).
   * **Three states remain unreliable**, all at `m = 1.5Λ` with large `J`:
     * `Δ_g = 100`: `J = 60` (`b̂ = 81`);
     * `Δ_g = 300`: `J = 112` (`b̂ = 150`, where `S/N_s` came out `−476`: plainly wrong);
     * `Δ_g = 1000`: `J = 112` (`b̂ = 150`).
   * **`J ≥ 140` not run.** `mack_mat` at `J = 150` needs about 6 GB and was OOM-killed (the integer matrices are not chunked).
     Combined with the `n`-sum problem at `m ≈ 1.5Λ`, **no reliable exact state has `b̂ > 100.5`**.
3. **`value_fast`.** The same double sum is evaluated in transposed-Horner order, `I(ν) = Σ_m g_m(ν)R_m`. It agrees with
   `Laplace3.value` to `≤ 10⁻⁹⁰` relative (`t_h13_fast.py`).
4. **Validation.**
   * Against the cached `e8` Arb records at `J = 16, 20` (`Δ_g = 300`; `t_h13_check.py`): identical (`5·10⁻¹⁷`) with the old
     weights; `≤ 2·10⁻⁹` with Arb weights and the new `N₀`.
   * `S`, `F_h` and `N_s` reproduce `e7_300_recs` to `10⁻¹⁵`.
   * The `p`-rule has 64 nodes for `b̂ ≤ 100` and `0.75b̂ + 32` beyond (Gauss error `~(e·1.35b̂/8n)^{2n}`). Quadrature checks (`t_h13_quad.py`):
     * `Δ_g = 300`, `1.5Λ`, `J = 74`, `b̂ = 99`: 64 and 128 nodes agree to `6·10⁻⁸` relative;
     * `Δ_g = 1000`, `2Λ`, `J = 100`, `b̂ = 100.5`: 107 and 214 nodes agree to `10⁻⁷`, 64 and 214 to `10⁻⁵`.
   * Odd `k`: quartic interpolation in `k` [M] as in `e8`. The cubic-vs-quartic change of `S` ("odd-k" column) is `≤ 1.9·10⁻⁵ N_s` on
     every reliable record, and `≤ 8·10⁻⁶ N_s` for `m ≥ 2Λ`.

### 3.2 Results

**The sample** (both skeletons, both channels, `J` even) is 25 states at `Δ_g = 300`, 25 at 1000 and 9 at 100 (24 / 24 / 8 reliable):
`m/Λ ∈ {1.5, 2, 3, 4, 5, 7, 10}`, `b̂ ≈ 20–100`, `J ≤ 120`.

**Table** (`cache/h13_exact_table.md`). Each entry is the worst over the reliable skeletons. "—" marks no reliable value.

**Δ_g = 100**

| `m/Λ` | `J` | `b̂` | `β` | `++` `F/N_s` | `+−` `F/N_s` | max `|S|/N_s` | max `|S|/F_hyb` | `n`-sum err/`N_s` | status |
|---|---|---|---|---|---|---|---|---|---|
| 1.50 | 30 | 40.7 | 0.410 | +1.0008 | +0.9777 | 3.1e-03 | 3.2e-03 | 2e-07 | ok (Δφ=3,2) |
| 2.00 | 40 | 40.5 | 0.408 | +1.0008 | +0.9863 | 2.5e-03 | 2.5e-03 | 1e-08 | ok (Δφ=3,2) |
| 1.50 | 46 | 62.0 | 0.620 | +1.0002 | +0.9673 | 1.1e-03 | 1.2e-03 | 2e-04 | ok (Δφ=3,2) |
| 1.50 | 60 | 80.7 | 0.798 | — | — | — | — | 1e-01 | **unreliable** (n-sum) |
| 3.00 | 60 | 40.3 | 0.407 | +1.0000 | +0.9938 | 2.8e-04 | 2.8e-04 | 5e-15 | ok (Δφ=3,2) |
| 2.00 | 70 | 70.5 | 0.701 | +0.9996 | +0.9758 | 1.3e-03 | 1.3e-03 | 5e-04 | ok (Δφ=3,2) |
| 5.00 | 100 | 40.2 | 0.405 | +0.9999 | +0.9977 | 1.3e-04 | 1.3e-04 | 0e+00 | ok (Δφ=3,2) |
| 3.00 | 104 | 69.7 | 0.693 | +1.0001 | +0.9900 | 1.9e-04 | 1.9e-04 | 1e-07 | ok (Δφ=3,2) |
| 4.00 | 120 | 60.2 | 0.603 | +1.0000 | +0.9953 | 6.4e-05 | 6.5e-05 | 1e-14 | ok (Δφ=3,2) |

**Δ_g = 1000**

| `m/Λ` | `J` | `b̂` | `β` | `++` `F/N_s` | `+−` `F/N_s` | max `|S|/N_s` | max `|S|/F_hyb` | `n`-sum err/`N_s` | status |
|---|---|---|---|---|---|---|---|---|---|
| 1.50 | 16 | 22.0 | 0.022 | +1.0000 | +0.9492 | 2.0e-05 | 2.0e-05 | 0e+00 | ok (Δφ=3,2) |
| 2.00 | 20 | 20.5 | 0.021 | +0.9936 | +0.9882 | 3.2e-05 | 3.3e-05 | 0e+00 | ok (Δφ=3,2) |
| 1.50 | 22 | 30.0 | 0.030 | +1.0000 | +0.9892 | 5.2e-05 | 5.3e-05 | 0e+00 | ok (Δφ=3,2) |
| 2.00 | 30 | 30.5 | 0.031 | +1.0000 | +0.9886 | 2.9e-05 | 2.9e-05 | 0e+00 | ok (Δφ=3,2) |
| 3.00 | 30 | 20.3 | 0.020 | +0.9936 | +1.0000 | 1.2e-05 | 1.2e-05 | 0e+00 | ok (Δφ=3,2) |
| 1.50 | 34 | 46.0 | 0.046 | +1.0000 | +0.9898 | 2.2e-05 | 2.2e-05 | 0e+00 | ok (Δφ=3,2) |
| 4.00 | 40 | 20.3 | 0.020 | +0.9950 | +1.0000 | 5.9e-06 | 5.9e-06 | 0e+00 | ok (Δφ=3,2) |
| 2.00 | 44 | 44.5 | 0.045 | +1.0000 | +0.9961 | 3.0e-05 | 3.0e-05 | 0e+00 | ok (Δφ=3,2) |
| 3.00 | 44 | 29.7 | 0.030 | +1.0000 | +0.9903 | 5.2e-06 | 5.2e-06 | 0e+00 | ok (Δφ=3,2) |
| 5.00 | 50 | 20.2 | 0.020 | +0.9963 | +1.0000 | 3.5e-06 | 3.5e-06 | 0e+00 | ok (Δφ=3,2) |
| 1.50 | 52 | 70.0 | 0.070 | +1.0000 | +0.9980 | 1.6e-05 | 1.6e-05 | 4e-11 | ok (Δφ=3,2) |
| 4.00 | 60 | 30.3 | 0.030 | +1.0000 | +0.9952 | 7.7e-07 | 7.7e-07 | 0e+00 | ok (Δφ=3,2) |
| 3.00 | 68 | 45.7 | 0.046 | +1.0000 | +0.9983 | 1.1e-05 | 1.1e-05 | 0e+00 | ok (Δφ=3,2) |
| 2.00 | 70 | 70.5 | 0.071 | +1.0000 | +0.9988 | 1.8e-05 | 1.8e-05 | 1e-11 | ok (Δφ=3,2) |
| 7.00 | 70 | 20.1 | 0.020 | +0.9978 | +1.0000 | 1.6e-06 | 1.6e-06 | 0e+00 | ok (Δφ=3,2) |
| 1.50 | 74 | 99.3 | 0.099 | +1.0000 | +0.9981 | 1.3e-05 | 1.3e-05 | 1e-07 | ok (Δφ=3,2) |
| 5.00 | 74 | 29.8 | 0.030 | +1.0000 | +0.9965 | 2.0e-06 | 2.0e-06 | 0e+00 | ok (Δφ=3,2) |
| 4.00 | 80 | 40.2 | 0.040 | +1.0000 | +0.9983 | 4.9e-07 | 4.9e-07 | 0e+00 | ok (Δφ=3,2) |
| 2.00 | 100 | 100.5 | 0.101 | +1.0000 | +0.9989 | 9.6e-06 | 9.6e-06 | 3e-08 | ok (Δφ=3,2) |
| 10.00 | 100 | 20.1 | 0.020 | +0.9988 | +1.0000 | 7.4e-07 | 7.4e-07 | 0e+00 | ok (Δφ=3,2) |
| 3.00 | 104 | 69.7 | 0.070 | +1.0000 | +0.9996 | 7.6e-06 | 7.6e-06 | 3e-13 | ok (Δφ=3,2) |
| 7.00 | 104 | 29.9 | 0.030 | +1.0000 | +0.9982 | 1.0e-06 | 1.0e-06 | 0e+00 | ok (Δφ=3,2) |
| 1.50 | 112 | 150.0 | 0.150 | — | — | — | — | 3e+07 | **unreliable** (n-sum) |
| 5.00 | 112 | 45.0 | 0.045 | +1.0000 | +0.9994 | 4.8e-06 | 4.8e-06 | 0e+00 | ok (Δφ=3,2) |
| 4.00 | 120 | 60.3 | 0.060 | +1.0000 | +0.9995 | 3.9e-06 | 3.9e-06 | 0e+00 | ok (Δφ=3,2) |

**Δ_g = 300**

| `m/Λ` | `J` | `b̂` | `β` | `++` `F/N_s` | `+−` `F/N_s` | max `|S|/N_s` | max `|S|/F_hyb` | `n`-sum err/`N_s` | status |
|---|---|---|---|---|---|---|---|---|---|
| 1.50 | 16 | 22.0 | 0.074 | +1.0002 | +0.9493 | 2.3e-04 | 2.3e-04 | 0e+00 | ok (Δφ=3,2) |
| 2.00 | 20 | 20.5 | 0.069 | +0.9939 | +0.9885 | 3.6e-04 | 3.7e-04 | 0e+00 | ok (Δφ=3,2) |
| 1.50 | 22 | 30.0 | 0.100 | +1.0000 | +0.9885 | 5.9e-04 | 5.9e-04 | 0e+00 | ok (Δφ=3,2) |
| 2.00 | 30 | 30.5 | 0.102 | +0.9999 | +0.9883 | 3.2e-04 | 3.2e-04 | 0e+00 | ok (Δφ=3,2) |
| 3.00 | 30 | 20.3 | 0.068 | +0.9937 | +1.0001 | 1.4e-04 | 1.4e-04 | 0e+00 | ok (Δφ=3,2) |
| 1.50 | 34 | 46.0 | 0.154 | +1.0000 | +0.9882 | 2.3e-04 | 2.4e-04 | 0e+00 | ok (Δφ=3,2) |
| 4.00 | 40 | 20.2 | 0.068 | +0.9951 | +1.0001 | 6.6e-05 | 6.6e-05 | 0e+00 | ok (Δφ=3,2) |
| 2.00 | 44 | 44.5 | 0.149 | +0.9998 | +0.9950 | 3.4e-04 | 3.4e-04 | 0e+00 | ok (Δφ=3,2) |
| 3.00 | 44 | 29.7 | 0.099 | +1.0001 | +0.9903 | 5.8e-05 | 5.8e-05 | 0e+00 | ok (Δφ=3,2) |
| 5.00 | 50 | 20.2 | 0.068 | +0.9963 | +1.0000 | 3.9e-05 | 3.9e-05 | 0e+00 | ok (Δφ=3,2) |
| 1.50 | 52 | 70.0 | 0.234 | +1.0000 | +0.9921 | 1.8e-04 | 1.8e-04 | 2e-15 | ok (Δφ=3,2) |
| 4.00 | 60 | 30.2 | 0.101 | +1.0000 | +0.9952 | 8.1e-06 | 8.1e-06 | 0e+00 | ok (Δφ=3,2) |
| 3.00 | 68 | 45.7 | 0.153 | +0.9999 | +0.9978 | 1.2e-04 | 1.2e-04 | 0e+00 | ok (Δφ=3,2) |
| 2.00 | 70 | 70.5 | 0.236 | +0.9999 | +0.9953 | 2.0e-04 | 2.0e-04 | 5e-19 | ok (Δφ=3,2) |
| 7.00 | 70 | 20.1 | 0.067 | +0.9978 | +1.0000 | 1.8e-05 | 1.8e-05 | 0e+00 | ok (Δφ=3,2) |
| 1.50 | 74 | 99.3 | 0.331 | +1.0000 | +0.9831 | 1.5e-04 | 1.5e-04 | 3e-08 | ok (Δφ=3,2) |
| 5.00 | 74 | 29.8 | 0.100 | +1.0000 | +0.9965 | 2.3e-05 | 2.3e-05 | 0e+00 | ok (Δφ=3,2) |
| 4.00 | 80 | 40.3 | 0.135 | +1.0000 | +0.9981 | 4.6e-06 | 4.6e-06 | 0e+00 | ok (Δφ=3,2) |
| 2.00 | 100 | 100.5 | 0.335 | +1.0000 | +0.9900 | 7.6e-05 | 7.7e-05 | 1e-07 | ok (Δφ=3,2) |
| 10.00 | 100 | 20.1 | 0.067 | +0.9988 | +1.0000 | 8.2e-06 | 8.2e-06 | 0e+00 | ok (Δφ=3,2) |
| 3.00 | 104 | 69.7 | 0.233 | +1.0000 | +0.9981 | 7.3e-05 | 7.3e-05 | 4e-16 | ok (Δφ=3,2) |
| 7.00 | 104 | 29.9 | 0.100 | +1.0000 | +0.9982 | 1.1e-05 | 1.1e-05 | 0e+00 | ok (Δφ=3,2) |
| 1.50 | 112 | 150.0 | 0.497 | — | — | — | — | 6e+02 | **unreliable** (n-sum) |
| 5.00 | 112 | 45.0 | 0.151 | +1.0000 | +0.9992 | 5.3e-05 | 5.3e-05 | 0e+00 | ok (Δφ=3,2) |
| 4.00 | 120 | 60.3 | 0.202 | +1.0000 | +0.9990 | 3.9e-05 | 3.9e-05 | 7e-19 | ok (Δφ=3,2) |


**Minima over the reliable exact records** (no non-positive record):

| `Δ_g` | min `F/N_s` (where) | max `|S|/F_hyb` |
|---|---|---|
| 300 | **0.9493** (`+−`, `1.5Λ`, `J = 16`, `b̂ = 22`) | `5.9·10⁻⁴` |
| 1000 | **0.9492** (same state) | `5.3·10⁻⁵` |
| 100 | **0.9673** (`+−`, `1.5Λ`, `J = 46`, `b̂ = 62`) | `3.2·10⁻³` |

**Scaling of the remainder (all reliable records):**
* `|S|/F_hyb · m² ≤ 131` on `b̂ ∈ [20,40)`, `≤ 121` on `[40,80)` and `≤ 38` on `[80,101]`, with `m` in AdS units.
* So `|S| ≲ 130 F_hyb/m²`: an `O(1/m²)` relative correction that does not grow with `b̂` in the sample.
* At fixed `(m/Λ, b̂)`, `S/F` falls like `Λ⁻²`. For example, `2.4·10⁻⁴ → 2.1·10⁻⁵` from `Δ_g = 300` to 1000 at `(2Λ, b̂ = 20.5)`.
* Since `m ≥ 1.5Λ ≥ 148`, this gives `|S|/F ≤ 6·10⁻³` at `Δ_g ≥ 100` on the sample [NE].

**The remainder cancels like the leading action.**
* The pointwise size of the remainder, `Σ_k∫|u^{k−2}φ_kE_kρ_k|` (column `sum|.|/N_s`), is up to `2·10²N_s` at `b̂ ≈ 100`.
* Its signed integral is `≤ 3·10⁻³N_s` (`≤ 6·10⁻⁴N_s` at `Δ_g ≥ 300`).
* So the exact `ρ_k` oscillates in `ν` with the phase of the leading Regge term and integrates against NF02 with the same
  forward-free cancellation.
* This is the property that `REFEREE_OMEGA_NF.md` M1 / `OMEGA_NF.md` I15 said "nothing argues". It is now observed on the sample [NE].
* It also explains why every pointwise envelope failed: they needed a relative pointwise correction `≲ 4·10⁻³` (`Δ_g = 100`) or
  `≲ 9·10⁻⁵` (`Δ_g = 1000`), while the true pointwise remainder is `O(1)` or larger.

## 4. Task 2: large-`b̂` envelope

### 4.1 Derivation [D]

Write `P_{−½+iν}(cosh β) = (√2/π)∫₀^β cos(νt)K(t)dt` with `K(t) = (cosh β − cosh t)^{−½}` (Mehler–Dirichlet). For a weight
`φ(p) = Σ_{n≥1}c_np^n` with `φ(0) = 0` (true for every NF02 `φ_k`, since `ψ₂, ψ₃ ∝ p(1−p)⁵`), the large-`Λ` endpoint expansion of
`I = ∫₀¹φ(p)P_{−½+iΛp}(cosh β)dp` at fixed `β` is

`I ≃ Σ_{n≥1} n! c_n M_n(β)/Λ^{n+1}` (+ `p = 1` terms `O(Λ⁻⁶)`, since `φ ∝ (1−p)⁵`),

`M_n(β) = (√2/π)∫₀^βK(t)[reg∫₀^∞ν^n cos(νt)dν]dt`.

* **Odd `n`: non-local.** `M_n = (√2/π)(−1)^{(n+1)/2}n!·FP∫₀^βK(t)t^{−n−1}dt`. These vanish in the flat limit, because
  `∫p^{odd}J₀(b̂p)` has no endpoint term.
* **Even `n`: local.** `M_n = (−1)^{n/2}K^{(n)}(0)/√2`. These reproduce the flat Hankel moments.

Hence the **leading Regge D-safe action** `FR = Σ_ku^{k−2}2ε_kσ^kI_k` satisfies

**`FR ≃ A₁(u,σ) M₁(β)/Λ² + A₂(u,σ) M₂(β)/(2Λ³) + A₃M₃/(6Λ⁴) + …`**, with `A_n = Σ_ku^{k−2}2ε_kσ^kφ_k^{(n)}(0)`,

and the pieces are:

* **`A₁ = 2X(0) + σuY(0) = 377.79 + 50.64σu`** (NF02 normalised to `k(0) = 1`). Only `k = 2, 3` contribute, since `φ_{k≥4} = O(p³)`.
  So `A₁ ≥ 355.3` for `u ≤ 1/2.25`, and `A₁ > 0` for all `u ≤ 1`.
* **`A₂ = 2φ₂″(0) + σuφ₃″(0) = −2047.09`**. `φ₃″(0) = 7·10⁻⁵ ≈ 0`, by NF02's construction.
* **`M₂(β) = −(cosh β − 1)^{−3/2}/(2√2) < 0`.** Its flat limit is `−1/β³`, so `A₂M₂/(2Λ³) → −A₂/(2b̂³)`: the known forward-free flat
  tail (`OMEGA_NF.md`; `NOFORWARD_FUNCTIONAL.md` Lemma M).
* **`M₁(β)` is the new, AdS-only term.** Equivalent forms:
  * `M₁(β) = reg∫₀^∞ν P_{−½+iν}(cosh β)dν = ∫₀^∞ 2ν P_{−½+iν}(cosh β)/(e^{2πν}+1) dν`;
  * `M₁(β) = (√2/π)∫₀^β K(t)[cosh(t/2)/(4sinh²(t/2)) − 1/t²]dt`.

  The first form uses the Mehler–Fock identity `∫ν tanh(πν)P dν = 0` for `β > 0`; the second follows from
  `∫₀^∞sin(νt)/(e^{2πν}+1)dν = 1/(2t) − 1/(4sinh(t/2))`. All three forms agree to 16 digits at `β = 0.3–14` (`h13_asym.py`).
  * `M₁(0) = 1/24`.
  * `M₁ > 0` on `β ∈ [10⁻³, 16]` (191 points) and at `β = 20, 30, 40` [NE].
  * Large `β`: `M₁ e^{β/2}β/√2 → 0.4–0.45`, slowly increasing [NE]. Analytically: the integral of the bracket over `[0, ∞)` vanishes,
    and the tail near `t = β` gives `M₁ ≈ √2e^{−β/2}(1/β − O(β⁻²)) > 0`.
  * Positivity on `(0, ∞)` is not certified: interval arithmetic on the second form would do it.
* **Signs.** Both leading terms of the asymptotic D-safe action are positive for every `u ∈ [0, 1]`, every `β > 0` and both
  helicities, because `A₁ > 0`, `M₁ > 0`, `A₂ < 0` and `M₂ < 0`.
* **Size.** The AdS term `A₁M₁/Λ² ≈ A₁/(24Λ²)` dominates the flat tail `−A₂/(2b̂³)` once `b̂ ≳ (12|A₂|/A₁)^{1/3}Λ^{2/3} ≈ 4Λ^{2/3}`.
  That is `b̂ ≈ 85 / 178 / 398` at `Δ_g = 100 / 300 / 1000`.
  * So at large `b̂` the D-safe action is **not** the flat `1/b̂³` tail: it is an `O(Λ⁻²)` AdS term, set by the moments at `ν = O(1)`.
  * This is why OMEGA_NF's flat-calibrated pointwise envelopes needed tolerances of `10⁻³–10⁻⁴`.

### 4.2 Check against the exact transforms and the full hybrid action [NE]

* `FR` against the 2-term asymptotic (`h13_asym.txt`): relative error `≤ 1.3·10⁻¹` at `b̂ = 20`, `≤ 7·10⁻³` at `b̂ = 50`, `≤ 8·10⁻⁴` at
  `b̂ = 100`, `≤ 8·10⁻⁵` at `b̂ = 200`, `≤ 1.4·10⁻⁵` at `b̂ = 500` (`Δ_g = 100, 300, 1000`; `m/Λ = 1.5, 3, 10`).
* `F_hyb` itself (flat Legendre part included) against the same 2-term form, over every record of §2:

  | `Δ_g` | `b̂ ≥ 50` | `b̂ ≥ 100` | `b̂ ∈ [200, 1000]` |
  |---|---|---|---|
  | 100 | [0.9982, 1.0111] | [0.9999, 1.0011] | [1.00001, 1.00004] |
  | 300 | [0.9978, 1.0116] | [0.9999, 1.0026] | [0.99999, 1.00034] |
  | 1000 | [0.9977, 1.0116] | [0.9998, 1.0029] | [0.99999, 1.00072] |
  (Including `m = 20, 50, 100Λ` changes the `Δ_g = 300, 1000` entries by `< 2·10⁻⁴`.)

### 4.3 The envelope statement and its grade

**E (hybrid part).** For every state with `m ≥ 1.5Λ` and `b̂ ≥ 100`:

`F_hyb ≥ (1 − 3·10⁻³)[A₁M₁(β)/Λ² + A₂M₂(β)/(2Λ³)] > 0`.

* **Derived [D]:**
  * the expansion and its coefficients;
  * `A₁ > 0`, `A₂ < 0`, `M₂ < 0`;
  * the closed forms of `M₁`.
* **Numerical [NE]:**
  * `M₁ > 0` (191 + 3 values of `β`, plus the large-`β` asymptotics);
  * the `3·10⁻³` remainder factor, which is verified on the grids of §2 (every `J`, `m ≤ 10Λ`, `b̂ ≤ 1000`; `m = 20, 50, 100Λ` sampled),
    not bounded.
* **Covers `m > 10Λ`:** `A_n` are polynomials in `u` with `A₁(0) = 377.8` and `A₂(0) = −2047`, and `β` depends on `b̂/Λ` only up to
  `O(1/m²)`.
* **Gap:** a rigorous bound on the asymptotic remainder (terms `n ≥ 3` and the Legendre-vs-`J₀` correction) is not derived. It would
  follow from bounding `Σ_{n≥3}` with the exact polynomial `φ_k`. This is the same kind of lemma as `NOFORWARD_FUNCTIONAL.md` Lemma M,
  with `M_n` in place of the Hankel moments.

**E′ (skeleton remainder).**

`F = F_hyb + S ≥ F_hyb(1 − C/m²)`, with `C = 131`, the maximum of `|S|m²/F_hyb` over all reliable exact records:
`b̂ ∈ [20, 100.5]`, `m ≤ 10Λ`, `J ≤ 120`, `Δ_g = 100, 300, 1000`, both skeletons.

* For `m ≥ 1.5Λ` and `Δ_g ≥ 100`, `C/m² ≤ 6·10⁻³`. So positivity holds with a factor `≥ 150` to spare **if** `C` stays below `~2·10⁴`
  (`m²` at `1.5Λ`, `Δ_g = 100`) outside the sample.
* **Grade:**
  * on the sampled states, `[NE]` (Arb arithmetic, with an `n`-sum truncation estimate rather than a bound; odd `k` [M]);
  * as an envelope, `[I]` (calibrated, not derived).
* **Structural support [I]:**
  * the exact moment is an entire function of `ν` (Barnes/Laplace form). If it has the same endpoint structure as `P`, then
    `S ≃ Σ_ku^{k−2}σ^kφ_k′(0)·M₁[E_kρ_k]/Λ² + …` scales like `F`;
  * `S/F` is then set by `ρ` at `ν = O(1)`, where `x = ν/m ≲ 1/(mβ)` and the remainder is the `O(1/m²)` correction of the heavy formula.
    That is what the data show. It is not proven.

**Combined grade for H13′, D-safe, `b̂ > 20`:**
* `[NE]` for `F_hyb` on the whole region `m ∈ [1.5, 10]Λ`, `b̂ ∈ [20, 200]` at `Δ_g = 300, 1000`, and `[40, 200]` at `Δ_g = 100`
  (every `J`, dense `m`);
* `[NE]` for the exact `F` on the sampled states;
* `[D + NE]` for the leading large-`b̂` envelope;
* `[I]` for `|S| < F` outside the sample.

The previous `[I]` (no evidence at all for `b̂ > 44`) becomes a calibrated `[I]` with a margin of `≥ 150×` (§0).

## 5. What is not done / caveats

* **`F_hyb` is the D-safe action minus the skeleton remainder.** Both are scalar-skeleton statements. H2/H4 (gravitons at `x ~ 1`,
  `+−` transfer, reading B) and H12′ (odd-`J` `+−` exchanges) are untouched.
* **Exact coverage is a sample, not a grid**, and it stops at `b̂ ≈ 100` (`m ≤ 2Λ`; 70 at `Δ_g = 100`; 60 at `4Λ`; 45 at `5Λ`;
  30 at `7Λ`; 20 at `10Λ`).
  * `b̂ ∈ (100, 200]` has **no** reliable exact state. There it rests on the hybrid scan plus envelope E′ [I].
  * The two blockers are the `n`-sum tail (§3.1 item 2) and the memory of `mack_mat` at `J ≥ 140`.
  * Covering these needs either an `n`-sum method that does not rely on Euler–Maclaurin at `J ~ τ`, or the analytic (eikonal) form of
    the remainder.
* **Rigour.**
  * The Arb radii (`≤ 10⁻⁴⁰`) cover arithmetic only. The `n`-tail truncation is estimated by spread, and the `p`-quadrature is checked
    by doubling, not bounded.
  * `M₁ > 0` and the asymptotic remainder are `[NE]`.
  * Nothing here is `[certified]`.
* **Odd `k`.** Interpolated [M]. The sensitivity (cubic vs quartic) is negligible (`≤ 1.9·10⁻⁵N_s` on reliable records).
* **Not run:** `Δφ = 2.5`; odd `J` in the exact sample (`+−` uses the even-`J` states; the hybrid scan has every `J`).

## 6. Files

| file | role |
|---|---|
| `h13_common.py` | conventions; flat rows for many `J`; hybrid action; `Δ ∈ (1/64)ℤ` states |
| `h13_hyb.py` | hybrid scan, every `J` → `cache/h13_hyb_<Dg>_<b0>_<b1>.pkl` |
| `h13_exact.py` | exact skeleton remainder on sampled states (3 `n`-sum configurations) → `cache/h13_ex_<Dg>.pkl` |
| `h13_laplace.py` | `Laplace3H`: `skelexact3.Laplace3` with Arb tail quadrature, `value_fast`, optional higher-order Euler–Maclaurin (worse; off) |
| `mack_flint.py`, `mack_mat.py` | exact Mack polynomials (FLINT), identical to `mack.mack_P_fast` |
| `h13_asym.py` | `M₁`, `M₂`, `M₃`; `A_n`; asymptotic vs exact `FR` → `cache/h13_asym.txt` |
| `h13_summary.py` | all tables → `cache/h13_summary.txt` |
| `t_h13_check.py`, `t_h13_fast.py`, `t_h13_nsum.py`, `t_h13_quad.py` | checks (§3.1) |

The runs used ≤ 2 cores. The chains are in `cache/h13_{hybchain,exA,exB,exC}.sh`; everything is resumable from the `.pkl` caches.
Superseded records, computed with double-precision tail weights, are kept as `cache/h13_ex_*_OLD_doubleweights.pkl` and are not used.


## 7. Follow-up (18 Nov 2026, later): items that do not depend on the functional

* **`M₁(β) > 0` for all `β > 0` [certified].** `h13_cert_m1.py`, `cache/h13_cert_m1.txt`. The proof has three parts:
  * **(a) `β ≤ 4.7`.** Here `c(t) > 0` on `(0, 4.7]`: by a series bound on `(0, 2]`, and by Arb on 2700 subintervals of `[2, 4.7]`.
    So the integrand is positive.
  * **(b) `β ∈ [4.7, 12]`.** Arb lower Riemann sums of `I(β) = ∫₀¹ g(v)dv`, with `t = β(1−v²)`, on 565 `β`-boxes. All are positive;
    the smallest lower bound is `2.25·10⁻⁴`.
  * **(c) `β ≥ 12`, analytic [D].** It uses the exact identity `∫₀^β c = 1/β − 1/(2 sinh(β/2))` and the bound
    `I/K(0) ≥ 1/β − 1/(2 sinh(β/2)) − √2(e^{−β/2}/4.7 + 4/β²) − 2√(2 ln 2)/(β − ln 2)²`, which is positive for `β ≥ 12`.
* **Two more exact identities [D; checked to 15 digits].**
  * `cosh β − 1 = b̂²/(2Λ²)` exactly. Hence `M₂/Λ³ = −1/b̂³`: the two-term envelope's `A₂M₂/(2Λ³)` *is* the flat tail `−A₂/(2b̂³)`.
  * `(2/π)∫₀^{π/2}[1/s² − w(s)/τ(s)²]_{s=b̂ sin θ}dθ = M₁(β)/Λ²`, with `τ = 2Λ asinh(s/2Λ)` and `w = (1+s²/4Λ²)^{−½}`.
* **Methods.**
  * `mack_mat.py` now uses only the pairs with `k₁₃ ≤ k₂₄` and accumulates the matrix product in chunks. It is bit-identical (24
    self-tests), 2× faster, and uses 0.13 GB at `J = 80` (previously about 6 GB at `J = 150`).
  * The `Δ_g = 1000` large-`J` failures of §3.1 were caused by the tail Gauss order: `n_q = 320` does not resolve the tail at `ν → Λ`,
    while `n_q ≥ 800` converges to `10⁻¹¹`. The new configurations are `(200, 960, 1)`, `(320, 1280, 0.35)`, `(320, 960, 0.5)`. With
    them `m = 1.5Λ`, `J = 112` at `Δ_g = 1000` (`b̂ = 150`) is reliable.
  * Near the twist gap (`Δ_g = 100`, `J = 60`; `Δ_g = 300`, `J = 112`, both at `1.5Λ`) the Euler–Maclaurin tail still fails.
