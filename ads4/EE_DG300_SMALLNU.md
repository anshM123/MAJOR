# `Ω^EE` at `Δ_g = 300`, `J ≥ 24` near the gap: the graviton diagonal at small `ν` (IN PROGRESS — hand-over at the end)

*19–20 Nov 2026, EE-Δ_g=300-small-ν agent. Unrefereed, not committed. No existing document or project code edited. New code
`ads4/code/spinfunc/z*_*.py` (+ `z2_queue*.sh`), caches `spinfunc/cache/z*`. Project modules imported read-only (in-process wrappers only).*

Grades: **[P]** proved, **[D]** derived, **[NE]** numerical evidence, **[M]** model, **[I]** open.

## 0. Status (updated as results come in)

* Fold of the symmetric crystal saddle located at large `J` (§1): `ν* ≈ 11.3–11.5`, independent of `J` and `Δ` [NE].
* The failure of the direct method at the heavy reference (`Δ_ref = 2·10⁴`, `J = 24`) is a floating-point overflow of the scalar
  prefactor `Q_J`, not a contour problem (§2); fixed in-process; the n3 contour then works at the reference with 1 % errors at `2⁸·16`
  points [NE]. No change of `Δ_ref` is needed.
* Production of the `Δ = 336` row nodes `ν ≈ 4–33` for `J = 24, 28` (state + reference, common random numbers): running.

## 1. Fold location at large `J` (`z1_fold.py`, `cache/z1_fold.txt`) [NE]

Symmetric saddle of the scalar flow function `Φ_ν` continued downwards from `ν = 16–22` (step 0.25–0.5), soft `4×4` determinant:

| block | `det_soft` at `ν = 12.5 / 12.0 / 11.5` | last symmetric saddle | `A±` pair |
|---|---|---|---|
| `Δ = 336`, `J = 24` | `1111 / 790 / 404` | `ν = 11.5` (Newton fails at 11.0) | exists at 11.2 (asym 0.19), not at 11.45 |
| `Δ = 336`, `J = 32` | identical to `J = 24` to all printed digits | 11.5 | exists at 11.2 (same `Φ₀`) |
| `Δ = 2·10⁴`, `J = 24` | `1113 / — / —` (`954` at 12.25) | scan stopped at 12.25 (same values as `Δ = 336`) | — |
| `Δ = 4000`, `J = 24` | `2329` at 14.25 (`Δ = 336`: `≈ 2330`) | stopped | — |

* `Φ_ν` (n1_direct.phi_w) is `J`-independent at the saddle family (the `(J/2) log Z²W²` term cancels the `J`-dependence of the `3pt`
  logs), and its soft block is `Δ`-independent to `≲ 0.2 %` between `Δ = 336` and `2·10⁴`. **So `ν* ≈ 11.3–11.5` for every `(Δ, J)` used
  here** (H0 §25.2 had 11.28 at `Δ = 100`, `J = 4`). The fold does not move with `J`; the large-`J` difficulty is in the prefactors
  (`P_J`, the spin-`J` structures), not in the contour.
* Production rule (as NEARGAP100 §8): n3 hybrid contour for `ν ≥ 11` (`ν_f = 13` below 13), `A₊ + A₋` thimbles (n4) below.

## 2. Why the direct method failed at the heavy reference, and the fix (`z2_points.py`) [NE]

* At `Δ = 2·10⁴`, `J = 24`, `ν = 20.85`, every point of the n3 contour (4090/4096; also with the reference's own proposal) was masked.
  Cause: `l1_crystal.build` evaluates the scalar prefactor `Q_J = Σ c_p (Z·W)^p (Z²W²)^{(J−p)/2} = |Z|^J|W|^J P_J(cos)` as a polynomial;
  at the heavy reference `|Z||W| ≈ 10¹³` at the crystal, so `Q_J` overflows double precision for `J ≳ 20` (`inf`/`nan`), and the
  contour's ok-mask (scalar) also gates the graviton kinds. At `J ≤ 10` (`Q ~ 10¹³⁰`) it does not overflow, which is why the
  `J ≤ 10` references of H0 §24–25 worked. (The real-contour failure of n1 at `ν ≲ 35` reported in EE_NEARGAP §5 may be this same
  overflow, or walls; not separated here [I].)
* Fix (in-process copy of `l1_crystal` in `n1_pw` only): `Z·W → Z·W/s`, `Z²W² → Z²W²/s²` with a per-block constant `s = |Z²W²|^{1/2}` at
  the block's saddle, i.e. the scalar integrand is multiplied by the constant `s^{−J}`; the stored scalar normalisation gets `+J log s`.
  The graviton integrands do not use `Q_J` and are unchanged.
* Test (`J = 24`, `ν = 20.85`, `2⁸·16` RQMC points, proposal adapted on the state and re-used): reference `Δ = 2·10⁴` now has 0 masked
  points and relative errors 1.0 % (EE) / 1.1 % (C3H) / 1.0 % (scalar), as the state `Δ = 336` (0.95 / 1.2 / 0.9 %).
* **Consequence for task 1's alternative reference:** the heavy reference is computable as it is; the normalisation `N̂` (state/reference
  ratio at `Δ_ref = 2·10⁴`, H0 §24.1, l9_analysis) need not be re-derived.

## 3. Error structure: common random numbers and the scalar control variate [NE]

At the test point (`J = 24`, `ν = 20.85`, `2⁸·16` points):
* The `2×2` graviton matrices are nearly rank-1 with a common fluctuating factor: element ratios within a matrix have relative errors
  `4·10⁻⁵` (state) / `4·10⁻⁸` (reference), while each element has 0.85 % (state) — all four elements share the same relative error.
* State/reference ratio of an element with common random numbers: 0.85 % (vs 1.4 % if independent). In `Γ` units (`m²×` relative,
  `m = 333.6`) that is `≈ 950` c₂ — too large; it falls as `N^{-1/2}`.
* Graviton/scalar ratio on the same samples: 0.33 % (state) / 0.23 % (reference) — the common factor is largely shared with the scalar.
  Since the exact scalar state/reference ratio is known (exact Mellin skeleton, n4_exact; `J` even), the estimator
  `G = [(g/s)_state/(g/s)_ref] × (s_state/s_ref)_exact` (scalar control variate) is ≈ 2.5× more precise. Both estimators are reported.

## 4. Interim (20 Nov, 00:45) — superseded by the final sections below when they are written

* Production (`z2_queueA.sh`: n3, `Z2_M = 11`, i.e. `2¹¹·2·8` points per block; `z2_queueB.sh`/`B2`: n4 thimbles, `Z2_M = 9`) on the
  `Δ = 336` rows, nodes `x = 0.012 … 0.1` (`ν ≈ 4.0 … 33.4`), each with `Δ_ref = 2·10⁴` on common random numbers.
* First two nodes, `J = 24`: scalar check `Re(d − d_ref) = −0.0024 ± 0.0017` (`ν = 13.34`), `+0.0018 ± 0.015` (`ν = 16.68`) [NE].
* **At `ν = 16.68` (`x = 0.05`, the lowest `K₂` node of the row) the direct `Γ_EE,++` is `+232` (`+−`: `+258`) against `K₂`'s `−1331`
  (`−1317`)**, and at `ν = 13.34` it is `+282` / `+297`. The `K₂` node at `x = 0.05` — the node from which the q3 spline extrapolates —
  has the wrong sign. With only these two direct nodes (NL weight below `ν = 25`, held below 13.34) the `J = 24` minimum is `+0.91`
  (at 1.3Λ; the 1.1Λ state is no longer the worst) [NE + M, preliminary].

* **Interim 2 ( `cache/z4_chain_int3.txt`, CV estimator, min-norm NL weight, `Γ_NL` held below the lowest direct node):**
  direct nodes so far `J = 24`: `ν = 11.68, 13.34, 16.68, 20.85`; `J = 28`: `13.33`. Scalar checks `|Δc₂| ≤ 99` (all ≤ 1.4σ).
  `Γ_EE,++` (complex `Γ_c`, real part) at `ν = 11.68 / 13.34 / 16.68 / 20.85`: `+440 / +358 / +223 / −30` (`J = 24`); `+280` at 13.33
  (`J = 28`). Chain: `J = 24` min `+0.853` (1.5Λ), `J = 28` min `+0.873` (1.5Λ) — the 1.1Λ states no longer bind; `n = 1` at 300
  (binding `J = 5`, `+0.2052`) [NE + M, preliminary: below `ν ≈ 11.7` (`J = 24`) / `13.3` (`J = 28`) still the held value].
* **Interim 3 ( `cache/z4_chain_i4_*.txt`):** `J = 24` n3 nodes complete (`ν = 11.68 … 33.4`), `J = 28` has `13.33,
  16.66`. Chain (`J = 24 / 28` min `λ/N_s`): CV-hold `+0.828 / +0.897`; direct-hold `+0.827 / +0.897`; CV-zero `+0.815 / +0.897`;
  CV-linear `+0.831 / +0.897`; CV-hold with `K₂` kept at the shared nodes `+0.839 / +0.873`. `J = 24` at 1.1Λ: `+0.905` (was `−0.053`).
  Stress test (uniform real `Γ_c = C` below the lowest direct node, `ν < 11.68` / `13.33`): `C = −4000` → `+0.199 / +0.309`;
  `−8000` → `−0.479 / −0.755`; `+8000` → pass. **Failure needs `Γ_c ≲ −5200` c₂ uniformly below the fold** (zero crossing by linear
  interpolation, both `J`) — the n4 nodes test exactly this.
* **Interim 4 (first n4 node, `J = 24`, `ν = 10.0`, `x = 0.03`):** scalar check `Re +0.0018 ± 0.0155` (`Δc₂ = +108 ± 908`),
  `Im +0.099 ± 0.046` (2.2σ phase offset, as NEARGAP100 §4 saw at 5.45). Element ratio errors: EE 4.7 % direct, **0.95 % with the scalar
  CV**; C3H 7.2 % / 5.2 %. `Γ_c` (CV): EE `−948 / −979` (`++/+−`), C3H `+6765 / −6021` (C3H error ≈ ±4000). The direct estimator
  carries the thimble phase error into `Im Γ_c` (`+7755i`); the CV removes it. Chain with this node held below 10 (CV): `J = 24`
  `+0.520` (1.1Λ), `J = 28` `+0.897`. The held value below `ν = 10` (where `X/φ ≈ 2–5` at `ν ≲ 4`) now drives `J = 24`; nodes 8.34,
  6.67, 5.0 are running.
* **Interim 5 (`J = 24` n4 nodes `ν = 10.0, 8.34`; `J = 28` n3 complete):** scalar checks at 8.34: `Re +0.012 ± 0.018`
  (`Δc₂ = +852 ± 1302`). `J = 24` `Γ_c,EE` (CV) at 8.34: `+3368 − 3290i` (next to the zero 6.64 `Γ_c` is large; the NL density
  `Re[(Y/φ)Γ_NL] = +388`). Chain `J = 24 / 28`: CV-hold `+0.807 / +0.812`; CV-zero `+0.726 / +0.788`; CV-lin `+0.847 / +0.818`;
  direct-hold `+0.822 / +0.806` (`J = 28`: n3 nodes only, lowest 11.66). `z5_queueC.sh` (`J = 32`) is set to start when `R1` ends.

## 9. Hand-over (kept current)

* **State (20 Nov, after the container restart at ≈ 00:25):** done and on disk: `J = 24` nodes `ν = 13.34, 16.68` (n3, state 336 +
  reference 2·10⁴; `cache/z2_pts_{EE,C3H}_J24.pkl`, `cache/z2_raw.pkl`). Lost: the in-flight `J = 24` n3 node and the `J = 28` n4 node
  (`ν = 10`, state only). Resumed with `z2_queueR1.sh` (`J = 24`: n3 `ν = 11.68, 20.85, 25.0, 33.4`, then n4 `10.0, 8.34, 6.67, 5.0`)
  and `z2_queueR2.sh` (`J = 28`: all n3 nodes, then n4); logs `cache/z2_queueR{1,2}.log`. `z2_points.py` skips finished nodes, so
  re-running a queue resumes it. The old `z2_queue{A,B,B2}.sh` are superseded. `z5_queueC.sh` (`J = 32`: `K₂` row on 336 + direct nodes)
  has not been started.
* **Code.** `z1_fold.py` (fold scan), `z2_points.py` (direct points, CRN, `Q_J` overflow fix, n3/n4 switch at `ν = 11`), `z3_nodes.py`
  (scalar validation vs exact Mellin `→ cache/z3_validate.txt`; conversion to chain records, estimators `dir`/`cv`, jackknife),
  `z4_chain.py` (w2 chain + direct nodes + non-local `k = 2` term below `ν = 25`; env in its docstring).
* **Re-run the analysis:** `python3 z3_nodes.py validate 24 28`;
  `Z4_JS=24,28 Z4_EST=dir|cv Z4_GAUGE=min|can Z4_EXTRAP=hold|zero|lin Z4_TAG=… Q3_JEVEN=24,28 Q3_JODD=5
  Q3_EXTRA_M=1.02,1.05,1.1,1.15,1.2,1.3,1.4,1.5,1.75,2,2.5,3,4,5 python3 z4_chain.py 300` (≈ 1 min);
  `Z4_DIRECT=0 Z4_NL=0` reproduces EE_NEARGAP's convention A (`−0.0534 / −0.0950`) exactly (checked).
