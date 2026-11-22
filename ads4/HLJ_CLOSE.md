# Closing `H_LJ(300)`? The spin-J engine with an analytic normalisation, its convergence structure, and the `ν < 8.3` windows

*22 Nov 2026, H_LJ-close agent. Unrefereed, not committed. No existing document, code file or cache was edited or overwritten. New code:
`ads4/code/spinfunc/lj_rs_*.py`, `lj_rs_*.sh`. New caches: `spinfunc/cache/lj_rs_*`, plus the unchanged `hc_chain.py`'s own outputs
`cache/hc_chain_z4_lj_rs*` (new names). Read-only imports: `bj_ee`, `bj_eeq`, `bj_chain` (functions that write nothing), `bj_gamma`,
`fw_gamma`, `bf_core`, `z3_nodes`, `bg_*`, `hc_chain` (run unchanged, in-process wrappers). ≤ 2 worker processes, `OMP/OPENBLAS = 1`,
`ulimit -v 3000000` on every job (one exploratory scratch script without the limit was killed; it wrote nothing). Nothing in
`code/resultq` was touched.*

Grades: **[P]** proved; **[D]** derived; **[D-sketch]**; **[DV]** derived and checked numerically; **[NE]** numerical evidence; **[M]** model;
**[A]** assumption; **[I]** open. `c₂` = `m² ×` relative correction. "Engine" = the closed-form spin-J analytic-radial `1/Δ` engine
(BELOW_FOLD_SPINJ, `bj_*`), images `I_ac(ν)` in the 2×2 conserved basis; order `aᵏ` = relative geometric order `k` (depth `k+1`).
"Direct" = the independent n3/z2 crystal records (state 336 and reference 2·10⁴ on common random numbers, CV estimator).

---

## 0. Summary

PLACEHOLDER-SUMMARY

---

## 1. What the obstruction actually is (diagnosis from the cached images) [DV + NE]

### 1.1 The large `a¹` term is a J-independent EE normalisation, not a `J²/Δ` effect
Per entry, `m²·P₁/P₀` (state, local ratio of the `a¹` to the `a⁰` image), ν = 2–10, away from the entries' common zeros (ν ≈ 1.07, 5.5–7, 13.3):

| state | `a¹` | `a²` | `a³` | `a⁴` | `a¹` coefficient `f₁ = Δ·(a¹ term)/m²` |
|---|---|---|---|---|---|
| J = 0, Δ = 300 (`bg_eeq_D300_J0_K4`) | +6230 ± 40 | +270 ± 40 | +9 | +0.3 | 20.8 |
| J = 28, Δ = 336 (`bj_eeq_D336_J28_d3`) | +4800 … +6500 | +160 … +305 | — | — | 15–19 |
| J = 100, Δ = 1100 (`bj_eeq_D1100_J100_d2_p`) | +1.7 … 2.1·10⁴ | — | — | — | 16–19 |
| reference Δ = 2·10⁴, J = 28 / 100 | 4.2·10⁵ (`m²` of 2·10⁴) | ≤ 0.6 % of it | — | — | 20.9 / 20.9 |

* The `a¹` term is (to 0.3 %) **the same for all four entries** and nearly ν-independent; its coefficient `f₁ ≈ 16–21` is independent of `J`
  and of `J/Δ` (it is present at J = 0). It is the radial-moment normalisation of the EE density (`∝ X₁^{9.25}`, `⟨X₁⟩ ≈ 10`), and in
  `c₂` units it is `≈ f₁ m`: `≈ 6·10³` at Δ ≈ 300, `≈ 2·10⁴` at Δ ≈ 1100. At J = 0 the series is geometric with ratio `≈ 1/25` through
  `a⁴`. The "`a¹ ≈ 2–4·10³` (J = 28), `8·10³` (1100, 100)" of BELOW_FOLD_SPINJ §2 are this term after subtracting its `ν < 1` mean, i.e.
  mostly the shift of the entries' zero at ν ≈ 1.07; they are not `J²/Δ`-enhanced terms.
* A common factor of this kind is a normalisation. With a fitted constant (BELOW_FOLD_SPINJ) it is absorbed differently at each order
  (the fitted ρ moves by `−137 → +203 c₂` from `a¹` to `a²`, §2.1), which is why the window `λ` moved by `+0.05` per order there.
  With the analytic normalisation of §2.1 it no longer moves (§4).

### 1.2 The genuinely J-dependent part: the helicity structure
Split the 2×2 state/reference log-ratio into a **common part** (entry mean) and a **structure part** (entry `(a,c)` minus entry `(0,0)`;
normalisation- and calibration-free; it is what separates `++` from `+−`). Per-order terms at J = 28 (`c₂`, `m²·Re log`):

| ν | 2 | 3 | 4 | 5 | 7 | 8.33 | 10 | 11.66 | 13.33 | 16.66 | 20.83 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| common `a¹` | +6305 | +5643 | +5331 | +4877 | +6220 | +5371 | +4676 | +3565 | (+464, zero) | +12099 | +7381 |
| common `a²` | +282 | +233 | +163 | +18 | +430 | +138 | −95 | −399 | (−1049) | +452 | −770 |
| structure `a¹` (01, 10, 11) | 6, 2, 3 | 7, 4, 3 | 9, 6, 3 | 15, 13, 7 | −28, −39, −15 | −8, −15, −5 | 19, 18, 9 | 99, 116, 49 | 425, 522, 212 | −1289, −1627, −645 | −1143, −1452, −570 |
| structure `a²` | −9, −13, −5 | −4, −7, −2 | 3, 2, 2 | 21, 24, 11 | −45, −59, −22 | −5, −8, −2 | 35, 42, 18 | 116, 145, 59 | 392, 493, 197 | −746, −950, −370 | −472, −602, −233 |

(ν = 6, the entries' zero: common `+318 / −1512`, structure `123–142 / 239–300`.)
* **Common part:** geometric, ratio `0.03–0.07` on ν ≤ 8.3 (outside the zero at 5.5–7); the ratio grows with ν (`≈ 0.1` at 11.7).
* **Structure part:** small in the window (`|a¹|, |a²| ≤ 25 c₂` on ν ≤ 5 and at 8.33; `≤ 60` at 7) but **not convergent**: `a²/a¹ ≈ 1.2–2`
  at ν = 5–11.7, `≈ 0.9` at 13.3, and it grows like `ν⁶–ν⁸`. At the reference (`J²/Δ = 0.04`) the same ratio is `≈ 0.06`. This is the
  `J`-enhanced non-uniformity; at J = 28 (`J²/Δ = 2.3`) the structure series is asymptotic, not convergent, at ν ≳ 5.
* At J = 32 (order `a¹` only): structure `a¹` = 7–43 `c₂` on ν = 2–5, 14–85 at 7–8 (zero at 6: local blow-up). At (1100, 100) (`J²/Δ = 9.1`):
  structure `a¹` = 15–150 (entries 01, 10) and 48–590 (entry 11) on ν = 2–10, 2.8·10³ at 11.7.


## 2. Method: analytic normalisation, bookkeeping, optimal truncation (and what was not achieved)

### 2.1 Analytic normalisation `ρ_an` (`lj_rs_norm.py`) [D + DV]
* The engine integrand (`bg_geom.assemble` of the `bj` pieces, with the `B`-factor, exactly what `bj_ee.images` integrates) is the
  `a`-expansion of the **same holomorphic integrand** as the direct one (w0 EE integrand of the z2/n3 records), written in the scaled
  variables of the real section `D` with every overall `a`-power dropped (log factors, `v′ − u′`, measure, `pref`) and its own branch
  constants. So point-wise on `D` (where both are defined and the series converges absolutely at small scaled radii)
  `f_eng(p; ν)/f_dir(p; ν) = c(Δ, J)` is `p`-independent up to the engine's truncation, with `f_dir` = `bf_core.grav_parts` (n1_pw
  point-wise w0 integrand, sheet `n₀ = +1`, `ψ = α − φ ∈ (0, π)`).
* By HRS_CLOSE Prop. RS (one populated cohomology class, `c_Γ(ν)` common to state and reference) the crystal integrals then obey
  `V_st/V_ref = (a_st/a_ref)⁶ (c_ref/c_st) I_st/I_ref`, i.e. `ρ_an ≡ κ_st/κ_ref = (a_st/a_ref)⁶ c_ref/c_st` (`bf_core` drops `a⁶`).
  **No parameter is fitted.**
* Checks (10 points, 4 entries, 4 ν):

| state | radial scale ×0.3 / ×1 / ×3: max (median) `|r/c − 1|` | `ρ_an` (×0.3 / ×1 / ×3) | ν-dependence |
|---|---|---|---|
| (336, 28), depth 3 | 1.8·10⁻⁴ (1.6·10⁻⁵) / 6.0·10⁻⁴ (6.4·10⁻⁵) / 1.7·10⁻³ (2.6·10⁻⁴) | 8.798991 / 8.798973 / 8.798927 ·10³¹ | `≤ 6·10⁻¹¹` |
| (336, 32), depth 2 | 7.5·10⁻⁴ (4·10⁻⁵) / 3.6·10⁻³ (1.8·10⁻⁴) / 2·10⁻² (2·10⁻³) | 8.79916 / 8.80107 / 8.81918 ·10³¹ | `≤ 6·10⁻¹¹` |
| references (2·10⁴, 28 / 32) | `≤ 10⁻⁵` | — | — |

  The spread shrinks with the radial scale as a truncation effect must. The ×0.3 value is the truncation-free limit to `≲ 2·10⁻⁶`
  (J = 28) / `≲ 2·10⁻⁵` (J = 32, depth 2). Production: the ×1 value for the J = 28 runs (0.2 `c₂` from the ×0.3 one) and the J = 32
  `a¹` runs (**24 `c₂`** from it, i.e. `≈ 0.004` in `λ/N_s`); the ×0.3 value for every later run (J = 32 `a²`, structure stresses).
* **Against the fitted constant** (bj_chain.calibrate, ±295 / ±540 `c₂` jackknife): `m²(ρ_fit/ρ_an − 1)` = **+203 − 115i** (J = 28, `a²`),
  **−137 − 114i** (J = 28, `a¹`), **−120** (J = 32, `a¹`). The analytic constant lies inside the fit's error at both orders; the fit
  moves by 340 `c₂` between orders, the analytic one does not move (it is order-independent by construction).

### 2.2 Bookkeeping at J = 28 (`lj_rs_book.py`): Prop. RS applies [DV]
On the full J = 28 engine integrand (`a²`, depth 3, N = 4; 5.4·10⁶ monomials per entry, all four entries):
* **B1** `E₁` (the ν-part of the leading exponents) is identical to the J = 0 EE2 one: `max|ΔE₁| = 0`;
* **B2** `E₀ − E₀(J = 0)` = 0 exactly (same lattice point);
* **B3** no negative power of `s, sf, ca, sa` (0 of 5.4·10⁶ monomials);
* **B4** the parity rules hold up to round-off (largest violating `|c|/max|c|` = `9·10⁻¹¹`).
So the J-resummed contraction (`D̂`-series, closed-form `B_k`) does not leave HRS_CLOSE's local system: the J = 28 engine is, order by
order, `c_Γ(ν)⁻¹ ×` the crystal-cycle integral, with the same `c_Γ` as the scalar and the reference, **below the fold as well**. This
answers the review's "test region disjoint from use region" objection (RU issue 2, second problem) at the level of Prop. RS:
[D given RD]. Not re-checked at J = 28: B5 (non-resonance; inherited, the exponents are identical) and B6 (no boundary layer) for the
new J-dependent factors (`√(V²W²)⁻¹`, `D̂`).

### 2.3 Resummation: what was attempted and what was achieved
* **Common part**: needs no resummation. Its large `a¹` is the J-independent normalisation of §1.1; the series is geometric with ratio
  `≤ 0.07` in the window (`≈ 1/25` at J = 0 through `a⁴`). An exponentiated (cumulant) form changes the `a²` image by
  `≈ ½(P₁/P₀)²m² ≈ 5 c₂` at (336, 28): irrelevant.
* **Structure part**: the `J`-enhanced terms are here. At J = 28 their series is asymptotic (ratio `≈ 1–2` at ν ≳ 5); the brief's
  suggested resummations (exponentiating a phase drift, or expanding at fixed `J²/Δ`) were not carried out: the J-dependence of the
  geometric coefficients enters through the conserved coefficients and the closed-form kernel in a way that is not a single phase, and a
  fixed-`J²/Δ` reorganisation needs all geometric orders of the `J²`-leading terms, which the engine does not produce. **A resummed
  engine was not obtained [I].**
* **What is used instead: optimal truncation with an error estimate, validated parameter-free against the direct nodes (§3).**
  Production: `a²` for both parts (the structure terms are `≤ 25 c₂` in the window); error band: `±1` and `±3.3 ×` the last computed
  term (whole image), plus separate structure-only stresses at the size of the unexplained structure offsets of §3.


## 3. Validation against the direct nodes, entry-resolved, no fitted constant (`lj_rs_val.py`) [NE]

Per node and common-seed replica (CV estimator), `m²Re[(V_st/V_ref)_ac/(ρ_an I_st/I_ref)_ac − 1]`; jackknife over the 8 replicas.
Errors of the **structure** combinations (entry − entry (0,0)) are computed on the differences themselves: they are `±1–9 c₂`, 20–500×
smaller than the common errors, because the MC noise is almost entirely common to the four entries (the review quoted per-entry errors
only, RU issue 2). `cache/lj_rs_val_J28_d3_N4_K2.txt`, `cache/lj_rs_val_J32.txt`.

### 3.1 (336, 28): residual direct − engine at each engine order
Common part (entry mean; `R_k` = residual with the engine truncated at `aᵏ`):

| ν | 8.33 | 10.00 | 11.66 | 13.33 | 16.66 | 20.83 |
|---|---|---|---|---|---|---|
| `R₁` (`a¹`) | −1512 ± 1466 | −373 ± 401 | **−75 ± 168** | **−208 ± 193** | **+229 ± 151** | (+1232 ± 18) |
| `R₂` (`a²`) | −1650 ± 1466 | −278 ± 401 | +324 ± 168 | +841 ± 193 | −223 ± 151 | (+2002 ± 18) |
| last term `a²` (common) | +138 | −95 | −399 | −1049 | +452 | −770 |

(20.83: the `D̂`-series is not converged at N = 4, `|N4 − N5| ≈ 10³`, BELOW_FOLD_SPINJ §2; the node tests nothing.) The engine with
**no free constant** reproduces the direct nodes on 10–16.7 within 1.5σ at its optimal truncation (`a¹`), and within the size of the
first omitted term at `a²` (`|R₂| ≤ 0.8 |a² term|` on 11.7–16.7). The common series is asymptotic above ν ≈ 11 (the `a²` term is
already larger than the agreement of `a¹`).

Structure part (entries (0,1), (1,0), (1,1) minus (0,0); `R₀`, `R₁`, `R₂`; first omitted terms `s₁`, `s₂` from §1.2):

| ν | `R₀` | `R₁` | `R₂` (± jk) | `s₁` | `s₂` |
|---|---|---|---|---|---|
| 8.33 | −45, −56, −23 | −37, −41, −18 | −32, −33, −16 (± 4, 5, 9) | −8, −15, −5 | −5, −8, −2 |
| 10.00 | −50, −63, −21 | −69, −81, −30 | −104, −123, −48 (± 3, 3, 6) | +19, +18, +9 | +35, +42, +18 |
| 11.66 | −60, −77, −23 | −159, −193, −72 | −275, −338, −131 (± 2, 2, 4) | +99, +116, +49 | +116, +145, +59 |
| 13.33 | −83, −105, −33 | −508, −627, −245 | −900, −1120, −442 (± 3, 4, 7) | +425, +522, +212 | +392, +493, +197 |
| 16.66 | −269, −337, −148 | +1020, +1290, +497 | +1766, +2240, +867 (± 4, 5, 8) | −1289, −1627, −645 | −746, −950, −370 |
| 20.83 | −425, −537, −217 | +718, +915, +353 | +1190, +1517, +586 (± 1, 1, 2) | −1143, −1452, −570 | −472, −602, −233 |

* The review's "3–16σ entry misses" are the structure part. They are **the engine's own truncation**: at every node ν ≥ 11.7 the
  residual of the order-`aᵏ` engine has the size and the sign of the omitted terms, the order-`a⁰` structure is the closest
  (`|R₀| < |s₁|` at 11.7–20.8), and `R₂/(a¹ − a²)` is the same number `2.2–2.5` for all three structure components and all of
  ν = 11.7–20.8 (a factor 20 in size, two sign changes): the residual tracks the engine's series term by term. That is what an
  asymptotic series at its optimal truncation looks like; it is not a cycle (H_RS) discrepancy, which would not follow the engine's
  own `a²` term.
* **Unexplained**: an offset of `−20 … −60 c₂` (J = 28) at ν = 8.33 and 10.0 that no engine order reproduces (`|R_k| > |s_k|` there),
  3–15σ on the precise structure errors. Same sign and pattern at J = 32 (§3.2), larger: `−40 … −130` at 8.32, `−80 … −220` at 10.0.
  Candidates: a systematic of the below-fold z2 nodes (8.33, 10.0 are below ν* ≈ 11.3), or a beyond-all-orders piece of the structure
  series [I]. It is carried as a structure stress in §4.

### 3.2 (336, 32), engine order `a¹` (ρ_an, no fit)
| ν | 8.32 | 9.99 | 11.65 |
|---|---|---|---|
| common `R₁` | −2602 ± 1659 | −533 ± 598 | +14 ± 400 |
| structure `R₁` (± jk) | −94, −37, −40 (± 6, 8, 14) | −224, −83, −82 (± 5, 7, 12) | −731, −278, −276 (± 5, 7, 13) |
| structure `a¹` term `s₁` | −36, −14, −9 | +58, +24, +30 | +500, +197, +207 |

Common part consistent at every node (≤ 1.6σ). Structure: same pattern as J = 28, `|R₁| ≈ 1.5–4 |s₁|`.
PLACEHOLDER-32D3


## 4. The `ν < 8.3` windows at `Δ_g = 300`: (336, 28) and (336, 32) through the unchanged chain (`lj_rs_chain.py`) [NE + M]

Protocol: exactly BELOW_FOLD_SPINJ §3.1 / `bj_chain.chain2 + chain_win` (two passes of the unchanged `hc_chain.py z4 300`; engine C2
records at every engine node `ν < 8.3`, the k = 2 NL density below 8.3 replaced by the engine's `Γ_NL`; everything above 8.3 is the
HYPOTHESES_CLOSE §5 protocol), **with `ρ := ρ_an` for every replica** (the per-replica reference constant still comes from the
reference records, which the reference engine matches to `10⁻⁷–10⁻⁵`). C3H below 8.3: `keep` (the n3/thimble nodes, as the main run)
or `drop`. Environments: `lj_rs_env28.sh`, `lj_rs_env32.sh` (identical to the `bj_chain` logs). `n = 1.0` in every run (the Δ_g = 300
minimum stays J = 5, `+0.2054`).

### 4.1 Engine values below the fold, (336, 28), `a²`, `ρ_an` (`Re Γ_c`, `c₂`; `cache/lj_rs_nl_J28_a2.pkl`)
| ν | 0.5 | 1.0 | 1.5 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|---|---|
| `++` | +808 | +854 | +689 | +705 | +661 | +586 | +472 | +204 | +892 | +478 |
| `+−` | −2020 | −5647 | +1944 | +401 | −216 | −562 | −1066 | −2656 | +3517 | +150 |
| `++`, `a¹` | +556 | +696 | +354 | +424 | +444 | +449 | +459 | +521 | +66 | +293 |

(ν ≈ 6–7 straddles the zero of `h₋ + h₊`, ν ≈ 1 the zero of `N̂`: the local `Γ_c` swings there, the NL density used by the chain is finite.)
BELOW_FOLD_SPINJ's fitted-constant `a²` values were `++ +903 … +1291`: the difference, `≈ 250–400`, is the fitted constant (+203 `c₂`
of ρ plus its jackknife).

### 4.2 `λ/N_s` per state (`hc_chain_z4_lj_rs*`, gap state is the minimum in every run)

| state | engine | C3H < 8.3 kept / dropped | at `m/Λ` |
|---|---|---|---|
| (336, 28) | **`a²`, ρ_an (production)** | **+0.7592 / +0.7652** | 1.091 (gap) |
| | `a¹`, ρ_an | +0.7568 / +0.7630 | gap |
| | `a² ± 1 ×` last term | +0.7616 / +0.7673 ; +0.7568 / +0.7630 | gap |
| | `a² ± 3.3 ×` last term | +0.7672 / +0.7723 ; +0.7512 / +0.7580 | gap |
| | uniform `Γ_c` shift −2000 / −4000 / −8000 (kept) | +0.431 / +0.103 / −0.553 → **failure at `≈ −4.6·10³`** | gap |
| | for comparison: BELOW_FOLD_SPINJ (fitted ρ), `a¹` / `a²` | +0.739 / +0.787 (kept) | gap / 1.3 |
| (336, 32) | `a¹`, ρ_an (×1) | **+0.8405 / +0.8405** | 1.3 |
| | uniform shift −2000 / −4000 | +0.540 / +0.087 → failure at `≈ −4.4·10³` | gap (1.102) |
| | for comparison: BELOW_FOLD_SPINJ (fitted ρ) | +0.834 | 1.3 |
PLACEHOLDER-32CHAIN
PLACEHOLDER-SSTRESS

**Reading.**
* With the normalisation fixed analytically the window result **no longer depends on the engine order** (J = 28: `a¹ → a²` moves
  `λ` by +0.002, against +0.048 with the fitted constant) and is stable under `±3.3×` the last term (`0.751 … 0.772`).
* The failure scale is an *additional* uniform `Γ_c ≈ −4.5·10³ c₂`; the engine's error in the window is bounded by the last term
  (`≲ 300 c₂` common, `≲ 25 c₂` structure, §1.2), by the unexplained structure offset of §3 (`≲ 10² c₂`, stressed above), and by the
  analytic constant (`≤ 2 c₂`). Margin `≥ 15`.
* For J = 32 the `±1 ×` stress at order `a¹` is not an error estimate: the last term there is the full `a¹` normalisation term
  (`≈ 5·10³ c₂`, §1.1), so `a¹ ± 1×(a¹ − a⁰)` spans `+0.23 … +0.96`; the `a²` run (above) replaces it.




## 8. Reproduce (`ads4/code/spinfunc/`; every job `ulimit -v 3000000; OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1`)

```
# analytic normalisation (needs the cached bj pieces; references: any depth)            ~25 min per J
python3 lj_rs_norm.py 28 336 3 4 2                    # -> cache/lj_rs_norm_J28_D336_d3_N4_K2.{txt,pkl}
BJ_ADAPT=1 BJ_PRUNE=1e-12 BJ_KTOL=1e-6 python3 lj_rs_eeq.py 20000 32 2 4     # reference pieces for J = 32 (3 min) -> cache/lj_rs_bj_*
python3 lj_rs_norm.py 32 336 2 5 1
BJ_PATCH=1 python3 lj_rs_norm.py 100 1100 2 4 1 p ; BJ_PATCH=1 python3 lj_rs_norm.py 64 364 2 4 1 p
# bookkeeping B1-B4 at J = 28 (12 min)
python3 lj_rs_book.py 28 336 3 4 2                    # -> cache/lj_rs_book_J28_D336_d3.txt
# validation, entry/structure resolved, no fit (2 min)
python3 lj_rs_val.py 28 3 4 2 ; python3 lj_rs_val.py 32 2 5 1 ; LJ_VMIN=4 LJ_VTAG=_win python3 lj_rs_val.py 28 3 4 2
# window chains (unchanged hc_chain z4; ~3 min per pass)
./lj_rs_queue28.sh                                    # J = 28: a2, a1, +-1, +-3.3 x last term, uniform shifts -> cache/lj_rs_chain28.log
./lj_rs_queueR.sh                                     # J = 32 reference pieces, norm, val, a1 chains     -> cache/lj_rs_chain32.log
./lj_rs_queue32d3.sh                                  # J = 32 depth 3 (a^2) images (~4 h), val, chains
./lj_rs_queueS.sh                                     # structure stresses                                 -> cache/lj_rs_chainS.log
# high J: (364, 64) depth-2 engine (fw basis), normalisation, fw_gamma-convention Gamma - kappa
./lj_rs_queueH.sh ; ./lj_rs_queueG.sh                 # -> cache/lj_rs_gamma.txt, cache/lj_rs_gamma_J*_K*.pkl
```
Single steps: `LJ_FEXT=f python3 lj_rs_chain.py nl J depth N Kr TAG` (engine `I_K + f(I_K − I_{K−1})`), `LJ_SSTRUCT=d01,d10,d11 …`
(structure stress), then `BJ_C3H=keep|drop [LJ_GSHIFT=g LJ_TAGSFX=_g] python3 lj_rs_chain.py win J TAG` with `source lj_rs_env<J>.sh`.
`LJ_KR=K BJ_NUCUT=ν FW_EST=dir python3 lj_rs_gamma.py J:Δ:depth:N:Nref lo hi`.
