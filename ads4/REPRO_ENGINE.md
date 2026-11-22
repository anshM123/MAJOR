# Independent reproduction of the engine ratio ρ(ν) behind Result Q's binding `s_d`

*22 Nov 2026, independent reproducer (engine). Unrefereed, not committed. No existing file was edited. New code is in
`ads4/code/repro_engine/`; outputs are in `ads4/code/repro_engine/cache/`, which is git-ignored by its own `.gitignore`. Compute was
≈ 1.6 CPU-h, at most 2 processes at a time, each under `ulimit -v 3000000` and `OMP_NUM_THREADS=1`. Nothing in `HLJ_CLOSE.md` or
`spinfunc/lj_*` was read or touched.*

Grades: **[P]** proved; **[D]** derived; **[DV]** derived and checked; **[NE]** numerical evidence; **[M]** model; **[I]** open.

## 0. Verdict

* **The engine ratio ρ(ν) = I_st/I_ref is reproduced by a different route, and so is the binding `s_d`.**
  * **Route.** I evaluated the EE `k = 4` image directly at finite Δ, with no 1/Δ expansion, no monomial algebra and no closed-form
    Γ/Beta integrals. The tool is a seven-dimensional tensor product-integration rule (own code) applied to my own point-wise
    implementation of the exact integrand on the real section D.
  * **Agreement in ρ.** On ν < 4, ρ_direct agrees with the engine's K = 5 ratio to ≤ 3.4·10⁻⁷ (Δ_g = 300) and ≤ 3·10⁻⁸
    (Δ_g = 1000). Both are normalised at ν = 0.02, as the chain uses them.
  * **Agreement per harmonic.** The difference is ≤ 0.009 `c₂` per harmonic on ν < 4. That is the size of the engine's `a⁵`/`a⁶`
    truncation, and 2000× below the δ = 0.2 tolerances (17–80).
* **`s_d` (H0).**
  * Δ_g = 300: direct **+0.50967**, against engine K = 5 +0.50962 (HRS_CLOSE +0.5096) and K = 4 +0.50924 (production, BF2 +0.509).
  * Δ_g = 1000: direct **+0.50776**, against engine +0.50774 / +0.50775.
  * **FF/N (δ = 0.2)**: +0.14265 (300) / +0.14282 (1000), against the claimed +0.1427 / +0.1428.
* **No discrepancy beyond stated precision.** The one difference I found is a record-keeping matter, not an error in any quoted number:
  * production K = 4 differs from both the direct result and K = 5 by up to 0.053 `c₂` per harmonic (ν ≈ 3, Δ_g = 300);
  * that is 5× HRS_CLOSE's quoted "`a⁵` term ≤ 0.0096";
  * it is explained by HRS_CLOSE §3.5's own finding (incomplete `a⁴` pieces in production; see §2).
* **Safety margin.** Moving `s_d` to the fail point 1.776 needs one of the following:
  * a worst-sign relative error of ρ of **2.1·10⁻⁴** (Δ_g = 300) or **1.9·10⁻⁵** (1000) on ν < 4;
  * a uniform relative error of 4.4·10⁻³ / 4.0·10⁻⁴ for ν ≥ 0.5;
  * scaling the whole beyond-leading-Regge content of ρ by **λ = 2.19 / 2.06**.

  The measured direct-vs-engine difference is about **600× below** the first threshold. In `s_d` terms, |Δs_d| ≲ 3·10⁻⁴ against a
  margin of 1.27.
* **What this does *not* test.**
  * Premise RD (crystal cycle ≡ real section D up to `e^{−O(Δ)}`). Both routes integrate over D.
  * The *definition* of the CLS-D EE integrand. I re-implemented it from the project's formulas, so this is a definition-level
    dependence. It was cross-checked point-wise against the project's independent TPS implementation.
  * The downstream data inputs `R_j` and `χ_d`.

  So this upgrades the engine's **numerics and algebra** for the binding state from "single implementation" to **[DV]**. `H_RS`/RD
  keep their grades.

## 1. Table

| quantity | claimed | reproduced | rel. diff | route |
|---|---|---|---|---|
| conserved ⟨TTO_Δ⟩ coefficients, J = 0 (V₁²V₂² : V₁V₂H₁₂ : H₁₂²) | numerical null space (`l5_spin.conserved`) | **closed form** Δ(Δ+2)/(2(Δ²−6Δ+6)) : 2Δ(Δ−4)/(Δ²−6Δ+6) : 1, 1↔2 symmetric | ≤ 1.8·10⁻¹⁴ at Δ = 300, 1000, 2·10⁴ | independent (sympy, Todorov operator) |
| point-wise EE k = 4 integrand on D | `bf_core.grav_parts` | own numpy (`integrand.py`, own traceless-tensor contraction instead of the null-circle Fourier contraction) | log-spread 1.6·10⁻¹¹ (300), 2.2·10⁻¹⁰ (1000), 1·10⁻⁷ (2·10⁴; the project's side loses digits in v′−u′, see §2), log A ≤ 4·10⁻¹² | independent code; shared definition |
| `A^{iν/2}` leading monomial | `E₁ = (−1, 1, 1, −1, −½, 1, −½, −½)`, "no a¹–a³ term" | `A = 4·Π atoms^{2E₁}(1 + O(a⁴))` (ratio 16 per doubling of Δ); constant 4 is Δ-independent, so no ν-dependent phase between state and reference | — | independent |
| ρ(ν)/ρ(0.02), Δ_g = 300, ν ∈ [0.3, 4] | engine K = 5 (HRS_CLOSE) | direct | ≤ 3.4·10⁻⁷ (ν < 4); ≤ 3·10⁻⁶ (ν < 12) | **independent** (direct integral) |
| ρ(ν)/ρ(0.02), Δ_g = 1000, ν ∈ [0.3, 4] | engine K = 5 | direct | ≤ 3·10⁻⁸ (ν < 4); ≤ 3·10⁻⁷ (ν < 12) | independent |
| Γ_K/Σ\|R_j\| at ν = 0.3/0.5/0.9/1.1/2.0/3.0/4.0/6.0/8.0 (300) | −0.3/−1.0/−3.6/−5.2/−10.4/−7.6/−2.1/+1.5/+0.2 (BF2 §4.2) | −0.307/−1.031/−3.624/−5.249/−10.447/−7.613/−2.151/+1.459/+0.190 | ≤ 0.009 `c₂` vs K = 5 (ν < 4); ≤ 0.053 vs K = 4 | independent |
| `s_d` by window (300, H0): [0,1) [1,2) [2,3) [3,4) [4,8) | +0.10, +0.57, −0.05, −0.10, −0.01 | +0.0989, +0.5702, −0.0458, −0.1021, −0.0116 | ≤ 1.4·10⁻⁴ abs. vs K = 5 per window | independent |
| **`s_d` (H0), Δ_g = 300** | +0.50924 (K = 4) / **+0.5096** (K = 5) | **+0.50967** (direct for ν < 12 = +0.50807; engine K = 5 tail ν ≥ 12 = +0.00159) | +8·10⁻⁴ vs K = 4; **+1·10⁻⁴** vs K = 5 | independent for ν < 12 (99.7 % of `s_d`); tail shared |
| **`s_d` (H0), Δ_g = 1000** | +0.50775 / **+0.5077** | **+0.50776** (direct +0.50642 + tail +0.00133) | 2·10⁻⁵ | as above |
| `s_d` gauges H2 / H4 (300) | +0.5098 / +0.5081 (HRS K = 5) | +0.50983 / +0.50822 | ≤ 2·10⁻⁴ | as above |
| `s_v` (300 / 1000) | ≈ −13.4 | −13.440 / −13.404 (engine K = 5: −13.435 / −13.401) | 4·10⁻⁴ | as above |
| order 0 alone (300): `s_d` restricted to ν < 4.1 | engine K = 0: +1.7448 (total +1.662) | +1.7448 (state at a → 0 with fixed ⟨TTO⟩(300) coefficients, Richardson in Δ = 10⁶/2·10⁶) | < 10⁻⁴; ≤ 0.028 `c₂` per harmonic | independent (ν < 4.1 only) |
| **FF/N, δ = 0.2, Δ_g = 300** | **+0.1427** (+0.14265, K = 5) | **+0.14265** | < 10⁻⁴ | ρ independent; HYB FF, N from REPRO's own code |
| **FF/N, δ = 0.2, Δ_g = 1000** | **+0.1428** | **+0.14282** | < 10⁻⁴ | as above |
| margin `s_fail/s_d` | 1.78/0.51 ≈ 3.5 | 1.7758/0.50967 = 3.48 | — | — |

## 2. Discrepancies and findings

* **E1. Production K = 4 vs K = 5 per harmonic** [NE, record-keeping, harmless].
  * The direct result sides with the K = 5 images. Production K = 4 differs from both by up to 0.053 `c₂` per harmonic, at ν ≈ 3
    (Δ_g = 300).
  * HRS_CLOSE §4 tabulates the "`a⁵` term" as ≤ 0.0096 on ν < 4. The difference K5 − K4(production) is ~5× larger.
  * Consistent explanation: HRS_CLOSE §3.5 already reports that production lacks the `a⁴` pieces of v′−u′ and R (≤ 0.045 `c₂` on
    [2,4)). So K5 − K4(production) = (`a⁴` completion) + (`a⁵` term).
  * The direct integral confirms that this completion is real and correctly signed. Its `s_d` effect is +0.0004, inside the quoted
    ±0.015.
  * **Recommendation:** in HRS_CLOSE §4's table, label the column "K5 − K4(complete)", or add the 0.045 completion explicitly.
* **E2. Precision of the project's point-wise integrand at Δ = 2·10⁴** [I-level note; not used by any quoted number].
  * `bf_core`/`w0_ee_common` compute v′−u′ as a difference of u′, v′ ~ a⁻⁴, which loses about (1/a)² of relative precision.
  * The effect is 10⁻⁷ at Δ = 2·10⁴ and 10⁻¹¹ at 300. My code uses the exact identity
    v′−u′ = 8(x₄·y)(x₄²+y²)/(m₁₃m₂₄), which follows from x₃ = −x₁ and x₂ = −x₄.
  * The engine (series) is unaffected. Only the BELOW_FOLD QMC used the point-wise form, and its noise was 10³× larger.
* **E3. The EE polarisation contraction cancels catastrophically near s → ±1** [noted for any future direct evaluation].
  * The contraction ⟨F, E⊗E, G⟩ is a residual of terms up to ~10⁹× larger at w = 1 − s² = 0.01.
  * In double precision this costs about 10⁻⁶ relative accuracy at such points.
  * I checked it by re-running in 80-bit precision on a full grid. Γ_K moved by ≤ 6·10⁻⁸ `c₂`; `s_d` was unchanged to 10⁻⁶.
* No disagreement with any quoted number.

## 3. Method (what was done)

### 3.1 Integrand [DV] (`integrand.py`, `tto_conserved.py`, `check_pointwise.py`)
* **Definition (read from `spinfunc/bg_geom.py`, `bg_ee.py`).** The integrand is
  ```
  A^{iν/2}·(−x²)⁻⁶ FP₁FP₂FP₃ N_n^{−½} m₂₄^{¾} m₂₀^{−¼} m₄₀^{−¼} c₁₃^{5/4} c₁₀^{−¼} c₃₀^{−¼} √(X₁X₂)
      · m₁₂⁻⁵ m₃₄^{−7/2} m₃₅^{−3/2} m₄₅^{−3/2} V₀₄₂⁻⁴ (v′−u′) pref R B^{Δ/2}
  ```
  with the stated geometry on D.
* **Re-implementation.** Everything was re-implemented in numpy, point-wise and exact in Δ.
* **⟨TTO⟩ coefficients.** These were re-derived symbolically from conservation (Todorov operator, d = 3) and 1↔2 symmetry. The result
  is a closed form not previously in the project.
* **Spin-2 contraction.** Done by explicit traceless projection and metric contraction. The project uses null-circle Fourier
  coefficients with weights κₙ instead.
* **Cross-check.** Against the project's point-wise builder: constant log offset, spread ≤ 2·10⁻¹⁰ at Δ ≤ 1000, identical signs.
* **Symmetries.** Two symmetries were verified point-wise and used:
  * s → −s averaging;
  * the transverse reflection (s, φ, α) → (−s, −φ, −α), which holds exactly.

  With them, the integrand divided by its real weight X₁^{2.25}X₂^{1.75}Y₁^{0.25}Y₂^{−0.25}e^{−S}w^{1.75}cos^{½}φ sin^{−½}α is
  analytic in all seven variables. The exponents were measured independently and agree with HRS_CLOSE's bookkeeping.

### 3.2 Quadrature [NE, converged] (`quad1d.py`, `direct_rho.py`)
* **Rules.** Product integration in each variable: Gauss nodes of the real weight, and ν-dependent weights ω_k(ν) = ∫W·atom^{iνE₁}ℓ_k.
  The weights were computed with fine tanh-sinh or exp-trapezoid rules, not with closed-form Γ functions. Rule self-test: ≤ 3·10⁻¹²
  for ν ≤ 12. The radial rules lose accuracy like e^{πν/4} beyond that (4·10⁻⁵ at ν = 40, a range not used here).
* **Remainder of A.** The remainder `A/(4Π atoms^{2E₁}) − 1` is O(a⁴), with |ε| ≤ 1.5·10⁻⁴ at the outermost nodes. It is kept to
  second order.
* **Convergence per dimension** (Δ_g = 300, change in Γ_K per harmonic, ν < 4):
  * X₁, X₂ and t = s² are converged at 8/8/6 nodes (changes ≤ 0.004);
  * the radial Y₁, Y₂ dominate: 8 nodes are off by 110 `c₂`, while 16 → 20 → 24 changes by < 5·10⁻⁴;
  * φ, α: 8 nodes are off by 0.6, while 16 → 20 changes by ≤ 0.004.
* **Production grid A** (10,10,16,16,8,16,16) against **grid B** (8,8,18,18,6,20,20):
  * total `s_d` differs by 2·10⁻⁴ (300) and 5·10⁻⁵ (1000);
  * per harmonic: ≤ 0.005 on ν < 4, growing to 0.13 at ν = 10, where the high-ν moments get harder.
* **Precision.** f64 against 80-bit agree (E3).
* **Order 0.** The state's a⁰ image was obtained from direct runs at Δ_geom = 10⁶ and 2·10⁶ with the ⟨TTO⟩(300) coefficients,
  Richardson-extrapolated. It reproduces the engine's K = 0 Γ_K on ν < 4.1 to 0.03 `c₂`. This confirms the size of the large order-0
  piece (+1.66 in total), from which the beyond-order-0 cancellation (−3.95, +2.52, …) starts.

### 3.3 Chain to s_d and FF/N (`sdchain.py`; definitions as in REPRO.md)
* Γ_exact = m²(c T_ref ρ − T_st), with c fixed at ν = 0.02.
* Γ_K = Γ_exact − 0.75Σ_j R_j(ν−ij)² − const·N̂₄.
* s_d = 0.04Σχ_d Re Γ_K.
* FF/N = (FF − 2a s_d − c(2/Λ²)s_v)/N, with FF, N from REPRO's own HYB code (`repro/q_binding.hyb_J0`).
* The direct ρ is used for ν < 12. Above that the engine's K = 5 ρ is spliced in; it contributes +0.0016 (300) / +0.0013 (1000).

## 4. What was shared with the project

* **Shared definitions (read, not imported):**
  * the CLS-D EE `k = 4` integrand on D, with its geometry, factor list, exponents, Lemma-B R and the EE prefactor recipe
    (`bg_geom.py`, `bg_ee.py`, `w0_ee_common.py`, `l5_spin.py`);
  * the real section D and the reference Δ_ref = 2·10⁴;
  * the Γ_K / s_d chain (`bg_qan.py`, `bf_qan.py`).
* **Shared data:**
  * `R_j` (`resultq/cache/stage4_hj.pkl`, EE3);
  * `χ_d`, `χ_v` (`resultq/cache/stage6_gauge.pkl`);
  * the engine's ρ for ν ≥ 12 (`spinfunc/cache/rs2_eeq_D*_J0_K5.pkl`), and K ≤ 4 (`bg_qan_D*_J0.pkl`) for comparison only.
* **Project code run** (read-only, comparison only): `bf_core.grav_parts` and `l5_spin.conserved` in `check_pointwise.py`.
* **Other reproducer's code:** REPRO's `q_binding.hyb_J0` (HYB FF, N; independent of the project).
* **Not shared:**
  * the series algebra (`bg_ser`), the a-expansion, the monomial bookkeeping and the closed-form integrals;
  * the ⟨TTO⟩ null-space solver;
  * the Fourier contraction;
  * all quadrature.

## 5. Sensitivity of FF/N to ρ (`sensitivity.py`)

| ρ error model (normalised at ν = 0.02; H0) | Δ_g = 300: error that brings `s_d` to 1.776 | Δ_g = 1000 | measured direct vs engine |
|---|---|---|---|
| worst-sign relative error, \|ε(ν)\| ≤ e₀ on ν < 4 | **e₀ = 2.1·10⁻⁴** | **1.9·10⁻⁵** | 3.4·10⁻⁷ / 3·10⁻⁸ (×600 margin) |
| same, on ν ∈ [1, 2) only | 4.6·10⁻⁴ | 4.2·10⁻⁵ | 1·10⁻⁷ / 1·10⁻⁸ |
| same, all ν < 41 | 2.0·10⁻⁴ | 1.8·10⁻⁵ | — |
| uniform ε = e₀ for ν ≥ 0.5 | 4.4·10⁻³ | 4.0·10⁻⁴ | — |
| scale the beyond-leading-Regge content, ρ = ρ_LR(1 + λ(ρ/ρ_LR − 1)) | **λ = 2.19** | **λ = 2.06** | λ − 1 ≈ 1·10⁻⁴ |
| pure phase error of ρ | ~4·10⁶ rad (insensitive: Re Γ_K) | ~1·10⁶ rad | — |

* FF/N moves by −0.1126 per unit `s_d`. FF/N = 0 at s_d = 1.7758, i.e. a margin of 1.266 (300) / 1.268 (1000).
* For scale: ρ deviates from the leading-Regge T_st/T_ref by up to 8.7·10⁻³ (300) / 7.7·10⁻⁴ (1000) on ν < 4. That is ≈ 770/m² at
  both Δ_g. A ≈ 2.5 % error *in that deviation's shape* (worst sign, both Δ_g) would already break δ = 0.2. So would a uniform doubling
  (λ ≈ 2.1).
* **The direct reproduction excludes both at ≥ 600×** [NE] — for the integral over D. RD (crystal cycle ≡ D beyond all orders) is a
  separate premise and is untouched.

## 6. Grades after this reproduction

| item | before | now |
|---|---|---|
| engine algebra and closed forms for ++ J = 0 at the gap (Δ_g = 300, 1000) | single implementation ([DV] only against the scalar skeleton) | **[DV]**: an independent direct integral at finite Δ agrees with K = 5 to the truncation level |
| truncation (K = 5 vs all orders) | [NE] (size of the last term) | **[DV] on D**: direct − K5 = +5·10⁻⁵ (300) / +2·10⁻⁵ (1000) in `s_d` |
| binding `s_d` = +0.510 (300) / +0.508 (1000) | [NE given H_RS] | [DV] on D, given the integrand definition; still conditional on RD for the crystal cycle |
| H_RS / RD | [D-sketch] (HRS_CLOSE) | unchanged (not tested) |

## 7. Reproduce

```
cd ads4/code/repro_engine
ulimit -v 3000000; export OMP_NUM_THREADS=1
python3 tto_conserved.py        # 10 s    -> cache/tto_conserved.txt   (closed-form conserved <TTO> coefficients)
python3 quad1d.py               # 5 s     (1-D product-integration rule self-test)
python3 check_pointwise.py      # 1 s     -> cache/check_pointwise.txt (runs project bf_core/l5_spin read-only; comparison only)
./dimscan.sh; ./yscan.sh; ./ascan.sh   # ~25 min (convergence scans; dimdiff.py prints them, e.g.
                                       #   python3 dimdiff.py 300 8-8-16-16-6-8-8 8-8-16-16-6-16-16 8-8-16-16-6-20-20)
./final_queue.sh                # ~35 min: grids A, B at Dg = 300, 1000, ref 2e4 (nu < 12); 80-bit check (nu < 4.1)
./order0_queue.sh               # ~11 min: order-0 diagnostic (Delta_geom = 1e6, 2e6, <TTO>(300))
python3 summary.py              # 30 s    -> cache/summary.txt   (tables of sections 1, 2, 5)
python3 sensitivity.py          # 20 s    -> cache/sensitivity.txt
python3 compare.py 300 10,10,16,16,8,16,16 f64   # optional detailed per-nu comparison -> cache/compare.txt
```
