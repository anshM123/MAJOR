# Independent reproduction of the headline numbers (UNIFIED_THEOREM, Results D and Q)

*22 Nov 2026, independent reproducer. Unrefereed, not committed. No existing file was edited and no project script was run or imported.
New files: `ads4/code/repro/bounds.py`, `q_binding.py`, `d_states.py`; outputs in `ads4/code/repro/cache/` (git-ignored by its own
`.gitignore`). Compute: under 2 CPU-min in total, one process at a time, `ulimit -v 3000000`, `OMP_NUM_THREADS=1`. Nothing in
`HLJ_CLOSE.md` or `spinfunc/lj_*` was read or touched.*

## 0. Verdict

* **The headline bound values stand.** Every one of them comes out of the stored functionals by my own exact-arithmetic code, to all
  quoted digits:
  * cubic `3.13·10⁻⁶ / 2.89·10⁻⁸`, `|ĝ₃|²Λ⁸ ≤ 1.07·10⁴ / 1.43·10⁴`;
  * quartic `72.04 / 92.81`, from `A = 17.2407`, `B_eff = −19.2197 / −19.2647`, `κ = 0.10991`;
  * `c₁` and `ε_{c₁}`, `H17`, `δk(0)`, the DG_SCAN asymptotic forms;
  * ABJM `3.34 % / 2.59 %`.
* **Q's binding value `FF/N = +0.1427 / +0.1428` is reproduced (+0.14265 / +0.14282).** My own code does everything downstream of the
  engine: the HYB action, the Legendre transforms `T(ν)`, `Γ_K`, `s_d`, and the `s_d → FF/N` coupling. The input it takes on trust is
  the engine's state/reference ratio `ρ(ν)`, which carries the whole physics of `s_d = +0.51`. That ratio was **not** reproduced.
* **D's binding values are reproduced at the model level.** The hybrid baselines `0.2043` (300, J = 5 `+−`, twist edge) and `0.2009`
  (1000) come out of my own code. The `−0.0005` crystal correction behind the adopted `+0.2038` was not reproduced.
* **No discrepancy beyond stated precision on any headline number.**
* There are nine minor findings in record-keeping and side claims (§2). The most substantive is that DG_SCAN's `s_d ≈ 0.50 + 1440/Δ²`
  does not fit its own anchor points. It errs in the pessimistic direction.

## 1. Table

"Fully" means own code from the definitions plus the stored functional JSONs, which are inputs. "Partly" names what was shared.
Where project code was read to settle a definition, §3 says so.

| quantity | claimed | reproduced | rel. diff | independent? |
|---|---|---|---|---|
| NF02 `k₀` | 27.48509817 | 27.48509817 | < 10⁻⁹ | fully |
| `A = 16ψ₂′(0)` (D) | 3022.34 | 3022.3377 | < 10⁻⁶ | fully |
| `B_AdS` / `B_∞` (D) | −6552.66 / −5327.21 | −6552.662 / −5327.2096 | < 10⁻⁶ | fully |
| `16L_A` (300 / 1000) | 10734.61 / 14339.59 | 10734.61 / 14339.59 | < 10⁻⁶ | fully (exact pole integral) |
| **cubic bound (300)** | **3.13·10⁻⁶** (3.13206·10⁻⁶) | **3.13206·10⁻⁶** | < 10⁻⁵ | fully |
| **cubic bound (1000)** | **2.89·10⁻⁸** (2.89126·10⁻⁸) | **2.89126·10⁻⁸** | < 10⁻⁵ | fully |
| `|ĝ₃|²Λ⁸ ≤` (300 / 1000) | 1.07·10⁴ / 1.43·10⁴ | 10734.6 / 14339.6 (with `δk(0)`: 10742.8 / 14340.6) | — | fully |
| `A_Q = 2ψ₂′(0)` | 17.241 | 17.240713 | < 10⁻⁷ | fully |
| `B_∞^Q = 2∫(ψ₂ − c₁p)/p²` | −19.2782 (stored B) | −19.278176 | < 10⁻⁷ | fully |
| `κ = (−2μ₂,₆ + 2μ₃,₄ + 4a)/8` | 0.1099 | 0.109908 | < 10⁻⁵ | fully (definition read from `common_q.py`) |
| **`2L_AdS` (300 / 1000)** | **72.04 / 92.81** (72.041 / 92.814) | **72.0406 / 92.8135** | 6·10⁻⁶ | fully |
| `B_eff` (300 / 1000) | −19.220 / −19.265 (−19.2647) | −19.2197 / −19.2647 | 2·10⁻⁵ | fully |
| DG_SCAN `b₁, b₂L, b₂` (D) | 19293.4, −13549, −14705 | 19293.4, −13549, −14705 | < 10⁻⁵ | fully (own expansion) |
| DG_SCAN `b₁, b₂L, b₂` (Q); `Λ₀^D`, `Λ₀^Q` | 10.314, 860.3, −2508.6; 8.7415, 4.5889 | 10.3137, 860.35, −2508.6; 8.7415, 4.5889 | < 10⁻⁴ | fully |
| DG_SCAN form errors (cubic, 300 / 1000) | −0.30 % / −0.067 %; +4.5·10⁻⁵ / +3.6·10⁻⁶ | −0.297 % / −0.067 %; +4.45·10⁻⁵ / +3.1·10⁻⁶ | — | fully |
| DG_SCAN form errors (quartic, 300) | −0.058, −0.024, −0.0009 | −0.0585, −0.0239, −0.0010 | — | fully |
| `c₁` (D), per unit | ∓0.22 % / ∓0.017 % (100: ∓2.37 %) | 0.2160 % / 0.0167 % (2.372 %) | — | fully (linearised formula re-derived) |
| `ε_{c₁}` (Q) with D's `t₄` | 0.57 / 0.059 `|c₁|` | 0.5671 / 0.0586 | — | fully |
| H17 band | `±212.5/Λ²`; "±0.2 % / ±0.02 %" | 212.48/Λ² = **±0.24 %** / ±0.021 % | see D2 | fully (heuristic read from `h3_H17_model.py`) |
| EE `δk(0)` first order | −68.08/Λ² | −68.076/Λ² | < 10⁻⁴ | fully |
| **Q binding FF/N (300)**, δ = 0.2, H0 | **+0.1427** (HRS: +0.14265, K = 5) | **+0.14270** (K = 4 from ρ); **+0.14265** with HRS's K = 5 `s_d` | < 10⁻⁴ | **partly**: shares the engine ratio ρ(ν) (BELOW_FOLD2 engine, cached), the de-shift weights χ (Stage 6), and the `R_j` rational functions (Stage 4) |
| **Q binding FF/N (1000)** | **+0.1428** | **+0.14282** | < 10⁻⁴ | partly (same) |
| `FF/N` at `s_d = 0` (HYB, ++ J = 0 gap) | 0.2001 (implied) | 0.200000 (FF = 0.601643, N = 3.008213) | — | fully |
| slope `−2a/N`; `s_fail` | −0.1126; 1.78 | −0.1126; 1.7758 | — | fully |
| `s_d` (K = 4, H0; 300 / 1000) | +0.50924 / +0.50775 (production cache) | +0.50924 / +0.50775 | < 10⁻⁶ | partly (given ρ, χ, `R_j`); own `T`, `P_{−½+iμ}`, `Γ_K` assembly agree with the cache to 3·10⁻¹³ |
| `s_d` by ν window (300) | +0.10, +0.57, −0.05, −0.10, −0.01 | +0.099, +0.569, −0.045, −0.101, −0.011 | — | partly (same) |
| weight shares of `|χ_d|Σ|R_j|` | 24.5 / 43.9 / 9.2 / 11.4 % | 24.5 / 43.9 / 9.2 / 11.4 % | — | partly (χ, `R_j`) |
| per-harmonic tolerances δ = 0.2 | 30.2 / 16.9 / 80.2 / 65.2 | 30.1 / 16.8 / 80.0 / 65.0 | ≤ 0.6 % | partly (χ, `R_j`) |
| flat Q functional, min FF/N | δ = 0.2 (LP margin) | **0.1966** (`+−`, `u = 0.0253`, J = 26); 0.19990 (`++` J = 2, `u ≈ 0.79`) | see D9 | fully |
| flat NF02 min `F/N` | 0.200 [certified] | 0.20000 (`u = 1`, `++` J = 4 / `+−` J = 5; grid) | — | fully (grid, not a certificate) |
| **D binding, 300** (J = 5 `+−`, twist edge 1.0166Λ) | **+0.2038** (hybrid +0.2043, then −0.0005 correction) | hybrid **+0.20427** at 1.01658Λ | +0.2 % vs 0.2038; 0 vs 0.2043 | fully for the hybrid; correction **not reproduced** |
| **D binding, 1000** (J = 5 `+−`, twist edge) | **+0.2009** at "1.0048Λ" | hybrid **+0.20090** at 1.00499Λ (+0.20086 at 1.0048) | < 10⁻³ | fully (hybrid; correction ≈ 0 per statement) |
| D, 1000: J = 4 / 6 / 7 / 11 / 13 | +0.2041 / +0.2296 / +0.2684 / +0.3482 / +0.4076 (corrected) | hybrid +0.2041 / +0.2295 / +0.2685 / +0.3492 / +0.4086 | ≤ 0.3 % | fully (hybrid only) |
| D, 300: J = 4 / 6 / 9 / 10 / 14 | +0.2137 / +0.3075 / +0.4506 / +0.6033 [2.0] / +0.7487 (corrected) | hybrid +0.2136 / +0.3092 / +0.4547 / +0.6047 [2.0] / +0.7499 | ≤ 1 % | fully (hybrid only) |
| ABJM `ML = 2^{5/4}π^{1/2}λ^{1/4}` | 4.22 λ^{1/4} | 4.21563 λ^{1/4} (from `M² = 4/α′`, `L²/α′ = π√(2λ)`) | — | fully |
| `g₄/(8πG) = 2ζ(3)/M⁶` | 2.404 | Virasoro–Shapiro expanded numerically: `VS + 1/(stu) → −2.404115` | — | fully |
| **ABJM placement** | **3.3 % / 2.6 %** (0.0334 / 0.0259) | **0.03337 / 0.02590** | < 10⁻³ | fully |
| `8πG = 256/c_T` | from `c_T = 32L²/(πG)` | `C_T^{Einstein}/C_T^{free scalar} = 32L²/(πG)` from the standard d = 3 formulas; then 256 | — | fully (formulas from the literature) |

## 2. Discrepancies and findings

None of these touches a headline number. They are ordered by how much they matter.

* **D1. DG_SCAN's law `s_d ≈ 0.500 + 1440/Δ²` (also in the UNIFIED update box) misses its own anchors.** [moderate; side claim]
  * It gives `0.516` (300) and `0.5014` (1000). The measured values are `0.5096` and `0.5077`: misses of `+0.007` and `−0.006`.
  * The second miss is 4× the stated `±0.0015` truncation at 1000.
  * A per-order `α_k + β_k/Δ` fit of the cached per-order contributions gives `s_d ≈ 0.507 + 0.6/Δ`. A two-point
    `s_∞ + C/Δ²` fit gives `0.5076 + 147/Δ²` (K = 4) or `0.5075 + 188/Δ²` (K = 5).
  * Likely cause: fitting the large, cancelling per-order terms (`+1.66 / −3.95 / +2.52 / …`) from rounded inputs. The rounded inputs
    sum to `0.516` at 300, while the unrounded cache values sum to `0.5092`.
  * The qualitative conclusion ("Δ_g-stable") holds more strongly than stated. The derived margins `×2.8 / 3.2` at Δ_g = 100 / 150
    are pessimistic: about `×3.45` with the corrected law. **Recommendation:** replace the law, or drop the coefficient.
* **D2. H17 band at 300.** The formula `212.5/Λ²` gives `±0.24 %`, and DG_SCAN itself says `±0.24 %`. The theorem and RESULT_D §0/§6.1
  say `±0.2 %`, which understates it slightly. At 1000, `±0.021 %` is correctly given as `±0.02 %`. **Recommendation:** write `±0.24 %`.
* **D3. Twist-edge mass of the 1000 binding state.**
  * `m_gap(J = 5) = √(998·1009)/998.5 = 1.00499Λ`, not "1.0048Λ" (RESULT_D §0, §5.2).
  * The value is unaffected: hybrid `0.20090` at 1.00499Λ, `0.20086` at 1.0048Λ, both `0.2009`.
* **D4. Result D, J = 15 at Δ_g = 1000.** The quoted `+0.869` is not a minimum over the physical mass range.
  * At its twist edge `1.0149Λ`, my hybrid baseline is `+0.695`. Elsewhere it is `+0.878` at 1.2Λ and `+0.948` at 1.05Λ.
  * The statement does flag "no twist-edge mass" and lists J = 15 as an `H_cont` exception, but the value column reads as a minimum.
  * Not binding (`0.695 ≫ 0.2009`). **Recommendation:** quote `≈ +0.69` (model), or mark the entry as "grid from 1.02Λ".
* **D5. Result D, J = 16 at Δ_g = 300.** The quoted `+0.884 ± 0.022` is above the hybrid `++` baseline at large mass.
  * That baseline is `+0.805` at 3.17Λ.
  * At 300 the `m ≥ 1.3Λ` corrections are negative (`−0.04 … −0.49`, H_R). So the true per-state minimum is plausibly `≤ 0.80`. The
    quoted value appears to be the near-gap minimum only.
  * The same `++` large-m dip gives `0.840` at 3.54Λ (J = 18), consistent with the quoted `+0.872 [4.0]`. Not binding.
* **D6. "Λ⁴ → Δ_g⁴ moves the cubic bound by −2.2 % at 300"** (DG_SCAN §3). I get `(298.5/300)⁴ − 1 = −1.99 %`.
* **D7. `δk(0)` at Δ_g = 100.** RESULT_D's `δk(0) = −7.245·10⁻³` gives `Λ²δk(0) = −70.29`. DG_SCAN quotes "exact −69.9" for the same
  quantity. At 300 and 1000 the two agree (`−68.31`, `−68.10`). One of the two at 100 is mis-transcribed. Δ_g = 100 is not claimed.
* **D8. The `Γ_K` cache reaches `|Γ_K| ~ 10¹²–10¹³` (raw units) near ν = 40.** My recomputation agrees with it to `≤ 10⁻⁹ c₂` per
  harmonic. Its weight there is nil (0.1 % of `|χ_d|Σ|R_j|` on [18, 41)). Noted only because the raw `max |ΔΓ_K| ≈ 10²–10³` looks
  alarming until it is expressed per harmonic.
* **D9. The δ = 0.2 functional `Q3-full-delta0.2` does not meet its own LP margin everywhere in the flat limit.** [minor]
  * On a dense flat grid, using `common_q`'s action and normalisation, min FF/N is `0.1966` at `+−`, `u = 0.0253` (`m ≈ 6.3M`),
    J = 26 (`b̂ ≈ 8.5`). It is `0.19990` at `++` J = 2, `u ≈ 0.79`.
  * This is LP-grid discretisation. Positivity is unaffected (margin ≈ 0.197), and in AdS such states are far above the binding
    `0.1427`. But "δ = 0.2" is nominal at the 2 % level.
  * For comparison, NF02 meets 0.200 on my grid, consistent with its certificate.

Things checked and found consistent, beyond the table:
* `c_T = 32L²/(πG)`, hence `8πG = 256/c_T`.
* `g₄M⁶/(8πG) = 2ζ(3)` is `λ`-independent, so `c_T` cannot enter the placement. The `∂²_m F` flavour-current caveat is correctly
  separated.
* The record correction "∓2.37 % per unit `c₁` at Δ_g = 100" is right.
* "log form underestimates by 0.6 % at 300" is right (−0.59 %).
* `ε_{c₁}` and the `c₁` relative shift follow from `16πG c₁t₄` integrated with `Λ⁴∫ψ₂`, as stated.
* The binding margin "`s_fail/s_d = 1.78/0.51 ≈ 3.5`" is `1.7758/0.5096 = 3.48`.
* The binding D values are, to `10⁻⁴`, the **flat NF02** values at the twist-edge mass. The leading-Regge curvature terms are
  `cv ≈ 10⁻⁴` (300) and `10⁻⁵` (1000). So the AdS input to D's minimum is entirely the beyond-leading-Regge correction, which here is
  `−0.0005`.

## 3. What was and was not reproduced, and what was shared

**Fully independent** (own code; inputs are the stored functional coefficients):
* all light-side constants and bounds of both results;
* the DG_SCAN closed forms;
* the `c₁`, `ε_{c₁}`, H17 and `δk(0)` terms;
* ABJM;
* the flat margins of NF02 and Q3;
* the Q HYB action at the binding state;
* the D hybrid per-state values. My Legendre `P_{−½+iν}(cosh β)` uses a Mehler–Dirichlet quadrature, checked against `mpmath.legenp`
  to `10⁻¹¹`.

**Definitions taken from project code (read only, nothing imported):**
* `κ = k/8`, `k = −2μ₂,₆ + 2μ₃,₄ + 4a`: `resultq/common_q.py`, `lp_g4.py`. The statements give the value, not the formula.
* The flat Stage-0 action, its normalisation `N`, and `D_{+−} = P^{(0,8)}_{J−4}` (stripped): `common_q.py`, `robust/rcommon.py`.
* The HYB model, and the coupling `FF → FF − 2a u³ s_d − c(2u³/Λ²)s_v` with `N` unchanged: `resultq/stage6_model.py`,
  `spinfunc/bg_qan.py`.
* `Γ_K = m²(c T_ref ρ − T_st) − 0.75Σ_j R_j(ν−ij)² − const·N̂₄`, with `Δ_ref = 2·10⁴`: `spinfunc/bg_qan.py`, `bf_qan.py`.
* D's hybrid state, i.e. `cv₂`, `cv₃`, `N_s`, the `+−` stripping `Dt = D/(1−y)⁴`: `spinfunc/ee1000_common.hyb`,
  `omegau/nf_common.flat_parts_nf`, `omegau/common.Dh_all`.
* The H17 heuristic: `spinfunc/h3_H17_model.py`.

**Shared data, not reproduced:**
* **The engine ratio `ρ(ν) = I_st/I_ref`** (BELOW_FOLD2's 7-variable analytic-radial `1/Δ` engine; cached in
  `spinfunc/cache/bg_qan_D{300,1000}_J0.pkl`). It is the whole physics of `s_d`:
  * with the pure leading-Regge ratio `T_st/T_ref`, `s_d` would be `−0.750`;
  * the engine's order 0 alone gives `+1.662`, from ratio deviations of only `~10⁻³` from leading Regge, amplified by `m² ≈ 9·10⁴`.
  * Reimplementing the engine (`~3·10⁵` monomials per order) was out of scope. `H_RS` and premise RD are therefore not tested here.
* The de-shift weights `χ_d`, `χ_v` (H0 min-norm, Stage 6; `resultq/cache/stage6_gauge.pkl`). These are Abel-equation solutions; I did
  not re-solve them.
* The rational `R_j(ν)` (`resultq/cache/stage4_hj.pkl`, EE3). I parsed the sympy strings and evaluated them myself.
* HRS_CLOSE's K = 5 `s_d` (`+0.5096` / `+0.5077`), which is used only for the second FF/N line. My own K = 4 result already gives
  `0.14270` / `0.14282`.

**Not reproduced at all:**
* the `−0.0005` (300) beyond-Regge correction to D's binding state, and every corrected per-state value (they need the near-gap crystal
  chains);
* any `H_*` hypothesis;
* the flat NF02 certificate (only a grid check);
* the Q full-set minimum over the 14 000-state Stage-6 scan (only the binding state);
* the `c_L` slope.

## 4. Reproduce

```
cd ads4/code/repro
ulimit -v 3000000; export OMP_NUM_THREADS=1
python3 bounds.py      # < 1 s  -> cache/bounds.txt     (constants, bounds, DG_SCAN forms, c1, H17, dk(0), ABJM, c_T)
python3 q_binding.py   # ~10 s  -> cache/q_binding.txt  (flat Q margin, HYB at the binding state, s_d from cached rho, FF/N)
python3 d_states.py    # ~45 s  -> cache/d_states.txt   (Legendre check, hybrid per-state values 300/1000, flat NF02 margin)
```
