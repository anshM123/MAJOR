# Below the fold, part 2: the analytic-radial 1/Δ expansion of the crystal integral, its validation, and a first measurement of the `++` EE `k = 4` correction

*21 Nov 2026, below-fold-2 agent. Unrefereed, not committed. No existing document or code file was edited; nothing of the `oo_*`
agent was touched or run. New code: `ads4/code/spinfunc/bg_*.py`, `bg_qrun.sh`. Caches `spinfunc/cache/bg_*` (git-ignored).
Read-only imports: `bf_core` (only for the pointwise checks), `bf_qan` (`mass`, `Rj`, `T`, `weights`), `bf_q` (QMC caches), `bf_exact`,
`n3_exact`, `l5_spin.conserved`, `l6_spin_action.shadow` (probed only), `l1_crystal.laplace` (scalar K₂ check), `resultq/stage4_hj`,
`stage6_model`, `stage6_stress`, `stage6_scan`, `common_q`. Compute ≈ 1.5 CPU-h, ≤ 2 cores, ≤ 3 GB (`ulimit -v`).*

Grades: **[P]** proved; **[D]** derived; **[DV]** derived and checked numerically; **[NE]** numerical evidence; **[M]** model; **[I]** open.
"c₂" = `m² ×` relative correction (Result D / Stage 7 convention). `s_d` = `0.04 Σ_ν χ_d(ν) Re Γ_K(ν)` (Stage 7 §2).

---

## 0. Answer in brief

1. **Construction (task 1) [DV].** BELOW_FOLD §5's design works, and better than proposed: on the real section `D` every factor of the
   CLS-D integrand (scalar and EE) is `a^d × (atom monomial) × (1 + ε)` with `ε` an **exact polynomial in `a`** (degree ≤ 4) in the
   radial variables, and the only non-monomial leading factors (`v′ − u′ ∝ X₁Y₂ + X₂Y₁ − 2X⊥Y⊥`, the EE prefactor, Lemma-B `R`,
   `P_J`) enter with integer powers. Hence every order in `a = 1/Δ` is a finite sum of monomials in
   `(√X₁, √X₂, √Y₁, √Y₂, s, 1−s², cos φ, sin φ, cos α, sin α, 1∓cos α)`, and **all seven integrals are closed form**
   (four radial Γ's, three Beta-type angular integrals; §1.3). No quadrature, no Monte Carlo, no fold, no `e^{−πν/2}` loss, no boost
   tails (they become Γ-function shifts). The EE polarisation prefactor (`w0_ee_common`: `<TTO>` conserved structures, null-circle
   Fourier contraction, `E ⊗ E`) is re-built term by term on this algebra (`bg_ee.py`) and **agrees pointwise with the project's exact
   integrand** (`bf_core.grav_parts`) to `0.11, 1.1·10⁻², 6.3·10⁻⁴, 3.6·10⁻⁵, 1.6·10⁻⁶` (orders a⁰…a⁴, Δ = 300) [DV].
2. **Validation (task 2).**
   * **(a) scalar vs the exact Mellin skeleton [DV/NE]:** J = 0, 2, 4 at Δ = 100, 300, 1000: `|m²ε| ≤ 0.03` (Δ ≥ 300) and `≤ 1.5`
     (Δ = 100) on **all of ν ∈ [0, 18]**; Δ = 400: `≤ 0.03` to ν = 40. Large spin needs the `J/Δ` resummation of `P_J` (§2.2):
     (336, 24/28/32), (300, 24), (1000, 24): `≤ 0.04` on ν ≤ 10, `≤ 0.23` to ν = 18; (100, 24): `≤ 1.2` (ν ≤ 4). BELOW_FOLD's `D`-QMC
     reached ν ≲ 2.6 (J > 0) with ±10–300; this reproduces it and extends it to every ν tested.
   * **truncation per order [NE]:** scalar Δ = 300: the a⁴ term is ≤ 0.15 c₂ (ν ≤ 5), 2.9 (ν = 18); a⁵: 0.18; a⁶: 0.01. EE k = 4,
     Δ = 300, J = 0: a⁴ term ≤ 0.4 (ν ≤ 5), 1.8 (ν = 8), ≤ 4 (ν = 18–30) except at the zeros of `N̂₄` (ill-conditioned log, 8 at
     13.1); Δ = 1000: ≤ 1 everywhere. Projected: the a⁴ term moves `s_d` by 0.016 (Δ = 300) / 0.002 (Δ = 1000). **Target
     (≲ 1 c₂ after a⁴ at Δ = 300) met for ν ≲ 18** (successive orders fall by ×12–17).
   * **(b) EE:** the reference (Δ = 2·10⁴) EE k = 4 image reproduces all four real zeros of Stage 4's exact `N̂₄^{EE}`
     (1.068, 4.546, 8.169, 13.094) with `N = I_ref/T_ref` smooth through them; k = 2 reproduces EE1's zeros [DV]. Consistent with the
     BELOW_FOLD QMC on `D` (χ²/n = 0.6, errors 150–2500) [NE, weak]. The K₂ Laplace rows (`w1_pts_EE_J0`, Δ = 400, ν = 20 → 40) differ
     by 512 c₂; the same Laplace machinery on the **scalar** is off by −107 c₂ there while this method is exact to 0.03 — attributed to the
     Laplace truncation, not proven [NE]. n3/n4 caches are C2/C3H kinds (J ≥ 4), not EE: not compared [I].
   * **(c) scalar large J → leading Regge [NE]:** at fixed `b̂ = (2J+1)/m` = 0.15, 0.3 the deviation `m² log(R/R_LR)` from Stage 4's
     `Σ_j h_j P_{−½+i(ν−ij)}(cosh β)` is a finite function of (ν, b̂), equal to 1–3 % for m = 333 → 2990 (J up to 450).
3. **Result Q (task 3) [NE + M].** `++` measured; `+−` not (§4.4).

   | state | Δ_g | `s_d` H0 | H2 | H3 | H4 | SB | XP | FF/N δ = 0.05 (K=0 → meas.) | δ = 0.2 |
   |---|---|---|---|---|---|---|---|---|---|
   | `++` J = 0, m = Λ | 300 | **+0.509 ± 0.016** | +0.510 | +0.509 | +0.508 | +0.525 | +0.515 | +0.050 → **+0.006** (SB +0.005) | +0.200 → **+0.143** (SB +0.141) |
   | `++` J = 0, m = Λ | 1000 | +0.508 ± 0.002 | +0.508 | +0.508 | +0.507 | +0.523 | +0.513 | +0.050 → +0.006 | +0.200 → +0.143 |
   | `++` J = 0, m = 1.05Λ | 300 | +0.509 ± 0.014 | +0.510 | +0.509 | +0.508 | +0.525 | +0.515 | +0.148 → +0.112 | +0.289 → +0.244 |
   | `++` J = 2, m = 1.05Λ | 300 | +0.511 ± 0.010 | +0.511 | +0.511 | +0.510 | +0.526 | +0.516 | +0.061 → +0.027 | +0.221 → +0.181 |
   | `++` J = 2, m = 1.05Λ | 1000 | +0.508 ± 0.001 | +0.508 | +0.508 | +0.507 | +0.523 | +0.513 | +0.061 → +0.027 | +0.221 → +0.181 |

   * errors = size of the a⁴ term (the omitted a⁵ is ≈ 1/15 of it); FF/N errors ±0.0013 (Δ_g = 300), ±0.0001 (1000).
   * **`s_d ≈ +0.51`, the same in every state, both Δ_g, every gauge to ±0.016** (gauge spread H0 ↔ SB = 0.016: the null projections
     on ν ∈ [2.5, 18] are small). It is **not** a covariant renormalisation of the HYB constant (fit: δγ = +0.014; the rest, s_d = 0.495,
     is a per-harmonic `Γ_K ≈ −1 / −5 / −10 / −2 / +1.5` c₂ at ν = 0.5 / 1.1 / 2 / 4 / 6; 70 % of `s_d` from ν ∈ [0.5, 2.5]).
   * **effective K** (`s_d / s`, Stage 7 derivative stress coefficients): structure (a) −0.57 (H0/H3/H4), −0.61 (H2), −0.89 (XP);
     structure (b) +1.38 (H0), −0.45 (SB), −3.0 (XP) — i.e. `|K| ≈ 0.5–0.9` in the (a)/(SB-b) reading, the natural scale.
   * Applying the measured `s_d` to all `++` J = 0, 2 states of Stage 6's set (`bg_ffscan`): min FF/N = **+0.0061 ± 0.0013 (δ = 0.05;
     SB +0.0048, XP +0.0056)** and **+0.1427 ± 0.0017 (δ = 0.2; SB +0.1410)**, Δ_g = 300 and 1000, binding state `++` J = 0 at the gap.
   * **δ verdict.** δ = 0.2: positive with margin +0.14 at every measured state and gauge → **quote δ = 0.2** [NE + M]. δ = 0.05: still
     positive at the measured binding state, but only by +0.005 ± 0.001 (≈ 10 % of δ) and the `+−` J = 4, 5 states (tolerance 4–6 c₂,
     Stage 7 §2.2) are unmeasured → **δ = 0.05 not established** [I].
4. **Result D (task 4) [I].** `HYPOTHESES_CLOSE2.md` does not exist. The **scalar** skeleton at (336, 28/32) is now reproduced to ≤ 0.04 c₂
   on the whole ν < 8.3 window (§2.2), but the EE engine samples the spin-J polarisation on 2J+1 null-circle points with `V₅^J`
   powers; at J = 28–32 (and (1150, 150)) this is out of reach. Needs a closed-form spin-J contraction (§5).

---

## 1. Construction (`bg_ser.py`, `bg_geom.py`, `bg_ee.py`)

### 1.1 Variables and exact geometry [P]
Section `D` of bf_core with `n₀ = +1` and ψ ∈ (−π, π] (equivalent to both sheets `n₀ = ±1`, ψ ∈ (0, π)), α = ψ + φ ∈ (−π, π]:
`x = a(p₁², p₂², s p₁p₂)`, `y = a(q₁², q₂², q₁q₂ sin φ)`, `u₀ = aY₁ cos(ψ−φ)`, `v₀ = −aY₂ cos α`, `y₀ = aρ sin ψ` (X_i = p_i², Y_i = q_i²).
Measure × CLS-D meas-factor × celestial Jacobian = `a⁵ p₁p₂ dX₁dX₂ ds dY₁dY₂ dφ dψ` exactly (the `|sin ψ|` cancels; the integrand is
analytic and 2π-periodic in ψ).

### 1.2 Structure of every factor [P (sympy) + DV]
| factor | leading | relative correction `ε` (exact) |
|---|---|---|
| `−x²`, `m₂₄`, `c₁₃`, `m₁₃`, `c₁₀ = aY₂(1 − cos α)`, `c₃₀ = aY₂(1 + cos α)` | atom monomial | 0 |
| `m₁₂, m₃₄, m₁₄, m₂₃`, FP₁,₂ | `1/(X₁X₂(1−s²))`, 1 | `∓a²C + a⁴ q_x q_y` |
| `m₁₅, m₂₅, m₃₅, m₄₅` | 1 or `1/(X₁X₂(1−s²))` | `a(Y₁+Y₂) − a² q_y` etc. (the only odd orders) |
| `m₂₀, m₄₀`, `N_n` | monomial | a², a⁴ polynomials (A has no a¹–a³ term) |
| `v′ − u′` | `C/(4X₁X₂Y₁Y₂(1−s²)cos²φ)`, `C = X₁Y₂+X₂Y₁−2s√(X₁X₂Y₁Y₂) sin φ` | power 1 |
| `B^{Δ/2}` | `e^{−(X₁+X₂+Y₁+Y₂)}` | `exp(Σ a^{2k} b_k)`, `b_k` polynomials |
| `t − 1` (spin angle) | — | `a²τ₂ + a⁴τ₄ + …`, `τ₂` a degree-2 polynomial (no inverse powers) |
Leading density (J = 0): exponents `(p₁, p₂, q₁, q₂, 1−s², cos φ, 1−cos α, 1+cos α) = (2.5, 1.5, 2.5, 1.5, 0.75, 2.5, −¼, −¼)
+ iν(−1, 1, 1, −1, −½, 1, −½, −½)` times `C` (BELOW_FOLD §1.2's density times `√(X₁X₂)`).

### 1.3 Algebra and integrals
* `bg_ser.Poly`: Laurent polynomial in 12 atoms (redundant, canonicalised by `s² = 1−w`, `sin² = 1 − cos²`); `ASer`: truncated Laurent
  series in `a`; division / logs / non-integer powers only through atom-monomial leading coefficients (recognised numerically by a
  log-linear fit, verified at independent points). `bg_geom.assemble` combines `Σ p_f log f` (with `p_f` possibly `∝ iν` for `A`) and
  the polynomial factors into `atoms^{E₀ + iνE₁} e^{−S} Σ_j (iν/2)^j Σ_k a^k corr_{j,k}`.
* `bg_ser.integrate_fast`: radial `Γ(1 + e/2)`, s: `Γ(½)Γ(G+1)/Γ(G+3/2)`, φ: `Γ(½)Γ((M+1)/2)/Γ(M/2+1)`, α:
  `2^{A+B+1} Σ_j C(r,j)(−2)^j B(A+j+½, B+½)`; each checked against mpmath quadrature to 10⁻¹² [DV].
* `bg_scalar`: the moments `M_{k,n,j}(ν)` are universal (Δ, J enter only through `Δ^{−k}` and `c_n(J) = P_J^{(n)}(1)/n!`); K = 6 in 10 s.
  `bg_scalarJ`: `P_J(1+δ) = Σ_n c_n(J) a^{2n} D̂ⁿ`, kept to n ≤ 12 at residual order R ≤ 4 (fixed J/Δ).
* `bg_ee`: `tto_samples`, `V`, `H`, `e_factor`, `_contract` of l5/w0 re-implemented on `ASer` (same mons/sol from `l5_spin.conserved`,
  same DFT/κ contraction). Pointwise check vs `bf_core.grav_parts` (10 random points, Δ = 300, ν = 0.5, 5): max |Δ log| per order
  0.11 / 1.1·10⁻² / 6.3·10⁻⁴ / 3.6·10⁻⁵ / 1.6·10⁻⁶; Lemma-B `R` to 5·10⁻¹² [DV]. Scalar: 2.6·10⁻⁴ / 2.5·10⁻⁷ / 2.6·10⁻¹⁰ (K = 2/4/6,
  J = 0); J = 4 likewise [DV].
* Cost: EE J = 0 state 2 min (depth 5 = K + 1; the prefactor loses one order to a leading cancellation), J = 2 14 min; ≤ 0.7 GB.

## 2. Scalar validation (`bg_validate.py`, `bg_validateJ.py` → `cache/bg_validate*.txt`)

### 2.1 Fixed order (K = 6), `m²ε = m²[log(R_bg/R_exact) − const]` (const on ν < 1.4)
| state | ν ≤ 3 | 3 < ν ≤ 10 | 10 < ν ≤ 18 | grade |
|---|---|---|---|---|
| Δ = 100, J = 0 / 2 / 4 | 0.000 / 0.003 / 0.03 | 0.16 / 0.22 / 0.54 | 0.15 / 0.03 / 1.5 | [DV] |
| Δ = 300, J = 0 / 2 / 4 | ≤ 0.005 | ≤ 0.01 | ≤ 0.005 | [DV] |
| Δ = 1000, J = 0 / 2 / 4 | ≤ 0.03 | ≤ 0.001 | ≤ 0.001 | [DV] |
| Δ = 400, J = 0 (ν = 20–40) | — | — | ≤ 0.03 (to ν = 40) | [DV] |
| Δ = 336, J = 24 / 28 / 32 (fixed order) | 0.4 / 1.1 / 2.3 | 11 / 28 / 65 | 10³ | fixed order fails at large J |
(the comparison is ill-conditioned at ν ≈ 5.5–6.6, zero of g + g̃; the exact object's own interpolation noise is ~0.04–0.16 at ν ≈ 3.)
Per-order sizes at Δ = 300, J = 0 (c₂, ν = 1 / 3 / 18): a²: 2.2 / 23 / 413; a³: 0.18 / 1.8 / 33; a⁴: 0.01 / 0.07 / 2.9; a⁵: 0 / 0 / 0.18;
a⁶: 0.01. The a¹ term is ν-independent (exactly: ⟨Y₁+Y₂⟩ has no ν-dependence) [DV].

### 2.2 J/Δ resummation (N = 12, R = 4) [DV]
max |m²ε| on ν ≤ 4 / 4–10 (excl. 5–7) / 10–18: (336, 24) 0.04 / 0.01 / 0.23; (336, 28) 0.04 / 0.01 / 0.16; (336, 32) 0.02 / 0.01 / 0.14;
(300, 24) 0.02 / 0.01 / 0.23; (1000, 24) 0.16 / 0.00 / 0.01; (100, 24) 1.2 / 0.13 / 3.8; N = 10 identical (n-series converged); R = 2
gives 3 / 6 / 50–80: the residual order matters, the n-order does not.

### 2.3 Large J → leading Regge (`bg_lrscalar.py`) [NE]
`m² log(R_bg/R_LR) − (ν₀)` with `R_LR = T_st/T_ref` (Stage-4 `sc1` kernel): at b̂ = 0.15, ν = 1 / 4 / 10 / 14: 2.96 / 63.9 / 129 / 151
(m = 333), 2.95 / 63.5 / 128 / 181 (995), 2.89 / 62.2 / 125 / 185 (2990): a finite b̂-function, approach O(1/m²). (Δ = 10⁴ rows drift by
20 %: double-precision limit of `m² log` at m² = 10⁸, not used.)

## 3. EE checks (`bg_lrcheck.py`, `bg_eetrunc.py`, `bg_k2check.py`)

### 3.1 Leading Regge [DV]
Reference image (Δ = 2·10⁴) vs Stage 4's exact kernels: sign changes of `Re I_ref` (phase-fixed) coincide with those of `T_ref` at
1.10, 4.58, 8.18, 13.10 (k = 4, EE3) and 1.18, 6.66 (k = 2, EE1); `|N|` is monotone through every zero (e.g. 3.94 / 3.78 / 3.63 at
1.03 / 1.07 / 1.11); median |Δ² log N| = 1.6·10⁻⁴. (The extra sign changes of Re I_ref are the phase of the ν-only normalisation.)

### 3.2 Order structure of the EE images [NE]
The state's conserved `<TTO>` coefficients differ from the reference's at O(1/Δ); the resulting a⁰ shape difference and the a¹ term are
each O(1) c₂ (not O(Δ)): Δ = 300 / 1000 give the same per-order sizes at k = 1, 2 (e.g. EE4 k = 1: +40.8 / +43.5 at ν = 3), k = 3 falls
like 1/Δ, k = 4 like 1/Δ². Per-order table: `cache/bg_eetrunc.txt`.

### 3.3 Other comparisons
* BELOW_FOLD's QMC on `D` (`bf_q` caches, Δ = 300 / 1000, J = 0, k = 4 and 2): χ²/n = 0.56–0.68 [NE, weak: errors 150–2500 c₂].
* K₂ rows (`w1_pts_EE_J0_k2_b`, Δ = 400): `m²[ΔD(20→40)]` differs by +512 c₂ (ours −599, converged to 2 between a³ and a⁴). The K₂ rows'
  own a⁴-type term is 12 028 c₂ at ν = 20 with no M6; for the scalar, the same Laplace machinery (`l1_crystal.laplace`, K₁/K₂) misses the
  exact ΔD(20→40) by −2 / −107 c₂ while this method is exact to 0.03. So the discrepancy is consistent with Laplace truncation; not proven [NE].
* n3/n4: their caches are C2/C3H kinds at J ≥ 4 (two structures); not compared [I].

## 4. Result Q (`bg_eeq.py`, `bg_qan.py`, `bg_covfit.py`, `bg_ffscan.py`)

### 4.1 Chain [D]
As bf_qan, but exact order by order: `Γ_exact = m²(c T_ref I_st/I_ref − T_st)` (no fit of `N`), `Γ_K = Γ_exact − 0.75 Σ_j R_j(ν−ij)²
− const·N̂₄` (`Γ_K(0.02) = 0`), projections on the whole de-shift grid (ν < 41; the restriction to ν < 18 changes `s_d` by ≤ 0.007).
FF/N at the state: Stage-6 HYB, `FF → FF − a·2u³ s_d − c·(2u³/Λ²) s_v` (`s_v ≈ −13.4`: Λ⁻²-suppressed, a_h-convention dependent).

### 4.2 Result [NE + M]
Table in §0.3. Per-order contributions to `s_d` (H0, Δ = 300, J = 0): +1.66 / −3.95 / +2.52 / +0.27 / +0.016 (k = 0…4); Δ = 1000:
+1.66 / −4.21 / +2.96 / +0.09 / +0.002. Decomposition (`bg_covfit`, fit ν < 12): covariant part δγ = +0.013 … +0.015; the remainder
carries `s_d(rem)` = +0.495 (H0) … +0.510 (SB); per harmonic it is ≤ 0.5 c₂ on ν ∈ [2.5, 18] and up to 5.7 c₂ on ν < 2.5.
The per-harmonic `Γ_K/Σ|R_j|` (Δ = 300, J = 0): −0.3 / −1.0 / −3.6 / −5.2 / −10.4 / −7.6 / −2.1 / +1.5 / +0.2 at ν = 0.3 / 0.5 / 0.9 / 1.1 /
2.0 / 3.0 / 4.0 / 6.0 / 8.0; the cumulative H0 projection peaks at +0.71 (ν = 2.5) and settles at +0.509 by ν ≈ 8.

### 4.3 Reading
* The measured correction is of the size Stage 7 called "natural" (`|K| ≈ 0.6–0.9` in its structures (a)/(b)) and is **almost
  gauge-independent** (spread 0.016 ≪ 0.5 tolerance): the three null projections are small for these states [NE].
* With it, δ = 0.05 survives at the measured states by +0.005 ± 0.001 (`++` J = 0 gap, all gauges, both Δ_g) and δ = 0.2 by +0.14.
  Both statements are [M] for the rest of the heavy action (HYB flat-exact + LR curvature; the EXA rows of Stage 6 are kept as they are)
  and [NE] for the EE correction.

### 4.4 Not measured [I]
`+−` J = 4 (gap 1.013Λ) and J = 5 (1.15–1.2Λ). J = 4 needs the 2×2 prefactor matrix (computable with `bg_eeq`, all four entries are kept)
contracted with the shadow matrix (`l6.shadow`, 16 s) **and a helicity basis**; at the a⁰ order the image is helicity-blind (Stage 4), so the
`++/+−` basis is fixed only by the x = ν/m dependence (order a²) — not attempted. J = 5 is parity-odd (other agent's sector).

## 5. Remaining open items [I]
1. `+−` J = 4, 5 EE k = 4 images (helicity basis at small x; odd structures for J = 5) → the δ = 0.05 decision.
2. Result D below the fold: an EE (and C2/C3H) engine for J = 28–150 needs the spin-J polarisation contraction in closed form (the scalar
   case is `P_J(t)`, resummed here); the sample-based `V₅^J` route scales like J³ in products.
3. An independent EE check at ν ≲ 10 (e.g. n3 contour for EE at J = 0, 2), and the source of the 512 c₂ K₂ discrepancy (§3.3).
4. The real-section ↔ crystal-cycle equivalence is established for the scalar at all J (to 0.04 c₂ with resummation) and for the EE only at
   leading order (zeros of N̂₄) and through the QMC consistency; for the EE beyond leading order it is assumed [D-sketch: the same
   analytic structure, Γ-function continuation of polynomial moments].

## 6. Graded claims
| claim | grade |
|---|---|
| every factor of the scalar and EE integrands on `D` = monomial × (1 + exact polynomial in a); all integrals closed form | [P] (sympy) + [DV] |
| series integrand = exact integrand pointwise to O(a^{K+1}) (scalar K ≤ 6, EE K ≤ 4) | [DV] |
| scalar J ≤ 4: `|m²ε| ≤ 0.03` (Δ ≥ 300), ≤ 1.5 (Δ = 100) on ν ≤ 18; Δ = 400 to ν = 40 | [DV] |
| scalar J = 24–32: resummed in J/Δ, ≤ 0.04 (ν ≤ 10), ≤ 0.23 (ν ≤ 18) | [DV] |
| scalar large-J approach to leading Regge, O(1/m²) | [NE] |
| EE reference image reproduces the four zeros of `N̂₄^{EE}` (and EE1's two) | [DV] |
| EE truncation: a⁴ term ≤ 0.4 c₂ (ν ≤ 5), ≲ 4 (ν ≤ 30) at Δ = 300; `s_d` error 0.016 / 0.002 | [NE] |
| `++` J = 0, 2 near-gap: `s_d = +0.51 ± 0.02` in every gauge, both Δ_g | [NE] (given D ≡ crystal for EE beyond LO) |
| min FF/N with it: +0.006 ± 0.001 (δ = 0.05), +0.143 ± 0.002 (δ = 0.2) | [NE + M] |
| δ = 0.2 recommended; δ = 0.05 not established (`+−` unmeasured, margin 10 %) | [M] / [I] |
| Result D points | [I] |

## 7. Files (`ads4/code/spinfunc/`)
| file | content | run |
|---|---|---|
| `bg_ser.py` | Laurent-polynomial / a-series algebra, atom recognition, closed-form integrals (`integrate_fast`) | library |
| `bg_geom.py` | exact geometry on D, scalar pieces (l1_crystal.build), B-correction, spin angle t, assembler | library |
| `bg_ee.py` | EE prefactor (l5_spin / w0_ee_common term by term), EE logs, Lemma-B R | library |
| `bg_scalar.py`, `bg_scalarJ.py` | scalar moments (fixed order K; J/Δ-resummed N, R) → `cache/bg_scalar_K6.pkl`, `bg_scalarJ_N*_R*.pkl` | 10 s; 1–5 min |
| `bg_validate.py`, `bg_validateJ.py`, `bg_lrscalar.py` | §2 tables → `cache/bg_validate.txt`, `bg_validateJ.txt`, `bg_lrscalar.txt` | 1 min each |
| `bg_eeq.py` | EE k = 4 / 2 images, state + reference, per order on the dν = 0.04 grid → `cache/bg_eeq_D*_J*_K4.pkl` | 2 / 14 min |
| `bg_lrcheck.py`, `bg_eetrunc.py`, `bg_k2check.py` | §3 → `cache/bg_lrcheck.txt`, `bg_eetrunc.txt`, `bg_k2check.txt` | < 5 min |
| `bg_qan.py`, `bg_covfit.py`, `bg_ffscan.py` | §4 → `cache/bg_qan.txt` (+ pkl), `bg_covfit.txt`, `bg_ffscan.txt` | 1 min each |
| `bg_qrun.sh` | the whole chain in order | ≈ 45 min |
