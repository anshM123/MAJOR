# Audit: U(2), N_f=8, W = V₊ + V₋ as a counterexample to c_S/C_T ≤ 16πζ(3)

Scripts:
- `u2_ct.py` builds the squashed-sphere matrix model.
- `candidate_u2.py` runs the audit. Its output is in `candidate_u2_output.txt` (about 10 minutes to run).

Convention: C_T is quoted in units of the free real scalar, C_T/C_s = (32/π²) ∂²_b F at b=1. Because the integrand is invariant under b → 1/b, the first b-derivative vanishes. So ∂²_b F = −⟨∂²_b log(integrand)⟩, and that is evaluated exactly.

## Verdict (revised after the SCAN_MONOPOLE.md follow-up)

**The first version of this audit concluded "existence FAIL". That was wrong; corrected in §3b.**

- **The number is right.** C_T = 92.01494393, F = 11.5448642825 > 9.7041, no operator violates unitarity at R_Q = 5/8, and the Higgs branch has dimension 28 with free chirals at its generic point.
- **The two simplest flows do not reach the fixed point (§3a):**
  - From plain SQCD, the monopole deformation is irrelevant: Δ(V±) = 3.43.
  - From the dual U(4)+M theory, V̂± are below the unitarity bound.
- **But BBP (1703.08460 §4.1) give a 3d UV completion through auxiliary Ising SCFTs, valid for N_f ≤ 3N_c+3. Our case N_f = 8 = 3N_c+2 is exactly their worked case.** I checked it numerically (§3b):
  - every deformation along the chain is relevant;
  - every intermediate fixed point respects unitarity;
  - F decreases strictly at each step (54.86 → 51.14 → 48.60 → 11.54);
  - BBP's alternative branch (Φ³ dropped) is excluded, because Φ³ is relevant at T_{ΦQQ̃}.
- **Consequence: this is a serious counterexample candidate, R ≥ 98/92.015 = 1.065.** It rests on standard assumptions: BBP's chain, localization with no accidental symmetry at T_M, and the thermal and moduli argument of §§4–5. None of my checks failed.

## 1. C_T at R_Q = 5/8 — computed, PASS (matches the claim)

**Method.** I used the Andreief determinant: Z_b = det_{j,k}∫dσ g_b(σ) e^{2π(bj + k/b − (N−1)Q/2)σ}. Its b-derivatives were taken analytically:
- The chiral pair term ∂²_b log[s_b(iy−σ)s_b(iy+σ)] at b=1 was written as a single w-integral and checked against finite differences to 1e-7.
- The vector term is ∂²_b log[sinh(πbx)sinh(πx/b)] = 2πx coth πx − 2π²x²/sinh²πx, with root x = σ₁ − σ₂.
- σ-integrals used the trapezoid rule on [−L, L], which converges exponentially here.

**Calibrations.**

| check | result | reference |
|---|---|---|
| free chiral | 4.0000000001 | 4 |
| XYZ | 8.72141 | 8.7214 |
| N=4 SQED_2 | 15.3029815148 | 15.3029815144 (1910.03600) |
| N=4 SQED_4 | 30.8910554574 | 30.8910554568 (1910.03600) |
| N=4 SU(2) with 3 flavours (tests the non-abelian vector term) | 30.8910554580 | must equal SQED_4 |

The last row relies on Sp(1) with N_f=3 being the same SCFT as U(1) with N_f=4. Both have so(6) ≅ su(4) symmetry, the Coulomb branches agree (D₃ = A₃), and the Higgs branches agree.

Side finding: the closed forms for SU(2)/SO(n) in 1910.03600 (its SO(5)/SO(6)/SO(8) table entries) reproduce my code only if the vector term uses root σ instead of 2σ. They give SO(6) = 39.51 ≠ SU(4) = 30.89. I believe those entries are wrong; with the corrected vector term the minimal SO(8) SCFT has C_T = 50.84 rather than 56.59. The MODULI.md conclusions are unchanged.

**Result.**

| quantity | value |
|---|---|
| U(2), N_f=8, R_Q=5/8 | **C_T = 92.0149439321** |
| grid change (L 8→10, h 0.01→0.005) | 0 to 1e-10 |
| independent 2D quadrature | identical |
| BBP dual U(4)+64 M at R_q=3/8 | **92.0149439346** |

The error is about 3e-9. The claim of 92.01494393 is confirmed, and so is the interval (92, 93).

**Symmetry check.** W = V₊ + V₋ carries U(1)_T charge ±1 and axial charge −N_f, so it breaks both. The diagonal baryonic U(1) is gauged. Only SU(8)×SU(8) and charge conjugation remain, so there is no abelian flavour symmetry to mix with and R_Q = 1 − 3/8 is fixed. PASS.

**Comparison without W.** Plain U(2) N_f=8 SQCD, after F-maximisation, has R_Q* = 0.44582, F* = 14.0710 and **C_T = 127.313**. The monopole term therefore lowers C_T by about 28%.

## 2. F-theorem test — PASS

- F(U(2), R_Q=5/8) = **11.5448642825**. The U(4)+M dual gives the same value; they differ by 9e-15, which confirms the BBP identity numerically.
- This exceeds F(28 free chirals) = 9.7040605.
- It is below F(plain SQCD) = 14.071, as required for any flow out of SQCD.

## 3. Accidental symmetries and existence of the fixed point

**Unitarity at the candidate point: PASS.**

| operator | R-charge |
|---|---|
| mesons QQ̃ | 1.25 |
| V± (flux (±1,0)) | 2 |
| flux (1,−1) | 4 |
| flux (1,1) | 6 |
| dual M | 1.25 |
| dual qq̃ | 0.75 (zero in the chiral ring) |
| dual V̂± | 2 |

There is no adjoint chiral, so there are no adjoint-dressed monopoles. U(2) has no baryons.

### 3a. The conservative flows fail (correct but not decisive)
- **Electric side.** At the plain SQCD fixed point, R[V±] = 8(1 − 0.44582) − 1 = **3.43 > 2**, so the deformation is irrelevant.
  - Every SU(8)²-invariant deformation that breaks U(1)_T involves monopoles with R ≥ 3.4.
  - Mesonic superpotentials such as (QQ̃)^k fix R_Q ≤ 1/2, which gives R[V] ≥ 3.
  - So no relevant flow from SQCD reaches R_Q = 5/8.
- **Dual side.** For the undeformed U(4) N_f=8 + M theory (W = Mqq̃):
  - Naive F-maximisation gives R_q* = 0.58283 and F = 26.58.
  - There R[V̂±] = 0.337 < 1/2, so V̂± are free and decoupled.
  - The decoupling-corrected F keeps increasing up to the convergence edge R_q → 5/8. Along the way the flux (1,−1,0,0) monopole also drops below 1/2, at R_q > 0.594.
  - The deforming operators V̂± are therefore free fields. W = V̂₊ + V̂₋ is then linear in free chirals, which signals SUSY breaking or a deformation that is not well defined, rather than a flow to an SCFT.
  - The relevance requirement R[V̂] < 2, i.e. R_q > 3/8, does hold, but it is moot.
- F-ordering is consistent in principle (26.6 > 11.54 and 14.07 > 11.54), but neither simple flow works.

### 3b. The BBP Ising-SCFT chain — PASS

Script `bbp_chain_u2.py`, output in `bbp_chain_u2_output.txt`. The chain is

T₀ ⊕ 64 Ising (Φᵢⱼ³) ⊕ 64 free σᵢⱼ → T_{ΦQQ̃+Φ³} ⊕ σ → T_{V+ΦQQ̃} ⊕ σ → T_M

**Step 1: couple Φ to the mesons.**
- ΦQQ̃ has R = 2/3 + 2(0.4458) = 1.558, so it is relevant.
- The IR has R_Q = R_Φ = 2/3. There V± has R = 8/3 − 1 = 5/3, and Φ has R = 2/3, so nothing is below 1/2.
- BBP's second option is T_{ΦQQ̃} with Φ³ dropped. F-maximising there gives R_Q* = 0.7163 and R_Φ = 0.567 < 2/3. So Φ³ is relevant at T_{ΦQQ̃}, and the flow ends at T_{ΦQQ̃+Φ³}. Option 2 is excluded.

**Step 2: add W = V₊ + V₋.**
- At R_Q = 2/3 this has R = 5/3, so it is relevant.
- The IR has R_Q = 5/8 and R_Φ = 3/4; Φ³ becomes irrelevant.

**Step 3: add the mass term σΦ.**
- σΦ has R = 5/4, so it is relevant.
- Integrating out Φ and σ leaves exactly T_M.

**F-theorem along the chain.**

| fixed point | F |
|---|---|
| T₀ ⊕ 64 Ising ⊕ 64 σ | 54.8623 |
| T_{ΦQQ̃+Φ³} ⊕ σ | 51.1443 |
| T_{V+ΦQQ̃} ⊕ σ | 48.6007 |
| T_M | 11.5449 |

F decreases strictly at every step.

**Remaining assumptions.** The flows end at the fixed points indicated (standard), and no accidental symmetry appears beyond those checked. BBP also give the dual version of the chain in their §4.2, which avoids the free-V̂ problem of the naive U(4)+M route.

## 4. Higgs branch — PASS, conditional on the fixed point existing

- **Electric side.** The Higgs branch is rank(M) ≤ 2, of dimension 2·2·8 − 4 = 28. At a generic point U(2) is completely Higgsed: 4 chirals are eaten and 28 free chirals remain. The monopole operators are heavy there.
- **Dual side.**
  - The F-term of M sets qq̃ = 0.
  - A vev of M of rank r gives mass to r flavours, leaving U(4) with 8−r flavours and W = V̂₊ + V̂₋.
  - That theory has a supersymmetric vacuum only if 8 − r ≥ 4 + 2, so r ≤ 2 and the dimension is 2(8+8−2) = 28.
  - At r = 2, the U(4) theory with 6 flavours and W is BBP-dual to U(0), i.e. 36 free mesons. I checked this numerically: log Z = −9.306556097384 on both sides. Those 36 mesons pair with the 6×6 block of M through W = M·(qq̃) and become massive.
  - What remains is the 28 tangent components of M: free, with no leftover gauge group or TQFT.

## 5. Thermal step — the argument holds

- **Setup.** Swap circles: the spatial circle has length L = β with antiperiodic fermions, so SUSY is broken and the vacuum energy density is ε₀ = −p.
- **Trial state.** On a large transverse torus, a broad zero-mode wavepacket centred at a Higgs-branch point v is a legitimate state. In finite volume there are no superselection sectors.
- **Energy of the trial state.** Its energy density is ε₂₈ + δ(v) + O(1/(A·Δv²)), where ε₂₈ is the energy density of 28 free chirals and A is the transverse area.
  - The SUSY-breaking Casimir potential on the moduli space has two sources:
    - massive modes with m ~ |v|², which contribute e^{−mL};
    - irrelevant sigma-model couplings, suppressed by powers of 1/(|v|²L).
  - So δ(v) → 0 as |v|²L → ∞.
- **Order of limits.** Take the transverse size to infinity first, which defines p and makes the kinetic term vanish. Then take v → ∞. This order is harmless.
- **Variational bound.** ε₀ ≤ ε_trial, hence p ≥ p₂₈.
  - If the Casimir potential pulls the state toward the origin, the bound still holds.
  - If it drives a runaway, p = p₂₈.
  - The only escape would be to restrict the conjecture to the metastable origin state.
- **The S² definition of 2610.33179 (eq. S3, the full spectral sum) gives no such escape.** The conformal mass costs π|v|²/T, which is negligible against the thermal gain of order (RT)² when 1 ≪ |v|²/T ≪ (RT)².

## What remains assumed

- The fixed point exists via BBP's Ising chain (§3b), under standard SUSY RG assumptions. All the checks I can run pass.
- Localization at R_Q = 5/8 computes the IR C_T, i.e. there is no hidden accidental symmetry. All gauge-invariant chiral operators I listed, on both sides, have R ≥ 1/2. BBP also state that T_M is completely interacting for N_f > 4(N_c+1)/3.
- BBP duality itself. Its partition-function identity is confirmed numerically to 1e-14 at this point.
- Given these, §§1, 4 and 5 give c_S ≥ 98 and C_T = 92.015, so **R ≥ 1.065**.
