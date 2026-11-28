# Scan: U(N_c), N_f with W = V₊ + V₋ (BBP 1703.08460) against c_S/C_T ≤ 16πζ(3)

## Files

Scripts:
- `scan_monopole.py` — the two-sided scan, plus the one-sided machinery.
- `onesided_scan.py` — the one-sided scan. It writes incrementally.
- `bbp_chain_general.py` — the Ising-chain checks.
- `u2_ct.py` — the matrix model.

Outputs:
- `scan_monopole_twosided_output.txt`
- `scan_monopole_onesided_output.txt`
- `bbp_chain_general_output.txt`
- `bbp_chain_u2_output.txt`
- `bbp_wz_check_output.txt`

## Conventions

- C_T is in real-scalar units, computed as (32/π²)∂²_b F through the Andreief determinant.
- Higgs branch: rank M ≤ N_c, so dim_H = 2N_cN_f − N_c². The free theory there has C_T = 4·dim_H and c_S = 3.5·dim_H.
- The conjecture is violated iff ratio := C_T(T_M)/(4·dim_H) < 7/8. In that case R ≥ 7/(8·ratio).
- T_M is the W = V₊ + V₋ fixed point, at R_Q = 1 − (N_c+1)/N_f.
- At N_f = N_c+2 I use BBP's IR identification:
  - N_c = 1: WZ model with W = det M, 9 chirals at R = 2/3.
  - N_c ≥ 2: N_f² free mesons.

  I verified the corresponding Z identity numerically for (1,3), (2,4), (3,5) and (5,7) to at least 1e-11.

## Calibrations (all pass)

| check | result |
|---|---|
| free chiral | 4 |
| XYZ | 8.7214 |
| N=4 SQED_n (1910.03600) | agrees to 1e-9 |
| Sp(1) N_f=3 vs U(1) N_f=4 (vector term) | equal |
| U(2) N_f=8 vs BBP dual U(4)+M | equal to 3e-9 in C_T, 1e-14 in F |
| topological-mixing machinery: U(1) N_f=2, W=V₊ vs its BBP WZ dual | F and C_T agree to 1e-10 |

## Two-sided results

Full table in `scan_monopole_twosided_output.txt`. Columns:
- R_Q* is the F-maximised R-charge of plain SQCD (T₀). At every point, all mesons and monopoles of T₀ have R ≥ 1/2, so no decoupling was needed.
- R[V]₀ is the monopole R-charge at T₀; criterion (i) needs R[V]₀ < 2.
- C_T(T_M) and the ratio are at the W fixed point.

| N_c | N_f | R_Q*(T₀) | R[V]₀ | (i) | C_T(T_M) | 7/8 threshold | ratio |
|---|---|---|---|---|---|---|---|
| 1 | 3 | 0.437 | 1.69 | yes | 26.16 | 17.5 | 1.308 |
| 1 | 5 | 0.461 | 2.69 | no | 32.93 | 31.5 | 0.915 |
| 1 | 6* | 0.467 | 3.20 | no | 34.00 | 38.5 | **0.773** |
| 2 | 4 | 0.385 | 1.46 | yes | 64 (free) | 42 | 1.333 |
| 2 | 5 | 0.410 | 1.95 | yes | 78.31 | 56 | 1.224 |
| 2 | 7 | 0.438 | 2.94 | no | 89.45 | 84 | 0.932 |
| 2 | 8 | 0.446 | 3.43 | no | 92.01 | 98 | **0.822** |
| 2 | 9* | 0.452 | 3.93 | no | 93.74 | 112 | **0.732** |
| 3 | 5, 6 | | 1.24, 1.72 | yes | 100, 138.65 | 73.5, 94.5 | 1.19, 1.28 |
| 3 | 10 | 0.433 | 3.68 | no | 173.17 | 178.5 | **0.849** |
| 3 | 11 | 0.439 | 4.17 | no | 176.26 | 199.5 | **0.773** |
| 4 | 6, 7, 8 | | 1.04, 1.50, 1.97 | yes | 144, 209.96, 235.46 | 112, 140, 168 | 1.13, 1.31, 1.23 |
| 4 | 12–14 | | 3.9–4.9 | no | 277.2–285.4 | 280–336 | **0.866–0.743** |
| 5 | 7, 8, 9 | | 0.84, 1.29, 1.75 | yes | 196, 291.0, 329.2 | 157.5, 192.5, 227.5 | 1.09, 1.32, 1.27 |
| 5 | 15–17 | | 4.7–5.6 | no | 410.1–419.3 | 437.5–507.5 | **0.820–0.723** |

\* marks N_f = 3(N_c+1). There BBP's existence argument is only a conformal-manifold proposal.

**Criteria (i) and (iii) never hold together.** At every point where the monopole deformation is relevant at T₀ (R[V]₀ < 2), the ratio is between 1.09 and 1.33, so there is no violation. Every violation sits at R[V]₀ ≥ 2.2. **Under the coordinator's strict criterion, no point passes all three.**

**Unitarity (ii) holds at every violating point.**
- Electric side: R_Q ≥ 0.58, so mesons have R ≥ 1.17. V± are in W with R = 2. Higher fluxes have R ≥ 4.
- Dual U(N_f−N_c−2) side: q has R = (N_c+1)/N_f; qq̃ vanishes by the M F-term; M has R = 2R_Q; V̂± have R = 2.
- N_f > 4(N_c+1)/3 at all of these points, so BBP say T_M is fully interacting there.

**For larger N_f the theories are inconsistent anyway.** For N_f ≥ 3N_c+4 the ratio keeps falling, but F(T_M) < F(free Higgs-branch theory); for example (1,7) has 4.359 < 4.506, and (2,10) has 12.10 < 12.48. That violates the F-theorem for the flow onto the moduli space, so T_M cannot exist there. This boundary coincides exactly with BBP's existence window N_f ≤ 3N_c+3. Inside the window, F(T_M) > F(Higgs-free) at every violating point.

## Existence through BBP's Ising-SCFT chain (§4.1 of 1703.08460)

Checked for N_f = 3N_c+2, N_c = 1…5, in `bbp_chain_general_output.txt`. The chain is

T₀ ⊕ N_f² Ising ⊕ N_f² σ → T_{ΦQQ̃+Φ³} → T_{V+ΦQQ̃} → T_M

- **Relevance at each step:**
  - ΦQQ̃ has R ≈ 1.53–1.59.
  - V at R_Q = 2/3 has R = 5/3.
  - σΦ has R ≈ 1.21–1.30.
- **F decreases strictly at every step.** For example, for (5,17) the sequence is 257.9 → 234.1 → 230.4 → 53.6.
- **BBP's alternative branch is excluded.** At T_{ΦQQ̃}, R_Φ* = 0.54–0.60 < 2/3, so Φ³ is relevant and option 1 is selected.
- **Smaller N_f.** These are reached from N_f = 3N_c+2 by complex masses. F(T_M) increases monotonically with N_f, as these flows require.

**Points that violate the bound and lie on established BBP-chain existence (N_f ≤ 3N_c+2):**

| point | ratio | lower bound on R |
|---|---|---|
| U(2) N_f=8 | 0.822 | 1.065 |
| U(3) N_f=10 | 0.849 | 1.031 |
| U(3) N_f=11 | 0.773 | 1.132 |
| U(4) N_f=12 | 0.866 | 1.010 |
| U(4) N_f=13 | 0.800 | 1.093 |
| U(4) N_f=14 | 0.743 | 1.177 |
| U(5) N_f=15 | 0.820 | 1.067 |
| U(5) N_f=16 | 0.769 | 1.138 |
| U(5) N_f=17 | 0.723 | **1.210** |

**Conformal-manifold proposal points (N_f = 3N_c+3, weaker existence argument):**

| point | ratio | lower bound on R |
|---|---|---|
| U(1) N_f=6 (abelian) | 0.773 | 1.13 |
| U(2) N_f=9 | 0.732 | 1.19 |
| U(3) N_f=12 | 0.709 | 1.23 |
| U(4) N_f=15 | 0.693 | 1.26 |
| U(5) N_f=18 | 0.682 | 1.28 |

At N_f = 3N_c+2 the violation grows with N_c: the ratio goes 0.915 → 0.822 → 0.773 → 0.743 → 0.723.

## Higgs branch (general version of the U(2), N_f=8 check)

- On the dual side, a rank-N_c vev of M leaves U(N') with N'+2 flavours and the monopole W.
- By BBP, that theory is WZ / free with (N_f−N_c)² mesons. I checked this numerically.
- Those mesons become massive with the corresponding block of M.
- What remains is N_f² − (N_f−N_c)² = 2N_cN_f − N_c² free chirals, with no leftover gauge group.
- This matches the electric count of complete Higgsing.
- Concrete case, U(1) N_f=6: 36 − 25 = 11 free chirals.

## One-sided W = V₊

R_Q and the topological mixing η are F-maximised along R[V₊] = 2. Results are in `scan_monopole_onesided_output.txt` for N_c = 1–4 and N_f = N_c+2 … 3N_c+4.

- The F-maximum is interior everywhere.
- Unitarity holds everywhere: mesons have R ≥ 0.56, and R[V₋] ≥ 0.61.
- **No violation.** The minimum ratio is 0.896, at (1,7). Ratios are 1.02–1.27 in the existence window N_f ≤ 3N_c+1, where existence comes from a real mass on T_M[N_f+1].
- As with the two-sided case, relevance at T₀ is the same condition and fails wherever the ratio approaches 7/8.

## Conclusion

- **Strict criterion (relevance at T₀):** no point passes all three conditions. Wherever the monopole W is relevant at plain SQCD, C_T(T_M) is 1.09–1.33 times the Higgs-branch free C_T, so there is no violation.
- **Existence via BBP's 3d UV completion (Ising-SCFT chain, N_f ≤ 3N_c+2):** nine (N_c, N_f) points violate c_S/C_T ≤ 16πζ(3), with R ≥ 1.01–1.21. Each of these points also has:
  - unitarity on both sides of the duality;
  - an F-theorem-consistent chain with F decreasing at every step;
  - F(T_M) > F(Higgs-free theory);
  - a 28-type Higgs-branch count, verified.

  The best established case is U(5) with N_f = 17, R ≥ 1.21. The simplest is U(2) with N_f = 8, or abelian U(1) with N_f = 6 on the conformal-manifold proposal.
- **What these rest on:**
  - BBP's chain flows terminate where indicated (standard SUSY RG reasoning, not a proof).
  - Localization C_T with no accidental symmetry at T_M.
  - The moduli-space thermal bound p ≥ p_Higgs-free (MODULI.md, CANDIDATE_U2.md §5).
