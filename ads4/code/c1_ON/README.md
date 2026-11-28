# c1_ON — sign of c₁ at order N⁰ in the 3D critical O(N) model (ledger)

Convention: log Z(S¹_β×S²_r) = 4πr² f T² − 8π c₁ + O(T⁻²); β = 1 throughout the code. Work by a research agent;
nothing committed.

## Formulation (derived here, [D] unless marked [C])

* log Z₁ = −½ Σ_{n,J} (2J+1) log Π_{n,J}, Π = G² bubble of the N = ∞ saddle (m² = m₀² + 1/(12r²) + …, conformal coupling).
* Split log Π = log Π^vac + δ, with Π^vac(ω,J) = (r/16)|Γ((J+1+iωr)/2)/Γ((J+2+iωr)/2)|² the exact R×S² vacuum
  (Δ = 1 two-point) kernel evaluated at Matsubara ω_n. Its full 3D determinant is exact:
  −½Σ(2J+1) log Π^vac = −S(β/r), S(x) = Σ_N σ₂(N) e^{−Nx}/N = ζ(3)/x² + (1/12) log x + ζ'(−1) + O(x²).
* Curvature expansion of the bubble (constant curvature K = 1/r², u = J(J+1)K, any rotation-invariant thermal kernel):
  **Π_n(u;K) = Π⁰_n(u) − (K/3) u ∂²_u Π⁰_n(u) + O(K²)**, Π⁰ the flat thermal bubble at mass m₀ = 2 log φ.
  (Reason: on S¹×M_K the propagator is G = (1 + K s²/12) G⁰ because the R/6 heat-kernel term cancels against
  ξR + δm² = K/4 + K/12; the area-element correction cancels the G² correction; what is left is the correction of the
  zonal spherical function, Φ = J₀(ks) + K[(s²/12)J₀ − (s/6k)J₁].)
* Euler–Maclaurin on S² (half-integer ℓ = J+½): Σ_J (2J+1) F(J(J+1)K) = (1/K)∫F du + F(0)/3 + O(K).

## Milestone 1 — bubble validated [C]

| check | result | file |
|---|---|---|
| flat thermal bubble (1D formula) vs direct 2D integral of closed-form Matsubara kernel | agree to 1e-15 | `flat.py`, `out` in run |
| analytic u-derivatives of Π⁰ vs numerical | agree to 30 digits | `check_flat_deriv.py` |
| large-p: 8pΠ⁰ − 1 → −128(h₁−h₂)(ω² − p²/3)/p⁵ (T_μν OPE term; no p⁻², no isotropic p⁻³) | ratio → 1 (1.0096 at p ≈ 70) | `check_flat_deriv.py` |
| exact S² 3j-sum bubble vs Π⁰ − (K/3)uΠ⁰'' | residual·r⁴/Π → const (e.g. n=0,J=0: 0.0138, 0.01308, 0.01292 at r = 4, 8, 16) | `check_s2_bubble.py`, `out_check_s2_bubble.txt` |
| S² saddle (m² − m₀²) r² | 0.08406, 0.08351, 0.08338 at r = 4, 8, 16 → 1/12 | same |
| vacuum: T=0, m=0 3j-sum bubble on R×S² vs Γ formula (normalisation r/16) | agree to 1e-13 | `check_vacuum.py`, `out_check_vacuum.txt` |
| curvature expansion of Π^vac at ω = 2π | (exact − 1/8p)/pred → 1 + O(K) | same |

## Milestone 2 — constant-term extraction (method B, semi-analytic) and preliminary sign [C]

Decomposition (β = 1): log Z₁ = −S(1/r) − ½ Σ_n S_n(r), S_n = Σ_J (2J+1) δ_{n,J}, δ = log(Π/Π^vac), summed n-outer.
* n ≠ 0: S_n = r²A_n + B_n + O(r⁻²), with A_n = ∫du δ⁰_n and B_n = δ⁰_n(0)/3 + ∫du D_n,
  D_n = −(u/3)[Π⁰''/Π⁰ − 3/(4(ω²+u)²)] (`c1_modes.py`, `out_modes.txt`, n ≤ 600, B_n ~ 0.006/n³).
  The per-sector expansion was checked against direct J-sums: the residual is ∝ r⁻² (`check_sector_EM.py`).
* n = 0: the reference is gapless, so the sector is done by direct J-sum plus Euler–Maclaurin boundary terms
  (`n0_sector.py`, `run_n0.py`, `out_n0.txt`). Result: S₀ − A₀r² − (1/6)log r → B₀ with a clean r⁻² approach:
  0.35260712, 0.35257637, 0.35256869, 0.35256678, 0.35256630, 0.35256618 at r = 12.5 … 400.
  **The log coefficient 1/6 is confirmed.** It cancels the reference's (1/12)log(r/β) exactly.
* Term 1 expansion checked (`check_term1.py`).

| quantity | value |
|---|---|
| −ζ'(−1) (reference constant) | +0.1654211437 |
| B₀ | 0.352566138 (±2e-9) |
| Σ_{n≥1} B_n | 0.0054749095 (±1e-12) |
| **b = constant of log Z₁** = −ζ'(−1) − ½[B₀ + 2ΣB_n] | **−0.0163368** |
| **c₁⁽⁰⁾ = −b/(8π)** | **+6.500e-4** |

**f₁ calibration and a regulator subtlety.**
* The n-outer sum of the conditionally convergent T_μν term (2ω² − k²)/p⁵ differs from a Lorentz-invariant
  regulator by a finite flat amount, 1/(3π²) per unit Σ∫ (cf. arXiv:2309.02347, eq. (2.43)).
* Corrected: f₁ = [−ζ(3) − ½(A₀ + 2ΣA_n)]/(4π) − c_T/(18π²) = −0.2294083 + 0.1654105 = **−0.0639978**. The
  literature value is −0.06399553 (2309.02347, eq. (1.11), their sign convention): agreement to 4e-5 relative.
* That this ordering difference only touches the r² term, and not the constant, was checked numerically on S¹×S²
  (`check_ordering.py`). With a Gaussian covariant cutoff, covariant − n-outer = 0.4244 r² + O(Λ⁻²), and the r⁰ part
  → 0 like Λ⁻². So c₁ is unaffected.

**Preliminary sign (go/no-go): c₁⁽⁰⁾ > 0.** The conjecture survives at O(N⁰). The result is small: it comes from a
cancellation between 0.1654 and 0.1818. Independent brute-force check (method A) pending.
