# Independent test of the claimed identity 𝒱(T₀) = −t₄/30 (CFT₃ on R × H²)

Status: 27 Nov 2026, verification agent. Units R = 1, T₀ = 1/2π. Grades: [P] proof, [D] derivation, [DV] numerically verified.

## 0. Verdict

| | 𝒱(T₀) computed here | claim −t₄/30 |
|---|---|---|
| free Dirac (or Majorana) fermion, t₄ = −4 | **+2/15** (exact, 40 digits) | +2/15 ✓ |
| Einstein gravity, t₄ = 0 | 0 | 0 ✓ |
| free conformal real scalar, t₄ = +4 (honest thermal ensemble on H²) | **+6/5** (exact, 40 digits) | −2/15 ✗ (off by 4/3) |

* **The ⟨TTT⟩-level formula is confirmed.** 𝒱 is exactly affine in t₄, not merely rational. Through its two ⟨TTT⟩-faithful points
  (Einstein and the fermion) the line is exactly −t₄/30. This agrees with the independent holographic slope from the Weyl³ row
  (a₀ = −48λ, t₄ = 1440λ).
* **"In ANY parity-even CFT₃" is false as stated.** The free conformal scalar is a counterexample: its true thermal 𝒱(T₀) is
  6/5, not −2/15.
  * The same scalar-only deviation shows up in the third cumulant ⟨HHH⟩ (Sec. 4). That is the known conformal-scalar Rényi
    anomaly: total-derivative (improvement) terms in the modular Hamiltonian; Lee–Lewkowycz–Perlmutter–Safdi,
    arXiv:1407.7816, and Hung–Myers–Smolkin, arXiv:1407.6429.
  * So the identity needs a hypothesis: no edge term in K = ∫x⁺T₊₊ inside ⟨K O_Y O_Y⟩.
* **Grade.** [D + DV] in the conditional form of Sec. 6. The unconditional "any CFT₃" version is refuted [DV].

## 1. Setup and checks of the inputs

* **Kubo normalisation.** η_Y = ρ_OO(ω)/(4 N_Y ω) as ω → 0, with ρ the commutator spectral function. In flat space this reduces to
  ρ_xy/2ω.
* **The shear tensor.** With ds² = 4|dz|²/(1−|z|²)², Y = Re(dz²) has
  √γ Y^{ij} ∂_iA ∂_jB = ½(1−|z|²)²(∂A∂B + ∂̄A∂̄B), and N_Y = π/6 ✓.
* **The scalar improvement term drops out of O_Y.**
  * Y is traceless, so the g_{ij} terms vanish.
  * G_{ij} = 0 on the spatial slices of R × H².
  * ∫Y^{ij}∇_i∇_jφ² = 0 by transversality. The boundary term is O(e^{−2ρ}).
  * So O_Y = ∫√γ Y^{ij}∂_iφ∂_jφ exactly.
* **The fermion in the conformal frame.** With χ = Ω^{1/2}ψ, T_zz = ½[χ₂*∂χ₁ − (∂χ₂*)χ₁]. All Ω-derivative terms cancel.
* **The entropy side depends only on C_T.** ∂_T s = c_H/T with c_H = β²⟨H;H⟩_c/V. That is a two-point (C_T-only) quantity;
  ⟨KKK⟩ enters only ∂_T c_H, not ∂_T s. So the denominator of 𝒱 is t₄-independent, and 𝒱 is **affine** in t₄ (not rational) for
  any ⟨TTT⟩-determined numerator.
* **Inputs reproduced** for both free fields, to 10⁻¹⁴ by direct quadrature and to 40 digits in closed form:
  * c_H(T₀) = π³C_T/12;
  * η_Y(T₀) = π²C_T/48 (Remark 14's normalisation);
  * c_H(T₀) = 4π η_Y(T₀).
  * The constants used are C_T = 3/32π² for the real scalar and Majorana fermion, and 3/16π² for the Dirac fermion.

## 2. Method: exact free-field η_Y(T) on H² (`spectral_kernels.py`)

* **Only particle–hole pairs contribute.** As ω → 0, only particle–hole terms survive (pair creation is ∝ ω³). So

  η_Y(T) = (π/2N_Y) ∫dE (−∂_E n_T) S(E,E)/norm,

  where S(E,E′) = Σ|⟨k|O_Y|l⟩|² is summed over the single-particle energy shells E and E′. S is T-independent: all the
  temperature dependence sits in n_T.
* **S(E,E′) by Schur's lemma.**
  * O_Y pairs the fields against the holomorphic quadratic differential dz². The L² holomorphic quadratic differentials form an
    irreducible (discrete-series) PSU(1,1) representation, with reproducing kernel 3/π(1−z w̄)⁻⁴.
  * Hence Σ|L_f|² = λ‖f‖², and the trace-per-area identity gives **λ = ∫d²w (1−|w|²)² 𝒦(0,w)**: a single radial integral.
  * Here 𝒦 is built from the H² spectral projectors:
    * **scalar:** μ(p)P_{−1/2+ip}(u), with μ(p) = p tanh(πp)/2π;
    * **component χ₁, χ₂ of the Dirac field:** F(σ)/(1 − w̄z), with F = (p coth πp/π)·cosh²(σ/2)·₂F₁(1+ip, 1−ip; 1; −sinh²(σ/2)).
      This kernel was derived here and checked against its ODE to 10⁻²⁶.
  * The positive-energy projector is ½(1 + h/E)P_E. Only its diagonal blocks enter T_zz.
* **Results.** The numerics fit the closed forms below to all digits. `out_spectral_kernels.txt` checks them at five values of p.

  * scalar: J(p) = p³ tanh(πp)(1+4p²)/(48π³)
  * Dirac: λ_F(p) = p coth(πp)(1+4p²)/(6π²)

  which gives the kinetic-theory-like exact formula [D + DV]

  **η_Y(T) = (1/24) Σ_species ∫₀^∞ dE ρ(E) (−∂n/∂E) (1 + 4E²R²).**

  Here ρ is the H² density of states per unit area:
  * scalar: ρ = E tanh(πE)/2π;
  * Dirac: ρ = 2·E coth(πE)/2π.
* **Sanity checks.**
  * At T₀ this gives exactly η_Y = π²C_T/48 for both fields.
  * As TR → ∞, η_Y ≃ ζ(3)T³R/2π for the scalar, so 4πη/s ∝ TR → ∞. This is the free-field infinite viscosity, with the mean free
    path cut off by the curvature.

## 3. Free-field values (`closed_form.py` → `out_closed_form.txt`, 40 digits)

Since c_H(T₀) = 4πη_Y(T₀), 𝒱(T₀) = d log η_Y/d log T |_{T₀} − 1.

| | d log η/d log T | 𝒱(T₀) | 𝒲 = d log c_H/d log T |
|---|---|---|---|
| real scalar | 11/5 | **6/5** | 34/15 |
| Dirac / Majorana fermion | 17/15 | **2/15** | 26/15 |
| Einstein | 1 | 0 | 7/4 |

* The direct quadrature from the hypergeometric kernels, which does not assume the closed forms, gives 11/5 and 17/15:
  * to 10⁻⁸ with N = 30 Gauss nodes (the saved output);
  * to 10⁻¹⁴ with N = 45 (an exploratory run).
* The closed-form kernels themselves agree with the hypergeometric integrals to 10⁻²¹.
* The Majorana value equals the Dirac value: both η and s halve.

## 4. ⟨TTT⟩ decomposition

* **Counting.** In d = 3, parity-even ⟨TTT⟩ = n_B·(scalar structure) + n_F·(Majorana structure). C_T ∝ n_B + n_F, and
  t₄ = 4(n_B − n_F)/(n_B + n_F). Here n_F counts **Majorana** fermions (C_T(Dirac) = 2 C_T(scalar)), which checks the formula
  quoted in the task.
* **What is affine.** Both 4π∂_Tη and ∂_Ts are additive over the structures, and ∂_Ts|_{T₀} ∝ C_T. So any ⟨TTT⟩-determined 𝒱 is
  𝒱 = ½(1 + t₄/4)v_B + ½(1 − t₄/4)v_F, where v_B and v_F are the structure values.
* **Fixing v_B and v_F.**
  * Einstein (t₄ = 0, 𝒱 = 0) forces v_B = −v_F.
  * The fermion gives v_F = 2/15.
  * Hence **𝒱_{⟨TTT⟩}(T₀) = −t₄/30** exactly, with boson-structure value v_B = −2/15.
* **The free scalar is off this line.** Its thermal value is 6/5. If it were ⟨TTT⟩-faithful, then together with the fermion it
  would force 𝒱(t₄ = 0) = 2/3 ≠ 0. That contradicts 𝒱 ≡ 0 for Einstein, which is exact (η = s/4π at every T).
* **The same pattern in an observable that does not involve η.** For 𝒲 (the ⟨HHH⟩ / Rényi S″_q datum), the line through
  Einstein and the fermion has slope 1/240 and predicts 53/30 for the boson structure. The thermal scalar gives 34/15, offset 1/2.
  This reproduces the known statement of LLPS 2014: the ⟨TTT⟩ expansion of Rényi entropies agrees with fermions and Einstein
  duals but not with the conformal scalar.
* **Not done.** I did not integrate the Osborn–Petkou boson structure through the CHM map directly. The boson-structure value
  −2/15 rests on the fermion, Einstein and Weyl³-slope agreement.

## 5. Subtleties

1. **Edge term (the scalar exception).**
   * T_tt^{imp} = ℋ_can − ξD²φ² (ξ = 1/8). Hence ∫T_tt^{imp} = H_can − ξ∮_∞ ∂_nφ².
   * H_can is the Fock Hamiltonian. It generates the unique KMS state, and at T₀ that state is exactly the vacuum (checked).
   * φ² has Δ = 1. Its connected thermal correlators decay like e^{−ρ}, so the boundary flux is O(1) and does not vanish. In mode
     language this comes from the soft Bose modes, n ≈ T/E at the gapless H² threshold.
   * So for the scalar, ∂_T at T₀ is **not** the insertion of K = ∫x⁺T₊₊^{imp}; an edge operator ∮φ² on the entangling circle is
     missing.
   * Fermions have no Δ ≤ 1 scalar in T and f(0) = ½ is regular, so there is no edge term. Holographic CFTs agree as well.
   * I did not compute the edge term directly. Its value, 4/3 in 𝒱, is inferred as the difference.
2. **Contact terms.** These are real polynomials in ω and do not enter Im G^R. The improvement term in O_Y vanishes identically
   (Sec. 1).
3. **Normalisation of η_Y.** The 2N_Y factor is right: it reproduces π²C_T/48 for two independent free fields.
4. **Order of limits.** Free-field ρ_T(ω)/ω is jointly smooth in (T, ω) near (T₀, 0), so the ω → 0 limit and ∂_T commute.
5. **Zero modes and volume.** There are no zero modes; both fields have E = p ≥ 0. Everything is per unit H² area. The O_Y matrix
   elements are L²-smeared, so no δ(ω) transport peak arises, unlike flat space.

## 6. Suggested statement for the paper (conditional form)

> **An exact CFT₃ relation between hyperbolic shear transport at the CHM temperature and ⟨TTT⟩.** Let a parity-even CFT₃ live on
> R × H²_R, with Y = Re(dz²) and η_Y the Kubo coefficient of O_Y = ∫√γ Y^{ij}T_{ij}. Suppose that the temperature derivative
> at T₀ = 1/2πR is the insertion of K = ∫x⁺T₊₊, with no edge operator on the entangling surface. Then
> 4π ∂_Tη_Y/∂_Ts|_{T₀} − 1 = −t₄/30. This is exact and affine in t₄, because ∂_Ts|_{T₀} = π⁴C_T/6R depends only on C_T. The
> free fermion (+2/15), Einstein gravity (0) and the Weyl³ coupling (slope −1/30) each check it independently. The free
> conformal scalar violates the hypothesis: its thermal value is +6/5. This is the edge-term anomaly familiar from its Rényi
> entropies.

## Files

| file | contents |
|---|---|
| `spectral_kernels.py` → `out_spectral_kernels.txt` | first-principles kernels, closed-form checks, direct T₀ quadrature (~4 min) |
| `closed_form.py` → `out_closed_form.txt` | 40-digit values of η, c_H, 𝒱, 𝒲; the affine-in-t₄ tests (seconds) |
