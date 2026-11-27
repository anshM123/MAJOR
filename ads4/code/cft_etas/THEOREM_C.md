# Theorem C: Euclidean (Matsubara) data and the Kubo limit. Exact identity, stability bound, sharpness

Status: 27 Nov 2026, research agent. Grades follow the paper: [P] proof, [D] derivation, [DV] numerically verified, [A] assumption.
Scripts: `einstein_ref.py`, `einstein_strip.py`, `pick_bound.py`, `sharp_opt.py`. Their outputs are in `out_*.txt`.

## 0. Setting and standing conditions (K)

We work with a CFT₃ in its KMS state on R² at temperature T and zero chemical potentials. Take G_R(ω) to be the parity-even part g
of the zero-momentum retarded correlator of T_xy (the notation of `sec_transport.tex`). Write δG = G − G_vac, and let ρ = Im G,
η(ω) = ρ(ω)/ω (even, total) and δη = η − η_vac. Thermodynamics gives ε = 2p and s_ent = 3p/T.

**(K)** [A; these are the conditions of the shear sum rule, Sec. B of the paper]:

* G_R is analytic in Im ω > 0.
* δG(ω) → J, a constant, as |ω| → ∞ in the closed upper half plane.
* δη ∈ L¹(R).
* δη is continuous at ω = 0. At c_T = ∞ this excludes a Drude delta function, which only free theories have.

Under (K),

  δG(ω) − J = (1/π) ∫ δρ(ω')/(ω' − ω) dω'.

J is the scheme-dependent contact constant. It never enters below: only differences of δG appear.

## 1. The identity [P given (K)]

For y > 0 the imaginary part of the dispersion integral vanishes, because δρ is odd. This gives

  [δG(0) − δG(iy)]/y = (1/π) ∫ δη(ω) · y/(ω² + y²) dω  =  (1/π) ∫_0^∞ Φ(t) e^{−yt} dt,

where Φ(t) := ∫ δη(ω) cos(ωt) dω is the vacuum-subtracted Kubo relaxation function. The second equality holds because the
Fourier transform of y/(π(ω²+y²)) is e^{−y|t|}.

Define the **Matsubara viscosity function**

  E(s) := 4π [δG(0) − δG(2πiTs)] / (2πT s · s_ent),   Re s > 0.

E(s) has three properties:

1. E(s) is a Poisson average of 4πδη/s_ent at width y = 2πTs.
2. E(0+) = 4πη/s, by Poisson-kernel convergence at a point where δη is continuous, and η_vac(0) = 0.
3. At s = n ∈ N, E(n) = 4π[G_E(0) − G_E(ω_n)]/(ω_n s_ent). This is a Euclidean S¹×R² correlator at Matsubara frequencies.
   * It is contact-term free.
   * s_ent comes from the free energy, so it is Euclidean too.
   * G_E(0) = G_R(0) requires no ω = 0 delta function, which is part of (K).

So

  **4πη/s − 1 = F(0+),   F := E − E_Einstein,**

and the Euclidean CFT data are the values F(n), n ≥ 1. The Kubo limit is a Laplace transform evaluated at s = 0, while the
Euclidean data are the same transform at s = 1, 2, 3, ….

*Relation to the Kubo–Mori formulation.* The fixed-ν quantity ⟨a, ν(Q+ν²)^{-1}a⟩ of the relative-form proposition is the
Poisson smear of the (regulated, total) η at width ν. After the vacuum part is subtracted, it equals E(ν/2πT)·s/4π. At
ν = ω_n it is therefore Euclidean.

**Check [DV]** (`out_einstein_ref.txt`). Einstein AdS₄ brane, units r₊ = 1 (2πT = 3/2, p = 1, η = 1, s = 4π).

* The Euclidean side, computed by solving the bulk equation at ω = iy, matches the real-frequency Poisson integral of δη to
  1e-10 for n = 1..8.
* E_E(n) = 0.4127522, 0.2381977, 0.1637545, 0.1239991, 0.0995750, … for n = 1, 2, 3, 4, 5.
* E_E(0+) = 1.
* δG(0) − δG(i∞) → 3/4 = ε S_E.

So the "smeared KSS ratio at resolution 2πT" for Einstein gravity is 4πη_{(1)}/s = 0.41275.

## 2. Theorem C (stability) [P]

**Hypotheses.**

* **(S_{γ,M})**: F extends analytically to H_γ = {Re s > −γ}, with sup_{H_γ} |F| ≤ M.
* **(E_{N,ε})**: |F(n)| ≤ ε for n = 1, …, N.

**Conclusion.**

  |4πη/s − 1| = |F(0)| ≤ ε Λ_N + (M + ε Λ_N) β_N,

  β_N := ∏_{n=1}^N n/(n+2γ) = Γ(N+1)Γ(1+2γ)/Γ(N+1+2γ),     Λ_N := Σ_{k=1}^N ∏_{j≠k, j≤N} (k+j+2γ)/|k−j|.

**Proof.** The Blaschke factor b_n(s) = (s−n)/(s+n+2γ) has |b_n| = 1 on Re s = −γ: the point s is equidistant from n and from
its mirror image −n−2γ. Inside H_γ, |b_n| < 1, and b_n → 1 at ∞. Define

* B = ∏_{n≤N} b_n and B_k = B/b_k;
* the interpolant L(s) = Σ_k F(k) B_k(s)/B_k(k).

The steps are:

1. F − L vanishes at s = 1..N, so (F − L)/B is analytic in H_γ: its singularities at the zeros of B are removable.
2. (F − L)/B is bounded in H_γ. Away from the N zeros, |B| is bounded below; near them the function is analytic.
3. Its boundary values satisfy limsup ≤ M + Σ_k ε/|B_k(k)| = M + εΛ_N, because |B| → 1 at the boundary and |B_k| ≤ 1 in H_γ.
4. The maximum principle for bounded analytic functions on a half plane (Phragmén–Lindelöf) gives
   |F − L| ≤ (M + εΛ_N)|B| on H_γ.
5. At s = 0, |B(0)| = β_N and |L(0)| ≤ εΛ_N. ∎

**Asymptotics [P for the β_N part; DV for Λ_N].**

* β_N ≤ Γ(1+2γ) N^{−2γ}.
* (1/N) log Λ_N → log(3+2√2) = 1.7627 (numerically 1.747 at N = 120, γ = 1).
* For γ = 0 there is a closed form, 1/|B_k(k)| = C(2k−1,k)·C(N+k,2k), which gives growth (3+2√2)^N via Σ_k C(N+k,2k)4^k ~ T_N(3).

Choosing N ≈ log(M/ε)/1.76 gives

  |4πη/s − 1| ≤ C_γ · M · (log(M/ε))^{−2γ}.

Numerically, upper·L^{2γ}/M ≈ 6.8–9 at γ = 1 for ε = 1e-10 … 1e-160, with L = log(M/ε) (`out_pick_bound.txt`).

**Physical reading.**

* **(S)**: G_R has no singularity in Im ω > −2πTγ, and the deviation of the normalized relaxation response from Einstein's is
  O(1)-bounded there.
  * A real-time sufficient condition: |Φ − Φ_E|(t) ≤ (s_ent/4)·2πT·m·e^{−2πTγ' t} with γ' > γ. This gives M ≤ m/(γ' − γ).
  * (S) does **not** assume Einstein. It holds whenever there is a uniform non-hydrodynamic relaxation gap and no near-axis
    resonance with a large residue.
  * It is exactly the hypothesis that the family of Remark 9 violates: there sup_{H_{γ/2}} |F_M| ~ 2^M.
* **Einstein values** (`out_einstein_strip.txt`):
  * The lowest k = 0 shear QNM is ω/2πT = 1.23295 − 1.77590 i, so γ_E = 1.776. The next is 2.1075 − 3.2776 i.
  * M_E(γ) = sup_{Re s = −γ} |E_E|:

    | γ | 0.25 | 0.5 | 0.75 | 0.9 | 1.1 | 1.25 | 1.5 |
    |---|---|---|---|---|---|---|---|
    | M_E(γ) | 1.35 | 1.89 | 2.74 | 3.48 | 4.81 | 6.39 | 14.5 |

    Hence M ≤ M_theory + M_E(γ).
* **(E)**: the Euclidean thermal ⟨TT⟩ at the first N Matsubara frequencies is within ε of Einstein, in the normalization of E. This
  is the unproved CFT input.

**Corollary C1.** If (E) holds with ε ≤ KΔ_g^{−4} for N up to (4/1.76)·log Δ_g, then

  |4πη/s − 1| ≤ C(γ, M, K) (log Δ_g)^{−2γ}.

This is an explicit rate.

Sample values (M = 3, from `out_pick_bound.txt`):

| γ | ε = 1e-12 | ε = 1e-24 |
|---|---|---|
| 0.5 | 0.209 | 0.102 |
| 1.0 | 0.0318 | 0.0075 |
| 1.5 | 0.0077 | 0.00087 |

## 3. Proposition C′ (sharpness) [P]

Take

  F_sh(s) = M · B_N(s) · ((a−γ)/(a+s))^k,   a > γ, k ≥ 2.

Its properties:

* **Bounded in the strip.** |F_sh| ≤ M on H_γ, since every factor has modulus ≤ 1 there (maximum principle).
* **Real and symmetric.** F_sh is real on the real axis, so G(−ω*) = G(ω)*.
* **OPE invisible.** F_sh = O(|s|^{−k}) at ∞, so the deviation δG_dev = −2πT s (s_ent/4π) F_sh decays at |ω| → ∞. It changes no
  OPE coefficient and no sum rule.
* **Small at every Matsubara point.** F_sh(n) = 0 for n ≤ N. For n > N, |F_sh(n)| ≤ M((a−γ)/(a+N+1))^k.
* **Size at s = 0.** |F_sh(0)| = M β_N (1−γ/a)^k.

With a = k = N = ⌈log₂(M/ε)⌉, every |F_sh(n)| ≤ ε and |F_sh(0)| ≥ M Γ(1+2γ) e^{−γ}(1+o(1)) (log 2/log(M/ε))^{2γ}. The upper and lower
bounds therefore have the same rate. The constants differ by about (1.76/log 2)^{2γ}e^{γ}.

* An optimized choice of (N, a, k) is computed in `sharp_opt.py`. It gives |F(0)|·L^{2γ} ≈ 1.29 at γ = 0.5, ε = 1e-12.
* **Positivity.** The deviation changes the normalized total viscosity 4πη(ω)/s by Re F_sh(−iω/2πT), and |Re F_sh| ≤ M on the
  real axis. The Einstein total viscosity is 4πη_E(ω)/s = 1 + O(ω²) and stays above its value at 0 (≥ 1, checked on the sampled
  grid in `out_sharpness_positivity.txt`). So M < 1 keeps ρ > 0.
* **Caveat.** As in Remarks 5 and 9 of the paper, this is a response function, not a constructed CFT correlator.

**Conclusion.** Euclidean closeness at every Matsubara frequency, together with any O(1) analytic control in a strip, gives
DC control only at a logarithmic rate, and that rate is optimal. A power-law rate (Δ_g^{−4}) needs a *small* real-time norm, as
in the relative-form proposition (ε, δ form bounds) or the disk-norm closure. Euclidean CFT data at precision ε ~ Δ_g^{−4} cannot
supply it. In terms of hypotheses: the form-comparison hypothesis controls F(s) for every real s ≥ 0, whereas Euclidean CFT
controls F only on N. Theorem C and Proposition C′ quantify exactly the price of that difference.

## 4. What remains open (honest list)

1. **(E) itself.** No CFT argument is known that bounds |E(n) − E_E(n)| by a gap power. Three facts are relevant:
   * The large-n tail is OPE-controlled: E(n) − E_E(n) = t₄/(24n) + O(n^{−4}). The O(n^{−4}) coefficient is set by the Δ = 6
     [TT] thermal one-point data, which is uncontrolled.
   * The first few n integrate the plane correlator beyond the radius of convergence of the thermal OPE (|x|² + τ² < β²).
   * In the bulk EFT, the quintic coupling I₂I₃, which is invisible in ⟨TTTT⟩, shifts 4πη/s by O(β₅). Theorem C then forces its
     Euclidean deviation to satisfy ε ≥ M·exp(−c|β₅|^{−1/2γ}) > 0. So (E) cannot follow from vacuum data up to four points plus
     the gap [D given the EFT deformation is a consistent CFT].
2. **(S) itself.** No CFT argument is known. The spectrum of the Liouvillian is not controlled by vacuum dimensions.
