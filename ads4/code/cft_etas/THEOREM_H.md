# The hyperbolic route: planar η/s as the n → 0 endpoint of a Rényi-type continuation from the vacuum (session 3)

Status: 27 Nov 2026, research agent. Grades: [P] proof, [D] derivation, [D-sketch], [DV] numerically verified, [A] assumption.
Units: L_AdS = R_{H²} = 1.

## 1. Setup

**The state.** Put the CFT₃ on R × H² with R = 1, at temperature T = T₀/n, where T₀ = 1/2π. At n = 1 the state is conformal to the
Minkowski vacuum restricted to a Rindler wedge (Casini–Huerta–Myers). In the bulk this is the massless hyperbolic black hole,
i.e. AdS–Rindler. The planar thermal state on R² is the limit T R → ∞, that is n → 0.

**The shear mode on H² [D].** In two dimensions a traceless, divergence-free symmetric tensor is a holomorphic quadratic
differential. Every TT tensor on H² therefore satisfies the same Laplace identity, which in Kodama–Ishibashi notation reads
k_T² + 2K = 0 with K = −1. So the TT sector carries a single "momentum", and the bulk tensor master equation reduces to the
**massless scalar** equation at zero momentum. This is the exact analogue of the planar k = 0 shear mode.

The explicit L² representative used in the code is H = sech²X (dX² − cosh²X dY²) in the metric dX² + cosh²X dY².
`hyp_test_einstein.py` checks, among other things, that α₀ = 0 for Einstein gravity.

**Kubo/membrane η on H².** Define η_H(T) = −lim_{ω→0} Im G_R^{TT}(ω)/ω for this mode.
* In Einstein gravity the membrane value is η = r_h²/16πG and s = r_h²/4G, so η/s = 1/4π at every T [D, standard].
* As T R → ∞ the tensor channel has no hydrodynamic pole. The fixed eigenvalue is ~1/R², which is zero momentum in units of T,
  so η_H(T)/s(T) → the planar η/s.

**The right observable: the GB-blind ratio [D].** Use

  𝒱(n) := 4π (∂_T η_H)/(∂_T s) − 1,   evaluated at T = T₀/n.

Why this ratio:
* **Gauss–Bonnet drops out.** The 4D Gauss–Bonnet term is topological. It shifts s on H² by a T-independent constant (its Wald
  entropy ∝ horizon R⁽²⁾ × area = const per unit H² area) and changes no correlator. So 𝒱 is exactly GB-blind.
* **Why not η/s itself.** η/s at n = 1 measures the global ratio c_T/F, because s(T₀) = F/(2πR²). GB shifts that ratio, and
  nothing local controls it.
* **Same planar limit.** In the planar limit η, s ∝ T², so 𝒱 → 4πη/s − 1 (planar).
* **Einstein is identically zero.** 𝒱 ≡ 0 for Einstein gravity, at all n.
* **Bulk formula.** At first order in couplings, with η = η₀(1 + e₁(r_h)) and s = s₀(1 + σ₁(r_h)), and η₀ = s₀/4π = r_h²/16πG on
  any static background:

  𝒱 = (e₁ − σ₁) + (r_h/2)(e₁ − σ₁)′ = δ(x) + (1 − x) δ′(x),

  where δ := e₁ − σ₁ is the hyperbolic η/s deviation and x = 1 − 1/r_h² is the horizon Weyl scale. The O(λ) corrections to the
  background and to T(r_h) cancel in the ratio of r_h-derivatives.

## 2. Taylor coefficients at n = 1 are vacuum correlators [D-sketch]

**β-derivatives are cumulants.** At n = 1 the H² Hamiltonian is the modular Hamiltonian K of the Rindler wedge. Since
∂_β Tr(e^{−βH} A B(t))/Z = −⟨H; A B(t)⟩_c, with H inserted on a slice away from A and B (H is conserved), we get:
* the k-th T-derivative of the Kubo correlator at T₀ is a connected vacuum (k+2)-point function of T, with k copies of the
  modular Hamiltonian K = 2π∫_{x⁺>0} x⁺ T_{++};
* the s-derivatives are ⟨K^{m}⟩_c.

**Regularity.**
* **UV.** In the S¹×H² frame different K's sit at different Euclidean times, at proper distance ≥ R·Δτ everywhere. There is no
  entangling-surface divergence for *densities* on H².
* **IR.** Connected correlators decay like e^{−Δ·dist} on H², and the volume grows like e^{dist}. The integrals converge because
  Δ_T = 3 > 1.

**Light-ray reading.** K is a weighted null-plane integral of T_{++}. So the Taylor data are **multi-point energy (light-ray)
correlators** of the vacuum. The first two:
* 𝒱(1) involves ∂_T η ~ ⟨K T T⟩ (t₄, α_T) and ∂_T s ~ ⟨KK⟩ (c_T);
* the x^k coefficient of the bulk response involves (k+3)-point vacuum data.

This matches the bulk count: an operator built from m Weyl tensors first enters at order x^{m−2}.

## 3. Theorem H (Cauchy / two-constants reconstruction) [P]

**Hypotheses.**
* **(b)** 𝒱 is analytic on a domain Ω ⊃ [0, 1] in the n-plane, with |𝒱| ≤ M on Ω.
* **(a)** The Taylor coefficients at n = 1, a_k = 𝒱^{(k)}(1)/k!, satisfy |a_k| ≤ ε C^k.
  * Here ε is the vacuum deviation scale. For example, ε ~ |t₄| ≲ Δ_g⁻⁴ for the leading term.
  * Recall that 𝒱 ≡ 0 for Einstein.

**(i) Disk version.** If Ω ⊇ D(1, ρ) with ρ > 1:

  |4πη/s − 1|_planar = |𝒱(0)| ≤ min_K [ ε Σ_{k<K} C^k + M ρ^{1−K}/(ρ − 1) ].

So |𝒱(0)| ≲ ε^{θ} M^{1−θ}, with θ = log ρ / log(ρ·max(C, 1)), up to logarithmic factors when C ≤ 1. Checked numerically
(`hyp_cauchy.py` → `out_hyp_cauchy.txt`; exponent of Δ_g with ε = Δ_g⁻⁴):

| (C, ρ) | (1, 1.1547) | (1.1, 1.1547) | (1.5, 1.1547) | (2, 2) |
|---|---|---|---|---|
| exponent | 3.93 | 2.41 | 1.05 | 2.00 |
| prediction 4θ | 4 | 2.406 | 1.047 | 2 |

**(ii) General domain.** If Cr₀ < 1, then (a) gives |𝒱| ≤ ε/(1 − Cr₀) on D(1, r₀). The two-constants theorem gives
|𝒱(0)| ≤ (ε/(1 − Cr₀))^{θ} M^{1−θ}, where θ = ω(0; D(1, r₀), Ω) is a harmonic measure. For Ω = D(1, ρ),
θ = log ρ / log(ρ/r₀).

**This is a power of Δ_g for every admissible Ω.** The target n = 0 lies in the **interior** of the analyticity domain, at finite
harmonic distance from the data. Compare:
* Theorem C (Euclidean data): the target is a boundary point, approached through a discrete set. Rate (log)^{−2γ}.
* Theorem R2 (early-time data): the target is at time → ∞. Rate exp(−c(log)^{0.4}).

## 4. The Einstein analytic structure [DV]

In Einstein gravity:
* r_h(n) = (1 + √(1 + 3n²))/(3n);
* x(n) = 1 − 9n²/(1 + √(1 + 3n²))², which is analytic at n = 0, with x(0) = 1;
* the only singularities are the branch points n = ±i/√3. These are the spinodal points dT/dr_h = 0, of Lee–Yang type.

Consequences:
* The Taylor radius about n = 1 is ρ = 2/√3 = 1.1547 > 1. So **n = 0 is inside the disk of convergence**, and every planar
  quantity that is analytic in x is a *convergent* series in vacuum Taylor data, with ratio √3/2 = 0.866 per order.
* Check: the Taylor series of x(n) about n = 1, summed at n = 0, has errors 4.7e-2, 5.7e-3, 8.3e-5, 1.7e-6 at K = 10, 20, 40, 60.
  This follows ρ^{−K}.
* At first order in couplings, 𝒱 = δ + (1 − x)δ′ with δ polynomial in x (see §5). So 𝒱 inherits exactly this domain.

## 5. Test on Einstein + λF (in progress)

`hyp_engine.py` and `hyp_test_einstein.py` compute the hyperbolic responses e₁, σ₁ for θ(r)I and θ(r)J and hence δ(x). The
results go in the README ledger.

## 6. Bottleneck (preliminary)

* **(a) is the bottleneck.** It asks for all-k bounds on vacuum multi-point energy (modular-Hamiltonian) correlators of T, with
  growth C^k and C < ρ for a useful exponent. Only k = 0 (c_T) and k = 1 (t₄, α_T, the paper's bounds) are controlled. k = 2
  (⟨TTTT⟩) is controlled up to the quartic positivity input. k ≥ 3 involves ≥ 5-point vacuum data (I₂I₃), which is uncontrolled.
* **(b) is O(1).** It is natural: no phase transition for n ∈ (0, 1], since hyperbolic black holes are stable for T ≥ T₀, and
  Rényi transitions occur at n > 1. Plus analyticity of the high-temperature curvature expansion at n = 0. But it requires
  continuing slightly past Re n = 0. The axiomatic domain of the thermal partition function is Re β > 0, i.e. Re n > 0, so n = 0
  is an axiomatic *boundary* point. At planar order a single saddle continues analytically (Einstein: up to |n| < 1/√3), but
  that is an [A].
