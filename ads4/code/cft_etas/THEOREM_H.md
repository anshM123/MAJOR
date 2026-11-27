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

## 5. General first-order response on the hyperbolic branch, and the Taylor data [DV]

**Engine.** `hyp_num.py` (sympy, fixed rational r_h) computes the membrane flux after the H² integral, in the det-preserving TT
gauge of the planar engine. It was run at r_h = 1, 5/4, 3/2, 2, 3, 5, 10 (outputs in `out_hyp_num_*.txt`). Results:
* There is no mass term at any r_h.
* η₀ ∝ r_h², so 4πη/s = 1 at every T.
* The fits below agree with all seven runs in **exact rational arithmetic**:

  e₁[θI] = 4t₀ − 2f′(r_h) t₁,    e₁[θJ] = x(−48 − 54x) t₀ + 6x f′(r_h) t₁,    f′(r_h) = 3r_h − 1/r_h = 4πT,

  where t₀ = θ(r_h) and t₁ = θ′(r_h).

**Wald entropy [D].** For L = R + θI, the Wald ratio is s₁/s₀ = 4θ C^{tr}{}_{tr}, and the horizon eigenvalue is C^{tr}{}_{tr} = x.
For θJ, use the type-D Weyl projection P(M²) = (A/2)M of the bivector square. Both are calibrated on the planar engine's 4θ and 6θ:

  σ₁[θI] = 4x t₀,   σ₁[θJ] = 6x² t₀.

**General response [DV + D].** At the horizon I = 12x², J = 12x³, r_h I′ = −72x², r_h J′ = −108x³ and f′/r_h = 2 + x. So, for
16πG𝓛 = R + 6 + F(I, J) at first order,

  δ(x) = 4(1−x) θ_I − 2(2+x) r_hθ_I′ + (−48x − 60x²) θ_J + 6x(2+x) r_hθ_J′,   with θ_y = F_y on the background,

  𝒱(x) = δ + (1−x) δ′.

`hyp_coeffs.py` → `out_hyp_coeffs.txt` gives these per-operator results:

| F | 𝒱(x)/λ | 𝒱 at x = 1 (planar) | Taylor coefficients a_k/λ at n = 1, k = 0, 1, 2, 3 |
|---|---|---|---|
| I | 0 (GB-blind) | 0 | 0, 0, 0, 0 |
| J | −48 − 120x + 60x² | −108 ✓ | −48, 120, 15, −165/4 |
| I² | 1344x − 96x² − 384x³ | 864 ✓ | 0, −1344, 408, 414 |
| IJ | −2880x² − 2016x³ + 2952x⁴ | −1944 ✓ | 0, 0, −2880, 4176 |
| J² | −14976x³ − 2448x⁴ + 10944x⁵ | −6480 ✓ | 0, 0, 0, 14976 |
| I³ | 89856x³ − 24192x⁴ − 34560x⁵ | 31104 ✓ | 0, 0, 0, −89856 |

All six planar values reproduce the paper's table. J² and I³ were **not** used in the fits, so they are independent checks.

**What the table shows.**
* **The counting.** An operator with m Weyl factors enters first at a_{m−3}. So a_k involves (k+3)-point vacuum data:
  * a₀ is ⟨TTT⟩ (t₄);
  * a₁ is ⟨TTTT⟩ (quartic contacts);
  * a₂ is the first 5-point coefficient, where the quintic IJ (I₂I₃) enters;
  * a₃ is 6-point.
* **Note.** x(n) = (1−n) + (3/8)(1−n)² − (1/32)(1−n)³ + …, so x ≈ 1 − n near T₀.
* **An exact CFT identity at T₀ [D].** ∂_T η at T₀ is the parity-even ⟨K; TT⟩. Since ∂_T s ∝ c_T and η(T₀) ∝ C_T,
  𝒱(T₀) is affine in t₄ in *any* CFT₃. It vanishes for Einstein, and its slope is fixed by the J row (λ = −α₃/6, t₄ = −240α₃):

  **𝒱(T₀) = 4π ∂_Tη/∂_Ts |_{T₀} − 1 = −t₄/30.**

  With the paper's bound this gives |𝒱(T₀)| ≤ τ_*(Δ_g)/30. This is an exact, gap-controlled k = 0 statement, analogous to the
  shear sum rule. It is testable on free fields (not done).

## 6. What the bound must look like, and why positivity cannot supply it

**(a) and (b) are one statement [P].** By Cauchy, |a_k| ≤ εC^k for all k holds if and only if |𝒱| ≲ ε on the disk
D(1, r) with r < 1/C. So hypothesis (a) is: *the deviation is uniformly small for complex temperatures near T₀*, a near-vacuum
statement. Theorem H then transports it to n = 0 using only the O(1) analyticity bound (b), at a cost of ε → ε^θ.

**Consistency with the collaborator's counterexample.** In F_N = ((1−q)/(1+q))^N, the point 0 lies on the boundary of
Re q ≥ 0. On any disk D(1, ρ) with ρ > 1, sup|F_N| grows like (ρ/(2−ρ))^N. That is exactly why Theorem H needs a *uniform* M on
an open domain containing n = 0.

**Required growth.**
* In the bulk, each operator contributes a polynomial in x(n). Its Taylor coefficients decay like (√3/2)^k, so C = √3/2 < 1 for
  any finite set of operators.
* Summed over operators, C ≤ 1 holds if the m-graviton couplings obey |λ_m| ≲ Δ_g^{−4} K^m with K below a fixed constant.
  Natural power counting, λ_m ~ Δ_g^{−(2m−4)}, satisfies this easily.
* With C ≤ 1 and the Einstein radius ρ = 2/√3, Theorem H gives |4πη/s − 1| ≲ Δ_g⁻⁴ log Δ_g. That is the rate in the target.

**What positivity gives [P/D].**
1. **Lehmann positivity.** Z(β) and Λ(β) = Z η_Y/β are Laplace transforms of positive measures on the H² spectrum: Λ ∝ Σ_E
   e^{−βE}|O_{EE}|². So both are completely monotone, analytic for Re β > 0, and bounded by |g(β)| ≤ g(Re β) (Bernstein). This gives
   *per-theory* geometric bounds on Taylor coefficients at β₀, with radius < β₀, i.e. |n − 1| < 1. Two problems remain:
   * n = 0 is on the boundary of that disk, not inside it.
   * The ratios in 𝒱 need zero-free regions (Fisher zeros), which positivity does not give.
2. **Multi-point energy-correlator positivity.** Commuting ANEC operators give positive k-point energy correlators, and energy
   conservation gives sum rules. Together these give *per-theory boundedness* |a_k| ≤ C^k, not smallness of deviations.
3. **Modular positivity.** Δ = e^{−2πK_full} ≥ 0 gives log-convexity in the modular parameter: one-sided inequalities, e.g.
   specific heat ≥ 0. K itself is not positive.

None of these produces *small* deviations. Positivity bounds each theory, not the difference between two theories.

## 7. The bottleneck

1. **Single open input: uniform-in-k smallness of the modular-cumulant deviations.**

   |δa_k| ≤ Δ_g⁻⁴ C^k for all k, with C < ρ (C ≤ 1 for the full Δ_g⁻⁴ log Δ_g rate).

   * Equivalently: the H² thermal free energy and Kubo correlator are Einstein-close to O(Δ_g⁻⁴) uniformly for complex T in a
     disk around T₀.
   * Status by order:
     * k = 0 is controlled exactly: 𝒱(T₀) = −t₄/30.
     * k = 1 is ⟨TTTT⟩ data, controlled under the paper's quartic-positivity input.
     * k ≥ 2 needs vacuum 5-point and higher data (I₂I₃ first), for which no bound exists.
   * The closest known analogue is the d = 2 sparseness theorem of Hartman (2013): large c plus a sparse light spectrum gives
     universal Rényi entropies. A d = 3 analogue for hyperbolic/Rényi quantities at all orders in (n − 1) is exactly what is
     needed.
2. **Secondary [A]: the O(1) analyticity domain (b).** n = 0 must be an interior point of a domain of uniform boundedness. It sits
   at the boundary of the axiomatic region Re β > 0. In Einstein gravity the domain is D(1, 2/√3); the singularities are spinodal
   (Fisher-zero) branch points at n = ±i/√3. On the CFT side this is a statement about Fisher zeros of the planar H² partition
   function.
