# Session 4: the k = 2 Taylor coefficient (the first five-point modular cumulant)

Status: 27 Nov 2026, research agent. Grades: [P] proof, [D] derivation, [D-sketch], [DV] numerically verified, [A] assumption.
Conventions: THEOREM_H.md. The observable is 𝒱(n) = 4π∂_Tη_Y/∂_Ts − 1 on R×H², with T = T₀/n. a_k are its Taylor coefficients at n = 1.

## 1. a₂ as a vacuum correlator

**Lemma M (left–right form of modular cumulants) [P, formal level for unbounded K].** Let K_R and K_L be the half-sided
(Rindler) modular Hamiltonians, so that the boost is B = K_R − K_L and B|Ω⟩ = 0. Then for every m and every A in the right-wedge
algebra (including boost-translates of A):

  ⟨Ω|K_R^m A|Ω⟩ = ⟨Ω|K_L^m A|Ω⟩ = ⟨Ω|A K_L^m|Ω⟩.

*Proof.*
1. ⟨Ω|B = 0, and [K_L, B] = 0, since K_L commutes with K_R and with itself.
2. By induction, ⟨Ω|K_R^m = ⟨Ω|(K_L + B)^m = ⟨Ω|K_L^m. Every word containing B can be reordered so that B stands leftmost.
3. K_L commutes with A, because the two are spacelike separated. ∎

**Consequence.** The temperature derivatives of H² quantities at T₀ are cumulants of *left-wedge* modular Hamiltonians with
*right-wedge* operators: correlators across the horizon of the entangled vacuum. In particular:
* the k = 2 coefficient involves ∂_T³η_Y. That is the vacuum 5-point cumulant ⟨O_Y O_Y K_L K_L K_L⟩_c, where O_Y is the H²
  shear strain, integrated along boost orbits for the Kubo limit;
* it also involves ∂_T³s ~ ⟨K⁴⟩_c (4-point);
* it also involves lower cumulants, from expanding the ratio and from the n ↔ T map.

K_L = 2π∫_{x⁺<0}|x⁺|T₊₊ on a null plane is a weighted light-ray integral. Under the Hofman–Maldacena map the null plane goes to
null infinity, so these are weighted energy-detector correlators [D-sketch]. The weight x⁺ means they are not the unweighted ℰ of
the collider literature. They are boost-charge light-ray operators.

## 2. Decomposition of δa₂ (from the bulk table, THEOREM_H §5) [DV + D]

At first order in the bulk couplings 16πG𝓛 = R + 6 + λ_J J + λ_{I²} I² + β₅ IJ + …:

  δa₂ = 15 λ_J + 408 λ_{I²} − 2880 β₅ + 0·(J², I³, …) + O(λ²).

The terms map to vacuum data as follows:

| term | vacuum data | status |
|---|---|---|
| 15 λ_J | t₄, a 3-point kinematic descendant (λ_J = −α₃/6, t₄ = −240α₃) | controlled by the paper |
| 408 λ_{I²} | quartic contact, ⟨TTTT⟩ | controlled under the paper's quartic-positivity input (G₄ ≲ 92.8 Λ⁻⁶) |
| O(λ²) | products of lower-point data, e.g. t₄², which are exchange-of-T pieces of the 5-point function | not computed; second-order hyperbolic engine needed |
| −2880 β₅ | **the genuinely five-point part**: the quintic I₂I₃ ≅ IJ coefficient, the object of Remarks 3, 4, 7, 7′ | uncontrolled |

**Exchange structure of β₅ [D].** At tree level a heavy scalar X coupled to g₂XI + g₃XJ gives β₅ = Σ_X g₂g₃/M_X². More generally,
the X-exchange contribution to the hyperbolic Taylor data is Σ_X Q_X(g₂(X), g₃(X); x), a fixed quadratic form per species. So the
k = 2 exchange piece is the **off-diagonal (2,3) entry of a positive Gram matrix** G_{ij} = Σ_X g_i g_j/M². Hence

  |β₅^{exch}| ≤ √(G₂₂ G₃₃),   G₂₂ = 2a ≲ G₄ ≲ Δ_g⁻⁶ (4-point),   G₃₃ = 2c = the six-point (I₃²) diagonal.

**Required size.** Theorem H needs only |δa₂| ≲ Δ_g⁻⁴ C². That means |β₅| ≲ Δ_g⁻⁴/2880, against the natural size Δ_g⁻⁸. Through the
Gram bound, it is enough that the **six-point diagonal satisfies c ≲ Δ_g⁻²**, against the natural Δ_g⁻¹⁰. This is a very weak
requirement, but it is still a 6-point input.

The planar exchange response, (432g₂² − 1944g₂g₃ − 3240g₃²)/M², is an indefinite quadratic form. So the exchange sector is not
sign-definite either, which is consistent with Proposition NS.

## 3. The three proposed mechanisms, and why each fails

**(a) Five-point light-ray superconvergence** [D-sketch]. Commutativity [ℰ(n₁), ℰ(n₂)] = 0 inside ⟨T K ℰ ℰ T⟩ gives *relations*
for 5-point data: β₅ equals a heavy-side sum of the interference products λ_{TTX}λ_{TTTX}. The heavy side is **not positive**.
* A positive functional that removes the spectator also removes this interference (Remark 7′, support lemma).
* Any bound therefore needs the diagonal Σλ_{TTTX}², i.e. G₃₃: six-point data.
* Failure: superconvergence turns β₅ into an interference sum. It does not bound it.

**(b) Modular / Tomita–Takesaki positivity** [P]. The modular operator of the wedge is Δ = e^{−2πB}, and B|Ω⟩ = 0. So
Δ^s O|Ω⟩ = (Δ^s O Δ^{−s})|Ω⟩ is the *complex-boosted single-operator state*, and energy-flux positivity in such states is the
collider positivity of ⟨O ℰ O⟩ with complexified boosts: **3-point data only**.
* The higher modular cumulants come from the half-sided K_R, which is not a symmetry. Its positivity content is the Lehmann /
  complete-monotonicity structure of the H² thermal state. That bounds each theory separately and gives no small deviations
  (THEOREM_H §6).
* Relative-entropy positivity S(ρ_n‖ρ₁) ≥ 0 gives one-sided convexity constraints, for example c_H ≥ 0.
* Failure: Δ sees only boosts. The genuine 5-point content enters only through K_R, where the only positivity is per-theory.

**(c) Coherent source / Rényi states** [D].
* By Lemma M, e^{−sK_R}|Ω⟩ = e^{−sK_L}|Ω⟩ is a normalizable purification of ρ_R^{1+s/π} (per unit H² volume). So ρ^{1−z} states exist.
* Energy-flux positivity in them, ⟨ℰ⟩_n ≥ 0, is one-sided and O(1).
* Energy conservation makes it two-sided, but still O(1), with no Δ_g: 0 ≤ ⟨ℰ(n̂)⟩_n ≤ energy.
* Expanding in s, the 5-point coefficient (the cubic response K₅ s³) appears inside a polynomial that is nonnegative on an interval.
  That bounds it only relative to the other O(1) coefficients.
* Failure: positivity of a sourced response bounds it at O(1). The gap smallness in the paper came from heavy-spectrum sum rules,
  not from positivity alone, and those sum rules exist only at 4 points.

## 4. The one handle that does see β₅: causality in the thermal state [D-sketch; conjectural]

Remark 3: type-N backgrounds (vacuum shocks, ANEC probes) are blind to I₂I₃ at quadratic order. **Type-D thermal backgrounds are
not.**
* On the brane, or on the hyperbolic branch with x > 0, the IJ term changes the quadratic graviton operator through background-
  Weyl factors (∝ x²–x³).
* Its four-derivative part grows with momentum relative to Einstein, by a factor of order β₅ x³ k².

CFT microcausality holds in *every* state, so the thermal commutator ⟨[T, T]⟩_β must vanish outside the light cone. The
boundary-to-boundary signal time must therefore be non-negative for every bounce depth and every polarization. This is the logic
of Brigante–Liu–Myers–Shenker–Yaida, now applied to all couplings at O(1) depth, with probe momenta up to the EFT cutoff k ≲ Λ ~ Δ_g.
* If the β₅-induced delay takes both signs across polarizations or depths, the constraints are two-sided. Naive power counting
  gives |β₅| ≲ Δ_g⁻²·(lower couplings, O(1)), with x³ = O(1) on the planar brane.
* On its own, Δ_g⁻² is *weaker* than the Δ_g⁻⁴/2880 that Theorem H needs at the full rate.
* But Theorem H works with any ε. If thermal causality gave |δa_k| ≲ Δ_g⁻² C^k uniformly in k, the result would be
  |4πη/s − 1| ≲ Δ_g⁻² log Δ_g: a power law.
* On the CFT side, this regime (probe momenta from T up to Δ_g at fixed ω/k < 1) is the thermal lightcone / bulk-cone regime. There
  the minimal-twist multi-stress coefficients [T^m] dominate, cf. Karlsson et al. 2206.05544. The causality statement would
  constrain those coefficients, in which β₅ appears at m = 3.
* What is not done:
  * the bulk computation of the β₅ delay and its sign structure (d = 3 has no TT channel at k ≠ 0, so the vector and scalar
    channels are needed);
  * the CFT derivation (thermal lightcone bootstrap with the gap);
  * uniformity in m.

## 5. Extension in k

The Gram structure extends to all orders: |G_{ij}| ≤ √(G_{ii}G_{jj}). So uniform geometric bounds on the *diagonal* heavy weights,
G_{jj} = Σ_X g_j(X)²/M_X² ≤ ε C^{2j} (heavy couplings to the j-th multi-stress composite), would control every exchange
off-diagonal. That would give the uniform-in-k input of Theorem H for the exchange sector, with C < 2/√3. Contact terms not produced
by tree exchange (stringy towers) need the same statement for their spectral representation.

**Status.** No diagonal bound beyond j = 2 is known. The bottleneck is now isolated as *inclusive heavy spectral weights in
multi-stress channels*, Σ_X |⟨X|T^j⟩|²/M², for j ≥ 3. The positive, natural CFT-side carrier for these is the diagonal
2j-point function of T in its (T^j)(T^j) channel. A dispersive bound on that (a 3 → 3 analogue of the paper's sum rules) is the
missing tool. Thermal causality (§4) is the one candidate found that sees β₅ with a gap-enhanced strength.
