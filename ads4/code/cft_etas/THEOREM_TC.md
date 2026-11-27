# Session 5: thermal causality ("thermal CEMZ") and the quintic β₅

Status: 27 Nov 2026, research agent. Grades: [P] proof, [D] derivation, [D-sketch], [DV] numerically verified, [A] assumption.
Script: `eikonal_check.py` → `out_eikonal_check.txt`.

## 1. Leading-eikonal propagation for 16πG𝓛 = R + 6 + F(I, J) in a type-D background [D; symbolic check]

**Setup.**
* At high momentum, k ≫ T (in AdS units, up to the cutoff Λ ~ Δ_g), a graviton is a local null plane wave h = ε e^{iS} with
  k = ∂S.
* At leading eikonal order its linearized Weyl tensor is type N. In the self-dual bivector block, δW = v vᵀ with v·v = 0 (and the
  conjugate in the anti-self-dual block).
* The background brane (and the hyperbolic branch) is type D: W = q·diag(−2, 1, 1). On the brane q(r) = −1/(2r³), so the exterior
  is q ∈ [−1/2, 0). On the hyperbolic branch q ∈ [−x/2, 0).

**The second variation.**

  δ²F = F_I δ²I + F_J δ²J + F_II (δI)² + 2F_IJ δI δJ + F_JJ (δJ)²,   δI = 2 tr(WδW) = −6q v₁²,   δJ = 3 tr(W²δW) = 9q² v₁².

The O(k⁴) parts of δ²I and δ²J are 2tr(δW²) and 3tr(WδW²). Both vanish, because δW² = 0 (type N). So the
**momentum-enhanced** (k⁴ relative to Einstein's k²) term is

  δ²F|_{k⁴} = 36 q² v₁⁴ · K(q) · k⁴,   K(q) = F_II − 3q F_IJ + (9/4) q² F_JJ.

Consequences:
* For real waves, v₁⁴ becomes a square of Re(·), which is ≥ 0. The polarization only changes the size, never the sign.
* The terms F_I δ²I and F_J δ²J change the cone only at O(k²), with no enhancement. This is how C³ (t₄) produces O(1),
  collider-type thermal-causality constraints and no Λ⁻² constraint. It agrees with the known equivalence between near-boundary
  thermal causality and the conformal-collider bounds.

**Calibration [D].**
* Pure I² (F_II = 2a) gives K = 2a. So causality ⇔ a ≥ 0: the paper's quartic positivity α₄ ≥ 0.
* Heavy-scalar exchange (g₂XI + g₃XJ) gives K = (2g₂ − 3qg₃)²/(4M²) ≥ 0: a **perfect square**, as unitarity requires.
* These two facts fix the sign convention: K ≥ 0 ⇔ subluminal, i.e. time delay.

## 2. Answer to question (2): β₅ is one-signed

For F = β₅ IJ alone, K(q) = −3qβ₅ = 3|q| β₅ in the exterior.
* This has **one sign for every polarization, every channel** (the eikonal analysis covers the tensor, vector and scalar channels
  at once), **and every bounce depth**, since q never changes sign outside the horizon.
* Causality therefore constrains β₅ from **one side only**: β₅ ≥ −(slack).

With the slack from the finite cutoff, the condition is K(q) ≳ −κ/(Λ² q² z_*³). For a = c = 0 this gives
β₅ ≥ −O(Λ⁻²) = −O(Δ_g⁻²).

**With the quartic and sextic diagonals** (F = aI² + bIJ + cJ², b = β₅): K = 2a + 3b s + (9/2) c s², with s = |q| ∈ (0, 1/2].
Causality ⇔ K ≥ 0 on the segment, which means:
* a ≥ 0 and c ≥ 0;
* b ≥ −2√(ac) when the minimizer s* = −b/(3c) lies inside the segment. That is exactly Remark 4's Gram/Cauchy–Schwarz cone,
  b² ≤ 4ac;
* b ≥ −(4/3)a − (3/4)c otherwise.

All of these hold up to the O(Λ⁻²) slack.

**So thermal causality at leading eikonal order is the positivity of the Hessian of F along the curve (δI, δJ) ∝ (−6q, 9q²).**
* It reproduces the heavy-exchange positivity cone, and nothing beyond it.
* Positive b is unconstrained. That is the sign that *lowers* η/s, since the planar shift from IJ is −1944b. So it is the
  KSS-violating direction.
* The negative side needs an upper bound on the six-point diagonal c.

## 3. Comparison with Theorem H (question 3)

* Theorem H needs |β₅| ≲ Δ_g⁻⁴/2880 from **both** sides.
* Thermal causality gives only β₅ ≥ −min(2√(ac), (4/3)a + (3/4)c) − O(Δ_g⁻²).
* **No rate follows.** The mechanism fails for the same structural reason as Remarks 4 and 7′: unitarity-type positivity bounds the
  interference term only through the diagonals, and the 6-point diagonal is uncontrolled.
* There is one by-product. For b < 0 (the direction that *raises* η/s), b ≥ −(4/3)a − (3/4)c with a ≲ Δ_g⁻⁶. This gives an upper
  bound on the IJ-induced *increase* of η/s, in terms of a and c. KSS needs the opposite direction.

## 4. What remains, and other operators (question 4)

* **Gradient operator Y = ∇I·∇J** (planar shift −34992λ). The leading-eikonal part of δ²Y is 2(k·k)δIδJ, which vanishes for null k.
  The remaining terms carry background gradients and are down by 1/k. **Thermal causality is blind to Y at the enhanced order**,
  just as type-N shocks are blind to I₂I₃. Every operator whose fluctuation derivatives are contracted with each other is suppressed
  the same way.
* **General algebraic F(I, J).** Only the projected Hessian K(q) is constrained. That is one functional per point q, against the
  three-jet functional that η/s needs. Directions in coupling space with K ≡ 0 on the exterior segment exist and are invisible to
  causality. Example: b = (2a + (9/2)cs²)/(−3s) cannot vanish for all s unless a = b = c = 0. But infinitely many higher monomials
  (I²J, J³, …) give enough freedom to keep K(q) ≡ 0 while shifting η/s. So **uniformity in the number of Weyl factors fails**
  [D-sketch].
* **What a CFT derivation would still need.**
  1. Thermal microcausality in the CFT. This is axiomatic: commutators vanish at spacelike separation in any state.
  2. The thermal lightcone bootstrap: show that the k ≫ T, fixed ω/k < 1 regime of ⟨TT⟩_β is governed by minimal-twist
     multi-stress data, and turn "no superluminal signal" into constraints on those coefficients, with the gap setting Λ.
  3. In any case, the outcome would only reproduce the positivity cone above. It cannot supply the two-sided bound Theorem H needs.

## 5. Verdict

The "thermal CEMZ" route is a precise negative result [D]:
* at leading eikonal order, thermal causality is equivalent to heavy-exchange positivity, K(q) ≥ 0, of the Hessian of F along a
  curve;
* it bounds β₅ only from one side, and only through the six-point diagonal;
* the KSS-violating sign is unconstrained;
* gradient operators are invisible to it.

The open input is unchanged: inclusive heavy weights in the multi-stress channels j ≥ 3 (the six-point diagonal first), or some
non-positivity mechanism that controls interference terms.
