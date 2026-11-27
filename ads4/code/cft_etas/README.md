# cft_etas: toward a CFT derivation of η/s. Running ledger

Research agent working alongside the human collaborator. Nothing outside this directory is edited, and nothing is committed.
Grades: [P] proof, [D] derivation, [D-sketch], [DV] numerically verified, [A] assumption.

## 2026-11-27, session 1

### Milestone 1: Route C. Euclidean data, the Kubo limit, and a sharp log-stability theorem (see THEOREM_C.md)

* **Identity [P given the sum-rule conditions (K)].** Define the Matsubara viscosity function
  E(s) := 4π[δG_R(0) − δG_R(2πiTs)]/(2πTs·s_ent).
  * E(s) is the Poisson smear of 4πδη/s at width 2πTs.
  * It is also the Laplace transform of the vacuum-subtracted Kubo relaxation function at 2πTs.
  * E(0+) = 4πη/s.
  * E(n) is a Euclidean S¹×R² Matsubara correlator, contact-term free; s_ent = 3p/T is Euclidean too.
  * Hence 4πη/s − 1 = F(0+) with F = E − E_Einstein, and the Euclidean CFT data are exactly F(1), F(2), ….
  * Link to the collaborator's Kubo–Mori proposition: their fixed-ν quantity, vacuum-subtracted, equals E(ν/2πT)·s/4π. At
    ν = ω_n it is Euclidean.
* **[DV] Einstein check.** The identity holds to 1e-10 for n = 1..8. E_E(n) = 0.41275, 0.23820, 0.16375, 0.12400, 0.09957. The lowest
  zero-momentum shear QNM is ω/2πT = 1.23295 − 1.77590i, so the relaxation gap is γ_E = 1.776.
* **Theorem C [P].** Hypotheses:
  * (S_{γ,M}): F analytic on Re s > −γ with |F| ≤ M there. This is a uniform relaxation gap plus an O(1) bound, and it does not
    assume Einstein.
  * (E_{N,ε}): |F(n)| ≤ ε for n ≤ N.

  Conclusion: |4πη/s − 1| ≤ εΛ_N + (M + εΛ_N)β_N, with β_N = ∏ n/(n+2γ) and Λ_N ~ (3+2√2)^N. Optimized over N this is
  |4πη/s − 1| ≤ C_γ M (log(M/ε))^{−2γ}. With ε ~ Δ_g⁻⁴ the rate is (log Δ_g)^{−2γ}, an explicit rate.
* **Proposition C′ (sharpness) [P].** There is an explicit response family F_sh = M·B_N(s)·((a−γ)/(a+s))^k with these properties:
  * it satisfies (S);
  * it is ≤ ε at every Matsubara point;
  * it leaves every OPE coefficient and sum rule unchanged;
  * it keeps ρ > 0 for M < 1;
  * it has |F(0)| ≈ c_γ M (log 1/ε)^{−2γ}.

  The numerically optimized constants are c ≈ 1.29 at γ = 0.5 and ≈ 3.36 at γ = 1. The upper bound has C ≈ 7 at γ = 1. Upper and
  lower bounds therefore agree within a factor of about 2 to 3.
* **Consequence (decisive for the program).** Euclidean CFT data, even at power-law precision and at all Matsubara frequencies,
  plus any O(1) real-time analytic control, yield at best a logarithmic rate for the DC ratio. A Δ_g⁻⁴ rate needs a *small*
  real-time norm: the collaborator's ε, δ form bounds, or the disk norm. This sharpens Remark 9: a bounded strip continuation is
  enough for a rate, and the rate can only be logarithmic.
* **Corollary [D].** Hypothesis (E) cannot follow from vacuum data up to four points plus the gap. The I₂I₃ deformation is invisible
  to ⟨TTTT⟩ but shifts η/s by O(β₅), so by Theorem C it must shift some Euclidean F(n) by at least about exp(−c β₅^{−1/2γ}). This
  assumes the bulk EFT deformation is a consistent CFT.
* Files:
  * `einstein_ref.py` → `out_einstein_ref.txt`: the identity check.
  * `einstein_strip.py` → `out_einstein_strip.txt`: the QNM, γ_E and sup norms M_E(γ).
  * `pick_bound.py` → `out_pick_bound.txt`: the Theorem C bound and the simple sharp family.
  * `sharp_opt.py` → `out_sharp_opt.txt`: the optimized sharp family.
  * `sharpness_positivity.py` → `out_sharpness_positivity.txt`: the Einstein total η(ω) ≥ η(0) on the sampled grid, used for
    positivity.

### Milestone 2: the bulk mechanism. Euclidean data are blind at the tip of the cigar; η is tip-local ([DV], `horizon_bump.py` → `out_horizon_bump.txt`)

* **Setup.** Use the Remark-8 shear class with a = b = 1 + A·χ(u), where χ is supported within u ≥ 1 − w of the horizon and is
  constant near it.
  * At the horizon a = b, which is the KMS condition: the Euclidean mode must behave as ρ^{|n|} at the tip.
  * Entropy is held fixed.
  * Exact result: 4πη/s − 1 = A.
* **Numerics** (A = 0.1, w = 0.4 → 0.025):
  * F(0+) = 0.1000 for every w.
  * F(1) ∝ w, F(2) ∝ w², F(3) ∝ w³. For example, at w = 0.025: F(1) = 2.4e-3, F(2) = 2.8e-5, F(3) = 2.5e-7.
  * The strip value F(−0.5) grows roughly as w^{−0.6}.
* **Reading.** Matsubara mode n vanishes as (1−u)^{n/2} at the tip. A deformation within w of the horizon therefore shifts E(n) by
  O(A wⁿ) and the strip norm M by O(A w^{−γ}), while it shifts η/s by A. Consequences:
  * This is the bulk content of Theorem C. Near-tip structure at depth w corresponds to a late-time tail of the relaxation function
    at t ~ log(1/w)/2πT, the Rindler redshift.
  * Hypothesis (S) is, in bulk language, "no structure below an O(1) depth near the horizon."
  * Static (n = 0) data are also only weakly sensitive: they scale as w² through A and w through the gradient coefficient.
* **Why power counting would still win.** In a bulk EFT with cutoff Λ ~ Δ_g, structure is smooth down to w_min ~ Δ_g⁻². For
  bump-type deformations, F(1) ~ A·w then gives A ≲ ε Δ_g², a power law. But the extremal (Theorem C′) families show that this
  needs *structural* input, such as a sign or the bump shape, beyond any analytic-norm hypothesis.
  * For example, if the deviation of the relaxation function is sign-definite and supported at t ≤ t_L, then
    |F(0)| ≤ e^{2πT t_L}|F(1)| [P, one line].
  * Without the sign, compact support at t_L only gives ε·exp(c·e^{πT t_L}) [D-sketch, Chebyshev on [e^{−t_L}, 1]].

### Milestone 3: rigor upgrades

* **Explicit constants [P].** β_N ≤ Γ(1+2γ)N^{−2γ}, from log-convexity of Γ. Λ_N ≤ (N+1)^{2γ}[0.427(3+2√2)^N + 1/2], from the
  Morgan–Voyce identity. Checked in `lambda_bound.py` → `out_lambda_bound.txt`.
* **Sharp family made fully OPE-invisible.** An extra factor exp(−λ√(s+γ)) makes the sharp family decay faster than any power in the
  closed strip, so every asymptotic coefficient is unchanged. See THEOREM_C.md §3.
* **Trivial but recorded [P].** The controlled data (c_T, t₄, α_T, the sum rules, positivity of ρ) do not fix even the
  finite-resolution value E(1). A positive bump in η near ω ~ T, compensated at high frequency and decaying faster than any power,
  changes E(1) at O(1) and changes none of them. So a finite-resolution KSS theorem *is* hypothesis (E) at n = 1, and it needs
  thermal input.

## Status summary (end of session 1)

**Proved.** The identity (§1 of THEOREM_C), Theorem C and Proposition C′ [P]; the bulk mechanism [DV].

**What this means for the program.**
* Any route through Euclidean thermal data is capped at the rate (log Δ_g)^{−2γ}. Even (E) itself needs vacuum data at five or
  more points (I₂I₃).
* A Δ_g⁻⁴ rate needs a small *real-time* (Liouvillian, or tip-local) norm, or structural sign/support information.

**Most promising next step.** Aim any all-n control (for example the collaborator's coherent-source phase bound) at a real-time
object. Candidates are the relaxation-function deviation on t ≲ log(Δ_g)/2πT, or the relative form bound on the low-frequency
(near-horizon) sector of L. Euclidean Matsubara data are not a good target.

## 2026-11-27, session 2: the real-time route (see THEOREM_R.md)

### Milestone 4: early-time Lorentzian control (Lemma R3) and Theorem R2

* **Lemma R3 [D-sketch].** The zero-momentum commutator at t < β/2 is fixed by the *convergent* thermal OPE. Two facts make this
  work: Lorentzian convergence for |t| + |x| < β, and causal support |x| ≤ t.
  * Heavy states contribute ≤ C(2t/β)^{Δ_g}.
  * Analytic (double-trace) terms give no commutator. So the non-factorized thermal data b₁ can drop out. This is proven for
    scalar structures and open for the spinning structures of T.
  * Holographic evidence [DV] (`early_time_check.py`): a horizon-layer deformation changes F(s) only by e^{−2 r_*(depth)·s}. The
    rates are 1.32 and 1.98, against tortoise predictions of 1.22 and 1.91, consistent once the width of the transition region is
    included. So no OPE coefficient changes, and the early commutator is untouched.
* **Theorem R2 [P].** Hypotheses: k(t) analytic and bounded by M in the sector |arg t| < α, with decay e^{−γt}, and |k| ≤ ε on
  (0, t₂].
  * Conclusion: |k(t)| ≤ M^{1−θ}ε^θ with the exact θ(t) = (2/π) arcsin((t₂/t)^{π/2α}).
  * Rate: exp(−c (log 1/ε)^{1/(1+p)}). For Einstein-like data (α → π/3, from the QNM angles 55.2°, 57.3°, 58.1°;
    `qnm_sector.py`) this is exp(−2.3 (log 1/ε)^{0.4}).
  * This is better than any power of log, but it is not a power of Δ_g.

### Milestone 5: the sign route is dead (Proposition NS)

* A one-signed δΦ turns the log into Δ_g^{−2.56} (math_transfer; my Theorem R1). But:
  * Einstein's own Φ_E oscillates at late times, with zero spacing 1.698–1.701 against π/Re ω₁ = 1.699;
  * a nonzero deviation is one-signed only if the first QNM in which the two theories differ is purely damped
    [P given the QNM expansion];
  * even monotone horizon-coupling deformations give strictly log-concave F(s), hence a non-one-signed φ [DV].
* **Conclusion.** No perturbative-around-Einstein principle can supply the sign.

### Chaos, pole-skipping, modular flow (THEOREM_R.md §4)

* The shear-channel pole-skipping points s = −n are horizon-*amplified*: F(−γ) ~ A w^{−γ}, while F(n) ~ A wⁿ and F(0) ~ A.
* Chaos is the energy-channel pole-skipping point at s = +1. It is CFT-accessible through the OTOC and tip-local at O(1), but it
  probes a different horizon coupling.
* Modular inclusions (Leutheusser–Liu, Borchers) fix the near-horizon *kinematics*, not the coupling normalization.
* **Concrete next test.** At linear order, compute δv_B for the θ_I I and θ_J J deformations on the brane, and check whether the
  functional is aligned with δ(4πη/s) = −6θ_I′ − 108θ_J + 18θ_J′. If aligned, chaos fixes η at linear order. If not, there is a
  precise obstruction.

### Milestone 6: map of the CFT data that are sensitive to the tip (THEOREM_R.md §6) [D-sketch]

* Only three kinds of CFT data are O(1)-sensitive to the tip: thermodynamics (s), chaos (λ_L, v_B) and hyperbolic T₀ data.
  * For planar R², the first two are uncontrolled thermal data.
  * The third is vacuum data, but it includes the global ratio c_T/F, which Gauss–Bonnet shifts.
* Planar η/s is the endpoint (x → 1) of an analytic continuation in temperature on H². The Taylor coefficients of that
  continuation are vacuum n-point data, so this is the CFT-native form of the all-n requirement.

## Status (end of session 2)

**Proved.**
* Theorem R2 [P]: from early-time real-time data, rate exp(−c(log 1/ε)^{0.4}).
* Proposition NS [P given the QNM expansion]: the sign route fails generically.
* Theorem R1 [P]: the sign gives a power law, but its hypothesis is false generically.

**Derived, with numerical support.**
* Lemma R3 [D-sketch]: the early-time commutator is OPE-controlled, b₁ may drop out, and heavy states are ≤ (2t/β)^{Δ_g}.
* Holographic evidence [DV]: deep (horizon) data never enter the OPE coefficients or the early commutator.

**Dead routes.**
* Every positivity structure: frequency positivity, Φ ≥ 0, one-signed δΦ.
* Modular kinematics.
* Euclidean data, which give at best a logarithmic rate.

**Most promising next step.** Bound a *small real-time norm* directly. The natural CFT carrier is the hyperbolic-temperature
continuation from T₀ (vacuum) to T = ∞ (planar). It needs:
* all-n vacuum gap bounds, for example from the collaborator's coherent-source phase bound;
* a uniform analyticity domain in the horizon Weyl scale x ∈ [0, 1].

With those, a Cauchy estimate gives a power law. This is the CFT image of Proposition R, without the bulk effective theory.

## 2026-11-27, session 3: the hyperbolic route (see THEOREM_H.md)

### Milestone 7: setup, observable and Theorem H

* **Setup.** The CFT lives on R × H² at T = T₀/n, with T₀ = 1/2π. At n = 1 this is the vacuum (Rindler); n → 0 is the planar limit.
* **Shear mode [D; DV].** TT tensors on H² all satisfy k_T² + 2K = 0, so the shear mode is a massless scalar at "zero momentum".
  * The engine confirms this: the reduced Einstein Lagrangian has no mass term (α₀ = 0), α₁ = −(π/4) r² f, and η₀ ∝ r_h², so
    4πη/s = 1 at every T in Einstein gravity.
  * Files: `hyp_engine.py`, `hyp_test_einstein.py` → `out_hyp_test_einstein.txt`.
* **Observable.** 𝒱(n) = 4π∂_Tη/∂_Ts − 1.
  * It is GB-blind: topological GB shifts s by a T-independent constant.
  * It tends to planar 4πη/s − 1 as n → 0.
  * It is ≡ 0 in Einstein gravity.
  * In bulk variables, 𝒱 = δ(x) + (1 − x)δ′(x), where x = 1 − 1/r_h² is the horizon Weyl scale.
* **Taylor data [D-sketch].** The coefficients of 𝒱 at n = 1 are connected vacuum correlators of T with insertions of the modular
  Hamiltonian K. Equivalently, they are multi-point light-ray/energy correlators.
  * They converge in the UV (no entangling surface in the H² frame) and in the IR (Δ_T = 3 > 1 on H²).
  * Order k involves up to (k+3)-point vacuum data.
* **Theorem H [P].** Hypotheses: (a) |a_k| ≤ εC^k; (b) 𝒱 analytic and bounded by M on a domain Ω ⊃ [0, 1].
  * Conclusion: |planar 4πη/s − 1| ≤ ε^θ M^{1−θ}, with θ a harmonic measure.
  * For a disk Ω = D(1, ρ), θ = log ρ / log(ρ max(C, 1)).
  * **This is a power of Δ_g.** The planar point is *interior* to the analyticity domain.
* **Einstein analytic structure [DV]** (`hyp_cauchy.py`).
  * x(n) has branch points only at n = ±i/√3 (spinodal points), so ρ = 2/√3 = 1.1547 > 1.
  * The Taylor series about n = 1 converges at n = 0 at the rate (√3/2)^K. Checked: error 1.7e-6 at K = 60.
  * Rates for ε = Δ_g⁻⁴: exponent 3.93 for C = 1, 2.41 for C = 1.1, 1.05 for C = 1.5.

### Milestone 8: the general response on the hyperbolic branch, the Taylor data and the bottleneck (THEOREM_H.md §§5–7)

* **Engine [DV].** `hyp_num.py` was run at r_h = 1, 5/4, 3/2, 2, 3, 5, 10 and agrees in exact rationals with
  * e₁[θI] = 4t₀ − 2f′t₁;
  * e₁[θJ] = x(−48 − 54x)t₀ + 6x f′ t₁.
* **Wald [D].** σ₁ = 4x t₀ for θI and 6x² t₀ for θJ.
* **General response.** For any F(I, J):
  δ(x) = 4(1−x)θ_I − 2(2+x)r_hθ_I′ + (−48x − 60x²)θ_J + 6x(2+x)r_hθ_J′, and 𝒱 = δ + (1−x)δ′.
  * All six tested operators (I, J, I², IJ, J², I³) reproduce the planar table at x = 1. J² and I³ are independent of the fit.
  * I gives 𝒱 ≡ 0, confirming GB-blindness.
* **Counting.** An operator with m Weyl factors enters first at Taylor order k = m − 3. So a_k is (k+3)-point vacuum data.
* **Exact CFT identity [D].** 𝒱(T₀) = 4π∂_Tη/∂_Ts|_{T₀} − 1 = −t₄/30 in any CFT₃. Hence |𝒱(T₀)| ≤ τ_*/30.
* **(a) and (b) collapse to one input [P].** Hypothesis (a) is the same as uniform smallness of the deviation for complex T near T₀.
  * Positivity (Lehmann/complete monotonicity, energy-correlator positivity, modular positivity) gives per-theory geometric bounds
    on the disk |n − 1| < 1. It does not reach n = 0 in the interior, and it never gives small deviations.
* **Bottleneck.** Uniform-in-k smallness of the modular-cumulant deviations, |δa_k| ≤ Δ_g⁻⁴ C^k with C < 2/√3:
  * k = 0 is controlled (t₄);
  * k = 1 is controlled (quartic, under the paper's quartic-positivity input);
  * k ≥ 2 is 5-point and higher, and uncontrolled.

  This would be a d = 3, all-order analogue of Hartman's d = 2 Rényi sparseness theorem.
  * Secondary [A]: the Fisher-zero-free analyticity domain around [0, 1].
  * With C ≤ 1, the rate would be Δ_g⁻⁴ log Δ_g.
* **Corrections recorded.** Lemma R3(iii) is withdrawn (THEOREM_R.md §1), following the collaborator's Remark 14(A).
* **Files (session 3).**
  * `hyp_engine.py`: helpers.
  * `hyp_test_einstein.py` → `out_hyp_test_einstein.txt`: Einstein massless check.
  * `hyp_num.py` → `out_hyp_num_*.txt`: the production engine, one run per r_h.
  * `hyp_coeffs.py` → `out_hyp_coeffs.txt`: the general response, 𝒱 and the Taylor data.
  * `hyp_cauchy.py` → `out_hyp_cauchy.txt`: the Einstein analytic structure and the Theorem H rates.
  * Superseded, too slow (unsimplified rational expressions): `hyp_full.py`, `hyp_fast.py`, and their `out_` files. `hyp_response.py`
    was never run.

## 2026-11-27, session 4: the k = 2 coefficient (THEOREM_K2.md)

* **Lemma M [P].** ⟨Ω|K_R^m A|Ω⟩ = ⟨Ω|A K_L^m|Ω⟩. So modular cumulants are left–right correlators across the Rindler horizon.
  * a₂ contains the vacuum 5-point cumulant ⟨O_Y O_Y K_L³⟩_c, a weighted light-ray (boost-charge) correlator.
* **Decomposition [DV + D].** δa₂ = 15λ_J (t₄) + 408λ_{I²} (quartic) − 2880β₅ (five-point) + O(λ²).
  * Only β₅ is uncontrolled.
  * Its exchange part is the (2,3) entry of a positive Gram matrix, so |β₅| ≤ √(G₂₂G₃₃).
  * Requirement: Theorem H needs only |β₅| ≲ Δ_g⁻⁴/2880. That follows from a **six-point diagonal c ≲ Δ_g⁻²** (natural size
    Δ_g⁻¹⁰), which is very weak but still 6-point.
* **Failures, with reasons.**
  * (a) 5-point light-ray superconvergence gives relations of the form β₅ = interference sum. The heavy side is not positive (support
    lemma).
  * (b) Tomita–Takesaki Δ^s only complex-boosts single operators, since B|Ω⟩ = 0, so it reaches 3-point data only [P]. The K_R
    positivity is per-theory only.
  * (c) Rényi/coherent states e^{−sK_R}|Ω⟩: energy-flux positivity there is one-sided or O(1), with no Δ_g.
* **Most promising handle [D-sketch, conjectural]: thermal microcausality ("thermal CEMZ").**
  * Type-D thermal backgrounds see I₂I₃, unlike type-N shocks. The momentum-enhanced delay is ∝ β₅ x³ k² for k up to Δ_g.
  * Two-sided causality constraints would give |β₅| ≲ Δ_g⁻². With uniformity in k, Theorem H would then give
    |4πη/s − 1| ≲ Δ_g⁻² log Δ_g.
  * Not done: the bulk delay computation in the vector/scalar channels, the CFT derivation through the thermal lightcone bootstrap,
    and uniformity.
* **Bottleneck.** Inclusive heavy weights in the multi-stress channels, G_{jj} = Σ_X |⟨X|T^j⟩|²/M², for j ≥ 3. Gram bounds then
  control every exchange off-diagonal.
