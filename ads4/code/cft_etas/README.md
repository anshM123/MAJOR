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
