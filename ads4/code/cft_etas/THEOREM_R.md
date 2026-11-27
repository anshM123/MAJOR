# Real-time route: relaxation function, early-time OPE control, sign structure (session 2)

Status: 27 Nov 2026, research agent. Grades: [P] proof, [D] derivation, [D-sketch], [DV] numerically verified, [A] assumption.
Units: 2πT = 1 unless stated, so β = 2π. In the numerics, r₊ = 1 and 2πT = 3/2.

## 0. Objects

* **Normalized deviation relaxation function φ(t).** Define it by

  F(s) = E(s) − E_E(s) = ∫_0^∞ φ(t) e^{−st} dt,   so   4πη/s − 1 = F(0) = ∫_0^∞ φ(t) dt.

  Explicitly, φ(t) = (4/s_ent)·(Φ − Φ_E)(t), with Φ(t) = ∫ δη(ω) cos ωt dω the vacuum-subtracted Kubo relaxation function.
* **Retarded-kernel deviation k(t).** k(t) = −φ′(t), the zero-momentum commutator deviation. Then
  4πη/s − 1 = ∫_0^∞ t k(t) dt (integrate by parts, assuming φ decays).
* **Causal support.** The commutator ⟨[T_xy(t,x), T_xy(0)]⟩ vanishes for |x| > t. So k(t) involves only |x| ≤ t.

## 1. Lemma R3: early-time real-time control by the thermal OPE [D-sketch]

**(i) Convergence.** The thermal OPE g = Σ_O a_O (zz̄)^{Δ/2}·(structures) converges for |z|, |z̄| < 1, with
z, z̄ = (τ ± i|x|)/β (Iliesiu–Koloğlu–Mahajan–Perlmutter–Simmons-Duffin, 1802.10266). In real time, τ = it gives
|z|, |z̄| = |t ± |x||/β, so the expansion converges for |t| + |x| < β. By causal support, k(t) for t < β/2 is therefore given
entirely by the convergent OPE.

**(ii) Heavy states.** A single trace of dimension Δ ≥ Δ_g contributes at most |a_O| (2t/β)^Δ (up to structure-dependent
constants). If the OPE converges absolutely at some radius ρ′ < 1 [A], the heavy tail on t ≤ ρβ/2 with ρ < ρ′ is
≤ C(ρ′)(ρ/ρ′)^{Δ_g}. This is exponentially small in Δ_g.

**(iii) Analytic terms drop out.** Thermal-OPE terms analytic in x^μ, such as double traces at leading order with scalar
structures, have no discontinuity across the light cone. So they give no commutator, and no large-frequency (OPE) term in G_R.
Through this mechanism the independent thermal one-point data b₁ of [TT] (Remark 2) can drop out of the early-time commutator.
The surviving inputs are:
* t₄ and α_T;
* the [TT] anomalous dimensions times the factorized ⟨T⟩², i.e. ⟨TTTT⟩ data;
* the odd-power multi-stress terms [T^k], i.e. (k+2)-point vacuum data times ⟨T⟩^k;
* the heavy tail.

Status of (iii):
* **Proven** for scalar tensor structures, where |x|^{2n+J}C_J(τ/|x|) is a polynomial.
* **Open** for the spinning structures of T, which carry x̂^μ factors.
* **Holographic check [DV]** (`early_time_check.py` → `out_early_time_check.txt`). A deformation supported at depth u ≥ 1 − w
  changes F(s) only by an amount that decays **exponentially** in s. The decay rate 1.32 at w = 0.6 and 1.98 at w = 0.4 matches
  twice the tortoise depth of the layer (1.22 to 1.91 up to the width of the transition region). Consequences:
  * no large-frequency coefficient changes;
  * by bulk causality the commutator is unchanged for t < 2r_*(1−w);
  * deep (horizon) data, which is what b₁ encodes, sit only in "analytic" OPE data.

**Consequence.** On t ∈ (0, ρβ/2], the deviation k(t) is controlled by vacuum data (t₄, ⟨TTTT⟩, higher-point vacuum data)
plus an exponentially small heavy tail. This is *unlike* the Euclidean E(n), which also needs the non-factorized thermal
data. The open part is (iii) for spinning structures, and the higher-point vacuum data (I₂I₃).

## 2. Theorem R2: from early-time real-time control to the Kubo slope [P]

**Hypotheses.**
* **(A_α)** k extends analytically to the sector S_α = {|arg t| < α}, with |k| ≤ M on S_α, and |k(t)| ≤ M e^{−γt} for real
  t > 0. This is the QNM gap γ together with a bounded QNM sum.
* **(B)** |k(t)| ≤ ε for t ∈ (0, t₂], with t₂ = ρβ/2. This is the output of Lemma R3.

**Conclusion.** For every t > 0 in S_α,

  |k(t)| ≤ M^{1−θ(t)} ε^{θ(t)},   θ(t) = (2/π) arcsin((t₂/t)^p) for t ≥ t₂,   p = π/(2α),

and |4πη/s − 1| ≤ ∫_0^∞ t·min(M^{1−θ(t)}ε^{θ(t)}, M e^{−γt}) dt.

**Proof.** The two-constants theorem gives |k| ≤ M^{1−ω}ε^{ω}, where ω is the harmonic measure of the slit (0, t₂] in
S_α \ (0, t₂]. Compute ω through two maps:
1. ζ = t^p sends S_α to the right half-plane and the slit to (0, b], with b = t₂^p.
2. ξ = iζ followed by ξ ↦ √(ξ² + b²) sends the upper half-plane minus the vertical slit [0, ib] onto the upper half-plane. The two
   sides of the slit go to [−b, b], and the real point c goes to i√(c² − b²).

The harmonic measure of [−b, b] at ih is (2/π) arctan(b/h). With h = √(c² − b²) this gives θ = (2/π) arcsin(b/c). ∎

**Rate [D].** Let L = log(M/ε). The minimum switches from the first bound to the second at
t_c = ((2/π) t₂^p L/γ)^{1/(1+p)}, which gives

  |4πη/s − 1| ≲ M t_c² exp(−γ^{p/(1+p)} ((2/π) t₂^p L)^{1/(1+p)}).

**Einstein values [DV]** (`qnm_sector.py` → `out_qnm_sector.txt`):
* The zero-momentum QNM angles below the real axis are 55.2°, 57.3°, 58.1°, tending to 60°.
* So k(t) = Σ r_n e^{−iω_n t} is analytic and bounded in |arg t| < α for any α < 55.2°, which means p ≥ 1.63.
* With α → π/3, p = 3/2, t₂ → π and γ = 1.776: the exponent is ≈ −2.34 L^{0.4}. With ε = Δ_g⁻⁴ this gives a rate of
  exp(−c (log Δ_g)^{0.4}).
* This beats every power of log, but it is not a power of Δ_g.

**Comparison with Theorem C.**
* Theorem C: Euclidean data on N, rate (log)^{−2γ}.
* Theorem R2: Lorentzian data on an interval, rate exp(−c (log)^{1/(1+p)}).
* What is new is that the input (B) is a **real-time, causal, OPE-convergent** quantity. By Lemma R3 it does not need the
  non-factorized thermal data b₁ (holographic evidence; still open for spinning structures).
* The remaining real-time hypothesis (A_α) is O(1). It is not small.
* Sharpness: two-constants bounds are attained by exp(log ε·(ω + iω̃)), so the rate cannot be improved under (A_α) + (B)
  alone [D-sketch].

## 3. Theorem R1: sign structure gives a power of Δ_g [P]

Let φ = φ₊ + r, with φ₊ ≥ 0 (or ≤ 0 throughout) and ∫|r| ≤ ϱ.

**(a) Support.** If supp φ₊ ⊂ [0, t_L], then

  |F(0)| ≤ e^{t_L}(|F(1)| + ϱ) + ϱ.

Proof: ∫φ₊ ≤ e^{t_L} ∫φ₊ e^{−t}. Equality holds for φ₊ = δ(t − t_L), so the factor e^{t_L} is sharp.

**(b) Strip (Lyapunov).** If ∫φ₊ e^{γt} ≤ M (the strip bound at s = −γ), then

  ∫φ₊ ≤ M^{1/(1+γ)} (∫φ₊ e^{−t})^{γ/(1+γ)}.

This is Hölder's inequality for the positive measure φ₊ dt. The math_transfer agent has the same rate, ε^{γ/(1+γ)}.

**Physical version.** Take t_L to be the "stringy time", t_s = log(Δ_g²) in units 1/2πT. This is the Rindler redshift time
of the depth w_min ~ Δ_g⁻² at which bulk effective field theory breaks down. If F(1) = O(Δ_g⁻⁴), then
|4πη/s − 1| ≤ C Δ_g⁻². Under the strip version with γ = γ_E = 1.776 the bound is C Δ_g^{−2.56}.

**Where the sign could come from.** These statements are D-level, from the bulk Remark-8 class.
* The first-order Euclidean kernel is positive. The Euclidean on-shell action is ∫[A φ′² + y² B φ²] with both terms positive,
  so δa ≥ 0 pointwise gives F(s) ≥ 0 for every s > 0. This is a sign in s, which is not enough.
* The *depth measure* dν = −dδa(u) is what decomposes φ into layer responses: F(s) = ∫dν(w) K_s(w), with K_0 ≡ 1 and
  K_1(w) ≈ w. So:
  * a monotone running coupling (ν ≥ 0) gives F(0) ≤ F(1)/K_1(w_min);
  * with w_min ~ Δ_g⁻², this is the Δ_g² loss of (a).
* In the bulk EFT, a single leading operator (C³, i.e. t₄) has a monotone profile. Mixed signs enter only through subleading
  operators.
* The CFT image of "monotone running of the horizon coupling" is not known [open].
* `relax_sign.py` tests whether layer responses are sign-definite in t; see §5.

## 4. Assessment of the CFT handles requested

1. **Pole-skipping.**
   * The shear-channel special points are s = −n (ω = −i2πTn), the mirror images of the Matsubara points.
   * In the layer model F(−γ) ~ A w^{−γ}, so they are *amplified* by near-horizon structure, just as Matsubara data
     (F(n) ~ A wⁿ) are suppressed. For a single layer, F(0)² ≈ F(1)F(−1): DC transport sits between the Euclidean and the
     pole-skipping data, and Hölder (Theorem R1(b)) makes this exact for positive depth measures.
   * So a CFT bound on F at s = −1 is exactly hypothesis (S) with γ = 1.
   * But: (i) at k = 0 the shear channel at ω = −i2πT is not a pole-skipping point (the residue is generic), and (ii) no CFT
     observable is known to access the lower half-plane. Grozdanov's univalence bound uses s = −1 data at the shear
     pole-skipping momentum; it is holographic input [D-sketch; literature].
2. **Chaos.**
   * The chaos pole-skipping point is in the *energy* channel at ω = +i2πT, k = i2πT/v_B: the first Matsubara frequency at
     imaginary momentum.
   * It is tip-local at O(1). The bulk pole-skipping condition is the near-horizon vv-equation, so v_B is set by horizon
     couplings. This makes chaos the one known CFT-accessible observable that is O(1)-sensitive to the tip.
   * It probes a *different* horizon coupling: the shock-wave equation involves the entropy-type coupling ∂L/∂R_{uvuv} and the
     transverse ∂L/∂R_{uxux}. η/s instead involves κ_xy/κ_rt (Brustein–Medved).
   * Dimension count [D-sketch]: at linear order, δ(4πη/s) depends on three horizon jets (θ_I′(1), θ_J(1), θ_J′(1)). One
     tip-local number v_B fixes one combination. So chaos can fix η only if the two linear functionals happen to be aligned.
   * Testable next step: compute δv_B for θ_I I and θ_J J on the brane and compare the functional with
     −6θ_I′ − 108θ_J + 18θ_J′.
   * The Δ_g⁻² stringy shift of λ_L is a UV (Regge) effect. It is not the EFT shift that governs η.
3. **QNEC, ANEC, light rays.**
   * These constrain one-point functions or entropies in states: ⟨T_kk⟩ ≥ (S″)/2π, ⟨E⟩ ≥ 0.
   * η is a two-point dissipative coefficient. In a homogeneous thermal state the ANEC and QNEC are saturated trivially, or are
     statements about ⟨T⟩.
   * The coherent-source version (the collaborator's) is the right way to make ANEC dynamical, and its causality input in
     coherent states is unproved. [No new result.]
4. **Modular flow / the equivalence principle.**
   * At N = ∞, the time-band algebras of the thermal GFF form half-sided modular inclusions (Leutheusser–Liu). Borchers then
     gives a positive generator P with Δ^{it}PΔ^{−it} = e^{2πt}P: the emergent near-horizon (Kruskal–Rindler) kinematics.
   * This structure is **kinematic**. It holds for every spectral density with full support, so it holds equally for the
     Einstein and the deformed (layer) theories. It therefore does not fix the horizon coupling normalization that η/s measures.
   * Precise failure: the CFT "equivalence principle" available from modular theory fixes boosts, not couplings.
   * What is missing is a *dynamical* statement: the low-frequency Liouvillian sector is close to the vacuum boost sector *with
     the same normalization as c_T and s*. This is the human's relative form bound restricted to that sector.

## 5. Sign tests [DV, `relax_sign.py` → `out_relax_sign.txt`]

(Filled in when the run completes; see the README ledger.)
