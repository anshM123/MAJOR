# FREE_NULLPLANE: the regulated null-plane identity, tested exactly in the d = 3 free scalar and free Dirac fermion

*26 Nov 2026. Research agent, about 4 h. Unreviewed and not committed. No existing file was edited; BATCH4_FOLD.md and
`code/b4fold/` were not touched.*

* **Code.** New code is in `ads4/code/freenp/`:
  * `np1_kernels.py`: exact sympy for the light-ray kernels (scalar with any ξ, Dirac), pair terms, and the algebra;
  * `lfcore.py`: shared light-front one-body machinery;
  * `np2_onepart.py`, `np2b_small_eta.py`: I_η in one-particle packets, cutoffs, m-scaling, η → 0;
  * `np3_tstate.py`, `np3b_small_eta.py`: the same in ψ_T = T(F)|0⟩, including the two-body quadratic-form pieces;
  * `np4_k0.py`: the author's K⁰ estimate in the original (celestial) frame;
  * `np5_spectral.py`: endpoint measure, endpoint rule and family bound.
* **Outputs.** They go to `code/freenp/cache/`, which ignores everything but its own `.gitignore`.
* **Compute.** One or two processes, `ulimit -v 3000000`, `OMP_NUM_THREADS=1`. Wall times: np1 4 s, np2 21 min, np2b about 25 min,
  np3 8 min, np3b about 4 min, np4 11 s, np5 under 1 s.

Grades are as in ANGULAR_MATCH: [P] exact or proved; [DV] derived and checked numerically; [NE] numerical evidence; [D] derived by hand.

---

## 0. Verdict

| question | answer | grade |
|---|---|---|
| **Does the regulated identity hold in free CFT₃?** | **Yes, exactly.** For every profile and state tested, I_η increases to ½⟨ψ, A((f′)²)ψ⟩ = 2‖H^{1/2}Aψ‖² − 2Re⟨Aψ, AHψ⟩. The approach is from below with error −c·η² (c > 0). | [P] algebra; [NE] limits to 1e-10 relative (one particle) and 1e-8 (ψ_T) |
| **Sign and factor** | [A, [H, A]] = **+½** A((f′)²), with H ≥ 0 the u-translation generator, O(u) = e^{iHu}Oe^{−iHu}, ds² = −du dv + dy². The ½ is the 1/(2·mass) of a Galilean particle of mass 2K. | [P] |
| **Extra contact terms in [𝓝, 𝓔]?** | **None.** The commutator is exactly Córdova–Shao's canonical form, (1/i)[𝓝(y), 𝓔(z)] = −δ(y−z)∂_z𝓔 + δ′(y−z)𝓔. There is no δ‴, δ″ or c-number (Schwinger) term, for the scalar at every ξ and for the Dirac fermion. The residual R of §18 is identically 0 (a = 0). | [P] |
| **m-scaling** | Null plane, fixed ψ: ‖H^{1/2}A(f_m)ψ‖² = m²⟨P_K⟩ + O(m), with ⟨P_K⟩ = ⟨A(1)⟩. For one particle this is an exact quadratic in m. Original frame (author's ψ_m of width ε/(1+m)): ‖(K⁰)^{1/2}E_mψ_m‖²/‖ψ_m‖² = κ(ε)(1+m)²/ε² with κ ≈ 1.50–1.55 (κ → 1.50 as ε → 0). So the author's **O(m²) holds**, with the constant set by the packet's own K⁰. | [P] null plane; [NE] celestial |
| **Spectral defect** | In any free CFT₃, E\|T⟩ lies in the two-particle space, which is a sum of **twist-1** conserved-current families (plus ψ̄ψ, at twist 2 and J = 0). The band τ ≳ √J is empty, so **μ((L,∞)) = 0 for every L and 𝒟 = 0**. The endpoint rule is satisfied exactly, **3/2 + 3/2 = 3**. The family bound Σ_OΔ²\|a_{O,N}\|² = O(N²) is **saturated**: it equals (N+3)². | [P] |

**What this does and does not establish.** In free theories every null-plane light-ray operator is a one-body operator. 𝓔 is the
mass density and −𝓝 the momentum density of a 1+1-dimensional Galilean system with particle mass 2K = 2(−p_v). The identity is that
system's **longitudinal f-sum rule**. Its only input is "momentum density = mass current", which is exactly
[N(g), A(f)] = iA(gf′). The free test therefore confirms:
* the conventions (sign, ½, Córdova–Shao form);
* that improvement terms and the constrained fermion component do not spoil the algebra in d = 3;
* that the zero-mode sector (p_v = 0) does no harm on physical states;
* that the regulator detects real non-domain situations (§4);
* that the author's O(m²) and O(N²) estimates are correct and **sharp**.

It cannot test the genuinely non-perturbative content. In an interacting CFT that content is whether the light-ray current algebra
survives, which in Galilean language is whether interactions are "density–density". Free theories satisfy this trivially. A
structural by-product (§3.3): for |m| ≤ 1 the profiles f_m are quadratic polynomials, A(f_m) is a global conformal charge, and the
identity holds non-perturbatively in **any** CFT. The open statement concerns |m| ≥ 2 only.

---

## 1. Conventions and the light-front dictionary

* **Coordinates.** u = t − x, v = t + x, so ds² = −du dv + dy², g_uv = −½, g_yy = 1. The null plane is u = 0, and its generators
  ∂_v point along +x.
* **Massless momentum p.** K := −p_v = (p⁰ − p^x)/2 ≥ 0 is the v-momentum. h := −p_u = (p⁰ + p^x)/2 = q²/(4K) ≥ 0, with q = p_y. A
  plane wave is e^{ip·x} = e^{−ihu − iKv + iqy}. H is the u-translation generator, H ≥ 0, and on one particle it equals h.
* **Measure.** The covariant measure d²p/(2p⁰) equals dK dq/(2K). There are no constrained modes to solve for: the null-plane fields
  are the covariant fields restricted to u = 0, so **this is light-front quantisation with no extra assumptions**.
* **Zero modes.** K = 0 at fixed q ≠ 0 means p⁰ → ∞ (UV, Schwartz-suppressed). K → 0 with q² ~ K means finite-energy particles moving
  **along the null generators** (the forward direction +x). These are the only "zero modes" a physical state can populate.

## 2. Commutators: exact one-body form, no contact terms (`np1_kernels.py`) [P]

**2.1 Kernels.** ⟨p₁|T_μν(x)|p₂⟩ = M_μν e^{i(p₂−p₁)·x}, and ∫dv sets K₁ = K₂.
* *Scalar*, T_μν = ∂φ∂φ − ½g(∂φ)² − ξ(∂∂ − g□)φ²:
  * ⟨T_vv⟩ = 2K² and ⟨T_vy⟩ = −K(q₁+q₂).
  * The **ξ-terms drop**: they are 2ξΔ_vΔ_μ with Δ_v = K₁ − K₂ = 0.
* *Dirac* (γ⁰ = iσ₂, γ¹ = σ₁, γ² = σ₃, ψ̄ = −ψ†γ⁰, u(K,q) = (√(2K), q/√(2K))):
  * ūγ_v u = −2√(K₁K₂), the good component only.
  * ūγ_y u = (K₁q₂ + K₂q₁)/√(K₁K₂), good × constrained.
  * Hence **⟨T_vv⟩ = 2K² and ⟨T_vy⟩ = −K(q₁+q₂), identical to the scalar.** In d = 3 the transverse space is one-dimensional, so the
    spin-density term ε_AB∂_B of d = 4 has no analogue.
* *Antiparticles*: the same, by C.

**2.2 One-body operators.** With normalised modes b and b̃(K,y) = (2π)^{−½}∫dq e^{iqy}b(K,q):
* 𝓔(y) = ∫dK K b̃†b̃(K,y);
* 𝓝(y) = (i/2)∫dK [b̃†∂_y b̃ − (∂_y b̃†)b̃];
* H = ∫dK dq (q²/4K) b†b.

So on each particle, **a(f) = K f(y), n(g) = −½{g(y), p} and h = p²/(4K)**, a Galilean particle of mass 2K in one transverse
dimension.

**2.3 Pair (vacuum) terms.** ∫dv puts the pair terms of both theories on δ(K₁+K₂), i.e. on exactly-forward collinear pairs. With
θ the angles from +x, K₁+K₂ = (r₁θ₁² + r₂θ₂²)/4 is positive definite. The pair matrix elements vanish there:

| operator | pair matrix element on forward collinear pairs |
|---|---|
| T_vv (scalar, every ξ) | O(θ⁴) |
| T_vy (scalar, every ξ) | O(θ³) |
| T_vv (Dirac) | O(θ⁴) |
| T_vy (Dirac) | O(θ³) |

Hence ∫d²p₁d²p₂/(p₁⁰p₂⁰) δ(K₁+K₂)Mψ ∝ ∫r dr δ(cr²)r^k = 0 (k ≥ 3).
* **𝓔\|0⟩ = 𝓝\|0⟩ = 0** on every state with a bounded wavefunction near the forward set.
* So **⟨0\|[N, A]\|0⟩ = 0: no Schwinger term.**
* A pair term would need a stress-tensor matrix element that is non-zero on exactly collinear pairs, i.e. k = 0.
* *Caveat* (not a contact term): finite-segment operators ∫_{−L}^{L}dv T_vv have **divergent** vacuum fluctuations at finite L
  (∫dq_rel is unbounded). Contact structure must therefore be studied on states, not by truncating the light ray.

**2.4 Algebra.** These are exact kernel identities for all (Q, k, l), checked in sympy with plane-wave profiles:
* [n(g), a(f)] = i a(gf′);
* [h, a(f)] = (i/2) n(f′);
* [a(f₁), [h, a(f₂)]] + (1↔2) = a(f₁′f₂′).

The commutator of two second-quantised one-body operators is the second quantisation of their one-body commutator, so in Fock space:
* **[N(g), A(f)] = iA(gf′)**;
* **[H, A(f)] = (i/2)N(f′)**, which is also the conservation law ∂_u𝓔 = ½∂_y𝓝 used in §18;
* **[A(f), [H, A(f)]] = ½A((f′)²)**;
* [A(f), A(g)] = 0.

In local form, (1/i)[𝓝(y), 𝓔(z)] = −δ(y−z)∂_z𝓔(z) + δ′(y−z)𝓔(z). The general contact residual R = aδ‴ − 2a′δ″ + a″δ′ of §18
is **identically zero**: the commutator kernel equals the canonical one at every (Q, k, l), so no contact term of any order is present.
This is consistent with the author's uniqueness argument, which forces a = 0. Free fields do not test that argument's Ward/Jacobi
hypotheses, but they are a model in which all of them hold.

**2.5 The f-sum-rule reading.** Take ρ_f = Σ_i K_i f(y_i) = A(f) and masses m_i = 2K_i. Then
[ρ_f, [H, ρ_f]] = Σ K_i²f′(y_i)²/m_i = ½A(f′²). This is the Thomas–Reiche–Kuhn (longitudinal f-sum) rule. In an interacting CFT
the analogous identity holds exactly when the light-ray algebra keeps j = p/m, i.e. [N, A] = iA(gf′). That, not anything about
free fields, is the content of §18's assumption 2.

## 3. The bound in concrete states

**Method.**
* *Profiles.* F = Re or Im of f_m(y) = (1+y²)e^{2im arctan y}, with f_m′ = 2(y + im)e^{2im arctan y}; optionally a C^∞ cutoff
  e^{−(y/R)⁸}, or a sharp cutoff.
* *Exact evolution.* G(t) = ⟨ψ, [A, A_t]ψ⟩/i with A_t = e^{iHt}Ae^{−iHt}, evaluated by FFT per K-slice. **The commutator algebra is
  not used.**
* *Regulator.* I_η = −∫ρ_η′G dt, with ρ_η Gaussian (even, ρ̂ > 0), evaluated by 80-point (ψ_T: 30-point) Gauss–Hermite.
* *Comparison values.* W := ½⟨A(F′²)⟩, 2‖H^{1/2}Aψ‖², 2Re⟨Aψ, AHψ⟩, and the author's spectral form
  2∫λρ̂(ηλ)dν_v − 2Re⟨v, AHψ⟩.

**Reading of the task's formula.**
* What converges is **lim I_η = ½⟨A((f′)²)⟩**, which equals 2‖H^{1/2}Aψ‖² − 2Re⟨Aψ, AHψ⟩.
* Equivalently, 2‖H^{1/2}Aψ‖² = 2Re⟨Aψ, AHψ⟩ + ½⟨A((f′)²)⟩, the §17 form.

### 3.1 One-particle packets (`np2_onepart.py`, `np2b_small_eta.py`) [NE]

**States** (Gaussian bumps in (K, q)): S1 has K₀ = 0.6, q₀ = 0.3; S2 has K₀ = 0.2 ("light", fast transverse motion); S3 is S1
displaced to y₀ = 1.5. There are 48 K-nodes, N_y = 16384 and dy ≈ 0.015.

| state, F | W = ½⟨A(F′²)⟩ | 2‖H^{1/2}Aψ‖² | 2Re⟨Aψ,AHψ⟩ | I_η: η = 0.1, 0.025, 0.0125 | (I−W)/η² |
|---|---|---|---|---|---|
| S1 Re f₀ | 26.6666666667 | 52.7470416667 | 26.0803750000 | 26.6666666667 (all η) | 0 |
| S1 Im f₁ | 1.2000000000 | 3.3000000000 | 2.1000000000 | 1.2000000000 (all η) | 0 |
| S1 Re f₂ | 26.0777347842 | 45.0208191095 | 18.9430843254 | 26.0519, 26.0759, 26.0773 | −2.88 |
| S1 Re f₄ | 25.0283578933 | 33.3709097942 | 8.3425519009 | 23.8907, 24.8810, 24.9891 | −257.3 (η → 0) |
| S1 Im f₄ | 20.8383087734 | 38.5761318725 | 17.7378230991 | 19.7002, 20.6909, 20.7990 | −251.5 |
| S2 Re f₂ | 4.8109240424 | 7.7177339957 | 2.9068099533 | 4.7168, 4.8020, 4.8086 | −12.5 |
| S2 Im f₄ | 6.3622147594 | 10.4622901710 | 4.1000754116 | 5.4009, 5.9883, 6.2258 | −873 (not yet asymptotic) |

* **Quadratic-form identity.** 2‖H^{1/2}Aψ‖² − 2Re⟨Aψ,AHψ⟩ − W ≤ 4e-14 in every case.
* **η → 0** (np2b, S1 Re f₄, η = 1/80 … 1/1280):
  * I − W = −3.93e-2, −9.99e-3, −2.51e-3, −6.28e-4, −1.57e-4;
  * (I − W)/η² → −257.3, a clean η² approach from below. S2 is in `cache/np2b_small_eta.txt`.
* **The author's spectral form** 2∫λρ̂(ηλ)dν_v − 2Re⟨v, AHψ⟩ tracks I_η to O(η²) and has the same limit, confirming §18's
  regulator theorem.
* **m = 0, ±1.** I_η is **exactly η-independent**: G(t) is exactly linear in t (§3.3).
* **m = 8** (np2b): I − W = −1.08e-2 (S1) and −4.30e-2 (S2) at η = 1/1280; the S1 η² coefficient is −1.77e4.
* **Onset scale.** For |m| ≥ 2 the approach sets in only for η ≲ η_* ~ K/m², because Aψ carries transverse momenta ~ m and so
  h ~ m²/K. The limit is **not uniform in m**. This is harmless for the author's argument, which takes η → 0 at fixed m, but any
  attempt to take m → ∞ jointly must respect it.

### 3.2 The stress-tensor state ψ_T = T(F)\|0⟩ (`np3_tstate.py`, `np3b_small_eta.py`) [NE]

**The state.**
* ψ(p₁,p₂) = F̂(p₁+p₂)·[−¾(ζ·p₁)(ζ·p₂) + ⅛((ζ·p₁)² + (ζ·p₂)²)], with ζ = (0,1,i)/2 (helicity +2, improved T; same polynomial as
  `code/kernels/fock_free.py`).
* F̂ = exp(−\|P−P₀\|²/2σ²), P₀ = (1,0,0), σ = 0.25: a T wavepacket near rest.

**Light-front subtlety.** ψ_T necessarily populates K → 0 (one particle moving along +x).
* Measured: the slice density of ⟨A(y²)⟩ is ∝ K^{−½}, with √K·density = 8.0150, 8.0157, 8.0196, 8.0308 at K = 8.7e-6, 2.4e-4, 1.4e-3,
  4.9e-3.
* The transverse spread grows like K^{−½}.
* The quadrature is K = s², 28 Gauss–Legendre nodes in s, with per-K grids up to N_y = 65536.
* For a profile growing like \|y\|^α, every term of the identity has K-density ∝ K^{3/2−α}. **f ~ y² (α = 2) is integrable with
  margin ½**; α ≥ 5/2 would diverge. The conformal weight (1+y²) is therefore safe, but not by much.

**Method.**
* *Two-body pieces*: Schmidt decomposition of each (K₁,K₂) block, with operator matrices in the q-support basis; exact grid sums.
* *I_η*: one-body reduced density matrix. This is exact, because [A, A_t] is one-body.

| F | 2‖H^{1/2}Aψ_T‖² | 2Re⟨Aψ_T,AHψ_T⟩ | difference | W (two-body) | W (ρ₁) | I_η: η = 0.1, 0.05, 0.025 |
|---|---|---|---|---|---|---|
| Re f₀ | 649.7330 | 618.3062 | 31.4268164 | 31.4268164 | 31.4268164 | 31.42681635 (all η) |
| Re f₁ | 611.2645 | 579.8376 | 31.4268164 | 31.4268164 | 31.4268164 | 31.42681635 (all η) |
| Im f₁ | 17.0211 | 15.9311 | 1.0900717 | 1.0900717 | 1.0900717 | 1.09007171 (all η) |
| Re f₂ | 528.3117 | 497.4062 | 30.9055509 | 30.9055509 | 30.9055509 | 30.8396, 30.8831, 30.8985 |
| Im f₂ | 58.9002 | 54.0187 | 4.8815523 | 4.8815523 | 4.8815523 | 4.8152, 4.8590, 4.8745 |
| Re f₄ | 353.0778 | 322.5332 | 30.5445657 | 30.5445657 | 30.5445657 | 29.0577, 29.7466, 30.2162 |
| Im f₄ | 156.1131 | 137.7897 | 18.3233980 | 18.3233980 | 18.3233980 | 16.8346, 17.5249, 17.9949 |

* Smaller η (np3b; `cache/np3b_small_eta.txt`):

  | F | I − W at η = 1/80 | at η = 1/160 | at η = 1/320 |
  |---|---|---|---|
  | Re f₂ | −2.11e-3 | −6.21e-4 | −1.80e-4 |
  | Re f₄ | −1.136e-1 | −3.60e-2 | −1.10e-2 |
  | Im f₄ | −1.136e-1 | −3.60e-2 | −1.10e-2 |

  * The successive ratios 3.40 and 3.45 (m = 2) and 3.15 and 3.28 (m = 4) mean I − W ≈ −c·η^{1.7–1.8}, slower than the clean η² of
    the compact-K packets (S1 and S2 approach η² cleanly; S2 Re f₈ is still pre-asymptotic at η = 1/1280).
  * The cause is the forward (K → 0) particles, where h ~ q²/K is unbounded. The limit is unaffected: I_η → W monotonically from
    below.
* ⟨P_K⟩ = ⟨A(1)⟩ = 0.54504, ⟨A(y²)⟩ = 15.713.
* **Conclusion.** Every term is finite, the identity holds to ≤ 1e-8 relative, and liminf I_η = lim I_η = W < ∞ in ψ_T.

### 3.3 Why m = 0, ±1 are exact in any CFT [P]

f₀ = 1+y² and f_{±1} = (1 ± iy)² are quadratic polynomials. For quadratic f, A(f) = ∫dv dy f(y)T_vv is the charge of a conformal
Killing vector tangent to the null plane: P_v, a null rotation (y), and a null special conformal transformation (y²).
* So [H, A(f)] and [A, [H, A]] are fixed by the conformal algebra, with no light-ray input at all.
* A_t is a degree-2 polynomial in t, G(t) is exactly linear, and I_η ≡ ½⟨A(f′²)⟩ for every η. This is seen to 1e-13.
* In an interacting CFT the regulated statement for |m| ≤ 1 is therefore a domain statement about conformal generators, not an open
  operator-algebra question.
* **The open content of §18's assumption 2 is the |m| ≥ 2 profiles**: the genuinely light-ray (supertranslation-type) operators.

### 3.4 Cutoffs: smooth is fine, sharp diverges (`np2_onepart.py`, `np2b_small_eta.py`) [NE]

**Smooth cutoff** f_{2,R} = f₂e^{−(y/R)⁸}, R = 1, 2, 4, 8 (S2): I_η → W at every R, from below. Example R = 8: W = 3.0330,
I_{0.0125} = 3.0305.

**Sharp cutoff** f₂·1_{\|y\|<2}: Aψ has a jump, so Aψ ∉ Dom(H^{1/2}).
* I_η = 0.257, 0.340, 0.389, 0.431, 0.480, 0.545 at η = 0.1 … 0.003; it grows without bound.
* The spectral form at very small η (dy = 0.005), η = 3.2e-3, 8e-4, 2e-4, 5e-5, 1.25e-5: 0.540, 0.748, 1.159, 1.988, 3.724.
  Successive ratios for η → η/4: 1.385, 1.550, 1.716, 1.873, tending to 2 = 4^{½}.
* Expected asymptotics I_η ~ η^{−½}: \|v̂(q)\|² ~ q⁻² up to q ~ (K/η)^{½}.
* **The regulator detects exactly the failure mode the author's Fatou argument is designed to exclude.** The cutoff f_{m,R} of
  §18 must be at least C¹; its uniform bound \|f′\|²/(1+y²) ≤ C(1+m²) presupposes that.

### 3.5 m-scaling [P/NE]

**Null plane, one particle** (exact; S1 fit residual 6e-14):

  ‖H^{1/2}A(f_m)ψ‖² = m²⟨A(1)⟩ + m⟨(1+y²)K j⟩ + ‖H^{1/2}A(f₀)ψ‖²

with j the transverse current. For S1: 0.6m² + 4.18m + 26.37. Also W_re + W_im = 2⟨A(y² + m²)⟩ exactly.

**Null plane, ψ_T** (complex A_m):

| m | 0 | 2 | 4 | 8 | 16 | 24 | 32 |
|---|---|---|---|---|---|---|---|
| ‖H^{1/2}A_mψ_T‖² | 324.87 | 293.61 | 254.60 | 253.16 | 389.13 | 560.68 | 801.90 |

* Successive slopes Δ‖·‖²/Δ(m²) = 0.536 (16→24) and 0.538 (24→32), tending to ⟨A(1)⟩ = 0.545. Two-particle cross terms decay in m.
* So **‖H^{1/2}A_mψ‖² = m²⟨P_K⟩ + O(m)**. The large constant (~ ⟨A(y²)⟩) comes from the forward, large-\|y\| region.

**Original frame, the author's statement** (`np4_k0.py`).
* *Setup:* ψ_m = T(F_m)\|0⟩ with width δ_m = ε/(1+m), ε = 0.15, and E_m = ∫e^{imθ}E(θ) → \|p\|e^{imθ_p} per particle. K⁰ = −\|p\|∇²_p
  on one-particle wavefunctions (measure d²p/\|p\|); this is fixed by [K⁰, P⁰] = −2iD and by the ground state e^{−\|p\|} having
  (P⁰+K⁰)/2 = ½ = Δ_φ.

| m | 0 | 2 | 4 | 8 | 16 | 32 |
|---|---|---|---|---|---|---|
| ‖(K⁰)^{½}E_mψ_m‖²/‖ψ_m‖² | 82.0 | 626.0 | 1709.5 | 5494.6 | 19536 | 73507 |
| ÷ (1+m)² | 82.0 | 69.56 | 68.38 | 67.83 | 67.60 | 67.50 |
| ⟨P⁰⟩⟨K⁰⟩/‖w‖⁴ | 80.2 | 627.6 | 1720 | 5542 | 19730 | 74283 |
| Δ² = (m+3)² | 9 | 25 | 49 | 121 | 361 | 1225 |
| target-family fraction | 1 | 0.999992 | 0.999986 | 0.999980 | 0.999975 | 0.999972 |

* **O(m²) confirmed.** At m = 8 the ε-dependence is 1.503/ε², 1.512/ε², 1.526/ε², 1.547/ε² for ε = 0.05, 0.10, 0.15, 0.20. The
  constant is dominated by the packet's own K⁰ ~ 1/δ_m².
* **Coercivity** ⟨P⁰⟩⟨K⁰⟩ ≥ Δ²‖w‖⁴ holds with a factor ≈ 1.5/ε² to spare.
* **Doppler/aberration loss** is ≈ 3e-5 at ε = 0.15, far inside the author's "≥ ½".
* Odd m: E_mψ ≈ 0 for back-to-back pairs (7e-3 at m = 3, from the packet width).

## 4. The spectral defect (`np5_spectral.py`) [P]

**Spectrum of E\|T⟩.** E(θ) is one-body and commutes with P, so E_n\|T(q),+2⟩ at rest is a two-particle state with angular
wavefunction ∝ e^{i(n+2)θ}(1+(−1)ⁿ). Odd n vanish. By Flato–Fronsdal, the two-particle states at rest with angular momentum j form
**one** family:
* *scalar:* J_{\|j\|} (Δ = \|j\|+1), or φ² (Δ = 1) at j = 0;
* *Dirac:* J_{\|j\|}, or ψ̄ψ (Δ = 2) at j = 0.

The normalised amplitudes are \|a_{O,n}\|² = 1 for the one allowed O and every even n, since free event shapes are δ-functions. This
matches `code/kernels/free_canonical.py` (\|â^{J_s}_{±s,2}\|² = 1).

| quantity | free scalar | free Dirac | author's statement |
|---|---|---|---|
| families in E\|T⟩ | J_s (twist 1), φ² | J_s (twist 1), ψ̄ψ (twist 2, J = 0) | — |
| R = T/(T−C) | −(J−1)/2 ≤ 0 | same; C = 0 at J = 0 | band τ ≳ √J needs R ≫ 1 |
| μ((L,∞)) | **0 for all L** | **0** | needs o(L⁻²), §16 |
| B_N | **0** (X_N = C − s(s+1) = 0 at every ladder top) | **0** | → 𝒟 |
| endpoint rule Σ[(T−b(b−1))/3 + 3/2]\|a_{O,b,0}\|² | J₄ (b = 4): 3/2; φ² (b = 0): 3/2; **sum 3** | J₄: 3/2; ψ̄ψ: 3/2; **sum 3** | = 3 + 𝒟 ⇒ **𝒟 = 0** |
| Σ_OΔ²\|a_{O,N}\|² | (N+3)² | (N+3)² | O(N²): **saturated**, constant 1 |
| fourth harmonic \|a_{J₆,4,2}\|² | 1 = \|c\|², t₄ = +4 | 1 = \|c\|², t₄ = −4 | ≥ \|c\|²: equality |

* **Endpoint rule.** The free check of its normalisation (the 1/3 and the RHS 3) is non-trivial. At the unitarity bound
  T_Δ = J(J−1) = b(b−1) for the endpoint family, so each endpoint contributes exactly 3/2.
* **Caveat.** This applies the long-multiplet formula at its degenerate boundary. The ladder itself is degenerate here: Y_k < 0, all
  lower helicities are null, and the whole amplitude sits at the ladder top.
* **The gravitational inequality fails** (\|c\|² = 1), as it must. The heavy sum is entirely twist-1 higher-spin currents, which §18's
  assumption 4 excludes, and φ² violates the Δ_φ > 1 assumption.

**What a non-trivial test of μ would take.** In **any** free CFT₃, E\|T⟩ stays in the two-particle sector, so μ on the band is always
empty. Populating τ ≳ √J needs multi-particle (multi-trace) components in E\|T⟩, i.e. interactions.

The minimal setting is the large-N critical O(N) model (or Gross–Neveu) at O(1/N).
* There E\|T⟩ acquires 4-particle components with twist ≥ 2 and unbounded twist.
* One needs ⟨TTTT⟩ at O(1/N) on the second sheet, decomposed into the [σσ]/[φφφφ]-type families, with endpoint amplitudes a_{O,b,0}.
* One then measures Σ_{τ≳√J, R>L}C\|a₀\|² against L⁻².

A cheaper proxy is the author's own sufficient condition, ‖(K⁰)^{½}E_mψ_m‖² = O(m²). In a free theory it can be computed for any
n-particle external state without a conformal-family decomposition, as in `np4_k0.py`. But the angular route needs the external state
to be T, which is two-particle in free theory.

## 5. Light-front subtleties: summary

1. **No constrained-mode or improvement contamination.** The ξ-improvement and the fermion's constrained component both enter
   T_vv and T_vy only through ∫dv ∂_v(…) or through K₁ = K₂, and drop out exactly.
2. **p_v = 0.** In physical states K → 0 means forward-moving particles (finite energy) or p⁰ → ∞ (Schwartz-suppressed).
   * Pair terms are supported on exactly-collinear forward pairs, where the matrix elements vanish as θ³ (𝓝) and θ⁴ (𝓔).
   * The forward region carries an integrable K^{−½} density, with transverse spread K^{−½}.
   * With the conformal weight y² this leaves margin ½ in the growth exponent.
   * A state singular at K → 0, i.e. with a non-vanishing collinear-forward amplitude at k = 0, would be needed to produce a pair or
     Schwinger term. Stress-tensor-generated states are not of that type.
3. **Finite light-ray segments** have divergent vacuum fluctuations. Regulate in t (as the author does), not in v.
4. **Non-uniformity in m.** The η → 0 limit sets in at η ≲ K/m².

## 6. Status for the record

* §18 assumption 2 (the regulated null-plane identity at finite c_T) is **verified exactly in both free CFT₃s**:
  * sign +, factor ½, Córdova–Shao canonical form;
  * no contact or Schwinger term;
  * liminf I_η = ½⟨A(f′²)⟩ < ∞ in one-particle packets and in ψ_T;
  * divergence only for non-C¹ profiles.
* The **|m| ≤ 1 part holds in any CFT** by the conformal algebra.
* The author's **O(m²)** (null-plane and K⁰ forms) and the **family bound O(N²)** are confirmed and **sharp**: free theories saturate
  them.
* **𝒟 = 0** in free theories, trivially: the band is empty. The endpoint rule's normalisation "3" is reproduced exactly.
* **Not tested** (impossible in free theory): whether interactions preserve [N(g), A(f)] = iA(gf′), i.e. the light-ray f-sum rule,
  at |m| ≥ 2; and any non-empty intermediate-band measure.
