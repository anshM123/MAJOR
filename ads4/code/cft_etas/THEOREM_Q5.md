# Session 6: can the quintic coupling be bounded two-sidedly? A negative answer, and what is left

Status: 27 Nov 2026, research agent. Grades: [P] proof, [D] derivation, [D-sketch], [DV] numerically verified, [A] assumption.
Script: `quintic_family.py` → `out_quintic_family.txt`.

## 1. Proposition NQ: four-point data plus every known positivity or causality cone leave β₅ and η/s free [D at tree level; DV]

**The family.** Take Einstein gravity plus heavy scalars X with couplings X(g₂I + g₃J) at tree level. Integrating out X gives
F = (g₂I + g₃J)²/(2M²), i.e. (a, b, c) = (g₂²/2M², g₂g₃/M², g₃²/2M²).

What fixed four-point data leave free:
* Hold the four-point data fixed: g₂ (equivalently the quartic a) and everything in ⟨TT⟩, ⟨TTT⟩, ⟨TTTT⟩. At tree level the graviton
  2 → 2 amplitude contains X only through g₂.
* Then g₃ is a free real parameter, and

  4πη/s − 1 (planar) = (432g₂² − 1944g₂g₃ − 3240g₃²)/M² = 216(2g₂² − 9g₂g₃ − 15g₃²)/M².

* This is bounded **above** by (3618/5)g₂²/M² ≈ 724 g₂²/M², attained at g₃ = −0.3g₂. That upper bound is controlled by four-point
  data.
* It is **unbounded below**, so KSS can be violated by any amount.
* On the hyperbolic branch the Taylor data are a₀ = 0 and a₁ = −672g₂²/M², which see only g₂. g₃ enters first at
  a₂ = 12g₂(17g₂ − 240g₃)/M² (5-point) and a₃ = 9(23g₂² + 464g₂g₃ + 832g₃²)/M² (6-point).

**Constraints the family satisfies.**
1. **2 → 2 positivity, null constraints and crossing**, including the paper's sum rules and Caron-Huot–Van Duong: these involve
   g₂ only [P]. Note that CHVD's "g₃" is a four-point Wilson coefficient (the s³ term of the 2 → 2 amplitude). Their two-sided
   confinement acts on the heavy 2 → 2 spectral moments. It never involves couplings to three-graviton states.
2. **Gram / Cauchy–Schwarz cone b² ≤ 4ac** (Remark 4): saturated, 4ac − b² = 0.
3. **Higher-point positivity.** The family is a genuine unitary tree exchange, so its 4-, 5- and 6-point coefficients form a PSD
   (rank-one) matrix. That is exactly the content of multipositivity, g₅² ≤ g₄g₆:
   * Cheung–Remmen, "Multipositivity bounds for scattering amplitudes", arXiv:2505.05553, PRD 112, 016017 (2025);
   * the six-point preprint arXiv:2609.17861, which reports that higher-point positivity is obstructed by sign-indefinite cubic
     couplings and holds only in a sector that forms PSD matrices;
   * the precursor Chandrasekaran–Remmen–Shahbazi-Moghaddam, arXiv:1804.03153.

   None of these bounds a five-point coefficient by four-point data alone.
4. **Thermal causality** (session 5): K(q) = (2g₂ − 3qg₃)²/4M² ≥ 0.
5. **CEMZ.** X does not change the three-graviton vertex at tree level.
6. **Planar CFT order.** All couplings are O(1) in planar normalization. Unitarity at fixed c_T bounds them only at O(c_T).

**Not checked.** The family has not been checked against Regge boundedness of the **five- and six-point** amplitudes.
* A finite set of scalars already fails four-point Regge boundedness (the 4-derivative vertex gives s³). So any consistent
  completion needs a tower, and the paper's sum rules assume one.
* Whether five- and six-point Regge consistency of such a tower confines g₃/g₂ uniformly is unknown. No theorem in the literature
  covers it, as far as these searches found.
* **This is the only escape left.**

**Conclusion.** Under every constraint now available, the KSS-violating direction is controlled only by the six-point diagonal
c = Σ_X g₃²/2M². The one-sided upper bound on η/s *is* controlled by four-point data.

## 2. Route (A), crossing confinement: why it does not transfer [P/D]

* The confinement mechanism in the literature (2 → 2 null constraints) acts on ratios of 2 → 2 heavy spectral data.
* The quintic needs the ratio ⟨X|T³⟩/⟨X|T²⟩. The XTTT coupling is absent from every tree-level 2 → 2 graviton amplitude, so no
  2 → 2 null constraint can see it [P].
* At five points, crossing of a tree exchange is automatic and the heavy spectral weights are the interference g₂g₃, which has no
  sign. Five-point null constraints are therefore linear *relations*: β₅ equals a heavy interference sum. They are not bounds.
* At six points the diagonal Σg₃² is positive, which leads to route (B).

## 3. Route (B), a 3 → 3 dispersion relation: obstructions [D-sketch]

The heavy side, Σ_X |⟨X|T³⟩|²/M^k, is positive (optical theorem). But four obstructions stand in the way:
1. **Disconnected and spectator terms.** Forward 3 → 3 contains 2 → 2 subprocesses with a spectator. These are larger by 1/G (in
   CFT, O(c_T⁻¹) against O(c_T⁻²)), and they carry graviton t-channel poles 1/t → ∞ in exactly forward kinematics.
   * A functional must annihilate them. By the support lemma it then kills the interference terms too, but it *can* keep the
     positive diagonal. So this obstruction is not fatal for the diagonal.
2. **Analyticity.** Fixed-angle or fixed-invariant analyticity in s₁₂₃ for 3 → 3 is not established. Multi-particle thresholds and
   Landau singularities intervene, unlike the 2 → 2 Martin/Mandelstam domain.
3. **Subtractions.** The low-energy side of the s₁₂₃ dispersion contains c itself. A bound needs an *improved* (superconvergent)
   sum rule in which c drops out and only gravity-exchange terms remain, as CHMRS do at four points. That requires a
   Regge bound on 3 → 3 amplitudes in the s₁₂₃ channel. Nothing like that is known, even in flat space.
4. **Partial-wave positivity.** Positivity for three-particle states is a matrix statement over the multiple 3-body angular
   structures at fixed J (cf. Degrande–Li–Xu, arXiv:2511.15524, for 2 → N unitarity). Positivity holds, but the functional must be
   designed in that larger space.

**Status.** A 3 → 3 CHMRS would give c ≲ G²Λ⁻¹⁰-type bounds. Then the Gram cone would give |b| ≤ 2√(ac) ≲ Δ_g⁻⁸, far below what
Theorem H needs (Δ_g⁻⁴/2880). The weakest useful version is c ≲ Δ_g⁻². This is the cleanest target left, and it is a genuine
multi-particle bootstrap problem.

## 4. Route (A′), mixed correlators ⟨TTTX⟩ and ⟨TXTX⟩ with heavy X [D-sketch]

**Per species.**
* The contact g₃(X) appears linearly in the TT → TX amplitude. Its s-channel threshold is s = M², so it is analytic on [0, M²)
  apart from graviton poles.
* If that amplitude is Regge-bounded (a Cauchy–Schwarz/Rindler argument bounds it by ⟨TTTT⟩ and ⟨XTTX⟩), a three-subtracted sum
  rule gives g₃(X) = Σ_Y λ_{TTY}λ_{TXY}/m_Y⁸.
* Cauchy–Schwarz then gives |g₃(X)|² ≤ A·W(X), with A = Σ_Y λ²_{TTY}/m⁸ (four-point, controlled) and W(X) = Σ_Y λ²_{TXY}/m⁸.
* W(X) is bounded by the CHMRS-type gravitational sum rule for TX → TX.
* So per X, |g₃(X)| is controlled at natural size [D-sketch].

**Obstruction.** β₅ needs the *inclusive* Σ_X g₃(X)²/M_X² ≤ A·Σ_X W(X)/M_X². That is a double sum Σ_{X,Y} λ²_{TXY} over heavy
pairs: heavy-heavy-light data whose growth is set by the density of heavy states (Hagedorn-like in string completions). It is not
controlled by light four-point data.

## 5. Route (C), structure of the hyperbolic continuation [P]

At first order in couplings, 𝒱(x) = Σ (coupling) × (fixed polynomial in x). So 𝒱 is entire in x for *every* coupling vector, and
analyticity on D(1, 2/√3) imposes no relation between a₁ and a₂. In the family, a₁ = −672g₂²/M² is fixed while a₂ ranges over ℝ.
Unitarity saturation and full S-matrix bounds give O(1), not Δ_g-small, constraints at planar order (§1, item 6).

## 6. Verdict and next target

* **[D] No combination of four-point vacuum data with positivity, crossing, causality or multipositivity bounds β₅, or bounds the
  KSS-violating direction of η/s.** The explicit tree-level counterexample family is §1.
* The only unexamined constraint is **Regge boundedness of the five- and six-point amplitudes of the heavy tower**.
* Two equivalent targets remain:
  * a 3 → 3 superconvergence/CHMRS bound on the six-point diagonal;
  * a bound on the inclusive heavy-heavy-light weights Σ_{X,Y} λ²_{TXY}, via route A′.
* Either would give |β₅| ≲ Δ_g^{−p}. The weakest useful input is c ≲ Δ_g⁻², i.e. Σ_X ⟨X|T³⟩²/M_X² ≲ Δ_g⁻².
* Corollary worth recording: within the exchange sector, the four-point data *do* bound the **increase** of η/s, by
  4πη/s − 1 ≤ (3618/5)g₂²/M² + O(Δ_g⁻⁴). The open problem is entirely the lower (KSS) side.
