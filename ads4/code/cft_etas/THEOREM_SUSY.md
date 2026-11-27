# Session 8: maximally supersymmetric holographic CFT₃s (3d N = 8, N = 6 ABJM)

Status: 27 Nov 2026, research agent. Grades: [P] proof, [D] derivation, [D-sketch], [DV] numerically verified, [Q] quoted from literature, [A] assumption.
Scripts: `weyl4_reduce.py` → `out_weyl4_reduce.txt`.

## 1. What supersymmetry fixes

1. **Three-point data: τ = 0 exactly [Q + D].**
   * Buchbinder–Kuzenko–Samsonov (arXiv:1503.04961): for 1 ≤ 𝒩 ≤ 3 the supercurrent three-point function is fixed up to one
     overall coefficient.
   * Hence ⟨TTT⟩ equals the free-multiplet value. A free 𝒩 = 2 chiral multiplet has 2 real scalars and 1 Dirac fermion, so
     t₄ ∝ (2·4 + 2·(−4))/4 = 0.
   * Any 𝒩 ≥ 2 theory contains the 𝒩 = 2 subalgebra and its stress multiplet, so **t₄ = 0** there.
   * The 𝒩 = 1 analysis (arXiv:2102.04827) finds no parity-odd structure, so α_T = 0. (Jain–Nizami 2205.11157 discuss parity-odd
     𝒩 = 2 invariants for general currents; for T I follow the 𝒩 = 1 statement.)
   * Hence τ = 0, and by the identity of THEOREM_H §5, **𝒱(T₀) = a₀ = −t₄/30 = 0 exactly**.
2. **Four-point data [Q].**
   * The 𝒩 = 8 stress-multiplet four-point Mellin amplitude, through the order that matters here, is supergravity + R⁴ + higher
     D^{2k}R⁴. These come from Chester–Pufu–Yin 1804.00949, Binder–Chester–Pufu 1808.10554 and 1906.07195, and Alday–Chester–Raj
     2107.10274, 2207.11138.
   * The R⁴ coefficient is fixed by supersymmetric localization (the mass-deformed S³ free energy). In ABJM this is known as a
     function of N/k⁵, and in the M-theory limit D⁴R⁴ is absent.
3. **Five- and six-point data.**
   * Flat space [Q]: maximal supersymmetry forbids R⁵, and D²R⁴ vanishes on shell (Elvang–Freedman–Kiermaier 1003.5018). Also
     Rⁿ, D²Rⁿ, D⁴Rⁿ and D⁶Rⁿ with n > 4 are excluded as independent MHV superinvariants.
   * So the I₂I₃ ≅ IJ quintic (10 derivatives) and the J² sextic (12) have **no independent supersymmetric completion** at their
     natural orders.
   * AdS₄ [D-sketch]: an 11D or 10D superinvariant of fixed derivative order, reduced on S⁷ or CP³, produces 4D terms with *equal
     or fewer* derivatives (internal curvature only lowers the count). So 4D 10-derivative five-point terms come only from 11D terms
     with ≥ 10 derivatives:
     * in M-theory, ℓ_p¹²D⁶R⁴, since D⁴R⁴ is absent;
     * in IIA, α′⁵ζ(5)D⁴R⁴ and up.

     Their five-point parts are completions of four-point superinvariants, so they are fixed by four-point data. Size:
     β₅ ≲ Δ_g⁻¹⁰ in the IIA planar regime, Δ_g ~ λ^{1/4}.
   * **No 3d 𝒩 = 8 five-point stress-multiplet correlator has been computed.** The only computation found is the AdS₅ analogue
     (Rastelli–Zhou 1906.05305). The CFT-side statement, that 𝒩 = 8 Ward identities fix the five-point contact data from four-point
     data, is therefore **open**.
4. **No scalar-hair sector.** 𝒩 = 8 SO(8) gauged supergravity has no SO(8)-singlet scalar: the 70 scalars sit in 35_v ⊕ 35_c. The
   KK towers on the round S⁷ are in SO(8) representations [n,0,0,0], with a singlet only at n = 0. So no light singlet scalar couples
   to TT, and the Remark 12(a) sector is absent [D]. The IIA/CP³ dilaton needs a separate check.

## 2. Check against the literature: the R⁴ correction to η/s for M2-branes [DV]

`weyl4_reduce.py`: in 4D, the t₈t₈-type Weyl⁴ invariant
W = C^h_{mnk}C_p^{mnq}C_h^{rsp}C_q^{rsk} + ½C^{hk}_{mn}C_{pq}^{mn}C_h^{rsp}C_q^{rsk}
reduces **exactly** to

  W = (1/32)(I² − P²)   (least squares over random 4D Weyl tensors, residual 3e-11).

Our planar response table assigns I² → +864 and P² → 0. So for 16πG𝓛 = R + 6 + γW in AdS₄ units:

  **4πη/s − 1 = 864/32 · γ = 27 γ / R_AdS⁶.**

* Garousi–Ghodsi (arXiv:0808.0411, §3.2, the 4D-Weyl scheme for AdS₄×S⁷/Z_k) find η/s = (1 + 1728γ)/4π, with their L the S⁷
  radius.
* In their metric P du² ∝ L²du²/(4u²f), so R_AdS = L/2. Then 1728γ/L⁶ = 1728γ/(64 R_AdS⁶) = **27γ/R_AdS⁶**. This agrees with ours
  exactly.
* They obtain the same 1728 in the 11D-Weyl scheme. Their type IIA result is likewise positive and ∝ ζ(3)λ^{−3/2} [Q].
* So our general response formula reproduces the known M2-brane R⁴ correction. The P² piece gives 0 to η, and to η_H as well,
  consistent with parity.

## 3. Theorem H for maximally supersymmetric theories [conditional]

**Taylor data of 𝒱.**
* δa₀ = 0 exactly (t₄ = α_T = 0).
* δa₁ = −1344 λ_{I²}, with λ_{I²} = γ/32, which is fixed by localization (R⁴).
* δa_k for k ≥ 2 = (kinematic descendants of R⁴, from V_{I²}(x(n)): 408λ_{I²}, 414λ_{I²}, …)
  + (D^{2k}R⁴ four-point superinvariants, O(Δ_g⁻¹⁰) and below)
  + (genuinely ≥ 5-point superinvariants, which first appear at much higher derivative order).

**Statement [D given (i) and (ii)].** For a maximally supersymmetric holographic CFT₃ in the planar, large-gap limit (IIA ABJM,
λ ≫ 1), assume:
* (i) the bulk is 4D 𝒩 = 8 gauged supergravity (or its 10D/11D lift) plus higher-derivative superinvariants;
* (ii) the superinvariant classification of §1.3 carries over to AdS₄.

Then

  **4πη/s − 1 = 27 γ_{R⁴} + O(Δ_g⁻¹⁰),   4πη_H/s = 0,**

where γ_{R⁴} ~ ζ(3)λ^{−3/2} ~ Δ_g⁻⁶ is the localization-fixed R⁴ coefficient, and the O(Δ_g⁻¹⁰) remainder is fixed by four-point
D⁴R⁴ data. The deviation is therefore **positive**: KSS is respected, and the theory sits on the Remark 4 cone at its natural
apex.

**In terms of Theorem H.** Supersymmetry supplies exactly the missing all-order input (a). Every Taylor coefficient is determined by
four-point (Mellin) data up to the order of the first independent higher-point superinvariant. (b), analyticity in n, holds in the
supergravity regime, with radius 2/√3 as in Einstein gravity.

## 4. What is still assumed (question 4)

* **Not purely CFT.** The step "five-point contact data are fixed by four-point data" uses the flat-space superinvariant
  classification, lifted to AdS₄ by a derivative-counting argument [D-sketch]. That is bulk EFT plus supersymmetry, with no bulk
  F(Weyl)-power counting assumption (Proposition R is not needed).
* **The purely CFT version would need:**
  1. 3d 𝒩 = 8 superconformal Ward identities for the five- and six-point stress-multiplet correlators, showing that their contact
     ambiguities at the relevant Mellin orders are absent or fixed. Not available in the literature found.
  2. Theorem H's analyticity input (b).
  3. The CFT definition of the H² observable and its Taylor data (THEOREM_H §§1–2).
* **Inputs that are genuinely CFT-proven:**
  * t₄ = α_T = 0, from superconformal Ward identities;
  * the R⁴ coefficient, from localization;
  * the identity 𝒱(T₀) = −t₄/30 [D], which is being independently checked in `code/hyp_identity`.

**Bottom line.** Under maximal supersymmetry the quintic obstruction disappears. The planar ratio becomes
4πη/s = 1 + 27γ_{R⁴} + O(Δ_g⁻¹⁰), with γ_{R⁴} fixed by localization. That is a Δ_g⁻⁶-accurate prediction controlled by
four-point data, conditional on (i) and (ii). For a fully CFT-native theorem, the precise remaining gap is the 𝒩 = 8 five-point
Ward-identity rigidity.
