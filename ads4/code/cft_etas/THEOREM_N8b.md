# Session 10: six points under 𝒩 = 8. The J² ↔ W₊³W₋³ sector, what it means for the Taylor coefficients, and the honest 𝒩 = 8 statement

Status: 27 Nov 2026, research agent. Grades: [P] proof, [D] derivation, [D-sketch], [DV] numerically verified, [Q] quoted from literature, [A] assumption.

## 1. MHV at n points: the bound D ≥ 2n + 8 [D; DV for n = 5]

**Setup.** Write the n-point local MHV superamplitude as δ¹⁶(Q)·f.
* f must be S_n-symmetric, with little-group weight t_i⁻⁴ on every leg.
* Every component ⟨ij⟩⁸f must be a polynomial.

**Argument.**
* f = P₁₂/⟨12⟩⁸ = P₃₄/⟨34⟩⁸, with both P polynomial.
* For n ≥ 5, generic kinematics with ⟨12⟩ = 0 still has ⟨34⟩ ≠ 0. So f has no pole at ⟨12⟩ = 0, and by symmetry no pole at all.
  Hence **f is a polynomial**.
* Weight t⁻⁴ on every leg means 4 net square brackets per leg. That forces dim f ≥ 2n, so D = 8 + dim f ≥ 2n + 8.
* At n = 4 the argument fails, because ⟨12⟩ = 0 forces ⟨34⟩ = 0. That is why R⁴ (D = 8) exists.

**Checks.**
* n = 5 gives D ≥ 18. The exact count of session 9 is 0, 0, 0, 0 at D = 10–16 and 1 at D = 18 [DV].
* In EFK language this is "D^{2k}Rⁿ excluded at MHV for k ≤ 3" [Q, arXiv:1003.5018], since D = 2n + 2k ≥ 2n + 8 means k ≥ 4.
* **n = 6 MHV: the first independent contact has D ≥ 20** (D⁸R⁶).

## 2. NMHV at six points, where J² lives [Q]

**Why NMHV.** The J² = (W₊³ + W₋³)² term contains 2W₊³W₋³, whose six-graviton component is (−−−+++). In η-degree counting that
is (8, 8, 8, 0, 0, 0) = 24, i.e. **NMHV**. In 3 → 3 crossing it is the **elastic** channel (+++ → +++). So the NMHV six-point
contact is exactly the six-point *diagonal* (THEOREM_Q5, THEOREM_RB).

**Literature result.** Elvang–Freedman–Kiermaier (arXiv:1003.5018, §4.1; review arXiv:1012.3401, §6.2.1) [Q]:
* NMHV **R⁶ (D = 12) and D²R⁶ (D = 14) have no 𝒩 = 8 supersymmetrization**. They construct a six-point NMHV superamplitude from 5
  basis matrix elements and show by a complex shift that the pole factor cannot cancel.
* NMHV **D⁴R⁶ (D = 16) exists**. They construct it explicitly; it is the 7-loop candidate.
* The authors note that the k = 1 exclusion (D²R⁶) may be where their argument stops rather than a sharp bound. D⁴R⁶ is
  constructed, so the first six-point NMHV contact is at **D ≤ 16**, and exactly 16 if their k ≤ 1 exclusions are complete.

**Not done.** I did not count the NMHV sector myself. A full superspace count needs the EFK solution of the NMHV Ward identities,
and it is beyond this session. The weight bound of §1 does not apply to NMHV, where f carries Grassmann structure.

**So:** the bare J² contact (R⁶, D = 12) is forbidden by 𝒩 = 8. The first independent six-point (elastic 3 → 3) datum is D⁴R⁶
(D = 16).

## 3. Translation to the Taylor coefficients of 𝒱 [D-sketch]

**Bookkeeping.** In gauged 𝒩 = 8, a superinvariant of dimension d contains terms ∂^{d−2j}g^{2j} with g ~ 1/L. An n-point invariant
contributes to Weyl monomials C^m with m ≥ n. At first order in the couplings, C^m enters the Taylor series of 𝒱 first at
a_{m−3}.

| datum | n | d | 4D Weyl term | first Taylor order | natural size (IIA, Δ_g ~ λ^{1/4}) | fixed by |
|---|---|---|---|---|---|---|
| R⁴ | 4 | 8 | C⁴ (I², P²) | a₁ | Δ_g⁻⁶ | localization |
| D⁴R⁴ (incl. its C⁵ at g²) | 4 | 12 | C⁴, C⁵ | a₁, a₂ | Δ_g⁻¹⁰ | 4-pt Mellin data (partly localization) |
| D⁶R⁴, D⁸R⁴ | 4 | 14, 16 | C⁴, C⁵, C⁶ (g-dressed) | a₁ … | Δ_g⁻¹², Δ_g⁻¹⁴ | 4-pt data [A: unique completion] |
| **D⁴R⁶ (NMHV)** | 6 | 16 | C⁶ (J²-type) at g⁴ | **a₃** | **Δ_g⁻¹⁴** | **independent** |
| **D⁸R⁵** | 5 | 18 | C⁵ (IJ-type) at g⁸ | **a₂** | **Δ_g⁻¹⁶** | **independent** |
| D⁸R⁶ (MHV) | 6 | 20 | C⁶ at g⁸ | a₃ | Δ_g⁻¹⁸ | independent |
| n ≥ 7 (NMHV D⁴Rⁿ?; N²MHV at n ≥ 8 not known) | ≥ 7 | ≥ 18 | C^{≥7} | a_{≥4} | ≤ Δ_g⁻¹⁶ | independent / unknown |

**Reading the table.**
* **a₀ = 0** exactly (t₄ = 0).
* **a₁** is fixed by four-point data: R⁴ (localization) + D⁴R⁴ + ….
* **a₂** is fixed by four-point data up to c₁₈ (D⁸R⁵), at Δ_g⁻¹⁶.
* **a₃** is fixed by four- and five-point data up to c₁₆⁽⁶⁾ (D⁴R⁶), at Δ_g⁻¹⁴. **That is the first genuinely independent datum,
  ordered by natural size.**

## 4. The 𝒩 = 8 statement, and how large p can be (question 3)

**With natural sizes for the independent superinvariant coefficients [A], plus assumptions (i)–(iii) below:**

  planar 4πη/s − 1 = 27γ_{R⁴}  [localization-fixed, ~Δ_g⁻⁶]
                     + (D⁴R⁴, D⁶R⁴, D⁸R⁴ contributions)  [four-point Mellin data, Δ_g⁻¹⁰ … Δ_g⁻¹⁴]
                     + O(c₁₆⁽⁶⁾ ~ Δ_g⁻¹⁴, c₁₈ ~ Δ_g⁻¹⁶, …).

* **Localization-fixed up to O(Δ_g⁻¹⁰):** this is p = 10 relative to the 27γ_{R⁴} term.
* **Four-point-fixed up to O(Δ_g⁻¹⁴):** p = 14.

**Without the natural-size assumption, no p is justified.** The independent coefficients c₁₆⁽⁶⁾, c₁₈, … are **unbounded** by any
proven CFT statement, and an O(1) value in AdS units would shift η/s at O(1). Supersymmetry moves the obstruction from D = 10
(β₅) to D = 16 (c₁₆⁽⁶⁾). It does not remove the need for a bound.

**Assumptions list.**
* (i) Bulk 𝒩 = 8 gauged supergravity (or the IIA/11D lift) plus local superinvariants, valid for the planar thermal state.
* (ii) The flat-space superinvariant counts carry over to AdS₄, with unique gauged completions of the four-point invariants.
* (iii) Theorem H's analyticity input (b). In the supergravity regime it holds with radius 2/√3.
* (iv) The counts: 5-point (session 9) [DV]; 6-point MHV bound [D]; 6-point NMHV [Q, EFK, not re-derived]; n ≥ 7 and N²MHV not
  counted.
* (v) A natural-size bound on the independent coefficients, **if** a value of p is to be quoted.

## 5. How one would bound the first independent coefficient (question 4) [D-sketch]

* **The six-point datum c₁₆⁽⁶⁾ (D⁴R⁶, NMHV) is the elastic 3 → 3 contact**, so it is the natural candidate for positivity.
  * At tree level, heavy exchange in the three-particle channel contributes Σ_X |⟨X|TTT⟩|²·(kernel) ≥ 0. So a forward 3 → 3
    dispersion relation would make it **sign-definite**, as unitarity makes α₄ ≥ 0 for R⁴.
  * An **upper** bound needs the 3 → 3 superconvergence (CHMRS) sum rule, whose obstructions are listed in THEOREM_Q5 §3. Maximal
    supersymmetry should help: tree amplitudes have improved large-z behaviour, and superamplitudes package the forward-limit
    polarization sums. That has not been worked out.
* **The five-point datum c₁₈ (D⁸R⁵, MHV) is a 2 → 3 interference.**
  * It is **not sign-definite**, for the same reason as β₅.
  * It is bounded by the Gram inequality |c₁₈| ≲ √(R⁴-type diagonal × six-point diagonal at the matching order). So it follows from
    the six-point bound.
* **Does the R⁴ positivity argument extend to D⁸R⁵?** No. R⁴'s sign comes from a 2 → 2 elastic forward amplitude. D⁸R⁵ is a
  five-point (inelastic) contact, and the analogous elastic statement lives at six points.

## 6. Summary

* Under 𝒩 = 8:
  * the first independent six-point (J²-sector) contact is D⁴R⁶ at D = 16 [Q];
  * the first independent five-point contact is D⁸R⁵ at D = 18 [DV];
  * MHV n-point contacts need D ≥ 2n + 8 [D].
* So the planar deviation is localization-fixed to O(Δ_g⁻¹⁰), and four-point-fixed to O(Δ_g⁻¹⁴), **if** the independent
  coefficients have natural size.
* Unconditionally, no p is justified. The single remaining input is a bound on the elastic 3 → 3 contact D⁴R⁶. That would also
  bound D⁸R⁵ through Gram, and then Theorem H closes.
