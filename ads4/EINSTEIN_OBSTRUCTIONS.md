# Emergence of Einstein gravity in AdS4 from CFT3: the chain of implications and its missing links

Audit date: 16 November 2026. Scope: the target "P1" of `PROBLEM_CANDIDATES.md` — for unitary CFT3s with
`c_T → ∞`, large-N factorisation and a large higher-spin gap, bound the non-Einstein parameters `(t_4, α_T)`
of ⟨TTT⟩ by `K (log Δ_gap)^a / Δ_gap^b`. This file records the exact chain of implications used by the best
existing derivations, the status of each link, and a ranked list of tractable sub-theorems.

**Source marks.**

* **[R]** — PDF fetched from arXiv in this session, converted with `pdftotext`, and the quoted passages read;
  quotations are verbatim (line breaks and typography of formulas normalised). Equation numbers are those of
  the arXiv PDF version fetched.
* **[A]** — arXiv abstract only.
* **[M]** — from memory, not checked in this session.
* **[D]** — derived in this note (not found in the literature searched); unrefereed. Short derivations are
  given so they can be checked.
* **[N]** — numerical check run in this session (scripts reproduced in the Appendix); a sanity check, never a proof.

Status vocabulary as in `LITERATURE.md`: **theorem / derivation / numerics / conjecture / missing**.
"Proved" in a physics paper's own words is recorded as *derivation*.

**Corrections to `LITERATURE.md` found during this audit** (that file is not modified):

1. `LITERATURE.md` §(b) says no d = 3 statement with explicit bounds was found. There is one:
   Afkhami-Jeddi–Kundu–Tajdini, arXiv:1805.07393, §7, eq. (7.7) (quoted in §1.3 below).
2. The AdS4 scalar bound quoted there as "eq. (236)" is eq. (5.64) of the arXiv v1 PDF of 2106.10274.
3. "Shocks, superconvergence and a stringy equivalence principle" *is* arXiv:1904.05905 (Kologlu–Kravchuk–
   Simmons-Duffin–Zhiboedov, "KKSZ"); it is the same paper that contains the ANEC-commutativity argument.
   There is no separate 2008.xxxxx paper with that title (arXiv title search). The second author of
   arXiv:2311.04271 is Yakov Landau.

Papers read [R]: 0907.0151 (§2.4), 1407.5597 (§3.5), 1610.09378 (§5), 1707.03007, 1707.07689 (§1, §4.2),
1708.05718 (§2.1), 1712.04861 (§1, §3.1), 1805.07393 (§1, §2.2–2.3, §7), 1904.05892 (§1, §4), 1904.05905
(§1, §2, §4.1, §4.4–4.5, §5.1.1, §5.3, §6.1, App. A.3, E), 2106.10274 (§1.1.1, §2 summary, §5.3.2, §5.5,
App. B, F), 2201.06602 (§1–§3, §4.2, App. A, conclusions), 2311.04271 (§1, §2, §3.5.2, §4.6, §5),
2512.09986 (§3.4, §6), 2608.19322 (§5), 2610.24496 (abstract, §6 opening).
Abstract only [A]: 2205.01495, 2610.19149, 2512.13780, 2501.17949, 2408.06440, 2502.14963, 2609.20924.

Searches: INSPIRE citation lists (to Nov 2026) of 2311.04271 (15 citing records), 1904.05905 (161),
2201.06602 (206), 2106.10274 (141), 1805.07393 (56), 1707.03007 (48), 1904.05892 (81), read by title;
arXiv abstract searches "superconvergence", "higher-spin gap", "energy correlator" + "stress tensor".
"Not found" means not found by these searches.

---

## 0. Summary of findings

* The quantitative d = 3 statement exists in two forms, neither a theorem: (i) a Regge/chaos-bound estimate
  with exponents `(a, b) = (1, 4)` (1805.07393 (7.7)); (ii) the expected AdS4 uplift of the sharp 4D
  flat-space bound on `|ĝ_3|²`, which gives `(a, b) = (1/2, 4)` for `(t_4² + α_T²)^{1/2}` (2201.06602 (4.4)
  and conclusions). The two differ by `(log Δ_gap)^{1/2}`.
* No paper found (to Nov 2026) completes HPPS for stress tensors in any dimension. By their own account,
  2311.04271 lack: the remaining sum rules, the Regge regime, and the light-contribution dictionary. Their
  dictionary is written for generic bulk dimension `D ≥ 8` tensor structures; the `D = 4` (d = 3) reduction,
  including the parity-odd structure, is not written.
* The ANEC-commutator (superconvergence) sum rule ⟨T|[E,E]|T⟩ = 0 is written explicitly only for general d
  with scalar exchange and for d = 4 (KKSZ (5.134), (1.6)). The general-d formulas degenerate at d = 3
  (poles `1/(d−3)`; the transverse group is SO(2)). **The d = 3 sum rule is not written in the literature
  searched**; 2512.09986 (Dec 2025) and 2608.19322 (Aug 2026) list stress-tensor-source event-shape blocks and
  their superconvergence relations as future work.
* [D] In d = 3 the structure is simple enough to write down (§3.3): with `c := (t_4 + iα_T)/4` the helicity-flip
  amplitude of the energy one-point function, the stress-tensor contribution to ⟨T|[E(θ_1),E(θ_2)]|T⟩ is
  `2i |c|² sin 4(θ_1−θ_2) σ_3`, i.e. the quadratic form is exactly `t_4² + α_T²`, and the sum rule reads
  `|c|² + Σ_{O≠T} |â^O_{−2}|² = Σ_{spin ≥ 6} |â^O_{6}|²` with every term a square. This is an inequality of the
  requested type ("non-Einstein ⟨TTT⟩ ⇒ single-trace operators of spin ≥ 6 with ⟨TTO⟩ ≠ 0"), conditional on
  the same convergence/no-contact-term assumptions as KKSZ. The kernels are not yet computed.
* [D] Superconvergence alone is homogeneous in (non-Einstein data, heavy data): it cannot give a decay rate
  in `Δ_gap`. A rate needs a Regge-spin-2 sum rule containing the graviton pole (strategy S-A), or the
  unproved "no outliers" assumption of KKSZ §6.1.
* For S-A, the curved-regime positivity statements for external gravitons are not given in any source. By
  inference from the flat-space structure (2201.06602 (3.5), (A.11)) they are ν-moments of
  `P^{−μ}_{−1/2+iν}(cosh β)` with `μ ∈ {0, 4, 8}`; for the sum rules that bound `|ĝ_3|²` only `μ = 0` enters at
  leading order, i.e. exactly `Q_1 > 0`, `Q_2 < 0` of `POSITIVITY_LEMMA.md`. The `μ = 4, 8` moments are
  easier (closed forms in §4.3, [D], [N]).

---

## 1. The conjecture(s) in d = 3 as stated in the literature

### 1.1 General statements

**HPPS**, arXiv:0907.0151 §2.4 [R]: "Conjecture: any CFT that has a large-N expansion, and in which all
single-trace operators of spin greater than two have parametrically large dimensions, has a local bulk dual."

**CEMZ**, arXiv:1407.5597 abstract [R]: the extra graviton three-point structures lead to causality violation
which "cannot be fixed by adding conventional particles with spins J ≤ 2. But, it can be fixed by adding an
infinite tower of extra massive particles with higher spins, J > 2." For D = 4 (§3.5) [R]: "in four
dimensions the α_2 structure is identically zero … Thus, we have only the Einstein term structure and the α_4
one", plus "the parity violating structure"; the Einstein term gives a `log(L/b)` time delay and "we will
later repeat the computations for AdS4 space and we will see that L → R_AdS4". Their eq. (3.20) is the 2 × 2
phase-shift matrix `δ = 2Gs[−1·log|β| + (3γ/β⁴)|−⟩⟨+| + (3γ*/β*⁴)|+⟩⟨−|]` with complex `γ` "the coefficient
of the +++ amplitude".

**Meltzer–Perlmutter**, arXiv:1712.04861 [R]: eq. (1.2) `⟨TTT⟩ = ⟨TTT⟩_Einstein + α_2 ⟨TTT⟩_{R²} + α_4 ⟨TTT⟩_{R³}`,
"(In d = 3, one of these is replaced by a parity-odd structure.)"; eq. (1.3) "CEMZ showed that, indeed,
`α_n ≲ Δ_gap^{−n}`"; eq. (1.11) "Proposal: Counting AdS derivatives ⟷ Counting CFT powers of Δ_gap". For
d = 3: eq. (3.13) `Π_{iν_0}(L ≪ 1) ≈ −log L/(2π)` and footnote 9: "We leave the logarithmic beahvior in d = 3,
shown in (3.13), understood."

**Costa–Hansen–Penedones**, arXiv:1707.07689 eq. (1.4) [R]: "We would like to prove that
`C^{(2)}_{TTT} ∼ α_2 ∼ 1/Δ_g²`, `C^{(3)}_{TTT} ∼ α_4 ∼ 1/Δ_g⁴`", with `Δ_g ≡ Δ(J = 4)` on the leading Regge
trajectory; shown "with mild assumptions on the behaviour of OPE coefficients in the large gap limit".
For d = 3 their (4.3x): `t_4(0, L) = −48 L⁻⁴ ln(L)⁻¹ r_3(0)/r_1(0) (1 + O(1/ln L))`.

### 1.2 The d = 3 parameters

**Chowdhury–David–Prakash**, arXiv:1707.03007 [R]. Eq. (1.1): `⟨TTT⟩ = n_s^T ⟨TTT⟩_free boson + n_f^T ⟨TTT⟩_free
fermion + p_T ⟨TTT⟩_parity odd`. Eq. (1.5): energy matrix for stress-tensor states (polarisations `ε_xy = ε_yx = 1`
and `ε'_xx = −ε'_yy = 1`, detector in the y direction)

`Ê(T) = [[ (E/4π)(1 − t_4/4), α_T E/(16π) ], [ α_T E/(16π), (E/4π)(1 + t_4/4) ]]`,

eq. (1.6): `t_4 = −4(n_f^T − n_s^T)/(n_f^T + n_s^T)`, `α_T = 8π⁴ p_T / (3(n_f^T + n_s^T))`; positivity ⇒
`t_4² + α_T² ≤ 16`; large-N Chern–Simons-matter theories: `t_4 = −4 cos θ`, `α_T = 4 sin θ` (eq. (1.10)).

[D] In the helicity basis `|±⟩ = (|ε'⟩ ± i|ε⟩)/√2` the matrix (1.5) is `(E/4π)[[1, c̄],[c, 1]]` with

`c = (t_4 + i α_T)/4`, `|c|² = (t_4² + α_T²)/16 ≤ 1`,

and rotating the detector by `θ` multiplies `c` by `e^{4iθ}` (checked with sympy, Appendix A.3). This is the
CFT3 image of the CEMZ matrix (3.20) and of the complex 4D cubic coupling below.

**Bulk side.** 1805.07393 eq. (7.8) [R]:
`S = M_Pl² ∫ d⁴x √g [R − 2Λ + α_4 W_{μναβ}W^{μνρσ}W_{ρσ}^{αβ} + α̃_4 W̃_{μναβ}W^{μνρσ}W_{ρσ}^{αβ}]`,
"Coupling constants α_4 and α̃_4 are related to the coefficients ñ_s − 4ñ_f and ñ_odd respectively"; "the
four-derivative terms do not contribute to ⟨TTT⟩". 2201.06602 eqs. (2.6), (2.9)–(2.11) [R]:
`M(1⁺2⁺3⁺) = (ĝ_3/2) √(8πG) ([12][13][23])²`, `ĝ_3 = α_3 + i α̃_3` ("the real and imaginary part representing
parity-even and parity-odd couplings"). So `(t_4, α_T) ↔ (Re, Im) ĝ_3 / R_AdS⁴` up to a real positive constant
and a phase convention. **The numerical constant relating `c` to `ĝ_3/R_AdS⁴` was not found in the sources
read** (needed to turn a bulk bound into a number for `t_4² + α_T²`).

Note on normalisations: the Einstein point is `t_4 = 0`, i.e. `n_s^T = n_f^T` in the basis of 1707.03007 and
`ñ_s = 4ñ_f` in the basis of 1805.07393.

### 1.3 Quantitative statements in d = 3

**(C1) Regge/chaos-bound estimate.** 1805.07393 §7 [R]: "we will show that in d = 3, this simple relationship
between bulk derivative and Δ_gap has a logarithmic violation." Eqs. (7.5)–(7.6): positivity of the holographic
null energy operator at `ρ → 1` gives `ñ_s = 4ñ_f` and `ñ_odd = 0` at infinite gap. Eq. (7.7):

`|ñ_s − 4ñ_f| / c_T ≲ ln Δ_gap / Δ_gap⁴`,  `|ñ_odd| / c_T ≲ ln Δ_gap / Δ_gap⁴`,

"Hence, causality constraints translate into `α_4 ∼ ln Δ_gap/Δ_gap⁴`, `α̃_4 ∼ ln Δ_gap/Δ_gap⁴`." Hypotheses (2.25):
`√c_T ≫ Δ_gap ≫ 1` (as printed), `Δ_gap` the dimension of "the lightest single trace operator with spin ℓ > 2".
Status: **derivation with estimates**; see chain S-C in §2.3. Exponents `(a, b) = (1, 4)`.

**(C2) Sharp dispersive form.** 2201.06602 [R], eq. (4.4): `|ĝ_3|² M⁸ ≤ 24.9 log(M/m_IR) − 27.6`; (4.5):
`g_4 M⁶/(8πG) ≤ 12.3 log(M/m_IR) − 13.5`; conclusions: "As explained in [36], dispersive bounds for scalar
amplitudes in flat space can be lifted to AdS and imply corresponding bounds on OPE coefficients of
holographic CFTs. Physically, the same must be true for graviton scattering as well. For example, (4.4)
translates to a bound on three-point OPE coefficients of stress tensors in holographic 3d CFTs. In AdS4,
infrared logarithms are simply replaced by log(M R_AdS)." With `M R_AdS = Δ_gap` this is the conjecture

`t_4² + α_T² ≤ K² log(κ Δ_gap) / Δ_gap⁸`, i.e. `(a, b) = (1/2, 4)` for `(t_4² + α_T²)^{1/2}`,

with `K²` = 24.9 × (the unknown conversion constant of §1.2)². Status: **conjecture** (flat-space bound:
derivation + numerics; uplift: not done). The parametric CEMZ-type form in the same paper, eq. (3.28), is
`|ĝ_3| ≲ (1/M⁴) log[1/(b m_IR)]`, which matches (C1); the sharp form is stronger by `√log`.

**(C3) ⟨TTO⟩ for light single-trace O.** 1805.07393 eq. (7.14) (d = 3) [R]:
`C^{Even}_{TTψ}/√c_T ≲ ln Δ_gap/Δ_gap²`, `C^{Odd}_{TTψ}/√c_T ≲ ln Δ_gap/Δ_gap²`, bulk couplings (7.15)
`α_{Ψhh} ΨW² + α̃_{Ψhh} ΨW̃W`. For non-conserved spin 2: "We believe causality of the four-point function
⟨TTTT⟩, requires that ⟨TTO_{ℓ=2}⟩ must be suppressed by the gap for all single trace O_{ℓ=2}. However, a
detailed analysis requires the computation of ⟨TTTT⟩ using the conformal Regge theory which we will not
attempt in this paper." Flat-space sharp analogue, 2201.06602 (3.10): `|ĝ_3|² + F_{g3 matter} ≤ (1/M⁸)(37.8
log(M/m_IR) − 45.4)` with `F_{g3 matter} ≥ 0` a positive quadratic form in the couplings `|g_0(m_ℓ)|²`,
`|g_2(m_ℓ)|²` of light spin-0 and spin-2 states (eq. (3.9)).

**(C4) Contact terms in ⟨TTTT⟩.** 2201.06602: (4.5) above and "all dimensionless ratios of the form
`g_k M^{2(k−4)}/g_4` satisfy two-sided bounds … We do not have a formal proof for all couplings". 2311.04271
eq. (1.3) states the general-D target `|α_2| ≤ #/Δ_gap²` (in D = 4 the `α_2` structure is absent). Expected
CFT3 form: `g_4 ≲ (8πG) log(Δ_gap)/Δ_gap⁶` in AdS units. Status: **conjecture**.

**Which gap.** HPPS/1805.07393: dimension of the lightest spin > 2 single-trace operator. 1707.07689,
1712.04861: `Δ_g = Δ(J = 4)` on the leading Regge trajectory. 2106.10274 footnote 2: "Δ_gap is conceptually a
twist gap (not a dimension gap) for us", applied to *all* single-trace exchanges other than the listed light
ones; light single-traces of spin ≤ 2 are allowed ("as long as they all have spin J ≤ 2"). The hypothesis of
P1 (dimension gap for spin > 2 only) is therefore weaker than what strategy S-A uses; a theorem along S-A
needs a twist gap, or an extra statement excluding large-spin, small-twist single-trace operators.

---

## 2. The chains of implication

Notation: `E(n)` the energy-flux (ANEC) operator at direction `n` on the celestial circle; `L[O]` the light
transform; `J_0` the Regge intercept. "T-block" = contribution of the stress-tensor conformal family.

### 2.1 Strategy S-A: dispersive sum rules for ⟨TTTT⟩ with the bulk-point / flat-space dictionary

| Step | Exact statement | Source | Status | What a rigorous proof needs |
|---|---|---|---|---|
| **L-A0** | Lorentzian axioms (temperedness, Wightman positivity, analyticity domains) for the four-point function of T | Kravchuk–Qiao–Rychkov 2001.08778, 2104.02090 treat scalar four-point functions [A, via `LITERATURE.md`] | **theorem for scalars; missing for spinning operators** | Extension of the KQR theorems to conserved spinning operators (expected, not written) |
| **L-A1** | Regge boundedness of spinning correlators: `J_0 ≤ 1` nonperturbatively | 1904.05905 §4.2.4, eq. (4.73): "Using Rindler positivity, we have shown that J_0 ≤ 1"; footnote 29: "This argument has a small loophole discussed, e.g., in section 4.2.6" | **derivation** | Close the lightcone-limit loophole (asymptotic lightcone expansion on the second sheet); causal configurations with 3, 4 timelike ("we have not found a way to rigorously bound the behavior of the correlator") |
| **L-A2** | Dispersive sum rules `Ψ_{k,ν,f}[G] = 0` for stress tensors, labelled by Regge spin `k`, momentum `ν`, transverse spin `j′`: "for the integral of the functional to be convergent, we must have k > J_0 for even j′ and k − 1 > J_0 for odd j′ … So, we should choose k ≥ 2 for even j′ and k ≥ 3 for odd j′." Unsubtracted `k = 2, 3` sum rules for `J_e = 2`, eq. (3.80): `f⁺_{k=2;J_e=2} = −ν⁻²(2m² − ν²)`, `f⁻_{k=3;J_e=2} = 1` at the saddle; "the action on all double-traces will remain zero" | 2311.04271 (2.8), (3.54), §3.5.2 | **derivation** (general d) | d = 3 specialisation: the stabiliser groups `SO(d−3)`, `SO(d−2)` are trivial/degenerate, `j′` reduces to a parity label; parity-odd three- and four-point structures (97 parity-even permutation-invariant four-point structures, 5 functions of one variable on the `z = z̄` line: 1708.05718 §2.2) |
| **L-A3** | Completeness of the family of sum rules: only the lowest-subtracted flat-space sum rules are obtained (7 of 19 with even `k`, 3 of 10 with odd `k`): "The full set of flat space sum rules for gravitons includes additional sum rules that cannot be expressed in terms of 'shock amplitudes' … We leave the problem of obtaining these additional sum rules using CFT techniques to future work." | 2311.04271 §1, §4.6, §5 | **missing** (for the full ⟨TTTT⟩ statement) | For `(t_4, α_T)` alone the two families `B_2^{(1)}`, `B_3^{(1)}` of 2201.06602 (2.26) suffice (their footnote-10 functional (3.7) proves `|ĝ_3|² ≤ 1598 log(M/m_IR) − 3211`); these are "anti-subtracted" and correspond to (3.80). The identification in D = 4 helicity variables (f, g, h of (2.1)–(2.3)) is **not written**: the dictionary (4.86) is in the generic-dimension basis (4.81), "(D ≥ 8)" |
| **L-A4** | Heavy action in the bulk-point regime `Δ, ν ≫ 1`, `ν < m`, `J/Δ ≪ 1`: CFT block action = flat-space partial-wave action, eq. (4.86) with matrices (4.88) | 2311.04271 §3–§4 | **derivation** (saddle point: "scattering crystal") | Error terms uniform in `(Δ, J, ν)`; extension to `ν = O(1)` ("We expect the same to be true in the spinning case"); D = 4 |
| **L-A5** | Heavy action in the Regge regime `Δ, J → ∞` at fixed AdS impact parameter `β` | scalar case: 2106.10274 (2.59) `C_{k,ν}[Δ,J] = 2 P_{(2−d)/2+iν}(η_AdS)/m^{2k} × (1 + O(ν²/m², 1/m²))` | **missing for spinning operators**: "We leave an exploration of the Regge regime for four-point functions of spinning operators to future work." (2311.04271 §2.2.2) | Spinning Regge-moment formula in d = 3 for the `k = 2, 3` sum rules, as a Hermitian form in the ⟨TTO⟩ structures (2⁺ + 1⁻ for even `J ≥ 4`, 1⁻ for odd `J ≥ 5`: 1708.05718 table in §2.1) |
| **L-A6** | Heavy action when `Δ, J, ν` are all large, and near threshold `τ ≈ Δ_gap` | scalar: 2106.10274 App. F, (F.3)–(F.8) | **missing for spinning** | as L-A5 |
| **L-A7** | Flat-space positive functional: e.g. (3.3a) `F_{g3} = ∫_{m_IR}^{M} dp [(1−p)³ p (65 + 155p + 47p²) B_2^{(1)} + 2(1−p)³(81 − 1927p + 4601p²) B_3^{(1)}] − 5 ∂_{p²} B_4^{(1)}|_{p=0}` (units M = 1); "Without the cutoff m_IR, the functionals are positive for all states with m > M"; "With the cutoff m_IR, they only become negative at some large b ∼ m_IR^{−2/3}" | 2201.06602 §3.2; optimal (4.4) from numerical optimisation (App. C) | **numerics**: "How do we ascertain that a functional is positive on all heavy states? A simple strategy is to plot the action …" (J ≤ 300, m ∈ [1, 16]) plus the `m → ∞` scaling limit (3.5) | A proof of the two-parameter family of inequalities (§5, ST3). [N] my transcription of (3.3a) is positive on a grid; the printed 6-digit functional (3.7) is *not* positive as printed for `M ≤ m ≤ 1.05 M` in the +− channel (Appendix A.4) — rounding (`ψ_3(1) = 0.0307` instead of 0), so (3.7) cannot be used verbatim |
| **L-A8** | Positivity in the curved regime `β ∼ 1` | scalar: 2106.10274 (5.63) `lim m⁴ ω_AdS[Δ,J] ∝ a_1 Q_1(β)/M² + a_2 Q_2(β)/M³ + …`, "We claim that ω_AdS is rigorously positive at large M"; `Q_1 > 0`, `−Q_2 > 0` "by checking numerical examples" — now `POSITIVITY_LEMMA.md` | scalar: **theorem for the two leading moments** (Lemma Q), **derivation** for the reduction to them; gravitons: **missing** (depends on L-A5) | (i) L-A5; (ii) non-asymptotic control of the remainder (§5, ST2); (iii) for off-diagonal sum rules the `μ = 4, 8` moments of §4.3 |
| **L-A9** | Light contribution: `−ω|_light` = graviton pole (regulated by AdS) + `|c|²` term + light matter | flat space: 2201.06602 (2.31a) `−B_2^{(1)}|_low = 16πG/p² + 2πG|ĝ_3|²p⁶`, (2.31b) `−B_3^{(1)}|_low = −2πG|ĝ_3|²p⁴`, matter (3.8); scalar AdS: 2106.10274 (5.62) with `8πG/(ν² + (d/2)²)` | **missing for stress tensors in AdS**: "we leave exploration of light states to future work" (2311.04271 footnote 2) | The T-block (and light spin-0, spin-2 blocks) action of the `k = 2, 3` functionals in d = 3, exactly or with controlled `1/ν` corrections |
| **L-A10** | Control of the EFT tail | 2106.10274 "technical assumption" (5.3): `Σ O(g/R²_AdS) = O(1/R²_AdS)`; "it would be nice to remove this assumption" | **conjecture (assumed)** for scalar bounds; for the `k = 2, 3` stress-tensor sum rules "the EFT series of the corresponding flat space sum rule should automatically truncate" (2311.04271 §3.5.2), and (2.31): "The absence of contact term contributions is a hallmark of superconvergent sum rules" | Likely avoidable for `(t_4, α_T)` if only `k = 2, 3` are used; curvature corrections `O(1/R²)` to L-A9 still need a bound |
| **L-A11** | Large-N: double-twist families drop out (double zeros), multi-trace and loop contributions subleading | 2106.10274 §1.1.1 and §4; 2311.04271 below (2.11) | **derivation** | Uniform-in-`c_T` error terms; a definition of "single-trace" for a sequence of finite-`c_T` theories |
| **L-A12** | Conclusion `t_4² + α_T² ≤ K² log(κΔ_gap)/Δ_gap⁸ (1 + o(1))` | 2201.06602 conclusions (expected) | **conjecture** | all of the above + the constant `c ↔ ĝ_3/R⁴` |

Variant with a mixed correlator ⟨TTφφ⟩: requires a light scalar single-trace φ in the spectrum (not implied by
the hypotheses of P1), and positivity then involves the 3 × 3 system ⟨φφφφ⟩, ⟨TTφφ⟩, ⟨TTTT⟩; the existing
derivations with scalar probes are of Regge type (S-C below), not dispersive.

### 2.2 Strategy S-B: superconvergence / ANEC-commutator sum rules in the planar limit

| Step | Exact statement | Source | Status | What a rigorous proof needs |
|---|---|---|---|---|
| **L-B1** | `[E(n_1), E(n_2)] = 0` for `n_1 ≠ n_2`. Sufficient conditions (4.71)–(4.72): `J_1 + J_2 > max(2 − Δ′_0, 1 + J_0)`, `J_1 + J_2 > 2 − τ′_0 + |Δ_1 − Δ_2|`; for `T`: (4.79) "a product of ANEC operators on the same null plane is well-defined and commutative if 3 > J_0 = J_0^{dDisc}". Planar: (4.80) `J_0^{planar} ≤ 2` from the chaos bound; "commutativity can be lost at higher orders in large-N perturbation theory (and only recovered nonperturbatively)" | 1904.05905 §4 | **derivation** (the abstract says "We prove"). Stated gaps: "our analysis establishes well-definedness and commutativity (when applicable) acting on this dense subspace. However, this subspace does not include important states like momentum eigenstates"; "In principle we have not excluded the possibility of contact terms at z_1 ∝ z_2 … we will work under this assumption" | Commutativity as a statement about matrix elements between `∫ f T|0⟩` states including momentum-space wavepackets; exclusion (or classification) of contact terms at `n_1 = n_2`; the planar statement needs the chaos bound (itself a derivation) |
| **L-B2** | t-channel expansion (5.2)–(5.3) of ⟨O_4(p)|L[O_1]L[O_2]|O_3(p)⟩ over primaries; each ordered term has a fixed sign: "in a unitary theory, due to the reality properties of OPE coefficients, the contributions of operators have fixed signs" | 1904.05905 §5.1, §5.3 | **derivation**: "it is not obvious that this expansion converges. In this paper we will assume that it does"; "we will assume that no subtractions are necessary if O_1 = O_2 = T" | Convergence of the harmonic-projected expansion; it is automatic "for applications to planar theories with a finite number of single-trace exchanges" (footnote 38) |
| **L-B3** | T-block of ⟨T|E E|T⟩ = product of one-point transfer matrices; commutator quadratic in the non-Einstein couplings | d = 4: 1904.05892 (4.16), (4.22)–(4.24), `[M_1, M_2] = t_2²[…] + t_4²[…] − t_2 t_4[…] − (1↔2)` ⇒ `t_2 = t_4 = 0`; general d, `j_1 = j_2 = 2` component: 1904.05905 (5.133) ∝ `C_T((4+d)t_2 + 4t_4)²` | **derivation** (d ≥ 4) | d = 3: [D] §3.3, `[M(θ_1), M(θ_2)] = 2i|c|² sin 4θ_12 σ_3` |
| **L-B4** | Planar limit: only single-trace operators contribute: "The factors `Γ(4 − Δ_φ/2)⁻²` in (1.6) ensure that contributions of double-trace operators are suppressed by O(1/N⁴)" (general d: `Γ(d − Δ/2)⁻²` in (5.132)); 1904.05892 §4.2.1: double-trace corrections enter as `γ²_{n,l}/N⁴` | 1904.05905 §1, 1904.05892 | **derivation** | Uniformity over the infinite double-trace towers |
| **L-B5** | Event-shape blocks for an exchanged primary `(Δ, J)`: scalars, general d: (5.131)–(5.132) with `q(Δ) ≥ 0` | 1904.05905 §5.1–5.3; algorithm for general exchange in §5.1.3–5.1.4 | **derivation**; **missing in d = 3 for stress-tensor sources** — 2512.09986 §6: "A straightforward extension of the present work is the calculation of conformal blocks for external sources with nonzero spin"; 2608.19322 §5: superconvergence with "the external state … created by the stress-energy tensor itself" listed as future work | Light transform + Fourier transform of the d = 3 ⟨TTO⟩ structures, parity-odd included (§5, ST1) |
| **L-B6** | The sum rule. General d (5.134): `−C_T 2^{−1−3d} π^{3−d/2} ((4+d)t_2+4t_4)² Γ(d+3) / [(d−1)² d (d+1)² (d+4) Γ(3+d/2) Γ(d/2)²] + (1/4) Σ_φ |λ_TTφ|² q(Δ_φ) + non-scalar = 0` and a second equation with relative factor `(d² − 9)` and `(2/(d+2))`; d = 4: (1.6) | 1904.05905 | **derivation**; "non-scalar" terms not computed | d = 3: the general-d decomposition is singular ((5.128) has `1/((d+3)(d−3))`, (E.4) has `1/(d−3)`, `1/(d²−5d+6)`); [D] §3.3 gives the d = 3 structure |
| **L-B7** | Qualitative emergence: "if no single-trace operators are present other than the stress tensor, then t_2 and t_4 must vanish in the planar limit"; 1904.05892 "Result 1: ⟨T|[E_1, E_2]|T⟩ = 0 ⟹ a = c" under "CFT with large N, large gap ⟹ Only O appears in the E O OPE" | 1904.05905 §1, 1904.05892 (1.12) | **derivation**, d = 4 written; the heavy contributions "give small corrections suppressed by 1/Δ_gap" is asserted, not derived | d = 3: [D] (★) of §3.3 gives it with spin ≥ 6 only |
| **L-B8** | Quantitative suppression of the heavy contribution: (6.6)–(6.7) `|α²_GB(s)| ≤ a⁺_{J=3}(s) ∼ 1/(Δ⁴_gap c_T)`; CFT form (6.8) `C⁻(Δ = d/2 + iν, J = 3) = 0` | 1904.05905 §6.1: "we argue (non-rigorously)"; "In what follows we assume that there are no outliers"; "A quantitative understanding of the superconvergence sum rule requires some extra analysis which we postpone for future work" | **missing** | A bound on the heavy side of the sum rule. [D] It cannot come from commutator sum rules alone (§3.4); it requires the Regge-spin-2 sum rule with the graviton pole, i.e. S-A |

### 2.3 Strategy S-C (for reference): Regge limit + chaos bound, the source of (C1)

L-C1: in the Regge limit of smeared correlators, T exchange dominates when `1/Δ_gap ≲ s = 1 − ρ`; a heavy
spin-ℓ operator contributes `δG ∼ i e^{−sΔ_gap/2}/(σ^{ℓ−1} s^a)` (1805.07393 (2.26)); footnote 13: "δG has
large numerical factors … we are making an additional assumption that OPE coefficients which appear in δG
are small enough to cancel these large numerical factors." — **assumption**.
L-C2: holographic null energy condition (2.21) `i lim_{B→∞} ⟨ε·O(B) E_{r=√ρ B}(0) ε·O(B)⟩ ≥ 0`, from the
chaos/causality bound; "we conjecture that the smearing projects out [OO] double-trace contributions for any
operator O" — **derivation with a conjectural step**.
L-C3: d = 3 computation (7.1)–(7.6) ⇒ `ñ_s = 4ñ_f`, `ñ_odd = 0` — **derivation**.
L-C4: setting `ρ ∼ 1 − 1/Δ_gap`, with the leading allowed term `−c_0 ln(1 − ρ)` in d = 3 (footnote 14) ⇒
(2.28) `t_n/c_0 ≲ ln Δ_gap/Δ_gap^n` ⇒ (7.7) — **estimate**.
The d = 4 original (1610.09378 (5.2)–(5.5)) states the same assumption: "There are also large numerical
factors in this expression, which we must assume are cancelled by small OPE coefficients."

---

## 3. Strategy S-B in d = 3: what is written, and the explicit structure

### 3.1 What is written (general d, d = 4, flat space)

* General d, components coupled to scalars: (5.134) above, with (5.132)
  `q(Δ) = 2^{1−3d} π^{3+d/2} (d−1)² Γ(d+1)² Γ(Δ) Γ(Δ+1−d/2) / [(d−2)² Γ(d−Δ/2)² Γ(2+Δ/2)⁴ Γ((d+Δ)/2)²] ≥ 0`
  ("is a non-negative function"). In d = 4 (1.6):
  `(t_4 + 2t_2)² = Σ_φ |λ_TTφ|² 15·2⁴π⁴ Γ(Δ_φ−1)Γ(Δ_φ) / [C_T Γ(4−Δ_φ/2)² Γ(2+Δ_φ/2)⁶] + non-scalar`,
  `(t_4 + 2t_2)² = −Σ_φ |λ_TTφ|² 360²π⁴ Γ(Δ_φ−1)Γ(Δ_φ) / [7 C_T Γ(4−Δ_φ/2)² Γ(2+Δ_φ/2)⁶] + non-scalar`.
  The scalar kernel is sign-definite, but the two equations have opposite relative signs, and the "non-scalar"
  kernels (massive spin 2 and higher) are not computed: "We see explicitly that there can be no cancellation
  between massive scalars." and "It is less obvious whether a spin-J operator can be completely canceled by
  spin-J operators … It might even be true that no finite set of spins is sufficient. We leave the
  investigation of this question to future work." "Other sum rules give expressions for other combinations of
  t_2 and t_4, but involve exclusively non-scalars" — not written out.
* Flat space, general D (KKSZ (2.60), (2.63), (A.13)): the shock commutator from the graviton is
  `(1/2)α_2² …  − α_2α_4 … + (α_4²(q_1·q_2)² − α_2²/(D−2)) [(e_3·q_1)²(e_4·q_2)² − (e_3·q_2)²(e_4·q_1)²]`, "this
  vanishes if and only if α_2 = α_4 = 0"; a scalar with `φR²` coupling contributes the same structure as `α_4²`
  with a positive coefficient; "massive spin-2 particles cannot cancel the contributions of higher-derivative
  gravitational interactions".
* Flat space, D = 4, all spins (2201.06602): the Regge-spin-3 sum rule is, from (2.31b), (2.35a), (2.28),

  `−2πG |ĝ_3|² p⁴ = ⟨ ( |c⁺⁺_{J,m²}|² d̃^J_{0,0}(x) − |c⁺⁻_{J,m²}|² d̃^J_{4,4}(x) ) / m⁶ ⟩`, `x = 1 − 2p²/m²`,

  (tree level, no light matter; `⟨·⟩` the positive heavy average (2.34); `c⁺⁻` exists for `J ≥ 4`, `c⁺⁺` for
  even `J`). This is the closest written analogue of the requested sum rule: quadratic form `|ĝ_3|²` on one side,
  a **difference of two sign-definite sums** on the other.
  [D] Expanding formally in `p²` with `d̃^J_{0,0} = ₂F₁(−J, J+1; 1; p²/m²)`, `d̃^J_{4,4} = ₂F₁(4−J, J+5; 1; p²/m²)`
  (their (A.3); coefficients checked with sympy):

  `p⁰`: `⟨(|c⁺⁺|² − |c⁺⁻|²)/m⁶⟩ = 0`; `p²`: `⟨(J(J+1)|c⁺⁺|² − (J−4)(J+5)|c⁺⁻|²)/m⁸⟩ = 0`;
  `p⁴`: `2πG|ĝ_3|² + ⟨(J−1)J(J+1)(J+2) |c⁺⁺|²/(4m¹⁰)⟩ = ⟨(J−5)(J−4)(J+5)(J+6) |c⁺⁻|²/(4m¹⁰)⟩`.

  So in flat space a non-zero `ĝ_3` requires helicity-violating couplings to states of **spin ≥ 6**, with a
  manifestly non-negative kernel. (Formal: forward expansions of gravitational sum rules are delicate,
  2201.06602 §2.5.)

### 3.2 Is the d = 3 CFT sum rule written?

No. Evidence: (i) KKSZ work in general d with `SO(d−1)` irreps `{j_1, j_2|λ|λ_3, λ_4}`; their explicit
equations have poles at d = 3 (see L-B6) because for `SO(2)` the two-row diagrams `(3,1)`, `(2,2)` do not
exist and `2 ⊗ 2` contains no spin-2 irrep, and they do not include the parity-odd ⟨TTT⟩ structure;
(ii) 1904.05892 is d = 4; (iii) 2512.09986 (Dec 2025) computes d = 3 energy-correlator blocks, parity-odd
included, but only "in the background of a scalar source", and notes "For a scalar external source, the
contribution of each conformal block to this commutator vanishes trivially"; (iv) 2608.19322 (Aug 2026)
repeats this as an outlook; (v) none of the 161 papers citing KKSZ has a title indicating it.

### 3.3 [D] The d = 3 structure

*Setting.* Unitary CFT3 with unique stress tensor. Fix timelike `q = (E, 0, 0)`. The states `T(q)|0⟩` span the
two helicities `h = ±2` of the little group SO(2) (conservation and tracelessness remove the rest). For a
primary `O` of dimension `Δ` above the unitarity bound and spin `J`, `O(q)|0⟩` decomposes into SO(2) charges
`m = −J, …, J`, each once, with positive norms `N_{O,m}`. The detector `E(θ)` commutes with translations, and
`E(θ) = U(R_θ) E(0) U(R_θ)†`, hence

`⟨O, m; q| E(θ) |T, h; q′⟩ = (2π)³ δ³(q − q′) e^{i(h−m)θ} a^O_{m,h}`.

Normalise `â^O_{m,h} := a^O_{m,h} (N_{T,h}/N_{O,m})^{1/2} / a^T_{h,h}`. From §1.2: `â^T_{2,2} = 1`,
`â^T_{−2,2} = c` (up to a convention-dependent phase), `|c|² = (t_4² + α_T²)/16`.

*Exact constraints from Poincaré charges.* `∫dθ E(θ) = H` and `∫dθ e^{∓iθ}E(θ) = P^x ∓ iP^y` give
`â^O_{2,2} = δ_{O,T}` and `â^O_{1,2} = â^O_{3,2} = 0` for every `O`.

*Sum rules.* Inserting a complete set of states,
`⟨T,h′|E(θ_1)E(θ_2)|T,h⟩ ∝ Σ_{O,m} conj(â^O_{m,h′}) â^O_{m,h} e^{−i(h′−m)θ_1} e^{i(h−m)θ_2}`. Put

`Â_m := Σ_O |â^O_{m,2}|² ≥ 0`, `B̂_m := Σ_O conj(â^O_{m,−2}) â^O_{m,2}`.

If `[E(θ_1), E(θ_2)] = 0` on these matrix elements with no contact term at `θ_1 = θ_2`, and the sums converge
mode by mode (hypotheses of L-B1, L-B2), then for all integers `m`

**(D)** `Â_m = Â_{4−m}`, **(O)** `B̂_m = B̂_{−m}`.

The stress tensor contributes `Â_2 ∋ 1`, `Â_{−2} ∋ |c|²` and `B̂_{±2} ∋ c` (so it drops out of (O)). The first
non-trivial relations are `Â_0 = Â_4`, `Â_{−1} = Â_5`, and

**(★)** `(t_4² + α_T²)/16 + Σ_{O≠T, J_O ≥ 2} |â^O_{−2,2}|² = Σ_{O: J_O ≥ 6} |â^O_{6,2}|²`.

Equivalently, the T-block of the commutator is `[M(θ_1), M(θ_2)] = 2i|c|² sin(4θ_12) diag(−1, +1)` with
`M(θ) = [[1, c̄e^{−4iθ}],[c e^{4iθ}, 1]]` (sympy, Appendix A.3) — the d = 3 counterpart of 1904.05892 (4.23).

*Consequences (conditional on the stated hypotheses; planar limit: sums over single-trace operators by L-B4).*

1. `t_4² + α_T² ≤ 16 Σ_{J_O ≥ 6} |â^O_{6,2}|²`: a non-Einstein ⟨TTT⟩ requires operators of spin ≥ 6 in `T × T`
   whose `m = 6` amplitude is non-zero. Operators of spin 4 and light operators of spin ≤ 2 cannot help: they
   enter (★) only on the left, with the same sign as `|c|²`. No assumption about light scalars is needed
   (scalars have only `m = 0` and do not enter (★)).
2. `Â_0 = Â_4`: any ⟨TTφ⟩ coupling to a scalar (and any `m = 0` amplitude) must be balanced by `m = 4`
   amplitudes of operators of spin ≥ 4 — the sum-rule form of (C3).
3. If the product `E(θ_1)E(θ_2)` is a positive operator-valued measure (commuting positive operators), then
   `{Â_{2−k}}_k` are the Fourier coefficients of a positive measure, so `0 ≤ Â_m ≤ Â_2 = 1`; with (★) this
   gives back the collider disc `t_4² + α_T² ≤ 16` and nothing stronger.

*Checks.* (i) Free scalar: `T(q)|0⟩` is a back-to-back two-particle state with wavefunction `∝ e^{ihθ}`;
`⟨h|E(θ_1)E(θ_2)|h⟩ ∝ δ(θ_12) + δ(θ_12 − π)`, so `Â_m = 1` for every even `m` and 0 for odd `m`; `|c| = 1`
(`t_4 = 4`); the bilinear primaries are the conserved currents `J_s` with helicities `±s` only; (★) reads
`1 + 0 = |â^{J_6}_{6,2}|² = 1`. ✓ (This example also shows that the mode sums converge only distributionally,
and that a symmetric contact term `δ(θ_12)` is present.) (ii) Flat space: the same "spin ≥ 6" threshold
appears in the `p⁴` relation of §3.1. (iii) The d → 3 degeneration of KKSZ (5.134): the stress-tensor
coefficient of the second equation carries `(d² − 9)`, consistent with T and scalars decoupling into
different relations ((★) and `Â_0 = Â_4`) at d = 3.

*What is not done.* The kernels: `â^O_{m,2} = Σ_a λ^{(a)}_{TTO} K^{(a)}_m(Δ, J)` with `a` running over the
conserved ⟨TTO⟩ structures. Until the `K^{(a)}_m` are computed, (★) is a structural identity, not a usable
inequality on OPE data.

### 3.4 [D] Why superconvergence alone gives no rate, and what would

Every relation (D), (O) is homogeneous: the Einstein part of the T-block is the identity matrix in helicity
space, contributes only to `Â_2`, and cancels identically. A linear combination `Σ_k w_k (Â_{2−k} − Â_{2+k}) = 0`
has the form `w_4|c|² + Σ_heavy (quadratic form in λ_TTO) = 0`. If the heavy form were non-negative for all
`Δ ≥ Δ_gap` with `w_4 > 0`, it would force `c = 0` exactly at finite gap, which is not expected (finite-gap
string backgrounds are expected to have non-zero higher-derivative cubic couplings [M]). Hence a decay rate
requires an *a priori* bound on the heavy side. Available inputs:

* `Â_m ≤ 1` (item 3 above) — gives only the collider disc, unless the kernels have a hierarchy
  `|K_6|² ≪ |K_m|²` at large `Δ` uniformly in `J` (unknown until L-B5 is done; no reason to expect it in the
  Regge regime `J ∼ Δ`).
* KKSZ §6.1: Regge-pole dominance with "no outliers" — an assumption.
* The Regge-spin-2 sum rule (flat space: `−B_2^{(1)}|_low = 16πG/p² + 2πG|ĝ_3|²p⁶`), in which the Einstein part
  does *not* cancel. In CFT this is the `k = 2` family of L-A2; its T-block is not a product of two one-point
  event shapes, and its convergence is marginal in the planar expansion (`J_0^{planar} = 2`), so it must be
  used at finite `c_T` as in 2106.10274. This is where S-B merges into S-A: in 2201.06602 the bound on
  `|ĝ_3|²` is obtained from `B_2` and `B_3` together (functional (3.7)).

**Answer to the question of item 3.** There is no sum rule in the literature of the form "quadratic form in
non-Einstein couplings = Σ |λ_TTO|² × (kernel of one sign)". What exists is "quadratic form + (sign-definite
sum) = (sign-definite sum)": in general d only the scalar kernel is explicit (KKSZ (5.132), (5.134)); in 4D
flat space all kernels are explicit (2201.06602 (2.35a)); in d = 3 CFT nothing is written, and (★) is the
form it takes. Each such relation gives directly "non-Einstein ⇒ higher-spin operators" as an inequality
(with spin ≥ 6 on the right in d = 3 and in 4D flat space), under: planarity or exact commutativity (Regge
intercept `< 3`), convergence of the mode-projected sums, absence of contact terms. No assumption on light
scalars is needed.

---

## 4. Strategy S-A: harmonic analysis needed in the AdS-curvature regime

### 4.1 What the sources give

* Scalars (2106.10274): Regge-regime action (2.59) with the `H^{d−1}` harmonic function
  `P_{(2−d)/2+iν}(cosh β)`; reduction (5.41)–(5.43) to `Q_n(β) = ∫_0^∞ dν ν^n P_{(2−d)/2+iν}(cosh β)`; AdS4
  (5.63) with `f^{d=3}(p) → 1310p − 1080p²` (5.56). "(Generically Q_n(β) ∼ 1/β^{n+1} at small β, but Q_1(β)
  approaches a constant as β → 0.)"
* Gravitons in flat space (2201.06602): heavy densities (2.33), (2.35) involve `d̃^J_{0,0}`, `d̃^J_{4,4}`
  (diagonal, coefficients `|c⁺⁺|²`, `|c⁺⁻|²`), `d̃^J_{4,0}` (coefficient `c⁺⁺c⁺⁻`) and `d̃^J_{4,−4}`; large-J limit
  (A.11): `lim d̃^J_{h,h′}(1 − 2p²/m²) = J_{h−h′}(bp)/(p/m)^{h−h′}`, `b = 2J/m`; scaling limit (3.5):
  `lim_{m→∞} F[m², bm/2] = (4/m⁴) ∫_0^1 dp ψ(p) (|c⁺⁺|² + |c⁺⁻|²) J_0(pb) + …`, "dominate[d] by the B_2 component",
  "In this limit, F⁺⁺ and F⁺⁻ coincide".
* Gravitons in AdS: nothing. 2311.04271 give no Regge-regime formula.

### 4.2 Inference (not in any source)

The flat transverse plane `R²` with Bessel kernels `J_μ(bp) e^{iμφ}` is replaced in AdS4 by `H²`, whose
harmonic functions of SO(2) weight `μ` are the associated conical functions. The expected substitution,
generalising `J̃(βp) ↦ P_{(2−d)/2+iν}(cosh β)` of 2106.10274 §5.3.2, is

`J_μ(bp)/p^μ  ↦  P^{−μ}_{−1/2+iν}(cosh β)`, `μ = |h − h′| ∈ {0, 4, 8}`,

consistent with the flat limit `P^{−μ}_{−1/2+iν}(cosh β) → ν^{−μ} J_μ(νβ)` ([N], Appendix A.1). Equivalently,
`P^{−μ}` is obtained from `P_{−1/2+iν}` by `μ` applications of the lowering operator, i.e. `β`-derivatives:
these are the "spinning harmonics on H²". The needed positivity statements would then be sign conditions on

`Q_n^{(μ)}(β) := ∫_0^∞ dν ν^n P^{−μ}_{−1/2+iν}(cosh β)`, `μ ∈ {0, 4, 8}`,

for the exponents `n` set by the small-`p` expansion of the weights, together with positive-semidefiniteness
of the 2 × 2 (parity-even) or 3 × 3 forms in `(c⁺⁺, c⁺⁻)` once off-diagonal sum rules are used.

For the bound on `|ĝ_3|²` only the MHV-type sum rules `B^{(1)}_k` are needed (2201.06602: "These use only
B^{(1)}"), whose heavy densities are diagonal; in the Regge regime the leading term is the `B_2` weight against
the `μ = 0` kernel. For (3.3a) the `B_2` weight is `(1−p)³p(65 + 155p + 47p²) = 65p − 40p² + …`, i.e.
`a_1 > 0`, `a_2 < 0`, the same pattern as (5.56). **So, granted the missing spinning Regge formula (L-A5), the
leading curved-regime positivity for the graviton functional reduces to exactly `Q_1 > 0` and `Q_2 < 0`, which
`POSITIVITY_LEMMA.md` proves.** Differences between the ++ and +− channels, and the `B_3`, `B_4` components,
appear at relative order `1/m` and are not controlled by any existing formula.

### 4.3 [D] The `μ = 4, 8` moments by the method of the lemma

Mehler–Dirichlet formula ([M]: DLMF §14.12; [N] verified to 20 digits for `μ = 0, 4, 8`, Appendix A.1): for
`μ > −1/2`,

`P^{−μ}_{−1/2+iν}(cosh β) = (2/π)^{1/2} (sinh β)^{−μ} Γ(μ+1/2)^{−1} ∫_0^β cos(νt) g_μ(t) dt`, `g_μ(t) := (cosh β − cosh t)^{μ−1/2}`.

For `μ > 1/2`, `g_μ` is positive, decreasing on `(0, β)`, vanishes at `t = β`, and
`|P^{−μ}_{−1/2+iν}| = O(ν^{−μ−1/2})`, so the moments converge absolutely for `n < μ − 1/2`. Fourier inversion
(even `n`) and the finite-part identity `∫_0^∞ ν cos(νt) dν = −1/t²` (odd `n = 1`; justified by Abel
regularisation as in Lemma Q) give, with `C_μ(β) := (2/π)^{1/2} (sinh β)^{−μ}/Γ(μ+1/2)`:

* `Q_0^{(μ)} = (π/2) C_μ (cosh β − 1)^{μ−1/2} > 0`;
* `Q_1^{(μ)} = C_μ [ ∫_0^β (g_μ(0) − g_μ(t)) t^{−2} dt + g_μ(0)/β ] > 0`;
* `Q_2^{(μ)} = (π/2) C_μ (μ − 1/2)(cosh β − 1)^{μ−3/2} > 0`.

[N] These agree with direct numerical integration for `μ ∈ {4, 8}`, `β ∈ {0.5, 2}` (Appendix A.2). Note the sign
of `Q_2` is opposite to the `μ = 0` case (`Q_2 = −1/(8 sinh³(β/2)) < 0`), as in flat space
(`∫_0^∞ p² J_μ(bp) p^{−μ} dp ∝ Γ(3/2)/Γ(μ − 1/2)`). Status: proof sketch + numerics; relevance conditional on §4.2.

### 4.4 What is not reducible to such moments

(i) The regime `m ∼ M`, `J ≫ 1` (scalar case: App. F of 2106.10274) and the threshold region `τ ≈ Δ_gap`,
where the functional is evaluated on blocks, not on harmonic functions. (ii) Finite `M = Δ_gap`: the moments
are the leading term of an expansion in `1/M`; oscillatory endpoint terms are only shown to be subleading.
(iii) The light side (graviton pole with AdS propagator, `|c|²` term).

---

## 5. Ranked sub-theorems (unproved, on the route, pen-and-paper + modest computation)

### ST1. The d = 3 ANEC-commutator sum rules with explicit kernels

* **Statement.** In the setting of §3.3, compute for each conserved ⟨TTO⟩ structure `a` (1⁺ + 1⁻ for `J = 0, 2`;
  2⁺ + 1⁻ for even `J ≥ 4`; 1⁻ for odd `J ≥ 5`) the functions `K^{(a)}_m(Δ, J)`, `|m| ≤ J`, such that
  `â^O_{m,2} = c_T^{−1/2} Σ_a λ^{(a)}_{TTO} K^{(a)}_m(Δ, J)` for unit-normalised `O`, and prove:
  (a) under the hypotheses of L-B1, L-B2, for every integer `m`: `Σ_O |â^O_{m,2}|² = Σ_O |â^O_{4−m,2}|²`, in
  particular (★) `(t_4² + α_T²)/16 + Σ_{O≠T} |â^O_{−2,2}|² = Σ_{J_O≥6} |â^O_{6,2}|²`;
  (b) `K^{(a)}_m(Δ, J)` has double zeros in the squared amplitude at the double-trace dimensions
  `Δ = 6 + J′ + 2n` appropriate to `[TT]` (planar reduction to single-trace operators);
  (c) for the free scalar and free Majorana fermion, `|â^{J_s}_{±s,2}|² = 1` for all even `s ≥ 4`;
  (d) decide whether `sup_{Δ ≥ Δ_gap, J ≥ 6, λ} |â_{6,2}|² / Σ_{m ∉ {1,2,3}} w_m |â_{m,2}|²` tends to 0 as
  `Δ_gap → ∞` for some fixed positive weights `w_m` (if yes, a gap-dependent bound follows from `Â_m ≤ 1`).
* **Obstruction removed.** L-B5, L-B6 in d = 3; the d = 3 analogue of KKSZ (1.6) with all spins, in a form
  where every term is a square; the first explicit inequality "non-Einstein ⟨TTT⟩ ⇒ spin ≥ 6 single-trace
  operators" for CFT3, parity-odd structure included.
* **Inputs.** d = 3 ⟨TTO⟩ structures (1708.05718 App. A; helicity basis 2102.08160 [A]); light transform and
  Fourier transform of three-point structures (1904.05905 §5.1.3–5.1.4; 2512.09986 §3 for the scalar-source
  d = 3 case, results in `₅F₄`); free-field ⟨TTJ_s⟩ for the checks.
* **Risk.** Computational (conventions, parity-odd structures): moderate. Conceptual: the result is
  conditional on commutativity without contact terms and on mode-wise convergence, as in KKSZ; (d) is
  probably negative.

### ST2. Uniform (non-asymptotic) curved-regime positivity in AdS4 — "Lemma Q⁺"

* **Statement.** For `f(p) = p(1−p)²(1310 + 1540p − 10318p² − 4697p³ + 15680p⁴)` (2106.10274 (5.55)) and for
  `f(p) = (1−p)³ p (65 + 155p + 47p²)` (2201.06602 (3.3a), weight of `B_2^{(1)}`), define
  `G_M[f](β) := ∫_0^1 dp f(p) P_{−1/2+iMp}(cosh β)`. Then there is an explicit `M_0` (target: `M_0 ≤ 5`) such
  that `G_M[f](β) > 0` for all `β > 0` and all `M ≥ M_0`, with
  `M² G_M[f](β) = a_1 Q_1(β) + a_2 Q_2(β)/M + R_M(β)`, `|R_M(β)| ≤ C M^{−2} (a_1 Q_1(β) + |a_2 Q_2(β)|/M)` for
  `β ≥ β_0/M`, and `G_M[f](β) ≥ (1 − ε) f̂(Mβ)` for `β ≤ β_0/M`, where `f̂(b) = ∫_0^1 f(p) J_0(bp) dp > 0`.
* **Obstruction removed.** L-A8 beyond the two leading moments: it replaces "rigorously positive at large M"
  and "dropping the oscillatory terms" by a statement valid for every `M ≥ M_0`, for the exact leading Regge
  kernel `2G_M/m⁴` of (2.59)/(5.61). It is needed to make the scalar AdS4 bound (5.64) a theorem and, granted
  L-A5, is the `B_2` part of the graviton functional.
* **Inputs.** Mehler's formula used in `POSITIVITY_LEMMA.md`; integration by parts at the endpoint `p = 1`
  (the factors `(1−p)²`, `(1−p)³`); interval arithmetic on a compact `(M, β)` window.
* **Evidence.** [N] `G_M > 0` on the grid `M ∈ {5, 20, 80}`, `β ∈ [0.02, 10]` for both `f`; `f̂(b) > 0` for
  `b ≤ 40` (Appendix A.2).
* **Risk.** Low–moderate. The lemma concerns the leading Regge kernel only; corrections `O(ν²/m², 1/m²)` in
  (2.59) are outside it.

### ST3. A proof of positivity of the 4D graviton functional on all heavy states

* **Statement.** With `ψ_2(p) = (1−p)³p(65 + 155p + 47p²)`, `ψ_3(p) = 2(1−p)³(81 − 1927p + 4601p²)` and `M = 1`,
  for all real `m ≥ 1`:
  `F⁺⁺_J(m) := m⁻⁶ ∫_0^1 dp [ψ_2(p)(2m² − p²) + ψ_3(p)] P_J(1 − 2p²/m²) + 5(2J(J+1) − 1) m⁻¹⁰ ≥ 0` for even `J ≥ 0`,
  `F⁺⁻_J(m) := m⁻⁶ ∫_0^1 dp [ψ_2(p)(2m² − p²) − ψ_3(p)] P^{(0,8)}_{J−4}(1 − 2p²/m²) + 5(2(J−4)(J+5) − 1) m⁻¹⁰ ≥ 0`
  for all integers `J ≥ 4`.
  (Heavy action of (3.3a) from (2.33), (2.35a); `d̃^J_{4,4}(x) = P^{(0,8)}_{J−4}(x)` [M, Jacobi form of (A.3)].)
* **Obstruction removed.** L-A7: the published evidence is plots for `J ≤ 300`, `m ≤ 16` and a scaling limit.
* **Inputs.** Uniform (Hilb-type) Bessel asymptotics of Jacobi polynomials; exact rational arithmetic for
  low `J`; interval arithmetic for the compact window.
* **Evidence.** [N] positive on a grid of 11 values of `m ∈ [1, 16]`, `J ≤ 100` (minimum `1.2·10⁻⁶`); negative
  without the `B_4` term.
* **Risk.** Moderate (two parameters, threshold region `m → 1`, `x → −1`). Two caveats found in this session:
  (i) by my reading of (2.32a) the low-energy side of (3.3a) contains `g_5` through `∂_{p²}B_4`; the printed
  bound (3.4a) is reproduced (`37.84 log − 45.45`) only when that term is dropped — I could not determine from
  the text how `g_5` is disposed of; (ii) the pure `B_2, B_3` functional (3.7) reproduces its printed bound
  (`1596.7 log − 3208` vs `1598 log − 3211`) but is not positive as printed near threshold. An exact-rational
  functional of type (3.7) would have to be re-derived (a small linear programme, not a large SDP).

### ST4. Spinning ν-moments on H²

* **Statement.** For `μ ∈ {4, 8}`, all `β > 0` and all real `0 ≤ n < μ − 1/2`: `Q_n^{(μ)}(β) > 0`, with the closed
  forms of §4.3 for `n = 0, 1, 2`.
* **Obstruction removed.** The off-diagonal part of L-A8 (needed for bounds beyond `|ĝ_3|`, e.g. `ĝ_4`), if
  §4.2 is right.
* **Inputs.** Mehler–Dirichlet formula; finite-part integrals. **Risk.** Low mathematically; relevance
  conditional on L-A5.

### ST5. The spinning Regge-moment formula in d = 3

* **Statement.** For the `k = 2` and `k = 3` stress-tensor sum rules of 2311.04271 in d = 3, and a heavy exchanged
  primary `(Δ, J)` with couplings `λ = (λ⁺⁺, λ⁺⁻, …)`: as `Δ, J → ∞` at fixed `β` and `ν`,
  `Ψ_{k,ν}[Δ, J; λ] = (2/m^{2k}) [ P_{−1/2+iν}(cosh β) (|λ⁺⁺|² + (−1)^k |λ⁺⁻|²) + (off-diagonal terms with
  P^{−4}, P^{−8}) ] (1 + O(ν²/m², 1/m²))`, in a basis of structures specified by the flat limit.
* **Obstruction removed.** L-A5 (and with it the applicability of Lemma Q, ST2, ST4 to gravitons).
* **Inputs.** Conformal Regge theory for spinning blocks; the helicity basis of 2102.08160; the "wavefunction
  trick" of 2106.10274 §3.1. **Risk.** High: this is the step 2311.04271 defer; it may need the d = 3
  spinning blocks in a regime where no closed form is available.

---

## 6. What cannot be done without large-scale numerics or new ideas

1. **A quantitative theorem `|c| ≤ K (log Δ_gap)^a/Δ_gap⁴`.** It needs L-A2/L-A3 in d = 3, L-A5, L-A6, L-A9
   together, plus positivity of one functional in *all* regimes at finite `Δ_gap`. None of ST1–ST4 delivers
   it; ST5 is one of three missing derivations. Without L-A5 there is no route by pen and paper.
2. **A rate from superconvergence alone** (§3.4): requires a new idea for bounding the heavy side (a
   replacement for "no outliers").
3. **Optimal constants.** `24.9` in (4.4) comes from numerical optimisation over improved sum rules; the AdS4
   optimum needs semidefinite programming with d = 3 spinning blocks, or at least the 4D graviton partial
   waves with many null constraints.
4. **Contact terms of ⟨TTTT⟩ (`g_4`, `ĝ_4`, …).** Needs the sum rules not obtainable from shock commutators
   (L-A3) and the full matrix positivity.
5. **Theorem-level foundations.** Commutativity `[E, E] = 0` on the relevant states without loopholes;
   Wightman-type theorems for spinning four-point functions; the chaos bound in the planar limit. These are
   long technical programmes.
6. **Uniformity in `c_T`.** Superconvergence fails order by order beyond the planar level ("at one-loop in the
   bulk … we expect J_0^{1-loop} = 3"; "it may be necessary to use eikonal techniques to re-sum gravitational
   exchanges", KKSZ §4.5, §6.2.1); a statement for a sequence of finite-`c_T` CFTs needs a definition of the
   single-trace gap at finite `c_T` and error terms in `1/c_T` relative to powers of `1/Δ_gap` (hierarchies
   such as `√c_T ≫ Δ_gap` in 1805.07393 (2.25)).
7. **The logarithm.** Whether the true exponent is `a = 1/2` (uplift of (4.4)), `a = 1` (C1), or `a = 0` cannot
   be decided by present methods: "If the flat-space divergences could be resolved using a purely flat-space
   mechanism, then it should be possible to remove the log term" (2106.10274 §5.5). Recent IR-finite
   proposals (2512.13780, 2610.19149, 2501.17949 [A]) concern scalar/EFT amplitudes in flat or half-flat
   space and have not been applied to `ĝ_3`.

---

## Appendix A. Numerical and symbolic checks run in this session (mpmath / sympy)

All are sanity checks. Working precision 15–30 digits.

**A.1 Mehler–Dirichlet formula and flat limit.**

```python
from mpmath import mp, mpf, mpc, legenp, quad, cos, cosh, sinh, sqrt, pi, gamma, besselj
mp.dps = 25
P  = lambda mu, nu, b: legenp(mpc(-0.5, nu), -mu, cosh(b), type=3)        # P^{-mu}_{-1/2+i nu}(cosh b)
MD = lambda mu, nu, b: sqrt(2/pi)*sinh(b)**(-mu)/gamma(mu+mpf(1)/2) * \
     quad(lambda t: cos(nu*t)*(cosh(b)-cosh(t))**(mu-mpf(1)/2), [0, b])
# mu in {0,4,8}, (nu,b) in {(0.7,0.5),(3,1.3),(6,2.5)}: agreement to >= 14 digits (20 digits for mu = 4, 8)
# flat limit nu = 400, b = 0.02: P^{-mu} vs nu^{-mu} J_mu(nu b): 0.171646/0.171651, -4.1163e-12/-4.1155e-12, 3.4085e-22/3.4097e-22
```

**A.2 Moments.** `Q_n^{(μ)}` by direct integration over `ν ∈ [0, 80]` against the closed forms of §4.3:
`μ = 4, β = 0.5`: `Q_0` 0.00108527/0.00108526, `Q_1` 0.00470561/0.00470512, `Q_2` 0.0298006/0.0297620;
`μ = 8, β = 2`: `Q_0` 6.08173e-6 (both), `Q_1` 8.24066e-6 (both), `Q_2` 1.651332e-5 (both). Differences are the
truncated tails. `G_M[f](β)` of ST2 via `P_{−1/2+iν}(cosh β) = (√2/π)∫_0^β cos(νt)(cosh β − cosh t)^{−1/2}dt`:
for `f` = (5.55): `M = 5`: 37.13 (β = 0.02) … 0.0191 (β = 10); `M = 20`: 37.13 … 0.00118; `M = 80`: 37.08 …
7.6e-5; for `f` = B_2 weight of (3.3a): `M = 5`: 6.17 … 9.4e-4; `M = 20`: 6.12 … 5.9e-5; `M = 80`: 5.49 … 3.7e-6;
all ten `β` values positive in each case. `min_b f̂(b)` on `b ∈ [0, 40]`: 0.0020 and 0.00064.

**A.3 Algebra.**

```python
import sympy as sp
t4, aT = sp.symbols('t4 alphaT', real=True)
Mlin = sp.Matrix([[1 - t4/4, aT/4], [aT/4, 1 + t4/4]])                    # 1707.03007 (1.5), basis (xy, xx-yy)
plus, minus = sp.Matrix([sp.I, 1])/sp.sqrt(2), sp.Matrix([-sp.I, 1])/sp.sqrt(2)
# <+|M|+> = 1 ; <-|M|+> = (t4 + i alphaT)/4 ; eigenvalues 1 +- sqrt(t4^2 + alphaT^2)/4
c, th1, th2 = sp.symbols('c theta1 theta2'); cb = sp.conjugate(c)
M = lambda t: sp.Matrix([[1, cb*sp.exp(-4*sp.I*t)], [c*sp.exp(4*sp.I*t), 1]])
# M(th1) M(th2) - M(th2) M(th1) = 2 i |c|^2 sin(4(th1 - th2)) diag(-1, 1)
# 2F1(-J, J+1; 1; y)  = 1 - J(J+1) y + J(J-1)(J+1)(J+2) y^2/4 + ...
# 2F1(4-J, J+5; 1; y) = 1 - (J-4)(J+5) y + (J-5)(J-4)(J+5)(J+6) y^2/4 + ...
```

**A.4 Transcription checks of 2201.06602.** Low-energy side, units `πG = 1`, `M = 1`, `L = log(1/m_IR)`:
(3.7) with `ψ_2 = (1−p)²(−154.020p⁷ + 2006.03p⁶ − 3443.55p⁵ − 2322.66p⁴ + 5294.87p³ − 2194.19p² + 798.968p)` on
`B_2` and `ψ_3 = 3701.82p⁶ − 18251.2p⁵ + 32684.9p⁴ − 25581.3p³ + 7603.72p² − 157.911p + 0.00171` on `B_3`:
`|ĝ_3|² ≤ 1596.7 L − 3208.0` (printed: `1598 L − 3211`). (3.3a) without the `g_5` term: `37.84 L − 45.45`
(printed: `37.8 L − 45.4`). Heavy side on the grid `m ∈ {1, 1.02, 1.05, 1.1, 1.2, 1.5, 2, 3, 5, 8, 16}`,
`J ≤ 30` and `J ∈ {40, 60(, 100)}`: (3.3a) all positive; (3.7) as printed negative in the +− channel at
`m = 1` (`J = 4`: −0.0053; `J = 12`: −9.1; growing with `J`), `m = 1.02` (`J = 12 … 60`: −2.6 … −19) and
`m = 1.05` (≥ −0.74); `ψ_3(1) = 0.0307`.

Fetched texts and scripts are in the session scratchpad (not part of the repository):
`/tmp/claude-0/-home-user-MAJOR/60eea3a7-07d6-5453-af55-c4346d86f9b5/scratchpad/hpps/` and
`/tmp/claude-0/-home-user-MAJOR/60eea3a7-07d6-5453-af55-c4346d86f9b5/scratchpad/hpps_scripts/`.
