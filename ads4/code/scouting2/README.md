# Scouting pass 2: one sharp 2024–2026 conjecture to settle (28 Nov 2026)

Written by a research agent. This file is the only one written in the repository, and nothing was committed. The probe
scripts ran from the session scratchpad. Their full text is in the Appendix, so they can be re-run.

**Source marks**
* **[R]** arXiv LaTeX source read (or grepped) this session.
* **[A]** abstract read this session (arXiv API).
* **[C]** computed here (Appendix).
* **[D]** derived here by hand, not machine-checked beyond what [C] says.
* **[I]** my inference.
* **[M]** memory, not checked.

"Open" means no resolution was found in the INSPIRE citation list (checked by recid, 28 Nov 2026) or in arXiv full-text
searches. That is a search-level statement, not a proof.

**Method of search.**
* I harvested all 2,410 arXiv abstracts (hep-th, cond-mat.str-el, cond-mat.stat-mech, math-ph; Jan 2024 – 28 Nov 2026)
  that contain "conjecture". I filtered them for bounds, inequalities and positivity, and read about 400 by hand.
* I ran targeted web searches and read the sources of the 12 most promising papers.
* Everything ruled out in `../scouting/README.md` was skipped: the 3D entropy bound and Kusuki–Tang, the 6D a-theorem,
  η/s, collider energy correlators, Bueno–Casini `C_T/F` in d = 3, and modular-bootstrap gaps.

---

## 0. Bottom line (blunt)

1. **Recommendation: the Allameh–Shaghoulian sign conjecture `c₁ ≥ 0`** for the curvature coefficient of the
   thermal effective action (arXiv:2402.13337, JHEP 01 (2025) 200).
   * **The test.** Compute the sign of `c₁` in the critical O(N) model at order `N⁰`.
   * **Why that test is decisive.** At `N = ∞`, `c₁` is *exactly zero*: the coefficient is proportional to the
     renormalised thermal `⟨φ²⟩`, which the gap equation sets to zero (§1.4; reproduced numerically here).
   * **Consequence.**
     * The sign of the single number `c₁⁽⁰⁾` decides the conjecture in an interacting unitary CFT. If `c₁⁽⁰⁾ < 0`, the
       conjecture fails for all large enough N.
     * The same holds for CP^{N−1}, whose `c₁` also vanishes at `N = ∞`.
     * The `1/N` correction is explicitly open: Mauro–Vichi 2606.05059 list it as future work, and 2606.30716 (Nov 2026
       v1) says "the c₁ ≥ 0 conjecture … remains open".
   * **Fit to our methods.** It is an exact large-N thermal computation on S¹×S² (and S¹×H² as a check), with certified
     arithmetic for the sign.
   * **Odds.**
     * P(a clean, cross-checked sign within 4 weeks) ≈ 0.6.
     * P(the sign is negative, i.e. a counterexample) ≈ 0.4. I have no real prior.
     * P(PRL | counterexample) ≈ 0.35. P(PRL | the conjecture is confirmed) ≈ 0.05: that outcome is a JHEP paper.
     * **Overall P(PRL) ≈ 0.1.** That is no better than the previous pass's best candidates, but the decision comes
       fast and the downside is cheap.
2. **Flag: adjacency to the excluded topic.** `c₁` is the next term after `f` (that is, `c_S/3`) in the same
   high-temperature expansion on S¹×H² and S¹×S².
   * In the Herzog–Stergiou variables, `c₁/C_T` is the constant term of `Ψ(x)` as `x → ∞` [I].
   * It is **not** a bound on the entropy and not a Cardy-growth bound. But if the user counts it as a "relative", drop
     it and take rank 2 or 3 below.
   * The computation would in any case feed the user's S¹×H² large-N work.
3. **No other candidate found has comparable odds × payoff.**
   * The best alternative with a real PRL payoff is the interface bound `c_cap ≥ c_eff` (Kusuki–Suzuyama–Tang
     2610.02309): high payoff, but a general proof within weeks is unlikely (≈ 5%).
   * The cheapest decisive alternative is the Pelissetto–Vicari exponent conjecture `Δ_ε ≥ 2Δ_φ` (2510.17637). It is
     *saturated at N = ∞* in every Gross–Neveu–Yukawa (GNY) model, so `1/N` scans could break it. But my probe of the
     chiral O(M) GNY family gives the wrong sign for a counterexample (§2.1), so P(counterexample) is about 5–10%.

---

## 1. Recommended candidate: `c₁ ≥ 0` for the thermal effective action

### 1.1 Statement and source [R]

* **Setting.** A CFT_d on `S¹_β × M_{d−1}` at high temperature is described by a local effective action on M:
  `log Z = −∫d^{d−1}x √g (−f T^{d−1} + c₁ T^{d−3} R + …)`. Here f is the free-energy density and `c₁` the first
  curvature coefficient (Benjamin–Lee–Ooguri–Simmons-Duffin 2306.08031).
* **d = 3 on a round S² of radius r** (this normalisation was checked numerically, §1.4):

      log Z(S¹_β × S²_r) = 4π r² f T² − 8π c₁ + O(T⁻²).

* **Conjecture.** Allameh–Shaghoulian 2402.13337, §"Evidence for c₁ ≥ 0", verbatim: *"Conjecture: The coefficient c₁ of
  the leading subextensive correction to the effective action coming from the gapped sector satisfies c₁ ≥ 0."* They
  state it for general d ≥ 3.
  * They write: *"A significant check of this constraint would be to calculate the leading subextensive correction to
    the Wilson-Fisher fixed point in d = 3. Since c₁ = 0 for the d = 3, 4 conformally coupled scalar, perturbative
    methods like the epsilon expansion or the large N expansion would suffice to probe a potential violation."*
  * Their evidence is free fields (fermions and vectors > 0, the 4D scalar = 0) and Einstein holography.
  * In d = 3 the free scalar has `c₁ = 0`. Its `log T` term comes from the gapless dimensionally reduced sector,
    which is why the conjecture is phrased for the gapped sector.
* **Why d = 3 is special** (Mauro–Vichi 2606.05059 §3 [R]). In d = 3, `∫√g R = 4πχ(M)` is topological. So `c₁` multiplies
  the Euler characteristic and cannot be extracted from flat-space correlators or perturbations around flat space. It
  is a genuinely global, "topological" thermal coefficient, fixed only by putting the theory on S² or on a higher-genus
  surface.

### 1.2 Status [R/A]

* INSPIRE lists 36 citing papers. None claims a proof or a counterexample.
* **2606.30716** ("CFTs on squashed spheres and the thermal effective action", v1, Nov 2026), conclusions, verbatim:
  *"All examples here and in previous work are … compatible with c₁≥0 and c₂≤0. This includes free fermions, the d=3
  Ising model, the large-N critical O(N) model, Einstein-gravity holographic CFTs … It would be valuable either to find a
  counterexample or to identify a general principle underlying these signs … In d=4, a proper counterexample to c₂ ≤ 0
  does exist: free Maxwell theory has c₂=1/144 >0 … The c₁≥0 conjecture … remains open."*
* **Large-N O(N): `c₁ = 0` at leading order.**
  * Mauro–Vichi 2606.05059, eq. (1.1): `f = 2ζ(3)N/5π`, `c₁ = 0`, `c₂ = −√5 log φ N/96π`.
  * 2606.30716 lens computation: "Note the absence of a β⁰ term".
  * Mauro–Vichi conclusions, verbatim: *"it would be interesting to study 1/N corrections, where the dependence of the
    saddle-point field on the background geometry can no longer be ignored"*. **So the decisive order is untouched.**
* **Large-N Gross–Neveu: `c₁ > 0` at leading order.** The thermal saddle has `m_* = 0`, so `c₁` = N × (free Dirac
  value) (2606.30716 [R]). Hence it is not decisive.
* **3D Ising: the sign of `c₁` is unknown.** 2306.08031, App. 3d Ising, verbatim: *"we are not yet able to reliably
  estimate the sign of c₁"*. So the "compatible" in 2606.30716 means only "not excluded".

### 1.3 Why a resolution would be broadly interesting

* **What `c₁` is.** It is the only universal O(T⁰) datum of a 3D CFT at high temperature. It multiplies χ(M), so it is
  a "thermal Euler term".
* **Links to other quantities** [I]:
  * the constant term of `log Z` on S¹×S²;
  * via S¹×H², the n-independent part of the universal disk Rényi term as n → 0;
  * the constant term of the Herzog–Stergiou `Ψ(x)` at large x.
* **What either outcome gives.**
  * **A counterexample in the O(N) model** would be clean: the canonical interacting unitary CFT, every N above some
    N₀. It would show that the thermal effective action has no sign rule at curvature order. In the language of
    Allameh–Shaghoulian's induced-gravity remark, it would also show that the sign of the induced Einstein–Hilbert
    term on the spatial slice is not universal.
  * **A positive sign** would be the first interacting-CFT check and would motivate a proof. It would not be a PRL.

### 1.4 Feasibility probe [C, D]

All four checks are in `probe/gap_s2.py` and `probe/c1_ninf.py` (Appendix), with β = 1, 30 digits.

| check | result |
|---|---|
| flat thermal mass `m₀ = 2 log φ` | reproduced (0.962423650…) |
| curvature shift of the N = ∞ saddle on S² | `(m² − m₀²) r² = 0.083445, 0.083361, 0.083340, 0.083335` at r = 10, 20, 40, 80 → **1/12** [C]; analytic `δm² = 1/(12r²)` from Euler–Maclaurin [D] |
| `log Z/N − (2ζ(3)/5π)·4πr²` at the saddle | `6.47e-5, 1.61e-5, 4.03e-6` at r = 10, 20, 40: **constant term = 0** [C], residual `∝ β²/r²` |
| the residual's coefficient, 6.455e-3 | equals `−16π c_{2,0}/N` with `c_{2,0} = −N/(2304π√5 log φ)` from David et al. and 2606.30716 [C] |
| free Dirac fermion, constant term of `log Z` on S² | `−0.2310474` at β/r = 0.02, against `−log2/3 = −0.2310491` [C]; hence `c₁^{Dirac} = log 2/(24π) > 0` |

**Structural reason for `c₁(N = ∞) = 0`** [D].
* Reduce each Matsubara mode to 2D. Conformal coupling leaves the term `(1/6 − 1/8)R = R/24` in the 2D heat kernel.
  Summing over the modes (ζ-regularised) gives `c₁/N = log(2 sinh(βm/2))/(96π) = −(β/48)⟨φ²⟩_β^{ren}/N`.
* The critical gap equation is exactly `⟨φ²⟩^{ren} = 0`, i.e. `2 sinh(βm/2) = 1`.
* So `c₁ = 0` at N = ∞ is not an accident of O(N). It holds for every bosonic large-N theory whose saddle is the
  critical-boson gap equation, including CP^{N−1} at N = ∞ [I]. In each of them the sign is decided entirely by the
  fluctuation determinant.

**What the decisive computation is.**
* At order N⁰, `log Z₁ = −½ Tr log Π` on S¹_β×S²_r. Here `Π = G²` is the bubble of the N = ∞ propagator, which has
  eigenvalues `ω_n² + (l+½)²/r² + m²(r)`.
* In the (n, J) basis this is exact:
  `Π_{n,J} = (4πr²)⁻¹ Σ_{l₁,l₂} (2l₁+1)(2l₂+1) (l₁ l₂ J; 0 0 0)² K_n(E_{l₁}, E_{l₂})`, where K_n is the closed-form
  Matsubara sum of two propagators.
* Extract the β⁰ term of `−½ Σ_{n,J} (2J+1) log Π_{n,J}` as r/β → ∞.
* There is also a semi-analytic form [D]. By Euler–Maclaurin for J (exactly as in the probe):
  `χ-term = (1/3) Σ_n log Π_n^{flat}(k=0) + Σ_n ∫d²k (…) Π_n^{(1)}(k)/Π_n^{flat}(k)`.
  Here `Π^{(1)}` is the O(1/r²) part of the S² bubble: the `δm² = 1/12r²` term plus the next-order 3j asymptotics.
  Both n-sums are ζ-regularisable.
* **The UV divergences cannot fake the answer.** They are local, `β∫√g(Λ³ + ΛR)`, and so linear in β. They therefore
  separate from the β-independent `c₁`.

**Time to a decisive result.** About 2–3 weeks for one method, plus 1 week for an independent cross-check (S¹×H² by
SL(2,ℝ) harmonic analysis, or the semi-analytic route). Extending to CP^{N−1} (adding the photon determinant) takes about
1 more week.

**Risks.**
1. **Precision.** `8πc₁⁽⁰⁾` must be extracted under a `4πr²fT²` term. The probe shows the N = ∞ analogue converges as
   β²/r², which is fine, but the N⁰ bubble sums converge slowly (∝ 1/l).
2. **Subtle regularisation of `Tr log Π`.** The checks are that the β → ∞ Casimir energy vanishes and that the known 1/N
   correction to f is reproduced.
3. **Scoop.** The groups behind 2411.18509 (David et al.), 2508.14872, 2608.31151 and 2610.09815: a run of 2024–26
   papers on large-N O(N) on exactly this geometry. Also Mauro–Vichi.
4. **"Conjecture holds" is the likelier publishable-but-small outcome.**

---

## 2. Other candidates (ranked after §1)

### 2.1 Pelissetto–Vicari: `Δ_ε ≥ 2Δ_φ`, i.e. `ν ≥ 1/(2−η)` and `γ ≥ 1` (2510.17637, Oct 2025) [R]

* **Statement.** For continuous transitions described by LGW Φ⁴ theories with a unique quadratic invariant, extended by
  them to GNY and Abelian-Higgs (AH) models, `Δ_ε ≥ 2Δ_φ`. With η ≥ 0, this implies ν ≥ 1/2 for unitary theories.
* **Status.** Three citing papers (lattice gauge models), all consistent. Open.
* **Interest.** A counterexample would be a continuous transition with γ < 1. That is PRL-level in statistical
  mechanics and condensed matter.
* **Key observation** [R + D]. In every GNY model `Σ ≡ Δ_ε − 2Δ_φ = 0` at N_f = ∞ (`Δ_φ = 1`, `Δ_{φ²} = 2`), so the
  sign is decided at O(1/N_f).
  * The authors checked Ising (`8/π²N_f`), XY, Heisenberg (`4/π²N_f`) and QED₃-GNY (`8/π²N_f`). All are positive.
* **Probe** [D, uses Gracey 1801.01320 [R]].
  * Write the O(1/N) exponents of the chiral O(M) GNY family. The fermion self-energy is ∝ M and the vertex is
    ∝ (2−M)/M (Clifford algebra). The relation `λ₁ = −(2μ−1)η₁` holds for both M = 1 and M = 3.
  * Then `Σ = 2μη₁[2 + (2−M)/(M(μ−1))]/N = 12η₁/(MN)` at d = 3. Since η₁ ∝ M, this is **positive and independent of M**.
  * So no counterexample comes from this family. It still has to be confirmed with the M-dependence of the
    φ²-insertion diagrams.
* **Other near-misses.**
  * In large-N CP^{N−1} the gauge-invariant bilinear Q has `Δ_ε − 2Δ_Q = −16/(π²N) < 0` (from `ν = 1−48/π²N`,
    `η_N = 1−32/π²N` [M]). But P–V define the order parameter as the nonlocal Dirac-dressed φ, for which
    `Σ = 1 − 28/π²N`, so Q is outside their class.
  * The DQCP's ν ≈ 0.45 is pre-empted by their "crossover to first order" remark.
* **Tools.** Exact large-N computations.
* **Time.** 1–2 weeks to scan unexplored GNY variants: QED₃-chiral-Heisenberg, Sp(N) and SU(2)-triplet couplings
  without cubic invariants.
* **Odds.** P(counterexample) ≈ 5–10%. A proof is impossible: the class is not defined axiomatically.

### 2.2 Kusuki–Suzuyama–Tang: `c_cap ≥ c_eff` at conformal interfaces (2610.02309, Oct 2026) [R]

* **Statement.** The log coefficient of the capacity of entanglement across a 2D conformal interface is at least the
  entropy coefficient. A homogeneous CFT saturates it.
* **Status.** New. It is proved for free bosons and fermions (all transmissions), for holographic thin branes, and to
  leading order in defect perturbation theory. DMRG supports it at the W₃ TCI–Ising interface. The authors' open goal
  is verbatim: "prove the bound for general unitary conformal interfaces and identify the principle".
* **Reformulation** [D]. Write `log Tr ρⁿ ≈ −W·E(n)`, with W = 2 log(L/ε) and E(n) the excess ground-state energy of a
  circle of length 2πn carrying 2n interfaces. The conjecture is then `d²E/d(log n)² ≤ 0` at n = 1. The homogeneous
  case gives `E = (c/6) sinh(log n)`, which is marginal.
* **Interest.** High (PRL format).
* **Tools.** Partial: exact free-field computations, and positivity of spectral densities if a representation of E(n)
  as a positive mixture exists.
* **Odds.** Proof within weeks ≈ 5%. A counterexample is unlikely: the perturbative proof covers weak reflection.

### 2.3 Freivogel–Stoffels: positivity of `L_{2n} = ∫dx⁺ (x⁺)^{2n} T₊₊` in "physically preparable" states (2501.08386) [R]

* **Status.** One citation. Open, but the conjecture is not sharp: "preparable" is defined only qualitatively.
* **Probe** [D].
  * In a chiral 2D CFT, L₂ = K is conjugate to P under inversion, so L₂ ≥ 0 exactly.
  * For L₄, the O(ε) interference term in `|0⟩ + ε|2⟩` is a moment of a product of positive-frequency wavefunctions. It
    vanishes for smooth (IR-decaying) packets, so the sign is decided at O(ε²).
* **Tools.** Light-ray and exact free-field computations.
* **Assessment.** A decisive result is possible only after fixing a definition. PRL odds are low.

### 2.4 Hubeny–Liu: every holographic entropy inequality has a "tripartite form" (2609.22504, Sep 2026) [A]

* **Test.** Check representability of the known 5-party facets and the partial 6-party facet list by unit-coefficient
  I₃ and conditional-I₃ terms. This is an integer or LP feasibility check, which fits our integer-certificate tools,
  and takes days.
* **Context.** Similar combinatorial conjectures were settled within months by Grimaldi–Headrick–Hubeny–Shteyner
  2601.09987 and Ju–Czech 2611.21726. The field moves fast.
* **PRL odds.** About 2%.

### 2.5 Sokolovs: `I₃^{spin} > 0` for free fermions, all z = k_F w (2603.14154, Mar 2026) [A]

* **The open part.** A positivity conjecture for one scaling function. Small z is proven; the minimum is 0.100 at
  z ≈ 1.5.
* **How to settle it.** Certified interval evaluation of Gaussian-state entropies, plus a large-z asymptotic lemma.
  Ideal for our certified-positivity tools: 1–2 weeks.
* **Interest.** Narrow. PRL odds are below 3%.

### 2.6 Fernández-Sarmiento–Penco–Rosen: leading symmetric irrelevant deformations raise the entropy (2508.14978) [A]

* **Probe** [I]. For an interacting IR CFT, e.g. 2+1D quantum Ising with leading irrelevant operator ε′, the first-order
  entropy shift is ∝ `g·a_{ε′}`. The sign of g is non-universal: corrections-to-scaling amplitudes change sign across
  Blume–Capel-type families.
* **Consequence.** Unless the conjecture is restricted to IR theories with S-matrix positivity (free Goldstones or
  photons, which are their examples), it fails trivially. That would be a comment, not a PRL.

### 2.7 Bueno et al.: `C_T/F₀ ≤ [free scalar] ≈ 0.314` in d = 5 (2604.01436, Apr 2026) [R]

* **Probe** [D].
  * Einstein holography gives 0.0936 (their table).
  * Gauss–Bonnet multiplies the ratio by `(1−2λf_∞)/(1−4λf_∞)`. The conformal-collider tensor bound `t₂ ≤ 4` gives
    `λf_∞ ≤ 3/16`, so the factor is at most 2.5 and the ratio at most 0.234 < 0.314. **No violation.**
* **Further routes.** Cubic gravities might reach further. But a holographic theory with higher-curvature terms and no
  UV completion is a weak counterexample.
* **Proof.** Same global obstruction as `C_T/F` in d = 3 (already ruled out).

### 2.8 Checked and discarded (already settled, or not our tools)

| conjecture | status or reason |
|---|---|
| Li's total positivity of 1D blocks | **resolved**: Ouimet–Richards 2609.21210 [A] |
| Krylov bound `λ_L ≤ λ_K ≤ 2πT` and the operator-growth hypothesis | **violated**: 2506.01957, with tension noted in 2508.05444 [A] |
| Rényi area law (α < 1) | **disproved**: 2610.15199 [A] |
| Planckian lower bound on τ_eq | **proved** in regulated form: 2602.17638 [A] |
| CFT distance conjecture | 2D bounds proved (2405.00674); SharpDC violated in class S (2608.28333) [A] |
| thermal-EFT `c₂ ≤ 0` | **violated** by 4D Maxwell (2606.30716) [R] |
| spatial positivity of 2D partition functions (2610.30234) | RCFT classification, not our tools |
| Latorre–Sierra `c_max ≤ d−1` (2403.17242) | lattice constructions, not our tools |
| positivity of Zhukovsky-variable coefficients (2611.09867, footnote) | "plays no role" in their theorem; low value |
| line-defect `h ≥ 0` | pre-2024 |

---

## 3. Ranking

| rank | target | P(decisive, ≤ 4 wk) | P(PRL \| decisive) | tool fit | scoop risk |
|---|---|---|---|---|---|
| **1** | **`c₁⁽⁰⁾` sign in large-N O(N) (+ CP^{N−1})** | 0.6 | 0.35 if negative (P ≈ 0.4); 0.05 if positive | excellent (exact large-N on S¹×S² / H², certified sign) | medium |
| 2 | P–V `Δ_ε ≥ 2Δ_φ`: 1/N scan of GNY variants | 0.8 (a scan); 0.05–0.1 (a counterexample) | 0.3 | good | low |
| 3 | KST `c_cap ≥ c_eff` general proof | 0.05 | 0.4 | partial | medium (Kusuki–Tang active) |
| 4 | Sokolovs `I₃^{spin} > 0` certified | 0.7 | < 0.03 | excellent | low |
| 5 | Hubeny–Liu tripartite form | 0.5 | 0.02 | good (integer LP) | high (fast field) |
| 6 | Freivogel–Stoffels `L_{2n}` | 0.2 | 0.1 | good | low |
| 7 | entropy positivity (2508.14978) | 0.8 (comment) | < 0.02 | fair | low |
| 8 | d = 5 `C_T/F₀` | 0.1 (counterexample) | 0.1 | holographic only | low |

Expected PRL value: about 0.1 for rank 1, about 0.02–0.03 for ranks 2–3, and negligible below that. Rank 1 wins
clearly, subject to the adjacency flag in §0.

---

## 4. Recommended plan, first week: the sign of `c₁⁽⁰⁾` in the critical O(N) model

**Target statement.** Compute `c₁ = c₁⁽⁰⁾ + O(1/N)` for the 3D critical O(N) model, with
`−8π c₁⁽⁰⁾ = β⁰-coefficient of [−½ Tr log Π]` on S¹_β×S²_r as r/β → ∞. If `c₁⁽⁰⁾ < 0` with certified error, then
Allameh–Shaghoulian's conjecture fails for all N > N₀.

**Day 1 (done in part here).**
* The N = ∞ machinery on S¹×S² reproduces `m₀`, `δm² = 1/12r²`, `f`, `c₁ = 0` and `c_{2,0}` (§1.4).
* Still to do: write the closed-form Matsubara kernel `K_n(E₁, E₂)` and the 3j recursion for `(l₁ l₂ J; 000)²`.

**Days 2–3: the bubble `Π_{n,J}(β, r)`.**
* Exact sums for `l₁, l₂ ≤ L`, with an analytic large-l tail from the flat 3D bubble `1/(8p)` and its 1/l corrections.
* Checks:
  * the flat limit, against the known flat thermal N = ∞ bubble `Π_n(k)`;
  * the vacuum limit β → ∞, against the exact Δ = 2 σ-propagator on R×S² (conformal-map formula).

**Days 3–5: `Tr log Π` and the β⁰ term, by two routes.**
* **(A) Direct.** Compute `F₁(β, r) = −½ Σ (2J+1) log Π_{n,J}`, regulate, and fit at fixed β = 1, r = 8…40:
  `a r² + b + c/r²`. The divergent pieces are linear in β, so check β-independence of b at two values of β.
* **(B) Semi-analytic.** Use the Euler–Maclaurin decomposition
  `b = −½[(1/3)Σ_n log Π_n^{flat}(0) + Σ_n ∫ Π^{(1)}_n/Π^{flat}_n]_ζ`, with `Π^{(1)}` from `δm²` and the next-order
  3j asymptotics.
* **Calibration:**
  * zero Casimir energy at O(N⁰) on R×S²;
  * the coefficient a, compared with the literature 1/N correction to f (Chubukov–Sachdev–Ye) [M; to be checked].

**Day 6.**
* If |b| exceeds the combined error by 5σ, certify the sign: interval bounds on the truncation tails and the fit
  remainder (our Bernstein / interval tools).
* Start the S¹×H² cross-check: the Plancherel measure `p tanh(πp)/2π` with `χ_reg = 1` [I]. It is independent of the 3j
  route.

**Day 7: go / no-go.**
* **Go (b > 0, i.e. `c₁⁽⁰⁾ < 0`).** Draft "No sign rule for the thermal Euler term of 3D CFTs". Weeks 2–3: CP^{N−1}
  (photon determinant), an estimate of N₀, and the Rényi and Ψ interpretation.
* **Confirm (b < 0).** Write it up as the first interacting check (JHEP). Hand the S¹×H² code to the user's
  Herzog–Stergiou work, where the same determinant enters Ψ(x). Then switch to rank 2: run the P–V 1/N scan for
  1–2 weeks, as a cheap lottery ticket.
* **Kill.** If the two routes disagree by more than the error after day 7, stop and report. Do not publish a sign
  that is not cross-checked.

**Honest odds.** About 60% that we have a cross-checked sign by the end of week 3. About 40% that it is negative. About
10% that the whole effort ends in a PRL.

---

## Appendix: probe scripts (run as `python3 gap_s2.py; python3 c1_ninf.py` in one folder; about 2 min with mpmath)

`gap_s2.py`

```python
from mpmath import mp, mpf, sqrt, exp, log, quad, findroot, nsum, inf, pi
mp.dps = 30
beta = mpf(1)
def G(m, r):   # renormalised <phi^2>/N on S^1 x S^2 (conformal coupling), Abel-Plana for half-integers
    M = m*r
    vac = (-M + 2*quad(lambda t: t/(sqrt(M**2-t**2)*(exp(2*pi*t)+1)), [0, M]))/(4*pi*r)
    th = nsum(lambda l: (2*l+1)/(exp(beta*sqrt((l+0.5)**2/r**2+m**2))-1)/sqrt((l+0.5)**2/r**2+m**2), [0, inf])/(4*pi*r**2)
    return vac + th
m0 = 2*log((1+sqrt(5))/2)
for r in [10, 20, 40, 80]:
    m = findroot(lambda mm: G(mm, mpf(r)), m0); print(r, m, (m**2-m0**2)*r**2)      # -> 1/12
for x in [mpf('0.05'), mpf('0.02')]:                                                 # free Dirac on S^2
    lz = 4*nsum(lambda k: k*log(1+exp(-x*k)), [1, inf]); print(x, lz - 3*mp.zeta(3)/x**2, -log(2)/3)
```

`c1_ninf.py` (it reuses `G` from `gap_s2.py`)

```python
def logZ(m, r):   # log Z / N at N = infinity, zeta-regularised Casimir sum via Abel-Plana
    M = m*r
    casimir = (-(mpf(2)/3)*M**3 + 4*quad(lambda t: t*sqrt(M**2-t**2)/(exp(2*pi*t)+1), [0, M]))/r
    th = nsum(lambda l: (2*l+1)*log(1-exp(-sqrt((l+0.5)**2/r**2+m**2))), [0, inf])
    return -casimir/2 - th
for r in [10, 20, 40]:
    r = mpf(r); m = findroot(lambda mm: G(mm, r), m0)
    print(r, logZ(m, r) - 2*mp.zeta(3)/(5*pi)*4*pi*r**2)   # -> 0 like beta^2/r^2: c_1 = 0; slope = -16 pi c_{2,0}/N
```

**Output.**

| check | values | expected |
|---|---|---|
| `(m² − m₀²)r²` | 0.0834456, 0.0833613, 0.0833403, 0.0833351 | 1/12 |
| Dirac constant term | −0.2310386 (β/r = 0.05), −0.2310474 (β/r = 0.02) | −log2/3 = −0.2310491 |
| O(N) constant term | 6.47e-5, 1.61e-5, 4.03e-6 (r = 10, 20, 40) | 0 |
