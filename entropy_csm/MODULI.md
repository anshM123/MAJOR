# Moduli spaces and the bound c_S/C_T ≤ 16πζ(3)

Files: `moduli_check.py`, with its output in `moduli_output.txt`.

Units used throughout:
- C_T is quoted in free-real-scalar units, C_T/C_s.
- c_S is quoted in real-scalar units.
- A free N=2 chiral multiplet has c_S = 3.5 and C_T = 4.
- A free hypermultiplet has c_S = 7 and C_T = 8.
- A free N=8 multiplet has c_S = 14 and C_T = 16.
- Every free supermultiplet therefore has c_S/C_T equal to 7/8 of the free-scalar value.

## 1. What 2610.33179 states

Line numbers refer to `thermal_envelope_arXiv.tex`, which I read in LaTeX.

**The only hypothesis is unitarity.**
- L23 (abstract): "We conjecture that every unitary three-dimensional conformal field theory obeys $c_S/C_T\leq16\pi\zeta(3)$".
- L30–33 state the boxed eq. (1) with no further qualification.
- Nowhere does the paper exclude moduli spaces, supersymmetric theories, or degenerate vacua. Nowhere does it require a unique or gapped vacuum.

**c_S is defined via S¹×S² with no vacuum selection.**
- L326: "Place the theory on $S^1_\beta\times S^2$ … In the high-temperature limit, $c_S=\frac{3}{4\pi}\lim_{\beta\to0}\beta^2\log\sum_{\mathcal O}e^{-\beta\Delta_{\mathcal O}}$".
- L330: "where the sum includes all local operators, descendants, and spin components, with multiplicity."
- This is the full radial-quantization trace. On S², the conformal coupling makes the vacuum unique, so nothing is projected out. In particular, nothing excludes runaway along flat directions.
- L65 and L69 give an equivalent Casimir-energy form: "Both quantities in Eq. (casimirbound) are vacuum data."

**Supersymmetric examples are not discussed.**
- The words supersymmetry, N=2/4/8, ABJM, BLG, Coulomb branch and Higgs branch do not occur in the paper.
- Holography is invoked only as "Einstein holography" (L50, L180, L218).

**Flat directions are mentioned only in the near-4d multiscalar supplement.** These passages explicitly flag thermal ordering and emergent flat directions as a loophole:
- L275: "persistent thermal symmetry breaking can contribute with the opposite sign already at $O(\varepsilon)$ … an emergent flat direction leaves the thermal coefficient equal to that of free scalars at leading order, though the CFT remains interacting".
- L508: "Ordering *increases* $c_S$ already at $O(\varepsilon)$, so the argument … does not cover ordered phases."
- L530: "whose zero-temperature flat direction $3\sigma^2=\boldsymbol\phi^2$, related to the large-rank moduli space of Ref. [Chai2020], lets the condensate cancel the leading thermal masses".
- L532: "None of these flat-space results yet establishes an inequality for the response on $S^1\times H^{d-1}$".
- L402 and L606 say saddles are assumed and that competing phases are not excluded.

**The thermal state far out on a moduli space is not discussed** for 3d CFTs.

## 2. (a) Is p(T) ≥ p_moduli valid?

**Yes, in flat space.**
- The thermal free energy density is independent of boundary conditions in the thermodynamic limit, and it equals the infimum of the effective free energy over homogeneous configurations: p = sup_v p_eff(v).
- A rigorous version uses the Gibbs variational principle on a finite torus. The trial state is a broad zero-mode wavepacket centred at v. In finite volume there are no superselection sectors, so this state lies in the same Hilbert space. Its zero-mode kinetic cost per unit volume goes to 0.
- By scale invariance, p_eff(v) = T³ h(|v|²/T).
- As |v|²/T → ∞, the massive sector has mass of order |v|² ≫ T, and the sigma-model curvature couplings are suppressed by T/|v|². Hence h(∞) is the free IR theory on the moduli space.
- So c_S ≥ c_S(IR on moduli).

**Yes, also on S² (the paper's definition).**
- The conformal mass |φ|²/(4R²) gives the zero mode an action πβ|v|². Write ν = |v|²/T, and set R = 1.
- The weight is then exp[4π(RT)² h(ν) − πν].
- At RT → ∞, the thermal term of order (RT)² beats the O(ν) conformal cost for any ν ≪ (RT)². Taking ν → ∞ slowly gives log Z ≥ 4π(RT)² h(∞)(1 − o(1)).
- So the spectral c_S is also at least h(∞). The conformal coupling lifts the moduli only by a subleading amount.
- In operator language: the S² spectrum of the origin CFT must contain the free-field towers of the moduli theory, as heavy states with large |v|.

**Consequence.** For an SCFT whose generic moduli point carries n free chirals, the conjecture implies a criterion: **C_T(origin) ≥ 3.5n = (7/8)·C_T(IR moduli theory).** Hypers and N=8 multiplets give the same 7/8 factor.

## 3. (b) Test of the criterion on known theories

Data sources:
- C_T for N=4 theories: 1910.03600 (eq. SQED_n, its CTCJ table, T[SU(3)] = 75.53). That paper's convention has a real scalar at 3/2, so I divided by 1.5.
- c_T for N=8 theories: Chester et al. 1412.0334, their table cTValues.
- XYZ: computed here from the squashed-sphere double sine, ∂²_b F. The code reproduces F(free chiral) = log2/2 and C_T(free chiral) = 4. It gives C_T(R=2/3 chiral)/C_T(free chiral) = 0.72679, which agrees with the W=Φ³ value I recall from the literature (about 0.7268). I did not re-fetch that reference.

| theory | C_T(origin) | IR on generic moduli point | C_T(IR) | ratio C_T(0)/C_T(IR) |
|---|---|---|---|---|
| XYZ (= N=2 SQED_1) | 8.721 | 1 chiral | 4 | 2.18 |
| N=4 SQED_2 = T[SU(2)] | 15.30 | 1 hyper | 8 | 1.91 |
| N=4 SQED_5 / SQED_10 / SQED_50 | 38.8 / 78.6 / 398.5 | n−1 hypers | 32 / 72 / 392 | 1.21 / 1.09 / 1.017 |
| T[SU(3)] | 50.36 | 3 hypers | 24 | 2.10 |
| minimal SO(8) SCFT | 56.59 | 5 hypers | 40 | 1.42 |
| T_3 (E6 minimal) | 106.83 | 11 hypers | 88 | 1.21 |
| U(1)_k×U(1)_{−k} ABJM | 16 | 1 N=8 multiplet | 16 | 1.000 (free) |
| U(2)_2×U(1)_{−2} ABJ | 21.33 | 1 N=8 multiplet | 16 | 1.33 |
| U(2)_1×U(2)_{−1} ABJM | 37.33 | 2 N=8 multiplets | 32 | 1.17 |
| U(2)_2×U(2)_{−2} ABJM | 42.67 | 2 | 32 | 1.33 |
| BLG_3 / BLG_4 / BLG_5 | 47.0 / 50.4 / 52.9 | 2 | 32 | 1.47 / 1.57 / 1.65 |
| large-N ABJM, T[SU(N)] | ∝ N^{3/2}, ∝ N² log N | N multiplets, ∝ N² hypers | | → ∞ |
| **candidate** | **92–93** | **28 chirals** | **112** | **0.821–0.830** |

- Every known theory has C_T(origin) ≥ C_T(IR moduli), with the minimum of 1 reached only by free theories. This is stronger than the needed 7/8.
- All of them therefore pass. The induced lower bound R ≥ (7/8)/ratio is at most 0.875, attained only by free U(1)×U(1) ABJM.
- **The candidate would be the first SCFT with C_T(origin) < C_T(IR moduli).** No theorem forbids this, because C_T is not monotonic in 3d. But it has no precedent among the examples here.
- If the premises hold, then R ≥ 98/93 = 1.054 (up to 98/92 = 1.065). That would be a genuine counterexample, and none of the paper's stated hypotheses would escape it.

## 4. (c) F-theorem check for the candidate

- Free fields give F_scalar = 0.0638071, F_Dirac = 0.2189595 and F_chiral = log2/2 = 0.3465736.
- F of the IR moduli theory is 28 × 0.3465736 = **9.7041**.
- If the F-theorem applies to moduli-space flows, the candidate needs **F(origin) > 9.704**. The entropic proof needs only a Poincaré-invariant vacuum interpolating between the two fixed points, which a vev provides.
- If the localized F of the candidate is ≤ 9.70, its premises are inconsistent. Possible causes are a wrong superconformal R-charge, quantum lifting of the moduli, or a misidentified C_T.
- For comparison, the Bueno–Casini C_T/F bound only requires F ≥ 92·0.0638 = 5.87, which is weaker.
- The candidate's C_T/F would be at most about 9.5. That is below the free chiral's 11.5 but above holographic ABJM's 64/π² ≈ 6.5, so it is not implausible by itself.

## 5. What must be triple-checked before claiming the counterexample

1. **Accidental symmetries.** C_T from τ_RR or ∂²_b F is valid only with the true superconformal R-symmetry. Every gauge-invariant chiral operator must have R ≥ 1/2. Any operator that violates this decouples as a free field, and both C_T and F must then be recomputed. Theories with large Higgs branches often have exactly this issue.
2. **The moduli space is not quantum-lifted.** N=2 Higgs and mixed branches can be modified by monopole superpotentials or quantum chiral-ring relations. The 28-dimensional count must hold in the IR, not only at the classical D- and F-flat level.
3. **C_T is the total C_T of the full theory at the origin**, including any decoupled free or topological sectors. The moduli count must refer to that same theory.
4. **F(origin) > 9.704**, as in §4.

If items 1–4 hold, the argument in §2 is robust in both flat space and on S²: R ≥ 1.054.

**Correction, added during the CANDIDATE_U2 audit.** I cross-checked 1910.03600's closed forms for SU(2) with k half-hypers (the SO(5) and SO(8) entries) using `u2_ct.py`. They reproduce only if the vector-multiplet term uses root σ instead of 2σ. With the correct root, SU(2) with 3 flavours gives C_T = 30.8911, exactly equal to SQED_4, as the Sp(1)_3 = U(1)_4 equivalence requires; the 1910 table instead gives 39.5 for SO(6) versus 30.9 for SU(4).

Corrected values:
- The minimal SO(8) SCFT has C_T = 50.84, not 56.59. Its ratio is then 1.27 instead of 1.42, which still passes.
- The SO(5) row (k = 5) is SU(2) with an odd number of half-hypers, which has a parity anomaly without a Chern–Simons term. That row should be discarded.

None of the conclusions above change.
