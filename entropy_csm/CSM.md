# c_S/C_T <= 16πζ(3) on planar Chern–Simons–matter theories

Files: `csm_ratio.py` (formulas, gap equations, consistency checks; `python3 csm_ratio.py`),
`scan.py` (family scans; output in `scan_output.txt`, ~10 min), `R_vs_lambda.png`.
Notation: R := (c_S/C_T)/(16πζ(3)), so R = 1 for free scalars and R = 3/4 for a free Dirac fermion.

## Verdict

**R < 1 everywhere on every planar CS-vector family tested. No counterexample.** The only point where
sup R = 1 is reached is the free-scalar corner (λ, λ_6) → (0, 0) of the regular-boson family.
This is a theorem about the planar R^2 saddles, and the primary paper (arXiv 2610.33179, Herzog–Stergiou,
Supplemental Material "Planar Chern–Simons-matter theories") already proves it for the four
single-species theories. I re-derived it independently from the source formulas and extended it to
the boson+fermion theories (N=2, and the N=1 / "w" family of 1305.7235).

## 1. What the primary paper already contains

- It states that at planar order the bound holds analytically for the regular and critical boson and
  fermion theories. The ratio is written as c_S/C_T = 16π G(μ,φ), where
  `G(μ,φ) = (1/12)∫_μ^∞ t(t²−μ²)/(cosh t − cos φ) dt`. G decreases in both μ and φ, and G(0,0) = ζ(3).
- C_T = C_T(0) sin(πλ)/(πλ).
- It uses the table-top ("uniform") holonomy saddle of 1211.4843.
- It scans the regular-boson λ_6 stability domain of Dey et al. 1808.04415, with benchmarks at λ = 1/2.
- It excludes a competing Higgsed saddle at the conformal point.
- Its Table 1 gives an RB range of 45.316–60.422 and "other CS-matter < 60.422".
- It does **not** treat the boson+fermion theories (N=2 or N=1), 1/N corrections, or the S^1×H^2 response for CS-matter.

## 2. Independent derivation (computed, not assumed)

**Inputs taken from the literature.**
- Leading-N free energies on R^2 with the table-top holonomy ρ(α) = 1/(2πλ) for |α| ≤ πλ. This is the
  correct V T² ≫ N limit (1211.4843 §2, 1301.6169), which is the right order of limits for the conjecture.
- Gap equations from the 1211.4843 summary and 1305.7235 eqs. (girfreesquare) and (bosfree).
- Planar C_T = Ñ C_s per species, with Ñ = 2N sin(πλ)/(πλ) (1207.4593 eq. (sumtwo); 1211.1866 for fermions).
  Ñ is invariant under bosonization and independent of λ_4 and λ_6.

**Key identity.** I checked this both symbolically and numerically. For any real thermal mass c, the
on-shell planar free energy of a fundamental species is
`F = −N V T³ · sin(πλ)/(π²λ) · G(c, φ)`, with φ = πλ for a boson and φ = π(1−λ) for a fermion.
- This follows from 1305.7235 eq. (free1). Only the definition of the holonomy-averaged tadpoles is used:
  - S(c) = c/2 − Im Li₂(e^{−c+iπλ})/(πλ)
  - C(c) = c/2 − Im Li₂(−e^{−c+iπλ})/(πλ)
- The specific gap equation is **not** needed.
- The factor sin(πλ)/(πλ) is exactly the C_T factor, so it cancels. Each species therefore contributes
  G(c,φ)/ζ(3) to R, and R is the C_T-weighted average.
- Because the integrand is positive and G is monotone in μ and φ, R ≤ 1 for any real masses. Equality
  holds only at μ = φ = 0.
- Consequently, R < 1 holds **for every value of λ_6** in the unHiggsed branch, including values outside
  the stability window. The window details therefore cannot create a violation.

**Normalisation checks** (`csm_ratio.py`, all passed to about 1e-15):
- f for a real scalar is ζ(3)/(2π), and f for a Dirac fermion is 3ζ(3)/(4π). Both were computed from
  scratch, with no CS input.
- R(free Dirac) = 3/4.
- G from polylogarithms, G from the integral, and 1211.4843's 𝓕(μ,λ) by direct quadrature all agree.
- S and C agree with brute-force holonomy averages.
- At λ → 0: G(0,0) = ζ(3) (R = 1 for a complex scalar), G(0,π) = (3/4)ζ(3), and the critical O(2N) model gives R = 4/5 (Sachdev).
- The 1211.4843 perturbative results are reproduced: μ_F = 2 ln2 λ + …, and F_F agrees through O(λ²).
- The primary-paper benchmarks are reproduced: at λ = 1/2 I get c_S/C_T = 49.0337677 for λ_6 = 6π² and 48.8231125 for λ_6 = 12π².
- The ungauged BMB endpoint gives μ = log 2 and R = 7/8.

**Bosonization checks.**
- RF(λ_f) = CB(1−λ_f) holds identically, for both μ and R.
- CF(λ_f, ĝ) = RB(1−λ_f, λ̂ = 2λ_b/ĝ) holds; this map comes from the λ_6 map, 1211.4843 eq. (lambda6map).
- N=2 is self-dual under λ → 1−λ.
- The w = 1 member of the N=1 family reproduces N=2, with c_B = c_F.

## 3. Results (planar, 0 < λ < 1)

| family | R range | sup / where | λ→0 | λ→1 |
|---|---|---|---|---|
| regular fermion (= CB at 1−λ) | 0.75 → 0.80, increasing | 4/5, not attained (λ_f→1) | 3/4 (free Dirac) | 4/5 (critical O(2N)) |
| critical boson | 0.80 → 0.75, decreasing | 4/5 at λ_b→0 | 4/5 | 3/4 |
| regular boson, stable λ_6 | (3/4, 1) | 1 only at (λ,λ_6)→(0,0) | [7/8, 1] (free → BMB) | 3/4 at tip λ̂=2 |
| critical fermion (= RB) | same as RB | 1 as λ_f→1, ĝ→large | 3/4 (GN) | [7/8, 1] |
| N=2, one chiral | [0.7723, 7/8) | 7/8 at both ends | 7/8 | 7/8 |
| N=1 w-family (sampled w∈[−2,4]) | ≤ 0.867 sampled | < 7/8 | → 7/8 | → 7/8 |
| any 1 boson + 1 fermion, any real masses | ≤ [G(0,πλ)+G(0,π−πλ)]/(2ζ3) | 7/8 (λ→0,1) | | |

**Regular boson.**
- At fixed λ, R decreases monotonically in λ̂ = (λ² + λ_6/8π²)^{1/2}. This was checked on a 250×40 grid
  of the stable domain λ̂ ∈ [2λ/(2−λ), 2].
- The maximum at each λ is therefore on the lower stability boundary (x6_Dey = φ₁), and the minimum is
  on the upper boundary λ̂ = 2.
- Examples of the R range: [0.841, 0.883] at λ = 0.3, [0.802, 0.816] at λ = 0.5, and [0.770, 0.773] at λ = 0.7.
- No sampled point has R > 1. The grid maximum is 0.99990 at λ = 0.002.

**Danger region near the free scalar (analytic).** Expanding the polylogarithms at small u = μ − iφ, with ρ = |u|:

`G(μ,φ) = ζ(3) − (μ²+φ²)/6 · ln(1/ρ) + O(ρ²)`

The O(ρ²) term is −μ²/4 at φ = 0 and −(11/36 − ζ(3)/6)φ² at μ = 0.
- The C_T fall-off (1 − π²λ²/6) cancels exactly between f and C_T, because f ∝ sin(πλ)/(πλ)·G.
- The net deficit is therefore governed by G alone, and it is **positive** and log-enhanced.
- On the RB lower boundary λ̂ ≈ λ, the gap equation gives μ ≈ λ ln(1/λ). Hence
  1 − R ≈ λ² ln³(1/λ)/(6ζ(3)) at leading log. Numerically 1 − R = 7.3e-7 at λ = 1e-4, and it is positive at every λ tested.
- On the regular-boson family, f never falls more slowly than C_T near the free scalar.
- By contrast, the regular fermion has R = 3/4 + 0.1002 λ² + …, so here f does fall more slowly than C_T.
  This is harmless because it starts from 3/4. The coefficient 0.1002 comes analytically from the
  1211.4843 O(λ²) term and matches the numerics.

**Holonomy matters.**
- With the naive zero-holonomy formulas of 1110.4386, the regular fermion would give a **spurious
  violation** R > 1 for λ_f > 0.4405, with R → ∞ as λ_f → 1. The cause is that the naive f does not
  vanish like C_T. This is the dotted curve in `R_vs_lambda.png`.
- With the correct table-top holonomy, f ∝ sin(πλ)/λ, exactly like C_T.

**Saddle dominance.**
- The lowest F corresponds to the largest Σ G, so dominance selects the largest R among the competing
  saddles. The G bound applies to every real-mass saddle, so this cannot produce a violation.
- In the w-family scan, only one real solution was found at each sampled point. At (λ, w) = (0.95, −2) there was none.
- For the RB Higgsed branch I rely on the primary paper's exclusion. I did not re-derive it.
  Under duality that branch maps to fermionic G-form free energies (Dey et al.), so it would also obey the bound.

## 4. Status: computed vs assumed

**Computed here.** All R values, limits, the bosonization and duality checks, the small-λ expansion, and
the naive-holonomy contrast.

**Assumed from the literature.**
- The planar free energies and gap equations (1211.4843, 1305.7235).
- The table-top holonomy as the R^2 thermodynamic saddle (1301.6169).
- The Ñ formula for C_T (1207.4593, 1211.1866).
- For N=2 and the N=1 family: planar C_T equals the sum of the boson and fermion Ñ. Mixed boson–fermion
  contributions are O(N^0), and 1211.4843 §7.3 states that spin ≥ 1 planar correlators are unchanged by λ_4 and λ_6.
- The RB stability window from Dey et al. 1808.04415. It was only needed for the "max on boundary" statement, not for the bound.

**Not done.**
- 1/N corrections. At finite N, λ_6 is not marginal, and the planar deficit near (0,0) is only
  O(λ² ln³(1/λ)), so a sign statement at λ ~ 1/N is not controlled by this analysis.
- The S^1×H^2 heat-capacity version of the conjecture.
- SO/Sp groups and N_f > 1 flavours. The same G-average argument applies at fixed N_f, but I did not compute them separately.

**Clean statement.** For any planar CS-vector model with fundamental bosons and/or fermions, at any
marginal couplings with real thermal masses, on the table-top R^2 saddle:
R = ⟨G(μ_i, φ_i)⟩/ζ(3) ≤ 1, with equality only for free scalars.
For one boson + one fermion, R ≤ 7/8 (numerically on a λ grid). With unequal numbers of bosons and fermions, the bound is the weighted average of these envelopes.
