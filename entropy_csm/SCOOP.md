# Literature and scoop check: c_S/C_T <= 16 pi zeta(3) versus BBP T_M

Checked 2026-11-28. Sources: INSPIRE API (`refersto:recid:3218785` is 2610.33179; the `refersto:arxiv:` form of the query returns junk), Semantic Scholar citations, arXiv API searches over hep-th from 2026-10-23 to 2026-11-28, WebSearch in extended mode, and full-text greps of the PDFs listed below.

## Verdict

- **No scoop.** 2610.33179 has exactly one citer, Kusuki–Tang 2610.37449, on both INSPIRE and Semantic Scholar. That paper reports no violation, and it does not discuss moduli spaces, Higgs branches or monopole superpotentials.
- **No hep-th paper since Oct 23 reports a violation.** I found nothing that claims a violation via supersymmetry, moduli spaces or flat directions.
- **No paper questions the existence of T_M.** None of BBP's 92 INSPIRE citers challenges T_M. This was screened by title and spot checks, not by reading all 92.
- **The likely referee objection is pre-figured, not pre-empted.** HS themselves name flat directions and thermal ordering as the main loophole. They do so only in d = 4 − ε, where a flat direction gives saturation, not violation.
- **Bonus: the candidate also violates Kusuki–Tang, and more strongly.** This is my deduction, so check it. KT's conjecture f <= f_vac is stronger than HS for any 3d N>=1 SCFT. Running the numbers gives a violation of at least 1.217 for U(2) N_f=8 (versus 1.065 against HS) and at least ~1.38 for U(5) N_f=17. The derivation is in section 1.

## 1. Papers citing or responding to 2610.33179

| arXiv | Authors | Relevance | Scoops/undercuts? |
|---|---|---|---|
| 2610.33179 | C. P. Herzog, A. Stergiou | The conjecture itself (v1 only). Tested on free fields, O(N) at N=∞, the sextic/BMB line, Einstein holography, planar CS-matter, GN/QED3 at large N_f, and Monte Carlo data for Ising/O(N). It never mentions supersymmetry, Higgs branches or moduli-space EFTs in d=3. Main text and SM: "persistent thermal symmetry breaking can contribute with the opposite sign already at O(ε)"; in the large-n biconical model "an emergent flat direction leaves the thermal coefficient equal to that of free scalars ... saturating the bound"; "No general inequality for thermally ordered phases ... has been established." In planar CS-matter they explicitly exclude a competing Higgsed thermal saddle. | Not a scoop. It flags flat directions and thermal order as the soft spot, but finds no violation, and only in 4−ε. Cite this passage as the authors' own caveat. |
| 2610.37449 (v2, Nov 1) | Y. Kusuki, Y. Tang | Independent and overlapping work (its note added cites HS). Conjectures f <= f_vac, where f_vac is fixed by parity-even ⟨TTT⟩: in 3d f_vac = f_s (C_T/C_T,s)(28+t4)/32, eq. (S12). Also conjectures f <= 2.998 F (using BCLM) and f <= 4π²a in 4d. Its proofs assume "the symmetric thermal phase with equilibrium scalar background φ_a = 0" (SM §S2.C). It checks a free N=2 chiral (f = f_vac) and weakly coupled 4d SQCD. | Not a scoop, but the counterexample also hits it. See the deduction below. |

**Deduction to check: the candidate violates KT by more than it violates HS.**
- Buchbinder–Kuzenko–Samsonov (1503.04961) find that the 3d N=1,2,3 supercurrent three-point function has a single structure. If that is right, t4 = 0 for 3d N>=1 SCFTs, as for a free chiral.
- KT's bound then reads c_S/C_T <= 14πζ(3) for SCFTs.
- It is violated whenever C_T(origin) < C_T(free Higgs-branch chirals), with no 7/8 factor.
- U(2) N_f=8: 112/92.015 = 1.217. U(5) N_f=17: about 1.21 × 8/7 ≈ 1.38.
- KT's F-form bound, f <= 2.998 F, is *not* threatened. The free-chiral ratio is f/F = 1.932, and F(origin) exceeds F(Higgs EFT); in U(2) the values are 11.545 and 9.704.

## 2. Prior literature on moduli, thermal state and C_T monotonicity

| arXiv | Authors | Relevance | Scoops/undercuts? |
|---|---|---|---|
| 2005.03676 | Chai, Chaudhuri, Choi, Komargodski, Rabinovici, Smolkin | Thermal order in CFTs. At large rank the vector model has a moduli space and a "deformed moduli space of vacua at finite temperature", and the thermal state sits away from the origin. | Precedent that the thermal state need not be at the origin. It is non-SUSY, with no c_S/C_T statement. Neither scoops nor undercuts. |
| 2608.02720 | Smolkin, Yung | Biconical model at large N in 2+1d. Z2 breaking persists to all T. HS's ref [38]. | Same as above. |
| 2609.00160 | Chaudhuri, Hawashin, Rabinovici, Scherer | Thermal order with multiple scalars coupled to O(N). Says persistent breaking "could, in principle, occur for any quantum potential with flat directions." | Supports the mechanism in general terms only. |
| hep-th/0504151 | Apreda, Erdmenger, Evans, Guralnik | Strongly coupled 4d N=2 (D3/D7). Finite T generates a potential on the Higgs branch, and "for T>Tc the theory is driven to a non-trivial point on the Higgs branch". | Closest precedent for a thermal state on a Higgs branch at strong coupling. It concerns the phase, not a bound. No scoop. |
| hep-th/9809032 | Tseytlin, Yankielowicz | N=4 SYM on the Coulomb branch at finite T. The probe D3 is attracted to the horizon, so the origin is preferred. | Typical holographic case: the origin has far more degrees of freedom than the EFT. Consistent with, not against, the candidate. |
| 2306.08031 | Benjamin, Lee, Ooguri, Simmons-Duffin | The thermal EFT assumes the theory is gapped on S¹. Says this "fails when a symmetry protects gapless modes ... supersymmetric compactifications ... spontaneously break a continuous symmetry". | The standard thermal-EFT/Cardy machinery that HS lean on excludes ungapped cases like this one. That helps the candidate's framing. |
| 2402.13337 | Allameh, Shaghoulian | Same gapped assumption. Notes free-scalar KK zero modes. | Context only. |
| 1508.02728 | Shaghoulian | The Casimir swap c_S ↔ ε_vac (NS sector), which HS use in eq. (7). | This is the identity the variational (Casimir-swapped) argument uses. Fine. |
| 1303.1522 | Nishioka, Yonekura | τ_RR ∝ C_T can *increase* along a 3d N=2 RG flow (a WZ example): "rules out ... C_T ... monotonically decreasing". | Established prior art that C_T(UV) < C_T(IR) is possible. Cite it. Not a scoop: it is a deformation flow, not a moduli-space flow, and makes no thermal claim. |
| hep-th/9708042 | Anselmi, Freedman, Grisaru, Johansen | In 4d, c is not monotonic; a is. | 4d analogue of the above. |
| 1910.03600 | Chang, Fluder, Lin, Shao, Wang | Localization C_T for 3d N=4 theories whose Higgs branch is a minimal nilpotent orbit (C_T(hyper)=12). Every entry exceeds the free Higgs-branch value: SQED_{N=2..5} gives 22.9/34.5/46.3/58.2 versus 12/24/36/48; SU(2) N_f=4 gives 84.9 versus 60; T3 gives 160.2 versus 132. | No prior case of C_T(origin) < C_T(moduli EFT) here. The candidate's monopole-W N=2 example appears to be new. |
| 0809.3238 | Shapere, Tachikawa | Claimed a_UV < a_IR on a 4d N=2 moduli space. v3 retracts it: the low-energy theory at the special point had decoupled sectors (see 1011.4568). | **Cautionary precedent.** A "moduli-space counterexample" died over a misidentified IR theory at the origin. The candidate's check for no accidental symmetry or decoupled sector at T_M (R_Q = 5/8) is exactly the point a referee will press. |
| hep-th/9901109 | Appelquist, Cohen, Schmaltz | Conjecture f_IR <= f_UV. | Same direction as the pressure inequality f(T_M) >= f(Higgs EFT). Supportive, with no C_T normalisation. |

I found no prior paper that computes the thermal free energy of a 3d N=2 gauge theory on its Higgs branch. I also found no published C_T or τ_RR for the BBP T_M theories.

On the 4d analogues: the f/a, f/c and f_vac bounds are in KT 2610.37449. The f <= 4π²a form cannot be violated via moduli, because a_origin > a_EFT (Komargodski–Schwimmer) and the EFT itself obeys the bound. I found no paper reporting a 4d violation via moduli.

## 3. Older bounds of the same type

| arXiv | Authors | Relevance | Scoops/undercuts? |
|---|---|---|---|
| 2307.05164 | Bueno, Casini, Lasso Andino, Moreno | C_T/F is maximised by the free scalar (BCLM). HS's Ψ(1) uses it, and KT's F-bound depends on it. None of its 16 INSPIRE citers tests a SUSY or moduli theory. | Not threatened by the candidate: C_T(origin) is small and F(origin) is large, which is the safe direction. |
| 0806.0110 | Kovtun, Ritz | Central charges versus thermal and transport data, with possible bounds suggested by holography. | Ancestor in spirit. No moduli discussion. |
| 0712.0743 | Kats, Petrov | Finite-N violation of the η/s >= 1/4π bound (c ≠ a). | Precedent for a "universal" thermal bound failing in a legitimate theory. Useful framing only. |
| hep-th/0101076 | Klemm, Petkou, Siopsis | Entropy/energy ratios bounded in free CFTs. | Historical context only. |
| 2412.04122, 2501.09886 | R.-X. Miao | Holographic lower bound on the Casimir effect, for BCFT strips. | Different (boundary) setup. Not relevant. |

None of these older bounds was treated for theories with moduli spaces in the papers I read.

## 4. Existence of the BBP fixed points T_M

| arXiv | Authors | Relevance | Scoops/undercuts? |
|---|---|---|---|
| 1703.08460 | Benini, Benvenuti, Pasquetti | W_mon is relevant at T0 only for 2–3 values of N_f per N_c. Existence for N_f up to 3N_c+3 is argued via the Ising-SCFT chain (4.9) at N_f = 3N_c+2, then mass deformations to lower N_f; 3(N_c+1) is a conformal manifold. Two caveats are in their own text. (i) The first flow could instead end at T_{ΦQQ̃} with R[Φ] > 2/3: "With enough Φij it might be possible to use similar arguments to rule out the second flow, but we will not do it here". (ii) With too few Ising copies, the first flow violates the F-theorem. They also note that T_M is interacting for N_f > 4(N_c+1)/3. | Both U(2) N_f=8 and U(5) N_f=17 sit exactly at N_f = 3N_c+2, the endpoint of the chain, so no mass deformation is needed. Existence there still rests on caveat (i), which BBP leave open. State this explicitly. |
| Citers of 1703.08460: 1806.01356, 1901.07559 (Amariti et al.), 2012.09773, 2012.08556 (Benvenuti et al.), 1808.09983 (Amariti, SUSY breaking in monopole quivers), 2507.04950, 2610.22054, 2610.34752, etc. | various | These use or extend monopole-deformed dualities. 1808.09983 is about SUSY breaking in other quivers, not U(N_c) SQCD with V₊+V₋. | None questions T_M for N_f <= 3N_c+2. No undercut found. |

## Recommendations for the write-up

1. Quote HS's own flat-direction caveat. Frame the result as the first d=3 realisation of that loophole, and as a *violation* rather than the saturation they found.
2. Add the KT corollary, assuming t4 = 0 for N=2 is confirmed. It is the stronger statement, and KT explicitly assume φ_a = 0.
3. Cite Nishioka–Yonekura for C_T non-monotonicity, and Chai et al. and Apreda et al. for thermal states off the origin.
4. Pre-empt the Shapere–Tachikawa-type objection: list all chiral-ring operators at T_M with R >= 1/2 and argue there is no decoupled sector.
5. Flag that BBP's chain leaves one flow ambiguity open.
