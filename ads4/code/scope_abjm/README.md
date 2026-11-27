# Scope study: can the quartic bound be inverted into an ABJM "string-gap" headline?

Script: `invert.py` (mpmath). Its output is in `invert_output.txt`. Labels: **[paper]** = value stated in letter/supplement,
**[cited]** = from the literature (equation numbers given), **[computed]** = computed here from those,
**[estimated]** = order-of-magnitude or O(1) statements not computed.

## 1. The bound B(Delta_g) and which value applies

Statement (supplement item 2a): `g4~ Lambda^6/8piG <= 2 L_AdS =: B(Lambda)`, with `Lambda = [(Delta_g-2)(Delta_g-1)]^{1/2}`.

| Delta_g | Lambda | B, 2a (scalar condition) | B, 2b (no scalar condition) | status in paper |
|---|---|---|---|---|
| 300 | 298.50 | 72.041 [paper] | 125.71 [paper] | quartic headline only if HLJclosed=true; one hypothesis open |
| 500 | 498.50 | 80.855 [paper, tab:dg] | 142.15 [computed from tab:dg] | value of the implication only |
| 1000 | 998.50 | **92.814** [paper] = QuartHead | **164.51** [paper] = QuartRep | the only claimed quartic gap (Pos still assumed) |
| 2000 | 1998.5 | 104.769 [paper, tab:dg] | 186.90 [computed] | value of the implication only |
| 5000 | 4998.5 | 120.571 [paper, tab:dg] | 216.53 [computed] | value of the implication only |

Asymptotic form [paper, eq:quarticL]:
`B(Lambda) = 17.2407 ln(2Lambda/3) - 19.2782 + 10.314/Lambda + (860.3 ln(2Lambda/3) - 2508.6)/Lambda^2 + ...`. It reproduces the table to 1e-3.
The 2b values at 500/2000/5000 use `2L_AdS + 0.005*(16 L_A)` (eps = 0.040, eq:repair). They are not in the paper, and eps was tuned only at 300/1000.

Correction to the brief: "164.5 vs 125.7" is not a comparison of with versus without the scalar condition. Both numbers are the
*without* value (2b), at Delta_g = 1000 and 300 respectively (the HLJclosed switch). The *with*-condition values are 92.8 (1000) and 72.0 (300).
For ABJM the paper argues the scalar condition holds [D-sketch, supplement sec:abjm], so 2a applies: **B = 92.81 at 1000**.
`tab:gap2000` concerns the cubic functional Omega^EE near the gap. It contains no quartic data. The quartic Omega^Q_up is examined only at 300 and 1000.

## 2. ABJM inputs (cited)

* **Dictionary** [cited]: `nu := L^4/l_s^4 = 2 pi^2 (lambda - 1/24) + ...`, `g_s^2 = 512 lambda^2/(3 c_T)`. Source: Chester-Hansen-Zhong,
  arXiv:2412.08689, eq.(1.1), citing ABJM 0806.1218, Bergman-Hirano 0902.1743 and Aharony-Hashimoto-Hirano-Ouyang 0906.2390.
  Equivalently `L^2/alpha' = pi sqrt(2 lambda)` [BCP 1906.07195, eq.(4.11): L^8/l_s^8 = 4 pi^4 lambda^2].
  Hence `M L = 2 nu^{1/4} = 2^{5/4} pi^{1/2} lambda^{1/4} = 4.2156 lambda^{1/4}` with `M^2 = 4/alpha'` [computed].
* **c_T and G** [cited/computed]: `8 pi G/L^2 = 256/c_T` (holographic C_T with free real scalar c_T = 1, BCP eq.(1.3)). Leading order:
  `c_T = (64/3pi) sqrt(2k) N^{3/2} = (64 sqrt2/3pi) N^2/sqrt(lambda)`. This is consistent with CHZ (1.1) and the ABJM g_s.
* **g4 from localization** [cited]: in the 't Hooft limit, BCP 1906.07195 eq.(1.2) give
  `M = (1/c_T)[(32/pi^2) M_SG + (3 sqrt2 zeta(3)/4pi^5)(35 M_4 - 72 M_3) lambda^{-3/2} + O(lambda^{-5/2})] + O(c_T^-2)`. Equivalently,
  CHZ eq.(2.8): `M_tree = M_SG + 3zeta(3)/(32 nu^{3/2}) [35 M^{4,1} - 72 M^2] + ...`. The coefficients B_3^4, B_4^4 of eq.(4.31) are fixed by
  the localization constraints (∂^4_m log Z, Fermi gas, eq.(3.43)) together with the superblock constraints. They agree with the flat limit
  `A_IIA = A_SG(1 + l_s^6 zeta(3)/32 stu + ...)` (eq.(1.4)). In the Letter's normalization this is
  `g4 M^6/8piG = 2 zeta(3)`, and in AdS units **`g4/(8 pi G L^6) = zeta(3)/(32 nu^{3/2}) = 4.283e-4 lambda^{-3/2}`** [computed].
  Note: localization fixes the <SSSS> contact term. Transferring it to <TTTT> uses N=6 superconformal symmetry and the flat limit. For N=6 that
  transfer is the paper's [D-sketch] (supplement sec:abjm). CPY 1804.00949 and ACP 1711.07343 / 1907.13222 are M-theory (k = 1, 2), which has no string scale.
* **First massive level** [cited]: CHZ eq.(3.3)-(3.4): `Delta(Delta-3) = R^2 m^2 = 4 delta R^2/alpha'`, so the twist is `2 nu^{1/4} sqrt(delta) + O(1)`.
  Leading Regge trajectory (CHZ eq.(1.6)): superprimary `Delta = -3/2 + sqrt(2(l+1)) nu^{1/4}(1 + O(nu^{-1/2}))` (odd l, Z+). The
  integrability checks are Gromov-Sizov 1403.1894 (slope function) and Bombardelli-Cavaglià-Conti-Tateo 1803.04748.
  The flat first massive level has maximal spin 4. Its residue in the (+-) graviton channel is pure d^4_{44} [computed from VS: the residue at S=1 is constant].
  So **the lightest spin>=4 single trace in TxT has twist `4.2156 lambda^{1/4} + O(1)`**. The O(1) term (paper: about -4; the superblock gives
  -3/2 + spin shifts) was not computed here [estimated]. It moves Y by less than 1% at lambda ~ 1e8.
* **Check of "Delta_g = 1000 <-> lambda ~ 3e9"** [computed]: Lambda = ML gives lambda = 3.147e9. Twist(J=4) = ML - 4 = 1000 gives 3.217e9. **Correct.**
  (At Delta_g = 300: 2.51e7 / 2.70e7.)

## 3. The inversion [computed]

For ABJM, `g4 Lambda^6/8piG = 2 zeta(3) (Lambda/ML)^6`. This grows like Lambda^6, while B grows like ln Lambda, so the bound fails at large Lambda
and therefore gives an **upper** bound on the gap.

**(A) Fixed gap. This is the only form the paper's claims support.** If ABJM's twist gap were >= 1000, the bound at 1000 would apply, and it is violated iff
`ML < 998.5 (2zeta3/92.81)^{1/6} = 543.1`, i.e. **lambda < lambda_* = 2.76e8**. So, for 1 << lambda < 2.76e8 (IIA regime), ABJM must contain a
spin>=4 single trace in TxT with twist < 1000. At lambda_*, the true twist is about 539, a ratio of **1.855**. For lambda << lambda_* the statement weakens:
at lambda = 1e6 the true twist is about 129 and the bound is still 1000, a factor 7.7.
Without the scalar condition (2b): lambda_* = 1.88e8, ratio 2.04. At 300 (HLJclosed): lambda_* = 2.6e6, ratio 1.81.

**(B) Running gap. This needs (Pos) for Omega^Q at the gap Delta_max(lambda) itself.** Solve `2zeta3 (Lambda/ML)^6 = B(Lambda)`:

`Delta_max(lambda) ≈ 4.2156 Y(lambda) lambda^{1/4}`, with `Y = [ (17.24 ln(2 Y ML/3) - 19.28)/(2 zeta(3)) ]^{1/6}`

| lambda | true twist ~ML-4 | Delta_max | Y = Lambda_max/ML | Delta_max/lambda^{1/4} |
|---|---|---|---|---|
| 1e6 | 129 | 234 | 1.74 | 7.40 |
| 1e7 | 233 | 425 | 1.79 | 7.56 |
| 2.75e8 | 539 | 1000 | 1.84 | 7.76 |
| 3.15e9 | 995 | 1872 | 1.87 | 7.90 |
| 1e14 | 13327 | 26531 | 1.99 | 8.39 |
| 1e20 | 4.2e5 | 8.9e5 | 2.11 | 8.88 |

(Rows below Lambda ~ 300 use eq:quarticL outside its validated range.) The rough estimate (1/0.026)^{1/6} ≈ 1.8 is **verified**:
Y = (92.81/2.404)^{1/6} = 1.838, or 1.855 in twist, and 2.02 without the scalar condition.
Y is not constant: the 4d infrared log makes it grow as (ln lambda)^{1/6}, so "X lambda^{1/4}" holds only with a slowly running X ≈ 7.4-8.4.
Sensitivity: Y = (B/2zeta3)^{1/6}. A factor 2 in B changes Y by 12%, and a factor 10 by 47%.

Comparison of constants [computed]: 10d flat window [0, 3.000] (CHMRS-swamp) gives Y_10d = 1.038. The AdS5/N=4 bound alpha/8piG <= 10.9/Delta_gap^6
(CHMRS, arXiv:2106.10274, eq.(5.49)) gives Y_AdS5 = 1.29. The AdS4 value 1.84 is the weakest of the three, because of the 4d infrared log.

**Validity.** IIA 't Hooft regime: lambda >> 1, g_s^2 = 512 lambda^2/(3c_T) << 1 (k >> lambda^{1/4}), c_T -> infinity at fixed lambda (A2).
The paper's 1/c_T corrections, of size Delta_g^p/c_T^2, are not estimated. The M-theory regime (k fixed) is excluded: the inversion would put
Lambda_max ~ 5.7 L/l_11, at the 11d Planck scale. There the KK-tower loop parameter (E l_11)^9 is O(1), so A2 fails, and no single-trace higher-spin
tower exists to be bounded.

## 4. Assumptions of the inverted statement and their grades in the paper

| assumption | grade in paper | needed where |
|---|---|---|
| (Pos) of Omega^Q_up on every heavy state | **[A]**: assumption at every gap. Error-controlled only on binding near-gap states and masses 1.3-3 Lambda at 1000; elsewhere modelled or interpolated. Examined only at 300 and 1000. | (A) at 1000 only. (B) at Delta_max(lambda) for each lambda, i.e. at every gap in the quoted range. The brief's "at all gaps above Delta_max" is not quite right: one gap per claim is enough. |
| g4~ = g4 (L7(b)/H_EFT(b)) | **[A]** | both |
| scalar condition (no scalar with m < 0.46 Lambda in TxT) | [D-sketch] for ABJM (exchange argument, N=6 transfer) | 2a. Otherwise use 2b (B = 164.5, Y = 2.02). It needs (Pos) for Omega^EE too, which is certified at 1000. |
| tau = 0 or kappa > 0 | ABJM: ghat_3 = 0 by 10d SUSY at leading order | both |
| A2: large c_T, tree level below the gap uniformly in spin | standard input; 1/c_T corrections not estimated | both |
| g_5 term and light single traces in the light side | [D-sketch] (I4) | both |
| SSSS -> TTTT transfer of the localization R^4 coefficient | [D-sketch] in the paper (N=6 fixes 2 of 5 structures) | input g4 |
| Lambda <-> ML | **not needed for the inverted statement.** Delta_max is a twist bound. ML enters only the comparison Y, with O(1/ML) ≈ 1% error. | comparison only |

**Logic.** Suppose twist gap >= Delta' >= Delta_max. A theory with gap Delta_true also has every smaller gap, so the bound at Delta' applies.
For Lambda' > Lambda_max, `2zeta3 (Lambda'/ML)^6 > B(Lambda')`, which is a contradiction. So the gap is < Delta'. Direction checks:
2zeta(3) > 0 lies inside [0, B] (same normalization as CHMRS). Upper bound on g4 means upper bound on the gap. Lambda enters with the 6th power.
**The logic is sound**, given the assumptions. Two caveats. First, the contradiction only shows that one of {gap >= Delta', (Pos), L7(b), scalar condition, A2} fails.
Concluding "gap < Delta'" treats the others as independent of the spectrum. That holds for (Pos), a property of the functional on blocks, but
L7(b) is an EFT-naturalness assumption. Second, the statement is conditional on (Pos) for Omega^Q, which is the least-established positivity in the paper.

## 5. Novelty [cited]

* **CHMRS, "AdS bulk locality from sharp CFT bounds", arXiv:2106.10274, JHEP 11 (2021) 164, Sec. 5.3**: for N=4 SYM they give
  `0 <= alpha/8piG <= 10.9/Delta_gap^6` and juxtapose it with `Delta_gap ≈ 2 lambda^{1/4}` from integrability (GKV) and the localization R^4 value
  (BCPW 1902.06263): "establishes bulk locality in the canonical example of holographic duality, with a sharp bound that depends on the coupling lambda."
  They do not write the inversion, but it follows in one line and gives Y = 1.29, sharper than ours. This is the closest precedent and largely
  removes the conceptual novelty of the ABJM version. Their bound is a flat-space uplift with a truncation assumption (their eq. 5.3), for the
  superconformal scalar correlator.
* **Chester-Dempsey-Pufu**, arXiv:2111.07989 and 2312.12576: localization plus numerical bootstrap give rigorous upper bounds on the dimension of the lowest
  unprotected scalar of N=4 SYM at finite N. This is a "localization + bootstrap upper bound on a gap", but the operator is the lowest
  unprotected scalar (double-trace-like at strong coupling), not the string state.
* **Caron-Huot-Coronado-Trinh-Zahraee**, 2207.01615, and **Caron-Huot-Coronado-Zahraee**, 2412.00249: these input the integrability spectrum and bound OPE
  coefficients (two-sided, with localization), the opposite direction.
* **Alday-Hansen-Silva**, 2204.07542, **Alday-Hansen** 2306.12786, and **Chester-Hansen-Zhong**, 2412.08689 (ABJM): dispersive sum rules +
  localization + a worldsheet ansatz *determine* the first-massive-level dimensions to O(1/nu) and match integrability. This is far sharper than a
  factor-1.84 upper bound, though it rests on extra assumptions.
* **Binder-Chester-Jerdee-Pufu**, 2011.05728: N=6 numerical bootstrap with localization inputs (islands; ABJ higher-spin, string, M-theory regimes).
  It gives no string-gap bound.
* **Agmon-Chester-Pufu** 1711.07343 and 1907.13222, and **Chester-Dempsey-Pufu** 2412.14094: N=8 M-theory regime, R^4 / D^4R^4 from bootstrap and
  localization. No string scale.
* **Guerrieri-Penedones-Vieira**, PRL 127 (2021) 081601, 2102.02847: a 10d S-matrix bootstrap *lower* bound on the R^4 coefficient in Planck units,
  with no gap. Orthogonal.
* **Kaplan-Kundu**, 2008.05477: a flat-space causality upper bound on the mass of the lightest higher-spin state in terms of non-gravitational data.
  Same spirit ("upper bound on the string scale"), different input.
* The paper's own Einstein-rigidity theorem already contains a "read in reverse" bound on the string scale from tau. For ABJM tau = 0, so that
  theorem says nothing there, and the quartic inversion would be the first ABJM-specific gap statement in the paper.
* No 2023-2026 paper was found that states an *upper bound on the string/higher-spin gap of ABJM* from <TTTT> sum rules plus localization. The search was
  at web-search level: CHMRS 2106.10274 was checked in full text; the others from abstracts or the relevant sections.

How our statement would differ: (i) <TTTT>, 4d graviton helicity sum rules, rather than the superconformal scalar <SSSS>; (ii) AdS4 with the
infrared log and an explicit constant; (iii) the bound uses no SUSY. The **input** does, because the localization R^4 coefficient lives in <SSSS> and
is transferred by N=6 Ward identities and the flat limit, so "no SUSY used" is true of the bound only.

## 6. Verdict

"ABJM must have a stringy spin-4 operator below Delta_max ≈ 7.8 lambda^{1/4}, within a factor Y ≈ 1.84 of the true value":

* **Correct?** Yes, as an implication. The arithmetic and the inequality directions check. The λ^{1/4} *law* is not supported: it needs (Pos) for
  Omega^Q at every Delta_max(lambda), and the paper examines it only at 300 and 1000, and assumes it even there. The defensible form is the fixed-gap one:
  *for lambda < 2.8e8, twist < 1000*, which is within a factor 1.85 of the truth only at lambda ≈ 2.8e8.
* **Non-vacuous?** Formally yes. It is a CFT-side existence statement. Physically it is weaker than the known answer (QSC/integrability and CHZ),
  and it is *the same information* as "2.6% of the bound", with Y = (1/0.026)^{1/6}. The sixth root compresses a factor 39 of slack into 1.84.
  A referee will see this immediately.
* **Novel?** Only marginally. CHMRS 2021 made the analogous juxtaposition for N=4 with a sharper constant (Y ≈ 1.29). The ABJM/<TTTT>/AdS4 version is new
  as far as the search shows, but it is an incremental corollary.
* **PRL-headline-worthy?** **No.** The conclusion rests on the paper's weakest input ((Pos) for Omega^Q, grade [A]) plus L7(b) [A]. It reproduces
  known physics more weakly than existing CFT-based determinations, and its Y degrades as (ln lambda)^{1/6}. Recommended use: one sentence in the
  ABJM paragraph, replacing or alongside "2.6%".

Draft sentence (if used, not as a headline):
> Read in reverse, the quartic bound and the localization value of g4 imply that for 1 ≪ λ ≲ 3×10^8 ABJM theory must contain a single trace of spin
> four or more in T×T with twist below 1000, at most 1.9 times the first massive string level at λ ≈ 3×10^8. This is conditional on positivity of the
> quartic functional. CHMRS made the analogous comparison for N=4 SYM [2106.10274].

Draft 3-sentence headline paragraph (provided as requested; **not recommended**):
> Supersymmetric localization fixes the R^4 coefficient of ABJM theory in AdS units, g4/(8πG L^6) = ζ(3)/[32(2π²λ)^{3/2}]. Combined with our bound, and
> assuming positivity of the quartic functional at the relevant gap, this forces ABJM to contain a single-trace operator of spin ≥ 4 in T×T with twist
> below Δ_max ≈ 7.8 λ^{1/4} (Y = [(17.24 ln(2Λ/3) − 19.28)/2ζ(3)]^{1/6} ≈ 1.8-1.9 for λ ~ 10^7-10^9), derived from CFT sum rules and
> localization alone. The first massive string level of AdS4×CP3 lies at twist 4.22 λ^{1/4}, a factor 1.84 below.
