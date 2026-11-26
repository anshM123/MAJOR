# Referee check of Theorem RP (ALL_GAP_NOTE.md) and of the mixed-profile strengthening (MIXED_PROFILE_THEOREM.md)

26 Nov 2026. This is an independent review. No other file was edited and nothing was committed.

The notes under review are:
* `ALL_GAP_NOTE.md`: the upload, also copied at `code/uniform_gap/rational_profile/`.
* `MIXED_PROFILE_THEOREM.md`: in the session scratchpad `mixed/`.

The paper inputs I read:
* `prl/sec_fullpos.tex`: fp-U, fp-target, fp-H, fp-single, fp-V, fp-LR.
* `prl/supplement.tex`: the definitions, eq:abel, sec:HK, sec:c1, item (1), sec:light (LS rows), and the allgap Remark.
* `prl/data/kernels.json` and `prl/data/nf02_psi.json`.
* `K_B2A4.md` and `LEMMAS_CLOSE.md` §3.

## Scripts (python3, mpmath/sympy; logs alongside)

| script | what it checks | result |
|---|---|---|
| `r1_symbolic.py` | product rule for D_b on z^k(1-z)^M (symbolic k, M, b); D_b z^k against the coefficients of fp-single; R_+(z)=R_-(-z) from kernels.json; P(z+i)/R_-(z); Fstar(z+i)/Fstar(z-i)=R_-/R_+; raw light identities R_-Q(ν+i)+R_+Q(ν-i)=P_K (K_B2A4 §4.5) and the Einstein one; uniqueness of Q; g_{k,3} built independently by a Pfaff/hypergeometric series; Q=-100g_2-36g_3; [ν^6]g_3=-1/36; the §6 J=0 and J=2 formulas | 16/16 pass |
| `r2_special.py` | Mehler–Fock for 1/(x+A) (DLMF 14.20.14 form) to 1e-30; the transform ∫P/(x+A)=πP(A)/cosh(πν); the resolvent ∫Q_1P=1/(ν²+9/4) to 1e-17; Φ>0 with Φ~2.22x^{23/2}e^{-πx}; B1=B2 in the strip; B real-even with a double zero at 0 and simple poles at ±i/2 and no others on a grid; \|B\|~x^{19/2}e^{-πx} uniformly in \|y\|≤1; the Mehler–Dirichlet bounds | all pass |
| `r3_endtoend.py` (M=3) | end-to-end at small M, at 60 digits. f by partial fractions; (2/Λ)∫E P=f; F_k=(2/Λ)∫E g_k (k≤4); -(2/Λ^7)∫EQ=36; I(z-form)=(2/Λ)∫ED; χ=Bq built from E, with its leading-Regge image =f(η); ∫χP_K=∫EQ; ∫χ·Ein=∫ED; the **paper's own moment route** (Z_k from χ, then L_k, then D_k²L_k) reproduces F_k | all agree to ≥12 digits |
| `r4_heavy.py` | at M=998, 800 digits, J up to 1500 (beyond the earlier J≤64). The lemma identity L_b[z^k(1-z)^{m0}] = U_k·E(1-B0)^X, compared against the direct alternating sum, including m0=M+1, a=M and fractional twists; the full V_b for F_strict and F at a∈{998, 998.37, 1000, 2500, 1e5, 1e9} | identity holds to 1e-480 or better; all strict values >0; non-strict values are 0 exactly for n≤1 and >0 otherwise |
| `r5_lightbound.py` | the integral I by quadrature against the analytic bound (ratio 0.62 at M=998); the τ bound against the 280 claim (ratio 0.79); the crossover with τ0(1000) | pass. **Crossover at Δ_g ≈ 1.79e5** |
| `r6_mixed.py` | the mixed profile at M=998, J≤1200, 800 digits. min H/(w³L(F)); tolerance c* of b=4 positivity to the k=3 heavy normalisation; KΛ^6 split; I; the 420 claim | min ratio 0.337 (claim ≥0.04); c*=1.49; KΛ^6=5.54e6 (k=2 part −3.58e4, k=3 part +5.58e6); τ/claim=0.74; **crossover Δ_g≈1.83e4** |
| `r7_k3light.py` | g_{k,4}=D_k²q_k (fp-cdh d=1) and 4g_{2,4}=ν⁴−(35/2)ν²+569/16 | pass |
| `r8_flatJ0.py` | the heavy formula continued to J=0 against the paper's flat light-scalar F_D^(0) for NF02 | 20.771893 and 11.317027 reproduced exactly |

## Findings: ALL_GAP_NOTE (Theorem RP and the 280 bound)

**Correct (confirmed):**
* **Positive-basis lemma (§2.2).** The shifted-ladder algebra is right: n'=n−k, N'=N−k, n'+s'=n+s+k and N'+s'=N+s+k.
* **Nonnegativity of p_j.** The p_j are nonnegative for real a≥M and m0≤M+1. The smallest factor of (a)_{m0−j} is a−M≥0. The factors of (N')_{m0} are positive because N'−m0+1≥1. Zeros of (n')_j truncate j at n'.
* **Distribution and moments.** Vandermonde and the factorial moments hold as polynomial identities. The Beta(n'+s', a) moments are correct, so E(1−B0)^X>0.
* **Dressing (§2.3).**
  * The product rule is exact.
  * The off-diagonal coefficients are ≥0. M(M−9+4b)/2>0 for both b at M≥998.
  * The diagonal bound −C<T/2 for k≤n is valid. It uses only k≤J, a≥M and J+a≥998.
  * At k=n+1 the diagonal ladder value is zero. For k>n+1 every term vanishes. Termwise summation is exact, because only k≤n+1 enters and each P_{k,M} starts at z^k.
* **Strict version (§6).** It is correct: n=0 gives ε/(2T).
* **Realisation (§1, §2.1).** It is correct, and end-to-end numerics confirm the de-shift, the leading-Regge image and the paper's moment route.
* **Light normalisation (§3).** It is consistent with the paper's conventions. With ψ2(p)=E(Λp), −(2/Λ^7)∫Eν^6 is −2∫ψ2p^6 and 16Λ∫ED is 16L_A. The inequality is homogeneous in E, so K=36/Λ^6 is simply the unnormalised value.
* **Constants.** 57600, 16 and 8 combine to give 12800 and 460800 as claimed.
* **No linear t_4 term.** c_1=0 is pointwise in ν (sec:c1), so it holds for any weight.
* **Robustness of K.** Because F_0=F_1=F_2=0, the moments ∫Eν^{2j} vanish for j≤2, and K depends only on the flat-calibrated ν^6 coefficient. The lower (AdS, H17-type) coefficients of Q and any t_4²/α_T² asymmetry at O(Λ^{-2}) drop out exactly.

**Problems:**
1. **[presentational, but must be fixed before stating] The bound adds nothing for Δ_g < 1.79·10^5.**
   * The paper already gives τ ≤ τ0(1000) = 2.89e-8 at every Δ_g ≥ 1000. This is Theorem FP at 1000 combined with the monotonicity Remark of sec:allgap.A.
   * 280√log/Δ_g² is weaker than that bound below Δ_g ≈ 1.79e5. At 1000 the new bound is 6.8e-4.
   * The note says only that the old bound is "much stronger at 1000".
2. **[fixable: grade and an explicit input] Interchange without the cutoff θ.**
   * The paper says that without θ "the interchange cannot be proved under J0≤1". Its fallback (sec:defs, LEMMAS §3.1 step 4) needs only a finite γ, because cut-off weights decay super-Gaussianly.
   * The new weight matches the heavy normalisation exactly: χ/N_com ~ ν^{-10}, and on the edges χR_-/(E_7 s) ~ ν^{-10}. The margin is therefore only polynomial.
   * Admissibility thus needs K2ν in its sharp form (γ=0 in position normalisation). That form rests on a numerical sign check of the dual bases (β≤16). A weight-side check of [D] grade is not the paper's [D] for its own weights.
   * State the sharp K2ν, together with K1′, K2, K3, LRC-a and LRC*, as retained inputs of the new functional.
   * The analytic structure (strip analyticity, the simple poles at ±i/2 cancelled by R_∓, the double zero at 0) is confirmed.
3. **[presentational] Size of the realised weight.**
   * The asymptotic O(ν^{-1}e^{-πν}) decay of χ, with E,E′=O(νe^{-πν}), is true, but its constants are astronomically large: Σ\|C_j\| ~ 10^{5700} at M=998.
   * There is a long pre-asymptotic regime in which E grows like ν^{2m−1}e^{-πν}; this was seen at M=3. At M=998, E peaks near ν≈2m/π≈640≈0.64Λ.
   * This is harmless for an exact existence theorem at leading order in 1/c_T. It does mean that nothing quantitative about how large c_T must be can be inherited.
   * Wording: "E_7 … grows as \|ν\|^11 e^{-π\|ν\|}" should say "behaves as".
4. **[presentational] The §6 light-spin check uses the heavy formula, whereas the paper treats light exchanges by their Regge polynomials (LS4, (E′)).**
   * Supporting evidence: the heavy formula continued to J=0 reproduces the paper's flat F_D^(0) for NF02 exactly (r8).
   * For F_strict, F_0=0. The J=0 value therefore depends only on F_1 (the flat top coefficient), and the μ² AdS shift multiplies zero, so the result is robust given (E′).
   * For the non-strict F, only the ν^6 moment survives at spin 2.
   * The conclusions are consistent with the paper's light treatment at grade [D-sketch, given (E′)]. They are not independent.

**Grades I would assign:**
* **Theorem RP (heavy positivity, both versions):** [P given Eq. fp-target], i.e. the same conditionality as Theorem FP, and analytic.
* **Realisation:** [P]. Mehler–Fock, the inverse and the 2i-periodicity are exact; DV end-to-end.
* **The no-cutoff interchange for this weight:** [D, new], conditional on the sharp K2ν (γ=0), K1′–K3, LRC-a and LRC*.
* **The 280 bound:** [D] at leading order in 1/c_T, with item (1)'s inputs plus that interchange. It must be stated together with the crossover at 1.79e5.
* **Low-spin light check:** [D-sketch given (E′)].

## Findings: MIXED_PROFILE_THEOREM (the 420 bound)

**Correct (confirmed):**
* **Positivity argument.**
  * The diagonal bound −C ≤ x²/4+(9/2)x+15/2 < 0.26x² < 0.26T, so 1+C/T > 37/50.
  * Domination L(G) ≤ A·L(F) holds, because F−G/A = εz(1−z)^M R is in the positive cone.
  * A/T ≤ 7/10 because T ≥ a² ≥ M².
  * Together these give the 1/25 margin, and n=0 gives w³ε/(2T).
  * The proof uses only \|σ_b\|=1, so it is robust to the σ assignment.
* **Numerics at M=998, J≤1200.** min H/(w³L(F)) = 0.337.
* **(b) k=3 light normalisation.**
  * g_{k,4}=D_k²q_k, and 4g_{2,4} is exactly K_B2A4's k=3 light polynomial. So K_{k=3} = 2∫ψ3[p⁴−(35/2)p²/Λ²+(569/16)/Λ⁴] = 8G_2/Λ^6, in the paper's convention (K_0 = −2∫ψ2p^6 + 2∫ψ3p^4).
  * Because G_0=G_1=0, only the flat top k=3 coefficient enters.
  * The Einstein term at k=3 is zero in the paper's sum rule (L_A contains ψ2 only). The k=3 functional is a commutator structure, so the diagonal Einstein one-point matrix does not contribute.
  * The k=2 part is −3.6e4 (≈100−36M), negative and O(M) against the k=3 +5.6e6.
* **(d) Averaged-R realisation.** The pushforward in A of ω(t)C_j(t)dt is a smooth density with compact support in A>1. Integrating by parts with ∂_A((A²−1)∂_A), whose eigenvalue on P_ν is −(ν²+1/4), and applying the Mehler bound gives O(ν^{-K}e^{-πν}). Both f_2 and f_3 are proper. Σρ_j(t) < 1/100.
* **(e) The 420 chain.** The harmonic differences are ≤ log(M+2) ≤ log Δ_g. The step (M+1)^{-3}M^{-2} ≤ (1+1/M)²a_g^{-5} is valid. The squared coefficient is 170103 < 176400, and the numerical τ is 0.74 of the claim.

**Problems:**
1. **[fixable, but must be an explicit assumption] (a) Admissibility of a non-compact k=3 weight.**
   * The paper's k=3 piece is "unsubtracted, not de-shifted". Its interchange is "plain Fubini" (sec:HK conclusion) with compact ψ3 on [0,1].
   * K2ν in LEMMAS_CLOSE is proved only for the j′=0 kernels (k=2, 4). The growth of 1/N̂_3 is not given in the paper.
   * A non-compact E_3 therefore needs two things that are new: a j′=1 kernel modulus bound, and ∫\|E_3\|/\|N̂_3\| < ∞, again with an exactly matched e^{πν} rate.
   * The note's C_A, C_B computation is supportive, but the note itself says it does not fix the conventions. Keep this as an explicit [A] or [D-sketch] input.
2. **[fixable: grade] The k=3 normalisation is now load-bearing twice.** It enters both as the sign-critical b=4 heavy term and as the leading light coefficient.
   * The paper grades the k=3 prefactor of fp-target as [D]; only the k=2 part is Einstein-calibrated in exact arithmetic.
   * K_B2A4 §7.3 lists the k=3 steps (j′=1 structure, N̂_3 normalisation) as lacking an end-to-end check.
   * Tolerance: the analytic proof survives a relative factor c ≤ 1.057 on the k=3 heavy term. The sampled states survive c ≤ 1.49.
   * The bound scales like c^{-1/2}. A sign error in the k=3 light coefficient would destroy the bound entirely.
   * In NF02 these were small corrections. Here they are the result.
3. **[presentational, must be stated] Crossover.** The 420 bound is below τ0(1000)=2.89e-8 only for Δ_g ≳ 1.83·10^4. The note's §4 sentence "the fixed-gap certificate … is not transferred … by monotonicity" is misleading. The positivity does not transfer, but the τ0(1000) *bound* does transfer to every Δ_g ≥ 1000 (sec:allgap.A Remark).
4. **[minor] Low-spin check.** Now F_2=1 and F_3≈−M, so the spin-2 sign relies on G_2 = 0.7M² beating −9M/2. Lower coefficients with O(1) AdS shifts multiply F_2=1 only, so the check is robust given (E′). It is still [D-sketch].

**Grades I would assign (mixed):**
* **Heavy positivity:** [P given fp-target], where the fp-target k=3 prefactor is [D].
* **The 420 bound:** [D-sketch]. It is conditional on the explicit k=3 admissibility input and on the k=3 normalisation chain (K_B2A4 k=3 [D] with unchecked steps), in addition to the inputs of the RP bound.
