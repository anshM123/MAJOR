# Independent review of the "uniform-gap continuation" research package

*Reviewed 26 Nov 2026. I read the package's `RESEARCH_NOTE.md` (Secs. 1–6) and the TeX (`prl/supplement.tex`, `prl/sec_fullpos.tex`),
then re-implemented the checks from the written statements. The package's scripts were read only to learn two definitions: what
its "normalized" lower bounds are normalised by, and which boxes the certificates use. None of them was run or imported.
The xcheck_ladder implementation of Eq. (fp-H) and of N_s (`code/xcheck_ladder/xl_common.py`, `xl_norm.py`) was imported read-only.
Nothing outside this directory was modified, and I made no commits. Commits `2c4c08b`, `ab9261e` ("Work in progress: review
outputs") and `695687c`, which contain in-progress files from this directory, were made by some other process.*

**Grades used here.** "confirmed" means an exact or symbolic check, or a derivation re-done by hand. "confirmed numerically"
means a high-precision numerical test, which is not a proof. "not checked" means what it says. "problem" means a statement
that is false or does not follow.

## Verdict

**I found no false statement.** All the numbers I could recompute agree, and none of the inequalities I tested is violated.

**Status of each part:**
* **Exact algebra (Secs. 1 and 2).** Confirmed: the Gamma/rational identities, the light polynomial identity and the four
  exact K(t) coefficients.
* **Analytic arguments (Sec. 1 inverse, Sec. 2 log bound).** Re-derived. They are correct and short.
* **Angular action formula (Sec. 3).** Confirmed numerically: the closed form equals the paper's Eq. (fp-H) applied to the new
  weight, with relative normalisation exactly 1.
* **Positivity statements (Sec. 3).** Consistent with dense sweeps in ball arithmetic. I did not re-run their Bernstein
  certificates.
* **Cutoff and residual estimates (Secs. 4–5).** Every elementary constant I tested holds, most with large room. The long
  chains of 2^N bounds were spot-checked, not re-derived line by line.
* **T0 transfer (Sec. 5).** Logically sound given T0_CERT. The quoted margins exist in T0_CERT.

## 1. Claims table

| # | claim (note section) | my check (script → log) | result |
|---|---|---|---|
| 1 | F(z+i)/F(z−i) = R₋/R₊ (§1) | Gamma/cosh recurrences + SymPy (exact zero); mpmath with the actual Γ at 5 complex points, 4·10⁻⁴¹ (`rv1_symbolic.py`) | **confirmed** |
| 2 | P(z+i)/R₋(z) = −z(2z−7i)²(2z−3i)²(2z+3i)(2z+7i)²/2048 (§1) | SymPy, exact zero | **confirmed** |
| 3 | B = F(z+i)/R₋ = F(z−i)/R₊ is even and real on ℝ (§1) | follows from 1; numerically 5·10⁻⁴¹ | **confirmed** |
| 4 | R₋(x−i)χ(x−i)+R₊(x+i)χ(x+i) = f(x) for χ = Bq, q the sech convolution (§1) | built q and χ in mpmath at Λ² = 900 and Λ² = 997002 (`rv2_deshift_check.py`). (A) limit from inside the strip, y = 1−ε → 1, Richardson; no PV formula used: agrees with f to 10⁻¹⁴–10⁻¹⁸ relative (ε-extrapolation-limited). (B) stated PV form h: 10⁻²³–10⁻²⁹. Beyond Λ, where f = 0, the two terms cancel to 10⁻¹⁹ of their size | **confirmed** (analytic argument re-derived; numerics agree) |
| 5 | same equation and normalisation as the paper (Eq. eq:abel; E = 2Re[R₋χ](ν−i) = ψ₂(ν/Λ)) | by hand: Fourier transform of the identity gives Eq. eq:abel exactly, with r±Y = (R₋±R₊)χ and ψ₂(ν/Λ) including 1/k0; R₊(x) = conj R₋(x) on ℝ (SymPy), so the two-term sum is 2Re[…] for real-analytic χ. Numerically, the **published** Δg=1000 weight gives 2Re[R₋χ](μ−i) = ψ₂(μ/Λ) to 1.85·10⁻⁷ of max ψ₂ (`rv3_paper_convention.py`) | **confirmed: identical conventions** (no factor or sign difference) |
| 6 | four exact K(t) coefficients, K(1/997002) = 0.999967523127696 (§2) | exact rationals from nf02_psi.json (`rv4_light.py`): identical. Paper/K_B2A4: 0.99996752 / 0.999967523128. Λ²·(t-coefficient) = −32.379 (paper "−32.4/Λ²"); k = 3 part −82.0165 (K_B2A4 −82.0165); K0−1 = 3.16·10⁻¹¹ | **confirmed** |
| 7 | raw k = 2 light function = R₋Q(z+i)+R₊Q(z−i) (§2) | SymPy against K_B2A4 §4.5's P = −(4ν²+1)N₁₀/(8(4ν²+9)(4ν²+49)): exact zero. Einstein: R₋D(z+i)+R₊D(z−i) equals the Lemma-R function iff D = 1/(z²+9/4) | **confirmed** |
| 8 | −511.773 < g ≤ −188.896, ∫g = −332.9506 (§2) | exact ∫g; dense sampling gives −511.77285 / −188.89611 | **confirmed** (range numerically) |
| 9 | 16L_A ≤ 3023 log(2Λ/3) − 5307 (§2) | 16L_A = 14339.589 at Λ = 998.5 (paper 14339.59) vs bound 14344.88 (slack 5.3); slack 12.0 / 18.8 / 24.1 / 27.4 / 30.6 at Λ = 1500 / 3000 / 10⁴ / 10⁵ / 10⁷. The error estimate Aε²/2 + (π/2)ε‖g‖ is correct | **confirmed** (derivation + numerics) |
| 10 | conditional rigidity √(t₄²+α_T²) ≤ (240/Λ⁴)√(E/K) (§2) | same form as the paper's τ₀; gives 2.8919·10⁻⁸ at 1000 (paper 2.8913·10⁻⁸, as an upper bound should) | consistent |
| 11 | recurrence for g[k,p] (§3) | against mpmath Taylor coefficients of (1−z)^−p P_{−1/2+iν}((1+z)/(1−z)), k ≤ 8, p = 3, 4: 2·10⁻³⁹ (`rv6_compare_fpH.py`) | **confirmed** |
| 12 | closed form = H_ang,b/w³ of Eq. (fp-H) (§3) | fed xcheck_ladder's `H_ang` with moments Z_k, Y_k computed from the **new** χ (double-precision GL measure, `rv5_build_chi.py`; spot-checked against mpmath to ≤ 3·10⁻¹³ (`rv5b_spot.py`)), with and without θ(Λ+2000,10). Δg = 1000: 18 states (J = 0–64, both b, twist edge and other masses), ≤ 4.3·10⁻¹² relative. Λ = 30: 10 states, ≤ 5·10⁻¹¹. Δg = 3000 (separate χ table): 11 states, ≤ 4.1·10⁻¹² (`out/rv6_Dg3000.log`). Normalisation factor exactly 1. Also Z₀/3 = 20.7718930133 = 2∫ψ₂ | **confirmed numerically** |
| 13 | J ≤ 64: every physical channel positive for all gaps ≥ 1000 (94 exact certificates; min normalised bound 2.72·10⁻⁸) (§3) | ball-arithmetic sweep of the closed form (`rv7_sweep.py low`): 94 channels × 12 gaps (1000, 1001, 1100, 1500, 2000, 3000, 5000, 10⁴, 3·10⁴, 10⁵, 10⁶, ∞) × 410 masses from the twist edge to m = ∞, plus golden-section refinement. **All positive.** Min S = H/w³ = 4.67·10⁻⁴ (J = 64, b = 4, Δg = ∞, twist edge); 5.94·10⁻⁴ at Δg = 1000. Normalised by N_s at the twist edge: min H/N_s = 0.20087 (1000), 0.20026 (3000), 0.200007 (10⁵), all at J = 5, +− (`rv14_HNs.py`). Full [0,1]² (q,u) rectangle: J ≤ 32 positive everywhere; non-positive points only for J ≥ 44 at u ≥ 0.99, i.e. non-physical m < m_gap, as the note says (`rv15_rectangle.py`) | **confirmed numerically**; the exact certificates themselves **not re-run** |
| 14 | J ≥ 65, b̂ ≤ 40: positive, normalised lower bound 0.01519546 (§3) | sweep (`rv7_sweep.py imp`, `imp2`): J ∈ {65…140, 200, 300, 500, 1000, 2000, 5000, 10⁵}, 6 gaps, b̂ step 0.25 or 0.5. **All positive.** Min S = 0.0159064 (J = 10⁵, +−, Δg = ∞, b̂ = 40); 0.015907 at J = 1000. Above their certified lower bound 0.0151955, ratio 0.955 | **confirmed numerically** (consistent); certificate **not re-run** |
| 15 | F, |F(x±i)|, |B| envelopes (§4) | `rv8_bounds.py`, x ∈ [0, 10⁷]: F/(1+x)^{23/2} ∈ [0.304, 14.6] (claim [2⁻¹¹, 2¹³]); |F(x±i)|/(1+x)^{23/2} ≤ 4.52 (claim 2¹⁴); |B|/(1+x)^{19/2} ≤ 0.449 (claim 512); limits at 10¹²: 1.1107 (= π/2^{3/2}), 0.0694 | **confirmed numerically**, large room |
| 16 | Binet bounds for cosh(πy)|Γ(a+iy)|² (§4) | re-derived from Binet's |R(z)| ≤ 1/(12 Re z) and 2y·arctan(y/a) = π|y| − 2|y|arctan(a/|y|); numerically for a = 3/4, 5/4, 7/4, room ≥ 2.2 | **confirmed** |
| 17 | shifted factor ≤ 16(1+x)^{3/2}; (1+s)¹² e^{−s/2} ≤ 24¹²; 26¹³ < 2⁶²; csch(πs/2) ≤ 3e^{−1.5s} (§4–5) | max 1.11; all true | **confirmed** |
| 18 | |χR±| ≤ 2⁹⁰ e^{−(x−L)₊}, |χ| ≤ 2⁹⁰(1+x)⁻² e^{−(x−L)₊} (§4) | on the χ tables (`rv12_chi_envelope.py`): 2⁸·⁰ and 2⁴·³ (Λ = 998.5). Actual decay rate beyond Λ is 1.56 ≈ π/2 > 1 | **confirmed numerically** |
| 19 | (k!)² t^k |g[k,3]((ν±i)²)| ≤ (1+y²)^k (k ≤ 64); t^k|g[k,3]| ≤ 2e^{86y}1600^{−k} (§4) | flint/acb, k ≤ 65 and k ≤ 250, 18 values of ν up to 10⁶, Λ² = 997002 and 10⁸: both hold. The first is **tight** at k = 1 for ν ≫ Λ (log-ratio −5·10⁻⁶, exactly true at k = 1 by hand); the second has room e^{−0.71} | **confirmed numerically**; the analytic "majorants sum to < .01" step was not re-derived |
| 20 | envelopes (2+y²)e^{130(1+y)} and 16(1+y)²e^{140(1+y)}; the 2⁻⁵⁰⁰ light and angular cutoff errors (§4) | plausibility only: the exponents are consistent with items 18–19 | **not checked in detail** |
| 21 | r = (Re θ(x−i)−1)f − 2Im θ(x−i)G, r(0) = 0 (§5) | re-derived (E_θ = 2F Re[(W/2+ih)θ]); numeric identity check with a deliberately close cutoff (C = L+30) agrees to 10 digits (`rv9_residual.py`) | **confirmed** |
| 22 | L(‖r‖₁+|r′(0+)|+‖r″‖₁) < 2⁻¹⁰⁰⁰ (§5) | computed r directly at Λ = 998.5 with the stated cutoff. max |r| ≈ 10⁻¹³⁴⁶, near x ≈ Λ+1900 (2⁻¹⁰⁰⁰ ≈ 10⁻³⁰¹). The dominant terms are far smaller than claimed. r″ was not computed separately; its prefactor is (C−x)/50-type, so the same order | **confirmed numerically** (‖r‖ part); the derivative bounds and the 2¹⁸¹/2²⁰⁰ chain **not re-derived** |
| 23 | Sec. 5 Sobolev-type constants: ψ₂^{(j)} = (1−p)^{5−j}Q_j, ‖Q_j‖₁ < 2⁴⁰; |(log F)^{(j)}| < 100; |F^{(j)}|/F ≤ 2²¹ | exact division (SymPy): ‖Q_j‖₁ ≤ 2.2·10⁹. Numerically, max |(log F)^{(j)}| = 5.4 / 19.2 / 49.3 and |F^{(j)}|/F ≤ 164 (`rv13_sec5_constants.py`) | **confirmed** |
| 24 | M₁ > 9·10⁻⁵ on β ∈ [7/4, 12]; M₁ ≥ (√2/π)K₀/(10β) for β ≥ 12; β·ε₁ ≤ 1 (§5) | M₁ from T0_CERT Lemma 2's definition, mpmath: min 1.146·10⁻⁴ (at β = 12); ratio ≥ 8.7 on β ∈ [12, 100]; max β·ε₁ = 0.0145 | **confirmed numerically** |
| 25 | relative Mehler bound |∫rP| ≤ (√2/π)K₀[B/β + 8B log2/β² + Aε₁] (TAIL_CERT R.2) | toy residuals with r(0) = 0 (`rv11_mehler_toy.py`): ratio LHS/bound ≤ 0.39 on β ∈ [12, 60]. The two integrations by parts giving min(A, B/t²) are correct | **confirmed numerically** (sanity); lemma not re-derived |
| 26 | T0_CERT margins: Part F slack ≥ 6·10⁻⁷, Part K ≥ 4.58/Λ², Part L ≥ 0.81·S_min M₁/Λ² with S_min > 200 (§5) | T0_CERT §3 table and §4.2. F: 6·10⁻⁷ (Δg ≥ 300 run; 9·10⁻⁷ in the ≥ 1000 run). K: 4.58 (≥ 300 run; 4.92 at ≥ 1000). L: max ratio 0.1905, so 0.81; S_min = 241.89. All for c = 0.199 and valid for every Λ ≥ Λ_min | **confirmed** (they exist and are quoted correctly) |
| 27 | transfer: H_hyb ≥ 0.199 N_s for every gap ≥ 1000 for the realised (cut) family (§5) | the logic is sound (see §3 below) | **confirmed (logic)**, given T0_CERT |

No "problem found" rows.

## 2. Details and caveats (none is an error)

* **Convention (item 5).** The note's identity, the paper's shifted density E(μ) = 2Re[χ(μ−i)R₋(μ−i)], Eq. eq:abel, and the
  published data file's `identity` field (W[χ](t) = ∫₀^∞ψ cos νt) are the same equation. Same target, ψ₂(ν/Λ) with 1/k0; same
  factor 2; same signs. The new family solves it exactly. The published weight solves it to ~2·10⁻⁷ of max ψ₂.
* **The new family is a different functional.** It differs from the published continuum weight at Δg = 1000 by a null-space
  element. At the twist edge, H/w³ differs by 3·10⁻¹⁰ (J = 4) to 1.4·10⁻⁵ (J = 32, +−) relative (`rv10_vs_published.py`). So, as
  the note says, Theorem FP's Δg = 1000 certificates do not transfer automatically, and the note's results are not about the
  published functional.
* **What the J ≤ 64 / b̂ ≤ 40 positivity is about.** It is about the ideal (uncut) closed form plus the Sec. 4 cutoff transfer,
  and only given Eq. (fp-target)/(fp-H) as operator identities, which are [D] in the paper. My checks validate the closed form
  against Eq. (fp-H). They do not test (fp-H) against the CFT.
* **Normalisations of the quoted minima.** The note's 2.72·10⁻⁸ and 0.01519 are lower bounds on S = H_ang/w³. For J ≤ 32 the
  package divides by the denominator bound (1+C/D)^{M+3}. For J ≥ 33 it also divides by (998+2J+1)^M, which is why 2.72·10⁻⁸ is
  ~10⁴ below the true minimum. My minima are actual values of S, so they must lie above the bounds, and they do.
* **Tight constant.** In item 19, the first inequality is tight (relative 5·10⁻⁶) at k = 1 for ν ≫ Λ. It holds there, as checked
  by hand, but a referee should check that the general-k proof treats ν → ∞ uniformly. The note's one-line justification
  ("lower-order majorants sum to < .01") does not show this explicitly.
* **Analyticity for H_K.** I checked only that χ = Bq is analytic in the open strip |Im ν| < 1 (B analytic for |Im| < 3/2,
  double zero at 0, which cancels the 1/ν pole of R±), with L¹ lines. I did not re-check the paper's LRC*/corner-count argument
  for this weight. The note retains it as an assumption ("same strip analyticity mechanism").
* **Numerical precision notes, for anyone re-running.** Two traps, both in my own code, both fixed. (i) The monomial form of
  ψ₂ loses all accuracy near p = 1 in double precision. Use the factored (1−p)⁵pX(p)/k0 form. (ii) mpmath `quad` over long
  intervals needs breakpoints about every ½ unit for the sech convolution. Without them, spot checks can be off by 10⁻²
  relative.

## 3. Sec. 5 transfer argument

The steps are:
1. The realised k = 2 leading-Regge image equals the target image (2/Λ)∫ψ₂P plus (2/Λ)∫₀^∞ r P_{−1/2+iν}(cosh β)dν.
2. T0 is a statement about the target action ("T0 is proved for the target (de-shifted) weights at every Δg ≥ 300").
3. H − cN_s with the realised component therefore moves by at most (1+c)·2/Λ·|∫rP|.
4. **Parts F and K and β ≤ 12 of Part L.** |P| ≤ P_{−1/2} ≤ 1 bounds this by 2.4‖r‖₁/Λ. That is ~10⁻¹³⁴⁰, against margins of
   ≥ 6·10⁻⁷, ≥ 4.58/Λ² and ≥ 0.81·242·9·10⁻⁵/Λ².
5. **β ≥ 12 in Part L.** The Mehler relative bound gives |∫rP|/[(√2/π)K₀/(10β)] ≤ 10[B(1+8 log2/12) + A]. This is uniform in
   β, so the floor dominates. This is exactly the β ≲ 200 cap that T0_CERT §4.2 could not remove for the published weight; the new
   weight removes it because r(0) = 0 exactly and r″ ∈ L¹.
6. The note correctly uses c = 0.199, not 0.1999. The 0.1999 run has zero Part-F slack.

**Inherited caveats.** T0's OO sector is [certified given helicity lemma (A)], and T0 inherits CERT_CLOSE's flat certificate,
whose corner code is marked "unreviewed". The conclusion bounds H_hyb only, not H_ang − H_hyb. The note states this.

## 4. Sweep minima

All sweeps evaluate the closed form in Arb at 600–6000 bits, with ball radii ≤ 10⁻⁹⁵ relative. Each point is rigorous; the
sweeps are not.

| region | points | min S = H/w³ | where | normalised min |
|---|---|---|---|---|
| J ≤ 64, all physical masses, 12 gaps incl. ∞ | 94 × 12 × ~410 | **4.67·10⁻⁴** (5.94·10⁻⁴ at 1000) | J = 64, +−, Δg = ∞, twist edge | H/N_s ≥ **0.20087** (1000), 0.20026 (3000), 0.200007 (10⁵), at J = 5 +− twist edge; the paper's +0.200868 at 1000 is reproduced |
| J ≥ 65, b̂ ≤ 40 (claimed region) | ~25 000 | **0.0159064** | J = 10⁵, +−, Δg = ∞, b̂ = 40 | claimed lower bound 0.0151955 ✓ |
| J = 65–300, 40 < b̂ ≤ 400 (**not covered** by the note) | 44 320 | **1.60·10⁻⁵**, all positive | J = 300, ++, Δg = ∞, b̂ = 400 | H/N_s ≥ **0.9016** at Δg = 1000 (J = 300 +− twist edge, b̂ = 475); H − H_hyb ≈ 10⁻⁵ N_s (`rv16_HNs_out.py`) |
| twist edge J = 65–1000 (uncovered) | 330 | 1.28·10⁻⁷, all positive | J = 1000, +−, Δg = ∞, b̂ = 2001 | H/N_s = 0.918 at J = 1000 +−, Δg = 1000 |

**For planning.** In the uncovered region the action is numerically far from zero in N_s units (≥ 0.90 on every point
sampled). The finite-mass correction H − H_hyb is ~10⁻⁵ N_s. The obstruction is purely one of proof technique: cancellations
of size e^{b̂} in the ladder sum.

## 5. Recommendation

**Safe to put in a paper now:**
* **Proved** (exact algebra plus short, correct analytic arguments):
  * Sec. 1 (exact inverse of the de-shift equation, same convention as the paper);
  * the Sec. 2 identities and exact K(t) (K ≥ 0.99996 on all gaps ≥ 1000);
  * the Sec. 2 log bound E ≤ 3023 log(2Λ/3) − 5307 for the uncut family.
* **Certified** (their exact integer certificates; not re-run by me, but consistent with my sweeps and resting on a closed form
  that I confirmed numerically against Eq. (fp-H)):
  * the J ≤ 64 fixed-spin positivity at all gaps;
  * J ≥ 65, b̂ ≤ 40.

  Before printing, an independent re-run of `check_uniform.py`, `check_physical.py` and `impact_tensor.py` is still needed
  (someone else is doing this).
* **Proved modulo routine verification.** The Secs. 4–5 cutoff and residual estimates. The constants hold with large room, except
  the tight k = 1 bound of item 19. A referee should still read the 2^N chains line by line; I spot-checked them numerically.
* **Certified given T0_CERT.** The realised hybrid floor H_hyb ≥ 0.199 N_s for every gap ≥ 1000.
* **Grading note.** Everything about H_ang carries the paper's grade for Eq. (fp-target) ([D] as an operator identity).

**What remains:**
* Full finite-mass positivity for J ≥ 65, b̂ > 40 at all gaps. It is positive on every point I sampled, with margin ≥ 0.90 N_s,
  but it is not proved.
* An independent re-run of the exact certificates.
* A line-by-line referee pass of the Sec. 4/5 bound chains, including the general-k proof of item 19.
* The H_K/LRC* analyticity check for the new weight.

## Files

| script | purpose | output |
|---|---|---|
| `rv_common.py` | NF02 (exact; factored form asserted equal to psi*_exact), R±, F, B in mpmath | — |
| `rv1_symbolic.py` | items 1–3, 5 (R₊ = conj R₋), 7 | `out/rv1_symbolic.log` |
| `rv2_deshift_check.py` | item 4 (args: Λ², dps) | `out/rv2_L30.log`, `out/rv2_L998.log` |
| `rv3_paper_convention.py` | item 5 on the published weight | `out/rv3_paper_convention.log` |
| `rv4_light.py` | items 6, 8, 9 | `out/rv4_light.log` |
| `rv5_build_chi.py`, `rv5b_spot.py` | double-precision χ table of the new family; mpmath spot check | `out/chi_*.npz`, `out/rv5*.log` |
| `rv_closed.py` | the note's Sec. 3 closed form (exact rationals) | — |
| `rv6_compare_fpH.py` | items 11, 12 | `out/rv6_L30.log`, `out/rv6_L998.log`, `out/rv6_Dg3000.log` |
| `rv7_sweep.py` | items 13–14 and the uncovered-region probes (`low`, `imp`, `imp2`, `out`, `edge`) | `out/rv7_*.json/log` |
| `rv8_bounds.py` | items 15–17, 19, 24 | `out/rv8_bounds.log` |
| `rv9_residual.py` | items 21–22 | `out/rv9_residual.log` |
| `rv10_vs_published.py` | new vs published functional at 1000 | `out/rv10_vs_published.log` |
| `rv11_mehler_toy.py` | item 25 | `out/rv11_mehler_toy.log` |
| `rv12_chi_envelope.py` | item 18 | `out/rv12_chi_envelope.log` |
| `rv13_sec5_constants.py` | item 23 | `out/rv13_sec5_constants.log` |
| `rv14_HNs.py`, `rv16_HNs_out.py` | H/N_s normalised margins (N_s from xcheck_ladder/xl_norm) | `out/rv14_*`, `out/rv16_*` |
| `rv15_rectangle.py` | full-rectangle scan (non-physical failures for J ≥ 44) | `out/rv15_rectangle.log` |

On `rv7_imp2.log`: "max truncation rel 7·10⁻⁶" is a diagnostic. It compares the full ladder with the ladder minus its last 20
terms, for J ≤ 176. For J > 176, where the ladder is truncated at k = 176, the change is ≤ 10⁻¹¹².
