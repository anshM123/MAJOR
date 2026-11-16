# Referee report: `ads4/FLAT_GRAVITON_BOUND.md` + `ads4/code/flatbound/`

Referee: independent, hostile. Date 16 Nov 2026. Sources checked: arXiv:2201.06602 **v1** (only arXiv version; PDF
and LaTeX source `gravity_4d.tex` from arxiv.org/e-print) and the **published version JHEP05(2023)122** (full text
via INSPIRE). Equation numbers below are arXiv v1 unless marked "JHEP" (JHEP shifts §2 numbers by −1).

## Summary verdicts

| aspect | verdict |
|---|---|
| Mathematical correctness (Theorem 1, regions P/C/M/S) | **Correct as far as checked.** No gap in coverage, no invalid inequality found. Minor write-up gaps only. |
| Computational soundness | **Sound.** `run_all.sh` re-run: all 7 PASS, 2 min 04 s. Certified paths use exact integers/rationals (FLINT, sympy) or Arb balls; floats only in print statements. Mutations fail as they should. Two numeric constants are assumed, not asserted (true; minor). |
| Fidelity to the paper | **v1 transcription correct** (sum rules, low sides, kernels, spins, Φ_h, measure, bracketing all verified against the LaTeX source). Finding (b) is **arithmetically correct for v1 (3.3a)**, but the note **ignores the version of record (JHEP)**, which prints a *different* (3.3a)/(3.3b), and it misses evidence that the authors' own numerics used the g5-free functional. Framing as "error in the paper" must be revised (Issues 2–4). |
| Physical applicability | **Theorem 3 is vacuous as stated**: H-IR (`⟨[F^ε]₋⟩ = 0`) contradicts H1–H4 with G > 0 (Issue 1). Theorems 1, 2 unaffected. |
| Novelty | A certified (interval/exact) proof of heavy-state positivity of F_g3 appears **new**: no rigorous follow-up found among the 207 INSPIRE citations of 2201.06602 or by arXiv/web search. Incremental relative to the paper, which already claims positivity from dense plots. No erratum found. |

## Issues

**1. [MAJOR, physics — Theorem 3 vacuous] H-IR is inconsistent with H1–H4.**
H-IR requires that no heavy spectral weight sit where `F^ε < 0`. But Theorem 2 holds for every ε, so the heavy
measure must reproduce the `16πG/p²` pole of `−B₂|low`. Differentiating the B₂ part in ε, in the scaling limit, gives
`∫ρ(b)J₀(εb)db = 16πG/ε²`, i.e. heavy density `ρ(b) ≃ 16πG·b db` at large b: the eikonal states, cf. arXiv:2410.21499.
That weight necessarily covers `ε^{−2/3} ≲ b ≲ ε^{−1}`, which is exactly where `F^ε < 0` (paper eq. (3.6)). So
`⟨[F^ε]₋⟩ = O(G) > 0` whenever G > 0. Theorem 3 is therefore true only vacuously. Its wording is also inconsistent:
"contribution non-negative" is not the same as "`⟨[F^ε]₋⟩ = 0`".
*Fix:* state the hypothesis that is actually needed, `⟨F^ε⟩ ≥ 0`. Alternatively assume that the heavy contribution from
`b > b*` equals its eikonal (EFT-computable) value, and bound or compute the resulting `Δ_IR`, which is O(1) in the constant
B under an eikonal density. Relabel §0, §4.3 and §5 accordingly. Do not call this "the paper's prescription": the paper only
says positivity "should not be crucial there".

**2. [MAJOR, fidelity] The published version is not addressed.** JHEP05(2023)122 prints
`F_g3 = ∫… − 5∂_{p²}(B₄^{(1)} + B₄^{(2)})|_{p=0}` and `F_g4 = ∫… − (7B₄^{(1)} + ∂_{p²}B₄^{(1)} + ∂_{p²}B₄^{(2)})|_{p=0}`, so the
authors changed the forward terms in publication. The surrounding text still says "These use only B^{(1)}". I computed
(sympy, own script) `−B₄^{(2)}|low ≡ 0` through O(g₆), including the graviton, ĝ₃ and matter terms, with B^{(2)} defined as in
(2.26)/JHEP (2.26). So the JHEP term does **not** remove the −5g₅ either. On my reading of `d̃_{4,−4}` (normalisation
binom(J+4,8) at x=1), its heavy action on `+−` states carries `(−1)^J binom(J+4,8)/m^{18}` terms, which look incompatible
with positivity. I have not verified that reading, so treat it as a flag, not a finding.
*Fix:* say explicitly which version is analysed, quote JHEP (3.3a)/(3.3b), and either analyse the JHEP functional or explain
why it is not the intended one. Theorem 1 currently covers v1-literal and the note's own `+5B₅^{(1)}|₀` completion, not the
JHEP-printed functional.

**3. [MINOR→MAJOR in presentation, fidelity] Finding (b) should be presented as a misprint, with evidence.**
Verified independently: the v1 (3.3a) low side is `16πG L − πG k|ĝ₃|² − 5g₅`. However, the paper's own Fig. 4
(`figures/func_g3_pm.pdf`, `func_g3_pp.pdf`, rendered at 600 dpi) agrees with the **completed** functional and not the literal one:
- in the `+−` panel, the point (m=1, b≈8), i.e. J=4, falls in the lowest (<0.1) band, matching completed 19/315 ≈ 0.060 rather than literal 1594/315 ≈ 5.06;
- in the `++` panel, (m=1, b=0), i.e. J=0, lies in the 10–100 band, matching completed 3904/315 ≈ 12.39 rather than literal 2329/315 ≈ 7.39, which would fall in the <10 band.

So (3.4a) and the plots were almost certainly computed with a g5-free combination; the defect is in the printed
formula. The note's sentence "(3.4a) does not follow from (3.3a) as printed" is literally true, but "error in the paper" overstates it.
*Fix:* add this evidence, and call it a typographical omission in v1 that is not repaired in JHEP (see Issue 2).

**4. [MINOR, fidelity] The same omission affects (3.3b)/(3.4b).** The low side of `F_g4` contains `−14g₄ − g₅ + …`
(coefficient of g₅ is −1; own computation). (3.4b) therefore also needs `+B₅^{(1)}|₀` or a caveat. The note says F_g4 is "not treated"
but should flag this. As a check, my computation reproduces the paper's 14.1 = 99/7 and 0.04 = 35249/840840.

**5. [MINOR] Bracketing.** The LaTeX source shows the bracket explicitly:
`(1-p)^3 \big(p(65+155p+47p^2)B_2^{(1)}+2(81-1927p+4601p^2)B_3^{(1)}\big)`. The question is settled directly, so §1.4
"PDF does not show the bracket glyph" can be replaced. An older commented-out version in the source used a different
normalisation (`−(5/2)∂B₄`, extra overall p). It is irrelevant but worth noting.

**6. [MINOR] Matter formula (3.8).** The paper prints `(m_ℓ² − 6m_ℓ p² + 6p⁴)` in (3.8a,b) (macro `\mlight = m_\ell`). It is dimensionally
inconsistent; the correct factor is `m_ℓ⁴ − 6m_ℓ²p² + 6p⁴`. The note silently uses the correct form, and its result reproduces (3.9) exactly.
The note should flag the misprint. "Up to the overall sign convention" is wrong: the signs agree exactly once `B|matter` is
read as `B|low`.

**7. [MINOR] Normalisation of the paper's plots.** The colour-bar label in Fig. 4 is `m⁸F`, while the caption says `m¹⁰`.
The near-zero at m ≈ 21.8 lies outside the plotted range m ≤ 16 anyway. In `m⁸F` it would be ≈ 6, not ≈ 2.9·10³. Rephrase "invisible
in the paper's normalisation".

**8. [MINOR, code] Unasserted constants.** `cert_S.tail_bound` assumes `|T0_n|, |T1_n| ≤ Tmax = 100`, and `small_b` assumes the
n ≥ 400 tail ratio is below 1/2 and |T| ≤ 100. Both are true: max over n of |T0|, |T1| is 12.34, and the ratio is about 10⁻⁷. But neither is
checked in code. *Fix:* add exact asserts. For n > NMAX, `|μ_j| ≤ ∫|ψ| ≤ 65` suffices.

**9. [MINOR, write-up] Lemma M.**
- The step `τ²|A″| ≤ (π/4)sup|w″|` is asserted without proof. It is true: for τ ≥ 1, `τ²∫₀^{arcsin(1/τ)}sin²t dt` is maximal at τ=1 with value π/4, which I checked numerically. Add one line.
- In the `+−` application the text says "value −2624.4/b̂₋⁴". It should read `−2624.4/(b̂² b̂₋²)` (`= −2624.4/(b̂⁴(1−δ)²)`). The code's final inequality is correct.
- Rest of Lemma M re-derived and correct: the Fubini step, `G(0) = πw(0)/2`, `G″(0) = πm²w″(0)/16`, the four integrations by parts using `cos Nπ = 0` and evenness about π, the Faà di Bruno coefficients, and the swap bounds (2/3), 1, 1, `π/(2√u)`.

**10. [MINOR] "This is exactly the paper's (3.4a)".** Only the ε→0 constants (37.8387, −45.4452) match. The uniform B = −45.2059 for
ε ≤ 1/100 is weaker. The note says so in §0, but "exactly" should be softened.

**11. [MINOR] Quoted classics "equation numbers from memory".** Verify the DLMF/Whittaker–Watson/Szegő numbers. I checked (Q2)
numerically, with agreement to 10⁻²², and the Szegő maximum `binom(n+8,n)` at x = −1.

## Point-by-point checks of the requested items

**Fidelity (v1), all confirmed against the LaTeX source:**
- (2.26a,b): kernels, `ds/(4πi)`, and the sign `[−f(s,−p²)+f(p²−s,−p²)]`.
- (2.33), (2.35a): powers `m^{k+4}(m²−p²)^{(k−2)/2}` and `m^{k+3}(m²−p²)^{(k−3)/2}`; σ₊₊ = +1, σ₊₋ = −1.
- (2.34): measure `16Σ(2J+1)∫dm²/m²`.
- (2.16): `d̃₀₀ = P_J`; `d̃₄₄ = ₂F₁(4−J, J+5; 1; y) = P^{(0,8)}_{J−4}` with N₄₄ = 1.
- Spins: `++` even J ≥ 0, `+−` J ≥ 4 (footnote 8).
- Forward derivatives `D′₊₊(1) = 𝒥/2`, `D′₊₋(1) = (J−4)(J+5)/2`, hence `Φ₊₊ = 10𝒥−5`, `Φ₊₋ = 10𝒥−205`. Completion with `+5σ_h`: `10𝒥` and `10𝒥−210`.
- (2.31), (2.32): re-derived independently. (3.9) polynomials reproduced exactly.

**Mathematics:**
- **(P)** Polynomial = `m⁴F` (verified against direct mpmath quadrature at (J,m) = (4,1), (37,3.3), (50,1.7), (245,21.79): agreement to 12 digits). Degree J+1 (`++`), J−3 (`+−`). Covers the closed interval [0,1]: `P(0) > 0`, leaves cover [0,1], endpoints exact-nonzero. Descartes after the Möbius map is correct. Coverage: 200 + 396 = 596 polynomials. The pessimistic Φ covers both variants because Φ enters with coefficient u³ > 0.
- **(C)** Valid. Uses `|d₄₄| ≤ 1` with `(1−up²)^{−4} ≤ (1−p²)^{−4}`, the Szegő bound, the split at p*, the `log B/4` estimate, and the monotonicity bound `16γΣ1/(J+i) ≤ 128γ/(J−3)`.
- **(M)** Valid (see Issue 9). Additional checks:
  - `w(0)` and `w″(0)` values for both channels;
  - `c_j ≥ 0`, `Σc_j = 1`, J+j ≥ 0;
  - `b_j ∈ [b̂(1−δ), b̂(1+δ)]` with δ taken at u₁;
  - minimisation over u;
  - monotonicity in b̂;
  - for u ∈ [u₁, 1/4], increase in u follows from `5b²u₁(1−η) > 324/b₋`.
- **(S)** Valid:
  - the Mehler substitution `φ = θ sin σ`, giving the `(a/sin a)^{1/2}(c/sin c)^{1/2}` factors, is correct;
  - `(ν+½)θ = b̂p s(up²)` makes the bound branch-free;
  - the nonnegative-coefficient bounds for s, arcsin and `(z/sin z)^{1/2}` hold (log z/sin z has positive coefficients);
  - the CG bound is monotone for |Z| ≥ Z_min ≥ 110;
  - the polydisc Cauchy tail formula is correct;
  - the a_k series is truncated with a ratio-checked geometric tail, carried to b̂-derivatives by Cauchy on radius r;
  - Φ is added exactly;
  - cells tile [1,120] contiguously, and `u₀ = ((b_c+h)/801)²` dominates every physical u with J ≥ 400.

  The continuous extension is legitimate: G_h is analytic and equals 𝔉 at the physical points (integer J, so the series
  terminates), and positivity on a superset of physical points is sufficient. Integrality is needed only for that coincidence and
  for Lemma C at integer J, which then extends identically.
- **Coverage:** J < 400 is (P). For J ≥ 400: u ≥ 1/4 is (C); u < 1/4 with b̂ ≥ 120 is (M); u < 1/4 with b̂ < 120 is (S), with b̂ ≤ 1 handled by `small_b`. There are no gaps.

**Computation, run by me:** `bash run_all.sh` returned FINAL VERDICT PASS (cert_P 60 s, 11 582 intervals; cert_M:
`80b−C_R = 2748`, `+−` margin 721; cert_S 604 cells per channel). Mutation tests on a scratch copy:
- `Φ₊₋ −1` makes cert_P FAIL at J=4;
- flipping σ in `+−` makes cert_P FAIL for 33 spins (J=5…);
- `cert_M` with B_A = 60 FAILs;
- cert_S on b̂ ∈ [20,25] with δ subtracted from 𝔉: δ = 10⁻⁴ passes, while δ = 2·10⁻³ and above FAIL. The true minimum there is 2.09·10⁻³ (J=400, b̂≈22.95), so the certificate is sharp to about 5%.

The σ-flip passes cert_S, which is legitimate because σ enters at O(u). No mpmath/numpy is imported in certified scripts.

**Independent spot checks (own sympy/mpmath scripts, no repo code):**
- `m¹⁰F(+−, J=4, m=1) = 1594/315 ≈ 5.0603` (literal) and `19/315 ≈ 0.0603` (completed). Confirmed.
- Further values at m=1, literal / completed:

  | channel, J | literal | completed |
  |---|---|---|
  | `+−`, J=5 | 166.69 | 161.69 |
  | `++`, J=0 | 7.394 | 12.394 |
  | `++`, J=2 | 34.59 | 39.59 |

- `k(0) = 330151/12012 = 27.48510`, `L(ε) = 65 log(1/ε) − 1171/15 + 40ε + 223ε²/2 − 259ε³/3 + 7ε⁴/2 + 47ε⁵/5`, `A₀ = 37.838686`, `B₀ = −45.445232`. All confirmed.
- The g₅ coefficient of the literal action is −5. Confirmed.
- Near-zero: `m⁴F(+−, 245, 21.79) = 2.6665·10⁻⁵`, `b̂³m⁴F = 0.305`. Confirmed.

**Hypotheses:** H1–H5 are honestly stated and Theorem 3 is labelled conditional, but H-IR is inconsistent (Issue 1). H3's use
of unsmeared forward limits `∂B₄|₀` and `B₅|₀` matches the paper's own stance (§2.5, §3.2).

**Novelty and literature search:**
- arXiv lists only v1 of 2201.06602. The published version is JHEP05(2023)122, with a modified (3.3) (Issue 2).
- No erratum was found (web search, INSPIRE record).
- I scanned all 207 INSPIRE citing records. None gives a rigorous or computer-assisted positivity proof of these functionals, and none discusses the g₅ term. Closest related work:
  - 2201.07177 (non-projective bounds, SDP/analytic);
  - 2205.01495 (D ≥ 5 graviton partial waves, same authors);
  - 2410.21499 (graviton pole made of eikonal states; relevant to Issue 1);
  - 2604.15235 (sampling the graviton pole);
  - 2610.19149 (IR-finite 4D positivity);
  - 2304.01259 and 2410.09794 (Hong–Wang–Zhou-type scalar-tensor causality bounds).
- arXiv API searches combining interval arithmetic or computer-assisted methods with positivity, bootstrap or graviton found nothing relevant.

---

## Round 2 (revised note, commit 4f92504)

### Verdicts

| aspect | verdict |
|---|---|
| Mathematical correctness | Theorem 1: unchanged and correct. **Proposition 1′: correct**, re-derived independently. **Lemma IR: algebra correct given E**; one rounding slip and one wrong order-of-magnitude statement (Issues R2.3, R2.4). |
| Computational soundness | **Sound.** `run_all.sh` re-run: 9/9 PASS, 2 min 04 s. The `Tmax` assertion and the precision restore in `cert_S.py` are correct. `cert_ir.py` and `cert_jhep.py` are rigorous (exact or Arb). |
| Physical applicability | **Theorem 3(i): honest**, essentially a rewriting. **Theorem 3(ii): MAJOR problem.** Hypothesis E1 is either vacuous (read literally) or unproved (read approximately), and at ε = 1/100 the negative set N_ε contains ordinary, non-eikonal heavy states (Issue R2.1). |
| Novelty | Unchanged. Proposition 1′ is a genuine finding about the version of record and should be reported to the authors as a misprint. |

### R2 issues

**R2.1 [MAJOR, applicability — Theorem 3(ii) and Lemma IR]: E1 is not an eikonal hypothesis at ε = 1/100.**
(E1) says that the states of N_ε (those with F^ε < 0) are in the scaling regime, with **exact** equality `F^ε = (2/m⁴)f̂_ε(2J/m)`.
- **Read literally**, no finite-mass state satisfies the equality. N_ε must then carry zero weight, which is the withdrawn H-IR again, and that is incompatible with G > 0.
- **Read as an approximation**, Lemma IR is no longer a proof: the finite-m corrections are not bounded.
- **In content**, E1 excludes states that are certainly not eikonal. My own float scan (Gauss–Legendre, scipy `eval_legendre`/`eval_jacobi`, completed functional) finds:
  - at ε = 1/100, `++` states with **m ≈ 3.8–6.5 M, J = 14–28, b̂ ≈ 7.4–9.6** have F^ε < 0. Examples: (J=26, m=6.15): m⁴F = +0.00195 but m⁴F^ε = −0.0418; (J=24, m=5.8): −0.0445; (J=14, m=3.8) is negative too.
  - `+−` states at b̂ ≈ 22.3–22.6, m ≈ 18–27 are negative too, e.g. (J=245, m=21.79): −0.0034.
  - the `++` states remain negative at ε = 10⁻³ and become positive at ε = 10⁻⁴.

  The cause is the term `162uε` in R_ε, of order uε rather than ε², which dominates the nearly marginal `++` island of the paper's Fig. 3a. So for m_IR near M/100, E(λ,η) forbids spectral weight on ordinary M-scale states at b ~ 8/M: exactly the states a cubic-coupling bound is meant to constrain. Calling this an "eikonal model" misdescribes it.

*Fix.* Use the exact inequality `[F^ε]₋ ≤ [R_ε]₊` (Theorem 1), with `m⁴R_ε = ∫₀^ε Ψ_h D_h dp`. Then either:
- (a) certify an ε* such that for ε ≤ ε* every state with `m⁴F < m⁴|R_ε|` has b̂ ≥ b₀ and u ≤ u₀ (large b₀, small u₀); this needs a quantitative lower bound on 𝔉, which cert_S and cert_M can supply. Then state E2 only there, with an explicit `O(u/ε)` correction to `[f̂_ε]₋`. Or
- (b) keep finite-m states in Lemma IR and bound their contribution under a stated measure hypothesis in (m, J), not in the b-marginal ν.

Restrict Theorem 3(ii) and the headline to ε ≤ ε*. For the `++` island this means ε ≲ 4·10⁻⁴ (estimate). Until then, label 3(ii) an estimate, not "[proved]".

**R2.2 [MINOR, applicability] E2 and the IR scale.** E2's density `16πG b log⁺(1/(μb))` with a hard support cutoff at b = 1/μ is a modelling choice: the IR regulator is unavoidable in D = 4 for a positive measure. It is not excluded by H1–H5, and it is consistent with Theorem 2, which constrains only ⟨F^ε⟩, not the N_ε part. So E is non-vacuous once E1 is repaired. The result then depends on λ = μ/m_IR, and the "best" constant −54.6 at λ = 1/2 is obtained by optimising over this model parameter. Say explicitly that μ is a model IR scale and not an observable, and that the comparison with the paper's −45.4 is model-dependent.

**R2.3 [MINOR] Rounding.** `B + Δ_IR = [−40.6008 ± 3·10⁻⁵]` (cert_ir log). The headline "−40.601" is *below* the certified value, so it is slightly stronger than what is proved; write −40.600. The other entry, −54.6 (from −54.63), is rounded in the safe direction.

**R2.4 [MINOR] Theorem 3(i).** It says Δ_IR "only receives contributions from states whose action without cutoff is O(ε²)". From the note's own bound `|R_ε| ≤ (1−ε²)^{−4}(65ε² + 162uε + 17uε⁴)`, the threshold is **O(ε² + uε)**, i.e. O(ε) at u ~ 1. This is exactly what produces the counterexamples in R2.1. Correct the statement.

**R2.5 [MINOR, code] cert_S.** `small_b` asserts the term ratio only for 400 ≤ n < 2000. Beyond that it is obvious, since `λ_{n+1}u/(n+1)² → u ≈ 1.6·10⁻⁶`, but say so. The `_tmax_ok` bound `|T1_n| ≤ ∫ψ₂ + ∫|ψ₃|` is correct, and restoring `set_prec(PREC)` after importing `cert_C` (200 bits) is correct; precision affects only tightness, not validity. In `cert_ir`, `count_roots` on (10⁻⁹, 1/100] leaves (0, 10⁻⁹) unchecked. That gap is trivially fine (`65 − φ = 40p + 223p² − 259p³ + …`, `d₂ ~ 259p³`) but should be stated.

### R2 checks done myself

- **(1) Proposition 1′: independently confirmed.**
  - **Heavy density.** From v1 (2.16a), Im f(m², u) is proportional to `|c⁺⁻|² d̃₄₄(1+2u/m²)`. At fixed t = −p² we have u = p² − m², so the argument is −x. The B₄^{(2)} kernel `(2s−p²)/[s(s−p²)]³` at s = m², with the same dm²-to-⟨⟩ normalisation as the B^{(1)} rules (I re-checked this against (2.33)), gives per state `(2m²−p²) P^{(0,8)}_{J−4}(−1+2p²/m²)/(m¹²(m²−p²)³)`. This agrees with JHEP (2.34b), since `(−1)^J d̃_{4,−4}(x) = d̃₄₄(−x)`.
  - **Forward derivative, by hand.** `∂_{p²}[(2m²−p²)/(m²−p²)³]₀ = 5/m⁶`, and `P^{(0,8)}_n(−1+2q) = (−1)^n binom(n+8,8) ₂F₁(−n, n+9; 9; q)`. Together these give `−5∂B₄^{(2)}|₀ = 5(−1)^J binom(J+4,8)(2𝒥−85)/(9m¹⁸)`.
  - **Exact values, sympy from scratch.**

    | J | `−5∂B₄^{(2)}` term at m = M | `m¹⁰F^JHEP(+−, m = M)` |
    |---|---|---|
    | 4 | −25 | −6281/315 |
    | 5 | +125 | |
    | 6 | −25 | |
    | 7 | −2475 | −28385719/13860 |
    | 9 | −67925 | −11495259817/171171 |

  - **Low side.** `−B₄^{(2)}|low ≡ 0` (EFT through g₆ plus matter). So −5g₅ survives.
  - **No harmless reading.** Any nonzero rescaling c of B₄^{(2)} adds `c·(−1)^J·O(J¹⁰)` at m = M against an action growing only like J². Both signs of c therefore fail for one parity of J. A Regge-convergent fixed-t rule needs the power (k+2)/2; with (k−2)/2 it diverges, since f(s, p²−s) ~ s. The t-channel cut does not enter at fixed t. So no consistent normalisation of B₄^{(2)} makes the printed term harmless, and the JHEP (3.3a) is a misprint. The note's wording ("indicates that the printed forward terms are not the ones used") is appropriate.
- **(2) Lemma IR.**
  - **Proof checked line by line, given E1–E2:**
    - `[f̂_ε]₋ ≤ [∫₀^εψ₂J₀]₊ ≤ 65ε[J₁(εb)]₊/b + 14ε³`, using `65 − φ = 40p + 223p² − 259p³ + …` and `40/3 + 223ε/4 ≤ 14`;
    - the substitution x = εb gives 65 I(λ), and `∫₀^{1/μ} b log(1/(μb)) db = 1/(4μ²)` gives 3.5ε/λ²;
    - the series for I(λ) and its alternating tail bound are right;
    - `∫₀^{R} b log(R/b) J₀(bp) db = (1−J₀(pR))/p²` is right.
  - **Numbers reproduced from the log:** Δ_IR(λ=1) = 4.6051, B + Δ = −40.6008, and `A log(1/2) − 28.4028 = −54.63` (needs μ ≤ M/200).
  - **Non-vacuity and Theorem 2:** see R2.1–R2.2.
- **(3) Overclaims.** Theorem 3(i) is honest, apart from R2.4. Finding 4 and §5 call E(λ,η) "non-vacuous" and Theorem 3(ii) "[proved]"; both should be qualified per R2.1. Findings 2–3 on the misprint are fine.
- **(4) Code.** Diff reviewed; run_all PASS (9/9). Nothing commits or alters data outside `logs/`.

### R2 summary

No fatal issue. One major issue: Theorem 3(ii)/E1 is vacuous if read literally, unproved if read approximately, and excludes M-scale `++` states for m_IR ≈ M/100. Everything else is correct, with minor fixes. Proposition 1′ (the JHEP-printed functional is not positive) is confirmed independently.
