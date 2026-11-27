# NF02 / Omega^EE positivity at every gap >= 1000 — working directory (status file, kept current)

Scope: the sharp NF02 cubic functional Omega^EE, the exact all-gap realisation of the unchanged NF02 target
(sech-convolution inverse, erf cutoff at Lambda+2000; `../uniform_gap/RESEARCH_NOTE.md` Sec. 1), and the target
inequality H_ang,b > 0 on every heavy state (a = Delta-J-2 >= Delta_g-2), b in {0,4}, every real Delta_g >= 1000.
Retained input throughout: Eq. (fp-target) of `prl/sec_fullpos.tex` [D] (only the sign of the prefactor is used).
Nothing outside this directory was modified; no commits.

**Status: final for this session (see Sec. 0).**

## 0. Bottom line

* **NOT closed. No counterexample found.** The open region is J > 200 with bhat > 40 (all gaps > 1000).
* New rigorous results:
  - **Theorem NG / NG'** (proved, Arb-verified for all Lambda >= 998.5): the NF02 generating functions F (k=2) and
    G (k=3) are not non-negative combinations of z^k(1-z)^d for any exponent bound D <= Lambda^2/2000 — in
    particular not for D = M+1 (Theorems RP/MP) nor for the twist-edge range a+1. Direction (1) is impossible, also
    with "controlled remainder" (the remainder must carry all of F's e^{0.88 Lambda} growth at z = -1/4).
  - **Lemma PB+**: Beta–Cauchy representation L(f) = E_B[[z^n](1-z)^{-a-1} f(Bz)]/binom(N,n), giving a two-line
    proof of the positive-basis lemma with the state's own range d <= a+1 (real twist).
  - **Lemma DIL**: H(state; Lambda) = 2 int psi2(p) K_st(Lambda p) dp + sigma Lambda^2 int psi3(p) K3_st(Lambda p) dp with
    gap-independent kernels; all-gap positivity = positivity of all dilations Lambda in [998.5, sqrt(a(a+1))].
  - **Certified fixed-spin positivity for 65 <= J <= 200** at every real gap >= 1000 and every mass, both channels
    (204 exact-integer Bernstein certificates with the package's method; ideal family. Transfer to the erf-cut
    realised family: J <= 140 with the package's cutoff constants, J <= 200 with a routine sharpening of the
    chi-tail rate — caveat C1, grade D).
* New numerical evidence (labelled): ~2150 valid points at gaps 1000, 1001.5, 1500, 3000 in the open region (J up to
  1.4e5, bhat up to 5.8e3, ell <= 1.24): all positive; per-frequency kernel study; scaling test showing the finite-mass
  deviation H/H_LR - 1 falls ~Lambda^-2 under joint scaling. Large ell (>~1.3) at large N was NOT reachable by my
  evaluator at any gap != 1000.
* The sharp all-gap rigidity bound tau <= (240/Lambda^4) sqrt((3023 log(2Lambda/3)-5307)/0.99995) therefore remains
  conditional on positivity in the open region.

## 1. Ledger

| region (all real Delta_g >= 1000 unless stated) | method | status | retained inputs | where |
|---|---|---|---|---|
| J <= 64, all masses | exact Bernstein in q (package) | certified (previous work) | fp-target | `../uniform_gap/` |
| J >= 65, bhat <= 40 | tensor Bernstein (package) | certified (previous work) | fp-target | `../uniform_gap/` |
| 65 <= J <= 200 (Sec. 4), all masses | package's `check_physical.py`, copied unchanged (sha256 in `cert_copy/`), run on new spins | ideal family: certified (exact integer Bernstein). Realised family: certified given the routine transfer extension of caveat C1 (J <= 140 with the note's constants; J <= 200 with the rate-3/2 tail) | fp-target | `cert_copy/phys_*.log, physical_*.json` |
| positive-cone representation of F, G (route 1) | Theorem NG | **proved impossible** | none | `nogo_cone.py`, Sec. 2 |
| J >= 65, bhat > 40, N=n+a <~ 2e4 (ell <~ 1) | direct Arb evaluation, gaps 1000/1001.5/1500/2000 | numerical only (truncation not certified); all positive | — | `hev.py`, `sweep1.py`, `out/` |
| J > 200, bhat > 40 (all gaps > 1000; at 1000 itself covered only for the published weight by Thm FP) | — | **open** | | |

Caveat C1 (ideal vs realised family). `check_physical.py` certifies the *uncut ideal* closed form. RESEARCH_NOTE
Sec. 4 transfers ideal -> realised (erf cutoff at Lambda+2000) with an envelope proved for J <= 64. Its two
ingredients extend as follows (my estimate, grade D = routine, not re-derived line by line):
 * coefficient envelope (k!)^2 t^k |g[k,3]((nu+-i)^2)| <= (1+y^2)^k: numerically true for all k <= 201 at 18 sample
   frequencies nu in [0,1e6] (max log-ratio -5e-6, tight at k=1 as the review noted) — `envelope_check.py`;
   the ladder sum then has envelope ~ exp(E_J(1+y)), E_J = 2J+12;
 * exterior cutoff piece x > Lambda+1000: with the note's tail |chi R| <= 2^90 e^{-(x-L)} the exponent is
   ~ -1000 + 3E_J, which stays below -100 only for J <= ~140. The true tail rate is pi/2 (it comes from the sech
   kernel of the inverse; the note weakened it to 1 to absorb (1+s)^12); keeping rate 3/2 costs a factor <= 2^72 and
   gives -1500 + 3E_J + 50 < -100 for J <= ~200.
So: J <= 140 transfers with the note's constants as written; 141 <= J <= 200 needs the (routine) rate-3/2 version of
the chi-tail bound. Until that is written out, the 141..200 certificates are for the ideal family only.

## 2. Theorem NG (no positive-cone representation) — proved

Let C_D be the convex cone of power series sum c_{k,d} z^k (1-z)^d with c_{k,d} >= 0 and 0 <= d <= D, convergent on
|z| < 1.

Lemma NG.1. For G in C_D and 0 < rho < 1: |G(-rho)| <= ((1+rho)/(1-rho))^D G(rho).
Proof: |(-rho)^k (1+rho)^d| = rho^k (1-rho)^d ((1+rho)/(1-rho))^d and d <= D. []

Theorem NG. For every real Lambda >= 998.5 the ideal NF02 generating functions
F(z) = 2(1-z)^-3 int_0^1 psi2(p) P_{-1/2+i Lambda p}((1+z)/(1-z)) dp and
G(z) = Lambda^2 (1-z)^-4 int_0^1 psi3(p) P_{-1/2+i Lambda p}((1+z)/(1-z)) dp satisfy
|F(-1/4)| > (5/3)^(Lambda+1) |F(1/4)|, |G(-1/4)| > (5/3)^(Lambda+1) |G(1/4)|. Hence F, G are not in C_D for any
D <= Lambda+1, at any gap. At Lambda = 998.5 the margin is e^345 for both.

Proof. At z = -1/4, x = 3/5 = cos(theta). Mehler–Dirichlet gives Phi(nu) := P_{-1/2+i nu}(3/5) =
(2/pi) int_0^theta cosh(nu t)(2(cos t - cos theta))^{-1/2} dt, positive and increasing in nu, with
sinh(nu theta)/(pi nu sin(theta/2)) <= Phi(nu) <= cosh(nu theta) (2/pi) K(sin(theta/2)).
psi2 (psi3) has exactly one sign change in (0,1), at p* = 0.91682 (0.92328) (Sturm count, exact), positive before,
negative after. Monotonicity of Phi gives
int psi Phi <= Phi(p_hi Lambda) int_0^{p_hi} psi^+ - Phi(0.95 Lambda) int_{0.95}^{0.99} |psi| < 0,
so |F(-1/4)| >= 2(4/5)^3 [Phi(0.95 Lambda) I_- - Phi(p_hi Lambda) I_+]. For x >= 1, |P_{-1/2+i nu}(x)| <= 1, so
|F(1/4)| <= 2(4/3)^3 int|psi2|. After dividing by exp(0.95 Lambda theta)/Lambda the left side increases and the right
side decreases for Lambda >= 998.5 (exponent gaps 0.0307 and 0.370 per unit Lambda), so the Arb check at
Lambda = 998.5 (`nogo_cone.py`, exact rational polynomial integrals, Arb for theta, K, sinh, cosh) proves all Lambda. []

Consequences. (i) Any proof that writes F (or G, or a positive mixture of RP/MP-type profiles) as a non-negative
combination of z^k(1-z)^d, d <= M+1, is impossible; so is the state-dependent version with d <= a+1 at twist-edge
states (a+1 = Delta_g-1 < Lambda+1/2+10^-3). (ii) The dressed functions (1+D_b/T)F + sigma_b G/T converge pointwise
on (-1,1) to F as T -> infinity at fixed b, and the inequality of NG.1 is closed under pointwise limits, so for the
gap-fixed cone C_{M+1} the dressed functions are also outside the cone for all sufficiently heavy states.
(iii) "Cone + controlled remainder" F = Phi_+ + R with |R| <= eps on [0,1/4] forces |R(-1/4)| >= |F(-1/4)| -
(5/3)^(Lambda+1)(|F(1/4)|+eps), i.e. R carries essentially all of the e^{0.88 Lambda} growth of F at negative z; its
Taylor coefficients are then exponentially large and bounding L(R) termwise is the original problem.
Theorem NG' (strong form; proved, `nogo_cone_strong.py`). The same conclusion holds for every D <= Lambda^2/2000, at
every Lambda >= 998.5: take rho = (600/Lambda)^2, so X = Lambda*theta in [1081, 1200] is bounded and the cone factor
((1+rho)/(1-rho))^D <= exp(2 rho D/(1-rho0)) <= exp(360/(1-rho0)) = e^{563.5} (rho0 = (600/998.5)^2); every bound is monotone in X or rho and is
evaluated at its worst end (log LHS >= 1005.6 vs log RHS <= 567.3 for F; 1005.3 vs 568.0 for G).
So the cone method is excluded for exponents up to Lambda^2/2000, while the basis lemma itself is false at
d ~ M^2/4 (mixed-profile counterexample): only the window Lambda^2/2000 < d < M^2/4 is not covered by either
statement, and there the basis lemma is not known to hold.
Mechanism: the cone allows growth e^{2 rho D} at -rho, while F grows like e^{2 Lambda sqrt(rho)} (Bessel type,
coefficients ~ Lambda^{2k}/(k!)^2). Heuristically the obstruction persists for D up to ~ Lambda^2/log Lambda; the
exact counterexample in `../uniform_gap/mixed_profile` shows the basis lemma itself fails at d ~ M^2.

## 3. Lemma PB+ and Lemma DIL — proved

Lemma PB+ (Beta–Cauchy form of the ladder). With n = J-b, real twist a > 0, N = n+a, s = 2b+1 and
B ~ Beta(n+s, a): U_{b,k} = E[B^k] (n)_k↓/(N)_k↓ and, for every power series f,
L(f) = E_B [ [z^n] (1-z)^{-a-1} f(Bz) ] / binom(N, n).
Consequently, for integers k >= 0 and 0 <= d <= a+1: L[z^k (1-z)^d] >= 0 (strict for k <= n), because
(1-Bz)^d (1-z)^{-a-1} = (1 + (1-B) z/(1-z))^d (1-z)^{-(a+1-d)} has non-negative coefficients.
Proof: (n+s)_k/(N+s)_k = E B^k; and [z^n](1-z)^{-a-1} z^k / [z^n](1-z)^{-a-1} = Gamma(N-k+1)Gamma(n+1)/
(Gamma(n-k+1)Gamma(N+1)) = (n)_k↓/(N)_k↓ (valid for real a). []
This is the positive-basis lemma of Sec. S11.K with M+1 replaced by the state's own a+1 (equal at the twist edge,
larger above it), with a two-line proof. Exact-rational sanity checks: `pb_lemma_check.py` (0 violations).

Lemma DIL (dilation form; exact). Let K_st(nu) := (1-r)^3 L_st[(1+D_b/T)(1-z)^-3 P_{-1/2+i nu}((1+z)/(1-z))] and
K3_st(nu) := (1-r)^3 L_st[(1-z)^-4 P_{-1/2+i nu}((1+z)/(1-z))]/T. These depend only on the state (J, Delta, b), not on
the gap. For the ideal family
H_ang,b(state; Lambda) = 2 int_0^1 psi2(p) K_st(Lambda p) dp + sigma_b Lambda^2 int_0^1 psi3(p) K3_st(Lambda p) dp,
and all-gap positivity is: for every state with twist a >= 998 and every Lambda in [998.5, sqrt(a(a+1))], this
dilation of the fixed NF02 profile is positive. (Delta_g enters only through the dilation Lambda and the constraint
Lambda^2 = (Delta_g-2)(Delta_g-1) <= a(a+1).) Proof: fp-H/fp-V with F = 2(1-z)^-3 int E P, E(nu) = psi2(nu/Lambda);
the moments are linear in E. []
Use: at Delta_g = 1000 only Lambda = 998.5 is tested on each state; all-gap positivity tests the whole dilation
interval on each state. Under (J, a, Lambda) -> (lam J, lam a, lam Lambda) the leading-Regge kernel P(eta) is invariant
up to O(1/a) (eta depends on r = C/T only), which suggests the scaling limit of Sec. 5.

## 4. Certified extension of fixed-spin positivity (all gaps)

`cert_copy/` holds byte-identical copies of the package scripts (`SHA256_of_copied_sources.txt`); outputs are
written only there. `check_physical.py JMAX JMIN` builds, for each spin-channel, the exact integer polynomial of the
uncut ideal action in (s, zeta) in [0,1]^2 (s = 998/(Delta_g-2), zeta = m_gap(J)^2/m^2: every real gap >= 1000 and
every mass above the twist edge) and proves positivity by exact Bernstein/de Casteljau. (b=0 is run for even J only,
as in the package: H_ang,0 is needed only for even J.)

| spins | channels | result | log |
|---|---|---|---|
| 65-105 | 61 | all PASS | `cert_copy/phys_65_105.log` |
| 106-130 | 38 | all PASS | `cert_copy/phys_106_130.log` |
| 131-160 | 45 | all PASS | `cert_copy/phys_131_160.log` |
| 161-180 | 30 | all PASS | `cert_copy/phys_161_180.log` |
| 181-192 | 18 | all PASS | `cert_copy/phys_181_192.log` |
| 193-200 | 12 | all PASS | `cert_copy/phys_193_200.log` |

Total: 204 new spin-channels (65 <= J <= 200), all PASS; with the package's 94 channels, every physical channel with
J <= 200 is certified for the ideal family at every real gap >= 1000 and every mass.

Normalised lower bounds decrease with J only through the normalising product (998+2J+1)^M (e.g. 1e-35 at J=193),
not through H; they are exact positive rationals.  The J=65..130 part independently confirms (by re-running the
package's own method, not by independent code) the collaborator's unreviewed claim of fixed-spin certificates for
J <= 130.

## 5. Numerical evidence (not proofs)

### 5.0 Sweeps (`sweep1.py`, `sweep2.py`, `summarize.py` -> `out/sweep_summary.md`)
Grid: J in 65..2e5 (log grid), twist a = a_g(1+x), x in {0, 1e-3, .01, .03, .1, .3, 1, 2, 4, 10, 30} (sweep1) and
large-J near-gap grid (sweep2), both b, only bhat > 40. A point is 'valid' if the ladder terminates before the
truncation K (then the Arb value is a rigorous enclosure of the ideal action at that state — GL quadrature is exact
for the polynomial moments) or if the last retained terms are < 1e-25 |S|. Invalid points (truncation not resolved)
are those with N = n+a large and ell >~ 1: the direct method cannot reach them (tail coefficients ~ e^{pi Lambda/2}).

| file | valid/total | nonpositive | min S=H/w^3 (state) | max ell | max bhat | max J | max abs(H/H_LR-1) (r<0.45) |
|---|---|---|---|---|---|---|---|
| sweep1 gap 1000 | 502 / 798 | 0 | 2.42e-05 (J=420,b=4,a=998) | 1.01 | 2.63e+03 | 87000 | 4.77e-03 |
| sweep1 gap 1001.5 | 502 / 798 | 0 | 2.41e-05 (J=420,b=4,a=999.5) | 1.01 | 2.63e+03 | 87000 | 4.77e-03 |
| sweep1 gap 1500 | 496 / 798 | 0 | 9.93e-06 (J=670,b=4,a=1498) | 0.658 | 2.89e+03 | 87000 | 5.57e-03 |
| sweep1 gap 3000 | 560 / 798 | 0 | 2.23e-06 (J=1070,b=4,a=2998) | 0.657 | 5.78e+03 | 140000 | 6.04e-03 |
| sweep2 gap 1000 (K=6000) | 90 / 350 | 0 | 2.52e-05 (J=1000,b=4,a=1996) | 1.24 | 3.14e+03 | 10541 | - |

The largest |H/H_LR - 1| always sits at the smallest spins (J = 65-70, bhat ~ 45-65), where H_LR (Legendre proxy) is
not the right comparison (finite-J Wigner effect), not at the finite-mass-critical states. The minima of S move with
the gap as the scaling (J,a) ~ Lambda predicts (S ~ bhat^-3 in the flat regime).

### 5.1 Per-frequency kernel (Lemma DIL), `pernu.py`, `out/pernu.log`
K_st(nu) computed exactly-in-Arb on nu = 0,5,...,1000 for (J,b,a) = (420,4,998), (420,0,998), (2000,0,998), (1000,4,1500):
* no exponential hump at the per-frequency level: |K_st(nu)| <= 0.98 for nu <= 900 in all four states; it grows only
  near the turning point nu ~ a+1 (|K| = 9.4 at nu = 1000, b = 4, a = 998), consistent with a
  [1-nu^2/(a+1)^2]^(-3/4)-type multiplier;
* the finite-mass deviation is tiny at small frequency (|K-P|/|P| = 3e-8 ... 3e-5 at nu = 0) but O(|P|) for
  nu >~ 200 (the finite-mass kernel is a dephased copy of P_{-1/2+i nu}(eta) at large nu): max |K-P| ~ 0.1 (b=0) to
  9.4 (b=4, at the turning point).
* Consequence for proof strategy: |H - H_LR| <= 2 int |psi2||K-P| is ~1, while H ~ 2e-5 for (420,4,998): the
  p-integral cancels by a factor ~3e4. Any uniform proof must keep the oscillation in nu (contour deformation /
  endpoint expansion at p = 0), exactly as the Delta_g = 1000 routes T', T'', signed transfer did. Termwise-in-nu
  bounds cannot work, just as termwise-in-k bounds cannot (fact (a) of S11) and positive-cone bounds cannot (Thm NG).
  The p = 1 endpoint is harmless away from extreme impact: psi2 has a fifth-order zero there and the turning point is
  at distance >= 1/2 beyond nu = Lambda ((a+1)^2 - Lambda^2 >= Delta_g - 1).

### 5.3 Scaling of the finite-mass deviation (`scaling.py`, `out/scaling_table.md`)
dev := H/H_LR - 1 at 16 base states at gap 1000 and at the doubled states (2J, 2a) at gap 2000:
ratios dev(1000)/dev(2000) = 3.3-6.2 for 13 of the 16 states (i.e. ~Lambda^-2), 2.5 and 31 for two more, and one
state ((401,4,998), dev = -1.3e-6 -> -3.5e-6) where the (tiny) deviation grows. So in scaled variables the
infinite-gap limit of the full action is the leading-Regge action (numerical evidence).

## 6. Files
* `hev.py` — Arb evaluator of the exact action for any gap (validated against `../review_uniform_gap/rv_closed.py`
  exact rationals to all printed digits: `validate_hev.py`).
* `sweep1.py`, `analyze1.py` — open-region sweeps.
* `scaling.py` — scaling test.
* `nogo_cone.py`, `nogo_cone_strong.py` — Theorems NG, NG'; `out/nogo_cone*.json`.
* `pernu.py` — per-frequency kernel study (`out/pernu.log`); `envelope_check.py` — caveat C1 envelope check.
* `pb_sharp.py` — attempted scan of the sharp d-range of the basis lemma at large n (inconclusive: Arb precision lost for d >= 5000; not used).
* `pb_lemma_check.py` — Lemma PB+ exact sanity checks.
* `cert_copy/` — certified fixed-spin runs.

## 7. What would be needed to finish (assessment)

Open region after this session: J > 200 (ideal family; J > 140 for the realised family until caveat C1 is written out) and bhat > 40, at
every real gap > 1000 (at Delta_g = 1000 itself Theorem FP covers the *published* weight, which differs from the
all-gap family by 1e-10..1e-5 relative; FP does not transfer automatically). In the ladder variable this is
every L in (~b̂^2/4Λ^2, infinity) for large spin, including the moving near-gap band a in [Delta_g-2, ~4.5 Lambda].

Why the obvious routes fail (now proved or quantified):
1. positive cone (route 1): impossible for every exponent bound up to Lambda^2/2000 (Theorems NG, NG');
2. termwise in the ladder index k: cancellation e^{2 Lambda sqrt(c)} (and e^{pi Lambda/2} in the tail) — fact (a) of S11;
3. termwise in frequency nu (Lemma DIL): per-frequency finite-mass deviation is O(|P|) for nu >~ 200, while H is
   a 1e-4..1e-5 fraction of int|psi2||K| — Sec. 5.1;
4. q-truncated tensor certificates: the q-degree truncation needs X K^2/Lambda^2 << 1 and the v-truncation
   (b+K)^2 v << 1; both fail beyond bhat ~ 60-80 at J ~ 65-200 (my X = 1600 attempt: truncation majorant 3e126).

What the numerics say the true structure is: H/H_LR - 1 is <= 0.6% wherever computed and decreases ~Lambda^-2 under
(J,a,Lambda) -> lam(J,a,Lambda) (Sec. 5.3), so in the scaling variables (J/Lambda, a/Lambda, ell) the infinite-gap
limit is exactly the leading-Regge action, for which H_hyb >= 0.199 N_s is already proved at every gap. The gap is a
proof-technique gap, with a factor >~30 of room in N_s units.

Most promising path (a q-uniform Regge-contour route, i.e. the Delta_g = 1000 routes T'/T''/transfer redone in
scaled variables):
 (R1) uniform per-power multiplier bounds for S(alpha) = L((1-z)^-alpha) = 3F2(alpha,-n,n+s;-N,N+s;1) on
      1/4 <= Re alpha <= 5, |Im alpha| <= Lambda <= a+1/2, uniform in n and in a >= Lambda-1/2, *including the
      deviation* |S(alpha) d^alpha - N_a(alpha)| <= C |alpha|^2 h-type bounds (the collaborator's eigenframe bound
      1024/[1-nu^2/(a+1)^2]^(3/4) is the boundedness half);
 (R2) localisation of the connection series in scaled variables (Lemmas Loc', Loc'' with Lambda kept symbolic);
 (R3) contour deformation in p for the two Hankel halves of P_{-1/2+i Lambda p}: main term = hybrid floor (proved
      for all gaps), remainders = R1 x R2 x e^{-tau beta}; the p = 1 endpoint is suppressed by psi2's fifth-order
      zero and the turning point lies >= 1/2 beyond nu = Lambda;
 (R4) the flat strip 40 < bhat < B0 (B0 where R3's asymptotic regime starts, plausibly a few hundred) by q-uniform
      tensor certificates on v-cells (engineering, not conceptual);
 (R5) extreme impact (floor e^{-ell/2}/(Lambda^2 ell)) — the collaborator's ell >= log Lambda + 28 result, if confirmed.
Each of R1-R3 has a Delta_g = 1000 instance in S11.F-G; the work is to keep Lambda as a variable (q in [0,1]) in every
inequality. I estimate this at the scale of the original §§100-129 effort, not hours.
