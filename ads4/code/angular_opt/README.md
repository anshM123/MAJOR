# angular_opt: optimising the positivity-free angular bound tau_*(Delta_g)

Research ledger (agent run, 28 Nov 2026). No file outside this directory was edited, and nothing was committed.
Grades: [P] proved (exact or hand proof); [P|X] proved given the input X; [C] computer-verified exactly at the stated
finite set; [NE] numerical evidence; [D] derived, not checked by machine; [open].

Notation (supplement Sec. sec:angular): T = T_Delta = (Delta-1)(Delta-2), C = J(J+1), and C_i = i(i+1).
r_k^(b)(Delta,J) = |a_{O,b,k}|^2/|a_{O,b,0}|^2 = D_k^2 prod_{i<k} (C - C_{b+i})/(T - C_{b+i}), with D_k = (k+1)(k+2)(k+3)/6
(the all-J ladder, exact). The endpoint weight is W^(b) = (T - b(b-1))/3 + 3/2, and |c|^2 = (t_4^2 + alpha_T^2)/16 = tau^2/16.
Heavy means twist Delta - J >= g := Delta_g. Ladder support: b=0 only for even J (J = 0, 2, 4, ...); b=4 for every J >= 4; odd
J >= 5 is b=4 only; J = 1, 3 have no structures (LADDER_ALLJ, KERNELS_D3 §5.5).

## Milestone 1 — the 1/Delta_g comes from discarding exact identities (28 Nov)

### 1a. Where 1/g comes from [P]
* The paper uses two rules: the endpoint rule  sum_{O,b} W^(b) |a_{O,b,0}|^2 = 3 + D,  and the k = 2 ANEC rule,
  sum_O |a_{O,4,2}|^2 >= |c|^2.
* Per heavy state, the ratio of (endpoint weight) to (fourth-harmonic weight) is ~ T/(300 x^2), with x = C/T. It is
  minimised at J ~ 2g, Delta ~ 3g, giving K ~ (729/4800) g^2, hence tau ~ 1/g.
* For these two rules alone the paper's K(g) is exactly optimal. A single b=4 state at the argmin is a primal point that
  saturates it. So 1/g is sharp **for that pair of rules**.

### 1b. The unused identities: ANEC harmonic balance R_k, for every k [P, same argument as the paper's (d)]
* The paper's argument for k = 2 works for every m = k+2. [E_m, E_{-m}] = 0 in |T,+2> (commuting POVM) gives
      (R_k)   sum_O |a_{O,4,k}|^2 - sum_O |a_{O,0,k}|^2 = delta_{k,2} |c|^2 ,   k = 0, 1, 2, ...
* The bookkeeping:
  * E_{k+2}|T,+2> has helicity k+4, which is reached only by the b=4 ladder at index k.
  * E_{-k-2}|T,+2> has helicity -k, which is reached by the b=0 ladder at index k, and by T itself only at k = 2 (helicity -2,
    amplitude c).
  * Heavy states have no helicities 1, 2, 3 (F_1 = F_2 = F_3 = 0).
* What the paper uses: only R_2, and it drops the b=0 term there.
* Why the other R_k matter: on a heavy state r_k^(4) = r_k^(0) (1 + O(k^2/T)). So the R_k force the two helicity ladders to
  balance harmonic by harmonic, and |c|^2 can only come from O(1/T) asymmetries or from low-J ladder ends.

### 1c. The LP and its exact value [P for the primal; NE for the dual optimum]
* **Primal.** Maximise |c|^2 over nonnegative heavy measures rho_b(Delta,J), subject to the endpoint rule (= 3, D = 0) and R_k
  for all k.
* **Dual.** Multipliers lambda_k on R_k, normalised to multiplier 1 on the endpoint rule. With L := -lambda_2:
      b=4, J>=4:      W^(4) + sum_k lambda_k r_k^(4) >= 0
      b=0, J even:    W^(0) - sum_k lambda_k r_k^(0) >= 0
  This gives |c|^2 <= (3+D)/L, i.e. tau^2 <= 16(3+D)/L (plus the same light-state term L_g as the paper).
* **Exact primal point** (`primal3.py`, exact rationals). Three states at the gap:
  * A = (b=0, J=2, Delta=g+2), with ladder entries k <= 2;
  * B = (b=4, J=4, Delta=g+4), with k = 0 only;
  * C = (b=4, J=6, Delta=g+6), with k <= 2.
  They satisfy R_0 and R_1 by choice of weights, and R_k for k >= 3 trivially. The endpoint rule is used with total 3. The
  result is |c|^2 = rho_C r_2^(4)(C) - rho_A r_2^(0)(A) > 0.
* **Consequence [P].** No bound derived from {endpoint rule, all R_k} can be stronger than
      tau_LP(g) = 4 sqrt(|c|^2_primal(g)),  tau_LP g^3 -> 4 sqrt(21600) = 587.88   (g -> infinity).
  Values: tau_LP = 0.2841 (g=10), 0.04965 (20), 3.981e-3 (50), 5.398e-4 (100), **5.828e-7 (g=1000)**.
  So within this framework the best possible scaling is **g^-3**, not g^-4.
* **LP check.** The full LP (HiGHS, `lp_primal.py`, `lp_dual.py`, `run2.py`) at g = 10, 20 reaches exactly this value:
  * the dual needs about n ~ 0.6 g multipliers;
  * the optimal dual equioscillates over all J, up to J >> g;
  * truncating J at 6g gives a fake dual that fails at larger J, which is fixed by including J up to 2000g.

### 1d. A closed-form dual: the Racah family [NE]
* The family: lambda_k D_k^2 = -W_B (-n)_k (n)_k/(k!)^2, with W_B = W^(4)(g+4).
  * Lambda^(0) is a balanced terminating 4F3: -W_B 4F3(-n, n, -J, J+1; 1, 2-Delta, Delta-1; 1).
  * Lambda^(4) has the same form with (4-J, J+5; 6-Delta, Delta+3).
  * Both are Racah polynomials, with (alpha, beta) = (0, -1) forced by the ladder.
  * The scaling limit is Legendre/Jacobi with c = n^2/T, i.e. Bessel J_0.
* It gives L = W_B n^2 (n^2-1)/400.
* Feasibility, scanned on J <= 6-8g and a Delta grid (`racah_explore.py`, `racah_nmax.py`):
  * n_max = 13 at g=20, so tau = 0.0638;
  * n_max = 28 at g=50, so tau = 5.84e-3;
  * at g=100, n = 58 already fails, so n_max < 58.
* The binding states are b=4, J = 6, 7 near the gap.
* In the scaling limit, feasibility reduces to P^(0,8)_m(1-2c) <= 1 for all m, which holds iff c <= 0.2550 (the m=3 root). So
  n_max/g -> 0.505, and tau -> 941/g^3.

### 1e. The positive-mixture lemma, partial [P for J <= twist - 2; NE beyond; fails for J >~ twist^2/3]
* The ladder sequence can be written as m_k = prod_{i<k}(C-C_i)/(T-C_i) = [(J)_k/(N)_k falling] x [(J+1)_k/(N+1)_k], with
  N = Delta - 2.
  * For integer N, the first factor is a hypergeometric factorial moment and the second a Beta moment. Binomial thinning then
    gives m_k = E[C(Y,k)]/C(J,k), with Y in {0..J} random. Hence F = E[Q_n(Y)], with Q_n a Hahn polynomial.
  * For real N this holds when J <= N - J, i.e. J <= twist - 2 [P].
  * Numerically it holds up to J ~ twist^2/3 and fails beyond (`mixture_test*.py`, `fk_test.py`).
* This is the remaining obstacle to a proof for all J with n ~ g: the large-spin band J >> twist^2, with exp(2n) cancellations.

## Milestone 2 — Lemma M and a rigorous tau_* ~ 85/g^2 (c_T -> infinity limit) (28 Nov)

### Lemma M (positive mixture representation of the ladder) [P]
* **Setup.** Take a state of twist h > 2 and write d := h - 2. For a ladder b in {0,4} set:
  * b = 0: J_b = J, N_b = Delta - 2, a_b = J + 1;
  * b = 4: J_b = J - 4, N_b = Delta - 6, a_b = J + 5.
* **Factorisation.** For both ladders N_b - J_b = d, and
      m_k^(b) := prod_{i<k} (C - C_{b+i})/(T - C_{b+i}) = [prod_{i<k} (J_b - i)/(N_b - i)] * [ (a_b)_k / (a_b + d)_k ].
* **The random variables.** Take any integer M >= 1 with M - 1 < d.
  * U is generalised hypergeometric: P(U=u) = C(J_b,u) C(d,M-u)/C(N_b,M), for u = 0..M.
    * This is >= 0 because C(d,m) > 0 for m <= M < d+1.
    * It sums to 1 by Chu–Vandermonde.
    * Its factorial moments are E C(U,k) = C(M,k) prod_{i<k}(J_b-i)/(N_b-i), a polynomial identity in d.
  * p ~ Beta(a_b, d) has moments E p^k = (a_b)_k/(a_b+d)_k.
  * Y | (U,p) ~ Binomial(U,p).
* **Conclusion.** E C(Y,k) = C(M,k) m_k^(b) for k <= M. So for every polynomial
      q(z) = sum_{k<=n} c_k C(z,k)/C(M,k), with n <= M,
  one has  sum_k c_k m_k^(b) = E q(Y),  with Y in {0,...,M}.
* **Scope.** Valid for every J and every real Delta, with no cancellation issue.
* **Check.** Exact rational check on 40 random cases (`mixture_lemma_exact.py`).

### Theorem N (new positivity-free bound) [P | paper's angular-route inputs (endpoint rule with D, all-J ladder, ANEC POVM) + Lemma M; finite certificate C]
* **Hypotheses.** Let g > 2 and M = ceil(g) - 2. Let q have degree n <= M, with q(0) = 1 and -theta_g <= q(z) <= 1 for
  z = 0..M, where theta_g = W^(0)(g+2)/W_B and W_B = W^(4)(g+4).
* **Dual feasibility.** The multipliers
      lambda_k = -W_B C(M,k) (Delta^k q)(0) / D_k^2   (k <= n)
  are dual-feasible on all heavy states:
  * b=4: F <= 1 <= W^(4)/W_B;
  * b=0: F >= -theta_g >= -W^(0)/W_B;
  * J = 0 and J = 4 are trivial, since there F = 1.
* **The bound.**
      tau^2 <= 3200 (3 + D + E_Lambda) / ( W_B M (M-1) S ),   S := q(2) - 2 q(1) + 1 <= 2(1+theta_g),
  with the light term E_Lambda = sum_{light O,b} rho (|Lambda| - W)_+, which tends to 0 under A2'.
* **Certificates** (`cert_lemmaM2.py`, exact integer arithmetic on the whole grid; `cache/certI_*.json`). Each covers every
  gap g with Mb >= ceil(g)-2, n <= ceil(g)-2 and theta_g >= theta_a, by restricting q:
  * (Mb, n, g_a) = (98, 25, 30) gives S = 3.7825;
  * (998, 80, 100) gives S = 3.9249;
  * (998, 80, 1000) gives S = 3.99205;
  * (9998, 250, 1000) gives S = 3.99205.
* **Values** (`tau_new.py`; D = 0 and light term 0):

  | g | tau_new | tau_new g^2 | paper tau_* | ratio |
  |---|---|---|---|---|
  | 32 | 8.60e-2 | 88.1 | 0.595 | 6.9 |
  | 100 | 8.575e-3 | 85.7 | 0.1815 | 21 |
  | 300 | 9.519e-4 | 85.7 | 5.97e-2 | 63 |
  | **1000** | **8.494e-5** | 84.9 | 1.781e-2 | **210** |
  | 10^4 | 8.494e-7 | 84.9 | 1.778e-3 | 2094 |

* **Asymptotics.** tau_new -> 40 sqrt(3/(W_B M(M-1)(1+theta))) -> 60 sqrt(2)/g^2 = 84.85/g^2.
* **Independent check.** The certificate's lambda_k were evaluated directly (no Lemma M) on about 4000 heavy states, at
  120 digits (`cert_direct_check.py`). All margins are >= 0, with equality only at the tight state (b=4, J=4, at the gap).
* **Limit of the method [P].** Lemma M cannot beat g^-2: c_2 = C(M,2) S <= C(M,2) 2(1+theta). The framework optimum is
  588/g^3 (Milestone 1).

### Caveat: light states [NE, important]
* **The light term.** Light states (2 < twist < g) enter through the light term. For the power-improved functional,
  |Lambda|/W on light states is astronomically large (`light_check.py`, g=1000):
  * 1e58 near twist 2;
  * 2.5e50 at twist 6;
  * 1e11 at twist 50;
  * 1e3 at twist 80;
  * about 1 near the gap.
* **Cause.** The partial sums of the Newton coefficients c_k = C(M,k) Delta^k q(0) are exponentially large.
* **Consequence.** The g^-2 gain holds in the same c_T -> infinity (A2') limit as the paper's tau_*. At physical c_T it needs
  E_light << 1e-58: an exponential c_T threshold. That is exactly the property the angular route was valued for avoiding.
* **Light-robust LP** (`lp_light.py`, `lp_light2.py`, `run_light_*.py`). Impose |Lambda| <= kappa W on light states. The
  bound is then
      tau^2 <= [48 + 16(kappa-1) E'_light]/L,  with E'_light the light share of the endpoint rule (<= 3).
* **Robust results.** At light coefficient 16 kappa/L ~ 460–1330 (the paper's own is ~16*70 = 1120):
  * the optimum needs only n ~ 12;
  * it gives L/g^2 = 7.5, 9.7, 12.0 at g = 36, 60, 100;
  * so tau g ~ 2.5 -> 2.0, against the paper's 17.8.
  The scaling is essentially 1/g with an ~8–9x better constant. Power gains in g appear only when kappa (hence c_T) grows
  faster than any fixed power of g.
