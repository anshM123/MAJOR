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
