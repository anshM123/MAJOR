# angular_g4: new positivity-free identities for tau = sqrt(t4^2 + alpha_T^2) beyond the k-rules

Research ledger (agent run, 28 Nov 2026). Only files in this directory were written. Nothing was committed.

**Grades.**

| grade | meaning |
|---|---|
| [P] | exact algebra or proof |
| [P\|X] | proved given X |
| [C] | exact computer check at the stated finite set |
| [NE] | numerical (float LP) evidence |
| [D] | derived, not independently re-derived |
| [open] | open |

**Inputs.** Everything here uses the same inputs as the paper's endpoint rule (supplement sec:angular):
* the all-J ladder (LADDER_ALLJ);
* the first-moment matrix M = (i/6)AH (LADDER_ALLJ §§2–3);
* the LADDER_ALLJ §4 mixed-commutator identity with its D-operator;
* the contact iδ(n,n')E(n) of the mixed commutator (KZ algebra; ANGULAR_ROUTE §15; ANGULAR_MATCH §3);
* the mirror relation â_{m,−2}(λ_e, λ_o) = â_{−m,2}(λ_e, −λ_o) (KERNELS_D3 §5.1; LADDER_ALLJ §1: "for h = −2 everything is mirrored,
  m → −m, times σ = ±1 by parity").

No positivity of any heavy functional is used. The only inequalities are |ξ| ≤ sqrt(ρ0 ρ4) (Cauchy–Schwarz on one state's own amplitudes)
and the endpoint budget.

**Notation.** As in `code/angular_opt/README.md`: T = (Δ−1)(Δ−2), C = J(J+1), r_j^{(b)}, W^{(b)} = (T − b(b−1))/3 + 3/2, and
|c|² = τ²/16. Also d_j = (Aa)_j/a_j = 3[T − (s−1)s]/(j+3) − 3[C − s(s+1)]/(j+1), with s = b + j. This d_j formula is checked against the
tridiagonal A to 1e−15 (`core.py`).

## Milestone 1 — the per-harmonic ("k-resolved") endpoint rules M_n (28 Nov)

### 1a. The identity [P | paper's inputs]
* **Contact fixes every Fourier coefficient.** The state-summed mixed commutator C(φ) = ⟨T,+2|[E(θ₁), E⁽¹⁾(θ₂)]|T,+2⟩ is c₀δ(φ)
  (ANGULAR_MATCH §3: c₂ = 0 by the Ward step), with Ĉ(0) = 1. Hence **Ĉ(n) = 1 for every n**, not only Σ_n q_n Ĉ(n) = 0.
* **Per-harmonic coefficients.** The coefficients are those of `code/angmatch/am3_ward.py`, which is exact. Per unit |a_0|², a ladder
  state contributes to Ĉ(n), n ≥ 2:

      h_n = ½ r_{n−2} (1 + d_{n−2}/3) + (n+1)/4 · r_{n−3} − (n−1)/4 · r_{n−1}          (both ladders, same form)

* **The flip entry of T.** |c|² adds γ_3 = −½, γ_4 = +½, γ_5 = +3/2.
* **The rules:**

      (M_n)   Σ_heavy ρ h_n + γ_n |c|² = 1 ,   n = 2, 3, 4, ...

* **Relation to the endpoint rule.** The endpoint rule is Σ_n 2q_n (M_n − 1), by telescoping. It is the only combination that needs the
  harmonic/state Fubini exchange, and its defect is 𝒟. **Each M_n is a finite-harmonic statement and carries no defect 𝒟**: each needs
  only the convergence of the state sum at fixed n.
* **Free-field check [C].** With all |a|² = 1, the free scalar (φ², all even currents, |c|² = 1) gives Ĉ(n) = 1 for n = 0..40
  (`core.py`).
  * Free fields have A·v = 0 on every state, so this tests the D-term and T-part, not the A-term. The A-term is the paper's.
* **Effect on the 3-state no-go.** M_n for n > 8 has zero heavy contribution from {A, B, C}, but RHS 1. So the configuration {A, B, C}
  is **excluded**: every admissible spectrum needs heavy support at arbitrarily high harmonics.
  * This exclusion is cheap. A high-J background can supply the M_n.

### 1b. Finite-harmonic obstruction [P, = ANGULAR_ROUTE §15]
* For large T at fixed x = C/T, h_n → (T/2) D_{n−2}² x^{n−2}[1/(n+1) − x/(n−1)], and ∫₀¹(1−x)(·)dx = 0.
* So **no finite combination of M_n alone** has a nonnegative leading heavy action. The M_n are useful only on top of the endpoint
  rule (multiplier ν > 0).

### 1c. LP results, heavy states only [NE] (`lp1.py` primal, `lp2.py` dual; they agree to 1%)
The rule set is {endpoint, R_k (k ≤ K), M_n (n ≤ N)}, with heavy states only (twist ≥ g) and no light constraint.

| g | old (no M) τ | with M_n τ | N | τ g³ new |
|---|---|---|---|---|
| 10 | 0.2841 | 0.1468 | 16–32 (saturated) | 147 |
| 20 | 0.04965 | 0.03092 | 40 (still falling slowly) | 247 |

* Ratio new/old: 0.52 at g = 10 and 0.62 at g = 20.
* The 3-state primal of the old LP is replaced by a primal that uses high-J support. The exponent stays consistent with **g⁻³**: τg³
  grows towards a constant, as in the old LP. The gain is a constant factor of about 0.5–0.6.

### 1d. LP results, light-robust [NE] (`lp2.py` mode theta, `escan.py`, `thetascan.py`)
* **Constraints.** Light states (2 < twist < g, all J) obey F ≥ −Θ W. The bound is τ² ≤ (48 + 16Θ E'_light)/L, with light coefficient
  16Θ/L.
* **Without M_n,** L(Θ = g^p) ∝ g^{e(p)}. The exponent is stable over g = 20…160:

  | p | e(p) | τ ∝ g^{−e/2} | light coefficient ∝ g^{p−e} |
  |---|---|---|---|
  | 3 | 2.42–2.45 | g^{−1.22} | g^{0.57} (130–300) |
  | 3.5 | 2.58–2.62 | g^{−1.30} | g^{0.9} (225–1700) |
  | 4 | 2.76–2.77 | g^{−1.38} | g^{1.24} (640–15000) |
  | 5 | 3.06–3.09 (to g = 60) | g^{−1.54} | g^{1.94} |

  The value at p = 5, g ≥ 100 is saturated by the K = 32 truncation.
* **With M_n (N = 24),** L is 10–20% larger and e(p) is unchanged (e.g. p = 4: 2.80 against 2.77). M_n are not a light-robust power
  gain.
* **Correction to "reverts to 1/g".** A light coefficient polynomial in g buys a power better than 1/g: τ ~ g^{−e(p)/2}, with e(p) > 2 for
  p > ~2.7. Each extra power of the bound costs more than one extra power of g in the light coefficient. There is no exponential
  threshold. The other agent's "≈1/g" corresponds to a light coefficient that is constant in g.

### 1e. Heavy-only exponent check at larger g [NE] (`heavyscan.py`, dual LP, J up to 3000 g)
* **Values.** τ g³ without and with M_n, and their ratio:

  | g | without M_n | with M_n | ratio |
  |---|---|---|---|
  | 30 | 455 | 349 | 0.77 |
  | 40 | 671 | 577 | 0.86 |

* **g = 40 is not trustworthy.** Its "without" value exceeds the 3-state asymptote 588. Double precision is not enough for the
  unregularised dual at g ≥ 40, so treat this point as indicative only.
* **Trend.** The ratio is 0.52, 0.62, 0.77, 0.86 at g = 10, 20, 30, 40, rising towards 1. **The M_n gain disappears asymptotically: no
  change of the g⁻³ exponent, and probably not of its constant.**
* **Heuristic reason.** The {A, B, C} core puts about ½ of its budget into M_2, which is more than its RHS of 1. A cheap high-J
  background with h_2/W → −1 (x = C/T → 1) absorbs the excess at a budget fraction that is O(1) at small g and falls at large g.
* The g = 50 primal run (`scan1.py`) violated tolerances and is discarded.

## Milestone 2 — the helicity-flip channel ⟨T,−2|[E, E⁽¹⁾]|T,+2⟩ (28 Nov)

### 2a. Derivation [P algebra (`flip_check.py`); D for the physics inputs, which are the paper's]
* **Setup.** In the LADDER_ALLJ §4 identity, replace the bra amplitudes by u_m = ⟨O,m|E|T,−2⟩, with the mirror relation
  u_m = η_O v_{−m} (η = ±1 by parity; parity-invariant theory). Set G' = (i/6)A F'. The D-operator is unchanged: it is the boost
  [K, E], with no dependence on the helicity of the external states.
* **Expansion.** Write LHS = −i e^{2i(θ₁+θ₂)} Σ_n X_n e^{inφ}. Then exactly (random-vector check to 8e−14):

      X_n = (1/6)(Q_{−n} + Q'_n) + ½[S_n + (n+3)/2 S_{n−1} − (n+1)/2 S_{n+1}],
      P_m = ū_m v_m,   S_m = P_m + P_{−m},   Q_m = conj((Au)_m) v_m,   Q'_m = ū_m (Av)_m .

* **The rules.** The contact iδE gives **X_n = c for every n**. The T family contributes c at n = −3, −2, −1, 2, 3c at n = 3, and −c at
  n = 1.
* **Free-field check [C].** Take the free scalar (φ², every even current with a_{J,+}a_{J,−} = +1, c = 1). It satisfies X_n = c for
  all |n| ≤ 28.
  * The free fermion is the mirror image: ψ̄ψ is parity-odd, the λ_{++} have opposite sign, and c = −1.
  * This tests the D-term, the T-part and the contact. As in 1a, it does not test the A-term.
* **Heavy content.** Let χ_p = ξ sqrt(r0_p r4_{p−4}), with ξ = η α β. The cross term needs a parity-even, even-J ≥ 4 state with both
  λ_{++} and λ_{+−}. Then:
  * **F_0:** Σ η ρ0 (1 + (T − 3C)/3) = c.
  * **F_1:** Σ η ρ0 (2 − 4C/3) = 2c. This uses the off-ladder first-moment component (Av)_1 = −sqrt(CT) v_0; the pipeline confirms
    M_1 ≠ 0 (LADDER_ALLJ §3).
  * **F_3:** Σ ξ[−(4/21)(T − 12) − 2] χ_4/ξ = −2c.
  * **F_{±p}, p ≥ 4:** the cross terms; see `flip.py`.
* **Corollary, a new exact parity sum rule** (2F_0 − F_1):

      Σ_{O ≠ T} η_O |a_{O,0}|² (T_Δ − J(J+1)) = 0 ,     t4/4 = Σ_O η_O |a_{O,0}|² (1 − 2J(J+1)/3)

  * Here |a_{O,0}|² is the b = 0 head (helicity-0 component) and η_O is the parity.
  * Twist > 2 means T − C > 0. So, without operators of twist below 2 (or scalars with Δ < 2), **the b = 0 couplings of the two
    parities must balance**.
  * t4 is then a signed, parity-weighted, linear sum over b = 0 heads.

### 2b. Flip channel alone [NE] (`fliplp.py`)
* **The bound.** |c| ≤ 3t, where |G_κ(s)| ≤ 2t sqrt(W0 W4) on heavy states (endpoint budget plus Cauchy–Schwarz).
* **Values.** τ g ≈ 14.9 (g = 10), 8.2 (20), 7.9 (40), and 14.3 (80, N = 48, not converged in N). So the flip channel alone gives about
  1/g.

### 2c. Exact obstruction for the flip channel [P, exact rationals, `unenh.py`]
* **Conversion identity.** On every state, the cross ladder satisfies
  (C − s(s+1)) χ_s = μ_s (T − s(s+1)) χ_{s+1}, with μ_s = (s+1)(s−3)/((s+4)s). So C-weights convert into T-weights.
* **Unenhanced combinations exist.** A combination Σκ_n F_n is "unenhanced" if, after conversion, it has no T-coefficient. These
  combinations exist: the null space has dimension N − 3.
* **But none of them sees c.** For every unenhanced combination, Σκ_n rhs_n = 0. Exactly, rhs = yᵀM with y_p = (p+3)(p−1)/2 (checked
  N = 6…14).
* **Consequence.** Every c-sensitive flip functional carries T-enhanced heavy weights (∝ T or C ~ g²). A linear bound
  |c| ≲ max(weight)/W can therefore not beat O(1)·x² from this mechanism. The power must come from cancellations, as in the diagonal
  channel.

### 2d. Combined diagonal + flip, light-robust [NE] (`combo2.py`, `combo_run.py`; `combo_p4.txt`, `combo_heavy.txt`)
* **Method.** This is an LP inner approximation of the per-state 2×2 PSD dual condition, |Ĝ| ≤ 2F̂0/s and ≤ 2sF̂4, with s updated
  iteratively; each LP gives a valid bound. Bisection is on c = t4/4 (parity-invariant theory).
* **Light-robust results, Θ = g⁴, K = 24, flip N_f = 24:**

  | g | diagonal only, τ | + flip, τ | gain |
  |---|---|---|---|
  | 10 | 0.3327 | 0.3223 | 3.1% |
  | 20 | 0.1097 | 0.1089 | 0.7% |
  | 36 | 0.04917 | (run not finished) | |

* **Verdict.** **The flip channel's gain vanishes with g. It does not change e(p) and does not break the g⁻³ wall.**
* **Why it is asymptotically free.** Any diagonal-feasible configuration can be made flip-feasible as follows:
  * split every b = 0-only state into two equal halves of opposite parity, which zeroes F_0, F_1 and Φ;
  * add a "c-generator" of parity-even states that carry both ladders.
  * The c-generator's endpoint cost is O(|c|·g) by 2b. It carries both ladders, and their leading terms cancel in R_k, so its R_k
    pollution is about ρ·k²/T ~ |c|/T² ~ g⁻⁷ at |c| ~ g⁻³. That is ≪ |c|² ~ g⁻⁶, so it is harmless.
* **Heavy-only run with flip (g = 10): an artefact, not a loss of information.** It gave τ = 0.216 against 0.1454 without flip. The
  dual multipliers are huge there, and the iterated PSD inner approximation (2 s-updates) is poor. The bound is valid but not
  optimal. A proper SOCP (Clarabel) failed numerically. The heavy-only flip effect is therefore covered by the 2d argument, not by
  this number.

## Milestone 3 — conclusions and structural no-go statements (28 Nov)

### 3a. What every identity considered here is [P]
* **Per-state data is complete.** By KERNELS_D3 §5.5 (the number of independent amplitudes equals the number of conserved structures),
  the full rest-frame transition data of one O is (λ_{++}, λ_{+−}, λ_odd) times the all-J ladders.
* **So every rule is linear in the per-state Gram data (ρ0^±, ρ4, ξ).** This holds for any T-state matrix element of a commutator of
  two light-ray moments, with any helicities, harmonics or moments, given the ladder. Its kernels are polynomials in (T, C) of degree
  equal to the total moment order:
  * degree 0: [E, E]; this gives R_k (diagonal), and in the flip channel Im-parts that are trivial for parity-invariant real data;
  * degree 1: [E, E⁽¹⁾]; this gives the endpoint rule, M_n, F_0, F_1, F_3, F_{±p}, and Φ.
* **Higher degree is not available.** It would need [E⁽¹⁾, E⁽¹⁾] (Regge J₀ < 1, marginal at the nonperturbative J₀ ≤ 1) or moments
  E⁽²⁾ (J₀ < 0, false).
* **Two-detector event shapes ⟨T|E[E, E]|T⟩** (candidate 2) bring in heavy–heavy ⟨X|E|Y⟩ data, which is not determined by T-state
  data. Without positivity they do not close.

### 3b. No-go statements
1. **[P] The flip channel has no unenhanced c-sensitive combination** (2c): rhs = yᵀM, y_p = (p+3)(p−1)/2.
2. **[P] Finite-harmonic obstruction for M_n** (1b, = ANGULAR_ROUTE §15).
3. **[NE + heuristic] Degree-1 rules do not change the heavy-only g⁻³ asymptote.**
   * They are "asymptotically free": their RHS can be met by a high-J background (M_n) or a c-generator (flip) whose budget share
     vanishes as g → ∞.
   * Evidence: the M_n gain ratio is 0.52, 0.62, 0.77, 0.86 at g = 10, 20, 30, 40, and the flip gain is 3% then 0.7%.
   * **So the only identities that can move the g⁻³ wall are new degree-0 (unenhanced) rules.** In T-states the complete degree-0 set
     is {R_k}, plus the trivial flip Im-parts. The other agent's 3-state no-go then applies to the whole class. Within positivity-free
     T-state commutator identities with moments ≤ 1, **g⁻³ is the wall, and g⁻⁴ needs positivity or heavy–heavy data**.
4. **[argument] ANEC positivity** (candidate 3), used through T-state data only, gives:
   * the POVM bound E(f)² ≤ ‖f‖∞E(f);
   * 2×2 minors with ⟨O|E|O⟩;
   * the conformal collider bound on c.

   All of these are unenhanced caps of the form Σ|a|² ≤ 1, O(1) per state. They cannot produce g-powers. A g⁻⁴ via positivity needs
   heavy functionals of dispersive type (the paper's Ω^EE route).
5. **[P] Parity and Bose structure** (candidate 4) do not forbid {A, B, C}:
   * A = (b0, J = 2) is allowed in either parity; B and C are parity-even S^{+−} states.
   * The flip channel excludes the literal configuration: it forces c = 0 unless cross-term states exist.
   * It does not exclude the family: halve A into two parities and add a cheap c-generator.

### 3c. Positive results worth keeping
* **M_n.** These are D-free, finite-harmonic, exact per-harmonic rules: Ĉ(n) = 1 for all n. They give a constant-factor tightening at
  finite g, which disappears asymptotically.
* **The exact parity sum rule Φ** and the linear t4 rule:
  * Σ_O η_O |a_{O,0}|²(T_Δ − J(J+1)) = 0;
  * t4/4 = Σ_O η_O |a_{O,0}|²(1 − 2J(J+1)/3).
  * These are new and positivity-free [P | the paper's inputs + mirror relation + (Av)_1]. They are not a bound on τ.
* **Light-robust tradeoff e(p)** (1d): a polynomial light coefficient g^{p−e(p)} buys τ ~ g^{−e(p)/2}. This is a polynomial
  threshold, not an exponential one.
