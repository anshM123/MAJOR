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
